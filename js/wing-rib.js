/* Fig. 2 — live finite-element analysis of a wing rib.
 * Geometry, meshing, assembly and solution all run in the browser:
 * plane-stress linear elasticity, constant-strain triangles (CST) and a
 * banded Cholesky solver. Units: mm, N, MPa. */
(function () {
  'use strict';

  // ---------- Model definition ----------
  const CHORD = 600;
  const NACA = { m: 0.02, p: 0.4, t: 0.18 };                       // NACA 2418
  const SPARS = [0.18, 0.68];                                        // x/c of front and rear spar
  const HOLES = [[0.07, 0.55], [0.30, 0.62], [0.44, 0.62], [0.57, 0.58], [0.80, 0.55]]; // x/c, radius / local half-thickness
  const MAT = { E: 71700, nu: 0.33, sy: 503, t: 1.6 };               // aluminium 7075-T6 web, 1.6 mm
  const H = 7.5;                                                     // target element size
  const PATCH = 40;                                                  // load spread along the edge
  const F_MAX = 12000;
  const F_DEFAULT = 3800;
  const BANDS = ['#2B4560', '#3E5E78', '#5A7A8E', '#8199A3', '#AEBDB6', '#D8D3BC', '#E3BC86', '#D4934F', '#B9652E', '#9C3D1F'];
  const BEYOND = '#5E1A0C';

  // ---------- Airfoil ----------
  function camber(x) {
    const { m, p } = NACA;
    return x < p ? m / (p * p) * (2 * p * x - x * x) : m / ((1 - p) * (1 - p)) * ((1 - 2 * p) + 2 * p * x - x * x);
  }
  function halfThickness(x) {
    return 5 * NACA.t * (0.2969 * Math.sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
  }
  const upper = x => (camber(x) + halfThickness(x)) * CHORD;
  const lower = x => (camber(x) - halfThickness(x)) * CHORD;

  function insidePoly(poly, x, y) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const [xi, yi] = poly[i], [xj, yj] = poly[j];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }
  function distToPoly(poly, x, y) {
    let best = Infinity;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const [ax, ay] = poly[j], [bx, by] = poly[i];
      const dx = bx - ax, dy = by - ay;
      const t = Math.max(0, Math.min(1, ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)));
      best = Math.min(best, Math.hypot(x - ax - t * dx, y - ay - t * dy));
    }
    return best;
  }

  // ---------- Geometry: outline, holes, spar lines and interior seeds ----------
  function buildGeometry() {
    // Dense outline (cosine spacing): trailing edge -> upper -> leading edge -> lower
    const M = 240, dense = [];
    for (let i = 0; i <= M; i++) {
      const x = 0.5 * (1 + Math.cos(Math.PI * i / M));
      dense.push({ x: x * CHORD, y: upper(x), forced: i === 0 });
    }
    for (let i = 1; i < M; i++) {
      const x = 0.5 * (1 - Math.cos(Math.PI * i / M));
      dense.push({ x: x * CHORD, y: lower(x), forced: false });
    }
    // Spar attachment points become exact outline vertices
    SPARS.forEach(xs => {
      const X = xs * CHORD;
      [upper(xs), lower(xs)].forEach(Y => {
        for (let k = 0; k < dense.length; k++) {
          const a = dense[k], b = dense[(k + 1) % dense.length];
          if ((a.x - X) * (b.x - X) < 0 && Math.abs((a.y + b.y) / 2 - Y) < Math.abs((a.y + b.y) / 2 - camber(xs) * CHORD)) {
            dense.splice(k + 1, 0, { x: X, y: Y, forced: true });
            break;
          }
        }
      });
    });

    // Resample at the element size between forced vertices
    const n = dense.length, s = [0];
    for (let i = 1; i <= n; i++) {
      const a = dense[i - 1], b = dense[i % n];
      s.push(s[i - 1] + Math.hypot(b.x - a.x, b.y - a.y));
    }
    const forced = [];
    dense.forEach((p, i) => { if (p.forced) forced.push(i); });
    const outline = [];
    forced.forEach((i0, f) => {
      const i1 = f + 1 < forced.length ? forced[f + 1] : n;
      const k = Math.max(1, Math.round((s[i1] - s[i0]) / H));
      let seg = i0;
      for (let j = 0; j < k; j++) {
        const target = s[i0] + (s[i1] - s[i0]) * j / k;
        while (s[seg + 1] < target) seg++;
        const a = dense[seg], b = dense[(seg + 1) % n], t = (target - s[seg]) / (s[seg + 1] - s[seg]);
        outline.push([a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)]);
      }
    });
    const densePoly = dense.map(p => [p.x, p.y]);

    const holes = HOLES.map(([xc, f]) => ({ cx: xc * CHORD, cy: camber(xc) * CHORD, r: f * halfThickness(xc) * CHORD }));
    const holePts = holes.map(h => {
      const k = Math.max(20, Math.ceil(2 * Math.PI * h.r / (0.75 * H)));
      return Array.from({ length: k }, (_, i) => [h.cx + h.r * Math.cos(2 * Math.PI * i / k), h.cy + h.r * Math.sin(2 * Math.PI * i / k)]);
    });

    const sparPts = [];
    SPARS.forEach(xs => {
      const yl = lower(xs), yu = upper(xs), k = Math.max(2, Math.round((yu - yl) / H));
      for (let j = 1; j < k; j++) sparPts.push([xs * CHORD, yl + (yu - yl) * j / k]);
    });

    // Interior seeds on a lightly jittered hexagonal lattice
    let seed = 7;
    const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    const ys = densePoly.map(p => p[1]), ymin = Math.min(...ys), ymax = Math.max(...ys);
    const dy = H * Math.sqrt(3) / 2, interior = [];
    for (let row = Math.floor(ymin / dy); row * dy <= ymax; row++) {
      for (let x = (row & 1) ? H / 2 : 0; x <= CHORD; x += H) {
        const px = x + (rnd() - 0.5) * 0.04 * H, py = row * dy + (rnd() - 0.5) * 0.04 * H;
        if (!insidePoly(densePoly, px, py) || distToPoly(densePoly, px, py) < 0.55 * H) continue;
        if (holes.some(h => Math.hypot(px - h.cx, py - h.cy) < h.r + 0.55 * H)) continue;
        if (SPARS.some(xs => Math.abs(px - xs * CHORD) < 0.5 * H)) continue;
        interior.push([px, py]);
      }
    }
    return { outline, densePoly, holes, holePts, sparPts, interior, ymin, ymax };
  }

  // ---------- Delaunay triangulation (Bowyer–Watson) ----------
  function circumcircle(P, a, b, c) {
    const [ax, ay] = P[a], [bx, by] = P[b], [cx, cy] = P[c];
    const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    const ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
    const uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
    return { a, b, c, ux, uy, r2: (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy) };
  }
  function delaunay(pts) {
    const n = pts.length, D = CHORD * 20;
    const P = pts.concat([[-D, -D], [CHORD / 2, D], [CHORD + D, -D]]);
    let tris = [circumcircle(P, n, n + 1, n + 2)];
    for (let i = 0; i < n; i++) {
      const [px, py] = P[i], keep = [], edges = new Map();
      for (const t of tris) {
        if ((px - t.ux) * (px - t.ux) + (py - t.uy) * (py - t.uy) < t.r2) {
          [[t.a, t.b], [t.b, t.c], [t.c, t.a]].forEach(([u, v]) => {
            const key = u < v ? u * 65536 + v : v * 65536 + u;
            if (edges.has(key)) edges.delete(key); else edges.set(key, [u, v]);
          });
        } else keep.push(t);
      }
      edges.forEach(([u, v]) => keep.push(circumcircle(P, u, v, i)));
      tris = keep;
    }
    return tris.filter(t => t.a < n && t.b < n && t.c < n).map(t => [t.a, t.b, t.c]);
  }

  // ---------- Mesh ----------
  function buildMesh(geo) {
    const pts = [...geo.outline, ...geo.holePts.flat(), ...geo.sparPts, ...geo.interior];
    const tris = delaunay(pts).filter(([a, b, c]) => {
      const cx = (pts[a][0] + pts[b][0] + pts[c][0]) / 3, cy = (pts[a][1] + pts[b][1] + pts[c][1]) / 3;
      return insidePoly(geo.densePoly, cx, cy) && !geo.holes.some(h => Math.hypot(cx - h.cx, cy - h.cy) < h.r);
    });

    // Renumber nodes along the chord so the stiffness matrix stays banded
    const used = new Uint8Array(pts.length);
    tris.forEach(t => t.forEach(i => { used[i] = 1; }));
    const order = [];
    pts.forEach((_, i) => { if (used[i]) order.push(i); });
    order.sort((i, j) => pts[i][0] - pts[j][0]);
    const map = new Int32Array(pts.length).fill(-1);
    order.forEach((old, k) => { map[old] = k; });

    const nn = order.length, X = new Float64Array(nn), Y = new Float64Array(nn);
    order.forEach((old, k) => { X[k] = pts[old][0]; Y[k] = pts[old][1]; });
    const T = new Int32Array(tris.length * 3);
    tris.forEach(([a, b, c], e) => {
      a = map[a]; b = map[b]; c = map[c];
      const twoA = (X[b] - X[a]) * (Y[c] - Y[a]) - (X[c] - X[a]) * (Y[b] - Y[a]);
      T.set(twoA > 0 ? [a, b, c] : [a, c, b], 3 * e);
    });

    const edgeSet = new Set();
    for (let e = 0; e < tris.length; e++) {
      for (let k = 0; k < 3; k++) {
        const u = T[3 * e + k], v = T[3 * e + (k + 1) % 3];
        edgeSet.add(u < v ? u * 65536 + v : v * 65536 + u);
      }
    }
    const edges = Int32Array.from([...edgeSet].flatMap(key => [Math.floor(key / 65536), key % 65536]));

    let offset = geo.outline.length;
    const outlineNodes = geo.outline.map((_, i) => map[i]);
    const holeNodes = geo.holePts.map(h => { const ids = h.map((_, i) => map[offset + i]); offset += h.length; return ids; });
    const fixed = new Uint8Array(nn);
    for (let i = 0; i < nn; i++) if (SPARS.some(xs => Math.abs(X[i] - xs * CHORD) < 1e-6)) fixed[i] = 1;

    // Arc length along the outer edge, used to spread point loads
    const outlineS = [0];
    for (let i = 1; i < outlineNodes.length; i++) {
      const a = outlineNodes[i - 1], b = outlineNodes[i];
      outlineS.push(outlineS[i - 1] + Math.hypot(X[b] - X[a], Y[b] - Y[a]));
    }
    const a0 = outlineNodes[outlineNodes.length - 1], b0 = outlineNodes[0];
    const perimeter = outlineS[outlineS.length - 1] + Math.hypot(X[b0] - X[a0], Y[b0] - Y[a0]);

    return { X, Y, T, edges, nn, ne: tris.length, outlineNodes, outlineS, perimeter, holeNodes, fixed };
  }

  // ---------- Assembly and factorisation ----------
  function assemble(mesh) {
    const { X, Y, T, nn, ne, fixed } = mesh, ndof = 2 * nn;
    let hb = 0;
    for (let e = 0; e < ne; e++) {
      const a = T[3 * e], b = T[3 * e + 1], c = T[3 * e + 2];
      hb = Math.max(hb, 2 * (Math.max(a, b, c) - Math.min(a, b, c)) + 1);
    }
    const w = hb + 1, K = new Float64Array(ndof * w);
    const { E, nu, t } = MAT, k0 = E / (1 - nu * nu);
    const D = [k0, k0 * nu, 0, k0 * nu, k0, 0, 0, 0, k0 * (1 - nu) / 2];
    const B = new Float64Array(ne * 18), area = new Float64Array(ne), DB = new Float64Array(18);

    for (let e = 0; e < ne; e++) {
      const n1 = T[3 * e], n2 = T[3 * e + 1], n3 = T[3 * e + 2];
      const twoA = (X[n2] - X[n1]) * (Y[n3] - Y[n1]) - (X[n3] - X[n1]) * (Y[n2] - Y[n1]);
      const b1 = Y[n2] - Y[n3], b2 = Y[n3] - Y[n1], b3 = Y[n1] - Y[n2];
      const c1 = X[n3] - X[n2], c2 = X[n1] - X[n3], c3 = X[n2] - X[n1];
      const Be = B.subarray(18 * e, 18 * e + 18);
      Be.set([b1, 0, b2, 0, b3, 0, 0, c1, 0, c2, 0, c3, c1, b1, c2, b2, c3, b3].map(v => v / twoA));
      area[e] = twoA / 2;
      for (let r = 0; r < 3; r++) {
        for (let s = 0; s < 6; s++) DB[6 * r + s] = D[3 * r] * Be[s] + D[3 * r + 1] * Be[6 + s] + D[3 * r + 2] * Be[12 + s];
      }
      const dofs = [2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1];
      for (let r = 0; r < 6; r++) {
        for (let s = 0; s < 6; s++) {
          const I = dofs[r], J = dofs[s];
          if (I < J) continue;
          K[I * w + (I - J)] += t * area[e] * (Be[r] * DB[s] + Be[6 + r] * DB[6 + s] + Be[12 + r] * DB[12 + s]);
        }
      }
    }

    // Spar attachments: clamp both displacement components
    for (let node = 0; node < nn; node++) {
      if (!fixed[node]) continue;
      for (const i of [2 * node, 2 * node + 1]) {
        for (let j = Math.max(0, i - hb); j <= i; j++) K[i * w + (i - j)] = 0;
        for (let r = i; r <= Math.min(ndof - 1, i + hb); r++) K[r * w + (r - i)] = 0;
        K[i * w] = 1;
      }
    }

    // In-place banded Cholesky: K = L L^T
    for (let i = 0; i < ndof; i++) {
      const j0 = Math.max(0, i - hb);
      for (let j = j0; j <= i; j++) {
        let sum = K[i * w + (i - j)];
        for (let k = Math.max(j0, j - hb); k < j; k++) sum -= K[i * w + (i - k)] * K[j * w + (j - k)];
        K[i * w + (i - j)] = i === j ? Math.sqrt(sum) : sum / K[j * w];
      }
    }
    return { L: K, w, hb, ndof, B, area, D };
  }

  function solve(sys, f) {
    const { L, w, hb, ndof } = sys, u = Float64Array.from(f);
    for (let i = 0; i < ndof; i++) {
      let sum = u[i];
      for (let k = Math.max(0, i - hb); k < i; k++) sum -= L[i * w + (i - k)] * u[k];
      u[i] = sum / L[i * w];
    }
    for (let i = ndof - 1; i >= 0; i--) {
      let sum = u[i];
      for (let k = i + 1; k <= Math.min(ndof - 1, i + hb); k++) sum -= L[k * w + (k - i)] * u[k];
      u[i] = sum / L[i * w];
    }
    return u;
  }

  // Area-weighted nodal average of element stresses (sx, sy, txy)
  function nodalStress(mesh, sys, u) {
    const { T, nn, ne } = mesh, { B, area, D } = sys;
    const S = new Float64Array(3 * nn), wsum = new Float64Array(nn), eps = [0, 0, 0];
    for (let e = 0; e < ne; e++) {
      const n = [T[3 * e], T[3 * e + 1], T[3 * e + 2]];
      const ue = [u[2 * n[0]], u[2 * n[0] + 1], u[2 * n[1]], u[2 * n[1] + 1], u[2 * n[2]], u[2 * n[2] + 1]];
      for (let r = 0; r < 3; r++) {
        let v = 0;
        for (let s = 0; s < 6; s++) v += B[18 * e + 6 * r + s] * ue[s];
        eps[r] = v;
      }
      for (let r = 0; r < 3; r++) {
        const sig = D[3 * r] * eps[0] + D[3 * r + 1] * eps[1] + D[3 * r + 2] * eps[2];
        for (const node of n) S[3 * node + r] += area[e] * sig;
      }
      for (const node of n) wsum[node] += area[e];
    }
    for (let i = 0; i < nn; i++) for (let r = 0; r < 3; r++) S[3 * i + r] /= wsum[i];
    return S;
  }

  // Unit load (1 N in x, 1 N in y) spread over the edge around outline node k
  function unitCases(mesh, sys, k) {
    const { outlineNodes, outlineS, perimeter, fixed } = mesh;
    const fx = new Float64Array(sys.ndof), fy = new Float64Array(sys.ndof);
    const weights = [];
    let total = 0;
    outlineNodes.forEach((node, j) => {
      if (fixed[node]) return;
      let d = Math.abs(outlineS[j] - outlineS[k]);
      d = Math.min(d, perimeter - d);
      if (d >= PATCH) return;
      const wgt = Math.cos(Math.PI * d / (2 * PATCH)) ** 2;
      weights.push([node, wgt]);
      total += wgt;
    });
    weights.forEach(([node, wgt]) => { fx[2 * node] = wgt / total; fy[2 * node + 1] = wgt / total; });
    const ux = solve(sys, fx), uy = solve(sys, fy);
    return { ux, uy, sx: nodalStress(mesh, sys, ux), sy: nodalStress(mesh, sys, uy) };
  }

  const core = { CHORD, MAT, SPARS, F_MAX, F_DEFAULT, BANDS, BEYOND, upper, lower, buildGeometry, buildMesh, assemble, solve, nodalStress, unitCases };
  if (typeof document === 'undefined') {
    if (typeof module !== 'undefined') module.exports = core;
    return;
  }

  // ---------- Figure ----------
  const figure = document.getElementById('ribFigure');
  if (!figure) return;
  const stage = figure.querySelector('.fe-stage');
  const canvas = document.getElementById('feCanvas');
  const ctx = canvas.getContext('2d');
  const logEl = document.getElementById('feLog');
  const out = {
    load: document.getElementById('feLoad'),
    smax: document.getElementById('feSmax'),
    dmax: document.getElementById('feDmax'),
    rf: document.getElementById('feRf')
  };
  const meshBtn = document.getElementById('feMesh');
  const resetBtn = document.getElementById('feReset');
  const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  const TEXT = {
    en: {
      meshing: 'Meshing…',
      meshed: (n, e) => `${n} nodes · ${e} CST elements`,
      factor: (d, hb) => `${d} DOF · banded Cholesky, half-bandwidth ${hb}`,
      done: (n, e, ms) => `${n} nodes · ${e} elements · meshed and solved in ${ms} ms`,
      front: 'FRONT SPAR', rear: 'REAR SPAR',
      scale: k => `Deformation ×${k}`,
      hint: 'Drag the edge to apply a load',
      aria: 'Interactive finite-element analysis of a wing rib. Drag on the edge to apply a load; arrow keys rotate and scale it.'
    },
    fr: {
      meshing: 'Maillage…',
      meshed: (n, e) => `${n} nœuds · ${e} éléments CST`,
      factor: (d, hb) => `${d} DDL · Cholesky en bande, demi-largeur ${hb}`,
      done: (n, e, ms) => `${n} nœuds · ${e} éléments · maillé et résolu en ${ms} ms`,
      front: 'LONGERON AVANT', rear: 'LONGERON ARRIÈRE',
      scale: k => `Déformée ×${k}`,
      hint: 'Faites glisser le bord pour charger',
      aria: 'Analyse par éléments finis interactive d’une nervure d’aile. Faites glisser le bord pour appliquer une charge ; les flèches du clavier la font tourner et varier.'
    }
  };
  const lang = () => (document.documentElement.lang === 'fr' ? 'fr' : 'en');
  const fmt = v => v.toLocaleString(lang() === 'fr' ? 'fr-FR' : 'en-GB');

  let geo = null;        // outline, holes and seed points
  let geoTime = 0;
  let model = null;      // mesh, factorised system, timings
  let view = null;       // canvas layout
  const cases = new Map();
  const state = {
    k: 0, Fx: 0, Fy: F_DEFAULT,
    showMesh: true, interacted: false, dragging: false,
    meshP: 1, contourP: 1, introDone: false, visible: false, raf: 0, t0: 0
  };

  function build() {
    const t0 = performance.now();
    const mesh = buildMesh(geo);
    const t1 = performance.now();
    const sys = assemble(mesh);
    const t2 = performance.now();
    model = { mesh, sys, tMesh: geoTime + t1 - t0, tSolve: t2 - t1 };
    // Default load: lift on the trailing-edge panel, aft of the rear spar
    const target = 0.88 * CHORD;
    let best = Infinity;
    mesh.outlineNodes.forEach((node, j) => {
      const d = Math.abs(mesh.X[node] - target);
      if (mesh.Y[node] > 0 && d < best && !mesh.fixed[node]) { best = d; state.k = j; }
    });
    model.defaultK = state.k;
    const t3 = performance.now();
    const c = getCase(state.k);
    model.tSolve += performance.now() - t3;
    // Fixed deformation magnification, picked so the default case reads clearly
    let dmax = 0;
    for (let i = 0; i < mesh.nn; i++) dmax = Math.max(dmax, Math.hypot(c.uy[2 * i], c.uy[2 * i + 1]) * F_DEFAULT);
    const raw = 20 / dmax, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    model.defScale = [1, 2, 5, 10].map(m => m * mag).filter(v => v <= raw).pop();
  }

  function getCase(k) {
    if (!cases.has(k)) cases.set(k, unitCases(model.mesh, model.sys, k));
    return cases.get(k);
  }

  function layout() {
    const cssW = stage.clientWidth;
    const pad = cssW < 640 ? 16 : 44;
    const s = (cssW - 2 * pad) / CHORD;
    const arrowMax = Math.max(70, Math.min(140, cssW * 0.12));
    const ribH = (geo.ymax - geo.ymin) * s;
    const cssH = Math.round(ribH + 2 * (arrowMax * 0.62 + 22));
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    canvas.style.height = cssH + 'px';
    view = { cssW, cssH, dpr, s, arrowMax, ox: pad, oy: cssH / 2 + (geo.ymax + geo.ymin) / 2 * s };
  }

  // Combine the two unit solutions for the current load
  function evaluate() {
    const { mesh } = model, c = getCase(state.k), { Fx, Fy } = state;
    const g = state.introDone ? 1 : state.contourP;
    const nn = mesh.nn, U = new Float64Array(2 * nn), vm = new Float64Array(nn);
    let smax = 0, dmax = 0;
    for (let i = 0; i < nn; i++) {
      const ux = g * (Fx * c.ux[2 * i] + Fy * c.uy[2 * i]);
      const uy = g * (Fx * c.ux[2 * i + 1] + Fy * c.uy[2 * i + 1]);
      U[2 * i] = ux; U[2 * i + 1] = uy;
      const sx = g * (Fx * c.sx[3 * i] + Fy * c.sy[3 * i]);
      const sy = g * (Fx * c.sx[3 * i + 1] + Fy * c.sy[3 * i + 1]);
      const txy = g * (Fx * c.sx[3 * i + 2] + Fy * c.sy[3 * i + 2]);
      vm[i] = Math.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy);
      smax = Math.max(smax, vm[i]);
      dmax = Math.max(dmax, Math.hypot(ux, uy));
    }
    return { U, vm, smax, dmax };
  }

  // Sutherland–Hodgman clip of a polygon against a level of a linear field
  function clip(poly, vals, level, above) {
    const p = [], v = [], n = vals.length;
    for (let i = 0; i < n; i++) {
      const j = (i + 1) % n, vi = vals[i], vj = vals[j];
      const inI = above ? vi >= level : vi <= level, inJ = above ? vj >= level : vj <= level;
      if (inI) { p.push(poly[2 * i], poly[2 * i + 1]); v.push(vi); }
      if (inI !== inJ) {
        const t = (level - vi) / (vj - vi);
        p.push(poly[2 * i] + t * (poly[2 * j] - poly[2 * i]), poly[2 * i + 1] + t * (poly[2 * j + 1] - poly[2 * i + 1]));
        v.push(level);
      }
    }
    return [p, v];
  }

  function draw() {
    if (!view) return;
    const { s, ox, oy, dpr, cssW, cssH } = view;
    const ink = 'rgba(28,27,24,', t = TEXT[lang()];
    const analysed = !!model && state.contourP > 0;

    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.lineJoin = 'round';

    // Undeformed shape: a solid blueprint until results are shown, then dashed
    if (analysed) ctx.setLineDash([4, 4]);
    ctx.strokeStyle = ink + (analysed ? '0.35)' : '0.85)');
    ctx.lineWidth = 1;
    ctx.beginPath();
    geo.densePoly.forEach(([x, y], i) => (i ? ctx.lineTo : ctx.moveTo).call(ctx, ox + x * s, oy - y * s));
    ctx.closePath();
    geo.holes.forEach(h => { ctx.moveTo(ox + (h.cx + h.r) * s, oy - h.cy * s); ctx.arc(ox + h.cx * s, oy - h.cy * s, h.r * s, 0, 2 * Math.PI); });
    ctx.stroke();
    ctx.setLineDash([]);

    let res = null, k = 1, px = null, py = null;
    if (model) {
      const { mesh } = model;
      res = evaluate();
      // Fixed magnification, reduced automatically for very large loads (as post-processors do)
      k = res.dmax > 0 ? Math.min(model.defScale, 28 / res.dmax) : model.defScale;
      px = i => ox + (mesh.X[i] + res.U[2 * i] * k) * s;
      py = i => oy - (mesh.Y[i] + res.U[2 * i + 1] * k) * s;

      // Banded von Mises contours
      if (analysed) {
        const levels = BANDS.map((_, b) => MAT.sy * b / BANDS.length).concat([MAT.sy, Infinity]);
        const paths = levels.slice(0, -1).map(() => new Path2D());
        for (let e = 0; e < mesh.ne; e++) {
          const n = [mesh.T[3 * e], mesh.T[3 * e + 1], mesh.T[3 * e + 2]];
          const vals = n.map(i => res.vm[i]);
          const poly = [px(n[0]), py(n[0]), px(n[1]), py(n[1]), px(n[2]), py(n[2])];
          const vmin = Math.min(...vals), vmax = Math.max(...vals);
          for (let b = 0; b < paths.length; b++) {
            const lo = levels[b], hi = levels[b + 1];
            if (vmax < lo || vmin > hi) continue;
            let p = poly, v = vals;
            if (lo > vmin) [p, v] = clip(p, v, lo, true);
            if (hi < vmax && p.length >= 6) [p, v] = clip(p, v, hi, false);
            if (p.length < 6) continue;
            paths[b].moveTo(p[0], p[1]);
            for (let q = 2; q < p.length; q += 2) paths[b].lineTo(p[q], p[q + 1]);
            paths[b].closePath();
          }
        }
        ctx.globalAlpha = state.introDone ? 1 : state.contourP;
        ctx.lineWidth = 0.6;
        paths.forEach((path, b) => {
          const color = b < BANDS.length ? BANDS[b] : BEYOND;
          ctx.fillStyle = color; ctx.strokeStyle = color;
          ctx.fill(path); ctx.stroke(path);
        });
        ctx.globalAlpha = 1;
      }

      // Mesh (swept in during the intro)
      if (state.showMesh || !state.introDone) {
        const limit = state.meshP * CHORD * 1.05;
        ctx.strokeStyle = analysed ? 'rgba(248,245,238,0.16)' : ink + '0.45)';
        ctx.lineWidth = 0.5;
        ctx.beginPath();
        for (let e = 0; e < mesh.edges.length; e += 2) {
          const a = mesh.edges[e], b = mesh.edges[e + 1];
          if ((mesh.X[a] + mesh.X[b]) / 2 > limit) continue;
          ctx.moveTo(px(a), py(a));
          ctx.lineTo(px(b), py(b));
        }
        ctx.stroke();
      }

      // Deformed edges
      ctx.strokeStyle = ink + '0.9)';
      ctx.lineWidth = 1.1;
      ctx.beginPath();
      [mesh.outlineNodes, ...mesh.holeNodes].forEach(loop => {
        loop.forEach((node, i) => (i ? ctx.lineTo : ctx.moveTo).call(ctx, px(node), py(node)));
        ctx.closePath();
      });
      ctx.stroke();
    }

    // Spars: clamped supports
    ctx.font = '500 10px "IBM Plex Mono", monospace';
    ctx.textAlign = 'center';
    SPARS.forEach((xs, i) => {
      const x = ox + xs * CHORD * s, yTop = oy - geo.ymax * s - 10, yBot = oy - geo.ymin * s + 10;
      ctx.strokeStyle = ink + '0.95)';
      ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.moveTo(x, yTop); ctx.lineTo(x, yBot); ctx.stroke();
      ctx.lineWidth = 1;
      [yTop, yBot].forEach((y, end) => {
        const dir = end ? 1 : -1;
        ctx.beginPath();
        ctx.moveTo(x - 9, y); ctx.lineTo(x + 9, y);
        for (let h = -9; h < 9; h += 4.5) { ctx.moveTo(x + h, y); ctx.lineTo(x + h + 4.5, y + dir * 5); }
        ctx.stroke();
      });
      ctx.fillStyle = ink + '0.55)';
      ctx.fillText(i ? t.rear : t.front, x, yTop - 10);
    });
    if (!model) return;

    // Load arrow
    const g = state.introDone ? 1 : state.contourP;
    const F = Math.hypot(state.Fx, state.Fy) * g;
    if (analysed && F > 1) {
      const node = model.mesh.outlineNodes[state.k];
      const x0 = px(node), y0 = py(node), len = F / F_MAX * view.arrowMax;
      const ux = state.Fx / Math.hypot(state.Fx, state.Fy), uy = -state.Fy / Math.hypot(state.Fx, state.Fy);
      const x1 = x0 + ux * len, y1 = y0 + uy * len;
      ctx.strokeStyle = '#9C3D1F'; ctx.fillStyle = '#9C3D1F';
      ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1 - ux * 8, y1 - uy * 8); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(x1, y1);
      ctx.lineTo(x1 - ux * 11 - uy * 5, y1 - uy * 11 + ux * 5);
      ctx.lineTo(x1 - ux * 11 + uy * 5, y1 - uy * 11 - ux * 5);
      ctx.closePath(); ctx.fill();
      ctx.beginPath(); ctx.arc(x0, y0, 3, 0, 2 * Math.PI); ctx.fill();
      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.textAlign = 'left';
      const label = `${(F / 1000).toFixed(2)} kN`, tw = ctx.measureText(label).width;
      const lx = Math.max(4, Math.min(cssW - tw - 4, ux >= 0 ? x1 + 8 : x1 - 8 - tw));
      const ly = Math.max(12, Math.min(cssH - 4, y1 + (uy > 0 ? 14 : -6)));
      ctx.fillText(label, lx, ly);
      view.head = [x1, y1];
    }

    // Labels
    ctx.font = '400 10.5px "IBM Plex Mono", monospace';
    ctx.fillStyle = ink + '0.55)';
    ctx.textAlign = 'left';
    ctx.fillText(t.scale(fmt(k >= 10 ? Math.round(k) : Math.round(k * 10) / 10)), ox, cssH - 10);
    if (!state.interacted && state.introDone) {
      ctx.textAlign = 'right';
      ctx.fillStyle = '#9C3D1F';
      ctx.fillText(t.hint, cssW - ox, cssH - 10);
    }

    // Readouts
    if (analysed) {
      const ang = Math.round((Math.atan2(state.Fy, state.Fx) * 180 / Math.PI + 360) % 360);
      out.load.textContent = `${(F / 1000).toFixed(2)} kN · ${ang}°`;
      out.smax.textContent = `${Math.round(res.smax)} MPa`;
      out.dmax.textContent = `${res.dmax.toFixed(3)} mm`;
      const rf = res.smax > 0 ? MAT.sy / res.smax : Infinity;
      out.rf.textContent = rf > 99 ? '> 99' : rf.toFixed(2) + (rf < 1 ? ' ✕' : '');
      out.rf.classList.toggle('fail', rf < 1);
    }
  }

  // ---------- Animation loop ----------
  function frame(now) {
    state.raf = 0;
    if (!state.introDone) {
      const t = (now - state.t0) / 1000;
      state.meshP = Math.min(1, t / 0.9);
      state.contourP = Math.max(0, Math.min(1, (t - 1.1) / 0.8));
      const { mesh } = model, lt = TEXT[lang()];
      const nodes = fmt(mesh.nn), els = fmt(mesh.ne);
      logEl.textContent = t < 0.9 ? lt.meshing : t < 1.6 ? lt.meshed(nodes, els) : lt.factor(fmt(model.sys.ndof), model.sys.hb);
      if (t >= 2.1) { state.introDone = true; setDoneLog(); }
    } else if (!state.interacted && !reduceMotion) {
      // Gentle "gust" until the visitor takes over
      const t = now / 1000, a = Math.PI / 2 + 0.12 * Math.sin(t * 0.7), F = F_DEFAULT * (0.85 + 0.15 * Math.sin(t * 1.3));
      state.Fx = F * Math.cos(a); state.Fy = F * Math.sin(a);
    }
    draw();
    if (state.visible && (!state.introDone || (!state.interacted && !reduceMotion))) request();
  }
  function request() { if (!state.raf) state.raf = requestAnimationFrame(frame); }

  function setDoneLog() {
    const ms = Math.max(1, Math.round(model.tMesh + model.tSolve));
    logEl.textContent = TEXT[lang()].done(fmt(model.mesh.nn), fmt(model.mesh.ne), fmt(ms));
  }

  // ---------- Interaction ----------
  function pointerPos(ev) {
    const r = canvas.getBoundingClientRect();
    return [ev.clientX - r.left, ev.clientY - r.top];
  }
  function nearestEdgeNode(x, y) {
    const { mesh } = model, { s, ox, oy } = view;
    let best = Infinity, k = state.k;
    mesh.outlineNodes.forEach((node, j) => {
      if (mesh.fixed[node]) return;
      const d = Math.hypot(ox + mesh.X[node] * s - x, oy - mesh.Y[node] * s - y);
      if (d < best) { best = d; k = j; }
    });
    return k;
  }
  function setForceFromPointer(x, y) {
    const node = model.mesh.outlineNodes[state.k];
    const dx = x - (view.ox + model.mesh.X[node] * view.s), dy = y - (view.oy - model.mesh.Y[node] * view.s);
    let F = Math.hypot(dx, dy) / view.arrowMax * F_MAX;
    if (F < 1) return;
    F = Math.min(F, F_MAX);
    const len = Math.hypot(dx, dy);
    state.Fx = F * dx / len;
    state.Fy = -F * dy / len;
  }
  function takeOver() {
    if (!state.introDone) { state.introDone = true; state.meshP = 1; state.contourP = 1; setDoneLog(); }
    state.interacted = true;
  }

  canvas.addEventListener('pointerdown', ev => {
    if (!model) return;
    takeOver();
    const [x, y] = pointerPos(ev);
    const onHead = view.head && Math.hypot(x - view.head[0], y - view.head[1]) < 18;
    if (!onHead) state.k = nearestEdgeNode(x, y);
    state.dragging = true;
    canvas.setPointerCapture(ev.pointerId);
    setForceFromPointer(x, y);
    draw();
  });
  canvas.addEventListener('pointermove', ev => {
    if (!model || !state.dragging) return;
    const [x, y] = pointerPos(ev);
    setForceFromPointer(x, y);
    draw();
  });
  ['pointerup', 'pointercancel'].forEach(type => canvas.addEventListener(type, () => { state.dragging = false; }));

  canvas.addEventListener('keydown', ev => {
    if (!model) return;
    const F = Math.hypot(state.Fx, state.Fy) || F_DEFAULT;
    let a = Math.atan2(state.Fy, state.Fx), f = F;
    if (ev.key === 'ArrowLeft') a += Math.PI / 18;
    else if (ev.key === 'ArrowRight') a -= Math.PI / 18;
    else if (ev.key === 'ArrowUp') f = Math.min(F_MAX, F + 500);
    else if (ev.key === 'ArrowDown') f = Math.max(0, F - 500);
    else return;
    ev.preventDefault();
    takeOver();
    state.Fx = f * Math.cos(a); state.Fy = f * Math.sin(a);
    draw();
  });

  meshBtn.addEventListener('click', () => {
    state.showMesh = !state.showMesh;
    meshBtn.setAttribute('aria-pressed', state.showMesh);
    draw();
  });
  resetBtn.addEventListener('click', () => {
    if (!model) return;
    state.interacted = false;
    state.k = model.defaultK;
    state.Fx = 0; state.Fy = F_DEFAULT;
    draw();
    request();
  });

  function setStaticText() {
    canvas.setAttribute('aria-label', TEXT[lang()].aria);
    if (model && state.introDone) setDoneLog();
    draw();
  }
  document.addEventListener('portfolio:lang', setStaticText);

  // Legend swatches
  const scaleEl = document.getElementById('feScale');
  BANDS.concat([BEYOND]).forEach((color, i) => {
    const sw = document.createElement('span');
    sw.style.background = color;
    if (i === BANDS.length) sw.className = 'beyond';
    scaleEl.appendChild(sw);
  });
  document.getElementById('feTicks').innerHTML =
    [0, 0.25, 0.5, 0.75, 1].map(f => `<span${f === 1 ? ' class="end"' : ''} style="left:${f * 100 * BANDS.length / (BANDS.length + 1)}%">${Math.round(MAT.sy * f)}</span>`).join('') +
    '<span class="beyond" style="left:100%">&gt;σ<sub>y</sub></span>';

  // Draw the geometry straight away; mesh and solve once the figure is in view
  const tg = performance.now();
  geo = buildGeometry();
  geoTime = performance.now() - tg;
  layout();
  draw();

  // Build when the figure approaches the viewport, play the intro once it is in view
  const io = new IntersectionObserver(entries => {
    entries.forEach(entry => {
      state.visible = entry.isIntersecting;
      if (!entry.isIntersecting) return;
      if (!model) {
        build();
        canvas.setAttribute('aria-label', TEXT[lang()].aria);
        if (reduceMotion) { state.introDone = true; setDoneLog(); }
        else { state.meshP = 0; state.contourP = 0; state.t0 = performance.now(); }
      }
      request();
    });
  }, { threshold: 0.25 });
  io.observe(stage);

  if ('ResizeObserver' in window) {
    new ResizeObserver(() => { if (geo) { layout(); draw(); } }).observe(stage);
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => draw());
})();
