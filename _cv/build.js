// Builds the English and French CV PDFs from the portfolio site's own content.
// Usage (from the repository root): node _cv/build.js
// Requires Playwright with Chromium: npm i -D playwright && npx playwright install chromium
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const ROOT = path.resolve(__dirname, '..');
const SITE_URL = 'https://rahulsg96-hub.github.io/rahulsubramanian.github.io/';
const OUTPUT = lang => path.join(ROOT, 'files', `Rahul_Subramanian_Girija_CV_${lang.toUpperCase()}.pdf`);

// Runs inside the portfolio page after the language has been applied
function extract() {
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
  const txt = el => (el ? el.textContent.replace(/\s+/g, ' ').trim() : '');
  const plain = markup => { const d = document.createElement('div'); d.innerHTML = markup; return d.textContent.replace(/\s+/g, ' ').trim(); };

  const section = id => {
    const sec = document.getElementById(id);
    return {
      heading: txt($('h2', sec)),
      entries: $$('.entry', sec).map(e => {
        const when = $('.when', e).cloneNode(true);
        const award = $('.award', when);
        if (award) award.remove();
        return {
          when: when.innerHTML.split(/<br\s*\/?>/i).map(plain).filter(Boolean),
          award: txt(award),
          title: txt($('h3', e)),
          org: txt($('.org', e)),
          sub: txt($('.sub', e)),
          bullets: $$('.bullets li:not([data-cv="omit"])', e).map(li => li.innerHTML.trim()),
          note: txt($('.note', e)),
          badge: txt($('.badge', e)),
          authors: txt($('.authors', e)),
          venue: txt($('.venue', e))
        };
      })
    };
  };

  const tel = $('a[href^="tel:"]');
  return {
    lang: document.documentElement.lang,
    name: $('h1').innerText.replace(/\s+/g, ' ').trim(),
    role: txt($('.kicker')),
    tagline: $('.lede').innerHTML.trim(),
    location: txt($('.masthead span:last-child')),
    email: $('a[href^="mailto:"]').getAttribute('href').slice('mailto:'.length),
    phone: txt(tel),
    phoneHref: tel.getAttribute('href'),
    linkedin: $('a[href*="linkedin.com"]').getAttribute('href'),
    facts: $$('.facts > div').map(d => ({ label: txt($('small', d)), value: txt($('span', d)) })),
    languages: {
      heading: txt($('.lang-list h3')),
      items: $$('.lang-item').map(i => ({ name: txt($('.name', i)), level: txt($('.level', i)) }))
    },
    skills: {
      heading: txt($('#skills h2')),
      groups: $$('.filter-btn').filter(b => b.dataset.filter !== 'all').map(b => ({
        key: b.dataset.filter,
        label: txt(b),
        items: $$(`.skill-tag[data-cat="${b.dataset.filter}"]`).map(el => txt(el))
      }))
    },
    experience: section('experience'),
    education: section('education'),
    projects: section('projects'),
    publications: section('publications'),
    certificates: section('certificates')
  };
}

const footer = (name, url) => `
  <div style="width:100%;padding:0 15mm;display:flex;justify-content:space-between;font-family:'IBM Plex Mono',monospace;font-size:6.5px;color:#736D5E;">
    <span>${name} — CV</span>
    <span>${url.replace(/^https?:\/\//, '').replace(/\/$/, '')}</span>
    <span><span class="pageNumber"></span> / <span class="totalPages"></span></span>
  </div>`;

async function build() {
  const browser = await chromium.launch();
  const site = await browser.newPage();
  await site.goto('file://' + path.join(ROOT, 'index.html'), { waitUntil: 'networkidle' });
  const cv = await browser.newPage();

  for (const lang of ['en', 'fr']) {
    await site.evaluate(l => applyLanguage(l), lang);
    const data = await site.evaluate(extract);
    await cv.goto('file://' + path.join(__dirname, 'template.html'), { waitUntil: 'networkidle' });
    await cv.evaluate(([d, url]) => window.renderCV(d, url), [data, SITE_URL]);

    const file = OUTPUT(lang);
    await cv.pdf({
      path: file,
      format: 'A4',
      printBackground: true,
      displayHeaderFooter: true,
      headerTemplate: '<span></span>',
      footerTemplate: footer(data.name, SITE_URL),
      margin: { top: '13mm', bottom: '15mm', left: '15mm', right: '15mm' }
    });
    const pages = (fs.readFileSync(file, 'latin1').match(/\/Type\s*\/Page(?!s)/g) || []).length;
    console.log(`${path.relative(ROOT, file)}: ${pages} page${pages === 1 ? '' : 's'}, ${Math.round(fs.statSync(file).size / 1024)} KB`);
    if (pages !== 2) console.warn(`  warning: expected 2 pages, got ${pages}`);
  }
  await browser.close();
}

if (require.main === module) build();
module.exports = { extract, ROOT, SITE_URL };
