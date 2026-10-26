# CV generator

Builds `files/Rahul_Subramanian_Girija_CV_EN.pdf` and `files/Rahul_Subramanian_Girija_CV_FR.pdf`
from the content of `index.html`, so the CVs always match the site.

```sh
npm i -D playwright && npx playwright install chromium   # once
node _cv/build.js
```

- Edit content on the site (`index.html`, including its EN/FR `translations`), then rebuild.
- `template.html` holds the A4 layout. The script warns if a CV is no longer exactly 2 pages.
- A bullet marked `data-cv="omit"` stays on the site but is left out of the CVs.
- This folder starts with `_`, so GitHub Pages does not publish it.
