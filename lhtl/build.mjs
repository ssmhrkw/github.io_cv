// Assemble self-contained single-file pages from src/:
//   dist/lh-tl-listening.html  both tools with a tab bar
//   dist/lh-listening.html     floor impact (LH) only
//   dist/tl-listening.html     airborne insulation (Dr / Rr / TL) only
// All three share the same code and the same saved settings.
import fs from 'node:fs';
const VERSION = '2.0.0';
const r = f => fs.readFileSync(new URL('./src/' + f, import.meta.url), 'utf8');
const parts = { '/*__CSS__*/': r('style.css'), '/*__DSP__*/': r('dsp.js'), '/*__SHELL__*/': r('shell.js'), '/*__LH__*/': r('lh.js'), '/*__TL__*/': r('tl.js') };
if (/<\/script>/i.test(Object.values(parts).join(''))) throw new Error('a source part contains </script>');
const variants = {
  both: { file: 'lh-tl-listening.html', title: '床衝撃音・遮音の聴き比べ（LH／Dr・Rr・TL）', h1: '床衝撃音・遮音の聴き比べ' },
  lh: { file: 'lh-listening.html', title: '床衝撃音の聴き比べ（LH）', h1: '床衝撃音の聴き比べ（LH）' },
  tl: { file: 'tl-listening.html', title: '遮音の聴き比べ（Dr／Rr・TL）', h1: '壁の遮音の聴き比べ（Dr／Rr・TL）' }
};
fs.mkdirSync(new URL('./dist/', import.meta.url), { recursive: true });
for (const [v, o] of Object.entries(variants)){
  let html = r('index.template.html').replace(/<!--IF:(\w+)-->([\s\S]*?)<!--\/IF:\1-->/g, (m, tag, body) => (v === 'both' || tag === v) ? body : '');
  for (const [k, text] of Object.entries(parts)) html = html.replace(k, () => text);
  html = html.split('__VERSION__').join(VERSION).replace('__TITLE__', o.title).replace('__H1__', o.h1);
  if (v !== 'both') html = html.replace(/ role="tabpanel" aria-labelledby="tab-\w+"/g, ' role="region" aria-label="' + o.h1 + '"');
  if (/__\w+__|\/\*__|<!--\/?IF/.test(html)) throw new Error('unresolved placeholder in ' + o.file);
  fs.writeFileSync(new URL('./dist/' + o.file, import.meta.url), html);
  console.log('built dist/' + o.file, (html.length / 1024).toFixed(1) + ' kB', 'v' + VERSION);
}
