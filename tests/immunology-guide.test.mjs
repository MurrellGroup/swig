import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, resolve, relative } from 'node:path';
import test from 'node:test';

import { renderMarkdown, readerUrl, safeUrl } from '../public/immunology/guide-markdown.js';
import { pages, findPage, pageUrl } from '../public/immunology/guide-pages.js';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../public/immunology');
const documents = new Map(readdirSync(root).filter((name) => name.endsWith('.md')).map((name) => [
  name, renderMarkdown(readFileSync(resolve(root, name), 'utf8')),
]));

test('every reader-facing Markdown file is in the navigation exactly once', () => {
  assert.equal(new Set(pages.map((page) => page.slug)).size, pages.length);
  assert.deepEqual(pages.map((page) => `${page.slug}.md`).sort(), [...documents.keys()].filter((name) => name !== 'MAINTAINERS.md').sort());
  assert.equal(pages.length, 16);
  assert.equal(findPage('../README'), null);
  assert.equal(findPage('https://example.com'), null);
  assert.equal(findPage('missing'), null);
  assert.equal(findPage().slug, 'index');
});

for (const [name, document] of documents) {
  test(`${name}: headings, links, anchors, and figures`, () => {
    assert.equal(document.headings.filter((heading) => heading.level === 1).length, 1);
    assert.equal(new Set(document.headings.map((heading) => heading.id)).size, document.headings.length);
    for (const link of document.links) {
      assert.ok(safeUrl(link), `${name}: unsafe link ${link}`);
      if (/^(https?:|mailto:)/.test(link)) continue;
      const [path, fragment] = link.split('#');
      const target = path ? path.replace(/^\.\//, '') : name;
      assert.ok(documents.has(target), `${name}: missing page ${link}`);
      if (target !== 'MAINTAINERS.md') assert.ok(findPage(target.replace(/\.md$/, '')), `${name}: unlisted reader page ${target}`);
      if (fragment) assert.ok(documents.get(target).headings.some((heading) => heading.id === decodeURIComponent(fragment)), `${name}: missing anchor ${link}`);
    }
    for (const image of document.images) {
      const path = resolve(root, image.url);
      assert.ok(!relative(root, path).startsWith('..'));
      assert.ok(existsSync(path), `${name}: missing figure ${image.url}`);
      assert.ok(image.alt.length > 20, `${name}: missing descriptive alt text`);
      const svg = readFileSync(path, 'utf8');
      assert.match(svg, /<svg[^>]*\bviewBox=/);
      assert.match(svg, /<title[ >]/);
      assert.match(svg, /<desc[ >]/);
      assert.doesNotMatch(svg, /<script\b|<foreignObject\b|\bon\w+\s*=/i);
    }
  });
}

test('collection includes all four SVGs and the expected table-bearing chapters', () => {
  assert.equal([...documents.values()].flatMap((document) => document.images).length, 4);
  for (const file of ['02-vdj-recombination.md', '06-vdj-assignment.md', '13-analysis-guide.md']) {
    assert.match(documents.get(file).html, /<table>/);
  }
});

test('inline concepts, code, emphasis, and parenthesized Wikipedia targets', () => {
  const document = renderMarkdown('# Receptors\n\n[**V(D)J**](https://en.wikipedia.org/wiki/V(D)J_recombination) uses `IGHV3-23*01`; *mutation* and **selection**.');
  assert.match(document.html, /<a href="https:\/\/en.wikipedia.org\/wiki\/V\(D\)J_recombination"><strong>V\(D\)J<\/strong><\/a>/);
  assert.match(document.html, /<code>IGHV3-23\*01<\/code>/);
  assert.match(document.html, /<em>mutation<\/em> and <strong>selection<\/strong>/);
});

test('GitHub-style heading IDs and duplicate suffixes', () => {
  const document = renderMarkdown('# B- and T-cell immunology\n\n## UCA: `V(D)J`\n\n## UCA: `V(D)J`\n\n## αβ T cells');
  assert.deepEqual(document.headings.map((heading) => heading.id), ['b--and-t-cell-immunology', 'uca-vdj', 'uca-vdj-1', 'αβ-t-cells']);
});

test('fenced sequences, rules, and lists are separate blocks', () => {
  const document = renderMarkdown('# Example\n\n```text\nA < T & G\n```\n\n---\n\n- One\n- Two\n\n3. Three\n4. Four');
  assert.match(document.html, /<pre><code>A &lt; T &amp; G\n<\/code><\/pre>/);
  assert.match(document.html, /<hr>/);
  assert.match(document.html, /<ul><li>One<\/li><li>Two<\/li><\/ul>/);
  assert.match(document.html, /<ol start="3"><li>Three<\/li><li>Four<\/li><\/ol>/);
});

test('tables retain inline markup, escaped pipes, alignment, and column headers', () => {
  const document = renderMarkdown('| Feature | Meaning |\n| :--- | ---: |\n| **V** | [Gene](glossary.md#gene) |\n| a\\|b | `x|y` |');
  assert.match(document.html, /<th scope="col" class="align-left">Feature<\/th>/);
  assert.match(document.html, /<th scope="col" class="align-right">Meaning<\/th>/);
  assert.match(document.html, /<strong>V<\/strong>/);
  assert.match(document.html, /href="\?page=glossary#gene"/);
  assert.match(document.html, /a\|b<\/td>/);
  assert.match(document.html, /<code>x\|y<\/code>/);
});

test('internal links retain fragments and work under a project base path', () => {
  const base = 'https://example.org/swig/immunology/';
  assert.equal(new URL(readerUrl('11-unmutated-common-ancestors.md#uncertainty'), base).href, `${base}?page=11-unmutated-common-ancestors#uncertainty`);
  assert.equal(readerUrl('#local'), '#local');
  assert.equal(readerUrl('https://example.org/file.md'), 'https://example.org/file.md');
  assert.equal(new URL(pageUrl('index'), base).pathname, '/swig/immunology/');
  assert.equal(new URL('../', base).pathname, '/swig/');
  assert.equal(new URL(pageUrl('glossary'), base).pathname, '/swig/immunology/');
});

test('raw HTML, executable URLs, and external image loads are rejected', () => {
  const document = renderMarkdown('<script>alert(1)</script>\n\n[bad](javascript:alert(1)) [data](data:text/html,x) [bad](//example.com)\n\n![remote](https://example.com/image.svg)');
  assert.match(document.html, /&lt;script&gt;/);
  assert.doesNotMatch(document.html, /<script|href=|<img/);
  for (const url of ['javascript:alert(1)', 'data:text/html,test', 'java\nscript:alert(1)', '//example.org', '\\example.org']) assert.equal(safeUrl(url), null);
  assert.equal(safeUrl('../image.svg', true), null);
  assert.equal(safeUrl('assets/01-vdj-recombination.svg', true), 'assets/01-vdj-recombination.svg');
});

test('image links preserve alt text and open the original SVG safely', () => {
  const document = renderMarkdown('![V, D, and J segments](assets/01-vdj-recombination.svg)');
  assert.match(document.html, /target="_blank" rel="noopener noreferrer"/);
  assert.match(document.html, /alt="V, D, and J segments"/);
  assert.match(document.html, /href="assets\/01-vdj-recombination.svg"/);
});

test('reader resources are self-contained and do not import analysis code', () => {
  const html = readFileSync(resolve(root, 'index.html'), 'utf8');
  assert.match(html, /Content-Security-Policy/);
  for (const name of ['guide.css', 'guide.js', 'guide-markdown.js', 'guide-pages.js']) assert.ok(existsSync(resolve(root, name)));
  assert.doesNotMatch(html, /(?:src|href)="https?:\/\/[^\"]+\.(?:js|css)/);
  const script = readFileSync(resolve(root, 'guide.js'), 'utf8');
  assert.doesNotMatch(script, /indexedDB|localStorage|swig-app|swiftig-runtime/);
});

test('application entry points preserve the analysis tab and use Vite BASE_URL', () => {
  const entry = readFileSync(resolve(root, '../../src/main.tsx'), 'utf8');
  assert.match(entry, /import\.meta\.env\.BASE_URL/);
  assert.equal((entry.match(/target="_blank" rel="noopener noreferrer"/g) || []).length, 2);
  assert.match(entry, /page=13-analysis-guide/);
  assert.match(entry, /<SwigApp \/>/);
});

test('README exposes the guide near the top with working repository links', () => {
  const readme = readFileSync(new URL('../README.md', import.meta.url), 'utf8');
  const intro = readme.slice(0, readme.indexOf('## Deploy on GitHub Pages'));
  for (const file of ['index.md', '13-analysis-guide.md', 'glossary.md', 'references.md', 'MAINTAINERS.md']) {
    assert.ok(intro.includes(`public/immunology/${file}`));
    assert.ok(existsSync(new URL(`../public/immunology/${file}`, import.meta.url)));
  }
});
