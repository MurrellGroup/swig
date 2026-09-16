import { escapeHtml, renderMarkdown } from './guide-markdown.js';
import { findPage, pages, pageUrl } from './guide-pages.js';

const chapter = document.getElementById('chapter');
const main = document.getElementById('content');
const status = document.getElementById('load-status');
const source = document.getElementById('markdown-source');
const pagination = document.getElementById('chapter-pagination');
const chapterLinks = document.getElementById('chapter-links');
const chapters = document.getElementById('chapters');
const requested = new URLSearchParams(location.search).get('page') || 'index';
const page = findPage(requested);

// Each page has an ordinary URL; reload, Back/Forward, and fragment links work
// without a client-side router or a host-specific fallback rule.
let group = '';
for (const item of pages) {
  if (item.group !== group) {
    group = item.group;
    const heading = document.createElement('p');
    heading.className = 'nav-group';
    heading.textContent = group;
    chapterLinks.append(heading);
  }
  const link = document.createElement('a');
  link.href = pageUrl(item.slug);
  link.textContent = item.title;
  if (item.slug === page?.slug) link.setAttribute('aria-current', 'page');
  chapterLinks.append(link);
}

const desktop = matchMedia('(min-width: 900px)');
function updateChapters() { chapters.open = desktop.matches; }
updateChapters();
desktop.addEventListener('change', updateChapters);

function showError(message) {
  status.textContent = message;
  status.setAttribute('role', 'alert');
  main.setAttribute('aria-busy', 'false');
  const help = document.createElement('p');
  const link = document.createElement('a');
  link.href = './';
  link.textContent = 'Return to the guide overview';
  help.append(link);
  chapter.replaceChildren(help);
}

async function loadChapter() {
  if (!page) {
    document.title = 'Chapter not found · Swig immunology';
    showError('This chapter is not in the guide. Choose a chapter from the navigation.');
    return;
  }
  // The manifest is an allow-list, so a query string cannot select an external
  // URL, traverse directories, or load arbitrary files from the application.
  const url = new URL(`${page.slug}.md`, import.meta.url);
  source.href = url.href;
  try {
    const response = await fetch(url, { credentials: 'omit' });
    if (!response.ok) throw new Error(`Chapter request failed (${response.status}).`);
    const text = await response.text();
    if (!/^# /m.test(text)) throw new Error('The server returned an unexpected chapter format.');
    const rendered = renderMarkdown(text);
    // The renderer creates a fixed set of tags, escapes all source text, and
    // rejects executable URLs and raw HTML. It has full-collection tests.
    chapter.innerHTML = rendered.html;
    document.title = `${rendered.headings.find((heading) => heading.level === 1)?.text || page.title} · Swig`;
    status.hidden = true;
    main.setAttribute('aria-busy', 'false');

    const sections = rendered.headings.filter((heading) => heading.level === 2);
    if (sections.length > 1) {
      const contents = document.createElement('details');
      contents.className = 'on-this-page';
      contents.innerHTML = `<summary>On this page</summary><ul>${sections.map((heading) => `<li><a href="#${escapeHtml(heading.id)}">${escapeHtml(heading.text)}</a></li>`).join('')}</ul>`;
      chapter.querySelector('h1')?.after(contents);
    }

    const position = pages.indexOf(page);
    for (const [offset, label] of [[-1, 'Previous'], [1, 'Next']]) {
      const neighbour = pages[position + offset];
      if (!neighbour) continue;
      const link = document.createElement('a');
      link.href = pageUrl(neighbour.slug);
      link.rel = offset < 0 ? 'prev' : 'next';
      const small = document.createElement('span');
      small.textContent = label;
      link.append(small, document.createTextNode(neighbour.title));
      pagination.append(link);
    }

    // Fragment targets only exist after the Markdown has loaded. Decoding is
    // guarded because a malformed fragment must not prevent reading a chapter.
    if (location.hash) {
      let id;
      try { id = decodeURIComponent(location.hash.slice(1)); } catch { id = location.hash.slice(1); }
      const target = document.getElementById(id);
      if (target) {
        target.tabIndex = -1;
        target.focus({ preventScroll: true });
        target.scrollIntoView();
      }
    }
  } catch (error) {
    console.error(error);
    showError('The chapter could not be loaded. Reload the page, use the Markdown source link, or choose another chapter.');
  }
}

loadChapter();
