/**
 * Renderer for the guide's deliberately small Markdown vocabulary: headings,
 * paragraphs, links/images, emphasis, code, lists, rules, and pipe tables.
 * Raw HTML is escaped. This is not used for uploaded sequences or arbitrary URLs.
 * Keep the full-collection regression tests in sync when adding Markdown syntax.
 */
export function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  })[char]);
}

export function headingId(text) {
  return text.toLowerCase().replace(/[^\p{L}\p{M}\p{N}_\-\s]/gu, '').replace(/\s/g, '-');
}

export function safeUrl(value, image = false) {
  if (/[\u0000-\u0020\u007f\\]/u.test(value) || value.startsWith('//')) return null;
  if (/^[a-z][a-z\d+.-]*:/i.test(value)) {
    return !image && /^(https?:|mailto:)/i.test(value) ? value : null;
  }
  if (image && !/^(?:\.\/)?assets\/[a-z\d-]+\.svg$/i.test(value)) return null;
  return value;
}

export function readerUrl(value) {
  const match = /^(?:\.\/)?([a-z\d_-]+)\.md(#[^\s]*)?$/i.exec(value);
  return match ? `?page=${encodeURIComponent(match[1])}${match[2] || ''}` : value;
}

function closeBracket(text, start, open, close) {
  let depth = 1;
  for (let i = start + 1; i < text.length; i += 1) {
    if (text[i] === '\\') { i += 1; continue; }
    if (text[i] === open) depth += 1;
    if (text[i] === close && --depth === 0) return i;
  }
  return -1;
}

/** Render inline syntax while retaining plain text for accessible headings. */
function inline(text, state) {
  let html = '';
  let plain = '';
  for (let i = 0; i < text.length;) {
    if (text[i] === '\\' && /[!"#$%&'()*+,\-./:;<=>?@[\\\]^_`{|}~]/.test(text[i + 1] || '')) {
      html += escapeHtml(text[i + 1]); plain += text[i + 1]; i += 2; continue;
    }
    if (text[i] === '`') {
      const fence = /^`+/.exec(text.slice(i))[0];
      const end = text.indexOf(fence, i + fence.length);
      if (end !== -1) {
        let code = text.slice(i + fence.length, end).replace(/\n/g, ' ');
        if (/^ .* $/.test(code) && /[^ ]/.test(code)) code = code.slice(1, -1);
        html += `<code>${escapeHtml(code)}</code>`; plain += code;
        i = end + fence.length; continue;
      }
    }
    const image = text.startsWith('![', i);
    if (image || text[i] === '[') {
      const start = i + (image ? 1 : 0);
      const end = closeBracket(text, start, '[', ']');
      const targetEnd = end !== -1 && text[end + 1] === '(' ? closeBracket(text, end + 1, '(', ')') : -1;
      if (targetEnd !== -1) {
        const label = inline(text.slice(start + 1, end), state);
        const original = text.slice(end + 2, targetEnd);
        const url = safeUrl(original, image);
        plain += label.plain;
        if (url === null) {
          html += label.html;
        } else if (image) {
          state.images.push({ url, alt: label.plain });
          html += `<a class="figure-link" href="${escapeHtml(url)}" target="_blank" rel="noopener noreferrer" aria-label="${escapeHtml(`Open figure at full size: ${label.plain} (new tab)`)}"><img src="${escapeHtml(url)}" alt="${escapeHtml(label.plain)}" decoding="async"></a>`;
        } else {
          state.links.push(url);
          html += `<a href="${escapeHtml(readerUrl(url))}">${label.html}</a>`;
        }
        i = targetEnd + 1; continue;
      }
    }
    const marker = text.startsWith('**', i) ? '**' : text.startsWith('__', i) ? '__' : /[*_]/.test(text[i]) ? text[i] : null;
    if (marker && !/\s/.test(text[i + marker.length] || ' ') && !(marker.includes('_') && /[\p{L}\p{N}]/u.test(text[i - 1] || ''))) {
      const end = text.indexOf(marker, i + marker.length);
      if (end > i + marker.length && !/\s/.test(text[end - 1])) {
        const inner = inline(text.slice(i + marker.length, end), state);
        const tag = marker.length === 2 ? 'strong' : 'em';
        html += `<${tag}>${inner.html}</${tag}>`; plain += inner.plain;
        i = end + marker.length; continue;
      }
    }
    html += escapeHtml(text[i]); plain += text[i]; i += 1;
  }
  return { html, plain };
}

function tableCells(line) {
  const value = line.trim().replace(/^\|/, '').replace(/(?<!\\)\|$/, '');
  const cells = [];
  let cell = '';
  let codeFence = 0;
  for (let i = 0; i < value.length; i += 1) {
    if (value[i] === '\\' && value[i + 1] === '|') { cell += '\\|'; i += 1; continue; }
    if (value[i] === '`') {
      const length = /^`+/.exec(value.slice(i))[0].length;
      if (!codeFence) codeFence = length;
      else if (codeFence === length) codeFence = 0;
      cell += '`'.repeat(length); i += length - 1; continue;
    }
    if (value[i] === '|' && !codeFence) { cells.push(cell.trim()); cell = ''; }
    else cell += value[i];
  }
  cells.push(cell.trim());
  return cells;
}

function isTableRule(line) {
  const cells = tableCells(line || '');
  return cells.length > 1 && cells.every((cell) => /^:?-{3,}:?$/.test(cell));
}

export function renderMarkdown(markdown) {
  const state = { links: [], images: [], headings: [] };
  const seenIds = new Set();
  const lines = markdown.replace(/\r\n?/g, '\n').split('\n');
  const blocks = [];
  const isRule = (line) => /^ {0,3}(?:(?:-\s*){3,}|(?:\*\s*){3,}|(?:_\s*){3,})$/.test(line);
  const blockStart = (index) => /^ {0,3}(?:#{1,6} |```|~~~|[-+*] |\d+\. )/.test(lines[index] || '') || isRule(lines[index] || '') || isTableRule(lines[index + 1]);
  const render = (text) => inline(text, state).html;

  for (let i = 0; i < lines.length;) {
    if (!lines[i].trim()) { i += 1; continue; }
    const fence = /^ {0,3}(`{3,}|~{3,})([^`]*)$/.exec(lines[i]);
    if (fence) {
      const code = [];
      const end = new RegExp(`^ {0,3}${fence[1][0]}{${fence[1].length},}\\s*$`);
      i += 1;
      while (i < lines.length && !end.test(lines[i])) code.push(lines[i++]);
      i += 1;
      blocks.push(`<pre><code>${escapeHtml(code.join('\n'))}\n</code></pre>`);
      continue;
    }
    const heading = /^ {0,3}(#{1,6})\s+(.+?)\s*#*\s*$/.exec(lines[i]);
    if (heading) {
      const content = inline(heading[2], state);
      const base = headingId(content.plain);
      let id = base;
      for (let suffix = 1; seenIds.has(id); suffix += 1) id = `${base}-${suffix}`;
      seenIds.add(id);
      const level = heading[1].length;
      state.headings.push({ level, id, text: content.plain });
      blocks.push(`<h${level} id="${escapeHtml(id)}">${content.html}</h${level}>`);
      i += 1; continue;
    }
    if (isRule(lines[i])) { blocks.push('<hr>'); i += 1; continue; }
    if (isTableRule(lines[i + 1])) {
      const headers = tableCells(lines[i]);
      const alignments = tableCells(lines[i + 1]).map((cell) => cell.endsWith(':') ? (cell.startsWith(':') ? 'center' : 'right') : 'left');
      const cell = (text, index, tag) => `<${tag}${tag === 'th' ? ' scope="col"' : ''} class="align-${alignments[index] || 'left'}">${render(text)}</${tag}>`;
      const rows = [`<thead><tr>${headers.map((text, index) => cell(text, index, 'th')).join('')}</tr></thead><tbody>`];
      i += 2;
      while (i < lines.length && lines[i].trim() && lines[i].includes('|')) {
        const cells = tableCells(lines[i++]);
        rows.push(`<tr>${headers.map((_, index) => cell(cells[index] || '', index, 'td')).join('')}</tr>`);
      }
      blocks.push(`<div class="table-scroll" role="region" aria-label="Table (scroll horizontally on small screens)" tabindex="0"><table>${rows.join('')}</tbody></table></div>`);
      continue;
    }
    const list = /^ {0,3}([-+*]|\d+\.)\s+(.+)$/.exec(lines[i]);
    if (list) {
      const ordered = /\d/.test(list[1]);
      const items = [];
      const pattern = ordered ? /^ {0,3}\d+\.\s+(.+)$/ : /^ {0,3}[-+*]\s+(.+)$/;
      while (i < lines.length) {
        const item = pattern.exec(lines[i]);
        if (!item) break;
        items.push(`<li>${render(item[1])}</li>`); i += 1;
      }
      const tag = ordered ? 'ol' : 'ul';
      const start = ordered && parseInt(list[1], 10) !== 1 ? ` start="${parseInt(list[1], 10)}"` : '';
      blocks.push(`<${tag}${start}>${items.join('')}</${tag}>`); continue;
    }
    const paragraph = [lines[i++]];
    while (i < lines.length && lines[i].trim() && !blockStart(i)) paragraph.push(lines[i++]);
    blocks.push(`<p>${render(paragraph.join('\n'))}</p>`);
  }
  return { html: blocks.join('\n'), ...state };
}
