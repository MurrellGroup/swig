import { useEffect, useMemo, useRef, useState } from "react";

export interface FacetPickerItem {
  value: string;
  count?: number;
}

function selectedValues(value: string): string[] {
  return value.split(/[\n,;]+/).map((item) => item.trim()).filter(Boolean);
}

export function FacetPicker({
  label,
  value,
  items,
  onChange,
  multiple = false,
  placeholder = "Any",
  className = "",
  allowCustom = false,
  help,
}: {
  label: string;
  value: string;
  items: FacetPickerItem[];
  onChange: (value: string) => void;
  multiple?: boolean;
  placeholder?: string;
  className?: string;
  allowCustom?: boolean;
  help?: string;
}) {
  const [query, setQuery] = useState("");
  const detailsRef = useRef<HTMLDetailsElement>(null);
  const selected = useMemo(() => selectedValues(value), [value]);
  const selectedSet = useMemo(() => new Set(selected), [selected]);
  const options = useMemo(() => {
    const normalized = query.trim().toLocaleLowerCase();
    return items
      .filter((item) => item.value && (!normalized || item.value.toLocaleLowerCase().includes(normalized)))
      .slice(0, 300);
  }, [items, query]);
  const summary = selected.length === 0
    ? placeholder
    : selected.length === 1
      ? selected[0]
      : `${selected.length.toLocaleString()} selected`;
  const explanatoryText = help ?? `Choose from ${label.toLocaleLowerCase()} values observed in this analysis. Search only narrows this menu; it does not run an analysis.`;

  useEffect(() => {
    const closeOnOutsidePointer = (event: PointerEvent) => {
      const details = detailsRef.current;
      if (details?.open && !details.contains(event.target as Node)) details.open = false;
    };
    const closeOnEscape = (event: KeyboardEvent) => {
      if (event.key === "Escape" && detailsRef.current?.open) {
        detailsRef.current.open = false;
        detailsRef.current.querySelector("summary")?.focus();
      }
    };
    document.addEventListener("pointerdown", closeOnOutsidePointer);
    document.addEventListener("keydown", closeOnEscape);
    return () => {
      document.removeEventListener("pointerdown", closeOnOutsidePointer);
      document.removeEventListener("keydown", closeOnEscape);
    };
  }, []);

  function choose(next: string) {
    if (!multiple) {
      onChange(next);
      if (detailsRef.current) detailsRef.current.open = false;
      setQuery("");
      return;
    }
    const updated = selectedSet.has(next)
      ? selected.filter((item) => item !== next)
      : [...selected, next];
    onChange(updated.join(", "));
  }

  const normalizedQuery = query.trim();
  const hasExactQuery = items.some((item) => item.value.toLocaleLowerCase() === normalizedQuery.toLocaleLowerCase());

  return <div className={`facet-picker ${selected.length ? "has-value" : ""} ${className}`.trim()} title={explanatoryText} data-field-help={explanatoryText}>
    <span className="facet-picker-label">{label}</span>
    <details ref={detailsRef} onToggle={(event) => { if (!(event.currentTarget as HTMLDetailsElement).open) setQuery(""); }}>
      <summary aria-label={`${label}: ${summary}`} title={selected.length ? `${explanatoryText}\nSelected: ${selected.join(", ")}` : explanatoryText}><span>{summary}</span><i aria-hidden="true">⌄</i></summary>
      <div className="facet-picker-popover">
        {(items.length > 8 || allowCustom) && <label className="facet-picker-search"><span className="visually-hidden">Search {label}</span><input type="search" value={query} autoComplete="off" placeholder={`Search ${label.toLocaleLowerCase()}…`} onChange={(event) => setQuery(event.target.value)} onKeyDown={(event) => { if (event.key === "Enter" && allowCustom && normalizedQuery) { event.preventDefault(); choose(normalizedQuery); } }} /></label>}
        <div className="facet-picker-options" role="listbox" aria-multiselectable={multiple || undefined}>
          <button type="button" className={!selected.length ? "selected" : ""} onClick={() => { onChange(""); if (detailsRef.current) detailsRef.current.open = false; setQuery(""); }}><span>{placeholder}</span>{!selected.length && <b>✓</b>}</button>
          {allowCustom && normalizedQuery && !hasExactQuery && <button type="button" className="facet-custom-value" onClick={() => choose(normalizedQuery)}><span>Use “{normalizedQuery}”</span><small>typed value</small></button>}
          {options.map((item) => <button type="button" role="option" aria-selected={selectedSet.has(item.value)} className={selectedSet.has(item.value) ? "selected" : ""} key={item.value} onClick={() => choose(item.value)}><span title={item.value}>{item.value}</span>{typeof item.count === "number" && <small>{item.count.toLocaleString()}</small>}{selectedSet.has(item.value) && <b>✓</b>}</button>)}
          {!options.length && <p>No matching values.</p>}
        </div>
        {multiple && selected.length > 0 && <footer><span>{selected.length.toLocaleString()} selected</span><button type="button" onClick={() => { if (detailsRef.current) detailsRef.current.open = false; }}>Done</button></footer>}
      </div>
    </details>
  </div>;
}

export function uniqueFacetItems(values: string[]): FacetPickerItem[] {
  const counts = new Map<string, number>();
  for (const raw of values) {
    const value = raw.trim();
    if (value) counts.set(value, (counts.get(value) ?? 0) + 1);
  }
  return [...counts.entries()]
    .map(([value, count]) => ({ value, count }))
    .sort((left, right) => left.value.localeCompare(right.value, undefined, { numeric: true }));
}
