import { useMemo, useState } from "react";

import { referenceFastaNames } from "./reference-fasta";

interface Props {
  fasta: string;
  excluded: string[];
  onChange: (excluded: string[]) => void;
  label?: string;
}

const DISPLAY_LIMIT = 300;

export function ReferenceAlleleExclusionEditor({ fasta, excluded, onChange, label = "Advanced reference-alignment options" }: Props) {
  const [search, setSearch] = useState("");
  const names = useMemo(() => referenceFastaNames(fasta), [fasta]);
  const excludedSet = useMemo(() => new Set(excluded), [excluded]);
  const normalizedSearch = search.trim().toLocaleUpperCase();
  const matches = useMemo(() => names.filter((name) => !normalizedSearch || name.toLocaleUpperCase().includes(normalizedSearch)), [names, normalizedSearch]);
  const visible = matches.slice(0, DISPLAY_LIMIT);
  const activeCount = names.reduce((count, name) => count + (excludedSet.has(name) ? 1 : 0), 0);

  const commit = (next: Set<string>) => onChange([...next].sort((left, right) => left.localeCompare(right, undefined, { numeric: true })));
  const toggle = (name: string, checked: boolean) => {
    const next = new Set(excludedSet);
    if (checked) next.add(name); else next.delete(name);
    commit(next);
  };
  const setMatches = (checked: boolean) => {
    const next = new Set(excludedSet);
    for (const name of matches) if (checked) next.add(name); else next.delete(name);
    commit(next);
  };

  return <details className="post-advanced reference-exclusion-editor">
    <summary>{label}{activeCount ? ` · ${activeCount.toLocaleString()} excluded` : ""}</summary>
    <div>
      <p>Exclusions use exact FASTA identifiers and are applied before CHMMAIRRa validates or builds the MSA. An assigned allele removed here will be reported as a missing reference rather than silently replaced.</p>
      <div className="reference-exclusion-toolbar">
        <label><span>Find allele</span><input type="search" value={search} onChange={(event) => setSearch(event.target.value)} placeholder="IGHV1-2*01…" /></label>
        <div className="result-actions"><button type="button" disabled={!matches.length} onClick={() => setMatches(true)}>Exclude {normalizedSearch ? "matches" : "all"}</button><button type="button" disabled={!matches.some((name) => excludedSet.has(name))} onClick={() => setMatches(false)}>Include {normalizedSearch ? "matches" : "all"}</button></div>
      </div>
      <div className="reference-exclusion-summary"><strong>{activeCount.toLocaleString()} excluded</strong><span>{Math.max(0, names.length - activeCount).toLocaleString()} of {names.length.toLocaleString()} reference records retained</span></div>
      {visible.length ? <div className="reference-exclusion-list">{visible.map((name) => <label key={name}><input type="checkbox" checked={excludedSet.has(name)} onChange={(event) => toggle(name, event.target.checked)} /><code>{name}</code><span>{excludedSet.has(name) ? "excluded" : "included"}</span></label>)}</div> : <div className="method-placeholder small"><span>∅</span><h4>{names.length ? "No allele matches this search" : "No reference FASTA is available"}</h4></div>}
      {matches.length > DISPLAY_LIMIT && <small>{DISPLAY_LIMIT.toLocaleString()} of {matches.length.toLocaleString()} matching identifiers are shown. Bulk include/exclude applies to every match.</small>}
    </div>
  </details>;
}
