# Swig 0.16.1

This release replaces the previous one-pass missing-V heuristic with a stricter two-pass, lineage-aware diagnostic.

## Missing-V evidence audit and correction

- Candidate discovery uses one deterministic lowest-SHM representative per assigned lineage. A clonal expansion can contribute only one discovery observation.
- Multi-site candidates require the linked alternate bases to co-occur on the same molecule and to recur in nearly the same set of independent lineages. Linked supersets suppress redundant singleton/subset warnings.
- Validation makes a second streaming pass over every retained lineage member. A lineage is vetoed if any covered member contains the parent-reference nucleotide at any candidate site. Third nucleotide states and inconsistent V assignments are separate vetoes.
- Supporting lineages must span configurable numbers of distinct CDR3 nucleotide sequences, CDR3 lengths, and J calls. A configurable number must be near germline after candidate sites are removed from their mutation count. Third nucleotide states across covered lineages must remain below a configurable fraction (2% by default).
- V coordinates now honor AIRR `v_germline_start`. Periods are alignment gaps rather than counted nucleotide positions.
- AID WRCY/RGYW context is marked explicitly. Hotspots receive a larger SHM null probability; the reported tail remains a screening/ranking statistic, not a posterior probability.
- Exact sequences already present in the composed V reference are suppressed.
- The results panel reports proposed patterns, reference-state vetoes, independent support, joint coverage, CDR3/J diversity, and near-germline evidence. CSV, TSV, JSON Lines, FASTA, and SVG exports include the expanded audit trail.

## Compatibility and scale

- Interactive and pipeline modes both run the same two streaming passes over the current committed working set.
- Candidate validation stores compact per-lineage bit flags. A 64-pattern-per-parent-V guard prevents an adversarial hotspot-rich repertoire from exhausting browser memory; any omitted preliminary patterns are reported.
- A synthetic 50,000-lineage, 220-nt two-pass benchmark completed in 3.55 seconds in the release environment and retained about 17.3 MiB of additional V8 heap after forced garbage collection. This is an engineering stress check, not a browser- or hardware-independent throughput claim.
- Sessions retain the new options and results. A saved dashboard produced by the removed one-pass algorithm is deliberately not restored; rerun the screen after loading that session.

## Regression coverage

Tests now verify absolute V coordinates through a truncated alignment with a pure-gap column, complete multi-SNP linkage, reference-state vetoes from non-representative lineage members, and resistance to support inflation by 100 descendants of one lineage.
