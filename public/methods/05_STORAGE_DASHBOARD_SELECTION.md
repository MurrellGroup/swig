# Result storage, dashboards, local filtering, and repertoire selection

## Browser-local AIRR store

Completed AIRR text is retained as bounded chunks in IndexedDB or referenced by byte ranges in a user-selected output file. A compact record index stores ordinals, identifiers, study fields, primary calls, identities, CDR3 fields, and QC flags. Full AIRR rows are decoded only for requested pages/scans. Browser-resident chunks are gzip-compressed when the platform supports it.

This design is why post-analysis can address a million-row table without keeping a million JavaScript objects in memory. It is also why loading a linked session must rebuild compact indexes from the linked AIRR table: the session deliberately does not duplicate the main table.

## Repertoire dashboard

Facet counts, V/D/J/C or isotype use, CDR3-length distributions, V-identity distributions, and V–J pair counts are accumulated as AIRR batches commit. Opening the dashboard uses those aggregates rather than rescanning full rows. For a comma-separated ambiguous call, `first call` assigns one count to the first label; `fractional` divides one record equally among its listed labels. Plot CSV and SVG contain the exact displayed values.

These are descriptive counts, not posterior allele frequencies and not genomic copy-number estimates.

## Sequence browser

Exact call/facet filters query compact indexes. Substring and combinations that cannot use one index scan only the candidate ordinals. Filters compose with logical AND. For V/D/J/C, `include ambiguous` permits a matching label within a comma-separated call; otherwise a multi-call record is excluded. Opening a record retrieves its one full AIRR row plus sparse Double-D evidence, if present.

## Metadata correction

Applying edited dataset/sample/donor/cohort/timepoint/compartment values updates all rows through an overlay keyed by stable dataset ID and rebuilds affected facets/index values in one IndexedDB transaction. The AIRR sequence and V(D)J calls are unchanged. Downstream state is cleared because its scope partitions may now be stale. Cancellation aborts the transaction and retains the previous overlay, facets, and downstream state.

## Cumulative post-analysis working set

Post-analysis never deletes assigned rows. It keeps a `Uint8Array` mask with one entry per AIRR ordinal. Applying collapse representatives, a chimera threshold, or repertoire selection composes a new mask. Later methods consume only ordinals marked active. `Reset to all records` restores the complete assigned population; it does not discard fitted models or files unless those results depend on an incompatible upstream call policy.

Every post-analysis card has an explicit **Skip step** control. Skipping bypasses that method in guided navigation without fabricating a result or changing the current mask; **Include step** puts it back. Repertoire selection (step 04) starts skipped because its untouched configuration would retain every record and add no information. This is navigation state, is saved with the analysis session, and can be changed before or after earlier steps have run.

## Repertoire selection

The selection block supports dataset/sample/donor/cohort/timepoint/compartment, identifier, locus, V/D1/D2/J/C/isotype calls, CDR3 nucleotide/amino-acid substring, literal/IUPAC/regular-expression motifs, productivity/completeness/frame/stop/D/CDR3 states, Double-D status, CDR3 lengths, V/J identity, and V mutation fraction. Non-empty conditions are ANDed. Comma/newline-separated values within one field use that field's documented any/all behavior.

Call matching accepts exact allele, gene-without-allele, or contained call text. Unless `include ambiguous` is enabled, a multi-call assignment fails a non-empty call filter. V mutation fraction is computed as `1 − v_identity`; a maximum mutation filter also requires a reported nonzero identity.

Motif modes are:

- **literal substring:** regex metacharacters are escaped;
- **IUPAC:** each ambiguity symbol expands to its nucleotide class;
- **regular expression:** the browser's case-insensitive JavaScript regex engine is used.

`Preview count` scans the inherited working set and returns a proposed mask but does not alter downstream state. `Apply` commits that exact preview. Editing a control invalidates the preview so an old count cannot be applied to new settings.

## Literature relationship

The storage, dashboard, filtering, and mask composition are **custom Swig engineering**. AIRR fields and terminology follow the [AIRR Rearrangement schema](https://docs.airr-community.org/en/latest/datarep/rearrangements.html). This is not an AIRR database server or ADC API implementation.

## Complexity and limitations

- Indexed equality/facet filters are near the number of indexed hits; arbitrary substrings/regexes are linear in the selected candidates.
- Dashboard aggregates are exact for committed rows, subject only to the chosen ambiguous-call counting policy.
- JavaScript regular expressions can be expensive for pathological user expressions.
- A working-set mask is selection state, not a new biological data object. Exports make the selected population explicit.
