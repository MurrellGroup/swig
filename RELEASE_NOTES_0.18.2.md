# Swig 0.18.2

## Double-D lineage germline reconstruction

- A supported Double-D lineage root is now constructed as V–D1–D2–J. Swig removes the unchanged baseline single-D junction template, projects both sparse D alignments through the combined AIRR query coordinates, and retains NP1, NP2, and NP3 as unresolved `N` states.
- Closest-member reconstruction ranks only safely projected VDDJ-aware members when a lineage contains usable Double-D evidence. This prevents a slightly cleaner V/J match with only the baseline single-D composite from silently dropping D2.
- Consensus reconstruction likewise votes only across safely projected VDDJ-aware members for a Double-D lineage. Ordinary VDJ lineages follow the previous path unchanged.
- The lineage root panel explicitly reports V–D1–D2–J construction and the selected D1/D2 calls. An incomplete imported D2 sidecar is reported as unresolved and is never partially projected.
- Regression tests cover the D1/D2 root sequence, N-masked insertion intervals, closest-member selection, consensus isolation, quick-view/FastTree outgroup input, incomplete sidecars, and unchanged ordinary VDJ reconstruction.
