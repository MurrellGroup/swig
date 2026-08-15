# Swig 0.25.2

This patch removes the redundant full-table metadata rewrite from portable-session and project restoration.

- Saved sample, donor, cohort, timepoint, and compartment values are now installed before the linked AIRR table is indexed.
- Each record receives its final searchable metadata during the existing import write instead of through a second IndexedDB cursor update.
- Linked-file fingerprinting still hashes the original AIRR bytes, unchanged by saved metadata overlays.
- Filtering, record details, downstream scans, and AIRR exports retain the restored study metadata.
- V(D)J assignment and all downstream scientific methods are unchanged.
