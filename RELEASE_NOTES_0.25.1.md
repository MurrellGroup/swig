# Swig 0.25.1

This patch fixes UCA settings that reverted immediately after a completed inference.

- A saved UCA result snapshot is now adopted only once. Clearing that result by editing a setting no longer causes the same stale snapshot to restore its old options.
- Every UCA setting, including **Additional-D probability**, can be changed after a run and used for **change and rerun**.
- Legitimate saved-session results and results belonging to a newly selected lineage/alignment are still restored.
- Additional-D probability continues to accept the full range from 0 through 1; its default remains 0.015.
