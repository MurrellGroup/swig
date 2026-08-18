# Swig 0.37.1

## One assignment default everywhere

- RIAT-MP is now the default for hand-written/normalized pipeline configs and for direct `swig-cli --vdj` runs, matching all browser workflows.
- `--assigner aer` and `--assigner standard` remain explicit overrides, as do corresponding config values.
- Direct V(D)J startup logging reports the selected assignment strategy.
- Restoring an older browser session that has no recorded assignment strategy now also falls back to RIAT-MP.

This changes only selection of the default strategy. The AER, RIAT-MP, and standard implementations and their outputs when explicitly selected are unchanged from 0.37.0.
