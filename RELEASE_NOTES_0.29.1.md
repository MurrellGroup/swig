# Swig 0.29.1 release notes

## Standalone CLI packaging fix

- The compiled Bun executable now opens its embedded `swig-worker.js` entry point instead of looking for an unbundled `swig-worker.mjs` file.
- The native CLI smoke test now completes annotation and lazy lineage-study export on the compiled macOS Apple-silicon executable.

No analysis method, default, or output schema changed from 0.29.0.
