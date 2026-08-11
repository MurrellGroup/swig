# SwiftIG WebAssembly source provenance

The algorithm sources in this directory are derived from MurrellGroup/swiftig commit
`21304be2499f8dce330fc8890f0d891311ef2012` (SwiftIG 0.2.0).

Browser-specific differences are deliberately narrow:

- `src/web_adapter.cpp` supplies the allocation, database-initialization, annotation,
  result, and error ABI exported to JavaScript.
- filesystem/index serialization code is excluded under `SWIG_WEB`;
- exception-only validation paths return an error-safe empty result because the browser
  core is built with exceptions disabled; and
- the CLI, native thread pool, CUDA backend, and process I/O layer are not linked.

Parallelism is supplied outside the module by independent browser Workers. Every Worker
instantiates this same core and owns its own immutable germline database.
