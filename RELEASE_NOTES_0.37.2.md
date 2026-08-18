# Swig 0.37.2

## Standalone CLI packaging fix

- Release executables retain Bun minification but no longer use Bun bytecode compilation.
- Bun 1.3.14 bytecode compilation failed to embed the multi-entry `swig-worker` module, so tagged `0.37.1` executables compiled but could not start their annotation workers.
- The release workflow continues to smoke-test the actual compiled Linux and macOS executables before publishing any assets.

This patch changes only standalone executable packaging. It retains the `0.37.1` RIAT-MP default and all `0.37.0` SwiftIG kernel and AIRR-output optimizations.
