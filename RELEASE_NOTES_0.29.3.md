# Swig 0.29.3 release notes

## Windows release cross-compilation

- The Windows x86-64 executable is cross-compiled on Linux because Bun 1.3.14 fails during compilation on the GitHub Windows runner.
- Native executable smoke tests remain enabled for Linux x86-64/ARM64 and macOS Intel/Apple silicon. The Windows artifact is compiled and uploaded but not executed in CI.

No analysis method, default, or output schema changed from 0.29.2.
