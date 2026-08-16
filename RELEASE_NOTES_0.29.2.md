# Swig 0.29.2 release notes

## Cross-platform release packaging

- The release workflow now creates its output directory with Node and uses each GitHub runner's native default shell.
- This avoids routing the Windows Bun compiler through Git Bash while retaining the native executable smoke test on every platform.

No analysis method, default, or output schema changed from 0.29.1.
