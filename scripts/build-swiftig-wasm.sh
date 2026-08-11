#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
wasi_sdk="${WASI_SDK:-}"

if [[ -z "$wasi_sdk" ]]; then
  echo "Set WASI_SDK to an extracted WASI SDK directory (WASI SDK 25 or newer)." >&2
  exit 2
fi

compiler="$wasi_sdk/bin/clang++"
if [[ ! -x "$compiler" ]]; then
  echo "No executable clang++ was found at $compiler." >&2
  exit 2
fi

build_dir="$project_dir/.build/wasm"
mkdir -p "$build_dir"
export TMPDIR="$build_dir"
export LD_LIBRARY_PATH="$wasi_sdk/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

sources=(
  "$project_dir/wasm/src/types.cpp"
  "$project_dir/wasm/src/index.cpp"
  "$project_dir/wasm/src/alignment.cpp"
  "$project_dir/wasm/src/engine.cpp"
  "$project_dir/wasm/src/airr.cpp"
  "$project_dir/wasm/src/web_adapter.cpp"
)

"$compiler" \
  --target=wasm32-wasip1 \
  -O3 \
  -DNDEBUG \
  -std=c++20 \
  -fno-exceptions \
  -fno-rtti \
  -DSWIG_WEB \
  -I"$project_dir/wasm/include" \
  "${sources[@]}" \
  -mexec-model=reactor \
  -Wl,--export-memory \
  -Wl,--strip-all \
  -o "$build_dir/swiftig.wasm"

install -m 755 "$build_dir/swiftig.wasm" "$project_dir/public/swiftig.wasm"
echo "Built public/swiftig.wasm"
