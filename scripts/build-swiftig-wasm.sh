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

wasm_opt="$project_dir/node_modules/.bin/wasm-opt"
if [[ ! -x "$wasm_opt" ]]; then
  echo "Run npm install before building so the pinned Binaryen wasm-opt is available." >&2
  exit 2
fi

build_dir="$project_dir/.build/wasm"
mkdir -p "$build_dir"
export TMPDIR="$build_dir"
export LD_LIBRARY_PATH="$wasi_sdk/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

extra_defines=("-DSWIG_AER_ROBUST_EXHAUSTIVE=0")
if [[ "${SWIG_AER_ROBUST_EXHAUSTIVE:-0}" == "1" ]]; then
  extra_defines=("-DSWIG_AER_ROBUST_EXHAUSTIVE=1")
fi
output_path="${SWIG_WASM_OUTPUT:-$project_dir/public/swiftig.wasm}"

sources=(
  "$project_dir/wasm/src/types.cpp"
  "$project_dir/wasm/src/index.cpp"
  "$project_dir/wasm/src/allele_tree.cpp"
  "$project_dir/wasm/src/alignment.cpp"
  "$project_dir/wasm/src/statistics.cpp"
  "$project_dir/wasm/src/engine.cpp"
  "$project_dir/wasm/src/double_d.cpp"
  "$project_dir/wasm/src/airr.cpp"
  "$project_dir/wasm/src/web_adapter.cpp"
)

"$compiler" \
  --target=wasm32-wasip1 \
  -O3 \
  -flto \
  -msimd128 \
  -mbulk-memory \
  -DNDEBUG \
  -std=c++20 \
  -fno-exceptions \
  -fno-rtti \
  -DSWIG_WEB \
  -DSWIG_V_TREE_ROOT_ALIGNMENTS=3 \
  -DSWIG_V_TREE_TRACEBACKS=2 \
  -DSWIG_V_TREE_TRACEBACK_TOLERANCE=4 \
  -DSWIG_V_TREE_TRACE_STATE_LIMIT=1024 \
  "${extra_defines[@]}" \
  -I"$project_dir/wasm/include" \
  "${sources[@]}" \
  -mexec-model=reactor \
  -Wl,--lto-O3 \
  -Wl,--gc-sections \
  -Wl,--export-memory \
  -Wl,--strip-all \
  -o "$build_dir/swiftig.raw.wasm"

"$wasm_opt" \
  "$build_dir/swiftig.raw.wasm" \
  -O4 \
  --converge \
  --enable-bulk-memory \
  --enable-simd \
  --strip-debug \
  --strip-producers \
  -o "$build_dir/swiftig.wasm"

install -m 755 "$build_dir/swiftig.wasm" "$output_path"
echo "Built optimized SIMD/LTO ${output_path#$project_dir/}"
