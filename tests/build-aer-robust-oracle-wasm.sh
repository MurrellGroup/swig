#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output="${1:-$project_dir/.build/aer-robust-oracle.wasm}"
mkdir -p "$(dirname "$output")"

SWIG_AER_ROBUST_EXHAUSTIVE=1 \
SWIG_WASM_OUTPUT="$output" \
  bash "$project_dir/scripts/build-swiftig-wasm.sh"
