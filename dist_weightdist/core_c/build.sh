#!/usr/bin/env bash
# Build the dwd_core hot kernel.  Produces a single self-contained binary that is
# scp-able to a compute node (no Python / no toolchain dependency at run time).
#
#   ./build.sh             # native build (-march=native, AVX2 if available)
#   ./build.sh portable    # -mavx2 (portable to any AVX2 x86-64 node)
#   ./build.sh generic     # no AVX2 (scalar fallback) -- still correct, slower
#
# The frozen design names Rust+rayon primary; this host has no Rust toolchain, so
# the documented C+OpenMP fallback IS the core (identical algorithm & file formats).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/dwd_core.c"
OUT="$HERE/dwd_core"
MODE="${1:-native}"
CC="${CC:-gcc}"
case "$MODE" in
  native)   ARCH="-march=native" ;;
  portable) ARCH="-mavx2 -mtune=generic" ;;
  generic)  ARCH="-mtune=generic" ;;       # no AVX2: scalar fallback path
  *) echo "unknown mode $MODE (native|portable|generic)"; exit 2 ;;
esac
set -x
$CC -O3 -funroll-loops $ARCH -fopenmp -o "$OUT" "$SRC"
set +x
echo "built $OUT  (mode=$MODE)"
"$OUT" 2>/dev/null || true
