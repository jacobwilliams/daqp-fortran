#!/usr/bin/env bash
#
# Build the upstream DAQP C library (the git submodule in upstream/daqp) into
# build/upstream/libdaqp.a, for the comparison in compare/.
#
# The package itself never needs this: fpm does not compile the C sources.
#
# Environment variables:
#   CC       C compiler (default: gcc if available, else cc)
#   CFLAGS   optimization flags (default: -O3). Use the same level as the
#            Fortran comparison build. Upstream's CMake build adds
#            -fassociative-math -fno-signed-zeros -fno-trapping-math
#            (FASTER_MATH); it is left out by default, so that the C code
#            sums in the same order as the Fortran port.
#
# Usage: tools/build_upstream.sh   (from the repository root)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/upstream/daqp"
OUT="$ROOT/build/upstream"

if [ ! -f "$SRC/src/daqp.c" ]; then
    echo "upstream/daqp is empty: run 'git submodule update --init' first" >&2
    exit 1
fi

if [ -z "${CC:-}" ]; then
    if command -v gcc >/dev/null 2>&1; then CC=gcc; else CC=cc; fi
fi
CFLAGS="${CFLAGS:--O3}"

mkdir -p "$OUT/obj"
rm -f "$OUT"/obj/*.o "$OUT/libdaqp.a"

# double precision, individual soft weights, and PROFILING (upstream's defaults:
# the time limit and the timing of setup and solve, as the port has them)
for f in "$SRC"/src/*.c; do
    "$CC" $CFLAGS -DPROFILING -I"$SRC/include" -c "$f" -o "$OUT/obj/$(basename "${f%.c}").o"
done
ar rcs "$OUT/libdaqp.a" "$OUT"/obj/*.o

{
    echo "upstream commit: $(git -C "$SRC" rev-parse HEAD) ($(git -C "$SRC" describe --tags 2>/dev/null || echo untagged))"
    echo "compiler: $("$CC" --version | head -1)"
    echo "flags: $CFLAGS -DPROFILING"
} > "$OUT/build_info.txt"
cat "$OUT/build_info.txt"
echo "built $OUT/libdaqp.a"
