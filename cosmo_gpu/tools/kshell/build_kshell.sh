#!/usr/bin/env bash
# Build KSHELL (GaffaSnobb/kshell fork) outside this repository, as a numerical oracle.
# KSHELL's license is unclear, so it is never vendored here; only its outputs are committed.
#
# usage: tools/kshell/build_kshell.sh          (KSHELL_DIR defaults to $HOME/kshell)
set -euo pipefail

KSHELL_DIR=${KSHELL_DIR:-$HOME/kshell}
KSHELL_REPO=https://github.com/GaffaSnobb/kshell
KSHELL_COMMIT=41ea03a6707fb588e176e906edea85be2a46ac79
FC=${FC:-gfortran}

if [[ ! -d $KSHELL_DIR/.git ]]; then
    git clone "$KSHELL_REPO" "$KSHELL_DIR"
fi
git -C "$KSHELL_DIR" fetch --quiet origin
git -C "$KSHELL_DIR" checkout --quiet "$KSHELL_COMMIT"

# Link against the BLAS/LAPACK next to the compiler (the conda toolchain CoSMo uses).
LIBDIR=$(dirname "$(dirname "$(command -v "$FC")")")/lib
cd "$KSHELL_DIR/src"
make clean >/dev/null 2>&1 || true
make -j8 FC="$FC" \
    FFLAGS="-O3 -fopenmp -fallow-argument-mismatch" \
    LIBS="-L$LIBDIR -Wl,-rpath,$LIBDIR -llapack -lblas -lm" \
    kshell.exe transit.exe count_dim.exe
cp kshell.exe transit.exe count_dim.exe "$KSHELL_DIR/bin/"

{
    echo "commit $KSHELL_COMMIT"
    echo "compiler $("$FC" --version | head -1)"
    echo "libdir $LIBDIR"
    echo "built $(date -Iseconds) on $(hostname)"
} >"$KSHELL_DIR/BUILD_INFO"
cat "$KSHELL_DIR/BUILD_INFO"
