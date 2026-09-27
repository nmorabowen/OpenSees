#!/usr/bin/env bash
# WP-130 / F18(c): one arm of F12's bearing deck (x10z8, leg h1.0_e0.6944),
# the same COMMON options as ../adr92_f12/run_bvp.sh. Run arms SEQUENTIALLY --
# the wall-clock column is the measurement.
#   run_arm.sh <name> <IntScheme> <wall_s> [--predictor] [-- <LadrunoSANISAND flags>]
# Needs LADRUNO_DIST_BIN (the WP-130 dist/bin) and PY (CPython 3.12).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
NAME=$1; SCHEME=$2; WALL=$3; shift 3
PRED=""
if [ "${1:-}" = "--predictor" ]; then PRED="--predictor"; shift; fi
mkdir -p bvp/$NAME logs
PY=${PY:-python3.12}
"$PY" -u f130_bvp.py $SCHEME --out bvp/$NAME --legs h1.0_e0.6944 --xlim 10 --zbot 8 \
    --wall $WALL --maxsubsteps 20000 --tantype 2 $PRED "$@" > logs/$NAME.log 2>&1
tail -4 logs/$NAME.log
