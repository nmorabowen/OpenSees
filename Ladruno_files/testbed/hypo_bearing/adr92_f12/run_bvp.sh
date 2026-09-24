#!/usr/bin/env bash
# WP-105 / F12 phase (b)/(c): the CP1/ADR-95 bearing deck, purpose-sized x10z8,
# leg h1.0_e0.6944, three arms. SEQUENTIAL -- concurrent legs would corrupt the
# wall-clock column, which is the measurement.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p bvp logs
WALL=1200
COMMON="--legs h1.0_e0.6944 --xlim 10 --zbot 8 --wall $WALL --maxsubsteps 20000"

echo "=== B1: IntScheme 1, TanType 0 (the parser default) ==="
mkdir -p bvp/s1_T0
python3.12 -u f12_bvp.py 1 --out bvp/s1_T0 $COMMON --tantype 0 > logs/bvp_s1_T0.log 2>&1
tail -6 logs/bvp_s1_T0.log

echo "=== B2: IntScheme 1, TanType 2 (the deck baseline) ==="
mkdir -p bvp/s1_T2
python3.12 -u f12_bvp.py 1 --out bvp/s1_T2 $COMMON --tantype 2 > logs/bvp_s1_T2.log 2>&1
tail -6 logs/bvp_s1_T2.log

echo "=== B3: IntScheme 2, TanType 2 ==="
mkdir -p bvp/s2_T2
python3.12 -u f12_bvp.py 2 --out bvp/s2_T2 $COMMON --tantype 2 > logs/bvp_s2_T2.log 2>&1
tail -6 logs/bvp_s2_T2.log

echo "=== DONE ==="
