#!/usr/bin/env bash
# WP-130 / F18(c): the arms, SEQUENTIALLY (wall clock is the measurement).
# Needs LADRUNO_DIST_BIN and PY. WALL defaults to F12's 1200 s.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
W=${WALL:-1200}
A="$HERE/run_arm.sh"
echo "=== s1_T2: IntScheme 1, TanType 2 (F12 baseline) ===";           bash $A s1_T2 1 $W
echo "=== s2_T2: IntScheme 2, TanType 2, vanilla CPPM (F12 refuted) ==="; bash $A s2_T2 2 $W
echo "=== s2_ref0: -cppmOnFail refuse -cppmHalvings 0 ===";              bash $A s2_ref0 2 $W -- -cppmOnFail refuse -cppmHalvings 0
echo "=== s2_ref0_se_ls: + -cppmStart explicit -cppmLineSearch on ===";  bash $A s2_ref0_se_ls 2 $W -- -cppmOnFail refuse -cppmHalvings 0 -cppmStart explicit -cppmLineSearch on
echo "=== s2_ref3_se_ls: -cppmHalvings 3 + start + LS ===";              bash $A s2_ref3_se_ls 2 $W -- -cppmOnFail refuse -cppmHalvings 3 -cppmStart explicit -cppmLineSearch on
echo "=== s2_se_ls: vanilla fallback + start + LS (no refusal) ===";     bash $A s2_se_ls 2 $W -- -cppmStart explicit -cppmLineSearch on
echo "=== DONE ==="
