#!/bin/sh
# WP-106 / ADR-93 II.1 -- the BVP leg the Gauss-point work could not decide:
# the strip footing's FREE-SURFACE RING at the campaign's own 7.65 kPa surcharge
# (PM-01 D20's minimum-embedment pressure), pRe = 0 vs pRe = 1 kPa.
#
# COARSENED ON PURPOSE, and said out loud: the COARSE leg only (h0 = 1.0 m,
# Gorini's calibrated e_init = 0.6944), `--wall 2000` (33 min) per arm, and
# `--maxsubsteps 20000` on BOTH arms. GATE U measured every uncapped leg of this
# deck seizing inside ModifiedEuler, so uncapped both arms would simply spend
# their wall budget and measure the budget. The cap is applied identically to
# both arms, so the substep census below compares equals.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
export LADRUNO_DIST_BIN="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/wp-106-sanisand-pre/dist/bin"
export LADRUNO_A2_EXPECT_BUILD=any
cd "$HERE"
for PRE in 0 1.0; do
  OUT="wp106_bvp/pRe$PRE"
  mkdir -p "$OUT"
  echo "=== pRe = $PRE ==="
  python3.12 -u sanisand_tau0_band.py --out "$OUT" --legs h1.0_e0.6944 \
      --surcharge 7.65 --maxsubsteps 20000 --wall 2000 --pRe "$PRE"
done
echo "WP106 BVP DONE"
