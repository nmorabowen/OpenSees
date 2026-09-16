#!/bin/sh
# WP-106 / ADR-93 II.1 -- capacity neutrality at p0 = 20 kPa.
# UNCAPPED (no -maxSubsteps): the cap FAILS the update, which under
# LadrunoBrick's sentinel forwarding stalls the push at ez ~ 0.7 % -- well short
# of the peak this test needs. Uncapped is vanilla's own behaviour and both arms
# run it.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
POST="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/wp-106-sanisand-pre/dist/bin"
cd "$HERE"
mkdir -p data
run() { echo "=== $*"; python3.12 -u probe_pre.py "$@"; }
N="--nstep 90 --ez-max 0.036"
run --bin "$POST" --p0 20 --pRe 0    $N --out data/cap_p20_pre0.csv
run --bin "$POST" --p0 20 --pRe 1.0  $N --out data/cap_p20_pre1.csv
run --bin "$POST" --p0 20 --pRe 0.1  $N --out data/cap_p20_pre0p1.csv
echo "CAPACITY DONE"
