#!/bin/sh
# WP-106 / ADR-93 II.1 -- the binary probe matrix.
#   PRE  = the pre-change binary at ladruno 634824e1f
#   POST = this worktree's build
# Every leg is a drained triaxial on one LadrunoBrick, IntScheme 1, TanType 2,
# at the SAME axial step du = 4e-4, so the substep census compares equals.
#
# `-maxSubsteps 20000`: uncapped, the low-confinement legs are the ADR-93
# seizure itself and do not terminate in a session. The cap is what the campaign
# runs (ADR-86b) and it is applied to EVERY arm.
#
# Step counts are deliberately short (the p0 = 100 leg's own global solve stalls
# at ez = 2.75 % anyway): every comparison below is arm-against-arm on the SAME
# number of steps of the SAME path, so a peak that the push has not reached yet
# is a shared property of both arms, not a bias in one.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
PRE="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/release-build-634824e1f/dist/bin"
POST="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/wp-106-sanisand-pre/dist/bin"
cd "$HERE"
mkdir -p data

run() { echo "=== $*"; python3.12 -u probe_pre.py "$@"; }

N20="--nstep 50 --ez-max 0.020 --max-substeps 20000"
NLO="--nstep 40 --ez-max 0.016 --max-substeps 20000"

# (2) the certificate: pRe = 0 must be byte-identical to the pre-change binary
run --bin "$PRE"  --p0 20  --pRe 0   $N20 --out data/cert_p20_PRE.csv
run --bin "$POST" --p0 20  --pRe 0   $N20 --out data/cert_p20_POST.csv

# (3) capacity neutrality at p0 = 20 kPa, pRe = 1 kPa
run --bin "$POST" --p0 20  --pRe 1.0 $N20 --out data/tx_p20_pre1.csv

# the decisive cost test: Gauss points that actually REACH the floor
run --bin "$PRE"  --p0 0.5 --pRe 0    $NLO --out data/cert_p0.5_PRE.csv
run --bin "$POST" --p0 0.5 --pRe 0    $NLO --out data/tx_p0.5_pre0.csv
run --bin "$POST" --p0 0.5 --pRe 0.1  $NLO --out data/tx_p0.5_pre0p1.csv
run --bin "$POST" --p0 0.5 --pRe 1.0  $NLO --out data/tx_p0.5_pre1.csv
run --bin "$POST" --p0 0.1 --pRe 0    $NLO --out data/tx_p0.1_pre0.csv
run --bin "$POST" --p0 0.1 --pRe 1.0  $NLO --out data/tx_p0.1_pre1.csv
echo "ALL PROBES DONE"
