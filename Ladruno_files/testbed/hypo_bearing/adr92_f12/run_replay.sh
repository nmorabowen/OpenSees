#!/usr/bin/env bash
# WP-105 / F12 phase (a), the CLEAN certificate: replay the scheme-1 baseline's
# own recorded strain path on a zero-free-DOF cube under both schemes. No global
# Newton, so no solver stall can confound the constitutive comparison, and wall
# time is pure constitutive cost. (ADR-92 G0's discipline, one level up.)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p logs data
export LADRUNO_OPENSEES_QUIET=1

run () {
  tag="$1"; shift
  echo "=== $tag ==="
  python3.12 -u f12_matpoint.py --tag "$tag" "$@" > "logs/$tag.log" 2>&1
  rc=$?
  grep -h "@@SUMMARY" "logs/$tag.log" || echo "NO SUMMARY (rc=$rc) for $tag"
}

for P in 100 20; do
  SRC="data/tx_p${P}_s1_nom.csv"
  for N in 1820 182 40; do
    for S in 1 2; do
      run "rp_p${P}_s${S}_n${N}" --kind replay --p0 $P --scheme $S --tan-type 2 \
          --nstep $N --path-csv "$SRC"
    done
  done
done

echo "=== DONE ==="
