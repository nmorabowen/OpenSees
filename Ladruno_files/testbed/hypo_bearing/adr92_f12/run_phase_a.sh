#!/usr/bin/env bash
# WP-105 / F12 phase (a) battery. One process per case (static latches in the
# material class must not leak between arms). Logs per case.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p logs data
export LADRUNO_OPENSEES_QUIET=1

run () {  # run <tag> <args...>
  tag="$1"; shift
  echo "=== $tag ==="
  python3.12 -u f12_matpoint.py --tag "$tag" "$@" > "logs/$tag.log" 2>&1
  rc=$?
  grep -h "@@SUMMARY" "logs/$tag.log" || echo "NO SUMMARY (rc=$rc) for $tag"
}

# --- triaxial certificate battery: p0 = 20 and 100 kPa, e = 0.6944 ---
for P in 100 20; do
  run "tx_p${P}_s1_ref"  --kind tx --p0 $P --scheme 1 --tan-type 2 --nstep 2000 --ez-max 0.02
  run "tx_p${P}_s2_ref"  --kind tx --p0 $P --scheme 2 --tan-type 2 --nstep 2000 --ez-max 0.02
  run "tx_p${P}_s1_nom"  --kind tx --p0 $P --scheme 1 --tan-type 2 --nstep 200  --ez-max 0.02
  run "tx_p${P}_s2_nom"  --kind tx --p0 $P --scheme 2 --tan-type 2 --nstep 200  --ez-max 0.02
  run "tx_p${P}_s1_c5"   --kind tx --p0 $P --scheme 1 --tan-type 2 --nstep 40   --ez-max 0.02
  run "tx_p${P}_s2_c5"   --kind tx --p0 $P --scheme 2 --tan-type 2 --nstep 40   --ez-max 0.02
done

# --- the p -> p_min floor path (ADR-93 ring regime), seeded at p0 = 5 kPa ---
for N in 160 40; do
  run "fl_p5_s1_n${N}" --kind floor --p0 5 --scheme 1 --tan-type 2 --nstep $N
  run "fl_p5_s2_n${N}" --kind floor --p0 5 --scheme 2 --tan-type 2 --nstep $N
done

echo "=== DONE ==="
