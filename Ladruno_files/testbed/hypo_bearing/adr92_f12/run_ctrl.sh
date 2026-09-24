#!/usr/bin/env bash
# F12 phase (a) controls:
#  (1) free-standing triaxial re-run with the patched wall accounting, at the
#      certificate's own NormDispIncr 1e-9 AND at a looser 1e-7, to separate a
#      GLOBAL-SOLVER stall from a constitutive one.
#  (2) -maxSubsteps on IntScheme 2: does the seam the class warns is inert on
#      scheme 2 actually fire through the CPPM's explicit fallback?
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p logs data
export LADRUNO_OPENSEES_QUIET=1

run () { tag="$1"; shift; echo "=== $tag ==="
  python3.12 -u f12_matpoint.py --tag "$tag" "$@" > "logs/$tag.log" 2>&1
  grep -h "@@SUMMARY" "logs/$tag.log" || echo "NO SUMMARY for $tag"; }

for P in 100 20; do
  for S in 1 2; do
    for N in 40 200; do
      run "ct_p${P}_s${S}_n${N}_t9" --kind tx --p0 $P --scheme $S --tan-type 2 \
          --nstep $N --ez-max 0.02 --push-tol 1e-9
      run "ct_p${P}_s${S}_n${N}_t7" --kind tx --p0 $P --scheme $S --tan-type 2 \
          --nstep $N --ez-max 0.02 --push-tol 1e-7
    done
  done
done

# -maxSubsteps seam on scheme 2, on the floor path (where subME reached 1282)
run "ms_fl_s2_n40_cap100"  --kind floor --p0 5 --scheme 2 --nstep 40 --max-substeps 100
run "ms_fl_s1_n40_cap100"  --kind floor --p0 5 --scheme 1 --nstep 40 --max-substeps 100

echo "=== DONE ==="
