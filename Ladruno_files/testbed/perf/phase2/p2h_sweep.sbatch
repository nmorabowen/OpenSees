#!/bin/bash
#SBATCH --job-name=adr75-blr
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ADR-75 item 1+2 — MUMPS BLR vs full-rank vs symmetric half-storage at SCALE.
# The prior desktop verdict (#626) was explicitly "small end only": BLR's payoff
# is governed by FRONT size (~N^(2/3) DOF, front memory ~N^(4/3)), so the whole
# question is whether the sign flips as N grows. This sweeps N.

set -uo pipefail
echo "HOST: $(hostname)"; date
export PATH=/opt/openmpi/bin:$PATH
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
NP=${NP:-16}

# (NX NY NZ) -> approx DOF = 3*(NX+1)*(NY+1)*(NZ+1)
SIZES=("14 14 28" "28 28 56" "40 40 80" "52 52 104")

for s in "${SIZES[@]}"; do
  set -- $s
  NX=$1; NY=$2; NZ=$3
  DOF=$(( 3 * (NX+1) * (NY+1) * (NZ+1) ))
  for m in full blr sym symblr; do
    echo "##### SIZE ${NX}x${NY}x${NZ} DOF~${DOF} MODE=$m NP=$NP"
    SMOKE_NX=$NX SMOKE_NY=$NY SMOKE_NZ=$NZ SMOKE_MODE=$m \
      timeout 2400 mpirun -n "$NP" "$WRAP" p2h_smoke.tcl 2>&1 \
      | grep -E "P2H_RESULT|P2H_TIP|INFOG\(9\)|INFOG\(21\)|INFOG\(22\)|RINFOG\(3\)|Warning|ERROR" \
      | sort -u
    echo "##### END ${NX}x${NY}x${NZ} $m rc=$?"
  done
done
echo "SWEEP_DONE"; date
