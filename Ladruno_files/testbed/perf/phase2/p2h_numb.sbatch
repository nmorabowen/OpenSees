#!/bin/bash
#SBATCH --job-name=adr75-numb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=02:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ADR-75 P2h -- CONFOUND CONTROL the first sweep missed.
# The sweep used vanilla `numberer ParallelRCM`; the fork's own mp_blr_smoke.py
# prefers `LadrunoParallelRCM` (ADR-74 exists because the vanilla parallel
# numberer degrades at scale). Numbering is additive and identical across solver
# modes, so it dilutes wall RATIOS toward 1.0 -- it cannot flip a verdict, but it
# can understate one. This measures how much.
set -uo pipefail
echo "HOST: $(hostname)"; date
export PATH=/opt/openmpi/bin:$PATH
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=/mnt/deadmanschest/nmorabowen/apegmsh_jobs/p2h_numb_$SLURM_JOB_ID
mkdir -p "$LOGDIR"

for nb in ParallelRCM LadrunoParallelRCM; do
  for m in full sym; do
    LOG="$LOGDIR/${nb}_${m}.log"
    echo "##### NUMBERER=$nb MODE=$m"
    SMOKE_NX=28 SMOKE_NY=28 SMOKE_NZ=56 SMOKE_MODE=$m SMOKE_NUMBERER=$nb \
      timeout 2400 mpirun -n 16 "$WRAP" p2h_smoke2.tcl > "$LOG" 2>&1
    rc=$?
    if grep -q "P2H_RESULT" "$LOG"; then
      grep -E "P2H_RESULT|P2H_NUMBERER|INFOG\(21\)" "$LOG" | sort -u
    else
      echo "!!! NO P2H_RESULT rc=$rc"; tail -15 "$LOG"
    fi
    echo "##### END $nb $m rc=$rc"
  done
done
echo "NUMB_TEST_DONE"; date
