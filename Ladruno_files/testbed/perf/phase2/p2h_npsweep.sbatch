#!/bin/bash
#SBATCH --job-name=adr75-npsweep
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=03:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ADR-75 P2h follow-up — CLOSE THE np CONFOUND.
#
# The phase2 sweep reported that symmetric's PER-PROCESS peak saving shrinks with
# rank count (-43.1% at np=16 -> -8.8% at np=32) while TOTAL factor memory stays
# ~-43%. But N and np moved together across those rows, so it was recorded as an
# observed pattern, not a law. This isolates np.
#
# DESIGN: N is FIXED at 143 811 DOF and every np runs on ONE node (32 cores), so
# there is no interconnect variable either -- the only thing changing is np.
# 143.8k is chosen because full-rank at np=16 was 993 MB/proc = ~16 GB total, so
# every np in {4,8,16,32} fits comfortably in 60 GB. (408.5k would OOM at np=32
# on one node: ~2 GB/proc x 32 = 64 GB.)
#
# Uses LadrunoParallelRCM throughout -- the numberer the first sweep should have
# used.

set -uo pipefail
echo "HOST: $(hostname)"; date
export PATH=/opt/openmpi/bin:$PATH
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=/mnt/deadmanschest/nmorabowen/apegmsh_jobs/p2h_np_$SLURM_JOB_ID
mkdir -p "$LOGDIR"
echo "LOGDIR: $LOGDIR"

for NP in 4 8 16 32; do
  for m in full sym; do
    LOG="$LOGDIR/np${NP}_${m}.log"
    echo "##### NP=$NP MODE=$m  (N fixed = 143811 DOF, 1 node)"
    SMOKE_NX=28 SMOKE_NY=28 SMOKE_NZ=56 SMOKE_MODE=$m \
    SMOKE_NUMBERER=LadrunoParallelRCM \
      timeout 3000 mpirun -n "$NP" "$WRAP" p2h_smoke2.tcl > "$LOG" 2>&1
    rc=$?
    if grep -q "P2H_RESULT" "$LOG"; then
      echo -n "  maxINFOG21="; grep "INFOG(21)" "$LOG" | awk '{print $5+0}' | sort -n | tail -1
      echo -n "  maxINFOG22="; grep "INFOG(22)" "$LOG" | awk '{print $6+0}' | sort -n | tail -1
      grep -E "P2H_RESULT" "$LOG"
    else
      echo "!!! NO P2H_RESULT (rc=$rc)"; tail -15 "$LOG"
    fi
    echo "##### END NP=$NP $m rc=$rc"
  done
done
echo "NPSWEEP_DONE"; date
