#!/bin/bash
#SBATCH --job-name=adr75-anom1n
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=04:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# Close the 884.8k anomaly: symmetric's per-proc saving was -8.8% there vs
# -52.8% at 143.8k/np=32. The np-sweep already proved np is NOT the cause, so the
# two remaining candidates were NODE COUNT (1 -> 2) and N (6.2x) -- confounded.
# This holds np=32 fixed and runs on TWO nodes at two N, so node count can be
# compared directly against the 1-node rows already measured.
set -uo pipefail
echo "HOST: $(hostname)"; date; echo "NODES: $SLURM_JOB_NODELIST"
export PATH=/opt/openmpi/bin:$PATH; export LD_LIBRARY_PATH="/opt/openmpi/lib:/mnt/nfshare/lib:${LD_LIBRARY_PATH:-}"
BASE=/mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=$BASE/anom1n_$SLURM_JOB_ID; mkdir -p "$LOGDIR"; cd "$BASE"
echo "LOGDIR: $LOGDIR"
for s in "40 40 80"; do
  set -- $s; NX=$1; NY=$2; NZ=$3
  DOF=$(( 3*(NX+1)*(NY+1)*(NZ+1) ))
  for m in full sym; do
    LOG="$LOGDIR/n${DOF}_${m}.log"
    echo "##### N~$DOF MODE=$m NP=32 NODES=1"
    SMOKE_NX=$NX SMOKE_NY=$NY SMOKE_NZ=$NZ SMOKE_MODE=$m \
    SMOKE_NUMBERER=LadrunoParallelRCM \
      timeout 5400 mpirun -n 32 \
        "$WRAP" p2h_smoke2.tcl > "$LOG" 2>&1
    rc=$?
    if grep -q P2H_RESULT "$LOG"; then
      echo -n "  maxINFOG21="; grep "INFOG(21)" "$LOG" | awk '{print $5+0}' | sort -n | tail -1
      echo -n "  maxINFOG22="; grep "INFOG(22)" "$LOG" | awk '{print $6+0}' | sort -n | tail -1
      grep P2H_RESULT "$LOG"
    else echo "!!! NO P2H_RESULT rc=$rc"; tail -15 "$LOG"; fi
    echo "##### END $DOF $m rc=$rc"
  done
done
echo "ANOM1N_DONE"; date
