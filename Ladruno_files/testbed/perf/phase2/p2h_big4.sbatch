#!/bin/bash
#SBATCH --job-name=adr75-big4
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --time=05:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ADR-75 items 1+2, LARGE end, take 3.
#
# Take 1 (144448): `mpirun` -> instant ORTE plm_slurm failure across nodes.
# Take 2 (144449): switched to srun, still nothing -- and I could not see WHY,
#   because the script piped srun's output straight into a grep whose patterns
#   did not match `srun: fatal:`. That is the banked rule "never let a bench
#   script grep away its own log", violated in the very run meant to honour it.
# Take 3: EVERY run's full output goes to its own file; grep reads the FILE.
# Also dropped #SBATCH --exclusive: srun steps inherit it and that is a known
# way for a step to fail to launch inside an exclusive allocation.

set -uo pipefail
echo "HOST: $(hostname)"; date
echo "NODES: $SLURM_JOB_NODELIST"
export LD_LIBRARY_PATH="/opt/openmpi/lib:/mnt/nfshare/lib:${LD_LIBRARY_PATH:-}"
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
NP=32
LOGDIR=/mnt/deadmanschest/nmorabowen/apegmsh_jobs/p2h_logs_$SLURM_JOB_ID
mkdir -p "$LOGDIR"
echo "LOGDIR: $LOGDIR"

SIZES=("40 40 80" "52 52 104")

for s in "${SIZES[@]}"; do
  set -- $s
  NX=$1; NY=$2; NZ=$3
  DOF=$(( 3 * (NX+1) * (NY+1) * (NZ+1) ))
  for m in full blr sym; do
    LOG="$LOGDIR/${NX}x${NY}x${NZ}_${m}.log"
    echo "##### SIZE ${NX}x${NY}x${NZ} DOF~${DOF} MODE=$m NP=$NP -> $LOG"
    SMOKE_NX=$NX SMOKE_NY=$NY SMOKE_NZ=$NZ SMOKE_MODE=$m \
      timeout 4200 /opt/slurm/bin/srun --ntasks="$NP" --cpu-bind=cores --mpi=pmix_v3 \
        "$WRAP" p2h_smoke.tcl > "$LOG" 2>&1
    rc=$?
    # assert on ARTIFACTS, not rc
    if grep -q "P2H_RESULT" "$LOG"; then
      grep -E "P2H_RESULT|P2H_TIP|INFOG\(9\)|INFOG\(21\)|INFOG\(22\)|RINFOG\(3\)" "$LOG" | sort -u
    else
      echo "!!! NO P2H_RESULT (rc=$rc) -- first 25 lines of $LOG:"
      head -25 "$LOG"
    fi
    echo "##### END ${NX}x${NY}x${NZ} $m rc=$rc"
  done
done
echo "BIG4_SWEEP_DONE"; date
