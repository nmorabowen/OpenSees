#!/bin/bash
#SBATCH --job-name=adr75-wall
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# ADR-75 item 2 — the CAPABILITY question, the analogue of P1c's
# "UmfPack OOM'd at 86 490 DOF; PARDISO solved it".
#
# At 884 835 DOF on ONE 60 GB node (np=16), full-rank FAILED (rc=1) and -BLR did
# not rescue it -- it drove the node to 323 MB free without completing. Here the
# same model gets TWO nodes (120 GB) and we ask whether `-matrixType 2`'s ~-43%
# clears a wall that full-rank does not.
#
# Order matters: sym FIRST, so that if the allocation dies on `full` we still
# have the answer that actually matters.

set -uo pipefail
echo "HOST: $(hostname)"; date
echo "NODES: $SLURM_JOB_NODELIST"
export LD_LIBRARY_PATH="/opt/openmpi/lib:/mnt/nfshare/lib:${LD_LIBRARY_PATH:-}"
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=/mnt/deadmanschest/nmorabowen/apegmsh_jobs/p2h_wall_$SLURM_JOB_ID
mkdir -p "$LOGDIR"
echo "LOGDIR: $LOGDIR"

NX=52; NY=52; NZ=104
DOF=$(( 3 * (NX+1) * (NY+1) * (NZ+1) ))

for m in sym full; do
  LOG="$LOGDIR/${m}.log"
  echo "##### 884835-class DOF~${DOF} MODE=$m NP=32 -> $LOG"
  SMOKE_NX=$NX SMOKE_NY=$NY SMOKE_NZ=$NZ SMOKE_MODE=$m \
    timeout 6000 /opt/slurm/bin/srun --ntasks=32 --cpu-bind=cores --mpi=pmix_v3 \
      "$WRAP" p2h_smoke.tcl > "$LOG" 2>&1
  rc=$?
  if grep -q "P2H_RESULT" "$LOG"; then
    echo "RESULT_OK $m"
    grep -E "P2H_RESULT|P2H_TIP|INFOG\(9\)|INFOG\(21\)|INFOG\(22\)" "$LOG" | sort -u
  else
    echo "!!! $m produced NO P2H_RESULT (rc=$rc). Tail of log:"
    tail -20 "$LOG"
  fi
  echo "##### END $m rc=$rc"
done
echo "WALL_TEST_DONE"; date
