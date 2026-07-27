#!/bin/bash
#SBATCH --job-name=adr75b-gl3t
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=03:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -uo pipefail
echo "HOST: $(hostname)"; date
export PATH=/opt/openmpi/bin:$PATH
BASE=/mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=$BASE/gl3t_$SLURM_JOB_ID; mkdir -p "$LOGDIR"; cd "$LOGDIR"
echo "LOGDIR: $LOGDIR"
# Same config as the 540 675-DOF G-L3 point (np=16, LOADF=0.35 plastic,
# -matrixType 2, LadrunoParallelRCM) at two smaller N, so the element-vs-solve
# split can be read as a TREND in N rather than a single number.
for s in "10 10 20" "16 16 32" "28 28 56"; do
  set -- $s; NX=$1; NY=$2; NZ=$3
  DOF=$(( 3*(NX+1)*(NY+1)*(NZ+1) ))
  LOG="$LOGDIR/n${DOF}.log"
  echo "##### N~$DOF"
  GL3_NX=$NX GL3_NY=$NY GL3_NZ=$NZ GL3_LOADF=0.35 GL3_NSTEP=2 GL3_PROF=1 \
  GL3_OUT="$LOGDIR/p${DOF}" \
    timeout 6000 mpirun -n 16 "$WRAP" $BASE/gl3_deck.tcl > "$LOG" 2>&1
  rc=$?
  if grep -q GL3_RESULT "$LOG"; then grep -E "GL3_RESULT|GL3_TIP" "$LOG" | head -2
  else echo "!!! NO GL3_RESULT rc=$rc"; tail -12 "$LOG"; fi
  echo "##### END $DOF rc=$rc"
done
echo "GL3T_DONE"; date
