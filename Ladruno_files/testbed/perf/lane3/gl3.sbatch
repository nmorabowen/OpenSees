#!/bin/bash
#SBATCH --job-name=adr75b-gl3
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=04:00:00
#SBATCH --partition=computes
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -uo pipefail
echo "HOST: $(hostname)"; date
export PATH=/opt/openmpi/bin:$PATH
cd /mnt/deadmanschest/nmorabowen/apegmsh_jobs
WRAP=/mnt/deadmanschest/nmorabowen/ladruno_build_test/openseesmp.sh
LOGDIR=/mnt/deadmanschest/nmorabowen/apegmsh_jobs/gl3_$SLURM_JOB_ID
mkdir -p "$LOGDIR"; cd "$LOGDIR"
echo "LOGDIR: $LOGDIR"
# 44x44x88 = 540 675 DOF (>= the 500k G-L3 threshold), moderate plasticity
# (LOADF 0.35 -> tip 22% above the elastic response at smoke scale, so elem.update
# does real constitutive work), coarse -perStep profiling, np=16 on ONE node.
LOG="$LOGDIR/gl3.log"
GL3_NX=44 GL3_NY=44 GL3_NZ=88 GL3_LOADF=0.35 GL3_NSTEP=2 GL3_PROF=1 \
GL3_OUT="$LOGDIR/gl3prof" \
  timeout 12000 mpirun -n 16 "$WRAP" /mnt/deadmanschest/nmorabowen/apegmsh_jobs/gl3_deck.tcl > "$LOG" 2>&1
rc=$?
if grep -q "GL3_RESULT" "$LOG"; then
  echo "GL3_OK"; grep -E "GL3_|INFOG\(21\)" "$LOG" | sort -u
else
  echo "!!! NO GL3_RESULT rc=$rc"; tail -25 "$LOG"
fi
ls -la "$LOGDIR"/*.h5 2>/dev/null | head
echo "GL3_DONE"; date
