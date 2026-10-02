#!/bin/bash
#SBATCH --job-name=wp144-p3smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --partition=computes
#SBATCH --output=/mnt/deadmanschest/nmorabowen/wp144/p3_smoke-%j.log
# WP-144 P3 harness gate (smoke_recovery.py) on a compute node. Python only (numpy/scipy, ~/wp144/venv).
# Submit from the head node:  sbatch ~/wp144/repo/Ladruno_files/testbed/norsand_calib/harness/run_smoke.sbatch
set -uo pipefail
export PATH=$HOME/wp144/venv/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $HOME/wp144/repo/Ladruno_files/testbed/norsand_calib || exit 1
echo "host=$(hostname) start=$(date -Is) job=${SLURM_JOB_ID:-none}"
python -u -m harness.smoke_recovery --starts ${STARTS:-8} --workers ${WORKERS:-8} ${SMOKE_ARGS:-}
rc=$?
echo "end=$(date -Is) rc=$rc"
exit $rc
