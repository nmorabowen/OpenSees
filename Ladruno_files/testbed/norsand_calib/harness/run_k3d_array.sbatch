#!/bin/bash
#SBATCH --job-name=wp144-k3d
#SBATCH --partition=computes,new
#SBATCH --exclude=node5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                  # Slurm 24.11 allocates CR_CORE_MEMORY with no default: without it one task holds the whole node (58 GB); a start uses ~125 MB
#SBATCH --array=0-16              # 16 Sobol starts + start 16 = --warm (K3 start 11, pi0_ratio 1, g_fac 1)
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%A_%a.log
# WP-144 K3c (NorSand's own CSL e0, lambda_c, xi free; owner 2026-11-09; else as K3) Toyoura fit, one Sobol start per array task (Python only: numpy/scipy, ~/wp144/venv).
#
# Why an array (2026-11-06): job 151170 (2026-10-04) ran all 16 starts in one 12 h job and was killed at the
# limit with no output. One full objective evaluation (--lit: 19 curves + 22 point tests, HAR) takes 42 s; a
# start takes up to max_nfev (k+1) = 1400 evaluations (~16 h), typically 44-294 (smoke). Each task now writes
# its start to $CKPT/start_<k>.json when it finishes; a resubmitted task skips a start already on disk.
#
# Submit from the K3 run directory (a snapshot of norsand_calib + norsand_oracle under snap/testbed):
#   jid=$(sbatch --parsable --export=ALL,RUN=$PWD harness_path/run_k3_array.sbatch)
#   sbatch --dependency=afterany:$jid --export=ALL,RUN=$PWD harness_path/run_k3_collect.sbatch
# WP-144 K3d (owner 2026-11-12): the six plastic constants plus the initial image-pressure ratio (pi0_ratio) and one
# factor on the elastic stiffness (g_fac), with the sr residuals up to 0.5 % axial strain weighted by 3.
set -euo pipefail
: "${RUN:?RUN (the K3 run directory) must be exported}"
export PATH=$HOME/wp144/venv/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RUN/snap/testbed/norsand_calib"
python -u -m harness.fit_tatsuoka --lit --energy HAR --p-min 0.505 --pi0-rule unified \
    --free chi,h,N,N_bar,rho,rho_bar,pi0_ratio,g_fac --small-strain 0.5 --small-weight 3 --starts 16 --seed 0 --workers 1 --warm "$RUN/warm.json" \
    --checkpoint "$RUN/ckpt" --only "${SLURM_ARRAY_TASK_ID}"
