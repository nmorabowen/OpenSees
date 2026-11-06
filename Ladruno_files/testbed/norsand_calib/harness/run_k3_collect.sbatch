#!/bin/bash
#SBATCH --job-name=wp144-k3-collect
#SBATCH --partition=computes
#SBATCH --exclude=node5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.log
# WP-144 K3: collect the array's checkpoints (run_k3_array.sbatch), then identifiability and the chi / h profiles
# from the best start, and write the final record to $RUN/out. Submit with --dependency=afterany:<array job>, so a
# partial array (a task at its time limit) is still collected; the record lists the missing starts.
set -euo pipefail
: "${RUN:?RUN (the K3 run directory) must be exported}"
export PATH=$HOME/wp144/venv/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RUN/snap/testbed/norsand_calib"
python -u -m harness.fit_tatsuoka --lit --energy HAR --p-min 0.505 --pi0-rule unified \
    --free chi,h,N,N_bar,rho,rho_bar --starts 16 --seed 0 --workers 16 --profiles \
    --checkpoint "$RUN/ckpt" --collect --out "$RUN/out" --tag K3_HAR_lit
