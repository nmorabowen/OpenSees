#!/bin/bash
# WP-238 biaxial demo on esmeralda: one single-threaded case as a 1-CPU srun step.
#   bash start_case.sh <case> <biaxial.py args...>
# Output: ~/ladruno_wp238/deck/out/<case>/{run.npz,log.txt}, stdout in out/<case>.txt
set -u
CASE=$1; shift
export PATH=/opt/slurm-24.11/bin:$PATH
ROOT=$HOME/ladruno_wp238
D=$ROOT/deck
M=/mnt/nfshare/software/spack/opt/spack/linux-zen3/intel-oneapi-mkl-2024.2.2-4olf4bnrh4yzaqbhq742doqvnim6yynl/mkl/2024.2
export LD_LIBRARY_PATH=$ROOT/lib:$M/lib:${LD_LIBRARY_PATH:-}   # lib/: libpython3.10, libgfortran, libquadmath copied from the head node (the nodes lack them)
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 LADRUNO_OPENSEES_QUIET=1
export BIAX_BIN=$ROOT/bin BIAX_SCRIPTS=$ROOT/OpenSees/Ladruno_scripts
export TMPDIR=$ROOT/tmp
PY=$HOME/ladruno_build_test/conan_venv/bin/python
mkdir -p $D/out
cd $D
setsid nohup srun --job-name=wp238_$CASE --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G \
  $PY -u $D/biaxial.py --out out/$CASE "$@" > $D/out/$CASE.txt 2>&1 < /dev/null &
sleep 3
echo "$(date -Is) $CASE: $(squeue -h -u $USER -n wp238_$CASE -o '%i %N %T' | head -1) args: $*" | tee -a $D/JOBS.txt
