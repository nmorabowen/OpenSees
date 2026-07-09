#!/bin/bash
# ADR-43 L2 profiler — SLURM batch submission for the esmeralda cluster.
# Runs the L3-only distributed-FEAST strong-scaling sweep on ONE exclusive
# compute node (node* = 32 CPUs / 60 GB — 2x the head node's RAM). The profiler
# spawns `mpirun -n {np}` internally per data point, using this allocation's
# slots. Submit from the repo's feast_l2_profile dir:
#     sbatch l2_sbatch.sh            # defaults: ne=32, np=2,4,8,16
#     NE=40 NPS=2,4,8,16,32 sbatch l2_sbatch.sh
# Multi-node (for the np=32-64 "real budget" regime), override the header:
#     sbatch --nodes=2 --ntasks-per-node=32 l2_sbatch.sh
#
#SBATCH --job-name=feast-l2
#SBATCH --partition=computes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --output=%x_%j.log
set -e

SPACK_MKL=/mnt/nfshare/software/spack/opt/spack/linux-zen3/intel-oneapi-mkl-2024.2.2-4olf4bnrh4yzaqbhq742doqvnim6yynl/mkl/2024.2/lib
export PATH=/opt/openmpi/bin:$PATH
# The profiler spawns `mpirun -n {np}` internally. When this script itself runs
# ON a compute node (launched via `srun ... bash l2_sbatch.sh`), force the inner
# mpirun to launch ranks LOCALLY on that node (rsh/local launcher) rather than
# opening a nested SLURM step, which conflicts with the enclosing srun step.
# Result: all np ranks run on the srun-assigned compute node. Override with
# OMPI_MCA_plm=slurm for a true multi-node sbatch run.
export OMPI_MCA_plm=${OMPI_MCA_plm:-rsh}
export OMPI_MCA_rmaps_base_oversubscribe=1
export TMPDIR=$HOME/ladruno_build_test/tmp
export L2_MODDIR=$HOME/ladruno_build_test/OpenSees/build/Release
export L2_MPIEXEC=/opt/openmpi/bin/mpirun
export L2_LDPATH=$SPACK_MKL:/opt/openmpi/lib:/mnt/nfshare/lib
export L2_TMPDIR=$HOME/ladruno_build_test/tmp
# inside a SLURM allocation OpenMPI maps -n ranks onto the allotted slots;
# --oversubscribe lets np exceed physical cores if needed; keep MPI tmp off /tmp.
# -x forwards the runtime env to REMOTE-node ranks (OpenMPI does NOT forward
# LD_LIBRARY_PATH/MKL vars by default — without this, compute-node ranks can't
# load MKL/openseesmp). Works for salloc (this node) and sbatch multi-node alike.
export L2_MPI_FLAGS="--oversubscribe --mca orte_tmpdir_base $HOME/ladruno_build_test/tmp -x LD_LIBRARY_PATH -x MKL_INTERFACE_LAYER -x MKL_THREADING_LAYER -x MKL_NUM_THREADS -x OMP_NUM_THREADS -x LADRUNO_FEAST_MPI_DEBUG -x LADRUNO_OPENSEES_QUIET -x LADRUNO_FEAST_PHI"
export L2_M0=${L2_M0:-16}

NE=${NE:-32}
NPS=${NPS:-2,4,8,16}
WANT=${WANT:-8}

cd "$HOME/ladruno_build_test/OpenSees/Ladruno_implementation/feast_l2_profile"
echo "host=$(hostname) ne=$NE nps=$NPS want=$WANT m0=$L2_M0"
date
python3.10 l2_profile.py "$NE" "$NPS" "$WANT"
date
