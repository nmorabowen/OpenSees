#!/bin/bash
# WP-253 (#1097): the esmeralda rows of the MPM benchmark -- mode F at N = 10^6 and the 24-thread
# rows (pre-registration Ladruno_implementation/253_mpm_performance.md sections 1, 3, 4).
#
# Submit (the --exclude is REQUIRED on every sbatch/srun, owner's rule 2026-11-20; the #SBATCH line
# below repeats it so a bare `sbatch` cannot forget it):
#
#   sbatch --exclude=node1,node2,node12,node13 --mem-per-cpu=2G --exclusive \
#          Ladruno_scripts/mpm_bench/esmeralda/wp253_bench.sbatch <CLONE> <PYD> <EXE> <OUT> [SET]
#
#   CLONE  a FRESH clone of the WP-253 branch (never ~/ladruno_build_test/OpenSees), built there
#   PYD    the opensees python module built from that clone (opensees.so)
#   EXE    the OpenSees executable built from that clone (its ladrunoBuild must equal the PYD's)
#   OUT    results directory (created)
#   SET    big (default: N = 10^6, threads 1 4 8 16 24) | t24 (N = 10^4 10^5, threads 16 24) | all
#
# 16 threads is an extra, unregistered row: the compute nodes have 16 physical cores (2 HT each),
# so the registered 24 runs on hyperthreads. It sets no target.
# Environment: WP253_VENV (a venv with numpy + h5py; its bin goes first on PATH) and
# TCL_LIBRARY (default: the first conan tcl8.6 under ~/.conan2; without it the exe cannot
# report ladrunoBuild and the harness refuses the run).
#
# Nothing here ssh's or submits by itself; the orchestrator submits it.
#SBATCH --job-name=wp253-mpm-bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2G   # 24 x 2G = 48G < the 58 GB nodes
#SBATCH --exclude=node1,node2,node12,node13
#SBATCH --time=12:00:00
#SBATCH --output=wp253-%j.out

set -euo pipefail

CLONE=${1:?usage: wp253_bench.sbatch CLONE PYD EXE OUT [big|t24|all]}
PYD=${2:?PYD}
EXE=${3:?EXE}
OUT=${4:?OUT}
SET=${5:-big}

case "$(readlink -f "$CLONE")" in
  */ladruno_build_test/OpenSees|*/ladruno_build_test/OpenSees/)
    echo "refusing to run from ~/ladruno_build_test/OpenSees: pass a fresh clone" >&2; exit 2 ;;
esac
[ -f "$CLONE/Ladruno_scripts/mpm_bench/bench.py" ] || { echo "no bench.py under $CLONE" >&2; exit 2; }
[ -f "$PYD" ] && [ -f "$EXE" ] || { echo "PYD/EXE missing" >&2; exit 2; }

# the machine and the toolchain, recorded next to the results
mkdir -p "$OUT"
{
  echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "slurm job: ${SLURM_JOB_ID:-none} on ${SLURM_JOB_NODELIST:-?}"
  lscpu; echo; free -h; echo; echo "clone head: $(git -C "$CLONE" rev-parse HEAD)"
  echo "gcc: $(gcc --version | head -1)"; echo "MKLROOT=${MKLROOT:-unset}"
} > "$OUT/machine.txt" 2>&1

# one thread per core, placed once; the harness sets OMP_NUM_THREADS / MKL_NUM_THREADS per run
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export MKL_DYNAMIC=FALSE
export LADRUNO_OPENSEES_QUIET=1
if [ -n "${WP253_VENV:-}" ]; then export PATH="$WP253_VENV/bin:$PATH"; fi
if [ -z "${TCL_LIBRARY:-}" ]; then
  TCL_INIT=$(find "$HOME/.conan2" -name init.tcl -path '*tcl8.6*' 2>/dev/null | head -1)
  [ -n "$TCL_INIT" ] && export TCL_LIBRARY=$(dirname "$TCL_INIT")
fi

BENCH="$CLONE/Ladruno_scripts/mpm_bench/bench.py"
COMMON=(--site esmeralda --pyd "$PYD" --exe "$EXE" --affinity none)

run() {   # label, extra args...
  local label=$1; shift
  srun --exclude=node1,node2,node12,node13 --mem-per-cpu=${SLURM_MEM_PER_CPU:-2048}M --exclusive --cpus-per-task=24 \
       --cpu-bind=none python3.12 "$BENCH" "${COMMON[@]}" --out "$OUT/$label" "$@"
}

case "$SET" in
  big|all) run big --modes F --N 1000000 --threads 1,4,8,16,24 ;;
esac
case "$SET" in
  t24|all) run t24 --modes F --N 10000,100000 --threads 16,24 ;;
esac

echo "done: $OUT"
