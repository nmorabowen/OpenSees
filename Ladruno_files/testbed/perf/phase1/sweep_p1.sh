#!/bin/bash
# ADR-75 P1 gate: Lane-B UmfPack (reference) vs threaded PARDISO, MKL thread sweep.
# Only these two solvers: SuperLU/SparseSYM are thread-invariant and already
# measured in RESULTS_laneB_baseline.md, so re-running them per thread count
# would triple the wall for no information.
PYD="C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\mumps-opensees-study-f833bf\dist\bin"
PY=/c/Users/nmb/AppData/Local/Python/bin/python3.12
for T in 1 2 4 8; do
  echo "########## MKL_NUM_THREADS=$T ##########"
  OPS_PYD="$PYD" MKL_NUM_THREADS=$T OMP_NUM_THREADS=$T LADRUNO_OPENSEES_QUIET=1 \
    OPS_BENCH_SOLVERS="UmfPack,Pardiso" OPS_BENCH_REPEATS=5 \
    "$PY" -S laneB_solver_bench.py 2>&1 | grep -E "probe |^UmfPack|^Pardiso|^solver"
  mv -f laneB_solver_baseline.json "laneB_p1_threads${T}.json" 2>/dev/null
done
echo "SWEEP-DONE"
