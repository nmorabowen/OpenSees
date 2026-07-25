#!/bin/bash
# ADR-75 P1d gate: does SYMMETRIC PARDISO (upper-triangle half-storage, mtype
# +/-2) actually beat the shipped UNSYMMETRIC PARDISO (mtype 11) on Lane B?
#
# This is the measure-don't-assume gate ADR-75 §3 demands. The only symmetric
# solver previously measurable in this fork (SparseSYM) was 2.10x SLOWER than
# unsymmetric UmfPack, so "symmetric halves the work" is a HYPOTHESIS here.
#
# UmfPack stays in the set as the anchor every ADR-75 number is quoted against,
# and as the correctness oracle: Lane B's LadrunoJ2 tangent is symmetric, so a
# correct half-storage path MUST reproduce its tip displacement.
#
# PYD deliberately points at THIS worktree (memory: never bench the shared
# checkout -- its branch moves under you).
PYD="C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr-75-session-handoff-29effa\dist\bin"
PY=/c/Users/nmb/AppData/Local/Python/bin/python3.12
for T in 1 4; do
  echo "########## MKL_NUM_THREADS=$T ##########"
  OPS_PYD="$PYD" MKL_NUM_THREADS=$T OMP_NUM_THREADS=$T LADRUNO_OPENSEES_QUIET=1 \
    OPS_BENCH_SOLVERS="UmfPack,Pardiso,PardisoSym,PardisoSPD" OPS_BENCH_REPEATS=5 \
    "$PY" -S laneB_solver_bench.py > "p1d_run_t${T}.log" 2>&1
  # Keep the FULL log: the grep-only form hid a silent ProfileSPDLinSOE fallback
  # (a mistyped option left theSOE null) and the run looked merely "slow".
  grep -E "probe |^UmfPack|^Pardiso|^solver|WARNING" "p1d_run_t${T}.log"
  # cp, NOT mv -- mv destroyed the locked baseline JSON twice in the P1 session.
  cp -f laneB_solver_baseline.json "laneB_p1d_threads${T}.json" 2>/dev/null
done
echo "SWEEP-DONE"
