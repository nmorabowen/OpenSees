#!/usr/bin/env bash
# ADR-77 T0b — does the transient step-anatomy flip with THREAD COUNT?
#
# T0 measured at 1 thread (chosen so the run is a bit-exact oracle: threaded
# PARDISO is not byte-reproducible, ADR-75 P1f) and found linearSolve STILL the
# largest loop under PARDISO (45.2% vs elem.tangent 40.8%) -- unlike the L3-0
# static lane's 74.9% element share. That comparison is not apples-to-apples on
# thread count, so no assembly-side conclusion may be drawn from T0 alone (G0).
#
# This sweeps MKL_NUM_THREADS. Each thread count is a SEPARATE PROCESS because
# MKL reads MKL_NUM_THREADS at library init -- an in-process loop would silently
# measure the first setting for every arm.
#
# UmfPack is swept too: it is a serial factorization, so its curve is the
# control. If UmfPack's split also moves with threads, the mover is BLAS inside
# the element kernels, not the solver, and the whole comparison needs rethinking.
#
# Usage:  ./t0b_thread_sweep.sh [n] [steps] [repeats]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
WT="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/pardisio-profiling-0a03b1"
N="${1:-15}"
STEPS="${2:-4}"
REPEATS="${3:-2}"
LADDER="${OPS_T0B_THREADS:-1 2 4 8 16}"

echo "ADR-77 T0b thread sweep | n=$N steps=$STEPS repeats=$REPEATS ladder=[$LADDER]"
for t in $LADDER; do
  echo "--- MKL_NUM_THREADS=$t ---"
  OPS_PYD="$WT/dist/bin" \
  MKL_NUM_THREADS="$t" OMP_NUM_THREADS="$t" \
  OPS_T0_N="$N" OPS_T0_STEPS="$STEPS" OPS_T0_REPEATS="$REPEATS" \
  OPS_T0_TAG="_t$t" \
  python3.12 -S "$HERE/t0_step_anatomy.py" 2>&1 | grep -vE "^WARNING|^$"
done
echo "done. Report with:  <venv>/python.exe t0b_report.py $N"
