#!/usr/bin/env bash
# ADR-77 V0 — /arch:AVX2 probe. Two questions, two measurements:
#
#   SPEED: does AVX2 move the element kernels? Metric = elem.tangent ms from
#   the deep profile (thread-invariant to 4T per T0b, so 4T is fair), plus step
#   wall. A/B INTERLEAVED across 3 rounds in separate processes -- the T3 noise
#   study measured +/-12% wall run-to-run, so uninterleaved wall comparisons of
#   single runs are worthless.
#
#   DRIFT: does FMA contraction change answers? Metric = uz at 1 THREAD, where
#   PARDISO is bit-reproducible (P1f), so ANY difference is the codegen, not
#   the solver. This is the number the exact-QA policy decision needs.
#
# Baseline pyd = dist/bin (no /arch). AVX2 pyd = dist/bin_avx2 (staged by the
# probe driver from the ninja build; dist/bin is never overwritten -- build.bat
# staging is deliberately bypassed).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
WT="C:/Users/nmb/Documents/Github/OpenSees/.claude/worktrees/pardisio-profiling-0a03b1"
BASE="$WT/dist/bin"
AVX2="$WT/dist/bin_avx2"
N="${1:-15}"

echo "== SPEED: 3 interleaved rounds, 4 threads, n=$N, PARDISO only =="
for r in 1 2 3; do
  for arm in base avx2; do
    d=$BASE; [ "$arm" = "avx2" ] && d=$AVX2
    OPS_PYD="$d" MKL_NUM_THREADS=4 OMP_NUM_THREADS=4 \
    OPS_T0_N="$N" OPS_T0_STEPS=4 OPS_T0_REPEATS=1 OPS_T0_SOLVERS=Pardiso \
    OPS_T0_TAG="_v0${arm}_r${r}" \
    python3.12 -S "$HERE/t0_step_anatomy.py" 2>&1 | grep -E "Pardiso.*wall" \
      | sed "s/^/  r$r $arm /"
  done
done

echo "== DRIFT: 1 thread (deterministic), both builds, same deck =="
for arm in base avx2; do
  d=$BASE; [ "$arm" = "avx2" ] && d=$AVX2
  OPS_PYD="$d" MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  OPS_T0_N=10 OPS_T0_STEPS=4 OPS_T0_REPEATS=1 OPS_T0_SOLVERS=Pardiso \
  OPS_T0_TAG="_v0drift_${arm}" \
  python3.12 -S "$HERE/t0_step_anatomy.py" 2>&1 | grep -E "Pardiso" \
    | sed "s/^/  $arm /"
done
echo "done — report with t0-style h5 readout"
