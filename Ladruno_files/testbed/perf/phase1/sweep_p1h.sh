#!/usr/bin/env bash
# ADR-75 P1h -- the PARDISO solver-phase split sweep.
#
# Closes the measurement half of the instrumentation item: PR #667 added the brackets and
# CI proves they EXIST, but nothing published what the split IS. That split is the number
# that bounds every reuse lever in the lane -- the `factored` gate, `-krylov`,
# ModifiedNewton and IMPL-EX all attack phase 22 (soe.factor) and do nothing for phase 33
# (soe.trisolve) -- and it is how you see where the 1->8 thread ceiling of 1.58x dies.
#
# COARSE by default, not deep. The headline split needs only coarse scopes (soe.factor
# and friends are OPS_PROFILE_SCOPE), and deep costs ~0.5 us per scope instance on an
# element loop -- a tax that lands on the assembly side and would skew the very fraction
# being reported. One extra DEEP round is run at the end for `soe.addA`, which is
# deep-gated by design.
#
# Interleaved rounds are mandatory: this box swings +-30% on background load
# (phase0/PHASE0_RECIPE.md).
#
# Usage:  ./sweep_p1h.sh [rounds]        (default 3)
# Env:    OPS_PYD -- dist/bin to test (default: this worktree's)
set -u

ROUNDS="${1:-3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE0="$(cd "$HERE/../phase0" && pwd)"
WT="$(cd "$HERE/../../../.." && pwd)"

export OPS_PYD="${OPS_PYD:-$(cygpath -w "$WT/dist/bin" 2>/dev/null || echo "$WT/dist/bin")}"
export LADRUNO_OPENSEES_QUIET=1

# python3.12 runs OpenSees (no h5py); opensees_env parses the h5 (has h5py). Two
# interpreters on purpose -- do not "simplify" this to one. (lane3/sweep_l3a.sh)
PY_RUN="python3.12"
PY_H5="/c/Users/nmb/venv/opensees_env/Scripts/python.exe"
SITE='C:\Users\nmb\venv\opensees_env\Lib\site-packages'

echo "=== ADR-75 P1h -- PARDISO phase-split sweep ==="
echo "pyd    : $OPS_PYD"
echo "rounds : $ROUNDS"
BUSY="$(tasklist 2>/dev/null | grep -icE '\b(cl|ninja|link)\.exe' || true)"
echo "busy compilers detected: ${BUSY:-0}  (must be 0 for publishable walls)"
if [ "${BUSY:-0}" != "0" ]; then
  echo "!! ABORT: a compile is running; every wall in this sweep would be invalid."
  exit 1
fi
# The brackets are the whole point -- refuse to produce a table of zeros from a stale pyd.
if ! grep -q "soe.cgs" "$WT/dist/bin/opensees.pyd" 2>/dev/null; then
  echo "!! ABORT: $OPS_PYD/opensees.pyd has no 'soe.cgs' string => it predates PR #667."
  echo "   Rebuild (Ladruno_scripts\\build.bat OpenSeesPy) before measuring."
  exit 1
fi
echo

run_one() {   # $1=tag  $2=OPS_SYSTEM  $3=MKL_NUM_THREADS  $4=round  $5=prof flags
  local tag="$1" sysopt="$2" thr="$3" r="$4" flags="$5"
  local h5="$HERE/p1h_${tag}_r${r}.h5"
  local log="$HERE/p1h_${tag}_r${r}.log"
  rm -f "$h5"
  OPS_PROF_OUT="$(cygpath -w "$h5" 2>/dev/null || echo "$h5")" \
  OPS_SYSTEM="$sysopt" OPS_PROF_FLAGS="$flags" \
  OMP_NUM_THREADS="$thr" MKL_NUM_THREADS="$thr" \
    "$PY_RUN" -S "$PHASE0/laneB_model.py" > "$log" 2>&1
  local rc=$?
  # NEVER grep the log away -- a silent ProfileSPDLinSOE fallback hides in the lines you
  # would have dropped (banked P1d trap). Match the REAL fallback strings
  # (OpenSeesCommands.cpp:689-690), not an invented pattern (banked L3-0 N6 finding).
  if grep -qiE "ProfileSPDLinSOE default|no LinearSOE specified|unknown system" "$log"; then
      echo "    !! SILENT SOLVER FALLBACK for $tag r$r -- the requested system was NOT used"
  fi
  # Positive assertion beats pattern-matching a warning.
  local want; want="$(echo "$sysopt" | awk '{print $1}')"
  if ! grep -qE "^system: \['${want}'" "$log"; then
      echo "    !! $tag r$r did not report 'system: ['${want}'...]' -- check $log"
  fi
  local wall; wall="$(grep -oE 'wall=[0-9.]+' "$log" | tail -1)"
  echo "    $tag r$r rc=$rc $wall  (MKL=$thr, flags='${flags:-coarse}')"
  [ -f "$h5" ] || echo "    !! no h5 produced for $tag r$r"
}

# ---- COARSE rounds: the headline split -------------------------------------------
# UmfPack is the reference row, not padding: it has carried soe.symbolic/factor/trisolve
# since ADR-40, so it validates the parser AND is the cross-solver comparison the shared
# scope names were chosen to enable.
for r in $(seq 1 "$ROUNDS"); do
  echo "--- coarse round $r ---"
  run_one umf1     "UmfPack"                          1 "$r" ""
  run_one unsym1   "Pardiso -matrixType 0"            1 "$r" ""
  run_one unsym4   "Pardiso -matrixType 0"            4 "$r" ""
  run_one sym4     "Pardiso -matrixType 2"            4 "$r" ""
  # -krylov on mtype 11: the CGS path is REFUSED for -matrixType 2 (mtype -2 has no
  # documented Intel CGS mode, P1e), so exercising soe.cgs requires 0 or 1.
  run_one krylov4  "Pardiso -matrixType 0 -krylov 6"  4 "$r" ""
done

# ---- one DEEP round: soe.addA per-element scatter attribution ---------------------
echo "--- deep round (soe.addA is deep-gated; expect a wall tax) ---"
run_one deep_unsym4 "Pardiso -matrixType 0" 4 1 "-deep"
run_one deep_sym4   "Pardiso -matrixType 2" 4 1 "-deep"

echo
echo "=== parse ==="
for f in "$HERE"/p1h_*.h5; do
  [ -e "$f" ] || continue
  PYTHONPATH="$SITE" "$PY_H5" -S "$HERE/p1h_parse_split.py" "$f" \
      --json "${f%.h5}.json" 2>&1 | sed 's/^/  /'
done
echo "done. h5 + json + logs in $HERE"
