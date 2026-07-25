#!/usr/bin/env bash
# ADR-75b L3-0 -- the Lane-3 gate sweep.
#
# Measures, per lane, the fraction of step wall inside each of the three threadable
# element loops (ADR-75b §2), split into KERNEL (threadable) and SCATTER (not), using
# the elem.update / elem.tangent / elem.residual scopes that already ship.
#
# Lane B is run on three solver configs INTERLEAVED, because the point is not the
# absolute wall but whether the Lane-1 solver win pushes the element fraction across
# ADR-75 §3's >40% gate. Interleaving is mandatory: this box swings +-30% on background
# load (phase0/PHASE0_RECIPE.md).
#
# Usage:  ./sweep_l3a.sh [rounds]        (default 3)
# Env:    OPS_PYD   -- dist/bin to test (default: this worktree's)
#         L3_LANES  -- space-separated subset of "B A D" (default: all)
set -u

ROUNDS="${1:-3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE0="$(cd "$HERE/../phase0" && pwd)"
WT="$(cd "$HERE/../../../.." && pwd)"

export OPS_PYD="${OPS_PYD:-$(cygpath -w "$WT/dist/bin" 2>/dev/null || echo "$WT/dist/bin")}"
export LADRUNO_OPENSEES_QUIET=1
LANES="${L3_LANES:-B A D}"

# python3.12 runs OpenSees (no h5py); opensees_env parses the h5 (has h5py). Two
# interpreters on purpose -- do not "simplify" this to one.
PY_RUN="python3.12"
PY_H5="/c/Users/nmb/venv/opensees_env/Scripts/python.exe"
SITE='C:\Users\nmb\venv\opensees_env\Lib\site-packages'

echo "=== ADR-75b L3-0 gate sweep ==="
echo "pyd    : $OPS_PYD"
echo "rounds : $ROUNDS   lanes: $LANES"
# Recipe caveat: a foreign compile eating 8 cores invalidates every wall in the sweep.
BUSY="$(tasklist 2>/dev/null | grep -icE '\b(cl|ninja|link)\.exe' || true)"
echo "busy compilers detected: ${BUSY:-0}  (must be 0 for publishable walls)"
echo

run_one() {   # $1=tag  $2=script  $3=OPS_SYSTEM  $4=MKL_NUM_THREADS  $5=round
  local tag="$1" script="$2" sysopt="$3" thr="$4" r="$5"
  local h5="$HERE/l3a_${tag}_r${r}.h5"
  rm -f "$h5"
  OPS_PROF_OUT="$(cygpath -w "$h5" 2>/dev/null || echo "$h5")" \
  OPS_SYSTEM="$sysopt" OMP_NUM_THREADS="$thr" MKL_NUM_THREADS="$thr" \
    "$PY_RUN" -S "$PHASE0/$script" > "$HERE/l3a_${tag}_r${r}.log" 2>&1
  local rc=$?
  # NEVER grep the log away -- a silent ProfileSPDLinSOE fallback hides in the lines
  # you would have dropped (banked P1d trap). The full log is kept; we only PEEK here.
  #
  # Adversarial review (N6): the ORIGINAL guard here matched "unknown system|WARNING.*fall|FAILED",
  # and VERIFIED BY EXECUTION it matches NONE of the three real fallback sites -- OpenSees emits
  #   "WARNING analysis Static - no LinearSOE specified, " / " ProfileSPDLinSOE default will be used"
  # (OpenSeesCommands.cpp:689-690, :863-864, :941-942), which contains none of those tokens. The one
  # guard written against the banked trap fired on ZERO of its targets. It also fired SPURIOUSLY on
  # every Lane-A run ("FAILED at step 86"), training the operator to ignore it -- the worst outcome.
  local log="$HERE/l3a_${tag}_r${r}.log"
  if grep -qiE "ProfileSPDLinSOE default|no LinearSOE specified|unknown system" "$log"; then
      echo "    !! SILENT SOLVER FALLBACK for $tag r$r -- the requested system was NOT used"
  fi
  # Positive assertion beats pattern-matching a warning: the h5 records the solver that actually ran.
  if [ -n "$sysopt" ]; then
      local want; want="$(echo "$sysopt" | awk '{print $1}')"
      if ! grep -qE "^system: \['${want}'" "$log"; then
          echo "    !! $tag r$r did not report 'system: ['${want}'...]' -- check $log"
      fi
  fi
  local wall
  wall="$(grep -oE 'wall[= ]+[0-9.]+' "$HERE/l3a_${tag}_r${r}.log" | tail -1)"
  echo "    $tag r$r rc=$rc $wall"
  [ -f "$h5" ] || echo "    !! no h5 produced for $tag r$r"
}

for r in $(seq 1 "$ROUNDS"); do
  echo "--- round $r ---"
  case " $LANES " in *" B "*)
    run_one laneB_umf1   laneB_model.py      "UmfPack"               1 "$r"
    run_one laneB_pard1  laneB_model.py      "Pardiso -matrixType 2" 1 "$r"
    run_one laneB_pard4  laneB_model.py      "Pardiso -matrixType 2" 4 "$r"
    # UNSYMMETRIC PARDISO is not optional padding: without it, "PARDISO's SOE" and
    # "half-storage" are confounded, and the first draft of RESULTS_l3a_update_scope.md
    # credited the addA win to the wrong one. It turns out unsym PARDISO's addA is the
    # SLOWEST of the three (1.28x UmfPack). Keep this row. (adversarial review, finding 1)
    run_one laneBunsym   laneB_model.py      "Pardiso -matrixType 0" 4 "$r"
  ;; esac
  case " $LANES " in *" A "*) run_one laneA laneA_fiber_frame.py "UmfPack" 1 "$r" ;; esac
  case " $LANES " in *" D "*) run_one laneD laneD_model.py       ""        1 "$r" ;; esac
done

echo
echo "=== parse ==="
for f in "$HERE"/l3a_*.h5; do
  [ -e "$f" ] || continue
  PYTHONPATH="$SITE" "$PY_H5" -S "$HERE/parse_lane3.py" "$f" \
      --json "${f%.h5}.json" 2>&1 | sed 's/^/  /'
  echo
done
echo "done. h5 + json + logs in $HERE"
