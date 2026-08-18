"""`printA -sparse` flag order in the CLASSIC Tcl engine (zone_a).

The Tcl twin of `test_printa_unsized_soe.py::test_printa_sparse_flag_order_is_free`.
`-sparse <baseIndex>` is optional, but the read was unconditional in BOTH
engines, so `printA -sparse -ret` consumed `-ret` as the index and failed with
"failed to read -sparse <baseIndex>". Only `printA -ret -sparse 0` worked, which
makes the flags read as positional when they are not.

WHY A SEPARATE TCL GATE AT ALL. `SRC/tcl/commands.cpp` is a SECOND, independent
implementation of `printA` -- it is the one `OpenSees.exe`, `OpenSeesSP` and
`OpenSeesMP` run, and it does its own argv walk. Fixing the openseespy path in
`OpenSeesCommands.cpp` does nothing for it. This fork has a documented history of
exactly that drift (see the `#729` ledger row: `ladrunoArcLength`/`ladrunoDR`
existed only in the Python wrappers, so a `.tcl` deck could start a dynamic
relaxation but never read its settling gate; and the ADR-75 P2h row, where a Tcl
`system Mumps` ladder silently DROPPED `-BLR`/`-stats`). A Python-only gate would
have let the Tcl half stay broken and reported green.

The classic engine's version was also slightly worse than the Python one:
`currentArg++` ran BEFORE the bounds test, so a bare trailing `-sparse` survived
only by luck -- the cursor landed on argc and the loop ended. The deck covers
that form too.

Skips when the exe has not been built (`build.bat OpenSees`), rather than
failing: the pyd-only build is a legitimate dev state.
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "printa_sparse_flag_order.tcl"


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(
    TCL_EXE is None,
    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)",
)


@needs_exe
def test_tcl_printa_sparse_accepts_any_flag_order():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL is required, not tidiness: under pytest's capture sys.stdin has
    # no inheritable Windows handle, and letting the child inherit it raises
    # `OSError: [WinError 6] The handle is invalid` before the exe ever starts.
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO),
                          stdin=subprocess.DEVNULL, errors="replace")
    out = proc.stdout + proc.stderr

    # Non-vacuity FIRST: a deck that died on line 1 also emits no FAIL lines, so
    # "no failures" is not evidence of anything until the verdict line is present.
    assert "SELF-TEST:" in out, f"deck did not reach its verdict:\n{out}"
    assert out.count("PASS ") + out.count("FAIL ") == 4, (
        f"expected 4 checks, got {out.count('PASS ')} pass / {out.count('FAIL ')} "
        f"fail -- the deck did not run to completion:\n{out}"
    )
    assert "SELF-TEST: 0 failure(s)" in out, (
        "a printA -sparse flag ordering was rejected by the classic Tcl engine "
        "-- the optional base index is eating the next flag again:\n" + out
    )
