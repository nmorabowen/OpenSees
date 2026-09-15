"""Classic-Tcl gate for OPS_GetStringFromAll's BUFFER contract (WP-103).

`elementAPI.h` documents OPS_GetStringFromAll as "does a strcpy", and the
openseespy backend (`PythonModule::getStringFromAll`) does. The classic-Tcl
backend did not -- it was

    extern "C" const char* OPS_GetStringFromAll(char *buffer, int len)
    { return OPS_GetString(); }          // never writes `buffer`

so the family idiom

    char tok[64];
    OPS_GetStringFromAll(tok, sizeof(tok));
    if (strcmp(tok, "auto") == 0) ...

read uninitialised stack under `OpenSees.exe`. Measured on the 9c2f964 release
build, before the fix::

    WARNING LadrunoKinematicCoupling: -k wants a number or 'auto', got '<garbage>'
    WARNING LadrunoEmbeddedNode: nHost must be >= 1 (or use -host eleTag); got '<garbage>'

WHY THIS TEST SHELLS OUT: the defect is classic-Tcl-ONLY. openseespy filled the
buffer, so the whole openseespy battery for these elements passed while the same
decks were unusable from `OpenSees.exe`. A Python test cannot reach the code
path -- only the exe can. (Same shape as
`tests/test_ladruno_solver_queries_tcl.py`, the #729 registration gate.)

Deck: tests/tcl/wp103_getstringfromall.tcl (also runnable by hand).
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "wp103_getstringfromall.tcl"


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
def test_getstringfromall_fills_the_buffer_under_classic_tcl():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL is required, not tidiness: under pytest's capture sys.stdin has
    # no inheritable Windows handle, and letting the child inherit it raises
    # `OSError: [WinError 6] The handle is invalid` before the exe ever starts.
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr

    # non-vacuity: a deck that died on line 1 also prints no FAIL lines, so the
    # absence of failure is not by itself evidence.
    assert "SELF-TEST:" in out, f"deck did not reach its verdict:\n{out}"
    n_pass = out.count("PASS ")
    assert n_pass >= 15, f"only {n_pass} checks ran; expected the full deck:\n{out}"

    # The signatures of the unwritten buffer, asserted BY NAME rather than by a
    # blanket "no element was refused": block A2 refuses an element on purpose
    # (`-k auto` with no -host), so a blanket check would fire on a healthy run.
    assert "got ''" not in out, (
        "an option token came back EMPTY -- OPS_GetStringFromAll is not filling "
        "the caller's buffer:\n" + out
    )
    assert "-dof needs at least one component" not in out, (
        "the -dof greedy reader strtol'd an unwritten buffer:\n" + out
    )
    assert "nHost must be >= 1" not in out, (
        "LadrunoEmbeddedNode's host spec came back as garbage -- the explicit "
        "<nHost> h1..hN form is unusable when the buffer is not filled:\n" + out
    )

    assert "SELF-TEST: PASS" in out and proc.returncode == 0, (
        f"WP-103 classic-Tcl buffer deck failed (rc={proc.returncode}):\n{out}"
    )
