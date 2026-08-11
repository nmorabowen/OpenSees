"""Classic-Tcl registration gate for `ladrunoDR` / `ladrunoArcLength` (#729).

The defect this guards is invisible to every other test in the suite: both
commands worked from openseespy and from the TclWrapper engine, and were simply
never registered in `SRC/tcl/commands.cpp` — the engine `OpenSees.exe`,
`OpenSeesSP` and `OpenSeesMP` actually run. So a `.tcl` deck could SET
`integrator LadrunoDynamicRelaxation` and then die on

    invalid command name "ladrunoDR"

with no way to read the settling gate that says when the relaxation is done. A
Python test cannot reach that code path, so this one shells out to the exe.

It also guards the *fix*, which has a silent failure mode of its own: bridging to
the no-arg `OPS_Ladruno*Cmd()` forms (the shape the contact family uses, #726)
compiles and links, but those read the `cmds` singleton that only the DL command
engine constructs — under `OpenSees.exe` they hit `cmds == 0`, return 0 (SUCCESS)
and write nothing, leaving a command that exists and answers an empty string. The
deck therefore fails on an EMPTY result, not merely on a Tcl error.

Deck: tests/tcl/ladruno_solver_queries.tcl (also runnable by hand).
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "ladruno_solver_queries.tcl"


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
def test_classic_tcl_answers_the_solver_queries():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL is required, not tidiness: under pytest's capture, sys.stdin has
    # no inheritable Windows handle, and letting the child inherit it raises
    # `OSError: [WinError 6] The handle is invalid` before the exe ever starts.
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr

    # non-vacuity: the deck must have actually run its checks. A deck that died on
    # line 1 also produces no FAIL lines, so absence of failure is not evidence.
    assert "SELF-TEST:" in out, f"deck did not reach its verdict:\n{out}"
    n_pass = out.count("PASS ")
    assert n_pass >= 22, f"only {n_pass} checks ran; expected the full deck:\n{out}"

    assert "invalid command name" not in out, (
        "a command is still unregistered in the classic Tcl engine:\n" + out
    )
    assert "EMPTY RESULT" not in out, (
        "a command answered nothing — this is the cmds==0 silent path, i.e. the "
        "bridge is calling the no-arg OPS_Ladruno*Cmd() form instead of ...CmdOn:\n"
        + out
    )
    assert "SELF-TEST: PASS" in out and proc.returncode == 0, (
        f"classic-Tcl solver-query deck failed (rc={proc.returncode}):\n{out}"
    )
