"""Classic-Tcl gate for WP-104: `wipe` zeroes LadrunoSANISAND's process-wide
IMPL-EX ledger under `OpenSees.exe` as well as under openseespy.

The openseespy half is `test_wp104_implex_refusals_wipe_reset.py`.  This one
exists because `wipe` in `SRC/tcl/commands.cpp` and `wipe` in
`SRC/interpreter/OpenSeesCommands.cpp` are two code paths; the fix lives in
`OPS_clearAllNDMaterial()` (`SRC/material/nD/NDMaterial.cpp`), which both
call, and this deck proves the classic-Tcl one really does reach it.  (Same
shape as `tests/test_wp103_getstringfromall_tcl.py`.)

Deck: tests/tcl/wp104_implex_refusals_wipe.tcl (also runnable by hand).
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "wp104_implex_refusals_wipe.tcl"


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
def test_wipe_zeroes_the_implex_ledger_under_classic_tcl():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL + utf-8/replace: the two Windows capture traps recorded in
    # tests/_testbed/subprocess_run.py -- without them a passing deck can look
    # like a crashed one under `pytest tests/`.
    proc = subprocess.run(
        [TCL_EXE, str(DECK)],
        cwd=str(DECK.parent),
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
        timeout=300,
    )
    out = proc.stdout

    assert "inherited the previous model's refusal ledger" not in out, (
        "a FRESH LadrunoSANISAND after `wipe` still reads the previous model's "
        "implexRefusals under OpenSees.exe -- the Tcl wipe path did not reach "
        "OPS_clearAllNDMaterial()'s reset:\n" + out
    )
    assert "SELF-TEST: PASS" in out and proc.returncode == 0, (
        f"WP-104 classic-Tcl wipe deck failed (rc={proc.returncode}):\n{out}"
    )
