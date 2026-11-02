"""WP-221: classic-Tcl (OpenSees.exe) gate for `ladrunoSymmetryGuard`.

Shells out to the exe, as test_ladruno_solver_queries_tcl.py does, because the
classic engine has its own `system` parser and its own `theSOE`, neither of which
a Python test can reach. Deck: tests/tcl/wp221_symmetry_guard.tcl.
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "wp221_symmetry_guard.tcl"


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_symmetry_guard():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see test_ladruno_solver_queries_tcl.py (WinError 6 under capture)
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert "SELF-TEST:" in out, f"deck did not reach its verdict:\n{out}"
    # 2 header checks + 5 per system; ProfileSPD x2 always, Pardiso x2 with MKL
    assert out.count("PASS ") >= 12, f"too few checks ran:\n{out}"
    assert "SELF-TEST: PASS" in out and proc.returncode == 0, (
        f"classic-Tcl symmetry-guard deck failed (rc={proc.returncode}):\n{out}")
