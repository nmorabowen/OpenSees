"""WP-173 -- an unrecognised `system` type no longer frees the live SOE.

Measured 2026-10-06 on ladruno f4a14761e (Windows, dist\\bin\\OpenSees.exe):
`system FullGeneral; ...; analysis Static; system Mumps; analyze 1` killed the
process with exit 127 and no message.

Root cause (classic Tcl `specifySOE`, SRC/tcl/commands.cpp): serial OpenSees.exe
is built WITHOUT `_MUMPS` (CMake adds `${MUMPS_FLAG}` only to the SP/MP/openseesmp
targets), so `Mumps` matches no branch -- exactly like a typo. The global
`theSOE` is not cleared on entry, so it still pointed at the SOE the analysis
owns, and the tail handed that SAME object to `setLinearSOE()`, whose first act
is `delete theSOE`. The analysis kept a dangling pointer -> crash at `analyze`.
With no analysis yet the same fall-through returned TCL_OK silently, and the
deck ran on the previous system while believing it had switched.

The FIX is WP-172's `specifySOE` wrapper (clear the global, restore it when
nothing was created, return TCL_ERROR). WP-173 is the regression gate for the
measured `Mumps` deck: on ladruno f4a14761e the `*-with-analysis` cases die
(exit 127) and `no-analysis` returns RC 0 with no warning; on the WP-172 build
all cases pass. openseespy is safe by construction (`OPS_System` starts from a
local `theSOE = 0` and returns -1 on an unknown type); the parity test pins it.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
U_EXPECTED = 0.1   # 1-DOF truss, k = 1, LoadControl 0.1, one step


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()

TCL_MODEL = """
model basic -ndm 1 -ndf 1
node 1 0.0; node 2 1.0; fix 1 1; mass 2 1.0
uniaxialMaterial Elastic 1 1.0; element Truss 1 1 2 1.0 1
timeSeries Linear 1; pattern Plain 1 1 { load 2 1.0 }
constraints Transformation; numberer Plain; system FullGeneral
test FixedNumIter 3; algorithm Newton
"""

STATIC_RUN = """
puts "AN [analyze 1]"
puts "U [nodeDisp 2 1]"
"""

# name -> (deck tail, must_refuse). must_refuse=None: `Mumps` may legitimately
# be compiled into some future serial build; then it must simply work.
TCL_CASES = {
    # the measured reproducer
    "mumps-static-with-analysis": ("""
integrator LoadControl 0.1
analysis Static
puts "RC [catch {system Mumps} m]"
""" + STATIC_RUN, None),
    "typo-static-with-analysis": ("""
integrator LoadControl 0.1
analysis Static
puts "RC [catch {system NoSuchSystem} m]"
""" + STATIC_RUN, True),
    # DirectIntegrationAnalysis::setLinearSOE has the same delete-first shape
    "typo-transient-with-analysis": ("""
integrator Newmark 0.5 0.25
analysis Transient
puts "RC [catch {system NoSuchSystem} m]"
puts "AN [analyze 3 0.01]"
puts "TIME [getTime]"
""", True),
    # no analysis yet: pre-fix this returned TCL_OK with no warning
    "typo-no-analysis": ("""
puts "RC [catch {system NoSuchSystem} m]"
integrator LoadControl 0.1
analysis Static
""" + STATIC_RUN, True),
    # premise: a recognised system still replaces the live one
    "accepted-swap-with-analysis": ("""
integrator LoadControl 0.1
analysis Static
puts "RC [catch {system BandGeneral} m]"
""" + STATIC_RUN, False),
}


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
@pytest.mark.parametrize("case", list(TCL_CASES))
def test_classic_tcl_unknown_system_keeps_live_soe(case, tmp_path):
    tail, must_refuse = TCL_CASES[case]
    deck = tmp_path / f"{case}.tcl"
    deck.write_text(TCL_MODEL + tail + 'puts "DECK-DONE"\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see tests/test_wp103_getstringfromall_tcl.py (WinError 6).
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=120, env=env, cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, (
        f"exit {proc.returncode} (127/139/0xC0000005 = the freed SOE was used):\n{out}")
    assert "DECK-DONE" in out, out
    assert "AN 0" in out, out

    refused = "RC 1" in out
    if must_refuse is None:
        must_refuse = refused   # Mumps: either outcome, but consistent below
    if must_refuse:
        assert refused, "unrecognised system must be a Tcl error:\n" + out
        assert ("is unknown or not installed" in out
                or "previous system left unchanged" in out), out
    else:
        assert "RC 0" in out, out

    if "TIME " in out:
        t = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("TIME ")]
        assert t == pytest.approx([0.03], rel=1e-12), out
    else:
        u = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("U ")]
        assert u == pytest.approx([U_EXPECTED], rel=1e-12), out


def test_openseespy_unknown_system_raises_and_keeps_soe():
    """Parity: the Python ladder already refused; pin it."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("FixedNumIter", 3)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", U_EXPECTED)
    ops.analysis("Static")
    with pytest.raises(Exception):
        ops.system("NoSuchSystem")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(2, 1) == pytest.approx(U_EXPECTED, rel=1e-12)
