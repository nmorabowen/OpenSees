"""WP-171 -- a refused `integrator` is an ERROR and keeps the previous integrator.

openseespy: `OPS_Integrator()` returned 0 even when the factory returned null,
so `Py_ops_integrator` raised nothing; the previous integrator silently stayed
in force. ADR-76 had fixed exactly this for `OPS_Algorithm` only.

classic Tcl (`specifyIntegrator`, commands.cpp): ~55 branches call
`theXAnalysis->setIntegrator(*theXIntegrator)` with no null check, so a refused
integrator SEGFAULTED when an analysis existed; without an analysis it returned
TCL_OK and nulled the global pointer, silently discarding the previous
integrator (`analysis` then fell back to a default).

Measured on the pre-fix build: every openseespy refusal case fails with
`DID NOT RAISE`; Tcl `transient-with-analysis` segfaults (0xC0000005) and
`transient-no-analysis` returns TCL_OK and falls back to the default
integrator. Tcl `static-no-analysis` is a premise case (stock LoadControl
already refused correctly) and only lacked the new message.

OBSERVABLE "the previous integrator was kept": a 1-DOF linear static model.
`LoadControl 0.1` then one step -> t = 0.1. Losing it to the `analysis Static`
default (`LoadControl 1.0`) gives t = 1.0.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DT_KEPT = 0.1


def _static_model():
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
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 10, 0)
    ops.algorithm("Newton")


def _dynamic_model():
    _static_model()
    ops.mass(2, 1.0)


# refused calls: (args, message the factory prints)
REFUSED_STATIC = [
    pytest.param(("LoadControl",), id="LoadControl-no-args"),
    pytest.param(("NoSuchIntegrator", 1.0), id="unknown-type"),
]
REFUSED_TRANSIENT = [
    pytest.param(("ExplicitBathe", 0.54, "-lnvd", 1.5), id="ExplicitBathe-lnvd-out-of-range"),
    pytest.param(("Newmark",), id="Newmark-no-args"),
]


@pytest.mark.parametrize("refused", REFUSED_STATIC)
@pytest.mark.parametrize("analysis_first", [False, True], ids=["no-analysis", "with-analysis"])
def test_refused_static_integrator_raises_and_keeps_previous(refused, analysis_first, capfd):
    _static_model()
    ops.integrator("LoadControl", DT_KEPT)
    if analysis_first:
        ops.analysis("Static")
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.integrator(*refused)
    assert "previous integrator left unchanged" in capfd.readouterr().err
    if not analysis_first:
        ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.getTime() == pytest.approx(DT_KEPT, rel=1e-12), (
        "the previous LoadControl was not kept")


@pytest.mark.parametrize("refused", REFUSED_TRANSIENT)
@pytest.mark.parametrize("analysis_first", [False, True], ids=["no-analysis", "with-analysis"])
def test_refused_transient_integrator_raises_and_keeps_previous(refused, analysis_first, capfd):
    _dynamic_model()
    ops.integrator("Newmark", 0.5, 0.25)
    if analysis_first:
        ops.analysis("Transient")
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.integrator(*refused)
    assert "previous integrator left unchanged" in capfd.readouterr().err
    if not analysis_first:
        ops.analysis("Transient")
    assert ops.analyze(3, 0.01) == 0
    assert ops.getTime() == pytest.approx(0.03, rel=1e-12)


def test_accepted_integrator_still_replaces(capfd):
    """Premise: the guard must not block a good integrator."""
    _static_model()
    ops.integrator("LoadControl", 0.5)
    ops.analysis("Static")
    ops.integrator("LoadControl", DT_KEPT)
    assert ops.analyze(1) == 0
    assert ops.getTime() == pytest.approx(DT_KEPT, rel=1e-12)


# --- classic Tcl ---------------------------------------------------------------
REPO = Path(__file__).resolve().parents[1]


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
constraints Plain; numberer Plain; system BandGeneral
test NormDispIncr 1e-10 10 0; algorithm Newton
"""

TCL_CASES = {
    # refused TRANSIENT integrator, no analysis yet: the pre-fix branch nulled
    # the global pointer and returned TCL_OK, so `analysis Transient` fell back
    # to its default and said so ("no Integrator specified").
    "transient-no-analysis": """
integrator Newmark 0.5 0.25
puts "RC [catch {integrator ExplicitBathe 0.54 -lnvd 1.5} m]"
analysis Transient
analyze 3 0.01
puts "TIME [getTime]"
""",
    # refused STATIC integrator whose stock branch already returned TCL_ERROR
    # (premise: the wrapper must not change an already-correct refusal)
    "static-no-analysis": """
integrator LoadControl 0.1
puts "RC [catch {integrator LoadControl} m]"
analysis Static
analyze 1
puts "TIME [getTime]"
""",
    # refused TRANSIENT integrator while an analysis exists: used to SEGFAULT
    "transient-with-analysis": """
integrator Newmark 0.5 0.25
analysis Transient
puts "RC [catch {integrator ExplicitBathe 0.54 -lnvd 1.5} m]"
analyze 3 0.01
puts "TIME [getTime]"
""",
}
TCL_EXPECTED_TIME = {"transient-no-analysis": 0.03, "static-no-analysis": 0.1,
                     "transient-with-analysis": 0.03}


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
@pytest.mark.parametrize("case", list(TCL_CASES))
def test_classic_tcl_refused_integrator(case, tmp_path):
    deck = tmp_path / f"{case}.tcl"
    deck.write_text(TCL_MODEL + TCL_CASES[case] + 'puts "DECK-DONE"\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see tests/test_wp103_getstringfromall_tcl.py (WinError 6).
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=120, env=env, cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, f"exit {proc.returncode} (139/-11/0xC0000005 = segfault):\n{out}"
    assert "DECK-DONE" in out, out
    assert "RC 1" in out, "refused integrator must be a Tcl error:\n" + out
    assert "previous integrator left unchanged" in out, out
    assert "no Integrator specified" not in out, "previous integrator was discarded:\n" + out
    t = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("TIME ")]
    assert t == pytest.approx([TCL_EXPECTED_TIME[case]], rel=1e-12), out
