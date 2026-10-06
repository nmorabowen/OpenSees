"""WP-172 -- a refused system/numberer/constraints/test/analysis is an ERROR
and keeps the previous object.

The ADR-76 (`OPS_Algorithm`) / WP-171 (`OPS_Integrator`, `specifyIntegrator`)
shape, audited for the remaining analysis-setup factory commands:

openseespy (`OPS_System`, `OPS_Numberer`, `OPS_ConstraintHandler`,
`OPS_CTest`): each called `cmds->setX(factory())` and returned 0 even when the
factory returned null. `Py_ops_*` raises only on `< 0`, so nothing was raised,
and `setX(0)` NULLED the global (deleting the previous object when no analysis
existed), so the next `analysis` silently built a default. `OPS_Analysis`
returned 0 for an unknown analysis type.

classic Tcl: `specifyCTest`, `specifyNumberer` and `specifyAnalysis` were
already correct (local temporaries, TCL_ERROR before touching a global).
`specifySOE` returned TCL_OK for an unknown type once ANY system had been set
(its final check reads the global, which still held the previous SOE), and
`specifyConstraintHandler`'s Auto / LadrunoProjection / LadrunoContact
branches assigned the factory result to the global BEFORE the null check, so
a refused handler returned TCL_ERROR with the previous handler discarded.

Observables for "the previous object was kept" (1-DOF linear static model):
  * system / numberer / constraints: re-issuing `analysis Static` prints
    `no LinearSOE specified` / `no Numberer specified` /
    `no ConstraintHandler yet specified` iff the global was nulled.
  * test: the previous test is `FixedNumIter 3` -> `testIter` == 3. The silent
    default (`NormUnbalance 1e-6 25`) converges this linear model in 1.
  * analysis: the previous Static analysis still runs (`getTime` == 0.1).
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DT = 0.1
DEFAULT_WARNING = {
    "system": "no LinearSOE specified",
    "numberer": "no Numberer specified",
    "constraints": "no ConstraintHandler yet specified",
}


def _model():
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
    ops.integrator("LoadControl", DT)


# (command, refused args). Every factory listed returns null for these args;
# `unknown-type` rows are premises (the dispatcher already returned -1).
REFUSED = [
    pytest.param("system", ("PFEM", "-noSuchOption"), id="system-PFEM-bad-option"),
    pytest.param("system", ("Mumps", "-ICNTL14"), id="system-Mumps-missing-value"),
    pytest.param("system", ("NoSuchSystem",), id="system-unknown-type"),
    pytest.param("numberer", ("ParallelRCM",), id="numberer-ParallelRCM-serial"),
    pytest.param("numberer", ("NoSuchNumberer",), id="numberer-unknown-type"),
    pytest.param("constraints", ("Penalty",), id="constraints-Penalty-no-args"),
    pytest.param("constraints", ("Auto", "-autoPenalty"), id="constraints-Auto-missing-value"),
    pytest.param("constraints", ("NoSuchHandler",), id="constraints-unknown-type"),
    pytest.param("test", ("NormUnbalance",), id="test-NormUnbalance-no-args"),
    pytest.param("test", ("EnergyIncr", 1e-8), id="test-EnergyIncr-missing-iter"),
    pytest.param("test", ("NoSuchTest", 1e-8, 10), id="test-unknown-type"),
]


@pytest.mark.parametrize("cmd,refused", REFUSED)
@pytest.mark.parametrize("analysis_first", [False, True], ids=["no-analysis", "with-analysis"])
def test_refused_factory_raises_and_keeps_previous(cmd, refused, analysis_first, capfd):
    _model()
    if analysis_first:
        ops.analysis("Static")
    capfd.readouterr()
    with pytest.raises(Exception):
        getattr(ops, cmd)(*refused)
    err = capfd.readouterr().err
    if not refused[0].startswith("NoSuch"):
        assert "left unchanged" in err, err
    # Re-issuing `analysis` rebuilds it from the globals: a nulled global would
    # be replaced by a default here.
    ops.analysis("Static")
    err = capfd.readouterr().err
    for c, warning in DEFAULT_WARNING.items():
        assert warning not in err, f"previous {c} was discarded:\n{err}"
    assert ops.analyze(1) == 0
    assert ops.getTime() == pytest.approx(DT, rel=1e-12)
    assert ops.testIter() == 3, "previous test (FixedNumIter 3) was discarded"


def test_unknown_analysis_type_raises_and_keeps_previous(capfd):
    _model()
    ops.analysis("Static")
    with pytest.raises(Exception):
        ops.analysis("NoSuchAnalysis")
    assert ops.analyze(1) == 0
    assert ops.getTime() == pytest.approx(DT, rel=1e-12)


@pytest.mark.parametrize("analysis_first", [False, True], ids=["no-analysis", "with-analysis"])
def test_accepted_objects_still_replace(analysis_first):
    """Premise: the guards must not block a good object."""
    _model()
    if analysis_first:
        ops.analysis("Static")
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.test("FixedNumIter", 2)
    if not analysis_first:
        ops.constraints("Plain")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.testIter() == 2


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
node 1 0.0; node 2 1.0; fix 1 1
uniaxialMaterial Elastic 1 1.0; element Truss 1 1 2 1.0 1
timeSeries Linear 1; pattern Plain 1 1 { load 2 1.0 }
constraints Transformation; numberer Plain; system FullGeneral
test FixedNumIter 3; algorithm Newton; integrator LoadControl 0.1
"""

# name -> (refused command, fixed by WP-172?)  Unfixed rows are premises that
# were already correct pre-fix and must stay so.
TCL_CASES = {
    "system-unknown-after-system": ("system NoSuchSystem", True),
    "constraints-Auto-missing-value": ("constraints Auto -autoPenalty", True),
    "numberer-unknown": ("numberer NoSuchNumberer", False),
    "test-NormUnbalance-no-args": ("test NormUnbalance", False),
    "test-unknown": ("test NoSuchTest 1e-8 10", False),
    "analysis-unknown": ("analysis NoSuchAnalysis", False),
}


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
@pytest.mark.parametrize("case", list(TCL_CASES))
def test_classic_tcl_refused_factory(case, tmp_path):
    refused, fixed = TCL_CASES[case]
    deck = tmp_path / f"{case}.tcl"
    deck.write_text(TCL_MODEL + f"""
puts "RC [catch {{{refused}}} m]"
analysis Static
analyze 1
puts "TIME [getTime]"
puts "ITER [testIter]"
puts "DECK-DONE"
""")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see tests/test_wp103_getstringfromall_tcl.py (WinError 6).
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=120, env=env, cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, f"exit {proc.returncode} (139/-11/0xC0000005 = segfault):\n{out}"
    assert "DECK-DONE" in out, out
    assert "RC 1" in out, f"refused `{refused}` must be a Tcl error:\n{out}"
    if fixed:
        assert "left unchanged" in out, out
    for c, warning in DEFAULT_WARNING.items():
        assert warning not in out, f"previous {c} was discarded:\n{out}"
    t = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("TIME ")]
    assert t == pytest.approx([DT], rel=1e-12), out
    assert "ITER 3" in out, "previous test (FixedNumIter 3) was discarded:\n" + out
