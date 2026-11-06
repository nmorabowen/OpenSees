"""ASDConcrete3D's per-tag hardening-law registry must not survive `wipe` (apeGmsh #2044).

`HardeningLawStorage` (ASDConcrete3DMaterial.cpp) keeps the ORIGINAL law per
material tag, stores only into an empty slot, and `HardeningLaw::regularize()`
copies it back over the instance (`deRegularize`). Before the fix nothing cleared
it, so the first ASDConcrete3D ever defined with tag N owned tag N for the life
of the process: a later model's same-tag material silently ran on the earlier
backbones, elastic modulus included. apeGmsh #2044 saw it as Newton diverging in
the elastic first step of a pytest process (sigma_zz = -720 for -7200).

Oracle: the uniaxial closed form on model 2 -- sigma_zz = -P/A at every Gauss
point, uz = -P L/(E A), ux = nu P/(E A), base reactions = P. Model 1 (same tag,
E = 3e6) is defined and wiped first. `-autoRegularization` with lch_ref != lch
is required: it is what makes regularize() reach the registry.

The same deck runs through the classic Tcl engine (tests/tcl/
asdconcrete3d_wipe_storage.tcl), whose `wipe` is a separate code path.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "asdconcrete3d_wipe_storage.tcl"

E, NU, P, L = 30.0e6, 0.2, 7200.0, 2.0
TAG = 1          # deliberately the common tag: the bug is a same-tag collision
TABLES = [
    "-Te", 0.0, 9e-05, 0.00015, 0.00507, 0.0250501, 0.250501,
    "-Ts", 0.0, 2700.0, 3000.0, 600.0, 3.0, 3.0,
    "-Td", 0.0, 0.0, 0.0, 0.960552268244576, 0.9999800399998403, 0.9999995660869531,
    "-Ce", 0.0, 0.0005, 0.001, 0.0015, 0.002, 0.18327272727272728, 0.18377272727272728,
    "-Cs", 0.0, 15000.0, 24852.8137423857, 28923.048454132637, 30000.0, 3000.0, 3000.0,
    "-Cd", 0.0, 0.0, 0.0, 0.0, 0.0, 0.9981618744961699, 0.9981786141748574,
    "-autoRegularization", 0.03112529146138255,
]


def _model_2():
    ops.model("BasicBuilder", "-ndm", 3, "-ndf", 3)
    for t, x, y, z in [(1, 0, 0, 2), (2, 0, 0, 0), (3, 0, 1, 2), (4, 0, 1, 0),
                       (5, 1, 0, 2), (6, 1, 0, 0), (7, 1, 1, 2), (8, 1, 1, 0)]:
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ASDConcrete3D", TAG, E, NU, *TABLES)
    ops.element("LadrunoBrick", 2, 3, 1, 2, 4, 7, 5, 6, 8, TAG)
    for n, f in [(2, (1, 1, 1)), (4, (1, 0, 1)), (6, (0, 1, 1)), (8, (0, 0, 1)),
                 (1, (1, 1, 0)), (3, (1, 0, 0)), (5, (0, 1, 0))]:
        ops.fix(n, *f)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (1, 3, 5, 7):
        ops.load(n, 0.0, 0.0, -P / 4.0)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def test_same_tag_after_wipe_gets_its_own_law():
    ops.wipe()
    ops.model("BasicBuilder", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ASDConcrete3D", TAG, 3.0e6, 0.2, "-fc", 3.0e3, "-ft", 3.0e2)
    ops.wipe()

    _model_2()
    try:
        assert ops.analyze(1) == 0, "Newton did not converge on the elastic step"
        for gp in range(1, 9):
            szz = ops.eleResponse(2, "material", gp, "stress")[2]
            assert szz == pytest.approx(-P, rel=1e-5), f"GP {gp}: sigma_zz = {szz}"
        for n in (1, 3, 5, 7):
            assert ops.nodeDisp(n, 3) == pytest.approx(-P * L / E, rel=1e-5)
        for n in (5, 7):
            assert ops.nodeDisp(n, 1) == pytest.approx(NU * P / E, rel=1e-5)
        ops.reactions()
        rz = sum(ops.nodeReaction(n, 3) for n in (2, 4, 6, 8))
        assert rz == pytest.approx(P, rel=1e-5)
    finally:
        ops.wipe()


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    name = "OpenSees.exe" if os.name == "nt" else "OpenSees"
    for cand in (REPO / "dist" / "bin" / name,           # build.bat output
                 REPO / "build" / "Release" / name):     # Zone-A CI build tree
        if cand.exists():
            return str(cand)
    return None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_wipe_clears_the_registry():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see test_ladruno_solver_queries_tcl.py (WinError 6 under capture)
    proc = subprocess.run([TCL_EXE, str(DECK)], capture_output=True, text=True,
                          timeout=300, env=env, cwd=str(REPO), stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    # the marker, not the exit code: OpenSees exits 0 after a Tcl error
    assert "ASDC3D_WIPE_2044:" in out, f"deck did not reach its verdict:\n{out}"
    assert out.count("PASS ") >= 16, f"too few checks ran:\n{out}"
    assert "ASDC3D_WIPE_2044: PASS" in out, out
