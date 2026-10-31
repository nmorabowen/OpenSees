"""ASDConcrete1D's per-tag hardening-law registry must not survive `wipe` (apeGmsh #2044).

The uniaxial sibling of test_asdconcrete3d_wipe_storage.py: ASDConcrete1D keeps its
own store-if-absent `HardeningLawStorage`, and `HardeningLaw::regularize()` copies
the stored law back over the instance (`deRegularize`). Before the fix a later
model's same-tag material ran on the FIRST model's laws, elastic modulus included:
E=3e6 then E=30e6 under tag 1 gave a singular truss (u = 0, stress 0).

Oracle: the closed form on model 2, one truss (L = 2, A = 1) under an axial
compression P = 7200 -- stress = -P/A, u = -P L/(E A), reaction = P. Model 1 (same
tag, E = 3e6) is defined and wiped first. `-autoRegularization` with lch_ref != lch
is required: it is what makes regularize() reach the registry.

The same deck runs through the classic Tcl engine (tests/tcl/
asdconcrete1d_wipe_storage.tcl), whose `wipe` is a separate code path.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "asdconcrete1d_wipe_storage.tcl"

E, P, L = 30.0e6, 7200.0, 2.0
TAG = 1          # deliberately the common tag: the bug is a same-tag collision


def test_same_tag_after_wipe_gets_its_own_law():
    ops.wipe()
    ops.model("BasicBuilder", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("ASDConcrete1D", TAG, 3.0e6, "-fc", 3.0e3, "-ft", 3.0e2)
    ops.wipe()

    ops.model("BasicBuilder", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("ASDConcrete1D", TAG, E, "-fc", 30.0e3, "-ft", 3.0e3,
                         "-autoRegularization", 0.03112529146138255)
    ops.element("truss", 1, 1, 2, 1.0, TAG)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -P, 0.0)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    try:
        assert ops.analyze(1) == 0, "Newton did not converge on the elastic step"
        stress = ops.eleResponse(1, "material", "stress")[0]
        assert stress == pytest.approx(-P, rel=1e-5), f"stress = {stress}"
        assert ops.nodeDisp(2, 1) == pytest.approx(-P * L / E, rel=1e-5)
        ops.reactions()
        assert ops.nodeReaction(1, 1) == pytest.approx(P, rel=1e-5)
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
    assert "ASDC1D_WIPE_2044:" in out, f"deck did not reach its verdict:\n{out}"
    assert out.count("PASS ") >= 4, f"too few checks ran:\n{out}"
    assert "ASDC1D_WIPE_2044: PASS" in out, out
