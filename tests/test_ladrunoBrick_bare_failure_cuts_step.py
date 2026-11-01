"""A failing material cuts the LadrunoBrick step through the ADR-94 sentinel (Ladruno C3b).

Contract (ADR-86b, pinned fork-wide by ADR-94 R5): every ``setTrialStrain`` call site in
``LadrunoBrick.cpp`` compares the return code with ``LADRUNO_MATERIAL_REFUSED`` and with nothing
else -- never a blanket ``< 0`` / ``!= 0`` (ASDConcrete3D's advisory codes -10 IMPL-EX error
control and -1000 eigen must keep being accepted; a blanket ``< 0`` broke
test_ladrunoBrick_asdconcrete_bend.py), and not a bare -1 either: the host is sentinel-only, so a
material that fails CONVERTS its failure into the sentinel.  Producers converted here:
``StagedStrain -maxStrain`` and ``LadrunoRCConcrete``'s crack-band failures (regularisation miss,
no-converge, sigma-33 condensation miss).

Forced failure: ``nDMaterial StagedStrain ... -maxStrain e`` returns the sentinel (without
touching the inner) once |eps - eps0|_inf > e.  Chain pinned: StagedStrain -maxStrain -> sentinel
-> LadrunoBrick cuts the step on every formulation -> a smaller step from the same committed
state reproduces the unguarded answer.
"""
import os
import re

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

E, NU = 30000.0, 0.2
CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
FORMS = [["-formulation", "std"], ["-formulation", "bbar"], ["-formulation", "ssp"],
         ["-formulation", "uri"], ["-formulation", "uri", "-hourglass", "stiffness"],
         ["-formulation", "eas"]]


def _build(form, guard=1.0e-3):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in CUBE.items():
        ops.node(t, *map(float, c))
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    if guard is None:
        mat = 1
    else:
        ops.nDMaterial("StagedStrain", 2, 1, "-noInit", "-maxStrain", guard)   # pass-through + guard
        mat = 2
    ops.element("LadrunoBrick", 1, *CUBE.keys(), 1 if guard is None else mat, *form)
    # bottom: uz = 0 (+ minimal in-plane restraint); top: prescribed uz; laterals free
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.sp(t, 3, 1.0)                 # lambda = the axial strain
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 10, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


@pytest.mark.parametrize("form", FORMS)
def test_sentinel_cuts_the_step_and_a_smaller_step_recovers(form):
    _build(form)
    ops.integrator("LoadControl", 5.0e-4)
    assert ops.analyze(1) == 0                         # 5e-4 < guard
    s_ok = ops.eleResponse(1, "stresses")[:6]
    ops.integrator("LoadControl", 1.0e-3)
    assert ops.analyze(1) != 0                          # 1.5e-3 > guard: sentinel -> step cut
    # the committed state is untouched: a smaller step from the SAME committed point converges
    # to exactly the unguarded material's answer
    ops.integrator("LoadControl", 4.0e-4)
    assert ops.analyze(1) == 0                          # 9e-4 < guard
    s9 = ops.eleResponse(1, "stresses")[:6]
    _build(form, guard=None)
    ops.integrator("LoadControl", 5.0e-4)
    assert ops.analyze(1) == 0
    ops.integrator("LoadControl", 4.0e-4)
    assert ops.analyze(1) == 0
    ref = ops.eleResponse(1, "stresses")[:6]
    assert s9 == pytest.approx(ref, rel=1e-10, abs=1e-12)
    assert s_ok != s9


def test_unguarded_material_is_unaffected():
    """No -maxStrain: StagedStrain is a pass-through and the element runs past the old guard."""
    _build(["-formulation", "bbar"], guard=None)
    ops.integrator("LoadControl", 1.0e-3)
    for _ in range(3):
        assert ops.analyze(1) == 0


def test_maxstrain_parser_refuses_non_positive():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    with pytest.raises(Exception):
        ops.nDMaterial("StagedStrain", 2, 1, "-maxStrain", 0.0)


# --- structural pins: the sentinel-only contract (ADR-94 R5) and its producers -------------------
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _src(*parts):
    with open(os.path.join(_REPO, *parts), encoding="utf-8", errors="replace") as fh:
        return fh.read()


def test_brick_is_sentinel_only_and_asd_advisory_codes_stay_accepted():
    """No blanket failure test on a setTrialStrain return (the ASDConcrete3D advisory -10/-1000 are
    accepted), and no bare -1 test either: the brick cuts on the sentinel only."""
    src = _src("SRC", "element", "ladrunoBrick", "LadrunoBrick.cpp")
    assert not re.findall(r"setTrialStrain\([^)]*\)\s*(?:\n\s*)?(?:<\s*0|!=\s*0|==\s*-1)(?!\w)", src)
    assert "ladrunoBrickMustCut" not in src and "lastTrialRc" not in src
    assert len(re.findall(r"setTrialStrain\([^)]*\)\s*(?:\n\s*)?==\s*LADRUNO_MATERIAL_REFUSED", src)) == 6


def test_producers_return_the_sentinel_not_a_bare_minus_one():
    staged = _src("SRC", "material", "nD", "StagedStrainNDMaterial.cpp")
    guard = re.findall(r"if \(m > maxStrain\) return ([^;]+);", staged)
    assert guard == ["LADRUNO_MATERIAL_REFUSED", "LADRUNO_MATERIAL_REFUSED"], guard
    rc = _src("SRC", "material", "nD", "LadrunoRCConcrete.cpp")
    i = rc.index("int LadrunoRCConcrete::setTrialStrain(const Vector& e)")
    body = rc[i:rc.index("int LadrunoRCConcrete::setTrialStrain(const Vector& v, const Vector&)", i)]
    assert "return -1" not in body and "? -1 :" not in body
    assert body.count("LADRUNO_MATERIAL_REFUSED") >= 3


def test_wrapping_hosts_forward_only_the_sentinel():
    assert "ladrunoMaterialMustCut" not in _src("SRC", "material", "LadrunoMaterialStatus.h")
    shell = _src("SRC", "element", "ladrunoSolidShell", "LadrunoSolidShell.cpp")
    assert "rcm == LADRUNO_MATERIAL_REFUSED" in shell
    logs = _src("SRC", "material", "nD", "LogStrainNDMaterial.cpp")
    assert "irc == LADRUNO_MATERIAL_REFUSED" in logs
