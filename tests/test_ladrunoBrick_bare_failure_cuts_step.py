"""LadrunoBrick honours a bare -1 from ``setTrialStrain`` (Ladruno C3b).

Before C3b every ``setTrialStrain`` call site in ``LadrunoBrick.cpp`` compared the return code
with ``LADRUNO_MATERIAL_REFUSED`` only, so a material that failed with the plain OpenSees code
-1 (LadrunoRCConcrete's loud crack-band failure, a condensation miss, ...) reached the analysis
as a SUCCESSFUL state determination -- unlike ``LadrunoQuad`` and ``TenNodeTetrahedron``, which
propagate it. The element now cuts the step on the sentinel OR -1; ASDConcrete3D's advisory
codes (-10 IMPL-EX error control, -1000 eigen) stay accepted (ADR-86b: a blanket ``< 0`` broke
test_ladrunoBrick_asdconcrete_bend.py).

Forced failure: ``nDMaterial StagedStrain ... -maxStrain e`` returns -1 (without touching the
inner) once |eps - eps0|_inf > e.
"""
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
def test_bare_minus_one_cuts_the_step_and_a_smaller_step_recovers(form):
    _build(form)
    ops.integrator("LoadControl", 5.0e-4)
    assert ops.analyze(1) == 0                         # 5e-4 < guard
    s_ok = ops.eleResponse(1, "stresses")[:6]
    ops.integrator("LoadControl", 1.0e-3)
    assert ops.analyze(1) != 0                          # 1.5e-3 > guard: -1 -> step cut
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
