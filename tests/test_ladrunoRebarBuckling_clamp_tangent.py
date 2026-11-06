"""LadrunoRebarBuckling: the tangent must follow the stress when g_law is clamped.

For a stocky bar (s/d <= 2) fStar exceeds the backbone onset stress, so r > 1 and the
buckling law sits BELOW bare (g_law = sigma_v1 - sigma_bare < 0). The raise is clamped to
g >= 0, so sigma = sigma_bare, and the tangent must be k_bare, not the law's k_v1. A scan
found no clamped state at s/d >= 3, including the recipe's s/d = 6.

setStrain commits, so each finite-difference stencil point gets its own fresh model
driven through the same committed history.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_FY, _E = 438.0, 200000.0
_H = 1.0e-7


def _model(wrapped, lsr):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Steel02", 1, _FY, _E, 0.01, 20.0, 0.925, 0.15)
    if wrapped:
        ops.uniaxialMaterial("LadrunoRebarBuckling", 2, 1, "-lsr", lsr, "-model", "dm",
                             "-alpha", 1.0, "-fy", _FY, "-E", _E, "-restraighten", "lambda")
    ops.testUniaxialMaterial(2 if wrapped else 1)


def _path(peaks, de=2.0e-4):
    out, cur = [], 0.0
    for p in peaks:
        n = max(1, int(abs(p - cur) / de))
        out += [cur + (p - cur) * (k + 1) / n for k in range(n)]
        cur = p
    return out


def _at(wrapped, lsr, committed, eps):
    _model(wrapped, lsr)
    for e in committed:
        ops.setStrain(e)
    ops.setStrain(eps)
    return ops.getStress(), ops.getTangent()


def _state(lsr, committed, eps):
    s_w, k_w = _at(True, lsr, committed, eps)
    s_b, k_b = _at(False, lsr, committed, eps)
    sp, _ = _at(True, lsr, committed, eps + _H)
    sm, _ = _at(True, lsr, committed, eps - _H)
    return dict(s_w=s_w, k_w=k_w, s_b=s_b, k_b=k_b, fd=(sp - sm) / (2.0 * _H))


# lsr, committed peaks, trial offset from the last committed strain
_HISTORIES = {
    "monotonic_lsr2": (2.0, [-0.0036], -2.0e-4),
    "monotonic_lsr1": (1.0, [-0.0036], -2.0e-4),
    "after_tension_lsr2": (2.0, [0.02, 0.0078], -2.0e-4),
    "reversed_lsr2": (2.0, [0.03, -0.02, 0.012, 0.0], -2.0e-4),
}


@pytest.mark.parametrize("name", sorted(_HISTORIES))
def test_tangent_matches_stress_at_the_clamp(name):
    lsr, peaks, de = _HISTORIES[name]
    committed = _path(peaks)
    st = _state(lsr, committed, committed[-1] + de)
    # the history really is at the clamp: stress equals bare past the onset
    assert st["s_w"] == pytest.approx(st["s_b"], abs=1e-9), st
    assert st["fd"] == pytest.approx(st["k_b"], rel=1e-3, abs=1.0), st
    assert st["k_w"] == pytest.approx(st["fd"], rel=1e-3, abs=1.0), st


def test_unclamped_buckling_tangent_is_unchanged():
    # s/d = 12, real softening (g > 0): the law's own tangent still applies.
    committed = _path([-0.02])
    st = _state(12.0, committed, committed[-1] - 2.0e-4)
    assert st["s_w"] > st["s_b"] + 1.0, st
    assert st["k_w"] == pytest.approx(st["fd"], rel=1e-3, abs=1.0), st
    assert abs(st["k_w"] - st["k_b"]) > 1.0, st


def _truss_reduction_history(peaks, de):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Steel02", 1, _FY, _E, 0.01, 20.0, 0.925, 0.15)
    ops.uniaxialMaterial("LadrunoRebarBuckling", 2, 1, "-lsr", 12.0, "-model", "dm",
                         "-alpha", 1.0, "-fy", _FY, "-E", _E, "-restraighten", "lambda")
    ops.element("Truss", 1, 1, 2, 1.0, 2)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    out = []
    cur = 0.0
    for e in _path(peaks, de=de):
        ops.integrator("DisplacementControl", 2, 1, e - cur)
        assert ops.analyze(1) == 0
        cur = e
        out.append(ops.eleResponse(1, "material", "reduction")[0])
    return out


def test_reduction_readout_is_bounded_where_bare_stress_crosses_zero():
    # sigma/sigma_bare divides by a stress that passes through zero while a raise is held,
    # so the raw ratio is a large number with no meaning. The fine step lands close to zero.
    hist = _truss_reduction_history([-0.02, 0.01], 1.0e-6)
    assert max(abs(r) for r in hist) < 250.0
    assert any(abs(r - 1.0) > 0.01 for r in hist)   # the readout still reports buckling
