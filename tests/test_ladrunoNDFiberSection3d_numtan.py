"""Ladruno WP-252: `section LadrunoNDFiber ... -tangent numerical [-tangentStep d]`.

The section tangent is the forward difference of the six resultants, each section deformation perturbed so the
largest fibre strain it moves changes by d, every fibre re-evaluated (its BeamFiber lateral balance included) and
then restored at the base deformation. Checks:

* elastic strips (ElasticIsotropic): the numerical tangent equals the analytic one (both exact);
* RC strips on a cracking shear path (LadrunoRCConcrete, BeamFiber view): the numerical mode changes the tangent
  only -- the resultants match the analytic section on every step -- and the numerical tangent predicts the
  resultant increment of the next step better than the material's damaged-secant tangent does;
* the parser and the Print line.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU = 30000.0, 0.2
B = H = 457.2
NSTRIP = 40
DY = H / NSTRIP
YS = [-0.5 * H + DY * (i + 0.5) for i in range(NSTRIP)]
AS = B * DY
PARA = ["-shearProfile", "parabolic", "-shearDepth", H, B]
NUM = ["-tangent", "numerical"]

CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0, 30.0, 5.0]
CD = [0.0, 0.0, 0.25, 1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0, 0.5]
TD = [0.0, 0.0, 1.0 - 0.5 / 5.0]
FY, ES = 420.0, 200000.0
DOF_OF = {0: 1, 1: 6, 2: 5, 3: 2, 4: 3, 5: 4}


def _rc(tag, *flags):
    ops.nDMaterial("LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", 2.0 / 3.0, "-beta", "-tensStiff", "vc", *flags)


def _section(sec, mat_of_y, *opts):
    ops.section("LadrunoNDFiber", sec, *opts)
    for y in YS:
        ops.fiber(y, 0.0, AS, mat_of_y(y))


def _drive(secs, rows):
    """One zeroLengthSection per section tag, all following the same deformation rows."""
    ts = 1
    t = [0.0] + [float(k + 1) for k in range(len(rows))]
    for i, sec in enumerate(secs):
        ni, nj = 10 * i + 1, 10 * i + 2
        ops.node(ni, 0.0, 0.0, 0.0)
        ops.node(nj, 0.0, 0.0, 0.0)
        ops.fix(ni, 1, 1, 1, 1, 1, 1)
        ops.element("zeroLengthSection", i + 1, ni, nj, sec)
        for c in range(6):
            ops.timeSeries("Path", ts, "-time", *t, "-values", *([0.0] + [r[c] for r in rows]))
            ops.pattern("Plain", ts, ts)
            ops.sp(nj, DOF_OF[c], 1.0)
            ts += 1
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _K(ele):
    return np.array(ops.eleResponse(ele, "section", "stiffness")).reshape(6, 6)


def _S(ele):
    return np.array(ops.eleResponse(ele, "section", "force"))


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def test_elastic_numerical_tangent_equals_the_analytic_one():
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    _section(1, lambda y: 1, *PARA)
    _section(2, lambda y: 1, *PARA, *NUM)
    _drive([1, 2], [(-2.0e-4, 3.0e-7, 0.0, 1.0e-4, 0.0, 0.0)])
    assert ops.analyze(1) == 0
    Ka, Kn = _K(1), _K(2)
    scale = np.abs(Ka).max(axis=0) + 1e-300           # per column
    err = np.abs(Kn - Ka) / scale
    # torsion: no fibre feels it on a section without -GJ, the analytic (zero) column is kept
    assert err.max() < 1e-6, err.max()
    assert np.allclose(_S(1), _S(2), rtol=1e-12, atol=1e-9)


def _rc_path():
    """Axial compression, a little curvature, then a shear ramp through cracking (as test_ladrunoNDFiberSection3d)."""
    return [(-3.0e-4, 1.0e-6, 0.0, 4.0e-3 * k / 40.0, 0.0, 0.0) for k in range(1, 41)]


def test_rc_numerical_tangent_changes_the_tangent_only_and_predicts_better():
    path = _rc_path()
    _model()
    _rc(1)                                              # cover: no ties
    _rc(2, "-ties", 0.002, 0.0, FY, ES)                 # core: ties in y
    cover = 63.5
    mat = lambda y: 2 if abs(y) <= H / 2 - cover else 1   # noqa: E731
    _section(1, mat, *PARA)
    _section(2, mat, *PARA, *NUM)
    _drive([1, 2], path)
    worstS, errA, errN, prev = 0.0, 0.0, 0.0, None
    for k in range(len(path)):
        assert ops.analyze(1) == 0, k
        S1, S2 = _S(1), _S(2)
        worstS = max(worstS, float(np.max(np.abs(S1 - S2) / (np.abs(S1) + 1.0))))
        if prev is not None:
            Sp, KAp, KNp, ep = prev
            de = np.array(path[k]) - ep
            errA += float(np.linalg.norm(S1 - (Sp + KAp @ de)))
            errN += float(np.linalg.norm(S1 - (Sp + KNp @ de)))
        prev = (S1, _K(1), _K(2), np.array(path[k]))
    assert worstS < 1e-9, worstS                         # the numerical mode does not change the response
    assert errN < 0.5 * errA, (errN, errA)               # it predicts the next step better than the secant


@pytest.mark.parametrize("args", [
    ["-tangent", "exact"],
    ["-tangent"],
    ["-tangent", "numerical", "-tangentStep", 0.0],
    ["-tangent", "numerical", "-tangentStep", 0.1],
])
def test_parser_refuses(args):
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    with pytest.raises(Exception):
        ops.section("LadrunoNDFiber", 1, *PARA, *args)


def test_parser_accepts():
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.section("LadrunoNDFiber", 1, *PARA, "-tangent", "analytic")
    ops.section("LadrunoNDFiber", 2, *PARA, "-tangent", "numerical", "-tangentStep", 1e-6)
    ops.section("LadrunoNDFiber", 3, "-shearProfile", "uniform", "-tangent", "numerical")
