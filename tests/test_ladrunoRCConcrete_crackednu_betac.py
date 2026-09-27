"""LadrunoRCConcrete (classTag 33015) — Ladruno C2: ``-crackedNu``, ``-betaC`` and the vc
tension-stiffening default. Zone-A, the material as integrated in OpenSees (stdBrick,
homogeneous strain).

* ``-crackedNu nu_c`` (default off): once the in-plane principal tensile strain reaches the
  cracking strain eps_cr (= ft/E unless -crackStrain), the elastic operator of the
  effective-stress predictor switches to nu_c (irreversible latch, response ``nuCracked``).
  The PV panels showed nu = 0.2 retained after cracking costs 8-10 % of the shear strength:
  the transverse tensile strain unloads the strut effective stress (sig2_eff ~ eps2 + nu eps1).
* ``-betaC C`` (default 170): beta = 1/(0.8 + C eps1). Vecchio-Collins 1986 is 0.34/|eps'c|.
* ``-tensStiff vc`` now defaults to c = 200 (Vecchio-Collins 1986); covered in
  test_ladrunoRCConcrete_tensstiff.py::test_vc_cm_default_equivalence.

With the new flags off the material is byte-identical to the pre-C2 build (the rest of the
test_ladrunoRCConcrete_*.py suite, unchanged, is that gate).
"""
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]
EPS_CR = 3.0 / E                         # = ft/E = 1e-4 (the kernel's default crack strain)
LAM = E * NU / ((1 + NU) * (1 - 2 * NU))  # Lame lambda at nu = 0.2


def _rc(tag, *flags):
    ops.nDMaterial("LadrunoRCConcrete", tag, E, NU,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC, *flags)


def _two_stage(flags, eyy, n1, exx, n2, probe=None):
    """Stage 1: ramp eyy (n1 equal steps), exx = ezz = 0. Stage 2: hold eyy, ramp exx.
    Every nodal displacement prescribed (homogeneous strain). Returns centroid stress."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    _rc(1, *flags)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 2, eyy * y)
        ops.sp(t, 3, 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n1)
    ops.analysis("Static")
    lat = []
    for _ in range(n1):
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
        if probe:
            lat.append(ops.eleResponse(1, "material", 1, probe)[0])
    ops.loadConst("-time", 0.0)
    ops.remove("loadPattern", 2)
    ops.pattern("Plain", 3, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)
    ops.integrator("LoadControl", 1.0 / n2)
    for _ in range(n2):
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
        if probe:
            lat.append(ops.eleResponse(1, "material", 1, probe)[0])
    return list(ops.eleResponse(1, "stresses"))[0:6], lat


@pytest.mark.t1
def test_cracked_nu_zero_decouples_the_strut_closed_form():
    """Stage 1: eyy = 2e-3 in 3e-5 steps -> the latch fires on step 4 (e1 = 1.2e-4 >= eps_cr).
    Stage 2: hold eyy, compress exx (elastic range: dc = 0, beta off). Once sig_xx is the
    compressive principal it is the unreduced compressive cone, so with nu_c = 0 the strut is
    decoupled from the transverse tension: d sig_xx / d exx = E exactly (path independent),
    and no eyy term enters. With nu kept at 0.2 the strut carries the Poisson term
    sig_xx_eff ~ LAM*eyy + (LAM+2mu)*exx: at exx = -3e-4 it is still in TENSION (the transverse
    tension erased the strut -- the PV-panel mechanism)."""
    n1 = int(round(2.0e-3 / 3.0e-5))
    s2, lat = _two_stage(["-crackedNu", 0.0], 3.0e-5 * n1, n1, -2.0e-4, 20, probe="nuCracked")
    s3, _ = _two_stage(["-crackedNu", 0.0], 3.0e-5 * n1, n1, -3.0e-4, 30)
    assert s2[0] < 0.0 and s3[0] < 0.0
    assert s3[0] - s2[0] == pytest.approx(E * -1.0e-4, rel=1e-9), (s2[0], s3[0])
    assert lat[:3] == [0.0, 0.0, 0.0] and all(v == 1.0 for v in lat[3:])
    # nu kept at 0.2: the strut stiffness is the coupled lambda + 2 mu (and the eyy Poisson
    # term sits in its effective stress) -- the 11 % stiffer, tension-contaminated strut.
    n2, _ = _two_stage([], 3.0e-5 * n1, n1, -2.0e-4, 20)
    n3, _ = _two_stage([], 3.0e-5 * n1, n1, -3.0e-4, 30)
    mu = E / (2.0 * (1.0 + NU))
    assert n3[0] - n2[0] == pytest.approx((LAM + 2.0 * mu) * -1.0e-4, rel=1e-9), (n2[0], n3[0])
    assert s3[0] < n3[0] - 0.5, (s3[0], n3[0])       # nu_c = 0 strut is the more compressive


def test_cracked_nu_inert_before_cracking_bit_identical():
    """A path that never reaches eps_cr (confined elastic compression, e1 = 0) is bit-identical
    with and without the flag; so is -betaC 170 (the old hard-wired value) with beta on."""
    a, _ = _two_stage([], 0.0, 5, -4.0e-4, 20)
    b, _ = _two_stage(["-crackedNu", 0.0], 0.0, 5, -4.0e-4, 20)
    assert a == b
    c, _ = _two_stage(["-beta", "-lublinerReduced"], 2.0e-3, 40, -1.5e-3, 40)
    d, _ = _two_stage(["-beta", "-lublinerReduced", "-betaC", 170.0], 2.0e-3, 40, -1.5e-3, 40)
    assert c == d


def test_cracked_nu_latch_is_irreversible():
    """Crack in y, then unload y to e1 = 1e-4 (below the 5e-4 peak): nuCracked stays 1."""
    ops.wipe()
    _, lat = _two_stage(["-crackedNu", 0.0], 5.0e-4, 10, 0.0, 2, probe="nuCracked")
    assert lat[-1] == 1.0
    # unload: a second call cannot unload inside the helper, so check via a fresh ramp down
    ops.loadConst("-time", 0.0)
    ops.remove("loadPattern", 3)
    ops.pattern("Plain", 4, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 2, -4.0e-4 * y)       # stage-1 pattern is constant (+5e-4) -> total 1e-4
    ops.test("NormDispIncr", 1.0e-9, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.05)
    # This unloading path sits right at a marginal corner of the Newton basin (confirmed locally:
    # KrylovNewton and NewtonLineSearch both converge here at 20 steps/1e-9, but plain Newton, more
    # steps, or a looser tolerance all fail even on Windows/MKL -- the step is genuinely tight, not
    # merely mistuned). The Zone-A Ubuntu runner failed this exact step with -3 (AcceleratedNewton,
    # KrylovNewton's base class) while it passes here, i.e. a platform BLAS/LAPACK-level difference in
    # how close the two algorithms land to that basin boundary. Retry each step with NewtonLineSearch
    # before failing so the assertion below (the physical claim: the latch is irreversible) is not at
    # the mercy of that platform-dependent margin.
    for _ in range(20):
        if ops.analyze(1) != 0:
            ops.algorithm("NewtonLineSearch")
            assert ops.analyze(1) == 0, "latch unload step failed to converge under both algorithms"
            ops.algorithm("KrylovNewton")
    ops.eleResponse(1, "forces")
    assert ops.eleResponse(1, "material", 1, "nuCracked")[0] == 1.0


def _homog_biax(flags, exx, eyy, nsteps=100):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    _rc(1, *flags)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 2, eyy * y)
        ops.sp(t, 3, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
    return list(ops.eleResponse(1, "stresses"))[0:6]


@pytest.mark.t1
@pytest.mark.parametrize("C", [170.0, 189.0, 300.0])
def test_betaC_scales_the_strut_exactly(C):
    """test_compression_softening_biaxial with a user C: beta-ON / beta-OFF compressive
    principal ratio == 1/(0.8 + C*e1) (lublinerReduced off, so beta touches only the cone)."""
    import numpy as np
    e1 = 3.0e-3
    off = _homog_biax([], -2.0e-3, e1)
    on = _homog_biax(["-beta", "-betaC", C], -2.0e-3, e1)
    mn = lambda s: float(np.min(np.linalg.eigvalsh(np.array(
        [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]]))))
    assert mn(on) / mn(off) == pytest.approx(1.0 / (0.8 + C * e1), rel=1e-3)


def test_serialization_carries_the_c2_params():
    """Schema v6: betaC, crackedNu and the committed nuCracked latch survive a database
    round trip (a dropped slot reverts to the defaults and changes the probed stress).
    The path ends on a TENSILE-UNLOADING stage-2 (exx > 0 would load; exx < 0 grows the
    Lubliner tensile measure through J2 and keeps damage active): re-evaluating a committed
    state whose last step was actively damaging is not idempotent in this kernel (the
    committed effective stress is the post-reduction one, a property inherited from
    ASDConcrete3D), so the probe must sit on an elastic step."""
    def build():
        _two_stage(["-beta", "-betaC", 250.0, "-crackedNu", 0.0], 1.0e-3, 10, 0.0, 1)
        ops.loadConst("-time", 0.0)
        ops.remove("loadPattern", 3)
        ops.pattern("Plain", 5, 1)
        for t, (x, y, z) in _CUBE.items():
            ops.sp(t, 2, -2.0e-4 * y)     # unload eyy 1e-3 -> 8e-4 (elastic unloading)
        ops.integrator("LoadControl", 0.5)
        for _ in range(2):
            assert ops.analyze(1) == 0
        ops.eleResponse(1, "forces")
    database_roundtrip(build, probe_nodes=[2], ndf=3, dbname="rc_c2",
                       probe_fn=lambda: list(ops.eleResponse(1, "stresses"))[0:6]
                       + list(ops.eleResponse(1, "material", 1, "nuCracked")))


@pytest.mark.parametrize("flags", [["-betaC", 0.0], ["-betaC", -5.0], ["-crackedNu", 0.5]])
def test_parser_refuses(flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    with pytest.raises(Exception):
        _rc(1, *flags)
