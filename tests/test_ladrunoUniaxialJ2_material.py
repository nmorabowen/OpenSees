"""LadrunoUniaxialJ2 (uniaxial combined isotropic Voce+linear + Chaboche AF
kinematic von Mises, MAT_TAG 33000) — Zone-A material battery.

Driver: a unit Truss (L=1, A=1) along x. Axial strain = nodeDisp(2,1); axial
stress = -nodeReaction(1,1) (A=1). A clean single-point uniaxial constitutive
probe driven by displacement control.

Oracles (derived in Ladruno_implementation/13_ladruno_uniaxial_j2_adr.md):
  * V0 elastic:        sigma = E * eps
  * V1 reduce-to-Hardening: linear iso+kin reproduces upstream HardeningMaterial
                       (both closed-form linear 1D plasticity) to ~1e-9
  * V2 monotonic:      linear kinematic C == linear isotropic Hiso=C (monotone)
  * V3 Bauschinger:    kin vs iso DIVERGE on reversal (kin loop narrower)
  * V4 AF saturation:  monotonic single AF -> sig0 + C/gamma
  * V7 THE oracle:     uniaxial reduction of the 3D nDMaterial LadrunoJ2 matches
                       this material's sigma(eps) under uniaxial stress, with the
                       SAME (C_k, gamma_k, iso) -- pins the 3D Chaboche calibration
  * state recording:   plasticStrain / equivalentPlasticStrain / backStress /
                       plasticWork are recordable through the element
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _kg(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


# ==========================================================================
#  1D truss driver (uniaxial stress, displacement-controlled)
# ==========================================================================
def _build_1d(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    mat_fn(1)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)                           # free in x, restrained in y
    ops.element("Truss", 1, 1, 2, 1.0, 1)      # A = 1, L = 1 -> strain=u2, stress=force
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_1d(mat_fn, segments):
    """segments = [(eps_target, nsteps), ...]. Returns [(eps, sig), ...]."""
    _build_1d(mat_fn)
    out = []
    for eps_target, nsteps in segments:
        cur = ops.nodeDisp(2, 1)
        d = (eps_target - cur) / nsteps
        ops.integrator("DisplacementControl", 2, 1, d)
        for _ in range(nsteps):
            assert ops.analyze(1) == 0, f"1D analyze failed heading to eps={eps_target}"
            ops.reactions()
            out.append((ops.nodeDisp(2, 1), -ops.nodeReaction(1, 1)))   # A=1
    return out


# ==========================================================================
#  3D uniaxial-stress driver (for V7): stdBrick unit cube, 1/8-symmetry
#  determinate restraints -> homogeneous uniaxial stress w/ free lateral
#  contraction. Same axial-strain path; read centroid GP sigma_xx.
# ==========================================================================
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def _build_3d(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 1, 0)
    ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.25, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_3d(mat_fn, segments):
    _build_3d(mat_fn)
    out = []
    for eps_target, nsteps in segments:
        cur = ops.nodeDisp(2, 1)
        d = (eps_target - cur) / nsteps
        ops.integrator("DisplacementControl", 2, 1, d)
        for _ in range(nsteps):
            assert ops.analyze(1) == 0, f"3D analyze failed heading to eps={eps_target}"
            ops.eleResponse(1, "forces")
            sig = list(ops.eleResponse(1, "stresses"))[0:6]
            out.append((ops.nodeDisp(2, 1), sig[0]))
    return out


# --------------------------------------------------------------------------
# V0 elastic
# --------------------------------------------------------------------------
@pytest.mark.t0m
def test_elastic_uniaxial():
    E = 1000.0

    def mat(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", 5.0, 0.0, 0.0, 0.0, "-kin", 0)

    eps = 1.0e-4   # below yield strain 5/1000
    e, s = _run_1d(mat, [(eps, 1)])[-1]
    assert s == pytest.approx(E * e, rel=1e-8), f"elastic sigma {s} != E*eps {E*e}"


# --------------------------------------------------------------------------
# V1: reduce to upstream HardeningMaterial (linear iso + linear kinematic).
# Both are closed-form linear 1D plasticity -> must agree to ~1e-9.
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_reduce_to_HardeningMaterial():
    E, sigmaY, Hiso, Hkin = 1000.0, 5.0, 40.0, 60.0
    segs = [(0.03, 60), (-0.03, 120)]   # push then fully reverse (exercises kinematic)

    def mat_ladruno(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sigmaY, 0.0, 0.0, Hiso, "-kin", 1, Hkin, 0.0)

    def mat_hard(tag):
        ops.uniaxialMaterial("Hardening", tag, E, sigmaY, Hiso, Hkin)

    ours = _run_1d(mat_ladruno, segs)
    ref = _run_1d(mat_hard, segs)
    assert len(ours) == len(ref)
    for (eo, so), (er, sr) in zip(ours, ref):
        assert eo == pytest.approx(er, abs=1e-12)
        assert so == pytest.approx(sr, rel=1e-9, abs=1e-9), f"sig {so} vs {sr} at eps={eo}"


# --------------------------------------------------------------------------
# V2: monotonic equivalence — linear kinematic C == linear isotropic Hiso=C
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_linear_kinematic_equals_isotropic_monotonic():
    E, sig0, H = 1000.0, 5.0, 120.0
    segs = [(0.03, 60)]

    def mat_kin(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, H, 0.0)

    def mat_iso(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, H, "-kin", 0)

    kin = _run_1d(mat_kin, segs)
    iso = _run_1d(mat_iso, segs)
    for (ek, sk), (ei, si) in zip(kin, iso):
        assert sk == pytest.approx(si, rel=1e-9, abs=1e-9), (
            f"monotonic kin {sk} != iso {si} at eps={ek}")


# --------------------------------------------------------------------------
# V3: Bauschinger — kinematic vs isotropic diverge on reversal
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_bauschinger_kinematic_vs_isotropic():
    E, sig0, H = 1000.0, 5.0, 120.0
    segs = [(0.03, 60), (-0.03, 120)]

    def mat_kin(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, H, 0.0)

    def mat_iso(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, H, "-kin", 0)

    kin = _run_1d(mat_kin, segs)
    iso = _run_1d(mat_iso, segs)

    assert kin[59][1] == pytest.approx(iso[59][1], rel=1e-7)   # forward peaks coincide
    sk_end, si_end = kin[-1][1], iso[-1][1]
    assert sk_end < 0.0 and si_end < 0.0
    assert abs(sk_end) < abs(si_end) - 1e-6, (
        f"kinematic reverse |sig|={abs(sk_end):.4f} should be < isotropic "
        f"|sig|={abs(si_end):.4f} (Bauschinger)")


# --------------------------------------------------------------------------
# V4: Armstrong-Frederick monotonic saturation -> sigma -> sig0 + C/gamma
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_af_monotonic_saturation():
    E, sig0, C, gamma = 1000.0, 5.0, 400.0, 80.0   # X_sat = C/gamma = 5

    def mat(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, C, gamma)

    s_end = _run_1d(mat, [(0.15, 300)])[-1][1]
    s_sat = sig0 + C / gamma                         # = 10.0
    assert s_end == pytest.approx(s_sat, rel=1e-3), (
        f"AF saturation sigma {s_end} != sig0+C/gamma {s_sat}")


# --------------------------------------------------------------------------
# V7: THE oracle — uniaxial reduction of the 3D nDMaterial LadrunoJ2 must match
# this uniaxial material's sigma(eps) under uniaxial stress, using IDENTICAL
# (C_k, gamma_k, iso). Full Chaboche N=2 over a push-pull-cycle path.
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_oracle_vs_3D_LadrunoJ2():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)
    sig0, Qinf, b, Hiso = 5.0, 3.0, 25.0, 8.0
    # two Chaboche backstresses (distinct recovery rates)
    C1, g1, C2, g2 = 300.0, 60.0, 120.0, 12.0
    segs = [(0.02, 50), (-0.02, 100), (0.025, 100)]

    def mat_1d(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, Qinf, b, Hiso,
                             "-kin", 2, C1, g1, C2, g2)

    def mat_3d(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, Qinf, b, Hiso,
                       "-kin", 2, C1, g1, C2, g2)

    one = _run_1d(mat_1d, segs)
    three = _run_3d(mat_3d, segs)
    assert len(one) == len(three)

    # must actually be in the plastic regime (guard against vacuous elastic match)
    assert max(abs(s) for _, s in one) > sig0 + 1.0, "expected plastic stress"
    for (e1, s1), (e3, s3) in zip(one, three):
        assert e1 == pytest.approx(e3, abs=1e-12)
        # both share the SAME hardening law (LadrunoHardening.h) + the AF reduction;
        # residual diff is the 3D element's lateral-equilibrium iteration only.
        assert s1 == pytest.approx(s3, rel=1e-7, abs=1e-7), (
            f"1D sigma {s1} != 3D-uniaxial sigma {s3} at eps={e1}")


# --------------------------------------------------------------------------
# state recording: plasticStrain / equivalentPlasticStrain / backStress /
# plasticWork must be recordable through the truss element's material response.
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_state_recording():
    E, sig0, C, gamma = 1000.0, 5.0, 400.0, 80.0

    def mat(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, E,
                             "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, C, gamma)

    _run_1d(mat, [(0.15, 300)])

    ebar = list(ops.eleResponse(1, "material", "equivalentPlasticStrain"))
    pstr = list(ops.eleResponse(1, "material", "plasticStrain"))
    back = list(ops.eleResponse(1, "material", "backStress"))
    wp = list(ops.eleResponse(1, "material", "plasticWork"))
    assert ebar and pstr and back and wp, "material responses not wired through element"

    assert ebar[0] > 0.1, f"equivalentPlasticStrain {ebar[0]} not accumulated"
    # single-AF axial backstress saturates to C/gamma (stress units)
    assert back[0] == pytest.approx(C / gamma, rel=0.02), (
        f"backstress {back[0]} != C/gamma {C/gamma}")
    assert wp[0] > 0.0, f"plasticWork {wp[0]} should be positive (dissipation)"
