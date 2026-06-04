"""LadrunoJ2 (combined isotropic Voce+linear + Chaboche AF kinematic von Mises,
classTag 33011) — Zone-A material battery.

Driver: a single ``stdBrick`` unit cube under statically-determinate
(1/8-symmetry) restraints develops a *homogeneous uniaxial stress* state with
free lateral contraction. We drive the x=1 face by displacement control and read
the centroid Gauss-point stress — a clean single-point constitutive probe.

Oracles (all derived in Ladruno_implementation/10_ladruno_j2_plasticity.md):
  * elastic:   sigma_xx = E * eps_xx
  * V1:        N=0 (pure Voce+linear isotropic) reproduces upstream J2Plasticity
  * monotonic equivalence (dSNPO eq 6.222): under *monotonic* loading a linear
    kinematic term C is identical to a linear isotropic modulus Hiso=C
    (this pins the (2/3)C numerator scaling).
  * Bauschinger: the same C as kinematic vs isotropic must DIVERGE on reversal,
    with the kinematic model yielding earlier in reverse (narrower loop).
  * AF saturation: monotonic single-AF stress -> sigma_y0 + C/gamma  (standard
    Chaboche C, gamma).

Plan: Ladruno_implementation/10_ladruno_j2_plasticity.md.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def _kg(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


def _build(mat_fn):
    """Single stdBrick unit cube, 1/8-symmetry determinate restraints, uniform
    x-traction reference load. Caller drives node 2 dof 1 by displacement."""
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
        ops.load(n, 0.25, 0.0, 0.0)        # unit x-traction reference
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run(mat_fn, segments):
    """segments = [(eps_target, nsteps), ...]. Returns list of (eps_xx, sig_xx)."""
    _build(mat_fn)
    out = []
    for eps_target, nsteps in segments:
        cur = ops.nodeDisp(2, 1)
        d = (eps_target - cur) / nsteps
        ops.integrator("DisplacementControl", 2, 1, d)
        for _ in range(nsteps):
            assert ops.analyze(1) == 0, f"analyze failed heading to eps={eps_target}"
            ops.eleResponse(1, "forces")                       # set lazy strain
            sig = list(ops.eleResponse(1, "stresses"))[0:6]    # centroid GP
            out.append((ops.nodeDisp(2, 1), sig[0]))
    return out


# --------------------------------------------------------------------------
# elastic
# --------------------------------------------------------------------------
@pytest.mark.t0m
def test_elastic_uniaxial():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)

    def mat(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", 5.0, 0.0, 0.0, 0.0, "-kin", 0)

    eps = 1.0e-4   # well below yield strain 5/1000
    res = _run(mat, [(eps, 1)])
    e, s = res[-1]
    assert s == pytest.approx(E * e, rel=1e-6), f"elastic sigma {s} != E*eps {E*e}"


# --------------------------------------------------------------------------
# V1: reduce to upstream J2Plasticity (isotropic-only, nonlinear Voce+linear)
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_reduce_to_J2Plasticity():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)
    sig0, Qinf, b, Hiso = 5.0, 4.0, 30.0, 10.0
    segs = [(0.03, 60)]

    def mat_ladruno(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, Qinf, b, Hiso, "-kin", 0)

    def mat_j2(tag):
        # J2Plasticity q(xi) = sigInf + (sig0-sigInf) e^{-delta xi} + H xi
        #            = sig0 + Qinf(1-e^{-b xi}) + Hiso xi   with sigInf=sig0+Qinf
        ops.nDMaterial("J2Plasticity", tag, K, G, sig0, sig0 + Qinf, b, Hiso, 0.0)

    ours = _run(mat_ladruno, segs)
    ref = _run(mat_j2, segs)
    assert len(ours) == len(ref)
    for (eo, so), (er, sr) in zip(ours, ref):
        assert eo == pytest.approx(er, abs=1e-12)
        # J2Plasticity carries an internal gamma*=(1-1e-8) fudge -> ~1e-7 match
        assert so == pytest.approx(sr, rel=2e-6, abs=2e-6), f"sig {so} vs {sr}"


# --------------------------------------------------------------------------
# monotonic equivalence: linear kinematic C == linear isotropic Hiso=C
# (dSNPO eq 6.222) -> pins the (2/3)C numerator scaling
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_linear_kinematic_equals_isotropic_monotonic():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)
    sig0, H = 5.0, 120.0
    segs = [(0.03, 60)]

    def mat_kin(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, H, 0.0)

    def mat_iso(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, H, "-kin", 0)

    kin = _run(mat_kin, segs)
    iso = _run(mat_iso, segs)
    for (ek, sk), (ei, si) in zip(kin, iso):
        assert sk == pytest.approx(si, rel=1e-7, abs=1e-7), (
            f"monotonic kin {sk} != iso {si} at eps={ek}")


# --------------------------------------------------------------------------
# Bauschinger: linear kinematic vs isotropic diverge on reversal; kinematic
# yields earlier in reverse (narrower hysteresis loop)
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_bauschinger_kinematic_vs_isotropic():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)
    sig0, H = 5.0, 120.0
    segs = [(0.03, 60), (-0.03, 120)]   # push then fully reverse

    def mat_kin(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, H, 0.0)

    def mat_iso(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, H, "-kin", 0)

    kin = _run(mat_kin, segs)
    iso = _run(mat_iso, segs)

    # forward peaks must coincide (monotonic equivalence)
    assert kin[59][1] == pytest.approx(iso[59][1], rel=1e-6)

    # at full reversal the curves must differ: kinematic loop is narrower, so the
    # compressive stress magnitude is SMALLER than the isotropic one.
    sk_end = kin[-1][1]
    si_end = iso[-1][1]
    assert sk_end < 0.0 and si_end < 0.0
    assert abs(sk_end) < abs(si_end) - 1e-6, (
        f"kinematic reverse |sig|={abs(sk_end):.4f} should be < isotropic "
        f"|sig|={abs(si_end):.4f} (Bauschinger)")


# --------------------------------------------------------------------------
# Armstrong-Frederick monotonic saturation: sigma -> sig0 + C/gamma
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_af_monotonic_saturation():
    E, nu = 1000.0, 0.3
    K, G = _kg(E, nu)
    sig0, C, gamma = 5.0, 400.0, 80.0    # X_sat (stress) = C/gamma = 5.0

    def mat(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, C, gamma)

    # push to large plastic strain so the AF backstress reaches its asymptote
    # (X(p) = (C/gamma)(1-e^{-gamma p}); need gamma*p >> 1)
    res = _run(mat, [(0.15, 200)])
    s_end = res[-1][1]
    s_sat = sig0 + C / gamma             # = 10.0
    assert s_end == pytest.approx(s_sat, rel=0.01), (
        f"AF saturation sigma {s_end} != sig0+C/gamma {s_sat}")


# --------------------------------------------------------------------------
# setResponse / getResponse: the internal state (backstress, equivalent
# plastic strain) must be recordable via the element's material response.
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_state_recording():
    K, G = _kg(1000.0, 0.3)
    sig0, C, gamma = 5.0, 400.0, 80.0    # single AF; X_sat(stress)=C/gamma=5

    def mat(t):
        ops.nDMaterial("LadrunoJ2", t, K, G,
                       "-iso", "voce", sig0, 0.0, 0.0, 0.0, "-kin", 1, C, gamma)

    _run(mat, [(0.15, 200)])             # push to near AF saturation (single stdBrick)

    ebar = list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))
    back = list(ops.eleResponse(1, "material", 1, "backStress"))
    assert len(ebar) == 1 and len(back) == 6, "responses not wired through the element"

    # accumulated plastic strain should be ~ (total axial strain - elastic)
    assert ebar[0] > 0.1, f"equivalentPlasticStrain {ebar[0]} not accumulated"
    # backstress is deviatoric (trace ~ 0)
    assert abs(back[0] + back[1] + back[2]) < 1e-6 * max(abs(back[0]), 1.0), "backstress not deviatoric"
    # axial backstress saturates to (2/3)(C/gamma)  (sigma_back=(3/2)alpha_axial=C/gamma)
    assert back[0] == pytest.approx((2.0/3.0) * (C / gamma), rel=0.03), (
        f"axial backstress {back[0]} != (2/3)(C/gamma) {(2.0/3.0)*(C/gamma)}")


# --------------------------------------------------------------------------
# 3D mixed (shear-inclusive) reduce-to-J2Plasticity — exercises the tensor
# machinery the uniaxial tests never touch: dotT's factor-2, the IIDEV6 shear
# block, and the n(x)n shear slots; validates the consistent tangent
# indirectly (Newton must converge on a full 3D plastic BVP).
# --------------------------------------------------------------------------
_NODES3D = {1: (0, 0, 0), 2: (1, 0.1, 0), 3: (1.1, 1, 0.1), 4: (0.05, 0.95, 0),
            5: (0, 0.05, 1), 6: (1, 0, 1.05), 7: (1.05, 1, 1.1), 8: (0, 1, 0.95)}
_TOPLOADS = {5: (8.0, -3.0, 4.0), 6: (-5.0, 6.0, -2.0),
             7: (3.0, 4.0, 7.0), 8: (-4.0, 5.0, -6.0)}   # mixed -> 3D stress w/ shear


def _solve_hex_plastic(mat_fn, nsteps=40, scale=0.9):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _NODES3D.items():
        ops.node(t, *c)
    mat_fn(1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _TOPLOADS.items():
        ops.load(n, fx * scale, fy * scale, fz * scale)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, "3D plastic solve did not converge"
    ops.eleResponse(1, "forces")
    stresses = list(ops.eleResponse(1, "stresses"))[0:6]   # full 6-comp at GP1
    disps = [ops.nodeDisp(n, d) for n in (5, 6, 7, 8) for d in (1, 2, 3)]
    return disps, stresses


@pytest.mark.t1
def test_reduce_to_J2Plasticity_3D_mixed_shear():
    K, G = _kg(1000.0, 0.3)
    sig0, Qinf, b, Hiso = 5.0, 4.0, 30.0, 10.0

    lad = _solve_hex_plastic(lambda t: ops.nDMaterial(
        "LadrunoJ2", t, K, G, "-iso", "voce", sig0, Qinf, b, Hiso, "-kin", 0))
    ref = _solve_hex_plastic(lambda t: ops.nDMaterial(
        "J2Plasticity", t, K, G, sig0, sig0 + Qinf, b, Hiso, 0.0))

    # guard against a vacuous test: the state must actually carry shear stress
    assert any(abs(s) > 1e-3 for s in lad[1][3:6]), "expected nonzero shear stress"
    # all 6 stress components (incl. shear) and all nodal disps must match J2Plasticity
    for a, c in zip(lad[1], ref[1]):
        assert a == pytest.approx(c, rel=2e-6, abs=2e-6), f"stress {a} vs {c}"
    for a, c in zip(lad[0], ref[0]):
        assert a == pytest.approx(c, rel=2e-6, abs=2e-6), f"disp {a} vs {c}"


# --------------------------------------------------------------------------
# dimensional views: a 2D quad reduce-to-J2Plasticity for PlaneStrain and
# PlaneStress (the condensation path). The element calls getCopy(type), so this
# exercises the reduced strain/stress mapping AND the sigma_22=0 nested solve.
# --------------------------------------------------------------------------
_QNODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}


def _solve_quad_plastic(mat_fn, qtype, nsteps=40, scale=0.8):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in _QNODES.items():
        ops.node(t, x, y)
    mat_fn(1)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, qtype, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 6.0 * scale, 0.0)
    ops.load(3, 6.0 * scale, 4.0 * scale)     # mixed -> in-plane shear too
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"{qtype} quad solve did not converge"
    disps = [ops.nodeDisp(n, d) for n in (2, 3, 4) for d in (1, 2)]
    sig = list(ops.eleResponse(1, "stress"))     # all GPs
    return disps, sig


@pytest.mark.t1
@pytest.mark.parametrize("qtype", ["PlaneStrain", "PlaneStress"])
def test_dimensional_view_reduce_to_J2Plasticity(qtype):
    K, G = _kg(1000.0, 0.3)
    sig0, Qinf, b, Hiso = 5.0, 4.0, 30.0, 10.0

    lad = _solve_quad_plastic(lambda t: ops.nDMaterial(
        "LadrunoJ2", t, K, G, "-iso", "voce", sig0, Qinf, b, Hiso, "-kin", 0), qtype)
    ref = _solve_quad_plastic(lambda t: ops.nDMaterial(
        "J2Plasticity", t, K, G, sig0, sig0 + Qinf, b, Hiso, 0.0), qtype)

    # must have actually yielded (guard against a vacuous elastic comparison)
    assert any(abs(d) > sig0 / 1000.0 for d in lad[0]), "expected plastic deformation"
    for a, c in zip(lad[0], ref[0]):
        assert a == pytest.approx(c, rel=2e-6, abs=2e-6), f"{qtype} disp {a} vs {c}"
    for a, c in zip(lad[1], ref[1]):
        assert a == pytest.approx(c, rel=2e-6, abs=2e-6), f"{qtype} stress {a} vs {c}"
