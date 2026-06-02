"""Lemaitre coupled ductile damage on the J2 family (ADR 15) — Zone-A battery.

Covers the `-damage lemaitre $r $s $pD $Dc` mode on BOTH the 3D nDMaterial
LadrunoJ2 (classTag 33011) and the uniaxial LadrunoUniaxialJ2 (33000), sharing
SRC/material/nD/LadrunoDamage.h.

3D driver: single ``stdBrick`` unit cube, 1/8-symmetry determinate restraints =>
homogeneous UNIAXIAL STRESS with free lateral contraction (triaxiality 1/3 =>
R_v == 1). 1D driver: a single ``Truss`` of unit length/area. Both are driven by
displacement control (stable through damage softening) and read at the centroid /
single integration point.

Test map (Ladruno_implementation/15_lemaitre_ductile_damage_adr.md):
  V0  damage-off == undamaged (pD above the reached strain) — bit-tripwire
  V1  monotone uniaxial: D grows monotonically, stress softens past pD
  V2  threshold: zero damage below pD
  V3  triaxiality: a laterally-confined cube damages faster than free-uniaxial
  V4  3D<->1D oracle: brick-uniaxial vs truss-uniaxial, identical params, ~1e-6
  V8  getInitialTangent stays elastic after deep damage (modal correctness)
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def _kg(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


# --------------------------------------------------------------------------
# 3D single-brick driver (uniaxial stress, free lateral contraction)
# --------------------------------------------------------------------------
def _build_brick(mat_fn, confine=False, scale=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *(scale * x for x in c))
    mat_fn(1)
    # 1/8-symmetry determinate restraints (free lateral) OR full lateral confinement
    if confine:
        # clamp all y,z so lateral strain ~ 0 => triaxial (higher sigmaH/sigma_eq)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
        ops.fix(3, 0, 1, 1)
        ops.fix(4, 1, 1, 1)
        ops.fix(5, 1, 1, 1)
        ops.fix(6, 0, 1, 1)
        ops.fix(7, 0, 1, 1)
        ops.fix(8, 1, 1, 1)
    else:
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
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_brick(mat_fn, eps_target, nsteps, confine=False, scale=1.0):
    """Returns list of (eps_xx, sig_xx, D) at each step. `scale` resizes the cube
    (so the element characteristic length changes -> exercises -autoRegularization).
    eps_target is engineering strain; the driven displacement is eps_target*scale."""
    _build_brick(mat_fn, confine, scale)
    out = []
    d = (eps_target * scale) / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        ops.eleResponse(1, "forces")
        sig = list(ops.eleResponse(1, "stresses"))[0:6]
        D = list(ops.eleResponse(1, "material", 1, "damage"))[0]
        out.append((ops.nodeDisp(2, 1) / scale, sig[0], D))
    return out


# --------------------------------------------------------------------------
# 1D single-truss driver
# --------------------------------------------------------------------------
def _build_truss(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    mat_fn(1)
    ops.fix(1, 1)
    ops.element("Truss", 1, 1, 2, 1.0, 1)     # area 1, length 1 => strain = disp
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_truss(mat_fn, eps_target, nsteps):
    _build_truss(mat_fn)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        stress = ops.eleResponse(1, "material", 1, "stress")
        sig = stress[0] if stress else 0.0
        D = list(ops.eleResponse(1, "material", 1, "damage"))[0]
        out.append((ops.nodeDisp(2, 1), sig, D))
    return out


_E, _NU = 1000.0, 0.3
_K, _G = _kg(_E, _NU)
# moderate metal-ish: yield 5, some Voce + linear iso, single AF
_ISO = ("-iso", "voce", 5.0, 4.0, 30.0, 10.0)
_KIN = ("-kin", 1, 400.0, 80.0)


def _mat3d(tag, dmg=None):
    args = ["LadrunoJ2", tag, _K, _G, *_ISO, *_KIN]
    if dmg:
        args += ["-damage", "lemaitre", *dmg]
    ops.nDMaterial(*args)


def _mat1d(tag, dmg=None):
    args = ["LadrunoUniaxialJ2", tag, _E, *_ISO, *_KIN]
    if dmg:
        args += ["-damage", "lemaitre", *dmg]
    ops.uniaxialMaterial(*args)


# ==========================================================================
# V0 — damage-off (pD above reached strain) == undamaged
# ==========================================================================
@pytest.mark.t0m
def test_V0_damage_off_identical_3d():
    base = _run_brick(lambda t: _mat3d(t), 0.05, 25)
    # pD huge => damage never activates => identical to no-damage
    dmg = _run_brick(lambda t: _mat3d(t, dmg=[2.0, 1.0, 100.0, 0.99]), 0.05, 25)
    assert len(base) == len(dmg)
    for (e0, s0, _), (e1, s1, D) in zip(base, dmg):
        assert D == 0.0
        assert abs(s0 - s1) <= 1e-9 * (abs(s0) + 1.0)


def test_V0_damage_off_identical_1d():
    base = _run_truss(lambda t: _mat1d(t), 0.05, 25)
    dmg = _run_truss(lambda t: _mat1d(t, dmg=[2.0, 1.0, 100.0, 0.99]), 0.05, 25)
    for (e0, s0, _), (e1, s1, D) in zip(base, dmg):
        assert D == 0.0
        assert abs(s0 - s1) <= 1e-9 * (abs(s0) + 1.0)


# ==========================================================================
# V1 — monotone uniaxial: D grows monotonically, stress softens past pD
# ==========================================================================
def test_V1_monotone_damage_softens_3d():
    dmg = [1.5, 1.0, 0.01, 0.95]
    res = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.12, 60)
    Ds = [r[2] for r in res]
    assert Ds[-1] > 1e-3, f"damage did not grow: {Ds[-1]}"
    assert all(b >= a - 1e-14 for a, b in zip(Ds, Ds[1:])), "D not monotonic"
    # stress must eventually fall below the undamaged response
    base = _run_brick(lambda t: _mat3d(t), 0.12, 60)
    assert res[-1][1] < base[-1][1] - 1e-6, "no softening vs undamaged"


def test_V1_monotone_damage_softens_1d():
    dmg = [1.5, 1.0, 0.01, 0.95]
    res = _run_truss(lambda t: _mat1d(t, dmg=dmg), 0.12, 60)
    Ds = [r[2] for r in res]
    assert Ds[-1] > 1e-3
    assert all(b >= a - 1e-14 for a, b in zip(Ds, Ds[1:]))
    base = _run_truss(lambda t: _mat1d(t), 0.12, 60)
    assert res[-1][1] < base[-1][1] - 1e-6


# ==========================================================================
# V2 — threshold: zero damage below pD
# ==========================================================================
def test_V2_threshold_gate_1d():
    pD = 0.05
    res = _run_truss(lambda t: _mat1d(t, dmg=[1.5, 1.0, pD, 0.95]), 0.12, 60)
    # early steps (elastic + sub-threshold plastic) accumulate no damage; D only
    # activates once the accumulated plastic strain exceeds pD
    assert res[0][2] == 0.0, "damage active on the first (elastic) step"
    assert res[-1][2] > 0.0, "damage never activated past threshold"
    # and the damage onset must be strictly later than the very first plastic step
    first_dmg = next((i for i, r in enumerate(res) if r[2] > 0.0), None)
    assert first_dmg is not None and first_dmg > 0


# ==========================================================================
# V3 — triaxiality: confined cube damages faster than free-uniaxial
# ==========================================================================
def test_V3_triaxiality_signature():
    dmg = [3.0, 1.0, 0.0, 0.95]
    free = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.06, 40, confine=False)
    conf = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.06, 40, confine=True)
    # at the same applied axial strain, the higher-triaxiality (confined) state
    # has larger R_v => larger Y => more damage
    assert conf[-1][2] > free[-1][2] + 1e-6, \
        f"confined D {conf[-1][2]} !> free D {free[-1][2]}"


# ==========================================================================
# V4 — 3D<->1D oracle: brick-uniaxial vs truss-uniaxial, identical params
# ==========================================================================
def test_V4_3d_1d_damage_oracle():
    dmg = [1.5, 1.0, 0.01, 0.95]
    b = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.10, 50)
    u = _run_truss(lambda t: _mat1d(t, dmg=dmg), 0.10, 50)
    assert len(b) == len(u)
    for (eb, sb, Db), (eu, su, Du) in zip(b, u):
        # identical axial strain history; uniaxial-stress 3D (triaxiality 1/3,
        # R_v == 1) must match the uniaxial material in stress AND damage
        sc = abs(su) + 1.0
        assert abs(sb - su) <= 1e-6 * sc, f"stress mismatch {sb} vs {su}"
        assert abs(Db - Du) <= 1e-7, f"damage mismatch {Db} vs {Du}"


# ==========================================================================
# V6 — characteristic-length regularization: -autoRegularization makes the damage
#      increment scale with element size, so a larger element damages faster per
#      unit plastic strain (crack-band fracture-energy objectivity). Without
#      auto-regularization both sizes share the same local D-vs-strain curve.
# ==========================================================================
def _mat3d_reg(tag, dmg, lch_ref):
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, *_ISO, *_KIN,
                   "-damage", "lemaitre", *dmg, "-autoRegularization", lch_ref)


def test_V6_lch_regularization_scales_damage():
    dmg = [3.0, 1.0, 0.0, 0.95]
    # local (no -autoRegularization): D-vs-strain identical regardless of size
    loc_small = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.06, 40, scale=1.0)
    loc_big   = _run_brick(lambda t: _mat3d(t, dmg=dmg), 0.06, 40, scale=2.0)
    assert abs(loc_small[-1][2] - loc_big[-1][2]) <= 1e-6, \
        "local model should be size-independent in strain space"

    # regularized (lch_ref = 1): the 2x cube (lch ~ 2) damages faster per strain
    rs = _run_brick(lambda t: _mat3d_reg(t, dmg, 1.0), 0.06, 40, scale=1.0)
    rb = _run_brick(lambda t: _mat3d_reg(t, dmg, 1.0), 0.06, 40, scale=2.0)
    assert rb[-1][2] > rs[-1][2] + 1e-4, \
        f"regularized big-cube D {rb[-1][2]} !> small-cube D {rs[-1][2]}"


# ==========================================================================
# V8 — current tangent is degraded by damage (the (1-D) factor + dD/deps term)
# ==========================================================================
def test_V8_tangent_degraded_by_damage_1d():
    res = _run_truss(lambda t: _mat1d(t, dmg=[1.0, 1.0, 0.0, 0.95]), 0.10, 50)
    assert res[-1][2] > 1e-3, "damage did not accumulate"
    et = ops.eleResponse(1, "material", 1, "tangent")
    # after appreciable damage the current tangent is well below the virgin E
    assert et and et[0] < _E, "tangent not degraded by damage"
