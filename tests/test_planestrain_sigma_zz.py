"""Plane-strain out-of-plane stress exposure — getStressZZ() + "stressesPlaneStrain".

Plane-strain continuum elements historically report only [sxx, syy, sxy]; the
constitutive sigma_zz (computed internally under eps_zz = 0) was discarded at
the element/recorder boundary, so full-tensor invariants (von Mises, J2,
triaxiality) were inexact downstream — unrecoverable for nonlinear materials
where sigma_zz != nu*(sxx+syy). This battery gates the fix:

  * material accessor NDMaterial::getStressZZ() (quiet-NaN default),
  * element response "stressesPlaneStrain" = [sxx, syy, sxy, szz] per GP on
    FourNodeQuad / Tri31 / LadrunoQuad / LadrunoCST / BezierTri6.

Gates:
  1. elastic: szz == nu*(sxx+syy) exactly, on every element type,
  2. the first 3 components match the legacy "stresses" response bitwise,
  3. J2 plastic: von Mises built WITH the true szz sits on the yield surface
     (radial return), while the elastic nu*(sxx+syy) reconstruction does NOT —
     the whole point of exposing it,
  4. the PlaneStrainMaterial wrapper (any 3D material under eps_zz=0) reports
     the wrapped 6-vector's exact sigma_zz,
  5. plane-STRESS materials (no override) report NaN, never a misleading 0.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU = 30_000.0, 0.3
DELTA = 1.0e-4  # prescribed top settlement (elastic runs)


def _mises(sxx, syy, sxy, szz):
    m = (sxx + syy + szz) / 3.0
    dxx, dyy, dzz = sxx - m, syy - m, szz - m
    return math.sqrt(1.5 * (dxx * dxx + dyy * dyy + dzz * dzz + 2.0 * sxy * sxy))


def _uniaxial_strain_column(nodes, fix_x_all=True):
    """Common rig: eps_xx = 0 (all x fixed), bottom y fixed, top y prescribed."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in nodes.items():
        ops.node(tag, x, y)
        ops.fix(tag, 1 if fix_x_all else 0, 1 if y == 0.0 else 0)


def _drive(top_nodes, delta, nsteps=10):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top_nodes:
        ops.sp(n, 2, -delta)
    # Penalty (not Transformation): the rig prescribes EVERY free DOF, and a
    # fully-condensed system gives FullGenLinSOE a fatal zero-size vectX
    ops.constraints("Penalty", 1.0e15, 1.0e15)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0


def _check_elastic_gps(ele_tag, ncomp_legacy, rel=1.0e-9):
    """Gates 1 + 2 on one element: 4-comp response vs legacy + elastic szz."""
    s3 = ops.eleResponse(ele_tag, "stresses")
    s4 = ops.eleResponse(ele_tag, "stressesPlaneStrain")
    ngp = ncomp_legacy // 3
    assert len(s3) == ncomp_legacy
    assert len(s4) == 4 * ngp
    for g in range(ngp):
        sxx, syy, sxy, szz = s4[4 * g : 4 * g + 4]
        # first three components are the SAME material query — bitwise equal
        assert (sxx, syy, sxy) == tuple(s3[3 * g : 3 * g + 3])
        # elastic plane strain: szz = nu*(sxx+syy), and nonzero here
        ref = NU * (sxx + syy)
        assert abs(szz - ref) <= rel * abs(ref)
        assert abs(szz) > 0.0


# ---------------------------------------------------------------------------
# elastic szz on every element type
# ---------------------------------------------------------------------------

def test_quad_elastic_szz():
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    _drive([3, 4], DELTA)
    _check_elastic_gps(1, 12)


def test_tri31_elastic_szz():
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("Tri31", 1, 1, 2, 3, 1.0, "PlaneStrain", 1)
    _drive([3], DELTA)
    _check_elastic_gps(1, 3)


def test_ladruno_quad_elastic_szz():
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0)
    _drive([3, 4], DELTA)
    _check_elastic_gps(1, 12)


def test_ladruno_cst_elastic_szz():
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoCST", 1, 1, 2, 3, 1, "-type", "PlaneStrain", "-thick", 1.0)
    _drive([3], DELTA)
    _check_elastic_gps(1, 3)


def test_bezier_tri6_elastic_szz():
    # straight-sided quadratic triangle: corners + edge midpoints
    nodes = {
        1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0),
        4: (0.5, 0.0), 5: (0.5, 0.5), 6: (0.0, 0.5),
    }
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("BezierTri6", 1, 1, 2, 3, 4, 5, 6, 1.0, "PlaneStrain", 1)
    _drive([3], DELTA)
    _check_elastic_gps(1, 9)


# ---------------------------------------------------------------------------
# nonlinear: the true szz is NOT the elastic reconstruction (the point)
# ---------------------------------------------------------------------------

def test_j2_plastic_szz_on_yield_surface():
    """Radial return keeps ||dev sigma|| on the yield surface — but only when
    von Mises is built with the TRUE sigma_zz. The elastic nu*(sxx+syy)
    reconstruction misses it once plastic flow starts."""
    K, G, SIG0 = 25_000.0, 11_500.0, 30.0
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    # perfect plasticity: sig0 == sigInf, H = 0 -> q == SIG0 always
    ops.nDMaterial("J2Plasticity", 1, K, G, SIG0, SIG0, 0.0, 0.0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    _drive([3, 4], 0.02, nsteps=40)  # oedometric push far past yield

    nu3d = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    s4 = ops.eleResponse(1, "stressesPlaneStrain")
    for g in range(4):
        sxx, syy, sxy, szz = s4[4 * g : 4 * g + 4]
        assert math.isfinite(szz)
        seq_true = _mises(sxx, syy, sxy, szz)
        seq_elastic_guess = _mises(sxx, syy, sxy, nu3d * (sxx + syy))
        # on the yield surface with the true szz ...
        assert seq_true == pytest.approx(SIG0, rel=1.0e-8)
        # ... and measurably OFF it with the elastic reconstruction
        assert abs(seq_elastic_guess - SIG0) > 0.05 * SIG0


def test_planestrain_wrapper_szz_matches_dedicated():
    """PlaneStrainMaterial wrapping the 3D J2 must expose the same exact szz
    as the dedicated J2PlaneStrain implementation (identical return mapping)."""
    K, G, SIG0 = 25_000.0, 11_500.0, 30.0

    def run(mat_setup):
        nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
        _uniaxial_strain_column(nodes)
        mat_tag = mat_setup()
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", mat_tag)
        _drive([3, 4], 0.02, nsteps=40)
        return ops.eleResponse(1, "stressesPlaneStrain")

    def direct():
        ops.nDMaterial("J2Plasticity", 1, K, G, SIG0, SIG0, 0.0, 0.0)
        return 1

    def wrapped():
        ops.nDMaterial("J2Plasticity", 1, K, G, SIG0, SIG0, 0.0, 0.0)
        ops.nDMaterial("PlaneStrain", 2, 1)
        return 2

    a, b = run(direct), run(wrapped)
    assert a == pytest.approx(b, rel=1.0e-10)
    assert all(math.isfinite(v) for v in b)


def test_ladruno_j2_szz():
    """Fork LadrunoJ2 in plane strain exposes its internal stress6[2]."""
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    K, G, SIGY = 25_000.0, 11_500.0, 30.0
    # perfect plasticity: Qinf = b = Hiso = 0 -> q == SIGY always
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", SIGY, 0.0, 0.0, 0.0)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0)
    _drive([3, 4], 0.02, nsteps=40)
    s4 = ops.eleResponse(1, "stressesPlaneStrain")
    for g in range(4):
        sxx, syy, sxy, szz = s4[4 * g : 4 * g + 4]
        assert math.isfinite(szz)
        assert _mises(sxx, syy, sxy, szz) == pytest.approx(SIGY, rel=1.0e-6)


# ---------------------------------------------------------------------------
# NaN contract: a material that doesn't expose szz reports NaN, not 0
# ---------------------------------------------------------------------------

def test_plane_stress_reports_nan():
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    _uniaxial_strain_column(nodes)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
    _drive([3, 4], DELTA)
    s4 = ops.eleResponse(1, "stressesPlaneStrain")
    assert len(s4) == 16
    for g in range(4):
        assert math.isnan(s4[4 * g + 3])
        # in-plane components still real
        assert all(math.isfinite(v) for v in s4[4 * g : 4 * g + 3])
