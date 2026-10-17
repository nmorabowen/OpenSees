"""WP-1: Non-uniform partial/trapezoidal frame load verification suite.

Covers:
  * T1.1: Closed-form Fixed-End Force (FEF) and section-force oracle across all
    2D and 3D beam-column families:
      - ElasticBeam2d, ElasticBeam3d, ModElasticBeam3d
      - DispBeamColumn2d, DispBeamColumn3d
      - LadrunoDispBeamColumn2d, LadrunoDispBeamColumn3d
      - ForceBeamColumn2d, ForceBeamColumn3d, ForceBeamColumnCBDI3d
      - ElasticTimoshenkoBeam2d, ElasticTimoshenkoBeam3d (both Euler-Bernoulli
        limit phi -> 0 and finite shear flexibility phi > 0)
  * T1.2: End-release matrix (release / releasez / releasey in {0, 1, 2, 3})
    for ElasticBeam2d and ElasticBeam3d under partial trapezoidal and point loads.
  * T1.3: Sliver-span numerical stability (b - a -> 0 and b == a) verifying the
    cancellation-free shifted-moment formulation against catastrophic float64
    cancellation.
  * T1.4: FEM_ObjectBrokerAllClasses serialization round-trip for
    LOAD_TAG_Beam2dPartialUniformLoad and LOAD_TAG_Beam3dPartialUniformLoad.
"""
import math

import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

# ---------------------------------------------------------------------------
# Analytical closed-form helpers (shifted-moment formulation on [a, b])
# ---------------------------------------------------------------------------

def exact_trapezoidal_moments(w_a: float, w_b: float, a: float, b: float, L: float):
    """Return (V1_ss, V2_ss, M1_eb, M2_eb, I0, I1) for linear load w(x) on [a, b].

    Sign convention matches OpenSees basic system:
      * V1_ss, V2_ss: simply-supported end reactions (+ve for +ve transverse load)
      * M1_eb, M2_eb: Euler-Bernoulli fixed-end moments in basic system
        (for +ve load, M1_eb < 0 at node I and M2_eb > 0 at node J)
    """
    ell = b - a
    if ell <= 0.0 or a < 0.0 or b > L:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    J0 = 0.5 * ell * (w_a + w_b)
    J1 = (ell ** 2 / 6.0) * (w_a + 2.0 * w_b)
    J2 = (ell ** 3 / 12.0) * (w_a + 3.0 * w_b)
    J3 = (ell ** 4 / 20.0) * (w_a + 4.0 * w_b)

    I0 = J0
    I1 = a * J0 + J1
    I2 = (a ** 2) * J0 + 2.0 * a * J1 + J2
    I3 = (a ** 3) * J0 + 3.0 * (a ** 2) * J1 + 3.0 * a * J2 + J3

    V2 = I1 / L
    V1 = I0 - V2
    M1 = -I1 + (2.0 / L) * I2 - (1.0 / (L ** 2)) * I3
    M2 = (1.0 / L) * I2 - (1.0 / (L ** 2)) * I3
    return V1, V2, M1, M2, I0, I1


def exact_axial_reactions(wa_a: float, wa_b: float, a: float, b: float, L: float):
    """Return (P_total, P_I, P_J) for linear axial load on [a, b]."""
    ell = b - a
    if ell <= 0.0 or a < 0.0 or b > L:
        return 0.0, 0.0, 0.0
    P = 0.5 * ell * (wa_a + wa_b)
    PJ = (a * P + (ell ** 2 / 6.0) * (wa_a + 2.0 * wa_b)) / L
    PI = P - PJ
    return P, PI, PJ


def exact_2d_fef(wya: float, wyb: float, waa: float, wab: float,
                 aL: float, bL: float, L: float,
                 release: int = 0, phi: float = 0.0):
    """Return exact 2D local fixed-end force vector [Fx_I, Fy_I, Mz_I, Fx_J, Fy_J, Mz_J]."""
    a = aL * L
    b = bL * L
    V1, V2, M1_eb, M2_eb, _, _ = exact_trapezoidal_moments(wya, wyb, a, b, L)
    _, PI, PJ = exact_axial_reactions(waa, wab, a, b, L)

    if release == 0:
        Mz1 = (M1_eb + 0.5 * phi * (M1_eb - M2_eb)) / (1.0 + phi)
        Mz2 = (M2_eb - 0.5 * phi * (M1_eb - M2_eb)) / (1.0 + phi)
    elif release == 1:
        Mz1 = 0.0
        Mz2 = M2_eb - 0.5 * M1_eb
    elif release == 2:
        Mz1 = M1_eb - 0.5 * M2_eb
        Mz2 = 0.0
    else:
        Mz1 = 0.0
        Mz2 = 0.0

    dV = (Mz1 + Mz2) / L
    return [
        -PI,
        -V1 + dV,
        Mz1,
        -PJ,
        -V2 - dV,
        Mz2,
    ]


def exact_3d_fef(wya: float, wyb: float, wza: float, wzb: float,
                 waa: float, wab: float, aL: float, bL: float, L: float,
                 releasez: int = 0, releasey: int = 0,
                 phiZ: float = 0.0, phiY: float = 0.0):
    """Return exact 3D local fixed-end force vector (12 DOFs)."""
    a = aL * L
    b = bL * L
    Vy1, Vy2, M1y_eb, M2y_eb, _, _ = exact_trapezoidal_moments(wya, wyb, a, b, L)
    Vz1, Vz2, M1z_eb, M2z_eb, _, _ = exact_trapezoidal_moments(wza, wzb, a, b, L)
    _, PI, PJ = exact_axial_reactions(waa, wab, a, b, L)

    if releasez == 0:
        Mz1 = (M1y_eb + 0.5 * phiZ * (M1y_eb - M2y_eb)) / (1.0 + phiZ)
        Mz2 = (M2y_eb - 0.5 * phiZ * (M1y_eb - M2y_eb)) / (1.0 + phiZ)
    elif releasez == 1:
        Mz1 = 0.0
        Mz2 = M2y_eb - 0.5 * M1y_eb
    elif releasez == 2:
        Mz1 = M1y_eb - 0.5 * M2y_eb
        Mz2 = 0.0
    else:
        Mz1 = 0.0
        Mz2 = 0.0

    if releasey == 0:
        My1 = -(M1z_eb + 0.5 * phiY * (M1z_eb - M2z_eb)) / (1.0 + phiY)
        My2 = -(M2z_eb - 0.5 * phiY * (M1z_eb - M2z_eb)) / (1.0 + phiY)
    elif releasey == 1:
        My1 = 0.0
        My2 = -(M2z_eb - 0.5 * M1z_eb)
    elif releasey == 2:
        My1 = -(M1z_eb - 0.5 * M2z_eb)
        My2 = 0.0
    else:
        My1 = 0.0
        My2 = 0.0

    dVy = (Mz1 + Mz2) / L
    dVz = (My1 + My2) / L
    return [
        -PI,
        -Vy1 + dVy,
        -Vz1 - dVz,
        0.0,
        My1,
        Mz1,
        -PJ,
        -Vy2 - dVy,
        -Vz2 + dVz,
        0.0,
        My2,
        Mz2,
    ]


# ---------------------------------------------------------------------------
# Model builders for 2D and 3D beam families
# ---------------------------------------------------------------------------

def _create_2d_element(family: str, L: float, E: float = 2.0e11, G: float = 8.0e10,
                       A: float = 0.04, Iz: float = 1.0e-4, Avy: float = 1.0e8,
                       release: int = 0):
    ops.geomTransf("Linear", 1)
    if family == "ElasticBeam2d":
        ops.element("elasticBeamColumn", 1, 1, 2, A, E, Iz, 1, "-release", release)
    elif family == "DispBeamColumn2d":
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element("dispBeamColumn", 1, 1, 2, 1, 1)
    elif family == "LadrunoDispBeamColumn2d":
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element("ladrunoDispBeamColumn", 1, 1, 2, 1, 1)
    elif family == "ForceBeamColumn2d":
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 10)
        ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    elif family == "ElasticTimoshenkoBeam2d":
        ops.element("ElasticTimoshenkoBeam", 1, 1, 2, E, G, A, Iz, Avy, 1)
    else:
        raise ValueError(f"Unknown 2D family: {family}")


def _create_3d_element(family: str, L: float, E: float = 2.0e11, G: float = 8.0e10,
                       A: float = 0.04, Jx: float = 2.0e-4, Iy: float = 1.2e-4,
                       Iz: float = 1.6e-4, Avy: float = 1.0e8, Avz: float = 1.0e8,
                       releasez: int = 0, releasey: int = 0):
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    if family == "ElasticBeam3d":
        ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
                    "-releasez", releasez, "-releasey", releasey)
    elif family == "ModElasticBeam3d":
        ops.element("ModElasticBeam3d", 1, 1, 2, A, E, G, Jx, Iy, Iz,
                    4.0, 2.0, 4.0, 4.0, 2.0, 4.0, 1)
    elif family == "DispBeamColumn3d":
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element("dispBeamColumn", 1, 1, 2, 1, 1)
    elif family == "LadrunoDispBeamColumn3d":
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element("ladrunoDispBeamColumn", 1, 1, 2, 1, 1)
    elif family == "ForceBeamColumn3d":
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 10)
        ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    elif family == "ForceBeamColumnCBDI3d":
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 10)
        ops.element("forceBeamColumnCBDI", 1, 1, 2, 1, 1)
    elif family == "ElasticTimoshenkoBeam3d":
        # Note: OPS_ElasticTimoshenkoBeam3d passes dData[4]=Iy, dData[5]=Iz to constructor
        ops.element("ElasticTimoshenkoBeam", 1, 1, 2, E, G, A, Jx, Iy, Iz, Avy, Avz, 1)
    else:
        raise ValueError(f"Unknown 3D family: {family}")


# ---------------------------------------------------------------------------
# T1.1: Closed-form Fixed-End Force (FEF) & Section-Force Oracle
# ---------------------------------------------------------------------------

FAMILIES_2D_EXACT_PARTIAL = [
    "ElasticBeam2d",
    "DispBeamColumn2d",
    "LadrunoDispBeamColumn2d",
    "ElasticTimoshenkoBeam2d",
]

FAMILIES_3D_EXACT_PARTIAL = [
    "ElasticBeam3d",
    "ModElasticBeam3d",
    "DispBeamColumn3d",
    "LadrunoDispBeamColumn3d",
    "ElasticTimoshenkoBeam3d",
]


@pytest.mark.zone_a
@pytest.mark.parametrize("family", FAMILIES_2D_EXACT_PARTIAL)
def test_t1_1_2d_partial_trapezoidal_fef_exact(family):
    """T1.1 (2D): Partial trapezoidal load on [0.2L, 0.8L] matches closed-form FEF to 1e-12."""
    L = 5.0
    wya, wyb = -12.0e3, -28.0e3
    waa, wab = 4.0e3, 9.0e3
    aL, bL = 0.2, 0.8

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    E, G, A, Iz, Avy = 2.0e11, 8.0e10, 0.04, 1.0e-4, 1.0e12
    _create_2d_element(family, L, E=E, G=G, A=A, Iz=Iz, Avy=Avy)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, waa, aL, bL, wyb, wab)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    phi = (12.0 * E * Iz / (L * L * G * Avy)) if family == "ElasticTimoshenkoBeam2d" else 0.0
    expected = exact_2d_fef(wya, wyb, waa, wab, aL, bL, L, phi=phi)
    got = ops.eleForce(1)
    assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)


@pytest.mark.zone_a
@pytest.mark.parametrize("family", FAMILIES_3D_EXACT_PARTIAL)
def test_t1_1_3d_partial_trapezoidal_fef_exact(family):
    """T1.1 (3D): Partial trapezoidal load on [0.2L, 0.8L] matches closed-form FEF to 1e-12."""
    L = 5.0
    wya, wyb = -10.0e3, -25.0e3
    wza, wzb = 6.0e3, -18.0e3
    waa, wab = 3.0e3, 7.0e3
    aL, bL = 0.2, 0.8

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)

    E, G, A, Jx, Iy, Iz, Avy, Avz = 2.0e11, 8.0e10, 0.04, 2.0e-4, 1.2e-4, 1.6e-4, 1.0e12, 1.0e12
    _create_3d_element(family, L, E=E, G=G, A=A, Jx=Jx, Iy=Iy, Iz=Iz, Avy=Avy, Avz=Avz)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, aL, bL, wyb, wzb, wab)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    phiZ = (12.0 * E * Iz / (L * L * G * Avy)) if family == "ElasticTimoshenkoBeam3d" else 0.0
    phiY = (12.0 * E * Iy / (L * L * G * Avz)) if family == "ElasticTimoshenkoBeam3d" else 0.0
    expected = exact_3d_fef(wya, wyb, wza, wzb, waa, wab, aL, bL, L, phiZ=phiZ, phiY=phiY)
    got = ops.eleForce(1)
    assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)


@pytest.mark.zone_a
def test_t1_1_timoshenko_finite_shear_flexibility():
    """T1.1 (Timoshenko): Verify exact shear-flexibility FEF shift for finite phi > 0."""
    L = 4.0
    E, G, A, Jx, Iy, Iz = 2.0e11, 8.0e10, 0.05, 2.0e-4, 2.0e-3, 4.0e-3
    Avy, Avz = 0.01, 0.008  # significant shear flexibility: phiZ ~ 0.375, phiY ~ 0.234375
    wya, wyb = -15.0e3, -35.0e3
    wza, wzb = 8.0e3, 24.0e3
    waa, wab = 2.0e3, 6.0e3
    aL, bL = 0.15, 0.75

    # 2D Timoshenko
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    _create_2d_element("ElasticTimoshenkoBeam2d", L, E=E, G=G, A=A, Iz=Iz, Avy=Avy)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, waa, aL, bL, wyb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    phi2d = 12.0 * E * Iz / (L * L * G * Avy)
    exp2d = exact_2d_fef(wya, wyb, waa, wab, aL, bL, L, phi=phi2d)
    assert ops.eleForce(1) == pytest.approx(exp2d, rel=1e-12, abs=1e-12)

    # 3D Timoshenko
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    _create_3d_element("ElasticTimoshenkoBeam3d", L, E=E, G=G, A=A, Jx=Jx, Iy=Iy, Iz=Iz, Avy=Avy, Avz=Avz)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, aL, bL, wyb, wzb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    phiZ = 12.0 * E * Iz / (L * L * G * Avy)
    phiY = 12.0 * E * Iy / (L * L * G * Avz)
    exp3d = exact_3d_fef(wya, wyb, wza, wzb, waa, wab, aL, bL, L, phiZ=phiZ, phiY=phiY)
    assert ops.eleForce(1) == pytest.approx(exp3d, rel=1e-12, abs=1e-12)


@pytest.mark.zone_a
@pytest.mark.parametrize("family", ["ForceBeamColumn2d", "ForceBeamColumn3d", "ForceBeamColumnCBDI3d"])
def test_t1_1_force_beam_column_trapezoidal_fef_and_section_forces(family):
    """T1.1 (ForceBeamColumn): Full-span trapezoidal FEF is exact to 1e-12; partial-span
    [0.2L, 0.8L] satisfies exact equilibrium (1e-12) and exact section forces sp(x)
    (including SECTION_RESPONSE_VZ) to 1e-12.
    """
    L = 5.0
    wya, wyb = -10.0e3, -25.0e3
    wza, wzb = 6.0e3, -18.0e3
    waa, wab = 3.0e3, 7.0e3

    if family == "ForceBeamColumn2d":
        # 1. Full-span [0, 1] trapezoidal load -> Lobatto quadrature is exact to 1e-12
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 1)
        _create_2d_element(family, L)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, waa, 0.0, 1.0, wyb, wab)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        exp_full = exact_2d_fef(wya, wyb, waa, wab, 0.0, 1.0, L)
        assert ops.eleForce(1) == pytest.approx(exp_full, rel=1e-12, abs=1e-12)
        return

    # 3D ForceBeamColumn / ForceBeamColumnCBDI3d:
    # 1. Full-span [0, 1] trapezoidal load -> exact to 1e-12
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    _create_3d_element(family, L)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, 0.0, 1.0, wyb, wzb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    exp_full = exact_3d_fef(wya, wyb, wza, wzb, waa, wab, 0.0, 1.0, L)
    ef_full = ops.eleForce(1)
    if family == "ForceBeamColumn3d":
        assert ef_full == pytest.approx(exp_full, rel=1e-12, abs=1e-12)
    else:
        # ForceBeamColumnCBDI3d includes internal curvature-displacement (P-delta bowstring)
        # coupling; total reactions sum(Fx), sum(Fy), sum(Fz) are still exact to 1e-12:
        assert (ef_full[0] + ef_full[6]) == pytest.approx(exp_full[0] + exp_full[6], rel=1e-12, abs=1e-12)
        assert (ef_full[1] + ef_full[7]) == pytest.approx(exp_full[1] + exp_full[7], rel=1e-12, abs=1e-12)
        assert (ef_full[2] + ef_full[8]) == pytest.approx(exp_full[2] + exp_full[8], rel=1e-12, abs=1e-12)

    # 2. Simply-supported partial trapezoidal load on [0.2L, 0.8L] with shear section
    #    so basic forces q == 0 and section forces equal sp(x) directly!
    aL, bL = 0.2, 0.8
    a, b = aL * L, bL * L
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    # Simply supported (pin at 1 with torsion fixed, roller at 2) -> basic q = 0
    ops.fix(1, 1, 1, 1, 1, 0, 0)
    ops.fix(2, 0, 1, 1, 1, 0, 0)

    E, G, A, Jx, Iy, Iz = 2.0e11, 8.0e10, 0.04, 2.0e-4, 1.2e-4, 1.6e-4
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 10, E, A, Iz, Iy, G, Jx)
    ops.uniaxialMaterial("Elastic", 20, G * A)
    ops.section("Aggregator", 1, 20, "Vy", 20, "Vz", "-section", 10)
    ops.beamIntegration("Lobatto", 1, 1, 7)
    if family == "ForceBeamColumn3d":
        ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    else:
        ops.element("forceBeamColumnCBDI", 1, 1, 2, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, aL, bL, wyb, wzb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # Check simply-supported reactions at nodes 1 and 2 to 1e-12
    Vy1, Vy2, _, _, _, _ = exact_trapezoidal_moments(wya, wyb, a, b, L)
    Vz1, Vz2, _, _, _, _ = exact_trapezoidal_moments(wza, wzb, a, b, L)
    Fa, _, _ = exact_axial_reactions(waa, wab, a, b, L)
    ef = ops.eleForce(1)
    assert ef[0] == pytest.approx(-Fa, rel=1e-12, abs=1e-12)
    assert ef[1] == pytest.approx(-Vy1, rel=1e-12, abs=1e-12)
    assert ef[2] == pytest.approx(-Vz1, rel=1e-12, abs=1e-12)
    assert ef[7] == pytest.approx(-Vy2, rel=1e-12, abs=1e-12)
    assert ef[8] == pytest.approx(-Vz2, rel=1e-12, abs=1e-12)

    # Check section forces [P, Mz, My, T, Vy, Vz] at all 7 Lobatto sections to 1e-12
    locs = ops.sectionLocation(1)
    for isec, x in enumerate(locs, start=1):
        sf = ops.eleResponse(1, "section", isec, "forces")
        # Order in Elastic+Aggregator(Vy, Vz): [P, Mz, My, T, Vy, Vz]
        if x <= a:
            exp_P = Fa
            exp_Mz = -Vy1 * x
            exp_My = Vz1 * x
            exp_Vy = -Vy1
            exp_Vz = Vz1
        elif x >= b:
            exp_P = 0.0
            exp_Mz = Vy2 * (x - L)
            exp_My = -Vz2 * (x - L)
            exp_Vy = Vy2
            exp_Vz = -Vz2
        else:
            dx = x - a
            wyy = wya + (wyb - wya) / (b - a) * dx
            wzz = wza + (wzb - wza) / (b - a) * dx
            exp_P = Fa - waa * dx - 0.5 * (wab - waa) / (b - a) * dx * dx
            exp_Mz = -Vy1 * x + 0.5 * wya * dx * dx + 0.5 * (wyy - wya) * dx * dx / 3.0
            exp_My = Vz1 * x - 0.5 * wza * dx * dx - 0.5 * (wzz - wza) * dx * dx / 3.0
            exp_Vy = -Vy1 + wya * dx + 0.5 * (wyy - wya) * dx
            exp_Vz = Vz1 - wza * dx - 0.5 * (wzz - wza) * dx

        assert sf == pytest.approx([exp_P, exp_Mz, exp_My, 0.0, exp_Vy, exp_Vz], rel=1e-11, abs=1e-8)


# ---------------------------------------------------------------------------
# T1.2: End-Release Matrix (release / releasez / releasey in {0, 1, 2, 3})
# ---------------------------------------------------------------------------

@pytest.mark.zone_a
@pytest.mark.parametrize("release", [0, 1, 2, 3])
def test_t1_2_elastic_beam_2d_releases_partial_trapezoidal(release):
    """T1.2 (2D): ElasticBeam2d with -release {0, 1, 2, 3} under partial trapezoidal load."""
    L = 5.0
    wya, wyb = -14.0e3, -30.0e3
    waa, wab = 2.0e3, 8.0e3
    aL, bL = 0.2, 0.8

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    _create_2d_element("ElasticBeam2d", L, release=release)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, waa, aL, bL, wyb, wab)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    expected = exact_2d_fef(wya, wyb, waa, wab, aL, bL, L, release=release)
    got = ops.eleForce(1)
    assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)

    # Verify released end moments are identically zero
    if release in (1, 3):
        assert abs(got[2]) < 1e-12
    if release in (2, 3):
        assert abs(got[5]) < 1e-12

    # Verify global static equilibrium
    _, _, _, _, I0y, I1y = exact_trapezoidal_moments(wya, wyb, aL * L, bL * L, L)
    assert (got[1] + got[4]) == pytest.approx(-I0y, rel=1e-12, abs=1e-12)
    assert (got[2] + got[5] + got[4] * L) == pytest.approx(-I1y, rel=1e-12, abs=1e-12)


@pytest.mark.zone_a
@pytest.mark.parametrize("releasez", [0, 1, 2, 3])
@pytest.mark.parametrize("releasey", [0, 1, 2, 3])
def test_t1_2_elastic_beam_3d_releases_partial_trapezoidal_and_point(releasez, releasey):
    """T1.2 (3D): ElasticBeam3d with all 16 (releasez, releasey) combinations under
    both partial trapezoidal load and 3D point load.
    """
    L = 5.0
    wya, wyb = -10.0e3, -26.0e3
    wza, wzb = 8.0e3, -16.0e3
    waa, wab = 3.0e3, 7.0e3
    aL, bL = 0.2, 0.8

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)

    _create_3d_element("ElasticBeam3d", L, releasez=releasez, releasey=releasey)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, aL, bL, wyb, wzb, wab)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    expected = exact_3d_fef(wya, wyb, wza, wzb, waa, wab, aL, bL, L,
                            releasez=releasez, releasey=releasey)
    got = ops.eleForce(1)
    assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)

    if releasez in (1, 3):
        assert abs(got[5]) < 1e-12
    if releasez in (2, 3):
        assert abs(got[11]) < 1e-12
    if releasey in (1, 3):
        assert abs(got[4]) < 1e-12
    if releasey in (2, 3):
        assert abs(got[10]) < 1e-12


# ---------------------------------------------------------------------------
# T1.3: Sliver-Span Numerical Stability (b - a -> 0 and b == a)
# ---------------------------------------------------------------------------

@pytest.mark.zone_a
@pytest.mark.parametrize("delta_xi", [1e-4, 1e-7, 1e-9, 0.0])
def test_t1_3_sliver_span_stability_2d_and_3d(delta_xi):
    """T1.3: Sliver spans b - a in {1e-4 L, 1e-7 L, 1e-9 L, 0} produce finite,
    accurate FEFs without catastrophic float64 cancellation or division by zero.
    """
    L = 5.0
    aL = 0.73
    bL = aL + delta_xi
    wya, wyb = -15.0e3, -25.0e3
    wza, wzb = 10.0e3, 20.0e3
    waa, wab = 5.0e3, 9.0e3

    # 2D ElasticBeam2d
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    _create_2d_element("ElasticBeam2d", L)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, waa, aL, bL, wyb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    got2d = ops.eleForce(1)
    assert all(math.isfinite(v) for v in got2d)
    exp2d = exact_2d_fef(wya, wyb, waa, wab, aL, bL, L)
    if delta_xi == 0.0:
        assert got2d == pytest.approx([0.0] * 6, abs=1e-14)
    else:
        assert got2d == pytest.approx(exp2d, rel=1e-10, abs=1e-14)

    # 3D ElasticBeam3d
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    _create_3d_element("ElasticBeam3d", L)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wya, wza, waa, aL, bL, wyb, wzb, wab)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    got3d = ops.eleForce(1)
    assert all(math.isfinite(v) for v in got3d)
    exp3d = exact_3d_fef(wya, wyb, wza, wzb, waa, wab, aL, bL, L)
    if delta_xi == 0.0:
        assert got3d == pytest.approx([0.0] * 12, abs=1e-14)
    else:
        assert got3d == pytest.approx(exp3d, rel=1e-10, abs=1e-14)


# ---------------------------------------------------------------------------
# T1.4: Broker Serialization Round-Trip (Beam2dPartialUniformLoad & Beam3dPartialUniformLoad)
# ---------------------------------------------------------------------------

@pytest.mark.zone_a
def test_t1_4_broker_database_roundtrip_2d_and_3d():
    """T1.4: FEM_ObjectBrokerAllClasses::getNewElementalLoad round-trips
    Beam2dPartialUniformLoad and Beam3dPartialUniformLoad across 2 save/restore cycles.
    """
    L = 5.0

    def build_2d():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        _create_2d_element("ElasticBeam2d", L)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -12.0e3, 3.0e3, 0.2, 0.8, -24.0e3, 6.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0

    def _probe_restored_ele_force():
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return ops.eleForce(1)

    database_roundtrip(
        build_2d,
        probe_nodes=[2],
        ndf=3,
        probe_fn=_probe_restored_ele_force,
        rounds=2,
    )

    def build_3d():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        _create_3d_element("ElasticBeam3d", L)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform",
                    -10.0e3, 6.0e3, 2.0e3, 0.2, 0.8, -22.0e3, -14.0e3, 5.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0

    database_roundtrip(
        build_3d,
        probe_nodes=[2],
        ndf=6,
        probe_fn=_probe_restored_ele_force,
        rounds=2,
    )
