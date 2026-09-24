"""Dynamic Rayleigh + inertia regression for LadrunoDispBeamColumn (2D/3D) -- WP-116.

LadrunoDispBeamColumn2d/3d accumulated inertia and the Rayleigh force into
the class-shared static P (the #562 pattern, LEDGER_quirks "MUST snapshot the
shared static `resid`"). WP-116 traced every path from
Element::getRayleighDampingForces() -- nothing reachable writes P, so it was
safe -- and converted both classes to a function-local snapshot anyway, with
the same operation order (bit-identical). No transient Rayleigh test existed.

Gate: differential against elasticBeamColumn. With an elastic section and
Legendre integration the displacement-based element IS the elastic beam, and
both use the same lumped (rho*L/2) or consistent (rho*L/420) mass, so every
displacement history must agree. Legs:
  {alphaM, betaK, betaK0} x
    lumped element -mass   x {step tip load, UniformExcitation}
    nodal mass             x {step, ground}
    consistent -cMass      x {step}
The lumped/nodal oracle carries the masses as NODAL masses (vanilla
ElasticBeam2d with element -mass subtracts the ground-motion Q twice --
LEDGER_quirks "ElasticBeam2d subtracts the ground-motion load Q TWICE").
The -cMass oracle is elasticBeamColumn -cMass on step loads only, for the
same reason. Fails in both directions: dropped inertia, dropped Rayleigh,
or a wrong -Q term break the agreement.

Break-on-purpose evidence: WP-116 plan doc.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

XI = 0.05
L_BEAM, A_SEC, E_MOD, I_SEC = 3.0, 0.02, 2.0e5, 1.5e-4
G_MOD, J_TOR = 8.0e4, 3.0e-4
RHO_L = 0.8
NEL = 4


def _rayleigh_for(kind, w1, xi):
    if kind == "alphaM":
        return (2.0 * xi * w1, 0.0, 0.0, 0.0)
    if kind == "betaK":
        return (0.0, 2.0 * xi / w1, 0.0, 0.0)
    if kind == "betaK0":
        return (0.0, 0.0, 2.0 * xi / w1, 0.0)
    raise ValueError(kind)


def _mesh(dim):
    ops.wipe()
    if dim == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E_MOD, A_SEC, I_SEC)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.section("Elastic", 1, E_MOD, A_SEC, I_SEC, I_SEC, G_MOD, J_TOR)
    ops.beamIntegration("Legendre", 1, 1, 3)
    h = L_BEAM / NEL
    for i in range(NEL + 1):
        ops.node(i + 1, *((i * h, 0.0) if dim == 2 else (i * h, 0.0, 0.0)))
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))
    return h


def _elements(dim, elem, mass):
    """mass: 'lumped' / 'consistent' element mass, or 'nodal' (none on the element)."""
    for e in range(1, NEL + 1):
        if elem == "ldbc":
            args = ["LadrunoDispBeamColumn", e, e, e + 1, 1, 1]
        else:
            sec = (A_SEC, E_MOD, I_SEC) if dim == 2 else (A_SEC, E_MOD, G_MOD, J_TOR, I_SEC, I_SEC)
            args = ["elasticBeamColumn", e, e, e + 1, *sec, 1]
        if mass in ("lumped", "consistent"):
            args += ["-mass", RHO_L]
        if mass == "consistent":
            args += ["-cMass"]
        ops.element(*args)


def _nodal_masses(dim, h):
    for i in range(2, NEL + 2):
        m = RHO_L * h * (0.5 if i == NEL + 1 else 1.0)
        ops.mass(i, *([m, m, 0.0] if dim == 2 else [m, m, m, 0.0, 0.0, 0.0]))


def _tip_load(dim):
    ops.load(NEL + 1, *([0.0, 1.0, 0.0] if dim == 2 else [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]))


def _history(dim, elem, mass, damping, load):
    h = _mesh(dim)
    _elements(dim, elem, mass)
    if mass == "nodal":
        _nodal_masses(dim, h)
    ops.timeSeries("Constant", 1)
    if load == "step":
        ops.pattern("Plain", 1, 1)
        _tip_load(dim)
    else:
        ops.pattern("UniformExcitation", 1, 2, "-accel", 1)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 50)
    ops.analysis("Transient")
    w1 = math.sqrt(ops.eigen("-fullGenLapack", 1)[0])
    ops.rayleigh(*_rayleigh_for(damping, w1, XI))
    nstep = 120
    dt = 1.5 * (2.0 * math.pi / w1) / nstep
    out = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0, f"{elem} transient step failed"
        out.append(ops.nodeDisp(NEL + 1, 2))
    return out


@pytest.mark.parametrize("dim", [2, 3])
def test_static_equals_elastic_beam(dim):
    """Premise: elastic section + Legendre IPs => the elastic beam."""
    tips = {}
    for elem in ("ldbc", "ebc"):
        _mesh(dim)
        _elements(dim, elem, "none")
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        _tip_load(dim)
        ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
        assert ops.analyze(1) == 0
        tips[elem] = ops.nodeDisp(NEL + 1, 2)
    exact = L_BEAM ** 3 / (3.0 * E_MOD * I_SEC)
    assert math.isclose(tips["ebc"], exact, rel_tol=1e-9)
    assert math.isclose(tips["ldbc"], tips["ebc"], rel_tol=1e-9), tips


LEGS = ([(m, ld) for m in ("lumped", "nodal") for ld in ("step", "ground")]
        + [("consistent", "step")])


@pytest.mark.parametrize("mass,load", LEGS)
@pytest.mark.parametrize("damping", ["alphaM", "betaK", "betaK0"])
@pytest.mark.parametrize("dim", [2, 3])
def test_dynamic_rayleigh_matches_elastic_beam(dim, damping, mass, load):
    ldbc = _history(dim, "ldbc", mass, damping, load)
    oracle_mass = "consistent" if mass == "consistent" else "nodal"
    ebc = _history(dim, "ebc", oracle_mass, damping, load)
    scale = max(abs(u) for u in ebc)
    assert scale > 0.0, "reference beam did not move -- rig is not driving the tip"
    worst = max(abs(a - b) for a, b in zip(ldbc, ebc))
    assert worst <= 1e-7 * scale, (
        f"LadrunoDispBeamColumn {dim}D [{damping}, {mass} mass, {load}] departs from "
        f"elasticBeamColumn by {worst / scale:.3e} of peak -- inertia, -Q or the "
        "Rayleigh force is being dropped or mis-signed")


def _k_elastic_2d(h):
    """Euler-Bernoulli element stiffness, local = global (beam along +x)."""
    a, b = E_MOD * A_SEC / h, E_MOD * I_SEC
    k = [[0.0] * 6 for _ in range(6)]
    k[0][0] = k[3][3] = a
    k[0][3] = k[3][0] = -a
    for (i, j, v) in ((1, 1, 12 / h ** 3), (1, 2, 6 / h ** 2), (1, 4, -12 / h ** 3), (1, 5, 6 / h ** 2),
                      (2, 2, 4 / h), (2, 4, -6 / h ** 2), (2, 5, 2 / h),
                      (4, 4, 12 / h ** 3), (4, 5, -6 / h ** 2), (5, 5, 4 / h)):
        k[i][j] = k[j][i] = b * v
    return k


@pytest.mark.parametrize("betaK", [0.0, 1.0e-3])
def test_damping_force_response_2d(betaK):
    """getResponse id 12 ('dampingForces', converted to its own buffer by WP-116)
    must be exactly betaK * K_e * v_e for an elastic element under betaK-only
    Rayleigh -- closed form from the analytical stiffness and the node velocities."""
    h = _mesh(2)
    _elements(2, "ldbc", "lumped")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    _tip_load(2)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25); ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 50); ops.analysis("Transient")
    ops.rayleigh(0.0, betaK, 0.0, 0.0)
    for _ in range(5):
        assert ops.analyze(1, 0.001) == 0
    k = _k_elastic_2d(h)
    for e in range(1, NEL + 1):
        v = list(ops.nodeVel(e)) + list(ops.nodeVel(e + 1))
        expect = [betaK * sum(k[i][j] * v[j] for j in range(6)) for i in range(6)]
        got = ops.eleResponse(e, "dampingForces")
        assert len(got) == 6, got
        scale = max(1e-30, max(abs(x) for x in expect))
        assert max(abs(g - x) for g, x in zip(got, expect)) <= 1e-9 * scale + 1e-300, (e, got, expect)
        if betaK == 0.0:
            assert all(x == 0.0 for x in got), got
