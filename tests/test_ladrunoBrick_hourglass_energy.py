"""LadrunoBrick — hourglass / stabilization energy reporting.

The element exposes a scalar `"hourglassEnergy"` response: the recoverable elastic
hourglass energy (uri 'stiffness': ½κ·Σq²) or condensed-stabilization energy (ssp:
½·u·Kstab·u). It is the GLSTAT-style spurious-mode diagnostic LS-DYNA reports in
MATSUM/GLSTAT — letting an explicit run check that hourglass energy stays a small
fraction of the internal energy.

Gates:
  * constant strain / rigid body  ->  Ehg = 0 (γ is orthogonal to the linear and
    rigid fields), for uri-stiffness AND ssp;
  * a genuinely hourglassing (bending) state -> Ehg > 0;
  * std / bbar / physical -> Ehg = 0 (no stored hourglass energy);
  * uri-viscous stores nothing either — in a static solve (no velocity) it reports
    0; under dynamics it instead reports its CUMULATIVE dissipation (covered in
    test_ladrunoBrick_viscous_dissipation.py);
  * Ehg >= 0 always, and on a single element Ehg <= the external work (it is a
    part of the strain energy, so it cannot exceed it).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


def _hg(tag=1):
    ops.eleResponse(tag, "forces")          # ensure trial strain/state is set
    return list(ops.eleResponse(tag, "hourglassEnergy"))[0]


# --------------------------------------------------------------------------
# constant strain -> zero hourglass energy (gamma orthogonality)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("extra", [
    ["-formulation", "std"],
    ["-formulation", "bbar"],
    ["-formulation", "uri", "-hourglass", "stiffness"],
    ["-formulation", "uri", "-hourglass", "viscous"],
    ["-formulation", "ssp"],
    ["-formulation", "uri", "-hourglass", "physical"],
])
def test_constant_strain_zero_hourglass_energy(extra):
    """A uniform uniaxial-stress state is constant strain: every formulation must
    report Ehg = 0 (the hourglass modes carry no energy in a linear field)."""
    E, nu, sig0 = 1000.0, 0.3, 2.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    # 1/8-symmetry restraints consistent with the exact uniaxial field
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, sig0 / 4.0, 0.0, 0.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0

    Ehg = _hg()
    assert abs(Ehg) < 1e-9, f"{extra}: constant-strain hourglass energy {Ehg:.3e} should be ~0"


# --------------------------------------------------------------------------
# bending -> positive hourglass/stabilization energy for the controlled forms
# --------------------------------------------------------------------------
def _cantilever_hg(form_args, nx=4, L=10.0, E=1000.0, nu=0.0, P=1.0):
    """Slender cantilever (1 element deep) under a transverse tip load; returns
    (tip_z, Ehg of the root element)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = L / nx

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k

    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * 1.0, k * 1.0)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    e = 1
    for i in range(nx):
        ops.element("LadrunoBrick", e,
                    nt(i, 0, 0), nt(i + 1, 0, 0), nt(i + 1, 1, 0), nt(i, 1, 0),
                    nt(i, 0, 1), nt(i + 1, 0, 1), nt(i + 1, 1, 1), nt(i, 1, 1),
                    1, *form_args)
        e += 1
    for k in range(2):
        for j in range(2):
            ops.fix(nt(0, j, k), 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tip = [nt(nx, j, k) for k in range(2) for j in range(2)]
    for n in tip:
        ops.load(n, 0.0, 0.0, P / len(tip))
    ops.system("FullGeneral"); ops.numberer("RCM"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    tip_z = sum(ops.nodeDisp(n)[2] for n in tip) / len(tip)
    return tip_z, _hg(1)


@pytest.mark.parametrize("form_args", [
    ["-formulation", "uri", "-hourglass", "stiffness"],
    ["-formulation", "ssp"],
])
def test_bending_has_positive_hourglass_energy(form_args):
    """In a bending-dominated state the controlled formulations store positive
    hourglass/stabilization energy (for uri this energy *is* the bending stiffness
    of the 1-pt element)."""
    _, Ehg = _cantilever_hg(form_args)
    assert Ehg > 0.0, f"{form_args}: bending hourglass energy {Ehg:.3e} should be > 0"


@pytest.mark.parametrize("form_args", [
    ["-formulation", "std"],
    ["-formulation", "bbar"],
    ["-formulation", "uri", "-hourglass", "physical"],
    ["-formulation", "uri", "-hourglass", "viscous"],
])
def test_fully_integrated_and_viscous_report_zero(form_args):
    """std/bbar (full integration) and physical (assumed strain folded into the
    strain energy) store no separable hourglass energy; uri-viscous dissipates
    rather than stores and accumulates nothing in this STATIC solve (no velocity).
    All report exactly 0 here."""
    _, Ehg = _cantilever_hg(form_args)
    assert Ehg == 0.0, f"{form_args}: expected 0 hourglass energy, got {Ehg:.3e}"


# --------------------------------------------------------------------------
# physical bound: on a single element Ehg <= external work (part of strain energy)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("form_args", [
    ["-formulation", "uri", "-hourglass", "stiffness"],
    ["-formulation", "ssp"],
])
def test_hourglass_energy_bounded_by_external_work(form_args):
    """For a single element under a static load, the elastic strain energy equals
    the external work ½·F·u. The hourglass/stabilization energy is one part of
    that strain energy, so 0 < Ehg <= ½·F·u (+tol)."""
    E, nu = 1000.0, 0.3
    loads = {5: (10.0, 0.0, 0.0), 6: (0.0, 8.0, 0.0),
             7: (0.0, 0.0, -12.0), 8: (5.0, 5.0, 5.0)}
    nodes = {1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.00), 3: (1.10, 1.00, 0.10),
             4: (0.05, 0.95, 0.00), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
             7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95)}
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in nodes.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in loads.items():
        ops.load(n, fx, fy, fz)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
    assert ops.analyze(1) == 0

    Wext = 0.0
    for n, (fx, fy, fz) in loads.items():
        u = ops.nodeDisp(n)
        Wext += 0.5 * (fx * u[0] + fy * u[1] + fz * u[2])
    Ehg = _hg(1)
    assert Ehg > 0.0, f"{form_args}: Ehg should be > 0 in this mixed-load state"
    assert Ehg <= Wext * (1.0 + 1e-6) + 1e-9, (
        f"{form_args}: Ehg {Ehg:.4e} should not exceed external work {Wext:.4e}"
    )
