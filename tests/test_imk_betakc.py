"""betaKc (committed-stiffness) Rayleigh damping on LadrunoIMKBeam 2D/3D -- WP-118.

LadrunoIMKBeam(2d)::commitState() never chained to Element::commitState(), so
Kc -- the committed stiffness behind betaKc -- was never refreshed. It was
captured once, correctly, when `rayleigh` ran (Domain::addElement calls
update(), so the tangent is valid by then) and then stayed frozen at the
element's INITIAL stiffness: betaKc silently behaved like betaK0. Invisible in
a linear elastic model; wrong as soon as the tangent changes (hinge yielding,
large rotation).

Gate: differential against elasticBeamColumn, which chains Element::commitState.
With no hinge materials IMK is the elastic beam, so under betaKc-only Rayleigh
every displacement history must agree. Legs:
  Corotational transf  -- large rotation changes the global tangent, so Kc must
                          be REFRESHED every commit. THE discriminating leg:
                          pre-fix it departs by ~45% of peak.
  Linear transf        -- constant tangent; pins that Kc is captured correctly
                          (a regression guard, not diagnostic of the refresh).
x {2D, 3D} x {element -mass, nodal mass}. The oracle carries nodal masses
(vanilla ElasticBeam2d double-counts ground-motion Q with element -mass; step
loads here, but keep the oracle uniform).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

L_BEAM, A_SEC, E_MOD, I_SEC = 3.0, 0.02, 2.0e5, 1.5e-4
G_MOD, J_TOR = 8.0e4, 3.0e-4
RHO_L = 0.8
NEL = 4
XI = 0.05
# first cantilever mode (continuous): w1 = 1.875^2 sqrt(EI / (m L^4)); used only
# to size betaKc and the time step, so no eigen call is needed
W1 = 1.875104 ** 2 * math.sqrt(E_MOD * I_SEC / (RHO_L * L_BEAM ** 4))
BETA_KC = 2.0 * XI / W1


def _history(dim, elem, mass, transf, tip_load, beta_kc=BETA_KC):
    ops.wipe()
    if dim == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf(transf, 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)
    h = L_BEAM / NEL
    for i in range(NEL + 1):
        ops.node(i + 1, *((i * h, 0.0) if dim == 2 else (i * h, 0.0, 0.0)))
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))
    for e in range(1, NEL + 1):
        if dim == 2:
            sec = (A_SEC, E_MOD, I_SEC)
        else:
            sec = (A_SEC, E_MOD, G_MOD, J_TOR, I_SEC, I_SEC)
        name = "LadrunoIMKBeam" if elem == "imk" else "elasticBeamColumn"
        args = [name, e, e, e + 1, *sec, 1]
        if mass == "element":
            args += ["-mass", RHO_L]
        ops.element(*args)
    if mass == "nodal":
        for i in range(2, NEL + 2):
            m = RHO_L * h * (0.5 if i == NEL + 1 else 1.0)
            ops.mass(i, *([m, m, 0.0] if dim == 2 else [m, m, m, 0.0, 0.0, 0.0]))

    # betaKc-only Rayleigh, in the usual script position (before the analysis)
    ops.rayleigh(0.0, 0.0, 0.0, beta_kc)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(NEL + 1, *([0.0, tip_load, 0.0] if dim == 2 else [0.0, tip_load, 0.0, 0.0, 0.0, 0.0]))
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 50)
    ops.analysis("Transient")
    nstep = 150
    dt = 2.0 * (2.0 * math.pi / W1) / nstep
    out = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0, f"{elem} transient step failed"
        out.append(ops.nodeDisp(NEL + 1, 2))
    return out


# small load (linear regime) for Linear; large load (tip rotation ~0.5 rad) for Corotational
LOAD = {"Linear": 1.0, "Corotational": 0.35 * E_MOD * I_SEC / L_BEAM ** 2 * 2.0}


@pytest.mark.parametrize("mass", ["element", "nodal"])
@pytest.mark.parametrize("transf", ["Linear", "Corotational"])
@pytest.mark.parametrize("dim", [2, 3])
def test_betakc_matches_elastic_beam(dim, transf, mass):
    imk = _history(dim, "imk", mass, transf, LOAD[transf])
    ebc = _history(dim, "ebc", "nodal", transf, LOAD[transf])
    scale = max(abs(u) for u in ebc)
    assert scale > 0.0
    worst = max(abs(a - b) for a, b in zip(imk, ebc))
    assert worst <= 1e-7 * scale, (
        f"IMK {dim}D [{transf}, {mass} mass] under betaKc departs from elasticBeamColumn by "
        f"{worst / scale:.3e} of peak -- Kc not refreshed at commit (commitState not "
        "chaining to Element::commitState) or captured wrong")


@pytest.mark.parametrize("dim", [2, 3])
def test_corotational_rig_actually_rotates(dim):
    """Premise of the Corotational legs: the tip rotates enough that a Kc frozen
    at the undeformed configuration differs from the refreshed one."""
    _history(dim, "ebc", "nodal", "Corotational", LOAD["Corotational"])
    rot_dof = 3 if dim == 2 else 6
    rot = max(abs(ops.nodeDisp(NEL + 1, rot_dof)), 0.0)
    assert rot > 0.1, f"tip rotation {rot:.3f} rad -- too small to exercise the Kc refresh"


@pytest.mark.parametrize("dim", [2, 3])
def test_betakc_actually_damps(dim):
    """Premise: the chosen betaKc visibly damps the response (else the gate is vacuous)."""
    damped = _history(dim, "ebc", "nodal", "Linear", 1.0)
    undamped = _history(dim, "ebc", "nodal", "Linear", 1.0, beta_kc=0.0)
    late = slice(len(damped) // 2, None)
    amp_d = max(damped[late]) - min(damped[late])
    amp_u = max(undamped[late]) - min(undamped[late])
    assert amp_d < 0.8 * amp_u, (amp_d, amp_u)
