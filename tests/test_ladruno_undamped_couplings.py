"""The four pure-penalty coupling elements carry NO Rayleigh damping (WP-123).

LadrunoDistributingCoupling (RBE3), LadrunoKinematicCoupling (RBE2), LadrunoEmbeddedNode
and LadrunoEmbeddedRebar refuse Rayleigh factors (a betaK must not shrink their explicit
dt_cr, ADR 28 §5 / ADR 20 §10.6) and report a zero damping matrix. Since WP-123 they get
that behaviour from one base, `LadrunoUndampedElement`; before, each carried its own copy,
and the implicit-transient crash of a missing getDamp override was fixed in them one after
another (#219 RBE3, #220 embedded; LEDGER_quirks "A no-op `setRayleighDampingFactors`
WITHOUT a `getDamp` override").

The existing per-element smoke tests only assert that Newmark does not crash. These assert
the contract itself, so a change to the shared code that lets damping leak in fails here:

  1. a damped transient (global rayleigh with all four factors, Newmark and HHT) is
     BIT-IDENTICAL to the same run with rayleigh(0,0,0,0). The models carry no nodal
     masses — all inertia is the bipenalty element mass — so the only way the factors can
     act is through the element;
  2. eleResponse 'dampingForce' is exactly zero in a moving state;
  3. the explicit critical step (element self-report, and the integrator's value under
     CentralDifferenceLadruno) does not depend on betaK;
  4. the element's damping is ZERO in absolute terms, not merely factor-independent:
     under Newmark average acceleration (gamma=1/2, beta=1/4) a linear undamped system
     conserves discrete energy exactly, so a constant load from rest gives a free
     vibration whose amplitude must not decay. Test 1 alone cannot see a getDamp that
     is nonzero regardless of the factors (both runs would carry it equally) -- WP-123
     mutation row C proved that gap before this case was added. The case uses
     `algorithm Linear`: these elements put no D·v in their residual, so a spurious C
     only reaches the tangent, which Newton would iterate away.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DT = 1.0e-4
RAYLEIGH_ON = (0.5, 1.0e-3, 1.0e-3, 1.0e-3)   # alphaM, betaK, betaK0, betaKc
RAYLEIGH_OFF = (0.0, 0.0, 0.0, 0.0)


# ------------------------------------------------------------------ models
def _face(kind):
    """RBE3 / RBE2: 6-DOF reference node 1 tied to 4 fixed corner nodes; the reference
    mass is the bipenalty m_p / I_p. Returns the free (node, dof) list."""
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    for i, c in enumerate([(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)]):
        ops.node(2 + i, c[0], c[1], 0.0)
        ops.fix(2 + i, 1, 1, 1)
    ops.element(kind, 1, 1, 4, 2, 3, 4, 5, "-k", 1.0e7, "-bipenalty", "-dtcr", 1.0e-3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 100.0, -50.0, 25.0, 10.0, -20.0, 1000.0)
    return [(1, d) for d in range(1, 7)]


def _embedded_node():
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)                              # constrained node, carries m_p
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1)                                     # host fixed => g = u_c
    ops.element("LadrunoEmbeddedNode", 1, 1, 1, 2, "-shape", 1.0, "-k", 1.0e6,
                "-bipenalty", "-dtcr", 1.0e-3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, -0.5, 0.25)
    return [(1, d) for d in range(1, 4)]


def _embedded_rebar():
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)                              # rebar node, carries m_p
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1)                                     # host fixed
    ops.element("LadrunoEmbeddedRebar", 1, 1, 1, 2, "-shape", 1.0, "-dir", 1.0, 0.0, 0.0,
                "-perfect", 1.0e4, "-kt", 1.0e4, "-bipenalty", "-dtcr", 1.0e-3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, -0.5, 0.25)
    return [(1, d) for d in range(1, 4)]


MODELS = {
    "LadrunoDistributingCoupling": lambda: _face("LadrunoDistributingCoupling"),
    "LadrunoKinematicCoupling": lambda: _face("LadrunoKinematicCoupling"),
    "LadrunoEmbeddedNode": _embedded_node,
    "LadrunoEmbeddedRebar": _embedded_rebar,
}
ELEMENTS = sorted(MODELS)
INTEGRATORS = {
    "Newmark": ("Newmark", 0.5, 0.25),
    "HHT": ("HHT", 0.9),
}


def _implicit_run(name, rayleigh, integrator, nsteps=20):
    """Build, apply `rayleigh`, run `nsteps` implicit steps; return the recorded
    (disp, vel) of every free DOF at every step, plus the final dampingForce."""
    ops.wipe()
    free = MODELS[name]()
    ops.rayleigh(*rayleigh)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 30)
    ops.algorithm("Newton")
    ops.integrator(*INTEGRATORS[integrator])
    ops.analysis("Transient")
    hist = []
    for _ in range(nsteps):
        assert ops.analyze(1, DT) == 0
        hist.append(tuple((ops.nodeDisp(n, d), ops.nodeVel(n, d)) for n, d in free))
    return hist, list(ops.eleResponse(1, "dampingForce"))


def _explicit_dtcr(name, rayleigh):
    ops.wipe()
    MODELS[name]()
    ops.rayleigh(*rayleigh)
    self_report = ops.eleResponse(1, "dtcr")[0]
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-5) == 0
    return self_report, ops.criticalTimeStep()


# ------------------------------------------------------------------ 1. no damping leaks in
@pytest.mark.parametrize("integrator", sorted(INTEGRATORS))
@pytest.mark.parametrize("name", ELEMENTS)
def test_damped_transient_is_bit_identical_to_undamped(name, integrator):
    on, _ = _implicit_run(name, RAYLEIGH_ON, integrator)
    off, _ = _implicit_run(name, RAYLEIGH_OFF, integrator)
    moving = max(abs(v) for step in off for _, v in step)
    assert moving > 0.0, "the probe must actually move, or equality proves nothing"
    assert on == off          # exact: the element must not feel alphaM/betaK/betaK0/betaKc


# ------------------------------------------------------------------ 2. zero damping force
@pytest.mark.parametrize("name", ELEMENTS)
def test_damping_force_is_zero_while_moving(name):
    hist, damp = _implicit_run(name, RAYLEIGH_ON, "Newmark")
    assert any(v != 0.0 for _, v in hist[-1]), "velocity must be nonzero at the probe"
    assert len(damp) > 0, "dampingForce must be answered (an empty reply would pass vacuously)"
    assert all(math.isfinite(x) and x == 0.0 for x in damp)


# ------------------------------------------------------------------ 4. no absolute damping
@pytest.mark.parametrize("name", ELEMENTS)
def test_free_vibration_does_not_decay(name):
    """Constant load from rest, Newmark average acceleration, no Rayleigh: the loaded DOFs
    oscillate about the static solution forever. Compare the oscillation amplitude in the
    first and last fifth of the run on every free DOF that moves."""
    ops.wipe()
    free = MODELS[name]()
    ops.remove("loadPattern", 1)                      # the model's Linear pattern -> Constant
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    if name in ("LadrunoDistributingCoupling", "LadrunoKinematicCoupling"):
        ops.load(1, 100.0, -50.0, 25.0, 10.0, -20.0, 1000.0)
    else:
        ops.load(1, 1.0, -0.5, 0.25)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 1)
    # Linear, not Newton: these elements add no D·v to their residual, so a spurious
    # getDamp only pollutes the TANGENT and Newton would iterate it away (same answer,
    # slower). One solve per step with the element's own tangent is exact for these
    # linear ties; a C-polluted tangent then shows up as a wrong (decaying) response.
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    nsteps = 1500                                     # >= 30 periods of the slowest DOF here
    hist = []
    for _ in range(nsteps):
        assert ops.analyze(1, DT) == 0
        hist.append([ops.nodeDisp(n, d) for n, d in free])
    k = nsteps // 5
    checked = 0
    for j in range(len(free)):
        series = [h[j] for h in hist]
        early = max(series[:k]) - min(series[:k])
        late = max(series[-k:]) - min(series[-k:])
        if early < 1.0e-12:
            continue                                  # a DOF the load does not excite
        checked += 1
        assert late > 0.95 * early, f"DOF {free[j]}: amplitude {early:.3e} -> {late:.3e} (damped)"
    assert checked > 0, "no DOF moved: the probe proves nothing"


# ------------------------------------------------------------------ 3. dt_cr ignores betaK
@pytest.mark.parametrize("name", ELEMENTS)
def test_explicit_critical_step_ignores_rayleigh(name):
    self_on, glob_on = _explicit_dtcr(name, (0.0, 1.0e-3, 0.0, 0.0))
    self_off, glob_off = _explicit_dtcr(name, RAYLEIGH_OFF)
    assert self_on > 0.0
    assert self_on == self_off
    assert glob_on == glob_off
