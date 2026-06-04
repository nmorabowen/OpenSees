"""LadrunoIndirectControl (Ladruno indirect / CMOD displacement-control static
integrator, classTag 33006) — Zone-A battery.

Generalizes stock DisplacementControl (a single nodal DOF) to a weighted
multi-DOF control quantity  zeta = c.U = sum_k coef_k * U(node_k, dof_k)  with
the per-step constraint  c.dU = Delta-zeta.  A relative control quantity (e.g.
CMOD = u_A - u_B) stays monotone through snap-back where the load and individual
DOFs reverse — the de Borst indirect-control answer to snap-back. ADR-20 §8 #3.

  IC-1  single-DOF, coef=1  ==  stock DisplacementControl, step-for-step
        (proves the predictor/corrector reduces exactly).
  IC-2  the constraint is honored exactly: c.dU == Delta-zeta each converged
        step, including a weighted coefficient (coef=2 => dU = Delta-zeta/2).
  IC-3  a relative two-DOF control quantity (CMOD-style) increments by the
        prescribed amount each step.
  IC-S  parser smoke: documented argument forms are accepted.

Canonical model: the shallow von Mises truss (corotational, snap-through limit
load ~3.8) shared with the LadrunoArcLength battery.

Theory / design: Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

LIC = "LadrunoIndirectControl"

# von Mises truss parameters (snap-through limit load ~3.80)
_H = 0.10
_SPAN = 1.0
_E = 1.0e4
_A = 1.0
_APEX = 2


def _build_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -_SPAN, 0.0)
    ops.node(_APEX, 0.0, _H)
    ops.node(3, _SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(_APEX, 1, 0)            # apex: X fixed by symmetry, Y free
    ops.uniaxialMaterial("Elastic", 1, _E)
    ops.element("corotTruss", 1, 1, _APEX, _A, 1)
    ops.element("corotTruss", 2, _APEX, 3, _A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)


def _setup_solver():
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")


def _trace(integrator_cmd, nsteps):
    """Run a static analysis; return the list of (lambda, uy_apex) per
    converged step."""
    _build_truss()
    _setup_solver()
    ops.integrator(*integrator_cmd)
    ops.analysis("Static")
    path = []
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        path.append((ops.getTime(), ops.nodeDisp(_APEX, 2)))
    return path


# --------------------------------------------------------------------------
# IC-1 : single-DOF indirect control reduces to stock DisplacementControl
# --------------------------------------------------------------------------
def test_ic1_reduces_to_displacement_control():
    incr, nsteps = -0.004, 40
    ref = _trace(("DisplacementControl", _APEX, 2, incr), nsteps)
    lic = _trace((LIC, incr, "-dof", _APEX, 2, 1.0), nsteps)

    assert len(ref) > 20, "reference DisplacementControl produced too few steps"
    assert len(lic) == len(ref), (
        f"step count differs: DisplacementControl {len(ref)} vs {LIC} {len(lic)}"
    )
    for k, ((lr, ur), (ll, ul)) in enumerate(zip(ref, lic)):
        assert abs(lr - ll) < 1.0e-8, f"step {k}: lambda {lr} vs {ll}"
        assert abs(ur - ul) < 1.0e-8, f"step {k}: uy {ur} vs {ul}"


# --------------------------------------------------------------------------
# IC-2 : the constraint c.dU = Delta-zeta is honored exactly (incl. weight)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("coef", [1.0, 2.0, 0.5])
def test_ic2_constraint_increment_exact(coef):
    incr, nsteps = -0.004, 15
    _build_truss()
    _setup_solver()
    ops.integrator(LIC, incr, "-dof", _APEX, 2, coef)
    ops.analysis("Static")
    prev = ops.nodeDisp(_APEX, 2)
    nconv = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        now = ops.nodeDisp(_APEX, 2)
        # c.dU = coef * dU(apex,2) must equal the prescribed control increment
        assert abs(coef * (now - prev) - incr) < 1.0e-9, (
            f"control increment not honored: coef*dU={coef * (now - prev)} vs {incr}"
        )
        prev = now
        nconv += 1
    assert nconv >= 10, f"too few converged steps ({nconv})"


# --------------------------------------------------------------------------
# IC-3 : a relative (CMOD-style) two-DOF control quantity is honored
# --------------------------------------------------------------------------
def test_ic3_relative_cmod_control():
    # Control the RELATIVE quantity  c.U = 1*u(apex,2) + (-1)*u(apex,2) ... needs
    # two distinct free DOFs. Use a 2-DOF cantilever-ish truss pair so the control
    # quantity mixes two independent displacements.
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 2.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)               # node 2: X free, Y fixed
    ops.fix(3, 0, 1)               # node 3: X free, Y fixed
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.element("truss", 2, 2, 3, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, 0.0)
    _setup_solver()

    # CMOD-style relative control: zeta = u(3,1) - u(2,1) (bar-2 elongation)
    incr, nsteps = 0.002, 10
    ops.integrator(LIC, incr, "-dof", 3, 1, 1.0, "-dof", 2, 1, -1.0)
    ops.analysis("Static")
    prev = ops.nodeDisp(3, 1) - ops.nodeDisp(2, 1)
    nconv = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        now = ops.nodeDisp(3, 1) - ops.nodeDisp(2, 1)
        assert abs((now - prev) - incr) < 1.0e-9, (
            f"relative control increment not honored: dZeta={now - prev} vs {incr}"
        )
        prev = now
        nconv += 1
    assert nconv >= 8, f"too few converged steps ({nconv})"


# --------------------------------------------------------------------------
# IC-S : parser smoke — documented argument forms are accepted
# --------------------------------------------------------------------------
def test_ics_parser_smoke():
    # multi -dof + -iter adaptation form
    _build_truss()
    _setup_solver()
    ops.integrator(LIC, -0.004, "-dof", _APEX, 2, 1.0, "-iter", 4, -0.02, -0.001)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "indirect control with -iter failed to take a step"
