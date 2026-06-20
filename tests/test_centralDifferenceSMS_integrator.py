"""CentralDifferenceSMS (Ladruno selective mass scaling, classTag 33007) —
Zone-A validation.

CentralDifferenceSMS subclasses CentralDifferenceLadruno and, at domainChanged,
inflates the mass of elements whose per-element stable step dt_e < dtTarget so the
whole model can run at a global step set by the bulk. The central claim (test
SMS-2): a mesh with one tiny stiff element is UNSTABLE under plain explicit CD at
the bulk timestep, but STABLE under CentralDifferenceSMS targeting that timestep.

Theory: Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

SMS = "CentralDifferenceSMS"
CDL = "CentralDifferenceLadruno"


def _chain_maxdisp(integrator_args, dt, nsteps, d0=1.0e-3):
    """Fixed–free 2-truss chain in 2D: element 1 is short+stiff+light (tiny dt_e),
    element 2 is the soft bulk (~100x larger dt_e). Seed node 3 with d0 and free-
    vibrate. Returns max |u_x(node 3)| over the run (inf if it blew up / diverged)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)
    ops.node(3, 1.1, 0.0); ops.fix(3, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e6)   # stiff short element -> tiny dt_e1
    ops.uniaxialMaterial("Elastic", 2, 1.0e4)   # soft bulk element  -> larger dt_e2
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    ops.element("Truss", 2, 2, 3, 1.0, 2, "-rho", 1.0)
    ops.setNodeDisp(3, 1, d0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    umax = abs(ops.nodeDisp(3, 1))
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return float("inf")
        u = abs(ops.nodeDisp(3, 1))
        if not math.isfinite(u):
            return float("inf")
        if u > umax:
            umax = u
    return umax


# dt_e1 ~ 1e-4 (stiff/light), dt_e2 ~ 1e-2 (soft bulk). Target the bulk: well above
# the tiny element's limit, so plain CD is unstable but SMS scales e1 to recover it.
DT_TARGET = 5.0e-3
N = 600


# --------------------------------------------------------------------------
# SMS-1: the integrator constructs and takes a step (parse + lifecycle smoke).
# --------------------------------------------------------------------------
def test_sms_constructs_and_runs():
    umax = _chain_maxdisp((SMS, DT_TARGET, "-verbose"), DT_TARGET, 5)
    assert math.isfinite(umax), "CentralDifferenceSMS failed to run a few steps"


# --------------------------------------------------------------------------
# SMS-2 (the central claim): at the bulk timestep, plain CD DIVERGES (the tiny
# element's dt_cr is ~50x smaller) while CentralDifferenceSMS stays STABLE.
# --------------------------------------------------------------------------
def test_sms_stabilizes_tiny_element():
    u_sms = _chain_maxdisp((SMS, DT_TARGET), DT_TARGET, N)
    u_cdl = _chain_maxdisp((CDL,), DT_TARGET, N)

    # SMS: bounded free vibration seeded at 1e-3 (added inertia shifts frequencies
    # but the response stays O(seed), not exploding).
    assert math.isfinite(u_sms) and u_sms < 1.0, (
        "SMS run not stable/bounded: max|u|=%r" % u_sms
    )
    # plain CD at the same dt is unstable on the tiny element -> blows up.
    assert (not math.isfinite(u_cdl)) or u_cdl > 1.0e2, (
        "plain CentralDifferenceLadruno was expected to DIVERGE at dt=%g "
        "(tiny-element dt_cr ~50x smaller), but max|u|=%r" % (DT_TARGET, u_cdl)
    )
    # and concretely: SMS is many orders of magnitude smaller than unscaled CD.
    assert u_sms < u_cdl, (u_sms, u_cdl)


# --------------------------------------------------------------------------
# SMS-3: dtTarget below the smallest element's dt_e is a no-op (nothing to scale)
# -> behaves like plain CDL (also stable here since dt is tiny).
# --------------------------------------------------------------------------
def test_sms_noop_when_target_below_min():
    # 5e-5 is below dt_e1 (~1e-4): no element exceeds target, nothing scaled.
    u = _chain_maxdisp((SMS, 5.0e-5), 5.0e-5, 200)
    assert math.isfinite(u) and u < 1.0, "SMS no-op run not stable: %r" % u


# --------------------------------------------------------------------------
# SMS-4 (review SMS-NO-RESTORE, critical): the injected fictitious mass must be
# RESTORED when the SMS integrator is torn down, so a later stage / integrator on
# the same Domain sees the ORIGINAL node masses. nodeMass() returns the node's own
# (mutated) mass matrix, so we can assert the round-trip directly.
# --------------------------------------------------------------------------
def _setup_chain():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)
    ops.node(3, 1.1, 0.0); ops.fix(3, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e6)
    ops.uniaxialMaterial("Elastic", 2, 1.0e4)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    ops.element("Truss", 2, 2, 3, 1.0, 2, "-rho", 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")


def test_sms_restores_mass_on_teardown():
    _setup_chain()
    base = ops.nodeMass(2)  # node's own mass (element mass is NOT nodal) -> ~0
    assert abs(base[0]) < 1e-9, "expected zero baseline nodeMass, got %r" % base

    ops.integrator(SMS, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    scaled = ops.nodeMass(2)
    assert scaled[0] > 1.0, "SMS did not inject nodal mass: %r" % scaled

    # tear down the SMS integrator by replacing it -> its destructor must restore.
    ops.wipeAnalysis()
    restored = ops.nodeMass(2)
    assert abs(restored[0]) < 1e-6, (
        "SMS did not restore node mass on teardown: baseline~0, after=%r "
        "(scaled was %r) -- the Domain is left corrupted (SMS-NO-RESTORE)"
        % (restored, scaled)
    )
