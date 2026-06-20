"""Mass-scaling + HRZ FIDELITY validation (Ladruno ADR 37).

The shipped SMS tests (test_centralDifferenceSMS_integrator.py) prove the feature
runs, stabilizes, and restores. These tests close the ACCURACY gap: that the
mass-scaled response actually MATCHES the unscaled reference (Tier 0), and that
the review-fix behaviors are asserted, not just warned (Tier 1).

Plan: Ladruno_implementation/37_ladruno_mass_scaling_validation_plan.md.

Validation finding (recorded here for provenance): SMS demonstrably raises the
global stable step — on the chain below, unscaled central difference diverges at
dt=0.011 while CentralDifferenceSMS(dtTarget=0.012) is stable through 0.012. The
fidelity tests below seed the FUNDAMENTAL MODE so the comparison reflects the
mass-scaling frequency shift (the quantity of interest), not the per-mode
central-difference period error (~(omega*dt)^2/24) that a point-seed's high modes
would otherwise inject.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

SMS = "CentralDifferenceSMS"
CDL = "CentralDifferenceLadruno"

FREE = (2, 3, 4)   # free nodes of the chain (x-dof)


# --------------------------------------------------------------------------
# numpy-free helpers
# --------------------------------------------------------------------------
def _interp(ts, us, t):
    if t <= ts[0]:
        return us[0]
    if t >= ts[-1]:
        return us[-1]
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t:
            lo = mid
        else:
            hi = mid
    w = (t - ts[lo]) / (ts[hi] - ts[lo])
    return us[lo] * (1.0 - w) + us[hi] * w


def _rms(xs):
    return math.sqrt(sum(v * v for v in xs) / len(xs)) if xs else 0.0


# --------------------------------------------------------------------------
# 3-element chain, throttling element e2 INTERIOR (both nodes free) so scaling
# fully raises its stable step:
#   node1(fix) - e1(bulk) - node2 - e2(tiny) - node3 - e3(bulk) - node4(free)
#   bulk e1/e3: L=1.0, E=1e4, rho=1  -> dt_e=0.0100
#   tiny e2   : L=0.1, E=138.4, rho=1 -> dt_e~0.0085, m(node)=0.05
# --------------------------------------------------------------------------
def _model():
    """Define the chain (no integrator / analysis yet)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1)
    ops.node(3, 1.1, 0.0); ops.fix(3, 0, 1)
    ops.node(4, 2.1, 0.0); ops.fix(4, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.uniaxialMaterial("Elastic", 2, 138.4)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    ops.element("Truss", 2, 2, 3, 1.0, 2, "-rho", 1.0)
    ops.element("Truss", 3, 3, 4, 1.0, 1, "-rho", 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")


def _mode1():
    """Return (omega1, {node: phi_x}) for the fundamental mode, normalized so the
    largest |component| is 1."""
    _model()
    ops.system("FullGeneral")
    lam = ops.eigen("-fullGenLapack", 1)
    w1 = math.sqrt(lam[0])
    phi = {n: ops.nodeEigenvector(n, 1, 1) for n in FREE}
    pmax = max(abs(v) for v in phi.values())
    phi = {n: v / pmax for n, v in phi.items()}
    return w1, phi


def _run_seeded(integrator_args, dt, T, phi, amp):
    """Seed amp*phi (fundamental mode) and free-vibrate to ~T. Return (t[], u4[])."""
    _model()
    for n in FREE:
        ops.setNodeDisp(n, 1, amp * phi[n], "-commit")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    t = [ops.getTime()]
    u = [ops.nodeDisp(4, 1)]
    for _ in range(int(round(T / dt))):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime())
        u.append(ops.nodeDisp(4, 1))
    return t, u


# --------------------------------------------------------------------------
# T-ACC: the scaled run (bulk dt) must MATCH the unscaled reference (fine dt) in
#        the fundamental mode. Seeding mode 1 isolates the mass-scaling shift from
#        high-mode discretization drift. This is the accuracy claim mass scaling
#        needs — usable, not merely bounded.
# --------------------------------------------------------------------------
def test_sms_reference_match():
    DT_E2 = 0.0085           # tiny element stable step (the throttle)
    DT_TARGET = 0.010        # modest target -> small added mass (~few %)
    AMP = 1.0e-3

    w1, phi = _mode1()
    T1 = 2.0 * math.pi / w1
    # 1.5 fundamental periods: long enough to characterize the response through a
    # full cycle, short enough that the (intended, modest) mode-1 frequency shift
    # Δf/f~0.7% accumulates to only ~2π·1.5·ε/√3 ~ 4% phase-drift RMS — i.e. the
    # tolerance is set by the physics of the shift, not chosen to pass.
    T = 1.5 * T1

    t_ref, u_ref = _run_seeded((CDL,), 0.9 * DT_E2, T, phi, AMP)
    t_sms, u_sms = _run_seeded((SMS, DT_TARGET), 0.9 * DT_TARGET, T, phi, AMP)

    assert all(math.isfinite(x) for x in u_sms) and _rms(u_sms) > 1e-6, "SMS run degenerate"
    assert all(math.isfinite(x) for x in u_ref) and _rms(u_ref) > 1e-6, "ref run degenerate"
    assert len(t_sms) < len(t_ref), "SMS did not take a larger step than the reference"

    t_end = min(t_ref[-1], t_sms[-1])
    diff, ref_on_grid = [], []
    for ts, us in zip(t_sms, u_sms):
        if ts > t_end:
            break
        ur = _interp(t_ref, u_ref, ts)
        diff.append(us - ur)
        ref_on_grid.append(ur)
    rel = _rms(diff) / _rms(ref_on_grid)
    assert rel < 0.05, (
        "SMS@bulk-dt does not match unscaled@fine-dt in mode 1: relative RMS "
        "%.3f%% (steps sms=%d vs ref=%d)" % (100 * rel, len(t_sms), len(t_ref))
    )


# --------------------------------------------------------------------------
# T-MODAL: SMS adds mass -> shifts frequencies. Assert the fundamental shift is
#          small at a modest dtTarget, and GROWS as dtTarget is pushed (the honest
#          fidelity tradeoff curve). Uses ops.eigen before/after the SMS injection.
# --------------------------------------------------------------------------
def _f1_after_sms(dtTarget):
    """Fundamental frequency after one SMS domainChanged injection. Returns f1.
    Inject by priming a transient step with SMS, THEN run eigen while the scaled
    nodal mass is committed (the SMS integrator is still alive, so the mass stands)."""
    _model()
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(SMS, dtTarget)
    ops.analysis("Transient")
    ops.analyze(1, 0.5 * dtTarget)          # triggers domainChanged -> mass injection
    ops.system("FullGeneral")               # eigen needs a full system
    lam = ops.eigen("-fullGenLapack", 1)
    return math.sqrt(lam[0]) / (2.0 * math.pi)


def test_sms_frequency_shift_grows():
    w1, _ = _mode1()
    f1_0 = w1 / (2.0 * math.pi)             # unscaled fundamental

    f1_modest = _f1_after_sms(0.0095)       # just above the tiny element -> small s
    f1_aggr = _f1_after_sms(0.020)          # well above -> large s

    d_modest = abs(f1_modest - f1_0) / f1_0
    d_aggr = abs(f1_aggr - f1_0) / f1_0

    assert math.isfinite(f1_modest) and math.isfinite(f1_aggr), (f1_modest, f1_aggr)
    # modest scaling barely moves f1
    assert d_modest < 0.01, "modest SMS shifted f1 by %.3f%% (expected <1%%)" % (100 * d_modest)
    # pushing dtTarget moves it MORE (the honest tradeoff: added mass lowers f1)
    assert d_aggr > d_modest, (
        "aggressive SMS (Δf1=%.3f%%) did not shift f1 more than modest (%.3f%%)"
        % (100 * d_aggr, 100 * d_modest)
    )
    assert f1_aggr < f1_0, "added mass should LOWER the fundamental frequency"


# --------------------------------------------------------------------------
# T-HRZ-END2END: on a CONSISTENT-mass beam (rotational DOFs), -lump rowsum is
#   DANGEROUS — it zeros the rotational mass, the dt_cr eigensolve discards those
#   massless modes, and the reported dt_cr is grossly INFLATED (unsafe). A run at
#   rowsum's reported step diverges; hrz gives a physical, safe estimate.
# --------------------------------------------------------------------------
def _beam(lump):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0); ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 200.0, 80.0, 1.0, 1.0, 1.0, 1,
                "-mass", 1.0, "-cMass")             # consistent mass: rotational DOFs
    ops.setNodeDisp(2, 2, 1.0e-3, "-commit")        # transverse tip displacement
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", "-lump", lump)
    ops.analysis("Transient")


def _beam_dtcr(lump):
    _beam(lump)
    ops.analyze(1, 1.0e-6)
    return ops.criticalTimeStep()


def _beam_maxtip(lump, dt, nsteps):
    _beam(lump)
    m = abs(ops.nodeDisp(2, 2))
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return float("inf")
        u = abs(ops.nodeDisp(2, 2))
        if not math.isfinite(u) or u > 1.0e3:
            return float("inf")
        m = max(m, u)
    return m


def test_hrz_safe_dtcr_vs_rowsum_on_beam():
    dt_rowsum = _beam_dtcr("rowsum")
    dt_hrz = _beam_dtcr("hrz")

    # rowsum INFLATES dt_cr on a rotational-DOF element (zeros rotational mass ->
    # the eigensolve drops those modes). hrz keeps it physical.
    assert dt_rowsum > 2.0 * dt_hrz, (
        "expected rowsum dt_cr grossly inflated vs hrz on a consistent beam: "
        "rowsum=%.5g hrz=%.5g" % (dt_rowsum, dt_hrz)
    )

    # a run at rowsum's reported step DIVERGES (its estimate is unsafe) ...
    u_unsafe = _beam_maxtip("hrz", 0.9 * dt_rowsum, 200)
    assert not math.isfinite(u_unsafe) or u_unsafe > 1.0e2, (
        "running at rowsum's reported dt_cr (%.5g) should be unstable, max|tip|=%r"
        % (0.9 * dt_rowsum, u_unsafe)
    )
    # ... while a step safely below the hrz/diagonal estimate is stable.
    u_safe = _beam_maxtip("hrz", 0.4 * dt_hrz, 200)
    assert math.isfinite(u_safe) and u_safe < 1.0, (
        "a step below the hrz dt_cr should be stable, got max|tip|=%r" % u_safe
    )
