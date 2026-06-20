"""Mass-scaling + HRZ FIDELITY validation (Ladruno ADR 37).

These close the ACCURACY gap the shipped SMS tests left: that mass scaling enables a
SUPRA-stable step (where plain central difference diverges) AND keeps the response
faithful. Reworked after a Tier-0 adversarial review found the first cut vacuous —
every test now has a control that FAILS on a broken / no-op / absent SMS.

Plan: Ladruno_implementation/37_ladruno_mass_scaling_validation_plan.md.

Empirical anchor (3-element chain below): plain CentralDifferenceLadruno diverges at
dt=0.011 (above the global stable step ~0.0105) while CentralDifferenceSMS(dtTarget
=0.012) stays stable through 0.012. That supra-stable regime is what makes the
control legs meaningful.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

SMS = "CentralDifferenceSMS"
CDL = "CentralDifferenceLadruno"

FREE = (2, 3, 4)
M_E2_NODE = 0.05          # tiny element e2 lumped nodal mass = 0.5*rho*A*L2 = 0.5*1*0.1
DT_SUPRA = 0.011          # > unscaled global stable step (~0.0105), < SMS-raised limit


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


def _rel_rms(t_a, u_a, t_b, u_b):
    """Relative RMS of (run a) vs the interpolated (run b reference), on a's grid."""
    t_end = min(t_a[-1], t_b[-1])
    diff, base = [], []
    for ts, ua in zip(t_a, u_a):
        if ts > t_end:
            break
        ub = _interp(t_b, u_b, ts)
        diff.append(ua - ub)
        base.append(ub)
    return _rms(diff) / _rms(base) if _rms(base) > 0 else float("inf")


# --------------------------------------------------------------------------
# node1(fix) - e1(bulk) - node2 - e2(tiny,interior) - node3 - e3(bulk) - node4(free)
#   bulk e1/e3: L=1.0, E=1e4, rho=1
#   tiny e2   : L=0.1, E=138.4, rho=1  -> m_node=0.05; its nodes (2,3) participate in
#               the governing high mode, so scaling them raises the global step.
# --------------------------------------------------------------------------
def _model():
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


def _eigen_f1(extra_mass=None):
    """Fundamental frequency f1 of the (optionally mass-augmented) model.
    extra_mass: dict {node: dm_x} added as nodal mass before the eigensolve."""
    _model()
    if extra_mass:
        for n, dm in extra_mass.items():
            ops.mass(n, dm, 0.0)
    ops.system("FullGeneral")
    lam = ops.eigen("-fullGenLapack", 1)
    return math.sqrt(lam[0]) / (2.0 * math.pi)


def _mode1():
    _model()
    ops.system("FullGeneral")
    lam = ops.eigen("-fullGenLapack", 2)
    assert lam[0] < lam[1], "eigen did not return ascending modes (LAPACK/harness issue)"
    w1 = math.sqrt(lam[0])
    phi = {n: ops.nodeEigenvector(n, 1, 1) for n in FREE}
    pmax = max(abs(v) for v in phi.values())
    return w1, {n: v / pmax for n, v in phi.items()}


def _seed_and_run(integrator_args, dt, T, seed):
    """seed: dict {node: u_x}. Returns (t[], u4[])."""
    _model()
    for n, u in seed.items():
        ops.setNodeDisp(n, 1, u, "-commit")
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


def _maxabs(integrator_args, dt, nsteps, seed):
    t, u = _seed_and_run(integrator_args, dt, nsteps * dt, seed)
    m = 0.0
    for v in u:
        if not math.isfinite(v) or abs(v) > 1.0e3:
            return float("inf")
        m = max(m, abs(v))
    return m


# --------------------------------------------------------------------------
# T-ACC: (A) NECESSITY — at a supra-stable step, plain CD diverges while SMS holds
#        (a broken/no-op/absent SMS diverges here too -> caught). (B) ACCURACY —
#        in the fundamental mode, SMS@supra-dt matches the fine-dt reference.
# --------------------------------------------------------------------------
def test_sms_reference_match():
    # (A) necessity + stability — point-seed excites the high mode so round-off
    # divergence is fast; this is the control that fails on a broken SMS.
    pt = {4: 1.0e-3}
    u_unscaled = _maxabs((CDL,), DT_SUPRA, 400, pt)
    u_sms = _maxabs((SMS, 0.012), DT_SUPRA, 400, pt)
    assert not math.isfinite(u_unscaled) or u_unscaled > 1.0e2, (
        "plain CD must DIVERGE at the supra-stable step dt=%g (else SMS adds nothing "
        "and this test is vacuous); got max|u4|=%r" % (DT_SUPRA, u_unscaled)
    )
    assert math.isfinite(u_sms) and u_sms < 1.0, (
        "SMS must stay stable at dt=%g; got max|u4|=%r" % (DT_SUPRA, u_sms)
    )

    # (B) fidelity at MODEST scaling — at a small added mass the fundamental-mode
    # response tracks the fine-dt reference. Seed mode 1 so the residual is the
    # mass-scaling frequency shift, not high-mode drift; run SMS comfortably INSIDE
    # its stability limit (not at the edge, which would ring the high modes). Part A
    # already proved the supra-stable win at aggressive scaling; this isolates the
    # accuracy cost of a modest, usable scaling.
    w1, phi = _mode1()
    T = 1.5 * (2.0 * math.pi / w1)
    seed = {n: 1.0e-3 * phi[n] for n in FREE}
    t_ref, u_ref = _seed_and_run((CDL,), 0.005, T, seed)             # fine-dt truth
    t_ctl, u_ctl = _seed_and_run((CDL,), 0.008, T, seed)            # SAME dt, NO scaling
    t_sms, u_s = _seed_and_run((SMS, 0.0095), 0.008, T, seed)       # modest scaling, same dt
    assert _rms(u_s) > 1e-6 and _rms(u_ref) > 1e-6, "degenerate run"

    floor = _rel_rms(t_ctl, u_ctl, t_ref, u_ref)   # pure discretization error at dt=0.008
    sms = _rel_rms(t_sms, u_s, t_ref, u_ref)        # discretization + mass-scaling shift
    # faithful: the scaled run stays within engineering tolerance of the converged truth.
    assert sms < 0.045, "SMS@modest-scaling not faithful: rel RMS %.3f%%" % (100 * sms)
    # real injection: the mass-scaling shift adds a MEASURABLE error beyond the no-scaling
    # discretization floor, so a no-op (which would tie the floor) fails this leg. The
    # injection MAGNITUDE is pinned separately + tightly by T-MODAL.
    assert sms > 1.3 * floor, (
        "SMS error %.3f%% not meaningfully above the no-scaling floor %.3f%% — injection "
        "not exercised (a no-op would tie the floor)" % (100 * sms, 100 * floor)
    )


# --------------------------------------------------------------------------
# T-MODAL: SMS must inject the RIGHT mass. Compare SMS's post-injection f1 to the
#          ANALYTIC prediction f1[(s-1)*m on e2's nodes], and show f1 falls
#          monotonically as dtTarget (hence added mass) grows.
# --------------------------------------------------------------------------
def _dt_e2():
    """The element-pencil critical step SMS sizes against (= e2's, the min)."""
    _model()
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 1e-6)
    return ops.criticalTimeStep()


def _f1_after_sms(dtTarget):
    _model()
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(SMS, dtTarget)
    ops.analysis("Transient")
    ops.analyze(1, 0.5 * dtTarget)          # domainChanged -> injection committed
    ops.system("FullGeneral")
    lam = ops.eigen("-fullGenLapack", 1)
    return math.sqrt(lam[0]) / (2.0 * math.pi)


def test_sms_injects_correct_mass_and_lowers_f1():
    dte2 = _dt_e2()
    f1_0 = _eigen_f1()

    # magnitude check: SMS's injected mass must MATCH (s-1)*m_e2 on nodes 2,3.
    # Use a target BETWEEN the tiny (dt_e2~0.0085) and bulk (dt_e~0.010) steps so ONLY
    # e2 is scaled (otherwise the bulk elements scale too and the e2-only prediction
    # undercounts).
    dtT = 0.0095
    s = (dtT / dte2) ** 2
    dm = (s - 1.0) * M_E2_NODE
    f1_pred = _eigen_f1({2: dm, 3: dm})       # analytic: add the predicted mass, eigensolve
    f1_sms = _f1_after_sms(dtT)
    rel = abs(f1_sms - f1_pred) / f1_pred
    assert rel < 0.02, (
        "SMS post-injection f1=%.5g != analytic (s-1)*m prediction f1=%.5g (%.2f%%) "
        "-> injection magnitude wrong" % (f1_sms, f1_pred, 100 * rel)
    )

    # monotone tradeoff: more aggressive dtTarget -> more added mass -> lower f1.
    f1s = [_f1_after_sms(dt) for dt in (0.0095, 0.011, 0.013, 0.016)]
    assert all(f1s[i + 1] < f1s[i] for i in range(len(f1s) - 1)), (
        "f1 must fall monotonically as dtTarget grows (added mass lowers it): %r" % f1s
    )
    assert f1s[0] < f1_0, "even modest scaling should lower f1 below the unscaled value"


# --------------------------------------------------------------------------
# T-HRZ-END2END: on a consistent-mass beam, -lump rowsum produces a NEGATIVE /
#   indefinite rotational lumped mass; the dt_cr eigensolve rejects those indefinite
#   pairs (DGGEV beta-threshold), grossly INFLATING the reported dt_cr (unsafe). hrz
#   and diagonal both give a physical estimate near the true diagonal-mass run limit.
#   (The run itself uses system Diagonal regardless of -lump, so -lump drives only the
#   dt_cr ESTIMATE; the point is that TRUSTING rowsum's number destabilizes a run.)
# --------------------------------------------------------------------------
def _beam(lump):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0); ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 200.0, 80.0, 1.0, 1.0, 1.0, 1,
                "-mass", 1.0, "-cMass")
    ops.setNodeDisp(2, 2, 1.0e-3, "-commit")
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


def _beam_maxtip(dt, nsteps):
    _beam("diagonal")                          # run mass is diagonal-of-consistent regardless
    m = abs(ops.nodeDisp(2, 2))
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return float("inf")
        u = abs(ops.nodeDisp(2, 2))
        if not math.isfinite(u) or u > 1.0e3:
            return float("inf")
        m = max(m, u)
    return m


def test_hrz_rowsum_unsafe_on_consistent_beam():
    dt_rowsum = _beam_dtcr("rowsum")
    dt_diag = _beam_dtcr("diagonal")
    dt_hrz = _beam_dtcr("hrz")

    # rowsum grossly INFLATES dt_cr (indefinite rotational mass -> eigensolve rejects
    # those pairs); hrz and diagonal are both in the physical ballpark.
    assert dt_rowsum > 2.0 * dt_diag, (
        "rowsum dt_cr should be grossly inflated vs the diagonal-mass run limit: "
        "rowsum=%.5g diag=%.5g" % (dt_rowsum, dt_diag)
    )
    ratio = dt_hrz / dt_diag
    assert 1.05 < ratio < 1.35, (
        "hrz/diagonal dt_cr ratio %.3f outside [1.05,1.35] (HRZ's mass-conserving rescale "
        "gives ~1.18 here; a broken HRZ moves it out). hrz=%.5g diag=%.5g"
        % (ratio, dt_hrz, dt_diag)
    )

    # A user who TRUSTS rowsum's dt_cr picks a step that destabilizes the run; one who
    # trusts diagonal's number is safe. (The run mass is diagonal-of-consistent either
    # way; rowsum's estimate is unsafe BECAUSE its indefinite rotational lump inflates it.)
    u_trust_rowsum = _beam_maxtip(0.9 * dt_rowsum, 200)
    u_trust_diag = _beam_maxtip(0.9 * dt_diag, 200)
    assert not math.isfinite(u_trust_rowsum) or u_trust_rowsum > 1.0e2, \
        "running at 0.9*rowsum dt_cr (%.5g) should diverge" % (0.9 * dt_rowsum)
    assert math.isfinite(u_trust_diag) and u_trust_diag < 1.0, \
        "running at 0.9*diagonal dt_cr (%.5g) should be stable, got %r" \
        % (0.9 * dt_diag, u_trust_diag)
