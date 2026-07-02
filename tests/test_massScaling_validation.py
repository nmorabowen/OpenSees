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


# ==========================================================================
# Tier 1 — robustness: convert the Review-3 fixes from one-time probes/warnings
# into asserted regression tests. Plan §"Tier 1". Each warning reaches stderr via
# opserr -> StandardStream -> std::cerr (fd 2), so `capfd` (NOT capsys, which only
# sees Python-level sys.stderr) is the capture fixture.
# ==========================================================================
import re


def _prime_sms(integrator_args, dt=1.0e-4):
    """Run ONE explicit step so CentralDifferenceSMS::domainChanged fires (the
    scaling pass + all its warnings happen there). Returns analyze()'s code."""
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    return ops.analyze(1, dt)


_ADDED_RE = re.compile(r"added mass\s+([0-9.eE+\-]+)% of model mass")


# --------------------------------------------------------------------------
# T-CAP: the -maxAddedMass cap actually fires, and the reported fraction is the
#   REAL (element-mass-denominator) value, not the dead 0% of the SMS-CAP-DEAD bug.
#   Control: the SAME over-scaled model under a LARGE cap must NOT warn -> proves
#   the cap compares (not always-warns) and the denominator is non-zero.
# --------------------------------------------------------------------------
def test_cap_warning_fires_with_real_fraction(capfd):
    # aggressive target -> large added mass (well past any small cap)
    _model()
    capfd.readouterr()                                   # drop setup chatter
    _prime_sms((SMS, 0.02, "-maxAddedMass", 0.001))      # 0.1% cap, ~exceeded
    err_small = capfd.readouterr().err

    _model()
    capfd.readouterr()
    # -verbose so the info line prints even when the cap is NOT tripped (the info line is
    # otherwise emitted only on over/self-report/mismatch) -> the control can positively
    # assert SMS ran and computed a sub-cap fraction.
    _prime_sms((SMS, 0.02, "-maxAddedMass", 50.0, "-verbose"))   # 5000% cap, never exceeded
    err_big = capfd.readouterr().err

    # NB the asserted substring carries the LEADING '%' ("% exceeds ...") on purpose: that
    # '%' straddles the X<<"% exceeds"<<cap message chunks and is exactly what the openseespy
    # PythonStream format-string bug used to eat. So this leg doubly guards (a) the cap fires
    # and (b) the '%'-roundtrip survives -> a PythonStream regression fails HERE, localized.
    assert "% exceeds -maxAddedMass cap" in err_small, (
        "the -maxAddedMass cap WARNING (with its leading '%%') must fire when added mass "
        "exceeds a 0.1%% cap; stderr was:\n%s" % err_small
    )
    # positive control: under a 5000%% cap SMS still RUNS and reports a (sub-cap) fraction
    # ('% of model mass' line present) but does NOT trip the cap. Proves the cap COMPARES
    # (not always-warns) AND that scaling happened (not a silent no-op that trivially passes
    # the 'not in' check).
    assert _ADDED_RE.search(err_big), (
        "under a 5000%% cap SMS must still run and report 'added mass X%% of model mass' "
        "(else the control is a vacuous silence); stderr:\n%s" % err_big
    )
    assert "exceeds -maxAddedMass cap" not in err_big, (
        "the cap must NOT warn under a 5000%% cap (else it always-warns and the test is "
        "vacuous); stderr was:\n%s" % err_big
    )
    # the reported added-mass % must be the REAL element-mass-denominator value, not the
    # dead 0% (SMS-CAP-DEAD: nodal getMass() is 0 on -rho models -> a 0/0 -> 0% cap).
    m = _ADDED_RE.search(err_small)
    assert m is not None, "no 'added mass X%% of model mass' line found:\n%s" % err_small
    frac_pct = float(m.group(1))
    assert frac_pct > 0.1, (
        "reported added-mass fraction %.4g%% must be the real (non-zero) value above the "
        "0.1%% cap; a 0%% here would be the SMS-CAP-DEAD bug (dead nodal denominator)"
        % frac_pct
    )


# --------------------------------------------------------------------------
# T-ENERGY: energy closes under scaling. The EnergyBalanceRecorder's KE sums BOTH
#   element mass AND nodal mass (EnergyBalanceKernel addNodeEnergy -> node->getMass),
#   and SMS injects its fictitious mass at the NODES -> the recorder's KE picks up
#   exactly the augmentation the leap-frog M^-1 uses. Two legs:
#   (1) undamped free vibration: total mechanical energy KE+IE stays conserved under
#       SMS (closure: recorder mass == integrator mass);
#   (2) the SMS run's conserved energy level is HIGHER than the unscaled run's, because
#       the recorder sees the SCALED nodal mass (a recorder blind to nodal mass would
#       tie them -> leg 2 fails: non-vacuous).
# --------------------------------------------------------------------------
def _energy_run(integrator_args, dt, nsteps, v0, efile):
    _model()
    for n in FREE:
        ops.setNodeVel(n, 1, v0, "-commit")              # seed KE in mode-rich state
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
    # injected nodal x-mass per free node (baseline is 0 on the -rho chain, so whatever
    # nodeMass reads here is purely the SMS injection). Read BEFORE remove/wipe.
    inj = {nn: ops.nodeMass(nn, 1) for nn in FREE}
    ops.remove("recorders")                              # flush/close
    rows = []
    with open(efile) as fh:
        for line in fh:
            v = line.split()
            if v:
                rows.append([float(x) for x in v])
    return rows, inj


def test_energy_closes_under_scaling(tmp_path):
    dt, n, v0 = 0.008, 300, 2.0
    r_sms, inj_sms = _energy_run((SMS, 0.016), dt, n, v0, str(tmp_path / "sms.txt"))
    r_cdl, inj_cdl = _energy_run((CDL,), dt, n, v0, str(tmp_path / "cdl.txt"))
    assert len(r_sms) > 20 and len(r_cdl) > 20, "too few energy rows recorded"
    assert all(v == 0.0 for v in inj_cdl.values()), "unscaled CDL must inject NO nodal mass"

    # (1) closure: KE+IE conserved under SMS (recorder mass consistent with the
    #     integrator's scaled M). cols (with -time): time KE IE DW ULW RES ERR%
    mech_sms = [row[1] + row[2] for row in r_sms]
    e_mean = sum(mech_sms) / len(mech_sms)
    drift = (max(mech_sms) - min(mech_sms)) / e_mean
    assert max(row[2] for row in r_sms) > 0.05 * r_sms[0][1], "no KE<->IE exchange (dead IE)"
    # < 5% is the plan's fidelity tolerance; this is ordinary explicit discretization, NOT
    # non-closure. A KE inconsistent with the integrator's (scaled) mass does not merely
    # drift — it accumulates without bound; the magnitude check below is what proves the
    # recorder actually reads the scaled nodal mass.
    assert drift < 0.05, (
        "SMS mechanical-energy drift %.3f%% (expect < 5%%) — gross non-closure (KE read "
        "against the wrong mass) would blow up here, not merely drift" % (100 * drift)
    )

    # (2) PINNED magnitude (not a loose inequality): undamped free vibration conserves KE+IE
    #     at its initial value E0 = 1/2 sum(m_i v_i^2). Seeding only x-velocity v0, the SMS
    #     run's conserved level exceeds the unscaled run's by EXACTLY the injected nodal KE,
    #     1/2 v0^2 sum(injected m_x), IFF the recorder reads the scaled nodal mass. A recorder
    #     blind to nodal mass predicts 0 uplift -> fails; a mis-scaled mass -> wrong magnitude.
    e_sms = sum(mech_sms) / len(mech_sms)
    e_cdl = sum(row[1] + row[2] for row in r_cdl) / len(r_cdl)
    pred_uplift = 0.5 * v0 * v0 * sum(inj_sms.values())
    assert pred_uplift > 0.05 * e_cdl, (
        "test weak: predicted injected-KE uplift %.4g is <5%% of the unscaled energy %.4g — "
        "pick a more aggressive dtTarget" % (pred_uplift, e_cdl)
    )
    rel = abs((e_sms - e_cdl) - pred_uplift) / pred_uplift
    assert rel < 0.15, (
        "SMS conserved energy uplift (e_sms-e_cdl)=%.4g must match the analytic injected KE "
        "%.4g within 15%% (got %.1f%%); off-magnitude => recorder mis-reads the scaled nodal "
        "mass, zero => blind to it" % (e_sms - e_cdl, pred_uplift, 100 * rel)
    )


# --------------------------------------------------------------------------
# T-SELFREP: an element whose stability bound is MASS-INDEPENDENT (self-reported,
#   e.g. a bipenalty LadrunoKinematicCoupling: dt = 2*sqrt(m_p/k)) cannot be helped
#   by adding nodal mass, so SMS must SKIP it (nSelfReport warning) while a normal
#   eigensolve element in the same model IS scaled.
# --------------------------------------------------------------------------
def _selfrep_model(coupling_dtcr):
    """3D model: a massed reference R + a massless slave tied by a bipenalty RBE2
    (self-reports coupling_dtcr), plus a separate tiny stiff truss that SMS scales."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    # self-reporting coupling (test-9 pattern: massed R, massless slave)
    ops.node(1, 0.0, 0.0, 0.0); ops.mass(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    ops.node(2, 1.0, 0.0, 0.0, "-ndf", 3)                # free massless 3-DOF slave
    ops.element("LadrunoKinematicCoupling", 1, 1, 1, 2,
                "-k", 1.0e7, "-bipenalty", "-dtcr", coupling_dtcr)
    # scalable normal element: a short 3-DOF truss (dt_e ~ 1e-3 << dtTarget) with -rho.
    # Use -ndf 3 nodes so the element mass is non-singular (a 6-DOF truss has zero
    # rotational mass -> elementCriticalDt's generalized eigensolve can't size it -> it
    # would be skipped, not scaled).
    ops.node(3, 0.0, 1.0, 0.0, "-ndf", 3); ops.fix(3, 1, 1, 1)
    ops.node(4, 0.1, 1.0, 0.0, "-ndf", 3); ops.fix(4, 0, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)            # c=100, L=0.1 -> dt_e~1e-3
    ops.element("Truss", 2, 3, 4, 1.0, 1, "-rho", 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")


def _truss_only_dt_e():
    """The truss's element-pencil step in isolation (no coupling), to GUARD that the
    fixture's truss really is sub-target on this platform (a drift above dtTarget would
    make it skipped, not scaled, and muddy the self-report assertions)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(3, 0.0, 1.0, 0.0); ops.fix(3, 1, 1, 1)
    ops.node(4, 0.1, 1.0, 0.0); ops.fix(4, 0, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("Truss", 2, 3, 4, 1.0, 1, "-rho", 1.0)
    ops.constraints("Transformation"); ops.numberer("Plain")
    ops.system("Diagonal"); ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl"); ops.analysis("Transient")
    ops.analyze(1, 1e-6)
    return ops.criticalTimeStep()


def test_selfreport_element_skipped_normal_scaled(capfd):
    dtcr = 1.0e-3
    dtTarget = 0.01                                       # > coupling bound AND > truss dt_e
    dt_e = _truss_only_dt_e()
    assert 0.0 < dt_e < dtTarget, (
        "fixture guard: the scalable truss's dt_e=%.4g must be in (0, dtTarget=%.4g) so SMS "
        "scales it; a platform drift here muddies the test, not the feature" % (dt_e, dtTarget)
    )
    _selfrep_model(dtcr)
    capfd.readouterr()
    _prime_sms((SMS, dtTarget, "-verbose"), dt=2.0e-4)    # safe sub-bound step
    err = capfd.readouterr().err

    assert "MASS-INDEPENDENT (self-reported)" in err, (
        "the self-report SKIP warning must fire for the bipenalty coupling; stderr:\n%s" % err
    )
    # the normal truss WAS scaled: SMS injects fictitious mass at its free node (baseline
    # nodal mass is 0 on a -rho truss, so any positive nodeMass is the injection).
    assert ops.nodeMass(4, 1) > 0.0, (
        "the normal truss element must be scaled (injected nodal mass > 0 on node 4); "
        "got %r" % ops.nodeMass(4, 1)
    )
    # the coupling was NOT scaled: its massless slave node 2 receives NO injected mass (SMS
    # skips self-reporting elements before injection; the coupling's own bipenalty lumping
    # lives in the element matrix, not the nodal mass). This is the direct skip detector — a
    # broken SMS that scaled the coupling anyway would push node-2 mass > 0.
    assert ops.nodeMass(2, 1) == 0.0, (
        "the self-reporting coupling must NOT be scaled (slave node 2 stays massless); "
        "got injected nodal mass %r" % ops.nodeMass(2, 1)
    )
    # and exactly one of the two elements was scaled (truss yes, coupling no).
    m = re.search(r"scaled\s+(\d+)/(\d+)\s+elements", err)
    assert m and int(m.group(2)) == 2 and int(m.group(1)) == 1, (
        "expected 'scaled 1/2 elements' (truss scaled, coupling skipped); line:\n%s" % err
    )


# --------------------------------------------------------------------------
# T-CONSTR: the constraint-exclusion guard (ADR-36 v1.1). A sub-target element that
#   touches an MP-constraint (equalDOF/rigidDiaphragm/rigidLink) SLAVE node is EXCLUDED — its
#   injected mass would be redistributed through the constraint's T^T M T and the dt
#   boost would not land. The guard must be TARGETED: it excludes the constrained
#   element (its slave node stays massless) AND reports it still governs, while a FREE
#   sub-target element in the same model is still scaled.
# --------------------------------------------------------------------------
def _constr_model():
    """Two tiny trusses: trussA touches an equalDOF SLAVE node (must be EXCLUDED), trussB
    is free (must be SCALED). Both dt_e ~ 1e-3 << dtTarget."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    # constrained leg: node 2 is the equalDOF slave -> trussA(1-2) excluded
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)
    ops.node(3, 0.2, 0.0); ops.fix(3, 0, 1)              # retained (master) node
    # free leg: node 5 is unconstrained -> trussB(4-5) scaled
    ops.node(4, 0.0, 1.0); ops.fix(4, 1, 1)
    ops.node(5, 0.1, 1.0); ops.fix(5, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)            # c=100, L=0.1 -> dt_e~1e-3 < dtTarget
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)   # touches slave node 2
    ops.element("Truss", 2, 4, 5, 1.0, 1, "-rho", 1.0)   # free
    ops.equalDOF(3, 2, 1)                                 # node-2 x slaved to node-3 (retained)
    ops.constraints("Transformation")
    ops.numberer("Plain")


def test_constrained_element_excluded_free_scaled(capfd):
    _constr_model()
    capfd.readouterr()
    _prime_sms((SMS, 0.01, "-verbose"), dt=2.0e-4)
    err = capfd.readouterr().err

    # the v1.1 disclosure: constrained elements are EXCLUDED and still govern
    assert "EXCLUDED from scaling" in err and "still GOVERN" in err, (
        "the constraint-exclusion warning (excluded + still-govern) must fire; stderr:\n%s" % err
    )
    # guard is TARGETED: the constrained slave node 2 gets NO injection (excluded)...
    assert ops.nodeMass(2, 1) == 0.0, (
        "trussA touches the equalDOF slave node 2 -> it must be EXCLUDED (slave stays "
        "massless); got injected nodal mass %r" % ops.nodeMass(2, 1)
    )
    # ...while the FREE truss's node 5 IS scaled (injection landed). A blanket exclusion
    # that skipped everything would fail this leg.
    assert ops.nodeMass(5, 1) > 0.0, (
        "the free trussB must still be scaled (injected nodal mass > 0 on node 5); got %r"
        % ops.nodeMass(5, 1)
    )
    # exactly one element scaled (the free one); the constrained one excluded, counted.
    m = re.search(r"scaled (\d+)/(\d+)\s+elements", err)
    assert m and int(m.group(1)) == 1 and int(m.group(2)) == 2, (
        "expected 'scaled 1/2 elements' (free trussB scaled, constrained trussA excluded); "
        "line:\n%s" % err
    )


# --------------------------------------------------------------------------
# T-BETAK: stiffness-proportional (betaK) Rayleigh damping shrinks the explicit stable
#   step, so SMS must size against the DAMPED step (closed form s = T^2 + 2*T*c,
#   c = betaK/dt_e) — injecting MORE mass than the undamped s = T^2. A betaK-blind SMS
#   (the pre-v1.1 behavior) injects only T^2 -> under-scales -> still unstable at dtTarget.
#   Single stiff truss so the injected nodal mass is purely this element's.
# --------------------------------------------------------------------------
M_BK_NODE = 0.05          # truss lumped nodal mass = 0.5*rho*A*L = 0.5*1*1*0.1


def _betaK_model(betaK):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)          # free in x
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)        # k=EA/L=1e5, m_node=0.05 -> tiny dt_e
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    if betaK > 0.0:
        ops.rayleigh(0.0, betaK, 0.0, 0.0)           # alphaM=0, betaK(current)=betaK
    ops.constraints("Transformation")
    ops.numberer("Plain")


def _betaK_inject(betaK, dtTarget):
    """Prime SMS one step; return the injected nodal mass on the free node (node 2, x)."""
    _betaK_model(betaK)
    _prime_sms((SMS, dtTarget), dt=2.0e-4)
    return ops.nodeMass(2, 1)


def _betaK_dt_e():
    _betaK_model(0.0)
    ops.system("Diagonal"); ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl"); ops.analysis("Transient")
    ops.analyze(1, 1e-6)
    return ops.criticalTimeStep()


def test_betaK_damped_sizing():
    dtTarget = 0.005
    dt_e = _betaK_dt_e()
    assert 0.0 < dt_e < dtTarget, "fixture: truss dt_e=%.4g must be sub-target %.4g" % (dt_e, dtTarget)
    betaK = 5.0e-4

    dm_u = _betaK_inject(0.0, dtTarget)              # undamped sizing (s = T^2)
    dm_d = _betaK_inject(betaK, dtTarget)            # betaK-damped sizing (s = T^2 + 2*T*c)
    assert dm_u > 0.0, "undamped injection must be positive (element sub-target)"

    # undamped baseline IS the textbook s = (dtTarget/dt_e)^2 (proves the no-damping path
    # is byte-unchanged): injected = (T^2 - 1)*m_node.
    T = dtTarget / dt_e
    su = T * T
    assert abs(dm_u - (su - 1.0) * M_BK_NODE) / ((su - 1.0) * M_BK_NODE) < 0.05, (
        "undamped injection %.5g != textbook (T^2-1)*m=%.5g" % (dm_u, (su - 1.0) * M_BK_NODE)
    )

    # betaK MUST inject MORE (a betaK-blind SMS would tie dm_d to dm_u -> fails here)...
    assert dm_d > 1.02 * dm_u, (
        "betaK damping must increase the injected mass over the undamped sizing (got dm_d=%.5g "
        "<= 1.02*dm_u=%.5g) — betaK term not folded into the scale" % (dm_d, dm_u)
    )
    # ...and by the CLOSED-FORM amount: s_d = T^2 + 2*T*c, c = betaK/dt_e.
    c = betaK / dt_e
    sd = T * T + 2.0 * T * c
    pred_ratio = (sd - 1.0) / (su - 1.0)
    got_ratio = dm_d / dm_u
    assert abs(got_ratio - pred_ratio) / pred_ratio < 0.05, (
        "betaK injected-mass ratio %.4f != closed-form (T^2+2Tc-1)/(T^2-1)=%.4f -> wrong "
        "betaK-damped scale" % (got_ratio, pred_ratio)
    )


# --------------------------------------------------------------------------
# T-BETAK-SLOTS: betaKinit / betaKcomm produce the SAME stiffness-proportional
#   damping force as current-tangent betaK (C = aM*M + bK*K + bK0*K0 + bKc*Kc), so
#   they shrink the explicit stable step identically and MUST enter the damped
#   sizing c = (bK + bK0 + bKc)/dt_e. A slot-blind SMS (pre-fix: only coefs(1) was
#   read) sizes rayleigh(0,0,bKinit,0) models with the UNDAMPED formula ->
#   under-scales -> still unstable at dtTarget. Two legs:
#   (1) equivalence: the injected mass is identical whichever slot carries beta
#       (same closed form, same numbers -> tight tolerance);
#   (2) skip-test: an element whose UNDAMPED step clears dtTarget but whose
#       betaKinit-damped step does not MUST be scaled (the damped skip test must
#       see slot-2/3 damping too).
# --------------------------------------------------------------------------
def _betaK_model_slot(beta, slot):
    """Same single-truss model as _betaK_model, with beta in rayleigh slot 1|2|3
    (1=betaK current, 2=betaKinit, 3=betaKcomm)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    args = [0.0, 0.0, 0.0, 0.0]
    args[slot] = beta
    ops.rayleigh(*args)
    ops.constraints("Transformation")
    ops.numberer("Plain")


def _betaK_inject_slot(beta, slot, dtTarget):
    _betaK_model_slot(beta, slot)
    _prime_sms((SMS, dtTarget), dt=2.0e-4)
    return ops.nodeMass(2, 1)


def test_betaKinit_betaKcomm_enter_damped_sizing():
    dtTarget = 0.005
    dt_e = _betaK_dt_e()
    assert 0.0 < dt_e < dtTarget, "fixture: truss dt_e=%.4g must be sub-target" % dt_e
    beta = 5.0e-4

    # (1) slot equivalence: identical injected mass whichever slot carries beta.
    # At the initial state K == K0 == Kc, so the three sizings are the SAME number
    # (identical code path after the fix) -> tight tolerance, not a loose band.
    dm = {slot: _betaK_inject_slot(beta, slot, dtTarget) for slot in (1, 2, 3)}
    assert dm[1] > 0.0, "current-betaK injection must be positive (fixture)"
    for slot, name in ((2, "betaKinit"), (3, "betaKcomm")):
        rel = abs(dm[slot] - dm[1]) / dm[1]
        assert rel < 1.0e-9, (
            "%s (slot %d) injected %.6g != current-betaK %.6g (rel %.2e) — the slot is "
            "invisible to the damped sizing (pre-fix behavior: undamped T^2 under-scale)"
            % (name, slot, dm[slot], dm[1], rel)
        )

    # (2) skip-test sees slot-2 damping: pick dtTarget2 slightly BELOW dt_e so the
    # element is undamped-safe (skip), then add betaKinit large enough that the
    # damped step dt_e*(sqrt(1+c^2)-c) falls below dtTarget2 -> must be scaled.
    dtTarget2 = 0.9 * dt_e
    assert _betaK_inject_slot(0.0, 1, dtTarget2) == 0.0, (
        "fixture: undamped element must be SKIPPED at dtTarget2=0.9*dt_e (its undamped "
        "step clears the target)"
    )
    beta_big = 0.30 * dt_e            # c = 0.30 -> sqrt(1+c^2)-c ~ 0.744 < 0.9
    dm_skip = _betaK_inject_slot(beta_big, 2, dtTarget2)
    assert dm_skip > 0.0, (
        "an element whose betaKinit-damped step (~%.4g) is below dtTarget2=%.4g MUST be "
        "scaled; got zero injection — the damped SKIP test is blind to betaKinit"
        % (dt_e * (math.sqrt(1.0 + 0.09) - 0.3), dtTarget2)
    )
