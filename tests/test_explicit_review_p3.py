"""Explicit-integrator review PR 3 (2026-07-02): diagnostics batch, scriptable legs.

Failing-first tests for the diagnostics fixes that a Zone-A model can reach:

1. Under selective mass scaling the step-1 dt_cr report compared dt against the
   PRE-SCALING element pencil and printed "expect INSTABILITY" on every correctly
   scaled run (dt <= dtTarget). Post-fix the report is labelled "[PRE-SCALING
   estimate]" and the warning fires only when dt exceeds the POST-SCALING
   effective limit (dtTarget capped by excluded/self-reported elements; the
   Noh-Bathe factor applies on the ExplicitBathe side).

2. The SMS alias parsers downgrade `-recompute` to report-only but did NOT
   consume its trailing N, which fell into the unknown-option branch
   ("unknown option 5" / garbage text under openseespy).

3. The -maxAddedMass cap warning said "% of model mass" while the denominator is
   the total element (translational) lumped mass (nodal getMass() is zero on
   element-rho models, SMS-CAP-DEAD) -- reworded to match.

4. The -divergence KE-proxy breaker compared each step's KE against the
   PREVIOUS step's, so plain free vibration false-tripped at velocity troughs
   (KE ~ 0 near a displacement extreme leaves prevKE ~ eps; the quadratic
   regrowth off that floor is an unbounded ratio -- phase luck decides).
   Post-fix the baseline is the RUNNING MAX of KE: identical for monotonic
   divergence (there the previous step IS the max), immune to troughs.

The kernel-guard hardenings of the same PR (consistent-builder ndf clamp, lumped
sum-ndf pre-walk) are NOT scriptable from a runnable model: any model that
assembles has sum(node ndf) == element numDOF by construction, so the guards
protect only against exotic reduced-basis element masses. They are covered by
the existing SMS batteries staying byte-identical (zero behavior change for
node-major elements) + inspection.

Fixture: 2-node bar, E=A=L=1, rho=2 => m_node = rho*A*L/2 = 1, k = EA/L = 1.
Element pencil lambda_max = 2k/m = 2 => dt_e = 2/sqrt(2) = 1.41421 (undamped,
lumping-invariant for a truss).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DT_E = 2.0 ** 0.5   # the bar element's pre-scaling stable step


def _bar_sms(integrator_args, dt, nsteps=3):
    """Bar fixture with an SMS integrator; returns ops.analyze rc."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 2.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    return ops.analyze(nsteps, dt)


# --------------------------------------------------------------------------
# 1. no false "expect INSTABILITY" on a correctly-scaled run
# --------------------------------------------------------------------------
def test_cd_sms_no_false_instability_warning(capfd):
    # dtTarget = 4 > dt_e = 1.414 -> the bar is scaled; dt = 2 sits ABOVE the
    # pre-scaling limit (old code warned) but BELOW dtTarget (the run is stable).
    capfd.readouterr()
    rc = _bar_sms(("CentralDifferenceSMS", 4.0), dt=2.0)
    err = capfd.readouterr().err
    assert rc == 0, "fixture: dt=2.0 < dtTarget=4.0 must run post-scaling"
    assert "expect INSTABILITY" not in err, (
        "CentralDifferenceSMS warned 'expect INSTABILITY' against the PRE-scaling "
        "pencil on a correctly-scaled run (dt=2.0 <= dtTarget=4.0); stderr:\n%s" % err
    )
    assert "[PRE-SCALING estimate]" in err, (
        "the step-1 dt_cr report must be labelled a PRE-SCALING estimate under "
        "mass scaling; stderr:\n%s" % err
    )


def test_eb_sms_no_false_instability_warning(capfd):
    # dt = 3.0 > pre-scaling Noh-Bathe limit 2*1.414 = 2.83 (old code warned) but
    # < dtTarget = 4.0 (stable post-scaling; post NB limit = 2*4 = 8).
    capfd.readouterr()
    rc = _bar_sms(("ExplicitBathe", 0.54, "-sms", 4.0), dt=3.0)
    err = capfd.readouterr().err
    assert rc == 0, "fixture: dt=3.0 < dtTarget=4.0 must run post-scaling"
    assert "expect INSTABILITY" not in err, (
        "ExplicitBathe -sms warned 'expect INSTABILITY' against the PRE-scaling "
        "Noh-Bathe limit on a correctly-scaled run; stderr:\n%s" % err
    )
    assert "[PRE-SCALING estimate]" in err, (
        "the dt_cr report must be labelled a PRE-SCALING estimate under -sms; "
        "stderr:\n%s" % err
    )


def test_cd_sms_still_warns_above_post_scaling_limit(capfd):
    # dt = 5.0 > dtTarget = 4.0: the run really is unstable -- the warning must
    # survive the reword (now against the POST-scaling effective limit).
    capfd.readouterr()
    _bar_sms(("CentralDifferenceSMS", 4.0), dt=5.0, nsteps=1)
    err = capfd.readouterr().err
    assert "expect INSTABILITY" in err and "POST-SCALING" in err, (
        "dt above dtTarget must still draw the instability warning (against the "
        "post-scaling effective limit); stderr:\n%s" % err
    )


# --------------------------------------------------------------------------
# 2. -recompute downgrade must CONSUME its trailing N
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "args",
    [
        ("CentralDifferenceSMS", 0.1, "-recompute", 5),
        ("CentralDifferenceSMSConsistent", 0.1, "-recompute", 5),
        ("ExplicitBatheSMS", 0.54, 0.1, "-recompute", 5),
        ("ExplicitBatheSMSConsistent", 0.54, 0.1, "-recompute", 5),
    ],
    ids=["cd-sms", "cd-sms-consistent", "eb-sms", "eb-sms-consistent"],
)
def test_sms_recompute_value_consumed(args, capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    capfd.readouterr()
    ops.integrator(*args)
    err = capfd.readouterr().err
    assert "downgraded to REPORT-ONLY" in err, (
        "fixture: -recompute must draw the W1-E3a downgrade note; stderr:\n%s" % err
    )
    assert "unknown option" not in err, (
        "the downgraded -recompute must consume its trailing N (it fell through "
        "to the unknown-option branch); stderr:\n%s" % err
    )


def test_sms_recompute_followed_by_flag(capfd):
    # regression lock: -recompute directly followed by another flag must not eat it.
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    capfd.readouterr()
    ops.integrator("CentralDifferenceSMS", 0.1, "-recompute", "-verbose")
    err = capfd.readouterr().err
    assert "unknown option" not in err, (
        "-recompute followed by a flag must un-read the peeked flag; stderr:\n%s" % err
    )


# --------------------------------------------------------------------------
# 3. -divergence KE breaker must not false-trip at free-vibration troughs
# --------------------------------------------------------------------------
def _sdof_divergence(integrator_args, nsteps=200, dt=0.005):
    """Undamped SDOF (k=100, m=1, omega=10, T=0.628) in free vibration; KE passes
    through ~0 every half period. Returns ops.analyze rc."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeVel(2, 1, 1.0, "-commit")
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    return ops.analyze(nsteps, dt)


@pytest.mark.parametrize(
    "integrator_args",
    [
        ("CentralDifferenceLadruno", "-divergence", 10.0),
        ("ExplicitBathe", 0.54, "-divergence", 10.0),
    ],
    ids=["cd", "eb"],
)
def test_divergence_breaker_no_trough_false_trip(integrator_args, capfd):
    capfd.readouterr()
    rc = _sdof_divergence(integrator_args)
    err = capfd.readouterr().err
    assert rc == 0 and "ABORT" not in err, (
        "the -divergence KE breaker false-tripped on plain free vibration (KE "
        "trough leaves prevKE ~ 0 and the regrowth ratio is unbounded); the "
        "baseline must be the running MAX of KE, not the previous step; "
        "rc=%s, stderr:\n%s" % (rc, err)
    )


def test_divergence_breaker_still_trips_on_growth(capfd):
    # dt above the CD stability limit (dt_cr = 0.2) -> genuine exponential blow-up;
    # the running-max breaker must still fire (before the NaN breaker sees Inf).
    capfd.readouterr()
    rc = _sdof_divergence(("CentralDifferenceLadruno", "-divergence", 10.0),
                          nsteps=200, dt=0.25)
    err = capfd.readouterr().err
    assert rc != 0 and "kinetic-energy proxy" in err, (
        "a genuinely unstable run (dt=0.25 > dt_cr=0.2) must still trip the "
        "-divergence KE breaker; rc=%s, stderr:\n%s" % (rc, err)
    )


# --------------------------------------------------------------------------
# 4. cap warning names the actual denominator
# --------------------------------------------------------------------------
def test_cap_warning_names_element_mass_denominator(capfd):
    # dtTarget=4 -> s = (4/1.414)^2 = 8 -> added mass 7x the element mass >> 5% cap.
    capfd.readouterr()
    _bar_sms(("CentralDifferenceSMS", 4.0), dt=1.0, nsteps=1)
    err = capfd.readouterr().err
    assert "exceeds -maxAddedMass cap" in err, (
        "fixture: 8x scaling must trip the default 5%% cap; stderr:\n%s" % err
    )
    assert "of total element (translational) mass exceeds -maxAddedMass cap" in err, (
        "the cap warning must name its true denominator (total element "
        "translational mass, not 'model mass'); stderr:\n%s" % err
    )
