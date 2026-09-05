"""Fingerprints for `nDMaterial LadrunoSANISAND` (ADR 86, PR-1).

`LadrunoSANISAND` is a thin subclass of the upstream UW `ManzariDafalias`
whose ONLY v1 difference is that the two low-stress constants that
``ManzariDafalias::initialize()`` hardcodes

    m_Pmin      = 1.0e-4 * m_P_atm;    //  0.0101 kPa at P_atm = 101
    m_Presidual = 1.0e-2 * m_P_atm;    //  1.01   kPa at P_atm = 101

become optional deck arguments (``-Presidual`` / ``-Pmin``), are carried on
the wire, and are echoed at construction.  Defaults are ``0.0`` and
``1.0e-3 * P_atm`` (NTUASand02's documented values: a cohesionless sand has
no cohesion).  ``m_Presidual`` is threaded through ~30 mean-stress
evaluations as ``p = one3*GetTrace(stress) + m_Presidual``, including the one
that feeds ``GetPSI()``, so it is an apparent cohesion that displaces the
state parameter psi and with it M^b, M^d and the dilatancy.

Every test below is the fingerprint for one specific SILENT revert-to-vanilla.
Each one fails only if that particular override is missing or wrong; none of
them is a general "does it run" smoke.

  1. ``test_vanilla_equivalence`` -- THE LOAD-BEARING GATE.  Hand
     `LadrunoSANISAND` exactly the two constants vanilla hardcodes and the
     committed stresses must match `ManzariDafalias` to ~1e-12 relative.
     WHAT IT PROVES: along the code paths this strain path exercises (elastic
     stage, the stage-0 -> stage-1 switch through ``Elastic2Plastic()``, and
     20 steps of the ModifiedEuler elastoplastic integrator with an active
     yield surface), the subclass is INERT -- it changes the two constants and
     nothing else, so a difference measured between the two class names is
     attributable to those constants alone.
     WHAT IT DOES NOT PROVE: anything about the paths this one does not visit
     -- other integration schemes (2/3/5/45), the analytical Jacobian branch
     (JacoType), the ``p < m_Pmin + m_Presidual`` clamp, cyclic reversal, or
     the plane-strain wrapper.  It is a strong inertness gate on one path, not
     a proof of inertness everywhere.
     WHAT THE READER MUST BE TOLD ABOUT ITS DECK: both legs run with the
     calibrated ``M_c = 1.3309`` silently raised to 1.99878 -- about +50 % --
     by ``Elastic2Plastic`` at the stage switch, 8 times per leg, one per Gauss
     point, on vanilla ``ManzariDafalias`` too.  This deck's plastic leg does
     not run Gorini's soil.  That is acceptable for THIS test and for no other:
     an inertness gate compares two class names on one deck, both legs take the
     same repair, and it cancels exactly -- inertness across a path that
     traverses the repair is if anything a stronger claim.  It is NOT
     acceptable for any test that reads a difference as "the constant did
     this", and every such test now runs on the confine-first deck instead.
     It carries a SENSITIVITY companion,
     ``test_vanilla_equivalence_is_not_vacuous``, which runs on the
     CONFINE-FIRST deck and makes two comparisons: 1.01 vs 0 (the parser-drop
     fingerprint -- without it, a parser that silently dropped the flags, the
     defect ADR 86 sec.7.1 documents in the sibling ``OPS_SAniSandMSMaterial``,
     would make this test compare vanilla with vanilla and pass while proving
     nothing) and 1.01 vs 0.5 (two interior values, both legs genuinely
     plastic, which isolates the CONSTANT rather than a branch flip).
  2. ``test_echo_and_defaults`` -- the material must SAY what it is running,
     once per construction and not latched behind a static bool (a latched
     message is observable only by whichever test in a session runs first,
     which is exactly how the p_residual defect stayed invisible).  Paired
     with ``test_tcl_subprocess_smoke``, which runs a 3-command deck through
     ``OpenSees.exe``: registration in this fork is FIVE sites and two gates, and a feature wired only into openseespy fails in classic Tcl
     with an error that reads like a broken install.  That deck uses the
     fullest argument form (18 positional + 5 positional optionals + all three
     -flags) with a ``-Presidual`` no default can produce, so it also proves
     the classic-Tcl backend PARSED the flags rather than accepting and
     dropping them.
  3. ``test_revert_to_start_preserves_settings`` -- the fingerprint for the
     ``revertToStart`` override, and the ONLY test that fails without it.
     ``ManzariDafalias::revertToStart`` calls ``this->initialize()``, which
     binds statically to the BASE ``initialize()`` inside the base TU and so
     restores ``p_r = 1.01 kPa`` MID-ANALYSIS.  See the comment in the test for
     why this deck would notice.
  4. ``test_db_roundtrip_carries_presidual`` -- ADR 86 sec.3, the sharp one.
     The base ``Vector(97)`` wire format has NO slot for ``m_Presidual``, and
     ``recvSelf`` never re-runs ``initialize()``, so a restored VANILLA
     material silently runs ``p_r = 0`` while the process beside it runs 1.01.
     This test demonstrates that defect live on `ManzariDafalias` and then
     shows `LadrunoSANISAND` surviving the identical round trip, in both
     directions (``p_r = 0`` and ``p_r = 1.01``), via the extra ``Vector(4)``.
  5. ``test_presidual_is_the_low_p_defect`` -- the physics claim: at
     p0 = 0.01*P_atm the residual pressure makes the model exceed its OWN
     bounding-surface identity eta_end == M^b(psi_end) by about +20 %, and a
     p_r = 0 run does not.  Measured here: +18.13 % vs +0.075 %, with the
     invariant `err * p_end == p_r` holding on both legs.  SLOW TIER.
  6. ``test_planestrain_lane_carries_the_settings`` -- the ONLY test that
     constructs a `LadrunoSANISANDPlaneStrain` at all.  Tests 1-5 all drive
     ndm=3 `stdBrick`, so until this one existed the whole plane-strain lane --
     registered as class tag 33021 in both interpreters and in the broker, i.e.
     reachable by any user today -- had nothing standing behind it.  It is the
     fingerprint for TWO silent reverts:
       * the geotechnical sign flip.  `setTrialStrain` does
         ``mEpsilon = -1.0 * strain_from_element``; changing that ``-1.0`` to
         ``+1.0`` (the defect ADR 86 sec.4.2 calls "catastrophic and silent")
         left the entire battery green.  See the test for the mirrored-deck
         trick that measures what the mutant would return WITHOUT touching C++.
       * a ``getCopy(const char*)`` "PlaneStrain" branch that quietly built a
         vanilla `ManzariDafaliasPlaneStrain`: the same deck is run at
         ``-Presidual`` 1.01 and 0 and the answers must differ materially
         (measured 1.091772e-01), which they cannot if the Gauss-point copy
         never received the constants.
     It runs on the CONFINE-FIRST plane-strain deck.  On the ramp deck it used
     before, the sign fingerprint's discriminating quantity -- the mirrored
     deck's stress deviator -- was +0.187 against the correct deck's -31.132,
     i.e. 0.6 % of it, because that deck's mirror leg spends its whole plastic
     history clamped at the low-p floor.  Confining first makes those -21.643
     and +19.157, a factor of 147 more margin, and drops the clamp warnings to
     zero.  The settings-reach comparison improved 3.178e-2 -> 1.091772e-01.

  6b. ``test_radial_ramp_with_pr0_never_yields`` -- not a fork fingerprint at
     all, but a CHARACTERISATION of upstream behaviour that several tests here
     now depend on.  A proportional-ramp deck run at ``p_residual = 0`` with a
     stage-0 (staged-gravity) leg in front of it gets silent pure elasticity:
     ``ManzariDafalias.cpp:974-977`` re-pins the back-stress ratio from the
     current stress at every elastic step, and on a radial path with p_r = 0
     that freezes alpha at the stress ratio and holds ``f < 0`` for ever.  The
     test pins the two-sided signature -- exactly zero plastic strain with the
     elastic leg, non-zero without it -- and doubles as the positive control
     for the "Outside Bounding" absence the two confined tests assert.
  7. ``test_print_states_what_it_ran`` -- the fingerprint for the ``Print``
     override, which nothing called before.  Deleting it outright left the
     battery green: `test_echo_and_defaults` exercises the CONSTRUCTOR echo
     (``echoLadrunoConstants``), a different function, and since the echo is now
     class-tag-guarded the Gauss-point copies do not echo at all.  Without
     ``Print`` a `printModel` record falls back to `ManzariDafalias::Print`,
     which emits a tag and a type and NOTHING about the two constants -- "a
     record that cannot state what it ran", which is what the override's own
     comment says it exists to prevent.

DETERMINISM / ISOLATION.  ``mElastFlag`` is a **static** member of
``ManzariDafalias`` (ManzariDafalias.cpp:58) and every Manzari-family
constructor writes it, so constructing ANY material of the family resets the
stage flag for EVERY instance in the process (ADR 86 risk 3).  Consequently
every test here starts from ``ops.wipe()``, rebuilds from scratch, and calls
``updateMaterialStage`` EXPLICITLY for its own tag at both stage boundaries;
no test relies on stage state surviving a construction, and no test depends on
execution order.  Test 1 proves this for itself by running its two legs in
BOTH orders and asserting the same answer.

The driver is a ZERO-FREE-DOF material-point cube: a single ``stdBrick`` with
the three negative faces rollered and every positive-face DOF prescribed by an
``sp`` constraint.  There are no free equations, so the answer is the
material's own stress response to a prescribed strain history and cannot be
contaminated by Newton tolerance -- which matters here, because a
displacement-norm convergence test on this material CAN report success without
equilibrium (see the note on test 5).  Test 6 uses the plane-strain TWIN of
that driver -- the same unit square, the same rollers-plus-``sp`` recipe, again
with zero free equations -- so the two lanes are compared on the same protocol
and not on two different analyses.

TWO DECKS, and which test uses which.  The driver above comes in a RAMP form
(``_build``/``_drive``, one proportional strain ramp from zero stress) and a
CONFINE-FIRST form (``_build_confined``/``_drive_confined``: isotropic strain
confinement at stage 0, stage flip at eta = 0, then isochoric deviatoric
loading).  Both are zero-free-DOF; both have plane-strain twins.  Tests 1, 3
and 4 are on the ramp deck and 2 and 6 on the confined one, for reasons each
states in its own docstring.  The short version: a test that asserts an
EQUALITY between two runs of the same deck is indifferent to the deck's
quirks, because they cancel exactly; a test that reads a DIFFERENCE and
attributes it to p_residual is not, and belongs on the confined deck.

CONFINEMENT LEVEL, and the ADR 86 risk-6 exposure it carries.  The confined
deck sits at p0 = 1.7175 kPa after its confinement stage -- that is
0.017*P_atm, squarely in the extreme-low-p territory ADR 86 sec.8 risk 6 says
no gate should depend on.  Three things mitigate it, and they are stated rather
than assumed:
  * the p_residual signal only EXISTS at low p.  Measured on this deck, holding
    everything else fixed and raising the confinement:
        p0 =   1.72 kPa   reldiff(p_r=1.01, p_r=0) = 9.325164e-02
        p0 =   5.73 kPa                              3.167012e-02
        p0 =  17.18 kPa                              1.424105e-02
        p0 =  57.25 kPa                              5.363259e-03
        p0 = 200.38 kPa                              1.451334e-03
    At 200 kPa the signal is barely above the 1e-3 sensitivity floor.  Running
    the sensitivity tests at a comfortable confinement would mean running them
    where there is almost nothing to measure.
  * there is no Newton anywhere in a zero-free-DOF deck, so the failure mode
    risk 6 is about -- an equilibrium iteration that cannot converge at
    vanishing p -- has no equation to fail in.  ``analyze(1)`` still returns 0
    or non-0 and every step is asserted.
  * measured: of the 80-plus confined legs run across the sweeps behind these
    numbers -- both lanes, p_residual 0 / 0.5 / 1.01, mirrored and not, n_dev
    20 / 40 / 80 / 160, p0 from 1.7 to 200 kPa -- not one failed a step.
The one gate that genuinely does depend on an extreme-low-p leg completing is
test 5, and that one is marked slow for exactly that reason.

STEP-SIZE CONVERGENCE -- neither deck has it, and it matters for one kind of
claim only.  Measured, halving the deviatoric step:
    confined 3D, p_r = 1.01, n_dev  20 ->  40   reldiff 4.467130e-03
                                    40 ->  80           1.580339e-02
                                    80 -> 160           1.147915e-02
    ramp 3D, p_r = 1.01, (N_EL,N_PL) x1 -> x2   reldiff 4.094500e-02
                                        x2 -> x4        2.189826e-02
                                        x4 -> x8        1.134525e-02
So the answers are still moving at the ~1e-2 level at the step counts shipped.
The EQUIVALENCE gates (tests 1, 3, 4) are completely unaffected: they compare
runs at identical step counts, so the discretisation error is common to both
sides and cancels exactly -- that is why they hold to 1e-12 and not to 1e-2.
What the wobble does bound is how much a SENSITIVITY number means.  On the
confined deck the 9.325164e-02 sensitivity is about 5.9x the 1.58e-2 wobble at
the nearest refinement, so it is a signal and not discretisation noise.  On the
ramp deck it is NOT comfortable: 5.579874e-02 against a 4.09e-2 wobble is a
factor of 1.4, which is a further reason the sensitivity tests moved.  Neither
number should be quoted as a converged physical quantity; both are quoted as
what this deck reproducibly returns.

SLOW TIER: ``test_presidual_is_the_low_p_defect`` is the only test marked
``@pytest.mark.slow``; measured wall time 28.5 s on the dev box (Windows,
2026-08-27) -- two 1200-step drained triaxial legs at p ~ 1 kPa.  Measured in
the shape CI actually runs (no ``--runslow``), the whole file is 10.32 s, of
which the confine-first decks are ~8.8 s: test 2 5.75 s (3 legs; 6.77 s run
alone) and test 6 3.05 s (3 legs; 4.10 s run alone).  Test 6b costs 0.08 s.
Before those two moved decks the same run was ~1.5 s, so the port costs about
8.8 s.  They are deliberately NOT marked slow: ``tests/conftest.py`` makes that
tier opt-in and no workflow passes ``--runslow`` or sets ``LADRUNO_RUN_SLOW``,
so a slow marker means the test never runs in CI at all -- which is exactly the
defect just fixed for the Tcl smoke.  ~9 s on a ~15-minute Zone-A job is the
right price for the two tests that carry the file's only plastic-regime
measurements.  Test 5 keeps the marker for the reason above -- ADR 86 risk 6,
not cost: 28 s is not minutes-scale, but its p_r = 0 leg fails outright at 400
steps where vanilla survives, and no gate should depend on that.
"""
import math
import os
import re
import subprocess
import tempfile

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_DIST_BIN = os.path.join(_ROOT, "dist", "bin")
_DIST_MP = os.path.join(_ROOT, "dist", "openseesmp")


# ---------------------------------------------------------------------------
#  Gorini's calibrated parameter set (ADR 86 sec.5).  Stress unit is kPa --
#  which matters: the D_factor dilatancy sigmoid ships kPa-dimensional in PR-1
#  (ADR 86 D5a/D5b), so these results are as intended only in kPa.
# ---------------------------------------------------------------------------
_PARAMS = [
    264.32,    # G0
    0.3129,    # nu
    0.6944,    # e_init
    1.33090,   # Mc
    0.71,      # c
    0.027,     # lambda_c
    0.83,      # e0
    0.45,      # ksi
    101.0,     # P_atm
    0.005,     # m
    1.3,       # h0
    0.968,     # ch
    3.5,       # nb
    0.05,      # A0
    5.75,      # nd
    12.5,      # z_max
    1100.0,    # cz
    2.0,       # Rho
]
_P_ATM = _PARAMS[8]

# The two constants ManzariDafalias::initialize() hardcodes.
_VANILLA_PR = 1.0e-2 * _P_ATM      # 1.01   kPa
_VANILLA_PMIN = 1.0e-4 * _P_ATM    # 0.0101 kPa

# LadrunoSANISAND given EXACTLY vanilla's two constants -- the equivalence leg.
_OPTS_VANILLA = ("-Presidual", _VANILLA_PR, "-Pmin", _VANILLA_PMIN, "-honorTolR", 0)
# Same deck, p_residual switched off, p_min PINNED to vanilla's so that exactly
# ONE variable moves (ADR 86 risk 4: two variables at low p if -Pmin is loose).
_OPTS_PR0 = ("-Presidual", 0.0, "-Pmin", _VANILLA_PMIN, "-honorTolR", 0)
# A THIRD, strictly-interior value.  0 is a special number on some paths (it is
# the branch that turns the `if (p > small)` alpha update off, and it is the
# class default), so an A/B that only ever compares 1.01 against 0 cannot tell
# "the constant is being used" from "a branch flipped".  0.5 kPa is neither
# endpoint and is producible by no default, so a leg run at it is a plain
# quantitative sample of the same continuous parameter.  Same `-Pmin` pin.
_OPTS_PR05 = ("-Presidual", 0.5, "-Pmin", _VANILLA_PMIN, "-honorTolR", 0)


# ---------------------------------------------------------------------------
#  Zero-free-DOF material-point driver
#
#  Unit stdBrick cube.  Negative faces rollered (12 DOFs fixed), every
#  positive-face DOF prescribed (12 DOFs sp-constrained) => 0 free equations.
#  Strain path: a single proportional ramp, lateral EXTENSION 0.25*a and axial
#  COMPRESSION a, so the volumetric part compresses (confinement builds) while
#  the deviatoric part grows.  One pattern and one time series only -- no
#  loadConst -- so the whole history can be replayed after a reset() (test 3).
#
#  _N_EL steps are taken with the material at stage 0 (elastic, the community
#  staged-gravity idiom), then the stage is flipped and _N_PL steps are taken
#  with the elastoplastic integrator running.
#
#  TWO MEASURED DEFECTS OF THIS DECK.  An earlier revision of this comment
#  claimed the plastic leg is "genuinely plastic ... so an elastic-only path
#  would show zero sensitivity and prove nothing".  Both halves of that are
#  wrong and the corrections are why the CONFINE-FIRST driver below exists:
#
#   (1) THE p_r = 0 LEG NEVER YIELDS.  Measured on this build,
#       `eleResponse(1,'material',1,'plasticstrains')` after the full history
#       is identically [-0,-0,-0, 0,-0,-0] on the `_OPTS_PR0` leg -- the norm
#       is exactly 0.0, not small.  The mechanism is the radial freeze
#       characterised by `test_radial_ramp_with_pr0_never_yields` below and it
#       is upstream (ManzariDafalias.cpp:974-977).  The `_OPTS_VANILLA` leg
#       DOES yield, but only barely: |eps_pl| = 3.349e-05, three orders below
#       the 5.67e-03 the confined deck reaches.  So the 5.58e-2
#       "sensitivity" this deck reports is a PLASTIC-vs-ELASTIC gap, not a
#       plastic A/B of the constant.  It is a real, reproducible difference and
#       it is fine as a NON-VACUITY witness for the tests that only need to
#       know the two legs are distinguishable (3, 4) -- but it does not isolate
#       p_residual, and the test that claimed to (test 2) has been moved to the
#       confined deck.
#   (2) EVERY LEG ENTERS STAGE 1 OUTSIDE THE BOUNDING SURFACE.  Measured 8
#       times per leg (once per Gauss point), on vanilla `ManzariDafalias` too:
#       "stage-switch stress ratio (1.81707) exceeded the bounding surface M_c
#       (1.3309); raising M_c to 1.99878 (Outside Bounding!)" -- and on the
#       p_r = 0 leg 1.3309 -> 2.3514.  `ManzariDafalias::Elastic2Plastic`
#       inflates the calibrated friction constant by 50-77 % before the plastic
#       leg starts, so this deck's plastic leg does not run Gorini's soil.
#
#  Neither defect touches the EQUIVALENCE gates (tests 1, 3, 4): those compare
#  two runs of the IDENTICAL deck, so whatever the deck does to M_c it does
#  identically to both legs and it cancels exactly.  What the defects break is
#  any claim that a difference measured here isolates the constant.
# ---------------------------------------------------------------------------
_XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
_LAT = 0.25          # lateral extension / axial compression
_E_AX = 3.0e-4       # total axial compressive strain
_N_EL = 5            # stage-0 (elastic) steps
_N_PL = 20           # stage-1 (elastoplastic) steps
_NTOT = _N_EL + _N_PL


def _analysis():
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _NTOT)
    ops.analysis('Static')


def _build(matcmd, tag, opts=()):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(matcmd, tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):                       # rollers on the three negative faces
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    ops.timeSeries('Linear', 1)              # ONE pattern, ONE series: replayable
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, _LAT * _E_AX)
            if y == 1.:
                ops.sp(n, 2, _LAT * _E_AX)
            if k == 1:
                ops.sp(n, 3, -_E_AX)
    _analysis()


def _stress():
    """Committed material stress at Gauss point 1 (element face sign convention)."""
    return list(ops.eleResponse(1, 'material', 1, 'stress'))


def _elastic_leg(tag):
    # EXPLICIT, every time: mElastFlag is static and the constructor of any
    # Manzari-family material just reset it (ADR 86 risk 3).
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(_N_EL):
        assert ops.analyze(1) == 0, f'elastic-stage step {step + 1} failed'


def _plastic_leg(tag):
    ops.updateMaterialStage('-material', tag, '-stage', 1)
    for step in range(_N_PL):
        assert ops.analyze(1) == 0, f'plastic-stage step {step + 1} failed'


def _drive(matcmd, tag, opts=()):
    """wipe -> build -> elastic leg -> stage flip -> plastic leg -> committed stress."""
    _build(matcmd, tag, opts)
    _elastic_leg(tag)
    _plastic_leg(tag)
    return _stress()


def _reldiff(u, v):
    """||u - v|| / ||u|| -- a norm, not a per-component ratio, so a near-zero
    shear component cannot manufacture a large 'relative' difference."""
    nu = math.sqrt(sum(x * x for x in u))
    assert nu > 0.0, 'reference stress is identically zero -- nothing was driven'
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(u, v))) / nu


# ---------------------------------------------------------------------------
#  CONFINE-FIRST zero-free-DOF driver  (ADR 86 PR-2 test-deck rebuild)
#
#  USED BY `test_vanilla_equivalence_is_not_vacuous` (3D) and, through its
#  plane-strain twin below, by `test_planestrain_lane_carries_the_settings`.
#  The other tests deliberately stay on `_drive` -- see the note at the end of
#  this block for which, and why.
#
#  WHY.  `_drive` above runs ONE proportional strain ramp from a zero stress
#  state, and that has two measured defects:
#
#   (1) THE p_r = 0 LEGS NEVER YIELD.  Every stage-0 step re-sets the
#       back-stress ratio (ManzariDafalias.cpp:974-977)
#           p = one3*GetTrace(NextStress) + m_Presidual;
#           if (p > small) NextAlpha = GetDevPart(NextStress)/p;
#       On an exactly radial path (a proportional ramp at fixed K/G) with
#       p_r = 0 the ratio s/p is constant along the ray, so `alpha` is frozen
#       AT the current stress ratio and f = ||s - p*alpha|| - root23*m*p
#       collapses to -root23*m*p < 0 for ever.  Measured on `_drive`:
#       `eleResponse(1,'material',1,'plasticstrains')` is identically
#       [0,0,0,0,0,0] on the p_r = 0 leg.  The 5.58e-2 gap that
#       `test_vanilla_equivalence_is_not_vacuous` USED to measure here was
#       therefore a PLASTIC-vs-ELASTIC gap, not a plastic A/B.  The mechanism
#       is characterised in its own right by
#       `test_radial_ramp_with_pr0_never_yields`.
#   (2) EVERY 3D LEG ENTERS STAGE 1 OUTSIDE THE BOUNDING SURFACE.  Measured 8
#       times per leg (once per Gauss point), on vanilla `ManzariDafalias`
#       too:  "stage-switch stress ratio (1.81707) exceeded the bounding
#       surface M_c (1.3309); raising M_c to 1.99878 (Outside Bounding!)".
#       `ManzariDafalias::Elastic2Plastic` inflates the calibrated friction
#       constant by 50 % before the plastic leg starts, so the deck's plastic
#       leg does not run Gorini's soil.
#       SCOPE, measured, and NARROWER than an earlier draft of this comment
#       claimed: this is a property of the 3D ramp deck only.  The PLANE-STRAIN
#       ramp deck `_drive_ps` fires ZERO "Outside Bounding" on its correct leg
#       -- eps_zz == 0 holds the stage-flip ratio below M_c there.  Its MIRROR
#       leg is the one that misbehaves, and differently: it fires 8 copies of
#       the low-p CLAMP warning ("mean stress p = 0.63108 is below the floor
#       m_Pmin + m_Presidual = 1.0201; CLAMPING ...").  Both go to zero on the
#       confined twin.
#
#  THE FIX is `_drained_triaxial`'s idea -- confine FIRST -- carried onto
#  `_drive`'s zero-free-DOF cube so that no Newton tolerance enters:
#
#    stage 1  isotropic strain confinement at material stage 0.  Isotropic
#             strain on an isotropic elastic material gives an exactly
#             isotropic stress, so GetDevPart(sigma) == 0 and alpha == 0.
#    stage 2  updateMaterialStage -> 1 with eta == 0 << M_c, so
#             Elastic2Plastic finds nothing to repair and M_c is untouched.
#    stage 3  ISOCHORIC (constant-volume, i.e. undrained-triaxial-compression)
#             stepped deviatoric loading: d_eps_zz = -de, d_eps_xx = d_eps_yy
#             = +de/2, hence `_C_LAT = 0.5`.  Isochoric and not the drained
#             0.25 ratio because with zero free DOFs there is no cell pressure
#             to hold p down: a net-compressive strain path drives p to
#             ~9000 kPa (measured) and washes out the p_residual signal, which
#             is a LOW-p effect by construction.  Constant volume keeps the
#             void ratio EXACTLY at e_init - (1+e_init)*tr(eps_conf) -- an
#             independent check that the path is what it claims -- and lets
#             dilatancy alone carry the state up to the bounding surface.
#
#  ZERO FREE DOFs IS PRESERVED.  Both stages are prescribed with `sp`, exactly
#  one per (node, dof) as in `_build`; the two stages differ only in the SHAPE
#  of the ramp, and a shape is a time series.  So the lateral DOFs and the
#  axial DOFs each get their own always-active `Plain` pattern carrying a
#  `Path` series, and NO `loadConst` is used -- the whole history stays a pure
#  function of pseudo-time and is therefore replayable after `ops.reset()`
#  exactly the way `_build`'s single Linear series is.  Verified by execution:
#  run, reset, replay gives a BIT-IDENTICAL committed stress.
#
#  MEASURED (dev box, 2026-08-27, this build), all three legs of test 2:
#
#     leg          p_end     eta/M^b   |eps_pl|    void ratio      "Outside
#                  (kPa)                                           Bounding"
#     p_r = 1.01   9.45882   0.974787  5.665e-03   0.694384750     0
#     p_r = 0.5    9.99822   0.935808  5.672e-03   0.694384750     0
#     p_r = 0      10.62689  0.902355  5.675e-03   0.694384750     0
#
#     p after confinement 1.7175 kPa;  eta at the stage flip 0.0 (exact);
#     void ratio identical to all 9 digits across the three legs, which is the
#     isochoric claim above checking itself.
#     reldiff(p_r=1.01, p_r=0)   = 9.325164e-02   (`_drive`: 5.579874e-02)
#     reldiff(p_r=1.01, p_r=0.5) = 4.241113e-02   (`_drive`: not meaningful --
#                                                  the p_r=0 leg there is
#                                                  elastic, so there is no
#                                                  plastic A/B to take)
#     wall time  ~1.7-3.0 s/leg  (on `_drive`: ~0.003-0.12 s/leg)
#  The cost is the price of the deck: it is spent inside the material's
#  adaptive ModifiedEuler substepping over 0.5 % of axial strain, and measuring
#  it against n_dev 20/40/80/160 shows it barely moves with the LOAD step
#  count -- so it cannot be bought back by taking fewer steps, only by
#  travelling less strain (and then eta/M^b stops short of 1).
#
#  WHICH TESTS MOVED, AND WHICH DID NOT.  Only the two that MEASURE a
#  difference and interpret it as "the constant did this" were ported:
#  `test_vanilla_equivalence_is_not_vacuous` and
#  `test_planestrain_lane_carries_the_settings`.  Tests 1, 3 and 4 stay on
#  `_drive` on purpose:
#    * test 1 (`test_vanilla_equivalence`) is an INERTNESS gate.  It asserts two
#      class names agree to 1e-12 on the same deck; that claim is regime-free,
#      and the ramp deck additionally drags the comparison through the
#      `Elastic2Plastic` M_c-repair path, which is strictly MORE code for the
#      subclass to be inert across, not less.  Porting it would cost ~5 s and
#      weaken it.
#    * tests 3 and 4 (`revertToStart`, the database wire) test a MECHANISM, not
#      a regime, and each needs only that its two reference answers be
#      distinguishable -- which the ramp deck's 5.4e-2 plastic-vs-elastic gap
#      supplies perfectly well.  Test 3 additionally reads `plasticstrains`
#      territory after a `reset()`, where that response is known to carry an
#      offset (see `_plastic_strains`), so moving it would add risk for no gain.
#  Every one of those three now says so in its own docstring rather than
#  leaving the reader to infer it.
# ---------------------------------------------------------------------------
_C_E_CONF = 3.0e-6      # isotropic compressive strain per direction, stage 1
_C_N_CONF = 10          # stage-0 confinement steps
_C_E_AX = 5.0e-3        # further axial compressive strain, stage 3
_C_LAT = 0.5            # lateral extension / axial compression -- 0.5 = isochoric
_C_N_DEV = 40           # stage-1 deviatoric steps
_C_LAT_PS = 1.0         # the plane-strain analogue: eps_zz == 0, so d_eps_xx =
                        # -d_eps_yy is what keeps the volume constant


def _c_series(n_conf=_C_N_CONF, n_dev=_C_N_DEV, e_conf=_C_E_CONF,
              e_ax=_C_E_AX, lat=_C_LAT):
    """The two ramp SHAPES, on a unit pseudo-time grid.

    Both are multiplied by the same `sp` magnitude `-e_conf`, so the shape
    carries the whole path:

        lateral   0 -> 1                    (confine, compression -e_conf)
                  1 -> 1 - lat*e_ax/e_conf  (deviator, EXTENSION)
        axial     0 -> 1                    (confine)
                  1 -> 1 + e_ax/e_conf      (deviator, further COMPRESSION)

    One `sp` per (node, dof) and one series per direction group: no DOF is
    constrained twice, which is what keeps the equation count at zero.

    TRAP FOR WHOEVER MEASURES THIS DECK NEXT.  `n_conf`, `n_dev`, `e_conf`,
    `e_ax` and `lat` are DEFAULT ARGUMENTS, bound to the module constants once
    at def time.  Monkeypatching `_C_N_DEV` or `_C_E_CONF` to sweep the deck
    therefore changes `_deviator_leg`'s loop count and `_build_confined`'s `sp`
    magnitude but NOT the series this function builds, and the two silently
    disagree: with `_C_N_DEV = 80` the analysis walks 80 steps along a 51-point
    Path, runs off its end where `PathSeries` returns 0, and unloads the model
    to a meaningless near-zero stress (observed: reldiff 1.6e+01).  Pass the
    values in as ARGUMENTS and thread them through a local copy of
    `_build_confined` instead.
    """
    s_lat = [i / n_conf for i in range(n_conf + 1)]
    s_ax = list(s_lat)
    r_lat = lat * e_ax / e_conf
    r_ax = e_ax / e_conf
    for i in range(1, n_dev + 1):
        s_lat.append(1.0 - r_lat * i / n_dev)
        s_ax.append(1.0 + r_ax * i / n_dev)
    # One extra HOLD point.  `PathSeries::getFactor` returns 0 at and beyond
    # its LAST time point, so without this the final load step would silently
    # unload the whole model back to zero displacement -- measured: the last
    # step returned sigma ~ 0 and a nonsense eta of 63.6.
    s_lat.append(s_lat[-1])
    s_ax.append(s_ax[-1])
    return s_lat, s_ax


def _analysis_confined():
    """`_analysis` with a unit time step, to match the unit-dt Path series."""
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _build_confined(matcmd, tag, opts=()):
    """The confine-first twin of `_build`: same cube, same rollers, same
    zero free equations, two Path-driven patterns instead of one Linear."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(matcmd, tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):                       # rollers on the three negative faces
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = _c_series()
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)               # lateral DOFs
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -_C_E_CONF)
            if y == 1.:
                ops.sp(n, 2, -_C_E_CONF)
    ops.pattern('Plain', 2, 2)               # axial DOFs
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -_C_E_CONF)
    _analysis_confined()


def _confine_leg(tag):
    # EXPLICIT, every time -- static mElastFlag, ADR 86 risk 3, as everywhere.
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(_C_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'


def _deviator_leg(tag):
    ops.updateMaterialStage('-material', tag, '-stage', 1)
    for step in range(_C_N_DEV):
        assert ops.analyze(1) == 0, f'deviatoric step {step + 1} failed'


def _drive_confined(matcmd, tag, opts=()):
    """wipe -> build -> confine (stage 0) -> stage flip -> shear -> stress."""
    _build_confined(matcmd, tag, opts)
    _confine_leg(tag)
    _deviator_leg(tag)
    return _stress()


def _plastic_strains():
    """Committed plastic strain at Gauss point 1.

    CAVEAT, measured: `revertToStart` restores the total strain but not the
    elastic-strain accumulator this response subtracts, so after an
    `ops.reset()` this reads an OFFSET (exactly minus the previous run's final
    elastic strain) even though the committed STRESS is bit-identical.  Present
    on `_drive` too -- it is upstream, not a property of this deck -- but any
    assertion on plastic strain across a reset has to account for it.
    """
    return list(ops.eleResponse(1, 'material', 1, 'plasticstrains'))


def _state_probe():
    """p, q, eta, M^b and eta/M^b from the committed state.  3D LANE ONLY.

    NOTE M^b is computed from the CALIBRATED `_PARAMS[3]` (M_c = 1.3309).  On
    any deck that fires "Outside Bounding" the material is no longer running
    that M_c -- `Elastic2Plastic` overwrote it -- so eta/M^b measured there is
    against a surface the material has stopped using.  That is precisely why
    the confined deck asserts zero such events.

    REFUSES the plane-strain lane rather than answering wrongly: there the
    'stress' response is the 3-vector (sigma_xx, sigma_yy, sigma_xy) with no
    sigma_zz in it, so p and q are not computable from it -- and read as if it
    were the 3D 6-vector it silently returns a CONSTANT eta of -1.5 on every
    leg, which is a number that looks like a measurement and is not one.
    """
    sig = ops.eleResponse(1, 'material', 1, 'stress')
    assert len(sig) == 6, (
        'p/q/eta are not computable from this response: got %d components, '
        'need the 3D 6-vector. The plane-strain response omits sigma_zz.'
        % len(sig))
    state = ops.eleResponse(1, 'material', 1, 'state')
    s = [-v for v in sig[:3]]                      # compression positive
    p = (s[0] + s[1] + s[2]) / 3.0
    q = s[2] - 0.5 * (s[0] + s[1])
    void_ratio = state[24]
    e_c = _PARAMS[6] - _PARAMS[5] * (max(p, 1e-12) / _P_ATM) ** _PARAMS[7]
    psi = void_ratio - e_c
    m_b = _PARAMS[3] * math.exp(-_PARAMS[12] * psi)
    eta = q / p if p != 0.0 else float('nan')
    return dict(p=p, q=q, eta=eta, e=void_ratio, psi=psi, Mb=m_b, ratio=eta / m_b)


# ===========================================================================
#  1. vanilla equivalence -- the load-bearing gate
# ===========================================================================

_EQ_TOL = 1.0e-12
_SENSITIVITY_FLOOR = 1.0e-3    # measured 9.3e-2 / 4.2e-2 on the confined deck
# A leg that yielded at all is ~5.7e-3; a leg that did not is EXACTLY 0.0, so
# any floor in between separates them.  1e-5 is two orders under the measured
# value and infinitely above the defect's, i.e. it is a "did plasticity happen"
# gate and not a calibration of how much.
_PLASTIC_FLOOR = 1.0e-5
# The diagnostic `ManzariDafalias::Elastic2Plastic` emits when it silently
# raises M_c at a stage switch.  Asserted ABSENT on every confined leg, and
# asserted PRESENT on the ramp deck by
# `test_radial_ramp_with_pr0_never_yields`, so the absence can never go vacuous
# by the message being reworded or removed.
_OB_MARKER = 'Outside Bounding'


def test_vanilla_equivalence():
    """LadrunoSANISAND handed vanilla's two constants IS vanilla, to 1e-12.

    Measured on the dev box: identically 0.0 (bit-for-bit), well inside 1e-12.

    STAYS ON `_drive`, THE RAMP DECK, DELIBERATELY -- and the reader is owed the
    disclosure that this deck runs BOTH legs with M_c inflated by
    `Elastic2Plastic`.  Measured, 8 events per leg (one per Gauss point) and on
    vanilla `ManzariDafalias` too: the stage-switch ratio is 1.81707 against the
    calibrated M_c = 1.3309, so M_c is raised to 1.99878 -- about +50 % -- before
    the plastic leg starts.  This deck's plastic leg therefore does NOT run
    Gorini's calibrated soil.

    That is acceptable HERE and nowhere else in this file, because this test is
    an INERTNESS gate, not a measurement: it asserts that two class names give
    the same answer on the same deck.  Both legs take the same repair, so the
    repair cancels exactly, and inertness across a path that traverses the
    repair is a strictly STRONGER statement than inertness across one that does
    not.  Any test that instead interprets a DIFFERENCE as "p_residual did this"
    has been moved to the confine-first deck, which fires zero such events.

    Order-independence is asserted explicitly because `mElastFlag` is static:
    building the second leg resets the stage flag that the first leg set, so a
    driver that did not re-assert the stage per leg would give an
    order-dependent answer.  Both legs are built in a freshly wipe()d model
    with their own explicit updateMaterialStage calls, and the pair is run in
    BOTH orders; all four results must agree.
    """
    md_first = _drive('ManzariDafalias', 1)
    sani_second = _drive('LadrunoSANISAND', 2, _OPTS_VANILLA)

    sani_first = _drive('LadrunoSANISAND', 3, _OPTS_VANILLA)
    md_second = _drive('ManzariDafalias', 4)

    assert md_first == md_second, (
        'ManzariDafalias gave a different answer depending on whether a '
        'LadrunoSANISAND was constructed before it -- static mElastFlag leaked '
        'across the two legs', md_first, md_second)
    assert sani_first == sani_second, (
        'LadrunoSANISAND gave an order-dependent answer', sani_first, sani_second)

    rel = _reldiff(md_first, sani_first)
    assert rel <= _EQ_TOL, (
        'LadrunoSANISAND given exactly vanilla\'s p_residual and p_min did NOT '
        'reproduce ManzariDafalias -- the subclass is not inert, so no '
        'A/B measurement against it isolates one constant',
        md_first, sani_first, rel)


def test_vanilla_equivalence_is_not_vacuous(capfd):
    """The companion to test 1: a PLASTIC path must be sensitive to p_residual.

    ON THE CONFINE-FIRST DECK (`_drive_confined`), not the ramp deck test 1
    uses, and the difference is the whole point of this test.  Its previous
    revision ran on `_drive`, asserted a 5.58e-2 gap, and its docstring claimed
    that gap proved the strain path "is sensitive to p_residual" on a plastic
    path.  Measured, it did not: on `_drive` the p_r = 0 leg never yields at all
    (|eps_pl| identically 0.0 -- see `test_radial_ramp_with_pr0_never_yields`),
    so 5.58e-2 was the distance between a barely-plastic answer and an ELASTIC
    one.  That is a real difference and it does prove the flags were parsed, but
    it does not isolate the constant, which is what the docstring said it did.

    THREE LEGS, TWO COMPARISONS, PROVING DIFFERENT THINGS.  Keep both:

      p_r = 1.01 vs p_r = 0    THE PARSER-DROP FINGERPRINT.  0 is the class
          default, so if the deck's `-Presidual` were accepted and dropped (ADR
          86 sec.7.1 documents exactly that bug in the sibling
          `OPS_SAniSandMSMaterial`) BOTH legs would run the default, this
          reldiff would collapse to 0, and `test_vanilla_equivalence` would be
          comparing vanilla with vanilla.  Measured 9.325164e-02.

      p_r = 1.01 vs p_r = 0.5  THE CONSTANT ITSELF.  Neither endpoint is a
          default and neither is the `p > small` branch boundary, and BOTH legs
          are genuinely, comparably plastic (|eps_pl| 5.665e-03 and 5.672e-03).
          So this one cannot be explained by a branch flip or by one leg simply
          not yielding: it is a quantitative response to moving a continuous
          material constant.  This is the assertion the old docstring claimed to
          be making and was not.  Measured 4.241113e-02.

    THE |eps_pl| > 0 ASSERTION IS THE GUARD.  It is asserted on every leg,
    and it is the single line that makes the whole class of defect above
    impossible to reintroduce silently: any future edit to the deck that lets a
    leg go elastic -- a different strain path, a smaller strain, a p_r that
    re-freezes alpha -- fails here instead of quietly turning this test back
    into an elastic-vs-plastic comparison wearing a plastic-A/B docstring.

    THE ZERO-M_c-INFLATION GUARD.  `Elastic2Plastic` silently raises M_c when
    the stage-switch stress ratio already exceeds the bounding surface; on the
    ramp deck that fires 8 times per leg and lifts the calibrated 1.3309 to
    1.99878.  Confining first puts eta = 0 (exact) at the flip, so it fires
    zero times -- measured 0 across all three legs -- and this test asserts that
    absence.  That assertion is what makes the confined deck's reason for
    existing something the suite ENFORCES rather than something a comment
    claims.  `test_radial_ramp_with_pr0_never_yields` supplies the matching
    positive control, so a reworded or deleted diagnostic cannot make this
    silently vacuous.

    Only `-Presidual` moves between legs: `-Pmin` is pinned to vanilla's
    1e-4*P_atm in all three (ADR 86 risk 4).  Measured wall time under pytest:
    5.75 s in a full-file run, 6.77 s run alone.
    """
    legs = {}
    for name, tag, opts in (('p_r = 1.01', 5, _OPTS_VANILLA),
                            ('p_r = 0', 6, _OPTS_PR0),
                            ('p_r = 0.5', 9, _OPTS_PR05)):
        legs[name] = (_drive_confined('LadrunoSANISAND', tag, opts),
                      _plastic_strains())

    err = capfd.readouterr().err

    # --- the deck ran the soil it was calibrated for --------------------------
    assert _OB_MARKER not in err, (
        'a confined leg fired the Elastic2Plastic M_c repair. The point of the '
        'confine-first deck is that it enters stage 1 at eta = 0, so nothing '
        'is outside the bounding surface and the calibrated M_c = %g survives '
        'into the plastic leg. If this fires, the numbers below are being '
        'measured on a soil with a different friction angle than the one the '
        'parameter set names' % _PARAMS[3], err)

    # --- every leg is genuinely plastic --------------------------------------
    for name, (_, eps) in legs.items():
        mag = math.sqrt(sum(x * x for x in eps))
        assert mag >= _PLASTIC_FLOOR, (
            f'the {name} leg accumulated |eps_pl| = {mag:.3e}, below the '
            f'{_PLASTIC_FLOOR:g} floor -- it did not meaningfully yield. Any '
            'difference measured against it is then a plastic-vs-elastic gap '
            'and does NOT isolate p_residual, which is exactly the defect this '
            'test was moved off the ramp deck to escape', eps)

    # --- (i) the parser-drop fingerprint: 1.01 vs the class default ----------
    rel_branch = _reldiff(legs['p_r = 1.01'][0], legs['p_r = 0'][0])
    assert rel_branch >= _SENSITIVITY_FLOOR, (
        'switching p_residual from 1.01 kPa to 0 changed the committed stress '
        f'by only {rel_branch:.3e}. Either the deck is not sensitive to '
        'p_residual, or -Presidual was accepted and dropped so both legs ran '
        'the class default -- in which case test_vanilla_equivalence is '
        'comparing vanilla with vanilla',
        legs['p_r = 1.01'][0], legs['p_r = 0'][0])

    # --- (ii) the constant itself: two interior, both plastic ----------------
    rel_const = _reldiff(legs['p_r = 1.01'][0], legs['p_r = 0.5'][0])
    assert rel_const >= _SENSITIVITY_FLOOR, (
        'moving p_residual between two INTERIOR values, 1.01 and 0.5 kPa, with '
        f'both legs plastic, changed the committed stress by only {rel_const:.3e}. '
        'The 1.01-vs-0 comparison above can be explained by a branch flip or by '
        'a default; this one cannot, and without it nothing here shows the '
        'constant is being USED rather than merely SWITCHED',
        legs['p_r = 1.01'][0], legs['p_r = 0.5'][0])


# ---------------------------------------------------------------------------
#  1b. the radial freeze -- why the tests above had to move decks
# ---------------------------------------------------------------------------

def test_radial_ramp_with_pr0_never_yields(capfd):
    """CHARACTERISATION of upstream behaviour a real user can hit.

    NOT a fingerprint for anything this fork wrote, and not a bug report against
    `LadrunoSANISAND`: it pins a property of `ManzariDafalias` that the meaning
    of several tests in this file now turns on, so that it cannot change under
    them silently.

    THE MECHANISM, ManzariDafalias.cpp:974-977.  Every ELASTIC (stage-0) step
    re-sets the back-stress ratio from the current stress:

        p = one3*GetTrace(NextStress) + m_Presidual;
        if (p > small) NextAlpha = GetDevPart(NextStress)/p;

    On an exactly RADIAL stress path -- which a proportional strain ramp at
    fixed K/G is -- the ratio s/p is constant along the ray.  With
    m_Presidual = 0 that `p` is the true mean stress, so `alpha` is pinned to
    the current stress ratio at every step, and the yield function

        f = ||s - p*alpha|| - root23*m*p

    collapses to -root23*m*p < 0 for the whole history.  The stage flip then
    hands the elastoplastic integrator a state that is interior by construction
    and it never yields.  A non-zero m_Presidual breaks the collinearity --
    s/(p_true + p_r) is NOT constant along the ray -- and the freeze lifts.

    WHAT A USER HITS.  A single-element calibration deck that runs staged
    gravity (stage 0) and then a proportional load ramp, with p_residual = 0,
    gets SILENT PURE ELASTICITY: no warning, no non-convergence, a perfectly
    plausible stress history, and no plastic strain anywhere in it.  This is not
    an exotic configuration -- `NTUASand02` (Gorini) ships p_residual = 0 with
    the identical alpha code, and is used for staged foundation problems.  Those
    are safe: the freeze is a property of PROPORTIONAL decks, not of the
    default, and a foundation problem's stress DIRECTION changes as the load
    spreads, which is all it takes to break the collinearity.

    QUALITATIVE ON PURPOSE.  An earlier draft of this check asserted the
    residual yield-function value sits in 0.023-0.024 kPa.  That number tracks
    `m`, the K/G ratio and the step count, so it would fail on any re-tuning of
    the deck for reasons that have nothing to do with the mechanism.  What is
    asserted instead is the mechanism's own two-sided signature, which is exact:
    with the radial elastic leg the plastic strain is EXACTLY zero, and removing
    that leg -- the only change -- makes it non-zero.

    ALSO THE POSITIVE CONTROL for the `_OB_MARKER` guard in the two confined
    tests.  Those assert the diagnostic is ABSENT; if it were reworded or
    removed in C++ their assertions would pass vacuously for ever.  Here the
    same marker is asserted PRESENT, 8 times, on the deck that is known to
    trigger it.  If this half fails, the marker is stale -- fix the string, do
    not relax the guards.

    Cost: 0.08 s for both halves.
    """
    # --- (i) radial ramp + elastic leg + p_r = 0 -> exactly zero plastic strain
    _build('LadrunoSANISAND', 10, _OPTS_PR0)
    _elastic_leg(10)
    _plastic_leg(10)
    frozen = _plastic_strains()
    assert all(v == 0.0 for v in frozen), (
        'the p_r = 0 radial ramp accumulated plastic strain. That is a BETTER '
        'outcome than the one characterised here -- if upstream changed the '
        'alpha update at ManzariDafalias.cpp:974-977 the freeze is gone -- but '
        'it invalidates the reasoning in this file about what the ramp deck '
        'measures, so read it as a notification and re-derive, not as a '
        'regression', frozen)

    ob_after_ramp = capfd.readouterr().err

    # --- (ii) the SAME deck with the elastic leg skipped -> it yields ---------
    # The only change is that stage 1 is asserted before the first step, so the
    # material never takes a stage-0 step and `alpha` is never re-pinned. If the
    # zero above were caused by the deck being too gentle to yield rather than
    # by the freeze, this half would be zero too.
    _build('LadrunoSANISAND', 11, _OPTS_PR0)
    ops.updateMaterialStage('-material', 11, '-stage', 1)
    for step in range(_NTOT):
        assert ops.analyze(1) == 0, f'no-elastic-leg step {step + 1} failed'
    thawed = _plastic_strains()
    mag = math.sqrt(sum(x * x for x in thawed))
    assert mag >= _PLASTIC_FLOOR, (
        f'with the stage-0 leg removed the same strain path still accumulated '
        f'only |eps_pl| = {mag:.3e}. Then the zero above is not the alpha '
        'freeze -- the deck simply never reaches yield -- and the explanation '
        'this file gives for moving tests onto the confined deck is wrong',
        thawed)

    # --- the positive control for the confined decks' zero-M_c-inflation guard
    n_ob = ob_after_ramp.count(_OB_MARKER)
    assert n_ob > 0, (
        'the ramp deck did not emit %r at its stage switch. Measured on this '
        'build it fires 8 times, once per Gauss point, raising M_c from 1.3309 '
        'to 2.3514 on this p_r = 0 leg. If the diagnostic was reworded or '
        'removed, the "not in err" guards in the confined tests are now vacuous '
        'and the marker must be updated there too' % _OB_MARKER,
        ob_after_ramp)


# ===========================================================================
#  2. echo + defaults, and the classic-Tcl wiring
# ===========================================================================

def test_echo_and_defaults(capfd):
    """Constructed with defaults, the material must SAY p_r = 0 and p_min = 1e-3*P_atm.

    The echo is deliberately NOT latched behind a static bool (ADR 86
    sec.4.4), so this test does not need to be first in the file or the
    session -- and it verifies the non-latching itself by constructing twice
    and counting.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 1, *_PARAMS)
    err = capfd.readouterr().err

    assert 'LadrunoSANISAND tag 1' in err, (
        'no per-construction echo at all -- a record that cannot state what it '
        'ran is exactly what the Print/echo overrides exist to prevent', err)
    # EXACT, not a substring prefix.  `'p_residual = 0' in err` also matches
    # `p_residual = 0.101` and `p_residual = 0.0101`, so a copy-paste slip that
    # echoed m_Pmin in m_Presidual's slot passed this assertion -- measured, by
    # mutation, during the pre-merge adversarial review.
    assert re.search(r'p_residual = 0(?![.\d])', err), (
        'the default p_residual is not exactly 0, or is not echoed', err)
    assert '(default, cohesionless)' in err, err
    # 1e-3 * 101.0 = 0.101 kPa.  Asserted as the number AND as the rule, so a
    # change to either the default or the P_atm it is derived from is caught.
    assert 'p_min = %g' % (1.0e-3 * _P_ATM) in err, (
        'the default p_min is not 1e-3*P_atm, or is not echoed', err)
    assert 'default = 1e-3*P_atm' in err, err
    # The echo must still declare where the two D_factor decisions stand.
    # PR-1's wording ("ships UNCHANGED and is still kPa-dimensional in this PR")
    # became FALSE the moment PR-2 non-dimensionalised the sigmoid in vanilla,
    # and sat in the echo unchallenged because this assertion only looked for
    # the substring 'D5a/D5b'.  PR-3 asserts the two claims separately and, for
    # D5b, asserts the CURRENT state rather than a label -- so the next time one
    # of them moves, this test moves with it instead of blessing a stale string.
    assert 'D5b' in err and 'D5a' in err, (
        'the echo no longer declares the D_factor decisions', err)
    assert 'non-dimensionalised in vanilla (PR-2)' in err, (
        'the echo does not state that D5b is DONE. If it still says the sigmoid '
        '"ships UNCHANGED and is still kPa-dimensional", that claim has been '
        'false since PR-2 merged', err)
    assert 'still OPEN' in err, (
        'the echo no longer says D5a is open', err)

    # --- not latched: a second construction must echo again ------------------
    ops.nDMaterial('LadrunoSANISAND', 2, *_PARAMS)
    ops.nDMaterial('LadrunoSANISAND', 3, *_PARAMS)
    err2 = capfd.readouterr().err
    # Count the echo's own LINE SIGNATURE, not a value substring.  This assertion
    # used to count 'p_residual = 0', which is a fragment of free text as well as
    # of the echo: when PR-3 corrected the stale D5b claim, the replacement
    # sentence contained the words "once p_residual = 0" and the count silently
    # went from 2 to 4.  A non-latching check must count MESSAGES, and
    # 'LadrunoSANISAND tag N: p_residual = ' occurs exactly once per echo.
    assert len(re.findall(r'LadrunoSANISAND tag \d+: p_residual = ', err2)) == 2, (
        'the echo is latched (once per process).  A latched message is '
        'observable only by whichever test in a session runs first, which is '
        'exactly how the p_residual defect stayed invisible', err2)

    # --- -honorTolR: 0 and 1 accepted, everything else refused (PR-3) --------
    # PR-1 refused every value but 0, because the base seam did not exist and a
    # flag that claims to have done something it did not do is the class of
    # defect this ADR exists to fix.  PR-2 opened the seam in vanilla
    # (`mHonorTolRInME`, read once in ModifiedEuler); PR-3 wires it, so 1 is now
    # a real request and the refusal moves to values that are neither.
    #
    # The base member the flag drives is a `bool`, so 2 and -1 would both mean
    # "1" silently.  Refusing them is not pedantry: it is the same rule the
    # PR-1 refusal encoded, applied to what is still unrepresentable.
    ops.nDMaterial('LadrunoSANISAND', 4, *_PARAMS, '-honorTolR', 1)
    err3 = capfd.readouterr().err
    assert 'LadrunoSANISAND tag 4' in err3, (
        '-honorTolR 1 was refused. PR-2 landed the seam and PR-3 wires it, so 1 '
        'is a supported value', err3)
    # The echo must name the tolerance the integrator will actually run with --
    # the deck default TolR is 1e-7, so honouring it is a 1000x change from
    # vanilla's hardcoded 1e-4 and must not be reported as merely "flag on".
    assert 'honorTolR = 1' in err3, err3
    assert "this deck's TolR" in err3, (
        'the echo does not say which tolerance ModifiedEuler will honour', err3)
    assert 'TolE = 1e-07' in err3 or 'TolE = 1e-007' in err3, (
        'the echo does not report the honoured TolE as a NUMBER; a flag state '
        'alone does not tell the reader what ran', err3)

    ops.nDMaterial('LadrunoSANISAND', 5, *_PARAMS, '-honorTolR', 0)
    err4 = capfd.readouterr().err
    assert 'honorTolR = 0' in err4 and '0.0001' in err4, (
        'with the flag off the echo must still name vanilla 1e-4 explicitly',
        err4)

    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 6, *_PARAMS, '-honorTolR', 2)
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 7, *_PARAMS, '-honorTolR', -1)
    capfd.readouterr()


# Deliberately the FULLEST deck form: 18 positional doubles + all five
# positional optionals + all three -flags.  That is the shape where the sibling
# `OPS_SAniSandMSMaterial` silently drops TolF/TolR (ADR 86 sec.7.1), and
# `-Presidual 0.5` is a value no default can produce, so the echo asserted
# below proves the CLASSIC-TCL backend really parsed the flag rather than
# merely accepting the command.  (`OPS_GetString` behaves differently in the
# two backends -- see LadrunoSANISAND.cpp:113-125 -- so flag parsing working
# under openseespy does not imply it works here.)
_TCL_PRESIDUAL = 0.5
_TCL_PMIN = 0.0101
# TWO materials, because the two -honorTolR values take different parser
# branches and `OPS_GetString` behaves differently in the classic-Tcl backend
# than in openseespy (LadrunoSANISAND.cpp:150-159).  Tag 2 also passes
# TolR = 1e-06 POSITIONALLY, a value no default produces, so the echoed
# "TolE = 1e-06" proves the backend carried BOTH the positional optional and
# the flag -- which is exactly the pair `OPS_SAniSandMSMaterial` silently
# drops (ADR 86 sec.7.1).
_TCL_DECK = """\
model basic -ndm 3 -ndf 3
nDMaterial LadrunoSANISAND 1 {params} 1 0 1 1e-07 1e-07 \\
    -Presidual {pr} -Pmin {pmin} -honorTolR 0
nDMaterial LadrunoSANISAND 2 {params} 1 0 1 1e-07 1e-06 \\
    -Presidual {pr} -Pmin {pmin} -honorTolR 1
puts "LADRUNO-SANISAND-TCL-OK"
"""


_TCL_INIT_MISSING = "Can't find a usable init.tcl"


def _find_tcl_library():
    """Locate a `tcl8.6` library directory containing `init.tcl`.

    Mirrors the discovery idiom WORKFLOW_GOTCHAS.md sec.7 uses in the CI shell
    (`find ~/.conan2 -name init.tcl -path '*/tcl8.6/*'`), but bounded: it only
    descends into directories whose name hints at Tcl, so it does not walk the
    whole conan cache.  Returns None if nothing is found.
    """
    roots = [os.path.join(os.path.expanduser('~'), '.conan2'),
             os.path.join(_ROOT, 'build')]
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            if 'init.tcl' in filenames and 'tcl8.6' in dirpath.replace('\\', '/'):
                return dirpath
            # prune: only follow branches that can still lead to a tcl runtime
            dirnames[:] = [d for d in dirnames
                           if 'tcl' in d.lower() or d in
                           ('p', 'b', 'lib', 'library', 'share', 'Release', 'res')]
    return None


def _opensees_exe():
    """Locate a classic-Tcl OpenSees binary.

    Searches the dev-box `dist/bin` layout AND the CI build tree.  The second
    location is load-bearing: `.github/workflows/ladruno.yml` builds the
    `OpenSees` target into `build/Release` and copies only `opensees.so` into
    `tests/` -- it never creates `dist/bin`.  Looking only there made this test
    skip on EVERY CI run, which silently disarmed the sole gate on the
    classic-Tcl registration path (the fork's most common silent failure) and
    the only exact assertion on the echoed p_residual/p_min values.  Found by
    mutation testing during the pre-merge adversarial review.
    """
    for root in (_DIST_BIN, os.path.join(_ROOT, "build", "Release")):
        for name in ('OpenSees.exe', 'OpenSees'):
            cand = os.path.join(root, name)
            if os.path.isfile(cand):
                return cand
    return None


def test_tcl_subprocess_smoke(tmp_path):
    """The classic-Tcl registration path must work too.

    Registration in this fork is FIVE sites and two gates; a material wired
    only into `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` works in
    openseespy and fails in `OpenSees.exe` with an error that reads like a
    broken install.  Every other test in this file drives the openseespy gate
    only, so without this one the classic-Tcl gate is untested.
    """
    exe = _opensees_exe()
    if exe is None:
        pytest.skip(f'no OpenSees executable in {_DIST_BIN} or build/Release')

    deck = tmp_path / 'sanisand_smoke.tcl'
    deck.write_text(_TCL_DECK.format(
        params=' '.join(repr(v) for v in _PARAMS),
        pr=_TCL_PRESIDUAL, pmin=_TCL_PMIN))

    def _run(env):
        p = subprocess.run([exe, str(deck)], cwd=str(tmp_path), env=env,
                           capture_output=True, text=True, timeout=180)
        return p, p.stdout + p.stderr

    proc, out = _run(None)

    # WORKFLOW_GOTCHAS.md sec.7: on the CI runner this binary aborts Tcl
    # initialization -- "Can't find a usable init.tcl" -- because TCL_LIBRARY is
    # unset and the conan Tcl runtime is not on its search list.  Retry ONCE with
    # a discovered TCL_LIBRARY rather than skipping: a skip here is what kept this
    # gate silently disarmed on every CI run in the first place.  Retry-on-failure
    # (not set-always) so a working environment is never perturbed.
    if _TCL_INIT_MISSING in out:
        tcl_lib = _find_tcl_library()
        if tcl_lib is None:
            pytest.skip(
                'the Tcl runtime is missing for this binary: it cannot find '
                'init.tcl and no tcl8.6 library dir was found in the conan '
                'cache to point TCL_LIBRARY at. This is an environment gap, '
                'not a defect in LadrunoSANISAND -- see '
                'Ladruno_internal/WORKFLOW_GOTCHAS.md sec.7.')
        env = dict(os.environ, TCL_LIBRARY=tcl_lib)
        proc, out = _run(env)
        assert _TCL_INIT_MISSING not in out, (
            f'TCL_LIBRARY={tcl_lib} did not satisfy this binary', out)

    assert proc.returncode == 0, (
        f'OpenSees.exe exited {proc.returncode} on a 3-command LadrunoSANISAND '
        'deck -- the classic-Tcl registration is broken', out)
    assert 'LADRUNO-SANISAND-TCL-OK' in out, (
        'the deck did not reach its last line', out)
    assert 'LadrunoSANISAND tag 1' in out, (
        'the material constructed but did not echo through the Tcl path', out)
    assert 'p_residual = %g (user)' % _TCL_PRESIDUAL in out, (
        'the classic-Tcl backend did not carry -Presidual through to the '
        'material -- the flag was accepted and silently dropped, which is '
        'ADR 86 sec.7.1 all over again', out)
    assert 'p_min = %g (user)' % _TCL_PMIN in out, (
        'the classic-Tcl backend did not carry -Pmin through', out)

    # ADR-86 PR-3: `-honorTolR 1` through the CLASSIC-TCL parser specifically.
    # Every other assertion on the wired seam drives the openseespy backend, and
    # the two take different paths through the flag loop -- so without this the
    # newly-accepted value is untested on the lane where this fork's
    # registration failures actually happen (ADR 86 sec.4.6).
    assert 'LadrunoSANISAND tag 2' in out, (
        'the second material (-honorTolR 1) did not construct under classic '
        'Tcl', out)
    assert 'honorTolR = 1' in out, (
        'the classic-Tcl backend did not carry -honorTolR 1 through', out)
    # Both exponent spellings accepted, exactly as the openseespy sibling above
    # does: this reads the SAME `opserr << double` path, and MSVC has historically
    # emitted 3-digit exponents (1e-006).  Without the fallback a printf cosmetic
    # would fail this test with a message accusing the parser of the ADR sec.7.1
    # dropped-TolR defect -- a misdiagnosis worse than the cosmetic.
    assert 'TolE = 1e-06' in out or 'TolE = 1e-006' in out, (
        'the honoured ModifiedEuler tolerance is not the deck TolR of 1e-06. '
        'Either the flag or the positional TolR was dropped on the classic-Tcl '
        'lane -- and a dropped TolR is ADR 86 sec.7.1 exactly', out)
    assert 'NO EFFECT' not in out, (
        'IntScheme 1 does reach ModifiedEuler, so the inert-scheme warning must '
        'not fire here', out)
    for bad in ('unknown nDMaterial', 'WARNING invalid', 'Want: nDMaterial'):
        assert bad not in out, (f'classic-Tcl parser rejected the deck ({bad})', out)


# ===========================================================================
#  3. revertToStart
# ===========================================================================

def test_revert_to_start_preserves_settings():
    """reset() must not restore p_r = 1.01 kPa mid-life.

    `ManzariDafalias::revertToStart` (:465-476) calls `this->initialize()`.
    `initialize()` is NOT virtual, so inside the base translation unit that
    call binds statically to `ManzariDafalias::initialize()`, whose :836 is
    `m_Presidual = 1.0e-2 * m_P_atm`.  A subclass that omitted the
    `revertToStart` override would therefore silently go back to running
    1.01 kPa of apparent cohesion from the reset onwards, with nothing
    warning.

    WHY THIS DECK WOULD NOTICE (a test that would pass with or without the
    override is worthless): the strain path ends at p ~ 13.6 kPa, and the
    identical path run with p_r = 1.01 instead of 0 gives a committed stress
    5.4e-2 (5.4 %) away -- see the `_reldiff(first, reverted)` assertion at the
    bottom, which measures that gap in this very test rather than trusting it.
    The pass/fail margin is therefore ~5e-2 against a 1e-12 tolerance.

    STAYS ON `_drive`, AND WHAT THAT GAP ACTUALLY IS.  Measured, this deck's
    p_r = 0 leg never yields -- |eps_pl| is identically 0.0, the radial freeze
    characterised by `test_radial_ramp_with_pr0_never_yields` -- and both legs
    enter stage 1 outside the bounding surface, so `Elastic2Plastic` raises M_c
    (1.3309 -> 2.3514 on the p_r = 0 leg, 1.99878 on the p_r = 1.01 one).  The
    5.4e-2 gap is therefore a plastic-vs-elastic distance measured on a soil
    with an inflated friction constant.

    That is FINE HERE, and the reason is that this test does not interpret the
    gap.  It needs one thing from it: that the two settings are DISTINGUISHABLE,
    so that "the same history replayed after reset() gives the same answer"
    is a claim with content.  Any reproducible difference serves.  What the test
    actually gates is `revertToStart` -- a mechanism, not a regime -- via the
    1e-12 equality above, and that comparison is between two runs of the
    identical deck, so the M_c repair cancels exactly.  Moving it to the
    confine-first deck would cost ~5 s and would additionally have to reason
    about `plasticstrains` reading an offset after a `reset()` (see
    `_plastic_strains`); it would buy nothing.

    The reset is taken MID-LIFE, after a full analysis, which is the case the
    override exists for.  The material stage must be re-asserted after the
    reset for the same static-mElastFlag reason as everywhere else.
    """
    tag = 7
    _build('LadrunoSANISAND', tag, _OPTS_PR0)
    _elastic_leg(tag)
    _plastic_leg(tag)
    first = _stress()

    ops.reset()                 # Domain::revertToStart -> material revertToStart
    ops.setTime(0.0)
    ops.wipeAnalysis()
    _analysis()
    _elastic_leg(tag)
    _plastic_leg(tag)
    second = _stress()

    assert _reldiff(first, second) <= _EQ_TOL, (
        'the same strain history gave a different answer after reset() -- '
        'revertToStart re-ran the BASE initialize() and restored '
        'p_residual = 1.0e-2*P_atm mid-life', first, second)

    # ... and prove the assertion above had something to catch.
    reverted = _drive('LadrunoSANISAND', 8, _OPTS_VANILLA)
    gap = _reldiff(first, reverted)
    assert gap >= _SENSITIVITY_FLOOR, (
        'this deck cannot tell p_r = 0 from p_r = 1.01 kPa, so the reset '
        'assertion above would pass with or without the revertToStart '
        f'override (gap {gap:.3e})')


# ===========================================================================
#  4. FileDatastore round trip -- ADR 86 sec.3
# ===========================================================================

_N_CUT = 12     # plastic steps taken before the save


def _roundtrip(matcmd, tag, opts):
    """Run to mid-path, save, wipe, rebuild, restore, finish.  Returns
    (stress at the save, stress right after the restore, final stress)."""
    with tempfile.TemporaryDirectory(prefix='ladruno_sanisand_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'sanisand_rt')

        _build(matcmd, tag, opts)
        _elastic_leg(tag)
        ops.updateMaterialStage('-material', tag, '-stage', 1)
        for step in range(_N_CUT):
            assert ops.analyze(1) == 0, f'pre-save plastic step {step + 1} failed'
        mid = _stress()

        try:
            ops.database('File', dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f'database() unsupported in this build: {exc}')
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip('database save returned failure on this build')

        _build(matcmd, tag, opts)                      # fresh, uncommitted skeleton
        ops.database('File', dbpath)
        ops.restore(1)
        after = _stress()

        ops.updateMaterialStage('-material', tag, '-stage', 1)
        ops.wipeAnalysis()
        _analysis()
        for step in range(_N_PL - _N_CUT):
            assert ops.analyze(1) == 0, f'post-restore plastic step {step + 1} failed'
        out = _stress()
        ops.wipe()                                     # release FileDatastore handles
        return mid, after, out


def test_db_roundtrip_carries_presidual():
    """A restored LadrunoSANISAND must keep running the soil it was saving.

    The assertion is BEHAVIOURAL, not "restore returned 0": the restored
    material must finish the remaining plastic steps on the p_r = 0 answer and
    NOT on the p_r = 1.01 answer.

    The vanilla leg is asserted too, and it is the fail-before evidence that
    makes the LadrunoSANISAND legs non-vacuous: `ManzariDafalias` restored
    through the identical code path CHANGES its stress the instant it is
    restored (measured 5.5e-2 at the restore, 4.2e-2 at the end of the path),
    because `Vector(97)` has no slot for `m_Presidual` and `recvSelf` never
    re-runs `initialize()`.  If that assertion ever fails, upstream fixed the
    wire format -- read it as a notification, not as a regression here.

    STAYS ON `_drive` for the same reason as test 3: this tests the WIRE, and
    every assertion in it is either an equality between two runs of the identical
    deck (where the deck's quirks cancel exactly) or a "these two are
    distinguishable" non-vacuity check, which any reproducible difference
    satisfies.  The reader should still know what that difference is on this
    deck: the p_r = 0 leg never yields (|eps_pl| identically 0.0, the radial
    freeze), and both legs run with M_c inflated by `Elastic2Plastic`.  So
    "the restored material finished on the p_r = 0 answer" means it finished on
    the ELASTIC answer -- which is a perfectly sharp thing to land on, and is
    arguably a harsher target than a nearby plastic one, but it is not evidence
    about how p_residual behaves in the plastic regime.  That evidence lives in
    `test_vanilla_equivalence_is_not_vacuous` and test 5.
    """
    ref_pr0 = _drive('LadrunoSANISAND', 20, _OPTS_PR0)
    ref_vanilla = _drive('LadrunoSANISAND', 21, _OPTS_VANILLA)
    sanity = _reldiff(ref_vanilla, ref_pr0)
    assert sanity >= _SENSITIVITY_FLOOR, (
        'the two reference answers are indistinguishable, so nothing below '
        'can prove which soil the restored material ran', sanity)

    # --- LadrunoSANISAND with p_r = 0 (the class default) --------------------
    mid, after, out = _roundtrip('LadrunoSANISAND', 22, _OPTS_PR0)
    assert _reldiff(mid, after) <= _EQ_TOL, (
        'the restored p_r = 0 material did not come back on the state it was '
        'saved at', mid, after)
    assert _reldiff(ref_pr0, out) <= _EQ_TOL, (
        'after the round trip the material no longer finishes on the p_r = 0 '
        'answer', ref_pr0, out)
    assert _reldiff(ref_vanilla, out) >= _SENSITIVITY_FLOOR, (
        'after the round trip the material is running the p_r = 1.01 soil')

    # --- LadrunoSANISAND with p_r = 1.01 -------------------------------------
    # This is the leg that actually exercises the new Vector(4): without it the
    # broker-built material would come back at p_r = 0 (P_atm is 0 in the null
    # constructor, so the base initialize() leaves 1.0e-2 * 0.0), i.e. it would
    # land on ref_pr0 instead.
    mid, after, out = _roundtrip('LadrunoSANISAND', 23, _OPTS_VANILLA)
    assert _reldiff(mid, after) <= _EQ_TOL, (mid, after)
    assert _reldiff(ref_vanilla, out) <= _EQ_TOL, (
        'the extra Vector(4) did not carry p_residual = 1.01 across the wire; '
        'the restored material fell back to the base default', ref_vanilla, out)
    assert _reldiff(ref_pr0, out) >= _SENSITIVITY_FLOOR

    # --- vanilla ManzariDafalias: the defect this class exists to fix --------
    mid, after, out = _roundtrip('ManzariDafalias', 24, ())
    assert _reldiff(mid, after) >= _SENSITIVITY_FLOOR, (
        'ManzariDafalias survived the database round trip unchanged. ADR 86 '
        'sec.3 says it cannot: Vector(97) carries m_Pmin at data(96) and '
        'nothing for m_Presidual. If upstream added the slot, ADR 86 sec.3 '
        'and this test both need updating', mid, after)
    assert _reldiff(ref_vanilla, out) >= _SENSITIVITY_FLOOR


def _mp_available():
    return (os.path.isfile(os.path.join(_DIST_MP, 'mpiexec.exe'))
            and os.path.isfile(os.path.join(_DIST_MP, 'openseesmp.pyd')))


def test_db_roundtrip_mp_two_rank_smoke():
    """The parallel half of ADR 86 sec.3, when the box can run it.

    The serial database round trip above already drives sendSelf/recvSelf
    through the same FileDatastore-as-Channel path an MP worker uses, so this
    is a complementary check, not the gate.  It is skipped -- loudly -- rather
    than omitted, so its absence never reads as a pass.
    """
    if not _mp_available():
        pytest.skip(
            'no MP build in this worktree: '
            f'{os.path.join(_DIST_MP, "mpiexec.exe")} and '
            f'{os.path.join(_DIST_MP, "openseesmp.pyd")} are both absent, so '
            'the 2-rank leg cannot run here. The serial FileDatastore round '
            'trip in test_db_roundtrip_carries_presidual covers the same '
            'sendSelf/recvSelf path.')
    # NOT self-arming: an MP build alone is not enough, a 2-rank driver deck
    # still has to be written.  Stated plainly because an earlier draft of the
    # PR text claimed this test would arm itself once an MP build appeared --
    # it will not, and a coverage claim that waits on a file nobody has written
    # is exactly the kind of overclaim this material exists to stop making.
    pytest.skip('MP build present, but the 2-rank LadrunoSANISAND driver deck '
                'has not been written -- this is a PLACEHOLDER, not a gate. '
                'Writing it is owed work, tracked in the manifest row.')


# ===========================================================================
#  5. the physics: p_residual is an apparent cohesion at low confinement
# ===========================================================================
#
#  Instrument: the model's OWN identity eta_end == M^b(psi_end), which the
#  formulation satisfies by construction at the bounding state, so any
#  departure is the instrument and no experiment is needed to see it.  In
#  triaxial compression ||r|| = sqrt(2/3)*eta with eta = q/p, and the bounding
#  condition ||r|| = sqrt(2/3)*M^b, so eta_end == M^b(psi_end) exactly.
#
#  Because every mean-stress evaluation inside the model is p_true + p_r while
#  the measurement uses p_true, the expected departure is
#
#       err = eta_true / M^b - 1 = p_r / p_end       (ADR 86 sec.2's invariant)
#
#  which is why err * p_end reproduces p_r across the whole pressure sweep.
#
#  DRIVER.  Unlike tests 1-4 this leg needs a STRESS boundary condition
#  (constant cell pressure), so it cannot use the zero-free-DOF cube: the
#  lateral faces carry constant nodal loads and the axial is driven by
#  DisplacementControl.  Two consequences, both measured, both ugly:
#
#   * A displacement-norm convergence test on this material CAN report
#     SUCCESS without equilibrium -- observed during development on a
#     stress-controlled variant of this deck, which "converged" every step
#     while the plastic integrator never engaged.  That observation is why
#     tests 1-4 use the zero-free-DOF cube (no Newton at all) instead of a
#     load-driven stack.  CAVEAT, measured 2026-08-26 during the pre-merge
#     adversarial review: the anecdote does NOT reproduce on the deck as
#     shipped below -- at 1e-10 this deck FAILS at axial step 3, and at 1e-6 it
#     completes with the correct +18.12 %.  So do not read the number 1e-8 as a
#     cliff edge; it is simply the value these results were measured at, and
#     the band in the assertions is what actually guards the answer.
#   * Prescribing the axial with an `sp` added mid-analysis forces the node
#     back to zero displacement at load factor 0 and destroys the step.
#     DisplacementControl is used precisely to avoid that.
#   * The p_r = 0 leg is markedly HARDER to converge than vanilla at this
#     pressure -- it needs ~3x the steps (it fails at axial step 1 of 400 and
#     completes at 450; 1200 is used for margin -- the cliff was bracketed
#     between 400 and 450 by sweep during the pre-merge review).  That is not a defect in this class; it is what the
#     residual pressure was silently buying, and it is ADR 86 risk 6 in the
#     flesh.  1200 steps is used for BOTH legs so the protocol is identical.
#
#  MEASURED (dev box, this build, p0 = 1.01 kPa, 2 % axial strain, 1200 steps):
#
#     leg                p_end     eta      M^b      err          err * p_end
#     vanilla (p_r=1.01) 5.4624   2.4456   2.0702   +18.1297 %   0.9903
#     LadrunoSANISAND    3.2936   2.0801   2.0786    +0.0747 %   0.0025
#       (p_r=0)
#
#  (bit-identical across repeated runs in the same process, and across
#  processes.)  err * p_end reproduces p_r on BOTH legs (0.9903 vs 1.0100,
#  1.9 %; and 0.0025 vs 0.0), which is what identifies the departure AS the
#  residual pressure rather than as "the leg had not reached the bounding
#  surface".
#  ADR 86's table reports +20.1 % and err * p_end = 1.0136 at p0 = 1 kPa on an
#  earlier build; +18.13 % here is the same phenomenon at p0 = 1.01 kPa.
#
#  NOTE, and it disagrees with ADR 86 sec.5: the ADR expects the p_r = 0 leg to
#  show a visible residual departure because "the D_factor sigmoid still
#  suppresses 59 % of dilatancy at that pressure".  It does not -- +0.075 % is
#  zero to measurement.  D_factor scales the DILATANCY D, which changes how the
#  state gets to the bounding surface (hence the very different p_end: 3.29 vs
#  5.46) but not the identity eta == M^b that holds once it is there.  The band
#  below is set from the measurement, not from the ADR's prose.

_LOWP_P0 = 0.01 * _P_ATM          # 1.01 kPa, ADR 86 sec.5
_LOWP_EPS_AX = 0.02
_LOWP_STEPS = 1200
_LOWP_TOL = 1.0e-8

# Band for the vanilla leg's departure.  Measured +18.1 % (+18.5 % at 400
# steps).  Lower edge 12 %: the invariant err = p_r / p_end means the number
# tracks p_end, which moves with the step count and the integration scheme --
# 12 % corresponds to p_end = 8.4 kPa, well beyond anything this path reaches.
# Upper edge 30 %: that corresponds to p_end = 3.4 kPa, i.e. the leg stopping
# short of the bounding state; ADR 86's own +52 % row is at p0 = 0.2 kPa, five
# times lower a confinement than this deck.
_LOWP_BAND = (0.12, 0.30)

# Ceiling for the p_r = 0 leg.  Measured +0.075 %; 3 % is 40x that, so this
# tolerates a large drift and still cannot be confused with the vanilla band.
_LOWP_PR0_CEILING = 0.03


def _lowp_probe():
    sig = ops.eleResponse(1, 'material', 1, 'stress')
    state = ops.eleResponse(1, 'material', 1, 'state')
    s = [-v for v in sig[:3]]                      # compression positive
    p = (s[0] + s[1] + s[2]) / 3.0
    q = s[2] - 0.5 * (s[0] + s[1])
    void_ratio = state[24]
    e_c = _PARAMS[6] - _PARAMS[5] * (max(p, 1e-12) / _P_ATM) ** _PARAMS[7]
    psi = void_ratio - e_c
    m_b = _PARAMS[3] * math.exp(-_PARAMS[12] * psi)
    return dict(p=p, q=q, eta=q / p, e=void_ratio, psi=psi, Mb=m_b,
                err=q / p / m_b - 1.0, sxx=s[0])


def _drained_triaxial(matcmd, tag, opts):
    """Isotropic confinement to _LOWP_P0 (elastic stage), then drained
    displacement-controlled axial compression.  Returns the probe dict, or
    None with the reason if a leg fails to complete."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(matcmd, tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    ops.constraints('Plain')          # only homogeneous fixes: no sp anywhere
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', _LOWP_TOL, 60, 0)
    # KrylovNewton, not Newton: same answer to 4 figures (+18.07 % vs +18.12 %
    # on the vanilla leg) and 6x faster, which is what makes this affordable at
    # all.  It DOES need the fine step: at 400 steps KrylovNewton's least-
    # squares acceleration goes degenerate on the p_r = 0 leg ("Intel oneMKL
    # ERROR: Parameter 5 was incorrect on entry to DGELS") and the step fails.
    ops.algorithm('KrylovNewton')

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    f = _LOWP_P0 / 4.0
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            ops.load(n, -f if x == 1. else 0.0, -f if y == 1. else 0.0,
                     -f if k == 1 else 0.0)
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')
    for step in range(10):
        if ops.analyze(1) != 0:
            return None, f'confinement step {step + 1} failed'

    ops.updateMaterialStage('-material', tag, '-stage', 1)
    ops.loadConst('-time', 0.0)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j in range(4):
        ops.load(4 + j + 1, 0.0, 0.0, -0.25)
    ops.integrator('DisplacementControl', 5, 3, -_LOWP_EPS_AX / _LOWP_STEPS)
    ops.analysis('Static')
    for step in range(_LOWP_STEPS):
        if ops.analyze(1) != 0:
            return None, (f'axial step {step + 1}/{_LOWP_STEPS} failed at '
                          f'eps_ax = {ops.nodeDisp(5, 3):.4e}')
    return _lowp_probe(), None


@pytest.mark.slow
def test_presidual_is_the_low_p_defect():
    """At p0 = 0.01*P_atm, p_residual makes the model beat its own M^b.

    Three assertions, in increasing sharpness:

      * the vanilla leg's departure lands in _LOWP_BAND (measured +18.13 %);
      * the p_r = 0 leg's departure is below _LOWP_PR0_CEILING (measured
        +0.075 %) -- "materially less", and in fact zero to measurement;
      * on BOTH legs the invariant `err * p_end == p_residual` holds.  That
        third one is what identifies the departure AS the residual pressure
        rather than as "the leg had not reached the bounding surface yet",
        and it is the only one of the three that a coincidence could not
        produce.

    ADR 86 risk 6 says no gate should depend on an extreme-low-p leg
    completing.  Both legs DO complete here at 1200 steps -- but the p_r = 0
    leg fails at 400 steps (which vanilla survives), so the margin is thin.
    If this becomes flaky, skip it with the measured numbers recorded rather
    than loosening the band into vacuity: the band's job is to distinguish
    +18 % from +0.07 %, and anything that cannot do that is not this test.

    Wall time: 21.4 s on the dev box (Windows, 2026-08-26).  Marked slow for
    the risk-6 reason given in the module docstring, not for the cost.
    """
    vanilla, why = _drained_triaxial('ManzariDafalias', 30, ())
    assert vanilla is not None, (
        'the vanilla low-p drained leg did not complete: ' + str(why) +
        ' -- ADR 86 risk 6 says no gate should depend on an extreme-low-p leg '
        'completing, so if this becomes chronic, skip rather than loosen')

    # p_r = 0, p_min PINNED to vanilla's, so only p_residual moves.
    cohesionless, why = _drained_triaxial('LadrunoSANISAND', 31, _OPTS_PR0)
    assert cohesionless is not None, (
        'the p_r = 0 low-p drained leg did not complete: ' + str(why) +
        ' -- at 400 steps this leg fails where vanilla succeeds, because the '
        'residual pressure is also a numerical regulariser at p ~ 1 kPa '
        '(ADR 86 risk 6). Raise _LOWP_STEPS or skip; do not compare a '
        'completed leg with a stalled one')

    lo, hi = _LOWP_BAND
    assert lo <= vanilla['err'] <= hi, (
        'vanilla ManzariDafalias departed from its own bounding-surface '
        f'identity by {100 * vanilla["err"]:+.2f} %, outside the measured band '
        f'[{100 * lo:.0f} %, {100 * hi:.0f} %]. The band was calibrated against '
        'a measurement on this build (+18.13 % at p_end = 5.462 kPa); if the '
        'real number has moved, move the band and say so -- do not tighten it '
        'around whatever it now reads', vanilla)

    assert abs(cohesionless['err']) <= _LOWP_PR0_CEILING, (
        'the p_residual = 0 leg departed from its own bounding-surface '
        f'identity by {100 * cohesionless["err"]:+.2f} %, above the '
        f'{100 * _LOWP_PR0_CEILING:.0f} % ceiling. Measured on this build: '
        '+0.075 %', cohesionless)

    # The invariant that names the constant: err * p_end == p_residual.
    inv_vanilla = vanilla['err'] * vanilla['p']
    assert abs(inv_vanilla - _VANILLA_PR) <= 0.10 * _VANILLA_PR, (
        f'err * p_end = {inv_vanilla:.4f} does not reproduce p_residual = '
        f'{_VANILLA_PR:.4f}. Without this the departure is not identified as '
        'the residual pressure at all', vanilla)
    inv_pr0 = cohesionless['err'] * cohesionless['p']
    assert abs(inv_pr0) <= 0.10 * _VANILLA_PR, (
        f'err * p_end = {inv_pr0:.4f} on the p_r = 0 leg, which should be 0 '
        'because that is the residual pressure it was given', cohesionless)


# ===========================================================================
#  6. the PLANE-STRAIN lane -- ADR 86 sec.4.2
# ===========================================================================
#
#  ELEMENT CHOICE.  This fork registers both `quad`/`stdQuad` (FourNodeQuad,
#  OpenSeesElementCommands.cpp:719-720) and `SSPquad`/`SSPQuad` (:839-840), and
#  BOTH construct the material with `getCopy("PlaneStrain")`, so either reaches
#  the lane.  Measured on this build, the two agree on the committed stress to
#  ~1e-13 relative on the deck below, so the choice is cosmetic for the physics.
#  `quad` is used because of one non-cosmetic difference: `FourNodeQuad::Print`
#  forwards to `materialPointers[0]->Print` (FourNodeQuad.cpp), whereas
#  `SSPquad::Print` prints only the element tag and its nodes and never touches
#  the material -- so on an SSPquad deck there is no route by which a
#  `printModel` record can state what the plane-strain Gauss point ran.  That
#  costs nothing here (4 Gauss points instead of 1, all seeing the identical
#  uniform strain) and keeps the two lanes' diagnostics symmetric.
#
#  BOUNDARY CONDITIONS.  The zero-free-DOF idea carries over EXACTLY: the unit
#  square's x = 0 edge is rollered in x and its y = 0 edge in y (4 DOFs fixed),
#  and the x = 1 and y = 1 edges have their normal DOF prescribed by `sp`
#  (4 DOFs).  8 of 8 DOFs are constrained, so there are no free equations, the
#  Gauss-point strain is exactly the prescribed uniform field, and no Newton
#  tolerance enters the answer.  No compromise was needed for 2D.
#
#  Strain path: `_build_ps` is the RAMP twin -- lateral EXTENSION 0.25*a in x,
#  axial COMPRESSION a in y, eps_zz = 0 by the kinematics, _N_EL elastic steps
#  then _N_PL elastoplastic ones.  `_build_ps_confined` (below) is the
#  CONFINE-FIRST twin and is what the test now uses; `_build_ps` is kept because
#  it is the deck the previous revision measured and the two are compared here.
#
#  MEASURED ON THE RAMP TWIN `_drive_ps` (dev box, this build;
#  `eleResponse(1,'material',1,'stress')` = sigma_xx, sigma_yy, sigma_xy at the
#  element face):
#
#     leg                                sigma_xx    sigma_yy   sigma_yy-sigma_xx
#     compressive, p_r = 1.01           -10.8022    -41.9341        -31.1318
#     compressive, p_r = 0               -9.9742    -43.0331        -33.0589
#     MIRRORED (== the +1.0 mutant)      -0.9172     -0.7306         +0.1866
#
#  reldiff(p_r = 1.01, p_r = 0) = 3.178e-2, and the mirror margin is +0.1866
#  against -31.13 -- i.e. the mutant is separated from the correct answer by
#  0.19 kPa of deviator, 0.6 % of the correct deck's own deviator.  The
#  confined twin moves that margin to +19.157 against -21.643: 103x more
#  deviator in absolute kPa, and 147x more relative to the correct deck's own.
#  That, plus the fact that this deck's mirror leg spends its whole plastic
#  history CLAMPED at the low-p floor (8 copies of "mean stress p = 0.63108 is
#  below the floor m_Pmin + m_Presidual = 1.0201; CLAMPING ..."), is why the
#  test moved.  Both numbers are recorded because a fingerprint whose margin
#  used to be 0.6 % is worth remembering as such.

_PS_TOL_SIGMA = 1.0e-9      # absolute kPa floor: "distinctly", not "by a hair"


def _build_ps(tag, opts=(), mirror=False):
    """The plane-strain twin of `_build`: same unit square, same
    rollers-plus-`sp` recipe, zero free equations.

    `mirror=True` negates the whole prescribed strain path.  That is not a
    second physics case -- it is how this test measures the ADR 86 sec.4.2
    mutant WITHOUT editing C++.  See `test_planestrain_lane_carries_the_settings`.
    """
    sgn = -1.0 if mirror else 1.0
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    for j, (x, y) in enumerate(_XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    # thickness 1.0, type string "PlaneStrain" -> getCopy("PlaneStrain")
    ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', tag)
    for j, (x, y) in enumerate(_XY):         # rollers on the two negative edges
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j, (x, y) in enumerate(_XY):
        n = j + 1
        if x == 1.:
            ops.sp(n, 1, sgn * _LAT * _E_AX)     # lateral extension
        if y == 1.:
            ops.sp(n, 2, sgn * -_E_AX)           # axial compression
    _analysis()


def _drive_ps(tag, opts=(), mirror=False):
    """wipe -> build 2D -> elastic leg -> stage flip -> plastic leg -> stress."""
    _build_ps(tag, opts, mirror)
    _elastic_leg(tag)
    _plastic_leg(tag)
    return _stress()


# ---------------------------------------------------------------------------
#  The plane-strain twin of the CONFINE-FIRST driver.  USED by
#  `test_planestrain_lane_carries_the_settings`.
#
#  ONE HONEST DIFFERENCE FROM THE 3D DECK, and it cannot be designed away:
#  plane strain pins eps_zz == 0, so an isotropic IN-PLANE strain is NOT an
#  isotropic stress.  With eps = (-e, -e, 0) the elastic answer carries a
#  deviator, and the stage-switch ratio works out to G/K = 3*(1-2*nu)/
#  (2*(1+nu)) = 0.4276 at nu = 0.3129 -- not 0 as it is in 3D.  That is an
#  ALGEBRAIC estimate, not a measurement: the plane-strain 'stress' response
#  has no sigma_zz in it, so the ratio cannot be read back off the deck.  What
#  IS measured is the consequence that matters -- ZERO "Outside Bounding"
#  events, i.e. `Elastic2Plastic` found nothing to repair and M_c was left at
#  its calibrated 1.3309.  "alpha == 0 exactly" is not available on this lane
#  and is not claimed.
#
#  For reference the 3D deck's flip state IS measured, and is exact:
#  sigma = (-1.7175383, -1.7175383, -1.7175383, 0, ~1e-17, ~4e-18), so
#  ||dev(sigma)|| = 1.8e-17 and the ratio `Elastic2Plastic` computes is
#  1.1e-17 -- zero to machine precision, against M_c = 1.3309.
#
#  Constant volume in plane strain means d_eps_xx = -d_eps_yy, i.e.
#  `_C_LAT_PS = 1.0` where the 3D deck uses 0.5.
#
#  THE MIRROR LEG ON THIS DECK, measured, because it is NOT simply the correct
#  deck with signs flipped and the difference has to be stated.  `mirror=True`
#  negates the prescribed path, so at the FACE the confinement stage becomes
#  hydrostatic in-plane TENSION: sigma = (+1.3082, +1.3082, 0) after the confine
#  leg, exactly minus the correct leg's (-1.3082, -1.3082, 0).  Internally that
#  is extension, so at the stage flip `Elastic2Plastic` takes its FIRST branch
#  (GetTrace(mSigma) < 3*m_Pmin) and RESETS the state: measured, the face stress
#  immediately after `updateMaterialStage ... 1` is (-0.0101, -0.0101, 0), i.e.
#  m_Pmin*I, with mAlpha zeroed.  The mirror leg therefore begins its deviatoric
#  stage at the p_min floor, not at 1.72 kPa of confinement.
#
#  THAT DOES NOT WEAKEN THE FINGERPRINT, and the reason is worth being exact
#  about.  The claim is never "the mirror leg is a physically sensible run"; it
#  is "the mirror leg's INTERNAL history is bit-for-bit the +1.0 mutant's on the
#  correct path".  That identity is purely a statement about the wrapper's two
#  negations and about the Gauss-point strain being exactly the prescribed one
#  (zero free DOFs), both of which hold here exactly as on `_build_ps`.  So
#  whatever branch the internal path takes -- including this reset -- is
#  precisely the branch the mutant would take, and exhibiting it is the point.
#  What HAD to be checked, and was: the mirror leg still CONVERGES on this deck
#  (all 10 + 40 steps return 0) and still reverses the deviator.
#
#  MEASURED (dev box, 2026-08-27, this build):
#     sigma after confinement    (-1.3082019, -1.3082019, 0)  isotropic in plane
#     final, p_r = 1.01          (-1.1484102, -22.790964, 0)  |eps_pl| 6.595e-3
#     final, p_r = 0             (-2.1092998, -25.089619, 0)  |eps_pl| 6.609e-3
#     final, p_r = 0.5           (-1.6334022, -23.900908, 0)  |eps_pl| 6.602e-3
#     reldiff(1.01, 0)            1.091772e-01  (on `_drive_ps`: 3.178e-2)
#     reldiff(1.01, 0.5)          5.307990e-02
#     "Outside Bounding"          0 on all four legs, mirror included
#     low-p CLAMP warnings        0 on all four legs (on `_drive_ps` the mirror
#                                 leg fires 8)
#     mirrored (the +1.0 mutant)  (-20.018304, -0.8613164, 0), deviator +19.157
#                                 vs -21.643 on the correct deck.  On
#                                 `_drive_ps` those were +0.187 vs -31.132.
#                                 The mutant/correct separation, measured as
#                                 |dev_mutant| / |dev_correct|, goes from 0.006
#                                 to 0.885 -- a factor of 147.
#     wall time                   ~0.6-1.3 s/leg  (on `_drive_ps`: ~0.01 s/leg)
# ---------------------------------------------------------------------------

def _build_ps_confined(tag, opts=(), mirror=False):
    sgn = -1.0 if mirror else 1.0
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    for j, (x, y) in enumerate(_XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', tag)
    for j, (x, y) in enumerate(_XY):         # rollers on the two negative edges
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    s_lat, s_ax = _c_series(lat=_C_LAT_PS)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)               # x: confine then EXTEND
    for j, (x, y) in enumerate(_XY):
        if x == 1.:
            ops.sp(j + 1, 1, sgn * -_C_E_CONF)
    ops.pattern('Plain', 2, 2)               # y: confine then COMPRESS further
    for j, (x, y) in enumerate(_XY):
        if y == 1.:
            ops.sp(j + 1, 2, sgn * -_C_E_CONF)
    _analysis_confined()


def _drive_ps_confined(tag, opts=(), mirror=False):
    _build_ps_confined(tag, opts, mirror)
    _confine_leg(tag)
    _deviator_leg(tag)
    return _stress()


def test_planestrain_lane_carries_the_settings(capfd):
    """The plane-strain lane must keep the sign convention AND the constants.

    TWO silent reverts, one deck -- now the CONFINE-FIRST deck
    (`_drive_ps_confined`).  Its previous revision ran on `_drive_ps`, the ramp
    twin, and both halves got measurably weaker for it:

      * the mutant margin.  The sign fingerprint separates the correct answer
        from the mutant's by the DEVIATOR.  On the ramp deck those were -31.132
        and +0.187 -- the mutant's deviator was 0.6 % of the correct one, i.e.
        the discriminating quantity was a rounding error away from zero, and it
        was small because that deck's mirror leg spends its whole plastic
        history CLAMPED at the low-p floor (8 copies of the "mean stress p =
        0.63108 is below the floor ...; CLAMPING" warning).  On the confined
        deck they are -21.643 and +19.157: 88.5 %, a factor of 147, and no
        clamp warnings on any leg.
      * the settings-reach comparison, 3.178e-2 -> 1.091772e-01.

    The structure is UNCHANGED -- same two halves, same order, same reasoning.
    Only the deck under them moved.

    (A) THE SIGN CONVENTION.  `LadrunoSANISANDPlaneStrain::setTrialStrain` is
    ported byte-for-byte from `ManzariDafaliasPlaneStrain::setTrialStrain`
    (UWmaterials/ManzariDafaliasPlaneStrain.cpp): the material is
    COMPRESSION-POSITIVE internally and TENSION-POSITIVE at the element face, so
    the wrapper negates on the way in (`mEpsilon(i) = -1.0 * strain(i)`, index
    map {0,1,3}) and negates again on the way out (`getStress`).  Turning the
    inbound `-1.0` into `+1.0` -- ADR 86 sec.4.2's "catastrophic and silent"
    defect -- feeds the model extension where the element applied compression.

    HOW THE MUTANT IS MEASURED WITHOUT TOUCHING C++.  Because the ONLY thing the
    flip changes is the sign of the strain handed to `integrate()`, driving the
    UNMUTATED material with the NEGATED strain path produces bit-for-bit the
    stress the mutant would return on the original path:

        mutant  on path  e :  eps_int = +e  ->  face stress = -sigma_int(+e)
        correct on path -e :  eps_int = +e  ->  face stress = -sigma_int(+e)

    and with zero free equations the Gauss-point strain IS the prescribed path,
    so negating the `sp` values negates it exactly, at every step of the
    history.  `mirror=True` is that run.  It is asserted against, so this test
    does not merely claim the flip would be caught -- it exhibits the flip's
    own output and shows the assertion rejecting it.

    THE ARGUMENT SURVIVES THE MOVE TO THE CONFINED DECK, and it had to be
    checked rather than assumed, because mirroring a CONFINEMENT stage means
    driving hydrostatic in-plane TENSION at the face.  Measured: it does.  The
    identity above depends only on the wrapper's two negations and on the
    Gauss-point strain being exactly the prescribed one, and both hold on
    `_build_ps_confined` exactly as on `_build_ps`.  What the tension does is
    send the internal state into `Elastic2Plastic`'s low-p branch at the stage
    flip, which resets the stress to m_Pmin*I -- measured, face sigma goes
    (+1.3082, +1.3082, 0) -> (-0.0101, -0.0101, 0) across the flip.  That is not
    a defect in the argument: it is a branch the MUTANT would also take, on the
    same internal history, so exhibiting it is exactly what this leg is for.
    The two things that could have broken -- convergence and the deviator
    reversal -- were both verified: all 10 + 40 steps return 0, and the deviator
    comes back +19.157.

    WHY A BARE "the stress is negative" ASSERTION IS NOT ENOUGH, measured: on
    the ramp deck the mutant returned sigma = (-0.9172, -0.7306, 0) and on the
    confined deck it returns (-20.0183, -0.8613, 0) -- STILL compressive at the
    face in both cases.  Fed extension, the sand loses its confinement, `p` pins
    near `m_Pmin + m_Presidual`, and the face sees compression rather than
    tension.  So the sign of sigma alone does not separate the two.  What does
    separate them is the DEVIATOR: under the correct convention the direction
    that was compressed carries the GREATER compression (sigma_yy - sigma_xx =
    -21.643 here), while the mutant reverses that ordering (+19.157) because
    internally it compressed x and extended y.  Both statements are asserted;
    the deviatoric one is the discriminating half.

    (B) THE CONSTANTS REACHING THE GAUSS POINT.  `LadrunoSANISAND::getCopy(const
    char*)` exists only because the base version hardcodes
    `new ManzariDafaliasPlaneStrain(...)`; a PlaneStrain branch that reverted to
    that would build a vanilla wrapper at every Gauss point, and the deck would
    echo the settings once and then run vanilla everywhere -- which the echo
    cannot catch, because `echoLadrunoConstants` is class-tag-guarded and
    Gauss-point copies deliberately stay quiet.  So the same deck is run at
    `-Presidual` 1.01 and at 0 with `-Pmin` PINNED in both (ADR 86 risk 4, one
    variable), and the two answers must differ materially: measured 1.091772e-01.

    GUARDS, as on the 3D confined test: zero `Elastic2Plastic` M_c inflation
    (measured 0 on all three legs run here, mirror included), and |eps_pl| above
    the floor on every leg (measured 6.595e-03 / 6.609e-03 / 6.258e-03), so this
    test can never regress into comparing an elastic leg either.

    Measured wall time under pytest: 3.05 s in a full-file run, 4.10 s alone.
    """
    # --- (B) the settings reach the plane-strain Gauss points ----------------
    ps_vanilla = _drive_ps_confined(40, _OPTS_VANILLA)
    eps_vanilla = _plastic_strains()
    ps_pr0 = _drive_ps_confined(41, _OPTS_PR0)
    eps_pr0 = _plastic_strains()

    rel = _reldiff(ps_vanilla, ps_pr0)
    assert rel >= _SENSITIVITY_FLOOR, (
        'on the PLANE-STRAIN lane, switching p_residual from 1.01 kPa to 0 '
        f'changed the committed stress by only {rel:.3e}. Either this path is '
        'not sensitive to p_residual (in which case nothing below proves the '
        'lane carries the settings), or getCopy("PlaneStrain") is handing the '
        'element a vanilla ManzariDafaliasPlaneStrain that never received them',
        ps_vanilla, ps_pr0, rel)

    # --- (A) the geotechnical sign convention --------------------------------
    sxx, syy, sxy = ps_vanilla

    # tension-positive AT THE FACE: a compressive strain path must not come back
    # as tension.  (True of the mutant too -- see the docstring -- so this is the
    # convention statement, not yet the fingerprint.)
    assert sxx < -_PS_TOL_SIGMA and syy < -_PS_TOL_SIGMA, (
        'the plane-strain lane returned a non-compressive stress for a '
        'compressive strain path -- the element-face sign convention is not '
        'tension-positive', ps_vanilla)
    assert abs(sxy) <= _PS_TOL_SIGMA, (
        'a shear stress appeared on an axis-aligned biaxial path; the {0,1,3} '
        'index map is wrong', ps_vanilla)

    # THE FINGERPRINT: the axis we compressed must carry the greater compression.
    dev = syy - sxx
    assert dev < -_PS_TOL_SIGMA, (
        f'sigma_yy - sigma_xx = {dev:+.4f} >= 0: the axis driven in COMPRESSION '
        '(y) is not the one carrying the greater compression, so the strain '
        'reached the model with its sign inverted -- ADR 86 sec.4.2',
        ps_vanilla)

    # ... and prove that assertion had something to catch, by running the deck
    # whose Gauss-point history is bit-identical to the +1.0 mutant's.
    mutant = _drive_ps_confined(42, _OPTS_VANILLA, mirror=True)
    eps_mutant = _plastic_strains()
    dev_mutant = mutant[1] - mutant[0]
    assert dev_mutant > _PS_TOL_SIGMA, (
        'the mirrored deck -- which reproduces exactly what a +1.0 sign flip in '
        'setTrialStrain would return -- did NOT reverse the deviator '
        f'(sigma_yy - sigma_xx = {dev_mutant:+.4f}), so the assertion above '
        'would pass with or without the flip and this test is not a '
        'fingerprint', mutant)
    assert dev * dev_mutant < 0.0, (
        'the correct and mirrored decks agree on the sign of the stress '
        'deviator, so the sign convention is not observable on this path',
        ps_vanilla, mutant)

    # --- the same two guards the 3D confined test carries --------------------
    err = capfd.readouterr().err
    assert _OB_MARKER not in err, (
        'a confined plane-strain leg fired the Elastic2Plastic M_c repair, so '
        'the numbers above were measured on a soil whose friction constant was '
        'silently raised above the calibrated %g' % _PARAMS[3], err)

    for name, eps in (('p_r = 1.01', eps_vanilla), ('p_r = 0', eps_pr0),
                      ('mirror', eps_mutant)):
        mag = math.sqrt(sum(x * x for x in eps))
        assert mag >= _PLASTIC_FLOOR, (
            f'the plane-strain {name} leg accumulated |eps_pl| = {mag:.3e}, '
            f'below the {_PLASTIC_FLOOR:g} floor -- it did not meaningfully '
            'yield, so the comparisons above are not between two plastic runs',
            eps)


# ===========================================================================
#  7. Print -- a record that can state what it ran
# ===========================================================================
#
#  ROUTE, established by probing rather than assumed.  An `nDMaterial` lives in
#  the material registry, NOT in the Domain, and no interpreter -- openseespy or
#  classic Tcl -- has a `print material <tag>` command, so the dimensionless
#  PROTOTYPE's `Print` is unreachable from a deck.  What IS reachable is the
#  Gauss-point copy: `printModel` walks the Domain, `Brick::Print`
#  (brick/Brick.cpp:437) forwards to `materialPointers[0]->Print(s, flag)`, and
#  that lands in `LadrunoSANISAND::Print` on the LadrunoSANISAND3D copy.  That is
#  the better target anyway -- it is the object that actually integrated the
#  soil.  Measured routes on this build:
#
#     ops.printModel('-file', path)            -> whole domain, reaches Print
#     ops.printModel('-file', path, '-ele')    -> elements only, reaches Print
#     ops.printModel('-ele')                   -> same, to stderr
#     ops.printModel('-ele', 1)                -> FAILS, "print ele failed to
#                                                 get integer" (upstream arg
#                                                 off-by-one, see the test)
#
#  The `-file` form is used so the assertion reads a file rather than racing the
#  C++ `opserr` stream against pytest's capture.
#
#  WHY THIS CANNOT PASS ON A VANILLA OBJECT OR ON DEFAULTS.  `ManzariDafalias::
#  Print` (UWmaterials/ManzariDafalias.cpp:730-734) is two lines -- a tag and a
#  type -- and mentions neither constant, so deleting the override drops both
#  asserted strings entirely.  And the two values asserted, 0.5 and 0.25, are
#  producible by no default: the defaults are p_residual = 0 and
#  p_min = 1e-3*P_atm = 0.101.

_PRINT_PR = 0.5      # no default produces this
_PRINT_PMIN = 0.25   # nor this (default is 1e-3*P_atm = 0.101)


def test_print_states_what_it_ran(tmp_path):
    """`printModel` must report the constants this material actually ran.

    This is NOT `test_echo_and_defaults` in another costume.  That test drives
    `echoLadrunoConstants()`, which fires once per `nDMaterial` command from the
    CONSTRUCTOR and -- since it became class-tag-guarded, to stop ~400k lines of
    stderr on a 50k-element model -- fires only on the deck-level prototype.
    `Print` is a different function, deliberately left UNGUARDED so that any
    instance, prototype or Gauss-point copy, can still be interrogated.  Nothing
    called it before this test: deleting the override outright left the whole
    battery green, which is how "a record that cannot state what it ran" -- the
    exact failure the override's own comment cites -- would have shipped.
    """
    tag = 50
    _build('LadrunoSANISAND', tag,
           ('-Presidual', _PRINT_PR, '-Pmin', _PRINT_PMIN))

    out = tmp_path / 'printmodel.out'
    # '-file <path>' then '-ele': OPS_printModel opens the file, then
    # printElement() with no tags left dumps every element to it.  A tag list
    # ('-ele', 1) is NOT usable from openseespy -- printElement() rewinds with
    # OPS_ResetCurrentInputArg(2), which assumes the Tcl argv where index 0 is
    # the command name, and in the Python backend that skips past the tag and
    # aborts with "print ele failed to get integer".  Upstream, and reported
    # rather than worked around.
    ops.printModel('-file', str(out), '-ele')
    txt = out.read_text()

    assert 'LadrunoSANISAND Material, tag: %d' % tag in txt, (
        'printModel produced no LadrunoSANISAND record at all. If the header '
        'says "ManzariDafalias Material" instead, the Print override is gone '
        'and the record cannot state what it ran', txt)
    assert 'ManzariDafalias Material, tag:' not in txt, (
        'the base ManzariDafalias::Print ran, so the override is not in effect',
        txt)
    # The Gauss-point copy, not the dimensionless prototype -- see the note
    # above for why the prototype is unreachable from any interpreter.
    assert 'Type: ThreeDimensional' in txt, txt

    # BOTH constants, at their non-default values.  Regexes, not `in`, so the
    # resolved value and the echoed user input are asserted TOGETHER: a Print
    # that reported the input but ran something else, or vice versa, fails here.
    assert re.search(r'p_residual = %g\s+\(input %g, user\)'
                     % (_PRINT_PR, _PRINT_PR), txt), (
        'printModel did not report p_residual = %g as a user value' % _PRINT_PR,
        txt)
    assert re.search(r'p_min\s+= %g\s+\(input %g, user\)'
                     % (_PRINT_PMIN, _PRINT_PMIN), txt), (
        'printModel did not report p_min = %g as a user value' % _PRINT_PMIN,
        txt)

    # ... and they are not the defaults wearing the same words.
    assert 'default: cohesionless' not in txt, (
        'Print reported the p_residual DEFAULT on a material given '
        '-Presidual %g' % _PRINT_PR, txt)
    assert 'default = 1e-3*P_atm' not in txt, (
        'Print reported the p_min DEFAULT on a material given -Pmin %g'
        % _PRINT_PMIN, txt)


# ===========================================================================
#  8. -honorTolR -- the PR-2 seam, wired (ADR 86 PR-3)
# ===========================================================================
#
#  WHAT THE SEAM IS.  `ManzariDafalias::ModifiedEuler()` opened with a bare
#  `double TolE = 1e-4` and ignored `mTolR` outright, so a deck passing a tight
#  TolR on IntScheme 1 silently ran substep error control at 1e-4.  RungeKutta45
#  and SAniSandMS both honour mTolR; ModifiedEuler was the outlier.  PR-2 made
#  the literal a seam -- `TolE = mHonorTolRInME ? mTolR : 1e-4` -- with
#  `mHonorTolRInME` a protected base member set false in all four
#  ManzariDafalias constructors and written NOWHERE in that class, so vanilla is
#  bit-identical.  PR-3 makes `LadrunoSANISAND::applyLadrunoConstants()` the one
#  writer, driven by the deck's `-honorTolR`.
#
#  TWO CLAIMS, AND ONLY ONE OF THEM IS A MEASUREMENT.
#
#  (a) STRUCTURAL, not measured: at TolR == 1e-4 the flag is a no-op on ANY deck
#      and ANY strain path, because both operands of the conditional are then the
#      same double and TolE is bit-identical.  `test_honor_tolr_is_inert_at_1e_4`
#      exhibits it on one deck; the generality comes from the source, and what
#      would falsify it is a SECOND read site for the flag appearing.  That is a
#      one-line check:
#          grep -n mHonorTolRInME SRC/material/nD/UWmaterials/ManzariDafalias.cpp
#      must show the flag DECLARED, set false in the four constructors, and READ
#      in ModifiedEuler -- nowhere else.
#
#  (b) MEASURED, and deck-specific: at TolR < 1e-4 the flag moves the answer.
#      On the confine-first deck, IntScheme 1:
#          TolR 1e-4  reldiff(flag 0, flag 1) = 0.000000e+00   (structural)
#          TolR 1e-5                            2.783e-04
#          TolR 1e-6                            2.641e-04
#          TolR 1e-7                            2.580e-04
#      Wall cost of turning it on: 1.9-2.3 s/leg against 1.8-2.8 s/leg off, i.e.
#      inside the run-to-run noise on this deck.  It is NOT free in general --
#      TolE is a substep tolerance and dT_min is 1e-6, so a tight TolR on a
#      harder path buys accuracy with substeps.
#
#  (c) MEASURED ONCE, and only HALF of it is gated -- read this before quoting it.
#          IntScheme  3 (RungeKutta4)   reldiff(flag 0, flag 1) = 0.000000e+00
#          IntScheme 45 (RungeKutta45)  reldiff(flag 0, flag 1) = 0.000000e+00
#      That is not a defect -- it is the seam's scope -- but a flag that is
#      accepted, stored, echoed and inert is the exact defect class this ADR
#      exists to fix, so the material WARNS at construction.  IntScheme 45
#      already honours TolR unconditionally, which is why the seam was needed for
#      ModifiedEuler and not for it.
#
#      THE GAP, named rather than papered over (found by adversarial review):
#      `test_honor_tolr_warns_when_the_scheme_cannot_honour_it` below constructs
#      IntScheme 1 and IntScheme 45 ONLY.  The IntScheme 3 half of the line above
#      is a one-off measurement plus a source argument (RungeKutta4 never calls
#      ModifiedEuler -- traced through `explicit_integrator`'s switch), NOT a
#      standing gate.  It is deliberately not gated here: `ManzariDafalias`'s
#      scheme-3/5 "no error control" warning is behind a C++ static latch that
#      fires once per PROCESS, and `test_manzari_safety_pack.py::
#      test_scheme3_warns_no_error_control` depends on being the first test in the
#      session to construct such a material -- its own comment says so.  Adding an
#      IntScheme 3 construction to this file would break that test in one file
#      order and not the other, which is a worse defect than the one it would
#      close.  A gate for it belongs in a SUBPROCESS, like the classic-Tcl smoke.

_HONOR_TOLR_TIGHT = 1.0e-6      # measured 2.641e-04 against flag 0
_HONOR_TOLR_FLOOR = 1.0e-5      # 26x under the measured 2.641e-04. NOT "two orders":
                                # an earlier comment said so and was wrong by a decade.
_HONOR_TOLR_VANILLA_TOLE = 1.0e-4


def _opts_tolr(tolr, honor):
    """Positional optionals + flags, with ONLY the honorTolR flag free.

    IntScheme 1 (ModifiedEuler -- the one function that reads the seam),
    TanType 0, JacoType 1, TolF 1e-7.  p_residual and p_min pinned to vanilla's
    so the honorTolR A/B moves exactly one variable (ADR 86 risk 4, applied to
    the third constant).
    """
    return [1, 0, 1, 1.0e-7, tolr,
            '-Presidual', _VANILLA_PR, '-Pmin', _VANILLA_PMIN,
            '-honorTolR', honor]


def test_honor_tolr_is_inert_at_1e_4():
    """TolR == 1e-4 => the flag selects the same double either way.

    This is the seam's INERTNESS half and it is structural: PR-2 wrote the
    conditional so that a false flag selects the literal `1e-4` verbatim, with
    no arithmetic.  Asking for TolR = 1e-4 makes the TRUE branch select the same
    value, so TolE -- and therefore all of ModifiedEuler -- is bit-identical.

    Bit-identical is asserted, not 'small': anything above 0.0 here means the
    flag is doing something other than choosing between mTolR and 1e-4.
    """
    a = _drive_confined('LadrunoSANISAND', 1, _opts_tolr(_HONOR_TOLR_VANILLA_TOLE, 0))
    b = _drive_confined('LadrunoSANISAND', 1, _opts_tolr(_HONOR_TOLR_VANILLA_TOLE, 1))
    assert _reldiff(a, b) == 0.0, (
        'at TolR = 1e-4 the two operands of `mHonorTolRInME ? mTolR : 1e-4` are '
        'the same double, so the flag cannot change anything. A nonzero '
        'difference means the seam is reading or scaling something else',
        a, b, _reldiff(a, b))


def test_honor_tolr_changes_the_answer_when_tolr_is_tight():
    """The load-bearing gate: `-honorTolR 1` is not a decoration.

    Same deck, same TolR, ONLY the flag moves.  If the wiring line in
    `applyLadrunoConstants()` is deleted, `mHonorTolRInME` stays false in every
    constructor, TolE stays 1e-4 on both legs and this reldiff collapses to
    exactly 0.0 -- which is what `test_honor_tolr_is_inert_at_1e_4` above
    asserts for the case where 0.0 is CORRECT.  The two tests are a pair: one
    pins where the flag must do nothing, the other where it must do something,
    and no single mutation satisfies both.
    """
    off = _drive_confined('LadrunoSANISAND', 1, _opts_tolr(_HONOR_TOLR_TIGHT, 0))
    on = _drive_confined('LadrunoSANISAND', 1, _opts_tolr(_HONOR_TOLR_TIGHT, 1))
    rd = _reldiff(off, on)
    assert rd > _HONOR_TOLR_FLOOR, (
        '-honorTolR 1 did not change the answer at TolR = %g. Either the seam '
        'is not wired (applyLadrunoConstants no longer sets mHonorTolRInME) or '
        'ModifiedEuler no longer reads it. Measured on this deck: 2.641e-04'
        % _HONOR_TOLR_TIGHT, rd, off, on)


def test_honor_tolr_warns_when_the_scheme_cannot_honour_it(capfd):
    """A flag that is accepted, stored and inert must SAY so.

    `mHonorTolRInME` is read at exactly one site, inside ModifiedEuler.  On
    IntScheme 45 the flag is therefore inert -- measured, reldiff exactly 0.0 --
    and RungeKutta45 already honours TolR unconditionally, so there is nothing
    for the flag to add there.  Accepting it silently would be a flag claiming
    to have done something it did not do: the defect ADR 86 was written about.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    capfd.readouterr()
    # IntScheme 45 (RungeKutta45), honorTolR 1.
    ops.nDMaterial('LadrunoSANISAND', 1, *_PARAMS, 45, 0, 1, 1.0e-7, 1.0e-6,
                   '-honorTolR', 1)
    err = capfd.readouterr().err
    assert 'NO EFFECT with IntScheme 45' in err, (
        '-honorTolR 1 on IntScheme 45 was accepted without warning. It is '
        'inert there (measured reldiff 0.0), so silence overstates what ran',
        err)
    assert 'RungeKutta45' in err, (
        'the warning does not name what the scheme does instead', err)

    # ... and it must NOT warn on the scheme that does honour it.
    ops.nDMaterial('LadrunoSANISAND', 2, *_PARAMS, 1, 0, 1, 1.0e-7, 1.0e-6,
                   '-honorTolR', 1)
    err2 = capfd.readouterr().err
    assert 'NO EFFECT' not in err2, (
        'the material warned that -honorTolR 1 is inert on IntScheme 1, which '
        'is the one scheme that DOES read the seam. A warning that fires '
        'everywhere is a warning nobody reads', err2)


# ===========================================================================
#  9. -Pmin is BEHAVIOURAL at low confinement (ADR 86 PR-3)
# ===========================================================================
#
#  THE OWED GATE.  Every test above this one pins `-Pmin` to vanilla's
#  1e-4*P_atm, deliberately, so that the p_residual A/B moves one variable
#  (ADR 86 risk 4).  That left the class's OWN default -- 1e-3*P_atm, TEN TIMES
#  vanilla's, inherited from NTUASand02 -- with no gate that TARGETS it.
#
#  AND THE OBVIOUS WAY TO SAY THAT IS WRONG, measured.  An earlier draft of this
#  comment claimed "no test could tell a build that honours -Pmin from one that
#  ignores it".  The PR-3 mutation campaign refuted it: with
#  `applyLadrunoConstants` mutated to ignore `mPminInput` and always resolve to
#  1e-3*P_atm, THREE tests fail -- this one,
#  `test_radial_ramp_with_pr0_never_yields` (which pins -Pmin only as a control
#  and is incidentally sensitive to it), and `test_print_states_what_it_ran`
#  (which catches the input/resolved-value mismatch in the record, not in the
#  behaviour).  So -Pmin was not invisible; what it lacked was a test that
#  measures it ON PURPOSE and identifies the MECHANISM, which is what the
#  deviator-wipe assertion below adds.  The incidental catches are fragile in a
#  way this one is not: change the ramp deck's constants and one of them stops
#  catching anything, silently.
#
#  WHAT IT DOES, MEASURED.  Confine-first deck, `-Presidual` PINNED at vanilla's
#  1.01 so only `-Pmin` moves, isochoric deviatoric leg, 10 + 40 steps:
#
#    e_conf    p at stage flip   -Pmin 0.0101      -Pmin 0.101       reldiff
#    3.0e-6    1.7175 kPa        p_end  9.4589     p_end  9.4589     0.0 exact
#    5.0e-7    0.2863 kPa        p_end  5.9791     p_end  5.9791     0.0 exact
#    4.0e-7    0.2290 kPa        p_end  5.5955     p_end 25.853      differs
#    2.0e-7    0.1145 kPa        p_end  4.2254     p_end 10.388      1.26
#
#  So it is a LOW-STRESS effect with a sharp boundary, and above it the two
#  values are bit-identical -- which is also why pinning `-Pmin` in the other
#  tests costs them nothing.
#
#  THE MECHANISM, and it is not the one this gate was expected to use.  The
#  handoff predicted the gate would assert on PR-2's throttled clamp diagnostic.
#  It cannot: on this deck ZERO clamp warnings are emitted on either leg.  What
#  actually happens is a WHOLE-TENSOR RESET at a site PR-2 did not instrument.
#  Traced per step at e_conf = 2e-7 (p at the flip 0.1145 kPa):
#
#    -Pmin 0.0101   step 1  p = 0.02389   |dev| = 0.5750   ... p rises, no reset
#                   min |dev| over the whole 40-step leg = 5.750e-01
#    -Pmin 0.101    step 1  p = -0.5647   |dev| = 9.5818   <- p goes NEGATIVE
#                   step 2  p =  0.1010   |dev| = 2.40e-17 <- RESET to m_Pmin*I
#
#  At step 2 the committed stress is exactly `m_Pmin * mI1`: p equals m_Pmin to
#  the last digit and the deviator is zero.  That is
#  `ManzariDafalias::explicit_integrator` (:1074/:1078) and/or `Stress_Correction`
#  (:2555/:2557/:2624), neither of which carries a warning -- PR-2's diagnostic
#  covers `ModifiedEuler` (:1378) and `RungeKutta45` (:1902) only, and BOTH of
#  those preserve the deviator (`GetDevPart(NextStress) + m_Pmin*mI1`) while the
#  uninstrumented ones zero it.  Recorded as owed work rather than fixed here;
#  see LEDGER_quirks.md.  (A draft of this note also cited `Stress_Correction`
#  :2608; that line is dead code, inside the `if (false)` opening at :2559.)
#
#  WHY THE DECK IS SOUND DESPITE THE REGIME.  ADR 86 risk 6 says no gate should
#  depend on an extreme-low-p leg COMPLETING.  This one does not depend on it --
#  it asserts it, and it holds: every one of the 50 steps converges on BOTH legs
#  (measured, 0 failed steps), with zero "Outside Bounding" events.  The
#  quantity asserted is a floor and a categorical wipe, never a calibrated
#  number: reldiff is 9.22 / 1.26 / 0.748 at n_dev 20 / 40 / 80, so the step
#  count moves it by an order of magnitude and the floor is set two orders below
#  the smallest.
#
#  COST: 3.2-3.5 s for the pair (measured), on a file that already costs ~10 s.
#  Deliberately NOT marked slow -- see the module docstring: no workflow passes
#  --runslow, so a slow marker deletes a test from CI.

_PMIN_VANILLA = _VANILLA_PMIN          # 1e-4 * P_atm = 0.0101 kPa
_PMIN_LADRUNO = 1.0e-3 * _P_ATM        # 1e-3 * P_atm = 0.101  kPa, the class default
_PMIN_E_CONF_LOW = 2.0e-7              # -> p = 0.1145 kPa at the stage flip
_PMIN_E_CONF_ORDINARY = _C_E_CONF      # -> p = 1.7175 kPa, the battery's own deck
_PMIN_N_DEV = 40
_PMIN_RELDIFF_FLOOR = 1.0e-2           # measured 0.748 at the LEAST sensitive n_dev
_PMIN_WIPE_TOL = 1.0e-10               # measured 2.4e-17 wiped vs 5.75e-01 not
_PMIN_NO_WIPE_FLOOR = 1.0e-2           # the vanilla leg's minimum is 5.75e-01


def _c_series_at(e_conf, n_dev, n_conf=_C_N_CONF, e_ax=_C_E_AX, lat=_C_LAT):
    """`_c_series` with the deck constants passed IN.

    `_c_series` binds n_conf/n_dev/e_conf/e_ax/lat as DEFAULT ARGUMENTS at def
    time, and its own docstring warns that monkeypatching the module constants
    desynchronises the Path series from the loop count -- the analysis then walks
    off the end of the series and unloads the model to a meaningless near-zero
    stress.  This is the local copy that docstring prescribes.
    """
    s_lat = [i / n_conf for i in range(n_conf + 1)]
    s_ax = list(s_lat)
    r_lat = lat * e_ax / e_conf
    r_ax = e_ax / e_conf
    for i in range(1, n_dev + 1):
        s_lat.append(1.0 - r_lat * i / n_dev)
        s_ax.append(1.0 + r_ax * i / n_dev)
    s_lat.append(s_lat[-1])            # the PathSeries hold point -- see _c_series
    s_ax.append(s_ax[-1])
    return s_lat, s_ax


def _drive_confined_at(e_conf, opts, n_dev=_PMIN_N_DEV, tag=1):
    """`_drive_confined` with the confinement magnitude and step count passed in.

    Returns (committed stress, [(p, |dev|) per deviatoric step]).  Every step is
    asserted to converge, on both legs -- this deck runs at p ~ 0.1 kPa and a
    silent non-convergence would look exactly like the effect being measured.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = _c_series_at(e_conf, n_dev)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -e_conf)
            if y == 1.:
                ops.sp(n, 2, -e_conf)
    ops.pattern('Plain', 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -e_conf)
    _analysis_confined()

    ops.updateMaterialStage('-material', tag, '-stage', 0)   # static mElastFlag
    for step in range(_C_N_CONF):
        assert ops.analyze(1) == 0, 'confinement step %d failed' % (step + 1)
    ops.updateMaterialStage('-material', tag, '-stage', 1)
    trace = []
    for step in range(n_dev):
        assert ops.analyze(1) == 0, 'deviatoric step %d failed' % (step + 1)
        sig = _stress()
        p = -(sig[0] + sig[1] + sig[2]) / 3.0
        dev = [-sig[0] - p, -sig[1] - p, -sig[2] - p, -sig[3], -sig[4], -sig[5]]
        trace.append((p, math.sqrt(sum(v * v for v in dev))))
    return _stress(), trace


def _opts_pmin(pmin):
    """ONLY -Pmin free: p_residual pinned to vanilla's, the seam flag off."""
    return ('-Presidual', _VANILLA_PR, '-Pmin', pmin, '-honorTolR', 0)


def test_pmin_is_inert_at_ordinary_confinement():
    """Above the floor, the class default and vanilla's are BIT-IDENTICAL.

    This is the half that licenses every other test in this file pinning
    `-Pmin`: on the battery's own confine-first deck (p = 1.7175 kPa at the
    stage flip, p_end 9.46) a tenfold change in `-Pmin` moves nothing at all.

    It is a property of THIS DECK'S STRAIN PATH, not a general claim -- what
    falsifies it is a path whose p dips below 0.101 kPa, which is exactly what
    the companion test below constructs.  Stated as a pair on purpose: an
    inertness measurement is only as general as the deck it was taken on.
    """
    a, _ = _drive_confined_at(_PMIN_E_CONF_ORDINARY, _opts_pmin(_PMIN_VANILLA))
    b, _ = _drive_confined_at(_PMIN_E_CONF_ORDINARY, _opts_pmin(_PMIN_LADRUNO))
    assert _reldiff(a, b) == 0.0, (
        'a 10x change in -Pmin moved the answer on a deck that never goes near '
        'the floor (p >= 1.72 kPa). Either the floor is being read somewhere it '
        'should not be, or this deck no longer stays above it', a, b)


def test_pmin_is_behavioural_at_low_confinement():
    """The owed gate: `-Pmin` is not inert, and the class default is not free.

    Only `-Pmin` moves -- p_residual is pinned at vanilla's 1.01 and the
    honorTolR seam is off.  Two independent assertions, because a difference
    alone does not identify a cause:

      1. the committed answer differs by a factor, not a rounding;
      2. the LadrunoSANISAND-default leg has its DEVIATOR WIPED -- one committed
         step at exactly `p = m_Pmin`, |dev| = 0 -- and the vanilla-default leg
         never comes within two orders of that.

    Assertion 2 is what makes this mutation-proof.  A build that ignores
    `mPminInput` and always resolves to 1e-3*P_atm runs BOTH legs at 0.101, so
    both wipe and the 'vanilla leg never wipes' assertion fails; a build that
    always resolves to 1e-4*P_atm wipes NEITHER and assertion 1 collapses to 0.
    """
    a, tr_a = _drive_confined_at(_PMIN_E_CONF_LOW, _opts_pmin(_PMIN_VANILLA))
    b, tr_b = _drive_confined_at(_PMIN_E_CONF_LOW, _opts_pmin(_PMIN_LADRUNO))

    rd = _reldiff(a, b)
    assert rd > _PMIN_RELDIFF_FLOOR, (
        '-Pmin moved from %g to %g and the answer did not change. Every other '
        'test in this file pins -Pmin, so this is the only one that varies it '
        'on purpose. Measured: 1.26 at n_dev=40'
        % (_PMIN_VANILLA, _PMIN_LADRUNO), rd, a, b)

    min_dev_a = min(d for _, d in tr_a)
    wiped_b = [(i + 1, p, d) for i, (p, d) in enumerate(tr_b)
               if d < _PMIN_WIPE_TOL]

    assert min_dev_a > _PMIN_NO_WIPE_FLOOR, (
        'the -Pmin = %g leg had its deviator wiped too, so the two legs are not '
        'being told apart by the floor. Measured minimum |dev| on this leg: '
        '5.75e-01' % _PMIN_VANILLA, min_dev_a)

    assert wiped_b, (
        'the -Pmin = %g leg never reset. The mechanism this gate identifies is a '
        'whole-tensor reset to m_Pmin*I at deviatoric step 2 (measured |dev| = '
        '2.4e-17); without it the reldiff above could be any low-p wobble'
        % _PMIN_LADRUNO, [d for _, d in tr_b])

    # ... and the reset lands exactly ON the floor, which is what identifies it
    # as `NextStress = m_Pmin * mI1` rather than as ordinary softening.
    for step, p, _d in wiped_b:
        assert abs(p - _PMIN_LADRUNO) < 1.0e-9 * _PMIN_LADRUNO, (
            'deviator wiped at step %d but p = %g is not m_Pmin = %g, so this is '
            'not the m_Pmin*I reset this gate claims to have found'
            % (step, p, _PMIN_LADRUNO))


# ===========================================================================
#  10. LadrunoSANISAND3D::getCopy(void) -- the last uncovered override
# ===========================================================================
#
#  WHY IT WAS UNCOVERED, AND WHY THAT IS NOT OBVIOUS.  A brick asks its material
#  for `getCopy("ThreeDimensional")`, never `getCopy()`.  So the no-argument
#  form -- `LadrunoSANISAND3D::getCopy()`, which clones via
#  `new LadrunoSANISAND3D(); *clone = *this;` -- has no caller on any ordinary
#  deck, and PR-1's mutation campaign proved it: the whole battery stayed green
#  with that function sabotaged.  The ADR log records it as knowingly open.
#
#  THE ROUTE, established by reading the callers rather than assumed.  Only one
#  reachable path calls `getCopy()` on an object that is ALREADY a
#  dimension-specific SANISAND:
#
#    nDMaterial LadrunoSANISAND 1 ...      -> the dimensionless PROTOTYPE
#                                             (classTag ND_TAG_LadrunoSANISAND)
#    nDMaterial InitStrain 2 1 0.0         -> InitStrainNDMaterial(2, *proto, 0.0)
#                                             ctor does proto.getCopy("ThreeDimensional")
#                                             (InitStrainNDMaterial.cpp:153)
#                                             => it now HOLDS a LadrunoSANISAND3D
#    element stdBrick ... 2                -> Brick asks the WRAPPER for
#                                             getCopy("ThreeDimensional")
#                                             -> InitStrainNDMaterial::getCopy(const char*)
#                                                returns getCopy()          (:293)
#                                             -> InitStrainNDMaterial::getCopy(void)
#                                                calls theMaterial->getCopy() (:285)
#                                             => LadrunoSANISAND3D::getCopy(void),
#                                                once per Gauss point, 8 times.
#
#  `nDMaterial InitStress` does NOT work for this: its 3-argument constructor
#  calls `material.getCopy()` on the PROTOTYPE, and the prototype's
#  `getCopy(void)` is the base's "subclass responsibility" + `exit(-1)` -- it
#  kills the process. That is upstream behaviour and true of vanilla
#  ManzariDafalias too; it is recorded here so the next person does not spend the
#  same hour on it.
#
#  eps0 = 0.0 makes the wrapper a pure pass-through, which is asserted rather
#  than assumed: measured BIT-IDENTICAL (reldiff exactly 0.0) at both p_r values.
#
#  MEASURED on the confine-first deck:
#      unwrapped p_r=0     vs wrapped p_r=0      reldiff 0.0000e+00
#      unwrapped p_r=1.01  vs wrapped p_r=1.01   reldiff 0.0000e+00
#      wrapped   p_r=0     vs wrapped p_r=1.01   reldiff 8.6042e-02

_GETCOPY_NONVACUITY_FLOOR = 1.0e-3     # measured 8.6042e-02


def _drive_confined_wrapped(opts, tag=1, wrap_tag=2):
    """The confine-first deck with an eps0 = 0 InitStrain wrapper in between.

    The wrapper is what forces the element through
    `LadrunoSANISAND3D::getCopy(void)` -- see the route note above.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.nDMaterial('InitStrain', wrap_tag, tag, 0.0)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, wrap_tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = _c_series()
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -_C_E_CONF)
            if y == 1.:
                ops.sp(n, 2, -_C_E_CONF)
    ops.pattern('Plain', 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -_C_E_CONF)
    _analysis_confined()
    # The stage parameter is keyed on the INNER material's tag: the copies the
    # wrapper holds carry `tag`, not `wrap_tag`, because
    # LadrunoSANISAND::getCopy(const char*) passes `this->getTag()` through.
    # InitStrainNDMaterial::setParameter forwards anything it does not recognise
    # to the wrapped material, and "updateMaterialStage" is not one of its own
    # (its parameters are all eps0*), so the request reaches the SANISAND.
    _confine_leg(tag)
    _deviator_leg(tag)
    return _stress()


def test_getcopy_void_carries_the_settings():
    """`LadrunoSANISAND3D::getCopy(void)` must clone the Ladruno constants.

    Two assertions, and the second is the one that cannot be faked:

      1. PASS-THROUGH. An eps0 = 0 InitStrain wrapper must not change the
         answer, at either p_residual. Asserted BIT-IDENTICAL, so a
         `getCopy(void)` that dropped or reset any part of the state -- not just
         the two constants -- shows up here.
      2. NON-VACUITY THROUGH THE WRAPPER. With the wrapper in place, the p_r = 0
         and p_r = 1.01 decks must still differ. Every Gauss-point material on
         this deck was built by `LadrunoSANISAND3D::getCopy(void)`, so a version
         of it that returned a vanilla ManzariDafalias3D, or that let
         p_residual fall back to 1.0e-2*P_atm, would make BOTH legs run the same
         soil and collapse this to ~0.

    The UNWRAPPED decks are the control: they never call `getCopy(void)` at all,
    so any wrapped-vs-unwrapped difference is attributable to that one function.
    Assertion 1 is therefore the primary kill, and the p_r = 0 leg is the sharp
    one -- a clone that let p_residual fall back to vanilla's 1.0e-2*P_atm would
    make the wrapped p_r = 0 deck run 1.01 while its unwrapped control runs 0,
    a measured 8.6e-2 apart.

    Assertion 2 is the non-vacuity companion, in the same spirit as G1's: it
    proves the wrapped deck is still a live experiment that RESPONDS to
    p_residual, so assertion 1's 'equal' is equality between two things that
    could have differed. G1 taught this the hard way -- under one mutation it
    passed by comparing vanilla with vanilla.
    """
    opts_pr0 = _OPTS_PR0
    opts_van = _OPTS_VANILLA

    plain_pr0 = _drive_confined('LadrunoSANISAND', 1, opts_pr0)
    plain_van = _drive_confined('LadrunoSANISAND', 1, opts_van)
    wrap_pr0 = _drive_confined_wrapped(opts_pr0)
    wrap_van = _drive_confined_wrapped(opts_van)

    assert _reldiff(plain_pr0, wrap_pr0) == 0.0, (
        'an eps0 = 0 InitStrain wrapper changed the p_r = 0 answer. The wrapper '
        'is a pass-through by construction, so this is the clone made by '
        'LadrunoSANISAND3D::getCopy(void) differing from the one made directly '
        'by getCopy("ThreeDimensional")', plain_pr0, wrap_pr0)
    assert _reldiff(plain_van, wrap_van) == 0.0, (
        'an eps0 = 0 InitStrain wrapper changed the p_r = 1.01 answer',
        plain_van, wrap_van)

    rd = _reldiff(wrap_pr0, wrap_van)
    assert rd > _GETCOPY_NONVACUITY_FLOOR, (
        'with every Gauss point built by LadrunoSANISAND3D::getCopy(void), the '
        'p_r = 0 and p_r = 1.01 decks became indistinguishable. That function '
        'is not carrying p_residual. Measured on this deck: 8.6042e-02', rd,
        wrap_pr0, wrap_van)


# ===========================================================================
#  11. LadrunoSANISANDPlaneStrain::getCopy(void) -- the plane-strain twin
# ===========================================================================
#
#  WHY THE 3D TRICK ABOVE DOES NOT PORT.  Section 10's route works only
#  because `InitStrainNDMaterial` special-cases the LITERAL STRING
#  "ThreeDimensional": its constructor always resolves `theMaterial` via
#  `material.getCopy("ThreeDimensional")` (InitStrainNDMaterial.cpp:135,153),
#  and `InitStrainNDMaterial::getCopy(const char*)` forwards that SAME string
#  straight to its own void `getCopy()` (:292-293) -- i.e. to
#  `theMaterial->getCopy()`.  `quad ... PlaneStrain` never asks a wrapper for
#  "ThreeDimensional" (FourNodeQuad.cpp:334 passes the literal type string the
#  element command was given), and InitStrain's OWN "PlaneStrain" branch (the
#  `ncomp` switch, InitStrainNDMaterial.cpp:302-322) calls
#  `theMaterial->getCopy("PlaneStrain")` -- the type-string constructor, NOT
#  the void clone.  So wrapping a SANISAND in InitStrain and driving it with a
#  `quad` element never reaches `LadrunoSANISANDPlaneStrain::getCopy(void)`;
#  it only re-exercises `LadrunoSANISAND::getCopy(const char*)`'s PlaneStrain
#  branch, already covered directly by `_build_ps`/`_build_ps_confined`.
#  `StagedStrainNDMaterial` (StagedStrainNDMaterial.cpp:293-309, the pattern
#  WP-C's wrapper is modelled on) has the identical shape: no special case for
#  any string but its own construction-time template.  Neither wrapper, nor
#  any nesting of the two, can substitute.
#
#  THE ROUTE THAT DOES WORK: `FluidSolidPorousMaterial`
#  (SRC/material/nD/soil/FluidSolidPorousMaterial.cpp), the pore-pressure
#  wrapper `test_fspm_over_manzari_family.py` already gates for the "subclass
#  responsibility" defect (ADR 86 follow-up).  Its constructor asks the soil
#  for a DIMENSION-MATCHED copy up front
#  (`soilMat.getCopy(nd == 3 ? "ThreeDimensional" : "PlaneStrain")`, :125), so
#  `theSoilMaterial` is ALREADY a `LadrunoSANISANDPlaneStrain` once `nd == 2`.
#  Its OWN `getCopy(const char*)` (:415-425) does NOT forward the type string
#  to the soil at all -- for both "PlaneStrain" and "ThreeDimensional" it just
#  runs its own copy constructor (`FluidSolidPorousMaterial(const
#  FluidSolidPorousMaterial&)`, :152-162), which clones the soil with the
#  dimension-FREE, no-argument `theSoilMaterial->getCopy()` (:156).  So:
#
#    nDMaterial LadrunoSANISAND 1 ...            -> the dimensionless PROTOTYPE
#    nDMaterial FluidSolidPorous 2 2 1 K atm      -> FSPM ctor asks 1 for
#                                                     getCopy("PlaneStrain")
#                                                     (the type-string path,
#                                                     already covered) ->
#                                                     theSoilMaterial is now a
#                                                     LadrunoSANISANDPlaneStrain
#    element quad 1 ... PlaneStrain 2             -> quad asks the WRAPPER
#                                                     (tag 2) for
#                                                     getCopy("PlaneStrain"),
#                                                     once per Gauss point ->
#                                                     FSPM's copy ctor ->
#                                                     theSoilMaterial->getCopy()
#                                                  => LadrunoSANISANDPlaneStrain
#                                                     ::getCopy(void), 4 times
#                                                     (once per Gauss point).
#
#  `loadStagex` (the wrapper's OWN undrained switch) is left at its
#  construction default of 0 for the whole deck below -- `updateMaterialStage`
#  is only ever aimed at the SOIL's tag, never the wrapper's -- so
#  `FluidSolidPorousMaterial::getStress()` (:296) takes its `loadStage == 0`
#  branch and returns `theSoilMaterial->getStress()` completely unmodified: no
#  pore-pressure term is ever added.  That is what licenses comparing the
#  wrapped deck against the bare `_drive_ps_confined` deck bit-for-bit below.
# ---------------------------------------------------------------------------

_FSPM_KCOMB = 1.0e6     # value is irrelevant here: loadStage never leaves 0


def _build_ps_fspm(tag, wrap_tag, opts=(), mirror=False):
    """`_build_ps_confined`'s twin, with the element pointed at a
    `FluidSolidPorous` wrapper (tag `wrap_tag`) instead of directly at the
    soil tag."""
    sgn = -1.0 if mirror else 1.0
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    for j, (x, y) in enumerate(_XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.nDMaterial('FluidSolidPorous', wrap_tag, 2, tag, _FSPM_KCOMB, _P_ATM)
    ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', wrap_tag)
    for j, (x, y) in enumerate(_XY):         # rollers on the two negative edges
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    s_lat, s_ax = _c_series(lat=_C_LAT_PS)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for j, (x, y) in enumerate(_XY):
        if x == 1.:
            ops.sp(j + 1, 1, sgn * -_C_E_CONF)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        if y == 1.:
            ops.sp(j + 1, 2, sgn * -_C_E_CONF)
    _analysis_confined()


def _drive_ps_fspm(tag, wrap_tag, opts=(), mirror=False):
    """wipe -> build (through FluidSolidPorous) -> confine -> shear -> stress.

    `updateMaterialStage` below is aimed at the SOIL's own tag, never the
    wrapper's -- exactly as the 3D InitStrain route requires.  FSPM forwards
    any "-material $tag" request whose tag is not its own straight to
    `theSoilMaterial->setParameter` (FluidSolidPorousMaterial.cpp:318-333),
    and every clone the wrapper holds carries the SOIL's tag (both
    `getCopy(type)` and `getCopy(void)` preserve it), so the request reaches
    the SANISAND clone, not the inert wrapper.
    """
    _build_ps_fspm(tag, wrap_tag, opts, mirror)
    _confine_leg(tag)
    _deviator_leg(tag)
    return _stress()


def test_getcopy_void_carries_the_settings_planestrain():
    """`LadrunoSANISANDPlaneStrain::getCopy(void)` must clone the Ladruno
    constants -- the plane-strain twin of `test_getcopy_void_carries_the_
    settings` (section 10), routed through `FluidSolidPorous` rather than
    `InitStrain` because the InitStrain route does not reach this override at
    all in 2D (see the block comment above).

    Same two assertions as the 3D test, and the second is the one that cannot
    be faked:

      1. PASS-THROUGH. The `FluidSolidPorous`-wrapped deck must be BIT
         IDENTICAL to the bare `_drive_ps_confined` deck, at both p_r = 0 and
         p_r = 1.01 kPa -- the wrapper's own undrained coupling never fires
         (`loadStage` stays at its construction default of 0 for the whole
         run), so any difference is attributable to
         `LadrunoSANISANDPlaneStrain::getCopy(void)` alone.
      2. NON-VACUITY. Through the wrapper, p_r = 0 and p_r = 1.01 must still
         differ materially. Every Gauss-point material on the wrapped deck was
         built by `LadrunoSANISANDPlaneStrain::getCopy(void)`
         (`FluidSolidPorousMaterial`'s copy constructor clones the wrapped
         soil with its dimension-free `getCopy()`), so a version of that
         override that dropped `mPresidualInput`/`mPminInput`/`mHonorTolR` --
         or let them fall back to compiler-default-initialized garbage, or to
         the base class's own values -- would make both legs collapse to the
         same answer, or diverge from the unwrapped control in assertion 1.

    The UNWRAPPED deck (`_drive_ps_confined`) is the control: it goes through
    `LadrunoSANISAND::getCopy(const char*)` directly (already covered by
    `test_planestrain_lane_carries_the_settings`) and never calls
    `getCopy(void)` at all, so a wrapped-vs-unwrapped difference is
    attributable to that one function.
    """
    opts_pr0 = _OPTS_PR0
    opts_van = _OPTS_VANILLA

    plain_pr0 = _drive_ps_confined(1, opts_pr0)
    plain_van = _drive_ps_confined(1, opts_van)
    wrap_pr0 = _drive_ps_fspm(1, 2, opts_pr0)
    wrap_van = _drive_ps_fspm(1, 2, opts_van)

    assert _reldiff(plain_pr0, wrap_pr0) == 0.0, (
        'a FluidSolidPorous wrapper with loadStage == 0 changed the p_r = 0 '
        'answer. The wrapper is a pass-through by construction at that stage, '
        'so this is the clone made by LadrunoSANISANDPlaneStrain::getCopy(void) '
        'differing from the one made directly by getCopy("PlaneStrain")',
        plain_pr0, wrap_pr0)
    assert _reldiff(plain_van, wrap_van) == 0.0, (
        'a FluidSolidPorous wrapper with loadStage == 0 changed the p_r = 1.01 '
        'answer', plain_van, wrap_van)

    rd = _reldiff(wrap_pr0, wrap_van)
    assert rd > _GETCOPY_NONVACUITY_FLOOR, (
        'with every Gauss point built by LadrunoSANISANDPlaneStrain::getCopy'
        '(void), the p_r = 0 and p_r = 1.01 decks became indistinguishable. '
        'That function is not carrying p_residual.', rd, wrap_pr0, wrap_van)
