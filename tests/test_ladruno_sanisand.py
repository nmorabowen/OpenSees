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
     It carries a SENSITIVITY companion
     (``test_vanilla_equivalence_is_not_vacuous``): the same comparison with
     ``-Presidual 0`` must give a MATERIALLY different answer.
     Without it, a parser that silently dropped the flags (exactly the defect
     ADR 86 sec.7.1 documents in the sibling ``OPS_SAniSandMSMaterial``) would
     make this test compare vanilla with vanilla and pass while proving
     nothing.
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
equilibrium (see the note on test 5).

SLOW TIER: ``test_presidual_is_the_low_p_defect`` is marked
``@pytest.mark.slow``; measured wall time 21.4 s on the dev box (Windows,
2026-08-26) -- two 1200-step drained triaxial legs at p ~ 1 kPa.  Everything
else in this file runs in about 0.2 s total.  22 s is not minutes-scale, so
the marker is NOT primarily about cost: it is there because ADR 86 risk 6 says
no gate should depend on an extreme-low-p leg completing, and this test does
(its p_r = 0 leg fails outright at 400 steps, where vanilla survives).  If the
project decides that gate is worth paying for on every push, dropping the
marker costs 22 s -- that is a scoping call, not a technical obstacle.
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
#  with the elastoplastic integrator running.  The plastic leg is genuinely
#  plastic: it ends at p ~ 13.7 kPa with an active yield surface, and the
#  p_residual sensitivity below (5.6 %) exists ONLY on plastic paths --
#  GetElasticModuli does not read m_Presidual at all, so an elastic-only path
#  would show zero sensitivity and prove nothing.
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


# ===========================================================================
#  1. vanilla equivalence -- the load-bearing gate
# ===========================================================================

_EQ_TOL = 1.0e-12
_SENSITIVITY_FLOOR = 1.0e-3    # measured 5.6e-2; see below


def test_vanilla_equivalence():
    """LadrunoSANISAND handed vanilla's two constants IS vanilla, to 1e-12.

    Measured on the dev box: identically 0.0 (bit-for-bit), well inside 1e-12.

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


def test_vanilla_equivalence_is_not_vacuous():
    """The companion to test 1: the SAME path must be sensitive to p_residual.

    If it were not -- e.g. because the parser silently dropped `-Presidual`
    (ADR 86 sec.7.1 documents exactly that bug in the sibling
    `OPS_SAniSandMSMaterial`), or because the strain path never left the
    elastic range -- then test_vanilla_equivalence would be comparing vanilla
    with vanilla and would pass while proving nothing.

    Only `-Presidual` moves between the two legs: `-Pmin` is pinned to
    vanilla's 1e-4*P_atm in both (ADR 86 risk 4).  Measured: 5.58e-2, i.e.
    5.6 % -- ten orders of magnitude above the 1e-12 equivalence tolerance.
    """
    vanilla_consts = _drive('LadrunoSANISAND', 5, _OPTS_VANILLA)
    no_cohesion = _drive('LadrunoSANISAND', 6, _OPTS_PR0)

    rel = _reldiff(vanilla_consts, no_cohesion)
    assert rel >= _SENSITIVITY_FLOOR, (
        'switching p_residual from 1.01 kPa to 0 changed the committed stress '
        f'by only {rel:.3e} -- this strain path is NOT sensitive to '
        'p_residual, so test_vanilla_equivalence proves nothing on it',
        vanilla_consts, no_cohesion)


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
    # PR-1 ships the D_factor sigmoid unchanged and kPa-dimensional (D5a/D5b).
    assert 'D5a/D5b' in err, (
        'the echo no longer declares the deferred D_factor decision', err)

    # --- not latched: a second construction must echo again ------------------
    ops.nDMaterial('LadrunoSANISAND', 2, *_PARAMS)
    ops.nDMaterial('LadrunoSANISAND', 3, *_PARAMS)
    err2 = capfd.readouterr().err
    assert err2.count('p_residual = 0') == 2, (
        'the echo is latched (once per process).  A latched message is '
        'observable only by whichever test in a session runs first, which is '
        'exactly how the p_residual defect stayed invisible', err2)

    # --- -honorTolR is a PR-1 refusal, not a no-op ---------------------------
    # A flag that claims to have done something it did not do is the class of
    # defect this ADR exists to fix, so PR-1 rejects any value but 0.
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 4, *_PARAMS, '-honorTolR', 1)
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
_TCL_DECK = """\
model basic -ndm 3 -ndf 3
nDMaterial LadrunoSANISAND 1 {params} 1 0 1 1e-07 1e-07 \\
    -Presidual {pr} -Pmin {pmin} -honorTolR 0
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
    override is worthless): the strain path ends at p ~ 13.7 kPa, and the
    identical path run with p_r = 1.01 instead of 0 gives a committed stress
    5.4e-2 (5.4 %) away -- see the `_reldiff(first, reverted)` assertion at the
    bottom, which measures that gap in this very test rather than trusting it.
    The pass/fail margin is therefore ~5e-2 against a 1e-12 tolerance.

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
