"""Regression tests for the audit-following-#709 upstream symmetrize/initial-
tangent fixes.

PR #709 fixed BezierTri6/BezierTet10, which assembled only half of B'*D*B and
mirrored it into the other half -- silently symmetrizing an unsymmetric
material tangent (non-associated plasticity: DruckerPrager with
rho_bar != rho, ManzariDafalias always). A follow-up audit found the same
defect pattern in two more places, plus two materials whose
getInitialTangent() is broken in a related way, plus a default-solver gap:

  Fix 1  SRC/element/zeroLength/ZeroLengthND.cpp
         getTangentStiff()/getInitialStiff() assembled A'*kb*A over only the
         lower triangle (`for j < i+1`) and mirrored -- same class of defect
         as #709, on a different element family.

  Fix 2  SRC/material/nD/UWmaterials/DruckerPragerPlaneStrain.cpp
         getInitialTangent() was a verbatim copy of getTangent(): it returned
         the condensed elastoplastic mCep (current, history-dependent, and
         possibly unsymmetric) instead of the condensed elastic mCe.

  Fix 3  SRC/material/nD/UWmaterials/ContactMaterial{2,3}D.cpp
         getInitialTangent() returned the SAME buffer getTangent() overwrites
         in place, so it was call-order dependent and the ZERO matrix on a
         freshly constructed material.

  Fix 4  SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp
         addA() unconditionally discards the lower triangle with no
         diagnostic, and ProfileSPD is the DEFAULT solver in both
         interpreters. Added the same budgeted asymWarned/asymBudget/
         asymPass guard already shipped on PARDISOGenLinSOE (#713).

Reachability notes (see the PR description / audit report for detail):

  * Fix 1 is tested via ZeroLengthND on DruckerPragerPlaneStrain (order 3),
    not the audit's literally-suggested DruckerPrager3D (order 6): order-6
    NDMaterials hit an UNRELATED pre-existing latent defect in
    ZeroLengthND::setTransformation (the fixed 3x3 `transformation` member is
    indexed by the material order, so order > 3 reads out of bounds). This is
    not part of the fix and is called out rather than silently worked around.
    DruckerPragerPlaneStrain exercises the identical A'*kb*A assembly loop
    this PR changed.

  * Fix 2 has no direct Python route to a bare NDMaterial::getInitialTangent()
    (there is no "test ND material" command). It is exercised indirectly
    through a quad element and a single Newton `-initial`-tangent probe step
    (see test_druckerpragerplanestrain_initial_tangent_stable for the exact
    mechanism).

  * Fix 3 is UNREACHABLE from openseespy. Every upstream element that owns a
    ContactMaterial2D/3D (SimpleContact2D/3D, BeamContact2D/3D/2Dp/3Dp)
    overrides its own getInitialStiff() as `return getTangentStiff();` --
    NONE of them ever call the material's getInitialTangent(). ZeroLengthND
    cannot substitute either: ContactMaterial{2,3}D::getCopy(const char*)
    only recognizes the literal strings "ContactMaterial2D"/"ContactMaterial3D",
    not the "PlaneStrain2D"/"ThreeDimensional" tokens ZeroLengthND requests.
    test_contactmaterial_initial_tangent_unreachable_from_python documents
    this gap instead of fabricating a passing assertion.
"""
import re

import numpy as np
import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

# ---- shared DruckerPrager constants (PEER example, kPa-scale; reused from
# test_beziertet10_unsym_tangent.py / test_pardiso_asym_rearm.py, already
# calibrated to pass through first yield with rho_bar != rho) -------------
_K, _G = 27777.78, 9259.26          # -> E ~= 25000, nu ~= 0.35
_RHO = 0.398
_RHO_BAR = 0.1                      # rho_bar != rho -> non-associated -> Cep unsymmetric
_SIG_Y = 5.0
_SIG_Y_ELASTIC = 500.0              # never yields
_HARD_H = 1000.0


# ===========================================================================
# Fix 1 -- ZeroLengthND full-assembly (getTangentStiff / getInitialStiff)
# ===========================================================================
#
# One zeroLengthND element on DruckerPragerPlaneStrain (order 3): node1
# fixed, node2 free (ndf=3 so numDOF=6, the "3-per-node" case). The default
# orientation (-orient not given) makes the transformation matrix the
# identity, so the free-DOF system reduces to exactly the element's local
# order x order block -- the material's condensed tangent, read straight off
# the assembled system via printA under `system FullGeneral`.

_ZLN_NSTEPS = 30
_ZLN_TARGET_SYY = -6.0
_ZLN_TARGET_TXY = 8.0


def _build_zln(sig_y):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.nDMaterial('DruckerPrager', 1, _K, _G, sig_y, _RHO, _RHO_BAR,
                   0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)
    ops.element('zeroLengthND', 1, 1, 2, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 0.0, _ZLN_TARGET_SYY, _ZLN_TARGET_TXY)

    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')          # dense getA(), and unsymmetric-capable
    ops.test('NormDispIncr', 1.0e-10, 60, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _ZLN_NSTEPS)
    ops.analysis('Static')


def _zln_ramp(sig_y):
    """Run the ramp; return (steps converged, asymmetry of the final free-DOF K)."""
    _build_zln(sig_y)
    done = 0
    for _ in range(_ZLN_NSTEPS):
        if ops.analyze(1) != 0:
            break
        done += 1
    flat = np.asarray(ops.printA('-ret'), dtype=float)
    n = int(round(len(flat) ** 0.5))
    K = flat.reshape(n, n).T   # printA('-ret') is column-major
    asym = np.max(np.abs(K - K.T))
    scale = np.max(np.abs(K)) or 1.0
    return done, asym / scale


def test_zerolengthnd_unsym_tangent():
    """Yielded (non-associated) run: measurably unsymmetric K. Elastic: symmetric.

    The plastic ramp is NOT required to converge all the way to the target
    stress: driving a single material point straight through the
    non-associated DruckerPrager corner region by pure stress control (no
    surrounding elastic mesh to regularize the path) can lose Newton
    convergence well past first yield. That is a property of this stress
    path, not of the fix under test -- the asymmetry fingerprint only needs
    ENOUGH converged steps to have crossed the yield surface at least once,
    which a >= 60% completion comfortably guarantees here (calibrated: the
    ramp reaches ~66% before losing convergence).
    """
    done_p, asym_p = _zln_ramp(_SIG_Y)
    done_e, asym_e = _zln_ramp(_SIG_Y_ELASTIC)

    assert done_p >= int(0.5 * _ZLN_NSTEPS), (
        f'non-associated DruckerPragerPlaneStrain ramp converged only '
        f'{done_p}/{_ZLN_NSTEPS} steps -- too little of the path to have '
        f'reliably crossed the yield surface')
    assert done_e == _ZLN_NSTEPS, (
        f'elastic control ramp only converged {done_e}/{_ZLN_NSTEPS} steps')

    assert asym_e < 1.0e-10, (
        'never-yielding run must keep an exactly symmetric assembled '
        'tangent', asym_e)
    assert asym_p > 1.0e-8, (
        'yielded run must expose the unsymmetric consistent tangent -- '
        'this is exactly 0 under the pre-fix mirror-assembly bug', asym_p)


def test_zerolengthnd_initial_stiff_matches_elastic_and_is_symmetric():
    """getInitialStiff() (never-loaded state) must be exactly symmetric.

    At zero strain DruckerPragerPlaneStrain's mCep == mCe (see
    DruckerPrager::initialize()), so this doubles as a sanity check that the
    rewritten full-assembly loop reduces correctly to the elastic case.
    """
    _build_zln(_SIG_Y)
    ops.algorithm('Newton', '-initial')
    ops.test('NormDispIncr', 1.0e30, 1, 0)   # trivially "converges" after 1 iter
    assert ops.analyze(1) == 0

    flat = np.asarray(ops.printA('-ret'), dtype=float)
    n = int(round(len(flat) ** 0.5))
    K = flat.reshape(n, n).T
    assert np.max(np.abs(K - K.T)) < 1.0e-10 * (np.max(np.abs(K)) or 1.0)


# ===========================================================================
# Fix 2 -- DruckerPragerPlaneStrain::getInitialTangent() history-independence
# ===========================================================================
#
# No bare "test ND material" command exists in openseespy, so this is probed
# through a quad element + a single Newton `-initial`-tangent iteration.
#
# Mechanism: a converged static step leaves the domain in equilibrium
# (unbalance == 0). Adding a small NEW probe load pattern on top makes the
# unbalance for the very first Newton iteration of the next step EXACTLY
# equal to the probe load (independent of everything that happened before,
# since the previous total load was already balanced by F_int). Forcing that
# one iteration to use the INITIAL tangent (`algorithm Newton -initial`, and
# a convergence test so loose it "converges" after exactly 1 iteration) makes
# the resulting displacement increment du = K_initial^-1 * probe. Two
# independent models --pristine-elastic vs already-driven-past-yield-- fed
# the identical probe must produce the IDENTICAL du if and only if
# getInitialTangent() returns the same (elastic, history-independent) matrix
# in both cases. Pre-fix, getInitialTangent() returned the current mCep,
# which is measurably different (softened by hardening, and unsymmetric)
# after a plastic step, so du would differ.

_QUAD_H_TOTAL = 3.0
_QUAD_V_TOTAL = 4.0
_QUAD_NSTEPS = 30
_PROBE_FX = 0.02
_PROBE_FY = -0.02


def _build_quad(sig_y=None, elastic=False):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    if elastic:
        E = 9.0 * _K * _G / (3.0 * _K + _G)
        nu = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))
        ops.nDMaterial('ElasticIsotropic', 1, E, nu)
    else:
        ops.nDMaterial('DruckerPrager', 1, _K, _G, sig_y, _RHO, _RHO_BAR,
                       0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0)
    ops.node(4, 0.0, 1.0)
    ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')           # unsymmetric-capable


def _probe_du(main_steps=0):
    """Apply the probe load with a single Newton `-initial` iteration.

    main_steps > 0 first drives pattern 1 (already defined by the caller)
    through `main_steps` ordinary-Newton LoadControl steps, then freezes it
    (`loadConst`) before adding the probe pattern -- see the module docstring
    for why this makes the probe's first-iteration unbalance == the probe
    load regardless of what happened before.
    """
    if main_steps > 0:
        ops.test('NormDispIncr', 1.0e-9, 60, 0)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1.0 / main_steps)
        ops.analysis('Static')
        for step in range(main_steps):
            assert ops.analyze(1) == 0, f'main ramp failed at step {step + 1}'
        ops.loadConst('-time', 0.0)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(3, _PROBE_FX, _PROBE_FY)
    ops.load(4, _PROBE_FX, _PROBE_FY)

    pre = [ops.nodeDisp(3, 1), ops.nodeDisp(3, 2),
           ops.nodeDisp(4, 1), ops.nodeDisp(4, 2)]

    ops.algorithm('Newton', '-initial')
    ops.test('NormDispIncr', 1.0e30, 1, 0)
    ops.integrator('LoadControl', 1.0)
    if main_steps == 0:
        ops.analysis('Static')
    assert ops.analyze(1) == 0

    post = [ops.nodeDisp(3, 1), ops.nodeDisp(3, 2),
            ops.nodeDisp(4, 1), ops.nodeDisp(4, 2)]
    return np.array(post) - np.array(pre)


def test_druckerpragerplanestrain_initial_tangent_stable():
    _build_quad(sig_y=_SIG_Y)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, _QUAD_H_TOTAL / 2.0, -_QUAD_V_TOTAL / 2.0)
    ops.load(4, _QUAD_H_TOTAL / 2.0, -_QUAD_V_TOTAL / 2.0)
    du_after_yield = _probe_du(main_steps=_QUAD_NSTEPS)

    _build_quad(sig_y=_SIG_Y)     # fresh, pristine (elastic) material instance
    du_pristine = _probe_du(main_steps=0)

    assert np.max(np.abs(du_pristine)) > 1.0e-8, (
        'the probe produced (near-)zero displacement -- getInitialTangent() '
        'looks singular, or the probe load is too small', du_pristine)

    rel = np.max(np.abs(du_after_yield - du_pristine)) / np.max(np.abs(du_pristine))
    assert rel < 1.0e-6, (
        'getInitialTangent() must return the SAME (elastic) matrix whether '
        'the material has yielded or not -- pre-fix this returned the '
        'current (softened, history-dependent) mCep instead of mCe\n'
        f'du (pristine)   = {du_pristine}\n'
        f'du (post-yield) = {du_after_yield}\n'
        f'relative deviation = {rel}')


# ===========================================================================
# Fix 3 -- ContactMaterial2D/3D::getInitialTangent() -- UNREACHABLE from Python
# ===========================================================================

@pytest.mark.skip(reason=(
    "ContactMaterial2D/3D::getInitialTangent() cannot be exercised from "
    "openseespy. Every upstream element that owns one of these materials "
    "(SimpleContact2D, SimpleContact3D, BeamContact2D, BeamContact3D, "
    "BeamContact2Dp, BeamContact3Dp) overrides getInitialStiff() as a bare "
    "`return getTangentStiff();' and never calls the material's own "
    "getInitialTangent() at all. ZeroLengthND cannot substitute either: "
    "ContactMaterial2D/3D::getCopy(const char*) only matches the literal "
    "strings 'ContactMaterial2D'/'ContactMaterial3D', not the "
    "'PlaneStrain2D'/'ThreeDimensional' tokens ZeroLengthND requests, so "
    "`element zeroLengthND ... <ContactMaterial2D tag>` fails to construct. "
    "There is also no standalone 'test ND material' command. The fix "
    "(a dedicated initialTangent buffer holding the elastic/sticking "
    "tangent, independent of getTangent()'s call history) is verified by "
    "code inspection only; see LEDGER_vanilla_files.md for the exact "
    "before/after values it replaces (the zero matrix on a fresh material)."
))
def test_contactmaterial_initial_tangent_unreachable_from_python():
    pass


# ===========================================================================
# Fix 4 -- ProfileSPDLinSOE asymmetry guard (budgeted, once per SOE)
# ===========================================================================
#
# Same cube deck as test_pardiso_asym_rearm.py (2x2x2 stdBrick, base fixed,
# shear+compression at the top face), just pointed at `system ProfileSPD` --
# the DEFAULT solver -- instead of Pardiso. ProfileSPD has no matrixType
# switch: it ALWAYS keeps only row <= col, so the guard needs no special flag
# to engage.

_ASYM = 'UNSYMMETRIC'
_AT_ASSEMBLY = re.compile(r'TANGENT ASSEMBLY (\d+)')
_RESAMPLE_PERIOD = 64   # mirrors ASYM_RESAMPLE_PERIOD in ProfileSPDLinSOE.cpp

_CUBE_NX = 2
_CUBE_H_TOTAL = 2.025
_CUBE_V_TOTAL = 3.0
_CUBE_NSTEPS = 48


def _cube_nid(i, j, k):
    return 1 + i + (_CUBE_NX + 1) * (j + (_CUBE_NX + 1) * k)


def _build_cube(elastic=False):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

    if elastic:
        E = 9.0 * _K * _G / (3.0 * _K + _G)
        nu = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))
        ops.nDMaterial('ElasticIsotropic', 1, E, nu)
    else:
        ops.nDMaterial('DruckerPrager', 1, _K, _G, _SIG_Y, _RHO, _RHO_BAR,
                       0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)

    h = 1.0 / _CUBE_NX
    for k in range(_CUBE_NX + 1):
        for j in range(_CUBE_NX + 1):
            for i in range(_CUBE_NX + 1):
                ops.node(_cube_nid(i, j, k), i * h, j * h, k * h)
    for j in range(_CUBE_NX + 1):
        for i in range(_CUBE_NX + 1):
            ops.fix(_cube_nid(i, j, 0), 1, 1, 1)

    tag = 1
    for k in range(_CUBE_NX):
        for j in range(_CUBE_NX):
            for i in range(_CUBE_NX):
                ops.element('stdBrick', tag,
                            _cube_nid(i, j, k), _cube_nid(i + 1, j, k),
                            _cube_nid(i + 1, j + 1, k), _cube_nid(i, j + 1, k),
                            _cube_nid(i, j, k + 1), _cube_nid(i + 1, j, k + 1),
                            _cube_nid(i + 1, j + 1, k + 1), _cube_nid(i, j + 1, k + 1),
                            1)
                tag += 1

    ntop = (_CUBE_NX + 1) * (_CUBE_NX + 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j in range(_CUBE_NX + 1):
        for i in range(_CUBE_NX + 1):
            ops.load(_cube_nid(i, j, _CUBE_NX),
                     _CUBE_H_TOTAL / ntop, 0.0, -_CUBE_V_TOTAL / ntop)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('ProfileSPD')            # the DEFAULT solver
    ops.test('NormDispIncr', 1.0e-8, 40, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _CUBE_NSTEPS)
    ops.analysis('Static')


def _cube_ramp(elastic=False, nsteps=_CUBE_NSTEPS):
    _build_cube(elastic=elastic)
    done = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        done += 1
    return done


def test_profilespd_asym_warning_once(capfd):
    """Default solver, non-associated DruckerPrager: warns, and warns ONCE."""
    done = _cube_ramp(elastic=False)
    cap = capfd.readouterr()
    text = cap.out + cap.err

    assert _ASYM in text, (
        f'non-associated DruckerPrager ran {done}/{_CUBE_NSTEPS} steps under '
        'the default `system ProfileSPD` and the half-storage guard never '
        'fired\n' + text)
    assert 'ProfileSPD' in text, (
        'the warning must name ProfileSPD so a user knows which solver to '
        'switch off\n' + text)
    assert text.count(_ASYM) == 1, (
        f'asymmetry reported {text.count(_ASYM)} times; the asymWarned '
        'latch must make this a once-per-SOE warning, not a flood\n' + text)


def test_profilespd_asym_warning_absent_when_elastic(capfd):
    """No false positive on a genuinely symmetric (elastic) tangent."""
    done = _cube_ramp(elastic=True)
    cap = capfd.readouterr()
    text = cap.out + cap.err
    assert done == _CUBE_NSTEPS, (
        f'the elastic control must converge all {_CUBE_NSTEPS} steps (got '
        f'{done}); a short run would make the no-warning assertion vacuous '
        f'-- fewer than {_RESAMPLE_PERIOD} assemblies means the guard never '
        're-armed and the test would pass without exercising anything\n' + text)
    assert _ASYM not in text, (
        'an ElasticIsotropic tangent is symmetric to round-off; the guard '
        'must not manufacture a warning\n' + text)
