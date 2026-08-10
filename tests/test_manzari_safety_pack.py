"""Regression tests for the ManzariDafalias/tet safety pack (TIMs report).

Four silent traps were identified in the upstream UW ``ManzariDafalias``
material and its ``TenNodeTetrahedron`` sibling, all fixed with ``// Ladruno``
markers in ``SRC/material/nD/UWmaterials/ManzariDafalias.cpp`` and
``SRC/element/tetrahedron/TenNodeTetrahedron.cpp``:

  1. (item 7) ``mElastFlag`` -- a STATIC, class-wide flag -- defaulted to
     ELASTOPLASTIC (1) in the two full constructors used by the Tcl/Python
     ``nDMaterial ManzariDafalias`` command, even though the community
     staged-gravity idiom (elastic gravity, then
     ``updateMaterialStage -material N -stage 1``) assumes the material
     starts elastic. Without an explicit stage-0 call, "elastic" gravity ran
     the full plastic model from zero stress (a singular regime) and stalled
     Newton with no diagnostic. Fixed by setting ``mElastFlag = 0`` in both
     full constructors, matching the parallel/classTag constructor that
     already did this.
  2. (item 8) The "Outside Bounding" diagnostic in ``Elastic2Plastic()`` was
     commented out upstream, so a stage-switch stress ratio above the
     bounding surface silently bumps ``M_c`` by 10% with no trace. Restored
     as an ``opserr`` warning (numerics unchanged).
  3. (item 3/4) Integration schemes 3 (RungeKutta4, whose adaptive branch is
     dead code) and 5 (ForwardEuler) take one uncorrected substep with no
     error estimate and no yield-drift correction. Added a once-per-process
     warning in both full constructors.
  4. (item 8) ``TenNodeTetrahedron::update()`` discarded the per-Gauss-point
     ``setTrialStrain`` return code and always returned 0, so the element
     could never report a constitutive failure (unlike its fork sibling
     BezierTet10, which accumulates the sum). Fixed to accumulate and return.

Test 1 is the load-bearing fingerprint for item 7: pre-fix, this exact deck
(no ``updateMaterialStage`` call) runs the plastic model from zero stress and
is at best highly nonlinear, at worst a Newton stall. Test 2 is a smoke test
that the staged-gravity idiom still works after the fix. Test 3 is the
fingerprint for item 3/4's warning text.
"""
import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

# ---- material: canonical Toyoura-sand ManzariDafalias parameter set ------
# (OpenSees docs example; units are self-consistent, P_atm in the same
# stress units as the applied loads below).

_MD_PARAMS = dict(
    G0=125, nu=0.05, e_init=0.8, Mc=1.25, c=0.712, lambda_c=0.019, e0=0.934,
    ksi=0.7, P_atm=100, m=0.01, h0=7.05, ch=0.968, nb=1.1, A0=0.704, nd=3.5,
    z_max=4, cz=600, Rho=1.42,
)
_MD_ORDER = ["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi",
             "P_atm", "m", "h0", "ch", "nb", "A0", "nd", "z_max", "cz", "Rho"]


def _md_args():
    return [_MD_PARAMS[k] for k in _MD_ORDER]


# ---- Test 3 (kept first): the scheme 3/5 "no error control" warning ------
#
# The warning latch (a C++ static bool) fires only ONCE per process. This
# must be the first test in the file/session to construct a scheme-3 or -5
# ManzariDafalias so it is the one that observes the warning text. No other
# test in this file (or, as of writing, anywhere else in tests/) constructs
# IntScheme 3 or 5.

def test_scheme3_warns_no_error_control(capfd):
    """IntScheme=3 (RungeKutta4) must emit the "no error control" warning."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ManzariDafalias', 1, *_md_args(), 3, 2, 1, 1e-7, 1e-7)
    err = capfd.readouterr().err
    assert 'IntScheme 3' in err and 'no error estimate' in err, (
        'expected the scheme-3/5 no-error-control warning on stderr', err)
    assert 'no yield-drift correction' in err
    assert 'Use 1 (ModifiedEuler' in err


# ---- Test 1: elastic default + linearity fingerprint (TIMs item 7) -------

_NSTEPS = 4


def _build_force_cube(mat_tag, fz_total):
    """Unit stdBrick cube: base fixed, top rollered in x/y, free/loaded in z."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

    xy = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
    for k in range(2):                       # k=0 bottom (z=0), k=1 top (z=1)
        for j, (x, y) in enumerate(xy):
            ops.node(4 * k + j + 1, x, y, float(k))

    ops.nDMaterial('ManzariDafalias', mat_tag, *_md_args())
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)

    for j in range(4):                       # bottom face: fully fixed
        ops.fix(j + 1, 1, 1, 1)
    for j in range(4):                       # top face: roller in x/y
        ops.fix(4 + j + 1, 1, 1, 0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j in range(4):
        ops.load(4 + j + 1, 0.0, 0.0, fz_total / 4.0)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-10, 15, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _NSTEPS)
    ops.analysis('Static')


def _run_force_cube(mat_tag, fz_total):
    _build_force_cube(mat_tag, fz_total)
    for step in range(_NSTEPS):
        rc = ops.analyze(1)
        assert rc == 0, (
            f'step {step + 1}/{_NSTEPS} failed to converge (fz_total={fz_total}) -- '
            'pre-fix this deck runs the plastic model from zero stress and stalls'
        )
        iters = ops.testIter()
        assert iters <= 5, (
            f'step {step + 1}/{_NSTEPS} took {iters} iterations -- an elastic-stage '
            'material should converge in a handful (the tangent equals the true '
            'response, so Newton is exact)'
        )
    return ops.nodeDisp(5, 3)


def test_elastic_default_converges_and_is_linear():
    """No updateMaterialStage call: the material must start elastic (mElastFlag=0).

    Convergence in a handful of iterations AND an exactly linear force-vs-
    displacement response are the fingerprint that the material never
    touched the plastic (explicit/implicit substepping) code path at all.
    """
    fz1 = -20.0
    ux1 = _run_force_cube(1, fz1)

    fz2 = 2.0 * fz1
    ux2 = _run_force_cube(2, fz2)

    assert ux1 != 0.0
    rel = abs(ux2 - 2.0 * ux1) / abs(ux2)
    assert rel < 1.0e-6, (
        'doubling the load did not double the displacement -- the material '
        'is not behaving as linear-elastic', ux1, ux2, rel)


# ---- Test 2: updateMaterialStage still reaches the plastic path -----------
#
# Confined + moderate deviatoric strain ratio, same idiom (and same 2-element
# stack topology, for the same reason) as
# tests/test_pdmy01_plastic_strain.py: rollers on the three negative faces,
# prescribed displacement on the three positive faces, applied on BOTH
# stdBrick elements of a 1x1x2 stack so the shared mid-face node z-DOFs stay
# free -- a single cube with all three positive/negative face pairs
# prescribed has zero free DOFs and never actually exercises Newton. DELTA is
# chosen so the elastic stress ratio stays comfortably below Mc (=1.25) even
# after the stage flip -- this is a smoke test that the plastic integrator
# code path (explicit_integrator / BackwardEuler_CPPM, entered via
# Elastic2Plastic() at the flip) runs without crashing, not a check of
# plastic-regime accuracy.

_E0 = 2.0e-4
_DELTA = 1.5
_NSTEP_STAGE = 8
_XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]


def _build_sp_stack(mat_tag):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, 0.5 * k)

    ops.nDMaterial('ManzariDafalias', mat_tag, *_md_args())
    for e in (1, 2):
        b = 4 * (e - 1)
        ops.element('stdBrick', e, b + 1, b + 2, b + 3, b + 4, b + 5, b + 6, b + 7, b + 8, mat_tag)

    # rollers on the three negative faces
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)

    # prescribed uniform strain on the three positive faces
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -_E0)
            if y == 1.:
                ops.sp(n, 2, -_E0)
            if k == 2:
                ops.sp(n, 3, -_E0 * (1.0 + _DELTA))

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-10, 60, 0)
    ops.algorithm('KrylovNewton')
    ops.integrator('LoadControl', 1.0 / _NSTEP_STAGE)
    ops.analysis('Static')


def test_update_material_stage_still_reaches_plastic_path():
    """Elastic gravity, flip stage, keep loading -- the staged idiom must work."""
    _build_sp_stack(3)

    for step in range(_NSTEP_STAGE):
        assert ops.analyze(1) == 0, f'elastic-stage step {step + 1} failed'

    ops.updateMaterialStage('-material', 3, '-stage', 1)

    for step in range(_NSTEP_STAGE):
        assert ops.analyze(1) == 0, (
            f'plastic-stage step {step + 1} failed after updateMaterialStage -- '
            'the staged-gravity idiom must remain reachable'
        )
