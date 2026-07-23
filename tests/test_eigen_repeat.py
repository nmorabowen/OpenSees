"""Regression: repeated `eigen` on the interpreter (openseespy) path.

Upstream bug (present in OpenSees/OpenSees master, fixed in this fork):
`OpenSeesCommands::eigen` builds an EPHEMERAL DirectIntegrationAnalysis when
the script has not defined an analysis, deletes it after the solve (no-op
destructor, so the cached theEigenSOE survives), and on the NEXT call builds
a fresh analysis whose internal EigenSOE is null.  With a same-classTag
cached SOE `eigenSOEUpdated` stays false and setEigenSOE is never re-invoked,
so every second consecutive `eigen` failed with
"DirectIntegrationAnalysis::eigen() - no EigenSOE has been set".

Classic Tcl re-attaches unconditionally and never had the bug.
"""

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

_REL_TOL = 1.0e-8


def _build_model():
    """Truss chain with lateral springs: 20 free DOFs, well-separated modes."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.uniaxialMaterial('Elastic', 1, 3000.0)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1)
    for i in range(1, 11):
        ops.node(i, 10.0 * i, 0.0)
        ops.mass(i, 5.0, 5.0)
        ops.element('Truss', i, i - 1, i, 10.0, 1)
        ops.node(100 + i, 10.0 * i, 10.0)
        ops.fix(100 + i, 1, 1)
        ops.element('Truss', 100 + i, 100 + i, i, 4.0, 1)


def _assert_close(a, b, tol=_REL_TOL):
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x > 0.0
        assert abs(x - y) <= tol * abs(x), (x, y)


def test_repeated_eigen_no_analysis_same_type():
    """The original failure: second consecutive default eigen call raised."""
    _build_model()
    first = ops.eigen(3)
    second = ops.eigen(3)
    third = ops.eigen(3)
    _assert_close(first, second)
    _assert_close(first, third)


def test_repeated_eigen_no_analysis_explicit_type():
    _build_model()
    first = ops.eigen('-genBandArpack', 3)
    second = ops.eigen('-genBandArpack', 3)
    _assert_close(first, second)


def test_repeated_eigen_no_analysis_switched_type():
    """Type switch exercises the SOE replacement path, then same-type reuse."""
    _build_model()
    arpack = ops.eigen('-genBandArpack', 3)
    lapack = ops.eigen('-fullGenLapack', 3)
    back = ops.eigen('-genBandArpack', 3)
    _assert_close(arpack, lapack, tol=1.0e-6)
    _assert_close(arpack, back)


def test_repeated_eigen_with_defined_transient_analysis():
    """Persistent-analysis flow must keep working (never had the bug)."""
    _build_model()
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-8, 10)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    first = ops.eigen(3)
    second = ops.eigen(3)
    _assert_close(first, second)


def test_eigen_then_define_analysis_then_eigen():
    """No-analysis eigen followed by a real analysis, then eigen again."""
    _build_model()
    before = ops.eigen(3)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-8, 10)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    after = ops.eigen(3)
    again = ops.eigen(3)
    _assert_close(before, after)
    _assert_close(before, again)
