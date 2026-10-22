"""WP-203: eigen -feast + modalProperties, and FEAST's stop criterion.

(1) The FEAST reconcile trimmed the domain eigenvalue list to the found count
    but left every node's eigenvector matrix at the analysis-loop cap (>= 128
    columns), so modalProperties aborted with "inconsistent eigenvector matrix
    and ID" under every constraint handler. Oracle: the same model through
    `eigen -fullGenLapack` + modalProperties (periods, participation masses).
(2) FEAST now stops on MKL's relative-residual test (fpm(6)=1); with the trace
    default every solve ran all refinement loops on a building model. The
    eigenpairs must not depend on the criterion (-traceStop restores it).
"""
import math
import sys

import pytest

from _testbed import ops

# ci-coverage: local-only -- needs FEAST (MKL, Windows/oneAPI build) (WP-143)
pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="FEAST requires MKL (Windows/oneAPI build)"),
]


def _building(handler):
    """3D 2x2-bay, 4-storey frame; rigid diaphragms (MP constraints) per floor."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    E, G = 25.0e6, 10.0e6
    tag = lambda k, i, j: 1 + i + 3 * j + 9 * k
    for k in range(5):
        for j in range(3):
            for i in range(3):
                ops.node(tag(k, i, j), 6.0 * i, 5.0 * j, 3.0 * k)
                if k == 0:
                    ops.fix(tag(k, i, j), 1, 1, 1, 1, 1, 1)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    ops.geomTransf('Linear', 2, 0.0, 0.0, 1.0)
    e = 0
    for k in range(1, 5):
        for j in range(3):
            for i in range(3):
                e += 1
                ops.element('elasticBeamColumn', e, tag(k - 1, i, j), tag(k, i, j),
                            0.25, E, G, 0.009, 0.0052, 0.0052, 1)
        for j in range(3):
            for i in range(2):
                e += 1
                ops.element('elasticBeamColumn', e, tag(k, i, j), tag(k, i + 1, j),
                            0.12, E, G, 0.002, 0.0016, 0.0036, 2)
        for i in range(3):
            for j in range(2):
                e += 1
                ops.element('elasticBeamColumn', e, tag(k, i, j), tag(k, i, j + 1),
                            0.12, E, G, 0.002, 0.0016, 0.0036, 2)
        master = 1000 + k
        ops.node(master, 6.0, 5.0, 3.0 * k)
        ops.fix(master, 0, 0, 1, 1, 1, 0)
        ops.mass(master, 50.0, 50.0, 0.0, 0.0, 0.0, 50.0 * (12.0 ** 2 + 10.0 ** 2) / 12.0)
        ops.rigidDiaphragm(3, master, *[tag(k, i, j) for j in range(3) for i in range(3)])
    if handler == 'Penalty':
        ops.constraints('Penalty', 1.0e14, 1.0e14)
    else:
        ops.constraints(handler)
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')


def _modal(eigen_args):
    lam = ops.eigen(*eigen_args)
    mp = ops.modalProperties('-return')
    return list(lam), mp


@pytest.mark.parametrize("handler", ['Transformation', 'Penalty'])
def test_feast_modal_properties_match_lapack(handler):
    # reference: the exact (Transformation) model through dense LAPACK. A
    # 1e14 penalty is ill-conditioned for the DENSE solver (fullGenLapack
    # drifts 7e-5 on mode 1 while ARPACK and FEAST agree with Transformation
    # to 1e-5), so the Penalty run is held to penalty accuracy.
    _building('Transformation')
    lam_ref, mp_ref = _modal(['-fullGenLapack', 12])
    f_cut = 0.5 * (math.sqrt(lam_ref[5]) + math.sqrt(lam_ref[6])) / (2.0 * math.pi)
    _building(handler)
    lam_f, mp_f = _modal(['-feast', 0.0, f_cut])          # used to abort here
    rel, absr = (1e-8, 1e-6) if handler == 'Transformation' else (2e-5, 1e-4)
    assert len(lam_f) == 6
    for a, b in zip(lam_f, lam_ref[:6]):
        assert a == pytest.approx(b, rel=rel)
    for key in ('partiMassRatiosMX', 'partiMassRatiosMY', 'partiMassRatiosRMZ'):
        for a, b in zip(mp_f[key], mp_ref[key][:6]):
            assert a == pytest.approx(b, abs=absr), key
    for n in ops.getNodeTags():
        # node eigenvector storage trimmed to the found modes
        with pytest.raises(Exception):
            ops.nodeEigenvector(n, 7)


def test_rsa_after_feast():
    _building('Transformation')
    lam = ops.eigen('-feast', 0.0, 5.0)
    ops.modalProperties()
    for k in range(1, len(lam) + 1):
        ops.responseSpectrumAnalysis(1, '-Tn', 0.1, 1.0, 4.0, '-Sa', 2.0, 2.0, 0.5, '-mode', k)
    assert any(abs(ops.nodeDisp(1004, 1)) > 0.0 for _ in [0])


def test_stop_criterion_does_not_change_eigenpairs():
    _building('Transformation')
    lam_r = ops.eigen('-feast', 0.0, 6.0)
    _building('Transformation')
    lam_t = ops.eigen('-feast', 0.0, 6.0, '-traceStop')
    assert len(lam_r) == len(lam_t) > 0
    for a, b in zip(lam_r, lam_t):
        assert a == pytest.approx(b, rel=1e-8)
