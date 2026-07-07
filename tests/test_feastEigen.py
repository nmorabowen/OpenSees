"""FeastEigenSOE/FeastEigenSolver (ADR 43 P1, tags 33022/33023) — G-B gate.

Serial MKL-FEAST band-targeted eigensolver: `eigen -feast fmin fmax` answers
"ALL modes in [fmin, fmax] Hz" via contour integration (MKL Extended
Eigensolver dfeast_scsrgv on the fork-assembled CSR K, M).

Gate G-B (umbrella ADR 45 §6): same eigenpairs as ARPACK on a medium model
(lambda <=1e-6 rel, MAC >=0.999 per mode) AND band-targeting returns exactly
the modes whose f lies in [f1, f2] against a fullGenLapack full-spectrum
oracle — interior, edge-straddling, and empty bands; m0 saturation is
flagged and auto-enlarged, never silently truncated.

  test_feast_matches_arpack     G-B parity: lambda + MAC vs ARPACK
  test_feast_band_interior      exactly modes 3..6 of the full spectrum
  test_feast_band_empty         no modes in band -> empty result, no crash
  test_feast_saturation_grows   -m0 too small -> warns + still complete
  test_feast_bad_input          malformed band rejected
"""
import math
import sys

import numpy as np
import pytest

from _testbed import ops

# FEAST lives in the MKL Extended Eigensolver: linked on the Windows/oneAPI
# build only (Zone-A Ubuntu CI uses reference LAPACK — the classes compile
# there but solve() refuses). The G-B gate is enforced on the Windows box.
pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="FEAST requires MKL (Windows/oneAPI build)"),
]

TWO_PI = 2.0 * math.pi


def _frame(nx=3, ny=10):
    """2D elastic frame grid (~(nx+1)*(ny+1)*3 DOF) with lumped masses."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    E, A, Iz, m = 30.0e9, 0.09, 6.75e-4, 5.0e3
    tag = 0
    for j in range(ny + 1):
        for i in range(nx + 1):
            tag += 1
            ops.node(tag, i * 5.0, j * 3.0)
            if j == 0:
                ops.fix(tag, 1, 1, 1)
            else:
                ops.mass(tag, m, m, 0.0)
    ops.geomTransf("Linear", 1)
    ele = 0
    def nid(i, j):
        return j * (nx + 1) + i + 1
    for j in range(1, ny + 1):        # columns
        for i in range(nx + 1):
            ele += 1
            ops.element("elasticBeamColumn", ele, nid(i, j - 1), nid(i, j),
                        A, E, Iz, 1)
    for j in range(1, ny + 1):        # beams
        for i in range(nx):
            ele += 1
            ops.element("elasticBeamColumn", ele, nid(i, j), nid(i + 1, j),
                        A, E, Iz, 1)


def _free_nodes(nx=3, ny=10):
    return [j * (nx + 1) + i + 1 for j in range(1, ny + 1)
            for i in range(nx + 1)]


def _mode_matrix(nodes, mode):
    """Stacked (ux, uy) eigenvector entries for MAC comparisons."""
    out = []
    for n in nodes:
        v = ops.nodeEigenvector(n, mode)
        out.extend(v[:2])
    return np.asarray(out)


def _mac(a, b):
    return float((a @ b) ** 2 / ((a @ a) * (b @ b)))


def _hz(lam):
    return math.sqrt(lam) / TWO_PI


def test_feast_matches_arpack():
    """G-B parity: FEAST over a band enclosing exactly the first 8 modes
    reproduces ARPACK's eigenpairs (lambda <=1e-6 rel, MAC >= 0.999)."""
    _frame()
    lam_ref = np.asarray(ops.eigen(10))
    nodes = _free_nodes()
    phi_ref = [_mode_matrix(nodes, k + 1) for k in range(8)]
    f_hi = 0.5 * (_hz(lam_ref[7]) + _hz(lam_ref[8]))  # between f8 and f9

    lam_feast = np.asarray(ops.eigen("-feast", 0.0, f_hi))
    assert lam_feast.shape == (8,)
    np.testing.assert_allclose(lam_feast, lam_ref[:8], rtol=1e-6)
    for k in range(8):
        phi_f = _mode_matrix(nodes, k + 1)
        assert _mac(phi_f, phi_ref[k]) >= 0.999, f"MAC mode {k + 1}"


def test_feast_band_interior():
    """Interior and edge-straddling bands against the fullGenLapack full
    spectrum: exactly the enclosed modes, no more, no fewer."""
    _frame(nx=1, ny=3)   # small: full spectrum is cheap
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 12))
    f = np.array([_hz(l) for l in lam_all])
    # interior band [f3, f6] with margins to the neighbours
    lo = 0.5 * (f[1] + f[2])
    hi = 0.5 * (f[5] + f[6])
    lam_band = np.asarray(ops.eigen("-feast", lo, hi))
    np.testing.assert_allclose(lam_band, lam_all[2:6], rtol=1e-6)
    # edge-straddling: cut through mode 4's frequency -> 3..4 only
    hi_cut = 0.5 * (f[3] + f[4])
    lam_cut = np.asarray(ops.eigen("-feast", lo, hi_cut))
    np.testing.assert_allclose(lam_cut, lam_all[2:4], rtol=1e-6)


def test_feast_band_empty():
    """A band between two modes holds nothing: graceful empty result."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 6))
    f = np.array([_hz(l) for l in lam_all])
    lo = f[2] * 1.02
    hi = f[3] * 0.98
    out = ops.eigen("-feast", lo, hi)
    empty = (out is None or out == 0 or
             (hasattr(out, "__len__") and len(out) == 0))
    assert empty, f"expected no modes in ({lo:.3f},{hi:.3f}) Hz, got {out}"


def test_feast_saturation_grows(capfd):
    """-m0 deliberately undersized: the solver must WARN and auto-enlarge
    (never silently truncate the band)."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 8))
    f = np.array([_hz(l) for l in lam_all])
    hi = 0.5 * (f[6] + f[7])          # band holds 7 modes
    lam = np.asarray(ops.eigen("-feast", 0.0, hi, "-m0", 3))
    err = capfd.readouterr().err.lower()
    assert lam.shape == (7,)
    np.testing.assert_allclose(lam, lam_all[:7], rtol=1e-6)
    assert "m0" in err or "subspace" in err or "enlarg" in err


def test_feast_mass_orthonormal():
    """ADR 43 §8.1: Phi^T M Phi == I. The masses are purely nodal here, so
    the reduced Gram assembles exactly from nodeEigenvector + the known m —
    this pins mass-normalization AND is far more sensitive to a CSR scatter
    error than eigenvalue parity on a symmetric grid (adversarial-gate
    proposal)."""
    _frame()
    m = 5.0e3
    lam_ref = np.asarray(ops.eigen(8))
    f_hi = _hz(lam_ref[6]) * 0.5 + _hz(lam_ref[7]) * 0.5
    lam = np.asarray(ops.eigen("-feast", 0.0, f_hi))
    p = lam.shape[0]
    assert p == 7
    nodes = _free_nodes()
    phi = np.array([[c for n in nodes
                     for c in ops.nodeEigenvector(n, k + 1)[:2]]
                    for k in range(p)])          # p x (2*nfree)
    gram = m * (phi @ phi.T)
    np.testing.assert_allclose(gram, np.eye(p), atol=1e-8)


def test_feast_bad_input():
    _frame(nx=1, ny=3)
    with pytest.raises(Exception):
        ops.eigen("-feast", 5.0)              # missing fmax
    with pytest.raises(Exception):
        ops.eigen("-feast", 10.0, 5.0)        # inverted band
