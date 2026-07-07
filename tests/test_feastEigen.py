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

P2 (-certify, the Sturm/inertia completeness certificate — ADR 43 §8.4):
count the band content independently via the negative-pivot inertia of the
LDL^T factorization of (K - sigma*M) at the two band edges (MKL PARDISO on
the SOE's own CSR); refuse on mismatch with FEAST's count.

  test_certify_matches          certified == FEAST count on the frame band
  test_certify_closely_spaced   near-degenerate pair straddling a band edge:
                                exactly one counted, one excluded (§8.2)
  test_certify_empty_band       0 == 0 certifies
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


# ---------------------------------------------------------------- P2 certify
def test_certify_matches(capfd):
    """-certify: the inertia count must equal FEAST's m on a normal band,
    and say so."""
    _frame()
    lam_ref = np.asarray(ops.eigen(8))
    f_hi = 0.5 * (_hz(lam_ref[6]) + _hz(lam_ref[7]))
    lam = np.asarray(ops.eigen("-feast", 0.0, f_hi, "-certify"))
    assert lam.shape == (7,)
    err = capfd.readouterr().err
    assert "certified complete" in err


def test_certify_closely_spaced(capfd):
    """ADR 43 §8.2/§8.4: two nearly-identical uncoupled oscillators; a band
    edge cutting BETWEEN the close pair must certify exactly one of them —
    the inertia count is exact where a loose contour could smear."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    k, m, split = 100.0, 1.0, 1.0e-3
    # two independent SDOFs: k and k*(1+split)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 1.0); ops.mass(2, m)
    ops.node(3, 2.0); ops.fix(3, 1)
    ops.node(4, 3.0); ops.mass(4, m)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.uniaxialMaterial("Elastic", 2, k * (1.0 + split))
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 3, 4, 1.0, 2)
    f1 = _hz(k / m)
    f2 = _hz(k * (1.0 + split) / m)
    f_cut = 0.5 * (f1 + f2)          # between the close pair
    lam = np.asarray(ops.eigen("-feast", 0.0, f_cut, "-certify"))
    assert lam.shape == (1,)
    assert lam[0] == pytest.approx(k / m, rel=1e-10)
    assert "certified complete" in capfd.readouterr().err
    # and the full band certifies both
    lam2 = np.asarray(ops.eigen("-feast", 0.0, f2 * 1.1, "-certify"))
    assert lam2.shape == (2,)


def test_certify_edge_on_eigenvalue(capfd):
    """Adversarial-gate fix: a band edge sitting numerically ON an
    eigenvalue makes (K - sigma*M) singular at the shift — the crossing
    pivot emerges near-zero in the ELIMINATION (a coupled model: an
    isolated diagonal zero would be normalized away by PARDISO's scaling),
    gets perturbed, and the inertia sign is no longer guaranteed. The
    certificate must detect the perturbed pivots and REFUSE rather than
    certify a possibly off-by-one count."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    # coupled 2-DOF chain (NOT uncoupled oscillators — see docstring)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 1.0); ops.mass(2, 2.0)
    ops.node(3, 2.0); ops.mass(3, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 400.0)
    ops.uniaxialMaterial("Elastic", 2, 150.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 2, 3, 1.0, 2)
    lam = np.asarray(ops.eigen("-fullGenLapack", 2))
    f1 = _hz(lam[0])                # edge EXACTLY on the first eigenvalue
    with pytest.raises(Exception):
        ops.eigen("-feast", 0.0, f1, "-certify")
    assert "perturbed pivot" in capfd.readouterr().err


def test_certify_empty_band(capfd):
    """0 == 0 certifies (no modes in the band, inertia agrees)."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 6))
    f = np.array([_hz(l) for l in lam_all])
    out = ops.eigen("-feast", f[2] * 1.02, f[3] * 0.98, "-certify")
    empty = (out is None or out == 0 or
             (hasattr(out, "__len__") and len(out) == 0))
    assert empty
    assert "certified complete" in capfd.readouterr().err


def test_feast_bad_input():
    _frame(nx=1, ny=3)
    with pytest.raises(Exception):
        ops.eigen("-feast", 5.0)              # missing fmax
    with pytest.raises(Exception):
        ops.eigen("-feast", 10.0, 5.0)        # inverted band


# ------------------------------------------------------------- P3b blockZGate
def test_blockzgate_frame(capfd):
    """P3b: the block-real complex-shift solve kernel (the FEAST-RCI inner
    solve of P3c) self-tests on the frame's own CSR — for a shift at the
    widest spectral gap it sweeps the contour flattening b/r down to 1e-6,
    recovering a KNOWN complex solution through the symmetric 2n block-real
    LDL^T path, and REFUSES if the error exceeds 1e-8. The native
    complex-symmetric PARDISO baseline runs beside it (the D2-review cost/
    accuracy record), and the a-on-eigenvalue degradation regime is reported.
    The eigen call succeeding IS the gate."""
    _frame()
    lam = np.asarray(ops.eigen("-feast", 0.0, 30.0, "-blockZGate"))
    assert lam.size >= 1
    err = capfd.readouterr().err
    assert "-blockZGate PASS" in err
    assert "degradation regime" in err


def test_blockzgate_small_model(capfd):
    """P3b on the small frame with a band that holds a handful of modes —
    exercises the gap-placement logic with interior eigenvalues present."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 8))
    f = np.array([_hz(l) for l in lam_all])
    hi = 0.5 * (f[5] + f[6])
    lam = np.asarray(ops.eigen("-feast", 0.0, hi, "-blockZGate"))
    assert lam.size == 6
    assert "-blockZGate PASS" in capfd.readouterr().err


def test_blockzgate_empty_band(capfd):
    """P3b on an empty band (info=1, m=0): the gate must still run — the
    empty-band path falls THROUGH to the hook rather than returning early —
    self-testing a band-midpoint shift with no interior eigenvalues."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 6))
    f = np.array([_hz(l) for l in lam_all])
    out = ops.eigen("-feast", f[2] * 1.02, f[3] * 0.98, "-blockZGate")
    empty = (out is None or out == 0 or
             (hasattr(out, "__len__") and len(out) == 0))
    assert empty
    assert "-blockZGate PASS" in capfd.readouterr().err


def test_blockzgate_composes_with_certify(capfd):
    """-certify and -blockZGate on the same run: both diagnostics fire and
    the eigenpairs are unchanged vs the plain band solve."""
    _frame(nx=1, ny=3)
    lam_plain = np.asarray(ops.eigen("-feast", 0.0, 50.0))
    lam_gated = np.asarray(ops.eigen("-feast", 0.0, 50.0, "-certify",
                                     "-blockZGate"))
    np.testing.assert_allclose(lam_gated, lam_plain, rtol=1e-12)
    err = capfd.readouterr().err
    assert "certified complete" in err
    assert "-blockZGate PASS" in err


# ------------------------------------------------------------------ P3c -rci
# dfeast_srci RCI orchestration with LadrunoBlockZKernel as the inner
# complex-shift solve (the seam the MPI rung distributes). Gate: the RCI
# path reproduces the dfeast_scsrgv driver-path eigenpairs (lambda <=1e-8
# rel, MAC >=0.999) on the frame battery, plus -certify composition.

def test_rci_matches_driver_frame():
    """P3c gate: -rci eigenpairs == driver-path eigenpairs on the frame
    (lambda <=1e-8 rel, MAC >=0.999 per mode)."""
    _frame()
    lam_ref = np.asarray(ops.eigen(10))
    f_hi = 0.5 * (_hz(lam_ref[7]) + _hz(lam_ref[8]))

    lam_drv = np.asarray(ops.eigen("-feast", 0.0, f_hi))
    nodes = _free_nodes()
    phi_drv = [_mode_matrix(nodes, k + 1) for k in range(8)]

    lam_rci = np.asarray(ops.eigen("-feast", 0.0, f_hi, "-rci"))
    assert lam_rci.shape == lam_drv.shape == (8,)
    np.testing.assert_allclose(lam_rci, lam_drv, rtol=1e-8)
    for k in range(8):
        phi_r = _mode_matrix(nodes, k + 1)
        assert _mac(phi_r, phi_drv[k]) >= 0.999, f"MAC mode {k + 1}"


def test_rci_band_interior():
    """-rci band-targeting: exactly the enclosed modes of the full spectrum
    (interior and edge-straddling cuts), like the driver path."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 12))
    f = np.array([_hz(l) for l in lam_all])
    lo = 0.5 * (f[1] + f[2])
    hi = 0.5 * (f[5] + f[6])
    lam_band = np.asarray(ops.eigen("-feast", lo, hi, "-rci"))
    np.testing.assert_allclose(lam_band, lam_all[2:6], rtol=1e-6)
    hi_cut = 0.5 * (f[3] + f[4])
    lam_cut = np.asarray(ops.eigen("-feast", lo, hi_cut, "-rci"))
    np.testing.assert_allclose(lam_cut, lam_all[2:4], rtol=1e-6)


def test_rci_empty_band():
    """-rci on an empty band: graceful empty result (info=1 path through
    the RCI termination)."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 6))
    f = np.array([_hz(l) for l in lam_all])
    out = ops.eigen("-feast", f[2] * 1.02, f[3] * 0.98, "-rci")
    empty = (out is None or out == 0 or
             (hasattr(out, "__len__") and len(out) == 0))
    assert empty, f"expected no modes, got {out}"


def test_rci_saturation_grows(capfd):
    """-rci with a deliberately undersized -m0: the saturation loop must
    enlarge (rebuilding the RCI workspace at the new m0) and return the
    complete band, exactly like the driver path."""
    _frame(nx=1, ny=3)
    lam_all = np.asarray(ops.eigen("-fullGenLapack", 8))
    f = np.array([_hz(l) for l in lam_all])
    hi = 0.5 * (f[6] + f[7])          # band holds 7 modes
    lam = np.asarray(ops.eigen("-feast", 0.0, hi, "-m0", 3, "-rci"))
    err = capfd.readouterr().err.lower()
    assert lam.shape == (7,)
    np.testing.assert_allclose(lam, lam_all[:7], rtol=1e-6)
    assert "m0" in err or "subspace" in err or "enlarg" in err


def test_rci_mass_orthonormal():
    """Phi^T M Phi == I through the RCI path — the CSR-scatter- and
    projector-accumulation-sensitive check (a wrong ijob=30/40 matvec or a
    wrong ijob=11 solve shows up here before it shows up in lambda)."""
    _frame()
    m = 5.0e3
    lam_ref = np.asarray(ops.eigen(8))
    f_hi = 0.5 * (_hz(lam_ref[6]) + _hz(lam_ref[7]))
    lam = np.asarray(ops.eigen("-feast", 0.0, f_hi, "-rci"))
    p = lam.shape[0]
    assert p == 7
    nodes = _free_nodes()
    phi = np.array([[c for n in nodes
                     for c in ops.nodeEigenvector(n, k + 1)[:2]]
                    for k in range(p)])
    gram = m * (phi @ phi.T)
    np.testing.assert_allclose(gram, np.eye(p), atol=1e-8)


def test_rci_composes_with_certify(capfd):
    """P3c gate composition: -rci -certify — the Sturm/inertia certificate
    must certify the RCI-found band content."""
    _frame()
    lam_ref = np.asarray(ops.eigen(8))
    f_hi = 0.5 * (_hz(lam_ref[6]) + _hz(lam_ref[7]))
    lam = np.asarray(ops.eigen("-feast", 0.0, f_hi, "-rci", "-certify"))
    assert lam.shape == (7,)
    assert "certified complete" in capfd.readouterr().err


def test_rci_composes_with_blockzgate(capfd):
    """All three diagnostics on one run (-rci -certify -blockZGate):
    eigenpairs unchanged vs the plain driver solve."""
    _frame(nx=1, ny=3)
    lam_plain = np.asarray(ops.eigen("-feast", 0.0, 50.0))
    lam_gated = np.asarray(ops.eigen("-feast", 0.0, 50.0, "-rci",
                                     "-certify", "-blockZGate"))
    np.testing.assert_allclose(lam_gated, lam_plain, rtol=1e-8)
    err = capfd.readouterr().err
    assert "certified complete" in err
    assert "-blockZGate PASS" in err
