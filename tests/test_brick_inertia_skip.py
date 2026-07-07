"""LadrunoBrick CDL inertia no-op skip (ADR-68 T7) — Zone-A gate.

Under CentralDifferenceLadruno the residual inertia pass (`formInertiaTerms(0)`)
runs at trial accel == 0 (the Azero solve): momentum = Sum_j N_j*a_j = 0 exactly,
so the whole block adds only zeros to resid -- yet still pays a geometry
recompute. T7 skips it when every nodal trial accel is exactly 0.0. The mass
matrix (tangFlag==1) is never skipped.

Gate = G-BYTE. The skip must be BIT-IDENTICAL to the un-skipped path
(`-noInertiaSkip`, the A/B escape):

  IS-1  CDL moving-wave bit-identity of the displacement trajectory (the claim);
  IS-2  SIGNBIT lock: the recorded nodal ACCEL is bit-for-bit identical (struct
        bytes), catching any -0.0/+0.0 divergence a value-== gate would miss --
        the adversarial-gate recommendation. Proves the skip is bit-identical,
        not merely ==-identical (resid's Zero()-prime + trailing bodyForce
        normalization means -0.0 never materializes; this test locks it);
  IS-3  implicit passthrough: under Newmark the trial accel is nonzero so the
        skip never fires -- identical to `-noInertiaSkip`;
  IS-4  Rayleigh (-alphaM) under CDL: the mass-proportional damping path
        (getMass re-forms the mass) is bit-identical.

Design: Ladruno_implementation/68_ladruno_state_determination_perf_adr.md (T7)
and 40b_phase0_dominance_report.md (T3 micro-drill addendum).
"""
import struct

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# 3D column of LadrunoBrick, fixed base, elastic (rate-independent so the CDL
# skip fires every step), lumped mass for the explicit Diagonal solve.
NX, NY, NZ = 2, 2, 6
L = 1.0
E, NU, RHO = 1.0e4, 0.25, 1.0
DT = 0.008                 # < dt_cr ~ L/sqrt(E/rho) = 1/100
FZ = -2.0                  # lateral-ish tip load driving a wave


def _nid(i, j, k):
    return 1 + i + (NX + 1) * (j + (NY + 1) * k)


def _build(ele_args=(), rayleigh_alpha=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                ops.node(_nid(i, j, k), i * L, j * L, k * L)
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
                ops.element("LadrunoBrick", e,
                            _nid(i, j, k), _nid(i+1, j, k), _nid(i+1, j+1, k), _nid(i, j+1, k),
                            _nid(i, j, k+1), _nid(i+1, j, k+1), _nid(i+1, j+1, k+1), _nid(i, j+1, k+1),
                            1, "-lumped", *ele_args)
                e += 1
    # top-face load in x (drives lateral vibration)
    top = [_nid(i, j, NZ) for j in range(NY + 1) for i in range(NX + 1)]
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for nd in top:
        ops.load(nd, FZ, 0.0, 0.0)
    if rayleigh_alpha != 0.0:
        ops.rayleigh(rayleigh_alpha, 0.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    return top


def _free_nodes():
    return [_nid(i, j, k) for k in range(1, NZ + 1)
            for j in range(NY + 1) for i in range(NX + 1)]


def _run_cdl(ele_args=(), rayleigh_alpha=0.0, nsteps=300):
    top = _build(ele_args, rayleigh_alpha)
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    tip = top[0]
    traj, accel = [], None
    for s in range(nsteps):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(tip, 1))
    # full free-DOF accel field at the final step (for the signbit lock)
    accel = [ops.nodeAccel(n, d) for n in _free_nodes() for d in (1, 2, 3)]
    return traj, accel


def _run_newmark(ele_args=(), nsteps=20):
    _build(ele_args)
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-12, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    tip = _nid(0, 0, NZ)
    traj = []
    for s in range(nsteps):
        assert ops.analyze(1, 0.02) == 0
        traj.append(ops.nodeDisp(tip, 1))
    return traj


def _assert_bit_identical(a, b, label):
    assert len(a) == len(b)
    bad = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert not bad, (f"{label}: skip vs -noInertiaSkip diverge at {bad[:5]} "
                     f"(first {a[bad[0]]!r} vs {b[bad[0]]!r})")


def _assert_bit_for_bit(a, b, label):
    # struct bytes: distinguishes -0.0 from +0.0 (which == would not)
    assert len(a) == len(b)
    bad = [i for i, (x, y) in enumerate(zip(a, b))
           if struct.pack("<d", x) != struct.pack("<d", y)]
    assert not bad, (f"{label}: byte-level divergence at {bad[:5]} "
                     f"(first {a[bad[0]]!r} vs {b[bad[0]]!r})")


def test_is1_cdl_wave_bit_identity():
    t_skip, _ = _run_cdl(())
    t_noskip, _ = _run_cdl(("-noInertiaSkip",))
    _assert_bit_identical(t_skip, t_noskip, "IS-1 CDL wave")
    assert max(abs(u) for u in t_skip) > 0.0        # the column actually moves


def test_is2_signbit_lock_on_accel():
    _, a_skip = _run_cdl(())
    _, a_noskip = _run_cdl(("-noInertiaSkip",))
    # bit-for-bit (not just ==): locks the resid Zero()+bodyForce normalization
    # invariant the skip's bit-identity depends on.
    _assert_bit_for_bit(a_skip, a_noskip, "IS-2 accel signbit")


def test_is3_implicit_passthrough():
    t_skip = _run_newmark(())
    t_noskip = _run_newmark(("-noInertiaSkip",))
    _assert_bit_identical(t_skip, t_noskip, "IS-3 Newmark passthrough")
    assert max(abs(u) for u in t_skip) > 0.0


def test_is4_rayleigh_alphaM_bit_identity():
    t_skip, _ = _run_cdl((), rayleigh_alpha=0.5)
    t_noskip, _ = _run_cdl(("-noInertiaSkip",), rayleigh_alpha=0.5)
    _assert_bit_identical(t_skip, t_noskip, "IS-4 Rayleigh alphaM")
