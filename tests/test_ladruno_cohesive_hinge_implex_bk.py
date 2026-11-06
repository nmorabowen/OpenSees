"""Zone-A regression (WP-217): the IMPL-EX clock and the B-K mode-mix elastic domain.

Found by the independent (paper/spec-derived) Tier-I oracles of frame-ladruno-validation
(Wave B, B3), which pin each defect with a strict xfail:

1. IMPL-EX was not explicit under DisplacementControl. LadrunoCohesiveHinge,
   LadrunoCohesiveHingeBiaxial, ModIMKPeakOriented and Bilin scaled the extrapolation by
   ops_Dt / dt_n with ops_Dt read live on every trial. In a static analysis ops_Dt is the
   load-factor increment, which under DisplacementControl is part of the solution, so the
   factor changed with every Newton iterate. Repro: Mc = 120, Gf = 1.5, Kpen = 48000
   (kappa0 = 2.5e-3), uniform DisplacementControl steps h = 1.40625e-4: step 18
   (theta = 1.0125 kappa0) returned 121.5 = Kpen theta (undamaged); the published scheme
   (Prazeres et al. 2016 Eq. 17) gives 119.6671. Fix: LadrunoImplexClock (step clock by
   default, -implexDt pseudo latches ops_Dt at the first trial of the step).

2./3. LadrunoCohesiveHingeBiaxial with a B-K mix (eta != 1) evaluated D at the CURRENT mode
   mix also while unloading, so inside the activation surface the secant changed with the
   jump direction: the diagonal unloading tangent was not the Jacobian (central-FD mismatch
   3.5e-4 at (1.5, 3) after a square path), a closed loop strictly inside did non-zero work
   (+1.15e-4 of sum|M.dalpha|) and D could decrease (healing). Fix (ADR 34, 2026-10-29
   amendment): D is the internal variable, frozen with the mix of its last growth; loading
   is D_law(alpha) >= D_c.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _uniform_path(vertices, h):
    """Polyline through `vertices` with every step of length exactly h (segments must be
    integer multiples of h). Returns the node values, starting at vertices[0]."""
    out = [float(vertices[0])]
    for a, b in zip(vertices[:-1], vertices[1:]):
        n = int(round(abs(b - a) / h))
        assert n > 0 and abs(n * h - abs(b - a)) < 1e-12 * max(1.0, abs(b - a))
        out += list(np.linspace(a, b, n + 1)[1:])
    return np.array(out)


# --------------------------------------------------------------------------- #
#  uniaxial drivers: DisplacementControl (zeroLength + unit spring + unit load)
#  and sp on a Path series under LoadControl(1) (ops_Dt = 1 on every step)
# --------------------------------------------------------------------------- #
def _drive_dispctrl(mat_args, theta, spring=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial(*mat_args)
    ops.element("zeroLength", 1, 1, 2, "-mat", mat_args[1], "-dir", 1)
    ops.uniaxialMaterial("Elastic", 99, spring)
    ops.element("zeroLength", 99, 1, 2, "-mat", 99, "-dir", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    out = np.zeros(theta.size)
    u = 0.0
    for i in range(1, theta.size):
        ops.integrator("DisplacementControl", 2, 1, float(theta[i] - u))
        assert ops.analyze(1) == 0, f"DisplacementControl failed at step {i}"
        u = ops.nodeDisp(2, 1)
        out[i] = ops.eleResponse(1, "basicForce")[0]
    ops.wipe()
    return out


def _drive_sp(mat_args, theta, dts=None):
    """sp on a Path series; `dts` (optional) are the LoadControl increments per step."""
    if dts is None:
        dts = np.ones(theta.size - 1)
    times = np.concatenate([[0.0], np.cumsum(dts)])
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial(*mat_args)
    ops.element("zeroLength", 1, 1, 2, "-mat", mat_args[1], "-dir", 1)
    ops.timeSeries("Path", 1, "-time", *[float(t) for t in times], "-values", *[float(x) for x in theta])
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 25, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    out = np.zeros(theta.size)
    for i in range(1, theta.size):
        ops.integrator("LoadControl", float(dts[i - 1]))
        assert ops.analyze(1) == 0, f"sp drive failed at step {i}"
        out[i] = ops.eleResponse(1, "basicForce")[0]
    ops.wipe()
    return out


_MC, _GF, _KPEN = 120.0, 1.5, 48000.0      # kappa0 = Mc/Kpen = 2.5e-3
_H = 1.40625e-4                             # = kappa0 / 17.777...


def _hinge(*extra):
    return ("LadrunoCohesiveHinge", 1, _MC, _GF, "-exp", "-penalty", _KPEN, "-implex", *extra)


def test_hinge_implex_dispctrl_step18_matches_published_scheme():
    """The B3 repro: step 18 crosses kappa0 by extrapolation. Prazeres Eq. 17 with a uniform
    step: kappa~ = kappa_17 + (kappa_17 - kappa_16) = 18 h, M = Menv(kappa~) theta/kappa~."""
    theta = np.arange(0, 25) * _H
    m = _drive_dispctrl(_hinge(), theta)
    k0 = _MC / _KPEN
    a = _MC / (_GF - _MC * _MC / (2.0 * _KPEN))
    expected = _MC * math.exp(-a * (18 * _H - k0))
    assert expected == pytest.approx(119.6671, abs=5e-5)
    assert m[18] == pytest.approx(expected, rel=1e-10), m[18]   # was 121.5 = Kpen*theta


@pytest.mark.parametrize("extra", [(), ("-implexDt", "step")])
def test_hinge_implex_independent_of_integrator(extra):
    """Same uniform cyclic path: DisplacementControl == sp + LoadControl(1). Before WP-217
    the DisplacementControl run used the iterate's load-factor increment as the clock."""
    # nodes at multiples of h = 0.05625 kappa0, so kappa0 falls inside a step (as in the repro)
    theta = _uniform_path([0.0, 40 * _H, 20 * _H, 80 * _H, -80 * _H, 120 * _H], _H)
    m_dc = _drive_dispctrl(_hinge(*extra), theta)
    m_sp = _drive_sp(_hinge(*extra), theta)
    assert np.max(np.abs(m_dc - m_sp)) / _MC < 1.0e-9


def test_hinge_implexdt_pseudo_keeps_the_time_ratio():
    """-implexDt pseudo: f = dt_{n+1}/dt_n from LoadControl, latched per step. With uniform
    dlambda it equals the step clock; when dlambda halves, the extrapolated kappa~ uses the
    ratio 1/2 (the step clock uses 1), which is visible on the first softening step."""
    k0 = _MC / _KPEN
    theta = _uniform_path([0.0, 1.0 * k0], 0.05 * k0)
    th = np.concatenate([theta, theta[-1] + 0.025 * k0 * np.arange(1, 9)])
    dts = np.concatenate([np.ones(theta.size - 1), 0.5 * np.ones(8)])
    m_step = _drive_sp(_hinge("-implexDt", "step"), th, dts)
    m_pseudo = _drive_sp(_hinge("-implexDt", "pseudo"), th, dts)
    m_uni = _drive_sp(_hinge("-implexDt", "pseudo"), th)
    assert np.max(np.abs(m_uni - _drive_sp(_hinge(), th))) == 0.0
    # the first half step (dlambda 1 -> 0.5) is the first one past kappa0:
    # kappa~ = 1.05 k0 (step clock) vs 1.025 k0 (pseudo, f = 0.5)
    i = theta.size                           # theta = 1.025 k0
    kn, knm1 = th[i - 1], th[i - 2]
    a = _MC / (_GF - _MC * _MC / (2.0 * _KPEN))
    for m, f in ((m_step, 1.0), (m_pseudo, 0.5)):
        kt = kn + f * (kn - knm1)
        env = _KPEN * kt if kt <= k0 else _MC * math.exp(-a * (kt - k0))
        assert m[i] == pytest.approx(env / kt * th[i], rel=1e-10), (f, m[i])
    assert m_step[i] != m_pseudo[i]


_BILIN_ARGS = (2.0e7, 0.05, 0.05, 100.0e3, -100.0e3,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
               0.01, 0.01, 0.02, 0.02, 0.3, 0.3, 0.4, 0.4, 1.0, 1.0)


@pytest.mark.parametrize("name", ["Bilin", "ModIMKPeakOriented"])
def test_imk_bilin_implex_independent_of_integrator(name):
    """The same per-step clock in Bilin and ModIMKPeakOriented -implex: a cyclic path through
    yield, capping and reversal gives the same moments under both drivers."""
    theta = _uniform_path([0.0, 0.012, -0.004, 0.02, -0.02], 2.0e-4)
    args = (name, 1, *_BILIN_ARGS, "-implex")
    m_dc = _drive_dispctrl(args, theta, spring=1.0e3)
    m_sp = _drive_sp(args, theta)
    assert np.max(np.abs(m_dc - m_sp)) / 100.0e3 < 1.0e-9


def test_implexdt_rejects_unknown_clock():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    with pytest.raises(Exception):
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 1, _MC, _GF, "-implex", "-implexDt", "wall")
    ops.wipe()


# --------------------------------------------------------------------------- #
#  biaxial material point: node-2 translations ARE alpha, each on its own Path
# --------------------------------------------------------------------------- #
_BMAT = (120.0, 1.5, 80.0, 1.0)            # Mcz, Gfz, Mcy, Gfy (the 00b slice's values)


def _a0(ratio=1000.0):
    Mcz, Gfz, Mcy, Gfy = _BMAT
    return np.array([Mcz / (ratio * Mcz * Mcz / (2 * Gfz)), Mcy / (ratio * Mcy * Mcy / (2 * Gfy))])


def _norm_path(vertices, h):
    """Polyline in normalized e = alpha/alpha0 with step h, returned in alpha units."""
    pts = [np.asarray(vertices[0], float)]
    for a, b in zip(vertices[:-1], vertices[1:]):
        a, b = np.asarray(a, float), np.asarray(b, float)
        n = max(1, int(round(np.linalg.norm(b - a) / h)))
        pts += [a + (b - a) * k / n for k in range(1, n + 1)]
    return np.array(pts) * _a0()


def _drive_biaxial(P, eta, implex=False, read_tangent=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 0, 1)
    args = ["LadrunoCohesiveHingeBiaxial", 7, *_BMAT, "-exp"]
    if eta != 1.0:
        args += ["-bk", eta]
    if implex:
        args += ["-implex"]
    ops.nDMaterial(*args)
    ops.element("zeroLengthND", 1, 1, 2, 7)
    t = [float(i) for i in range(P.shape[0])]
    for j in (0, 1):
        ops.timeSeries("Path", j + 1, "-time", *t, "-values", *[float(x) for x in P[:, j]])
        ops.pattern("Plain", j + 1, j + 1)
        ops.sp(2, j + 1, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    M = np.zeros_like(P)
    D = np.zeros(P.shape[0])
    for i in range(1, P.shape[0]):
        assert ops.analyze(1) == 0, f"biaxial drive failed at step {i}"
        M[i] = [ops.eleResponse(1, "material", "momentZ")[0], ops.eleResponse(1, "material", "momentY")[0]]
        D[i] = ops.eleResponse(1, "material", "damage")[0]
    K = np.array(ops.eleResponse(1, "material", "tangent")).reshape(2, 2) if read_tangent else None
    ops.wipe()
    return M, D, K


_SQUARE_TO_INTERIOR = [[0, 0], [3, 0], [3, 3], [1.5, 3]]


@pytest.mark.parametrize("eta", [1.0, 1.5])
def test_biaxial_unloading_tangent_is_the_jacobian(eta):
    """After the square path the point (1.5, 3) is unloaded and off-axis. The fork's tangent
    there must equal the central FD of its own M(alpha). eta = 1.5 was 3.5e-4 off."""
    P = _norm_path(_SQUARE_TO_INTERIOR, 0.05)
    _, _, K = _drive_biaxial(P, eta, read_tangent=True)
    a, h = P[-1], 1.0e-6 * _a0()
    Kfd = np.zeros((2, 2))
    for j in (0, 1):
        e = np.zeros(2)
        e[j] = h[j]
        Mp, _, _ = _drive_biaxial(np.vstack([P, a + e]), eta)
        Mm, _, _ = _drive_biaxial(np.vstack([P, a - e]), eta)
        Kfd[:, j] = (Mp[-1] - Mm[-1]) / (2.0 * h[j])
    err = np.max(np.abs(K - Kfd)) / np.max(np.abs(K))
    assert err < 1.0e-6, (K, Kfd, err)


@pytest.mark.parametrize("eta", [1.0, 1.5, 0.5])
def test_biaxial_interior_closed_loop_is_conservative(eta):
    """Load to r_max = sqrt(18), unload to (1.5, 3), run a closed loop strictly inside:
    the loop work vanishes (M is linear in alpha there) and D is constant on it.
    eta = 1.5 used to give +1.15e-4 of sum|M.dalpha|."""
    pre = _norm_path(_SQUARE_TO_INTERIOR, 0.01)
    loop = _norm_path([[1.5, 3], [0.2, 3], [0.2, 0.5], [1.5, 0.5], [1.5, 3]], 0.01)[1:]
    P = np.vstack([pre, loop])
    M, D, _ = _drive_biaxial(P, eta)
    i0 = pre.shape[0] - 1
    dW = np.sum(0.5 * (M[i0 + 1:] + M[i0:-1]) * np.diff(P[i0:], axis=0), axis=1)
    assert abs(np.sum(dW)) / np.sum(np.abs(dW)) < 1.0e-9
    assert np.max(np.abs(D[i0:] - D[i0])) < 1.0e-15   # eta = 1 recomputes D at the current mix: round-off


@pytest.mark.parametrize("eta", [1.5, 0.5])
def test_biaxial_damage_never_heals_and_stress_is_continuous(eta):
    """Full square path (off-radial loading, unloading through a pure axis, re-entry at the
    opposite corner): D is non-decreasing and M has no jump larger than the step's elastic
    increment (the damage criterion D_law >= D_c is continuous by construction)."""
    h = 0.02
    P = _norm_path([[0, 0], [3, 0], [3, 3], [-3, 3], [-3, -3], [0, 0]], h)
    M, D, _ = _drive_biaxial(P, eta)
    assert np.all(np.diff(D) >= 0.0), np.min(np.diff(D))
    Kp = np.array(_BMAT[0::2]) / _a0()          # Kpen_i = Mc_i / a0_i
    jump = np.max(np.abs(np.diff(M, axis=0)) / (Kp * h * _a0()))
    assert jump < 1.0 + 1.0e-9, jump


def test_biaxial_eta1_frozen_equals_current_mix():
    """eta = 1: D is mix-independent, so the amended law is the ADR 34 law to round-off."""
    P = _norm_path([[0, 0], [3, 0], [3, 3], [-3, 3], [1.5, 0.5]], 0.05)
    M, _, _ = _drive_biaxial(P, 1.0)
    a0 = _a0()
    Mcz, Gfz, Mcy, Gfy = _BMAT
    Kz, Ky = Mcz / a0[0], Mcy / a0[1]
    rmax = 0.0
    for i in range(1, P.shape[0]):
        e = P[i] / a0
        r = math.hypot(*e)
        rmax = max(rmax, r)
        wz, wy = (e[0] / r) ** 2, (e[1] / r) ** 2
        S = wz * Mcz ** 2 / Kz + wy * Mcy ** 2 / Ky
        Es = wz * Gfz + wy * Gfy - 0.5 * S
        Dd = 0.0 if rmax <= 1.0 else 1.0 - math.exp(-(S / Es) * (rmax - 1.0)) / rmax
        ref = (1.0 - Dd) * np.array([Kz * P[i, 0], Ky * P[i, 1]])
        assert np.max(np.abs(M[i] - ref) / np.array([Mcz, Mcy])) < 1.0e-12, i


def test_biaxial_implex_independent_of_integrator_on_axis():
    """IMPL-EX clock in the biaxial law: a pure-z cyclic path under DisplacementControl equals
    the sp + LoadControl(1) drive (and the scalar hinge on that axis)."""
    a0z = _a0()[0]
    th = _uniform_path([0.0, 3.0, 1.0, 6.0, -6.0], 0.05) * a0z
    P = np.column_stack([th, np.zeros_like(th)])
    M_sp, _, _ = _drive_biaxial(P, 1.5, implex=True)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 7, *_BMAT, "-exp", "-bk", 1.5, "-implex")
    ops.element("zeroLengthND", 1, 1, 2, 7)
    ops.uniaxialMaterial("Elastic", 99, 1.0)
    ops.element("zeroLength", 2, 1, 2, "-mat", 99, "-dir", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    u = 0.0
    Mz = np.zeros(th.size)
    for i in range(1, th.size):
        ops.integrator("DisplacementControl", 2, 1, float(th[i] - u))
        assert ops.analyze(1) == 0, i
        u = ops.nodeDisp(2, 1)
        Mz[i] = ops.eleResponse(1, "material", "momentZ")[0]
    ops.wipe()
    assert np.max(np.abs(Mz - M_sp[:, 0])) / _BMAT[0] < 1.0e-9
