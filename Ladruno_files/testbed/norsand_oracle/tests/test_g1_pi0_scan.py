"""WP-144 round 3b, G1 (Zone B): the unified pi_i0 rule (S.53) and the gated scan-step refusal (S.56).  Equation sheet 144a 5.4, 10.2, 13.15, 16.3.

TEST RULE (plan 5): every expected value comes from the sheet -- the closed form (S.53) / (S.56) coded in g1_hf_common.py, or the
number PRINTED in 13.15 / 10.2 -- and was fixed before the oracles ran on these paths.  A failing assertion is a finding.

Gates
  K1.15   pi_i0 on the K2 set at p_init = -100: eta* = c2 M = 0.18 -> -50.995881 (smooth cap c2 = 0.15, and the planar cap c1 = c2 = 0.15),
          eta* = 0 -> -46.475800 (no cap: the apex), eta_init = 0.75 -> -71.554175, eta* = M -> -100 (the image point).  Both oracles
          (O2 `initial_state` default = the unified rule; O1 `pi_rule="S53"`).  Closed-form identities: eta(p_init, pi_i0) = eta* in both
          branches (N > 0 and N = 0), N -> 0 limit O(N), F = |p|(eta_init - c2 M) < 0 when the state is inside, F = 0 when it is not, the
          lode dependence eta_init = zeta(theta) q/|p|, the floored p_init (the floor acts first), the refusal eta* >= M/N, `-pi0` overrides
  yield   the consequence printed in 13.15: first yield on a drained TXC path from the isotropic state with pi_i0 = -50.9959 at q =
          11.3469, p = -103.7823, eta = 0.1093, w = 0.337 (cap partly active); on a constant-p path at eta = 0.18 exactly
  S.56    the gated validate() refusal PI_SCAN_REL <= W_ramp / 10, cap = smooth ONLY (round 3b A2, A3): W_ramp = 0.0606 on the K2 defaults;
          c1 = 0.05 with c2 = 0.07 (0.0122) is admissible, c2 = 0.06 (0.0061) is refused; planar (c1 = c2) and no cap are never refused

MUTANTS (named for the mutation gate, plan 5.3).
  pi_i0 rule:   the pre-round-3 apex default kept for a capped model (K1.15 eta* = c2 M); c2 taken as c1 or 1 - c2; eta_init computed without zeta
                (theta ignored, TXE/off-corner case); the floor applied after the rule instead of before (floored-p case); the N = 0 branch
                evaluated with the power formula (division by zero) or the power branch used for N = 0; eta* = max replaced by min or by a sum
                (the image-point and eta_init = 0.75 cases); no refusal at eta* >= M/N; an explicit -pi0 ignored.
  S.56 gate:    refusal absent; refusal ungated (planar / none refused: W_ramp = 0 there); W_ramp with the round-3 inverted ratio (-0.0645:
                every smooth cap refused, including the K2 default); the factor 10 dropped or inverted; PI_SCAN_REL changed (a contract constant).
"""
import math
import warnings

import numpy as np
import pytest

import g1_hf_common as C
from conftest import ORACLES

O1, O2 = ORACLES["O1"], ORACLES["O2"]
I3 = np.eye(3)
M_, N_ = 1.2, 0.4


def init_state(oname, kw, sigma0, v0=1.70, pi_i0=None):
    """The unified rule: O2's default; O1 with pi_rule = 'S53'.  pi_i0 given -> the explicit override."""
    P = C.make_hf(oname, **kw)
    if oname == "O1":
        return P, O1.initial_state(P, sigma0, v0, pi_i0, pi_rule="S53")
    return P, O2.initial_state(P, sigma0, v0, pi_i0)


def k2(cap="smooth", **over):
    kw = C.ba06_kw(p_min=0.0, cap=cap)
    if cap == "smooth":
        kw.update(c1=0.05, c2=0.15)
    elif cap == "planar":
        kw.update(c1=0.15, c2=0.15)
    kw.update(over)
    return kw


def txc(p, eta):
    """Compression-meridian (theta = pi/3, zeta = 1) principal stress with q = eta |p|."""
    return np.diag(C.sigma_pq(p, eta * abs(p), np.array([1.0, 1.0, -2.0])))


# ======================================================================================================================
# K1.15  printed values and the closed form
# ======================================================================================================================
@pytest.mark.parametrize("cap,sigma_eta,printed", [
    ("smooth", 0.0, -50.995881),            # eta* = c2 M = 0.18
    ("planar", 0.0, -50.995881),            # planar: c2 := c1 = chi_cap = 0.15
    ("none", 0.0, -46.475800),              # eta* = 0: the apex through p_init
    ("smooth", 0.75, -71.554175),           # eta_init = 0.75 > c2 M
    ("none", 0.75, -71.554175),
    ("planar", 0.75, -71.554175),
    ("none", 1.2, -100.0),                  # eta* = M: the image point, pi_i0 = p_init
    ("smooth", 1.2, -100.0)])
def test_k1_15_pi_i0_unified_rule_printed_values_both_oracles(oracle_name, cap, sigma_eta, printed):
    """13.15, K2 set (M = 1.2, N = 0.4, c1 = 0.05, c2 = 0.15 smooth; planar c1 = c2 = 0.15), p_init = -100: eta* = max(eta_init, c2 M).  Isotropic
    start (eta_init = 0): smooth / planar -> eta* = 0.18, pi_i0 = -50.995881; no cap -> eta* = 0 -> -46.475800 (the apex, the G1 value of
    9.6 (F)).  Compression-meridian start with eta_init = 0.75 (zeta(pi/3) = 1, q = 75): -71.554175 for every cap.  eta_init = M: pi_i0 =
    p_init = -100.  The closed form (S.53) of the common module reproduces the printed digits to 6e-7 (sheet check inside the test) and
    both oracles are compared to the PRINTED values (1e-8 relative: 8 printed digits) and to the closed form (1e-13).
    KILLS: the pre-round-3 apex default kept for a capped model (-46.4758 instead of -50.9959), c2 mapped to c1 / (1 - c2), eta* = min or a sum
    of eta_init and c2 M, the planar cap's c2 taken as 0, the image-point case off by the (1-N) power."""
    kw = k2(cap)
    sig = -100.0 * I3 if sigma_eta == 0.0 else txc(-100.0, sigma_eta)
    c2 = {"none": 0.0, "planar": 0.15, "smooth": 0.15}[cap]
    eta_star = max(sigma_eta, c2 * M_)
    cf = C.pi_S53(M_, N_, -100.0, eta_star)
    assert abs(cf - printed) <= 6e-7 * abs(printed) + 1e-6, (cf, printed)
    P, st = init_state(oracle_name, kw, sig)
    assert abs(st.pi_i - printed) <= 1e-8 * abs(printed) + 6e-7, (st.pi_i, printed)
    assert abs(st.pi_i - cf) <= 1e-12 * abs(cf), (st.pi_i, cf)
    # eta(p_init, pi_i0) = eta* exactly (the inverse of S.12)
    assert abs(C.eta_S12(M_, N_, -100.0, st.pi_i) - eta_star) <= 1e-11


@pytest.mark.parametrize("N", [0.4, 0.2, 0.0])
def test_pi_i0_inverse_identity_both_branches_and_f_at_the_initial_stress(oracle_name, N):
    """(S.53) / 5.4: for every N the surface through (p_init, eta*) has eta(p_init, pi_i0) = eta* (power branch for N > 0, the exp branch
    for N = 0), and the yield function at the initial stress is F = |p|(eta_init - c2 M) < 0 when eta_init < c2 M (the state is INSIDE:
    the whole cap ramp lies inside the initial surface), F = 0 when eta_init >= c2 M (the surface passes through the state: the K0 deck
    states).  Three starts: isotropic, an off-corner state with eta_init = 0.09 (< 0.18) and one with eta_init = 0.75 at theta = 0.5.
    rho = rho_bar = 0.7, N_bar = N, smooth cap (0.05, 0.15), seeded v0 with B > 0.
    KILLS: the N = 0 branch evaluated with the power formula; eta* not entering the identity; F evaluated against a surface that does not
    pass through the state although eta_init >= c2 M; eta_init computed without zeta (the F identity breaks at theta != pi/3)."""
    kw = k2("smooth", N=N, N_bar=N, rho=0.7, rho_bar=0.7)
    for th, eta_i in ((None, 0.0), (0.5, 0.09), (0.5, 0.75)):
        if th is None:
            sig = -100.0 * I3
        else:
            nh = C.dir_for_theta(th)
            z = C.zeta_ww(th, 0.7)
            sig = np.diag(C.sigma_pq(-100.0, eta_i * 100.0 / z, nh))
        P, st = init_state(oracle_name, kw, sig)
        eta_star = max(eta_i, 0.15 * M_)
        assert abs(st.pi_i - C.pi_S53(M_, N, -100.0, eta_star)) <= 1e-12 * abs(st.pi_i)
        assert abs(C.eta_S12(M_, N, -100.0, st.pi_i) - eta_star) <= 1e-11
        F = C.F_yield(M_, N, 0.7, np.diag(sig), st.pi_i)
        if eta_i >= 0.15 * M_:
            assert abs(F) <= 1e-10 * 100.0, F
        else:
            assert abs(F - 100.0 * (eta_i - 0.15 * M_)) <= 1e-10 * 100.0 and F < 0.0


def test_pi_i0_n_to_zero_limit_is_continuous_order_n():
    """5.4: 'the N -> 0 limit of the power branch is the exp branch (gap O(N): 2.4e-10 at N = 1e-9)'.  The printed 2.4e-10 is a 30-digit
    evaluation at an unstated set (on the K2 numbers, eta* = 0.18, the gap is 0.36 N); in double precision the power branch at N = 1e-9 loses
    ~1e-7 to the exponent 1/N = 1e9, so the gate is the O(N) LAW, where round-off is 1e-11: gap(1e-5)/gap(1e-4) = 0.1 within 5 % and gap/N in
    [0.1, 1] (order one), from the common closed form AND from O2's initial_state (N_bar = 0, rho = rho_bar = 0.7), the latter agreeing with the
    closed form to 1e-9.  KILLS: a discontinuous branch switch (a different formula at N = 0), an O(1) gap."""
    cf0 = C.pi_S53(M_, 0.0, -100.0, 0.18)
    gaps = {}
    for N in (1e-4, 1e-5):
        g_cf = abs(C.pi_S53(M_, N, -100.0, 0.18) / cf0 - 1.0)
        assert 0.1 <= g_cf / N <= 1.0, (N, g_cf / N)
        gaps[N] = g_cf
        out = []
        for NN in (0.0, N):
            P = C.make_hf("O2", **k2("smooth", N=NN, N_bar=0.0, rho=0.7, rho_bar=0.7))
            out.append(O2.initial_state(P, -100.0 * I3, 1.70, None).pi_i)
        assert abs(out[0] - cf0) <= 1e-12 * abs(cf0)
        assert abs(out[1] / out[0] - C.pi_S53(M_, N, -100.0, 0.18) / cf0) <= 1e-9, out
    assert abs(gaps[1e-5] / gaps[1e-4] - 0.1) <= 5e-3


def test_pi_i0_refuses_eta_star_at_or_above_M_over_N(oracle_name):
    """(S.53): refuse at initialState if eta* >= M/N (N > 0): no surface passes through the state.  K2: M/N = 3.  A compression-meridian
    start at p = -100 with eta_init = 3.3 and 3.0 is refused (ValueError); eta_init = 2.0 is accepted and has pi_i0 = p [(1-N)/(1 - 2N/M)]^((1-N)/N)
    (closed form).  KILLS: no refusal (a NaN / complex pi_i0 or the apex rule applied), the refusal at > instead of >= (eta* = M/N exactly)."""
    kw = k2("none")
    for eta in (3.3, 3.0):
        with pytest.raises(ValueError):
            init_state(oracle_name, kw, txc(-100.0, eta))
    P, st = init_state(oracle_name, kw, txc(-100.0, 2.0), v0=1.60)
    assert abs(st.pi_i - C.pi_S53(M_, N_, -100.0, 2.0)) <= 1e-12 * abs(st.pi_i)


def test_pi_i0_explicit_value_overrides_the_rule(oracle_name):
    """16.3: 'an explicit -pi0 overrides' the rule (the paper-mode K2 benchmark needs pi_i0 = p_c 0.6^1.5, 14).  pi_i0 = -80 is used as
    given for a smooth-cap model whose rule would give -50.995881.  KILLS: an override ignored."""
    P, st = init_state(oracle_name, k2("smooth"), -100.0 * I3, pi_i0=-80.0)
    assert st.pi_i == -80.0


def test_pi_i0_uses_the_floored_initial_pressure(oracle_name):
    """5.4 / 9.7: 'the pi_i0 rule uses the floored p_init'.  BA06 K2, p_min = 0.5: sigma0 = isotropic -0.2 floors to p_init = -0.5 (counted), and
    the rule is applied THEN: no cap -> pi_i0 = p_init (1-N)^((1-N)/N) = -0.5 x 0.464758 = -0.232379 (the apex through the floored point);
    smooth cap -> pi_i0 = -0.5 x 0.50995881 = -0.2549794 (eta* = c2 M scales with p_init).  KILLS: the rule applied to the deck's raw p =
    -0.2 (the order floor -> rule reversed)."""
    for cap, ref in (("none", 0.464758001544890), ("smooth", 0.50995881)):
        kw = k2(cap, p_min=0.5)
        P, st = init_state(oracle_name, kw, -0.2 * I3)
        pm = float(np.trace(st.sigma)) / 3.0
        assert abs(pm + 0.5) <= 1e-12
        assert abs(st.pi_i - (-0.5 * ref)) <= (1e-7 if cap == "smooth" else 1e-12) * 0.5 * ref


def test_pi_i0_lode_dependence_eta_init_uses_zeta_of_theta(oracle_name):
    """5.4: eta_init := zeta(theta, rho) q/|p| (0 on the axis).  Two states with the SAME q = 40 at p = -100 on the extension (theta = 0,
    zeta = 1/rho = 1/0.7) and the compression (theta = pi/3, zeta = 1) meridians give eta_init = 0.5714 and 0.40 and therefore different
    pi_i0 (no cap): pi_S53(p, 0.5714) and pi_S53(p, 0.40).  KILLS: eta_init = q/|p| without zeta (the two pi_i0 coincide)."""
    kw = k2("none", rho=0.7, rho_bar=0.8)
    out = []
    for nh, th in ((np.array([-1.0, -1.0, 2.0]), 0.0), (np.array([1.0, 1.0, -2.0]), math.pi / 3)):
        sig = np.diag(C.sigma_pq(-100.0, 40.0, nh))
        P, st = init_state(oracle_name, kw, sig)
        out.append(st.pi_i)
        z = C.zeta_ww(th, 0.7)
        assert abs(C.theta_of_dir(nh) - th) <= 1e-9
        assert abs(st.pi_i - C.pi_S53(M_, N_, -100.0, z * 0.4)) <= 1e-11 * abs(st.pi_i), (th, st.pi_i)
    assert abs(out[0] - out[1]) > 1.0, "extension and compression starts gave the same pi_i0: zeta ignored"


# ======================================================================================================================
# the consequence printed in 13.15: where the first yield happens
# ======================================================================================================================
def smooth_w(eta, c1=0.05, c2=0.15, M=M_):
    """(S.35) quintic blend weight."""
    t = min(1.0, max(0.0, (eta - c1 * M) / ((c2 - c1) * M)))
    return t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)


def test_first_yield_on_a_drained_txc_path_from_the_unified_start():
    """13.15 consequence: from the isotropic start p = -100 with pi_i0 = -50.9959 (smooth cap 0.05 / 0.15) the response is elastic up to a
    finite q: first yield at q = 11.3469 kPa, p = -103.7823 kPa, eta = 0.1093 (inside the ramp: w = 0.337), not plastic from the first
    increment at w = 0 (the apex start).  O2 `triaxial` (drained), 3000 equal axial increments to -1.5e-3: the first plastic state brackets
    the printed point by one increment (dq = 7e-3 kPa): q within 2e-2, p within 1e-2, eta within 3e-4, w = S(t) of (S.35) within 3e-3.
    CONTROL: the apex start (no cap, pi_i0 = -46.4758) is plastic from the FIRST increment.
    KILLS: the pre-round-3 apex start under a smooth cap (plastic from increment 1 at w = 0, no ramp), c2 taken wrong (the first-yield q moves)."""
    P = C.make_hf("O2", **k2("smooth"))
    st0 = O2.initial_state(P, -100.0 * I3, 1.70, None)
    assert abs(st0.pi_i + 50.995881) <= 1e-5
    sts = O2.triaxial(P, st0, "drained", -1.5e-3, 3000)
    first = next(i for i, s in enumerate(sts) if s.flags["plastic"])
    assert first >= 1 and not any(s.flags["refused"] for s in sts[:first + 1])
    s = sts[first]
    p, q, _ = C.invariants_pq(s.sigma)
    eta = C.eta_S12(M_, N_, p, s.pi_i)
    assert abs(q - 11.3469) <= 2e-2 and abs(p + 103.7823) <= 1e-2, (q, p)
    assert abs(eta - 0.1093) <= 3e-4 and abs(smooth_w(eta) - 0.337) <= 3e-3, (eta, smooth_w(eta))
    # control: the apex start is plastic from the first increment
    P2 = C.make_hf("O2", **k2("none"))
    a0 = O2.initial_state(P2, -100.0 * I3, 1.70, None)
    assert abs(a0.pi_i + 46.4758) <= 1e-4
    a1 = O2.triaxial(P2, a0, "drained", -1.5e-3, 3000)
    assert a1[0].flags["plastic"]


def test_first_yield_on_a_constant_p_path_is_at_eta_c2_M_exactly():
    """13.15: 'on a constant-p path eta = 0.18 exactly' (w = 1, the cap inactive at the first yield).  A pure shear path (tr = 0, p0 constant
    under alpha0 = 0) along the compression meridian in steps of eps_s = 1e-5 (dq = 0.162 kPa, d eta = 1.62e-3): the last elastic state has
    eta <= 0.18 (1e-9), the first plastic state eta in (0.18, 0.18 + 1.7e-3], p = -100 throughout (1e-9).  KILLS: a surface that does not pass
    through (p, c2 M) at p_init (first yield elsewhere), c2 / M mixed up."""
    P = C.make_hf("O2", **k2("smooth"))
    st = O2.initial_state(P, -100.0 * I3, 1.70, None)
    d = np.diag(C.SQ32 * 1e-5 * np.array([-1.0, -1.0, 2.0]) / math.sqrt(6.0)) * -1.0       # compression meridian, tr = 0
    assert abs(np.trace(d)) <= 1e-18
    last_eta = None
    for k in range(200):
        nxt = O2.step(P, st, d)
        assert not nxt.flags["refused"]
        p, q, _ = C.invariants_pq(nxt.sigma)
        eta = q / abs(p)
        if not nxt.flags["plastic"]:
            assert abs(p + 100.0) <= 1e-9, "p must stay constant on the elastic part of the shear path (alpha0 = 0)"
        if nxt.flags["plastic"]:
            assert last_eta is not None and last_eta <= 0.18 + 1e-9 and 0.18 < eta <= 0.18 + 1.7e-3, (last_eta, eta)
            return
        last_eta, st = eta, nxt
    raise AssertionError("no yield within 200 shear increments")


# ======================================================================================================================
# S.56  the gated scan-step refusal
# ======================================================================================================================
def _validated(oname, **kw):
    return C.make_hf(oname, **kw).validate()


@pytest.mark.sheet_check
def test_w_ramp_closed_form_printed_values_and_sign():
    """(S.56) / round 3b A2: W_ramp = 1 - pi_i(eta_1)/pi_i(eta_2) = 1 - [(1 - c2 N)/(1 - c1 N)]^((1-N)/N) (N > 0), 1 - exp(-(c2 - c1)) (N = 0): 0.0606
    on the K2 defaults (N = 0.4, c1 = 0.05, c2 = 0.15), 0.0122 for (0.05, 0.07), 0.0061 for (0.05, 0.06); POSITIVE (the round-3 text's inverted ratio
    is -0.0645 on K2).  Also by the definition at fixed p: the ratio of the two (S.53) image pressures at eta_1 and eta_2.  Sheet check."""
    assert abs(C.w_ramp(0.4, 0.05, 0.15) - 0.0606) <= 5e-5
    assert abs(C.w_ramp(0.4, 0.05, 0.07) - 0.0122) <= 5e-5 and abs(C.w_ramp(0.4, 0.05, 0.06) - 0.0061) <= 5e-5
    for N, c1, c2 in ((0.4, 0.05, 0.15), (0.2, 0.05, 0.1), (0.0, 0.05, 0.15)):
        by_def = 1.0 - C.pi_S53(M_, N, -100.0, c1 * M_) / C.pi_S53(M_, N, -100.0, c2 * M_)
        assert abs(C.w_ramp(N, c1, c2) - by_def) <= 1e-13 and by_def > 0.0
    inverted = 1.0 - C.pi_S53(M_, 0.4, -100.0, 0.15 * M_) / C.pi_S53(M_, 0.4, -100.0, 0.05 * M_)
    assert abs(inverted - (-0.0645)) <= 5e-5                           # the withdrawn round-3 words


@pytest.mark.parametrize("N", [0.4, 0.0])
def test_validate_gate_smooth_cap_refused_when_narrower_than_ten_scan_steps(oracle_name, N):
    """(S.56): PI_SCAN_REL (1e-3) <= W_ramp / 10, cap = smooth ONLY.  N = 0.4: (0.05, 0.15) W = 0.0606 accepted; (0.05, 0.07) W = 0.0122 accepted;
    (0.05, 0.06) W = 0.0061 REFUSED (W/10 = 6.1e-4 < 1e-3).  N = 0: W = 1 - exp(-(c2 - c1)): (0.05, 0.07) 0.0198 accepted, (0.05, 0.055) 0.00499 refused.
    The Params W_ramp property is the closed form.  Boundary: the narrowest admissible c2 is where W_ramp = 0.01 (one scan step per 10): c2 just
    above it accepted, just below refused (both energies' common Params).
    KILLS: refusal absent; the factor 10 dropped (0.0061 accepted at W >= 1e-3) or inverted; W_ramp with the inverted ratio (-0.0645: even the K2
    default refused); PI_SCAN_REL changed."""
    base = dict(N=N, N_bar=N, rho=0.7, rho_bar=0.7) if N == 0.0 else {}
    ok = [(0.05, 0.15), (0.05, 0.07)] if N > 0 else [(0.05, 0.15), (0.05, 0.07)]
    bad = [(0.05, 0.06)] if N > 0 else [(0.05, 0.055)]
    for c1, c2 in ok:
        P = _validated(oracle_name, **k2("smooth", c1=c1, c2=c2, **base))
        assert abs(P.W_ramp - C.w_ramp(N, c1, c2)) <= 1e-13
        assert C.w_ramp(N, c1, c2) / 10.0 >= 1e-3
    for c1, c2 in bad:
        assert C.w_ramp(N, c1, c2) / 10.0 < 1e-3
        with pytest.raises(ValueError):
            _validated(oracle_name, **k2("smooth", c1=c1, c2=c2, **base))
    # the boundary W_ramp = 0.01: bisection on c2 with the common closed form
    lo, hi = 0.051, 0.15
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if C.w_ramp(N, 0.05, mid) < 0.01:
            lo = mid
        else:
            hi = mid
    c2_edge = hi
    _validated(oracle_name, **k2("smooth", c1=0.05, c2=c2_edge + 1e-6, **base))
    with pytest.raises(ValueError):
        _validated(oracle_name, **k2("smooth", c1=0.05, c2=c2_edge - 1e-6, **base))


@pytest.mark.parametrize("cap", ["planar", "none"])
def test_validate_gate_never_refuses_planar_or_no_cap(oracle_name, cap):
    """(S.56) round 3b A3: planar (c1 = c2, W_ramp == 0) and no cap have no ramp and are NOT subject to the refusal -- an ungated check would
    reject every planar and no-cap model.  Planar chi_cap in {0.01, 0.10 (BA06's value), 0.15} and the no-cap default are accepted (and W_ramp is
    0).  KILLS: the ungated refusal (every planar / none model refused: W_ramp = 0 <  10 PI_SCAN_REL)."""
    if cap == "none":
        P = _validated(oracle_name, **k2("none"))
        assert P.W_ramp == 0.0
    else:
        for chi in (0.01, 0.10, 0.15):
            P = _validated(oracle_name, **k2("planar", c1=chi, c2=chi))
            assert P.W_ramp == 0.0


def test_scan_contract_constants_are_the_sheets():
    """(S.56) / 10.2: PI_SCAN_REL = 1e-3 and PI_SCAN_MAX = 1000 (the nested scan travels at most |pi_i,n|); the withdrawn G1 value is 1e-4.
    The refusal in both oracles uses the same constant.  KILLS: a constant changed without the sheet (the 10 x contract moves)."""
    assert O2.kernel.PI_SCAN_REL == 1.0e-3 and O2.kernel.PI_SCAN_MAX == 1000
    from o1_rate import params as P1
    assert P1.PI_SCAN_REL == 1.0e-3
