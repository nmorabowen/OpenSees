"""WP-144 round 3b, G1 (Zone B): the p' floor Pi_f (S.48) on both oracles.  Equation sheet 144a 9.7, 13.12 - 13.14b, 10.2 (S.56 chains), 16.3.

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form of the sheet (coded in
g1_hf_common.py from the sheet's formulas, section cited at each use), (b) a number PRINTED in the sheet (compared to the printed
digits), or (c) a stated convergence / truncation argument; all of it was written BEFORE the oracles were run on these paths.
No oracle output is an expected value.  A failing assertion is a finding (against the sheet or an oracle), never a reason to
edit an oracle or loosen a bound.  Owner decisions in force (2026-10-03): the floor is the strain-space projection Pi_f (S.48),
applied at the trial and after convergence, NEVER inside the local Newton, counted (W_f), never refuses; the kernel keeps the EXACT
consistent tangent (zero bulk stiffness at a floored state, no regularisation); E_f is an on-demand diagnostic from the closed-form Psi.

Gates (ids follow the sheet's K1 numbering)
  K1.12   BA06 floor: eps_v,f = 0.0529831737, a trial at p = -p_min/2 (eps_v 0.0599146455) projects to p = -0.5 with
          d eps^f_v = kappa ln 2, W_f = 3.465735903e-3, E_f = 2.5e-3; q, n_hat, eps_s, pi_i, v unchanged; the floored tangent is
          2 mu0 (I_sym - delta delta / 3) with delta:C = 0 (full 4th-order tensor, non-coaxial state); idempotence; p_min = 0 off
  K1.13   HAR floor and the out-of-domain trial: eps_v,f(0) = 9.83645033e-4, G(p_min) = 5769.156 (2G = 11538.31), the in-domain
          (d eps_v = +1e-4) and out-of-domain (+1.1e-4) isotropic increments floor to the same eps_v,f with d eps^f_v = 6.9522805e-5 /
          7.9522805e-5, no refusal; the same increments are ordinary elastic steps under BA06 (M-F4); p_min = 0 refuses (HAR file)
  K1.14   HAR floored state under shear: eps_s = 2e-4 -> x = 0.09185198375, eps_v,f = 1.04102893e-3, q_f = 14.8362051, eps'_f =
          0.08679793, eps_s unchanged, tangent a^e Phi with delta:C = 0 and delta:C/max|C| = 1.45e-2 without the eps' term
  K1.14b  HAR dry-side FPf (the printed state, increment and every printed intermediate), the chains FPf,FPf and -P-,-Pf
  FD      (S.32f) / (S.54) against the central FD of the O2 stress update: BA06 cases (A) elastic + trial floor, (B) plastic + trial
          floor, (C) plastic + post floor, (D) m = 2 chain, alpha0 = 0 and 5 (the D12 != 0 stand-in), and the HAR cases above
  S.50    the floor closed forms (S.49), (S.50) for n in {0, 0.3, 1/2, 0.7} against p(eps_v,f, eps_s) = -p_min in the sheet's own
          forward map and eps'_f = -D12/D11; E_f = Psi difference = the quadrature of |p| d eps_v, in [0, W_f] (S.52)
  count   counters: n_f_tr / n_f_post / n_f_init, eps^f_v, W_f = p_min eps^f_v, at_floor, sub-increments sum, init projection,
          no refusal and the invariant p <= -p_min on every committed state of a random census, D >= 0 under the floor
  conv    O2 -> O1 first order on a floored path (S.55), alpha0 = 0 and 5: the sheet's printed error ladder

MUTANTS (named for the mutation gate, plan 5.3; the numbers are the sheet's).
  M-F1  floor not applied (pass-through or refusal): K1.12 / K1.13 expect a committed p = -p_min and no refusal.
  M-F2  applied but not counted: every counter assertion.            M-F3a unprojected tangent: delta:C != 0, 0.67 .. 2.1 off the FD.
  M-F3b eps'_f dropped: K1.14 delta:C (HAR) and the FD of the HAR cases.        M-F3c v-column tied to the floored trial: FD (B), 2.6e-3.
  M-F3d Phi operators omitted from the chain: FD (D) and FPf,FPf / -P-,-Pf, 0.98 .. 2.1.   M-F4 HAR floor with the BA06 inverse: K1.13
  (0.0530 vs 9.84e-4).  M-F5 trial floor skipped (post only): K1.13 out-of-domain trial.  M-F6 projection direction wrong (q kept):
  K1.14 eps_s unchanged and q_f.  M-F7 pi_i or v altered by the projection: K1.12 equalities.  M-F8 default reference wrong: the
  default-p_min test of test_g1_har.py.  M-F9 floor inside the local Newton: the printed K1.12 (C) and K1.14b values (Delta lambda, pi_i,
  q_c, the FPf pattern itself).
"""
import math

import numpy as np
import pytest

import g1_hf_common as C
from conftest import ORACLES

O1, O2 = ORACLES["O1"], ORACLES["O2"]
KER = O2.kernel
O2api = O2.api
I3 = np.eye(3)
ONES = np.ones(3)
SQ23, SQ32 = C.SQ23, C.SQ32

PI_FAR = -5000.0
K2_P0, K2_KH, K2_MU0 = -100.0, 0.01, 5400.0
PMIN_K2 = 0.5          # 5e-3 |p0| (sheet 1.3)


# ----------------------------------------------------------------------------------------------------------------------
# plumbing
# ----------------------------------------------------------------------------------------------------------------------
def make_state(oname, kw, sigma0, pi0, v0):
    P = C.make_hf(oname, **kw)
    return P, ORACLES[oname].initial_state(P, sigma0, v0, pi0)


def run1(oname, P, st, deps):
    """One increment (a 3 x 3 tensor); O1 integrates at rtol 1e-10."""
    deps = np.asarray(deps, float)
    if oname == "O1":
        return O1.run_path(P, st, deps[None], rtol=1e-10)[-1]
    return O2.step(P, st, deps)


def inv_e(st):
    return C.invariants_eps(st.eps_e)


def pressure(st):
    return float(np.trace(st.sigma)) / 3.0


def dev_norm(st):
    return float(np.linalg.norm(st.sigma - pressure(st) * I3))


def f_counters(oname, st):
    """(eps^f_v, W_f, at_floor) from either oracle's state."""
    if oname == "O2":
        return st.eps_f_v, st.W_f, bool(st.flags["at_floor"])
    return st.flags.get("eps_f_v", 0.0), st.flags.get("W_f", 0.0), bool(st.flags.get("at_floor", False))


def iso(p):
    return p * I3


def sym_identity4():
    return 0.5 * (np.einsum("ik,jl->ijkl", I3, I3) + np.einsum("il,jk->ijkl", I3, I3))


def iso_dev4():
    """I_sym - delta delta / 3 (the deviatoric projector)."""
    return sym_identity4() - np.einsum("ij,kl->ijkl", I3, I3) / 3.0


def tan_block(P, st):
    C4 = O2.tangent(P, st)
    return np.array([[C4[a, a, b, b] for b in range(3)] for a in range(3)])


def rel(A, B):
    return float(np.abs(A - B).max() / np.abs(B).max())


def fd_block(P, st0, deps, h, stepper):
    """Central FD of the diagonal stress over the three normal strain components, all 6 evaluation points on ONE branch
    pattern (fpattern + plastic); returns (block, pattern set)."""
    cols, pats = [], set()
    for b in range(3):
        E = np.zeros((3, 3))
        E[b, b] = h
        sp, sm = stepper(P, st0, deps + E), stepper(P, st0, deps - E)
        for s in (sp, sm):
            assert not s.flags["refused"], s.flags["reason"]
            pats.add((s.flags["fpattern"], bool(s.flags["plastic"])))
        cols.append((np.diag(sp.sigma) - np.diag(sm.sigma)) / (2.0 * h))
    return np.array(cols).T, pats


def single(P, st, d):
    return O2.step(P, st, d)


def mutants(stn):
    """The sheet's negative controls, built from the returned blocks (res.ae at the committed strain, res.chain, res.Phi_*):
    'none'    both floor operators dropped (the plain a^ep, or a^e when elastic),
    'no_post' Phi^post dropped, 'no_tr' Phi^tr dropped."""
    res = stn.cache["res"]
    if res.chain is None:                                           # elastic: C = a^e Phi^tr
        return dict(none=res.ae, no_post=res.ae @ (res.Phi_tr if res.Phi_tr is not None else I3), no_tr=res.ae)
    ch = res.chain
    vcol = np.outer(ch.u[:3] * ch.Pi_v * stn.v, ONES)
    Ptr = res.Phi_tr if res.Phi_tr is not None else I3
    Ppo = res.Phi_post if res.Phi_post is not None else I3
    return dict(none=res.ae @ (ch.b[:3, :3] - vcol), no_post=res.ae @ (ch.b[:3, :3] @ Ptr - vcol),
                no_tr=res.ae @ Ppo @ (ch.b[:3, :3] - vcol))


def fractions(fr):
    return lambda P, st, d: O2.step_fractions(P, st, d, fr)


def best_fd_error(P, st0, deps, stepper, hs):
    """{h: relative error of the returned tangent against the central FD}, the tangent of the state the stepper returns at deps."""
    stn = stepper(P, st0, deps)
    assert not stn.flags["refused"], stn.flags["reason"]
    Cb = tan_block(P, stn)
    out = {}
    for h in hs:
        fd, pats = fd_block(P, st0, deps, h, stepper)
        assert len(pats) == 1, f"FD points left the branch of the base increment: {pats}"
        out[h] = rel(Cb, fd)
    return stn, Cb, out


# ======================================================================================================================
# S.49 / S.50 closed forms: the sheet's own forward map is the check (not the oracle)
# ======================================================================================================================
@pytest.mark.sheet_check
@pytest.mark.parametrize("n", [0.0, 0.3, 0.5, 0.7])
@pytest.mark.parametrize("es", [0.0, 5e-5, 2e-4, 1e-3])
def test_har_floor_closed_form_solves_p_equals_minus_pmin_and_eps_prime_is_minus_D12_over_D11(n, es):
    """(S.50): eps_v,f(eps_s) is the solution of p(eps_v, eps_s) = -p_min of the forward map (S.5h), q_f = 3 g p_a eps_s x^n, and
    eps'_f = -D12/D11 at p = -p_min (the implicit derivative; D from (S.5h')).  TIMs k, g, p_a = 101, p_min = 5e-3 p_a, n in
    {0, 0.3, 0.5, 0.7}; coded in g1_hf_common from the sheet (n = 1/2 closed, general n the exact bracket).  This is a check of the
    sheet and of the common closed form against the sheet's other formulas; the oracle is compared to it below.
    KILLS (as a reference for the oracle tests): a wrong x-equation (a, b, exponent 2n), a wrong eps'_f formula."""
    pmin = 5e-3 * C.PA_T
    f = C.har_floor(es, pmin, n=n)
    p, q = C.har_pq(f["ev_f"], es, n=n)
    assert abs(p + pmin) <= 1e-11 * pmin, (p, pmin)
    assert abs(q - f["q_f"]) <= 1e-10 * max(q, 1e-12)
    D11, D12, D22, _ = C.har_D(p, q, n=n)
    assert abs(f["epsp"] - (-D12 / D11)) <= 1e-9 * max(abs(f["epsp"]), 1e-12), (f["epsp"], -D12 / D11)
    # finite-difference confirmation of the implicit derivative (independent of the formula): p(ev_f(es + h), es + h) = -pmin
    h = 1e-8
    fp_, fm_ = C.har_floor(es + h, pmin, n=n), C.har_floor(max(es - h, 0.0) if es > h else es, pmin, n=n)
    if es > 1e-6:
        assert abs((fp_["ev_f"] - fm_["ev_f"]) / (2 * h) - f["epsp"]) <= 1e-5 * max(abs(f["epsp"]), 1e-9)


@pytest.mark.sheet_check
def test_sheet_printed_floor_numbers_k1_13_and_k1_14_reproduce_from_the_closed_forms():
    """The printed digits of 13.13 / 13.14 against the closed forms of (S.50) at the TIMs set (n = 1/2): eps_v,f(0) =
    9.83645033e-4 (domain edge 1.05849170e-3), G(p_min) = g p_a (p_min/p_a)^(1/2) = 5769.156 and 2G = 11 538.31; at eps_s = 2e-4:
    x = 0.09185198375, eps_v,f = 1.04102893e-3, q_f = 14.8362051 kPa, eta_f = 29.38, eps'_f = 0.08679793.  Pure sheet check."""
    pmin = C.PMIN_T
    f0 = C.har_floor(0.0, pmin)
    assert abs(f0["ev_f"] - 9.83645033e-4) <= 6e-12 and abs(C.EDGE_T - 1.05849170e-3) <= 1e-11
    G = C.G_T * C.PA_T * math.sqrt(pmin / C.PA_T)
    assert abs(G - 5769.156) <= 6e-4 and abs(2 * G - 11538.31) <= 6e-3
    f = C.har_floor(2e-4, pmin)
    assert abs(f["x"] - 0.09185198375) <= 6e-12 and abs(f["ev_f"] - 1.04102893e-3) <= 6e-12
    assert abs(f["q_f"] - 14.8362051) <= 6e-8 and abs(f["q_f"] / pmin - 29.38) <= 5e-3 and abs(f["epsp"] - 0.08679793) <= 6e-9


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_floor_ev_of_the_oracle_matches_the_sheet_closed_forms_ba06_and_har(a0):
    """The O2 kernel operator floor_ev (the only direct kernel call of this file) against (S.49) and (S.50) coded in the common
    module: BA06 K2 (alpha0 = 0 and 5) at 4 values of eps_s, HAR n in {0, 0.3, 1/2, 0.7} at 4 values.  1e-12 relative.
    KILLS: a wrong (S.49) alpha0 factor or log argument; the BA06 inverse used under HAR (M-F4); a wrong HAR exponent or bracket."""
    for es in (0.0, 3e-3, 1e-2, 3e-2):
        P = C.make_hf("O2", **C.ba06_kw(p_min=PMIN_K2, alpha0=a0))
        ev_f, epsp = KER.floor_ev(P, es)
        ref = C.ba06_floor(es, PMIN_K2, a0=a0)
        assert abs(ev_f - ref["ev_f"]) <= 1e-12 * abs(ref["ev_f"]) and abs(epsp - ref["epsp"]) <= 1e-12 * max(abs(ref["epsp"]), 1.0)
        p, _ = C.ba06_pq(ev_f, es, a0=a0)
        assert abs(p + PMIN_K2) <= 1e-12 * PMIN_K2
    if a0 == 0.0:
        assert abs(KER.floor_ev(C.make_hf("O2", **C.ba06_kw(p_min=PMIN_K2)), 0.0)[0] - 0.0529831737) <= 6e-11      # printed (13.12)
        for n in (0.0, 0.3, 0.5, 0.7):
            for es in (0.0, 5e-5, 2e-4, 1e-3):
                P = C.make_hf("O2", **C.har_kw(p_min="default", n_e=n))
                ev_f, epsp = KER.floor_ev(P, es)
                ref = C.har_floor(es, C.PMIN_T, n=n)
                assert abs(ev_f - ref["ev_f"]) <= 1e-12 * abs(ref["ev_f"]), (n, es, ev_f, ref["ev_f"])
                assert abs(epsp - ref["epsp"]) <= 1e-9 * max(abs(ref["epsp"]), 1e-9), (n, es)


# ======================================================================================================================
# K1.12  BA06 floor (sheet 13.12)
# ======================================================================================================================
def k112_setup(oname, pmin=PMIN_K2, shear=True, theta=0.45):
    """K2 BA06 (alpha0 = 0), p_min given; start at p = -1 kPa with (optionally) a deviatoric stress along a direction of Lode
    angle `theta`; the surface is far (pi_i0 = -5000: the whole path is elastic).  Returns (P, st0, deps) with the pure
    expansion of d eps_v = kappa ln 4 that takes p from -1 to -0.25 = -p_min/2 (13.12)."""
    kw = C.ba06_kw(p_min=pmin)
    es0 = 2.0e-5 if shear else 0.0       # q = 0.324 kPa: F < 0 at p = -0.5 on the huge surface (q = 1.62 would leave it at the floor)
    nh = C.dir_for_theta(theta)
    q0 = 3.0 * K2_MU0 * es0
    sig0 = C.sigma_pq(-1.0, q0, nh)
    P, st0 = make_state(oname, kw, np.diag(sig0), PI_FAR, 1.70)
    dev = K2_KH * math.log(4.0)                 # p = p0 exp(-eps_v/kappa): -1 -> -0.25
    return P, st0, np.eye(3) * dev / 3.0


@pytest.mark.parametrize("shear", [False, True])
def test_k1_12_ba06_trial_floor_closed_forms_counters_and_unchanged_fields(oracle_name, shear):
    """13.12, K2 set, alpha0 = 0, p_min = 5e-3 |p0| = 0.5 kPa.  A trial at p^tr = -p_min/2 = -0.25 (eps_v,tr = 0.0599146455) projects
    to eps_v,f = -kappa ln(p_min/|p0|) = 0.0529831737 (printed), d eps^f_v = kappa ln 2 = 6.931471806e-3, W_f = p_min d eps^f_v =
    3.465735903e-3, committed p = -0.5 exactly.  The pure expansion from p = -1 (with and without a deviatoric stress of 0.324 kPa:
    alpha0 = 0 => q = 3 mu0 eps_s depends on eps_s alone) leaves q, n_hat, eps_s, pi_i unchanged (M-F6, M-F7) and v = v_n exp(tr
    d eps) built from the RAW trace (the floor does not touch v).  O2: pattern 'FE-', counters (n_f_tr, n_f_post, n_f_init) = (1, 0,
    0), at_floor; the E_f response (api.floor_energy) = kappa (p_min - |p^tr|) = 2.5e-3 and 0 <= E_f <= W_f (S.52).  O1: the rate
    form (S.55) with lam_f = tr eps. gives the same final state (BA06 alpha0 = 0; sheet 9.7 O1 paragraph) to the integration tolerance.
    KILLS: M-F1 (no floor: p = -0.25), M-F2 (not counted), M-F6 (q or eps_s changed), M-F7 (pi_i or v altered), a floor evaluated at
    the wrong p (the committed p is not -p_min), E_f mis-signed."""
    P, st0, deps = k112_setup(oracle_name, shear=shear)
    q0, ev0 = C.invariants_pq(st0.sigma)[1], inv_e(st0)[0]
    nhat0 = (st0.sigma - pressure(st0) * I3) / max(dev_norm(st0), 1e-300)
    es0 = inv_e(st0)[1]
    st = run1(oracle_name, P, st0, deps)
    assert C.ok_state(oracle_name, st), C.why_state(oracle_name, st)
    tol = 1e-13 if oracle_name == "O2" else 1e-7
    assert abs(pressure(st) + PMIN_K2) <= tol * PMIN_K2 * (1 if oracle_name == "O2" else 10), pressure(st)
    ev, es, _ = inv_e(st)
    assert abs(ev - 0.0529831737) <= 6e-11, ev                                    # printed eps_v,f (alpha0 = 0: independent of eps_s)
    ev_f_closed = C.ba06_floor(es, PMIN_K2)["ev_f"]
    assert abs(ev - ev_f_closed) <= (1e-12 if oracle_name == "O2" else 1e-7) * abs(ev_f_closed)
    assert abs(es - es0) <= (1e-14 if oracle_name == "O2" else 1e-9) + 1e-12 * es0, "eps_s changed by the projection (M-F6)"
    assert abs(C.invariants_pq(st.sigma)[1] - q0) <= (1e-12 if oracle_name == "O2" else 1e-7) * max(q0, 1e-3)
    if shear:
        nh = (st.sigma - pressure(st) * I3) / dev_norm(st)
        assert np.linalg.norm(nh - nhat0) <= (1e-13 if oracle_name == "O2" else 1e-7), "n_hat changed (M-F6)"
    assert st.pi_i == PI_FAR or abs(st.pi_i - PI_FAR) <= 1e-9 * abs(PI_FAR), "pi_i altered by the projection (M-F7)"
    assert abs(st.v - st0.v * math.exp(float(np.trace(deps)))) <= (1e-13 if oracle_name == "O2" else 1e-8) * st.v, "v not v_n exp(raw tr) (M-F7)"
    dfv = (math.log(2.0)) * K2_KH
    efv, Wf, atf = f_counters(oracle_name, st)
    # the floor volumetric strain is trial minus floored: tr d eps - (ev0 -> ev_f); with the start at p = -1: eps_v0 = kappa ln 100
    assert abs(efv - dfv) <= (1e-13 if oracle_name == "O2" else 1e-7) * dfv, (efv, dfv)
    assert abs(Wf - 3.465735903e-3) <= 6e-13 + (0 if oracle_name == "O2" else 1e-9), Wf        # printed
    assert abs(Wf - PMIN_K2 * efv) <= 1e-15 + 1e-12 * Wf
    assert atf
    if oracle_name == "O2":
        assert (st.n_f_tr, st.n_f_post, st.n_f_init) == (1, 0, 0), "floor events not counted exactly (M-F2)"
        assert st.flags["fpattern"] == "FE-" and (st.flags["floor_tr"], st.flags["floor_post"]) == (1, 0)
        Ef = O2api.floor_energy(P, st)
        assert abs(Ef - K2_KH * (PMIN_K2 - 0.25)) <= 1e-13, Ef                             # kappa (p_min - |p^tr|) = 2.5e-3
        assert 0.0 <= Ef <= Wf
        # the raw trial strain is eps_v,tr = eps_v(-0.25) = kappa ln(100/0.25) = 0.0599146455 (printed)
        assert abs(ev0 + float(np.trace(deps)) - 0.0599146455) <= 6e-11


def test_k1_12_floored_tangent_is_2mu0_deviator_with_zero_bulk_response_noncoaxial(oracle_name):
    """13.12 / 9.7: the tangent of an elastic trial-floored step under BA06 alpha0 = 0 is a^e(eps_f) Phi^tr with a^e = K delta delta
    + 2 mu0 (I - delta delta/3) and Phi = I - delta/3 (eps' = 0): C = 2 mu0 (I_sym - delta delta/3) = 10 800 (I_sym - delta delta/3),
    delta:C = 0 -- the FULL 4th-order tensor (all 81 components, shear columns included), for a NON-COAXIAL start and a non-coaxial
    increment (so the spin terms of (S.33) are exercised; Phi has unit spin).  The O1 continuum tangent at the floored state (floor
    branch of S.55) is the same operator.  No bulk stiffness is faked (owner decision (c)).
    KILLS: M-F3a (unprojected a^e: delta:C = K = 100/0.01 x ... != 0), a bulk regularisation added to the tangent, the spin of Pi_f
    wrong (a shear block off by the (sigma_a - sigma_b)/(eps_a - eps_b) ratio), a floor tangent used on an unfloored step."""
    kw = C.ba06_kw(p_min=PMIN_K2)
    sig0 = np.array([[-1.0, 0.03, 0.01], [0.03, -1.12, 0.02], [0.01, 0.02, -0.95]])      # small deviator: F < 0 at the floor (|p| = 0.5)
    P, st0 = make_state(oracle_name, kw, sig0, PI_FAR, 1.70)
    deps = np.array([[0.0140 + 1.5e-5, 4e-6, -2e-6], [4e-6, 0.0140 - 1.0e-5, 3e-6], [-2e-6, 3e-6, 0.0140 - 0.5e-5]])    # expansion + small shear
    st = run1(oracle_name, P, st0, deps)
    assert C.ok_state(oracle_name, st)
    if oracle_name == "O2":
        assert st.flags["fpattern"] == "FE-", st.flags["fpattern"]
    else:
        assert st.flags.get("floor", False), "O1 did not end on the floor"
    Ct = O2.tangent(P, st) if oracle_name == "O2" else O1.tangent(P, st)
    want = 2.0 * K2_MU0 * iso_dev4()
    assert np.abs(Ct - want).max() <= 1e-9 * 2.0 * K2_MU0, np.abs(Ct - want).max()
    assert np.abs(np.einsum("iikl->kl", Ct)).max() <= 1e-9 * 2.0 * K2_MU0, "delta:C != 0 at a floored state"


def test_floor_is_idempotent_and_a_second_zero_increment_is_not_a_floor_event(oracle_name):
    """9.7 determinism: after a projection p = -p_min to round-off and the activation test p > -p_min (1 - 1e-12) is false, so a
    second (zero) increment does nothing: the pattern is '-E-', the counters do not move, p stays -p_min, the state is unchanged.
    KILLS: re-projection chatter (n_f_tr increments on a state already on the floor), an activation test with the wrong tolerance
    sign (p = -p_min counts as above the floor), a floor that drifts."""
    if oracle_name == "O1":
        pytest.skip("O2 counters; the O1 rate form has no per-increment event count (sheet 9.7 'Counted' is per Gauss point of the kernel)")
    P, st0, deps = k112_setup("O2")
    st1 = O2.step(P, st0, deps)
    st2 = O2.step(P, st1, np.zeros((3, 3)))
    assert st1.flags["fpattern"] == "FE-" and st2.flags["fpattern"] == "-E-"
    assert (st2.n_f_tr, st2.n_f_post, st2.eps_f_v) == (st1.n_f_tr, st1.n_f_post, st1.eps_f_v)
    assert np.abs(st2.sigma - st1.sigma).max() <= 1e-14 and np.abs(st2.eps_e - st1.eps_e).max() <= 1e-15
    # a small compressive step leaves the floor with no event: the full volumetric stiffness is back (sheet 9.7: one-sided)
    st3 = O2.step(P, st2, -1e-4 * I3)
    assert st3.flags["fpattern"] == "-E-" and not st3.flags["at_floor"] and pressure(st3) < -PMIN_K2 * (1.0 + 1e-6)
    K3 = float(np.einsum("iijj->", O2.tangent(P, st3))) / 9.0
    assert abs(K3 - (-pressure(st3) / K2_KH)) <= 1e-8 * K3, "the bulk stiffness K = -p/kappa did not return after leaving the floor"


def test_floor_off_restores_the_unfloored_response_and_counts_nothing(oracle_name):
    """-pmin 0 = off (9.7, 1.3): the very increment of K1.12 is an ordinary elastic step to p = p0 exp(-eps_v/kappa) = -0.25, no
    counter moves, at_floor false.  KILLS: M-F1 inverted (floor applied at p_min = 0 / a hard-coded floor), a counter incremented
    with the floor off, the pre-round-3 behaviour not restored by `p_min = 0`."""
    P, st0, deps = k112_setup(oracle_name, pmin=0.0)
    st = run1(oracle_name, P, st0, deps)
    assert C.ok_state(oracle_name, st)
    assert abs(pressure(st) + 0.25) <= 1e-9, pressure(st)
    efv, Wf, atf = f_counters(oracle_name, st)
    assert efv == 0.0 and Wf == 0.0 and not atf
    if oracle_name == "O2":
        assert (st.n_f_tr, st.n_f_post, st.n_f_init) == (0, 0, 0) and st.flags["fpattern"] == "-E-"


# ======================================================================================================================
# K1.13  HAR floor and the out-of-domain trial (sheet 13.13)
# ======================================================================================================================
EV_P1 = 9.53167838e-4          # the isotropic state at p = -1 kPa (printed, 13.13)


def har_iso_start(oname, pmin="default"):
    kw = C.har_kw(p_min=pmin)
    P, st = make_state(oname, kw, iso(-1.0), PI_FAR, 1.70)
    return P, st


@pytest.mark.parametrize("dev,dfv_printed,wf_printed,out_of_domain", [(1.0e-4, 6.9522805e-5, 3.5109017e-5, False),
                                                                      (1.1e-4, 7.9522805e-5, None, True)])
def test_k1_13_har_isotropic_floor_in_and_out_of_the_domain(oracle_name, dev, dfv_printed, wf_printed, out_of_domain):
    """13.13, TIMs HAR set, p_min = 5e-3 p_a = 0.505 kPa: an isotropic state at p = -1 kPa (eps_v = 9.53167838e-4, the inverse map
    (S.5h'')) with d eps_v = +1e-4 is still in the domain (p^tr = -2.56e-3 kPa) and floors to eps_v,f = 9.83645033e-4 with d eps^f_v
    = 6.9522805e-5, W_f = 3.5109017e-5; with d eps_v = +1.1e-4 the trial (eps_v = 1.0632e-3 > 1.0585e-3 = the domain edge) is OUT of
    the domain and floors to the SAME eps_v,f with d eps^f_v = 7.9522805e-5: NO REFUSAL (M-F1, M-F5).  Committed p = -0.505, q = 0
    (isotropic), eps_s = 0, pi_i unchanged.  O2: pattern 'FE-', n_f_tr = 1 (the out-of-domain trial counts as a trial event).  O1: the
    rate form floors after the floor-hit point and ends on the same state (HAR05 isotropic axis: D12 = 0).  The shear stiffness at
    the floor: G(p_min) = 5769.156 kPa, the tangent is 2G (I_sym - delta delta/3), delta:C = 0 (volumetric response 0).
    KILLS: M-F1 (refusal or pass-through of the out-of-domain trial), M-F4 (the BA06 inverse: eps_v,f = 0.0530 vs 9.84e-4), M-F5
    (trial floor skipped: the out-of-domain trial refuses instead of flooring), a floor computed at the wrong p_ref (p_min = 5e-3 |p0|
    with p0 None), the shear stiffness regularised (2G != 11538.31)."""
    P, st0 = har_iso_start(oracle_name)
    assert abs(inv_e(st0)[0] - EV_P1) <= 1e-10 and abs(pressure(st0) + 1.0) <= 1e-12
    assert (C.EDGE_T < inv_e(st0)[0] + dev) == out_of_domain
    st = run1(oracle_name, P, st0, np.eye(3) * dev / 3.0)
    assert C.ok_state(oracle_name, st), ("a floor event was refused", C.why_state(oracle_name, st))
    ev, es, _ = inv_e(st)
    tol = 1e-12 if oracle_name == "O2" else 1e-7
    assert abs(ev - 9.83645033e-4) <= 6e-12 + tol * 1e-3, ev                                  # printed eps_v,f(0)
    assert abs(pressure(st) + C.PMIN_T) <= tol * C.PMIN_T * (1 if oracle_name == "O2" else 10), pressure(st)
    assert es <= 1e-14 and dev_norm(st) <= 1e-12 and abs(st.pi_i - PI_FAR) <= 1e-9 * abs(PI_FAR)
    efv, Wf, atf = f_counters(oracle_name, st)
    assert abs(efv - dfv_printed) <= 6e-13 + (0 if oracle_name == "O2" else 1e-10), (efv, dfv_printed)
    if wf_printed is not None:
        assert abs(Wf - wf_printed) <= 6e-13 + (0 if oracle_name == "O2" else 1e-10), (Wf, wf_printed)
    assert abs(Wf - C.PMIN_T * efv) <= 1e-14 + 1e-12 * Wf and atf
    if oracle_name == "O2":
        assert (st.n_f_tr, st.n_f_post, st.n_f_init) == (1, 0, 0) and st.flags["fpattern"] == "FE-"
        Ef = O2api.floor_energy(P, st)
        assert 0.0 <= Ef <= Wf + 1e-18
        if not out_of_domain:
            # in the domain: E_f = Psi(eps_f) - Psi(pre) (S.52); the pre-floor state is eps_v = 9.53167838e-4 + 1e-4
            want = C.har_psi(C.har_floor(0.0, C.PMIN_T)["ev_f"], 0.0) - C.har_psi(inv_e(st0)[0] + dev, 0.0)
            assert abs(Ef - want) <= 1e-9 * abs(want) + 1e-18, (Ef, want)
        else:
            assert abs(Ef - Wf) <= 1e-18, "out-of-domain trial: only the bound W_f is stated (S.52)"
    Ct = O2.tangent(P, st) if oracle_name == "O2" else O1.tangent(P, st)
    G = C.G_T * C.PA_T * math.sqrt(C.PMIN_T / C.PA_T)
    assert abs(G - 5769.156) <= 6e-4
    assert np.abs(Ct - 2.0 * G * iso_dev4()).max() <= 1e-8 * 2.0 * G, "floored HAR tangent is not 2G(p_min) (I - delta delta/3)"
    assert np.abs(np.einsum("iikl->kl", Ct)).max() <= 1e-8 * 2.0 * G


def test_k1_13_control_the_same_increments_are_ordinary_elastic_steps_under_ba06(oracle_name):
    """13.13 (M-F4): under BA06 (K2, p_min = 0.5) the same two isotropic increments from p = -1 (+1e-4, +1.1e-4) are ordinary
    elastic steps: p = p0 exp(-eps_v/kappa) = -0.99005 / -0.98911..., no floor event.  Makes the HAR assertions above discriminate.
    KILLS: a floor computed with the wrong energy (events where there are none)."""
    for dev in (1.0e-4, 1.1e-4):
        P, st0 = make_state(oracle_name, C.ba06_kw(p_min=PMIN_K2), iso(-1.0), PI_FAR, 1.70)
        st = run1(oracle_name, P, st0, np.eye(3) * dev / 3.0)
        assert abs(pressure(st) - (-1.0 * math.exp(-dev / K2_KH))) <= 1e-9
        if oracle_name == "O2":
            assert (st.n_f_tr, st.n_f_post) == (0, 0) and st.eps_f_v == 0.0


# ======================================================================================================================
# K1.14  HAR floored state under shear (sheet 13.14)
# ======================================================================================================================
EPS_V_TR_14 = 1.050e-3      # a trial expansion into the floor with eps_s = 2e-4 (in the domain: p^tr = -0.245 kPa)


def k114_build(es_target=2e-4):
    """O2 on the TIMs HAR set: from the isotropic state at p = -100 a single elastic increment to (eps_v, eps_s) = (1.050e-3,
    2e-4), n_hat at theta = 0.5 with ASCENDING principal components (so the principal block equals C[a,a,b,b])."""
    kw = C.har_kw(p_min="default", M=50.0)           # M = 50: the huge surface stays above eta_f = 29 (every NorSand surface has eta <= M/N)
    P, st0 = make_state("O2", kw, iso(-100.0), PI_FAR, 1.70)
    ev0 = inv_e(st0)[0]
    nh = np.sort(C.dir_for_theta(0.5))
    deps = np.diag((EPS_V_TR_14 - ev0) / 3.0 * ONES + SQ32 * es_target * nh)
    return P, st0, deps, nh


def ae_har(p, q, nh):
    """(S.3) a^e in the principal basis for HAR (D12 != 0, q/eps_s != D22), from the stress-form (S.5h') in the common module."""
    D11, D12, D22, ratio = C.har_D(p, q)
    return (D11 * np.outer(ONES, ONES) + SQ23 * D12 * (np.outer(ONES, nh) + np.outer(nh, ONES)) + (2.0 / 3.0) * D22 * np.outer(nh, nh)
            + (2.0 * ratio / 3.0) * (I3 - np.outer(ONES, ONES) / 3.0 - np.outer(nh, nh)))


def test_k1_14_har_floor_under_shear_printed_values_and_tangent():
    """13.14, TIMs HAR (M = 50 keeps the whole path elastic: eta_f = 29 lies above any physical surface, only the elastic floor operator
    is under test), p_min = 0.505: a trial with eps_s = 2e-4 and eps_v = 1.050e-3 (in the domain, p^tr = -0.245 > -p_min) floors to
    x = varpi_f/p_a = 0.09185198375, eps_v,f = 1.04102893e-3, q_f = 14.8362051 kPa (eta_f = 29.38); eps_s^e UNCHANGED by Pi_f (M-F6:
    the projection keeps the deviatoric elastic strain, q is NOT kept); committed p = -0.505.  Tangent: C_f = a^e(eps_f) Phi with Phi =
    I - 1/3 + (1/3) eps'_f sqrt(2/3) 1 n_hat^T, eps'_f = 0.08679793 (printed): a^e from (S.3) with the HAR D of (S.5h') coded in the
    common module; delta:C_f = 0 with the eps' term; without it (Phi = I - 1/3) delta:C = sqrt(6) D12 n_hat (M-F3b; the sheet prints delta:C/max|C| = 1.45e-2 without a direction).
    KILLS: M-F3b (eps'_f dropped: delta:C = 1.45e-2 max|C|), M-F6 (q kept: q_f != 14.836 and eps_s changes), M-F4 (BA06 floor value),
    a^e built without the t4 / t2 terms of (S.3), the HAR floor evaluated at q = 0."""
    P, st0, deps, nh = k114_build()
    st = O2.step(P, st0, deps)
    assert st.flags["fpattern"] == "FE-", st.flags["fpattern"]
    ev, es, _ = inv_e(st)
    assert abs(es - 2e-4) <= 1e-13, "eps_s changed by the projection (M-F6)"
    assert abs(ev - 1.04102893e-3) <= 6e-12, ev
    p, q, _ = C.invariants_pq(st.sigma)
    assert abs(p + C.PMIN_T) <= 1e-12 and abs(q - 14.8362051) <= 6e-8 * 1 and abs(q / C.PMIN_T - 29.38) <= 5e-3
    vp = math.sqrt(p * p + C.K_T * (1 - C.N_T) * q * q / (3.0 * C.G_T))
    assert abs(vp / C.PA_T - 0.09185198375) <= 6e-12
    # tangent
    epsp = 0.08679793
    Phi = I3 - 1.0 / 3.0 + (1.0 / 3.0) * epsp * SQ23 * np.outer(ONES, nh)
    ae = ae_har(p, q, nh)
    Cexp = ae @ Phi
    Cb = tan_block(P, st)
    assert np.abs(Cb - Cexp).max() <= 5e-8 * np.abs(Cexp).max(), np.abs(Cb - Cexp).max() / np.abs(Cexp).max()      # eps'_f printed to 7 digits
    assert np.abs(Cb.sum(axis=0)).max() <= 1e-12 * np.abs(Cb).max(), "delta:C_f != 0 at a floored HAR state"
    # M-F3b control: without the eps' term C = a^e (I - 1/3), and by (S.3) (n_hat.1 = 0) delta:C = sqrt(6) D12 n_hat EXACTLY; the ratio
    # delta:C/max|C| must exceed the sheet's printed 1.45e-2 (the printed digits themselves are not reproduced: the sheet gives no
    # direction n_hat and every direction from the compression to the extension meridian gives 0.075-0.078 here -- reported)
    Cmut = ae @ (I3 - 1.0 / 3.0)
    D12 = C.har_D(p, q)[1]
    dm = Cmut.sum(axis=0)
    assert np.abs(dm - math.sqrt(6.0) * D12 * nh).max() <= 1e-9 * np.abs(dm).max()
    assert np.abs(dm).max() / np.abs(Cmut).max() >= 1.45e-2, "M-F3b would go unseen"
    assert rel(Cb, Cmut) >= 1e-3, "the eps'_f term is invisible in the tangent: the gate has no power"


# ======================================================================================================================
# FD of (S.32f) / (S.54): BA06 cases (A) - (D), alpha0 = 0 and 5
# ======================================================================================================================
PMIN_FD = 50.0


def ba06_fd_params(a0):
    return C.make_hf("O2", **C.ba06_kw(p_min=PMIN_FD, alpha0=a0))


def case_A(a0):
    P = ba06_fd_params(a0)
    st0 = O2.initial_state(P, np.diag([-57.0, -61.0, -64.0]), 1.70, -300.0)
    return P, st0, np.diag([2.0e-3, 1.5e-3, 2.5e-3])


def _off_corner_state(P, p, eta, direction, v):
    sig, th, q = C.sigma_on_surface(1.2, 0.4, 0.7, p, eta, direction)
    pi = C.pi_S53(1.2, 0.4, p, eta)
    return O2.initial_state(P, np.diag(sig), v, pi), th


def case_B(a0):
    """dry side (eta = 1.4 > M), expansion + shear along n_hat until the first 'FP-' pattern of the sheet's search list."""
    P = ba06_fd_params(a0)
    st0, th = _off_corner_state(P, -55.0, 1.4, C.OFFC_DIR, 1.70)
    for sdev in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3):
        deps = np.diag(1.0e-3 * ONES + sdev * C.OFFC_DIR)
        s = O2.step(P, st0, deps)
        if not s.flags["refused"] and s.flags["fpattern"] == "FP-":
            return P, st0, deps
    raise AssertionError("case (B) FP- not reached on the sheet's search list")


def case_C(a0):
    """wet side (eta = 0.5 M), the post floor ('-Pf') of 13.12 (C)."""
    P = ba06_fd_params(a0)
    st0, th = _off_corner_state(P, -51.0, 0.5 * 1.2, C.WET_DIR, 1.70)
    for amp in (2e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1.2e-2):
        deps = np.diag(amp * np.array([0.55, 0.40, -0.95]) - 1e-5)
        s = O2.step(P, st0, deps)
        if not s.flags["refused"] and s.flags["fpattern"] == "-Pf":
            return P, st0, deps
    raise AssertionError("case (C) -Pf not reached on the sheet's search list")


# sheet 9.7 FD record, 3x the printed best-h error (h = 1e-6 / 1e-7), floored at 1e-9
GATE_A = 1e-9
GATE_B = 6e-8
GATE_C = 1.5e-7
GATE_D = 1e-8


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_floor_fd_A_elastic_with_trial_floor(a0):
    """(S.32f) elastic: C_f = a^e(eps_f) Phi^tr against the central FD of the O2 step, K2 BA06 with p_min scaled to 50 kPa (the
    algebra is scale-free), alpha0 = 0 and 5 (D12 != 0: the eps' term live), off the WW corners.  Sheet 13.12 / 9.7: p^tr = -33.3 ->
    -50.00 (elastic, 'FE-'), d eps^f_v = 4.066e-3; printed FD record 8.1e-13 / 1.6e-11 (h = 1e-6 / 1e-7); gate 1e-9; delta:C = 0.
    Negative control: the unprojected a^e (the return map's own atilde) is 0.67-0.69 off the FD.
    KILLS: M-F3a (plain a^e), M-F3b at alpha0 = 5 (eps' dropped, delta:C = 1.1e-2), the wrong spin of Phi."""
    P, st0, deps = case_A(a0)
    stn, Cb, errs = best_fd_error(P, st0, deps, single, (1e-6, 1e-7))
    assert stn.flags["fpattern"] == "FE-"
    ev_tr = C.invariants_eps(st0.eps_e + deps)[0]
    es_tr = C.invariants_eps(st0.eps_e + deps)[1]
    p_tr, _ = C.ba06_pq(ev_tr, es_tr, a0=a0)
    assert abs(p_tr - (-33.3)) <= 0.06 and abs(pressure(stn) + PMIN_FD) <= 1e-12 * PMIN_FD
    if a0 == 0.0:
        assert abs(stn.eps_f_v - 4.066e-3) <= 6e-7, stn.eps_f_v                                     # printed
    assert min(errs.values()) <= GATE_A, errs
    assert np.abs(Cb.sum(axis=0)).max() <= 1e-13 * np.abs(Cb).max(), "delta:C != 0"
    assert rel(mutants(stn)["none"], Cb) >= 0.3, "the unprojected a^e is not far from C_f: the gate has no power"


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_floor_fd_B_plastic_with_trial_floor_v_column_on_the_raw_trial(a0):
    """(S.32f) plastic + trial floor: C_f = a^e(eps_c) [b Phi^tr - u Pi_v v_{n+1} 1^T], the v-column on the RAW trial strain (sheet
    9.7; v = v_n exp(tr d eps) is built from the total strain and the floor does not touch it).  Dry side (eta = 1.4 > M), pattern
    'FP-', K2 BA06 alpha0 = 0 / 5, p_min = 50.  Printed (13.12 B, alpha0 = 0): p_c = -50.70 (the return from a floored trial compresses
    p; no post floor).  FD record 5.2e-9 / 1.7e-8 (alpha0 = 0), 5.2e-9 / 5.1e-9 (alpha0 = 5); gate 6e-8.  Negative controls: plain
    atilde 0.93 off, the v-column tied to the floored trial (a~^ep Phi^tr) 2.6e-3 / 2.5e-3 off (so the gate separates the mutant by
    4e4 x).
    KILLS: M-F3a, M-F3c (v-column tied to the floored trial), M-F3d at m = 1 (Phi^tr omitted), M-F9 (floor inside the Newton: p_c is
    not the unconstrained value), Phi^post applied to an inactive post floor."""
    P, st0, deps = case_B(a0)
    stn, Cb, errs = best_fd_error(P, st0, deps, single, (1e-6, 1e-7))
    assert stn.flags["fpattern"] == "FP-" and stn.flags["plastic"]
    if a0 == 0.0:
        assert abs(pressure(stn) + 50.70) <= 6e-3, pressure(stn)                                    # printed -50.70 (-50.698)
    assert min(errs.values()) <= GATE_B, errs
    assert stn.n_f_post == 0, "case (B) must not post-floor (BA06 dry side)"
    mu = mutants(stn)
    assert rel(mu["none"], Cb) >= 0.3 and rel(mu["no_tr"], Cb) >= 0.3


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_floor_fd_C_plastic_with_post_floor_wet_side(a0):
    """(S.32f) plastic + post floor, wet side (eta = 0.5 M, F_p < 0): the return RELAXES p (p_c = -48.76 -> -50.00 by the post floor,
    pattern '-Pf', eta_c = 0.953; printed for alpha0 = 0), delta:C_f = 0 (the last operator is an active Pi_f).  FD record 4.6e-8 /
    5.0e-10, 4.2e-8 / 4.4e-10 (alpha0 = 5); gate 1.5e-7.  Negative control: dropping Phi^post (plain atilde) is 2.0-2.1 off.  F(sigma_f)
    is NOT asserted <= F_tol at the floored state (9.7 item 2-3): the wet-side post floor leaves F(sigma_f) > 0 by O(p_min).
    KILLS: M-F3a, M-F3d at m = 1 (Phi^post omitted), M-F9, M-F2 (post event not counted: n_f_post = 1)."""
    P, st0, deps = case_C(a0)
    stn, Cb, errs = best_fd_error(P, st0, deps, single, (1e-6, 1e-7))
    assert stn.flags["fpattern"] == "-Pf" and stn.flags["plastic"]
    assert (stn.n_f_tr, stn.n_f_post) == (0, 1), "post floor event not counted (M-F2)"
    pre = stn.cache["floor_events"][-1][0]
    ev, es, _ = C.invariants_eps(pre)
    p_c, _ = C.ba06_pq(ev, es, a0=a0)
    if a0 == 0.0:
        assert abs(p_c + 48.76) <= 6e-3 and abs(stn.eta - 0.953) <= 6e-4, (p_c, stn.eta)         # printed
    assert p_c > -PMIN_FD, "the unconstrained return must end ABOVE the floor (wet side)"
    assert min(errs.values()) <= GATE_C, errs
    assert np.abs(Cb.sum(axis=0)).max() <= 1e-12 * np.abs(Cb).max(), "delta:C != 0 after a post floor"
    assert rel(mutants(stn)["no_post"], Cb) >= 0.5


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_floor_fd_D_chain_m2_trial_floored_plastic_substeps(a0):
    """(S.54) / (S.46f): the m = 2 chain with the trial floor and a plastic return in each half (pattern 'FP-,FP-', fractions held
    fixed) against the FD of the whole increment; K2 BA06 alpha0 = 0 / 5, increment of case (B).  FD record 2.9e-9 / 1.8e-9 (alpha0
    = 0), 2.9e-9 / 5.5e-10 (alpha0 = 5); gate 1e-8.  Counters sum over the sub-increments: n_f_tr = 2.  Negative control: the
    last-sub-increment tangent (tangent_last_substep, which carries neither the chain nor the Phi operators of the earlier
    sub-increment) is 0.5-0.99 off.
    KILLS: M-F3d (Phi operators omitted from the chain: 0.98-0.99 off), a chain that drops the floored state of the first
    sub-increment, counters not summed over sub-increments."""
    P, st0, deps = case_B(a0)
    fr = (0.5, 0.5)
    stn, Cb, errs = best_fd_error(P, st0, deps, fractions(fr), (1e-6, 1e-7))
    assert stn.flags["fpattern"] == "FP-,FP-", stn.flags["fpattern"]
    assert (stn.n_f_tr, stn.n_f_post) == (2, 0) and (stn.flags["floor_tr"], stn.flags["floor_post"]) == (2, 0)
    assert min(errs.values()) <= GATE_D, errs
    last = np.array([[O2api.tangent_last_substep(P, stn)[a, a, b, b] for b in range(3)] for a in range(3)])
    assert rel(last, Cb) >= 1e3 * GATE_D, "the last-sub-increment tangent coincides with the chain: the gate has no power"


# ======================================================================================================================
# K1.14b  HAR dry-side FPf and the chained patterns
# ======================================================================================================================
def fpf_build():
    """The printed state of 13.14b: TIMs HAR + M 1.3309, rho = rho_bar 0.71 (WW), fork CSL at the TIMs values, K2 plastic constants,
    p_min = 0.505, state on the surface at p = -0.6 kPa, eta = 1.2 M, off-corner direction (theta = 0.271), psi_i = -0.10."""
    kw = C.har_kw(p_min="default")
    M, N, rho = kw["M"], kw["N"], kw["rho"]
    sig, th, q = C.sigma_on_surface(M, N, rho, -0.6, 1.2 * M, C.OFFC_DIR)
    pi = C.pi_S53(M, N, -0.6, 1.2 * M)
    v = C.v_for_psi(kw, pi, -0.10)
    P = C.make_hf("O2", **kw)
    st0 = O2.initial_state(P, np.diag(sig), v, pi)
    return P, st0, dict(th=th, q=q, pi=pi, v=v, sig=sig, M=M, N=N, rho=rho)


def inc(tr, sh):
    """(2e-5/3) 1 + 2e-5 n_hat: tr d eps = `tr`, shear amplitude ||dev|| = `sh` along the state's n_hat (diagonal tensor)."""
    return np.diag(tr / 3.0 * ONES + sh * C.OFFC_DIR)


def test_k1_14b_printed_state_and_the_fpf_step_every_printed_intermediate():
    """13.14b, the dry-side post floor under HAR.  Premise (printed): p = -0.6, q = 0.6999, eta = 1.2 M = 1.5971, theta = 0.271, pi_i =
    -0.74366, v = 1.72704 (psi_i = -0.10), eps_v = 9.851e-4, eps_s = 3.335e-5, D11 = 13 525, D12 = -6 235, D22 = 28 255 kPa.  Step
    d eps = (2e-5/3) 1 + 2e-5 n_hat: pattern FPf; p^tr = -0.4607 -> floored trial d eps^f_v,tr = 3.82e-6 -> d lambda = 1.157e-5, eta_c =
    1.8200, pi_i = -0.74335, d eps^p_v = +7.07e-6 (dilative), d eps^p_s = 1.57e-5, p_c = -0.4877 (> -p_min: the return ends ABOVE
    the floor), post floor p = -0.505, q 0.6619 -> 0.6705, F(sigma_c)/p_min = 6e-13 (compared < 1e-9), F(sigma_f)/p_min = -0.0041.
    Counters (n_f_tr, n_f_post) = (1, 1).  Computed here from the common closed forms on the O2 pre-floor events, not from O2's p.
    KILLS: M-F9 (floor inside the Newton: d lambda, pi_i, q_c differ and the FPf pattern itself disappears), M-F1/M-F5 (no trial
    floor: p^tr stays -0.46), M-F2, M-F7 (pi_i altered by the projection), the BA06 sign rule applied under HAR ('dry side FP-')."""
    P, st0, I = fpf_build()
    M, N, rho = I["M"], I["N"], I["rho"]
    p0, q0, th0 = C.invariants_pq(st0.sigma)
    ev0, es0, _ = inv_e(st0)
    assert abs(p0 + 0.6) <= 1e-12 and abs(q0 - 0.6999) <= 5e-5 and abs(th0 - 0.271) <= 5e-4 and abs(I["pi"] + 0.74366) <= 5e-6
    assert abs(I["v"] - 1.72704) <= 5e-6 and abs(ev0 - 9.851e-4) <= 6e-8 and abs(es0 - 3.335e-5) <= 6e-9
    assert abs(C.F_yield(M, N, rho, st0.sigma, st0.pi_i)) <= 1e-12
    D11, D12, D22, _ = C.har_D(p0, q0)
    assert abs(D11 - 13525) <= 1 and abs(D12 + 6235) <= 1 and abs(D22 - 28255) <= 1
    deps = inc(2e-5, 2e-5)
    # the raw trial pressure and the first floor event
    ev_t, es_t, _ = C.invariants_eps(st0.eps_e + deps)
    p_tr, _ = C.har_pq(ev_t, es_t)
    assert abs(p_tr + 0.4607) <= 6e-5, p_tr
    st = O2.step(P, st0, deps)
    assert st.flags["fpattern"] == "FPf" and st.flags["plastic"], st.flags["fpattern"]
    assert (st.n_f_tr, st.n_f_post, st.n_f_init) == (1, 1, 0)
    ev1, ev2 = st.cache["floor_events"]
    assert abs(ev1[2] - 3.82e-6) <= 6e-9, ev1[2]                                                     # d eps^f_v,tr
    pre = ev2[0]
    evc, esc, _ = C.invariants_eps(pre)
    p_c, q_c = C.har_pq(evc, esc)
    assert abs(p_c + 0.4877) <= 6e-5 and p_c > -C.PMIN_T, p_c
    assert abs(q_c - 0.6619) <= 6e-5
    sig_c = C.sigma_pq(p_c, q_c, C.invariants_eps(pre)[2])          # the converged stress, co-axial with the pre-floor strain
    p_f, q_f, _ = C.invariants_pq(st.sigma)
    assert abs(p_f + C.PMIN_T) <= 1e-13 and abs(q_f - 0.6705) <= 6e-5
    assert abs(st.dlam - 1.157e-5) <= 1e-3 * 1.157e-5 + 6e-9, st.dlam
    assert abs(st.eta - 1.8200) <= 6e-5 and abs(st.pi_i + 0.74335) <= 6e-6
    # printed d eps^p_v = +7.07e-6, d eps^p_s = 1.57e-5 are the changes of the ELASTIC strain invariants over the return (floored trial ->
    # converged, pre-floor), -(ev_c - ev_trf) and -(es_c - es_trf); the accumulated plastic volumetric strain of the state agrees
    evt, est_, _ = C.invariants_eps(ev1[1])
    assert abs(-(evc - evt) - 7.07e-6) <= 6e-9 and abs(-(esc - est_) - 1.57e-5) <= 6e-8 and -(evc - evt) > 0.0
    assert abs(st.eps_p_v - 7.07e-6) <= 6e-9 and st.eps_p_v > 0.0
    Fc = C.F_yield(M, N, rho, sig_c, st.pi_i) / C.PMIN_T
    assert abs(Fc) <= 1e-9, Fc                                                                       # printed 6e-13
    Ff = C.F_yield(M, N, rho, st.sigma, st.pi_i) / C.PMIN_T
    assert abs(Ff - (-0.0041)) <= 6e-5, Ff
    # the competing terms of the return (mid-point D's from the common closed form): D12 d eps^e_s beats |D11 d eps^e_v| (the HAR dry side)
    assert pressure(st) == pytest.approx(-C.PMIN_T, abs=1e-13)


def test_k1_14b_printed_pattern_map_around_the_state():
    """Sheet 9.7 item 2: 'the floored trial returns above the floor (FPf) for shear amplitudes up to 8e-5 at any expansion >= 2e-5 and
    below it (FP-) from 1.6e-4 on, where the larger d lambda lets the dilative volumetric term win; a pure expansion floors elastically
    (FE-)'.  Asserted (cells inside those statements; the full map is fpf_region_har.log): pure expansion tr in {2e-5 ... 3.2e-4} ->
    'FE-'; FPf for tr in {4e-5, 8e-5, 1.6e-4} x shear in {5e-6, 1e-5, 2e-5} and tr = 2e-5 x the same; FP- for tr >= 8e-5 with shear in
    {1.6e-4, 3.2e-4}.  From the evidence log (not the sheet text): tr = 1e-5 with shear >= 1e-5 is '-P-'.
    KILLS: a post-floor test that follows the BA06 sign rule (dry side never post-floors), a missing post floor under HAR, a floor
    applied inside the Newton (the FP- / FPf boundary moves)."""
    P, st0, _ = fpf_build()
    for tr in (2e-5, 4e-5, 8e-5, 1.6e-4, 3.2e-4):
        assert O2.step(P, st0, inc(tr, 0.0)).flags["fpattern"] == "FE-", tr
        for sh in (5e-6, 1e-5, 2e-5):
            assert O2.step(P, st0, inc(tr, sh)).flags["fpattern"] == "FPf", (tr, sh)
        if tr >= 8e-5:
            for sh in (1.6e-4, 3.2e-4):
                assert O2.step(P, st0, inc(tr, sh)).flags["fpattern"] == "FP-", (tr, sh)
    for sh in (1e-5, 2e-5, 4e-5):
        assert O2.step(P, st0, inc(1e-5, sh)).flags["fpattern"] == "-P-", sh


HAR_FD_H = (1e-7, 1e-8)


def test_k1_14b_fpf_tangent_matches_fd_o_h2():
    """(S.32f) at the dry-side FPf step: C_f = a^e(eps_f) Phi^post [b Phi^tr - u Pi_v v_{n+1} 1^T] against the central FD of the O2
    step over the three principal strains.  Printed: 2.2e-5 / 2.2e-7 / 2.2e-9 at h = 1e-6 / 1e-7 / 1e-8 (ratios 100, 101: O(h^2); the
    increment is 2e-5, so h = 1e-6 is 5 % of it); delta:C_f = 2e-16; committed p = -p_min to 2e-15.  Gate: 3x the printed value at each of
    h = 1e-7, 1e-8 (7e-7, 7e-9) and the O(h^2) ratio err(1e-7)/err(1e-8) in [30, 300].  Negative controls (printed): Phi^post dropped
    0.14, Phi^tr dropped 0.12, both dropped 2.2 -- the first two are computed here from the returned blocks (res.chain): the gate
    separates them by > 1e5.
    KILLS: M-F3a, M-F3b under HAR (eps'_f dropped: Phi without the n_hat term), M-F3d at m = 1, M-F9, the BA06 sign rule."""
    P, st0, _ = fpf_build()
    deps = inc(2e-5, 2e-5)
    stn, Cb, errs = best_fd_error(P, st0, deps, single, (1e-6,) + HAR_FD_H)
    assert stn.flags["fpattern"] == "FPf"
    assert errs[1e-7] <= 7e-7 and errs[1e-8] <= 7e-9, errs
    assert 30.0 <= errs[1e-7] / errs[1e-8] <= 300.0, f"not O(h^2): {errs}"
    assert np.abs(Cb.sum(axis=0)).max() <= 1e-13 * np.abs(Cb).max()
    assert abs(pressure(stn) + C.PMIN_T) <= 1e-13
    mu = mutants(stn)
    assert rel(mu["none"], Cb) >= 1.0, "both operators dropped must be ~2.2 off"          # printed 2.2
    assert 0.05 <= rel(mu["no_post"], Cb) <= 0.5 and 0.05 <= rel(mu["no_tr"], Cb) <= 0.5   # printed 0.14, 0.12


@pytest.mark.parametrize("name,deps_fn,fr,pattern,gates", [
    ("FPf,FPf", lambda: inc(4e-5, 4e-5), (0.5, 0.5), "FPf,FPf", {1e-7: 7e-7, 1e-8: 7e-9}),
    ("-P-,-Pf", lambda: inc(1.4e-5, 2e-5), (0.5, 0.5), "-P-,-Pf", {1e-7: 3.6e-6, 1e-8: 3.6e-8}),
])
def test_k1_14b_chained_patterns_match_fd(name, deps_fn, fr, pattern, gates):
    """(S.54) / (S.46f) under HAR, m = 2, fractions (1/2, 1/2) held fixed: twice the printed increment gives FPf,FPf (the printed
    increment itself halves to -P-,FPf), tr d eps = 1.4e-5 with shear 2e-5 gives -P-,-Pf (unsplit: -P- with p_c = -0.5072, p^tr =
    -0.5311).  Printed FD records: 2.3e-7 / 2.3e-9 (FPf,FPf) and 1.2e-6 / 1.2e-8 (-P-,-Pf) at h = 1e-7 / 1e-8; gate 3x each, O(h^2).
    Counters sum over the two sub-increments: FPf,FPf -> (n_f_tr, n_f_post) = (2, 2); -P-,-Pf -> (0, 1).
    KILLS: M-F3d (Phi operators omitted from the chain: 2.1 off), a chain that applies the post operator of the wrong sub-increment,
    counters not summed, the floored state of sub-increment 1 not carried into sub-increment 2."""
    P, st0, _ = fpf_build()
    deps = deps_fn()
    stn, Cb, errs = best_fd_error(P, st0, deps, fractions(fr), (1e-7, 1e-8))
    assert stn.flags["fpattern"] == pattern, stn.flags["fpattern"]
    for h, g in gates.items():
        assert errs[h] <= g, (name, h, errs)
    assert 30.0 <= errs[1e-7] / errs[1e-8] <= 300.0, f"not O(h^2): {errs}"
    want = (2, 2) if name == "FPf,FPf" else (0, 1)
    assert (stn.n_f_tr, stn.n_f_post) == want
    last = np.array([[O2api.tangent_last_substep(P, stn)[a, a, b, b] for b in range(3)] for a in range(3)])
    assert rel(last, Cb) >= 1e3 * gates[1e-7], "the last-sub-increment tangent coincides with the chain: the gate has no power"
    if name == "FPf,FPf":
        assert abs(pressure(stn) + C.PMIN_T) <= 1e-12


def test_k1_14b_unsplit_and_halved_patterns_printed():
    """Sheet 13.14b: twice the increment unsplit is FPf; the printed increment halves to -P-,FPf (m = 2); tr 1.4e-5 + shear 2e-5 is
    unsplit '-P-' with p_c = -0.5072 and p^tr(floored-not) = -0.5311.  KILLS: a ladder or fractions path that skips the floor test of
    a sub-increment."""
    P, st0, _ = fpf_build()
    assert O2.step(P, st0, inc(4e-5, 4e-5)).flags["fpattern"] == "FPf"
    assert O2.step_fractions(P, st0, inc(2e-5, 2e-5), (0.5, 0.5)).flags["fpattern"] == "-P-,FPf"
    s = O2.step(P, st0, inc(1.4e-5, 2e-5))
    assert s.flags["fpattern"] == "-P-"
    assert abs(pressure(s) + 0.5072) <= 6e-5
    ev_t, es_t, _ = C.invariants_eps(st0.eps_e + inc(1.4e-5, 2e-5))
    assert abs(C.har_pq(ev_t, es_t)[0] + 0.5311) <= 6e-5


# ======================================================================================================================
# counters, initial projection, census
# ======================================================================================================================
def test_initial_state_is_projected_counted_and_the_stress_replaced_ba06_and_har(oracle_name):
    """9.7 `initialState`: eps^e := Pi_f(invert(sigma0)), counted (n_f_init = 1, eps^f_v, W_f = p_min eps^f_v), and the deck's sigma0 is
    REPLACED by sigma(eps^e_f): p = -p_min; under BA06 alpha0 = 0 q is unchanged (q = 3 mu0 eps_s), under HAR q is scaled by
    (varpi_f/varpi)^n.  BA06 K2: sigma0 = isotropic -0.2 -> p = -0.5, eps^f_v = kappa ln(0.5/0.2) = 9.162907319e-3.  HAR TIMs: (p, q) =
    (-0.3, 0.5) -> p = -0.505, q_f = q (varpi_f/varpi)^n with varpi from (S.5h'), eps_s unchanged, eps^f_v = eps_v - eps_v,f.
    KILLS: an initial state left at p = -0.2 (not floored: a deck's first row of Gauss points would sit below p_min), an uncounted
    projection, q not rescaled under HAR (the BA06 rule 'keep q' applied to HAR), a refusal instead of a projection."""
    P = C.make_hf(oracle_name, **C.ba06_kw(p_min=PMIN_K2))
    st = ORACLES[oracle_name].initial_state(P, iso(-0.2), 1.70, -5000.0)
    assert abs(pressure(st) + 0.5) <= 1e-13 * 0.5
    want = K2_KH * math.log(0.5 / 0.2)
    efv, Wf, atf = f_counters(oracle_name, st)
    n_init = st.n_f_init if oracle_name == "O2" else st.flags["n_f_init"]
    assert n_init == 1 and abs(efv - want) <= 1e-12 and abs(Wf - 0.5 * want) <= 1e-12 and atf
    # HAR with shear
    nh = C.dir_for_theta(0.5)
    p0_, q0_ = -0.3, 0.5
    kw = C.har_kw(p_min="default")
    P = C.make_hf(oracle_name, **kw)
    st = ORACLES[oracle_name].initial_state(P, np.diag(C.sigma_pq(p0_, q0_, nh)), 1.70, -5000.0)
    ev_i, es_i, vp_i = C.har_inverse(p0_, q0_)
    f = C.har_floor(es_i, C.PMIN_T)
    q_f = q0_ * (f["x"] * C.PA_T / vp_i) ** C.N_T
    p_f, q_c, _ = C.invariants_pq(st.sigma)
    assert abs(p_f + C.PMIN_T) <= 1e-12 and abs(q_c - q_f) <= 1e-10 * q_f, (q_c, q_f)
    assert abs(C.invariants_eps(st.eps_e)[1] - es_i) <= 1e-13
    efv, Wf, atf = f_counters(oracle_name, st)
    n_init = st.n_f_init if oracle_name == "O2" else st.flags["n_f_init"]
    assert n_init == 1 and abs(efv - (ev_i - f["ev_f"])) <= 1e-12 and abs(Wf - C.PMIN_T * efv) <= 1e-15 and atf
    # an initial state ABOVE the floor is not projected
    st = ORACLES[oracle_name].initial_state(P, iso(-5.0), 1.70, -5000.0)
    assert abs(pressure(st) + 5.0) <= 1e-12
    assert (st.n_f_init if oracle_name == "O2" else st.flags.get("n_f_init", 0)) == 0


def test_substepped_floor_events_are_summed_and_the_projection_is_deterministic():
    """9.7 counting: a substepped increment sums its sub-increments' events.  The K1.12 increment taken as 2 and 4 equal sub-increments:
    first sub-increment floors at the trial, later ones are on the floor, so n_f_tr is the number of sub-increments that START from or
    reach p > -p_min(1 - 1e-12) (every one after the first reaches it again: the trial of an expansion at the floor is above it), the
    final state equals the unsplit one (isotropic, volumetric: the floor is path independent for alpha0 = 0) and eps^f_v = kappa ln 2
    + ...  Gate: final p = -0.5, total eps^f_v = tr d eps - (eps_v,f - eps_v0) = kappa ln 2 for any m, counters >= 1 and W_f = p_min eps^f_v.
    KILLS: an event count not summed over the sub-increments, W_f not cumulative, a floor path-dependent for the isotropic axis."""
    P, st0, deps = k112_setup("O2", shear=False)
    ref = O2.step(P, st0, deps)
    for m in (2, 4):
        st = O2.step_fractions(P, st0, deps, tuple([1.0 / m] * m))
        assert abs(pressure(st) + PMIN_K2) <= 1e-13 and abs(st.eps_f_v - ref.eps_f_v) <= 1e-13
        assert st.n_f_tr >= 1 and st.flags["floor_tr"] == st.n_f_tr and abs(st.W_f - PMIN_K2 * st.eps_f_v) <= 1e-15
        assert np.abs(st.sigma - ref.sigma).max() <= 1e-12 and abs(st.v - ref.v) <= 1e-13


@pytest.mark.parametrize("kind", ["ba06", "har"])
def test_floor_census_never_refuses_for_the_floor_keeps_p_below_pmin_counts_and_keeps_dissipation_nonnegative(kind):
    """9.7 contract on 150 seeded random increments from a near-surface start: (i) every ACCEPTED committed state has p <= -p_min (1 +
    1e-10) -- the floor's invariant; (ii) no refusal carries a floor-related reason ('p_or_pi_nonneg', 'elastic_domain', 'trial_elastic'
    ... ): the floor introduces no refusal; (iii) the counters never decrease, eps^f_v >= 0 and W_f = p_min eps^f_v exactly; (iv) D >= 0
    at every plastic step ('D >= 0 holds for the plastic flow, and the floor's violation is exactly E_f', sheet 9.7 item 1); (v)
    E_f of every step in [0, p_min d eps^f_v of that step].  BA06: K2, p_min = 0.5, start p = -2; HAR: TIMs, p_min = 0.505, start p =
    -1.5, increments up to 6e-4 / 1.5e-4 per component with an expansion bias so the floor engages (the census asserts it engages on
    >= 10 steps).
    KILLS: M-F1 (a refusal where a floor event is due, a committed p above -p_min), M-F2, D leaking negative work through the floor, a
    floor that fires on the wrong side (p more compressive than p_min)."""
    rng = np.random.default_rng(1443 if kind == "ba06" else 1444)
    if kind == "ba06":
        kw = C.ba06_kw(p_min="default")
        amp, bias, pmin = 6e-3, 1.5e-3, 0.5
        P, st = make_state("O2", kw, iso(-2.0), -3.0, 1.70)
    else:
        kw = C.har_kw(p_min="default")
        amp, bias, pmin = 1.5e-4, 3e-5, C.PMIN_T
        P, st = make_state("O2", kw, iso(-1.5), -2.2, C.v_for_psi(kw, -2.2, -0.05))
    n_f = n_pl = n_ref = 0
    prev_ntr = prev_npost = 0
    prev_efv = st.eps_f_v
    for _ in range(150):
        a = rng.uniform(-amp, amp, size=(3, 3))
        d = 0.5 * (a + a.T) + bias * I3 * rng.uniform(-1.0, 1.0)
        nxt = O2.step(P, st, d)
        if nxt.flags["refused"]:
            n_ref += 1
            reason = str(nxt.flags["reason"])
            assert not any(s in reason for s in ("p_or_pi_nonneg", "elastic_domain", "trial_elastic", "p_nonneg")), reason
            continue
        assert pressure(nxt) <= -pmin * (1.0 - 1e-10), f"committed state above the floor: p = {pressure(nxt)}"
        assert nxt.n_f_tr >= prev_ntr and nxt.n_f_post >= prev_npost and nxt.eps_f_v >= prev_efv - 1e-18
        assert abs(nxt.W_f - pmin * nxt.eps_f_v) <= 1e-15 + 1e-12 * nxt.W_f
        step_ev = nxt.eps_f_v - prev_efv
        Ef = O2api.floor_energy(P, nxt)
        assert -1e-18 <= Ef <= pmin * step_ev * (1.0 + 1e-9) + 1e-18, (Ef, pmin * step_ev)
        if nxt.flags["plastic"]:
            n_pl += 1
            assert nxt.D >= 0.0, f"negative dissipation {nxt.D}"
        n_f += int(nxt.flags["floor_tr"] + nxt.flags["floor_post"] > 0)
        prev_ntr, prev_npost, prev_efv = nxt.n_f_tr, nxt.n_f_post, nxt.eps_f_v
        st = nxt
    assert n_f >= 10, f"the floor engaged on only {n_f} of the accepted steps: the census proves little"
    assert n_pl >= 10


@pytest.mark.sheet_check
@pytest.mark.parametrize("energy", ["ba06", "har"])
def test_floor_energy_is_the_psi_difference_equals_the_quadrature_and_lies_in_zero_to_Wf(energy):
    """(S.52): E_f = Psi(eps_f) - Psi(eps) = integral from eps_v,f to eps_v of |p(eps'_v, eps_s)| d eps'_v in [0, p_min d eps^f_v] for a
    pre-floor state in dom Psi; BA06 alpha0 = 0 exactly E_f = kappa (p_min - |p|) with d eps^f_v = kappa ln(p_min/|p|).  A grid of 4
    eps_s x 5 |p|/p_min (t = 0.1 .. 0.999): the Psi difference of the common closed forms equals a 40-point Gauss-Legendre quadrature of
    |p| (1e-10 relative), the bound holds with margin >= 5e-4 at t = 0.999, and the printed ranges hold: E_f / (p_min d eps^f_v) in
    [0.21, 1) (BA06) and [0.37, 1) (HAR TIMs).  A sheet check: it exercises the closed forms the oracle's floor_energy is compared to.
    KILLS (as a reference): a Psi that is not the potential of (p, q) (the quadrature would differ), a sign error of E_f."""
    xg, wg = np.polynomial.legendre.leggauss(40)
    pmin = 0.5 if energy == "ba06" else C.PMIN_T
    lo_ratio = 0.21 if energy == "ba06" else 0.37
    for es in ((0.0, 5e-3, 1e-2, 3e-2) if energy == "ba06" else (0.0, 5e-5, 2e-4, 5e-4)):
        if energy == "ba06":
            ev_f = C.ba06_floor(es, pmin)["ev_f"]
            pfun = lambda ev: abs(C.ba06_pq(ev, es)[0])           # noqa: E731
            psi = lambda ev: C.ba06_psi(ev, es)                   # noqa: E731
        else:
            ev_f = C.har_floor(es, pmin)["ev_f"]
            pfun = lambda ev: abs(C.har_pq(ev, es)[0])            # noqa: E731
            psi = lambda ev: C.har_psi(ev, es)                    # noqa: E731
        for t in (0.1, 0.3, 0.6, 0.9, 0.999):
            # pre-floor state with |p| = t p_min at this eps_s: bisection on the forward map (p is monotone in eps_v)
            a, b = ev_f, (ev_f + 0.05 if energy == "ba06" else C.EDGE_T - 1e-12)
            for _ in range(200):
                mid = 0.5 * (a + b)
                if pfun(mid) > t * pmin:
                    a = mid
                else:
                    b = mid
            ev = 0.5 * (a + b)
            assert abs(pfun(ev) - t * pmin) <= 1e-9 * pmin
            dfv = ev - ev_f
            Ef = psi(ev_f) - psi(ev)
            quad = 0.5 * (ev - ev_f) * sum(w * pfun(0.5 * (ev - ev_f) * x + 0.5 * (ev + ev_f)) for x, w in zip(xg, wg))
            assert abs(Ef - quad) <= 1e-9 * abs(quad) + 1e-16, (es, t, Ef, quad)
            assert 0.0 <= Ef < pmin * dfv
            assert Ef / (pmin * dfv) >= lo_ratio - 1e-3, (es, t, Ef / (pmin * dfv))
            if energy == "ba06" and es == 0.0:
                assert abs(Ef - K2_KH * (pmin - t * pmin)) <= 1e-12 and abs(dfv - K2_KH * math.log(1.0 / t)) <= 1e-12
            if t == 0.999:
                margin = 1.0 - Ef / (pmin * dfv)
                assert margin > 0.0
                if energy == "ba06":
                    assert abs(margin - 5e-4) <= 2e-5, margin                 # printed 'margin 5e-4 at t = 0.999'


# ======================================================================================================================
# O2 -> O1 first order on a floored path (S.55)
# ======================================================================================================================
def floored_path_setup(oname, a0):
    """The sheet's floor_rate.py path: K2 BA06 (alpha0 = a0), p_min = 50, state ON the floor and the surface (p = -50, eta = 0.6 M, wet
    side, off-corner direction (-1.2, -0.1, 1.3)), loading along the stress deviator with slight expansion: d eps = (n_hat + 0.02) T,
    T = 2e-3 -- both mechanisms (the NorSand flow and the floor) stay active (sheet: min lambda 0.38, min lambda_f 0.20)."""
    kw = C.ba06_kw(p_min=PMIN_FD, alpha0=a0)
    P = C.make_hf(oname, **kw)
    sig, th, q = C.sigma_on_surface(1.2, 0.4, 0.7, -PMIN_FD, 0.6 * 1.2, C.WET_DIR)
    pi0 = C.pi_S53(1.2, 0.4, -PMIN_FD, 0.6 * 1.2)
    st0 = ORACLES[oname].initial_state(P, np.diag(sig), 1.70, pi0)
    edot = C.WET_DIR + 0.02
    return P, st0, np.diag(edot * 2.0e-3)


@pytest.mark.parametrize("a0", [0.0, 5.0])
def test_o2_converges_to_o1_at_first_order_on_a_floored_path(a0):
    """Plan 5.1 / 16.3 under the floor: the projection splitting O2 (trial floor -> return -> post floor) converges to the rate
    solution of the Koiter system (S.55) at FIRST order.  The sheet's printed ladder (floor_rate.py, Radau rtol 1e-10): m = 1, 2, 4, 8,
    16, 32, 64 sub-increments miss the rate solution by 7.2e-3 ... 1.4e-4 in sigma, 2.2e-2 ... 4.1e-4 in pi_i, observed orders 0.87, 0.92,
    0.96, 0.98, 0.99, 0.99 (alpha0 = 0; alpha0 = 5: 7.5e-3 ... 1.4e-4 and 2.4e-2 ... 4.5e-4).  Gates: the O1 truth (one increment)
    completes on the floor, err strictly decreasing, log-log slope of the last four ratios in [0.8, 1.3] (first order), err(m) within
    1.5x of the printed values at m = 1 and m = 64.
    KILLS: an O1 floor form (S.55) that is not the Koiter system (the gap does not close), a first-order scheme degraded to O(1) by a
    floor applied inside the Newton, a wrong E_f/W_f feed-back."""
    P1, st01, d = floored_path_setup("O1", a0)
    sts = O1.run_path(P1, st01, d[None], rtol=1e-10)
    assert sts[-1].flags["status"] == "ok", sts[-1].flags["status"]
    sig1, pi1 = sts[-1].sigma, sts[-1].pi_i
    P2, st02, _ = floored_path_setup("O2", a0)
    errs, pierrs = [], []
    ms = (1, 2, 4, 8, 16, 32, 64)
    for m in ms:
        out = O2.run_path(P2, st02, np.array([d / m] * m))
        assert not any(s.flags["refused"] for s in out), f"O2 m={m} refused"
        errs.append(float(np.abs(out[-1].sigma - sig1).max() / np.abs(sig1).max()))
        pierrs.append(abs(out[-1].pi_i - pi1) / abs(pi1))
    print(f"\n[floor conv a0={a0}] m = {ms}\n  sigma {['%.3e' % e for e in errs]}\n  pi_i  {['%.3e' % e for e in pierrs]}")
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), errs
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(2.0) for i in range(len(errs) - 1)]
    assert all(0.8 <= o <= 1.3 for o in orders[2:]), orders
    printed = {0.0: (7.2e-3, 1.4e-4, 2.2e-2, 4.1e-4), 5.0: (7.5e-3, 1.4e-4, 2.4e-2, 4.5e-4)}[a0]
    assert errs[0] <= 1.5 * printed[0] and errs[-1] <= 1.5 * printed[1], (errs[0], errs[-1])
    assert pierrs[0] <= 1.5 * printed[2] and pierrs[-1] <= 1.5 * printed[3], (pierrs[0], pierrs[-1])
