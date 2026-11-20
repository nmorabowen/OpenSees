"""WP-144 gate G1 (Zone B): Gudehus-Argyris plastic-path census, both oracles (G1 review item N4).

Every existing plastic-path census (K1.5-K1.9, I.1, I.3) runs on Willam-Warnke. Here the SAME three scenarios (drained TXC,
drained TXE, undrained TXC) run on the Gudehus-Argyris shape (S.10) with rho = 0.8, rho_bar = 0.85, in both CSL modes (paper,
fork), and the two path-wise identities that depend on the shape are asserted:

  GA.D   D >= b_min min(|p_n|, |p_{n+1}|) Delta lambda at every plastic increment (sheet 11.1/11.2, I.3 of test_g1_identities),
         b_min the minimum of the (S.38) bracket over the surface, both oracles;
  GA.H   Delta pi_i = sqrt(2/3) h Delta lambda (pi_i* - pi_i) Omega (S.26) at every unsubstepped plastic O2 step, and the
         pointwise rate law d pi_i / d eps^p_s = h (pi_i* - pi_i) for O1 (Richardson probe), drained paths, with Omega of (S.20)
         evaluated with the GA corner values (zeta_bar(0) = 1/rho_bar, zeta_bar(pi/3) = 1, zeta_bar_y = sqrt6 (1 - rho_bar)/(2 rho_bar)
         at both corners, S.9 with S.10).

TEST RULE (plan 5): every expected value and tolerance below is a closed form of the sheet (section cited) or a stated round-off /
convergence argument, written before any oracle output was read for these paths; no oracle output is an expected value. A failing
assertion is a finding, not a reason to edit an oracle or a bound. The tolerances are those of the WW tests they mirror (the
argument is the same and is not repeated beyond what changes).
"""
import functools
import math

import numpy as np
import pytest

from test_g1_k1_closed_forms import (CASES, FORK, PAPER, SQ23, SIG0, assert_path_ok, case_v0, init_state, is_plastic,  # noqa: F401
                                     ok, pq_theta, psi_of, run_staged, run_tx, why)
from test_g1_identities import pistar_closed

RHO, RHO_BAR = 0.8, 0.85            # GA range [7/9, 1] (sheet 4.1); rho/rho_bar = 0.941 >= beta = 0.75 (S.39): condition A holds
GA_PAPER = dict(PAPER, zeta="GA", rho=RHO, rho_bar=RHO_BAR)
GA_FORK = dict(FORK, zeta="GA", rho=RHO, rho_bar=RHO_BAR)

GA_CASES = {
    "TXC_paper": dict(CASES["TXC_paper"], kw=GA_PAPER),
    "TXC_fork": dict(CASES["TXC_fork"], kw=GA_FORK),
    "TXE_paper": dict(CASES["TXE_paper"], kw=GA_PAPER),
    "TXE_fork": dict(CASES["TXC_fork"], kw=GA_FORK, sign=+1, corner="E"),
    "UND_paper": dict(CASES["UND_paper"], kw=GA_PAPER),
    "UND_fork": dict(CASES["UND_fork"], kw=GA_FORK),
}
DRAINED = ["TXC_paper", "TXC_fork", "TXE_paper", "TXE_fork"]


def ga_zeta(theta, rho):
    """(S.10): zeta = [(1 + rho) + (1 - rho) cos 3 theta]/(2 rho)."""
    return ((1.0 + rho) + (1.0 - rho) * np.cos(3.0 * theta)) / (2.0 * rho)


def b_min_ga(kw):
    """Minimum over the surface of the bracket of (S.38), sheet 11.1: linear in eta on [0, M/N], so at an end point; zeta_bar/zeta is a
    Moebius function of cos 3 theta (monotone), minimum min(rho/rho_bar, 1) at theta = 0 or pi/3. The numeric min over a
    (eta x theta) grid of the exact GA zeta is asserted equal (1e-9) in test_ga_b_min_closed_form."""
    M, N, Nb = kw["M"], kw["N"], kw["N_bar"]
    r = min(kw["rho"] / kw["rho_bar"], 1.0)
    return min(M / (1.0 - Nb), (M - M / N) / (1.0 - Nb) + r * M / N)


def omega_corner_ga(kw, sig, corner):
    """Omega of (S.20) at an axisymmetric stress with the GA corner values: zeta_bar = 1 (C) or 1/rho_bar (E);
    zeta_bar_y = -sqrt6 zeta''(0)/9 = +sqrt6 (1 - rho_bar)/(2 rho_bar) at E, +sqrt6 zeta''(pi/3)/9 = the same at C (S.9, S.10)."""
    p, q, _, w = pq_theta(sig)
    rb = kw["rho_bar"]
    zb = 1.0 if corner == "C" else 1.0 / rb
    zby = math.sqrt(6.0) * (1.0 - rb) / (2.0 * rb)
    xi = w - p
    R = float(np.linalg.norm(xi))
    ya = 3.0 * xi ** 2 / R ** 3 - 3.0 * float((xi ** 3).sum()) * xi / R ** 5 - 1.0 / R
    return math.sqrt(1.5 * zb ** 2 + (zby * q) ** 2 * float((ya ** 2).sum()))


@functools.lru_cache(maxsize=None)
def run_ga(oname, case):
    spec = GA_CASES[case]
    P, st0 = init_state(oname, spec["kw"], SIG0, case_v0(spec), spec["pi0"])
    sts, cum = run_staged(oname, P, st0, spec["kind"], spec["sign"], spec["sched"])
    return P, spec["kw"], st0, sts, cum


def test_ga_b_min_closed_form():
    """Sheet 11.1: b_min = 0.5735 for GA 0.8 / 0.85 on the K2 constants; equals the numeric min of the bracket over a
    (eta x theta) grid of the exact GA zeta (both end points in the grid), 1e-9."""
    kw = GA_PAPER
    eta = np.linspace(0.0, kw["M"] / kw["N"], 401)[:, None]
    th = np.linspace(0.0, math.pi / 3.0, 301)[None, :]
    br = (kw["M"] - eta) / (1.0 - kw["N_bar"]) + ga_zeta(th, kw["rho_bar"]) / ga_zeta(th, kw["rho"]) * eta
    assert abs(float(br.min()) - b_min_ga(kw)) <= 1e-9


@pytest.mark.parametrize("case", sorted(GA_CASES))
def test_ga_dissipation_lower_bound_every_plastic_step(oracle_name, case):
    """GA.D. S.38 on the GA surface: D >= b_min min(|p_n|, |p_{n+1}|) Delta lambda on every plastic increment with Delta eps^p_s > 0,
    Delta lambda = Delta eps^p_s/(sqrt(2/3) Omega) (S.20) with Omega the GA corner value recomputed here (the paths are
    axisymmetric). Same argument and tolerances as I.3 (O2: exact at the converged state; O1: the min over the increment is taken
    over the end points, the slack of the bound covers the second-order difference); round-off allowance 1e-9 relative; D > 0
    separately; >= 8 plastic increments. b_min = b_min_ga (0.5735).
    Kills: the WW zeta used in place of the GA one in either oracle's flow or dissipation (the bracket differs at theta = 0 by
    zeta_bar/zeta of the wrong shape function); a dropped zeta_bar_y q term (Omega and the flow change)."""
    P, kw, st0, sts, cum = run_ga(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    corner = GA_CASES[case]["corner"]
    bmin = b_min_ga(kw)
    prev, n_chk = st0, 0
    for k, st in enumerate(sts):
        dEs = st.eps_p_s - prev.eps_p_s
        if is_plastic(st) and dEs > 1e-12:
            Om = omega_corner_ga(kw, st.sigma, corner)
            dlam = dEs / (SQ23 * Om)
            pmin = min(abs(pq_theta(prev.sigma)[0]), abs(pq_theta(st.sigma)[0]))
            bound = bmin * pmin * dlam
            assert st.D > 0.0, f"{case} step {k}: D = {st.D:.3e} on a plastic increment"
            assert st.D >= bound * (1.0 - 1e-9), \
                f"{case} step {k}: D = {st.D:.6e} < b_min |p| dlam = {bound:.6e} (ratio {st.D / bound:.4f})"
            n_chk += 1
        prev = st
    assert n_chk >= 8, f"only {n_chk} plastic increments with Delta eps^p_s > 0"


def _hardening_gap(kw, st, corner):
    p = pq_theta(st.sigma)[0]
    Om = omega_corner_ga(kw, st.sigma, corner)
    psi = psi_of(kw, st.v, st.pi_i)
    return pistar_closed(kw, p, Om, psi) - st.pi_i, Om


@pytest.mark.parametrize("case", sorted(GA_CASES))
def test_ga_hardening_identity_o2_every_plastic_step(case):
    """GA.H (O2). (S.26), backward Euler: Delta pi_i = sqrt(2/3) h Delta lambda (pi_i* - pi_i)_{n+1} Omega_{n+1} per unsubstepped
    plastic step with q > 0; Delta lambda = Delta eps^p_s/(sqrt(2/3) Omega) with the GA Omega; pi_i* from (S.23) via the inverse of
    S.12 (pistar_closed of the I.1 test), psi_i from S.22 in the case's CSL mode. Tolerance 1e-8 |pred| + 1e-11 |pi_{i,n}|
    (nested solve |r| <= 1e-12 |pi_{i,n}|, as I.1). >= 20 steps checked and >= 8 with |pred| >= 1e-4 |pi_i| (signal).
    Kills: the WW Omega (no zeta_bar_y q term of the GA surface, or zeta_bar(0) = 1/rho_bar missing) inside the hardening rate."""
    P, kw, st0, sts, cum = run_ga("O2", case)
    assert_path_ok("O2", sts, cum, case)
    corner = GA_CASES[case]["corner"]
    prev, n_chk, n_signal = st0, 0, 0
    for k, st in enumerate(sts):
        if is_plastic(st) and pq_theta(st.sigma)[1] > 0.0 and st.flags.get("substeps", 1) == 1:
            dEs = st.eps_p_s - prev.eps_p_s
            assert dEs > 0.0, (k, dEs)
            a, Om = _hardening_gap(kw, st, corner)
            dlam = dEs / (SQ23 * Om)
            pred = SQ23 * kw["h"] * dlam * a * Om
            obs = st.pi_i - prev.pi_i
            tol = 1e-8 * abs(pred) + 1e-11 * abs(prev.pi_i)
            assert abs(obs - pred) <= tol, \
                f"{case} step {k}: d pi_i = {obs:.12e}, identity {pred:.12e}, miss {abs(obs - pred):.2e} > {tol:.2e}"
            n_chk += 1
            n_signal += abs(pred) >= 1e-4 * abs(prev.pi_i)
        prev = st
    assert n_chk >= 20, f"only {n_chk} unsubstepped plastic steps checked"
    assert n_signal >= 8, f"only {n_signal} steps with |d pi_i| >= 1e-4 |pi_i|"


@pytest.mark.parametrize("case", DRAINED)
def test_ga_hardening_identity_o1_pointwise_richardson(case):
    """GA.H (O1), drained paths (the drained set of I.1; the undrained paths are covered per step in O2 only). The rate law
    d pi_i / d eps^p_s = h (pi_i* - pi_i) holds pointwise; probed exactly as in I.1: chords of two tiny drained-triaxial increments
    (2e-6, 4e-6 axial), point value 2 r(d/2) - r(d), expected value h (pi_i* - pi_i) at the start state with the GA Omega; tolerance
    1e-4 |expected| + 0.05 (argument as I.1); every plastic state with |expected| >= 5 probed, >= 8 of them."""
    P, kw, st0, sts, cum = run_ga("O1", case)
    assert_path_ok("O1", sts, cum, case)
    spec = GA_CASES[case]
    n_chk = 0
    for k, st in enumerate(sts):
        if not (is_plastic(st) and pq_theta(st.sigma)[1] > 0.0):
            continue
        a, _ = _hardening_gap(kw, st, spec["corner"])
        expected = kw["h"] * a
        if abs(expected) < 5.0:
            continue

        def chord(d):
            s2 = run_tx("O1", P, st, spec["kind"], spec["sign"] * d, 1)[-1]
            assert ok("O1", s2) and is_plastic(s2), (k, d, why("O1", s2))
            return (s2.pi_i - st.pi_i) / (s2.eps_p_s - st.eps_p_s)

        point = 2.0 * chord(2e-6) - chord(4e-6)
        tol = 1e-4 * abs(expected) + 0.05
        assert abs(point - expected) <= tol, \
            f"{case} state {k}: point rate {point:.6f}, h(pi*-pi) = {expected:.6f}, tol {tol:.3g}"
        n_chk += 1
    assert n_chk >= 8, f"only {n_chk} states probed"
