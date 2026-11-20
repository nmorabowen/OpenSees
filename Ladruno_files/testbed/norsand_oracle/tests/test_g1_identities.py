"""WP-144 gate G1 (Zone B): four identity gates, both oracles.

  I.1  hardening-rate identity (S.25)/(S.26) and (S.40), sheet 8 -- TXC paper, TXC fork, TXE paper
  I.2  the vertex rule of sheet 3.2 as a POSITIVE test: pure isotropic compression past pi_c (K2 set, cap none and planar)
  I.3  lower bound on the plastic dissipation: D >= b_min |p| Delta lambda with b_min the minimum printed in sheet 11.2
  I.4  elastic-energy convexity table (plan 2.5 / 6.2 G1 line; sheet 2.2, 2.3, S.5), alpha0 in {0, 0.2}

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form in the equation sheet
144a (the section is cited at each use), (b) a number printed in the sheet, or (c) a stated convergence argument.  All of it
was written BEFORE the oracles were run on these paths.  No oracle output is used as an expected value.  A tolerance is
never loosened to make a test pass; a failing assertion is a finding (against the sheet or an oracle), not a reason to edit
an oracle or a bound.

How the two oracles are read (same adapter as test_g1_k1_closed_forms; conftest hides the naming differences):
  * stress p, q, theta are recomputed HERE from State.sigma; pi_i, v, eps_p_v, eps_p_s, D are State fields common to both;
  * O2 is backward Euler: a step quantity is a point value at the converged state and is EQUATED to the sheet's closed form;
  * O1 is the continuum rate integrator: a step ratio over an increment is an eps^p_s-weighted mean of the point value
    (sheet 13.5), so O1 is read POINTWISE by Richardson probing, never by chord.
Helpers (paths, closed forms, corner Omega) are imported from test_g1_k1_closed_forms, which is not edited.
"""
import functools
import math
import os

import numpy as np
import pytest

from conftest import K2_BASE, ORACLES, make_params
from test_g1_k1_closed_forms import (CASES, PAPER, PAPER_2INV, SQ23, SQ32, assert_path_ok, ok, omega_corner,  # noqa: F401
                                     pq_theta, psi_of, run_case, run_deps, run_tx, why, is_plastic)

SQ2 = math.sqrt(2.0)
I3 = np.eye(3)
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


# ==============================================================================================
# closed forms from the sheet (written here, not taken from the oracles)
# ==============================================================================================
def pistar_closed(kw, p, Om, psi):
    """pi_i* of (S.23), computed the way the gate text says: D* = chi psi_i gives the limit stress ratio
    eta* = M + sqrt(2/3) chi (1 - N_bar) psi_i Omega, and pi_i* is the inverse of S.12 at (p, eta*)
    (BA06 2.8, sheet 5.1: pi/p = [(1-N)/(1 - eta N/M)]^((1-N)/N), or exp(eta/M - 1) for N = 0).
    The power form B^((N-1)/N) printed in (S.23), B = 1 - sqrt(2/3) chi_bar psi Omega N / M, chi_bar = chi/beta, is computed
    as well and the two are asserted equal (an algebra check of the sheet, 1e-12 relative)."""
    M, N, Nb, chi = kw["M"], kw["N"], kw["N_bar"], kw["chi"]
    eta_star = M + SQ23 * chi * (1.0 - Nb) * psi * Om
    if N == 0.0:
        inv = p * math.exp(eta_star / M - 1.0)
        direct = p * math.exp(SQ23 * chi * (1.0 - Nb) * psi * Om / M)
    else:
        d = 1.0 - eta_star * N / M
        assert d > 0.0, f"B <= 0 guard of S.23: base {d}"
        inv = p * ((1.0 - N) / d) ** ((1.0 - N) / N)
        beta = (1.0 - N) / (1.0 - Nb)
        B = 1.0 - SQ23 * (chi / beta) * psi * Om * N / M
        direct = p * B ** ((N - 1.0) / N)
    assert abs(inv - direct) <= 1e-12 * abs(p), (inv, direct)
    return inv


def hardening_state_quantities(kw, st, corner):
    """(pi_i* - pi_i, Omega) at a state ON a corner path: p from the stress, Omega from (S.20) (corner values, S.9),
    psi_i from the state's v and pi_i (S.22), pi_i* from (S.23)."""
    p = pq_theta(st.sigma)[0]
    Om = omega_corner(kw, st.sigma, corner)
    psi = psi_of(kw, st.v, st.pi_i)
    return pistar_closed(kw, p, Om, psi) - st.pi_i, Om


# ==============================================================================================
# I.1  hardening-rate identity (sheet 8, S.25/S.26, S.40)
# ==============================================================================================
H_CASES = ["TXC_paper", "TXC_fork", "TXE_paper"]


@pytest.mark.parametrize("case", H_CASES)
def test_i1_hardening_identity_o2_every_plastic_step(case):
    """Sheet 8 (S.26), backward Euler: pi_i = pi_{i,n} + sqrt(2/3) h Dlam (pi_i* - pi_i) Omega, all of pi_i*, Omega at n+1.
    With Dlam = Deps^p_s/(sqrt(2/3) Omega) (S.20) this reads, per plastic step,
        Delta pi_i = sqrt(2/3) h Dlam (pi_i* - pi_i)_{n+1} Omega_{n+1},   to 1e-8 relative.
    Everything on the right is recomputed here: Omega from S.20 (corner values of zeta_bar, zeta_bar_y, S.9: the drained
    triaxial paths are axisymmetric), psi_i from S.22, pi_i* from S.23 (via the inverse of S.12), h from the parameter set;
    Deps^p_s and pi_i are the oracle's State increments.  Tolerance argument: the nested solve of (S.27) converges to
    |r| <= 1e-12 |pi_{i,n}| (README), so the identity holds to ~1e-12 |pi_{i,n}| absolute; the tolerance is
    1e-8 |Delta pi_i| + 1e-11 |pi_{i,n}| (the floor only matters at the steps where H -> 0 and Delta pi_i -> 0).
    Steps that the oracle substepped are skipped (the identity is per backward-Euler step); at least 20 must remain.
    MUTANTS KILLED: h scaled by 2 in the oracle (the identity misses by 100 %); Omega dropped from or doubled in the hardening
    rate; pi_i* evaluated with chi instead of chi_bar or at the start of the step (O(step) = 1e-2 miss); the CSL psi_i
    evaluated at pi_{i,n} instead of pi_{i,n+1}."""
    P, kw, st0, sts, cum = run_case("O2", case)
    assert_path_ok("O2", sts, cum, case)
    corner = CASES[case]["corner"]
    prev, n_chk, n_skip, n_signal, worst = st0, 0, 0, 0, 0.0
    for k, st in enumerate(sts):
        if is_plastic(st) and pq_theta(st.sigma)[1] > 0.0:
            if st.flags.get("substeps", 1) != 1:
                n_skip += 1
            else:
                dEs = st.eps_p_s - prev.eps_p_s
                assert dEs > 0.0, (k, dEs)
                a, Om = hardening_state_quantities(kw, st, corner)
                dlam = dEs / (SQ23 * Om)
                pred = SQ23 * kw["h"] * dlam * a * Om
                obs = st.pi_i - prev.pi_i
                err = abs(obs - pred)
                tol = 1e-8 * abs(pred) + 1e-11 * abs(prev.pi_i)
                assert err <= tol, f"{case} step {k}: d pi_i = {obs:.12e}, identity {pred:.12e}, miss {err:.2e} > {tol:.2e}"
                worst = max(worst, err / max(abs(pred), 1e-300))
                n_chk += 1
                if abs(pred) >= 1e-4 * abs(prev.pi_i):
                    n_signal += 1
        prev = st
    assert n_chk >= 20, f"only {n_chk} unsubstepped plastic steps checked ({n_skip} substepped)"
    assert n_signal >= 8, f"only {n_signal} steps with |d pi_i| >= 1e-4 |pi_i|: the identity is not exercised"


@pytest.mark.parametrize("case", H_CASES)
def test_i1_hardening_identity_o1_pointwise_richardson(case):
    """Sheet 8 (S.25)/(S.40) and 13.5 (revised): the rate law d pi_i / d eps^p_s = h (pi_i* - pi_i) holds POINTWISE.
    Over a finite increment the ratio Delta pi_i / Delta eps^p_s is an eps^p_s-weighted mean of the point value, so it is probed:
    from a plastic state, two tiny drained-triaxial increments d/2 = 2e-6 and d = 4e-6 (axial strain), chords r(d/2), r(d),
    point value = 2 r(d/2) - r(d) (first-order Richardson, chord error O(d) cancelled; sheet 13.5 measured 6.8e-10 at d = 2e-6).
    Expected value: h (pi_i* - pi_i) at the START state, with pi_i* from S.23 (Omega from S.20, psi_i from S.22) written here.
    Tolerance argument: the chord error is (1/2) d(ln r)/d eps^p_s Delta eps^p_s <= (1/2) h Delta eps^p_s ~ 4e-4 relative at
    d = 4e-6 and is O(d^2) ~ 2e-7 after Richardson; ODE noise: the oracle runs rtol 1e-10, so Delta pi_i is good to a few 1e-9
    absolute and Delta eps^p_s (~2e-6) to ~1e-10, i.e. the ratio to <~ 5e-2 absolute after the factor-3 Richardson
    amplification.  Tolerance 1e-4 |expected| + 0.05 (absolute in kPa per unit eps^p_s).  The hardening relaxes pi_i onto pi_i* at the rate
    h = 280 per unit eps^p_s, so only the early plastic states carry a signal: every plastic state with |expected| >= 5
    (|pi_i* - pi_i| >= 0.018 kPa; tolerance <= 1.1e-2 of the value, against 1.0 for the mutant below) is probed, at least 8 of them.
    MUTANTS KILLED: h scaled by 2 in the oracle (the point value misses by 100 %); Omega missing from the rate; chi_bar vs chi."""
    P, kw, st0, sts, cum = run_case("O1", case)
    assert_path_ok("O1", sts, cum, case)
    spec = CASES[case]
    corner = spec["corner"]
    n_pl, n_chk, worst = 0, 0, 0.0
    for k, st in enumerate(sts):
        if not (is_plastic(st) and pq_theta(st.sigma)[1] > 0.0):
            continue
        n_pl += 1
        a, _ = hardening_state_quantities(kw, st, corner)
        expected = kw["h"] * a
        if abs(expected) < 5.0:
            continue

        def chord(d):
            s2 = run_tx("O1", P, st, spec["kind"], spec["sign"] * d, 1)[-1]
            assert ok("O1", s2) and is_plastic(s2), (k, d, why("O1", s2))
            return (s2.pi_i - st.pi_i) / (s2.eps_p_s - st.eps_p_s)

        point = 2.0 * chord(2e-6) - chord(4e-6)
        tol = 1e-4 * abs(expected) + 0.05
        worst = max(worst, abs(point - expected) / abs(expected))
        assert abs(point - expected) <= tol, f"{case} state {k}: point rate {point:.6f}, h(pi*-pi) = {expected:.6f}, tol {tol:.3g}"
        n_chk += 1
    assert n_chk >= 8, f"only {n_chk} states probed (|h (pi* - pi)| >= 5)"


# ==============================================================================================
# I.2  vertex rule as a positive test (sheet 3.2; 10.1 for the planar cap)
# ==============================================================================================
VERTEX_CAPS = {"none": dict(cap="none"), "planar": dict(cap="planar", c1=0.10, c2=0.10)}   # BA06 chi_cap = 0.10 (sheet 10.1)
PI_C_CASES = {"apex_start": -100.0, "elastic_then_plastic": -120.0}     # target pi_c (kPa); the energy reference is p0 = -100
DELTA = 1.0e-3                    # isotropic strain increment: deps = -DELTA I/3, tr deps = -DELTA
N_VERTEX = 8
V0_K2 = 1.59                      # K2 initial specific volume (sheet 14)
TOL_P = {"O1": 1e-8, "O2": 1e-9}  # relative tolerance on p = pi_c, argued in the docstring of the test
TOL_Q = {"O1": 1e-9, "O2": 1e-12}  # q/|p| on the axis, argued in the docstring of the test


def pi_i_for_pic(kw, pic):
    """Inverse of the apex relation pi_c = pi_i/(1-N)^((1-N)/N) (sheet 5.1, AB06 p.1538)."""
    N = kw["N"]
    return pic * (1.0 - N) ** ((1.0 - N) / N)


def run_vertex(oname, cap, pic):
    kw = dict(PAPER, **VERTEX_CAPS[cap])
    P = make_params(oname, **kw)
    pi0 = pi_i_for_pic(kw, pic)
    st0 = ORACLES[oname].initial_state(P, -100.0 * I3, V0_K2, pi0)
    sts = run_deps(oname, P, st0, np.array([-DELTA / 3.0 * I3] * N_VERTEX))
    return kw, P, st0, sts, pi0


@pytest.mark.parametrize("cap", ["none", "planar"])
@pytest.mark.parametrize("case", ["apex_start", "elastic_then_plastic"])
def test_i2_isotropic_compression_past_pi_c_is_1d_axial_return_with_frozen_pi_i(oracle_name, case, cap):
    """Sheet 3.2 (vertex rule, revised): an isotropic strain path leaves q = 0 at every step; F := p eta there, so F = 0 means
    eta = 0, i.e. p = pi_c = pi_i/(1-N)^((1-N)/N) (S.12 apex); flow is purely volumetric, Omega := 0, so pi_i is NOT hardened
    and hydrostatic compression beyond pi_c is perfectly plastic.  Path: sigma0 = -100 I (energy reference p0 = -100, eps^e_v0 = 0),
    K2 set (rho 0.7 / rho_bar 0.8), pi_c = -100 (state AT the apex) or -120 (elastic first), deps = -DELTA I/3 for 8 steps,
    cap none and planar (c1 = c2 = 0.10).  Closed forms (all written here):
      p_k    = -min(|p0| exp(k DELTA/kappa), |pi_c|)           (K1.1 elastic branch, clipped at the apex)
      eps^p_v(k) = min(0, -k DELTA - eps^e_v(pi_c)),  eps^e_v(pi_c) = eps_v0 - kappa ln(pi_c/p0)     (K1.1; < 0 = compaction)
      D_step = pi_c Delta eps^p_v         (D = sigma:Delta eps^p with an isotropic eps^p and p = pi_c; sheet 3.2 prints 0.300 per
                                           tr = -3e-3 at pi_c = -100)
      q = 0, pi_i == pi_i0, eps^p_s == 0, deviatoric Delta eps^p == 0, plastic iff the trial |p| exceeds |pi_c|.
    Tolerances.  p: the backward-Euler Newton leaves |F| <= 1e-12 |p0| and dF/dp = M/(1-N) = 2 at the apex, so |p - pi_c|
    <~ 5e-11 kPa (relative 5e-13); 1e-9 relative (O2).  O1: F is held to the oracle's surface tolerance F_tol = 1e-8 M|p|
    (sheet 12), dF/dp = M/(1-N): 1e-8 relative.  q: the strain path is exactly isotropic, so O2 (algebraic, symmetric) can only
    break the isotropy by round-off: 1e-12 |p|.  O1 is an ODE solve whose deviatoric elastic strain is controlled only to its
    absolute tolerance atol_eps = 0.1 rtol min(kappa, |p|/(3 mu0)) = 6e-13 (rtol 1e-10, integrator.py), i.e. q <~ 3 mu0 sqrt(3/2)
    atol_eps ~ 1e-8 kPa = 1e-10 |p| at p = -100: 1e-9 |p| (factor 10).  pi_i: Omega = 0 makes the hardening rate exactly 0, so pi_i is
    bit-identical for O2; for O1 the Radau error test bounds the state error per internal step by atol + rtol |pi_i| with rtol = 1e-10
    and an increment takes ~5 internal steps: 1e-9 |pi_i| (factor 2).  eps^p_v: p pinned at pi_c fixes eps^e_v, error
    kappa x (relative error of p) -> 2 kappa TOL_P.  D: relative error of p times the strain: 10 TOL_P.
    MUTANTS KILLED: the vertex rule removed (Omega := sqrt(3/2) zeta_bar at the axis: pi_i hardens, eps^p_s > 0, q leaves 0);
    a return to the wrong point on the axis (not pi_c) or no return at all (elastic through the apex); a volumetric
    flow of the wrong size or sign (eps^p_v, D); hardening of pi_i by the volumetric strain (pi_i must stay frozen)."""
    pic = PI_C_CASES[case]
    kw, P, st0, sts, pi0 = run_vertex(oracle_name, cap, pic)
    kap, p0, ev0 = kw["kappa_hat"], kw["p0"], kw["eps_v0"]
    assert abs(pi0 / (1.0 - kw["N"]) ** ((1.0 - kw["N"]) / kw["N"]) - pic) <= 1e-12 * abs(pic)   # construction sanity (S.12 apex)
    ev_e_pic = ev0 - kap * math.log(pic / p0)                                                      # K1.1 at p = pi_c
    tolp = TOL_P[oracle_name]
    prev_epv, n_plastic = 0.0, 0
    for k, st in enumerate(sts, start=1):
        assert ok(oracle_name, st), (k, why(oracle_name, st))
        trial_abs = abs(p0) * math.exp(k * DELTA / kap)
        p_cf = -min(trial_abs, abs(pic))
        p, q, _, _ = pq_theta(st.sigma)
        assert abs(p - p_cf) <= tolp * abs(pic), f"step {k}: p = {p:.12f}, closed form {p_cf:.12f}"
        assert q <= TOL_Q[oracle_name] * abs(p), f"step {k}: q = {q:.3e} (isotropic path must stay on the axis)"
        epv_cf = min(0.0, -k * DELTA - ev_e_pic)                    # eps^p_v = tr eps - eps^e_v(pi_c) once yielded (compaction < 0)
        assert abs(st.eps_p_v - epv_cf) <= 2.0 * kap * tolp + 1e-14, f"step {k}: eps^p_v = {st.eps_p_v:.12e}, closed form {epv_cf:.12e}"
        plastic_cf = trial_abs > abs(pic) * (1.0 + 1e-7)
        elastic_cf = trial_abs < abs(pic) * (1.0 - 1e-7)
        if plastic_cf:
            n_plastic += 1
            assert is_plastic(st), f"step {k}: trial beyond pi_c but the step is elastic"
            if oracle_name == "O2":
                assert st.pi_i == pi0, f"step {k}: pi_i changed {pi0!r} -> {st.pi_i!r}"
                assert st.flags["vertex"], f"step {k}: the vertex rule did not fire"
                dev = st.eps_p - np.trace(st.eps_p) / 3.0 * I3
                assert np.linalg.norm(dev) <= 1e-14, f"step {k}: deviatoric plastic strain {np.linalg.norm(dev):.3e}"
            else:
                assert abs(st.pi_i - pi0) <= 1e-9 * abs(pi0), f"step {k}: pi_i changed {pi0!r} -> {st.pi_i!r}"
            assert st.eps_p_s <= 1e-14, f"step {k}: eps^p_s = {st.eps_p_s:.3e} (no deviatoric flow on the axis)"
            d_epv = st.eps_p_v - prev_epv
            assert d_epv < 0.0
            assert abs(st.D - pic * d_epv) <= 10.0 * tolp * abs(pic * d_epv), f"step {k}: D = {st.D:.8e}, pi_c Deps^p_v = {pic * d_epv:.8e}"
        elif elastic_cf:
            assert not is_plastic(st), f"step {k}: trial inside the surface but the step is plastic"
            assert abs(st.pi_i - pi0) <= 1e-9 * abs(pi0)
        prev_epv = st.eps_p_v
    assert n_plastic >= (N_VERTEX if case == "apex_start" else N_VERTEX - 1)


@pytest.mark.parametrize("cap", ["none", "planar"])
def test_i2_o1_stops_at_the_vertex_when_the_increment_is_not_isotropic(cap):
    """Sheet 3.2 (revised, 'What the two integrators do there'): on the axis the yield function has no deviatoric gradient
    (f = F_p delta/3) and Q has a vertex, so the rate problem has NO consistent solution for a non-isotropic strain rate: the
    deviatoric part of a^e:eps' is invisible to f yet moves the stress off the axis.  The rate oracle must therefore STOP and
    say so, not complete.  Documented statuses (o1_rate README): vertex_reached (the state arrives at q = 0 in finite time,
    sheet: measured on the near-isotropic path), vertex_nonisotropic (an increment that starts on the axis with a deviatoric
    part).  This test starts AT the apex (pi_c = -100, F = 0) and applies deps = -DELTA I/3 + 5e-5 diag(1, 0, -1): the stop must be
    reported at the first increment, with a status from that vocabulary and not 'ok'.  The pure isotropic increment of the same
    state is the positive test above (O1 rides the axis, README: 'A pure hydrostatic path is fine').
    MUTANT KILLED: an O1 that integrates through the vertex with some arbitrary deviatoric direction (the stop is the contract)."""
    kw = dict(PAPER, **VERTEX_CAPS[cap])
    P = make_params("O1", **kw)
    st0 = ORACLES["O1"].initial_state(P, -100.0 * I3, V0_K2, pi_i_for_pic(kw, -100.0))
    deps = np.array([-DELTA / 3.0 * I3 + 5e-5 * np.diag([1.0, 0.0, -1.0])])
    sts = run_deps("O1", P, st0, deps)
    status = sts[-1].flags["status"]
    assert status in ("vertex_reached", "vertex_nonisotropic"), f"O1 status {status!r}"
    assert len(sts) == 1


def test_i2_o1_vertex_reached_on_the_near_isotropic_path_status():
    """Sheet 3.2 (revised, 'What the two integrators do there'): on the near-isotropic G1 path (no cap, K2 set rho 0.7 / 0.8, paper
    CSL, deviator/volumetric ratio 0.2/3 -- AMP_STOP of test_g1_convergence_tangents, n = 40) the rate oracle must STOP with status
    vertex_reached BEFORE completing the path (fewer than 40 increments returned), at q = 0: the contract is that the state arrives
    at the vertex in finite time. q <= 1e-6 |p|: the event fires at R = 1e-7 |p|, i.e. q = 1.2e-7 |p| (sheet 3.2: R_tol, a specified
    constant, not a measurement). No number measured from O1 output is asserted here (see the regression pin below).
    MUTANT KILLED: an O1 without the vertex event (it would run on through q = 0 and complete or chatter)."""
    from test_g1_convergence_tangents import AMP_STOP, N_CAP, _run_cap
    P, sts = _run_cap("O1", "none", AMP_STOP)
    st = sts[-1]
    assert st.flags["status"] == "vertex_reached", st.flags["status"]
    assert len(sts) < N_CAP, (len(sts), N_CAP)
    p, q, _, _ = pq_theta(st.sigma)
    assert q <= 1e-6 * abs(p), q


def test_i2_o1_vertex_stop_location_regression_pin():
    """Regression pin of measured O1 output, not a gate. Sheet 3.2 prints the measured stop of the rate oracle on the same path
    (AMP_STOP, n = 40, no cap): increment 13 of 40, last state p = -236.7 kPa. These two numbers were MEASURED from O1, not derived
    from the sheet; they are pinned only so that a change in the integrator's event handling or step control shows up as a diff
    to be looked at. p within 0.05 kPa (rounding of the printed digits). A failure here is a prompt to re-measure and
    update the sheet, not by itself evidence of a defect; the vertex contract is gated by the status test above."""
    from test_g1_convergence_tangents import AMP_STOP, N_CAP, _run_cap
    P, sts = _run_cap("O1", "none", AMP_STOP)
    st = sts[-1]
    assert len(sts) == 13 and N_CAP == 40, (len(sts), N_CAP)
    p, _, _, _ = pq_theta(st.sigma)
    assert abs(p - (-236.7)) <= 0.05, p


# ==============================================================================================
# I.3  lower bound on the dissipation (sheet 11.1, 11.2)
# ==============================================================================================
# b_min = minimum over (theta, eta) of the bracket of S.38, (M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta, as printed in the table
# of sheet 11.2.  K2 set (M 1.2, N .4, N_bar .2, rho .7, rho_bar .8): +0.375.  rho = rho_bar = 1 (zeta_bar/zeta = 1): the bracket
# is linear in eta, smallest at eta = M/N: (M - M/N)/(1 - N_bar) + M/N = 0.750 = the table's third row (same M, N, N_bar; rho >= rho_bar
# keeps zeta_bar/zeta >= 1 with equality at theta = pi/3).
B_MIN = {"TXC_paper": 0.375, "TXC_fork": 0.375, "TXE_paper": 0.375, "UND_paper": 0.375, "UND_fork": 0.375,
         "OC_TXE_ok": 0.375, "TXC_2inv": 0.750}


@pytest.mark.parametrize("case", sorted(B_MIN))
def test_i3_dissipation_lower_bound(oracle_name, case):
    """S.38: D^p = -lambda' p [(M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta] >= b_min |p| lambda' on the yield surface, b_min from the
    table of sheet 11.2 (see B_MIN).  On a plastic increment with Delta eps^p_s > 0:
        D >= b_min min(|p_n|, |p_{n+1}|) Delta lambda,     Delta lambda = Delta eps^p_s / (sqrt(2/3) Omega)  (S.20).
    Omega is constant along an axisymmetric (corner) path (it depends on theta and rho_bar only), so Delta lambda is exact for both
    oracles from the corner Omega recomputed here (omega_corner).  O2: D = Delta lambda sigma_{n+1}:q_{n+1} exactly, so the
    bound holds at the converged state.  O1: D = int lambda' sigma:q dt >= b_min int |p| lambda' dt >= b_min min_t|p(t)| Delta lambda; the
    min over the increment is taken over the end points (p is monotone along these paths except across a q-peak, where the 30 %
    slack of the bound at its tightest path, OC_TXE_ok, covers the second-order difference).  The inequality is asserted with a
    round-off allowance of 1e-9 relative; D > 0 is asserted separately.  At least 8 plastic increments per path.
    MUTANTS KILLED: an oracle (O1 or O2) that reports D == 0 or D < 0 on plastic increments; a dissipation of the wrong sign or
    missing the volumetric term p beta F_p (the bracket would go below b_min at low eta); a dissipation that omits q zeta_bar."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    corner = CASES[case]["corner"]
    bmin = B_MIN[case]
    prev, n_chk, worst = st0, 0, np.inf
    for k, st in enumerate(sts):
        dEs = st.eps_p_s - prev.eps_p_s
        if is_plastic(st) and dEs > 1e-12:
            Om = omega_corner(kw, st.sigma, corner)
            dlam = dEs / (SQ23 * Om)
            pmin = min(abs(pq_theta(prev.sigma)[0]), abs(pq_theta(st.sigma)[0]))
            bound = bmin * pmin * dlam
            assert st.D > 0.0, f"{case} step {k}: D = {st.D:.3e} on a plastic increment"
            assert st.D >= bound * (1.0 - 1e-9), f"{case} step {k}: D = {st.D:.6e} < b_min |p| dlam = {bound:.6e} (ratio {st.D / bound:.4f})"
            worst = min(worst, st.D / bound)
            n_chk += 1
        prev = st
    assert n_chk >= 8, f"only {n_chk} plastic increments with Delta eps^p_s > 0"


# ==============================================================================================
# I.4  elastic-energy convexity table (plan 2.5 / 6.2 G1; sheet 2.2, 2.3, S.5)
# ==============================================================================================
# What the sheet claims (2.2): K = D11 = -p/kappa > 0 always; for alpha0 = 0 the Hessian determinant is -3 mu0 p0 e^omega/kappa > 0
# "always"; with alpha0 != 0 "it must be tabulated".  2.3: the full 6x6 Hessian additionally needs q/eps_s > 0.  Written here from
# S.5, with u := -p0 e^omega > 0 and a := 3 alpha0/(2 kappa):
#   D11 = (u/kappa)(1 + a eps_s^2) > 0,   D22 = 3(mu0 + alpha0 u),   D12 = -3 u alpha0 eps_s/kappa,   q/eps_s = 3 mu^e = 3(mu0 + alpha0 u),
#   det D = (3u/kappa) [ mu0 (1 + a eps_s^2) + alpha0 u (1 - a eps_s^2) ].                                   (E.1)
# 6x6: a^e = D11 1x1 + sqrt(2/3) D12 (1xn + nx1) + (2/3) D22 nxn + (2q/(3 eps_s)) (IDEV - nxn) has the eigenvalues
# 2 mu^e (x4) and those of [[3 D11, sqrt2 D12],[sqrt2 D12, (2/3) D22]] (det = 2 det D)  -- PD iff D11 > 0, det D > 0, mu^e > 0.  (E.2)
# Consequences: for a eps_s^2 < 1 both terms of (E.1) are positive for every alpha0 >= 0, i.e. PD everywhere for eps_s < sqrt(2 kappa/(3 alpha0))
# (0.183 for the K2 kappa, 0.316 for BA06 Table 1) -- the whole footing range eps_s <= 0.08.  Loss needs BOTH a eps_s^2 > 1 and
# alpha0 u > mu0, i.e. u > mu0/alpha0 (eps_v < -kappa ln(mu0/(alpha0 |p0|)) = -0.0560 for K2 at alpha0 = 0.2).
E_SETS = {"K2": dict(kappa_hat=0.01, mu0=5400.0), "BA06T1": dict(kappa_hat=0.03, mu0=2000.0)}   # sheet 14 / 14 K2b (BA06 Table 1)
E_ALPHAS = [0.0, 0.2]
EV_GRID = np.round(np.arange(-0.06, 0.02 + 1e-12, 0.002), 6)          # 41 points, the task's eps_v range
ES_FOOT = np.round(np.arange(0.0, 0.08 + 1e-12, 0.005), 6)            # 17 points, the task's eps_s range
ES_EXT = np.round(np.arange(0.0, 0.80 + 1e-12, 0.02), 6)              # 41 points, beyond the footing range (location of loss)
NHAT = np.array([1.0, 0.0, -1.0]) / SQ2          # deviatoric unit direction, distinct principal values (theta = pi/6): spin terms active
MANDEL_IDX = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]
MANDEL_W = [1.0, 1.0, 1.0, SQ2, SQ2, SQ2]


def e_params(setname, alpha0):
    return dict(PAPER, p0=-100.0, eps_v0=0.0, alpha0=alpha0, **E_SETS[setname])


def S5(kw, ev, es):
    """(S.5): p, q, mu^e, D11, D12, D22."""
    k, p0, a0, mu0 = kw["kappa_hat"], kw["p0"], kw["alpha0"], kw["mu0"]
    E = math.exp(-(ev - kw["eps_v0"]) / k)
    p = p0 * E * (1.0 + 1.5 * a0 * es * es / k)
    mu = mu0 - a0 * p0 * E
    return p, 3.0 * mu * es, mu, -p / k, 3.0 * p0 * a0 * es * E / k, 3.0 * mu0 - 3.0 * a0 * p0 * E


def eig_closed(D11, D12, D22, mu):
    """Eigenvalues of the 6x6 Mandel tangent (E.2) and of the 2x2 Hessian (S.5), ascending."""
    e2 = np.linalg.eigvalsh(np.array([[3.0 * D11, SQ2 * D12], [SQ2 * D12, 2.0 / 3.0 * D22]]))
    return np.sort(np.concatenate([[2.0 * mu] * 4, e2])), np.linalg.eigvalsh(np.array([[D11, D12], [D12, D22]]))


def mandel6(C):
    M = np.empty((6, 6))
    for a, (i, j) in enumerate(MANDEL_IDX):
        for b, (k, l) in enumerate(MANDEL_IDX):
            M[a, b] = MANDEL_W[a] * MANDEL_W[b] * C[i, j, k, l]
    return M


def oracle_point(oname, P, ev, es):
    """(el, C4): the oracle's energy-level quantities at eps^e = diag(eps_v/3 + sqrt(3/2) eps_s nhat)  (|e| = sqrt(3/2) eps_s)."""
    eps = ev / 3.0 * np.ones(3) + SQ32 * es * NHAT
    if oname == "O1":
        from o1_rate.model import energy
        el = energy(np.diag(eps), P)
        return el, el.a4
    from o2_algo import kernel as K
    el = K.elastic(P, eps)
    return el, K.tangent_small(el.ae, el.sig, eps, np.eye(3))


@functools.lru_cache(maxsize=None)
def grid_eval(oname, setname, alpha0, extended):
    kw = e_params(setname, alpha0)
    P = make_params(oname, **kw)
    es_grid = ES_EXT if extended else ES_FOOT
    rows = []
    for ev in EV_GRID:
        for es in es_grid:
            p, q, mu, D11, D12, D22 = S5(kw, float(ev), float(es))
            e6c, e2c = eig_closed(D11, D12, D22, mu)
            el, C4 = oracle_point(oname, P, float(ev), float(es))
            M6 = mandel6(C4)
            sym = float(np.max(np.abs(M6 - M6.T)) / np.max(np.abs(M6)))
            e6 = np.linalg.eigvalsh(0.5 * (M6 + M6.T))
            e2 = np.linalg.eigvalsh(np.array([[el.D11, el.D12], [el.D12, el.D22]]))
            rows.append(dict(ev=float(ev), es=float(es), p=p, q=q, mu=mu, D=(D11, D12, D22), e6c=e6c, e2c=e2c,
                             o_p=float(el.p), o_q=float(el.q), o_D=(float(el.D11), float(el.D12), float(el.D22)),
                             e6=e6, e2=e2, sym=sym, es_o=float(el.eps_s)))
    return rows


def predicted_pd(rows):
    """+1 = PD predicted by the closed form (min eigenvalue > 1e-8 lambda_max), -1 = not PD predicted, 0 = ambiguous band (skipped)."""
    out = []
    for r in rows:
        lmax = float(r["e6c"][-1])
        m6, m2 = float(r["e6c"][0]), float(r["e2c"][0])
        if m6 > 1e-8 * lmax and m2 > 1e-8 * float(r["e2c"][-1]):
            out.append(1)
        elif m6 < -1e-8 * lmax or m2 < -1e-8 * float(r["e2c"][-1]):
            out.append(-1)
        else:
            out.append(0)
    return np.array(out)


@pytest.mark.parametrize("alpha0", E_ALPHAS)
@pytest.mark.parametrize("setname", sorted(E_SETS))
def test_i4_sheet_closed_forms_det_D_and_eigen_structure(setname, alpha0):
    """Checks of the sheet itself (no oracle): (a) alpha0 = 0: det D = D11 D22 - D12^2 = -3 mu0 p0 e^omega/kappa (sheet 2.2) to
    1e-12 relative, and > 0 on the whole grid; (b) (E.1): det D = (3u/kappa)[mu0 (1 + a eps_s^2) + alpha0 u (1 - a eps_s^2)] from S.5, 1e-11
    relative; (c) the 6x6 eigenvalues (E.2) against a direct 6x6 Mandel assembly of the closed-form a^e = (S.3) built here with
    numpy (an algebra check of the derivation used below).  Footing grid and extended grid."""
    kw = e_params(setname, alpha0)
    k, p0, a0, mu0 = kw["kappa_hat"], kw["p0"], alpha0, kw["mu0"]
    a = 3.0 * a0 / (2.0 * k)
    for es_grid in (ES_FOOT, ES_EXT):
        for ev in EV_GRID:
            for es in es_grid:
                p, q, mu, D11, D12, D22 = S5(kw, float(ev), float(es))
                det = D11 * D22 - D12 ** 2
                u = -p0 * math.exp(-float(ev) / k)
                e1 = (3.0 * u / k) * (mu0 * (1.0 + a * es ** 2) + a0 * u * (1.0 - a * es ** 2))
                assert abs(det - e1) <= 1e-11 * abs(e1) + 1e-12 * abs(D11 * D22), (ev, es, det, e1)
                if alpha0 == 0.0:
                    e0 = -3.0 * mu0 * p0 * math.exp(-float(ev) / k) / k
                    assert abs(det - e0) <= 1e-12 * e0, (ev, es, det, e0)
                    assert det > 0.0 and D11 > 0.0 and mu > 0.0
    # (c) direct assembly of S.3 at one generic point
    ev, es = -0.03, 0.05
    p, q, mu, D11, D12, D22 = S5(kw, ev, es)
    n = np.diag(NHAT)
    IxI = np.einsum("ij,kl->ijkl", I3, I3)
    ISYM = 0.5 * (np.einsum("ik,jl->ijkl", I3, I3) + np.einsum("il,jk->ijkl", I3, I3))
    nxn = np.einsum("ij,kl->ijkl", n, n)
    C = (D11 * IxI + SQ23 * D12 * (np.einsum("ij,kl->ijkl", I3, n) + np.einsum("ij,kl->ijkl", n, I3)) + 2.0 / 3.0 * D22 * nxn
         + (2.0 * q / (3.0 * es)) * (ISYM - IxI / 3.0 - nxn))
    e_dir = np.linalg.eigvalsh(mandel6(C))
    e_cf = eig_closed(D11, D12, D22, mu)[0]
    assert np.max(np.abs(e_dir - e_cf)) <= 1e-10 * e_cf[-1], (e_dir, e_cf)


@pytest.mark.parametrize("alpha0", E_ALPHAS)
@pytest.mark.parametrize("setname", sorted(E_SETS))
def test_i4_oracle_hessian_and_6x6_tangent_match_the_sheet(oracle_name, setname, alpha0):
    """The oracle's elastic quantities at eps^e = diag(eps_v/3 + sqrt(3/2) eps_s nhat) against the sheet closed forms (S.5), (E.2):
    p, q, D11, D12, D22, eps_s to 1e-12 relative (to the largest of the quantities at the point); the 6x6 Mandel tangent symmetric
    to 1e-12 (it is a Hessian of Psi, sheet 2.1) and its six eigenvalues equal to (E.2) to 1e-10 lambda_max (round-off of the 6x6
    assembly, O2 builds it spectrally with g_ab = (sigma_a - sigma_b)/(eps_a - eps_b)); the oracle's 2x2 Hessian eigenvalues to
    1e-11.  Footing grid, eps_v in [-0.06, 0.02] x eps_s in [0, 0.08] (697 points).
    MUTANTS KILLED: dropping D12 or the (2/3)(D22 - q/eps_s) nxn term of S.3 (a 5e-4 relative shift of the lowest eigenvalue at
    eps_s = 0.08); K = D11 not proportional to p; the shear modulus not mu0 - alpha0 p0 e^omega; a non-symmetric 6x6."""
    rows = grid_eval(oracle_name, setname, alpha0, False)
    assert len(rows) == len(EV_GRID) * len(ES_FOOT)
    for r in rows:
        scale = max(abs(r["p"]), abs(r["q"]), r["D"][0], r["D"][2])
        tag = (oracle_name, setname, alpha0, r["ev"], r["es"])
        assert abs(r["es_o"] - r["es"]) <= 1e-12, tag
        assert abs(r["o_p"] - r["p"]) <= 1e-12 * abs(r["p"]), tag
        assert abs(r["o_q"] - r["q"]) <= 1e-12 * max(abs(r["q"]), 1e-3 * abs(r["p"])), tag
        for got, want in zip(r["o_D"], r["D"]):
            assert abs(got - want) <= 1e-12 * scale, (tag, got, want)
        assert r["sym"] <= 1e-12, (tag, r["sym"])
        assert np.max(np.abs(r["e6"] - r["e6c"])) <= 1e-10 * r["e6c"][-1], (tag, r["e6"], r["e6c"])
        assert np.max(np.abs(r["e2"] - r["e2c"])) <= 1e-11 * r["e2c"][-1], (tag, r["e2"], r["e2c"])


@pytest.mark.parametrize("alpha0", E_ALPHAS)
@pytest.mark.parametrize("setname", sorted(E_SETS))
def test_i4_convexity_on_the_footing_grid(oracle_name, setname, alpha0):
    """Where the sheet says positive definiteness must hold: sheet 2.2: alpha0 = 0 -> det D > 0 always, K > 0, mu^e = mu0 > 0, so the
    2x2 Hessian AND the 6x6 tangent are PD at every point of eps_v in [-0.06, 0.02] x eps_s in [0, 0.08] (also the 'q/eps_s > 0' of 2.3).
    alpha0 = 0.2: (E.1)/(E.2) are tabulated (the sheet's instruction) and give a eps_s^2 <= 0.19 (K2) / 0.064 (BA06 T1) < 1, hence both terms
    of det D positive and D11, mu^e > 0: PD at every grid point too -- asserted as a statement about the closed form AND about the
    oracle (smallest eigenvalue of the oracle's 6x6 and 2x2 > 0 at every point, in the direction the closed form predicts).  The
    sign agreement oracle <-> closed form is asserted point by point in both directions (generic 'assert the location' form).
    MUTANT KILLED: an energy with the wrong-sign shear coupling, mu^e = mu0 + alpha0 p0 e^omega (negative for u > mu0/alpha0, i.e. at
    eps_v <= -0.056 with the K2 kappa, inside this grid): the oracle's 6x6 and 2x2 would lose definiteness there."""
    rows = grid_eval(oracle_name, setname, alpha0, False)
    pred = predicted_pd(rows)
    assert np.all(pred == 1), f"closed form predicts loss of convexity on the footing grid: {[(r['ev'], r['es']) for r, c in zip(rows, pred) if c != 1]}"
    for r, c in zip(rows, pred):
        assert r["e6"][0] > 0.0 and r["e2"][0] > 0.0, (oracle_name, setname, alpha0, r["ev"], r["es"], r["e6"][0], r["e2"][0])


def test_i4_loss_of_convexity_location_beyond_the_footing_range(oracle_name):
    """The branch 'if the sheet predicts loss of convexity somewhere, assert its location instead', exercised where it exists.
    K2 kappa = 0.01, mu0 = 5400, alpha0 = 0.2, eps_s extended to 0.8: (E.1) is negative iff mu0 (1 + a eps_s^2) + alpha0 u (1 - a eps_s^2) < 0, which
    needs a eps_s^2 > 1 (eps_s > 0.1826) AND alpha0 u > mu0 (eps_v < -0.0560): a sliver at the strongest compression and large shear.
    Asserted: (1) the closed form predicts a non-empty not-PD set on this grid; (2) every predicted not-PD point lies in
    {eps_v < -0.0559, eps_s > 0.1826} (the necessary conditions above, with 1e-3 allowance for the rounding of those two numbers);
    (3) the oracle's 6x6 and 2x2 smallest eigenvalues are negative exactly at the predicted not-PD points and positive at the predicted
    PD points (points within 1e-8 lambda_max of zero are skipped); (4) for alpha0 = 0 the same grid is PD everywhere (sheet 2.2).
    MUTANT KILLED: any energy whose Hessian stays PD where S.5 says it is lost (e.g. a symmetrised or 'fixed' D12), and any that loses it elsewhere."""
    rows = grid_eval(oracle_name, "K2", 0.2, True)
    pred = predicted_pd(rows)
    nonpd = [r for r, c in zip(rows, pred) if c == -1]
    assert nonpd, "the closed form predicts no loss of convexity on the extended grid (expected a sliver at eps_v ~ -0.06, eps_s > 0.18)"
    k, mu0, a0, p0 = 0.01, 5400.0, 0.2, -100.0
    es_c = math.sqrt(2.0 * k / (3.0 * a0))
    ev_c = -k * math.log(mu0 / (a0 * abs(p0)))
    for r in nonpd:
        assert r["es"] > es_c * (1.0 - 1e-3) and r["ev"] < ev_c * (1.0 - 1e-3) + 1e-12, (r["ev"], r["es"], es_c, ev_c)
    for r, c in zip(rows, pred):
        if c == 1:
            assert r["e6"][0] > 0.0 and r["e2"][0] > 0.0, (r["ev"], r["es"], r["e6"][0], r["e2"][0])
        elif c == -1:
            assert r["e6"][0] < 0.0 and r["e2"][0] < 0.0, (r["ev"], r["es"], r["e6"][0], r["e2"][0])
    rows0 = grid_eval(oracle_name, "K2", 0.0, True)
    assert all(r["e6"][0] > 0.0 and r["e2"][0] > 0.0 for r in rows0)


@pytest.mark.parametrize("alpha0", E_ALPHAS)
def test_i4_public_api_tangent_matches_the_closed_form(oracle_name, alpha0):
    """The same eigenvalue check through the PUBLIC route (initial_state -> tangent), BA06 Table 1 set, 9 states eps_v in {-0.04, -0.01,
    0.01} x eps_s in {0.01, 0.04, 0.08}.  The state is built from the stress (p, q, nhat) of the closed form S.5, so the oracle's own inverse
    of the energy must return (eps_v, eps_s) (1e-9), and the oracle's 4th-order tangent, converted to the 6x6 Mandel matrix, must have
    the eigenvalues of (E.2) to 1e-9 lambda_max (inverse tolerance 1e-13 relative, amplified by <~ 1e3)."""
    kw = e_params("BA06T1", alpha0)
    P = make_params(oracle_name, **kw)
    ora = ORACLES[oracle_name]
    for ev in (-0.04, -0.01, 0.01):
        for es in (0.01, 0.04, 0.08):
            p, q, mu, D11, D12, D22 = S5(kw, ev, es)
            sig = np.diag(p * np.ones(3) + SQ23 * q * NHAT)
            st = ora.initial_state(P, sig, 1.7, -1.0e4)
            ee = np.linalg.eigvalsh(st.eps_e)
            ev_o = float(ee.sum())
            es_o = SQ23 * float(np.linalg.norm(ee - ev_o / 3.0))
            assert abs(ev_o - ev) <= 1e-9 and abs(es_o - es) <= 1e-9, (oracle_name, alpha0, ev, es, ev_o, es_o)
            M6 = mandel6(ora.tangent(P, st))
            e6 = np.linalg.eigvalsh(0.5 * (M6 + M6.T))
            e6c = eig_closed(D11, D12, D22, mu)[0]
            assert np.max(np.abs(e6 - e6c)) <= 1e-9 * e6c[-1], (oracle_name, alpha0, ev, es, e6, e6c)


def _bin_table(rows_by_oracle, key):
    """min over the oracles and over the points in each (|p|, eta = q/|p|) bin of the smallest eigenvalue (6x6 or 2x2)."""
    pb = [10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, float("inf")]
    eb = [0.0, 0.5, 1.0, 1.5, 2.5, 5.0, float("inf")]
    cell = [[None] * (len(eb) - 1) for _ in range(len(pb) - 1)]
    for rows in rows_by_oracle:
        for r in rows:
            ap, eta = abs(r["p"]), r["q"] / abs(r["p"])
            if ap < pb[0]:
                continue
            i = next(j for j in range(len(pb) - 1) if pb[j] <= ap < pb[j + 1])
            j = next(j for j in range(len(eb) - 1) if eb[j] <= eta < eb[j + 1])
            m = float(r[key][0])
            cell[i][j] = m if cell[i][j] is None else min(cell[i][j], m)
    hdr = "| \\|p\\| (kPa) \\ q/\\|p\\| | " + " | ".join(f"[{eb[j]:g}, {eb[j + 1]:g})" for j in range(len(eb) - 1)) + " |"
    sep = "|---|" + "---|" * (len(eb) - 1)
    lines = [hdr, sep]
    for i in range(len(pb) - 1):
        lines.append(f"| [{pb[i]:g}, {pb[i + 1]:g}) | " + " | ".join("-" if c is None else f"{c:.4g}" for c in cell[i]) + " |")
    return "\n".join(lines)


def test_i4_write_energy_convexity_table():
    """Writes tests/out/energy_convexity.md: for each parameter set and alpha0, the smallest eigenvalue of the oracle 6x6 tangent and of
    the 2x2 invariant Hessian binned over (|p|, q/|p|) on the footing grid, the p and q/|p| ranges the grid covers, and the closed-form
    prediction of the not-PD set on the extended grid.  Asserted: the file is written and (both oracles, all sets) the footing-grid
    minima are positive for alpha0 in {0, 0.2}, the table's own claim."""
    os.makedirs(OUT_DIR, exist_ok=True)
    out = ["# Elastic-energy convexity table (WP-144 G1, test_g1_identities.py::test_i4_*)", "",
           "BA06 energy (S.4)-(S.5), alpha0 in {0, 0.2}; grid eps_v in [-0.06, 0.02] (step 0.002) x eps_s in [0, 0.08] (step 0.005), "
           "elastic strain eps^e = diag(eps_v/3 + sqrt(3/2) eps_s nhat), nhat = (1, 0, -1)/sqrt2.  p0 = -100 kPa, eps_v0 = 0.",
           "The task's 'p ~ -10 ... -1000 kPa' cannot be hit by one kappa over this eps_v range (p = p0 exp(-eps_v/kappa)); the two kappa of the "
           "sheet are tabulated and their actual p ranges stated.  q/|p| is the stress ratio sqrt(3/2)|s|/|p| (not the yield-surface eta).",
           "Eigenvalues in kPa.  6x6 = Mandel form of the 4th-order tangent (equivalent to Voigt for definiteness).  Minimum over O1 and O2.", ""]
    for setname in sorted(E_SETS):
        for alpha0 in E_ALPHAS:
            both = [grid_eval(o, setname, alpha0, False) for o in ("O1", "O2")]
            rows = both[0]
            ps = [abs(r["p"]) for r in rows]
            etas = [r["q"] / abs(r["p"]) for r in rows]
            m6 = min(float(r["e6"][0]) for rs in both for r in rs)
            m2 = min(float(r["e2"][0]) for rs in both for r in rs)
            assert m6 > 0.0 and m2 > 0.0, (setname, alpha0, m6, m2)
            out += [f"## {setname} (kappa = {E_SETS[setname]['kappa_hat']}, mu0 = {E_SETS[setname]['mu0']:g}), alpha0 = {alpha0}", "",
                    f"|p| from {min(ps):.1f} to {max(ps):.1f} kPa; q/|p| from {min(etas):.3f} to {max(etas):.2f}; points {len(rows)}; "
                    f"smallest eigenvalue over the grid: 6x6 {m6:.5g}, 2x2 {m2:.5g}; closed-form not-PD points: "
                    f"{int(np.sum(predicted_pd(rows) == -1))}.", "",
                    "Smallest eigenvalue of the 6x6 tangent:", "", _bin_table(both, "e6"), "",
                    "Smallest eigenvalue of the 2x2 Hessian D:", "", _bin_table(both, "e2"), ""]
    rows = grid_eval("O2", "K2", 0.2, True)
    pred = predicted_pd(rows)
    bad = [(r["ev"], r["es"]) for r, c in zip(rows, pred) if c == -1]
    out += ["## Loss of convexity (K2, alpha0 = 0.2), beyond the footing range, eps_s up to 0.8", "",
            "(E.1): det D = (3u/kappa)[mu0 (1 + a eps_s^2) + alpha0 u (1 - a eps_s^2)], u = -p0 e^omega, a = 3 alpha0/(2 kappa); negative only for "
            "a eps_s^2 > 1 (eps_s > 0.1826) and alpha0 u > mu0 (eps_v < -0.0560).  Predicted not-PD grid points "
            f"({len(bad)}): " + (", ".join(f"({ev:g}, {es:g})" for ev, es in bad) if bad else "none") + "", "",
            "alpha0 = 0 is PD everywhere (det D = -3 mu0 p0 e^omega/kappa > 0, sheet 2.2).", ""]
    with open(os.path.join(OUT_DIR, "energy_convexity.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    assert os.path.getsize(os.path.join(OUT_DIR, "energy_convexity.md")) > 500
