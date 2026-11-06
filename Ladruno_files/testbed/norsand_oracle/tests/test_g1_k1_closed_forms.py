"""WP-144 gate G1 (Zone B): K1 closed-form identities, K1.1 - K1.10 (plan 5.2, sheet 13), both oracles.

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form in the
equation sheet 144a (the section is cited at each use), (b) a number printed in the sheet, or (c) a stated
convergence / quadrature argument.  All of it was written BEFORE the oracles were run on these paths.
No oracle output is used as an expected value.  A tolerance is never loosened to make a test pass; a failing
assertion is a finding (against the sheet or an oracle), not a reason to edit an oracle or a bound.

Gate ids (plan 5.2 K1.n):
  K1.1  isotropic hyperelastic compression p(eps_v)                         sheet 2.2, S.5, 13.1
  K1.2  closed non-coaxial elastic loop: zero net work, state returns       sheet 2.2, 13.2
  K1.3  zeta corners, zeta' = 0 at the corners, convexity ranges, refusals  sheet 4, S.10-S.11, 13.3
  K1.4  eta = M at p = pi_i (image point), q = M|pi_i|/zeta(theta)          sheet 5.1, 13.4
  K1.5  flow-rule identity at plastic steps with q > 0                      sheet 5.3, S.20, 13.5
  K1.6  peak identity: D = chi psi_i where H changes sign                   sheet 7, 8, 13.6
  K1.7  undrained critical-state endpoint (paper and fork CSL)              sheet 6, 13.7
  K1.8  drained critical-state asymptote                                    sheet 13.8
  K1.9  dissipation D >= 0; refusal rule (condition A, S.39); forced counterexample   sheet 11, 13.9
  K1.10 specific volume v = v0 exp(tr eps) on a path (exponential update, G2 decision)    sheet 1.2, 13.10

How the two oracles are read (they differ only in names / diagnostics; conftest hides the rest):
  * stress p, q, theta are always recomputed HERE from State.sigma (principal values), never taken from an oracle;
  * pi_i, v, eps_e, D, eps_p_v, eps_p_s are State fields common to both;
  * pi_i* is State.flags['pistar'] (O1) or State.pi_star (O2) -- used only to LOCATE H = 0, never as an expected
    value: what is asserted at the located state is D = chi psi_i, with D and psi_i recomputed here from the sheet.
  * O1 is a continuum integrator: a step quantity such as eps_p_v/eps_p_s over an increment is a weighted mean of the
    point value (weights eps_p_s'), so for O1 it is bounded, not equated (argument at each use).  O2 is backward Euler:
    its step quantities ARE the point values at the converged state (sheet 13.5), so they are equated.
"""
import functools
import math
import warnings

import numpy as np
import pytest
import sympy as sp

from conftest import K2_BASE, ORACLES, make_params

SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)

# ----------------------------------------------------------------------------------------------
# parameter sets (sheet 14 / 15).  Condition A (S.39) holds for all of them: N_bar <= N and rho/rho_bar >= beta.
# ----------------------------------------------------------------------------------------------
PAPER = dict(K2_BASE, rho=0.7, rho_bar=0.8)              # AB06 6.1 case 2; 0.875 >= beta = 0.75
PAPER_2INV = dict(K2_BASE, rho=1.0, rho_bar=1.0)         # 2-invariant: zeta = 1, Omega = sqrt(3/2) exactly (S.20)
FORK = dict(PAPER, csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45, p_a=101.325)   # task: fork CSL values


# ----------------------------------------------------------------------------------------------
# sheet closed forms, written here from the sheet (not from the oracles)
# ----------------------------------------------------------------------------------------------
def v0_for(kw, pi0, psi0):
    """Specific volume giving image state parameter psi0 at pi_i = pi0 (sheet 6, S.22)."""
    if kw["csl_mode"] == "paper":
        return psi0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-pi0)
    return 1.0 + kw["e0"] + psi0 - kw["lambda_c"] * (-pi0 / kw["p_a"]) ** kw["xi"]


def psi_of(kw, v, pi):
    """psi_i = v - v_c(pi_i)  (sheet 6, S.22)."""
    if kw["csl_mode"] == "paper":
        return v - kw["v_c0"] + kw["lambda_tilde"] * math.log(-pi)
    return (v - 1.0) - kw["e0"] + kw["lambda_c"] * (-pi / kw["p_a"]) ** kw["xi"]


def eta_S12(kw, p, pi):
    """eta(p, pi_i) of S.12."""
    M, N = kw["M"], kw["N"]
    if N == 0.0:
        return M * (1.0 + math.log(pi / p))
    return (M / N) * (1.0 - (1.0 - N) * (p / pi) ** (N / (1.0 - N)))


def ratio_from_eta(kw, eta):
    """x = p/pi_i on the surface for a given eta (inverse of S.12, BA06 2.8 as in sheet 5.1)."""
    M, N = kw["M"], kw["N"]
    if N == 0.0:
        return 1.0 / math.exp(eta / M - 1.0)
    return (((1.0 - eta * N / M) / (1.0 - N)) ** ((1.0 - N) / N))


def pq_theta(sig):
    """p, q, theta, principal values from a stress tensor: S.1 (theta = 0 tension corner, pi/3 compression corner)."""
    w = np.linalg.eigvalsh(0.5 * (sig + sig.T))
    p = float(w.sum()) / 3.0
    xi = w - p
    R = float(np.linalg.norm(xi))
    if R <= 0.0:
        return p, 0.0, float("nan"), w
    y = float((xi ** 3).sum()) / R ** 3
    return p, SQ32 * R, math.acos(max(-1.0, min(1.0, SQ6 * y))) / 3.0, w


@functools.lru_cache(maxsize=None)
def _ww_sym():
    """Willam-Warnke (S.11) and its theta-derivatives, exact via sympy, lambdified.
    returns (zeta, zeta', zeta'', convexity functional r^2 + 2 r'^2 - r r'' with r = 1/zeta)."""
    th, rho = sp.symbols("theta rho", positive=True)
    c = sp.cos(th)
    A = 4 * (1 - rho ** 2)
    B = 2 * rho - 1
    z = (A * c ** 2 + B ** 2) / (2 * (1 - rho ** 2) * c + B * sp.sqrt(A * c ** 2 + 5 * rho ** 2 - 4 * rho))
    z1 = sp.diff(z, th)
    z2 = sp.diff(z1, th)
    r = 1 / z
    r1 = sp.diff(r, th)
    r2 = sp.diff(r1, th)
    conv = r ** 2 + 2 * r1 ** 2 - r * r2
    lam = lambda e: sp.lambdify((th, rho), e, "numpy")          # noqa: E731
    return lam(z), lam(z1), lam(z2), lam(conv)


@functools.lru_cache(maxsize=None)
def _ga_sym():
    """Gudehus-Argyris (S.10), same quantities."""
    th, rho = sp.symbols("theta rho", positive=True)
    z = ((1 + rho) + (1 - rho) * sp.cos(3 * th)) / (2 * rho)
    z1 = sp.diff(z, th)
    z2 = sp.diff(z1, th)
    r = 1 / z
    r1 = sp.diff(r, th)
    r2 = sp.diff(r1, th)
    conv = r ** 2 + 2 * r1 ** 2 - r * r2
    lam = lambda e: sp.lambdify((th, rho), e, "numpy")          # noqa: E731
    return lam(z), lam(z1), lam(z2), lam(conv)


def zeta_sheet(kind, theta, rho):
    return np.asarray((_ww_sym() if kind == "WW" else _ga_sym())[0](theta, rho), dtype=float)


def zeta_corner(rho, corner):
    """(zeta, zeta_y) of S.9 at a corner: 'C' = theta pi/3 (compression), 'E' = theta 0 (extension), WW.
    zeta_y(0) = -sqrt6 zeta''(0)/9, zeta_y(pi/3) = +sqrt6 zeta''(pi/3)/9, with zeta'' exact from S.11."""
    z2 = _ww_sym()[2]
    if corner == "E":
        return 1.0 / rho, -SQ6 * float(z2(0.0, rho)) / 9.0
    return 1.0, SQ6 * float(z2(math.pi / 3.0, rho)) / 9.0


def omega_corner(kw, sig, corner):
    """Omega of S.20 at an axisymmetric (corner) stress: Omega^2 = (3/2) zeta_b^2 + (zeta_b_y q)^2 sum_a y_a^2,
    y_a from S.6 with delta_a = 1."""
    p, q, _, w = pq_theta(sig)
    zb, zby = zeta_corner(kw["rho_bar"], corner)
    xi = w - p
    R = float(np.linalg.norm(xi))
    ya = 3.0 * xi ** 2 / R ** 3 - 3.0 * float((xi ** 3).sum()) * xi / R ** 5 - 1.0 / R
    return math.sqrt(1.5 * zb ** 2 + (zby * q) ** 2 * float((ya ** 2).sum()))


def dilatancy_closed(kw, sig, pi, corner):
    """D = sqrt(3/2) beta F_p / Omega (sheet 5.3 / AB06 39) with beta F_p = (eta - M)/(1 - N_bar), eta from S.12."""
    p, _, _, _ = pq_theta(sig)
    eta = eta_S12(kw, p, pi)
    return SQ32 * (eta - kw["M"]) / ((1.0 - kw["N_bar"]) * omega_corner(kw, sig, corner))


def dissipation_per_lambda(kw, sig, pi, corner):
    """sigma:q = p beta F_p + q zeta_bar  (sheet 11.1, before the substitution q = -p eta/zeta), so D_step = dlam * this."""
    p, q, _, _ = pq_theta(sig)
    eta = eta_S12(kw, p, pi)
    return p * (eta - kw["M"]) / (1.0 - kw["N_bar"]) + q * zeta_corner(kw["rho_bar"], corner)[0]


# ----------------------------------------------------------------------------------------------
# oracle-agnostic runners
# ----------------------------------------------------------------------------------------------
def init_state(oname, kw, sig0, v0, pi0, unchecked=False):
    P = make_params(oname, unchecked=unchecked, **kw)
    return P, ORACLES[oname].initial_state(P, sig0, v0, pi0)


def run_deps(oname, P, st0, deps):
    ora = ORACLES[oname]
    deps = np.asarray(deps, float)
    return ora.run_path(P, st0, deps, rtol=1e-10) if oname == "O1" else ora.run_path(P, st0, deps)


def run_tx(oname, P, st0, kind, total, n):
    ora = ORACLES[oname]
    if oname == "O1":
        return ora.triaxial(P, st0, kind, total, n, rtol=1e-10)
    return ora.triaxial(P, st0, kind, total, n)


def ok(oname, st):
    f = st.flags
    return (f.get("status", "ok") in ("ok", "initial")) if oname == "O1" else (not f.get("refused", False))


def why(oname, st):
    f = st.flags
    return f.get("status") if oname == "O1" else f.get("reason")


def is_plastic(st):
    return bool(st.flags.get("plastic", False))


def pistar_of(oname, st):
    return st.flags["pistar"] if oname == "O1" else st.pi_star


def run_staged(oname, P, st0, kind, sign, sched):
    """Concatenate triaxial stages; sched = [(|axial increment of stage|, n_steps)].  Returns (states, cumulative signed
    axial strain after each state).  Stops at the first non-ok state (which is included)."""
    sts, cum, st, c = [], [], st0, 0.0
    for tot, n in sched:
        seg = run_tx(oname, P, st, kind, sign * tot, n)
        for s in seg:
            c += sign * tot / n
            sts.append(s)
            cum.append(c)
        st = seg[-1]
        if not ok(oname, st):
            break
    return sts, np.array(cum)


# ----------------------------------------------------------------------------------------------
# path library (initial state, stage schedule); cached per (oracle, case) and shared by K1.5 / 1.6 / 1.8 / 1.9
# ----------------------------------------------------------------------------------------------
SIG0 = -100.0 * I3
TX_SCHED = [(0.04, 20), (0.46, 46), (2.5, 50)]            # step 2e-3, 1e-2, 5e-2; total |eps_ax| = 3.0
U_SCHED = [(0.04, 20), (0.96, 48), (4.0, 80)]             # step 2e-3, 2e-2, 5e-2; total |eps_ax| = 5.0
P_CS_TARGET = -250.0                                      # undrained CS pressure the initial v0 is built to give

CASES = {
    "TXC_paper": dict(kw=PAPER, kind="drained", sign=-1, sched=TX_SCHED, pi0=-80.0, psi0=-0.05, corner="C"),
    "TXC_fork": dict(kw=FORK, kind="drained", sign=-1, sched=TX_SCHED, pi0=-80.0, psi0=-0.05, corner="C"),
    "TXE_paper": dict(kw=PAPER, kind="drained", sign=+1, sched=TX_SCHED, pi0=-80.0, psi0=-0.05, corner="E"),
    "TXC_2inv": dict(kw=PAPER_2INV, kind="drained", sign=-1, sched=TX_SCHED[:2], pi0=-80.0, psi0=-0.05, corner="C"),
    "UND_paper": dict(kw=PAPER, kind="undrained", sign=-1, sched=U_SCHED, pi0=-60.0, p_cs=P_CS_TARGET, corner="C"),
    "UND_fork": dict(kw=FORK, kind="undrained", sign=-1, sched=U_SCHED, pi0=-60.0, p_cs=P_CS_TARGET, corner="C"),
    # heavily over-consolidated extension: pi_i >> |p|, so yield is reached at the extension corner with |p|/|pi_i| ~ 0.06
    # (sheet 11 "state reaches theta ~ 0 near the apex").  _ok: condition A holds (K2 set); _bad: N_bar = N, unchecked.
    "OC_TXE_ok": dict(kw=PAPER, kind="drained", sign=+1, sched=[(0.03, 60)], pi0=-1000.0, psi0=-0.03, corner="E"),
    "OC_TXE_bad": dict(kw=dict(PAPER, N_bar=0.4), kind="drained", sign=+1, sched=[(0.03, 60)], pi0=-1000.0,
                       psi0=-0.03, corner="E", unchecked=True),
}


def case_v0(spec):
    kw = spec["kw"]
    if "p_cs" in spec:     # invert the closed-form undrained endpoint (sheet 13.7) for v0
        pcs = spec["p_cs"]
        if kw["csl_mode"] == "paper":
            return kw["v_c0"] - kw["lambda_tilde"] * math.log(-pcs)
        return 1.0 + kw["e0"] - kw["lambda_c"] * (-pcs / kw["p_a"]) ** kw["xi"]
    return v0_for(kw, spec["pi0"], spec["psi0"])


@functools.lru_cache(maxsize=None)
def run_case(oname, case):
    spec = CASES[case]
    P, st0 = init_state(oname, spec["kw"], SIG0, case_v0(spec), spec["pi0"], spec.get("unchecked", False))
    sts, cum = run_staged(oname, P, st0, spec["kind"], spec["sign"], spec["sched"])
    return P, spec["kw"], st0, sts, cum


def assert_path_ok(oname, sts, cum, case):
    bad = [k for k, s in enumerate(sts) if not ok(oname, s)]
    assert not bad, f"{oname}/{case}: path stopped at step {bad[0]} ({why(oname, sts[bad[0]])})"


# ==============================================================================================
# K1.1  isotropic hyperelastic compression (sheet 2.2, S.5 at eps_s = 0; 13.1)
# ==============================================================================================
K11_SETS = [
    # (overrides, sigma0 pressure, pi_i0).  pi_i0 is chosen so the whole path stays elastic (asserted below).
    (dict(), -100.0, -200.0),
    (dict(eps_v0=0.002, alpha0=0.2, kappa_hat=0.02), -120.0, -300.0),
]


@pytest.mark.parametrize("over,p_init,pi0", K11_SETS)
def test_k1_1_isotropic_compression_matches_closed_form(oracle_name, over, p_init, pi0):
    """p(eps_v) = p0 exp(-(eps_v - eps_v0)/kappa) for eps_s = 0, for any alpha0 (S.5).  Elastic path, so tolerance is
    round-off: the map is an exponential of a linear function; 1e-10 relative is the O1 ODE tolerance (rtol 1e-10)."""
    kw = dict(PAPER_2INV, **over)
    p0, kap, ev0 = kw["p0"], kw["kappa_hat"], kw["eps_v0"]
    _, st0 = init_state(oracle_name, kw, p_init * I3, 1.7, pi0)
    P = make_params(oracle_name, **kw)
    # the reference state of the energy: tr(eps^e) at p_init follows from inverting S.5 at eps_s = 0
    ev_init = ev0 - kap * math.log(p_init / p0)
    assert abs(np.trace(st0.eps_e) - ev_init) <= 1e-12, (np.trace(st0.eps_e), ev_init)
    dv = 0.001
    steps = [-dv] * 10 + [+dv] * 20                        # eps_v: 0 -> -0.01 -> +0.01
    deps = np.array([np.eye(3) * s / 3.0 for s in steps])
    sts = run_deps(oracle_name, P, st0, deps)
    cum = np.cumsum(steps)
    for k, (s, c) in enumerate(zip(sts, cum)):
        assert ok(oracle_name, s), (k, why(oracle_name, s))
        assert not is_plastic(s), f"step {k} yielded: not an elastic path (raise |pi_i0|)"
        p_closed = p0 * math.exp(-((ev_init + c) - ev0) / kap)
        p, q, _, _ = pq_theta(s.sigma)
        assert abs(p - p_closed) <= 1e-10 * abs(p_closed), (k, p, p_closed)
        assert q <= 1e-9 * abs(p), (k, q)                      # isotropic strain -> no deviatoric stress (S.5)
    if over == {}:
        # the number printed in sheet 13.1: eps_v = -0.01 -> p = -100 e = -271.828 kPa
        assert abs(pq_theta(sts[9].sigma)[0] - (-100.0 * math.e)) <= 1e-10 * 271.828


# ==============================================================================================
# K1.2  closed non-coaxial elastic loop (sheet 2.2, 13.2)
# ==============================================================================================
def _sym(a):
    return np.array(a, float)


# vertices of a closed loop in total-strain space; every segment changes the principal axes (shear terms on all pairs)
LOOP = [
    np.zeros((3, 3)),
    _sym([[0.0020, 0.0015, 0.0], [0.0015, -0.0010, 0.0010], [0.0, 0.0010, 0.0005]]),
    _sym([[-0.0010, -0.0020, 0.0015], [-0.0020, 0.0015, 0.0], [0.0015, 0.0, -0.0015]]),
    _sym([[0.0005, 0.0, -0.0020], [0.0, 0.0020, 0.0015], [-0.0020, 0.0015, -0.0005]]),
    np.zeros((3, 3)),
]


@pytest.mark.parametrize("alpha0", [0.0, 0.2])
def test_k1_2_closed_nonconaxial_loop_zero_work_and_state_returns(oracle_name, alpha0):
    """W = oint sigma:d(eps) = 0 because sigma = dPsi/d(eps^e) (BA06 2.1; sheet 2.2), to 1e-12 of oint |sigma:d(eps)|.

    Quadrature argument: along each segment the oracle stress at eps(t) = E_j + t Delta_j is obtained by ONE oracle step
    of size t Delta_j from the segment start (elastic, hence path independent); sigma:Delta_j is analytic in t
    (exp of a linear function, range ~0.4 in the exponent), so 12-point Gauss-Legendre is exact to ~1e-15.
    State return: sigma, eps^e, pi_i equal their initial values to 1e-12 (relative to |p0|; strain absolute)."""
    kw = dict(PAPER_2INV, alpha0=alpha0)
    P, st0 = init_state(oracle_name, kw, SIG0, 1.7, -250.0)
    xg, wg = np.polynomial.legendre.leggauss(12)
    tg, wg = 0.5 * (xg + 1.0), 0.5 * wg
    W = Wabs = 0.0
    st = st0
    for j in range(len(LOOP) - 1):
        D = LOOP[j + 1] - LOOP[j]
        for t, w in zip(tg, wg):
            s = run_deps(oracle_name, P, st, [t * D])[-1]
            assert ok(oracle_name, s) and not is_plastic(s), "loop left the elastic domain"
            d = float(np.sum(s.sigma * D))
            W += w * d
            Wabs += w * abs(d)
        st = run_deps(oracle_name, P, st, [D])[-1]
        assert ok(oracle_name, st) and not is_plastic(st)
    assert Wabs > 0.0
    assert abs(W) <= 1e-12 * Wabs, f"net work {W:.3e}, sum|.| {Wabs:.3e}, ratio {abs(W) / Wabs:.2e}"
    assert np.linalg.norm(st.sigma - st0.sigma) <= 1e-12 * abs(kw["p0"]), np.linalg.norm(st.sigma - st0.sigma)
    assert np.linalg.norm(st.eps_e - st0.eps_e) <= 1e-12, np.linalg.norm(st.eps_e - st0.eps_e)
    assert st.pi_i == st0.pi_i


# ==============================================================================================
# K1.3  zeta corners, zeta' = 0 at the corners, convexity ranges, refusals (sheet 4, 13.3)
# ==============================================================================================
WW_RHOS = [0.55, 0.7, 0.85, 1.0]
GA_RHOS = [7.0 / 9.0, 0.85, 1.0]
ZETA_CASES = [("WW", r) for r in WW_RHOS] + [("GA", r) for r in GA_RHOS]


def oracle_zeta(oname, theta, rho, kind):
    """(zeta, zeta_y) from the oracle's own shape function (the corner branch S.9 is used exactly at the corners)."""
    if oname == "O1":
        from o1_rate.model import zeta_fun
        z, zy = zeta_fun(theta, math.cos(3.0 * theta) / SQ6, rho, kind)
        return float(z), float(zy)
    from o2_algo.kernel import zeta_y
    z, zy, _ = zeta_y(theta, rho, kind)
    return float(z), float(zy)


@pytest.mark.sheet_check
@pytest.mark.parametrize("kind,rho", ZETA_CASES)
def test_k1_3_zeta_corner_values_and_zero_slope(kind, rho):
    """[sheet_check: never calls an oracle; not oracle coverage.]
    Sheet 4: zeta(0) = 1/rho, zeta(pi/3) = 1, zeta'(0) = zeta'(pi/3) = 0 -- on the sheet's own closed form
    (this is the statement the oracles are then held to below).  Tolerance 1e-12: exact in floating point up to round-off."""
    z, z1, _, _ = (_ww_sym() if kind == "WW" else _ga_sym())
    assert abs(float(z(0.0, rho)) - 1.0 / rho) <= 1e-12
    assert abs(float(z(math.pi / 3.0, rho)) - 1.0) <= 1e-12
    assert abs(float(z1(0.0, rho))) <= 1e-12
    assert abs(float(z1(math.pi / 3.0, rho))) <= 1e-12


@pytest.mark.parametrize("kind,rho", ZETA_CASES)
def test_k1_3_oracle_zeta_corners_interior_and_zeta_y(oracle_name, kind, rho):
    """Oracle shape function vs the sheet closed forms.
    corners: zeta(0) = 1/rho, zeta(pi/3) = 1 (1e-13).
    interior: zeta(theta) against S.10 / S.11 implemented here (1e-12).
    zeta_y = d zeta / d y, y = cos3theta/sqrt6 (S.1), away from the corners (sheet 3.1 rule ii: not within 1e-4 of a corner):
        central difference of the sheet's zeta in y with d theta = 1e-5; truncation O(d theta^2) ~ 1e-10, round-off ~ 1e-11;
        tolerance 1e-7 relative.
    at the corners: the limits S.9, zeta_y(0) = -sqrt6 zeta''(0)/9 and zeta_y(pi/3) = +sqrt6 zeta''(pi/3)/9, with zeta''
        exact from S.11 by sympy (GA: zeta_y = sqrt6 (1-rho)/(2 rho), S.10); tolerance 1e-9 relative (complex-step or
        closed-form constants in the oracles)."""
    z0, zy0 = oracle_zeta(oracle_name, 0.0, rho, kind)
    zc, zyc = oracle_zeta(oracle_name, math.pi / 3.0, rho, kind)
    assert abs(z0 - 1.0 / rho) <= 1e-13 and abs(zc - 1.0) <= 1e-13
    if kind == "GA":
        lim0 = limc = SQ6 * (1.0 - rho) / (2.0 * rho)
    else:
        z2 = _ww_sym()[2]
        lim0 = -SQ6 * float(z2(0.0, rho)) / 9.0
        limc = SQ6 * float(z2(math.pi / 3.0, rho)) / 9.0
    for got, lim in ((zy0, lim0), (zyc, limc)):
        assert abs(got - lim) <= 1e-9 * max(1.0, abs(lim)), (got, lim)
    y = lambda th: math.cos(3.0 * th) / SQ6                              # noqa: E731
    for th in (0.05, 0.3, 0.5, 0.8, 1.0):
        z, zy = oracle_zeta(oracle_name, th, rho, kind)
        assert abs(z - float(zeta_sheet(kind, th, rho))) <= 1e-12, th
        d = 1e-5
        fd = (float(zeta_sheet(kind, th + d, rho)) - float(zeta_sheet(kind, th - d, rho))) / (y(th + d) - y(th - d))
        assert abs(zy - fd) <= 1e-7 * max(1.0, abs(fd)), (th, zy, fd)


@pytest.mark.sheet_check
def test_k1_3_sheet_corner_second_derivative_table():
    """[sheet_check: never calls an oracle; not oracle coverage.]
    The table of sheet 4.2 (zeta''(0), zeta''(pi/3), zeta_y(0), zeta_y(pi/3)) against zeta'' exact from S.11 (sympy):
    a check of the sheet itself.  Tolerance 1e-9 relative (the table carries 10-11 digits)."""
    table = {   # rho: (zeta''(0), zeta''(pi/3), zeta_y(0), zeta_y(pi/3))   as printed in sheet 4.2
        0.7: (-1.1208791209, 9.5625, 0.3050646566, 2.6025828517),
        0.71: (-1.0828693089, 8.4336734694, 0.2947196961, 2.2953551841),
        0.8: (-0.75, 3.0, 0.2041241452, 0.8164965809),
        1.0: (0.0, 0.0, 0.0, 0.0),
    }
    z2 = _ww_sym()[2]
    for rho, (a0, a3, y0, y3) in table.items():
        assert abs(float(z2(0.0, rho)) - a0) <= 1e-9 * max(1.0, abs(a0)), rho
        assert abs(float(z2(math.pi / 3.0, rho)) - a3) <= 1e-9 * max(1.0, abs(a3)), rho
        assert abs(zeta_corner(rho, "E")[1] - y0) <= 1e-9, rho
        assert abs(zeta_corner(rho, "C")[1] - y3) <= 1e-9, rho


@pytest.mark.sheet_check
def test_k1_3_convexity_ranges_ww_half_to_one_ga_seven_ninths_to_one():
    """[sheet_check: never calls an oracle; not oracle coverage.]
    Sheet 4: the deviatoric section r = 1/zeta is convex iff r^2 + 2 r'^2 - r r'' >= 0 on [0, pi/3].
    WW: convex for rho in [1/2, 1] (min = -8.9e-15 at 1/2 printed in the sheet -> tolerance -1e-9), not at 0.45.
    rho = 1/2 is the CONVEXITY boundary only (the measure is identically 0 there, zeta = 2cos(theta)); it is not
    admissible: the parser refuses it (owner decision 2026-10-01, sheet 4.2), see the refusal tests below.
    GA: convex for [7/9, 1], not at 0.77.  4001-point grid, exact derivatives (sympy)."""
    th = np.linspace(0.0, math.pi / 3.0, 4001)
    ww, ga = _ww_sym()[3], _ga_sym()[3]
    for rho in (0.5, 0.55, 0.7, 0.85, 1.0):
        assert float(np.min(ww(th, rho))) >= -1e-9, rho
    assert float(np.min(ww(th, 0.45))) <= -1e-3
    for rho in (7.0 / 9.0, 0.85, 1.0):
        assert float(np.min(ga(th, rho))) >= -1e-9, rho
    assert float(np.min(ga(th, 0.77))) <= -1e-4


def build(oname, **kw):
    """Params as the parser would see them: O1 validates in the constructor, O2 in validate()."""
    P = make_params(oname, **kw)
    P.validate()
    return P


@pytest.mark.parametrize("label,kw", [
    ("WW rho below 1/2", dict(PAPER, rho=0.45, rho_bar=0.45)),
    ("WW rho_bar below 1/2", dict(PAPER, rho=0.7, rho_bar=0.45)),
    ("GA rho below 7/9", dict(PAPER, zeta="GA", rho=0.77, rho_bar=0.77)),
    ("rho above 1", dict(PAPER, rho=1.05, rho_bar=1.05)),
])
def test_k1_3_refuses_rho_outside_convex_range(oracle_name, label, kw):
    """K1.3: the parser refuses rho outside (1/2, 1] (WW) / [7/9, 1] (GA).  rho = rho_bar in every case so that
    condition A (S.39) holds and the refusal can only be the convexity range."""
    with pytest.raises(ValueError):
        build(oracle_name, **kw)


@pytest.mark.parametrize("label,kw", [
    ("GA rho = 7/9 (boundary)", dict(PAPER, zeta="GA", rho=7.0 / 9.0, rho_bar=7.0 / 9.0)),
    # WW just above 1/2 is admissible (sheet 4.2: "a value like 0.51 is admissible but a poor choice"); with the
    # refusal tests below it pins the boundary at exactly 1/2 and not at some larger value.
    ("WW rho = rho_bar = 0.51", dict(PAPER, rho=0.51, rho_bar=0.51)),
    ("WW rho = 0.51, rho_bar = 0.6", dict(PAPER, rho=0.51, rho_bar=0.6)),
    ("WW rho = 0.7, rho_bar = 0.51", dict(PAPER, rho=0.7, rho_bar=0.51)),
])
def test_k1_3_accepts_the_admissible_ends_of_the_ranges(oracle_name, label, kw):
    """GA 7/9 stays accepted (closed range [7/9, 1], sheet 4.1); WW accepts every rho in (1/2, 1], here 0.51 for rho and
    for rho_bar.  Condition A (S.39) holds for all three WW sets (0.51/0.6 = 0.85, 0.7/0.51, 1 >= beta = 0.75)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # rho > rho_bar is a warning, not a refusal (sheet 11.2)
        P = build(oracle_name, **kw)
    assert P.rho == kw["rho"] and P.rho_bar == kw["rho_bar"]


@pytest.mark.parametrize("label,kw", [
    ("WW rho = rho_bar = 1/2", dict(PAPER, rho=0.5, rho_bar=0.5)),
    # isolate each ellipticity: condition A (S.39) and the convexity range of the OTHER one are satisfied, so the only
    # possible reason for the refusal is the one named in the label (0.5/0.6 = 0.833 >= 0.75; 0.7/0.5 = 1.4 >= 0.75)
    ("WW rho = 1/2 alone (rho_bar = 0.6)", dict(PAPER, rho=0.5, rho_bar=0.6)),
    ("WW rho_bar = 1/2 alone (rho = 0.7)", dict(PAPER, rho=0.7, rho_bar=0.5)),
])
def test_k1_3_refuses_ww_rho_one_half_exactly(oracle_name, label, kw):
    """Owner decision 2026-10-01 (sheet 4.2, 13.3, 16.1.6): WW rho = 1/2 EXACTLY is refused, for rho and for rho_bar,
    because there B = 2 rho - 1 = 0, zeta = 2 cos(theta), zeta'(pi/3) = -sqrt3 != 0 and the compression corner is a vertex
    of the deviatoric section (so zeta o sigma is only C^0 there).  The same set with 0.5 -> 0.51 is accepted (test above),
    so the refusal is the range and not condition A.
    Mutants killed: the pre-decision parser (accept 1/2, refuse only rho < 1/2); a range check applied to rho but not to
    rho_bar (or the reverse); a range check applied only to the first of the two."""
    with pytest.raises(ValueError):
        build(oracle_name, **kw)
    ok_kw = dict(kw, rho=kw["rho"] + (0.01 if kw["rho"] == 0.5 else 0.0),
                 rho_bar=kw["rho_bar"] + (0.01 if kw["rho_bar"] == 0.5 else 0.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        build(oracle_name, **ok_kw)                    # control: 0.51 is accepted (not condition A)


# ==============================================================================================
# K1.4  image point: eta = M at p = pi_i (sheet 5.1, 13.4)
# ==============================================================================================
def sigma_from_pq_theta(p, q, theta):
    """Principal stresses sigma_k = p + (2/3) q cos(theta - 2 pi k/3), k = 0,1,2: sum xi = 0, |xi| = sqrt(2/3) q, and
    sum xi^3 / R^3 = cos3theta / sqrt6, i.e. the theta of S.1 (theta = pi/3: two equal, one more compressive than them)."""
    k = np.arange(3)
    return np.diag(p + (2.0 / 3.0) * q * np.cos(theta - 2.0 * math.pi * k / 3.0))


IMAGE_THETAS = [0.0, math.pi / 12.0, math.pi / 6.0, math.pi / 4.0, math.pi / 3.0]


@pytest.mark.parametrize("N,N_bar", [(0.4, 0.2), (0.0, 0.0)])
@pytest.mark.parametrize("kind,rho", [("WW", 0.6), ("WW", 0.7), ("WW", 1.0), ("GA", 0.8), ("GA", 1.0)])
def test_k1_4_image_point_eta_equals_M_at_p_equals_pi_i(oracle_name, kind, rho, N, N_bar):
    """F = zeta q + p eta = 0 with p = pi_i and eta(pi_i, pi_i) = M  <=>  q = M|pi_i| / zeta(theta): q = M|pi_i| in TXC,
    q = rho M |pi_i| in TXE, and the sheet's zeta(theta) in between (S.10 / S.11 written here).
    The oracle is given a stress ON that point and asked to place pi_i on the surface through it (initial_state with
    pi_i0 = None, the inverse of S.12): it must return pi_i = p.  The two branches N > 0 and N = 0 of S.12 are both hit.
    Tolerance 1e-9 relative: inversion of S.12 amplifies a zeta error by d ln(pi/p)/d eta <= 1/(M(1-N)) ~ 1."""
    kw = dict(PAPER, zeta=kind, rho=rho, rho_bar=rho, N=N, N_bar=N_bar)
    M, p = kw["M"], -100.0
    for th in IMAGE_THETAS:
        q = M * abs(p) / float(zeta_sheet(kind, th, rho))
        sig = sigma_from_pq_theta(p, q, th)
        tt = pq_theta(sig)[2]
        assert abs(tt - th) <= (1e-7 if th in (0.0, math.pi / 3.0) else 1e-9), (tt, th)          # construction sanity
        P, st = init_state(oracle_name, kw, sig, 1.7, None)
        assert abs(st.pi_i / p - 1.0) <= 1e-9, (th, st.pi_i)
        if oracle_name == "O1":
            # O1 reports F/(M|p|) of its own yield function at the state: pi_i0 = p on this stress is F = 0
            _, st2 = init_state(oracle_name, kw, sig, 1.7, p)
            assert abs(st2.flags["F_rel"]) <= 1e-9, (th, st2.flags["F_rel"])
    # corner values stated explicitly (the gate text): q = M|pi_i| (TXC), q = rho M |pi_i| (TXE)
    for th, qexp in ((math.pi / 3.0, M * abs(p)), (0.0, rho * M * abs(p))):
        assert abs(M * abs(p) / float(zeta_sheet(kind, th, rho)) - qexp) <= 1e-12 * qexp


# ==============================================================================================
# K1.5  flow-rule identity at plastic steps with q > 0 (sheet 5.3, S.20, 13.5)
# ==============================================================================================
@pytest.mark.parametrize("case", ["TXC_2inv", "TXC_paper", "TXC_fork", "TXE_paper"])
def test_k1_5_flow_rule_identity_at_plastic_steps(oracle_name, case):
    """Delta eps^p_v / Delta eps^p_s = sqrt(3/2) beta F_p / Omega  (AB06 39, sheet 5.3) at the converged state.
    Right side recomputed here: F_p from S.12 at (p, pi_i), Omega from S.20 with y_a of S.6 at the state's stress, with the
    corner values of zeta_bar, zeta_bar_y (S.9; sympy zeta'').  The paths are axisymmetric so every state IS on a corner.
    2-invariant case (rho = rho_bar = 1): Omega = sqrt(3/2) exactly, so D = (eta - M)/(1 - N_bar) (BA06 2.25).
    O2: backward Euler, Delta eps^p = Delta lambda q(sigma_{n+1}, pi_{i,n+1}) exactly, so the step ratio IS the point value at
        the converged state: equal to 1e-8 relative (local Newton tolerance 1e-12 scaled, amplification <~ 1e3), every plastic step.
    O1: a continuum integrator reports increments, and the ratio over an increment is the eps^p_s'-weighted MEAN of D(t), not its
        point value (D is not monotone across the drained peak, so the end-point values do not bracket it).  The point value
        is therefore measured by PROBING from a state with two tiny increments d = 2e-6 and d/2 (axial strain) and a first-order
        Richardson extrapolation, point = 2 r(d/2) - r(d), which cancels the O(d) chord error (remaining O(d^2 D'') ~ 1e-10).
        Noise: the ODE tolerance is rtol 1e-10 on eps^p (~5e-2), i.e. <~ 5e-12 absolute on an increment of ~2e-6, relative <~ 3e-6,
        so r is good to ~1e-6 and the extrapolation to ~5e-6; tolerance 5e-5 absolute (factor 10), every 4th plastic state."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    spec = CASES[case]
    corner = spec["corner"]
    prev, n_checked, n_pl = st0, 0, 0
    for k, st in enumerate(sts):
        if is_plastic(st) and pq_theta(st.sigma)[1] > 0.0:
            d_end = dilatancy_closed(kw, st.sigma, st.pi_i, corner)
            if oracle_name == "O2":
                ds = st.eps_p_s - prev.eps_p_s
                if ds > 1e-12:
                    ratio = (st.eps_p_v - prev.eps_p_v) / ds
                    assert abs(ratio - d_end) <= 1e-8 * max(1.0, abs(d_end)), (k, ratio, d_end)
                    n_checked += 1
            else:
                n_pl += 1
                if n_pl % 4 == 1:
                    def chord(d):
                        s2 = run_tx("O1", P, st, spec["kind"], spec["sign"] * d, 1)[-1]
                        assert ok("O1", s2) and is_plastic(s2), (k, d, why("O1", s2))
                        return (s2.eps_p_v - st.eps_p_v) / (s2.eps_p_s - st.eps_p_s)
                    point = 2.0 * chord(1e-6) - chord(2e-6)
                    assert abs(point - d_end) <= 5e-5, (k, point, d_end)
                    n_checked += 1
        prev = st
    assert n_checked >= 10, f"only {n_checked} plastic states checked: the path is not a flow-rule test"


# ==============================================================================================
# K1.6  peak identity: D = chi psi_i where H changes sign (sheet 7, 8, 13.6)
# ==============================================================================================
def _sign_changes(oname, sts):
    """Indices k with both states k-1, k plastic and a = pi_i* - pi_i changing sign.
    H = -M (p/pi_i)^{1/(1-N)} sqrt(2/3) h (pi_i* - pi_i) Omega (sheet 8) with every factor but (pi_i* - pi_i) positive, so
    sgn H = -sgn a: H > 0 (hardening) iff a < 0."""
    out = []
    for k in range(1, len(sts)):
        if is_plastic(sts[k]) and is_plastic(sts[k - 1]):
            a0 = pistar_of(oname, sts[k - 1]) - sts[k - 1].pi_i
            a1 = pistar_of(oname, sts[k]) - sts[k].pi_i
            if a0 * a1 < 0.0:
                out.append(k)
    return out


@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper"])
def test_k1_6_peak_identity_D_equals_chi_psi_at_H_zero(oracle_name, case):
    """At any state with pi_i = pi_i*: H = 0 and D = chi psi_i (S.23, sheet 13.6).  It is a STATE identity, so it does not
    depend on the path or step size that produced the state.  But O2's steps are finite, so the state with H = 0 is not one
    of the step end points: it is located by bracketing the sign change of a = pi_i* - pi_i between states k-1 and k of the
    coarse path and bisecting with single backward-Euler (O2) / Radau (O1) sub-steps of size f*Delta from state k-1 (the same
    constitutive call, so the found state is a genuine state of the oracle), until |a|/|pi_i| <= 1e-10.
    At that state D is recomputed here (closed form of 5.3, S.20) and compared with chi psi_i(v, pi_i) (S.22).
    Tolerance argument: D(eta) - chi psi_i vanishes identically at pi_i = pi_i* and is Lipschitz in a/|pi_i| with constant
    L <~ 10 (dD/d eta = sqrt(3/2)/((1-N_bar) Omega) ~ 1, d eta / d ln pi_i <= M/(1-N) ~ 2, chi dpsi/dln pi_i = chi Lambda ~ 0.05);
    with |a|/|pi_i| <= 1e-10 the residual is ~1e-9, and 1e-6 leaves a factor 1e3 for the oracle solver tolerance (1e-10).
    Sign: dense drained, so H > 0 before the crossing (a < 0) and H < 0 after (peak then softening), sheet 8.
    Secondary, coarse-path check (no re-solving): linear interpolation of b = D - chi psi_i to the zero of a between the two
    bracketing states.  b = c a + d a^2 + ... with b(a=0) = 0, so the interpolant is off by |d a_{k-1} a_k| (quadratic in the
    step) while the bracket variation |b_k - b_{k-1}| ~ |c (a_k - a_{k-1})| is linear: require |b_interp| <= 0.1 |Delta b| + 1e-9."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    spec = CASES[case]
    corner = spec["corner"]
    ks = _sign_changes(oracle_name, sts)
    assert ks, "H never changes sign on this path: not a peak-identity test (dense drained path must soften)"
    k = ks[0]
    lo_st, hi_st = sts[k - 1], sts[k]
    a_lo = pistar_of(oracle_name, lo_st) - lo_st.pi_i
    assert a_lo < 0.0, "expected hardening (H > 0, a < 0) before the crossing on a dense drained path"

    def bfun(st):
        return dilatancy_closed(kw, st.sigma, st.pi_i, corner) - kw["chi"] * psi_of(kw, st.v, st.pi_i)

    # secondary: linear interpolation on the coarse path
    a0, a1 = a_lo, pistar_of(oracle_name, hi_st) - hi_st.pi_i
    b0, b1 = bfun(lo_st), bfun(hi_st)
    s = a0 / (a0 - a1)
    b_interp = b0 + s * (b1 - b0)
    assert abs(b_interp) <= 0.1 * abs(b1 - b0) + 1e-9, (b0, b1, b_interp)

    # primary: bisection on the sub-step fraction
    dax = cum[k] - (cum[k - 1])
    lo, hi, found = 0.0, 1.0, None
    for _ in range(80):
        f = 0.5 * (lo + hi)
        st = run_tx(oracle_name, P, lo_st, spec["kind"], f * dax, 1)[-1]
        assert ok(oracle_name, st) and is_plastic(st), (f, why(oracle_name, st))
        a = pistar_of(oracle_name, st) - st.pi_i
        found = st
        if abs(a) <= 1e-10 * abs(st.pi_i):
            break
        if a * a_lo > 0.0:
            lo = f
        else:
            hi = f
    a = pistar_of(oracle_name, found) - found.pi_i
    assert abs(a) <= 1e-10 * abs(found.pi_i), f"bisection did not converge: a/|pi_i| = {abs(a / found.pi_i):.2e}"
    resid = bfun(found)
    assert abs(resid) <= 1e-6, (resid, dilatancy_closed(kw, found.sigma, found.pi_i, corner),
                                kw["chi"] * psi_of(kw, found.v, found.pi_i))
    assert pq_theta(found.sigma)[1] > 0.0


# ==============================================================================================
# K1.7  undrained critical-state endpoint (sheet 6, 13.7)
# ==============================================================================================
def p_cs_closed(kw, v):
    """Isochoric => v constant; CS <=> psi_i = 0 with pi_i = p (sheet 13.7)."""
    if kw["csl_mode"] == "paper":
        return -math.exp((kw["v_c0"] - v) / kw["lambda_tilde"])
    e = v - 1.0
    return -kw["p_a"] * ((kw["e0"] - e) / kw["lambda_c"]) ** (1.0 / kw["xi"])


@pytest.mark.parametrize("case", ["UND_paper", "UND_fork"])
def test_k1_7_undrained_triaxial_ends_on_the_csl(oracle_name, case):
    """An isochoric triaxial compression ends at p_cs = -exp((v_c0 - v)/lambda_tilde) (paper) or
    -p_a((e0 - e)/lambda_c)^(1/xi) (fork), with q_cs = M|p_cs| (theta = pi/3, zeta = 1), pi_i = p, psi_i = 0, D = 0.
    v is constant (sheet 1.2 / 13.7: v = v0 exp(tr eps), tr eps = 0, so v = v0 exactly, unchanged by the
    G2 exponential-update decision): asserted to 1e-12.
    Convergence argument: the CS is a fixed point (D = 0 => plastic volumetric strain rate 0 => p stationary) and an
    attractor; the approach is exponential in plastic shear strain with rate ~ |chi| per unit eps^p_s (psi_i relaxes as D = chi psi_i),
    so at |eps_ax| = 5 (eps^p_s ~ 5) the residual is ~ e^(-|chi| 5) ~ 1e-8; require p within 1e-6, the derived quantities
    within 1e-5 (they carry the same residual times O(1..10)), and the error at |eps_ax| = 2, 3.5, 5 non-increasing."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    v0 = case_v0(CASES[case])
    pcs = p_cs_closed(kw, v0)
    assert abs(pcs - P_CS_TARGET) <= 1e-9 * abs(P_CS_TARGET), "v0 inversion of the closed form is inconsistent"
    assert abs(cum[-1]) >= 5.0 - 1e-9, f"path ended at eps_ax = {cum[-1]}"
    assert max(abs(s.v - v0) for s in sts) <= 1e-12 * v0, "specific volume changed on an isochoric path"

    def p_err(at):
        k = int(np.argmin(np.abs(np.abs(cum) - at)))
        return abs(pq_theta(sts[k].sigma)[0] - pcs) / abs(pcs)

    e2, e35, e5 = p_err(2.0), p_err(3.5), p_err(5.0)
    assert e35 <= e2 + 1e-12 and e5 <= e35 + 1e-12, (e2, e35, e5)
    assert e5 <= 1e-6, f"p_cs relative error at eps_ax = 5: {e5:.3e}"
    end = sts[-1]
    p, q, th, _ = pq_theta(end.sigma)
    assert abs(th - math.pi / 3.0) <= 1e-7
    assert abs(q - kw["M"] * abs(pcs)) <= 1e-5 * kw["M"] * abs(pcs), (q, kw["M"] * abs(pcs))
    assert abs(end.pi_i / p - 1.0) <= 1e-5, end.pi_i / p
    assert abs(psi_of(kw, end.v, end.pi_i)) <= 1e-5
    assert abs(dilatancy_closed(kw, end.sigma, end.pi_i, "C")) <= 1e-5


# ==============================================================================================
# K1.8  drained critical-state asymptote (sheet 13.8)
# ==============================================================================================
@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper"])
def test_k1_8_drained_critical_state_asymptote(oracle_name, case):
    """At large shear strain: psi_i -> 0, pi_i -> p, eta -> M, i.e. -zeta(theta) q/p -> M, q/|p| -> M/zeta(theta)
    (M in TXC, rho M in TXE), D -> 0, H -> 0.  Closed forms only; approach argument: the plastic shear strain at
    |eps_ax| = 3 is ~3 and psi_i relaxes at the rate |chi| ~ 3.5 per unit eps^p_s (D = chi psi_i), i.e. e^(-10); tolerances
    2e-3 on eta/M - 1, 1e-3 on |psi_i|, 5e-3 on |pi_i/p - 1| and |D|, and the eta error at the end must be at most half
    of that at |eps_ax| = 1.5 (a stalled approach fails).  A miss is a finding (slower than the rate argument), not a tolerance to tune."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    assert abs(cum[-1]) >= 3.0 - 1e-9, f"path ended at eps_ax = {cum[-1]}"
    corner = CASES[case]["corner"]
    zc = 1.0 if corner == "C" else 1.0 / kw["rho"]
    M = kw["M"]

    def eta_err(st):
        p, q, _, _ = pq_theta(st.sigma)
        return abs(-zc * q / p / M - 1.0)

    kmid = int(np.argmin(np.abs(np.abs(cum) - 1.5)))
    end = sts[-1]
    p, q, th, _ = pq_theta(end.sigma)
    assert abs(th - (math.pi / 3.0 if corner == "C" else 0.0)) <= 1e-7
    assert eta_err(end) <= 2e-3, (eta_err(end), q / abs(p), M / zc)
    assert eta_err(end) <= 0.5 * eta_err(sts[kmid]), (eta_err(end), eta_err(sts[kmid]))
    assert abs(psi_of(kw, end.v, end.pi_i)) <= 1e-3
    assert abs(end.pi_i / p - 1.0) <= 5e-3
    assert abs(dilatancy_closed(kw, end.sigma, end.pi_i, corner)) <= 5e-3
    assert abs(pistar_of(oracle_name, end) / p - 1.0) <= 5e-3    # pi_i* -> p as psi_i -> 0 (S.23)


# ==============================================================================================
# K1.10  specific volume on a path: v = v0 exp(tr eps)  (sheet 1.2, 13.10; G2 owner decision 2026-10-01)
# ==============================================================================================
V_IDENTITY_TOL = {"O2": 1.0e-12, "O1": 1.0e-8}
V_IDENTITY_POWER = 1.0e-6


@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper", "TXC_2inv", "UND_paper", "UND_fork"])
def test_k1_10_specific_volume_is_v0_exp_of_the_total_strain_trace(oracle_name, case):
    """Sheet 13.10: at every committed state v = v0 exp(tr eps), eps = eps^e + eps^p the TOTAL strain (additivity,
    sheet 1.2), for any increment sequence and any substepping (exp of a sum).  Here tr eps = tr(State.eps_e) +
    State.eps_p_v is read from the state's own strain fields, v0 is the initial specific volume the case was built
    with (sheet 6 / 13.7 inversion), and the right-hand side uses no oracle v: the check ties the v-update to the
    independently integrated strains.  Drained paths (|tr eps| up to O(0.1 .. 1)) and undrained paths (tr eps = 0,
    v = v0: the control) are all covered.
    Tolerances.  O2: 1e-12 relative to v0 (sheet 13.10: round-off; exp of a sum, <= 148 steps of 1e-16).  O1: 1e-8,
    because tr eps^e + eps^p_v come from an ODE integrated at rtol = 1e-10 per increment (<= 148 increments, |eps^p_v|
    = O(0.1): accumulated <= 148 x 1e-10 x 0.5 ~ 7e-9); O1's v itself is algebraic in the total strain.
    The superseded linear update v0 (1 + x) would miss by v0 (1 + x - e^x) ~ -v0 x^2 / 2: asserted >= 1e-6 relative on
    every drained path, so the gate has power (sheet 13.10: any path with |tr eps| >~ 1e-3 .. 1e-6).
    Kills: v += v0 tr d_eps (linear), v_n (1 + tr d_eps), v += v tr d_eps; a v that is not carried across a substep or
    a plastic step; v reset to v0 at a branch switch."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    spec = CASES[case]
    v0 = case_v0(spec)
    assert abs(np.trace(st0.eps_e) + st0.eps_p_v) <= 1e-12, "premise: the start is the energy reference, tr eps = 0"
    worst, signal, xmax = 0.0, 0.0, 0.0
    for s in sts:
        x = float(np.trace(s.eps_e)) + float(s.eps_p_v)
        want = v0 * math.exp(x)
        worst = max(worst, abs(s.v - want) / v0)
        signal = max(signal, abs(v0 * (1.0 + x) - want) / v0)
        xmax = max(xmax, abs(x))
    print(f"\n[{oracle_name}/{case}] max |v - v0 exp(tr eps)| / v0 = {worst:.2e} over {len(sts)} states; max |tr eps| = {xmax:.3e}; "
          f"the linear rule would miss by {signal:.2e}")
    assert worst <= V_IDENTITY_TOL[oracle_name], f"v != v0 exp(tr eps): {worst:.3e}"
    if spec["kind"] == "undrained":
        assert xmax <= 1e-8 and signal <= 1e-12        # isochoric: x = 0 and v = v0 (O1: ODE-level trace)
    else:
        assert signal >= V_IDENTITY_POWER, f"gate has no power on {case}: max |tr eps| = {xmax:.2e}"


# ==============================================================================================
# K1.9  dissipation, refusal rule, forced counterexample (sheet 11, 13.9)
# ==============================================================================================
def step_scale(cum, k):
    return abs(cum[k] - (cum[k - 1] if k else 0.0))


@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper", "UND_paper", "UND_fork", "OC_TXE_ok"])
def test_k1_9_dissipation_nonnegative_every_step(oracle_name, case):
    """D >= 0 at every plastic step and D = 0 at every elastic step (K1.9; S.38/S.39, condition A holds for these sets).
    Roundoff floor: 1e-9 |p0| |Delta eps_step| (D is stress x strain; Delta eps_step = axial step)."""
    P, kw, st0, sts, cum = run_case(oracle_name, case)
    assert_path_ok(oracle_name, sts, cum, case)
    n_pl = 0
    for k, st in enumerate(sts):
        tol = 1e-9 * abs(kw["p0"]) * step_scale(cum, k)
        if is_plastic(st):
            n_pl += 1
            assert st.D >= -tol, f"step {k}: D = {st.D:.3e}"
        else:
            assert abs(st.D) <= 1e-15 * abs(kw["p0"]), f"elastic step {k}: D = {st.D:.3e}"
    assert n_pl >= 10, "path is mostly elastic: not a dissipation census"


@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper", "UND_paper", "OC_TXE_ok", "OC_TXE_bad"])
def test_k1_9_o2_step_dissipation_equals_the_S38_closed_form(case):
    """O2 only (backward Euler: D_step = Delta lambda * sigma_{n+1}:q_{n+1} exactly).  Sheet 11.1 before the substitution at
    yield: sigma:q = p beta F_p + q zeta_bar, with Delta lambda = Delta eps^p_s / (sqrt(2/3) Omega) (S.20).  Everything recomputed
    here from the state.  Tolerance 1e-8 relative (local Newton tolerance)."""
    oname = "O2"
    P, kw, st0, sts, cum = run_case(oname, case)
    assert_path_ok(oname, sts, cum, case)
    corner = CASES[case]["corner"]
    prev, n = st0, 0
    for k, st in enumerate(sts):
        if is_plastic(st):
            ds = st.eps_p_s - prev.eps_p_s
            dlam = ds / (SQ23 * omega_corner(kw, st.sigma, corner))
            expect = dlam * dissipation_per_lambda(kw, st.sigma, st.pi_i, corner)
            assert abs(st.D - expect) <= 1e-8 * max(abs(expect), 1e-12), (k, st.D, expect)
            n += 1
        prev = st
    assert n >= 5


# ---- the refusal rule -------------------------------------------------------------------------
REFUSED = [
    ("N_bar > N", dict(PAPER, N=0.3, N_bar=0.35, rho=0.8, rho_bar=0.8)),
    ("N_bar > N by 1%", dict(PAPER, N=0.4, N_bar=0.41, rho=0.8, rho_bar=0.8)),
    ("owner counterexample N_bar=N, rho=0.7, rho_bar=0.8", dict(PAPER, N=0.4, N_bar=0.4, rho=0.7, rho_bar=0.8)),
    ("rho/rho_bar = 0.74 < beta = 0.75", dict(PAPER, N=0.4, N_bar=0.2, rho=0.592, rho_bar=0.8)),
    ("AB06 typical: N=.3, N_bar=.2, rho=.71, rho_bar=.9", dict(PAPER, M=1.33, N=0.3, N_bar=0.2, rho=0.71, rho_bar=0.9)),
]
ACCEPTED = [
    ("K2 case 2 (rho < rho_bar)", PAPER),
    ("N_bar = N with rho = rho_bar", dict(PAPER, N=0.4, N_bar=0.4, rho=0.7, rho_bar=0.7)),
    ("rho/rho_bar = 0.76 > beta = 0.75", dict(PAPER, N=0.4, N_bar=0.2, rho=0.608, rho_bar=0.8)),
    ("fork CSL", FORK),
]


@pytest.mark.parametrize("label,kw", REFUSED)
def test_k1_9_params_validation_refuses(oracle_name, label, kw):
    """Owner-approved rule (sheet 11.2, G0 decision 2): refuse N_bar > N and rho/rho_bar < (1-N)/(1-N_bar)."""
    with pytest.raises(ValueError):
        build(oracle_name, **kw)


@pytest.mark.parametrize("label,kw", ACCEPTED)
def test_k1_9_params_validation_accepts_condition_A(oracle_name, label, kw):
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # no warning either: rho <= rho_bar
        build(oracle_name, **kw)


def test_k1_9_rho_greater_than_rho_bar_is_only_a_warning(oracle_name):
    """rho > rho_bar violates AB06's psi_c <= phi_c reading (sheet 11.2) but not dissipation under reading A: a warning,
    not a refusal.  (N = .4, N_bar = .2: rho/rho_bar = 1.14 >= beta = 0.75.)"""
    kw = dict(PAPER, rho=0.8, rho_bar=0.7)
    with pytest.warns(UserWarning):
        P = build(oracle_name, **kw)
    assert P.rho == 0.8 and P.rho_bar == 0.7


@pytest.mark.parametrize("row", [
    # (M, N, N_bar, rho, rho_bar, condition A holds, shape) -- the four WW rows of the table in sheet 11.2, then four GA rows
    # (rho, rho_bar >= 7/9, the GA range): the same four situations (pass; N_bar = N with rho < rho_bar; rho > rho_bar; rho/rho_bar < beta)
    (1.2, 0.4, 0.2, 0.7, 0.8, True, "WW"),
    (1.2, 0.4, 0.4, 0.7, 0.8, False, "WW"),
    (1.2, 0.4, 0.2, 0.8, 0.7, True, "WW"),
    (1.33, 0.3, 0.2, 0.71, 0.9, False, "WW"),
    (1.2, 0.4, 0.2, 0.8, 0.85, True, "GA"),
    (1.2, 0.4, 0.4, 0.8, 0.85, False, "GA"),
    (1.2, 0.4, 0.2, 0.85, 0.8, True, "GA"),
    (1.33, 0.3, 0.2, 0.78, 0.95, False, "GA"),
])
def test_k1_9_refusal_is_exactly_the_sign_of_the_dissipation_bracket(oracle_name, row):
    """S.38: D^p = -lambda' p [ (M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta ], -lambda' p >= 0.  The bracket is linear in eta on
    [0, M/N] and zeta_bar/zeta(theta) is monotone (WW and GA: for GA it is a Moebius function of cos 3 theta); its minimum over the
    grid (eta x theta, exact zeta from S.11 / S.10, both end points included) is negative iff condition A fails (S.39).  WW rows:
    the sheet's printed minima (+0.375, -0.375, +0.750, -0.382) are checked to 1e-3.  GA rows (nothing printed): the minimum is the
    closed form of sheet 11.1, min(M/(1 - N_bar), (M - M/N)/(1 - N_bar) + min(rho/rho_bar, 1) M/N), exact to 1e-9 (the grid contains
    both end points and the bracket is linear in eta).  The oracle must refuse exactly the failing rows."""
    M, N, Nb, rho, rhob, cond, shape = row
    eta = np.linspace(0.0, M / N, 401)[:, None]
    th = np.linspace(0.0, math.pi / 3.0, 301)[None, :]
    ratio = zeta_sheet(shape, th, rhob) / zeta_sheet(shape, th, rho)
    bracket = (M - eta) / (1.0 - Nb) + ratio * eta
    mn = float(bracket.min())
    printed = {(1.2, 0.4, 0.2, 0.7, 0.8): 0.375, (1.2, 0.4, 0.4, 0.7, 0.8): -0.375,
               (1.2, 0.4, 0.2, 0.8, 0.7): 0.750, (1.33, 0.3, 0.2, 0.71, 0.9): -0.382}
    if shape == "WW":
        assert abs(mn - printed[row[:5]]) <= 1e-3, (mn, printed[row[:5]])
    else:
        closed = min(M / (1.0 - Nb), (M - M / N) / (1.0 - Nb) + min(rho / rhob, 1.0) * M / N)
        assert abs(mn - closed) <= 1e-9, (mn, closed)
    assert (mn >= -1e-12) == cond
    kw = dict(PAPER, M=M, N=N, N_bar=Nb, rho=rho, rho_bar=rhob, zeta=shape)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if cond:
            build(oracle_name, **kw)
        else:
            with pytest.raises(ValueError):
                build(oracle_name, **kw)


def test_k1_9_forced_counterexample_negative_dissipation_at_the_extension_apex(oracle_name):
    """Sheet 11: with N_bar = N = 0.4, rho = 0.7, rho_bar = 0.8 (refused by the parser) and validation bypassed, the bracket of
    S.38 at the extension corner (zeta_bar/zeta = rho/rho_bar) is (M - eta)/(1 - N_bar) + (rho/rho_bar) eta, negative for
    eta > eta_c = [M/(1-N_bar)] / [1/(1-N_bar) - rho/rho_bar] = 2.526 (< M/N = 3), i.e. for |p|/|pi_i| < x_c with
    x_c from the inverse of S.12 at eta_c (0.135).  A generic TXE path never gets there (O1 reports); the path used is a drained
    extension from an isotropic state that is heavily over-consolidated (pi_i0 = -1000, p0 = -100): yield occurs at x ~ 0.06, at
    theta = 0, at eta ~ 2.7 > eta_c, and pi_i then collapses toward |p| at rate h.
    Asserted: (1) negative dissipation appears; (2) every step where the closed-form bracket (from the state) is negative has
    D < 0 and vice versa (O2: at the end state, exact; O1: where the bracket has that sign at both ends of the increment);
    (3) every negative step is at theta = 0 with x < x_c (O1: at the start or the end of the increment); (4) the same path with condition A satisfied (N_bar = .2, case OC_TXE_ok,
    K1.9 census above) has D >= 0 throughout."""
    P, kw, st0, sts, cum = run_case(oracle_name, "OC_TXE_bad")
    assert_path_ok(oracle_name, sts, cum, "OC_TXE_bad")
    r = kw["rho"] / kw["rho_bar"]
    eta_c = (kw["M"] / (1.0 - kw["N_bar"])) / (1.0 / (1.0 - kw["N_bar"]) - r)
    assert eta_c < kw["M"] / kw["N"]
    x_c = ratio_from_eta(kw, eta_c)
    corner = "E"
    prev, neg = st0, []
    for k, st in enumerate(sts):
        if is_plastic(st):
            g_end = dissipation_per_lambda(kw, st.sigma, st.pi_i, corner)
            tol = 1e-9 * abs(kw["p0"]) * step_scale(cum, k)
            if oracle_name == "O2":
                if abs(g_end) > 1e-9:
                    assert (st.D < 0.0) == (g_end < 0.0), (k, st.D, g_end)
            elif is_plastic(prev):
                g_start = dissipation_per_lambda(kw, prev.sigma, prev.pi_i, corner)
                if g_start < -1e-9 and g_end < -1e-9:
                    assert st.D < 0.0, (k, st.D, g_start, g_end)
                if g_start > 1e-9 and g_end > 1e-9:
                    assert st.D > -tol, (k, st.D, g_start, g_end)
            if st.D < -tol:
                p, q, th, _ = pq_theta(st.sigma)
                neg.append(k)
                assert abs(th) <= 1e-6, (k, th)
                # O2: the end state is the state.  O1: D is the integral over the increment, which may start in the negative zone
                # (x rises as pi_i collapses), so the smaller x of the two bracketing states must be below x_c.
                xs = [abs(p / st.pi_i)] + ([abs(pq_theta(prev.sigma)[0] / prev.pi_i)] if oracle_name == "O1" else [])
                assert min(xs) < x_c * (1.0 + 1e-6), (k, xs, x_c)
        prev = st
    assert neg, "no negative dissipation appeared: the forced counterexample did not reproduce (sheet 11)"
