"""WP-144 gate G1 (Zone B): the Q-cap closed forms (equation sheet 144a 3.2, 10, 11; plan 5.1-5.2).

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form in the
(revised) sheet, with the section cited at the use, or (b) a stated convergence / round-off argument, and was
written BEFORE the oracles were run on these paths.  No oracle output is used as an expected value.  A tolerance
is never loosened to make a test pass; a failing assertion is a finding (against the sheet or an oracle), not a
reason to edit an oracle or a bound.

What is tested (all on the K2 paper set, WW, rho = 0.7 / rho_bar = 0.8, near-isotropic compression path of
test_g1_convergence_tangents.e_iso_total):
  (a) smooth cap, O2 : Delta eps^p = Delta lambda q_a(S.35) at every plastic step, with w and eta RECOMPUTED here
  (b) planar cap, O2 : on the cap q_a = -delta_a/3, eps^p_s = 0, pi_i frozen exactly
  (c) capped dissipation closed form, O2 : D = dlam [-p(1-w) + w D_u/lam]   (sheet 10.2, 11.1, S.38)
  (d) O2 consistent tangent vs central FD of O2's own stress update, cap ACTIVE (0 < w < 1), off the theta corners
  (e) near-isotropic smooth-cap paths (deviator 1e-4, 5e-4, 2e-3; 40 steps) complete in BOTH oracles, D >= bound

Every flow quantity is re-derived here from the sheet in numpy / sympy (zeta from S.10-S.11, y_a from S.6, q_a from
S.17 and S.35, F_p from S.13); the oracles supply only the state (stress, pi_i, accumulated plastic strain, dlam, D).
"""
import functools
import math
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from conftest import K2_BASE, ORACLES, make_params

O1 = ORACLES["O1"]
O2 = ORACLES["O2"]

SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)
ONES = np.ones(3)

SIG0 = -100.0 * I3          # isotropic start at p0, eps^e = 0
PI0 = -80.0                 # initial image pressure (apex of the surface at pi_c = -80/0.6^1.5 = -172 kPa)
PSI0 = -0.05                # initial image state parameter (dense): v0 = psi0 + v_c0 - lambda_tilde ln(-pi0)
P0 = 100.0
AMP_STOP = 2.0e-3           # the near-isotropic path on which the planar cap and no cap stop (sheet 3.2, 10.1)
N_CAP = 40
R_TOL_REL = 1.0e-8          # R_tol of the vertex rule, sheet 3.2

CAPS = {
    "smooth": dict(cap="smooth", c1=0.05, c2=0.15),       # sheet 10.2 recommended defaults
    "planar": dict(cap="planar", c1=0.10, c2=0.10),       # BA06's chi_cap = 0.10 (sheet 10.1)
    "none": dict(cap="none"),
}


def kw_for(cap):
    """AB06 6.1 / K2 case 2 (rho = 0.7, rho_bar = 0.8; condition A holds: 0.875 >= beta = 0.75) with the cap mode."""
    return dict(K2_BASE, rho=0.7, rho_bar=0.8, **CAPS[cap])


def v0_for(kw):
    """Specific volume giving psi_i = PSI0 at pi_i = PI0 (sheet 6, paper CSL)."""
    return PSI0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-PI0)


def e_iso_total(amp):
    """Total strain of the near-isotropic path: -0.01 on every axis (tr = -0.03) plus the deviator amp*diag(1, 0, -1)."""
    return -0.01 * I3 + amp * np.diag([1.0, 0.0, -1.0])


# ----------------------------------------------------------------------------------------------
# the sheet's closed forms, written here (numpy / sympy), independent of both oracles
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _zeta_fns(kind):
    """zeta(theta, rho) and d zeta/d theta: WW (S.11) / GA (S.10), exact via sympy, lambdified."""
    th, rho = sp.symbols("theta rho", positive=True)
    if kind == "WW":
        c = sp.cos(th)
        A = 4 * (1 - rho ** 2)
        B = 2 * rho - 1
        z = (A * c ** 2 + B ** 2) / (2 * (1 - rho ** 2) * c + B * sp.sqrt(A * c ** 2 + 5 * rho ** 2 - 4 * rho))
    else:
        z = ((1 + rho) + (1 - rho) * sp.cos(3 * th)) / (2 * rho)
    return (sp.lambdify((th, rho), z, "math"), sp.lambdify((th, rho), sp.diff(z, th), "math"))


def eta_S12(kw, p, pi):
    """eta(p, pi_i) of (S.12)."""
    M, N = kw["M"], kw["N"]
    if N == 0.0:
        return M * (1.0 + math.log(pi / p))
    return (M / N) * (1.0 - (1.0 - N) * (p / pi) ** (N / (1.0 - N)))


def cap_weight_S35(kw, eta):
    """w(eta) of (S.35): quintic blend between eta_1 = c1 M and eta_2 = c2 M; planar = step with w = 1 at eta = eta_1
    (the convention of sheet 10.1); no cap = 1."""
    if kw["cap"] == "none":
        return 1.0
    e1, e2 = kw["c1"] * kw["M"], kw["c2"] * kw["M"]
    if kw["cap"] == "planar":
        return 1.0 if eta >= e1 else 0.0
    t = min(1.0, max(0.0, (eta - e1) / (e2 - e1)))
    return t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)


def theta_of(sig):
    """Lode angle from (S.1): y = tr(xi^3)/R^3 = cos3theta/sqrt6, theta = acos(sqrt6 y)/3 in [0, pi/3]."""
    lam = np.linalg.eigvalsh(0.5 * (sig + sig.T))
    xi = lam - lam.mean()
    R = float(np.linalg.norm(xi))
    y = float(np.sum(xi ** 3)) / R ** 3
    return math.acos(max(-1.0, min(1.0, SQ6 * y))) / 3.0


def sheet_flow(kw, sigma, pi):
    """Everything the cap flow needs, from a stress tensor and pi_i, by the sheet's closed forms.
    Principal values in eigh (ascending) order, V the eigenvectors; q_a is (S.35) with q^u_a = (S.17):
        q^u_a = beta F_p / 3 + sqrt(3/2) zeta_bar n_a + zeta_bar_y q y_a,   beta F_p = (eta - M)/(1 - N_bar)   (S.13, S.16)
        q_a   = -1/3 + w (q^u_a + 1/3)
    Vertex rule (sheet 3.2): R < R_tol  =>  n_a = y_a = 0, q = 0.   D_u / lambda = sigma:q^u = p beta F_p + q zeta_bar,
    and, on the surface (F = 0, q = -p eta/zeta), D_u / lambda = -p [(M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta]  (S.38)."""
    sig = 0.5 * (sigma + sigma.T)
    lam, V = np.linalg.eigh(sig)
    p = float(lam.mean())
    xi = lam - p
    R = float(np.linalg.norm(xi))
    q = SQ32 * R
    M, Nb = kw["M"], kw["N_bar"]
    eta = eta_S12(kw, p, pi)
    betaFp = (eta - M) / (1.0 - Nb)
    w = cap_weight_S35(kw, eta)
    out = SimpleNamespace(lam=lam, V=V, p=p, q=q, R=R, eta=eta, betaFp=betaFp, w=w, vertex=R < R_TOL_REL * abs(p),
                          theta=None, zeta=None, zeta_bar=None, Du=None, Du_yield=None)
    if out.vertex:
        qu = betaFp / 3.0 * ONES
        out.q = 0.0
        out.Du = p * betaFp
    else:
        nh = xi / R
        S3 = float(np.sum(xi ** 3))
        y = S3 / R ** 3
        theta = math.acos(max(-1.0, min(1.0, SQ6 * y))) / 3.0
        assert 1e-3 < theta < math.pi / 3.0 - 1e-3, \
            f"state within 1e-3 of a theta corner (theta = {theta}): zeta_y must not be tested there (sheet 3.1 rule ii)"
        y_a = 3.0 * xi ** 2 / R ** 3 - 3.0 * S3 * xi / R ** 5 - 1.0 / R                                  # (S.6)
        zfun, zfun_theta = _zeta_fns(kw["zeta"])
        z = zfun(theta, kw["rho"])
        zb = zfun(theta, kw["rho_bar"])
        zb_theta = zfun_theta(theta, kw["rho_bar"])
        zb_y = -(2.0 / SQ6) / math.sin(3.0 * theta) * zb_theta                                           # (S.8)
        qu = betaFp / 3.0 * ONES + SQ32 * zb * nh + zb_y * q * y_a                                       # (S.17)
        out.theta, out.zeta, out.zeta_bar = theta, z, zb
        out.Du = p * betaFp + q * zb
        out.Du_yield = -p * ((M - eta) / (1.0 - Nb) + (zb / z) * eta)                                    # (S.38)
    out.qu = qu
    out.q_a = -ONES / 3.0 + w * (qu + ONES / 3.0)                                                        # (S.35)
    return out


def dep_from_q(f, dlam):
    """Delta eps^p = dlam sum_a q_a n^a (x) n^a, with n^a the eigenvectors of sigma (coaxial, sheet 5.3)."""
    return dlam * (f.V * f.q_a) @ f.V.T


def bracket_min(kw):
    """Lower bound of the (S.38) bracket B(eta, theta) = (M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta over the whole
    surface (sheet 11.1): B is LINEAR in eta on [0, M/N], so its minimum is at an end point, and min_theta zeta_bar/zeta =
    min(rho/rho_bar, 1) (monotone in theta, minimum rho/rho_bar at theta = 0 when rho < rho_bar).  K2: min(1.5, 0.375)."""
    M, N, Nb = kw["M"], kw["N"], kw["N_bar"]
    r = min(kw["rho"] / kw["rho_bar"], 1.0)
    e_max = M / N
    return min(M / (1.0 - Nb), (M - e_max) / (1.0 - Nb) + r * e_max)


def check_o2_step(kw, prev, st, rel=1e-8):
    """One plastic O2 step against the sheet: Delta eps^p = dlam q_a(S.35) (relative 'rel' of ||Delta eps^p||) and, where
    the effective Omega = w Omega^u is exactly zero (w = 0 or the vertex rule), pi_i and eps^p_s are exactly frozen (S.25,
    sheet 10.1).  Returns (flow quantities, relative error)."""
    f = sheet_flow(kw, st.sigma, st.pi_i)
    dep = st.eps_p - prev.eps_p
    exp = dep_from_q(f, st.dlam)
    scale = max(float(np.linalg.norm(dep)), float(np.linalg.norm(exp)), 1e-300)
    err = float(np.linalg.norm(dep - exp)) / scale
    assert err <= rel, (f"Delta eps^p != dlam q_a: rel err {err:.3e} (w = {f.w:.6g}, eta = {f.eta:.6g}, "
                        f"dlam = {st.dlam:.6e})")
    if f.w == 0.0 or f.vertex:
        assert st.pi_i == prev.pi_i, f"pi_i moved with Omega = 0: {prev.pi_i} -> {st.pi_i} (w = {f.w}, vertex {f.vertex})"
        assert st.eps_p_s == prev.eps_p_s, "eps^p_s grew with Omega = 0"
    return f, err


# ----------------------------------------------------------------------------------------------
# runners
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def run_o2(cap, amp, n):
    kw = kw_for(cap)
    P = make_params("O2", **kw)
    st0 = O2.initial_state(P, SIG0, v0_for(kw), PI0)
    deps = np.array([e_iso_total(amp) / n] * n)
    return kw, P, st0, O2.run_path(P, st0, deps)


def iter_plastic(st0, sts):
    """(k, previous state, state) for every non-refused plastic step with dlam > 0 that was NOT substepped (on a
    substepped increment State.dlam is the last sub-step's and eps_p / D are sums: the per-step closed forms below
    are per backward-Euler step).  Stops at the first refusal."""
    prev = st0
    for k, st in enumerate(sts):
        if st.flags["refused"]:
            break
        if st.flags["plastic"] and st.flags.get("substeps", 1) == 1 and st.dlam > 0.0:
            yield k, prev, st
        prev = st


def refused_at(sts):
    r = [s.flags["refused"] for s in sts]
    return r.index(True) if any(r) else None


# ==============================================================================================
# (a) smooth cap, O2: Delta eps^p = dlam q_a(S.35) at every plastic step, w and eta recomputed
# ==============================================================================================
@pytest.mark.parametrize("amp,n", [(AMP_STOP, 160), (6.0e-3, 160), (1.0e-2, 160)])
def test_g1_cap_smooth_o2_plastic_flow_is_the_S35_blend(amp, n):
    """Sheet 10.2 (S.35): q_a = -delta_a/3 + w (q^u_a + delta_a/3), w = S(t), t = clamp((eta - eta_1)/(eta_2 - eta_1)),
    S(t) = t^3 (10 - 15 t + 6 t^2), eta = eta(p, pi_i) of (S.12), q^u_a = (S.17).  At every plastic step the accumulated
    plastic strain increment must equal Delta lambda q_a evaluated at the CONVERGED state (backward Euler:
    Delta eps^p = Delta lambda q(sigma_{n+1}, pi_{i,n+1})), to 1e-8 relative to ||Delta eps^p||: the two sides are
    the same algebra on the same converged state, so only round-off (1e-12) and the 1e-12 nested tolerance on pi_i
    enter; 1e-8 is the gate.  w, eta, zeta_bar, zeta_bar_y, y_a are recomputed HERE (sympy zeta, S.6, S.8), none is read
    from the oracle.  Where w = 0 (or the vertex) pi_i and eps^p_s must be exactly frozen (S.25 with Omega = 0).
    Path choice (scenario, not a value): n = 160 increments so that the steps are not substepped (the per-step closed
    form needs one backward-Euler step per increment; at n = 40 most plastic steps are substepped); deviator 2e-3 sits at
    the lower edge of the ramp, 6e-3 crosses its middle (w from 0.8 to 0.2), 1e-2 comes down from the uncapped branch
    (w = 1) into the ramp.  Path guard: the path must really visit the ramp, and its mixed region (min(w, 1-w) >= 0.05),
    or the test cannot tell the blend from the planar flow or from the uncapped flow.
    Mutants killed: q_a = q^u (cap ignored on the ramp); q_a = -delta/3 on the ramp (planar); a cubic 3t^2 - 2t^3 or a
    linear blend; the (1-w) term dropped from the volumetric part; the blend evaluated at the trial eta instead of the
    converged one."""
    kw, P, st0, sts = run_o2("smooth", amp, n)
    k_ref = refused_at(sts)
    assert k_ref is None, f"O2 refused the smooth-cap path at step {k_ref + 1}: {sts[k_ref].flags['reason']}"
    ws, errs = [], []
    for k, prev, st in iter_plastic(st0, sts):
        f, err = check_o2_step(kw, prev, st, rel=1e-8)
        errs.append(err)
        if 0.0 < f.w < 1.0:
            ws.append(f.w)
    print(f"\n[smooth a amp={amp} n={n}] plastic steps checked {len(errs)}, ramp steps {len(ws)}, "
          f"max w {max(ws) if ws else float('nan'):.4f}, max rel err {max(errs):.3e}")
    assert len(ws) >= 1, "test bug / finding: the path never enters the ramp 0 < w < 1"
    assert max(min(w, 1.0 - w) for w in ws) >= 0.05, f"the ramp steps all sit at the ramp edge: max w = {max(ws)}"


# ==============================================================================================
# (b) planar cap, O2: on the cap q_a = -delta_a/3, eps^p_s = 0, pi_i frozen exactly
# ==============================================================================================
@pytest.mark.parametrize("amp", [1.0e-4, AMP_STOP], ids=["dev1e-4", "dev2e-3"])
def test_g1_cap_planar_o2_on_the_cap_flow_is_volumetric_and_pi_frozen(amp):
    """Sheet 10.1 (BA06 2.76): for eta < chi_cap M the flow potential is Q = -p: q_a = -delta_a/3, q_ab = 0, Omega = 0,
    hence eps^p_s does not grow, pi_i does not evolve (S.25: pi_i' ~ Omega), and eps^p_v = -lambda (compaction).
    Asserted at every plastic, non-substepped step before any refusal whose CONVERGED eta (recomputed from S.12) is
    < c1 M:
      * Delta eps^p = -(Delta lambda/3) 1: relative 1e-12 (the eigenvector basis V V^T = 1 to round-off);
      * Delta eps^p_v = -Delta lambda (1e-12 Delta lambda);
      * eps^p_s unchanged and pi_i == pi_{i,n}, both EXACT (Omega = w Omega^u is exactly 0.0 there: no tolerance);
    and where eta >= c1 M (w = 1, convention of sheet 10.1) the uncapped flow (S.17) to 1e-8.
    n = 160 increments (not substepped: the per-step relations need one backward-Euler step per increment).
    Guard: at least 3 cap steps are checked.
    Mutants killed: deviatoric flow kept on the cap; pi_i hardening not frozen; q_a = -delta (missing 1/3) or +delta/3
    (sign); the planar switch inverted (eta > c1 M); the convention at the switch changed to w = 0 for eta >= c1 M."""
    kw, P, st0, sts = run_o2("planar", amp, 160)
    e1 = kw["c1"] * kw["M"]
    n_cap, n_unc = 0, 0
    for k, prev, st in iter_plastic(st0, sts):
        f = sheet_flow(kw, st.sigma, st.pi_i)
        if f.eta < e1:
            assert f.w == 0.0
            dep = st.eps_p - prev.eps_p
            dl = st.dlam
            assert np.linalg.norm(dep + dl / 3.0 * I3) <= 1e-12 * np.linalg.norm(dep), (k, dep, dl)
            assert abs((st.eps_p_v - prev.eps_p_v) + dl) <= 1e-12 * dl, (k, st.eps_p_v - prev.eps_p_v, dl)
            assert st.eps_p_s == prev.eps_p_s, f"step {k + 1}: eps^p_s grew on the cap"
            assert st.pi_i == prev.pi_i, f"step {k + 1}: pi_i moved on the cap: {prev.pi_i} -> {st.pi_i}"
            n_cap += 1
        else:
            check_o2_step(kw, prev, st, rel=1e-8)
            n_unc += 1
    print(f"\n[planar b amp={amp}] cap steps {n_cap}, uncapped steps {n_unc}, refused at {refused_at(sts)}")
    assert n_cap >= 3, f"only {n_cap} plastic steps on the cap: the path does not exercise the planar cap"


# ==============================================================================================
# (c) capped dissipation closed form, O2
# ==============================================================================================
@pytest.mark.parametrize("cap,amp,n", [("smooth", AMP_STOP, 160), ("smooth", 6.0e-3, 160), ("smooth", 1.0e-2, 160),
                                       ("planar", AMP_STOP, 160), ("none", AMP_STOP, 160)])
def test_g1_cap_o2_step_dissipation_is_the_capped_closed_form(cap, amp, n):
    """Sheet 10.2 / 11.1: D^p = lambda [ -p (1 - w) + w D_u / lambda ], with the uncapped dissipation per unit multiplier
    D_u / lambda = sigma:q^u = p beta F_p + q zeta_bar and, ON the yield surface (F = 0), by (S.38)
    D_u / lambda = -p [ (M - eta)/(1 - N_bar) + (zeta_bar/zeta) eta ].  Per backward-Euler step D = Delta lambda x that,
    evaluated at the converged state with w, eta, p, q recomputed here.  Two checks per plastic non-substepped step:
      * the two forms of D_u agree (that is F = 0 at the converged state): |difference| <= 1e-9 |p|  (F residual is
        <= 1e-12 |p0| = 1e-10 kPa times zeta_bar/zeta);
      * State.D == Delta lambda [ -p(1-w) + w D_u/lambda ]: 1e-9 Delta lambda |p| (round-off of sigma:q against the
        identities sum_a sigma_a y_a = 0, sum_a sigma_a n_a = R, sheet 3.1).
    D = 0 EXACTLY on every elastic step.  Guard: at least 3 plastic steps; the smooth runs must include a ramp step.
    Mutants killed: D from the uncapped flow (w ignored); the planar term -p(1-w) dropped or sign-flipped; zeta_bar and
    zeta swapped in S.38; D reported with q^u at the trial state."""
    kw, P, st0, sts = run_o2(cap, amp, n)
    n_pl, n_ramp = 0, 0
    for k, st in enumerate(sts):
        if st.flags["refused"]:
            break
        if not st.flags["plastic"]:
            assert st.D == 0.0, f"step {k + 1}: D = {st.D} on an elastic step"
    for k, prev, st in iter_plastic(st0, sts):
        f = sheet_flow(kw, st.sigma, st.pi_i)
        if not f.vertex:
            assert abs(f.Du - f.Du_yield) <= 1e-9 * abs(f.p), \
                f"step {k + 1}: D_u forms differ ({f.Du} vs {f.Du_yield}): the converged state is off the surface"
        Du = f.Du
        D_exp = st.dlam * (-f.p * (1.0 - f.w) + f.w * Du)
        assert abs(st.D - D_exp) <= 1e-9 * st.dlam * abs(f.p), \
            f"step {k + 1}: D = {st.D:.12e}, closed form {D_exp:.12e} (w = {f.w:.5g}, dlam = {st.dlam:.4e})"
        n_pl += 1
        n_ramp += 0.0 < f.w < 1.0
    print(f"\n[D closed form {cap} n={n}] plastic steps {n_pl}, ramp steps {n_ramp}, refused at {refused_at(sts)}")
    assert n_pl >= 3
    if cap == "smooth":
        assert n_ramp >= 1, "the smooth path never enters the ramp"


# ==============================================================================================
# (d) O2 CTO vs central FD, cap ACTIVE (0 < w < 1), off the theta corners
# ==============================================================================================
def _col_errors(P, A, deps, h):
    """Per-column relative error of the consistent tangent at O2.step(A, deps) against the central FD of the O2 stress
    update, columns (k, l), k <= l, perturbation E_kl = E_lk = 1/2 (d sigma/dh = C_ij(kl) by the minor symmetry of the
    small-strain tangent, sheet 9.4).  Returns ({(k, l): err}, the unperturbed step), or (None, the unperturbed step) when
    an FD point is not a single plain backward-Euler step (plastic, not refused, not substepped): the FD of a substepped
    increment is the FD of a different map (the chain of sub-steps), not of the one whose tangent is being tested."""
    stn = O2.step(P, A, deps)
    assert stn.flags["plastic"] and not stn.flags["refused"] and stn.flags.get("substeps", 1) == 1
    C = O2.tangent(P, stn)
    errs = {}
    for k in range(3):
        for l in range(k, 3):
            E = np.zeros((3, 3))
            if k == l:
                E[k, k] = 1.0
            else:
                E[k, l] = E[l, k] = 0.5
            sp_, sm_ = O2.step(P, A, deps + h * E), O2.step(P, A, deps - h * E)
            for s in (sp_, sm_):
                if not (s.flags["plastic"] and not s.flags["refused"] and s.flags.get("substeps", 1) == 1):
                    return None, stn
            fd = (sp_.sigma - sm_.sigma) / (2.0 * h)
            errs[(k, l)] = float(np.linalg.norm(C[:, :, k, l] - fd) / np.linalg.norm(C[:, :, k, l]))
    return errs, stn


@pytest.mark.parametrize("amp", [6.0e-3, 1.0e-2], ids=["dev6e-3", "dev1e-2"])
def test_g1_cap_o2_cto_matches_central_fd_with_the_cap_active(amp):
    """Plan 5.1 / sheet 9.3, 10.2 (S.36-S.37): the closed-form consistent tangent, built from q_ab = w q^u_ab +
    w_eta eta_p g_a/3 delta_b, q_{a,pi} = w q^u_{a,pi} + w_eta eta_pi g_a, Omega_a = w Omega^u_a + w_eta eta_p Omega^u/3 and
    the capped nested-solve derivative r' of (S.37), must match the central difference of O2's own stress update to 1e-6
    per column (best step of {1e-5, 1e-6, 1e-7}, the same gate as G1.fd), including the three shear columns.  A step size
    counts only if EVERY FD point is a single plain backward-Euler step (plastic, not refused, not substepped): a central
    difference across a substepped increment differentiates the chain of sub-steps, not the map whose tangent is tested.
    At h = 1e-5 (16% of the 6e-5 increment) the compressive diagonal perturbations push the loop gain G of the nested solve
    past 1 on several states and O2 substeps them (first attempt 'local_linesearch:pi_fold', sheet 10.2 non-uniqueness); those
    h are reported and not used, and h = 1e-6 and 1e-7 are required to be valid on every state.
    Gate states: the smooth-cap near-isotropic path of deviator 6e-3 (w from 0.8 to 0.2, across the middle of the ramp where
    w_eta is largest) and 1e-2 (down from w = 1) at n = 160 (a fine step: the loop gain G of the nested pi_i solve falls with
    the step, sheet 10.2, so r' stays away from the fold where a central difference would straddle the 0.05 kPa dip width).
    Selection rule, fixed before the run: every step k whose
    CONVERGED state has 0.02 <= w <= 0.98 (recomputed here), is plastic, single-step, with theta in [0.1, pi/3 - 0.1]
    (the sheet 4.3 rule for FD tangent tests), then up to 6 of them spread evenly; the FD base state is the state before
    the step and the increment is the path increment.  Truncation argument: eta changes by ~2e-3 per 1e-5 of strain
    (eta_p dp, dp = K d eps = 0.17 kPa) against a ramp width of 0.12, so (h/scale)^2/6 ~ 1e-7 at h = 1e-5 and 1e-9 at
    1e-6; round-off eps sigma/h ~ 1e-11 relative at 1e-7.
    Mutants killed: any w_eta term dropped from q_ab, q_{a,pi}, Omega_a or r' (sheet 10.2 says a w_eta-dropping mutant
    is caught at w = 0.19, 0.51); the (1-w)/w weighting of q_ab; a stale pi_i derivative P_a (Omega_pi missing)."""
    kw, P, st0, sts = run_o2("smooth", amp, 160)
    assert refused_at(sts) is None, f"O2 refused the n = 160 smooth path: {sts[refused_at(sts)].flags['reason']}"
    states = [st0] + list(sts)
    cands = []
    for i, nxt in enumerate(sts):
        if not (nxt.flags["plastic"] and nxt.flags.get("substeps", 1) == 1 and nxt.dlam > 0.0):
            continue
        f = sheet_flow(kw, nxt.sigma, nxt.pi_i)
        th = theta_of(nxt.sigma)
        if 0.02 <= f.w <= 0.98 and 0.1 < th < math.pi / 3.0 - 0.1:
            cands.append((i, f.w))
    assert len(cands) >= 2, f"only {len(cands)} usable ramp states on the path"
    idx = np.unique(np.linspace(0, len(cands) - 1, min(6, len(cands))).round().astype(int))
    deps = e_iso_total(amp) / 160
    results = {}
    for j in idx:
        i, w = cands[j]
        per_h, off_branch = {}, []
        for h in (1e-5, 1e-6, 1e-7):
            errs, stn = _col_errors(P, states[i], deps, h)
            if errs is None:
                off_branch.append(h)
            else:
                per_h[h] = max(errs.values())
        print(f"\n[cto cap] step {i + 1} (w = {w:.3f}): max-column err per h {per_h}; FD points off the single-step "
              f"branch (not used) at h = {off_branch}")
        results[(i, round(w, 3))] = per_h
    for key, per_h in results.items():
        assert 1e-6 in per_h and 1e-7 in per_h,             f"state {key}: FD points left the single-step branch even at h = 1e-6 / 1e-7 (valid h: {sorted(per_h)})"
        best = min(per_h.values())
        assert best <= 1e-6, f"state {key}: best-h max column error {best:.3e} > 1e-6; all: {per_h}"
    sel_w = [cands[j][1] for j in idx]
    assert max(min(w, 1.0 - w) for w in sel_w) >= 0.02, f"all selected states sit at the ramp edge: w = {sel_w}"


# ==============================================================================================
# (e) near-isotropic smooth-cap paths complete in BOTH oracles, D >= the sheet 11 bound
# ==============================================================================================
@pytest.mark.parametrize("amp", [1.0e-4, 5.0e-4, AMP_STOP], ids=["dev1e-4", "dev5e-4", "dev2e-3"])
def test_g1_cap_smooth_near_isotropic_paths_complete_with_the_dissipation_bound(oracle_name, amp):
    """Sheet 10.2 + 11.1 + 16.3 [G1]: with the smooth cap (quintic, c1 = 0.05, c2 = 0.15) the near-isotropic compression paths
    of deviator 1e-4, 5e-4, 2e-3 run to completion in 40 increments in BOTH oracles (O1: every increment status ok; O2: no
    refusal), after the O2 root-selection fix of sheet 8 / 10.2.  Dissipation bound (sheet 11.1, S.38-S.39 with the cap 10.2):
    per step D = dlam [-p(1-w) + w D_u/lam] with D_u/lam = -p B, B >= B_min = min over the surface (bracket_min, 0.375 for K2),
    so  D >= dlam |p| [(1 - w) + w B_min]  >= 0.375 dlam |p| > 0   at every plastic non-substepped O2 step (w, p recomputed
    here; checked at n = 40 and n = 160), D >= 0 at every step (and = 0 exactly on elastic O2 steps).  O1 integrates the same non-negative rate, so it may
    undershoot 0 by its ODE tolerance only: D >= -1e-9 |p0| eps_max (the tolerance of the existing census).
    Mutants killed: an O2 return map that refuses (or loops) in the ramp; a w that is not clamped to [0, 1]; a dissipation
    with the wrong cap weighting (D below the bound); a nested pi_i solve that jumps to the far root (path ends or D < bound)."""
    kw = kw_for("smooth")
    P = make_params(oracle_name, **kw)
    ora = ORACLES[oracle_name]
    st0 = ora.initial_state(P, SIG0, v0_for(kw), PI0)
    deps = np.array([e_iso_total(amp) / N_CAP] * N_CAP)
    sts = ora.run_path(P, st0, deps, rtol=1e-10) if oracle_name == "O1" else ora.run_path(P, st0, deps)
    if oracle_name == "O1":
        stat = sts[-1].flags["status"]
        assert len(sts) == N_CAP and all(s.flags["status"] == "ok" for s in sts), \
            f"O1 did not complete: {stat} after {len(sts)}/{N_CAP}"
    else:
        k = refused_at(sts)
        assert k is None, f"O2 refused at step {k + 1}: {sts[k].flags['reason']}"
    D = np.array([s.D for s in sts])
    plastic = np.array([bool(s.flags["plastic"]) for s in sts])
    assert plastic.any(), "the path never yields"
    if oracle_name == "O1":
        tol = -1e-9 * P0 * float(np.abs(e_iso_total(amp) / N_CAP).max())
        assert np.all(D >= tol), f"min D {D.min():.3e}"
        return
    assert np.all(D >= 0.0), f"negative dissipation {D.min():.3e}"
    assert np.all(D[~plastic] == 0.0), "D must be exactly 0 on an elastic step"
    bmin = bracket_min(kw)
    assert abs(bmin - 0.375) <= 1e-12, f"test bug: B_min = {bmin}, the K2 closed form is 0.375"
    # The bound is per backward-Euler step (needs dlam of ONE step), so it is checked on the non-substepped steps: at the
    # n = 40 gate most plastic steps are substepped (few remain), hence also on the same path at n = 160 (no substepping).
    for n_run in (N_CAP, 160):
        kw_, P_, st0_, sts_ = run_o2("smooth", amp, n_run)
        assert refused_at(sts_) is None, f"O2 n = {n_run}: refused at step {refused_at(sts_) + 1}"
        n_chk, ratio = 0, []
        for k, prev, st in iter_plastic(st0_, sts_):
            f = sheet_flow(kw_, st.sigma, st.pi_i)
            bound = st.dlam * abs(f.p) * ((1.0 - f.w) + f.w * bmin)
            assert st.D >= bound * (1.0 - 1e-9),                 f"n={n_run} step {k + 1}: D = {st.D:.6e} < bound {bound:.6e} (w = {f.w:.4g})"
            n_chk += 1
            ratio.append(st.D / bound)
        print(f"\n[smooth e O2 amp={amp} n={n_run}] plastic steps checked {n_chk}, "
              f"min D/bound {min(ratio) if ratio else float('nan'):.4f}")
        if n_run == 160:
            assert n_chk >= 3


# ==============================================================================================
# (f) the cap on the parameter sets the fork actually ships: fork CSL + smooth cap (sheet 15 TIMs mapping), and rho = rho_bar = 1
# ==============================================================================================
# Sheet 15: M = M_c = 1.3309 direct, rho = c = 0.71, fork CSL (e0 = 0.83, lambda_c = 0.027, xi = 0.45, p_a = 101.325 kPa),
# rho_bar = 0.75 (sheet 15 constraint: rho/rho_bar = 0.947 >= beta = (1-N)/(1-N_bar) = 0.75, condition A holds); smooth cap with the
# sheet 10.2 defaults. K2 base for the remaining BA06 constants (mu0, kappa_hat, h, chi, N, N_bar, lambda_tilde, v_c0 unused in fork mode).
FORK_CAP_CASES = {
    "fork_default": dict(K2_BASE, M=1.3309, rho=0.71, rho_bar=0.75, csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45,
                         p_a=101.325, cap="smooth", c1=0.05, c2=0.15),
    "rho1_rhobar1": dict(K2_BASE, rho=1.0, rho_bar=1.0, cap="smooth", c1=0.05, c2=0.15),
}


def v0_case(kw):
    """Specific volume giving psi_i = PSI0 at pi_i = PI0 in the case's own CSL mode (sheet 6, S.22)."""
    if kw["csl_mode"] == "paper":
        return PSI0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-PI0)
    return 1.0 + kw["e0"] + PSI0 - kw["lambda_c"] * (-PI0 / kw["p_a"]) ** kw["xi"]


@functools.lru_cache(maxsize=None)
def run_o2_case(label, amp, n):
    kw = FORK_CAP_CASES[label]
    P = make_params("O2", **kw)
    st0 = O2.initial_state(P, SIG0, v0_case(kw), PI0)
    deps = np.array([e_iso_total(amp) / n] * n)
    return kw, P, st0, O2.run_path(P, st0, deps)


@pytest.mark.parametrize("amp", [AMP_STOP, 6.0e-3, 1.0e-2], ids=["dev2e-3", "dev6e-3", "dev1e-2"])
@pytest.mark.parametrize("label", sorted(FORK_CAP_CASES))
def test_g1_cap_fork_sets_o2_flow_is_S35_and_dissipation_is_the_capped_closed_form(label, amp):
    """Sheet 10.2 (S.35), 11.1 (S.38), 15, on the parameter sets the fork ships: (i) fork_default = fork CSL + smooth cap (quintic,
    c1 = 0.05, c2 = 0.15), M = 1.3309, rho = 0.71, rho_bar = 0.75 (the sheet 15 TIMs mapping; the CSL does not enter the flow or the
    dissipation argument, sheet 11.4, but it does enter psi_i, pi_i* and the path); (ii) rho = rho_bar = 1 (zeta = zeta_bar = 1, so
    zeta_bar_y = 0 and the flow has no Lode term), smooth cap, paper CSL. Same two gates as the paper-set tests (a) and (c), same
    tolerances and for the same reasons: at every plastic non-substepped O2 step
      * Delta eps^p = Delta lambda q_a (S.35) with w, eta, zeta, zeta_bar, y_a recomputed here from the converged state
        (check_o2_step; relative 1e-8; exact freezing of pi_i and eps^p_s where w = 0);
      * the two forms of D_u agree (F = 0 at the converged state, 1e-9 |p|) and State.D = Delta lambda [-p (1 - w) + w D_u/lambda]
        (1e-9 Delta lambda |p|); D = 0 exactly on elastic steps.
    Paths: the three near-isotropic compression paths of (a) at n = 160 (not substepped); each must visit the ramp
    (min(w, 1 - w) >= 0.05 on some plastic step, otherwise the blend cannot be told from the planar flow or from no cap).
    Mutants killed: those of (a) and (c), plus a cap blend or dissipation that silently assumes the paper CSL or rho_bar != 1
    (a nonzero Lode term at rho_bar = 1 would break the S.35 identity)."""
    kw, P, st0, sts = run_o2_case(label, amp, 160)
    assert refused_at(sts) is None, f"O2 refused the {label} path at step {refused_at(sts) + 1}: {sts[refused_at(sts)].flags['reason']}"
    assert kw["N_bar"] <= kw["N"] and kw["rho"] / kw["rho_bar"] >= (1.0 - kw["N"]) / (1.0 - kw["N_bar"])   # condition A (S.39)
    ws, errs = [], []
    for k, prev, st in iter_plastic(st0, sts):
        f, err = check_o2_step(kw, prev, st, rel=1e-8)
        errs.append(err)
        if 0.0 < f.w < 1.0:
            ws.append(f.w)
        if not f.vertex:
            assert abs(f.Du - f.Du_yield) <= 1e-9 * abs(f.p), \
                f"step {k + 1}: D_u forms differ ({f.Du} vs {f.Du_yield}): the converged state is off the surface"
        D_exp = st.dlam * (-f.p * (1.0 - f.w) + f.w * f.Du)
        assert abs(st.D - D_exp) <= 1e-9 * st.dlam * abs(f.p), \
            f"step {k + 1}: D = {st.D:.12e}, closed form {D_exp:.12e} (w = {f.w:.5g}, dlam = {st.dlam:.4e})"
    for k, st in enumerate(sts):
        if not st.flags["plastic"]:
            assert st.D == 0.0, f"step {k + 1}: D = {st.D} on an elastic step"
    print(f"\n[cap f {label} amp={amp}] plastic steps checked {len(errs)}, ramp steps {len(ws)}, max rel err {max(errs):.3e}")
    assert len(errs) >= 3, f"only {len(errs)} plastic non-substepped steps"
    assert ws and max(min(w, 1.0 - w) for w in ws) >= 0.05, f"the path does not visit the ramp (ramp w values: {ws[:5]})"
