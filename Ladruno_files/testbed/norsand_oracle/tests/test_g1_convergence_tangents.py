"""WP-144 gate G1 (Zone B): O2 -> O1 convergence, consistent tangents, dissipation census, cap modes.

Plan 144 section 5.1 (gates between the oracles) and 5.2 K1.9; equation sheet 144a sections 4.3, 9, 11, 13.

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form in the
sheet, (b) a published/derived number, or (c) a stated convergence argument, and was written BEFORE the
oracles were run on these paths.  No oracle output is used as an expected value.  A tolerance is never
loosened to make a test pass; a failing assertion is a finding (against the sheet or an oracle), not a
reason to edit an oracle or a bound.

Gate ids used below
  G1.conv     O2(n) -> O1 endpoint convergence, first order (plan 5.1, first bullet of "Gates between them")
  G1.vclosed  specific volume is a closed form of the total strain (sheet 1.2), both oracles
  G1.fd       O2 consistent tangent vs central FD of O2's own stress update, off the theta corners (plan 5.1)
  G1.corner   the same at the exact WW compression corner: O(h) (sheet 4.3); GA control regular
  G1.cto      O2 CTO -> O1 continuum tangent (S.42) as the step -> 0
  G1.elastic  elastic tangent and stress closed forms (S.4, S.5 with alpha0 = 0)
  G1.diss     dissipation census D >= 0 every step, D == 0 on elastic steps (K1.9, S.38/S.39)
  G1.cap      cap modes: smooth completes; planar / no cap near-isotropic stop is REPORTED with the sheet's named status
              (O1) and O2 returns onto the cap / along the axis with pi_i frozen (sheet 3.2, 10.1); closed forms in test_g1_cap
"""
import functools
import math

import numpy as np
import pytest

from conftest import K2_BASE, ORACLES, make_params

O1 = ORACLES["O1"]
O2 = ORACLES["O2"]

SIG0 = -100.0 * np.eye(3)          # isotropic start, p = p0 = -100 kPa, eps^e = 0 (eps_v0 = 0 at p0)
P0 = 100.0
KAPPA = K2_BASE["kappa_hat"]       # 0.01
MU0 = K2_BASE["mu0"]               # 5400
H_HARD = K2_BASE["h"]              # 280
NS = (25, 50, 100, 200, 400)       # O2 increments (plan 5.1)
I3 = np.eye(3)


# ----------------------------------------------------------------------------------------------
# parameter sets, initial states, paths
# ----------------------------------------------------------------------------------------------
def common_params(mode, **over):
    """Parameter set in the sheet's common names.

    paper: AB06 6.1 set with rho = 0.7, rho_bar = 0.8 (K2 case 2); WW.
    fork : the TIMs DM04 mapping of sheet 15 (M = 1.3309, rho = 0.71, N = 0.3, N_bar = 0.2, e0/lambda_c/xi/p_a),
           rho_bar = 0.75 so that the Lode scaling is non-associative too.
    Both satisfy condition A (S.39): N_bar <= N and rho/rho_bar >= (1-N)/(1-N_bar)
    (paper 0.875 >= 0.75; fork 0.9467 >= 0.875).
    """
    kw = dict(K2_BASE)
    if mode == "paper":
        kw.update(rho=0.7, rho_bar=0.8)
    elif mode == "fork":
        kw.update(csl_mode="fork", M=1.3309, N=0.3, N_bar=0.2, rho=0.71, rho_bar=0.75,
                  e0=0.83, lambda_c=0.027, xi=0.45, p_a=101.325)
    else:
        raise KeyError(mode)
    kw.update(over)
    return kw


def v0_for(kw, pi0, psi0):
    """Initial specific volume that gives image state parameter psi0 at pi_i0 (sheet section 6, S.22)."""
    if kw["csl_mode"] == "paper":
        return psi0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-pi0)
    return 1.0 + kw["e0"] + psi0 - kw["lambda_c"] * (-pi0 / kw["p_a"]) ** kw["xi"]


def theta_of(sig):
    """Lode angle from (S.1): y = tr(s^3)/R^3 = cos3theta/sqrt6, theta = acos(sqrt6 y)/3 in [0, pi/3]."""
    p = np.trace(sig) / 3.0
    s = sig - p * I3
    R = np.linalg.norm(s)
    y = np.trace(s @ s @ s) / R ** 3
    return math.acos(max(-1.0, min(1.0, math.sqrt(6.0) * y))) / 3.0


# dense start (psi_i0 = -0.05, pi_i0 = -80: q_yield ~ 91 kPa at p = -100) for the drained paths;
# nearly-critical loose start (psi_i0 = +0.01, pi_i0 = -105) for the undrained and non-coaxial paths.
E1_NC = np.diag([0.003, 0.003, -0.010])                              # triaxial-like, tr = -0.004
# Segment 2 = tensor shear (rotates the principal axes) PLUS continued axial loading.  Pure shear on a coaxial
# yielded state is exactly neutral loading (f : a^e : d_eps = 0 to first order, f is coaxial with sigma), on which
# O1's mode decision chatters (observed: status 'chatter'); the axial part keeps the numerator strictly positive.
E2_NC = np.array([[0.0005, 0.006, 0.0], [0.006, 0.0005, 0.004], [0.0, 0.004, -0.002]])
KNOT = 0.6                                                           # 0.6*25 = 15: a knot for every n in NS

PATHS = {
    "TXC_drained": dict(kind="tx", tx="drained", ax=-0.06, pi0=-80.0, psi0=-0.05, n_truth=20),
    "TXC_undrained": dict(kind="tx", tx="undrained", ax=-0.05, pi0=-105.0, psi0=0.01, n_truth=20),
    "TXE_drained": dict(kind="tx", tx="drained", ax=+0.06, pi0=-80.0, psi0=-0.05, n_truth=20),
    "NONCOAXIAL": dict(kind="gen", pi0=-105.0, psi0=0.01, n_truth=25),
}
MODES = ("paper", "fork")


def noncoax_deps(n):
    """Piecewise-linear total-strain path: E1 over t in [0, KNOT], then E2 over [KNOT, 1]; n equal steps."""
    def eps(s):
        return E1_NC * min(s, KNOT) / KNOT + E2_NC * max(s - KNOT, 0.0) / (1.0 - KNOT)
    t = np.arange(n + 1) / n
    return np.array([eps(t[k + 1]) - eps(t[k]) for k in range(n)])


def run_case(oname, path, mode, n, rtol=1e-10, **over):
    ora = ORACLES[oname]
    spec = PATHS[path]
    kw = common_params(mode, **over)
    P = make_params(oname, **kw)
    v0 = v0_for(kw, spec["pi0"], spec["psi0"])
    st0 = ora.initial_state(P, SIG0, v0, spec["pi0"])
    if spec["kind"] == "tx":
        if oname == "O1":
            sts = ora.triaxial(P, st0, spec["tx"], spec["ax"], n, rtol=rtol)
        else:
            sts = ora.triaxial(P, st0, spec["tx"], spec["ax"], n)
    else:
        deps = noncoax_deps(n)
        sts = ora.run_path(P, st0, deps, rtol=rtol) if oname == "O1" else ora.run_path(P, st0, deps)
    return P, v0, st0, sts


@functools.lru_cache(maxsize=None)
def truth(path, mode):
    """O1 at rtol 1e-10, PATHS[path]['n_truth'] increments (the endpoint does not depend on the chunking)."""
    return run_case("O1", path, mode, PATHS[path]["n_truth"], rtol=1e-10)


@functools.lru_cache(maxsize=None)
def o2run(path, mode, n):
    return run_case("O2", path, mode, n)


def status(oname, sts):
    f = sts[-1].flags
    if oname == "O1":
        return f["status"]
    return f"refused:{f['reason']}" if f["refused"] else "ok"


def complete(oname, sts, n):
    return len(sts) == n and status(oname, sts) == "ok"


def endpoint(oname, path, st):
    """(stress vector or tensor, pi_i, v) at the end of a run.  Triaxial: (sigma_axial, sigma_lateral); the axial
    axis is x in O1 and z in O2 (their drivers differ), so triaxial paths are compared in principal terms."""
    S = st.sigma
    if PATHS[path]["kind"] == "tx":
        ax, lat = (0, 1) if oname == "O1" else (2, 0)
        return np.array([S[ax, ax], S[lat, lat]]), st.pi_i, st.v
    return S.copy(), st.pi_i, st.v


def observed_order(ns, errs):
    """Least-squares slope of log(err) against log(h), h = 1/n (or h itself when ns are step sizes)."""
    x = np.log(1.0 / np.asarray(ns, float))
    y = np.log(np.asarray(errs, float))
    return float(np.polyfit(x, y, 1)[0])


def order_h(hs, errs):
    return float(np.polyfit(np.log(np.asarray(hs, float)), np.log(np.asarray(errs, float)), 1)[0])


def table(names, rows):
    return "\n".join(f"  {nm:>10s}: " + "  ".join(f"{v:.3e}" for v in row) for nm, row in zip(names, rows))


# ----------------------------------------------------------------------------------------------
# G1.conv / G1.vclosed : O2(n) -> O1 convergence (plan 5.1)
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("path", list(PATHS))
def test_g1_o2_converges_to_o1_first_order(path, mode):
    """Plan 5.1: O2 -> O1 with first-order convergence as d_eps -> 0.

    Expected behaviour (stated convergence argument, not oracle output):
      * backward Euler is first order: err(n) ~ C/n, so err strictly decreases along n = 25..400 (h halves each step)
        and the least-squares log-log slope is in [0.8, 1.3];
      * a priori size of the finest error.  BE global error <= (h/2) * (variation of the state)/(strain scale of the
        fastest ODE) with two scales in the model: the elastic scale kappa_hat (0.01, K = -p/kappa_hat) and the
        hardening scale 1/h_hard = 1/280 = 3.6e-3 (pi_i' = h (pi_i* - pi_i) eps_s').  Taking the smaller scale and a
        safety factor 2:  err_rel(n = 400) <= 280 * h_step with h_step the largest strain-increment norm.
      The O1 truth must itself complete (status ok); a truth that stops is a path/finding, not a pass.
    """
    P1, v01, st01, sts1 = truth(path, mode)
    n_truth = PATHS[path]["n_truth"]
    assert complete("O1", sts1, n_truth), f"O1 truth did not complete: {status('O1', sts1)} after {len(sts1)}/{n_truth}"
    sig1, pi1, v1 = endpoint("O1", path, sts1[-1])
    if PATHS[path]["kind"] == "gen":
        # the path must actually rotate the principal axes: a stress state with shear components (the elastic shear
        # stress alone would be 2 mu0 * 0.006 = 65 kPa; require at least 1 kPa of off-diagonal stress survives)
        off = sig1 - np.diag(np.diag(sig1))
        assert np.abs(off).max() > 1.0, f"test bug: final stress is (nearly) coaxial: {sig1}"

    e_sig, e_pi, e_v = [], [], []
    for n in NS:
        P2, v02, st02, sts2 = o2run(path, mode, n)
        assert complete("O2", sts2, n), f"O2 n={n}: {status('O2', sts2)} after {len(sts2)} steps"
        sig2, pi2, v2 = endpoint("O2", path, sts2[-1])
        e_sig.append(np.linalg.norm(sig2 - sig1) / np.linalg.norm(sig1))
        e_pi.append(abs(pi2 - pi1) / abs(pi1))
        e_v.append(abs(v2 - v1) / abs(v1))

    drained = PATHS[path]["kind"] == "tx" and PATHS[path]["tx"] == "drained"
    quantities = [("sigma", e_sig), ("pi_i", e_pi)] + ([("v", e_v)] if drained else [])
    print(f"\n[{path}/{mode}] n = {NS}\n" + table([q[0] for q in quantities], [q[1] for q in quantities]))

    if PATHS[path]["kind"] == "tx":
        h_step = abs(PATHS[path]["ax"]) / NS[-1]
    else:
        h_step = np.linalg.norm(noncoax_deps(NS[-1]), axis=(1, 2)).max()
    bound = H_HARD * h_step

    for name, errs in quantities:
        assert all(e > 0.0 for e in errs), f"{name}: zero error is not a first-order signal: {errs}"
        assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), \
            f"{name}: error not strictly decreasing over n={NS}: {errs}"
        p = observed_order(NS, errs)
        assert 0.8 <= p <= 1.3, f"{name}: observed order {p:.3f} outside [0.8, 1.3]; errs {errs}"
        assert errs[-1] <= bound, f"{name}: err(n=400) = {errs[-1]:.3e} > a priori bound {bound:.3e}"

    if not drained:
        # G1.vclosed: v = v0 (1 + tr eps_total) exactly (sheet 1.2, small strain); strain-controlled paths.
        tr = float(np.trace(E1_NC + E2_NC)) if PATHS[path]["kind"] == "gen" else 0.0
        sts2 = o2run(path, mode, NS[-1])
        for oname, sts, v0, tol in (("O1", sts1, v01, 1e-9), ("O2", sts2[3], sts2[1], 1e-12)):
            v_end = sts[-1].v
            assert abs(v_end - v0 * (1.0 + tr)) <= tol * v0, \
                f"{oname}: v = {v_end}, closed form v0(1+tr eps) = {v0 * (1.0 + tr)}"


# ----------------------------------------------------------------------------------------------
# G1.fd : O2 consistent tangent vs central FD of O2's own stress update, off the corners
# ----------------------------------------------------------------------------------------------
SHEAR_PRE = np.array([[0.0, 2.0e-4, 1.0e-4], [2.0e-4, 0.0, 0.0], [1.0e-4, 0.0, 0.0]])
E_PRE = np.diag([4.0e-4, -1.0e-3, 0.0]) + SHEAR_PRE                  # three distinct principal strains, theta mid-range
DEPS_FD = np.diag([1.0e-4, -6.0e-4, 2.0e-4]) + 0.5 * SHEAR_PRE       # continues the loading direction
NPRE = 10


def _preloaded_o2(mode, **over):
    kw = common_params(mode, **over)
    P = make_params("O2", **kw)
    v0 = v0_for(kw, -80.0, -0.05)
    st0 = O2.initial_state(P, SIG0, v0, -80.0)
    sts = O2.run_path(P, st0, np.array([E_PRE] * NPRE))
    assert not any(s.flags["refused"] for s in sts)
    return P, v0, sts[-1]


def _fd_and_cto(P, st, deps, h):
    """(C, per-column relative error of C against the central FD of sigma(deps), all steps plastic)."""
    stn = O2.step(P, st, deps)
    assert stn.flags["plastic"] and not stn.flags["refused"]
    C = O2.tangent(P, stn)
    errs = {}
    for k in range(3):
        for l in range(k, 3):
            E = np.zeros((3, 3))
            if k == l:
                E[k, k] = 1.0
            else:
                E[k, l] = E[l, k] = 0.5          # dsigma/dh = C_ijkl E_kl = C_ij(kl) by minor symmetry
            sp = O2.step(P, st, deps + h * E)
            sm = O2.step(P, st, deps - h * E)
            assert sp.flags["plastic"] and sm.flags["plastic"], "FD points must stay on the plastic branch"
            fd = (sp.sigma - sm.sigma) / (2.0 * h)
            errs[(k, l)] = np.linalg.norm(C[:, :, k, l] - fd) / np.linalg.norm(C[:, :, k, l])
    return C, stn, errs


@pytest.mark.parametrize("mode", MODES)
def test_g1_o2_cto_matches_central_fd_off_corners(mode):
    """Plan 5.1: 'O2's tangent matches finite differences to ~1e-7 (central differences, several step sizes)';
    the task gate is <= 1e-6 relative, per column, INCLUDING the three shear columns.

    Argument for the step sizes: central-difference truncation ~ (h/kappa_hat)^2/6 relative (kappa_hat = 0.01 is the
    strain scale on which K = -p/kappa_hat changes), i.e. 1.7e-7 at h = 1e-5, 1.7e-9 at h = 1e-6; round-off grows as
    eps_sigma/h.  The gate is the best step size of {1e-5, 1e-6, 1e-7}, max over columns, plus the check that the
    truncation-dominated branch really converges (err(1e-5) < err(1e-4)/10, O(h^2) would give 1/100).
    The state is asserted off both corners (sheet 4.3): theta in [0.1, pi/3 - 0.1] from (S.1).
    """
    P, v0, st = _preloaded_o2(mode)
    stn = O2.step(P, st, DEPS_FD)
    th = theta_of(stn.sigma)
    assert 0.1 < th < math.pi / 3 - 0.1, f"test bug: theta = {th} is not off the corners"
    assert stn.flags["plastic"]
    worst = {}
    for h in (1e-4, 1e-5, 1e-6, 1e-7):
        _, _, errs = _fd_and_cto(P, st, DEPS_FD, h)
        worst[h] = max(errs.values())
        shear = max(errs[(0, 1)], errs[(0, 2)], errs[(1, 2)])
        print(f"\n[{mode}] h={h:.0e}: max column err {worst[h]:.3e}  (shear columns {shear:.3e}; theta {th:.3f})")
    best = min(worst[h] for h in (1e-5, 1e-6, 1e-7))
    assert best <= 1e-6, f"best-h max column error {best:.3e} > 1e-6; all: {worst}"
    assert worst[1e-5] < worst[1e-4] / 10.0, f"FD truncation branch not converging: {worst}"
    # shear columns individually at the best step (the gate names them)
    hb = min((1e-5, 1e-6, 1e-7), key=lambda h: worst[h])
    _, _, errs = _fd_and_cto(P, st, DEPS_FD, hb)
    for kl in ((0, 1), (0, 2), (1, 2)):
        assert errs[kl] <= 1e-6, f"shear column {kl}: {errs[kl]:.3e} at h = {hb}"


# ----------------------------------------------------------------------------------------------
# G1.corner : the exact WW compression corner (sheet 4.3)
# ----------------------------------------------------------------------------------------------
E_AX = np.diag([5.0e-4, 5.0e-4, -2.0e-3])      # axisymmetric (triaxial-compression) loading step
EAX_DIR = np.diag([-0.5, -0.5, 1.0])           # axisymmetric perturbation direction


def _corner_state(**over):
    kw = common_params("paper", **over)
    P = make_params("O2", **kw)
    v0 = v0_for(kw, -80.0, -0.05)
    st0 = O2.initial_state(P, SIG0, v0, -80.0)
    sts = O2.run_path(P, st0, np.array([E_AX] * 6))
    st = sts[-1]
    assert not st.flags["refused"] and st.flags["plastic"]
    assert abs(theta_of(st.sigma) - math.pi / 3) < 1e-6, "test bug: state is not on the TXC corner"
    return P, st


def _corner_errors(P, st, h):
    """(max column relative error over the six columns [symmetry-breaking], relative error of C:E_ax [axisymmetric])."""
    _, stn, errs = _fd_and_cto(P, st, E_AX, h)
    C = O2.tangent(P, stn)
    assert abs(theta_of(stn.sigma) - math.pi / 3) < 1e-6
    sp = O2.step(P, st, E_AX + h * EAX_DIR)
    sm = O2.step(P, st, E_AX - h * EAX_DIR)
    fd = (sp.sigma - sm.sigma) / (2.0 * h)
    ca = np.einsum("ijkl,kl->ij", C, EAX_DIR)
    return max(errs.values()), np.linalg.norm(ca - fd) / np.linalg.norm(ca)


def test_g1_corner_ww_is_first_order_for_symmetry_breaking_and_second_order_axisymmetric():
    """Sheet 4.3: at the exact WW compression corner zeta o sigma is only C^2 in stress, the tangent has a |phi|-kink,
    and a central-difference tangent test converges O(h) for symmetry-breaking perturbations and O(h^2) for
    axisymmetric ones.  Asserted as the ORDER (log-log slope), not as a tight tolerance:
      symmetry-breaking (all six columns, incl. shear): slope in [0.8, 1.3] over h = 1e-4, 1e-5, 1e-6;
      axisymmetric: slope in [1.7, 2.3] over h = 3e-4, 1e-4, 3e-5 (kept above the round-off floor)."""
    P, st = _corner_state()
    hs = (1e-4, 1e-5, 1e-6)
    sb = [_corner_errors(P, st, h)[0] for h in hs]
    p_sb = order_h(hs, sb)
    print(f"\n[corner WW] symmetry-breaking errs at h={hs}: {sb}  order {p_sb:.3f}")
    assert 0.8 <= p_sb <= 1.3, f"symmetry-breaking order {p_sb:.3f}; errs {sb}"
    ha = (3e-4, 1e-4, 3e-5)
    ax = [_corner_errors(P, st, h)[1] for h in ha]
    p_ax = order_h(ha, ax)
    print(f"[corner WW] axisymmetric errs at h={ha}: {ax}  order {p_ax:.3f}")
    assert 1.7 <= p_ax <= 2.3, f"axisymmetric order {p_ax:.3f}; errs {ax}"


def test_g1_corner_ga_control_is_regular():
    """Sheet 4.3: 'GA has no such kink'.  Control: with Gudehus-Argyris (rho 0.8, rho_bar 0.85; condition A holds:
    0.941 >= 0.75) the same corner FD test is O(h^2)-accurate for every column, so the O(h) of WW is the kink of
    zeta, not the FD or the tangent code.  Gate: best h of {1e-5, 1e-6, 1e-7} <= 1e-6 (same as G1.fd)."""
    P, st = _corner_state(zeta="GA", rho=0.8, rho_bar=0.85)
    best = min(_corner_errors(P, st, h)[0] for h in (1e-5, 1e-6, 1e-7))
    print(f"\n[corner GA] best max column err {best:.3e}")
    assert best <= 1e-6


# ----------------------------------------------------------------------------------------------
# G1.cto : O2's CTO -> O1's continuum tangent as the step -> 0
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("mode", MODES)
def test_g1_o2_cto_tends_to_o1_continuum_tangent(mode):
    """Plan 5.1 / sheet 12: the continuum a^ep of (S.42) at state A (O1) is the h -> 0 limit of the consistent
    tangent of a backward-Euler step of size h started at A (O2, same eps^e, pi_i, v).  Not equal at finite h:
    both the state shift and the BE Delta-lambda term change the tangent by O(h).
      * err(h) = |C_O2(h) - C_O1(A)| / |C_O1(A)| is strictly decreasing over h = 1e-4 ... 1e-6 (five values);
      * least-squares order in [0.8, 1.3] (first order);
      * a priori size: the tangent varies on the strain scale kappa_hat (K = -p/kappa_hat, so dC/C ~ d eps / kappa_hat)
        and on 1/h_hard; err(h) <= 10 h / kappa_hat  (safety factor 10) -> 1e-3 at h = 1e-6.
    Loading direction = DEPS_FD normalised to unit Frobenius norm; state A is O1's own after NPRE non-coaxial steps."""
    kw = common_params(mode)
    P1 = make_params("O1", **kw)
    P2 = make_params("O2", **kw)
    v0 = v0_for(kw, -80.0, -0.05)
    st0 = O1.initial_state(P1, SIG0, v0, -80.0)
    sts = O1.run_path(P1, st0, np.array([E_PRE] * NPRE), rtol=1e-10)
    assert len(sts) == NPRE and sts[-1].flags["status"] == "ok", sts[-1].flags["status"]
    A = sts[-1]
    assert A.flags["plastic"], "test bug: state A must be on the plastic branch"
    C1 = O1.tangent(P1, A, plastic_branch=True)
    C1n = np.linalg.norm(C1)

    Ehat = DEPS_FD / np.linalg.norm(DEPS_FD)
    hs = (1e-4, 3e-5, 1e-5, 3e-6, 1e-6)
    errs = []
    for h in hs:
        s2 = O2.initial_state(P2, A.sigma, A.v, A.pi_i)
        s2.v0 = v0                                   # dv/d(eps) = v0 (initial), not the current v (sheet 1.2)
        stn = O2.step(P2, s2, h * Ehat)
        assert stn.flags["plastic"] and not stn.flags["refused"]
        errs.append(np.linalg.norm(O2.tangent(P2, stn) - C1) / C1n)
    print(f"\n[cto->O1 {mode}] h = {hs}\n  err = {errs}")
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), f"not decreasing: {errs}"
    p = order_h(hs, errs)
    assert 0.8 <= p <= 1.3, f"order {p:.3f}; errs {errs}"
    assert errs[-1] <= 10.0 * hs[-1] / KAPPA, f"err(1e-6) = {errs[-1]:.3e} > {10.0 * hs[-1] / KAPPA:.3e}"


# ----------------------------------------------------------------------------------------------
# G1.elastic : closed forms of the elastic tangent and stress (S.4, S.5 with alpha0 = 0)
# ----------------------------------------------------------------------------------------------
def _closed_elastic(eps, p0=-100.0):
    """sigma = p(eps_v) 1 + 2 mu0 dev(eps), p = p0 exp(-eps_v/kappa) (S.5, alpha0 = 0, eps_v0 = 0);
    C = K 1(x)1 + 2 mu0 (I_sym - 1/3 1(x)1), K = -p/kappa."""
    ev = float(np.trace(eps))
    p = p0 * math.exp(-ev / KAPPA)
    sig = p * I3 + 2.0 * MU0 * (eps - ev / 3.0 * I3)
    K = -p / KAPPA
    Isym = 0.5 * (np.einsum("ik,jl->ijkl", I3, I3) + np.einsum("il,jk->ijkl", I3, I3))
    C = K * np.einsum("ij,kl->ijkl", I3, I3) + 2.0 * MU0 * (Isym - np.einsum("ij,kl->ijkl", I3, I3) / 3.0)
    return sig, C


def test_g1_elastic_tangent_and_stress_closed_forms(oracle_name):
    """Sheet 2.2 (S.4, S.5): with alpha0 = 0 the energy gives p = p0 exp(-eps_v/kappa_hat), sigma_dev = 2 mu0 e,
    a^e = K 1(x)1 + 2 mu0 (I - 1/3 1(x)1) with K = -p/kappa_hat.  Checked (1) at the isotropic start, (2) after a small
    non-coaxial ELASTIC step (F < 0 by the margin q_yield ~ 91 kPa >> q ~ 5 kPa).
    Tolerances from the arithmetic: O2's spectral formula (S.33) divides sigma differences (~1e-14 relative round-off,
    ~100 kPa) by strain differences >= 1e-4: error ~ 1e-16*100/(1e-4*5400) ~ 2e-14 relative, so 1e-10; O1 stress comes
    from an ODE at rtol 1e-10, so 1e-9."""
    ora = ORACLES[oracle_name]
    kw = common_params("paper")
    P = make_params(oracle_name, **kw)
    v0 = v0_for(kw, -80.0, -0.05)
    st0 = ora.initial_state(P, SIG0, v0, -80.0)
    tol = 1e-10 if oracle_name == "O2" else 1e-9

    _, C0 = _closed_elastic(np.zeros((3, 3)))
    assert np.linalg.norm(ora.tangent(P, st0) - C0) <= tol * np.linalg.norm(C0)

    deps = np.diag([1.0e-4, -2.0e-4, 5.0e-5]) + np.array([[0, 1e-4, 0], [1e-4, 0, 2e-5], [0, 2e-5, 0]])
    if oracle_name == "O1":
        st = ora.run_path(P, st0, deps[None], rtol=1e-10)[0]
        assert st.flags["status"] == "ok" and not st.flags["plastic"]
    else:
        st = ora.step(P, st0, deps)
        assert not st.flags["plastic"] and not st.flags["refused"]
    sig, C = _closed_elastic(deps)
    assert np.linalg.norm(st.sigma - sig) <= tol * np.linalg.norm(sig), (st.sigma, sig)
    assert np.linalg.norm(ora.tangent(P, st) - C) <= tol * np.linalg.norm(C)
    assert abs(st.pi_i - (-80.0)) == 0.0, "pi_i must not move on an elastic step"


# ----------------------------------------------------------------------------------------------
# G1.diss : dissipation census (K1.9)
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("path", list(PATHS))
def test_g1_dissipation_census_every_step(oracle_name, path, mode):
    """K1.9 / (S.38-S.39): D = Delta-lambda sum_a sigma_a q_a >= 0 at every plastic step and D = 0 at every elastic
    step, for both oracles on every convergence path (O2 at n = 100, O1 at its truth increments), under condition A
    (which both parameter sets satisfy).  O2: exact (D = dlam * positive bracket, dlam >= 0 by KKT).  O1: D is an ODE
    integral of a non-negative rate, so it may undershoot 0 by the integration tolerance only: the absolute tolerance
    on D is 1e-2 * rtol * |p| * eps_max = 1e-12 |p| eps_max, and 1e-9 |p0| eps_max (1000x that) is allowed."""
    if oracle_name == "O1":
        P, v0, st0, sts = truth(path, mode)
        n = PATHS[path]["n_truth"]
        emax = abs(PATHS[path]["ax"]) / n if PATHS[path]["kind"] == "tx" else np.abs(noncoax_deps(n)).max()
        tolD = -1e-9 * P0 * emax
    else:
        P, v0, st0, sts = o2run(path, mode, 100)
        n = 100
        tolD = 0.0
    assert complete(oracle_name, sts, n), f"{oracle_name}: {status(oracle_name, sts)} after {len(sts)}/{n}"
    D = np.array([s.D for s in sts])
    plastic = np.array([bool(s.flags["plastic"]) for s in sts])
    print(f"\n[{oracle_name} {path}/{mode}] steps {n}, plastic {int(plastic.sum())}, min D {D.min():.3e}, "
          f"min D plastic {D[plastic].min() if plastic.any() else float('nan'):.3e}")
    assert plastic.any(), "the path never yields: not a dissipation test"
    assert np.all(D >= tolD), f"negative dissipation: min D = {D.min():.3e}, steps {np.where(D < tolD)[0]}"
    assert np.all(D[~plastic] == 0.0), "D must be exactly 0 on elastic steps"


# ----------------------------------------------------------------------------------------------
# G1.cap : cap modes on a near-isotropic compression path
# ----------------------------------------------------------------------------------------------
def e_iso_total(amp):
    """Near-isotropic compression: -0.01 on every axis (tr = -0.03, well beyond the hydrostatic yield at
    p = pi_c = pi_i/(1-N)^((1-N)/N) = -172 kPa) plus a small deviator amp*diag(1, 0, -1)."""
    return -0.01 * I3 + amp * np.diag([1.0, 0.0, -1.0])


AMP_STOP = 2.0e-3            # deviator/volumetric ratio 0.2/3: near-isotropic, drives q/|p| to ~0.1 at the end
N_CAP = 40
CAP_KW = {
    "smooth": dict(cap="smooth", c1=0.05, c2=0.15),
    "planar": dict(cap="planar", c1=0.10, c2=0.10),     # BA06's chi_cap = 0.10 (sheet 10.1)
    "none": dict(cap="none"),
}
VERTEX_FAMILY = ("vertex_reached", "vertex_nonisotropic")


def _run_cap(oname, cap, amp, n=N_CAP):
    ora = ORACLES[oname]
    kw = common_params("paper", **CAP_KW[cap])
    P = make_params(oname, **kw)
    v0 = v0_for(kw, -80.0, -0.05)
    st0 = ora.initial_state(P, SIG0, v0, -80.0)
    deps = np.array([e_iso_total(amp) / n] * n)
    sts = ora.run_path(P, st0, deps, rtol=1e-10) if oname == "O1" else ora.run_path(P, st0, deps)
    return P, sts


@pytest.mark.parametrize("amp", [1.0e-4, AMP_STOP], ids=["dev1e-4", "dev2e-3"])
def test_g1_cap_smooth_runs_to_completion(oracle_name, amp):
    """Sheet 10.2: the smooth cap (quintic blend, c1 = 0.05, c2 = 0.15) removes the corner, so the near-isotropic
    compression path must run to completion in both oracles (O1: every increment status ok; O2: no refusal), with
    D >= 0 at every step (K1.9 with the cap: D = lambda[-p(1-w) + w D_u] >= 0, sheet 10.2).  Two deviator sizes: the
    weakly non-isotropic 1e-4 and the AMP_STOP path on which the planar cap and no cap stop (below)."""
    P, sts = _run_cap(oracle_name, "smooth", amp)
    assert complete(oracle_name, sts, N_CAP), f"{oracle_name}: {status(oracle_name, sts)} after {len(sts)}/{N_CAP}"
    D = np.array([s.D for s in sts])
    assert any(s.flags["plastic"] for s in sts), "the path never yields"
    tol = -1e-9 * P0 * np.abs(e_iso_total(amp) / N_CAP).max() if oracle_name == "O1" else 0.0
    assert np.all(D >= tol), f"min D {D.min():.3e}"


def test_g1_cap_smooth_o2_completes_and_converges_under_refinement():
    """Sheet 10.2 + plan 5.1 (convergence argument, not oracle output).  The strict gate above is N_CAP = 40 increments.
    If O2 refuses there, the next question is whether the refusal is a step-size artefact of backward Euler or a
    failure of the return map at a state the continuum path reaches.  A backward-Euler map that is consistent with the
    rate problem must, for a fine enough step, complete the same path and converge to the O1 truth at first order.
    Gate: O2 completes at n = 320 and n = 640 increments of the AMP_STOP smooth-cap path; endpoint stress and pi_i errors
    against O1 (rtol 1e-10, endpoint independent of the chunking) decrease and the two-point observed order is in
    [0.8, 1.3].  (Failure here = the smooth cap is not usable even on a fine grid: a sheet/O2 finding.)"""
    P1, sts1 = _run_cap("O1", "smooth", AMP_STOP)
    assert complete("O1", sts1, N_CAP), f"O1 truth did not complete: {status('O1', sts1)}"
    S1, pi1 = sts1[-1].sigma, sts1[-1].pi_i
    ns, e_sig, e_pi = (320, 640), [], []
    for n in ns:
        P2, sts2 = _run_cap("O2", "smooth", AMP_STOP, n=n)
        assert complete("O2", sts2, n), f"O2 n={n}: {status('O2', sts2)} after {len(sts2)}/{n}"
        e_sig.append(np.linalg.norm(sts2[-1].sigma - S1) / np.linalg.norm(S1))
        e_pi.append(abs(sts2[-1].pi_i - pi1) / abs(pi1))
    print(f"\n[smooth cap refinement] n = {ns}\n{table(['sigma', 'pi_i'], [e_sig, e_pi])}")
    for name, errs in (("sigma", e_sig), ("pi_i", e_pi)):
        assert errs[1] < errs[0], f"{name}: error not decreasing: {errs}"
        p = observed_order(ns, errs)
        assert 0.8 <= p <= 1.3, f"{name}: observed order {p:.3f}; errs {errs}"


@pytest.mark.parametrize("cap", ["planar", "none"])
def test_g1_cap_planar_and_nocap_stops_are_reported_not_hidden(oracle_name, cap):
    """Sheet 3.2 (vertex) and 10.1 (planar cap, Filippov sliding), as revised at G1.
    O1 (continuum rate oracle) must STOP, with the named status, strictly before the end of the path, and the stop state
    must be the one the sheet names:
      none  : status 'vertex_reached'.  Near-isotropic compression drives q -> 0 in finite time and the vertex has no flow
              direction: the stop state has q <= 1e-6 |p| (the oracle's event is R/|p| = 10 R_tol = 1e-7 with the
              R_tol = 1e-8 |p| of sheet 3.2; q = sqrt(3/2) R; margin 10x).
      planar: status 'cap_sliding'.  The state is attracted to the Q-corner eta = c1 M = 0.12 (sheet 10.1) and the stop
              state lies ON it: |eta(p, pi_i) - c1 M| <= 1e-6 (eta recomputed here from S.12), and pi_i is frozen at its
              initial value on EVERY state of the run (Omega = 0 on the capped side, S.25; ODE tolerance 1e-10 -> 1e-9 |pi_i|).
    O2 (backward Euler) has no rate problem: at every plastic, non-substepped step before a refusal the plastic strain
    increment is Delta lambda q_a with q_a the closed form of (S.17) / (S.35) / the vertex rule recomputed here
    (test_g1_cap.check_o2_step), to 1e-8; where Omega = 0 (planar cap: eta < c1 M; or R < R_tol) pi_i and eps^p_s are exactly
    frozen, and for the planar cap at least one step is on the cap (q_a = -delta_a/3).  If O2 refuses, the refusal carries the
    reason the sheet documents (local_linesearch after 2^8 substeps) and every remaining step repeats the refused state; an
    O2 that completes the path is also accepted (the closed forms were checked at every step).  The axial return to the
    apex with pi_i frozen is tested exactly in test_g1_o2_hydrostatic_apex_returns_along_the_axis_with_pi_frozen.
    Mutants killed: O1 completing 'ok' on these paths (the stop hidden) or stopping with another status; O1 stopping away
    from the corner; O1 hardening pi_i on the capped side; O2 returning a state whose plastic flow is not q_a; O2 evolving
    pi_i where Omega = 0; an O2 refusal with an empty or different reason."""
    P, sts = _run_cap(oracle_name, cap, AMP_STOP)
    st = status(oracle_name, sts)
    if oracle_name == "O1":
        want = "cap_sliding" if cap == "planar" else "vertex_reached"
        print(f"\n[O1 {cap}] status {st} after {len(sts)}/{N_CAP} increments")
        assert st != "ok", f"O1 completed silently on the {cap} near-isotropic path (expected an explicit stop)"
        assert st == want, f"O1 status {st!r}, the sheet names {want!r}"
        assert len(sts) < N_CAP, "a non-ok status must end the run (run_path breaks at the first non-ok increment)"
        last = sts[-1]
        lam = np.linalg.eigvalsh(0.5 * (last.sigma + last.sigma.T))
        p = float(lam.mean())
        if cap == "none":
            q = math.sqrt(1.5) * float(np.linalg.norm(lam - p))
            assert q <= 1e-6 * abs(p), f"vertex_reached with q/|p| = {q / abs(p):.3e}"
        else:
            M, N = K2_BASE["M"], K2_BASE["N"]
            eta = (M / N) * (1.0 - (1.0 - N) * (p / last.pi_i) ** (N / (1.0 - N)))
            assert abs(eta - CAP_KW["planar"]["c1"] * M) <= 1e-6, f"cap_sliding off the corner: eta = {eta}"
            assert all(abs(s.pi_i - (-80.0)) <= 1e-9 * 80.0 for s in sts), \
                f"pi_i moved on the capped side: {[s.pi_i for s in sts]}"
        return
    # O2: the closed forms live in test_g1_cap (imported here, not at module level, to keep the two files decoupled)
    from test_g1_cap import check_o2_step, iter_plastic, kw_for, refused_at
    kw = kw_for(cap)
    st0 = O2.initial_state(P, SIG0, v0_for(kw, -80.0, -0.05), -80.0)
    print(f"\n[O2 {cap}] status {st} after {len(sts)}/{N_CAP} steps")
    k = refused_at(sts)
    if k is not None:
        assert all(s.flags["refused"] for s in sts[k:]), "after a refusal every remaining step must repeat the refused state"
        reason = sts[k].flags["reason"]
        assert reason and "local_linesearch" in reason and "substeps exhausted at 2^8" in reason, \
            f"refusal reason {reason!r} is not the documented local_linesearch after 2^8 substeps"
    n_pl, n_on_cap = 0, 0
    for kk, prev, s in iter_plastic(st0, sts):
        f, _ = check_o2_step(kw, prev, s, rel=1e-8)
        n_pl += 1
        n_on_cap += f.w == 0.0
    print(f"[O2 {cap}] plastic steps checked {n_pl}, on the cap {n_on_cap}")
    assert n_pl >= 3
    if cap == "planar":
        assert n_on_cap >= 1, "the planar path never reached the cap"
    assert all(s.D >= 0.0 for s in sts if not s.flags["refused"]), "negative dissipation"


@pytest.mark.parametrize("cap", ["none", "planar", "smooth"])
def test_g1_o2_hydrostatic_apex_returns_along_the_axis_with_pi_frozen(cap):
    """Sheet 3.2 (vertex rule, 'what the two integrators do there') and 10.1: a hydrostatic trial state beyond the
    compression apex pi_c is returned ALONG THE AXIS to p = pi_c with pi_i frozen (hydrostatic compression beyond pi_c is
    perfectly plastic: no volumetric hardening in this model).  Start ON the apex: sigma = p 1 with p = -100 kPa, and
    pi_i = pi_c (1-N)^((1-N)/N) (sheet 5.1: pi_c = pi_i/(1-N)^((1-N)/N)), i.e. -100 * 0.6^1.5 = -46.4758 kPa; then three
    hydrostatic steps of tr(Delta eps) = -3e-3.  Closed forms, for the three cap modes alike (at eta = 0 the cap weight is
    w = 0 for the planar and the smooth cap, and the vertex rule gives a purely volumetric flow without a cap):
      * p = pi_c = -100 after every step  (F = p eta = 0 <=> eta = 0 <=> p = pi_c): 1e-9 |p| (Newton residual 1e-12);
      * pi_i unchanged, EXACTLY; eps^p_s unchanged, EXACTLY (Omega = 0);
      * Delta eps^p_v = tr(Delta eps) (the elastic strain returns to its apex value): 1e-8 |tr Delta eps|, and
        Delta eps^p is isotropic (deviatoric part <= 1e-12 ||Delta eps^p||);
      * D = p Delta eps^p_v = -p |tr Delta eps| = 0.300 kPa per step (sigma:Delta eps^p with sigma = p 1; in the no-cap case
        Delta eps^p_v = Delta lambda beta F_p and D = Delta lambda p beta F_p, in the cap case Delta eps^p_v = -Delta lambda and
        D = -p Delta lambda): 1e-8 relative.
    Mutants killed: Omega not zeroed at the vertex (pi_i would harden, eps^p_s would grow); the axis return not along the
    axis (p drifts off pi_c); a vertex rule that returns elastically or with a deviatoric flow; a cap weight that is not 0 at
    eta = 0 (the planar flow would not be the axial one)."""
    kw = common_params("paper", **CAP_KW[cap])
    P = make_params("O2", **kw)
    N = kw["N"]
    p_apex = -100.0
    pi0 = p_apex * (1.0 - N) ** ((1.0 - N) / N)
    st = O2.initial_state(P, p_apex * I3, 1.6, pi0)
    assert abs(st.pi_i - pi0) == 0.0
    dtr = -3.0e-3
    for step_no in range(1, 4):
        prev = st
        st = O2.step(P, prev, (dtr / 3.0) * I3)
        assert not st.flags["refused"], st.flags["reason"]
        assert st.flags["plastic"], f"step {step_no} did not yield (hydrostatic compression beyond the apex)"
        p = float(np.trace(st.sigma)) / 3.0
        assert abs(p - p_apex) <= 1e-9 * abs(p_apex), f"step {step_no}: p = {p} (apex {p_apex})"
        assert np.linalg.norm(st.sigma - p * I3) <= 1e-9 * abs(p_apex), "stress left the hydrostatic axis"
        assert st.pi_i == prev.pi_i, f"step {step_no}: pi_i moved {prev.pi_i} -> {st.pi_i}"
        assert st.eps_p_s == prev.eps_p_s, f"step {step_no}: eps^p_s grew"
        dep = st.eps_p - prev.eps_p
        dpv = float(np.trace(dep))
        assert abs(dpv - dtr) <= 1e-8 * abs(dtr), f"step {step_no}: d eps^p_v = {dpv}, tr d eps = {dtr}"
        assert np.linalg.norm(dep - dpv / 3.0 * I3) <= 1e-12 * np.linalg.norm(dep)
        assert abs(st.D - (-p_apex * abs(dtr))) <= 1e-8 * (-p_apex * abs(dtr)), f"step {step_no}: D = {st.D}"
