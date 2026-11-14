"""O2 self-checks (sheet 144a). Run:  python -m o2_algo.selfcheck   from norsand_oracle/.
Prints numbers; it is NOT the gate test suite (P0d writes that)."""
from __future__ import annotations

import math
import warnings

import numpy as np

from . import kernel as K
from .acoustic import acoustic_min_det, acoustic_principal, acoustic_tensor
from .api import (State, floor_energy, initial_state, k2_path, run_path, step, step_fractions, tangent, tangent_finite,
                  tangent_last_substep, triaxial)
from .params import Params

np.set_printoptions(precision=6, linewidth=140)
I3 = np.eye(3)


def k2_params(rho=0.7, rho_bar=0.8, **kw):
    d = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0, M=1.2, N=0.4, N_bar=0.2,
             rho=rho, rho_bar=rho_bar, zeta="WW", chi=-3.5, h=280.0, csl_mode="paper", lam_tilde=0.0135,
             v_c0=1.81, cap="none")
    d.update(kw)
    return Params(**d).validate()


def fork_params(**kw):
    d = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0, M=1.3309, N=0.3, N_bar=0.2,
             rho=0.71, rho_bar=0.71, zeta="WW", chi=-3.5, h=280.0, csl_mode="fork", e0=0.83, lam_c=0.027,
             xi=0.45, p_a=101.325, cap="none")
    d.update(kw)
    return Params(**d).validate()


def tims_params(p_a=101.0, **kw):
    """HAR TIMs set (sheet §2.3, §13.14b): n = 1/2, g = G0 f(e_ref), k from nu; fork CSL at the TIMs values, M 1.3309,
    rho = rho_bar 0.71 (WW), and the K2 plastic constants N 0.4, N_bar 0.2, chi -3.5, h 280 (TIMs' own are a P3 refit).
    p_a = 101 kPa, the campaign's `Patm 101` (round 3b, A4: 101.325 is WITHDRAWN as the TIMs value; it stays only as
    O2's inactive Params default). The sheet's K1.1h/K1.11/K1.13/K1.14/K1.14b values are at 101."""
    d = dict(energy="HAR", k=1889.48104361, g=807.80387674, n_e=0.5, p_a=p_a, M=1.3309, N=0.4, N_bar=0.2,
             rho=0.71, rho_bar=0.71, zeta="WW", chi=-3.5, h=280.0, csl_mode="fork", e0=0.83, lam_c=0.027, xi=0.45,
             cap="none")
    d.update(kw)
    return Params(**d).validate()


def sym(a):
    return 0.5 * (a + a.T)


def off_corner_sig(P, p_s, eta_s, direction=(-1.0, -0.35, 1.35)):
    """principal stress at (p, eta) on F = 0 with a deviatoric direction away from both WW corners."""
    xi = np.array(direction, float)
    xi -= xi.mean()
    nh = xi / np.linalg.norm(xi)
    inv = K.invariants(p_s + nh)
    z, _, _ = K.zeta_y(inv.theta, P.rho, P.zeta)
    q_s = eta_s * abs(p_s) / z
    return p_s + K.SQ23 * q_s * nh, inv.theta, nh


def v_for_psi(P, pi, psi_target):
    """v with psi_i = e - e_c(pi) = psi_target (fork) or v - v_c0 + lam ln(-pi) = psi_target (paper)."""
    if P.csl_mode == "fork":
        return 1.0 + P.e0 - P.lam_c * (-pi / P.p_a) ** P.xi + psi_target
    return psi_target + P.v_c0 - P.lam_tilde * math.log(-pi)


def relmax(A, B):
    return float(np.abs(A - B).max() / np.abs(B).max())


def ok(name, err, tol):
    print(f"  [{'OK' if err <= tol else 'FAIL'}] {name}: {err:.3e} (tol {tol:.0e})")
    return err <= tol


# ------------------------------------------------------------------ 1. Jacobian vs FD
def jacobian_fd(P, st, deps, h=1e-6):
    """returns (rel err, local iters) at the converged plastic state of the step st -> st+deps."""
    tr = float(np.trace(deps))
    v = st.v * math.exp(tr)                     # (S.26) exponential update, vfac = v (G2)
    eps_tr = st.eps_e + sym(deps)
    w, V = np.linalg.eigh(eps_tr)
    res = K.return_map(P, w, st.pi_i, v, v)
    assert res.plastic and not res.refused, res.reason
    x = np.append(res.eps_e, res.dlam)
    pe = K.evaluate(P, x[:3], x[3], w, v, st.pi_i)
    J = K.jacobian(P, pe)
    Jfd = np.empty((4, 4))
    for j in range(4):
        e = np.zeros(4)
        e[j] = h
        rp = K.evaluate(P, (x + e)[:3], (x + e)[3], w, v, st.pi_i).r
        rm = K.evaluate(P, (x - e)[:3], (x - e)[3], w, v, st.pi_i).r
        Jfd[:, j] = (rp - rm) / (2 * h)
    S = np.diag([1, 1, 1, 1 / P.p_ref])
    return np.linalg.norm(S @ (J - Jfd)) / np.linalg.norm(S @ J), res


def check_jacobian():
    print("\n=== 1. 4x4 Jacobian (S.30) vs central FD of the residual (S.29) ===")
    shear = np.array([[0, 2e-4, 1e-4], [2e-4, 0, 0], [1e-4, 0, 0]])
    cases = []
    for name, P in (("paper", k2_params()), ("fork", fork_params()), ("fork-cap", fork_params(cap="smooth", c1=0.3, c2=0.6))):
        s0 = initial_state(P, -100 * I3, 1.59 if P.csl_mode == "paper" else 1.65, -60.4 if "cap" not in name else -80.0)
        path = [np.diag([4e-4, -1e-3, 0.0]) + shear] * 6 + [np.diag([1e-4, -4e-4, 1e-3]) + 2 * shear] * 6
        sts = run_path(P, s0, np.array(path))
        for k in (4, 8, 11):
            errs = [jacobian_fd(P, sts[k - 1], path[k], h)[0] for h in (1e-6, 1e-7, 1e-8)]
            res = jacobian_fd(P, sts[k - 1], path[k])[1]
            cases.append(errs[2])
            print(f"  {name:9s} state {k:2d}: rel err h=1e-6/1e-7/1e-8: {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e}  "
                  f"(dlam {res.dlam:.3e}, iters {res.local_iters}, theta {K.invariants(res.sig).theta:.4f}, "
                  f"cap_active {res.cap_active})")
    print(f"  max rel err at h=1e-8 = {max(cases):.2e}  (target <= 1e-7; the h=1e-6 values are O(h^2) FD truncation "
          "of r4 in dlam, J_44 ~ 1e5)")


# ------------------------------------------------------------------ 2. CTO vs FD
def cto_fd(P, st, deps, h=1e-6, finite=False):
    stn = step(P, st, deps)
    C = tangent_finite(P, stn) if finite else tangent(P, stn)
    Cfd = np.zeros((3, 3, 3, 3))
    for k in range(3):
        for l in range(k, 3):
            E = np.zeros((3, 3))
            E[k, l] = E[l, k] = 0.5 if k != l else 1.0
            sp = step(P, st, deps + h * E).sigma
            sm = step(P, st, deps - h * E).sigma
            d = (sp - sm) / (2 * h)
            Cfd[:, :, k, l] = d
            Cfd[:, :, l, k] = d
    return np.linalg.norm(C - Cfd) / np.linalg.norm(C), stn, C


def check_cto():
    print("\n=== 2. consistent tangent (S.33) vs central FD of the converged stress ===")
    shear = np.array([[0, 3e-4, 1e-4], [3e-4, 0, 2e-4], [1e-4, 2e-4, 0]])
    for name, P in (("paper WW", k2_params()), ("fork WW", fork_params()), ("paper GA", k2_params(zeta="GA", rho=0.8, rho_bar=0.85)),
                    ("fork cap", fork_params(cap="smooth", c1=0.3, c2=0.6))):
        s0 = initial_state(P, -100 * I3, 1.59 if P.csl_mode == "paper" else 1.65, -60.4 if "cap" not in name else -80.0)
        pre = [np.diag([4e-4, -1e-3, 0.0]) + shear] * 5
        sts = run_path(P, s0, np.array(pre))
        st = sts[-1]
        for h in (1e-5, 1e-6):
            d = np.diag([1e-4, -6e-4, 2e-4]) + 0.5 * shear
            err, stn, _ = cto_fd(P, st, d, h)
            th = K.invariants(np.linalg.eigvalsh(stn.sigma)).theta
            print(f"  {name:9s} h={h:.0e}: rel err {err:.2e}  (plastic {stn.flags['plastic']}, theta {th:.4f}, "
                  f"cap_active {stn.flags['cap_active']}, sym err {np.linalg.norm(_.reshape(9,9)-_.reshape(9,9).T)/np.linalg.norm(_):.3f})")
    # elastic step
    P = k2_params()
    s0 = initial_state(P, -100 * I3, 1.59, -60.4)
    err, stn, _ = cto_fd(P, s0, np.diag([1e-4, -2e-4, 0.5e-4]) + shear, 1e-6)
    print(f"  elastic non-coaxial step (repeated-eigenvalue start): rel err {err:.2e}, plastic {stn.flags['plastic']}")
    # exact TXC corner
    P = k2_params()
    s0 = initial_state(P, -100 * I3, 1.59, -60.4)
    pre = [np.diag([5e-4, 5e-4, -2e-3])] * 5
    st = run_path(P, s0, np.array(pre))[-1]
    for h in (1e-4, 1e-5, 1e-6):
        err, stn, _ = cto_fd(P, st, np.diag([5e-4, 5e-4, -2e-3]), h)
        print(f"  TXC corner (theta=pi/3, WW) h={h:.0e}: rel err {err:.2e}   [O(h) expected, sheet §4.3]")
    P = k2_params(zeta="GA", rho=0.8, rho_bar=0.85)
    st = run_path(P, initial_state(P, -100 * I3, 1.59, -60.4), np.array(pre))[-1]
    err, stn, _ = cto_fd(P, st, np.diag([5e-4, 5e-4, -2e-3]), 1e-6)
    print(f"  TXC corner (GA) h=1e-6: rel err {err:.2e}   [regular, no kink]")


# ------------------------------------------------------------------ 3. quadratic convergence
def check_quadratic():
    print("\n=== 3. local Newton convergence on a hard step ===")
    P = k2_params()
    s0 = initial_state(P, -100 * I3, 1.59, -60.4)
    st = run_path(P, s0, np.array([np.diag([4e-4, -1e-3, 0.0])] * 3))[-1]
    big = np.diag([2e-3, -8e-3, 1e-3]) + np.array([[0, 2e-3, 0], [2e-3, 0, 1e-3], [0, 1e-3, 0]])
    tr = float(np.trace(big))
    w, V = np.linalg.eigh(st.eps_e + big)
    v = st.v * math.exp(tr)                     # (S.26) exponential update, vfac = v (G2)
    res = K.return_map(P, w, st.pi_i, v, v)
    print("  scaled residual history:", " ".join(f"{r:.2e}" for r in res.res_hist), "| refused:", res.refused, res.reason)
    print(f"  dlam {res.dlam:.4e}, iters {res.local_iters}, nested pi iters total {res.pi_iters}")


# ------------------------------------------------------------------ 4. K1
def check_k1():
    print("\n=== 4. K1 closed forms (sheet §13) ===")
    P = k2_params()
    # K1.1
    el = K.elastic(P, np.array([-0.01 / 3] * 3))
    print(f"  K1.1 p(eps_v=-0.01) = {el.p:.6f}  (closed form {P.p0 * math.e:.6f})")
    # K1.2 closed loop
    s0 = initial_state(P, -100 * I3, 1.59, -2000.0)   # far inside the surface
    verts = [np.zeros((3, 3)), np.diag([1e-3, -2e-3, 5e-4]) + np.array([[0, 5e-4, 0], [5e-4, 0, 2e-4], [0, 2e-4, 0]]),
             np.diag([-1e-3, -1e-3, -2e-3]) + np.array([[0, 0, 4e-4], [0, 0, 0], [4e-4, 0, 0]]),
             np.diag([5e-4, 0, -1e-3]) + np.array([[0, -3e-4, 0], [-3e-4, 0, 0], [0, 0, 0]]), np.zeros((3, 3))]
    xg, wg = np.polynomial.legendre.leggauss(8)
    W, Wabs = 0.0, 0.0
    for a, b in zip(verts[:-1], verts[1:]):
        d = b - a
        for xi_, wi in zip(xg, wg):
            e = a + 0.5 * (xi_ + 1) * d
            w, V = np.linalg.eigh(s0.eps_e + e)
            sig = (V * K.elastic(P, w).sig) @ V.T
            W += 0.5 * wi * float(np.sum(sig * d))
            Wabs += 0.5 * wi * abs(float(np.sum(sig * d)))
    sts = run_path(P, s0, np.array([b - a for a, b in zip(verts[:-1], verts[1:])]))
    print(f"  K1.2 closed non-coaxial elastic loop: W/sum|W| = {W / Wabs:.2e}; all elastic {all(not s.flags['plastic'] for s in sts)}; "
          f"state return |d sigma| {np.abs(sts[-1].sigma - s0.sigma).max():.2e}, |d eps_e| {np.abs(sts[-1].eps_e - s0.eps_e).max():.2e}")
    # K1.3
    for kind, rho in (("WW", 0.7), ("WW", 0.5), ("GA", 7 / 9), ("GA", 0.9)):
        z0, z10, _ = K.zeta_theta(0.0, rho, kind)
        z1, z11, _ = K.zeta_theta(math.pi / 3, rho, kind)
        print(f"  K1.3 {kind} rho={rho:.4f}: zeta(0)-1/rho = {z0 - 1 / rho:.1e}, zeta(pi/3)-1 = {z1 - 1:.1e}, zeta'(0) = {z10:.1e}, zeta'(pi/3) = {z11:.1e}")
    for bad in (dict(zeta="GA", rho=0.7, rho_bar=0.8), dict(rho=0.45, rho_bar=0.5), dict(rho=0.5, rho_bar=0.5),
                dict(rho=0.6, rho_bar=0.5), dict(N_bar=0.5), dict(N=0.4, N_bar=0.4, rho=0.7, rho_bar=0.8)):
        try:
            k2_params(**bad)
            print("  K1.3/K1.9 refusal MISSED for", bad)
        except ValueError as e:
            print(f"  K1.3/K1.9 refused {bad}: {str(e)[:60]}")
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        k2_params(rho=0.8, rho_bar=0.7)
        print(f"  K1.9 rho > rho_bar warns: {len(wl) == 1 and 'rho' in str(wl[0].message)}")
    # K1.4 image point
    pi = -150.0
    for th_name, sig in (("TXC", np.array([pi, pi, pi]) + np.array([1, 1, -2]) * P.M * abs(pi) / 3),
                         ("TXE", np.array([pi, pi, pi]) + np.array([-1, -1, 2]) * P.rho * P.M * abs(pi) / 3)):
        inv = K.invariants(sig)
        print(f"  K1.4 {th_name}: p = {inv.p:.6f} (= pi_i {pi}), theta = {inv.theta:.6f}, F = {K.flow(P, inv, pi).F:.2e}")
    # K1.5, K1.6, K1.8, K1.9 on a drained triaxial (fork, dense)
    for name, P2, v0 in (("fork", fork_params(), 1.70), ("paper", k2_params(), 1.72)):
        s0 = initial_state(P2, -100 * I3, v0, None)
        sts = triaxial(P2, s0, "drained", -0.25, 500)
        flow_err, Dmin, Hprev, peak = 0.0, np.inf, None, None
        for i, s in enumerate(sts):
            if s.flags["plastic"] and not s.flags["vertex"]:
                dv = s.eps_p_v - (sts[i - 1].eps_p_v if i else 0.0)
                ds = s.eps_p_s - (sts[i - 1].eps_p_s if i else 0.0)
                flow_err = max(flow_err, abs(dv / ds - math.sqrt(1.5) * P2.beta * s.F_p / s.Omega))
                Dmin = min(Dmin, s.D)
                w, _ = np.linalg.eigh(s.sigma)
                inv = K.invariants(w)
                H = -P2.M * (inv.p / s.pi_i) ** (1 / (1 - P2.N)) * K.SQ23 * P2.h * (s.pi_star - s.pi_i) * s.Omega
                Dil = math.sqrt(1.5) * P2.beta * s.F_p / s.Omega
                if Hprev is not None and Hprev * H < 0 and peak is None:
                    peak = (i + 1, Dil - P2.chi * s.psi_i, abs(s.pi_star - s.pi_i), s.eta)
                Hprev = H
            else:
                Dmin = min(Dmin, s.D)
        last = sts[-1]
        w, _ = np.linalg.eigh(last.sigma)
        inv = K.invariants(w)
        q_over_p = -inv.q / inv.p
        refused = any(s.flags["refused"] for s in sts)
        print(f"  drained TXC [{name}] {len(sts)} steps, refused {refused}: K1.5 max|flow-rule residual| = {flow_err:.2e}; "
              f"K1.9 min D = {Dmin:.3e}")
        if peak:
            print(f"     K1.6 H sign change at step {peak[0]}: D - chi psi_i = {peak[1]:+.3e}, |pi* - pi_i| = {peak[2]:.2e}, eta = {peak[3]:.4f}")
        print(f"     K1.8 at eps_a = -0.25: eta = {last.eta:.5f} (M = {P2.M}), q/|p| = {q_over_p:.5f}, psi_i = {last.psi_i:+.2e}, "
              f"|pi*-pi_i|/|pi_i| = {abs(last.pi_star - last.pi_i) / abs(last.pi_i):.2e}")
    # K1.7 undrained CS
    for name, P2, v0 in (("fork", fork_params(), 1.78), ("paper", k2_params(), 1.72)):
        s0 = initial_state(P2, -100 * I3, v0, None)
        sts = triaxial(P2, s0, "undrained", -2.0, 4000)   # asymptotic approach: -4.0/8000 gives 1e-8 (fork)
        last = sts[-1]
        w, _ = np.linalg.eigh(last.sigma)
        inv = K.invariants(w)
        if P2.csl_mode == "fork":
            e = v0 - 1
            p_cs = -P2.p_a * ((P2.e0 - e) / P2.lam_c) ** (1 / P2.xi)
        else:
            p_cs = -math.exp((P2.v_c0 - v0) / P2.lam_tilde)
        refused = any(s.flags["refused"] for s in sts)
        print(f"  K1.7 undrained [{name}] at eps_a = -2.0, refused {refused}: p = {inv.p:.4f} vs p_cs = {p_cs:.4f} (rel {abs(inv.p - p_cs) / abs(p_cs):.2e}); "
              f"q/|p| = {-inv.q / inv.p:.5f} (M {P2.M}); pi_i = {last.pi_i:.4f}; psi_i = {last.psi_i:+.2e}; v const {abs(last.v - v0):.1e}; "
              f"min D {min(s.D for s in sts):.2e}")


# ------------------------------------------------------------------ 5. finite-strain pieces + K2
def check_finite_and_k2(full_table=False):
    print("\n=== 5. finite-strain tangent pieces and K2 (sheet §14) ===")
    P = k2_params()
    s0 = initial_state(P, -100 * I3, 1.59, -60.4, finite=True)
    f1 = np.log(np.array([1 + 4e-4, 1 - 1e-3, 1.0]))
    sts = run_path(P, s0, np.array([np.diag(f1)] * 12))
    st = sts[-1]
    # a~^ep (finite, v-factor) vs FD of tau_a w.r.t. log-stretch increments (diagonal)
    at = st.cache["atilde"]
    Vc = st.cache["nvec"]                      # eigh ordering -> permutation of the diagonal axes
    perm = [int(np.argmax(np.abs(Vc[:, a]))) for a in range(3)]
    h = 1e-6
    at_fd = np.zeros((3, 3))
    for b in range(3):
        e = np.zeros(3)
        e[perm[b]] = h
        sp = np.diag(step(P, sts[-2], np.diag(f1 + e)).sigma)[perm]
        sm = np.diag(step(P, sts[-2], np.diag(f1 - e)).sigma)[perm]
        at_fd[:, b] = (sp - sm) / (2 * h)
    print(f"  a~^ep (S.32, vfac = v) vs FD, finite mode, plastic: rel err {np.linalg.norm(at - at_fd) / np.linalg.norm(at):.2e}")
    # (S.44) transcription vs generic contraction of (S.34)
    a4 = tangent_finite(P, st)
    tau = st.cache["sig"]
    lam = np.exp(st.cache["eps_tr"])
    c_ab = at - 2 * np.diag(tau)
    gam = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a != b:
                gam[a, b] = (tau[b] * lam[a] ** 2 - tau[a] * lam[b] ** 2) / (lam[b] ** 2 - lam[a] ** 2)
    V = st.cache["nvec"]
    worst = 0.0
    for th, ph in ((0.3, 1.1), (1.2, 2.0), (2.5, 0.4)):
        alpha = np.array([np.sin(th) * np.sin(ph), np.cos(ph), np.cos(th) * np.sin(ph)])
        n = V @ alpha
        A1 = V.T @ acoustic_tensor(a4, n) @ V
        A2 = acoustic_principal(c_ab, gam, tau, alpha)
        worst = max(worst, np.abs(A1 - A2).max() / np.abs(A1).max())
    print(f"  (S.44) principal-basis acoustic tensor vs n.a.n contraction: max rel diff {worst:.2e}")
    # K2 nominal
    print("  K2 nominal: pi_i0 = -60.4, chi = -3.5, v_c0 = 1.81, WW, paper CSL, no cap, grid 91x91 + Nelder-Mead")
    results = {}
    for tag, (rho, rb) in (("0.7/0.8", (0.7, 0.8)), ("1.0/1.0", (1.0, 1.0))):
        Pk = k2_params(rho=rho, rho_bar=rb)
        s0 = initial_state(Pk, -100 * I3, 1.59, -60.4, finite=True)
        r = k2_path(Pk, s0, 60)
        results[tag] = r
        nd = r["min_det_normalised"]
        nv = r["n_vec"][-1]
        refused = any(s.flags["refused"] for s in r["states"])
        print(f"    rho/rho_bar {tag}: n_first = {r['n_first']}, n_interp = {r['n_interp']:.3f}" if r["n_first"] else
              f"    rho/rho_bar {tag}: no localization within n_max")
        print(f"      normalised min det from step 10: " + " ".join(f"{x:+.3f}" for x in nd[9:]))
        print(f"      n at localization (global = principal basis, diagonal path): {np.round(nv, 4)}, refused {refused}, "
              f"first plastic step {next((i + 1 for i, s in enumerate(r['states']) if s.flags['plastic']), None)}")
    if results["0.7/0.8"]["n_first"] and results["1.0/1.0"]["n_first"]:
        print(f"    ordering rho=0.7 before rho=1: {results['0.7/0.8']['n_first'] < results['1.0/1.0']['n_first']}, "
              f"gap = {results['1.0/1.0']['n_first'] - results['0.7/0.8']['n_first']} (paper: 22 / 26, gap 4)")
    if full_table:
        print("  K2 sensitivity table (n_first / n_interp) [pi_i0, chi, v_c0]:")
        for pi0 in (-60.0, -80.0, -100.0):
            for chi in (-3.0, -3.5, -4.0):
                for vc0 in (1.80, 1.81, 1.82):
                    row = []
                    for rho, rb in ((0.7, 0.8), (1.0, 1.0)):
                        Pk = k2_params(rho=rho, rho_bar=rb, chi=chi, v_c0=vc0)
                        r = k2_path(Pk, initial_state(Pk, -100 * I3, 1.59, pi0, finite=True), 60,
                                    acoustic_kwargs=dict(grid_n=46, n_refine=2))
                        row.append((r["n_first"], r["n_interp"]))
                    print(f"    pi_i0 {pi0:6.1f} chi {chi:4.1f} v_c0 {vc0:.2f}: 0.7 -> {row[0][0]} ({row[0][1]:.2f})  1.0 -> {row[1][0]} ({row[1][1]:.2f})  gap {row[1][0] - row[0][0]}")


# ------------------------------------------------------------------ 6. smooth cap: root selection + substepping
def check_cap():
    """Near-isotropic compression through the smooth-cap ramp (the G1.cap path: -0.01 on every axis
    plus amp*diag(1,0,-1), c1 = 0.05, c2 = 0.15, paper set, pi_i0 = -80, psi_i0 = -0.05). Reports
    completion, the substep census and the nested-solve work; prints the ramp-entry steps of the
    n = 40 / amp = 2e-3 case (the step the factor-2 bracket used to send to the far root)."""
    print("\n=== 6. smooth cap: nested root selection (pi_fold guard) and substepping ===")
    P = k2_params(cap="smooth", c1=0.05, c2=0.15)
    v0 = -0.05 + P.v_c0 - P.lam_tilde * math.log(80.0)
    e1 = P.c1 * P.M
    for amp in (1e-4, 5e-4, 2e-3):
        for n in (40, 80, 160, 320, 640):
            s0 = initial_state(P, -100 * I3, v0, -80.0)
            deps = (-0.01 * I3 + amp * np.diag([1.0, 0.0, -1.0])) / n
            sts = run_path(P, s0, np.array([deps] * n))
            ref = [s.flags["refused"] for s in sts]
            k = ref.index(True) if any(ref) else None
            subs = [s.flags.get("substeps", 1) for s in sts]
            nsub = sum(1 for s in subs if s > 1)
            Dmin = min(s.D for s in sts)
            npl = sum(1 for s in sts if s.flags["plastic"])
            pit = sum(s.flags["pi_iters"] for s in sts)
            lit = sum(s.flags["local_iters"] for s in sts)
            w_end = K.cap_weight(P, sts[-1].eta)[0] if sts[-1].flags["plastic"] else float("nan")
            print(f"  amp {amp:g} n {n:4d}: " + (f"REFUSED at step {k}: {sts[k].flags['reason']}" if k is not None else "completed")
                  + f"; plastic steps {npl}; substepped increments {nsub} (max substeps {max(subs)}); local iters {lit}, "
                  f"nested iters {pit}; min D {Dmin:+.3e}; end eta {sts[-1].eta:.4f} w {w_end:.3e} pi_i {sts[-1].pi_i:.4f}")
    print("  ramp entry, amp 2e-3 n 40 (eta_1 = c1 M = %.3f): step, dlam, pi_i, eta, w, substeps" % e1)
    s0 = initial_state(P, -100 * I3, v0, -80.0)
    deps = (-0.01 * I3 + 2e-3 * np.diag([1.0, 0.0, -1.0])) / 40
    sts = run_path(P, s0, np.array([deps] * 40))
    for i, s in enumerate(sts):
        if s.flags["plastic"] and (abs(s.eta - e1) < 0.03 or s.flags.get("substeps", 1) > 1):
            print(f"    {i + 1:3d}: dlam {s.dlam:.4e}  pi_i {s.pi_i:.4f}  eta {s.eta:.5f}  w {K.cap_weight(P, s.eta)[0]:.3e}  "
                  f"substeps {s.flags.get('substeps', 1)}  local iters {s.flags['local_iters']}  nested iters {s.flags['pi_iters']}")


# ------------------------------------------------------------------ 7. §9.6 chained tangent across substeps
def chain_vs_fd(P, st, deps, fractions, h):
    """Chained tangent (S.47) of the increment deps taken with the given fractions, against the central FD
    of the WHOLE increment (fractions held fixed at every FD point, same branch pattern required).
    Returns (max over the 6 kernel columns of ||C_J - FD_J||/||C_J||, the same for the last-sub-increment
    CTO, the state, its branch pattern)."""
    stn = step_fractions(P, st, deps, fractions)
    assert not stn.flags["refused"], stn.flags["reason"]
    C, Cl = tangent(P, stn), tangent_last_substep(P, stn)
    pat = stn.flags["pattern"]
    fpat = stn.flags["fpattern"]
    e, el = 0.0, 0.0
    for J in range(6):
        E = K.CHAIN_E[J]
        sp = step_fractions(P, st, deps + h * E, fractions)
        sm = step_fractions(P, st, deps - h * E, fractions)
        assert sp.flags["pattern"] == pat == sm.flags["pattern"], (pat, sp.flags["pattern"], sm.flags["pattern"])
        assert sp.flags["fpattern"] == fpat == sm.flags["fpattern"], (fpat, sp.flags["fpattern"], sm.flags["fpattern"])
        fd = (sp.sigma - sm.sigma) / (2.0 * h)
        col, coll = np.einsum("ijkl,kl->ij", C, E), np.einsum("ijkl,kl->ij", Cl, E)
        e = max(e, np.linalg.norm(col - fd) / np.linalg.norm(col))
        el = max(el, np.linalg.norm(coll - fd) / np.linalg.norm(coll))
    return e, el, stn, pat


def _s33(st):
    c = st.cache
    return K.tangent_small(c["atilde"], c["sig"], c["eps_tr"], c["nvec"])


def check_chain():
    print("\n=== 7. chained consistent tangent across substeps (sheet §9.6) vs central FD of the whole increment ===")
    H = (1e-6, 1e-7, 1e-8)
    # (A) AMP_STOP smooth cap, n = 40: every substepped increment of the ladder, FD with the ladder's 1/m
    P = k2_params(cap="smooth", c1=0.05, c2=0.15)
    v0 = -0.05 + P.v_c0 - P.lam_tilde * math.log(80.0)
    s0 = initial_state(P, -100 * I3, v0, -80.0)
    deps = (-0.01 * I3 + 2e-3 * np.diag([1.0, 0.0, -1.0])) / 40
    sts = run_path(P, s0, np.array([deps] * 40))
    print("  (A) AMP_STOP smooth cap n = 40 (ladder substeps; FD with the same fractions): step m pattern | "
          "chain err h=1e-6 1e-7 1e-8 | last-sub CTO err")
    amp_states, rows = {}, []
    for i, s in enumerate(sts):
        m = s.flags.get("substeps", 1)
        if m == 1 or s.flags["refused"]:
            continue
        prev = sts[i - 1] if i else s0
        amp_states[i + 1] = prev
        errs, last = [], None
        for h in H:
            e, el, stn, pat = chain_vs_fd(P, prev, deps, [1.0 / m] * m, h)
            errs.append(e)
            last = el
        assert np.array_equal(stn.sigma, s.sigma) and stn.pi_i == s.pi_i, "ladder state != fixed-fraction state"
        assert np.array_equal(tangent(P, stn), tangent(P, s)), "ladder chain != fixed-fraction chain"
        rows.append((i + 1, m, pat, errs, last))
        print(f"    {i + 1:3d} {m} {pat:5s} | {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} | {last:.2e}")
    if rows:
        for hi, h in enumerate(H):
            v = [r[3][hi] for r in rows]
            print(f"      h = {h:.0e}: chain err min {min(v):.2e} max {max(v):.2e} over {len(rows)} substepped increments")
        v = [r[4] for r in rows]
        print(f"      last-sub-increment CTO err: min {min(v):.2f} max {max(v):.2f}")
    # (B) generic plastic increment with all three shears (no cap), forced m = 8 and m = 2
    shear = np.array([[0, 3e-4, 1e-4], [3e-4, 0, 2e-4], [1e-4, 2e-4, 0]])
    Pb = fork_params()
    sb0 = initial_state(Pb, -100 * I3, 1.65, -60.4)
    stb = run_path(Pb, sb0, np.array([np.diag([4e-4, -1e-3, 0.0]) + shear] * 5))[-1]
    db = np.diag([1e-4, -6e-4, 2e-4]) + 0.5 * shear
    print("  (B) generic non-coaxial plastic increment (fork WW, no cap), forced uniform fractions:")
    for m in (8, 2):
        errs, last = [], None
        for h in H:
            e, el, stn, pat = chain_vs_fd(Pb, stb, db, [1.0 / m] * m, h)
            errs.append(e)
            last = el
        th = K.invariants(np.linalg.eigvalsh(stn.sigma)).theta
        print(f"    m = {m} pattern {pat}: chain err {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} | last-sub CTO err {last:.2e} "
              f"(theta {th:.3f})")
    # (C) m = 1: the chain must reduce to (S.33); the ladder's m = 1 path returns (S.33) itself
    s1 = step_fractions(Pb, stb, db, [1.0])
    dC = np.linalg.norm(tangent(Pb, s1) - _s33(s1)) / np.linalg.norm(_s33(s1))
    e1, _, _, pat1 = chain_vs_fd(Pb, stb, db, [1.0], 1e-7)
    s1e = step_fractions(Pb, sb0, np.diag([1e-4, -2e-4, 0.5e-4]) + shear, [1.0])
    dCe = np.linalg.norm(tangent(Pb, s1e) - _s33(s1e)) / np.linalg.norm(_s33(s1e))
    s1l = step(Pb, stb, db)
    ident = ("C_chain" not in s1l.cache) and np.array_equal(tangent(Pb, s1l), _s33(s1))
    print(f"  (C) m = 1: chain vs (S.33) {dC:.1e} (plastic, pattern {pat1}, FD err {e1:.2e}); elastic step {dCe:.1e} "
          f"(plastic {s1e.flags['plastic']}); ladder m = 1 returns (S.33) bit-identically: {ident}")
    # (E) non-uniform fractions (recursive-halving shapes, sum = 1)
    print("  (E) non-uniform fractions:")
    cases = [("(B) increment", Pb, stb, db, (0.5, 0.25, 0.125, 0.125))]
    if 20 in amp_states:
        cases.append(("AMP step 20", P, amp_states[20], deps, (0.25, 0.25, 0.25, 0.125, 0.125)))
    if 11 in amp_states:
        cases.append(("AMP step 11", P, amp_states[11], deps, (0.5, 0.25, 0.25)))
    for name, Pc, sc, dc, fr in cases:
        errs, last = [], None
        for h in H:
            e, el, stn, pat = chain_vs_fd(Pc, sc, dc, fr, h)
            errs.append(e)
            last = el
        print(f"    {name:14s} alpha = {fr}: pattern {pat}: chain err {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} "
              f"| last-sub {last:.2e}")
    # (F) vertex branch: hydrostatic plastic step from the apex (no cap): pi_i frozen, C:1 = 0, FD along 1
    print("  (F) vertex branch, hydrostatic plastic step from the apex (no cap):")
    Pv = k2_params()
    sv0 = initial_state(Pv, -100 * I3, 1.59, None)
    dv = -1e-3 * I3
    for fr in ((1.0,), (0.5, 0.25, 0.25)):
        stn = step_fractions(Pv, sv0, dv, fr)
        C = tangent(Pv, stn)
        ones = np.einsum("ijkl,kl->ij", C, I3)
        h = 1e-7
        fd = (step_fractions(Pv, sv0, dv + h * I3, fr).sigma - step_fractions(Pv, sv0, dv - h * I3, fr).sigma) / (2 * h)
        msg = (f"    alpha = {fr}: pattern {stn.flags['pattern']}, vertex {stn.flags['vertex']}, p = {stn.sigma.trace() / 3:.4f}, "
               f"pi_i = {stn.pi_i:.4f} (frozen: {stn.pi_i == sv0.pi_i}), max|C:1|/max|C| = {np.abs(ones).max() / np.abs(C).max():.1e}, "
               f"max|FD along 1|/max|C| = {np.abs(fd).max() / np.abs(C).max():.1e}")
        if len(fr) == 1:
            msg += f", chain vs (S.33) {np.linalg.norm(C - _s33(stn)) / np.linalg.norm(_s33(stn)):.1e}"
        print(msg)


# ------------------------------------------------------------------ 8. HAR energy option (sheet §2.3-§2.4)
def check_har():
    print("\n=== 8. HAR energy option (sheet §2.3, §2.4, K1.1h, K1.11): closed forms, Hessian, inverse, FD re-runs with D12 != 0 ===")
    allok = True
    # (a) K1.1h / K1.11 closed forms at the TIMs p_a = 101 (A4; the sheet's values) and, as a p_a-scaling check, at 101.325
    for pa in (101.0, 101.325):
        P = tims_params(p_a=pa)
        kn = P.k * (1 - P.n_e)
        el = K.elastic(P, np.array([-1e-3 / 3] * 3))
        el2 = K.elastic(P, np.array([5e-4 / 3] * 3))
        edge = 1.0 / kn
        es = 8.66546968e-4
        e_cv = K.elastic(P, math.sqrt(1.5) * es * np.array([1, 0, -1]) / math.sqrt(2))     # eps_v = 0, eps_s = es
        print(f"  p_a = {pa}: K1.1h p(eps_v=-1e-3) = {el.p:.6f} (sheet @101: -381.983585), K = D11 = {el.D11:.3f} (371129.585), "
              f"p(+5e-4) = {el2.p:.6f} (-28.117707), domain edge eps_v = {edge:.9e} (1.058491699e-3)")
        print(f"            K1.11 constant-volume shear eps_s = {es:.8e}: eta = {e_cv.q / abs(e_cv.p):.9f} (3 g eps_s = {3 * P.g * es:.9f}), "
              f"p = {e_cv.p:.6f} (-166.548671), q = {e_cv.q:.6f} (349.752210)")
        allok &= ok("K1.11 eta = 3 g eps_s exactly", abs(e_cv.q / abs(e_cv.p) - 3 * P.g * es) / (3 * P.g * es), 1e-13)
        # ring spot inverse (S.5h''): p = -3.5, eta = 2.1
        sig = -3.5 * K.ONES + K.SQ23 * 2.1 * 3.5 * np.array([1, 0, -1]) / math.sqrt(2)
        eps = K.invert_elastic(P, sig)
        ev, es2, _ = K._split(eps)
        print(f"            inverse map p = -3.5, eta = 2.1: eps_v = {ev:.7e} (9.0504735e-4), eps_s = {es2:.7e} (1.2561906e-4); "
              f"round trip |sigma - sigma(eps)| = {np.abs(K.elastic(P, eps).sig - sig).max():.1e}")
        allok &= ok("inverse map round trip", np.abs(K.elastic(P, eps).sig - sig).max() / 3.5, 1e-13)
    # (b) Hessian (S.5h') stress form vs FD of (p, q) in (eps_v, eps_s); (S.3) vs FD of sigma_a; det D > 0; Psi derivative
    P = tims_params()
    rng = np.random.default_rng(3)
    worst = dict(D=0.0, ae=0.0, psi=0.0)
    for _ in range(12):
        ev = rng.uniform(-4e-3, 9e-4)
        es = rng.uniform(0.0, 2e-3)
        nh = rng.normal(size=3)
        nh -= nh.mean()
        nh /= np.linalg.norm(nh)
        eps = ev / 3 * K.ONES + math.sqrt(1.5) * es * nh
        el = K.elastic(P, eps)
        h = 1e-7

        def pq(dv, ds):
            e_ = (ev + dv) / 3 * K.ONES + math.sqrt(1.5) * (es + ds) * nh
            x = K.elastic(P, e_)
            return x.p, x.q
        D11 = (pq(h, 0)[0] - pq(-h, 0)[0]) / (2 * h)
        D22 = (pq(0, h)[1] - pq(0, -h)[1]) / (2 * h)
        D12 = (pq(0, h)[0] - pq(0, -h)[0]) / (2 * h)
        D21 = (pq(h, 0)[1] - pq(-h, 0)[1]) / (2 * h)
        Dfd = np.array([[D11, D12], [D21, D22]])
        Dcf = np.array([[el.D11, el.D12], [el.D12, el.D22]])
        worst["D"] = max(worst["D"], relmax(Dcf, Dfd))
        assert el.D11 * el.D22 - el.D12 ** 2 > 0 and el.D12 < 0
        ae_fd = np.zeros((3, 3))
        for b in range(3):
            e_ = np.zeros(3)
            e_[b] = h
            ae_fd[:, b] = (K.elastic(P, eps + e_).sig - K.elastic(P, eps - e_).sig) / (2 * h)
        worst["ae"] = max(worst["ae"], relmax(el.ae, ae_fd))
        gpsi = np.array([(K.energy_psi(P, eps + h * np.eye(3)[b]) - K.energy_psi(P, eps - h * np.eye(3)[b])) / (2 * h) for b in range(3)])
        worst["psi"] = max(worst["psi"], relmax(gpsi, el.sig))
    allok &= ok("(S.5h') Hessian D vs FD of (p, q), 12 random states", worst["D"], 1e-6)
    allok &= ok("(S.3) a^e vs FD of sigma_a (t2/t4 terms live: D12 != 0, D22 != q/eps_s)", worst["ae"], 1e-6)
    allok &= ok("sigma = dPsi/deps (S.4h)", worst["psi"], 1e-6)
    el = K.elastic(P, np.array([-1e-3, -1.2e-3, -0.3e-3]))
    print(f"      sample: D22/(q/eps_s) = {el.D22 / (el.q / el.eps_s):.4f} (in [1, 1/(1-n)) = [1, 2)), |D12|/sqrt(D11 D22) = "
          f"{abs(el.D12) / math.sqrt(el.D11 * el.D22):.3f}")
    # (c) K1.2 closed non-coaxial elastic loop under HAR (the HAR law is stiffer: the K1.2 loop strains are scaled by
    #     0.3 to stay elastic inside pi_i = -2e4, and |p|^n is not polynomial, so 40 Gauss points per edge)
    s0 = initial_state(P, -100 * I3, 1.65, -2.0e4)
    verts = [0.3 * x for x in (np.zeros((3, 3)), np.diag([1e-3, -2e-3, 5e-4]) + np.array([[0, 5e-4, 0], [5e-4, 0, 2e-4], [0, 2e-4, 0]]),
             np.diag([-1e-3, -1e-3, -2e-3]) + np.array([[0, 0, 4e-4], [0, 0, 0], [4e-4, 0, 0]]),
             np.diag([5e-4, 0, -1e-3]) + np.array([[0, -3e-4, 0], [-3e-4, 0, 0], [0, 0, 0]]), np.zeros((3, 3)))]
    xg, wg = np.polynomial.legendre.leggauss(40)
    W, Wabs = 0.0, 0.0
    for a, b in zip(verts[:-1], verts[1:]):
        d = b - a
        for xi_, wi in zip(xg, wg):
            e = a + 0.5 * (xi_ + 1) * d
            w, V = np.linalg.eigh(s0.eps_e + e)
            sig = (V * K.elastic(P, w).sig) @ V.T
            W += 0.5 * wi * float(np.sum(sig * d))
            Wabs += 0.5 * wi * abs(float(np.sum(sig * d)))
    sts = run_path(P, s0, np.array([b - a for a, b in zip(verts[:-1], verts[1:])]))
    print(f"  K1.2 closed loop (HAR): W/sum|W| = {W / Wabs:.2e}; all elastic {all(not s.flags['plastic'] for s in sts)}; "
          f"state return |d sigma| {np.abs(sts[-1].sigma - s0.sigma).max():.2e}")
    allok &= ok("K1.2 (HAR) W/sum|W|", abs(W / Wabs), 1e-12)
    # (d) Jacobian (S.30) vs FD under HAR (D12 != 0 live)
    shear = np.array([[0, 2e-4, 1e-4], [2e-4, 0, 0], [1e-4, 0, 0]])
    s0 = initial_state(P, -100 * I3, 1.65, -60.4)
    path = [np.diag([4e-4, -1e-3, 0.0]) + shear] * 6 + [np.diag([1e-4, -4e-4, 1e-3]) + 2 * shear] * 6
    sts = run_path(P, s0, np.array(path))
    jmax = 0.0
    for k in (4, 8, 11):
        errs = [jacobian_fd(P, sts[k - 1], path[k], h)[0] for h in (1e-6, 1e-7, 1e-8)]
        res = jacobian_fd(P, sts[k - 1], path[k])[1]
        jmax = max(jmax, errs[2])
        print(f"  Jacobian (HAR) state {k:2d}: rel err h=1e-6/1e-7/1e-8: {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e}  "
              f"(dlam {res.dlam:.3e}, iters {res.local_iters}, theta {K.invariants(res.sig).theta:.4f}, D12 {res.ae[0, 1]:.1f})")
    allok &= ok("Jacobian (S.30) vs FD under HAR, max at h=1e-8", jmax, 1e-7)
    # (e) CTO (S.33) vs FD under HAR: plastic non-coaxial off corners, elastic with q > 0 (t4 live), finite mode
    shear3 = np.array([[0, 3e-4, 1e-4], [3e-4, 0, 2e-4], [1e-4, 2e-4, 0]])
    st = run_path(P, s0, np.array([np.diag([4e-4, -1e-3, 0.0]) + shear3] * 5))[-1]
    d = np.diag([1e-4, -6e-4, 2e-4]) + 0.5 * shear3
    for h in (1e-5, 1e-6, 1e-7):
        err, stn, C = cto_fd(P, st, d, h)
        th = K.invariants(np.linalg.eigvalsh(stn.sigma)).theta
        print(f"  CTO (HAR) plastic non-coaxial h={h:.0e}: rel err {err:.2e} (theta {th:.4f}, sym err "
              f"{np.linalg.norm(C.reshape(9, 9) - C.reshape(9, 9).T) / np.linalg.norm(C):.3f})  [O(h^2) expected]")
    allok &= ok("CTO (S.33) vs FD under HAR, plastic, h=1e-7", err, 1e-7)
    # elastic steps (the HAR law is ~15x stiffer in shear than BA06's mu0 = 5400: smaller increments, pi_i far inside)
    sfar = initial_state(P, -100 * I3, 1.65, -2.0e4)
    for h in (1e-6, 1e-7):
        err, stn, _ = cto_fd(P, sfar, 0.1 * (np.diag([1e-4, -2e-4, 0.5e-4]) + shear3), h)
        print(f"  CTO (HAR) elastic non-coaxial from the isotropic start h={h:.0e}: rel err {err:.2e}, plastic {stn.flags['plastic']} "
              "[O(h^2): the HAR moduli vary on the eps* ~ 1e-3 scale]")
    assert not stn.flags["plastic"]
    allok &= ok("CTO elastic (HAR), repeated-eigenvalue start, h=1e-7", err, 1e-7)
    se = run_path(P, sfar, np.array([0.1 * (np.diag([1e-4, -2e-4, 0.5e-4]) + shear3)] * 2))[-1]
    for h in (1e-6, 1e-7):
        err, stn, _ = cto_fd(P, se, 0.1 * (np.diag([0.5e-4, -1e-4, 0.3e-4]) + 0.3 * shear3), h)
        print(f"  CTO (HAR) elastic with q > 0 (t2 and t4 of (S.3) live) h={h:.0e}: rel err {err:.2e}, plastic {stn.flags['plastic']}")
    assert not stn.flags["plastic"]
    allok &= ok("CTO elastic (HAR) at q > 0, h=1e-7", err, 1e-7)
    # finite mode a~^ep vs FD of tau_a (diagonal protocol)
    s0f = initial_state(P, -100 * I3, 1.65, -60.4, finite=True)
    f1 = np.log(np.array([1 + 4e-4, 1 - 1e-3, 1.0]))
    stsf = run_path(P, s0f, np.array([np.diag(f1)] * 12))
    stf = stsf[-1]
    at = stf.cache["atilde"]
    Vc = stf.cache["nvec"]
    perm = [int(np.argmax(np.abs(Vc[:, a]))) for a in range(3)]
    for h in (1e-6, 1e-7):
        at_fd = np.zeros((3, 3))
        for b in range(3):
            e = np.zeros(3)
            e[perm[b]] = h
            sp = np.diag(step(P, stsf[-2], np.diag(f1 + e)).sigma)[perm]
            sm = np.diag(step(P, stsf[-2], np.diag(f1 - e)).sigma)[perm]
            at_fd[:, b] = (sp - sm) / (2 * h)
        print(f"  finite mode (HAR) a~^ep (S.32) vs FD h={h:.0e}: rel err {relmax(at, at_fd):.2e} (plastic {stf.flags['plastic']})")
    allok &= ok("finite mode (HAR): a~^ep (S.32) vs FD, h=1e-7", relmax(at, at_fd), 1e-7)
    # (f) §9.6 chain under HAR: (B) m = 8 / 2, (C) m = 1 reduction, (E) non-uniform, (F) vertex
    H = (1e-6, 1e-7, 1e-8)
    print("  chain (HAR), generic non-coaxial plastic increment:")
    for m in (8, 2):
        errs, last = [], None
        for h in H:
            e, el_, stn, pat = chain_vs_fd(P, st, d, [1.0 / m] * m, h)
            errs.append(e)
            last = el_
        print(f"    (B) m = {m} pattern {pat}: chain err {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} | last-sub CTO err {last:.2e}")
        allok &= ok(f"chain (B) m = {m} under HAR, h=1e-7", errs[1], 1e-7)
    s1 = step_fractions(P, st, d, [1.0])
    dC = np.linalg.norm(tangent(P, s1) - _s33(s1)) / np.linalg.norm(_s33(s1))
    print(f"    (C) m = 1: chain vs (S.33) {dC:.1e}")
    allok &= ok("chain (C) m = 1 reduces to (S.33) under HAR", dC, 1e-12)
    errs = [chain_vs_fd(P, st, d, (0.5, 0.25, 0.125, 0.125), h)[0] for h in H]
    print(f"    (E) alpha = (1/2, 1/4, 1/8, 1/8): chain err {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e}")
    allok &= ok("chain (E) non-uniform under HAR, h=1e-7", errs[1], 1e-7)
    sv0 = initial_state(P, -100 * I3, 1.65, None)
    dv = -1e-3 * I3
    for fr in ((1.0,), (0.5, 0.25, 0.25)):
        stn = step_fractions(P, sv0, dv, fr)
        C = tangent(P, stn)
        ones = np.einsum("ijkl,kl->ij", C, I3)
        h = 1e-7
        fd = (step_fractions(P, sv0, dv + h * I3, fr).sigma - step_fractions(P, sv0, dv - h * I3, fr).sigma) / (2 * h)
        print(f"    (F) vertex alpha = {fr}: pattern {stn.flags['pattern']}, vertex {stn.flags['vertex']}, pi_i frozen "
              f"{stn.pi_i == sv0.pi_i}, max|C:1|/max|C| = {np.abs(ones).max() / np.abs(C).max():.1e}, "
              f"max|FD along 1|/max|C| = {np.abs(fd).max() / np.abs(C).max():.1e}")
    print("  HAR group:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ------------------------------------------------------------------ 9. the p' floor (sheet §9.7, K1.12-K1.14b)
def floor_fd(P, st, deps, fractions, h):
    """(S.32f)/(S.54) vs central FD of the whole increment over the six kernel columns, branch pattern (incl. the floor
    pattern) held fixed. Returns (chain err, direct (S.32f) CTO err of the last sub-increment, state)."""
    e, el, stn, _ = chain_vs_fd(P, st, deps, fractions, h)
    return e, el, stn


def check_floor():
    print("\n=== 9. the p' floor Pi_f (sheet §9.7, (S.48)-(S.54), K1.12-K1.14b, round-3b A1/A5) ===")
    allok = True
    # (a) K1.12 BA06, K2 set, default p_min = 5e-3 |p0| = 0.5 kPa
    P = k2_params()
    print(f"  (a) K1.12 BA06 K2 set: default p_min = {P.p_min} (5e-3 |p0|); eps_v,f(0) = {K.floor_ev(P, 0.0)[0]:.10f} (0.0529831737)")
    s0 = initial_state(P, -100 * I3, 1.59, -2.0e4)                        # far inside the surface (elastic at p = -0.5)
    ev_tr = K.floor_ev(P, 0.0)[0] + P.kappa_hat * math.log(2.0)            # p^tr = -p_min/2 (= 0.0599146455)
    shear = 1e-6 * np.array([[0, 1, 0.5], [1, 0, 0], [0.5, 0, 0]])         # a little non-coaxial shear: q, n^ to check
    d = ev_tr / 3 * I3 + shear
    st = step(P, s0, d)
    w_tr, V = np.linalg.eigh(s0.eps_e + d)
    el_tr = K.elastic(P, w_tr)
    el_c = K.elastic(P, np.linalg.eigvalsh(st.eps_e))
    print(f"      trial p = {el_tr.p:.6f} (-0.25) -> committed p = {np.trace(st.sigma) / 3:.10f} (-0.5); pattern {st.flags['fpattern']}; "
          f"d eps^f_v = {st.eps_f_v:.10e} (6.931471806e-3), W_f = {st.W_f:.10e} (3.465735903e-3), E_f = {floor_energy(P, st):.6e} (2.5e-3); "
          f"q {el_tr.q:.6f} -> {el_c.q:.6f}, pi_i {st.pi_i} (unchanged {st.pi_i == s0.pi_i}), v = v0 exp(tr): {abs(st.v - 1.59 * math.exp(ev_tr)):.1e}")
    assert st.flags["fpattern"] == "FE-", st.flags["fpattern"]
    allok &= ok("K1.12 committed p = -p_min", abs(np.trace(st.sigma) / 3 + 0.5) / 0.5, 1e-13)
    allok &= ok("K1.12 d eps^f_v = kappa ln 2", abs(st.eps_f_v - 0.01 * math.log(2)), 1e-15)
    allok &= ok("K1.12 E_f = kappa (p_min - |p_tr|)", abs(floor_energy(P, st) - 0.01 * 0.25) / 2.5e-3, 1e-9)
    allok &= ok("K1.12 E_f <= W_f", max(0.0, floor_energy(P, st) - st.W_f), 1e-18)
    C = tangent(P, st)
    at = st.cache["atilde"]
    allok &= ok("K1.12 a~_f = 2 mu0 (delta - 1/3)", relmax(at, 2 * P.mu0 * (I3 - 1.0 / 3.0)), 1e-13)
    allok &= ok("K1.12 delta : C_f = 0", np.abs(np.einsum("iikl->kl", C)).max() / np.abs(C).max(), 1e-14)
    allok &= ok("K1.12 q, n^ unchanged by Pi_f (BA06 alpha0 = 0)", abs(el_c.q - el_tr.q) / el_tr.q, 1e-13)
    # idempotence
    fl1 = K.floor_project(P, w_tr)
    fl2 = K.floor_project(P, fl1.eps_f)
    print(f"      idempotent: second projection active = {fl2.active}, p = {K.elastic(P, fl1.eps_f).p:.15f}")
    assert fl1.active and not fl2.active
    # (b) K1.13 HAR floor and the out-of-domain trial (sheet values at p_a = 101 (A4); 101.325 only shows they move with p_a)
    for pa in (101.0, 101.325):
        P = tims_params(p_a=pa)
        kn = P.k * (1 - P.n_e)
        evf0 = K.floor_ev(P, 0.0)[0]
        G = P.g * pa * (P.p_min / pa) ** 0.5
        print(f"  (b) K1.13 HAR p_a = {pa}: p_min = {P.p_min:.6f} (5e-3 p_a), eps_v,f(0) = {evf0:.8e} (sheet @101: 9.83645033e-4), "
              f"edge 1/(k(1-n)) = {1 / kn:.8e}, G(p_min) = {G:.3f} (5769.156), 2G = {2 * G:.2f}, eps*_f(0) = {1 / kn - evf0:.4e}")
        s1 = initial_state(P, -1.0 * I3, 1.70, -50.0)
        ev1 = float(np.trace(s1.eps_e))
        rows = []
        for dv in (1e-4, 1.1e-4):
            stn = step(P, s1, dv / 3 * I3)
            try:
                ptr = K.elastic(P, np.linalg.eigvalsh(s1.eps_e) + dv / 3).p
            except K.EvalError:
                ptr = float("nan")
            rows.append((dv, ptr, stn))
            print(f"      isotropic p = -1 kPa (eps_v = {ev1:.8e}; 9.53167838e-4 @101) + d eps_v = {dv:+.1e}: p^tr = {ptr:.3e}, "
                  f"in domain {ev1 + dv < 1 / kn}, refused {stn.flags['refused']}, pattern {stn.flags['fpattern']}, "
                  f"d eps^f_v = {stn.eps_f_v:.7e} ({'6.9522805e-5' if dv == 1e-4 else '7.9522805e-5'} @101), W_f = {stn.W_f:.7e}, "
                  f"E_f = {floor_energy(P, stn):.4e} (bound {P.p_min * stn.eps_f_v:.4e}), committed p = {np.trace(stn.sigma) / 3:.6f}")
            allok &= ok(f"K1.13 d eps_v = {dv:.1e}: not refused, committed p = -p_min", abs(np.trace(stn.sigma) / 3 + P.p_min) / P.p_min, 1e-12)
        allok &= ok("K1.13 the out-of-domain trial floors to the SAME eps_v,f (d eps^f_v differs by 1e-5)",
                    abs((rows[1][2].eps_f_v - rows[0][2].eps_f_v) - 1e-5), 1e-15)
        P0 = tims_params(p_a=pa, p_min=0.0)
        s1 = initial_state(P0, -1.0 * I3, 1.70, -50.0)
        stn = step(P0, s1, 1.1e-4 / 3 * I3)
        print(f"      p_min = 0: refused {stn.flags['refused']} reason '{stn.flags['reason']}' (the pre-round-3 refusal, M-F5)")
        assert stn.flags["refused"] and "elastic_domain" in stn.flags["reason"]
    Pb = k2_params()
    s1 = initial_state(Pb, -1.0 * I3, 1.59, -50.0)
    stn = step(Pb, s1, 1.1e-4 / 3 * I3)
    print(f"      the same +1.1e-4 under BA06 (K2 set): pattern {stn.flags['fpattern']}, plastic {stn.flags['plastic']}, "
          f"p = {np.trace(stn.sigma) / 3:.6f} (ordinary elastic step, M-F4)")
    # (c) K1.14 HAR floored state under shear: eps_s = 2e-4 at the floor
    for pa in (101.0, 101.325):
        P = tims_params(p_a=pa)
        es = 2e-4
        evf, epsp = K.floor_ev(P, es)
        kn = P.k * (1 - P.n_e)
        a_, b_ = 3 * kn * P.g * es * es, (P.p_min / pa) ** 2
        x = 0.5 * (a_ + math.sqrt(a_ * a_ + 4 * b_))
        nh = np.array([1, 0, -1]) / math.sqrt(2)
        eps_f = evf / 3 * K.ONES + math.sqrt(1.5) * es * nh
        el = K.elastic(P, eps_f)
        hh = 1e-7
        epsp_fd = (K.floor_ev(P, es + hh)[0] - K.floor_ev(P, es - hh)[0]) / (2 * hh)
        print(f"  (c) K1.14 HAR p_a = {pa}, eps_s = 2e-4 at the floor: x = {x:.11f} (0.09185198375 @101), eps_v,f = {evf:.8e} (1.04102893e-3), "
              f"q_f = {el.q:.7f} (14.8362051), eta_f = {el.q / abs(el.p):.2f}, eps'_f = {epsp:.8f} (0.08679793); -D12/D11 = {-el.D12 / el.D11:.8f}, "
              f"FD {epsp_fd:.8f}; p = {el.p:.9f}")
        allok &= ok("K1.14 eps'_f = -D12/D11 at p = -p_min", abs(epsp + el.D12 / el.D11) / epsp, 1e-10)
        allok &= ok("K1.14 p(eps_v,f(eps_s), eps_s) = -p_min", abs(el.p + P.p_min) / P.p_min, 1e-13)
    # the elastic floored tangent a^e Phi (K1.14): delta:C_f = 0 with the eps'_f term, 1.45e-2 without (M-F3b)
    #     The K1.14 state has eta_f = 29 >> M/N (past the tension apex of every surface), so it is never an elastic STEP of
    #     the model: the sheet's K1.14 is a check of the elastic floor OPERATOR a^e(eps_f) Phi (S.51a) itself (floor_sympy (d)).
    P = tims_params(p_a=101.0)
    es = 2e-4
    nh = np.array([1, 0, -1]) / math.sqrt(2)
    eps_tr = (K.floor_ev(P, es)[0] + 2e-5) / 3 * K.ONES + math.sqrt(1.5) * es * nh        # above the floor (in domain)
    fl = K.floor_project(P, eps_tr)
    el = K.elastic(P, fl.eps_f)
    Cop = el.ae @ fl.Phi
    Cmut = el.ae @ (I3 - 1.0 / 3.0)
    hh = 1e-7
    Cfd = np.zeros((3, 3))
    for b in range(3):
        e_ = np.zeros(3)
        e_[b] = hh
        Cfd[:, b] = (K.elastic(P, K.floor_project(P, eps_tr + e_).eps_f).sig - K.elastic(P, K.floor_project(P, eps_tr - e_).eps_f).sig) / (2 * hh)
    print(f"      K1.14 operator a^e Phi at eps_s = 2e-4 (active {fl.active}, eps_s after Pi_f {el.eps_s:.6e}, q_f = {el.q:.7f}): "
          f"vs FD of sigma(Pi_f(eps)) {relmax(Cop, Cfd):.2e}; delta:C/max|C| = {np.abs(Cop.sum(axis=0)).max() / np.abs(Cop).max():.1e}; "
          f"mutant eps'_f dropped (M-F3b): {np.abs(Cmut.sum(axis=0)).max() / np.abs(Cmut).max():.3e} (sheet 1.45e-2)")
    allok &= ok("K1.14 elastic floor operator a^e Phi vs FD", relmax(Cop, Cfd), 1e-7)
    allok &= ok("K1.14 delta : C_f = 0 (operator, HAR)", np.abs(Cop.sum(axis=0)).max() / np.abs(Cop).max(), 1e-13)
    allok &= ok("K1.14 eps_s unchanged by Pi_f (M-F6)", abs(el.eps_s - es) / es, 1e-12)
    # an actual elastic trial-floored STEP under HAR (low eta, pi_i far inside): FE- with the (S.32f) elastic tangent
    es2 = 2e-6
    eps0 = (K.floor_ev(P, es2)[0] - 2e-5) / 3 * K.ONES + math.sqrt(1.5) * es2 * nh
    s0 = initial_state(P, np.diag(K.elastic(P, eps0).sig), 1.70, -2.0e4)
    stn = step(P, s0, 3e-5 / 3 * I3 + 1e-7 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]))     # p ~ -0.81 -> trial -0.38 kPa
    C = tangent(P, stn)
    el2 = K.elastic(P, np.linalg.eigvalsh(stn.eps_e))
    print(f"      elastic trial-floored step (HAR, eta_f = {el2.q / abs(el2.p):.3f}): pattern {stn.flags['fpattern']}, "
          f"delta:C/max|C| = {np.abs(np.einsum('iikl->kl', C)).max() / np.abs(C).max():.1e}, p = {np.trace(stn.sigma) / 3:.6f}")
    assert stn.flags["fpattern"] == "FE-", stn.flags["fpattern"]
    allok &= ok("HAR elastic FE- step: delta : C_f = 0", np.abs(np.einsum("iikl->kl", C)).max() / np.abs(C).max(), 1e-13)
    # (d) FD of (S.32f) under BA06, K2 set, floor scaled to 50 kPa (the §9.7 record): (A) elastic + trial floor,
    #     (B) plastic + trial floor (FP-, dry side), (C) plastic + post floor (-Pf, wet side); alpha0 = 0 and 5
    print("  (d) (S.32f) vs FD, BA06 K2 set, p_min = 50 kPa (chain m = 1 | direct CTO), six kernel columns:")
    for alpha0 in (0.0, 5.0):
        P = k2_params(alpha0=alpha0, p_min=50.0)
        sA = initial_state(P, np.diag([-57.0, -61.0, -64.0]), 1.70, -300.0)
        dA = np.diag([2.0e-3, 1.5e-3, 2.5e-3]) + 2e-4 * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        for h in (1e-6, 1e-7):
            e, el_, stn = floor_fd(P, sA, dA, [1.0], h)
            print(f"      alpha0 {alpha0}: (A) {stn.flags['fpattern']} h={h:.0e}: chain {e:.2e} | direct {el_:.2e}; p = {np.trace(stn.sigma) / 3:.4f}")
        assert stn.flags["fpattern"] == "FE-", stn.flags["fpattern"]
        allok &= ok(f"(A) elastic + trial floor, alpha0 {alpha0}", max(e, el_), 1e-8)
        sigB, thB, nhB = off_corner_sig(P, -55.0, 1.4)
        piB = K.pi_of_eta(P, -55.0, 1.4)
        sB = initial_state(P, np.diag(sigB), 1.70, piB)
        dB = None
        for sdev in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3):
            cand = np.diag(1.0e-3 * K.ONES + sdev * nhB) + 1e-4 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
            stn = step(P, sB, cand)
            if stn.flags["fpattern"] == "FP-":
                dB = cand
                break
        assert dB is not None, "case (B) not found"
        for h in (1e-6, 1e-7):
            e, el_, stn = floor_fd(P, sB, dB, [1.0], h)
            print(f"      alpha0 {alpha0}: (B) {stn.flags['fpattern']} h={h:.0e}: chain {e:.2e} | direct {el_:.2e}; p_c = {np.trace(stn.sigma) / 3:.4f}, "
                  f"eta {stn.eta:.3f}, dlam {stn.dlam:.3e}")
        allok &= ok(f"(B) plastic + trial floor, alpha0 {alpha0}", max(e, el_), 1e-6)
        sigC, _, _ = off_corner_sig(P, -51.0, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
        piC = K.pi_of_eta(P, -51.0, 0.5 * P.M)
        sC = initial_state(P, np.diag(sigC), 1.70, piC)
        dC = None
        for amp in (2e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1.2e-2):
            cand = np.diag(amp * np.array([0.55, 0.40, -0.95]) - 1e-5) + 1e-4 * np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
            stn = step(P, sC, cand)
            if stn.flags["fpattern"] == "-Pf":
                dC = cand
                break
        assert dC is not None, "case (C) not found"
        for h in (1e-6, 1e-7):
            e, el_, stn = floor_fd(P, sC, dC, [1.0], h)
            print(f"      alpha0 {alpha0}: (C) {stn.flags['fpattern']} h={h:.0e}: chain {e:.2e} | direct {el_:.2e}; p = {np.trace(stn.sigma) / 3:.4f}, "
                  f"eta {stn.eta:.3f}; delta:C/max|C| = {np.abs(np.einsum('iikl->kl', tangent(P, stn))).max() / np.abs(tangent(P, stn)).max():.1e}")
        allok &= ok(f"(C) plastic + post floor, alpha0 {alpha0}", max(e, el_), 1e-6)
        # (D) m = 2 chain FP-,FP- on the (B) increment
        for h in (1e-6, 1e-7):
            e, el_, stn = floor_fd(P, sB, dB, [0.5, 0.5], h)
            print(f"      alpha0 {alpha0}: (D) m = 2 '{stn.flags['fpattern']}' h={h:.0e}: chain {e:.2e} | last-sub {el_:.2e}")
        allok &= ok(f"(D) chain (S.54) m = 2, alpha0 {alpha0}", e, 1e-6)
    # (e) HAR dry-side FPf (K1.14b, A1): TIMs set, p_min 0.505, surface state p = -0.6 kPa, eta = 1.2 M, off-corner, expansion + shear
    print("  (e) K1.14b HAR FPf (A1): TIMs set (p_a 101, A4), p_min = 0.505 (= the default 5e-3 p_a):")
    P = tims_params(p_min=0.505)
    sigF, thF, nhF = off_corner_sig(P, -0.6, 1.2 * P.M)
    piF = K.pi_of_eta(P, -0.6, 1.2 * P.M)
    vF = v_for_psi(P, piF, -0.10)
    sF = initial_state(P, np.diag(sigF), vF, piF)
    elF = K.elastic(P, np.linalg.eigvalsh(sF.eps_e))
    print(f"      state: p = {np.trace(sF.sigma) / 3:.4f}, eta = {1.2 * P.M:.4f}, theta = {thF:.3f}, pi_i = {piF:.5f}, v = {vF:.5f}, "
          f"eps_v = {elF.eps_v:.6e}, eps_s = {elF.eps_s:.6e}, D11 = {elF.D11:.1f}, D12 = {elF.D12:.1f}, D22 = {elF.D22:.1f}; "
          f"eps*_f(eps_s) = {1 / (P.k * (1 - P.n_e)) - K.floor_ev(P, elF.eps_s)[0]:.3e}")
    dF = np.diag((2e-5 / 3) * K.ONES + 2e-5 * nhF)
    stn = step(P, sF, dF)
    res = stn.cache["res"]
    assert stn.flags["fpattern"] == "FPf", stn.flags["fpattern"]
    try:
        ptr = K.elastic(P, res.eps_tr_raw).p
    except K.EvalError:
        ptr = float("nan")
    el_c = K.elastic(P, res.eps_c)
    inv_f = K.invariants(stn.cache["sig"])
    Ff = K.flow(P, inv_f, stn.pi_i).F
    dpv = stn.eps_p_v - sF.eps_p_v
    print(f"      pattern {stn.flags['fpattern']}: p^tr = {ptr:.4f} -> floored {K.elastic(P, res.eps_tr).p:.4f} (d eps^f_v,tr {res.dfv_tr:.3e}) -> "
          f"return p_c = {el_c.p:.4f} (ABOVE the floor) -> post floor p = {np.trace(stn.sigma) / 3:.4f} (d eps^f_v,post {res.dfv_post:.3e}); "
          f"dlam = {stn.dlam:.3e}, eta_c = {stn.eta:.4f}, pi_i = {stn.pi_i:.5f}, d eps^p_v = {dpv:+.3e} (dilative; |.| < eps*_f), "
          f"q_c {el_c.q:.4f} -> q_f {inv_f.q:.4f}, F(sigma_f)/p_min = {Ff / P.p_min:+.4f} (A5: no F <= F_tol claimed), "
          f"E_f = {floor_energy(P, stn):.3e} <= W_f share {P.p_min * (res.dfv_tr + res.dfv_post):.3e}")
    allok &= ok("FPf committed p = -p_min", abs(np.trace(stn.sigma) / 3 + P.p_min) / P.p_min, 1e-12)
    errs = {}
    for h in (1e-7, 1e-8):
        e, el_, stn = floor_fd(P, sF, dF, [1.0], h)
        errs[h] = max(e, el_)
        print(f"      (S.32f) FPf vs FD h={h:.0e}: chain {e:.2e} | direct {el_:.2e}")
    allok &= ok("(S.32f) at FPf under HAR, h=1e-8", errs[1e-8], 1e-7)
    C = tangent(P, stn)
    allok &= ok("FPf delta : C_f = 0", np.abs(np.einsum("iikl->kl", C)).max() / np.abs(C).max(), 1e-13)
    # mutants at FPf (the direct CTO rebuilt by hand from the chain data)
    ch = res.chain
    v_ = stn.v
    ae_f = res.ae
    dx_full = ch.b[:3, :3] @ res.Phi_tr - np.outer(ch.u[:3] * ch.Pi_v * v_, K.ONES)
    dx_notr = ch.b[:3, :3] - np.outer(ch.u[:3] * ch.Pi_v * v_, K.ONES)
    E6 = K.CHAIN_E
    Cfd_cols = []
    for J in range(6):
        sp = step_fractions(P, sF, dF + 1e-8 * E6[J], [1.0]).sigma
        sm = step_fractions(P, sF, dF - 1e-8 * E6[J], [1.0]).sigma
        Cfd_cols.append((sp - sm) / 2e-8)
    def col_err(at_):
        C_ = K.tangent_small(at_, stn.cache["sig"], stn.cache["eps_tr"], stn.cache["nvec"])
        return max(np.linalg.norm(np.einsum("ijkl,kl->ij", C_, E6[J]) - Cfd_cols[J]) / np.linalg.norm(Cfd_cols[J]) for J in range(6))
    print(f"      mutants vs FD: Phi_post dropped {col_err(ae_f @ dx_full):.2e}; Phi_tr dropped {col_err(ae_f @ res.Phi_post @ dx_notr):.2e}; "
          f"plain a~^ep {col_err(K.elastic(P, res.eps_c).ae @ dx_notr):.2e}; exact {col_err(res.atilde):.2e}")
    # Jacobian vs FD at the FPf converged iterate (the floored trial is what the residual sees)
    x = np.append(res.eps_c, res.dlam)
    pe = K.evaluate(P, x[:3], x[3], res.eps_tr, v_, sF.pi_i)
    J = K.jacobian(P, pe)
    Jfd = np.empty((4, 4))
    hj = 1e-8
    for j in range(4):
        e_ = np.zeros(4)
        e_[j] = hj
        Jfd[:, j] = (K.evaluate(P, (x + e_)[:3], (x + e_)[3], res.eps_tr, v_, sF.pi_i).r
                     - K.evaluate(P, (x - e_)[:3], (x - e_)[3], res.eps_tr, v_, sF.pi_i).r) / (2 * hj)
    S = np.diag([1, 1, 1, 1 / P.p_ref])
    allok &= ok("Jacobian (S.30) vs FD at the FPf iterate (HAR)", np.linalg.norm(S @ (J - Jfd)) / np.linalg.norm(S @ J), 1e-6)
    # (f) chains under HAR: 'FPf,FPf' (2x the increment, halved) and '-P-,-Pf'
    print("  (f) chains (S.54) under HAR, m = 2, fractions 1/2, 1/2:")
    for h in (1e-7, 1e-8):
        e, el_, stn = floor_fd(P, sF, 2.0 * dF, [0.5, 0.5], h)
        print(f"      2x the FPf increment: pattern '{stn.flags['fpattern']}' h={h:.0e}: chain {e:.2e} | last-sub {el_:.2e}; "
              f"floor_tr {stn.flags['floor_tr']}, floor_post {stn.flags['floor_post']}")
    assert stn.flags["fpattern"] == "FPf,FPf", stn.flags["fpattern"]
    allok &= ok("chain 'FPf,FPf' under HAR, h=1e-8", e, 1e-7)
    d2 = None
    for amp_v in np.linspace(1.0e-5, 2.2e-5, 13):
        for amp_s in np.linspace(1.5e-5, 9.0e-5, 16):
            cand = np.diag((amp_v / 3) * K.ONES + amp_s * nhF)
            s2 = step_fractions(P, sF, cand, [0.5, 0.5])
            if not s2.flags["refused"] and s2.flags["fpattern"] == "-P-,-Pf":
                d2 = cand
                break
        if d2 is not None:
            break
    assert d2 is not None, "'-P-,-Pf' not found"
    for h in (1e-7, 1e-8):
        e, el_, stn = floor_fd(P, sF, d2, [0.5, 0.5], h)
        print(f"      tr deps = {np.trace(d2):.3e}: pattern '{stn.flags['fpattern']}' h={h:.0e}: chain {e:.2e} | last-sub {el_:.2e}; "
              f"unsplit pattern {step(P, sF, d2).flags['fpattern']}")
    allok &= ok("chain '-P-,-Pf' under HAR, h=1e-8 (O(h^2) from h=1e-7)", e, 3e-7)
    # (g) counting along a path through the floor (W_f always counted, E_f on demand; substepped sums)
    sts = run_path(P, sF, np.array([dF] * 4))
    print("  (g) counters along 4 FPf increments: " + "; ".join(
        f"{i + 1}: {s.flags['fpattern']} n_f_tr {s.n_f_tr} n_f_post {s.n_f_post} eps_f_v {s.eps_f_v:.3e} W_f {s.W_f:.3e} "
        f"E_f {floor_energy(P, s):.2e} at_floor {s.flags['at_floor']}" for i, s in enumerate(sts)))
    assert [s.n_f_tr for s in sts] == [1, 2, 3, 4] and [s.n_f_post for s in sts] == [1, 2, 3, 4]
    assert all(s.flags["at_floor"] for s in sts) and all(floor_energy(P, s) <= P.p_min * (s.eps_f_v - (sts[i - 1].eps_f_v if i else sF.eps_f_v)) + 1e-18 for i, s in enumerate(sts))
    s2 = step_fractions(P, sF, 2.0 * dF, [0.5, 0.5])
    print(f"      substepped (m = 2) increment: floor_tr {s2.flags['floor_tr']}, floor_post {s2.flags['floor_post']} (sums), "
          f"fpattern {s2.flags['fpattern']}, n_f_tr {s2.n_f_tr}, state.floor = {s2.floor}")
    assert s2.flags["floor_tr"] == 2 and s2.flags["floor_post"] == 2 and s2.n_f_tr == 2
    # init floor: a deck sigma0 above the floor is projected and counted
    si = initial_state(P, -0.2 * I3, 1.70, None)
    print(f"      initial_state at p = -0.2 kPa: n_f_init {si.n_f_init}, p = {np.trace(si.sigma) / 3:.4f}, eps_f_v {si.eps_f_v:.3e}, "
          f"at_floor {si.flags['at_floor']}, pi_i0 = {si.pi_i:.5f} (the (S.53) apex through the FLOORED p)")
    assert si.n_f_init == 1 and abs(np.trace(si.sigma) / 3 + P.p_min) < 1e-12
    # (h) general-n HAR floor (bracketed scalar solve) and n = 1/2: p(eps_v,f(eps_s), eps_s) = -p_min, eps'_f vs FD
    worst = 0.0
    for n in (0.3, 0.5, 0.7, 0.0):
        Pn = tims_params(n_e=n)
        for es in (0.0, 2e-4, 3e-3):
            evf, epsp = K.floor_ev(Pn, es)
            eps = evf / 3 * K.ONES + math.sqrt(1.5) * es * np.array([1, 0, -1]) / math.sqrt(2)
            el = K.elastic(Pn, eps)
            hh = 1e-7
            fd = (K.floor_ev(Pn, es + hh)[0] - K.floor_ev(Pn, es - hh)[0]) / (2 * hh)   # eps_v,f is even in eps_s
            worst = max(worst, abs(el.p + Pn.p_min) / Pn.p_min, abs(epsp - fd) / max(1e-3, abs(epsp)))
    allok &= ok("(S.50) general n in {0, 0.3, 0.5, 0.7}: p(eps_v,f) = -p_min and eps'_f vs FD", worst, 1e-6)
    # (i) defaults leave the default path unchanged: p_min = 0 vs the default floor on the K2 / fork drained paths (bit-identical)
    same = True
    for mk, v0 in ((k2_params, 1.72), (fork_params, 1.70)):
        Pa_, Pb_ = mk(p_min=0.0), mk()
        a_ = triaxial(Pa_, initial_state(Pa_, -100 * I3, v0, None), "drained", -0.05, 100)
        b_ = triaxial(Pb_, initial_state(Pb_, -100 * I3, v0, None), "drained", -0.05, 100)
        same &= all(np.array_equal(x.sigma, y.sigma) and x.pi_i == y.pi_i and np.array_equal(tangent(Pa_, x), tangent(Pb_, y))
                    for x, y in zip(a_, b_)) and len(a_) == len(b_) == 100
    print(f"  (i) default p_min = 5e-3 p_ref vs p_min = 0 on the K2 and fork drained paths (100 steps): bit-identical {same}")
    allok &= same
    # (j) NEAR-floor but not floored (HAR, the K1.14b state p = -0.6 kPa, p_min 0.505): a compressive + shear increment that
    #     stays above the floor, pattern -P-: the plain (S.33) CTO, bit-identical to the no-floor path, vs FD; and the
    #     sibling that only just floors at the trial (FP-/FPf): the one-sided switch is a branch decision (sheet §9.7)
    P = tims_params(p_min=0.505)
    sF = initial_state(P, np.diag(sigF), vF, piF)
    dN = np.diag((-1.5e-5 / 3) * K.ONES + 1.5e-5 * nhF) + 2e-6 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    stn = step(P, sF, dN)
    print(f"  (j) near-floor, not floored (HAR): pattern {stn.flags['fpattern']}, p_c = {np.trace(stn.sigma) / 3:.4f} (floor -0.505), "
          f"eta {stn.eta:.4f}, at_floor {stn.flags['at_floor']}")
    assert stn.flags["fpattern"] == "-P-", stn.flags["fpattern"]
    P0 = tims_params(p_min=0.0)
    st0 = step(P0, initial_state(P0, np.diag(sigF), vF, piF), dN)
    allok &= ok("(j) -P- near the floor: bit-identical to p_min = 0", float(np.abs(stn.sigma - st0.sigma).max()
                + np.abs(tangent(P, stn) - tangent(P0, st0)).max() + abs(stn.pi_i - st0.pi_i)), 0.0)
    for h in (1e-7, 1e-8):
        e, el_, stn = floor_fd(P, sF, dN, [1.0], h)
        print(f"      (S.33) -P- near the floor vs FD h={h:.0e}: chain {e:.2e} | direct {el_:.2e}")
    allok &= ok("(j) CTO (S.33) at the near-floor -P- state (HAR), h=1e-8", max(e, el_), 1e-7)
    # (k) finite mode (eps~ = log stretches, diagonal protocol): a~^ep_f of (S.32f) vs FD of tau_a over the principal log
    #     stretches, BA06 K2 set with p_min = 50 kPa, patterns FE- (trial floor, elastic) and FP- (trial floor, plastic);
    #     tangent_finite (S.34) takes a~^ep_f with the RAW trial stretches (sheet §9.7 "LogStrain provider")
    P = k2_params(p_min=50.0)
    print("  (k) finite mode, BA06 K2 set, p_min = 50 kPa: a~^ep_f vs FD of tau_a (diagonal log-stretch protocol):")
    sigB, _, nhB = off_corner_sig(P, -55.0, 1.4)
    piB = K.pi_of_eta(P, -55.0, 1.4)
    sfB = initial_state(P, np.diag(sigB), 1.70, piB, finite=True)
    ddB = None
    for sdev in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3):
        cand = 1.0e-3 * K.ONES + sdev * nhB
        if step(P, sfB, np.diag(cand)).flags["fpattern"] == "FP-":
            ddB = cand
            break
    assert ddB is not None, "finite-mode case FP- not found"
    for label, sig0, pi0, dd in (("FE-", np.diag([-57.0, -61.0, -64.0]), -300.0, np.array([2.0e-3, 1.5e-3, 2.5e-3])),
                                 ("FP-", np.diag(sigB), piB, ddB)):
        sf = initial_state(P, sig0, 1.70, pi0, finite=True)
        stf = step(P, sf, np.diag(dd))
        at = stf.cache["atilde"]
        Vc = stf.cache["nvec"]
        perm = [int(np.argmax(np.abs(Vc[:, a]))) for a in range(3)]
        errs = []
        for h in (1e-6, 1e-7):
            at_fd = np.zeros((3, 3))
            for b in range(3):
                e_ = np.zeros(3)
                e_[perm[b]] = h
                sp_ = step(P, sf, np.diag(dd + e_))
                sm_ = step(P, sf, np.diag(dd - e_))
                assert sp_.flags["fpattern"] == stf.flags["fpattern"] == sm_.flags["fpattern"]
                at_fd[:, b] = (np.diag(sp_.sigma)[perm] - np.diag(sm_.sigma)[perm]) / (2 * h)
            errs.append(relmax(at, at_fd))
        a4 = tangent_finite(P, stf)
        print(f"      {label}: pattern {stf.flags['fpattern']} (finite {stf.finite}), p = {np.trace(stf.sigma) / 3:.4f}, "
              f"a~_f vs FD h=1e-6/1e-7: {errs[0]:.2e} {errs[1]:.2e}; delta:a~_f/max = {np.abs(at.sum(axis=0)).max() / np.abs(at).max():.1e}; "
              f"(S.34) assembled, |a4| {np.abs(a4).max():.3e}")
        assert stf.flags["fpattern"] == label, (label, stf.flags["fpattern"])
        allok &= ok(f"(k) finite mode {label}: a~^ep_f (S.32f) vs FD of tau_a, h=1e-7", errs[1], 1e-7)
        if label == "FE-":      # p pinned only when the LAST operator applied is an active Pi_f (FE-, -Pf, FPf; not FP-)
            allok &= ok(f"(k) finite mode {label}: delta : a~_f = 0 (p pinned)", np.abs(at.sum(axis=0)).max() / np.abs(at).max(), 1e-13)
    print("  floor group:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ------------------------------------------------------------------ 10. pi_i0 rule (S.53), (S.56), parser refusals
def check_pi0():
    print("\n=== 10. unified pi_i0 rule (S.53, K1.15), the (S.56) refusal, the §2.4 parser refusals ===")
    allok = True
    Pc = k2_params(cap="smooth", c1=0.05, c2=0.15)
    Pn = k2_params()
    v0 = 1.59
    vals = dict(ramp_end=initial_state(Pc, -100 * I3, v0, None).pi_i, apex=initial_state(Pn, -100 * I3, v0, None).pi_i,
                legacy_cap=initial_state(Pc, -100 * I3, v0, None, pi0_rule="legacy").pi_i)
    sig75 = np.diag(-100 * K.ONES + np.array([1, 1, -2]) * 75.0 / 3)                 # TXC, q = 75, eta = 0.75
    vals["eta_0.75"] = initial_state(Pn, sig75, v0, None).pi_i
    sigM = np.diag(-100 * K.ONES + np.array([1, 1, -2]) * 120.0 / 3)                 # eta = M: the image point
    vals["eta_M"] = initial_state(Pn, sigM, v0, None).pi_i
    print(f"  K1.15 K2 set p_init = -100: eta* = c2 M = 0.18 -> {vals['ramp_end']:.6f} (-50.995881); eta* = 0 (apex, cap none) -> "
          f"{vals['apex']:.6f} (-46.475800); legacy rule with the cap -> {vals['legacy_cap']:.6f} (= apex); eta_init = 0.75 -> "
          f"{vals['eta_0.75']:.6f} (-71.554175); eta* = M -> {vals['eta_M']:.6f} (-100)")
    allok &= ok("K1.15 ramp_end", abs(vals["ramp_end"] + 50.995881), 1e-6)
    allok &= ok("K1.15 apex", abs(vals["apex"] + 46.475800), 1e-6)
    allok &= ok("K1.15 eta 0.75", abs(vals["eta_0.75"] + 71.554175), 1e-6)
    allok &= ok("K1.15 eta = M", abs(vals["eta_M"] + 100.0), 1e-9)
    assert vals["legacy_cap"] == vals["apex"]
    st0 = initial_state(Pc, -100 * I3, v0, None)
    sts = triaxial(Pc, st0, "drained", -0.004, 800)
    k = next(i for i, s in enumerate(sts) if s.flags["plastic"])
    w, _ = np.linalg.eigh(sts[k].sigma)
    inv = K.invariants(w)
    print(f"      first yield on the drained TXC path (800 x 5e-6): step {k + 1}, q = {inv.q:.4f} (11.3469), p = {inv.p:.4f} (-103.7823), "
          f"eta = {sts[k].eta:.4f} (0.1093), w = {K.cap_weight(Pc, sts[k].eta)[0]:.3f} (0.337); elastic before: {not sts[k - 1].flags['plastic']}")
    # (S.56)
    print(f"  (S.56) W_ramp: defaults {Pc.W_ramp:.4f} (0.0606); c2 = 0.07 -> {k2_params(cap='smooth', c1=0.05, c2=0.07).W_ramp:.4f} (0.0122, admissible)")
    for bad, label in ((dict(cap="smooth", c1=0.05, c2=0.06), "c2 = 0.06 (0.0061 < 10 PI_SCAN_REL)"),):
        try:
            k2_params(**bad)
            print(f"  (S.56) refusal MISSED for {label}")
            allok = False
        except ValueError as e:
            print(f"  (S.56) refused {label}: {str(e)[:70]}")
    for okset in (dict(cap="planar", c1=0.1, c2=0.1), dict(cap="none"), dict(cap="smooth", c1=0.05, c2=0.07)):
        k2_params(**okset)
    print("  (S.56) planar c1 = c2 = 0.1, none, smooth (0.05, 0.07): accepted (A3: gated to cap = smooth)")
    # parser refusals (§2.4)
    cases = [("HAR with p0 given", dict(p0=-100.0)), ("HAR with alpha0 given", dict(alpha0=0.0)), ("HAR k <= 0", dict(k=0.0)),
             ("HAR n_e = 1", dict(n_e=1.0)), ("HAR n_e < 0", dict(n_e=-0.1)), ("HAR p_a <= 0", dict(p_a=0.0)),
             ("p_min < 0", dict(p_min=-1.0))]
    for label, over in cases:
        try:
            tims_params(**over)
            print(f"  §2.4 refusal MISSED: {label}")
            allok = False
        except ValueError as e:
            print(f"  §2.4 refused {label}: {str(e)[:70]}")
    for label, over in (("BA06 with k given", dict(k=1.0)), ("BA06 with n_e given", dict(n_e=0.5)), ("energy typo", dict(energy="HAR05"))):
        try:
            k2_params(**over)
            print(f"  §2.4 refusal MISSED: {label}")
            allok = False
        except ValueError as e:
            print(f"  §2.4 refused {label}: {str(e)[:70]}")
    try:
        tims_params(k=None, g=None, n_e=None)
        allok = False
    except ValueError as e:
        print(f"  §2.4 refused HAR without k/g/n_e: {str(e)[:60]}")
    try:
        initial_state(Pn, 1.0 * I3, v0, None)
        allok = False
    except ValueError as e:
        print(f"  initial p >= 0 refused: {e}")
    print(f"  defaults: Params() p_min = {Params().p_min} (0.5), energy {Params().energy}, p0 {Params().p0}; "
          f"tims p_min = {tims_params().p_min:.6f} (5e-3 x 101 = 0.505), p_ref = {tims_params().p_ref}")
    print("  pi0 group:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ------------------------------------------------------------------ 11. D2: the initial tangent at a FLOORED initial state
def _rot3(seed=7):
    """a fixed proper rotation (QR of a seeded random matrix): sigma0 with a non-trivial spectral basis."""
    Q, R = np.linalg.qr(np.random.default_rng(seed).standard_normal((3, 3)))
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def _elastic_tangent_like_kernel(P, st):
    """the kernel's elasticTangent(): eigen-decompose the COMMITTED eps^e, a^e (S.3) there, (S.33) with the spin from
    sigma(eps^e) -- rebuilt here from the state alone (no cache), as an independent reference."""
    w, V = np.linalg.eigh(sym(st.eps_e))
    el = K.elastic(P, w)
    return K.tangent_small(el.ae, el.sig, w, V)


def _d2_init_cases():
    """(label, P, deck sigma0, v0, pi_i0) of the floored-init cases (+ one unfloored control)."""
    R = _rot3()
    nh = np.array([-1.0, -0.35, 1.35]); nh -= nh.mean(); nh /= np.linalg.norm(nh)
    PK = k2_params()                                                       # BA06 K2, p_min 0.5
    PT = tims_params(p_min=0.505)                                          # HAR TIMs, p_min 0.505
    sig_k_sh = R @ np.diag(-0.3 + K.SQ23 * 0.2 * nh) @ R.T                # BA06 (p, q) = (-0.3, 0.2), rotated
    sig_t_sh = R @ np.diag(-0.3 + K.SQ23 * 0.5 * nh) @ R.T                # HAR (p, q) = (-0.3, 0.5) (the G1 test case), rotated
    return [("BA06 K2 iso p = -0.2", PK, -0.2 * I3, 1.70, -5000.0),
            ("BA06 K2 (p, q) = (-0.3, 0.2) rotated", PK, sig_k_sh, 1.70, -5000.0),
            ("HAR TIMs iso p = -0.2, pi_i0 = -50", PT, -0.2 * I3, 1.70, -50.0),
            ("HAR TIMs iso p = -0.2, (S.53)", PT, -0.2 * I3, 1.70, None),
            ("HAR TIMs (p, q) = (-0.3, 0.5) rotated", PT, sig_t_sh, 1.70, -5000.0),
            ("control: BA06 K2 p = -100 (not floored)", PK, -100.0 * I3, 1.72, None)]


def _d2_step_cases():
    """(label, P, floored or plain initial state, increment, fractions or None) covering every floor pattern of the
    existing groups, a floored INIT followed by a step, a plain plastic step, and substepped increments."""
    cases = []
    PK = k2_params()
    s0 = initial_state(PK, -100 * I3, 1.59, -2.0e4)
    ev_tr = K.floor_ev(PK, 0.0)[0] + PK.kappa_hat * math.log(2.0)                      # p^tr = -p_min/2 (floor group (a))
    cases.append(("K1.12 FE- (BA06)", PK, s0, ev_tr / 3 * I3 + 1e-6 * np.array([[0, 1, 0.5], [1, 0, 0], [0.5, 0, 0]]), None))
    P = k2_params(p_min=50.0)
    sigB, _, nhB = off_corner_sig(P, -55.0, 1.4)
    sB = initial_state(P, np.diag(sigB), 1.70, K.pi_of_eta(P, -55.0, 1.4))
    for sdev in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3):
        cand = np.diag(1.0e-3 * K.ONES + sdev * nhB) + 1e-4 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        if step(P, sB, cand).flags["fpattern"] == "FP-":
            cases.append(("p50 (B) FP- (BA06)", P, sB, cand, None))
            cases.append(("p50 (D) FP-,FP- m = 2 (BA06)", P, sB, cand, [0.5, 0.5]))
            break
    sigC, _, _ = off_corner_sig(P, -51.0, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
    sC = initial_state(P, np.diag(sigC), 1.70, K.pi_of_eta(P, -51.0, 0.5 * P.M))
    for amp in (2e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1.2e-2):
        cand = np.diag(amp * np.array([0.55, 0.40, -0.95]) - 1e-5) + 1e-4 * np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
        if step(P, sC, cand).flags["fpattern"] == "-Pf":
            cases.append(("p50 (C) -Pf (BA06)", P, sC, cand, None))
            break
    PT = tims_params(p_min=0.505)
    sigF, _, nhF = off_corner_sig(PT, -0.6, 1.2 * PT.M)
    piF = K.pi_of_eta(PT, -0.6, 1.2 * PT.M)
    sF = initial_state(PT, np.diag(sigF), v_for_psi(PT, piF, -0.10), piF)
    dF = np.diag((2e-5 / 3) * K.ONES + 2e-5 * nhF)
    cases.append(("K1.14b FPf (HAR)", PT, sF, dF, None))
    cases.append(("K1.14b FPf,FPf m = 2 (HAR)", PT, sF, 2.0 * dF, [0.5, 0.5]))
    # a FLOORED initial state followed by a step (the kernel_parity FLOOR_HAR_init_above_floor paths)
    si = initial_state(PT, -0.2 * I3, 1.70, -50.0)
    cases.append(("floored init (HAR, p = -0.2) + compressive/shear step", PT, si,
                  -2e-5 / 3 * I3 + 1e-6 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), None))
    si2 = initial_state(PT, -0.2 * I3, 1.70, None)
    cases.append(("floored init (HAR, (S.53)) + hydrostatic step", PT, si2, -2e-5 / 3 * I3, None))
    si3 = initial_state(PK, -0.2 * I3, 1.70, -5000.0)
    cases.append(("floored init (BA06) + expansion step (stays floored, FE-)", PK, si3, 1e-4 / 3 * I3, None))
    # plain plastic step, no floor anywhere (the pre-round-3 model)
    Pp = k2_params()
    sp = initial_state(Pp, -100 * I3, v_for_psi(Pp, -80.0, -0.05), -80.0)
    cases.append(("plain plastic TXC step (BA06 K2)", Pp, sp, np.diag([3e-3, 3e-3, -1e-2]), None))
    return cases


def check_d2():
    print("\n=== 11. D2 (owner decision 2026-10-04): the initial tangent at a FLOORED initial state is the full elastic a^e ===")
    allok = True
    # (a) tangent(initial_state) == a^e at the floored eps^e_f (the kernel's elasticTangent()); a^e Phi_init NOT returned
    print("  (a) tangent() of a fresh initial state vs the full a^e (S.33 on the committed eps^e) and vs a^e Phi_init:")
    init_cases = _d2_init_cases()
    for label, P, sig0, v0, pi0 in init_cases:
        st = initial_state(P, sig0, v0, pi0)
        c = st.cache
        C = tangent(P, st)
        C_cache = K.tangent_small(K.elastic(P, c["eps_tr"]).ae, c["sig"], c["eps_tr"], c["nvec"])   # same spectral data
        C_indep = _elastic_tangent_like_kernel(P, st)                                              # re-diagonalised
        K_C = float(np.einsum("iijj->", C)) / 9.0
        D11 = K.elastic(P, c["eps_tr"]).D11
        print(f"      {label}: n_f_init {st.n_f_init}, p = {np.trace(st.sigma) / 3:.4f}, at_floor {st.flags['at_floor']}; "
              f"bit-identical to a^e: {np.array_equal(C, C_cache)}; vs re-diagonalised {relmax(C, C_indep):.1e}; "
              f"delta:C:delta/9 = {K_C:.6f} = D11(p) {D11:.6f}")
        allok &= ok(f"(a) {label}: tangent == a^e (bit-identical)", float(np.abs(C - C_cache).max()), 0.0)
        allok &= ok(f"(a) {label}: tangent vs a^e re-diagonalised", relmax(C, C_indep), 1e-12)
        allok &= ok(f"(a) {label}: bulk stiffness delta:C:delta/9 = D11 (compressive branch)", abs(K_C - D11) / D11, 1e-12)
        if st.n_f_init:
            assert c["Phi_init"] is not None
            C_phi = K.tangent_small(K.elastic(P, c["eps_tr"]).ae @ c["Phi_init"], c["sig"], c["eps_tr"], c["nvec"])
            gap = relmax(C, C_phi)
            dCphi = np.abs(np.einsum("iikl->kl", C_phi)).max() / np.abs(C_phi).max()
            K_gap = float(np.einsum("iijj->", C - C_phi)) / 9.0          # the discarded bulk: a^e (1 - Phi) 1/9 = D11 exactly
            print(f"          a^e Phi_init would give delta:C/max|C| = {dCphi:.1e} (zero bulk); it is {gap:.3f} away from the returned "
                  f"tangent (max-norm; 2 mu0 dominates under BA06), bulk part of the difference {K_gap:.6f} = D11")
            allok &= ok(f"(a) {label}: a^e Phi_init is NOT returned (bulk of C - a^e Phi = D11)", abs(K_gap - D11) / D11, 1e-12)
            allok &= ok(f"(a) {label}: a^e Phi_init has no bulk stiffness (the discarded operator)", dCphi, 1e-13)
        else:
            assert c["Phi_init"] is None
            allok &= ok(f"(a) {label}: unfloored init, no Phi_init", 0.0, 0.0)
    # (b) every STEP tangent is independent of the initial-state tangent cache: the same step from the same state with the
    #     cache carrying the old a^e Phi_init (or garbage) is bit-identical in sigma, pi_i, counters and tangent
    print("  (b) step tangents do not read the initial tangent cache (bit-identical with the old a^e Phi_init / garbage in it):")
    for label, P, s0, d, fr in _d2_step_cases():
        def run(s):
            return step(P, s, d) if fr is None else step_fractions(P, s, d, fr)
        st1 = run(s0)
        assert not st1.flags["refused"], (label, st1.flags["reason"])
        same = True
        for variant in ("old a^e Phi", "garbage"):
            sm = s0.copy()
            if variant == "old a^e Phi":
                Phi = s0.cache.get("Phi_init")
                if Phi is None:
                    continue
                sm.cache["atilde"] = K.elastic(P, s0.cache["eps_tr"]).ae @ Phi
            else:
                sm.cache["atilde"] = np.full((3, 3), np.nan)
                sm.cache["sig"] = np.full(3, np.nan)
            st2 = run(sm)
            same &= (np.array_equal(st1.sigma, st2.sigma) and st1.pi_i == st2.pi_i and st1.v == st2.v
                     and np.array_equal(tangent(P, st1), tangent(P, st2)) and st1.flags["fpattern"] == st2.flags["fpattern"]
                     and (st1.n_f_tr, st1.n_f_post, st1.eps_f_v, st1.W_f) == (st2.n_f_tr, st2.n_f_post, st2.eps_f_v, st2.W_f))
            # ... and along a 3-increment path
            a_ = run_path(P, s0, np.array([d] * 3)) if fr is None else None
            b_ = run_path(P, sm, np.array([d] * 3)) if fr is None else None
            if a_ is not None:
                same &= all(np.array_equal(x.sigma, y.sigma) and np.array_equal(tangent(P, x), tangent(P, y)) and x.pi_i == y.pi_i
                            for x, y in zip(a_, b_))
        print(f"      {label}: pattern '{st1.flags['fpattern']}', substeps {st1.flags['substeps']}, init floored {bool(s0.n_f_init)}: "
              f"bit-identical {same}")
        allok &= ok(f"(b) {label}: step tangent / state unchanged", 0.0 if same else 1.0, 0.0)
    # (c) the C++ kernel's elasticTangent() at the SAME floored initial states (kernel_parity shim, optional: needs g++ and
    #     the kernel header next to the repo; skipped with a note elsewhere)
    print("  (c) kernel elasticTangent() at the floored initial states vs O2 tangent() (kernel_parity shim):")
    kern, NK, note = _kernel_shim()
    if kern is None:
        print(f"      SKIPPED: {note}")
    else:
        worst = 0.0
        for label, P, sig0, v0, pi0 in init_cases:
            st = initial_state(P, sig0, v0, pi0)
            rc, ks, msg = kern.initial_state(P, sig0, v0, pi0)
            assert rc == 0, (label, rc, msg)
            Ck = kern.elastic_tangent(P, ks)
            Co = NK.c4_to_c6(tangent(P, st))
            err = float(np.abs(Ck - Co).max() / np.abs(Co).max())
            worst = max(worst, err)
            gap = ""
            if st.n_f_init:
                c = st.cache
                Cphi = NK.c4_to_c6(K.tangent_small(K.elastic(P, c["eps_tr"]).ae @ c["Phi_init"], c["sig"], c["eps_tr"], c["nvec"]))
                gap = f"; kernel vs the old O2 a^e Phi_init: {np.abs(Ck - Cphi).max() / np.abs(Ck).max():.3f}"
            print(f"      {label}: kernel n_f_init {int(ks[16])}, |C_k - C_o2|/max|C_o2| = {err:.2e}{gap}")
            assert int(ks[16]) == st.n_f_init
        allok &= ok("(c) kernel elasticTangent() == O2 tangent() at every floored initial state (gate 1e-10)", worst, 1e-10)
    print("  d2 group:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ------------------------------------------------------------------ 12. round 4: the image-pressure rule at a post floor
# Sheet §9.1 step 5g, §9.7 (S.57)-(S.59), §13.12 (C), §9.7 item 3 (wet-floor table), §13.16 (K1.16). Every expected
# number below is copied from the sheet (round 4b, 6ae9fe5ab) and was written here BEFORE the first run.
def k3_params(**kw):
    """K3 set of sheet §13.16 (K1.16): HAR, WW, fork CSL, smooth cap, p_min 0.505 kPa."""
    d = dict(energy="HAR", k=298.34244992965336, g=383.5831499095543, n_e=0.5, p_a=101.0, M=1.25,
             N=0.4686136860640983, N_bar=0.16034970350694513, rho=0.5828981510111576, rho_bar=0.7213710496292179,
             zeta="WW", chi=-2.8923322243122183, h=308.3223484879116, csl_mode="fork", e0=0.934, lam_c=0.019, xi=0.7,
             cap="smooth", c1=0.05, c2=0.15, p_min=0.505)
    d.update(kw)
    return Params(**d).validate()


def c_N(P):
    """c_N = (1-N)^((1-N)/N) (N > 0), e^-1 (N = 0): |pi_i|/|p| at the tip eta = 0 (S.58)."""
    return (1.0 - P.N) ** ((1.0 - P.N) / P.N) if P.N > 0.0 else math.exp(-1.0)


def F_of(P, sigma, pi):
    """F(sigma, pi_i) of a 3x3 stress, as the trial test evaluates it."""
    return K.flow(P, K.invariants(np.linalg.eigvalsh(sym(sigma))), pi).F


def pi_f_of(P, sigma):
    """pi_f(sigma): the inverse of (S.12) at (p, eta = zeta q/|p|) (S.57), recomputed outside the kernel path."""
    inv = K.invariants(np.linalg.eigvalsh(sym(sigma)))
    eta = 0.0 if inv.vertex else K.zeta_y(inv.theta, P.rho, P.zeta)[0] * inv.q / abs(inv.p)
    return K.pi_of_eta(P, inv.p, eta), eta


class rule_off:
    """Context: the round-3b return map (step 5g absent: pi_i := pi_i^BE whatever F(sigma_f)); the M-R4 mutant."""

    def __enter__(self):
        self._orig = K.floor_pi_rule

        def off(P, sig_f, ae_f, pi_be):
            return dict(pi=pi_be, pi_be=pi_be)
        K.floor_pi_rule = off

    def __exit__(self, *a):
        K.floor_pi_rule = self._orig


H19_DIR = None


def h19_records():
    """The dumped H19 Gauss-point states: out/h19/s185.txt and stall.txt. Line format (header of h19_replay.cpp):
    id ee11 ee22 eg12 s11 s22 s12 s33 pi_i v v0 epv eps at_floor n_f_tr n_f_post eps_f_v W_f d11 d22 dg12
    (plane strain; eg12 and dg12 are engineering shears; eps^e_33 is not dumped and is recovered from s33)."""
    import os
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                        "norsand_calib", "out", "h19")
    recs = {}
    for fn in ("s185.txt", "stall.txt"):
        with open(os.path.join(base, fn)) as fh:
            for line in fh:
                t = line.split()
                if not t or t[0].startswith("#"):
                    continue
                x = [float(a) for a in t[1:]]
                recs[t[0]] = dict(ee=x[0:3], s=x[3:7], pi=x[7], v=x[8], v0=x[9], epv=x[10], eps=x[11], at_floor=int(x[12]),
                                  ntr=int(x[13]), npost=int(x[14]), efv=x[15], wf=x[16], d=x[17:20])
    return recs


def h19_state(P, rec, norm=2e-4):
    """(State, increment tensor scaled to |d eps| = (d11^2 + d22^2 + dg12^2/2)^(1/2) = norm, eps^e_33, stress
    reproduction error). eps^e_33 solves sigma_33(eps^e) = s33 (secant on O2's own elastic law)."""
    e11, e22, g12 = rec["ee"]

    def sig_of(e33):
        E = np.array([[e11, 0.5 * g12, 0.0], [0.5 * g12, e22, 0.0], [0.0, 0.0, e33]])
        w, V = np.linalg.eigh(E)
        return E, (V * K.elastic(P, w).sig) @ V.T
    a, b = 0.5 * (e11 + e22), 0.5 * (e11 + e22) + 1e-5
    fa, fb = sig_of(a)[1][2, 2] - rec["s"][3], sig_of(b)[1][2, 2] - rec["s"][3]
    for _ in range(60):
        if fb == fa:
            break
        c = b - fb * (b - a) / (fb - fa)
        a, fa = b, fb
        b, fb = c, sig_of(c)[1][2, 2] - rec["s"][3]
        if abs(b - a) <= 1e-17:
            break
    E, sig = sig_of(b)
    err = max(abs(sig[0, 0] - rec["s"][0]), abs(sig[1, 1] - rec["s"][1]), abs(sig[0, 1] - rec["s"][2]),
              abs(sig[2, 2] - rec["s"][3]))
    st = State(sig, E, rec["pi"], rec["v"], v0=rec["v0"], eps_p_v=rec["epv"], eps_p_s=rec["eps"],
               eps_f_v=rec["efv"], W_f=rec["wf"], n_f_tr=rec["ntr"], n_f_post=rec["npost"])
    st.flags["at_floor"] = bool(rec["at_floor"])
    d11, d22, dg = rec["d"]
    sc = norm / math.sqrt(d11 * d11 + d22 * d22 + 0.5 * dg * dg)
    deps = sc * np.array([[d11, 0.5 * dg, 0.0], [0.5 * dg, d22, 0.0], [0.0, 0.0, 0.0]])
    return st, deps, b, err


def s4(st):
    return np.array([st.sigma[0, 0], st.sigma[1, 1], st.sigma[2, 2], st.sigma[0, 1]])


# ---- sheet §13.12 (C) with the round-4 rule
K112C = dict(pi0=-33.125471695, theta=0.454427, dlam=1.546466e-3, p_c=-48.764495, q_c=46.421589, eta_c=0.952607,
             pi_be=-40.198569107, sig_f=(-34.8599764962, -34.1946812138, -80.9453422900), eta_f=0.929068,
             F_over_pmin=0.010899, pi=-40.516314098, cN_pmin=23.2379)


def k112c_case():
    P = k2_params(p_min=50.0)
    sig, th, nh = off_corner_sig(P, -51.0, 0.6, direction=(-1.2, -0.1, 1.3))
    pi0 = K.pi_of_eta(P, -51.0, 0.6)
    st0 = initial_state(P, np.diag(sig), 1.70, pi0)
    deps = np.diag(6e-3 * np.array([0.55, 0.40, -0.95]) - 1e-5)
    return P, st0, deps, th


def check_r4_k112c():
    print("\n=== 12a. K1.12 (C) with the image-pressure rule (S.57) (sheet §13.12; BA06 K2 set, p_min = 50 kPa) ===")
    allok = True
    X = K112C
    P, st0, deps, th = k112c_case()
    st = step(P, st0, deps)
    res = st.cache["res"]
    elc = K.elastic(P, res.eps_c)
    sf = np.diag(st.sigma)
    print(f"  start pi_i = {st0.pi_i:.9f} ({X['pi0']}), theta = {th:.6f} ({X['theta']}); pattern {st.flags['rpattern']} (-Pf*), "
          f"substeps {st.flags['substeps']}")
    print(f"  return: dlam = {res.dlam:.6e} ({X['dlam']:.6e}), p_c = {elc.p:.6f} ({X['p_c']}), q_c = {elc.q:.6f} ({X['q_c']}), "
          f"eta_c = {res.eta:.6f} ({X['eta_c']}), pi_BE = {res.pi_be:.9f} ({X['pi_be']})")
    print(f"  post floor: sigma_f = ({sf[0]:.10f}, {sf[1]:.10f}, {sf[2]:.10f}) (sheet {X['sig_f']}), eta_f = {res.eta_f:.6f} "
          f"({X['eta_f']}), F(sigma_f, pi_BE)/p_min = {res.F_post / P.p_min:+.6f} (+{X['F_over_pmin']})")
    Fc = F_of(P, st.sigma, st.pi_i)
    print(f"  committed pi_i = {st.pi_i:.9f} (sheet {X['pi']}; round 3b {X['pi_be']}), F(sigma_f, pi_i)/p_min = {Fc / P.p_min:.1e} "
          f"(0 to 3e-16), n_f_pi = {st.n_f_pi}, floor = {st.floor}; c_N p_min = {c_N(P) * P.p_min:.4f} ({X['cN_pmin']})")
    assert st.flags["rpattern"] == "-Pf*" and st.flags["fpattern"] == "-Pf" and st.flags["substeps"] == 1
    assert (st.flags["floor_pi"], st.flags["floor_pi_guard"], st.flags["floor_pi_band"], st.n_f_pi) == (1, 0, 0, 1)
    allok &= ok("K1.12 (C) start pi_i", abs(st0.pi_i - X["pi0"]), 2e-9)
    allok &= ok("K1.12 (C) dlam", abs(res.dlam - X["dlam"]), 1e-9)
    allok &= ok("K1.12 (C) p_c", abs(elc.p - X["p_c"]), 1e-6)
    allok &= ok("K1.12 (C) q_c", abs(elc.q - X["q_c"]), 1e-6)
    allok &= ok("K1.12 (C) eta_c", abs(res.eta - X["eta_c"]), 1e-6)
    allok &= ok("K1.12 (C) pi_i^BE", abs(res.pi_be - X["pi_be"]), 2e-9)
    allok &= ok("K1.12 (C) sigma_f", float(np.abs(sf - np.array(X["sig_f"])).max()), 2e-10)
    allok &= ok("K1.12 (C) eta_f", abs(res.eta_f - X["eta_f"]), 1e-6)
    allok &= ok("K1.12 (C) F(sigma_f, pi_BE)/p_min", abs(res.F_post / P.p_min - X["F_over_pmin"]), 1e-6)
    allok &= ok("K1.12 (C) NEW pi_i = -40.516314098", abs(st.pi_i - X["pi"]), 2e-9)
    allok &= ok("K1.12 (C) F(sigma_f, pi_i)/p_min = 0", abs(Fc) / P.p_min, 1e-13)
    allok &= ok("K1.12 (C) c_N p_min", abs(c_N(P) * P.p_min - X["cN_pmin"]), 1e-4)
    # the rule touches pi_i alone: everything else of the step is the round-3b step, bit for bit (M-R9 included)
    with rule_off():
        s3 = step(P, st0, deps)
    same = (np.array_equal(s3.sigma, st.sigma) and np.array_equal(s3.eps_e, st.eps_e) and s3.v == st.v and s3.D == st.D
            and s3.dlam == st.dlam and s3.eps_f_v == st.eps_f_v and s3.W_f == st.W_f and s3.n_f_post == st.n_f_post
            and np.array_equal(tangent(P, s3), tangent(P, st)) and np.array_equal(s3.eps_p, st.eps_p)
            and floor_energy(P, s3) == floor_energy(P, st) and s3.psi_i == st.psi_i and s3.pi_star == st.pi_star
            and s3.eta == st.eta)
    print(f"  rule off (round 3b): pi_i = {s3.pi_i:.9f}; sigma_f, eps^e, v, dlam, D, eps^p, d eps^f_v, W_f, E_f, psi_i, pi*, eta and "
          f"the (S.32f) tangent bit-identical with the rule on: {same}")
    assert same and s3.pi_i == res.pi_be
    # cases (A), (B) of K1.12 at p_min = 50 do not change (no post floor): see `floor` group (d)
    print("  r4 K1.12 (C):", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ---- sheet §9.7 item 3, the wet-floor table with the rule: psi_i -> (F/p_min, pi_BE, eta_f, pi_i) at d gamma 1e-4/1e-3/1e-2
WET = {0.0: ((0.2699, 0.3496, 0.3572), (-0.328889, -0.329228, -0.329262), (0.940109, 1.021599, 1.029388),
             (-0.386215057, -0.406870322, -0.408940969)),
       0.04: ((0.2717, 0.3670, 0.5322), (-0.328346, -0.324406, -0.287273), (0.938960, 1.012743, 0.951647),
              (-0.385936497, -0.404537473, -0.389032266)),
       0.08: ((0.2732, 0.3820, 0.6866), (-0.327887, -0.320380, -0.258876), (0.937988, 1.005324, 0.897090),
              (-0.385701048, -0.402600177, -0.376005318)),
       0.12: ((0.2745, 0.3950, 0.8217), (-0.327494, -0.316976, -0.238442), (0.937157, 0.999032, 0.856623),
              (-0.385499896, -0.400969027, -0.366805331))}
WET_START = dict(q=0.2583, theta=0.4544, pi=-0.328007)


def check_r4_wet():
    print("\n=== 12b. wet-floor table with the rule (sheet §9.7 item 3; HAR TIMs set, p_min 0.505, eta = 0.5M, pure shear) ===")
    allok = True
    P = tims_params()
    sig, th, nh = off_corner_sig(P, -0.505, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
    pi0 = K.pi_of_eta(P, -0.505, 0.5 * P.M)
    q0 = K.invariants(sig).q
    print(f"  start: p = -0.505, q = {q0:.4f} ({WET_START['q']}), theta = {th:.4f} ({WET_START['theta']}), pi_i = {pi0:.6f} "
          f"({WET_START['pi']})")
    allok &= ok("wet table start (q, theta, pi_i)", max(abs(q0 - WET_START["q"]), abs(th - WET_START["theta"]),
                                                       abs(pi0 - WET_START["pi"]) * 100), 1e-4)
    eF = eB = eE = eP = eZ = 0.0
    for psi, (Fs, Bs, Es, Ps) in WET.items():
        st0 = initial_state(P, np.diag(sig), v_for_psi(P, pi0, psi), pi0)
        for j, dg in enumerate((1e-4, 1e-3, 1e-2)):
            st = step(P, st0, np.diag(dg / math.sqrt(2.0) * nh))
            res = st.cache["res"]
            Fc = F_of(P, st.sigma, st.pi_i) / P.p_min
            print(f"  psi_i {psi:+.2f} dgamma {dg:.0e}: {st.flags['rpattern']} m={st.flags['substeps']}  F/p_min {res.F_post / P.p_min:+.4f} "
                  f"({Fs[j]:+.4f})  pi_BE {res.pi_be:.6f} ({Bs[j]:.6f})  eta_f {res.eta_f:.6f} ({Es[j]:.6f})  pi_i {st.pi_i:.9f} "
                  f"({Ps[j]:.9f})  F after {Fc:.1e}")
            assert st.flags["rpattern"] == "-Pf*" and st.flags["substeps"] == 1 and st.n_f_pi == 1, st.flags
            eF = max(eF, abs(res.F_post / P.p_min - Fs[j]))
            eB = max(eB, abs(res.pi_be - Bs[j]))
            eE = max(eE, abs(res.eta_f - Es[j]))
            eP = max(eP, abs(st.pi_i - Ps[j]))
            eZ = max(eZ, abs(Fc))
    allok &= ok("wet table F(sigma_f, pi_BE)/p_min (12 cases)", eF, 1e-4)
    allok &= ok("wet table pi_i^BE [kPa]", eB, 1e-6)
    allok &= ok("wet table eta_f", eE, 1e-6)
    allok &= ok("wet table pi_i by (S.57) [kPa]", eP, 2e-9)
    allok &= ok("wet table F(sigma_f, pi_i)/p_min = 0 after the rule (sheet 7e-16)", eZ, 1e-13)
    print("  r4 wet-floor table:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ---- sheet §13.16 (K1.16): rows (steps -> s11, s22, s33, s12, pi_i), eps^e_33, the scaled increment, step-1 detail
K116 = {
    "e7367g1_s185": dict(
        e33=2.099254865382e-3, inc=(4.965401387252e-5, -6.100176746158e-6, -2.738512981529e-4), pi0=-0.3224962955, F0=0.2315,
        rows={1: (-0.5101907038, -0.5976456953, -0.4071636009, -0.2151586912, -0.3989652918),
              2: (-0.5082384097, -0.6104039680, -0.3963576223, -0.2510036052, -0.4461726222),
              60: (-0.5114970385, -0.6473686614, -0.3561343001, -0.3336845167, -0.5923861407)},
        dlam=1.427733871e-4, p_c=-0.365744965, eta_c=1.092426, pi_be=-0.3245384951, F_post=0.175979368, eta_f=0.922543205,
        dfv=7.223598e-5),
    "e7367g2_s185": dict(
        e33=2.100404994193e-3, inc=(6.298503508381e-5, -1.654199187512e-5, -2.674294219427e-4), pi0=-0.3666004637, F0=0.1472,
        rows={1: (-0.4889744155, -0.6273143992, -0.3987111852, -0.2332250428, -0.4276028897),
              2: (-0.4850426738, -0.6393433615, -0.3906139647, -0.2596017804, -0.4672749361),
              60: (-0.4814251232, -0.6764194559, -0.3571554210, -0.3278584883, -0.5976134600)},
        dlam=1.401240162e-4, p_c=-0.388932212, eta_c=1.177716, pi_be=-0.3676032012, F_post=0.118226999, eta_f=1.026016850,
        dfv=5.957758e-5),
    "e7367g1_m53": dict(
        e33=2.088911465473e-3, inc=(4.980341391149e-5, -6.390050089011e-6, -2.737838096843e-4), pi0=-0.1116066186, F0=1.6596,
        rows={1: (-0.5071972881, -0.5355837599, -0.4722189520, -0.0691517825, -0.2815026983),
              2: (-0.5031271070, -0.5765889736, -0.4352839194, -0.1789591553, -0.3671148499)},
        dlam=1.528227673e-4, p_c=-0.187753246, eta_c=0.484384, pi_be=-0.1150543895, F_post=1.439692306, eta_f=0.294256201,
        dfv=1.856529e-4, q_c=0.080655979)}
K116_CN_PMIN = 0.246558167
# round 3b (rule off) after 60 steps from e7367g1_s185: s11, s22, s33, s12, pi_i, F [kPa]
K116_R3B_60 = (-0.5098218813, -0.6109439554, -0.3942341633, -0.2483437655, -0.3894275650, 0.0929)


def check_r4_k116():
    print("\n=== 12c. K1.16: the rule on the H19 states (sheet §13.16; K3 HAR set, |d eps| = 2e-4 repeated) ===")
    allok = True
    P = k3_params()
    F_tol = K.F_TRIAL_TOL_REL * P.p_ref
    recs = h19_records()
    print(f"  c_N p_min = {c_N(P) * P.p_min:.9f} ({K116_CN_PMIN})")
    allok &= ok("c_N p_min", abs(c_N(P) * P.p_min - K116_CN_PMIN), 2e-9)
    for sid, X in K116.items():
        st, deps, e33, err = h19_state(P, recs[sid])
        F0 = F_of(P, st.sigma, st.pi_i)
        inc = (deps[0, 0], deps[1, 1], 2.0 * deps[0, 1])
        print(f"  {sid}: dumped stresses reproduced to {err:.1e} kPa; eps^e_33 = {e33:.12e} ({X['e33']:.12e}); pi_i = {st.pi_i:.10f} "
              f"({X['pi0']}), F = {F0:+.4f} kPa = {F0 / P.p_min:+.3f} p_min (sheet {X['F0']:+.4f}); legacy start, F > F_tol")
        allok &= ok(f"{sid} stress reproduction [kPa]", err, 1e-12)
        allok &= ok(f"{sid} eps^e_33", abs(e33 - X["e33"]), 2e-15)
        allok &= ok(f"{sid} scaled increment (d11, d22, dgamma12)", float(np.abs(np.array(inc) - np.array(X["inc"])).max()), 2e-16)
        allok &= ok(f"{sid} start F [kPa]", abs(F0 - X["F0"]), 1e-4)
        nmax = max(X["rows"])
        stall = sid.endswith("m53")
        cur, all_m1 = st, True
        for n in range(1, nmax + 1):
            prev = cur
            cur = step(P, cur, deps)
            assert not cur.flags["refused"], (sid, n, cur.flags["reason"])
            res = cur.cache["res"]
            m = cur.flags["substeps"]
            all_m1 &= (m == 1)
            Fc = F_of(P, cur.sigma, cur.pi_i)
            p = float(np.trace(cur.sigma)) / 3.0
            pif, _ = pi_f_of(P, cur.sigma)
            # (S.58) and the K1.16 assertions at EVERY accepted step
            assert cur.flags["floor_post"] >= 1 and cur.flags["floor_pi"] >= 1 and res.floor_post and res.floor_pi, (sid, n, cur.flags)
            assert cur.flags["floor_pi_guard"] == 0 and cur.flags["floor_pi_band"] == 0
            assert Fc <= F_tol and abs(Fc) <= 1e-13 * P.p_min * 100, (sid, n, Fc)
            assert abs(cur.pi_i) >= c_N(P) * P.p_min and abs(p + P.p_min) <= 1e-12 * P.p_min
            assert abs(cur.pi_i - pif) <= 1e-12 * abs(pif) and cur.pi_i < res.pi_be
            assert cur.n_f_pi == prev.n_f_pi + cur.flags["floor_pi"]
            if n == 1:
                elc = K.elastic(P, np.linalg.eigvalsh(sym(_from_pr(res.eps_c, cur.cache["nvec"]))))
                print(f"    step 1 (m = {m}, {cur.flags['rpattern']}, local iters {cur.flags['local_iters']}): dlam {res.dlam:.9e} ({X['dlam']:.9e}), "
                      f"p_c {elc.p:.9f} ({X['p_c']}), eta_c {res.eta:.6f} ({X['eta_c']}), pi_BE {res.pi_be:.10f} ({X['pi_be']}), "
                      f"F(sigma_f, pi_BE) {res.F_post:+.9f} ({X['F_post']:+.9f}) kPa, eta_f {res.eta_f:.9f} ({X['eta_f']}), "
                      f"d eps^f_v {res.dfv_post:.6e} ({X['dfv']:.6e})")
                if m == 1:
                    d1 = max(abs(res.dlam - X["dlam"]) / X["dlam"] * 1e-1, abs(elc.p - X["p_c"]), abs(res.eta - X["eta_c"]) * 1e-3,
                             abs(res.pi_be - X["pi_be"]), abs(res.F_post - X["F_post"]), abs(res.eta_f - X["eta_f"]),
                             abs(res.dfv_post - X["dfv"]) * 1e1,
                             abs(elc.q - X.get("q_c", elc.q)))
                    allok &= ok(f"{sid} step-1 detail (dlam, p_c, eta_c, pi_BE, F, eta_f, d eps^f_v)", d1, 2e-9)
                    # the stresses of step 1 are those of round 3b; round 3b commits pi_BE with F = F_post
                    with rule_off():
                        s3 = step(P, prev, deps)
                    assert np.array_equal(s3.sigma, cur.sigma) and s3.pi_i == res.pi_be and s3.D == cur.D
                    print(f"    step 1 rule off (round 3b, M-R4): same sigma bit for bit, pi_i = {s3.pi_i:.10f}, "
                          f"F = {F_of(P, s3.sigma, s3.pi_i):+.4f} kPa")
            if n in X["rows"]:
                got = np.append(s4(cur), cur.pi_i)
                e = float(np.abs(got - np.array(X["rows"][n])).max())
                print(f"    step {n:2d}: m = {m} {cur.flags['rpattern']}  s11 {got[0]:.10f} s22 {got[1]:.10f} s33 {got[2]:.10f} s12 {got[3]:.10f} "
                      f"pi_i {got[4]:.10f}  |F| = {abs(Fc):.1e} kPa  q = {K.invariants(np.linalg.eigvalsh(cur.sigma)).q:.3f}")
                print(f"             sheet      s11 {X['rows'][n][0]:.10f} s22 {X['rows'][n][1]:.10f} s33 {X['rows'][n][2]:.10f} "
                      f"s12 {X['rows'][n][3]:.10f} pi_i {X['rows'][n][4]:.10f}   max diff {e:.1e}")
                if stall and not all_m1:
                    print("             (stall row at m > 1: asserted by the invariants only, sheet §13.16 round 4b)")
                else:
                    assert all_m1, (sid, n, "the sheet digits need m = 1 at every step")
                    allok &= ok(f"{sid} step {n} (sigma, pi_i) [kPa]", e, 2e-10)
        print(f"    {nmax} steps: every step post-floored with the rule fired, F <= F_tol, pi_i = pi_f(sigma_f), p = -p_min, "
              f"|pi_i| >= c_N p_min; m = 1 throughout: {all_m1}; n_f_pi {st.n_f_pi} -> {cur.n_f_pi}, n_f_post {st.n_f_post} -> {cur.n_f_post}")
        assert cur.n_f_pi == nmax or not all_m1
        if sid == "e7367g1_s185":
            with rule_off():
                c3 = st
                for n in range(60):
                    c3 = step(P, c3, deps)
            got = np.append(np.append(s4(c3), c3.pi_i), F_of(P, c3.sigma, c3.pi_i))
            print(f"    rule off, 60 steps (round 3b): s11 {got[0]:.10f} s22 {got[1]:.10f} s33 {got[2]:.10f} s12 {got[3]:.10f} pi_i {got[4]:.10f} "
                  f"F {got[5]:+.4f}; sheet {K116_R3B_60}; q = {K.invariants(np.linalg.eigvalsh(c3.sigma)).q:.3f} (0.469; 0.631 with the rule)")
            allok &= ok("round-3b 60-step values with the rule off (sigma, pi_i)", float(np.abs(got[:5] - np.array(K116_R3B_60[:5])).max()), 2e-10)
            allok &= ok("round-3b 60-step F [kPa]", abs(got[5] - K116_R3B_60[5]), 1e-4)
    # a non-converging return from the stall state is refused, state frozen, nothing counted (never F > F_tol committed)
    st, deps, _, _ = h19_state(P, recs["e7367g1_m53"], norm=1e-6)
    r = step(P, st, deps)
    print(f"  stall state at |d eps| = 1e-6: refused {r.flags['refused']} '{r.flags['reason']}' (sheet: refused, local_linesearch:pi_fold), "
          f"pi_i {r.pi_i:.10f} (-0.1116066186), n_f_pi {r.n_f_pi}, n_f_post {r.n_f_post} (unchanged {st.n_f_post})")
    if r.flags["refused"]:
        assert r.pi_i == st.pi_i and (r.n_f_pi, r.n_f_post, r.n_f_tr, r.eps_f_v) == (st.n_f_pi, st.n_f_post, st.n_f_tr, st.eps_f_v)
    else:
        assert F_of(P, r.sigma, r.pi_i) <= F_tol
    print("  r4 K1.16:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


def _from_pr(vals, V):
    return (V * vals) @ V.T


# ---- (S.59): the chain row, its gate fixture (case (b) of the sheet's FD table) and the named mutants
S59_GATE = dict(pi_start=-0.3989652918, F_over_pmin=(0.162, 0.099), pi=-0.479306004, C_fd=(1.0e-5, 1.0e-7, 1.3e-9),
                Spi_fd=(4.7e-6, 4.7e-8, 9.4e-10), gate_h7=3e-7,
                mutants={"M-R1": 0.20, "M-R2": 0.19, "M-R3": 0.020, "M-R5": 0.016, "M-R6": 0.26})


def r4_chain(P, st, deps, fractions, mutant=None):
    """The increment with prescribed fractions through api._step_once and kernel.chain_propagate, with one named
    mutant of the (S.59) row applied on every sub-increment where the rule fired (mutant None = the shipped chain).
    Returns (C (S.47), S^pi (6,), final state, [StepResult per sub-increment])."""
    from dataclasses import replace
    from .api import _step_once
    S_eps, S_pi, cum = K.chain_start()
    cur, rs = st, []
    for a in fractions:
        cur = _step_once(P, cur, a * deps)
        assert not cur.flags["refused"], cur.flags["reason"]
        res, nv = cur.cache["res"], cur.cache["nvec"]
        rs.append(res)
        if mutant is None or not res.floor_pi:
            S_eps, S_pi, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, res, nv)
            continue
        inv = K.invariants(res.sig)
        fl = K.flow(P, inv, res.pi)
        G = -fl.f_a / fl.Y.F_pi
        if mutant == "M-R1":        # row not replaced: the (S.54) S^pi row kept
            S_eps, S_pi, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, replace(res, floor_pi=False), nv)
        elif mutant == "M-R2":      # zeta_y q y_a term dropped from f_a
            Gm = -(fl.Y.F_p / 3.0 * K.ONES + K.SQ32 * fl.zeta * inv.nhat) / fl.Y.F_pi
            S_eps, S_pi, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, replace(res, gamma_pi=Gm @ res.ae), nv)
        elif mutant == "M-R3":      # (S.59) row ADDED to the (S.54) row
            _, S54, _ = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, replace(res, floor_pi=False), nv)
            S_eps, S59, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, res, nv)
            S_pi = S54 + S59
        elif mutant == "M-R5":      # a^e of the PRE-floor converged state in gamma_b
            S_eps, S_pi, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v,
                                                 replace(res, gamma_pi=G @ K.elastic(P, res.eps_c).ae), nv)
        elif mutant == "M-R6":      # the row built from the strain sensitivity BEFORE Phi^post
            Spre, _, _ = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, replace(res, floor_post=False, floor_pi=False), nv)
            S_eps, _, cum = K.chain_propagate(S_eps, S_pi, cum, a, cur.v, res, nv)
            S_pi = np.einsum("ia,Jij,ja->Ja", nv, Spre, nv) @ res.gamma_pi
        else:
            raise ValueError(mutant)
    C = K.chain_assemble(P, cur.cache["res"], cur.cache["nvec"], S_eps)
    return C, S_pi, cur, rs


def r4_fd(P, st, deps, fractions, h):
    """Central FD over the six kernel columns of the final stress and of the final pi_i; fractions fixed; the SAME
    branch pattern, floor pattern and reset flags asserted at every FD point; no FD point in the noise band."""
    ref = step_fractions(P, st, deps, fractions)
    fdS, fdP = [], np.empty(6)
    for J in range(6):
        sp = step_fractions(P, st, deps + h * K.CHAIN_E[J], fractions)
        sm = step_fractions(P, st, deps - h * K.CHAIN_E[J], fractions)
        for s_ in (sp, sm):
            assert not s_.flags["refused"] and s_.flags["rpattern"] == ref.flags["rpattern"], (J, s_.flags)
            assert s_.flags["pattern"] == ref.flags["pattern"] and s_.flags["floor_pi_band"] == 0
        fdS.append((sp.sigma - sm.sigma) / (2.0 * h))
        fdP[J] = (sp.pi_i - sm.pi_i) / (2.0 * h)
    return np.array(fdS), fdP


def _c_err(C, fdS):
    return max(np.linalg.norm(np.einsum("ijkl,kl->ij", C, K.CHAIN_E[J]) - fdS[J]) /
               np.linalg.norm(np.einsum("ijkl,kl->ij", C, K.CHAIN_E[J])) for J in range(6))


def check_r4_s59():
    print("\n=== 12d. (S.59) chain row: the forced (1/2, 1/2) HAR gate, FD at h = 1e-6/1e-7/1e-8, mutants M-R1..M-R6 ===")
    allok = True
    X = S59_GATE
    P = k3_params()
    st0, d1, _, _ = h19_state(P, h19_records()["e7367g1_s185"])
    st1 = step(P, st0, d1)                                   # ONE K1.16 step
    deps, fr = 2.0 * d1, (0.5, 0.5)
    stn = step_fractions(P, st1, deps, fr)
    C, S_pi, cur, rs = r4_chain(P, st1, deps, fr)
    Fp = [r.F_post / P.p_min for r in rs]
    print(f"  fixture: start pi_i = {st1.pi_i:.10f} ({X['pi_start']}), |d eps| = {math.sqrt(deps[0, 0] ** 2 + deps[1, 1] ** 2 + 2 * deps[0, 1] ** 2):.1e}, "
          f"fractions {fr}, pattern {stn.flags['rpattern']} (-Pf*,-Pf*); F(sigma_f, pi_BE)/p_min = {Fp[0]:+.3f}, {Fp[1]:+.3f} "
          f"(+{X['F_over_pmin'][0]}, +{X['F_over_pmin'][1]}); band limit F_tol/p_min = {K.F_BAND_TOL_REL * P.p_ref / P.p_min:.1e}")
    print(f"  committed pi_i = {stn.pi_i:.9f} ({X['pi']}), p = {np.trace(stn.sigma) / 3:.12f} (-0.505); floor_pi = {stn.flags['floor_pi']}, "
          f"n_f_pi {st1.n_f_pi} -> {stn.n_f_pi}")
    assert stn.flags["rpattern"] == "-Pf*,-Pf*" and stn.flags["floor_pi"] == 2 and stn.flags["floor_pi_band"] == 0
    assert stn.n_f_pi == st1.n_f_pi + 2 and rs[0].floor_pi, "the rule must fire in the NON-FINAL sub-increment"
    allok &= ok("gate start pi_i", abs(st1.pi_i - X["pi_start"]), 2e-10)
    allok &= ok("gate F(sigma_f, pi_BE)/p_min of the two sub-increments", max(abs(Fp[0] - X["F_over_pmin"][0]), abs(Fp[1] - X["F_over_pmin"][1])), 1e-3)
    allok &= ok("gate committed pi_i", abs(stn.pi_i - X["pi"]), 2e-9)
    allok &= ok("gate committed p = -p_min", abs(np.trace(stn.sigma) / 3 + P.p_min) / P.p_min, 1e-12)
    same = np.array_equal(C, tangent(P, stn)) and np.array_equal(S_pi, stn.cache["S_pi_chain"]) and cur.pi_i == stn.pi_i
    print(f"  selfcheck chain driver (mutant None) == api.step_fractions chain, bit for bit: {same}")
    assert same
    tolC = (3e-5, X["gate_h7"], 3e-8)
    tolP = (3e-5, 3e-7, 3e-8)
    fd7 = None
    for i, h in enumerate((1e-6, 1e-7, 1e-8)):
        fdS, fdP = r4_fd(P, st1, deps, fr, h)
        eC = _c_err(C, fdS)
        eP = float(np.abs(S_pi - fdP).max() / np.abs(S_pi).max())
        print(f"  h = {h:.0e}: chained C (S.47) vs FD {eC:.2e} (sheet {X['C_fd'][i]:.1e})   S^pi (S.59) vs FD {eP:.2e} (sheet {X['Spi_fd'][i]:.1e})")
        allok &= ok(f"(S.59) gate: C vs FD at h = {h:.0e}", eC, tolC[i])
        allok &= ok(f"(S.59) gate: S^pi vs FD at h = {h:.0e}", eP, tolP[i])
        if h == 1e-7:
            fd7 = (fdS, fdP)
    print(f"  delta:C / max|C| = {np.abs(np.einsum('iikl->kl', C)).max() / np.abs(C).max():.1e} (sheet <= 4e-15: the last sub-increment is post-floored)")
    allok &= ok("delta : C = 0", np.abs(np.einsum("iikl->kl", C)).max() / np.abs(C).max(), 1e-13)
    print("  mutants of the row, C vs FD at h = 1e-7 (pass threshold 3e-7) and S^pi vs FD:")
    for name, sheet in X["mutants"].items():
        Cm, Sm, curm, _ = r4_chain(P, st1, deps, fr, mutant=name)
        assert np.array_equal(curm.sigma, stn.sigma) and curm.pi_i == stn.pi_i          # a mutant of the ROW only
        eC = _c_err(Cm, fd7[0])
        eP = float(np.abs(Sm - fd7[1]).max() / np.abs(Sm).max())
        print(f"    {name}: C {eC:.3f} (sheet {sheet})   S^pi {eP:.3f}")
        allok &= ok(f"{name} is far from FD (>= 5e-3, i.e. 1/err <= 200)", 1.0 / eC, 200.0)
        allok &= ok(f"{name} reproduces the sheet's figure (relative)", abs(eC - sheet) / sheet, 0.1)
    # m = 1: the single-step tangent is unchanged by the rule (M-R9); the chain at m = 1 reduces to (S.32f)/(S.33)
    s1 = step_fractions(P, st1, d1, [1.0])
    with rule_off():
        s1off = step_fractions(P, st1, d1, [1.0])
    e1 = relmax(tangent(P, s1), _s33(s1))
    print(f"  m = 1 ({s1.flags['rpattern']}): chain vs closed form (S.32f)/(S.33) {e1:.1e} (sheet 5e-16); tangent with the rule == rule off, "
          f"bit for bit: {np.array_equal(tangent(P, s1), tangent(P, s1off))}; pi_i {s1.pi_i:.10f} vs {s1off.pi_i:.10f}")
    assert s1.flags["rpattern"] == "-Pf*" and np.array_equal(tangent(P, s1), tangent(P, s1off)) and s1.pi_i != s1off.pi_i
    allok &= ok("m = 1 chain == (S.33) with the rule fired (M-R9)", e1, 1e-13)
    # the reset sub-increment feeding an ordinary (S.54) row: BA06 K1.12 (C), 2x the increment + shear 2e-4, (1/2, 1/2) (case (d2))
    Pb, sb, db, _ = k112c_case()
    dd = 2.0 * db + 2e-4 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    sd = step_fractions(Pb, sb, dd, fr)
    Cd, Sd, _, _ = r4_chain(Pb, sb, dd, fr)
    print(f"  BA06 K1.12 (C) 2x + shear 2e-4, (1/2, 1/2): pattern {sd.flags['rpattern']} (sheet (d2): -Pf*,-P-)")
    for i, h in enumerate((1e-6, 1e-7, 1e-8)):
        fdS, fdP = r4_fd(Pb, sb, dd, fr, h)
        eC, eP = _c_err(Cd, fdS), float(np.abs(Sd - fdP).max() / np.abs(Sd).max())
        print(f"    h = {h:.0e}: C vs FD {eC:.2e} (sheet {(3.5e-8, 8.4e-9, 7.8e-9)[i]:.1e})   S^pi vs FD {eP:.2e} (sheet {(7.4e-9, 1.4e-8, 1.3e-8)[i]:.1e})")
        if h == 1e-7:
            allok &= ok("(d2) BA06 chain with a reset feeding an (S.54) row: C vs FD, h = 1e-7", eC, 3e-7)
            allok &= ok("(d2) S^pi vs FD, h = 1e-7", eP, 3e-7)
            if sd.flags["rpattern"] == "-Pf*,-P-":
                Cm = r4_chain(Pb, sb, dd, fr, mutant="M-R1")[0]
                print(f"    M-R1 on (d2): C {_c_err(Cm, fdS):.3f} (sheet 0.089)")
                allok &= ok("M-R1 on (d2) far from FD", 1.0 / _c_err(Cm, fdS), 200.0)
    print("  r4 (S.59):", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ---- (S.58) on a random census of post-floor events, K1.14b "never shrink" (M-R7)
def _census_set(name, P, rng, n_paths, n_steps, stats):
    F_tol = K.F_TRIAL_TOL_REL * P.p_ref
    dF = F_tol / P.p_min
    cn = c_N(P)
    fac = (1.0 + dF * P.N / P.M) ** (-(1.0 - P.N) / P.N) if P.N > 0.0 else math.exp(-dF / P.M)
    c2 = 0.0 if P.cap == "none" else P.c2
    worst = dict(a=0.0, b=-1e300, c=1e300, Fz=0.0, pif=0.0)
    n = dict(steps=0, post=0, fired=0, notfired=0, band=0, guard=0, refused=0, elastic=0, sub=0, tr=0)
    for _ in range(n_paths):
        p0 = -P.p_min * rng.uniform(1.0, 2.5)
        eta0 = rng.uniform(0.0, 1.1) * P.M
        dvec = rng.normal(size=3)
        sig, _, _ = off_corner_sig(P, p0, eta0, direction=tuple(dvec))
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        pi_g = K.pi_of_eta(P, p0, max(eta0, c2 * P.M))
        try:
            st = initial_state(P, Q @ np.diag(sig) @ Q.T, v_for_psi(P, pi_g, rng.uniform(-0.05, 0.15)), None)   # (S.53): admissible
        except ValueError:
            continue
        assert F_of(P, st.sigma, st.pi_i) <= F_tol
        A = rng.normal(size=(3, 3))
        A = sym(A)
        A -= np.trace(A) / 3.0 * I3
        A /= np.linalg.norm(A)
        mag = 10.0 ** rng.uniform(-5.0, -3.0)
        trv = rng.uniform(-0.3, 1.0)
        for _k in range(n_steps):
            B = rng.normal(size=(3, 3))
            deps = mag * (A + 0.3 * sym(B) / np.linalg.norm(B) + trv / 3.0 * I3)
            with warnings.catch_warnings():
                warnings.simplefilter("error")                 # the (S.57) guard warning would be an error here
                new = step(P, st, deps)
            if new.flags["refused"]:
                n["refused"] += 1
                break
            n["steps"] += 1
            n["sub"] += new.flags["substeps"] > 1
            n["post"] += new.flags["floor_post"]
            n["tr"] += new.flags["floor_tr"]
            n["fired"] += new.flags["floor_pi"]
            n["notfired"] += new.flags["floor_post"] - new.flags["floor_pi"]
            n["band"] += new.flags["floor_pi_band"]
            n["guard"] += new.flags["floor_pi_guard"]
            n["elastic"] += not new.flags["plastic"]
            p = float(np.trace(new.sigma)) / 3.0
            Fc = F_of(P, new.sigma, new.pi_i)
            worst["a"] = max(worst["a"], (P.p_min - abs(p)) / P.p_min)                    # (a) |p| >= p_min (1 - 1e-12)
            worst["b"] = max(worst["b"], Fc / F_tol)                                      # (b) F <= F_tol
            worst["c"] = min(worst["c"], abs(new.pi_i) / (cn * abs(p) * fac))             # (c) >= 1
            assert new.n_f_pi == st.n_f_pi + new.flags["floor_pi"] and new.n_f_post == st.n_f_post + new.flags["floor_post"]
            res = new.cache["res"]                                                        # the LAST sub-increment
            if res.floor_post:
                if res.floor_pi:
                    pif, _e = pi_f_of(P, new.sigma)
                    worst["Fz"] = max(worst["Fz"], abs(Fc) / P.p_min)
                    worst["pif"] = max(worst["pif"], abs(new.pi_i - pif) / abs(pif))
                    assert res.F_post > 0.0 and new.pi_i < res.pi_be and res.gamma_pi is not None
                else:
                    assert new.pi_i == res.pi_be and (res.F_post <= 0.0 or res.floor_pi_band or res.floor_pi_guard)
            else:
                assert not res.floor_pi and res.gamma_pi is None
                if not res.plastic:
                    assert new.flags["substeps"] > 1 or new.pi_i == st.pi_i
            st = new
    stats[name] = (n, worst)
    print(f"  {name}: {n['steps']} accepted steps ({n['sub']} substepped, {n['elastic']} elastic, {n['refused']} paths refused), "
          f"trial floors {n['tr']}, post floors {n['post']}: rule fired {n['fired']}, not fired {n['notfired']} (in the band {n['band']}, "
          f"guard {n['guard']})")
    print(f"      (a) max (p_min - |p|)/p_min = {worst['a']:.1e} (<= 1e-12)   (b) max F/F_tol = {worst['b']:.3e} (<= 1)   "
          f"(c) min |pi_i| / (c_N |p| (1 + dF N/M)^(-(1-N)/N)) = {worst['c']:.6f} (>= 1)   fired: max |F(sigma_f, pi_i)|/p_min = "
          f"{worst['Fz']:.1e}, max |pi_i - pi_f(sigma_f)|/|pi_f| = {worst['pif']:.1e}")
    return n, worst


def check_r4_s58():
    print("\n=== 12e. (S.58) invariants on a random census of post-floor events (admissible starts by (S.53)); never shrink ===")
    allok = True
    rng = np.random.default_rng(20261111)
    stats = {}
    sets = [("K3 HAR smooth cap p_min 0.505", k3_params(), 60, 40),
            ("TIMs HAR no cap p_min 0.505", tims_params(), 60, 40),
            ("BA06 K2 p_min 50", k2_params(p_min=50.0), 60, 40),
            ("BA06 K2 alpha0 5 p_min 50", k2_params(p_min=50.0, alpha0=5.0), 40, 40),
            ("BA06 N = 0 (N_bar 0, rho = rho_bar 0.8) p_min 50", k2_params(rho=0.8, rho_bar=0.8, p_min=50.0, N=0.0, N_bar=0.0), 40, 40)]
    tot_fired = tot_not = 0
    for name, P, npth, nst in sets:
        n, w = _census_set(name, P, rng, npth, nst, stats)
        tot_fired += n["fired"]
        tot_not += n["notfired"]
        allok &= ok(f"{name}: (S.58a) |p| >= p_min (1 - 1e-12)", w["a"], 1e-12)
        allok &= ok(f"{name}: (S.58b) F <= F_tol at every committed state (F/F_tol)", w["b"], 1.0)
        allok &= ok(f"{name}: (S.58c) |pi_i| bound (1/ratio <= 1)", 1.0 / w["c"], 1.0)
        allok &= ok(f"{name}: F(sigma_f, pi_i)/p_min = 0 where the rule fired", w["Fz"], 1e-11)
        allok &= ok(f"{name}: pi_i = pi_f(sigma_f) where the rule fired", w["pif"], 1e-12)
        assert n["guard"] == 0
    print(f"  census total: rule fired at {tot_fired} post-floor events, not fired at {tot_not}")
    assert tot_fired >= 100 and tot_not >= 10, "the census must contain both outcomes"
    # K1.14b, the regression point of "never shrink" (sheet §13.14b round 4; M-R7 would commit -0.740398976)
    P = tims_params()
    sig, th, nh = off_corner_sig(P, -0.6, 1.2 * P.M, direction=(-1.0, -0.35, 1.35))
    pi0 = K.pi_of_eta(P, -0.6, 1.2 * P.M)
    st0 = initial_state(P, np.diag(sig), v_for_psi(P, pi0, -0.10), pi0)
    st = step(P, st0, np.diag(2e-5 / 3.0 * K.ONES + 2e-5 * nh))
    res = st.cache["res"]
    pif, _ = pi_f_of(P, st.sigma)
    print(f"  K1.14b single {st.flags['rpattern']} (FPf, no star): F(sigma_f, pi_BE)/p_min = {res.F_post / P.p_min:+.6f} (-0.004098), "
          f"pi_i = {st.pi_i:.9f} (-0.743350807 unchanged), sigma_f = {np.diag(st.sigma)} (sheet -0.8497174515, -0.5790682544, "
          f"-0.0862142940); the always-reset mutant M-R7 would commit pi_f = {pif:.9f} (-0.740398976)")
    assert st.flags["rpattern"] == "FPf" and st.n_f_pi == 0 and st.pi_i == res.pi_be
    allok &= ok("K1.14b F(sigma_f, pi_BE)/p_min = -0.004098", abs(res.F_post / P.p_min + 0.004098), 1e-6)
    allok &= ok("K1.14b pi_i unchanged = -0.743350807", abs(st.pi_i + 0.743350807), 2e-9)
    allok &= ok("K1.14b sigma_f", float(np.abs(np.diag(st.sigma) - np.array([-0.8497174515, -0.5790682544, -0.0862142940])).max()), 2e-10)
    allok &= ok("K1.14b M-R7 value pi_f(sigma_f) = -0.740398976", abs(pif + 0.740398976), 2e-9)
    for inc, want, pis in ((np.diag(4e-5 / 3.0 * K.ONES + 4e-5 * nh), "FPf,FPf", -0.743038708),
                           (np.diag(1.4e-5 / 3.0 * K.ONES + 2e-5 * nh), "-P-,-Pf", -0.743573074)):
        s2 = step_fractions(P, st0, inc, (0.5, 0.5))
        print(f"  K1.14b chain {s2.flags['rpattern']} ({want}): final pi_i = {s2.pi_i:.9f} ({pis}), floor_pi = {s2.flags['floor_pi']}")
        assert s2.flags["rpattern"] == want and s2.n_f_pi == 0
        allok &= ok(f"K1.14b chain {want} final pi_i", abs(s2.pi_i - pis), 2e-9)
    # the two branches of step 5g no march reaches, called directly on kernel.floor_pi_rule
    for Pv in (k3_params(), k2_params(p_min=50.0), k2_params(rho=0.8, rho_bar=0.8, p_min=50.0, N=0.0, N_bar=0.0)):
        sv = -Pv.p_min * K.ONES                                  # vertex (q = 0) at the floor, pi_BE inside the tip: F = p eta > 0
        aev = K.elastic(Pv, K.invert_elastic(Pv, sv)).ae
        o = K.floor_pi_rule(Pv, sv, aev, -0.3 * c_N(Pv) * Pv.p_min)
        assert o["floor_pi"] and o["eta_f"] == 0.0 and not np.any(o["gamma_pi"]) and o["F_post"] > 0.0
        allok &= ok(f"vertex (N = {Pv.N:.3f}, {Pv.energy}): pi_f = c_N p_f, S^pi row zero", abs(o["pi"] + c_N(Pv) * Pv.p_min) / Pv.p_min, 1e-15)
        allok &= ok("vertex: F(sigma_f, pi_f) = 0", abs(K.flow(Pv, K.invariants(sv), o["pi"]).F) / Pv.p_min, 1e-14)
        o = K.floor_pi_rule(Pv, sv, aev, -2.0 * c_N(Pv) * Pv.p_min)       # F < 0: nothing is done
        assert not o.get("floor_pi", False) and o["pi"] == -2.0 * c_N(Pv) * Pv.p_min and o["F_post"] < 0.0
        if Pv.N > 0.0:                                           # guard: eta_f >= M/N (unreachable from a return; forced here)
            nh_ = np.array([1.0, 0.2, -1.2])
            nh_ /= np.linalg.norm(nh_)
            sg = -Pv.p_min * K.ONES + K.SQ23 * (1.05 * Pv.M / Pv.N / 0.5) * Pv.p_min * nh_      # zeta >= rho > 0.5
            with warnings.catch_warnings(record=True) as wl:
                warnings.simplefilter("always")
                o = K.floor_pi_rule(Pv, sg, aev, -Pv.p_min)
            print(f"  guard (N = {Pv.N:.3f}, {Pv.energy}): eta_f = {o['eta_f']:.3f} >= M/N = {Pv.M / Pv.N:.3f}: pi_i unchanged {o['pi'] == -Pv.p_min}, "
                  f"floor_pi_guard {o.get('floor_pi_guard', False)}, floor_pi {o.get('floor_pi', False)}, warnings {len(wl)}")
            assert o["pi"] == -Pv.p_min and o["floor_pi_guard"] and not o.get("floor_pi", False) and len(wl) == 1
    # p_min = 0: nothing of §9.7 acts
    P0 = k3_params(p_min=0.0)
    assert not K.floor_project(P0, np.array([1e-3, 1e-3, 1e-3])).active
    print("  r4 (S.58):", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


# ---- D >= 0 and the floor energy accounting with the rule (sheet §9.7 "Dissipation, energy")
R4_MARCH = dict(n=600, minD=6.2e-5, sumD=4.384e-2, sumEf=2.920e-3, Wf=2.987e-3, resid=0.012, resid_1200=0.0060,
                off_sumEf=1.52e-2, off_Wf=1.71e-2, off_resid=0.086)


def _march(P, st, deps, n):
    sD = sE = sW = wext = 0.0
    minD, post, fired, cur = 1e300, 0, 0, st
    psi0 = K.energy_psi(P, np.linalg.eigvalsh(st.eps_e))
    for _ in range(n):
        new = step(P, cur, deps)
        assert not new.flags["refused"], new.flags["reason"]
        res = new.cache["res"]
        Ef = floor_energy(P, new)
        dW = new.W_f - cur.W_f
        # per event: D of the return formed at the PRE-floor converged state, E_f in [0, p_min d eps^f_v]
        if new.flags["substeps"] == 1 and res.plastic:
            sc = K.elastic(P, res.eps_c).sig
            assert new.D == res.dlam * float(np.dot(sc, res.q_a))
        assert -1e-18 <= Ef <= dW * (1.0 + 1e-12) + 1e-18 and abs(dW - P.p_min * (new.eps_f_v - cur.eps_f_v)) <= 1e-15
        sD += new.D
        sE += Ef
        sW += dW
        minD = min(minD, new.D)
        post += new.flags["floor_post"]
        fired += new.flags["floor_pi"]
        wext += 0.5 * float(np.sum((cur.sigma + new.sigma) * deps))
        cur = new
    dpsi = K.energy_psi(P, np.linalg.eigvalsh(cur.eps_e)) - psi0
    return dict(sD=sD, sE=sE, sW=sW, minD=minD, post=post, fired=fired, resid=(wext - (dpsi + sD - sE)) / sD, st=cur)


def check_r4_energy():
    print("\n=== 12f. D >= 0 and the floor energy accounting with the rule (sheet §9.7; 600 steps of 2e-4 from e7367g1_s185) ===")
    allok = True
    X = R4_MARCH
    P = k3_params()
    st, deps, _, _ = h19_state(P, h19_records()["e7367g1_s185"])
    m = _march(P, st, deps, 600)
    print(f"  with the rule: min D = {m['minD']:.2e} ({X['minD']:.1e}) kPa, sum D = {m['sD']:.4e} ({X['sumD']:.3e}), sum E_f = {m['sE']:.4e} "
          f"({X['sumEf']:.3e}) <= W_f = {m['sW']:.4e} ({X['Wf']:.3e}) kPa, post floors {m['post']} (600), resets {m['fired']} (600), "
          f"residual of W_ext - (dPsi + sum D - sum E_f) = {100 * m['resid']:.2f} % of sum D (1.2 %)")
    allok &= ok("march: D > 0 at every step (min D vs sheet 6.2e-5, relative)", abs(m["minD"] - X["minD"]) / X["minD"], 0.01)
    assert m["minD"] > 0.0 and m["post"] == 600 and m["fired"] == 600
    allok &= ok("march: sum D (relative)", abs(m["sD"] - X["sumD"]) / X["sumD"], 2e-4)
    allok &= ok("march: sum E_f (relative)", abs(m["sE"] - X["sumEf"]) / X["sumEf"], 2e-4)
    allok &= ok("march: W_f (relative)", abs(m["sW"] - X["Wf"]) / X["Wf"], 2e-4)
    allok &= ok("march: sum E_f <= W_f", max(0.0, m["sE"] - m["sW"]), 0.0)
    allok &= ok("march: work-balance residual / sum D = 1.2 %", abs(abs(m["resid"]) - X["resid"]), 5e-4)
    m2 = _march(P, st, 0.5 * deps, 1200)
    print(f"  1200 steps of 1e-4: residual {100 * m2['resid']:.2f} % (0.60 %: O(d), the truncation of the trapezoid W_ext); min D = {m2['minD']:.2e}")
    allok &= ok("march 1200: residual / sum D = 0.60 %", abs(abs(m2["resid"]) - X["resid_1200"]), 5e-4)
    assert m2["minD"] > 0.0
    with rule_off():
        m3 = _march(P, st, deps, 600)
    print(f"  rule off (round 3b) on the same path: sum E_f = {m3['sE']:.3e} ({X['off_sumEf']:.2e}), W_f = {m3['sW']:.3e} ({X['off_Wf']:.2e}), "
          f"residual {100 * m3['resid']:.1f} % (8.6 %), min D = {m3['minD']:.2e}")
    allok &= ok("rule off: sum E_f, W_f (relative)", max(abs(m3["sE"] - X["off_sumEf"]) / X["off_sumEf"], abs(m3["sW"] - X["off_Wf"]) / X["off_Wf"]), 5e-3)
    allok &= ok("rule off: residual 8.6 %", abs(abs(m3["resid"]) - X["off_resid"]), 1e-3)
    print("  r4 energy:", "ALL OK" if allok else "SOME CHECK FAILED")
    return allok


def check_r4():
    res = [check_r4_k112c(), check_r4_wet(), check_r4_k116(), check_r4_s59(), check_r4_s58(), check_r4_energy()]
    print("\n  r4 group:", "ALL OK" if all(res) else "SOME CHECK FAILED")
    return all(res)


def _kernel_shim():
    """(Kernel, ns_kernel module, note) from kernel_parity/ when g++ and the header are available, else (None, None, why).
    O2 itself never imports the kernel; this is only the (c) cross-check of the D2 self-check."""
    import os
    import sys
    import tempfile
    if sys.platform == "win32":
        return None, None, "Windows: the kernel_parity shim is a Zone B (Linux / Esmeralda, g++) build"
    try:
        kp = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "kernel_parity")
        if kp not in sys.path:
            sys.path.insert(0, kp)
        import ns_kernel as NK
        if NK.compiler() is None:
            return None, None, "no g++ on PATH (Zone B: run on Esmeralda)"
        if not os.path.exists(NK.KERNEL_H):
            return None, None, f"kernel header not found: {NK.KERNEL_H}"
        lib = NK.build(tempfile.mkdtemp(prefix="ns_kernel_d2_"))
        return NK.Kernel(lib), NK, ""
    except Exception as e:                                   # noqa: BLE001 - the shim is optional here
        return None, None, f"{type(e).__name__}: {e}"


if __name__ == "__main__":
    import sys
    groups = dict(jac=check_jacobian, cto=check_cto, newton=check_quadratic, k1=check_k1,
                  k2=lambda: check_finite_and_k2(full_table="--table" in sys.argv), cap=check_cap,
                  chain=check_chain, har=check_har, floor=check_floor, pi0=check_pi0, d2=check_d2,
                  r4=check_r4, r4c=check_r4_k112c, r4wet=check_r4_wet, r4k116=check_r4_k116, r4s59=check_r4_s59,
                  r4s58=check_r4_s58, r4energy=check_r4_energy)
    # no verb: every group once (the r4 sub-verbs r4c ... r4energy are the parts of r4)
    sel = [a for a in sys.argv[1:] if a in groups] or [g for g in groups if not (g.startswith('r4') and g != 'r4')]
    for g in sel:
        groups[g]()
