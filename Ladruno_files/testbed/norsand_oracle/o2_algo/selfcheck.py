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


if __name__ == "__main__":
    import sys
    groups = dict(jac=check_jacobian, cto=check_cto, newton=check_quadratic, k1=check_k1,
                  k2=lambda: check_finite_and_k2(full_table="--table" in sys.argv), cap=check_cap,
                  chain=check_chain, har=check_har, floor=check_floor, pi0=check_pi0)
    sel = [a for a in sys.argv[1:] if a in groups] or list(groups)
    for g in sel:
        groups[g]()
