"""O2 self-checks (sheet 144a). Run:  python -m o2_algo.selfcheck   from norsand_oracle/.
Prints numbers; it is NOT the gate test suite (P0d writes that)."""
from __future__ import annotations

import math
import warnings

import numpy as np

from . import kernel as K
from .acoustic import acoustic_min_det, acoustic_principal, acoustic_tensor
from .api import (State, initial_state, k2_path, run_path, step, step_fractions, tangent, tangent_finite,
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


def sym(a):
    return 0.5 * (a + a.T)


# ------------------------------------------------------------------ 1. Jacobian vs FD
def jacobian_fd(P, st, deps, h=1e-6):
    """returns (rel err, local iters) at the converged plastic state of the step st -> st+deps."""
    tr = float(np.trace(deps))
    v = st.v + st.v0 * tr
    eps_tr = st.eps_e + sym(deps)
    w, V = np.linalg.eigh(eps_tr)
    res = K.return_map(P, w, st.pi_i, v, st.v0)
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
    S = np.diag([1, 1, 1, 1 / abs(P.p0)])
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
    res = K.return_map(P, w, st.pi_i, st.v + st.v0 * tr, st.v0)
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
    print(f"  a~^ep (S.32 with v0->v) vs FD, finite mode, plastic: rel err {np.linalg.norm(at - at_fd) / np.linalg.norm(at):.2e}")
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
    e, el = 0.0, 0.0
    for J in range(6):
        E = K.CHAIN_E[J]
        sp = step_fractions(P, st, deps + h * E, fractions)
        sm = step_fractions(P, st, deps - h * E, fractions)
        assert sp.flags["pattern"] == pat == sm.flags["pattern"], (pat, sp.flags["pattern"], sm.flags["pattern"])
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


if __name__ == "__main__":
    import sys
    groups = dict(jac=check_jacobian, cto=check_cto, newton=check_quadratic, k1=check_k1,
                  k2=lambda: check_finite_and_k2(full_table="--table" in sys.argv), cap=check_cap,
                  chain=check_chain)
    sel = [a for a in sys.argv[1:] if a in groups] or list(groups)
    for g in sel:
        groups[g]()
