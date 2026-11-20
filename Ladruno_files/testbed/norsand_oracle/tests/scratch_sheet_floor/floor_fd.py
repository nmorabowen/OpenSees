"""Sheet round 2026-10-03: FD check of the p-floor operator inside the backward-Euler step (sheet S.48-S.51,
S.32f, S.46f) using the shipped O2 kernel (BA06 energy, K2 set, alpha0 = 0 and alpha0 = 5 so that D12 != 0
exercises the eps' term). The floor is implemented HERE (O2 has none yet): Pi_f in strain space at fixed
deviatoric elastic strain, applied to the trial and after convergence, p_min scaled up to 50 kPa so that the
floor engages at K2 stress levels (the algebra is scale-free). All states are OFF the Willam-Warnke corners
(theta ~ 0.4-0.6 rad): at the exact TXC corner a symmetry-breaking central FD is O(h) (sheet 4.3), which the
first draft of this script ran into.

  (A) elastic step with the TRIAL floored:   C_f = a^e(eps_f) DPi_f^tr                         vs central FD
  (B) plastic step with the trial floored:   C_f = a^e(eps_c) [b Phi_tr - u Pi_v v 1^T]           vs FD
      (the v column stays attached to the RAW trial strain: the floor does not touch v; the mutant that
       ties it to the floored trial, a~^ep Phi_tr, is reported)
  (C) plastic step with the POST floor:      C_f = a^e(eps_f) Phi_post [b - u Pi_v v 1^T]          vs FD
  (D) m = 2 chain, sub-increment 1 trial-floored and plastic: (S.46f) vs FD of the whole increment.
Every FD point must stay in the same branch (floor active on both sides); asserted.
The O2 oracle is imported read-only from the WP-144 worktree (ORACLE below).
Run:  python -u floor_fd.py
"""
import math
import os
import sys
import warnings

import numpy as np

ORACLE = r"C:/Users/nmora/Documents/Github/OpenSees/.claude/worktrees/ladrunonorsand-implementation-review-7bbd75/Ladruno_files/testbed/norsand_oracle"
sys.path.insert(0, ORACLE)
sys.path.insert(0, os.path.join(ORACLE, "tests"))
warnings.simplefilter("ignore")
from conftest import K2_BASE, make_params          # noqa: E402
from o2_algo import kernel as K                   # noqa: E402

I3 = np.eye(3)
ONES = np.ones(3)
SQ23 = math.sqrt(2.0 / 3.0)
np.set_printoptions(precision=4, linewidth=150)
ok_all = True


def rep(name, err, tol):
    global ok_all
    ok = err <= tol
    ok_all &= ok
    print(f"  [{'OK' if ok else 'FAIL'}] {name}: {err:.3e} (tol {tol:.0e})")


def floor_op(P, eps_e, pmin):
    """Pi_f (sheet S.48-S.49, BA06): returns (eps_f, active, Phi_block, deps_f_v).
    eps_v,f = eps_v0 - kappa ln[p_min / (|p0| (1 + 3 alpha0 eps_s^2/(2 kappa)))], eps' = 3 alpha0 eps_s/(1 + ...)."""
    ev = float(eps_e.sum())
    e = eps_e - ev / 3.0
    ne = float(np.linalg.norm(e))
    es = SQ23 * ne
    nh = e / ne if es > K.EPS_S_TOL else np.zeros(3)
    el = K.elastic(P, eps_e)
    if el.p <= -pmin * (1.0 - 1e-12):
        return eps_e.copy(), False, I3.copy(), 0.0
    fac = 1.0 + 1.5 * P.alpha0 * es * es / P.kappa_hat
    ev_f = P.eps_v0 - P.kappa_hat * math.log(pmin / (abs(P.p0) * fac))
    epsp = 3.0 * P.alpha0 * es / fac
    eps_f = eps_e - (ev - ev_f) / 3.0
    Phi = I3 - 1.0 / 3.0 + (1.0 / 3.0) * epsp * SQ23 * np.outer(ONES, nh)
    return eps_f, True, Phi, ev - ev_f


def step_floor(P, eps_e_n, pi_n, v_n, deps3, pmin):
    """One sub-increment with the floor at the trial and after convergence (principal, co-axial case)."""
    v = v_n * math.exp(float(deps3.sum()))
    eps_tr = eps_e_n + deps3
    eps_trf, act_tr, Phi_tr, dfv_tr = floor_op(P, eps_tr, pmin)
    res = K.return_map(P, eps_trf, pi_n, v, v)
    assert not res.refused, res.reason
    eps_c = res.eps_e
    eps_f, act_post, Phi_post, dfv_post = floor_op(P, eps_c, pmin)
    el = K.elastic(P, eps_f)
    if res.plastic:
        ch = res.chain
        dxde = ch.b[:3, :3] @ Phi_tr - np.outer(ch.u[:3] * ch.Pi_v * v, ONES)   # (S.32f): v column on the RAW trial
        dxde_mut = np.linalg.solve(res.ae, res.atilde) @ Phi_tr                 # mutant: v column tied to the floored trial
    else:
        dxde = Phi_tr
        dxde_mut = Phi_tr
    C = el.ae @ Phi_post @ dxde
    C_mut_v = el.ae @ Phi_post @ dxde_mut
    return dict(sig=el.sig, eps_f=eps_f, pi=res.pi, v=v, res=res, C=C, C_mut_v=C_mut_v, act_tr=act_tr,
                act_post=act_post, dfv=dfv_tr + dfv_post, p_tr=K.elastic(P, eps_tr).p, p_c=K.elastic(P, eps_c).p,
                p_f=el.p)


def fd_tangent(P, eps_e_n, pi_n, v_n, deps3, pmin, h):
    C = np.zeros((3, 3))
    acts = set()
    for b in range(3):
        dp = deps3.copy()
        dp[b] += h
        dm = deps3.copy()
        dm[b] -= h
        sp_ = step_floor(P, eps_e_n, pi_n, v_n, dp, pmin)
        sm_ = step_floor(P, eps_e_n, pi_n, v_n, dm, pmin)
        acts.add((sp_["act_tr"], sp_["act_post"], sp_["res"].plastic))
        acts.add((sm_["act_tr"], sm_["act_post"], sm_["res"].plastic))
        C[:, b] = (sp_["sig"] - sm_["sig"]) / (2 * h)
    return C, acts


def relerr(A, B):
    return float(np.abs(A - B).max() / np.abs(B).max())


def off_corner_sig(P, p_s, eta_s, direction=(-1.0, -0.35, 1.35)):
    """principal stress at (p, eta) with a deviatoric direction away from both WW corners."""
    xi = np.array(direction)
    xi -= xi.mean()
    nh = xi / np.linalg.norm(xi)
    inv = K.invariants(p_s + nh)           # theta of the direction
    z, _, _ = K.zeta_y(inv.theta, P.rho, P.zeta)
    q_s = eta_s * abs(p_s) / z             # F = 0: zeta q = -p eta
    sig = p_s + SQ23 * q_s * nh
    return sig, inv.theta


PMIN = 50.0
for alpha0 in (0.0, 5.0):
    print(f"\n================ BA06 K2 set, alpha0 = {alpha0}, p_min = {PMIN} kPa ================")
    kw = dict(K2_BASE, rho=0.7, rho_bar=0.8, alpha0=alpha0)
    P = make_params("O2", **kw)
    # (A) elastic + trial floor: off-corner state at p ~ -60 inside a large surface (pi_i = -300), expansion step
    eps_n = K.invert_elastic(P, np.array([-57.0, -61.0, -64.0]))
    pi_n, v_n = -300.0, 1.70
    deps = np.array([2.0e-3, 1.5e-3, 2.5e-3])
    s = step_floor(P, eps_n, pi_n, v_n, deps, PMIN)
    th = K.invariants(s["sig"]).theta
    print(f"(A) theta = {th:.3f}, p_tr = {s['p_tr']:.3f}, p_conv = {s['p_c']:.3f}, p_committed = {s['p_f']:.3f}, plastic = {s['res'].plastic}, "
          f"floor tr/post = {s['act_tr']}/{s['act_post']}, deps^f_v = {s['dfv']:.3e}")
    assert s['act_tr'] and not s['res'].plastic
    for h in (1e-6, 1e-7):
        Cfd, acts = fd_tangent(P, eps_n, pi_n, v_n, deps, PMIN, h)
        assert len(acts) == 1, acts
        rep(f"(A) elastic+trial floor: C_f vs FD, h = {h:.0e}", relerr(s["C"], Cfd), 3e-7 if h == 1e-6 else 3e-8)
    rep("(A) delta : C_f = 0", float(np.abs(s["C"].sum(axis=0)).max() / np.abs(s["C"]).max()), 1e-14)
    rep("(A) committed p = -p_min", abs(s["p_f"] + PMIN) / PMIN, 1e-13)
    print(f"      mutant 'unprojected a^e' vs FD: {relerr(K.elastic(P, s['eps_f']).ae, Cfd):.3e}")
    if alpha0 != 0.0:
        C_mut2 = K.elastic(P, s["eps_f"]).ae @ (I3 - 1.0 / 3.0)
        print(f"      mutant 'eps'' term dropped' vs FD: {relerr(C_mut2, Cfd):.3e}; its delta:C/max|C| = "
              f"{np.abs(C_mut2.sum(axis=0)).max()/np.abs(C_mut2).max():.3e}")

    # (B) plastic + trial floor: on the surface at p = -55, eta = 1.4 > M (the return compresses p), expansion + shear
    sig_B, thB = off_corner_sig(P, -55.0, 1.4)
    pi_B = K.pi_of_eta(P, -55.0, 1.4)
    eps_B = K.invert_elastic(P, sig_B)
    vB = 1.70
    xiB = np.array([-1.0, -0.35, 1.35])
    nhB = (xiB - xiB.mean()) / np.linalg.norm(xiB - xiB.mean())
    deps_B = None
    for sdev in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3):
        cand = 1.0e-3 * ONES + sdev * nhB            # expansion (trial p above the floor) + shear along the state's n^
        s = step_floor(P, eps_B, pi_B, vB, cand, PMIN)
        if s['act_tr'] and s['res'].plastic and not s['act_post']:
            deps_B = cand
            break
    assert deps_B is not None, "case (B) not in the intended branch"
    print(f"(B) theta = {thB:.3f} -> {K.invariants(s['sig']).theta:.3f}, p_tr = {s['p_tr']:.3f}, p_conv = {s['p_c']:.3f}, "
          f"p_committed = {s['p_f']:.3f}, plastic = {s['res'].plastic}, eta = {s['res'].eta:.3f}, "
          f"floor tr/post = {s['act_tr']}/{s['act_post']}, dlam = {s['res'].dlam:.3e}")
    for h in (1e-6, 1e-7):
        Cfd, acts = fd_tangent(P, eps_B, pi_B, vB, deps_B, PMIN, h)
        assert len(acts) == 1, acts
        rep(f"(B) plastic+trial floor: C_f vs FD, h = {h:.0e}", relerr(s["C"], Cfd), 3e-6 if h == 1e-6 else 3e-7)
    print(f"      mutant 'v column tied to the floored trial (a~^ep Phi_tr)' vs FD: {relerr(s['C_mut_v'], Cfd):.3e}")
    print(f"      mutant 'DPi_tr dropped (plain a~^ep)' vs FD: {relerr(s['res'].atilde, Cfd):.3e}")

    # (C) post floor: wet side (eta < M, F_p < 0): the return makes p LESS compressive; start on the surface at
    #     p = -51, eta = 0.5 M off-corner; deviatoric increment with slight compression (trial below the floor)
    sig_C, thC = off_corner_sig(P, -51.0, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
    pi_C = K.pi_of_eta(P, -51.0, 0.5 * P.M)
    eps_C = K.invert_elastic(P, sig_C)
    found = None
    for amp in (2e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1.2e-2):
        deps_C = amp * np.array([0.55, 0.40, -0.95]) - 1e-5
        s = step_floor(P, eps_C, pi_C, 1.70, deps_C, PMIN)
        if (not s['act_tr']) and s['res'].plastic and s['act_post']:
            found = deps_C
            break
    if found is None:
        print("(C) could not build a post-floor case (reported, not fatal); last try:",
              s['act_tr'], s['res'].plastic, s['act_post'], s['p_tr'], s['p_c'])
    else:
        print(f"(C) theta = {thC:.3f} -> {K.invariants(s['sig']).theta:.3f}, p_tr = {s['p_tr']:.3f}, p_conv = {s['p_c']:.3f}, "
              f"p_committed = {s['p_f']:.3f}, plastic = {s['res'].plastic}, eta = {s['res'].eta:.3f}, "
              f"floor tr/post = {s['act_tr']}/{s['act_post']}, dlam = {s['res'].dlam:.3e}")
        for h in (1e-6, 1e-7):
            Cfd, acts = fd_tangent(P, eps_C, pi_C, 1.70, found, PMIN, h)
            assert len(acts) == 1, acts
            rep(f"(C) plastic+post floor: C_f vs FD, h = {h:.0e}", relerr(s["C"], Cfd), 3e-6 if h == 1e-6 else 3e-7)
        rep("(C) delta : C_f = 0", float(np.abs(s["C"].sum(axis=0)).max() / np.abs(s["C"]).max()), 1e-13)
        print(f"      mutant 'post floor absent from the tangent (plain a~^ep)' vs FD: {relerr(s['res'].atilde, Cfd):.3e}")

    # (D) m = 2 chain with the floor in sub-increment 1 (trial-floored, plastic).
    def chain_floor(P, eps_n, pi_n, v_n, deps3, fr, pmin, with_ops=True):
        """(S.46f) eigenvalue blocks (co-axial case). with_ops=False: the mutant that leaves the DPi operators
        out of the chain while the states are still floored."""
        S_eps = np.zeros((3, 3))
        S_pi = np.zeros(3)
        cum = 0.0
        eps_k, pi_k, v_k = eps_n.copy(), pi_n, v_n
        pattern = []
        for a in fr:
            T = S_eps + a * I3
            v_new = v_k * math.exp(a * float(deps3.sum()))
            eps_trf, act_tr, Phi_tr, _ = floor_op(P, eps_k + a * deps3, pmin)
            res = K.return_map(P, eps_trf, pi_k, v_new, v_new)
            assert not res.refused
            eps_f, act_post, Phi_post, _ = floor_op(P, res.eps_e, pmin)
            if not with_ops:
                Phi_tr, Phi_post = I3, I3
            cum += a
            S_v = v_new * cum * ONES                       # tr E_J = 1 on the three normal columns (raw trial)
            Tf = Phi_tr @ T
            if res.plastic:
                ch = res.chain
                S_new = ch.b[:3, :3] @ Tf - np.outer(ch.u[:3] / ch.c, S_pi) - np.outer(ch.u[:3] * ch.Pi_v, S_v)
                S_pi = ch.w @ Tf + ((1 - ch.kappa) / ch.c) * S_pi + (1 - ch.kappa) * ch.Pi_v * S_v
            else:
                S_new = Tf
            S_eps = Phi_post @ S_new
            eps_k, pi_k, v_k = eps_f, res.pi, v_new
            pattern.append(("F" if act_tr else "-") + ("P" if res.plastic else "E") + ("f" if act_post else "-"))
        el = K.elastic(P, eps_k)
        return el.ae @ S_eps, el.sig, pattern

    fr = [0.5, 0.5]
    deps_D = deps_B.copy()
    C_ch, sig_D, pat = chain_floor(P, eps_B, pi_B, 1.70, deps_D, fr, PMIN)
    print(f"(D) m = 2 chain, branch pattern {pat}  (F = trial floored, P/E plastic/elastic, f = post floored)")
    for h in (1e-6, 1e-7):
        Cfd = np.zeros((3, 3))
        pats = set()
        for b in range(3):
            dp = deps_D.copy()
            dp[b] += h
            dm = deps_D.copy()
            dm[b] -= h
            _, sp_, pp = chain_floor(P, eps_B, pi_B, 1.70, dp, fr, PMIN)
            _, sm_, pm = chain_floor(P, eps_B, pi_B, 1.70, dm, fr, PMIN)
            pats.add(tuple(pp))
            pats.add(tuple(pm))
            Cfd[:, b] = (sp_ - sm_) / (2 * h)
        assert len(pats) == 1, pats
        rep(f"(D) chained C with the floor operators vs FD of the whole increment, h = {h:.0e}", relerr(C_ch, Cfd),
            3e-6 if h == 1e-6 else 3e-7)
    C_mut, _, _ = chain_floor(P, eps_B, pi_B, 1.70, deps_D, fr, PMIN, with_ops=False)
    print(f"      mutant 'DPi operators omitted from the chain' vs FD: {relerr(C_mut, Cfd):.3e}")

print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
