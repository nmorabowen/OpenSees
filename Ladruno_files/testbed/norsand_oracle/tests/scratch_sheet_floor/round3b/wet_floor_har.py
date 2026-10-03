"""Round 3b (2026-10-03), Adversary item 5: size of F(sigma_f) after a WET-side post floor under HAR, and the
self-limiting of the plastic contraction against the domain edge.

State: AT the floor (p = -p_min = -0.505 kPa, TIMs HAR set), on the surface at eta = 0.5 M (wet side), off the WW
corners, psi_i in {0, +0.04, +0.08, +0.12} (loose of critical; set through v in fork mode). Increment: a pure
deviatoric step along the state's n^ of engineering shear strain Delta_gamma := sqrt(2) ||dev deps|| (the gamma of a
simple-shear increment; deps_s = Delta_gamma/sqrt(3)) in {1e-4, 1e-3, 1e-2}, one backward-Euler step with the floor
at the trial and after convergence (no substepping; a refusal is reported as such, then the increment is retried
as 2, 4, ... equal sub-increments as O2's api.step does).
Reported per case: pattern, p_c (converged, above the floor), F(sigma_f)/p_min after the post floor (the 'slightly
outside F' of sheet §9.7 item 3), deps^p_v of the return (contraction < 0 ... in sheet signs the elastic eps_v moves
UP by |deps^p_v| toward the domain edge), the margin eps*_f(eps_s) = edge - eps_v,f at the committed state, and the
in-Newton minimum of eps* over the accepted iterates (how close the local solve came to leaving dom Psi).
Run:  python -u wet_floor_har.py
"""
import math

import numpy as np

import har_patch as HP
from har_patch import K, ONES, SQ23

from conftest import make_params          # noqa: E402

H = HP.HarParams(**HP.TIMS_HAR)
HP.install(H)
PMIN = 5.0e-3 * H.p_a
kw = dict(p0=-H.p_a, M=1.3309, N=0.4, N_bar=0.2, chi=-3.5, h=280.0, rho=0.71, rho_bar=0.71, zeta="WW",
          csl_mode="fork", e0=0.83, lam_c=0.027, xi=0.45, p_a=H.p_a, cap="none")
P = make_params("O2", **kw)

# instrument the local Newton: record the smallest eps* seen at any EVALUATED iterate (accepted or backtracked)
_orig_eval = K.evaluate
min_est = [float("inf")]
n_domain_fail = [0]


def eval_logged(P_, eps_e, dlam, eps_tr, v, pi_n):
    est = H.edge - float(eps_e.sum())
    min_est[0] = min(min_est[0], est)
    if est <= 0.0:
        n_domain_fail[0] += 1
    return _orig_eval(P_, eps_e, dlam, eps_tr, v, pi_n)


K.evaluate = eval_logged


def v_for_psi(pi, psi):
    return 1.0 + P.e0 - P.lam_c * (-pi / P.p_a) ** P.xi + psi


def F_of(sig, pi):
    return K.flow(P, K.invariants(sig), pi).F


def one_step(eps_n, pi_n, v_n, d):
    v = v_n * math.exp(float(d.sum()))
    etf, at, _, dfv_t = HP.floor_op(H, eps_n + d, PMIN)
    res = K.return_map(P, etf, pi_n, v, v)
    if res.refused:
        return None, res.reason
    ef, ap, _, dfv_p = HP.floor_op(H, res.eps_e, PMIN)
    pat = ("F" if at else "-") + ("P" if res.plastic else "E") + ("f" if ap else "-")
    el_c = K.elastic(P, res.eps_e)
    el_f = K.elastic(P, ef)
    depv_p = float((etf + 0.0).sum()) - float(res.eps_e.sum())      # deps^p_v = eps_v,tr(floored) - eps_v,c
    return dict(pat=pat, p_c=el_c.p, p_f=el_f.p, q_c=el_c.q, q_f=el_f.q, F_f=F_of(el_f.sig, res.pi),
                F_c=F_of(el_c.sig, res.pi), depv_p=depv_p, dlam=res.dlam, eta=res.eta, eps_f=ef, pi=res.pi, v=v,
                est_f=H.edge - float(ef.sum()), eps_s=el_f.eps_s, dfv=dfv_t + dfv_p), ""


def run(eps_n, pi_n, v_n, d, max_halvings=8):
    """one step; on refusal retry as 2^k equal sub-increments (floor at every sub-increment)."""
    for kk in range(max_halvings + 1):
        m = 2 ** kk
        eps_k, pi_k, v_k = eps_n.copy(), pi_n, v_n
        pats, worst_F, max_contr, p_c_min = [], -1e9, 0.0, -1e9
        ok = True
        for _ in range(m):
            out, reason = one_step(eps_k, pi_k, v_k, d / m)
            if out is None:
                ok = False
                last = reason
                break
            pats.append(out["pat"])
            worst_F = max(worst_F, out["F_f"] / PMIN)
            # contraction: deps^p_v < 0 in sheet signs; its size is eps_v,c - eps_v,tr(floored) = -deps^p_v
            max_contr = max(max_contr, -out["depv_p"]) if out["pat"][1] == "P" else max_contr
            p_c_min = max(p_c_min, out["p_c"]) if out["pat"][1] == "P" else p_c_min
            eps_k, pi_k, v_k = out["eps_f"], out["pi"], out["v"]
        if ok:
            return m, pats, worst_F, max_contr, p_c_min, out
    return None, [last], None, None, None, None


print(f"HAR TIMs, p_min = {PMIN}, domain edge eps_v = {H.edge:.6e}; eps*_f(eps_s = 0) = "
      f"{HP.har_floor_evf(H, 0.0, PMIN)[4]:.3e}")
sig0, th0, nh0 = HP.off_corner_sig(P, -PMIN, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
pi0 = K.pi_of_eta(P, -PMIN, 0.5 * P.M)
eps0 = K.invert_elastic(P, sig0)
el0 = K.elastic(P, eps0)
print(f"state at the floor: p = {el0.p:.4f}, q = {el0.q:.4f}, eta = {0.5 * P.M:.4f}, theta = {th0:.3f}, pi_i = {pi0:.5f}, "
      f"eps_v = {el0.eps_v:.6e} (margin eps* = {H.edge - el0.eps_v:.3e}), eps_s = {el0.eps_s:.3e}")
print(f"{'psi_i':>6} {'dgamma':>7} {'m':>3} {'pattern(s)':<22} {'max F(sig_f)/p_min':>18} {'max|deps^p_v| (contr.)':>23} "
      f"{'p_c max':>8} {'eps*_f end':>10} {'min eps* in Newton':>18} {'dom fails':>9}")
hdr_done = True
summary = {}
for psi in (0.0, 0.04, 0.08, 0.12):
    v0 = v_for_psi(pi0, psi)
    for dg in (1e-4, 1e-3, 1e-2):
        min_est[0] = float("inf")
        n_domain_fail[0] = 0
        d = (dg / math.sqrt(2.0)) * nh0                      # ||dev d|| = dgamma/sqrt2
        m, pats, wF, mc, pcm, out = run(eps0, pi0, v0, d)
        if m is None:
            print(f"{psi:6.2f} {dg:7.0e} {'-':>3} {pats[0]:<22} REFUSED")
            continue
        pat_s = ",".join(pats) if m <= 4 else f"{pats[0]} x{m} ({''.join(sorted(set(pats)))})"
        print(f"{psi:6.2f} {dg:7.0e} {m:3d} {pat_s:<22} {wF:18.4f} {mc:23.3e} {pcm:8.4f} {out['est_f']:10.3e} "
              f"{min_est[0]:18.3e} {n_domain_fail[0]:9d}")
        summary[(psi, dg)] = (wF, mc, m)
# the (p, q) picture for one case: the contributions to F after the post floor
psi, dg = 0.12, 1e-2
v0 = v_for_psi(pi0, psi)
m, pats, wF, mc, pcm, out = run(eps0, pi0, v0, (dg / math.sqrt(2.0)) * nh0)
if out is not None:
    Y = K.yield_p(P, out["p_f"], out["pi"])
    zeta_f = K.zeta_y(K.invariants(K.elastic(P, out["eps_f"]).sig).theta, P.rho, P.zeta)[0]
    est1 = abs(Y.F_p) * (PMIN - abs(out["p_c"]))
    est2 = zeta_f * (out["q_f"] - out["q_c"])
    print(f"\ncase psi_i = {psi}, dgamma = {dg}, m = {m}: last sub-increment p_c = {out['p_c']:.4f} -> p_f = {out['p_f']:.4f}, "
          f"q_c = {out['q_c']:.4f} -> q_f = {out['q_f']:.4f}, eta_c = {out['eta']:.4f}, F_c/p_min = {out['F_c'] / PMIN:+.2e}, "
          f"F_f/p_min = {out['F_f'] / PMIN:+.4f}; first-order estimate of sheet 9.7 item 2: |F_p| (p_min - |p_c|) = {est1:.4f} kPa "
          f"+ zeta (q_f - q_c) = {est2:.4f} kPa = {(est1 + est2) / PMIN:.3f} p_min (F_p = {Y.F_p:+.4f}, zeta = {zeta_f:.4f}); "
          f"the exact F differs from this first-order size at second order in dp (dp = {PMIN - abs(out['p_c']):.3f} kPa is "
          f"half of p_min here)")
K.evaluate = _orig_eval
HP.uninstall()
print("\nmax over the table: F(sigma_f)/p_min = %.3f, converged contraction |deps^p_v| = %.3e (vs the margin at the floor "
      "eps*_f(0) = %.3e, eps*_f(eps_s of the state) = %.3e): the converged return never leaves dom Psi; the full Newton "
      "steps do overshoot it (dom fails > 0) and the line search backtracks them." % (
          max(v[0] for v in summary.values()), max(v[1] for v in summary.values()), HP.har_floor_evf(H, 0.0, PMIN)[4],
          H.edge - el0.eps_v))
print("done")
