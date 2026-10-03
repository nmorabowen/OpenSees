"""Round 3b (2026-10-03), Adversary item 1: the HAR dry-side 'FPf' floor case and the chained patterns, FD-checked.

The O2 return map runs on the HAR law through har_patch (the plastic part is energy-blind, sheet §2.1/§2.4); the
floor Pi_f (S.48) with the HAR closed form (S.50) is wrapped around `return_map` exactly as floor_fd.py did for BA06.
TIMs elastic set (n 1/2, g 807.80, k 1889.48, p_a 101), p_min = 5e-3 p_a = 0.505 kPa, fork CSL at the TIMs values
(e0 0.83, lambda_c 0.027, xi 0.45), M 1.3309, rho = rho_bar 0.71 (WW); the plastic constants N 0.4, N_bar 0.2,
chi -3.5, h 280 are the K2 values (TIMs' are a P3 refit). All states OFF the WW corners.

  (B-HAR) dry side, trial floored, plastic, POST floor active ('FPf'): surface state p = -0.6 kPa, eta = 1.2 M,
          expansion + shear. The return from the floored trial (p = -p_min) ends ABOVE the floor because the D12
          coupling of the plastic shear strain (+) beats the volumetric compression of the dilative return (-).
          C_f = a^e(eps_f) Phi_post [b Phi_tr - u Pi_v v 1^T]  (S.32f) vs central FD, h = 1e-6 / 1e-7.
  (D-HAR) m = 2 chains: pattern FPf,FPf (the same increment halved) and -P-,-Pf (a smaller increment whose trial
          stays below the floor and whose second return crosses it): (S.54)/(S.46f) vs FD of the whole increment.
  Also: F(sigma_f)/p_min after the dry-side post floor (sign and size), and the two competing dp terms.
Every FD point must stay in the same branch; asserted.
Run:  python -u floor_fd_har.py
"""
import math

import numpy as np

import har_patch as HP
from har_patch import K, I3, ONES, SQ23

from conftest import make_params          # noqa: E402  (ORACLE/tests on sys.path via har_patch)

np.set_printoptions(precision=5, linewidth=150)
ok_all = True


def rep(name, err, tol):
    global ok_all
    ok = err <= tol
    ok_all &= ok
    print(f"  [{'OK' if ok else 'FAIL'}] {name}: {err:.3e} (tol {tol:.0e})")


H = HP.HarParams(**HP.TIMS_HAR)
HP.install(H)
PMIN = 5.0e-3 * H.p_a
kw = dict(p0=-H.p_a, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0,          # BA06 fields unused (p0 := -p_a)
          M=1.3309, N=0.4, N_bar=0.2, chi=-3.5, h=280.0, rho=0.71, rho_bar=0.71, zeta="WW",
          csl_mode="fork", e0=0.83, lam_c=0.027, xi=0.45, p_a=H.p_a, cap="none")
P = make_params("O2", **kw)
print(f"HAR TIMs set, p_min = {PMIN} kPa, domain edge eps_v = {H.edge:.6e}, p0 := -p_a = {P.p0}")


def step_floor(eps_e_n, pi_n, v_n, deps3):
    v = v_n * math.exp(float(deps3.sum()))
    eps_tr = eps_e_n + deps3
    eps_trf, act_tr, Phi_tr, dfv_tr = HP.floor_op(H, eps_tr, PMIN)
    res = K.return_map(P, eps_trf, pi_n, v, v)
    assert not res.refused, res.reason
    eps_c = res.eps_e
    eps_f, act_post, Phi_post, dfv_post = HP.floor_op(H, eps_c, PMIN)
    el = K.elastic(P, eps_f)
    if res.plastic:
        ch = res.chain
        dxde = ch.b[:3, :3] @ Phi_tr - np.outer(ch.u[:3] * ch.Pi_v * v, ONES)
    else:
        dxde = Phi_tr
    C = el.ae @ Phi_post @ dxde
    p_tr = None
    try:
        p_tr = K.elastic(P, eps_tr).p
    except K.EvalError:
        p_tr = float("nan")                      # out of domain
    return dict(sig=el.sig, eps_f=eps_f, pi=res.pi, v=v, res=res, C=C, act_tr=act_tr, act_post=act_post,
                dfv=dfv_tr + dfv_post, p_tr=p_tr, p_c=K.elastic(P, eps_c).p, p_f=el.p, q_c=K.elastic(P, eps_c).q,
                q_f=el.q, eps_c=eps_c, eps_trf=eps_trf, el=el, Phi_post=Phi_post)


def fd_tangent(eps_e_n, pi_n, v_n, deps3, h):
    C = np.zeros((3, 3))
    acts = set()
    for b in range(3):
        dp = deps3.copy()
        dp[b] += h
        dm = deps3.copy()
        dm[b] -= h
        sp_ = step_floor(eps_e_n, pi_n, v_n, dp)
        sm_ = step_floor(eps_e_n, pi_n, v_n, dm)
        acts.add((sp_["act_tr"], sp_["act_post"], sp_["res"].plastic))
        acts.add((sm_["act_tr"], sm_["act_post"], sm_["res"].plastic))
        C[:, b] = (sp_["sig"] - sm_["sig"]) / (2 * h)
    return C, acts


def relerr(A, B):
    return float(np.abs(A - B).max() / np.abs(B).max())


def pattern(s):
    return ("F" if s["act_tr"] else "-") + ("P" if s["res"].plastic else "E") + ("f" if s["act_post"] else "-")


def F_of(sig, pi):
    inv = K.invariants(sig)
    return K.flow(P, inv, pi).F


def v_for_psi(pi, psi_target):
    """v such that psi_i = e - e_c(pi) = psi_target in fork mode."""
    e_c = P.e0 - P.lam_c * (-pi / P.p_a) ** P.xi
    return 1.0 + e_c + psi_target


# ------------------------------------------------------------------------------------------------ (B-HAR)
print("\n================ (B-HAR) dry side, trial floored, plastic, post floor active ('FPf') ================")
p_s, eta_s = -0.6, 1.2 * P.M
sig_B, thB, nhB = HP.off_corner_sig(P, p_s, eta_s)
pi_B = K.pi_of_eta(P, p_s, eta_s)
eps_B = K.invert_elastic(P, sig_B)
elB = K.elastic(P, eps_B)
vB = v_for_psi(pi_B, -0.10)                      # dense of critical at the image pressure
print(f"  state: p = {elB.p:.4f}, q = {elB.q:.4f}, eta = {eta_s:.4f} (= 1.2 M), theta = {thB:.3f}, pi_i = {pi_B:.5f}, "
      f"v = {vB:.5f} (psi_i = -0.10), F = {F_of(sig_B, pi_B):.2e}, eps_v = {elB.eps_v:.6e}, eps_s = {elB.eps_s:.6e}")
print(f"  D11 = {elB.D11:.3f}, D12 = {elB.D12:.3f}, D22 = {elB.D22:.3f} kPa")
found = None
for amp_v in (2e-5, 4e-5, 6e-5, 8e-5):
    for amp_s in (2e-5, 5e-5, 1e-4, 2e-4, 4e-4):
        cand = (amp_v / 3.0) * ONES + amp_s * nhB       # expansion (trial p above the floor) + shear along n^
        s = step_floor(eps_B, pi_B, vB, cand)
        if s["act_tr"] and s["res"].plastic and s["act_post"]:
            found = cand
            break
    if found is not None:
        break
assert found is not None, "FPf not found"
deps_B = found
s = step_floor(eps_B, pi_B, vB, deps_B)
res = s["res"]
# decompose p_c - (-p_min) into the D11 and D12 contributions at the converged state (secant: use converged D's)
eps_trf = s["eps_trf"]
el_trf = K.elastic(P, eps_trf)
dEv = float((s["eps_c"] - eps_trf).sum())
e_c = s["eps_c"] - float(s["eps_c"].sum()) / 3
e_t = eps_trf - float(eps_trf.sum()) / 3
dEs = SQ23 * (np.linalg.norm(e_c) - np.linalg.norm(e_t))
el_c = K.elastic(P, s["eps_c"])
D11m, D12m = 0.5 * (el_trf.D11 + el_c.D11), 0.5 * (el_trf.D12 + el_c.D12)
print(f"  increment: tr deps = {deps_B.sum():.3e} (expansion), shear amp {np.linalg.norm(deps_B - deps_B.mean()):.3e} "
      f"along n^; pattern {pattern(s)}")
print(f"  p_tr = {s['p_tr']:.4f} -> floored trial p = {el_trf.p:.4f} (deps^f_v,tr = {float((eps_B + deps_B - eps_trf).sum()):.3e}) "
      f"-> return p_c = {s['p_c']:.4f} (ABOVE the floor) -> post floor p = {s['p_f']:.4f}; "
      f"dlam = {res.dlam:.3e}, eta_c = {res.eta:.4f}, pi_i = {res.pi:.5f}")
print(f"  return decomposition (mid-point D's): D11 deps^e_v = {D11m * dEv:+.4f} kPa (deps^p_v = {-dEv:+.3e}, dilative), "
      f"D12 deps^e_s = {D12m * dEs:+.4f} kPa (deps^p_s = {-dEs:+.3e}); sum {D11m * dEv + D12m * dEs:+.4f} vs "
      f"p_c - p_trf = {s['p_c'] - el_trf.p:+.4f}")
Ff = F_of(s["sig"], res.pi)
Fc = F_of(el_c.sig, res.pi)
print(f"  yield after the post floor: F(sigma_c)/p_min = {Fc / PMIN:+.3e} (converged, ~0), F(sigma_f)/p_min = {Ff / PMIN:+.4f}; "
      f"q_c = {s['q_c']:.4f} -> q_f = {s['q_f']:.4f} (HAR: q moves with the floored p), dp = {s['p_f'] - s['p_c']:+.4f}")
errs = {}
for h in (1e-6, 1e-7, 1e-8):
    Cfd, acts = fd_tangent(eps_B, pi_B, vB, deps_B, h)
    assert len(acts) == 1, acts
    errs[h] = relerr(s["C"], Cfd)
    rep(f"(B-HAR) FPf: C_f (S.32f) vs FD, h = {h:.0e}", errs[h], {1e-6: 3e-5, 1e-7: 3e-7, 1e-8: 3e-8}[h])
print(f"      O(h^2): err(1e-6)/err(1e-7) = {errs[1e-6]/errs[1e-7]:.1f}, err(1e-7)/err(1e-8) = {errs[1e-7]/errs[1e-8]:.1f} "
      f"(the increment is 2e-5: h = 1e-6 is 5 % of it)")
rep("(B-HAR) delta : C_f = 0 (last operator an active Pi_f)", float(np.abs(s["C"].sum(axis=0)).max() / np.abs(s["C"]).max()), 1e-13)
rep("(B-HAR) committed p = -p_min", abs(s["p_f"] + PMIN) / PMIN, 1e-12)
ch = res.chain
Phi_tr_B = HP.floor_op(H, eps_B + deps_B, PMIN)[2]
dxde_B = ch.b[:3, :3] @ Phi_tr_B - np.outer(ch.u[:3] * ch.Pi_v * s["v"], ONES)
print(f"      mutant 'Phi_post dropped' vs FD: {relerr(s['el'].ae @ dxde_B, Cfd):.3e}")
print(f"      mutant 'Phi_tr dropped' vs FD: {relerr(s['el'].ae @ s['Phi_post'] @ (ch.b[:3, :3] - np.outer(ch.u[:3] * ch.Pi_v * s['v'], ONES)), Cfd):.3e}")
print(f"      mutant 'both operators dropped (plain a~^ep)' vs FD: {relerr(res.atilde, Cfd):.3e}")
print(f"      mutant 'v-column tied to the floored trial (a~^ep Phi_tr, then Phi_post)' vs FD: "
      f"{relerr(s['el'].ae @ s['Phi_post'] @ np.linalg.solve(res.ae, res.atilde) @ Phi_tr_B, Cfd):.3e}")
# the dry-side rule: a 'pure volumetric' floored trial under BA06 cannot do this; show the sign of the two terms
print(f"  -> dry side under HAR: dp_D12 = {D12m * dEs:+.4f} > |dp_D11| = {abs(D11m * dEv):.4f}: post floor on the DRY side.")

# ------------------------------------------------------------------------------------------------ (D-HAR)
print("\n================ (D-HAR) m = 2 chains ================")


def chain_floor(eps_n, pi_n, v_n, deps3, fr, with_ops=True):
    S_eps = np.zeros((3, 3))
    S_pi = np.zeros(3)
    cum = 0.0
    eps_k, pi_k, v_k = eps_n.copy(), pi_n, v_n
    pat = []
    for a in fr:
        T = S_eps + a * I3
        v_new = v_k * math.exp(a * float(deps3.sum()))
        eps_trf, act_tr, Phi_tr, _ = HP.floor_op(H, eps_k + a * deps3, PMIN)
        res = K.return_map(P, eps_trf, pi_k, v_new, v_new)
        assert not res.refused, res.reason
        eps_f, act_post, Phi_post, _ = HP.floor_op(H, res.eps_e, PMIN)
        if not with_ops:
            Phi_tr, Phi_post = I3, I3
        cum += a
        S_v = v_new * cum * ONES
        Tf = Phi_tr @ T
        if res.plastic:
            ch = res.chain
            S_new = ch.b[:3, :3] @ Tf - np.outer(ch.u[:3] / ch.c, S_pi) - np.outer(ch.u[:3] * ch.Pi_v, S_v)
            S_pi = ch.w @ Tf + ((1 - ch.kappa) / ch.c) * S_pi + (1 - ch.kappa) * ch.Pi_v * S_v
        else:
            S_new = Tf
        S_eps = Phi_post @ S_new
        eps_k, pi_k, v_k = eps_f, res.pi, v_new
        pat.append(("F" if act_tr else "-") + ("P" if res.plastic else "E") + ("f" if act_post else "-"))
    el = K.elastic(P, eps_k)
    return el.ae @ S_eps, el.sig, pat, el.p


def chain_fd(eps_n, pi_n, v_n, deps3, fr, h):
    Cfd = np.zeros((3, 3))
    pats = set()
    for b in range(3):
        dp = deps3.copy()
        dp[b] += h
        dm = deps3.copy()
        dm[b] -= h
        _, sp_, pp, _ = chain_floor(eps_n, pi_n, v_n, dp, fr)
        _, sm_, pm, _ = chain_floor(eps_n, pi_n, v_n, dm, fr)
        pats.add(tuple(pp))
        pats.add(tuple(pm))
        Cfd[:, b] = (sp_ - sm_) / (2 * h)
    return Cfd, pats


fr = [0.5, 0.5]
# (D1) twice the (B-HAR) increment, halved: each half floors at the trial, returns, floors again: FPf, FPf
deps_D1 = 2.0 * deps_B
C_ch, sig_D, pat, pD = chain_floor(eps_B, pi_B, vB, deps_D1, fr)
s1 = step_floor(eps_B, pi_B, vB, deps_D1)
print(f"(D1) 2x the (B-HAR) increment (tr deps = {deps_D1.sum():.1e}), m = 2: pattern {pat}, final p = {pD:.4f} "
      f"(unsplit: {pattern(s1)}); the (B-HAR) increment itself halves to {chain_floor(eps_B, pi_B, vB, deps_B, fr)[2]}")
assert pat == ["FPf", "FPf"], pat
for h in (1e-7, 1e-8):
    Cfd, pats = chain_fd(eps_B, pi_B, vB, deps_D1, fr, h)
    assert len(pats) == 1, pats
    rep(f"(D1) chained C (S.54) 'FPf,FPf' vs FD, h = {h:.0e}", relerr(C_ch, Cfd), 3e-7 if h == 1e-7 else 3e-8)
C_mut, _, _, _ = chain_floor(eps_B, pi_B, vB, deps_D1, fr, with_ops=False)
print(f"      mutant 'Phi operators omitted from the chain' vs FD: {relerr(C_mut, Cfd):.3e}")

# (D2) a smaller increment from the same state: trial stays below the floor, the second return crosses it: -P-, -Pf
found2, fallback = None, None
for amp_v in np.linspace(1.0e-5, 2.2e-5, 13):
    for amp_s in np.linspace(1.5e-5, 9.0e-5, 16):
        cand = (amp_v / 3.0) * ONES + amp_s * nhB
        try:
            _, _, pat2, _ = chain_floor(eps_B, pi_B, vB, cand, fr)
        except AssertionError:
            continue
        if pat2 == ["-P-", "-Pf"]:
            found2 = cand
            break
        if pat2 == ["-P-", "FPf"] and fallback is None:
            fallback = cand
    if found2 is not None:
        break
if found2 is None:
    print("      pattern -P-,-Pf NOT found on the search grid from this state; using the -P-,FPf increment instead")
    found2 = fallback
assert found2 is not None
deps_D2 = found2
C_ch2, sig_D2, pat2, pD2 = chain_floor(eps_B, pi_B, vB, deps_D2, fr)
print(f"(D2) increment tr deps = {deps_D2.sum():.3e}, shear amp {np.linalg.norm(deps_D2 - deps_D2.mean()):.3e}: pattern {pat2}, final p = {pD2:.4f}")
errs2 = {}
for h in (1e-7, 1e-8):
    Cfd2, pats = chain_fd(eps_B, pi_B, vB, deps_D2, fr, h)
    assert len(pats) == 1, pats
    errs2[h] = relerr(C_ch2, Cfd2)
    rep(f"(D2) chained C (S.54) '{','.join(pat2)}' vs FD, h = {h:.0e}", errs2[h], 3e-6 if h == 1e-7 else 3e-8)
print(f"      O(h^2): err(1e-7)/err(1e-8) = {errs2[1e-7]/errs2[1e-8]:.1f} (the increment is {deps_D2.sum():.1e} + shear 2e-5)")
C_mut2, _, _, _ = chain_floor(eps_B, pi_B, vB, deps_D2, fr, with_ops=False)
print(f"      mutant 'Phi operators omitted from the chain' vs FD: {relerr(C_mut2, Cfd2):.3e}")
# the same increment unsplit: what single-step pattern is it?
s2 = step_floor(eps_B, pi_B, vB, deps_D2)
print(f"      the same increment in one step: pattern {pattern(s2)}, p_tr = {s2['p_tr']:.4f}, p_c = {s2['p_c']:.4f}, p_f = {s2['p_f']:.4f}")

HP.uninstall()
print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
