"""Sheet round 2026-10-03: the continuum (rate) form of the p-floor for O1 (sheet S.55) and the first-order
convergence of the O2/kernel split (trial floor -> NorSand return -> post floor) to it.

Rate form while the floor is active (p = -p_min and the unconstrained rate would raise p): a second mechanism
eps^f. = lam_f. P, P = delta/3, with the two consistency conditions F. = 0 and p. = 0 (Koiter):
    [ f:a^e:q + H   f:a^e:P ] [lam  ]   [ f:a^e:eps. ]
    [ P:a^e:q       P:a^e:P ] [lam_f] = [ P:a^e:eps. ]      lam, lam_f >= 0
Built from O2's primitives (BA06 energy, K2 set, alpha0 = 0 and 5), principal co-axial path on which BOTH
mechanisms stay active (wet side eta < M, deviatoric + slight expansion; checked a posteriori: min lam, min lam_f
over every RHS evaluation are > 0, so the KKT branch never switches and the RHS is smooth), integrated with Radau
at rtol 1e-10; the split backward-Euler with m = 1 ... 64 sub-increments is compared at the end of the path.
Expected: |sigma_BE(m) - sigma_rate| = O(1/m) (first order, like the plain NorSand backward Euler vs O1).
Also checked: with only the floor active (elastic + floor) and BA06 alpha0 = 0, lam_f = tr eps. exactly
(P:a^e:P = K, P:a^e:eps. = K tr eps.).
Run:  python -u floor_rate.py
"""
import math
import os
import sys
import warnings

import numpy as np
from scipy.integrate import solve_ivp

ORACLE = r"C:/Users/nmora/Documents/Github/OpenSees/.claude/worktrees/ladrunonorsand-implementation-review-7bbd75/Ladruno_files/testbed/norsand_oracle"
sys.path.insert(0, ORACLE)
sys.path.insert(0, os.path.join(ORACLE, "tests"))
warnings.simplefilter("ignore")
from conftest import K2_BASE, make_params          # noqa: E402
from o2_algo import kernel as K                   # noqa: E402

I3 = np.eye(3)
ONES = np.ones(3)
SQ23 = math.sqrt(2.0 / 3.0)
Pv = ONES / 3.0
PMIN = 50.0
ok_all = True


def floor_op(P, eps_e, pmin):
    """Pi_f (sheet S.48-S.49, BA06): same operator as floor_fd.py."""
    ev = float(eps_e.sum())
    e = eps_e - ev / 3.0
    ne = float(np.linalg.norm(e))
    es = SQ23 * ne
    el = K.elastic(P, eps_e)
    if el.p <= -pmin * (1.0 - 1e-12):
        return eps_e.copy()
    fac = 1.0 + 1.5 * P.alpha0 * es * es / P.kappa_hat
    ev_f = P.eps_v0 - P.kappa_hat * math.log(pmin / (abs(P.p0) * fac))
    return eps_e - (ev - ev_f) / 3.0


def off_corner_sig(P, p_s, eta_s, direction):
    xi = np.array(direction)
    xi -= xi.mean()
    nh = xi / np.linalg.norm(xi)
    inv = K.invariants(p_s + nh)
    z, _, _ = K.zeta_y(inv.theta, P.rho, P.zeta)
    q_s = eta_s * abs(p_s) / z
    return p_s + SQ23 * q_s * nh, inv.theta


def rhs_factory(P, track):
    def rhs(t, y, edot):
        eps_e, pi, v = y[:3], y[3], y[4]
        el = K.elastic(P, eps_e)
        inv = K.invariants(el.sig)
        fl = K.flow(P, inv, pi)
        psi, _ = K.csl(P, v, pi)
        ps, _, _ = K.pistar(P, inv.p, fl.Om, psi)
        H = -fl.Y.F_pi * SQ23 * P.h * (ps - pi) * fl.Om
        fA = fl.f_a @ el.ae
        PA = Pv @ el.ae
        G = np.array([[fA @ fl.q_a + H, fA @ Pv], [PA @ fl.q_a, PA @ Pv]])
        lam, lam_f = np.linalg.solve(G, np.array([fA @ edot, PA @ edot]))
        track.append((lam, lam_f, el.p, fl.F))
        de = edot - lam * fl.q_a - lam_f * Pv
        dpi = SQ23 * P.h * lam * (ps - pi) * fl.Om
        dv = v * float(edot.sum())
        return np.concatenate([de, [dpi, dv]])
    return rhs


def be_split(P, eps_e, pi, v, deps, m, pmin):
    for _ in range(m):
        d = deps / m
        v = v * math.exp(float(d.sum()))
        eps_trf = floor_op(P, eps_e + d, pmin)
        res = K.return_map(P, eps_trf, pi, v, v)
        assert not res.refused, res.reason
        assert res.plastic
        eps_e = floor_op(P, res.eps_e, pmin)
        pi = res.pi
    return K.elastic(P, eps_e).sig, pi


for alpha0 in (0.0, 5.0):
    print(f"\n===== BA06 K2 set, alpha0 = {alpha0}, p_min = {PMIN} =====")
    P = make_params("O2", **dict(K2_BASE, rho=0.7, rho_bar=0.8, alpha0=alpha0))
    sig0, th = off_corner_sig(P, -PMIN, 0.6 * P.M, direction=(-1.2, -0.1, 1.3))
    pi0 = K.pi_of_eta(P, -PMIN, 0.6 * P.M)
    eps0 = K.invert_elastic(P, sig0)
    v0 = 1.70
    xi = np.array([-1.2, -0.1, 1.3])
    nh = (xi - xi.mean()) / np.linalg.norm(xi - xi.mean())
    edot = nh + 0.02                                       # LOADING (along the stress deviator) with slight expansion
    T = 2.0e-3
    track = []
    rhs = rhs_factory(P, track)
    sol = solve_ivp(rhs, (0.0, T), np.concatenate([eps0, [pi0, v0]]), args=(edot,), method="Radau",
                    rtol=1e-10, atol=1e-14)
    assert sol.success, sol.message
    y = sol.y[:, -1]
    sig_rate = K.elastic(P, y[:3]).sig
    lam_min = min(t_[0] for t_ in track)
    lamf_min = min(t_[1] for t_ in track)
    p_dev = max(abs(t_[2] + PMIN) for t_ in track) / PMIN
    F_dev = max(abs(t_[3]) for t_ in track) / abs(P.p0)
    print(f"  rate (Radau rtol 1e-10, nfev {sol.nfev}): p_end = {sig_rate.mean():.6f}, pi_end = {y[3]:.6f}; over all RHS evaluations: "
          f"min lam = {lam_min:.3e}, min lam_f = {lamf_min:.3e}, max |p + p_min|/p_min = {p_dev:.1e}, max |F|/|p0| = {F_dev:.1e}")
    both = lam_min > 0 and lamf_min > 0
    print(f"  [{'OK' if both else 'FAIL'}] both mechanisms active along the whole path (KKT branch never switches)")
    errs = []
    for m in (1, 2, 4, 8, 16, 32, 64):
        sig_be, pi_be = be_split(P, eps0, pi0, v0, edot * T, m, PMIN)
        e = float(np.abs(sig_be - sig_rate).max() / np.abs(sig_rate).max())
        errs.append(e)
        print(f"  split BE m = {m:3d}: |sigma - sigma_rate|/|sigma| = {e:.3e}, |pi - pi_rate|/|pi| = {abs(pi_be - y[3])/abs(y[3]):.3e}, p = {sig_be.mean():.5f}")
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(2.0) for i in range(len(errs) - 1)]
    print("  observed orders (m -> 2m):", ", ".join(f"{o:.2f}" for o in orders))
    first_order = all(0.8 <= o <= 1.3 for o in orders[2:])
    ok_all &= first_order and both
    print(f"  [{'OK' if first_order else 'FAIL'}] split converges at first order to the simultaneous (Koiter) rate solution")
    if alpha0 == 0.0:
        el = K.elastic(P, eps0)
        PA = Pv @ el.ae
        lam_f = (PA @ edot) / (PA @ Pv)
        good = abs(lam_f - edot.sum()) < 1e-12
        print(f"  [{'OK' if good else 'FAIL'}] floor-only branch, alpha0 = 0: lam_f = tr eps. ({lam_f:.6e} vs {edot.sum():.6e}); "
              f"P:a^e:P = {PA @ Pv:.3f} = K = {-el.p / P.kappa_hat:.3f}")
        ok_all &= good

print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
