"""ADR-159 oracle -- the smoothed mortar contact law (`-smoothN <g0>`, `-smoothT <r>`) and its tangent.

Pure numpy, no OpenSees. Run: python proto_adr159_smooth_normal.py  (exit 0 = every gate passed).

NORMAL LAW (per slave node I of a mortar pair; x = -(lambda_I + epsN*gbar_I) = epsN*z with z the
penetration; S = epsN*g0 is the band in pressure units, g0 > 0 the band half-width, a LENGTH):

    x <= -S        P = 0                       (open: no pressure and never tension)
    |x| <  S       P = (x + S)^2 / (4 S)       (quadratic onset over the gap band -g0 < gbar < g0)
    x >= S         P = x                       (the shipped penalty ramp, unchanged)

  p_I = -P (the shipped sign). P is C1 with a Lipschitz P' = (x+S)/(2S), P >= 0, and
  P -> max(0, x) uniformly as g0 -> 0 (max error S/4 at x = 0).

FRICTION ONSET (part of -smoothN): the tangential traction is scaled by the C1 weight
    chi(x) = t^2 (3 - 2 t),  t = clamp((x + S)/(2S), 0, 1)
so a cohesive bond (cap = mu*P + c, c > 0) fades out over the band instead of jumping from c to 0
when the node lifts off (the shipped law's residual JUMP). chi = 1 for x >= S.

STICK/SLIP CORNER (-smoothT r, 0 < r < 1): the return-map magnitude min(rho, cap), rho = ||tT*||,
is rounded over |rho - cap| < delta = r*cap by phi = rho - (rho - cap + delta)^2/(4 delta) (C1),
with the slip that keeps T = kt (gT - gpT_trial).

TANGENT (tang = -dR/du with the shipped frozen-D,M,n mortar operator, b = [D, -M], dz/du_B = -b_IB n/a):
  normal:   kn * P'(x) * b_IA b_IB / a  n (x) n                    (P' replaces the 0/1 active mask)
  friction: chi * K_ss(dN/dz = kn P')  b_IA b_IB / a               (K_ss the shipped/smoothed block)
          + kn * chi'(x) * b_IA b_IB / a  tF (x) n                  (the onset coupling; NON-symmetric,
                                                                    assembled only under -consistanttan)
  K_ss (rounded corner) = kt[phi'(rho) nh (x) nh + phi/rho (P_t - nh (x) nh)] - (dphi/dcap) mu kn P' nh (x) n

Gates:
  G1  scalar law: C1 at the band edges, P >= 0, P == ramp outside the band, sup|P - ramp| = S/4.
  G2  one-node residual vs analytic tangent, central FD at 4000 random states (open / band / closed,
      stick / blend / slip, cohesion / Coulomb / Tresca cap), with and without -smoothT; relative
      error <= 1e-6 (states within FD reach of the UNSMOOTHED stick/slip switch are skipped when
      r = 0: that switch is a genuine kink, which is what -smoothT removes).
  G3  Newton on a shared crease node (two facet pairs, normals +-20 deg, lateral load ramped to full
      slip with lift-off of one flank, force control): the smoothed law converges in <= 8 iterations
      per step for Coulomb and mixed cones at g0 = 1e-6..1e-4 (pure cohesion reported, see gate_g3).
      A one-node toy does NOT reproduce the shipped stall (shipped converges here too): the stall
      needs the mesh coupling of the pile deck, so the binary gates are on the R3 deck itself.
  G4  g0 -> 0 on a closed case: the smoothed answer converges to the shipped one (error O(g0), exact
      once the band no longer reaches the penetration).
"""
import sys

import numpy as np


# ----------------------------------------------------------------------------------- the law
def law(x, S):
    """P, dP/dx, chi, dchi/dx; S <= 0 => the shipped ramp (chi = 1 iff x > 0)."""
    if S <= 0.0:
        return (x, 1.0, 1.0, 0.0) if x > 0.0 else (0.0, 0.0, 0.0, 0.0)
    if x <= -S:
        return 0.0, 0.0, 0.0, 0.0
    if x >= S:
        return x, 1.0, 1.0, 0.0
    t = (x + S) / (2.0 * S)
    return (x + S) ** 2 / (4.0 * S), t, t * t * (3.0 - 2.0 * t), 6.0 * t * (1.0 - t) / (2.0 * S)


def corner(rho, cap, r):
    """phi, dphi/drho, dphi/dcap of min(rho, cap) rounded over |rho - cap| < delta = r*cap (r = 0 =>
    the shipped min). dphi/dcap is TOTAL: delta moves with cap."""
    delta = r * cap
    if delta <= 0.0:
        return (rho, 1.0, 0.0) if rho <= cap else (cap, 0.0, 1.0)
    if rho <= cap - delta:
        return rho, 1.0, 0.0
    if rho >= cap + delta:
        return cap, 0.0, 1.0
    e = rho - cap + delta
    dphi_dcap = e / (2 * delta) + r * (e * e / (4 * delta * delta) - e / (2 * delta))
    return rho - e * e / (4 * delta), 1.0 - e / (2 * delta), dphi_dcap


def cap_of(N, mu, c, tmax):
    capC = mu * N + c if N > 0.0 else 0.0
    capped = tmax > 0.0 and tmax < capC
    return (tmax if capped else capC), capped


def return_map(gT, gpT, N, kt, mu, c, tmax, r=0.0):
    """Applied (negated) traction tF and the trial slip (LadrunoFrictionKernel, shipped or smooth)."""
    tr = kt * (gT - gpT)
    cap, _ = cap_of(N, mu, c, tmax)
    rho = np.linalg.norm(tr)
    if cap <= 0.0:
        return np.zeros(3), gT.copy()
    phi, _, _ = corner(rho, cap, r)
    if rho <= cap - r * cap or rho == 0.0:
        return -tr, gpT.copy()
    nh = tr / rho
    return -phi * nh, gpT + (rho - phi) / kt * nh


def kss_block(gT, gpT, n, N, dNdz, kt, mu, c, tmax, r=0.0, consistent=True):
    tr = kt * (gT - gpT)
    cap, capped = cap_of(N, mu, c, tmax)
    Pt = np.eye(3) - np.outer(n, n)
    if cap <= 0.0:
        return np.zeros((3, 3))
    rho = np.linalg.norm(tr)
    if rho <= cap - r * cap or rho == 0.0:
        return kt * Pt
    phi, dr, dc = corner(rho, cap, r)
    nh = tr / rho
    K = phi * kt / rho * (Pt - np.outer(nh, nh)) + kt * dr * np.outer(nh, nh)
    if consistent:
        dcap = 0.0 if (capped or N <= 0.0) else mu
        K += np.outer(-dc * dcap * dNdz * nh, n)
    return K


# ------------------------------------------------- one slave node vs a rigid master (b = [a])
def node_residual(u, n, a, kn, S, gap0, gpT, kt, mu, c, tmax, r=0.0):
    """The contact force APPLIED to one slave node (area a): +a(P n + chi tF); tangent = -dR/du."""
    x = -kn * (gap0 + n @ u)
    P, _, chi, _ = law(x, S)
    gT = u - (n @ u) * n
    R = a * P * n                                    # -(D p) n with p = -P, D = a
    if P > 0.0:
        tF, _ = return_map(gT, gpT, P, kt, mu, c, tmax, r)
        R = R + a * chi * tF
    return R


def node_tangent(u, n, a, kn, S, gap0, gpT, kt, mu, c, tmax, r=0.0):
    x = -kn * (gap0 + n @ u)
    P, dP, chi, dchi = law(x, S)
    gT = u - (n @ u) * n
    K = a * kn * dP * np.outer(n, n)                 # b b / a = a
    if P > 0.0:
        tF, _ = return_map(gT, gpT, P, kt, mu, c, tmax, r)
        K += a * chi * kss_block(gT, gpT, n, P, kn * dP, kt, mu, c, tmax, r)
        K += a * kn * dchi * np.outer(tF, n)
    return K


def node_commit(u, n, kn, S, gap0, gpT, kt, mu, c, tmax, r=0.0):
    x = -kn * (gap0 + n @ u)
    P = law(x, S)[0]
    if P <= 0.0:
        return gpT
    return return_map(u - (n @ u) * n, gpT, P, kt, mu, c, tmax, r)[1]


# ----------------------------------------------------------------------------------- gates
def gate_g1():
    S = 2.0
    xs = np.linspace(-3 * S, 3 * S, 20001)
    P = np.array([law(x, S)[0] for x in xs])
    dP = np.array([law(x, S)[1] for x in xs])
    ch = np.array([law(x, S)[2] for x in xs])
    ok = P.min() >= 0.0
    ok &= np.allclose(P[xs >= S], xs[xs >= S]) and np.all(P[xs <= -S] == 0.0)
    for e in (-S, S):                                # C1 at the edges (left/right limits)
        h = 1e-9
        for k, tol in ((0, 1e-8), (1, 1e-8), (2, 1e-8), (3, 1e-6)):
            ok &= abs(law(e - h, S)[k] - law(e + h, S)[k]) < tol
    for e in (0.7, 1.3):                             # rounded corner: C1 at cap -+ delta
        h = 1e-9
        for k in range(3):
            ok &= abs(corner(e - h, 1.0, 0.3)[k] - corner(e + h, 1.0, 0.3)[k]) < 1e-7
    err = np.max(np.abs(P - np.maximum(xs, 0.0)))
    ok &= abs(err - S / 4.0) < 1e-6
    ok &= np.all(np.diff(P) >= -1e-15) and np.all(np.diff(ch) >= -1e-15)   # monotone
    ok &= np.max(np.abs(np.gradient(P, xs) - dP)) < 1e-3
    print("G1 scalar law: C1 edges (normal + corner), P>=0, ramp outside band, sup|P-ramp| = %.4f"
          " (S/4 = %.4f): %s" % (err, S / 4, "PASS" if ok else "FAIL"))
    return ok


def gate_g2(rng):
    worst, nchk, nblend = 0.0, 0, 0
    for k in range(4000):
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        a, kn, kt = rng.uniform(0.2, 2.0), 10 ** rng.uniform(5, 7), 10 ** rng.uniform(4, 7)
        g0 = 10 ** rng.uniform(-5, -3)
        S = kn * g0
        r = rng.choice([0.0, 0.05, 0.3])
        gap0 = rng.uniform(-2.5, 2.5) * g0           # open / in band / closed
        mu = rng.choice([0.0, 0.3, 1.0]); c = rng.choice([0.0, 0.0, 0.3 * S])
        if mu == 0.0 and c == 0.0:
            c = 0.5 * S
        tmax = rng.choice([0.0, 0.0, 0.5 * (mu * S + c) + 1e-9])
        gpT = rng.normal(size=3) * g0 * 0.3; gpT -= (gpT @ n) * n
        u = rng.normal(size=3) * g0 * 0.4
        args = (n, a, kn, S, gap0, gpT, kt, mu, c, tmax, r)
        x = -kn * (gap0 + n @ u)
        P = law(x, S)[0]
        rho = np.linalg.norm(kt * ((u - (n @ u) * n) - gpT))
        capC = mu * P + c if P > 0 else 0.0
        cap = tmax if (tmax > 0 and tmax < capC) else capC
        h = 1e-7 * g0
        reach = 1e3 * h * (kt + mu * kn)
        edges = [cap] if r == 0.0 else [cap * (1 - r), cap * (1 + r)]
        if P > 0 and min(abs(rho - e) for e in edges) < reach:
            continue                                 # FD straddles a corner / blend edge
        if tmax > 0 and abs(capC - tmax) < 1e3 * h * mu * kn:
            continue                                 # the min(mu N + c, tmax) switch (not smoothed)
        if P > 0 and r > 0 and abs(rho - cap) < r * cap:
            nblend += 1
        K = node_tangent(u, *args)
        Kfd = np.zeros((3, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = h
            Kfd[:, j] = -(node_residual(u + e, *args) - node_residual(u - e, *args)) / (2 * h)
        scale = max(np.abs(Kfd).max(), a * kn * 1e-3)
        worst = max(worst, np.abs(K - Kfd).max() / scale)
        nchk += 1
    ok = worst <= 1e-6 and nchk > 3000 and nblend > 50
    print("G2 node tangent vs central FD: %d states (%d in the rounded corner), worst relative "
          "error %.2e: %s" % (nchk, nblend, worst, "PASS" if ok else "FAIL"))
    return ok


AL = np.radians(20.0)
CREASE = [np.array([np.sin(AL), 0.0, np.cos(AL)]), np.array([-np.sin(AL), 0.0, np.cos(AL)])]


def crease_newton(S, r, mu=0.0, c=200.0, kd=1e4, kn=1e7, kt=1e6, steps=20, maxit=40):
    """One slave node shared by two facet pairs of a crease (area 0.5 each), on a spring kd, pressed
    by 1e3 and pushed laterally to 3e3 (past full slip on the up-slope flank, which lifts off)."""
    u = np.zeros(3); gp = [np.zeros(3), np.zeros(3)]; its = []
    for s in range(1, steps + 1):
        f = np.array([s / steps * 3e3, 0.0, -1e3])
        for it in range(1, maxit + 1):
            R = f - kd * u; K = kd * np.eye(3)
            for k in range(2):
                R = R + node_residual(u, CREASE[k], 0.5, kn, S, 0.0, gp[k], kt, mu, c, 0.0, r)
                K = K + node_tangent(u, CREASE[k], 0.5, kn, S, 0.0, gp[k], kt, mu, c, 0.0, r)
            if np.linalg.norm(R) < 1e-4:
                break
            u = u + np.linalg.solve(K, R)
        else:
            its.append(None)
            return its
        its.append(it - 1)
        gp = [node_commit(u, CREASE[k], kn, S, 0.0, gp[k], kt, mu, c, 0.0, r) for k in range(2)]
    return its


def gate_g3():
    """Gated: Coulomb (mu = 0.3) and mixed (mu = 0.5, c = 100) cones at three bands, r = 0.1.
    Reported, not gated: pure cohesion (c = 200). Losing a cohesive bond across the band is a local
    SOFTENING of stiffness ~ c*a*chi'_max = 0.75*c*a/g0; when that rivals the normal stiffness
    epsN*a the force-controlled path has a limit point inside the step and no Newton variant (or
    law) can follow it -- a brittle debond, physics rather than numerics. Shipped is reported too."""
    rows, ok = [], True
    for g0 in (1e-6, 1e-5, 1e-4):
        for mu, c in ((0.3, 0.0), (0.5, 100.0)):
            its = crease_newton(1e7 * g0, 0.1, mu=mu, c=c)
            ok &= None not in its and max(its) <= 8
            rows.append(max(its) if None not in its else None)
    coh = [crease_newton(1e7 * g0, 0.1, mu=0.0, c=200.0) for g0 in (1e-6, 1e-5, 1e-4)]
    coh = [max(i) if None not in i else "snap@%d" % len(i) for i in coh]
    shp = [crease_newton(0.0, 0.0, mu=mu, c=c) for mu, c in ((0.3, 0.0), (0.5, 100.0), (0.0, 200.0))]
    shp = [max(i) if None not in i else "fail" for i in shp]
    print("G3 crease node Newton to full slip (-smoothN g0 -smoothT 0.1): max its per step %s (gated "
          "<= 8); pure cohesion %s (reported); shipped %s (reported): %s"
          % (rows, coh, shp, "PASS" if ok else "FAIL"))
    return ok


def gate_g4():
    """Closed case (pressed only, penetration 1e-5): smoothed -> shipped as g0 -> 0."""
    def solve(S):
        n = np.array([0.0, 0.0, 1.0]); u = np.zeros(3); kd, f = 1e5, np.array([10.0, 0.0, -100.0])
        for _ in range(60):
            R = f - kd * u + node_residual(u, n, 1.0, 1e7, S, 0.0, np.zeros(3), 1e6, 0.3, 0.0, 0.0)
            if np.linalg.norm(R) < 1e-10:
                break
            K = kd * np.eye(3) + node_tangent(u, n, 1.0, 1e7, S, 0.0, np.zeros(3), 1e6, 0.3, 0.0, 0.0)
            u = u + np.linalg.solve(K, R)
        return u
    u0 = solve(0.0)
    errs = [np.linalg.norm(solve(1e7 * g) - u0) / np.linalg.norm(u0) for g in (1e-4, 3e-5, 1e-5, 3e-6)]
    ok = all(e2 < e1 for e1, e2 in zip(errs, errs[1:])) and errs[-1] < 1e-12 and errs[0] > 1e-3
    print("G4 g0 = 1e-4, 3e-5, 1e-5, 3e-6 (penetration 1e-5): relative error vs shipped %s: %s"
          % (", ".join("%.1e" % e for e in errs), "PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    rng = np.random.default_rng(159)
    res = [gate_g1(), gate_g2(rng), gate_g3(), gate_g4()]
    sys.exit(0 if all(res) else 1)
