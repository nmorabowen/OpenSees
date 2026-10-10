"""ADR-155 ORACLE -- pile-contact R0.5: -augment, -maxGap, -gapOffset/-adjust (numpy only).

Build-free closed forms that tests/test_adr155_pile_contact_r05.py checks in compiled OpenSees.

  T1  interference pressure (-gapOffset -delta) on two elastic bars in series with the mortar
      penalty spring: sigma = delta / (2L/E + 1/epsN); Uzawa (held-load augmentation) -> delta*E/(2L).
  T2  N-2: a LINEAR penalty tie solved in n load steps. With the commit-cycle Uzawa update after
      every step (the shipped default) the answer depends on n; with the update off (-augment never)
      it does not, and the error against the exact bond is O(1/eps).
  T3  -adjust arithmetic: the shift gbar_shift = g/a + (-g0/a0) is an EXACT 0.0 at the reference
      (g == g0, a == a0 bitwise), for any doubles -- the zero-pressure start is not a round-off
      residue. -gapOffset alone is a pure additive shift of gbar.
  T4  -maxGap geometry on a closed circle (the test's 16-facet tube in a 24-facet ring): the plane
      distance of every NEAR pair (angular overlap) is < 0.01*R, of every ANTIPODAL pair whose
      projection overlaps (the N-1 defect -- the clip accepts it) ~ 2R; any maxGap in between keeps
      exactly the near pairs. The window is two orders of magnitude wide.

Run: python proto_adr155_r05.py
"""
import math
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------------------------------- T1
def t1_interference():
    print("T1 interference pressure (two bars + penalty, gap shifted by -delta)")
    E, L, A, delta, eps = 2.0e4, 1.0, 1.0, 1.0e-3, 1.0e6
    kb = E * A / L
    # u1 = bottom bar's top (z = 1-), u2 = top bar's bottom (z = 1+); the far ends are clamped.
    # shifted gap g = u2 - u1 - delta (< 0 = penetration); ALM pressure p = lam + eps*g (<= 0);
    # compressive force t = -A*p. Equilibrium: kb*u1 + t = 0, kb*u2 - t = 0 (linear while closed).

    def solve(lam):
        K = np.array([[kb + A * eps, -A * eps], [-A * eps, kb + A * eps]])
        f = np.array([A * lam - A * eps * delta, -A * lam + A * eps * delta])
        u = np.linalg.solve(K, f)
        g = u[1] - u[0] - delta
        return -A * (lam + eps * g), g              # (t, g)

    t, _ = solve(0.0)
    closed = delta * A / (2 * L / E + 1.0 / eps)
    check("penalty sigma == delta/(2L/E+1/eps)", abs(t - closed) < 1e-12 * closed,
          f"{t:.12g} vs {closed:.12g}")
    lam = 0.0
    for _ in range(60):                               # one Uzawa update per (held-load) commit
        t, g = solve(lam)
        lam = min(0.0, lam + eps * g)
    t, _ = solve(lam)
    exact = delta * E * A / (2 * L)
    check("Uzawa limit sigma == delta*E/2L", abs(t - exact) < 1e-9 * exact, f"{t:.12g} vs {exact:.12g}")


# ---------------------------------------------------------------------------------------- T2
def t2_step_count():
    print("T2 N-2: linear penalty tie, n load steps, Uzawa per commit on/off")
    k1, k2, P = 2.0e4, 2.0e4, 20.0          # bar below the tie (k1, clamped), bar above (k2, loaded)

    def run(eps, n, augment):
        # dofs: a (master side of the tie), b (slave side), c (tip). tie: r = b - a -> 0.
        lam = 0.0
        for i in range(1, n + 1):
            f = np.array([0.0, 0.0, P * i / n])
            K = np.array([[k1 + eps, -eps, 0.0], [-eps, eps + k2, -k2], [0.0, -k2, k2]])
            f = f + np.array([lam, -lam, 0.0])   # tie traction t = lam + eps*r
            u = np.linalg.solve(K, f)
            if augment:
                lam += eps * (u[1] - u[0])        # no clamp (equality bond)
        return u[2]

    exact = P / k1 + P / k2
    eps = 1.0e5
    never = [run(eps, n, False) for n in (1, 2, 5)]
    commit = [run(eps, n, True) for n in (1, 2, 5)]
    check("never: n-independent", max(abs(x - never[0]) for x in never) < 1e-15 * abs(never[0]),
          f"{never}")
    check("commit: n-dependent", abs(commit[2] - commit[0]) > 1e-4 * abs(commit[0]), f"{commit}")
    errs = [abs(run(e, 1, False) - exact) for e in (1e5, 1e6, 1e7, 1e8)]
    ratios = [a / b for a, b in zip(errs, errs[1:])]
    check("never: error O(1/eps)", all(9.5 < r < 10.5 for r in ratios), f"ratios {ratios}")


# ---------------------------------------------------------------------------------------- T3
def t3_adjust_exact_zero():
    print("T3 -adjust: the reference shift is an exact 0.0; -gapOffset is additive")
    rng = np.random.default_rng(155)
    ok = True
    for _ in range(100000):
        g0 = float(rng.normal()) * 10.0 ** float(rng.uniform(-8, 1))
        a0 = float(rng.uniform(1e-6, 10.0))
        gbar = g0 / a0 + (-(g0 / a0))        # what LadrunoContactFE::mortarActive computes at u=0
        if (gbar * a0) / a0 != 0.0 or gbar != 0.0:
            ok = False
            break
    check("100000 random (g0, a0): shifted gap == 0.0 bitwise", ok)
    g, a, off = 0.37, 1.9, -1e-3
    shifted = (g / a + 0.0 + off) * a
    check("offset only: gbar -> gbar + g0", abs(shifted / a - (g / a + off)) < 1e-15)


# ---------------------------------------------------------------------------------------- T4
def t4_maxgap_window():
    print("T4 -maxGap: near vs antipodal plane distances on the test's closed cylinder")
    R, nti, nto = 0.5, 16, 24

    def facets(n):
        out = []
        for j in range(n):
            t0, t1 = 2 * math.pi * j / n, 2 * math.pi * (j + 1) / n
            p0 = np.array([R * math.cos(t0), R * math.sin(t0)])
            p1 = np.array([R * math.cos(t1), R * math.sin(t1)])
            c = 0.5 * (p0 + p1)
            nrm = c / np.linalg.norm(c)
            out.append((p0, p1, c, nrm))
        return out

    M, S = facets(nti), facets(nto)
    near, far = [], []
    for (m0, m1, mc, mn) in M:
        t = np.array([-mn[1], mn[0]])
        lo, hi = sorted((float((m0 - mc) @ t), float((m1 - mc) @ t)))
        for (s0, s1, sc, sn) in S:
            # the kernel clip: project the slave chord onto the master line along the master normal
            a, b = sorted((float((s0 - mc) @ t), float((s1 - mc) @ t)))
            if min(hi, b) - max(lo, a) <= 1e-12:
                continue                     # no overlap -> the clip rejects it anyway
            d = abs(float((sc - mc) @ mn))
            (near if float(sn @ mn) > 0.0 else far).append(d)
    check("antipodal overlapping pairs exist (the N-1 defect)", len(far) > 0, f"{len(far)} pairs")
    check("near pairs are within 0.01 R", max(near) < 0.01 * R, f"max near {max(near):.3e}")
    check("antipodal pairs are ~2R away", min(far) > 1.9 * R, f"min far {min(far):.3f}")
    check("the test's maxGap = 0.1 splits them", max(near) < 0.1 < min(far))


if __name__ == "__main__":
    t1_interference()
    t2_step_count()
    t3_adjust_exact_zero()
    t4_maxgap_window()
    print(f"\n{'ALL PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
