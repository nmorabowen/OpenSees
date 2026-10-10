"""ADR-157 ORACLE -- mortar friction path state per (slave node, facet PAIR), not per slave node.

The kernel math is UNCHANGED (LadrunoFrictionKernel::frictionReturnMap, the C3.1 weighted slip,
the C3.3 lambda_T offset trick). What changes is WHERE the path state lives. ADR-41 C3.1-C3.3
stored gpT / gT0 / engaged / lambdaT once per GLOBAL slave node; every (slave facet, master facet)
adapter that integrates the node runs its OWN return map on its OWN local slip

    gbarT_I^pair = P_n(pair) [ sum_J D_IJ u_s,J - sum_K M_IK u_m,K ] / a_I^pair

and wrote the node's trial state, last writer wins (LEDGER_quirks "Mortar friction committed slip is
last-writer-wins", the C3.1 gate MAJOR-1). This oracle mirrors LadrunoContactFE::addMortarFriction
line for line on a 1-D mortar interface embedded in 3-D (segments in the x-z plane) and runs BOTH
storage layouts:

    'node'  key = slave node                  (the shipped C3.1-C3.3 layout)
    'pair'  key = (slave node, sf, mf)        (ADR-157)

Gates:
  T1  FLAT non-matching interface, UNIFORM slip, 3 slave/master alignments x 2 sweep orders:
      stick force = -kt*s*L, full-slip force = -cap*L (cohesion, Coulomb, tauMax cap) for BOTH
      layouts (the pairs agree when the field is uniform -- this is why the fence held on the
      matched battery). Pinned so the new layout is not a regression on the analytic patch.
  T2  FLAT non-matching, NON-UNIFORM (linear) slip, step 1 then a step-2 increment, both sweep
      orders: 'pair' is order-INDEPENDENT and in full slip gives f_K = -cap*a_K exactly;
      'node' is order-DEPENDENT (the committed slip is the last writer's).
  T3  CURVED interface (a regular polygon, the slave polygon rotated half a facet -- every slave
      facet straddles a master vertex, the R1 pile geometry in 2-D), rigid LATERAL slave slip in
      two steps, cohesion-only: 'pair' friction tractions stay in each pair's tangent plane
      (|t.n| = 0) and the total lateral force equals the analytic sum over master facets of
      cap*L_f*dir_f; 'node' leaks a normal component |t.n| ~ cap*sin(facet angle) at step 2 (a
      neighbour's lambda_T / gpT read in the wrong plane) and misses the analytic force.
  T4  one facet pair per node (the matched single-facet battery): the two layouts are BIT-identical.

Run: python proto_adr157_mortar_pair_friction.py   (numpy only; exit 1 on any FAIL)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' -- ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ---------------------------------------------------------------- the shipped friction kernel
def friction_cap(N, mu, c, tmax):
    capC = mu * N + c if N > 0.0 else 0.0
    return tmax if (tmax > 0.0 and tmax < capC) else capC


def return_map(gTeff, gpT, N, kt, mu, c=0.0, tmax=0.0):
    """LadrunoFrictionKernel::frictionReturnMap: returns (tFric = -tT, gpTtrial)."""
    tTtr = kt * (gTeff - gpT)
    cap = friction_cap(N, mu, c, tmax)
    nrm = np.linalg.norm(tTtr)
    if cap <= 0.0:
        return np.zeros(3), gTeff.copy()
    if nrm <= cap:
        return -tTtr, gpT.copy()
    nh = tTtr / nrm
    return -cap * nh, gpT + (nrm - cap) / kt * nh


# ---------------------------------------------------------------- 1-D mortar pair integration
GP, GW = np.polynomial.legendre.leggauss(6)


def pair_ops(Xs, Xm, outward):
    """Clip slave segment Xs (2x3) against master segment Xm (2x3) by normal projection onto the
    master line; return (D, M, n) or None. n = the master normal oriented toward `outward`."""
    e = Xm[1] - Xm[0]
    Lm = np.linalg.norm(e)
    e = e / Lm
    n = np.array([-e[2], 0.0, e[0]])
    if n @ outward < 0.0:
        n = -n
    xi = [(Xs[i] - Xm[0]) @ e / Lm for i in range(2)]      # master param of each slave node
    # xi(eta) = xi0 + (xi1-xi0) eta must lie in [0,1]
    a, b = 0.0, 1.0
    d = xi[1] - xi[0]
    if abs(d) < 1e-14:
        if not (0.0 <= xi[0] <= 1.0):
            return None
    else:
        e0, e1 = (0.0 - xi[0]) / d, (1.0 - xi[0]) / d
        a, b = max(a, min(e0, e1)), min(b, max(e0, e1))
    if b - a <= 1e-12:
        return None
    Ls = np.linalg.norm(Xs[1] - Xs[0])
    D = np.zeros((2, 2))
    M = np.zeros((2, 2))
    for g, w in zip(GP, GW):
        eta = a + (b - a) * (g + 1.0) / 2.0
        wj = w * (b - a) / 2.0 * Ls
        Ns = np.array([1.0 - eta, eta])
        x = xi[0] + d * eta
        Nm = np.array([1.0 - x, x])
        D += wj * np.outer(Ns, Ns)
        M += wj * np.outer(Ns, Nm)
    return D, M, n


class Store:
    """The Domain-owned friction slots, keyed per 'node' or per 'pair'."""

    def __init__(self, layout):
        self.layout, self.s = layout, {}

    def slot(self, node, sf, mf):
        k = node if self.layout == "node" else (node, sf, mf)
        if k not in self.s:
            z = np.zeros(3)
            self.s[k] = dict(gpT=z.copy(), gpTt=z.copy(), gT0=z.copy(), eng=False,
                             lamT=z.copy(), lamTt=z.copy())
        return self.s[k]

    def commit(self, augment=True):
        for st in self.s.values():
            st["gpT"] = st["gpTt"].copy()
            if augment:
                st["lamT"] = st["lamTt"].copy()


def sweep(pairs, us, um, store, p, kt, mu, c, tmax, order):
    """One residual sweep (addMortarFriction per pair, in `order`). Returns slave nodal forces
    (nS x 3), master nodal forces (nM x 3), and the per-pair friction tractions [(n, tF[2])]."""
    fs = np.zeros_like(us)
    fm = np.zeros_like(um)
    trac = []
    for ip in order:
        sf, mf, sn, mn, D, M, n = pairs[ip]
        tF = np.zeros((2, 3))
        for I in range(2):
            a = D[I].sum()
            if a <= 1e-300:
                continue
            r = D[I] @ us[sn] - M[I] @ um[mn]
            gbarT = (r - (r @ n) * n) / a
            st = store.slot(sn[I], sf, mf)
            if not st["eng"]:
                st["gT0"], st["eng"] = gbarT.copy(), True
            gTeff = (gbarT - st["gT0"]) + st["lamT"] / kt
            tf, gpt = return_map(gTeff, st["gpT"], p, kt, mu, c, tmax)
            st["gpTt"], st["lamTt"] = gpt, -tf
            tF[I] = tf
        for K in range(2):
            fs[sn[K]] += D[K] @ tF
        for L in range(2):
            fm[mn[L]] -= M[:, L] @ tF
        trac.append((n, tF))
    return fs, fm, trac


def engage(pairs, nS, nM, store, kt, mu, c, tmax):
    """Step 0: the as-built configuration (u = 0) is the first in-contact evaluation, so every
    slot captures gT0 = 0 there -- exactly what the first Newton iterate does in the C++."""
    sweep(pairs, np.zeros((nS, 3)), np.zeros((nM, 3)), store, P, kt, mu, c, tmax,
          range(len(pairs)))
    store.commit()


def build_pairs(Xs, segS, Xm, segM, outward):
    pairs = []
    for sf, (i, j) in enumerate(segS):
        for mf, (k, l) in enumerate(segM):
            r = pair_ops(Xs[[i, j]], Xm[[k, l]], outward)
            if r is not None:
                pairs.append((sf, mf, np.array([i, j]), np.array([k, l])) + r)
    return pairs


def flat(ns, nm, off, L=1.0):
    """Slave nodes on z=0 over [0,L] with ns segments; master over [-off, L+off'] clipped so the
    master covers the slave; master nodes shifted by `off` (non-matching)."""
    xs = np.linspace(0.0, L, ns + 1)
    xm = np.concatenate([[-0.5], np.linspace(off, L - 0.37 * off, nm + 1)[1:-1], [L + 0.5]])
    Xs = np.array([[x, 0.0, 0.0] for x in xs])
    Xm = np.array([[x, 0.0, 0.0] for x in xm])
    segS = [(i, i + 1) for i in range(ns)]
    segM = [(i, i + 1) for i in range(len(xm) - 1)]
    return Xs, segS, Xm, segM


KT, P = 1.0e4, 100.0
CONES = {"cohesion": (0.0, 50.0, 0.0), "Coulomb": (0.4, 0.0, 0.0),
         "Coulomb+c capped": (0.4, 30.0, 45.0)}

# ------------------------------------------------------------------------------------------ T1
print("T1  flat non-matching, UNIFORM slip: analytic stick/slip, both layouts, any alignment")
for cone, (mu, c, tmax) in CONES.items():
    cap = friction_cap(P, mu, c, tmax)
    worst = 0.0
    for (ns, nm, off) in [(5, 3, 0.11), (7, 4, 0.23), (4, 4, 0.0)]:
        Xs, segS, Xm, segM = flat(ns, nm, off)
        pairs = build_pairs(Xs, segS, Xm, segM, np.array([0, 0, 1.0]))
        for layout in ("node", "pair"):
            for order in (list(range(len(pairs))), list(range(len(pairs)))[::-1]):
                store = Store(layout)
                um = np.zeros((len(Xm), 3))
                engage(pairs, len(Xs), len(Xm), store, KT, mu, c, tmax)
                for s, expect in [(0.2 * cap / KT, -KT * 0.2 * cap / KT),   # stick
                                  (3.0 * cap / KT, -cap), (5.0 * cap / KT, -cap)]:   # slip, slip
                    us = np.zeros((len(Xs), 3)); us[:, 0] = s
                    fs, fm, _ = sweep(pairs, us, um, store, P, KT, mu, c, tmax, order)
                    store.commit(augment=(s > cap / KT))   # no Uzawa on the stick step
                    worst = max(worst, abs(fs[:, 0].sum() - expect * 1.0) / cap,
                                abs(fs[:, 0].sum() + fm[:, 0].sum()) / cap)
    check(f"{cone}: sum f_x == -kt*s*L (stick) / -cap*L (slip), self-equilibrated", worst < 1e-12,
          f"max rel err {worst:.1e}")

# ------------------------------------------------------------------------------------------ T2
print("T2  flat non-matching, LINEAR slip field, two steps: order (in)dependence")
mu, c, tmax = 0.0, 50.0, 0.0
cap = c
Xs, segS, Xm, segM = flat(6, 4, 0.13)
pairs = build_pairs(Xs, segS, Xm, segM, np.array([0, 0, 1.0]))
aK = np.zeros(len(Xs))
for (_, _, sn, _, D, _, _) in pairs:
    aK[sn] += D.sum(axis=1)
res = {}
for layout in ("node", "pair"):
    for oname, order in [("fwd", list(range(len(pairs)))), ("rev", list(range(len(pairs)))[::-1])]:
        store = Store(layout)
        um = np.zeros((len(Xm), 3))
        engage(pairs, len(Xs), len(Xm), store, KT, mu, c, tmax)
        # step 1: a steep slip field (all in slip); step 2: a much smaller uniform increment
        us = np.zeros((len(Xs), 3)); us[:, 0] = (2.0 + 400.0 * Xs[:, 0]) * cap / KT
        sweep(pairs, us, um, store, P, KT, mu, c, tmax, order)
        store.commit()
        us[:, 0] += 1.5 * cap / KT
        fs, _, _ = sweep(pairs, us, um, store, P, KT, mu, c, tmax, order)
        res[(layout, oname)] = fs[:, 0].copy()
d_pair = np.abs(res[("pair", "fwd")] - res[("pair", "rev")]).max() / cap
d_node = np.abs(res[("node", "fwd")] - res[("node", "rev")]).max() / cap
e_pair = np.abs(res[("pair", "fwd")] + cap * aK).max() / cap
e_node = max(np.abs(res[("node", o)] + cap * aK).max() for o in ("fwd", "rev")) / cap
check("'pair': step-2 forces independent of the sweep order", d_pair < 1e-13, f"{d_pair:.1e}")
check("'pair': full slip => f_K = -cap*a_K exactly", e_pair < 1e-12, f"{e_pair:.1e}")
check("'node' (shipped): order-DEPENDENT and off the analytic slip force (the defect)",
      d_node > 1e-3 and e_node > 1e-3, f"order diff {d_node:.2e}, analytic err {e_node:.2e}")

# ------------------------------------------------------------------------------------------ T3
print("T3  CURVED interface (polygon, slave rotated half a facet), lateral slip, cohesion only")
NF, R = 24, 0.5
th_m = np.arange(NF) * 2 * np.pi / NF
th_s = th_m + np.pi / NF
Xm = np.array([[R * np.cos(t), 0.0, R * np.sin(t)] for t in th_m])     # the pile skin (master)
Rs = R                                                                    # slave vertices on R too
Xs = np.array([[Rs * np.cos(t), 0.0, Rs * np.sin(t)] for t in th_s])    # the soil hole (slave)
segM = [(i, (i + 1) % NF) for i in range(NF)]
segS = [(i, (i + 1) % NF) for i in range(NF)]
pairs = []
for sf, (i, j) in enumerate(segS):
    mid = 0.5 * (Xs[i] + Xs[j])
    for mf, (k, l) in enumerate(segM):
        mm = 0.5 * (Xm[k] + Xm[l])
        if np.linalg.norm(mid - mm) > 2.5 * R * np.sin(np.pi / NF) * 2:   # neighbours only
            continue
        r = pair_ops(Xs[[i, j]], Xm[[k, l]], -mm)   # outward: slave lies outside the skin
        if r is not None:
            pairs.append((sf, mf, np.array([i, j]), np.array([k, l])) + r)
check("every slave node is shared by >= 3 facet pairs (the curved shared-node case)",
      min(sum(1 for pr in pairs if I in pr[2]) for I in range(NF)) >= 3)
out = {}
for layout in ("node", "pair"):
    store = Store(layout)
    um = np.zeros((NF, 3))
    engage(pairs, NF, NF, store, KT, 0.0, c, 0.0)
    for step, s in enumerate([40.0 * c / KT, 48.0 * c / KT]):   # every pair slips (min |sin| = sin 7.5deg)
        us = np.zeros((NF, 3)); us[:, 0] = s              # rigid lateral slip of the hole ring
        fs, fm, trac = sweep(pairs, us, um, store, P, KT, 0.0, c, 0.0, range(len(pairs)))
        store.commit()
        leak = max(abs(tF[I] @ n) for n, tF in trac for I in range(2)) / c
        out[(layout, step)] = (fs[:, 0].sum(), fs[:, 2].sum(), leak)
# analytic: each pair slips; its traction is c along -P_n x_hat/|P_n x_hat|, integrated over the
# pair's overlap (sum of D rows = the overlap length) -- a sum over master facets.
Fx = Fz = 0.0
for (_, _, _, _, D, _, n) in pairs:
    tdir = np.array([1.0, 0, 0]) - n[0] * n
    tdir /= np.linalg.norm(tdir)
    Fx -= c * D.sum() * tdir[0]
    Fz -= c * D.sum() * tdir[2]
for step in (0, 1):
    fxp, fzp, leakp = out[("pair", step)]
    fxn, fzn, leakn = out[("node", step)]
    check(f"step {step + 1} 'pair': |t.n| == 0 and F == analytic sum_f c*L_f*dir_f",
          leakp < 1e-12 and abs(fxp - Fx) < 1e-9 * abs(Fx) and abs(fzp - Fz) < 1e-9 * abs(Fx),
          f"leak {leakp:.1e}, Fx {fxp:.6f} vs {Fx:.6f}")
    print(f"        step {step + 1} 'node' (shipped): leak |t.n|/c = {leakn:.3e}, "
          f"Fx {fxn:.6f} vs analytic {Fx:.6f} (rel {abs(fxn - Fx) / abs(Fx):.2e})")
check("'node' (shipped) leaks a normal friction component at step 2 (the defect)",
      out[("node", 1)][2] > 1e-2, f"{out[('node', 1)][2]:.3e} (sin(360/NF) = {np.sin(2*np.pi/NF):.3f})")

# ------------------------------------------------------------------------------------------ T4
print("T4  one facet pair per node: the two layouts are bit-identical")
Xs, segS, Xm, segM = flat(1, 1, 0.0)
Xm = np.array([[0.0, 0, 0], [1.0, 0, 0]])
segM = [(0, 1)]
pairs = build_pairs(Xs, segS, Xm, segM, np.array([0, 0, 1.0]))
got = {}
for layout in ("node", "pair"):
    store = Store(layout)
    um = np.zeros((2, 3))
    engage(pairs, 2, 2, store, KT, 0.3, c, 0.0)
    hist = []
    for s in (0.3, 1.7, 2.9, 2.0):
        us = np.array([[s * c / KT, 0, 0], [1.1 * s * c / KT, 0, 0]])
        fs, _, _ = sweep(pairs, us, um, store, P, KT, 0.3, c, 0.0, [0])
        store.commit()
        hist.append(fs.copy())
    got[layout] = np.array(hist)
check("single pair: node == pair bit for bit", np.array_equal(got["node"], got["pair"]))

print(f"\n{'ALL PASS' if _fails == 0 else f'{_fails} FAIL(S)'}")
sys.exit(1 if _fails else 0)
