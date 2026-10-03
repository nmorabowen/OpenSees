"""numpy replica of one mortar pair's slave-side force (normal + Coulomb/cohesion), committed state
gT0=gpT=lamT=lamN=0 (step 1 from u=0), its FD and the shipped analytic tangent block."""
import numpy as np
import os
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "proto_c1_mortar.py"), encoding="utf-8").read()
_ns = {}
exec(compile(_src[: _src.index("def grid_nodes")], "proto_c1_mortar", "exec"), _ns)
mortar_pair, facet_normal = _ns["mortar_pair"], _ns["facet_normal"]
RD = np.array([0.0, 0.0, 1.0])


def rmap(gT, N, kt, mu, coh):
    t = kt * gT; cap = mu * N + coh if N > 0 else 0.0; nr = np.linalg.norm(t)
    if nr <= cap: return -t, False
    return -cap * t / nr, True


def force(Xs, Xm, kn, kt, mu, coh, us):
    r = mortar_pair(Xs, 4, Xm, 4, RD)
    if r is None: return np.zeros((4, 3)), None
    n = facet_normal(Xm, 4, RD); D, g = r["D"], r["g"]; a = D.sum(1)
    p = np.zeros(4); tf = np.zeros((4, 3)); slip = [None] * 4
    for I in range(4):
        if a[I] <= 1e-300: continue
        pr = kn * g[I] / a[I]; p[I] = min(pr, 0.0)
        if p[I] < 0:
            rr = D[I] @ us              # master fixed: r = sum_J D_IJ u_s,J
            gT = (rr - (rr @ n) * n) / a[I]
            tf[I], slip[I] = rmap(gT, -p[I], kt, mu, coh)
    f = -(D @ p)[:, None] * n[None, :] + D @ tf
    return f, dict(D=D, a=a, p=p, n=n, us=us, slip=slip, area=r["area"])


def kan(info, kn, kt, mu, coh, consistent=False):
    D, a, p, n, us = info["D"], info["a"], info["p"], info["n"], info["us"]
    K = np.zeros((12, 12)); Pt = np.eye(3) - np.outer(n, n)
    for I in range(4):
        if p[I] >= 0: continue
        bb = np.outer(D[I], D[I]) / a[I]
        K += kn * np.kron(bb, np.outer(n, n))
        rr = D[I] @ us; gT = (rr - (rr @ n) * n) / a[I]; t = kt * gT
        N = -p[I]; cap = mu * N + coh; nr = np.linalg.norm(t)
        if nr <= cap: Kss = kt * Pt
        else:
            nh = t / nr; Kss = cap * kt / nr * (Pt - np.outer(nh, nh))
            if consistent: Kss = Kss + np.outer(-mu * kn * nh, n)
        K += np.kron(bb, Kss)
    return K
