"""Reference (oracle) for the LadrunoCSTPair macro-element (ADR 70 P4a).

Pure-numpy assembly of the disjoint 2-triangle F-bar-Patch macro-element
(dSNPO §15.1.9): four nodes n1..n4, triangles T1=(n1,n2,n3) / T2=(n1,n3,n4)
(split along the n1-n3 diagonal), one centroid GP per triangle, and the
shared patch dilatation

    J-bar = v_patch / V_patch = (J1*V1 + J2*V2) / (V1 + V2)      (eq 15.36)
    F-bar_e = (J-bar / J_e)^{1/2} F_e                            (eq 15.34, n=2)

The material is driven at F-bar_e; the spatial configuration (gradients g_e,
volume weight dv_e = J_e*A_e*t) stays on the ACTUAL F_e (dSNPO Remark 15.2).

Exact consistent tangent (the P4a novelty). Two algebraically identical
forms are implemented and cross-checked against each other:

  (a) the substitute-row form the C++ element uses through the shared kernel:
        K += g_e^T a g_e dv_e  +  p_e (g-bar - g_e) dv_e,
        p_e[a,i] = sum_j g_e[a,j] q_ij,
        q_ij = (1/2) sum_p a_ijpp - (1/2) sigma_ij           (eq 15.10, n=2)
        g-bar = sum_s (v_s / v_patch) g_s   (8-DOF padded rows — the patch
        volume-weighted spatial dilatation-gradient operator)
  (b) the dSNPO eq 15.37/15.38 split:
        K^(ee) += g_e^T a g_e dv_e + (v_e/v_p - 1) g_e^T q g_e dv_e
        K^(es) += (v_s/v_p) g_e^T q g_s dv_e            (s = the other triangle)

Both are generally UNSYMMETRIC and stress-proportional in the coupling
(q -> (1/2)c:I != 0 at F=I, but the (g-bar - g_e) row is displacement-
dependent; the FD gate must run at finite stress to exercise everything).

No OpenSees dependency — numpy + the LogStrain2D material oracle only.
"""
from __future__ import annotations

import numpy as np

import logstrain2d_reference as l2

I2 = np.eye(2)
I3 = np.eye(3)

# canonical mildly-distorted reference quad (CCW), split along n1-n3
REF4 = np.array([[0.0, 0.0], [2.0, 0.1], [2.2, 1.8], [0.1, 1.9]])
TRIS = [(0, 1, 2), (0, 2, 3)]


def _t3_dNdX(Xtri):
    """Constant reference gradients dN/dX [2x3] and area for one T3."""
    d = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])   # dN/dxi
    J = d @ Xtri                                          # dX/dxi (2x2)
    detJ = np.linalg.det(J)
    assert detJ > 0.0, "triangle winding"
    dNdX = np.linalg.inv(J) @ d                           # [i][a]
    return dNdX, 0.5 * detJ


def _material(F2, K, G, kind="elastic", plane="strain"):
    step = l2.plane_stress_step if plane == "stress" else l2.plane_strain_step
    res = step(F2, I3.copy(), I3.copy(), kind, K=K, G=G)
    s = res["sigma_voigt"]
    sig = np.array([[s[0], s[2]], [s[2], s[1]]])
    return sig, res["c2"]


def _tri_state(X4, U4, tri):
    """Per-triangle kinematics: (padded dNdX [2x4], A, F [2x2], J)."""
    idx = list(tri)
    dNdX, A = _t3_dNdX(X4[idx])
    dN4 = np.zeros((2, 4))
    dN4[:, idx] = dNdX                                    # zero-padded 4-node rows
    F = I2 + U4.T @ dN4.T                                 # F_ij = d_ij + sum_a U_ai dN4_ja
    return dN4, A, F, np.linalg.det(F)


def assemble_pair(u, K, G, thickness=1.0, X=None, form="row"):
    """Return (f_int [8], Kmat [8x8]) for the pair at nodal displacement u [8].

    form="row"   : substitute-row coupling (what the C++ element assembles)
    form="dsnpo" : explicit eq 15.37/15.38 (ve/vp - 1) / (vs/vp) split
    """
    if X is None:
        X = REF4
    U = u.reshape(4, 2)
    states = [_tri_state(X, U, t) for t in TRIS]

    Vp = sum(A for (_, A, _, _) in states)
    vp = sum(J * A for (_, A, _, J) in states)
    Jbar = vp / Vp
    assert Jbar > 0.0 and all(J > 0.0 for (_, _, _, J) in states)

    # per-triangle spatial rows g[a,j] (4-node padded) + deformed volumes
    gs, dvs, sigs, a4s = [], [], [], []
    for (dN4, A, F, J) in states:
        Fbar = np.sqrt(Jbar / J) * F
        sig, c2 = _material(Fbar, K, G)
        g = dN4.T @ np.linalg.inv(F)                      # g[a][j], actual F
        a4 = c2 - np.einsum("il,jk->ijkl", sig, I2)
        gs.append(g); dvs.append(J * A * thickness)
        sigs.append(sig); a4s.append(a4)

    vtot = sum(dvs)
    gbar = sum((dvs[e] / vtot) * gs[e] for e in range(2))  # patch row operator

    f = np.zeros(8)
    Km = np.zeros((8, 8))
    for e in range(2):
        g, dv, sig, a4 = gs[e], dvs[e], sigs[e], a4s[e]
        for a in range(4):
            for i in range(2):
                f[2 * a + i] += np.dot(sig[i], g[a]) * dv
        Km += np.einsum("aj,ijkl,bl->aibk", g, a4, g).reshape(8, 8) * dv

        q = 0.5 * np.einsum("ijpp->ij", a4) - 0.5 * sig
        p = g @ q.T                                       # p[a][i]
        if form == "row":
            Km += np.einsum("ai,bk->aibk", p, gbar - g).reshape(8, 8) * dv
        elif form == "dsnpo":
            s = 1 - e
            Km += (dvs[e] / vtot - 1.0) * np.einsum(
                "ai,bk->aibk", p, g).reshape(8, 8) * dv
            Km += (dvs[s] / vtot) * np.einsum(
                "ai,bk->aibk", p, gs[s]).reshape(8, 8) * dv
        else:
            raise ValueError(form)
    return f, Km


def internal_force(u, K, G, **kw):
    return assemble_pair(u, K, G, **kw)[0]


def fd_tangent(u, K, G, h=1e-7, **kw):
    """d f_int / d u by central differences (unsymmetric-aware)."""
    Kfd = np.zeros((8, 8))
    for d in range(8):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (internal_force(up, K, G, **kw)
                     - internal_force(um, K, G, **kw)) / (2 * h)
    return Kfd
