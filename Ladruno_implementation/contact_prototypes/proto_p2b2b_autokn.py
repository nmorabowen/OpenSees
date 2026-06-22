"""ADR-39 P2b-2b ORACLE — `-kn auto` penalty-stiffness auto-sizing, frictionless.

Oracle-first (the fork discipline): pin the reduction arithmetic + the scaling
laws of the auto-kn heuristic in pure numpy BEFORE the C++ resolves it inside
LadrunoContactHandler::handle(). No OpenSees, no build.

DESIGN (user-gated iteration 18): kn is sourced GENERICALLY from the owning solid
element's *initial stiffness matrix* — NOT a LadrunoBrick-specific material call
(materialPointers[] is private). The handler, at handle() time, has the Domain, so
per (slave node, master segment) pair it:
  1. auto-detects the owning solid element (the one whose node set ⊇ the segment
     face nodes),
  2. pulls Element::getInitialStiff() (24×24 for a hex; base virtual, any solid),
  3. reduces to a scalar interface stiffness via the segment normal n:

        kn = f_si · mean_over_face_nodes( nᵀ K_block_node n )                (★)

     where K_block_node is the 3×3 diagonal block of getInitialStiff() at that
     node's translational DOFs. f_si = 0.10 (LS-DYNA SLSFAC default).

WHY (★) is dimensionally the LS-DYNA 26.14a form `f·K·A²/V`: for a 3D solid of
characteristic size L and modulus E, a nodal diagonal stiffness K_diag ~ E·L, and
A²/V ~ L⁴/L³ = L, so E·A²/V ~ E·L ~ K_diag. The element-stiffness route absorbs
the A²/V geometry into the assembled matrix exactly and the modulus into D, so (★)
is the same scaling with the constant folded into f_si — and it needs no modulus
extraction, no element-type coupling, no vanilla edit. The SOFT Courant floor
(26.15) is DEFERRED to P2b-2c.

This oracle validates (non-circularly):
  T1  the reduction (★) on a synthetic SPD block matrix == hand arithmetic;
  T2  on a REAL assembled trilinear hex: kn ∝ E   (modulus scaling);
  T3  on a REAL assembled hex: kn ∝ L   (size scaling, K_diag ~ E·L);
  T4  penetration g = P/kn is a SMALL fraction of the element size for a
      physically-scaled nodal load (the heuristic is usable, not lock/blow-up);
  T5  normal-direction sensitivity: an oblique normal reduces a different block
      combination (n picks the directional stiffness), still positive & bounded;
  T6  mesh-refinement: halving L keeps g/L bounded (size-objective contact).

Run: <pythoncore-3.12> proto_p2b2b_autokn.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)

F_SI = 0.10   # penalty scale factor (LS-DYNA SLSFAC default); MUST match the C++


# ---------------------------------------------------------------- hex stiffness
# Standard isoparametric trilinear 8-node hexahedron, 2x2x2 Gauss, isotropic
# linear elastic. This is exactly `LadrunoBrick -formulation std` getInitialStiff(),
# so the C++ test can compare kn to this oracle to ~1% on a std-formulation brick.
# (bbar/eas differ in the matrix but the SCALING laws T2/T3/T6 hold regardless.)

# natural-coord corner ordering (the conventional brick ordering)
_HEX_XI = np.array([
    [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
    [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
], dtype=float)
_G = 1.0 / np.sqrt(3.0)   # 2-pt Gauss
_GP = np.array([[a, b, c] for a in (-_G, _G) for b in (-_G, _G) for c in (-_G, _G)])


def _shape_grad(xi):
    """trilinear N and dN/dxi (8x3) at natural coord xi=[r,s,t]."""
    r, s, t = xi
    sgn = _HEX_XI
    N = 0.125 * (1 + sgn[:, 0] * r) * (1 + sgn[:, 1] * s) * (1 + sgn[:, 2] * t)
    dN = np.zeros((8, 3))
    dN[:, 0] = 0.125 * sgn[:, 0] * (1 + sgn[:, 1] * s) * (1 + sgn[:, 2] * t)
    dN[:, 1] = 0.125 * sgn[:, 1] * (1 + sgn[:, 0] * r) * (1 + sgn[:, 2] * t)
    dN[:, 2] = 0.125 * sgn[:, 2] * (1 + sgn[:, 0] * r) * (1 + sgn[:, 1] * s)
    return N, dN


def _D_iso(E, nu):
    """6x6 isotropic elastic tangent (Voigt [xx,yy,zz,xy,yz,zx])."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    return D


def hex_stiffness(L=1.0, E=1.0, nu=0.3, x0=(0.0, 0.0, 0.0)):
    """assemble the 24x24 initial stiffness of an LxLxL hex at corner x0.
    returns (K[24x24], X[8x3] node coords in the conventional ordering)."""
    x0 = np.asarray(x0, float)
    X = x0 + 0.5 * L * (_HEX_XI + 1.0)   # map [-1,1]^3 corners to [0,L]^3 + x0
    D = _D_iso(E, nu)
    K = np.zeros((24, 24))
    for gp in _GP:
        N, dN = _shape_grad(gp)
        J = dN.T @ X                       # 3x3 Jacobian (dx/dxi)
        detJ = np.linalg.det(J)
        dNx = dN @ np.linalg.inv(J)        # 8x3 dN/dx
        B = np.zeros((6, 24))
        for a in range(8):
            bx, by, bz = dNx[a]
            c = 3 * a
            B[0, c + 0] = bx
            B[1, c + 1] = by
            B[2, c + 2] = bz
            B[3, c + 0] = by; B[3, c + 1] = bx
            B[4, c + 1] = bz; B[4, c + 2] = by
            B[5, c + 0] = bz; B[5, c + 2] = bx
        K += (B.T @ D @ B) * detJ           # gauss weight = 1 for 2-pt
    return K, X


# --------------------------------------------------------------- the reduction
def auto_kn(K, faceNodeIdx, n, f_si=F_SI):
    """(★) kn = f_si · mean_over_face_nodes( nᵀ K_block_node n ).
    K            : element initial stiffness (3*nNode square)
    faceNodeIdx  : list of element-local node indices on the contact segment
    n            : unit segment normal (3,)
    """
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    acc = 0.0
    for a in faceNodeIdx:
        c = 3 * a
        Kb = K[c:c + 3, c:c + 3]            # 3x3 diagonal block
        acc += n @ Kb @ n                   # directional nodal stiffness (force/len)
    return f_si * acc / len(faceNodeIdx)


def top_face(X):
    """element-local node indices of the +z (top) face of the conventional hex."""
    zmax = X[:, 2].max()
    return [a for a in range(len(X)) if abs(X[a, 2] - zmax) < 1e-9 * max(1.0, zmax + 1)]


# ----------------------------------------------------------------------- tests
def test_reduction_arithmetic():
    """T1 — (★) on a synthetic block-diagonal SPD matrix == hand value."""
    # 2 face nodes, hand-chosen 3x3 blocks; n = +z picks the (2,2) entry.
    K = np.zeros((6, 6))
    K[0:3, 0:3] = np.diag([10.0, 20.0, 30.0])
    K[3:6, 3:6] = np.diag([40.0, 50.0, 70.0])
    n = np.array([0.0, 0.0, 1.0])
    kn = auto_kn(K, [0, 1], n, f_si=0.10)
    expect = 0.10 * (30.0 + 70.0) / 2.0     # mean of the zz entries
    assert abs(kn - expect) < 1e-12, (kn, expect)
    # off-diagonal coupling: rotate n into x-y, nᵀKn must pick the rotated combo
    K2 = np.zeros((3, 3)); K2[0:3, 0:3] = np.array([[10., 2., 0.], [2., 20., 0.], [0., 0., 5.]])
    Kc = np.zeros((3, 3)); Kc[:] = K2
    n2 = np.array([1.0, 1.0, 0.0]); n2 = n2 / np.linalg.norm(n2)
    kn2 = auto_kn(Kc, [0], n2, f_si=1.0)
    assert abs(kn2 - (n2 @ K2 @ n2)) < 1e-12, kn2
    print(f"[T1 reduction] kn={kn:.4f} (expect {expect:.4f}); oblique={kn2:.4f} OK")


def test_modulus_scaling():
    """T2 — kn ∝ E on a real assembled hex (top face, n=+z)."""
    K1, X = hex_stiffness(L=1.0, E=100.0, nu=0.3)
    K2, _ = hex_stiffness(L=1.0, E=300.0, nu=0.3)
    face = top_face(X)
    n = np.array([0, 0, 1.0])
    kn1 = auto_kn(K1, face, n)
    kn2 = auto_kn(K2, face, n)
    assert abs(kn2 / kn1 - 3.0) < 1e-9, (kn1, kn2, kn2 / kn1)
    assert kn1 > 0
    print(f"[T2 E-scaling] kn(E=100)={kn1:.3f} kn(E=300)={kn2:.3f} ratio={kn2/kn1:.4f} OK")


def test_size_scaling():
    """T3 — kn ∝ L (K_diag ~ E·L): doubling L doubles kn."""
    n = np.array([0, 0, 1.0])
    Ka, Xa = hex_stiffness(L=1.0, E=100.0, nu=0.3)
    Kb, Xb = hex_stiffness(L=2.0, E=100.0, nu=0.3)
    kna = auto_kn(Ka, top_face(Xa), n)
    knb = auto_kn(Kb, top_face(Xb), n)
    assert abs(knb / kna - 2.0) < 1e-9, (kna, knb, knb / kna)
    print(f"[T3 L-scaling] kn(L=1)={kna:.3f} kn(L=2)={knb:.3f} ratio={knb/kna:.4f} OK")


def test_penetration_fraction():
    """T4 — g = P/kn is a SMALL fraction of L for a physically-scaled load.
    Load a slave with the force one element would transmit at strain eps:
    P ~ sigma·A = (E·eps)·L². Then g/L = P/(kn·L). With kn ~ f_si·E·L and
    P ~ E·eps·L², g/L ~ eps/f_si — bounded, mesh-independent, and small for eps«f_si."""
    E, nu, L = 200.0, 0.3, 1.5
    K, X = hex_stiffness(L=L, E=E, nu=nu)
    n = np.array([0, 0, 1.0])
    kn = auto_kn(K, top_face(X), n)
    eps = 1e-3                       # 0.1% strain — a working contact pressure
    P = E * eps * L * L              # nodal-scale contact force
    g = P / kn
    frac = g / L
    assert 0 < frac < 0.1, frac      # penetration < 10% of element size
    print(f"[T4 penetration] kn={kn:.2f} P={P:.4f} g={g:.3e} g/L={frac:.3e} (<0.1) OK")


def test_oblique_normal():
    """T5 — an oblique normal yields a positive, bounded kn (picks directional stiff)."""
    K, X = hex_stiffness(L=1.0, E=100.0, nu=0.3)
    face = top_face(X)
    n_axis = np.array([0, 0, 1.0])
    n_obl = np.array([1.0, 0.0, 2.0]); n_obl = n_obl / np.linalg.norm(n_obl)
    kn_axis = auto_kn(K, face, n_axis)
    kn_obl = auto_kn(K, face, n_obl)
    assert kn_obl > 0, kn_obl
    # the oblique value mixes shear+normal nodal stiffness — must stay same order
    assert 0.2 * kn_axis < kn_obl < 2.0 * kn_axis, (kn_axis, kn_obl)
    print(f"[T5 oblique] kn(axis)={kn_axis:.3f} kn(oblique)={kn_obl:.3f} OK")


def test_mesh_objectivity():
    """T6 — refine the mesh (halve L) → g/L stays bounded (size-objective).
    Repeats T4's physically-scaled load at two element sizes; the penetration
    FRACTION g/L is invariant (it depends on eps/f_si, not on L)."""
    E, nu, eps = 200.0, 0.3, 1e-3
    fracs = []
    for L in (2.0, 1.0, 0.5):
        K, X = hex_stiffness(L=L, E=E, nu=nu)
        kn = auto_kn(K, top_face(X), n=np.array([0, 0, 1.0]))
        P = E * eps * L * L
        fracs.append((P / kn) / L)
    spread = max(fracs) - min(fracs)
    assert spread < 1e-6, fracs       # g/L identical across refinement
    print(f"[T6 mesh-obj] g/L over L=[2,1,0.5] = {[f'{f:.3e}' for f in fracs]} spread={spread:.1e} OK")


if __name__ == "__main__":
    print("=== ADR-39 P2b-2b `-kn auto` oracle ===")
    test_reduction_arithmetic()
    test_modulus_scaling()
    test_size_scaling()
    test_penetration_fraction()
    test_oblique_normal()
    test_mesh_objectivity()
    print("ALL P2b-2b AUTO-KN ORACLE CHECKS PASSED")
