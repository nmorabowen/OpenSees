# hex20_reference.py — sympy/numpy oracle for the H20 serendipity kernel
# (SRC/element/ladrunoBrick/LadrunoHex20Shape.h), ADR 72 P0.
#
# INDEPENDENCE CONTRACT (R2-auditable): the shape functions here are built
# SYMBOLICALLY from the serendipity DEFINITION (ADR 72 §3.1; Felippa AFEM
# §18.3) and all derivatives come from sympy.diff — nothing is transcribed
# from the C++ header. The node/GP ordering tables are hardcoded from the
# PRIMARY sources (SRC/element/UP-ucsd/shp3dv.cpp:51-79 L/M/N node pattern,
# :251-277 brcshl RA/SA/TA; SRC/recorder/Ladruno_ElementResults.h:1261-1289),
# not from the header.
#
# Run standalone (all P0 gates, exact-rational where possible):
#     python tests/hex20_reference.py
# Imported by tests/test_hex20_kernel_cpp.py for the C++ cross-check.
#
# Gates implemented (ADR 72 §6 P0):
#   G1 partition of unity, Kronecker delta, sum of gradients = 0   (symbolic)
#   G2 degree-2 completeness on the reference element              (symbolic)
#   G3 quadrature degree: 27-pt and Irons 14-pt exact to degree 5,
#      8-pt exact to degree 3 (random-polynomial vs exact integral) (symbolic)
#   G4 Jacobian/detJ on a distorted straight-edged hex, detJ > 0   (numeric)
#   G5 K rank: 54 @ 27-pt; 48 @ 8-pt with EXACTLY 6 spurious modes
#      (catalogued to hex20_spurious_modes.txt); 14-pt rank verdict (numeric)
#   G6 consistent-mass diagonals > 0                               (numeric)
#   G7 EXACT lumping fractions on the reference cube:
#      row-sum corner = -1/8, mid-edge = +1/6 (the FORBIDDEN lump);
#      HRZ corner = 7/248, mid-edge = 2/31 (all positive)          (exact)

import sys
import numpy as np
import sympy as sp

# ----------------------------------------------------------------- orderings --
# Node pattern (shp3dv.cpp:51-67 comment + L/M/N tables :75-79): corners 0-3
# lower CCW from (-1,-1,-1), 4-7 upper; mid-edges 8-11 lower ring (0-1,1-2,2-3,
# 3-0), 12-15 upper ring, 16-19 vertical (0-4,1-5,2-6,3-7).
NODE_XI = [
    (-1, -1, -1), (+1, -1, -1), (+1, +1, -1), (-1, +1, -1),
    (-1, -1, +1), (+1, -1, +1), (+1, +1, +1), (-1, +1, +1),
    (0, -1, -1), (+1, 0, -1), (0, +1, -1), (-1, 0, -1),
    (0, -1, +1), (+1, 0, +1), (0, +1, +1), (-1, 0, +1),
    (-1, -1, 0), (+1, -1, 0), (+1, +1, 0), (-1, +1, 0),
]

# element edges (corner pairs) in mid-edge node order 8..19
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
         (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]

# 27-pt rule, brcshl order (shp3dv.cpp:251-277: natural coords = 2*{RA,SA,TA}*g)
_G = sp.sqrt(sp.Rational(3, 5))
_W5, _W8 = sp.Rational(5, 9), sp.Rational(8, 9)


def _w27(pt):
    return sp.prod([_W8 if c == 0 else _W5 for c in pt])


_GP27_PATTERN = [
    (-1, -1, -1), (+1, -1, -1), (+1, +1, -1), (-1, +1, -1),   # corners lower
    (-1, -1, +1), (+1, -1, +1), (+1, +1, +1), (-1, +1, +1),   # corners upper
    (0, -1, -1), (+1, 0, -1), (0, +1, -1), (-1, 0, -1),       # lower edge mids
    (0, -1, +1), (+1, 0, +1), (0, +1, +1), (-1, 0, +1),       # upper edge mids
    (-1, -1, 0), (+1, -1, 0), (+1, +1, 0), (-1, +1, 0),       # vertical mids
    (+1, 0, 0), (0, +1, 0), (0, 0, +1),                       # faces +x +y +z
    (-1, 0, 0), (0, -1, 0), (0, 0, -1),                       # faces -x -y -z
    (0, 0, 0),                                                 # centroid
]
GP27 = [(p[0] * _G, p[1] * _G, p[2] * _G, _w27(p)) for p in _GP27_PATTERN]

# 8-pt rule (recorder Hexahedron_GaussLegendre_2: lexicographic, zeta fastest)
_A = 1 / sp.sqrt(3)
GP8 = [(sx * _A, sy * _A, sz * _A, sp.Integer(1))
       for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]

# Irons 14-pt degree-5 rule (Irons 1971): 6 face points (+-b,0,0)... weight B6,
# 8 corner points (+-c,+-c,+-c) weight C8. Values verified below by the G3
# degree-5 exactness gate (self-validating, so a typo here cannot pass).
_B14 = sp.sqrt(sp.Rational(19, 30))
_C14 = sp.sqrt(sp.Rational(19, 33))
_WB14 = sp.Rational(320, 361)
_WC14 = sp.Rational(121, 361)
GP14 = ([(s * _B14, 0, 0, _WB14) for s in (+1, -1)]
        + [(0, s * _B14, 0, _WB14) for s in (+1, -1)]
        + [(0, 0, s * _B14, _WB14) for s in (+1, -1)]
        + [(sx * _C14, sy * _C14, sz * _C14, _WC14)
           for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])

# ------------------------------------------------- symbolic shape functions --
XI, ETA, ZETA = sp.symbols("xi eta zeta")
_XYZ = (XI, ETA, ZETA)


def _shape_symbolic():
    """Serendipity H20 from the definition (ADR 72 §3.1 / Felippa AFEM §18.3)."""
    N = []
    for (xa, ya, za) in NODE_XI:
        if 0 not in (xa, ya, za):  # corner
            n = sp.Rational(1, 8) * (1 + XI * xa) * (1 + ETA * ya) * (1 + ZETA * za) \
                * (XI * xa + ETA * ya + ZETA * za - 2)
        elif xa == 0:
            n = sp.Rational(1, 4) * (1 - XI**2) * (1 + ETA * ya) * (1 + ZETA * za)
        elif ya == 0:
            n = sp.Rational(1, 4) * (1 - ETA**2) * (1 + XI * xa) * (1 + ZETA * za)
        else:
            n = sp.Rational(1, 4) * (1 - ZETA**2) * (1 + XI * xa) * (1 + ETA * ya)
        N.append(sp.expand(n))
    return N


N_SYM = _shape_symbolic()
DN_SYM = [[sp.diff(n, v) for v in _XYZ] for n in N_SYM]  # derivatives via sympy

# fast numeric evaluators (numpy)
_N_FUN = sp.lambdify(_XYZ, N_SYM, "numpy")
_DN_FUN = sp.lambdify(_XYZ, DN_SYM, "numpy")


def shape_numeric(xi):
    """N[20], dN[20][3] at natural point xi (floats)."""
    N = np.array(_N_FUN(*xi), dtype=float)
    dN = np.array(_DN_FUN(*xi), dtype=float)
    return N, dN


# ------------------------------------------------------- geometry / assembly --
def straight_edge_nodes(corners):
    """20-node coordinate array from 8 corners: mid-edge nodes at TRUE midpoints."""
    X = np.zeros((20, 3))
    X[:8] = np.asarray(corners, dtype=float)
    for k, (i, j) in enumerate(EDGES):
        X[8 + k] = 0.5 * (X[i] + X[j])
    return X


# distorted straight-edged hex shared with the C++ driver (keep in sync!)
CORNERS_DISTORT = [
    (0.00, 0.00, 0.00), (2.10, -0.10, 0.15), (2.30, 2.05, -0.10), (-0.15, 1.90, 0.10),
    (0.10, 0.20, 2.20), (2.00, 0.00, 1.90), (2.20, 2.30, 2.35), (0.05, 2.10, 2.00),
]
X_DISTORT = straight_edge_nodes(CORNERS_DISTORT)
X_CUBE = straight_edge_nodes(NODE_XI[:8])  # reference biunit cube

# shared natural-point samples for the C++ cross-check (keep in sync!)
XI_SAMPLES = [
    (0.0, 0.0, 0.0),
    (0.3, -0.2, 0.7),
    (-0.9, 0.55, -0.35),
    (1.0, 1.0, 1.0),
    (0.123, 0.456, -0.789),
]


def jacobian(X, dN):
    """J[i][j] = dx_i/dxi_j; returns (J, detJ)."""
    J = X.T @ dN  # (3x20)(20x3)
    return J, float(np.linalg.det(J))


def cart_grad(X, dN):
    J, detJ = jacobian(X, dN)
    dNdx = dN @ np.linalg.inv(J)  # dNdx[a,i] = dN[a,j] Jinv[j,i]
    return dNdx, detJ


def fill_B(dNdx):
    """6x60 Voigt {xx,yy,zz,xy,yz,zx}, engineering shear."""
    B = np.zeros((6, 60))
    for a in range(20):
        bx, by, bz = dNdx[a]
        c = 3 * a
        B[0, c] = bx
        B[1, c + 1] = by
        B[2, c + 2] = bz
        B[3, c], B[3, c + 1] = by, bx
        B[4, c + 1], B[4, c + 2] = bz, by
        B[5, c], B[5, c + 2] = bz, bx
    return B


def iso_D(E=100.0, nu=0.25):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * G
        D[3 + i, 3 + i] = G
    return D


def gauss_table(rule):
    return np.array([[float(sp.N(c, 20)) for c in row] for row in rule])


GP27_F, GP8_F, GP14_F = (gauss_table(r) for r in (GP27, GP8, GP14))


def stiffness(X, D, gp):
    K = np.zeros((60, 60))
    for x, y, z, w in gp:
        _, dN = shape_numeric((x, y, z))
        dNdx, detJ = cart_grad(X, dN)
        B = fill_B(dNdx)
        K += w * detJ * (B.T @ D @ B)
    return K


def consistent_mass(X, rho, gp):
    M = np.zeros((60, 60))
    for x, y, z, w in gp:
        N, dN = shape_numeric((x, y, z))
        _, detJ = jacobian(X, dN)
        NN = rho * w * detJ * np.outer(N, N)
        for d in range(3):
            M[d::3, d::3] += NN
    return M


def hrz_lump(M):
    """Hinton-Rock-Zienkiewicz diagonal scaling (per translational direction)."""
    mdiag = np.zeros(60)
    for d in range(3):
        idx = np.arange(d, 60, 3)
        S = M[idx, idx].sum()
        m = M[idx, :].sum()
        mdiag[idx] = (m / S) * M[idx, idx]
    return mdiag


# ------------------------------------------------------------------- gates --
def gate_G1_basics():
    assert sp.simplify(sum(N_SYM) - 1) == 0, "partition of unity FAILED"
    for j, v in enumerate(_XYZ):
        assert sp.simplify(sum(dn[j] for dn in DN_SYM)) == 0, "grad-sum FAILED"
    for b, xb in enumerate(NODE_XI):
        sub = dict(zip(_XYZ, xb))
        for a, n in enumerate(N_SYM):
            val = n.subs(sub)
            assert val == (1 if a == b else 0), f"Kronecker FAILED N_{a}({b})={val}"
    print("G1 PASS  partition of unity / Kronecker / grad-sum")


def gate_G2_completeness():
    import random
    rng = random.Random(72)
    mono = [1, XI, ETA, ZETA, XI * ETA, ETA * ZETA, ZETA * XI, XI**2, ETA**2, ZETA**2]
    u = sum(sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)) * m for m in mono)
    interp = sum(N_SYM[a] * u.subs(dict(zip(_XYZ, NODE_XI[a]))) for a in range(20))
    assert sp.simplify(interp - u) == 0, "degree-2 completeness FAILED"
    print("G2 PASS  degree-2 completeness (reference element)")


def gate_G3_quadrature():
    import random
    rng = random.Random(272)
    # Exactness notion differs by rule TYPE:
    #  - 27-pt (3x3x3) and 8-pt (2x2x2) are TENSOR-PRODUCT Gauss rules: exact for
    #    each variable's polynomial degree up to `deg` INDEPENDENTLY (product=True),
    #    so the probe spans monomials with max(i,j,k) <= deg (e.g. x^5 y^5 z^5).
    #  - Irons 14-pt is a NON-product rule exact only to TOTAL degree 5
    #    (product=False): the probe must satisfy i+j+k <= deg. Probing it with a
    #    per-variable term like x^4 y^4 (total degree 8) would spuriously "fail" a
    #    rule that is genuinely degree-5 exact (verified by hand: 2B b^4 + 8C c^4 =
    #    8/5, 8C c^4 = 8/9, weight sum 6B + 8C = 8).
    for rule, deg, product, name in ((GP27, 5, True, "27-pt"),
                                     (GP8, 3, True, "8-pt"),
                                     (GP14, 5, False, "Irons 14-pt")):
        p = sum(sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
                * XI**i * ETA**j * ZETA**k
                for i in range(deg + 1) for j in range(deg + 1) for k in range(deg + 1)
                if (max(i, j, k) <= deg if product else i + j + k <= deg)
                and rng.random() < 0.15)
        # R2 hardening: guarantee a nonzero high-degree EVEN distinguisher
        # survives the random filter (else an unlucky seed could drop every
        # degree>=4 even term and a mere degree-3 rule would pass vacuously).
        # Top even order each rule must integrate exactly: per-variable for the
        # tensor-product rules (xi^te eta^te zeta^te), TOTAL for Irons (xi^te).
        te = 2 * (deg // 2)  # 4 for deg 5, 2 for deg 3
        p += (XI**te * ETA**te * ZETA**te) if product else XI**te
        exact = sp.integrate(p, (XI, -1, 1), (ETA, -1, 1), (ZETA, -1, 1))
        approx = sum(w * p.subs(dict(zip(_XYZ, (x, y, z)))) for x, y, z, w in rule)
        assert sp.simplify(exact - approx) == 0, f"{name} degree-{deg} exactness FAILED"
        assert sp.simplify(sum(r[3] for r in rule) - 8) == 0, f"{name} weight sum FAILED"
    print("G3 PASS  quadrature degrees (27-pt/8-pt per-variable deg 5/3, "
          "Irons 14-pt total deg-5) + weight sums")


def gate_G4_jacobian():
    for smp in XI_SAMPLES:
        _, dN = shape_numeric(smp)
        _, detJ = jacobian(X_DISTORT, dN)
        assert detJ > 0, f"detJ <= 0 at {smp}"
    for x, y, z, w in GP27_F:
        _, dN = shape_numeric((x, y, z))
        _, detJ = jacobian(X_DISTORT, dN)
        assert detJ > 0, "detJ <= 0 at a 27-pt GP"
    print("G4 PASS  distorted-hex Jacobian positive at samples + all GPs")


def _rank(K, tol_ratio=1e-9):
    s = np.linalg.svd(K, compute_uv=False)
    return int((s > tol_ratio * s[0]).sum()), s


def gate_G5_ranks(log="hex20_spurious_modes.txt"):
    D = iso_D()
    r27, _ = _rank(stiffness(X_CUBE, D, GP27_F))
    assert r27 == 54, f"27-pt rank {r27} != 54"
    K8 = stiffness(X_CUBE, D, GP8_F)
    r8, s8 = _rank(K8)
    assert r8 == 48, f"8-pt rank {r8} != 48 (deficiency != 6)"
    r14, _ = _rank(stiffness(X_CUBE, D, GP14_F))
    assert r14 == 54, f"Irons 14-pt rank {r14} != 54 (NOT rank-sufficient)"
    # catalogue: 12-dim nullspace of K8 = 6 RBM + 6 spurious; project out RBMs
    w, V = np.linalg.eigh(K8)
    null = V[:, w < 1e-9 * w[-1]]
    assert null.shape[1] == 12, f"K8 nullspace dim {null.shape[1]} != 12"
    X = X_CUBE
    rbm = np.zeros((60, 6))
    for a in range(20):
        rbm[3 * a:3 * a + 3, :3] = np.eye(3)
        x, y, z = X[a]
        rbm[3 * a:3 * a + 3, 3] = (0, -z, y)   # rot x
        rbm[3 * a:3 * a + 3, 4] = (z, 0, -x)   # rot y
        rbm[3 * a:3 * a + 3, 5] = (-y, x, 0)   # rot z
    Q, _ = np.linalg.qr(rbm)
    spurious = null - Q @ (Q.T @ null)
    U, sv, _ = np.linalg.svd(spurious, full_matrices=False)
    n_spur = int((sv > 1e-8).sum())
    assert n_spur == 6, f"spurious mode count {n_spur} != 6"
    with open(log, "w") as f:
        f.write("# 6 spurious zero-energy modes of H20 @ 2x2x2 (cube), ADR 72 P0\n")
        f.write("# rows = 60 dofs (node-major x,y,z); columns = modes 1..6\n")
        np.savetxt(f, U[:, :6], fmt="%+.6e")
    print(f"G5 PASS  ranks 54/48/54 (27/8/14-pt); 6 spurious modes -> {log}")


def gate_G6_G7_mass():
    M = consistent_mass(X_CUBE, 1.0, GP27_F)
    assert (np.diag(M) > 0).all(), "consistent-mass diagonal not positive"
    # exact fractions via sympy (rho = 1, V = 8, per-direction mass = 8)
    ints = [sp.integrate(n, (XI, -1, 1), (ETA, -1, 1), (ZETA, -1, 1)) for n in N_SYM]
    ints2 = [sp.integrate(n**2, (XI, -1, 1), (ETA, -1, 1), (ZETA, -1, 1)) for n in N_SYM]
    V = sp.Integer(8)
    S = sum(ints2)
    for a in range(20):
        corner = 0 not in NODE_XI[a]
        rs_frac = sp.nsimplify(ints[a] / V)
        hrz_frac = sp.nsimplify(ints2[a] / S)
        if corner:
            assert rs_frac == sp.Rational(-1, 8), f"row-sum corner {a}: {rs_frac} != -1/8"
            assert hrz_frac == sp.Rational(7, 248), f"HRZ corner {a}: {hrz_frac} != 7/248"
        else:
            assert rs_frac == sp.Rational(1, 6), f"row-sum mid-edge {a}: {rs_frac} != 1/6"
            assert hrz_frac == sp.Rational(2, 31), f"HRZ mid-edge {a}: {hrz_frac} != 2/31"
    # numeric HRZ path agrees with the exact fractions and conserves mass
    md = hrz_lump(M)
    assert np.allclose(md[0::3].sum(), 8.0, rtol=1e-12), "HRZ mass not conserved"
    assert (md > 0).all(), "HRZ produced a non-positive mass"
    assert np.allclose(md[0], 8 * 7 / 248, rtol=1e-10), "numeric HRZ corner != 7/248 * M"
    assert np.allclose(md[24], 8 * 2 / 31, rtol=1e-10), "numeric HRZ mid-edge != 2/31 * M"
    rs = M.sum(axis=1)
    assert rs[0] < 0, "row-sum corner mass unexpectedly positive"
    print("G6/G7 PASS  mass diag > 0; EXACT row-sum -1/8 & +1/6; HRZ 7/248 & 2/31 (positive, conserving)")


def main():
    gate_G1_basics()
    gate_G2_completeness()
    gate_G3_quadrature()
    gate_G4_jacobian()
    gate_G5_ranks()
    gate_G6_G7_mass()
    print("ALL P0 ORACLE GATES PASS (hex20_reference.py)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
