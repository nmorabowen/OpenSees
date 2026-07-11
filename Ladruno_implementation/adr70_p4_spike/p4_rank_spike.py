# ADR-70 P4 design spike — rank tests for triangle volumetric-cure candidates.
# Small-strain limit (F = I) rank analysis: the F-bar B-bar operator is
#   eps_bar = eps + ((theta_tilde - theta)/2) * m,  m = [1,1,0]^T (Voigt {xx,yy,xy})
# where theta = tr(eps) and theta_tilde is the PROJECTED dilatation.
# Zero-energy modes of any SPD-D stiffness = null space of the stacked
# sqrt(w)*B_bar rows.  RBM count in 2D = 3.
import numpy as np

np.set_printoptions(precision=3, suppress=True)

# ---------- T6 machinery ----------
# area coords L1 = 1-xi-eta, L2 = xi, L3 = eta
# nodes: 3 corners then midsides n4=(1-2), n5=(2-3), n6=(3-1)
def t6_shape_derivs(xi, eta):
    L1, L2, L3 = 1.0 - xi - eta, xi, eta
    # dN/d(xi,eta): dL1 = (-1,-1), dL2 = (1,0), dL3 = (0,1)
    dL = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    dN = np.zeros((6, 2))
    L = [L1, L2, L3]
    for i in range(3):
        dN[i] = (4.0 * L[i] - 1.0) * dL[i]
    dN[3] = 4.0 * (L2 * dL[0] + L1 * dL[1])   # 4*L1*L2
    dN[4] = 4.0 * (L3 * dL[1] + L2 * dL[2])   # 4*L2*L3
    dN[5] = 4.0 * (L1 * dL[2] + L3 * dL[0])   # 4*L3*L1
    return dN

GP3 = [(2/3, 1/6), (1/6, 2/3), (1/6, 1/6)]
W3 = [1/6, 1/6, 1/6]

def t6_nodes(corners):
    c = np.asarray(corners, dtype=float)
    mids = np.array([(c[0]+c[1])/2, (c[1]+c[2])/2, (c[2]+c[0])/2])
    return np.vstack([c, mids])

def element_B_rows(xy6):
    """per-GP: (B(3x12), b_dil(12,), w*detJ) for the 3-pt interior rule."""
    out = []
    for (xi, eta), w in zip(GP3, W3):
        dN = t6_shape_derivs(xi, eta)          # 6x2 in (xi,eta)
        J = xy6.T @ dN                          # 2x2  dx/dxi
        detJ = np.linalg.det(J)
        assert detJ > 0, "winding"
        dNx = dN @ np.linalg.inv(J)             # 6x2 in (x,y)
        B = np.zeros((3, 12))
        for a in range(6):
            B[0, 2*a]   = dNx[a, 0]
            B[1, 2*a+1] = dNx[a, 1]
            B[2, 2*a]   = dNx[a, 1]
            B[2, 2*a+1] = dNx[a, 0]
        b = B[0] + B[1]                         # dilatation row
        out.append((B, b, w * detJ))
    return out

M_VOIGT = np.array([1.0, 1.0, 0.0])

def stack_Bbar(elem_rows_list, dofmaps, ndof, projector):
    """elem_rows_list[e] = list of (B, b, wdet) per GP; dofmaps[e] maps local->global.
    projector(all_gp) -> b_tilde rows (global-dof length) per GP index.
    all_gp = list of (e, B, b, wdet) in order."""
    all_gp = []
    for e, rows in enumerate(elem_rows_list):
        for (B, b, wdet) in rows:
            Bg = np.zeros((3, ndof)); bg = np.zeros(ndof)
            for a in range(len(dofmaps[e]) // 2):
                ga = dofmaps[e][2*a] // 2
                Bg[:, 2*ga:2*ga+2] = B[:, 2*a:2*a+2]
                bg[2*ga:2*ga+2] = b[2*a:2*a+2]
            all_gp.append((e, Bg, bg, wdet))
    btil = projector(all_gp)
    rows = []
    for k, (e, Bg, bg, wdet) in enumerate(all_gp):
        Bbar = Bg + 0.5 * np.outer(M_VOIGT, btil[k] - bg)
        rows.append(np.sqrt(wdet) * Bbar)
    return np.vstack(rows)

def zero_modes(S, ndof, tol=1e-10):
    sv = np.linalg.svd(S, compute_uv=False)
    return ndof - int((sv > tol * sv[0]).sum())

# projectors -----------------------------------------------------------
def proj_none(all_gp):
    return [bg for (_, _, bg, _) in all_gp]

def proj_const_per_element(all_gp):
    btil = []
    for (e0, _, _, _) in all_gp:
        num = sum(wd * bg for (e, _, bg, wd) in all_gp if e == e0)
        den = sum(wd for (e, _, _, wd) in all_gp if e == e0)
        btil.append(num / den)
    return btil

def proj_const_patch(all_gp):
    num = sum(wd * bg for (_, _, bg, wd) in all_gp)
    den = sum(wd for (_, _, _, wd) in all_gp)
    return [num / den] * len(all_gp)

def make_proj_P1_patch(gp_xy):
    """L2 projection of the dilatation field onto GLOBAL P1 {1,x,y} over the patch,
    using the same quadrature. gp_xy[k] = physical (x,y) of GP k."""
    def proj(all_gp):
        P = np.array([[1.0, x, y] for (x, y) in gp_xy])          # nGP x 3
        Wd = np.diag([wd for (_, _, _, wd) in all_gp])
        Bmat = np.vstack([bg for (_, _, bg, _) in all_gp])       # nGP x ndof
        Mm = P.T @ Wd @ P
        coef = np.linalg.solve(Mm, P.T @ Wd @ Bmat)              # 3 x ndof
        Btil = P @ coef                                          # nGP x ndof
        return [Btil[k] for k in range(len(all_gp))]
    return proj

def gp_phys(xy6):
    pts = []
    for (xi, eta) in GP3:
        L1, L2, L3 = 1 - xi - eta, xi, eta
        N = np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                      4*L1*L2, 4*L2*L3, 4*L3*L1])
        pts.append(N @ xy6)
    return pts

# ---------- geometry: unit square split along the diagonal ----------
A, B_, C, D = (0, 0), (1, 0), (1, 1), (0, 1)
T1c, T2c = [A, B_, C], [A, C, D]
xy1, xy2 = t6_nodes(T1c), t6_nodes(T2c)

# global node table for the pair (9 unique nodes)
def build_pair():
    nodes = []
    def nid(p):
        for i, q in enumerate(nodes):
            if np.allclose(p, q): return i
        nodes.append(np.array(p, float)); return len(nodes) - 1
    m1 = [nid(p) for p in xy1]
    m2 = [nid(p) for p in xy2]
    return nodes, m1, m2

nodes, map1, map2 = build_pair()
ndof_pair = 2 * len(nodes)
dof1 = [2*n + d for n in map1 for d in (0, 1)]
dof2 = [2*n + d for n in map2 for d in (0, 1)]
rows1, rows2 = element_B_rows(xy1), element_B_rows(xy2)

print(f"pair: {len(nodes)} nodes, {ndof_pair} DOFs (expect 9 / 18)")

# --- Case 1: single T6 std -------------------------------------------
S = stack_Bbar([rows1], [list(range(12))], 12, proj_none)
print("C1  single T6, std                      : zero modes =", zero_modes(S, 12), "(expect 3)")

# --- Case 2: single T6, element-constant projection (P3 refutation) --
S = stack_Bbar([rows1], [list(range(12))], 12, proj_const_per_element)
print("C2  single T6, const-mean B-bar (P3)    : zero modes =", zero_modes(S, 12), "(expect 5 = P3 refutation)")

# --- Case 3: T6 pair, per-element constant projection ----------------
S = stack_Bbar([rows1, rows2], [dof1, dof2], ndof_pair, proj_const_per_element)
print("C3  T6 pair, per-elem const B-bar       : zero modes =", zero_modes(S, ndof_pair), "(predict 4)")

# --- Case 4: T6 pair, PATCH-constant projection (F-bar-Patch on T6) --
S = stack_Bbar([rows1, rows2], [dof1, dof2], ndof_pair, proj_const_patch)
print("C4  T6 pair, patch-const (F-bar-Patch)  : zero modes =", zero_modes(S, ndof_pair), "(claim A: predict 5 -> REFUTED)")

# --- Case 5: T6 pair, patch-P1 projected dilatation ------------------
gpxy = gp_phys(xy1) + gp_phys(xy2)
S = stack_Bbar([rows1, rows2], [dof1, dof2], ndof_pair, make_proj_P1_patch(gpxy))
print("C5  T6 pair, patch-P1 projection        : zero modes =", zero_modes(S, ndof_pair), "(claim D: predict 3 -> rank OK)")

# --- Case 6: single T6, element-local P1disc projection with 3-pt rule
prj = make_proj_P1_patch(gp_phys(xy1))
all_gp = []
for (Bm, b, wdet) in rows1:
    all_gp.append((0, None, b, wdet))
btil = prj([(e, None, bg, wd) for (e, _, bg, wd) in all_gp])
diff = max(np.max(np.abs(btil[k] - all_gp[k][2])) for k in range(3))
print(f"C6  single T6, elem-local P1disc @3GP   : max|b_tilde - b| = {diff:.2e} (claim B: identity -> no-op)")

# --- Case 7: T3 pair, dSNPO 15.1.9 patch-constant --------------------
def t3_rows(corners):
    c = np.asarray(corners, float)
    area = 0.5 * np.linalg.det(np.array([c[1]-c[0], c[2]-c[0]]))
    assert area > 0
    Bm = np.zeros((3, 6))
    for a in range(3):
        j, k = (a+1) % 3, (a+2) % 3
        bx = (c[j][1] - c[k][1]) / (2*area)
        by = (c[k][0] - c[j][0]) / (2*area)
        Bm[0, 2*a] = bx; Bm[1, 2*a+1] = by
        Bm[2, 2*a] = by; Bm[2, 2*a+1] = bx
    return [(Bm, Bm[0] + Bm[1], area)]

t3nodes = [np.array(p, float) for p in (A, B_, C, D)]
t3map1, t3map2 = [0, 1, 2], [0, 2, 3]
t3dof1 = [2*n+d for n in t3map1 for d in (0, 1)]
t3dof2 = [2*n+d for n in t3map2 for d in (0, 1)]
r31, r32 = t3_rows([A, B_, C]), t3_rows([A, C, D])
S = stack_Bbar([r31, r32], [t3dof1, t3dof2], 8, proj_const_patch)
print("C7  T3 pair, 15.1.9 patch-const         : zero modes =", zero_modes(S, 8), "(predict 3 -> OK)")
S = stack_Bbar([r31, r32], [t3dof1, t3dof2], 8, proj_none)
print("C7b T3 pair, std (reference)            : zero modes =", zero_modes(S, 8), "(expect 3)")

# --- Case 8: constraint counting sanity (volumetric rows rank) -------
# count independent volumetric constraints seen by each formulation on the T6 pair
def vol_rank(projector):
    all_gp = []
    for e, rows, dm in ((0, rows1, dof1), (1, rows2, dof2)):
        for (Bm, b, wdet) in rows:
            bg = np.zeros(ndof_pair)
            for a in range(6):
                ga = dm[2*a] // 2
                bg[2*ga:2*ga+2] = b[2*a:2*a+2]
            all_gp.append((e, None, bg, wdet))
    btil = projector(all_gp)
    Vm = np.vstack([np.sqrt(all_gp[k][3]) * btil[k] for k in range(len(all_gp))])
    sv = np.linalg.svd(Vm, compute_uv=False)
    return int((sv > 1e-10 * sv[0]).sum())

print("C8  T6-pair volumetric constraint count : std =", vol_rank(proj_none),
      "| patch-const =", vol_rank(proj_const_patch),
      "| patch-P1 =", vol_rank(make_proj_P1_patch(gpxy)))
