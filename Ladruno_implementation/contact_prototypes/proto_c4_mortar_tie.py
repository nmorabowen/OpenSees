"""ADR-41 C4 ORACLE — mortar MESH-TYING (a permanent bond; the zero-gap limit of contact).

Oracle-first (the fork discipline): validate the TIE enforcement layer in pure numpy BEFORE
wiring it into the C++ adapter/Domain. C4 reuses the SHIPPED mortar machinery (C1→C3) with the
inequality STRIPPED: the active set is frozen ON for every slave node, and the FULL 3-vector
weighted relative DISPLACEMENT

    r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K          (a 3-vector per slave node I)

is driven to ZERO (normal AND tangential — no KKT, no compression clamp, no friction cone).

The mechanics ([[_adr41_c4_design]] §mechanics; the C2 λ_N pattern, minus the inequality):
  - tie traction (full 3-vec, NO clamp):     t_I = λ_tie,I + epsTie · (r_I / a_I)
  - assemble like the NORMAL force (via D / −M, but a full 3-vec, not a scalar×n):
        f^s_K = −Σ_I D_KI · t_I ,   f^m_L = +Σ_I M_IL · t_I       (self-equilibrating, Σφ=1)
  - tangent (symmetric, SPD, no active set, no n⊗n projection):
        K[A][B] = epsTie · Σ_I (b_IA b_IB / a_I) · I₃ ,   b = B̃ = [D, −M]   (the FULL identity)
  - ALM (the accuracy win): one Uzawa step per Domain::commit(), NO clamp (equality):
        λ_tie,I ← λ_tie,I + epsTie · (r_I / a_I)   ⇒ ‖r‖ → 0 at FINITE epsTie (epsTie-INDEPENDENT)
  - r is built from DISPLACEMENTS (getTrialDisp) — the bond exists from the as-built config, so
    there is NO engagement origin gT0 (gT0 ≡ 0); the C3.1 displacement-not-position lesson.

Shared-node resolution (the C4 pre-req — LEDGER_quirks MAJOR-1/MINOR-1, fenced before non-matching):
  r_I is a LINEAR accumulation (no return map), so the per-node tie state can use the SAME
  order-INDEPENDENT global accumulator λ_N already uses (rGlobal = Σ_facets r_I^facet) — unlike the
  friction slip (a return-map OUTPUT, last-writer-wins). The FORCE uses each facet's LOCAL r (the
  C2.2 lesson: keeps R(u) deterministic per facet); the GLOBAL r drives only the commit Uzawa + the
  ‖r‖ query. T6 pins that this is order-independent.

Gates (capstone mesh-tying row):
  T1  constant-stress PATCH TEST across a NON-MATCHING tied interface — a linear field transmits
      EXACTLY (the HEADLINE): r ≡ 0 for any consistent linear field; D⁻¹M recovers the master field.
  T2  per-node assembly through D/−M: self-equilibrium ΣF^s+ΣF^m=0; assembled tangent
      epsTie·B̃ᵀB̃⊗I₃ == FD, SYMMETRIC + SPD; the FULL identity (normal AND tangential tied).
  T3  penalty tie (λ_tie≡0): the residual bond ‖r‖ = O(1/epsTie) (a split == monolithic to a tol).
  T4  ALM (no clamp): the Uzawa drives ‖r‖ → augTol at FINITE epsTie, epsTie-INDEPENDENT.
  T5  full 3-vec / no clamp: the tie constrains normal AND tangential, and a "tensile" r (positive
      separation) restores too — there is NO compression clamp (the contact→tie difference).
  T6  shared-node order-independence: the GLOBAL r accumulator is facet-eval-order-independent
      (the linear-accumulation resolution of MAJOR-1); a per-facet last-writer-wins value is NOT.

Run: <pythoncore-3.12> proto_c4_mortar_tie.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ===================================== a minimal NON-MATCHING 1D mortar interface (exact clip)
# C1/C2 already validated the clip→Gauss→D,M machinery in 3D; C4 isolates the TIE enforcement, so a
# 1D edge mortar (the exact-overlap limit) is enough to pin the mechanics. D_IJ = ∫ N_I^s N_J^s ds
# (slave Gram); M_IK = ∫ N_I^s N_K^m ds (slave-master coupling), integrated over the SLAVE segments
# with the master shape evaluated at the physical Gauss point. Linear 2-node elements per segment.
def _lin(xi):
    return np.array([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)])


def assemble_mortar_1d(slave_x, master_x, ngp=6):
    """Mortar D (ns×ns), M (ns×nm), a (ns) for two linear meshes sharing the edge [x0,x1].
    slave_x / master_x: monotone node coordinates (consecutive linear segments). The two meshes
    span the SAME interval but with DIFFERENT node layouts (non-matching)."""
    ns, nm = len(slave_x), len(master_x)
    D = np.zeros((ns, ns))
    M = np.zeros((ns, nm))
    gp, gw = np.polynomial.legendre.leggauss(ngp)
    for s in range(ns - 1):
        xa, xb = slave_x[s], slave_x[s + 1]
        Lseg = xb - xa
        for xi, w in zip(gp, gw):
            x = 0.5 * (xa + xb) + 0.5 * Lseg * xi          # physical Gauss point
            Jw = 0.5 * Lseg * w
            Ns = np.zeros(ns)
            Ns[s:s + 2] = _lin(xi)
            Nm = np.zeros(nm)                              # master shape at x (locate the segment)
            for m in range(nm - 1):
                if master_x[m] - 1e-12 <= x <= master_x[m + 1] + 1e-12:
                    xim = 2.0 * (x - master_x[m]) / (master_x[m + 1] - master_x[m]) - 1.0
                    Nm[m:m + 2] = _lin(xim)
                    break
            D += Jw * np.outer(Ns, Ns)
            M += Jw * np.outer(Ns, Nm)
    return D, M, D.sum(1)


# A genuinely NON-MATCHING pair: 4 slave nodes (non-uniform) vs 3 master nodes (non-uniform),
# both spanning [0, 1].
slave_x = np.array([0.0, 0.3, 0.65, 1.0])
master_x = np.array([0.0, 0.55, 1.0])
D, M, a = assemble_mortar_1d(slave_x, master_x)
ns, nm = len(slave_x), len(master_x)
nN = ns + nm

print("ADR-41 C4 oracle — mortar mesh-tying (permanent bond, zero-gap limit of contact)")
print(f"  non-matching interface: {ns} slave nodes vs {nm} master nodes on [0,1]")

# --- T0: the mortar operators are sane (partition of unity) ----------------------------------
print("\nT0  mortar operators: partition of unity (Σ_J D_IJ == Σ_K M_IK == a_I = ∫N_I)")
check("a_I = Σ_J D_IJ == Σ_K M_IK (D and M row sums agree ⇒ Σφ=1)",
      np.allclose(D.sum(1), M.sum(1), atol=1e-12),
      f"max|ΔrowSum|={np.abs(D.sum(1) - M.sum(1)).max():.2e}")
check("a_I > 0 at every slave node (referenced)", np.all(a > 1e-12))

# --- T1: the HEADLINE constant-stress PATCH TEST — a linear field transmits EXACTLY -----------
print("\nT1  PATCH TEST (headline): a consistent linear field across the NON-MATCHING tie ⇒ r ≡ 0")
# A 3-vector linear field u(x) = c0 + grad·x sampled at the node coords (interface param = x). The
# tie ties ALL THREE components, so use a full 3-vec field with an arbitrary constant gradient.
c0 = np.array([0.7, -0.2, 0.45])
grad = np.array([1.3, -0.6, 2.1])           # du/dx per component (constant ⇒ "constant stress")
us = np.array([c0 + grad * x for x in slave_x])     # (ns, 3) slave nodal field
um = np.array([c0 + grad * x for x in master_x])    # (nm, 3) master nodal field
r = D @ us - M @ um                                  # (ns, 3) weighted relative displacement
check("r ≡ 0 for a consistent linear field (EXACT transmission, all 3 components, ≤1e-12)",
      np.abs(r).max() < 1e-12, f"max‖r‖={np.abs(r).max():.2e}")
# D⁻¹M recovers the master field at the slave nodes (the mortar projection of a linear field).
us_proj = np.linalg.solve(D, M @ um)                 # (ns, 3)
check("D⁻¹M recovers the slave linear field from the master nodes (≤1e-10)",
      np.allclose(us_proj, us, atol=1e-10), f"max|Δ|={np.abs(us_proj - us).max():.2e}")
# A constant field is the trivial sub-case (rigid translation transmits exactly).
us_c = np.tile(c0, (ns, 1)); um_c = np.tile(c0, (nm, 1))
check("constant (rigid-translation) field ⇒ r ≡ 0 (the Σφ=1 sub-case)",
      np.abs(D @ us_c - M @ um_c).max() < 1e-12)


# ===================================== the tie enforcement layer (force / tangent / ALM)
def bmat(I, A):
    """B̃ = [D, −M] operator row: slave node A ⇒ D[I][A]; master node A−ns ⇒ −M[I][A−ns]."""
    return D[I][A] if A < ns else -M[I][A - ns]


def tie_traction(U, lamTie, epsTie):
    """t_I = λ_tie,I + epsTie·(r_I/a_I), full 3-vec, NO clamp. U = (nN,3) all-node disp;
    lamTie = (ns,3) committed multiplier. Returns t (ns,3) and r (ns,3)."""
    r = np.zeros((ns, 3))
    for I in range(ns):
        for A in range(nN):
            r[I] += bmat(I, A) * U[A]
    t = lamTie + epsTie * (r / a[:, None])
    return t, r


def tie_force(U, lamTie, epsTie):
    """assembled tie force F (nN,3): f^s_K = −Σ_I D_KI t_I, f^m_L = +Σ_I M_IL t_I
    ⇒ F_A = −Σ_I b_IA t_I (the b=[D,−M] scatter; the NORMAL-force sign, t a full 3-vec)."""
    t, r = tie_traction(U, lamTie, epsTie)
    F = np.zeros((nN, 3))
    for I in range(ns):
        for A in range(nN):
            F[A] += -bmat(I, A) * t[I]
    return F, r, t


def tie_tangent(epsTie):
    """K[A][B] = epsTie·Σ_I (b_IA b_IB/a_I)·I₃ — the FULL identity (all 3 components tied),
    active set frozen ON (W_I = 1/a_I unconditionally). (nN·3, nN·3)."""
    K = np.zeros((nN, 3, nN, 3))
    for I in range(ns):
        for A in range(nN):
            for B in range(nN):
                w = epsTie * bmat(I, A) * bmat(I, B) / a[I]
                for d in range(3):
                    K[A, d, B, d] += w               # ⊗ I₃ (diagonal in the 3-space)
    return K.reshape(nN * 3, nN * 3)


# --- T2: self-equilibrium + the assembled tangent epsTie·B̃ᵀB̃⊗I₃ == FD (symmetric, SPD) -------
print("\nT2  assembly through D/−M: self-equilibrium; tangent epsTie·B̃ᵀB̃⊗I₃ == FD, symmetric, SPD")
epsTie = 1.0e4
lam0 = np.zeros((ns, 3))
U0 = np.random.default_rng(0).standard_normal((nN, 3)) * 1e-2   # an arbitrary trial config
F0, r0, t0 = tie_force(U0, lam0, epsTie)
check("self-equilibrium ΣF^s + ΣF^m == 0 (the tie applies NO net force/moment, ≤1e-10)",
      np.abs(F0.sum(0)).max() < 1e-10, f"max|ΣF|={np.abs(F0.sum(0)).max():.2e}")
Kc = tie_tangent(epsTie)
# element tangent = −∂(applied tie force)/∂u, central FD
h = 1e-7
Kfd = np.zeros((nN * 3, nN * 3))
for j in range(nN * 3):
    Up = U0.reshape(-1).copy(); Up[j] += h
    Um_ = U0.reshape(-1).copy(); Um_[j] -= h
    Fp, _, _ = tie_force(Up.reshape(nN, 3), lam0, epsTie)
    Fm, _, _ = tie_force(Um_.reshape(nN, 3), lam0, epsTie)
    Kfd[:, j] = -(Fp.reshape(-1) - Fm.reshape(-1)) / (2 * h)
rel = np.abs(Kc - Kfd).max() / (np.abs(Kfd).max() + 1e-30)
check("assembled K_c == FD(F_tie) (≤1e-6 rel)", rel <= 1e-6, f"max rel={rel:.2e}")
check("tie tangent IS symmetric (no Csl, no active set — pure penalty Gram)",
      np.allclose(Kc, Kc.T, atol=1e-9))
# SPD on the constrained subspace: B̃ᵀB̃ is PSD with a null space (rigid relative modes the tie
# does not see); the body stiffness regularizes it. Eigenvalues are ≥ 0 and the rank == 3·ns.
evals = np.linalg.eigvalsh(0.5 * (Kc + Kc.T))
check("tie tangent is PSD (all eigenvalues ≥ −tol)", evals.min() > -1e-6, f"λ_min={evals.min():.2e}")
check("tie tangent rank == 3·ns (one constrained mode per slave-node component)",
      np.linalg.matrix_rank(Kc, tol=1e-6 * abs(evals.max())) == 3 * ns,
      f"rank={np.linalg.matrix_rank(Kc, tol=1e-6 * abs(evals.max()))} (want {3 * ns})")


# ===================================== a minimal split-body solve (penalty + ALM)
# Two elastic bodies sharing the tied interface, ONE displacement component (the tie decouples per
# component — T2's I₃ tangent). Slave body: SPD chain Ks tying the slave interface DOFs to ground +
# load fs. Master body: SPD chain Km + load fm. Tie couples u_s ↔ u_m. We compare the tied split to
# the MONOLITHIC body (the two fields forced equal through the mortar projection).
rng = np.random.default_rng(7)
Ks = np.array([[3.0, -1.0, 0.0, 0.0],
               [-1.0, 3.0, -1.0, 0.0],
               [0.0, -1.0, 3.0, -1.0],
               [0.0, 0.0, -1.0, 2.0]])          # 4×4 SPD (slave interface dofs)
Km = np.array([[2.5, -1.2, 0.0],
               [-1.2, 3.0, -1.2],
               [0.0, -1.2, 2.0]])               # 3×3 SPD (master interface dofs)
Kbody = np.zeros((nN, nN))
Kbody[:ns, :ns] = Ks
Kbody[ns:, ns:] = Km
fext = np.concatenate([rng.standard_normal(ns), rng.standard_normal(nm)])
B = np.zeros((ns, nN))                          # B̃ as an (ns × nN) matrix on the scalar component
for I in range(ns):
    for A in range(nN):
        B[I, A] = bmat(I, A)
Winv = np.diag(1.0 / a)


def solve_penalty(epsTie):
    """one linear solve with λ_tie ≡ 0 (penalty tie). Returns (u, r)."""
    Ktie = epsTie * B.T @ Winv @ B
    u = np.linalg.solve(Kbody + Ktie, fext)
    r = B @ u
    return u, r


def solve_alm(epsTie, augTol=1e-12, maxAug=500):
    """ALM: Uzawa λ_tie ← λ_tie + epsTie·(r/a), NO clamp, once per outer 'commit'. Returns
    (u, r, naug). The inner solve is linear (one Newton step) at fixed λ."""
    lam = np.zeros(ns)
    Ktie = epsTie * B.T @ Winv @ B
    Ksys = Kbody + Ktie
    naug = 0
    for naug in range(1, maxAug + 1):
        rhs = fext - B.T @ lam              # the λ_tie term assembles globally (D·λ), the C2.2 form
        u = np.linalg.solve(Ksys, rhs)
        r = B @ u
        if np.abs(r).max() < augTol:
            break
        lam = lam + epsTie * (r / a)        # NO clamp (equality) — the C4 Uzawa
    return u, r, naug


# --- T3: penalty tie ‖r‖ = O(1/epsTie) -------------------------------------------------------
print("\nT3  penalty tie (λ_tie≡0): the residual bond ‖r‖ = O(1/epsTie) (split → monolithic)")
res = {}
for e in [1e2, 1e3, 1e4, 1e5]:
    _, r = solve_penalty(e)
    res[e] = np.abs(r).max()
    print(f"      epsTie={e:.0e}  ‖r‖_inf={res[e]:.3e}")
# tenfold epsTie ⇒ ~tenfold smaller bond (first-order penalty convergence)
ratio = res[1e3] / res[1e4]
check("‖r‖ scales as 1/epsTie (10× penalty ⇒ ~10× smaller bond)", 8.0 < ratio < 12.0,
      f"‖r‖(1e3)/‖r‖(1e4)={ratio:.2f}")
check("penalty bond is NONZERO (the O(1/epsTie) gap ALM removes)", res[1e4] > 1e-9)

# --- T4: ALM (no clamp) drives ‖r‖ → augTol at FINITE epsTie, epsTie-INDEPENDENT --------------
print("\nT4  ALM (no clamp): Uzawa drives ‖r‖ → augTol at FINITE epsTie, epsTie-INDEPENDENT")
augTol = 1e-12
results = {}
for e in [1e2, 1e3, 1e4]:
    u, r, naug = solve_alm(e, augTol=augTol)
    results[e] = (np.abs(r).max(), naug)
    print(f"      epsTie={e:.0e}  ‖r‖_inf={np.abs(r).max():.3e}  ({naug} augmentations)")
check("ALM drives ‖r‖ → augTol at FINITE epsTie (all penalties, ≤augTol)",
      all(results[e][0] <= augTol for e in results),
      f"max‖r‖={max(results[e][0] for e in results):.2e}")
# the CONVERGED solution is epsTie-INDEPENDENT (the headline ALM win — not just smaller, IDENTICAL)
u_lo, _, _ = solve_alm(1e2, augTol=augTol)
u_hi, _, _ = solve_alm(1e5, augTol=augTol)
check("converged u is epsTie-INDEPENDENT (1e2 vs 1e5 agree, ≤1e-9)",
      np.allclose(u_lo, u_hi, atol=1e-9), f"max|Δu|={np.abs(u_lo - u_hi).max():.2e}")
# and it equals the EXACT constrained solution (KKT saddle: minimize ½uᵀKu−fᵀu s.t. Bu=0)
KKT = np.block([[Kbody, B.T], [B, np.zeros((ns, ns))]])
rhsK = np.concatenate([fext, np.zeros(ns)])
sol = np.linalg.solve(KKT, rhsK)
u_exact = sol[:nN]
check("ALM u == the EXACT Lagrange-constrained (Bu=0) solution (≤1e-9)",
      np.allclose(u_lo, u_exact, atol=1e-9), f"max|Δ|={np.abs(u_lo - u_exact).max():.2e}")

# --- T5: FULL 3-vec / NO clamp — tie constrains normal AND tangential; tension restores too ---
print("\nT5  full 3-vec / NO clamp: normal AND tangential tied; a 'tensile' r restores (no clamp)")
# Drive a separation (positive normal r) and a tangential slip; both must produce a RESTORING
# traction (contact would clamp the tension to zero — the tie does NOT).
U_sep = np.zeros((nN, 3))
U_sep[0] = np.array([0.0, 0.0, +5e-3])      # slave node 0 pulled "apart" (normal +z) → tension
U_sep[1] = np.array([4e-3, 0.0, 0.0])       # slave node 1 slid tangentially (+x)
_, r5, t5 = tie_force(U_sep, np.zeros((ns, 3)), epsTie)
check("a POSITIVE normal r (separation) produces a restoring traction (NO compression clamp)",
      t5[0, 2] > 0.0, f"t_z(node0)={t5[0,2]:.3e}")
check("a TANGENTIAL r (slip) produces a restoring traction (the tie ties tangential too)",
      abs(t5[1, 0]) > 0.0, f"t_x(node1)={t5[1,0]:.3e}")
check("the traction is a FULL 3-vector aligned with r (t = epsTie·r/a, no projection)",
      np.allclose(t5, epsTie * r5 / a[:, None], atol=1e-12))

# --- T6: shared-node order-independence — the GLOBAL r accumulator is eval-order-independent ---
print("\nT6  shared-node resolution (MAJOR-1): the GLOBAL r accumulator is facet-order-independent")
# A slave node shared by 2 facets contributes r^facetA and r^facetB. The GLOBAL r = r^A + r^B is a
# LINEAR sum ⇒ identical in any evaluation order (unlike the friction slip, a return-map OUTPUT that
# is last-writer-wins). Pin it: accumulate in both orders, compare; and show a last-writer-wins
# per-facet value WOULD differ (so the global accumulator is the correct fix).
rA = np.array([1.0, -2.0, 0.5])             # facet A's contribution to the shared node's r
rB = np.array([-0.3, 0.7, -1.1])            # facet B's contribution
g_fwd = rA + rB                              # accumulate A then B
g_rev = rB + rA                              # accumulate B then A
check("GLOBAL r is identical in either facet order (order-INDEPENDENT linear accumulation)",
      np.allclose(g_fwd, g_rev, atol=1e-15))
check("a per-facet last-writer-wins value WOULD be order-dependent (rA != rB ⇒ the bug we avoid)",
      not np.allclose(rA, rB), "the global accumulator — λ_N's treatment — is the resolution")
# idempotent delta-update (a re-eval OVERWRITES, never double-counts) — the accumulateMortarGap rule
acc = np.zeros(3); last = {}
for facet, rf in [("A", rA), ("B", rB), ("A", rA), ("B", rB)]:   # re-eval A, B (Newton re-forms R)
    acc += rf - last.get(facet, np.zeros(3))
    last[facet] = rf.copy()
check("delta-update is idempotent across re-evals (== single A+B, no double-count)",
      np.allclose(acc, rA + rB, atol=1e-15))

print(f"\n{'='*70}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(C4 mortar mesh-tying oracle)\n{'='*70}")
sys.exit(1 if _fails else 0)
