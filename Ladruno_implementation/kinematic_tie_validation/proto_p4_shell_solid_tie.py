"""ADR-64 (LadrunoTie P4) ORACLE — shell-to-solid plane-section rigid-arm tie.

Ties an ndf-6 shell EDGE (MASTER) to an ndf-3 solid FACE (SLAVES) — direction b-B,
signed off (OQ-1): the solid face nodes are the SLAVES, the shell edge nodes the
MASTERS, and the tie is the linearized rigid plane-section map (Abaqus shell-to-solid
coupling, TG §6.6.2/§6.6.6, kinematic/RBE2-flavored form). Per solid face node s,
projected onto the master edge polyline at ξ (closest point, per-node FROZEN arm
d = x_s − x̄(ξ)), THREE translation rows:

    u_{s,i} = Σ_{j=A,B} N_j(ξ) · [ u_{j,i} + (θ_j × d)_i ]        i = 1..3

i.e. per master node j: N_j on the translation columns and N_j·(−[d]×) on the
rotation columns — the exact nk×6 mixed-triple shape ltEmitMixedRow already consumes
(LadrunoTie.cpp:154), including its 1e-12 near-zero coefficient drop, which is
LOAD-BEARING here: the arm d runs ~along the shell normal n, so the θ·n (drilling)
column of θ×d is ~zero and gets DROPPED ⇒ shell drilling is auto-free (the Abaqus
behavior), no rigid-body failure, no local-frame gymnastics.

v1 CONTRACT, signed off (OQ-2): statically-exact couple with a FROZEN small-rotation
arm (same status as rigidLink -beam; the handler's 0.1-rad flagRotMonitor auto-arms),
and through-thickness Poisson stretch SUPPRESSED at the seam (RBE2-style St-Venant
boundary layer, error ∝ ν·t·strain, exact at ν=0). Both honest limits are GATED
below (T6 Poisson, T7 Timoshenko shear), not engineered away. Per BLOCKER-4: if
either gate shows anything other than the predicted scaling, STOP — the operator is
wrong, not the tolerance.

Gates (falsifiable, numeric):
  T1  Rigid body, 6 modes, incl. drilling-falls-out-free: the θ·n column is DROPPED
      by the 1e-12 filter (not clamped — predictions invariant under arbitrary
      master drilling), all 6 modes strain-free; on-edge slave (|d|≈0) degenerates
      to identity translation rows (θ columns vanish naturally).
  T2  Constant-stress membrane patch (ν=0), non-conforming through-thickness slave
      layout: any in-plane constant strain crosses exactly.
  T3  HEADLINE: constant-MOMENT patch between dissimilar element types — shell side
      θ=κx, w=½κx² ↔ solid side u_a=−κxz linear-in-thickness, ν=0: EXACT by
      construction; the gate proves it.
  T4  Work-conjugacy / static equivalence: Cᵀλ reproduces the exact edge resultants
      (linear momentum by PoU, angular momentum incl. the transferred couples); the
      θ columns deliver a PURE COUPLE (zero net force); reduced stiffness TᵀKT
      symmetric to 1e-14.
  T5  3D orientation invariance: a randomly rotated frame (skew edge, skew arm)
      reproduces rigid + bending exactly; drilling stays free via the NULL SPACE of
      −[d]× (the filter can't drop a skew column — freedom is structural); a kinked
      polyline with per-node frozen arms still crosses rigid modes exactly (OQ-7).
  T6  Honest limit, Poisson (OQ-2): thickness-stretch suppression misfit = ν·ε·|z|,
      linear in ν AND in t, EXACT at ν=0; in-plane components stay exact.
  T7  Honest limit, Timoshenko shear (OQ-2): the parabolic-shear warping f(z) =
      (γ/2)(z − z³/c²) is unrepresentable by a plane section ⇒ misfit O(γ·t)
      (max = γc/(3√3)), linear in γ, linear in t, NOT shrinking with in-plane
      refinement — while the resultant transfer stays exact.
  T8  Contrast: θ columns zeroed (= the translations-only -dof 3 1 2 3 tie available
      today) still passes T2 but FAILS T3 (hinge) — the θ columns are the
      load-bearing delta of this rung.

Run: <pythoncore-3.12> proto_p4_shell_solid_tie.py   (numpy only, no build)
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


# ============================================================ the b-B operator
def skew(v):
    """[v]× : skew(v) @ w = v × w."""
    return np.array([[0.0, -v[2], v[1]],
                     [v[2], 0.0, -v[0]],
                     [-v[1], v[0], 0.0]])


def lt_project_edge(edge, xs, tol=1e-9):
    """Closest-point projection of xs onto the master edge polyline (mirrors the
    planned ltProjectEdge): closed-form per segment, brute force over segments,
    in-bounds (ξ ∈ [0,1] within tol) required — out-of-polyline is a REFUSAL.
    Returns (seg, ξ, footpoint)."""
    best = None
    for k in range(len(edge) - 1):
        a, b = edge[k], edge[k + 1]
        ab = b - a
        L2 = ab @ ab
        if L2 <= 0.0:
            continue
        xi = ((xs - a) @ ab) / L2
        if xi < -tol or xi > 1.0 + tol:
            continue
        xic = min(max(xi, 0.0), 1.0)
        foot = a + xic * ab
        dd = np.linalg.norm(xs - foot)
        if best is None or dd < best[0]:
            best = (dd, k, xic, foot)
    if best is None:
        raise ValueError("slave does not project in-bounds onto the master edge polyline")
    return best[1], best[2], best[3]


DROP = 1e-12          # ltEmitMixedRow's near-zero coefficient filter (LadrunoTie.cpp:161)


def plane_section_rows(edge, xs, zero_theta=False):
    """The 3 translation rows for one solid slave node, built exactly as
    generateShellSolid will: N_j on translations, N_j·(−[d]×) on rotations, frozen
    arm d = xs − footpoint, 1e-12 filter applied. Returns (C, d, dropped_mask):
    C is 3 × (n_m·6) over the stacked master (u,θ) DOFs.
    zero_theta=True = the T8 contrast (translations-only -dof 3 1 2 3 tie)."""
    seg, xi, foot = lt_project_edge(edge, xs)
    d = xs - foot
    B = -skew(d)                                   # ∂(θ×d)/∂θ
    C = np.zeros((3, len(edge) * 6))
    for j, Nj in ((seg, 1.0 - xi), (seg + 1, xi)):
        if abs(Nj) <= DROP:
            continue
        C[:, j * 6:j * 6 + 3] += Nj * np.eye(3)
        if not zero_theta:
            C[:, j * 6 + 3:j * 6 + 6] += Nj * B
    dropped = np.abs(C) <= DROP                    # mirror the emission-time drop
    C = np.where(dropped, 0.0, C)
    return C, d, dropped


def predict(rows, U):
    """Stack per-slave predictions u_s = C @ U (U = master (u,θ) DOFs, flat)."""
    return np.array([C @ U for C, _, _ in rows])


def master_state(edge, ufun, tfun):
    """Flat master DOF vector: per edge node, 3 translations + 3 rotations."""
    U = np.zeros(len(edge) * 6)
    for j, x in enumerate(edge):
        U[j * 6:j * 6 + 3] = ufun(x)
        U[j * 6 + 3:j * 6 + 6] = tfun(x)
    return U


# ============================================================ geometry
# Shell mid-surface in the xy-plane (normal n = ẑ), edge line along y at x = X0;
# solid face = the plane x = X0 with nodes spread through the shell thickness
# (z ∈ [−c, c], t = 2c). Master polyline deliberately NON-UNIFORM; slave (y, z)
# layout deliberately NON-CONFORMING (no slave y matches a master y).
X0 = 0.8
T_SHELL = 0.1
C_HALF = T_SHELL / 2.0

MY = np.array([0.0, 0.6, 1.0])                     # master edge y (2 segments, non-uniform)
EDGE = [np.array([X0, y, 0.0]) for y in MY]

SY = np.array([0.15, 0.3, 0.45, 0.7, 0.9])         # slave y (non-conforming)
SZ = np.array([-C_HALF, -0.4 * C_HALF, 0.2 * C_HALF, C_HALF])   # through-thickness
SLAVES = [np.array([X0, y, z]) for y in SY for z in SZ]
SLAVES.append(np.array([X0, 0.3, 0.0]))            # on-edge slave: |d| = 0 degeneracy (ξ=0.5 exact)

ROWS = [plane_section_rows(EDGE, xs) for xs in SLAVES]
XS = np.array(SLAVES)

print("ADR-64 P4 oracle — shell-to-solid plane-section rigid-arm tie (b-B, signed off)")
print(f"\n     edge x={X0}, master y={MY}, {len(SLAVES)} solid slaves over z∈[±{C_HALF}]"
      f" (t={T_SHELL}), non-conforming")

# --- T1: rigid body (6 modes) + drilling falls out FREE + on-edge degeneracy ------
print("\nT1  rigid body: 6 modes strain-free, drilling DROPPED by the filter (not clamped)")
rng = np.random.default_rng(64)
P0 = np.array([0.3, -0.2, 0.05])                   # rotation reference point
err = 0.0
for mode in range(6):
    u0 = np.zeros(3); om = np.zeros(3)
    if mode < 3:
        u0[mode] = 0.017
    else:
        om[mode - 3] = 0.013
    U = master_state(EDGE, lambda x: u0 + np.cross(om, x - P0), lambda x: om)
    exact = u0 + np.cross(om, XS - P0)
    err = max(err, np.abs(predict(ROWS, U) - exact).max())
check("all 6 rigid modes cross exactly", err < 1e-15, f"err={err:.2e}")

# drilling: the arm d = (0,0,z) ∥ shell normal ⇒ the θ_z column of −[d]× is exactly
# zero and must be DROPPED (this is what keeps the shell's drilling DOF out of the
# handler's union-find groups — free, not clamped).
drill_dropped = all(dropped[:, 3 + 2::6].all() for _, d, dropped in ROWS)
check("θ·n (drilling) columns all ≤1e-12 and dropped, every slave", drill_dropped)
Urand = rng.normal(size=len(EDGE) * 6) * 0.01
Udrill = Urand.copy()
Udrill[5::6] += rng.normal(size=len(EDGE)) * 10.0   # huge arbitrary drilling
e1 = np.abs(predict(ROWS, Urand) - predict(ROWS, Udrill)).max()
check("predictions INVARIANT under arbitrary master drilling (free, not clamped)",
      e1 == 0.0, f"delta={e1:.2e}")
# on-edge slave (|d| = 0): identity translation rows, θ columns vanish naturally
Con, don, dron = ROWS[-1]
theta_all_dropped = dron[:, [j * 6 + k for j in range(len(EDGE)) for k in (3, 4, 5)]].all()
row_sum = np.abs(Con[:, 0::6].sum(axis=1) + Con[:, 1::6].sum(axis=1) +
                 Con[:, 2::6].sum(axis=1) - 1.0).max()   # each row's N-weights sum to 1
check("on-edge slave (|d|≈0): θ columns vanish, identity-PoU translation rows",
      np.linalg.norm(don) == 0.0 and theta_all_dropped and row_sum < 1e-15,
      f"|d|={np.linalg.norm(don):.1e}, PoU err={row_sum:.1e}")

# --- T2: constant-stress membrane patch (ν=0, non-conforming layout) --------------
print("\nT2  constant-stress membrane patch (ν=0): any in-plane strain crosses exactly")
EPS = np.array([[3.1e-4, 1.2e-4, 0.0],
                [1.2e-4, -2.4e-4, 0.0],
                [0.0, 0.0, 0.0]])                  # constant in-plane strain, ν=0
U = master_state(EDGE, lambda x: EPS @ x, lambda x: np.zeros(3))
exact = XS @ EPS.T
e1 = np.abs(predict(ROWS, U) - exact).max()
check("membrane patch exact through the tie", e1 < 1e-15, f"err={e1:.2e}")

# --- T3: HEADLINE — constant-MOMENT patch between dissimilar elements -------------
print("\nT3  HEADLINE: constant moment crosses shell↔solid (θ=κx / u_a=−κxz, ν=0)")
KAPPA = 0.05
# shell side at the interface x=X0: w=½κx², section rotation about y with
# u_x = −z·w'(x) ⇒ θ_y = −κx (u = θ×d arm kinematics); solid side: u_x = −κxz,
# u_z = ½κx² (exact pure-bending elasticity at ν=0).
U = master_state(EDGE,
                 lambda x: np.array([0.0, 0.0, 0.5 * KAPPA * x[0] ** 2]),
                 lambda x: np.array([0.0, -KAPPA * x[0], 0.0]))
exact = np.stack([-KAPPA * X0 * XS[:, 2],
                  np.zeros(len(XS)),
                  np.full(len(XS), 0.5 * KAPPA * X0 ** 2)], axis=1)
e1 = np.abs(predict(ROWS, U) - exact).max()
check("constant-moment patch EXACT (linear-in-thickness u_a reproduced)",
      e1 < 1e-15, f"err={e1:.2e}")

# --- T4: work-conjugacy / static equivalence --------------------------------------
print("\nT4  work-conjugacy: Cᵀλ momentum-exact, θ columns a pure couple, TᵀKT symmetric")
LAM = rng.normal(size=(len(SLAVES), 3)) * 2.7      # arbitrary slave reactions
Fm = np.zeros(len(EDGE) * 6)
for (Crow, _, _), lam in zip(ROWS, LAM):
    Fm += Crow.T @ lam                             # f_master = Cᵀλ (both handlers)
F = Fm.reshape(-1, 6)[:, :3]                       # nodal forces
M = Fm.reshape(-1, 6)[:, 3:]                       # nodal couples (from the θ columns)
e1 = np.abs(F.sum(axis=0) - LAM.sum(axis=0)).max()
check("linear momentum exact (PoU on the translation columns)", e1 < 1e-13,
      f"err={e1:.2e}")
ang_m = np.cross(np.array(EDGE), F).sum(axis=0) + M.sum(axis=0)
ang_s = np.cross(XS, LAM).sum(axis=0)
e2 = np.abs(ang_m - ang_s).max()
check("angular momentum exact (θ columns carry the transferred couple N_j·(d×λ))",
      e2 < 1e-13, f"err={e2:.2e}")
# the θ columns deliver a PURE couple: they land on rotation DOFs only — the nodal
# FORCE block is independent of the arms (rebuild with θ columns zeroed ⇒ same F).
Fm0 = np.zeros(len(EDGE) * 6)
for xs, lam in zip(SLAVES, LAM):
    C0, _, _ = plane_section_rows(EDGE, xs, zero_theta=True)
    Fm0 += C0.T @ lam
e3 = np.abs(Fm0.reshape(-1, 6)[:, :3] - F).max()
check("θ columns contribute ZERO net force (pure couple)", e3 == 0.0, f"delta={e3:.2e}")
# reduced stiffness: eliminate slaves via u_s = C u_m ⇒ T = [I; Call], K_red = TᵀKT
nm, ns = len(EDGE) * 6, len(SLAVES) * 3
Call = np.vstack([C for C, _, _ in ROWS])
A = rng.normal(size=(nm + ns, nm + ns))
K = A.T @ A + np.eye(nm + ns)                      # synthetic SPD "model" stiffness
T = np.vstack([np.eye(nm), Call])
K_red = T.T @ K @ T
e4 = np.abs(K_red - K_red.T).max() / np.abs(K_red).max()
check("reduced stiffness TᵀKT symmetric", e4 < 1e-14, f"asym={e4:.2e}")

# --- T5: 3D orientation invariance + kinked polyline (OQ-7) -----------------------
print("\nT5  rotated frame (skew edge/arm): exactness holds, drilling free via NULL SPACE")
th = np.array([0.4, -0.75, 0.3])
Rm = (np.eye(3) + np.sin(np.linalg.norm(th)) * skew(th / np.linalg.norm(th)) +
      (1 - np.cos(np.linalg.norm(th))) * skew(th / np.linalg.norm(th)) @
      skew(th / np.linalg.norm(th)))               # Rodrigues
edge_r = [Rm @ x for x in EDGE]
slaves_r = [Rm @ x for x in SLAVES]
rows_r = [plane_section_rows(edge_r, xs) for xs in slaves_r]
XSr = np.array(slaves_r)
# rotated bending state: X' = RX, u'(X') = R u(X), θ' = R θ
U = master_state(edge_r,
                 lambda xr: Rm @ np.array([0.0, 0.0, 0.5 * KAPPA * (Rm.T @ xr)[0] ** 2]),
                 lambda xr: Rm @ np.array([0.0, -KAPPA * (Rm.T @ xr)[0], 0.0]))
exact = (np.stack([-KAPPA * X0 * XS[:, 2],
                   np.zeros(len(XS)),
                   np.full(len(XS), 0.5 * KAPPA * X0 ** 2)], axis=1)) @ Rm.T
e1 = np.abs(predict(rows_r, U) - exact).max()
check("rotated constant-moment state exact (global-DOF rows, skew arm)",
      e1 < 1e-14, f"err={e1:.2e}")
# drilling in a skew frame: no global column is near-zero (the filter keeps all),
# but −[d]× annihilates its own axis ⇒ θ += α·d̂ leaves every prediction unchanged.
n_r = Rm @ np.array([0.0, 0.0, 1.0])               # rotated shell normal = arm direction
Urand = rng.normal(size=len(EDGE) * 6) * 0.01
Udrill = Urand.copy()
for j in range(len(EDGE)):
    Udrill[j * 6 + 3:j * 6 + 6] += 7.7 * n_r
e2 = np.abs(predict(rows_r, Urand) - predict(rows_r, Udrill)).max()
check("skew drilling free via the null space of −[d]× (structural, not the filter)",
      e2 < 1e-14, f"delta={e2:.2e}")
# kinked polyline, per-node frozen arms (OQ-7): rigid modes must still be exact.
edge_k = [np.array([X0, 0.0, 0.0]), np.array([X0, 0.5, 0.0]),
          np.array([X0 + 0.25, 1.0, 0.0])]         # 2nd segment kinks in x
slaves_k = [np.array([X0, 0.2, z]) for z in SZ] + \
           [np.array([X0 + 0.15, 0.8 + 0.075, z]) for z in SZ]  # off both segments
rows_k = [plane_section_rows(edge_k, xs) for xs in slaves_k]
err = 0.0
for mode in range(6):
    u0 = np.zeros(3); om = np.zeros(3)
    (u0 if mode < 3 else om)[mode % 3] = 0.011
    U = master_state(edge_k, lambda x: u0 + np.cross(om, x - P0), lambda x: om)
    exact = u0 + np.cross(om, np.array(slaves_k) - P0)
    err = max(err, np.abs(predict(rows_k, U) - exact).max())
check("kinked polyline + per-node frozen arms: all 6 rigid modes exact (OQ-7)",
      err < 1e-15, f"err={err:.2e}")

# --- T6: honest limit — Poisson thickness-stretch suppression (OQ-2) --------------
print("\nT6  honest limit: through-thickness Poisson stretch suppressed — misfit = ν·ε·|z|")
EPS0 = 4.0e-4                                      # axial stretch along the edge (y)


def poisson_misfit(nu, c):
    edge = [np.array([X0, y, 0.0]) for y in MY]
    slaves = [np.array([X0, y, z]) for y in SY for z in (-c, 0.5 * c, c)]
    rows = [plane_section_rows(edge, xs) for xs in slaves]
    # solid field u = (−νε x, ε y, −νε z); shell masters at z=0, membrane ⇒ θ=0
    U = master_state(edge, lambda x: np.array([-nu * EPS0 * x[0], EPS0 * x[1], 0.0]),
                     lambda x: np.zeros(3))
    Xs = np.array(slaves)
    exact = np.stack([-nu * EPS0 * Xs[:, 0], EPS0 * Xs[:, 1], -nu * EPS0 * Xs[:, 2]],
                     axis=1)
    mis = predict(rows, U) - exact
    return np.abs(mis[:, 2]).max(), np.abs(mis[:, :2]).max()


m0, ip0 = poisson_misfit(0.0, C_HALF)
m1, ip1 = poisson_misfit(0.15, C_HALF)
m2, ip2 = poisson_misfit(0.30, C_HALF)
m2t, _ = poisson_misfit(0.30, 2.0 * C_HALF)
check("EXACT at ν=0 (in-plane always exact)", m0 < 1e-15 and max(ip0, ip1, ip2) < 1e-15,
      f"ν=0 misfit={m0:.2e}")
check("misfit linear in ν (suppression = ν·ε·|z|)",
      abs(m2 / m1 - 2.0) < 1e-6 and abs(m2 - 0.30 * EPS0 * C_HALF) < 1e-15,
      f"ν-ratio={m2/m1:.3f} (→2), max={m2:.2e} vs νεc={0.30*EPS0*C_HALF:.2e}")
check("misfit linear in t (a boundary layer, not growing with the model)",
      abs(m2t / m2 - 2.0) < 1e-6, f"t-ratio={m2t/m2:.3f} (→2)")

# --- T7: honest limit — Timoshenko shear warping O(γ·t) (OQ-2) --------------------
print("\nT7  honest limit: parabolic shear warping f(z)=(γ/2)(z−z³/c²) — misfit O(γ·t)")
# Timoshenko section state at the interface: u_x = −z·θ_b + f(z), u_z = w0. The
# plane-section tie reproduces −z·θ_b exactly (T3) and CANNOT represent f(z);
# misfit = |f(z)|, analytic max γc/(3√3) at z = c/√3 — independent of in-plane (y)
# mesh, linear in γ AND in c, while the RESULTANT transfer stays exact (T4 identity).
W0, THB = 2.3e-3, 1.7e-3


def shear_misfit(gam, c, ny):
    edge = [np.array([X0, y, 0.0]) for y in MY]
    ys = np.linspace(0.05, 0.95, ny)
    zs = np.linspace(-c, c, 21)                    # dense in z to catch the extremum
    slaves = [np.array([X0, y, z]) for y in ys for z in zs]
    rows = [plane_section_rows(edge, xs) for xs in slaves]
    U = master_state(edge, lambda x: np.array([0.0, 0.0, W0]),
                     lambda x: np.array([0.0, -THB, 0.0]))
    Xs = np.array(slaves)
    f = 0.5 * gam * (Xs[:, 2] - Xs[:, 2] ** 3 / c ** 2)
    exact = np.stack([-THB * Xs[:, 2] + f, np.zeros(len(Xs)), np.full(len(Xs), W0)],
                     axis=1)
    return np.abs(predict(rows, U) - exact).max()


GAM = 8.0e-4
s1 = shear_misfit(GAM, C_HALF, 5)
s2 = shear_misfit(2 * GAM, C_HALF, 5)
s3 = shear_misfit(GAM, 2 * C_HALF, 5)
s4 = shear_misfit(GAM, C_HALF, 20)                 # 4× in-plane refinement
bound = GAM * C_HALF / (3.0 * np.sqrt(3.0))
check("misfit linear in γ", abs(s2 / s1 - 2.0) < 1e-3, f"γ-ratio={s2/s1:.4f} (→2)")
check("misfit linear in t (O(γ·t))", abs(s3 / s1 - 2.0) < 1e-3,
      f"t-ratio={s3/s1:.4f} (→2)")
check("misfit does NOT shrink with in-plane refinement (a thickness-local layer)",
      abs(s4 / s1 - 1.0) < 1e-12, f"refine-ratio={s4/s1:.6f} (→1)")
check("misfit ≈ the analytic bound γc/(3√3), ≪ the section displacement scale",
      s1 <= bound * 1.001 and s1 > 0.95 * bound and s1 < 0.05 * abs(THB) * X0,
      f"misfit={s1:.2e} vs bound={bound:.2e}")
# resultants stay exact even for shear-consistent reactions (parabolic τ(z) lumping)
edge = EDGE
zs = np.linspace(-C_HALF, C_HALF, 9)
slaves7 = [np.array([X0, 0.5, z]) for z in zs]
rows7 = [plane_section_rows(edge, xs) for xs in slaves7]
tau = (1.0 - (zs / C_HALF) ** 2)                   # parabolic shear traction shape
lam7 = np.stack([np.zeros(9), np.zeros(9), tau], axis=1)
Fm = np.zeros(len(edge) * 6)
for (Crow, _, _), lam in zip(rows7, lam7):
    Fm += Crow.T @ lam
F = Fm.reshape(-1, 6)[:, :3]; M = Fm.reshape(-1, 6)[:, 3:]
e1 = np.abs(F.sum(axis=0) - lam7.sum(axis=0)).max()
e2 = np.abs(np.cross(np.array(edge), F).sum(axis=0) + M.sum(axis=0)
            - np.cross(np.array(slaves7), lam7).sum(axis=0)).max()
check("resultant transfer stays exact under the shear distribution",
      max(e1, e2) < 1e-13, f"errs={e1:.2e},{e2:.2e}")

# --- T8: contrast — θ columns zeroed = the hinge ----------------------------------
print("\nT8  contrast: translations-only tie (-dof 3 1 2 3) passes T2 but FAILS T3 (hinge)")
rows_h = [plane_section_rows(EDGE, xs, zero_theta=True) for xs in SLAVES]
U = master_state(EDGE, lambda x: EPS @ x, lambda x: np.zeros(3))
e1 = np.abs(predict(rows_h, U) - XS @ EPS.T).max()
check("membrane patch still exact without θ columns", e1 < 1e-15, f"err={e1:.2e}")
U = master_state(EDGE,
                 lambda x: np.array([0.0, 0.0, 0.5 * KAPPA * x[0] ** 2]),
                 lambda x: np.array([0.0, -KAPPA * x[0], 0.0]))
exact = np.stack([-KAPPA * X0 * XS[:, 2], np.zeros(len(XS)),
                  np.full(len(XS), 0.5 * KAPPA * X0 ** 2)], axis=1)
e2 = np.abs(predict(rows_h, U) - exact).max()
check("constant-moment patch FAILS (u_a=−κxz unreachable ⇒ hinge; misfit=κ·x·c)",
      abs(e2 - KAPPA * X0 * C_HALF) < 1e-15, f"misfit={e2:.2e} = κ·x0·c={KAPPA*X0*C_HALF:.2e}")

# ===================================================================== summary
print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'} — "
      "P4 shell-to-solid plane-section tie (b-B): rigid+membrane+constant-moment "
      "exact, drilling auto-free (filter axis-aligned / null-space skew), momentum "
      "identities exact, honest limits scale as ν·t and γ·t as predicted (OQ-2), "
      "translations-only contrast is a hinge.")
sys.exit(0 if _fails == 0 else 1)
