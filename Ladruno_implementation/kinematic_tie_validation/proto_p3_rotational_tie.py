"""ADR-62 P3 ORACLE — shell / ROTATIONAL kinematic mesh-tie (ndf-6 shell-to-shell).

P1/P2 tie the TRANSLATIONAL DOFs of a non-conforming interface: each slave node's
displacement is condensed onto master nodes as  u_s = Σ_k P_sk u_{m,k}  (P = the
collocation shape weights N_i, or the mortar row D⁻¹M). P3 extends the SAME per-slave
transfer P to the ROTATIONAL DOFs of an ndf-6 shell node:

        θ_s = Σ_k P_sk θ_{m,k}          (same P, applied per rotational DOF)

so the generator emits three extra EQ_Constraint rows per slave node (DOFs 4,5,6),
each with the identical P-row coefficients it already uses for DOFs 1,2,3. The shipped
LadrunoProjectionHandler (ADR-30) enforces them unchanged — it treats every (node,dof)
as an independent union-find vertex, so a rotational row is just more master-only
vertices (verified in LadrunoProjectionHandler.cpp:326-366).

WHY THIS IS CORRECT — and its PRECISE LIMIT (the honest scope of v1):

  A constant-curvature (constant-moment) shell state has a rotation field that is
  LINEAR in position and a transverse-displacement field that is QUADRATIC:
        θ(x) = θ0 + κ·x            (linear   — curvature κ = ∂θ/∂x constant)
        w(x) = w0 + θ0·x + ½κ·x²   (quadratic — zero transverse shear: ∂w/∂x = θ)

  The transfer P is LINEARLY complete (P reproduces a + b·x exactly — proven in P2 T3)
  but NOT quadratically complete. Therefore:

   • The ROTATION field (linear) crosses ANY non-matching tie EXACTLY.           ← NEW
   • The TRANSVERSE DISPLACEMENT (quadratic) is reproduced exactly ONLY where it is
     at most linear ALONG the interface line. For CYLINDRICAL bending whose curvature
     axis is PARALLEL to the (straight) interface — the shell-strip patch — both w and
     θ are CONSTANT along the interface, so the linear-complete tie reproduces them
     exactly and the constant-moment state crosses the tie EXACTLY. This is the v1
     scope (and the FE test in tests/test_ladrunoTie_shell.py).
   • If the curvature makes w vary QUADRATICALLY *along* the interface, the tie leaves
     an O(h²) residual on w (never on θ) — same order as the shell discretization
     error, vanishing on refinement and for matching meshes. Removing it needs a
     Hermite / rotation-consistent w–θ transfer (couple w to θ across the interface
     tangent) — deferred as P3.1, NOT a straight reuse of P.

  Scope v1 = shell-to-shell, same ndf, same P on rotations. Shell-to-SOLID (which needs
  θ = ½∇×u to synthesize a rotation from the solid translation field, a DIFFERENT
  operator) is a later rung — not covered here.

Gates (falsifiable, numeric):
  T1  P is linearly complete (precondition inherited from P1/P2): P·1 = 1 and
      P·(linear) = linear at the slave nodes.
  T2  RIGID ROTATION: a constant rotation field θ ≡ c transfers exactly (partition of
      unity on the rotational DOFs) ⇒ a rigid-body rotation of the master crosses.
  T3  CONSTANT CURVATURE — rotation: a linear rotation field θ = θ0 + κ·x is reproduced
      EXACTLY at the slave nodes by the SAME P (this is the headline P3 result).
  T4  ALIGNED bending patch (curvature axis ∥ interface): w and θ are CONSTANT along
      the interface ⇒ BOTH cross exactly ⇒ the constant-moment strip patch is exact.
  T5  THE LIMIT (honest): curvature varying ALONG the interface ⇒ quadratic w leaves a
      NONZERO residual that is O(h²) (4× smaller when the interface is refined h→h/2),
      while θ stays exact. Proves we bound the limit rather than hide it.

Run: <pythoncore-3.12> proto_p3_rotational_tie.py   (numpy only, no build)
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


# ============================================================ collocation transfer P
# The shell interface of a flat strip IS a 1D line; the tie condenses each slave node
# onto the master edge it projects into. We build the P1 COLLOCATION P on that line
# (linear master edge, shape N0=1-s, N1=s). This is the exact operator the generator
# emits for a shell-to-shell collocation tie; its 2D-mortar sibling is proven linearly
# complete in proto_p2 (T3) and inherits every result below identically.
def collocation_P(master_s, slave_s):
    """P[i,k] = master shape weight of master node k at slave node i's projection.
    master_s: sorted master-node arc positions (defines the linear master segments).
    slave_s : slave-node arc positions (each lands in one master segment)."""
    master_s = np.asarray(master_s, float)
    P = np.zeros((len(slave_s), len(master_s)))
    for i, s in enumerate(slave_s):
        # locate the master segment [a,b] containing s (clamp at the ends)
        seg = np.searchsorted(master_s, s) - 1
        seg = min(max(seg, 0), len(master_s) - 2)
        a, b = master_s[seg], master_s[seg + 1]
        xi = 0.0 if b == a else (s - a) / (b - a)
        P[i, seg] += 1.0 - xi
        P[i, seg + 1] += xi
    return P


# a genuinely NON-MATCHING interface: master edge split at s=0,1 (2 nodes); slave nodes
# at s=0, 0.5, 1 — the s=0.5 slave is non-conforming (0.5/0.5 master weights).
MASTER_S = [0.0, 1.0]
SLAVE_S = [0.0, 0.5, 1.0]
P = collocation_P(MASTER_S, SLAVE_S)

print("ADR-62 P3 oracle — shell/rotational tie (θ_s = Σ P_sk θ_{m,k}, same P)")
print(f"\n     non-matching interface  master s={MASTER_S}  slave s={SLAVE_S}")
print(f"     collocation P =\n{P}")

# --- T1: P linearly complete (precondition) -----------------------------------
print("\nT1  P linearly complete (inherited P1/P2 precondition)")
check("Σ_k P_ik = 1 (partition of unity)", np.allclose(P.sum(1), 1.0, atol=1e-14),
      f"max|ΣP-1|={np.abs(P.sum(1) - 1).max():.2e}")
lin_m = 0.7 + 2.3 * np.asarray(MASTER_S)
lin_s = 0.7 + 2.3 * np.asarray(SLAVE_S)
check("P·(linear on masters) = linear at slaves", np.allclose(P @ lin_m, lin_s, atol=1e-14),
      f"err={np.abs(P @ lin_m - lin_s).max():.2e}")

# --- T2: rigid rotation (constant θ field transfers exactly) ------------------
print("\nT2  rigid rotation: a CONSTANT rotation field crosses exactly")
theta_const = 0.031          # rad, same at every master node
tm = np.full(len(MASTER_S), theta_const)
ts = P @ tm
check("θ_s ≡ θ_master (partition of unity on rotations)",
      np.allclose(ts, theta_const, atol=1e-14), f"max|θ_s-θ|={np.abs(ts - theta_const).max():.2e}")

# --- T3: constant curvature — rotation reproduced EXACTLY (the headline) -------
print("\nT3  constant curvature: linear rotation field θ=θ0+κ·x reproduced EXACTLY")
theta0, kappa = -0.004, 0.05
tm = theta0 + kappa * np.asarray(MASTER_S)
ts_exact = theta0 + kappa * np.asarray(SLAVE_S)
ts_tie = P @ tm
err_theta = np.abs(ts_tie - ts_exact).max()
check("θ_s (tie) == θ_s (exact linear field)", err_theta < 1e-13,
      f"max nodal θ err={err_theta:.2e}  (incl. the 0.5/0.5 non-conforming node)")

# --- T4: ALIGNED bending patch — w and θ CONSTANT along the interface ---------
print("\nT4  aligned bending (curvature axis ∥ interface): w and θ CONSTANT along it")
# strip bends about y (∥ interface); along the interface line s (the y-direction) BOTH
# the transverse displacement w and the rotation θ are position-independent constants.
w_iface, th_iface = 0.5 * kappa * 1.0**2, kappa * 1.0    # values at the interface x=1
wm = np.full(len(MASTER_S), w_iface)
thm = np.full(len(MASTER_S), th_iface)
ew = np.abs(P @ wm - w_iface).max()
et = np.abs(P @ thm - th_iface).max()
check("transverse w crosses exactly (constant along interface)", ew < 1e-14, f"err={ew:.2e}")
check("rotation θ crosses exactly (constant along interface)", et < 1e-14, f"err={et:.2e}")
print("      → the constant-moment shell-strip patch crosses the non-matching tie EXACTLY")

# --- T5: THE LIMIT — quadratic w ALONG the interface leaves an O(h²) residual --
print("\nT5  the honest limit: curvature ALONG the interface ⇒ quadratic w ⇒ O(h²) residual")


def w_residual(n_master):
    """Refine the master edge into n_master segments (h = 1/n_master); slave nodes at the
    segment midpoints (maximally non-conforming). Return the worst |P·w_quad - w_quad|
    for a field w = ½·c·s² that varies QUADRATICALLY *along* the interface, and the
    matching worst θ error for the LINEAR θ = c·s field (must stay ~0)."""
    ms = np.linspace(0.0, 1.0, n_master + 1)
    ss = 0.5 * (ms[:-1] + ms[1:])            # midpoints — non-conforming
    Ph = collocation_P(ms, ss)
    c = 0.05
    wq_m, wq_s = 0.5 * c * ms**2, 0.5 * c * ss**2      # quadratic
    th_m, th_s = c * ms, c * ss                        # linear
    return np.abs(Ph @ wq_m - wq_s).max(), np.abs(Ph @ th_m - th_s).max()


rw1, rt1 = w_residual(2)      # h = 1/2
rw2, rt2 = w_residual(4)      # h = 1/4
rw3, rt3 = w_residual(8)      # h = 1/8
check("linear θ ALONG the interface stays EXACT under refinement",
      max(rt1, rt2, rt3) < 1e-13, f"max θ err={max(rt1, rt2, rt3):.2e}")
check("quadratic w ALONG the interface leaves a NONZERO residual (expected limit)",
      rw1 > 1e-6, f"w residual(h=1/2)={rw1:.3e}")
ratio1, ratio2 = rw1 / rw2, rw2 / rw3
check("w residual is O(h²) (halving h ⇒ ~4× smaller)",
      abs(ratio1 - 4.0) < 0.4 and abs(ratio2 - 4.0) < 0.4,
      f"ratios ≈ {ratio1:.2f}, {ratio2:.2f} (→4)")
print("      → θ is exact for ALL curvatures; only w-along-interface is O(h²)-approximate.")
print("        v1 scope = aligned strip patch (T4, EXACT). Hermite w–θ transfer = P3.1.")

# ===================================================================== summary
print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'} — "
      "P3 rotational tie: the SAME linear-complete P transfers rotations "
      "(θ_s=Σ P_sk θ_{m,k}); rigid rotation + constant-curvature rotation cross EXACTLY; "
      "the aligned bending patch is exact; along-interface w is O(h²) (P3.1).")
sys.exit(0 if _fails == 0 else 1)
