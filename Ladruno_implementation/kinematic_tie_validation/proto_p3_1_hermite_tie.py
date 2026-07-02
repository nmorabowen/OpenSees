"""ADR-62 P3.1 ORACLE — Hermite (rotation-consistent) w–θ shell transfer.

P3 ties every DOF of an ndf-6 shell node with the SAME linearly-complete transfer P.
Its proven limit (P3 oracle T5): a curvature varying ALONG the interface makes the
transverse displacement w QUADRATIC along the interface line, and a linear-complete
P leaves an O(h²) residual on w (never on θ). P3.1 removes it: reconstruct w from a
CUBIC HERMITE combination of the master nodal w AND the master interface-tangent
slope φ = dw/ds, instead of the straight linear interpolation.

THE TRANSFER (1D along the straight interface, master segment [a,b], h=b−a, ξ=(s−a)/h):

    w_s  = H1(ξ)·w_a + H2(ξ)·φ_a + H3(ξ)·w_b + H4(ξ)·φ_b          (w row)
    φ_s  = H1'·w_a  + H2'·φ_a  + H3'·w_b  + H4'·φ_b               (tangent-slope row)
    H1 = 1−3ξ²+2ξ³   H2 = h(ξ−2ξ²+ξ³)   H3 = 3ξ²−2ξ³   H4 = h(ξ³−ξ²)

where φ is the SLOPE OF w ALONG THE INTERFACE. For a shell with normal n and rotation
vector θ (right-hand rule, u = θ×r), the transverse-displacement gradient obeys
∇w = n×θ, so the slope along the unit tangent t is

    φ = dw/ds = t·(n×θ) = θ·(t×n) = θ·e ,   e ≡ t×n  (in-plane, ⊥ interface).

So φ_k is a fixed LINEAR COMBINATION of the master node's three global rotation DOFs
(coefficients e_x,e_y,e_z) — and the Hermite row is exactly an EQ_Constraint: a list
of (masterNode, dof, coef) triples. The remaining slave DOF rows split by projector:

    u_s = (I−n⊗n)·Σ P_k u_k  +  n·[Σ_k H_w,k (n·u_k) + H_φ,k (e·θ_k)]
    θ_s = (I−e⊗e)·Σ P_k θ_k  +  e·[Σ_k H'_w,k (n·u_k) + H'_φ,k (e·θ_k)]

(in-plane translations + non-slope rotations keep the linear P; only the n-component
of u and the e-component of θ go Hermite). KIRCHHOFF assumption: the slope equals the
rotation (zero transverse shear). A shear-deformable (Mindlin) state has φ = dw/ds − γ;
the Hermite row then carries an O(γ·h) error — bounded honestly in T6.

WHY THIS STAYS A LadrunoTie-LEVEL TRANSFORM (the OQ, answered):
  • EQ_Constraint natively stores per-retained (node, dof, coef) triples — a mixed
    w←(w,θ) row needs NO constraint-class change (EQ_Constraint.h:37-45).
  • LadrunoProjectionHandler::buildGroups builds one union-find vertex per (node,dof)
    and walks arbitrary retained-DOF lists — NO retained-dof==constrained-dof
    assumption (LadrunoProjectionHandler.cpp EQ path). Mixed rows are just more edges.
  • COLLOCATION mode: the closest-point projection already yields the master edge +
    parameter ξ; evaluating H1..H4 is emission-time arithmetic. NO kernel change.
  • MORTAR mode is the exception: a weak-form Hermite needs the Hermite basis inside
    the M integral, but LadrunoMortarKernel::integratePair only returns ACCUMULATED
    D/M (PairResult has no per-GP hook) ⇒ mortar-Hermite = kernel extension ⇒ defer
    (P3.1 scope = collocation/edge ties, where the along-interface bending case lives).

Gates (falsifiable, numeric):
  T1  Rigid modes: constant w (φ=0) and rigid rotation (w linear, φ const) cross the
      Hermite rows EXACTLY (PoU + linear completeness preserved — no regression vs P).
  T2  HEADLINE: constant curvature ALONG the interface (w=½κs² quadratic, φ=κs) is
      EXACT through the Hermite rows — the SAME field the linear P misses at O(h²)
      (shown side by side).
  T3  Cubic completeness: w cubic + φ=w' EXACT (both rows) — full Hermite order.
  T4  Rates: for a QUARTIC w the Hermite residual is O(h⁴) vs the linear P's O(h²)
      (strict dominance, not a wash).
  T5  3D orientation invariance: an INCLINED interface (general t, n) with the rows
      built at the GLOBAL-DOF level from the projector split above reproduces the
      rotated constant-curvature bending state exactly on all 6 DOFs.
  T6  Honest limit (Kirchhoff assumption): Mindlin data (φ = dw/ds − γ) leaves an
      O(γ·h) error on the w row — linear in the shear angle γ AND in h (thin-shell /
      refinement safe), while the linear-P rows are γ-insensitive.

Run: <pythoncore-3.12> proto_p3_1_hermite_tie.py   (numpy only, no build)
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


# ============================================================ transfer operators
def linear_P(master_s, slave_s):
    """P1/P3 linear collocation transfer (per-DOF row of master shape weights)."""
    master_s = np.asarray(master_s, float)
    P = np.zeros((len(slave_s), len(master_s)))
    for i, s in enumerate(slave_s):
        seg = min(max(np.searchsorted(master_s, s) - 1, 0), len(master_s) - 2)
        a, b = master_s[seg], master_s[seg + 1]
        xi = 0.0 if b == a else (s - a) / (b - a)
        P[i, seg] += 1.0 - xi
        P[i, seg + 1] += xi
    return P


def hermite_rows(master_s, slave_s):
    """P3.1 Hermite transfer. Returns (Hw, Hphi, dHw, dHphi):
       w_s   = Hw  @ w_m  + Hphi  @ phi_m      (transverse-displacement row)
       phi_s = dHw @ w_m  + dHphi @ phi_m      (tangent-slope row)
    Each row touches ONE master segment (4 nonzeros) — local, EQ-triple-shaped."""
    master_s = np.asarray(master_s, float)
    ns, nm = len(slave_s), len(master_s)
    Hw, Hp = np.zeros((ns, nm)), np.zeros((ns, nm))
    dHw, dHp = np.zeros((ns, nm)), np.zeros((ns, nm))
    for i, s in enumerate(slave_s):
        seg = min(max(np.searchsorted(master_s, s) - 1, 0), nm - 2)
        a, b = master_s[seg], master_s[seg + 1]
        h = b - a
        xi = 0.0 if h == 0.0 else (s - a) / h
        # cubic Hermite basis on [a,b] (slope DOFs carry the segment length h)
        Hw[i, seg]      = 1 - 3*xi**2 + 2*xi**3
        Hp[i, seg]      = h * (xi - 2*xi**2 + xi**3)
        Hw[i, seg + 1]  = 3*xi**2 - 2*xi**3
        Hp[i, seg + 1]  = h * (xi**3 - xi**2)
        # d/ds = (1/h) d/dξ
        dHw[i, seg]     = (-6*xi + 6*xi**2) / h
        dHp[i, seg]     = 1 - 4*xi + 3*xi**2
        dHw[i, seg + 1] = (6*xi - 6*xi**2) / h
        dHp[i, seg + 1] = 3*xi**2 - 2*xi
    return Hw, Hp, dHw, dHp


def grids(n):
    """n-segment master on [0,1]; slaves at segment midpoints (maximally non-conforming)."""
    ms = np.linspace(0.0, 1.0, n + 1)
    return ms, 0.5 * (ms[:-1] + ms[1:])


MS, SS = grids(2)
P = linear_P(MS, SS)
HW, HP, DHW, DHP = hermite_rows(MS, SS)
print("ADR-62 P3.1 oracle — Hermite w–θ transfer (collocation/edge ties)")
print(f"\n     master s={MS}  slave s={SS} (midpoints, non-conforming)")

# --- T1: rigid modes preserved --------------------------------------------------
print("\nT1  rigid modes cross the Hermite rows exactly (no regression vs linear P)")
w_c = np.full(len(MS), 0.37)                       # constant w, zero slope
e1 = np.abs(HW @ w_c + HP @ np.zeros(len(MS)) - 0.37).max()
e2 = np.abs(DHW @ w_c).max()                       # slope of a constant = 0
check("constant w: w row PoU + slope row → 0", max(e1, e2) < 1e-14,
      f"errs={e1:.2e},{e2:.2e}")
c = 0.041                                          # rigid rotation: w=c·s, φ=c
w_m, p_m = c * MS, np.full(len(MS), c)
e1 = np.abs(HW @ w_m + HP @ p_m - c * SS).max()
e2 = np.abs(DHW @ w_m + DHP @ p_m - c).max()
check("rigid rotation (w linear, φ const): both rows exact", max(e1, e2) < 1e-14,
      f"errs={e1:.2e},{e2:.2e}")

# --- T2: HEADLINE — constant curvature ALONG the interface ----------------------
print("\nT2  headline: quadratic w ALONG the interface — Hermite EXACT, linear P O(h²)")
kappa = 0.05
w_m, p_m = 0.5 * kappa * MS**2, kappa * MS         # w=½κs², φ=κs (Kirchhoff-consistent)
w_s, p_s = 0.5 * kappa * SS**2, kappa * SS
err_H = np.abs(HW @ w_m + HP @ p_m - w_s).max()
err_P = np.abs(P @ w_m - w_s).max()
err_slope = np.abs(DHW @ w_m + DHP @ p_m - p_s).max()
check("Hermite w row reproduces the quadratic w EXACTLY", err_H < 1e-15,
      f"err={err_H:.2e}")
check("Hermite slope row reproduces φ=κs EXACTLY", err_slope < 1e-15,
      f"err={err_slope:.2e}")
check("the SAME field through the linear P is O(h²)-wrong (the P3 T5 limit)",
      err_P > 1e-6, f"linear-P err={err_P:.3e}")

# --- T3: cubic completeness -----------------------------------------------------
print("\nT3  cubic completeness (full Hermite order)")
co = np.array([0.3, -0.7, 0.4, 0.9])               # w = c0 + c1 s + c2 s² + c3 s³
w_m = co[0] + co[1]*MS + co[2]*MS**2 + co[3]*MS**3
p_m = co[1] + 2*co[2]*MS + 3*co[3]*MS**2
w_s = co[0] + co[1]*SS + co[2]*SS**2 + co[3]*SS**3
p_s = co[1] + 2*co[2]*SS + 3*co[3]*SS**2
e1 = np.abs(HW @ w_m + HP @ p_m - w_s).max()
e2 = np.abs(DHW @ w_m + DHP @ p_m - p_s).max()
check("cubic w + consistent φ: both rows exact", max(e1, e2) < 1e-13,
      f"errs={e1:.2e},{e2:.2e}")

# --- T4: rates — Hermite O(h⁴) vs linear O(h²) on a quartic ---------------------
print("\nT4  refinement rates on a quartic w (strict dominance)")


def rates(n):
    ms, ss = grids(n)
    Ph = linear_P(ms, ss)
    hw, hp, _, _ = hermite_rows(ms, ss)
    w_m, p_m, w_s = ms**4, 4*ms**3, ss**4
    return (np.abs(hw @ w_m + hp @ p_m - w_s).max(),
            np.abs(Ph @ w_m - w_s).max())


# the quartic's curvature-max sits at the LAST segment and drifts outward with n, so
# the linear rate approaches 4 from below — measure at n large enough to be asymptotic
(h1, l1), (h2, l2), (h3, l3) = rates(8), rates(16), rates(32)
rH1, rH2 = h1 / h2, h2 / h3
rL1, rL2 = l1 / l2, l2 / l3
check("Hermite residual is O(h⁴) (halving h ⇒ ~16×)",
      abs(rH1 - 16) < 3 and abs(rH2 - 16) < 3, f"ratios≈{rH1:.1f},{rH2:.1f} (→16)")
check("linear-P residual is O(h²) (halving h ⇒ ~4×)",
      abs(rL1 - 4) < 0.5 and abs(rL2 - 4) < 0.5, f"ratios≈{rL1:.2f},{rL2:.2f} (→4)")

# --- T5: 3D orientation invariance (global-DOF row construction) ----------------
print("\nT5  inclined interface: GLOBAL-DOF rows (projector split) reproduce bending")
# unit tangent t (interface direction), shell normal n, in-plane transverse e = t×n
t = np.array([2.0, 1.0, 0.5]); t /= np.linalg.norm(t)
n = np.array([-1.0, 2.0, 0.0]); n -= (n @ t) * t; n /= np.linalg.norm(n)
e = np.cross(t, n)
# bending along the interface: u = w(s)·n, θ = φ(s)·e with w=½κs², φ=κs  (∇w = n×θ ✓)
Xm = np.outer(MS, t)                                # master nodes on the interface line
Xs = np.outer(SS, t)                                # slave nodes on the same line
Um = 0.5 * kappa * MS[:, None]**2 * n               # master translations (3 per node)
Tm = kappa * MS[:, None] * e                        # master rotations   (3 per node)
# build the slave DOFs through the P3.1 global-DOF transfer:
#   u_s = (I−n⊗n)(P u) + n·[Hw (u·n) + Hp (θ·e)] ;  θ_s = (I−e⊗e)(P θ) + e·[dHw (u·n) + dHp (θ·e)]
w_nod, phi_nod = Um @ n, Tm @ e
Us = (P @ Um) @ (np.eye(3) - np.outer(n, n)).T + np.outer(HW @ w_nod + HP @ phi_nod, n)
Ts = (P @ Tm) @ (np.eye(3) - np.outer(e, e)).T + np.outer(DHW @ w_nod + DHP @ phi_nod, e)
Us_exact = 0.5 * kappa * SS[:, None]**2 * n
Ts_exact = kappa * SS[:, None] * e
e1 = np.abs(Us - Us_exact).max()
e2 = np.abs(Ts - Ts_exact).max()
check("all 3 slave translations exact (incl. the in-plane projector part)",
      e1 < 1e-14, f"err={e1:.2e}")
check("all 3 slave rotations exact (incl. the non-slope projector part)",
      e2 < 1e-14, f"err={e2:.2e}")

# --- T6: honest limit — Mindlin shear breaks slope=rotation at O(γ·h) -----------
print("\nT6  honest limit: Mindlin data (φ = dw/ds − γ) ⇒ O(γ·h) w-row error")


def shear_err(n, gamma):
    # slaves at ξ=0.3 of each segment — NOT the midpoint: the shear error enters the w
    # row through (H2+H4) ∝ ξ−3ξ²+2ξ³, which has a root exactly at ξ=0.5, so midpoint
    # slaves would (spuriously) hide the Kirchhoff-assumption error entirely
    ms = np.linspace(0.0, 1.0, n + 1)
    ss = ms[:-1] + 0.3 * (ms[1] - ms[0])
    hw, hp, _, _ = hermite_rows(ms, ss)
    w_m, p_m = 0.5 * kappa * ms**2, kappa * ms - gamma   # rotation lags the slope by γ
    w_s = 0.5 * kappa * ss**2
    return np.abs(hw @ w_m + hp @ p_m - w_s).max()


g1, g2 = shear_err(2, 1e-3), shear_err(2, 2e-3)          # linear in γ
h1, h2 = shear_err(2, 1e-3), shear_err(4, 1e-3)          # linear in h
check("w-row error is linear in the shear angle γ", abs(g2 / g1 - 2.0) < 0.05,
      f"γ-ratio≈{g2/g1:.3f} (→2)")
check("w-row error is ~linear in h (refinement-safe)", 1.5 < h1 / h2 < 3.0,
      f"h-ratio≈{h1/h2:.2f}")
bound = 1e-3 * 0.5 * 0.0962          # analytic: γ·h·max_ξ|ξ−3ξ²+2ξ³| (max ≈ 0.0962)
check("error obeys the analytic bound γ·h·0.0962 and is ≪ the w-scale",
      g1 <= bound * 1.01 and g1 < 1e-2 * 0.5 * kappa,
      f"err={g1:.2e} ≤ bound={bound:.2e}, w-scale={0.5*kappa:.2e}")
print("      → the Hermite row assumes Kirchhoff (slope=rotation); shear-deformable")
print("        states carry an O(γ·h) tie error — vanishing for thin shells AND on refinement.")

# ===================================================================== summary
print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'} — "
      "P3.1 Hermite w–θ transfer: rigid modes + quadratic/cubic w cross EXACTLY "
      "(linear P is O(h²) on the same field), O(h⁴) on quartics, orientation-invariant "
      "at the global-DOF level, Kirchhoff-limit bounded at O(γ·h).")
sys.exit(0 if _fails == 0 else 1)
