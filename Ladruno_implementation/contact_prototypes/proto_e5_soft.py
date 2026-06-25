"""ADR-57 E5 ORACLE — edge-edge explicit Courant-stable SOFT penalty.

Oracle-first (the fork discipline): validate the SOFT=edge enforcement in pure numpy BEFORE
transcribing it into the LadrunoContactFE EDGE_EDGE mode. No OpenSees, no build. E5 is the
edge-edge analogue of the shipped B1 (NTS) / B2 (mortar) SOFT lane: under the explicit
CentralDifferenceLadruno the per-pair penalty εN is REPLACED by the Courant-stable
k_soft = SOFSCL·4·m_eff/dt² so the edge contact's own central-difference mode ω·dt = 2√SOFSCL ≤ 2
NEVER throttles dt_cr — edge-on impact runs at the STRUCTURAL dt (the progressive-collapse enabler,
now for the cos_t→0 case). Inert under implicit / SOFT-absent ⇒ byte-identical (the configured εN).

What E5 owns (the ADR-57 §6 SOFT tier + the E5 gate row):
  - the 4-node edge gap-mode generalized mass m_eff = 1/(B M⁻¹ Bᵀ), B = [(1−s)n, s n, −(1−t)n, −t n]
    the E2 B-operator. For a diagonal lumped mass this is the CLOSED FORM
        m_eff = 1/( (1−s)²·invMproj_a0 + s²·invMproj_a1 + (1−t)²·invMproj_b0 + t²·invMproj_b1 )
    invMproj = Σ_d n_d²/m_d the nodal mass projected on the gap direction n; a FIXED/massless node
    ⇒ ∞ mass ⇒ 0 contribution (the B1 ladrunoInvMassProj rule, generalized n→4-node edge operator).
  - k_soft = SOFSCL·4·m_eff/dt² ⇒ the edge mode ω·dt = √(k_soft/m_eff)·dt = 2√SOFSCL ≤ 2 (the
    central-difference stability bound), INDEPENDENT of dt — so the contact never sets dt_cr.
  - STIFF (configured εN) DIVERGES vs SOFT bounded: a central-difference impact at the structural dt
    with the configured εN has ω·dt ≫ 2 and blows up; the SAME dt with k_soft stays bounded (energy
    restitution — the rebound speed ≈ the impact speed, no spurious injection).
  - friction softKt follows the B1-kt n→t rule: the tangential STICK penalty
    k_soft_t = SOFSCL·4·m_eff_t/dt² with m_eff_t the WORST-CASE (smallest) edge gap-mode mass over
    the two basis tangents t1,t2 to n, so a stiff friction kt never throttles dt_cr via ω_t.
  - SOFSCL≤0 / implicit (no CDL) ⇒ the configured εN is returned ⇒ byte-identical to E2.

Run: <pythoncore-3.12> proto_e5_soft.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9
_fails = 0

OK, PARALLEL, DEGENERATE = 0, -1, -2
TAU_PAR = np.sin(np.deg2rad(15.0)) ** 2
MARGIN = 1e-6


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ============================================ E0 kernel (compact, mirrors LadrunoEdgeKernel)
def closest_pt_segseg(p1, q1, p2, q2, tau_par=TAU_PAR):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    out = {"status": OK, "s": 0.0, "t": 0.0, "c1": p1.copy(), "c2": p2.copy(), "interior": False}
    if a <= 1e-9 or e <= 1e-9:
        out["status"] = DEGENERATE
        return out
    c, b, f = d1 @ r, d1 @ d2, d2 @ r
    denom = a * e - b * b
    s = np.clip((b * f - c * e) / denom, 0.0, 1.0) if denom > 1e-300 else 0.0
    t = (b * s + f) / e
    if t < 0.0:
        t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
    elif t > 1.0:
        t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)
    out["s"], out["t"] = s, t
    out["c1"], out["c2"] = p1 + s * d1, p2 + t * d2
    out["interior"] = (MARGIN <= s <= 1 - MARGIN) and (MARGIN <= t <= 1 - MARGIN)
    if denom < tau_par * a * e:
        out["status"] = PARALLEL
    return out


def edge_normal(d1, d2, sign_ref):
    m = np.cross(d1, d2)
    n = m / np.linalg.norm(m)
    return n if (n @ sign_ref) >= 0 else -n


# ============================================ the E5 SOFT layer (mirror the C++ to-be)
def inv_mass_proj(m, n):
    """invMproj = Σ_d n_d²/m_d over the 3 translational DOFs; a DOF with m_d ≤ 0 ⇒ ∞ mass ⇒ 0
    (the B1 ladrunoInvMassProj rule — a fixed/massless node drops out)."""
    s = 0.0
    for d in range(3):
        if m[d] > 0.0:
            s += n[d] * n[d] / m[d]
    return s


def edge_gap_mode_inv_mass(masses, dir_, s, t):
    """B M⁻¹ Bᵀ for the 4-node edge operator B = [(1−s)dir, s dir, −(1−t)dir, −t dir] and diagonal
    nodal masses. The signs square away ⇒ the closed form Σ w_i²·invMproj_i, w=[(1−s),s,(1−t),t]."""
    w = np.array([1 - s, s, 1 - t, t])
    return sum(w[i] ** 2 * inv_mass_proj(masses[i], dir_) for i in range(4))


def edge_gap_mode_inv_mass_matrix(masses, dir_, s, t):
    """the SAME quantity by the FULL B M⁻¹ Bᵀ matrix product (independent reference for the closed
    form). B is 1×12 over [a0,a1,b0,b1]; M⁻¹ is the 12×12 diagonal inverse nodal mass."""
    w = np.array([1 - s, s, -(1 - t), -t])           # the signed edge weights
    B = np.concatenate([wi * dir_ for wi in w])      # 1×12
    Minv = np.zeros(12)
    for i in range(4):
        for d in range(3):
            Minv[3 * i + d] = 1.0 / masses[i][d] if masses[i][d] > 0.0 else 0.0
    return float(B @ (Minv * B))


def soft_kn(masses, n, s, t, dt, sofscl):
    """k_soft = SOFSCL·4·m_eff/dt², m_eff = 1/(B M⁻¹ Bᵀ). ≤0 invMass / dt ⇒ caller falls back to εN."""
    inv = edge_gap_mode_inv_mass(masses, n, s, t)
    if dt <= 0.0 or inv <= 0.0:
        return None
    return sofscl * 4.0 * (1.0 / inv) / (dt * dt)


def basis_tangents(n):
    """two orthonormal tangents to n (the softKt construction: t1 from the least-aligned axis)."""
    k = int(np.argmin(np.abs(n)))
    e = np.zeros(3); e[k] = 1.0
    t1 = e - (e @ n) * n
    t1 = t1 / np.linalg.norm(t1)
    t2 = np.cross(n, t1)
    return t1, t2


def soft_kt(masses, n, s, t, dt, sofscl):
    """k_soft_t = SOFSCL·4·m_eff_t/dt², m_eff_t = the WORST-CASE (smallest) edge gap-mode mass over
    the two basis tangents to n (largest inverse mass ⇒ the binding stick-mode stability bound)."""
    t1, t2 = basis_tangents(n)
    inv1 = edge_gap_mode_inv_mass(masses, t1, s, t)
    inv2 = edge_gap_mode_inv_mass(masses, t2, s, t)
    inv = max(inv1, inv2)
    if dt <= 0.0 or inv <= 0.0:
        return None
    return sofscl * 4.0 * (1.0 / inv) / (dt * dt)


# ============================================================================ tests
print("ADR-57 E5 oracle — edge-edge explicit Courant-stable SOFT penalty")

# canonical crossed bars: slave edge ∥ x at z=zc, master edge ∥ y at z=0, n≈+z.
def bars(zc=-0.02):
    return (np.array([-1, 0, zc]), np.array([1, 0, zc]),
            np.array([0, -1, 0.]), np.array([0, 1, 0.]))
ORIENT = np.array([0, 0, 1.])

# --- T1: closed-form m_eff == the full B M⁻¹ Bᵀ matrix product (the generalization is exact) ---
print("\nT1  m_eff = 1/(BM⁻¹Bᵀ): closed form == full matrix product, anisotropic + off-center")
p1, q1, p2, q2 = bars(-0.02)
k = closest_pt_segseg(p1, q1, p2, q2)
n = edge_normal(q1 - p1, q2 - p2, ORIENT)
s, t = k["s"], k["t"]
# anisotropic, distinct nodal masses on all 4 nodes (master nodes NOT fixed here)
masses = [np.array([2.0, 3.0, 1.5]), np.array([1.0, 1.0, 4.0]),
          np.array([2.5, 0.5, 2.0]), np.array([1.5, 2.0, 0.8])]
inv_closed = edge_gap_mode_inv_mass(masses, n, s, t)
inv_matrix = edge_gap_mode_inv_mass_matrix(masses, n, s, t)
check("closed form == matrix product (B M⁻¹ Bᵀ)", abs(inv_closed - inv_matrix) < 1e-12,
      f"closed={inv_closed:.6e} matrix={inv_matrix:.6e}")
# off-center crossing (s,t away from 0.5) still matches
ps1, qs1 = np.array([-1, 0.3, -0.02]), np.array([1, -0.2, -0.02])
pm2, qm2 = np.array([0.4, -1, 0.0]), np.array([-0.3, 1, 0.0])
ko = closest_pt_segseg(ps1, qs1, pm2, qm2)
no = edge_normal(qs1 - ps1, qm2 - pm2, ORIENT)
ic = edge_gap_mode_inv_mass(masses, no, ko["s"], ko["t"])
im = edge_gap_mode_inv_mass_matrix(masses, no, ko["s"], ko["t"])
check("off-center (s,t≠0.5): closed == matrix", abs(ic - im) < 1e-12,
      f"s={ko['s']:.3f} t={ko['t']:.3f}")

# --- T2: a FIXED (massless) master node drops out (∞ mass) — m_eff → the slave-only mass ---
print("\nT2  fixed/massless master node ⇒ ∞ mass ⇒ 0 contribution (the B1 invMassProj rule)")
m_fixed = [np.array([1.0, 1.0, 1.0]), np.array([1.0, 1.0, 1.0]),
           np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])]   # master fixed
inv_fixed = edge_gap_mode_inv_mass(m_fixed, n, s, t)
inv_slave_only = (1 - s) ** 2 * inv_mass_proj(m_fixed[0], n) + s ** 2 * inv_mass_proj(m_fixed[1], n)
check("fixed master drops out (m_eff = slave-only)", abs(inv_fixed - inv_slave_only) < 1e-12,
      f"inv_fixed={inv_fixed:.6e}")
# with the symmetric crossing s=t=0.5, equal unit slave masses: invMproj on n=ẑ is 1/1=1 each ⇒
# inv = 0.25·1 + 0.25·1 = 0.5 ⇒ m_eff = 2.0 (the harmonic half-mass of the two slave nodes)
check("symmetric unit slave m_eff = 2.0 (1/(0.25+0.25))", abs(1.0 / inv_fixed - 2.0) < 1e-9,
      f"m_eff={1.0/inv_fixed:.4f}")

# --- T3: ω·dt = 2√SOFSCL (the Courant identity) — INDEPENDENT of dt and εN ---
print("\nT3  k_soft = SOFSCL·4·m_eff/dt² ⇒ ω·dt = 2√SOFSCL ≤ 2 (independent of dt)")
masses_u = [np.array([1.0, 1.0, 1.0])] * 2 + [np.array([0.0, 0.0, 0.0])] * 2   # fixed master
inv = edge_gap_mode_inv_mass(masses_u, n, s, t)
m_eff = 1.0 / inv
for sofscl in (0.05, 0.1, 0.25, 1.0):
    for dt in (1e-3, 1e-4, 5e-5):
        ks = soft_kn(masses_u, n, s, t, dt, sofscl)
        omega_dt = np.sqrt(ks / m_eff) * dt
        check(f"SOFSCL={sofscl:<5} dt={dt:.0e}: ω·dt = 2√SOFSCL = {2*np.sqrt(sofscl):.4f}",
              abs(omega_dt - 2 * np.sqrt(sofscl)) < 1e-9, f"ω·dt={omega_dt:.6f}")

# --- T4: STIFF εN DIVERGES vs SOFT bounded — a central-difference impact at the structural dt ---
print("\nT4  central-difference impact: stiff εN DIVERGES, SOFT bounded + energy-restituting")
# 1-DOF gap-mode oscillator: x'' = -k/m_eff·⟨−x⟩ (one-sided penalty contact), unit slave mass.
# structural dt sized from the BODY (not the contact): here just a fixed dt the contact must not throttle.
def central_diff_contact(k, m_eff, dt, x0, v0, nsteps):
    """one-sided penalty contact: accel = +(k/m_eff)·⟨−x⟩ pushes OUT of penetration (x<0). The
    penalty traction tN = εN⟨−x⟩ acts toward +x — the restoring sign. Returns (max_penetration,
    v_out) where v_out is the post-rebound separation speed (energy-restitution probe)."""
    x = x0
    a0 = (k / m_eff) * max(0.0, -x0)              # restoring accel at the seed (0 outside contact)
    xm = x0 - v0 * dt + 0.5 * a0 * dt * dt        # x_{-1} seed
    pen = 0.0
    for _ in range(nsteps):
        a = (k / m_eff) * max(0.0, -x)            # penalty accel (Macaulay one-sided, OUTWARD)
        xp = 2 * x - xm + a * dt * dt             # central difference
        xm, x = x, xp
        if -x > pen:
            pen = -x                              # deepest penetration reached
        if not np.isfinite(x):
            return np.inf, np.inf
    v_out = (x - xm) / dt                          # final velocity (rebound, body separated)
    return pen, v_out
dt = 1e-3
m_eff = 2.0
EPS_STIFF = 1.0e8                                  # the configured (stiff) penalty
k_soft = soft_kn(masses_u, n, s, t, dt, 0.1)       # SOFSCL=0.1
# impact: start just outside contact moving in (x0>0, v0<0). The STRUCTURAL dt (1e-3) is set by the
# body, NOT the contact — SOFT must run stably at it; the stiff εN must NOT.
x0, v0, ns = 0.01, -1.0, 300
pen_stiff, vout_stiff = central_diff_contact(EPS_STIFF, m_eff, dt, x0, v0, ns)
pen_soft, vout_soft = central_diff_contact(k_soft, m_eff, dt, x0, v0, ns)
omega_dt_stiff = np.sqrt(EPS_STIFF / m_eff) * dt
omega_dt_soft = np.sqrt(k_soft / m_eff) * dt
# analytic SOFT penetration: d = |v0|·√(m_eff/k_soft) (energy balance ½m v² = ½k d²)
d_analytic = abs(v0) * np.sqrt(m_eff / k_soft)
check("stiff εN: ω·dt ≫ 2 (unstable) ⇒ the integration BLOWS UP (energy injected)",
      omega_dt_stiff > 2.0 and (not np.isfinite(vout_stiff) or abs(vout_stiff) > 2.0 * abs(v0)),
      f"ω·dt={omega_dt_stiff:.1f}, pen={pen_stiff:.2e}, v_out={vout_stiff:.2e}")
check("SOFT: ω·dt = 2√0.1 ≈ 0.632 ≤ 2 (stable) ⇒ penetration bounded ≈ analytic",
      omega_dt_soft < 2.0 and abs(pen_soft - d_analytic) < 0.05 * d_analytic,
      f"ω·dt={omega_dt_soft:.4f}, pen={pen_soft:.4e} vs analytic {d_analytic:.4e}")
check("SOFT energy-restituting: rebound speed ≈ impact speed (|v_out|≈|v0|, no injection/dissipation)",
      vout_soft > 0 and abs(abs(vout_soft) - abs(v0)) < 0.05 * abs(v0),
      f"v_out={vout_soft:.4f} vs v0={v0}")

# --- T5: softKt (the B1-kt n→t rule) — worst-case over the two tangents; isotropic ⇒ == softKn mass ---
print("\nT5  softKt: k_soft_t = SOFSCL·4·m_eff_t/dt², m_eff_t worst-case over t1,t2 ⊥ n")
# ISOTROPIC lumped mass ⇒ invMproj is direction-independent ⇒ m_eff_t == m_eff_n exactly
kt_soft = soft_kt(masses_u, n, s, t, dt, 0.1)
check("isotropic mass: softKt == softKn (invMproj direction-independent)", abs(kt_soft - k_soft) < 1e-6,
      f"softKt={kt_soft:.4e} softKn={k_soft:.4e}")
omega_t_dt = np.sqrt(kt_soft / m_eff) * dt
check("stick mode ω_t·dt = 2√SOFSCL ≤ 2 (kt never throttles dt_cr)", abs(omega_t_dt - 2 * np.sqrt(0.1)) < 1e-9,
      f"ω_t·dt={omega_t_dt:.4f}")
# anisotropic mass ⇒ softKt picks the WORST (smallest m_eff_t = largest inverse mass) ⇒ conservative
t1, t2 = basis_tangents(n)
i1 = edge_gap_mode_inv_mass(masses, t1, s, t)
i2 = edge_gap_mode_inv_mass(masses, t2, s, t)
kt_aniso = soft_kt(masses, n, s, t, dt, 0.1)
kt_from_worst = 0.1 * 4.0 * (1.0 / max(i1, i2)) / (dt * dt)
check("anisotropic: softKt uses the WORST-CASE (max inverse mass) tangent",
      abs(kt_aniso - kt_from_worst) < 1e-6, f"softKt={kt_aniso:.4e}")

# --- T6: SOFSCL≤0 / no dt (implicit) ⇒ fall back to the configured εN (byte-identical contract) ---
print("\nT6  SOFSCL≤0 / no CDL dt ⇒ the configured εN (byte-identical to E2)")
# the C++ contract: softScale<=0 OR a failed CDL dynamic_cast ⇒ return kn unchanged. Mirror the GATE:
def soft_kn_gated(masses, n, s, t, dt, sofscl, kn, is_cdl):
    if sofscl <= 0.0 or not is_cdl:
        return kn                                  # off / implicit ⇒ configured εN (byte-identical)
    ks = soft_kn(masses, n, s, t, dt, sofscl)
    return kn if ks is None else ks
KN = 1.0e5
check("SOFSCL=0 ⇒ kn (byte-identical)", soft_kn_gated(masses_u, n, s, t, dt, 0.0, KN, True) == KN)
check("implicit (not CDL) ⇒ kn even with SOFSCL>0", soft_kn_gated(masses_u, n, s, t, dt, 0.1, KN, False) == KN)
check("massless slave / dt≤0 ⇒ kn fallback", soft_kn_gated([np.zeros(3)] * 4, n, s, t, dt, 0.1, KN, True) == KN)
check("SOFT active (CDL, mass, dt) ⇒ k_soft ≠ kn", soft_kn_gated(masses_u, n, s, t, dt, 0.1, KN, True) != KN)

print(f"\n{'='*66}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E5 edge-edge explicit SOFT oracle)\n{'='*66}")
sys.exit(1 if _fails else 0)
