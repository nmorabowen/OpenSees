"""B2 follow-up ORACLE — frictional SOFT=2 (Courant-stable explicit Coulomb/Tresca on the
segment-based soft penalty).

B2 (#406) ships the SOFT=2 segment-based NORMAL penalty (oracle proto_b2_soft2_segment.py): per
slave node I the area-normalised gap ḡ_I = g̃_I/a_I is restrained by a Courant-stable penalty
k_soft,I = SOFSCL·4·m_eff,I/dt² so explicit dt_cr is never throttled. This adds friction on the SAME
segment overlap. Two pieces are NEW (everything else reuses the shipped LadrunoFrictionKernel +
LadrunoMortarKernel, already oracle-validated):

  (A) the TANGENTIAL gap-mode mass m_eff_t,I — the B2 normal m_eff with the gap normal n replaced by a
      surface tangent t. The tangential slip ḡT_I = B_t,I·u is linear in the facet displacements with
          B_t,I : slave J → (D_IJ/a_I) t ,  master K → −(M_IK/a_I) t ,
      so m_eff_t,I = 1/(B_t,I M⁻¹ B_t,Iᵀ) = 1/(Σ_J (D_IJ/a_I)² invMproj(m_sJ,t) + Σ_K (M_IK/a_I)²
      invMproj(m_mK,t)). The Courant-stable STICK penalty is k_soft_t,I = SOFSCL·4·m_eff_t,I/dt² (the
      B1-kt rule generalized from the NTS tangent operator to the segment operator), worst-cased over
      the two basis tangents t1,t2 — so a stiff friction kt never throttles dt_cr via the stick mode
      ω_t = √(kt/m_eff_t) (continuous stick at a high cone diverges when ω_t·dt > 2).
  (B) the cone over the SOFT normal pressure — N_I = |p_I| = |k_soft,I·ḡ_I| is the soft penalty (not a
      configured epsN), feeding the SAME return map cap = min(μ N_I + c, τmax). Penalty friction
      (λ_T ≡ 0; explicit single-pass, like the SOFT=2 normal).

GATES:
  T1 — m_eff_t closed form == assembled B_t M⁻¹ B_t ᵀ; isotropic lumped mass ⇒ m_eff_t == m_eff_n
       (direction-independent); fixed-master ⇒ slave-side mass.
  T2 — Courant STICK stability: k_soft_t = SOFSCL·4·m_eff_t/dt² ⇒ ω_t·dt = 2√SOFSCL ≤ 2; a stiff kt
       (ω_t·dt ≫ 2) is unstable at the same dt.
  T3 — the cone over the soft pressure: stick (|tT^tr| ≤ cap) keeps the trial; slip returns to cap,
       direction opposes the slip; Coulomb/Tresca/cohesion via min(μN+c, τmax); μ=0 ⇒ zero friction.
  T4 — D/M scatter self-equilibrium (Σφ=1): Σ_K f^s_K + Σ_L f^m_L = 0 for a constant nodal traction.

Run: <pythoncore-3.12> proto_soft2_friction.py   (numpy only)
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


def invMproj(m, d):
    """inverse mass projected on unit direction d: Σ_k d_k²/m_k; a DOF with no mass ⇒ ∞ mass ⇒ 0."""
    return sum((d[k] * d[k] / m[k]) for k in range(3) if m[k] > 0.0)


def m_eff_seg(D_row, M_row, aI, ms, mm, d):
    """segment-based gap-mode generalized mass for node I in direction d: 1/(B M⁻¹ Bᵀ),
    B : slave J → (D_IJ/aI) d, master K → −(M_IK/aI) d."""
    s = 0.0
    for J in range(len(D_row)):
        c = D_row[J] / aI
        s += c * c * invMproj(ms[J], d)
    for K in range(len(M_row)):
        c = M_row[K] / aI
        s += c * c * invMproj(mm[K], d)
    return 1.0 / s if s > 0.0 else np.inf


print("B2 follow-up — frictional SOFT=2 oracle\n")

# ---- a matched unit facet pair: D = (area/4)·I on a bilinear quad's lumped rows, M = D (matched) ----
# Use the lumped (row-summed) mortar matrices for a matched [0,1]² facet: a_I = 1/4 per node.
area = 1.0
aI = area / 4.0
D_row = np.array([aI, 0.0, 0.0, 0.0])          # lumped (diagonal) D row for node 0
M_row = np.array([aI, 0.0, 0.0, 0.0])          # matched ⇒ M row == D row (master node 0 carries it)
ms = [np.array([2.0, 2.0, 2.0])] * 4            # isotropic slave nodal mass
mm = [np.array([3.0, 3.0, 3.0])] * 4            # isotropic master nodal mass

n = np.array([0.0, 0.0, 1.0])
t1 = np.array([1.0, 0.0, 0.0])
t2 = np.array([0.0, 1.0, 0.0])

# -- T1: m_eff_t closed form == assembled; isotropic ⇒ direction-independent; fixed-master limit ----
# assembled B_t M⁻¹ B_tᵀ for node 0 (only node-0 entries nonzero in the lumped rows)
def assembled_invmass(D_row, M_row, aI, ms, mm, d):
    inv = 0.0
    for J in range(4):
        c = D_row[J] / aI
        inv += c * c * invMproj(ms[J], d)
    for K in range(4):
        c = M_row[K] / aI
        inv += c * c * invMproj(mm[K], d)
    return inv

me_t1 = m_eff_seg(D_row, M_row, aI, ms, mm, t1)
me_t2 = m_eff_seg(D_row, M_row, aI, ms, mm, t2)
me_n = m_eff_seg(D_row, M_row, aI, ms, mm, n)
inv_assem = assembled_invmass(D_row, M_row, aI, ms, mm, t1)
check("T1 m_eff_t == 1/(assembled B_t M⁻¹ B_tᵀ)", abs(me_t1 - 1.0 / inv_assem) < 1e-12,
      f"m_eff_t1={me_t1:.6f}")
check("T1 isotropic lumped mass ⇒ m_eff_t == m_eff_n (direction-independent)",
      abs(me_t1 - me_n) < 1e-12 and abs(me_t2 - me_n) < 1e-12, f"t1={me_t1:.4f}, n={me_n:.4f}")
# matched lumped: m_eff = 1/(1/m_s + 1/m_m) (harmonic) = 1/(1/2 + 1/3) = 1.2
check("T1 matched harmonic m_eff = 1/(1/m_s + 1/m_m)", abs(me_n - 1.0 / (1.0/2.0 + 1.0/3.0)) < 1e-12,
      f"m_eff={me_n:.4f}")
# fixed master (∞ mass ⇒ 0 contribution) ⇒ m_eff = m_s
mm_fixed = [np.array([0.0, 0.0, 0.0])] * 4
me_fixed = m_eff_seg(D_row, M_row, aI, ms, mm_fixed, t1)
check("T1 fixed-master limit ⇒ m_eff_t == m_slave", abs(me_fixed - 2.0) < 1e-12, f"m_eff={me_fixed:.4f}")

# -- T2: Courant STICK stability k_soft_t = SOFSCL·4·m_eff_t/dt² ⇒ ω_t·dt = 2√SOFSCL ----------------
SOFSCL, dt = 0.10, 1.0e-3
k_soft_t = SOFSCL * 4.0 * me_t1 / (dt * dt)
omega_t = np.sqrt(k_soft_t / me_t1)
check("T2 ω_t·dt = 2√SOFSCL ≤ 2 (Courant-stable stick)", abs(omega_t * dt - 2.0 * np.sqrt(SOFSCL)) < 1e-9,
      f"ω_t·dt={omega_t*dt:.4f}")
# central-difference stick oscillator: stable when ω·dt < 2, diverges when > 2
def cd_stick_stable(k, m, dt, nsteps=4000):
    x, xp = 1.0e-3, 1.0e-3
    for _ in range(nsteps):
        xn = 2.0 * x - xp - (dt * dt) * (k / m) * x
        xp, x = x, xn
        if abs(x) > 1e3:
            return False
    return True
check("T2 soft stick stable; stiff kt (ω_t·dt≫2) diverges at the same dt",
      cd_stick_stable(k_soft_t, me_t1, dt) and not cd_stick_stable(1.0e9, me_t1, dt))

# -- T3: the cone over the SOFT normal pressure (the shipped return map, re-shown) ------------------
def return_map(gTeff, gpT, N, kt, mu, coh=0.0, tauMax=-1.0):
    tT_tr = kt * (gTeff - gpT)            # trial tangential traction (1-D here)
    cap = mu * N + coh
    if tauMax >= 0.0:
        cap = min(cap, tauMax)
    if abs(tT_tr) <= cap:                 # STICK
        return tT_tr, gpT
    tT = cap * np.sign(tT_tr)             # SLIP — return to the cone
    gpT_new = gTeff - tT / kt
    return tT, gpT_new

N_I = 5.0                                 # soft normal pressure magnitude |p_I|
kt = k_soft_t
# stick: small slip stays elastic
tT_s, _ = return_map(0.1 * (0.3*N_I)/kt, 0.0, N_I, kt, mu=0.3)
check("T3 stick keeps the elastic trial (|tT| < μN)", abs(tT_s) < 0.3 * N_I + 1e-12)
# slip: large slip caps at μN, opposes motion
tT_l, _ = return_map(10.0, 0.0, N_I, kt, mu=0.3)
check("T3 slip caps at μN and opposes motion", abs(tT_l - 0.3 * N_I) < 1e-9 and tT_l > 0.0)
# Tresca cap dominates when τmax < μN
tT_tr, _ = return_map(10.0, 0.0, N_I, kt, mu=0.3, tauMax=1.0)
check("T3 Tresca cap min(μN, τmax) dominates", abs(tT_tr - 1.0) < 1e-9)
# cohesion adds an intercept
tT_c, _ = return_map(10.0, 0.0, 0.0, kt, mu=0.3, coh=0.7)
check("T3 cohesion gives a no-pressure shear strength c", abs(tT_c - 0.7) < 1e-9)
# μ=0, c=0 ⇒ no friction
tT_0, _ = return_map(10.0, 0.0, N_I, kt, mu=0.0)
check("T3 μ=0, c=0 ⇒ zero friction (byte-identical frictionless SOFT=2)", abs(tT_0) < 1e-12)

# -- T4: D/M scatter self-equilibrium (Σφ=1 ⇒ Σ D = Σ M = a_I) -------------------------------------
# constant nodal tangential traction tF over the matched facet: f^s_K = Σ_I D_KI tF_I, f^m_L = −Σ_I M_IL tF_I
D = np.diag([aI, aI, aI, aI])
M = np.diag([aI, aI, aI, aI])             # matched
tF = np.array([0.4, -0.2, 0.1, 0.3])      # arbitrary per-node tangential traction (1 component)
fs = D @ tF
fm = -(M.T @ tF)
check("T4 D/M scatter self-equilibrates (Σf^s + Σf^m = 0)", abs(fs.sum() + fm.sum()) < 1e-12,
      f"Σ={fs.sum()+fm.sum():.2e}")

print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAIL(S)'}  (frictional SOFT=2)")
sys.exit(1 if _fails else 0)
