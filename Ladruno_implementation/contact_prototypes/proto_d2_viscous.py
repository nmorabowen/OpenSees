"""ADR-41 D2 ORACLE — viscous contact stabilization (p_visc = μ_c·v_rel).

Oracle-first (the fork discipline): validate the velocity-proportional NORMAL contact damper that
suppresses contact chatter / snap-through in the pounding/rocking/uplift regime, BEFORE wiring it
into the C++ adapter. It reuses the shipped normal-penalty operator driven by VELOCITY:

    gap rate ġ = B·v        (B = ∂g/∂u : RIGID_PLANE B=n ; SEGMENT B=[n | −Nᵢ n])
    viscous force  f_visc = −μ_c·(B·v)·B = −μ_c·ġ·B      (into the residual, from getTrialVel)
    viscous damping matrix C_visc = μ_c·B Bᵀ             (into C, via addCtoTang)

active ONLY while in contact (gap<0); v≡0 in statics ⇒ identically inert (static byte-identical).

The mechanics are a Kelvin–Voigt (linear spring + dashpot) contact, so a normal impact has the
EXACT coefficient of restitution

    e = exp(−ζπ / √(1−ζ²)),   ζ = μ_c / (2·√(kₙ·m))      (half-period damped decay over contact)

— a strong analytic gate the numerical integration must reproduce.

Gates:
  T1  viscous operator: C_visc = μ_c·B Bᵀ == −∂f_visc/∂v (FD), SPD; f_visc opposes the gap rate.
  T2  static-inert: v ≡ 0 ⇒ f_visc ≡ 0, C_visc ≠ 0 but contributes no force (byte-identical statics).
  T3  HEADLINE impact restitution: μ_c=0 ⇒ e≈1 (elastic penalty); μ_c>0 ⇒ e = exp(−ζπ/√(1−ζ²))
      (numerics vs analytic, monotone decreasing in μ_c) — chatter/snap-through energy bleed.
  T4  active mask: a SEPARATED pair (gap≥0) ⇒ zero viscous force (never damps a free-flying node).
  T5  energy balance: the KE lost across the impact == the dashpot dissipation ∫ μ_c·ġ² dt.

Run: <pythoncore-3.12> proto_d2_viscous.py   (numpy only)
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


print("ADR-41 D2 oracle — viscous contact stabilization (p_visc = μ_c·v_rel)")

# ===================================== T1: the viscous operator (SEGMENT-like B = [n | −Nᵢ n])
print("\nT1  viscous operator: C_visc = μ_c·B Bᵀ == −∂f_visc/∂v (FD), SPD; force opposes gap rate")
n = np.array([0.0, 0.0, 1.0])                      # contact normal (+z)
N = np.array([0.4, 0.35, 0.25])                    # master shape weights at the projection (Σ=1)
mu_c = 50.0
# B (the gap-gradient row) over [slave xyz | m0 xyz | m1 xyz | m2 xyz] = [n | −N0 n | −N1 n | −N2 n]
ndof = 3 * (1 + 3)
B = np.zeros(ndof)
B[0:3] = n
for i in range(3):
    B[3 * (1 + i):3 * (1 + i) + 3] = -N[i] * n


def f_visc(v):                                     # f_visc = −μ_c·(B·v)·B  (the residual term)
    return -mu_c * (B @ v) * B


C_visc = mu_c * np.outer(B, B)                     # the damping matrix
# C_visc == −∂f_visc/∂v (central FD)
h = 1e-6
Cfd = np.zeros((ndof, ndof))
for j in range(ndof):
    vp = np.zeros(ndof); vp[j] = h
    vm = np.zeros(ndof); vm[j] = -h
    Cfd[:, j] = -(f_visc(vp) - f_visc(vm)) / (2 * h)
rel = np.abs(C_visc - Cfd).max() / (np.abs(Cfd).max() + 1e-30)
check("C_visc == −∂f_visc/∂v (FD ≤1e-6)", rel <= 1e-6, f"rel={rel:.2e}")
check("C_visc is symmetric", np.allclose(C_visc, C_visc.T, atol=1e-12))
evals = np.linalg.eigvalsh(C_visc)
check("C_visc is PSD (rank-1 dashpot ⇒ one +eig, rest 0)", evals.min() > -1e-9 and evals.max() > 0,
      f"λ=[{evals.min():.2e}, {evals.max():.3f}]")
# self-equilibrium: a rigid approach velocity (slave & master move together) ⇒ zero gap rate ⇒ no force
v_rigid = np.tile([0.1, -0.2, 0.3], 1 + 3)
check("rigid-body velocity ⇒ ġ=0 ⇒ f_visc=0 (self-equilibrating, frame-invariant)",
      np.abs(f_visc(v_rigid)).max() < 1e-12, f"max|f|={np.abs(f_visc(v_rigid)).max():.2e}")
# sign: slave approaching the master (slave −z velocity) ⇒ ġ<0 ⇒ viscous force pushes slave +z (apart)
v_approach = np.zeros(ndof); v_approach[2] = -1.0          # slave moving down (into master)
f = f_visc(v_approach)
check("viscous force OPPOSES the approach (slave gets +z restoring force)", f[2] > 0,
      f"f_slave_z={f[2]:.3f}")

# ===================================== T2: static-inert
print("\nT2  static-inert: v ≡ 0 ⇒ f_visc ≡ 0 (statics byte-identical regardless of μ_c)")
check("v=0 ⇒ f_visc=0", np.abs(f_visc(np.zeros(ndof))).max() == 0.0)
check("C_visc ≠ 0 (the matrix is assembled, but multiplies v=0 ⇒ no static force)", np.abs(C_visc).max() > 0)

# ===================================== T3: HEADLINE — impact restitution vs the Kelvin–Voigt analytic
print("\nT3  HEADLINE: 1-DOF normal impact restitution e = exp(−ζπ/√(1−ζ²)) (numerics vs analytic)")
m = 1.0
kn = 1.0e4
v_in = -2.0                                        # initial downward velocity (approaching the plane)


def impact_restitution(muc, dt=None, x0=1.0e-3):
    """Integrate a point mass dropping onto a rigid plane at x=0 (contact when x<0) with a penalty
    spring kn + viscous dashpot μ_c (BOTH active only while x<0). Central difference. Returns the
    coefficient of restitution e = |v_out|/|v_in| (measured once the mass separates, x>0, v>0)."""
    # stable dt: well below the penalty period T=2π√(m/kn); halve for the damped case margin.
    T = 2 * np.pi * np.sqrt(m / kn)
    if dt is None:
        dt = T / 400.0
    x = x0
    v = v_in
    a = 0.0
    # leapfrog/central-difference: x_{n+1} = x_n + v_n dt + 0.5 a_n dt²; a from force at x_{n+1};
    # v_{n+1} = v_n + 0.5(a_n + a_{n+1}) dt  (velocity-Verlet — symplectic, exact energy for μ_c=0)
    def force(xx, vv):
        if xx < 0.0:                               # in contact: penalty + viscous (both gated on gap<0)
            return -kn * xx - muc * vv             # f = −kn·g − μ_c·ġ  (g=x, ġ=v)
        return 0.0
    in_contact_seen = False
    for _ in range(200000):
        a = force(x, v) / m
        x_new = x + v * dt + 0.5 * a * dt * dt
        # velocity-Verlet with a velocity-dependent force ⇒ one fixed-point pass on a_{n+1}
        a_new = force(x_new, v + a * dt) / m
        v_new = v + 0.5 * (a + a_new) * dt
        x, v = x_new, v_new
        if x < 0.0:
            in_contact_seen = True
        if in_contact_seen and x > 0.0 and v > 0.0:    # separated and moving away
            return abs(v) / abs(v_in)
    return None


def e_analytic(muc):
    zeta = muc / (2.0 * np.sqrt(kn * m))
    if zeta >= 1.0:
        return 0.0                                 # overdamped: no rebound
    return float(np.exp(-zeta * np.pi / np.sqrt(1.0 - zeta * zeta)))


# μ_c = 0 ⇒ elastic (e ≈ 1)
e0 = impact_restitution(0.0)
check("μ_c=0 ⇒ elastic penalty rebound e≈1 (energy conserved)", abs(e0 - 1.0) < 5e-3, f"e={e0:.5f}")
# μ_c > 0 ⇒ e matches the Kelvin–Voigt analytic restitution
es = []
# underdamped sweep (ζ<1 ⇒ the exp(−ζπ/√(1−ζ²)) formula applies; ζ≥1 is overdamped — the
# Kelvin–Voigt model has no clean restitution there and is not how the stabilizer is used).
for muc in [20.0, 50.0, 120.0, 150.0]:
    e_num = impact_restitution(muc)
    e_an = e_analytic(muc)
    es.append(e_num if e_num is not None else 0.0)
    ok = (e_num is not None) and abs(e_num - e_an) < 0.02
    print(f"      μ_c={muc:6.1f}  ζ={muc/(2*np.sqrt(kn*m)):.4f}  "
          f"e_num={e_num if e_num is None else f'{e_num:.4f}'}  e_analytic={e_an:.4f}")
    check(f"μ_c={muc}: e_num == exp(−ζπ/√(1−ζ²)) (≤0.02)", ok)
check("restitution e DECREASES monotonically with μ_c (more damping ⇒ less rebound)",
      all(es[i] > es[i + 1] for i in range(len(es) - 1)), f"e={[f'{x:.3f}' for x in es]}")

# ===================================== T4: active mask (separated ⇒ no viscous)
print("\nT4  active mask: a SEPARATED pair (gap≥0) ⇒ zero viscous force (never damps free flight)")
# the force() above returns 0 for x>=0 regardless of v — pin it directly.
def contact_force(x, v, muc):
    return (-kn * x - muc * v) if x < 0.0 else 0.0
check("gap>0 (separated), any velocity ⇒ contact force = 0 (no spurious drag in free flight)",
      contact_force(+1e-3, -5.0, 100.0) == 0.0 and contact_force(+1e-3, +5.0, 100.0) == 0.0)
check("gap<0 (in contact), approaching ⇒ viscous opposes (force has a −μ_c·v term)",
      contact_force(-1e-4, -5.0, 100.0) > contact_force(-1e-4, 0.0, 100.0))

# ===================================== T5: energy balance (KE loss == dashpot dissipation)
print("\nT5  energy balance: KE lost across the impact == dashpot dissipation ∫ μ_c·ġ² dt")
muc = 120.0
T = 2 * np.pi * np.sqrt(m / kn)
dt = T / 800.0
x, v = 1.0e-3, v_in
diss = 0.0
in_contact_seen = False
e_out = None
def force(xx, vv, mc):
    return (-kn * xx - mc * vv) if xx < 0.0 else 0.0
for _ in range(400000):
    a = force(x, v, muc) / m
    x_new = x + v * dt + 0.5 * a * dt * dt
    a_new = force(x_new, v + a * dt, muc) / m
    v_new = v + 0.5 * (a + a_new) * dt
    # dashpot power = μ_c·ġ² (only in contact); trapezoidal over the step
    if x < 0.0:
        diss += muc * v * v * dt                    # ∫ μ_c ġ² dt (ġ=v)
    x, v = x_new, v_new
    if x < 0.0:
        in_contact_seen = True
    if in_contact_seen and x > 0.0 and v > 0.0:
        e_out = abs(v) / abs(v_in)
        break
ke_in = 0.5 * m * v_in ** 2
ke_out = 0.5 * m * (e_out * abs(v_in)) ** 2
ke_loss = ke_in - ke_out
check("KE lost == dashpot dissipation ∫μ_c ġ² dt (≤2% — viscous removes exactly the dissipated energy)",
      abs(ke_loss - diss) / ke_in < 0.02, f"ΔKE={ke_loss:.4f} diss={diss:.4f} rel={abs(ke_loss-diss)/ke_in:.3%}")

print(f"\n{'='*70}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(D2 viscous-stabilization oracle)\n{'='*70}")
sys.exit(1 if _fails else 0)
