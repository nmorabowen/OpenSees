"""ADR-41 A1 ORACLE — LadrunoFrictionKernel: unified cone min(μN+c, τmax) return map
+ consistent tangent (Css/Csl), FD-checked.

A1 extracts the P3/P3.5 Coulomb friction into LadrunoFrictionKernel.h and GENERALIZES
the cone from cap=μN to the unified:

    cap = min( μ·N + c ,  τmax )        N = kn·<−gap>₊ (contact pressure), c = cohesion

The pressure-coupling block Csl = ∂tT/∂gN is the min()-selected coupling: ∂cap/∂N = μ on
the Coulomb branch, 0 when the τmax cap binds (pressure-independent) or separated. This
oracle:
  (1) BIT-FOR-BIT: with c=0, τmax≤0 the unified map/tangent == the shipped pure-Coulomb
      cap=μN forms (the extraction must not change existing behavior);
  (2) Css/Csl vs central-difference FD to 1e-6 on Coulomb / Coulomb+cohesion (non-sym,
      d_TN=μ) and τmax-capped / Tresca (SYMMETRIC, d_TN=0);
  (3) return-map values (stick = elastic, slip = capped) and the Tresca/μ=0 corner.

Mirrors proto_p35_implicit_tangent.py's assembled-slave-block FD with n FROZEN (∂n/∂u is
the deferred P2b-2c curvature block, excluded). Run: <pythoncore-3.12> proto_a1_friction.py
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


# ---------------------------------------------------- unified cone (mirror C++)
def friction_cap(N, mu, c, tmax):
    capC = (mu * N + c) if N > 0.0 else 0.0
    return tmax if (tmax > 0.0 and tmax < capC) else capC


def tract_3d(gT, gpT, n, N, kt, mu, c=0.0, tmax=0.0):
    """traction tT (global tangent-plane vector). Mirrors frictionReturnMap (returns tT,
    not the negated applied force)."""
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT, float))
    cap = friction_cap(N, mu, c, tmax)
    nrm = np.linalg.norm(tT_tr)
    if cap <= 0.0:
        return 0.0 * tT_tr                            # zero cone radius: FREE SLIP (review P3)
    if nrm <= cap:
        return tT_tr                                  # stick
    return cap * (tT_tr / nrm)                         # slip: capped at the cone


def slave_block_analytic(gT, gpT, n, N, kn, kt, mu, c=0.0, tmax=0.0):
    """K_ss = D_TT·P_t + d_TN⊗n  with the min()-selected d_TN (∂cap/∂N = μ on the Coulomb
    branch, 0 when τmax binds or separated)."""
    n = np.asarray(n, float)
    Pt = np.eye(3) - np.outer(n, n)
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT, float))
    capC = (mu * N + c) if N > 0.0 else 0.0
    capped = (tmax > 0.0 and tmax < capC)
    cap = tmax if capped else capC
    nrm = np.linalg.norm(tT_tr)
    if cap <= 0.0:                                    # zero cone: FREE SLIP, K_ss = 0 (review P3)
        return 0.0 * Pt
    if nrm <= cap:                                    # STICK
        return kt * Pt
    nh = tT_tr / nrm
    D_TT = (cap * kt / nrm) * (Pt - np.outer(nh, nh))
    dCap_dN = 0.0 if (capped or N <= 0.0) else mu
    d_TN = -dCap_dN * kn * nh
    return D_TT @ Pt + np.outer(d_TN, n)


def fd_slave_block(X_s, gpT, n, p0, kn, kt, mu, c=0.0, tmax=0.0, h=1e-7):
    """CENTRAL-difference K_ss = −∂(tFric)/∂u_s with n FROZEN (flat master through p0)."""
    n = np.asarray(n, float)
    Pt = np.eye(3) - np.outer(n, n)

    def tfric(u):
        d = (X_s + u) - p0
        gT = Pt @ d
        gN = n @ d
        N = kn * max(0.0, -gN)
        return -tract_3d(gT, gpT, n, N, kt, mu, c, tmax)   # applied force on slave

    K = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        K[:, j] = -(tfric(e) - tfric(-e)) / (2 * h)        # K = −∂r/∂u (central)
    return K


def relerr(Ka, Kfd, kt):
    return np.max(np.abs(Ka - Kfd)) / max(np.max(np.abs(Ka)), kt)


# ============================================================ tests
print("ADR-41 A1 oracle — LadrunoFrictionKernel unified cone min(μN+c, τmax)")
n = np.array([0.0, 0.0, 1.0]); p0 = np.zeros(3)
kn, kt, mu = 1e6, 1e5, 0.5

# --- T1: BIT-FOR-BIT — unified (c=0, τmax≤0) reduces to pure Coulomb cap=μN ---
print("\nT1  bit-for-bit: unified (c=0, τmax=0) == shipped pure-Coulomb cap=μN")
for (gT, lbl) in [([0.03, 0.04, 0.0], "slip"), ([1e-3, 0.0, 0.0], "stick")]:
    gT = np.array(gT); N = kn * 1e-3
    t_uni = tract_3d(gT, np.zeros(3), n, N, kt, mu, 0.0, 0.0)
    t_cou = (lambda: (kt*gT) if kt*np.linalg.norm(gT) <= mu*N
             else mu*N*(kt*gT)/np.linalg.norm(kt*gT))()
    Kuni = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu, 0.0, 0.0)
    Kcou = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu)  # default == pure Coulomb
    check(f"{lbl}: traction identical", np.allclose(t_uni, t_cou, atol=1e-12))
    check(f"{lbl}: tangent identical", np.allclose(Kuni, Kcou, atol=1e-12))

# --- T2: Coulomb (slip) — Css/Csl vs FD to 1e-6, NON-symmetric (d_TN=μ) ---
print("\nT2  Coulomb slip: cap=μN; Css/Csl vs central-FD ≤1e-6; NON-symmetric")
Xs = np.array([0.03, 0.04, -1e-3])                    # penetrated + tangential ⇒ slip
d = Xs - p0; Pt = np.eye(3) - np.outer(n, n)
gT = Pt @ d; N = kn * max(0.0, -(n @ d))
Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu)
Kfd = fd_slave_block(Xs, np.zeros(3), n, p0, kn, kt, mu)
e = relerr(Ka, Kfd, kt); asym = np.max(np.abs(Ka - Ka.T)) / kt
check("Css/Csl vs FD ≤ 1e-6", e < 1e-6, f"relerr={e:.2e}")
check("non-symmetric (d_TN=μ active)", asym > 1e-3, f"asym/kt={asym:.2e}")

# --- T3: Coulomb + COHESION (slip) — cap=μN+c; still NON-sym (∂cap/∂N=μ) ---
print("\nT3  Coulomb+cohesion slip: cap=μN+c; FD ≤1e-6; still non-symmetric")
c = 0.3 * mu * N                                      # cohesion intercept
Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu, c=c)
Kfd = fd_slave_block(Xs, np.zeros(3), n, p0, kn, kt, mu, c=c)
e = relerr(Ka, Kfd, kt); asym = np.max(np.abs(Ka - Ka.T)) / kt
check("cap == μN + c", abs(friction_cap(N, mu, c, 0.0) - (mu * N + c)) < 1e-9)
check("Css/Csl vs FD ≤ 1e-6", e < 1e-6, f"relerr={e:.2e}")
check("non-symmetric (cohesion is pressure-independent ⇒ d_TN still μ)", asym > 1e-3)

# --- T4: τmax CAP binds (slip) — cap=τmax; Csl=0 ⇒ SYMMETRIC ---
print("\nT4  τmax cap binds: cap=τmax<μN; Csl=0 ⇒ SYMMETRIC; FD ≤1e-6")
tmax = 0.4 * mu * N                                   # below μN ⇒ cap clamps to τmax
Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu, tmax=tmax)
Kfd = fd_slave_block(Xs, np.zeros(3), n, p0, kn, kt, mu, tmax=tmax)
e = relerr(Ka, Kfd, kt); asym = np.max(np.abs(Ka - Ka.T)) / kt
check("cap == τmax (cap binds)", abs(friction_cap(N, mu, c, tmax) - tmax) < 1e-9)
check("Css vs FD ≤ 1e-6", e < 1e-6, f"relerr={e:.2e}")
check("SYMMETRIC (d_TN=0 on the τmax branch)", asym < 1e-9, f"asym/kt={asym:.2e}")

# --- T5: TRESCA (μ=0, c=τ) — constant cap; Csl=0 ⇒ SYMMETRIC ---
print("\nT5  Tresca (μ=0, cohesion=τ): constant cap; Csl=0 ⇒ SYMMETRIC; FD ≤1e-6")
tau = 0.5 * mu * N                                    # a constant shear strength
Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, 0.0, c=tau)
Kfd = fd_slave_block(Xs, np.zeros(3), n, p0, kn, kt, 0.0, c=tau)
e = relerr(Ka, Kfd, kt); asym = np.max(np.abs(Ka - Ka.T)) / kt
check("cap == τ (pressure-independent)", abs(friction_cap(N, 0.0, tau, 0.0) - tau) < 1e-9)
check("Css vs FD ≤ 1e-6", e < 1e-6, f"relerr={e:.2e}")
check("SYMMETRIC (μ=0 ⇒ d_TN=0)", asym < 1e-9, f"asym/kt={asym:.2e}")

# --- T6: STICK — elastic traction, tangent kt·P_t symmetric (cone-agnostic) ---
print("\nT6  stick (‖tT*‖≤cap): tFric=-tT*, tangent=kt·P_t symmetric")
Xst = np.array([1e-4, 0.0, -1e-3])                    # tiny slip ⇒ stick
d2 = Xst - p0; gTs = Pt @ d2; Ns = kn * max(0.0, -(n @ d2))
t = tract_3d(gTs, np.zeros(3), n, Ns, kt, mu)
Kst = slave_block_analytic(gTs, np.zeros(3), n, Ns, kn, kt, mu)
check("stick ⇒ tT == elastic kt·gT", np.allclose(t, kt * gTs, atol=1e-12))
check("stick tangent == kt·P_t", np.allclose(Kst, kt * Pt, atol=1e-9))
check("stick tangent symmetric", np.max(np.abs(Kst - Kst.T)) / kt < 1e-12)

# --- T7: return-map slip value — tT = cap·n̂, magnitude == cap ---
print("\nT7  slip return-map value: ‖tT‖ == cap, direction == slip dir")
tslip = tract_3d(gT, np.zeros(3), n, N, kt, mu)
check("‖tT‖ == μN (cap)", abs(np.linalg.norm(tslip) - mu * N) < 1e-6 * mu * N,
      f"‖tT‖={np.linalg.norm(tslip):.3f} cap={mu*N:.3f}")

print(f"\n{'='*60}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(A1 friction-cone oracle)\n{'='*60}")
sys.exit(1 if _fails else 0)
