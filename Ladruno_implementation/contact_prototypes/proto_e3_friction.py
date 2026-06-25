"""ADR-57 E3 ORACLE — edge-edge FRICTION (Coulomb/Tresca; reuse LadrunoFrictionKernel).

Oracle-first (the fork discipline): validate the tangential friction layer in pure numpy BEFORE
transcribing it into the LadrunoContactFE EDGE_EDGE mode. No OpenSees, no build. E3 adds the
tangential traction on top of the E2 normal penalty — the edge-edge contact is the FOURTH consumer
of the ONE shipped friction return map (NTS, mortar, SOFT=2, now edge-edge), with NO new friction
code: only the slip MEASURE and the SCATTER are edge-specific.

What E3 owns (the ADR-57 E3 gate row + §4):
  - the edge tangent frame is t₁=ê_s, t₂=n×ê_s (n=ê_s×ê_m ⊥ ê_s by construction, so the draft's
    ê_s⊥ projection was a no-op); but the return map operates on the full 3-vec tangential slip, so
    the kernel call is identical to NTS/mortar — only n is needed (tangentPart projects out n).
  - SLIP FROM DISPLACEMENT, NOT POSITION (the C3.1 / NTS-SEGMENT trap): the closest-point
    construction makes the relative POSITION at (s*,t*) purely NORMAL (w ∥ n), so positions carry
    NO tangential information; the slip is the relative DISPLACEMENT of the two closest points,
    tangential part — gT = tangentPart((1−s)u_a0+s u_a1 − (1−t)u_b0 − t u_b1, n).
  - the return map tT* = kt(gTeff−gpT), STICK ‖tT*‖≤cap / SLIP ‖tT*‖=cap; cap = min(μN+c, τmax),
    N = εN⟨−gN⟩ the normal pressure; tFric = −tT (APPLIED, motion-opposing — the kernel BLOCKER-1).
  - the incline criterion a = g(sinθ−μcosθ): the cone threshold ‖drive‖>μN IS the slide criterion
    tanθ>μ; the net tangential accel is g(sinθ−μcosθ) when sliding, 0 when sticking.
  - self-equilibrating SCATTER via the SAME edge weights as the normal force (f_a0=(1−s)tFric,
    f_a1=s tFric, f_b0=−(1−t)tFric, f_b1=−t tFric ⇒ Σf=0, since Σw = (1−s)+s−(1−t)−t = 0).
  - the friction TANGENT K[A][B]=w_A w_B·K_ss (K_ss = frictionTangentBlock), SYMMETRIC by default
    (the Coulomb Csl pressure-coupling behind -edgeConsistentTan); FD-checked.
  - μ=0 ∧ c=0 ∧ τmax=0 ⇒ cap=0 ⇒ tFric=0 ⇒ BYTE-IDENTICAL to the E2 frictionless path.

Run: <pythoncore-3.12> proto_e3_friction.py   (numpy only)
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


# ============================================ LadrunoFrictionKernel (mirror the C++ verbatim)
def tangent_part(v, n):
    return v - (v @ n) * n


def friction_cap(N, mu, cohesion, tauMax):
    capC = (mu * N + cohesion) if N > 0.0 else 0.0
    return tauMax if (tauMax > 0.0 and tauMax < capC) else capC


def friction_return_map(gTeff, gpT, N, kt, mu, cohesion=0.0, tauMax=0.0):
    """returns (tFric, gpTtrial, slipping). tFric = APPLIED (negated) traction opposing motion."""
    tTtr = kt * (gTeff - gpT)
    cap = friction_cap(N, mu, cohesion, tauMax)
    norm = np.linalg.norm(tTtr)
    if cap <= 0.0 or norm <= cap:
        return -tTtr, gpT.copy(), False
    nh = tTtr / norm
    dlam = (norm - cap) / kt
    return -cap * nh, gpT + dlam * nh, True


def friction_tangent_block(gTeff, gpT, n, N, kn, kt, mu, consistent, cohesion=0.0, tauMax=0.0):
    tTtr = kt * (gTeff - gpT)
    capC = (mu * N + cohesion) if N > 0.0 else 0.0
    capped = (tauMax > 0.0 and tauMax < capC)
    cap = tauMax if capped else capC
    nrm = np.linalg.norm(tTtr)
    Pt = np.eye(3) - np.outer(n, n)
    if cap <= 0.0 or nrm <= cap:
        return kt * Pt
    nh = tTtr / nrm
    s = cap * kt / nrm
    Kss = s * (Pt - np.outer(nh, nh))
    if consistent:
        dCap_dN = 0.0 if (capped or N <= 0.0) else mu
        Kss = Kss + np.outer(-dCap_dN * kn * nh, n)
    return Kss


# ============================================ the E3 edge-edge friction layer
def friction_active(mu, cohesion, tauMax):
    """the CALLER's guard: the friction block is assembled ONLY when the cone is real. Byte-identity
    to the frictionless E2 path comes from SKIPPING the block here (the kernel's cap<=0 branch
    returns the RAW elastic traction, NOT zero — so the guard, not the kernel, gives byte-identity;
    exactly the NTS `mu>0` / mortar `mu>0 || c>0 || τmax>0` short-circuit)."""
    return mu > 0.0 or cohesion > 0.0 or tauMax > 0.0


def edge_weights(s, t):
    """the 4 nodal scatter weights for [a0,a1,b0,b1] (the E2 B-operator shape)."""
    return np.array([1 - s, s, -(1 - t), -t])


def edge_slip(p1, q1, p2, q2, ua0, ua1, ub0, ub1, n):
    """tangential slip from DISPLACEMENT at the closest point (the C3.1 trap)."""
    k = closest_pt_segseg(p1, q1, p2, q2)
    s, t = k["s"], k["t"]
    d_disp = ((1 - s) * ua0 + s * ua1) - ((1 - t) * ub0 + t * ub1)   # rel. DISPLACEMENT at the cp
    return tangent_part(d_disp, n), s, t


def edge_friction_force(s, t, tFric):
    """scatter tFric to the 4 nodes with the edge weights; returns 4x3."""
    w = edge_weights(s, t)
    return np.array([wi * tFric for wi in w])


def edge_friction_tangent(s, t, Kss):
    """K[A][B] = w_A w_B K_ss (12x12)."""
    w = edge_weights(s, t)
    K = np.zeros((12, 12))
    for A in range(4):
        for B in range(4):
            K[3*A:3*A+3, 3*B:3*B+3] = w[A] * w[B] * Kss
    return K


# ============================================================================ tests
print("ADR-57 E3 oracle — edge-edge Coulomb/Tresca friction (reuse LadrunoFrictionKernel)")

# canonical crossed bars: slave edge ∥ x at z=zc, master edge ∥ y at z=0, n≈+z (toward slave).
def bars(zc=-0.02):
    return (np.array([-1, 0, zc]), np.array([1, 0, zc]),
            np.array([0, -1, 0.]), np.array([0, 1, 0.]))
ORIENT = np.array([0, 0, 1.])
EPS, KT, MU = 1.0e3, 1.0e3, 0.3

# --- T1: slip is from DISPLACEMENT, not POSITION (the C3.1 trap) ---
print("\nT1  slip from DISPLACEMENT: the closest-point relative POSITION is purely normal")
p1, q1, p2, q2 = bars()
k = closest_pt_segseg(p1, q1, p2, q2)
n = edge_normal(q1 - p1, q2 - p2, ORIENT)
# relative POSITION at the closest point ∥ n ⇒ tangential part ≈ 0 (carries NO slip)
s, t = k["s"], k["t"]
relpos = ((1 - s) * p1 + s * q1) - ((1 - t) * p2 + t * q2)
check("relative POSITION tangential part ≈ 0 (positions purely normal)",
      np.linalg.norm(tangent_part(relpos, n)) < TOL, f"|posT|={np.linalg.norm(tangent_part(relpos,n)):.2e}")
# a tangential DISPLACEMENT of the slave edge (slide along x) DOES register as slip
ua = np.array([0.05, 0, 0.])                      # slave edge slides +x
gT, s, t = edge_slip(p1, q1, p2, q2, ua, ua, np.zeros(3), np.zeros(3), n)
check("displacement slide +x ⇒ slip ≈ (0.05,0,0) (tangential, from DISPLACEMENT)",
      np.allclose(gT, [0.05, 0, 0], atol=1e-6), f"gT={gT}")

# --- T2: return map stick vs slip on the edge tangent plane; friction OPPOSES motion ---
print("\nT2  return map: STICK (small slip) vs SLIP (cap), friction opposes the motion")
N = EPS * 0.02                                    # normal pressure tn = εN|gN|, gN=-0.02
small = np.array([1e-4, 0, 0.])                   # tiny slip ⇒ stick
tF_s, gp_s, slip_s = friction_return_map(small, np.zeros(3), N, KT, MU)
check("small slip ⇒ STICK (‖tF‖<cap, elastic)", (not slip_s) and np.linalg.norm(tF_s) < MU * N,
      f"|tF|={np.linalg.norm(tF_s):.3f} cap={MU*N:.3f}")
check("stick traction opposes motion (tF·slip < 0)", tF_s @ small < 0)
big = np.array([1.0, 0, 0.])                      # large slip ⇒ slip, capped at μN
tF_b, gp_b, slip_b = friction_return_map(big, np.zeros(3), N, KT, MU)
check("large slip ⇒ SLIP at the cap ‖tF‖=μN", slip_b and abs(np.linalg.norm(tF_b) - MU * N) < TOL,
      f"|tF|={np.linalg.norm(tF_b):.4f} μN={MU*N:.4f}")
check("slip traction opposes motion (tF·slip < 0)", tF_b @ big < 0)

# --- T3: the incline criterion a = g(sinθ−μcosθ) (the cone threshold IS the slide criterion) ---
print("\nT3  incline a = g(sinθ−μcosθ): the friction cone threshold == the slide criterion tanθ>μ")
g = 1.0
for thd in (10.0, 16.7, 30.0, 45.0):
    th = np.deg2rad(thd)
    N_inc, drive = g * np.cos(th), g * np.sin(th)     # block on a θ-incline (unit mass)
    cap = MU * N_inc
    a = max(0.0, drive - cap)                          # net tangential accel (0 if stuck)
    sliding = np.tan(th) > MU
    a_pred = g * (np.sin(th) - MU * np.cos(th)) if sliding else 0.0
    check(f"θ={thd:>4}°: a={'slide' if sliding else 'stick'} == g(sinθ−μcosθ)",
          abs(a - max(0.0, a_pred)) < TOL and (a > 0) == sliding, f"a={a:.4f} (tanθ={np.tan(th):.3f} vs μ={MU})")

# --- T4: Tresca cap (μ=0, c=τmax) — pressure-independent constant ---
print("\nT4  Tresca cap (μ=0, cohesion=c) — traction capped at c, independent of N")
C = 5.0
tF_tr, _, slip_tr = friction_return_map(np.array([1.0, 0, 0.]), np.zeros(3), 999.0, KT, 0.0, cohesion=C)
check("Tresca: ‖tF‖ == c (pressure-independent)", abs(np.linalg.norm(tF_tr) - C) < TOL,
      f"|tF|={np.linalg.norm(tF_tr):.4f} c={C}")

# --- T5: μ=0 ∧ c=0 ∧ τmax=0 ⇒ the CALLER short-circuits ⇒ byte-identical to frictionless E2 ---
print("\nT5  μ=0 ∧ c=0 ∧ τmax=0 ⇒ friction lane INACTIVE (caller short-circuit ⇒ byte-identical to E2)")
check("frictionless ⇒ friction lane INACTIVE (C++ skips the block — no kernel call, no slot touch)",
      not friction_active(0.0, 0.0, 0.0))
check("any friction param > 0 ⇒ ACTIVE (Coulomb / cohesion / Tresca cap)",
      friction_active(0.3, 0, 0) and friction_active(0, 5.0, 0) and friction_active(0, 0, 2.0))
# NB: the kernel's cap<=0 branch returns the RAW elastic traction (NOT zero) — byte-identity is the
# guard's job, not the kernel's. Demonstrate the trap (why the C++ MUST short-circuit, like NTS/mortar):
tF0, _, _ = friction_return_map(np.array([0.3, 0.2, 0.]), np.zeros(3), N, KT, 0.0)
check("trap: kernel with cap=0 returns NONZERO elastic (⇒ the block MUST be guarded out)",
      np.linalg.norm(tF0) > TOL, f"|tF(cap=0)|={np.linalg.norm(tF0):.3f}")

# --- T6: self-equilibrium + symmetric tangent + FD-check the assembled friction tangent ---
print("\nT6  scatter self-equilibrium Σf=0; tangent K=w⊗w·K_ss symmetric + FD-checked")
# a skew penetrating pair with a tangential slip (stick regime ⇒ a smooth tangent to FD)
p1, q1 = np.array([-1, 0.1, -0.02]), np.array([1, -0.05, -0.02])
p2, q2 = np.array([0.1, -1, 0.0]), np.array([-0.1, 1, 0.0])
n = edge_normal(q1 - p1, q2 - p2, ORIENT)
k = closest_pt_segseg(p1, q1, p2, q2)
s, t = k["s"], k["t"]
ua0, ua1 = np.array([2e-4, 1e-4, 0.]), np.array([1e-4, -1e-4, 0.])    # small ⇒ stick
ub0, ub1 = np.zeros(3), np.zeros(3)
N6 = EPS * 0.02
gT, _, _ = edge_slip(p1, q1, p2, q2, ua0, ua1, ub0, ub1, n)
tF, gpt, slipping = friction_return_map(gT, np.zeros(3), N6, KT, MU)
f = edge_friction_force(s, t, tF)
check("Σf == 0 (self-equilibrating scatter, Σw=0)", np.allclose(f.sum(axis=0), 0.0, atol=TOL))
Kss = friction_tangent_block(gT, np.zeros(3), n, N6, EPS, KT, MU, consistent=False)
K = edge_friction_tangent(s, t, Kss)
check("symmetric tangent (default, solver-safe)", np.allclose(K, K.T, atol=TOL),
      f"max asym={np.abs(K-K.T).max():.2e}")
# FD: freeze n,s,t; r(u)=scatter(tFric(gTeff(u))); K_solver = −∂r/∂u == w⊗w·K_ss
U0 = np.concatenate([ua0, ua1, ub0, ub1]).astype(float)
def fric_resid(U):
    u = U.reshape(4, 3)
    d = ((1 - s) * u[0] + s * u[1]) - ((1 - t) * u[2] + t * u[3])
    gTe = tangent_part(d, n)                       # frozen n
    tf, _, _ = friction_return_map(gTe, np.zeros(3), N6, KT, MU)
    return edge_friction_force(s, t, tf).reshape(-1)
h = 1e-8
Kfd = np.zeros((12, 12))
for j in range(12):
    Up, Um = U0.copy(), U0.copy(); Up[j] += h; Um[j] -= h
    Kfd[:, j] = (fric_resid(Up) - fric_resid(Um)) / (2 * h)
check("K == −∂r/∂u FD (the consistent stick tangent)", np.abs(-Kfd - K).max() < 1e-4,
      f"max err={np.abs(-Kfd - K).max():.2e}")
# the consistent (non-symmetric Coulomb Csl) branch exists and differs from symmetric on SLIP
gT_slip, _, _ = edge_slip(p1, q1, p2, q2, np.array([1.0, 0, 0.]), np.array([1.0, 0, 0.]), ub0, ub1, n)
Kss_sym = friction_tangent_block(gT_slip, np.zeros(3), n, N6, EPS, KT, MU, consistent=False)
Kss_con = friction_tangent_block(gT_slip, np.zeros(3), n, N6, EPS, KT, MU, consistent=True)
check("consistent Csl branch is NON-symmetric on slip (≠ the symmetric default)",
      not np.allclose(Kss_con, Kss_con.T) and not np.allclose(Kss_sym, Kss_con))

print(f"\n{'='*66}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E3 edge-edge friction oracle)\n{'='*66}")
sys.exit(1 if _fails else 0)
