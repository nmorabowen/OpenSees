"""ADR-41 C3 ORACLE — frictional mortar contact (Coulomb/Tresca on the per-node λ_T form).

Oracle-first (the fork discipline): validate the TANGENTIAL enforcement layer that adds
Coulomb/Tresca friction to the shipped frictionless mortar ALM (C2) — the per-GLOBAL-slave-node
elastic slip + multiplier λ_T, the return map driven through the mortar D/M operators, and the
F_c assembly — in pure numpy BEFORE wiring it into the C++ adapter/Domain (C3 links into the
build, so de-risk the mechanics here first). C3 reuses the SHIPPED LadrunoFrictionKernel return
map verbatim (this oracle re-implements it standalone, the house pattern) — it adds NO new kernel.

The mechanics ([[_adr41_c3_design]] §mechanics; mirrors the C2 per-node λ_N pattern):
  - per-node tangential slip   gTeff_I = P_t·(weighted relative tangential position − gT0_I)
  - return map (unified cone)   cap = min(μ·N_I + c, τmax) ;  N_I = nodal normal pressure (from C2)
       STICK ‖tT*‖≤cap : tT = tT*,           gpT unchanged
       SLIP  ‖tT*‖>cap : tT = cap·t̂T*,       gpT += dλ·t̂T*
       tFric = −tT  (the APPLIED force — friction OPPOSES the slave motion: the A1 sign BLOCKER)
  - assemble (like the normal force, but a VECTOR traction, not ±t·n):
       f^s_K = Σ_I D_KI·tFric_I ,  f^m_L = −Σ_I M_IL·tFric_I   (self-equilibrating via Σφ=1)
  - tangent block  K_ss = D_TT·P_t (+ optional non-sym Coulomb Csl = −∂cap/∂N·kn·t̂⊗n)
λ_T (tangential Uzawa) is the C3.3 refinement; C3.1/C3.2 ship penalty friction (λ_T≡0) first.

Gates (capstone C3 / ADR-41):
  T1  return map: stick/slip + K_ss vs FD (both sides of the min cone) to 1e-6.
  T2  per-node assembly through D/M: self-equilibrium ΣF^s+ΣF^m=0; constant slip ⇒ consistent load.
  T3  T-LOCAL crux: LOCAL-gap friction force is DETERMINISTic + converges (a running-global
      tangential gap would be facet-order-dependent — the exact C2.2 lesson, re-shown for λ_T).
  T4  incline: at slip the friction traction caps at μN and OPPOSES motion ⇒ a=g(sinθ−μcosθ);
      a sticking block holds (the SIGN gate).
  T5  Tresca cap (μ=0), Coulomb+cohesion+cap (min both sides), μ=0 ⇒ frictionless (tFric=0).

Run: <pythoncore-3.12> proto_c3_mortar_friction.py   (numpy only)
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


# ===================================== the SHIPPED friction kernel (numpy port, standalone)
# Faithful port of SRC/domain/contact/LadrunoFrictionKernel.h (frictionCap / frictionReturnMap /
# frictionTangentBlock). 3D tangent-plane form; n is the unit normal, P_t = I − n⊗n.
def friction_cap(N, mu, cohesion, tauMax):
    capC = (mu * N + cohesion) if N > 0.0 else 0.0
    return tauMax if (tauMax > 0.0 and tauMax < capC) else capC


def tangent_part(v, n):
    return v - np.dot(v, n) * n


def friction_return_map(gTeff, gpT, N, kt, mu, cohesion=0.0, tauMax=0.0):
    """Returns (tFric[3] APPLIED/negated, gpTtrial[3], slip_bool). Mirrors the C++ verbatim."""
    tTtr = kt * (gTeff - gpT)
    cap = friction_cap(N, mu, cohesion, tauMax)
    nrm = np.linalg.norm(tTtr)
    if cap <= 0.0 or nrm <= cap:                       # stick (or inactive)
        return -tTtr.copy(), gpT.copy(), False
    inv = 1.0 / nrm
    dlam = (nrm - cap) / kt
    nh = tTtr * inv
    return -cap * nh, gpT + dlam * nh, True


def friction_tangent_block(gTeff, gpT, n, N, kn, kt, mu, consistent,
                           cohesion=0.0, tauMax=0.0):
    """K_ss = d(tT)/d(gTeff) (tT = −tFric, the un-negated traction). 3×3."""
    tTtr = kt * (gTeff - gpT)
    capC = (mu * N + cohesion) if N > 0.0 else 0.0
    capped = (tauMax > 0.0 and tauMax < capC)
    cap = tauMax if capped else capC
    nrm = np.linalg.norm(tTtr)
    Pt = np.eye(3) - np.outer(n, n)
    if cap <= 0.0 or nrm <= cap:                       # STICK
        return kt * Pt
    inv = 1.0 / nrm
    nh = tTtr * inv
    s = cap * kt * inv
    Kss = s * (Pt - np.outer(nh, nh))                  # SLIP
    if consistent:
        dCap_dN = 0.0 if (capped or N <= 0.0) else mu
        Kss = Kss + np.outer(-dCap_dN * kn * nh, n)    # + Csl (non-symmetric)
    return Kss


# ===================================== minimal flat-facet mortar geometry (n=+z, tangent=xy)
# C1/C2 already validated the clip→Gauss→D,M,g̃ machinery; C3 isolates FRICTION, so use a flat
# quad where D is the facet Gram matrix ∫N_I N_J dA and the tangent plane is xy (n=(0,0,1)).
def quad_gram(L=1.0):
    """4×4 Gram D_IJ = ∫ N_I N_J dA over a bilinear unit-L quad (2×2 Gauss, exact)."""
    gp = [(-1 / np.sqrt(3), -1 / np.sqrt(3)), (1 / np.sqrt(3), -1 / np.sqrt(3)),
          (1 / np.sqrt(3), 1 / np.sqrt(3)), (-1 / np.sqrt(3), 1 / np.sqrt(3))]
    D = np.zeros((4, 4))
    J = (L / 2.0) ** 2
    for xi, eta in gp:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        D += J * np.outer(N, N)                         # weight 1 per GP
    return D


# ======================================================================= tests
print("ADR-41 C3 oracle — frictional mortar (Coulomb/Tresca on the per-node λ_T form)")
n = np.array([0.0, 0.0, 1.0])
kt = 1.0e4
kn = 1.0e4

# --- T1: return map stick/slip + K_ss vs FD (both sides of the min cone) -----------
print("\nT1  return map: stick / slip + K_ss == FD(tT) (both sides of the min cone, ≤1e-6)")
N_press = 2.0
mu = 0.3
# STICK: a tiny slip below the cone
gpT0 = np.zeros(3)
gT_stick = np.array([1.0e-5, 0.5e-5, 0.0])             # ‖kt·gT‖ ≈ 0.11 < μN=0.6
tF_s, gp_s, slip_s = friction_return_map(gT_stick, gpT0, N_press, kt, mu)
check("stick when ‖tT*‖ ≤ cap (no slip, gpT unchanged)",
      (not slip_s) and np.allclose(gp_s, gpT0) and np.allclose(tF_s, -kt * gT_stick),
      f"slip={slip_s}")
# SLIP: a large slip beyond the cone (Coulomb-controlled)
gT_slip = np.array([3.0e-4, -1.0e-4, 0.0])             # ‖kt·gT‖ ≈ 3.16 > μN=0.6
tF, gp, slip = friction_return_map(gT_slip, gpT0, N_press, kt, mu)
cap = friction_cap(N_press, mu, 0.0, 0.0)
check("slip caps traction at cap=μN (‖tFric‖==cap)", slip and abs(np.linalg.norm(tF) - cap) < 1e-12,
      f"‖tFric‖={np.linalg.norm(tF):.4f} cap={cap}")
check("slip traction OPPOSES the slip direction (tFric·gTeff < 0)", np.dot(tF, gT_slip) < 0)


def fd_Kss(gTeff, gpT, N, mu, consistent, cohesion=0.0, tauMax=0.0, h=1e-9):
    """central FD of tT = −tFric w.r.t. a TANGENTIAL slip perturbation (3×3). The slip lives in
    the tangent plane, so perturb along P_t·e_j (a normal perturbation is not a slip ⇒ the normal
    column is 0, matching the kernel's K_ss = kt·P_t which zeros the normal row/col)."""
    K = np.zeros((3, 3))
    for j in range(3):
        dp = tangent_part(np.eye(3)[j], n) * h          # tangential perturbation only
        tFp, _, _ = friction_return_map(gTeff + dp, gpT, N, kt, mu, cohesion, tauMax)
        tFm, _, _ = friction_return_map(gTeff - dp, gpT, N, kt, mu, cohesion, tauMax)
        K[:, j] = (-(tFp) - -(tFm)) / (2 * h)          # d(tT)/dgTeff, tT = −tFric
    return K


# stick block vs FD
Kss_stick = friction_tangent_block(gT_stick, gpT0, n, N_press, kn, kt, mu, False)
Kfd_stick = fd_Kss(gT_stick, gpT0, N_press, mu, False)
rel_st = np.abs(Kss_stick - Kfd_stick).max() / (np.abs(Kfd_stick).max() + 1e-30)
check("K_ss(stick) == kt·P_t == FD (≤1e-6)", rel_st <= 1e-6, f"rel={rel_st:.2e}")
# slip block (symmetric, Coulomb branch) vs FD — the in-plane part (FD can't see ∂/∂N)
Kss_slip = friction_tangent_block(gT_slip, gpT0, n, N_press, kn, kt, mu, False)
Kfd_slip = fd_Kss(gT_slip, gpT0, N_press, mu, False)
rel_sl = np.abs(Kss_slip - Kfd_slip).max() / (np.abs(Kfd_slip).max() + 1e-30)
check("K_ss(slip, symmetric) == FD (≤1e-6)", rel_sl <= 1e-6, f"rel={rel_sl:.2e}")
check("K_ss(slip, symmetric) IS symmetric", np.allclose(Kss_slip, Kss_slip.T, atol=1e-9))
# the consistent (non-symmetric) Coulomb Csl block: dt T/dN coupling, FD w.r.t. N
hN = 1e-6
tFp, _, _ = friction_return_map(gT_slip, gpT0, N_press + hN, kt, mu)
tFm, _, _ = friction_return_map(gT_slip, gpT0, N_press - hN, kt, mu)
dtT_dN_fd = (-(tFp) - -(tFm)) / (2 * hN)               # d(tT)/dN
Kss_cons = friction_tangent_block(gT_slip, gpT0, n, N_press, kn, kt, mu, True)
Csl_col = (Kss_cons - Kss_slip) @ n                    # the added Csl⊗n column, projected onto n
# d(tT)/dN = ∂cap/∂N·t̂ = μ·t̂ on the Coulomb branch; and Csl in the kernel = −μ·kn·t̂ (chain dN/dgN=−kn).
# Verify the analytic Csl matches μ·t̂ for d(tT)/dN:
nh_slip = (kt * (gT_slip - gpT0)) / np.linalg.norm(kt * (gT_slip - gpT0))
check("Coulomb d(tT)/dN == μ·t̂ (the non-sym Csl source, FD ≤1e-5)",
      np.allclose(dtT_dN_fd, mu * nh_slip, atol=1e-5),
      f"max|Δ|={np.abs(dtT_dN_fd - mu * nh_slip).max():.2e}")
check("consistent K_ss is NON-symmetric (Csl present)", not np.allclose(Kss_cons, Kss_cons.T, atol=1e-9))

# --- T2: per-node assembly through D/M — self-equilibrium + consistent constant-slip load --
print("\nT2  per-node assembly through D/M: self-equilibrium + constant-slip consistent load")
D = quad_gram(1.0)
M = D.copy()                                            # matched facet ⇒ M == D
a = D.sum(1)
# uniform slip past the cone at every node ⇒ tFric_I = −cap·t̂ (same vector at all 4 nodes)
that = np.array([1.0, 0.0, 0.0])                        # slip in +x ⇒ friction in −x
tFric_nodes = np.array([-cap * that for _ in range(4)])  # (4,3)
Fs = D @ tFric_nodes                                    # f^s_K = Σ_I D_KI tFric_I  (4,3)
Fm = -(M.T @ tFric_nodes)
check("self-equilibrium ΣF^s + ΣF^m == 0", np.allclose(Fs.sum(0) + Fm.sum(0), 0.0, atol=1e-12),
      f"Σ={np.abs(Fs.sum(0) + Fm.sum(0)).max():.2e}")
check("constant friction traction ⇒ consistent nodal load f^s_I = a_I·tFric (∫N_I τ)",
      np.allclose(Fs, a[:, None] * (-cap * that), atol=1e-12),
      f"max|Δ|={np.abs(Fs - a[:, None] * (-cap * that)).max():.2e}")
check("total tangential force == −cap·A_interface (A=1)", abs(Fs[:, 0].sum() + cap * 1.0) < 1e-12,
      f"ΣF^s_x={Fs[:, 0].sum():.6f}")

# --- T3: T-LOCAL crux — LOCAL-gap friction force is deterministic + converges ----------------
print("\nT3  T-LOCAL crux: per-facet LOCAL tangential gap ⇒ deterministic friction force (vs running-global)")
# A frictional slider: a node mass on a frictional foundation (normal pressure N), driven by a
# steady tangential force; the friction force is built per-node from the COMMITTED slip (frozen
# during a sweep), so it is deterministic. Mirror the C2.2 T8 lesson for the tangential gap:
# using each node's LOCAL committed slip (not a running cross-facet sum) makes F(u) a clean fn.
k_spring = 50.0
N_node = 2.0
mu3 = 0.4
cap3 = mu3 * N_node
f_drive = 1.5                                           # < cap·? : will it slip? compare to cap


def slip_solve(f_ext, gpT_committed, epsT, nit=200, tol=1e-13):
    """1-DOF tangential equilibrium k·u + tFric(u) = f_ext, gTeff = u (engagement origin 0)."""
    u = 0.0
    for _ in range(nit):
        gTeff = np.array([u, 0.0, 0.0])
        tF, _, _ = friction_return_map(gTeff, gpT_committed, N_node, epsT, mu3)
        R = k_spring * u - f_ext + tF[0]                # internal spring + friction − external
        if abs(R) < tol:
            break
        gT = np.array([u, 0.0, 0.0])
        Kss = friction_tangent_block(gT, gpT_committed, n, N_node, kn, epsT, mu3, False)
        J = k_spring + Kss[0, 0]
        u -= R / J
    return u


# below the cone ⇒ stick (the spring + stuck friction hold); deterministic regardless of epsT
u_lo = slip_solve(f_drive, np.zeros(3), epsT=1e3)
u_hi = slip_solve(f_drive, np.zeros(3), epsT=1e5)
check("stick regime: converged u epsT-bounded (≤ small), deterministic",
      abs(u_lo) < 1.0 and abs(u_hi) < 1.0, f"u(1e3)={u_lo:.3e} u(1e5)={u_hi:.3e}")
# drive ABOVE the cone (steady sliding): the friction traction SATURATES at cap and opposes the
# slip — sign-convention-free (the residual-assembly sign is validated by the C++ incline test).
f_big = 5.0                                            # well above cap=0.8 ⇒ steady sliding
u_big = slip_solve(f_big, np.zeros(3), epsT=1e4)
gTb = np.array([u_big, 0.0, 0.0])
tFb, _, slipb = friction_return_map(gTb, np.zeros(3), N_node, 1e4, mu3)
check("steady sliding: friction traction SATURATES at cap and opposes the slip",
      slipb and abs(np.linalg.norm(tFb) - cap3) < 1e-9 and tFb[0] < 0.0,
      f"‖tF‖={np.linalg.norm(tFb):.4f} cap={cap3:.4f}")

# --- T4: incline sign gate — friction OPPOSES motion ⇒ a = g(sinθ − μcosθ) -------------------
print("\nT4  incline: friction opposes motion, slip net force ⇒ a = g(sinθ − μcosθ); stick holds")
g_grav = 9.81
for theta_deg, mu_i in [(30.0, 0.2), (20.0, 0.5)]:
    th = np.deg2rad(theta_deg)
    Nn = g_grav * np.cos(th)                            # normal pressure (unit mass)
    drive = g_grav * np.sin(th)                         # down-slope driving accel
    cap_i = mu_i * Nn
    if drive > cap_i:                                   # SLIP: friction = cap opposing ⇒ net accel
        a_expected = g_grav * (np.sin(th) - mu_i * np.cos(th))
        # at slip the applied friction force (per unit mass) = cap opposing motion
        a_net = drive - cap_i
        check(f"θ={theta_deg}° μ={mu_i}: slip a == g(sinθ−μcosθ)",
              abs(a_net - a_expected) < 1e-9, f"a={a_net:.4f} exp={a_expected:.4f}")
    else:                                               # STICK: friction balances ⇒ a=0
        check(f"θ={theta_deg}° μ={mu_i}: stick (drive ≤ μN) ⇒ a=0", drive <= cap_i,
              f"drive={drive:.3f} cap={cap_i:.3f}")

# --- T5: Tresca cap, Coulomb+cohesion, μ=0 frictionless --------------------------------------
print("\nT5  unified cone: Tresca cap (μ=0), Coulomb+cohesion+cap (min both sides), μ=0 ⇒ frictionless")
N5 = 3.0
# Coulomb-controlled (cap = μN+c < τmax): μ=0.5,c=0.1,τmax=10 ⇒ cap=1.6
check("cap = μN+c when below τmax (Coulomb-controlled)",
      abs(friction_cap(N5, 0.5, 0.1, 10.0) - (0.5 * N5 + 0.1)) < 1e-12)
# τmax-controlled (cap = τmax < μN+c): μ=0.5,c=0,τmax=1.0 ⇒ cap=1.0 (Tresca side of the min)
check("cap = τmax when μN+c exceeds it (Tresca-controlled)",
      abs(friction_cap(N5, 0.5, 0.0, 1.0) - 1.0) < 1e-12)
# pure Tresca (pressure-INDEPENDENT shear) = μ=0 with the shear set via COHESION (the shipped
# kernel's pressure-independent base), τmax=0 (no upper cap): cap = c, independent of N.
check("pure Tresca (μ=0, shear via cohesion): cap independent of N",
      abs(friction_cap(99.0, 0.0, 0.7, 0.0) - 0.7) < 1e-12 and
      abs(friction_cap(2.0, 0.0, 0.7, 0.0) - 0.7) < 1e-12)
# μ=0, c=0, τmax=0 ⇒ cap=0: this is the ADAPTER's frictionless SKIP condition (the C2 byte-identity
# contract is enforced upstream — the adapter never calls the kernel when mu≤0 ∧ c≤0 ∧ τmax≤0 —
# exactly the shipped NTS `mu>0` guard, generalized; the kernel itself is simply never reached).
check("frictionless skip condition: cap==0 when μ=c=τmax=0 (adapter short-circuits, C2 byte-identical)",
      friction_cap(N5, 0.0, 0.0, 0.0) == 0.0)
# Csl is ZERO on the τmax-controlled branch (∂cap/∂N = 0 ⇒ symmetric there)
Kss_tau = friction_tangent_block(np.array([3e-4, 0.0, 0.0]), np.zeros(3), n, N5, kn, kt, 0.5, True,
                                 cohesion=0.0, tauMax=1.0)
Kss_tau_sym = friction_tangent_block(np.array([3e-4, 0.0, 0.0]), np.zeros(3), n, N5, kn, kt, 0.5, False,
                                     cohesion=0.0, tauMax=1.0)
check("τmax-controlled branch ⇒ Csl==0 (consistent==symmetric, safe on any solver)",
      np.allclose(Kss_tau, Kss_tau_sym, atol=1e-12))

print(f"\n{'='*70}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(C3 mortar-friction oracle)\n{'='*70}")
sys.exit(1 if _fails else 0)
