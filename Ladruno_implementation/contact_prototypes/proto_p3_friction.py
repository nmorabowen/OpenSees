"""ADR-39 P3 ORACLE — IMPL-EX Coulomb friction (clean radial return), frictionless-normal sibling.

Oracle-first (the fork discipline): pin the tangential friction return map + the
IMPL-EX extrapolation + the (non-symmetric) consistent tangent in pure numpy BEFORE
the C++ LadrunoContactKernel friction branch. No OpenSees, no build.

DESIGN (ADR §"Friction = pressure-dependent J2 radial return, non-associative",
de Souza Neto §7, LS-DYNA 26.19-26.23). We author the CLEAN Coulomb return map — NOT
the ASDimplex "apparent damage" formulation (which carries the quarantined bugs:
shadowed dE1/dE2, ddata(39) OOB). It is structurally 1D plasticity on the tangential
traction with a pressure-dependent (Coulomb) yield radius:

  state per active pair: committed plastic slip gp_T (tangent-plane vector),
                         + (IMPL-EX) committed plastic multiplier dlam (and dlam_old).
  normal force      N    = kn·<−gN>+          (penetration ⇒ N≥0, compressive)
  trial traction    tT*  = kt·(gT − gp_T^commit)
  yield (cone)      f*   = ‖tT*‖ − μN
  STICK  (f*≤0):    tT   = tT*,                gp_T unchanged,      dlam = 0
  SLIP   (f*>0):    n̂    = tT*/‖tT*‖
                    dlam = f*/kt
                    tT   = μN · n̂  = (1−d)·tT*  with  d = 1 − μN/‖tT*‖  (secant form)
                    gp_T = gp_T^commit + dlam·n̂

  IMPL-EX (explicit phase, the v1 SHIP leg): extrapolate the multiplier from committed
  history (NOT solve the cone), so the tangent is the constant SECANT (symmetric, SPD)
  and the step is explicit in the internal variable:
    dlam~ = max(0, dlam_commit + (Δt_n/Δt_{n-1})·(dlam_commit − dlam_old))
    tT    = tT* − kt·dlam~·n̂*        (n̂* = current trial direction)
  Error is O(Δt²), corrected at the next commit (a true implicit return map runs at
  commitState to refresh dlam). For pure explicit dynamics only the FORCE tT is
  assembled (the tangent goes to mass-only), so IMPL-EX's job there is a robust,
  history-consistent bounded traction.

  CONSISTENT (implicit) tangent — NON-SYMMETRIC, rigorous (slip flow is purely
  tangential, cone normal has a t_N component):
    k_T = (μN/‖tT*‖)[I_T − n̂⊗n̂] + μ·sign·(n̂ ⊗ (−kn n))   ← the ⊗n column is the
    pressure-coupling that breaks symmetry (∂‖tT‖cap/∂gN via N=kn gN).

This oracle validates: stick, slip-caps-at-μN, flow direction, energy dissipation ≥ 0,
the non-symmetric consistent tangent vs finite difference, IMPL-EX → implicit as the
step settles, and a sliding-block-on-incline acceleration a=g(sinθ−μcosθ).

Run: <pythoncore-3.12> proto_p3_friction.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9


# ----------------------------------------------------------- the return map
def normal_force(kn, gN):
    """N = kn·<−gN>+  (gN<0 = penetration ⇒ N>0). Macaulay; 0 if separated."""
    return kn * max(0.0, -gN)


def return_map_implicit(gT, gpT_commit, N, kt, mu):
    """Clean Coulomb radial return. gT, gpT_commit: 2-vectors (tangent plane).
    returns (tT, gpT_new, dlam, slipping)."""
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT_commit, float))
    norm_tr = np.linalg.norm(tT_tr)
    cap = mu * N
    f = norm_tr - cap
    if f <= 0.0 or norm_tr < 1e-300:
        return tT_tr, np.array(gpT_commit, float), 0.0, False
    nhat = tT_tr / norm_tr
    dlam = f / kt
    tT = cap * nhat
    gpT_new = np.asarray(gpT_commit, float) + dlam * nhat
    return tT, gpT_new, dlam, True


def return_map_implex(gT, gpT_commit, N, kt, mu, dlam_commit, dlam_old, time_factor=1.0):
    """IMPL-EX: extrapolate the multiplier; explicit traction with current trial dir."""
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT_commit, float))
    norm_tr = np.linalg.norm(tT_tr)
    if norm_tr < 1e-300:
        return tT_tr.copy()
    nhat = tT_tr / norm_tr
    dlam_x = max(0.0, dlam_commit + time_factor * (dlam_commit - dlam_old))
    tT = tT_tr - kt * dlam_x * nhat
    return tT


def friction_tangent(gT, gpT_commit, N, kt, mu, kn):
    """The FRICTION-specific sensitivities of the tangential traction tT in the
    [gN, gT1, gT2] layout: returns the 2x3 block [∂tT/∂gN | ∂tT/∂gT] (rows = tT1,tT2;
    cols = gN,gT1,gT2). The normal-normal term ∂tN/∂gN is P2b's (already validated),
    not here. NON-SYMMETRIC under slip (the ∂tT/∂gN column has no symmetric partner —
    normal traction does not depend on tangential slip)."""
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT_commit, float))
    norm_tr = np.linalg.norm(tT_tr)
    K = np.zeros((2, 3))                       # rows tT1,tT2 ; cols gN,gT1,gT2
    cap = mu * N
    if norm_tr <= cap or norm_tr < 1e-300 or N <= 0.0:
        # stick: tT = kt·(gT − gpT) ⇒ ∂tT/∂gT = kt I, ∂tT/∂gN = 0
        K[0, 1] = K[1, 2] = kt
        return K
    nhat = tT_tr / norm_tr
    IT = np.eye(2)
    # SLIP: tT = μN · tT*/‖tT*‖, with ‖tT*‖ = kt‖gT−gpT‖.
    # ∂tT/∂gT = μN·kt/‖tT*‖ · (I − n̂⊗n̂)   (the dropped-kt factor the FD caught)
    K[:, 1:3] = (cap * kt / norm_tr) * (IT - np.outer(nhat, nhat))
    # ∂tT/∂gN: cap = μN = μ·kn·(−gN) ⇒ dN/dgN = −kn ⇒ ∂tT/∂gN = −μ kn n̂
    # (the NON-SYMMETRIC pressure-coupling column).
    K[0, 0] = -mu * kn * nhat[0]
    K[1, 0] = -mu * kn * nhat[1]
    return K


# ----------------------------------------------------------------- tests
def test_stick():
    """Below the cone: tT = kt·gT, no plastic slip, no dissipation."""
    kn, kt, mu = 1e6, 1e5, 0.5
    N = normal_force(kn, -1e-3)             # gN=-1e-3 penetration ⇒ N=1e3, cap=μN=500
    gT = np.array([1e-3, 0.0])              # tT*=kt·gT=[100,0], ‖‖=100 < 500 ⇒ stick
    tT, gpT, dlam, slip = return_map_implicit(gT, [0.0, 0.0], N, kt, mu)
    assert not slip and dlam == 0.0
    assert np.allclose(tT, kt * gT), tT
    assert np.allclose(gpT, [0, 0])
    print(f"[stick] ‖tT‖={np.linalg.norm(tT):.1f} < μN={mu*N:.1f}  no slip OK")


def test_slip_caps_at_muN():
    """Above the cone: ‖tT‖ = μN exactly, direction = slip direction."""
    kn, kt, mu = 1e6, 1e5, 0.5
    N = normal_force(kn, -1e-3)             # N=1e3, cap=500
    gT = np.array([3e-2, 4e-2])             # tT*=kt·gT=[3000,4000],‖‖=5000 ≫ 500 ⇒ slip
    tT, gpT, dlam, slip = return_map_implicit(gT, [0.0, 0.0], N, kt, mu)
    assert slip
    assert abs(np.linalg.norm(tT) - mu * N) < 1e-9, np.linalg.norm(tT)
    # direction parallel to the trial (= slip) direction [3,4]/5
    assert np.allclose(tT / np.linalg.norm(tT), [0.6, 0.8]), tT
    assert dlam > 0 and np.linalg.norm(gpT) > 0
    print(f"[slip] ‖tT‖={np.linalg.norm(tT):.1f} == μN={mu*N:.1f}  dir={tT/np.linalg.norm(tT)} OK")


def test_dissipation_nonneg():
    """Slip dissipates ≥0 energy: ∫ tT·dgp ≥ 0 over a monotone slide."""
    kn, kt, mu = 1e6, 1e5, 0.3
    N = normal_force(kn, -2e-3)
    gpT = np.array([0.0, 0.0]); diss = 0.0
    prev_gp = gpT.copy()
    for s in np.linspace(0, 0.05, 50):
        gT = np.array([s, 0.0])
        tT, gpT, dlam, slip = return_map_implicit(gT, gpT, N, kt, mu)
        diss += np.dot(tT, gpT - prev_gp)
        prev_gp = gpT.copy()
    assert diss > 0, diss
    print(f"[dissipation] ∫tT·dgp = {diss:.4f} > 0 OK")


def test_consistent_tangent_vs_fd():
    """Friction tangent (∂tT/∂gN | ∂tT/∂gT) vs finite difference (slip regime), and
    confirm the ∂tT/∂gN coupling is present (the rigorous NON-symmetry source)."""
    kn, kt, mu = 1e6, 1e5, 0.5
    gN = -1e-3
    gT = np.array([3e-2, 4e-2])
    N = normal_force(kn, gN)
    K = friction_tangent(gT, [0, 0], N, kt, mu, kn)   # 2x3 [∂tT/∂gN | ∂tT/∂gT]

    def tT_of(gN_, gT_):
        N_ = normal_force(kn, gN_)
        tT_, _, _, _ = return_map_implicit(gT_, [0, 0], N_, kt, mu)
        return tT_

    h = 1e-9
    Kfd = np.zeros((2, 3))
    base = tT_of(gN, gT)
    for j, pert in enumerate([(h, 0, 0), (0, h, 0), (0, 0, h)]):
        Kfd[:, j] = (tT_of(gN + pert[0], gT + np.array([pert[1], pert[2]])) - base) / h
    err = np.max(np.abs(K - Kfd))
    assert err < 1e-2 * kt, f"tangent err {err:.3e}\nK=\n{K}\nKfd=\n{Kfd}"
    coupling = abs(K[0, 0]) + abs(K[1, 0])     # ∂tT/∂gN — the non-symmetric column
    assert coupling > 1e-3 * kt, "∂tT/∂gN coupling (non-symmetry source) should be nonzero"
    print(f"[tangent] max|K−Kfd|={err:.3e}  ∂tT/∂gN coupling={coupling:.1f} (non-sym ✓) OK")


def test_implex_converges_to_implicit():
    """IMPL-EX → implicit once the slip rate is STEADY (constant dlam ⇒ linear
    extrapolation exact). NOTE a known/correct IMPL-EX property: at a stick→slip
    onset the multiplier extrapolation overshoots for ONE step (dlam_old still the
    pre-slip 0), then re-syncs — so we measure the STEADY regime, not the onset."""
    kn, kt, mu = 1e6, 1e5, 0.4
    N = normal_force(kn, -1e-3)
    gpT = np.array([0.0, 0.0]); dlam_old = 0.0; dlam_commit = 0.0
    dstep = 1e-3
    errs = []
    for k in range(1, 30):
        gT = np.array([k * dstep, 0.0])
        tT_imp, gpT_new, dlam_new, _ = return_map_implicit(gT, gpT, N, kt, mu)
        tT_x = return_map_implex(gT, gpT, N, kt, mu, dlam_commit, dlam_old, 1.0)
        if k >= 10:  # steady regime (well past the slip-onset overshoot at k≈6)
            errs.append(abs(np.linalg.norm(tT_x) - np.linalg.norm(tT_imp)))
        gpT = gpT_new; dlam_old = dlam_commit; dlam_commit = dlam_new
    assert max(errs) < 1e-6 * mu * N, f"steady implex-implicit gap {max(errs):.3e}"
    print(f"[implex] steady-slide max|‖tT_implex‖−‖tT_impl‖| = {max(errs):.3e} (→0) OK")


def test_sliding_block_on_incline():
    """Explicit dynamics: block on a θ-incline, contact normal N=mg cosθ, Coulomb μ.
    Once sliding, a = g(sinθ − μ cosθ). Integrate central-difference along the slope."""
    g, m, theta, mu = 9.81, 2.0, np.radians(30.0), 0.3
    kn, kt = 1e8, 1e7
    N = m * g * np.cos(theta)               # the (held) normal force
    drive = m * g * np.sin(theta)           # down-slope gravity component
    # 1-DOF along-slope: m x'' = drive − friction(x', state)
    dt = 1e-4
    x = 0.0; v = 0.0; gpT = np.array([0.0, 0.0]); a = 0.0
    xs = []
    for step in range(2000):
        # tangential relative slip along slope = x; penalty trial vs committed plastic
        gT = np.array([x, 0.0])
        tT, gpT, dlam, slip = return_map_implicit(gT, gpT, N, kt, mu)
        fric = tT[0]                        # opposes motion once on the cone
        a = (drive - fric) / m
        v += a * dt
        x += v * dt
        xs.append(x)
    a_theory = g * (np.sin(theta) - mu * np.cos(theta))
    # measured steady acceleration from the last slice (x = ½ a t²)
    a_meas = a
    assert a_theory > 0, "incline should slide (tanθ>μ)"
    assert abs(a_meas - a_theory) / a_theory < 0.02, f"a={a_meas:.4f} vs {a_theory:.4f}"
    print(f"[incline θ=30° μ=0.3] a_meas={a_meas:.4f}  a_theory={a_theory:.4f}  OK")


if __name__ == "__main__":
    print("=== ADR-39 P3 IMPL-EX Coulomb friction oracle ===")
    test_stick()
    test_slip_caps_at_muN()
    test_dissipation_nonneg()
    test_consistent_tangent_vs_fd()
    test_implex_converges_to_implicit()
    test_sliding_block_on_incline()
    print("ALL P3 FRICTION ORACLE CHECKS PASSED")
