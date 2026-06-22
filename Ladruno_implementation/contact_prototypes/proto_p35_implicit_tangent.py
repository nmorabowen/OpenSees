"""ADR-39 P3.5 ORACLE — the ASSEMBLED implicit friction tangent (3D, global coords).

P3 shipped friction FORCE-only (explicit). P3.5 makes it work IMPLICITLY by assembling
the consistent friction tangent so Newton converges (a free tangential DOF is singular
without it). The per-traction blocks ∂tT/∂gT and ∂tT/∂gN are already validated in
proto_p3_friction.py (friction_tangent vs FD). Here we pin the 3D ASSEMBLED slave-block
tangent that the FE addKtToTang will build, against a finite difference of the friction
RESIDUAL — the new transcription risk.

Assembly (OpenSees convention K = −∂r/∂q, established in P2/P3):
  rel-disp operator   G       = [ I | −N_1 I | ... | −N_nps I ]   (3 × ndof)
  tangent projector   P_t     = I − n⊗n
  slip                gT      = P_t · G · u   (from engagement; gT0/gpT committed-frozen)
  normal gap          gN      = n · G · u
  friction force      r_fric  = −Gᵀ tT        (applied; slave += −tT, master_i += N_i tT)
  ASSEMBLED tangent   K_fric  = Gᵀ [ D_TT·P_t + d_TN⊗n ] G
                       with D_TT = ∂tT/∂gT (3×3), d_TN = ∂tT/∂gN (3-vec, slip only)

The SLAVE 3×3 block (G_slave = I) is  K_ss = D_TT·P_t + d_TN⊗n  — the dominant block and
the one carrying the NON-SYMMETRY (d_TN⊗n ≠ its transpose under slip). This oracle:
 (1) FD-validates K_ss against −∂(tFric)/∂u_slave with n FROZEN (the main-term tangent the
     FE assembles; ∂n/∂u is the deferred P2b-2c curvature block, excluded here);
 (2) confirms STICK ⇒ symmetric, SLIP ⇒ non-symmetric (needs FullGeneral/UmfPack);
 (3) checks the committed-N variant (freeze N ⇒ d_TN=0 ⇒ SYMMETRIC) — the robust knob.

Run: <pythoncore-3.12> proto_p35_implicit_tangent.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)


# ---------------------------------------------------------- traction + tangent (3D)
def tract_3d(gT, gpT, n, N, kt, mu):
    """3D Coulomb traction tT (global, tangent-plane vector). gT, gpT, n: 3-vectors.
    Mirrors frictionReturnMap but returns tT (not the negated force)."""
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT, float))
    cap = mu * N
    nrm = np.linalg.norm(tT_tr)
    if cap <= 0.0 or nrm <= cap:                 # stick (or inactive)
        return tT_tr
    return cap * (tT_tr / nrm)                    # slip: capped at the cone


def slave_block_analytic(gT, gpT, n, N, kn, kt, mu, committed_N=False):
    """K_ss = D_TT·P_t + d_TN⊗n  (3×3 slave block of the assembled friction tangent).
    committed_N: freeze N for the cone ⇒ d_TN=0 ⇒ symmetric (the robust variant)."""
    n = np.asarray(n, float)
    Pt = np.eye(3) - np.outer(n, n)              # tangent-plane projector
    tT_tr = kt * (np.asarray(gT, float) - np.asarray(gpT, float))
    cap = mu * N
    nrm = np.linalg.norm(tT_tr)
    if cap <= 0.0 or nrm <= cap:                 # STICK: D_TT = kt·P_t, d_TN = 0
        D_TT = kt * Pt
        d_TN = np.zeros(3)
    else:                                        # SLIP
        nh = tT_tr / nrm                         # slip direction (in tangent plane)
        D_TT = (cap * kt / nrm) * (Pt - np.outer(nh, nh))
        # ∂tT/∂gN = μ·(dN/dgN)·n̂ ; dN/dgN = −kn (N = kn·<−gN>) ⇒ d_TN = −μ·kn·n̂.
        # committed_N freezes the cap within the step ⇒ d_TN = 0 ⇒ symmetric tangent.
        d_TN = np.zeros(3) if committed_N else (-mu * kn * nh)
    # K_ss = D_TT·P_t + d_TN⊗n   (D_TT already carries P_t on the left for slip; the
    # ∂gT/∂u_s = P_t and ∂gN/∂u_s = nᵀ chain rule gives this exact form).
    return D_TT @ Pt + np.outer(d_TN, n)


# ----------------------------------------------------------------- tests
def _fd_slave_block(X_s, gpT, n, p0_master, kn, kt, mu, committed_N=False, h=1e-9):
    """FD of K_ss = −∂(tFric)/∂u_s with n FROZEN (flat master, ∂n/∂u excluded). The
    slave sits at X_s; the master is the flat plane through p0_master with normal n.
    gT = P_t·(x_s − p0) [engagement gT0=0 here]; gN = n·(x_s − p0)."""
    n = np.asarray(n, float)
    Pt = np.eye(3) - np.outer(n, n)
    N0 = kn * max(0.0, -(n @ (X_s - p0_master)))   # committed N for the frozen-cap variant

    def tfric(u):
        xs = X_s + u
        d = xs - p0_master
        gT = Pt @ d
        gN = n @ d
        N = N0 if committed_N else kn * max(0.0, -gN)
        return -tract_3d(gT, gpT, n, N, kt, mu)    # applied friction force on the slave

    base = tfric(np.zeros(3))
    K = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        K[:, j] = -(tfric(e) - base) / h           # K = −∂r/∂u
    return K


def test_slip_assembled_tangent_vs_fd():
    """SLIP regime: analytic K_ss == FD, and it is NON-SYMMETRIC (the d_TN⊗n column)."""
    n = np.array([0.0, 0.0, 1.0]); p0 = np.zeros(3)
    kn, kt, mu = 1e6, 1e5, 0.5
    X_s = np.array([0.03, 0.04, -1e-3])            # penetrated (gN<0) + tangential offset
    d = X_s - p0; Pt = np.eye(3) - np.outer(n, n)
    gT = Pt @ d; gN = n @ d; N = kn * max(0.0, -gN)
    # confirm we are in slip: ‖kt·gT‖ > μN
    assert kt * np.linalg.norm(gT) > mu * N, "test must be in the slip regime"
    Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu)
    Kfd = _fd_slave_block(X_s, np.zeros(3), n, p0, kn, kt, mu)
    err = np.max(np.abs(Ka - Kfd)) / kt
    asym = np.max(np.abs(Ka - Ka.T)) / kt
    assert err < 1e-3, f"slip tangent err {err:.3e}\nKa=\n{Ka}\nKfd=\n{Kfd}"
    assert asym > 1e-3, f"slip tangent should be NON-symmetric, asym={asym:.3e}"
    print(f"[slip] max|Ka−Kfd|/kt = {err:.2e}  asym/kt = {asym:.2e} (non-sym ✓)")


def test_stick_assembled_tangent_vs_fd():
    """STICK regime: analytic K_ss == FD == kt·P_t, and SYMMETRIC."""
    n = np.array([0.0, 0.0, 1.0]); p0 = np.zeros(3)
    kn, kt, mu = 1e6, 1e5, 0.5
    X_s = np.array([1e-3, 0.0, -1e-3])             # small tangential ⇒ stick (kt·gT<μN)
    d = X_s - p0; Pt = np.eye(3) - np.outer(n, n)
    gT = Pt @ d; gN = n @ d; N = kn * max(0.0, -gN)
    assert kt * np.linalg.norm(gT) < mu * N, "test must be in the stick regime"
    Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu)
    Kfd = _fd_slave_block(X_s, np.zeros(3), n, p0, kn, kt, mu)
    err = np.max(np.abs(Ka - Kfd)) / kt
    asym = np.max(np.abs(Ka - Ka.T)) / kt
    assert err < 1e-4, f"stick tangent err {err:.3e}"
    assert asym < 1e-10, f"stick tangent must be symmetric, asym={asym:.3e}"
    assert np.allclose(Ka, kt * Pt), f"stick tangent should be kt·P_t\n{Ka}"
    print(f"[stick] max|Ka−Kfd|/kt = {err:.2e}  symmetric ✓ (K=kt·P_t)")


def test_committed_N_is_symmetric():
    """The committed-N cap variant freezes N ⇒ d_TN=0 ⇒ the slip tangent is SYMMETRIC
    (the robust knob: works with a symmetric solver, at the cost of exact quadratic
    convergence when the normal force changes within a step). FD must match too."""
    n = np.array([0.0, 0.0, 1.0]); p0 = np.zeros(3)
    kn, kt, mu = 1e6, 1e5, 0.5
    X_s = np.array([0.03, 0.04, -1e-3])
    d = X_s - p0; Pt = np.eye(3) - np.outer(n, n)
    gT = Pt @ d; gN = n @ d; N = kn * max(0.0, -gN)
    Ka = slave_block_analytic(gT, np.zeros(3), n, N, kn, kt, mu, committed_N=True)
    Kfd = _fd_slave_block(X_s, np.zeros(3), n, p0, kn, kt, mu, committed_N=True)
    err = np.max(np.abs(Ka - Kfd)) / kt
    asym = np.max(np.abs(Ka - Ka.T)) / kt
    assert err < 1e-3, f"committed-N tangent err {err:.3e}"
    assert asym < 1e-9, f"committed-N tangent must be symmetric, asym={asym:.3e}"
    print(f"[committed-N] max|Ka−Kfd|/kt = {err:.2e}  symmetric ✓ (d_TN=0)")


def test_full_assembled_self_equilibrium():
    """The FULL ndof×ndof assembled K_fric = Gᵀ[D_TT·P_t+d_TN⊗n]G must be in self-
    equilibrium: a rigid tangential translation of the whole pair (slave + all master
    nodes by the same δ) produces ZERO net tangential force ⇒ rows/cols sum to 0 over
    the [slave, master_i] blocks weighted by [1, N_i]. (ΣN_i = 1 ⇒ G·1_rigid = 0.)"""
    n = np.array([0.0, 0.0, 1.0])
    kn, kt, mu = 1e6, 1e5, 0.5
    N = kn * 1e-3; gT = np.array([0.03, 0.04, 0.0]); gpT = np.zeros(3)
    Kss = slave_block_analytic(gT, gpT, n, N, kn, kt, mu)
    nps = 4
    Nsh = np.array([0.25, 0.25, 0.25, 0.25])       # shape weights at the quad center
    # G = [I | −N_i I]; K_fric = Gᵀ Kss G  (Kss is the 3×3 traction-block; G scatters it)
    G = np.zeros((3, 3 * (1 + nps)))
    G[:, 0:3] = np.eye(3)
    for i in range(nps):
        G[:, 3*(1+i):3*(2+i)] = -Nsh[i] * np.eye(3)
    Kfric = G.T @ Kss @ G
    # rigid tangential translation u = [δ ; δ ; ... ; δ] (slave + all masters by δ)
    delta = np.array([1.0, 1.0, 0.0])
    u_rigid = np.tile(delta, 1 + nps)
    f = Kfric @ u_rigid
    assert np.max(np.abs(f)) < 1e-6 * kt, f"rigid-translation force not zero: {f}"
    print(f"[self-eq] ‖K_fric·u_rigid‖ = {np.max(np.abs(f)):.2e} (≈0 ✓)")


if __name__ == "__main__":
    print("=== ADR-39 P3.5 implicit friction ASSEMBLED-tangent oracle ===")
    test_stick_assembled_tangent_vs_fd()
    test_slip_assembled_tangent_vs_fd()
    test_committed_N_is_symmetric()
    test_full_assembled_self_equilibrium()
    print("ALL P3.5 ASSEMBLED-TANGENT ORACLE CHECKS PASSED")
