"""
ADR 44 P2 (FRF / SSD) oracle — pins the modal frequency-response algebra, sign
convention, and normalization BEFORE the C++ build.

Modal transfer function (ADR 44 §4.4), in the (w, d=c/m) form that P1a already
uses and validated against direct Newmark:

    H_a(Om) = 1 / (w_a^2 - Om^2 + i*Om*d_a),    d_a = c_a/m_a = 2*xi_a*w_a.

For a uniform base acceleration u_g''(t) = e^{+i Om t} along influence vector R
(relative formulation, ADR 44 §4.3), the modal load amplitude is f_a = -Gamma_a
with Gamma_a = (phi_a^T M R)/m_a, and the physical *relative* displacement FRF is

    u(Om) = sum_a  psi_a * (-Gamma_a) * H_a(Om),     psi_a = phi_a * Vscale_a.

The DEFINITIVE end-to-end check: with all modes retained and CLASSICAL damping
C = a0 M + a1 K, this modal sum must equal the direct complex solve

    u_direct(Om) = (K - Om^2 M + i Om C)^{-1} (-M R).

(Same convention e^{+iOmt}; velocity FRF = iOm*u, relative-accel FRF = -Om^2*u.)

Sign convention PINNED here: e^{+i Om t}, so the +i*Om*d_a imaginary part gives
the physically correct -90 deg response phase lag at resonance.  A flipped sign
(e^{-iOmt}) would put +90 deg and is a classic silent FRF bug (ADR 44 §9).
"""
import numpy as np


def modal_eig(M, K):
    """Undamped modes, mass-normalized: phi^T M phi = I, K phi = M phi diag(w^2)."""
    # generalized symmetric eig via M^{-1/2}
    Lam, Q = np.linalg.eigh(M)
    Minvsqrt = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    A = Minvsqrt @ K @ Minvsqrt
    w2, U = np.linalg.eigh(A)
    Phi = Minvsqrt @ U                      # columns mass-normalized
    w = np.sqrt(np.clip(w2, 0.0, None))
    return w, Phi


def H_modal(w_a, d_a, Om):
    return 1.0 / (w_a * w_a - Om * Om + 1j * Om * d_a)


def frf_modal(M, K, a0, a1, R, Om, resp="disp"):
    """Physical relative-response FRF (all modes) at circular frequency Om."""
    w, Phi = modal_eig(M, K)
    n = M.shape[0]
    u = np.zeros(n, dtype=complex)
    for a in range(Phi.shape[1]):
        phi = Phi[:, a]
        wa = w[a]
        da = a0 + a1 * wa * wa               # = 2*xi*wa; rigid (wa=0) -> a0
        Gamma = phi @ (M @ R)                # m_a = 1 (mass-normalized)
        u += phi * (-Gamma) * H_modal(wa, da, Om)
    if resp == "disp":
        return u
    if resp == "vel":
        return 1j * Om * u
    if resp == "accel":                      # relative acceleration
        return -Om * Om * u
    raise ValueError(resp)


def frf_direct(M, K, a0, a1, R, Om, resp="disp"):
    C = a0 * M + a1 * K
    D = K - Om * Om * M + 1j * Om * C
    u = np.linalg.solve(D, -M @ R)
    if resp == "disp":
        return u
    if resp == "vel":
        return 1j * Om * u
    if resp == "accel":
        return -Om * Om * u
    raise ValueError(resp)


# ============================================================================
# 1) SDOF closed-form checks (gate a): resonant peak height + phase
# ============================================================================
print("[1] SDOF closed form")
m, k, xi = 2.0, 800.0, 0.05
w0 = np.sqrt(k / m)
M = np.array([[m]]); K = np.array([[k]]); R = np.array([1.0])
a1 = 0.0
a0 = 2.0 * xi * w0                            # pure mass-proportional -> xi at w0

# at resonance Om=w0: |u| = |H| = 1/(w0*d) with d=2*xi*w0  -> 1/(2*xi*w0^2)
u_res = frf_modal(M, K, a0, a1, R, w0)[0]
peak_analytic = 1.0 / (2.0 * xi * w0 * w0)
print(f"  |u(w0)| = {abs(u_res):.8e}  analytic 1/(2 xi w0^2) = {peak_analytic:.8e}")
assert abs(abs(u_res) - peak_analytic) / peak_analytic < 1e-10

# response phase at resonance: u = phi(-Gamma)/(i w0 d); phi=1/sqrt(m), Gamma=sqrt(m)
# u = (1/sqrt m)(-sqrt m)/(i w0 d) = -1/(i w0 d) = +i/(w0 d)  -> +90 deg
phase_deg = np.degrees(np.angle(u_res))
print(f"  angle(u(w0)) = {phase_deg:.4f} deg (expect +90)")
assert abs(phase_deg - 90.0) < 1e-8

# static limit Om->0: u(0) = -1/w0^2 = -m/k
u0 = frf_modal(M, K, a0, a1, R, 1e-9)[0]
print(f"  u(0) = {u0.real:.8e}  analytic -m/k = {-m/k:.8e}")
assert abs(u0.real - (-m / k)) < 1e-10 and abs(u0.imag) < 1e-12

# ============================================================================
# 2) MDOF modal-sum FRF == direct complex solve (THE end-to-end pin)
#    3-DOF shear chain, classical Rayleigh damping, full frequency sweep.
# ============================================================================
print("[2] MDOF modal FRF vs direct complex solve")
mm = np.array([2.0, 3.0, 1.5])
kk = np.array([1200.0, 900.0, 600.0])        # story stiffnesses
M = np.diag(mm)
K = np.zeros((3, 3))
for i in range(3):
    K[i, i] += kk[i]
    if i + 1 < 3:
        K[i, i] += kk[i + 1]
        K[i, i + 1] -= kk[i + 1]
        K[i + 1, i] -= kk[i + 1]
R = np.ones(3)
a0, a1 = 0.30, 2.0e-3

w, _ = modal_eig(M, K)
print("  modal freqs (Hz):", np.round(w / (2 * np.pi), 4))

worst = 0.0
for resp in ("disp", "vel", "accel"):
    for f in np.linspace(0.01, 5.0, 400):
        Om = 2 * np.pi * f
        um = frf_modal(M, K, a0, a1, R, Om, resp)
        ud = frf_direct(M, K, a0, a1, R, Om, resp)
        worst = max(worst, np.max(np.abs(um - ud)))
print(f"  worst |modal - direct| over disp/vel/accel sweep = {worst:.2e}")
assert worst < 1e-9, "MODAL-vs-DIRECT MISMATCH — sign/normalization bug"

# ============================================================================
# 3) FRF magnitude peaks at the modal frequencies (gate b)
# ============================================================================
print("[3] FRF peaks at modal freqs")
fgrid = np.linspace(0.01, 5.0, 20000)
mag = np.array([abs(frf_modal(M, K, a0, a1, R, 2 * np.pi * f)[0]) for f in fgrid])
# local maxima
peaks = [fgrid[i] for i in range(1, len(fgrid) - 1)
         if mag[i] > mag[i - 1] and mag[i] > mag[i + 1]]
modal_hz = np.sort(w / (2 * np.pi))
print("  detected peaks (Hz):", np.round(peaks, 4))
print("  modal freqs   (Hz):", np.round(modal_hz, 4))
# each detected damped peak should be within 1% of a modal frequency
for pf in peaks:
    assert min(abs(pf - modal_hz) / modal_hz) < 0.01

# ============================================================================
# 4) static limit for MDOF: u(0) == -K^{-1} M R
# ============================================================================
print("[4] MDOF static limit")
u0 = frf_modal(M, K, a0, a1, R, 1e-9)
ustat = np.linalg.solve(K, -M @ R)
print("  u(0)   =", np.round(u0.real, 8))
print("  static =", np.round(ustat, 8))
assert np.max(np.abs(u0.real - ustat)) < 1e-8
assert np.max(np.abs(u0.imag)) < 1e-9

print("\nALL OK — modal FRF algebra, sign convention (e^{+iOmt}), and "
      "normalization verified against the direct complex solve.")
