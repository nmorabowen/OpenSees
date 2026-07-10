"""
ADR 44 P3 (randomResponse) oracle — pins the PSD convention, the RMS integral,
the spectral-moment statistics, and the factor-of-2 traps BEFORE the C++ build.

CONVENTION (pinned here, tested empirically in check 3):

  * The input PSD is the ONE-SIDED PSD of the base acceleration, G(f), as a
    function of frequency f in Hz, units (accel)^2/Hz:

        sigma_ug^2 = integral_0^inf G(f) df.

    This is the engineering-spec convention (wind / floor-vibration / equipment
    qualification curves are one-sided in Hz).  The classic silent bug is mixing
    it with the two-sided rad/s PSD S(Om) of the random-vibration textbooks:

        G(f) = 4*pi * S(Om = 2*pi*f)      (x2 one-sided, x2*pi rad/s -> Hz).

  * Response PSD through the (validated, ADR 44 P2) modal FRF H_x(f):

        G_xx(f) = |H_x(f)|^2 * G(f),

    where H_x is the P2 complex FRF of the requested response quantity
    (relative disp / vel / accel), i.e. exactly what `frequencyResponse` emits.

  * RMS = sqrt( integral_{fmin}^{fmax} G_xx(f) df ), trapezoidal on the sweep
    grid (the -biased grid resolves sharp resonant peaks).

  * Spectral moments (one-sided, Hz):  m_k = integral f^k G_xx(f) df.
    Zero-upcrossing rate  nu_0 = sqrt(m2/m0)  [Hz]  (narrow-band: nu_0 ~ f_n).
    Expected-peak factor over duration T (Davenport):
        g = sqrt(2 ln(nu0 T)) + 0.5772/sqrt(2 ln(nu0 T)),   E[peak] = g*sigma.

ANALYTIC ANCHOR (check 1): SDOF  x'' + d x' + w^2 x = -ug''(t) (relative form,
Gamma=1), white base-accel PSD.  With TWO-SIDED S0 in rad/s the textbook result
is sigma_x^2 = pi*S0/(d*w^2) = pi*S0/(2*xi*w^3).  Converting S0 = G0/(4*pi):

        sigma_x^2 = G0 / (8 * xi * w^3)          <-- the Hz/one-sided form
        sigma_v^2 = G0 / (8 * xi * w)            (velocity: extra (2*pi*f)^2 == w^2 at
                                                  the narrow-band peak; exact result
                                                  sigma_v = w*sigma_x for white noise)

Run:  python psd_rms_oracle.py   -> all checks must print OK.
"""
import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# modal machinery (identical to the P2 spike — the P3 layer rides on it)
# ---------------------------------------------------------------------------
def modal_eig(M, K):
    Lam, Q = np.linalg.eigh(M)
    Mis = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    w2, U = np.linalg.eigh(Mis @ K @ Mis)
    return np.sqrt(np.clip(w2, 0.0, None)), Mis @ U


def frf_modal(M, K, dfun, R, freqs, dof, resp="disp"):
    """P2 modal-sum complex FRF (relative response to unit base accel)."""
    w, Phi = modal_eig(M, K)
    out = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        uc = 0.0 + 0.0j
        for a in range(Phi.shape[1]):
            Gam = Phi[:, a] @ (M @ R)
            H = 1.0 / (w[a] ** 2 - Om ** 2 + 1j * Om * dfun(w[a]))
            uc += Phi[dof, a] * (-Gam) * H
        if resp == "vel":
            uc *= 1j * Om
        elif resp == "accel":
            uc *= -(Om ** 2)
        out[i] = uc
    return out


def frf_direct(M, K, C, R, freqs, dof, resp="disp"):
    """Independent direct complex solve (K - Om^2 M + i Om C)^{-1}(-M R)."""
    out = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        u = np.linalg.solve(K - Om ** 2 * M + 1j * Om * C, -(M @ R))[dof]
        if resp == "vel":
            u *= 1j * Om
        elif resp == "accel":
            u *= -(Om ** 2)
        out[i] = u
    return out


def biased_grid(fmin, fmax, nf, f_modes, nclust=15, halfw_frac=0.05):
    """The P2 -biased grid: LIN base + +/-5% cluster around each in-band mode."""
    g = list(np.linspace(fmin, fmax, nf))
    for fa in f_modes:
        if fmin <= fa <= fmax:
            hw = halfw_frac * fa
            g += [fa + hw * k / nclust for k in range(-nclust, nclust + 1)
                  if fmin <= fa + hw * k / nclust <= fmax]
    return np.unique(np.array(g))


def rms_from_psd(freqs, Gxx):
    return np.sqrt(np.trapezoid(Gxx, freqs))


def moments(freqs, Gxx):
    m0 = np.trapezoid(Gxx, freqs)
    m2 = np.trapezoid(freqs ** 2 * Gxx, freqs)
    return m0, m2


# ---------------------------------------------------------------------------
# 1. analytic anchor: white-noise SDOF, one-sided-Hz form
# ---------------------------------------------------------------------------
def check_white_noise_sdof():
    m, k, xi = 2.0, 800.0, 0.05
    w = np.sqrt(k / m)
    fn = w / TWO_PI
    G0 = 0.13  # (accel)^2/Hz, flat

    M = np.array([[m]]); K = np.array([[k]]); R = np.ones(1)
    f_modes = [fn]
    # band wide enough that the tail contributes ~nothing (integrand ~ f^-4)
    freqs = biased_grid(1e-4, 40.0 * fn, 6000, f_modes, nclust=400, halfw_frac=0.30)
    H = frf_modal(M, K, lambda wa: 2.0 * xi * wa, R, freqs, dof=0)
    sig = rms_from_psd(freqs, np.abs(H) ** 2 * G0)

    sig_analytic = np.sqrt(G0 / (8.0 * xi * w ** 3))
    err = abs(sig - sig_analytic) / sig_analytic
    assert err < 2e-4, f"white-noise SDOF RMS err {err:.2e}"
    print(f"OK  1. white-noise SDOF: sigma={sig:.6e} vs analytic "
          f"sqrt(G0/(8 xi w^3))={sig_analytic:.6e}  (rel err {err:.1e})")

    # velocity: sigma_v = w * sigma_x for white noise (exact)
    Hv = frf_modal(M, K, lambda wa: 2.0 * xi * wa, R, freqs, dof=0, resp="vel")
    sig_v = rms_from_psd(freqs, np.abs(Hv) ** 2 * G0)
    err_v = abs(sig_v - w * sig_analytic) / (w * sig_analytic)
    assert err_v < 2e-3, f"white-noise SDOF vel RMS err {err_v:.2e}"
    print(f"OK  1b. velocity RMS = w*sigma_x to {err_v:.1e}")


# ---------------------------------------------------------------------------
# 2. MDOF: modal-|H|^2 PSD propagation == direct-solve PSD propagation
# ---------------------------------------------------------------------------
def check_mdof_vs_direct():
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    n = len(masses)
    M = np.diag(masses)
    K = np.zeros((n, n))
    for i, kk in enumerate(ks):
        lo = i - 1
        K[i, i] += kk
        if lo >= 0:
            K[lo, lo] += kk; K[lo, i] -= kk; K[i, lo] -= kk
    C = a0 * M + a1 * K
    R = np.ones(n)
    w, _ = modal_eig(M, K)
    f_modes = w / TWO_PI

    # a NON-flat input PSD (band-limited, sloped) to exercise G(f) sampling
    def Gin(f):
        return 0.05 * (1.0 + 0.3 * f) * ((f > 0.2) & (f < 6.0))

    freqs = biased_grid(0.05, 8.0, 800, f_modes, nclust=60, halfw_frac=0.10)
    Gf = Gin(freqs)
    for resp in ("disp", "vel", "accel"):
        # Rayleigh damping law in (w, d) form: d(w) = a0 + a1 w^2
        Hm = frf_modal(M, K, lambda wa: a0 + a1 * wa ** 2, R, freqs, 2, resp)
        Hd = frf_direct(M, K, C, R, freqs, 2, resp)
        s_m = rms_from_psd(freqs, np.abs(Hm) ** 2 * Gf)
        s_d = rms_from_psd(freqs, np.abs(Hd) ** 2 * Gf)
        err = abs(s_m - s_d) / s_d
        assert err < 1e-10, f"{resp}: modal vs direct RMS err {err:.2e}"
        print(f"OK  2. MDOF {resp:5s}: RMS modal {s_m:.6e} == direct {s_d:.6e} "
              f"({err:.1e})")


# ---------------------------------------------------------------------------
# 3. THE factor-of-2 empirical pin: Monte-Carlo synthetic time history
# ---------------------------------------------------------------------------
def check_monte_carlo():
    """Generate a stationary realization of the target ONE-SIDED-Hz PSD, run it
    through the exact SDOF integrator, and compare the time-domain RMS to the
    spectral prediction.  This is the check that catches a wrong one-sided /
    two-sided or Hz/rad-per-s factor (they'd show up as sqrt(2) or sqrt(2*pi))."""
    rng = np.random.default_rng(12345)
    m, k, xi = 2.0, 800.0, 0.03
    w = np.sqrt(k / m); fn = w / TWO_PI; d = 2.0 * xi * w

    G0, f1, f2 = 0.02, 0.2, 12.0          # band-limited white, one-sided Hz
    # synthesis: ug''(t) = sum_k sqrt(2 G0 df) cos(2 pi f_k t + phi_k)
    df = 0.001953125                        # dyadic: fits an exact FFT length
    fk = np.arange(f1, f2, df) + 0.5 * df
    amp = np.sqrt(2.0 * G0 * df)
    phi = rng.uniform(0.0, TWO_PI, len(fk))

    dt = 0.01
    T = 1.0 / df                            # one full period of the synthesis
    t = np.arange(0.0, T, dt)
    ug = (amp * np.cos(TWO_PI * np.outer(t, fk) + phi)).sum(axis=1)

    # sanity: input variance == integral of the one-sided PSD == G0*(f2-f1)
    var_in_target = G0 * (f2 - f1)
    var_in = ug.var()
    err_in = abs(var_in - var_in_target) / var_in_target
    assert err_in < 5e-3, f"synthetic input variance err {err_in:.2e}"

    # exact PWL SDOF march (the P1a recurrence, trusted at 1e-14)
    from numpy.linalg import solve as lsolve
    A = np.array([[0.0, 1.0], [-w * w, -d]])
    import scipy.linalg as sla
    E = sla.expm(A * dt)
    # first-order-hold blocks via the van-Loan augmented matrix
    Z = np.zeros((2, 2))
    Aug = np.block([[A, np.eye(2), Z],
                    [Z, Z, np.eye(2)],
                    [Z, Z, Z]])
    Eaug = sla.expm(Aug * dt)
    P1 = Eaug[:2, 2:4]; P2 = Eaug[:2, 4:6]
    B0 = (P1 - P2 / dt) @ np.array([0.0, 1.0])
    B1 = (P2 / dt) @ np.array([0.0, 1.0])

    x = np.zeros(2)
    q = np.empty(len(t)); qd = np.empty(len(t))
    fload = -ug                              # Gamma = 1 for the SDOF
    for i in range(len(t)):
        q[i], qd[i] = x
        if i + 1 < len(t):
            x = E @ x + B0 * fload[i] + B1 * fload[i + 1]

    # discard the start-up transient (~10 damped periods)
    n0 = int(10.0 / (xi * fn) / dt)
    sig_mc = q[n0:].std()
    sigv_mc = qd[n0:].std()

    # spectral prediction on a biased grid over the excitation band
    freqs = biased_grid(f1, f2, 3000, [fn], nclust=500, halfw_frac=0.25)
    M1 = np.array([[m]]); K1 = np.array([[k]])
    H = frf_modal(M1, K1, lambda wa: 2.0 * xi * wa, np.ones(1), freqs, 0)
    sig_sp = rms_from_psd(freqs, np.abs(H) ** 2 * G0)
    Hv = frf_modal(M1, K1, lambda wa: 2.0 * xi * wa, np.ones(1), freqs, 0, "vel")
    sigv_sp = rms_from_psd(freqs, np.abs(Hv) ** 2 * G0)

    err = abs(sig_mc - sig_sp) / sig_sp
    errv = abs(sigv_mc - sigv_sp) / sigv_sp
    # one realization of a narrow-band process: sampling error ~ few %
    assert err < 0.05, f"MC disp RMS err {err:.2e} (factor-of-2 bug?)"
    assert errv < 0.05, f"MC vel RMS err {errv:.2e}"
    print(f"OK  3. Monte-Carlo: disp RMS time-domain {sig_mc:.5e} vs spectral "
          f"{sig_sp:.5e} ({err * 100:.1f}%); vel {errv * 100:.1f}%")
    # a wrong one-sided/two-sided factor would give err ~ 41% (sqrt(2)); a
    # Hz/rad mixup ~ 150% (sqrt(2 pi)); both far above the 5% gate.

    # 3b. zero-upcrossing rate: count actual upcrossings vs nu0 = sqrt(m2/m0)
    m0, m2 = moments(freqs, np.abs(H) ** 2 * G0)
    nu0 = np.sqrt(m2 / m0)
    qz = q[n0:]
    ups = np.sum((qz[:-1] < 0.0) & (qz[1:] >= 0.0))
    nu0_mc = ups / (len(qz) * dt)
    err_nu = abs(nu0 - nu0_mc) / nu0_mc
    assert err_nu < 0.05, f"nu0 err {err_nu:.2e}"
    assert abs(nu0 - fn) / fn < 0.02, "narrow-band nu0 should sit near fn"
    print(f"OK  3b. nu0 spectral {nu0:.4f} Hz vs counted {nu0_mc:.4f} Hz "
          f"({err_nu * 100:.1f}%); fn={fn:.4f}")

    # 3c. Davenport expected peak vs the realized max (advisory, wide tol)
    Tdur = len(qz) * dt
    lnt = 2.0 * np.log(nu0 * Tdur)
    g = np.sqrt(lnt) + 0.5772 / np.sqrt(lnt)
    peak_pred = g * sig_sp
    peak_mc = np.abs(qz).max()
    print(f"ADV 3c. Davenport peak {peak_pred:.5e} vs realized {peak_mc:.5e} "
          f"(ratio {peak_mc / peak_pred:.3f}; one realization, no assert)")


# ---------------------------------------------------------------------------
# 4. grid-resolution honesty: a plain LIN grid under-integrates a sharp peak
# ---------------------------------------------------------------------------
def check_grid_resolution():
    m, k, xi = 2.0, 800.0, 0.005          # sharp: half-power width ~2 xi fn
    w = np.sqrt(k / m); fn = w / TWO_PI
    G0 = 0.05
    M1 = np.array([[m]]); K1 = np.array([[k]]); R = np.ones(1)

    ref = np.sqrt(G0 / (8.0 * xi * w ** 3))

    lin = np.linspace(0.1, 10.0, 200)
    Hl = frf_modal(M1, K1, lambda wa: 2 * xi * wa, R, lin, 0)
    sig_lin = rms_from_psd(lin, np.abs(Hl) ** 2 * G0)

    bia = biased_grid(0.1, 10.0, 200, [fn], nclust=40, halfw_frac=0.05)
    Hb = frf_modal(M1, K1, lambda wa: 2 * xi * wa, R, bia, 0)
    sig_bia = rms_from_psd(bia, np.abs(Hb) ** 2 * G0)

    err_lin = abs(sig_lin - ref) / ref
    err_bia = abs(sig_bia - ref) / ref
    assert err_bia < 0.02, f"biased grid err {err_bia:.2e}"
    print(f"OK  4. sharp peak (xi=0.5%): LIN-200 err {err_lin * 100:.0f}% vs "
          f"BIASED-200 err {err_bia * 100:.2f}%  -> -biased is the documented "
          f"default recommendation")


if __name__ == "__main__":
    check_white_noise_sdof()
    check_mdof_vs_direct()
    check_monte_carlo()
    check_grid_resolution()
    print("\nALL CHECKS PASS — conventions pinned:")
    print("  input  = ONE-SIDED base-accel PSD G(f), Hz;  sigma_ug^2 = int G df")
    print("  output = G_xx(f) = |H_x(f)|^2 G(f);  RMS = sqrt(trapz);")
    print("  nu0 = sqrt(m2/m0) Hz;  Davenport peak factor with -duration.")
