"""
ADR 44 follow-up (-load nodal-force FRF/SSD/random) oracle — pins the modal-load
normalization BEFORE the C++ build.

For a nodal force pattern P (spatial shape) driven harmonically P e^{+iOm t},
the steady physical response is

    u(Om) = sum_a  psi_a * (psi_a^T P / m~_a) * H_a(Om),
    H_a   = 1/(w_a^2 - Om^2 + i Om d_a),
    m~_a  = psi_a^T M psi_a.

NORMALIZATION PIN (the load-bearing fact, from DomainModalProperties.cpp): the
`-unorm` scale factors are applied to the local V IN PLACE (DMP ~:636-655)
BEFORE the generalized masses (:771-772) and participation factors (:815) are
computed. So with psi_a = evec_a * Vscale_a (exactly what RSA/P1a/P2 recover
with), m~_a == mp.generalizedMasses()(a) DIRECTLY — no extra Vscale^2.

The formula is normalization-invariant by construction (psi psi^T/m~ is scale
free); this oracle exercises it with a DELIBERATELY non-mass-normalized basis
to mimic the OpenSees path, against the direct solve (K - Om^2 M + i Om C)^{-1} P.

Consistency identity also pinned here: -baseAccel -dir k is the special case
P = -M R_k (then psi^T P/m~ == -Gamma), so the two excitation channels must
agree exactly on the same model.

Random (-load) anchor: SDOF m q'' + c q' + k q = F(t) with one-sided-Hz white
force PSD G_F: |H_uF| = (1/m)|H|, so sigma_u^2 = G_F/(8 xi w^3 m^2).

Run:  python load_frf_oracle.py  -> all checks must print OK.
"""
import numpy as np

TWO_PI = 2.0 * np.pi


def modal_eig_unnormalized(M, K, seed=7):
    """Modes with a DELIBERATE per-mode random scaling (non-mass-normalized),
    mimicking arbitrary eigensolver + -unorm normalization."""
    Lam, Q = np.linalg.eigh(M)
    Mis = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    w2, U = np.linalg.eigh(Mis @ K @ Mis)
    Phi = Mis @ U                                   # mass-normalized
    rng = np.random.default_rng(seed)
    scales = rng.uniform(0.3, 3.0, Phi.shape[1])    # arbitrary re-scaling
    return np.sqrt(np.clip(w2, 0.0, None)), Phi * scales


def frf_modal_load(M, K, dfun, P, freqs, dof, resp="disp"):
    w, Psi = modal_eig_unnormalized(M, K)
    out = np.zeros(len(freqs), dtype=complex)
    gm = np.array([Psi[:, a] @ (M @ Psi[:, a]) for a in range(Psi.shape[1])])
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        uc = 0.0 + 0.0j
        for a in range(Psi.shape[1]):
            H = 1.0 / (w[a] ** 2 - Om ** 2 + 1j * Om * dfun(w[a]))
            uc += Psi[dof, a] * (Psi[:, a] @ P / gm[a]) * H
        if resp == "vel":
            uc *= 1j * Om
        elif resp == "accel":
            uc *= -(Om ** 2)
        out[i] = uc
    return out


def frf_direct(M, K, C, P, freqs, dof, resp="disp"):
    out = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        u = np.linalg.solve(K - Om ** 2 * M + 1j * Om * C, P)[dof]
        if resp == "vel":
            u *= 1j * Om
        elif resp == "accel":
            u *= -(Om ** 2)
        out[i] = u
    return out


def chain_MK(masses, ks):
    n = len(masses)
    M = np.diag(np.asarray(masses, float))
    K = np.zeros((n, n))
    for i, kk in enumerate(ks):
        lo = i - 1
        K[i, i] += kk
        if lo >= 0:
            K[lo, lo] += kk; K[lo, i] -= kk; K[i, lo] -= kk
    return M, K


def check_load_frf_vs_direct():
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    M, K = chain_MK(masses, ks)
    C = a0 * M + a1 * K
    P = np.array([0.0, 5.0, -2.0])                   # arbitrary multi-node load
    freqs = np.linspace(0.05, 6.0, 300)
    for resp in ("disp", "vel", "accel"):
        um = frf_modal_load(M, K, lambda w: a0 + a1 * w * w, P, freqs, 2, resp)
        ud = frf_direct(M, K, C, P, freqs, 2, resp)
        err = np.max(np.abs(um - ud)) / np.max(np.abs(ud))
        assert err < 1e-12, f"{resp}: modal-load FRF vs direct {err:.2e}"
        print(f"OK  1. load-FRF {resp:5s} vs direct: {err:.2e} "
              f"(non-mass-normalized basis)")

    # static limit
    u0 = frf_modal_load(M, K, lambda w: a0 + a1 * w * w, P, [1e-9], 2)[0]
    ustat = np.linalg.solve(K, P)[2]
    assert abs(u0.real - ustat) / abs(ustat) < 1e-9
    print(f"OK  2. static limit u(0) = K^-1 P  ({u0.real:.6e} vs {ustat:.6e})")


def check_baseaccel_equivalence():
    """-baseAccel -dir k == -load with P = -M R_k (unit amplitude)."""
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    xi = 0.04
    M, K = chain_MK(masses, ks)
    R = np.ones(3)
    P = -(M @ R)
    freqs = np.linspace(0.05, 6.0, 200)

    # base-accel channel in the P2 form: u = sum psi (-Gamma) H,  Gamma = psi^T M R/m~
    w, Psi = modal_eig_unnormalized(M, K)
    gm = np.array([Psi[:, a] @ (M @ Psi[:, a]) for a in range(Psi.shape[1])])
    ub = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        for a in range(Psi.shape[1]):
            H = 1.0 / (w[a] ** 2 - Om ** 2 + 1j * Om * 2 * xi * w[a])
            Gam = Psi[:, a] @ (M @ R) / gm[a]
            ub[i] += Psi[2, a] * (-Gam) * H

    ul = frf_modal_load(M, K, lambda w_: 2 * xi * w_, P, freqs, 2)
    err = np.max(np.abs(ub - ul)) / np.max(np.abs(ub))
    assert err < 1e-13, f"baseAccel != load(P=-MR): {err:.2e}"
    print(f"OK  3. -baseAccel == -load(P=-M R): {err:.2e}")


def check_force_psd_rms():
    """White force PSD through an SDOF: sigma_u^2 = G_F/(8 xi w^3 m^2)."""
    m, k, xi = 2.0, 800.0, 0.05
    w = np.sqrt(k / m)
    GF = 0.4
    fn = w / TWO_PI
    # dense biased grid over a wide band
    g = list(np.linspace(1e-4, 40 * fn, 6000))
    hw = 0.30 * fn
    g += [fn + hw * kk / 400 for kk in range(-400, 401)]
    freqs = np.unique(np.array(g))
    M1 = np.array([[m]]); K1 = np.array([[k]])
    H = frf_modal_load(M1, K1, lambda w_: 2 * xi * w_, np.ones(1), freqs, 0)
    sig = np.sqrt(np.trapezoid(np.abs(H) ** 2 * GF, freqs))
    ref = np.sqrt(GF / (8.0 * xi * w ** 3 * m * m))
    err = abs(sig - ref) / ref
    assert err < 2e-4, f"force-PSD RMS err {err:.2e}"
    print(f"OK  4. white force-PSD SDOF: sigma_u={sig:.6e} vs "
          f"sqrt(G_F/(8 xi w^3 m^2))={ref:.6e} ({err:.1e})")


if __name__ == "__main__":
    check_load_frf_vs_direct()
    check_baseaccel_equivalence()
    check_force_psd_rms()
    print("\nALL CHECKS PASS — modal-load normalization pinned:")
    print("  f_a = psi_a^T P / m~_a,  m~_a = generalizedMasses()(a),")
    print("  psi_a = evec_a * Vscale_a (the SAME pair RSA/P1a/P2 recover with);")
    print("  -baseAccel is the special case P = -M R.")
