"""Reference (oracle) for the LadrunoHS / LadrunoHSSmall double-hardening soil model.

Pure-numpy reference implementation of SRC/material/nD/LadrunoHSKernel.h:
  * Schanz, Vermeer & Bonnier (1999) Hardening Soil (HS) 4-surface model:
      - Surface 1 (f^s): Hyperbolic shear-hardening cone with Rowe (1962)
        stress-dilatancy and optional e >= e_max dilatancy cutoff.
      - Surface 2 (f^f): Ultimate failure truncation at q_tilde = M_c * (p + a).
      - Surface 3 (f^c): Elliptical volumetric cap with associated flow in
        (p, q_tilde) space and exact constructor auto-calibration of (alpha, beta)
        from (E_oed^ref, K_0^nc).
      - Surface 4 (f^t): Rankine/apex tension cutoff at p + a >= p_cut.
  * Dual Lode-angle formulations (-lodeShape mc | mn):
      - 'mc': Mohr-Coulomb hexagon (Schanz et al. 1999 / PLAXIS HS / FLAC PH)
        with C^infinity endpoint-preserving corner smoothing at sin(3*theta) = +-1.
      - 'mn': Panteghini & Lagioia (2014) strictly convex Matsuoka-Nakai locus
        (I1*I2/I3 = m_MN), coinciding identically with 'mc' at triaxial
        compression (theta = +pi/6) and triaxial extension (theta = -pi/6) while
        strictly circumscribing 'mc' at intermediate Lode angles (PlaneStrain).
  * Benz (2007) / Benz, Vermeer & Schwab (2009) HSsmall overlay (-smallStrain):
      - 6-component deviatoric strain-history tensor H updated via spectral
        projection in the principal frame of d(e), with exact step-integrated
        Hardin-Drnevich secant/tangent modulus G_t(gamma_hist) in [G_ur, G_0]
        and Masing's 2nd rule upon strain reversal.

Sign convention inside the kernel:
  * Compression-positive effective stress (p = tr(sigma)/3 > 0 in compression)
    and compression-positive strain (eps_v = tr(eps) > 0 in compression).
  * External OpenSees interface negates stress and strain at the wrapper boundary
    (leaving the 6x6 tangent matrix C^alg unchanged in sign).
  * Voigt order (6-component): [11, 22, 33, 12, 23, 31], with engineering shear
    strains [eps_11, eps_22, eps_33, gamma_12, gamma_23, gamma_31].

No OpenSees / openseespy dependency — numpy only.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Tuple

import numpy as np


# --------------------------------------------------------------------------- #
#  Voigt <-> 3x3 symmetric tensor helpers (OpenSees order: 11,22,33,12,23,31) #
# --------------------------------------------------------------------------- #
_VOIGT_I = (0, 1, 2, 0, 1, 2)
_VOIGT_J = (0, 1, 2, 1, 2, 0)
_I3 = np.eye(3, dtype=float)
_I6_M = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=float)
_IDEV6 = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]) - np.outer(_I6_M, _I6_M) / 3.0


def stress_voigt_to_tensor(sig6: np.ndarray) -> np.ndarray:
    """Convert 6-component stress [s11, s22, s33, s12, s23, s31] to 3x3 matrix."""
    T = np.zeros((3, 3), dtype=float)
    for k in range(6):
        i, j = _VOIGT_I[k], _VOIGT_J[k]
        T[i, j] = sig6[k]
        T[j, i] = sig6[k]
    return T


def stress_tensor_to_voigt(T: np.ndarray) -> np.ndarray:
    """Convert symmetric 3x3 stress tensor to [s11, s22, s33, s12, s23, s31]."""
    return np.array([T[_VOIGT_I[k], _VOIGT_J[k]] for k in range(6)], dtype=float)


def strain_voigt_to_tensor(eps6: np.ndarray) -> np.ndarray:
    """Convert 6-component engineering strain [e11, e22, e33, g12, g23, g31] to 3x3 tensor."""
    T = np.zeros((3, 3), dtype=float)
    for k in range(3):
        T[k, k] = eps6[k]
    for k in range(3, 6):
        i, j = _VOIGT_I[k], _VOIGT_J[k]
        T[i, j] = 0.5 * eps6[k]
        T[j, i] = 0.5 * eps6[k]
    return T


def strain_tensor_to_voigt(T: np.ndarray) -> np.ndarray:
    """Convert symmetric 3x3 strain tensor to engineering Voigt [e11, e22, e33, g12, g23, g31]."""
    return np.array(
        [T[0, 0], T[1, 1], T[2, 2], 2.0 * T[0, 1], 2.0 * T[1, 2], 2.0 * T[2, 0]],
        dtype=float,
    )


# --------------------------------------------------------------------------- #
#  Parameters & Constructor Auto-Calibration of (alpha, beta_cap)             #
# --------------------------------------------------------------------------- #
@dataclass
class HSParams:
    E50_ref: float
    Eoed_ref: float
    Eur_ref: float
    m: float
    c: float
    phi_deg: float
    psi_deg: float
    nu_ur: float = 0.20
    p_ref: float = 100.0
    Rf: float = 0.90
    K0nc: float = -1.0          # <= 0 means auto: 1 - sin(phi)
    sigma_t: float = 0.0        # Rankine tensile strength >= 0
    OCR: float = 1.0
    pop: float = 0.0            # Pre-overburden pressure >= 0
    e0: float = 0.50
    emax: float = -1.0          # <= 0 means dilatancy cutoff inactive
    use_small_strain: bool = False
    G0_ref: float = 0.0
    gamma_07: float = 1.0e-4
    lode_shape: str = "mc"      # "mc" (Mohr-Coulomb) or "mn" (Matsuoka-Nakai)
    stiff_stress: str = "sigma3"  # "sigma3" (PLAXIS HS) or "p" (mean stress GHS)
    alpha_cap: float = -1.0     # <= 0 means auto-calibrate from (Eoed_ref, K0nc)
    beta_cap: float = -1.0      # <= 0 means auto-calibrate from (Eoed_ref, K0nc)
    p_min_ratio: float = 1.0e-3

    # Derived constants populated in __post_init__
    phi: float = field(init=False)
    psi: float = field(init=False)
    sin_phi: float = field(init=False)
    cos_phi: float = field(init=False)
    sin_psi: float = field(init=False)
    sin_phicv: float = field(init=False)
    a_apex: float = field(init=False)
    p_hat_ref: float = field(init=False)
    p_min: float = field(init=False)
    M_c: float = field(init=False)
    delta_lode: float = field(init=False)
    gamma_mn: float = field(init=False)
    c0_mn: float = field(init=False)
    gamma_mc: float = field(init=False)
    th_max_mc: float = field(init=False)
    Ei_ref: float = field(init=False)
    Gur_ref: float = field(init=False)
    Kur_ref: float = field(init=False)

    def __post_init__(self) -> None:
        if self.E50_ref <= 0.0 or self.Eoed_ref <= 0.0 or self.Eur_ref <= 0.0:
            raise ValueError("Reference moduli E50_ref, Eoed_ref, Eur_ref must be > 0.")
        if not (0.0 <= self.m <= 1.0):
            raise ValueError("Power exponent m must be in [0, 1].")
        if not (0.0 < self.phi_deg < 89.0):
            raise ValueError("Friction angle phi_deg must be in (0, 89) deg.")
        if not (0.0 <= self.psi_deg <= self.phi_deg + 1.0e-12):
            raise ValueError("Dilatancy angle psi_deg must be in [0, phi_deg].")
        if not (-1.0 < self.nu_ur < 0.495):
            raise ValueError("Poisson ratio nu_ur must be in (-1, 0.495).")
        if not (0.0 < self.Rf < 1.0):
            raise ValueError("Failure ratio Rf must be in (0, 1).")

        self.phi = math.radians(self.phi_deg)
        self.psi = math.radians(self.psi_deg)
        self.sin_phi = math.sin(self.phi)
        self.cos_phi = math.cos(self.phi)
        self.sin_psi = math.sin(self.psi)
        self.sin_phicv = max(
            0.0, (self.sin_phi - self.sin_psi) / (1.0 - self.sin_phi * self.sin_psi)
        )
        if self.K0nc <= 0.0:
            self.K0nc = 1.0 - self.sin_phi
        self.a_apex = self.c * (self.cos_phi / self.sin_phi) if self.sin_phi > 1e-12 else 0.0
        self.p_hat_ref = self.p_ref + self.a_apex
        self.p_min = max(1.0e-6, self.p_min_ratio * self.p_ref)

        self.M_c = 6.0 * self.sin_phi / (3.0 - self.sin_phi)
        self.delta_lode = (3.0 + self.sin_phi) / (3.0 - self.sin_phi)

        # Panteghini & Lagioia (2014) Eq. (18) exact Matsuoka-Nakai constants:
        s2 = self.sin_phi * self.sin_phi
        self.gamma_mn = (9.0 - s2) * self.sin_phi / ((3.0 + s2) ** 1.5)
        self.c0_mn = math.cos(math.acos(-self.gamma_mn) / 3.0)

        # Endpoint-preserving C^infinity Mohr-Coulomb smoothing constant:
        self.gamma_mc = 1.0 - 1.0e-7
        self.th_max_mc = math.asin(self.gamma_mc) / 3.0

        self.Ei_ref = 2.0 * self.E50_ref / (2.0 - self.Rf)
        if self.Eur_ref < 1.001 * self.Ei_ref:
            raise ValueError(
                f"Eur_ref ({self.Eur_ref}) must exceed Ei_ref = 2*E50_ref/(2-Rf) ({self.Ei_ref})."
            )
        self.Gur_ref = self.Eur_ref / (2.0 * (1.0 + self.nu_ur))
        self.Kur_ref = self.Eur_ref / (3.0 * (1.0 - 2.0 * self.nu_ur))

        self.lode_shape = self.lode_shape.lower()
        if self.lode_shape not in ("mc", "mn"):
            raise ValueError(f"Unknown lode_shape '{self.lode_shape}' (expected 'mc' or 'mn').")
        self.stiff_stress = self.stiff_stress.lower()
        if self.stiff_stress not in ("sigma3", "p"):
            raise ValueError(f"Unknown stiff_stress '{self.stiff_stress}' (expected 'sigma3' or 'p').")

        if self.alpha_cap <= 0.0 or self.beta_cap <= 0.0:
            alpha_cal, beta_cal = calibrate_cap_parameters(self)
            if self.alpha_cap <= 0.0:
                self.alpha_cap = alpha_cal
            if self.beta_cap <= 0.0:
                self.beta_cap = beta_cal


def calibrate_cap_parameters(p: HSParams) -> Tuple[float, float]:
    """Closed-form calibration of (alpha_cap, beta_cap) from (Eoed_ref, K0nc).

    Guarantees that along virgin 1D oedometric compression (d_eps2 = d_eps3 = 0)
    on the K0nc normal-consolidation ray (hat_sig3 = K0nc * hat_sig1):
      1. Lateral stress ratio remains hat_sig3 / hat_sig1 = K0nc at all stresses, and
      2. Axial tangent stiffness satisfies d(sig1)/d(eps1) = Eoed_ref * (hat_sig1 / p_hat_ref)^m.
    """
    K0 = p.K0nc
    kp = (1.0 + 2.0 * K0) / 3.0
    kq = 1.0 - K0
    eta_K0 = kq / kp

    r_K0 = p.Rf * ((1.0 - K0) / K0) * ((1.0 - p.sin_phi) / (2.0 * p.sin_phi))
    r_K0 = min(r_K0, 0.99)

    s_phim_K0 = (1.0 - K0) / (1.0 + K0)
    s_psim_K0 = (s_phim_K0 - p.sin_phicv) / (1.0 - s_phim_K0 * p.sin_phicv)

    s_K0 = (K0 ** p.m) if p.stiff_stress == "sigma3" else (kp ** p.m)

    # Rate of plastic shear strain d(gamma^p)/d(sig1) along the virgin K0nc ray:
    d_gam_p = (1.0 - p.m) * kq * (
        2.0 / (p.Ei_ref * s_K0 * (1.0 - r_K0)) - 2.0 / (p.Eur_ref * s_K0)
    )

    Av = 1.0 / p.Eoed_ref - kp / (p.Kur_ref * s_K0) + s_psim_K0 * d_gam_p
    Aq = 2.0 / (3.0 * p.Eoed_ref) - kq / (3.0 * p.Gur_ref * s_K0) - 0.5 * d_gam_p

    # Guard admissible positive cap compliances:
    Av_eff = max(Av, 1.0e-3 / p.Eoed_ref)
    Aq_eff = max(Aq, 1.0e-3 / p.Eoed_ref)

    alpha = math.sqrt(eta_K0 * Av_eff / Aq_eff)
    alpha = max(0.20, min(5.0, alpha))
    kpp = kp * math.sqrt(1.0 + (eta_K0 / alpha) ** 2)
    beta_cap = Av_eff * (kpp ** (p.m - 1.0))
    return alpha, beta_cap


# --------------------------------------------------------------------------- #
#  Lode Shape Function Theta(x) and its derivative dTheta/dx                  #
#  where x = sin(3 * theta_c) = 13.5 * J3 / q^3 in [-1, +1]                   #
# --------------------------------------------------------------------------- #
def lode_theta_and_deriv(p: HSParams, x: float) -> Tuple[float, float]:
    """Return (Theta(x), dTheta/dx) normalized so Theta(+1)=1 (TXC) and Theta(-1)=delta (TXE)."""
    x_clamped = max(-1.0, min(1.0, float(x)))
    if p.lode_shape == "mn":
        # Panteghini & Lagioia (2014) Eq. (17)-(19):
        arg = -p.gamma_mn * x_clamped
        psi = math.acos(arg) / 3.0
        Theta = math.cos(psi) / p.c0_mn
        denom = math.sqrt(max(1.0e-30, 1.0 - arg * arg))
        dpsi_dx = p.gamma_mn / (3.0 * denom)
        dTheta_dx = -math.sin(psi) * dpsi_dx / p.c0_mn
        return Theta, dTheta_dx
    else:
        # Endpoint-exact C^infinity Mohr-Coulomb hexagon:
        arg = p.gamma_mc * x_clamped
        th = math.asin(arg) / 3.0
        c_max = math.cos(p.th_max_mc)
        s_max = math.sin(p.th_max_mc)
        A_c = 3.0 / c_max
        B_s = p.sin_phi / s_max
        inv_denom = 1.0 / (3.0 - p.sin_phi)
        Theta = (A_c * math.cos(th) - B_s * math.sin(th)) * inv_denom
        dth_dx = p.gamma_mc / (3.0 * math.sqrt(max(1.0e-30, 1.0 - arg * arg)))
        dTheta_dx = (-A_c * math.sin(th) - B_s * math.cos(th)) * inv_denom * dth_dx
        return Theta, dTheta_dx


def invariants_and_lode_grad(
    p: HSParams, sig6: np.ndarray
) -> Tuple[float, float, float, float, np.ndarray, np.ndarray]:
    """Compute (p_mean, q, x_sin3, Theta, dq_dsig6, dTheta_dsig6) in Voigt space.

    Note: gradient wrt Voigt stress [s11, s22, s33, s12, s23, s31] includes the
    factor of 2 on shear components so d(f)/d(sig6) . d(sig6) = df.
    """
    p_mean = (sig6[0] + sig6[1] + sig6[2]) / 3.0
    s0 = sig6[0] - p_mean
    s1 = sig6[1] - p_mean
    s2 = sig6[2] - p_mean
    s3, s4, s5 = sig6[3], sig6[4], sig6[5]

    J2 = 0.5 * (s0 * s0 + s1 * s1 + s2 * s2) + s3 * s3 + s4 * s4 + s5 * s5
    q = math.sqrt(3.0 * J2)
    if q < 1.0e-12:
        return p_mean, 0.0, 1.0, 1.0, np.zeros(6), np.zeros(6)

    # Voigt gradient of q = sqrt(3*J2): dq/dsig = (3/(2*q)) * s_voigt_weight
    s_vec_w = np.array([s0, s1, s2, 2.0 * s3, 2.0 * s4, 2.0 * s5], dtype=float)
    dq_dsig = (1.5 / q) * s_vec_w

    # J3 = det(s) and its Voigt gradient dJ3/dsig = dev(s^2)
    J3 = (
        s0 * s1 * s2
        + 2.0 * s3 * s4 * s5
        - s0 * s4 * s4
        - s1 * s5 * s5
        - s2 * s3 * s3
    )
    # s^2 components:
    t00 = s0 * s0 + s3 * s3 + s5 * s5
    t11 = s3 * s3 + s1 * s1 + s4 * s4
    t22 = s5 * s5 + s4 * s4 + s2 * s2
    t01 = s0 * s3 + s3 * s1 + s5 * s4
    t12 = s3 * s5 + s1 * s4 + s4 * s2
    t20 = s5 * s0 + s4 * s3 + s2 * s5
    tr_t = (t00 + t11 + t22) / 3.0  # == (2/3) * J2
    dJ3_dsig = np.array(
        [
            t00 - tr_t,
            t11 - tr_t,
            t22 - tr_t,
            2.0 * t01,
            2.0 * t12,
            2.0 * t20,
        ],
        dtype=float,
    )

    q3 = q * q * q
    x_raw = 13.5 * J3 / q3
    x_sin3 = max(-1.0, min(1.0, x_raw))
    Theta, dTheta_dx = lode_theta_and_deriv(p, x_sin3)

    # dx/dsig = (13.5 / q^3) * dJ3/dsig - (3 * x / q) * dq/dsig
    dx_dsig = (13.5 / q3) * dJ3_dsig - (3.0 * x_sin3 / q) * dq_dsig
    dTheta_dsig = dTheta_dx * dx_dsig
    return p_mean, q, x_sin3, Theta, dq_dsig, dTheta_dsig


# --------------------------------------------------------------------------- #
#  Internal State                                                             #
# --------------------------------------------------------------------------- #
@dataclass
class HSState:
    sig: np.ndarray = field(default_factory=lambda: np.zeros(6, dtype=float))
    eps: np.ndarray = field(default_factory=lambda: np.zeros(6, dtype=float))
    gamma_p: float = 0.0
    p_p: float = 100.0
    void_ratio: float = 0.50
    # HSsmall deviatoric strain history state:
    e_dev: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=float))
    H_hist: np.ndarray = field(default_factory=lambda: np.zeros((3, 3), dtype=float))
    is_virgin: bool = True

    def copy(self) -> "HSState":
        return HSState(
            sig=self.sig.copy(),
            eps=self.eps.copy(),
            gamma_p=float(self.gamma_p),
            p_p=float(self.p_p),
            void_ratio=float(self.void_ratio),
            e_dev=self.e_dev.copy(),
            H_hist=self.H_hist.copy(),
            is_virgin=bool(self.is_virgin),
        )


def initialize_state_from_stress(
    p: HSParams, sig0_comp: np.ndarray, ocr: float = -1.0, pop: float = -1.0
) -> HSState:
    """Initialize a consistent committed HSState from an initial compressive stress sig0_comp."""
    st = HSState(
        sig=np.array(sig0_comp, dtype=float).copy(),
        eps=np.zeros(6, dtype=float),
        void_ratio=p.e0,
    )
    p_mean, q, _, Theta, _, _ = invariants_and_lode_grad(p, st.sig)
    q_tilde = Theta * q
    p_hat = max(p.p_min, p_mean + p.a_apex)

    # Equivalent preconsolidation pressure on the cap passing through (p_mean, q_tilde):
    p_eq = math.sqrt(p_hat * p_hat + (q_tilde / p.alpha_cap) ** 2) - p.a_apex
    ocr_use = p.OCR if ocr <= 0.0 else ocr
    pop_use = p.pop if pop < 0.0 else pop
    st.p_p = max(p.p_min, max(p_eq * ocr_use, p_eq + pop_use))

    # Mobilized shear hardening parameter gamma_p on the initial stress state:
    if q_tilde > 1.0e-12:
        st.gamma_p = max(0.0, _phi_shear(p, p_mean, q_tilde)[0])
    else:
        st.gamma_p = 0.0
    return st


# --------------------------------------------------------------------------- #
#  Constitutive Helper Functions in (p, q_tilde) Space                        #
# --------------------------------------------------------------------------- #
def _stiffness_scale_and_derivs(
    p: HSParams, p_mean: float, q_tilde: float
) -> Tuple[float, float, float]:
    """Return (s, ds/dp, ds/dq_tilde) for stress-dependent stiffness scaling."""
    p_hat = p_mean + p.a_apex
    if p.stiff_stress == "sigma3":
        sig3_hat = p_hat - (1.0 / 3.0) * q_tilde
        if sig3_hat <= p.p_min:
            s = (p.p_min / p.p_hat_ref) ** p.m
            return s, 0.0, 0.0
        ratio = sig3_hat / p.p_hat_ref
        s = ratio ** p.m
        ds_dsig3 = (p.m / p.p_hat_ref) * (ratio ** (p.m - 1.0)) if p.m > 0.0 else 0.0
        return s, ds_dsig3, -(1.0 / 3.0) * ds_dsig3
    else:
        if p_hat <= p.p_min:
            s = (p.p_min / p.p_hat_ref) ** p.m
            return s, 0.0, 0.0
        ratio = p_hat / p.p_hat_ref
        s = ratio ** p.m
        ds_dp = (p.m / p.p_hat_ref) * (ratio ** (p.m - 1.0)) if p.m > 0.0 else 0.0
        return s, ds_dp, 0.0


def _phi_shear(
    p: HSParams, p_mean: float, q_tilde: float
) -> Tuple[float, float, float]:
    """Return (Phi_s, dPhi_s/dp, dPhi_s/dq_tilde) for the hyperbolic shear cone f^s."""
    p_hat = max(p.p_min, p_mean + p.a_apex)
    sig3_hat = max(p.p_min, p_hat - (1.0 / 3.0) * q_tilde)
    dsig3_dp = 1.0 if (p_hat - (1.0 / 3.0) * q_tilde) > p.p_min else 0.0
    dsig3_dq = -(1.0 / 3.0) if (p_hat - (1.0 / 3.0) * q_tilde) > p.p_min else 0.0

    kf = 2.0 * p.sin_phi / (1.0 - p.sin_phi)
    qa = (kf / p.Rf) * sig3_hat
    dqa_dp = (kf / p.Rf) * dsig3_dp
    dqa_dq = (kf / p.Rf) * dsig3_dq

    r = min(0.9999, max(0.0, q_tilde / qa))
    if q_tilde / qa < 0.9999:
        dr_dp = -(q_tilde / (qa * qa)) * dqa_dp
        dr_dq = (1.0 / qa) - (q_tilde / (qa * qa)) * dqa_dq
    else:
        dr_dp = 0.0
        dr_dq = 0.0

    one_m_r = 1.0 - r
    c_hyp = 2.0 / (p.Ei_ref * one_m_r) - 2.0 / p.Eur_ref
    dc_dr = 2.0 / (p.Ei_ref * one_m_r * one_m_r)
    dc_dp = dc_dr * dr_dp
    dc_dq = dc_dr * dr_dq

    s, ds_dp, ds_dq = _stiffness_scale_and_derivs(p, p_mean, q_tilde)
    Phi_s = c_hyp * q_tilde / s
    dPhi_dp = (dc_dp * q_tilde) / s - (Phi_s / s) * ds_dp
    dPhi_dq = (dc_dq * q_tilde + c_hyp) / s - (Phi_s / s) * ds_dq
    return Phi_s, dPhi_dp, dPhi_dq


def _rowe_dilatancy(
    p: HSParams, p_mean: float, q_tilde: float, void_ratio: float
) -> Tuple[float, float, float]:
    """Return (sin_psi_m, d(sin_psi_m)/dp, d(sin_psi_m)/dq_tilde) per Rowe (1962)."""
    p_hat = max(p.p_min, p_mean + p.a_apex)
    denom = max(p.p_min, 2.0 * p_hat + (1.0 / 3.0) * q_tilde)
    s_phim = min(0.995, max(0.0, q_tilde / denom))
    if 0.0 < q_tilde / denom < 0.995 and (2.0 * p_hat + (1.0 / 3.0) * q_tilde) > p.p_min:
        dsphim_dp = -2.0 * q_tilde / (denom * denom)
        dsphim_dq = 2.0 * p_hat / (denom * denom)
    else:
        dsphim_dp = 0.0
        dsphim_dq = 0.0

    scv = p.sin_phicv
    denom_r = 1.0 - s_phim * scv
    s_psim = (s_phim - scv) / denom_r
    dspsim_dsphim = (1.0 - scv * scv) / (denom_r * denom_r)

    # Dilatancy cutoff when void_ratio >= emax and soil attempts to dilate (s_psim > 0):
    if p.emax > 0.0 and void_ratio >= p.emax and s_psim > 0.0:
        return 0.0, 0.0, 0.0

    return s_psim, dspsim_dsphim * dsphim_dp, dspsim_dsphim * dsphim_dq


def _cap_pp_and_deriv(p: HSParams, p_p_n: float, deps_v_pc: float) -> Tuple[float, float]:
    """Return (p_p_new, dp_p_new / d(deps_v_pc)) from exact integrated cap hardening law."""
    pp_hat_n = max(p.p_min, p_p_n + p.a_apex)
    if deps_v_pc <= 0.0:
        return p_p_n, 0.0

    if abs(1.0 - p.m) > 1.0e-8:
        om = 1.0 - p.m
        base_n = (pp_hat_n / p.p_hat_ref) ** om
        inc = om * deps_v_pc / (p.beta_cap * p.p_hat_ref)
        val = max(1.0e-12, base_n + inc)
        pp_hat_new = p.p_hat_ref * (val ** (1.0 / om))
        dpp_ddeps = (1.0 / p.beta_cap) * ((pp_hat_new / p.p_hat_ref) ** p.m)
        return pp_hat_new - p.a_apex, dpp_ddeps
    else:
        # m == 1 logarithmic Cam-Clay limit:
        arg = deps_v_pc / (p.beta_cap * p.p_hat_ref)
        pp_hat_new = pp_hat_n * math.exp(min(50.0, arg))
        dpp_ddeps = pp_hat_new / (p.beta_cap * p.p_hat_ref)
        return pp_hat_new - p.a_apex, dpp_ddeps


# --------------------------------------------------------------------------- #
#  Benz (2007) HSsmall Paraelastic Predictor & Exact Fourth-Order Tangent     #
# --------------------------------------------------------------------------- #
def _phi_hssmall_integral(
    G0: float, Gur: float, c_h: float, gam: float
) -> Tuple[float, float]:
    """Return (Phi(gam) = int_0^gam G_t(s) ds, G_t(gam))."""
    if G0 <= Gur or c_h <= 0.0:
        return Gur * gam, Gur
    gam_cut = (math.sqrt(G0 / Gur) - 1.0) / c_h
    if gam <= gam_cut:
        op = 1.0 + c_h * gam
        return G0 * gam / op, G0 / (op * op)
    else:
        op_c = 1.0 + c_h * gam_cut
        val = G0 * gam_cut / op_c + Gur * (gam - gam_cut)
        return val, Gur


def _hssmall_predictor(
    p: HSParams, st_n: HSState, eps_new: np.ndarray, s_n: float
) -> Tuple[float, np.ndarray, np.ndarray, bool, np.ndarray, np.ndarray]:
    """Evaluate Benz (2007) HSsmall shear modulus G_e, updated history (e_dev, H_hist, is_virgin),
    6x6 paraelastic tangent C_e_alg = d(sig_tr)/d(eps_new), and 6-vector dGe/d(eps_new).
    """
    Gur = p.Gur_ref * s_n
    Kur = p.Kur_ref * s_n
    e_new = strain_voigt_to_tensor(eps_new)
    eps_v_new = np.trace(e_new)
    e_dev_new = e_new - (eps_v_new / 3.0) * _I3
    de_dev = e_dev_new - st_n.e_dev

    if (not p.use_small_strain) or p.G0_ref <= p.Gur_ref:
        C_e = Kur * np.outer(_I6_M, _I6_M) + 2.0 * Gur * _IDEV6
        return Gur, e_dev_new, st_n.H_hist.copy(), st_n.is_virgin, C_e, np.zeros(6)

    G0 = p.G0_ref * s_n
    vals, vecs = np.linalg.eigh(de_dev)
    H_hat = vecs.T @ st_n.H_hist @ vecs
    w = np.ones(3, dtype=float)
    reversed_any = False
    for a in range(3):
        if vals[a] * H_hat[a, a] < -1.0e-16:
            w[a] = 0.0
            if abs(vals[a]) > 1.0e-12 and abs(H_hat[a, a]) > 1.0e-12:
                reversed_any = True

    virgin_new = st_n.is_virgin and (not reversed_any)
    W_mat = np.sqrt(np.outer(w, w))
    H0_hat = W_mat * H_hat
    H1_hat = H0_hat + np.diag(vals)
    H0 = vecs @ H0_hat @ vecs.T
    H1 = vecs @ H1_hat @ vecs.T

    gam0 = math.sqrt(max(0.0, 2.0 * float(np.sum(H0 * H0))))
    gam1 = math.sqrt(max(0.0, 2.0 * float(np.sum(H1 * H1))))
    m_fac = 1.0 if virgin_new else 2.0
    c_h = 0.385 / (m_fac * max(1.0e-12, p.gamma_07))

    phi1, gt1 = _phi_hssmall_integral(G0, Gur, c_h, gam1)
    phi0, _ = _phi_hssmall_integral(G0, Gur, c_h, gam0)

    dgam = gam1 - gam0
    if abs(dgam) > 1.0e-12 and gam1 > 1.0e-14:
        G_sec = (phi1 - phi0) / dgam
        dGsec_dgam1 = (gt1 - G_sec) / dgam
        H1_voigt = stress_tensor_to_voigt(H1 - (np.trace(H1) / 3.0) * _I3)
        de_voigt = stress_tensor_to_voigt(de_dev)
        dgam1_deps = (2.0 / gam1) * H1_voigt
        dGe_deps = dGsec_dgam1 * dgam1_deps
        C_e = (
            Kur * np.outer(_I6_M, _I6_M)
            + 2.0 * G_sec * _IDEV6
            + 2.0 * np.outer(de_voigt, dGe_deps)
        )
        return G_sec, e_dev_new, H1, virgin_new, C_e, dGe_deps
    else:
        C_e = Kur * np.outer(_I6_M, _I6_M) + 2.0 * gt1 * _IDEV6
        return gt1, e_dev_new, H1, virgin_new, C_e, np.zeros(6)


# --------------------------------------------------------------------------- #
#  2D (p, q) Return-Map Solvers + Analytical Sensitivity d(p,q)/d(p_tr,q_tr,Th,Ge)#
# --------------------------------------------------------------------------- #
def _solve_regime1_shear(
    p: HSParams,
    p_tr: float,
    q_tr: float,
    Theta: float,
    Ke: float,
    Ge: float,
    gamma_p_n: float,
    void_ratio: float,
) -> Tuple[bool, float, float, float, np.ndarray]:
    """Solve Regime 1 (Hyperbolic shear cone f^s = 0 only, dlambda_c = 0).

    Unknowns y = [p, q].
    R0(p, q) = p - p_tr - Ke * sin_psi_m(p, Th*q) * (2*(q_tr - q)/(3*Ge)) = 0
    R1(p, q) = Phi_s(p, Th*q) - gamma_p_n - 2*(q_tr - q)/(3*Ge) = 0
    """
    inv_3Ge2 = 2.0 / (3.0 * Ge)
    p_cur = max(p.p_min - p.a_apex, p_tr)
    q_max_fail = p.M_c * max(p.p_min, p_cur + p.a_apex) / Theta
    q_cur = min(q_tr, 0.95 * q_max_fail)

    for _ in range(35):
        q_tilde = Theta * q_cur
        p_hat = p_cur + p.a_apex
        if p_hat <= p.p_min or q_tilde >= p.M_c * p_hat:
            return False, p_cur, q_cur, 0.0, np.zeros((2, 4))

        dgam = inv_3Ge2 * (q_tr - q_cur)
        s_psim, dspsi_dp, dspsi_dqtilde = _rowe_dilatancy(p, p_cur, q_tilde, void_ratio)
        Phi_s, dPhi_dp, dPhi_dqtilde = _phi_shear(p, p_cur, q_tilde)

        R0 = p_cur - p_tr - Ke * s_psim * dgam
        R1 = Phi_s - gamma_p_n - dgam
        if max(abs(R0) / max(1.0, p.p_ref), abs(R1)) < 1.0e-13:
            if dgam < -1.0e-12:
                return False, p_cur, q_cur, 0.0, np.zeros((2, 4))
            # Analytical Jacobian J_y and J_tr = dR / d(p_tr, q_tr, Theta, Ge):
            J00 = 1.0 - Ke * dspsi_dp * dgam
            J01 = -Ke * (dspsi_dqtilde * Theta * dgam - s_psim * inv_3Ge2)
            J10 = dPhi_dp
            J11 = dPhi_dqtilde * Theta + inv_3Ge2
            J_y = np.array([[J00, J01], [J10, J11]], dtype=float)
            J_tr = np.array(
                [
                    [-1.0, -Ke * s_psim * inv_3Ge2, -Ke * dspsi_dqtilde * q_cur * dgam, Ke * s_psim * dgam / Ge],
                    [0.0, -inv_3Ge2, dPhi_dqtilde * q_cur, dgam / Ge],
                ],
                dtype=float,
            )
            sens = -np.linalg.solve(J_y, J_tr)
            return True, p_cur, q_cur, max(0.0, dgam), sens

        J00 = 1.0 - Ke * dspsi_dp * dgam
        J01 = -Ke * (dspsi_dqtilde * Theta * dgam - s_psim * inv_3Ge2)
        J10 = dPhi_dp
        J11 = dPhi_dqtilde * Theta + inv_3Ge2
        det = J00 * J11 - J01 * J10
        if abs(det) < 1.0e-30:
            return False, p_cur, q_cur, 0.0, np.zeros((2, 4))
        dp = -(J11 * R0 - J01 * R1) / det
        dq = -(-J10 * R0 + J00 * R1) / det

        alpha_ls = 1.0
        for _ in range(12):
            p_next = p_cur + alpha_ls * dp
            q_next = max(0.0, min(q_tr, q_cur + alpha_ls * dq))
            if (p_next + p.a_apex) > p.p_min and (Theta * q_next) < 0.9995 * p.M_c * (p_next + p.a_apex):
                p_cur, q_cur = p_next, q_next
                break
            alpha_ls *= 0.5
        else:
            return False, p_cur, q_cur, 0.0, np.zeros((2, 4))

    return False, p_cur, q_cur, 0.0, np.zeros((2, 4))


def _solve_regime2_failure(
    p: HSParams,
    p_tr: float,
    q_tr: float,
    Theta: float,
    Ke: float,
    Ge: float,
    void_ratio: float,
) -> Tuple[bool, float, float, float, np.ndarray]:
    """Closed-form solution for Regime 2 (Ultimate failure truncation f^f = 0 only)."""
    s_psi = 0.0 if (p.emax > 0.0 and void_ratio >= p.emax and p.sin_psi > 0.0) else p.sin_psi
    denom = 1.5 * Ge * Theta + p.M_c * Ke * s_psi
    dgam = (Theta * q_tr - p.M_c * (p_tr + p.a_apex)) / denom
    if dgam < -1.0e-12:
        return False, p_tr, q_tr, 0.0, np.zeros((2, 4))
    dgam = max(0.0, dgam)
    p_new = p_tr + Ke * s_psi * dgam
    q_new = max(0.0, q_tr - 1.5 * Ge * dgam)
    if (p_new + p.a_apex) <= p.p_min:
        return False, p_new, q_new, dgam, np.zeros((2, 4))

    # Exact derivatives of (p_new, q_new) wrt (p_tr, q_tr, Theta, Ge):
    ddgam_dptr = -p.M_c / denom
    ddgam_dqtr = Theta / denom
    ddgam_dTh = (q_tr - 1.5 * Ge * dgam) / denom
    ddgam_dGe = -1.5 * Theta * dgam / denom
    sens = np.array(
        [
            [
                1.0 + Ke * s_psi * ddgam_dptr,
                Ke * s_psi * ddgam_dqtr,
                Ke * s_psi * ddgam_dTh,
                Ke * s_psi * ddgam_dGe,
            ],
            [
                -1.5 * Ge * ddgam_dptr,
                1.0 - 1.5 * Ge * ddgam_dqtr,
                -1.5 * Ge * ddgam_dTh,
                -1.5 * dgam - 1.5 * Ge * ddgam_dGe,
            ],
        ],
        dtype=float,
    )
    return True, p_new, q_new, dgam, sens


def _solve_regime3_cap(
    p: HSParams,
    p_tr: float,
    q_tr: float,
    Theta: float,
    Ke: float,
    Ge: float,
    p_p_n: float,
) -> Tuple[bool, float, float, float, float, np.ndarray]:
    """Solve Regime 3 (Elliptical volumetric cap f^c = 0 only, dgam_p = 0)."""
    alpha2 = p.alpha_cap * p.alpha_cap
    p_cur = max(p_p_n * 0.5, min(p_tr, p_p_n))
    q_cur = q_tr

    for _ in range(35):
        p_hat = max(p.p_min, p_cur + p.a_apex)
        deps_v_pc = max(0.0, (p_tr - p_cur) / Ke)
        c_q = 1.0 + (3.0 * Ge * Theta * Theta / alpha2) * (deps_v_pc / p_hat)
        q_cur = q_tr / c_q
        pp_new, dpp_ddeps = _cap_pp_and_deriv(p, p_p_n, deps_v_pc)
        pp_hat = pp_new + p.a_apex

        q_tilde = Theta * q_cur
        R = (q_tilde * q_tilde) / alpha2 + p_hat * p_hat - pp_hat * pp_hat
        if abs(R) / (p.p_hat_ref * p.p_hat_ref) < 1.0e-13:
            dcq_dp = (3.0 * Ge * Theta * Theta / alpha2) * (
                -(1.0 / (Ke * p_hat)) - (deps_v_pc / (p_hat * p_hat))
            )
            dcq_dptr = (3.0 * Ge * Theta * Theta / alpha2) * (1.0 / (Ke * p_hat))
            dcq_dTh = (6.0 * Ge * Theta / alpha2) * (deps_v_pc / p_hat)
            dcq_dGe = (3.0 * Theta * Theta / alpha2) * (deps_v_pc / p_hat)

            J00 = q_cur * dcq_dp
            J01 = c_q
            J10 = 2.0 * p_hat + 2.0 * pp_hat * dpp_ddeps * (1.0 / Ke)
            J11 = 2.0 * Theta * Theta * q_cur / alpha2
            J_y = np.array([[J00, J01], [J10, J11]], dtype=float)
            J_tr = np.array(
                [
                    [q_cur * dcq_dptr, -1.0, q_cur * dcq_dTh, q_cur * dcq_dGe],
                    [-2.0 * pp_hat * dpp_ddeps * (1.0 / Ke), 0.0, 2.0 * Theta * q_cur * q_cur / alpha2, 0.0],
                ],
                dtype=float,
            )
            sens = -np.linalg.solve(J_y, J_tr)
            return True, p_cur, q_cur, pp_new, deps_v_pc, sens

        dcq_dp = (3.0 * Ge * Theta * Theta / alpha2) * (
            -(1.0 / (Ke * p_hat)) - (deps_v_pc / (p_hat * p_hat))
        )
        dq_dp = -(q_cur / c_q) * dcq_dp
        dR_dp = (
            (2.0 * Theta * Theta * q_cur / alpha2) * dq_dp
            + 2.0 * p_hat
            + 2.0 * pp_hat * dpp_ddeps * (1.0 / Ke)
        )
        dp = -R / dR_dp
        p_cur = max(p.p_min - p.a_apex, min(p_tr, p_cur + dp))

    return False, p_cur, q_cur, p_p_n, 0.0, np.zeros((2, 4))


def _solve_regime4_corner(
    p: HSParams,
    p_tr: float,
    q_tr: float,
    Theta: float,
    Ke: float,
    Ge: float,
    gamma_p_n: float,
    p_p_n: float,
    void_ratio: float,
    failure_mode: bool = False,
) -> Tuple[bool, float, float, float, float, np.ndarray]:
    """Solve Regime 4 (Simultaneous Shear [f^s or f^f] + Cap [f^c] corner)."""
    alpha2 = p.alpha_cap * p.alpha_cap
    p_cur = max(p.p_min - p.a_apex, min(p_tr, max(p_p_n, 0.7 * p_tr)))
    q_cur = min(q_tr, 0.85 * p.M_c * (p_cur + p.a_apex) / Theta)

    def eval_residuals(
        p_val: float, q_val: float, p_tr_v: float, q_tr_v: float, Th_v: float, Ge_v: float
    ):
        p_hat = max(p.p_min, p_val + p.a_apex)
        q_tilde = Th_v * q_val
        deps_v_p = (p_tr_v - p_val) / Ke
        deps_q_p = (q_tr_v - q_val) / (3.0 * Ge_v)
        if failure_mode:
            s_psim = 0.0 if (p.emax > 0.0 and void_ratio >= p.emax and p.sin_psi > 0.0) else p.sin_psi
        else:
            s_psim, _, _ = _rowe_dilatancy(p, p_val, q_tilde, void_ratio)

        k_cap = (Th_v * Th_v * q_val) / (alpha2 * p_hat)
        denom = max(0.05, 0.5 + k_cap * s_psim)
        dgam_p = (deps_q_p - k_cap * deps_v_p) / denom
        deps_v_pc = deps_v_p + s_psim * dgam_p

        pp_new, _ = _cap_pp_and_deriv(p, p_p_n, max(0.0, deps_v_pc))
        pp_hat = pp_new + p.a_apex

        if failure_mode:
            R0 = q_tilde - p.M_c * p_hat
        else:
            Phi_s, _, _ = _phi_shear(p, p_val, q_tilde)
            R0 = Phi_s - gamma_p_n - dgam_p
        R1 = (q_tilde * q_tilde) / alpha2 + p_hat * p_hat - pp_hat * pp_hat
        return R0, R1, dgam_p, deps_v_pc, pp_new

    for _ in range(40):
        R0, R1, dgam_p, deps_v_pc, pp_new = eval_residuals(p_cur, q_cur, p_tr, q_tr, Theta, Ge)
        if max(abs(R0) / max(1.0, p.p_ref if failure_mode else 1.0), abs(R1) / (p.p_hat_ref ** 2)) < 1.0e-12:
            if dgam_p < -1.0e-10 or deps_v_pc < -1.0e-10:
                return False, p_cur, q_cur, 0.0, p_p_n, np.zeros((2, 4))
            J_y = np.zeros((2, 2), dtype=float)
            J_tr = np.zeros((2, 4), dtype=float)
            hp = max(1.0e-8, 1.0e-8 * abs(p_cur))
            hq = max(1.0e-8, 1.0e-8 * max(1.0, abs(q_cur)))
            hth = 1.0e-8
            hge = max(1.0e-6, 1.0e-8 * Ge)
            r0p, r1p, _, _, _ = eval_residuals(p_cur + hp, q_cur, p_tr, q_tr, Theta, Ge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur - hp, q_cur, p_tr, q_tr, Theta, Ge)
            J_y[0, 0] = (r0p - r0m) / (2.0 * hp)
            J_y[1, 0] = (r1p - r1m) / (2.0 * hp)

            r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur + hq, p_tr, q_tr, Theta, Ge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur - hq, p_tr, q_tr, Theta, Ge)
            J_y[0, 1] = (r0p - r0m) / (2.0 * hq)
            J_y[1, 1] = (r1p - r1m) / (2.0 * hq)

            r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur, p_tr + hp, q_tr, Theta, Ge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur, p_tr - hp, q_tr, Theta, Ge)
            J_tr[0, 0] = (r0p - r0m) / (2.0 * hp)
            J_tr[1, 0] = (r1p - r1m) / (2.0 * hp)

            r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr + hq, Theta, Ge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr - hq, Theta, Ge)
            J_tr[0, 1] = (r0p - r0m) / (2.0 * hq)
            J_tr[1, 1] = (r1p - r1m) / (2.0 * hq)

            r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr, Theta + hth, Ge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr, Theta - hth, Ge)
            J_tr[0, 2] = (r0p - r0m) / (2.0 * hth)
            J_tr[1, 2] = (r1p - r1m) / (2.0 * hth)

            r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr, Theta, Ge + hge)
            r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur, p_tr, q_tr, Theta, Ge - hge)
            J_tr[0, 3] = (r0p - r0m) / (2.0 * hge)
            J_tr[1, 3] = (r1p - r1m) / (2.0 * hge)

            sens = -np.linalg.solve(J_y, J_tr)
            return True, p_cur, q_cur, max(0.0, dgam_p), pp_new, sens

        hp = max(1.0e-7, 1.0e-7 * abs(p_cur))
        hq = max(1.0e-7, 1.0e-7 * max(1.0, abs(q_cur)))
        r0p, r1p, _, _, _ = eval_residuals(p_cur + hp, q_cur, p_tr, q_tr, Theta, Ge)
        r0m, r1m, _, _, _ = eval_residuals(p_cur - hp, q_cur, p_tr, q_tr, Theta, Ge)
        J00 = (r0p - r0m) / (2.0 * hp)
        J10 = (r1p - r1m) / (2.0 * hp)
        r0p, r1p, _, _, _ = eval_residuals(p_cur, q_cur + hq, p_tr, q_tr, Theta, Ge)
        r0m, r1m, _, _, _ = eval_residuals(p_cur, q_cur - hq, p_tr, q_tr, Theta, Ge)
        J01 = (r0p - r0m) / (2.0 * hq)
        J11 = (r1p - r1m) / (2.0 * hq)
        det = J00 * J11 - J01 * J10
        if abs(det) < 1.0e-30:
            break
        dp = -(J11 * R0 - J01 * R1) / det
        dq = -(-J10 * R0 + J00 * R1) / det

        alpha_ls = 1.0
        for _ in range(10):
            p_next = p_cur + alpha_ls * dp
            q_next = max(0.0, q_cur + alpha_ls * dq)
            if (p_next + p.a_apex) > p.p_min and (
                failure_mode or (Theta * q_next) < 0.9995 * p.M_c * (p_next + p.a_apex)
            ):
                p_cur, q_cur = p_next, q_next
                break
            alpha_ls *= 0.5
        else:
            break

    return False, p_cur, q_cur, 0.0, p_p_n, np.zeros((2, 4))


# --------------------------------------------------------------------------- #
#  Full 3D Return Map & 6x6 Consistent Algorithmic Tangent                    #
# --------------------------------------------------------------------------- #
def return_map_3d(
    p: HSParams, st_n: HSState, eps_new: np.ndarray
) -> Tuple[HSState, np.ndarray, str]:
    """Perform the 4-surface Hardening Soil + HSsmall return map for strain eps_new.

    Returns:
      (st_new, C_alg_6x6, regime_label)
    where regime_label is one of:
      'elastic', 'shear', 'failure', 'cap', 'shear_cap', 'failure_cap', 'tension'
    """
    st_new = st_n.copy()
    st_new.eps = np.array(eps_new, dtype=float).copy()
    deps = st_new.eps - st_n.eps
    deps_v = deps[0] + deps[1] + deps[2]

    # Void ratio update (compression-positive: deps_v > 0 decreases void ratio):
    st_new.void_ratio = st_n.void_ratio - (1.0 + st_n.void_ratio) * deps_v

    # Committed reference stress scale s_n:
    if np.linalg.norm(st_n.sig) < 1.0e-12:
        s_n = 1.0
    else:
        p_n, q_n, _, Th_n, _, _ = invariants_and_lode_grad(p, st_n.sig)
        s_n, _, _ = _stiffness_scale_and_derivs(p, p_n, Th_n * q_n)

    Ke = p.Kur_ref * s_n
    Ge, e_dev_new, H_new, virgin_new, C_e_alg, dGe_deps = _hssmall_predictor(
        p, st_n, st_new.eps, s_n
    )
    st_new.e_dev = e_dev_new
    st_new.H_hist = H_new
    st_new.is_virgin = virgin_new

    deps_dev_sig = np.array(
        [
            deps[0] - deps_v / 3.0,
            deps[1] - deps_v / 3.0,
            deps[2] - deps_v / 3.0,
            0.5 * deps[3],
            0.5 * deps[4],
            0.5 * deps[5],
        ],
        dtype=float,
    )
    sig_tr = st_n.sig + Ke * deps_v * _I6_M + 2.0 * Ge * deps_dev_sig
    p_tr, q_tr, _, Theta, dq_dsig_tr, dTheta_dsig_tr = invariants_and_lode_grad(p, sig_tr)

    # 1. Check Rankine / Apex Tension Cutoff (Surface 4):
    p_cut = max(p.p_min - p.a_apex, -p.sigma_t)
    if p_tr <= p_cut:
        st_new.sig = p_cut * _I6_M + 1.0e-3 * (sig_tr - p_cut * _I6_M)
        C_floor = 1.0e-3 * C_e_alg
        return st_new, C_floor, "tension"

    # 2. Evaluate trial yield functions at (p_tr, q_tr, Theta):
    q_tilde_tr = Theta * q_tr
    p_hat_tr = p_tr + p.a_apex
    f_fail_tr = q_tilde_tr - p.M_c * p_hat_tr
    if f_fail_tr < 0.0:
        Phi_s_tr, _, _ = _phi_shear(p, p_tr, q_tilde_tr)
        f_shear_tr = Phi_s_tr - st_n.gamma_p
    else:
        f_shear_tr = 1.0

    pp_hat_n = st_n.p_p + p.a_apex
    f_cap_tr = (q_tilde_tr * q_tilde_tr) / (p.alpha_cap ** 2) + p_hat_tr * p_hat_tr - pp_hat_n * pp_hat_n

    if f_shear_tr <= 1.0e-12 and f_fail_tr <= 1.0e-12 and f_cap_tr <= 1.0e-10:
        st_new.sig = sig_tr
        return st_new, C_e_alg, "elastic"

    # 3. Koiter Active-Set Selection & Return Map in (p, q):
    regime = "elastic"
    p_new, q_new = p_tr, q_tr
    sens = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]], dtype=float)
    solved = False

    # Try single-surface Shear (Regime 1) or Failure (Regime 2) first if shear is violated:
    # Note: e_max dilatancy cutoff is gated on st_n.void_ratio for C^1 Newton-Raphson smoothness.
    if f_shear_tr > 1.0e-12 or f_fail_tr > 1.0e-12:
        ok1, p1, q1, dgam1, sens1 = _solve_regime1_shear(
            p, p_tr, q_tr, Theta, Ke, Ge, st_n.gamma_p, st_n.void_ratio
        )
        if ok1:
            fc1 = ((Theta * q1) ** 2) / (p.alpha_cap ** 2) + (p1 + p.a_apex) ** 2 - pp_hat_n ** 2
            if fc1 <= 1.0e-8:
                p_new, q_new, st_new.gamma_p, sens = p1, q1, st_n.gamma_p + dgam1, sens1
                regime = "shear"
                solved = True
        if not solved:
            ok2, p2, q2, dgam2, sens2 = _solve_regime2_failure(
                p, p_tr, q_tr, Theta, Ke, Ge, st_n.void_ratio
            )
            if ok2:
                fc2 = ((Theta * q2) ** 2) / (p.alpha_cap ** 2) + (p2 + p.a_apex) ** 2 - pp_hat_n ** 2
                if fc2 <= 1.0e-8:
                    p_new, q_new, st_new.gamma_p, sens = p2, q2, st_n.gamma_p + dgam2, sens2
                    regime = "failure"
                    solved = True

    # Try single-surface Cap (Regime 3) if not yet solved and cap is violated:
    if (not solved) and f_cap_tr > 1.0e-10:
        ok3, p3, q3, pp3, _, sens3 = _solve_regime3_cap(
            p, p_tr, q_tr, Theta, Ke, Ge, st_n.p_p
        )
        if ok3:
            q_tilde3 = Theta * q3
            ff3 = q_tilde3 - p.M_c * (p3 + p.a_apex)
            fs3 = (_phi_shear(p, p3, q_tilde3)[0] - st_n.gamma_p) if ff3 < 0.0 else 1.0
            if fs3 <= 1.0e-13 and ff3 <= 1.0e-13:
                p_new, q_new, st_new.p_p, sens = p3, q3, pp3, sens3
                regime = "cap"
                solved = True

    # Otherwise, solve Regime 4 (Simultaneous Shear + Cap Corner):
    if not solved:
        ok4, p4, q4, dgam4, pp4, sens4 = _solve_regime4_corner(
            p, p_tr, q_tr, Theta, Ke, Ge, st_n.gamma_p, st_n.p_p, st_n.void_ratio, failure_mode=False
        )
        if ok4:
            p_new, q_new, st_new.gamma_p, st_new.p_p, sens = (
                p4,
                q4,
                st_n.gamma_p + dgam4,
                pp4,
                sens4,
            )
            regime = "shear_cap"
            solved = True
        else:
            ok4f, p4f, q4f, dgam4f, pp4f, sens4f = _solve_regime4_corner(
                p, p_tr, q_tr, Theta, Ke, Ge, st_n.gamma_p, st_n.p_p, st_n.void_ratio, failure_mode=True
            )
            if ok4f:
                p_new, q_new, st_new.gamma_p, st_new.p_p, sens = (
                    p4f,
                    q4f,
                    st_n.gamma_p + dgam4f,
                    pp4f,
                    sens4f,
                )
                regime = "failure_cap"
                solved = True

    # 4. Reconstruct 3D stress tensor & 6x6 Consistent Algorithmic Tangent:
    s_tr = sig_tr - p_tr * _I6_M
    if q_tr > 1.0e-12:
        scale_q = q_new / q_tr
        st_new.sig = p_new * _I6_M + scale_q * s_tr

        dp_tr_dsig = _I6_M / 3.0
        grad_p = (
            sens[0, 0] * dp_tr_dsig
            + sens[0, 1] * dq_dsig_tr
            + sens[0, 2] * dTheta_dsig_tr
        )
        grad_q = (
            sens[1, 0] * dp_tr_dsig
            + sens[1, 1] * dq_dsig_tr
            + sens[1, 2] * dTheta_dsig_tr
        )
        I_dev_stress = np.eye(6, dtype=float) - np.outer(_I6_M, _I6_M) / 3.0
        n_s = s_tr / q_tr
        dsig_dsig_tr = (
            np.outer(_I6_M, grad_p)
            + np.outer(n_s, grad_q)
            + scale_q * (I_dev_stress - np.outer(n_s, dq_dsig_tr))
        )
        C_alg = (
            dsig_dsig_tr @ C_e_alg
            + sens[0, 3] * np.outer(_I6_M, dGe_deps)
            + sens[1, 3] * np.outer(n_s, dGe_deps)
        )
    else:
        st_new.sig = p_new * _I6_M
        dp_tr_dsig = _I6_M / 3.0
        grad_p = sens[0, 0] * dp_tr_dsig
        dsig_dsig_tr = np.outer(_I6_M, grad_p)
        C_alg = dsig_dsig_tr @ C_e_alg + sens[0, 3] * np.outer(_I6_M, dGe_deps)

    return st_new, C_alg, regime
