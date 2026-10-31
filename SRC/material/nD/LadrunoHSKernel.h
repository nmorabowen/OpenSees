/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno: LadrunoHSKernel.h -- OpenSees-free C++17 constitutive kernel for the
// Hardening Soil (HS, Schanz, Vermeer & Bonnier 1999) and Hardening Soil with
// Small-Strain Stiffness (HSsmall, Benz 2007) models with selectable Lode-angle
// cross-section (-lodeShape mc vs -lodeShape mn, Panteghini & Lagioia 2014).
//
// Sign convention inside the kernel (`return_map_3d_comp`):
//   Compression-positive soil mechanics:
//     p = (sig11 + sig22 + sig33) / 3 > 0 in compression,
//     eps_v = eps11 + eps22 + eps33 > 0 in compression.
//   Voigt ordering (order 6):
//     sig = [sig11, sig22, sig33, tau12, tau23, tau31]
//     eps = [eps11, eps22, eps33, gam12, gam23, gam31]  (engineering shear)
//   OpenSees continuum wrapper (`return_map_3d_opensees`) flips stress & strain
//   signs on entry/exit (`sig_ops = -sig_comp`, `eps_ops = -eps_comp`), while
//   the 6x6 tangent `d(sig_ops)/d(eps_ops) == d(sig_comp)/d(eps_comp)` is
//   invariant under simultaneous sign flip.
//
// See Ladruno_implementation/210_ladruno_hs_plan.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoHSKernel_h
#define LadrunoHSKernel_h

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace LadrunoHSKernel {

static constexpr double kPi = 3.141592653589793238462643383279502884;

enum LodeShape : int {
  LODE_MC = 0,  // Mohr-Coulomb (Schanz et al. 1999 / PLAXIS HS / FLAC PH)
  LODE_MN = 1   // Matsuoka-Nakai (Panteghini & Lagioia 2014 / numgeo)
};

enum StiffStress : int {
  STIFF_SIGMA3 = 0,  // minor principal stress sigma'_3 (PLAXIS HS default)
  STIFF_MEAN_P = 1   // mean effective stress p' (GHS / isotropic variant)
};

enum ActiveRegime : int {
  REGIME_REFUSED     = -1,
  REGIME_ELASTIC     = 0,
  REGIME_SHEAR       = 1,
  REGIME_FAILURE     = 2,
  REGIME_CAP         = 3,
  REGIME_SHEAR_CAP   = 4,
  REGIME_FAILURE_CAP = 5,
  REGIME_TENSION     = 6
};

struct HSParams {
  // Primary user parameters:
  double E50_ref  = 30000.0;
  double Eoed_ref = 30000.0;
  double Eur_ref  = 90000.0;
  double m        = 0.50;
  double c        = 0.0;
  double phi_deg  = 35.0;
  double psi_deg  = 5.0;
  double nu_ur    = 0.20;
  double p_ref    = 100.0;
  double Rf       = 0.90;
  double K0nc     = -1.0;     // <= 0 -> auto Jaky (1 - sin(phi))
  double sigma_t  = 0.0;
  double OCR      = 1.0;
  double pop      = 0.0;
  double e0       = 0.50;
  double emax     = -1.0;     // <= 0 -> dilatancy cutoff disabled

  // Optional Benz (2007) HSsmall parameters:
  bool   use_small_strain = false;
  double G0_ref           = 0.0;
  double gamma_07         = 1.0e-4;

  // Formulation switches:
  int    lode_shape   = LODE_MC;
  int    stiff_stress = STIFF_SIGMA3;

  // Cap parameters (<= 0 -> closed-form auto-calibration from Eoed_ref & K0nc):
  double alpha_cap      = -1.0;
  double beta_cap       = -1.0;
  double user_alpha_cap = -1.0;
  double user_beta_cap  = -1.0;
  double p_min_ratio    = 1.0e-3;

  // Derived constants populated by finalize():
  double phi        = 0.0;
  double psi        = 0.0;
  double sin_phi    = 0.0;
  double cos_phi    = 1.0;
  double sin_psi    = 0.0;
  double sin_phicv  = 0.0;
  double a_apex     = 0.0;
  double p_hat_ref  = 100.0;
  double p_min      = 0.1;
  double M_c        = 1.0;
  double delta_lode = 1.0;
  double gamma_mn   = 0.0;
  double c0_mn      = 1.0;
  double gamma_mc   = 1.0 - 1.0e-7;
  double th_max_mc  = 0.0;
  double Ei_ref     = 60000.0;
  double Gur_ref    = 37500.0;
  double Kur_ref    = 50000.0;

  bool finalize() {
    if (E50_ref <= 0.0 || Eoed_ref <= 0.0 || Eur_ref <= 0.0) return false;
    if (Eur_ref < E50_ref) return false;
    if (m < 0.0 || m > 1.0) return false;
    if (phi_deg <= 0.0 || phi_deg >= 89.0) return false;
    if (nu_ur <= -0.99 || nu_ur >= 0.499) return false;
    if (p_ref <= 0.0 || Rf <= 0.0 || Rf >= 1.0) return false;
    if (use_small_strain && (G0_ref <= 0.0 || gamma_07 <= 0.0)) return false;

    phi     = phi_deg * (kPi / 180.0);
    psi     = psi_deg * (kPi / 180.0);
    sin_phi = std::sin(phi);
    cos_phi = std::cos(phi);
    sin_psi = std::sin(psi);

    // Rowe (1962) critical-state friction ratio sin(phi_cv):
    double num_cv = sin_phi - sin_psi;
    double den_cv = 1.0 - sin_phi * sin_psi;
    sin_phicv = std::max(0.0, num_cv / std::max(1.0e-12, den_cv));

    // Apex shift a = c * cot(phi'):
    a_apex    = (sin_phi > 1.0e-12) ? (c * cos_phi / sin_phi) : 0.0;
    p_hat_ref = p_ref + a_apex;
    p_min     = std::max(1.0e-6, p_min_ratio * p_ref);

    if (K0nc <= 0.0) {
      K0nc = std::max(0.15, 1.0 - sin_phi);
    }

    // Triaxial compression slope M_c and extension ratio delta:
    M_c        = 6.0 * sin_phi / (3.0 - sin_phi);
    delta_lode = (3.0 + sin_phi) / (3.0 - sin_phi);

    // Panteghini & Lagioia (2014) Eq. (18) Matsuoka-Nakai constants:
    // m_MN = (9 - sin^2(phi)) / (1 - sin^2(phi)) => m_MN - 3 = 2*(3 + sin^2(phi)) / (1 - sin^2(phi)).
    double s2   = sin_phi * sin_phi;
    double d_mn = 3.0 + s2;
    gamma_mn    = (9.0 - s2) * sin_phi / (d_mn * std::sqrt(d_mn));
    gamma_mn    = std::max(0.0, std::min(1.0 - 1.0e-14, gamma_mn));
    c0_mn       = std::cos(std::acos(-gamma_mn) / 3.0);

    // Endpoint-exact C^inf Mohr-Coulomb smoothing constants:
    gamma_mc  = 1.0 - 1.0e-7;
    th_max_mc = std::asin(gamma_mc) / 3.0;

    // Elastic and hyperbolic reference moduli (Schanz, Vermeer & Bonnier 1999):
    Ei_ref  = 2.0 * E50_ref / (2.0 - Rf);
    Gur_ref = Eur_ref / (2.0 * (1.0 + nu_ur));
    Kur_ref = Eur_ref / (3.0 * (1.0 - 2.0 * nu_ur));

    if (user_alpha_cap > 0.0) alpha_cap = user_alpha_cap;
    if (user_beta_cap > 0.0)  beta_cap  = user_beta_cap;
    if (alpha_cap <= 0.0 || beta_cap <= 0.0) {
      calibrate_cap_closed_form();
    }
    return true;
  }

  void calibrate_cap_closed_form() {
    // 1D virgin oedometer path: sig_1 + a = p_hat_ref, sig_3 + a = K0nc * p_hat_ref
    double eta_k0 = 3.0 * (1.0 - K0nc) / (1.0 + 2.0 * K0nc);
    double s_k0   = (stiff_stress == STIFF_SIGMA3)
                        ? std::pow(K0nc, m)
                        : std::pow((1.0 + 2.0 * K0nc) / 3.0, m);

    double kf = 2.0 * sin_phi / (1.0 - sin_phi);
    double r_k0 = std::min(0.99, Rf * (1.0 - K0nc) / (K0nc * kf));

    // Rowe dilatancy on the K0nc line:
    double s_phim_k0 = (1.0 - K0nc) / (1.0 + K0nc);
    double s_psim_k0 = (s_phim_k0 - sin_phicv) / (1.0 - s_phim_k0 * sin_phicv);

    // Elastic strain rates per unit d(sig_1):
    double deps_v_e_d1 = (1.0 + 2.0 * K0nc) / (3.0 * Kur_ref * s_k0);
    double deps_q_e_d1 = (1.0 - K0nc) / (3.0 * Gur_ref * s_k0);

    // Shear-cone plastic multiplier rate d(gamma^p) / d(sig_1) along proportional K0nc path:
    double c_hyp = 2.0 / (Ei_ref * (1.0 - r_k0)) - 2.0 / Eur_ref;
    double dgam_p_d1 = std::max(0.0, (1.0 - m) * c_hyp * (1.0 - K0nc) / s_k0);

    double deps_q_ps_d1 = 0.5 * dgam_p_d1;
    double deps_v_ps_d1 = -s_psim_k0 * dgam_p_d1;

    // Total 1D oedometric strain rate per unit d(sig_1) is 1 / Eoed_ref:
    double deps1_d1     = 1.0 / Eoed_ref;
    double deps_v_pc_d1 = deps1_d1 - deps_v_e_d1 - deps_v_ps_d1;
    double deps_q_pc_d1 = (2.0 / 3.0) * deps1_d1 - deps_q_e_d1 - deps_q_ps_d1;

    double Av_eff = std::max(deps_v_pc_d1, 1.0e-3 / Eoed_ref);
    double Aq_eff = std::max(deps_q_pc_d1, 1.0e-3 / Eoed_ref);

    if (alpha_cap <= 0.0) {
      alpha_cap = std::sqrt(eta_k0 * Av_eff / Aq_eff);
      alpha_cap = std::max(0.20, std::min(5.0, alpha_cap));
    }

    if (beta_cap <= 0.0) {
      double f_k0 = std::sqrt(1.0 + (eta_k0 / alpha_cap) * (eta_k0 / alpha_cap));
      double kpp  = ((1.0 + 2.0 * K0nc) / 3.0) * f_k0;
      beta_cap    = Av_eff * std::pow(kpp, m - 1.0);
    }
  }
};

struct HSState {
  double sig[6]         = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eps[6]         = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double gamma_p        = 0.0;
  double p_p            = 100.0;
  double void_ratio     = 0.50;
  double e_dev[3][3]    = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double H_hist[3][3]   = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  bool   is_virgin      = true;
  double gamma_hist_max = 0.0;
};

// ---------------------------------------------------------------------------
//  Symmetric 3x3 Jacobi Eigensolver (for Benz 2007 HSsmall strain history)
// ---------------------------------------------------------------------------
inline void sym_eigen_3x3(const double A[3][3], double vals[3], double vecs[3][3]) {
  double D[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      D[i][j]    = A[i][j];
      vecs[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int iter = 0; iter < 50; ++iter) {
    // Find largest off-diagonal entry:
    int p = 0, q = 1;
    double max_off = std::abs(D[0][1]);
    if (std::abs(D[0][2]) > max_off) { max_off = std::abs(D[0][2]); p = 0; q = 2; }
    if (std::abs(D[1][2]) > max_off) { max_off = std::abs(D[1][2]); p = 1; q = 2; }
    if (max_off < 1.0e-16 * (std::abs(D[0][0]) + std::abs(D[1][1]) + std::abs(D[2][2]) + 1.0e-30)) {
      break;
    }

    double app = D[p][p];
    double aqq = D[q][q];
    double apq = D[p][q];
    double tau = (aqq - app) / (2.0 * apq);
    double t   = (tau >= 0.0) ? (1.0 / (tau + std::sqrt(1.0 + tau * tau)))
                              : (-1.0 / (-tau + std::sqrt(1.0 + tau * tau)));
    double c   = 1.0 / std::sqrt(1.0 + t * t);
    double s   = t * c;

    D[p][p] = app - t * apq;
    D[q][q] = aqq + t * apq;
    D[p][q] = D[q][p] = 0.0;

    int r = 3 - p - q;
    double dpr = D[p][r];
    double dqr = D[q][r];
    D[p][r] = D[r][p] = c * dpr - s * dqr;
    D[q][r] = D[r][q] = s * dpr + c * dqr;

    for (int k = 0; k < 3; ++k) {
      double vkp = vecs[k][p];
      double vkq = vecs[k][q];
      vecs[k][p] = c * vkp - s * vkq;
      vecs[k][q] = s * vkp + c * vkq;
    }
  }

  vals[0] = D[0][0];
  vals[1] = D[1][1];
  vals[2] = D[2][2];

  // Sort eigenvalues ascending for deterministic ordering:
  for (int i = 0; i < 2; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      if (vals[j] < vals[i]) {
        std::swap(vals[i], vals[j]);
        for (int k = 0; k < 3; ++k) {
          std::swap(vecs[k][i], vecs[k][j]);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  Lode-Angle Shape Functions Theta(x), k3(x), Theta_sig(x) and Derivatives
// ---------------------------------------------------------------------------
inline void lode_theta_and_deriv(const HSParams& p, double x_sin3, double& Theta, double& dTheta_dx) {
  double x = std::max(-1.0, std::min(1.0, x_sin3));
  if (p.lode_shape == LODE_MN) {
    double u     = -p.gamma_mn * x;
    double phi3  = std::acos(u) / 3.0;
    Theta        = std::cos(phi3) / p.c0_mn;
    double denom = std::sqrt(std::max(1.0e-30, 1.0 - u * u));
    dTheta_dx    = -(p.gamma_mn / (3.0 * p.c0_mn)) * std::sin(phi3) / denom;
  } else {
    double u          = p.gamma_mc * x;
    double th_bar     = std::asin(u) / 3.0;
    double dth_dx     = (p.gamma_mc / 3.0) / std::sqrt(std::max(1.0e-30, 1.0 - u * u));
    double s_norm     = std::sin(th_bar) / std::sin(p.th_max_mc);
    double ds_dth     = std::cos(th_bar) / std::sin(p.th_max_mc);
    double c_ratio    = std::cos(th_bar) / std::cos(p.th_max_mc);
    double dc_dth     = -std::sin(th_bar) / std::cos(p.th_max_mc);
    double d          = p.delta_lode;
    Theta             = 0.5 * ((d + 1.0) * c_ratio - (d - 1.0) * s_norm);
    double dTheta_dth = 0.5 * ((d + 1.0) * dc_dth  - (d - 1.0) * ds_dth);
    dTheta_dx         = dTheta_dth * dth_dx;
  }
}

inline void lode_factors_and_derivs(
    const HSParams& p,
    double x_sin3,
    double& Theta,
    double& dTheta_dx,
    double& k3,
    double& dk3_dx,
    double& Theta_sig,
    double& dTheta_sig_dx) {
  lode_theta_and_deriv(p, x_sin3, Theta, dTheta_dx);
  double x       = std::max(-1.0, std::min(1.0, x_sin3));
  double u       = p.gamma_mc * x;
  double th_bar  = std::asin(u) / 3.0;
  double dth_dx  = (p.gamma_mc / 3.0) / std::sqrt(std::max(1.0e-30, 1.0 - u * u));
  double s_norm  = std::sin(th_bar) / std::sin(p.th_max_mc);
  double ds_dth  = std::cos(th_bar) / std::sin(p.th_max_mc);
  double c_ratio = std::cos(th_bar) / std::cos(p.th_max_mc);
  double dc_dth  = -std::sin(th_bar) / std::cos(p.th_max_mc);

  // Endpoint-exact C^infty minor principal stress factor:
  //   sig_3 = p - k3(x) * q, where k3(+1) = 1/3 (TXC) and k3(-1) = 2/3 (TXE)
  k3     = 0.5 * (c_ratio - (1.0 / 3.0) * s_norm);
  dk3_dx = 0.5 * (dc_dth  - (1.0 / 3.0) * ds_dth) * dth_dx;

  // Principal-stress-difference factor q_sig = Theta_sig(x) * q:
  //   For LODE_MC, Theta_sig(x) == c_ratio(x) so q_sig == sig_1 - sig_3 identically.
  //   For LODE_MN, Theta_sig(+1) = Theta_sig(-1) = 1 and q_sig = q_f(sig_3) on f^f = 0.
  double alpha_phi = 3.0 * (1.0 - p.sin_phi) / (3.0 - p.sin_phi);
  Theta_sig        = (Theta     - p.M_c * k3)     / alpha_phi;
  dTheta_sig_dx    = (dTheta_dx - p.M_c * dk3_dx) / alpha_phi;
}

inline void invariants_and_lode_grad(
    const HSParams& p,
    const double sig[6],
    double& p_mean,
    double& q,
    double& x_sin3,
    double& Theta,
    double dq_dsig[6],
    double dTheta_dsig[6],
    double* dx_dsig_out = nullptr) {
  p_mean = (sig[0] + sig[1] + sig[2]) / 3.0;
  double s0 = sig[0] - p_mean;
  double s1 = sig[1] - p_mean;
  double s2 = sig[2] - p_mean;
  double s3 = sig[3];
  double s4 = sig[4];
  double s5 = sig[5];

  double J2 = 0.5 * (s0 * s0 + s1 * s1 + s2 * s2) + s3 * s3 + s4 * s4 + s5 * s5;
  q = std::sqrt(std::max(0.0, 3.0 * J2));
  if (q < 1.0e-12) {
    x_sin3 = 1.0;
    Theta  = 1.0;
    for (int i = 0; i < 6; ++i) {
      dq_dsig[i]     = 0.0;
      dTheta_dsig[i] = 0.0;
      if (dx_dsig_out) dx_dsig_out[i] = 0.0;
    }
    return;
  }

  double inv_2q = 1.5 / q;
  dq_dsig[0] = inv_2q * s0;
  dq_dsig[1] = inv_2q * s1;
  dq_dsig[2] = inv_2q * s2;
  dq_dsig[3] = 3.0 * s3 / q;
  dq_dsig[4] = 3.0 * s4 / q;
  dq_dsig[5] = 3.0 * s5 / q;

  double J3 = s0 * (s1 * s2 - s4 * s4)
            - s3 * (s3 * s2 - s4 * s5)
            + s5 * (s3 * s4 - s1 * s5);

  double t00  = s0 * s0 + s3 * s3 + s5 * s5;
  double t11  = s3 * s3 + s1 * s1 + s4 * s4;
  double t22  = s5 * s5 + s4 * s4 + s2 * s2;
  double t01  = s0 * s3 + s3 * s1 + s5 * s4;
  double t12  = s3 * s5 + s1 * s4 + s4 * s2;
  double t20  = s5 * s0 + s4 * s3 + s2 * s5;
  double tr_t = (t00 + t11 + t22) / 3.0;

  double dJ3_dsig[6] = {
      t00 - tr_t,
      t11 - tr_t,
      t22 - tr_t,
      2.0 * t01,
      2.0 * t12,
      2.0 * t20
  };

  double q3    = q * q * q;
  double x_raw = 13.5 * J3 / q3;
  x_sin3       = std::max(-1.0, std::min(1.0, x_raw));
  double dTheta_dx = 0.0;
  lode_theta_and_deriv(p, x_sin3, Theta, dTheta_dx);

  for (int i = 0; i < 6; ++i) {
    double dx_dsig_i = (13.5 / q3) * dJ3_dsig[i] - (3.0 * x_sin3 / q) * dq_dsig[i];
    dTheta_dsig[i]   = dTheta_dx * dx_dsig_i;
    if (dx_dsig_out) dx_dsig_out[i] = dx_dsig_i;
  }
}

// ---------------------------------------------------------------------------
//  Constitutive Helper Functions in 3D Haigh-Westergaard Space (p, q, x)
// ---------------------------------------------------------------------------
inline void stiffness_scale_and_derivs(
    const HSParams& p, double p_mean, double q,
    double& s, double& ds_dp, double& ds_dq,
    double k3 = 1.0 / 3.0, double* ds_dk3_out = nullptr) {
  double p_hat = p_mean + p.a_apex;
  if (p.stiff_stress == STIFF_SIGMA3) {
    double sig3_hat = p_hat - k3 * q;
    if (sig3_hat <= p.p_min) {
      s     = std::pow(p.p_min / p.p_hat_ref, p.m);
      ds_dp = 0.0;
      ds_dq = 0.0;
      if (ds_dk3_out) *ds_dk3_out = 0.0;
      return;
    }
    double ratio  = sig3_hat / p.p_hat_ref;
    s             = std::pow(ratio, p.m);
    double ds_ds3 = (p.m > 0.0) ? ((p.m / p.p_hat_ref) * std::pow(ratio, p.m - 1.0)) : 0.0;
    ds_dp         = ds_ds3;
    ds_dq         = -k3 * ds_ds3;
    if (ds_dk3_out) *ds_dk3_out = -q * ds_ds3;
  } else {
    if (p_hat <= p.p_min) {
      s     = std::pow(p.p_min / p.p_hat_ref, p.m);
      ds_dp = 0.0;
      ds_dq = 0.0;
      if (ds_dk3_out) *ds_dk3_out = 0.0;
      return;
    }
    double ratio = p_hat / p.p_hat_ref;
    s            = std::pow(ratio, p.m);
    ds_dp        = (p.m > 0.0) ? ((p.m / p.p_hat_ref) * std::pow(ratio, p.m - 1.0)) : 0.0;
    ds_dq        = 0.0;
    if (ds_dk3_out) *ds_dk3_out = 0.0;
  }
}

inline void phi_shear(
    const HSParams& p, double p_mean, double q,
    double& Phi_s, double& dPhi_dp, double& dPhi_dq,
    double k3 = 1.0 / 3.0, double Theta_sig = 1.0,
    double dk3_dx = 0.0, double dTheta_sig_dx = 0.0,
    double* dPhi_dx_out = nullptr) {
  double p_hat     = std::max(p.p_min, p_mean + p.a_apex);
  double s3_raw    = p_hat - k3 * q;
  double sig3_hat  = std::max(p.p_min, s3_raw);
  double dsig3_dp  = (s3_raw > p.p_min && (p_mean + p.a_apex) > p.p_min) ? 1.0 : 0.0;
  double dsig3_dq  = (s3_raw > p.p_min) ? -k3 : 0.0;
  double dsig3_dk3 = (s3_raw > p.p_min) ? -q  : 0.0;

  double q_sig      = Theta_sig * q;
  double dqsig_dq   = Theta_sig;
  double dqsig_dThs = q;

  double kf      = 2.0 * p.sin_phi / (1.0 - p.sin_phi);
  double qa      = (kf / p.Rf) * sig3_hat;
  double dqa_ds3 = kf / p.Rf;

  double r_raw    = q_sig / qa;
  double r        = std::min(0.9999, std::max(0.0, r_raw));
  double dr_ds3   = 0.0;
  double dr_dqsig = 0.0;
  if (r_raw > 0.0 && r_raw < 0.9999) {
    dr_ds3   = -(q_sig / (qa * qa)) * dqa_ds3;
    dr_dqsig = 1.0 / qa;
  }

  double one_m_r   = 1.0 - r;
  double c_hyp_raw = 2.0 / (p.Ei_ref * one_m_r) - 2.0 / p.Eur_ref;
  double c_hyp     = std::max(0.0, c_hyp_raw);
  double dc_dr     = (c_hyp_raw > 0.0) ? (2.0 / (p.Ei_ref * one_m_r * one_m_r)) : 0.0;

  double s = 1.0, ds_dp = 0.0, ds_dq = 0.0, ds_dk3 = 0.0;
  stiffness_scale_and_derivs(p, p_mean, q, s, ds_dp, ds_dq, k3, &ds_dk3);

  Phi_s = c_hyp * q_sig / s;

  double dc_dp = dc_dr * (dr_ds3 * dsig3_dp);
  double dc_dq = dc_dr * (dr_ds3 * dsig3_dq + dr_dqsig * dqsig_dq);
  dPhi_dp      = (dc_dp * q_sig) / s - (Phi_s / s) * ds_dp;
  dPhi_dq      = (dc_dq * q_sig + c_hyp * dqsig_dq) / s - (Phi_s / s) * ds_dq;

  if (dPhi_dx_out) {
    double dsig3_dx = dsig3_dk3 * dk3_dx;
    double dqsig_dx = dqsig_dThs * dTheta_sig_dx;
    double ds_dx    = ds_dk3 * dk3_dx;
    double dc_dx    = dc_dr * (dr_ds3 * dsig3_dx + dr_dqsig * dqsig_dx);
    *dPhi_dx_out    = (dc_dx * q_sig + c_hyp * dqsig_dx) / s - (Phi_s / s) * ds_dx;
  }
}

inline void rowe_dilatancy(
    const HSParams& p, double p_mean, double q, double void_ratio,
    double& s_psim, double& dspsim_dp, double& dspsim_dq,
    double k3 = 1.0 / 3.0, double Theta_sig = 1.0,
    double dk3_dx = 0.0, double dTheta_sig_dx = 0.0,
    double* dspsim_dx_out = nullptr) {
  double p_hat     = std::max(p.p_min, p_mean + p.a_apex);
  double s3_raw    = p_hat - k3 * q;
  double sig3_hat  = std::max(p.p_min, s3_raw);
  double dsig3_dp  = (s3_raw > p.p_min && (p_mean + p.a_apex) > p.p_min) ? 1.0 : 0.0;
  double dsig3_dq  = (s3_raw > p.p_min) ? -k3 : 0.0;
  double dsig3_dk3 = (s3_raw > p.p_min) ? -q  : 0.0;

  double q_sig      = Theta_sig * q;
  double dqsig_dq   = Theta_sig;
  double dqsig_dThs = q;

  double denom_raw    = 2.0 * sig3_hat + q_sig;
  double denom        = std::max(p.p_min, denom_raw);
  double sp_raw       = q_sig / denom;
  double s_phim       = std::min(0.995, std::max(0.0, sp_raw));
  double dsphim_ds3   = 0.0;
  double dsphim_dqsig = 0.0;
  if (sp_raw > 0.0 && sp_raw < 0.995 && denom_raw > p.p_min) {
    dsphim_ds3   = -2.0 * q_sig / (denom * denom);
    dsphim_dqsig =  2.0 * sig3_hat / (denom * denom);
  }

  double scv           = p.sin_phicv;
  double denom_r       = 1.0 - s_phim * scv;
  s_psim               = (s_phim - scv) / denom_r;
  double dspsim_dsphim = (1.0 - scv * scv) / (denom_r * denom_r);

  if (p.emax > 0.0 && void_ratio >= p.emax && s_psim > 0.0) {
    s_psim    = 0.0;
    dspsim_dp = 0.0;
    dspsim_dq = 0.0;
    if (dspsim_dx_out) *dspsim_dx_out = 0.0;
    return;
  }
  dspsim_dp = dspsim_dsphim * (dsphim_ds3 * dsig3_dp);
  dspsim_dq = dspsim_dsphim * (dsphim_ds3 * dsig3_dq + dsphim_dqsig * dqsig_dq);
  if (dspsim_dx_out) {
    double dsig3_dx = dsig3_dk3 * dk3_dx;
    double dqsig_dx = dqsig_dThs * dTheta_sig_dx;
    *dspsim_dx_out  = dspsim_dsphim * (dsphim_ds3 * dsig3_dx + dsphim_dqsig * dqsig_dx);
  }
}

inline void cap_pp_and_deriv(
    const HSParams& p, double p_p_n, double deps_v_pc,
    double& p_p_new, double& dpp_ddeps) {
  double pp_hat_n = std::max(p.p_min, p_p_n + p.a_apex);
  if (deps_v_pc <= 0.0) {
    p_p_new   = p_p_n;
    dpp_ddeps = 0.0;
    return;
  }
  if (std::abs(1.0 - p.m) > 1.0e-8) {
    double om         = 1.0 - p.m;
    double base_n     = std::pow(pp_hat_n / p.p_hat_ref, om);
    double inc        = om * deps_v_pc / (p.beta_cap * p.p_hat_ref);
    double val        = std::max(1.0e-12, base_n + inc);
    double pp_hat_new = p.p_hat_ref * std::pow(val, 1.0 / om);
    dpp_ddeps         = (1.0 / p.beta_cap) * std::pow(pp_hat_new / p.p_hat_ref, p.m);
    p_p_new           = pp_hat_new - p.a_apex;
  } else {
    double arg        = deps_v_pc / (p.beta_cap * p.p_hat_ref);
    double pp_hat_new = pp_hat_n * std::exp(std::min(50.0, arg));
    dpp_ddeps         = pp_hat_new / (p.beta_cap * p.p_hat_ref);
    p_p_new           = pp_hat_new - p.a_apex;
  }
}

inline void initialize_state_from_stress(
    const HSParams& p, const double sig0_comp[6], HSState& st,
    double ocr = -1.0, double pop = -1.0) {
  for (int i = 0; i < 6; ++i) {
    st.sig[i] = sig0_comp[i];
    st.eps[i] = 0.0;
  }
  st.void_ratio = p.e0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      st.e_dev[i][j]  = 0.0;
      st.H_hist[i][j] = 0.0;
    }
  }
  st.is_virgin      = true;
  st.gamma_hist_max = 0.0;

  double p_mean = 0.0, q = 0.0, x_sin3 = 1.0, Theta = 1.0, dTheta_dx = 0.0;
  double k3 = 1.0 / 3.0, dk3_dx = 0.0, Theta_sig = 1.0, dTheta_sig_dx = 0.0;
  double dq_dsig[6], dTheta_dsig[6];
  invariants_and_lode_grad(p, st.sig, p_mean, q, x_sin3, Theta, dq_dsig, dTheta_dsig);
  lode_factors_and_derivs(p, x_sin3, Theta, dTheta_dx, k3, dk3_dx, Theta_sig, dTheta_sig_dx);
  double q_tilde = Theta * q;
  double p_hat   = std::max(p.p_min, p_mean + p.a_apex);

  double p_eq    = std::sqrt(p_hat * p_hat + (q / p.alpha_cap) * (q / p.alpha_cap)) - p.a_apex;
  double ocr_use = (ocr <= 0.0) ? p.OCR : ocr;
  double pop_use = (pop <  0.0) ? p.pop : pop;
  st.p_p         = std::max(p.p_min, std::max(p_eq * ocr_use, p_eq + pop_use));

  if (q > 1.0e-12) {
    double Phi_s = 0.0, dPhi_dp = 0.0, dPhi_dq = 0.0;
    phi_shear(p, p_mean, q, Phi_s, dPhi_dp, dPhi_dq, k3, Theta_sig);
    st.gamma_p = std::max(0.0, Phi_s);
  } else {
    st.gamma_p = 0.0;
  }
}

// ---------------------------------------------------------------------------
//  Benz (2007) HSsmall Paraelastic Predictor & Exact Fourth-Order Tangent
// ---------------------------------------------------------------------------
inline void phi_hssmall_integral(
    double G0, double Gur, double c_h, double gam,
    double& phi_val, double& gt_val) {
  if (G0 <= Gur || c_h <= 0.0) {
    phi_val = Gur * gam;
    gt_val  = Gur;
    return;
  }
  double gam_cut = (std::sqrt(G0 / Gur) - 1.0) / c_h;
  if (gam <= gam_cut) {
    double op = 1.0 + c_h * gam;
    phi_val   = G0 * gam / op;
    gt_val    = G0 / (op * op);
  } else {
    double op_c = 1.0 + c_h * gam_cut;
    phi_val     = G0 * gam_cut / op_c + Gur * (gam - gam_cut);
    gt_val      = Gur;
  }
}

inline void build_isotropic_elastic_tangent(double Kur, double Gur, double C_e[6][6]) {
  static constexpr double I6_M[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
  static constexpr double IDEV6_DIAG[6] = {2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 0.5, 0.5, 0.5};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double id_ij = 0.0;
      if (i == j) {
        id_ij = IDEV6_DIAG[i];
      } else if (i < 3 && j < 3) {
        id_ij = -1.0 / 3.0;
      }
      C_e[i][j] = Kur * I6_M[i] * I6_M[j] + 2.0 * Gur * id_ij;
    }
  }
}

inline void hssmall_predictor(
    const HSParams& p,
    const HSState& st_n,
    const double eps_new[6],
    double s_n,
    double& G_sec,
    double e_dev_new[3][3],
    double H_new[3][3],
    bool& virgin_new,
    double C_e[6][6],
    double dGe_deps[6],
    double* gamma_hist_max_out = nullptr) {
  static constexpr double I6_M[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
  double Gur      = p.Gur_ref * s_n;
  double Kur      = p.Kur_ref * s_n;
  double kappa_KG = p.Kur_ref / p.Gur_ref;

  double eps_v_new = eps_new[0] + eps_new[1] + eps_new[2];
  double eps_v_n   = st_n.eps[0] + st_n.eps[1] + st_n.eps[2];
  double deps_v    = eps_v_new - eps_v_n;

  e_dev_new[0][0] = eps_new[0] - eps_v_new / 3.0;
  e_dev_new[1][1] = eps_new[1] - eps_v_new / 3.0;
  e_dev_new[2][2] = eps_new[2] - eps_v_new / 3.0;
  e_dev_new[0][1] = e_dev_new[1][0] = 0.5 * eps_new[3];
  e_dev_new[1][2] = e_dev_new[2][1] = 0.5 * eps_new[4];
  e_dev_new[0][2] = e_dev_new[2][0] = 0.5 * eps_new[5];

  for (int i = 0; i < 6; ++i) dGe_deps[i] = 0.0;

  if (!p.use_small_strain || p.G0_ref <= p.Gur_ref) {
    G_sec      = Gur;
    virgin_new = st_n.is_virgin;
    if (gamma_hist_max_out) *gamma_hist_max_out = st_n.gamma_hist_max;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        H_new[i][j] = st_n.H_hist[i][j];
    build_isotropic_elastic_tangent(Kur, Gur, C_e);
    return;
  }

  double de_dev[3][3];
  double sum_edev_sq = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      de_dev[i][j] = e_dev_new[i][j] - st_n.e_dev[i][j];
      sum_edev_sq += e_dev_new[i][j] * e_dev_new[i][j];
    }
  }
  double gam_dev_new = std::sqrt(std::max(0.0, 2.0 * sum_edev_sq));

  double G0 = p.G0_ref * s_n;
  double vals[3];
  double vecs[3][3];
  sym_eigen_3x3(de_dev, vals, vecs);

  // H_hat = vecs^T * H_hist * vecs:
  double H_hat[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      double sum = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          sum += vecs[i][a] * st_n.H_hist[i][j] * vecs[j][b];
      H_hat[a][b] = sum;
    }
  }

  double w[3] = {1.0, 1.0, 1.0};
  bool reversed_any = false;
  for (int a = 0; a < 3; ++a) {
    if (vals[a] * H_hat[a][a] < -1.0e-16) {
      w[a] = 0.0;
      if (std::abs(vals[a]) > 1.0e-12 && std::abs(H_hat[a][a]) > 1.0e-12) {
        reversed_any = true;
      }
    }
  }

  virgin_new = st_n.is_virgin && (!reversed_any);
  double H0_hat[3][3], H1_hat[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      H0_hat[a][b] = std::sqrt(w[a] * w[b]) * H_hat[a][b];
      H1_hat[a][b] = H0_hat[a][b] + ((a == b) ? vals[a] : 0.0);
    }
  }

  double H0[3][3];
  double sum_H0_sq = 0.0;
  double sum_H1_sq = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double s0 = 0.0, s1 = 0.0;
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          s0 += vecs[i][a] * H0_hat[a][b] * vecs[j][b];
          s1 += vecs[i][a] * H1_hat[a][b] * vecs[j][b];
        }
      }
      H0[i][j]    = s0;
      H_new[i][j] = s1;
      sum_H0_sq  += s0 * s0;
      sum_H1_sq  += s1 * s1;
    }
  }

  double gam0 = std::sqrt(std::max(0.0, 2.0 * sum_H0_sq));
  double gam1 = std::sqrt(std::max(0.0, 2.0 * sum_H1_sq));
  double dgam = gam1 - gam0;

  double tr_H1 = (H_new[0][0] + H_new[1][1] + H_new[2][2]) / 3.0;
  double H1_dev_voigt[6] = {
      H_new[0][0] - tr_H1,
      H_new[1][1] - tr_H1,
      H_new[2][2] - tr_H1,
      H_new[0][1],
      H_new[1][2],
      H_new[0][2]
  };
  double edev_voigt[6] = {
      e_dev_new[0][0], e_dev_new[1][1], e_dev_new[2][2],
      e_dev_new[0][1], e_dev_new[1][2], e_dev_new[0][2]
  };

  // Extended Masing rule (Masing's 2nd rule): when reloading exceeds the historical
  // maximum deviatoric strain norm (gam_dev_new > gamma_hist_max), split the step at
  // gamma_hist_max so the loop closes on Masing (m_fac=2) and rejoins the virgin
  // backbone (m_fac=1) beyond gamma_hist_max:
  if (!virgin_new && st_n.gamma_hist_max > 0.0 &&
      gam_dev_new > st_n.gamma_hist_max + 1.0e-12 && dgam > 1.0e-12 && gam1 > 1.0e-14) {
    double c_h_m = 0.385 / (2.0 * std::max(1.0e-12, p.gamma_07));
    double c_h_v = 0.385 / std::max(1.0e-12, p.gamma_07);

    double gam_v0 = std::max(st_n.gamma_hist_max, gam_dev_new - dgam);
    double dgam_v = std::max(0.0, gam_dev_new - gam_v0);
    double dgam_m = std::max(0.0, dgam - dgam_v);
    double gam_m1 = gam0 + dgam_m;

    double phi_m1 = 0.0, gt_m1 = Gur, phi_m0 = 0.0, gt_m0 = Gur;
    double phi_v1 = 0.0, gt_v1 = Gur, phi_v0 = 0.0, gt_v0 = Gur;
    phi_hssmall_integral(G0, Gur, c_h_m, gam_m1,      phi_m1, gt_m1);
    phi_hssmall_integral(G0, Gur, c_h_m, gam0,        phi_m0, gt_m0);
    phi_hssmall_integral(G0, Gur, c_h_v, gam_dev_new, phi_v1, gt_v1);
    phi_hssmall_integral(G0, Gur, c_h_v, gam_v0,      phi_v0, gt_v0);

    double dphi = (phi_m1 - phi_m0) + (phi_v1 - phi_v0);
    G_sec       = std::max(Gur, std::min(G0, dphi / dgam));

    double coef_gam1 = 0.0, coef_gcur = 0.0;
    if (gam_dev_new - dgam < st_n.gamma_hist_max) {
      coef_gam1 = (gt_m1 - G_sec) / dgam;
      coef_gcur = (gt_v1 - gt_m1) / dgam;
    } else {
      coef_gam1 = (gt_v0 - G_sec) / dgam;
      coef_gcur = (gt_v1 - gt_v0) / dgam;
    }
    for (int i = 0; i < 6; ++i) {
      dGe_deps[i] = coef_gam1 * (2.0 / gam1) * H1_dev_voigt[i]
                  + coef_gcur * (2.0 / gam_dev_new) * edev_voigt[i];
    }

    virgin_new = true;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        H_new[i][j] = e_dev_new[i][j];
    if (gamma_hist_max_out) *gamma_hist_max_out = std::max(st_n.gamma_hist_max, gam_dev_new);

    double de_sig[6] = {
        de_dev[0][0], de_dev[1][1], de_dev[2][2],
        de_dev[0][1], de_dev[1][2], de_dev[0][2]
    };
    double Ke_sec = kappa_KG * G_sec;
    build_isotropic_elastic_tangent(Ke_sec, G_sec, C_e);
    for (int i = 0; i < 6; ++i) {
      double dsig_tr_dGe_i = kappa_KG * deps_v * I6_M[i] + 2.0 * de_sig[i];
      for (int j = 0; j < 6; ++j) {
        C_e[i][j] += dsig_tr_dGe_i * dGe_deps[j];
      }
    }
    return;
  }

  if (gamma_hist_max_out) {
    *gamma_hist_max_out = std::max(st_n.gamma_hist_max, gam_dev_new);
  }

  double m_fac = virgin_new ? 1.0 : 2.0;
  double c_h   = 0.385 / (m_fac * std::max(1.0e-12, p.gamma_07));

  double phi1 = 0.0, gt1 = Gur, phi0 = 0.0, gt0 = Gur;
  phi_hssmall_integral(G0, Gur, c_h, gam1, phi1, gt1);
  phi_hssmall_integral(G0, Gur, c_h, gam0, phi0, gt0);

  if (std::abs(dgam) > 1.0e-12 && gam1 > 1.0e-14) {
    G_sec = std::max(Gur, std::min(G0, (phi1 - phi0) / dgam));
    double dGsec_dgam1 = (gt1 - G_sec) / dgam;

    for (int i = 0; i < 6; ++i) {
      dGe_deps[i] = dGsec_dgam1 * (2.0 / gam1) * H1_dev_voigt[i];
    }

    double de_sig[6] = {
        de_dev[0][0], de_dev[1][1], de_dev[2][2],
        de_dev[0][1], de_dev[1][2], de_dev[0][2]
    };
    double Ke_sec = kappa_KG * G_sec;
    build_isotropic_elastic_tangent(Ke_sec, G_sec, C_e);
    for (int i = 0; i < 6; ++i) {
      double dsig_tr_dGe_i = kappa_KG * deps_v * I6_M[i] + 2.0 * de_sig[i];
      for (int j = 0; j < 6; ++j) {
        C_e[i][j] += dsig_tr_dGe_i * dGe_deps[j];
      }
    }
  } else {
    G_sec = std::max(Gur, std::min(G0, gt0));
    build_isotropic_elastic_tangent(kappa_KG * G_sec, G_sec, C_e);
  }
}

// ---------------------------------------------------------------------------
//  2x2 Linear Solver Helper: sens (2x4) = - J_y^{-1} * J_tr
// ---------------------------------------------------------------------------
inline bool solve_2x4_neg(const double J_y[2][2], const double J_tr[2][4], double sens[2][4]) {
  double det = J_y[0][0] * J_y[1][1] - J_y[0][1] * J_y[1][0];
  if (std::abs(det) < 1.0e-30) return false;
  double inv_det = -1.0 / det;
  for (int col = 0; col < 4; ++col) {
    sens[0][col] = inv_det * ( J_y[1][1] * J_tr[0][col] - J_y[0][1] * J_tr[1][col]);
    sens[1][col] = inv_det * (-J_y[1][0] * J_tr[0][col] + J_y[0][0] * J_tr[1][col]);
  }
  return true;
}

// ---------------------------------------------------------------------------
//  Regime 1: Single-Surface Hyperbolic Shear Hardening (f^s = 0, q_tilde < q_f)
// ---------------------------------------------------------------------------
inline bool solve_regime1_shear(
    const HSParams& p,
    double p_tr, double q_tr,
    double Theta, double k3, double dk3_dx, double Theta_sig, double dTheta_sig_dx,
    double Ke, double Ge, double gamma_p_n, double void_ratio,
    double& p_out, double& q_out, double& dgam_out, double sens[2][4]) {
  double inv_3Ge2   = 2.0 / (3.0 * Ge);
  double p_cur      = std::max(p.p_min - p.a_apex, p_tr);
  double q_max_fail = p.M_c * std::max(p.p_min, p_cur + p.a_apex) / Theta;
  double q_cur      = std::min(q_tr, 0.95 * q_max_fail);

  for (int iter = 0; iter < 35; ++iter) {
    double q_tilde = Theta * q_cur;
    double p_hat   = p_cur + p.a_apex;
    if (p_hat <= p.p_min || q_tilde >= p.M_c * p_hat) {
      return false;
    }

    double dgam = inv_3Ge2 * (q_tr - q_cur);
    double s_psim = 0.0, dspsi_dp = 0.0, dspsi_dq = 0.0, dspsi_dx = 0.0;
    rowe_dilatancy(p, p_cur, q_cur, void_ratio, s_psim, dspsi_dp, dspsi_dq,
                   k3, Theta_sig, dk3_dx, dTheta_sig_dx, &dspsi_dx);
    double Phi_s = 0.0, dPhi_dp = 0.0, dPhi_dq = 0.0, dPhi_dx = 0.0;
    phi_shear(p, p_cur, q_cur, Phi_s, dPhi_dp, dPhi_dq,
              k3, Theta_sig, dk3_dx, dTheta_sig_dx, &dPhi_dx);

    double R0 = p_cur - p_tr - Ke * s_psim * dgam;
    double R1 = Phi_s - gamma_p_n - dgam;

    double J_y[2][2] = {
        {1.0 - Ke * dspsi_dp * dgam, -Ke * (dspsi_dq * dgam - s_psim * inv_3Ge2)},
        {dPhi_dp,                    dPhi_dq + inv_3Ge2}
    };

    if (std::max(std::abs(R0) / std::max(1.0, p.p_ref), std::abs(R1)) < 1.0e-13) {
      if (dgam < -1.0e-12) return false;
      // Note: Since Ke = (Kur_ref / Gur_ref) * Ge scales proportionally with Ge,
      // the product Ke * dgam = (2/3)*(Kur_ref / Gur_ref)*(q_tr - q_cur) in R0
      // has zero partial derivative wrt Ge at fixed (p_tr, q_tr, x_tr).
      double J_tr[2][4] = {
          {-1.0, -Ke * s_psim * inv_3Ge2, -Ke * dspsi_dx * dgam, 0.0},
          { 0.0, -inv_3Ge2,                dPhi_dx,              dgam / Ge}
      };
      if (!solve_2x4_neg(J_y, J_tr, sens)) return false;
      p_out    = p_cur;
      q_out    = q_cur;
      dgam_out = std::max(0.0, dgam);
      return true;
    }

    double det = J_y[0][0] * J_y[1][1] - J_y[0][1] * J_y[1][0];
    if (std::abs(det) < 1.0e-30) return false;
    double dp = -( J_y[1][1] * R0 - J_y[0][1] * R1) / det;
    double dq = -(-J_y[1][0] * R0 + J_y[0][0] * R1) / det;

    double alpha_ls = 1.0;
    bool accepted   = false;
    for (int ls = 0; ls < 12; ++ls) {
      double p_next = p_cur + alpha_ls * dp;
      double q_next = std::max(0.0, std::min(q_tr, q_cur + alpha_ls * dq));
      if ((p_next + p.a_apex) > p.p_min &&
          (Theta * q_next) < 0.9995 * p.M_c * (p_next + p.a_apex)) {
        p_cur    = p_next;
        q_cur    = q_next;
        accepted = true;
        break;
      }
      alpha_ls *= 0.5;
    }
    if (!accepted) return false;
  }
  return false;
}

// ---------------------------------------------------------------------------
//  Regime 2: Ultimate Failure Truncation (f^f = 0 only)
// ---------------------------------------------------------------------------
inline bool solve_regime2_failure(
    const HSParams& p,
    double p_tr, double q_tr, double Theta, double dTheta_dx,
    double Ke, double Ge, double void_ratio,
    double& p_out, double& q_out, double& dgam_out, double sens[2][4]) {
  double s_psi = (p.emax > 0.0 && void_ratio >= p.emax && p.sin_psi > 0.0) ? 0.0 : p.sin_psi;
  double denom = 1.5 * Ge * Theta + p.M_c * Ke * s_psi;
  double dgam  = (Theta * q_tr - p.M_c * (p_tr + p.a_apex)) / denom;
  if (dgam < -1.0e-12) return false;
  dgam  = std::max(0.0, dgam);
  p_out = p_tr + Ke * s_psi * dgam;
  q_out = std::max(0.0, q_tr - 1.5 * Ge * dgam);
  if ((p_out + p.a_apex) <= p.p_min) return false;
  dgam_out = dgam;

  double ddgam_dptr = -p.M_c / denom;
  double ddgam_dqtr = Theta / denom;
  double ddgam_dTh  = (q_tr - 1.5 * Ge * dgam) / denom;

  // Since Ke = (Kur_ref / Gur_ref) * Ge, both Ke*dgam and Ge*dgam are independent of Ge
  // at fixed (p_tr, q_tr, Theta), so sens[0][3] = sens[1][3] = 0.0.
  sens[0][0] = 1.0 + Ke * s_psi * ddgam_dptr;
  sens[0][1] = Ke * s_psi * ddgam_dqtr;
  sens[0][2] = (Ke * s_psi * ddgam_dTh) * dTheta_dx;
  sens[0][3] = 0.0;

  sens[1][0] = -1.5 * Ge * ddgam_dptr;
  sens[1][1] = 1.0 - 1.5 * Ge * ddgam_dqtr;
  sens[1][2] = (-1.5 * Ge * ddgam_dTh) * dTheta_dx;
  sens[1][3] = 0.0;
  return true;
}

// ---------------------------------------------------------------------------
//  Regime 3: Elliptical Volumetric Cap Only (f^c = 0 only, Associative)
// ---------------------------------------------------------------------------
inline bool solve_regime3_cap(
    const HSParams& p,
    double p_tr, double q_tr, double Theta, double dTheta_dx,
    double Ke, double Ge, double p_p_n,
    double& p_out, double& q_out, double& pp_out, double sens[2][4]) {
  (void)Theta;
  (void)dTheta_dx;
  double p_cur  = std::max(p.p_min - p.a_apex, std::min(p_tr, p_p_n));
  double alpha2 = p.alpha_cap * p.alpha_cap;
  double c_q    = (3.0 * Ge) / alpha2;

  for (int iter = 0; iter < 35; ++iter) {
    double deps_v_pc = (p_tr - p_cur) / Ke;
    if (deps_v_pc < -1.0e-12) return false;
    double deps_pos = std::max(0.0, deps_v_pc);
    double p_hat    = std::max(p.p_min, p_cur + p.a_apex);

    double pp_new = p_p_n, dpp_ddeps = 0.0;
    cap_pp_and_deriv(p, p_p_n, deps_pos, pp_new, dpp_ddeps);
    double pp_hat = pp_new + p.a_apex;

    // Note: Since Ke = (Kur_ref / Gur_ref) * Ge, c_q * deps_pos = (3*Ge/(Ke*alpha^2))*(p_tr - p_cur)
    // is independent of Ge, while d(deps_pos)/dGe = -deps_pos / Ge in the cap hardening law.
    double denom_q = p_hat + c_q * deps_pos;
    double q_cur   = q_tr * p_hat / denom_q;
    double dq_dp   = q_tr * c_q * (p_tr + p.a_apex) / (Ke * denom_q * denom_q);
    double dq_dptr = -q_tr * p_hat * c_q / (Ke * denom_q * denom_q);
    double dq_dqtr = p_hat / denom_q;
    double dq_dGe  = 0.0;

    double R = (q_cur * q_cur) / alpha2 + p_hat * p_hat - pp_hat * pp_hat;
    double dR_dp = (2.0 * q_cur / alpha2) * dq_dp
                 + 2.0 * p_hat
                 + (2.0 * pp_hat / Ke) * dpp_ddeps;

    if (std::abs(R) / (p.p_hat_ref * p.p_hat_ref) < 1.0e-13) {
      double dR_dptr = (2.0 * q_cur / alpha2) * dq_dptr - (2.0 * pp_hat / Ke) * dpp_ddeps;
      double dR_dqtr = (2.0 * q_cur / alpha2) * dq_dqtr;
      double dR_dGe  = (2.0 * q_cur / alpha2) * dq_dGe
                     + (2.0 * pp_hat * deps_pos / Ge) * dpp_ddeps;

      sens[0][0] = -dR_dptr / dR_dp;
      sens[0][1] = -dR_dqtr / dR_dp;
      sens[0][2] = 0.0;
      sens[0][3] = -dR_dGe  / dR_dp;

      sens[1][0] = dq_dptr + dq_dp * sens[0][0];
      sens[1][1] = dq_dqtr + dq_dp * sens[0][1];
      sens[1][2] = 0.0;
      sens[1][3] = dq_dGe  + dq_dp * sens[0][3];

      p_out  = p_cur;
      q_out  = q_cur;
      pp_out = pp_new;
      return true;
    }

    if (std::abs(dR_dp) < 1.0e-30) return false;
    double dp = -R / dR_dp;
    p_cur = std::max(p.p_min - p.a_apex, std::min(p_tr, p_cur + dp));
  }
  return false;
}

// ---------------------------------------------------------------------------
//  Regime 4: Two-Surface Corner Return Map (Shear/Failure + Cap Active)
// ---------------------------------------------------------------------------
inline void eval_regime4_residuals(
    const HSParams& p,
    double p_val, double q_val,
    double p_tr_v, double q_tr_v,
    double Th_v, double k3_v, double Th_sig_v,
    double Ke_v, double Ge_v,
    double gamma_p_n, double p_p_n, double void_ratio, bool failure_mode,
    double& R0, double& R1, double& dgam_p, double& deps_v_pc, double& pp_new) {
  double alpha2   = p.alpha_cap * p.alpha_cap;
  double p_hat    = std::max(p.p_min, p_val + p.a_apex);
  double q_tilde  = Th_v * q_val;
  double deps_v_p = (p_tr_v - p_val) / Ke_v;
  double deps_q_p = (q_tr_v - q_val) / (3.0 * Ge_v);

  double s_psim = 0.0, d1 = 0.0, d2 = 0.0;
  if (failure_mode) {
    s_psim = (p.emax > 0.0 && void_ratio >= p.emax && p.sin_psi > 0.0) ? 0.0 : p.sin_psi;
  } else {
    rowe_dilatancy(p, p_val, q_val, void_ratio, s_psim, d1, d2, k3_v, Th_sig_v);
  }

  double k_cap = q_val / (alpha2 * p_hat);
  double denom = std::max(0.05, 0.5 + k_cap * s_psim);
  dgam_p       = (deps_q_p - k_cap * deps_v_p) / denom;
  deps_v_pc    = deps_v_p + s_psim * dgam_p;

  double dpp_ddeps = 0.0;
  cap_pp_and_deriv(p, p_p_n, std::max(0.0, deps_v_pc), pp_new, dpp_ddeps);
  double pp_hat = pp_new + p.a_apex;

  if (failure_mode) {
    R0 = q_tilde - p.M_c * p_hat;
  } else {
    double Phi_s = 0.0;
    phi_shear(p, p_val, q_val, Phi_s, d1, d2, k3_v, Th_sig_v);
    R0 = Phi_s - gamma_p_n - dgam_p;
  }
  R1 = (q_val * q_val) / alpha2 + p_hat * p_hat - pp_hat * pp_hat;
}

inline bool solve_regime4_corner(
    const HSParams& p,
    double p_tr, double q_tr,
    double Theta, double dTheta_dx, double k3, double dk3_dx, double Theta_sig, double dTheta_sig_dx,
    double Ke, double Ge, double gamma_p_n, double p_p_n, double void_ratio,
    bool failure_mode,
    double& p_out, double& q_out, double& dgam_s_out, double& pp_out, double sens[2][4]) {
  double p_cur    = std::max(p.p_min - p.a_apex, std::min(p_tr, std::max(p_p_n, 0.7 * p_tr)));
  double q_cur    = std::min(q_tr, 0.85 * p.M_c * (p_cur + p.a_apex) / Theta);
  double kappa_KG = p.Kur_ref / p.Gur_ref;

  for (int iter = 0; iter < 40; ++iter) {
    double R0 = 0.0, R1 = 0.0, dgam_p = 0.0, deps_v_pc = 0.0, pp_new = p_p_n;
    eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge,
                           gamma_p_n, p_p_n, void_ratio, failure_mode,
                           R0, R1, dgam_p, deps_v_pc, pp_new);

    double scale0 = failure_mode ? std::max(1.0, p.p_ref) : 1.0;
    double scale1 = p.p_hat_ref * p.p_hat_ref;
    if (std::max(std::abs(R0) / scale0, std::abs(R1) / scale1) < 1.0e-12) {
      if (dgam_p < -1.0e-10 || deps_v_pc < -1.0e-10) return false;
      double J_y[2][2]  = {{0.0, 0.0}, {0.0, 0.0}};
      double J_tr[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
      double hp  = std::max(1.0e-8, 1.0e-8 * std::abs(p_cur));
      double hq  = std::max(1.0e-8, 1.0e-8 * std::max(1.0, std::abs(q_cur)));
      double hx  = 1.0e-8;
      double hge = std::max(1.0e-6, 1.0e-8 * Ge);
      double r0p, r1p, r0m, r1m, g_tmp, e_tmp, pp_tmp;

      eval_regime4_residuals(p, p_cur + hp, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur - hp, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_y[0][0] = (r0p - r0m) / (2.0 * hp);
      J_y[1][0] = (r1p - r1m) / (2.0 * hp);

      eval_regime4_residuals(p, p_cur, q_cur + hq, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur, q_cur - hq, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_y[0][1] = (r0p - r0m) / (2.0 * hq);
      J_y[1][1] = (r1p - r1m) / (2.0 * hq);

      eval_regime4_residuals(p, p_cur, q_cur, p_tr + hp, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur, q_cur, p_tr - hp, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_tr[0][0] = (r0p - r0m) / (2.0 * hp);
      J_tr[1][0] = (r1p - r1m) / (2.0 * hp);

      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr + hq, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr - hq, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_tr[0][1] = (r0p - r0m) / (2.0 * hq);
      J_tr[1][1] = (r1p - r1m) / (2.0 * hq);

      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr, Theta + hx * dTheta_dx, k3 + hx * dk3_dx, Theta_sig + hx * dTheta_sig_dx, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr, Theta - hx * dTheta_dx, k3 - hx * dk3_dx, Theta_sig - hx * dTheta_sig_dx, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_tr[0][2] = (r0p - r0m) / (2.0 * hx);
      J_tr[1][2] = (r1p - r1m) / (2.0 * hx);

      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, kappa_KG * (Ge + hge), Ge + hge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
      eval_regime4_residuals(p, p_cur, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, kappa_KG * (Ge - hge), Ge - hge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
      J_tr[0][3] = (r0p - r0m) / (2.0 * hge);
      J_tr[1][3] = (r1p - r1m) / (2.0 * hge);

      if (!solve_2x4_neg(J_y, J_tr, sens)) return false;
      p_out      = p_cur;
      q_out      = q_cur;
      dgam_s_out = std::max(0.0, dgam_p);
      pp_out     = pp_new;
      return true;
    }

    double hp = std::max(1.0e-7, 1.0e-7 * std::abs(p_cur));
    double hq = std::max(1.0e-7, 1.0e-7 * std::max(1.0, std::abs(q_cur)));
    double r0p, r1p, r0m, r1m, g_tmp, e_tmp, pp_tmp;
    eval_regime4_residuals(p, p_cur + hp, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
    eval_regime4_residuals(p, p_cur - hp, q_cur, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
    double J00 = (r0p - r0m) / (2.0 * hp);
    double J10 = (r1p - r1m) / (2.0 * hp);

    eval_regime4_residuals(p, p_cur, q_cur + hq, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0p, r1p, g_tmp, e_tmp, pp_tmp);
    eval_regime4_residuals(p, p_cur, q_cur - hq, p_tr, q_tr, Theta, k3, Theta_sig, Ke, Ge, gamma_p_n, p_p_n, void_ratio, failure_mode, r0m, r1m, g_tmp, e_tmp, pp_tmp);
    double J01 = (r0p - r0m) / (2.0 * hq);
    double J11 = (r1p - r1m) / (2.0 * hq);

    double det = J00 * J11 - J01 * J10;
    if (std::abs(det) < 1.0e-30) return false;
    double dp = -( J11 * R0 - J01 * R1) / det;
    double dq = -(-J10 * R0 + J00 * R1) / det;

    double alpha_ls = 1.0;
    bool accepted   = false;
    for (int ls = 0; ls < 12; ++ls) {
      double p_next = p_cur + alpha_ls * dp;
      double q_next = std::max(0.0, std::min(q_tr, q_cur + alpha_ls * dq));
      if ((p_next + p.a_apex) > p.p_min) {
        if (failure_mode || (Theta * q_next) < 0.9995 * p.M_c * (p_next + p.a_apex)) {
          p_cur    = p_next;
          q_cur    = q_next;
          accepted = true;
          break;
        }
      }
      alpha_ls *= 0.5;
    }
    if (!accepted) return false;
  }
  return false;
}

// ---------------------------------------------------------------------------
//  Master 3D Constitutive Return Map (Compression-Positive Mechanics Sign)
// ---------------------------------------------------------------------------
inline ActiveRegime return_map_3d_comp(
    const HSParams& p,
    const HSState& st_n,
    const double eps_new[6],
    HSState& st_new,
    double C_alg[6][6]) {
  static constexpr double I6_M[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

  st_new = st_n;
  double deps[6];
  double norm_deps_sq = 0.0;
  for (int i = 0; i < 6; ++i) {
    if (!std::isfinite(eps_new[i])) {
      return REGIME_REFUSED;
    }
    st_new.eps[i] = eps_new[i];
    deps[i]       = eps_new[i] - st_n.eps[i];
    norm_deps_sq += deps[i] * deps[i];
  }
  // Refuse wild trial strain increments (|deps| > 5.0 = 500% strain in a single step):
  if (norm_deps_sq > 25.0) {
    return REGIME_REFUSED;
  }
  double deps_v = deps[0] + deps[1] + deps[2];

  // Void ratio update (compression-positive: deps_v > 0 decreases void ratio):
  st_new.void_ratio = st_n.void_ratio - (1.0 + st_n.void_ratio) * deps_v;

  // Committed reference stress scale s_n using true 3D minor principal stress factor k3_n:
  double norm_sig_n = 0.0;
  for (int i = 0; i < 6; ++i) norm_sig_n += st_n.sig[i] * st_n.sig[i];
  double s_n = 1.0;
  if (norm_sig_n > 1.0e-24) {
    double p_n = 0.0, q_n = 0.0, x_n = 1.0, Th_n = 1.0, dTh_dx_n = 0.0;
    double k3_n = 1.0 / 3.0, dk3_dx_n = 0.0, Th_sig_n = 1.0, dTh_sig_dx_n = 0.0;
    double dq_dummy[6], dTh_dummy[6];
    invariants_and_lode_grad(p, st_n.sig, p_n, q_n, x_n, Th_n, dq_dummy, dTh_dummy);
    lode_factors_and_derivs(p, x_n, Th_n, dTh_dx_n, k3_n, dk3_dx_n, Th_sig_n, dTh_sig_dx_n);
    double ds_dp = 0.0, ds_dq = 0.0;
    stiffness_scale_and_derivs(p, p_n, q_n, s_n, ds_dp, ds_dq, k3_n);
  }

  double Ge = p.Gur_ref * s_n;
  double C_e_alg[6][6];
  double dGe_deps[6];
  hssmall_predictor(
      p, st_n, st_new.eps, s_n, Ge,
      st_new.e_dev, st_new.H_hist, st_new.is_virgin,
      C_e_alg, dGe_deps, &st_new.gamma_hist_max);
  double Ke = (p.Kur_ref / p.Gur_ref) * Ge;

  double deps_dev_sig[6] = {
      deps[0] - deps_v / 3.0,
      deps[1] - deps_v / 3.0,
      deps[2] - deps_v / 3.0,
      0.5 * deps[3],
      0.5 * deps[4],
      0.5 * deps[5]
  };
  double sig_tr[6];
  for (int i = 0; i < 6; ++i) {
    sig_tr[i] = st_n.sig[i] + Ke * deps_v * I6_M[i] + 2.0 * Ge * deps_dev_sig[i];
  }

  double p_tr = 0.0, q_tr = 0.0, x_tr = 1.0, Theta = 1.0, dTheta_dx = 0.0;
  double k3 = 1.0 / 3.0, dk3_dx = 0.0, Theta_sig = 1.0, dTheta_sig_dx = 0.0;
  double dq_dsig_tr[6], dTheta_dsig_tr[6], dx_dsig_tr[6];
  invariants_and_lode_grad(p, sig_tr, p_tr, q_tr, x_tr, Theta, dq_dsig_tr, dTheta_dsig_tr, dx_dsig_tr);
  lode_factors_and_derivs(p, x_tr, Theta, dTheta_dx, k3, dk3_dx, Theta_sig, dTheta_sig_dx);

  // 1. Check Rankine / Apex Tension Cutoff (Surface 4):
  double p_cut = std::max(p.p_min - p.a_apex, -p.sigma_t);
  if (p_tr <= p_cut) {
    for (int i = 0; i < 6; ++i) {
      st_new.sig[i] = p_cut * I6_M[i] + 1.0e-3 * (sig_tr[i] - p_cut * I6_M[i]);
      for (int j = 0; j < 6; ++j) {
        C_alg[i][j] = 1.0e-3 * C_e_alg[i][j];
      }
    }
    return REGIME_TENSION;
  }

  // 2. Evaluate trial yield functions at (p_tr, q_tr, x_tr):
  double q_tilde_tr = Theta * q_tr;
  double p_hat_tr   = p_tr + p.a_apex;
  double f_fail_tr  = q_tilde_tr - p.M_c * p_hat_tr;
  double f_shear_tr = 1.0;
  if (f_fail_tr < 0.0) {
    double Phi_s_tr = 0.0, dPhi_dp = 0.0, dPhi_dq = 0.0;
    phi_shear(p, p_tr, q_tr, Phi_s_tr, dPhi_dp, dPhi_dq, k3, Theta_sig);
    f_shear_tr = Phi_s_tr - st_n.gamma_p;
  }

  double pp_hat_n = st_n.p_p + p.a_apex;
  double f_cap_tr = (q_tr * q_tr) / (p.alpha_cap * p.alpha_cap)
                  + p_hat_tr * p_hat_tr - pp_hat_n * pp_hat_n;

  if (f_shear_tr <= 1.0e-12 && f_fail_tr <= 1.0e-12 && f_cap_tr <= 1.0e-10) {
    for (int i = 0; i < 6; ++i) {
      st_new.sig[i] = sig_tr[i];
      for (int j = 0; j < 6; ++j) {
        C_alg[i][j] = C_e_alg[i][j];
      }
    }
    return REGIME_ELASTIC;
  }

  // 3. Koiter Active-Set Selection & Return Map in (p, q):
  ActiveRegime regime = REGIME_ELASTIC;
  double p_new = p_tr, q_new = q_tr;
  double sens[2][4] = {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}};
  bool solved = false;
  double alpha2 = p.alpha_cap * p.alpha_cap;

  // Note: e_max dilatancy cutoff is intentionally gated on the committed void ratio
  // st_n.void_ratio so that the step return map remains C^1-smooth with respect to
  // trial volumetric strain during global Newton-Raphson iterations.
  if (f_shear_tr > 1.0e-12 || f_fail_tr > 1.0e-12) {
    double p1 = 0.0, q1 = 0.0, dgam1 = 0.0, sens1[2][4];
    if (solve_regime1_shear(p, p_tr, q_tr, Theta, k3, dk3_dx, Theta_sig, dTheta_sig_dx,
                            Ke, Ge, st_n.gamma_p, st_n.void_ratio,
                            p1, q1, dgam1, sens1)) {
      double fc1 = (q1 * q1) / alpha2 + (p1 + p.a_apex) * (p1 + p.a_apex) - pp_hat_n * pp_hat_n;
      if (fc1 <= 1.0e-8) {
        p_new = p1; q_new = q1; st_new.gamma_p = st_n.gamma_p + dgam1;
        for (int r = 0; r < 2; ++r) for (int c = 0; c < 4; ++c) sens[r][c] = sens1[r][c];
        regime = REGIME_SHEAR;
        solved = true;
      }
    }
    if (!solved) {
      double p2 = 0.0, q2 = 0.0, dgam2 = 0.0, sens2[2][4];
      if (solve_regime2_failure(p, p_tr, q_tr, Theta, dTheta_dx, Ke, Ge, st_n.void_ratio,
                                p2, q2, dgam2, sens2)) {
        double fc2 = (q2 * q2) / alpha2 + (p2 + p.a_apex) * (p2 + p.a_apex) - pp_hat_n * pp_hat_n;
        if (fc2 <= 1.0e-8) {
          p_new = p2; q_new = q2; st_new.gamma_p = st_n.gamma_p + dgam2;
          for (int r = 0; r < 2; ++r) for (int c = 0; c < 4; ++c) sens[r][c] = sens2[r][c];
          regime = REGIME_FAILURE;
          solved = true;
        }
      }
    }
  }

  if (!solved && f_cap_tr > 1.0e-10) {
    double p3 = 0.0, q3 = 0.0, pp3 = 0.0, sens3[2][4];
    if (solve_regime3_cap(p, p_tr, q_tr, Theta, dTheta_dx, Ke, Ge, st_n.p_p,
                          p3, q3, pp3, sens3)) {
      double q_tilde3 = Theta * q3;
      double ff3      = q_tilde3 - p.M_c * (p3 + p.a_apex);
      double fs3      = 1.0;
      if (ff3 < 0.0) {
        double Phi3 = 0.0, d1 = 0.0, d2 = 0.0;
        phi_shear(p, p3, q3, Phi3, d1, d2, k3, Theta_sig);
        fs3 = Phi3 - st_n.gamma_p;
      }
      if (fs3 <= 1.0e-13 && ff3 <= 1.0e-13) {
        p_new = p3; q_new = q3; st_new.p_p = pp3;
        for (int r = 0; r < 2; ++r) for (int c = 0; c < 4; ++c) sens[r][c] = sens3[r][c];
        regime = REGIME_CAP;
        solved = true;
      }
    }
  }

  if (!solved) {
    double p4 = 0.0, q4 = 0.0, dgam4 = 0.0, pp4 = 0.0, sens4[2][4];
    if (solve_regime4_corner(p, p_tr, q_tr, Theta, dTheta_dx, k3, dk3_dx, Theta_sig, dTheta_sig_dx,
                             Ke, Ge, st_n.gamma_p, st_n.p_p, st_n.void_ratio,
                             false, p4, q4, dgam4, pp4, sens4)) {
      p_new = p4; q_new = q4; st_new.gamma_p = st_n.gamma_p + dgam4; st_new.p_p = pp4;
      for (int r = 0; r < 2; ++r) for (int c = 0; c < 4; ++c) sens[r][c] = sens4[r][c];
      regime = REGIME_SHEAR_CAP;
      solved = true;
    } else if (solve_regime4_corner(p, p_tr, q_tr, Theta, dTheta_dx, k3, dk3_dx, Theta_sig, dTheta_sig_dx,
                                    Ke, Ge, st_n.gamma_p, st_n.p_p, st_n.void_ratio,
                                    true, p4, q4, dgam4, pp4, sens4)) {
      p_new = p4; q_new = q4; st_new.gamma_p = st_n.gamma_p + dgam4; st_new.p_p = pp4;
      for (int r = 0; r < 2; ++r) for (int c = 0; c < 4; ++c) sens[r][c] = sens4[r][c];
      regime = REGIME_FAILURE_CAP;
      solved = true;
    }
  }

  if (!solved) {
    return REGIME_REFUSED;
  }

  // 4. Reconstruct 3D stress tensor & 6x6 Consistent Algorithmic Tangent:
  double s_tr[6];
  for (int i = 0; i < 6; ++i) s_tr[i] = sig_tr[i] - p_tr * I6_M[i];

  if (q_tr > 1.0e-12) {
    double scale_q = q_new / q_tr;
    double n_s[6];
    for (int i = 0; i < 6; ++i) {
      st_new.sig[i] = p_new * I6_M[i] + scale_q * s_tr[i];
      n_s[i]        = s_tr[i] / q_tr;
    }

    double grad_p[6], grad_q[6];
    for (int j = 0; j < 6; ++j) {
      double dp_tr_j = I6_M[j] / 3.0;
      grad_p[j] = sens[0][0] * dp_tr_j + sens[0][1] * dq_dsig_tr[j] + sens[0][2] * dx_dsig_tr[j];
      grad_q[j] = sens[1][0] * dp_tr_j + sens[1][1] * dq_dsig_tr[j] + sens[1][2] * dx_dsig_tr[j];
    }

    double dsig_dsig_tr[6][6];
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double i_dev_ij = ((i == j) ? 1.0 : 0.0) - (I6_M[i] * I6_M[j]) / 3.0;
        dsig_dsig_tr[i][j] = I6_M[i] * grad_p[j]
                           + n_s[i]  * grad_q[j]
                           + scale_q * (i_dev_ij - n_s[i] * dq_dsig_tr[j]);
      }
    }

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) sum += dsig_dsig_tr[i][k] * C_e_alg[k][j];
        C_alg[i][j] = sum
                    + sens[0][3] * I6_M[i] * dGe_deps[j]
                    + sens[1][3] * n_s[i]  * dGe_deps[j];
      }
    }
  } else {
    for (int i = 0; i < 6; ++i) st_new.sig[i] = p_new * I6_M[i];
    double grad_p[6];
    for (int j = 0; j < 6; ++j) grad_p[j] = sens[0][0] * (I6_M[j] / 3.0);
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 6; ++k) sum += (I6_M[i] * grad_p[k]) * C_e_alg[k][j];
        C_alg[i][j] = sum + sens[0][3] * I6_M[i] * dGe_deps[j];
      }
    }
  }

  return regime;
}

}  // namespace LadrunoHSKernel

#endif  // LadrunoHSKernel_h
