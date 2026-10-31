// Standalone C++17 verification harness for SRC/material/nD/LadrunoHSKernel.h.
// Compiles with zero OpenSees dependencies:
//   g++ -O2 -std=c++17 tests/ladruno_hs_kernel_check.cpp -I SRC/material/nD -o ladruno_hs_kernel_check
// or:
//   cl /EHsc /O2 /std:c++17 tests\ladruno_hs_kernel_check.cpp /I SRC\material\nD

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "LadrunoHSKernel.h"

using namespace LadrunoHSKernel;

static double frob_norm_diff(const double A[6][6], const double B[6][6]) {
  double s = 0.0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      double d = A[i][j] - B[i][j];
      s += d * d;
    }
  return std::sqrt(s);
}

static double frob_norm(const double A[6][6]) {
  double s = 0.0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      s += A[i][j] * A[i][j];
  return std::sqrt(s);
}

static bool check_fd_tangent(
    const HSParams& p,
    const HSState& st_base,
    const double eps_target[6],
    ActiveRegime expected_reg,
    const char* label) {
  HSState st_alg;
  double C_alg[6][6];
  ActiveRegime reg = return_map_3d_comp(p, st_base, eps_target, st_alg, C_alg);
  if (reg != expected_reg) {
    std::fprintf(stderr, "FAIL [%s]: regime=%d != expected=%d\n",
                 label, static_cast<int>(reg), static_cast<int>(expected_reg));
    return false;
  }

  const double delta = 1.0e-7;
  double C_fd[6][6];
  for (int j = 0; j < 6; ++j) {
    double ep[6], em[6];
    for (int i = 0; i < 6; ++i) {
      ep[i] = eps_target[i];
      em[i] = eps_target[i];
    }
    ep[j] += delta;
    em[j] -= delta;
    HSState sp, sm;
    double Cp[6][6], Cm[6][6];
    return_map_3d_comp(p, st_base, ep, sp, Cp);
    return_map_3d_comp(p, st_base, em, sm, Cm);
    for (int i = 0; i < 6; ++i) {
      C_fd[i][j] = (sp.sig[i] - sm.sig[i]) / (2.0 * delta);
    }
  }

  double rel_err = frob_norm_diff(C_alg, C_fd) / std::max(1.0, frob_norm(C_fd));
  if (rel_err > 1.0e-7) {
    std::fprintf(stderr, "FAIL [%s]: FD rel_err=%.3e > 1e-7\n", label, rel_err);
    return false;
  }
  return true;
}

int main() {
  for (int lode : {LODE_MC, LODE_MN}) {
    HSParams p;
    p.E50_ref    = 37000.0;
    p.Eoed_ref   = 29600.0;
    p.Eur_ref    = 90000.0;
    p.m          = 0.50;
    p.c          = 5.0;
    p.phi_deg    = 38.0;
    p.psi_deg    = 10.0;
    p.K0nc       = 0.38;
    p.lode_shape = lode;
    if (!p.finalize()) return 1;

    HSState st_ocr, st_fail, st_iso, st_k0, st_fail_cap;
    const double sig_ocr[6]  = {150.0, 95.0, 65.0, 12.0, -8.0, 5.0};
    const double sig_fail[6] = {415.0, 95.0, 90.0, 12.0, -8.0, 5.0};
    const double sig_iso[6]  = {100.0, 100.0, 100.0, 0.0, 0.0, 0.0};
    const double sig_k0[6]   = {100.0, 38.0, 38.0, 0.0, 0.0, 0.0};
    initialize_state_from_stress(p, sig_ocr,  st_ocr,  15.0);
    initialize_state_from_stress(p, sig_fail, st_fail, 15.0);
    initialize_state_from_stress(p, sig_iso,  st_iso,  1.0);
    st_iso.gamma_p = 0.05;
    initialize_state_from_stress(p, sig_k0,   st_k0,   1.0);
    initialize_state_from_stress(p, sig_fail, st_fail_cap, 1.0);

    const double e_el[6]   = {-1.0e-4, 5.0e-5, 5.0e-5, -2.0e-5, 1.0e-5, -1.0e-5};
    const double e_sh[6]   = {1.2e-3, -4.0e-4, -5.0e-4, 3.0e-4, -1.0e-4, 1.0e-4};
    const double e_fl[6]   = {2.5e-2, -1.2e-2, -1.2e-2, 3.0e-3, -1.0e-3, 1.0e-3};
    const double e_cp[6]   = {2.0e-3, 1.9e-3, 1.8e-3, 1.0e-4, -5.0e-5, 5.0e-5};
    const double e_sc[6]   = {1.5e-3, 0.0, 0.0, 1.0e-4, -5.0e-5, 3.0e-5};
    const double e_fc[6]   = {2.5e-2, -9.0e-3, -9.0e-3, 3.0e-3, -1.0e-3, 1.0e-3};
    const double e_tn[6]   = {-2.0e-2, -2.0e-2, -2.0e-2, 0.0, 0.0, 0.0};

    if (!check_fd_tangent(p, st_ocr,      e_el, REGIME_ELASTIC,     "elastic"))     return 2;
    if (!check_fd_tangent(p, st_ocr,      e_sh, REGIME_SHEAR,       "shear"))       return 3;
    if (!check_fd_tangent(p, st_fail,     e_fl, REGIME_FAILURE,     "failure"))     return 4;
    if (!check_fd_tangent(p, st_iso,      e_cp, REGIME_CAP,         "cap"))         return 5;
    if (!check_fd_tangent(p, st_k0,       e_sc, REGIME_SHEAR_CAP,   "shear_cap"))   return 6;
    if (!check_fd_tangent(p, st_fail_cap, e_fc, REGIME_FAILURE_CAP, "failure_cap")) return 10;
    if (!check_fd_tangent(p, st_ocr,      e_tn, REGIME_TENSION,     "tension"))     return 7;

    // 3D cap normality check at general Lode angle (theta != +/- pi/6):
    {
      HSState st_cap;
      double C_cap[6][6];
      ActiveRegime reg_cap = return_map_3d_comp(p, st_iso, e_cp, st_cap, C_cap);
      if (reg_cap != REGIME_CAP) return 13;
      double p_n = 0.0, q_n = 0.0, x_n = 1.0, Th_n = 1.0, dTh_dx_n = 0.0;
      double k3_n = 1.0 / 3.0, dk3_dx_n = 0.0, Th_sig_n = 1.0, dTh_sig_dx_n = 0.0;
      double dq_n[6], dTh_n[6], ds_dp = 0.0, ds_dq = 0.0, s_n = 1.0;
      invariants_and_lode_grad(p, st_iso.sig, p_n, q_n, x_n, Th_n, dq_n, dTh_n);
      lode_factors_and_derivs(p, x_n, Th_n, dTh_dx_n, k3_n, dk3_dx_n, Th_sig_n, dTh_sig_dx_n);
      stiffness_scale_and_derivs(p, p_n, q_n, s_n, ds_dp, ds_dq, k3_n);
      double Ke = p.Kur_ref * s_n;
      double Ge = p.Gur_ref * s_n;
      double dsig[6];
      for (int i = 0; i < 6; ++i) dsig[i] = st_cap.sig[i] - st_iso.sig[i];
      double dp = (dsig[0] + dsig[1] + dsig[2]) / 3.0;
      double deps_e[6] = {
          dp / (3.0 * Ke) + (dsig[0] - dp) / (2.0 * Ge),
          dp / (3.0 * Ke) + (dsig[1] - dp) / (2.0 * Ge),
          dp / (3.0 * Ke) + (dsig[2] - dp) / (2.0 * Ge),
          dsig[3] / Ge,
          dsig[4] / Ge,
          dsig[5] / Ge
      };
      double deps_p[6];
      for (int i = 0; i < 6; ++i) deps_p[i] = e_cp[i] - deps_e[i];
      double p_c = 0.0, q_c = 0.0, x_c = 1.0, Th_c = 1.0;
      double dq_c[6], dTh_c[6];
      invariants_and_lode_grad(p, st_cap.sig, p_c, q_c, x_c, Th_c, dq_c, dTh_c);
      double p_hat_c = p_c + p.a_apex;
      double grad_fc[6];
      for (int i = 0; i < 6; ++i) {
        double m_i = (i < 3) ? (1.0 / 3.0) : 0.0;
        grad_fc[i] = 2.0 * p_hat_c * m_i + (2.0 * q_c / (p.alpha_cap * p.alpha_cap)) * dq_c[i];
      }
      static const double w_t[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};
      double dot = 0.0, np2 = 0.0, ng2 = 0.0;
      for (int i = 0; i < 6; ++i) {
        dot += w_t[i] * deps_p[i] * grad_fc[i];
        np2 += w_t[i] * deps_p[i] * deps_p[i];
        ng2 += w_t[i] * grad_fc[i] * grad_fc[i];
      }
      double cos_ang = dot / std::sqrt(np2 * ng2);
      if (std::abs(1.0 - cos_ang) > 1.0e-13) {
        std::fprintf(stderr, "FAIL [cap_normality]: 1 - cos_ang = %.3e\n", 1.0 - cos_ang);
        return 14;
      }
    }

    // HSsmall elastoplastic tangent check:
    HSParams p_hss = p;
    p_hss.E50_ref          = 30000.0;
    p_hss.Eoed_ref         = 30000.0;
    p_hss.Eur_ref          = 90000.0;
    p_hss.m                = 0.55;
    p_hss.c                = 0.0;
    p_hss.phi_deg          = 38.0;
    p_hss.psi_deg          = 8.0;
    p_hss.use_small_strain = true;
    p_hss.G0_ref           = 112500.0;
    p_hss.gamma_07         = 2.0e-4;
    p_hss.alpha_cap        = -1.0;
    p_hss.beta_cap         = -1.0;
    if (!p_hss.finalize()) return 8;

    HSState st_hss0, st_hss1, st_hss_unl, st_hss_iso;
    double C_tmp[6][6];
    initialize_state_from_stress(p_hss, sig_ocr, st_hss0, 5.0);
    const double e_pre[6]  = {2.0e-4, -1.0e-4, -1.0e-4, 1.5e-4, 0.0, 0.0};
    const double e_next[6] = {1.4e-3, -6.0e-4, -5.0e-4, 4.0e-4, -1.0e-4, 1.0e-4};
    return_map_3d_comp(p_hss, st_hss0, e_pre, st_hss1, C_tmp);
    if (!check_fd_tangent(p_hss, st_hss1, e_next, REGIME_SHEAR, "hss_shear")) return 9;

    const double e_unl[6]   = {-2.0e-4, 1.0e-4, 1.0e-4, -1.5e-4, 0.0, 0.0};
    const double e_cross[6] = {3.2e-4, -1.4e-4, -1.5e-4, 2.5e-4, -2.0e-5, 1.0e-5};
    return_map_3d_comp(p_hss, st_hss1, e_unl, st_hss_unl, C_tmp);
    if (!check_fd_tangent(p_hss, st_hss_unl, e_cross, REGIME_SHEAR, "hss_masing_cross")) return 11;

    initialize_state_from_stress(p_hss, sig_iso, st_hss_iso, 1.0);
    st_hss_iso.gamma_p = 0.05;
    if (!check_fd_tangent(p_hss, st_hss_iso, e_cp, REGIME_CAP, "hss_cap")) return 12;

    // Verify 1D fallback for failure_mode == false (REGIME_SHEAR_CAP) in solve_regime4_corner:
    {
      double Th_f = 1.0, dTh_f = 0.0, k3_f = 1.0 / 3.0, dk3_f = 0.0, Ths_f = 1.0, dThs_f = 0.0;
      lode_factors_and_derivs(p, 0.35, Th_f, dTh_f, k3_f, dk3_f, Ths_f, dThs_f);
      double p_2d = 0.0, q_2d = 0.0, dg_2d = 0.0, pp_2d = 0.0, sens_2d[2][4];
      double p_1d = 0.0, q_1d = 0.0, dg_1d = 0.0, pp_1d = 0.0, sens_1d[2][4];
      bool ok2d = solve_regime4_corner(
          p, 110.0, 138.5, Th_f, dTh_f, k3_f, dk3_f, Ths_f, dThs_f,
          p.Kur_ref, p.Gur_ref, st_k0.gamma_p, st_k0.p_p, st_k0.void_ratio,
          false, p_2d, q_2d, dg_2d, pp_2d, sens_2d, false);
      bool ok1d = solve_regime4_corner(
          p, 110.0, 138.5, Th_f, dTh_f, k3_f, dk3_f, Ths_f, dThs_f,
          p.Kur_ref, p.Gur_ref, st_k0.gamma_p, st_k0.p_p, st_k0.void_ratio,
          false, p_1d, q_1d, dg_1d, pp_1d, sens_1d, true);
      if (!ok2d || !ok1d || std::abs(p_2d - p_1d) > 1.0e-10 || std::abs(q_2d - q_1d) > 1.0e-10) {
        std::fprintf(stderr, "FAIL [regime4_1d_fallback]: p_2d=%.12f p_1d=%.12f q_2d=%.12f q_1d=%.12f\n",
                     p_2d, p_1d, q_2d, q_1d);
        return 15;
      }

      // High-confinement / small-p_ref ULP scaling check for solve_regime4_corner:
      HSParams p_hi = p;
      p_hi.p_ref     = 0.1;
      p_hi.alpha_cap = -1.0;
      p_hi.beta_cap  = -1.0;
      if (!p_hi.finalize()) return 16;
      double p_hi2d = 0.0, q_hi2d = 0.0, dg_hi2d = 0.0, pp_hi2d = 0.0, sens_hi2d[2][4];
      double p_hi1d = 0.0, q_hi1d = 0.0, dg_hi1d = 0.0, pp_hi1d = 0.0, sens_hi1d[2][4];
      bool ok_hi2d = solve_regime4_corner(
          p_hi, 65000.0, 30000.0, Th_f, dTh_f, k3_f, dk3_f, Ths_f, dThs_f,
          p_hi.Kur_ref, p_hi.Gur_ref, 1.0e-5, 50000.0, 0.65,
          false, p_hi2d, q_hi2d, dg_hi2d, pp_hi2d, sens_hi2d, false);
      bool ok_hi1d = solve_regime4_corner(
          p_hi, 65000.0, 30000.0, Th_f, dTh_f, k3_f, dk3_f, Ths_f, dThs_f,
          p_hi.Kur_ref, p_hi.Gur_ref, 1.0e-5, 50000.0, 0.65,
          false, p_hi1d, q_hi1d, dg_hi1d, pp_hi1d, sens_hi1d, true);
      if (!ok_hi2d || !ok_hi1d ||
          std::abs(p_hi2d - p_hi1d) / p_hi2d > 1.0e-11 ||
          std::abs(q_hi2d - q_hi1d) / q_hi2d > 1.0e-11) {
        std::fprintf(stderr, "FAIL [regime4_hi_ulp]: ok2d=%d ok1d=%d p_2d=%.12f p_1d=%.12f\n",
                     static_cast<int>(ok_hi2d), static_cast<int>(ok_hi1d), p_hi2d, p_hi1d);
        return 17;
      }
    }

    // Emit machine-precision comparison row for Python parity test:
    HSState st_out;
    double C_out[6][6];
    ActiveRegime reg_out = return_map_3d_comp(p, st_k0, e_sc, st_out, C_out);
    std::cout << std::setprecision(17)
              << "PARITY " << lode << " " << static_cast<int>(reg_out) << " "
              << p.alpha_cap << " " << p.beta_cap << " "
              << st_out.sig[0] << " " << st_out.sig[1] << " " << st_out.sig[2] << " "
              << st_out.sig[3] << " " << st_out.sig[4] << " " << st_out.sig[5] << " "
              << st_out.gamma_p << " " << st_out.p_p << " "
              << C_out[0][0] << " " << C_out[0][1] << " " << C_out[3][3] << "\n";
  }
  std::cout << "LADRUNO_HS_KERNEL_OK\n";
  return 0;
}
