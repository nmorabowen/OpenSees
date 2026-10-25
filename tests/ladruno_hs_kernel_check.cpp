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

    HSState st_hss0, st_hss1;
    double C_tmp[6][6];
    initialize_state_from_stress(p_hss, sig_ocr, st_hss0, 5.0);
    const double e_pre[6]  = {2.0e-4, -1.0e-4, -1.0e-4, 1.5e-4, 0.0, 0.0};
    const double e_next[6] = {1.4e-3, -6.0e-4, -5.0e-4, 4.0e-4, -1.0e-4, 1.0e-4};
    return_map_3d_comp(p_hss, st_hss0, e_pre, st_hss1, C_tmp);
    if (!check_fd_tangent(p_hss, st_hss1, e_next, REGIME_SHEAR, "hss_shear")) return 9;

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
