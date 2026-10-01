// WP-144 P1b: confirms the Voigt <-> tensor mapping used by LadrunoNorSand.cpp against the REAL kernel.
//   strain  eps_tensor = e (normal), e/2 (shear);  stress voigt = tensor;  T[a][b] = C[a][b] * (b>=3 ? 0.5 : 1)
// by central finite differences of the kernel stress with respect to the ENGINEERING strain vector,
// at a plastic, non-coaxial state. Build/run on Esmeralda (login node, seconds):
//   g++ -std=c++17 -O2 -I SRC/material/nD -o /tmp/voigt_check Ladruno_files/testbed/norsand_oracle/kernel_parity/tools/voigt_check.cpp && /tmp/voigt_check
#include "LadrunoNorSandKernel.h"
#include <cstdio>
#include <cmath>
using namespace ladruno_norsand;

static StepInfo g_info{};
static int stepEng(const Params& p, const State& n, const double e[6], State& np1, double sig[6], double C[6][6]) {
  double d[6]; for (int i = 0; i < 6; i++) d[i] = (i >= 3) ? 0.5 * e[i] : e[i];
  int rc = step(p, n, d, np1, sig, C, g_info);
  return rc != 0 ? rc : g_info.refusal;
}

int main() {
  Params p{};
  p.p0 = -100; p.kappa_hat = 0.01; p.eps_v0 = 0; p.mu0 = 5400; p.alpha0 = 0; p.M = 1.2; p.N = 0.4; p.N_bar = 0.2;
  p.rho = 0.7; p.rho_bar = 0.8; p.chi = -3.5; p.h = 280; p.lambda_tilde = 0.0135; p.v_c0 = 1.81; p.e0 = 0.83;
  p.lambda_c = 0.027; p.xi = 0.45; p.p_a = 101.325; p.c1 = 0.05; p.c2 = 0.15; p.csl_mode = 0; p.zeta = 0; p.cap = 0;
  std::string msg; bool w = false;
  if (validate(p, msg, w)) { std::printf("validate refused: %s\n", msg.c_str()); return 2; }
  double sig0[6] = {-100, -100, -100, 0, 0, 0};
  State s0{}; if (initialState(p, sig0, 1.70, std::nan(""), s0, msg)) { std::printf("init refused: %s\n", msg.c_str()); return 2; }
  // drive to a plastic non-coaxial state with engineering-shear strains
  State a = s0, b{}; double sig[6], C[6][6];
  const double path[3][6] = {{-0.004, 0.0015, 0.0008, 0.0030, 0.0, 0.0012}, {-0.003, 0.001, 0.0005, 0.0025, 0.0010, 0.0}, {-0.002, 0.0, 0.0004, 0.0020, 0.0005, 0.0008}};
  for (int k = 0; k < 3; k++) { if (stepEng(p, a, path[k], b, sig, C)) { std::printf("path refused at %d\n", k); return 3; } a = b; }
  // base increment + tangent at the plastic state a
  const double e0[6] = {-0.0010, 0.0004, 0.0002, 0.0012, 0.0006, 0.0003};
  State r0{}; if (stepEng(p, a, e0, r0, sig, C)) { std::printf("base step refused\n"); return 3; }
  const StepInfo info = g_info;
  double T[6][6]; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) T[i][j] = C[i][j] * (j >= 3 ? 0.5 : 1.0);
  double maxrel = 0, maxabs = 0, scale = 0; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) scale = std::fmax(scale, std::fabs(T[i][j]));
  // symmetric-shear-consistency: vary the ENGINEERING component j, compare stress change to column j of T
  const double h = 1e-7;
  for (int j = 0; j < 6; j++) {
    double ep[6], em[6]; for (int i = 0; i < 6; i++) { ep[i] = e0[i]; em[i] = e0[i]; } ep[j] += h; em[j] -= h;
    State sp{}, sm{}; double sgp[6], sgm[6], Cd[6][6];
    if (stepEng(p, a, ep, sp, sgp, Cd) || stepEng(p, a, em, sm, sgm, Cd)) { std::printf("FD step refused\n"); return 3; }
    for (int i = 0; i < 6; i++) {
      double fd = (sgp[i] - sgm[i]) / (2 * h), err = std::fabs(fd - T[i][j]);
      maxabs = std::fmax(maxabs, err); maxrel = std::fmax(maxrel, err / scale);
    }
  }
  // the NON-symmetry the contract promises (never a symmetric solver)
  double asym = 0; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) asym = std::fmax(asym, std::fabs(T[i][j] - T[j][i]));
  std::printf("plastic=%d substeps=%d  |T_fd - T| max abs %.3e, rel-to-max|T| %.3e  (max|T| %.3e), max|T-T^T| %.3e\n",
              info.plastic, info.substeps, maxabs, maxrel, scale, asym);
  // v0 != v: the state a carries v0=1.70 while v has moved
  std::printf("state a: v=%.6f v0=%.6f (v0 != v: %s)\n", a.v, a.v0, (a.v != a.v0) ? "yes" : "NO");
  return (maxrel < 1e-5) ? 0 : 1;
}
