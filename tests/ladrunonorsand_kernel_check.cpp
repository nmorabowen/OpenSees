// Standalone driver for the LadrunoNORSAND return-map kernel (WP-144 P1a).
//
// Includes ONLY SRC/material/nD/LadrunoNorSandKernel.h (pure, OpenSees-free), runs a
// set of prescribed strain PATHS threading the committed state, and prints per step the
// stress, the state, the StepInfo and the 6x6 tangent in a parseable format:
//
//   SCENARIO <name>
//   INIT <eps_e x6> <pi_i> <v> <v0>
//   STEP <k> INFO <refusal> <finest> <finest_sub> <plastic> <vertex> <cap_active> <local_iters> <pi_iters> <substeps>
//   STEP <k> SIGMA <6 tensor comps {00,11,22,01,12,02}>
//   STEP <k> STATE <eps_e x6> <pi_i> <v> <v0> <eps_p_v> <eps_p_s> <D_last>
//   STEP <k> TAN <36 entries, row-major, C[I][J] = d sigma_I / d eps_J (tensor shear)>
//   CHECK <name> <value> <tol> PASS|FAIL
//
// It then runs self-checks (CHECK lines): FD of the tangent (pins the shear convention
// of C), the symmetric elastic / non-symmetric plastic engineering tangent, the
// validate() refusals, the frozen state on a refused step, bounded work on wild trial
// increments, and an elastic closed loop. Exit status 1 if any CHECK fails.
// Parity against the O2 oracle is run by
//   Ladruno_files/testbed/norsand_oracle/kernel_parity/test_kernel_parity.py
//
// Build (one line):
//   g++ -std=c++17 -O2 -Wall -Wextra -Werror -I SRC/material/nD tests/ladrunonorsand_kernel_check.cpp -o nsk
//   ./nsk

#include <LadrunoNorSandKernel.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace ladruno_norsand;

static int g_fail = 0;

static void check(const char* name, double value, double tol, bool pass)
{
  printf("CHECK %s %.6e %.6e %s\n", name, value, tol, pass ? "PASS" : "FAIL");
  if (!pass) g_fail++;
}

// ---- parameter sets (sheet §14 K2 / G1 common_params) ------------------------ //
static Params paper()
{
  Params P;
  P.p0 = -100.0; P.kappa_hat = 0.01; P.eps_v0 = 0.0; P.mu0 = 5400.0; P.alpha0 = 0.0;
  P.M = 1.2; P.N = 0.4; P.N_bar = 0.2; P.rho = 0.7; P.rho_bar = 0.8; P.chi = -3.5; P.h = 280.0;
  P.lambda_tilde = 0.0135; P.v_c0 = 1.81; P.e0 = 0.83; P.lambda_c = 0.027; P.xi = 0.45;
  P.p_a = 101.325; P.c1 = 0.05; P.c2 = 0.15; P.csl_mode = 0; P.zeta = 0; P.cap = 0;
  return P;
}

static Params fork()
{
  Params P = paper();
  P.csl_mode = 1; P.M = 1.3309; P.N = 0.3; P.N_bar = 0.2; P.rho = 0.71; P.rho_bar = 0.75;
  return P;
}

static double v0_for(const Params& P, double pi0, double psi0)
{
  if (P.csl_mode == 0) return psi0 + P.v_c0 - P.lambda_tilde * std::log(-pi0);
  return 1.0 + P.e0 + psi0 - P.lambda_c * std::pow(-pi0 / P.p_a, P.xi);
}

static void print6(const char* tag, int k, const double* x, int n)
{
  printf("STEP %d %s", k, tag);
  for (int i = 0; i < n; ++i) printf(" %.17g", x[i]);
  printf("\n");
}

static State init_state(const Params& P, double p, double pi0, double psi0, double* v0_out = nullptr)
{
  const double sig0[6] = {p, p, p, 0.0, 0.0, 0.0};
  const double v0 = std::isnan(pi0) ? 1.7 : v0_for(P, pi0, psi0);
  State s;
  std::string msg;
  const int rc = initialState(P, sig0, v0, pi0, s, msg);
  if (rc) { printf("INIT FAILED %d %s\n", rc, msg.c_str()); g_fail++; }
  if (v0_out) *v0_out = v0;
  return s;
}

// Runs a prescribed path; prints everything; stops after the first refusal.
static State run(const char* name, const Params& P, State s, const std::vector<std::vector<double>>& path)
{
  printf("SCENARIO %s\n", name);
  printf("INIT");
  for (int i = 0; i < 6; ++i) printf(" %.17g", s.eps_e[i]);
  printf(" %.17g %.17g %.17g\n", s.pi_i, s.v, s.v0);
  for (size_t k = 0; k < path.size(); ++k) {
    State np1;
    double sig[6], C[6][6];
    StepInfo info;
    int finest = 0, sub = 0;
    detail::step_ex(P, s, path[k].data(), np1, sig, C, info, finest, sub);
    printf("STEP %d INFO %d %d %d %d %d %d %d %d %d\n", (int)k, info.refusal, finest, sub, info.plastic,
           info.vertex, info.cap_active, info.local_iters, info.pi_iters, info.substeps);
    print6("SIGMA", (int)k, sig, 6);
    const double st[12] = {np1.eps_e[0], np1.eps_e[1], np1.eps_e[2], np1.eps_e[3], np1.eps_e[4], np1.eps_e[5],
                           np1.pi_i, np1.v, np1.v0, np1.eps_p_v, np1.eps_p_s, np1.D_last};
    print6("STATE", (int)k, st, 12);
    print6("TAN", (int)k, &C[0][0], 36);
    s = np1;
    if (info.refusal != OK) break;
  }
  return s;
}

static std::vector<std::vector<double>> repeat(const double d[6], int n)
{
  return std::vector<std::vector<double>>(n, std::vector<double>(d, d + 6));
}

// G1 NONCOAXIAL path: E1 = diag(0.003, 0.003, -0.010) over t in [0, 0.6], then E2 (with
// tensor shear 01 = 0.006, 12 = 0.004) over [0.6, 1]; n equal steps (G1 noncoax_deps).
static std::vector<std::vector<double>> noncoax(int n)
{
  const double E1[6] = {0.003, 0.003, -0.010, 0.0, 0.0, 0.0};
  const double E2[6] = {0.0005, 0.0005, -0.002, 0.006, 0.004, 0.0};
  const double KNOT = 0.6;
  auto eps = [&](double s, double* out) {
    for (int i = 0; i < 6; ++i)
      out[i] = E1[i] * std::fmin(s, KNOT) / KNOT + E2[i] * std::fmax(s - KNOT, 0.0) / (1.0 - KNOT);
  };
  std::vector<std::vector<double>> path;
  for (int k = 0; k < n; ++k) {
    double a[6], b[6];
    eps(double(k) / n, a);
    eps(double(k + 1) / n, b);
    std::vector<double> d(6);
    for (int i = 0; i < 6; ++i) d[i] = b[i] - a[i];
    path.push_back(d);
  }
  return path;
}

static double maxabs(const double* x, int n)
{
  double m = 0.0;
  for (int i = 0; i < n; ++i) m = std::fmax(m, std::fabs(x[i]));
  return m;
}

int main()
{
  const double nan = std::nan("");

  // ---- prescribed paths --------------------------------------------------------
  {   // undrained TXC (paper), loose start
    const Params P = paper();
    const int n = 25; const double da = -0.05 / n;
    const double d[6] = {-0.5 * da, -0.5 * da, da, 0, 0, 0};
    run("txc_undrained_paper", P, init_state(P, -100.0, -105.0, 0.01), repeat(d, n));
  }
  {   // non-coaxial (fork CSL)
    const Params P = fork();
    run("noncoaxial_fork", P, init_state(P, -100.0, -105.0, 0.01), noncoax(25));
  }
  {   // non-coaxial (paper, Gudehus-Argyris)
    Params P = paper(); P.zeta = 1; P.rho = 0.8; P.rho_bar = 0.9;
    run("noncoaxial_paper_GA", P, init_state(P, -100.0, -105.0, 0.01), noncoax(25));
  }
  {   // smooth cap, AMP_STOP near-isotropic path, n = 40 (substepping is normal here)
    Params P = paper(); P.cap = 2; P.c1 = 0.05; P.c2 = 0.15;
    const int n = 40;
    const double d[6] = {(-0.01 + 2e-3) / n, -0.01 / n, (-0.01 - 2e-3) / n, 0, 0, 0};
    run("cap_smooth_ampstop", P, init_state(P, -100.0, -80.0, -0.05), repeat(d, n));
  }
  {   // apex: on the surface at p = -100 (hydrostatic), three hydrostatic steps (vertex rule)
    const Params P = paper();
    const double d[6] = {-1e-3, -1e-3, -1e-3, 0, 0, 0};
    run("apex_hydrostatic", P, init_state(P, -100.0, nan, 0.0), repeat(d, 3));
  }

  // ---- CHECK: tangent vs central FD (non-coaxial plastic step, off the corners) ----
  {
    const Params P = paper();
    State s = init_state(P, -100.0, -80.0, -0.05);
    const double EPRE[6] = {4.0e-4, -1.0e-3, 0.0, 2.0e-4, 0.0, 1.0e-4};     // G1 E_PRE
    for (int k = 0; k < 10; ++k) {
      State np1; double sig[6], C[6][6]; StepInfo info;
      step(P, s, EPRE, np1, sig, C, info);
      s = np1;
    }
    const double DFD[6] = {1.0e-4, -6.0e-4, 2.0e-4, 1.0e-4, 0.0, 0.5e-4};    // G1 DEPS_FD
    State np1; double sig[6], C[6][6]; StepInfo info;
    step(P, s, DFD, np1, sig, C, info);
    const bool base_ok = info.refusal == OK && info.plastic && info.substeps == 1;
    const double h = 1e-7;
    double err = 0.0;
    bool all_ok = base_ok;
    for (int J = 0; J < 6; ++J) {
      double dp[6], dm[6];
      for (int i = 0; i < 6; ++i) { dp[i] = DFD[i]; dm[i] = DFD[i]; }
      dp[J] += h; dm[J] -= h;
      State a, b; double sp[6], sm[6], Ca[6][6], Cb[6][6]; StepInfo ia, ib;
      step(P, s, dp, a, sp, Ca, ia);
      step(P, s, dm, b, sm, Cb, ib);
      all_ok = all_ok && ia.refusal == OK && ib.refusal == OK && ia.substeps == 1 && ib.substeps == 1;
      for (int I = 0; I < 6; ++I) err = std::fmax(err, std::fabs((sp[I] - sm[I]) / (2.0 * h) - C[I][J]));
    }
    err /= maxabs(&C[0][0], 36);
    check("fd_tangent_noncoaxial_plastic", err, 1e-6, all_ok && err <= 1e-6);

    // engineering form D (shear columns halved): non-symmetric when plastic
    double D[6][6], asym = 0.0;
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) D[I][J] = (J < 3) ? C[I][J] : 0.5 * C[I][J];
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) asym = std::fmax(asym, std::fabs(D[I][J] - D[J][I]));
    asym /= maxabs(&D[0][0], 36);
    check("plastic_tangent_is_nonsymmetric", asym, 1e-4, asym > 1e-4);
    // elastic tangent at the same state: symmetric in engineering form
    double Ce[6][6], De[6][6], easym = 0.0;
    elasticTangent(P, np1, Ce);
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) De[I][J] = (J < 3) ? Ce[I][J] : 0.5 * Ce[I][J];
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) easym = std::fmax(easym, std::fabs(De[I][J] - De[J][I]));
    easym /= maxabs(&De[0][0], 36);
    check("elastic_tangent_engineering_symmetric", easym, 1e-12, easym <= 1e-12);
  }

  // ---- CHECK: validate() refusals (owner decisions in force) ----------------------
  {
    std::string msg; bool warn = false;
    Params P = paper(); P.N_bar = P.N;                       // the counterexample: rho < rho_bar but beta = 1
    int rc = validate(P, msg, warn);
    check("refuse_counterexample_Nbar_eq_N", rc, 13, rc == 13);
    P = paper(); P.N_bar = 0.5;
    rc = validate(P, msg, warn);
    check("refuse_Nbar_gt_N", rc, 12, rc == 12);
    P = paper(); P.rho = 0.5;
    rc = validate(P, msg, warn);
    check("refuse_WW_rho_half", rc, 11, rc == 11);
    P = paper(); P.rho_bar = 0.5; P.rho = 0.5;
    rc = validate(P, msg, warn);
    check("refuse_WW_rhobar_half", rc, 11, rc == 11);
    P = paper(); P.zeta = 1; P.rho = 0.75; P.rho_bar = 0.9;
    rc = validate(P, msg, warn);
    check("refuse_GA_rho_below_7_9", rc, 11, rc == 11);
    P = paper(); P.rho = 0.9; P.rho_bar = 0.8;
    rc = validate(P, msg, warn);
    check("warn_only_rho_gt_rhobar", rc, 0, rc == 0 && warn);
    P = paper();
    rc = validate(P, msg, warn);
    check("accept_K2_set", rc, 0, rc == 0 && !warn);
  }

  // ---- CHECK: refusal freezes the state (no-cap AMP_STOP path; O2 refuses step 13) ----
  {
    const Params P = paper();
    State s = init_state(P, -100.0, -80.0, -0.05);
    const int n = 40;
    const double d[6] = {(-0.01 + 2e-3) / n, -0.01 / n, (-0.01 - 2e-3) / n, 0, 0, 0};
    int refused_at = -1;
    bool frozen = false;
    auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < n; ++k) {
      State np1; double sig[6], C[6][6], sn[6]; StepInfo info;
      step(P, s, d, np1, sig, C, info);
      if (info.refusal != OK) {
        refused_at = k;
        stress(P, s, sn);
        frozen = info.refusal == SUBSTEPS_EXHAUSTED && info.substeps == 256
              && np1.pi_i == s.pi_i && np1.v == s.v && np1.v0 == s.v0 && np1.eps_p_v == s.eps_p_v
              && np1.eps_p_s == s.eps_p_s && np1.D_last == s.D_last;
        for (int i = 0; i < 6; ++i) frozen = frozen && np1.eps_e[i] == s.eps_e[i] && sig[i] == sn[i];
        break;
      }
      s = np1;
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("INFO nocap_ampstop refused_at_step %d (zero-based) in %.3f s\n", refused_at, secs);
    check("refusal_freezes_state", refused_at, 0, refused_at >= 0 && frozen);
  }

  // ---- CHECK: bounded work on wild trial increments (material checklist) ----------
  {
    const Params P = paper();
    const State s = init_state(P, -100.0, -80.0, -0.05);
    const double wild[4][6] = {{1e4, -1e4, 5e3, 2e3, -3e3, 1e3},     // |deps| ~ 1e4: elastic overflow
                               {-1.0, -1.0, -1.0, 0.3, 0.0, 0.0},      // huge finite compression
                               {1.0, 1.0, 1.0, 0.0, 0.0, 0.0},         // tension: p > 0 at the trial
                               {0.2, -0.3, 0.1, 0.25, -0.1, 0.05}};    // large deviatoric
    for (int w = 0; w < 4; ++w) {
      State np1; double sig[6], C[6][6]; StepInfo info;
      auto t0 = std::chrono::steady_clock::now();
      step(P, s, wild[w], np1, sig, C, info);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      char nm[64];
      std::snprintf(nm, sizeof(nm), "wild_increment_%d_bounded_seconds", w);
      printf("INFO wild %d refusal %d substeps %d plastic %d time %.3f s\n", w, info.refusal, info.substeps,
             info.plastic, secs);
      // bound: refuse (or accept a finite answer) within 10 s; never NaN on accept
      bool fin = true;
      for (int i = 0; i < 6; ++i) fin = fin && std::isfinite(sig[i]);
      check(nm, secs, 10.0, secs <= 10.0 && fin);
    }
    State np1; double sig[6], C[6][6]; StepInfo info;
    const double dn[6] = {nan, 0, 0, 0, 0, 0};
    step(P, s, dn, np1, sig, C, info);
    check("nan_increment_refused", info.refusal, 0, info.refusal != OK);
  }

  // ---- CHECK: elastic closed loop returns the stress, D = 0 -----------------------
  {
    const Params P = paper();
    const State s0 = init_state(P, -100.0, -80.0, -0.05);
    double sig0[6];
    stress(P, s0, sig0);
    const double d[6] = {1e-4, -5e-5, -5e-5, 2e-5, 0.0, -1e-5};
    const double dm[6] = {-1e-4, 5e-5, 5e-5, -2e-5, 0.0, 1e-5};
    State s1, s2; double sg1[6], sg2[6], C[6][6]; StepInfo i1, i2;
    step(P, s0, d, s1, sg1, C, i1);
    step(P, s1, dm, s2, sg2, C, i2);
    double e = 0.0;
    for (int i = 0; i < 6; ++i) e = std::fmax(e, std::fabs(sg2[i] - sig0[i]));
    e /= maxabs(sig0, 6);
    check("elastic_closed_loop", e, 1e-12,
          e <= 1e-12 && !i1.plastic && !i2.plastic && s2.D_last == 0.0 && s2.pi_i == s0.pi_i);
  }

  printf("SUMMARY %s (%d failed)\n", g_fail ? "FAIL" : "PASS", g_fail);
  return g_fail ? 1 : 0;
}
