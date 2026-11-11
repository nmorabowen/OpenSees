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
// increments, the exponential specific-volume law v = v0 exp(tr eps), an elastic closed loop, and the chained tangent of substepped increments
// (sheet §9.6) against the central FD of the whole increment (ladder, uniform and
// recursive-halving fractions, m = 1 reduction, vertex), and (round 4, sheet §9.7 (S.57)-(S.59)) the image-pressure
// rule at a post floor: K1.12 (C), the §9.7 wet-floor table, K1.14b, K1.16 (the H19 states), the (S.58) census and
// the forced (1/2, 1/2) gate of the chained row (S.59). Exit status 1 if any CHECK fails.
// Parity against the O2 oracle is run by
//   Ladruno_files/testbed/norsand_oracle/kernel_parity/test_kernel_parity.py
//
// Build (one line):
//   g++ -std=c++17 -O2 -Wall -Wextra -Werror -Wshadow -I SRC/material/nD tests/ladrunonorsand_kernel_check.cpp -o nsk
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

// Chained tangent (S.47) of the increment d taken with the fixed fractions fr, against the central FD
// of the WHOLE increment (same fractions at every FD point; O2 selfcheck chain_vs_fd). e: max over the
// six columns of ||C_J - FD_J|| / ||C_J||; el: the same for the last sub-increment's (S.33) CTO.
// Returns false if any of the 13 evaluations refused or the branch (plastic flag) changed.
static bool chain_fd(const Params& P, const State& s, const double d[6], const std::vector<double>& fr, double h,
                     double& e, double& el)
{
  const int m = (int)fr.size();
  State np1; double sig[6], C[6][6], Cl[6][6]; StepInfo info, il;
  detail::step_fractions(P, s, d, fr.data(), m, true, np1, sig, C, info);
  detail::step_fractions(P, s, d, fr.data(), m, false, np1, sig, Cl, il);
  bool ok = info.refusal == OK && il.refusal == OK && info.plastic;
  e = 0.0; el = 0.0;
  for (int J = 0; J < 6; ++J) {
    double dp[6], dm[6];
    for (int i = 0; i < 6; ++i) { dp[i] = d[i]; dm[i] = d[i]; }
    dp[J] += h; dm[J] -= h;
    State a, b; double sp[6], sm[6], Ca[6][6], Cb[6][6]; StepInfo ia, ib;
    detail::step_fractions(P, s, dp, fr.data(), m, true, a, sp, Ca, ia);
    detail::step_fractions(P, s, dm, fr.data(), m, true, b, sm, Cb, ib);
    ok = ok && ia.refusal == OK && ib.refusal == OK && ia.plastic == info.plastic && ib.plastic == info.plastic;
    double n2 = 0.0, d2 = 0.0, nl2 = 0.0, dl2 = 0.0;
    for (int I = 0; I < 6; ++I) {
      const double fd = (sp[I] - sm[I]) / (2.0 * h);
      // norm of the 3x3 column tensor: shear rows count twice
      const double wgt = I < 3 ? 1.0 : 2.0;
      n2 += wgt * (C[I][J] - fd) * (C[I][J] - fd);
      d2 += wgt * C[I][J] * C[I][J];
      nl2 += wgt * (Cl[I][J] - fd) * (Cl[I][J] - fd);
      dl2 += wgt * Cl[I][J] * Cl[I][J];
    }
    e = std::fmax(e, std::sqrt(n2 / d2));
    el = std::fmax(el, std::sqrt(nl2 / dl2));
  }
  return ok;
}

// =====================================================================================================
// Round 4 (sheet 144a §9.1 step 5g, §9.7 (S.57)-(S.59), §13.12 (C), §13.14b, §13.16; owner decision O13 A):
// the image-pressure rule at a post floor. Every expected number below is the SHEET's (round 4b), written
// in before the kernel was run; tolerances are stated beside them.
// =====================================================================================================
static double fdiff(double a, double b) { return std::fabs(a - b); }

// principal stresses of a tensor sigma (eig_sym3) and F(sigma, pi) = zeta q + p eta(p, pi) (vertex: p eta)
static bool Fval(const Params& P, const double sig[6], double pi, double& F)
{
  double S[3][3], w[3], V[3][3];
  detail::t6_to_m(sig, S);
  if (!detail::eig_sym3(S, w, V)) return false;
  detail::Invariants inv; detail::invariants(w, inv);
  detail::Flow fl;
  if (detail::flow(P, inv, pi, fl) != detail::EE_NONE) return false;
  F = fl.F;
  return true;
}

// pi_f(sigma): the image pressure of the surface through sigma (inverse of (S.12) at (p, zeta q/|p|))
static bool pif_of(const Params& P, const double sig[6], double& pif, double& eta)
{
  double S[3][3], w[3], V[3][3];
  detail::t6_to_m(sig, S);
  if (!detail::eig_sym3(S, w, V)) return false;
  detail::Invariants inv; detail::invariants(w, inv);
  eta = 0.0;
  if (!inv.vertex) { double z, zy, zyy; detail::zeta_y(inv.theta, P.rho, P.zeta, z, zy, zyy); eta = z * inv.q / -inv.p; }
  return detail::pi_of_eta(P, inv.p, eta, pif);
}

static double pmean(const double s[6]) { return ((s[0] + s[1]) + s[2]) / 3.0; }
static double c_N(const Params& P) { return P.N == 0.0 ? std::exp(-1.0) : std::pow(1.0 - P.N, (1.0 - P.N) / P.N); }

// (S.58) at a committed state: (a) |p| >= p_min (1 - 1e-12), (b) F <= F_tol, (c) |pi_i| >= c_N p_min (1 - 1e-12 - d_F (1-N)/M)
static bool s58(const Params& P, const State& s, double& Fover, double& piratio)
{
  double sig[6]; stress(P, s, sig);
  double F = 0.0;
  if (!Fval(P, sig, s.pi_i, F)) return false;
  const double Ftol = detail::F_TRIAL_TOL_REL * pRef(P), dF = Ftol / P.p_min;
  Fover = F / P.p_min;
  piratio = std::fabs(s.pi_i) / (c_N(P) * P.p_min);
  return std::fabs(pmean(sig)) >= P.p_min * (1.0 - 1e-12) && F <= Ftol
      && piratio >= 1.0 - 1e-12 - dF * (1.0 - P.N) / P.M;
}

// K3 set of K1.16 (sheet §13.16; = Ladruno_files/testbed/norsand_calib/out/h19/h19_replay.cpp make_params)
static Params k3()
{
  Params P = paper();
  P.M = 1.25; P.N = 0.4686136860640983; P.N_bar = 0.16034970350694513;
  P.rho = 0.5828981510111576; P.rho_bar = 0.7213710496292179;
  P.chi = -2.8923322243122183; P.h = 308.3223484879116;
  P.e0 = 0.934; P.lambda_c = 0.019; P.xi = 0.7; P.p_a = 101.0; P.c1 = 0.05; P.c2 = 0.15;
  P.csl_mode = 1; P.zeta = 0; P.cap = 2; P.energy = 1;
  P.k = 298.34244992965336; P.g = 383.5831499095543; P.n_e = 0.5; P.p_min = 0.505;
  return P;
}

// TIMs HAR set of the §9.7 item 3 wet-floor table and of K1.14b (sheet "Inputs to rebuild", round 4b)
static Params tims_har()
{
  Params P = paper();
  P.energy = 1; P.k = 1889.48104361; P.g = 807.80387674; P.n_e = 0.5; P.p_a = 101.0;
  P.M = 1.3309; P.N = 0.4; P.N_bar = 0.2; P.chi = -3.5; P.h = 280.0;
  P.zeta = 0; P.rho = 0.71; P.rho_bar = 0.71;
  P.csl_mode = 1; P.e0 = 0.83; P.lambda_c = 0.027; P.xi = 0.45; P.cap = 0; P.p_min = 0.505;
  return P;
}

// A principal (diagonal) state on the surface: n^ = dev(d)/|dev(d)|, sigma_a = p + sqrt(2/3) q n^_a with
// q = eta |p| / zeta(theta(n^), rho), pi_i the inverse of (S.12) at (p, eta), eps^e by initialState's inversion.
static bool surface_state(const Params& P, const double d[3], double p, double eta, double v, State& s, double nh[3],
                          double* pi_out = nullptr)
{
  const double m = (d[0] + d[1] + d[2]) / 3.0;
  double nn = 0.0;
  for (int a = 0; a < 3; ++a) { nh[a] = d[a] - m; nn += nh[a] * nh[a]; }
  for (int a = 0; a < 3; ++a) nh[a] /= std::sqrt(nn);
  double t[3];
  for (int a = 0; a < 3; ++a) t[a] = -1.0 + std::sqrt(2.0 / 3.0) * nh[a];
  detail::Invariants inv; detail::invariants(t, inv);
  double z, zy, zyy; detail::zeta_y(inv.theta, P.rho, P.zeta, z, zy, zyy);
  const double q = eta * std::fabs(p) / z;
  double pi = 0.0;
  if (!detail::pi_of_eta(P, p, eta, pi)) return false;
  const double vv = std::isnan(v) ? 1.0 : v;
  const double sig0[6] = {p + std::sqrt(2.0 / 3.0) * q * nh[0], p + std::sqrt(2.0 / 3.0) * q * nh[1],
                          p + std::sqrt(2.0 / 3.0) * q * nh[2], 0.0, 0.0, 0.0};
  std::string msg;
  if (initialState(P, sig0, vv, pi, s, msg) != 0) { printf("INFO surface_state initialState: %s\n", msg.c_str()); return false; }
  if (pi_out) *pi_out = pi;
  return true;
}

// H19 dumped Gauss-point states (out/h19/s185.txt, stall.txt; format of out/h19/h19_replay.cpp):
// id ee11 ee22 eg12 s11 s22 s12 s33 pi_i v v0 epv eps at_floor n_f_tr n_f_post eps_f_v W_f d11 d22 dg12
struct H19 { const char* id; double ee[3], s[4], pi, v, v0, epv, eps; int at_floor, ntr, npost; double efv, wf, d[3]; };
static const H19 H19_STATES[] = {
  {"e7367g1_s185", {0.0020794563675880407, 0.0020623708231096217, -8.455253430988486e-05},
   {-0.5100767879440233, -0.6059940516753852, -0.2373365315574276, -0.39892916038058873}, -0.32249629546987785,
   1.9352522944653925, 1.69162, 0.11983998212133348, 0.41572595288512826, 1, 95, 9, 0.014565473339402796,
   0.007355564036398412, {0.0017300526378768999, -0.00021254327793723826, -0.009541568219878838}},
  {"e7367g2_s185", {0.0020837249895674934, 0.0020581170534215276, -8.706109638000796e-05},
   {-0.488250857597471, -0.6323752725912393, -0.24499494049350506, -0.3943738698112941}, -0.3666004636846257,
   1.9352522944653987, 1.69162, 0.12335248849525761, 0.41751582667471926, 1, 27, 5, 0.011051801984596942,
   0.005581160002221456, {0.0020580128620412974, -0.0005405035021015248, -0.008738158029346854}},
  {"e7367g1_m53", {0.002077423879239629, 0.00206727896553894, -4.942769003255963e-05},
   {-0.5074726482337255, -0.5635208707163949, -0.13653808447681343, -0.44400648104988016}, -0.11160661856019498,
   1.956532170893081, 1.69162, 0.12300707620431346, 0.4564006454713661, 1, 102, 24, 0.022341749899901803,
   0.011282583699450411, {0.0018391544802891235, -0.00023597356742177666, -0.01011036554856204}},
  {"e7367g2_m53", {0.00208013703803662, 0.002064364748965215, -5.1105036850271606e-05},
   {-0.4933360656283694, -0.5805613064996943, -0.14131268862820037, -0.4411026278719341}, -0.12495702194105537,
   1.9565321708930865, 1.69162, 0.126685149934657, 0.45547398396194977, 1, 37, 19, 0.01866320667248494,
   0.009424919369604896, {0.002214282120011929, -0.0006111012071446653, -0.00915474717917164}},
  {"e7775g1_m53", {0.0020801370380366226, 0.0020643647489652114, 5.1105036850273775e-05},
   {-0.49333606562835575, -0.5805613064997139, 0.1413126886282059, -0.4411026278719296}, -0.12495702194107475,
   1.9565321708930883, 1.69162, 0.12668514993465835, 0.45547398396195443, 1, 37, 19, 0.018663206672486155,
   0.009424919369605508, {0.002214282120011929, -0.0006111012071452482, 0.009154747179170308}},
  {"e7775g2_m53", {0.0020774238792396283, 0.002067278965538941, 4.942769003255269e-05},
   {-0.5074726482337295, -0.5635208707163885, 0.13653808447679222, -0.44400648104988627}, -0.11160661856016911,
   1.95653217089309, 1.69162, 0.12300707620431166, 0.4564006454713701, 1, 102, 24, 0.02234174989990671,
   0.011282583699452889, {0.0018391544802889015, -0.0002359735674218877, 0.01011036554856104}},
};

// the H19 state as a kernel State (eps^e_33 from sigma_33 by bisection on the kernel's stress(), h19_replay.cpp)
static bool h19_state(const Params& P, const H19& r, State& n)
{
  n = State();
  n.eps_e[0] = r.ee[0]; n.eps_e[1] = r.ee[1]; n.eps_e[2] = 0.0; n.eps_e[3] = 0.5 * r.ee[2]; n.eps_e[4] = 0.0; n.eps_e[5] = 0.0;
  n.pi_i = r.pi; n.v = r.v; n.v0 = r.v0; n.eps_p_v = r.epv; n.eps_p_s = r.eps; n.D_last = 0.0;
  n.eps_f_v = r.efv; n.W_f = r.wf; n.n_f_tr = r.ntr; n.n_f_post = r.npost; n.n_f_init = 0; n.at_floor = r.at_floor; n.n_f_pi = 0;
  auto f = [&](double e33) { n.eps_e[2] = e33; double sg[6]; stress(P, n, sg); return sg[2] - r.s[3]; };
  double lo = -0.2, hi = 0.2, flo = f(lo), fhi = f(hi);
  if (!(std::isfinite(flo) && std::isfinite(fhi)) || flo * fhi > 0) {
    bool ok = false; double a = -0.2, fa = f(a);
    for (int i = 1; i <= 4000 && !ok; i++) {
      double b = -0.2 + 0.4 * i / 4000.0, fb = f(b);
      if (std::isfinite(fa) && std::isfinite(fb) && fa * fb <= 0) { lo = a; hi = b; flo = fa; fhi = fb; ok = true; }
      a = b; fa = fb;
    }
    if (!ok) return false;
  }
  for (int it = 0; it < 200; it++) {
    double mid = 0.5 * (lo + hi), fm = f(mid);
    if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
  }
  n.eps_e[2] = 0.5 * (lo + hi);
  return true;
}

// the K1.16 increment: the listed (d11, d22, dg12) scaled to |de| = (d11^2 + d22^2 + dg12^2/2)^(1/2) = norm
static void h19_incr(const H19& r, double norm, double de[6])
{
  const double dn = std::sqrt(r.d[0] * r.d[0] + r.d[1] * r.d[1] + 0.5 * r.d[2] * r.d[2]);
  const double sc = norm / dn;
  de[0] = r.d[0] * sc; de[1] = r.d[1] * sc; de[2] = 0.0; de[3] = 0.5 * r.d[2] * sc; de[4] = 0.0; de[5] = 0.0;
}

// chained tangent of d with fractions fr against the central FD of the whole increment (chain_fd's measure),
// asserting at every FD point the base pattern counts (floor_tr, floor_post, floor_pi, plastic) of the base
static bool chain_fd_r4(const Params& P, const State& s, const double d[6], const std::vector<double>& fr, double h,
                        double& e)
{
  const int m = (int)fr.size();
  State np1; double sig[6], C[6][6]; StepInfo info;
  detail::step_fractions(P, s, d, fr.data(), m, true, np1, sig, C, info);
  bool ok = info.refusal == OK;
  e = 0.0;
  for (int J = 0; J < 6; ++J) {
    double dp[6], dm[6];
    for (int i = 0; i < 6; ++i) { dp[i] = d[i]; dm[i] = d[i]; }
    dp[J] += h; dm[J] -= h;
    State a, b; double sp[6], sm[6], Ca[6][6], Cb[6][6]; StepInfo ia, ib;
    detail::step_fractions(P, s, dp, fr.data(), m, true, a, sp, Ca, ia);
    detail::step_fractions(P, s, dm, fr.data(), m, true, b, sm, Cb, ib);
    for (const StepInfo* x : {&ia, &ib})
      ok = ok && x->refusal == OK && x->plastic == info.plastic && x->floor_tr == info.floor_tr
              && x->floor_post == info.floor_post && x->floor_pi == info.floor_pi;
    double n2 = 0.0, d2 = 0.0;
    for (int I = 0; I < 6; ++I) {
      const double fd = (sp[I] - sm[I]) / (2.0 * h);
      const double wgt = I < 3 ? 1.0 : 2.0;
      n2 += wgt * (C[I][J] - fd) * (C[I][J] - fd);
      d2 += wgt * C[I][J] * C[I][J];
    }
    e = std::fmax(e, std::sqrt(n2 / d2));
  }
  return ok;
}

static void round4_checks()
{
  // ---- K1.12 (C) (sheet §13.12 (C), round 4): BA06 wet side -Pf*, pi_i -40.198569107 -> -40.516314098 -------------
  {
    Params P = paper(); P.p_min = 50.0;                 // K2 set (mu0 5400, alpha0 0, WW rho 0.7 / 0.8, paper CSL), cap none
    const double dd[3] = {-1.2, -0.1, 1.3};
    State s; double nh[3], pi0 = 0.0;
    const bool built = surface_state(P, dd, -51.0, 0.6, 1.70, s, nh, &pi0);
    // sheet: pi_i(start) = -33.125471695 (tol 2e-9 kPa)
    check("r4_K1_12C_start_pi", fdiff(pi0, -33.125471695), 2e-9, built && fdiff(pi0, -33.125471695) <= 2e-9);
    const double de[6] = {6e-3 * 0.55 - 1e-5, 6e-3 * 0.40 - 1e-5, 6e-3 * -0.95 - 1e-5, 0.0, 0.0, 0.0};
    detail::OnceOut o;
    detail::step_once(P, s, de, o);
    detail::Elastic elc; detail::elastic(P, o.res.eps_conv, elc);
    detail::Invariants ic; detail::invariants(elc.sig, ic);
    double zc, zyc, zyyc; detail::zeta_y(ic.theta, P.rho, P.zeta, zc, zyc, zyyc);
    const double etac = zc * ic.q / -ic.p;
    double Fbe = 0.0, Fnew = 0.0, pif = 0.0, etaf = 0.0;
    Fval(P, o.sigma, o.res.pi_be, Fbe);
    Fval(P, o.sigma, o.st.pi_i, Fnew);
    pif_of(P, o.sigma, pif, etaf);
    // sheet values (tol): dlam 1.546466e-3 (rel 1e-6); p_c -48.764495, q_c 46.421589 (1e-6 kPa); eta_c 0.952607,
    // eta_f 0.929068 (1e-6); pi_BE -40.198569107 (2e-9); sigma_f (-34.8599764962, -34.1946812138, -80.9453422900)
    // (2e-9 kPa); F(sigma_f, pi_BE)/p_min +0.010899 (1e-6); pi_i -40.516314098 (2e-9); F(sigma_f, pi_i) = 0
    // (|F| <= 1e-12 p_min); pattern -Pf* (floor_tr 0, plastic, floor_post 1, floor_pi 1, m = 1); v unchanged.
    const double sf[3] = {-34.8599764962, -34.1946812138, -80.9453422900};
    double es = 0.0;
    for (int a = 0; a < 3; ++a) es = std::fmax(es, fdiff(o.sigma[a], sf[a]));
    printf("INFO r4 K1.12(C): dlam %.9e p_c %.9f q_c %.9f eta_c %.9f pi_BE %.12f eta_f %.9f F_BE/pmin %.9e pi_i %.12f "
           "F/pmin %.3e sigma_f (%.10f %.10f %.10f) pattern %d%d%d%d\n", o.res.dlam, elc.p, elc.q, etac, o.res.pi_be,
           etaf, Fbe / P.p_min, o.st.pi_i, Fnew / P.p_min, o.sigma[0], o.sigma[1], o.sigma[2], (int)o.res.floor_tr,
           (int)o.res.plastic, (int)o.res.floor_post, (int)o.res.floor_pi);
    const bool pat = !o.res.refused && !o.res.floor_tr && o.res.plastic && o.res.floor_post && o.res.floor_pi
                  && o.st.n_f_pi == 1 && o.st.v == s.v * std::exp(de[0] + de[1] + de[2]);
    check("r4_K1_12C_pattern_-Pf*", pat ? 0.0 : 1.0, 0.0, pat);
    check("r4_K1_12C_dlam", fdiff(o.res.dlam, 1.546466e-3) / 1.546466e-3, 1e-6, fdiff(o.res.dlam, 1.546466e-3) <= 1e-6 * 1.546466e-3);
    check("r4_K1_12C_p_c_q_c", std::fmax(fdiff(elc.p, -48.764495), fdiff(elc.q, 46.421589)), 1e-6,
          fdiff(elc.p, -48.764495) <= 1e-6 && fdiff(elc.q, 46.421589) <= 1e-6);
    check("r4_K1_12C_eta_c_eta_f", std::fmax(fdiff(etac, 0.952607), fdiff(etaf, 0.929068)), 1e-6,
          fdiff(etac, 0.952607) <= 1e-6 && fdiff(etaf, 0.929068) <= 1e-6);
    check("r4_K1_12C_pi_BE", fdiff(o.res.pi_be, -40.198569107), 2e-9, fdiff(o.res.pi_be, -40.198569107) <= 2e-9);
    check("r4_K1_12C_sigma_f", es, 2e-9, es <= 2e-9);
    check("r4_K1_12C_F_BE_over_pmin", fdiff(Fbe / P.p_min, 0.010899), 1e-6, fdiff(Fbe / P.p_min, 0.010899) <= 1e-6);
    check("r4_K1_12C_pi_i_new", fdiff(o.st.pi_i, -40.516314098), 2e-9, fdiff(o.st.pi_i, -40.516314098) <= 2e-9);
    check("r4_K1_12C_F_zero_after", std::fabs(Fnew) / P.p_min, 1e-12, std::fabs(Fnew) <= 1e-12 * P.p_min);
    // (S.58c) c_N p_min = 23.2379 kPa (1e-4)
    check("r4_K1_12C_cN_pmin", fdiff(c_N(P) * P.p_min, 23.2379), 1e-4, fdiff(c_N(P) * P.p_min, 23.2379) <= 1e-4);
  }

  // ---- §9.7 item 3, the wet-floor table with the rule (round 4): HAR TIMs set, one BE step, every case -Pf*, m = 1 ----
  {
    const Params P = tims_har();
    const double dd[3] = {-1.2, -0.1, 1.3};
    const double psis[4] = {0.0, 0.04, 0.08, 0.12};
    const double dgs[3] = {1e-4, 1e-3, 1e-2};
    // sheet table: F(sigma_f, pi_BE)/p_min (1e-4), pi_BE (1e-6 kPa), eta_f (1e-6), pi_i by (S.57) (2e-9 kPa)
    const double FB[4][3] = {{0.2699, 0.3496, 0.3572}, {0.2717, 0.3670, 0.5322}, {0.2732, 0.3820, 0.6866}, {0.2745, 0.3950, 0.8217}};
    const double PB[4][3] = {{-0.328889, -0.329228, -0.329262}, {-0.328346, -0.324406, -0.287273},
                             {-0.327887, -0.320380, -0.258876}, {-0.327494, -0.316976, -0.238442}};
    const double EF[4][3] = {{0.940109, 1.021599, 1.029388}, {0.938960, 1.012743, 0.951647},
                             {0.937988, 1.005324, 0.897090}, {0.937157, 0.999032, 0.856623}};
    const double PN[4][3] = {{-0.386215057, -0.406870322, -0.408940969}, {-0.385936497, -0.404537473, -0.389032266},
                             {-0.385701048, -0.402600177, -0.376005318}, {-0.385499896, -0.400969027, -0.366805331}};
    double eF = 0.0, eP = 0.0, eE = 0.0, eN = 0.0, eZ = 0.0, pi_start = 0.0;
    bool ok = true;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 3; ++j) {
        State s; double nh[3];
        double pi0 = 0.0;
        ok = ok && surface_state(P, dd, -0.505, 0.5 * P.M, std::nan(""), s, nh, &pi0);
        s.v = s.v0 = 1.0 + P.e0 - P.lambda_c * std::pow(-pi0 / P.p_a, P.xi) + psis[i];
        pi_start = pi0;
        const double a = dgs[j] / std::sqrt(2.0);
        const double de[6] = {a * nh[0], a * nh[1], a * nh[2], 0.0, 0.0, 0.0};
        detail::OnceOut o;
        detail::step_once(P, s, de, o);
        double Fbe = 0.0, Fn = 0.0, pif = 0.0, etaf = 0.0;
        Fval(P, o.sigma, o.res.pi_be, Fbe); Fval(P, o.sigma, o.st.pi_i, Fn); pif_of(P, o.sigma, pif, etaf);
        ok = ok && !o.res.refused && !o.res.floor_tr && o.res.plastic && o.res.floor_post && o.res.floor_pi;
        eF = std::fmax(eF, fdiff(Fbe / P.p_min, FB[i][j])); eP = std::fmax(eP, fdiff(o.res.pi_be, PB[i][j]));
        eE = std::fmax(eE, fdiff(etaf, EF[i][j])); eN = std::fmax(eN, fdiff(o.st.pi_i, PN[i][j]));
        eZ = std::fmax(eZ, std::fabs(Fn) / P.p_min);
        printf("INFO r4 wet-floor psi %.2f dg %.0e: F_BE/pmin %.4f pi_BE %.6f eta_f %.6f pi_i %.9f F/pmin %.1e\n", psis[i],
               dgs[j], Fbe / P.p_min, o.res.pi_be, etaf, o.st.pi_i, Fn / P.p_min);
      }
    // sheet: start pi_i = -0.328007 (1e-6)
    check("r4_wetfloor_start_pi", fdiff(pi_start, -0.328007), 1e-6, fdiff(pi_start, -0.328007) <= 1e-6);
    check("r4_wetfloor_all_-Pf*", ok ? 0.0 : 1.0, 0.0, ok);
    check("r4_wetfloor_F_BE", eF, 1e-4, eF <= 1e-4);
    check("r4_wetfloor_pi_BE", eP, 1e-6, eP <= 1e-6);
    check("r4_wetfloor_eta_f", eE, 1e-6, eE <= 1e-6);
    check("r4_wetfloor_pi_S57", eN, 2e-9, eN <= 2e-9);
    check("r4_wetfloor_F_zero_after", eZ, 1e-12, eZ <= 1e-12);
  }

  // ---- K1.14b (sheet §13.14b, round 4: no value changes; F(sigma_f, pi_BE) <= 0, the rule does NOT fire; M-R7) ----
  {
    const Params P = tims_har();
    const double dd[3] = {-1.0, -0.35, 1.35};
    State s; double nh[3], pi0 = 0.0;
    bool ok = surface_state(P, dd, -0.6, 1.2 * P.M, std::nan(""), s, nh, &pi0);
    s.v = s.v0 = 1.0 + P.e0 - P.lambda_c * std::pow(-pi0 / P.p_a, P.xi) - 0.10;
    // sheet: pi_i(start) -0.743658 (1e-6), v 1.72704 (1e-5)
    check("r4_K1_14b_start", std::fmax(fdiff(pi0, -0.743658), fdiff(s.v, 1.72704)), 1e-5,
          ok && fdiff(pi0, -0.743658) <= 1e-6 && fdiff(s.v, 1.72704) <= 1e-5);
    const double one[1] = {1.0};
    const double half[2] = {0.5, 0.5};
    double d1[6], d2[6], d3[6];
    for (int a = 0; a < 3; ++a) {
      d1[a] = 2e-5 / 3.0 + 2e-5 * nh[a];
      d2[a] = 2.0 * d1[a];
      d3[a] = 1.4e-5 / 3.0 + 2e-5 * nh[a];
    }
    for (int a = 3; a < 6; ++a) d1[a] = d2[a] = d3[a] = 0.0;
    State n1, n2, n3; double s1[6], s2[6], s3[6], C[6][6]; StepInfo i1, i2, i3;
    detail::step_fractions(P, s, d1, one, 1, true, n1, s1, C, i1);
    detail::step_fractions(P, s, d2, half, 2, true, n2, s2, C, i2);
    detail::step_fractions(P, s, d3, half, 2, true, n3, s3, C, i3);
    // sheet (tol 2e-9 kPa on pi_i, 2e-10 on sigma_f): single FPf pi -0.743350807, sigma_f (-0.8497174515,
    // -0.5790682544, -0.0862142940); chain FPf,FPf -0.743038708; chain -P-,-Pf -0.743573074; floor_pi = 0 in all
    const double sf[3] = {-0.8497174515, -0.5790682544, -0.0862142940};
    double es = 0.0;
    for (int a = 0; a < 3; ++a) es = std::fmax(es, fdiff(s1[a], sf[a]));
    printf("INFO r4 K1.14b: single pi %.12f (floor %d/%d/%d) sigma (%.10f %.10f %.10f); FPf,FPf pi %.12f (%d/%d/%d); "
           "-P-,-Pf pi %.12f (%d/%d/%d)\n", n1.pi_i, i1.floor_tr, i1.floor_post, i1.floor_pi, s1[0], s1[1], s1[2], n2.pi_i,
           i2.floor_tr, i2.floor_post, i2.floor_pi, n3.pi_i, i3.floor_tr, i3.floor_post, i3.floor_pi);
    const bool pat = i1.refusal == OK && i2.refusal == OK && i3.refusal == OK && i1.floor_tr == 1 && i1.floor_post == 1
                  && i2.floor_tr == 2 && i2.floor_post == 2 && i3.floor_tr == 0 && i3.floor_post == 1
                  && i1.floor_pi == 0 && i2.floor_pi == 0 && i3.floor_pi == 0 && n1.n_f_pi == 0 && n2.n_f_pi == 0;
    check("r4_K1_14b_patterns_rule_silent", pat ? 0.0 : 1.0, 0.0, pat);
    const double ep = std::fmax(fdiff(n1.pi_i, -0.743350807), std::fmax(fdiff(n2.pi_i, -0.743038708), fdiff(n3.pi_i, -0.743573074)));
    check("r4_K1_14b_pi_never_shrunk", ep, 2e-9, ep <= 2e-9);
    check("r4_K1_14b_sigma_f", es, 2e-10, es <= 2e-10);
  }

  // ---- K1.16 (sheet §13.16): the H19 regression points. K3 set, |de| = 2e-4 repeated, every step m = 1, -Pf* ----
  const Params P = k3();
  const double cNp = c_N(P) * P.p_min;
  check("r4_K1_16_cN_pmin", fdiff(cNp, 0.246558167), 1e-9, fdiff(cNp, 0.246558167) <= 1e-9);
  {
    struct Row { int step; double s[4], pi; };
    // sheet table (tol 5e-10 kPa = 1e-9 p_min on sigma and pi_i); eps^e_33 (tol 1e-14); |F| <= F_tol at every step
    const Row R1[3] = {{1, {-0.5101907038, -0.5976456953, -0.4071636009, -0.2151586912}, -0.3989652918},
                       {2, {-0.5082384097, -0.6104039680, -0.3963576223, -0.2510036052}, -0.4461726222},
                       {60, {-0.5114970385, -0.6473686614, -0.3561343001, -0.3336845167}, -0.5923861407}};
    const Row R2[3] = {{1, {-0.4889744155, -0.6273143992, -0.3987111852, -0.2332250428}, -0.4276028897},
                       {2, {-0.4850426738, -0.6393433615, -0.3906139647, -0.2596017804}, -0.4672749361},
                       {60, {-0.4814251232, -0.6764194559, -0.3571554210, -0.3278584883}, -0.5976134600}};
    const Row R3[2] = {{1, {-0.5071972881, -0.5355837599, -0.4722189520, -0.0691517825}, -0.2815026983},
                       {2, {-0.5031271070, -0.5765889736, -0.4352839194, -0.1789591553}, -0.3671148499}};
    const Row* rows[3] = {R1, R2, R3};
    const int nrows[3] = {3, 3, 2};
    const double e33[3] = {2.099254865382e-3, 2.100404994193e-3, 2.088911465473e-3};
    // sheet: F of the legacy starts +0.458 / (0.1472 kPa) / +3.29 p_min; step-1 detail of e7367g1_s185 below
    const char* names[3] = {"r4_K1_16_s185g1", "r4_K1_16_s185g2", "r4_K1_16_m53g1_stall"};
    for (int c = 0; c < 3; ++c) {
      State n;
      bool ok = h19_state(P, H19_STATES[c], n);
      check((std::string(names[c]) + "_eps_e33").c_str(), fdiff(n.eps_e[2], e33[c]), 1e-14, ok && fdiff(n.eps_e[2], e33[c]) <= 1e-14);
      double de[6]; h19_incr(H19_STATES[c], 2e-4, de);
      double err = 0.0, Fmax = 0.0, pmax = 0.0, pirel = 0.0;
      bool inv = ok, m1 = true;
      int r = 0;
      const int last = rows[c][nrows[c] - 1].step;
      for (int k = 1; k <= last && ok; ++k) {
        State np1; double sg[6], C[6][6]; StepInfo info;
        const int rc = step(P, n, de, np1, sg, C, info);
        double F = 0.0, pif = 0.0, etaf = 0.0;
        Fval(P, sg, np1.pi_i, F); pif_of(P, sg, pif, etaf);
        inv = inv && rc == OK && info.plastic && info.floor_post >= 1 && info.floor_pi >= 1
                  && np1.n_f_pi == n.n_f_pi + info.floor_pi
                  && F <= detail::F_TRIAL_TOL_REL * pRef(P) && std::fabs(np1.pi_i) >= cNp * (1.0 - 1e-9);
        m1 = m1 && info.substeps == 1 && info.floor_tr == 0 && info.floor_post == 1 && info.floor_pi == 1;
        Fmax = std::fmax(Fmax, std::fabs(F)); pmax = std::fmax(pmax, std::fabs(pmean(sg) + P.p_min));
        pirel = std::fmax(pirel, fdiff(np1.pi_i, pif) / std::fabs(pif));
        if (c == 0 && k == 1) {
          // sheet step-1 detail (tol): dlam 1.427733871e-4 (rel 1e-8), p_c -0.365744965 (1e-9), eta_c 1.092426 (1e-6),
          // pi_BE -0.3245384951 (2e-10), F(sigma_f, pi_BE) +0.175979368 kPa (1e-9), eta_f 0.922543205 (1e-9),
          // deps^f_v 7.223598e-5 (rel 1e-6)
          detail::OnceOut o; detail::step_once(P, n, de, o);
          detail::Elastic elc; detail::elastic(P, o.res.eps_conv, elc);
          detail::Invariants ic; detail::invariants(elc.sig, ic);
          double zc, zyc, zyyc; detail::zeta_y(ic.theta, P.rho, P.zeta, zc, zyc, zyyc);
          double Fbe = 0.0; Fval(P, o.sigma, o.res.pi_be, Fbe);
          printf("INFO r4 K1.16 s185g1 step 1: dlam %.10e p_c %.10f eta_c %.7f pi_BE %.11f F_BE %.10f eta_f %.10f dfv %.7e\n",
                 o.res.dlam, elc.p, zc * ic.q / -ic.p, o.res.pi_be, Fbe, etaf, o.res.dfv_post);
          const double e1 = std::fmax(std::fmax(fdiff(o.res.dlam, 1.427733871e-4) / 1.427733871e-4 / 1e-8, fdiff(elc.p, -0.365744965) / 1e-9),
                            std::fmax(std::fmax(fdiff(zc * ic.q / -ic.p, 1.092426) / 1e-6, fdiff(o.res.pi_be, -0.3245384951) / 2e-10),
                            std::fmax(std::fmax(fdiff(Fbe, 0.175979368) / 1e-9, fdiff(etaf, 0.922543205) / 1e-9),
                                      fdiff(o.res.dfv_post, 7.223598e-5) / 7.223598e-5 / 1e-6)));
          check("r4_K1_16_s185g1_step1_detail_(units_of_tol)", e1, 1.0, e1 <= 1.0);
        }
        if (r < nrows[c] && rows[c][r].step == k) {
          const double got[4] = {sg[0], sg[1], sg[2], sg[3]};
          for (int a = 0; a < 4; ++a) err = std::fmax(err, fdiff(got[a], rows[c][r].s[a]));
          err = std::fmax(err, fdiff(np1.pi_i, rows[c][r].pi));
          printf("INFO r4 %s step %d: s11 %.10f s22 %.10f s33 %.10f s12 %.10f pi %.10f F %.2e m %d pattern %d%d%d%d\n",
                 names[c], k, sg[0], sg[1], sg[2], sg[3], np1.pi_i, F, info.substeps, info.floor_tr, info.plastic,
                 info.floor_post, info.floor_pi);
          ++r;
        }
        n = np1;
      }
      check((std::string(names[c]) + "_invariants_S58").c_str(), Fmax, detail::F_TRIAL_TOL_REL * pRef(P),
            inv && pmax <= 1e-12 * P.p_min && pirel <= 1e-12);
      // the stall row (c = 2): digits only if the kernel took every step at m = 1 (sheet §13.16 round 4b)
      if (c < 2 || m1) check((std::string(names[c]) + "_digits").c_str(), err, 5e-10, m1 && err <= 5e-10);
      else printf("INFO r4 %s: the kernel substepped; digits not asserted (invariants only, sheet 13.16)\n", names[c]);
    }
  }

  // ---- (S.58) census: every committed state after an accepted post-floored step, from the six H19 states ---------
  // (legacy starts: (S.58) applies after their first accepted step). Marches of 2e-4 along the listed direction
  // (600 steps from e7367g1_s185, 60 from the others) and 48 in-plane directions x |de| in {1e-3, 2e-4, 1e-5}, two
  // steps each. Sheet: F(sigma_f, pi_i) = 0 to 2.5e-14 p_min over 1800 reset steps; min |pi_i| 0.381 ... 0.399 kPa.
  {
    int nev = 0, nfire = 0, nquiet = 0, nbad = 0, nref = 0, nsubst = 0;
    double Fmax = -1e300, prmin = 1e300, pimin = 1e300;
    auto census = [&](const State& s, const StepInfo& info) {
      if (info.refusal != OK || info.floor_post == 0) return;
      ++nev; nfire += info.floor_pi > 0 ? 1 : 0; nquiet += info.floor_pi == 0 ? 1 : 0; nsubst += info.substeps > 1;
      double Fo = 0.0, pr = 0.0;
      if (!s58(P, s, Fo, pr)) ++nbad;
      Fmax = std::fmax(Fmax, Fo); prmin = std::fmin(prmin, pr); pimin = std::fmin(pimin, std::fabs(s.pi_i));
    };
    for (int c = 0; c < 6; ++c) {
      State n0;
      if (!h19_state(P, H19_STATES[c], n0)) { ++nbad; continue; }
      double de[6]; h19_incr(H19_STATES[c], 2e-4, de);
      State n = n0;
      const int steps = c == 0 ? 600 : 60;
      for (int k = 0; k < steps; ++k) {
        State np1; double sg[6], C[6][6]; StepInfo info;
        step(P, n, de, np1, sg, C, info);
        if (info.refusal != OK) { ++nref; break; }
        if (k > 0) census(np1, info);
        n = np1;
      }
      for (int i = 0; i < 12; ++i)
        for (int j = -1; j <= 2; ++j)
          for (double nr : {1e-3, 2e-4, 1e-5}) {
            const double PI = 3.14159265358979323846;
            const double az = 2 * PI * i / 12.0, el = PI * j / 6.0;
            const double u[3] = {std::cos(el) * std::cos(az), std::cos(el) * std::sin(az), std::sin(el)};
            const double d6[6] = {u[0] * nr, u[1] * nr, 0.0, u[2] * nr / std::sqrt(2.0), 0.0, 0.0};
            State a = n0;
            for (int k = 0; k < 2; ++k) {
              State np1; double sg[6], C[6][6]; StepInfo info;
              step(P, a, d6, np1, sg, C, info);
              if (info.refusal != OK) { ++nref; break; }
              if (k > 0 || a.n_f_pi > 0) census(np1, info);   // a legacy start: (S.58) after the first accepted reset
              a = np1;
            }
          }
    }
    printf("INFO r4 (S.58) census: %d post-floored committed states (%d with the rule fired, %d without; %d substepped), "
           "%d refused increments; max F/p_min %.3e, min |pi_i|/(c_N p_min) %.6f, min |pi_i| %.6f kPa, violations %d\n",
           nev, nfire, nquiet, nsubst, nref, Fmax, prmin, pimin, nbad);
    check("r4_S58_census_F_le_Ftol_pi_ge_cN_pmin", (double)nbad, 0.0, nbad == 0 && nev >= 1000 && nfire >= 600);
  }

  // ---- (S.59) gate (sheet §9.7 "Gate that makes (S.59) testable", round 4b) ---------------------------------------
  // start = e7367g1_s185 after ONE K1.16 step (pi_i -0.3989652918); de = 2x the K1.16 increment (|de| 4e-4); fractions
  // (1/2, 1/2); pattern -Pf*,-Pf* with F(sigma_f, pi_BE) = +0.162 and +0.099 p_min (outside the noise band, 1e-3);
  // committed pi_i -0.479306004 (2e-9), p -0.505. Chained C vs central FD (sheet): 1.0e-5 / 1.0e-7 / 1.3e-9 at
  // h = 1e-6 / 1e-7 / 1e-8; PASS <= 3e-7 at h = 1e-7 (the (S.54) row kept: 0.20). Also case (a): m = 1, |de| 2e-4:
  // 3.0e-5 / 3.0e-7 / 3.0e-9 (pass <= 3e-6 at 1e-7), and the chained m = 1 tangent = (S.32f) to round-off.
  {
    State n;
    bool ok = h19_state(P, H19_STATES[0], n);
    double d1[6]; h19_incr(H19_STATES[0], 2e-4, d1);
    { State np1; double sg[6], C[6][6]; StepInfo info; step(P, n, d1, np1, sg, C, info); ok = ok && info.refusal == OK; n = np1; }
    check("r4_S59_fixture_start_pi", fdiff(n.pi_i, -0.3989652918), 5e-10, ok && fdiff(n.pi_i, -0.3989652918) <= 5e-10);
    double d2[6], dh[6];
    for (int i = 0; i < 6; ++i) { d2[i] = 2.0 * d1[i]; dh[i] = 0.5 * d2[i]; }
    detail::OnceOut o1, o2;
    detail::step_once(P, n, dh, o1);
    detail::step_once(P, o1.st, dh, o2);
    double F1 = 0.0, F2 = 0.0;
    Fval(P, o1.sigma, o1.res.pi_be, F1); Fval(P, o2.sigma, o2.res.pi_be, F2);
    const std::vector<double> fr = {0.5, 0.5};
    State nf; double sf[6], Cf[6][6]; StepInfo inf;
    detail::step_fractions(P, n, d2, fr.data(), 2, true, nf, sf, Cf, inf);
    printf("INFO r4 S.59 fixture: F_BE/pmin %.4f %.4f, pattern floor %d/%d/%d plastic %d, pi_i %.10f p %.12f\n",
           F1 / P.p_min, F2 / P.p_min, inf.floor_tr, inf.floor_post, inf.floor_pi, inf.plastic, nf.pi_i, pmean(sf));
    const bool pat = inf.refusal == OK && inf.floor_tr == 0 && inf.floor_post == 2 && inf.floor_pi == 2
                  && o1.res.floor_pi && o2.res.floor_pi && F1 >= 1e3 * detail::F_TRIAL_TOL_REL * pRef(P)
                  && F2 >= 1e3 * detail::F_TRIAL_TOL_REL * pRef(P);
    check("r4_S59_fixture_pattern_-Pf*,-Pf*_outside_band", pat ? 0.0 : 1.0, 0.0, pat);
    check("r4_S59_fixture_F_BE_0.162_0.099", std::fmax(fdiff(F1 / P.p_min, 0.162), fdiff(F2 / P.p_min, 0.099)), 1e-3,
          fdiff(F1 / P.p_min, 0.162) <= 1e-3 && fdiff(F2 / P.p_min, 0.099) <= 1e-3);
    check("r4_S59_fixture_pi_i", fdiff(nf.pi_i, -0.479306004), 2e-9, fdiff(nf.pi_i, -0.479306004) <= 2e-9 &&
          fdiff(pmean(sf), -P.p_min) <= 1e-12 * P.p_min);
    double e[3] = {0.0, 0.0, 0.0};
    const double hs[3] = {1e-6, 1e-7, 1e-8};
    bool okfd = true;
    for (int i = 0; i < 3; ++i) okfd = chain_fd_r4(P, n, d2, fr, hs[i], e[i]) && okfd;
    printf("INFO r4 S.59 chain (1/2,1/2) C vs FD: %.2e / %.2e / %.2e at h = 1e-6 / 1e-7 / 1e-8 (sheet 1.0e-5 / 1.0e-7 / 1.3e-9)\n",
           e[0], e[1], e[2]);
    check("r4_S59_chain_vs_fd_h1e-7", e[1], 3e-7, okfd && e[1] <= 3e-7);
    check("r4_S59_chain_vs_fd_O(h2)", e[0] / e[1], 30.0, okfd && e[0] / e[1] >= 30.0);
    // case (a): m = 1
    double ea[3] = {0.0, 0.0, 0.0};
    const std::vector<double> f1 = {1.0};
    bool oka = true;
    for (int i = 0; i < 3; ++i) oka = chain_fd_r4(P, n, d1, f1, hs[i], ea[i]) && oka;
    printf("INFO r4 S.59 case (a) m = 1 C vs FD: %.2e / %.2e / %.2e (sheet 3.0e-5 / 3.0e-7 / 3.0e-9)\n", ea[0], ea[1], ea[2]);
    check("r4_S59_case_a_m1_vs_fd_h1e-7", ea[1], 3e-6, oka && ea[1] <= 3e-6);
    State a, b; double sa[6], sb[6], Ca[6][6], Cb[6][6]; StepInfo ia, ib;
    detail::step_fractions(P, n, d1, f1.data(), 1, true, a, sa, Ca, ia);
    detail::step_fractions(P, n, d1, f1.data(), 1, false, b, sb, Cb, ib);
    double dC = 0.0;
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) dC = std::fmax(dC, std::fabs(Ca[I][J] - Cb[I][J]));
    dC /= maxabs(&Cb[0][0], 36);
    check("r4_S59_case_a_chain_m1_equals_S32f", dC, 1e-12, ia.floor_pi == 1 && ib.floor_pi == 1 && dC <= 1e-12);
  }
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
    const double wild[5][6] = {{1e4, -1e4, 5e3, 2e3, -3e3, 1e3},     // |deps| ~ 1e4: elastic overflow
                               {-1.0, -1.0, -1.0, 0.3, 0.0, 0.0},      // huge finite compression
                               {1.0, 1.0, 1.0, 0.0, 0.0, 0.0},         // tension: p > 0 at the trial
                               {0.2, -0.3, 0.1, 0.25, -0.1, 0.05},     // large deviatoric
                               {300.0, 300.0, 300.0, 0.0, 0.0, 0.0}};  // tr = 900: exp(tr) overflows in v
    for (int w = 0; w < 5; ++w) {
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
      if (w == 4) check("v_overflow_increment_refused", info.refusal, 0, info.refusal != OK && np1.v == s.v);
    }
    State np1; double sig[6], C[6][6]; StepInfo info;
    const double dn[6] = {nan, 0, 0, 0, 0, 0};
    step(P, s, dn, np1, sig, C, info);
    check("nan_increment_refused", info.refusal, 0, info.refusal != OK);
  }

  // ---- CHECK: exponential specific-volume law (sheet §1.2 (S.26), G2 owner decision 2026-10-01) --
  // v_{n+1} = v_n exp(tr deps) on every accepted increment (bit-exact on a non-substepped one, a product of
  // m sub-increment exponentials on a substepped one), hence v = v0 exp(tr eps) along the path. A linear
  // law v_n + v0 tr deps is off by ~v0 (tr eps)^2 / 2 ~ 1e-4 relative on these paths.
  {
    double worst_step = 0.0, worst_path = 0.0, lin_gap = 0.0;
    int nsub = 0, nacc = 0;
    bool ok = true;
    for (int c = 0; c < 2; ++c) {
      Params P = c == 0 ? fork() : paper();
      if (c == 1) { P.cap = 2; P.c1 = 0.05; P.c2 = 0.15; }
      State s = c == 0 ? init_state(P, -100.0, -105.0, 0.01) : init_state(P, -100.0, -80.0, -0.05);
      const int n = c == 0 ? 25 : 40;
      std::vector<std::vector<double>> path = c == 0 ? noncoax(n) : std::vector<std::vector<double>>();
      if (c == 1) {
        const double d[6] = {(-0.01 + 2e-3) / n, -0.01 / n, (-0.01 - 2e-3) / n, 0, 0, 0};
        path = repeat(d, n);
      }
      double trsum = 0.0;
      for (int k = 0; k < n; ++k) {
        const double* d = path[k].data();
        State np1; double sig[6], C[6][6]; StepInfo info;
        step(P, s, d, np1, sig, C, info);
        if (info.refusal != OK) { ok = false; break; }
        ++nacc;
        const double tr = (d[0] + d[1]) + d[2];
        trsum += tr;
        const double es = std::fabs(np1.v - s.v * std::exp(tr)) / s.v;
        if (info.substeps == 1) ok = ok && np1.v == s.v * std::exp(tr);
        else ++nsub;
        worst_step = std::fmax(worst_step, es);
        worst_path = std::fmax(worst_path, std::fabs(np1.v - np1.v0 * std::exp(trsum)) / np1.v0);
        lin_gap = std::fmax(lin_gap, std::fabs(np1.v - np1.v0 * (1.0 + trsum)) / np1.v0);
        ok = ok && np1.v0 == s.v0;
        s = np1;
      }
    }
    printf("INFO v-law: %d accepted increments (%d substepped): per-step |v - v_n exp(tr)|/v max %.3e, "
           "path |v - v0 exp(tr eps)|/v0 max %.3e, linear-law gap %.3e\n", nacc, nsub, worst_step, worst_path, lin_gap);
    check("v_exponential_per_step", worst_step, 1e-14, ok && nsub >= 20 && worst_step <= 1e-14);
    check("v_equals_v0_exp_tr_eps", worst_path, 1e-13, ok && worst_path <= 1e-13 && lin_gap > 1e-6);
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

  // ---- CHECK: chained consistent tangent across sub-increments (sheet §9.6) --------------
  // (A) AMP_STOP smooth cap n = 40: every substepped ladder increment. The ladder tangent must be
  //     bit-identical to step_fractions(1/m, chain) and match the central FD of the WHOLE increment
  //     (fractions held fixed); the last sub-increment's CTO must be far from it (O2: 0.51-0.90).
  {
    Params P = paper(); P.cap = 2; P.c1 = 0.05; P.c2 = 0.15;
    State s = init_state(P, -100.0, -80.0, -0.05);
    const int n = 40;
    const double d[6] = {(-0.01 + 2e-3) / n, -0.01 / n, (-0.01 - 2e-3) / n, 0, 0, 0};
    int nsub = 0;
    double worst = 0.0, last_min = 1e300;
    bool ident = true, ok = true;
    for (int k = 0; k < n; ++k) {
      State np1; double sig[6], C[6][6]; StepInfo info;
      step(P, s, d, np1, sig, C, info);
      if (info.refusal != OK) { ok = false; break; }
      if (info.substeps > 1) {
        ++nsub;
        const int m = info.substeps;
        std::vector<double> fr(m, 1.0 / m);
        double e = 0.0, el = 0.0;
        if (!chain_fd(P, s, d, fr, 1e-7, e, el)) ok = false;
        worst = std::fmax(worst, e);
        last_min = std::fmin(last_min, el);
        State b; double sb[6], Cb[6][6]; StepInfo ib;
        detail::step_fractions(P, s, d, fr.data(), m, true, b, sb, Cb, ib);
        for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) ident = ident && Cb[I][J] == C[I][J];
      }
      s = np1;
    }
    printf("INFO chain AMP_STOP: %d substepped increments, chain vs FD (h = 1e-7) max %.3e, last-sub CTO min %.3f\n",
           nsub, worst, last_min);
    check("chain_ampstop_vs_fd", worst, 1e-6, ok && nsub >= 20 && worst <= 1e-6);
    check("chain_ladder_equals_step_fractions", ident ? 0.0 : 1.0, 0.0, ok && ident);
    check("chain_last_substep_cto_is_far", last_min, 0.1, ok && last_min > 0.1);
  }
  // (B, E) generic non-coaxial plastic increment (fork WW, three shears): forced uniform m = 8, 2 and
  //     the non-uniform (recursive-halving) fractions (1/2, 1/4, 1/8, 1/8); (C) m = 1 chain = (S.33).
  {
    Params P = fork(); P.rho_bar = 0.71;
    const double sh[6] = {0, 0, 0, 3e-4, 2e-4, 1e-4};
    State s = init_state(P, -100.0, -60.4, 0.0);
    s.v = s.v0 = 1.65;
    const double pre[6] = {4e-4 + sh[0], -1e-3, 0.0, sh[3], sh[4], sh[5]};
    for (int k = 0; k < 5; ++k) {
      State np1; double sig[6], C[6][6]; StepInfo info;
      step(P, s, pre, np1, sig, C, info);
      s = np1;
    }
    const double db[6] = {1e-4, -6e-4, 2e-4, 0.5 * sh[3], 0.5 * sh[4], 0.5 * sh[5]};
    const std::vector<std::vector<double>> frs = {std::vector<double>(8, 0.125), {0.5, 0.5},
                                                  {0.5, 0.25, 0.125, 0.125}};
    const char* names[3] = {"chain_generic_m8_vs_fd", "chain_generic_m2_vs_fd", "chain_generic_nonuniform_vs_fd"};
    for (int c = 0; c < 3; ++c) {
      double e = 0.0, el = 0.0;
      const bool okc = chain_fd(P, s, db, frs[c], 1e-7, e, el);
      printf("INFO %s: chain err %.3e, last-sub CTO err %.3e\n", names[c], e, el);
      check(names[c], e, 1e-6, okc && e <= 1e-6 && el > 1e-3);
    }
    // (C) m = 1: chain vs (S.33) (distinct trial eigenvalues: equal to round-off)
    State a, b; double sa[6], sb[6], Ca[6][6], Cb[6][6]; StepInfo ia, ib;
    const double one = 1.0;
    detail::step_fractions(P, s, db, &one, 1, true, a, sa, Ca, ia);
    step(P, s, db, b, sb, Cb, ib);
    double dC = 0.0;
    for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) dC = std::fmax(dC, std::fabs(Ca[I][J] - Cb[I][J]));
    dC /= maxabs(&Cb[0][0], 36);
    check("chain_m1_equals_S33", dC, 1e-12, ia.refusal == OK && ib.substeps == 1 && ib.plastic && dC <= 1e-12);
  }
  // (F) vertex: hydrostatic plastic step from the apex, fractions (1/2, 1/4, 1/4): C:1 = 0
  {
    const Params P = paper();
    State s = init_state(P, -100.0, nan, 0.0);
    s.v = s.v0 = 1.59;
    const double dv[6] = {-1e-3, -1e-3, -1e-3, 0, 0, 0};
    const double fr[3] = {0.5, 0.25, 0.25};
    State np1; double sig[6], C[6][6]; StepInfo info;
    detail::step_fractions(P, s, dv, fr, 3, true, np1, sig, C, info);
    double c1 = 0.0;
    for (int I = 0; I < 6; ++I) c1 = std::fmax(c1, std::fabs(C[I][0] + C[I][1] + C[I][2]));
    c1 /= maxabs(&C[0][0], 36);
    check("chain_vertex_C_on_1_is_zero", c1, 1e-12,
          info.refusal == OK && info.vertex && np1.pi_i == s.pi_i && c1 <= 1e-12);
  }

  // ---- round 4: the image-pressure rule at a post floor (sheet §9.7 (S.57)-(S.59), §13.12 (C), §13.14b, §13.16) --
  round4_checks();

  printf("SUMMARY %s (%d failed)\n", g_fail ? "FAIL" : "PASS", g_fail);
  return g_fail ? 1 : 0;
}
