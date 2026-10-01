// extern "C" shim around the LadrunoNORSAND kernel for the ctypes parity harness (WP-144 P1a).
// Test-bed code (not part of the OpenSees build). Compiled by ns_kernel.py:
//   g++ -std=c++17 -O2 -fPIC -shared -I <repo>/SRC/material/nD ns_shim.cpp -o libns_kernel.so
//
// Flat layouts (all double unless stated):
//   params d[20] = p0, kappa_hat, eps_v0, mu0, alpha0, M, N, N_bar, rho, rho_bar, chi, h,
//                  lambda_tilde, v_c0, e0, lambda_c, xi, p_a, c1, c2;   params i[3] = csl_mode, zeta, cap
//   state[12]    = eps_e[6] {00,11,22,01,12,02}, pi_i, v, v0, eps_p_v, eps_p_s, D_last
//   info int[9]  = refusal, plastic, vertex, cap_active, local_iters, pi_iters, substeps, finest, finest_sub
//   C[36]        = row-major C[I][J] = d sigma_I / d eps_J (tensor shear, see the kernel header)
#include <LadrunoNorSandKernel.h>

#include <cstring>
#include <string>

using namespace ladruno_norsand;

namespace {

Params mk(const double* d, const int* i)
{
  Params P;
  P.p0 = d[0]; P.kappa_hat = d[1]; P.eps_v0 = d[2]; P.mu0 = d[3]; P.alpha0 = d[4]; P.M = d[5];
  P.N = d[6]; P.N_bar = d[7]; P.rho = d[8]; P.rho_bar = d[9]; P.chi = d[10]; P.h = d[11];
  P.lambda_tilde = d[12]; P.v_c0 = d[13]; P.e0 = d[14]; P.lambda_c = d[15]; P.xi = d[16];
  P.p_a = d[17]; P.c1 = d[18]; P.c2 = d[19];
  P.csl_mode = i[0]; P.zeta = i[1]; P.cap = i[2];
  return P;
}

void s2a(const State& s, double* a)
{
  for (int k = 0; k < 6; ++k) a[k] = s.eps_e[k];
  a[6] = s.pi_i; a[7] = s.v; a[8] = s.v0; a[9] = s.eps_p_v; a[10] = s.eps_p_s; a[11] = s.D_last;
}

State a2s(const double* a)
{
  State s;
  for (int k = 0; k < 6; ++k) s.eps_e[k] = a[k];
  s.pi_i = a[6]; s.v = a[7]; s.v0 = a[8]; s.eps_p_v = a[9]; s.eps_p_s = a[10]; s.D_last = a[11];
  return s;
}

void copy_msg(const std::string& m, char* out, int len)
{
  if (!out || len <= 0) return;
  std::strncpy(out, m.c_str(), static_cast<size_t>(len - 1));
  out[len - 1] = '\0';
}

}  // namespace

extern "C" {

int ns_validate(const double* d, const int* i, char* msg, int msglen, int* warn)
{
  std::string m;
  bool w = false;
  const int rc = validate(mk(d, i), m, w);
  copy_msg(m, msg, msglen);
  *warn = w ? 1 : 0;
  return rc;
}

int ns_initial_state(const double* d, const int* i, const double* sigma0, double v0, double pi0,
                     double* st, char* msg, int msglen)
{
  State s;
  std::string m;
  const int rc = initialState(mk(d, i), sigma0, v0, pi0, s, m);
  copy_msg(m, msg, msglen);
  if (rc == 0) s2a(s, st);
  return rc;
}

int ns_step(const double* d, const int* i, const double* st_n, const double* deps, double* st_np1,
            double* sigma, double* C, int* info)
{
  const Params P = mk(d, i);
  const State n = a2s(st_n);
  State np1;
  double Cm[6][6];
  StepInfo inf;
  int finest = 0, sub = 0;
  const int rc = detail::step_ex(P, n, deps, np1, sigma, Cm, inf, finest, sub);
  s2a(np1, st_np1);
  for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) C[6 * I + J] = Cm[I][J];
  info[0] = inf.refusal; info[1] = inf.plastic; info[2] = inf.vertex; info[3] = inf.cap_active;
  info[4] = inf.local_iters; info[5] = inf.pi_iters; info[6] = inf.substeps; info[7] = finest; info[8] = sub;
  // StepInfo.finest / finest_sub must be the same numbers as step_ex's out-parameters
  if (inf.finest != finest || inf.finest_sub != sub) info[7] = info[8] = -1;
  return rc;
}

// detail::step_fractions (O2 api.step_fractions): m sub-increments fr[k] * deps, no ladder;
// chain != 0: the chained tangent (S.47) for every m (m = 1 included).
int ns_step_fractions(const double* d, const int* i, const double* st_n, const double* deps, const double* fr,
                      int m, int chain, double* st_np1, double* sigma, double* C, int* info)
{
  const Params P = mk(d, i);
  const State n = a2s(st_n);
  State np1;
  double Cm[6][6];
  StepInfo inf;
  const int rc = detail::step_fractions(P, n, deps, fr, m, chain != 0, np1, sigma, Cm, inf);
  s2a(np1, st_np1);
  for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) C[6 * I + J] = Cm[I][J];
  info[0] = inf.refusal; info[1] = inf.plastic; info[2] = inf.vertex; info[3] = inf.cap_active;
  info[4] = inf.local_iters; info[5] = inf.pi_iters; info[6] = inf.substeps; info[7] = inf.finest;
  info[8] = inf.finest_sub;
  return rc;
}

void ns_stress(const double* d, const int* i, const double* st, double* sigma)
{
  stress(mk(d, i), a2s(st), sigma);
}

void ns_elastic_tangent(const double* d, const int* i, const double* st, double* C)
{
  double Cm[6][6];
  elasticTangent(mk(d, i), a2s(st), Cm);
  for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) C[6 * I + J] = Cm[I][J];
}

}  // extern "C"
