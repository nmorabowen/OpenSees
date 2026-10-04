// extern "C" shim around the LadrunoNORSAND kernel for the ctypes parity harness (WP-144 P1a).
// Test-bed code (not part of the OpenSees build). Compiled by ns_kernel.py:
//   g++ -std=c++17 -O2 -fPIC -shared -I <repo>/SRC/material/nD ns_shim.cpp -o libns_kernel.so
//
// Flat layouts (all double unless stated):
//   params d[24] = p0, kappa_hat, eps_v0, mu0, alpha0, M, N, N_bar, rho, rho_bar, chi, h,
//                  lambda_tilde, v_c0, e0, lambda_c, xi, p_a, c1, c2, k, g, n_e, p_min;
//   params i[4]  = csl_mode, zeta, cap, energy
//   state[18]    = eps_e[6] {00,11,22,01,12,02}, pi_i, v, v0, eps_p_v, eps_p_s, D_last,
//                  eps_f_v, W_f, n_f_tr, n_f_post, n_f_init, at_floor     (round 3b: the floor counters)
//   info int[12] = refusal, plastic, vertex, cap_active, local_iters, pi_iters, substeps, finest, finest_sub,
//                  floor_tr, floor_post, at_floor
//   infod[3]     = deps_f_v, W_f, E_f
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
  P.k = d[20]; P.g = d[21]; P.n_e = d[22]; P.p_min = d[23];
  P.csl_mode = i[0]; P.zeta = i[1]; P.cap = i[2]; P.energy = i[3];
  return P;
}

void s2a(const State& s, double* a)
{
  for (int k = 0; k < 6; ++k) a[k] = s.eps_e[k];
  a[6] = s.pi_i; a[7] = s.v; a[8] = s.v0; a[9] = s.eps_p_v; a[10] = s.eps_p_s; a[11] = s.D_last;
  a[12] = s.eps_f_v; a[13] = s.W_f; a[14] = s.n_f_tr; a[15] = s.n_f_post; a[16] = s.n_f_init; a[17] = s.at_floor;
}

State a2s(const double* a)
{
  State s;
  for (int k = 0; k < 6; ++k) s.eps_e[k] = a[k];
  s.pi_i = a[6]; s.v = a[7]; s.v0 = a[8]; s.eps_p_v = a[9]; s.eps_p_s = a[10]; s.D_last = a[11];
  s.eps_f_v = a[12]; s.W_f = a[13];
  s.n_f_tr = static_cast<int>(a[14]); s.n_f_post = static_cast<int>(a[15]);
  s.n_f_init = static_cast<int>(a[16]); s.at_floor = static_cast<int>(a[17]);
  return s;
}

void info2a(const StepInfo& inf, int* info, double* infod)
{
  info[0] = inf.refusal; info[1] = inf.plastic; info[2] = inf.vertex; info[3] = inf.cap_active;
  info[4] = inf.local_iters; info[5] = inf.pi_iters; info[6] = inf.substeps; info[7] = inf.finest;
  info[8] = inf.finest_sub; info[9] = inf.floor_tr; info[10] = inf.floor_post; info[11] = inf.at_floor;
  infod[0] = inf.deps_f_v; infod[1] = inf.W_f; infod[2] = inf.E_f;
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

int ns_initial_state(const double* d, const int* i, const double* sigma0, double v0, double pi0, int pi0_rule,
                     double* st, char* msg, int msglen)
{
  State s;
  std::string m;
  const int rc = initialState(mk(d, i), sigma0, v0, pi0, s, m, pi0_rule);
  copy_msg(m, msg, msglen);
  if (rc == 0) s2a(s, st);
  return rc;
}

int ns_step(const double* d, const int* i, const double* st_n, const double* deps, double* st_np1,
            double* sigma, double* C, int* info, double* infod)
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
  info2a(inf, info, infod);
  // StepInfo.finest / finest_sub must be the same numbers as step_ex's out-parameters
  if (inf.finest != finest || inf.finest_sub != sub) info[7] = info[8] = -1;
  return rc;
}

// detail::step_fractions (O2 api.step_fractions): m sub-increments fr[k] * deps, no ladder;
// chain != 0: the chained tangent (S.47) for every m (m = 1 included).
int ns_step_fractions(const double* d, const int* i, const double* st_n, const double* deps, const double* fr,
                      int m, int chain, double* st_np1, double* sigma, double* C, int* info, double* infod)
{
  const Params P = mk(d, i);
  const State n = a2s(st_n);
  State np1;
  double Cm[6][6];
  StepInfo inf;
  const int rc = detail::step_fractions(P, n, deps, fr, m, chain != 0, np1, sigma, Cm, inf);
  s2a(np1, st_np1);
  for (int I = 0; I < 6; ++I) for (int J = 0; J < 6; ++J) C[6 * I + J] = Cm[I][J];
  info2a(inf, info, infod);
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

// detail::elastic at principal elastic strains: out[17] = p, q, D11, D12, D22, sig[3], ae[9] (row-major).
// Returns the EvalErr code (EE_ELASTIC_DOMAIN outside the HAR domain; out untouched then).
int ns_elastic(const double* d, const int* i, const double* eps3, double* out)
{
  detail::Elastic el;
  const int rc = detail::elastic(mk(d, i), eps3, el);
  if (rc) return rc;
  out[0] = el.p; out[1] = el.q; out[2] = el.D11; out[3] = el.D12; out[4] = el.D22;
  for (int a = 0; a < 3; ++a) out[5 + a] = el.sig[a];
  for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) out[8 + 3 * a + b] = el.ae[a][b];
  return 0;
}

// detail::floor_ev: out[2] = eps_v,f(eps_s), eps'_f.
void ns_floor_ev(const double* d, const int* i, double es, double* out)
{
  detail::floor_ev(mk(d, i), es, out[0], out[1]);
}

// detail::floor_project at principal strains: out[14] = eps_f[3], active, dfv, in_domain, epsp, Phi[9].
void ns_floor_project(const double* d, const int* i, const double* eps3, double* out)
{
  detail::FloorResult fr;
  detail::floor_project(mk(d, i), eps3, fr);
  for (int a = 0; a < 3; ++a) out[a] = fr.eps_f[a];
  out[3] = fr.active ? 1.0 : 0.0; out[4] = fr.dfv; out[5] = fr.in_domain ? 1.0 : 0.0; out[6] = fr.epsp;
  for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) out[7 + 3 * a + b] = fr.Phi[a][b];
}

// pRef(P), defaultPmin(P) (sheet §2.4, §9.7).
void ns_pref(const double* d, const int* i, double* out)
{
  const Params P = mk(d, i);
  out[0] = pRef(P);
  out[1] = defaultPmin(P);
}

// detail::energy_psi at principal strains (rc: EvalErr code).
int ns_energy_psi(const double* d, const int* i, const double* eps3, double* psi)
{
  return detail::energy_psi(mk(d, i), eps3, *psi);
}

}  // extern "C"
