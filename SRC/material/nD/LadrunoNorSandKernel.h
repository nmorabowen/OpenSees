/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 10/2026 (WP-144 P1a)
//
// LadrunoNorSandKernel -- the PURE numerical core of LadrunoNORSAND (NorSand in the
// Andrade & Borja 2006 three-invariant form). Header-only, std-only, NO OpenSees
// dependency (the kernel-oracle doctrine, as LadrunoJ2Kernel.h): it is tested
// standalone against the O2 algorithmic oracle
//     Ladruno_files/testbed/norsand_oracle/o2_algo/{kernel.py, api.py, params.py}
// through tests/ladrunonorsand_kernel_check.cpp and the ctypes parity harness
// Ladruno_files/testbed/norsand_oracle/kernel_parity/.
//
// CONTRACT. This file reproduces O2's ALGORITHM, branches and constants exactly
// (equation sheet Ladruno_implementation/144a_norsand_equation_sheet.md, sections
// 1-10, and plan 144 section 2.8 "Kernel contract learned at G1"). Every function in
// namespace detail mirrors the O2 function of the same name, expression by
// expression, so a parity diff is line-by-line traceable:
//   ww_theta, ga_theta, zeta_y          kernel.py §4   (S.8)-(S.11), corner branch (S.9)
//   elastic, invert_elastic             kernel.py §2   BA06 energy (S.3)-(S.5)
//   invariants                          kernel.py §3   (S.1),(S.6),(S.7), vertex rule §3.2
//   eta_of, pi_of_eta, yield_p          kernel.py §5   (S.12)-(S.13)
//   csl                                 kernel.py §6   (S.22), paper / fork CSL
//   pistar                              kernel.py §7   (S.23)-(S.24), B > 0 guard
//   cap_weight, omega_only, flow        kernel.py §5.3/§10 (S.14)-(S.21), cap (S.35)-(S.36)
//   pi_residual, solve_pi               kernel.py §8   nested pi_i solve (S.26)-(S.27),(S.37)
//   evaluate, scaled_norm, jacobian     kernel.py §9   (S.28)-(S.30)
//   atilde_ep                           kernel.py §9.3 (S.31)-(S.32)
//   return_map                          kernel.py §9.1 AB06 Box 2
//   spectral, tangent_small             kernel.py §9.4 (S.33)
//   step_once, step_ex                  api.py _step_once / step (substepping)
// Deliberate differences from O2 (none algorithmic; all at the API boundary):
//   * on a refused step() the returned state is the committed n (frozen) with
//     sigma = stress(n) and C = elasticTangent(n); O2 returns the trial-elastic state
//     of the whole increment. The refusal decision and its reason are O2's.
//   * where O2 would raise an uncaught Python exception (math.exp OverflowError on a
//     wild trial strain, a non-finite iterate) the kernel turns it into an
//     evaluation failure (detail::EE_NONFINITE): a refusal at the trial state, a
//     line-search backtrack at an iterate. Never reached on O2-defined paths.
//   * step() refuses a non-finite deps outright (no substepping), refusal
//     LOCAL_NOCONV, substeps = 0.
//   * initialState() also refuses pi_i0 >= 0 and a non-positive / non-finite v0
//     (O2 would refuse at the first step instead); validate() also refuses a
//     non-finite parameter (after O2's own checks, which run first and in O2's order).
//
// NUMERICAL CONTRACT (O2 kernel.py constants, reproduced exactly): RES_TOL 1e-12,
// MAX_LOCAL_ITERS 30, MAX_LINESEARCH 10, PI_TOL_REL 1e-12, MAX_PI_ITERS 50,
// PI_SCAN_REL 1e-3, PI_SCAN_MAX 1000, MAX_SUBSTEP_HALVINGS 8, R_TOL_REL 1e-8,
// CORNER_SIN3T 1e-8, F_TRIAL_TOL_REL 1e-10, EPS_S_TOL 1e-14, REPEATED_EIG_TOL 1e-10,
// PI_MAX_NEG 0. Bounded work: one step() is at most 2^9 - 1 = 511 backward-Euler
// solves, each at most 1 + 30 x 11 residual evaluations, each with a nested solve of
// at most 1000 scan + 50 safeguarded-Newton evaluations of r(pi_i); hitting any bound
// REFUSES (returns a nonzero refusal), it never force-accepts.
//
// Spectral decomposition: a cyclic Jacobi eigen-solver for symmetric 3x3 (Numerical
// Recipes rotation scheme). Tolerance: sweeps continue until every off-diagonal entry
// is exactly zero; from the 5th sweep on an off-diagonal a_pq is set to zero when
// 100|a_pq| is below the rounding unit of both |a_pp| and |a_qq| (|a_pq| < ~1e-18
// max|a_pp|), i.e. eigenvalues and eigenvectors to machine precision
// (||A V - V Lambda|| ~ 1e-16 ||A||). At most 50 sweeps (a 3x3 converges in <= ~6);
// non-convergence is a refusal. Eigenvalues are sorted ascending, as numpy.eigh.
//
// TENSOR CONVENTION. Compression negative throughout. Symmetric tensors are stored as
// 6 TENSOR components in the order {00,11,22,01,12,02} (true tensor shear, NOT
// engineering), as LadrunoJ2Kernel.h. The tangent C[6][6] is the true derivative
//     C[I][J] = d sigma_I / d eps_J,
// where eps_J for J >= 3 is the TENSOR shear component eps_kl varied symmetrically
// (eps_kl and eps_lk move together). In terms of the 4th-order tangent C4_ijkl:
//     C[I][J] = C4_ijkk             (J = kk normal),
//     C[I][J] = C4_ijkl + C4_ijlk   (J = kl shear; = 2 C4_ijkl when minor-symmetric).
// An OpenSees (engineering shear, gamma = 2 eps_kl) caller therefore uses
//     D_eng[I][J] = C[I][J]          (J normal),   D_eng[I][J] = 0.5 * C[I][J]  (J shear),
// and passes deps_tensor[J] = 0.5 * gamma_J for the shear slots.
// C is NON-SYMMETRIC (non-associated flow, the pi_i* and v terms, sheet §9.3):
// never use a symmetric solver with this material.

#ifndef LadrunoNorSandKernel_h
#define LadrunoNorSandKernel_h

#include <cmath>
#include <cstdio>
#include <string>

namespace ladruno_norsand {

// ---- public API (fixed; the wiring author codes to exactly this) ------------ //

// Parameter names follow sheet 144a §1.3. csl_mode 0 = paper (log), 1 = fork (power);
// zeta 0 = Willam-Warnke, 1 = Gudehus-Argyris; cap 0 = none, 1 = planar, 2 = smooth.
struct Params {
  double p0, kappa_hat, eps_v0, mu0, alpha0, M, N, N_bar, rho, rho_bar, chi, h,
         lambda_tilde, v_c0, e0, lambda_c, xi, p_a, c1, c2;
  int csl_mode; /*0 paper, 1 fork*/
  int zeta;     /*0 WW, 1 GA*/
  int cap;      /*0 none, 1 planar, 2 smooth*/
};

// eps_e: tensor components {00,11,22,01,12,02}. pi_i < 0 image pressure; v specific
// volume; v0 the INITIAL specific volume (a separate committed datum, plan §2.8);
// eps_p_v / eps_p_s accumulated plastic volumetric / deviatoric (sum dlam sqrt(2/3) Omega)
// strains; D_last the plastic dissipation of the last step (sum over its sub-steps).
struct State {
  double eps_e[6];
  double pi_i;
  double v;
  double v0;
  double eps_p_v;
  double eps_p_s;
  double D_last;
};

enum Refusal { OK = 0, LOCAL_NOCONV, LOCAL_LINESEARCH, PI_NOBRACKET, PI_NOCONV, B_NONPOS,
               P_OR_PI_NONNEG, NEGATIVE_DLAMBDA, SUBSTEPS_EXHAUSTED };

// refusal: OK, or SUBSTEPS_EXHAUSTED for every refusal that went down the substep ladder
// (by O2's contract every refused increment does: O2's reason suffix "(substeps
// exhausted at 2^8)"); the finest-level cause is available from detail::step_ex.
// plastic: any sub-step plastic; vertex / cap_active: the last sub-step's;
// local_iters / pi_iters: sums over the sub-steps; substeps: 1 = none, 2^k used,
// 256 on exhaustion. On a refusal the flags are those of the whole-increment attempt.
struct StepInfo {
  int refusal;
  int plastic;
  int vertex;
  int cap_active;
  int local_iters;
  int pi_iters;
  int substeps;
};

inline int validate(const Params& P, std::string& msg, bool& warn_rho_gt_rhobar);
inline int initialState(const Params& P, const double sigma0[6], double v0, double pi_i0,
                        State& out, std::string& msg);
inline int step(const Params& P, const State& n, const double deps[6], State& np1,
                double sigma[6], double C[6][6], StepInfo& info);
inline void stress(const Params& P, const State& s, double sigma[6]);
inline void elasticTangent(const Params& P, const State& s, double C[6][6]);

// ---- internals: O2 kernel.py / api.py, function by function ----------------- //
namespace detail {

// NUMERICAL CONTRACT (o2_algo/kernel.py, README)
constexpr double RES_TOL = 1.0e-12;
constexpr int    MAX_LOCAL_ITERS = 30;
constexpr int    MAX_LINESEARCH = 10;
constexpr double PI_TOL_REL = 1.0e-12;
constexpr int    MAX_PI_ITERS = 50;
constexpr double PI_SCAN_REL = 1.0e-3;
constexpr int    PI_SCAN_MAX = 1000;
constexpr int    MAX_SUBSTEP_HALVINGS = 8;
constexpr double R_TOL_REL = 1.0e-8;
constexpr double CORNER_SIN3T = 1.0e-8;
constexpr double F_TRIAL_TOL_REL = 1.0e-10;
constexpr double EPS_S_TOL = 1.0e-14;
constexpr double REPEATED_EIG_TOL = 1.0e-10;
constexpr double PI_MAX_NEG = 0.0;
constexpr int    JACOBI_MAX_SWEEPS = 50;
constexpr double PI_ = 3.14159265358979323846;   // == Python math.pi

inline double SQ23() { return std::sqrt(2.0 / 3.0); }
inline double SQ32() { return std::sqrt(1.5); }
inline double SQ6() { return std::sqrt(6.0); }

// Evaluation failures (O2's EvalError strings). Raised inside the residual evaluation
// only; turned into a line-search backtrack or a refusal by return_map.
enum EvalErr { EE_NONE = 0, EE_P_OR_PI_NONNEG, EE_PI_NONNEG, EE_B_NONPOS, EE_PI_FOLD,
               EE_PI_NOBRACKET, EE_PI_NOCONV, EE_NONFINITE, EE_SINGULAR_J };

inline const char* evalErrName(int e)
{
  switch (e) {
    case EE_NONE: return "";
    case EE_P_OR_PI_NONNEG: return "p_or_pi_nonneg";
    case EE_PI_NONNEG: return "pi_nonneg";
    case EE_B_NONPOS: return "B_nonpos";
    case EE_PI_FOLD: return "pi_fold";
    case EE_PI_NOBRACKET: return "pi_nobracket";
    case EE_PI_NOCONV: return "pi_noconv";
    case EE_NONFINITE: return "nonfinite";
    case EE_SINGULAR_J: return "singular_J";
    default: return "?";
  }
}

inline const char* refusalName(int r)
{
  switch (r) {
    case OK: return "OK";
    case LOCAL_NOCONV: return "LOCAL_NOCONV";
    case LOCAL_LINESEARCH: return "LOCAL_LINESEARCH";
    case PI_NOBRACKET: return "PI_NOBRACKET";
    case PI_NOCONV: return "PI_NOCONV";
    case B_NONPOS: return "B_NONPOS";
    case P_OR_PI_NONNEG: return "P_OR_PI_NONNEG";
    case NEGATIVE_DLAMBDA: return "NEGATIVE_DLAMBDA";
    case SUBSTEPS_EXHAUSTED: return "SUBSTEPS_EXHAUSTED";
    default: return "?";
  }
}

// O2 reason prefix "local_<err>" (failure at the first local evaluation) -> code.
inline int refusalOfLocalErr(int e)
{
  switch (e) {
    case EE_P_OR_PI_NONNEG:
    case EE_PI_NONNEG: return P_OR_PI_NONNEG;
    case EE_B_NONPOS: return B_NONPOS;
    case EE_PI_NOBRACKET: return PI_NOBRACKET;
    case EE_PI_NOCONV: return PI_NOCONV;
    case EE_PI_FOLD: return LOCAL_LINESEARCH;   // no own code; unreachable at dlam = 0
    default: return LOCAL_NOCONV;               // EE_NONFINITE, EE_SINGULAR_J
  }
}

inline double beta(const Params& P) { return (1.0 - P.N) / (1.0 - P.N_bar); }
inline double chi_bar(const Params& P) { return P.chi / beta(P); }

// --------------------------------------------------------------------------------
// small dense linear algebra
// --------------------------------------------------------------------------------
// LU with partial pivoting (first max-|.| pivot, as LAPACK idamax), solve A x = b for
// nrhs right-hand sides stored as columns of B (n <= 4). Returns false on an exactly
// zero (or non-finite) pivot, the case where LAPACK dgesv reports info > 0 and numpy
// raises LinAlgError.
inline bool lu_solve(int n, const double Ain[4][4], double B[4][4], int nrhs)
{
  double A[4][4];
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A[i][j] = Ain[i][j];
  for (int k = 0; k < n; ++k) {
    int piv = k;
    double amax = std::fabs(A[k][k]);
    for (int i = k + 1; i < n; ++i)
      if (std::fabs(A[i][k]) > amax) { amax = std::fabs(A[i][k]); piv = i; }
    if (!(amax > 0.0) || !std::isfinite(amax)) return false;
    if (piv != k) {
      for (int j = 0; j < n; ++j) { double t = A[k][j]; A[k][j] = A[piv][j]; A[piv][j] = t; }
      for (int j = 0; j < nrhs; ++j) { double t = B[k][j]; B[k][j] = B[piv][j]; B[piv][j] = t; }
    }
    for (int i = k + 1; i < n; ++i) {
      const double l = A[i][k] / A[k][k];
      A[i][k] = l;
      for (int j = k + 1; j < n; ++j) A[i][j] -= l * A[k][j];
      for (int j = 0; j < nrhs; ++j) B[i][j] -= l * B[k][j];
    }
  }
  for (int c = 0; c < nrhs; ++c)
    for (int i = n - 1; i >= 0; --i) {
      double s = B[i][c];
      for (int j = i + 1; j < n; ++j) s -= A[i][j] * B[j][c];
      B[i][c] = s / A[i][i];
    }
  return true;
}

// Cyclic Jacobi for a symmetric 3x3 (see the header comment for the tolerance).
// w ascending, V columns = eigenvectors. Returns false on non-convergence / non-finite.
inline void jrot(double a[3][3], double s, double tau, int i, int j, int k, int l)
{
  const double g = a[i][j];
  const double h = a[k][l];
  a[i][j] = g - s * (h + g * tau);
  a[k][l] = h + s * (g - h * tau);
}

inline bool eig_sym3(const double Ain[3][3], double w[3], double V[3][3])
{
  double a[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      a[i][j] = 0.5 * (Ain[i][j] + Ain[j][i]);       // numpy: eigh(0.5 (A + A^T))
      if (!std::isfinite(a[i][j])) return false;
    }
  double b[3], z[3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) V[i][j] = (i == j) ? 1.0 : 0.0;
    b[i] = w[i] = a[i][i];
    z[i] = 0.0;
  }
  bool done = false;
  for (int sweep = 1; sweep <= JACOBI_MAX_SWEEPS; ++sweep) {
    const double sm = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
    if (sm == 0.0) { done = true; break; }
    const double tresh = (sweep < 4) ? 0.2 * sm / 9.0 : 0.0;
    for (int ip = 0; ip < 2; ++ip) {
      for (int iq = ip + 1; iq < 3; ++iq) {
        const double g = 100.0 * std::fabs(a[ip][iq]);
        if (sweep > 4 && std::fabs(w[ip]) + g == std::fabs(w[ip])
                      && std::fabs(w[iq]) + g == std::fabs(w[iq])) {
          a[ip][iq] = 0.0;
        } else if (std::fabs(a[ip][iq]) > tresh) {
          double hh = w[iq] - w[ip];
          double t;
          if (std::fabs(hh) + g == std::fabs(hh)) {
            t = a[ip][iq] / hh;
          } else {
            const double th = 0.5 * hh / a[ip][iq];
            t = 1.0 / (std::fabs(th) + std::sqrt(1.0 + th * th));
            if (th < 0.0) t = -t;
          }
          const double c = 1.0 / std::sqrt(1.0 + t * t);
          const double s = t * c;
          const double tau = s / (1.0 + c);
          hh = t * a[ip][iq];
          z[ip] -= hh; z[iq] += hh;
          w[ip] -= hh; w[iq] += hh;
          a[ip][iq] = 0.0;
          for (int j = 0; j < ip; ++j) jrot(a, s, tau, j, ip, j, iq);
          for (int j = ip + 1; j < iq; ++j) jrot(a, s, tau, ip, j, j, iq);
          for (int j = iq + 1; j < 3; ++j) jrot(a, s, tau, ip, j, iq, j);
          for (int j = 0; j < 3; ++j) jrot(V, s, tau, j, ip, j, iq);
        }
      }
    }
    for (int i = 0; i < 3; ++i) { b[i] += z[i]; w[i] = b[i]; z[i] = 0.0; }
  }
  if (!done) return false;
  // ascending sort (insertion; stable), columns follow
  for (int i = 1; i < 3; ++i) {
    for (int k = i; k > 0 && w[k] < w[k - 1]; --k) {
      double t = w[k]; w[k] = w[k - 1]; w[k - 1] = t;
      for (int r = 0; r < 3; ++r) { t = V[r][k]; V[r][k] = V[r][k - 1]; V[r][k - 1] = t; }
    }
  }
  for (int i = 0; i < 3; ++i) if (!std::isfinite(w[i])) return false;
  return true;
}

// 6 tensor comps {00,11,22,01,12,02} <-> 3x3
inline void t6_to_m(const double t[6], double m[3][3])
{
  m[0][0] = t[0]; m[1][1] = t[1]; m[2][2] = t[2];
  m[0][1] = m[1][0] = t[3];
  m[1][2] = m[2][1] = t[4];
  m[0][2] = m[2][0] = t[5];
}

// api._from_principal: (V * vals) @ V.T
inline void from_principal(const double vals[3], const double V[3][3], double t[6])
{
  static const int I6[6] = {0, 1, 2, 0, 1, 0};
  static const int J6[6] = {0, 1, 2, 1, 2, 2};
  for (int k = 0; k < 6; ++k) {
    const int i = I6[k], j = J6[k];
    t[k] = (V[i][0] * vals[0]) * V[j][0] + (V[i][1] * vals[1]) * V[j][1]
         + (V[i][2] * vals[2]) * V[j][2];
  }
}

inline bool all_finite(const double* x, int n)
{
  for (int i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
  return true;
}

// --------------------------------------------------------------------------------
// §4 zeta(theta, rho): shape functions in the y-form of (S.8)-(S.9)
// --------------------------------------------------------------------------------
// Willam-Warnke (S.11): zeta, zeta'(theta), zeta''(theta) by the quotient rule.
inline void ww_theta(double theta, double rho, double& z, double& z1, double& z2)
{
  const double A = 4.0 * (1.0 - rho * rho);
  const double B = 2.0 * rho - 1.0;
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const double Nm = A * c * c + B * B;
  const double Nm1 = 2.0 * A * c;
  const double Nm2 = 2.0 * A;
  const double two = 2.0 * (1.0 - rho * rho);
  double Dn, Dn1, Dn2;
  if (std::fabs(B) < 1e-15) {     // rho = 1/2 exactly: unreachable through validated Params
    Dn = two * c; Dn1 = two; Dn2 = 0.0;
  } else {
    const double S = std::sqrt(A * c * c + 5.0 * rho * rho - 4.0 * rho);
    const double S1 = A * c / S;
    const double S2 = A / S - A * A * c * c / std::pow(S, 3.0);
    Dn = two * c + B * S;
    Dn1 = two + B * S1;
    Dn2 = B * S2;
  }
  z = Nm / Dn;
  const double num1 = Nm1 * Dn - Nm * Dn1;
  const double zc = num1 / std::pow(Dn, 2.0);
  const double zcc = (Nm2 * Dn - Nm * Dn2) / std::pow(Dn, 2.0) - 2.0 * num1 * Dn1 / std::pow(Dn, 3.0);
  z1 = -zc * s;                    // d zeta / d theta
  z2 = zcc * s * s - zc * c;       // d2 zeta / d theta2
}

// Gudehus-Argyris (S.10) in theta.
inline void ga_theta(double theta, double rho, double& z, double& z1, double& z2)
{
  z = ((1.0 + rho) + (1.0 - rho) * std::cos(3.0 * theta)) / (2.0 * rho);
  z1 = -3.0 * (1.0 - rho) * std::sin(3.0 * theta) / (2.0 * rho);
  z2 = -9.0 * (1.0 - rho) * std::cos(3.0 * theta) / (2.0 * rho);
}

// zeta, zeta_y, zeta_yy of (S.8)-(S.10) at theta. sin3theta, cos3theta are taken from
// the SAME theta as zeta', zeta'' (rule (i) of sheet §3.1). kind: 0 WW, 1 GA.
inline void zeta_y(double theta, double rho, int kind, double& z, double& zy, double& zyy)
{
  double z1, z2;
  if (kind == 1) {
    ga_theta(theta, rho, z, z1, z2);
    zy = SQ6() * (1.0 - rho) / (2.0 * rho);
    zyy = 0.0;
    return;
  }
  ww_theta(theta, rho, z, z1, z2);
  const double s3 = std::sin(3.0 * theta);
  if (std::fabs(s3) < CORNER_SIN3T) {                 // corner branch (S.9)
    double zc0, zc1, z2c;
    if (theta < PI_ / 6.0) {
      ww_theta(0.0, rho, zc0, zc1, z2c);
      zy = -SQ6() * z2c / 9.0;
      zyy = 0.0;
      return;
    }
    ww_theta(PI_ / 3.0, rho, zc0, zc1, z2c);
    zy = SQ6() * z2c / 9.0;
    zyy = 0.0;
    return;
  }
  const double c3 = std::cos(3.0 * theta);
  zy = -(2.0 / SQ6()) / s3 * z1;
  zyy = (2.0 / 3.0) / (s3 * s3) * (z2 - 3.0 * z1 * c3 / s3);
}

// --------------------------------------------------------------------------------
// §2 BA06 energy in principal elastic strains
// --------------------------------------------------------------------------------
struct Elastic {
  double eps_v, eps_s, nhat_e[3], p, q, D11, D12, D22, sig[3], ae[3][3];
};

inline void elastic(const Params& P, const double eps_e[3], Elastic& el)
{
  const double ev = (eps_e[0] + eps_e[1]) + eps_e[2];
  double e[3];
  for (int a = 0; a < 3; ++a) e[a] = eps_e[a] - ev / 3.0;
  const double ne = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  const double es = SQ23() * ne;
  double nh[3];
  for (int a = 0; a < 3; ++a) nh[a] = (es > EPS_S_TOL) ? e[a] / ne : 0.0;
  const double om = -(ev - P.eps_v0) / P.kappa_hat;
  const double E = std::exp(om);
  const double p = P.p0 * E * (1.0 + 1.5 * P.alpha0 / P.kappa_hat * es * es);
  const double q = 3.0 * (P.mu0 - P.alpha0 * P.p0 * E) * es;
  const double D11 = -p / P.kappa_hat;
  const double D22 = 3.0 * P.mu0 - 3.0 * P.alpha0 * P.p0 * E;
  const double D12 = 3.0 * P.p0 * P.alpha0 * es / P.kappa_hat * E;
  el.eps_v = ev; el.eps_s = es; el.p = p; el.q = q; el.D11 = D11; el.D12 = D12; el.D22 = D22;
  for (int a = 0; a < 3; ++a) {
    el.nhat_e[a] = nh[a];
    el.sig[a] = p + SQ23() * q * nh[a];
  }
  const double ratio = (es > EPS_S_TOL) ? q / es : D22;     // eps_s -> 0 limit (S.3 note)
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      const double t1 = D11;
      const double t2 = SQ23() * D12 * (nh[b] + nh[a]);
      const double t3 = (2.0 / 3.0) * D22 * (nh[a] * nh[b]);
      const double t4 = (2.0 * ratio / 3.0) * (((a == b ? 1.0 : 0.0) - 1.0 / 3.0) - nh[a] * nh[b]);
      el.ae[a][b] = ((t1 + t2) + t3) + t4;
    }
}

inline bool elastic_finite(const Elastic& el)
{
  if (!std::isfinite(el.p) || !std::isfinite(el.q)) return false;
  if (!all_finite(el.sig, 3)) return false;
  for (int a = 0; a < 3; ++a) if (!all_finite(el.ae[a], 3)) return false;
  return true;
}

// Principal elastic strains from principal stresses (Newton on (eps_v, eps_s) with the
// 2x2 Hessian; closed form when alpha0 = 0). Used by initialState only.
// Returns 0, or 1 if p >= 0, or 2 if the 2x2 Hessian is singular / non-finite.
inline int invert_elastic(const Params& P, const double sig[3], double out[3])
{
  const double p = ((sig[0] + sig[1]) + sig[2]) / 3.0;
  double xi[3];
  for (int a = 0; a < 3; ++a) xi[a] = sig[a] - p;
  const double R = std::sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
  const double q = SQ32() * R;
  double nh[3];
  for (int a = 0; a < 3; ++a) nh[a] = (R > 0.0) ? xi[a] / R : 0.0;
  if (!(p < 0.0)) return 1;
  double ev = P.eps_v0 - P.kappa_hat * std::log(p / P.p0);
  double es = q / (3.0 * P.mu0);
  if (P.alpha0 != 0.0) {
    for (int it = 0; it < 100; ++it) {
      double e3[3];
      for (int a = 0; a < 3; ++a) e3[a] = (es > 0.0) ? ev / 3.0 + SQ32() * es * nh[a] : ev / 3.0;
      Elastic el;
      elastic(P, e3, el);
      const double r0 = el.p - p, r1 = el.q - q;
      if (std::sqrt(r0 * r0 + r1 * r1) <= 1e-13 * std::fabs(p)) break;
      double H[4][4] = {{el.D11, el.D12, 0, 0}, {el.D12, el.D22, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
      double Bm[4][4] = {{r0, 0, 0, 0}, {r1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
      if (!lu_solve(2, H, Bm, 1)) return 2;
      ev -= Bm[0][0];
      es -= Bm[1][0];
    }
  }
  for (int a = 0; a < 3; ++a) out[a] = ev / 3.0 + SQ32() * es * nh[a];
  if (!all_finite(out, 3)) return 2;
  return 0;
}

// --------------------------------------------------------------------------------
// §3 invariants and their derivatives in principal space (S.1),(S.6),(S.7); vertex §3.2
// --------------------------------------------------------------------------------
struct Invariants {
  double p, q, R;
  bool vertex;
  double nhat[3], y, theta, y_a[3], nhat_ab[3][3], y_ab[3][3];
};

inline void invariants(const double sig[3], Invariants& I)
{
  I.p = ((sig[0] + sig[1]) + sig[2]) / 3.0;
  double xi[3];
  for (int a = 0; a < 3; ++a) xi[a] = sig[a] - I.p;
  I.R = std::sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
  if (I.R < R_TOL_REL * std::fabs(I.p)) {
    // vertex rule (§3.2): n_hat, y_a, n_hat_ab, y_ab := 0 ; q := 0 in F (F = p eta)
    I.q = 0.0; I.vertex = true; I.y = 0.0; I.theta = PI_ / 3.0;
    for (int a = 0; a < 3; ++a) {
      I.nhat[a] = 0.0; I.y_a[a] = 0.0;
      for (int b = 0; b < 3; ++b) { I.nhat_ab[a][b] = 0.0; I.y_ab[a][b] = 0.0; }
    }
    return;
  }
  const double R = I.R;
  I.vertex = false;
  I.q = SQ32() * R;
  for (int a = 0; a < 3; ++a) I.nhat[a] = xi[a] / R;
  const double S3 = (std::pow(xi[0], 3.0) + std::pow(xi[1], 3.0)) + std::pow(xi[2], 3.0);
  const double R2 = std::pow(R, 2.0), R3 = std::pow(R, 3.0), R5 = std::pow(R, 5.0);
  I.y = S3 / R3;
  const double arg = std::fmin(1.0, std::fmax(-1.0, SQ6() * I.y));
  I.theta = std::acos(arg) / 3.0;
  for (int a = 0; a < 3; ++a)
    I.y_a[a] = 3.0 * xi[a] * xi[a] / R3 - 3.0 * S3 * xi[a] / R5 - 1.0 / R;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      const double P1 = (a == b ? 1.0 : 0.0) - 1.0 / 3.0;
      I.nhat_ab[a][b] = (P1 - I.nhat[a] * I.nhat[b]) / R;
      const double T1 = 6.0 * (a == b ? xi[a] : 0.0) / R3;
      const double T2 = 3.0 * S3 / R5 * (P1 - 5.0 * (xi[a] * xi[b]) / R2);
      const double T3 = (xi[b] + xi[a]) / R3;
      const double T4 = 9.0 * (xi[a] * (xi[b] * xi[b]) + (xi[a] * xi[a]) * xi[b]) / R5;
      I.y_ab[a][b] = ((T1 - T2) + T3) - T4;
    }
}

// --------------------------------------------------------------------------------
// §5 yield function F, its p / pi_i derivatives (S.12)-(S.13)
// --------------------------------------------------------------------------------
inline double eta_of(const Params& P, double p, double pi)
{
  if (P.N == 0.0) return P.M * (1.0 + std::log(pi / p));
  return (P.M / P.N) * (1.0 - (1.0 - P.N) * std::pow(p / pi, P.N / (1.0 - P.N)));
}

// Inverse of (S.12): pi_i on the surface through (p, eta). Returns false if eta >= M/N.
inline bool pi_of_eta(const Params& P, double p, double eta, double& pi)
{
  if (P.N == 0.0) { pi = p * std::exp(eta / P.M - 1.0); return true; }
  const double d = 1.0 - eta * P.N / P.M;
  if (d <= 0.0) return false;
  pi = p * std::pow((1.0 - P.N) / d, (1.0 - P.N) / P.N);
  return true;
}

struct Yield { double eta, F_p, F_pp, F_pi, F_ppi, eta_p, eta_pi; };

inline int yield_p(const Params& P, double p, double pi, Yield& Y)
{
  if (p >= 0.0 || pi >= PI_MAX_NEG) return EE_P_OR_PI_NONNEG;
  const double N = P.N;
  const double eta = eta_of(P, p, pi);
  const double r = p / pi;
  const double F_p = (eta - P.M) / (1.0 - N);                 // = M ln(pi/p) for N = 0
  const double F_pp = -(P.M / (1.0 - N)) / p * std::pow(r, N / (1.0 - N));
  const double F_pi = P.M * std::pow(r, 1.0 / (1.0 - N));
  const double F_ppi = (P.M / ((1.0 - N) * p)) * std::pow(r, 1.0 / (1.0 - N));
  Y.eta = eta; Y.F_p = F_p; Y.F_pp = F_pp; Y.F_pi = F_pi; Y.F_ppi = F_ppi;
  Y.eta_p = (F_p - eta) / p;
  Y.eta_pi = F_pi / p;
  return EE_NONE;
}

// --------------------------------------------------------------------------------
// §6 CSL: psi_i(v, pi_i) and Lambda(pi_i) (S.22)
// --------------------------------------------------------------------------------
inline void csl(const Params& P, double v, double pi, double& psi, double& Lam)
{
  if (P.csl_mode == 0) {
    psi = v - P.v_c0 + P.lambda_tilde * std::log(-pi);
    Lam = P.lambda_tilde;
    return;
  }
  const double t = std::pow(-pi / P.p_a, P.xi);
  psi = (v - 1.0) - P.e0 + P.lambda_c * t;
  Lam = P.lambda_c * P.xi * t;
}

// --------------------------------------------------------------------------------
// §7 pi_i* = Pi(p, Omega, psi_i) and its derivatives (S.23)-(S.24)
// --------------------------------------------------------------------------------
inline int pistar(const Params& P, double p, double Om, double psi,
                  double& ps, double& Ppsi, double& POm)
{
  const double cb = chi_bar(P);
  const double a = SQ23() * cb;
  if (P.N == 0.0) {
    ps = p * std::exp(a * psi * Om / P.M);
    Ppsi = ps * a * Om / P.M;
    POm = ps * a * psi / P.M;
    return EE_NONE;
  }
  const double B = 1.0 - a * psi * Om * P.N / P.M;
  if (B <= 0.0) return EE_B_NONPOS;
  ps = p * std::pow(B, (P.N - 1.0) / P.N);
  const double den = P.M - a * psi * Om * P.N;
  Ppsi = a * Om * (1.0 - P.N) * ps / den;
  POm = a * psi * (1.0 - P.N) * ps / den;
  return EE_NONE;
}

// --------------------------------------------------------------------------------
// §5.3 + §10 flow vector, Hessian, Omega (S.17)-(S.21), cap (S.35)-(S.36)
// --------------------------------------------------------------------------------
struct Flow {
  double F, f_a[3], q_a[3], q_ab[3][3], q_api[3], Om, Om_a[3], Om_pi, w, zeta, zeta_bar;
  Yield Y;
};

// w(eta), dw/deta of (S.35). none: w = 1; planar: step at eta_1 with w = 1 at equality.
inline void cap_weight(const Params& P, double eta, double& w, double& w_eta)
{
  if (P.cap == 0) { w = 1.0; w_eta = 0.0; return; }
  const double e1 = P.c1 * P.M;
  const double e2 = P.c2 * P.M;
  if (P.cap == 1) { w = (eta >= e1) ? 1.0 : 0.0; w_eta = 0.0; return; }
  const double t = (eta - e1) / (e2 - e1);
  if (t <= 0.0) { w = 0.0; w_eta = 0.0; return; }
  if (t >= 1.0) { w = 1.0; w_eta = 0.0; return; }
  const double S = std::pow(t, 3.0) * (10.0 - 15.0 * t + 6.0 * t * t);
  const double Sp = 30.0 * t * t * std::pow(1.0 - t, 2.0);
  w = S;
  w_eta = Sp / (e2 - e1);
}

// Omega (capped), Omega_pi, Y -- the pieces the nested pi_i loop needs (cheap path).
inline int omega_only(const Params& P, const Invariants& inv, double pi,
                      double& Om, double& Om_pi, Yield& Y)
{
  const int e = yield_p(P, inv.p, pi, Y);
  if (e) return e;
  if (inv.vertex) { Om = 0.0; Om_pi = 0.0; return EE_NONE; }
  double zb, zby, zbyy;
  zeta_y(inv.theta, P.rho_bar, P.zeta, zb, zby, zbyy);
  const double Sy = inv.y_a[0] * inv.y_a[0] + inv.y_a[1] * inv.y_a[1] + inv.y_a[2] * inv.y_a[2];
  const double Omu = std::sqrt(1.5 * zb * zb + std::pow(zby * inv.q, 2.0) * Sy);
  double w, w_eta;
  cap_weight(P, Y.eta, w, w_eta);
  Om = w * Omu;
  Om_pi = w_eta * Y.eta_pi * Omu;
  return EE_NONE;
}

inline int flow(const Params& P, const Invariants& inv, double pi, Flow& fl)
{
  Yield Y;
  const int e = yield_p(P, inv.p, pi, Y);
  if (e) return e;
  const double bt = beta(P);
  double F, f_a[3], qu_a[3], qu_ab[3][3], qu_api[3], Omu, Omu_a[3], z, zb;
  if (inv.vertex) {
    // §3.2: purely volumetric flow; F = p eta; Omega = 0
    F = inv.p * Y.eta;
    for (int a = 0; a < 3; ++a) {
      f_a[a] = Y.F_p / 3.0;
      qu_a[a] = bt * Y.F_p / 3.0;
      qu_api[a] = bt * Y.F_ppi / 3.0;
      Omu_a[a] = 0.0;
      for (int b = 0; b < 3; ++b) qu_ab[a][b] = bt * Y.F_pp / 9.0;
    }
    Omu = 0.0;
    z = zb = 1.0;
  } else {
    double zy, zyy_unused, zby, zbyy;
    zeta_y(inv.theta, P.rho, P.zeta, z, zy, zyy_unused);
    zeta_y(inv.theta, P.rho_bar, P.zeta, zb, zby, zbyy);
    const double q = inv.q;
    const double* nh = inv.nhat;
    const double* y_a = inv.y_a;
    F = z * q + inv.p * Y.eta;
    for (int a = 0; a < 3; ++a) {
      f_a[a] = Y.F_p / 3.0 + SQ32() * z * nh[a] + zy * q * y_a[a];                  // (S.14)
      qu_a[a] = bt * Y.F_p / 3.0 + SQ32() * zb * nh[a] + zby * q * y_a[a];           // (S.17)
      qu_api[a] = bt * Y.F_ppi / 3.0;                                                // (S.19)
    }
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)                                                    // (S.18)
        qu_ab[a][b] = bt * Y.F_pp / 9.0 + SQ32() * zb * inv.nhat_ab[a][b]
                    + zby * q * inv.y_ab[a][b] + zbyy * q * (y_a[a] * y_a[b])
                    + SQ32() * zby * (nh[a] * y_a[b] + y_a[a] * nh[b]);
    const double Sy = y_a[0] * y_a[0] + y_a[1] * y_a[1] + y_a[2] * y_a[2];
    Omu = std::sqrt(1.5 * zb * zb + std::pow(zby * q, 2.0) * Sy);                    // (S.20)
    for (int a = 0; a < 3; ++a) {                                                    // (S.21)
      const double yaby = (inv.y_ab[a][0] * y_a[0] + inv.y_ab[a][1] * y_a[1]) + inv.y_ab[a][2] * y_a[2];
      Omu_a[a] = (1.5 * zb * zby * y_a[a]
                  + zby * q * (zbyy * q * y_a[a] + SQ32() * zby * nh[a]) * Sy
                  + std::pow(zby * q, 2.0) * yaby) / Omu;
    }
  }
  double w, w_eta;
  cap_weight(P, Y.eta, w, w_eta);
  fl.F = F; fl.w = w; fl.zeta = z; fl.zeta_bar = zb; fl.Y = Y;
  for (int a = 0; a < 3; ++a) fl.f_a[a] = f_a[a];
  if (w == 1.0 && w_eta == 0.0) {
    for (int a = 0; a < 3; ++a) {
      fl.q_a[a] = qu_a[a]; fl.q_api[a] = qu_api[a]; fl.Om_a[a] = Omu_a[a];
      for (int b = 0; b < 3; ++b) fl.q_ab[a][b] = qu_ab[a][b];
    }
    fl.Om = Omu;
    fl.Om_pi = 0.0;
    return EE_NONE;
  }
  double g_a[3];
  for (int a = 0; a < 3; ++a) g_a[a] = qu_a[a] + 1.0 / 3.0;
  for (int a = 0; a < 3; ++a) {
    fl.q_a[a] = -1.0 / 3.0 + w * g_a[a];
    for (int b = 0; b < 3; ++b)                                                      // (S.36)
      fl.q_ab[a][b] = w * qu_ab[a][b] + w_eta * Y.eta_p / 3.0 * g_a[a];
    fl.q_api[a] = w * qu_api[a] + w_eta * Y.eta_pi * g_a[a];
    fl.Om_a[a] = w * Omu_a[a] + w_eta * Y.eta_p / 3.0 * Omu;
  }
  fl.Om = w * Omu;
  fl.Om_pi = w_eta * Y.eta_pi * Omu;
  return EE_NONE;
}

// --------------------------------------------------------------------------------
// §8 nested scalar solve for pi_i (S.26)-(S.27), (S.37)
// --------------------------------------------------------------------------------
inline int pi_residual(const Params& P, const Invariants& inv, double dlam, double v,
                       double pi_n, double pi, double& r, double& rp)
{
  if (pi >= PI_MAX_NEG) return EE_PI_NONNEG;
  const double k = SQ23() * P.h;
  double Om, Om_pi;
  Yield Y;
  int e = omega_only(P, inv, pi, Om, Om_pi, Y);
  if (e) return e;
  double psi, Lam;
  csl(P, v, pi, psi, Lam);
  double ps, Ppsi, POm;
  e = pistar(P, inv.p, Om, psi, ps, Ppsi, POm);
  if (e) return e;
  r = pi - pi_n - k * dlam * (ps - pi) * Om;
  rp = 1.0 - k * dlam * ((POm * Om_pi + Ppsi * Lam / pi - 1.0) * Om + (ps - pi) * Om_pi);  // (S.37)
  return EE_NONE;
}

// Nested scalar solve of r(pi_i) = 0, the root CONTINUOUS WITH pi_i,n (O2 solve_pi):
// 1. |r(pi_i,n)| <= tol -> pi_i,n, 0 iterations;
// 2. scan from pi_i,n in fixed steps PI_SCAN_REL*|pi_i,n| (away from 0 if r > 0, toward
//    0 if r < 0) to the FIRST sign change; |r| growing first -> pi_fold; PI_SCAN_MAX
//    steps without a sign change -> pi_nobracket;
// 3. safeguarded Newton/bisection inside that bracket, MAX_PI_ITERS -> pi_noconv.
inline int solve_pi(const Params& P, const Invariants& inv, double dlam, double v, double pi_n,
                    double& pi_out, double& c_out, int& iters)
{
  const double tol = PI_TOL_REL * std::fabs(pi_n);
  double a = pi_n, ra, rpa;
  int e = pi_residual(P, inv, dlam, v, pi_n, a, ra, rpa);
  if (e) return e;
  if (std::fabs(ra) <= tol) { pi_out = a; c_out = rpa; iters = 0; return EE_NONE; }
  // 2. scan for the first bracket
  const double d = (ra > 0.0) ? -1.0 : 1.0;          // pi_i < 0: -1 = away from zero
  const double hstep = PI_SCAN_REL * std::fabs(pi_n);
  int it = 0;
  double b = 0.0, rb = 0.0, rpb = 0.0;
  bool bracketed = false;
  for (int s = 0; s < PI_SCAN_MAX; ++s) {
    b = a + d * hstep;
    e = pi_residual(P, inv, dlam, v, pi_n, b, rb, rpb);
    if (e) return e;
    it += 1;
    if (std::fabs(rb) <= tol) { pi_out = b; c_out = rpb; iters = it; return EE_NONE; }
    if (rb * ra < 0.0) { bracketed = true; break; }
    if (std::fabs(rb) > std::fabs(ra)) return EE_PI_FOLD;
    a = b; ra = rb; rpa = rpb;
  }
  if (!bracketed) return EE_PI_NOBRACKET;
  // 3. safeguarded Newton inside [lo, hi]
  double lo, hi, rlo;
  if (a < b) { lo = a; hi = b; rlo = ra; } else { lo = b; hi = a; rlo = rb; }
  double x, r, rp;
  if (std::fabs(ra) < std::fabs(rb)) { x = a; r = ra; rp = rpa; }
  else { x = b; r = rb; rp = rpb; }
  for (int k = 0; k < MAX_PI_ITERS; ++k) {
    double xn = (rp != 0.0) ? x - r / rp : lo;
    if (!(lo < xn && xn < hi)) xn = 0.5 * (lo + hi);
    e = pi_residual(P, inv, dlam, v, pi_n, xn, r, rp);
    if (e) return e;
    x = xn;
    it += 1;
    if (std::fabs(r) <= tol) { pi_out = x; c_out = rp; iters = it; return EE_NONE; }
    if (r * rlo > 0.0) { lo = x; rlo = r; }
    else { hi = x; }
  }
  return EE_PI_NOCONV;
}

// --------------------------------------------------------------------------------
// §9 residual, Jacobian, tangent
// --------------------------------------------------------------------------------
struct PointEval {
  double eps_e[3];
  double dlam;
  Elastic el;
  Invariants inv;
  Flow fl;
  double pi, c;
  int pi_iters;
  double psi, Lam, ps, Ppsi, POm;
  double r[4];
  double P_a[3], Pi_b[3], Pi_lam, Pi_v;   // sensitivities of the converged pi_i (S.28)
};

inline int evaluate(const Params& P, const double eps_e[3], double dlam, const double eps_tr[3],
                    double v, double pi_n, PointEval& pe)
{
  for (int a = 0; a < 3; ++a) pe.eps_e[a] = eps_e[a];
  pe.dlam = dlam;
  if (!all_finite(eps_e, 3) || !std::isfinite(dlam)) return EE_NONFINITE;
  elastic(P, eps_e, pe.el);
  if (!elastic_finite(pe.el)) return EE_NONFINITE;
  invariants(pe.el.sig, pe.inv);
  int e = solve_pi(P, pe.inv, dlam, v, pi_n, pe.pi, pe.c, pe.pi_iters);
  if (e) return e;
  e = flow(P, pe.inv, pe.pi, pe.fl);
  if (e) return e;
  csl(P, v, pe.pi, pe.psi, pe.Lam);
  e = pistar(P, pe.inv.p, pe.fl.Om, pe.psi, pe.ps, pe.Ppsi, pe.POm);
  if (e) return e;
  for (int a = 0; a < 3; ++a) pe.r[a] = eps_e[a] - eps_tr[a] + dlam * pe.fl.q_a[a];
  pe.r[3] = pe.fl.F;
  const double k = SQ23() * P.h;
  const double p = pe.inv.p, ps = pe.ps, pi = pe.pi, c = pe.c, Om = pe.fl.Om;
  for (int a = 0; a < 3; ++a) {
    const double pistar_a = (ps / (3.0 * p)) + pe.POm * pe.fl.Om_a[a];                 // (S.24)
    pe.P_a[a] = (k * dlam / c) * (Om * pistar_a + (ps - pi) * pe.fl.Om_a[a]);          // (S.28)/(S.37)
  }
  for (int b = 0; b < 3; ++b)
    pe.Pi_b[b] = (pe.P_a[0] * pe.el.ae[0][b] + pe.P_a[1] * pe.el.ae[1][b]) + pe.P_a[2] * pe.el.ae[2][b];
  pe.Pi_lam = (k * Om / c) * (ps - pi);
  pe.Pi_v = (k * dlam * Om / c) * pe.Ppsi;
  if (!all_finite(pe.r, 4)) return EE_NONFINITE;
  return EE_NONE;
}

inline double scaled_norm(const Params& P, const double r[4])
{
  return std::sqrt(std::pow(r[0], 2.0) + std::pow(r[1], 2.0) + std::pow(r[2], 2.0)
                   + std::pow(r[3] / std::fabs(P.p0), 2.0));
}

// (S.30)
inline void jacobian(const PointEval& pe, double J[4][4])
{
  const Flow& fl = pe.fl;
  const Elastic& el = pe.el;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      const double qae = (fl.q_ab[a][0] * el.ae[0][b] + fl.q_ab[a][1] * el.ae[1][b])
                       + fl.q_ab[a][2] * el.ae[2][b];
      J[a][b] = (a == b ? 1.0 : 0.0) + pe.dlam * (qae + fl.q_api[a] * pe.Pi_b[b]);
    }
    J[a][3] = fl.q_a[a] + pe.dlam * fl.q_api[a] * pe.Pi_lam;
  }
  for (int b = 0; b < 3; ++b) {
    const double fae = (fl.f_a[0] * el.ae[0][b] + fl.f_a[1] * el.ae[1][b]) + fl.f_a[2] * el.ae[2][b];
    J[3][b] = fae + fl.Y.F_pi * pe.Pi_b[b];
  }
  J[3][3] = fl.Y.F_pi * pe.Pi_lam;
}

// (S.31)-(S.32): a~^ep_ab = d sigma_a / d eps~_b. vfac = v0 (small strain).
inline bool atilde_ep(const PointEval& pe, const double J[4][4], double vfac, double at[3][3])
{
  const Flow& fl = pe.fl;
  double s[4];
  for (int a = 0; a < 3; ++a) s[a] = pe.dlam * fl.q_api[a] * pe.Pi_v * vfac;
  s[3] = fl.Y.F_pi * pe.Pi_v * vfac;
  double b[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  if (!lu_solve(4, J, b, 4)) return false;               // b = J^{-1}
  double bs[4];
  for (int i = 0; i < 4; ++i)
    bs[i] = ((b[i][0] * s[0] + b[i][1] * s[1]) + b[i][2] * s[2]) + b[i][3] * s[3];
  double dxde[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) dxde[i][j] = b[i][j] - bs[i];
  for (int a = 0; a < 3; ++a)
    for (int j = 0; j < 3; ++j)
      at[a][j] = (pe.el.ae[a][0] * dxde[0][j] + pe.el.ae[a][1] * dxde[1][j]) + pe.el.ae[a][2] * dxde[2][j];
  for (int a = 0; a < 3; ++a) if (!all_finite(at[a], 3)) return false;
  return true;
}

// O2 StepResult (the fields the kernel needs). reason: Refusal code of the refusal at
// THIS backward-Euler solve (OK if accepted); sub: the EvalErr behind it (O2's
// "trial_<e>", "local_<e>", "local_linesearch:<e>", "local_singular_J").
struct ReturnResult {
  double eps_e[3], sig[3], pi, dlam, q_a[3], Om, D, atilde[3][3];
  bool plastic, vertex, cap_active, refused;
  int reason, sub;
  int local_iters, pi_iters;
};

inline void rr_refuse(ReturnResult& R, const double eps_tr[3], const Elastic& el, double pi_n,
                      bool plastic, bool vertex, int reason, int sub, int it, int pit)
{
  for (int a = 0; a < 3; ++a) {
    R.eps_e[a] = eps_tr[a]; R.sig[a] = el.sig[a]; R.q_a[a] = 0.0;
    for (int b = 0; b < 3; ++b) R.atilde[a][b] = el.ae[a][b];
  }
  R.pi = pi_n; R.dlam = 0.0; R.Om = 0.0; R.D = 0.0;
  R.plastic = plastic; R.vertex = vertex; R.cap_active = false; R.refused = true;
  R.reason = reason; R.sub = sub; R.local_iters = it; R.pi_iters = pit;
}

// One backward-Euler step in principal space (sheet §9.1, O2 return_map). Never throws;
// refusals come back in R.refused / R.reason / R.sub.
inline void return_map(const Params& P, const double eps_tr[3], double pi_n, double v, double vfac,
                       ReturnResult& R)
{
  // 1-2. trial
  Elastic el;
  elastic(P, eps_tr, el);
  Invariants inv;
  if (!elastic_finite(el)) {
    // O2: math.exp OverflowError (uncaught). Kernel: refuse at the trial state.
    for (int a = 0; a < 3; ++a) { el.sig[a] = 0.0; for (int b = 0; b < 3; ++b) el.ae[a][b] = 0.0; }
    rr_refuse(R, eps_tr, el, pi_n, false, false, LOCAL_NOCONV, EE_NONFINITE, 0, 0);
    return;
  }
  invariants(el.sig, inv);
  Flow fl0;
  int e = flow(P, inv, pi_n, fl0);
  if (e) {   // "trial_<e>"
    rr_refuse(R, eps_tr, el, pi_n, false, inv.vertex,
              (e == EE_NONFINITE ? LOCAL_NOCONV : P_OR_PI_NONNEG), e, 0, 0);
    return;
  }
  double psi0, Lam0;
  csl(P, v, pi_n, psi0, Lam0);
  (void)psi0; (void)Lam0;
  if (fl0.F <= F_TRIAL_TOL_REL * std::fabs(P.p0)) {
    for (int a = 0; a < 3; ++a) {
      R.eps_e[a] = eps_tr[a]; R.sig[a] = el.sig[a]; R.q_a[a] = 0.0;
      for (int b = 0; b < 3; ++b) R.atilde[a][b] = el.ae[a][b];
    }
    R.pi = pi_n; R.dlam = 0.0; R.Om = 0.0; R.D = 0.0;
    R.plastic = false; R.vertex = inv.vertex; R.cap_active = fl0.w < 1.0; R.refused = false;
    R.reason = OK; R.sub = EE_NONE; R.local_iters = 0; R.pi_iters = 0;
    return;
  }
  // 3-4. local Newton on x = (eps_e, dlam)
  double x[4] = {eps_tr[0], eps_tr[1], eps_tr[2], 0.0};
  int pi_total = 0;
  PointEval pe, pen;
  e = evaluate(P, x, x[3], eps_tr, v, pi_n, pe);
  if (e) {   // "local_<e>"
    rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, refusalOfLocalErr(e), e, 0, 0);
    return;
  }
  double rn = scaled_norm(P, pe.r);
  pi_total += pe.pi_iters;
  int it = 0;
  bool converged = rn <= RES_TOL;
  while (!converged && it < MAX_LOCAL_ITERS) {
    double J[4][4];
    jacobian(pe, J);
    double dxB[4][4] = {{pe.r[0], 0, 0, 0}, {pe.r[1], 0, 0, 0}, {pe.r[2], 0, 0, 0}, {pe.r[3], 0, 0, 0}};
    if (!lu_solve(4, J, dxB, 1)) {   // "local_singular_J"
      rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, LOCAL_NOCONV, EE_SINGULAR_J, it, pi_total);
      return;
    }
    const double dx[4] = {dxB[0][0], dxB[1][0], dxB[2][0], dxB[3][0]};
    double alpha = 1.0;
    bool accepted = false;
    int last_err = EE_NONE;
    double xn[4];
    double rnn = 0.0;
    for (int ls = 0; ls < MAX_LINESEARCH + 1; ++ls) {
      for (int i = 0; i < 4; ++i) xn[i] = x[i] - alpha * dx[i];
      const int ee = evaluate(P, xn, xn[3], eps_tr, v, pi_n, pen);
      if (ee == EE_NONE) {
        rnn = scaled_norm(P, pen.r);
        if (rnn < rn || rnn <= RES_TOL) { accepted = true; break; }
      } else {
        // a nested failure or an inadmissible iterate: rejected, the step (dlam
        // included) is backtracked exactly like a non-decreasing residual.
        last_err = ee;
      }
      alpha *= 0.5;
    }
    if (!accepted) {   // "local_linesearch[:<e>]"
      rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, LOCAL_LINESEARCH, last_err, it, pi_total);
      return;
    }
    for (int i = 0; i < 4; ++i) x[i] = xn[i];
    pe = pen;
    rn = rnn;
    pi_total += pe.pi_iters;
    it += 1;
    converged = rn <= RES_TOL;
  }
  if (!converged) {
    rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, LOCAL_NOCONV, EE_NONE, it, pi_total);
    return;
  }
  if (pe.dlam < 0.0) {
    rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, NEGATIVE_DLAMBDA, EE_NONE, it, pi_total);
    return;
  }
  // 5. tangent and diagnostics
  double J[4][4];
  jacobian(pe, J);
  if (!atilde_ep(pe, J, vfac, R.atilde)) {   // O2: LinAlgError from np.linalg.inv (uncaught)
    rr_refuse(R, eps_tr, el, pi_n, true, inv.vertex, LOCAL_NOCONV, EE_SINGULAR_J, it, pi_total);
    return;
  }
  for (int a = 0; a < 3; ++a) {
    R.eps_e[a] = pe.eps_e[a]; R.sig[a] = pe.el.sig[a]; R.q_a[a] = pe.fl.q_a[a];
  }
  R.pi = pe.pi; R.dlam = pe.dlam; R.Om = pe.fl.Om;
  R.D = pe.dlam * ((pe.el.sig[0] * pe.fl.q_a[0] + pe.el.sig[1] * pe.fl.q_a[1]) + pe.el.sig[2] * pe.fl.q_a[2]);
  R.plastic = true; R.vertex = pe.inv.vertex; R.cap_active = pe.fl.w < 1.0; R.refused = false;
  R.reason = OK; R.sub = EE_NONE; R.local_iters = it; R.pi_iters = pi_total;
}

// --------------------------------------------------------------------------------
// §9.4 spectral assembly of the 4th-order tangent (S.33), compressed to 6x6
// --------------------------------------------------------------------------------
// C4 = sum_ab diag_ab m^a (x) m^b + half sum_{a!=b} spin_ab (m^ab (x) m^ab + m^ab (x) m^ba)
inline void spectral(const double V[3][3], const double diag_ab[3][3], const double spin_ab[3][3],
                     double half, double C4[3][3][3][3])
{
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) C4[i][j][k][l] = 0.0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) {
          C4[i][j][k][l] += diag_ab[a][b] * ((V[i][a] * V[j][a]) * (V[k][b] * V[l][b]));
          if (a != b) {
            const double mab_ij = V[i][a] * V[j][b];
            C4[i][j][k][l] += half * spin_ab[a][b] * (mab_ij * (V[k][a] * V[l][b])
                                                      + mab_ij * (V[k][b] * V[l][a]));
          }
        }
}

// C[I][J] = d sigma_I / d eps_J with symmetric variation of the shear slots (header).
inline void compress_c4(const double C4[3][3][3][3], double C[6][6])
{
  static const int I6[6] = {0, 1, 2, 0, 1, 0};
  static const int J6[6] = {0, 1, 2, 1, 2, 2};
  for (int I = 0; I < 6; ++I)
    for (int Jc = 0; Jc < 6; ++Jc) {
      const int i = I6[I], j = J6[I], k = I6[Jc], l = J6[Jc];
      C[I][Jc] = (Jc < 3) ? C4[i][j][k][k] : C4[i][j][k][l] + C4[i][j][l][k];
    }
}

// (S.33): small-strain consistent tangent from a~^ep, converged sigma_a, trial eps~_a.
inline void tangent_small(const double atilde[3][3], const double sig[3], const double eps_tr[3],
                          const double V[3][3], double C[6][6])
{
  double g[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      if (a == b) continue;
      const double d = eps_tr[a] - eps_tr[b];
      if (std::fabs(d) < REPEATED_EIG_TOL) g[a][b] = atilde[a][a] - atilde[a][b];
      else g[a][b] = (sig[a] - sig[b]) / d;
    }
  double C4[3][3][3][3];
  spectral(V, atilde, g, 0.5, C4);
  compress_c4(C4, C);
}

// --------------------------------------------------------------------------------
// api.py: _step_once and step (substepping)
// --------------------------------------------------------------------------------
struct OnceOut {
  State st;
  double sigma[6];
  ReturnResult res;
  double w[3], V[3][3];
};

// One backward-Euler increment without substepping (O2 api._step_once, small strain).
inline void step_once(const Params& P, const State& n, const double deps[6], OnceOut& o)
{
  const double tr = (deps[0] + deps[1]) + deps[2];
  const double v = n.v + n.v0 * tr;
  const double vfac = n.v0;
  double et6[6], et[3][3];
  for (int i = 0; i < 6; ++i) et6[i] = n.eps_e[i] + deps[i];
  t6_to_m(et6, et);
  o.st = n;
  o.st.v = v;
  if (!eig_sym3(et, o.w, o.V)) {
    for (int a = 0; a < 3; ++a) { o.w[a] = 0.0; for (int b = 0; b < 3; ++b) o.V[a][b] = (a == b); }
    Elastic el;
    for (int a = 0; a < 3; ++a) { el.sig[a] = 0.0; for (int b = 0; b < 3; ++b) el.ae[a][b] = 0.0; }
    rr_refuse(o.res, o.w, el, n.pi_i, false, false, LOCAL_NOCONV, EE_NONFINITE, 0, 0);
    for (int i = 0; i < 6; ++i) o.sigma[i] = 0.0;
    return;
  }
  return_map(P, o.w, n.pi_i, v, vfac, o.res);
  from_principal(o.res.sig, o.V, o.sigma);
  from_principal(o.res.eps_e, o.V, o.st.eps_e);
  o.st.pi_i = o.res.pi;
  o.st.D_last = 0.0;
  if (o.res.refused) return;
  o.st.D_last = o.res.D;
  if (o.res.plastic) {
    o.st.eps_p_v = n.eps_p_v + o.res.dlam * ((o.res.q_a[0] + o.res.q_a[1]) + o.res.q_a[2]);
    o.st.eps_p_s = n.eps_p_s + o.res.dlam * SQ23() * o.res.Om;
  }
}

// O2 api.step: the increment is attempted whole; a refused increment is retried as
// 2, 4, ..., 2^MAX_SUBSTEP_HALVINGS equal sub-increments, each a full BE step chained on
// the previous one. The tangent of a substepped increment is the CTO of the LAST
// sub-increment. finest / finest_sub: the reason of the refusal at the finest level
// (OK / EE_NONE on success).
inline int step_ex(const Params& P, const State& n, const double deps[6], State& np1,
                   double sigma[6], double C[6][6], StepInfo& info, int& finest, int& finest_sub)
{
  finest = OK;
  finest_sub = EE_NONE;
  if (!all_finite(deps, 6)) {
    np1 = n;
    stress(P, n, sigma);
    elasticTangent(P, n, C);
    info.refusal = LOCAL_NOCONV; info.plastic = 0; info.vertex = 0; info.cap_active = 0;
    info.local_iters = 0; info.pi_iters = 0; info.substeps = 0;
    finest = LOCAL_NOCONV; finest_sub = EE_NONFINITE;
    return info.refusal;
  }
  StepInfo first = {0, 0, 0, 0, 0, 0, 0};
  int last_reason = OK, last_sub = EE_NONE;
  OnceOut cur, nxt;
  for (int j = 0; j <= MAX_SUBSTEP_HALVINGS; ++j) {
    const int m = 1 << j;
    double dsub[6];
    for (int i = 0; i < 6; ++i) dsub[i] = deps[i] / static_cast<double>(m);
    cur.st = n;
    double D = 0.0;
    int iters = 0, piters = 0;
    bool plastic = false, ok = true;
    for (int k = 0; k < m; ++k) {
      step_once(P, cur.st, dsub, nxt);
      cur = nxt;
      if (cur.res.refused) { ok = false; break; }
      D += cur.st.D_last;
      iters += cur.res.local_iters;
      piters += cur.res.pi_iters;
      plastic = plastic || cur.res.plastic;
    }
    if (ok) {
      np1 = cur.st;
      np1.D_last = D;
      for (int i = 0; i < 6; ++i) sigma[i] = cur.sigma[i];
      tangent_small(cur.res.atilde, cur.res.sig, cur.w, cur.V, C);
      info.refusal = OK; info.plastic = plastic ? 1 : 0;
      info.vertex = cur.res.vertex ? 1 : 0; info.cap_active = cur.res.cap_active ? 1 : 0;
      info.local_iters = iters; info.pi_iters = piters; info.substeps = m;
      return OK;
    }
    if (j == 0) {
      first.plastic = cur.res.plastic ? 1 : 0; first.vertex = cur.res.vertex ? 1 : 0;
      first.cap_active = cur.res.cap_active ? 1 : 0;
      first.local_iters = cur.res.local_iters; first.pi_iters = cur.res.pi_iters;
    }
    last_reason = cur.res.reason;
    last_sub = cur.res.sub;
  }
  // exhausted: frozen at n (API contract), O2's flags of the whole-increment attempt
  np1 = n;
  stress(P, n, sigma);
  elasticTangent(P, n, C);
  info = first;
  info.refusal = SUBSTEPS_EXHAUSTED;
  info.substeps = 1 << MAX_SUBSTEP_HALVINGS;
  finest = last_reason;
  finest_sub = last_sub;
  return info.refusal;
}

}  // namespace detail

// ---- public API implementation --------------------------------------------- //

// O2 params.validate(): 0 ok, nonzero = refusal code (msg explains). On success msg
// holds the warnings (rho > rho_bar; chi > 0), if any.
inline int validate(const Params& P, std::string& msg, bool& warn_rho_gt_rhobar)
{
  char buf[512];
  msg.clear();
  warn_rho_gt_rhobar = false;
  if (P.zeta != 0 && P.zeta != 1) { msg = "zeta must be 0 (WW) or 1 (GA)"; return 1; }
  if (P.csl_mode != 0 && P.csl_mode != 1) { msg = "csl_mode must be 0 (paper) or 1 (fork)"; return 2; }
  if (P.cap < 0 || P.cap > 2) { msg = "cap must be 0 (none), 1 (planar) or 2 (smooth)"; return 3; }
  if (!(P.p0 < 0.0)) { msg = "p0 must be negative (compression negative)"; return 4; }
  if (!(P.kappa_hat > 0.0)) { msg = "kappa_hat must be > 0"; return 5; }
  if (!(P.mu0 > 0.0)) { msg = "mu0 must be > 0"; return 6; }
  if (!(P.M > 0.0)) { msg = "M must be > 0"; return 7; }
  if (!(0.0 <= P.N && P.N < 1.0)) { msg = "N must satisfy 0 <= N < 1"; return 8; }
  if (!(0.0 <= P.N_bar && P.N_bar < 1.0)) { msg = "N_bar must satisfy 0 <= N_bar < 1"; return 9; }
  if (P.h < 0.0) { msg = "h must be >= 0"; return 10; }
  // shape-function admissibility, both rho and rho_bar. GA: [7/9, 1]. WW: (1/2, 1] --
  // rho = 1/2 EXACTLY is refused (owner decision 2026-10-01: the compression corner is a vertex).
  const char* names[2] = {"rho", "rho_bar"};
  const double vals[2] = {P.rho, P.rho_bar};
  for (int i = 0; i < 2; ++i) {
    const double r = vals[i];
    if (P.zeta == 1) {
      if (!(7.0 / 9.0 - 1e-15 <= r && r <= 1.0 + 1e-15)) {
        std::snprintf(buf, sizeof(buf), "%s=%.17g outside the convex range [7/9, 1] of zeta='GA'", names[i], r);
        msg = buf; return 11;
      }
    } else {
      if (!(0.5 < r && r <= 1.0 + 1e-15)) {
        std::snprintf(buf, sizeof(buf), "%s=%.17g outside the admissible range (1/2, 1] of zeta='WW' "
                      "(rho = 1/2 exactly is refused: the compression corner is a vertex)", names[i], r);
        msg = buf; return 11;
      }
    }
  }
  // dissipation condition A (sheet S.39), owner-approved hard refusal
  if (P.N_bar > P.N + 1e-15) {
    std::snprintf(buf, sizeof(buf), "dissipation refusal: N_bar=%.17g > N=%.17g (sheet S.39)", P.N_bar, P.N);
    msg = buf; return 12;
  }
  const double bt = detail::beta(P);
  if (P.rho / P.rho_bar < bt - 1e-15) {
    std::snprintf(buf, sizeof(buf), "dissipation refusal: rho/rho_bar=%.6g < (1-N)/(1-N_bar)=%.6g (sheet S.39)",
                  P.rho / P.rho_bar, bt);
    msg = buf; return 13;
  }
  std::string warn;
  if (P.rho > P.rho_bar + 1e-15) {
    warn_rho_gt_rhobar = true;
    std::snprintf(buf, sizeof(buf), "warning: rho=%.17g > rho_bar=%.17g violates AB06's psi_c <= phi_c reading "
                  "(dissipation under reading A is still guaranteed)", P.rho, P.rho_bar);
    warn += buf;
  }
  if (P.chi > 0.0) {
    std::snprintf(buf, sizeof(buf), "%swarning: chi=%.17g > 0: the model expects chi < 0 (D* = chi psi_i)",
                  warn.empty() ? "" : "; ", P.chi);
    warn += buf;
  }
  if (P.csl_mode == 0) {
    if (!(P.lambda_tilde > 0.0)) { msg = "lambda_tilde must be > 0"; return 14; }
  } else {
    if (!(P.lambda_c > 0.0 && P.xi > 0.0 && P.p_a > 0.0)) { msg = "fork CSL needs lambda_c > 0, xi > 0, p_a > 0"; return 15; }
  }
  if (P.cap != 0) {
    if (!(0.0 <= P.c1 && P.c1 <= P.c2 && P.c2 < 1.0)) { msg = "cap needs 0 <= c1 <= c2 < 1"; return 16; }
    if (P.cap == 1 && P.c1 != P.c2) { msg = "planar cap needs c1 == c2 (= chi_cap)"; return 17; }
    if (P.cap == 2 && !(P.c2 > P.c1)) { msg = "smooth cap needs c2 > c1"; return 18; }
  }
  // kernel addition (after O2's checks): every parameter finite
  const double all[20] = {P.p0, P.kappa_hat, P.eps_v0, P.mu0, P.alpha0, P.M, P.N, P.N_bar, P.rho, P.rho_bar,
                          P.chi, P.h, P.lambda_tilde, P.v_c0, P.e0, P.lambda_c, P.xi, P.p_a, P.c1, P.c2};
  if (!detail::all_finite(all, 20)) { msg = "every parameter must be finite"; return 19; }
  msg = warn;
  return 0;
}

// O2 api.initial_state: state at sigma0 (tensor comps) with eps^p = 0, v = v0.
// pi_i0 NaN: pi_i placed so that F(sigma0, pi_i) = 0 (inverse of (S.12); hydrostatic
// sigma0 => the apex through p). Returns 0, or nonzero with msg:
// 1-19 validate refusal codes + 100; 101 p >= 0; 102 elastic inversion failed;
// 103 eta >= M/N (no surface through sigma0); 104 pi_i0 >= 0; 105 v0 not > 0;
// 106 eigen-decomposition failed / non-finite sigma0. On success msg holds the
// validate() warnings, if any.
inline int initialState(const Params& P, const double sigma0[6], double v0, double pi_i0,
                        State& out, std::string& msg)
{
  bool warn = false;
  const int vr = validate(P, msg, warn);
  if (vr) return 100 + vr;
  const std::string warnings = msg;   // kept in msg on success
  msg.clear();
  if (!detail::all_finite(sigma0, 6)) { msg = "sigma0 must be finite"; return 106; }
  if (!(v0 > 0.0) || !std::isfinite(v0)) { msg = "v0 must be finite and > 0"; return 105; }
  double S[3][3], w[3], V[3][3];
  detail::t6_to_m(sigma0, S);
  if (!detail::eig_sym3(S, w, V)) { msg = "eigen-decomposition of sigma0 failed"; return 106; }
  double eps_p[3];
  const int ie = detail::invert_elastic(P, w, eps_p);
  if (ie == 1) { msg = "initial stress must have p < 0"; return 101; }
  if (ie) { msg = "elastic inversion of sigma0 failed"; return 102; }
  detail::Invariants inv;
  detail::invariants(w, inv);
  if (std::isnan(pi_i0)) {
    double eta = 0.0;
    if (!inv.vertex) {
      double z, zy, zyy;
      detail::zeta_y(inv.theta, P.rho, P.zeta, z, zy, zyy);
      eta = -z * inv.q / inv.p;
    }
    if (!detail::pi_of_eta(P, inv.p, eta, pi_i0)) {
      msg = "eta >= M/N: no yield surface through this stress"; return 103;
    }
  }
  if (!(pi_i0 < 0.0) || !std::isfinite(pi_i0)) { msg = "pi_i0 must be finite and < 0"; return 104; }
  detail::from_principal(eps_p, V, out.eps_e);
  out.pi_i = pi_i0;
  out.v = v0;
  out.v0 = v0;
  out.eps_p_v = 0.0;
  out.eps_p_s = 0.0;
  out.D_last = 0.0;
  msg = warnings;
  return 0;
}

inline int step(const Params& P, const State& n, const double deps[6], State& np1,
                double sigma[6], double C[6][6], StepInfo& info)
{
  int finest = 0, finest_sub = 0;
  return detail::step_ex(P, n, deps, np1, sigma, C, info, finest, finest_sub);
}

// sigma(eps_e) by the BA06 energy, in the eigenbasis of eps_e. Non-finite on failure.
inline void stress(const Params& P, const State& s, double sigma[6])
{
  double E[3][3], w[3], V[3][3];
  detail::t6_to_m(s.eps_e, E);
  if (!detail::eig_sym3(E, w, V)) {
    for (int i = 0; i < 6; ++i) sigma[i] = std::nan("");
    return;
  }
  detail::Elastic el;
  detail::elastic(P, w, el);
  detail::from_principal(el.sig, V, sigma);
}

// Elastic tangent a^e at s (S.3) assembled by (S.33) with a~ = a^e (the tangent of an
// elastic step, and of a fresh initial state). Same C convention as step().
inline void elasticTangent(const Params& P, const State& s, double C[6][6])
{
  double E[3][3], w[3], V[3][3];
  detail::t6_to_m(s.eps_e, E);
  if (!detail::eig_sym3(E, w, V)) {
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) C[i][j] = std::nan("");
    return;
  }
  detail::Elastic el;
  detail::elastic(P, w, el);
  detail::tangent_small(el.ae, el.sig, w, V, C);
}

}  // namespace ladruno_norsand

#endif  // LadrunoNorSandKernel_h
