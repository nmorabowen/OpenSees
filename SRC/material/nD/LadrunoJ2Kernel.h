/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// LadrunoJ2Kernel — the PURE numerical core of the combined-hardening von Mises
// (J2) return map. Plain doubles + <math.h>, NO OpenSees dependency, so it can be
// unit-tested standalone against a numpy oracle (tests/ladrunoj2_reference.py)
// WITHOUT building OpenSees, and so the SAME verified return map serves BOTH:
//   * the small-strain nDMaterial LadrunoJ2 (which delegates integrate() here), and
//   * the finite-strain path — LadrunoJ2 wrapped by LogStrainNDMaterial (the
//     Hencky log-strain adaptor reuses the inner small-strain map verbatim, Box
//     14.3/14.4 of de Souza Neto, Perić & Owen 2008). This is the kernel-sharing
//     that gives `LadrunoBrick -geom finite` large-strain combined-hardening J2.
//
// Model:
//   * Isotropic: Voce saturation + linear  sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso*p
//   * Kinematic: Chaboche superposed Armstrong-Frederick, alpha = sum_k alpha_k
//                alpha_k,{n+1} = (alpha_k,n + (2/3) C_k dG n) / (1 + sqrt(2/3) g_k dG)
//                (g_k = 0 => linear Prager term; N=1 => single AF)
//   * Integration: implicit backward-Euler reduced to a SINGLE scalar Newton on
//     dGamma (Kobayashi-Ohno 2002): direction n = M(dG)/||M(dG)|| with
//     M(dG) = s_tr - sum_k alpha_k,n/(1 + sqrt(2/3) g_k dG); ||xi|| = ||M|| - theta(dG).
//   * Analytic consistent tangent (dSNPO 2008 eq 7.213 structure + AF dn/deps term).
//
// Tensor convention: symmetric tensors stored as 6 TENSOR components in the
// canonical OpenSees order {00,11,22,01,12,02} (the off-diagonal slots hold TRUE
// tensor components, NOT engineering — the engineering<->tensor conversion is the
// caller's responsibility, exactly as in LadrunoJ2.cpp). The 6x6 tangent Dtan is
// the engineering "J2ThreeDimensional convention" matrix (deviatoric projector
// IIDEV6 has 0.5 on the shear diagonal), so it maps engineering-shear strain
// increments to true-stress increments — the same matrix J2ThreeDimensional and
// LogStrainKernel consume.

#ifndef LadrunoJ2Kernel_h
#define LadrunoJ2Kernel_h

#include <math.h>

namespace ladruno_j2_kernel {

// Max number of Chaboche backstress terms (matches LadrunoJ2::MAXBACK). Sizes
// the local scratch arrays; nBack must satisfy 0 <= nBack <= MAXBACK.
const int MAXBACK = 8;

// Double contraction A:B of two symmetric tensors stored as 6 tensor components,
// with the factor-2 weights on the off-diagonal pairs.
inline double dotT(const double* a, const double* b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
       + 2.0*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]);
}

// 6x6 deviatoric projector (Belytschko IIdev) in the J2ThreeDimensional mapping:
// normal block (2/3,-1/3,...), shear-diagonal 0.5.
const double IIDEV6[6][6] = {
  { 2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0 },
  {-1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0 },
  {-1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0 },
  { 0.0, 0.0, 0.0, 0.5, 0.0, 0.0 },
  { 0.0, 0.0, 0.0, 0.0, 0.5, 0.0 },
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5 }
};

// ---- material parameters (POD; copied from the LadrunoJ2 instance) ---------- //
struct Params {
  double K;       // bulk modulus
  double G;       // shear modulus
  double sig0;    // initial yield stress
  double Qinf;    // Voce saturation increment
  double bIso;    // Voce saturation rate
  double Hiso;    // linear isotropic modulus
  int    nBack;   // number of backstress terms
  double C[MAXBACK];    // AF kinematic moduli C_k
  double gam[MAXBACK];  // AF recall constants gamma_k
};

// ---- hardening law ---------------------------------------------------------- //
inline double yieldStress(const Params& p, double pbar)
{
  // Voce saturation + linear; reduces to perfect (Qinf=Hiso=0) or linear (Qinf=0)
  return p.sig0 + p.Qinf*(1.0 - exp(-p.bIso*pbar)) + p.Hiso*pbar;
}

inline double yieldSlope(const Params& p, double pbar)
{
  return p.Qinf*p.bIso*exp(-p.bIso*pbar) + p.Hiso;
}

// ---- elastic tangent (engineering 6x6, J2-3D convention) -------------------- //
inline void elasticTangent(const Params& p, double Kt[6][6])
{
  const double G = p.G, K = p.K;
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) {
      double oneone = (I < 3 && J < 3) ? 1.0 : 0.0;
      Kt[I][J] = K*oneone + 2.0*G*IIDEV6[I][J];
    }
}

// ===========================================================================
//  The return map (scalar Newton, Kobayashi-Ohno) + consistent tangent.
//
//  Inputs:
//    p          material parameters
//    strain6    trial TOTAL strain (tensor components {00,11,22,01,12,02})
//    epsP_n     committed plastic strain (deviatoric, tensor comps)
//    ebarP_n    committed accumulated equivalent plastic strain
//    alpha_n    committed backstress terms [nBack][6] (tensor comps)
//
//  Outputs (all caller-allocated):
//    stress6    Cauchy stress (tensor comps; shear = true stress)
//    Dtan       algorithmic tangent (engineering 6x6, J2-3D convention)
//    epsP       updated plastic strain
//    ebarP      updated accumulated equivalent plastic strain
//    alpha      updated backstress terms [nBack][6]
//    dGamma     plastic multiplier increment (IMPL-EX hook / element record)
//
//  Returns a status code (the caller may turn it into an opserr warning):
//    STATUS_OK (0)           converged / elastic
//    STATUS_SINGULAR (1)     singular local Jacobian (dR == 0) — loop aborted
//    STATUS_NO_CONVERGE (2)  scalar Newton hit maxIter without |R| <= tol
// ===========================================================================
enum { STATUS_OK = 0, STATUS_SINGULAR = 1, STATUS_NO_CONVERGE = 2 };

// outResidual (optional): on the plastic path, receives the final |R| of the
// scalar Newton — lets the caller reprint the original "|R|=..." diagnostic.
inline int returnMap(const Params& p,
                     const double strain6[6],
                     const double epsP_n[6], double ebarP_n,
                     const double alpha_n[][6],
                     double stress6[6], double Dtan[6][6],
                     double epsP[6], double& ebarP,
                     double alpha[][6], double& dGamma,
                     double* outResidual = 0)
{
  const double G = p.G, K = p.K;
  const int nBack = p.nBack;
  const double root23 = sqrt(2.0/3.0);

  // copy committed history into the trial outputs (elastic-step default)
  for (int i = 0; i < 6; i++) epsP[i] = epsP_n[i];
  ebarP = ebarP_n;
  dGamma = 0.0;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha[k][i] = alpha_n[k][i];

  // deviatoric strain (tensor comps)
  double tr = strain6[0] + strain6[1] + strain6[2];
  double edev[6];
  for (int i = 0; i < 6; i++) edev[i] = strain6[i];
  for (int i = 0; i < 3; i++) edev[i] -= tr/3.0;

  // trial deviatoric stress: s_tr = 2G (edev - epsP_n)
  double s_tr[6];
  for (int i = 0; i < 6; i++) s_tr[i] = 2.0*G*(edev[i] - epsP_n[i]);

  // trial relative stress M_tr = s_tr - sum_k alpha_k,n  (dG = 0 => denominators 1)
  double M[6];
  for (int i = 0; i < 6; i++) {
    M[i] = s_tr[i];
    for (int k = 0; k < nBack; k++) M[i] -= alpha_n[k][i];
  }
  double normM = sqrt(dotT(M, M));
  double sy0   = yieldStress(p, ebarP_n);
  double f_tr  = normM - root23*sy0;

  // Tolerances scaled to the operative stress magnitude (unit-aware; never
  // collapses to a fixed absolute when sig0 == 0, e.g. kinematic-only). The
  // normFloor guards the 1/||M|| divisions below — the relative-stress norm
  // ||M(dG)|| can collapse toward zero at a sharp reversal with backstress
  // (J2Plasticity guards the analogous norm_tau division; we restore that).
  const double stressScale = fmax(normM, root23*sy0);
  const double tolY      = 1.0e-12 * fmax(stressScale, 1.0e-300);
  const double normFloor = 1.0e-10 * fmax(stressScale, 1.0e-300);

  if (f_tr <= tolY) {
    // ---- elastic step: state frozen (already copied above) ----
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    elasticTangent(p, Dtan);
    return STATUS_OK;
  }

  // ---- plastic corrector: scalar Newton on dG ----
  int status = STATUS_NO_CONVERGE;
  double dG = 0.0;
  double Dk[MAXBACK];
  double Mp[6];            // d M / d dG
  double theta = 0.0, dtheta = 0.0, R = 1.0;
  int iter = 0;
  const int maxIter = 50;

  while (iter < maxIter) {
    // denominators
    for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*p.gam[k]*dG;

    // M(dG) and its derivative wrt dG
    for (int i = 0; i < 6; i++) {
      M[i]  = s_tr[i];
      Mp[i] = 0.0;
      for (int k = 0; k < nBack; k++) {
        double inv = 1.0/Dk[k];
        M[i]  -= alpha_n[k][i]*inv;
        Mp[i] += alpha_n[k][i]*root23*p.gam[k]*inv*inv;
      }
    }
    normM = sqrt(dotT(M, M));

    // theta(dG) = 2G dG + sum_k (2/3) C_k dG / Dk ; dtheta/ddG
    theta = 2.0*G*dG;
    dtheta = 2.0*G;
    for (int k = 0; k < nBack; k++) {
      double inv = 1.0/Dk[k];
      theta  += (2.0/3.0)*p.C[k]*dG*inv;
      dtheta += (2.0/3.0)*p.C[k]*inv*inv;   // d(dG/Dk)/ddG = 1/Dk^2
    }

    double xiNorm = normM - theta;                      // ||xi_{n+1}||
    double pbar   = ebarP_n + root23*dG;
    R = xiNorm - root23*yieldStress(p, pbar);

    // dR/ddG = (M:Mp)/||M|| - dtheta - (2/3) sig_y'   (guard 1/||M||)
    double dnormM = (normM > normFloor) ? dotT(M, Mp)/normM : 0.0;
    double dR = dnormM - dtheta - (2.0/3.0)*yieldSlope(p, pbar);

    if (fabs(R) <= tolY) { status = STATUS_OK; break; }
    if (dR == 0.0) { status = STATUS_SINGULAR; break; }  // singular local Jacobian
    dG -= R/dR;
    if (dG < 0.0) dG = 0.0;   // keep admissible
    iter++;
  }
  if (outResidual) *outResidual = fabs(R);   // final scalar-Newton residual

  // recompute final quantities at converged dG
  for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*p.gam[k]*dG;
  for (int i = 0; i < 6; i++) {
    M[i]  = s_tr[i];
    Mp[i] = 0.0;
    for (int k = 0; k < nBack; k++) {
      double inv = 1.0/Dk[k];
      M[i]  -= alpha_n[k][i]*inv;
      Mp[i] += alpha_n[k][i]*root23*p.gam[k]*inv*inv;
    }
  }
  normM = sqrt(dotT(M, M));
  if (normM <= normFloor) {
    // Degenerate relative-stress direction (||M|| -> 0): no well-defined plastic
    // flow. Fall back to the elastic predictor at the committed state (history
    // already copied above) rather than dividing by ~0.
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    elasticTangent(p, Dtan);
    return status;
  }
  double n[6];
  for (int i = 0; i < 6; i++) n[i] = M[i]/normM;
  double pbar = ebarP_n + root23*dG;

  // ---- state update ----
  ebarP = pbar;
  dGamma = dG;
  for (int i = 0; i < 6; i++) epsP[i] = epsP_n[i] + dG*n[i];
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++)
      alpha[k][i] = (alpha_n[k][i] + (2.0/3.0)*p.C[k]*dG*n[i]) / Dk[k];

  // ---- stress ----
  for (int i = 0; i < 6; i++) stress6[i] = s_tr[i] - 2.0*G*dG*n[i];
  for (int i = 0; i < 3; i++) stress6[i] += K*tr;

  // ---- consistent tangent (dSNPO 7.213 structure + AF dn term) ----
  dtheta = 2.0*G;
  for (int k = 0; k < nBack; k++) {
    double inv = 1.0/Dk[k];
    dtheta += (2.0/3.0)*p.C[k]*inv*inv;
  }
  double nMp = dotT(n, Mp);
  double h   = dtheta + (2.0/3.0)*yieldSlope(p, pbar) - nMp;   // = -df/ddG > 0

  double Mperp[6];
  for (int i = 0; i < 6; i++) Mperp[i] = Mp[i] - nMp*n[i];

  double G2 = G*G;
  double beta1   = 1.0 - 2.0*G*dG/normM;          // coeff on 2G * IIdev
  double betaNN  = 4.0*G2*dG/normM - 4.0*G2/h;    // coeff on n (x) n
  double betaMpN = -4.0*G2*dG/(h*normM);          // coeff on Mperp (x) n

  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) {
      double oneone = (I < 3 && J < 3) ? 1.0 : 0.0;
      Dtan[I][J] = K*oneone
                 + 2.0*G*beta1*IIDEV6[I][J]
                 + betaNN*n[I]*n[J]
                 + betaMpN*Mperp[I]*n[J];
    }

  return status;
}

} // namespace ladruno_j2_kernel

#endif
