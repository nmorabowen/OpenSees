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

//
// Scalar type (WP-262): the undamaged return map is written ONCE as
// returnMapT<T>. T = double is the production path (returnMap() below is the
// unchanged double-only entry point every caller uses). T = ladruno_dual::Dual<N>
// pushes forward-mode derivatives through the SAME code: seeding the strain
// gives the exact algorithmic tangent, the oracle the analytic Dtan is checked
// against (tests/ladrunoj2_dual_check.cpp). Branches read value(x); the value
// part of a Dual run follows the double run bit for bit.

#include <math.h>
#include <type_traits>
#include <LadrunoDual.h>        // WP-262: dual numbers (forward-mode AD)
#include <LadrunoHardening.h>   // shared Voce+linear sig_y (oracle contract w/ LadrunoUniaxialJ2)
#include <LadrunoDamage.h>      // shared Lemaitre ductile-damage law (oracle contract w/ LadrunoUniaxialJ2)

namespace ladruno_j2_kernel {

using ladruno_dual::value;      // identity on double, the value part of a Dual

// Max number of Chaboche backstress terms (matches LadrunoJ2::MAXBACK). Sizes
// the local scratch arrays; nBack must satisfy 0 <= nBack <= MAXBACK.
const int MAXBACK = 8;

// Double contraction A:B of two symmetric tensors stored as 6 tensor components,
// with the factor-2 weights on the off-diagonal pairs.
template <class T>
inline T dotT(const T* a, const T* b)
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
  // Lemaitre coupled ductile damage (default-constructed OFF; existing aggregate
  // initializers that stop at gam leave this off => returnMap is byte-identical).
  Ladruno::DamageParams dmg;
};

// ---- hardening law (delegated to the SHARED Ladruno::* backbone) ------------ //
// Using LadrunoHardening.h keeps sig_y(pbar) byte-identical to the uniaxial
// LadrunoUniaxialJ2 (the V7 1e-12 oracle contract) AND across the finite-strain
// path that reuses this kernel.
template <class T>
inline T yieldStress(const Params& p, const T& pbar)
{
  return Ladruno::yieldStressVoceLinear(pbar, p.sig0, p.Qinf, p.bIso, p.Hiso);
}

template <class T>
inline T yieldSlope(const Params& p, const T& pbar)
{
  return Ladruno::yieldSlopeVoceLinear(pbar, p.sig0, p.Qinf, p.bIso, p.Hiso);
}

// ---- elastic tangent (engineering 6x6, J2-3D convention) -------------------- //
template <class T>
inline void elasticTangent(const Params& p, T Kt[6][6])
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

// Plastic-corrector residual of the scalar Newton at a given dG (shared by the
// Newton loop and, for a Dual T, the implicit-function-theorem step). Writes the
// denominators Dk, M(dG), dM/ddG and ||M||; returns R(dG) and dR/ddG.
//   R(dG) = ||M(dG)|| - theta(dG) - sqrt(2/3) sig_y(ebarP_n + sqrt(2/3) dG)
template <class T>
inline void plasticResidual(const Params& p, const T s_tr[6], const T alpha_n[][6],
                            const T& ebarP_n, const T& dG, double normFloor,
                            T Dk[], T M[6], T Mp[6], T& normM, T& R, T& dR)
{
  const double G = p.G;
  const int nBack = p.nBack;
  const double root23 = sqrt(2.0/3.0);

  // denominators
  for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*p.gam[k]*dG;

  // M(dG) and its derivative wrt dG
  for (int i = 0; i < 6; i++) {
    M[i]  = s_tr[i];
    Mp[i] = 0.0;
    for (int k = 0; k < nBack; k++) {
      T inv = 1.0/Dk[k];
      M[i]  -= alpha_n[k][i]*inv;
      Mp[i] += alpha_n[k][i]*root23*p.gam[k]*inv*inv;
    }
  }
  normM = sqrt(dotT(M, M));

  // theta(dG) = 2G dG + sum_k (2/3) C_k dG / Dk ; dtheta/ddG
  T theta = 2.0*G*dG;
  T dtheta = 2.0*G;
  for (int k = 0; k < nBack; k++) {
    T inv = 1.0/Dk[k];
    theta  += (2.0/3.0)*p.C[k]*dG*inv;
    dtheta += (2.0/3.0)*p.C[k]*inv*inv;   // d(dG/Dk)/ddG = 1/Dk^2
  }

  T xiNorm = normM - theta;                           // ||xi_{n+1}||
  T pbar   = ebarP_n + root23*dG;
  R = xiNorm - root23*yieldStress(p, pbar);

  // dR/ddG = (M:Mp)/||M|| - dtheta - (2/3) sig_y'   (guard 1/||M||)
  T dnormM = (value(normM) > normFloor) ? T(dotT(M, Mp)/normM) : T(0.0);
  dR = dnormM - dtheta - (2.0/3.0)*yieldSlope(p, pbar);
}

// outResidual (optional): on the plastic path, receives the final |R| of the
// scalar Newton — lets the caller reprint the original "|R|=..." diagnostic.
//
// returnMapT<T>: the return map for any scalar T (double, or a Dual to get
// derivatives). Every input and output is T so a caller can seed any of them;
// the material parameters stay double (a parameter-sensitivity extension would
// template Params as well).
template <class T>
inline int returnMapT(const Params& p,
                      const T strain6[6],
                      const T epsP_n[6], const T& ebarP_n,
                      const T alpha_n[][6],
                      T stress6[6], T Dtan[6][6],
                      T epsP[6], T& ebarP,
                      T alpha[][6], T& dGamma,
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
  T tr = strain6[0] + strain6[1] + strain6[2];
  T edev[6];
  for (int i = 0; i < 6; i++) edev[i] = strain6[i];
  for (int i = 0; i < 3; i++) edev[i] -= tr/3.0;

  // trial deviatoric stress: s_tr = 2G (edev - epsP_n)
  T s_tr[6];
  for (int i = 0; i < 6; i++) s_tr[i] = 2.0*G*(edev[i] - epsP_n[i]);

  // trial relative stress M_tr = s_tr - sum_k alpha_k,n  (dG = 0 => denominators 1)
  T M[6];
  for (int i = 0; i < 6; i++) {
    M[i] = s_tr[i];
    for (int k = 0; k < nBack; k++) M[i] -= alpha_n[k][i];
  }
  T normM = sqrt(dotT(M, M));
  T sy0   = yieldStress(p, ebarP_n);
  T f_tr  = normM - root23*sy0;

  // Tolerances scaled to the operative stress magnitude (unit-aware; never
  // collapses to a fixed absolute when sig0 == 0, e.g. kinematic-only). The
  // normFloor guards the 1/||M|| divisions below — the relative-stress norm
  // ||M(dG)|| can collapse toward zero at a sharp reversal with backstress
  // (J2Plasticity guards the analogous norm_tau division; we restore that).
  // Tolerances are control flow only: they read the VALUE of a Dual.
  const double stressScale = fmax(value(normM), root23*value(sy0));
  const double tolY      = 1.0e-12 * fmax(stressScale, 1.0e-300);
  const double normFloor = 1.0e-10 * fmax(stressScale, 1.0e-300);

  if (value(f_tr) <= tolY) {
    // ---- elastic step: state frozen (already copied above) ----
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    elasticTangent(p, Dtan);
    return STATUS_OK;
  }

  // ---- plastic corrector: scalar Newton on dG ----
  int status = STATUS_NO_CONVERGE;
  T dG = 0.0;
  T Dk[MAXBACK];
  T Mp[6];                 // d M / d dG
  T R = 1.0, dR = 0.0;
  int iter = 0;
  const int maxIter = 50;

  while (iter < maxIter) {
    plasticResidual(p, s_tr, alpha_n, ebarP_n, dG, normFloor, Dk, M, Mp, normM, R, dR);

    if (fabs(value(R)) <= tolY) { status = STATUS_OK; break; }
    if (value(dR) == 0.0) { status = STATUS_SINGULAR; break; }  // singular local Jacobian
    dG -= R/dR;
    if (value(dG) < 0.0) dG = 0.0;   // keep admissible
    iter++;
  }
  if (outResidual) *outResidual = fabs(value(R));   // final scalar-Newton residual

  // Dual T only: the derivative of an iterate lags its value by one Newton step
  // (it is exact for the PREVIOUS iterate). One more Newton step at the converged
  // dG, with its value part zeroed, makes d(dG) the implicit-function-theorem
  // derivative -R_eps/R_x at the converged point while leaving value(dG)
  // bit-identical to the double run. Compiled out for T = double.
  if constexpr (ladruno_dual::is_dual<T>::value) {
    if (status == STATUS_OK && value(dR) != 0.0) {
      T corr = R/dR;
      ladruno_dual::zero_value(corr);
      dG -= corr;
    }
  }

  // recompute final quantities at converged dG
  for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*p.gam[k]*dG;
  for (int i = 0; i < 6; i++) {
    M[i]  = s_tr[i];
    Mp[i] = 0.0;
    for (int k = 0; k < nBack; k++) {
      T inv = 1.0/Dk[k];
      M[i]  -= alpha_n[k][i]*inv;
      Mp[i] += alpha_n[k][i]*root23*p.gam[k]*inv*inv;
    }
  }
  normM = sqrt(dotT(M, M));
  if (value(normM) <= normFloor) {
    // Degenerate relative-stress direction (||M|| -> 0): no well-defined plastic
    // flow. Fall back to the elastic predictor at the committed state (history
    // already copied above) rather than dividing by ~0.
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    elasticTangent(p, Dtan);
    return status;
  }
  T n[6];
  for (int i = 0; i < 6; i++) n[i] = M[i]/normM;
  T pbar = ebarP_n + root23*dG;

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
  T dtheta = 2.0*G;
  for (int k = 0; k < nBack; k++) {
    T inv = 1.0/Dk[k];
    dtheta += (2.0/3.0)*p.C[k]*inv*inv;
  }
  T nMp = dotT(n, Mp);
  T h   = dtheta + (2.0/3.0)*yieldSlope(p, pbar) - nMp;   // = -df/ddG > 0

  T Mperp[6];
  for (int i = 0; i < 6; i++) Mperp[i] = Mp[i] - nMp*n[i];

  double G2 = G*G;
  T beta1   = 1.0 - 2.0*G*dG/normM;          // coeff on 2G * IIdev
  T betaNN  = 4.0*G2*dG/normM - 4.0*G2/h;    // coeff on n (x) n
  T betaMpN = -4.0*G2*dG/(h*normM);          // coeff on Mperp (x) n

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

// The production entry point (unchanged signature; every caller uses this).
inline int returnMap(const Params& p,
                     const double strain6[6],
                     const double epsP_n[6], double ebarP_n,
                     const double alpha_n[][6],
                     double stress6[6], double Dtan[6][6],
                     double epsP[6], double& ebarP,
                     double alpha[][6], double& dGamma,
                     double* outResidual = 0)
{
  return returnMapT<double>(p, strain6, epsP_n, ebarP_n, alpha_n,
                            stress6, Dtan, epsP, ebarP, alpha, dGamma, outResidual);
}

// ===========================================================================
//  Lemaitre coupled ductile damage layer over the (undamaged) effective return
//  map above.  dSNPO 2008 Ch.12 sec 12.3/12.4.
//
//  KEY STRUCTURE: with strain-equivalence + effective-stress von Mises, the
//  effective trial stress s~_tr = 2G(e - e^p_n) carries NO damage, so the plastic
//  multiplier is solved on the UNDAMAGED effective problem (returnMap, verbatim)
//  and D is recovered from the converged effective state -- EXACT, not an
//  operator-split approximation. The state update decouples; the (Delta-gamma, D)
//  coupling lives entirely in the CONSISTENT TANGENT:
//
//      sigma = (1 - D) sigma~                                       (nominal stress)
//      dsigma/deps = (1 - D) dsigma~/deps  -  sigma~ (x) dD/deps    (consistent)
//
//  with dD/deps = ddD_dp * d(pbar)/deps + ddD_dY * dY/deps, assembled here from the
//  flow direction n, hardening modulus h (recomputed at the converged dGamma), and
//  the effective algorithmic tangent dsigma~/deps returned by returnMap.
//
//  Extra args:  D_n   committed damage (in)
//               D     updated damage (out, clamped to [0, Dc]) -- ALWAYS the
//                     IMPLICIT value, so the caller commits a true backward-Euler D.
//  When p.dmg.on == false this is a pure pass-through of returnMap (D = 0), so a
//  damage-OFF material is byte-identical to the undamaged LadrunoJ2 (regression V0).
//
//  IMPL-EX (Oliver-Huespe-Cante 2008): pass dScaleOverride >= 0 to degrade the
//  stress/tangent by a FROZEN extrapolated damage D~ instead of the implicit D.
//  The plastic return stays implicit (it is hardening => SPD in effective space;
//  ONLY damage softens), so freezing D~ gives a CONSTANT SPD tangent
//  (1 - D~) Dtan_eff -- no -sigma~ (x) dD/deps softening term, hence no Newton
//  stall on the indefinite implicit tangent. D (out) is still the implicit value,
//  so the caller updates the extrapolation from a true backward-Euler history.
//  dScaleOverride < 0 (default) = fully implicit (the consistent tangent above).
// ===========================================================================
template <class T>
inline int returnMapDamagedT(const Params& p,
                             const T strain6[6],
                             const T epsP_n[6], const T& ebarP_n,
                             const T alpha_n[][6], const T& D_n,
                             T stress6[6], T Dtan[6][6],
                             T epsP[6], T& ebarP,
                             T alpha[][6], T& dGamma, T& D,
                             double* outResidual = 0, double dScaleOverride = -1.0)
{
  // 1) effective (undamaged) return map -- the proven kernel, untouched
  int status = returnMapT<T>(p, strain6, epsP_n, ebarP_n, alpha_n,
                             stress6, Dtan, epsP, ebarP, alpha, dGamma, outResidual);

  if (!p.dmg.on) { D = 0.0; return status; }   // pure pass-through (bit-identical)

  const double G = p.G;
  const int nBack = p.nBack;
  const double root23 = sqrt(2.0/3.0);

  // 2) effective deviatoric stress s~ and hydrostatic stress from stress6 (effective)
  T sigEff[6];
  for (int i = 0; i < 6; i++) sigEff[i] = stress6[i];
  const T sigH = (sigEff[0] + sigEff[1] + sigEff[2]) / 3.0;   // effective pressure
  T sdev[6];
  for (int i = 0; i < 6; i++) sdev[i] = sigEff[i];
  for (int i = 0; i < 3; i++) sdev[i] -= sigH;
  const T qEff = sqrt(1.5 * dotT(sdev, sdev));                // effective von Mises

  // 3) damage energy release rate + backward-Euler damage update
  const double E  = Ladruno::youngsFromKG(p.K, p.G);
  const double nu = Ladruno::poissonFromKG(p.K, p.G);
  const T Y  = Ladruno::damageReleaseRate(qEff, sigH, nu, E);
  const T dDraw = Ladruno::damageIncrement(Y, ebarP_n, ebarP, p.dmg);

  T Dnew = D_n + dDraw;
  bool capped = false;
  if (value(Dnew) >= p.dmg.Dc) { Dnew = p.dmg.Dc; capped = true; }
  if (value(Dnew) < 0.0) Dnew = 0.0;
  D = Dnew;                          // ALWAYS the implicit damage (for the caller to commit)

  // IMPL-EX: degrade by the frozen extrapolated D~ (SPD tangent) but still return
  // the implicit D above. Fully-implicit (default) degrades by the implicit D.
  const bool implex = (dScaleOverride >= 0.0);
  T dScale = implex ? T(dScaleOverride) : Dnew;
  if (value(dScale) >= p.dmg.Dc) dScale = p.dmg.Dc;
  if (value(dScale) < 0.0) dScale = 0.0;
  const T omd = 1.0 - dScale;

  // 4) dD/deps for the consistent tangent (only where damage is actively growing).
  //    Skipped entirely in IMPL-EX: D~ is frozen => the rank-one softening term
  //    vanishes => tangent (1 - D~) Dtan_eff is SPD.
  T dDdeps[6] = {0,0,0,0,0,0};
  if (!implex && !capped && value(dDraw) > 0.0 && value(dGamma) > 0.0) {
    // recompute the converged flow direction n and hardening modulus h (the
    // quantities returnMap consumes internally but does not return)
    T tr = strain6[0] + strain6[1] + strain6[2];
    T edev[6];
    for (int i = 0; i < 6; i++) edev[i] = strain6[i];
    for (int i = 0; i < 3; i++) edev[i] -= tr/3.0;
    T s_tr[6];
    for (int i = 0; i < 6; i++) s_tr[i] = 2.0*G*(edev[i] - epsP_n[i]);

    T Dk[MAXBACK];
    for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*p.gam[k]*dGamma;
    T M[6], Mp[6];
    for (int i = 0; i < 6; i++) {
      M[i] = s_tr[i]; Mp[i] = 0.0;
      for (int k = 0; k < nBack; k++) {
        T inv = 1.0/Dk[k];
        M[i]  -= alpha_n[k][i]*inv;
        Mp[i] += alpha_n[k][i]*root23*p.gam[k]*inv*inv;
      }
    }
    T normM = sqrt(dotT(M, M));
    const double normFloor = 1.0e-10 * fmax(fmax(value(normM), root23*value(yieldStress(p, ebarP))), 1.0e-300);
    if (value(normM) > normFloor) {
      T n[6];
      for (int i = 0; i < 6; i++) n[i] = M[i]/normM;
      T dtheta = 2.0*G;
      for (int k = 0; k < nBack; k++) { T inv = 1.0/Dk[k]; dtheta += (2.0/3.0)*p.C[k]*inv*inv; }
      T h = dtheta + (2.0/3.0)*yieldSlope(p, ebarP) - dotT(n, Mp);   // = -df/ddG > 0

      // d(pbar)/deps : pbar = ebarP_n + sqrt(2/3) dGamma ; d(dGamma)/deps_eng = (2G/h) n w
      //   (w = 1 normal, 1/2 shear: tensor->engineering strain on the differentiation var)
      // dY/deps : Y(sigma~) contracted with the effective algorithmic tangent Dtan
      const double a = (2.0/3.0)*(1.0 + nu);
      const double b = 3.0*(1.0 - 2.0*nu);
      T gY[6];                                        // grad_Y in dotT convention
      for (int i = 0; i < 6; i++) {
        T base = 3.0*a*sdev[i];
        if (i < 3) base += (2.0*b/3.0)*sigH;
        gY[i] = base/(2.0*E);
      }
      T ddD_dY, ddD_dp;
      Ladruno::damageIncrementPartials(Y, ebarP_n, ebarP, p.dmg, ddD_dY, ddD_dp);

      const double wt[6] = {1,1,1,2,2,2};             // dotT weights (off-diagonal pairs)
      for (int j = 0; j < 6; j++) {
        // dY/deps_eng,j = sum_I wt_I gY_I Dtan[I][j]   (Dtan already eng-convention)
        T dY_dej = 0.0;
        for (int I = 0; I < 6; I++) dY_dej += wt[I]*gY[I]*Dtan[I][j];
        // d(pbar)/deps_eng,j = sqrt(2/3) (2G/h) n_j.  d||M||/dM carries the dotT weight
        // wt_j (2 on shear); the tensor->engineering strain factor is 1/wt_j, so they
        // cancel and the engineering derivative is just n_j (no shear factor).
        T dpbar_dej = root23*(2.0*G/h)*n[j];
        dDdeps[j] = ddD_dp*dpbar_dej + ddD_dY*dY_dej;
      }
    }
  }

  // 5) degrade stress and assemble the consistent tangent
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] = omd*Dtan[I][J] - sigEff[I]*dDdeps[J];
  for (int i = 0; i < 6; i++) stress6[i] = omd*sigEff[i];

  return status;
}

// The production entry point (unchanged signature; every caller uses this).
inline int returnMapDamaged(const Params& p,
                            const double strain6[6],
                            const double epsP_n[6], double ebarP_n,
                            const double alpha_n[][6], double D_n,
                            double stress6[6], double Dtan[6][6],
                            double epsP[6], double& ebarP,
                            double alpha[][6], double& dGamma, double& D,
                            double* outResidual = 0, double dScaleOverride = -1.0)
{
  return returnMapDamagedT<double>(p, strain6, epsP_n, ebarP_n, alpha_n, D_n,
                                   stress6, Dtan, epsP, ebarP, alpha, dGamma, D,
                                   outResidual, dScaleOverride);
}

} // namespace ladruno_j2_kernel

#endif
