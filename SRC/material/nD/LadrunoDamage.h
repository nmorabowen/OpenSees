/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: shared Lemaitre coupled-ductile-damage law for the J2 family.
//
// The SAME damage kinematics drive BOTH the 3D nDMaterial LadrunoJ2 (via the
// shared return-map kernel LadrunoJ2Kernel.h) and the uniaxial LadrunoUniaxialJ2.
// Keeping it in ONE header is the damage half of the "oracle contract": the
// uniaxial Lemaitre material must reduce to the 3D model under a uniaxial stress
// state to ~1e-12, which only holds if the damage backbone is byte-identical in
// both. (Companion of LadrunoHardening.h, which shares the isotropic law.)
//
// Model (de Souza Neto, Peric & Owen 2008, Ch.12 sec 12.3 / 12.4; Lemaitre 1985):
//   * Strain-equivalence + effective stress  sigma~ = sigma / (1 - D).
//   * Damage energy release rate
//       Y = sigma~_eq^2 / (2E) * R_v ,  R_v = (2/3)(1+nu) + 3(1-2nu)(sigmaH/sigma_eq)^2
//     the triaxiality function R_v carries the stress-triaxiality dependence; it
//     equals 1 at uniaxial tension (sigmaH/sigma_eq = 1/3), so Y reduces to the
//     classical sigma~_eq^2/(2E) there. Using sigmaH^2 = (sigmaH/sigma_eq)^2 q^2
//     gives the clean equivalent form  Y = (a*q^2 + b*sigmaH^2)/(2E)  with
//       a = (2/3)(1+nu),  b = 3(1-2nu),  q = sigma~_eq.
//     Triaxiality ratio is identical for nominal and effective stress (both scale
//     by 1-D), so Y is computed directly from the EFFECTIVE stress.
//   * Damage evolution (backward Euler), accumulated only past a plastic-strain
//     threshold pD:
//       D_{n+1} = D_n + (Y/r)^s * dp_dmg ,  dp_dmg = max(0, p_{n+1} - max(p_n, pD))
//     clamped at the critical value Dc (rupture). r (energy denominator, Lemaitre's
//     "S"), s (exponent), pD (threshold), Dc (critical) are the four constants.
//
// For the Lemaitre strain-equivalence model the effective return map is
// D-INDEPENDENT (sigma~_tr = D^e:(eps - eps^p_n) carries no D), so the plastic
// multiplier is solved on the undamaged effective problem (the existing kernel,
// verbatim) and D is recovered from the converged effective state -- exact, not an
// operator-split approximation. The state update therefore decouples; the
// CONSISTENT TANGENT is where the (Delta-gamma, D) coupling lives, through the
// dD/deps sensitivity assembled by the caller from the partials exposed here.
//
// Header-only, no OpenSees dependency. Written: N. Mora-Bowen (Ladruno), 2026.
// See Ladruno_implementation/15_lemaitre_ductile_damage_adr.md.

#ifndef LadrunoDamage_h
#define LadrunoDamage_h

#include <math.h>

namespace Ladruno {

// ---- damage parameters (POD; copied from the material instance) ------------ //
struct DamageParams {
  bool   on   = false;   // master switch (off => exact undamaged behaviour)
  double r    = 0.0;     // energy denominator S
  double s    = 1.0;     // exponent
  double pD   = 0.0;     // accumulated-plastic-strain damage threshold
  double Dc   = 1.0;     // critical (rupture) damage
  // Crack-band (characteristic-length) regularization factor. Scaling the energy
  // denominator r by (lch_ref/lch)^(1/s) is identically a linear scale lch/lch_ref
  // on the damage increment dD; we carry the latter directly. 1.0 = no
  // regularization (local model). Set per-step by the material from the active
  // element's getCharacteristicLength() when -autoRegularization is given.
  double lchScale = 1.0; // = lch / lch_ref
};

// Recover (E, nu) from the kernel's (K, G) pair so the same damage law serves a
// material parameterised in bulk/shear (3D) or directly in E (uniaxial nu=const).
inline double youngsFromKG(double K, double G) { return 9.0*K*G/(3.0*K + G); }
inline double poissonFromKG(double K, double G) { return (3.0*K - 2.0*G)/(2.0*(3.0*K + G)); }

// Triaxiality function R_v(eta), eta = sigmaH/sigma_eq. R_v(1/3) = 1.
inline double triaxialityRv(double eta, double nu)
{
  const double a = (2.0/3.0)*(1.0 + nu);
  const double b = 3.0*(1.0 - 2.0*nu);
  return a + b*eta*eta;
}

// Damage energy release rate Y = (a q^2 + b sigmaH^2)/(2E), with q = effective von
// Mises stress and sigmaH = effective hydrostatic stress. (Equivalent to
// sigma~_eq^2/(2E)*R_v; this form avoids the eta = sigmaH/q division and its q->0
// singularity.) Returns 0 if E <= 0.
inline double damageReleaseRate(double qEff, double sigHeff, double nu, double E)
{
  if (E <= 0.0) return 0.0;
  const double a = (2.0/3.0)*(1.0 + nu);
  const double b = 3.0*(1.0 - 2.0*nu);
  return (a*qEff*qEff + b*sigHeff*sigHeff)/(2.0*E);
}

// Backward-Euler damage increment dD = (Y/r)^s * dp_dmg over a step, with the
// plastic-strain threshold pD applied to the accumulated plastic strain p.
//   pN, pNp1 : accumulated plastic strain at start / end of step
// Returns the RAW increment (NOT clamped to Dc -- the caller clamps so it can also
// flag rupture). Guards: r > 0, Y >= 0, dp_dmg >= 0.
inline double damageIncrement(double Y, double pN, double pNp1,
                              const DamageParams& dp)
{
  if (!dp.on || dp.r <= 0.0 || Y <= 0.0) return 0.0;
  const double pStart = (pN > dp.pD) ? pN : dp.pD;   // max(pN, pD)
  const double dpDmg  = pNp1 - pStart;
  if (dpDmg <= 0.0) return 0.0;
  return dp.lchScale * pow(Y/dp.r, dp.s) * dpDmg;
}

// Partial derivatives of the (unclamped) damage increment for the consistent
// tangent.  dD = g(Y) * dp_dmg,  g(Y) = (Y/r)^s.
//   ddD_dY   = s (Y/r)^(s-1) (1/r) dp_dmg = s * dD / Y      (Y > 0)
//   ddD_dp   = (Y/r)^s                                       (dp_dmg > 0 region)
// Both zero outside the damaging region (matches damageIncrement returning 0).
inline void damageIncrementPartials(double Y, double pN, double pNp1,
                                    const DamageParams& dp,
                                    double& ddD_dY, double& ddD_dp)
{
  ddD_dY = 0.0;
  ddD_dp = 0.0;
  if (!dp.on || dp.r <= 0.0 || Y <= 0.0) return;
  const double pStart = (pN > dp.pD) ? pN : dp.pD;
  const double dpDmg  = pNp1 - pStart;
  if (dpDmg <= 0.0) return;
  const double g = dp.lchScale * pow(Y/dp.r, dp.s);  // lchScale (Y/r)^s
  ddD_dp = g;
  ddD_dY = dp.s * g / Y * dpDmg;                     // lchScale s (Y/r)^(s-1)(1/r) dp_dmg
}

} // namespace Ladruno

#endif
