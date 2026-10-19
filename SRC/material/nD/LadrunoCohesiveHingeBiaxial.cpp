/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Ladruno: coupled biaxial cohesive moment-rotation law (Mz-My interaction
// surface) for the LadrunoDispBeamColumn3d Tier-2 embedded hinge.
// See LadrunoCohesiveHingeBiaxial.h and Ladruno_implementation/
// 34_ladruno_cohesive_hinge_biaxial_adr.md for the derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoCohesiveHingeBiaxial.h>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

#include <stdlib.h>

static bool lchbTryOptionalDouble(double& out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  char tok[64];
  const char* s = OPS_GetStringFromAll(tok, sizeof(tok));
  OPS_ResetCurrentInputArg(-1);
  if (s == 0 || s[0] == '\0') return false;
  char* endp = 0;
  (void)strtod(s, &endp);
  if (endp == s || *endp != '\0') return false;
  int one = 1;
  return (OPS_GetDoubleInput(&one, &out) == 0);
}

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoCohesiveHingeBiaxial tag Mcz Gfz Mcy Gfy
//        <-exp | -linear> <-penaltyRatio r> <-bk eta> <-pm P0 <alpha>> <-implex>
// ===========================================================================
void* OPS_LadrunoCohesiveHingeBiaxial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoCohesiveHingeBiaxial tag? Mcz? Gfz? Mcy? Gfy? "
           << "<-exp | -linear> <-penaltyRatio r?> <-bk eta?> <-pm P0? <alpha?>> <-implex>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoCohesiveHingeBiaxial tag\n";
    return 0;
  }

  double mp[4];   // Mcz, Gfz, Mcy, Gfy
  numData = 4;
  if (OPS_GetDoubleInput(&numData, mp) < 0) {
    opserr << "WARNING invalid LadrunoCohesiveHingeBiaxial Mcz/Gfz/Mcy/Gfy\n";
    return 0;
  }
  if (mp[0] <= 0.0 || mp[1] <= 0.0 || mp[2] <= 0.0 || mp[3] <= 0.0) {
    opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag
           << ": Mc and Gf (both axes) must be > 0 (got Mcz=" << mp[0] << " Gfz=" << mp[1]
           << " Mcy=" << mp[2] << " Gfy=" << mp[3] << ")\n";
    return 0;
  }

  double penaltyRatio = 1000.0;
  double bkEta        = 1.0;
  int    shape = LadrunoCohesiveHingeBiaxial::EXP;
  double P0 = 0.0;
  double alphaPM = 0.5;
  bool   implex = false;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-exp") == 0 || strcmp(flag, "-exponential") == 0) {
      shape = LadrunoCohesiveHingeBiaxial::EXP;
    }
    else if (strcmp(flag, "-linear") == 0 || strcmp(flag, "-lin") == 0) {
      shape = LadrunoCohesiveHingeBiaxial::LINEAR;
    }
    else if (strcmp(flag, "-penaltyRatio") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &penaltyRatio) < 0) {
        opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag << ": -penaltyRatio needs a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-bk") == 0 || strcmp(flag, "-eta") == 0 || strcmp(flag, "-modemix") == 0) {
      // Benzeggagh-Kenane mode-mix exponent: Gf_mix = Gf_z + (Gf_y - Gf_z) w_y^eta (eta=1 -> linear)
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &bkEta) < 0) {
        opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag << ": -bk needs an exponent value\n";
        return 0;
      }
      if (bkEta <= 0.0) {
        opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag
               << ": -bk eta must be > 0 (got " << bkEta << "); reset to 1\n";
        bkEta = 1.0;
      }
    }
    else if (strcmp(flag, "-pm") == 0 || strcmp(flag, "-PM") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &P0) < 0 || P0 <= 0.0) {
        opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag << ": -pm needs P0 > 0\n";
        return 0;
      }
      double optAlpha = alphaPM;
      if (lchbTryOptionalDouble(optAlpha)) {
        alphaPM = optAlpha;
      }
    }
    else if (strcmp(flag, "-implex") == 0) {
      implex = true;
    }
    else {
      opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag
             << ": unknown option '" << flag << "'\n";
      return 0;
    }
  }

  return new LadrunoCohesiveHingeBiaxial(tag, mp[0], mp[1], mp[2], mp[3],
                                         shape, penaltyRatio, bkEta, P0, alphaPM, implex);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial(int tag, double mcz, double gfz,
                                                         double mcy, double gfy,
                                                         int shp, double pratio, double eta,
                                                         double p0, double apm, bool impl)
  : NDMaterial(tag, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(mcz), Gfz(gfz), Mcy(mcy), Gfy(gfy), shape(shp), penaltyRatio(pratio), bkEta(eta),
    P0(p0), alphaPM(apm), implex(impl), P_axial(0.0),
    Kpenz(0.0), Kpeny(0.0), a0z(0.0), a0y(0.0), Kinit(2,2),
    CrMax(0.0), CrMax_old(0.0), Cstrain(2), Cstress(2), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(2), Tstress(2), Ttangent(2,2), TrMax(0.0), Twork(0.0)
{
  computeDerived();
}

LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial()
  : NDMaterial(0, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(0.0), Gfz(0.0), Mcy(0.0), Gfy(0.0), shape(EXP), penaltyRatio(1000.0), bkEta(1.0),
    P0(0.0), alphaPM(0.5), implex(false), P_axial(0.0),
    Kpenz(0.0), Kpeny(0.0), a0z(0.0), a0y(0.0), Kinit(2,2),
    CrMax(0.0), CrMax_old(0.0), Cstrain(2), Cstress(2), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(2), Tstress(2), Ttangent(2,2), TrMax(0.0), Twork(0.0)
{
  // parameters arrive via recvSelf -> computeDerived
}

LadrunoCohesiveHingeBiaxial::~LadrunoCohesiveHingeBiaxial()
{
}

double LadrunoCohesiveHingeBiaxial::getPMMultiplier(void) const
{
  if (P0 <= 0.0) return 1.0;
  // OpenSees beam basic axial force q(0): positive in tension, negative in compression.
  // Parabolic P-M scaling: f(n) = clamp(1 - alphaPM * n - (1 - alphaPM) * n^2, 0.05, 2.5)
  double n = P_axial / P0;
  double f = 1.0 - alphaPM * n - (1.0 - alphaPM) * n * n;
  if (f < 0.05) f = 0.05;
  else if (f > 2.5) f = 2.5;
  return f;
}

const Matrix& LadrunoCohesiveHingeBiaxial::getInitialTangent(void)
{
  if (P0 > 0.0) {
    double f = getPMMultiplier();
    Kinit(0,0) = Kpenz * f;
    Kinit(1,1) = Kpeny * f;
  } else {
    Kinit(0,0) = Kpenz;
    Kinit(1,1) = Kpeny;
  }
  return Kinit;
}

// ===========================================================================
//  calibration: per-axis guarded penalty floor + activation jumps
// ===========================================================================
void LadrunoCohesiveHingeBiaxial::computeDerived(void)
{
  Kinit.Zero();
  if (Mcz <= 0.0 || Gfz <= 0.0 || Mcy <= 0.0 || Gfy <= 0.0) {
    Kpenz = Kpeny = a0z = a0y = 0.0;
    return;
  }
  double r = (penaltyRatio > 1.0) ? penaltyRatio : 1000.0;
  Kpenz = r * Mcz*Mcz / (2.0*Gfz);   // = penaltyRatio * floor; always above the floor
  Kpeny = r * Mcy*Mcy / (2.0*Gfy);
  a0z   = Mcz / Kpenz;
  a0y   = Mcy / Kpeny;
  Kinit(0,0) = Kpenz;
  Kinit(1,1) = Kpeny;
}

// ===========================================================================
//  effective 1D law in the elliptical norm r at the current mode mix (wz,wy):
//  returns S (peak-scale slope), T_eff(rMax), dT_eff/dr(rMax). int T_eff dr == Gf_mix.
// ===========================================================================
void LadrunoCohesiveHingeBiaxial::effectiveLaw(double wz, double wy, double rMax,
                                               double& S, double& Teff, double& dTeff) const
{
  S = wz*Mcz*Mcz/Kpenz + wy*Mcy*Mcy/Kpeny;
  // Benzeggagh-Kenane fracture-energy mode mix (eta=1 -> the original linear interpolation,
  // computed via the identical expression so the default is bit-for-bit unchanged):
  double Gfmix = (bkEta == 1.0) ? (wz*Gfz + wy*Gfy)
                                : (Gfz + (Gfy - Gfz)*pow(wy, bkEta));
  // Esoft > 0 for EVERY mode mix: S = (2/penaltyRatio)*(wz Gf_z + wy Gf_y) <= (2/pr)*max(Gf), tiny vs
  // Gfmix (which lies in [min Gf, max Gf] for any eta>0), so Esoft = Gfmix - S/2 > 0 strictly and the
  // divisions by S and Esoft below are safe. (For eta=1 this is exactly the old Esoft = Gf_mix(1-1/pr).)
  double Esoft = Gfmix - 0.5*S;

  if (rMax <= 1.0) {                       // pre-peak (closed): T_eff = S r, D = 0
    Teff  = S*rMax;
    dTeff = S;
    return;
  }
  double s = rMax - 1.0;
  if (shape == LINEAR) {
    double rf = 2.0*Esoft / S;
    if (s >= rf) { Teff = 0.0; dTeff = 0.0; }       // fully broken
    else         { Teff = S*(1.0 - s/rf); dTeff = -S/rf; }
  } else {                                 // EXP
    double A = S / Esoft;
    double e = exp(-A*s);
    Teff  =  S*e;
    dTeff = -A*S*e;
  }
}

// ===========================================================================
//  trial state
// ===========================================================================
int LadrunoCohesiveHingeBiaxial::setTrialStrain(const Vector& v)
{
  double az = v(0), ay = v(1);
  Tstrain(0) = az;
  Tstrain(1) = ay;

  double ez = (a0z > 0.0) ? az/a0z : 0.0;
  double ey = (a0y > 0.0) ? ay/a0y : 0.0;
  double r  = sqrt(ez*ez + ey*ey);

  bool loading = (r >= CrMax);
  TrMax = loading ? r : CrMax;

  // mode mix from the current jump direction (instantaneous; radial-path exact)
  double wz, wy;
  if (r > 1.0e-300) { wz = ez*ez/(r*r); wy = ey*ey/(r*r); }
  else              { wz = 1.0;         wy = 0.0; }

  if (implex) {
    double time_factor = 1.0;
    if (ops_Dt > 0.0 && dtime_n_commit > 0.0) {
      time_factor = ops_Dt / dtime_n_commit;
      if (time_factor > 2.0) time_factor = 2.0;
    } else if (implexCommittedSteps == 0) {
      time_factor = 0.0;
    }
    double rMax_tilde = CrMax + time_factor * (CrMax - CrMax_old);
    if (rMax_tilde < CrMax) rMax_tilde = CrMax;

    double S_t, Teff_t, dTeff_t;
    effectiveLaw(wz, wy, rMax_tilde, S_t, Teff_t, dTeff_t);
    double D_t = (rMax_tilde > 1.0e-300 && S_t > 0.0) ? (1.0 - Teff_t/(S_t*rMax_tilde)) : 0.0;
    if (D_t < 0.0) D_t = 0.0;
    if (D_t > 1.0) D_t = 1.0;

    double fpm = (P0 > 0.0) ? getPMMultiplier() : 1.0;
    double secFactor = (1.0 - D_t) * fpm;
    double tanFactor = (1.0 - D_t > 1.0e-6 ? (1.0 - D_t) : 1.0e-6) * fpm;

    Tstress(0) = secFactor * Kpenz * az;
    Tstress(1) = secFactor * Kpeny * ay;
    Ttangent.Zero();
    Ttangent(0,0) = tanFactor * Kpenz;
    Ttangent(1,1) = tanFactor * Kpeny;
    Twork = Cwork + 0.5*((Cstress(0)+Tstress(0))*(Tstrain(0)-Cstrain(0))
                       + (Cstress(1)+Tstress(1))*(Tstrain(1)-Cstrain(1)));
    return 0;
  }

  double S, Teff, dTeff;
  effectiveLaw(wz, wy, TrMax, S, Teff, dTeff);

  // isotropic secant damage D = 1 - T_eff(rMax)/(S rMax)  (0 pre-peak; in [0,1])
  double D = (TrMax > 1.0e-300 && S > 0.0) ? (1.0 - Teff/(S*TrMax)) : 0.0;
  if (D < 0.0) D = 0.0;
  if (D > 1.0) D = 1.0;

  // secant-damaged per-axis penalty
  Tstress(0) = (1.0 - D)*Kpenz*az;
  Tstress(1) = (1.0 - D)*Kpeny*ay;

  // 2x2 tangent: (1-D) diagonal penalty + the cohesive damage-gradient coupling
  Ttangent.Zero();
  Ttangent(0,0) = (1.0 - D)*Kpenz;
  Ttangent(1,1) = (1.0 - D)*Kpeny;
  if (loading && TrMax > 1.0 && S > 0.0) {
    // CONSISTENT off-radial tangent (ADR 34 PR-4b). dD/dalpha_j carries BOTH the
    // elliptical-norm sensitivity  dD/dr * dr/dalpha_j  AND the mode-mix sensitivity
    // dD/dw_z * dw_z/dalpha_j. The mode-mix term VANISHES on the pure axes (e_y or e_z = 0)
    // and along the radial direction (dw_z/dalpha . alpha = 0), so radial/pure-axis paths
    // are bit-unchanged; off-radial it makes Ttangent the TRUE Jacobian of M(alpha). Under the
    // LINEAR mix (bkEta=1) the term is identically zero (D is mix-independent — see below), so it
    // only matters for the Benzeggagh-Kenane nonlinear mix; FD-gated in the material-point battery.
    const double r = TrMax;                                // = current norm while loading
    const double s = r - 1.0;
    // (a) norm sensitivity:  dD/dr (frozen mix) and  dr/dalpha_j = alpha_j/(a0_j^2 r)
    const double dDdr  = (Teff - dTeff*r) / (S*r*r);
    const double drdaz = az / (a0z*a0z*r);
    const double drday = ay / (a0y*a0y*r);
    double dDaz = dDdr*drdaz;
    double dDay = dDdr*drday;
    // (b) mode-mix sensitivity.  For the LINEAR Gf interpolation (eta=1) the BK calibration makes
    //     A = S/Esoft mix-INDEPENDENT, so dD/dw == 0 exactly and the frozen-mix tangent is already
    //     the consistent Jacobian (we add nothing -> the eta=1 default stays bit-identical). For a
    //     NONLINEAR interpolation (eta!=1) Gf_mix is no longer proportional to S, so A and hence D
    //     depend on the mix and the off-radial term is live: the softening ratio rho = T_eff/S
    //     depends on w_z only through A = S/Esoft, so dD/dw_z = -(1/r) d(rho)/dA * dA/dw_z (fixed r).
    if (bkEta != 1.0 && wz > 1.0e-300 && wy > 1.0e-300) {
      const double cz = Mcz*Mcz/Kpenz, cy = Mcy*Mcy/Kpeny; // S = wz*cz + wy*cy (peak-scale, linear)
      const double Sp = cz - cy;                           // dS/dwz
      const double Gfmix = Gfz + (Gfy - Gfz)*pow(wy, bkEta);
      const double Esoft = Gfmix - 0.5*S;
      const double Gp = bkEta*(Gfz - Gfy)*pow(wy, bkEta - 1.0);   // dGf_mix/dwz (chain via wy=1-wz)
      const double Ep = Gp - 0.5*Sp;                       // dEsoft/dwz
      const double Ap = (Sp*Esoft - S*Ep) / (Esoft*Esoft); // d(S/Esoft)/dwz
      double dDdw;
      if (shape == LINEAR) dDdw = (s/(2.0*r)) * Ap;        // rho = 1 - (s/2)*A
      else                 dDdw = (s/r) * (Teff/S) * Ap;   // rho = exp(-A s) = T_eff/S
      // dw_z/dalpha_j from w_z = e_z^2/r^2 (e_i = alpha_i/a0_i)
      const double r4 = r*r*r*r;
      const double dwz_daz =  2.0*ez*ey*ey/(r4*a0z);
      const double dwz_day = -2.0*ez*ez*ey/(r4*a0y);
      dDaz += dDdw*dwz_daz;
      dDay += dDdw*dwz_day;
    }
    // ASYMMETRIC BY DESIGN: row i carries its own Kpen_i, so Ttangent(0,1) != Ttangent(1,0) when
    // Gfz != Gfy. This is the correct Jacobian of a history-dependent (non-potential) damage model,
    // NOT a bug — do not "fix" it into symmetry. The element symmetrizes (Czy = ½(Kc01+Kc10)) before
    // its eigenvalue-floored 2x2 inverse, which is what needs a symmetric matrix.
    Ttangent(0,0) -= Kpenz*az*dDaz;
    Ttangent(0,1) -= Kpenz*az*dDay;
    Ttangent(1,0) -= Kpeny*ay*dDaz;
    Ttangent(1,1) -= Kpeny*ay*dDay;
  }

  if (P0 > 0.0) {
    double f = getPMMultiplier();
    Tstress  *= f;
    Ttangent *= f;
  }

  // path work int M . dalpha (2D trapezoidal from the committed anchor)
  Twork = Cwork + 0.5*((Cstress(0)+Tstress(0))*(Tstrain(0)-Cstrain(0))
                     + (Cstress(1)+Tstress(1))*(Tstrain(1)-Cstrain(1)));
  return 0;
}

int LadrunoCohesiveHingeBiaxial::setTrialStrain(const Vector& v, const Vector&)
{
  return this->setTrialStrain(v);
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoCohesiveHingeBiaxial::commitState(void)
{
  CrMax_old  = CrMax;
  CrMax      = TrMax;
  if (implex) {
    double az = Tstrain(0), ay = Tstrain(1);
    double ez = (a0z > 0.0) ? az/a0z : 0.0;
    double ey = (a0y > 0.0) ? ay/a0y : 0.0;
    double r  = sqrt(ez*ez + ey*ey);
    double wz = (r > 1.0e-300) ? (ez*ez/(r*r)) : 1.0;
    double wy = (r > 1.0e-300) ? (ey*ey/(r*r)) : 0.0;
    double S, Teff, dTeff;
    effectiveLaw(wz, wy, CrMax, S, Teff, dTeff);
    double D = (CrMax > 1.0e-300 && S > 0.0) ? (1.0 - Teff/(S*CrMax)) : 0.0;
    if (D < 0.0) D = 0.0;
    if (D > 1.0) D = 1.0;
    double fpm = (P0 > 0.0) ? getPMMultiplier() : 1.0;
    double Mz_imp = (1.0 - D) * fpm * Kpenz * az;
    double My_imp = (1.0 - D) * fpm * Kpeny * ay;
    Cwork = Cwork + 0.5*((Cstress(0)+Mz_imp)*(az-Cstrain(0))
                       + (Cstress(1)+My_imp)*(ay-Cstrain(1)));
    Cstress(0) = Mz_imp;
    Cstress(1) = My_imp;
    Twork = Cwork;
  } else {
    Cstress(0) = Tstress(0); Cstress(1) = Tstress(1);
    Cwork      = Twork;
  }
  Cstrain(0) = Tstrain(0); Cstrain(1) = Tstrain(1);
  dtime_n_commit = ops_Dt;
  implexCommittedSteps++;
  return 0;
}

int LadrunoCohesiveHingeBiaxial::revertToLastCommit(void)
{
  TrMax = CrMax;
  Tstrain(0) = Cstrain(0); Tstrain(1) = Cstrain(1);
  Tstress(0) = Cstress(0); Tstress(1) = Cstress(1);
  Twork = Cwork;
  // restore a consistent secant tangent at the committed point (diagonal)
  this->setTrialStrain(Cstrain);
  return 0;
}

int LadrunoCohesiveHingeBiaxial::revertToStart(void)
{
  CrMax = 0.0; CrMax_old = 0.0; Cwork = 0.0; Cstrain.Zero(); Cstress.Zero();
  dtime_n_commit = 0.0; implexCommittedSteps = 0;
  TrMax = 0.0; Twork = 0.0; Tstrain.Zero(); Tstress.Zero();
  Ttangent = Kinit;
  if (P0 > 0.0) Ttangent *= getPMMultiplier();
  return 0;
}

// ===========================================================================
//  copy
// ===========================================================================
NDMaterial* LadrunoCohesiveHingeBiaxial::getCopy(void)
{
  LadrunoCohesiveHingeBiaxial* theCopy =
    new LadrunoCohesiveHingeBiaxial(this->getTag(), Mcz, Gfz, Mcy, Gfy,
                                    shape, penaltyRatio, bkEta, P0, alphaPM, implex);
  theCopy->P_axial              = P_axial;
  theCopy->CrMax                = CrMax;
  theCopy->CrMax_old            = CrMax_old;
  theCopy->Cstrain              = Cstrain;
  theCopy->Cstress              = Cstress;
  theCopy->Cwork                = Cwork;
  theCopy->dtime_n_commit       = dtime_n_commit;
  theCopy->implexCommittedSteps = implexCommittedSteps;
  theCopy->TrMax                = TrMax;
  theCopy->Tstrain              = Tstrain;
  theCopy->Tstress              = Tstress;
  theCopy->Ttangent             = Ttangent;
  theCopy->Twork                = Twork;
  return theCopy;
}

NDMaterial* LadrunoCohesiveHingeBiaxial::getCopy(const char*)
{
  // order-2 hinge law: only one representation; the element uses getCopy(void)
  return this->getCopy();
}

// ===========================================================================
//  parallel / database
// ===========================================================================
int LadrunoCohesiveHingeBiaxial::sendSelf(int cTag, Channel& theChannel)
{
  static thread_local Vector data(25);   // tag + 7 params + CrMax + Cstrain(2) + Cstress(2) + Cwork + 11 extended
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = Mcz;  data(c++) = Gfz;  data(c++) = Mcy;  data(c++) = Gfy;
  data(c++) = (double)shape;
  data(c++) = penaltyRatio;
  data(c++) = bkEta;
  data(c++) = CrMax;
  data(c++) = Cstrain(0); data(c++) = Cstrain(1);
  data(c++) = Cstress(0); data(c++) = Cstress(1);
  data(c++) = Cwork;
  data(c++) = P0;
  data(c++) = alphaPM;
  data(c++) = implex ? 1.0 : 0.0;
  data(c++) = P_axial;
  data(c++) = CrMax_old;
  data(c++) = dtime_n_commit;
  data(c++) = (double)implexCommittedSteps;
  data(c++) = Tstress(0); data(c++) = Tstress(1);
  data(c++) = Ttangent(0, 0); data(c++) = Ttangent(1, 1);

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHingeBiaxial::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoCohesiveHingeBiaxial::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static thread_local Vector data(25);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHingeBiaxial::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  Mcz = data(c++); Gfz = data(c++); Mcy = data(c++); Gfy = data(c++);
  shape = (int)data(c++);
  penaltyRatio = data(c++);
  bkEta = data(c++);
  CrMax = data(c++);
  Cstrain(0) = data(c++); Cstrain(1) = data(c++);
  Cstress(0) = data(c++); Cstress(1) = data(c++);
  Cwork = data(c++);
  P0 = data(c++);
  alphaPM = data(c++);
  implex = (data(c++) != 0.0);
  P_axial = data(c++);
  CrMax_old = data(c++);
  dtime_n_commit = data(c++);
  implexCommittedSteps = (int)data(c++);
  double savedTstress0 = data(c++); double savedTstress1 = data(c++);
  double savedTtan00   = data(c++); double savedTtan11   = data(c++);

  computeDerived();
  if (implex) {
    Tstrain = Cstrain;
    TrMax   = CrMax;
    Twork   = Cwork;
    Tstress(0) = savedTstress0;
    Tstress(1) = savedTstress1;
    Ttangent.Zero();
    Ttangent(0, 0) = savedTtan00;
    Ttangent(1, 1) = savedTtan11;
  } else {
    this->revertToLastCommit();   // sync trial to committed
  }
  return 0;
}

// ===========================================================================
//  print
// ===========================================================================
void LadrunoCohesiveHingeBiaxial::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoCohesiveHingeBiaxial\", ";
    s << "\"Mcz\": " << Mcz << ", \"Gfz\": " << Gfz << ", ";
    s << "\"Mcy\": " << Mcy << ", \"Gfy\": " << Gfy << ", ";
    s << "\"shape\": \"" << (shape == LINEAR ? "linear" : "exp") << "\", ";
    s << "\"bkEta\": " << bkEta << "}";
    return;
  }
  s << "LadrunoCohesiveHingeBiaxial (coupled Mz-My cohesive interaction hinge)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  Mcz/Gfz: " << Mcz << " / " << Gfz << "   (strong-axis capacity / fracture energy)" << endln;
  s << "  Mcy/Gfy: " << Mcy << " / " << Gfy << "   (weak-axis capacity / fracture energy)" << endln;
  s << "  shape  : " << (shape == LINEAR ? "linear" : "exponential") << endln;
  s << "  bkEta  : " << bkEta << "   (Benzeggagh-Kenane mode-mix exponent; 1 = linear)" << endln;
  if (P0 > 0.0)
    s << "  P0     : " << P0 << "   alphaPM: " << alphaPM << endln;
  if (implex)
    s << "  implex : on" << endln;
  s << "  Kpenz  : " << Kpenz << "   Kpeny: " << Kpeny << endln;
}

// ===========================================================================
//  recordable responses (scalars — sufficient for the element gates)
// ===========================================================================
Response* LadrunoCohesiveHingeBiaxial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "energy") == 0 || strcmp(a, "work") == 0 || strcmp(a, "dissipatedEnergy") == 0)
    return new MaterialResponse(this, 200, 0.0);
  if (strcmp(a, "momentZ") == 0 || strcmp(a, "Mz") == 0)
    return new MaterialResponse(this, 201, 0.0);
  if (strcmp(a, "momentY") == 0 || strcmp(a, "My") == 0)
    return new MaterialResponse(this, 202, 0.0);
  if (strcmp(a, "jumpZ") == 0 || strcmp(a, "alphaZ") == 0)
    return new MaterialResponse(this, 203, 0.0);
  if (strcmp(a, "jumpY") == 0 || strcmp(a, "alphaY") == 0)
    return new MaterialResponse(this, 204, 0.0);
  if (strcmp(a, "maxNorm") == 0 || strcmp(a, "rMax") == 0)
    return new MaterialResponse(this, 205, 0.0);
  if (strcmp(a, "damage") == 0 || strcmp(a, "D") == 0)
    return new MaterialResponse(this, 206, 0.0);
  if (strcmp(a, "tangent") == 0 || strcmp(a, "stiffness") == 0)   // raw 2x2 [K00,K01,K10,K11]
    return new MaterialResponse(this, 207, Vector(4));
  if (strcmp(a, "pmMultiplier") == 0 || strcmp(a, "fpm") == 0)
    return new MaterialResponse(this, 208, 0.0);
  if (strcmp(a, "implexError") == 0 || strcmp(a, "ImplexError") == 0)
    return new MaterialResponse(this, 209, 0.0);

  return NDMaterial::setResponse(argv, argc, s);
}

double LadrunoCohesiveHingeBiaxial::getImplexError(void) const
{
  if (!implex) return 0.0;
  double az = Tstrain(0), ay = Tstrain(1);
  double ez = (a0z > 0.0) ? az/a0z : 0.0;
  double ey = (a0y > 0.0) ? ay/a0y : 0.0;
  double r  = sqrt(ez*ez + ey*ey);
  double wz = (r > 1.0e-300) ? (ez*ez/(r*r)) : 1.0;
  double wy = (r > 1.0e-300) ? (ey*ey/(r*r)) : 0.0;
  double S, Teff, dTeff;
  effectiveLaw(wz, wy, TrMax, S, Teff, dTeff);
  double D = (TrMax > 1.0e-300 && S > 0.0) ? (1.0 - Teff/(S*TrMax)) : 0.0;
  if (D < 0.0) D = 0.0;
  if (D > 1.0) D = 1.0;
  double fpm = (P0 > 0.0) ? getPMMultiplier() : 1.0;
  double Mz_imp = (1.0 - D) * fpm * Kpenz * az;
  double My_imp = (1.0 - D) * fpm * Kpeny * ay;
  double dMz = Tstress(0) - Mz_imp;
  double dMy = Tstress(1) - My_imp;
  double diffNorm = sqrt(dMz*dMz + dMy*dMy);
  double impNorm  = sqrt(Mz_imp*Mz_imp + My_imp*My_imp);
  return diffNorm / (impNorm + 1.0e-12);
}

int LadrunoCohesiveHingeBiaxial::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 200: matInfo.theDouble = Twork;      return 0;
    case 201: matInfo.theDouble = Tstress(0); return 0;
    case 202: matInfo.theDouble = Tstress(1); return 0;
    case 203: matInfo.theDouble = Tstrain(0); return 0;
    case 204: matInfo.theDouble = Tstrain(1); return 0;
    case 205: matInfo.theDouble = TrMax;      return 0;
    case 206: {                                          // isotropic damage at r_max
      double ez = (a0z > 0.0) ? Tstrain(0)/a0z : 0.0;
      double ey = (a0y > 0.0) ? Tstrain(1)/a0y : 0.0;
      double r  = sqrt(ez*ez + ey*ey);
      double wz, wy;
      if (r > 1.0e-300) { wz = ez*ez/(r*r); wy = ey*ey/(r*r); }
      else              { wz = 1.0;         wy = 0.0; }
      double S, Teff, dTeff;
      effectiveLaw(wz, wy, TrMax, S, Teff, dTeff);
      double D = (TrMax > 1.0e-300 && S > 0.0) ? (1.0 - Teff/(S*TrMax)) : 0.0;
      matInfo.theDouble = (D < 0.0) ? 0.0 : (D > 1.0 ? 1.0 : D);
      return 0;
    }
    case 207: {                                            // raw (un-symmetrized) 2x2 tangent
      static thread_local Vector t(4);
      t(0) = Ttangent(0,0); t(1) = Ttangent(0,1);
      t(2) = Ttangent(1,0); t(3) = Ttangent(1,1);
      return matInfo.setVector(t);
    }
    case 208:
      matInfo.theDouble = getPMMultiplier();
      return 0;
    case 209:
      matInfo.theDouble = getImplexError();
      return 0;
    default:
      return NDMaterial::getResponse(responseID, matInfo);
  }
}

int LadrunoCohesiveHingeBiaxial::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "Mcz") == 0 || strcmp(argv[0], "Mc") == 0) return param.addObject(1, this);
  if (strcmp(argv[0], "Gfz") == 0 || strcmp(argv[0], "Gf") == 0) return param.addObject(2, this);
  if (strcmp(argv[0], "Mcy") == 0)                               return param.addObject(3, this);
  if (strcmp(argv[0], "Gfy") == 0)                               return param.addObject(4, this);
  if (strcmp(argv[0], "P0") == 0)                                return param.addObject(5, this);
  if (strcmp(argv[0], "alphaPM") == 0)                           return param.addObject(6, this);
  if (strcmp(argv[0], "bkEta") == 0 || strcmp(argv[0], "eta") == 0) return param.addObject(7, this);
  if (strcmp(argv[0], "penaltyRatio") == 0)                      return param.addObject(8, this);
  return -1;
}

int LadrunoCohesiveHingeBiaxial::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 1: Mcz = info.theDouble; computeDerived(); return 0;
    case 2: Gfz = info.theDouble; computeDerived(); return 0;
    case 3: Mcy = info.theDouble; computeDerived(); return 0;
    case 4: Gfy = info.theDouble; computeDerived(); return 0;
    case 5: P0 = info.theDouble; return 0;
    case 6: alphaPM = info.theDouble; return 0;
    case 7: bkEta = info.theDouble; return 0;
    case 8: penaltyRatio = info.theDouble; computeDerived(); return 0;
    default: return -1;
  }
}
