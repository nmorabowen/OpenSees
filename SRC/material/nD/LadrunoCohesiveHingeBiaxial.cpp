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
#include <string.h>
#include <math.h>
#include <elementAPI.h>

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoCohesiveHingeBiaxial tag Mcz Gfz Mcy Gfy
//        <-exp | -linear> <-penaltyRatio r>
// ===========================================================================
void* OPS_LadrunoCohesiveHingeBiaxial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoCohesiveHingeBiaxial tag? Mcz? Gfz? Mcy? Gfy? "
           << "<-exp | -linear> <-penaltyRatio r?>\n";
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
    else {
      opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag
             << ": unknown option '" << flag << "' (ignored)\n";
    }
  }

  return new LadrunoCohesiveHingeBiaxial(tag, mp[0], mp[1], mp[2], mp[3], shape, penaltyRatio, bkEta);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial(int tag, double mcz, double gfz,
                                                         double mcy, double gfy,
                                                         int shp, double pratio, double eta)
  : NDMaterial(tag, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(mcz), Gfz(gfz), Mcy(mcy), Gfy(gfy), shape(shp), penaltyRatio(pratio), bkEta(eta),
    Kpenz(0.0), Kpeny(0.0), a0z(0.0), a0y(0.0), Kinit(2,2),
    CrMax(0.0), Cstrain(2), Cstress(2), Cwork(0.0),
    Tstrain(2), Tstress(2), Ttangent(2,2), TrMax(0.0), Twork(0.0)
{
  computeDerived();
}

LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial()
  : NDMaterial(0, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(0.0), Gfz(0.0), Mcy(0.0), Gfy(0.0), shape(EXP), penaltyRatio(1000.0), bkEta(1.0),
    Kpenz(0.0), Kpeny(0.0), a0z(0.0), a0y(0.0), Kinit(2,2),
    CrMax(0.0), Cstrain(2), Cstress(2), Cwork(0.0),
    Tstrain(2), Tstress(2), Ttangent(2,2), TrMax(0.0), Twork(0.0)
{
  // parameters arrive via recvSelf -> computeDerived
}

LadrunoCohesiveHingeBiaxial::~LadrunoCohesiveHingeBiaxial()
{
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
  CrMax      = TrMax;
  Cstrain(0) = Tstrain(0); Cstrain(1) = Tstrain(1);
  Cstress(0) = Tstress(0); Cstress(1) = Tstress(1);
  Cwork      = Twork;
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
  CrMax = 0.0; Cwork = 0.0; Cstrain.Zero(); Cstress.Zero();
  TrMax = 0.0; Twork = 0.0; Tstrain.Zero(); Tstress.Zero();
  Ttangent = Kinit;
  return 0;
}

// ===========================================================================
//  copy
// ===========================================================================
NDMaterial* LadrunoCohesiveHingeBiaxial::getCopy(void)
{
  LadrunoCohesiveHingeBiaxial* theCopy =
    new LadrunoCohesiveHingeBiaxial(this->getTag(), Mcz, Gfz, Mcy, Gfy, shape, penaltyRatio, bkEta);
  theCopy->CrMax   = CrMax;
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Cwork   = Cwork;
  theCopy->TrMax   = TrMax;
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  theCopy->Twork   = Twork;
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
  static Vector data(14);   // tag + 7 params + CrMax + Cstrain(2) + Cstress(2) + Cwork
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

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHingeBiaxial::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoCohesiveHingeBiaxial::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(14);
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

  computeDerived();
  this->revertToLastCommit();   // sync trial to committed
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

  return NDMaterial::setResponse(argv, argc, s);
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
      static Vector t(4);
      t(0) = Ttangent(0,0); t(1) = Ttangent(0,1);
      t(2) = Ttangent(1,0); t(3) = Ttangent(1,1);
      return matInfo.setVector(t);
    }
    default:
      return NDMaterial::getResponse(responseID, matInfo);
  }
}
