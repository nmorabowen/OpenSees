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
    else {
      opserr << "WARNING LadrunoCohesiveHingeBiaxial tag " << tag
             << ": unknown option '" << flag << "' (ignored)\n";
    }
  }

  return new LadrunoCohesiveHingeBiaxial(tag, mp[0], mp[1], mp[2], mp[3], shape, penaltyRatio);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial(int tag, double mcz, double gfz,
                                                         double mcy, double gfy,
                                                         int shp, double pratio)
  : NDMaterial(tag, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(mcz), Gfz(gfz), Mcy(mcy), Gfy(gfy), shape(shp), penaltyRatio(pratio),
    Kpenz(0.0), Kpeny(0.0), a0z(0.0), a0y(0.0), Kinit(2,2),
    CrMax(0.0), Cstrain(2), Cstress(2), Cwork(0.0),
    Tstrain(2), Tstress(2), Ttangent(2,2), TrMax(0.0), Twork(0.0)
{
  computeDerived();
}

LadrunoCohesiveHingeBiaxial::LadrunoCohesiveHingeBiaxial()
  : NDMaterial(0, ND_TAG_LadrunoCohesiveHingeBiaxial),
    Mcz(0.0), Gfz(0.0), Mcy(0.0), Gfy(0.0), shape(EXP), penaltyRatio(1000.0),
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
  double Gfmix = wz*Gfz + wy*Gfy;
  // Esoft > 0 for EVERY mode mix because computeDerived forces Kpen_i = penaltyRatio*Mc_i^2/(2 Gf_i)
  // with penaltyRatio >= 1000, so S = sum w_i*Mc_i^2/Kpen_i = (2/penaltyRatio)*sum w_i*Gf_i = (2/pr)*Gf_mix
  // and Esoft = Gf_mix*(1 - 1/pr) > 0 strictly. Hence the divisions by S and Esoft below are safe.
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
    // dD/dr at frozen mix; dr/dalpha_j = alpha_j/(a0_j^2 r)
    double dDdr  = (Teff - dTeff*TrMax) / (S*TrMax*TrMax);
    double drdaz = az / (a0z*a0z*TrMax);
    double drday = ay / (a0y*a0y*TrMax);
    // ASYMMETRIC BY DESIGN: row i carries its own Kpen_i, so Ttangent(0,1) != Ttangent(1,0) when
    // Gfz != Gfy. This is the correct Jacobian of a history-dependent (non-potential) damage model,
    // NOT a bug — do not "fix" it into symmetry. The element symmetrizes (Czy = ½(Kc01+Kc10)) before
    // its eigenvalue-floored 2x2 inverse, which is what needs a symmetric matrix.
    Ttangent(0,0) -= Kpenz*az*dDdr*drdaz;
    Ttangent(0,1) -= Kpenz*az*dDdr*drday;
    Ttangent(1,0) -= Kpeny*ay*dDdr*drdaz;
    Ttangent(1,1) -= Kpeny*ay*dDdr*drday;
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
    new LadrunoCohesiveHingeBiaxial(this->getTag(), Mcz, Gfz, Mcy, Gfy, shape, penaltyRatio);
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
  static Vector data(13);   // tag + 6 params + CrMax + Cstrain(2) + Cstress(2) + Cwork
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = Mcz;  data(c++) = Gfz;  data(c++) = Mcy;  data(c++) = Gfy;
  data(c++) = (double)shape;
  data(c++) = penaltyRatio;
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
  static Vector data(13);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHingeBiaxial::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  Mcz = data(c++); Gfz = data(c++); Mcy = data(c++); Gfy = data(c++);
  shape = (int)data(c++);
  penaltyRatio = data(c++);
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
    s << "\"shape\": \"" << (shape == LINEAR ? "linear" : "exp") << "\"}";
    return;
  }
  s << "LadrunoCohesiveHingeBiaxial (coupled Mz-My cohesive interaction hinge)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  Mcz/Gfz: " << Mcz << " / " << Gfz << "   (strong-axis capacity / fracture energy)" << endln;
  s << "  Mcy/Gfy: " << Mcy << " / " << Gfy << "   (weak-axis capacity / fracture energy)" << endln;
  s << "  shape  : " << (shape == LINEAR ? "linear" : "exponential") << endln;
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
    default:
      return NDMaterial::getResponse(responseID, matInfo);
  }
}
