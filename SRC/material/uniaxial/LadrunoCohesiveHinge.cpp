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

// Ladruno: discrete cohesive moment–rotation law M([[theta]]) for the
// LadrunoDispBeamColumn Tier-2 embedded strong-discontinuity hinge.
// See LadrunoCohesiveHinge.h and Ladruno_implementation/
// 32_ladruno_dispbeamcolumn_regularization_adr.md for the derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoCohesiveHinge.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoCohesiveHinge tag Mc Gf
//        <-penalty K | -penaltyRatio r> <-exp | -linear>
// ===========================================================================
void* OPS_LadrunoCohesiveHinge(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoCohesiveHinge tag? Mc? Gf? "
           << "<-penalty K? | -penaltyRatio r?> <-exp | -linear>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoCohesiveHinge tag\n";
    return 0;
  }

  double mp[2];   // Mc, Gf
  numData = 2;
  if (OPS_GetDoubleInput(&numData, mp) < 0) {
    opserr << "WARNING invalid LadrunoCohesiveHinge Mc/Gf\n";
    return 0;
  }
  double Mc = mp[0];
  double Gf = mp[1];
  if (Mc <= 0.0 || Gf <= 0.0) {
    opserr << "WARNING LadrunoCohesiveHinge tag " << tag
           << ": Mc and Gf must be > 0 (got Mc=" << Mc << " Gf=" << Gf << ")\n";
    return 0;
  }

  // defaults: exponential softening, near-rigid default penalty
  double Kpen = -1.0;            // <0 -> constructor builds penaltyRatio*floor
  double penaltyRatio = 1000.0;
  int    shape = LadrunoCohesiveHinge::EXP;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-penalty") == 0 || strcmp(flag, "-Kpen") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &Kpen) < 0) {
        opserr << "WARNING LadrunoCohesiveHinge tag " << tag << ": -penalty needs a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-penaltyRatio") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &penaltyRatio) < 0) {
        opserr << "WARNING LadrunoCohesiveHinge tag " << tag << ": -penaltyRatio needs a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-exp") == 0 || strcmp(flag, "-exponential") == 0) {
      shape = LadrunoCohesiveHinge::EXP;
    }
    else if (strcmp(flag, "-linear") == 0 || strcmp(flag, "-lin") == 0) {
      shape = LadrunoCohesiveHinge::LINEAR;
    }
    else {
      opserr << "WARNING LadrunoCohesiveHinge tag " << tag
             << ": unknown option '" << flag << "' (ignored)\n";
    }
  }

  return new LadrunoCohesiveHinge(tag, Mc, Gf, Kpen, shape, penaltyRatio);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHinge::LadrunoCohesiveHinge(int tag, double mc, double gf,
                                           double kpen, int shp, double pratio)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoCohesiveHinge),
    Mc(mc), Gf(gf), Kpen(kpen), shape(shp), penaltyRatio(pratio),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0)
{
  computeDerived();
  Ttangent = Kpen;
}

LadrunoCohesiveHinge::LadrunoCohesiveHinge()
  : UniaxialMaterial(0, MAT_TAG_LadrunoCohesiveHinge),
    Mc(0.0), Gf(0.0), Kpen(0.0), shape(EXP), penaltyRatio(1000.0),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0)
{
  // parameters arrive via recvSelf -> computeDerived
}

LadrunoCohesiveHinge::~LadrunoCohesiveHinge()
{
}

// ===========================================================================
//  calibration: enforce the guarded penalty floor, build kappa0/aExp/kappaf
//  so the monotonic envelope integrates to EXACTLY Gf.
// ===========================================================================
void LadrunoCohesiveHinge::computeDerived(void)
{
  if (Mc <= 0.0 || Gf <= 0.0) { kappa0 = aExp = kappaf = 0.0; return; }

  const double KpenFloor = Mc*Mc / (2.0*Gf);   // pre-peak energy == Gf at the floor

  // default penalty (Kpen<=0): penaltyRatio multiples above the floor
  if (Kpen <= 0.0) {
    double r = (penaltyRatio > 1.0) ? penaltyRatio : 1000.0;
    Kpen = r * KpenFloor;
  }

  // guarded floor: Kpen must exceed the floor or the softening branch is non-physical
  if (Kpen <= KpenFloor) {
    double rSafe = (penaltyRatio > 1.0) ? penaltyRatio : 1000.0;
    double KpenSafe = rSafe * KpenFloor;
    opserr << "WARNING LadrunoCohesiveHinge tag " << this->getTag()
           << ": penalty Kpen=" << Kpen << " <= floor Mc^2/(2Gf)=" << KpenFloor
           << " (pre-peak energy would exceed Gf); clamping to " << KpenSafe << "\n";
    Kpen = KpenSafe;
  }

  kappa0 = Mc / Kpen;                       // jump at activation
  double Esoft = Gf - Mc*Mc/(2.0*Kpen);     // energy left for the softening branch (>0 by floor)
  aExp   = Mc / Esoft;                       // exp: area Mc/aExp == Esoft
  kappaf = 2.0*Esoft / Mc;                   // lin: area Mc*kappaf/2 == Esoft
}

// ===========================================================================
//  envelope moment Menv(km) and consistent slope dMenv/dkm at km >= 0
// ===========================================================================
void LadrunoCohesiveHinge::envelope(double km, double& Menv, double& dMenv) const
{
  if (km <= kappa0) {                        // pre-peak penalty (hinge closed)
    Menv  = Kpen * km;
    dMenv = Kpen;
    return;
  }
  if (shape == LINEAR) {
    double d = km - kappa0;
    if (d >= kappaf) {                       // fully broken: zero residual moment
      Menv = 0.0; dMenv = 0.0;
    } else {
      Menv  = Mc * (1.0 - d/kappaf);
      dMenv = -Mc/kappaf;
    }
  } else {                                   // EXP
    double e = exp(-aExp * (km - kappa0));
    Menv  =  Mc * e;
    dMenv = -aExp * Mc * e;
  }
}

// ===========================================================================
//  trial state
// ===========================================================================
int LadrunoCohesiveHinge::setTrialStrain(double strain, double)
{
  Tstrain = strain;
  double akappa = fabs(strain);

  // irreversible damage driver: max |jump| ever seen
  bool loading = (akappa >= CkappaMax);
  TkappaMax = loading ? akappa : CkappaMax;

  double Menv, dMenv;
  envelope(TkappaMax, Menv, dMenv);

  // secant-to-origin (isotropic damage); Ksec = Menv(kappaMax)/kappaMax
  double Ksec = (TkappaMax > 0.0) ? (Menv / TkappaMax) : Kpen;

  Tstress = Ksec * strain;                   // sign carried by `strain`

  if (loading && akappa > 0.0)
    Ttangent = dMenv;                        // consistent envelope slope (may be < 0)
  else
    Ttangent = Ksec;                         // unload/reload secant (>= 0)

  // path work int M d[[theta]] (trapezoidal from the committed anchor; exact for the
  // piecewise-linear LINEAR envelope, and the standard accumulation for EXP)
  Twork = Cwork + 0.5*(Cstress + Tstress)*(Tstrain - Cstrain);

  return 0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoCohesiveHinge::commitState(void)
{
  CkappaMax = TkappaMax;
  Cstrain   = Tstrain;
  Cstress   = Tstress;
  Cwork     = Twork;
  return 0;
}

int LadrunoCohesiveHinge::revertToLastCommit(void)
{
  TkappaMax = CkappaMax;
  Tstrain   = Cstrain;
  Tstress   = Cstress;
  Twork     = Cwork;
  // recover a consistent trial tangent at the committed point
  double Menv, dMenv;
  envelope(TkappaMax, Menv, dMenv);
  Ttangent = (TkappaMax > 0.0) ? (Menv / TkappaMax) : Kpen;
  return 0;
}

int LadrunoCohesiveHinge::revertToStart(void)
{
  CkappaMax = 0.0; Cstrain = 0.0; Cstress = 0.0; Cwork = 0.0;
  TkappaMax = 0.0; Tstrain = 0.0; Tstress = 0.0; Twork = 0.0;
  Ttangent = Kpen;
  return 0;
}

// ===========================================================================
//  copy
// ===========================================================================
UniaxialMaterial* LadrunoCohesiveHinge::getCopy(void)
{
  LadrunoCohesiveHinge* theCopy =
    new LadrunoCohesiveHinge(this->getTag(), Mc, Gf, Kpen, shape, penaltyRatio);

  theCopy->CkappaMax = CkappaMax;
  theCopy->Cstrain   = Cstrain;
  theCopy->Cstress   = Cstress;
  theCopy->Cwork     = Cwork;
  theCopy->TkappaMax = TkappaMax;
  theCopy->Tstrain   = Tstrain;
  theCopy->Tstress   = Tstress;
  theCopy->Ttangent  = Ttangent;
  theCopy->Twork     = Twork;
  return theCopy;
}

// ===========================================================================
//  parallel / database
// ===========================================================================
int LadrunoCohesiveHinge::sendSelf(int cTag, Channel& theChannel)
{
  // tag + 4 params (Mc,Gf,Kpen,shape) + penaltyRatio + 4 committed history
  static Vector data(1 + 4 + 1 + 4);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = Mc;
  data(c++) = Gf;
  data(c++) = Kpen;
  data(c++) = (double)shape;
  data(c++) = penaltyRatio;
  data(c++) = CkappaMax;
  data(c++) = Cstrain;
  data(c++) = Cstress;
  data(c++) = Cwork;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHinge::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoCohesiveHinge::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(1 + 4 + 1 + 4);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHinge::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  Mc           = data(c++);
  Gf           = data(c++);
  Kpen         = data(c++);
  shape        = (int)data(c++);
  penaltyRatio = data(c++);
  CkappaMax    = data(c++);
  Cstrain      = data(c++);
  Cstress      = data(c++);
  Cwork        = data(c++);

  computeDerived();           // Kpen already valid (was floor-checked on send side)
  this->revertToLastCommit(); // sync trial to committed
  return 0;
}

// ===========================================================================
//  print
// ===========================================================================
void LadrunoCohesiveHinge::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoCohesiveHinge\", ";
    s << "\"Mc\": " << Mc << ", \"Gf\": " << Gf << ", \"Kpen\": " << Kpen << ", ";
    s << "\"shape\": \"" << (shape == LINEAR ? "linear" : "exp") << "\"}";
    return;
  }
  s << "LadrunoCohesiveHinge (discrete cohesive moment-rotation hinge)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  Mc     : " << Mc << "   (cohesive moment capacity)" << endln;
  s << "  Gf     : " << Gf << "   (fracture energy per hinge)" << endln;
  s << "  Kpen   : " << Kpen << "   (penalty; floor Mc^2/(2Gf)=" << Mc*Mc/(2.0*Gf) << ")" << endln;
  s << "  shape  : " << (shape == LINEAR ? "linear" : "exponential") << endln;
  s << "  kappa0 : " << kappa0 << "   activation jump" << endln;
  if (shape == LINEAR) s << "  kappaf : " << kappaf << "   softening span" << endln;
  else                 s << "  aExp   : " << aExp   << "   softening rate"  << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoCohesiveHinge::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return UniaxialMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "rotationJump") == 0 || strcmp(a, "jump") == 0 || strcmp(a, "kappa") == 0)
    return new MaterialResponse(this, 100, 0.0);
  if (strcmp(a, "maxJump") == 0 || strcmp(a, "kappaMax") == 0)
    return new MaterialResponse(this, 101, 0.0);
  if (strcmp(a, "dissipatedEnergy") == 0 || strcmp(a, "energy") == 0 ||
      strcmp(a, "work") == 0)
    return new MaterialResponse(this, 102, 0.0);
  if (strcmp(a, "damage") == 0 || strcmp(a, "D") == 0)
    return new MaterialResponse(this, 103, 0.0);
  if (strcmp(a, "failed") == 0 || strcmp(a, "open") == 0)
    return new MaterialResponse(this, 104, 0.0);

  return UniaxialMaterial::setResponse(argv, argc, s);
}

int LadrunoCohesiveHinge::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 100: matInfo.theDouble = Tstrain;   return 0;   // rotation jump
    case 101: matInfo.theDouble = TkappaMax; return 0;   // max jump
    case 102: matInfo.theDouble = Twork;     return 0;   // path work int M d[[theta]]
    case 103: {                                           // isotropic damage 1 - Ksec/Kpen
      double Menv, dMenv;
      envelope(TkappaMax, Menv, dMenv);
      double Ksec = (TkappaMax > 0.0) ? (Menv/TkappaMax) : Kpen;
      matInfo.theDouble = (Kpen > 0.0) ? (1.0 - Ksec/Kpen) : 0.0;
      return 0;
    }
    case 104: {                                           // fully-open (broken) flag
      double open = 0.0;
      if (shape == LINEAR && TkappaMax >= kappa0 + kappaf) open = 1.0;
      matInfo.theDouble = open;
      return 0;
    }
    default:
      return UniaxialMaterial::getResponse(responseID, matInfo);
  }
}
