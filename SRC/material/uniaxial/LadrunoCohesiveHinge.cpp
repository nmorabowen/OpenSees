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
//        <-pm P0 <alpha>> <-implex>
// ===========================================================================
static bool lchTryOptionalDouble(double& out)
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

void* OPS_LadrunoCohesiveHinge(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoCohesiveHinge tag? Mc? Gf? "
           << "<-penalty K? | -penaltyRatio r?> <-exp | -linear> <-pm P0? <alpha?>> <-implex>\n";
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
  double P0 = 0.0;
  double alphaPM = 0.5;
  bool   implex = false;

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
    else if (strcmp(flag, "-pm") == 0 || strcmp(flag, "-PM") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &P0) < 0 || P0 <= 0.0) {
        opserr << "WARNING LadrunoCohesiveHinge tag " << tag << ": -pm needs P0 > 0\n";
        return 0;
      }
      double optAlpha = alphaPM;
      if (lchTryOptionalDouble(optAlpha)) {
        alphaPM = optAlpha;
      }
    }
    else if (strcmp(flag, "-implex") == 0) {
      implex = true;
    }
    else {
      opserr << "WARNING LadrunoCohesiveHinge tag " << tag
             << ": unknown option '" << flag << "'\n";
      return 0;
    }
  }

  return new LadrunoCohesiveHinge(tag, Mc, Gf, Kpen, shape, penaltyRatio, P0, alphaPM, implex);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHinge::LadrunoCohesiveHinge(int tag, double mc, double gf,
                                           double kpen, int shp, double pratio,
                                           double p0, double apm, bool impl)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoCohesiveHinge),
    Mc(mc), Gf(gf), Kpen(kpen), userKpen(kpen > 0.0), shape(shp), penaltyRatio(pratio),
    P0(p0), alphaPM(apm), implex(impl), P_axial(0.0),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), CkappaMax_old(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0)
{
  computeDerived();
  Ttangent = Kpen;
}

LadrunoCohesiveHinge::LadrunoCohesiveHinge()
  : UniaxialMaterial(0, MAT_TAG_LadrunoCohesiveHinge),
    Mc(0.0), Gf(0.0), Kpen(0.0), userKpen(false), shape(EXP), penaltyRatio(1000.0),
    P0(0.0), alphaPM(0.5), implex(false), P_axial(0.0),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), CkappaMax_old(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0)
{
  // parameters arrive via recvSelf -> computeDerived
}

LadrunoCohesiveHinge::~LadrunoCohesiveHinge()
{
}

double LadrunoCohesiveHinge::getPMMultiplier(void) const
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

// ===========================================================================
//  calibration: enforce the guarded penalty floor, build kappa0/aExp/kappaf
//  so the monotonic envelope integrates to EXACTLY Gf.
// ===========================================================================
void LadrunoCohesiveHinge::computeDerived(void)
{
  if (Mc <= 0.0 || Gf <= 0.0) { kappa0 = aExp = kappaf = 0.0; return; }

  const double KpenFloor = Mc*Mc / (2.0*Gf);   // pre-peak energy == Gf at the floor

  // default penalty (!userKpen or Kpen<=0): penaltyRatio multiples above the floor
  if (!userKpen || Kpen <= 0.0) {
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
  }
  else if (shape == LINEAR) {
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
  if (P0 > 0.0) {
    double f = getPMMultiplier();
    Menv  *= f;
    dMenv *= f;
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

  if (implex) {
    // Oliver et al. (2008) IMPL-EX step extrapolation of the damage driver:
    // extrapolate kappaMax from the previous committed increment and evaluate the
    // secant stiffness explicitly so the step tangent is strictly positive-definite
    // and constant across global Newton iterations.
    double time_factor = 1.0;
    if (ops_Dt > 0.0 && dtime_n_commit > 0.0) {
      time_factor = ops_Dt / dtime_n_commit;
      if (time_factor > 2.0) time_factor = 2.0;
    } else if (implexCommittedSteps == 0) {
      time_factor = 0.0;
    }
    double dKappaMax_n = CkappaMax - CkappaMax_old;
    double kappaMax_tilde = CkappaMax + time_factor * dKappaMax_n;
    if (kappaMax_tilde < CkappaMax) kappaMax_tilde = CkappaMax;

    double Menv_tilde = 0.0, dMenv_tilde = 0.0;
    envelope(kappaMax_tilde, Menv_tilde, dMenv_tilde);
    double fpm = (P0 > 0.0) ? getPMMultiplier() : 1.0;
    double KpenEff = Kpen * fpm;
    double Ksec_tilde = (kappaMax_tilde > 0.0) ? (Menv_tilde / kappaMax_tilde) : KpenEff;
    double Kfloor = 1.0e-6 * KpenEff;
    if (Ksec_tilde < Kfloor) Ksec_tilde = Kfloor;

    Tstress  = Ksec_tilde * strain;
    Ttangent = Ksec_tilde;
    Twork    = Cwork + 0.5*(Cstress + Tstress)*(Tstrain - Cstrain);
    return 0;
  }

  double Menv, dMenv;
  envelope(TkappaMax, Menv, dMenv);

  // secant-to-origin (isotropic damage); Ksec = Menv(kappaMax)/kappaMax
  double KpenEff = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
  double Ksec = (TkappaMax > 0.0) ? (Menv / TkappaMax) : KpenEff;

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
  CkappaMax_old = CkappaMax;
  CkappaMax     = TkappaMax;
  if (implex) {
    // Commit the implicit stress at the true trial damage state for next-step anchor
    double Menv, dMenv;
    envelope(CkappaMax, Menv, dMenv);
    double KpenEff = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
    double Ksec = (CkappaMax > 0.0) ? (Menv / CkappaMax) : KpenEff;
    double stressImplicit = Ksec * Tstrain;
    Cwork   = Cwork + 0.5*(Cstress + stressImplicit)*(Tstrain - Cstrain);
    Cstress = stressImplicit;
    Twork   = Cwork;
  } else {
    Cstress = Tstress;
    Cwork   = Twork;
  }
  Cstrain        = Tstrain;
  dtime_n_commit = ops_Dt;
  implexCommittedSteps++;
  return 0;
}

int LadrunoCohesiveHinge::revertToLastCommit(void)
{
  return this->setTrialStrain(Cstrain);
}

int LadrunoCohesiveHinge::revertToStart(void)
{
  CkappaMax = 0.0; CkappaMax_old = 0.0;
  Cstrain = 0.0; Cstress = 0.0; Cwork = 0.0;
  dtime_n_commit = 0.0; implexCommittedSteps = 0;
  TkappaMax = 0.0; Tstrain = 0.0; Tstress = 0.0; Twork = 0.0;
  Ttangent = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
  return 0;
}

// ===========================================================================
//  copy
// ===========================================================================
UniaxialMaterial* LadrunoCohesiveHinge::getCopy(void)
{
  LadrunoCohesiveHinge* theCopy =
    new LadrunoCohesiveHinge(this->getTag(), Mc, Gf, userKpen ? Kpen : -1.0, shape, penaltyRatio, P0, alphaPM, implex);

  theCopy->userKpen             = userKpen;
  theCopy->P_axial              = P_axial;
  theCopy->CkappaMax            = CkappaMax;
  theCopy->CkappaMax_old        = CkappaMax_old;
  theCopy->Cstrain              = Cstrain;
  theCopy->Cstress              = Cstress;
  theCopy->Cwork                = Cwork;
  theCopy->dtime_n_commit       = dtime_n_commit;
  theCopy->implexCommittedSteps = implexCommittedSteps;
  theCopy->TkappaMax            = TkappaMax;
  theCopy->Tstrain              = Tstrain;
  theCopy->Tstress              = Tstress;
  theCopy->Ttangent             = Ttangent;
  theCopy->Twork                = Twork;
  return theCopy;
}

// ===========================================================================
//  parallel / database
// ===========================================================================
int LadrunoCohesiveHinge::sendSelf(int cTag, Channel& theChannel)
{
  // tag + 4 params (Mc,Gf,Kpen,shape) + penaltyRatio + 4 committed history + 10 extended
  static thread_local Vector data(20);
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
  data(c++) = P0;
  data(c++) = alphaPM;
  data(c++) = implex ? 1.0 : 0.0;
  data(c++) = P_axial;
  data(c++) = CkappaMax_old;
  data(c++) = dtime_n_commit;
  data(c++) = (double)implexCommittedSteps;
  data(c++) = userKpen ? 1.0 : 0.0;
  data(c++) = Tstress;
  data(c++) = Ttangent;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHinge::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoCohesiveHinge::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static thread_local Vector data(20);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHinge::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  Mc                   = data(c++);
  Gf                   = data(c++);
  Kpen                 = data(c++);
  shape                = (int)data(c++);
  penaltyRatio         = data(c++);
  CkappaMax            = data(c++);
  Cstrain              = data(c++);
  Cstress              = data(c++);
  Cwork                = data(c++);
  P0                   = data(c++);
  alphaPM              = data(c++);
  implex               = (data(c++) != 0.0);
  P_axial              = data(c++);
  CkappaMax_old        = data(c++);
  dtime_n_commit       = data(c++);
  implexCommittedSteps = (int)data(c++);
  userKpen             = (data(c++) != 0.0);
  double savedTstress  = data(c++);
  double savedTtangent = data(c++);

  computeDerived();           // Kpen already valid (was floor-checked on send side)
  if (implex) {
    Tstrain   = Cstrain;
    TkappaMax = CkappaMax;
    Twork     = Cwork;
    Tstress   = savedTstress;
    Ttangent  = savedTtangent;
  } else {
    this->revertToLastCommit(); // sync trial to committed
  }
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
  if (P0 > 0.0)
    s << "  P0     : " << P0 << "   alphaPM: " << alphaPM << endln;
  if (implex)
    s << "  implex : on" << endln;
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
  if (strcmp(a, "pmMultiplier") == 0 || strcmp(a, "fpm") == 0)
    return new MaterialResponse(this, 105, 0.0);
  if (strcmp(a, "implexError") == 0 || strcmp(a, "ImplexError") == 0)
    return new MaterialResponse(this, 106, 0.0);

  return UniaxialMaterial::setResponse(argv, argc, s);
}

double LadrunoCohesiveHinge::getImplexError(void) const
{
  if (!implex) return 0.0;
  double Menv_imp = 0.0, dMenv_imp = 0.0;
  envelope(TkappaMax, Menv_imp, dMenv_imp);
  double KpenEff = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
  double Ksec_imp = (TkappaMax > 0.0) ? (Menv_imp / TkappaMax) : KpenEff;
  double M_imp = Ksec_imp * Tstrain;
  return fabs(Tstress - M_imp) / (fabs(M_imp) + 1.0e-12);
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
      double KpenEff = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
      double Ksec = (TkappaMax > 0.0) ? (Menv/TkappaMax) : KpenEff;
      matInfo.theDouble = (KpenEff > 0.0) ? (1.0 - Ksec/KpenEff) : 0.0;
      return 0;
    }
    case 104: {                                           // fully-open (broken) flag
      double open = 0.0;
      if (shape == LINEAR && TkappaMax >= kappa0 + kappaf) open = 1.0;
      matInfo.theDouble = open;
      return 0;
    }
    case 105:
      matInfo.theDouble = getPMMultiplier();
      return 0;
    case 106:
      matInfo.theDouble = getImplexError();
      return 0;
    default:
      return UniaxialMaterial::getResponse(responseID, matInfo);
  }
}

int LadrunoCohesiveHinge::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "Mc") == 0)      return param.addObject(1, this);
  if (strcmp(argv[0], "Gf") == 0)      return param.addObject(2, this);
  if (strcmp(argv[0], "Kpen") == 0)    return param.addObject(3, this);
  if (strcmp(argv[0], "P0") == 0)      return param.addObject(4, this);
  if (strcmp(argv[0], "alphaPM") == 0) return param.addObject(5, this);
  return -1;
}

int LadrunoCohesiveHinge::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 1: Mc = info.theDouble; computeDerived(); return 0;
    case 2: Gf = info.theDouble; computeDerived(); return 0;
    case 3: Kpen = info.theDouble; userKpen = (Kpen > 0.0); computeDerived(); return 0;
    case 4: P0 = info.theDouble; return 0;
    case 5: alphaPM = info.theDouble; return 0;
    default: return -1;
  }
}

