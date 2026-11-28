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
//        <-pm P0 <alpha>> <-implex> <-implexDt step|pseudo>
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
           << "<-penalty K? | -penaltyRatio r?> <-exp | -linear> <-pm P0? <alpha?>> <-implex> <-implexDt step|pseudo>\n";
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
  double P0t = 0.0;
  double Pbal = 0.0;
  double muBal = 1.0;
  bool   implex = false;
  int    implexClock = LadrunoImplexClock::STEP;

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
      P0t = 0.0;
      Pbal = 0.0;
      double optAlpha = alphaPM;
      if (lchTryOptionalDouble(optAlpha)) {
        alphaPM = optAlpha;
      }
    }
    else if (strcmp(flag, "-pmBal") == 0 || strcmp(flag, "-PMBal") == 0) {
      double pmb[4];
      numData = 4;
      if (OPS_GetNumRemainingInputArgs() < 4 || OPS_GetDoubleInput(&numData, pmb) < 0 ||
          pmb[0] <= 0.0 || pmb[1] <= 0.0 || pmb[2] <= 0.0 || pmb[2] >= pmb[0] || pmb[3] <= 0.0) {
        opserr << "WARNING LadrunoCohesiveHinge tag " << tag
               << ": -pmBal needs P0c>0 P0t>0 0<Pbal<P0c muBal>0\n";
        return 0;
      }
      P0    = pmb[0];
      P0t   = pmb[1];
      Pbal  = pmb[2];
      muBal = pmb[3];
    }
    else if (strcmp(flag, "-implex") == 0) {
      implex = true;
    }
    else if (strcmp(flag, "-implexDt") == 0) {   // Ladruno WP-217
      const char* src = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if (!LadrunoImplexClock::parse(src, implexClock)) {
        opserr << "WARNING LadrunoCohesiveHinge tag " << tag << ": -implexDt needs step or pseudo\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoCohesiveHinge tag " << tag
             << ": unknown option '" << flag << "'\n";
      return 0;
    }
  }

  LadrunoCohesiveHinge* mat = new LadrunoCohesiveHinge(tag, Mc, Gf, Kpen, shape, penaltyRatio,
                                                       P0, alphaPM, implex, P0t, Pbal, muBal);
  mat->setImplexClock(implexClock);
  return mat;
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoCohesiveHinge::LadrunoCohesiveHinge(int tag, double mc, double gf,
                                           double kpen, int shp, double pratio,
                                           double p0, double apm, bool impl,
                                           double p0t, double pbal, double mubal)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoCohesiveHinge),
    Mc(mc), Gf(gf), Kpen(kpen), userKpen(kpen > 0.0), shape(shp), penaltyRatio(pratio),
    P0(p0), alphaPM(apm), P0t(p0t), Pbal(pbal), muBal(mubal), implex(impl), P_axial(0.0),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), CkappaMax_old(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0),
    parameterID(0)
{
  computeDerived();
  Ttangent = Kpen;
}

LadrunoCohesiveHinge::LadrunoCohesiveHinge()
  : UniaxialMaterial(0, MAT_TAG_LadrunoCohesiveHinge),
    Mc(0.0), Gf(0.0), Kpen(0.0), userKpen(false), shape(EXP), penaltyRatio(1000.0),
    P0(0.0), alphaPM(0.5), P0t(0.0), Pbal(0.0), muBal(1.0), implex(false), P_axial(0.0),
    kappa0(0.0), aExp(0.0), kappaf(0.0),
    CkappaMax(0.0), CkappaMax_old(0.0), Cstrain(0.0), Cstress(0.0), Cwork(0.0),
    dtime_n_commit(0.0), implexCommittedSteps(0),
    Tstrain(0.0), Tstress(0.0), Ttangent(0.0), TkappaMax(0.0), Twork(0.0),
    parameterID(0)
{
  // parameters arrive via recvSelf -> computeDerived
}

LadrunoCohesiveHinge::~LadrunoCohesiveHinge()
{
}

double LadrunoCohesiveHinge::getPMMultiplier(void) const
{
  // WP-267: the law lives in LadrunoCohesiveHingeKernel.h (pmMultiplier); this is
  // its double instantiation (bit-identity gate: tests/ladruno_cohesivehinge_dual_check.cpp).
  ladruno_cohesive_hinge::ParamsT<double> p;
  fillKernelParams(p);
  return ladruno_cohesive_hinge::pmMultiplier(p);
}

double LadrunoCohesiveHinge::getPMDerivative(void) const
{
  if (P0 <= 0.0) return 0.0;
  if (Pbal > 0.0 && P0 > Pbal && P0t > 0.0) {
    double f = 0.0, df = 0.0;
    if (P_axial >= 0.0 && P_axial <= P0t) {
      f  = 1.0 - P_axial / P0t;
      df = -1.0 / P0t;
    } else if (P_axial < 0.0 && P_axial >= -Pbal) {
      double u = (P_axial + Pbal) / Pbal;
      f  = 1.0 + (muBal - 1.0) * (1.0 - u * u);
      df = -2.0 * (muBal - 1.0) * (P_axial + Pbal) / (Pbal * Pbal);
    } else if (P_axial < -Pbal && P_axial >= -P0) {
      double denom = P0 - Pbal;
      double u = (P_axial + Pbal) / denom;
      f  = muBal * (1.0 - u * u);
      df = -2.0 * muBal * (P_axial + Pbal) / (denom * denom);
    }
    if (f <= 0.05 || f >= 3.0) return 0.0;
    return df;
  }
  double n = P_axial / P0;
  double f = 1.0 - alphaPM * n - (1.0 - alphaPM) * n * n;
  if (f <= 0.05 || f >= 2.5) return 0.0;
  return (-alphaPM - 2.0 * (1.0 - alphaPM) * n) / P0;
}

// ===========================================================================
//  calibration: enforce the guarded penalty floor, build kappa0/aExp/kappaf
//  so the monotonic envelope integrates to EXACTLY Gf.
// ===========================================================================
void LadrunoCohesiveHinge::computeDerived(void)
{
  // WP-267: LadrunoCohesiveHingeKernel.h::computeDerived (double instantiation)
  ladruno_cohesive_hinge::ParamsT<double> p;
  fillKernelParams(p);
  double KpenIn = 0.0, KpenFloor = 0.0, KpenSafe = 0.0;
  if (ladruno_cohesive_hinge::computeDerived(p, &KpenIn, &KpenFloor, &KpenSafe))
    opserr << "WARNING LadrunoCohesiveHinge tag " << this->getTag()
           << ": penalty Kpen=" << KpenIn << " <= floor Mc^2/(2Gf)=" << KpenFloor
           << " (pre-peak energy would exceed Gf); clamping to " << KpenSafe << "\n";
  Kpen   = p.Kpen;
  kappa0 = p.kappa0;
  aExp   = p.aExp;
  kappaf = p.kappaf;
}

void LadrunoCohesiveHinge::fillKernelParams(ladruno_cohesive_hinge::ParamsT<double>& p) const
{
  p.Mc = Mc; p.Gf = Gf; p.Kpen = Kpen; p.userKpen = userKpen;
  p.shape = shape; p.penaltyRatio = penaltyRatio;
  p.P0 = P0; p.alphaPM = alphaPM; p.P0t = P0t; p.Pbal = Pbal; p.muBal = muBal;
  p.P_axial = P_axial;
  p.kappa0 = kappa0; p.aExp = aExp; p.kappaf = kappaf;
}

// ===========================================================================
//  envelope moment Menv(km) and consistent slope dMenv/dkm at km >= 0
// ===========================================================================
void LadrunoCohesiveHinge::envelope(double km, double& Menv, double& dMenv) const
{
  ladruno_cohesive_hinge::ParamsT<double> p;
  fillKernelParams(p);
  ladruno_cohesive_hinge::envelope(p, km, Menv, dMenv);
}

// ===========================================================================
//  trial state
// ===========================================================================
int LadrunoCohesiveHinge::setTrialStrain(double strain, double)
{
  Tstrain = strain;

  // WP-267: the trial state lives in the templated kernel
  // (LadrunoCohesiveHingeKernel.h::trialState); this is its double instantiation,
  // which performs exactly the operations of the former in-class code.
  // IMPL-EX (Oliver et al. 2008): the step ratio comes from LadrunoImplexClock
  // (WP-217): frozen per step, never the live ops_Dt, which under
  // DisplacementControl is the iterate's load-factor increment.
  ladruno_cohesive_hinge::ParamsT<double> p;
  fillKernelParams(p);
  ladruno_cohesive_hinge::StateT<double> C;
  C.kappaMax = CkappaMax;
  C.kappaMaxOld = CkappaMax_old;
  double time_factor = 0.0;
  if (implex)
    time_factor = (implexCommittedSteps == 0) ? 0.0 : implexClock.factor(dtime_n_commit);
  ladruno_cohesive_hinge::trialState(p, implex, time_factor, strain, C,
                                     TkappaMax, Tstress, Ttangent);

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
  dtime_n_commit = implexClock.commitDt();
  implexCommittedSteps++;
  return 0;
}

int LadrunoCohesiveHinge::revertToLastCommit(void)
{
  int res = this->setTrialStrain(Cstrain);
  implexClock.reopen();   // the retried step latches its own dt
  return res;
}

int LadrunoCohesiveHinge::revertToStart(void)
{
  SHV.clear();                                   // WP-267: history sensitivities
  CkappaMax = 0.0; CkappaMax_old = 0.0;
  Cstrain = 0.0; Cstress = 0.0; Cwork = 0.0;
  dtime_n_commit = 0.0; implexCommittedSteps = 0; implexClock.reopen();
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
    new LadrunoCohesiveHinge(this->getTag(), Mc, Gf, userKpen ? Kpen : -1.0, shape, penaltyRatio,
                             P0, alphaPM, implex, P0t, Pbal, muBal);

  theCopy->userKpen             = userKpen;
  theCopy->P_axial              = P_axial;
  theCopy->CkappaMax            = CkappaMax;
  theCopy->CkappaMax_old        = CkappaMax_old;
  theCopy->Cstrain              = Cstrain;
  theCopy->Cstress              = Cstress;
  theCopy->Cwork                = Cwork;
  theCopy->dtime_n_commit       = dtime_n_commit;
  theCopy->implexCommittedSteps = implexCommittedSteps;
  theCopy->implexClock          = implexClock;
  theCopy->TkappaMax            = TkappaMax;
  theCopy->Tstrain              = Tstrain;
  theCopy->Tstress              = Tstress;
  theCopy->Ttangent             = Ttangent;
  theCopy->Twork                = Twork;
  theCopy->parameterID          = parameterID;   // WP-267: DDM
  theCopy->SHV                  = SHV;
  return theCopy;
}

// ===========================================================================
//  parallel / database
// ===========================================================================
int LadrunoCohesiveHinge::sendSelf(int cTag, Channel& theChannel)
{
  // tag + 4 params (Mc,Gf,Kpen,shape) + penaltyRatio + 4 committed history + 13 extended
  static thread_local Vector data(23);
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
  data(c++) = LadrunoImplexClock::pack(implex, implexClock.mode);   // WP-217: bit 1 = clock
  data(c++) = P_axial;
  data(c++) = CkappaMax_old;
  data(c++) = dtime_n_commit;
  data(c++) = (double)implexCommittedSteps;
  data(c++) = userKpen ? 1.0 : 0.0;
  data(c++) = Tstress;
  data(c++) = Ttangent;
  data(c++) = P0t;
  data(c++) = Pbal;
  data(c++) = muBal;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoCohesiveHinge::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoCohesiveHinge::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static thread_local Vector data(23);
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
  implexClock.mode     = LadrunoImplexClock::unpackMode(data(c));
  implex               = LadrunoImplexClock::unpackImplex(data(c++));
  P_axial              = data(c++);
  CkappaMax_old        = data(c++);
  dtime_n_commit       = data(c++);
  implexCommittedSteps = (int)data(c++);
  userKpen             = (data(c++) != 0.0);
  double savedTstress  = data(c++);
  double savedTtangent = data(c++);
  P0t                  = data(c++);
  Pbal                 = data(c++);
  muBal                = data(c++);

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
  if (P0 > 0.0 && Pbal > 0.0)
    s << "  pmBal  : P0c=" << P0 << " P0t=" << P0t << " Pbal=" << Pbal << " muBal=" << muBal << endln;
  else if (P0 > 0.0)
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
  if (strcmp(a, "pmDerivative") == 0 || strcmp(a, "dfpm") == 0)
    return new MaterialResponse(this, 107, 0.0);

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
    case 107:
      matInfo.theDouble = getPMDerivative();
      return 0;
    default:
      return UniaxialMaterial::getResponse(responseID, matInfo);
  }
}

int LadrunoCohesiveHinge::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "Mc") == 0)      { param.setValue(Mc);      return param.addObject(1, this); }
  if (strcmp(argv[0], "Gf") == 0)      { param.setValue(Gf);      return param.addObject(2, this); }
  if (strcmp(argv[0], "Kpen") == 0)    { param.setValue(Kpen);    return param.addObject(3, this); }
  if (strcmp(argv[0], "P0") == 0)      { param.setValue(P0);      return param.addObject(4, this); }
  if (strcmp(argv[0], "alphaPM") == 0) { param.setValue(alphaPM); return param.addObject(5, this); }
  if (strcmp(argv[0], "P0t") == 0)     { param.setValue(P0t);     return param.addObject(6, this); }
  if (strcmp(argv[0], "Pbal") == 0)    { param.setValue(Pbal);    return param.addObject(7, this); }
  if (strcmp(argv[0], "muBal") == 0)   { param.setValue(muBal);   return param.addObject(8, this); }
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
    case 6: P0t = info.theDouble; return 0;
    case 7: Pbal = info.theDouble; return 0;
    case 8: muBal = info.theDouble; return 0;
    default: return -1;
  }
}

int LadrunoCohesiveHinge::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

// ===========================================================================
//  WP-267: DDM sensitivities by dual numbers (the WP-262 LadrunoUniaxialJ2 recipe)
//
//  OpenSees DDM contract (computeSensitivities BEFORE commit):
//  getStressSensitivity(g) is dM/dtheta at the CURRENT trial jump held fixed,
//  given the committed history and its stored sensitivities;
//  commitSensitivity(d jump/d theta, g, n) advances the history sensitivities
//  with the converged jump derivative. Both are one Dual<1> pass of the kernel.
//
//  The seed is the active parameter (parameterID: 1 Mc, 2 Gf, 3 Kpen, 4 P0,
//  5 alphaPM, 6 P0t, 7 Pbal, 8 muBal). The derived constants (kappa0, aExp,
//  kappaf and a default Kpen = penaltyRatio Mc^2/(2Gf)) are re-derived in Dual,
//  exactly as updateParameter re-derives them, so the result is the derivative
//  of what updateParameter + the run would do. Seeding Kpen marks it user-set,
//  as updateParameter(3) does. The history kappaMax and kappaMax_old (IMPL-EX
//  extrapolation) carry their stored sensitivities. Not differentiated: P_axial
//  (the parent beam's axial force: its sensitivity is the element's business).
// ===========================================================================
void LadrunoCohesiveHinge::seededParams(ladruno_cohesive_hinge::ParamsT<ladruno_dual::Dual<1> >& p) const
{
  ladruno_cohesive_hinge::ParamsT<double> pd;
  fillKernelParams(pd);
  p.Mc = pd.Mc; p.Gf = pd.Gf; p.Kpen = pd.Kpen; p.userKpen = pd.userKpen;
  p.shape = pd.shape; p.penaltyRatio = pd.penaltyRatio;
  p.P0 = pd.P0; p.alphaPM = pd.alphaPM; p.P0t = pd.P0t; p.Pbal = pd.Pbal; p.muBal = pd.muBal;
  p.P_axial = pd.P_axial;
  p.kappa0 = pd.kappa0; p.aExp = pd.aExp; p.kappaf = pd.kappaf;
  switch (parameterID) {
    case 1: p.Mc.d[0]      = 1.0; break;
    case 2: p.Gf.d[0]      = 1.0; break;
    case 3: p.Kpen.d[0]    = 1.0; p.userKpen = true; break;
    case 4: p.P0.d[0]      = 1.0; break;
    case 5: p.alphaPM.d[0] = 1.0; break;
    case 6: p.P0t.d[0]     = 1.0; break;
    case 7: p.Pbal.d[0]    = 1.0; break;
    case 8: p.muBal.d[0]   = 1.0; break;
    default: break;                        // not this material's parameter
  }
  ladruno_cohesive_hinge::computeDerived(p);   // values unchanged; derivatives of the derived constants
}

void LadrunoCohesiveHinge::sensitivityPass(int gradIndex, double strainGradient,
                                           double& dStress, double& dTangent, double* dTrialState)
{
  typedef ladruno_dual::Dual<1> D1;
  ladruno_cohesive_hinge::ParamsT<D1> p;
  seededParams(p);

  ladruno_cohesive_hinge::StateT<D1> C;
  C.kappaMax = CkappaMax;
  C.kappaMaxOld = CkappaMax_old;
  const size_t base = (size_t)gradIndex * SHV_N;
  if (gradIndex >= 0 && SHV.size() >= base + SHV_N) {
    C.kappaMax.d[0]    = SHV[base];
    C.kappaMaxOld.d[0] = SHV[base + 1];
  }

  double time_factor = 0.0;
  if (implex)
    time_factor = (implexCommittedSteps == 0) ? 0.0 : implexClock.factor(dtime_n_commit);

  D1 strain(Tstrain);
  strain.d[0] = strainGradient;
  D1 kmax, stress, tangent;
  ladruno_cohesive_hinge::trialState(p, implex, time_factor, strain, C, kmax, stress, tangent);

  dStress  = stress.d[0];
  dTangent = tangent.d[0];
  if (dTrialState) {
    dTrialState[0] = kmax.d[0];            // commitState: CkappaMax     = TkappaMax
    dTrialState[1] = C.kappaMax.d[0];      //              CkappaMax_old = CkappaMax
  }
}

double LadrunoCohesiveHinge::getStressSensitivity(int gradIndex, bool)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dS;
}

double LadrunoCohesiveHinge::getTangentSensitivity(int gradIndex)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dT;
}

double LadrunoCohesiveHinge::getInitialTangentSensitivity(int)
{
  ladruno_cohesive_hinge::ParamsT<ladruno_dual::Dual<1> > p;
  seededParams(p);
  return ladruno_cohesive_hinge::initialTangent(p).d[0];
}

int LadrunoCohesiveHinge::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
  if (gradIndex < 0 || numGrads <= gradIndex) return -1;
  if (SHV.size() < (size_t)numGrads * SHV_N) SHV.resize((size_t)numGrads * SHV_N, 0.0);
  double dTr[SHV_N];
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, strainGradient, dS, dT, dTr);
  for (int i = 0; i < SHV_N; i++) SHV[(size_t)gradIndex * SHV_N + i] = dTr[i];
  return 0;
}
