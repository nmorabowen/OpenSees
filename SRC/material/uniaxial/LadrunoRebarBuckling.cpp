/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: rebar-buckling WRAPPER UniaxialMaterial (Dhakal-Maekawa 2002).
// See LadrunoRebarBuckling.h and Ladruno_implementation/14_ladruno_rebar_buckling_adr.md.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoRebarBuckling.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

const double LadrunoRebarBuckling::LSR_MIN = 1.0e-12;

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoRebarBuckling tag matTag -lsr s/d
//        <-fy fy> <-E E> <-model dm|ga> <-alpha a>
// ===========================================================================
void* OPS_LadrunoRebarBuckling(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoRebarBuckling tag? matTag? "
           << "-lsr lsr? <-fy fy?> <-E E?> <-model dm|ga> <-alpha a?>\n";
    return 0;
  }

  int idata[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, idata) < 0) {
    opserr << "WARNING invalid LadrunoRebarBuckling tag/matTag\n";
    return 0;
  }

  double lsr = 0.0, fy = 0.0, E = 0.0, alpha = 1.0;
  int model = LadrunoRebarBuckling::MODEL_DM;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    numData = 1;
    if (strcmp(flag, "-lsr") == 0) {
      if (OPS_GetDoubleInput(&numData, &lsr) < 0) { opserr << "WARNING LadrunoRebarBuckling: -lsr wants a value\n"; return 0; }
    }
    else if (strcmp(flag, "-fy") == 0) {
      if (OPS_GetDoubleInput(&numData, &fy) < 0) { opserr << "WARNING LadrunoRebarBuckling: -fy wants a value\n"; return 0; }
    }
    else if (strcmp(flag, "-E") == 0) {
      if (OPS_GetDoubleInput(&numData, &E) < 0) { opserr << "WARNING LadrunoRebarBuckling: -E wants a value\n"; return 0; }
    }
    else if (strcmp(flag, "-alpha") == 0) {
      if (OPS_GetDoubleInput(&numData, &alpha) < 0) { opserr << "WARNING LadrunoRebarBuckling: -alpha wants a value\n"; return 0; }
    }
    else if (strcmp(flag, "-model") == 0) {
      const char* m = OPS_GetString();
      if (strcmp(m, "dm") == 0 || strcmp(m, "DM") == 0) model = LadrunoRebarBuckling::MODEL_DM;
      else if (strcmp(m, "ga") == 0 || strcmp(m, "GA") == 0) model = LadrunoRebarBuckling::MODEL_GA;
      else { opserr << "WARNING LadrunoRebarBuckling: -model wants dm|ga\n"; return 0; }
    }
    else {
      opserr << "WARNING LadrunoRebarBuckling: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* bar = OPS_getUniaxialMaterial(idata[1]);
  if (bar == 0) {
    opserr << "WARNING LadrunoRebarBuckling: wrapped material " << idata[1] << " not found\n";
    return 0;
  }

  if (model == LadrunoRebarBuckling::MODEL_GA) {
    opserr << "WARNING LadrunoRebarBuckling: -model ga not implemented in v1 (use dm)\n";
    return 0;
  }
  if (lsr > LadrunoRebarBuckling::LSR_MIN && fy <= 0.0) {
    opserr << "WARNING LadrunoRebarBuckling: -fy required (>0) when -lsr enables buckling\n";
    return 0;
  }

  UniaxialMaterial* mat =
    new LadrunoRebarBuckling(idata[0], *bar, lsr, fy, E, model, alpha);
  if (mat == 0) {
    opserr << "WARNING LadrunoRebarBuckling: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoRebarBuckling::LadrunoRebarBuckling(int tag, UniaxialMaterial& bar,
                                           double lsr_, double fy_, double E_,
                                           int model_, double alpha_)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), lsr(lsr_), fy(fy_), E(E_), model(model_), alpha(alpha_),
    CmaxTenStrain(0.0), CfStarL(0.0), Conset(0),
    TmaxTenStrain(0.0), TfStarL(0.0), Tonset(0),
    Tstress(0.0), Ttangent(0.0)
{
  theBar = bar.getCopy();
  if (theBar == 0)
    opserr << "LadrunoRebarBuckling::LadrunoRebarBuckling - failed to copy wrapped material\n";
  if (E <= 0.0 && theBar != 0)
    E = theBar->getInitialTangent();
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : 0.0;
}

LadrunoRebarBuckling::LadrunoRebarBuckling()
  : UniaxialMaterial(0, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), lsr(0.0), fy(0.0), E(0.0), model(MODEL_DM), alpha(1.0),
    CmaxTenStrain(0.0), CfStarL(0.0), Conset(0),
    TmaxTenStrain(0.0), TfStarL(0.0), Tonset(0),
    Tstress(0.0), Ttangent(0.0)
{

}

LadrunoRebarBuckling::~LadrunoRebarBuckling()
{
  if (theBar != 0)
    delete theBar;
}

// ===========================================================================
//  Dhakal-Maekawa buckled-average stress (+ analytic consistent tangent)
//   e = eps - e_cross (relative compressive strain). Mirrors
//   ReinforcingSteel::Buckled_stress_Dhakal; fStarL is the bare backbone stress
//   anchored at the onset strain eStar.
// ===========================================================================
double LadrunoRebarBuckling::eStarStrain(void) const
{
  double eyp = fy / E;
  double f = 55.0 - 2.3 * sqrt(fy / E * 2000.0) * lsr;
  if (f < 7.0) f = 7.0;
  return -f * eyp;                 // compressive onset strain (relative to e_cross)
}

double LadrunoRebarBuckling::dmStress(double e, double sBare, double kBare,
                                      double fStarL, double& tanOut) const
{
  double eyp   = fy / E;
  double eStar = eStarStrain();

  if (e >= -eyp || fStarL == 0.0) {     // not buckled (or no anchor yet)
    tanOut = kBare;
    return sBare;
  }

  double fStar = fStarL * alpha * (1.1 - 0.016 * sqrt(fy / E * 2000.0) * lsr);
  if (fStar > -0.2 * fy) fStar = -0.2 * fy;     // residual floor

  double red, dRedDe;
  if (e >= eStar) {                     // intermediate branch  [eStar, -eyp)
    red    = 1.0 - (1.0 - fStar / fStarL) * (e + eyp) / (eStar + eyp);
    dRedDe = -(1.0 - fStar / fStarL) / (eStar + eyp);
  } else {                              // deep branch  e < eStar
    red    = (fStar - 0.02 * E * (e - eStar)) / fStarL;
    dRedDe = -0.02 * E / fStarL;
  }

  double ave = red * sBare;
  if (ave > -0.2 * fy) {                // residual-stress floor (constant => zero tangent)
    tanOut = 0.0;
    return -0.2 * fy;
  }
  // analytic consistent tangent: d(red*sBare)/deps = red*kBare + sBare*(dRed/de)
  tanOut = red * kBare + sBare * dRedDe;
  return ave;
}

// ===========================================================================
//  strain interface
// ===========================================================================
int LadrunoRebarBuckling::setTrialStrain(double strain, double rate)
{
  theBar->setTrialStrain(strain, rate);
  double sBare = theBar->getStress();
  double kBare = theBar->getTangent();

  // track the tensile-strain anchor e_cross = maxTenStrain - sigma(maxTen)/E
  TmaxTenStrain = CmaxTenStrain;
  double maxTenStress = 0.0;
  if (strain > TmaxTenStrain) { TmaxTenStrain = strain; maxTenStress = (sBare > 0.0 ? sBare : 0.0); }
  double e_cross = (TmaxTenStrain > 0.0) ? (TmaxTenStrain - maxTenStress / E) : 0.0;

  double eyp = fy / E;
  double e   = strain - e_cross;

  if (lsr <= LSR_MIN || E <= 0.0 || fy <= 0.0 || e >= -eyp) {
    // pass-through (tension / pre-onset / buckling disabled): bare response, reset onset
    Tonset = 0;
    TfStarL = 0.0;
    Tstress = sBare;
    Ttangent = kBare;
    return 0;
  }

  // --- buckling branch ---
  // anchor fStarL = bare backbone stress at eStar, once per compression excursion,
  // by a TRIAL probe (setTrialStrain is memoryless given committed state, so the
  // restore is exact).
  if (Conset == 0 && Tonset == 0) {
    double eStarTot = eStarStrain() + e_cross;
    theBar->setTrialStrain(eStarTot, rate);
    TfStarL = theBar->getStress();
    theBar->setTrialStrain(strain, rate);          // restore
    sBare = theBar->getStress();
    kBare = theBar->getTangent();
    Tonset = 1;
  } else {
    Tonset = 1;
    TfStarL = (Conset != 0) ? CfStarL : TfStarL;
  }

  double tan;
  Tstress = dmStress(e, sBare, kBare, TfStarL, tan);
  Ttangent = tan;
  return 0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRebarBuckling::commitState(void)
{
  CmaxTenStrain = TmaxTenStrain;
  CfStarL = TfStarL;
  Conset  = Tonset;
  return theBar->commitState();
}

int LadrunoRebarBuckling::revertToLastCommit(void)
{
  TmaxTenStrain = CmaxTenStrain;
  TfStarL = CfStarL;
  Tonset  = Conset;
  return theBar->revertToLastCommit();
}

int LadrunoRebarBuckling::revertToStart(void)
{
  CmaxTenStrain = 0.0; CfStarL = 0.0; Conset = 0;
  TmaxTenStrain = 0.0; TfStarL = 0.0; Tonset = 0;
  Tstress = 0.0;
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : 0.0;
  return (theBar != 0) ? theBar->revertToStart() : 0;
}

// ===========================================================================
//  copy
// ===========================================================================
UniaxialMaterial* LadrunoRebarBuckling::getCopy(void)
{
  LadrunoRebarBuckling* theCopy =
    new LadrunoRebarBuckling(this->getTag(), *theBar, lsr, fy, E, model, alpha);
  theCopy->CmaxTenStrain = CmaxTenStrain;
  theCopy->CfStarL = CfStarL;
  theCopy->Conset  = Conset;
  theCopy->TmaxTenStrain = TmaxTenStrain;
  theCopy->TfStarL = TfStarL;
  theCopy->Tonset  = Tonset;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  return theCopy;
}

// ===========================================================================
//  parallel (nested-material idiom, mirrors FatigueMaterial)
// ===========================================================================
int LadrunoRebarBuckling::sendSelf(int cTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theBar->getClassTag();
  int matDbTag = theBar->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theBar->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send ID\n";
    return -1;
  }

  static Vector data(8);
  data(0) = lsr;
  data(1) = fy;
  data(2) = E;
  data(3) = model;
  data(4) = alpha;
  data(5) = CmaxTenStrain;
  data(6) = CfStarL;
  data(7) = Conset;
  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send data\n";
    return -2;
  }

  if (theBar->sendSelf(cTag, theChannel) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send wrapped material\n";
    return -3;
  }
  return 0;
}

int LadrunoRebarBuckling::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv ID\n";
    return -1;
  }
  this->setTag((int)dataID(0));

  int matClassTag = dataID(1);
  if (theBar == 0 || theBar->getClassTag() != matClassTag) {
    if (theBar != 0) delete theBar;
    theBar = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theBar == 0) {
      opserr << "LadrunoRebarBuckling::recvSelf - failed to create wrapped material\n";
      return -2;
    }
  }
  theBar->setDbTag(dataID(2));

  static Vector data(8);
  if (theChannel.recvVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv data\n";
    return -3;
  }
  lsr   = data(0);
  fy    = data(1);
  E     = data(2);
  model = (int)data(3);
  alpha = data(4);
  CmaxTenStrain = data(5);
  CfStarL = data(6);
  Conset  = (int)data(7);

  if (theBar->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv wrapped material\n";
    return -4;
  }
  this->revertToLastCommit();
  Tstress = theBar->getStress();
  Ttangent = theBar->getTangent();
  return 0;
}

// ===========================================================================
//  print / responses
// ===========================================================================
void LadrunoRebarBuckling::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoRebarBuckling\", ";
    s << "\"lsr\": " << lsr << ", \"fy\": " << fy << ", \"E\": " << E << ", ";
    s << "\"model\": \"" << (model == MODEL_GA ? "ga" : "dm") << "\", ";
    s << "\"material\": \"" << (theBar ? theBar->getTag() : -1) << "\"}";
    return;
  }
  s << "LadrunoRebarBuckling (Dhakal-Maekawa rebar-buckling wrapper)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  lsr    : " << lsr << "  fy: " << fy << "  E: " << E
    << "  model: " << (model == MODEL_GA ? "ga" : "dm") << endln;
  if (theBar) { s << "  wraps  : "; theBar->Print(s, flag); }
}

Response* LadrunoRebarBuckling::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc > 0 && (strcmp(argv[0], "buckled") == 0 || strcmp(argv[0], "isBuckled") == 0))
    return new MaterialResponse(this, 100, 0.0);
  // delegate everything else (stress/strain/tangent + wrapped-material responses)
  if (argc > 1 && strcmp(argv[0], "material") == 0 && theBar != 0)
    return theBar->setResponse(&argv[1], argc - 1, s);
  return UniaxialMaterial::setResponse(argv, argc, s);
}

int LadrunoRebarBuckling::getResponse(int responseID, Information& matInfo)
{
  if (responseID == 100) { matInfo.theDouble = (double)Tonset; return 0; }
  return UniaxialMaterial::getResponse(responseID, matInfo);
}
