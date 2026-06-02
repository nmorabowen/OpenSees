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

// Ladruno: rebar-buckling wrapper UniaxialMaterial. See LadrunoRebarBuckling.h
// and Ladruno_implementation/14_ladruno_rebar_buckling_adr.md for the model.
//
// The Dhakal-Maekawa buckled-average degradation is ported verbatim from
// ReinforcingSteel::Buckled_stress_Dhakal(ess, fss) but works in engineering
// strain on a generic wrapped bar (ReinforcingSteel works in log strain on its
// own monolithic backbone), and re-implements the cyclic bookkeeping with its
// own minimal state (an ε_max anchor) instead of the ReinforcingSteel
// branch-state machine. The consistent tangent is analytic (product rule)
// rather than finite-differenced.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoRebarBuckling.h>
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

static const double LRB_INVALID = 1.0e30;   // "fStarL not yet computed" sentinel

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoRebarBuckling tag matTag -lsr s/d
//                    <-model dm|ga> <-alpha a> <-fy fy> <-E E>
// ===========================================================================
void* OPS_LadrunoRebarBuckling(void)
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoRebarBuckling tag? matTag? -lsr s/d? "
           << "<-model dm|ga> <-alpha a?> <-fy fy?> <-E E?>\n";
    return 0;
  }

  int idata[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, idata) < 0) {
    opserr << "WARNING invalid LadrunoRebarBuckling tag/matTag\n";
    return 0;
  }

  double lsr   = 0.0;          // off by default (identity gate)
  double alpha = 1.0;          // DM residual-shape factor (ReinforcingSteel beta)
  double fy    = 0.0;          // queried/required when lsr>0
  double E     = 0.0;          // <=0 => take theBar->getInitialTangent()
  int    model = LadrunoRebarBuckling::MODEL_DM;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-lsr") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &lsr) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -lsr wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-alpha") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &alpha) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -alpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-fy") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &fy) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -fy wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-E") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &E) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -E wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-model") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "dm") == 0 || strcmp(law, "DM") == 0) {
        model = LadrunoRebarBuckling::MODEL_DM;
      } else if (strcmp(law, "ga") == 0 || strcmp(law, "GA") == 0) {
        opserr << "WARNING LadrunoRebarBuckling: -model ga (Gomes-Appleton) "
               << "is not implemented in v1; use -model dm\n";
        return 0;
      } else {
        opserr << "WARNING LadrunoRebarBuckling: unknown -model '" << law
               << "' (only 'dm' in v1)\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoRebarBuckling: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* bar = OPS_getUniaxialMaterial(idata[1]);
  if (bar == 0) {
    opserr << "WARNING LadrunoRebarBuckling: wrapped material " << idata[1]
           << " does not exist (for material " << idata[0] << ")\n";
    return 0;
  }

  if (E <= 0.0)
    E = bar->getInitialTangent();

  if (lsr > 0.0 && (fy <= 0.0 || E <= 0.0)) {
    opserr << "WARNING LadrunoRebarBuckling: with -lsr>0 a positive -fy is "
           << "required (and E from -E or the wrapped material's initial "
           << "tangent must be >0); material " << idata[0] << "\n";
    return 0;
  }

  UniaxialMaterial* theMat =
    new LadrunoRebarBuckling(idata[0], *bar, lsr, alpha, fy, E, model);
  if (theMat == 0) {
    opserr << "WARNING LadrunoRebarBuckling: failed to allocate material\n";
    return 0;
  }
  return theMat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoRebarBuckling::LadrunoRebarBuckling(int tag, UniaxialMaterial& bar,
                                           double lsr_, double alpha_,
                                           double fy_, double E_, int model_)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), model(model_), lsr(lsr_), alpha(alpha_), fy(fy_), E(E_),
    parameterID(0)
{
  theBar = bar.getCopy();
  if (theBar == 0) {
    opserr << "LadrunoRebarBuckling::LadrunoRebarBuckling -- failed to get a "
           << "copy of the wrapped material\n";
  }
  if (E <= 0.0 && theBar != 0)
    E = theBar->getInitialTangent();
  // Initialize OUR state inline -- do NOT call revertToStart() here: it would
  // forward to theBar->revertToStart() and wipe the state that bar.getCopy()
  // just preserved (the classic wrapper getCopy footgun; FatigueMaterial /
  // MinMaxMaterial avoid it the same way). getCopy() restores C*/T* afterward.
  CmaxTenStrain = CmaxTenStress = CfStarL = 0.0; CfStarLcross = LRB_INVALID;
  TmaxTenStrain = TmaxTenStress = TfStarL = 0.0; TfStarLcross = LRB_INVALID;
  Tstrain = Tstress = 0.0; Tr = 1.0;
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : E;
}

LadrunoRebarBuckling::LadrunoRebarBuckling()
  : UniaxialMaterial(0, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), model(MODEL_DM), lsr(0.0), alpha(1.0), fy(0.0), E(0.0),
    parameterID(0)
{
  // theBar is built by recvSelf via the broker.
  CmaxTenStrain = CmaxTenStress = CfStarL = 0.0; CfStarLcross = LRB_INVALID;
  TmaxTenStrain = TmaxTenStress = TfStarL = 0.0; TfStarLcross = LRB_INVALID;
  Tstrain = Tstress = Ttangent = 0.0; Tr = 1.0;
}

LadrunoRebarBuckling::~LadrunoRebarBuckling()
{
  if (theBar)
    delete theBar;
}

// ===========================================================================
//  committed-clone backbone probe -- the bare stress at an arbitrary strain,
//  reached monotonically from the last committed state. Material-agnostic.
// ===========================================================================
double LadrunoRebarBuckling::backboneStressAt(double engStrain)
{
  if (theBar == 0)
    return 0.0;
  UniaxialMaterial* probe = theBar->getCopy();   // copies theBar's trial state
  if (probe == 0)
    return 0.0;
  probe->revertToLastCommit();                   // back to the committed anchor
  probe->setTrialStrain(engStrain);              // push to the onset strain
  double s = probe->getStress();
  delete probe;
  return s;
}

// ===========================================================================
//  Dhakal-Maekawa buckled-average overlay (sets Tstress / Ttangent / Tr)
//  Ported from ReinforcingSteel::Buckled_stress_Dhakal(ess, fss); the clean
//  r*sigma_bare branch (the monotonic / TBranchNum%4>1 form).
// ===========================================================================
void LadrunoRebarBuckling::applyBuckling(double eps, double sBare, double kBare)
{
  const double eY = fy / E;                        // yield strain
  const double e_cross = TmaxTenStrain - TmaxTenStress / E;   // tensile-unload anchor
  const double e = eps - e_cross;                  // strain from the anchor

  // pre-onset (incl. tension): exact pass-through
  if (e >= -eY) {
    Tstress = sBare; Ttangent = kBare; Tr = 1.0;
    return;
  }

  const double kfac = sqrt(eY * 2000.0);           // = sqrt(fy/E*2000)
  double eStarMag = 55.0 - 2.3 * kfac * lsr;
  if (eStarMag < 7.0) eStarMag = 7.0;
  const double eStar = -eStarMag * eY;             // onset strain (negative)

  // backbone stress at the onset strain (cached per anchor). Mirrors upstream
  // ReinforcingSteel::Buckled_stress_Dhakal, which evaluates Backbone_f(eStar)
  // with eStar treated as the absolute backbone strain (it equals the onset's
  // absolute strain when e_cross==0, the monotonic-from-virgin case).
  if (TfStarLcross != e_cross) {
    TfStarL = backboneStressAt(eStar);
    TfStarLcross = e_cross;
  }
  const double fStarL = TfStarL;
  if (fabs(fStarL) < 1.0e-12) {                    // degenerate bare backbone => no overlay
    Tstress = sBare; Ttangent = kBare; Tr = 1.0;
    return;
  }

  double fStar = fStarL * alpha * (1.1 - 0.016 * kfac * lsr);
  if (fStar > -0.2 * fy) fStar = -0.2 * fy;        // residual floor

  if (e >= eStar) {
    // intermediate branch:  -eY > e >= eStar
    const double ratio = 1.0 - (1.0 - fStar / fStarL) * (e + eY) / (eStar + eY);
    const double dratio = -(1.0 - fStar / fStarL) / (eStar + eY);   // d(ratio)/d(eps)
    Tstress  = sBare * ratio;
    Ttangent = kBare * ratio + sBare * dratio;
    Tr = ratio;
  } else {
    // deep branch:  e < eStar  (ramp toward the -0.2 fy floor)
    const double num = fStar - 0.02 * E * (e - eStar);
    double ave = sBare * num / fStarL;
    if (ave > -0.2 * fy) {
      // floored: flat residual => (near) zero tangent
      Tstress  = -0.2 * fy;
      Ttangent = 0.0;
      Tr = (sBare != 0.0) ? Tstress / sBare : 1.0;
    } else {
      Tstress  = ave;
      Ttangent = (kBare * num + sBare * (-0.02 * E)) / fStarL;
      Tr = (sBare != 0.0) ? ave / sBare : 1.0;
    }
  }
}

// ===========================================================================
//  state setting
// ===========================================================================
int LadrunoRebarBuckling::setTrialStrain(double strain, double strainRate)
{
  if (theBar == 0)
    return -1;

  int res = theBar->setTrialStrain(strain, strainRate);
  double sBare = theBar->getStress();
  double kBare = theBar->getTangent();
  Tstrain = strain;

  // update the tensile-strain anchor (invalidate the fStarL cache if it moves)
  if (strain > TmaxTenStrain) {
    TmaxTenStrain = strain;
    TmaxTenStress = sBare;
    TfStarLcross  = LRB_INVALID;
  }

  if (lsr <= 0.0 || fy <= 0.0 || E <= 0.0) {
    // identity gate (B0): exact pass-through
    Tstress = sBare; Ttangent = kBare; Tr = 1.0;
    return res;
  }

  this->applyBuckling(strain, sBare, kBare);
  return res;
}

double LadrunoRebarBuckling::getStrainRate(void)
{
  return (theBar != 0) ? theBar->getStrainRate() : 0.0;
}

double LadrunoRebarBuckling::getDampTangent(void)
{
  return (theBar != 0) ? theBar->getDampTangent() : 0.0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRebarBuckling::commitState(void)
{
  CmaxTenStrain = TmaxTenStrain;
  CmaxTenStress = TmaxTenStress;
  CfStarL       = TfStarL;
  CfStarLcross  = TfStarLcross;
  return (theBar != 0) ? theBar->commitState() : 0;
}

int LadrunoRebarBuckling::revertToLastCommit(void)
{
  TmaxTenStrain = CmaxTenStrain;
  TmaxTenStress = CmaxTenStress;
  TfStarL       = CfStarL;
  TfStarLcross  = CfStarLcross;
  return (theBar != 0) ? theBar->revertToLastCommit() : 0;
}

int LadrunoRebarBuckling::revertToStart(void)
{
  CmaxTenStrain = CmaxTenStress = CfStarL = 0.0; CfStarLcross = LRB_INVALID;
  TmaxTenStrain = TmaxTenStress = TfStarL = 0.0; TfStarLcross = LRB_INVALID;
  Tstrain = Tstress = 0.0;
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : E;
  Tr = 1.0;
  return (theBar != 0) ? theBar->revertToStart() : 0;
}

UniaxialMaterial* LadrunoRebarBuckling::getCopy(void)
{
  LadrunoRebarBuckling* theCopy =
    new LadrunoRebarBuckling(this->getTag(), *theBar, lsr, alpha, fy, E, model);

  theCopy->CmaxTenStrain = CmaxTenStrain;
  theCopy->CmaxTenStress = CmaxTenStress;
  theCopy->CfStarL       = CfStarL;
  theCopy->CfStarLcross  = CfStarLcross;
  theCopy->TmaxTenStrain = TmaxTenStrain;
  theCopy->TmaxTenStress = TmaxTenStress;
  theCopy->TfStarL       = TfStarL;
  theCopy->TfStarLcross  = TfStarLcross;
  theCopy->Tstrain       = Tstrain;
  theCopy->Tstress       = Tstress;
  theCopy->Ttangent      = Ttangent;
  theCopy->Tr            = Tr;

  return theCopy;
}

// ===========================================================================
//  parallel
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

  static Vector data(10);
  data(0) = model;
  data(1) = lsr;
  data(2) = alpha;
  data(3) = fy;
  data(4) = E;
  data(5) = CmaxTenStrain;
  data(6) = CmaxTenStress;
  data(7) = CfStarL;
  data(8) = CfStarLcross;
  data(9) = this->getTag();
  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send Vector\n";
    return -2;
  }

  if (theBar->sendSelf(cTag, theChannel) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send the bar\n";
    return -3;
  }
  return 0;
}

int LadrunoRebarBuckling::recvSelf(int cTag, Channel& theChannel,
                                   FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  if (theBar == 0) {
    int matClassTag = int(dataID(1));
    theBar = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theBar == 0) {
      opserr << "LadrunoRebarBuckling::recvSelf - failed to create bar with "
             << "classTag " << matClassTag << endln;
      return -2;
    }
  }
  theBar->setDbTag(dataID(2));

  static Vector data(10);
  if (theChannel.recvVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv Vector\n";
    return -3;
  }
  model = int(data(0));
  lsr   = data(1);
  alpha = data(2);
  fy    = data(3);
  E     = data(4);
  CmaxTenStrain = data(5);
  CmaxTenStress = data(6);
  CfStarL       = data(7);
  CfStarLcross  = data(8);
  // trial = committed (recomputed on next setTrialStrain)
  TmaxTenStrain = CmaxTenStrain;
  TmaxTenStress = CmaxTenStress;
  TfStarL       = CfStarL;
  TfStarLcross  = CfStarLcross;
  Tstrain = 0.0; Tstress = 0.0; Ttangent = E; Tr = 1.0;

  if (theBar->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv the bar\n";
    return -4;
  }
  return 0;
}

// ===========================================================================
//  Print
// ===========================================================================
void LadrunoRebarBuckling::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoRebarBuckling\", ";
    s << "\"material\": \"" << theBar->getTag() << "\", ";
    s << "\"model\": \"" << (model == MODEL_DM ? "dm" : "ga") << "\", ";
    s << "\"lsr\": " << lsr << ", ";
    s << "\"alpha\": " << alpha << ", ";
    s << "\"fy\": " << fy << ", ";
    s << "\"E\": " << E << "}";
    return;
  }

  s << "LadrunoRebarBuckling, tag: " << this->getTag() << endln;
  s << "  wrapped material tag: " << theBar->getTag() << endln;
  s << "  model: " << (model == MODEL_DM ? "Dhakal-Maekawa" : "Gomes-Appleton") << endln;
  s << "  lsr (s/d): " << lsr << "  alpha: " << alpha << "  fy: " << fy
    << "  E: " << E << endln;
  s << "  wrapped material:" << endln;
  theBar->Print(s, flag);
}

// ===========================================================================
//  parameters
// ===========================================================================
int LadrunoRebarBuckling::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "lsr") == 0) {
    param.setValue(lsr);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "alpha") == 0) {
    param.setValue(alpha);
    return param.addObject(2, this);
  }
  // forward everything else to the wrapped bar (so element setParameter
  // can still reach the core material's parameters)
  return (theBar != 0) ? theBar->setParameter(argv, argc, param) : -1;
}

int LadrunoRebarBuckling::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case 1: lsr   = info.theDouble; return 0;
  case 2: alpha = info.theDouble; return 0;
  default: return -1;
  }
}

// ===========================================================================
//  responses
// ===========================================================================
Response* LadrunoRebarBuckling::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc == 0)
    return 0;

  // Responses we don't recognize are forwarded to the wrapped bar so core
  // quantities (plasticStrain, backStress, ...) stay reachable. Decide FIRST,
  // so forwarding does not emit a spurious empty UniaxialMaterialOutput block.
  const bool mine =
    (strcmp(argv[0], "stress") == 0) ||
    (strcmp(argv[0], "tangent") == 0) ||
    (strcmp(argv[0], "strain") == 0) ||
    (strcmp(argv[0], "stressStrain") == 0) ||
    (strcmp(argv[0], "stressANDstrain") == 0) ||
    (strcmp(argv[0], "reduction") == 0) ||
    (strcmp(argv[0], "buckling") == 0);

  if (!mine)
    return (theBar != 0) ? theBar->setResponse(argv, argc, s) : 0;

  Response* theResponse = 0;

  s.tag("UniaxialMaterialOutput");
  s.attr("matType", this->getClassType());
  s.attr("matTag", this->getTag());

  if (strcmp(argv[0], "stress") == 0) {
    s.tag("ResponseType", "sigma11");
    theResponse = new MaterialResponse(this, 1, this->getStress());
  }
  else if (strcmp(argv[0], "tangent") == 0) {
    s.tag("ResponseType", "C11");
    theResponse = new MaterialResponse(this, 2, this->getTangent());
  }
  else if (strcmp(argv[0], "strain") == 0) {
    s.tag("ResponseType", "eps11");
    theResponse = new MaterialResponse(this, 3, this->getStrain());
  }
  else if ((strcmp(argv[0], "stressStrain") == 0) ||
           (strcmp(argv[0], "stressANDstrain") == 0)) {
    s.tag("ResponseType", "sig11");
    s.tag("ResponseType", "eps11");
    theResponse = new MaterialResponse(this, 4, Vector(2));
  }
  else if ((strcmp(argv[0], "reduction") == 0) ||
           (strcmp(argv[0], "buckling") == 0)) {
    s.tag("ResponseType", "r");
    theResponse = new MaterialResponse(this, 5, Tr);
  }

  s.endTag();
  return theResponse;
}

int LadrunoRebarBuckling::getResponse(int responseID, Information& matInfo)
{
  static Vector stressStrain(2);
  switch (responseID) {
  case 1: matInfo.setDouble(this->getStress());  return 0;
  case 2: matInfo.setDouble(this->getTangent()); return 0;
  case 3: matInfo.setDouble(this->getStrain());  return 0;
  case 4:
    stressStrain(0) = this->getStress();
    stressStrain(1) = this->getStrain();
    matInfo.setVector(stressStrain);
    return 0;
  case 5: matInfo.setDouble(Tr); return 0;
  default: return -1;
  }
}
