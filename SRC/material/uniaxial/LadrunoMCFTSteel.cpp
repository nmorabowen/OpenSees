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

// Ladruno WP-254 (program #1131): see LadrunoMCFTSteel.h.

#include "LadrunoMCFTSteel.h"

#include <Vector.h>
#include <Channel.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <LadrunoOptSpec.h>

#include <string.h>
#include <cmath>

using ladruno_mcft::Steel;

void* OPS_LadrunoMCFTSteel()
{
  static const char* OWNER = "LadrunoMCFTSteel";
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "ERROR uniaxialMaterial LadrunoMCFTSteel: want  uniaxialMaterial LadrunoMCFTSteel $tag "
              "$fy $fu [-Es $Es] [-esh $esh] [-esu $esu]\n";
    return 0;
  }
  int tag, one = 1;
  if (OPS_GetIntInput(&one, &tag) < 0) {
    opserr << "ERROR uniaxialMaterial LadrunoMCFTSteel: invalid $tag\n";
    return 0;
  }
  double d2[2]; int two = 2;
  if (OPS_GetDoubleInput(&two, d2) < 0) {
    opserr << "ERROR uniaxialMaterial LadrunoMCFTSteel " << tag << ": want $fy $fu\n";
    return 0;
  }
  Steel S;
  S.fy = d2[0]; S.fu = d2[1];
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    const ladruno_opt::OptSpec* o = ladruno_opt::find(ladruno_opt::kLadrunoMCFTSteelOptions, opt);
    if (o == 0) { ladruno_opt::reportUnknown("LadrunoMCFTSteel", tag, opt, ladruno_opt::kLadrunoMCFTSteelOptions); return 0; }
    if (!ladruno_opt::haveValues(OWNER, tag, *o)) return 0;
    double v;
    if (OPS_GetDoubleInput(&one, &v) < 0) {
      opserr << "ERROR uniaxialMaterial LadrunoMCFTSteel " << tag << ": " << opt << " wants a number\n";
      return 0;
    }
    if      (strcmp(opt, "-Es") == 0)  S.Es = v;
    else if (strcmp(opt, "-esh") == 0) S.esh = v;
    else if (strcmp(opt, "-esu") == 0) S.esu = v;
    else { ladruno_opt::reportUnhandled(OWNER, tag, opt); return 0; }
  }
  if (!(S.fy > 0.0) || !(S.fu >= S.fy) || !(S.Es > 0.0) || !(S.esh > S.fy / S.Es) || !(S.esu > S.esh)) {
    opserr << "ERROR uniaxialMaterial LadrunoMCFTSteel " << tag << ": want fy > 0, fu >= fy, Es > 0 and "
              "fy/Es < esh < esu\n";
    return 0;
  }
  return new LadrunoMCFTSteel(tag, S);
}

LadrunoMCFTSteel::LadrunoMCFTSteel(int tag, const Steel& s)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoMCFTSteel), steel(s),
    epC(0.0), epsC(0.0), epsT(0.0), sigT(0.0), EtT(s.Es)
{
}

LadrunoMCFTSteel::LadrunoMCFTSteel()
  : UniaxialMaterial(0, MAT_TAG_LadrunoMCFTSteel), steel(),
    epC(0.0), epsC(0.0), epsT(0.0), sigT(0.0), EtT(0.0)
{
}

LadrunoMCFTSteel::~LadrunoMCFTSteel()
{
}

int
LadrunoMCFTSteel::setTrialStrain(double strain, double)
{
  if (!std::isfinite(strain)) return -1;
  epsT = strain;
  sigT = ladruno_mcft::steelStress(steel, epsT, epC, &EtT);
  return 0;
}

int
LadrunoMCFTSteel::commitState(void)
{
  epC = epsT - sigT / steel.Es;
  epsC = epsT;
  return 0;
}

int
LadrunoMCFTSteel::revertToLastCommit(void)
{
  return this->setTrialStrain(epsC);
}

int
LadrunoMCFTSteel::revertToStart(void)
{
  epC = epsC = 0.0;
  return this->setTrialStrain(0.0);
}

UniaxialMaterial*
LadrunoMCFTSteel::getCopy(void)
{
  LadrunoMCFTSteel* c = new LadrunoMCFTSteel(this->getTag(), steel);
  c->epC = epC; c->epsC = epsC;
  c->epsT = epsT; c->sigT = sigT; c->EtT = EtT;
  return c;
}

int
LadrunoMCFTSteel::sendSelf(int commitTag, Channel& theChannel)
{
  Vector d(9);
  d(0) = this->getTag();
  d(1) = steel.fy; d(2) = steel.fu; d(3) = steel.Es; d(4) = steel.esh; d(5) = steel.esu;
  d(6) = epC; d(7) = epsC; d(8) = epsT;
  if (theChannel.sendVector(this->getDbTag(), commitTag, d) < 0) {
    opserr << "LadrunoMCFTSteel::sendSelf - failed to send data\n";
    return -1;
  }
  return 0;
}

int
LadrunoMCFTSteel::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  Vector d(9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, d) < 0) {
    opserr << "LadrunoMCFTSteel::recvSelf - failed to recv data\n";
    return -1;
  }
  this->setTag((int)d(0));
  steel.fy = d(1); steel.fu = d(2); steel.Es = d(3); steel.esh = d(4); steel.esu = d(5);
  epC = d(6); epsC = d(7);
  return this->setTrialStrain(d(8));
}

void
LadrunoMCFTSteel::Print(OPS_Stream& s, int)
{
  s << "LadrunoMCFTSteel, tag: " << this->getTag() << "  fy " << steel.fy << "  fu " << steel.fu
    << "  Es " << steel.Es << "  esh " << steel.esh << "  esu " << steel.esu
    << "  (committed plastic strain " << epC << ")" << endln;
}
