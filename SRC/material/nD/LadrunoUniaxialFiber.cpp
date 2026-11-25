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

// Ladruno WP-245 (program #1131, rcbf-3): see LadrunoUniaxialFiber.h.

#include "LadrunoUniaxialFiber.h"

#include <UniaxialMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Response.h>
#include <MaterialResponse.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <LadrunoOptSpec.h>

#include <string.h>

void* OPS_LadrunoUniaxialFiber()
{
  static const char* OWNER = "LadrunoUniaxialFiber";
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber: want  nDMaterial LadrunoUniaxialFiber "
              "$tag $uniTag [-rho $rho] [-dowel $fc $fy $db $ls [-dowelK $kd] [-dowelSlip $su] [-dowelDecay $Cd]]\n";
    return 0;
  }
  int idata[2], two = 2;
  if (OPS_GetIntInput(&two, idata) < 0) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber: invalid $tag $uniTag\n";
    return 0;
  }
  const int tag = idata[0];
  double rho = -1.0;
  bool dowel = false;
  ladruno_dowel::Params dp;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    const ladruno_opt::OptSpec* o = ladruno_opt::find(ladruno_opt::kLadrunoUniaxialFiberOptions, opt);
    if (o == 0) { ladruno_opt::reportUnknown("LadrunoUniaxialFiber", tag, opt, ladruno_opt::kLadrunoUniaxialFiberOptions); return 0; }
    if (!ladruno_opt::haveValues(OWNER, tag, *o)) return 0;
    if (strcmp(opt, "-rho") == 0) {
      int one = 1;
      if (OPS_GetDoubleInput(&one, &rho) < 0 || rho < 0.0) {
        opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": -rho wants a value >= 0\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-dowel") == 0) {               // WP-257
      double v[4]; int four = 4;
      if (OPS_GetDoubleInput(&four, v) < 0) {
        opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": -dowel wants $fc $fy $db $ls\n";
        return 0;
      }
      dp.fc = v[0]; dp.fy = v[1]; dp.db = v[2]; dp.ls = v[3];
      dowel = true;
    }
    else if (strcmp(opt, "-dowelK") == 0 || strcmp(opt, "-dowelSlip") == 0 || strcmp(opt, "-dowelDecay") == 0) {
      double v; int one = 1;
      const bool decay = (opt[6] == 'D');
      if (OPS_GetDoubleInput(&one, &v) < 0 || !(decay ? v >= 0.0 : v > 0.0)) {
        opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": " << opt << " wants a value "
               << (decay ? ">= 0" : "> 0") << "\n";
        return 0;
      }
      if (opt[6] == 'K') dp.kd = v; else if (decay) dp.cd = v; else dp.su = v;
    }
    else { ladruno_opt::reportUnhandled(OWNER, tag, opt); return 0; }
  }
  if (!dowel && (dp.kd != 1.3 || dp.su > 0.0 || dp.cd > 0.0)) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": -dowelK / -dowelSlip / -dowelDecay need -dowel\n";
    return 0;
  }
  if (dowel && !dp.valid()) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": -dowel $fc $fy $db $ls must all be > 0\n";
    return 0;
  }
  UniaxialMaterial* uni = OPS_getUniaxialMaterial(idata[1]);
  if (uni == 0) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": uniaxial material "
           << idata[1] << " not found\n";
    return 0;
  }
  LadrunoUniaxialFiber* m = new LadrunoUniaxialFiber(tag, *uni, 3, rho);
  if (dowel) m->setDowel(dp);
  return m;
}

LadrunoUniaxialFiber::LadrunoUniaxialFiber(int tag, UniaxialMaterial& u, int ord, double r)
  : NDMaterial(tag, ND_TAG_LadrunoUniaxialFiber), theUni(0), order(ord), rho(r),
    eps3(3), sig3(3), eps2(2), sig2(2), D3(3, 3), D2(2, 2), dowel(false), dp()
{
  smaxC[0] = smaxC[1] = scC[0] = scC[1] = DC[0] = DC[1] = 0.0;
  theUni = u.getCopy();
  if (theUni == 0) {
    opserr << "LadrunoUniaxialFiber " << tag << ": failed to copy the uniaxial material\n";
  }
}

LadrunoUniaxialFiber::LadrunoUniaxialFiber()
  : NDMaterial(0, ND_TAG_LadrunoUniaxialFiber), theUni(0), order(3), rho(-1.0),
    eps3(3), sig3(3), eps2(2), sig2(2), D3(3, 3), D2(2, 2), dowel(false), dp()
{
  smaxC[0] = smaxC[1] = scC[0] = scC[1] = DC[0] = DC[1] = 0.0;
}

void
LadrunoUniaxialFiber::setDowel(const ladruno_dowel::Params& P)
{
  dowel = true;
  dp = P;
}

// WP-257: the dowel response of shear component k (0 = gxy, 1 = gxz) at the trial strain
ladruno_dowel::Out
LadrunoUniaxialFiber::dowelOut(int k)
{
  const Vector& e = (order == 3) ? eps3 : eps2;
  return ladruno_dowel::shear(dp, e(1 + k), theUni->getStress(), smaxC[k], scC[k], DC[k]);
}

LadrunoUniaxialFiber::~LadrunoUniaxialFiber()
{
  if (theUni != 0) delete theUni;
}

int
LadrunoUniaxialFiber::setTrialStrain(const Vector& v)
{
  if (v.Size() != order) {
    opserr << "LadrunoUniaxialFiber " << this->getTag() << ": strain of size " << v.Size()
           << " for an order-" << order << " view\n";
    return -1;
  }
  if (order == 3) eps3 = v; else eps2 = v;
  return theUni->setTrialStrain(v(0));
}

int LadrunoUniaxialFiber::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int
LadrunoUniaxialFiber::setTrialStrainIncr(const Vector& dv)
{
  static Vector v3(3), v2(2);
  if (order == 3) { v3 = eps3; v3.addVector(1.0, dv, 1.0); return this->setTrialStrain(v3); }
  v2 = eps2; v2.addVector(1.0, dv, 1.0); return this->setTrialStrain(v2);
}

int LadrunoUniaxialFiber::setTrialStrainIncr(const Vector& dv, const Vector&) { return this->setTrialStrainIncr(dv); }

const Matrix&
LadrunoUniaxialFiber::getTangent(void)
{
  Matrix& D = (order == 3) ? D3 : D2;
  D.Zero();
  D(0, 0) = theUni->getTangent();
  if (dowel) {                                   // WP-257: tau_k(gamma_k, sigma(eps))
    for (int k = 0; k < order - 1; k++) {
      const ladruno_dowel::Out o = dowelOut(k);
      D(1 + k, 1 + k) = o.dg;
      D(1 + k, 0) = o.dsig * D(0, 0);
    }
  }
  return D;
}

const Matrix&
LadrunoUniaxialFiber::getInitialTangent(void)
{
  Matrix& D = (order == 3) ? D3 : D2;
  D.Zero();
  D(0, 0) = theUni->getInitialTangent();
  if (dowel) {                                   // WP-257: the virgin slope at zero stress
    const double g0 = dp.tauU() * dp.ls / dp.s0();
    for (int k = 0; k < order - 1; k++) D(1 + k, 1 + k) = g0;
  }
  return D;
}

const Vector&
LadrunoUniaxialFiber::getStress(void)
{
  Vector& s = (order == 3) ? sig3 : sig2;
  s.Zero();
  s(0) = theUni->getStress();
  if (dowel)                                     // WP-257
    for (int k = 0; k < order - 1; k++) s(1 + k) = dowelOut(k).tau;
  return s;
}

const Vector&
LadrunoUniaxialFiber::getStrain(void)
{
  return (order == 3) ? eps3 : eps2;
}

double
LadrunoUniaxialFiber::getRho(void)
{
  return (rho >= 0.0) ? rho : theUni->getRho();
}

int
LadrunoUniaxialFiber::commitState(void)
{
  if (dowel) {                                   // WP-257: the largest slip so far
    const Vector& e = (order == 3) ? eps3 : eps2;
    for (int k = 0; k < order - 1; k++) {
      if (dp.cd > 0.0) {                         // the re-slip the converged trial used, then the new slip
        DC[k] = ladruno_dowel::commitReslip(dp, e(1 + k), smaxC[k], scC[k], DC[k]);
        scC[k] = std::fabs(e(1 + k)) * dp.ls;
      }
      smaxC[k] = ladruno_dowel::commitSlip(dp, e(1 + k), smaxC[k]);
    }
  }
  return theUni->commitState();
}

int LadrunoUniaxialFiber::revertToLastCommit(void) { return theUni->revertToLastCommit(); }

int
LadrunoUniaxialFiber::revertToStart(void)
{
  eps3.Zero(); eps2.Zero();
  smaxC[0] = smaxC[1] = scC[0] = scC[1] = DC[0] = DC[1] = 0.0;
  return theUni->revertToStart();
}

NDMaterial*
LadrunoUniaxialFiber::getCopy(void)
{
  LadrunoUniaxialFiber* c = new LadrunoUniaxialFiber(this->getTag(), *theUni, order, rho);
  c->eps3 = eps3; c->eps2 = eps2;
  c->dowel = dowel; c->dp = dp;
  for (int k = 0; k < 2; k++) { c->smaxC[k] = smaxC[k]; c->scC[k] = scC[k]; c->DC[k] = DC[k]; }
  return c;
}

NDMaterial*
LadrunoUniaxialFiber::getCopy(const char* type)
{
  LadrunoUniaxialFiber* c = 0;
  if (strcmp(type, "BeamFiber") == 0 || strcmp(type, "TimoshenkoFiber") == 0)
    c = new LadrunoUniaxialFiber(this->getTag(), *theUni, 3, rho);
  else if (strcmp(type, "BeamFiber2d") == 0 || strcmp(type, "TimoshenkoFiber2d") == 0)
    c = new LadrunoUniaxialFiber(this->getTag(), *theUni, 2, rho);
  if (c != 0) {
    if (dowel) c->setDowel(dp);
    return c;
  }
  opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << this->getTag() << ": no '" << type
         << "' view -- a bar is a beam fibre (BeamFiber, BeamFiber2d) only (refused)\n";
  return 0;
}

const char*
LadrunoUniaxialFiber::getType(void) const
{
  return (order == 3) ? "BeamFiber" : "BeamFiber2d";
}

int LadrunoUniaxialFiber::getOrder(void) const { return order; }

Response*
LadrunoUniaxialFiber::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc > 0 && strcmp(argv[0], "dowel") == 0) {   // WP-257: (tau_y, tau_z, smax_y, smax_z, cap, D_y, D_z)
    static Vector v(7);
    return new MaterialResponse(this, 101, v);
  }
  Response* r = NDMaterial::setResponse(argv, argc, s);
  if (r != 0) return r;
  return theUni->setResponse(argv, argc, s);
}

int
LadrunoUniaxialFiber::getResponse(int responseID, Information& info)
{
  if (responseID == 101) {
    static Vector v(7);
    v.Zero();
    if (dowel) {
      for (int k = 0; k < order - 1; k++) { v(k) = dowelOut(k).tau; v(2 + k) = smaxC[k]; v(5 + k) = DC[k]; }
      v(4) = ladruno_dowel::cap(dp, theUni->getStress(), 0);
    }
    return info.setVector(v);
  }
  return NDMaterial::getResponse(responseID, info);
}

int
LadrunoUniaxialFiber::sendSelf(int commitTag, Channel& theChannel)
{
  if (theUni == 0) return -1;
  static ID idata(5);
  idata(0) = this->getTag();
  idata(4) = dowel ? 1 : 0;                      // WP-257
  idata(1) = order;
  idata(2) = theUni->getClassTag();
  int uniDb = theUni->getDbTag();
  if (uniDb == 0) {
    uniDb = theChannel.getDbTag();
    if (uniDb != 0) theUni->setDbTag(uniDb);
  }
  idata(3) = uniDb;
  int res = theChannel.sendID(this->getDbTag(), commitTag, idata);
  if (res < 0) { opserr << "LadrunoUniaxialFiber::sendSelf - failed to send ID\n"; return res; }
  static Vector ddata(5);
  ddata.Zero();
  ddata(0) = rho;
  for (int i = 0; i < order; i++) ddata(1 + i) = (order == 3) ? eps3(i) : eps2(i);
  res = theChannel.sendVector(this->getDbTag(), commitTag, ddata);
  if (res < 0) { opserr << "LadrunoUniaxialFiber::sendSelf - failed to send data\n"; return res; }
  if (dowel) {                                   // WP-257: a size-13 vector, never the 5 above
    static Vector w(13);
    w(0) = dp.fc; w(1) = dp.fy; w(2) = dp.db; w(3) = dp.ls; w(4) = dp.kd; w(5) = dp.su;
    w(6) = smaxC[0]; w(7) = smaxC[1]; w(8) = dp.cd; w(9) = scC[0]; w(10) = scC[1]; w(11) = DC[0]; w(12) = DC[1];
    res = theChannel.sendVector(this->getDbTag(), commitTag, w);
    if (res < 0) { opserr << "LadrunoUniaxialFiber::sendSelf - failed to send the dowel data\n"; return res; }
  }
  return theUni->sendSelf(commitTag, theChannel);
}

int
LadrunoUniaxialFiber::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static ID idata(5);
  int res = theChannel.recvID(this->getDbTag(), commitTag, idata);
  if (res < 0) { opserr << "LadrunoUniaxialFiber::recvSelf - failed to recv ID\n"; return res; }
  this->setTag(idata(0));
  order = idata(1);
  if (theUni == 0 || theUni->getClassTag() != idata(2)) {
    if (theUni != 0) delete theUni;
    theUni = theBroker.getNewUniaxialMaterial(idata(2));
    if (theUni == 0) {
      opserr << "LadrunoUniaxialFiber::recvSelf - no uniaxial material of classTag " << idata(2) << "\n";
      return -1;
    }
  }
  theUni->setDbTag(idata(3));
  static Vector ddata(5);
  res = theChannel.recvVector(this->getDbTag(), commitTag, ddata);
  if (res < 0) { opserr << "LadrunoUniaxialFiber::recvSelf - failed to recv data\n"; return res; }
  rho = ddata(0);
  for (int i = 0; i < 3; i++) eps3(i) = ddata(1 + i);
  for (int i = 0; i < 2; i++) eps2(i) = ddata(1 + i);
  dowel = (idata(4) == 1);
  smaxC[0] = smaxC[1] = scC[0] = scC[1] = DC[0] = DC[1] = 0.0;
  if (dowel) {                                   // WP-257
    static Vector w(13);
    res = theChannel.recvVector(this->getDbTag(), commitTag, w);
    if (res < 0) { opserr << "LadrunoUniaxialFiber::recvSelf - failed to recv the dowel data\n"; return res; }
    dp.fc = w(0); dp.fy = w(1); dp.db = w(2); dp.ls = w(3); dp.kd = w(4); dp.su = w(5);
    smaxC[0] = w(6); smaxC[1] = w(7); dp.cd = w(8);
    scC[0] = w(9); scC[1] = w(10); DC[0] = w(11); DC[1] = w(12);
  }
  return theUni->recvSelf(commitTag, theChannel, theBroker);
}

void
LadrunoUniaxialFiber::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoUniaxialFiber, tag: " << this->getTag() << " (" << this->getType()
    << " bar fibre, " << (dowel ? "dowel shear" : "no shear") << "), uniaxial material: ";
  if (theUni != 0) s << theUni->getTag() << " (" << theUni->getClassType() << ")";
  s << endln;
  if (dowel)
    s << "  dowel: fc " << dp.fc << " fy " << dp.fy << " db " << dp.db << " ls " << dp.ls
      << " kd " << dp.kd << " su " << dp.slipU() << " decay " << dp.cd << " (tau_u " << dp.tauU() << ")" << endln;
  if (theUni != 0) theUni->Print(s, flag);
}
