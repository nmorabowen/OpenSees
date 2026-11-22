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
              "$tag $uniTag [-rho $rho]\n";
    return 0;
  }
  int idata[2], two = 2;
  if (OPS_GetIntInput(&two, idata) < 0) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber: invalid $tag $uniTag\n";
    return 0;
  }
  const int tag = idata[0];
  double rho = -1.0;
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
    else { ladruno_opt::reportUnhandled(OWNER, tag, opt); return 0; }
  }
  UniaxialMaterial* uni = OPS_getUniaxialMaterial(idata[1]);
  if (uni == 0) {
    opserr << "ERROR nDMaterial LadrunoUniaxialFiber " << tag << ": uniaxial material "
           << idata[1] << " not found\n";
    return 0;
  }
  return new LadrunoUniaxialFiber(tag, *uni, 3, rho);
}

LadrunoUniaxialFiber::LadrunoUniaxialFiber(int tag, UniaxialMaterial& u, int ord, double r)
  : NDMaterial(tag, ND_TAG_LadrunoUniaxialFiber), theUni(0), order(ord), rho(r),
    eps3(3), sig3(3), eps2(2), sig2(2), D3(3, 3), D2(2, 2)
{
  theUni = u.getCopy();
  if (theUni == 0) {
    opserr << "LadrunoUniaxialFiber " << tag << ": failed to copy the uniaxial material\n";
  }
}

LadrunoUniaxialFiber::LadrunoUniaxialFiber()
  : NDMaterial(0, ND_TAG_LadrunoUniaxialFiber), theUni(0), order(3), rho(-1.0),
    eps3(3), sig3(3), eps2(2), sig2(2), D3(3, 3), D2(2, 2)
{
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
  return D;
}

const Matrix&
LadrunoUniaxialFiber::getInitialTangent(void)
{
  Matrix& D = (order == 3) ? D3 : D2;
  D.Zero();
  D(0, 0) = theUni->getInitialTangent();
  return D;
}

const Vector&
LadrunoUniaxialFiber::getStress(void)
{
  Vector& s = (order == 3) ? sig3 : sig2;
  s.Zero();
  s(0) = theUni->getStress();
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

int LadrunoUniaxialFiber::commitState(void)        { return theUni->commitState(); }
int LadrunoUniaxialFiber::revertToLastCommit(void) { return theUni->revertToLastCommit(); }

int
LadrunoUniaxialFiber::revertToStart(void)
{
  eps3.Zero(); eps2.Zero();
  return theUni->revertToStart();
}

NDMaterial*
LadrunoUniaxialFiber::getCopy(void)
{
  LadrunoUniaxialFiber* c = new LadrunoUniaxialFiber(this->getTag(), *theUni, order, rho);
  c->eps3 = eps3; c->eps2 = eps2;
  return c;
}

NDMaterial*
LadrunoUniaxialFiber::getCopy(const char* type)
{
  if (strcmp(type, "BeamFiber") == 0 || strcmp(type, "TimoshenkoFiber") == 0)
    return new LadrunoUniaxialFiber(this->getTag(), *theUni, 3, rho);
  if (strcmp(type, "BeamFiber2d") == 0 || strcmp(type, "TimoshenkoFiber2d") == 0)
    return new LadrunoUniaxialFiber(this->getTag(), *theUni, 2, rho);
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
  Response* r = NDMaterial::setResponse(argv, argc, s);
  if (r != 0) return r;
  return theUni->setResponse(argv, argc, s);
}

int
LadrunoUniaxialFiber::getResponse(int responseID, Information& info)
{
  return NDMaterial::getResponse(responseID, info);
}

int
LadrunoUniaxialFiber::sendSelf(int commitTag, Channel& theChannel)
{
  if (theUni == 0) return -1;
  static ID idata(4);
  idata(0) = this->getTag();
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
  return theUni->sendSelf(commitTag, theChannel);
}

int
LadrunoUniaxialFiber::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static ID idata(4);
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
  return theUni->recvSelf(commitTag, theChannel, theBroker);
}

void
LadrunoUniaxialFiber::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoUniaxialFiber, tag: " << this->getTag() << " (" << this->getType()
    << " bar fibre, no shear), uniaxial material: ";
  if (theUni != 0) s << theUni->getTag() << " (" << theUni->getClassType() << ")";
  s << endln;
  if (theUni != 0) theUni->Print(s, flag);
}
