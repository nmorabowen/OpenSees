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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// LadrunoIMKBeam2d implementation. See LadrunoIMKBeam2d.h and
// Ladruno_implementation/14_ladruno_imk_beam.md for the formulation. The
// per-bending-axis hinge kernel is shared with the 3D element through
// LadrunoIMKHinge.h.

#include "LadrunoIMKBeam2d.h"
#include "LadrunoIMKHinge.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <CrdTransf.h>
#include <UniaxialMaterial.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <classTags.h>
#include <OPS_Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

Matrix LadrunoIMKBeam2d::K(6, 6);
Vector LadrunoIMKBeam2d::P(6);

// ---------------------------------------------------------------------------
// constructors / destructor
// ---------------------------------------------------------------------------
LadrunoIMKBeam2d::LadrunoIMKBeam2d()
  : Element(0, ELE_TAG_LadrunoIMKBeam2d),
    A(0.0), E(0.0), Iz(0.0), rho(0.0),
    q(3), kb(3, 3), Q(6),
    connectedExternalNodes(2), theCoordTransf(0)
{
  for (int i = 0; i < 2; i++) {
    theMat[i] = 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }
  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam2d::LadrunoIMKBeam2d(int tag, int Nd1, int Nd2,
                                   double a, double e, double iz,
                                   CrdTransf &theTransf,
                                   UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                                   double r)
  : Element(tag, ELE_TAG_LadrunoIMKBeam2d),
    A(a), E(e), Iz(iz), rho(r),
    q(3), kb(3, 3), Q(6),
    connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  UniaxialMaterial *src[2] = { mZi, mZj };
  for (int i = 0; i < 2; i++) {
    theMat[i] = (src[i] != 0) ? src[i]->getCopy() : 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }

  theCoordTransf = theTransf.getCopy2d();
  if (theCoordTransf == 0) {
    opserr << "LadrunoIMKBeam2d::LadrunoIMKBeam2d -- failed to get copy of CrdTransf\n";
    exit(-1);
  }

  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam2d::~LadrunoIMKBeam2d()
{
  for (int i = 0; i < 2; i++)
    if (theMat[i] != 0)
      delete theMat[i];
  if (theCoordTransf != 0)
    delete theCoordTransf;
}

// ---------------------------------------------------------------------------
// trivial accessors
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::getNumExternalNodes(void) const { return 2; }
const ID &LadrunoIMKBeam2d::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoIMKBeam2d::getNodePtrs(void) { return theNodes; }
int LadrunoIMKBeam2d::getNumDOF(void) { return 6; }

void LadrunoIMKBeam2d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }

  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

  if (theNodes[0] == 0 || theNodes[1] == 0) {
    opserr << "LadrunoIMKBeam2d::setDomain -- node not found for element "
           << this->getTag() << endln;
    return;
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  if (dofNd1 != 3 || dofNd2 != 3) {
    opserr << "LadrunoIMKBeam2d::setDomain -- both nodes must have 3 DOF (element "
           << this->getTag() << ")\n";
    return;
  }

  this->DomainComponent::setDomain(theDomain);

  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "LadrunoIMKBeam2d::setDomain -- CrdTransf failed to initialize\n";
    return;
  }

  double L = theCoordTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "LadrunoIMKBeam2d::setDomain -- element " << this->getTag()
           << " has zero length\n";
    return;
  }
}

// ---------------------------------------------------------------------------
// state cycle
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::commitState(void)
{
  int ok = 0;
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->commitState();
    thetaHcommit[i] = thetaH[i];
  }
  ok += theCoordTransf->commitState();
  return ok;
}

int LadrunoIMKBeam2d::revertToLastCommit(void)
{
  int ok = 0;
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToLastCommit();
    thetaH[i] = thetaHcommit[i];
  }
  ok += theCoordTransf->revertToLastCommit();
  return ok;
}

int LadrunoIMKBeam2d::revertToStart(void)
{
  int ok = 0;
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToStart();
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }
  ok += theCoordTransf->revertToStart();
  return ok;
}

// ---------------------------------------------------------------------------
// update: full basic-system state determination, caches q and kb
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::update(void)
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;
  const double EAoverL = E * A * oneOverL;

  kb.Zero();
  q.Zero();

  // axial: linear elastic, no hinge
  kb(0, 0) = EAoverL;  q(0) = EAoverL * v(0);

  double k2[2][2];

  // strong-axis bending: basic dofs 1,2 ; materials 0,1 ; EI = E*Iz
  double qi, qj;
  ladrunoIMKSolveAxis(theMat[0], theMat[1], v(1), v(2), L, E * Iz,
                      thetaH[0], thetaH[1], qi, qj, k2);
  q(1) = qi;  q(2) = qj;
  kb(1, 1) = k2[0][0];  kb(1, 2) = k2[0][1];
  kb(2, 1) = k2[1][0];  kb(2, 2) = k2[1][1];

  return ok;
}

const Matrix &LadrunoIMKBeam2d::getTangentStiff(void)
{
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &LadrunoIMKBeam2d::getInitialStiff(void)
{
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;

  static Matrix kbi(3, 3);
  kbi.Zero();
  kbi(0, 0) = E * A * oneOverL;

  double k2[2][2];
  ladrunoIMKInitBlock(theMat[0], theMat[1], L, E * Iz, k2);
  kbi(1, 1) = k2[0][0];  kbi(1, 2) = k2[0][1];
  kbi(2, 1) = k2[1][0];  kbi(2, 2) = k2[1][1];

  return theCoordTransf->getInitialGlobalStiffMatrix(kbi);
}

const Matrix &LadrunoIMKBeam2d::getMass(void)
{
  K.Zero();
  if (rho > 0.0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;
    K(0, 0) = K(1, 1) = m;   // node i translations
    K(3, 3) = K(4, 4) = m;   // node j translations
  }
  return K;
}

// ---------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------
void LadrunoIMKBeam2d::zeroLoad(void) { Q.Zero(); }

int LadrunoIMKBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "LadrunoIMKBeam2d::addLoad -- element loads not supported in v1 "
            "(element " << this->getTag() << "); use nodal loads\n";
  return -1;
}

int LadrunoIMKBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  double L = theCoordTransf->getInitialLength();
  double m = 0.5 * rho * L;

  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  Q(3) -= m * Raccel2(0);
  Q(4) -= m * Raccel2(1);

  return 0;
}

// ---------------------------------------------------------------------------
// resisting force
// ---------------------------------------------------------------------------
const Vector &LadrunoIMKBeam2d::getResistingForce(void)
{
  static Vector p0Vec(3);
  p0Vec.Zero();

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // subtract element nodal loads (inertia unbalance accumulated in Q)
  P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector &LadrunoIMKBeam2d::getResistingForceIncInertia(void)
{
  P = this->getResistingForce();

  // Rayleigh damping forces
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  if (rho == 0.0)
    return P;

  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();

  double L = theCoordTransf->getInitialLength();
  double m = 0.5 * rho * L;

  P(0) += m * accel1(0);
  P(1) += m * accel1(1);
  P(3) += m * accel2(0);
  P(4) += m * accel2(1);

  return P;
}

// ---------------------------------------------------------------------------
// parallel
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();

  // integer data: tag, nodes, transf class+db, 2x (mat class+db)
  static ID idData(9);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = theCoordTransf->getClassTag();
  int crdDb = theCoordTransf->getDbTag();
  if (crdDb == 0) {
    crdDb = theChannel.getDbTag();
    if (crdDb != 0)
      theCoordTransf->setDbTag(crdDb);
  }
  idData(4) = crdDb;

  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0) {
      idData(5 + 2 * i) = theMat[i]->getClassTag();
      int mDb = theMat[i]->getDbTag();
      if (mDb == 0) {
        mDb = theChannel.getDbTag();
        if (mDb != 0)
          theMat[i]->setDbTag(mDb);
      }
      idData(6 + 2 * i) = mDb;
    } else {
      idData(5 + 2 * i) = 0;
      idData(6 + 2 * i) = 0;
    }
  }

  res += theChannel.sendID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::sendSelf -- failed to send ID\n";
    return res;
  }

  static Vector dData(6);
  dData(0) = A;  dData(1) = E;  dData(2) = Iz;  dData(3) = rho;
  dData(4) = thetaHcommit[0];  dData(5) = thetaHcommit[1];

  res += theChannel.sendVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::sendSelf -- failed to send Vector\n";
    return res;
  }

  res += theCoordTransf->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::sendSelf -- failed to send CrdTransf\n";
    return res;
  }

  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0) {
      res += theMat[i]->sendSelf(cTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoIMKBeam2d::sendSelf -- failed to send material " << i << endln;
        return res;
      }
    }
  }

  return res;
}

int LadrunoIMKBeam2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static ID idData(9);
  res += theChannel.recvID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::recvSelf -- failed to recv ID\n";
    return res;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);

  static Vector dData(6);
  res += theChannel.recvVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::recvSelf -- failed to recv Vector\n";
    return res;
  }
  A = dData(0);  E = dData(1);  Iz = dData(2); rho = dData(3);
  for (int i = 0; i < 2; i++) {
    thetaHcommit[i] = dData(4 + i);
    thetaH[i] = thetaHcommit[i];
  }

  // CrdTransf
  int crdClass = idData(3);
  int crdDb = idData(4);
  if (theCoordTransf == 0 || theCoordTransf->getClassTag() != crdClass) {
    if (theCoordTransf != 0)
      delete theCoordTransf;
    theCoordTransf = theBroker.getNewCrdTransf(crdClass);
    if (theCoordTransf == 0) {
      opserr << "LadrunoIMKBeam2d::recvSelf -- could not create CrdTransf\n";
      return -1;
    }
  }
  theCoordTransf->setDbTag(crdDb);
  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::recvSelf -- failed to recv CrdTransf\n";
    return res;
  }

  // materials
  for (int i = 0; i < 2; i++) {
    int mClass = idData(5 + 2 * i);
    int mDb = idData(6 + 2 * i);
    if (mClass == 0) {
      if (theMat[i] != 0) {
        delete theMat[i];
        theMat[i] = 0;
      }
      continue;
    }
    if (theMat[i] == 0 || theMat[i]->getClassTag() != mClass) {
      if (theMat[i] != 0)
        delete theMat[i];
      theMat[i] = theBroker.getNewUniaxialMaterial(mClass);
      if (theMat[i] == 0) {
        opserr << "LadrunoIMKBeam2d::recvSelf -- could not create material " << i << endln;
        return -1;
      }
    }
    theMat[i]->setDbTag(mDb);
    res += theMat[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "LadrunoIMKBeam2d::recvSelf -- failed to recv material " << i << endln;
      return res;
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// output
// ---------------------------------------------------------------------------
void LadrunoIMKBeam2d::Print(OPS_Stream &s, int flag)
{
  s << "LadrunoIMKBeam2d, tag: " << this->getTag() << endln;
  s << "  nodes: " << connectedExternalNodes(0) << " "
    << connectedExternalNodes(1) << endln;
  s << "  A: " << A << " E: " << E << " Iz: " << Iz << endln;
  const char *lbl[2] = { "Mz@i", "Mz@j" };
  for (int i = 0; i < 2; i++)
    s << "  hinge " << lbl[i] << ": "
      << (theMat[i] != 0 ? "material " : "elastic")
      << (theMat[i] != 0 ? theMat[i]->getTag() : 0) << endln;
}

// ---------------------------------------------------------------------------
// recorders
// ---------------------------------------------------------------------------
Response *LadrunoIMKBeam2d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  Response *theResponse = 0;

  s.tag("ElementOutput");
  s.attr("eleType", "LadrunoIMKBeam2d");
  s.attr("eleTag", this->getTag());
  s.attr("node1", connectedExternalNodes(0));
  s.attr("node2", connectedExternalNodes(1));

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "localForce") == 0) {
    theResponse = new ElementResponse(this, 1, P);

  } else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {
    theResponse = new ElementResponse(this, 2, Vector(3));

  } else if (strcmp(argv[0], "basicDeformation") == 0 ||
             strcmp(argv[0], "deformation") == 0 ||
             strcmp(argv[0], "deformations") == 0) {
    theResponse = new ElementResponse(this, 3, Vector(3));

  } else if (strcmp(argv[0], "hingeRotation") == 0 ||
             strcmp(argv[0], "hingeRotations") == 0 ||
             strcmp(argv[0], "plasticRotation") == 0) {
    theResponse = new ElementResponse(this, 4, Vector(2));

  } else if (strcmp(argv[0], "hingeMoment") == 0 ||
             strcmp(argv[0], "hingeMoments") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));

  } else if ((strcmp(argv[0], "material") == 0 || strcmp(argv[0], "hinge") == 0) &&
             argc > 2) {
    int idx = atoi(argv[1]);
    if (idx >= 0 && idx < 2 && theMat[idx] != 0)
      theResponse = theMat[idx]->setResponse(&argv[2], argc - 2, s);
  }

  s.endTag();
  return theResponse;
}

int LadrunoIMKBeam2d::getResponse(int responseID, Information &info)
{
  switch (responseID) {
  case 1:
    return info.setVector(this->getResistingForce());
  case 2:
    return info.setVector(q);
  case 3:
    return info.setVector(theCoordTransf->getBasicTrialDisp());
  case 4: {
    static Vector th(2);
    for (int i = 0; i < 2; i++)
      th(i) = thetaH[i];
    return info.setVector(th);
  }
  case 5: {
    static Vector m(2);
    m(0) = q(1);  m(1) = q(2);
    return info.setVector(m);
  }
  default:
    return -1;
  }
}
