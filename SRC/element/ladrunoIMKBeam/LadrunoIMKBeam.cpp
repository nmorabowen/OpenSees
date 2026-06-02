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
// LadrunoIMKBeam implementation. See LadrunoIMKBeam.h and
// Ladruno_implementation/14_ladruno_imk_beam.md for the formulation.

#include "LadrunoIMKBeam.h"

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

Matrix LadrunoIMKBeam::K(12, 12);
Vector LadrunoIMKBeam::P(12);

// ---------------------------------------------------------------------------
// constructors / destructor
// ---------------------------------------------------------------------------
LadrunoIMKBeam::LadrunoIMKBeam()
  : Element(0, ELE_TAG_LadrunoIMKBeam),
    A(0.0), E(0.0), G(0.0), Jx(0.0), Iy(0.0), Iz(0.0), rho(0.0),
    q(6), kb(6, 6), Q(12),
    connectedExternalNodes(2), theCoordTransf(0)
{
  for (int i = 0; i < 4; i++) {
    theMat[i] = 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }
  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam::LadrunoIMKBeam(int tag, int Nd1, int Nd2,
                               double a, double e, double g, double jx,
                               double iy, double iz, CrdTransf &theTransf,
                               UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                               UniaxialMaterial *mYi, UniaxialMaterial *mYj,
                               double r)
  : Element(tag, ELE_TAG_LadrunoIMKBeam),
    A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), rho(r),
    q(6), kb(6, 6), Q(12),
    connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  UniaxialMaterial *src[4] = { mZi, mZj, mYi, mYj };
  for (int i = 0; i < 4; i++) {
    theMat[i] = (src[i] != 0) ? src[i]->getCopy() : 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }

  theCoordTransf = theTransf.getCopy3d();
  if (theCoordTransf == 0) {
    opserr << "LadrunoIMKBeam::LadrunoIMKBeam -- failed to get copy of CrdTransf\n";
    exit(-1);
  }

  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam::~LadrunoIMKBeam()
{
  for (int i = 0; i < 4; i++)
    if (theMat[i] != 0)
      delete theMat[i];
  if (theCoordTransf != 0)
    delete theCoordTransf;
}

// ---------------------------------------------------------------------------
// trivial accessors
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::getNumExternalNodes(void) const { return 2; }
const ID &LadrunoIMKBeam::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoIMKBeam::getNodePtrs(void) { return theNodes; }
int LadrunoIMKBeam::getNumDOF(void) { return 12; }

void LadrunoIMKBeam::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }

  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

  if (theNodes[0] == 0 || theNodes[1] == 0) {
    opserr << "LadrunoIMKBeam::setDomain -- node not found for element "
           << this->getTag() << endln;
    return;
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  if (dofNd1 != 6 || dofNd2 != 6) {
    opserr << "LadrunoIMKBeam::setDomain -- both nodes must have 6 DOF (element "
           << this->getTag() << ")\n";
    return;
  }

  this->DomainComponent::setDomain(theDomain);

  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "LadrunoIMKBeam::setDomain -- CrdTransf failed to initialize\n";
    return;
  }

  double L = theCoordTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "LadrunoIMKBeam::setDomain -- element " << this->getTag()
           << " has zero length\n";
    return;
  }
}

// ---------------------------------------------------------------------------
// state cycle
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::commitState(void)
{
  int ok = 0;
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->commitState();
    thetaHcommit[i] = thetaH[i];
  }
  ok += theCoordTransf->commitState();
  return ok;
}

int LadrunoIMKBeam::revertToLastCommit(void)
{
  int ok = 0;
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToLastCommit();
    thetaH[i] = thetaHcommit[i];
  }
  ok += theCoordTransf->revertToLastCommit();
  return ok;
}

int LadrunoIMKBeam::revertToStart(void)
{
  int ok = 0;
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToStart();
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
  }
  ok += theCoordTransf->revertToStart();
  return ok;
}

// ---------------------------------------------------------------------------
// per-axis internal state determination (2x2 Newton on hinge rotations)
// ---------------------------------------------------------------------------
void LadrunoIMKBeam::solveAxis(double vi, double vj, int mi, int mj,
                               double L, double EI,
                               double &thi, double &thj,
                               double &qi, double &qj, double k2[2][2])
{
  const double a = 4.0 * EI / L;   // elastic stiffness block  K_el = [[a,b],[b,a]]
  const double b = 2.0 * EI / L;
  const double ktRef = a;          // stiffness scale for guards
  const double ktFloor = 1.0e-8 * ktRef;

  const bool hi = (theMat[mi] != 0);
  const bool hj = (theMat[mj] != 0);

  // warm-start from the incoming trial values; non-hinged ends stay at 0
  double th_i = hi ? thi : 0.0;
  double th_j = hj ? thj : 0.0;

  double Mi = 0.0, Mj = 0.0, kti = 0.0, ktj = 0.0;

  if (hi || hj) {
    const int maxIter = 25;
    for (int it = 0; it < maxIter; it++) {
      if (hi) {
        theMat[mi]->setTrialStrain(th_i);
        Mi = theMat[mi]->getStress();
        kti = theMat[mi]->getTangent();
      }
      if (hj) {
        theMat[mj]->setTrialStrain(th_j);
        Mj = theMat[mj]->getStress();
        ktj = theMat[mj]->getTangent();
      }

      // basic forces from elastic interior: q = K_el (v - thetaH)
      const double qi_ = a * (vi - th_i) + b * (vj - th_j);
      const double qj_ = b * (vi - th_i) + a * (vj - th_j);

      // residual: hinged-end basic force must equal the hinge moment
      const double gi = hi ? (qi_ - Mi) : 0.0;
      const double gj = hj ? (qj_ - Mj) : 0.0;

      const double gnorm = sqrt(gi * gi + gj * gj);
      const double ref = fabs(qi_) + fabs(qj_) + fabs(Mi) + fabs(Mj);
      if (gnorm <= 1.0e-10 * ref + 1.0e-12)
        break;

      // guarded tangents for the Jacobian (force residual stays exact)
      double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
      double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;

      if (hi && hj) {
        // J = -(K_el_sub + diag(kt))
        const double J11 = -a - ki, J12 = -b;
        const double J21 = -b,       J22 = -a - kj;
        double det = J11 * J22 - J12 * J21;
        if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
        // dth = -J^{-1} g
        const double dthi = -(J22 * gi - J12 * gj) / det;
        const double dthj = -(-J21 * gi + J11 * gj) / det;
        th_i += dthi;
        th_j += dthj;
      } else if (hi) {
        const double J11 = -a - ki;
        th_i += -gi / J11;
      } else { // hj only
        const double J22 = -a - kj;
        th_j += -gj / J22;
      }
    }
  }

  thi = th_i;
  thj = th_j;

  // exact basic forces at the converged hinge rotations
  qi = a * (vi - th_i) + b * (vj - th_j);
  qj = b * (vi - th_i) + a * (vj - th_j);

  // condensed tangent block:  K = (F_el + F_h)^{-1}
  const double f = L / (6.0 * EI);     // F_el = f * [[2,-1],[-1,2]]
  double F11 = 2.0 * f, F12 = -f;
  double F21 = -f,       F22 = 2.0 * f;
  if (hi) {
    double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
    F11 += 1.0 / ki;
  }
  if (hj) {
    double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;
    F22 += 1.0 / kj;
  }
  double det = F11 * F22 - F12 * F21;
  if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
  k2[0][0] = F22 / det;
  k2[0][1] = -F12 / det;
  k2[1][0] = -F21 / det;
  k2[1][1] = F11 / det;
}

void LadrunoIMKBeam::initBlock(int mi, int mj, double L, double EI, double k2[2][2])
{
  const double a = 4.0 * EI / L;
  const double ktRef = a;
  const double ktFloor = 1.0e-8 * ktRef;
  const double f = L / (6.0 * EI);
  double F11 = 2.0 * f, F12 = -f;
  double F21 = -f,       F22 = 2.0 * f;
  if (theMat[mi] != 0) {
    double k0 = theMat[mi]->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F11 += 1.0 / k0;
  }
  if (theMat[mj] != 0) {
    double k0 = theMat[mj]->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F22 += 1.0 / k0;
  }
  double det = F11 * F22 - F12 * F21;
  if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
  k2[0][0] = F22 / det;
  k2[0][1] = -F12 / det;
  k2[1][0] = -F21 / det;
  k2[1][1] = F11 / det;
}

// ---------------------------------------------------------------------------
// update: full basic-system state determination, caches q and kb
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::update(void)
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;
  const double EAoverL = E * A * oneOverL;
  const double GJoverL = G * Jx * oneOverL;

  kb.Zero();
  q.Zero();

  // axial + torsion: linear elastic, no hinge
  kb(0, 0) = EAoverL;  q(0) = EAoverL * v(0);
  kb(5, 5) = GJoverL;  q(5) = GJoverL * v(5);

  double k2[2][2];

  // strong-axis bending: basic dofs 1,2 ; materials 0,1 ; EI = E*Iz
  double qi, qj;
  solveAxis(v(1), v(2), 0, 1, L, E * Iz, thetaH[0], thetaH[1], qi, qj, k2);
  q(1) = qi;  q(2) = qj;
  kb(1, 1) = k2[0][0];  kb(1, 2) = k2[0][1];
  kb(2, 1) = k2[1][0];  kb(2, 2) = k2[1][1];

  // weak-axis bending: basic dofs 3,4 ; materials 2,3 ; EI = E*Iy
  solveAxis(v(3), v(4), 2, 3, L, E * Iy, thetaH[2], thetaH[3], qi, qj, k2);
  q(3) = qi;  q(4) = qj;
  kb(3, 3) = k2[0][0];  kb(3, 4) = k2[0][1];
  kb(4, 3) = k2[1][0];  kb(4, 4) = k2[1][1];

  return ok;
}

const Matrix &LadrunoIMKBeam::getTangentStiff(void)
{
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &LadrunoIMKBeam::getInitialStiff(void)
{
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;

  static Matrix kbi(6, 6);
  kbi.Zero();
  kbi(0, 0) = E * A * oneOverL;
  kbi(5, 5) = G * Jx * oneOverL;

  double k2[2][2];
  initBlock(0, 1, L, E * Iz, k2);
  kbi(1, 1) = k2[0][0];  kbi(1, 2) = k2[0][1];
  kbi(2, 1) = k2[1][0];  kbi(2, 2) = k2[1][1];
  initBlock(2, 3, L, E * Iy, k2);
  kbi(3, 3) = k2[0][0];  kbi(3, 4) = k2[0][1];
  kbi(4, 3) = k2[1][0];  kbi(4, 4) = k2[1][1];

  return theCoordTransf->getInitialGlobalStiffMatrix(kbi);
}

const Matrix &LadrunoIMKBeam::getMass(void)
{
  K.Zero();
  if (rho > 0.0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;
    K(0, 0) = K(1, 1) = K(2, 2) = m;
    K(6, 6) = K(7, 7) = K(8, 8) = m;
  }
  return K;
}

// ---------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------
void LadrunoIMKBeam::zeroLoad(void) { Q.Zero(); }

int LadrunoIMKBeam::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "LadrunoIMKBeam::addLoad -- element loads not supported in v1 "
            "(element " << this->getTag() << "); use nodal loads\n";
  return -1;
}

int LadrunoIMKBeam::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  double L = theCoordTransf->getInitialLength();
  double m = 0.5 * rho * L;

  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  Q(2) -= m * Raccel1(2);
  Q(6) -= m * Raccel2(0);
  Q(7) -= m * Raccel2(1);
  Q(8) -= m * Raccel2(2);

  return 0;
}

// ---------------------------------------------------------------------------
// resisting force
// ---------------------------------------------------------------------------
const Vector &LadrunoIMKBeam::getResistingForce(void)
{
  static Vector p0Vec(5);
  p0Vec.Zero();

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // subtract element nodal loads (inertia unbalance accumulated in Q)
  P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector &LadrunoIMKBeam::getResistingForceIncInertia(void)
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
  P(2) += m * accel1(2);
  P(6) += m * accel2(0);
  P(7) += m * accel2(1);
  P(8) += m * accel2(2);

  return P;
}

// ---------------------------------------------------------------------------
// parallel
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();

  // integer data: tag, nodes, transf class+db, 4x (mat class+db)
  static ID idData(13);
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

  for (int i = 0; i < 4; i++) {
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
    opserr << "LadrunoIMKBeam::sendSelf -- failed to send ID\n";
    return res;
  }

  static Vector dData(11);
  dData(0) = A;  dData(1) = E;  dData(2) = G;  dData(3) = Jx;
  dData(4) = Iy; dData(5) = Iz; dData(6) = rho;
  dData(7) = thetaHcommit[0];  dData(8) = thetaHcommit[1];
  dData(9) = thetaHcommit[2];  dData(10) = thetaHcommit[3];

  res += theChannel.sendVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::sendSelf -- failed to send Vector\n";
    return res;
  }

  res += theCoordTransf->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::sendSelf -- failed to send CrdTransf\n";
    return res;
  }

  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0) {
      res += theMat[i]->sendSelf(cTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoIMKBeam::sendSelf -- failed to send material " << i << endln;
        return res;
      }
    }
  }

  return res;
}

int LadrunoIMKBeam::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static ID idData(13);
  res += theChannel.recvID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::recvSelf -- failed to recv ID\n";
    return res;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);

  static Vector dData(11);
  res += theChannel.recvVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::recvSelf -- failed to recv Vector\n";
    return res;
  }
  A = dData(0);  E = dData(1);  G = dData(2);  Jx = dData(3);
  Iy = dData(4); Iz = dData(5); rho = dData(6);
  for (int i = 0; i < 4; i++) {
    thetaHcommit[i] = dData(7 + i);
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
      opserr << "LadrunoIMKBeam::recvSelf -- could not create CrdTransf\n";
      return -1;
    }
  }
  theCoordTransf->setDbTag(crdDb);
  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::recvSelf -- failed to recv CrdTransf\n";
    return res;
  }

  // materials
  for (int i = 0; i < 4; i++) {
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
        opserr << "LadrunoIMKBeam::recvSelf -- could not create material " << i << endln;
        return -1;
      }
    }
    theMat[i]->setDbTag(mDb);
    res += theMat[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "LadrunoIMKBeam::recvSelf -- failed to recv material " << i << endln;
      return res;
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// output
// ---------------------------------------------------------------------------
void LadrunoIMKBeam::Print(OPS_Stream &s, int flag)
{
  s << "LadrunoIMKBeam, tag: " << this->getTag() << endln;
  s << "  nodes: " << connectedExternalNodes(0) << " "
    << connectedExternalNodes(1) << endln;
  s << "  A: " << A << " E: " << E << " G: " << G
    << " Jx: " << Jx << " Iy: " << Iy << " Iz: " << Iz << endln;
  const char *lbl[4] = { "Mz@i", "Mz@j", "My@i", "My@j" };
  for (int i = 0; i < 4; i++)
    s << "  hinge " << lbl[i] << ": "
      << (theMat[i] != 0 ? "material " : "elastic")
      << (theMat[i] != 0 ? theMat[i]->getTag() : 0) << endln;
}

// ---------------------------------------------------------------------------
// recorders
// ---------------------------------------------------------------------------
Response *LadrunoIMKBeam::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  Response *theResponse = 0;

  s.tag("ElementOutput");
  s.attr("eleType", "LadrunoIMKBeam");
  s.attr("eleTag", this->getTag());
  s.attr("node1", connectedExternalNodes(0));
  s.attr("node2", connectedExternalNodes(1));

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "localForce") == 0) {
    theResponse = new ElementResponse(this, 1, P);

  } else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {
    theResponse = new ElementResponse(this, 2, Vector(6));

  } else if (strcmp(argv[0], "basicDeformation") == 0 ||
             strcmp(argv[0], "deformation") == 0 ||
             strcmp(argv[0], "deformations") == 0) {
    theResponse = new ElementResponse(this, 3, Vector(6));

  } else if (strcmp(argv[0], "hingeRotation") == 0 ||
             strcmp(argv[0], "hingeRotations") == 0 ||
             strcmp(argv[0], "plasticRotation") == 0) {
    theResponse = new ElementResponse(this, 4, Vector(4));

  } else if (strcmp(argv[0], "hingeMoment") == 0 ||
             strcmp(argv[0], "hingeMoments") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(4));

  } else if ((strcmp(argv[0], "material") == 0 || strcmp(argv[0], "hinge") == 0) &&
             argc > 2) {
    int idx = atoi(argv[1]);
    if (idx >= 0 && idx < 4 && theMat[idx] != 0)
      theResponse = theMat[idx]->setResponse(&argv[2], argc - 2, s);
  }

  s.endTag();
  return theResponse;
}

int LadrunoIMKBeam::getResponse(int responseID, Information &info)
{
  switch (responseID) {
  case 1:
    return info.setVector(this->getResistingForce());
  case 2:
    return info.setVector(q);
  case 3:
    return info.setVector(theCoordTransf->getBasicTrialDisp());
  case 4: {
    static Vector th(4);
    for (int i = 0; i < 4; i++)
      th(i) = thetaH[i];
    return info.setVector(th);
  }
  case 5: {
    static Vector m(4);
    m(0) = q(1);  m(1) = q(2);  m(2) = q(3);  m(3) = q(4);
    return info.setVector(m);
  }
  default:
    return -1;
  }
}
