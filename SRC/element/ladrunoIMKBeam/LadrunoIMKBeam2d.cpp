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
#include <LadrunoResponseTokens.h>   // Ladruno — shared recorder-token aliases
#include <ElementalLoad.h>
#include <Renderer.h>
#include <classTags.h>
#include <OPS_Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static thread_local Matrix K(6, 6);
static thread_local Vector P(6);

// ---------------------------------------------------------------------------
// constructors / destructor
// ---------------------------------------------------------------------------
LadrunoIMKBeam2d::LadrunoIMKBeam2d()
  : Element(0, ELE_TAG_LadrunoIMKBeam2d),
    A(0.0), E(0.0), Iz(0.0), rho(0.0), cMass(0),
    GAvy(0.0), P0(0.0), pmAlpha(0.5),
    dI(0.0), dJ(0.0), useImplex(false), implexErr(0.0),
    dtime_n_commit(0.0),
    q(3), kb(3, 3), Q(6), q0Dirty(false), committedAtCurrentV(false),
    connectedExternalNodes(2), theCoordTransf(0)
{
  for (int i = 0; i < 2; i++) {
    theMat[i] = 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
    dThetaHprev[i] = 0.0;
    thetaHp[i] = 0.0;
    thetaHpCommit[i] = 0.0;
    dThetaHpPrev[i] = 0.0;
  }
  q0[0] = q0[1] = q0[2] = 0.0;
  p0[0] = p0[1] = p0[2] = 0.0;
  vCommit[0] = vCommit[1] = vCommit[2] = 0.0;
  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam2d::LadrunoIMKBeam2d(int tag, int Nd1, int Nd2,
                                   double a, double e, double iz,
                                   CrdTransf &theTransf,
                                   UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                                   double r, int cm,
                                   double gavy, double p0_pm, double alpha_pm,
                                   double di, double dj, bool implex)
  : Element(tag, ELE_TAG_LadrunoIMKBeam2d),
    A(a), E(e), Iz(iz), rho(r), cMass(cm),
    GAvy(gavy), P0(p0_pm), pmAlpha(alpha_pm),
    dI(di), dJ(dj), useImplex(implex), implexErr(0.0),
    dtime_n_commit(0.0),
    q(3), kb(3, 3), Q(6), q0Dirty(false), committedAtCurrentV(false),
    connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  UniaxialMaterial *src[2] = { mZi, mZj };
  for (int i = 0; i < 2; i++) {
    if (src[i] != 0) {
      theMat[i] = src[i]->getCopy();
      if (theMat[i] == 0) {
        opserr << "LadrunoIMKBeam2d::LadrunoIMKBeam2d -- failed to copy material " << i << endln;
        exit(-1);
      }
    } else {
      theMat[i] = 0;
    }
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
    dThetaHprev[i] = 0.0;
    thetaHp[i] = 0.0;
    thetaHpCommit[i] = 0.0;
    dThetaHpPrev[i] = 0.0;
  }

  theCoordTransf = theTransf.getCopy2d();
  if (theCoordTransf == 0) {
    opserr << "LadrunoIMKBeam2d::LadrunoIMKBeam2d -- failed to get copy of CrdTransf\n";
    exit(-1);
  }

  q0[0] = q0[1] = q0[2] = 0.0;
  p0[0] = p0[1] = p0[2] = 0.0;
  vCommit[0] = vCommit[1] = vCommit[2] = 0.0;

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

  this->update();
}

// ---------------------------------------------------------------------------
// state cycle
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::commitState(void)
{
  // Base-class commit first: it refreshes Kc, the committed stiffness behind
  // betaKc Rayleigh damping. Skipping it left Kc frozen at the tangent captured
  // when `rayleigh` ran, so betaKc silently behaved like betaK0 once hinges
  // yielded or the geometry rotated (WP-118). Same order as ElasticBeam2d.  // Ladruno
  int ok = this->Element::commitState();
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->commitState();
    dThetaHprev[i] = thetaH[i] - thetaHcommit[i];
    thetaHcommit[i] = thetaH[i];
    dThetaHpPrev[i] = thetaHp[i] - thetaHpCommit[i];
    thetaHpCommit[i] = thetaHp[i];
  }
  dtime_n_commit = ops_Dt;
  if (theCoordTransf != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    const Vector &vc = theCoordTransf->getBasicTrialDisp();
    vCommit[0] = vc(0);
    vCommit[1] = vc(1);
    vCommit[2] = vc(2);
  }
  ok += theCoordTransf->commitState();
  return ok;
}

int LadrunoIMKBeam2d::revertToLastCommit(void)
{
  int ok = 0;
  committedAtCurrentV = false;
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToLastCommit();
    thetaH[i] = thetaHcommit[i];
    thetaHp[i] = thetaHpCommit[i];
  }
  ok += theCoordTransf->revertToLastCommit();
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0)
    ok += this->update();
  return ok;
}

int LadrunoIMKBeam2d::revertToStart(void)
{
  int ok = 0;
  committedAtCurrentV = false;
  vCommit[0] = vCommit[1] = vCommit[2] = 0.0;
  for (int i = 0; i < 2; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToStart();
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
    dThetaHprev[i] = 0.0;
    thetaHp[i] = 0.0;
    thetaHpCommit[i] = 0.0;
    dThetaHpPrev[i] = 0.0;
  }
  dtime_n_commit = 0.0;
  implexErr = 0.0;
  ok += theCoordTransf->revertToStart();
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0)
    ok += this->update();
  return ok;
}

// ---------------------------------------------------------------------------
// update: full basic-system state determination, caches q and kb
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::update(void)
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();
  if (committedAtCurrentV) {
    if (!q0Dirty &&
        v(0) == vCommit[0] && v(1) == vCommit[1] && v(2) == vCommit[2]) {
      return ok;
    }
    committedAtCurrentV = false;
  }
  const double L = theCoordTransf->getInitialLength();
  if (L <= 0.0)
    return -1;
  const double oneOverL = 1.0 / L;
  const double EAoverL = E * A * oneOverL;

  kb.Zero();
  q.Zero();

  // axial: linear elastic, no hinge (+ fixed-end axial force q0[0])
  kb(0, 0) = EAoverL;  q(0) = EAoverL * v(0) + q0[0];

  const double fPM = ladrunoIMKScalePM(q(0), P0, pmAlpha);
  if (P0 <= 0.0) {
    for (int i = 0; i < 2; i++)
      if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
        ch->setAxialForce(q(0));
  }
  double k2[2][2];

  double dtScale = 1.0;
  if (dtime_n_commit > 1.0e-16 && ops_Dt > 0.0) {
    dtScale = ops_Dt / dtime_n_commit;
    if (dtScale > 2.0) dtScale = 2.0;
  }

  // strong-axis bending: basic dofs 1,2 ; materials 0,1 ; EI = E*Iz
  double qi, qj;
  ladrunoIMKSolveAxis(theMat[0], theMat[1], v(1), v(2), L, E * Iz,
                      thetaH[0], thetaH[1], qi, qj, k2, q0[1], q0[2],
                      GAvy, fPM, dI, dJ, useImplex,
                      thetaHpCommit[0], thetaHpCommit[1],
                      dThetaHpPrev[0], dThetaHpPrev[1], &implexErr,
                      dtScale, &thetaHp[0], &thetaHp[1]);
  q(1) = qi;  q(2) = qj;
  kb(1, 1) = k2[0][0];  kb(1, 2) = k2[0][1];
  kb(2, 1) = k2[1][0];  kb(2, 2) = k2[1][1];
  q0Dirty = false;

  return ok;
}

const Matrix &LadrunoIMKBeam2d::getTangentStiff(void)
{
  if (q0Dirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &LadrunoIMKBeam2d::getInitialStiff(void)
{
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;

  static thread_local Matrix kbi(3, 3);
  kbi.Zero();
  kbi(0, 0) = E * A * oneOverL;

  double k2[2][2];
  ladrunoIMKInitBlock(theMat[0], theMat[1], L, E * Iz, k2, GAvy, 1.0, dI, dJ);
  kbi(1, 1) = k2[0][0];  kbi(1, 2) = k2[0][1];
  kbi(2, 1) = k2[1][0];  kbi(2, 2) = k2[1][1];

  return theCoordTransf->getInitialGlobalStiffMatrix(kbi);
}

const Matrix &LadrunoIMKBeam2d::getMass(void)
{
  K.Zero();
  if (rho <= 0.0)
    return K;

  double L = theCoordTransf->getInitialLength();
  if (cMass == 0) {
    double m = 0.5 * rho * L;
    K(0, 0) = K(1, 1) = m;   // node i translations
    K(3, 3) = K(4, 4) = m;   // node j translations
  } else {
    static thread_local Matrix ml(6, 6);
    ml.Zero();
    double m = rho * L / 420.0;
    ml(0, 0) = ml(3, 3) = m * 140.0;
    ml(0, 3) = ml(3, 0) = m * 70.0;

    ml(1, 1) = ml(4, 4) = m * 156.0;
    ml(1, 4) = ml(4, 1) = m * 54.0;
    ml(2, 2) = ml(5, 5) = m * 4.0 * L * L;
    ml(2, 5) = ml(5, 2) = -m * 3.0 * L * L;
    ml(1, 2) = ml(2, 1) = m * 22.0 * L;
    ml(4, 5) = ml(5, 4) = -ml(1, 2);
    ml(1, 5) = ml(5, 1) = -m * 13.0 * L;
    ml(2, 4) = ml(4, 2) = -ml(1, 5);

    K = theCoordTransf->getGlobalMatrixFromLocal(ml);
  }
  return K;
}

// ---------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------
void LadrunoIMKBeam2d::zeroLoad(void)
{
  Q.Zero();
  if (q0[0] != 0.0 || q0[1] != 0.0 || q0[2] != 0.0)
    q0Dirty = true;
  q0[0] = q0[1] = q0[2] = 0.0;
  p0[0] = p0[1] = p0[2] = 0.0;
}

int LadrunoIMKBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0) * loadFactor;  // Transverse (+ve upward)
    double wa = data(1) * loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5 * wt * L;
    double M = V * L / 6.0;  // wt*L*L/12
    double P = wa * L;

    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    q0[0] -= 0.5 * P;
    q0[1] -= M;
    q0[2] += M;
    q0Dirty = true;
  }
  else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
    double wta = data(0) * loadFactor;  // Transverse at start
    double wtb = data(1) * loadFactor;  // Transverse at end
    double waa = data(2) * loadFactor;  // Axial at start
    double wab = data(3) * loadFactor;  // Axial at end
    double a   = data(4) * L;
    double b   = data(5) * L;

    double ell = b - a;
    if (ell <= 0.0 || a < 0.0 || b > L)
      return 0;

    double P  = 0.5 * ell * (waa + wab);
    double PJ = (a * P + (ell * ell / 6.0) * (waa + 2.0 * wab)) / L;

    double J0 = 0.5 * ell * (wta + wtb);
    double J1 = (ell * ell / 6.0) * (wta + 2.0 * wtb);
    double J2 = (ell * ell * ell / 12.0) * (wta + 3.0 * wtb);
    double J3 = (ell * ell * ell * ell / 20.0) * (wta + 4.0 * wtb);

    double I0 = J0;
    double I1 = a * J0 + J1;
    double I2 = a * a * J0 + 2.0 * a * J1 + J2;
    double I3 = a * a * a * J0 + 3.0 * a * a * J1 + 3.0 * a * J2 + J3;

    double V2 = I1 / L;
    double V1 = I0 - V2;
    double M1 = -I1 + (2.0 / L) * I2 - (1.0 / (L * L)) * I3;
    double M2 = (1.0 / L) * I2 - (1.0 / (L * L)) * I3;
    if (GAvy > 0.0 && L > 0.0) {
      const double Phi = 12.0 * E * Iz / (GAvy * L * L);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1 - 0.5 * Phi * M2) * inv;
      const double m2 = (-0.5 * Phi * M1 + (1.0 + 0.5 * Phi) * M2) * inv;
      M1 = m1;
      M2 = m2;
    }

    p0[0] -= P;
    p0[1] -= V1;
    p0[2] -= V2;

    q0[0] -= PJ;
    q0[1] += M1;
    q0[2] += M2;
    q0Dirty = true;
  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0) * loadFactor;
    double N = data(1) * loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL * L;
    double b = L - a;

    p0[0] -= N;
    double V1 = P * (1.0 - aOverL);
    double V2 = P * aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0 / (L * L);
    double a2 = a * a;
    double b2 = b * b;

    q0[0] -= N * aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    if (GAvy > 0.0 && L > 0.0) {
      const double Phi = 12.0 * E * Iz / (GAvy * L * L);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1 - 0.5 * Phi * M2) * inv;
      const double m2 = (-0.5 * Phi * M1 + (1.0 + 0.5 * Phi) * M2) * inv;
      M1 = m1;
      M2 = m2;
    }
    q0[1] += M1;
    q0[2] += M2;
    q0Dirty = true;
  }
  else if (type == LOAD_TAG_BeamUniformMoment) {
    double mz = data(2) * loadFactor;
    p0[1] += mz;
    p0[2] -= mz;
  }
  else {
    opserr << "LadrunoIMKBeam2d::addLoad -- load type unknown for element with tag: "
           << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int LadrunoIMKBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "LadrunoIMKBeam2d::addInertiaLoadToUnbalance -- incompatible node DOF size\n";
    return -1;
  }

  if (cMass == 0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(3) -= m * Raccel2(0);
    Q(4) -= m * Raccel2(1);
  } else {
    static thread_local Vector Raccel(6);
    for (int i = 0; i < 3; i++) {
      Raccel(i)     = Raccel1(i);
      Raccel(i + 3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

// ---------------------------------------------------------------------------
// resisting force
// ---------------------------------------------------------------------------
const Vector &LadrunoIMKBeam2d::getResistingForce(void)
{
  // Re-evaluate basic state only when elemental loads q0 were modified since the
  // last update() (e.g. under Linear algorithm before Domain::update() runs).
  if (q0Dirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();

  Vector p0Vec(p0, 3);

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // subtract element nodal loads (inertia unbalance accumulated in Q)
  if (rho != 0.0)
    P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector &LadrunoIMKBeam2d::getResistingForceIncInertia(void)
{
  // SNAPSHOT into a function-local buffer BEFORE getRayleighDampingForces():
  // betaK Rayleigh re-enters getTangentStiff(); returning `res` keeps the result
  // independent of anything that path writes, and of the static P that
  // getResistingForce() returns (Element::getResponse 444444 subtracts both).
  // Order is unchanged -- ((f - Q) + R) + m*a -- so results are bit-identical
  // (WP-115; LEDGER_quirks "MUST snapshot the shared static `resid`").  // Ladruno
  static thread_local Vector res(6);
  res = this->getResistingForce();

  // Rayleigh damping forces (added BEFORE inertia -- keep this order)
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  if (rho == 0.0)
    return res;

  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();

  if (cMass == 0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    res(0) += m * accel1(0);
    res(1) += m * accel1(1);
    res(3) += m * accel2(0);
    res(4) += m * accel2(1);
  } else {
    static thread_local Vector accel(6);
    for (int i = 0; i < 3; i++) {
      accel(i)     = accel1(i);
      accel(i + 3) = accel2(i);
    }
    res.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }

  return res;
}

// ---------------------------------------------------------------------------
// parallel
// ---------------------------------------------------------------------------
int LadrunoIMKBeam2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();

  // integer data: tag, nodes, transf class+db, 2x (mat class+db), cMass
  static thread_local ID idData(10);
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
  idData(9) = cMass;

  res += theChannel.sendID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::sendSelf -- failed to send ID\n";
    return res;
  }

  static thread_local Vector dData(45);
  dData(0) = A;  dData(1) = E;  dData(2) = Iz;  dData(3) = rho;
  dData(4) = thetaHcommit[0];  dData(5) = thetaHcommit[1];
  dData(6) = alphaM;  dData(7) = betaK;  dData(8) = betaK0;  dData(9) = betaKc;
  dData(10) = GAvy;   dData(11) = P0;   dData(12) = pmAlpha;
  dData(13) = dI;     dData(14) = dJ;   dData(15) = useImplex ? 1.0 : 0.0;
  dData(16) = thetaHpCommit[0];  dData(17) = thetaHpCommit[1];
  dData(18) = dThetaHpPrev[0];   dData(19) = dThetaHpPrev[1];
  dData(20) = dThetaHprev[0];    dData(21) = dThetaHprev[1];
  dData(22) = dtime_n_commit;
  dData(23) = implexErr;
  for (int i = 0; i < 3; i++)
    dData(24 + i) = q(i);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      dData(27 + 3 * i + j) = kb(i, j);
  for (int i = 0; i < 3; i++)
    dData(36 + i) = vCommit[i];
  for (int i = 0; i < 3; i++) {
    dData(39 + i) = q0[i];
    dData(42 + i) = p0[i];
  }

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

  static thread_local ID idData(10);
  res += theChannel.recvID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam2d::recvSelf -- failed to recv ID\n";
    return res;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  cMass = idData(9);

  static thread_local Vector dData(45);
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
  alphaM = dData(6);  betaK = dData(7);  betaK0 = dData(8);  betaKc = dData(9);
  GAvy = dData(10);   P0 = dData(11);    pmAlpha = dData(12);
  dI = dData(13);     dJ = dData(14);    useImplex = (dData(15) != 0.0);
  for (int i = 0; i < 2; i++) {
    thetaHpCommit[i] = dData(16 + i);
    thetaHp[i] = thetaHpCommit[i];
    dThetaHpPrev[i] = dData(18 + i);
    dThetaHprev[i] = dData(20 + i);
  }
  dtime_n_commit = dData(22);
  implexErr = dData(23);
  for (int i = 0; i < 3; i++)
    q(i) = dData(24 + i);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      kb(i, j) = dData(27 + 3 * i + j);
  for (int i = 0; i < 3; i++)
    vCommit[i] = dData(36 + i);
  for (int i = 0; i < 3; i++) {
    q0[i] = dData(39 + i);
    p0[i] = dData(42 + i);
  }
  q0Dirty = false;
  committedAtCurrentV = true;

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

  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0 && theCoordTransf != 0)
    theCoordTransf->update();

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
  s << "  A: " << A << " E: " << E << " Iz: " << Iz
    << " rho: " << rho << " cMass: " << cMass;
  if (GAvy > 0.0)
    s << " GAvy: " << GAvy;
  if (P0 > 0.0)
    s << " P0: " << P0 << " pmAlpha: " << pmAlpha;
  s << endln;
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

  if (argc < 1) return 0;

  s.tag("ElementOutput");
  s.attr("eleType", "LadrunoIMKBeam2d");
  s.attr("eleTag", this->getTag());
  s.attr("node1", connectedExternalNodes(0));
  s.attr("node2", connectedExternalNodes(1));

  // Ladruno — vocabulary realigned on the DispBeamColumn family: force /
  // globalForce are the GLOBAL end forces, localForce the element-frame ones
  // (it used to be a third spelling of the global vector), chordRotation ==
  // basicDeformation == deformation, plasticRotation == hingeRotation.
  if (LadrunoResp::isGlobalForce(argv[0])) {
    s.tag("ResponseType", "Px_1");
    s.tag("ResponseType", "Py_1");
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "Px_2");
    s.tag("ResponseType", "Py_2");
    s.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 1, P);

  } else if (LadrunoResp::isLocalForce(argv[0])) {
    s.tag("ResponseType", "N_1");
    s.tag("ResponseType", "V_1");
    s.tag("ResponseType", "M_1");
    s.tag("ResponseType", "N_2");
    s.tag("ResponseType", "V_2");
    s.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 6, Vector(6));

  } else if (LadrunoResp::isBasicForce(argv[0])) {
    s.tag("ResponseType", "N");
    s.tag("ResponseType", "M_1");
    s.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 2, Vector(3));

  } else if (strcmp(argv[0], "basicStiffness") == 0) {
    s.tag("ResponseType", "N");
    s.tag("ResponseType", "M1");
    s.tag("ResponseType", "M2");
    theResponse = new ElementResponse(this, 19, Matrix(3, 3));

  } else if (LadrunoResp::isBasicDef(argv[0])) {
    s.tag("ResponseType", "eps");
    s.tag("ResponseType", "theta1");
    s.tag("ResponseType", "theta2");
    theResponse = new ElementResponse(this, 3, Vector(3));

  } else if (LadrunoResp::isPlasticDef(argv[0])) {
    s.tag("ResponseType", "theta1P");
    s.tag("ResponseType", "theta2P");
    theResponse = new ElementResponse(this, 4, Vector(2));

  } else if (LadrunoResp::is(argv[0], "hingeMoment")) {
    s.tag("ResponseType", "M_1");
    s.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 5, Vector(2));

  } else if (LadrunoResp::is(argv[0], "hingeSummary")) {
    s.tag("ResponseType", "theta_i");
    s.tag("ResponseType", "M_i");
    s.tag("ResponseType", "theta_j");
    s.tag("ResponseType", "M_j");
    theResponse = new ElementResponse(this, 8, Vector(4));

  } else if (LadrunoResp::isRayleighForce(argv[0])) {
    theResponse = new ElementResponse(this, 12, P);

  } else if (LadrunoResp::isEnergy(argv[0])) {
    theResponse = new ElementResponse(this, 10, 0.0);

  } else if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
    theResponse = new ElementResponse(this, 32, 0.0);

  } else if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
    theResponse = new ElementResponse(this, 33, 0.0);

  } else if (LadrunoResp::isLocalAxes(argv[0])) {
    theResponse = new ElementResponse(this, 30, Vector(9));

  } else if (strcmp(argv[0], "section") == 0 && argc > 2) {
    int secNum = atoi(argv[1]);
    if ((secNum == 1 || secNum == 2) && LadrunoResp::is(argv[2], "force")) {
      theResponse = new ElementResponse(this, (secNum == 1) ? 71 : 72, Vector(2));
    }
  } else if ((LadrunoResp::is(argv[0], "material") ||
              strcmp(argv[0], "hinge") == 0) && argc > 2) {
    int idx = atoi(argv[1]);
    if (idx >= 0 && idx < 2 && theMat[idx] != 0)
      theResponse = theMat[idx]->setResponse(&argv[2], argc - 2, s);
  }

  if (theResponse == 0 && theCoordTransf != 0)
    theResponse = theCoordTransf->setResponse(argv, argc, s);

  s.endTag();

  // Ladruno — base vocabulary (dampingForce, dynamicForce, inertialForce);
  // Element::setResponse opens its own ElementOutput tag, so this MUST come
  // after endTag().
  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, s);
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
    static thread_local Vector th(2);
    for (int i = 0; i < 2; i++)
      th(i) = thetaH[i];
    return info.setVector(th);
  }
  case 5: {
    static thread_local Vector m(2);
    m(0) = q(1);  m(1) = q(2);
    return info.setVector(m);
  }
  case 6: {
    // Ladruno — local-frame end forces from the basic q = [N, Mz_i, Mz_j],
    // same mapping as DispBeamColumn2d including p0 reactions.
    static thread_local Vector Pl(6);
    const double L = theCoordTransf->getInitialLength();
    const double V = (L > 0.0) ? (q(1) + q(2)) / L : 0.0;
    Pl(0) = -q(0) + p0[0];  Pl(1) =  V + p0[1];  Pl(2) = q(1);
    Pl(3) =  q(0);          Pl(4) = -V + p0[2];  Pl(5) = q(2);
    return info.setVector(Pl);
  }
  case 8: {
    static thread_local Vector hs(4);
    hs(0) = thetaH[0];
    hs(1) = q(1);
    hs(2) = thetaH[1];
    hs(3) = q(2);
    return info.setVector(hs);
  }
  case 10: {
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    const double L = theCoordTransf->getInitialLength();
    double energy = 0.0;
    if (L > 0.0) {
      const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
      const double Lc = hasEZ ? (L - dI - dJ) : L;
      const double T00 = hasEZ ? (L - dJ) / Lc : 1.0;
      const double T01 = hasEZ ? dJ / Lc : 0.0;
      const double T10 = hasEZ ? dI / Lc : 0.0;
      const double T11 = hasEZ ? (L - dI) / Lc : 1.0;
      const double vic = T00 * v(1) + T01 * v(2);
      const double vjc = T10 * v(1) + T11 * v(2);
      const double ea = E * A / L;
      double a, b, F11, F12;
      ladrunoIMKCoreCoeffs(Lc, E * Iz, GAvy, a, b, F11, F12);
      const double ve0 = v(0);
      const double ve1 = vic - thetaH[0];
      const double ve2 = vjc - thetaH[1];
      energy += 0.5 * (ea * ve0 * ve0 + ve1 * (a * ve1 + b * ve2) + ve2 * (b * ve1 + a * ve2));
    }
    const double fPM = ladrunoIMKScalePM(q(0), P0, pmAlpha);
    for (int i = 0; i < 2; i++)
      if (theMat[i] != 0)
        energy += fPM * theMat[i]->getEnergy();
    return info.setDouble(energy);
  }
  case 12: {
    static thread_local Vector damp(6);
    damp.Zero();
    damp.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return info.setVector(damp);
  }
  case 19:
    return info.setMatrix(kb);
  case 30: {
    static thread_local Vector la(9);
    static thread_local Vector vx(3), vy(3), vz(3);
    theCoordTransf->getLocalAxes(vx, vy, vz);
    for (int i = 0; i < 3; i++) {
      la(i)     = vx(i);
      la(i + 3) = vy(i);
      la(i + 6) = vz(i);
    }
    return info.setVector(la);
  }
  case 32:
    return info.setDouble(implexErr);
  case 33: {
    int cnt = (theMat[0] != 0 ? 1 : 0) + (theMat[1] != 0 ? 1 : 0);
    return info.setDouble(cnt > 0 ? implexErr / cnt : 0.0);
  }
  case 71:
  case 72: {
    static thread_local Vector sf(2);
    if (responseID == 71) {
      sf(0) = q(0) - p0[0];
      sf(1) = -q(1);
    } else {
      sf(0) = q(0);
      sf(1) = q(2);
    }
    return info.setVector(sf);
  }
  default:
    return this->Element::getResponse(responseID, info);
  }
}

int LadrunoIMKBeam2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "A") == 0) {
    param.setValue(A);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0], "Iz") == 0 || strcmp(argv[0], "I") == 0) {
    param.setValue(Iz);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(rho);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0], "GAvy") == 0 || strcmp(argv[0], "GAy") == 0 || strcmp(argv[0], "GAv") == 0) {
    param.setValue(GAvy);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0], "P0") == 0) {
    param.setValue(P0);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0], "pmAlpha") == 0) {
    param.setValue(pmAlpha);
    return param.addObject(7, this);
  }

  if (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hinge1") == 0 ||
      strcmp(argv[0], "hingeZi") == 0 || strcmp(argv[0], "hingeZI") == 0 ||
      strcmp(argv[0], "hingeMzI") == 0 || strcmp(argv[0], "matZi") == 0) {
    if (argc < 2 || theMat[0] == 0)
      return -1;
    return theMat[0]->setParameter(&argv[1], argc - 1, param);
  }
  if (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hinge2") == 0 ||
      strcmp(argv[0], "hingeZj") == 0 || strcmp(argv[0], "hingeZJ") == 0 ||
      strcmp(argv[0], "hingeMzJ") == 0 || strcmp(argv[0], "matZj") == 0) {
    if (argc < 2 || theMat[1] == 0)
      return -1;
    return theMat[1]->setParameter(&argv[1], argc - 1, param);
  }
  if (strcmp(argv[0], "matZ") == 0 || strcmp(argv[0], "hingeZ") == 0) {
    if (argc < 2)
      return -1;
    int result = -1;
    for (int i = 0; i < 2; i++) {
      if (theMat[i] != 0) {
        int r = theMat[i]->setParameter(&argv[1], argc - 1, param);
        if (r != -1)
          result = r;
      }
    }
    return result;
  }
  if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "hinge") == 0) {
    if (argc == 2) {
      int result = -1;
      for (int i = 0; i < 2; i++) {
        if (theMat[i] != 0) {
          int r = theMat[i]->setParameter(&argv[1], 1, param);
          if (r != -1)
            result = r;
        }
      }
      return result;
    } else if (argc > 2) {
      int idx = -1;
      if (strcmp(argv[1], "i") == 0 || strcmp(argv[1], "0") == 0)
        idx = 0;
      else if (strcmp(argv[1], "j") == 0 || strcmp(argv[1], "1") == 0)
        idx = 1;
      else
        idx = atoi(argv[1]);
      if (idx >= 0 && idx < 2 && theMat[idx] != 0)
        return theMat[idx]->setParameter(&argv[2], argc - 2, param);
    }
    return -1;
  }

  return -1;
}

int LadrunoIMKBeam2d::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    E = info.theDouble;
    break;
  case 2:
    A = info.theDouble;
    break;
  case 3:
    Iz = info.theDouble;
    break;
  case 4:
    rho = info.theDouble;
    break;
  case 5:
    GAvy = info.theDouble;
    break;
  case 6:
    P0 = info.theDouble;
    break;
  case 7:
    pmAlpha = info.theDouble;
    break;
  default:
    return -1;
  }
  committedAtCurrentV = false;
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  return 0;
}

bool LadrunoIMKBeam2d::ladrunoThreadSafeUpdate(void) const
{
  for (int i = 0; i < 2; i++)
    if (theMat[i] != 0 && !theMat[i]->ladrunoThreadSafeUpdate())
      return false;
  return true;
}

bool LadrunoIMKBeam2d::ladrunoThreadSafeForm(void) const
{
  for (int i = 0; i < 2; i++)
    if (theMat[i] != 0 && !theMat[i]->ladrunoThreadSafeForm())
      return false;
  return true;
}
