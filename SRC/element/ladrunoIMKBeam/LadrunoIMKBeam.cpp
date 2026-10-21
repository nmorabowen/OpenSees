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

static thread_local Matrix K(12, 12);
static thread_local Vector P(12);

// ---------------------------------------------------------------------------
// constructors / destructor
// ---------------------------------------------------------------------------
LadrunoIMKBeam::LadrunoIMKBeam()
  : Element(0, ELE_TAG_LadrunoIMKBeam),
    A(0.0), E(0.0), G(0.0), Jx(0.0), Iy(0.0), Iz(0.0), rho(0.0), cMass(0),
    GAvy(0.0), GAvz(0.0), P0(0.0), pmAlpha(0.5),
    dI(0.0), dJ(0.0), useImplex(false), implexErr(0.0),
    dtime_n_commit(0.0),
    q(6), kb(6, 6), Q(12), q0Dirty(false), committedAtCurrentV(false),
    connectedExternalNodes(2), theCoordTransf(0)
{
  for (int i = 0; i < 4; i++) {
    theMat[i] = 0;
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
    dThetaHprev[i] = 0.0;
    thetaHp[i] = 0.0;
    thetaHpCommit[i] = 0.0;
    dThetaHpPrev[i] = 0.0;
  }
  for (int i = 0; i < 2; i++) {
    theBiaxialMat[i] = 0;
    theShearMat[i] = 0;
    gammaS[i] = 0.0;
    gammaScommit[i] = 0.0;
    dGammaSprev[i] = 0.0;
    gammaSp[i] = 0.0;
    gammaSpCommit[i] = 0.0;
    dGammaSpPrev[i] = 0.0;
  }
  for (int i = 0; i < 5; i++) {
    q0[i] = 0.0;
    p0[i] = 0.0;
  }
  for (int i = 0; i < 6; i++)
    vCommit[i] = 0.0;
  theNodes[0] = 0;
  theNodes[1] = 0;
}

LadrunoIMKBeam::LadrunoIMKBeam(int tag, int Nd1, int Nd2,
                               double a, double e, double g, double jx,
                               double iy, double iz, CrdTransf &theTransf,
                               UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                               UniaxialMaterial *mYi, UniaxialMaterial *mYj,
                               double r, int cm,
                               double gavy, double gavz,
                               double p0_pm, double alpha_pm,
                               double di, double dj, bool implex,
                               UniaxialMaterial *mShearY, UniaxialMaterial *mShearZ,
                               NDMaterial *mBiI, NDMaterial *mBiJ)
  : Element(tag, ELE_TAG_LadrunoIMKBeam),
    A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), rho(r), cMass(cm),
    GAvy(gavy), GAvz(gavz), P0(p0_pm), pmAlpha(alpha_pm),
    dI(di), dJ(dj), useImplex(implex), implexErr(0.0),
    dtime_n_commit(0.0),
    q(6), kb(6, 6), Q(12), q0Dirty(false), committedAtCurrentV(false),
    connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  UniaxialMaterial *src[4] = { mZi, mZj, mYi, mYj };
  for (int i = 0; i < 4; i++) {
    if (src[i] != 0) {
      theMat[i] = src[i]->getCopy();
      if (theMat[i] == 0) {
        opserr << "LadrunoIMKBeam::LadrunoIMKBeam -- failed to get copy of material " << i << endln;
        exit(-1);
      }
      if (P0 > 0.0) {
        if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
          ch->setPM(P0, pmAlpha);
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

  NDMaterial *srcBi[2] = { mBiI, mBiJ };
  UniaxialMaterial *srcShear[2] = { mShearY, mShearZ };
  for (int i = 0; i < 2; i++) {
    if (srcBi[i] != 0) {
      theBiaxialMat[i] = srcBi[i]->getCopy();
      if (theBiaxialMat[i] == 0) {
        opserr << "LadrunoIMKBeam::LadrunoIMKBeam -- failed to get copy of biaxial material " << i << endln;
        exit(-1);
      }
      if (P0 > 0.0) {
        if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
          cb->setPM(P0, pmAlpha);
      }
    } else {
      theBiaxialMat[i] = 0;
    }
    if (srcShear[i] != 0) {
      theShearMat[i] = srcShear[i]->getCopy();
      if (theShearMat[i] == 0) {
        opserr << "LadrunoIMKBeam::LadrunoIMKBeam -- failed to get copy of shear material " << i << endln;
        exit(-1);
      }
    } else {
      theShearMat[i] = 0;
    }
    gammaS[i] = 0.0;
    gammaScommit[i] = 0.0;
    dGammaSprev[i] = 0.0;
    gammaSp[i] = 0.0;
    gammaSpCommit[i] = 0.0;
    dGammaSpPrev[i] = 0.0;
  }

  for (int i = 0; i < 5; i++) {
    q0[i] = 0.0;
    p0[i] = 0.0;
  }
  for (int i = 0; i < 6; i++)
    vCommit[i] = 0.0;

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
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0)
      delete theBiaxialMat[i];
    if (theShearMat[i] != 0)
      delete theShearMat[i];
  }
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

  this->update();
}

// ---------------------------------------------------------------------------
// state cycle
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::commitState(void)
{
  // Base-class commit first: it refreshes Kc, the committed stiffness behind
  // betaKc Rayleigh damping. Skipping it left Kc frozen at the tangent captured
  // when `rayleigh` ran, so betaKc silently behaved like betaK0 once hinges
  // yielded or the geometry rotated (WP-118). Same order as ElasticBeam3d.  // Ladruno
  int ok = this->Element::commitState();
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->commitState();
    dThetaHprev[i] = thetaH[i] - thetaHcommit[i];
    thetaHcommit[i] = thetaH[i];
    dThetaHpPrev[i] = thetaHp[i] - thetaHpCommit[i];
    thetaHpCommit[i] = thetaHp[i];
  }
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0)
      ok += theBiaxialMat[i]->commitState();
    if (theShearMat[i] != 0)
      ok += theShearMat[i]->commitState();
    dGammaSprev[i] = gammaS[i] - gammaScommit[i];
    gammaScommit[i] = gammaS[i];
    dGammaSpPrev[i] = gammaSp[i] - gammaSpCommit[i];
    gammaSpCommit[i] = gammaSp[i];
  }
  dtime_n_commit = ops_Dt;
  if (theCoordTransf != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    const Vector &vc = theCoordTransf->getBasicTrialDisp();
    for (int i = 0; i < 6; i++)
      vCommit[i] = vc(i);
  }
  ok += theCoordTransf->commitState();
  return ok;
}

int LadrunoIMKBeam::revertToLastCommit(void)
{
  int ok = 0;
  committedAtCurrentV = false;
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToLastCommit();
    thetaH[i] = thetaHcommit[i];
    thetaHp[i] = thetaHpCommit[i];
  }
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0)
      ok += theBiaxialMat[i]->revertToLastCommit();
    if (theShearMat[i] != 0)
      ok += theShearMat[i]->revertToLastCommit();
    gammaS[i] = gammaScommit[i];
    gammaSp[i] = gammaSpCommit[i];
  }
  ok += theCoordTransf->revertToLastCommit();
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0)
    ok += this->update();
  return ok;
}

int LadrunoIMKBeam::revertToStart(void)
{
  int ok = 0;
  committedAtCurrentV = false;
  for (int i = 0; i < 6; i++)
    vCommit[i] = 0.0;
  for (int i = 0; i < 4; i++) {
    if (theMat[i] != 0)
      ok += theMat[i]->revertToStart();
    thetaH[i] = 0.0;
    thetaHcommit[i] = 0.0;
    dThetaHprev[i] = 0.0;
    thetaHp[i] = 0.0;
    thetaHpCommit[i] = 0.0;
    dThetaHpPrev[i] = 0.0;
  }
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0)
      ok += theBiaxialMat[i]->revertToStart();
    if (theShearMat[i] != 0)
      ok += theShearMat[i]->revertToStart();
    gammaS[i] = 0.0;
    gammaScommit[i] = 0.0;
    dGammaSprev[i] = 0.0;
    gammaSp[i] = 0.0;
    gammaSpCommit[i] = 0.0;
    dGammaSpPrev[i] = 0.0;
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
int LadrunoIMKBeam::update(void)
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();
  if (committedAtCurrentV) {
    if (!q0Dirty &&
        v(0) == vCommit[0] && v(1) == vCommit[1] && v(2) == vCommit[2] &&
        v(3) == vCommit[3] && v(4) == vCommit[4] && v(5) == vCommit[5]) {
      return ok;
    }
    committedAtCurrentV = false;
  }
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;
  const double EAoverL = E * A * oneOverL;
  const double GJoverL = G * Jx * oneOverL;

  kb.Zero();
  q.Zero();

  // axial + torsion: linear elastic, no hinge
  kb(0, 0) = EAoverL;  q(0) = EAoverL * v(0) + q0[0];
  kb(5, 5) = GJoverL;  q(5) = GJoverL * v(5);

  const double fPM = ladrunoIMKScalePM(q(0), P0, pmAlpha);
  const double dLogPM = (fPM > 0.0) ? (ladrunoIMKScalePMDeriv(q(0), P0, pmAlpha) / fPM) : 0.0;
  for (int i = 0; i < 4; i++)
    if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
      ch->setAxialForce(q(0));
  for (int i = 0; i < 2; i++) {
    if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
      cb->setAxialForce(q(0));
    if (LadrunoCohesiveHinge *chs = dynamic_cast<LadrunoCohesiveHinge *>(theShearMat[i]))
      chs->setAxialForce(q(0));
  }

  double dtScale = 1.0;
  if (dtime_n_commit > 1.0e-16 && ops_Dt > 0.0) {
    dtScale = ops_Dt / dtime_n_commit;
    if (dtScale > 2.0) dtScale = 2.0;
  }

  if (theBiaxialMat[0] != 0 || theBiaxialMat[1] != 0) {
    const double v4[4] = { v(1), v(2), v(3), v(4) };
    const double q0_4[4] = { q0[1], q0[2], q0[3], q0[4] };
    double q4[4] = { 0.0, 0.0, 0.0, 0.0 };
    double dq4_dN[4] = { 0.0, 0.0, 0.0, 0.0 };
    double k4[4][4];
    ladrunoIMKSolveBiaxial(theBiaxialMat[0], theBiaxialMat[1],
                           v4, L, E * Iz, E * Iy,
                           thetaH, q4, k4, q0_4,
                           GAvy, GAvz, fPM, dI, dJ, useImplex,
                           thetaHpCommit, dThetaHpPrev, &implexErr,
                           dtScale, thetaHp,
                           theShearMat[0], theShearMat[1],
                           &gammaS[0], &gammaS[1],
                           gammaSpCommit[0], gammaSpCommit[1],
                           dGammaSpPrev[0], dGammaSpPrev[1],
                           &gammaSp[0], &gammaSp[1],
                           dLogPM, dq4_dN);
    for (int i = 0; i < 4; i++) {
      q(1 + i) = q4[i];
      kb(1 + i, 0) = EAoverL * dq4_dN[i];
      for (int j = 0; j < 4; j++)
        kb(1 + i, 1 + j) = k4[i][j];
    }
    q0Dirty = false;
    return ok;
  }

  double k2[2][2];

  // strong-axis bending: basic dofs 1,2 ; materials 0,1 ; EI = E*Iz
  double qi, qj;
  double errZ = 0.0, errY = 0.0;
  double dqZ_dN[2] = { 0.0, 0.0 }, dqY_dN[2] = { 0.0, 0.0 };
  ladrunoIMKSolveAxis(theMat[0], theMat[1], v(1), v(2), L, E * Iz,
                      thetaH[0], thetaH[1], qi, qj, k2, q0[1], q0[2],
                      GAvy, fPM, dI, dJ, useImplex,
                      thetaHpCommit[0], thetaHpCommit[1],
                      dThetaHpPrev[0], dThetaHpPrev[1], &errZ,
                      dtScale, &thetaHp[0], &thetaHp[1],
                      theShearMat[0], &gammaS[0],
                      gammaSpCommit[0], dGammaSpPrev[0], &gammaSp[0],
                      dLogPM, dqZ_dN);
  q(1) = qi;  q(2) = qj;
  kb(1, 0) = EAoverL * dqZ_dN[0];
  kb(2, 0) = EAoverL * dqZ_dN[1];
  kb(1, 1) = k2[0][0];  kb(1, 2) = k2[0][1];
  kb(2, 1) = k2[1][0];  kb(2, 2) = k2[1][1];

  // weak-axis bending: basic dofs 3,4 ; materials 2,3 ; EI = E*Iy
  ladrunoIMKSolveAxis(theMat[2], theMat[3], v(3), v(4), L, E * Iy,
                      thetaH[2], thetaH[3], qi, qj, k2, q0[3], q0[4],
                      GAvz, fPM, dI, dJ, useImplex,
                      thetaHpCommit[2], thetaHpCommit[3],
                      dThetaHpPrev[2], dThetaHpPrev[3], &errY,
                      dtScale, &thetaHp[2], &thetaHp[3],
                      theShearMat[1], &gammaS[1],
                      gammaSpCommit[1], dGammaSpPrev[1], &gammaSp[1],
                      dLogPM, dqY_dN);
  implexErr = errZ + errY;
  q(3) = qi;  q(4) = qj;
  kb(3, 0) = EAoverL * dqY_dN[0];
  kb(4, 0) = EAoverL * dqY_dN[1];
  kb(3, 3) = k2[0][0];  kb(3, 4) = k2[0][1];
  kb(4, 3) = k2[1][0];  kb(4, 4) = k2[1][1];
  q0Dirty = false;

  return ok;
}

const Matrix &LadrunoIMKBeam::getTangentStiff(void)
{
  if (q0Dirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &LadrunoIMKBeam::getInitialStiff(void)
{
  const double L = theCoordTransf->getInitialLength();
  const double oneOverL = 1.0 / L;

  static thread_local Matrix kbi(6, 6);
  kbi.Zero();
  kbi(0, 0) = E * A * oneOverL;
  kbi(5, 5) = G * Jx * oneOverL;

  if (theBiaxialMat[0] != 0 || theBiaxialMat[1] != 0) {
    double k4[4][4];
    ladrunoIMKInitBiaxial(theBiaxialMat[0], theBiaxialMat[1],
                          L, E * Iz, E * Iy, k4,
                          GAvy, GAvz, 1.0, dI, dJ,
                          theShearMat[0], theShearMat[1]);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        kbi(1 + i, 1 + j) = k4[i][j];
    return theCoordTransf->getInitialGlobalStiffMatrix(kbi);
  }

  double k2[2][2];
  ladrunoIMKInitBlock(theMat[0], theMat[1], L, E * Iz, k2, GAvy, 1.0, dI, dJ, theShearMat[0]);
  kbi(1, 1) = k2[0][0];  kbi(1, 2) = k2[0][1];
  kbi(2, 1) = k2[1][0];  kbi(2, 2) = k2[1][1];
  ladrunoIMKInitBlock(theMat[2], theMat[3], L, E * Iy, k2, GAvz, 1.0, dI, dJ, theShearMat[1]);
  kbi(3, 3) = k2[0][0];  kbi(3, 4) = k2[0][1];
  kbi(4, 3) = k2[1][0];  kbi(4, 4) = k2[1][1];

  return theCoordTransf->getInitialGlobalStiffMatrix(kbi);
}

const Matrix &LadrunoIMKBeam::getMass(void)
{
  K.Zero();
  if (rho <= 0.0)
    return K;

  double L = theCoordTransf->getInitialLength();
  if (cMass == 0) {
    double m = 0.5 * rho * L;
    K(0, 0) = K(1, 1) = K(2, 2) = m;
    K(6, 6) = K(7, 7) = K(8, 8) = m;
  } else {
    static thread_local Matrix ml(12, 12);
    ml.Zero();
    double m = rho * L / 420.0;
    ml(0, 0) = ml(6, 6) = m * 140.0;
    ml(0, 6) = ml(6, 0) = m * 70.0;

    ml(2, 2) = ml(8, 8) = m * 156.0;
    ml(2, 8) = ml(8, 2) = m * 54.0;
    ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
    ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
    ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
    ml(8, 10) = ml(10, 8) = -ml(2, 4);
    ml(2, 10) = ml(10, 2) = m * 13.0 * L;
    ml(4, 8) = ml(8, 4) = -ml(2, 10);

    ml(1, 1) = ml(7, 7) = m * 156.0;
    ml(1, 7) = ml(7, 1) = m * 54.0;
    ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
    ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
    ml(1, 5) = ml(5, 1) = m * 22.0 * L;
    ml(7, 11) = ml(11, 7) = -ml(1, 5);
    ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
    ml(5, 7) = ml(7, 5) = -ml(1, 11);

    K = theCoordTransf->getGlobalMatrixFromLocal(ml);
  }
  return K;
}

// ---------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------
void LadrunoIMKBeam::zeroLoad(void)
{
  Q.Zero();
  for (int i = 0; i < 5; i++) {
    if (q0[i] != 0.0)
      q0Dirty = true;
    q0[i] = 0.0;
    p0[i] = 0.0;
  }
}

int LadrunoIMKBeam::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();
  const double Lc = ((dI > 0.0 || dJ > 0.0) && (dI + dJ < L)) ? (L - dI - dJ) : L;

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0) * loadFactor;  // Transverse Y
    double wz = data(1) * loadFactor;  // Transverse Z
    double wx = data(2) * loadFactor;  // Axial (+ve from node I to J)

    double Vy = 0.5 * wy * L;
    double Mz = Vy * L / 6.0; // wy*L*L/12
    double Vz = 0.5 * wz * L;
    double My = Vz * L / 6.0; // wz*L*L/12
    double P  = wx * L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5 * P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;
    q0Dirty = true;
  }
  else if (type == LOAD_TAG_BeamUniformMoment) {
    double my = data(1) * loadFactor;  // About y
    double mz = data(2) * loadFactor;  // About z

    // Reactions in basic system
    p0[1] += mz;
    p0[2] -= mz;
    p0[3] -= my;
    p0[4] += my;
  }
  else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
    double wya = data(0) * loadFactor;  // Transverse Y at start
    double wza = data(1) * loadFactor;  // Transverse Z at start
    double waa = data(2) * loadFactor;  // Axial at start
    double a   = data(3) * L;
    double b   = data(4) * L;
    double wyb = data(5) * loadFactor;  // Transverse Y at end
    double wzb = data(6) * loadFactor;  // Transverse Z at end
    double wab = data(7) * loadFactor;  // Axial at end

    double ell = b - a;
    if (ell <= 0.0 || a < 0.0 || b > L)
      return 0;

    // Axial
    double P  = 0.5 * ell * (waa + wab);
    double PJ = (a * P + (ell * ell / 6.0) * (waa + 2.0 * wab)) / L;

    // Transverse Y (bending about z)
    double J0y = 0.5 * ell * (wya + wyb);
    double J1y = (ell * ell / 6.0) * (wya + 2.0 * wyb);
    double J2y = (ell * ell * ell / 12.0) * (wya + 3.0 * wyb);
    double J3y = (ell * ell * ell * ell / 20.0) * (wya + 4.0 * wyb);
    double I0y = J0y;
    double I1y = a * J0y + J1y;
    double I2y = a * a * J0y + 2.0 * a * J1y + J2y;
    double I3y = a * a * a * J0y + 3.0 * a * a * J1y + 3.0 * a * J2y + J3y;
    double Vy2 = I1y / L;
    double Vy1 = I0y - Vy2;
    double M1y = -I1y + (2.0 / L) * I2y - (1.0 / (L * L)) * I3y;
    double M2y = (1.0 / L) * I2y - (1.0 / (L * L)) * I3y;
    if (GAvy > 0.0 && Lc > 0.0) {
      const double Phi = 12.0 * E * Iz / (GAvy * Lc * Lc);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1y - 0.5 * Phi * M2y) * inv;
      const double m2 = (-0.5 * Phi * M1y + (1.0 + 0.5 * Phi) * M2y) * inv;
      M1y = m1;
      M2y = m2;
    }

    // Transverse Z (bending about y)
    double J0z = 0.5 * ell * (wza + wzb);
    double J1z = (ell * ell / 6.0) * (wza + 2.0 * wzb);
    double J2z = (ell * ell * ell / 12.0) * (wza + 3.0 * wzb);
    double J3z = (ell * ell * ell * ell / 20.0) * (wza + 4.0 * wzb);
    double I0z = J0z;
    double I1z = a * J0z + J1z;
    double I2z = a * a * J0z + 2.0 * a * J1z + J2z;
    double I3z = a * a * a * J0z + 3.0 * a * a * J1z + 3.0 * a * J2z + J3z;
    double Vz2 = I1z / L;
    double Vz1 = I0z - Vz2;
    double M1z = -I1z + (2.0 / L) * I2z - (1.0 / (L * L)) * I3z;
    double M2z = (1.0 / L) * I2z - (1.0 / (L * L)) * I3z;
    if (GAvz > 0.0 && Lc > 0.0) {
      const double Phi = 12.0 * E * Iy / (GAvz * Lc * Lc);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1z - 0.5 * Phi * M2z) * inv;
      const double m2 = (-0.5 * Phi * M1z + (1.0 + 0.5 * Phi) * M2z) * inv;
      M1z = m1;
      M2z = m2;
    }

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;

    // Fixed end forces in basic system
    q0[0] -= PJ;
    q0[1] += M1y;
    q0[2] += M2y;
    q0[3] -= M1z;
    q0[4] -= M2z;
    q0Dirty = true;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py     = data(0) * loadFactor;
    double Pz     = data(1) * loadFactor;
    double N      = data(2) * loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL * L;
    double b = L - a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = Py * (1.0 - aOverL);
    double V2 = Py * aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz * (1.0 - aOverL);
    V2 = Pz * aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0 / (L * L);
    double a2 = a * a;
    double b2 = b * b;

    // Fixed end forces in basic system
    q0[0] -= N * aOverL;
    double M1 = -a * b2 * Py * L2;
    double M2 =  a2 * b * Py * L2;
    if (GAvy > 0.0 && Lc > 0.0) {
      const double Phi = 12.0 * E * Iz / (GAvy * Lc * Lc);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1 - 0.5 * Phi * M2) * inv;
      const double m2 = (-0.5 * Phi * M1 + (1.0 + 0.5 * Phi) * M2) * inv;
      M1 = m1;
      M2 = m2;
    }
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 =  a2 * b * Pz * L2;
    if (GAvz > 0.0 && Lc > 0.0) {
      const double Phi = 12.0 * E * Iy / (GAvz * Lc * Lc);
      const double inv = 1.0 / (1.0 + Phi);
      const double m1 = ((1.0 + 0.5 * Phi) * M1 - 0.5 * Phi * M2) * inv;
      const double m2 = (-0.5 * Phi * M1 + (1.0 + 0.5 * Phi) * M2) * inv;
      M1 = m1;
      M2 = m2;
    }
    q0[3] -= M1;
    q0[4] -= M2;
    q0Dirty = true;
  }
  else {
    opserr << "LadrunoIMKBeam::addLoad -- load type unknown for element with tag: "
           << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int LadrunoIMKBeam::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "LadrunoIMKBeam::addInertiaLoadToUnbalance -- incompatible node DOF size\n";
    return -1;
  }

  if (cMass == 0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);
    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);
  } else {
    static thread_local Vector Raccel(12);
    for (int i = 0; i < 6; i++) {
      Raccel(i)     = Raccel1(i);
      Raccel(i + 6) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

// ---------------------------------------------------------------------------
// resisting force
// ---------------------------------------------------------------------------
const Vector &LadrunoIMKBeam::getResistingForce(void)
{
  if (q0Dirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();

  Vector p0Vec(p0, 5);

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // subtract element nodal loads (inertia unbalance accumulated in Q)
  P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector &LadrunoIMKBeam::getResistingForceIncInertia(void)
{
  // SNAPSHOT into a function-local buffer BEFORE getRayleighDampingForces():
  // betaK Rayleigh re-enters getTangentStiff(); returning `res` keeps the result
  // independent of anything that path writes, and of the static P that
  // getResistingForce() returns (Element::getResponse 444444 subtracts both).
  // Order is unchanged -- ((f - Q) + R) + m*a -- so results are bit-identical
  // (WP-115; LEDGER_quirks "MUST snapshot the shared static `resid`").  // Ladruno
  static thread_local Vector res(12);
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
    res(2) += m * accel1(2);
    res(6) += m * accel2(0);
    res(7) += m * accel2(1);
    res(8) += m * accel2(2);
  } else {
    static thread_local Vector accel(12);
    for (int i = 0; i < 6; i++) {
      accel(i)     = accel1(i);
      accel(i + 6) = accel2(i);
    }
    res.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }

  return res;
}

// ---------------------------------------------------------------------------
// parallel
// ---------------------------------------------------------------------------
int LadrunoIMKBeam::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();

  // integer data: tag, nodes, transf class+db, 4x (mat class+db), cMass,
  //               2x (shearMat class+db), 2x (biaxialMat class+db)
  static thread_local ID idData(22);
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
  idData(13) = cMass;
  for (int i = 0; i < 2; i++) {
    if (theShearMat[i] != 0) {
      idData(14 + 2 * i) = theShearMat[i]->getClassTag();
      int sDb = theShearMat[i]->getDbTag();
      if (sDb == 0) {
        sDb = theChannel.getDbTag();
        if (sDb != 0)
          theShearMat[i]->setDbTag(sDb);
      }
      idData(15 + 2 * i) = sDb;
    } else {
      idData(14 + 2 * i) = 0;
      idData(15 + 2 * i) = 0;
    }
    if (theBiaxialMat[i] != 0) {
      idData(18 + 2 * i) = theBiaxialMat[i]->getClassTag();
      int bDb = theBiaxialMat[i]->getDbTag();
      if (bDb == 0) {
        bDb = theChannel.getDbTag();
        if (bDb != 0)
          theBiaxialMat[i]->setDbTag(bDb);
      }
      idData(19 + 2 * i) = bDb;
    } else {
      idData(18 + 2 * i) = 0;
      idData(19 + 2 * i) = 0;
    }
  }

  res += theChannel.sendID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::sendSelf -- failed to send ID\n";
    return res;
  }

  static thread_local Vector dData(102);
  dData(0) = A;  dData(1) = E;  dData(2) = G;  dData(3) = Jx;
  dData(4) = Iy; dData(5) = Iz; dData(6) = rho;
  dData(7) = thetaHcommit[0];  dData(8) = thetaHcommit[1];
  dData(9) = thetaHcommit[2];  dData(10) = thetaHcommit[3];
  dData(11) = alphaM;  dData(12) = betaK;  dData(13) = betaK0;  dData(14) = betaKc;
  dData(15) = GAvy;    dData(16) = GAvz;   dData(17) = P0;      dData(18) = pmAlpha;
  dData(19) = dI;      dData(20) = dJ;     dData(21) = useImplex ? 1.0 : 0.0;
  for (int i = 0; i < 4; i++) {
    dData(22 + i) = thetaHpCommit[i];
    dData(26 + i) = dThetaHpPrev[i];
    dData(30 + i) = dThetaHprev[i];
  }
  dData(34) = dtime_n_commit;
  dData(35) = implexErr;
  for (int i = 0; i < 6; i++)
    dData(36 + i) = q(i);
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      dData(42 + 6 * i + j) = kb(i, j);
  for (int i = 0; i < 6; i++)
    dData(78 + i) = vCommit[i];
  for (int i = 0; i < 5; i++) {
    dData(84 + i) = q0[i];
    dData(89 + i) = p0[i];
  }
  for (int i = 0; i < 2; i++) {
    dData(94 + i) = gammaScommit[i];
    dData(96 + i) = dGammaSprev[i];
    dData(98 + i) = gammaSpCommit[i];
    dData(100 + i) = dGammaSpPrev[i];
  }

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
  for (int i = 0; i < 2; i++) {
    if (theShearMat[i] != 0) {
      res += theShearMat[i]->sendSelf(cTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoIMKBeam::sendSelf -- failed to send shear material " << i << endln;
        return res;
      }
    }
    if (theBiaxialMat[i] != 0) {
      res += theBiaxialMat[i]->sendSelf(cTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoIMKBeam::sendSelf -- failed to send biaxial material " << i << endln;
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

  static thread_local ID idData(22);
  res += theChannel.recvID(dbTag, cTag, idData);
  if (res < 0) {
    opserr << "LadrunoIMKBeam::recvSelf -- failed to recv ID\n";
    return res;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  cMass = idData(13);

  static thread_local Vector dData(102);
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
  alphaM = dData(11);  betaK = dData(12);  betaK0 = dData(13);  betaKc = dData(14);
  GAvy = dData(15);    GAvz = dData(16);   P0 = dData(17);      pmAlpha = dData(18);
  dI = dData(19);      dJ = dData(20);     useImplex = (dData(21) != 0.0);
  for (int i = 0; i < 4; i++) {
    thetaHpCommit[i] = dData(22 + i);
    thetaHp[i] = thetaHpCommit[i];
    dThetaHpPrev[i] = dData(26 + i);
    dThetaHprev[i] = dData(30 + i);
  }
  dtime_n_commit = dData(34);
  implexErr = dData(35);
  for (int i = 0; i < 6; i++)
    q(i) = dData(36 + i);
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      kb(i, j) = dData(42 + 6 * i + j);
  for (int i = 0; i < 6; i++)
    vCommit[i] = dData(78 + i);
  for (int i = 0; i < 5; i++) {
    q0[i] = dData(84 + i);
    p0[i] = dData(89 + i);
  }
  for (int i = 0; i < 2; i++) {
    gammaScommit[i] = dData(94 + i);
    gammaS[i] = gammaScommit[i];
    dGammaSprev[i] = dData(96 + i);
    gammaSpCommit[i] = dData(98 + i);
    gammaSp[i] = gammaSpCommit[i];
    dGammaSpPrev[i] = dData(100 + i);
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
    if (P0 > 0.0) {
      if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
        ch->setPM(P0, pmAlpha);
    }
  }

  for (int i = 0; i < 2; i++) {
    int sClass = idData(14 + 2 * i);
    int sDb = idData(15 + 2 * i);
    if (sClass == 0) {
      if (theShearMat[i] != 0) {
        delete theShearMat[i];
        theShearMat[i] = 0;
      }
    } else {
      if (theShearMat[i] == 0 || theShearMat[i]->getClassTag() != sClass) {
        if (theShearMat[i] != 0)
          delete theShearMat[i];
        theShearMat[i] = theBroker.getNewUniaxialMaterial(sClass);
        if (theShearMat[i] == 0) {
          opserr << "LadrunoIMKBeam::recvSelf -- could not create shear material " << i << endln;
          return -1;
        }
      }
      theShearMat[i]->setDbTag(sDb);
      res += theShearMat[i]->recvSelf(cTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoIMKBeam::recvSelf -- failed to recv shear material " << i << endln;
        return res;
      }
    }

    int bClass = idData(18 + 2 * i);
    int bDb = idData(19 + 2 * i);
    if (bClass == 0) {
      if (theBiaxialMat[i] != 0) {
        delete theBiaxialMat[i];
        theBiaxialMat[i] = 0;
      }
    } else {
      if (theBiaxialMat[i] == 0 || theBiaxialMat[i]->getClassTag() != bClass) {
        if (theBiaxialMat[i] != 0)
          delete theBiaxialMat[i];
        theBiaxialMat[i] = theBroker.getNewNDMaterial(bClass);
        if (theBiaxialMat[i] == 0) {
          opserr << "LadrunoIMKBeam::recvSelf -- could not create biaxial material " << i << endln;
          return -1;
        }
      }
      theBiaxialMat[i]->setDbTag(bDb);
      res += theBiaxialMat[i]->recvSelf(cTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoIMKBeam::recvSelf -- failed to recv biaxial material " << i << endln;
        return res;
      }
      if (P0 > 0.0) {
        if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
          cb->setPM(P0, pmAlpha);
      }
    }
  }

  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0 && theCoordTransf != 0)
    theCoordTransf->update();

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
    << " Jx: " << Jx << " Iy: " << Iy << " Iz: " << Iz
    << " rho: " << rho << " cMass: " << cMass;
  if (GAvy > 0.0 || GAvz > 0.0)
    s << " GAvy: " << GAvy << " GAvz: " << GAvz;
  if (P0 > 0.0)
    s << " P0: " << P0 << " pmAlpha: " << pmAlpha;
  s << endln;
  if (theBiaxialMat[0] != 0 || theBiaxialMat[1] != 0) {
    s << "  biaxialHinge@i: " << (theBiaxialMat[0] ? theBiaxialMat[0]->getTag() : 0)
      << " biaxialHinge@j: " << (theBiaxialMat[1] ? theBiaxialMat[1]->getTag() : 0) << endln;
  } else {
    const char *lbl[4] = { "Mz@i", "Mz@j", "My@i", "My@j" };
    for (int i = 0; i < 4; i++)
      s << "  hinge " << lbl[i] << ": "
        << (theMat[i] != 0 ? "material " : "elastic")
        << (theMat[i] != 0 ? theMat[i]->getTag() : 0) << endln;
  }
  if (theShearMat[0] != 0 || theShearMat[1] != 0) {
    s << "  shearHinge Vy: " << (theShearMat[0] ? theShearMat[0]->getTag() : 0)
      << " Vz: " << (theShearMat[1] ? theShearMat[1]->getTag() : 0) << endln;
  }
}

// ---------------------------------------------------------------------------
// recorders
// ---------------------------------------------------------------------------
Response *LadrunoIMKBeam::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  Response *theResponse = 0;

  if (argc < 1) return 0;

  s.tag("ElementOutput");
  s.attr("eleType", "LadrunoIMKBeam");
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
    s.tag("ResponseType", "Pz_1");
    s.tag("ResponseType", "Mx_1");
    s.tag("ResponseType", "My_1");
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "Px_2");
    s.tag("ResponseType", "Py_2");
    s.tag("ResponseType", "Pz_2");
    s.tag("ResponseType", "Mx_2");
    s.tag("ResponseType", "My_2");
    s.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 1, P);

  } else if (LadrunoResp::isLocalForce(argv[0])) {
    s.tag("ResponseType", "N_1");
    s.tag("ResponseType", "Vy_1");
    s.tag("ResponseType", "Vz_1");
    s.tag("ResponseType", "T_1");
    s.tag("ResponseType", "My_1");
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "N_2");
    s.tag("ResponseType", "Vy_2");
    s.tag("ResponseType", "Vz_2");
    s.tag("ResponseType", "T_2");
    s.tag("ResponseType", "My_2");
    s.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 6, Vector(12));

  } else if (LadrunoResp::isBasicForce(argv[0])) {
    s.tag("ResponseType", "N");
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "Mz_2");
    s.tag("ResponseType", "My_1");
    s.tag("ResponseType", "My_2");
    s.tag("ResponseType", "T");
    theResponse = new ElementResponse(this, 2, Vector(6));

  } else if (strcmp(argv[0], "basicStiffness") == 0) {
    s.tag("ResponseType", "N");
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "Mz_2");
    s.tag("ResponseType", "My_1");
    s.tag("ResponseType", "My_2");
    s.tag("ResponseType", "T");
    theResponse = new ElementResponse(this, 19, Matrix(6, 6));

  } else if (LadrunoResp::isBasicDef(argv[0])) {
    s.tag("ResponseType", "eps");
    s.tag("ResponseType", "thetaZ_1");
    s.tag("ResponseType", "thetaZ_2");
    s.tag("ResponseType", "thetaY_1");
    s.tag("ResponseType", "thetaY_2");
    s.tag("ResponseType", "thetaX");
    theResponse = new ElementResponse(this, 3, Vector(6));

  } else if (LadrunoResp::isPlasticDef(argv[0])) {
    s.tag("ResponseType", "thetaZP_1");
    s.tag("ResponseType", "thetaZP_2");
    s.tag("ResponseType", "thetaYP_1");
    s.tag("ResponseType", "thetaYP_2");
    theResponse = new ElementResponse(this, 4, Vector(4));

  } else if (LadrunoResp::is(argv[0], "hingeMoment")) {
    s.tag("ResponseType", "Mz_1");
    s.tag("ResponseType", "Mz_2");
    s.tag("ResponseType", "My_1");
    s.tag("ResponseType", "My_2");
    theResponse = new ElementResponse(this, 5, Vector(4));

  } else if (LadrunoResp::is(argv[0], "hingeSummary")) {
    s.tag("ResponseType", "thetaZ_i");
    s.tag("ResponseType", "Mz_i");
    s.tag("ResponseType", "thetaZ_j");
    s.tag("ResponseType", "Mz_j");
    s.tag("ResponseType", "thetaY_i");
    s.tag("ResponseType", "My_i");
    s.tag("ResponseType", "thetaY_j");
    s.tag("ResponseType", "My_j");
    theResponse = new ElementResponse(this, 8, Vector(8));

  } else if (LadrunoResp::is(argv[0], "shearDeformation") ||
             LadrunoResp::is(argv[0], "shearStrain") ||
             LadrunoResp::is(argv[0], "shearDistortion")) {
    s.tag("ResponseType", "gamma_y");
    s.tag("ResponseType", "gamma_z");
    theResponse = new ElementResponse(this, 14, Vector(2));

  } else if (LadrunoResp::is(argv[0], "shearForce")) {
    s.tag("ResponseType", "Vy");
    s.tag("ResponseType", "Vz");
    theResponse = new ElementResponse(this, 15, Vector(2));

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
      theResponse = new ElementResponse(this, (secNum == 1) ? 71 : 72, Vector(4));
    }
  } else if (strcmp(argv[0], "shearHinge") == 0) {
    if (argc > 1) {
      if (theShearMat[0] != 0)
        theResponse = theShearMat[0]->setResponse(&argv[1], argc - 1, s);
    } else {
      s.tag("ResponseType", "gamma_y");
      s.tag("ResponseType", "Vy");
      s.tag("ResponseType", "gamma_z");
      s.tag("ResponseType", "Vz");
      theResponse = new ElementResponse(this, 16, Vector(4));
    }
  } else if (strcmp(argv[0], "shearHingeY") == 0) {
    if (argc > 1) {
      if (theShearMat[0] != 0)
        theResponse = theShearMat[0]->setResponse(&argv[1], argc - 1, s);
    } else {
      s.tag("ResponseType", "gamma_y");
      s.tag("ResponseType", "Vy");
      theResponse = new ElementResponse(this, 17, Vector(2));
    }
  } else if (strcmp(argv[0], "shearHingeZ") == 0) {
    if (argc > 1) {
      if (theShearMat[1] != 0)
        theResponse = theShearMat[1]->setResponse(&argv[1], argc - 1, s);
    } else {
      s.tag("ResponseType", "gamma_z");
      s.tag("ResponseType", "Vz");
      theResponse = new ElementResponse(this, 18, Vector(2));
    }
  } else if (strcmp(argv[0], "hingeBiaxial") == 0) {
    if (argc > 1) {
      if (theBiaxialMat[0] != 0)
        theResponse = theBiaxialMat[0]->setResponse(&argv[1], argc - 1, s);
    } else {
      theResponse = new ElementResponse(this, 8, Vector(8));
    }
  } else if (strcmp(argv[0], "hingeBiaxialI") == 0) {
    if (argc > 1) {
      if (theBiaxialMat[0] != 0)
        theResponse = theBiaxialMat[0]->setResponse(&argv[1], argc - 1, s);
    } else {
      theResponse = new ElementResponse(this, 20, Vector(4));
    }
  } else if (strcmp(argv[0], "hingeBiaxialJ") == 0) {
    if (argc > 1) {
      if (theBiaxialMat[1] != 0)
        theResponse = theBiaxialMat[1]->setResponse(&argv[1], argc - 1, s);
    } else {
      theResponse = new ElementResponse(this, 21, Vector(4));
    }
  } else if ((LadrunoResp::is(argv[0], "material") ||
              strcmp(argv[0], "hinge") == 0) && argc > 2) {
    if ((strcmp(argv[1], "shear") == 0 || strcmp(argv[1], "shearY") == 0) && theShearMat[0] != 0) {
      theResponse = theShearMat[0]->setResponse(&argv[2], argc - 2, s);
    } else if (strcmp(argv[1], "shearZ") == 0 && theShearMat[1] != 0) {
      theResponse = theShearMat[1]->setResponse(&argv[2], argc - 2, s);
    } else if ((strcmp(argv[1], "biaxialI") == 0 || strcmp(argv[1], "i") == 0 || strcmp(argv[1], "0") == 0)
               && theBiaxialMat[0] != 0) {
      theResponse = theBiaxialMat[0]->setResponse(&argv[2], argc - 2, s);
    } else if ((strcmp(argv[1], "biaxialJ") == 0 || strcmp(argv[1], "j") == 0 || strcmp(argv[1], "1") == 0)
               && theBiaxialMat[1] != 0) {
      theResponse = theBiaxialMat[1]->setResponse(&argv[2], argc - 2, s);
    } else {
      int idx = atoi(argv[1]);
      if (idx >= 0 && idx < 4 && theMat[idx] != 0)
        theResponse = theMat[idx]->setResponse(&argv[2], argc - 2, s);
    }
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
    static thread_local Vector th(4);
    for (int i = 0; i < 4; i++)
      th(i) = thetaH[i];
    return info.setVector(th);
  }
  case 5: {
    static thread_local Vector m(4);
    m(0) = q(1);  m(1) = q(2);  m(2) = q(3);  m(3) = q(4);
    return info.setVector(m);
  }
  case 6: {
    // Ladruno — local-frame end forces from the basic
    // q = [N, Mz_i, Mz_j, My_i, My_j, T], same mapping as DispBeamColumn3d
    // including p0 reactions.
    static thread_local Vector Pl(12);
    const double L = theCoordTransf->getInitialLength();
    const double oneOverL = (L > 0.0) ? 1.0 / L : 0.0;

    Pl(0)  = -q(0) + p0[0];   Pl(6)  = q(0);          // axial
    Pl(3)  = -q(5);           Pl(9)  = q(5);          // torsion

    double V = (q(1) + q(2)) * oneOverL;              // Mz -> shear along y
    Pl(5)  = q(1);            Pl(11) = q(2);
    Pl(1)  =  V + p0[1];      Pl(7)  = -V + p0[2];

    V = (q(3) + q(4)) * oneOverL;                     // My -> shear along z
    Pl(4)  = q(3);            Pl(10) = q(4);
    Pl(2)  = -V + p0[3];      Pl(8)  =  V + p0[4];

    return info.setVector(Pl);
  }
  case 8: {
    static thread_local Vector hs(8);
    hs(0) = thetaH[0];
    hs(1) = q(1);
    hs(2) = thetaH[1];
    hs(3) = q(2);
    hs(4) = thetaH[2];
    hs(5) = q(3);
    hs(6) = thetaH[3];
    hs(7) = q(4);
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
      const double vzic = T00 * v(1) + T01 * v(2);
      const double vzjc = T10 * v(1) + T11 * v(2);
      const double vyic = T00 * v(3) + T01 * v(4);
      const double vyjc = T10 * v(3) + T11 * v(4);
      const double ea = E * A / L;
      const double gj = G * Jx / L;
      double az, bz, F11z, F12z;
      double ay, by, F11y, F12y;
      ladrunoIMKCoreCoeffs(Lc, E * Iz, GAvy, az, bz, F11z, F12z);
      ladrunoIMKCoreCoeffs(Lc, E * Iy, GAvz, ay, by, F11y, F12y);
      double gamRotY = 0.0, gamRotZ = 0.0;
      if (theShearMat[0] != 0 && Lc > 0.0) {
        const double Ksy0 = theShearMat[0]->getInitialTangent();
        const double Fsy0 = (Ksy0 > 0.0) ? (1.0 / (Lc * Lc * Ksy0)) : 0.0;
        const double q0zic = hasEZ ? (((L - dI) * q0[1] - dI * q0[2]) / L) : q0[1];
        const double q0zjc = hasEZ ? ((-dJ * q0[1] + (L - dJ) * q0[2]) / L) : q0[2];
        const double gamRot0Y = -Fsy0 * (q0zic + q0zjc);
        gamRotY = (gammaS[0] / Lc) - gamRot0Y;
      }
      if (theShearMat[1] != 0 && Lc > 0.0) {
        const double Ksz0 = theShearMat[1]->getInitialTangent();
        const double Fsz0 = (Ksz0 > 0.0) ? (1.0 / (Lc * Lc * Ksz0)) : 0.0;
        const double q0yic = hasEZ ? (((L - dI) * q0[3] - dI * q0[4]) / L) : q0[3];
        const double q0yjc = hasEZ ? ((-dJ * q0[3] + (L - dJ) * q0[4]) / L) : q0[4];
        const double gamRot0Z = -Fsz0 * (q0yic + q0yjc);
        gamRotZ = (gammaS[1] / Lc) - gamRot0Z;
      }
      const double ve0 = v(0);
      const double ve1 = vzic - thetaH[0] + gamRotY;
      const double ve2 = vzjc - thetaH[1] + gamRotY;
      const double ve3 = vyic - thetaH[2] + gamRotZ;
      const double ve4 = vyjc - thetaH[3] + gamRotZ;
      const double ve5 = v(5);
      energy += 0.5 * (ea * ve0 * ve0 + gj * ve5 * ve5 +
                       ve1 * (az * ve1 + bz * ve2) + ve2 * (bz * ve1 + az * ve2) +
                       ve3 * (ay * ve3 + by * ve4) + ve4 * (by * ve3 + ay * ve4));
    }
    const double fPM = ladrunoIMKScalePM(q(0), P0, pmAlpha);
    for (int i = 0; i < 4; i++)
      if (theMat[i] != 0) {
        const double f = (dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]) == 0) ? fPM : 1.0;
        energy += f * theMat[i]->getEnergy();
      }
    for (int i = 0; i < 2; i++) {
      if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
        energy += cb->getEnergy();
      if (theShearMat[i] != 0)
        energy += theShearMat[i]->getEnergy();
    }
    return info.setDouble(energy);
  }
  case 12: {
    static thread_local Vector damp(12);
    damp.Zero();
    damp.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return info.setVector(damp);
  }
  case 14: {
    static thread_local Vector sd(2);
    sd(0) = gammaS[0];
    sd(1) = gammaS[1];
    return info.setVector(sd);
  }
  case 15: {
    static thread_local Vector sf(2);
    const double L = theCoordTransf->getInitialLength();
    sf(0) = (L > 0.0) ? -(q(1) + q(2)) / L : 0.0;
    sf(1) = (L > 0.0) ? -(q(3) + q(4)) / L : 0.0;
    return info.setVector(sf);
  }
  case 16: {
    static thread_local Vector sh(4);
    const double L = theCoordTransf->getInitialLength();
    sh(0) = gammaS[0];
    sh(1) = (theShearMat[0] != 0) ? theShearMat[0]->getStress() : ((L > 0.0) ? -(q(1) + q(2)) / L : 0.0);
    sh(2) = gammaS[1];
    sh(3) = (theShearMat[1] != 0) ? theShearMat[1]->getStress() : ((L > 0.0) ? -(q(3) + q(4)) / L : 0.0);
    return info.setVector(sh);
  }
  case 17: {
    static thread_local Vector shy(2);
    const double L = theCoordTransf->getInitialLength();
    shy(0) = gammaS[0];
    shy(1) = (theShearMat[0] != 0) ? theShearMat[0]->getStress() : ((L > 0.0) ? -(q(1) + q(2)) / L : 0.0);
    return info.setVector(shy);
  }
  case 18: {
    static thread_local Vector shz(2);
    const double L = theCoordTransf->getInitialLength();
    shz(0) = gammaS[1];
    shz(1) = (theShearMat[1] != 0) ? theShearMat[1]->getStress() : ((L > 0.0) ? -(q(3) + q(4)) / L : 0.0);
    return info.setVector(shz);
  }
  case 20: {
    static thread_local Vector hbi(4);
    hbi(0) = thetaH[0];
    hbi(1) = q(1);
    hbi(2) = thetaH[2];
    hbi(3) = q(3);
    return info.setVector(hbi);
  }
  case 21: {
    static thread_local Vector hbj(4);
    hbj(0) = thetaH[1];
    hbj(1) = q(2);
    hbj(2) = thetaH[3];
    hbj(3) = q(4);
    return info.setVector(hbj);
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
    int cnt = 0;
    for (int i = 0; i < 4; i++)
      if (theMat[i] != 0) cnt++;
    for (int i = 0; i < 2; i++) {
      if (theBiaxialMat[i] != 0) cnt++;
      if (theShearMat[i] != 0) cnt++;
    }
    return info.setDouble(cnt > 0 ? implexErr / cnt : 0.0);
  }
  case 71:
  case 72: {
    static thread_local Vector sf(4);
    if (responseID == 71) {
      sf(0) = q(0) - p0[0];
      sf(1) = -q(1);
      sf(2) = -q(3);
      sf(3) = q(5);
    } else {
      sf(0) = q(0);
      sf(1) = q(2);
      sf(2) = q(4);
      sf(3) = q(5);
    }
    return info.setVector(sf);
  }
  default:
    return this->Element::getResponse(responseID, info);
  }
}

int LadrunoIMKBeam::setParameter(const char **argv, int argc, Parameter &param)
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
  if (strcmp(argv[0], "Iz") == 0) {
    param.setValue(Iz);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0], "Iy") == 0) {
    param.setValue(Iy);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0], "G") == 0) {
    param.setValue(G);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0], "Jx") == 0 || strcmp(argv[0], "J") == 0) {
    param.setValue(Jx);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(rho);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0], "GAvy") == 0 || strcmp(argv[0], "GAy") == 0) {
    param.setValue(GAvy);
    return param.addObject(8, this);
  }
  if (strcmp(argv[0], "GAvz") == 0 || strcmp(argv[0], "GAz") == 0) {
    param.setValue(GAvz);
    return param.addObject(9, this);
  }
  if (strcmp(argv[0], "P0") == 0) {
    param.setValue(P0);
    return param.addObject(10, this);
  }
  if (strcmp(argv[0], "pmAlpha") == 0) {
    param.setValue(pmAlpha);
    return param.addObject(11, this);
  }

  int res = -1;
  if (strcmp(argv[0], "shearHinge") == 0 || strcmp(argv[0], "shearHingeY") == 0 ||
      strcmp(argv[0], "matShear") == 0 || strcmp(argv[0], "matShearY") == 0) {
    if (argc < 2 || theShearMat[0] == 0)
      return -1;
    res = theShearMat[0]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "shearHingeZ") == 0 || strcmp(argv[0], "matShearZ") == 0) {
    if (argc < 2 || theShearMat[1] == 0)
      return -1;
    res = theShearMat[1]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeBiaxialI") == 0 || strcmp(argv[0], "hingeCI") == 0) {
    if (argc < 2 || theBiaxialMat[0] == 0)
      return -1;
    res = theBiaxialMat[0]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeBiaxialJ") == 0 || strcmp(argv[0], "hingeCJ") == 0) {
    if (argc < 2 || theBiaxialMat[1] == 0)
      return -1;
    res = theBiaxialMat[1]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeBiaxial") == 0 || strcmp(argv[0], "hingeC") == 0) {
    if (argc < 2)
      return -1;
    for (int i = 0; i < 2; i++) {
      if (theBiaxialMat[i] != 0) {
        int r = theBiaxialMat[i]->setParameter(&argv[1], argc - 1, param);
        if (r != -1)
          res = r;
      }
    }
  } else if (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hinge1") == 0 ||
      strcmp(argv[0], "hingeZi") == 0 || strcmp(argv[0], "hingeZI") == 0 ||
      strcmp(argv[0], "hingeMzI") == 0 || strcmp(argv[0], "matZi") == 0) {
    if (argc < 2) return -1;
    if (theBiaxialMat[0] != 0)
      res = theBiaxialMat[0]->setParameter(&argv[1], argc - 1, param);
    else if (theMat[0] != 0)
      res = theMat[0]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hinge2") == 0 ||
      strcmp(argv[0], "hingeZj") == 0 || strcmp(argv[0], "hingeZJ") == 0 ||
      strcmp(argv[0], "hingeMzJ") == 0 || strcmp(argv[0], "matZj") == 0) {
    if (argc < 2) return -1;
    if (theBiaxialMat[1] != 0)
      res = theBiaxialMat[1]->setParameter(&argv[1], argc - 1, param);
    else if (theMat[1] != 0)
      res = theMat[1]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeYi") == 0 || strcmp(argv[0], "hingeYI") == 0 ||
      strcmp(argv[0], "hingeMyI") == 0 || strcmp(argv[0], "matYi") == 0) {
    if (argc < 2 || theMat[2] == 0)
      return -1;
    res = theMat[2]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "hingeYj") == 0 || strcmp(argv[0], "hingeYJ") == 0 ||
      strcmp(argv[0], "hingeMyJ") == 0 || strcmp(argv[0], "matYj") == 0) {
    if (argc < 2 || theMat[3] == 0)
      return -1;
    res = theMat[3]->setParameter(&argv[1], argc - 1, param);
  } else if (strcmp(argv[0], "matZ") == 0 || strcmp(argv[0], "hingeZ") == 0) {
    if (argc < 2)
      return -1;
    for (int i = 0; i < 2; i++) {
      if (theMat[i] != 0) {
        int r = theMat[i]->setParameter(&argv[1], argc - 1, param);
        if (r != -1)
          res = r;
      }
    }
  } else if (strcmp(argv[0], "matY") == 0 || strcmp(argv[0], "hingeY") == 0) {
    if (argc < 2)
      return -1;
    for (int i = 2; i < 4; i++) {
      if (theMat[i] != 0) {
        int r = theMat[i]->setParameter(&argv[1], argc - 1, param);
        if (r != -1)
          res = r;
      }
    }
  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "hinge") == 0) {
    if (argc == 2) {
      for (int i = 0; i < 4; i++) {
        if (theMat[i] != 0) {
          int r = theMat[i]->setParameter(&argv[1], 1, param);
          if (r != -1)
            res = r;
        }
      }
      for (int i = 0; i < 2; i++) {
        if (theBiaxialMat[i] != 0) {
          int r = theBiaxialMat[i]->setParameter(&argv[1], 1, param);
          if (r != -1)
            res = r;
        }
        if (theShearMat[i] != 0) {
          int r = theShearMat[i]->setParameter(&argv[1], 1, param);
          if (r != -1)
            res = r;
        }
      }
    } else if (argc > 2) {
      if ((strcmp(argv[1], "shear") == 0 || strcmp(argv[1], "shearY") == 0) && theShearMat[0] != 0)
        res = theShearMat[0]->setParameter(&argv[2], argc - 2, param);
      else if (strcmp(argv[1], "shearZ") == 0 && theShearMat[1] != 0)
        res = theShearMat[1]->setParameter(&argv[2], argc - 2, param);
      else if ((strcmp(argv[1], "biaxialI") == 0 || strcmp(argv[1], "i") == 0) && theBiaxialMat[0] != 0)
        res = theBiaxialMat[0]->setParameter(&argv[2], argc - 2, param);
      else if ((strcmp(argv[1], "biaxialJ") == 0 || strcmp(argv[1], "j") == 0) && theBiaxialMat[1] != 0)
        res = theBiaxialMat[1]->setParameter(&argv[2], argc - 2, param);
      else {
        int idx = -1;
        if (strcmp(argv[1], "Zi") == 0 || strcmp(argv[1], "i") == 0 || strcmp(argv[1], "0") == 0)
          idx = 0;
        else if (strcmp(argv[1], "Zj") == 0 || strcmp(argv[1], "j") == 0 || strcmp(argv[1], "1") == 0)
          idx = 1;
        else if (strcmp(argv[1], "Yi") == 0 || strcmp(argv[1], "2") == 0)
          idx = 2;
        else if (strcmp(argv[1], "Yj") == 0 || strcmp(argv[1], "3") == 0)
          idx = 3;
        else
          idx = atoi(argv[1]);
        if (idx >= 0 && idx < 4 && theMat[idx] != 0)
          res = theMat[idx]->setParameter(&argv[2], argc - 2, param);
      }
    }
  }

  if (res != -1)
    param.addObject(99, this);
  return res;
}

int LadrunoIMKBeam::updateParameter(int parameterID, Information &info)
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
    Iy = info.theDouble;
    break;
  case 5:
    G = info.theDouble;
    break;
  case 6:
    Jx = info.theDouble;
    break;
  case 7:
    rho = info.theDouble;
    break;
  case 8:
    GAvy = info.theDouble;
    break;
  case 9:
    GAvz = info.theDouble;
    break;
  case 10:
    P0 = info.theDouble;
    if (P0 > 0.0) {
      for (int i = 0; i < 4; i++)
        if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
          ch->setPM(P0, pmAlpha);
      for (int i = 0; i < 2; i++)
        if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
          cb->setPM(P0, pmAlpha);
    }
    break;
  case 11:
    pmAlpha = info.theDouble;
    if (P0 > 0.0) {
      for (int i = 0; i < 4; i++)
        if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(theMat[i]))
          ch->setPM(P0, pmAlpha);
      for (int i = 0; i < 2; i++)
        if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(theBiaxialMat[i]))
          cb->setPM(P0, pmAlpha);
    }
    break;
  case 99:
    break;
  default:
    return -1;
  }
  committedAtCurrentV = false;
  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  return 0;
}

bool LadrunoIMKBeam::ladrunoThreadSafeUpdate(void) const
{
  for (int i = 0; i < 4; i++)
    if (theMat[i] != 0 && !theMat[i]->ladrunoThreadSafeUpdate())
      return false;
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0 && !theBiaxialMat[i]->ladrunoThreadSafeUpdate())
      return false;
    if (theShearMat[i] != 0 && !theShearMat[i]->ladrunoThreadSafeUpdate())
      return false;
  }
  return true;
}

bool LadrunoIMKBeam::ladrunoThreadSafeForm(void) const
{
  for (int i = 0; i < 4; i++)
    if (theMat[i] != 0 && !theMat[i]->ladrunoThreadSafeForm())
      return false;
  for (int i = 0; i < 2; i++) {
    if (theBiaxialMat[i] != 0 && !theBiaxialMat[i]->ladrunoThreadSafeForm())
      return false;
    if (theShearMat[i] != 0 && !theShearMat[i]->ladrunoThreadSafeForm())
      return false;
  }
  return true;
}
