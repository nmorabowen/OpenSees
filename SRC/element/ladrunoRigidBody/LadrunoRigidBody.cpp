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

// Ladruno: RigidBody element (ADR 58) — P1 implementation (translation/ballistic).
// See Ladruno_implementation/58_ladruno_rigid_body_adr.md and the class note in
// LadrunoRigidBody.h. Realized as a zero-stiffness Element (ADR 58 D1(b)) that owns
// a private internal 6-DOF CoM Node carrying the condensed mass and ties its slave
// set to that node by rigid-link MP_Constraints (Joint3D pattern). The element
// assembles NOTHING (K = M = P = C = 0); the dynamics live on the internal node.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoRigidBody.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <math.h>
#include <string.h>

LadrunoRigidBody::LadrunoRigidBody(int tag, int ndmIn, const ID& slaveNodes,
                                   double mUserIn, int internalNodeTag)
  : Element(tag, ELE_TAG_LadrunoRigidBody),
    ndm(ndmIn), nSlave(slaveNodes.Size()), nDOF(0),
    slaveTags(slaveNodes), theSlaves(0),
    intNodeTag(internalNodeTag), theIntNode(0), theMPs(0), theDomainCache(0),
    mUser(mUserIn), mBody(0.0), xc(3), Ibody(3, 3), slaveMass0(0),
    valid(false), K0(0), M0(0), P0(0), C0(0), dampF(0)
{
  xc.Zero();
  Ibody.Zero();
}

LadrunoRigidBody::LadrunoRigidBody()
  : Element(0, ELE_TAG_LadrunoRigidBody),
    ndm(3), nSlave(0), nDOF(0),
    slaveTags(0), theSlaves(0),
    intNodeTag(-1), theIntNode(0), theMPs(0), theDomainCache(0),
    mUser(-1.0), mBody(0.0), xc(3), Ibody(3, 3), slaveMass0(0),
    valid(false), K0(0), M0(0), P0(0), C0(0), dampF(0)
{
  xc.Zero();
  Ibody.Zero();
}

LadrunoRigidBody::~LadrunoRigidBody()
{
  // remove the owned internal node + MP_Constraints and restore slave masses.
  // Use the cached Domain (the base getDomain() is nulled on detach).
  this->teardown(theDomainCache);

  if (theSlaves != 0)   delete[] theSlaves;
  if (theMPs != 0)      delete[] theMPs;
  if (slaveMass0 != 0) {
    for (int i = 0; i < nSlave; i++)
      if (slaveMass0[i] != 0) delete slaveMass0[i];
    delete[] slaveMass0;
  }
  if (K0 != 0)    delete K0;
  if (M0 != 0)    delete M0;
  if (P0 != 0)    delete P0;
  if (C0 != 0)    delete C0;
  if (dampF != 0) delete dampF;
}

// ---- domain ---------------------------------------------------------------

int LadrunoRigidBody::getNumExternalNodes(void) const { return nSlave; }
const ID& LadrunoRigidBody::getExternalNodes(void)    { return slaveTags; }
Node** LadrunoRigidBody::getNodePtrs(void)            { return theSlaves; }
int LadrunoRigidBody::getNumDOF(void)                 { return nDOF; }

void LadrunoRigidBody::setDomain(Domain* theDomain)
{
  // detach (element removed from a domain): drop resolved pointers only
  if (theDomain == 0) {
    if (theSlaves != 0)
      for (int i = 0; i < nSlave; i++) theSlaves[i] = 0;
    theIntNode = 0;
    this->DomainComponent::setDomain(0);
    return;
  }

  this->DomainComponent::setDomain(theDomain);
  theDomainCache = theDomain;          // cache for teardown (getDomain() nulls on detach)

  // idempotent: a second setDomain on an already-built body is a no-op
  if (valid)
    return;

  if (ndm != 3) {
    opserr << "WARNING LadrunoRigidBody " << this->getTag()
           << ": only ndm==3 is supported in v1 (P1)\n";
    return;
  }
  if (nSlave < 1) {
    opserr << "WARNING LadrunoRigidBody " << this->getTag()
           << ": needs at least one slave node\n";
    return;
  }

  // resolve slaves + count element DOFs (slaves must be 3- or 6-DOF in 3D)
  theSlaves = new Node*[nSlave];
  nDOF = 0;
  for (int i = 0; i < nSlave; i++) {
    theSlaves[i] = theDomain->getNode(slaveTags(i));
    if (theSlaves[i] == 0) {
      opserr << "WARNING LadrunoRigidBody " << this->getTag()
             << ": slave node " << slaveTags(i) << " not in domain\n";
      return;
    }
    int nd = theSlaves[i]->getNumberDOF();
    if (nd != 3 && nd != 6) {
      opserr << "WARNING LadrunoRigidBody " << this->getTag()
             << ": slave node " << slaveTags(i) << " has " << nd
             << " DOFs; v1 supports 3 or 6 (3D)\n";
      return;
    }
    nDOF += nd;
  }

  // allocate the zero returns now that nDOF is known, so getMass/getTangentStiff/
  // getResistingForce never dereference a null pointer even if a later step fails
  // (a failed body then assembles inert zeros instead of segfaulting)
  this->allocateZeroReturns();

  // CoM + total mass + body-frame inertia (ADR 58 D3), and zero the slaves'
  // nodal mass so it is not double-counted once condensed onto the CoM node
  if (this->resolveCentroidAndMass() != 0)   // warns on failure
    return;

  if (this->createInternalNode(theDomain) < 0) return;
  if (this->buildSlaveConstraints(theDomain) < 0) return;

  valid = true;
}

// ---- mass / centroid condensation (ADR 58 D3) -----------------------------

int LadrunoRigidBody::resolveCentroidAndMass(void)
{
  // pass 1: total mass + first moment → CoM
  double M = 0.0;
  Vector S(3); S.Zero();
  for (int i = 0; i < nSlave; i++) {
    const Vector& xi = theSlaves[i]->getCrds();
    double mi = (mUser >= 0.0)
                  ? mUser / nSlave                       // explicit total mass, split evenly
                  : theSlaves[i]->getMass()(0, 0);       // slave lumped translational mass
    M += mi;
    for (int k = 0; k < 3; k++) S(k) += mi * xi(k);
  }
  if (M <= 0.0) {
    opserr << "WARNING LadrunoRigidBody " << this->getTag()
           << ": total mass is zero — give the slaves nodal mass or pass -mass m\n";
    return -1;
  }
  mBody = M;
  for (int k = 0; k < 3; k++) xc(k) = S(k) / M;

  // pass 2: inertia tensor about the CoM, I = Σ m_i (|d|² I₃ − d dᵀ)
  Ibody.Zero();
  double Lc2 = 0.0;
  for (int i = 0; i < nSlave; i++) {
    const Vector& xi = theSlaves[i]->getCrds();
    double mi = (mUser >= 0.0) ? mUser / nSlave : theSlaves[i]->getMass()(0, 0);
    double d0 = xi(0) - xc(0), d1 = xi(1) - xc(1), d2 = xi(2) - xc(2);
    double r2 = d0 * d0 + d1 * d1 + d2 * d2;
    if (r2 > Lc2) Lc2 = r2;
    Ibody(0, 0) += mi * (r2 - d0 * d0);
    Ibody(1, 1) += mi * (r2 - d1 * d1);
    Ibody(2, 2) += mi * (r2 - d2 * d2);
    Ibody(0, 1) -= mi * d0 * d1;
    Ibody(0, 2) -= mi * d0 * d2;
    Ibody(1, 2) -= mi * d1 * d2;
  }
  Ibody(1, 0) = Ibody(0, 1);
  Ibody(2, 0) = Ibody(0, 2);
  Ibody(2, 1) = Ibody(1, 2);

  // floor the inertia so a collinear/coincident slave set does not leave a
  // singular rotational DOF on the internal node (rotation is unused in P1, but
  // the central-difference solve still inverts the node mass). P2 replaces this
  // with the body-frame side-channel integrator (ADR 58 D3).
  double floorI = 1.0e-9 * mBody * Lc2;
  double absFloor = 1.0e-12 * mBody;
  if (floorI < absFloor) floorI = absFloor;
  for (int k = 0; k < 3; k++) Ibody(k, k) += floorI;

  // zero the slaves' nodal mass (captured for restore at teardown), so the mass
  // now lives only at the CoM node — only when we read it from the slaves.
  slaveMass0 = new Matrix*[nSlave];
  for (int i = 0; i < nSlave; i++) slaveMass0[i] = 0;
  if (mUser < 0.0) {
    for (int i = 0; i < nSlave; i++) {
      const Matrix& mm = theSlaves[i]->getMass();
      slaveMass0[i] = new Matrix(mm);                 // capture original
      Matrix zero(mm.noRows(), mm.noCols());
      zero.Zero();
      theSlaves[i]->setMass(zero);
    }
  }
  return 0;
}

// ---- internal node + slave constraints (Joint3D pattern) ------------------

int LadrunoRigidBody::createInternalNode(Domain* theDomain)
{
  // pick a free internal-node tag (user-supplied or auto, bumped on collision)
  int base = (intNodeTag >= 0) ? intNodeTag : (9000000 + this->getTag());
  while (theDomain->getNode(base) != 0) base++;
  intNodeTag = base;

  theIntNode = new Node(intNodeTag, 6, xc(0), xc(1), xc(2));
  if (theIntNode == 0 || theDomain->addNode(theIntNode) == false) {
    opserr << "WARNING LadrunoRigidBody " << this->getTag()
           << ": could not create/add internal CoM node " << intNodeTag << "\n";
    return -1;
  }

  // condensed mass on the internal node: diag(m,m,m) ⊕ Ibody
  Matrix Mc(6, 6); Mc.Zero();
  Mc(0, 0) = Mc(1, 1) = Mc(2, 2) = mBody;
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++) Mc(3 + a, 3 + b) = Ibody(a, b);
  theIntNode->setMass(Mc);
  return 0;
}

int LadrunoRigidBody::buildSlaveConstraints(Domain* theDomain)
{
  theMPs = new MP_Constraint*[nSlave];
  for (int i = 0; i < nSlave; i++) theMPs[i] = 0;

  for (int i = 0; i < nSlave; i++) {
    int nd = theSlaves[i]->getNumberDOF();
    const Vector& xi = theSlaves[i]->getCrds();
    double d0 = xi(0) - xc(0), d1 = xi(1) - xc(1), d2 = xi(2) - xc(2);

    // U_slave = Ccr · U_internal :  u = u_R + θ_R × d  (rotations tied if nd==6)
    Matrix Ccr(nd, 6); Ccr.Zero();
    Ccr(0, 0) = 1.0; Ccr(1, 1) = 1.0; Ccr(2, 2) = 1.0;
    Ccr(0, 4) = d2;  Ccr(0, 5) = -d1;
    Ccr(1, 3) = -d2; Ccr(1, 5) = d0;
    Ccr(2, 3) = d1;  Ccr(2, 4) = -d0;
    if (nd >= 6) { Ccr(3, 3) = 1.0; Ccr(4, 4) = 1.0; Ccr(5, 5) = 1.0; }

    ID constrDOF(nd);
    for (int k = 0; k < nd; k++) constrDOF(k) = k;
    ID retainDOF(6);
    for (int k = 0; k < 6; k++) retainDOF(k) = k;

    MP_Constraint* mp = new MP_Constraint(intNodeTag, slaveTags(i), Ccr,
                                          constrDOF, retainDOF);
    if (mp == 0 || theDomain->addMP_Constraint(mp) == false) {
      opserr << "WARNING LadrunoRigidBody " << this->getTag()
             << ": could not add rigid-link MP for slave " << slaveTags(i) << "\n";
      if (mp != 0) delete mp;
      return -1;
    }
    theMPs[i] = mp;
  }
  return 0;
}

void LadrunoRigidBody::teardown(Domain* theDomain)
{
  if (theDomain == 0) return;

  if (theMPs != 0) {
    for (int i = 0; i < nSlave; i++) {
      if (theMPs[i] != 0) {
        int t = theMPs[i]->getTag();
        MP_Constraint* mp = theDomain->getMP_Constraint(t);
        if (mp != 0) { theDomain->removeMP_Constraint(t); delete mp; }
        theMPs[i] = 0;
      }
    }
  }
  // restore the slaves' captured nodal mass — fetch by tag (theSlaves may dangle
  // during a Domain wipe; a tag lookup returns 0 if the node is already gone)
  if (slaveMass0 != 0) {
    for (int i = 0; i < nSlave; i++) {
      if (slaveMass0[i] != 0) {
        Node* sn = theDomain->getNode(slaveTags(i));
        if (sn != 0) sn->setMass(*slaveMass0[i]);
      }
    }
  }
  if (theIntNode != 0) {
    Node* n = theDomain->removeNode(intNodeTag);
    if (n != 0) delete n;
    theIntNode = 0;
  }
}

void LadrunoRigidBody::allocateZeroReturns(void)
{
  K0 = new Matrix(nDOF, nDOF); K0->Zero();
  M0 = new Matrix(nDOF, nDOF); M0->Zero();
  C0 = new Matrix(nDOF, nDOF); C0->Zero();
  P0 = new Vector(nDOF);       P0->Zero();
  dampF = new Vector(nDOF);    dampF->Zero();
}

// ---- state (P1: trivial; P2 advances the body-frame SO(3) state here) ------
//
// P1 owns no element-level state: the internal Node commits/reverts itself, so
// these are no-ops. P2 ADDS body-frame SO(3) state (orientation quaternion +
// angular velocity) and MUST save it in commitState and restore it in
// revertToLastCommit/revertToStart (TODO P2 — do not leave these as no-ops once
// rotational state exists).

int LadrunoRigidBody::commitState(void)        { return 0; }
int LadrunoRigidBody::revertToLastCommit(void) { return 0; }
int LadrunoRigidBody::revertToStart(void)      { return 0; }
int LadrunoRigidBody::update(void)             { return 0; }

// ---- matrices / forces: a rigid body assembles none of its own -------------

const Matrix& LadrunoRigidBody::getTangentStiff(void) { return *K0; }
const Matrix& LadrunoRigidBody::getInitialStiff(void) { return *K0; }
const Matrix& LadrunoRigidBody::getMass(void)         { return *M0; }
const Vector& LadrunoRigidBody::getResistingForce(void)            { return *P0; }
const Vector& LadrunoRigidBody::getResistingForceIncInertia(void)  { return *P0; }

// no physical damping (a βK must not shrink dt_cr; the base lazy-damping slot is
// never allocated, so getDamp/getRayleighDampingForces must be overridden — same
// landmine the RBE2 sibling documents).
int LadrunoRigidBody::setRayleighDampingFactors(double, double, double, double) { return 0; }
const Matrix& LadrunoRigidBody::getDamp(void)                 { return *C0; }
const Vector& LadrunoRigidBody::getRayleighDampingForces(void){ return *dampF; }

// ---- parallel (serial-only v1 — ADR 58 §6) --------------------------------

int LadrunoRigidBody::sendSelf(int, Channel&)
{
  opserr << "LadrunoRigidBody::sendSelf - not implemented (serial-only v1)\n";
  return -1;
}
int LadrunoRigidBody::recvSelf(int, Channel&, FEM_ObjectBroker&)
{
  opserr << "LadrunoRigidBody::recvSelf - not implemented (serial-only v1)\n";
  return -1;
}

// ---- recorder / output (ADR 58 D9): CoM kinematics off the internal node ---

Response* LadrunoRigidBody::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  if (strcmp(argv[0], "comDisp")  == 0 || strcmp(argv[0], "centroidDisp") == 0)
    return new ElementResponse(this, 1, Vector(3));
  if (strcmp(argv[0], "comVel")   == 0)
    return new ElementResponse(this, 2, Vector(3));
  if (strcmp(argv[0], "comAccel") == 0)
    return new ElementResponse(this, 3, Vector(3));
  if (strcmp(argv[0], "mass")     == 0)
    return new ElementResponse(this, 4, Vector(1));
  return 0;
}

int LadrunoRigidBody::getResponse(int responseID, Information& eleInfo)
{
  Vector out3(3);
  if (theIntNode == 0) return -1;
  switch (responseID) {
  case 1: {
    const Vector& u = theIntNode->getTrialDisp();
    for (int k = 0; k < 3; k++) out3(k) = u(k);
    return eleInfo.setVector(out3);
  }
  case 2: {
    const Vector& v = theIntNode->getTrialVel();
    for (int k = 0; k < 3; k++) out3(k) = v(k);
    return eleInfo.setVector(out3);
  }
  case 3: {
    const Vector& a = theIntNode->getTrialAccel();
    for (int k = 0; k < 3; k++) out3(k) = a(k);
    return eleInfo.setVector(out3);
  }
  case 4: {
    Vector m(1); m(0) = mBody;
    return eleInfo.setVector(m);
  }
  default:
    return -1;
  }
}

void LadrunoRigidBody::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoRigidBody, tag: " << this->getTag() << "\n";
  s << "  slaves: " << slaveTags;
  s << "  internal CoM node: " << intNodeTag
    << "  CoM: " << xc(0) << " " << xc(1) << " " << xc(2)
    << "  m_body: " << mBody << "\n";
}
