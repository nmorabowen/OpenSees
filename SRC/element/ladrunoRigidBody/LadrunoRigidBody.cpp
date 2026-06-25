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
#include <Element.h>
#include <ElementIter.h>
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

// Rotate a vector by a Versor: v' = v + 2 qs (qv×v) + 2 qv×(qv×v). Implemented
// locally because Versor::rotate uses left-scalar multiplication (2.0*vec /
// scalar*vec) which VectorND<3> does not provide (only vec*double + one-arg cross).
static Vector3D rotateVec(const OpenSees::Versor& q, const Vector3D& v)
{
  Vector3D c1 = q.vector.cross(v);     // qv × v
  Vector3D c2 = q.vector.cross(c1);    // qv × (qv × v)
  Vector3D r;
  for (int i = 0; i < 3; i++)
    r[i] = v[i] + 2.0 * q.scalar * c1[i] + 2.0 * c2[i];
  return r;
}

// a slave's lumped translational mass = average of the 3 translational diagonals
// (robust to anisotropic nodal mass; for the usual isotropic mass they are equal)
static double lumpedTransMass(Node* n)
{
  const Matrix& m = n->getMass();
  if (m.noRows() < 3) return 0.0;
  return (m(0, 0) + m(1, 1) + m(2, 2)) / 3.0;
}

LadrunoRigidBody::LadrunoRigidBody(int tag, int ndmIn, const ID& slaveNodes,
                                   double mUserIn, int internalNodeTag)
  : Element(tag, ELE_TAG_LadrunoRigidBody),
    ndm(ndmIn), nSlave(slaveNodes.Size()), nDOF(0),
    slaveTags(slaveNodes), theSlaves(0),
    intNodeTag(internalNodeTag), theIntNode(0), theMPs(0), theDomainCache(0),
    mUser(mUserIn), mBody(0.0), xc(3), Ibody(3, 3), slaveMass0(0),
    slaveOffset0(0), nIncident(0), incidentTag(0), incidentOff(0), incidentSlaveIdx(0),
    haveOmega0(false), IbodyInv(3, 3),
    lastTimeCommit(0.0), lastTimeTrial(0.0), started(false),
    valid(false), K0(0), M0(0), P0(0), C0(0), dampF(0)
{
  xc.Zero();
  Ibody.Zero();
  qCommit.vector[0] = qCommit.vector[1] = qCommit.vector[2] = 0.0;
  qCommit.scalar = 1.0;
  qTrial = qCommit;
  for (int k = 0; k < 3; k++) { Lcommit[k] = Ltrial[k] = L0[k] = omega0[k] = 0.0; }
  IbodyInv.Zero();
}

LadrunoRigidBody::LadrunoRigidBody()
  : Element(0, ELE_TAG_LadrunoRigidBody),
    ndm(3), nSlave(0), nDOF(0),
    slaveTags(0), theSlaves(0),
    intNodeTag(-1), theIntNode(0), theMPs(0), theDomainCache(0),
    mUser(-1.0), mBody(0.0), xc(3), Ibody(3, 3), slaveMass0(0),
    slaveOffset0(0), nIncident(0), incidentTag(0), incidentOff(0), incidentSlaveIdx(0),
    haveOmega0(false), IbodyInv(3, 3),
    lastTimeCommit(0.0), lastTimeTrial(0.0), started(false),
    valid(false), K0(0), M0(0), P0(0), C0(0), dampF(0)
{
  xc.Zero();
  Ibody.Zero();
  qCommit.vector[0] = qCommit.vector[1] = qCommit.vector[2] = 0.0;
  qCommit.scalar = 1.0;
  qTrial = qCommit;
  for (int k = 0; k < 3; k++) { Lcommit[k] = Ltrial[k] = L0[k] = omega0[k] = 0.0; }
  IbodyInv.Zero();
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
  if (slaveOffset0 != 0)    delete[] slaveOffset0;
  if (incidentTag != 0)     delete[] incidentTag;       // tags, not owned elements
  if (incidentOff != 0)     delete[] incidentOff;
  if (incidentSlaveIdx != 0) delete[] incidentSlaveIdx;
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

  // idempotent: a second setDomain on an already-built body keeps its state, but
  // REBUILDS the incident-element cache so a domainChanged after late-added toe/
  // contact elements still picks them up (M1 — the ordering precondition is then a
  // best-effort, not the only chance).
  if (valid) {
    this->buildIncidentCache(theDomain);
    return;
  }

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

  // P2 SO(3) seed: invert the constant body-frame inertia for the side channel,
  // and set the initial spatial angular momentum from -omega (body frame, q0=I ⇒
  // L0 = Ibody·omega0). The orientation starts at identity.
  if (!invert3x3(Ibody, IbodyInv)) {
    opserr << "WARNING LadrunoRigidBody " << this->getTag()
           << ": body-frame inertia is not invertible\n";
    return;
  }
  for (int i = 0; i < 3; i++) {
    double Li = 0.0;
    for (int j = 0; j < 3; j++) Li += Ibody(i, j) * omega0[j];
    L0[i] = Li; Lcommit[i] = Li; Ltrial[i] = Li;
  }
  lastTimeCommit = lastTimeTrial = theDomain->getCurrentTime();
  started = false;

  // P2 Stage 2: capture each slave's spatial lever arm d_i^0 = x_i^0 − xc (the
  // body-frame offset at the reference config), and cache the external elements
  // attached to the slaves so their reactions can drive the spin (moment gather).
  slaveOffset0 = new Vector3D[nSlave];
  for (int i = 0; i < nSlave; i++) {
    const Vector& xi = theSlaves[i]->getCrds();
    for (int k = 0; k < 3; k++) slaveOffset0[i][k] = xi(k) - xc(k);
  }
  this->buildIncidentCache(theDomain);

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
                  : lumpedTransMass(theSlaves[i]);       // slave lumped translational mass
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
    double mi = (mUser >= 0.0) ? mUser / nSlave : lumpedTransMass(theSlaves[i]);
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

// ---- state: body-frame SO(3) integration (ADR 58 D2/D3, P2 Stage 1) --------
//
// Once per converged step, advance the rotation OFF the global solve:
//   ω_body = Ibody⁻¹·(qᵀ·L)        (current body-frame angular velocity)
//   q      ← q · exp(ω_body·dt)    (orientation, incremental exp-map compose)
//   L      ← L + m·dt              (spatial angular momentum; torque-free ⇒ L
//                                   unchanged ⇒ ‖L‖ conserved to machine precision)
// The Dzhanibekov flip emerges from ω_body re-derived as q rotates, even though L
// is constant. dt is the domain-time advance since the last commit.

int LadrunoRigidBody::commitState(void)
{
  if (valid && theDomainCache != 0) {
    double t = theDomainCache->getCurrentTime();
    if (!started) {                  // first commit just seeds the clock (no giant first dt)
      lastTimeTrial = t;
      started = true;
    }
    double dt = t - lastTimeTrial;
    if (started && dt > 0.0) {
      // P2 Stage 2: gather the external CoM torque at the CURRENT (start-of-step)
      // orientation — the configuration the slaves sit at and at which the toe/contact
      // force was just evaluated. Done BEFORE the q advance so the lever arm R·d_i^0
      // matches the force's config (an explicit forward torque; O(dt) consistent).
      Vector3D m = this->gatherCoMTorque();

      // 2nd-order midpoint orientation update (Krysl–Endres-style): sample ω at the
      // half-step orientation so energy stays bounded (a forward sample at step-start
      // is dissipative and spirals a free body to its major axis instead of flipping).
      Vector3D w0 = this->omegaBodyFrom(qTrial, Ltrial);
      Vector3D dthH;
      for (int k = 0; k < 3; k++) dthH[k] = w0[k] * (0.5 * dt);
      OpenSees::Versor qHalf = qTrial * OpenSees::Versor::from_vector(dthH);
      qHalf.normalize();

      Vector3D wH = this->omegaBodyFrom(qHalf, Ltrial);      // ω at the midpoint
      Vector3D dth;
      for (int k = 0; k < 3; k++) dth[k] = wH[k] * dt;
      qTrial = qTrial * OpenSees::Versor::from_vector(dth);  // q · exp(ω_mid·dt)
      qTrial.normalize();

      // kick the spatial angular momentum with the gathered torque. Free spin ⇒ m≈0
      // ⇒ ‖L‖ conserved (the Dzhanibekov gate). A toe reaction / applied moment now
      // drives the spin (P2 Stage 2, no longer ignored).
      for (int k = 0; k < 3; k++) Ltrial[k] += m[k] * dt;
      lastTimeTrial = t;
    }
  }
  qCommit = qTrial;
  Lcommit = Ltrial;
  lastTimeCommit = lastTimeTrial;
  return 0;
}

int LadrunoRigidBody::revertToLastCommit(void)
{
  qTrial = qCommit;
  Ltrial = Lcommit;
  lastTimeTrial = lastTimeCommit;     // keep the clock in sync with the rolled-back state
  return 0;
}

int LadrunoRigidBody::revertToStart(void)
{
  qCommit.vector[0] = qCommit.vector[1] = qCommit.vector[2] = 0.0;
  qCommit.scalar = 1.0;
  qTrial = qCommit;
  Lcommit = L0; Ltrial = L0;
  lastTimeTrial = lastTimeCommit;
  started = false;                    // re-baseline the clock on the next commit
  return 0;
}

// update() runs inside Domain::update(), which the integrator calls AFTER the
// constraint handler imposes the linear MP part (u_i = u_R) and BEFORE the residual
// is formed (CentralDifference::newStep → applyLoad/enforceSPs → Domain::update).
// We overwrite the slaves with the exact finite-rotation kinematics so the toe
// contact sees the rotated geometry when it forms its resisting force.
int LadrunoRigidBody::update(void)
{
  this->imposeSlaveKinematics();
  return 0;
}

// impose u_i = u_R + (R−I)·d_i^0 on every slave (ADR 58 §3, the exact finite-rotation
// slave map). The per-slave MP already carries the linear u_R; here we ADD the
// nonlinear (R−I)d_i^0 rotation term that the linear/homogeneous MP cannot represent
// (the orientation lives in the side channel q, not in any retained DOF).
void LadrunoRigidBody::imposeSlaveKinematics(void)
{
  if (!valid || theIntNode == 0) return;
  const Vector& uR = theIntNode->getTrialDisp();      // internal CoM node, 6 DOF
  const Vector& vR = theIntNode->getTrialVel();        // CoM translational velocity (0..2)
  // spatial angular velocity at the imposed orientation (M2): ω_s = R·(Ibody⁻¹ qᵀL)
  Vector3D wSpat = rotateVec(qTrial, this->omegaBodyFrom(qTrial, Ltrial));

  // 6-DOF slaves: body rotation vector from the quaternion log (set once — same body q)
  double nv = sqrt(qTrial.vector[0]*qTrial.vector[0]
                 + qTrial.vector[1]*qTrial.vector[1]
                 + qTrial.vector[2]*qTrial.vector[2]);
  double angle = 2.0 * atan2(nv, qTrial.scalar);

  for (int i = 0; i < nSlave; i++) {
    Vector3D r = rotateVec(qTrial, slaveOffset0[i]);   // R·d_i^0 (current spatial lever arm)
    for (int k = 0; k < 3; k++) {
      double u = uR(k) + (r[k] - slaveOffset0[i][k]);  // u_R + (R−I)·d_i^0
      theSlaves[i]->setTrialDisp(u, k);
    }

    // M2 — consistent slave VELOCITY v_i = v_R + ω × (R·d_i^0), so a velocity-dependent
    // incident element (dashpot / viscous contact) and a slave vel recorder stay
    // coherent with the imposed finite-rotation disp (the handler would otherwise leave
    // the linearized T·v_retained). setTrialVel needs the full nodal vector.
    int nd = theSlaves[i]->getNumberDOF();
    Vector3D wxr = wSpat.cross(r);
    Vector vel(nd); vel.Zero();
    for (int k = 0; k < 3; k++) vel(k) = vR(k) + wxr[k];

    // 6-DOF slaves: also drive the rotational DOFs (disp = body rotation vector, vel =
    // spatial ω). The MP rotation block tied them to the off-channel internal-node
    // rotation (~0); not exercised by the 3-DOF toe of the Housner gate.
    if (nd >= 6) {
      for (int k = 0; k < 3; k++) {
        double th = (nv > 1.0e-12) ? angle * qTrial.vector[k] / nv : 0.0;
        theSlaves[i]->setTrialDisp(th, 3 + k);
        vel(3 + k) = wSpat[k];
      }
    }
    theSlaves[i]->setTrialVel(vel);
  }
}

// cache the external elements attached to the slaves (ADR 58 D3, the dual of slaving):
// the toe/contact reaction on a slave is read off these elements' getResistingForce()
// in commitState to drive the spin. Built ONCE here, in the safe setDomain context —
// NOT lazily from commitState, which runs inside Domain::commit()'s element loop over
// the shared SingleDomEleIter (a reentrant getElements()→reset() there would abort the
// commit loop, silently skipping commitState on every later element).
void LadrunoRigidBody::buildIncidentCache(Domain* theDomain)
{
  // free any prior cache (this is re-runnable on a second setDomain)
  if (incidentTag != 0)      { delete[] incidentTag;      incidentTag = 0; }
  if (incidentOff != 0)      { delete[] incidentOff;      incidentOff = 0; }
  if (incidentSlaveIdx != 0) { delete[] incidentSlaveIdx; incidentSlaveIdx = 0; }
  nIncident = 0;

  for (int pass = 0; pass < 2; pass++) {
    int count = 0;
    ElementIter& eiter = theDomain->getElements();
    Element* e;
    while ((e = eiter()) != 0) {
      if (e->getTag() == this->getTag()) continue;     // skip self
      const ID& enodes = e->getExternalNodes();
      int off = 0;
      for (int j = 0; j < enodes.Size(); j++) {
        Node* nj = theDomain->getNode(enodes(j));
        int ndj = (nj != 0) ? nj->getNumberDOF() : 0;
        for (int s = 0; s < nSlave; s++) {
          if (enodes(j) == slaveTags(s)) {
            // one record per (element, slave) incidence: if an element spans two
            // slaves of THIS body it is recorded twice (once per slave), and the
            // gather sums both r_i × f_i — the correct total moment from that element.
            if (pass == 1) {
              incidentTag[count]      = e->getTag();
              incidentOff[count]      = off;
              incidentSlaveIdx[count] = s;
            }
            count++;
            break;
          }
        }
        off += ndj;
      }
    }
    if (pass == 0) {
      nIncident = count;
      if (nIncident == 0) return;
      incidentTag      = new int[nIncident];
      incidentOff      = new int[nIncident];
      incidentSlaveIdx = new int[nIncident];
    }
  }
}

// net spatial moment about the CoM: m = Σ (R·d_i^0) × f_i, with f_i the EXTERNAL force
// on slave i (toe/contact reaction). The reaction at a constrained massless slave is
// the force the rigid-link carries; the force the slave exerts on the BODY is its
// negative (Newton's third law vs the calculateNodalReactions sign convention, where
// reaction = Σ getResistingForce). Lever arm uses the CURRENT orientation qTrial
// (the config at which the slaves were imposed and the toe force evaluated) — a frozen
// d_i^0 would keep a constant restoring moment and destroy the rocking equilibrium at
// θ=α. The translational part of f_i is ALREADY gathered to the CoM by the MP
// transpose (drives CoM translation); only this moment is lost to the side channel.
Vector3D LadrunoRigidBody::gatherCoMTorque(void)
{
  Vector3D m; m[0] = m[1] = m[2] = 0.0;
  if (theIntNode == 0 || theDomainCache == 0) return m;   // detached — no force source
  for (int rr = 0; rr < nIncident; rr++) {
    // re-resolve by tag so a removed incident element (this fork's element-removal
    // feature) is skipped, never dereferenced after free.
    Element* e = theDomainCache->getElement(incidentTag[rr]);
    if (e == 0) continue;
    const Vector& P = e->getResistingForce();    // committed-config resisting force
    int off = incidentOff[rr];
    if (off + 2 >= P.Size()) continue;
    Vector3D f;
    for (int k = 0; k < 3; k++) f[k] = -P(off + k);    // force ON the body from the slave
    Vector3D lev = rotateVec(qTrial, slaveOffset0[incidentSlaveIdx[rr]]);  // R·d_i^0
    Vector3D mi = lev.cross(f);
    for (int k = 0; k < 3; k++) m[k] += mi[k];
  }
  return m;
}

// body-frame angular velocity from orientation + spatial momentum: ω = Ibody⁻¹·(qᵀ·L)
Vector3D LadrunoRigidBody::omegaBodyFrom(const OpenSees::Versor& q, const Vector3D& L) const
{
  Vector3D piBody = rotateVec(q.conjugate(), L);
  Vector3D w;
  for (int i = 0; i < 3; i++) {
    double wi = 0.0;
    for (int j = 0; j < 3; j++) wi += IbodyInv(i, j) * piBody[j];
    w[i] = wi;
  }
  return w;
}

void LadrunoRigidBody::setOmega0(double wx, double wy, double wz)
{
  omega0[0] = wx; omega0[1] = wy; omega0[2] = wz;
  haveOmega0 = true;
}

// analytic inverse of a 3×3 (symmetric) matrix; false if ~singular
bool LadrunoRigidBody::invert3x3(const Matrix& A, Matrix& Ainv)
{
  double a = A(0,0), b = A(0,1), c = A(0,2);
  double d = A(1,0), e = A(1,1), f = A(1,2);
  double g = A(2,0), h = A(2,1), i = A(2,2);
  double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  if (det == 0.0 || !(det == det)) return false;   // singular / NaN
  double inv = 1.0 / det;
  Ainv(0,0) = (e*i - f*h) * inv;
  Ainv(0,1) = (c*h - b*i) * inv;
  Ainv(0,2) = (b*f - c*e) * inv;
  Ainv(1,0) = (f*g - d*i) * inv;
  Ainv(1,1) = (a*i - c*g) * inv;
  Ainv(1,2) = (c*d - a*f) * inv;
  Ainv(2,0) = (d*h - e*g) * inv;
  Ainv(2,1) = (b*g - a*h) * inv;
  Ainv(2,2) = (a*e - b*d) * inv;
  return true;
}

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
  if (strcmp(argv[0], "orientation") == 0 || strcmp(argv[0], "quaternion") == 0)
    return new ElementResponse(this, 5, Vector(4));   // (qx, qy, qz, qs)
  if (strcmp(argv[0], "angularVel") == 0 || strcmp(argv[0], "omega") == 0)
    return new ElementResponse(this, 6, Vector(3));   // spatial ω
  if (strcmp(argv[0], "angularMom") == 0 || strcmp(argv[0], "L") == 0)
    return new ElementResponse(this, 7, Vector(3));   // spatial angular momentum
  if (strcmp(argv[0], "omegaBody") == 0 || strcmp(argv[0], "angularVelBody") == 0)
    return new ElementResponse(this, 8, Vector(3));   // body-frame ω (flip metric)
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
  case 5: {                                   // orientation quaternion (qx,qy,qz,qs)
    Vector q4(4);
    q4(0) = qCommit.vector[0]; q4(1) = qCommit.vector[1];
    q4(2) = qCommit.vector[2]; q4(3) = qCommit.scalar;
    return eleInfo.setVector(q4);
  }
  case 6: {                                   // spatial angular velocity ω = q·(I⁻¹ qᵀL)
    Vector3D wBody = this->omegaBodyFrom(qCommit, Lcommit);
    Vector3D wSpat = rotateVec(qCommit, wBody);
    for (int k = 0; k < 3; k++) out3(k) = wSpat[k];
    return eleInfo.setVector(out3);
  }
  case 7: {                                   // spatial angular momentum L
    for (int k = 0; k < 3; k++) out3(k) = Lcommit[k];
    return eleInfo.setVector(out3);
  }
  case 8: {                                   // BODY-frame ω (the Dzhanibekov flip metric)
    Vector3D wBody = this->omegaBodyFrom(qCommit, Lcommit);
    for (int k = 0; k < 3; k++) out3(k) = wBody[k];
    return eleInfo.setVector(out3);
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
