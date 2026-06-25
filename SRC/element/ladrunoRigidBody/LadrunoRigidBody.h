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

// Ladruno: RigidBody element (ADR 58). See
// Ladruno_implementation/58_ladruno_rigid_body_adr.md.
//
//   A genuine 6-DOF RIGID BODY: a part of the model that moves as one rigid
//   object (CoM translation + finite rotation) instead of being meshed as
//   deformable. Motivating use case = rocking foundations (a footing/wall that
//   tips about its toe on a compression-only contact), base isolation, rigid pile
//   caps, multi-body. = LS-DYNA *CONSTRAINED_NODAL_RIGID_BODY / Abaqus *RIGID BODY.
//
//   ARCHITECTURE — zero-stiffness Element (ADR 58 D1(b), the fired fallback):
//   the body is realized as an Element that OWNS a PRIVATE INTERNAL 6-DOF Node at
//   the center of mass (the Joint3D pattern, SRC/element/joint/Joint3D.cpp) and
//   ties its slave node set to that internal node by internal MP_Constraints. The
//   element contributes NO stiffness, NO resisting force, and NO mass of its own:
//       getTangentStiff = getInitialStiff = 0,  getResistingForce = 0,
//       getMass = 0   (the CONDENSED diagonal translational mass lives on the
//                      internal CoM Node; slaves are pure kinematic followers).
//   Chosen over a standalone DomainComponent because (a) the framework footprint
//   is 6 files vs ~15+ (Joint3D touches no Domain/broker quartet), and (b) — the
//   decider — Domain::commit()/revertToLastCommit() iterate Nodes+Elements ONLY
//   (Domain.cpp:2183-2237), so a DomainComponent would get no per-step state
//   callback; an Element inherits commitState()/revertToLastCommit()/update() for
//   free, which is exactly where the P2 SO(3) self-integration lands.
//
//   MASS — condensation to the CoM (ADR 58 D3). The slaves' translational masses
//   are lumped to a diagonal translational mass m_body on the internal node, and
//   their inertia about the CoM to a 3×3 tensor Ibody (m_body = Σ m_i, CoM =
//   Σ m_i x_i / Σ m_i). The TRANSLATIONAL block is diagonal (DiagonalSOE-valid). In
//   P1 the full 6×6 (incl. the dense Ibody block) is set on the internal node;
//   under system('Diagonal') the off-diagonal products of inertia are harmlessly
//   dropped (rotation is unused in P1). P2 moves the rotational integration into the
//   commitState side channel so a non-principal-axis body integrates correctly
//   regardless of the global SOE (ADR 58 D3). dt_cr: the body carries no element
//   stiffness, so it never enters the CFL eigensolve.
//
//   KINEMATICS — exact finite-rotation slave map (ADR 58 §3, corrected):
//       u_i = u_R + (R − I) d_i ,   d_i = x_i − x_R .
//   P1 (this milestone) freezes R = I (translation/ballistic only): the map
//   collapses to u_i = u_R, i.e. a rigid-link MP per slave, and the milestone
//   proves the internal-node + MP-slaving + condensed-mass + load/IC/recorder
//   plumbing end to end (free-fall = ballistic; followers exact). P2 adds the
//   exp-map SO(3) update (reuse SRC/matrix/GroupSO3.h + Versor.h) advanced in
//   commitState, re-forming the transport block of each slave MP as R evolves.
//
//   EXTERNAL nodes = the slave set [ S_1 .. S_N ]. The internal CoM node is added
//   to the Domain in setDomain and removed (with its MP_Constraints) in the
//   destructor / when detached — guarded against already-gone nodes (Joint3D).
//
//   classTag ELE_TAG_LadrunoRigidBody = 33015 (next free ELE slot after 33014).
//   Explicit-only v1 (no implicit tangent — ADR 58 §3). Serial-only v1.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoRigidBody_h
#define LadrunoRigidBody_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>     // OpenSees::Versor — unit-quaternion orientation (P2 SO(3))
#include <Vector3D.h>   // body-frame angular velocity / momentum (P2)

class Node;
class Domain;
class MP_Constraint;
class Channel;
class FEM_ObjectBroker;
class Response;

class LadrunoRigidBody : public Element
{
 public:
  // tag, spatial dim, the slave node set, an optional explicit mass/inertia
  // override (mUser>=0 ⇒ use it instead of condensing the slaves' nodal mass).
  LadrunoRigidBody(int tag, int ndm, const ID& slaveNodes,
                   double mUser = -1.0, int internalNodeTag = -1);
  LadrunoRigidBody();
  ~LadrunoRigidBody();

  const char* getClassType(void) const { return "LadrunoRigidBody"; }

  // domain
  int getNumExternalNodes(void) const;
  const ID& getExternalNodes(void);
  Node** getNodePtrs(void);
  int getNumDOF(void);
  void setDomain(Domain* theDomain);

  // initial body-frame angular velocity (the -omega IC; sets up free-spin/rocking)
  void setOmega0(double wx, double wy, double wz);

  // state — the Element hooks the DomainComponent route would NOT get for free.
  // P1: trivial (R=I). P2: commitState advances the body-frame SO(3) state and
  // re-forms the slave-MP transport blocks; update() refreshes follower kinematics.
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  // matrices / forces — a rigid body adds none of its own (see class note).
  const Matrix& getTangentStiff(void);
  const Matrix& getInitialStiff(void);
  const Matrix& getMass(void);
  const Vector& getResistingForce(void);
  const Vector& getResistingForceIncInertia(void);

  // no physical Rayleigh damping (a βK must not shrink dt_cr); getDamp /
  // getRayleighDampingForces return zero — same base lazy-damping-slot landmine
  // the RBE2 sibling documents (LadrunoKinematicCoupling, ADR 29 §6).
  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc);
  const Matrix& getDamp(void);
  const Vector& getRayleighDampingForces(void);

  // parallel (serial-only v1 — stubs, ADR 58 §6)
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  // recorder/output (ADR 58 D9): CoM {disp,vel,accel}, orientation (quaternion),
  // angular velocity/momentum. Most ride the internal CoM node directly.
  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

 private:
  int ndm;                 // spatial dimension (3 in v1)
  int nSlave;              // number of slave nodes N
  int nDOF;                // Σ_i ndf_i over the slaves (element DOF count), set at setDomain
  ID slaveTags;            // external node tags [ S_1 .. S_N ]
  Node** theSlaves;        // resolved slave pointers (size N)

  // private internal master node at the CoM (6 DOF: 3 translation + 3 rotation)
  int intNodeTag;          // <0 ⇒ auto-pick a free tag at setDomain
  Node* theIntNode;        // owned; added to / removed from the Domain here

  // slave-tying constraints (one per slave), owned + added/removed here
  MP_Constraint** theMPs;  // size N

  // cached Domain for teardown: the base DomainComponent::getDomain() is nulled on
  // detach (setDomain(0)), so the destructor cannot rely on it to remove the owned
  // internal node + MPs. Cache our own (Joint3D pattern); never nulled on detach.
  Domain* theDomainCache;

  // condensed body properties (ADR 58 D3)
  double mUser;            // user-supplied total mass (<0 ⇒ condense from slaves' nodal mass)
  double mBody;            // resolved total mass Σ m_i
  Vector xc;               // center of mass (size ndm), resolved at setDomain
  Matrix Ibody;            // 3×3 body-frame inertia about the CoM (floored; P2 uses it)

  // original slave nodal mass, captured so teardown can restore it (the body
  // ZEROES slave mass at setDomain to avoid double-counting once it is condensed
  // onto the internal CoM node — slaves become massless kinematic followers).
  Matrix** slaveMass0;     // size N (captured copies; 0 if a slave had no mass)

  // body-frame SO(3) rotational state (ADR 58 D2, P2 Stage 1). The body integrates
  // its own rotation in commitState OFF the global solve (D3): the primary state is
  // the SPATIAL angular momentum L (conserved exactly under zero torque ⇒ the
  // Dzhanibekov ‖L‖ gate is met by construction), plus the orientation quaternion.
  // Each step  ω_body = Ibody⁻¹·(qᵀ·L),  q ← q · exp(ω_body·dt),  L ← L + m·dt.
  OpenSees::Versor qCommit, qTrial;   // orientation (body→spatial); identity = {{0,0,0},1}
  Vector3D Lcommit, Ltrial;           // spatial angular momentum
  Vector3D L0;                        // initial spatial L (for revertToStart)
  Vector3D omega0;                    // initial body-frame angular velocity (-omega IC)
  bool haveOmega0;                    // -omega supplied
  Matrix IbodyInv;                    // 3×3 constant body-frame inertia inverse
  double lastTimeCommit, lastTimeTrial; // domain time of the last advance (commit/trial, for dt)
  bool started;                       // first commit seeds lastTime (no giant first-dt on a time jump)

  bool valid;              // setDomain succeeded (well-posed body)

  // zero returns (a rigid body assembles no stiffness/mass/force of its own)
  Matrix* K0;              // nDOF × nDOF, always zero
  Matrix* M0;              // nDOF × nDOF, always zero (mass is on the internal node)
  Vector* P0;              // nDOF, always zero
  Matrix* C0;              // damping, always zero (getDamp bypass)
  Vector* dampF;           // Rayleigh damping force, always zero

  // SO(3) helpers (P2)
  Vector3D gatherCoMTorque(void);      // net spatial moment about the CoM (Stage 1: 0 — free spin)
  Vector3D omegaBodyFrom(const OpenSees::Versor& q, const Vector3D& L) const;  // ω_body = Ibody⁻¹·(qᵀ·L)
  static bool invert3x3(const Matrix& A, Matrix& Ainv);  // analytic 3×3 inverse

  // lifecycle helpers (Joint3D pattern)
  int  resolveCentroidAndMass(void);   // CoM + m_body from slave coords/masses (or mUser); 0 ok, -1 fail
  int  createInternalNode(Domain* theDomain);   // build + addNode the CoM node
  int  buildSlaveConstraints(Domain* theDomain);// one rigid-link MP per slave
  void teardown(Domain* theDomain);    // remove + delete MPs and internal node (guarded)
  void allocateZeroReturns(void);      // size K0/M0/P0/C0/dampF to nDOF and zero
};

#endif
