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

// Ladruno: RBE2 / kinematic-coupling element (ADR 29). See
// Ladruno_implementation/29_ladruno_kinematic_coupling_rbe2_adr.md.
//
//   Kinematic coupling: a REFERENCE node R (translations + rotations) RIGIDLY
//   drives a set of SLAVE nodes S_i — each slave follows R's rigid-body motion
//   u_i = u_R + θ_R × d_i  (and θ_i = θ_R where the slave carries rotations and the
//   rotation is tied), d_i = x_i − x_R. A load/displacement at R is transmitted to
//   the set WITH rigid stiffness. = Abaqus *COUPLING,kinematic / LS-DYNA
//   *CONSTRAINED_NODAL_RIGID_BODY / Nastran RBE2. Contrast RBE3
//   (LadrunoDistributingCoupling): R is the WEIGHTED AVERAGE of the set, adding no
//   stiffness. In RBE2 the single node R is the MASTER; in RBE3 it is the dependent.
//
//   External nodes: [ R, S_1, ..., S_N ].  Constraint gaps (per slave, per tied DOF):
//       g_i^t = u_i − ( u_R + θ_R × d_i )        (translation rows; transport block)
//       g_i^r = θ_i − θ_R                         (rotation rows; identity blocks)
//   Penalty: t = D g, P = Bᵀt, K = Σ Bᵢᵀ D_i Bᵢ, D_i = diag(K_t I, K_r I) over the
//   tied DOFs. The transport block ∂g_i^t/∂θ_R = +[d_i]_× is a SIGN FLIP of RBE3's
//   transOp (RBE2 carries −θ_R×d_i). Self-equilibrated for any K.
//
//   -dof selects the dependent components on each slave (1..ndm translations,
//   ndm+1..ndm+nrot rotations; default = every DOF the slave possesses). The gap
//   vector is RAGGED (per-slave tied count varies for mixed 3-/6-DOF sets) — the
//   layout (gapNode/gapDof/gapIsRot) is resolved ONCE at setDomain.
//
//   Explicit safety (ADR 29 §6): R is the MASTER and is frequently massed (platen /
//   footing / equipment), so -bipenalty DEFAULTS OFF and, when on, lumps a penalty
//   mass only on tied DOFs that are ACTUALLY massless (scan over R AND every slave —
//   a massless slave is the RBE2-specific hazard), sized from the Gershgorin row-sum
//   of the assembled penalty tangent (≥ λ_max ⇒ conservative). getDamp /
//   getRayleighDampingForces are overridden to return zero (the no-op
//   setRayleighDampingFactors never allocates the base Element's lazy damping slot,
//   so the base path would dereference theMatrices[-1] when a TRANSIENT integrator
//   forms the C-tangent → crash; D ≡ 0 is physically correct).
//
//   classTag ELE_TAG_LadrunoKinematicCoupling = 33012 (next free after RBE3=33011;
//   33009/33010 reserved VEM/SBFEM). Sibling of LadrunoDistributingCoupling.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoKinematicCoupling_h
#define LadrunoKinematicCoupling_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;

class LadrunoKinematicCoupling : public Element
{
 public:
  LadrunoKinematicCoupling(int tag, int ndm, int refNode, const ID& slaveNodes,
                           const ID& dofSel, double kt, double kr, bool krUser,
                           int enforce, bool bipenalty, int bpMode, double bpDt,
                           double bpBeta, double kAlpha, int hostEleTag, bool ktAuto,
                           bool initGapCapture = true);
  LadrunoKinematicCoupling();
  ~LadrunoKinematicCoupling();

  const char* getClassType(void) const { return "LadrunoKinematicCoupling"; }

  // domain
  int getNumExternalNodes(void) const;
  const ID& getExternalNodes(void);
  Node** getNodePtrs(void);
  int getNumDOF(void);
  void setDomain(Domain* theDomain);

  // state
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  // matrices / forces
  const Matrix& getTangentStiff(void);
  const Matrix& getInitialStiff(void);
  const Matrix& getMass(void);
  const Vector& getResistingForce(void);
  const Vector& getResistingForceIncInertia(void);

  // a pure penalty coupling carries no physical Rayleigh damping (refuse the
  // factors so a βK can't spuriously shrink dt_cr). getDamp / getRayleighDampingForces
  // are ALSO overridden to return zero — see the class note + ADR 29 §6 (the base
  // lazy-damping-slot index-landmine that crashes a TRANSIENT run).
  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc);
  const Matrix& getDamp(void);
  const Vector& getRayleighDampingForces(void);
  // self-reported explicit critical step: min over lumped (massless) DOFs of
  // 2√(m_p/k_dof). −1 when bipenalty off / no massless DOF.
  double getExplicitCriticalTimeStep(void);

  // parallel
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

 private:
  int ndm;                  // 2 or 3
  int nrot;                 // rotation DOFs: 3 (3D) or 1 (2D drilling)
  int nSlave;               // number of slave nodes N
  ID connectedNodes;        // [ refNode, slave_1 .. slave_N ]
  ID dofSel;                // requested 1-based dependent components (size 0 = default all)

  double Kt;                // translational penalty (resolved value)
  double Kr;                // rotational-tie penalty (derived K_t·ℓ² or user -kr)
  bool krUser;              // -kr given numerically (else derive)
  bool ktAuto;              // -k auto: resolve K_t from a representative -host element
  double kAlpha;            // multiplier for the auto K_t
  int hostEleTag;           // representative host element (>=0) for -k auto / -wcap; -1 else
  bool ktResolved;          // transient: auto K_t / derived K_r resolved this run
  double ell2;              // floored rotation length scale; default K_r = K_t·ℓ² (ADR 29 §6 D1)

  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  Vector lambdaAL;          // per-gap-row AL multiplier (size nGap), Uzawa-updated

  // bipenalty (ADR 29 §6): default OFF; when on, lump a penalty mass on every tied DOF
  // that is ACTUALLY massless (R AND slaves), sized from the Gershgorin row-sum of K.
  bool bipenalty;
  int bpMode;               // 0 = -dtcr budget, 1 = -wcap β (needs -host)
  double bpDt, bpBeta;
  bool bpResolved;

  bool hasRefRot;           // R carries rotation DOFs (ndf_R >= ndm+nrot view)

  // geometry + ragged gap layout, resolved ONCE at setDomain (coords known there):
  bool valid;               // setDomain succeeded (well-posed)
  Matrix dvec;              // d_i = x_i − x_R (nSlave × ndm)
  int nGap;                 // total tied DOFs (ragged sum over slaves)
  ID gapNode;               // element node slot (1..N) of the dependent slave for each row
  ID gapDof;                // node-local DOF index (component − 1) for each row
  ID gapIsRot;              // 0 = translation row, 1 = rotation row

  // per-node DOF bookkeeping (R may be 6-DOF, slaves 3- or 6-DOF):
  int nDOF;                 // Σ_i ndf_i
  ID nodeNdf;               // ndf of each external node (size 1+N)
  ID dofOffset;             // element-DOF offset of each node (size 1+N)

  Matrix* B;                // constant gap operator nGap × nDOF (built at setDomain)

  // initial-gap (offset) capture for stress-free staged activation. g0 is per-row.
  bool initGapCapture;      // -absolute ⇒ false
  bool g0Computed;
  Vector g0;                // captured per-row offset (size nGap)

  Node** theNodes;          // size 1 + N
  Matrix* K;                // nDOF × nDOF
  Vector* P;                // nDOF
  Matrix* M0;               // diagonal lumped bipenalty mass (nDOF × nDOF; zero unless bipenalty)
  Matrix* C0;               // damping (nDOF × nDOF; ALWAYS zero — getDamp bypass)
  Vector* dampF;            // Rayleigh damping force (nDOF; ALWAYS zero)

  void allocate(void);
  void resolveGeometry(void);    // d_i, ragged layout, ℓ², refuse checks
  void buildB(void);             // fill *B from the layout + per-slave d_i
  void resolveAutoKt(void);      // -k auto (needs -host) + derive K_r = K_t·ℓ²
  void resolveBipenalty(void);   // Gershgorin per-DOF massless-scan lumping into *M0
  double effectiveCouplingStiffness(void) const { return Kt > 0.0 ? Kt : 0.0; }
  void computeGap(Vector& g);    // full gap (nGap), incl −g0 if captured
  void captureInitialGap(void);
  double rowPenalty(int row) const { return gapIsRot(row) ? Kr : Kt; }

  // transport operator: (T_i θ)_k = (θ × d_i)_k ; returns ∂(θ×d_i)_k/∂θ_r. The g_i^t
  // B-block is the NEGATIVE of this (sign flip vs RBE3 — ADR 29 §2.4 M1).
  double transOp(int i, int k, int r) const;
};

#endif
