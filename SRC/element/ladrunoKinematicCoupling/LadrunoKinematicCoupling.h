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
//   of the assembled penalty tangent (≥ λ_max ⇒ conservative). D ≡ 0: the Rayleigh
//   factors are refused and getDamp is zero, inherited from LadrunoUndampedElement
//   (WP-123; that header explains the transient-crash history).
//
//   classTag ELE_TAG_LadrunoKinematicCoupling = 33012 (next free after RBE3=33011;
//   33009/33010 reserved VEM/SBFEM). Sibling of LadrunoDistributingCoupling.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoKinematicCoupling_h
#define LadrunoKinematicCoupling_h

#include <LadrunoUndampedElement.h>   // WP-123: refuses Rayleigh, zero getDamp
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;

class LadrunoKinematicCoupling : public LadrunoUndampedElement
{
 public:
  LadrunoKinematicCoupling(int tag, int ndm, int refNode, const ID& slaveNodes,
                           const ID& dofSel, double kt, double kr, bool krUser,
                           int enforce, bool bipenalty, int bpMode, double bpDt,
                           double bpBeta, double kAlpha, int hostEleTag, bool ktAuto,
                           bool initGapCapture = true, int alUpdate = 0);
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

  // setRayleighDampingFactors (refused) / getDamp (zero): LadrunoUndampedElement (WP-123).
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
  // Ladruno (WP-101 r1): resolveAutoKt() must NOT latch ktResolved when the -host element is
  // not in the domain yet. Domain::addElement() calls update() at DECLARATION time, so a deck
  // that declares the coupling BEFORE its host used to latch on that failed lookup: -k auto
  // silently stayed at 1e12 and the conditioning warning never fired. Retry instead; warn once
  // if the host is still missing after an analysis exists.
  bool ktHostMissWarned;    // transient
  double ell2;              // floored rotation length scale; default K_r = K_t·ℓ² (ADR 29 §6 D1)

  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  // Ladruno (WP-101 / ADR 29 §4.2b): WHERE the Uzawa recursion advances.
  //
  //   alUpdate = 0 ("commit", DEFAULT) — one update in commitState: λ ← λ + D g once per
  //     COMMITTED step, a first-order Uzawa ACROSS steps. The residual is then a genuine
  //     function of u within a step (λ is frozen), which is what every algorithm's Jacobian
  //     model assumes. This is the ONLY cadence that is safe with every algorithm and
  //     integrator, and it is the cadence the ADR-41 D1 held-load augmentation sweep
  //     (`ladrunoBeginAugment` / `LoadControl 0.0` / `ladrunoEndAugment`) turns into a proper
  //     OUTER Uzawa loop: each held-load analyze is an inner solve at FIXED λ, and the
  //     Domain::commit() at its end does the outer update. Measured: the rigidity gate closes
  //     to ≤1e-9 in 4-5 passes for all of Newton / ModifiedNewton / KrylovNewton / BFGS /
  //     Broyden × LoadControl / DisplacementControl. THIS is the supported within-step route.
  //
  //   alUpdate = 1 ("iter", OPT-IN, EXPERT) — λ ← λ + D g inside update(), once per
  //     equilibrium iteration on the current trial displacements. The fixed point still has
  //     g ≡ 0 exactly (Δu = 0 ⇒ Δλ = 0 ⇒ D g = 0) and full Newton under LoadControl reaches
  //     it in one step, but the scheme is NOT generally safe: update() advances λ BEFORE the
  //     force is formed, so the tie force carries λ_k + 2·D·g(u_k) against a tangent that
  //     linearises a single D·g, and λ_k is path-dependent — the residual is NOT a function
  //     of u. Every secant / accelerated / re-solving method is then fed (du, dr) pairs that
  //     describe no Jacobian. MEASURED FAILURES on a linear 2×2×2 elastic gate:
  //       DisplacementControl  — 5/5 steps fail with EVERY algorithm (it re-solves dLambda
  //                              each iterate against a residual that moves independently of u)
  //       KrylovNewton / BFGS / Broyden + LoadControl — fail or diverge (Broyden to 2.5e275)
  //       ModifiedNewton       — survives on this linear model, but 10/10 fail on a
  //                              LadrunoBrick bbar + LadrunoJ2 host
  //     So `iter` is REFUSED at the first update() unless the active algorithm is full Newton
  //     AND the active static integrator is LoadControl (checked via OPS_GetAlgorithm /
  //     OPS_GetStaticIntegrator). Full Newton survives only because its contraction here is
  //     ~0.008. Prefer the augment sweep above.
  //
  // Either way λ must be restorable: revertToLastCommit rolls back to lambdaCommitted, or a
  // failed/retried step would inherit the multipliers of a discarded trial state.
  int alUpdate;             // 0 = per-commit (DEFAULT), 1 = per-iteration (opt-in, guarded)
  bool alGuardWarned;       // the `iter` refusal has been printed once (transient)
  bool alWarnedTransient;   // the AL-under-transient note has been printed once (transient)
  Vector lambdaAL;          // per-gap-row AL multiplier (size nGap), Uzawa-updated
  Vector lambdaCommitted;   // λ as of the last commitState (size nGap) — revert target
  // Domain::revertToLastCommit() and Domain::revertToStart() BOTH end with
  // `return this->update();` (Domain.cpp), so every element's update() is invoked once
  // more on the just-reverted state. An element that mutates state in update() — which is
  // exactly what the per-iteration Uzawa does — would therefore advance λ on the state it
  // just rolled back to, ratcheting a little further on every failed step. One-shot latch:
  // the revert arms it, the induced update() consumes it. Transient (not serialized).
  bool alSkipUpdate;

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

  void allocate(void);
  void resolveGeometry(void);    // d_i, ragged layout, ℓ², refuse checks
  void buildB(void);             // fill *B from the layout + per-slave d_i
  void resolveAutoKt(void);      // -k auto (needs -host) + derive K_r = K_t·ℓ² + conditioning warn
  // Ladruno (WP-101 r1): true ⇒ the active algorithm/integrator cannot carry -alUpdate iter,
  // and update() must abort. Reads OPS_GetAlgorithm / OPS_GetStaticIntegrator /
  // OPS_GetTransientIntegrator; silent (returns false) when no analysis is configured yet.
  bool refuseIterCadence(void);
  void warnAlUnderTransient(void);  // one-time note: -enforce al is unusable under transient
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
