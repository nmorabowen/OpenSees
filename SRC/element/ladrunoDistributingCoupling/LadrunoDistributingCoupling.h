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

// Ladruno: RBE3 / distributing-coupling element (ADR 28). See
// Ladruno_implementation/28_ladruno_distributing_coupling_rbe3_adr.md.
//
//   Distributing coupling: a REFERENCE node R (translations + rotations) whose
//   motion is the WEIGHTED LEAST-SQUARES rigid-body fit of a set of INDEPENDENT
//   nodes I_i (translations, weights w_i) — equivalently, a load applied at R is
//   DISTRIBUTED to the independents as a statically-equivalent force set, adding
//   NO stiffness to them. = Abaqus *COUPLING,distributing / LS-DYNA
//   *CONSTRAINED_INTERPOLATION / Nastran RBE3. Driver: the ndf-mismatch
//   moment-transfer interface (a 6-DOF beam/shell node into a face of 3-DOF
//   solid nodes — the moment enters the continuum as a distributed force couple,
//   the face stays free to deform). Contrast RBE2 (kinematic coupling): a
//   reference RIGIDLY drives a slave set (adds rigid stiffness).
//
//   External nodes: [ R, I_1, ..., I_N ].  Weighted centroid x_c = Σw_i x_i / W,
//   r_i = x_i − x_c, weighted position-inertia I_c = Σ w_i[(r_i·r_i)1 − r_i⊗r_i].
//   Constraint gaps (kinematically complete for an OFFSET reference node):
//       g_t = ( u_R + θ_R × (x_c − x_R) ) − (1/W) Σ w_i u_i        (size ndm)
//       g_r =   P θ_R − I_c⁺ Σ w_i (r_i × u_i)                      (size nrot)
//   with I_c⁺ the spectral pseudo-inverse (drop collinear/degenerate axes) and
//   P its range projector. Penalty: t_t = K_t g_t, t_r = K_r g_r; P = Bᵀt,
//   K = BᵀDB, D = diag(K_t I, K_r I). The force split f_i = (w_i/W)t_t +
//   w_i(α×r_i), α = I_c⁺ t_r, is the work-conjugate transpose — net force AND
//   moment balance exact for any K.
//
//   The TRANSLATION block reuses the shared LadrunoEmbeddedKernel (weight-agnostic);
//   the ROTATION + transport blocks are RBE3-specific (cross-DOF lever-arm coupling
//   the per-node kernel cannot assemble) and are built here.
//
//   classTag ELE_TAG_LadrunoDistributingCoupling = 33011 (33009/33010 left to the
//   VEM/SBFEM frontier reservations). Shares
//   SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.{h,cpp}.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoDistributingCoupling_h
#define LadrunoDistributingCoupling_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;

class LadrunoDistributingCoupling : public Element
{
 public:
  LadrunoDistributingCoupling(int tag, int ndm, int refNode, const ID& indepNodes,
                              const Vector& weights, double kt,
                              double kr, bool krUser, int enforce,
                              bool bipenalty, int bpMode, double bpDt, double bpBeta,
                              double kAlpha, int hostEleTag, bool ktAuto,
                              bool initGapCapture = true);
  LadrunoDistributingCoupling();
  ~LadrunoDistributingCoupling();

  const char* getClassType(void) const { return "LadrunoDistributingCoupling"; }

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
  // factors so a βK can't spuriously shrink dt_cr — ADR 28 §5).
  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc);
  // self-reported explicit critical step: min over the (translational, rotational)
  // penalty classes of 2√(m/k) (ADR 28 §5). −1 when bipenalty off / unbounded.
  double getExplicitCriticalTimeStep(void);

  // parallel
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

 private:
  int ndm;                  // 2 or 3
  int nrot;                 // rotation DOFs: 3 (3D) or 1 (2D, drilling)
  int nIndep;               // number of independent nodes N
  ID connectedNodes;        // [ refNode, indep_1 .. indep_N ]
  Vector weights;           // user weights w_i (size N)

  double Kt;                // translational penalty (resolved value)
  double Kr;                // rotational penalty (derived K_t·ℓ² or user -kr)
  bool krUser;              // -kr given numerically (else derive)
  bool ktAuto;              // -k auto: resolve K_t from a representative -host element
  double kAlpha;            // multiplier for the auto K_t
  int hostEleTag;           // representative host element (>=0) for -k auto / -wcap; -1 else
  bool ktResolved;          // transient: auto K_t / derived K_r resolved this run

  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  Vector lambda;            // AL translational multiplier (ndm), Uzawa-updated
  Vector lambda_r;          // AL rotational multiplier (nrot), Uzawa-updated

  // bipenalty (ADR 28 §5): m_p on the reference node's translational DOFs, I_p on
  // its rotation DOFs (the massless reference node is the explicit hazard).
  bool bipenalty;
  int bpMode;               // 0 = -dtcr budget, 1 = -wcap β (needs -host)
  double bpDt, bpBeta;
  double mPenalty, iPenalty;
  bool bpResolved;
  bool bpWarned;

  // geometry + RBE3 operators, resolved ONCE at setDomain (coords known there):
  bool valid;               // setDomain succeeded (geometry well-posed)
  double W;                 // Σ w_i
  double ell2;              // Σ w_i|r_i|²/W — the rotation stiffness scale (K_r = K_t·ℓ²)
  Vector dRef;              // x_c − x_R (the transport lever, size ndm)
  Matrix rvec;              // r_i = x_i − x_c (nIndep × ndm)
  Matrix Icplus;            // spectral pseudo-inverse of I_c (nrot × nrot)
  Matrix Pproj;             // range projector onto the KEPT rotation subspace (nrot × nrot)
  int nKept;                // resolvable reference-rotation axes (nKept ≤ nrot)

  // per-node DOF bookkeeping (the reference may be 6-DOF, independents 3-DOF):
  int nDOF;                 // Σ_i ndf_i
  ID nodeNdf;               // ndf of each external node (size 1+N)
  ID dofOffset;             // element-DOF offset of each node (size 1+N)

  Matrix* B;                // constant gap operator (ndm+nrot) × nDOF (built at setDomain)

  // initial-gap (offset) capture for stress-free staged activation (ADR 28; mirrors
  // ASDEmbeddedNodeElement m_U0 / LadrunoEmbeddedNode g0). g_r captured through the
  // SAME projected operators ⇒ gr0 lives in the kept subspace automatically.
  bool initGapCapture;      // -absolute ⇒ false (keep the absolute tie)
  bool g0Computed;          // capture-once guard
  Vector g0;                // captured translational offset (ndm)
  Vector gr0;               // captured rotational offset (nrot)

  Node** theNodes;          // size 1 + N
  Matrix* K;                // nDOF × nDOF
  Vector* P;                // nDOF
  Matrix* M0;               // mass (nDOF × nDOF; zero unless bipenalty)

  void allocate(void);
  void resolveGeometry(void);    // x_c, r_i, dRef, W, ℓ², I_c eigen → I_c⁺/P, nKept
  void buildB(void);             // fill *B from the resolved operators
  void resolveAutoKt(void);      // -k auto (needs -host) + derive K_r = K_t·ℓ²
  void resolveBipenalty(void);
  double effectiveCouplingStiffness(void) const { return Kt > 0.0 ? Kt : 0.0; }
  void computeGap(Vector& g);    // full gap (ndm+nrot), incl −g0/−gr0 if captured
  void captureInitialGap(void);

  // cross operator C_i (nrot × ndm): (C_i u)_r = (r_i × u)_r.
  double crossOp(int i, int r, int j) const;
  // transport operator T (ndm × nrot): (T θ)_k = (θ × (x_c−x_R))_k.
  double transOp(int k, int r) const;
};

#endif
