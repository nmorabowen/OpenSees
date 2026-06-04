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

// Ladruno: general node-to-host COUPLING element (ADR 23). Phase 1 = the isotropic
// translational tie (U).
//
//   Ties ONE constrained node's TRANSLATIONS to a host element's nodes through
//   precomputed shape-function weights N_i, so an arbitrary node can be embedded in
//   a NON-matching host mesh (mesh stitching, point embedding, SSI). It is the
//   isotropic sibling of LadrunoEmbeddedRebar: same kinematic tie + penalty/AL +
//   `-kt auto` + bipenalty machinery (shared via LadrunoEmbeddedKernel), but with an
//   ISOTROPIC penalty (no bar axis, no axial/transverse split, no bond-slip).
//
//   External nodes: [ cNode, hostNode_1, ..., hostNode_M ].
//   Gap (translational, current-config relative disp of the constrained point):
//       g = u_c − Σ_i N_i u_host,i        (ndm-vector, the first ndm DOFs of each node)
//   Isotropic coupling traction and tangent:
//       t = K_u · g  (+ λ for -enforce al),   D_u = K_u · I_{ndm}.
//   Resisting force  P = Bᵀt,  tangent  K = Bᵀ D_u B,  B = [ I | −N_1 I … −N_M I ]
//   (the shared LadrunoEmbeddedKernel assembles these in a compact ndm-packed space;
//   the element SCATTERS them into the full per-node-ndf layout via DOF offsets, so
//   the constrained / host nodes may carry more DOFs than ndm — only the first ndm
//   translational DOFs are coupled in Phase 1; rotation (UR) and pressure (UP) are
//   ADR 23 Phase 2 / Phase 1b).
//
//   The isotropic tie is ALREADY frame-objective under rigid host rotation (g and
//   N_i(ξ) both transform with the host), so — unlike the anisotropic rebar — there
//   is NO `-corot` (LEDGER_quirk). Penalty adds STIFFNESS ONLY ⇒ implicit + explicit.
//
//   classTag ELE_TAG_LadrunoEmbeddedNode = 33006 (Ladruno private element band).
//   Shares SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.{h,cpp}.
//   See Ladruno_implementation/23_ladruno_embedded_node_adr.md.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoEmbeddedNode_h
#define LadrunoEmbeddedNode_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;
class UniaxialMaterial;

class LadrunoEmbeddedNode : public Element
{
 public:
  LadrunoEmbeddedNode(int tag, int ndm, int cNode, const ID& hostNodes,
                      const Vector& shape, double ku,
                      int hostEleTag = -1, bool ktAuto = false,
                      double ktAlpha = 0.0, int enforce = 0,
                      bool bipenalty = false, int bpMode = 0,
                      double bpDt = 0.0, double bpBeta = 0.0,
                      bool pressure = false, double kp = 0.0,
                      bool rot = false, double kr = 0.0, bool krAuto = false,
                      double krAlpha = 0.0, const Matrix* gradN = 0,
                      bool matMode = false, const Vector* normalDir = 0,
                      const Vector* orientDir = 0,
                      UniaxialMaterial* matN = 0, UniaxialMaterial* matT1 = 0,
                      UniaxialMaterial* matT2 = 0);
  LadrunoEmbeddedNode();
  ~LadrunoEmbeddedNode();

  const char* getClassType(void) const { return "LadrunoEmbeddedNode"; }

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

  // ADR 23 D5 / ADR 20 §10.6 (D-bp-5) — a pure penalty COUPLING carries NO physical
  // Rayleigh damping; refuse the factors so a βK can't spuriously shrink dt_cr.
  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc);

  // ADR 23 D5 / ADR 20 §10.6.1 — self-reported explicit critical step 2√(m_p/k_eff)
  // (the per-element eigensolve sees this coupling as λ_max=0). −1 when bipenalty off.
  double getExplicitCriticalTimeStep(void);

  // parallel
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

 private:
  int ndm;                  // 2 or 3
  int nHost;                // number of host nodes M
  ID connectedNodes;        // [cNode, host_1..host_M]
  Vector Nshape;            // host shape-function weights N_i (size M)

  double Ku;                // isotropic translational penalty (resolved value)

  // ADR 23 Phase 1b — pressure tie (UP). Opt-in via -pressure (pflag). Couples the
  // constrained node's pressure DOF (index ndm, the u-p convention) to the host's
  // interpolated pressure:  g_p = p_c − Σ N_i p_host,i;  t_p = K_p·g_p (+ λ_p for AL).
  // Activated in setDomain only if all coupled nodes are u-p (ndf ≥ ndm+1); else the
  // element degrades to U-only (warn). The pressure mode is NOT bipenalty-bounded —
  // the pressure tie is implicit-recommended (ADR 23 D5 / M1 "pressure → implicit").
  double Kp;                // pressure penalty (numeric; -kp). Auto-scale deferred.
  bool upActive;            // resolved in setDomain: pflag && all nodes u-p
  double lambda_p;          // AL pressure multiplier (scalar); committed, Uzawa-updated
  double computeGapP(void); // g_p = p_c − Σ N_i p_host,i (pressure DOF = index ndm)

  // ADR 23 Phase 2 — rotation tie (UR). Opt-in via -rot (rflag). Ties the constrained
  // node's ROTATION DOFs to the host's CONTINUUM rotation at the embedded point,
  // θ_host = ½ curl(u) = skew(∇u)|_ξ, read from the host's cartesian shape gradients
  // ∂N_i/∂x (gradN, supplied by -dNdx or queried via Element::getInterpolationGradients).
  //   3D (nrot=3): θ_host = ½ Σ_i (∇N_i × u_host,i), rotation DOFs at indices ndm..ndm+2.
  //   2D (nrot=1): θ_z   = ½ Σ_i (∂N_i/∂x·u_y − ∂N_i/∂y·u_x), drilling DOF at index ndm.
  // gap g_r = θ_c − θ_host; traction t_r = K_r·g_r (+ λ_r for AL); D_r = K_r·I_{nrot}.
  // Activated in setDomain only if -rot was set AND the cNode carries the rotation DOFs
  // (ndf ≥ ndm+nrot); else degrade to no-UR (warn). UR is APPROXIMATE / mesh-limited:
  // on a CST/TET4 host ∂N/∂x is element-constant ⇒ a single rigid-spin tie; moment-
  // critical embedments need a higher-order host (BezierTet10) — ADR 23 §3 UR-4.
  int rflag;                // 0 = no UR; 1 = -rot requested
  int nrot;                 // rotation DOFs coupled: 3 (3D) or 1 (2D); set from ndm
  bool rActive;             // resolved in setDomain: rflag && cNode has rotation DOFs
  Matrix gradN;             // host cartesian gradients ∂N_i/∂x_j (nHost × ndm), for UR
  double Kr;                // rotational penalty (resolved value)
  bool krAuto;              // -kr auto: K_r = krAlpha·K_u·lch² (needs -host for lch)
  double krAlpha;           // multiplier for the auto K_r (ADR 23 D3)
  bool krResolved;          // transient: auto K_r already resolved this run
  Vector lambda_r;          // AL rotation multiplier (size nrot); Uzawa-updated
  void resolveKr(void);     // resolve auto K_r (after K_u; queries host lch)
  void computeGapR(Vector& gr);  // g_r = θ_c − θ_host(ξ) (size nrot)
  // host-rotation operator: row r, host-node i, local trans dof j → ∂θ_host_r/∂u_i,j.
  double rotOper(int r, int i, int j) const;

  // ADR 23 Phase 2b (D9) — material-driven interface mode (node ↔ interpolated host
  // point). Opt-in by supplying any of -matN/-matT1/-matT2: the ISOTROPIC translational
  // block D_u=K_u·I is replaced by a local frame {e_0=normal, e_1,e_2=tangents} where
  // each direction carries a UNIAXIAL material (penalty fallback K_u where no material).
  // Decompose the gap g into local comps g_d = g·e_d; t = Σ_d f_d(g_d)·e_d, D = Σ_d
  // k_d·e_d⊗e_d with (f_d,k_d) = mat_d stress/tangent in FORCE units (D9-4 — the gap is
  // in metres and the material returns force, so NO bondScale). Models: cohesive
  // (softening uniaxial), unilateral gap (ENT/ElasticPPGap on the normal), elastic
  // bedding (Elastic), bond. Coulomb friction is only APPROXIMATE (uncoupled per-
  // direction — a fixed ElasticPP slip force, not μ·N; rigorous → LadrunoContact).
  // v1 uses the REFERENCE frame (-normal/-orient); the -corot frame co-rotation (D9.1,
  // for a large-rotation contact normal) is DEFERRED to v2. AL augments only the PENALTY
  // directions and re-projects λ off the material directions each commit (M4/D9-3).
  int matMode;                 // 0 = isotropic K_u·I; 1 = material-frame interface
  Vector normalDir;            // -normal: unit frame normal e_0 (size ndm)
  Vector orientDir;            // -orient: hint for the first tangent e_1 (size ndm; opt)
  bool haveOrient;             // -orient supplied
  Matrix frame;                // ndm×ndm orthonormal local frame (cols e_0..e_{ndm-1})
  UniaxialMaterial* matDir[3]; // per-direction materials (null ⇒ penalty K_u); [0]=normal
  bool hasMat[3];              // direction d is material-driven (true) vs penalty (false)
  void buildFrame(void);       // construct the orthonormal frame from normalDir/orientDir
  // translational traction+tangent: isotropic (matMode 0) or material-frame (matMode 1).
  void formTransTraction(const Vector& g, Vector& t, Matrix& D);

  // ADR 23 D4 — constraint-enforcement strategy. 0 = penalty (default), 1 = AL.
  // AL adds a per-element multiplier λ (translational, size ndm) with the SAME
  // tangent K = BᵀD_uB; per-step Uzawa λ += K_u·g in commitState. The isotropic tie
  // has NO preferred axis, so (unlike the rebar) there is NO directional re-projection
  // of λ (ADR 23 D4 / ES-3).
  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  Vector lambda;            // AL multiplier (size ndm); committed, Uzawa-updated

  // ADR 23 D3 / M2 — auto-scaled translational penalty from the host's initial
  // stiffness:  K_u = ktAlpha * max_i|K_host(i,i)|  (getInitialStiff max-abs-diagonal;
  // NO lch factor — the diagonal already carries ~E·lch). Lazy (first assembly).
  int hostEleTag;           // -host element tag (>=0) for auto-Ku / wcap; -1 if none
  bool ktAuto;              // -k auto: resolve K_u from the host stiffness
  double ktAlpha;           // dimensionless multiplier for the auto K_u
  bool ktResolved;          // transient: auto K_u already resolved this run
  void resolveAutoKt(void);

  // ADR 23 D5 / ADR 20 §10.6 — bipenalty. For Phase 1 (U only) the coupling is
  // dimensionally homogeneous (translational), so the scalar k_eff = K_u with the
  // mass penalty m_p lumped on the cNode's first ndm DOFs is correct (the per-DOF-class
  // (stiffness,inertia) generalization of M1/ES-1 arrives with UR/UP in Phase 1b/2).
  bool bipenalty;           // -bipenalty: add the mass penalty m_p
  int bpMode;               // 0 = -dtcr budget, 1 = -wcap β (auto ω_host)
  double bpDt;              // -dtcr target step (bpMode 0)
  double bpBeta;            // -wcap penalty-frequency ratio β (bpMode 1)
  double mPenalty;          // resolved mass penalty m_p (per cNode translational DOF)
  // ADR 23 D5 / M1·ES-1 — per-DOF-class (stiffness, inertia) pair for the ROTATION
  // mode: a translational-only m_p cannot bound the rotation mode (different units).
  // I_p is lumped on the cNode rotation DOFs (same -dtcr/-wcap budget formula but
  // with K_r), so dt_r = 2√(I_p/K_r) = dt_u — the rotation mode self-bounds.
  double iPenalty;          // resolved rotational mass penalty I_p (per cNode rot DOF)
  bool bpResolved;          // transient: m_p / I_p already resolved this run
  double effectiveCouplingStiffness(void); // k_eff = K_u (translational class)
  void resolveBipenalty(void);

  // per-node DOF bookkeeping (resolved in setDomain): the constrained / host nodes
  // may carry more DOFs than ndm (e.g. u-p / u-r), so the compact ndm-packed coupling
  // is SCATTERED into the full layout via these node offsets.
  int nDOF;                 // Σ_i ndf_i (0 until setDomain)
  ID nodeNdf;               // ndf of each external node (size 1+M)
  ID dofOffset;             // global element-DOF offset of each node (size 1+M)
  void allocate(void);      // (re)allocate K/P/M0 to nDOF (after setDomain)

  // DOF-class request flag: 0 = U only, 1 = U+P requested (-pressure). The UR (rotation)
  // bit is reserved for Phase 2.
  int pflag;                // 0 = U; 1 = -pressure requested (UP, Phase 1b)

  Node** theNodes;          // size 1 + M
  Matrix* K;                // nDOF x nDOF
  Vector* P;                // nDOF
  Matrix* M0;               // mass (nDOF x nDOF; zero unless bipenalty)

  void computeGap(Vector& g);  // g = u_c − Σ N_i u_host,i (translational, size ndm)
};

#endif
