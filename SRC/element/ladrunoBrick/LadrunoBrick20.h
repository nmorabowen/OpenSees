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

// LadrunoBrick20 — 20-node serendipity quadratic hexahedron (Ladruno fork).
//
// One element class, one classTag (ELE_TAG_LadrunoBrick20 = 33018), 60 DOF,
// one NDMaterial per Gauss point. The 3D sibling of the ADR-70 quadratic
// triangle (LadrunoLST) and the second-order member of the fork solid family
// (H8 = LadrunoBrick, quadratic Bezier simplex = BezierTet10). ADR 72.
//
//   element('LadrunoBrick20', tag, n1..n20, matTag,
//           ['-formulation', <std|uri>]     # default std (v1: std only live)
//           ['-lumped']                      # HRZ lumped mass (lands P3)
//           ['-b', bx, by, bz]               # body force
//           ['-damp', dampTag])              # Damping objects, per GP
//
// v1 (P1) scope — small strain, geometrically linear:
//   std  — full 3x3x3 (27-pt) Gauss displacement element. Reduce-to anchor:
//          reproduces upstream Twenty_Node_Brick stiffness / resisting force /
//          consistent mass to ~1e-12 on a distorted hex.
//   uri  — uniform 2x2x2 (8-pt) reduced integration (the C3D20R analog). The
//          enum ordinal (URI=1) is DECLARED and SERIALIZED now (append-only
//          forever), but the formulation is parser-REJECTED until P2.
//
// Deliberately NOT here (vs LadrunoBrick): no geom seam (linear only, no
// SolidTransformation), no hourglass/ssp/eas machinery, no bulk viscosity.
// -hourglass is a HARD parser error on this element by design (H20@2x2x2 has 6
// non-communicable modes; no control exists here — ADR 72 §2.2/§3.2).
//
// All shape / GP / B / mass / volume math comes from the pure P0 kernel
// Ladruno::hex20::* (LadrunoHex20Shape.h); materialPointers[L] pairs with
// Ladruno::hex20::GP27[L] (brcshl order) so the recorder rule + the upstream
// per-GP comparison line up.  // Ladruno

#ifndef LadrunoBrick20_h
#define LadrunoBrick20_h

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Damping.h>

class Response;   // Ladruno — fwd-decl for setResponse / the U1 advisory probe

class LadrunoBrick20 : public Element {

 public:

  // Formulation selector. Ordinals are SERIALIZED (packed into idData by
  // sendSelf) — append-only forever: STD keeps ordinal 0, URI keeps ordinal 1.
  // v1 ships STD only; URI is parser-rejected until P2 but the ordinal is
  // reserved now so a future stream loads back as the same element.  // Ladruno
  enum class Formulation { STD, URI };

  // null constructor (broker)
  LadrunoBrick20();

  // full constructor
  LadrunoBrick20(int tag,
                 int node1,  int node2,  int node3,  int node4,
                 int node5,  int node6,  int node7,  int node8,
                 int node9,  int node10, int node11, int node12,
                 int node13, int node14, int node15, int node16,
                 int node17, int node18, int node19, int node20,
                 NDMaterial &theMaterial,
                 Formulation formulation = Formulation::STD,
                 double b1 = 0.0, double b2 = 0.0, double b3 = 0.0,
                 int massType = 0,
                 Damping *theDamping = 0);

  virtual ~LadrunoBrick20();

  const char *getClassType(void) const { return "LadrunoBrick20"; };

  // domain
  void setDomain(Domain *theDomain);
  int  setDamping(Domain *theDomain, Damping *theDamping);

  // topology
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);

  // Reference-config element size for crack-band / regularized-softening
  // materials: lch = cbrt(V), V by the 27-pt rule (the LadrunoBrick / hex
  // convention). Degenerate V<=0 falls back to the Element base.  // Ladruno
  double getCharacteristicLength(void);

  // state
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  void Print(OPS_Stream &s, int flag);

  // matrices / vectors
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // parallel / database
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  // output
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  // plotting
  int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);

 private:

  // -------- sizes (single source of truth; loop bounds derive from these) --
  static const int NEN  = 20;   // nodes
  static const int NGP  = 27;   // Gauss points (std, full 27-pt)
  static const int NDF  = 3;    // dofs per node
  static const int NDOF = 60;   // element dofs
  static const int NSTR = 6;    // stress/strain components (Voigt)

  // -------- private attributes --------
  ID connectedExternalNodes;          // 20 node numbers
  Node *nodePointers[NEN];            // pointers to 20 nodes
  NDMaterial *materialPointers[NGP];  // pointers to 27 materials (one per Gauss point)

  Formulation formulation;            // the -formulation selector
  double b[3];                        // body forces
  double appliedB[3];                 // body forces applied with load
  int applyLoad;

  Vector *load;
  Matrix *Ki;
  Matrix *M0;                         // cached consistent mass (integrated once)
  int massType;                       // 0 consistent, 1 lumped (HRZ; lands P3)

  Damping *theDamping[NGP];

  // -------- geometry cache (reference-config; ~13 KB per instance) --------
  // Filled once at setDomain (nodal coords are reference-constant), consumed by
  // every hot path (update / residual / stiffness / mass / volume) — no kernel
  // shape/jacobian/invert calls in the Gauss loops. The cache builder validates
  // detJ > 0 at all 27 GPs (upstream Twenty_Node_Brick::Jacobian3d analog); a
  // non-positive detJ marks the element bad: update() returns -1 and the
  // tangent/mass/residual paths emit-once + zero, so the analysis fails loudly
  // instead of silently integrating a sign-flipped dv.  // Ladruno
  bool geomCached;                    // cache valid (rebuilt on every setDomain)
  bool badGeom;                       // dead element: non-positive detJ (or failed
                                      // material clone — same kill switch)
  bool warnedBadUse;                  // emit-once flag for use-after-bad paths
  bool hasMass;                       // any GP rho != 0 (cached at setDomain)
  double cachedN[NGP][NEN];           // shape functions at the 27 GPs
  double cachedDNdx[NGP][NEN][3];     // cartesian gradients at the 27 GPs
  double cachedDV[NGP];               // w_L * detJ_L

  // U1 advisory (ADR 72 §3.7, RESOLVED BOTH): softening/crack-band materials on
  // quadratic elements are theory-shaky (lch = cbrt(V) mis-sizes a sub-cell
  // band). At setDomain, probe materialPointers[0] for a dual-scalar "damage"
  // response (the cached-Response pattern from LadrunoBrick::damageResponse);
  // if present, emit a PROCESS-once opserr advisory and proceed (a 100k-element
  // mesh prints it once, not 100k times). Never fires for materials with no
  // dual-scalar damage channel. Not serialized.  // Ladruno
  static bool advisedDamage;

  // -lumped honesty guard: massType==1 is parsed + serialized (API stable) but
  // the HRZ path lands P3 — getMass errors clearly ONCE (per process) and falls
  // back to the consistent mass so the run is never silently wrong-but-quiet.
  // Not serialized.  // Ladruno
  static bool warnedLumped;

  // URI ordinal defense (C++/wire surface): the parser rejects 'uri' until P2,
  // but the public constructor and recvSelf coerce any non-STD ordinal to STD
  // with a process-once notice, so Print/recorder metadata can never claim uri
  // while computing std. Not serialized.  // Ladruno
  static bool warnedUriCoerce;

  // -------- static workspace (sized 60) --------
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;

  // -------- private methods --------
  void buildGeometryCache(void);                 // fill N/dNdx/dv + detJ>0 gate (F1/F11)
  bool cacheUsable(const char *where);           // guard: emit-once + false when bad
  void coerceFormulationToStd(const char *where);// URI ordinal defense (F8)
  void formResidual(void);                       // B^T sigma (+ damping + body) Gauss loop
  void formStiffness(int initialFlag);           // shared B^T D B assembly (tangent/initial)
  void computeConsistentMass(void);              // 27-pt consistent mass into `mass`
  void formInertiaResidual(void);                // rho N_a N_b accel momentum pass (residual only)
  double computeVolume(void);                    // sum of cached dv
  void gatherCoords(double X[NEN][3]);           // reference nodal coords -> X[20][3]

  // string <-> enum helpers
  static const char *formulationName(Formulation f);
};

#endif
