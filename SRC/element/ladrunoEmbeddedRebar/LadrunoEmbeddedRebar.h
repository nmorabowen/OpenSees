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

// Ladruno: embedded-reinforcement COUPLING element (ADR 20, v2 slice 2).
//
//   Ties ONE rebar node to a host solid element's nodes through precomputed
//   shape-function weights N_i, so a rebar mesh can be embedded in a NON-matching
//   concrete mesh. The element is a pure COUPLING (bond): the rebar's own axial
//   stiffness lives on a separate rebar element (corotTruss/beam) along the bar;
//   this element only enforces the bar<->concrete attachment.
//
//   External nodes: [ rebarNode, hostNode_1, ..., hostNode_M ].
//   Gap (relative displacement of the rebar point):
//       g = u_rebar - sum_i N_i u_host_i        (ndm-vector)
//   Decompose along the bar axis `dir` (unit): s = g . dir  (axial slip),
//   g_t = g - s*dir (transverse). Coupling traction:
//       t = F_axial(s) * dir + kt * g_t
//   F_axial = perfect: kAxial*s ; bond: bondMat->stress(s) * bondScale
//   (bondScale = perimeter * L_trib converts the tau-s law to a nodal force).
//   Resisting force  P = B^T t,  tangent  K = B^T D B,  with
//       B = [ I | -N_1 I | ... | -N_M I ],
//       D = (dF_axial/ds) dir(x)dir + kt (I - dir(x)dir)   (anisotropic, ADR D6).
//
//   MODE (ADR D2): this is MODE P (penalty coupling) — works in BOTH implicit and
//   explicit (it adds stiffness only, mass stays where it is). Perfect bond is
//   Mode P with a large axial penalty; bond-slip plugs LadrunoBondSlip into the
//   axial slot. MODE T (transformation / DOF-elimination, no penalty) is DEFERRED:
//   it needs a MULTI-retained constraint (u_rebar = sum N_i u_host_i), which the
//   stock single-retained MP_Constraint cannot express — a separate larger lift.
//
//   Element type of the REBAR itself (truss vs beam, ADR D6) is the user's choice
//   on the rebar mesh; this coupling element is identical for both (it couples
//   translations only). The host shape functions N_i are precomputed by the
//   generator (apeGmsh) from the host type + reference-config natural coords, so
//   the element is HOST-AGNOSTIC (it never needs to know the host element type).
//
//   classTag ELE_TAG_LadrunoEmbeddedRebar = 33005 (Ladruno private element band).
//   See Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoEmbeddedRebar_h
#define LadrunoEmbeddedRebar_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class UniaxialMaterial;
class Channel;
class FEM_ObjectBroker;
class Response;

class LadrunoEmbeddedRebar : public Element
{
 public:
  LadrunoEmbeddedRebar(int tag, int ndm, int rebarNode, const ID& hostNodes,
                       const Vector& shape, const Vector& dir,
                       double kt, double bondScale,
                       UniaxialMaterial* bondMat, double kAxialPerfect,
                       int hostEleTag = -1, bool ktAuto = false,
                       double ktAlpha = 0.0, bool corot = false,
                       const Vector* shapeB = 0);
  LadrunoEmbeddedRebar();
  ~LadrunoEmbeddedRebar();

  const char* getClassType(void) const { return "LadrunoEmbeddedRebar"; }

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

  // parallel
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

 private:
  int ndm;                  // 2 or 3
  int nHost;                // number of host nodes M
  int nDOF;                 // (1 + M) * ndm
  ID connectedNodes;        // [rebar, host_1..host_M]
  Vector Nshape;            // host shape-function weights N_i (size M)
  Vector dir;               // REFERENCE unit bar-axis direction (size ndm, frozen)

  // ADR 20 §10.5 — co-rotated bar axis. With `-corot`, the axial/transverse split
  // is taken against the CURRENT bar direction instead of the frozen reference
  // `dir`, so the coupling traction stays frame-objective under large host
  // rotation (the only true large-rotation defect — the gap g and the weights
  // N_i(ξ) are already objective). Direction = secant between the embed point
  // (weights Nshape) and a second point B along the bar (weights NshapeB), using
  // CURRENT host node positions: dirCur = normalize(Σ NshapeB_i x_i − Σ Nshape_i x_i).
  // v1 drops the ∂dir/∂u consistent-tangent term (EICR practice: exact for
  // explicit, converges for moderate per-step rotation under step-halving).
  bool corot;               // -corot: recompute the bar axis each step
  Vector NshapeB;           // weights at point B along the bar (size M; corot only)
  Vector dirCur;            // last co-rotated unit bar axis (reported by "dir")
  void currentBarAxis(Vector& d);  // fill d with the working bar axis (dir or secant)
  double kt;                // transverse penalty stiffness (resolved value)
  double bondScale;         // perimeter * L_trib (tau -> axial force)
  double kAxialPerfect;     // perfect-bond axial penalty (used when bondMat==0)
  UniaxialMaterial* bondMat;// optional axial bond-slip law (owns a getCopy)

  // ADR 20 §10.2a — auto-scaled transverse penalty. With `-kt auto`, kt is
  // resolved (lazily, on first stiffness/force assembly) from the host's own
  // initial stiffness so conditioning is mesh/material-independent:
  //     kt = ktAlpha * max_i |K_host(i,i)|   ( ~ ktAlpha * E_host * lch ).
  // Lazy (not at setDomain) to dodge element setDomain-ordering: by first
  // assembly the host has setDomain'd and getInitialStiff() is valid.
  int hostEleTag;           // -host element tag (>=0) for auto-kt; -1 if none
  bool ktAuto;              // -kt auto: resolve kt from the host stiffness
  double ktAlpha;           // dimensionless multiplier for the auto kt
  bool ktResolved;          // transient: auto kt already resolved this run
  void resolveAutoKt(void); // resolve kt from the host element if needed

  Node** theNodes;          // size 1 + M
  Matrix* K;                // nDOF x nDOF
  Vector* P;                // nDOF
  Matrix* M0;               // zero mass (nDOF x nDOF)
  Response* bondEnergyResp; // cached bondMat "energy" sub-response (ADR §10.2b)

  // committed/trial scalar slip is held by bondMat; nothing else is path-dep.
  void formBandTraction(Vector& g, double& s, Vector& gt,
                        double& Faxial, double& kAxial); // helper
};

#endif
