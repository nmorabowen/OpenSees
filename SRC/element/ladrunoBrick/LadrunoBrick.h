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

// LadrunoBrick — unified 8-node hexahedral solid element (Ladruno fork).
//
// One element class, one classTag (ELE_TAG_LadrunoBrick = 33002), with the
// anti-locking formulation chosen at construction via a single selector:
//
//   element('LadrunoBrick', tag, n1..n8, matTag, '-formulation', <std|bbar|uri|ssp|eas>)
//
// v1 (small strain, geometrically linear):
//   std  — full 2x2x2 Gauss displacement  (reproduces upstream Brick)
//   bbar — mean-dilatation B-bar          (reproduces upstream bbarBrick)
//   uri  — 1-pt reduced + hourglass        (cheap explicit hex; new)
//   ssp  — stabilized single-point: bbar + statically-condensed stabilization
//          (SSPbrick port). NB this is NOT enhanced assumed strain.
//   eas  — true Simo-Rifai enhanced assumed strain (9 internal params, inner
//          Newton + static condensation; full 2x2x2, 8 live GPs; ADR 19). Robust
//          on smooth/homogeneous response; for notched/localization inelasticity
//          use ssp/bbar (EAS stabilization is ADR 20).
//
// The kernel is carved into three seams (kinematics ledger / geometry method /
// material adaptor) so corotational (v2) and finite-strain (v3) drop in without
// a rewrite — see Ladruno_implementation/09_ladruno_brick.md and
// solid_transformation_wrapper.md.
//
// Baseline (std/bbar kernels) adapted from Ed Love's Brick / BbarBrick. The
// upstream Brick::setParameter damping loop bug (i<4 over an 8-entry array) is
// fixed here.  // Ladruno

#ifndef LadrunoBrick_h
#define LadrunoBrick_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Damping.h>

class SolidTransformation;   // seam 2/3: geometry-method layer (linear/corot/finite)
class Response;              // Ladruno — cached material "damage" query (Tier-A Kstab)

class LadrunoBrick : public Element {

 public:

  // Formulation selector (the single -formulation axis).
  // Ordinals are serialized (packed into idData(28) by sendSelf), so the order is
  // load-bearing for DB/parallel back-compat: SSP keeps ordinal 3 (the slot the
  // old single-point "EAS" used), so legacy streams reload as the same element.
  // EAS is true Simo-Rifai enhanced assumed strain (ADR 19) and takes the NEW
  // ordinal 4 — distinct from the renamed single-point SSP=3.  // Ladruno
  enum class Formulation { STD, BBAR, URI, SSP, EAS };

  // Hourglass stabilization flavour (uri only).
  enum class Hourglass { VISCOUS, STIFFNESS, PHYSICAL };

  // null constructor (broker)
  LadrunoBrick();

  // full constructor
  LadrunoBrick(int tag,
               int node1, int node2, int node3, int node4,
               int node5, int node6, int node7, int node8,
               NDMaterial &theMaterial,
               Formulation formulation = Formulation::STD,
               double b1 = 0.0, double b2 = 0.0, double b3 = 0.0,
               int massType = 0,
               Hourglass hgType = Hourglass::PHYSICAL,
               double hgCoeff = 0.0,
               Damping *theDamping = 0,
               int geomMethodID = 0,   // 0=linear, 2=finite (SolidTransformation method id)
               double stabBeta = 0.0); // Ladruno — ADR 20 eas tangent-regularization beta

  virtual ~LadrunoBrick();

  const char *getClassType(void) const { return "LadrunoBrick"; };

  // domain
  void setDomain(Domain *theDomain);
  int  setDamping(Domain *theDomain, Damping *theDamping);

  // topology
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);

  // Reference-config element size for crack-band / regularized-softening
  // materials (ASDConcrete3D, Lemaitre LCH_REF). The Element base default
  // returns the MIN inter-node distance, which under-sizes the band on a
  // distorted hex and over-softens. We return the edge of an equal-volume
  // cube, lch = cbrt(V) — geometry-true (the hex analogue of BezierTet10's
  // cbrt(6V) and BezierTri6's sqrt(2A)). Degenerate V<=0 falls back to base.  // Ladruno
  double getCharacteristicLength(void);

  // Ladruno (ADR 20 §9): trilinear 8-node hex shape weights at natural coord
  // xi = (ξ,η,ζ) ∈ [-1,1]³, for embedded-reinforcement coupling. N sized to 8.  // Ladruno
  int getInterpolationWeights(const Vector &xi, Vector &N);

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

  // -------- private attributes --------
  ID connectedExternalNodes;          // 8 node numbers
  Node *nodePointers[8];              // pointers to 8 nodes
  NDMaterial *materialPointers[8];    // pointers to 8 materials (one per Gauss point)

  Formulation formulation;            // the -formulation selector
  Hourglass hourglassType;            // uri only
  double hourglassCoeff;              // uri only

  // Cumulative viscous-hourglass dissipation (uri + Hourglass::VISCOUS only).
  // The FB viscous hourglass force damps the spurious modes and stores NO energy,
  // so hourglassEnergy() cannot report it instantaneously. Instead we integrate
  // the work done against that damping force over committed steps:
  //   ΔE = c_visc·Σ_{a,i} q̇_aι·Δq_aι   (q̇ from getTrialVel, Δq from Δu).
  // hgDissipated accumulates it (committed, serialized); uPrevCommit holds the
  // last committed nodal displacement (the work baseline); hgPrevValid guards the
  // first commit / a post-recv reseed so no spurious increment is booked.  // Ladruno
  double hgDissipated;                // cumulative dissipated viscous-hourglass energy
  double uPrevCommit[24];             // last committed nodal displacement (work baseline)
  bool   hgPrevValid;                 // false until uPrevCommit is first seeded

  double b[3];                        // body forces
  double appliedB[3];                 // body forces applied with load
  int applyLoad;

  Vector *load;
  Matrix *Ki;
  int massType;                       // 0 consistent, 1 lumped

  Damping *theDamping[8];

  // Tier-A damage-scaled hourglass stabilization (softening support). For the
  // single-point STIFFNESS-stabilized formulations (ssp, uri+stiffness) the
  // constant elastic Kstab over-stiffens a cracked element and blocks crack
  // localization under a softening material (ASDConcrete3D). damageResponse is a
  // cached query of materialPointers[0]'s "damage" channel (built in setDomain,
  // for ASDConcrete3D = getAvgDamage() = [d_tension, d_compression]); damageScale()
  // reads it and returns max(floor, 1 - max(d_i)) to degrade Kstab. 0 (null) when
  // the material has no "damage" channel => damageScale() falls back to 1.0 (the
  // original constant elastic Kstab). Not serialized — rebuilt in setDomain on the
  // receive side.  // Ladruno
  Response *damageResponse;

  // Geometry-method layer (seam 2/3). v1 = SolidTransformationLinear (identity):
  // localizeDisp / globalizeForce / globalizeStiff are pass-throughs, so routing
  // through it changes nothing. corot (v2) / finite (v3) override on this same
  // interface. Owned by the element; rebuilt in recvSelf from a method id.  // Ladruno
  SolidTransformation *theGeom;

  // -------- static workspace --------
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;

  // quadrature data (full 2x2x2)
  static const double root3;
  static const double one_over_root3;
  static const double sg[2];
  static const double wg[8];

  // local nodal coordinates
  static double xl[3][8];

  // -------- private methods --------
  void formInertiaTerms(int tangFlag);
  void formResidAndTangent(int tang_flag);
  void computeBasis(void);

  // Reference-config element volume by 2x2x2 Gauss integration of |J|
  // (V = Σ wg·detJ over the 8 GPs). Formulation-independent — does NOT reuse
  // sspVol, which only exists after buildSSP and only for ssp. Backs
  // getCharacteristicLength().  // Ladruno
  double computeVolume(void);

  // -geom finite (v3, updated-Lagrangian). isFinite() is true when theGeom
  // reports a DeformationGradient strain measure: the element computes the full
  // F per GP, drives the material via setTrialF(F), and assembles the spatial
  // internal force + material tangent + geometric (initial-stress) stiffness.  // Ladruno
  bool isFinite(void) const;

  // True for the single-integration-point formulations (ssp, uri with stiffness
  // or viscous hourglass): the constitutive response is evaluated ONCE at the
  // centroid (material slot 0) and the other 7 slots are output mirrors — so we
  // skip 7 redundant (and, for materials like ASDConcrete3D, expensive) return
  // maps in update() and mirror slot 0 in the per-GP output. std/bbar/uri-physical
  // and -geom finite genuinely use all 8 Gauss points, so this is false.  // Ladruno
  bool isSinglePoint(void) const;

  int  updateFinite(void);                       // per GP: F = I + Σ uⱼ⊗∇ₓNⱼ → setTrialF
  void formResidAndTangentFinite(int tang_flag); // ∫Bᵀσ dv + ∫BᵀcB dv + ∫GᵀΣG dv
  // deformation gradient F (row-major [9]) and det at GP `gp` from reference
  // shape gradients shpRef[0..2][a] = ∂Nₐ/∂Xᵢ and the nodal trial displacements.
  double deformationGradient(const double shpRef[4][8], double F[9]);

  // F-bar (bbar + finite, dSNPO eq 15.5/15.10). Assumes computeBasis() has set
  // xl. Returns J0 = det F0 with F0 the deformation gradient at the element
  // centroid (natural coords 0,0,0). If G0 != 0, also fills the centroid spatial
  // gradient operator G0[k][b] = ∂N_b/∂x_k|_centroid (from F0⁻¹) for the eq 15.10
  // F-bar coupling term.  // Ladruno
  double centroidFbar(double (*G0)[8] = 0);

  // Seam 0+2: refresh theGeom from current geometry (reference + current nodal
  // coords) and return the localized (core-frame) 24-dof trial displacement.
  // For -geom linear this equals the global trial displacement.  // Ladruno
  const Vector &computeLocalDisp(void);

  // Seam 1 (kinematics ledger) + core formulation B selection.
  // computeB: the linear small-strain B at a node (std / uri / shear rows of bbar).
  const Matrix &computeB(int node, const double shp[4][8]);
  // computeBbar: the mean-dilatation B-bar (bbar formulation).
  const Matrix &computeBbar(int node, const double shp[4][8], const double shpBar[4][8]);

  // bbar helper: element-averaged shape-function gradients (mean dilatation).
  void computeShapeBar(double shpBar[4][8], double Shape[4][8][8],
                       const double dvol[8], double volume);

  // uri (1-point reduced integration + hourglass control) helpers.
  // formUri assembles stiff (+ resid when !useInitialTangent) at the centroid
  // plus Flanagan-Belytschko stiffness hourglass control.
  void formUri(int tang_flag, bool useInitialTangent);
  // FB-corrected hourglass vectors (orthogonal to the linear field => the
  // hourglass control is consistent: zero for any constant-strain/rigid state).
  // beta is the normalization (1.0 for stiffness scaling; 1/8 for the physical
  // assumed strain, where gamma must equal the true hourglass-mode coefficient).
  void computeGammaHourglass(const double b[3][8], double gamma[4][8], double beta);
  // the 4 raw hourglass base vectors (nodal values of {yz, zx, xy, xyz} per
  // Belytschko h_alpha = [eta*zeta, zeta*xi, xi*eta, xi*eta*zeta], 8.7.25).
  static const double hg0[4][8];

  // physical (uri + Belytschko-Bindeman isochoric assumed-strain hourglass,
  // sec 8.7.6/8.7.8) integrated at the full 2x2x2 rule (8.7.20). Assembles
  // stiff (+ resid when !useInitialTangent) with the assumed-strain B-bar.
  void formPhysical(int tang_flag, bool useInitialTangent);

  // ssp — Stabilized Single-Point: bbar (mean-dilatation Bnot) + a statically
  // condensed stabilization. Ported from UWelements/SSPbrick. NB despite the
  // SSPbrick derivation this is NOT enhanced assumed strain — there is no
  // per-step alpha state; the stabilization modes are condensed analytically at
  // element setup using the INITIAL material tangent, so sspBnot/sspKstab are
  // CONSTANT and cure both shear AND volumetric locking across all nu. buildSSP
  // computes them once (setDomain); formSSP assembles K = Kstab + V*Bnot^T C Bnot
  // and f = Kstab*u + V*Bnot^T*sigma - bodyForce. (True Simo-Rifai EAS, with live
  // alpha state and per-iteration condensation, is ADR 19's '-formulation eas'.)  // Ladruno
  void buildSSP(void);
  void formSSP(int tang_flag, bool useInitialTangent);
  Matrix *sspBnot;     // 6x24 mean-dilatation B (constant; strain = Bnot*u)
  Matrix *sspKstab;    // 24x24 condensed stabilization (constant)
  double  sspVol;      // element volume (8*Jo + higher-order terms)

  // eas — TRUE Simo-Rifai enhanced assumed strain (ADR 19). Unlike ssp this is a
  // mixed element: the compatible strain B*u is enriched by an enhanced field
  // M(xi)*alpha with 9 element-internal parameters alpha, solved each form pass by
  // an inner Newton enforcing int M^T sigma = 0, then statically condensed
  // (K* = Kdd - Kda Kaa^-1 Kad). Full 2x2x2 integration (8 LIVE material points).
  // Ported from the 2-D EnhancedQuad (Wilson incompatible-mode lineage):
  //   M_i(xi) = sym[ (j0/j(xi)) J0^-T E_i(xi) ]   (ADR 19 eq E.8)
  // = 3 natural-direction bubbles x 3 dofs. j0/J0inv are the centroid Jacobian
  // det/inverse (cached in setDomain via buildEAStrue). alpha is committed state
  // (alphaCommit); commit/revert/sendSelf carry it. v1 = small strain only
  // (-geom corot/finite + eas are parser-reserved).  // Ladruno
  void buildEAStrue(void);                                 // cache centroid J0inv/j0
  void formEAStrue(int tang_flag, bool useInitialTangent); // inner-Newton + condensation
  void computeMenh(const double gp[3], double jdet, Matrix &M);  // 6x9 enhanced operator
  Vector alpha;        // 9 enhanced parameters (trial; solved each form pass)
  Vector alphaCommit;  // committed enhanced parameters (serialized)
  Matrix easJ0inv;     // 3x3 centroid Jacobian inverse (mode map; cached)
  double easJ0det;     // centroid Jacobian determinant j0 (cached)
  // ADR 20 — optional Kaa tangent regularization for inelastic localization.
  // Kaa_stab = Kaa + beta*easKaa0, applied ONLY at the .Solve() operator (never to
  // the accumulated Kaa) so the residual int M^T sigma is untouched => the converged
  // state is beta-independent in the convex regime (a modified-Newton, not a physics
  // change). easKaa0 = int M^T C0 M dV is the elastic enhanced stiffness (geometry-
  // only in small strain), cached once in buildEAStrue alongside easJ0inv.  // Ladruno
  double easStabBeta;  // -stab beta (default 0 = bit-identical to bare eas)
  Matrix easKaa0;      // 9x9 elastic enhanced stiffness int M^T C0 M dV (cached)
  // Fill the assumed-strain B (Bbar[node][6][3], Voigt {xx,yy,zz,xy,yz,zx}) at a
  // Gauss point; returns |J|. gamma/bC are the (precomputed) hourglass vectors
  // and centroid gradients. Implements eq 8.7.26.
  double assumedStrainB(const double gp[3], const double gamma[4][8],
                        const double bC[3][8], double Bbar[8][6][3]);
  // node natural coordinates (Brick order), for analytic dN/dxi and dh/dxi.
  static const double natCoord[8][3];

  Matrix transpose(int dim1, int dim2, const Matrix &M);

  // Recoverable elastic hourglass / stabilization energy at the current trial
  // state (the GLSTAT spurious-mode diagnostic): for uri 'stiffness' it is the
  // FB perturbation energy ½κ·Σ q_aι², for ssp it is the condensed-stabilization
  // energy ½·u_core·Kstab·u_core. Zero for std/bbar (fully integrated, no
  // hourglass) and for physical (assumed strain folded into the strain energy).
  // Stateless; reported via the "hourglassEnergy" response. NOTE: uri 'viscous'
  // DISSIPATES rather than stores energy — for that flavour hourglassEnergy()
  // returns the committed accumulator hgDissipated instead of an instantaneous
  // stored value.  // Ladruno
  double hourglassEnergy(void);

  // Tier-A multiplier degrading the constant elastic hourglass stabilization
  // (formSSP / formUri stiffness) with the material's current damage:
  //   s = max(floor, 1 - max(d_i)),  d_i = the material's "damage" response
  // (ASDConcrete3D: [d_tension, d_compression]). Returns 1.0 when the material
  // reports no "damage" channel (damageResponse == 0) => the original elastic
  // Kstab. ASDConcrete3D returns the secant by default (+IMPLEX) so d_i is
  // positive & monotone and s stays in [floor, 1] — never negative.  // Ladruno
  double damageScale(void);

  // Work done against the FB viscous-hourglass damping force over the step now
  // being committed (uri + Hourglass::VISCOUS only): ΔE = c_visc·Σ q̇_aι·Δq_aι,
  // clamped ≥ 0. Mirrors the force assembled in formResidAndTangent's viscous
  // branch; accumulated into hgDissipated by commitState.  // Ladruno
  double viscousHourglassIncrement(void);

  // string <-> enum helpers (shared by OPS factory + Print + sendSelf)
  static const char *formulationName(Formulation f);
};

#endif
