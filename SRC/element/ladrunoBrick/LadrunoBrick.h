/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LadrunoBrick — unified 8-node hexahedral solid element (Ladruno fork).
//
// One element class, one classTag (ELE_TAG_LadrunoBrick = 33002), with the
// anti-locking formulation chosen at construction via a single selector:
//
//   element('LadrunoBrick', tag, n1..n8, matTag, '-formulation', <std|bbar|uri|eas>)
//
// v1 (small strain, geometrically linear):
//   std  — full 2x2x2 Gauss displacement  (reproduces upstream Brick)
//   bbar — mean-dilatation B-bar          (reproduces upstream bbarBrick)
//   uri  — 1-pt reduced + hourglass        (cheap explicit hex; new)
//   eas  — bbar + statically-condensed enhanced assumed strain (SSPbrick port)
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

class LadrunoBrick : public Element {

 public:

  // Formulation selector (the single -formulation axis).
  enum class Formulation { STD, BBAR, URI, EAS };

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
               int geomMethodID = 0);   // 0=linear, 2=finite (SolidTransformation method id)

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

  // -geom finite (v3, updated-Lagrangian). isFinite() is true when theGeom
  // reports a DeformationGradient strain measure: the element computes the full
  // F per GP, drives the material via setTrialF(F), and assembles the spatial
  // internal force + material tangent + geometric (initial-stress) stiffness.  // Ladruno
  bool isFinite(void) const;

  // True for the single-integration-point formulations (eas, uri with stiffness
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

  // eas (v2) — Stabilized Single-Point: bbar (mean-dilatation Bnot) + statically
  // condensed enhanced assumed strain. Ported from UWelements/SSPbrick. The
  // enhanced modes are condensed analytically at element setup using the INITIAL
  // material tangent, so easBnot/easKstab are CONSTANT (no per-step alpha state)
  // and cure both shear AND volumetric locking across all nu. buildEAS computes
  // them once (setDomain); formEAS assembles K = Kstab + V*Bnot^T C Bnot and
  // f = Kstab*u + V*Bnot^T*sigma - bodyForce.  // Ladruno
  void buildEAS(void);
  void formEAS(int tang_flag, bool useInitialTangent);
  Matrix *easBnot;     // 6x24 mean-dilatation B (constant; strain = Bnot*u)
  Matrix *easKstab;    // 24x24 condensed enhanced-strain stabilization (constant)
  double  easVol;      // element volume (8*Jo + higher-order terms)
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
  // FB perturbation energy ½κ·Σ q_aι², for eas it is the condensed-stabilization
  // energy ½·u_core·Kstab·u_core. Zero for std/bbar (fully integrated, no
  // hourglass) and for physical (assumed strain folded into the strain energy).
  // Stateless; reported via the "hourglassEnergy" response. NOTE: uri 'viscous'
  // DISSIPATES rather than stores energy — for that flavour hourglassEnergy()
  // returns the committed accumulator hgDissipated instead of an instantaneous
  // stored value.  // Ladruno
  double hourglassEnergy(void);

  // Work done against the FB viscous-hourglass damping force over the step now
  // being committed (uri + Hourglass::VISCOUS only): ΔE = c_visc·Σ q̇_aι·Δq_aι,
  // clamped ≥ 0. Mirrors the force assembled in formResidAndTangent's viscous
  // branch; accumulated into hgDissipated by commitState.  // Ladruno
  double viscousHourglassIncrement(void);

  // string <-> enum helpers (shared by OPS factory + Print + sendSelf)
  static const char *formulationName(Formulation f);
};

#endif
