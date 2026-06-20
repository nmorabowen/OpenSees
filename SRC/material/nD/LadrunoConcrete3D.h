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

// Ladruno: CDPM2-grade solid-concrete plastic-damage nDMaterial (3D). A thin OpenSees
// wrapper over the header-only, OpenSees-free LadrunoConcrete3DKernel.h (the SAME verified
// core that the g++ oracle-numeric-dump self-check pins to the numpy oracle). The kernel
// carries: an effective-stress Menetrey-Willam 3-invariant surface (friction m0, eccentricity
// e from Kupfer fcc/fc, Willam-Warnke Lode r), confinement-aware hardening (qh1/qh2/kappa_p +
// ductility xh), non-associated flow (dilatancy Df), a semi-implicit return map, and the dual
// scalar damage omega_t/omega_c on the spectral tension/compression split with crack-band Gf/Gc
// regularization + automatic unilateral crack-closure. See
// Ladruno_implementation/31_ladruno_concrete3d_adr.md. classTag 33017.
//
// SCOPE: all dimensional views ("ThreeDimensional" + the Phase-2 reduced PlaneStrain /
// AxiSymmetric / PlateFiber / PlaneStress) are served by this ONE class through a `dim` mode
// (the LadrunoJ2 "one core, many views" doctrine): the kernel return map always runs on the full
// 6-component tensor, and the element-facing strain/stress/tangent are mapped to the reduced
// ordering. The element requests a view via getCopy(type); the OPS parser always builds the 3D
// prototype. PlaneStress / PlateFiber enforce the out-of-plane sigma_22 = 0 by a nested Newton on
// eps_22 + static condensation of the 33 dof (de Souza Neto et al. sec 9.2.3). The finite-strain
// view is also free via nDMaterial LogStrain wrapping the 3D material. Tier-1 implicit only
// (IMPL-EX / Duvaut-Lions viscosity are P3). The CDPM2 consistent tangent is NON-SYMMETRIC
// (non-associated flow + the semi-implicit theta-freeze even associated) => use an UNSYMMETRIC
// solver (e.g. `system FullGeneral` / `UmfPack`); the parser warns. NB unconfined plane-STRESS
// softening is snap-backy and the sigma_22 = 0 nested Newton can stall on the post-peak limit
// point (Tier-3 explicit territory) — the reduced views are robust pre-peak / under confinement.
//
// Tensor convention (LadrunoJ2 lineage): the kernel stores symmetric tensors as 6 components
// {00,11,22,01,12,02} with TRUE tensor off-diagonals. The element passes ENGINEERING shear
// strain (gamma = 2 eps), so the wrapper halves the shear strains on the way in and doubles
// them on the way out; stresses are true components (no factor). The kernel tangent is in the
// tensor convention (dsig_ij = 2G deps_ij), so its shear COLUMNS are halved to give the
// element's d(sigma)/d(engineering strain).
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoConcrete3D_h
#define LadrunoConcrete3D_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LadrunoConcrete3D : public NDMaterial {
 public:
  // dimensional views (element-facing ordering, engineering shear)
  enum { DIM_3D = 0,         // {00,11,22,01,12,02}      order 6
         DIM_PSTRAIN,        // {00,11,01}  (eps22=0)    order 3
         DIM_AXISYM,         // {00,11,22,01}            order 4
         DIM_PLATEFIBER,     // {00,11,01,12,02} sig22=0 order 5
         DIM_PSTRESS };      // {00,11,01}  sig22=0      order 3

  LadrunoConcrete3D();
  LadrunoConcrete3D(int tag, double E, double nu, double fc, double ft, double Gf, double Gc,
                    double e, double Df, double As,
                    double qh0, double Hp, double Ah, double Bh, double Ch, double Dh,
                    double rho, double lch, bool autoReg, bool implex = false, int dimMode = DIM_3D);
  ~LadrunoConcrete3D();

  const char* getClassType(void) const { return "LadrunoConcrete3D"; }

  int setTrialStrain(const Vector& strain);
  int setTrialStrain(const Vector& v, const Vector& r);
  int setTrialStrainIncr(const Vector& v);
  int setTrialStrainIncr(const Vector& v, const Vector& r);

  const Matrix& getTangent(void);
  const Matrix& getInitialTangent(void);
  const Vector& getStress(void);
  const Vector& getStrain(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial* getCopy(void);
  NDMaterial* getCopy(const char* type);
  const char* getType(void) const;          // dim-dependent (see .cpp)
  int getOrder(void) const { return ncomp; }

  double getRho(void) { return rho; }

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  // ---- material parameters (POSITIVE strengths; compression negative) ----
  double E, nu, fc, ft;        // elasticity + uniaxial strengths
  double Gf, Gc;               // tensile / compressive fracture energies (crack-band)
  double ecc, m0;              // deviatoric eccentricity + derived friction
  double Df;                   // dilatancy (non-associated flow)
  double As;                   // confinement-ductility amplitude (Eq.56)
  double qh0, Hp;              // hardening (Eq.30-31)
  double Ah, Bh, Ch, Dh;       // ductility-measure params (Eq.33-36)
  double rho;                  // mass density
  double lchFixed;             // characteristic length when -autoRegularization is OFF
  bool   autoReg;              // pull lch from the active element each step
  bool   implex;               // Tier-2 IMPL-EX (P3 robustness; default OFF = Tier-1 implicit)

  // ---- dimensional view (element-facing ordering; the kernel is always 3D) ----
  int    dim;                  // DIM_*
  int    ncomp;                // element vector order (3..6)
  int    vmap[6];              // reduced index a -> full 6-comp tensor index
  bool   condense;            // enforce sigma_22 = 0 (PlaneStress / PlateFiber)
  double cEps22;               // committed out-of-plane tensor strain (condensed modes)

  // ---- committed kernel state (tensor components; shear = TRUE tensor) ----
  double eps_n[6];             // committed total strain
  double sig_n[6];             // committed NOMINAL stress
  double sigEff_n[6];          // committed EFFECTIVE stress (drives the next return)
  double kp_n;                 // committed kappa_p
  double etmax_n;              // committed max equiv-strain history (Eq.43)
  double kdt1_n, kdt2_n;       // committed tensile damage histories (Eq.44-45)
  double kdc_n, kdc1_n, kdc2_n;// committed compressive damage histories (Eq.47-49)
  // committed IMPL-EX bookkeeping (the implicit damage + per-variable increments + dt, for the next
  // step's extrapolation; unused when !implex)
  double wt_n, wc_n;           // committed IMPLICIT dual damage
  double dwt_n, dwc_n;         // committed implicit damage increments
  double depl_n[6];            // committed implicit plastic-strain increment (tensor)
  double dtn_n;                // committed time step

  // ---- trial state (recomputed each setTrialStrain) ----
  double strain6[6];           // trial total strain (tensor)
  double stress6[6];           // trial NOMINAL stress (true tensor components)
  double sigEff6[6];           // trial effective stress
  double kp_t;
  double etmax_t, kdt1_t, kdt2_t, kdc_t, kdc1_t, kdc2_t;
  double wt_t, wc_t, dwt_t, dwc_t, depl_t[6], dtn_t;   // trial IMPL-EX bookkeeping
  double Dtan6[6][6];          // trial damaged tangent (kernel TENSOR convention)
  double omegaT, omegaC;       // trial damage variables (for recorders)
  int    lastStatus;           // 0 converged, 2 no-converge (honest off-surface flag)

  // helpers
  void integrate(bool doTangent);
  void setupDim(void);             // fill ncomp/vmap/condense + size the return buffers
  void condenseTangent(void);      // static condensation of the 33 dof (sigma_22 = 0)

  // element-facing return buffers
  Vector stressOut;
  Vector strainOut;
  Matrix tangentOut;
};

#endif // LadrunoConcrete3D_h
