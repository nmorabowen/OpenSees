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

// LadrunoNorSand: NorSand in the Andrade & Borja 2006 (IJNME 65, 3-invariant) form, as an nDMaterial.
//
//   - ALL the constitutive mathematics lives in the header-only, OpenSees-free kernel
//     LadrunoNorSandKernel.h (namespace ladruno_norsand): the 4x4 spectral return in principal
//     elastic strains, the nested pi_i solve, Delta-lambda backtracking, substepping, and the
//     closed-form NON-symmetric consistent tangent. This class is the thin shell: it holds the
//     Params and the committed/trial kernel State, maps between the OpenSees Voigt vectors
//     (engineering shear) and the kernel's tensor components, and signals a refused step the
//     fork's way (LADRUNO_MATERIAL_REFUSED at the trial, the WP-99 commit-refusal seam at commit).
//   - LadrunoNorSand3D (33024) and LadrunoNorSandPlaneStrain (33025) are thin dimension wrappers
//     (LadrunoSANISAND3D / PlaneStrain pattern). The base class (33023) is itself a complete
//     3D material; its getCopy(type) hands out the wrappers.
//   - The tangent is NON-symmetric (non-associated flow + pi_i*(p, Omega) + the v term): pair it
//     with an unsymmetric solver. This class never selects one.
//
// Sign convention: kernel and OpenSees agree (compression NEGATIVE, tension positive), so there is
// NO sign flip anywhere here (unlike LadrunoSANISAND, which is compression-positive internally).
//
// Voigt <-> tensor mapping (the only place shear factors appear):
//   OpenSees 3D order {11,22,33,12,23,13}, engineering shear strain gamma = 2*eps_tensor.
//   Kernel order      {00,11,22,01,12,02}, TENSOR components (LadrunoJ2Kernel.h convention).
//   The index orders coincide (12<->01, 23<->12, 13<->02), so the map is the identity on indices:
//     strain  eps_tensor[i] = e[i]            (i < 3)        eps_tensor[i] = e[i]/2   (i >= 3)
//     stress  sigma_voigt[i] = sigma_tensor[i]               (shear stress is the same number)
//     tangent T[a][b] = C[a][b] * w_b,  w_b = 1 (b < 3), 1/2 (b >= 3)       (rows unchanged)
//   where C = d sigma_tensor / d eps_tensor is the kernel's 6x6 derivative taken with respect to
//   the six independent tensor components (eps_01 an independent variable). PlaneStrain picks the
//   index set {0,1,3} (eps_22 = eps_12 = eps_02 = 0); sigma_22 is kept and exposed by getStressZZ.
//
// Refusal (WP-99 pattern, mirrors LadrunoSANISAND):
//   - setTrialStrain returns LADRUNO_MATERIAL_REFUSED when the kernel refuses (committed state
//     untouched; the trial stress/tangent are those of the last committed state). A forwarding
//     element cuts the step. NO latch is set at the trial (a refused trial may be retried smaller).
//   - If a refused trial then reaches commitState() (a host element that DISCARDS the code), the
//     commit is refused out of band (ladrunoNoteCommitRefusal) and this point LATCHES: every later
//     trial and commit is refused until revertToStart(). getCopy propagates the latch.
//   - Bounded work: the kernel's caps (local iterations, line search, pi_i scan, 2^8 substeps).
//   - The refusal reason is recorded at two levels: the kernel's `refusal` (SUBSTEPS_EXHAUSTED for every
//     refusal that went down the substep ladder) and `finest` / `finest_sub`, the cause at the finest
//     level (StepInfo, from detail::step_ex). Both reach the one-time WARNING and the `refusal` /
//     `stepInfo` responses.
//
// Tangent: a substepped increment returns the CHAINED consistent tangent, d sigma_final / d (TOTAL
// strain increment), propagated through every sub-increment (owner decision 2026-10-01; plan 2.8).
// It is a derivative with respect to the six INDEPENDENT tensor components (eps_01 etc. independent),
// so an engineering-shear caller halves the shear COLUMNS (done in getTangent; unlike LadrunoJ2Kernel).
//
// Initial state: -pi0 is REQUIRED (no default onto the yield surface). F(sigma0, pi_i0) is computed at
// construction: a start OUTSIDE the surface (F > 1e-6 |p0|; |F| <= 1e-6 |p0| is ON the surface, warned) is refused; an on-surface start
// is accepted with a warning (the first loading step is then plastic).
//
// See Ladruno_implementation/144_ladruno_norsand_plan.md and 144a_norsand_equation_sheet.md.
// classTags 33023 / 33024 / 33025. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoNorSand_h
#define LadrunoNorSand_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <classTags.h>
#include "LadrunoNorSandKernel.h"

class LadrunoNorSand : public NDMaterial {
 public:
  // dimensional views (element-facing ordering, engineering shear)
  enum { DIM_3D = 0,        // {11,22,33,12,23,13}  order 6
         DIM_PSTRAIN = 1 }; // {11,22,12}  (eps_33 = 0)  order 3

  // sendSelf / recvSelf wire layout: ONE Vector of WIRE_LEN doubles (offsets in LadrunoNorSand.cpp,
  // namespace nsw). Test-friendly: tests/ci can rebuild a material from the raw vector. The base
  // NDMaterial sends nothing, so this is the only Vector under the dbTag (FE_Datastore keys vectors by
  // size: nothing to collide with).
  static const int WIRE_LEN = 128;

  // null constructor (broker / recvSelf): parameters are filled by recvSelf
  LadrunoNorSand();

  // full constructor. sig0 is the initial stress (OpenSees order 11,22,33,12,23,13; tension
  // positive = kernel order and sign), v0 the initial specific volume (a separate committed state
  // variable from v), pi0 the initial image pressure (REQUIRED, < 0; NaN is refused: initOK() false).
  LadrunoNorSand(int tag, const ladruno_norsand::Params& p, const double sig0[6],
                 double v0, double pi0, double dens = 0.0);

  ~LadrunoNorSand();

  const char* getClassType(void) const { return "LadrunoNorSand"; }

  int setTrialStrain(const Vector& strain);
  int setTrialStrain(const Vector& v, const Vector& r);
  int setTrialStrainIncr(const Vector& v);
  int setTrialStrainIncr(const Vector& v, const Vector& r);

  const Matrix& getTangent(void);        // NON-symmetric consistent tangent of the last step
  const Matrix& getInitialTangent(void); // hyperelastic tangent at the INITIAL state
  const Vector& getStress(void);
  const Vector& getStrain(void);
  double getStressZZ(void);              // plane-strain sigma_zz (NaN in other dims)

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);               // back to the INITIAL state (sigma0, v0, pi0), clears the latch

  NDMaterial* getCopy(void);
  NDMaterial* getCopy(const char* type);
  const char* getType(void) const;
  int getOrder(void) const;
  double getRho(void) { return density; }

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

  // ---- construction diagnostics (the parser) ----
  bool initOK(void) const { return initOk; }
  double initF0(void) const { return F0init; }   // F(sigma0, pi_i0) of the initial state (NaN if not built)
  const char* initMessage(void) const { return initMsg.c_str(); }
  void echoParameters(OPS_Stream& s) const;      // parameter echo at construction (like LadrunoSANISAND)

 protected:
  // derived-class constructors: fixed classTag and dimension
  LadrunoNorSand(int clsTag, int dimMode);
  LadrunoNorSand(int tag, int clsTag, const ladruno_norsand::Params& p, const double sig0[6],
                 double v0, double pi0, double dens, int dimMode);
  // copy EVERYTHING (parameters, both states, strains, latch, counters) from another instance;
  // the clone shares no storage with the source (history isolation). The dbTag is NOT copied.
  void copyFrom(const LadrunoNorSand& o);

 private:
  // parameters and initial-state inputs
  ladruno_norsand::Params kp;
  double density;                  // mass density (the element mass; NOT the ellipticity rho)
  double sigma0[6];                // initial stress, kernel order
  double v0init;                   // initial specific volume
  double pi0init;                  // initial image pressure (required, < 0)
  double F0init;                   // F(sigma0, pi_i0) at construction (derived; NaN until built)
  bool   initOk;
  std::string initMsg;

  // dimensional view
  int dim;
  int ncomp;
  int vmap[6];

  // kernel states: initial, committed, trial
  ladruno_norsand::State s0, sC, sT;
  double epsC[6], epsT[6];         // total strain, TENSOR components, kernel order
  double sigT[6];                  // trial stress (kernel order = OpenSees Voigt order)
  double CT[6][6];                 // trial consistent tangent, kernel convention (see header note)

  // refusal / latch / census (per instance)
  bool trialRefused;               // the last setTrialStrain was refused (cleared by a good trial / revert)
  bool latched;                    // a refused trial reached commitState: refuse everything until revertToStart
  bool warnedTrial, warnedLatch;
  int  lastRefusal;                // kernel Refusal code of the last refused trial; -1 = non-finite input/output
  int  lastFinest, lastFinestSub;  // finest-level cause (StepInfo.finest / finest_sub) of that refusal; 0 if none
  int  nRefusals;                  // trials refused since revertToStart
  int  nSubstepped;                // COMMITTED steps that needed substepping (counted in commitState)
  ladruno_norsand::StepInfo lastInfo;

  // helpers
  void   setupDim(void);
  int    integrate(void);
  void   restoreTrialFromCommitted(void);
  void   buildInitialState(void);

  // element-facing return buffers (sized to ncomp)
  Vector stressOut;
  Vector strainOut;
  Matrix tangentOut;
};

#endif
