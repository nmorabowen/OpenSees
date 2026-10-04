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

// LadrunoNorSand: NorSand in the Andrade & Borja 2006 (IJNME 67, 3-invariant) form, as an nDMaterial.
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
//     trial and commit is refused until revertToStart(). sendSelf/recvSelf keeps the latch (the same point,
//     restarted); getCopy does NOT: a clone starts with no latch and an empty refusal record (see copyFrom).
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
// Initial state: -pi0 is REQUIRED unless -pi0_auto is given (no silent default onto the yield surface).
// F(sigma0, pi_i0) is computed at construction, at the stress the state really has (the deck's sigma0, or
// sigma(eps^e_f) when the p' floor projected the initial state): a start OUTSIDE the surface (F > 1e-6 p_ref,
// p_ref = |p0| (BA06) or p_a (HAR); |F| <= 1e-6 p_ref is ON the surface, warned) is refused; an on-surface start
// is accepted with a warning (the first loading step is then plastic). -pi0_auto selects the unified
// pi_i0 rule of sheet 5.4 (S.53) (the surface through (p_init, max(eta_init, c2 M)), owner decision (d)).
//
// ENERGY OPTION and p' FLOOR (WP-144 round 3, sheet 2.3-2.4 and 9.7; the mathematics is the kernel's):
//   -energy BA06 (default, paper mode: -p0 -kappa_hat -mu0 [-eps_v0] [-alpha0]) | HAR (Houlsby-Amorosi-Rojas 2005:
//   -k -g [-n 0.5] -p_a, or the DM04 mapping -G0 -nu [-e_ref]); HAR REPLACES the BA06 constants, and giving one
//   of them (or a HAR one with BA06) is refused, never ignored. -p_a is ONE flag, shared by the HAR energy and the
//   fork CSL. -pmin p: the floor (default 5e-3 p_ref; 0 = off). The floor never refuses; it is COUNTED (responses
//   `floor` / `floored`: [at_floor, n_f_tr, n_f_post, eps_f_v, W_f]; `stepInfo` [9], [10] = floor_tr, floor_post of
//   the last step; `floorEnergy`: E_f). The consistent tangent of a floored step is the EXACT one (zero bulk
//   stiffness at a floored state; no regularisation, owner decision (c)).
//
// See Ladruno_implementation/144_ladruno_norsand_plan.md and 144a_norsand_equation_sheet.md.
// classTags 33023 / 33024 / 33025. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoNorSand_h
#define LadrunoNorSand_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <classTags.h>
#include "LadrunoElasticStrainProvider.h"   // LogStrain v2 elastic-strain mixin (WP-144 G2)
#include "LadrunoNorSandKernel.h"

class LadrunoNorSand : public NDMaterial, public LadrunoElasticStrainProvider {
 public:
  // dimensional views (element-facing ordering, engineering shear)
  enum { DIM_3D = 0,        // {11,22,33,12,23,13}  order 6
         DIM_PSTRAIN = 1 }; // {11,22,12}  (eps_33 = 0)  order 3

  // sendSelf / recvSelf wire layout: ONE Vector of WIRE_LEN doubles (offsets in LadrunoNorSand.cpp,
  // namespace nsw). Test-friendly: tests/ci can rebuild a material from the raw vector. The base
  // NDMaterial sends nothing, so this is the only Vector under the dbTag (FE_Datastore keys vectors by
  // size: nothing to collide with). The trial stress and tangent are NOT on the wire (Domain::recvSelf's
  // update() recomputes them); recvSelf rebuilds them from the restored trial State.
  static const int WIRE_LEN = 112;

  // null constructor (broker / recvSelf): parameters are filled by recvSelf
  LadrunoNorSand();

  // full constructor. sig0 is the initial stress (OpenSees order 11,22,33,12,23,13; tension
  // positive = kernel order and sign), v0 the initial specific volume (a separate committed state
  // variable from v; the kernel evolves v = v0 exp(tr eps), plan 2.8: the shell never computes v itself),
  // pi0 the initial image pressure (< 0). pi0 = NaN selects the UNIFIED rule of sheet 5.4 (S.53) (the parser's
  // -pi0_auto; the resolved value is then stored in place of the NaN). The parser itself refuses a deck with
  // neither -pi0 nor -pi0_auto; the kernel's "NaN = on the surface" apex rule is never reachable from a deck.
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

  // LadrunoElasticStrainProvider (WP-144 G2, owner decision 2026-10-01): the kernel TRIAL elastic strain
  // sT.eps_e (tensor, {00,11,22,01,12,02}) as engineering Voigt {11,22,33,12,23,13} (shear doubled), in
  // the material's own frame. Always the full 6 components, whatever the dimensional view. Returns false
  // if the point was never built (initOK() false) or the state is non-finite; epsE is then untouched.
  bool ladrunoGetElasticStrain(Vector& epsE) const;

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
  // copy the parameters, the initial/committed/trial states, the strains and the substep census from another
  // instance; the clone shares no storage with the source (history isolation). The dbTag is NOT copied.
  // The latch, the refused-trial flag and the refusal record are NOT copied: a clone starts un-refused.
  void copyFrom(const LadrunoNorSand& o);

 private:
  // parameters and initial-state inputs
  ladruno_norsand::Params kp;
  double density;                  // mass density (the element mass; NOT the ellipticity rho)
  double sigma0[6];                // initial stress, kernel order
  double v0init;                   // initial specific volume
  double pi0init;                  // initial image pressure (required, < 0)
  double F0init;                   // F(sigma0, pi_i0) at construction (derived; NaN until built)
  bool   pi0Auto;                  // pi_i0 came from the unified rule (S.53), not from the deck (echo only)
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
  double efC, efT;                 // cumulative floor energy E_f (S.52) over the committed history / incl. the trial step
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
