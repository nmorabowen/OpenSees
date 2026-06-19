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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 06/2026
//
// Description: LadrunoRCFiniteStrain — a finite-strain-NATIVE view of the RC
// plastic-damage material (ADR 19_ladruno_rc_shell_adr.md, Phase 4b). It is a
// FiniteStrainNDMaterial subclass driven by setTrialF(F) that reuses the SAME
// verified small-strain return map as nDMaterial LadrunoRCConcrete
// (LadrunoRCKernel.h) — the kernel is frame-agnostic, so the only new ingredient
// is the Hencky (logarithmic-strain) seam.
//
// WHY a NATIVE subclass and not the generic LogStrain wrapper: the RC spine is a
// plastic-DAMAGE law (stiffness degrades: tau = (1-d) Dᵉ:εᵉ), and the generic
// LogStrainNDMaterial commits its elastic left Cauchy–Green by recovering the
// elastic strain through the inner's INITIAL tangent (εᵉ = Cᵉ:tau, exact only for
// a non-degrading inner). For a damage inner that recovery is shrunk by (1-d), so
// the committed bᵉ drifts. This native view avoids it: the RC kernel carries NO
// tensorial plastic strain (its inelasticity is the damage thresholds xt/xc + the
// effective-stress recompose, returnMap3D is a pure function of the TOTAL strain),
// so the multiplicative split is trivial — the elastic left Cauchy–Green is just
// the TOTAL B = F Fᵀ, recomputed each step from the element's F; there is no bᵉ to
// track or mis-commit. (returnMap3D uses an incremental elastic predictor from the
// committed stress_eff + committed strain, so the committed Hencky strain stored in
// RCHist IS the de-facto committed elastic reference — and it is serialized. There is
// simply no SEPARATE bᵉ/Fn beyond RCHist.) The committed state is exactly
// LadrunoRCConcrete's RCHist; the finite view adds no committed state of its own.
//
// THE SEAM (de Souza Neto, Perić & Owen 2008, Box 14.3, §14.5):
//   B   = F Fᵀ                 (total left Cauchy–Green; J = det F)
//   εᵉ  = ½ ln B               (Hencky strain, engineering-shear Voigt)
//   τ,D = returnMap3D(εᵉ)      (nominal/Kirchhoff stress + small-strain tangent)
//   σ   = τ / J                (Cauchy stress)
//   c   = (1/2J)[D:L:B]        (spatial constitutive modulus; LogStrainKernel)
//
// OBJECTIVITY (§14.11). The ISOTROPIC-damage spine is objective by construction
// (B = F Fᵀ co-rotates; the spectral split of B and the scalar damage thresholds
// are frame-indifferent). The fixed-crack normal + aggregate-interlock direction
// stored in RCHist are DIRECTIONAL internal variables in a FIXED frame — this
// setTrialF material view does NOT co-rotate them, so under LARGE rotation a
// cracked point is NOT frame-indifferent (pinned as a strict xfail, the same
// boundary as LadrunoJ2Finite's kinematic backstress before its v2 co-rotation).
// The supported objective large-rotation route is the corotational ELEMENT (it
// feeds the de-rotated strain), proven for the small-strain material. A
// co-rotating-crack finite-native view (to flip the xfail) is deferred.
//
// IMPL-EX (`-implex`): carried through unchanged from LadrunoRCConcrete — the
// EXPLICIT extrapolation freezes the damage thresholds in setTrialF (secant
// tangent), commitState re-integrates IMPLICITLY to advance the true thresholds
// and measure implexError. Default off ⇒ implicit.
//
// State (per Gauss point): the RC Params + RCHist (no finite-specific committed
// state — B is recomputed from the element's F). Reuses LadrunoRCKernel.h (return
// map) + LogStrainKernel.h (Hencky kinematics + spatial tangent).
//
// Companion: nDMaterial LadrunoRCConcrete (small-strain, classTag 33015);
// element LadrunoBrick -geom finite (classTag 33002). classTag 33018.

#ifndef LadrunoRCFiniteStrain_h
#define LadrunoRCFiniteStrain_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include "LadrunoRCKernel.h"

class LadrunoRCFiniteStrain : public FiniteStrainNDMaterial {
 public:
  LadrunoRCFiniteStrain();
  LadrunoRCFiniteStrain(int tag, const ladruno_rc_kernel::Params& P, double rho);
  ~LadrunoRCFiniteStrain();

  const char* getClassType(void) const override { return "LadrunoRCFiniteStrain"; }

  // --- the finite-strain seam (FiniteStrainNDMaterial contract) ---
  int setTrialF(const Matrix& F) override;                     // total deformation gradient
  const Vector& getStress(void) override;                      // Cauchy σ (6)
  const Matrix& getTangent(void) override;                     // spatial constitutive c (6×6, lossy)
  int getSpatialTangentTensor(double c[3][3][3][3]) override;  // FULL c_ijkl (element channel)
  const Matrix& getInitialTangent(void) override;
  double getJ(void) const override { return Jdet; }

  const Vector& getStrain(void) override;                      // Hencky strain εᵉ (6)
  double getRho(void) override { return rho; }

  int commitState(void) override;
  int revertToLastCommit(void) override;
  int revertToStart(void) override;

  NDMaterial* getCopy(void) override;
  NDMaterial* getCopy(const char* type) override;
  const char* getType(void) const override { return "ThreeDimensional"; }
  int getOrder(void) const override { return 6; }

  int sendSelf(int commitTag, Channel& theChannel) override;
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) override;
  void Print(OPS_Stream& s, int flag = 0) override;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s) override;
  int getResponse(int responseID, Information& matInfo) override;

 private:
  ladruno_rc_kernel::Params P;       // material constants + backbones (kernel-side)
  ladruno_rc_kernel::RCHist histN;   // committed history
  ladruno_rc_kernel::RCHist histTr;  // trial history (overwritten each integrate)
  double rho;

  // finite-strain working state (per Gauss point)
  double Ftrial9[9];                 // trial total F (row-major)
  double Btotal9[9];                 // B = F Fᵀ (kept for the §14.5 spatial tangent)
  double Jdet;                       // det F
  double strain6[6];                 // trial Hencky strain εᵉ (engineering shear)
  double stress6[6];                 // kernel nominal stress = Kirchhoff τ (engineering)
  double Dtan[6][6];                 // kernel small-strain algorithmic tangent ∂τ/∂εᵉ
  int    status;

  // IMPL-EX runtime (config lives in P.implex* so getCopy propagates it)
  double implexError;
  double dtime_n, dtime_n_commit, dtime_0;
  bool   commitDone;

  // Phase 3b: crack-band regularization latch (config in P.autoReg/P.lchRef; getCopy propagates)
  bool   regularizationDone;
  double regLch;
  bool   regWarned;
  int    regularizeIfNeeded(void);

  // helpers
  void   integrate(bool do_implex, double tfac);   // Hencky -> returnMap3D -> stress6/Dtan
  double implexTimeFactor(void) const;

  // element-facing return buffers
  Vector sigmaCauchy;                // Cauchy σ (6)
  Vector henckyStrain;               // εᵉ (6, engineering shear)
  Matrix aTangent;                   // spatial constitutive 6×6 (lossy)
  Matrix tangentOut;                 // initial-tangent buffer
};

#endif
