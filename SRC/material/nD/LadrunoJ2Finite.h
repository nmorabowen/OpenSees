/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// Description: LadrunoJ2Finite — a finite-strain-NATIVE combined-hardening von
// Mises (J2) material (de Souza Neto, Perić & Owen 2008, §14.11). It is a
// FiniteStrainNDMaterial subclass driven by setTrialF(F) that OWNS its
// finite-strain state and reuses the SAME verified small-strain return map as
// nDMaterial LadrunoJ2 (LadrunoJ2Kernel.h) — the kernel is frame-agnostic, so the
// only genuinely new ingredient is the spatial CO-ROTATION of the kinematic
// backstress.
//
// WHY (the §14.11 boundary this crosses): the v1 path
// LadrunoBrick -geom finite → LogStrain → LadrunoJ2 lifts the UNCHANGED
// small-strain material, which stores the Chaboche backstress α in a FIXED frame.
// The log-strain (MATISU) algorithm co-rotates the elastic state automatically
// (bᵉᵗʳ = f_Δ bᵉ_n f_Δᵀ, so the deviatoric stress s emerges in the CURRENT frame),
// but α does NOT co-rotate ⇒ ‖s − α‖ subtracts two tensors in different frames ⇒
// combined hardening is NOT frame-indifferent under large rotation (pinned as a
// strict xfail on the v1 wrapper). Isotropic hardening, which depends only on ‖s‖
// and the scalar ε̄ᵖ, is unaffected.
//
// THE FIX: own the backstresses α_{n,k} as SPATIAL tensors and push them forward
// by the incremental material rotation each step BEFORE the return map:
//
//     α̃_{n,k} = R_Δ · α_{n,k} · R_Δᵀ ,   R_Δ = polar(f_Δ),   f_Δ = F·F_n⁻¹
//
// Then run the UNCHANGED return map in the rotated frame. For a pure rigid
// increment f_Δ = Q this gives α̃ = Q α Qᵀ, exactly cancelling the rotation of s;
// as f_Δ → I, R_Δ → I, so it reduces to the small-strain / non-rotating case
// exactly. Verified objective to ~1e-12 in tests/ladrunoj2_finite_native_reference.py.
//
// CONSISTENT TANGENT (two channels — see Ladruno_implementation/16_finite_native_j2_adr.md):
//   (A) the trial elastic strain εᵉᵗʳ — the standard log-strain spatial tangent
//       c = (1/2J)[D:L:B] from LogStrainKernel (analytic);
//   (B) the co-rotated backstress α̃(R(F)) — a small (~1e-4..2e-3 relative) but
//       real, frame-objective constitutive sensitivity that straddles the
//       tangent-FD gate. It is recovered EXACTLY-to-FD by perturbing ONLY R
//       (≈9 cheap return-map calls, plastic steps only), then mapped into the
//       element's cmat via Ḟ = L F:  cmatB_ijkl = Σ_m (∂σ_ij/∂F_km) F_lm. No
//       kernel re-derivation. (Validated in tests/test_ladrunoJ2_finite_native_tangent.py.)
//
// State (per Gauss point): F_n, bᵉ_n, ε̄ᵖ_n, α_{n,k} (SPATIAL). Reuses
// LadrunoJ2Kernel.h (return map), LogStrainKernel.h (Hencky kinematics + spatial
// tangent + eigenvalue-degeneracy branches), LadrunoHardening.h (shared σ_y).
//
// Companion: nDMaterial LadrunoJ2 (small-strain, classTag 33011);
// element LadrunoBrick -geom finite (classTag 33002).

#ifndef LadrunoJ2Finite_h
#define LadrunoJ2Finite_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LadrunoJ2Finite : public FiniteStrainNDMaterial
{
 public:
  static const int MAXBACK = 8;   // max Chaboche backstress terms (matches kernel)

  LadrunoJ2Finite(int tag, double K, double G,
                  double sig0, double Qinf, double bIso, double Hiso,
                  int nBack, const double *C, const double *gam, double rho = 0.0);
  LadrunoJ2Finite();
  ~LadrunoJ2Finite();

  const char *getClassType(void) const override { return "LadrunoJ2Finite"; }

  // --- the finite-strain seam (FiniteStrainNDMaterial contract) ---
  // `override` on every contract method is deliberate: a signature drift would be
  // a COMPILE error rather than a silent fallback to the base default (notably
  // getSpatialTangentTensor's base returns -1 → the element would fail loudly, but
  // a silent shadow on getStress/getTangent could corrupt results).
  int setTrialF(const Matrix &F) override;             // total deformation gradient
  const Vector &getStress(void) override;              // Cauchy σ (6-vector)
  const Matrix &getTangent(void) override;             // spatial constitutive c (6×6, lossy)
  int getSpatialTangentTensor(double c[3][3][3][3]) override;  // FULL c_ijkl (element channel)
  const Matrix &getInitialTangent(void) override;
  double getJ(void) const override { return Jdet; }

  const Vector &getStrain(void) override;              // Hencky strain εᵉᵗʳ (6-vector)
  double getRho(void) override { return rho; }

  int commitState(void) override;
  int revertToLastCommit(void) override;
  int revertToStart(void) override;

  NDMaterial *getCopy(void) override;
  NDMaterial *getCopy(const char *type) override;
  const char *getType(void) const override { return "ThreeDimensional"; }
  int getOrder(void) const override { return 6; }

  int sendSelf(int commitTag, Channel &theChannel) override;
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

  void Print(OPS_Stream &s, int flag = 0) override;

  int setParameter(const char **argv, int argc, Parameter &param) override;
  int updateParameter(int parameterID, Information &info) override;
  Response *setResponse(const char **argv, int argc, OPS_Stream &output) override;
  int getResponse(int responseID, Information &matInfo) override;

 private:
  // material constants
  double bulk, shear;                 // K, G
  double sig0, Qinf, bIso, Hiso;      // Voce+linear isotropic hardening
  double rho;                         // mass density
  int nBack;                          // number of Chaboche backstress terms
  double Ckin[MAXBACK], gKin[MAXBACK];

  // committed finite-strain state (per Gauss point)
  double Fn[9];                       // committed total deformation gradient (row-major)
  double Be_n[9];                     // committed elastic left Cauchy–Green bᵉ_n
  double ebarP_n;                     // committed accumulated equiv. plastic strain
  double alpha_n[MAXBACK][6];         // committed SPATIAL backstress (tensor comps)

  // trial state staged for commit
  double Ftrial9[9];                  // trial F
  double BeTrial9[9];                 // trial bᵉᵗʳ (used for the §14.5 spatial tangent)
  double Be_trialUpd[9];             // bᵉ to become bᵉ_n on commit = exp[2 εᵉ_{n+1}]
  double ebarP_trial;
  double alpha_trial[MAXBACK][6];    // updated backstress (current frame), staged
  double dGammaTrial;                // plastic multiplier increment (record / plastic flag)

  // cached trial outputs
  Vector sigmaCauchy;                // Cauchy stress (6)
  Vector henckyStrain;               // εᵉᵗʳ (6, engineering shear)
  Matrix aTangent;                   // spatial constitutive 6×6 (channel A, lossy)
  Matrix K0init;                     // per-instance elastic 6×6 for getInitialTangent
  double Jdet;                       // det F
  double DtanTrial[6][6];            // kernel algorithmic tangent at the converged dG

  void setIdentity(double M[9]);
};

#endif
