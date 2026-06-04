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
// Description: InitDefGradNDMaterial — the MULTIPLICATIVE, finite-strain sibling
// of InitStrainNDMaterial / InitStressNDMaterial. It is the "seam-2" staged-
// activation wrapper from Ladruno_implementation/staged_deformation_gradiend.md:
// it makes a finite-strain continuum element BORN STRESS-FREE at the current
// deformed geometry when it is appended mid-stage (staged construction — a new
// member, a concrete lift, a backfill layer).
//
// THE PROBLEM (why InitStrain is not enough). OpenSees stores displacement u from
// the original mesh X_ref, and that frame NEVER re-zeros. An element appended at
// u = u0 computes its strain from X_ref and is born carrying the full accumulated
// deformation as spurious stress. Subtracting u0 (additively, as InitStrainND-
// Material does) is EXACT for small strain but WRONG for finite strain: it is
// non-objective (a post-birth rigid rotation Q injects fake strain) and it adds in
// displacement/strain space where finite kinematics compose MULTIPLICATIVELY.
// Proof: subtracting u0 builds F_naive = F − F0 + I; under a pure post-birth
// rotation F → Q·F0 this gives Q·F0 − F0 + I ≠ I for Q ≠ I → spurious stress.
//
// THE FIX (multiplicative split, stays inside the one X_ref frame). At activation
// the wrapper captures the birth deformation gradient
//
//     F0 = F at first setTrialF   ( = I + Grad_X(u0), the committed deformation )
//
// and thereafter feeds the inner finite-strain material the RELATIVE gradient
//
//     F_rel = F · F0⁻¹            ( deformation measured from the stress-free birth )
//
// This is the multiplicative split F = F_rel · F0 (Lee), identical in form to
// F = Fe·Fp in finite plasticity or F = Fe·Fg in growth. Under a post-birth
// rotation F → Q·F: F_rel → Q·F_rel, so C_rel = F_relᵀF_rel is UNCHANGED → zero
// spurious stress. It degenerates correctly: appended at t = 0 ⇒ F0 = I ⇒
// F_rel = F (ordinary behaviour, zero overhead) — the finite-strain generalization
// of Truss's initialDisp-folded-into-L (a 1-D F0).
//
// WHY THE SPATIAL SEAM MAKES THIS A THIN WRAPPER. The shipped finite-strain seam
// (FiniteStrainNDMaterial: setTrialF(F) in, CAUCHY σ out, element integrates
// ∫bᵀσ dv on the current config with its own det F) is SPATIAL, not PK1/total-
// Lagrangian. Two facts collapse the F0 bookkeeping the (PK1-dialect) design note
// described:
//   1. The incremental gradient is invariant to F0:
//        F_Δ = F_rel · F_rel,n⁻¹ = (F·F0⁻¹)(F_n·F0⁻¹)⁻¹ = F · F_n⁻¹.
//      So F0 changes ONLY the stress-free initial state, never the step response.
//   2. Cauchy stress is reference-independent: the inner returns σ = τ_rel/J_rel,
//      which IS the true Cauchy stress. No P = P_rel·F0⁻ᵀ push-back and no J0
//      integration-weight are needed — the element already integrates on the
//      current configuration with its own det F. (The inner's getJ() returns
//      J_rel = det F / det F0; the element does NOT use it for assembly.)
// So the wrapper is a near pass-through: capture F0 once, forward F_rel, and
// delegate stress / tangent / state to the inner FiniteStrainNDMaterial.
//
// CONTRACT. Wraps an inner FiniteStrainNDMaterial (e.g. LogStrainNDMaterial). IS-A
// FiniteStrainNDMaterial, so LadrunoBrick -geom finite (and any future F-based
// solid) consumes it with ZERO element edits. getStress() = inner Cauchy σ (the
// true stress, relative to birth = relative to current, identical object);
// getStrain() = inner RELATIVE Hencky strain (from birth) by default; the captured
// F0 is exposed via setResponse("F0").
//
// See Ladruno_implementation/staged_deformation_gradiend.md (the design note) and
// constraints_reference_position.md (the SP/MP/element reference-frame trichotomy
// that motivates it).

#ifndef InitDefGradNDMaterial_h
#define InitDefGradNDMaterial_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class InitDefGradNDMaterial : public FiniteStrainNDMaterial
{
 public:
  // active=false  -> pure pass-through (F0 ≡ I, never re-references): the wrapper
  //                  is bit-identical to the bare inner (the -noInitF opt-out).
  // F0given!=0    -> explicit birth gradient (9 row-major components); skips the
  //                  auto-capture. Otherwise F0 is captured on the FIRST setTrialF.
  InitDefGradNDMaterial(int tag, NDMaterial &inner, bool active = true,
                        const double *F0given = 0);
  InitDefGradNDMaterial();
  ~InitDefGradNDMaterial();

  const char *getClassType(void) const { return "InitDefGradNDMaterial"; }

  // --- the finite-strain seam (FiniteStrainNDMaterial contract) ---
  int setTrialF(const Matrix &F);                     // total F in; forwards F·F0⁻¹
  const Vector &getStress(void);                      // Cauchy σ (delegated)
  const Matrix &getTangent(void);                     // spatial c 6×6 (delegated)
  int getSpatialTangentTensor(double c[3][3][3][3]);  // full c_ijkl (delegated)
  const Matrix &getInitialTangent(void);
  double getJ(void) const;                            // J_rel = det(F·F0⁻¹)

  const Vector &getStrain(void);                      // relative Hencky (delegated)
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const { return "ThreeDimensional"; }
  int getOrder(void) const { return 6; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

 private:
  FiniteStrainNDMaterial *theMaterial;  // inner finite-strain material (deep copy)

  bool   active;        // false ⇒ pass-through (no re-referencing)
  bool   captured;      // F0 has been snapshot (true once birth is seen / given)
  bool   f0Explicit;    // F0 supplied at construction (revertToStart keeps it;
                        // auto-capture mode re-arms instead)
  double F0[9];         // birth deformation gradient (row-major), captured at birth
  double invF0[9];      // F0⁻¹ (cached; rebuilt from F0 on construct / recvSelf)

  Matrix Frel;          // scratch 3×3 forwarded to the inner setTrialF
  Vector birthF0;       // setResponse("F0") channel (9, row-major)

  // 3×3 row-major helpers (self-contained; no external linear-algebra dep)
  static void   setIdentity9(double M[9]);
  static double det9(const double M[9]);
  static int    inv9(const double M[9], double Mi[9]);  // returns -1 if singular
  static void   mul9(const double A[9], const double B[9], double C[9]); // C = A·B
};

#endif
