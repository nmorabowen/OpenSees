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
// Description: SolidTransformationCorot — the COROTATIONAL geometry method
// (v2). Large rigid-body ROTATION + large displacement, but SMALL material
// strain. It strips the element rigid rotation R so the constitutive core sees
// only small strain in a co-rotating frame, then rotates the core force /
// stiffness back to global and adds the geometric (load-) stiffness implied by
// the configuration-dependence of R. It REUSES the full small-strain material
// library (setTrialStrain) — no finite-strain material is involved.
//
//   getStrainMeasure() == SmallStrain   (like linear; the element runs its
//                                        normal std/bbar/uri/eas B-matrix path
//                                        on the de-rotated displacement)
//
// THEORY (see Ladruno_implementation/10_solid_corotational_adr.md and the two
// independent derivations folded into it). For a continuum element with
// TRANSLATIONAL DOFs only:
//
//   * Rotation R = polar(H), Procrustes best fit, from the nodal cloud:
//         H = Σ_a (x_a − x_c)(X_a − X_c)ᵀ ,   R = H S⁻¹ ,  S = (HᵀH)^½.
//     Chosen over polar-of-volume-averaged-F because H needs only nodal
//     coordinates — exactly what the seam update() is handed — whereas F̄ would
//     need shape-function gradients the wrapper never receives.
//
//   * Deformational (corotated) displacement, fed to the small-strain core:
//         u_d,a = Rᵀ (x_a − x_c) − (X_a − X_c).
//     Under any rigid body motion x_a = Q X_a + c the polar gives R = Q exactly,
//     so u_d ≡ 0 ⇒ ε = 0 ⇒ σ = 0 ⇒ f = 0. The rigid-rotation / zero-stress gate
//     is therefore an ALGEBRAIC IDENTITY (its only failure mode is a non-exact
//     R — so the polar iteration converges to ~machine precision). No EICR
//     rigid-body projector is needed: centroid-centering removes the 3 rigid
//     translations and the Rᵀ pull-back removes the 3 rigid rotations.
//
//   * Spin map Γ = ∂θ/∂u (3 × 3·numNodes), the SPATIAL rotation sensitivity
//     (δR = spin(Γ δu) R), used for the geometric stiffness. Compute the body
//     spin first, then push to spatial:
//         ω_body = (tr(S) I − S)⁻¹ · 2 · axial( skew( Rᵀ δH ) ),  δH = e_k ⊗ X_b⁰,
//         Γ_{:,col} = R · ω_body.
//     (Verified against the U=I limit AND a Γ-vs-finite-difference gate — the
//     R· push to the spatial frame is essential; omitting it puts the geometric
//     terms in the wrong frame.)
//
//   * Internal force (EXACT gradient of the corotated energy):
//         f_global,b = R f_d,b − (1/n) R Σ_a f_d,a − Γ_{:,b}ᵀ m ,
//         m = Σ_a (x_a−x_c) × (R f_d,a).
//     The −(1/n)RΣf_d centroid back-reaction vanishes for the self-equilibrated
//     internal force but is kept for the body/applied-load case.
//
//   * Tangent (v2.0): K = R̄ k_d R̄ᵀ + K_geo , R̄ = blockdiag(R), and the
//     geometric term assembled symmetrically as G1 + G1ᵀ with G1[3b+i,c] =
//     −[spin(R f_d,b) Γ]_{i,c} (the G2 lever-arm term equals G1ᵀ at equilibrium,
//     so it is INCLUDED). TWO higher-order geometric pieces remain DEFERRED:
//       (a) the polar-Hessian   ∂Γ/∂u · m  (curvature of the rotation map), and
//       (b) the material-frame coupling  R k_d (∂Rᵀ/∂u) x⁰  (the k_d-bearing
//           part of ∂f_d/∂u through u_d = Rᵀx⁰ − X⁰).
//     Both are O(strain·force)/O(∂Γ·force) — they vanish at small deformation
//     and are exactly what the (looser) FD-tangent gate tolerates; both must be
//     added to reach quadratic Newton in v2.1. Because the RESIDUAL is the exact
//     gradient, Newton converges to the correct equilibrium regardless — only
//     the asymptotic rate is affected. See the ADR §5/§7.
//
// STATELESS. R is a pure function of the current nodal positions (continuum
// nodes carry no rotational DOFs, so there is no committed triad to advance —
// unlike CorotCrdTransf3d). commitState/revert*/serialization therefore need no
// corot-specific state; the base no-ops suffice and the element rebuilds the
// method from its id alone.

#ifndef SolidTransformationCorot_h
#define SolidTransformationCorot_h

#include <SolidTransformation.h>

class SolidTransformationCorot : public SolidTransformation
{
 public:
  SolidTransformationCorot();
  virtual ~SolidTransformationCorot();

  int getMethodID() const { return METHOD_COROT; }
  const char *getName() const { return "corot"; }
  SolidTransformation *getCopy() const;

  StrainMeasure getStrainMeasure() const { return StrainMeasure::SmallStrain; }

  int update(int numNodes, const Matrix &refCrds, const Matrix &curCrds);

  int localizeDisp(const Vector &uGlobal, Vector &uCore) const;

  int globalizeForce(const Vector &fCore, Vector &fGlobal) const;
  int globalizeStiff(const Matrix &kCore, const Vector &fCore,
                     Matrix &kGlobal) const;

  bool frameTimeVarying() const { return true; }
  void getCurrentFrame(Matrix &R) const;

 private:
  // -- maximum nodes supported. Consumers: LadrunoBrick (8-node hex) and
  //    BezierTet10 (10-node quadratic tet). Guarded in update(). All corot
  //    logic runs on the runtime numNodes/ndof; MAX_NODES only sizes the
  //    fixed scratch arrays (Xrel/xrel, Gamma, the Rf buffer in
  //    globalizeForce), so raising it is allocation-only — no logic change
  //    and the existing 8-node consumer is unaffected.  // Ladruno (Tet10 corot)
  enum { MAX_NODES = 10, NDF = 3, MAX_DOF = MAX_NODES * NDF };

  int numNodes;                 // set by update(); 0 until first update()
  double Rmat[9];               // current rotation R (row-major 3x3)
  double xrel[MAX_NODES][3];    // current centred coords x_a − x_c
  double Xrel[MAX_NODES][3];    // reference centred coords X_a − X_c
  double Gamma[3][MAX_DOF];     // spin map ∂θ/∂u (3 × 3·numNodes)
};

#endif
