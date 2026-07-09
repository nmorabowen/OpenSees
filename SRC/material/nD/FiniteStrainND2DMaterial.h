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
// Created: 07/2026
//
// Description: Abstract base class for 2D (PLANE) finite-strain material adaptors
// — the plane sibling of FiniteStrainNDMaterial (seam 3). It is the contract
// between an F-delivering plane finite-strain element (LadrunoQuad / LadrunoCST
// with -geom finite, ADR 25 P5) and an UNCHANGED small-strain material.
//
// This class adds ONE thing to the NDMaterial interface: a 2×2 deformation-
// gradient entry point setTrialF(F). A concrete subclass (LogStrain2D) presents
// an order-3 PlaneStrain / PlaneStress face over a 3D small-strain inner material
// via the same spatial multiplicative (logarithmic / Hencky) technique as the 3D
// adaptor — de Souza Neto, Perić & Owen (2008), Box 14.3 (MATISU), §14.7 (plane
// stress via thickness-stretch condensation).
//
// CONTRACT (locked for the ADR 25 P5 element wiring):
//
//   * The ELEMENT calls setTrialF(F) with the TOTAL 2×2 deformation gradient
//     F = I + ∂u/∂X (in-plane) at the Gauss point. The adaptor lifts F to 3×3
//     internally (F₃₃ = 1 for plane strain; F₃₃ solved from σ₃₃ = 0 for plane
//     stress) and NEVER has setTrialStrain(...) called on it. Returns < 0 on
//     failure (det F ≤ 0, non-convergence) so the element can cut the step.
//
//   * getStress() returns the in-plane CAUCHY stress σ as a 3-vector in the plane
//     Voigt order {11, 22, 12} with ENGINEERING shear (γ₁₂ = 2 ε₁₂). The element
//     assembles the SPATIAL internal force ∫ bᵀ σ t dA on the current config.
//
//   * The spatial CONSTITUTIVE modulus is c_ijkl = the in-plane block that the
//     element combines with the geometric term:  a_ijkl = c_ijkl − σ_il δ_jk
//     (i,j,k,l ∈ {0,1}), exactly as in 3D. As in 3D, c is NOT minor-symmetric in
//     (k,l), so getTangent() (3×3) is a LOSSY projection for output only; a
//     consistent element tangent MUST use getSpatialTangentTensor2D(c). For plane
//     stress, c is the STATICALLY CONDENSED modulus (the out-of-plane direction is
//     eliminated enforcing δσ₃₃ = 0), so a_ijkl = c_ijkl − σ_il δ_jk still gives
//     the correct consistent plane-stress tangent.
//
//   * getType() reports "PlaneStrain" or "PlaneStress"; getOrder() returns 3.
//
//   * STATE OWNERSHIP: the concrete adaptor owns the (3D) committed elastic left
//     Cauchy–Green state and, for plane stress, the committed thickness stretch
//     λ = F₃₃ (warm-starts the next Newton). It wraps the inner material's
//     commitState / revertToLastCommit / revertToStart / sendSelf / recvSelf.
//
// See Ladruno_implementation/25_ladruno_plane_elements_adr.md (P5) and
// FiniteStrainNDMaterial.h for the 3D contract this parallels.

#ifndef FiniteStrainND2DMaterial_h
#define FiniteStrainND2DMaterial_h

#include <NDMaterial.h>
#include <OPS_Globals.h>

class Matrix;
class Vector;

class FiniteStrainND2DMaterial : public NDMaterial
{
 public:
  FiniteStrainND2DMaterial(int tag, int classTag) : NDMaterial(tag, classTag) {}
  FiniteStrainND2DMaterial() : NDMaterial() {}
  virtual ~FiniteStrainND2DMaterial() {}

  // ===================================================================== //
  //  THE PLANE FINITE-STRAIN SEAM (seam 3, 2D)                            //
  // ===================================================================== //

  // Set the trial state from the TOTAL 2×2 deformation gradient F at the GP.
  // Returns 0 on success, < 0 on failure (det F ≤ 0, plane-stress non-convergence).
  virtual int setTrialF(const Matrix &F) = 0;

  // ===================================================================== //
  //  SPATIAL RESPONSE — what the plane element consumes                   //
  // ===================================================================== //

  // In-plane CAUCHY stress σ (3-vector, Voigt {11,22,12}, engineering shear).
  virtual const Vector &getStress(void) = 0;
  // In-plane spatial CONSTITUTIVE modulus as a 3×3 — the LOSSY minor-symmetric
  // projection (c is not (k,l)-symmetric). Output only; a consistent element
  // tangent MUST use getSpatialTangentTensor2D() below.
  virtual const Matrix &getTangent(void) = 0;

  // FULL 4th-order in-plane spatial constitutive modulus c_ijkl (i,j,k,l ∈ {0,1})
  // — the element's tangent channel. The element forms a_ijkl = c_ijkl − σ_il δ_jk
  // and assembles K = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) t dA. For plane stress this
  // c is already statically condensed (δσ₃₃ = 0). Default −1 = unsupported.
  virtual int getSpatialTangentTensor2D(double c[2][2][2][2]) { return -1; }

  // In-plane Jacobian J = det(2×2 F) (> 0), the area change the element uses for
  // the current-configuration mapping. The returned σ is already true Cauchy.
  virtual double getJ(void) const = 0;

  // ===================================================================== //
  //  The small-strain entry point is DISABLED for finite-strain materials //
  // ===================================================================== //
  virtual int setTrialStrain(const Vector & /*v*/) override {
    opserr << "FiniteStrainND2DMaterial: setTrialStrain() is invalid — this is a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
  virtual int setTrialStrain(const Vector & /*v*/, const Vector & /*r*/) override {
    opserr << "FiniteStrainND2DMaterial: setTrialStrain() is invalid — this is a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }

 protected:
 private:
};

#endif
