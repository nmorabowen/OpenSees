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
// Description: SolidTransformationFinite — the FINITE-STRAIN geometry method
// (v3), updated-Lagrangian / spatial multiplicative.
//
// Unlike linear and corot (whose material core sees a SMALL strain), finite
// declares getStrainMeasure() == DeformationGradient: it tells the consuming
// element to compute the full deformation gradient F per Gauss point and drive
// a FiniteStrainNDMaterial via setTrialF(F) (NOT setTrialStrain), receiving the
// Cauchy stress σ and the spatial CONSTITUTIVE modulus c.
//
// DIVISION OF LABOUR (per the LOCKED seam-3 tangent contract, see
// FiniteStrainNDMaterial.h). The element forms the complete consistent tangent
//
//        K = ∫ Bᵀ c B dv  +  ∫ Gᵀ Σ G dv
//
// where the geometric / initial-stress term ∫ Gᵀ Σ G dv is built from the Cauchy
// stress σ and the element's own (spatial) shape-function gradients G — data the
// transformation does not have. Therefore, for finite, the seam methods are all
// IDENTITY pass-throughs: the element computes F, pushes its gradients to the
// current configuration, and assembles BOTH the material and geometric parts in
// its own Gauss loop. The transformation's sole finite-specific job is to flip
// the strain measure so one element supports linear / corot / finite by swapping
// theGeom — no element rewrite.
//
//   getStrainMeasure() == DeformationGradient   (the F-path signal)
//   update()            — no-op (the element computes F per GP from ∇N)
//   localizeDisp()      — identity (F is built from the global displacements)
//   globalizeForce()    — identity (∫ Bᵀ σ dv is already in global coords)
//   globalizeStiff()    — identity (the element assembles K_mat + K_geo itself)
//   frameTimeVarying()  == false   (Cauchy σ is reported in the global frame)

#ifndef SolidTransformationFinite_h
#define SolidTransformationFinite_h

#include <SolidTransformation.h>

class SolidTransformationFinite : public SolidTransformation
{
 public:
  SolidTransformationFinite();
  virtual ~SolidTransformationFinite();

  int getMethodID() const { return METHOD_FINITE; }
  const char *getName() const { return "finite"; }
  SolidTransformation *getCopy() const;

  StrainMeasure getStrainMeasure() const { return StrainMeasure::DeformationGradient; }

  int update(int numNodes, const Matrix &refCrds, const Matrix &curCrds);

  int localizeDisp(const Vector &uGlobal, Vector &uCore) const;

  int globalizeForce(const Vector &fCore, Vector &fGlobal) const;
  int globalizeStiff(const Matrix &kCore, const Vector &fCore,
                     Matrix &kGlobal) const;

  bool frameTimeVarying() const { return false; }
  void getCurrentFrame(Matrix &R) const;
};

#endif
