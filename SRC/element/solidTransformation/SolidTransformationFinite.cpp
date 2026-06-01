/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// SolidTransformationFinite — finite-strain (updated-Lagrangian) geometry
// method. See the header: every seam is identity; the element owns the F
// computation and the full K = ∫BᵀcB + ∫GᵀΣG assembly. The only finite-specific
// signal is getStrainMeasure() == DeformationGradient.

#include <SolidTransformationFinite.h>
#include <Matrix.h>
#include <Vector.h>

SolidTransformationFinite::SolidTransformationFinite()
  : SolidTransformation()
{
}

SolidTransformationFinite::~SolidTransformationFinite()
{
}

SolidTransformation *
SolidTransformationFinite::getCopy() const
{
  return new SolidTransformationFinite();
}

// The element computes F per Gauss point from the shape-function gradients (data
// the transformation does not hold), so there is no element-level state to
// refresh here.
int
SolidTransformationFinite::update(int /*numNodes*/, const Matrix & /*refCrds*/,
                                  const Matrix & /*curCrds*/)
{
  return 0;
}

// Identity: the deformation gradient is built directly from the global nodal
// displacements; no rigid rotation is stripped (that is corot, not finite).
int
SolidTransformationFinite::localizeDisp(const Vector &uGlobal, Vector &uCore) const
{
  if (&uCore != &uGlobal)
    uCore = uGlobal;
  return 0;
}

// Identity: ∫ Bᵀ σ dv is assembled on the current configuration in global
// coordinates — there is nothing to rotate.
int
SolidTransformationFinite::globalizeForce(const Vector &fCore, Vector &fGlobal) const
{
  if (&fGlobal != &fCore)
    fGlobal = fCore;
  return 0;
}

// Identity: the element assembles BOTH the material tangent ∫ Bᵀ c B dv and the
// geometric/initial-stress stiffness ∫ Gᵀ Σ G dv in its own Gauss loop (it owns
// the spatial gradients G and the Cauchy stress σ). See the LOCKED seam-3
// contract in FiniteStrainNDMaterial.h.
int
SolidTransformationFinite::globalizeStiff(const Matrix &kCore, const Vector & /*fCore*/,
                                          Matrix &kGlobal) const
{
  if (&kGlobal != &kCore)
    kGlobal = kCore;
  return 0;
}

// Cauchy stress is reported in the global frame.
void
SolidTransformationFinite::getCurrentFrame(Matrix &R) const
{
  R.Zero();
  R(0, 0) = 1.0;
  R(1, 1) = 1.0;
  R(2, 2) = 1.0;
}
