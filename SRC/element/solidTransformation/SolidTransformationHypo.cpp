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
// SolidTransformationHypo — hypoelastic rate-form updated-Lagrangian geometry
// method (ADR 79). See the header: every seam is identity; the element owns
// the per-GP objective-increment integration (LadrunoHypoKernel.h) and the
// full current-configuration assembly. The only hypo-specific signal is
// getMethodID() == METHOD_HYPO.

#include <SolidTransformationHypo.h>
#include <Matrix.h>
#include <Vector.h>

SolidTransformationHypo::SolidTransformationHypo()
  : SolidTransformation()
{
}

SolidTransformationHypo::~SolidTransformationHypo()
{
}

SolidTransformation *
SolidTransformationHypo::getCopy() const
{
  return new SolidTransformationHypo();
}

// The element computes the midpoint increment per Gauss point from its own
// shape-function gradients and committed/trial nodal displacements — there is
// no element-level state to refresh here.
int
SolidTransformationHypo::update(int /*numNodes*/, const Matrix & /*refCrds*/,
                                const Matrix & /*curCrds*/)
{
  return 0;
}

// Identity: the increment is built directly from the global nodal
// displacements; no whole-element rigid rotation is stripped (per-GP polar
// rotations handle objectivity inside the element).
int
SolidTransformationHypo::localizeDisp(const Vector &uGlobal, Vector &uCore) const
{
  if (&uCore != &uGlobal)
    uCore = uGlobal;
  return 0;
}

// Identity: ∫ Bᵀ σ dv is assembled on the current configuration in global
// coordinates — there is nothing to rotate.
int
SolidTransformationHypo::globalizeForce(const Vector &fCore, Vector &fGlobal) const
{
  if (&fGlobal != &fCore)
    fGlobal = fCore;
  return 0;
}

// Identity: the element assembles the pushed-forward material tangent
// ∫ Bᵀ c B dv AND the initial-stress geometric term ∫ Gᵀ Σ G dv in its own
// Gauss loop (it owns the spatial gradients and the pushed-forward Cauchy
// stress).
int
SolidTransformationHypo::globalizeStiff(const Matrix &kCore, const Vector & /*fCore*/,
                                        Matrix &kGlobal) const
{
  if (&kGlobal != &kCore)
    kGlobal = kCore;
  return 0;
}

// Stress responses are reported in the (fixed) unrotated material frame; the
// per-GP rotations do not define an element frame, so the recorder frame is
// the identity.
void
SolidTransformationHypo::getCurrentFrame(Matrix &R) const
{
  R.Zero();
  R(0, 0) = 1.0;
  R(1, 1) = 1.0;
  R(2, 2) = 1.0;
}
