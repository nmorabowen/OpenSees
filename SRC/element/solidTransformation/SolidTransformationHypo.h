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
// Description: SolidTransformationHypo — the hypoelastic RATE-FORM
// updated-Lagrangian geometry method (v4, ADR 79). Large strain WITH the
// existing small-strain material library: the element integrates the
// Hughes–Winget midpoint objective strain increment in the UNROTATED
// (Green–Naghdi) material frame (LadrunoHypoKernel.h), accumulates it, and
// drives an UNCHANGED small-strain NDMaterial via setTrialStrain — then
// pushes the material stress/tangent forward with R = polar(F) and assembles
// on the current configuration with the initial-stress geometric term.
//
// Like SolidTransformationFinite this class is a pure MARKER: the hypo path
// needs PER-GAUSS-POINT machinery (midpoint gradients, per-GP polar rotation,
// per-GP accumulated strain) that the seam interface — whole-element vectors
// only — cannot see, so all of it lives in the consuming element. Every seam
// method here is identity; the only signals are getMethodID() == METHOD_HYPO
// and the name. getStrainMeasure() is honestly SmallStrain: the material IS
// driven through setTrialStrain with a small-strain-convention Voigt vector
// (the accumulated objective increments).
//
// The material and its whole tensorial state live in the FIXED unrotated
// frame (the corot convention — kinematic hardening stays frame-consistent
// with no state rotation ever; ADR 79 §2). Stress responses are reported in
// that material frame; frameTimeVarying() stays false because the rotation is
// per-GP, not per-element (a global-frame push-forward recorder channel is a
// demand-driven follow-up).

#ifndef SolidTransformationHypo_h
#define SolidTransformationHypo_h

#include <SolidTransformation.h>

class SolidTransformationHypo : public SolidTransformation
{
 public:
  SolidTransformationHypo();
  virtual ~SolidTransformationHypo();

  int getMethodID() const { return METHOD_HYPO; }
  const char *getName() const { return "hypo"; }
  SolidTransformation *getCopy() const;

  StrainMeasure getStrainMeasure() const { return StrainMeasure::SmallStrain; }

  int update(int numNodes, const Matrix &refCrds, const Matrix &curCrds);

  int localizeDisp(const Vector &uGlobal, Vector &uCore) const;

  int globalizeForce(const Vector &fCore, Vector &fGlobal) const;
  int globalizeStiff(const Matrix &kCore, const Vector &fCore,
                     Matrix &kGlobal) const;

  bool frameTimeVarying() const { return false; }
  void getCurrentFrame(Matrix &R) const;
};

#endif
