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

// LadrunoElasticStrainProvider: a mixin (WP-144, plan 2.8, owner decision 2026-10-01, G2 option c).
//
// An inner NDMaterial that carries its OWN elastic strain tensor implements this, so that a
// finite-strain wrapper (LogStrainNDMaterial) can read the updated elastic strain straight from the
// inner instead of recovering it as inv(D0) : tau. That recovery is exact only for a LINEAR-elastic
// inner (constant D0); for a nonlinear hyperelastic inner (LadrunoNorSand: pressure-dependent K, mu)
// it is wrong, and the committed b^e drifts (the wrapper is then not objective).
//
// The wrapper does  dynamic_cast<LadrunoElasticStrainProvider*>(inner)  and, if non-null and the call
// returns true, builds b^e = exp[2 eps^e] from the PROVIDED strain; otherwise it keeps its original
// D0-inversion path unchanged (every other inner material stays byte-identical).
//
// Contract of ladrunoGetElasticStrain:
//   - the inner's CURRENT TRIAL elastic strain (the one belonging to the stress getStress() returns),
//   - engineering Voigt order (11,22,33,12,23,13): shear entries are gamma = 2 eps_tensor,
//   - expressed in the inner's own (log-strain) frame: no rotation or push-forward is applied,
//   - epsE must be sized 6 by the caller; return false (epsE untouched) if unavailable.

#ifndef LadrunoElasticStrainProvider_h
#define LadrunoElasticStrainProvider_h

#include <Vector.h>

class LadrunoElasticStrainProvider {
 public:
  virtual ~LadrunoElasticStrainProvider() {}
  // The inner material's CURRENT TRIAL elastic strain, engineering Voigt (11,22,33,12,23,13) in its own
  // (log-strain) frame. Return false if unavailable.
  virtual bool ladrunoGetElasticStrain(Vector &epsE) const = 0;
};

#endif
