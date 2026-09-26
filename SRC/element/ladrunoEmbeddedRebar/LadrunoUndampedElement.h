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

// Ladruno (WP-123): base for elements that carry NO Rayleigh damping.
//
// Users: the four pure-penalty couplings LadrunoDistributingCoupling (RBE3),
// LadrunoKinematicCoupling (RBE2), LadrunoEmbeddedNode and LadrunoEmbeddedRebar. Each
// used to carry its own copy of this code, and the implicit-transient crash of a
// missing getDamp override was fixed in them one after another (#219, then #220).
// A new coupling / tie / constraint element that must ignore Rayleigh damping
// derives from this class instead of writing the overrides again.
//
// * setRayleighDampingFactors refuses the factors: it returns 0 WITHOUT storing
//   them, so Element's alphaM/betaK/betaK0/betaKc stay 0 and a betaK can never
//   shrink the reported explicit dt_cr (CriticalTimeStep reads the element's
//   factors; ADR 28 §5, ADR 20 §10.6).
//
// * getDamp returns an element-owned zero matrix sized getNumDOF(). Since the
//   2026-07-28 Element.cpp fix the base getDamp would also answer zero here
//   (LEDGER_quirks "makes 11 `Element` methods dereference `theMatrices[-1]`"); the
//   override keeps these elements independent of that vanilla edit surviving an
//   upstream sync (LEDGER_quirks "A no-op `setRayleighDampingFactors` WITHOUT a
//   `getDamp` override").
//
// * getRayleighDampingForces is deliberately NOT here: it is not virtual in Element,
//   so a same-named method would only SHADOW it -- the four pre-WP-123 copies were
//   never called. Element's own version already answers zero (stored factors are 0),
//   and that is what eleResponse 'dampingForce' reaches.

#ifndef LadrunoUndampedElement_h
#define LadrunoUndampedElement_h

#include <Element.h>
#include <Matrix.h>

class LadrunoUndampedElement : public Element
{
 public:
  LadrunoUndampedElement(int tag, int classTag) : Element(tag, classTag) {}
  virtual ~LadrunoUndampedElement() {}

  int setRayleighDampingFactors(double, double, double, double) { return 0; }

  const Matrix &getDamp(void)
  {
    const int n = this->getNumDOF();
    if (zeroDamp.noRows() != n || zeroDamp.noCols() != n)
      zeroDamp.resize(n, n);
    zeroDamp.Zero();
    return zeroDamp;
  }

 private:
  Matrix zeroDamp;          // per-instance, never shared: D == 0
};

#endif
