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

// Ladruno: shared artificial diagonal-mass machinery for the fork's two
// fictitious-mass static/quasi-static integrators (LadrunoDynamicRelaxation and
// LadrunoArcLength). Both build an integrator-OWNED diagonal "M*" by a Gershgorin
// row-sum of the element tangent stiffness (or lumped real mass, or unity) and
// poke it onto the diagonal SOE — deliberately NOT via Element::getMass(), which
// is zero on the zero-density models these solvers target (the ADR-20 BLOCKER).
//
// Keeping the row-sum builder and the diagonal poke in ONE header is the same
// "extraction" precedent as LadrunoJ2Kernel.h / LadrunoHardening.h: the two
// integrators carried near-identical private copies; sharing them keeps the M*
// definition byte-identical across both. The integrator-specific scaling stays
// at the call site (DR's dt²/4 Gershgorin prefactor and viscous rescale; the
// ArcLength cOverDt diagonal weight), passed in as plain doubles so the shared
// code is a pure, bit-faithful lift.
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.
// See Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md and
//     Ladruno_implementation/21_ladruno_dynamic_relaxation_adr.md.

#ifndef LadrunoFictitiousMass_h
#define LadrunoFictitiousMass_h

#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Element.h>
#include <math.h>

namespace Ladruno {

// Build the integrator-owned diagonal M* by an SOE-free element probe:
//   mode 0 (gershgorin): m_i = dt2Quarter * sum_j |K_ij|   (row-sum of element
//          tangent stiffness). dt2Quarter <= 0 means "no prefactor" (factor 1),
//          which is how LadrunoArcLength keeps the raw row-sum (its scale lives
//          in cOverDt); LadrunoDynamicRelaxation passes 0.25*dt*dt (already
//          divided by s^2 for the viscous case) so omega*dt ~ 2 by construction.
//   mode 1 (lumped):     m_i = massScale * sum_j |M_ij|     (row-sum of the real
//          element mass; needs a nonzero-density model).
//   mode 2 (unity):      m_i = 1.
// The mmax*1e-8 floor lifts any zero / non-positive diagonal (free DOF with no
// stiffness/mass) to keep M* invertible. Mstar must already be sized to nEqn.
// Returns 0 on success, -1 if no AnalysisModel.
inline int
buildGershgorinDiagonal(AnalysisModel *theModel, Vector &Mstar, int mode,
                        double massScale, double dt2Quarter)
{
  if (theModel == 0) return -1;
  int size = Mstar.Size();
  Mstar.Zero();

  if (mode == 2) {                       // unity
    for (int i = 0; i < size; i++) Mstar(i) = 1.0;
    return 0;
  }

  // gershgorin / lumped: accumulate |row-sums| of the element K (or M) into the
  // global diagonal at the element's equation numbers.
  FE_EleIter &theEles = theModel->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEles()) != 0) {
    Element *e = elePtr->getElement();
    if (e == 0) continue;                // subdomain / no backing element
    const ID &id = elePtr->getID();
    const Matrix &Ke = (mode == 1) ? e->getMass() : e->getTangentStiff();
    int n = id.Size();
    if (Ke.noRows() != n) continue;      // defensive
    for (int i = 0; i < n; i++) {
      int loc = id(i);
      if (loc < 0) continue;
      double rs = 0.0;
      for (int j = 0; j < n; j++) rs += fabs(Ke(i, j));
      Mstar(loc) += rs;
    }
  }

  if (mode == 0) {                       // gershgorin stability scaling
    double f = (dt2Quarter > 0.0) ? dt2Quarter : 1.0;
    if (f != 1.0)
      for (int i = 0; i < size; i++) Mstar(i) *= f;
  } else {                               // lumped real mass * scale
    for (int i = 0; i < size; i++) Mstar(i) *= massScale;
  }

  // floor any zero / non-positive diagonal (free DOF with no stiffness/mass)
  double mmax = 0.0;
  for (int i = 0; i < size; i++) if (Mstar(i) > mmax) mmax = Mstar(i);
  double floor = (mmax > 0.0) ? mmax * 1.0e-8 : 1.0;
  for (int i = 0; i < size; i++)
    if (!(Mstar(i) > 0.0)) Mstar(i) = floor;

  return 0;
}

// Poke the diagonal M* onto the SOE: A_ii <- (zeroFirst ? 0 : A_ii) + fac*M*_i.
//   zeroFirst=true  (LadrunoDynamicRelaxation): replace A with fac*M* so the
//                   solve degenerates to a M*^{-1} apply (matrix-free LHS).
//   zeroFirst=false (LadrunoArcLength): add fac*M* to the real K already
//                   assembled by the base formTangent (diagonal regularization).
inline void
addDiagonalToSOE(LinearSOE *theSOE, const Vector &Mstar, int size, double fac,
                 bool zeroFirst)
{
  if (theSOE == 0) return;
  if (zeroFirst) theSOE->zeroA();
  static Matrix m1(1, 1);
  static ID id1(1);
  for (int i = 0; i < size; i++) {
    m1(0, 0) = fac * Mstar(i);
    id1(0) = i;
    theSOE->addA(m1, id1);
  }
}

} // namespace Ladruno

#endif
