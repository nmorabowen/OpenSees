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

// Ladruno (WP-124) -- shared Element-contract helpers for the continuum element
// shells: LadrunoQuad, LadrunoCST, LadrunoLST, LadrunoCSTPair, LadrunoBrick,
// LadrunoBrick20, BezierTri6, BezierTet10.
//
// WHY. The eight elements re-implement the same Element-contract shell, and
// eight PRs had to fix one defect in several copies (#228 Ki cache, #562
// Rayleigh P-clobber, #670 recvSelf cache invalidation, #852 ground sign, ...).
// A fix made HERE reaches every element. Plan and inventory:
// Ladruno_implementation/124_continuum_shell_helpers.md.
//
// FREE FUNCTIONS, NOT A BASE CLASS (owner decision): the elements differ in
// node count, DOF layout and caches, and the Rayleigh tail needs Element's
// PROTECTED members, so the snapshot + getRayleighDampingForces() step stays in
// each element. Every helper is a byte-for-byte refactor of the code it
// replaced (proven by Ladruno_implementation/wp123_undamped/fingerprint.py
// --suite shells) and is pinned by tests/test_ladruno_element_shell_helpers.py.
//
// NEW ELEMENT? Use these instead of copying another element's shell
// (.claude/skills/ladruno-new-element/SKILL.md).

#ifndef LadrunoElementShell_h
#define LadrunoElementShell_h

#include <Matrix.h>
#include <Vector.h>
#include <Node.h>
#include <OPS_Globals.h>

namespace LadrunoShell {

// ---- initial-stiffness cache (the #228 contract) ---------------------------
// getInitialStiff() forms K0 into scratch (often a CLASS-STATIC Matrix that the
// next getTangentStiff of ANY instance overwrites), then
//     return LadrunoShell::cacheKi(Ki, formed);
// The cache owns a byte-copy; return *Ki, never the scratch. Ki is formed once
// (the vanilla convention); anything that invalidates the reference state
// (recvSelf into a live element) calls dropKi(Ki).
inline const Matrix &cacheKi(Matrix *&Ki, const Matrix &formed)
{
  if (Ki != 0)
    delete Ki;
  Ki = new Matrix(formed);
  return *Ki;
}

inline void dropKi(Matrix *&Ki)
{
  if (Ki != 0) {
    delete Ki;
    Ki = 0;
  }
}


// ---- ground-motion inertia load (addInertiaLoadToUnbalance) ----------------
// The OpenSees convention: the load vector Q accumulates -M R a_g and the
// residual SUBTRACTS Q, so the unbalance gains +M R a_g. `+M a` shakes the mesh
// the wrong way, silently (BezierTri6/Tet10 until WP-117 #852).
//   diagOnly = true : Q(i) += -M(i,i) ra(i)      (lumped/diagonal mass; plane four)
//   diagOnly = false: Q    += -M ra              (full M; Bezier, bricks)
// The massless early-out stays in the element (its predicate differs per
// element: material/element rho vs a probe of M's diagonal). ra is the
// caller's scratch (size nen*ndf), filled with the nodal R a_g. checkSize=false
// keeps the bricks' historical no-check behaviour (C13, owner decision).
inline int addGroundInertia(Vector &Q, const Matrix &M, Node **nodes, int nen,
                            int ndf, bool diagOnly, const Vector &accel,
                            Vector &ra, const char *who, bool checkSize = true)
{
  for (int a = 0; a < nen; a++) {
    const Vector &Raccel = nodes[a]->getRV(accel);
    if (checkSize && Raccel.Size() != ndf) {
      opserr << who << "::addInertiaLoadToUnbalance - matrix and vector sizes incompatible\n";
      return -1;
    }
    for (int j = 0; j < ndf; j++)
      ra(a * ndf + j) = Raccel(j);
  }
  if (diagOnly) {
    const int n = nen * ndf;
    for (int i = 0; i < n; i++)
      Q(i) += -M(i, i) * ra(i);
  } else
    Q.addMatrixVector(1.0, M, ra, -1.0);
  return 0;
}

} // namespace LadrunoShell

#endif
