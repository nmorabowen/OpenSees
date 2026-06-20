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

// Ladruno: shared mass-conserving (HRZ) diagonal-lumping utility. See
// Ladruno_implementation/35_ladruno_hrz_lumped_mass_adr.md.
//
// Given a CONSISTENT element mass matrix Mc, return the Hinton-Rock-Zienkiewicz
// diagonal lump: per translational spatial direction d, scale the consistent
// diagonal entries so their sum equals the element's total mass in that
// direction,
//        alpha_d        = m_d / S_d,   S_d = sum_{a in d} Mc_aa,
//                                       m_d = sum_{a in d} rowsum_a,
//        M^HRZ_aa,d     = alpha_d * Mc_aa.
// m_d is the sum of FULL row sums over direction-d DOFs; it equals the element's
// total mass in direction d when the consistent mass has no cross-coupling into
// direction d (true for global-frame translational blocks and for symmetric
// element mass where the uy<->rz coupling is antisymmetric and cancels in the
// sum). It is the correct conserved target in every shipped/tested case here.
// This is simultaneously (a) positive (Mc_aa > 0 for a properly integrated
// consistent mass), (b) mass-conserving (sum_a M^HRZ_aa,d = m_d), and (c)
// rotation-aware: a rotational DOF (no single direction) is scaled by the MEAN
// of the element's translational alphas, mean_d(alpha_d) -- well-defined even
// when alpha_x != alpha_y (e.g. a beam: axial 140 vs transverse 156).
//
// Reduces to row-sum EXACTLY on geometrically regular elements (square quad,
// parallelepiped hex, equal-length bar) where all Mc_aa in a direction are
// equal; diverges per-node on distorted/higher-order elements (where row-sum is
// wrong: zero/negative or zero-rotational) but STILL conserves total mass.
//
// The core math (hrzLumpRaw) is OpenSees-free so it can be g++-verified against
// a numpy oracle standalone (the LadrunoJ2Kernel.h / LadrunoHardening.h
// precedent): compile any consumer / the standalone driver with
// -DLADRUNO_MASSLUMPING_STANDALONE to drop the OpenSees Matrix/ID wrapper.
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoMassLumping_h
#define LadrunoMassLumping_h

#ifndef LADRUNO_MASSLUMPING_STANDALONE
#include <Matrix.h>
#include <ID.h>
#endif

namespace Ladruno {

// Return codes from the HRZ lumper.
enum HRZStatus {
  HRZ_OK            = 0,  // HRZ lump applied
  HRZ_FELLBACK_DIAG = 1,  // a guard failed -> mdiag holds the diagonal-of-consistent
                          //   (caller should warn; Diagonal is non-conserving)
  HRZ_BAD_INPUT     = -1  // n <= 0 / null buffers
};

// Core math, OpenSees-free. `Mc` is the n x n consistent mass packed ROW-MAJOR
// (Mc[i*n + j]); `dir[i]` is the spatial direction 0..ndm-1 for a translational
// DOF, or -1 for a rotational DOF. Fills `mdiag[i]` with the HRZ diagonal.
//
// Guards (per ADR, fall back to diagonal-of-consistent + warn if any fails):
//   - every translational Mc_aa > 0,
//   - every per-direction S_d > 0 and m_d > 0,
//   - every produced M^HRZ_aa > 0.
inline int hrzLumpRaw(const double *Mc, const int *dir, int n, double *mdiag)
{
  if (n <= 0 || Mc == 0 || dir == 0 || mdiag == 0)
    return HRZ_BAD_INPUT;

  // infer ndm = max translational direction index + 1 (capped at 3)
  int ndm = 0;
  for (int i = 0; i < n; ++i) {
    int d = dir[i];
    if (d >= 0 && d + 1 > ndm) ndm = d + 1;
  }
  if (ndm > 3) ndm = 3;

  // accumulate per-direction S_d (sum of consistent diagonal) and
  // m_d (sum of full row sums = total element mass in that direction)
  double S[3]  = {0.0, 0.0, 0.0};
  double Mt[3] = {0.0, 0.0, 0.0};
  bool ok = true;

  for (int i = 0; i < n && ok; ++i) {
    int d = dir[i];
    if (d < 0 || d >= ndm) continue;          // rotational handled below
    double mii = Mc[i * n + i];
    if (!(mii > 0.0)) { ok = false; break; }  // guard: positive consistent diagonal
    S[d] += mii;
    double rs = 0.0;
    for (int j = 0; j < n; ++j) rs += Mc[i * n + j];
    Mt[d] += rs;
  }

  if (ok)
    for (int d = 0; d < ndm; ++d)
      if (!(S[d] > 0.0) || !(Mt[d] > 0.0)) { ok = false; break; }

  if (ok) {
    double alpha[3] = {0.0, 0.0, 0.0};
    double alphaSum = 0.0;
    for (int d = 0; d < ndm; ++d) { alpha[d] = Mt[d] / S[d]; alphaSum += alpha[d]; }
    double alphaMean = (ndm > 0) ? alphaSum / ndm : 0.0;  // for rotational DOFs

    for (int i = 0; i < n; ++i) {
      int d = dir[i];
      double mii = Mc[i * n + i];
      if (d >= 0 && d < ndm) {
        // translational DOF: must stay strictly positive, scale by its
        // direction's HRZ factor.
        mdiag[i] = alpha[d] * mii;
        if (!(mdiag[i] > 0.0)) { ok = false; break; }
      } else {
        // rotational / non-translational DOF (shell drilling with ZERO
        // rotational mass; u-p pressure with a NEGATIVE compressibility term):
        // scale a genuinely positive inertia by the mean alpha (the HRZ
        // rotational benefit), but pass a legitimately non-positive entry
        // THROUGH UNCHANGED rather than failing the whole element. The
        // downstream eigensolve already tolerates a non-positive M diagonal
        // (DSYGVX -> DGGEV fallback), exactly as the Diagonal lump would here.
        mdiag[i] = (mii > 0.0) ? alphaMean * mii : mii;
      }
    }
  }

  if (!ok) {
    // fallback: diagonal-of-consistent (non-conserving last resort; caller warns)
    for (int i = 0; i < n; ++i) mdiag[i] = Mc[i * n + i];
    return HRZ_FELLBACK_DIAG;
  }
  return HRZ_OK;
}

#ifndef LADRUNO_MASSLUMPING_STANDALONE
// OpenSees wrapper: copies the consistent Matrix M into a row-major buffer and
// calls the raw kernel. `dofDir(i)` = direction 0..ndm-1 (translational) or -1
// (rotational); build it from the element's node layout
// (ele->getNodePtrs()[k]->getNumberDOF(): first ndm DOFs translational).
inline int hrzLump(const Matrix &M, const ID &dofDir, double *mdiag, int n)
{
  if (n <= 0 || M.noRows() != n || M.noCols() != n || dofDir.Size() != n)
    return HRZ_BAD_INPUT;
  double *Mc  = new double[n * n];
  int    *dir = new int[n];
  for (int i = 0; i < n; ++i) {
    dir[i] = dofDir(i);
    for (int j = 0; j < n; ++j) Mc[i * n + j] = M(i, j);
  }
  int rc = hrzLumpRaw(Mc, dir, n, mdiag);
  delete[] Mc;
  delete[] dir;
  return rc;
}
#endif

} // namespace Ladruno

#endif
