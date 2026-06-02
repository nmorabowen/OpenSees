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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// LadrunoIMKHinge: dimension-agnostic per-bending-axis hinge kernel shared by
// the 2D and 3D LadrunoIMKBeam macro elements. Each function operates on ONE
// bending axis: an elastic interior (4EI/L) in series with up to two end
// rotational hinges (one independent uniaxial law per end). A NULL material
// pointer means that end is elastic (no hinge). The 3D element calls these
// twice (strong axis Mz, weak axis My); the 2D element calls them once (Mz).
//
// These are header-only free functions (no Element coupling) so both elements
// share one kernel and cannot drift -- the same extraction pattern used by
// LadrunoJ2Kernel.h / LadrunoHardening.h. See LadrunoIMKBeam.cpp and
// Ladruno_implementation/14_ladruno_imk_beam.md for the formulation.

#ifndef LadrunoIMKHinge_h
#define LadrunoIMKHinge_h

#include <UniaxialMaterial.h>
#include <math.h>

// Per-axis internal state determination: solves the 2x2 internal Newton for the
// two hinge rotations of one bending axis, returns the basic forces (qi, qj at
// the two ends) and the condensed 2x2 tangent block k2. thi/thj are in/out: warm
// -started from the incoming trial hinge rotations and overwritten with the
// converged values. mi/mj are the end materials (NULL => that end is elastic).
inline void ladrunoIMKSolveAxis(UniaxialMaterial *mi, UniaxialMaterial *mj,
                                double vi, double vj, double L, double EI,
                                double &thi, double &thj,
                                double &qi, double &qj, double k2[2][2])
{
  const double a = 4.0 * EI / L;   // elastic stiffness block  K_el = [[a,b],[b,a]]
  const double b = 2.0 * EI / L;
  const double ktRef = a;          // stiffness scale for guards
  const double ktFloor = 1.0e-8 * ktRef;

  const bool hi = (mi != 0);
  const bool hj = (mj != 0);

  // warm-start from the incoming trial values; non-hinged ends stay at 0
  double th_i = hi ? thi : 0.0;
  double th_j = hj ? thj : 0.0;

  double Mi = 0.0, Mj = 0.0, kti = 0.0, ktj = 0.0;

  if (hi || hj) {
    const int maxIter = 25;
    for (int it = 0; it < maxIter; it++) {
      if (hi) {
        mi->setTrialStrain(th_i);
        Mi = mi->getStress();
        kti = mi->getTangent();
      }
      if (hj) {
        mj->setTrialStrain(th_j);
        Mj = mj->getStress();
        ktj = mj->getTangent();
      }

      // basic forces from elastic interior: q = K_el (v - thetaH)
      const double qi_ = a * (vi - th_i) + b * (vj - th_j);
      const double qj_ = b * (vi - th_i) + a * (vj - th_j);

      // residual: hinged-end basic force must equal the hinge moment
      const double gi = hi ? (qi_ - Mi) : 0.0;
      const double gj = hj ? (qj_ - Mj) : 0.0;

      const double gnorm = sqrt(gi * gi + gj * gj);
      const double ref = fabs(qi_) + fabs(qj_) + fabs(Mi) + fabs(Mj);
      if (gnorm <= 1.0e-10 * ref + 1.0e-12)
        break;

      // guarded tangents for the Jacobian (force residual stays exact)
      double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
      double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;

      if (hi && hj) {
        // J = -(K_el_sub + diag(kt))
        const double J11 = -a - ki, J12 = -b;
        const double J21 = -b,       J22 = -a - kj;
        double det = J11 * J22 - J12 * J21;
        if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
        // dth = -J^{-1} g
        const double dthi = -(J22 * gi - J12 * gj) / det;
        const double dthj = -(-J21 * gi + J11 * gj) / det;
        th_i += dthi;
        th_j += dthj;
      } else if (hi) {
        const double J11 = -a - ki;
        th_i += -gi / J11;
      } else { // hj only
        const double J22 = -a - kj;
        th_j += -gj / J22;
      }
    }
  }

  thi = th_i;
  thj = th_j;

  // exact basic forces at the converged hinge rotations
  qi = a * (vi - th_i) + b * (vj - th_j);
  qj = b * (vi - th_i) + a * (vj - th_j);

  // condensed tangent block:  K = (F_el + F_h)^{-1}
  const double f = L / (6.0 * EI);     // F_el = f * [[2,-1],[-1,2]]
  double F11 = 2.0 * f, F12 = -f;
  double F21 = -f,       F22 = 2.0 * f;
  if (hi) {
    double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
    F11 += 1.0 / ki;
  }
  if (hj) {
    double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;
    F22 += 1.0 / kj;
  }
  double det = F11 * F22 - F12 * F21;
  if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
  k2[0][0] = F22 / det;
  k2[0][1] = -F12 / det;
  k2[1][0] = -F21 / det;
  k2[1][1] = F11 / det;
}

// Initial (small-strain) 2x2 tangent block for one bending axis, built from the
// hinge initial tangents (NULL material => elastic end).
inline void ladrunoIMKInitBlock(UniaxialMaterial *mi, UniaxialMaterial *mj,
                                double L, double EI, double k2[2][2])
{
  const double a = 4.0 * EI / L;
  const double ktRef = a;
  const double ktFloor = 1.0e-8 * ktRef;
  const double f = L / (6.0 * EI);
  double F11 = 2.0 * f, F12 = -f;
  double F21 = -f,       F22 = 2.0 * f;
  if (mi != 0) {
    double k0 = mi->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F11 += 1.0 / k0;
  }
  if (mj != 0) {
    double k0 = mj->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F22 += 1.0 / k0;
  }
  double det = F11 * F22 - F12 * F21;
  if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
  k2[0][0] = F22 / det;
  k2[0][1] = -F12 / det;
  k2[1][0] = -F21 / det;
  k2[1][1] = F11 / det;
}

#endif
