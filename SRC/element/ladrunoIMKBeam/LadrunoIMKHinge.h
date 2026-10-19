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
#include <LadrunoCohesiveHinge.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <math.h>

inline bool ladrunoIMKMatIsImplex(UniaxialMaterial *m)
{
  if (m == 0) return false;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(m))
    return ch->isImplex();
  if (Bilin *bi = dynamic_cast<Bilin *>(m))
    return bi->isImplex();
  if (ModIMKPeakOriented *mk = dynamic_cast<ModIMKPeakOriented *>(m))
    return mk->isImplex();
  return false;
}

inline double ladrunoIMKMatImplexError(UniaxialMaterial *m)
{
  if (m == 0) return 0.0;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(m))
    return fabs(ch->getImplexError());
  if (Bilin *bi = dynamic_cast<Bilin *>(m))
    return fabs(bi->getImplexError());
  if (ModIMKPeakOriented *mk = dynamic_cast<ModIMKPeakOriented *>(m))
    return fabs(mk->getImplexError());
  return 0.0;
}

// Elastic core stiffness (a, b) and flexibility (F11, F12) coefficients for one
// bending axis, with optional Timoshenko shear rigidity GAv (Phi = 12*EI/(GAv*L^2)).
inline void ladrunoIMKCoreCoeffs(double L, double EI, double GAv,
                                 double &a, double &b, double &F11, double &F12)
{
  const double f = L / (6.0 * EI);
  if (GAv > 0.0 && L > 0.0) {
    const double Phi = 12.0 * EI / (GAv * L * L);
    const double inv = 1.0 / (1.0 + Phi);
    a = (EI / L) * (4.0 + Phi) * inv;
    b = (EI / L) * (2.0 - Phi) * inv;
    F11 = f * (2.0 + 0.5 * Phi);
    F12 = -f * (1.0 - 0.5 * Phi);
  } else {
    a = 4.0 * EI / L;
    b = 2.0 * EI / L;
    F11 = 2.0 * f;
    F12 = -f;
  }
}

// Parabolic axial-flexural (P-M) moment-capacity scaling factor:
//   f_PM(P) = clamp(1 - alpha * (P / P0) - (1 - alpha) * (P / P0)^2, 0.05, 2.5)
// where basic axial force P = q(0) is positive in tension (-P > 0 in compression).
inline double ladrunoIMKScalePM(double P, double P0, double alpha)
{
  if (P0 <= 0.0)
    return 1.0;
  const double n = P / P0;
  double f = 1.0 - alpha * n - (1.0 - alpha) * n * n;
  if (f < 0.05) f = 0.05;
  if (f > 2.50) f = 2.50;
  return f;
}

// Initial (small-strain) 2x2 tangent block for one bending axis, built from the
// hinge initial tangents (NULL material => elastic end), including optional
// Timoshenko shear GAv and rigid end offsets (dI, dJ).
inline void ladrunoIMKInitBlock(UniaxialMaterial *mi, UniaxialMaterial *mj,
                                double L, double EI, double k2[2][2],
                                double GAv = 0.0, double fPM = 1.0,
                                double dI = 0.0, double dJ = 0.0)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;

  double a, b, F11, F12;
  ladrunoIMKCoreCoeffs(Lc, EI, GAv, a, b, F11, F12);
  const double ktRef = a;
  const double ktFloor = 1.0e-8 * ktRef;
  double F21 = F12, F22 = F11;
  if (mi != 0) {
    double k0 = fPM * mi->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F11 += 1.0 / k0;
  }
  if (mj != 0) {
    double k0 = fPM * mj->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F22 += 1.0 / k0;
  }
  double det = F11 * F22 - F12 * F21;
  if (fabs(det) < 1.0e-30) det = (det < 0 ? -1.0e-30 : 1.0e-30);
  double kc[2][2];
  kc[0][0] = F22 / det;
  kc[0][1] = -F12 / det;
  kc[1][0] = -F21 / det;
  kc[1][1] = F11 / det;

  if (hasEZ) {
    const double T00 = (L - dJ) / Lc, T01 = dJ / Lc;
    const double T10 = dI / Lc,       T11 = (L - dI) / Lc;
    const double KT00 = kc[0][0] * T00 + kc[0][1] * T10;
    const double KT01 = kc[0][0] * T01 + kc[0][1] * T11;
    const double KT10 = kc[1][0] * T00 + kc[1][1] * T10;
    const double KT11 = kc[1][0] * T01 + kc[1][1] * T11;
    k2[0][0] = T00 * KT00 + T10 * KT10;
    k2[0][1] = T00 * KT01 + T10 * KT11;
    k2[1][0] = T01 * KT00 + T11 * KT10;
    k2[1][1] = T01 * KT01 + T11 * KT11;
  } else {
    k2[0][0] = kc[0][0];  k2[0][1] = kc[0][1];
    k2[1][0] = kc[1][0];  k2[1][1] = kc[1][1];
  }
}

// Per-axis internal state determination: solves the 2x2 internal Newton for the
// two hinge rotations of one bending axis, returns the basic forces (qi, qj at
// the two ends) and the condensed 2x2 tangent block k2. thi/thj are in/out: warm
// -started from the incoming trial hinge rotations and overwritten with the
// converged values. mi/mj are the end materials (NULL => that end is elastic).
inline void ladrunoIMKSolveAxis(UniaxialMaterial *mi, UniaxialMaterial *mj,
                                double vi, double vj, double L, double EI,
                                double &thi, double &thj,
                                double &qi, double &qj, double k2[2][2],
                                double q0i = 0.0, double q0j = 0.0,
                                double GAv = 0.0, double fPM = 1.0,
                                double dI = 0.0, double dJ = 0.0,
                                bool useImplex = false,
                                double thpCommitI = 0.0, double thpCommitJ = 0.0,
                                double dThpPrevI = 0.0, double dThpPrevJ = 0.0,
                                double *implexErrOut = 0,
                                double dtScale = 1.0,
                                double *thpOutI = 0, double *thpOutJ = 0)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;
  const double T00 = hasEZ ? (L - dJ) / Lc : 1.0;
  const double T01 = hasEZ ? dJ / Lc       : 0.0;
  const double T10 = hasEZ ? dI / Lc       : 0.0;
  const double T11 = hasEZ ? (L - dI) / Lc : 1.0;

  const double vic = T00 * vi + T01 * vj;
  const double vjc = T10 * vi + T11 * vj;
  const double q0ic = hasEZ ? (((L - dI) * q0i - dI * q0j) / L) : q0i;
  const double q0jc = hasEZ ? ((-dJ * q0i + (L - dJ) * q0j) / L) : q0j;

  double a, b, F11, F12;
  ladrunoIMKCoreCoeffs(Lc, EI, GAv, a, b, F11, F12);
  const double F11_core = F11;
  const double F12_core = F12;
  const double ktRef = a;          // stiffness scale for guards
  const double ktFloor = 1.0e-8 * ktRef;

  const bool hi = (mi != 0);
  const bool hj = (mj != 0);

  // warm-start from the incoming trial values; non-hinged ends stay at 0
  double th_i = hi ? thi : 0.0;
  double th_j = hj ? thj : 0.0;

  double Mi = 0.0, Mj = 0.0, kti = 0.0, ktj = 0.0;

  if (hi || hj) {
    const int maxIter = 50;
    double damp = 1.0;
    double gnormPrev = 1.0e300;
    for (int it = 0; it < maxIter; it++) {
      if (hi) {
        mi->setTrialStrain(th_i);
        Mi = fPM * mi->getStress();
        kti = fPM * mi->getTangent();
      }
      if (hj) {
        mj->setTrialStrain(th_j);
        Mj = fPM * mj->getStress();
        ktj = fPM * mj->getTangent();
      }

      // basic forces from elastic interior + fixed-end moments: q = K_el (v - thetaH) + q0
      const double qi_ = a * (vic - th_i) + b * (vjc - th_j) + q0ic;
      const double qj_ = b * (vic - th_i) + a * (vjc - th_j) + q0jc;

      // residual: hinged-end basic force must equal the hinge moment
      const double gi = hi ? (qi_ - Mi) : 0.0;
      const double gj = hj ? (qj_ - Mj) : 0.0;

      const double gnorm = sqrt(gi * gi + gj * gj);
      const double ref = fabs(qi_) + fabs(qj_) + fabs(Mi) + fabs(Mj) + fabs(q0ic) + fabs(q0jc);
      if (gnorm <= 1.0e-10 * ref + 1.0e-12)
        break;

      if (gnorm > gnormPrev) damp *= 0.5;
      else                   damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
      gnormPrev = gnorm;

      // guarded tangents for the Jacobian (force residual stays exact)
      double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
      double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;

      if (hi && hj) {
        // J = -(K_el_sub + diag(kt)) — symmetric 2x2; invert with eigenvalue floor
        const double J11 = -a - ki, J12 = -b;
        const double J22 = -a - kj;
        const double half = 0.5 * (J11 + J22);
        const double diff = 0.5 * (J11 - J22);
        const double rad  = sqrt(diff * diff + J12 * J12);
        const double lam1 = half + rad;
        const double lam2 = half - rad;
        double q1x, q1y;
        if (fabs(J12) > 1.0e-300) {
          q1x = lam1 - J22; q1y = J12;
        } else {
          if (J11 >= J22) { q1x = 1.0; q1y = 0.0; }
          else            { q1x = 0.0; q1y = 1.0; }
        }
        double n1 = sqrt(q1x * q1x + q1y * q1y);
        if (n1 < 1.0e-300) { q1x = 1.0; q1y = 0.0; n1 = 1.0; }
        q1x /= n1; q1y /= n1;
        const double q2x = -q1y, q2y = q1x;
        const double jFloor = 1.0e-8 * (fabs(a) + fabs(ki) + fabs(kj)) + 1.0e-300;
        const double f1 = (fabs(lam1) < jFloor) ? ((lam1 < 0.0) ? -jFloor : jFloor) : lam1;
        const double f2 = (fabs(lam2) < jFloor) ? ((lam2 < 0.0) ? -jFloor : jFloor) : lam2;
        const double i1 = 1.0 / f1, i2 = 1.0 / f2;
        const double i11 = i1 * q1x * q1x + i2 * q2x * q2x;
        const double i12 = i1 * q1x * q1y + i2 * q2x * q2y;
        const double i22 = i1 * q1y * q1y + i2 * q2y * q2y;
        const double dthi = -(i11 * gi + i12 * gj);
        const double dthj = -(i12 * gi + i22 * gj);
        th_i += damp * dthi;
        th_j += damp * dthj;
      } else if (hi) {
        double J11 = -a - ki;
        const double jFloor = 1.0e-8 * (fabs(a) + fabs(ki)) + 1.0e-300;
        if (fabs(J11) < jFloor) J11 = (J11 < 0.0) ? -jFloor : jFloor;
        th_i += damp * (-gi / J11);
      } else { // hj only
        double J22 = -a - kj;
        const double jFloor = 1.0e-8 * (fabs(a) + fabs(kj)) + 1.0e-300;
        if (fabs(J22) < jFloor) J22 = (J22 < 0.0) ? -jFloor : jFloor;
        th_j += damp * (-gj / J22);
      }
    }
  }

  thi = th_i;
  thj = th_j;

  // exact basic forces at the converged hinge rotations (+ fixed-end moments)
  double qic = a * (vic - th_i) + b * (vjc - th_j) + q0ic;
  double qjc = b * (vic - th_i) + a * (vjc - th_j) + q0jc;

  // Implicit plastic hinge rotations: thp = th - q_c / k0
  double Fh0_i = 0.0, Fh0_j = 0.0;
  if (hi) {
    double k0i = fPM * mi->getInitialTangent();
    if (fabs(k0i) < ktFloor) k0i = (k0i < 0 ? -ktFloor : ktFloor);
    Fh0_i = 1.0 / k0i;
  }
  if (hj) {
    double k0j = fPM * mj->getInitialTangent();
    if (fabs(k0j) < ktFloor) k0j = (k0j < 0 ? -ktFloor : ktFloor);
    Fh0_j = 1.0 / k0j;
  }
  const double thp_i = hi ? (th_i - Fh0_i * qic) : 0.0;
  const double thp_j = hj ? (th_j - Fh0_j * qjc) : 0.0;
  if (thpOutI != 0) *thpOutI = thp_i;
  if (thpOutJ != 0) *thpOutJ = thp_j;

  const bool matImplex = ladrunoIMKMatIsImplex(mi) || ladrunoIMKMatIsImplex(mj);

  if (useImplex && !matImplex) {
    const double thpExI = hi ? (thpCommitI + dtScale * dThpPrevI) : 0.0;
    const double thpExJ = hj ? (thpCommitJ + dtScale * dThpPrevJ) : 0.0;
    const double F11_0 = F11_core + Fh0_i;
    const double F22_0 = F11_core + Fh0_j;
    double det0 = F11_0 * F22_0 - F12_core * F12_core;
    if (fabs(det0) < 1.0e-30) det0 = (det0 < 0 ? -1.0e-30 : 1.0e-30);
    const double kc0_00 =  F22_0 / det0;
    const double kc0_01 = -F12_core / det0;
    const double kc0_10 = -F12_core / det0;
    const double kc0_11 =  F11_0 / det0;

    const double rhs_i = vic - thpExI + F11_core * q0ic + F12_core * q0jc;
    const double rhs_j = vjc - thpExJ + F12_core * q0ic + F11_core * q0jc;
    const double qic_ex = kc0_00 * rhs_i + kc0_01 * rhs_j;
    const double qjc_ex = kc0_10 * rhs_i + kc0_11 * rhs_j;
    if (implexErrOut != 0)
      *implexErrOut = fabs(qic - qic_ex) + fabs(qjc - qjc_ex);
    qic = qic_ex;
    qjc = qjc_ex;
    qi = T00 * qic + T10 * qjc;
    qj = T01 * qic + T11 * qjc;
    ladrunoIMKInitBlock(mi, mj, L, EI, k2, GAv, fPM, dI, dJ);
    return;
  }

  if (implexErrOut != 0) {
    *implexErrOut = matImplex
      ? (ladrunoIMKMatImplexError(mi) + ladrunoIMKMatImplexError(mj))
      : 0.0;
  }

  qi = T00 * qic + T10 * qjc;
  qj = T01 * qic + T11 * qjc;

  // condensed tangent block:  K = (F_el + F_h)^{-1}
  double F21 = F12, F22 = F11;
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
  double kc[2][2];
  kc[0][0] = F22 / det;
  kc[0][1] = -F12 / det;
  kc[1][0] = -F21 / det;
  kc[1][1] = F11 / det;

  if (hasEZ) {
    const double KT00 = kc[0][0] * T00 + kc[0][1] * T10;
    const double KT01 = kc[0][0] * T01 + kc[0][1] * T11;
    const double KT10 = kc[1][0] * T00 + kc[1][1] * T10;
    const double KT11 = kc[1][0] * T01 + kc[1][1] * T11;
    k2[0][0] = T00 * KT00 + T10 * KT10;
    k2[0][1] = T00 * KT01 + T10 * KT11;
    k2[1][0] = T01 * KT00 + T11 * KT10;
    k2[1][1] = T01 * KT01 + T11 * KT11;
  } else {
    k2[0][0] = kc[0][0];  k2[0][1] = kc[0][1];
    k2[1][0] = kc[1][0];  k2[1][1] = kc[1][1];
  }
}

#endif
