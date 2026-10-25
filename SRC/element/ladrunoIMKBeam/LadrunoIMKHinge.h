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
#include <NDMaterial.h>
#include <LadrunoCohesiveHinge.h>
#include <LadrunoCohesiveHingeBiaxial.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <Vector.h>
#include <Matrix.h>
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

inline void ladrunoIMKMatSetImplex(UniaxialMaterial *m, bool flag)
{
  if (m == 0) return;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(m))
    ch->setImplex(flag);
  else if (Bilin *bi = dynamic_cast<Bilin *>(m))
    bi->setImplex(flag);
  else if (ModIMKPeakOriented *mk = dynamic_cast<ModIMKPeakOriented *>(m))
    mk->setImplex(flag);
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

inline bool ladrunoIMKBiaxialIsImplex(NDMaterial *m)
{
  if (m == 0) return false;
  if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(m))
    return cb->isImplex();
  return false;
}

inline void ladrunoIMKBiaxialSetImplex(NDMaterial *m, bool flag)
{
  if (m == 0) return;
  if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(m))
    cb->setImplex(flag);
}

inline double ladrunoIMKBiaxialImplexError(NDMaterial *m)
{
  if (m == 0) return 0.0;
  if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(m))
    return fabs(cb->getImplexError());
  return 0.0;
}

// Spectral 2x2 symmetric inversion with eigenvalue floor (handles softening).
inline void ladrunoIMKInvert2x2SymFloored(double A00, double A01, double A11,
                                          double floorVal,
                                          double &I00, double &I01, double &I11)
{
  const double half = 0.5 * (A00 + A11);
  const double diff = 0.5 * (A00 - A11);
  const double rad  = sqrt(diff * diff + A01 * A01);
  const double lam1 = half + rad;
  const double lam2 = half - rad;
  double q1x, q1y;
  if (fabs(A01) > 1.0e-300) {
    q1x = lam1 - A11; q1y = A01;
  } else {
    if (A00 >= A11) { q1x = 1.0; q1y = 0.0; }
    else            { q1x = 0.0; q1y = 1.0; }
  }
  double n1 = sqrt(q1x * q1x + q1y * q1y);
  if (n1 < 1.0e-300) { q1x = 1.0; q1y = 0.0; n1 = 1.0; }
  q1x /= n1; q1y /= n1;
  const double q2x = -q1y, q2y = q1x;
  const double f1 = (fabs(lam1) < floorVal) ? ((lam1 < 0.0) ? -floorVal : floorVal) : lam1;
  const double f2 = (fabs(lam2) < floorVal) ? ((lam2 < 0.0) ? -floorVal : floorVal) : lam2;
  const double i1 = 1.0 / f1, i2 = 1.0 / f2;
  I00 = i1 * q1x * q1x + i2 * q2x * q2x;
  I01 = i1 * q1x * q1y + i2 * q2x * q2y;
  I11 = i1 * q1y * q1y + i2 * q2y * q2y;
}

// General 2x2 inversion with determinant floor (for non-symmetric mode-mix tangents).
inline void ladrunoIMKInvert2x2Floored(double K00, double K01, double K10, double K11,
                                       double floorVal,
                                       double &F00, double &F01, double &F10, double &F11)
{
  if (fabs(K01 - K10) <= 1.0e-12 * (fabs(K00) + fabs(K11) + fabs(K01) + 1.0e-300)) {
    double I00, I01, I11;
    ladrunoIMKInvert2x2SymFloored(K00, 0.5 * (K01 + K10), K11, floorVal, I00, I01, I11);
    F00 = I00; F01 = I01; F10 = I01; F11 = I11;
    return;
  }
  double det = K00 * K11 - K01 * K10;
  const double detFloor = floorVal * floorVal;
  if (fabs(det) < detFloor)
    det = (det < 0.0) ? -detFloor : detFloor;
  const double invDet = 1.0 / det;
  F00 =  K11 * invDet;
  F01 = -K01 * invDet;
  F10 = -K10 * invDet;
  F11 =  K00 * invDet;
}

// 4x4 Gauss-Jordan inversion with pivot floor.
inline void ladrunoIMKInvert4x4(const double A[4][4], double Ainv[4][4], double pivotFloor)
{
  double M[4][8];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = A[i][j];
      M[i][j + 4] = (i == j) ? 1.0 : 0.0;
    }
  }
  for (int col = 0; col < 4; col++) {
    int piv = col;
    double maxVal = fabs(M[col][col]);
    for (int row = col + 1; row < 4; row++) {
      if (fabs(M[row][col]) > maxVal) {
        maxVal = fabs(M[row][col]);
        piv = row;
      }
    }
    if (piv != col) {
      for (int j = 0; j < 8; j++) {
        double tmp = M[col][j];
        M[col][j] = M[piv][j];
        M[piv][j] = tmp;
      }
    }
    double p = M[col][col];
    if (fabs(p) < pivotFloor)
      p = (p < 0.0) ? -pivotFloor : pivotFloor;
    const double invP = 1.0 / p;
    for (int j = 0; j < 8; j++)
      M[col][j] *= invP;
    for (int row = 0; row < 4; row++) {
      if (row == col) continue;
      const double fac = M[row][col];
      if (fac != 0.0) {
        for (int j = 0; j < 8; j++)
          M[row][j] -= fac * M[col][j];
      }
    }
  }
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      Ainv[i][j] = M[i][j + 4];
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

inline double ladrunoIMKScalePMDeriv(double P, double P0, double alpha)
{
  if (P0 <= 0.0)
    return 0.0;
  const double n = P / P0;
  const double f = 1.0 - alpha * n - (1.0 - alpha) * n * n;
  if (f <= 0.05 || f >= 2.50)
    return 0.0;
  return (-alpha - 2.0 * (1.0 - alpha) * n) / P0;
}

// 4-parameter RC balance-point P-M scaling factor (DR-3):
inline double ladrunoIMKScalePMBal(double P, double P0c, double P0t, double Pbal, double muBal)
{
  if (P0c <= 0.0 || P0t <= 0.0 || Pbal <= 0.0 || Pbal >= P0c)
    return 1.0;
  double f = 1.0;
  if (P >= 0.0) {
    f = 1.0 - P / P0t;
  } else if (P >= -Pbal) {
    const double u = (P + Pbal) / Pbal;
    f = 1.0 + (muBal - 1.0) * (1.0 - u * u);
  } else {
    const double denom = P0c - Pbal;
    const double u = (P + Pbal) / denom;
    f = muBal * (1.0 - u * u);
  }
  if (f < 0.05) f = 0.05;
  if (f > 3.0)  f = 3.0;
  return f;
}

inline double ladrunoIMKScalePMBalDeriv(double P, double P0c, double P0t, double Pbal, double muBal)
{
  if (P0c <= 0.0 || P0t <= 0.0 || Pbal <= 0.0 || Pbal >= P0c)
    return 0.0;
  double f = 1.0;
  double df = 0.0;
  if (P >= 0.0) {
    f = 1.0 - P / P0t;
    df = -1.0 / P0t;
  } else if (P >= -Pbal) {
    const double u = (P + Pbal) / Pbal;
    f = 1.0 + (muBal - 1.0) * (1.0 - u * u);
    df = -2.0 * (muBal - 1.0) * (P + Pbal) / (Pbal * Pbal);
  } else {
    const double denom = P0c - Pbal;
    const double u = (P + Pbal) / denom;
    f = muBal * (1.0 - u * u);
    df = -2.0 * muBal * (P + Pbal) / (denom * denom);
  }
  if (f <= 0.05 || f >= 3.0)
    return 0.0;
  return df;
}

inline double ladrunoIMKMatDLogPM(UniaxialMaterial *m, double dLogPMFallback)
{
  if (m == 0) return 0.0;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge *>(m)) {
    if (ch->hasPM()) {
      const double f = ch->getPMMultiplier();
      return (f > 0.0) ? (ch->getPMDerivative() / f) : 0.0;
    }
  }
  return dLogPMFallback;
}

inline double ladrunoIMKBiaxialDLogPM(NDMaterial *m, double dLogPMFallback)
{
  if (m == 0) return 0.0;
  if (LadrunoCohesiveHingeBiaxial *cb = dynamic_cast<LadrunoCohesiveHingeBiaxial *>(m)) {
    if (cb->hasPM()) {
      const double f = cb->getPMMultiplier();
      return (f > 0.0) ? (cb->getPMDerivative() / f) : 0.0;
    }
  }
  return dLogPMFallback;
}

// Initial (small-strain) 2x2 tangent block for one bending axis, built from the
// hinge initial tangents (NULL material => elastic end), including optional
// Timoshenko shear GAv, rigid end offsets (dI, dJ), series shear hinge mShear,
// and optional Powell active-set compliance subtraction (usePowell).
inline void ladrunoIMKInitBlock(UniaxialMaterial *mi, UniaxialMaterial *mj,
                                double L, double EI, double k2[2][2],
                                double GAv = 0.0, double fPM = 1.0,
                                double dI = 0.0, double dJ = 0.0,
                                UniaxialMaterial *mShear = 0,
                                bool usePowell = false)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;

  double a, b, F11, F12;
  ladrunoIMKCoreCoeffs(Lc, EI, GAv, a, b, F11, F12);
  const double ktRef = a;
  const double ktFloor = 1.0e-8 * ktRef;
  double F21 = F12, F22 = F11;
  const double fPMi = (dynamic_cast<LadrunoCohesiveHinge *>(mi) == 0) ? fPM : 1.0;
  const double fPMj = (dynamic_cast<LadrunoCohesiveHinge *>(mj) == 0) ? fPM : 1.0;
  if (!usePowell && mi != 0) {
    double k0 = fPMi * mi->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F11 += 1.0 / k0;
  }
  if (!usePowell && mj != 0) {
    double k0 = fPMj * mj->getInitialTangent();
    if (fabs(k0) < ktFloor) k0 = (k0 < 0 ? -ktFloor : ktFloor);
    F22 += 1.0 / k0;
  }
  if (mShear != 0 && Lc > 0.0) {
    double ks0 = mShear->getInitialTangent() * Lc * Lc;
    if (fabs(ks0) < ktFloor) ks0 = (ks0 < 0 ? -ktFloor : ktFloor);
    const double Fs0 = 1.0 / ks0;
    F11 += Fs0;
    F12 += Fs0;
    F21 += Fs0;
    F22 += Fs0;
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

// Per-axis internal state determination: solves the internal Newton for the
// two hinge rotations (+ optional series shear hinge deformation gammaS) of one
// bending axis, returns the basic forces (qi, qj at the two ends) and the
// condensed 2x2 tangent block k2.
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
                                double *thpOutI = 0, double *thpOutJ = 0,
                                UniaxialMaterial *mShear = 0,
                                double *gammaSInOut = 0,
                                double gammaSpCommit = 0.0,
                                double dGammaSpPrev = 0.0,
                                double *gammaSpOut = 0,
                                double dLogPM = 0.0,
                                double *dq_dN = 0,
                                bool usePowell = false)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;
  const double T00 = hasEZ ? (L - dJ) / Lc : 1.0;
  const double T01 = hasEZ ? dJ / Lc       : 0.0;
  const double T10 = hasEZ ? dI / Lc       : 0.0;
  const double T11 = hasEZ ? (L - dI) / Lc : 1.0;

  const double vic = T00 * vi + T01 * vj;
  const double vjc = T10 * vi + T11 * vj;
  // Clear-span fixed-end moments (already integrated over Lc in addLoad)
  const double q0ic = q0i;
  const double q0jc = q0j;

  double a, b, F11, F12;
  ladrunoIMKCoreCoeffs(Lc, EI, GAv, a, b, F11, F12);
  const double F11_core = F11;
  const double F12_core = F12;
  const double ktRef = a;          // stiffness scale for guards
  const double ktFloor = 1.0e-8 * ktRef;

  const bool hi = (mi != 0);
  const bool hj = (mj != 0);
  const bool hs = (mShear != 0 && Lc > 0.0);
  const double fPMi = (dynamic_cast<LadrunoCohesiveHinge *>(mi) == 0) ? fPM : 1.0;
  const double fPMj = (dynamic_cast<LadrunoCohesiveHinge *>(mj) == 0) ? fPM : 1.0;

  double Fh0_i = 0.0, Fh0_j = 0.0, Fs0 = 0.0, k0i_raw = 0.0, k0j_raw = 0.0, ks0_raw = 0.0;
  if (hi) {
    k0i_raw = fPMi * mi->getInitialTangent();
    if (fabs(k0i_raw) < ktFloor) k0i_raw = (k0i_raw < 0 ? -ktFloor : ktFloor);
    Fh0_i = 1.0 / k0i_raw;
  }
  if (hj) {
    k0j_raw = fPMj * mj->getInitialTangent();
    if (fabs(k0j_raw) < ktFloor) k0j_raw = (k0j_raw < 0 ? -ktFloor : ktFloor);
    Fh0_j = 1.0 / k0j_raw;
  }
  if (hs) {
    ks0_raw = mShear->getInitialTangent();
    double ks0_rot = ks0_raw * Lc * Lc;
    if (fabs(ks0_rot) < ktFloor) ks0_rot = (ks0_rot < 0 ? -ktFloor : ktFloor);
    Fs0 = 1.0 / ks0_rot;
  }
  // Elastic shear rotation under fixed-end moments at v = 0:
  // subtracting gam_rot0 ensures q == q0 at v == 0 in the elastic range under asymmetric span loads.
  const double gam_rot0 = hs ? (-Fs0 * (q0ic + q0jc)) : 0.0;
  const double powellScale = usePowell ? 1.0 : 0.0;

  const bool matImplexI = ladrunoIMKMatIsImplex(mi);
  const bool matImplexJ = ladrunoIMKMatIsImplex(mj);
  const bool matImplexS = ladrunoIMKMatIsImplex(mShear);
  const bool matImplex = matImplexI || matImplexJ || matImplexS;

  // warm-start from the incoming trial values; non-hinged ends stay at 0
  double th_i = hi ? thi : 0.0;
  double th_j = hj ? thj : 0.0;
  double gam_s = (hs && gammaSInOut != 0) ? *gammaSInOut : 0.0;

  double Mi = 0.0, Mj = 0.0, Vs = 0.0, kti = 0.0, ktj = 0.0, kts = 0.0;

  if (hi || hj || hs) {
    // FB-10: when element-level -implex (useImplex) is active, disable any
    // redundant material-level -implex on the element's private hinge clones so
    // the inner series equilibrium is solved implicitly and extrapolated once
    // at the element level (also prevents Bilin::commitState from overwriting
    // Tangent = k0).
    if (useImplex) {
      if (matImplexI) ladrunoIMKMatSetImplex(mi, false);
      if (matImplexJ) ladrunoIMKMatSetImplex(mj, false);
      if (matImplexS) ladrunoIMKMatSetImplex(mShear, false);
    }

    const int maxIter = 50;
    double damp = 1.0;
    double gnormPrev = 1.0e300;
    for (int it = 0; it < maxIter; it++) {
      if (hi) {
        mi->setTrialStrain(th_i);
        Mi = fPMi * mi->getStress();
        kti = fPMi * mi->getTangent();
      }
      if (hj) {
        mj->setTrialStrain(th_j);
        Mj = fPMj * mj->getStress();
        ktj = fPMj * mj->getTangent();
      }
      if (hs) {
        mShear->setTrialStrain(gam_s);
        Vs = mShear->getStress();
        kts = mShear->getTangent();
      }

      const double thCore_i = hi ? (th_i - powellScale * Fh0_i * Mi) : 0.0;
      const double thCore_j = hj ? (th_j - powellScale * Fh0_j * Mj) : 0.0;
      const double gam_rot = hs ? (gam_s / Lc - gam_rot0) : 0.0;
      const double vel_i = vic - thCore_i + gam_rot;
      const double vel_j = vjc - thCore_j + gam_rot;

      // basic forces from elastic interior + fixed-end moments: q = K_el (v - thetaH + (gam_s - gam_s0)/Lc) + q0
      const double qi_ = a * vel_i + b * vel_j + q0ic;
      const double qj_ = b * vel_i + a * vel_j + q0jc;

      // residual: hinged-end basic force must equal the hinge moment;
      // shear residual scaled by Lc so all residuals have moment units:
      //   gs_mom = Lc * (-(qi_ + qj_) / Lc - Vs) = -(qi_ + qj_) - Lc * Vs
      const double gi = hi ? (qi_ - Mi) : 0.0;
      const double gj = hj ? (qj_ - Mj) : 0.0;
      const double gs_mom = hs ? (-(qi_ + qj_) - Lc * Vs) : 0.0;

      const double gnorm = sqrt(gi * gi + gj * gj + gs_mom * gs_mom);
      const double ref = fabs(qi_) + fabs(qj_) + fabs(Mi) + fabs(Mj)
                       + (hs ? Lc * fabs(Vs) : 0.0) + fabs(q0ic) + fabs(q0jc);
      if (gnorm <= 1.0e-10 * ref + 1.0e-12)
        break;

      if (gnorm > gnormPrev) damp *= 0.5;
      else                   damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
      gnormPrev = gnorm;

      // guarded tangents for the Jacobian (force residual stays exact)
      double ki = (fabs(kti) < ktFloor) ? (kti < 0 ? -ktFloor : ktFloor) : kti;
      double kj = (fabs(ktj) < ktFloor) ? (ktj < 0 ? -ktFloor : ktFloor) : ktj;
      const double eta_i = hi ? (1.0 - powellScale * Fh0_i * kti) : 0.0;
      const double eta_j = hj ? (1.0 - powellScale * Fh0_j * ktj) : 0.0;

      double a_eff = a;
      double b_eff = b;
      double gi_eff = gi;
      double gj_eff = gj;
      double Ds = 1.0;
      const double ab = a + b;
      if (hs) {
        const double ks_rot = kts * Lc * Lc;
        const double ks = (fabs(ks_rot) < ktFloor) ? (ks_rot < 0 ? -ktFloor : ktFloor) : ks_rot;
        Ds = 2.0 * ab + ks;
        const double dsFloor = 1.0e-8 * (2.0 * fabs(ab) + fabs(ks)) + 1.0e-300;
        if (fabs(Ds) < dsFloor) Ds = (Ds < 0.0) ? -dsFloor : dsFloor;
        a_eff = a - (ab * ab) / Ds;
        b_eff = b - (ab * ab) / Ds;
        if (hi) gi_eff = gi + (ab / Ds) * gs_mom;
        if (hj) gj_eff = gj + (ab / Ds) * gs_mom;
      }

      double dthi = 0.0, dthj = 0.0;
      if (hi && hj) {
        const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(ki) + fabs(kj)) + 1.0e-300;
        if (!usePowell) {
          // J = -(K_el_sub + diag(kt)) — symmetric 2x2; invert with eigenvalue floor
          const double J11 = -a_eff - ki, J12 = -b_eff;
          const double J22 = -a_eff - kj;
          double i11, i12, i22;
          ladrunoIMKInvert2x2SymFloored(J11, J12, J22, jFloor, i11, i12, i22);
          dthi = -(i11 * gi_eff + i12 * gj_eff);
          dthj = -(i12 * gi_eff + i22 * gj_eff);
        } else {
          double i00, i01, i10, i11;
          ladrunoIMKInvert2x2Floored(a_eff * eta_i + ki, b_eff * eta_j,
                                     b_eff * eta_i,      a_eff * eta_j + kj,
                                     jFloor, i00, i01, i10, i11);
          dthi = i00 * gi_eff + i01 * gj_eff;
          dthj = i10 * gi_eff + i11 * gj_eff;
        }
      } else if (hi) {
        double J11 = -a_eff * eta_i - ki;
        const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(ki)) + 1.0e-300;
        if (fabs(J11) < jFloor) J11 = (J11 < 0.0) ? -jFloor : jFloor;
        dthi = -gi_eff / J11;
      } else if (hj) { // hj only
        double J22 = -a_eff * eta_j - kj;
        const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(kj)) + 1.0e-300;
        if (fabs(J22) < jFloor) J22 = (J22 < 0.0) ? -jFloor : jFloor;
        dthj = -gj_eff / J22;
      }

      th_i += damp * dthi;
      th_j += damp * dthj;
      if (hs) {
        const double dgam_rot = (gs_mom + ab * (eta_i * dthi + eta_j * dthj)) / Ds;
        gam_s += damp * (dgam_rot * Lc);
      }
    }
  }

  thi = th_i;
  thj = th_j;
  if (gammaSInOut != 0) *gammaSInOut = gam_s;

  // exact basic forces at the converged hinge rotations (+ fixed-end moments)
  const double thCore_i = hi ? (th_i - powellScale * Fh0_i * Mi) : 0.0;
  const double thCore_j = hj ? (th_j - powellScale * Fh0_j * Mj) : 0.0;
  const double gam_rot = hs ? (gam_s / Lc - gam_rot0) : 0.0;
  double qic = a * (vic - thCore_i + gam_rot) + b * (vjc - thCore_j + gam_rot) + q0ic;
  double qjc = b * (vic - thCore_i + gam_rot) + a * (vjc - thCore_j + gam_rot) + q0jc;

  // Implicit plastic hinge rotations: thp = th - M / k0 (exact 0.0 in elastic regime)
  const double thp_i = hi ? (th_i - Mi / k0i_raw) : 0.0;
  const double thp_j = hj ? (th_j - Mj / k0j_raw) : 0.0;
  const double gam_sp = (hs && fabs(ks0_raw) > 1.0e-30) ? (gam_s - Vs / ks0_raw) : 0.0;
  if (thpOutI != 0) *thpOutI = thp_i;
  if (thpOutJ != 0) *thpOutJ = thp_j;
  if (gammaSpOut != 0) *gammaSpOut = gam_sp;

  if (useImplex) {
    const double thpExI = hi ? (thpCommitI + dtScale * dThpPrevI) : 0.0;
    const double thpExJ = hj ? (thpCommitJ + dtScale * dThpPrevJ) : 0.0;
    const double gamSpEx = hs ? (gammaSpCommit + dtScale * dGammaSpPrev) : 0.0;
    const double F11_0 = F11_core + (usePowell ? 0.0 : Fh0_i) + Fs0;
    const double F22_0 = F11_core + (usePowell ? 0.0 : Fh0_j) + Fs0;
    const double F12_0 = F12_core + Fs0;
    double det0 = F11_0 * F22_0 - F12_0 * F12_0;
    if (fabs(det0) < 1.0e-30) det0 = (det0 < 0 ? -1.0e-30 : 1.0e-30);
    const double kc0_00 =  F22_0 / det0;
    const double kc0_01 = -F12_0 / det0;
    const double kc0_10 = -F12_0 / det0;
    const double kc0_11 =  F11_0 / det0;

    const double rhs_i = vic - thpExI + (hs ? (gamSpEx / Lc - gam_rot0) : 0.0) + F11_core * q0ic + F12_core * q0jc;
    const double rhs_j = vjc - thpExJ + (hs ? (gamSpEx / Lc - gam_rot0) : 0.0) + F12_core * q0ic + F11_core * q0jc;
    const double qic_ex = kc0_00 * rhs_i + kc0_01 * rhs_j;
    const double qjc_ex = kc0_10 * rhs_i + kc0_11 * rhs_j;
    const double dthp_i = (thp_i - thpExI) - (hs ? (gam_sp - gamSpEx) / Lc : 0.0);
    const double dthp_j = (thp_j - thpExJ) - (hs ? (gam_sp - gamSpEx) / Lc : 0.0);
    const double dqic_imp = kc0_00 * dthp_i + kc0_01 * dthp_j;
    const double dqjc_imp = kc0_10 * dthp_i + kc0_11 * dthp_j;
    if (implexErrOut != 0)
      *implexErrOut = fabs(dqic_imp) + fabs(dqjc_imp);
    qic = qic_ex;
    qjc = qjc_ex;
    qi = T00 * qic + T10 * qjc;
    qj = T01 * qic + T11 * qjc;
    if (dq_dN != 0) {
      if (!usePowell) {
        const double dLogI = hi ? ladrunoIMKMatDLogPM(mi, dLogPM) : 0.0;
        const double dLogJ = hj ? ladrunoIMKMatDLogPM(mj, dLogPM) : 0.0;
        const double drhs_i = dLogI * Fh0_i * qic_ex;
        const double drhs_j = dLogJ * Fh0_j * qjc_ex;
        const double dqic_dN = kc0_00 * drhs_i + kc0_01 * drhs_j;
        const double dqjc_dN = kc0_10 * drhs_i + kc0_11 * drhs_j;
        dq_dN[0] = T00 * dqic_dN + T10 * dqjc_dN;
        dq_dN[1] = T01 * dqic_dN + T11 * dqjc_dN;
      } else {
        dq_dN[0] = 0.0;
        dq_dN[1] = 0.0;
      }
    }
    ladrunoIMKInitBlock(mi, mj, L, EI, k2, GAv, fPM, dI, dJ, mShear, usePowell);
    return;
  }

  if (implexErrOut != 0) {
    *implexErrOut = matImplex
      ? (ladrunoIMKMatImplexError(mi) + ladrunoIMKMatImplexError(mj) + ladrunoIMKMatImplexError(mShear))
      : 0.0;
  }

  qi = T00 * qic + T10 * qjc;
  qj = T01 * qic + T11 * qjc;

  // Stiffness Schur-complement condensed tangent block:
  //   First condense series shear hinge into K_eff, then condense rotational hinges
  //   via K_c = K_eff - K_eff * E * (K_eff * E + K_H)^{-1} * K_eff (never divides by kti/ktj/kts).
  double a_eff = a, b_eff = b;
  if (hs) {
    const double ab = a + b;
    const double ks_rot = kts * Lc * Lc;
    double Ds = 2.0 * ab + ks_rot;
    const double dsFloor = 1.0e-8 * (2.0 * fabs(ab) + fabs(ks_rot)) + 1.0e-300;
    if (fabs(Ds) < dsFloor) Ds = (Ds < 0.0) ? -dsFloor : dsFloor;
    a_eff = a - (ab * ab) / Ds;
    b_eff = b - (ab * ab) / Ds;
  }
  double kc[2][2] = {
    { a_eff, b_eff },
    { b_eff, a_eff }
  };
  const double dMi_dN = hi ? (ladrunoIMKMatDLogPM(mi, dLogPM) * Mi) : 0.0;
  const double dMj_dN = hj ? (ladrunoIMKMatDLogPM(mj, dLogPM) * Mj) : 0.0;
  const double eta_i = hi ? (1.0 - powellScale * Fh0_i * kti) : 0.0;
  const double eta_j = hj ? (1.0 - powellScale * Fh0_j * ktj) : 0.0;
  double dqic_dN = 0.0, dqjc_dN = 0.0;
  if (hi && hj) {
    const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(kti) + fabs(ktj)) + 1.0e-300;
    if (!usePowell) {
      const double A11 = a_eff + kti;
      const double A12 = b_eff;
      const double A22 = a_eff + ktj;
      double i11, i12, i22;
      ladrunoIMKInvert2x2SymFloored(A11, A12, A22, jFloor, i11, i12, i22);
      const double M00 = a_eff * i11 + b_eff * i12;
      const double M01 = a_eff * i12 + b_eff * i22;
      const double M10 = b_eff * i11 + a_eff * i12;
      const double M11 = b_eff * i12 + a_eff * i22;
      kc[0][0] -= M00 * a_eff + M01 * b_eff;
      kc[0][1] -= M00 * b_eff + M01 * a_eff;
      kc[1][0] -= M10 * a_eff + M11 * b_eff;
      kc[1][1] -= M10 * b_eff + M11 * a_eff;
      dqic_dN = M00 * dMi_dN + M01 * dMj_dN;
      dqjc_dN = M10 * dMi_dN + M11 * dMj_dN;
    } else {
      double i00, i01, i10, i11;
      ladrunoIMKInvert2x2Floored(a_eff * eta_i + kti, b_eff * eta_j,
                                 b_eff * eta_i,       a_eff * eta_j + ktj,
                                 jFloor, i00, i01, i10, i11);
      const double w00 = eta_i * i00;
      const double w01 = 0.5 * (eta_i * i01 + eta_j * i10);
      const double w11 = eta_j * i11;
      const double M00 = a_eff * w00 + b_eff * w01;
      const double M01 = a_eff * w01 + b_eff * w11;
      const double M10 = b_eff * w00 + a_eff * w01;
      const double M11 = b_eff * w01 + a_eff * w11;
      kc[0][0] -= M00 * a_eff + M01 * b_eff;
      kc[0][1] -= M00 * b_eff + M01 * a_eff;
      kc[1][0] -= M10 * a_eff + M11 * b_eff;
      kc[1][1] -= M10 * b_eff + M11 * a_eff;
      dqic_dN = M00 * dMi_dN + M01 * dMj_dN;
      dqjc_dN = M10 * dMi_dN + M11 * dMj_dN;
    }
  } else if (hi) {
    double A11 = a_eff * eta_i + kti;
    const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(kti)) + 1.0e-300;
    if (fabs(A11) < jFloor) A11 = (A11 < 0.0) ? -jFloor : jFloor;
    const double w00 = eta_i / A11;
    kc[0][0] -= (a_eff * a_eff) * w00;
    kc[0][1] -= (a_eff * b_eff) * w00;
    kc[1][0] -= (b_eff * a_eff) * w00;
    kc[1][1] -= (b_eff * b_eff) * w00;
    dqic_dN = a_eff * w00 * dMi_dN;
    dqjc_dN = b_eff * w00 * dMi_dN;
  } else if (hj) {
    double A22 = a_eff * eta_j + ktj;
    const double jFloor = 1.0e-8 * (fabs(a_eff) + fabs(ktj)) + 1.0e-300;
    if (fabs(A22) < jFloor) A22 = (A22 < 0.0) ? -jFloor : jFloor;
    const double w11 = eta_j / A22;
    kc[0][0] -= (b_eff * b_eff) * w11;
    kc[0][1] -= (b_eff * a_eff) * w11;
    kc[1][0] -= (a_eff * b_eff) * w11;
    kc[1][1] -= (a_eff * a_eff) * w11;
    dqic_dN = b_eff * w11 * dMj_dN;
    dqjc_dN = a_eff * w11 * dMj_dN;
  }
  if (dq_dN != 0) {
    dq_dN[0] = T00 * dqic_dN + T10 * dqjc_dN;
    dq_dN[1] = T01 * dqic_dN + T11 * dqjc_dN;
  }

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

// Initial 4x4 bending tangent block for 3D coupled biaxial NDMaterial hinges
// at ends I and J (indices [0,1,2,3] = [Mz_i, Mz_j, My_i, My_j]).
inline void ladrunoIMKInitBiaxial(NDMaterial *biI, NDMaterial *biJ,
                                  double L, double EIz, double EIy,
                                  double k4[4][4],
                                  double GAvy = 0.0, double GAvz = 0.0,
                                  double fPM = 1.0,
                                  double dI = 0.0, double dJ = 0.0,
                                  UniaxialMaterial *mShearY = 0,
                                  UniaxialMaterial *mShearZ = 0,
                                  bool usePowell = false)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;

  double az, bz, F11z, F12z;
  double ay, by, F11y, F12y;
  ladrunoIMKCoreCoeffs(Lc, EIz, GAvy, az, bz, F11z, F12z);
  ladrunoIMKCoreCoeffs(Lc, EIy, GAvz, ay, by, F11y, F12y);
  const double ktRef = 0.5 * (az + ay);
  const double ktFloor = 1.0e-8 * ktRef;

  double Fc[4][4] = {
    { F11z, F12z, 0.0,  0.0  },
    { F12z, F11z, 0.0,  0.0  },
    { 0.0,  0.0,  F11y, F12y },
    { 0.0,  0.0,  F12y, F11y }
  };

  if (mShearY != 0 && Lc > 0.0) {
    double ks0 = mShearY->getInitialTangent() * Lc * Lc;
    if (fabs(ks0) < ktFloor) ks0 = (ks0 < 0 ? -ktFloor : ktFloor);
    const double Fs0 = 1.0 / ks0;
    Fc[0][0] += Fs0; Fc[0][1] += Fs0;
    Fc[1][0] += Fs0; Fc[1][1] += Fs0;
  }
  if (mShearZ != 0 && Lc > 0.0) {
    double ks0 = mShearZ->getInitialTangent() * Lc * Lc;
    if (fabs(ks0) < ktFloor) ks0 = (ks0 < 0 ? -ktFloor : ktFloor);
    const double Fs0 = 1.0 / ks0;
    Fc[2][2] += Fs0; Fc[2][3] += Fs0;
    Fc[3][2] += Fs0; Fc[3][3] += Fs0;
  }

  const double fPMI = (dynamic_cast<LadrunoCohesiveHingeBiaxial *>(biI) == 0) ? fPM : 1.0;
  const double fPMJ = (dynamic_cast<LadrunoCohesiveHingeBiaxial *>(biJ) == 0) ? fPM : 1.0;

  if (!usePowell && biI != 0) {
    const Matrix &K0 = biI->getInitialTangent();
    double F00, F01, F10, F11;
    ladrunoIMKInvert2x2Floored(fPMI * K0(0, 0), fPMI * K0(0, 1),
                               fPMI * K0(1, 0), fPMI * K0(1, 1),
                               ktFloor, F00, F01, F10, F11);
    Fc[0][0] += F00; Fc[0][2] += F01;
    Fc[2][0] += F10; Fc[2][2] += F11;
  }
  if (!usePowell && biJ != 0) {
    const Matrix &K0 = biJ->getInitialTangent();
    double F00, F01, F10, F11;
    ladrunoIMKInvert2x2Floored(fPMJ * K0(0, 0), fPMJ * K0(0, 1),
                               fPMJ * K0(1, 0), fPMJ * K0(1, 1),
                               ktFloor, F00, F01, F10, F11);
    Fc[1][1] += F00; Fc[1][3] += F01;
    Fc[3][1] += F10; Fc[3][3] += F11;
  }

  double Kc[4][4];
  ladrunoIMKInvert4x4(Fc, Kc, 1.0e-30);

  if (hasEZ) {
    const double T00 = (L - dJ) / Lc, T01 = dJ / Lc;
    const double T10 = dI / Lc,       T11 = (L - dI) / Lc;
    const double T[4][4] = {
      { T00, T01, 0.0, 0.0 },
      { T10, T11, 0.0, 0.0 },
      { 0.0, 0.0, T00, T01 },
      { 0.0, 0.0, T10, T11 }
    };
    double KT[4][4];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        KT[i][j] = 0.0;
        for (int m = 0; m < 4; m++)
          KT[i][j] += Kc[i][m] * T[m][j];
      }
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        k4[i][j] = 0.0;
        for (int m = 0; m < 4; m++)
          k4[i][j] += T[m][i] * KT[m][j];
      }
  } else {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        k4[i][j] = Kc[i][j];
  }
}

// Coupled 4x4 internal state determination for 3D biaxial NDMaterial hinges
// at ends I and J (with optional series shear hinges mShearY, mShearZ).
inline void ladrunoIMKSolveBiaxial(NDMaterial *biI, NDMaterial *biJ,
                                   const double v4[4], double L,
                                   double EIz, double EIy,
                                   double th4[4], double q4[4], double k4[4][4],
                                   const double q0_4[4],
                                   double GAvy = 0.0, double GAvz = 0.0,
                                   double fPM = 1.0,
                                   double dI = 0.0, double dJ = 0.0,
                                   bool useImplex = false,
                                   const double thpCommit[4] = 0,
                                   const double dThpPrev[4] = 0,
                                   double *implexErrOut = 0,
                                   double dtScale = 1.0,
                                   double thpOut[4] = 0,
                                   UniaxialMaterial *mShearY = 0,
                                   UniaxialMaterial *mShearZ = 0,
                                   double *gammaSYInOut = 0,
                                   double *gammaSZInOut = 0,
                                   double gammaSpYCommit = 0.0,
                                   double gammaSpZCommit = 0.0,
                                   double dGammaSpYPrev = 0.0,
                                   double dGammaSpZPrev = 0.0,
                                   double *gammaSpYOut = 0,
                                   double *gammaSpZOut = 0,
                                   double dLogPM = 0.0,
                                   double *dq4_dN = 0,
                                   bool usePowell = false)
{
  const bool hasEZ = (dI > 0.0 || dJ > 0.0) && (dI + dJ < L);
  const double Lc = hasEZ ? (L - dI - dJ) : L;
  const double T00 = hasEZ ? (L - dJ) / Lc : 1.0;
  const double T01 = hasEZ ? dJ / Lc       : 0.0;
  const double T10 = hasEZ ? dI / Lc       : 0.0;
  const double T11 = hasEZ ? (L - dI) / Lc : 1.0;

  const double vc[4] = {
    T00 * v4[0] + T01 * v4[1],
    T10 * v4[0] + T11 * v4[1],
    T00 * v4[2] + T01 * v4[3],
    T10 * v4[2] + T11 * v4[3]
  };
  // Clear-span fixed-end moments (already integrated over Lc in addLoad)
  const double q0c[4] = { q0_4[0], q0_4[1], q0_4[2], q0_4[3] };

  double az, bz, F11z, F12z;
  double ay, by, F11y, F12y;
  ladrunoIMKCoreCoeffs(Lc, EIz, GAvy, az, bz, F11z, F12z);
  ladrunoIMKCoreCoeffs(Lc, EIy, GAvz, ay, by, F11y, F12y);
  const double ktRef = 0.5 * (az + ay);
  const double ktFloor = 1.0e-8 * ktRef;

  const bool hI = (biI != 0);
  const bool hJ = (biJ != 0);
  const bool hsY = (mShearY != 0 && Lc > 0.0);
  const bool hsZ = (mShearZ != 0 && Lc > 0.0);
  const double fPMI = (dynamic_cast<LadrunoCohesiveHingeBiaxial *>(biI) == 0) ? fPM : 1.0;
  const double fPMJ = (dynamic_cast<LadrunoCohesiveHingeBiaxial *>(biJ) == 0) ? fPM : 1.0;

  // Initial flexibility for plastic hinge rotation, Powell active-set, and IMPL-EX
  double Fh0I[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double Fh0J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  if (hI) {
    const Matrix &K0 = biI->getInitialTangent();
    ladrunoIMKInvert2x2Floored(fPMI * K0(0, 0), fPMI * K0(0, 1),
                               fPMI * K0(1, 0), fPMI * K0(1, 1),
                               ktFloor, Fh0I[0][0], Fh0I[0][1], Fh0I[1][0], Fh0I[1][1]);
  }
  if (hJ) {
    const Matrix &K0 = biJ->getInitialTangent();
    ladrunoIMKInvert2x2Floored(fPMJ * K0(0, 0), fPMJ * K0(0, 1),
                               fPMJ * K0(1, 0), fPMJ * K0(1, 1),
                               ktFloor, Fh0J[0][0], Fh0J[0][1], Fh0J[1][0], Fh0J[1][1]);
  }

  double Fs0Y = 0.0, Fs0Z = 0.0, ks0Y_raw = 0.0, ks0Z_raw = 0.0;
  if (hsY) {
    ks0Y_raw = mShearY->getInitialTangent();
    double ks0_rot = ks0Y_raw * Lc * Lc;
    if (fabs(ks0_rot) < ktFloor) ks0_rot = (ks0_rot < 0 ? -ktFloor : ktFloor);
    Fs0Y = 1.0 / ks0_rot;
  }
  if (hsZ) {
    ks0Z_raw = mShearZ->getInitialTangent();
    double ks0_rot = ks0Z_raw * Lc * Lc;
    if (fabs(ks0_rot) < ktFloor) ks0_rot = (ks0_rot < 0 ? -ktFloor : ktFloor);
    Fs0Z = 1.0 / ks0_rot;
  }
  const double gamRot0Y = hsY ? (-Fs0Y * (q0c[0] + q0c[1])) : 0.0;
  const double gamRot0Z = hsZ ? (-Fs0Z * (q0c[2] + q0c[3])) : 0.0;
  const double powellScale = usePowell ? 1.0 : 0.0;

  const bool matImplexI = ladrunoIMKBiaxialIsImplex(biI);
  const bool matImplexJ = ladrunoIMKBiaxialIsImplex(biJ);
  const bool matImplexSY = ladrunoIMKMatIsImplex(mShearY);
  const bool matImplexSZ = ladrunoIMKMatIsImplex(mShearZ);
  const bool matImplex = matImplexI || matImplexJ || matImplexSY || matImplexSZ;

  double th[4] = {
    hI ? th4[0] : 0.0,
    hJ ? th4[1] : 0.0,
    hI ? th4[2] : 0.0,
    hJ ? th4[3] : 0.0
  };
  double gamSY = (hsY && gammaSYInOut != 0) ? *gammaSYInOut : 0.0;
  double gamSZ = (hsZ && gammaSZInOut != 0) ? *gammaSZInOut : 0.0;

  double M[4] = {0.0, 0.0, 0.0, 0.0};
  double KI[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double KJ[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double VsY = 0.0, VsZ = 0.0, ktsY = 0.0, ktsZ = 0.0;
  Vector eps2(2);

  if (hI || hJ || hsY || hsZ) {
    if (useImplex) {
      if (matImplexI) ladrunoIMKBiaxialSetImplex(biI, false);
      if (matImplexJ) ladrunoIMKBiaxialSetImplex(biJ, false);
      if (matImplexSY) ladrunoIMKMatSetImplex(mShearY, false);
      if (matImplexSZ) ladrunoIMKMatSetImplex(mShearZ, false);
    }

    const int maxIter = 50;
    double damp = 1.0;
    double gnormPrev = 1.0e300;
    for (int it = 0; it < maxIter; it++) {
      if (hI) {
        eps2(0) = th[0]; eps2(1) = th[2];
        biI->setTrialStrain(eps2);
        const Vector &sig = biI->getStress();
        const Matrix &kt  = biI->getTangent();
        M[0] = fPMI * sig(0); M[2] = fPMI * sig(1);
        KI[0][0] = fPMI * kt(0, 0); KI[0][1] = fPMI * kt(0, 1);
        KI[1][0] = fPMI * kt(1, 0); KI[1][1] = fPMI * kt(1, 1);
      }
      if (hJ) {
        eps2(0) = th[1]; eps2(1) = th[3];
        biJ->setTrialStrain(eps2);
        const Vector &sig = biJ->getStress();
        const Matrix &kt  = biJ->getTangent();
        M[1] = fPMJ * sig(0); M[3] = fPMJ * sig(1);
        KJ[0][0] = fPMJ * kt(0, 0); KJ[0][1] = fPMJ * kt(0, 1);
        KJ[1][0] = fPMJ * kt(1, 0); KJ[1][1] = fPMJ * kt(1, 1);
      }
      if (hsY) {
        mShearY->setTrialStrain(gamSY);
        VsY = mShearY->getStress();
        ktsY = mShearY->getTangent();
      }
      if (hsZ) {
        mShearZ->setTrialStrain(gamSZ);
        VsZ = mShearZ->getStress();
        ktsZ = mShearZ->getTangent();
      }

      const double thCore[4] = {
        hI ? (th[0] - powellScale * (Fh0I[0][0] * M[0] + Fh0I[0][1] * M[2])) : 0.0,
        hJ ? (th[1] - powellScale * (Fh0J[0][0] * M[1] + Fh0J[0][1] * M[3])) : 0.0,
        hI ? (th[2] - powellScale * (Fh0I[1][0] * M[0] + Fh0I[1][1] * M[2])) : 0.0,
        hJ ? (th[3] - powellScale * (Fh0J[1][0] * M[1] + Fh0J[1][1] * M[3])) : 0.0
      };
      const double gamRotY = hsY ? (gamSY / Lc - gamRot0Y) : 0.0;
      const double gamRotZ = hsZ ? (gamSZ / Lc - gamRot0Z) : 0.0;
      const double vel0 = vc[0] - thCore[0] + gamRotY;
      const double vel1 = vc[1] - thCore[1] + gamRotY;
      const double vel2 = vc[2] - thCore[2] + gamRotZ;
      const double vel3 = vc[3] - thCore[3] + gamRotZ;

      const double qc[4] = {
        az * vel0 + bz * vel1 + q0c[0],
        bz * vel0 + az * vel1 + q0c[1],
        ay * vel2 + by * vel3 + q0c[2],
        by * vel2 + ay * vel3 + q0c[3]
      };

      const double g[4] = {
        hI ? (qc[0] - M[0]) : 0.0,
        hJ ? (qc[1] - M[1]) : 0.0,
        hI ? (qc[2] - M[2]) : 0.0,
        hJ ? (qc[3] - M[3]) : 0.0
      };
      const double gsY_mom = hsY ? (-(qc[0] + qc[1]) - Lc * VsY) : 0.0;
      const double gsZ_mom = hsZ ? (-(qc[2] + qc[3]) - Lc * VsZ) : 0.0;

      const double gnorm = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]
                              + gsY_mom*gsY_mom + gsZ_mom*gsZ_mom);
      const double ref = fabs(qc[0]) + fabs(qc[1]) + fabs(qc[2]) + fabs(qc[3])
                       + fabs(M[0]) + fabs(M[1]) + fabs(M[2]) + fabs(M[3])
                       + (hsY ? Lc * fabs(VsY) : 0.0) + (hsZ ? Lc * fabs(VsZ) : 0.0)
                       + fabs(q0c[0]) + fabs(q0c[1]) + fabs(q0c[2]) + fabs(q0c[3]);
      if (gnorm <= 1.0e-10 * ref + 1.0e-12)
        break;

      if (gnorm > gnormPrev) damp *= 0.5;
      else                   damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
      gnormPrev = gnorm;

      double az_eff = az, bz_eff = bz, DsY = 1.0;
      double ay_eff = ay, by_eff = by, DsZ = 1.0;
      const double abz = az + bz;
      const double aby = ay + by;
      double g_eff[4] = { g[0], g[1], g[2], g[3] };

      if (hsY) {
        const double ks_rot = ktsY * Lc * Lc;
        const double ks = (fabs(ks_rot) < ktFloor) ? (ks_rot < 0 ? -ktFloor : ktFloor) : ks_rot;
        DsY = 2.0 * abz + ks;
        const double dsFloor = 1.0e-8 * (2.0 * fabs(abz) + fabs(ks)) + 1.0e-300;
        if (fabs(DsY) < dsFloor) DsY = (DsY < 0.0) ? -dsFloor : dsFloor;
        az_eff = az - (abz * abz) / DsY;
        bz_eff = bz - (abz * abz) / DsY;
        if (hI) g_eff[0] += (abz / DsY) * gsY_mom;
        if (hJ) g_eff[1] += (abz / DsY) * gsY_mom;
      }
      if (hsZ) {
        const double ks_rot = ktsZ * Lc * Lc;
        const double ks = (fabs(ks_rot) < ktFloor) ? (ks_rot < 0 ? -ktFloor : ktFloor) : ks_rot;
        DsZ = 2.0 * aby + ks;
        const double dsFloor = 1.0e-8 * (2.0 * fabs(aby) + fabs(ks)) + 1.0e-300;
        if (fabs(DsZ) < dsFloor) DsZ = (DsZ < 0.0) ? -dsFloor : dsFloor;
        ay_eff = ay - (aby * aby) / DsZ;
        by_eff = by - (aby * aby) / DsZ;
        if (hI) g_eff[2] += (aby / DsZ) * gsZ_mom;
        if (hJ) g_eff[3] += (aby / DsZ) * gsZ_mom;
      }

      const double EI_mat[2][2] = {
        { hI ? (1.0 - powellScale * (Fh0I[0][0] * KI[0][0] + Fh0I[0][1] * KI[1][0])) : 0.0,
          hI ? (-powellScale * (Fh0I[0][0] * KI[0][1] + Fh0I[0][1] * KI[1][1]))      : 0.0 },
        { hI ? (-powellScale * (Fh0I[1][0] * KI[0][0] + Fh0I[1][1] * KI[1][0]))      : 0.0,
          hI ? (1.0 - powellScale * (Fh0I[1][0] * KI[0][1] + Fh0I[1][1] * KI[1][1])) : 0.0 }
      };
      const double EJ_mat[2][2] = {
        { hJ ? (1.0 - powellScale * (Fh0J[0][0] * KJ[0][0] + Fh0J[0][1] * KJ[1][0])) : 0.0,
          hJ ? (-powellScale * (Fh0J[0][0] * KJ[0][1] + Fh0J[0][1] * KJ[1][1]))      : 0.0 },
        { hJ ? (-powellScale * (Fh0J[1][0] * KJ[0][0] + Fh0J[1][1] * KJ[1][0]))      : 0.0,
          hJ ? (1.0 - powellScale * (Fh0J[1][0] * KJ[0][1] + Fh0J[1][1] * KJ[1][1])) : 0.0 }
      };

      // Build 4x4 system A * dth = g_eff, where A = Keff * E + KH
      double A[4][4] = {
        { az_eff * EI_mat[0][0] + KI[0][0], bz_eff * EJ_mat[0][0],            az_eff * EI_mat[0][1] + KI[0][1], bz_eff * EJ_mat[0][1]            },
        { bz_eff * EI_mat[0][0],            az_eff * EJ_mat[0][0] + KJ[0][0], bz_eff * EI_mat[0][1],            az_eff * EJ_mat[0][1] + KJ[0][1] },
        { ay_eff * EI_mat[1][0] + KI[1][0], by_eff * EJ_mat[1][0],            ay_eff * EI_mat[1][1] + KI[1][1], by_eff * EJ_mat[1][1]            },
        { by_eff * EI_mat[1][0],            ay_eff * EJ_mat[1][0] + KJ[1][0], by_eff * EI_mat[1][1],            ay_eff * EJ_mat[1][1] + KJ[1][1] }
      };
      if (!hI) {
        for (int m = 0; m < 4; m++) {
          A[0][m] = A[m][0] = 0.0;
          A[2][m] = A[m][2] = 0.0;
        }
        A[0][0] = 1.0; A[2][2] = 1.0;
        g_eff[0] = 0.0; g_eff[2] = 0.0;
      }
      if (!hJ) {
        for (int m = 0; m < 4; m++) {
          A[1][m] = A[m][1] = 0.0;
          A[3][m] = A[m][3] = 0.0;
        }
        A[1][1] = 1.0; A[3][3] = 1.0;
        g_eff[1] = 0.0; g_eff[3] = 0.0;
      }

      double Ainv[4][4];
      ladrunoIMKInvert4x4(A, Ainv, ktFloor);
      double dth[4] = {0.0, 0.0, 0.0, 0.0};
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          dth[i] += Ainv[i][j] * g_eff[j];

      for (int i = 0; i < 4; i++)
        th[i] += damp * dth[i];

      const double dthCore0 = EI_mat[0][0] * dth[0] + EI_mat[0][1] * dth[2];
      const double dthCore1 = EJ_mat[0][0] * dth[1] + EJ_mat[0][1] * dth[3];
      const double dthCore2 = EI_mat[1][0] * dth[0] + EI_mat[1][1] * dth[2];
      const double dthCore3 = EJ_mat[1][0] * dth[1] + EJ_mat[1][1] * dth[3];
      if (hsY) {
        const double dgamRotY = (gsY_mom + abz * (dthCore0 + dthCore1)) / DsY;
        gamSY += damp * (dgamRotY * Lc);
      }
      if (hsZ) {
        const double dgamRotZ = (gsZ_mom + aby * (dthCore2 + dthCore3)) / DsZ;
        gamSZ += damp * (dgamRotZ * Lc);
      }
    }
  }

  for (int i = 0; i < 4; i++)
    th4[i] = th[i];
  if (gammaSYInOut != 0) *gammaSYInOut = gamSY;
  if (gammaSZInOut != 0) *gammaSZInOut = gamSZ;

  const double thCore[4] = {
    hI ? (th[0] - powellScale * (Fh0I[0][0] * M[0] + Fh0I[0][1] * M[2])) : 0.0,
    hJ ? (th[1] - powellScale * (Fh0J[0][0] * M[1] + Fh0J[0][1] * M[3])) : 0.0,
    hI ? (th[2] - powellScale * (Fh0I[1][0] * M[0] + Fh0I[1][1] * M[2])) : 0.0,
    hJ ? (th[3] - powellScale * (Fh0J[1][0] * M[1] + Fh0J[1][1] * M[3])) : 0.0
  };
  const double gamRotY = hsY ? (gamSY / Lc - gamRot0Y) : 0.0;
  const double gamRotZ = hsZ ? (gamSZ / Lc - gamRot0Z) : 0.0;
  double qc[4] = {
    az * (vc[0] - thCore[0] + gamRotY) + bz * (vc[1] - thCore[1] + gamRotY) + q0c[0],
    bz * (vc[0] - thCore[0] + gamRotY) + az * (vc[1] - thCore[1] + gamRotY) + q0c[1],
    ay * (vc[2] - thCore[2] + gamRotZ) + by * (vc[3] - thCore[3] + gamRotZ) + q0c[2],
    by * (vc[2] - thCore[2] + gamRotZ) + ay * (vc[3] - thCore[3] + gamRotZ) + q0c[3]
  };

  const double thp[4] = {
    hI ? (th[0] - (Fh0I[0][0] * M[0] + Fh0I[0][1] * M[2])) : 0.0,
    hJ ? (th[1] - (Fh0J[0][0] * M[1] + Fh0J[0][1] * M[3])) : 0.0,
    hI ? (th[2] - (Fh0I[1][0] * M[0] + Fh0I[1][1] * M[2])) : 0.0,
    hJ ? (th[3] - (Fh0J[1][0] * M[1] + Fh0J[1][1] * M[3])) : 0.0
  };
  if (thpOut != 0) {
    for (int i = 0; i < 4; i++)
      thpOut[i] = thp[i];
  }
  const double gamSpY = (hsY && fabs(ks0Y_raw) > 1.0e-30) ? (gamSY - VsY / ks0Y_raw) : 0.0;
  const double gamSpZ = (hsZ && fabs(ks0Z_raw) > 1.0e-30) ? (gamSZ - VsZ / ks0Z_raw) : 0.0;
  if (gammaSpYOut != 0) *gammaSpYOut = gamSpY;
  if (gammaSpZOut != 0) *gammaSpZOut = gamSpZ;

  if (useImplex && thpCommit != 0 && dThpPrev != 0) {
    const double thpEx[4] = {
      hI ? (thpCommit[0] + dtScale * dThpPrev[0]) : 0.0,
      hJ ? (thpCommit[1] + dtScale * dThpPrev[1]) : 0.0,
      hI ? (thpCommit[2] + dtScale * dThpPrev[2]) : 0.0,
      hJ ? (thpCommit[3] + dtScale * dThpPrev[3]) : 0.0
    };
    const double gamSpYEx = hsY ? (gammaSpYCommit + dtScale * dGammaSpYPrev) : 0.0;
    const double gamSpZEx = hsZ ? (gammaSpZCommit + dtScale * dGammaSpZPrev) : 0.0;

    const double f0Scale = usePowell ? 0.0 : 1.0;
    double F0[4][4] = {
      { F11z + f0Scale * Fh0I[0][0] + Fs0Y, F12z + Fs0Y,                        f0Scale * Fh0I[0][1],               0.0                                 },
      { F12z + Fs0Y,                        F11z + f0Scale * Fh0J[0][0] + Fs0Y, 0.0,                                f0Scale * Fh0J[0][1]                },
      { f0Scale * Fh0I[1][0],               0.0,                                F11y + f0Scale * Fh0I[1][1] + Fs0Z, F12y + Fs0Z                         },
      { 0.0,                                f0Scale * Fh0J[1][0],               F12y + Fs0Z,                        F11y + f0Scale * Fh0J[1][1] + Fs0Z  }
    };
    double Kc0[4][4];
    ladrunoIMKInvert4x4(F0, Kc0, 1.0e-30);

    const double rhs[4] = {
      vc[0] - thpEx[0] + (hsY ? (gamSpYEx / Lc - gamRot0Y) : 0.0) + F11z * q0c[0] + F12z * q0c[1],
      vc[1] - thpEx[1] + (hsY ? (gamSpYEx / Lc - gamRot0Y) : 0.0) + F12z * q0c[0] + F11z * q0c[1],
      vc[2] - thpEx[2] + (hsZ ? (gamSpZEx / Lc - gamRot0Z) : 0.0) + F11y * q0c[2] + F12y * q0c[3],
      vc[3] - thpEx[3] + (hsZ ? (gamSpZEx / Lc - gamRot0Z) : 0.0) + F12y * q0c[2] + F11y * q0c[3]
    };
    const double dthp[4] = {
      (thp[0] - thpEx[0]) - (hsY ? (gamSpY - gamSpYEx) / Lc : 0.0),
      (thp[1] - thpEx[1]) - (hsY ? (gamSpY - gamSpYEx) / Lc : 0.0),
      (thp[2] - thpEx[2]) - (hsZ ? (gamSpZ - gamSpZEx) / Lc : 0.0),
      (thp[3] - thpEx[3]) - (hsZ ? (gamSpZ - gamSpZEx) / Lc : 0.0)
    };
    double err = 0.0;
    for (int i = 0; i < 4; i++) {
      double qex = 0.0;
      double dq_imp = 0.0;
      for (int j = 0; j < 4; j++) {
        qex += Kc0[i][j] * rhs[j];
        dq_imp += Kc0[i][j] * dthp[j];
      }
      err += fabs(dq_imp);
      qc[i] = qex;
    }
    if (implexErrOut != 0) *implexErrOut = err;
    q4[0] = T00 * qc[0] + T10 * qc[1];
    q4[1] = T01 * qc[0] + T11 * qc[1];
    q4[2] = T00 * qc[2] + T10 * qc[3];
    q4[3] = T01 * qc[2] + T11 * qc[3];
    if (dq4_dN != 0) {
      if (!usePowell) {
        const double dLogI = hI ? ladrunoIMKBiaxialDLogPM(biI, dLogPM) : 0.0;
        const double dLogJ = hJ ? ladrunoIMKBiaxialDLogPM(biJ, dLogPM) : 0.0;
        const double drhs[4] = {
          dLogI * (Fh0I[0][0] * qc[0] + Fh0I[0][1] * qc[2]),
          dLogJ * (Fh0J[0][0] * qc[1] + Fh0J[0][1] * qc[3]),
          dLogI * (Fh0I[1][0] * qc[0] + Fh0I[1][1] * qc[2]),
          dLogJ * (Fh0J[1][0] * qc[1] + Fh0J[1][1] * qc[3])
        };
        double dqc_dN[4] = { 0.0, 0.0, 0.0, 0.0 };
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            dqc_dN[i] += Kc0[i][j] * drhs[j];
        dq4_dN[0] = T00 * dqc_dN[0] + T10 * dqc_dN[1];
        dq4_dN[1] = T01 * dqc_dN[0] + T11 * dqc_dN[1];
        dq4_dN[2] = T00 * dqc_dN[2] + T10 * dqc_dN[3];
        dq4_dN[3] = T01 * dqc_dN[2] + T11 * dqc_dN[3];
      } else {
        for (int i = 0; i < 4; i++)
          dq4_dN[i] = 0.0;
      }
    }
    ladrunoIMKInitBiaxial(biI, biJ, L, EIz, EIy, k4, GAvy, GAvz, fPM, dI, dJ, mShearY, mShearZ, usePowell);
    return;
  }

  if (implexErrOut != 0) {
    *implexErrOut = matImplex
      ? (ladrunoIMKBiaxialImplexError(biI) + ladrunoIMKBiaxialImplexError(biJ)
         + ladrunoIMKMatImplexError(mShearY) + ladrunoIMKMatImplexError(mShearZ))
      : 0.0;
  }

  q4[0] = T00 * qc[0] + T10 * qc[1];
  q4[1] = T01 * qc[0] + T11 * qc[1];
  q4[2] = T00 * qc[2] + T10 * qc[3];
  q4[3] = T01 * qc[2] + T11 * qc[3];

  // Stiffness Schur-complement condensed tangent:
  //   K_c = K_eff - (K_eff * E) * (K_eff * E + K_H)^{-1} * (P_H * K_eff)
  // This never inverts K_I, K_J, or ktsY/ktsZ directly, so it remains exact
  // when a hinge tangent is zero (complete fracture) or negative (softening).
  double az_eff = az, bz_eff = bz;
  double ay_eff = ay, by_eff = by;
  if (hsY) {
    const double abz = az + bz;
    const double ks_rot = ktsY * Lc * Lc;
    double DsY = 2.0 * abz + ks_rot;
    const double dsFloor = 1.0e-8 * (2.0 * fabs(abz) + fabs(ks_rot)) + 1.0e-300;
    if (fabs(DsY) < dsFloor) DsY = (DsY < 0.0) ? -dsFloor : dsFloor;
    az_eff = az - (abz * abz) / DsY;
    bz_eff = bz - (abz * abz) / DsY;
  }
  if (hsZ) {
    const double aby = ay + by;
    const double ks_rot = ktsZ * Lc * Lc;
    double DsZ = 2.0 * aby + ks_rot;
    const double dsFloor = 1.0e-8 * (2.0 * fabs(aby) + fabs(ks_rot)) + 1.0e-300;
    if (fabs(DsZ) < dsFloor) DsZ = (DsZ < 0.0) ? -dsFloor : dsFloor;
    ay_eff = ay - (aby * aby) / DsZ;
    by_eff = by - (aby * aby) / DsZ;
  }

  const double Keff[4][4] = {
    { az_eff, bz_eff, 0.0,    0.0    },
    { bz_eff, az_eff, 0.0,    0.0    },
    { 0.0,    0.0,    ay_eff, by_eff },
    { 0.0,    0.0,    by_eff, ay_eff }
  };

  double Kc[4][4];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      Kc[i][j] = Keff[i][j];

  double dqc_dN[4] = { 0.0, 0.0, 0.0, 0.0 };
  if (hI || hJ) {
    const double EI_mat[2][2] = {
      { hI ? (1.0 - powellScale * (Fh0I[0][0] * KI[0][0] + Fh0I[0][1] * KI[1][0])) : 0.0,
        hI ? (-powellScale * (Fh0I[0][0] * KI[0][1] + Fh0I[0][1] * KI[1][1]))      : 0.0 },
      { hI ? (-powellScale * (Fh0I[1][0] * KI[0][0] + Fh0I[1][1] * KI[1][0]))      : 0.0,
        hI ? (1.0 - powellScale * (Fh0I[1][0] * KI[0][1] + Fh0I[1][1] * KI[1][1])) : 0.0 }
    };
    const double EJ_mat[2][2] = {
      { hJ ? (1.0 - powellScale * (Fh0J[0][0] * KJ[0][0] + Fh0J[0][1] * KJ[1][0])) : 0.0,
        hJ ? (-powellScale * (Fh0J[0][0] * KJ[0][1] + Fh0J[0][1] * KJ[1][1]))      : 0.0 },
      { hJ ? (-powellScale * (Fh0J[1][0] * KJ[0][0] + Fh0J[1][1] * KJ[1][0]))      : 0.0,
        hJ ? (1.0 - powellScale * (Fh0J[1][0] * KJ[0][1] + Fh0J[1][1] * KJ[1][1])) : 0.0 }
    };
    const double E_mat[4][4] = {
      { EI_mat[0][0], 0.0,          EI_mat[0][1], 0.0          },
      { 0.0,          EJ_mat[0][0], 0.0,          EJ_mat[0][1] },
      { EI_mat[1][0], 0.0,          EI_mat[1][1], 0.0          },
      { 0.0,          EJ_mat[1][0], 0.0,          EJ_mat[1][1] }
    };
    double KeffE[4][4];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        KeffE[i][j] = 0.0;
        for (int m = 0; m < 4; m++)
          KeffE[i][j] += Keff[i][m] * E_mat[m][j];
      }
    double A[4][4] = {
      { KeffE[0][0] + KI[0][0], KeffE[0][1],            KeffE[0][2] + KI[0][1], KeffE[0][3]            },
      { KeffE[1][0],            KeffE[1][1] + KJ[0][0], KeffE[1][2],            KeffE[1][3] + KJ[0][1] },
      { KeffE[2][0] + KI[1][0], KeffE[2][1],            KeffE[2][2] + KI[1][1], KeffE[2][3]            },
      { KeffE[3][0],            KeffE[3][1] + KJ[1][0], KeffE[3][2],            KeffE[3][3] + KJ[1][1] }
    };
    const double ph[4] = {
      hI ? 1.0 : 0.0,
      hJ ? 1.0 : 0.0,
      hI ? 1.0 : 0.0,
      hJ ? 1.0 : 0.0
    };
    if (!hI) {
      for (int m = 0; m < 4; m++) {
        A[0][m] = A[m][0] = 0.0;
        A[2][m] = A[m][2] = 0.0;
      }
      A[0][0] = 1.0; A[2][2] = 1.0;
    }
    if (!hJ) {
      for (int m = 0; m < 4; m++) {
        A[1][m] = A[m][1] = 0.0;
        A[3][m] = A[m][3] = 0.0;
      }
      A[1][1] = 1.0; A[3][3] = 1.0;
    }
    double Ainv[4][4];
    ladrunoIMKInvert4x4(A, Ainv, ktFloor);
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        double sub = 0.0;
        for (int m = 0; m < 4; m++) {
          if (ph[m] == 0.0) continue;
          for (int n = 0; n < 4; n++) {
            if (ph[n] == 0.0) continue;
            sub += KeffE[i][m] * Ainv[m][n] * Keff[n][j];
          }
        }
        Kc[i][j] -= sub;
      }
    }
    const double dLogI = ladrunoIMKBiaxialDLogPM(biI, dLogPM);
    const double dLogJ = ladrunoIMKBiaxialDLogPM(biJ, dLogPM);
    const double dM_dN[4] = {
      hI ? dLogI * M[0] : 0.0,
      hJ ? dLogJ * M[1] : 0.0,
      hI ? dLogI * M[2] : 0.0,
      hJ ? dLogJ * M[3] : 0.0
    };
    for (int i = 0; i < 4; i++) {
      double sum = 0.0;
      for (int m = 0; m < 4; m++) {
        if (ph[m] == 0.0) continue;
        for (int n = 0; n < 4; n++) {
          if (ph[n] == 0.0) continue;
          sum += KeffE[i][m] * Ainv[m][n] * dM_dN[n];
        }
      }
      dqc_dN[i] = sum;
    }
  }
  if (dq4_dN != 0) {
    dq4_dN[0] = T00 * dqc_dN[0] + T10 * dqc_dN[1];
    dq4_dN[1] = T01 * dqc_dN[0] + T11 * dqc_dN[1];
    dq4_dN[2] = T00 * dqc_dN[2] + T10 * dqc_dN[3];
    dq4_dN[3] = T01 * dqc_dN[2] + T11 * dqc_dN[3];
  }

  if (hasEZ) {
    const double T[4][4] = {
      { T00, T01, 0.0, 0.0 },
      { T10, T11, 0.0, 0.0 },
      { 0.0, 0.0, T00, T01 },
      { 0.0, 0.0, T10, T11 }
    };
    double KT[4][4];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        KT[i][j] = 0.0;
        for (int m = 0; m < 4; m++)
          KT[i][j] += Kc[i][m] * T[m][j];
      }
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        k4[i][j] = 0.0;
        for (int m = 0; m < 4; m++)
          k4[i][j] += T[m][i] * KT[m][j];
      }
  } else {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        k4[i][j] = Kc[i][j];
  }
}

#endif


