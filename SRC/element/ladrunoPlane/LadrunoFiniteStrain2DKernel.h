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
// LadrunoFiniteStrain2DKernel — the PURE numerical core of the shared 2D (plane)
// finite-strain, updated-Lagrangian element assembly (ADR 70). Plain doubles +
// <math.h>, NO OpenSees dependency, so it is unit-tested standalone against the
// numpy oracle (tests/finitestrain2d_reference.py) WITHOUT building OpenSees, and
// reused by every 2D continuum element (LadrunoQuad / LadrunoCST / LadrunoLST /
// BezierTri6) — each supplies only its shape gradients + GP rule; the finite-strain
// mechanics live here, written once.
//
// Per Gauss point (indices i,j,k,l ∈ {0,1}; a,b node indices; the material returns
// the in-plane Cauchy σ (2×2) and the FULL in-plane spatial modulus c2[2][2][2][2]
// via FiniteStrainND2DMaterial — that dispatch is the ELEMENT's job, not the
// kernel's, so this header stays OpenSees-free):
//
//   F      = I + Σ_a u_a ⊗ ∂N_a/∂X                    (deformationGradient2D)
//   F-bar  : F̄ = (J₀/J)^{1/2} F   (n=2, NOT ⅓)         (fbarScale2D)
//   g_j^a  = Σ_k ∂N_a/∂X_k (F⁻¹)_kj                    (spatialGradients2D)
//   a_ijkl = c_ijkl − σ_il δ_jk   (dSNPO eq 14.99; NOT minor-sym in (k,l))
//   f_{a,i} += σ_ij g_j^a · dv                         (addInternalForce2D)
//   K_{(a,i)(b,k)} += g_j^a a_ijkl g_l^b · dv          (addTangent2D)
//     + [F-bar] (Σ_j g_j^a q_ij)(g0_k^b − g_k^b) dv,   (addFbarCoupling2D)
//       q_ij = (1/n) Σ_p a_ijpp − ((n-1)/n) σ_ij       (dSNPO eq 15.10/15.11)
//
// GUARDS ARE THE ELEMENT'S JOB (this header stays OpenSees-free and silent on
// failure): the consuming element MUST reject det F ≤ 0 at every GP AND J₀ ≤ 0 at
// the centroid (F-bar) — returning < 0 so the solver cuts the step — BEFORE calling
// these routines, exactly as LadrunoBrick::updateFinite does. inv2() on a singular F
// returns an all-zero inverse (⇒ zero gradients ⇒ spuriously zero residual AND
// tangent, no error), and fbarScale2D() returns NaN for J₀/J < 0 — so an element
// that omits the guards inherits a silent-NaN / zero-stiffness trap.
//
// References: de Souza Neto, Perić & Owen (2008) §14.5 (14.99 spatial tangent),
// §15.1 (F-bar). Faithful i,j,k,l∈{0,1} restriction of LadrunoBrick's proven 3D
// finite path (formResidAndTangentFinite).

#ifndef LadrunoFiniteStrain2DKernel_h
#define LadrunoFiniteStrain2DKernel_h

#include <math.h>

namespace ladruno_fs2d {

// ---- 2×2 helpers (row-major double[4]) --------------------------------- //
inline double det2(const double F[4]) { return F[0] * F[3] - F[1] * F[2]; }

inline double inv2(const double F[4], double Fi[4]) {
  double d = det2(F);
  double id = (d != 0.0) ? 1.0 / d : 0.0;
  Fi[0] =  F[3] * id;  Fi[1] = -F[1] * id;
  Fi[2] = -F[2] * id;  Fi[3] =  F[0] * id;
  return d;
}

// F = I + Σ_a u_a ⊗ ∂N_a/∂X.  gradRef[k*nN + a] = ∂N_a/∂X_k ; u[2a+i] = (u_a)_i.
// Returns J = det F. F row-major [F00,F01,F10,F11].
inline double deformationGradient2D(const double *gradRef, const double *u,
                                    int nN, double F[4]) {
  F[0] = 1.0; F[1] = 0.0; F[2] = 0.0; F[3] = 1.0;
  for (int a = 0; a < nN; a++) {
    double gx = gradRef[0 * nN + a], gy = gradRef[1 * nN + a];
    double ux = u[2 * a + 0], uy = u[2 * a + 1];
    F[0] += ux * gx;  F[1] += ux * gy;
    F[2] += uy * gx;  F[3] += uy * gy;
  }
  return det2(F);
}

// F-bar dilatation scale (plane, n=2): (J0/J)^{1/2}. Returns F̄ = s·F in Fbar[4].
inline double fbarScale2D(double J0, double J) {
  return sqrt(J0 / J);
}

// spatial gradients g_j^a = Σ_k ∂N_a/∂X_k (F⁻¹)_kj.  gSpatial[j*nN + a].
inline void spatialGradients2D(const double *gradRef, const double Finv[4],
                               int nN, double *gSpatial) {
  for (int a = 0; a < nN; a++) {
    double gx = gradRef[0 * nN + a], gy = gradRef[1 * nN + a];
    // (F⁻¹)_kj: Finv row-major [00,01,10,11]; g_j = Σ_k gradRef_k Finv_kj
    gSpatial[0 * nN + a] = gx * Finv[0] + gy * Finv[2];   // j=0: k=0->Finv[0], k=1->Finv[2]
    gSpatial[1 * nN + a] = gx * Finv[1] + gy * Finv[3];   // j=1
  }
}

// a_ijkl = c_ijkl − σ_il δ_jk  (full 2×2×2×2 spatial consistent modulus).
inline void spatialTangent2D(const double c2[2][2][2][2], const double sig[2][2],
                             double a4[2][2][2][2]) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
          a4[i][j][k][l] = c2[i][j][k][l] - sig[i][l] * (j == k ? 1.0 : 0.0);
}

// f_{a,i} += σ_ij g_j^a · dv   (accumulate into resid[2a+i]).
inline void addInternalForce2D(const double sig[2][2], const double *g, int nN,
                               double dv, double *resid) {
  for (int a = 0; a < nN; a++)
    for (int i = 0; i < 2; i++) {
      double fi = 0.0;
      for (int j = 0; j < 2; j++) fi += sig[i][j] * g[j * nN + a];
      resid[2 * a + i] += fi * dv;
    }
}

// K_{(a,i)(b,k)} += g_j^a a_ijkl g_l^b · dv   (accumulate into K, ld = 2*nN).
inline void addTangent2D(const double a4[2][2][2][2], const double *g, int nN,
                         double dv, double *K, int ld) {
  for (int a = 0; a < nN; a++)
    for (int b = 0; b < nN; b++)
      for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++) {
          double s = 0.0;
          for (int j = 0; j < 2; j++)
            for (int l = 0; l < 2; l++)
              s += g[j * nN + a] * a4[i][j][k][l] * g[l * nN + b];
          K[(2 * a + i) * ld + (2 * b + k)] += s * dv;
        }
}

// F-bar unsymmetric G₀ coupling (dSNPO eq 15.10/15.11, n=2):
//   q_ij = (1/2) Σ_p a4_ijpp − (1/2) σ_ij
//   K_{(a,i)(b,k)} += (Σ_j g_j^a q_ij)(g0_k^b − g_k^b) dv
// g = actual GP spatial gradients, g0 = centroid spatial gradients (both [j*nN+a]).
inline void addFbarCoupling2D(const double a4[2][2][2][2], const double sig[2][2],
                              const double *g, const double *g0, int nN,
                              double dv, double *K, int ld) {
  double q[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      q[i][j] = 0.5 * (a4[i][j][0][0] + a4[i][j][1][1]) - 0.5 * sig[i][j];

  for (int a = 0; a < nN; a++)
    for (int i = 0; i < 2; i++) {
      double p = 0.0;                             // p_i^a = Σ_j g_j^a q_ij
      for (int j = 0; j < 2; j++) p += g[j * nN + a] * q[i][j];
      for (int b = 0; b < nN; b++)
        for (int k = 0; k < 2; k++)
          K[(2 * a + i) * ld + (2 * b + k)]
            += p * (g0[k * nN + b] - g[k * nN + b]) * dv;
    }
}

} // namespace ladruno_fs2d

#endif
