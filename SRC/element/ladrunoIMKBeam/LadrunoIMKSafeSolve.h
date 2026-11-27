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

// Author: N. Mora-Bowen (Ladruno), 11/2026 (WP-261)
//
// LadrunoIMKSafeSolve.h -- the safeguarded local solve of the LadrunoIMKBeam
// series hinges, kept free of OpenSees types (only <math.h>) so a standalone
// g++ harness can exercise it (tests/test_ladrunoIMK_safe_solve_cpp.py).
// LadrunoIMKHinge.h includes it and calls it when the plain local Newton
// exits unconverged.

#ifndef LadrunoIMKSafeSolve_h
#define LadrunoIMKSafeSolve_h

#include <math.h>

// Dense n x n solve (n <= 6) by Gaussian elimination with partial pivoting and
// a pivot floor. A and r are overwritten; the solution is returned in r.
inline void ladrunoIMKSolveDense(int n, double A[6][6], double r[6], double pivotFloor)
{
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int i = c + 1; i < n; i++)
      if (fabs(A[i][c]) > fabs(A[piv][c])) piv = i;
    if (piv != c) {
      for (int j = 0; j < n; j++) { const double t = A[c][j]; A[c][j] = A[piv][j]; A[piv][j] = t; }
      const double t = r[c]; r[c] = r[piv]; r[piv] = t;
    }
    if (fabs(A[c][c]) < pivotFloor) A[c][c] = (A[c][c] < 0.0) ? -pivotFloor : pivotFloor;
    for (int i = c + 1; i < n; i++) {
      const double f = A[i][c] / A[c][c];
      if (f == 0.0) continue;
      for (int j = c; j < n; j++) A[i][j] -= f * A[c][j];
      r[i] -= f * r[c];
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    double s = r[i];
    for (int j = i + 1; j < n; j++) s -= A[i][j] * r[j];
    r[i] = s / A[i][i];
  }
}

// Ladruno (WP-261): safeguarded fallback for the local hinge solve.
//
// Why: the plain local Newton cycles when a hinge must UNLOAD from a plastic
// committed state. With a bilinear law (Steel01, b = 1e-4, k0 ~ 150 a) the
// hinge residual g = a (v - th) - M(th) is monotone but S-shaped: slope ~a on
// the plastic branches, a + k0 on the narrow elastic band. A Newton step taken
// with the flat plastic slope overshoots the band onto the opposite plastic
// branch and the flat slope there sends it back: the arctan 2-cycle.
//
// What: Newton on G(x) = 0 with a BRACKETED LINE SEARCH on the directional
// residual h(t) = G(x + t d) . d. For hinge laws with a nonnegative tangent
// (more generally while K_el + K_H stays positive definite) G is minus the
// gradient of the convex series potential
//   Pi(x) = 1/2 e.K.e + q0.e + sum Phi_i(th_i) + Psi(gamma_s),
// so h(t) = -dPi/dt is nonincreasing in t. When the full Newton step passes
// the line minimum (h(1) < 0) the root of h in (0, 1) is bracketed and found
// by a safeguarded 1D solve (Newton in t when it stays inside the bracket,
// else Illinois regula falsi, else bisection): on a piecewise-linear law it
// lands on the elastic band in a few material calls instead of cycling. If
// the Newton direction is not a descent direction (indefinite tangent from a
// softening law) the steepest-descent direction d = G is used instead. Only
// the trial strain of each material is used, so any UniaxialMaterial or
// NDMaterial law works; there is nothing bilinear-specific here.
//
// Contract: ev(x, G, H) sets every material's trial state at x, fills the
// residual G (moment units) and the stiffness H = -dG/dx (inactive unknowns:
// G = 0, identity row), and returns true when x satisfies the convergence
// test. The LAST call to ev is always at the returned x, so the materials'
// trial states, moments and tangents are consistent with it whether or not
// it converged. Returns true when converged.
template <class Eval>
inline bool ladrunoIMKSafeguardedSolve(int n, double x[6], Eval &ev, double pivotFloor)
{
  const int maxOuter = 100;
  const int maxInner = 50;
  const double eta = 0.01;   // inner stop: |h(t)| <= eta h(0)
  double G[6], H[6][6], A[6][6], d[6], xt[6], Gt[6], Ht[6][6];
  bool conv = ev(x, G, H);
  for (int outer = 0; outer < maxOuter && !conv; outer++) {
    for (int i = 0; i < n; i++) {
      d[i] = G[i];
      for (int j = 0; j < n; j++) A[i][j] = H[i][j];
    }
    ladrunoIMKSolveDense(n, A, d, pivotFloor);
    double h0 = 0.0;
    for (int i = 0; i < n; i++) h0 += G[i] * d[i];
    if (!(h0 > 0.0)) {                 // not a descent direction: steepest descent
      h0 = 0.0;
      for (int i = 0; i < n; i++) { d[i] = G[i]; h0 += G[i] * G[i]; }
      if (!(h0 > 0.0)) break;          // G == 0 or NaN: nothing more to do
    }
    // full Newton step
    for (int i = 0; i < n; i++) xt[i] = x[i] + d[i];
    conv = ev(xt, Gt, Ht);
    double h1 = 0.0;
    for (int i = 0; i < n; i++) h1 += Gt[i] * d[i];
    if (!conv && h1 < 0.0) {
      // overshoot: the line minimum lies in (0, 1); bracket and solve h(t) = 0
      double tL = 0.0, hL = h0, tR = 1.0, hR = h1, tc = 1.0, hc = h1;
      double dpc = 0.0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) dpc -= d[i] * Ht[i][j] * d[j];
      int side = 0;
      for (int inner = 0; inner < maxInner; inner++) {
        double tn = (dpc < 0.0) ? (tc - hc / dpc) : -1.0;
        if (!(tn > tL && tn < tR)) {
          tn = (tL * hR - tR * hL) / (hR - hL);        // regula falsi
          if (!(tn > tL && tn < tR)) tn = 0.5 * (tL + tR);
        }
        for (int i = 0; i < n; i++) xt[i] = x[i] + tn * d[i];
        conv = ev(xt, Gt, Ht);
        double hn = 0.0;
        for (int i = 0; i < n; i++) hn += Gt[i] * d[i];
        tc = tn; hc = hn;
        if (conv || fabs(hn) <= eta * h0) break;
        if (hn > 0.0) {
          tL = tn; hL = hn;
          if (side == +1) hR *= 0.5;                   // Illinois
          side = +1;
        } else {
          tR = tn; hR = hn;
          if (side == -1) hL *= 0.5;
          side = -1;
        }
        if (tR - tL <= 1.0e-15) break;
        dpc = 0.0;
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++) dpc -= d[i] * Ht[i][j] * d[j];
      }
    }
    // accept the last evaluated point (keeps the materials consistent with x)
    for (int i = 0; i < n; i++) {
      x[i] = xt[i];
      G[i] = Gt[i];
      for (int j = 0; j < n; j++) H[i][j] = Ht[i][j];
    }
  }
  return conv;
}

#endif
