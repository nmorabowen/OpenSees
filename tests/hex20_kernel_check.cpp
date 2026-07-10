// Standalone driver for the pure H20 serendipity kernel — NO OpenSees build.
// Includes ONLY the pure headers LadrunoHex20Shape.h (shape/quadrature/B/mass/K)
// and, with LADRUNO_MASSLUMPING_STANDALONE, LadrunoMassLumping.h (the shared HRZ
// lumper). Prints %.17g tagged lines that the companion pytest
// (test_hex20_kernel_cpp.py) diffs against the sympy/numpy oracle
// (hex20_reference.py) — catching any transcription bug in the hand-coded C++
// serendipity algebra without building OpenSees. Mirrors the mechanics of
// tests/finitestrain2d_kernel_check.cpp.
//
// Sample points, distorted-hex corners and the elastic D block are DUPLICATED
// verbatim from hex20_reference.py (XI_SAMPLES / CORNERS_DISTORT / EDGES /
// iso_D(E=100, nu=0.25)) — keep in sync.

#include <cstdio>
#include <cmath>
#include <LadrunoHex20Shape.h>
#include <LadrunoMassLumping.h>

using namespace Ladruno::hex20;

// element edges (corner pairs) in mid-edge node order 8..19 (== oracle EDGES)
static const int EDGES[12][2] = {
  {0, 1}, {1, 2}, {2, 3}, {3, 0},
  {4, 5}, {5, 6}, {6, 7}, {7, 4},
  {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

// 20-node coord array from 8 corners: mid-edge nodes at TRUE edge midpoints
// (== oracle straight_edge_nodes).
static void straightEdge(const double corners[8][3], double X[20][3]) {
  for (int i = 0; i < 8; ++i)
    for (int c = 0; c < 3; ++c) X[i][c] = corners[i][c];
  for (int k = 0; k < 12; ++k) {
    const int i = EDGES[k][0], j = EDGES[k][1];
    for (int c = 0; c < 3; ++c) X[8 + k][c] = 0.5 * (X[i][c] + X[j][c]);
  }
}

// distorted straight-edged hex — DUPLICATE of oracle CORNERS_DISTORT
static const double CORNERS_DISTORT[8][3] = {
  {0.00, 0.00, 0.00}, {2.10, -0.10, 0.15}, {2.30, 2.05, -0.10}, {-0.15, 1.90, 0.10},
  {0.10, 0.20, 2.20}, {2.00, 0.00, 1.90}, {2.20, 2.30, 2.35}, {0.05, 2.10, 2.00}
};

// shared natural-point samples — DUPLICATE of oracle XI_SAMPLES
static const double XI_SAMPLES[5][3] = {
  {0.0, 0.0, 0.0},
  {0.3, -0.2, 0.7},
  {-0.9, 0.55, -0.35},
  {1.0, 1.0, 1.0},
  {0.123, 0.456, -0.789}
};

// isotropic elastic D (6x6, engineering-Voigt) — DUPLICATE of oracle iso_D
static void isoD(double D[6][6], double E, double nu) {
  const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double G = E / (2.0 * (1.0 + nu));
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) D[i][j] = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) D[i][j] = lam;
  for (int i = 0; i < 3; ++i) {
    D[i][i] += 2.0 * G;
    D[3 + i][3 + i] = G;
  }
}

int main() {
  // ---- node & GP tables (as the header defines them) ---------------------- //
  for (int a = 0; a < 20; ++a)
    for (int c = 0; c < 3; ++c)
      printf("NODEXI %d %d %.17g\n", a, c, NODE_XI[a][c]);
  for (int L = 0; L < 27; ++L)
    for (int c = 0; c < 4; ++c)
      printf("GP27 %d %d %.17g\n", L, c, GP27[L][c]);
  for (int L = 0; L < 8; ++L)
    for (int c = 0; c < 4; ++c)
      printf("GP8 %d %d %.17g\n", L, c, GP8[L][c]);

  // ---- N and dN at the 5 shared samples ----------------------------------- //
  for (int s = 0; s < 5; ++s) {
    double N[20], dN[20][3];
    shape(XI_SAMPLES[s], N, dN);
    for (int a = 0; a < 20; ++a) {
      printf("N %d %d %.17g\n", s, a, N[a]);
      for (int c = 0; c < 3; ++c)
        printf("DN %d %d %d %.17g\n", s, a, c, dN[a][c]);
    }
  }

  // ---- distorted hex: detJ at samples, B at samples 0 & 1, K27 & K8 -------- //
  double Xd[20][3];
  straightEdge(CORNERS_DISTORT, Xd);
  for (int s = 0; s < 5; ++s) {
    double N[20], dN[20][3], J[3][3];
    shape(XI_SAMPLES[s], N, dN);
    const double detJ = jacobian(Xd, dN, J);
    printf("DETJ %d %.17g\n", s, detJ);
    if (s == 0 || s == 1) {
      double Jinv[3][3], dNdx[20][3], B[6][60];
      invert3(J, detJ, Jinv);
      cartGrad(dN, Jinv, dNdx);
      fillB(dNdx, B);
      for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 60; ++c)
          printf("B %d %d %d %.17g\n", s, r, c, B[r][c]);
    }
  }

  double D[6][6];
  isoD(D, 100.0, 0.25);
  double K27[60 * 60], K8[60 * 60];
  stiffnessElastic(Xd, D, GP27, 27, K27);
  stiffnessElastic(Xd, D, GP8, 8, K8);
  for (int i = 0; i < 60; ++i)
    for (int j = 0; j < 60; ++j) {
      printf("K27 %d %d %.17g\n", i, j, K27[i * 60 + j]);
      printf("K8 %d %d %.17g\n", i, j, K8[i * 60 + j]);
    }

  // ---- reference cube: consistent M, row-sum, HRZ, volumes ---------------- //
  double Xc[8][3];
  for (int i = 0; i < 8; ++i)
    for (int c = 0; c < 3; ++c) Xc[i][c] = NODE_XI[i][c];   // cube corners = first 8 nodes
  double Xcube[20][3];
  straightEdge(Xc, Xcube);

  double M[60 * 60];
  consistentMass(Xcube, 1.0, GP27, 27, M);
  for (int i = 0; i < 60; ++i)
    for (int j = 0; j < 60; ++j)
      printf("M %d %d %.17g\n", i, j, M[i * 60 + j]);

  double rs[60];
  rowSumLump(M, rs);
  for (int i = 0; i < 60; ++i)
    printf("ROWSUM %d %.17g\n", i, rs[i]);

  int dir[60];
  dofDirs(dir);
  double hrz[60];
  Ladruno::hrzLumpRaw(M, dir, 60, hrz);
  for (int i = 0; i < 60; ++i)
    printf("HRZ %d %.17g\n", i, hrz[i]);

  printf("VOL27 %.17g\n", volume(Xcube, GP27, 27));
  printf("VOL8 %.17g\n", volume(Xcube, GP8, 8));
  return 0;
}
