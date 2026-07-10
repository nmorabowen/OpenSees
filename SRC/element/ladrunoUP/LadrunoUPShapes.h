/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoUPShapes — stateless shape providers for the Biot u-p family (ADR 71,
// companion to LadrunoUPKernel.h). Pure doubles + <math.h>, NO OpenSees deps, so
// the whole coupled kernel is unit-tested standalone (tests/ladruno_up_kernel_
// check.cpp) against the numpy oracle. One struct per geometry — T3, Q4, H8,
// BT6 (Bernstein-quadratic triangle), BTET10 (Bernstein-quadratic tetrahedron) —
// each exposing static-const GP data plus stateless eval / evalP / elemSizeH.
//
// Shape functions and Jacobians are LIFTED (and pinned to) the fork's own donor
// elements, cited per struct:
//   T3      LadrunoCST.cpp:655–688      (N0=ξ, N1=η, N2=1−ξ−η)
//   Q4      LadrunoQuad.cpp:1549–1605   (standard CCW bilinear)
//   H8      SRC/element/brick/shp3d.cpp (standard CCW trilinear)
//   BT6     BezierTri6.cpp:851–896 + GP3 (~:109)
//   BTET10  BezierTet10.cpp:1177–1236 + GP4 (~:101)
//
// PINNED GP rules & basis orders (ADR-71 plan 0.A addendum §2–4 — verified here
// against the donors, all agree): T3 = 1-pt (⅓,⅓) w=½; Q4 = 2×2 Gauss (±1/√3)
// w=1; H8 = 2×2×2; BT6 = 3-pt interior {(⅙,⅙),(⅔,⅙),(⅙,⅔)} w=⅙; BTET10 = 4-pt
// perms of a=0.585410196624968 / b=0.138196601125011 w=1/24.
//
// LAYOUTS: nodal coords xy[a*ndm + i] (row-major, node a component i). Outputs:
//   Nu[a], dNu[a*ndm + i] = ∂N_a/∂x_i (cartesian, node-major), *detJ (scalar).
//   Np[b], dNp[b*ndm + i] — the pressure basis:
//     equal-order (taylorHood=false) => Np ≡ Nu (nNp = nNodes);
//     Taylor–Hood (taylorHood=true)  => barycentric-linear L_i on the VERTEX
//     carriers, COMPACTED to nNp = #carriers, ordered by element node order
//     restricted to carriers (BT6: [ξ₃,ξ₁,ξ₂]; BTET10: [L1,L2,L3,L4]). This is
//     the linear barycentric basis, NOT the Bernstein vertex subset (ADR §3.3
//     ⟨UP-8⟩). Linear shapes ignore the flag (all nodes are vertices).
//
// PRECONDITION (Bézier providers): STRAIGHT SIDES / affine geometry map — the
// mid-edge control points at the edge midpoints, so the quadratic Bernstein map
// reduces to a constant-Jacobian affine map and ∇L_i are exact constants (the
// TH-stability requirement, ADR §3.3). The element-level straight-side guard
// arrives in P3; these providers ASSUME it and are silent on violation.
//
// GUARDS: none here (the *Kernel.h idiom) — a singular / inverted Jacobian
// yields garbage gradients silently; the consuming element must reject detJ ≤ 0.

#ifndef LadrunoUPShapes_h
#define LadrunoUPShapes_h

#include <math.h>

namespace ladruno_up {
namespace shapes {

// ---- small dense-inverse helpers (row-major) ------------------------------ //
// 2×2: J[0..3] = {J00,J01,J10,J11}; returns det, fills Jm1 = inv(J) row-major.
inline double inv2x2(const double J[4], double Jm1[4]) {
  double det = J[0] * J[3] - J[1] * J[2];
  double id = (det != 0.0) ? 1.0 / det : 0.0;
  Jm1[0] =  J[3] * id;  Jm1[1] = -J[1] * id;
  Jm1[2] = -J[2] * id;  Jm1[3] =  J[0] * id;
  return det;
}
// 3×3: J[0..8] row-major; returns det, fills Jm1 = inv(J) row-major.
inline double inv3x3(const double J[9], double Jm1[9]) {
  double c00 = J[4] * J[8] - J[5] * J[7];
  double c01 = J[5] * J[6] - J[3] * J[8];
  double c02 = J[3] * J[7] - J[4] * J[6];
  double det = J[0] * c00 + J[1] * c01 + J[2] * c02;
  double id = (det != 0.0) ? 1.0 / det : 0.0;
  // Jm1 = adjugate^T / det
  Jm1[0] = c00 * id;
  Jm1[1] = (J[2] * J[7] - J[1] * J[8]) * id;
  Jm1[2] = (J[1] * J[5] - J[2] * J[4]) * id;
  Jm1[3] = c01 * id;
  Jm1[4] = (J[0] * J[8] - J[2] * J[6]) * id;
  Jm1[5] = (J[2] * J[3] - J[0] * J[5]) * id;
  Jm1[6] = c02 * id;
  Jm1[7] = (J[1] * J[6] - J[0] * J[7]) * id;
  Jm1[8] = (J[0] * J[4] - J[1] * J[3]) * id;
  return det;
}

// largest vertex-pair distance over the listed vertex node indices (ADR §3.3).
inline double maxVertexPairDist(const double* xy, int ndm,
                                const int* vtx, int nvtx) {
  double hmax = 0.0;
  for (int p = 0; p < nvtx; p++)
    for (int q = p + 1; q < nvtx; q++) {
      double d2 = 0.0;
      for (int i = 0; i < ndm; i++) {
        double dd = xy[vtx[p] * ndm + i] - xy[vtx[q] * ndm + i];
        d2 += dd * dd;
      }
      if (d2 > hmax) hmax = d2;
    }
  return sqrt(hmax);
}

// =========================================================================== //
//  T3 — linear triangle (donor: LadrunoCST.cpp:655–688)                       //
//  N0=ξ, N1=η, N2=1−ξ−η. 1-pt rule at (⅓,⅓), w=½. Equal-order p only.         //
// =========================================================================== //
struct T3 {
  static const int nNodes = 3, ndm = 2, nGP = 1;
  static constexpr double gpCoord[nGP * ndm] = {1.0 / 3.0, 1.0 / 3.0};
  static constexpr double gpWeight[nGP] = {0.5};
  static constexpr int pCarrierTH[nNodes] = {1, 1, 1};

  static void basis(double xi, double eta, double N[3],
                    double dxi[3], double deta[3]) {
    N[0] = xi;  N[1] = eta;  N[2] = 1.0 - xi - eta;
    dxi[0]  = 1.0;  dxi[1]  = 0.0;  dxi[2]  = -1.0;
    deta[0] = 0.0;  deta[1] = 1.0;  deta[2] = -1.0;
  }
  static double jacobian(const double* xy, const double dxi[3],
                         const double deta[3], double Jm1[4]) {
    double J[4] = {0, 0, 0, 0};   // {∂x/∂ξ, ∂x/∂η, ∂y/∂ξ, ∂y/∂η}
    for (int a = 0; a < nNodes; a++) {
      J[0] += dxi[a]  * xy[a * ndm + 0];
      J[1] += deta[a] * xy[a * ndm + 0];
      J[2] += dxi[a]  * xy[a * ndm + 1];
      J[3] += deta[a] * xy[a * ndm + 1];
    }
    return inv2x2(J, Jm1);
  }
  static void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ) {
    double xi = gpCoord[gp * ndm + 0], eta = gpCoord[gp * ndm + 1];
    double N[3], dxi[3], deta[3];
    basis(xi, eta, N, dxi, deta);
    double Jm1[4];
    *detJ = jacobian(xy, dxi, deta, Jm1);
    for (int a = 0; a < nNodes; a++) {
      Nu[a] = N[a];
      dNu[a * ndm + 0] = Jm1[0] * dxi[a] + Jm1[2] * deta[a];   // ∂N/∂x
      dNu[a * ndm + 1] = Jm1[1] * dxi[a] + Jm1[3] * deta[a];   // ∂N/∂y
    }
  }
  static void evalP(int gp, const double* xy, bool taylorHood,
                    double* Np, double* dNp) {
    (void)taylorHood;                       // linear: all nodes are vertices
    double detJ;
    eval(gp, xy, Np, dNp, &detJ);           // p ≡ u
  }
  static double elemSizeH(const double* xy) {
    static const int vtx[3] = {0, 1, 2};
    return maxVertexPairDist(xy, ndm, vtx, 3);
  }
};

// =========================================================================== //
//  Q4 — bilinear quad (donor: LadrunoQuad.cpp:1549–1605)                      //
//  standard CCW; 2×2 Gauss (±1/√3), w=1. Equal-order p.                       //
// =========================================================================== //
struct Q4 {
  static const int nNodes = 4, ndm = 2, nGP = 4;
  static constexpr double G = 0.5773502691896258;   // 1/√3
  static constexpr double gpCoord[nGP * ndm] = {
      -G, -G,   G, -G,   G, G,   -G, G};
  static constexpr double gpWeight[nGP] = {1.0, 1.0, 1.0, 1.0};
  static constexpr int pCarrierTH[nNodes] = {1, 1, 1, 1};

  static void basis(double xi, double eta, double N[4],
                    double dxi[4], double deta[4]) {
    double xm = 1.0 - xi, xp = 1.0 + xi, em = 1.0 - eta, ep = 1.0 + eta;
    N[0] = 0.25 * xm * em;  N[1] = 0.25 * xp * em;
    N[2] = 0.25 * xp * ep;  N[3] = 0.25 * xm * ep;
    dxi[0]  = -0.25 * em;  dxi[1]  =  0.25 * em;  dxi[2]  =  0.25 * ep;  dxi[3]  = -0.25 * ep;
    deta[0] = -0.25 * xm;  deta[1] = -0.25 * xp;  deta[2] =  0.25 * xp;  deta[3] =  0.25 * xm;
  }
  static double jacobian(const double* xy, const double dxi[4],
                         const double deta[4], double Jm1[4]) {
    double J[4] = {0, 0, 0, 0};
    for (int a = 0; a < nNodes; a++) {
      J[0] += dxi[a]  * xy[a * ndm + 0];
      J[1] += deta[a] * xy[a * ndm + 0];
      J[2] += dxi[a]  * xy[a * ndm + 1];
      J[3] += deta[a] * xy[a * ndm + 1];
    }
    return inv2x2(J, Jm1);
  }
  static void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ) {
    double xi = gpCoord[gp * ndm + 0], eta = gpCoord[gp * ndm + 1];
    double N[4], dxi[4], deta[4];
    basis(xi, eta, N, dxi, deta);
    double Jm1[4];
    *detJ = jacobian(xy, dxi, deta, Jm1);
    for (int a = 0; a < nNodes; a++) {
      Nu[a] = N[a];
      dNu[a * ndm + 0] = Jm1[0] * dxi[a] + Jm1[2] * deta[a];
      dNu[a * ndm + 1] = Jm1[1] * dxi[a] + Jm1[3] * deta[a];
    }
  }
  static void evalP(int gp, const double* xy, bool taylorHood,
                    double* Np, double* dNp) {
    (void)taylorHood;
    double detJ;
    eval(gp, xy, Np, dNp, &detJ);
  }
  static double elemSizeH(const double* xy) {
    static const int vtx[4] = {0, 1, 2, 3};
    return maxVertexPairDist(xy, ndm, vtx, 4);
  }
};

// =========================================================================== //
//  H8 — trilinear hexahedron (donor: SRC/element/brick/shp3d.cpp)             //
//  standard CCW; 2×2×2 Gauss, w=1. Equal-order p.                             //
// =========================================================================== //
struct H8 {
  static const int nNodes = 8, ndm = 3, nGP = 8;
  static constexpr double G = 0.5773502691896258;
  static constexpr double gpCoord[nGP * ndm] = {
      -G, -G, -G,    G, -G, -G,    G, G, -G,    -G, G, -G,
      -G, -G,  G,    G, -G,  G,    G, G,  G,    -G, G,  G};
  static constexpr double gpWeight[nGP] = {1, 1, 1, 1, 1, 1, 1, 1};
  static constexpr int pCarrierTH[nNodes] = {1, 1, 1, 1, 1, 1, 1, 1};
  // natural corner coords (shp3d node order)
  static constexpr double nc[nNodes * 3] = {
      -1, -1, -1,   1, -1, -1,   1, 1, -1,   -1, 1, -1,
      -1, -1,  1,   1, -1,  1,   1, 1,  1,   -1, 1,  1};

  static void basis(double xi, double eta, double zeta, double N[8],
                    double d[3][8]) {
    for (int a = 0; a < nNodes; a++) {
      double xa = nc[a * 3 + 0], ya = nc[a * 3 + 1], za = nc[a * 3 + 2];
      double gx = 1.0 + xi * xa, gy = 1.0 + eta * ya, gz = 1.0 + zeta * za;
      N[a]    = 0.125 * gx * gy * gz;
      d[0][a] = 0.125 * xa * gy * gz;   // ∂N/∂ξ
      d[1][a] = 0.125 * ya * gx * gz;   // ∂N/∂η
      d[2][a] = 0.125 * za * gx * gy;   // ∂N/∂ζ
    }
  }
  static double jacobian(const double* xy, const double d[3][8], double Jm1[9]) {
    double J[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   // J[i*3+j] = ∂x_i/∂ξ_j
    for (int a = 0; a < nNodes; a++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          J[i * 3 + j] += d[j][a] * xy[a * ndm + i];
    return inv3x3(J, Jm1);
  }
  static void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ) {
    double xi = gpCoord[gp * ndm + 0], eta = gpCoord[gp * ndm + 1],
           zeta = gpCoord[gp * ndm + 2];
    double N[8], d[3][8];
    basis(xi, eta, zeta, N, d);
    double Jm1[9];
    *detJ = jacobian(xy, d, Jm1);
    for (int a = 0; a < nNodes; a++) {
      Nu[a] = N[a];
      for (int i = 0; i < 3; i++)      // ∂N_a/∂x_i = Σ_j dNdξ_j · Jm1[j][i]
        dNu[a * ndm + i] = d[0][a] * Jm1[0 * 3 + i]
                         + d[1][a] * Jm1[1 * 3 + i]
                         + d[2][a] * Jm1[2 * 3 + i];
    }
  }
  static void evalP(int gp, const double* xy, bool taylorHood,
                    double* Np, double* dNp) {
    (void)taylorHood;
    double detJ;
    eval(gp, xy, Np, dNp, &detJ);
  }
  static double elemSizeH(const double* xy) {
    static const int vtx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    return maxVertexPairDist(xy, ndm, vtx, 8);
  }
};

// =========================================================================== //
//  BT6 — Bernstein-quadratic triangle (donor: BezierTri6.cpp:851–896)         //
//  N0=ξ₃², N1=ξ₁², N2=ξ₂², N3=2ξ₁ξ₃, N4=2ξ₁ξ₂, N5=2ξ₂ξ₃  (ξ₃=1−ξ₁−ξ₂).       //
//  vertices = nodes 0,1,2 (node0↔ξ₃, node1↔ξ₁, node2↔ξ₂). 3-pt rule.          //
//  TH pressure carriers = nodes 0,1,2 → Np = [ξ₃, ξ₁, ξ₂] (linear barycentric).//
//  STRAIGHT-SIDE PRECONDITION (affine map). gpCoord = (ξ₁,ξ₂).                //
// =========================================================================== //
struct BT6 {
  static const int nNodes = 6, ndm = 2, nGP = 3;
  static constexpr double gpCoord[nGP * ndm] = {
      1.0 / 6.0, 1.0 / 6.0,   2.0 / 3.0, 1.0 / 6.0,   1.0 / 6.0, 2.0 / 3.0};
  static constexpr double gpWeight[nGP] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
  static constexpr int pCarrierTH[nNodes] = {1, 1, 1, 0, 0, 0};

  static void basisU(double x1, double x2, double N[6],
                     double dx1[6], double dx2[6]) {
    double x3 = 1.0 - x1 - x2;
    N[0] = x3 * x3;   N[1] = x1 * x1;   N[2] = x2 * x2;
    N[3] = 2.0 * x1 * x3;   N[4] = 2.0 * x1 * x2;   N[5] = 2.0 * x2 * x3;
    dx1[0] = -2.0 * x3;   dx1[1] = 2.0 * x1;   dx1[2] = 0.0;
    dx1[3] = 2.0 * (x3 - x1);   dx1[4] = 2.0 * x2;   dx1[5] = -2.0 * x2;
    dx2[0] = -2.0 * x3;   dx2[1] = 0.0;   dx2[2] = 2.0 * x2;
    dx2[3] = -2.0 * x1;   dx2[4] = 2.0 * x1;   dx2[5] = 2.0 * (x3 - x2);
  }
  static double jacobian(const double* xy, const double dx1[6],
                         const double dx2[6], double Jm1[4]) {
    double J[4] = {0, 0, 0, 0};
    for (int a = 0; a < nNodes; a++) {
      J[0] += dx1[a] * xy[a * ndm + 0];
      J[1] += dx2[a] * xy[a * ndm + 0];
      J[2] += dx1[a] * xy[a * ndm + 1];
      J[3] += dx2[a] * xy[a * ndm + 1];
    }
    return inv2x2(J, Jm1);
  }
  static void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ) {
    double x1 = gpCoord[gp * ndm + 0], x2 = gpCoord[gp * ndm + 1];
    double N[6], dx1[6], dx2[6];
    basisU(x1, x2, N, dx1, dx2);
    double Jm1[4];
    *detJ = jacobian(xy, dx1, dx2, Jm1);
    for (int a = 0; a < nNodes; a++) {
      Nu[a] = N[a];
      dNu[a * ndm + 0] = Jm1[0] * dx1[a] + Jm1[2] * dx2[a];
      dNu[a * ndm + 1] = Jm1[1] * dx1[a] + Jm1[3] * dx2[a];
    }
  }
  static void evalP(int gp, const double* xy, bool taylorHood,
                    double* Np, double* dNp) {
    double x1 = gpCoord[gp * ndm + 0], x2 = gpCoord[gp * ndm + 1];
    double N[6], dx1[6], dx2[6];
    basisU(x1, x2, N, dx1, dx2);
    double Jm1[4];
    (void)jacobian(xy, dx1, dx2, Jm1);        // affine: constant J (straight sides)
    if (!taylorHood) {                         // equal-order: p ≡ u
      for (int a = 0; a < nNodes; a++) {
        Np[a] = N[a];
        dNp[a * ndm + 0] = Jm1[0] * dx1[a] + Jm1[2] * dx2[a];
        dNp[a * ndm + 1] = Jm1[1] * dx1[a] + Jm1[3] * dx2[a];
      }
      return;
    }
    // Taylor–Hood: linear barycentric p on vertices, compacted [ξ₃, ξ₁, ξ₂].
    double x3 = 1.0 - x1 - x2;
    double L[3]  = {x3, x1, x2};
    double dL1[3] = {-1.0, 1.0, 0.0};          // ∂L/∂ξ₁
    double dL2[3] = {-1.0, 0.0, 1.0};          // ∂L/∂ξ₂
    for (int b = 0; b < 3; b++) {
      Np[b] = L[b];
      dNp[b * ndm + 0] = Jm1[0] * dL1[b] + Jm1[2] * dL2[b];
      dNp[b * ndm + 1] = Jm1[1] * dL1[b] + Jm1[3] * dL2[b];
    }
  }
  static double elemSizeH(const double* xy) {
    static const int vtx[3] = {0, 1, 2};       // vertex nodes only
    return maxVertexPairDist(xy, ndm, vtx, 3);
  }
};

// =========================================================================== //
//  BTET10 — Bernstein-quadratic tetrahedron (donor: BezierTet10.cpp:1177–1236)//
//  N0..3 = L1²..L4²; N4=2L1L2, N5=2L2L3, N6=2L1L3, N7=2L1L4, N8=2L3L4,        //
//  N9=2L2L4  (L4=1−L1−L2−L3). vertices = nodes 0–3 (node k ↔ L_{k+1}).        //
//  4-pt rule; TH carriers = nodes 0–3 → Np = [L1,L2,L3,L4]. gpCoord=(L1,L2,L3).//
//  STRAIGHT-SIDE PRECONDITION (affine map).                                   //
// =========================================================================== //
struct BTET10 {
  static const int nNodes = 10, ndm = 3, nGP = 4;
  static constexpr double A = 0.585410196624968, B = 0.138196601125011;
  static constexpr double gpCoord[nGP * ndm] = {
      A, B, B,   B, A, B,   B, B, A,   B, B, B};
  static constexpr double gpWeight[nGP] = {
      1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0};
  static constexpr int pCarrierTH[nNodes] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0};

  static void basisU(double L1, double L2, double L3, double N[10],
                     double d[3][10]) {
    double L4 = 1.0 - L1 - L2 - L3;
    N[0] = L1 * L1;  N[1] = L2 * L2;  N[2] = L3 * L3;  N[3] = L4 * L4;
    N[4] = 2.0 * L1 * L2;  N[5] = 2.0 * L2 * L3;  N[6] = 2.0 * L1 * L3;
    N[7] = 2.0 * L1 * L4;  N[8] = 2.0 * L3 * L4;  N[9] = 2.0 * L2 * L4;
    // ∂N/∂L1
    d[0][0] = 2.0 * L1;  d[0][1] = 0.0;  d[0][2] = 0.0;  d[0][3] = -2.0 * L4;
    d[0][4] = 2.0 * L2;  d[0][5] = 0.0;  d[0][6] = 2.0 * L3;
    d[0][7] = 2.0 * (L4 - L1);  d[0][8] = -2.0 * L3;  d[0][9] = -2.0 * L2;
    // ∂N/∂L2
    d[1][0] = 0.0;  d[1][1] = 2.0 * L2;  d[1][2] = 0.0;  d[1][3] = -2.0 * L4;
    d[1][4] = 2.0 * L1;  d[1][5] = 2.0 * L3;  d[1][6] = 0.0;
    d[1][7] = -2.0 * L1;  d[1][8] = -2.0 * L3;  d[1][9] = 2.0 * (L4 - L2);
    // ∂N/∂L3
    d[2][0] = 0.0;  d[2][1] = 0.0;  d[2][2] = 2.0 * L3;  d[2][3] = -2.0 * L4;
    d[2][4] = 0.0;  d[2][5] = 2.0 * L2;  d[2][6] = 2.0 * L1;
    d[2][7] = -2.0 * L1;  d[2][8] = 2.0 * (L4 - L3);  d[2][9] = -2.0 * L2;
  }
  static double jacobian(const double* xy, const double d[3][10], double Jm1[9]) {
    double J[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   // J[i*3+j] = ∂x_i/∂L_j
    for (int a = 0; a < nNodes; a++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          J[i * 3 + j] += d[j][a] * xy[a * ndm + i];
    return inv3x3(J, Jm1);
  }
  static void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ) {
    double L1 = gpCoord[gp * ndm + 0], L2 = gpCoord[gp * ndm + 1],
           L3 = gpCoord[gp * ndm + 2];
    double N[10], d[3][10];
    basisU(L1, L2, L3, N, d);
    double Jm1[9];
    *detJ = jacobian(xy, d, Jm1);
    for (int a = 0; a < nNodes; a++) {
      Nu[a] = N[a];
      for (int i = 0; i < 3; i++)
        dNu[a * ndm + i] = d[0][a] * Jm1[0 * 3 + i]
                         + d[1][a] * Jm1[1 * 3 + i]
                         + d[2][a] * Jm1[2 * 3 + i];
    }
  }
  static void evalP(int gp, const double* xy, bool taylorHood,
                    double* Np, double* dNp) {
    double L1 = gpCoord[gp * ndm + 0], L2 = gpCoord[gp * ndm + 1],
           L3 = gpCoord[gp * ndm + 2];
    double N[10], d[3][10];
    basisU(L1, L2, L3, N, d);
    double Jm1[9];
    (void)jacobian(xy, d, Jm1);               // affine (straight sides)
    if (!taylorHood) {
      for (int a = 0; a < nNodes; a++) {
        Np[a] = N[a];
        for (int i = 0; i < 3; i++)
          dNp[a * ndm + i] = d[0][a] * Jm1[0 * 3 + i]
                           + d[1][a] * Jm1[1 * 3 + i]
                           + d[2][a] * Jm1[2 * 3 + i];
      }
      return;
    }
    // Taylor–Hood: linear barycentric p on vertices, compacted [L1,L2,L3,L4].
    double L4 = 1.0 - L1 - L2 - L3;
    double L[4] = {L1, L2, L3, L4};
    // ∂L/∂(L1,L2,L3) parametric (L4 = 1−L1−L2−L3)
    double dpar[4][3] = {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, -1, -1}};
    for (int b = 0; b < 4; b++) {
      Np[b] = L[b];
      for (int i = 0; i < 3; i++)
        dNp[b * ndm + i] = dpar[b][0] * Jm1[0 * 3 + i]
                         + dpar[b][1] * Jm1[1 * 3 + i]
                         + dpar[b][2] * Jm1[2 * 3 + i];
    }
  }
  static double elemSizeH(const double* xy) {
    static const int vtx[4] = {0, 1, 2, 3};    // vertex nodes only
    return maxVertexPairDist(xy, ndm, vtx, 4);
  }
};

} // namespace shapes
} // namespace ladruno_up

#endif
