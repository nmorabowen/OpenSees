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

// LadrunoCPMLShapes.h — Standalone header-only C++17 shape functions and
// quadrature rules for the 9 LadrunoCPML continuum topologies:
//   2D (ndm=2): T3, Q4, T6, BT6
//   3D (ndm=3): TET4, H8, TET10, BTET10, H20
// Zero OpenSees Element/Domain/Matrix/Vector dependencies.

#ifndef LadrunoCPMLShapes_h
#define LadrunoCPMLShapes_h

#include <cmath>
#include <array>
#include <algorithm>
#include <cassert>

namespace LadrunoCPML {

enum class CPMLTopology {
  T3 = 0,
  Q4 = 1,
  T6 = 2,
  BT6 = 3,
  TET4 = 4,
  H8 = 5,
  TET10 = 6,
  BTET10 = 7,
  H20 = 8
};

namespace shapes {

static constexpr int MAX_NODES = 20;
static constexpr int MAX_GP = 8;
static constexpr int MAX_NDM = 3;

// ---- 2x2 and 3x3 Jacobian inverse helpers -------------------------------- //
// J[p * ndm + q] = d x_q / d xi_p (row p = reference coord, col q = physical coord).
// Then dN_a/dxi_p = sum_q J[p, q] * dN_a/dx_q, so:
//   dN_a/dx_i = sum_p Jinv[i, p] * dN_a/dxi_p.
// Scale-free degeneracy guard (WP-115): normalize |detJ| by maxRow^ndm so
// micro-scale elements (L = 1e-6) are accepted and collapsed macro-scale
// elements (L = 100) are rejected.
inline double invertJacobian2D(const double J[4], double Jinv[4]) {
  const double r0 = std::hypot(J[0], J[1]);
  const double r1 = std::hypot(J[2], J[3]);
  const double maxRow = std::max(r0, r1);
  const double scale = maxRow * maxRow;
  const double det = J[0] * J[3] - J[1] * J[2];
  if (!(scale > 0.0) || !std::isfinite(det) || std::fabs(det) <= 1.0e-10 * scale) {
    Jinv[0] = Jinv[1] = Jinv[2] = Jinv[3] = 0.0;
    return 0.0;
  }
  const double id = 1.0 / det;
  Jinv[0] =  J[3] * id;
  Jinv[1] = -J[1] * id;
  Jinv[2] = -J[2] * id;
  Jinv[3] =  J[0] * id;
  return det;
}

inline double invertJacobian3D(const double J[9], double Jinv[9]) {
  const double r0 = std::sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2]);
  const double r1 = std::sqrt(J[3] * J[3] + J[4] * J[4] + J[5] * J[5]);
  const double r2 = std::sqrt(J[6] * J[6] + J[7] * J[7] + J[8] * J[8]);
  const double maxRow = std::max({r0, r1, r2});
  const double scale = maxRow * maxRow * maxRow;
  const double c00 = J[4] * J[8] - J[5] * J[7];
  const double c01 = J[5] * J[6] - J[3] * J[8];
  const double c02 = J[3] * J[7] - J[4] * J[6];
  const double det = J[0] * c00 + J[1] * c01 + J[2] * c02;
  if (!(scale > 0.0) || !std::isfinite(det) || std::fabs(det) <= 1.0e-10 * scale) {
    for (int i = 0; i < 9; ++i) Jinv[i] = 0.0;
    return 0.0;
  }
  const double id = 1.0 / det;
  Jinv[0] = c00 * id;
  Jinv[1] = (J[2] * J[7] - J[1] * J[8]) * id;
  Jinv[2] = (J[1] * J[5] - J[2] * J[4]) * id;
  Jinv[3] = c01 * id;
  Jinv[4] = (J[0] * J[8] - J[2] * J[6]) * id;
  Jinv[5] = (J[2] * J[3] - J[0] * J[5]) * id;
  Jinv[6] = c02 * id;
  Jinv[7] = (J[1] * J[6] - J[0] * J[7]) * id;
  Jinv[8] = (J[0] * J[4] - J[1] * J[3]) * id;
  return det;
}

// Generic 2D evaluation from parametric basis and derivatives dN_dxi[a][2]
template <int NNodes>
inline void eval2DFromBasis(const double* coords,
                            const double N[NNodes],
                            const double dN_dxi[NNodes][2],
                            double weight,
                            double* N_out,
                            double* dNdx_out,
                            double* x_gp_out,
                            double* detJ_out,
                            double* w_vol_out) {
  double J[4] = {0.0, 0.0, 0.0, 0.0};
  double x_gp[2] = {0.0, 0.0};
  for (int a = 0; a < NNodes; ++a) {
    const double xa = coords[a * 2 + 0];
    const double ya = coords[a * 2 + 1];
    x_gp[0] += N[a] * xa;
    x_gp[1] += N[a] * ya;
    J[0] += dN_dxi[a][0] * xa; // dx/dxi_0
    J[1] += dN_dxi[a][0] * ya; // dy/dxi_0
    J[2] += dN_dxi[a][1] * xa; // dx/dxi_1
    J[3] += dN_dxi[a][1] * ya; // dy/dxi_1
  }
  double Jinv[4];
  const double detJ = invertJacobian2D(J, Jinv);
  for (int a = 0; a < NNodes; ++a) {
    N_out[a] = N[a];
    dNdx_out[a * 2 + 0] = Jinv[0] * dN_dxi[a][0] + Jinv[1] * dN_dxi[a][1];
    dNdx_out[a * 2 + 1] = Jinv[2] * dN_dxi[a][0] + Jinv[3] * dN_dxi[a][1];
  }
  if (x_gp_out) {
    x_gp_out[0] = x_gp[0];
    x_gp_out[1] = x_gp[1];
  }
  if (detJ_out) *detJ_out = detJ;
  if (w_vol_out) *w_vol_out = (detJ == 0.0) ? 0.0 : (weight * detJ);
}

// Generic 3D evaluation from parametric basis and derivatives dN_dxi[a][3].
// When UseAbsDetJ == true (used for TET10 and BTET10 whose barycentric
// parameterization zeta_4 = 1 - zeta_1 - zeta_2 - zeta_3 is an odd permutation
// on standard right-handed tets), w_vol uses |detJ| while detJ_out retains the
// signed detJ so orientation uniformity across GPs can still be verified.
template <int NNodes, bool UseAbsDetJ = false>
inline void eval3DFromBasis(const double* coords,
                            const double N[NNodes],
                            const double dN_dxi[NNodes][3],
                            double weight,
                            double* N_out,
                            double* dNdx_out,
                            double* x_gp_out,
                            double* detJ_out,
                            double* w_vol_out) {
  double J[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double x_gp[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < NNodes; ++a) {
    const double xa = coords[a * 3 + 0];
    const double ya = coords[a * 3 + 1];
    const double za = coords[a * 3 + 2];
    x_gp[0] += N[a] * xa;
    x_gp[1] += N[a] * ya;
    x_gp[2] += N[a] * za;
    for (int p = 0; p < 3; ++p) {
      J[p * 3 + 0] += dN_dxi[a][p] * xa;
      J[p * 3 + 1] += dN_dxi[a][p] * ya;
      J[p * 3 + 2] += dN_dxi[a][p] * za;
    }
  }
  double Jinv[9];
  const double detJ = invertJacobian3D(J, Jinv);
  for (int a = 0; a < NNodes; ++a) {
    N_out[a] = N[a];
    for (int i = 0; i < 3; ++i) {
      dNdx_out[a * 3 + i] = Jinv[i * 3 + 0] * dN_dxi[a][0]
                          + Jinv[i * 3 + 1] * dN_dxi[a][1]
                          + Jinv[i * 3 + 2] * dN_dxi[a][2];
    }
  }
  if (x_gp_out) {
    x_gp_out[0] = x_gp[0];
    x_gp_out[1] = x_gp[1];
    x_gp_out[2] = x_gp[2];
  }
  if (detJ_out) *detJ_out = detJ;
  if (w_vol_out) *w_vol_out = (detJ == 0.0) ? 0.0 : (weight * (UseAbsDetJ ? std::fabs(detJ) : detJ));
}

// =========================================================================== //
// 1. T3 (ndm=2, nNodes=3, nGP=1)
//    Reference coords (r, s) with L_0 = 1 - r - s, L_1 = r, L_2 = s.
//    1 GP: (r, s) = (1/3, 1/3), w = 0.5, w_vol = 0.5 * detJ.
// =========================================================================== //
struct T3 {
  static constexpr int ndm = 2;
  static constexpr int nNodes = 3;
  static constexpr int nGP = 1;
  static constexpr double gpCoord[nGP * ndm] = {1.0 / 3.0, 1.0 / 3.0};
  static constexpr double gpWeight[nGP] = {0.5};

  static void basis(double r, double s, double N[3], double dN[3][2]) {
    N[0] = 1.0 - r - s;
    N[1] = r;
    N[2] = s;
    dN[0][0] = -1.0; dN[0][1] = -1.0;
    dN[1][0] =  1.0; dN[1][1] =  0.0;
    dN[2][0] =  0.0; dN[2][1] =  1.0;
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[3], dN[3][2];
    basis(gpCoord[gp * 2 + 0], gpCoord[gp * 2 + 1], N, dN);
    eval2DFromBasis<3>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                       x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 2. Q4 (ndm=2, nNodes=4, nGP=4)
//    Nodes 0..3 at (-1,-1), (1,-1), (1,1), (-1,1).
//    2x2 Gauss rule with g = 0.57735026918962576451, w_g = 1.0.
//    Loop order: j (eta -g..+g) outer, i (xi -g..+g) inner:
//      GP0=(-g,-g), GP1=(+g,-g), GP2=(-g,+g), GP3=(+g,+g).
// =========================================================================== //
struct Q4 {
  static constexpr int ndm = 2;
  static constexpr int nNodes = 4;
  static constexpr int nGP = 4;
  static constexpr double G = 0.57735026918962576451;
  static constexpr double gpCoord[nGP * ndm] = {
    -G, -G,
    +G, -G,
    -G, +G,
    +G, +G
  };
  static constexpr double gpWeight[nGP] = {1.0, 1.0, 1.0, 1.0};
  static constexpr double nodeXi[4][2] = {
    {-1.0, -1.0}, {+1.0, -1.0}, {+1.0, +1.0}, {-1.0, +1.0}
  };

  static void basis(double xi, double eta, double N[4], double dN[4][2]) {
    for (int a = 0; a < 4; ++a) {
      const double xa = nodeXi[a][0];
      const double ea = nodeXi[a][1];
      N[a] = 0.25 * (1.0 + xa * xi) * (1.0 + ea * eta);
      dN[a][0] = 0.25 * xa * (1.0 + ea * eta);
      dN[a][1] = 0.25 * ea * (1.0 + xa * xi);
    }
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[4], dN[4][2];
    basis(gpCoord[gp * 2 + 0], gpCoord[gp * 2 + 1], N, dN);
    eval2DFromBasis<4>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                       x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 3. T6 (ndm=2, nNodes=6, nGP=3 — Lagrange 6-node triangle)
//    Reference coords (s, t) where Node 0=(1,0), Node 1=(0,1), Node 2=(0,0),
//    L_3 = 1 - s - t:
//      N_0 = s(2s - 1), N_1 = t(2t - 1), N_2 = L_3(2L_3 - 1)
//      N_3 = 4st, N_4 = 4tL_3, N_5 = 4sL_3
//    3 GPs: (1/6, 1/6), (2/3, 1/6), (1/6, 2/3) with w_g = 1/6 each.
// =========================================================================== //
struct T6 {
  static constexpr int ndm = 2;
  static constexpr int nNodes = 6;
  static constexpr int nGP = 3;
  static constexpr double gpCoord[nGP * ndm] = {
    1.0 / 6.0, 1.0 / 6.0,
    2.0 / 3.0, 1.0 / 6.0,
    1.0 / 6.0, 2.0 / 3.0
  };
  static constexpr double gpWeight[nGP] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

  static void basis(double s, double t, double N[6], double dN[6][2]) {
    const double L3 = 1.0 - s - t;
    N[0] = s * (2.0 * s - 1.0);
    N[1] = t * (2.0 * t - 1.0);
    N[2] = L3 * (2.0 * L3 - 1.0);
    N[3] = 4.0 * s * t;
    N[4] = 4.0 * t * L3;
    N[5] = 4.0 * s * L3;

    dN[0][0] = 4.0 * s - 1.0;      dN[0][1] = 0.0;
    dN[1][0] = 0.0;                dN[1][1] = 4.0 * t - 1.0;
    dN[2][0] = 1.0 - 4.0 * L3;     dN[2][1] = 1.0 - 4.0 * L3;
    dN[3][0] = 4.0 * t;            dN[3][1] = 4.0 * s;
    dN[4][0] = -4.0 * t;           dN[4][1] = 4.0 * (L3 - t);
    dN[5][0] = 4.0 * (L3 - s);     dN[5][1] = -4.0 * s;
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[6], dN[6][2];
    basis(gpCoord[gp * 2 + 0], gpCoord[gp * 2 + 1], N, dN);
    eval2DFromBasis<6>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                       x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 4. BT6 (ndm=2, nNodes=6, nGP=3 — Bezier 6-node triangle)
//    Reference coords (l_1, l_2) with l_0 = 1 - l_1 - l_2:
//      N_0 = l_0^2, N_1 = l_1^2, N_2 = l_2^2,
//      N_3 = 2 l_0 l_1, N_4 = 2 l_1 l_2, N_5 = 2 l_2 l_0.
//    3 GPs: (2/3, 1/6), (1/6, 2/3), (1/6, 1/6) with w_g = 1/6 each.
// =========================================================================== //
struct BT6 {
  static constexpr int ndm = 2;
  static constexpr int nNodes = 6;
  static constexpr int nGP = 3;
  static constexpr double gpCoord[nGP * ndm] = {
    2.0 / 3.0, 1.0 / 6.0,
    1.0 / 6.0, 2.0 / 3.0,
    1.0 / 6.0, 1.0 / 6.0
  };
  static constexpr double gpWeight[nGP] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

  static void basis(double l1, double l2, double N[6], double dN[6][2]) {
    const double l0 = 1.0 - l1 - l2;
    N[0] = l0 * l0;
    N[1] = l1 * l1;
    N[2] = l2 * l2;
    N[3] = 2.0 * l0 * l1;
    N[4] = 2.0 * l1 * l2;
    N[5] = 2.0 * l2 * l0;

    dN[0][0] = -2.0 * l0;          dN[0][1] = -2.0 * l0;
    dN[1][0] =  2.0 * l1;          dN[1][1] =  0.0;
    dN[2][0] =  0.0;               dN[2][1] =  2.0 * l2;
    dN[3][0] =  2.0 * (l0 - l1);   dN[3][1] = -2.0 * l1;
    dN[4][0] =  2.0 * l2;          dN[4][1] =  2.0 * l1;
    dN[5][0] = -2.0 * l2;          dN[5][1] =  2.0 * (l0 - l2);
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[6], dN[6][2];
    basis(gpCoord[gp * 2 + 0], gpCoord[gp * 2 + 1], N, dN);
    eval2DFromBasis<6>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                       x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 5. TET4 (ndm=3, nNodes=4, nGP=1)
//    Reference coords (r, s, t) with L_0 = 1 - r - s - t, L_1 = r, L_2 = s, L_3 = t.
//    1 GP at (0.25, 0.25, 0.25) with w = 1/6, w_vol = (1/6) detJ.
// =========================================================================== //
struct TET4 {
  static constexpr int ndm = 3;
  static constexpr int nNodes = 4;
  static constexpr int nGP = 1;
  static constexpr double gpCoord[nGP * ndm] = {0.25, 0.25, 0.25};
  static constexpr double gpWeight[nGP] = {1.0 / 6.0};

  static void basis(double r, double s, double t, double N[4], double dN[4][3]) {
    N[0] = 1.0 - r - s - t;
    N[1] = r;
    N[2] = s;
    N[3] = t;
    dN[0][0] = -1.0; dN[0][1] = -1.0; dN[0][2] = -1.0;
    dN[1][0] =  1.0; dN[1][1] =  0.0; dN[1][2] =  0.0;
    dN[2][0] =  0.0; dN[2][1] =  1.0; dN[2][2] =  0.0;
    dN[3][0] =  0.0; dN[3][1] =  0.0; dN[3][2] =  1.0;
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[4], dN[4][3];
    basis(gpCoord[gp * 3 + 0], gpCoord[gp * 3 + 1], gpCoord[gp * 3 + 2], N, dN);
    eval3DFromBasis<4, false>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                              x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 6. H8 (ndm=3, nNodes=8, nGP=8)
//    Nodes 0..7 at (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
//                  (-1,-1,1),  (1,-1,1),  (1,1,1),  (-1,1,1).
//    2x2x2 Gauss rule with g = 0.57735026918962576451, w_g = 1.0.
//    Loop order: k (zeta -g..+g) outer, j (eta -g..+g) middle, i (xi -g..+g) inner.
// =========================================================================== //
struct H8 {
  static constexpr int ndm = 3;
  static constexpr int nNodes = 8;
  static constexpr int nGP = 8;
  static constexpr double G = 0.57735026918962576451;
  static constexpr double gpCoord[nGP * ndm] = {
    -G, -G, -G,
    +G, -G, -G,
    -G, +G, -G,
    +G, +G, -G,
    -G, -G, +G,
    +G, -G, +G,
    -G, +G, +G,
    +G, +G, +G
  };
  static constexpr double gpWeight[nGP] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  static constexpr double nodeXi[8][3] = {
    {-1.0, -1.0, -1.0}, {+1.0, -1.0, -1.0}, {+1.0, +1.0, -1.0}, {-1.0, +1.0, -1.0},
    {-1.0, -1.0, +1.0}, {+1.0, -1.0, +1.0}, {+1.0, +1.0, +1.0}, {-1.0, +1.0, +1.0}
  };

  static void basis(double xi, double eta, double zeta, double N[8], double dN[8][3]) {
    for (int a = 0; a < 8; ++a) {
      const double xa = nodeXi[a][0];
      const double ea = nodeXi[a][1];
      const double za = nodeXi[a][2];
      const double gx = 1.0 + xa * xi;
      const double ge = 1.0 + ea * eta;
      const double gz = 1.0 + za * zeta;
      N[a]     = 0.125 * gx * ge * gz;
      dN[a][0] = 0.125 * xa * ge * gz;
      dN[a][1] = 0.125 * ea * gx * gz;
      dN[a][2] = 0.125 * za * gx * ge;
    }
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[8], dN[8][3];
    basis(gpCoord[gp * 3 + 0], gpCoord[gp * 3 + 1], gpCoord[gp * 3 + 2], N, dN);
    eval3DFromBasis<8, false>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                              x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 7. TET10 (ndm=3, nNodes=10, nGP=4 — Lagrange 10-node tetrahedron)
//    Reference coords (z1, z2, z3) with z4 = 1 - z1 - z2 - z3:
//      N_0 = z1(2z1-1), N_1 = z2(2z2-1), N_2 = z3(2z3-1), N_3 = z4(2z4-1)
//      N_4 = 4z1z2, N_5 = 4z2z3, N_6 = 4z3z1, N_7 = 4z1z4, N_8 = 4z3z4, N_9 = 4z2z4
//    4 GPs with a = 0.5854101966249685, b = 0.1381966011250105, w_g = 1/24.
//    Note: Uses |detJ| for w_vol so both right-handed and left-handed vertex
//    orderings work seamlessly while detJ_out reports signed detJ.
// =========================================================================== //
struct TET10 {
  static constexpr int ndm = 3;
  static constexpr int nNodes = 10;
  static constexpr int nGP = 4;
  static constexpr double A = 0.5854101966249685;
  static constexpr double B = 0.1381966011250105;
  static constexpr double gpCoord[nGP * ndm] = {
    A, B, B,
    B, A, B,
    B, B, A,
    B, B, B
  };
  static constexpr double gpWeight[nGP] = {
    1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0
  };

  static void basis(double z1, double z2, double z3, double N[10], double dN[10][3]) {
    const double z4 = 1.0 - z1 - z2 - z3;
    N[0] = z1 * (2.0 * z1 - 1.0);
    N[1] = z2 * (2.0 * z2 - 1.0);
    N[2] = z3 * (2.0 * z3 - 1.0);
    N[3] = z4 * (2.0 * z4 - 1.0);
    N[4] = 4.0 * z1 * z2;
    N[5] = 4.0 * z2 * z3;
    N[6] = 4.0 * z3 * z1;
    N[7] = 4.0 * z1 * z4;
    N[8] = 4.0 * z3 * z4;
    N[9] = 4.0 * z2 * z4;

    dN[0][0] = 4.0 * z1 - 1.0;       dN[0][1] = 0.0;                  dN[0][2] = 0.0;
    dN[1][0] = 0.0;                  dN[1][1] = 4.0 * z2 - 1.0;       dN[1][2] = 0.0;
    dN[2][0] = 0.0;                  dN[2][1] = 0.0;                  dN[2][2] = 4.0 * z3 - 1.0;
    dN[3][0] = 1.0 - 4.0 * z4;       dN[3][1] = 1.0 - 4.0 * z4;       dN[3][2] = 1.0 - 4.0 * z4;
    dN[4][0] = 4.0 * z2;             dN[4][1] = 4.0 * z1;             dN[4][2] = 0.0;
    dN[5][0] = 0.0;                  dN[5][1] = 4.0 * z3;             dN[5][2] = 4.0 * z2;
    dN[6][0] = 4.0 * z3;             dN[6][1] = 0.0;                  dN[6][2] = 4.0 * z1;
    dN[7][0] = 4.0 * (z4 - z1);      dN[7][1] = -4.0 * z1;            dN[7][2] = -4.0 * z1;
    dN[8][0] = -4.0 * z3;            dN[8][1] = -4.0 * z3;            dN[8][2] = 4.0 * (z4 - z3);
    dN[9][0] = -4.0 * z2;            dN[9][1] = 4.0 * (z4 - z2);      dN[9][2] = -4.0 * z2;
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[10], dN[10][3];
    basis(gpCoord[gp * 3 + 0], gpCoord[gp * 3 + 1], gpCoord[gp * 3 + 2], N, dN);
    eval3DFromBasis<10, true>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                              x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 8. BTET10 (ndm=3, nNodes=10, nGP=4 — Bezier 10-node tetrahedron)
//    Reference coords (z1, z2, z3) with z4 = 1 - z1 - z2 - z3:
//      N_0 = z1^2, N_1 = z2^2, N_2 = z3^2, N_3 = z4^2
//      N_4 = 2z1z2, N_5 = 2z2z3, N_6 = 2z3z1, N_7 = 2z1z4, N_8 = 2z3z4, N_9 = 2z2z4
//    Same 4 GPs and weights w_g = 1/24 as TET10.
//    Note: Uses |detJ| for w_vol so both right-handed and left-handed vertex
//    orderings work seamlessly while detJ_out reports signed detJ.
// =========================================================================== //
struct BTET10 {
  static constexpr int ndm = 3;
  static constexpr int nNodes = 10;
  static constexpr int nGP = 4;
  static constexpr double A = 0.5854101966249685;
  static constexpr double B = 0.1381966011250105;
  static constexpr double gpCoord[nGP * ndm] = {
    A, B, B,
    B, A, B,
    B, B, A,
    B, B, B
  };
  static constexpr double gpWeight[nGP] = {
    1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0
  };

  static void basis(double z1, double z2, double z3, double N[10], double dN[10][3]) {
    const double z4 = 1.0 - z1 - z2 - z3;
    N[0] = z1 * z1;
    N[1] = z2 * z2;
    N[2] = z3 * z3;
    N[3] = z4 * z4;
    N[4] = 2.0 * z1 * z2;
    N[5] = 2.0 * z2 * z3;
    N[6] = 2.0 * z3 * z1;
    N[7] = 2.0 * z1 * z4;
    N[8] = 2.0 * z3 * z4;
    N[9] = 2.0 * z2 * z4;

    dN[0][0] = 2.0 * z1;             dN[0][1] = 0.0;                  dN[0][2] = 0.0;
    dN[1][0] = 0.0;                  dN[1][1] = 2.0 * z2;             dN[1][2] = 0.0;
    dN[2][0] = 0.0;                  dN[2][1] = 0.0;                  dN[2][2] = 2.0 * z3;
    dN[3][0] = -2.0 * z4;            dN[3][1] = -2.0 * z4;            dN[3][2] = -2.0 * z4;
    dN[4][0] = 2.0 * z2;             dN[4][1] = 2.0 * z1;             dN[4][2] = 0.0;
    dN[5][0] = 0.0;                  dN[5][1] = 2.0 * z3;             dN[5][2] = 2.0 * z2;
    dN[6][0] = 2.0 * z3;             dN[6][1] = 0.0;                  dN[6][2] = 2.0 * z1;
    dN[7][0] = 2.0 * (z4 - z1);      dN[7][1] = -2.0 * z1;            dN[7][2] = -2.0 * z1;
    dN[8][0] = -2.0 * z3;            dN[8][1] = -2.0 * z3;            dN[8][2] = 2.0 * (z4 - z3);
    dN[9][0] = -2.0 * z2;            dN[9][1] = 2.0 * (z4 - z2);      dN[9][2] = -2.0 * z2;
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[10], dN[10][3];
    basis(gpCoord[gp * 3 + 0], gpCoord[gp * 3 + 1], gpCoord[gp * 3 + 2], N, dN);
    eval3DFromBasis<10, true>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                              x_gp_out, detJ_out, w_vol_out);
  }
};

// =========================================================================== //
// 9. H20 (ndm=3, nNodes=20, nGP=8 — 20-node serendipity brick with GP8 rule)
//    Loop order: k (zeta -g..+g) outer, j (eta -g..+g) middle, i (xi -g..+g) inner.
// =========================================================================== //
struct H20 {
  static constexpr int ndm = 3;
  static constexpr int nNodes = 20;
  static constexpr int nGP = 8;
  static constexpr double G = 0.57735026918962576451;
  static constexpr double gpCoord[nGP * ndm] = {
    -G, -G, -G,
    +G, -G, -G,
    -G, +G, -G,
    +G, +G, -G,
    -G, -G, +G,
    +G, -G, +G,
    -G, +G, +G,
    +G, +G, +G
  };
  static constexpr double gpWeight[nGP] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  static constexpr double nodeXi[20][3] = {
    {-1.0, -1.0, -1.0}, {+1.0, -1.0, -1.0}, {+1.0, +1.0, -1.0}, {-1.0, +1.0, -1.0},
    {-1.0, -1.0, +1.0}, {+1.0, -1.0, +1.0}, {+1.0, +1.0, +1.0}, {-1.0, +1.0, +1.0},
    { 0.0, -1.0, -1.0}, {+1.0,  0.0, -1.0}, { 0.0, +1.0, -1.0}, {-1.0,  0.0, -1.0},
    { 0.0, -1.0, +1.0}, {+1.0,  0.0, +1.0}, { 0.0, +1.0, +1.0}, {-1.0,  0.0, +1.0},
    {-1.0, -1.0,  0.0}, {+1.0, -1.0,  0.0}, {+1.0, +1.0,  0.0}, {-1.0, +1.0,  0.0}
  };

  static void basis(double xi, double eta, double zeta, double N[20], double dN[20][3]) {
    for (int a = 0; a < 20; ++a) {
      const double xa = nodeXi[a][0];
      const double ya = nodeXi[a][1];
      const double za = nodeXi[a][2];
      const double px = 1.0 + xi * xa;
      const double py = 1.0 + eta * ya;
      const double pz = 1.0 + zeta * za;

      if (xa != 0.0 && ya != 0.0 && za != 0.0) {
        const double s = xi * xa + eta * ya + zeta * za;
        N[a]     = 0.125 * px * py * pz * (s - 2.0);
        dN[a][0] = 0.125 * xa * py * pz * (2.0 * xi * xa + eta * ya + zeta * za - 1.0);
        dN[a][1] = 0.125 * ya * px * pz * (xi * xa + 2.0 * eta * ya + zeta * za - 1.0);
        dN[a][2] = 0.125 * za * px * py * (xi * xa + eta * ya + 2.0 * zeta * za - 1.0);
      } else if (xa == 0.0) {
        const double b = 1.0 - xi * xi;
        N[a]     =  0.25 * b * py * pz;
        dN[a][0] = -0.50 * xi * py * pz;
        dN[a][1] =  0.25 * b * ya * pz;
        dN[a][2] =  0.25 * b * py * za;
      } else if (ya == 0.0) {
        const double b = 1.0 - eta * eta;
        N[a]     =  0.25 * b * px * pz;
        dN[a][0] =  0.25 * b * xa * pz;
        dN[a][1] = -0.50 * eta * px * pz;
        dN[a][2] =  0.25 * b * px * za;
      } else {
        const double b = 1.0 - zeta * zeta;
        N[a]     =  0.25 * b * px * py;
        dN[a][0] =  0.25 * b * xa * py;
        dN[a][1] =  0.25 * b * px * ya;
        dN[a][2] = -0.50 * zeta * px * py;
      }
    }
  }

  static void eval(int gp, const double* coords, double* N_out, double* dNdx_out,
                   double* detJ_out = nullptr, double* x_gp_out = nullptr,
                   double* w_vol_out = nullptr) {
    assert(gp >= 0 && gp < nGP);
    double N[20], dN[20][3];
    basis(gpCoord[gp * 3 + 0], gpCoord[gp * 3 + 1], gpCoord[gp * 3 + 2], N, dN);
    eval3DFromBasis<20, false>(coords, N, dN, gpWeight[gp], N_out, dNdx_out,
                               x_gp_out, detJ_out, w_vol_out);
  }
};

// ---- Topology metadata & runtime dispatch helpers ------------------------ //
inline int topologyNdm(CPMLTopology topo) {
  switch (topo) {
    case CPMLTopology::T3:
    case CPMLTopology::Q4:
    case CPMLTopology::T6:
    case CPMLTopology::BT6:
      return 2;
    case CPMLTopology::TET4:
    case CPMLTopology::H8:
    case CPMLTopology::TET10:
    case CPMLTopology::BTET10:
    case CPMLTopology::H20:
      return 3;
  }
  return 0;
}

inline int topologyNumNodes(CPMLTopology topo) {
  switch (topo) {
    case CPMLTopology::T3:     return 3;
    case CPMLTopology::Q4:     return 4;
    case CPMLTopology::T6:     return 6;
    case CPMLTopology::BT6:    return 6;
    case CPMLTopology::TET4:   return 4;
    case CPMLTopology::H8:     return 8;
    case CPMLTopology::TET10:  return 10;
    case CPMLTopology::BTET10: return 10;
    case CPMLTopology::H20:    return 20;
  }
  return 0;
}

inline int topologyNumGP(CPMLTopology topo) {
  switch (topo) {
    case CPMLTopology::T3:     return 1;
    case CPMLTopology::Q4:     return 4;
    case CPMLTopology::T6:     return 3;
    case CPMLTopology::BT6:    return 3;
    case CPMLTopology::TET4:   return 1;
    case CPMLTopology::H8:     return 8;
    case CPMLTopology::TET10:  return 4;
    case CPMLTopology::BTET10: return 4;
    case CPMLTopology::H20:    return 8;
  }
  return 0;
}

inline int numVertices(CPMLTopology topo) noexcept {
  switch (topo) {
    case CPMLTopology::T3:     return 3;
    case CPMLTopology::Q4:     return 4;
    case CPMLTopology::T6:     return 3;
    case CPMLTopology::BT6:    return 3;
    case CPMLTopology::TET4:   return 4;
    case CPMLTopology::H8:     return 8;
    case CPMLTopology::TET10:  return 4;
    case CPMLTopology::BTET10: return 4;
    case CPMLTopology::H20:    return 8;
  }
  return 0;
}

inline int topologyNumVertices(CPMLTopology topo) noexcept {
  return numVertices(topo);
}

inline double elemSizeH(CPMLTopology topo, const double* coords) noexcept {
  if (coords == nullptr) return 0.0;
  static constexpr int edgesTri[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  static constexpr int edgesQ4[4][2]  = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static constexpr int edgesTet[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
  static constexpr int edgesHex[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
  };

  const int (*edges)[2] = nullptr;
  int nEdges = 0;
  const int ndm = topologyNdm(topo);
  switch (topo) {
    case CPMLTopology::T3:
    case CPMLTopology::T6:
    case CPMLTopology::BT6:
      edges = edgesTri;
      nEdges = 3;
      break;
    case CPMLTopology::Q4:
      edges = edgesQ4;
      nEdges = 4;
      break;
    case CPMLTopology::TET4:
    case CPMLTopology::TET10:
    case CPMLTopology::BTET10:
      edges = edgesTet;
      nEdges = 6;
      break;
    case CPMLTopology::H8:
    case CPMLTopology::H20:
      edges = edgesHex;
      nEdges = 12;
      break;
  }
  double hmax2 = 0.0;
  for (int e = 0; e < nEdges; ++e) {
    const int a = edges[e][0];
    const int b = edges[e][1];
    double d2 = 0.0;
    for (int i = 0; i < ndm; ++i) {
      const double diff = coords[a * ndm + i] - coords[b * ndm + i];
      d2 += diff * diff;
    }
    if (d2 > hmax2) hmax2 = d2;
  }
  return std::sqrt(hmax2);
}

inline void evalTopologyGP(CPMLTopology topo, int gp, const double* coords,
                           double* N_out, double* dNdx_out,
                           double* detJ_out = nullptr,
                           double* x_gp_out = nullptr,
                           double* w_vol_out = nullptr) {
  switch (topo) {
    case CPMLTopology::T3:
      T3::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::Q4:
      Q4::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::T6:
      T6::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::BT6:
      BT6::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::TET4:
      TET4::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::H8:
      H8::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::TET10:
      TET10::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::BTET10:
      BTET10::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
    case CPMLTopology::H20:
      H20::eval(gp, coords, N_out, dNdx_out, detJ_out, x_gp_out, w_vol_out);
      break;
  }
}

// Evaluate pressure basis Np[0..nP_out-1] and physical gradients
// dNp_dx[b * ndm + i] at parent Gauss point `gp`:
//   - pOrder == 0 (equal-order): nP_out = topologyNumNodes(topo), Np = Nu, dNp_dx = dNu_dx.
//   - pOrder == 1 (linear / Taylor-Hood): nP_out = numVertices(topo), vertex linear
//     reference basis mapped through the parent element's physical Jacobian J_g^{-1}.
inline bool evalShapeAndGradP(CPMLTopology topo, int gp, const double* coords,
                              int pOrder, double* Np, double* dNp_dx,
                              int& nP_out) noexcept {
  if (coords == nullptr || Np == nullptr || dNp_dx == nullptr) return false;
  if (pOrder == 0) {
    nP_out = topologyNumNodes(topo);
    double detJ = 0.0;
    evalTopologyGP(topo, gp, coords, Np, dNp_dx, &detJ, nullptr, nullptr);
    return (detJ != 0.0) && std::isfinite(detJ);
  }
  if (pOrder != 1) {
    nP_out = 0;
    return false;
  }

  nP_out = numVertices(topo);
  switch (topo) {
    case CPMLTopology::T3:
    case CPMLTopology::Q4:
    case CPMLTopology::TET4:
    case CPMLTopology::H8: {
      double detJ = 0.0;
      evalTopologyGP(topo, gp, coords, Np, dNp_dx, &detJ, nullptr, nullptr);
      return (detJ != 0.0) && std::isfinite(detJ);
    }
    case CPMLTopology::T6: {
      assert(gp >= 0 && gp < T6::nGP);
      const double s = T6::gpCoord[gp * 2 + 0];
      const double t = T6::gpCoord[gp * 2 + 1];
      double Nu[6], dNu[6][2];
      T6::basis(s, t, Nu, dNu);
      double J[4] = {0.0, 0.0, 0.0, 0.0};
      for (int a = 0; a < 6; ++a) {
        const double xa = coords[a * 2 + 0];
        const double ya = coords[a * 2 + 1];
        J[0] += dNu[a][0] * xa;
        J[1] += dNu[a][0] * ya;
        J[2] += dNu[a][1] * xa;
        J[3] += dNu[a][1] * ya;
      }
      double Jinv[4];
      const double detJ = invertJacobian2D(J, Jinv);
      Np[0] = s;
      Np[1] = t;
      Np[2] = 1.0 - s - t;
      const double dNp_dxi[3][2] = {
        { 1.0,  0.0},
        { 0.0,  1.0},
        {-1.0, -1.0}
      };
      for (int b = 0; b < 3; ++b) {
        dNp_dx[b * 2 + 0] = Jinv[0] * dNp_dxi[b][0] + Jinv[1] * dNp_dxi[b][1];
        dNp_dx[b * 2 + 1] = Jinv[2] * dNp_dxi[b][0] + Jinv[3] * dNp_dxi[b][1];
      }
      return (detJ != 0.0) && std::isfinite(detJ);
    }
    case CPMLTopology::BT6: {
      assert(gp >= 0 && gp < BT6::nGP);
      const double l1 = BT6::gpCoord[gp * 2 + 0];
      const double l2 = BT6::gpCoord[gp * 2 + 1];
      double Nu[6], dNu[6][2];
      BT6::basis(l1, l2, Nu, dNu);
      double J[4] = {0.0, 0.0, 0.0, 0.0};
      for (int a = 0; a < 6; ++a) {
        const double xa = coords[a * 2 + 0];
        const double ya = coords[a * 2 + 1];
        J[0] += dNu[a][0] * xa;
        J[1] += dNu[a][0] * ya;
        J[2] += dNu[a][1] * xa;
        J[3] += dNu[a][1] * ya;
      }
      double Jinv[4];
      const double detJ = invertJacobian2D(J, Jinv);
      Np[0] = 1.0 - l1 - l2;
      Np[1] = l1;
      Np[2] = l2;
      const double dNp_dxi[3][2] = {
        {-1.0, -1.0},
        { 1.0,  0.0},
        { 0.0,  1.0}
      };
      for (int b = 0; b < 3; ++b) {
        dNp_dx[b * 2 + 0] = Jinv[0] * dNp_dxi[b][0] + Jinv[1] * dNp_dxi[b][1];
        dNp_dx[b * 2 + 1] = Jinv[2] * dNp_dxi[b][0] + Jinv[3] * dNp_dxi[b][1];
      }
      return (detJ != 0.0) && std::isfinite(detJ);
    }
    case CPMLTopology::TET10:
    case CPMLTopology::BTET10: {
      assert(gp >= 0 && gp < TET10::nGP);
      const double z1 = TET10::gpCoord[gp * 3 + 0];
      const double z2 = TET10::gpCoord[gp * 3 + 1];
      const double z3 = TET10::gpCoord[gp * 3 + 2];
      double Nu[10], dNu[10][3];
      if (topo == CPMLTopology::TET10) {
        TET10::basis(z1, z2, z3, Nu, dNu);
      } else {
        BTET10::basis(z1, z2, z3, Nu, dNu);
      }
      double J[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      for (int a = 0; a < 10; ++a) {
        const double xa = coords[a * 3 + 0];
        const double ya = coords[a * 3 + 1];
        const double za = coords[a * 3 + 2];
        for (int p = 0; p < 3; ++p) {
          J[p * 3 + 0] += dNu[a][p] * xa;
          J[p * 3 + 1] += dNu[a][p] * ya;
          J[p * 3 + 2] += dNu[a][p] * za;
        }
      }
      double Jinv[9];
      const double detJ = invertJacobian3D(J, Jinv);
      Np[0] = z1;
      Np[1] = z2;
      Np[2] = z3;
      Np[3] = 1.0 - z1 - z2 - z3;
      const double dNp_dxi[4][3] = {
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0},
        {-1.0, -1.0, -1.0}
      };
      for (int b = 0; b < 4; ++b) {
        for (int i = 0; i < 3; ++i) {
          dNp_dx[b * 3 + i] = Jinv[i * 3 + 0] * dNp_dxi[b][0]
                            + Jinv[i * 3 + 1] * dNp_dxi[b][1]
                            + Jinv[i * 3 + 2] * dNp_dxi[b][2];
        }
      }
      return (detJ != 0.0) && std::isfinite(detJ);
    }
    case CPMLTopology::H20: {
      assert(gp >= 0 && gp < H20::nGP);
      const double xi   = H20::gpCoord[gp * 3 + 0];
      const double eta  = H20::gpCoord[gp * 3 + 1];
      const double zeta = H20::gpCoord[gp * 3 + 2];
      double Nu[20], dNu[20][3];
      H20::basis(xi, eta, zeta, Nu, dNu);
      double J[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      for (int a = 0; a < 20; ++a) {
        const double xa = coords[a * 3 + 0];
        const double ya = coords[a * 3 + 1];
        const double za = coords[a * 3 + 2];
        for (int p = 0; p < 3; ++p) {
          J[p * 3 + 0] += dNu[a][p] * xa;
          J[p * 3 + 1] += dNu[a][p] * ya;
          J[p * 3 + 2] += dNu[a][p] * za;
        }
      }
      double Jinv[9];
      const double detJ = invertJacobian3D(J, Jinv);
      double dNp_dxi[8][3];
      H8::basis(xi, eta, zeta, Np, dNp_dxi);
      for (int b = 0; b < 8; ++b) {
        for (int i = 0; i < 3; ++i) {
          dNp_dx[b * 3 + i] = Jinv[i * 3 + 0] * dNp_dxi[b][0]
                            + Jinv[i * 3 + 1] * dNp_dxi[b][1]
                            + Jinv[i * 3 + 2] * dNp_dxi[b][2];
        }
      }
      return (detJ != 0.0) && std::isfinite(detJ);
    }
  }
  return false;
}

} // namespace shapes

using shapes::numVertices;
using shapes::elemSizeH;
using shapes::evalShapeAndGradP;

} // namespace LadrunoCPML

namespace ladruno_cpml = LadrunoCPML;

#endif // LadrunoCPMLShapes_h
