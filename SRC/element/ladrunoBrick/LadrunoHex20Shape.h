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

// LadrunoHex20Shape.h — pure 20-node serendipity hexahedron (H20) shape-function
// + quadrature kernel (Ladruno fork, ADR 72 P0).
//
// PURE header: plain doubles, ZERO OpenSees includes — unit-testable with a
// standalone g++ driver against the sympy/numpy oracle (tests/hex20_reference.py
// + tests/hex20_kernel_check.cpp), the LadrunoMassLumping.h / LogStrainKernel.h
// idiom. Consumed by LadrunoBrick20 (ELE_TAG 33018, P1) and by any future
// embedded-host / basisInfo work.
//
// ============================ ORDERING MEMO (task 0.1) =======================
// Node ordering — upstream Twenty_Node_Brick / shap3dv "Local Node Pattern"
// (SRC/element/UP-ucsd/shp3dv.cpp:51-67, L/M/N selector tables :75-79; consumed
// by Twenty_Node_Brick.cpp via compuLocalShapeFunction; independently tabulated
// by the fork recorder, SRC/recorder/Ladruno_ElementResults.h ~:1418). This is
// the Abaqus C3D20 / LS-DYNA *ELEMENT_SOLID_H20 (Fig. 19-31) convention:
//   k 0- 3 : corners, LOWER face (zeta=-1), counterclockwise:
//            (-1,-1,-1) (+1,-1,-1) (+1,+1,-1) (-1,+1,-1)
//   k 4- 7 : corners, UPPER face (zeta=+1), same pattern
//   k 8-11 : mid-edges of the LOWER ring: edges 0-1, 1-2, 2-3, 3-0
//            ( 0,-1,-1) (+1, 0,-1) ( 0,+1,-1) (-1, 0,-1)
//   k12-15 : mid-edges of the UPPER ring: edges 4-5, 5-6, 6-7, 7-4
//            ( 0,-1,+1) (+1, 0,+1) ( 0,+1,+1) (-1, 0,+1)
//   k16-19 : VERTICAL mid-edges: edges 0-4, 1-5, 2-6, 3-7
//            (-1,-1, 0) (+1,-1, 0) (+1,+1, 0) (-1,+1, 0)
// (shp3dv is the variable-node 8..27 generator — nodes 21-27 are the face/centroid
// slots of the 27-node variant, absent here; with only the 20 serendipity nodes
// present its hierarchic corrections reduce to the standard serendipity set.
// Exact equivalence to upstream is pinned by the P1 reduce-to gate.)
//
// 27-pt Gauss ordering — the brcshl pattern (shp3dv.cpp brcshl, RA/SA/TA arrays
// :251-277; natural coords = 2*{RA,SA,TA}[L]*g, g = sqrt(3/5); matches the fork
// recorder rule Hexahedron_GaussLegendre_3, Ladruno_ElementResults.h:1266-1289):
//   L 0- 7 : corners  (+-g,+-g,+-g), node-pattern order (lower CCW, upper CCW)
//   L 8-19 : edge mids, node-pattern order (lower ring, upper ring, vertical)
//   L20-25 : face centres +xi, +eta, +zeta, -xi, -eta, -zeta
//   L26    : centroid (0,0,0)
//   weights: tensor products of {5/9, 8/9}: corner (5/9)^3, edge (5/9)^2(8/9),
//            face (5/9)(8/9)^2, centroid (8/9)^3; sum = 8.
//
// 8-pt Gauss ordering — the 8-node Brick rule (recorder Hexahedron_GaussLegendre_2,
// Ladruno_ElementResults.h:1261-1265: "i,j,k lexicographic, z fastest"):
//   (-a,-a,-a) (-a,-a,+a) (-a,+a,-a) (-a,+a,+a)
//   (+a,-a,-a) (+a,-a,+a) (+a,+a,-a) (+a,+a,+a),  a = 1/sqrt(3), w = 1.
// =============================================================================
//
// Voigt convention (matches LadrunoBrick::computeB / OpenSees NDMaterial 3D):
//   strain = {eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx}
//   (ENGINEERING shear: gamma_ij = u_i,j + u_j,i).
//
// HRZ lumping: this kernel provides the consistent mass, the row-sum lump (for
// the P0 negativity demonstration ONLY — row-sum H20 corner masses are -M/8,
// FORBIDDEN in production) and the dofDirs() table; the positive HRZ lump is
// obtained by the caller via the shared Ladruno::hrzLumpRaw
// (SRC/analysis/integrator/LadrunoMassLumping.h):
//     double M[60*60]; int dir[60]; double mdiag[60];
//     Ladruno::hex20::consistentMass(X, rho, Ladruno::hex20::GP27, 27, M);
//     Ladruno::hex20::dofDirs(dir);
//     Ladruno::hrzLumpRaw(M, dir, 60, mdiag);   // positive by construction
//
// Written: N. Mora-Bowen (Ladruno), 2026. ADR 72 P0.  // Ladruno

#ifndef LadrunoHex20Shape_h
#define LadrunoHex20Shape_h

namespace Ladruno {
namespace hex20 {

// ---------------------------------------------------------------- constants --
static const int NEN = 20;   // nodes
static const int NDF = 3;    // dofs per node
static const int NDOF = 60;  // element dofs

// node natural coordinates, brcshl order (see ORDERING MEMO)
static const double NODE_XI[20][3] = {
  {-1.0, -1.0, -1.0}, {+1.0, -1.0, -1.0}, {+1.0, +1.0, -1.0}, {-1.0, +1.0, -1.0},
  {-1.0, -1.0, +1.0}, {+1.0, -1.0, +1.0}, {+1.0, +1.0, +1.0}, {-1.0, +1.0, +1.0},
  { 0.0, -1.0, -1.0}, {+1.0,  0.0, -1.0}, { 0.0, +1.0, -1.0}, {-1.0,  0.0, -1.0},
  { 0.0, -1.0, +1.0}, {+1.0,  0.0, +1.0}, { 0.0, +1.0, +1.0}, {-1.0,  0.0, +1.0},
  {-1.0, -1.0,  0.0}, {+1.0, -1.0,  0.0}, {+1.0, +1.0,  0.0}, {-1.0, +1.0,  0.0}
};

// 27-pt rule (brcshl order): {xi, eta, zeta, w}
static const double GP_G  = 0.774596669241483377;              // sqrt(3/5)
static const double GP_WC = (5.0/9.0)*(5.0/9.0)*(5.0/9.0);     // corner  125/729
static const double GP_WE = (5.0/9.0)*(5.0/9.0)*(8.0/9.0);     // edge    200/729
static const double GP_WF = (5.0/9.0)*(8.0/9.0)*(8.0/9.0);     // face    320/729
static const double GP_W0 = (8.0/9.0)*(8.0/9.0)*(8.0/9.0);     // centre  512/729

static const double GP27[27][4] = {
  {-GP_G, -GP_G, -GP_G, GP_WC}, {+GP_G, -GP_G, -GP_G, GP_WC},   // 0-7 corners
  {+GP_G, +GP_G, -GP_G, GP_WC}, {-GP_G, +GP_G, -GP_G, GP_WC},
  {-GP_G, -GP_G, +GP_G, GP_WC}, {+GP_G, -GP_G, +GP_G, GP_WC},
  {+GP_G, +GP_G, +GP_G, GP_WC}, {-GP_G, +GP_G, +GP_G, GP_WC},
  {  0.0, -GP_G, -GP_G, GP_WE}, {+GP_G,   0.0, -GP_G, GP_WE},   // 8-11 lower edges
  {  0.0, +GP_G, -GP_G, GP_WE}, {-GP_G,   0.0, -GP_G, GP_WE},
  {  0.0, -GP_G, +GP_G, GP_WE}, {+GP_G,   0.0, +GP_G, GP_WE},   // 12-15 upper edges
  {  0.0, +GP_G, +GP_G, GP_WE}, {-GP_G,   0.0, +GP_G, GP_WE},
  {-GP_G, -GP_G,   0.0, GP_WE}, {+GP_G, -GP_G,   0.0, GP_WE},   // 16-19 vertical edges
  {+GP_G, +GP_G,   0.0, GP_WE}, {-GP_G, +GP_G,   0.0, GP_WE},
  {+GP_G,   0.0,   0.0, GP_WF}, {  0.0, +GP_G,   0.0, GP_WF},   // 20-25 faces +x,+y,+z,-x,-y,-z
  {  0.0,   0.0, +GP_G, GP_WF}, {-GP_G,   0.0,   0.0, GP_WF},
  {  0.0, -GP_G,   0.0, GP_WF}, {  0.0,   0.0, -GP_G, GP_WF},
  {  0.0,   0.0,   0.0, GP_W0}                                  // 26 centroid
};

// 8-pt rule (Brick / Hexahedron_GaussLegendre_2 order: lexicographic, zeta fastest)
static const double GP_A = 0.577350269189625765;               // 1/sqrt(3)
static const double GP8[8][4] = {
  {-GP_A, -GP_A, -GP_A, 1.0}, {-GP_A, -GP_A, +GP_A, 1.0},
  {-GP_A, +GP_A, -GP_A, 1.0}, {-GP_A, +GP_A, +GP_A, 1.0},
  {+GP_A, -GP_A, -GP_A, 1.0}, {+GP_A, -GP_A, +GP_A, 1.0},
  {+GP_A, +GP_A, -GP_A, 1.0}, {+GP_A, +GP_A, +GP_A, 1.0}
};

// ------------------------------------------------------------ shape functions --
// Serendipity H20 at natural point xi[3] (ADR 72 §3.1; Felippa AFEM §18.3):
//   corner a:   N = 1/8 (1+xi xi_a)(1+eta eta_a)(1+zeta zeta_a)(xi xi_a + eta eta_a + zeta zeta_a - 2)
//   mid-edge a (xi_a = 0 class): N = 1/4 (1-xi^2)(1+eta eta_a)(1+zeta zeta_a)  (+ cyclic)
// Fills N[20] and dN[20][3] = dN_a/d(xi,eta,zeta).
inline void shape(const double xi[3], double N[20], double dN[20][3])
{
  for (int a = 0; a < 20; ++a) {
    const double xa = NODE_XI[a][0], ya = NODE_XI[a][1], za = NODE_XI[a][2];
    const double px = 1.0 + xi[0] * xa;
    const double py = 1.0 + xi[1] * ya;
    const double pz = 1.0 + xi[2] * za;

    if (xa != 0.0 && ya != 0.0 && za != 0.0) {
      // corner node
      const double s = xi[0] * xa + xi[1] * ya + xi[2] * za;
      N[a]     = 0.125 * px * py * pz * (s - 2.0);
      dN[a][0] = 0.125 * xa * py * pz * (2.0 * xi[0] * xa + xi[1] * ya + xi[2] * za - 1.0);
      dN[a][1] = 0.125 * ya * px * pz * (xi[0] * xa + 2.0 * xi[1] * ya + xi[2] * za - 1.0);
      dN[a][2] = 0.125 * za * px * py * (xi[0] * xa + xi[1] * ya + 2.0 * xi[2] * za - 1.0);
    }
    else if (xa == 0.0) {
      // mid-edge, xi-direction bubble
      const double b = 1.0 - xi[0] * xi[0];
      N[a]     =  0.25 * b * py * pz;
      dN[a][0] = -0.50 * xi[0] * py * pz;
      dN[a][1] =  0.25 * b * ya * pz;
      dN[a][2] =  0.25 * b * py * za;
    }
    else if (ya == 0.0) {
      // mid-edge, eta-direction bubble
      const double b = 1.0 - xi[1] * xi[1];
      N[a]     =  0.25 * b * px * pz;
      dN[a][0] =  0.25 * b * xa * pz;
      dN[a][1] = -0.50 * xi[1] * px * pz;
      dN[a][2] =  0.25 * b * px * za;
    }
    else {
      // mid-edge, zeta-direction bubble (za == 0)
      const double b = 1.0 - xi[2] * xi[2];
      N[a]     =  0.25 * b * px * py;
      dN[a][0] =  0.25 * b * xa * py;
      dN[a][1] =  0.25 * b * px * ya;
      dN[a][2] = -0.50 * xi[2] * px * py;
    }
  }
}

// ------------------------------------------------------------------ geometry --
// Jacobian J[i][j] = d x_i / d xi_j = sum_a X[a][i] dN[a][j]; returns det J.
inline double jacobian(const double X[20][3], const double dN[20][3], double J[3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int a = 0; a < 20; ++a) s += X[a][i] * dN[a][j];
      J[i][j] = s;
    }
  return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
       - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
       + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

// inverse of a 3x3 given its determinant (adjugate / detJ)
inline void invert3(const double J[3][3], double detJ, double Jinv[3][3])
{
  const double d = 1.0 / detJ;
  Jinv[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * d;
  Jinv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * d;
  Jinv[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * d;
  Jinv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * d;
  Jinv[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * d;
  Jinv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * d;
  Jinv[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * d;
  Jinv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * d;
  Jinv[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * d;
}

// cartesian gradients dNdx[a][i] = dN_a/dx_i = sum_j dN[a][j] Jinv[j][i]
inline void cartGrad(const double dN[20][3], const double Jinv[3][3], double dNdx[20][3])
{
  for (int a = 0; a < 20; ++a)
    for (int i = 0; i < 3; ++i)
      dNdx[a][i] = dN[a][0] * Jinv[0][i] + dN[a][1] * Jinv[1][i] + dN[a][2] * Jinv[2][i];
}

// ------------------------------------------------------------------------ B --
// 6x60 small-strain B, Voigt {xx,yy,zz,xy,yz,zx} with engineering shear
// (LadrunoBrick::computeB per-node layout, tiled over 20 nodes).
inline void fillB(const double dNdx[20][3], double B[6][60])
{
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 60; ++j) B[i][j] = 0.0;
  for (int a = 0; a < 20; ++a) {
    const int c = 3 * a;
    const double bx = dNdx[a][0], by = dNdx[a][1], bz = dNdx[a][2];
    B[0][c    ] = bx;
    B[1][c + 1] = by;
    B[2][c + 2] = bz;
    B[3][c    ] = by;  B[3][c + 1] = bx;
    B[4][c + 1] = bz;  B[4][c + 2] = by;
    B[5][c    ] = bz;  B[5][c + 2] = bx;
  }
}

// -------------------------------------------------------------------- volume --
inline double volume(const double X[20][3], const double gp[][4], int ngp)
{
  double N[20], dN[20][3], J[3][3];
  double V = 0.0;
  for (int L = 0; L < ngp; ++L) {
    const double xi[3] = { gp[L][0], gp[L][1], gp[L][2] };
    shape(xi, N, dN);
    V += gp[L][3] * jacobian(X, dN, J);
  }
  return V;
}

// ---------------------------------------------------------------------- mass --
// consistent mass M[60*60] row-major: M(3a+d, 3b+d) = int rho N_a N_b dV.
// ALWAYS integrate mass with the FULL rule (GP27) regardless of the stiffness
// formulation (ADR 72 §2.3 / Abaqus AUG §28.1.1 convention).
inline void consistentMass(const double X[20][3], double rho,
                           const double gp[][4], int ngp, double *M /*60*60*/)
{
  for (int k = 0; k < 60 * 60; ++k) M[k] = 0.0;
  double N[20], dN[20][3], J[3][3];
  for (int L = 0; L < ngp; ++L) {
    const double xi[3] = { gp[L][0], gp[L][1], gp[L][2] };
    shape(xi, N, dN);
    const double dv = gp[L][3] * jacobian(X, dN, J) * rho;
    for (int a = 0; a < 20; ++a)
      for (int b = 0; b < 20; ++b) {
        const double m = dv * N[a] * N[b];
        for (int d = 0; d < 3; ++d)
          M[(3 * a + d) * 60 + (3 * b + d)] += m;
      }
  }
}

// row-sum lump (P0 negativity DEMONSTRATION ONLY — H20 corner row-sums are
// -M/8: never use in production; use hrzLumpRaw instead, see header comment).
inline void rowSumLump(const double *M /*60*60*/, double mdiag[60])
{
  for (int i = 0; i < 60; ++i) {
    double s = 0.0;
    for (int j = 0; j < 60; ++j) s += M[i * 60 + j];
    mdiag[i] = s;
  }
}

// dof direction table for Ladruno::hrzLumpRaw (all-translational: 0,1,2 repeating)
inline void dofDirs(int dir[60])
{
  for (int a = 0; a < 20; ++a) {
    dir[3 * a] = 0; dir[3 * a + 1] = 1; dir[3 * a + 2] = 2;
  }
}

// ---------------------------------------------------- elastic stiffness (P0) --
// K[60*60] row-major += B^T D B detJ w over the given rule, constant 6x6 D.
// P0 oracle/driver helper (rank + cross-check); the element assembles with the
// per-GP material tangent in its own loop.
inline void stiffnessElastic(const double X[20][3], const double D[6][6],
                             const double gp[][4], int ngp, double *K /*60*60*/)
{
  for (int k = 0; k < 60 * 60; ++k) K[k] = 0.0;
  double N[20], dN[20][3], J[3][3], Jinv[3][3], dNdx[20][3], B[6][60], DB[6][60];
  for (int L = 0; L < ngp; ++L) {
    const double xi[3] = { gp[L][0], gp[L][1], gp[L][2] };
    shape(xi, N, dN);
    const double detJ = jacobian(X, dN, J);
    invert3(J, detJ, Jinv);
    cartGrad(dN, Jinv, dNdx);
    fillB(dNdx, B);
    const double dv = gp[L][3] * detJ;
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 60; ++j) {
        double s = 0.0;
        for (int m = 0; m < 6; ++m) s += D[i][m] * B[m][j];
        DB[i][j] = s;
      }
    for (int i = 0; i < 60; ++i)
      for (int j = 0; j < 60; ++j) {
        double s = 0.0;
        for (int m = 0; m < 6; ++m) s += B[m][i] * DB[m][j];
        K[i * 60 + j] += dv * s;
      }
  }
}

} // namespace hex20
} // namespace Ladruno

#endif // LadrunoHex20Shape_h
