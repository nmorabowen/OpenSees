// Standalone driver for the shared 2D finite-strain assembly kernel — NO OpenSees
// deps. Reimplements a single-element updated-Lagrangian assembly using ONLY the
// pure kernels (SRC/element/ladrunoPlane/LadrunoFiniteStrain2DKernel.h +
// SRC/material/nD/LogStrainKernel.h for an inline elastic-Hencky plane-strain
// material) and prints f_int and K for a Q4 (std + F-bar) and a T3. The companion
// pytest (test_finitestrain2d_kernel_cpp.py) diffs this against the numpy oracle
// (finitestrain2d_reference.py) — catching any C++ transcription bug in the finite-
// strain element math (F-bar G₀ coupling included) without building OpenSees.

#include <cstdio>
#include <cmath>
#include <LogStrainKernel.h>
#include <LadrunoFiniteStrain2DKernel.h>

using namespace ladruno_fs2d;

static const double Kmod = 1500.0, Gmod = 700.0;

// inline isotropic elastic small-strain law (engineering-Voigt), τ = D:ε
static void elastic(const double eps6[6], double tau6[6], double D6[36]) {
  double tr = eps6[0] + eps6[1] + eps6[2];
  for (int i = 0; i < 3; i++) tau6[i] = Kmod * tr + 2.0 * Gmod * (eps6[i] - tr / 3.0);
  for (int i = 3; i < 6; i++) tau6[i] = Gmod * eps6[i];
  for (int k = 0; k < 36; k++) D6[k] = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      D6[6 * i + j] = Kmod + 2.0 * Gmod * ((i == j ? 1.0 : 0.0) - 1.0 / 3.0);
  D6[6 * 3 + 3] = Gmod; D6[6 * 4 + 4] = Gmod; D6[6 * 5 + 5] = Gmod;
}

// plane-strain elastic-Hencky material: F̄(2×2) → in-plane Cauchy σ + full c2 block
static void material2D(const double F2[4], double sig[2][2], double c2[2][2][2][2]) {
  using namespace logstrain_kernel;
  double I9[9] = {1,0,0, 0,1,0, 0,0,1};
  double F3[9] = {F2[0], F2[1], 0.0,  F2[2], F2[3], 0.0,  0.0, 0.0, 1.0};  // F33=1
  double Betr[9]; trial_Be(F3, I9, I9, Betr);
  double eps6[6]; hencky_voigt(Betr, eps6);
  double tau6[6], D6[36]; elastic(eps6, tau6, D6);
  double J3 = mat3_det(F3);
  double sig6[6]; for (int k = 0; k < 6; k++) sig6[k] = tau6[k] / J3;
  double c3[3][3][3][3]; spatial_tangent_full(Betr, D6, J3, c3);
  sig[0][0] = sig6[0]; sig[1][1] = sig6[1]; sig[0][1] = sig[1][0] = sig6[3];
  for (int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<2;k++) for(int l=0;l<2;l++)
    c2[i][j][k][l] = c3[i][j][k][l];
}

// ---- shape-function reference gradients (must match finitestrain2d_reference) --- //
// Fill gradRef[i*nN + a] = ∂N_a/∂X_i and return detJ0, given natural derivs dNdxi.
static double refGrads(const double *dNdxi /*2*nN*/, const double *X /*nN*2*/,
                       int nN, double *gradRef) {
  double J[4] = {0,0,0,0};                       // J[a][i] = Σ_n dN_n/dξ_α X_n,i
  for (int n = 0; n < nN; n++) {
    double dx = dNdxi[0 * nN + n], de = dNdxi[1 * nN + n];
    J[0] += dx * X[2 * n + 0];  J[1] += dx * X[2 * n + 1];
    J[2] += de * X[2 * n + 0];  J[3] += de * X[2 * n + 1];
  }
  double det = J[0] * J[3] - J[1] * J[2];
  double Ji[4] = { J[3] / det, -J[1] / det, -J[2] / det, J[0] / det };  // inv(J), row-major
  // ∂N_a/∂X_i = Σ_α Jinv[i][α] dNdξ[α][a]  (Jinv[i][α] = ∂ξ_α/∂X_i)
  for (int a = 0; a < nN; a++) {
    double dx = dNdxi[0 * nN + a], de = dNdxi[1 * nN + a];
    gradRef[0 * nN + a] = Ji[0] * dx + Ji[1] * de;   // i=0: Jinv[0][0], Jinv[0][1]
    gradRef[1 * nN + a] = Ji[2] * dx + Ji[3] * de;   // i=1: Jinv[1][0], Jinv[1][1]
  }
  return det;
}

static void q4_dN(double xi, double eta, double *dNdxi /*2*4*/) {
  dNdxi[0*4+0] = -0.25*(1-eta); dNdxi[0*4+1] = 0.25*(1-eta);
  dNdxi[0*4+2] = 0.25*(1+eta);  dNdxi[0*4+3] = -0.25*(1+eta);
  dNdxi[1*4+0] = -0.25*(1-xi);  dNdxi[1*4+1] = -0.25*(1+xi);
  dNdxi[1*4+2] = 0.25*(1+xi);   dNdxi[1*4+3] = 0.25*(1-xi);
}
static void t3_dN(double, double, double *dNdxi /*2*3*/) {
  dNdxi[0*3+0] = -1.0; dNdxi[0*3+1] = 1.0; dNdxi[0*3+2] = 0.0;
  dNdxi[1*3+0] = -1.0; dNdxi[1*3+1] = 0.0; dNdxi[1*3+2] = 1.0;
}
static void t6_dN(double xi, double eta, double *dNdxi /*2*6*/) {
  double L0 = 1.0 - xi - eta, L1 = xi, L2 = eta;
  // ∂/∂ξ  (dL0=-1, dL1=1, dL2=0)
  dNdxi[0*6+0] = 1 - 4*L0;   dNdxi[0*6+1] = 4*L1 - 1;   dNdxi[0*6+2] = 0.0;
  dNdxi[0*6+3] = 4*(L0-L1);  dNdxi[0*6+4] = 4*L2;       dNdxi[0*6+5] = -4*L2;
  // ∂/∂η  (dL0=-1, dL1=0, dL2=1)
  dNdxi[1*6+0] = 1 - 4*L0;   dNdxi[1*6+1] = 0.0;        dNdxi[1*6+2] = 4*L2 - 1;
  dNdxi[1*6+3] = -4*L1;      dNdxi[1*6+4] = 4*L1;       dNdxi[1*6+5] = 4*(L0-L2);
}
static void shapeDN(int elem, double xi, double eta, double *dN) {
  if (elem == 4) q4_dN(xi, eta, dN);
  else if (elem == 6) t6_dN(xi, eta, dN);
  else t3_dN(xi, eta, dN);
}

struct GP { double xi, eta, w; };

static void assembleElement(const char *name, int elem /*4=Q4,3=T3*/,
                            const double *X, const double *u, int nN,
                            const GP *gps, int ngp, double cx, double ceta, bool fbar) {
  int nd = 2 * nN, ld = nd;
  double resid[16] = {0}, K[16 * 16] = {0};

  double g0[2 * 8];  // centroid spatial gradients (F-bar)
  if (fbar) {
    double dN0[2 * 8], gr0[2 * 8];
    shapeDN(elem, cx, ceta, dN0);
    refGrads(dN0, X, nN, gr0);
    double F0[4]; deformationGradient2D(gr0, u, nN, F0);
    double F0i[4]; inv2(F0, F0i);
    spatialGradients2D(gr0, F0i, nN, g0);
  }

  for (int q = 0; q < ngp; q++) {
    double dN[2 * 8], gradRef[2 * 8];
    shapeDN(elem, gps[q].xi, gps[q].eta, dN);
    double detJ0 = refGrads(dN, X, nN, gradRef);
    double F[4]; double J = deformationGradient2D(gradRef, u, nN, F);
    double Fmat[4] = { F[0], F[1], F[2], F[3] };
    if (fbar) {                                       // F-bar scale uses centroid J0
      double dN0[2 * 8], gr0[2 * 8];
      shapeDN(elem, cx, ceta, dN0);
      refGrads(dN0, X, nN, gr0);
      double F0[4]; double J0 = deformationGradient2D(gr0, u, nN, F0);
      double s = fbarScale2D(J0, J);
      for (int m = 0; m < 4; m++) Fmat[m] = s * F[m];
    }
    double sig[2][2], c2[2][2][2][2];
    material2D(Fmat, sig, c2);

    double Fi[4]; inv2(F, Fi);
    double g[2 * 8]; spatialGradients2D(gradRef, Fi, nN, g);
    double dv = J * detJ0 * gps[q].w;                 // thickness = 1

    addInternalForce2D(sig, g, nN, dv, resid);
    double a4[2][2][2][2]; spatialTangent2D(c2, sig, a4);
    addTangent2D(a4, g, nN, dv, K, ld);
    if (fbar) addFbarCoupling2D(a4, sig, g, g0, nN, dv, K, ld);
  }

  printf("CASE %s\n", name);
  printf("F");
  for (int i = 0; i < nd; i++) printf(" %.15g", resid[i]);
  printf("\nK");
  for (int i = 0; i < nd; i++) for (int j = 0; j < nd; j++) printf(" %.15g", K[i * ld + j]);
  printf("\n");
}

int main() {
  // ref coords + disp must match finitestrain2d_reference.py / the pytest.
  double Xq4[8] = {0.0,0.0, 2.0,0.1, 2.2,1.8, 0.1,1.9};
  double Xt3[6] = {0.0,0.0, 2.0,0.2, 0.3,1.6};
  double Xt6[12] = {0.0,0.0, 2.0,0.0, 0.0,2.0, 1.0,-0.05, 1.05,1.0, -0.05,1.0};
  // deterministic small displacement (matches the pytest's fixed vector)
  double uq4[8] = { 0.010,-0.006, 0.004, 0.012,-0.008, 0.003, 0.007,-0.011};
  double ut3[6] = { 0.009,-0.005, 0.006, 0.008,-0.004, 0.010};
  double ut6[12] = { 0.008,-0.005, 0.006,0.009, -0.004,0.007,
                     0.003,-0.006, 0.005,0.004, -0.007,0.008};

  GP q4gp[4]; double g = 1.0 / sqrt(3.0);
  q4gp[0] = {-g,-g,1}; q4gp[1] = {g,-g,1}; q4gp[2] = {g,g,1}; q4gp[3] = {-g,g,1};
  GP t3gp[1] = {{1.0/3, 1.0/3, 0.5}};
  GP t6gp[3] = {{2.0/3, 1.0/6, 1.0/6}, {1.0/6, 2.0/3, 1.0/6}, {1.0/6, 1.0/6, 1.0/6}};

  assembleElement("Q4_std",  4, Xq4, uq4, 4, q4gp, 4, 0.0, 0.0, false);
  assembleElement("Q4_fbar", 4, Xq4, uq4, 4, q4gp, 4, 0.0, 0.0, true);
  assembleElement("T3_std",  3, Xt3, ut3, 3, t3gp, 1, 1.0/3, 1.0/3, false);
  assembleElement("T6_std",  6, Xt6, ut6, 6, t6gp, 3, 1.0/3, 1.0/3, false);
  assembleElement("T6_fbar", 6, Xt6, ut6, 6, t6gp, 3, 1.0/3, 1.0/3, true);
  return 0;
}
