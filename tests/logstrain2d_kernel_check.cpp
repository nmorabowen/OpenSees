// Standalone driver for the 2D (plane) log-strain algorithm — NO OpenSees deps.
// Reproduces the math of LogStrain2D::setTrialF / buildOutputs using ONLY the
// pure kernel (SRC/material/nD/LogStrainKernel.h) plus an inline isotropic elastic
// inner law: plane strain lifts F to 3×3 with F₃₃ = 1 and reads back the in-plane
// block; plane stress solves F₃₃ from σ₃₃ = 0 (local Newton) and statically
// condenses the out-of-plane direction. The companion pytest
// (test_logstrain2d_kernel_cpp.py) diffs σ (3) and the full in-plane modulus c2
// (16) against the numpy oracle (logstrain2d_reference.py) — catching any C++
// transcription bug without building OpenSees.

#include <cstdio>
#include <cmath>
#include <LogStrainKernel.h>

using namespace logstrain_kernel;

static const double K = 1500.0;
static const double G = 700.0;

// inline isotropic elastic small-strain law (engineering-Voigt), τ = D:ε
static void elastic(const double eps6[6], double tau6[6], double D6[36]) {
  double tr = eps6[0] + eps6[1] + eps6[2];
  for (int i = 0; i < 3; i++) tau6[i] = K*tr + 2.0*G*(eps6[i] - tr/3.0);
  for (int i = 3; i < 6; i++) tau6[i] = G*eps6[i];
  for (int k = 0; k < 36; k++) D6[k] = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      D6[6*i + j] = K + 2.0*G*((i==j ? 1.0 : 0.0) - 1.0/3.0);
  D6[6*3+3] = G; D6[6*4+4] = G; D6[6*5+5] = G;
}

// One lifted evaluation: F2 (row-major 2×2) + F33 → Cauchy σ6, J3, full c3d[3][3][3][3].
static void lifted(const double F2[4], double F33, double sig6[6],
                   double c3d[3][3][3][3]) {
  double I9[9] = {1,0,0, 0,1,0, 0,0,1};
  double F3[9] = {F2[0], F2[1], 0.0,  F2[2], F2[3], 0.0,  0.0, 0.0, F33};
  double Betr[9]; trial_Be(F3, I9, I9, Betr);
  double eps6[6]; hencky_voigt(Betr, eps6);
  double tau6[6], D6[36]; elastic(eps6, tau6, D6);
  double J3 = mat3_det(F3);
  for (int k = 0; k < 6; k++) sig6[k] = tau6[k] / J3;
  spatial_tangent_full(Betr, D6, J3, c3d);
}

// build the in-plane material modulus c2 exactly as LogStrain2D::buildOutputs
static void build_c2(const double sig6[6], const double c3d[3][3][3][3],
                     bool condense, double c2[2][2][2][2]) {
  if (!condense) {
    for (int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<2;k++) for(int l=0;l<2;l++)
      c2[i][j][k][l] = c3d[i][j][k][l];
    return;
  }
  double sig[3][3];
  sig[0][0]=sig6[0]; sig[1][1]=sig6[1]; sig[2][2]=sig6[2];
  sig[0][1]=sig[1][0]=sig6[3];
  sig[1][2]=sig[2][1]=sig6[4];
  sig[2][0]=sig[0][2]=sig6[5];
  double a3[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    a3[i][j][k][l] = c3d[i][j][k][l] - sig[i][l]*(j==k ? 1.0 : 0.0);
  double denom = a3[2][2][2][2];
  for (int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<2;k++) for(int l=0;l<2;l++) {
    double a2 = a3[i][j][k][l] - a3[i][j][2][2]*a3[2][2][k][l]/denom;
    c2[i][j][k][l] = a2 + sig[i][l]*(j==k ? 1.0 : 0.0);
  }
}

static void run(const char *name, int planeStress, const double F2[4]) {
  double sig6[6], c3d[3][3][3][3];
  double F33 = 1.0;
  if (planeStress) {                          // Newton: solve σ₃₃(λ) = 0
    double lam = 1.0;
    for (int it = 0; it < 30; it++) {
      lifted(F2, lam, sig6, c3d);
      double r = sig6[2];
      double sref = 1.0;
      if (fabs(sig6[0]) > sref) sref = fabs(sig6[0]);
      if (fabs(sig6[1]) > sref) sref = fabs(sig6[1]);
      if (fabs(r) <= 1e-12*sref) break;
      double h = 1e-8*((fabs(lam)>1.0)?fabs(lam):1.0);
      double s2[6], cc[3][3][3][3];
      lifted(F2, lam + h, s2, cc);
      double drdl = (s2[2] - r)/h;
      double lamNew = lam - r/drdl;
      if (lamNew <= 0.0) lamNew = 0.5*lam;
      lam = lamNew;
    }
    lifted(F2, lam, sig6, c3d);
    F33 = lam;
  } else {
    lifted(F2, 1.0, sig6, c3d);
  }
  double c2[2][2][2][2];
  build_c2(sig6, c3d, planeStress != 0, c2);

  printf("CASE %s\n", name);
  printf("SIGMA %.15g %.15g %.15g\n", sig6[0], sig6[1], sig6[3]);   // {σ11,σ22,σ12}
  printf("C2");
  for (int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<2;k++) for(int l=0;l<2;l++)
    printf(" %.15g", c2[i][j][k][l]);
  printf("\n");
}

int main() {
  // in-plane 2×2 F cases (row-major); names must match the pytest.
  double generic[4] = {1.15, 0.08, -0.04, 0.93};
  double shear[4]   = {1.0, 0.20, 0.0, 1.0};
  double equibi[4]  = {1.10, 0.0, 0.0, 1.10};   // 2-fold in-plane (kernel degeneracy)
  double compress[4]= {0.90, 0.05, 0.03, 0.88};

  run("ps_generic",  0, generic);   // plane strain
  run("ps_shear",    0, shear);
  run("ps_equibi",   0, equibi);
  run("pss_generic", 1, generic);   // plane stress
  run("pss_shear",   1, shear);
  run("pss_equibi",  1, equibi);
  run("pss_compress",1, compress);
  return 0;
}
