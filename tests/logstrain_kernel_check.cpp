// Standalone driver for the pure LogStrain numerical kernel — NO OpenSees deps.
// Compiles with plain g++ against SRC/material/nD/LogStrainKernel.h and prints
// the Cauchy stress σ and material spatial tangent c for a set of deformation
// gradients, using an inline isotropic elastic Hencky inner law. The companion
// pytest (test_logstrain_kernel_cpp.py) diffs this output against the numpy
// oracle (logstrain_reference.py) — catching any C++ transcription bug without
// building OpenSees.

#include <cstdio>
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

static void run(const char *name, const double F[9]) {
  double I9[9] = {1,0,0, 0,1,0, 0,0,1};
  double Betr[9]; trial_Be(F, I9, I9, Betr);
  double eps6[6]; hencky_voigt(Betr, eps6);
  double tau6[6], D6[36]; elastic(eps6, tau6, D6);
  double J = mat3_det(F);
  double sig6[6], c6[36];
  assemble_material(Betr, D6, tau6, J, sig6, c6);

  printf("CASE %s\n", name);
  printf("SIGMA");
  for (int k = 0; k < 6; k++) printf(" %.15g", sig6[k]);
  printf("\nC");
  for (int k = 0; k < 36; k++) printf(" %.15g", c6[k]);
  printf("\n");
}

int main() {
  double generic[9]  = {1.20, 0.10, -0.05,  0.02, 0.90, 0.08,  0.03, -0.04, 1.10};
  double uniaxial[9] = {1.25, 0.0, 0.0,  0.0, 0.95, 0.0,  0.0, 0.0, 0.90};
  double rot[9];                       // rotation about z by 0.5 rad (det=1)
  {
    double cz = 0.87758256189037271612, sz = 0.47942553860420300027;
    double R[9] = {cz, -sz, 0.0,  sz, cz, 0.0,  0.0, 0.0, 1.0};
    for (int i = 0; i < 9; i++) rot[i] = R[i];
  }
  double dbl[9] = {1.2, 0.0, 0.0,  0.0, 1.2, 0.0,  0.0, 0.0, 0.8}; // 2-fold b
  double tri[9] = {1.1, 0.0, 0.0,  0.0, 1.1, 0.0,  0.0, 0.0, 1.1}; // 3-fold b

  run("generic",  generic);
  run("uniaxial", uniaxial);
  run("rotation", rot);
  run("double",   dbl);
  run("triple",   tri);
  return 0;
}
