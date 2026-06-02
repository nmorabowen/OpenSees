// Standalone driver for the extracted LadrunoJ2 return-map kernel.
//
// Includes ONLY SRC/material/nD/LadrunoJ2Kernel.h (pure, OpenSees-free), runs a
// set of prescribed strain PATHS threading the committed history, and prints the
// per-step Cauchy stress, internal state, and algorithmic tangent so the numpy
// oracle (tests/ladrunoj2_reference.py, driven from test_ladrunoJ2_kernel_cpp.py)
// can cross-check the C++ transcription WITHOUT building OpenSees.
//
// Build:  g++ -O2 -std=c++17 -I SRC/material/nD ladrunoj2_kernel_check.cpp -o lj2k

#include <LadrunoJ2Kernel.h>
#include <cstdio>

using namespace ladruno_j2_kernel;

// A symmetric strain direction (6 TENSOR components {00,11,22,01,12,02}) scaled
// by a per-step multiplier — the same construction the Python side mirrors.
struct Scenario {
  const char* name;
  Params p;
  double dir[6];
  const double* mult;
  int nstep;
};

static void run(const Scenario& sc)
{
  printf("SCENARIO %s\n", sc.name);

  double epsP_n[6] = {0,0,0,0,0,0};
  double ebarP_n = 0.0;
  double alpha_n[MAXBACK][6];
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = 0.0;

  for (int s = 0; s < sc.nstep; s++) {
    double eps[6];
    for (int i = 0; i < 6; i++) eps[i] = sc.mult[s] * sc.dir[i];

    double stress[6], Dtan[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG;
    returnMap(sc.p, eps, epsP_n, ebarP_n, alpha_n,
              stress, Dtan, epsP, ebarP, alpha, dG);

    printf("STEP %d STRESS", s);
    for (int i = 0; i < 6; i++) printf(" %.17g", stress[i]);
    printf("\n");

    printf("STEP %d STATE %.17g %.17g", s, ebarP, dG);
    for (int i = 0; i < 6; i++) printf(" %.17g", epsP[i]);
    for (int k = 0; k < sc.p.nBack; k++)
      for (int i = 0; i < 6; i++) printf(" %.17g", alpha[k][i]);
    printf("\n");

    printf("STEP %d DTAN", s);
    for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) printf(" %.17g", Dtan[I][J]);
    printf("\n");

    // commit
    for (int i = 0; i < 6; i++) epsP_n[i] = epsP[i];
    ebarP_n = ebarP;
    for (int k = 0; k < sc.p.nBack; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha[k][i];
  }
}

// Multiplier sequences (mirrored exactly in Python).
static const double MONO[6] = {0.002, 0.005, 0.01, 0.02, 0.03, 0.05};
static const double CYC[5]  = {0.03, 0.0, -0.03, 0.0, 0.03};

int main()
{
  // K, G chosen as a typical metal (E≈1820, ν≈0.3 here is irrelevant — K,G direct)
  const double K = 1500.0, G = 700.0;

  Scenario scen[4];

  // 1) isotropic-only Voce+linear, uniaxial monotonic (elastic → plastic)
  scen[0] = {"elastic_iso",
             Params{K, G, 5.0, 4.0, 30.0, 10.0, 0, {0}, {0}},
             {1,0,0, 0,0,0}, MONO, 6};

  // 2) single Armstrong-Frederick, uniaxial cyclic (Bauschinger / saturation)
  scen[1] = {"single_af_cyclic",
             Params{K, G, 5.0, 0.0, 0.0, 0.0, 1, {400.0}, {80.0}},
             {1,0,0, 0,0,0}, CYC, 5};

  // 3) three-term Chaboche + Voce+linear, 3D mixed (shear-inclusive) monotonic
  scen[2] = {"chaboche3_3d",
             Params{K, G, 8.0, 3.0, 20.0, 5.0, 3, {500.0,300.0,150.0}, {100.0,50.0,10.0}},
             {0.5,-0.2,-0.3, 0.4, 0.1, 0.2}, MONO, 6};

  // 4) linear kinematic (Prager, gamma=0), uniaxial cyclic
  scen[3] = {"kin_linear_cyclic",
             Params{K, G, 5.0, 0.0, 0.0, 0.0, 1, {120.0}, {0.0}},
             {1,0,0, 0,0,0}, CYC, 5};

  for (int s = 0; s < 4; s++) run(scen[s]);
  return 0;
}
