// hypo_kernel_check — standalone g++ cross-check of LadrunoHypoKernel.h
// (ADR 79 P0). No OpenSees dependency: includes ONLY the kernel header.
//
//   g++ -O2 -std=c++17 -Wall -Wextra -I SRC/element/solidTransformation
//       tests/hypo_kernel_check.cpp -o hypo_check ; ./hypo_check tests/hypo_cases.txt
//
// (a) Always runs standalone self-tests — the kernel's contract invariants:
//     rigid-rotation increment => deps EXACTLY ~0 at ANY angle (the HW midpoint
//     property), polar orthogonality, the pure-stretch midpoint closed form
//     2(l-1)/(l+1), bond6(I) = identity, tangent-push symmetry preservation.
// (b) If the oracle-emitted cases file (tests/hypo_cases.txt, from
//     tests/hypo_reference.py — SVD polar, independent algorithm) is given and
//     present, diffs deps/R1/J1/sout/Dout to <= 1e-9 per case.
//
// Exit 0 + "RESULT: PASS" on success.

#include <LadrunoHypoKernel.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace ladruno_hypo;

static int g_fail = 0;

static void check(bool ok, const char *what)
{
  if (!ok) {
    std::printf("FAIL: %s\n", what);
    g_fail++;
  }
}

static double maxAbsDiff(const double *a, const double *b, int n)
{
  double m = 0.0;
  for (int i = 0; i < n; i++) {
    double d = std::fabs(a[i] - b[i]);
    if (d > m) m = d;
  }
  return m;
}

// Rodrigues rotation about a (unit) axis by angle th, row-major.
static void rodrigues(const double axis[3], double th, double Q[9])
{
  double c = std::cos(th), s = std::sin(th), v = 1.0 - c;
  double x = axis[0], y = axis[1], z = axis[2];
  Q[0] = c + x*x*v;   Q[1] = x*y*v - z*s; Q[2] = x*z*v + y*s;
  Q[3] = y*x*v + z*s; Q[4] = c + y*y*v;   Q[5] = y*z*v - x*s;
  Q[6] = z*x*v - y*s; Q[7] = z*y*v + x*s; Q[8] = c + z*z*v;
}

// The frozen H8 reference geometry of the oracle: corner cube [-1,1]^3,
// centroid gradients dN_a/dX = corner_a / 8 (linear-complete).
static const double CORNERS[8][3] = {
  {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
  {-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1},
};

static void centroidGrads(double dNdX[24])
{
  for (int a = 0; a < 8; a++)
    for (int j = 0; j < 3; j++)
      dNdX[a*3+j] = CORNERS[a][j] / 8.0;
}

// nodal displacements of the affine map x = A X : u_a = (A - I) X_a
static void affineDisp(const double A[9], double u[24])
{
  for (int a = 0; a < 8; a++)
    for (int i = 0; i < 3; i++) {
      double s = 0.0;
      for (int J = 0; J < 3; J++) s += A[3*i+J] * CORNERS[a][J];
      u[a*3+i] = s - CORNERS[a][i];
    }
}

// ---------------------------------------------------------------------------
// (a) standalone self-tests
// ---------------------------------------------------------------------------
static void selfTests()
{
  double dNdX[24];
  centroidGrads(dNdX);

  // -- 1. rigid-rotation increment: deps == 0 at ANY per-step angle -----------
  {
    const double ax1[3] = {0.267261241912424, 0.534522483824849, 0.801783725737273};
    const double ax2[3] = {0.577350269189626, -0.577350269189626, 0.577350269189626};
    double Q1[9], Q2m[9], Q2[9];
    rodrigues(ax1, 0.7, Q1);            // committed: rotated 0.7 rad
    rodrigues(ax2, 1.9, Q2m);           // increment: FURTHER 1.9 rad — huge step
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double s = 0.0;
        for (int k = 0; k < 3; k++) s += Q2m[3*i+k] * Q1[3*k+j];
        Q2[3*i+j] = s;
      }
    double un[24], u1[24], deps[6], R1[9], J1;
    affineDisp(Q1, un);
    affineDisp(Q2, u1);
    check(hypoIncrement(8, dNdX, un, u1, deps, R1, J1), "rigid: hypoIncrement ok");
    double zero6[6] = {0, 0, 0, 0, 0, 0};
    check(maxAbsDiff(deps, zero6, 6) < 1e-13, "rigid: deps == 0 at 1.9 rad/step");
    check(std::fabs(J1 - 1.0) < 1e-13, "rigid: J1 == 1");
    check(maxAbsDiff(R1, Q2, 9) < 1e-12, "rigid: R1 == Q2");
  }

  // -- 2. pure stretch from rest, one step: deps_ii = 2(l-1)/(l+1) ------------
  {
    double A[9] = {1.30, 0, 0,  0, 0.85, 0,  0, 0, 1.10};
    double un[24] = {0}, u1[24], deps[6], R1[9], J1;
    affineDisp(A, u1);
    check(hypoIncrement(8, dNdX, un, u1, deps, R1, J1), "stretch: ok");
    double lam[3] = {1.30, 0.85, 1.10};
    for (int i = 0; i < 3; i++) {
      double want = 2.0 * (lam[i] - 1.0) / (lam[i] + 1.0);
      check(std::fabs(deps[i] - want) < 1e-14, "stretch: midpoint closed form");
    }
    check(std::fabs(deps[3]) + std::fabs(deps[4]) + std::fabs(deps[5]) < 1e-14,
          "stretch: no shear");
    double I9[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    check(maxAbsDiff(R1, I9, 9) < 1e-12, "stretch: R1 == I");
    check(std::fabs(J1 - 1.30 * 0.85 * 1.10) < 1e-13, "stretch: J1 = prod(l)");
  }

  // -- 3. polar orthogonality on a rotated stretch ----------------------------
  {
    const double ax[3] = {0.6, -0.64, 0.48};
    double Q[9];
    rodrigues(ax, 2.4, Q);
    double U[9] = {1.2, 0.1, -0.05,  0.1, 0.9, 0.08,  -0.05, 0.08, 1.05};
    double F[9];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double s = 0.0;
        for (int k = 0; k < 3; k++) s += Q[3*i+k] * U[3*k+j];
        F[3*i+j] = s;
      }
    double R[9];
    check(polar3(F, R), "polar: ok");
    double RtR[9];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double s = 0.0;
        for (int k = 0; k < 3; k++) s += R[3*k+i] * R[3*k+j];
        RtR[3*i+j] = s;
      }
    double I9[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    check(maxAbsDiff(RtR, I9, 9) < 1e-13, "polar: R^T R == I");
    check(std::fabs(det3(R) - 1.0) < 1e-13, "polar: det R == 1");
    check(maxAbsDiff(R, Q, 9) < 1e-12, "polar: R == Q for F = Q U (U symm PD)");
  }

  // -- 4. bond6(I) = identity; tangent push preserves symmetry ----------------
  {
    double I9[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double M[36];
    bond6(I9, M);
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        check(std::fabs(M[6*r+c] - (r == c ? 1.0 : 0.0)) < 1e-15, "bond6(I)");

    const double ax[3] = {0.0, 0.0, 1.0};
    double Q[9];
    rodrigues(ax, 0.83, Q);
    double D[36], Dout[36];
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        D[6*r+c] = (r == c ? 5.0 : 0.0) + 0.1 * (r + c);   // symmetric
    pushTangent6(Q, D, Dout);
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        check(std::fabs(Dout[6*r+c] - Dout[6*c+r]) < 1e-12, "tangent push symm");

    // isotropy invariance: an isotropic 6x6 (eng-shear convention) is
    // rotation-invariant — the sharpest check of the bond convention.
    double lamK = 2.0, muK = 0.7, Diso[36] = {0}, DisoOut[36];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Diso[6*i+j] = lamK + (i == j ? 2.0 * muK : 0.0);
    for (int i = 3; i < 6; i++) Diso[6*i+i] = muK;   // GAMMA convention: G, not 2G
    const double ax2[3] = {0.48, 0.6, -0.64};
    double Q2[9];
    rodrigues(ax2, 1.2, Q2);
    pushTangent6(Q2, Diso, DisoOut);
    check(maxAbsDiff(Diso, DisoOut, 36) < 1e-12, "isotropic tangent invariant");
  }
}

// ---------------------------------------------------------------------------
// (b) oracle cases
// ---------------------------------------------------------------------------
static bool readVals(std::istringstream &ss, double *dst, int n)
{
  for (int i = 0; i < n; i++)
    if (!(ss >> dst[i])) return false;
  return true;
}

static int runCases(const char *path)
{
  std::ifstream in(path);
  if (!in) {
    std::printf("[cases] no cases file at %s — skipping oracle diff\n", path);
    std::printf("[cases] parsed 0 cases\n");
    return 0;
  }

  int nCases = 0;
  double dNdX[24], un[24], u1[24], D[36], s6[6];
  double eDeps[6], eR1[9], eJ1, eSout[6], eDout[36];
  bool have[10] = {false};

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string tag;
    ss >> tag;
    if (tag == "case") {
      std::memset(have, 0, sizeof(have));
    } else if (tag == "dNdX") { have[0] = readVals(ss, dNdX, 24); }
    else if (tag == "un")   { have[1] = readVals(ss, un, 24); }
    else if (tag == "u1")   { have[2] = readVals(ss, u1, 24); }
    else if (tag == "D")    { have[3] = readVals(ss, D, 36); }
    else if (tag == "s")    { have[4] = readVals(ss, s6, 6); }
    else if (tag == "deps") { have[5] = readVals(ss, eDeps, 6); }
    else if (tag == "R1")   { have[6] = readVals(ss, eR1, 9); }
    else if (tag == "J1")   { have[7] = readVals(ss, &eJ1, 1); }
    else if (tag == "sout") { have[8] = readVals(ss, eSout, 6); }
    else if (tag == "Dout") { have[9] = readVals(ss, eDout, 36); }
    else if (tag == "endcase") {
      bool all = true;
      for (int i = 0; i < 10; i++) all = all && have[i];
      if (!all) {
        std::printf("FAIL: case %d incomplete\n", nCases);
        g_fail++;
        continue;
      }
      double deps[6], R1[9], J1, sout[6], Dout[36];
      if (!hypoIncrement(8, dNdX, un, u1, deps, R1, J1)) {
        std::printf("FAIL: case %d hypoIncrement rejected\n", nCases);
        g_fail++;
        nCases++;
        continue;
      }
      pushStress6(R1, s6, sout);
      pushTangent6(R1, D, Dout);
      const double tol = 1e-9;
      if (maxAbsDiff(deps, eDeps, 6) > tol) { std::printf("FAIL: case %d deps (%.3e)\n", nCases, maxAbsDiff(deps, eDeps, 6)); g_fail++; }
      if (maxAbsDiff(R1, eR1, 9) > tol)     { std::printf("FAIL: case %d R1 (%.3e)\n", nCases, maxAbsDiff(R1, eR1, 9)); g_fail++; }
      if (std::fabs(J1 - eJ1) > tol)        { std::printf("FAIL: case %d J1\n", nCases); g_fail++; }
      if (maxAbsDiff(sout, eSout, 6) > tol) { std::printf("FAIL: case %d sout (%.3e)\n", nCases, maxAbsDiff(sout, eSout, 6)); g_fail++; }
      if (maxAbsDiff(Dout, eDout, 36) > tol){ std::printf("FAIL: case %d Dout (%.3e)\n", nCases, maxAbsDiff(Dout, eDout, 36)); g_fail++; }
      nCases++;
    }
  }
  std::printf("[cases] parsed %d cases\n", nCases);
  return nCases;
}

int main(int argc, char **argv)
{
  selfTests();
  runCases(argc > 1 ? argv[1] : "tests/hypo_cases.txt");
  if (g_fail == 0) {
    std::printf("RESULT: PASS\n");
    return 0;
  }
  std::printf("RESULT: FAIL (%d)\n", g_fail);
  return 1;
}
