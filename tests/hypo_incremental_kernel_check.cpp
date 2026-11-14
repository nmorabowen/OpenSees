// hypo_incremental_kernel_check — standalone g++ check of the incremental
// hypo-kernel entry (WP-232, tests/mpm_interface.md §8): hypoIncrementFromF and
// hypoIncrementFromFPlaneStrain in LadrunoHypoKernel.h. No OpenSees dependency.
//
//   g++ -O2 -std=c++17 -Wall -Wextra -I SRC/element/solidTransformation
//       tests/hypo_incremental_kernel_check.cpp -o hinc ; ./hinc tests/hypo_incremental_cases.txt
//
// (a) Self-tests: the contract invariants (rigid rotation => deps = 0 at any
//     angle below pi; the pure-stretch midpoint closed form 2e/(2+e); R1 =
//     Q·polar(Fn); refusals leave outputs untouched; the plane-strain entry
//     equals the 3D entry on the embedded F) and an AGREEMENT check against the
//     ADR-79 nodal-displacement entry hypoIncrement() on affine H8 fields — the
//     two entries are the same algebra by different arithmetic.
// (b) Every case in tests/hypo_incremental_cases.txt (WP-232b, written without
//     reading SRC/): status, deps, R1, J1 and, when given, the stress push
//     scauchy = R1 shat R1^T. A value marked `pending-oracle` is not compared
//     and is counted.
//
// Exit 0 + "RESULT: PASS" on success.

#include <LadrunoHypoKernel.h>

#include <cmath>
#include <cstdio>
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

static void matmul3(const double A[9], const double B[9], double C[9])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += A[3*i+k] * B[3*k+j];
      C[3*i+j] = s;
    }
}

static void rodrigues(const double axis[3], double th, double Q[9])
{
  double c = std::cos(th), s = std::sin(th), v = 1.0 - c;
  double x = axis[0], y = axis[1], z = axis[2];
  Q[0] = c + x*x*v;   Q[1] = x*y*v - z*s; Q[2] = x*z*v + y*s;
  Q[3] = y*x*v + z*s; Q[4] = c + y*y*v;   Q[5] = y*z*v - x*s;
  Q[6] = z*x*v - y*s; Q[7] = z*y*v + x*s; Q[8] = c + z*z*v;
}

static const double CORNERS[8][3] = {
  {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
  {-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1},
};

// u_a = (A - I) X_a for the affine map x = A X on the [-1,1]^3 cube; centroid
// gradients dN_a/dX = corner_a / 8 reproduce F = A exactly.
static void affineDisp(const double A[9], double u[24])
{
  for (int a = 0; a < 8; a++)
    for (int i = 0; i < 3; i++) {
      double s = 0.0;
      for (int J = 0; J < 3; J++) s += A[3*i+J] * CORNERS[a][J];
      u[a*3+i] = s - CORNERS[a][i];
    }
}

static void selfTests()
{
  const double I9[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  const double zero6[6] = {0, 0, 0, 0, 0, 0};

  // 1. rigid rotation from a stretched state: deps = 0, R1 = Q polar(Fn) = Q U-less
  {
    const double ax[3] = {0.267261241912424, 0.534522483824849, 0.801783725737273};
    const double Fn[9] = {1.2, 0.1, 0.0,  -0.05, 0.9, 0.02,  0.0, 0.03, 1.1};
    double Rn[9];
    check(polar3(Fn, Rn), "rigid: polar(Fn)");
    for (double th : {0.1, 1.0, 2.5, 3.0}) {
      double Q[9], F1[9], QRn[9];
      rodrigues(ax, th, Q);
      matmul3(Q, Fn, F1);
      matmul3(Q, Rn, QRn);
      double deps[6], R1[9], J1;
      check(hypoIncrementFromF(Fn, F1, deps, R1, &J1) == 0, "rigid: status 0");
      check(maxAbsDiff(deps, zero6, 6) < 1e-12, "rigid: deps == 0 below pi");
      check(maxAbsDiff(R1, QRn, 9) < 1e-12, "rigid: R1 == Q polar(Fn)");
      check(std::fabs(J1 - det3(Fn)) < 1e-12, "rigid: J1 == det Fn");
    }
  }

  // 2. pure stretch from rest: deps_ii = 2(l-1)/(l+1)
  {
    const double F1[9] = {1.30, 0, 0,  0, 0.85, 0,  0, 0, 1.10};
    double deps[6], R1[9], J1;
    check(hypoIncrementFromF(I9, F1, deps, R1, &J1) == 0, "stretch: status 0");
    const double lam[3] = {1.30, 0.85, 1.10};
    for (int i = 0; i < 3; i++)
      check(std::fabs(deps[i] - 2.0 * (lam[i] - 1.0) / (lam[i] + 1.0)) < 1e-14,
            "stretch: midpoint closed form");
    check(std::fabs(deps[3]) + std::fabs(deps[4]) + std::fabs(deps[5]) < 1e-15, "stretch: no shear");
    check(maxAbsDiff(R1, I9, 9) < 1e-14, "stretch: R1 == I");
  }

  // 3. refusals leave the outputs untouched
  {
    const double Fbad[9] = {-1, 0, 0,  0, 1, 0,  0, 0, 1};          // det F1 < 0
    const double Fflip[9] = {-1, 0, 0,  0, -1, 0,  0, 0, 1};        // det Fmid = 0
    for (const double *F1 : {Fbad, Fflip}) {
      double deps[6] = {7, 7, 7, 7, 7, 7}, R1[9] = {7, 7, 7, 7, 7, 7, 7, 7, 7}, J1 = 7.0;
      check(hypoIncrementFromF(I9, F1, deps, R1, &J1) == -1, "refuse: status -1");
      bool untouched = (J1 == 7.0);
      for (int i = 0; i < 6; i++) untouched = untouched && deps[i] == 7.0;
      for (int i = 0; i < 9; i++) untouched = untouched && R1[i] == 7.0;
      check(untouched, "refuse: outputs untouched");
    }
  }

  // 4. agreement with the ADR-79 nodal-displacement entry on affine fields
  {
    double dNdX[24];
    for (int a = 0; a < 8; a++)
      for (int j = 0; j < 3; j++) dNdX[a*3+j] = CORNERS[a][j] / 8.0;
    const double Fn[9] = {1.05, 0.12, -0.03,  0.02, 0.97, 0.08,  -0.06, 0.04, 1.01};
    const double ax[3] = {0.6, -0.64, 0.48};
    double Q[9], S[9], F1[9];
    rodrigues(ax, 0.7, Q);
    const double D[9] = {1.1, 0.05, 0.0,  0.0, 0.92, -0.04,  0.03, 0.0, 1.07};
    matmul3(D, Fn, S);
    matmul3(Q, S, F1);
    double un[24], u1[24];
    affineDisp(Fn, un);
    affineDisp(F1, u1);
    double dA[6], RA[9], JA, dB[6], RB[9], JB;
    check(hypoIncrement(8, dNdX, un, u1, dA, RA, JA), "agree: nodal entry ok");
    check(hypoIncrementFromF(Fn, F1, dB, RB, &JB) == 0, "agree: F entry ok");
    check(maxAbsDiff(dA, dB, 6) < 1e-13, "agree: deps matches hypoIncrement");
    check(maxAbsDiff(RA, RB, 9) < 1e-13, "agree: R1 matches hypoIncrement");
    check(std::fabs(JA - JB) < 1e-13, "agree: J1 matches hypoIncrement");
  }

  // 5. plane strain equals the 3D entry on the embedded F; push agrees too
  {
    const double Fn[4] = {1.1, 0.2, -0.05, 0.95};
    const double F1[4] = {0.9, 0.45, -0.3, 1.02};
    const double Fn3[9] = {Fn[0], Fn[1], 0, Fn[2], Fn[3], 0, 0, 0, 1};
    const double F13[9] = {F1[0], F1[1], 0, F1[2], F1[3], 0, 0, 0, 1};
    double d3[3], R2[4], J2, d6[6], R3[9], J3;
    check(hypoIncrementFromFPlaneStrain(Fn, F1, d3, R2, &J2) == 0, "ps: status 0");
    check(hypoIncrementFromF(Fn3, F13, d6, R3, &J3) == 0, "ps: 3D status 0");
    check(d6[2] == 0.0 && d6[4] == 0.0 && d6[5] == 0.0, "ps: out-of-plane deps exactly 0");
    check(d3[0] == d6[0] && d3[1] == d6[1] && d3[2] == d6[3], "ps: deps == embedded 3D");
    check(R2[0] == R3[0] && R2[1] == R3[1] && R2[2] == R3[3] && R2[3] == R3[4], "ps: R1 block");
    check(R3[2] == 0.0 && R3[5] == 0.0 && R3[6] == 0.0 && R3[7] == 0.0 && R3[8] == 1.0,
          "ps: R1 is a rotation about z exactly");
    check(J2 == J3, "ps: J1");
    const double s3[3] = {1.5, -2.0, 0.8};
    const double s6[6] = {1.5, -2.0, 0.0, 0.8, 0.0, 0.0};
    double o3[3], o6[6];
    pushStressPlaneStrain(R2, s3, o3);
    pushStress6(R3, s6, o6);
    check(std::fabs(o3[0] - o6[0]) + std::fabs(o3[1] - o6[1]) + std::fabs(o3[2] - o6[3]) < 1e-14,
          "ps: pushStressPlaneStrain == pushStress6 embedded");
    // WP-232f: the plane-strain tangent push equals the {xx,yy,xy} block of
    // pushTangent6 on a 6x6 tangent whose out-of-plane rows/columns are arbitrary.
    const double D3[9] = {7.0, 2.5, 0.3,  2.1, 6.0, -0.4,  0.2, -0.6, 1.9};
    double D6[36];
    for (int k = 0; k < 36; k++) D6[k] = 0.37 * (k % 7) - 0.5;   // out-of-plane junk
    const int ix[3] = {0, 1, 3};
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++) D6[6 * ix[r] + ix[c]] = D3[3 * r + c];
    double O3[9], O6[36];
    pushTangentPlaneStrain(R2, D3, O3);
    pushTangent6(R3, D6, O6);
    double dmax = 0.0;
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        dmax = std::max(dmax, std::fabs(O3[3 * r + c] - O6[6 * ix[r] + ix[c]]));
    check(dmax < 1e-13, "ps: pushTangentPlaneStrain == pushTangent6 embedded block");
  }
}

// ---------------------------------------------------------------------------
// (b) the WP-232b cases
// ---------------------------------------------------------------------------
struct Field {
  bool present = false;
  bool pending = false;
  std::vector<double> v;
};

static int g_pendingSkipped = 0;

static Field readField(std::istringstream &ss)
{
  Field f;
  f.present = true;
  std::string tok;
  while (ss >> tok) {
    if (tok == "pending-oracle") { f.pending = true; continue; }
    f.v.push_back(std::stod(tok));
  }
  return f;
}

static bool within(const double *got, const std::vector<double> &exp, int n,
                   double atol, double rtol, double &worst)
{
  worst = 0.0;
  bool ok = (int)exp.size() == n;
  for (int i = 0; ok && i < n; i++) {
    double d = std::fabs(got[i] - exp[i]);
    if (d > worst) worst = d;
    if (d > atol + rtol * std::fabs(exp[i])) ok = false;
  }
  return ok;
}

static void compareField(const std::string &name, const char *key, const Field &f,
                         const double *got, int n, double atol, double rtol)
{
  if (!f.present) return;
  if (f.pending) { g_pendingSkipped++; return; }
  double worst;
  if (!within(got, f.v, n, atol, rtol, worst)) {
    std::printf("FAIL: case %s %s (max |diff| %.3e, %d values expected %d)\n",
                name.c_str(), key, worst, (int)f.v.size(), n);
    g_fail++;
  }
}

static int runCases(const char *path)
{
  std::ifstream in(path);
  if (!in) {
    std::printf("FAIL: no cases file at %s\n", path);
    g_fail++;
    return 0;
  }
  int nCases = 0;
  std::string line, name;
  int dim = 0, status = 0;
  double atol = 1e-12, rtol = 1e-12;
  Field Fn, F1, deps, R1, J1, shat, scauchy;
  auto resetCase = [&]() {
    name.clear(); dim = 0; status = 0; atol = 1e-12; rtol = 1e-12;
    Fn = F1 = deps = R1 = J1 = shat = scauchy = Field();
  };
  resetCase();
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key)) continue;
    if (key == "case") { resetCase(); continue; }
    if (key == "name") { ss >> name; continue; }
    if (key == "note") continue;
    if (key == "dim") { ss >> dim; continue; }
    if (key == "status") { ss >> status; continue; }
    if (key == "tol") { ss >> atol >> rtol; continue; }
    if (key == "Fn") { Fn = readField(ss); continue; }
    if (key == "F1") { F1 = readField(ss); continue; }
    if (key == "deps") { deps = readField(ss); continue; }
    if (key == "R1") { R1 = readField(ss); continue; }
    if (key == "J1") { J1 = readField(ss); continue; }
    if (key == "shat") { shat = readField(ss); continue; }
    if (key == "scauchy") { scauchy = readField(ss); continue; }
    if (key != "endcase") {
      std::printf("FAIL: unknown key '%s'\n", key.c_str());
      g_fail++;
      continue;
    }
    // ---- endcase: run it ----
    const int nn = dim * dim, nv = (dim == 3) ? 6 : 3;
    if ((dim != 2 && dim != 3) || (int)Fn.v.size() != nn || (int)F1.v.size() != nn) {
      std::printf("FAIL: case %s malformed\n", name.c_str());
      g_fail++;
      continue;
    }
    double gd[6] = {9, 9, 9, 9, 9, 9}, gR[9] = {9, 9, 9, 9, 9, 9, 9, 9, 9}, gJ = 9.0;
    int st = (dim == 3) ? hypoIncrementFromF(Fn.v.data(), F1.v.data(), gd, gR, &gJ)
                        : hypoIncrementFromFPlaneStrain(Fn.v.data(), F1.v.data(), gd, gR, &gJ);
    if (st != status) {
      std::printf("FAIL: case %s status %d, expected %d\n", name.c_str(), st, status);
      g_fail++;
    } else if (st == 0) {
      compareField(name, "deps", deps, gd, nv, atol, rtol);
      compareField(name, "R1", R1, gR, nn, atol, rtol);
      compareField(name, "J1", J1, &gJ, 1, atol, rtol);
      if (shat.present && !shat.pending && (int)shat.v.size() == nv) {
        double sc[6];
        if (dim == 3) pushStress6(gR, shat.v.data(), sc);
        else pushStressPlaneStrain(gR, shat.v.data(), sc);
        compareField(name, "scauchy", scauchy, sc, nv, atol, rtol);
      }
    } else {
      bool untouched = (gJ == 9.0);
      for (int i = 0; i < 6; i++) untouched = untouched && gd[i] == 9.0;
      for (int i = 0; i < 9; i++) untouched = untouched && gR[i] == 9.0;
      if (!untouched) {
        std::printf("FAIL: case %s refused but wrote its outputs\n", name.c_str());
        g_fail++;
      }
    }
    nCases++;
  }
  std::printf("[cases] parsed %d cases\n", nCases);
  std::printf("[cases] pending-oracle values skipped: %d\n", g_pendingSkipped);
  return nCases;
}

int main(int argc, char **argv)
{
  selfTests();
  runCases(argc > 1 ? argv[1] : "tests/hypo_incremental_cases.txt");
  if (g_fail == 0) {
    std::printf("RESULT: PASS\n");
    return 0;
  }
  std::printf("RESULT: FAIL (%d)\n", g_fail);
  return 1;
}
