// Standalone g++ verification of Ladruno::hrzLumpRaw (OpenSees-free).
// Build: g++ -std=c++14 -I../../SRC/analysis/integrator -DLADRUNO_MASSLUMPING_STANDALONE hrz_standalone.cpp -o hrz_standalone
// Run in CI via tests/test_hrz_lumped_mass.py::test_hrz_standalone_kernel
// (compiles + runs this file). Exit 0 = all HRZ invariants hold.

#ifndef LADRUNO_MASSLUMPING_STANDALONE
#define LADRUNO_MASSLUMPING_STANDALONE
#endif
#include <LadrunoMassLumping.h>

#include <cstdio>
#include <cmath>
#include <vector>

using std::vector;

static int g_fail = 0;
static void check(bool cond, const char *msg) {
  std::printf("  [%s] %s\n", cond ? "PASS" : "FAIL", msg);
  if (!cond) g_fail++;
}
static bool close(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(b)); }

static double rowsum(const vector<double> &Mc, int n, int i) {
  double s = 0.0; for (int j = 0; j < n; ++j) s += Mc[i * n + j]; return s;
}
static double dirTotal(const vector<double> &Mc, const vector<int> &dir, int n, int d) {
  double s = 0.0; for (int i = 0; i < n; ++i) if (dir[i] == d) s += rowsum(Mc, n, i); return s;
}
static double dirSumMdiag(const vector<double> &md, const vector<int> &dir, int n, int d) {
  double s = 0.0; for (int i = 0; i < n; ++i) if (dir[i] == d) s += md[i]; return s;
}

int main() {
  // --- Case 1: regular 1D bar, HRZ must EQUAL row-sum -----------------------
  {
    int n = 2; vector<double> Mc = {2,1, 1,2}; vector<int> dir = {0,0};
    vector<double> md(n);
    int rc = Ladruno::hrzLumpRaw(Mc.data(), dir.data(), n, md.data());
    std::printf("Case1 regular bar: rc=%d md=[%g,%g]\n", rc, md[0], md[1]);
    check(rc == Ladruno::HRZ_OK, "rc==OK");
    check(close(md[0], rowsum(Mc, n, 0)) && close(md[1], rowsum(Mc, n, 1)), "HRZ == row-sum on regular bar");
    check(close(md[0] + md[1], dirTotal(Mc, dir, n, 0)), "conserves total mass");
  }
  // --- Case 2: distorted 1D, HRZ != row-sum but conserves -------------------
  {
    int n = 2; vector<double> Mc = {2,1, 1,4}; vector<int> dir = {0,0};
    vector<double> md(n);
    int rc = Ladruno::hrzLumpRaw(Mc.data(), dir.data(), n, md.data());
    std::printf("Case2 distorted: rc=%d md=[%g,%g]\n", rc, md[0], md[1]);
    check(rc == Ladruno::HRZ_OK, "rc==OK");
    check(!close(md[0], rowsum(Mc, n, 0)), "HRZ != row-sum (per-node) on distorted");
    check(close(md[0] + md[1], dirTotal(Mc, dir, n, 0)), "conserves total mass (=8)");
    check(close(md[0] + md[1], 8.0), "total == 8");
  }
  // --- Case 3: regular square quad (4 nodes x/y), HRZ == row-sum ------------
  {
    int n = 8; vector<int> dir = {0,1,0,1,0,1,0,1};
    // per-direction 4x4 consistent quad block, replicated into x and y DOFs
    double q[4][4] = {{4,2,1,2},{2,4,2,1},{1,2,4,2},{2,1,2,4}};
    vector<double> Mc(n * n, 0.0);
    for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) {
      Mc[(2*a + 0) * n + (2*b + 0)] = q[a][b];   // x DOFs
      Mc[(2*a + 1) * n + (2*b + 1)] = q[a][b];   // y DOFs
    }
    vector<double> md(n);
    int rc = Ladruno::hrzLumpRaw(Mc.data(), dir.data(), n, md.data());
    std::printf("Case3 regular quad: rc=%d md0..3=[%g,%g,%g,%g]\n", rc, md[0], md[1], md[2], md[3]);
    check(rc == Ladruno::HRZ_OK, "rc==OK");
    bool eqRowsum = true; for (int i = 0; i < n; ++i) if (!close(md[i], rowsum(Mc, n, i))) eqRowsum = false;
    check(eqRowsum, "HRZ == row-sum on regular quad (all DOFs)");
    check(close(dirSumMdiag(md, dir, n, 0), dirTotal(Mc, dir, n, 0)), "conserves x-mass");
    check(close(dirSumMdiag(md, dir, n, 1), dirTotal(Mc, dir, n, 1)), "conserves y-mass");
  }
  // --- Case 4: 2D beam, rotational DOFs + direction-dependent alpha ---------
  {
    int n = 6; vector<int> dir = {0,1,-1, 0,1,-1};   // ux,uy,rz per 2 nodes
    vector<double> Mc(n * n, 0.0);
    // x block (DOFs 0,3): alpha_x driver
    Mc[0*n+0]=140; Mc[0*n+3]=70; Mc[3*n+0]=70; Mc[3*n+3]=140;
    // y block (DOFs 1,4): different -> alpha_y != alpha_x
    Mc[1*n+1]=156; Mc[1*n+4]=54; Mc[4*n+1]=54; Mc[4*n+4]=156;
    // rot block (DOFs 2,5): positive diagonal
    Mc[2*n+2]=4;   Mc[2*n+5]=-3; Mc[5*n+2]=-3; Mc[5*n+5]=4;
    vector<double> md(n);
    int rc = Ladruno::hrzLumpRaw(Mc.data(), dir.data(), n, md.data());
    std::printf("Case4 beam: rc=%d md=[%g,%g,%g,%g,%g,%g]\n", rc, md[0],md[1],md[2],md[3],md[4],md[5]);
    double ax = 420.0/280.0, ay = 420.0/312.0, amean = (ax + ay) / 2.0;
    check(rc == Ladruno::HRZ_OK, "rc==OK");
    check(close(md[0], ax*140) && close(md[3], ax*140), "x DOFs scaled by alpha_x");
    check(close(md[1], ay*156) && close(md[4], ay*156), "y DOFs scaled by alpha_y");
    check(close(md[2], amean*4) && close(md[5], amean*4), "rot DOFs scaled by mean(alpha)");
    check(md[2] > 0 && md[5] > 0, "rotational inertia strictly positive (row-sum would give 1)");
    check(close(md[0]+md[3], dirTotal(Mc,dir,n,0)) && close(md[1]+md[4], dirTotal(Mc,dir,n,1)), "conserves x and y mass");
  }
  // --- Case 5: guard -> fallback to Diagonal on non-positive consistent diag -
  {
    int n = 2; vector<double> Mc = {0,1, 1,2}; vector<int> dir = {0,0};  // Mc_00 = 0
    vector<double> md(n);
    int rc = Ladruno::hrzLumpRaw(Mc.data(), dir.data(), n, md.data());
    std::printf("Case5 guard: rc=%d md=[%g,%g]\n", rc, md[0], md[1]);
    check(rc == Ladruno::HRZ_FELLBACK_DIAG, "rc==FELLBACK_DIAG when consistent diagonal <= 0");
    check(close(md[0], 0.0) && close(md[1], 2.0), "fallback returns diagonal-of-consistent");
  }

  std::printf("\n%s (%d failures)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
  return g_fail == 0 ? 0 : 1;
}
