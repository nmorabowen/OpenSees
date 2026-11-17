// WP-232j: the header-level checks of anti-locking (tests/mpm_antilocking_interface.md J1-J5;
// gate rows L1 a, b, c of wp-briefs.md "232j"). Standalone: g++ -std=c++17 -I <ladrunoMPM>
// -I <solidTransformation>; run by tests/test_ladruno_mpm_antilocking_kinematics.py.
//
//  L1a  an affine field du_a = A x_a + b: every particle's dFbar equals dF, rel 1e-12
//       (|dFbar - dF|_max / |dF - I|_max), every basis, 2D and 3D, untruncated supports;
//  L1b  a 30 deg rigid rotation A = Q - I: dFbar = Q to rounding, the kernel's deps at the
//       frozen rigid_rotation tolerance (tests/hypo_incremental_cases.txt: 1e-12);
//  L1c  plane strain: the 2x2 scaling, embedded with F_zz = 1, gives the 3D kernel's
//       increment with deps_zz = deps_yz = deps_zx = 0 exactly;
//  J3   the corners of a particle's cell lie in its support for every basis, including
//       positions on grid lines and on the max face; the corner gradients reproduce I;
//  J1   det(dFbar) = det(dF_c); a non-positive determinant refuses with dFbar untouched;
//       a non-affine field is changed (the averaging is not the identity).

#include <LadrunoMPMCore.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace ladruno_mpm;

static int g_fail = 0, g_checks = 0;

static void check(bool ok, const std::string &what)
{
  g_checks++;
  if (!ok) {
    g_fail++;
    std::printf("FAIL: %s\n", what.c_str());
  }
}

static const char *BN[4] = {"gimp", "linear", "bspline2", "bspline3"};

static bool setup(Engine &E, int nd, int basis, int ppc)
{
  std::string err;
  const int n[3] = {6, 6, 6};
  const double x0[3] = {-1.0, 0.5, 2.0};
  if (!E.setGrid(nd, x0, 0.5, n, err)) return false;
  E.par.basis = (Basis)basis;
  E.par.antiLocking = 1;
  std::vector<double> geom;
  for (int d = 0; d < nd; d++) geom.push_back(x0[d] + 1.0);
  for (int d = 0; d < nd; d++) geom.push_back(x0[d] + 2.0);
  int64_t next = 0;
  return E.seedLattice(0, geom.data(), ppc, 1, 1, 1000.0, next, err) > 0;
}

static void rotation(int nd, double deg, double Q[9])
{
  const double t = deg * 3.14159265358979323846 / 180.0, c = std::cos(t), s = std::sin(t);
  for (int k = 0; k < 9; k++) Q[k] = 0.0;
  if (nd == 2) { Q[0] = c; Q[1] = -s; Q[2] = s; Q[3] = c; return; }
  // about the axis (1, 2, 2)/3 (Rodrigues)
  const double a[3] = {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
  const double K[9] = {0, -a[2], a[1], a[2], 0, -a[0], -a[1], a[0], 0};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double kk = 0.0;
      for (int m = 0; m < 3; m++) kk += K[3 * i + m] * K[3 * m + j];
      Q[3 * i + j] = (i == j ? 1.0 : 0.0) + s * K[3 * i + j] + (1.0 - c) * kk;
    }
}

int main()
{
  // ---- L1a / L1b through Engine::affineCheck, every basis, 2D and 3D -----------
  for (int nd = 2; nd <= 3; nd++)
    for (int basis = 0; basis < 4; basis++)
      for (int ppc = 2; ppc <= 3; ppc++) {
        Engine E;
        const std::string tag = std::string(BN[basis]) + " " + std::to_string(nd) + "D ppc " +
                                std::to_string(ppc);
        if (!setup(E, nd, basis, ppc)) { check(false, tag + ": setup"); continue; }
        const int N = E.P.size();
        double A[9], b[3] = {0.3, -0.2, 0.1};
        std::vector<double> out;
        // deviatoric (shear + stretch, trace 0)
        for (int k = 0; k < 9; k++) A[k] = 0.0;
        A[0] = 1e-2; A[nd + 1] = -1e-2; A[1] = 7e-3; A[nd] = -3e-3;
        E.affineCheck(A, b, out);
        check(out[0] <= 1e-12 && out[2] == N && out[3] == 0,
              tag + ": L1a deviatoric rel " + std::to_string(out[0]));
        // volumetric
        for (int k = 0; k < 9; k++) A[k] = 0.0;
        for (int d = 0; d < nd; d++) A[(nd + 1) * d] = 1.2e-2;
        E.affineCheck(A, b, out);
        check(out[0] <= 1e-12 && out[2] == N && out[3] == 0,
              tag + ": L1a volumetric rel " + std::to_string(out[0]));
        // general affine (deviatoric + volumetric + spin)
        for (int k = 0; k < nd * nd; k++) A[k] = 1e-2 * std::sin(1.3 * k + 0.4);
        E.affineCheck(A, b, out);
        check(out[0] <= 1e-12 && out[2] == N && out[3] == 0,
              tag + ": L1a general rel " + std::to_string(out[0]));
        // L1b: 30 deg rigid rotation
        double Q[9];
        rotation(nd, 30.0, Q);
        for (int i = 0; i < nd; i++)
          for (int j = 0; j < nd; j++) A[nd * i + j] = Q[nd * i + j] - (i == j ? 1.0 : 0.0);
        E.affineCheck(A, b, out);
        check(out[0] <= 1e-12 && out[1] <= 1e-12 && out[2] == N,
              tag + ": L1b rotation rel " + std::to_string(out[0]) + " deps " + std::to_string(out[1]));
      }

  // ---- J3: corners inside the support; gradients reproduce I -----------------------
  for (int nd = 2; nd <= 3; nd++)
    for (int basis = 0; basis < 4; basis++) {
      Grid g;
      g.ndm = nd; g.h = 0.25;
      for (int d = 0; d < 3; d++) { g.x0[d] = 0.1 * d; g.n[d] = (d < nd) ? 5 : 0; }
      int bad = 0, badI = 0;
      for (int t = 0; t < 4000; t++) {
        double x[3], lp[3];
        for (int d = 0; d < nd; d++) {
          double f = std::fmod(0.6180339887 * (t + 1) * (d + 1.7), 1.0);
          if (t % 10 == 0) f = std::floor(f * 5.0) / 5.0;      // on a grid line
          if (t % 37 == 0) f = 1.0;                              // the max face
          x[d] = g.x0[d] + f * g.n[d] * g.h;
          lp[d] = 0.25 * g.h;
        }
        std::vector<long> node;
        std::vector<double> S, dS;
        particleSupport(g, (Basis)basis, x, lp, node, S, dS);
        long ce[8];
        double gN[8][3];
        const int nc = cellCorners(g, x, ce, gN);
        for (int k = 0; k < nc; k++) {
          bool in = false;
          for (long nn : node) in = in || nn == ce[k];
          if (!in) bad++;
        }
        double M[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int k = 0; k < nc; k++) {
          int ijk[3];
          g.extIJK(ce[k], ijk);
          for (int a = 0; a < nd; a++)
            for (int c = 0; c < nd; c++) M[nd * a + c] += (g.x0[a] + ijk[a] * g.h) * gN[k][c];
        }
        for (int a = 0; a < nd; a++)
          for (int c = 0; c < nd; c++)
            if (std::fabs(M[nd * a + c] - (a == c ? 1.0 : 0.0)) > 1e-12) badI++;
      }
      check(bad == 0, std::string("J3 corners in support, ") + BN[basis] + " " + std::to_string(nd) +
                      "D: " + std::to_string(bad) + " misses");
      check(badI == 0, std::string("J1 corner gradients reproduce I, ") + BN[basis] + " " +
                       std::to_string(nd) + "D");
    }

  // ---- L1c, J1: plane strain embedding, det, refusal, non-identity --------------------
  for (int t = 0; t < 200; t++) {
    double dF[4], dFc[4], dFb[4];
    for (int k = 0; k < 4; k++) {
      dF[k] = (k % 3 == 0 ? 1.0 : 0.0) + 0.03 * std::sin(2.1 * t + k);
      dFc[k] = (k % 3 == 0 ? 1.0 : 0.0) + 0.03 * std::cos(1.7 * t + 0.5 * k);
    }
    if (fbarIncrement(2, dF, dFc, dFb) != 0) { check(false, "L1c: refused a valid pair"); continue; }
    check(std::fabs(detN(2, dFb) - detN(2, dFc)) <= 1e-14 * detN(2, dFc), "J1: det(dFbar) = det(dFc), 2D");
    // embed with F_zz = 1 and run both kernel entries from a plane-strain F_n
    const double Fn2[4] = {1.02, 0.01, -0.02, 0.97};
    double F12[4], F13[9], Fn3[9] = {Fn2[0], Fn2[1], 0, Fn2[2], Fn2[3], 0, 0, 0, 1};
    for (int a = 0; a < 2; a++)
      for (int c = 0; c < 2; c++) F12[2 * a + c] = dFb[2 * a] * Fn2[c] + dFb[2 * a + 1] * Fn2[2 + c];
    for (int k = 0; k < 9; k++) F13[k] = 0.0;
    F13[0] = F12[0]; F13[1] = F12[1]; F13[3] = F12[2]; F13[4] = F12[3]; F13[8] = 1.0;
    double d2[3] = {0, 0, 0}, R2[4], J2 = 0.0, d3[6] = {0, 0, 0, 0, 0, 0}, R3[9], J3 = 0.0;
    const int s2 = ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn2, F12, d2, R2, &J2);
    const int s3 = ladruno_hypo::hypoIncrementFromF(Fn3, F13, d3, R3, &J3);
    check(s2 == 0 && s3 == 0 && d3[2] == 0.0 && d3[4] == 0.0 && d3[5] == 0.0 &&
          std::fabs(d2[0] - d3[0]) <= 1e-15 && std::fabs(d2[1] - d3[1]) <= 1e-15 &&
          std::fabs(d2[2] - d3[3]) <= 1e-15, "L1c: plane strain = the 3D entry with F_zz = 1");
  }
  for (int t = 0; t < 200; t++) {
    double dF[9], dFc[9], dFb[9];
    for (int k = 0; k < 9; k++) {
      dF[k] = (k % 4 == 0 ? 1.0 : 0.0) + 0.03 * std::sin(2.1 * t + k);
      dFc[k] = (k % 4 == 0 ? 1.0 : 0.0) + 0.03 * std::cos(1.7 * t + 0.5 * k);
    }
    if (fbarIncrement(3, dF, dFc, dFb) != 0) { check(false, "J1: refused a valid 3D pair"); continue; }
    check(std::fabs(detN(3, dFb) - detN(3, dFc)) <= 1e-14 * detN(3, dFc), "J1: det(dFbar) = det(dFc), 3D");
  }
  {
    double dF[4] = {1, 0, 0, 1}, dFc[4] = {-1, 0, 0, 1}, dFb[4] = {7, 7, 7, 7};
    check(fbarIncrement(2, dF, dFc, dFb) == -1 && dFb[0] == 7.0, "J1: det(dFc) <= 0 refused, output untouched");
    double dG[4] = {0, 1, 1, 0};
    check(fbarIncrement(2, dG, dF, dFb) == -1 && dFb[0] == 7.0, "J1: det(dF) <= 0 refused, output untouched");
    double same[4] = {1.01, 0.002, -0.003, 0.995};
    check(fbarIncrement(2, same, same, dFb) == 0 && dFb[0] == same[0] && dFb[1] == same[1] &&
          dFb[2] == same[2] && dFb[3] == same[3], "J1: an equal pair gives dFbar = dF bit for bit");
  }
  // a non-affine field changes the increment (the averaging is not the identity)
  {
    double dF[4] = {1.01, 0.0, 0.0, 1.0}, dFc[4] = {1.0, 0.0, 0.0, 1.0}, dFb[4];
    fbarIncrement(2, dF, dFc, dFb);
    check(std::fabs(dFb[0] - dF[0]) > 1e-4, "J1: unequal volume changes are modified");
  }

  std::printf("[checks] %d checks, %d failed\n", g_checks, g_fail);
  std::printf("RESULT: %s\n", g_fail == 0 ? "PASS" : "FAIL");
  return g_fail == 0 ? 0 : 1;
}
