// WP-261 harness: the LadrunoIMKBeam safeguarded local hinge solve
// (SRC/element/ladrunoIMKBeam/LadrunoIMKSafeSolve.h), without OpenSees.
//
// Two series hinges on an elastic interior (a = 4EI/L, b = 2EI/L), each hinge a
// Steel01 clip evaluated from its COMMITTED state (as Steel01 does), solved for
// a trial chord deformation v. Prints one "key value" line per check; the
// python driver (tests/test_ladrunoIMK_safe_solve_cpp.py) asserts on them.
//
//   plain_*   the plain local Newton of ladrunoIMKSolveAxis (50 iterations,
//             residual damping), transcribed: it must FAIL on the unloading case
//   safe_*    the safeguarded solve: converges, matches region enumeration
//   noroot_*  a law with an upward moment jump has no root: the solve must
//             return false, with its last evaluation at the returned point

#include "LadrunoIMKSafeSolve.h"
#include <stdio.h>
#include <float.h>

struct Clip {
  double k0, kp, cb, thn, Mn, ct;
  void eval(double th, double &M, double &kt) const {
    const double d = th - thn;
    if (fabs(d) <= DBL_EPSILON) { M = Mn; kt = ct; return; }
    const double c = Mn + k0 * d;
    double T = kp * th + cb; if (c < T) T = c;
    const double lo = kp * th - cb; if (T < lo) T = lo;
    M = T; kt = (fabs(T - c) < DBL_EPSILON) ? k0 : kp;
  }
};

struct Pair {
  double a, b, v[2];
  Clip h[2];
  int nEval;
  double kRes() const { return a + b + h[0].k0 + h[1].k0; }
  double lastX[2];
  // residual G = q - M, H = -dG/dth, tight test of ladrunoIMKSolveAxis
  bool operator()(const double x[6], double G[6], double H[6][6]) {
    nEval++;
    lastX[0] = x[0]; lastX[1] = x[1];
    double M[2], kt[2];
    for (int k = 0; k < 2; k++) h[k].eval(x[k], M[k], kt[k]);
    const double qi = a * (v[0] - x[0]) + b * (v[1] - x[1]);
    const double qj = b * (v[0] - x[0]) + a * (v[1] - x[1]);
    G[0] = qi - M[0]; G[1] = qj - M[1];
    for (int r = 0; r < 6; r++) for (int c = 0; c < 6; c++) H[r][c] = 0.0;
    H[0][0] = a + kt[0]; H[1][1] = a + kt[1]; H[0][1] = H[1][0] = b;
    const double gn = sqrt(G[0] * G[0] + G[1] * G[1]);
    const double rf = fabs(qi) + fabs(qj) + fabs(M[0]) + fabs(M[1]);
    const double sc = (a + b) * (fabs(v[0]) + fabs(v[1]) + fabs(x[0]) + fabs(x[1]))
                    + fabs(kt[0] * x[0]) + fabs(kt[1] * x[1]);
    return gn <= 1.0e-10 * rf + 1.0e-12 + 1.0e-12 * sc + kRes() * 1.0e-14;
  }
};

// the plain loop of ladrunoIMKSolveAxis (hi && hj, no shear), transcribed
static bool plainNewton(Pair &p, double th[2])
{
  double damp = 1.0, gnPrev = 1.0e300;
  const double ktFloor = 1.0e-8 * p.a;
  for (int it = 0; it < 50; it++) {
    double M[2], kt[2];
    for (int k = 0; k < 2; k++) p.h[k].eval(th[k], M[k], kt[k]);
    const double qi = p.a * (p.v[0] - th[0]) + p.b * (p.v[1] - th[1]);
    const double qj = p.b * (p.v[0] - th[0]) + p.a * (p.v[1] - th[1]);
    const double gi = qi - M[0], gj = qj - M[1];
    const double gn = sqrt(gi * gi + gj * gj);
    if (gn <= 1.0e-10 * (fabs(qi) + fabs(qj) + fabs(M[0]) + fabs(M[1])) + 1.0e-12) return true;
    if (gn > gnPrev) damp *= 0.5; else damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
    gnPrev = gn;
    double ki = (fabs(kt[0]) < ktFloor) ? ktFloor : kt[0];
    double kj = (fabs(kt[1]) < ktFloor) ? ktFloor : kt[1];
    const double A00 = p.a + ki, A01 = p.b, A11 = p.a + kj;
    const double det = A00 * A11 - A01 * A01;
    th[0] += damp * ( A11 * gi - A01 * gj) / det;
    th[1] += damp * (-A01 * gi + A00 * gj) / det;
  }
  return false;
}

// exact solution by enumerating the 3x3 Steel01 regions (unique by monotonicity)
static bool exactPair(const Pair &p, double th[2])
{
  for (int r0 = -1; r0 <= 1; r0++)
    for (int r1 = -1; r1 <= 1; r1++) {
      const int r[2] = { r0, r1 };
      double kk[2], m0[2];
      for (int k = 0; k < 2; k++) {
        if (r[k] == 0) { kk[k] = p.h[k].k0; m0[k] = p.h[k].Mn - p.h[k].k0 * p.h[k].thn; }
        else           { kk[k] = p.h[k].kp; m0[k] = r[k] * p.h[k].cb; }
      }
      // (K + diag kk) th = K v - m0
      const double A00 = p.a + kk[0], A01 = p.b, A11 = p.a + kk[1];
      const double r0_ = p.a * p.v[0] + p.b * p.v[1] - m0[0];
      const double r1_ = p.b * p.v[0] + p.a * p.v[1] - m0[1];
      const double det = A00 * A11 - A01 * A01;
      const double t0 = (A11 * r0_ - A01 * r1_) / det;
      const double t1 = (-A01 * r0_ + A00 * r1_) / det;
      double M0, M1, k_;
      p.h[0].eval(t0, M0, k_); p.h[1].eval(t1, M1, k_);
      const double qi = p.a * (p.v[0] - t0) + p.b * (p.v[1] - t1);
      const double qj = p.b * (p.v[0] - t0) + p.a * (p.v[1] - t1);
      if (fabs(qi - M0) + fabs(qj - M1) <= 1.0e-9 * (fabs(qi) + fabs(qj) + fabs(M0) + fabs(M1) + 1.0)) {
        th[0] = t0; th[1] = t1; return true;
      }
    }
  return false;
}

// committed Steel01 state of a hinge pair after monotonic loading to vn
static void commitMonotonic(Pair &p, double k0, double My, double bb, const double vn[2])
{
  for (int k = 0; k < 2; k++) {
    p.h[k].k0 = k0; p.h[k].kp = bb * k0; p.h[k].cb = My * (1.0 - bb);
    p.h[k].thn = 0.0; p.h[k].Mn = 0.0; p.h[k].ct = k0;
  }
  p.v[0] = vn[0]; p.v[1] = vn[1];
  double th[2];
  exactPair(p, th);
  for (int k = 0; k < 2; k++) {
    double M, kt;
    p.h[k].eval(th[k], M, kt);
    p.h[k].thn = th[k]; p.h[k].Mn = M; p.h[k].ct = kt;
  }
}

// a law with an upward moment jump J at th = thStar: g = q - M jumps from
// positive to negative there and has no root when the jump straddles q
struct JumpPair {
  double a, b, v[2], k, thStar, J;
  double lastX[2];
  bool operator()(const double x[6], double G[6], double H[6][6]) {
    lastX[0] = x[0]; lastX[1] = x[1];
    double M[2], kt[2];
    for (int i = 0; i < 2; i++) {
      M[i] = k * x[i] + ((x[i] >= thStar) ? J : 0.0);
      kt[i] = k;
    }
    const double qi = a * (v[0] - x[0]) + b * (v[1] - x[1]);
    const double qj = b * (v[0] - x[0]) + a * (v[1] - x[1]);
    G[0] = qi - M[0]; G[1] = qj - M[1];
    for (int r = 0; r < 6; r++) for (int c = 0; c < 6; c++) H[r][c] = 0.0;
    H[0][0] = a + kt[0]; H[1][1] = a + kt[1]; H[0][1] = H[1][0] = b;
    const double gn = sqrt(G[0] * G[0] + G[1] * G[1]);
    return gn <= 1.0e-10 * (fabs(qi) + fabs(qj) + fabs(M[0]) + fabs(M[1])) + 1.0e-12;
  }
};

int main()
{
  // element of tests/test_ladrunoIMKBeam2d_unload.py: E 2e11, Iz 2e-5, L 3
  const double E = 2.0e11, Iz = 2.0e-5, L = 3.0;
  const double a = 4.0 * E * Iz / L, b = 2.0 * E * Iz / L;
  const double k0 = 100.0 * 6.0 * E * Iz / L, My = 1.0e5, bb = 1.0e-4;

  // --- 1. the deterministic unloading 2-cycle --------------------------------
  {
    Pair p; p.a = a; p.b = b;
    const double vn[2] = { 0.0, 0.02 };
    commitMonotonic(p, k0, My, bb, vn);
    p.v[0] = 0.0; p.v[1] = 0.018;                 // unload by 0.002 rad
    double thP[2] = { p.h[0].thn, p.h[1].thn };
    const bool okPlain = plainNewton(p, thP);
    double x[6] = { p.h[0].thn, p.h[1].thn, 0, 0, 0, 0 };
    p.nEval = 0;
    const bool okSafe = ladrunoIMKSafeguardedSolve(2, x, p, 1.0e-8 * a);
    double thE[2];
    const bool okE = exactPair(p, thE);
    const double thy = My / k0;
    const double err = (fabs(x[0] - thE[0]) + fabs(x[1] - thE[1])) / thy;
    printf("plain_unload_converged %d\n", okPlain ? 1 : 0);
    printf("safe_unload_converged %d\n", okSafe ? 1 : 0);
    printf("safe_unload_exact_found %d\n", okE ? 1 : 0);
    printf("safe_unload_err %.3e\n", err);
    printf("safe_unload_evals %d\n", p.nEval);
    printf("safe_unload_last_at_x %d\n", (p.lastX[0] == x[0] && p.lastX[1] == x[1]) ? 1 : 0);
  }

  // --- 2. a sweep of loading / unloading steps from plastic states -----------
  {
    int nPlainFail = 0, nSafeFail = 0, nCases = 0, maxEval = 0;
    double maxErr = 0.0;
    unsigned s = 12345u;
    for (int n = 0; n < 4000; n++) {
      s = s * 1664525u + 1013904223u; const double u1 = (s >> 8) / 16777216.0;
      s = s * 1664525u + 1013904223u; const double u2 = (s >> 8) / 16777216.0;
      s = s * 1664525u + 1013904223u; const double u3 = (s >> 8) / 16777216.0;
      s = s * 1664525u + 1013904223u; const double u4 = (s >> 8) / 16777216.0;
      Pair p; p.a = a; p.b = b;
      const double vn[2] = { 0.03 * (2.0 * u1 - 1.0), 0.03 * (2.0 * u2 - 1.0) };
      commitMonotonic(p, k0, My, bb, vn);
      p.v[0] = vn[0] + 0.01 * (2.0 * u3 - 1.0);
      p.v[1] = vn[1] + 0.01 * (2.0 * u4 - 1.0);
      double thP[2] = { p.h[0].thn, p.h[1].thn };
      if (!plainNewton(p, thP)) nPlainFail++;
      double x[6] = { p.h[0].thn, p.h[1].thn, 0, 0, 0, 0 };
      p.nEval = 0;
      const bool ok = ladrunoIMKSafeguardedSolve(2, x, p, 1.0e-8 * a);
      if (!ok) { nSafeFail++; continue; }
      if (p.nEval > maxEval) maxEval = p.nEval;
      double thE[2];
      if (exactPair(p, thE)) {
        const double e = (fabs(x[0] - thE[0]) + fabs(x[1] - thE[1])) / (My / k0);
        if (e > maxErr) maxErr = e;
      }
      nCases++;
    }
    printf("sweep_plain_fail %d\n", nPlainFail);
    printf("sweep_safe_fail %d\n", nSafeFail);
    printf("sweep_safe_max_err %.3e\n", maxErr);
    printf("sweep_safe_max_evals %d\n", maxEval);
    printf("sweep_cases %d\n", nCases);
  }

  // --- 2b. Steel01's absolute DBL_EPSILON dead zone ------------------------
  // a virgin hinge pair at v ~ 50 eps: M(th) = 0 for |th| <= eps, then k0 th,
  // so the hinge equation has no root. The residual left (~k0 eps) is the
  // law's resolution, not a failure: the solve must accept it.
  {
    Pair p; p.a = a; p.b = b;
    const double vn[2] = { 0.0, 0.0 };
    commitMonotonic(p, k0, My, bb, vn);
    p.v[0] = 0.0; p.v[1] = 50.0 * DBL_EPSILON;
    double x[6] = { 0, 0, 0, 0, 0, 0 };
    p.nEval = 0;
    const bool ok = ladrunoIMKSafeguardedSolve(2, x, p, 1.0e-8 * a);
    printf("deadzone_converged %d\n", ok ? 1 : 0);
    printf("deadzone_evals %d\n", p.nEval);
  }

  // --- 3. no root: the solve reports failure, consistently ------------------
  {
    JumpPair p; p.a = a; p.b = b; p.k = 0.5 * a; p.thStar = 1.0e-3; p.J = 1.0e7;
    // with x0 = x1 = thStar and v symmetric, q - k*thStar lies inside (0, J)
    p.v[0] = p.v[1] = p.thStar + (p.k * p.thStar + 0.5 * p.J) / (a + b);
    double x[6] = { 0, 0, 0, 0, 0, 0 };
    const bool ok = ladrunoIMKSafeguardedSolve(2, x, p, 1.0e-8 * a);
    printf("noroot_converged %d\n", ok ? 1 : 0);
    printf("noroot_last_at_x %d\n", (p.lastX[0] == x[0] && p.lastX[1] == x[1]) ? 1 : 0);
    printf("noroot_x_finite %d\n", (x[0] == x[0] && x[1] == x[1]) ? 1 : 0);
  }
  return 0;
}
