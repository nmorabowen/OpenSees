// mpm_bspline_core_check — standalone g++ check of the WP-232i additions to the
// explicit MPM core (SRC/domain/pattern/ladrunoMPM/LadrunoMPMCore.h) against the
// frozen tests/mpm_bspline_interface.md and the "232i" gate table. No OpenSees.
//
//   g++ -O2 -std=c++17 -Wall -Wextra -I SRC/element/solidTransformation
//       -I SRC/domain/pattern/ladrunoMPM tests/mpm_bspline_core_check.cpp -o bsc
//   ./bsc            (or ./bsc mms  : also the mode-F MMS mini-ladder, INFO only)
//
// G1     basis properties per B-spline degree, 2D and 3D, >= 1e4 random admissible
//        positions including every face's 2h band; cardinal values; C1 across knots.
// G1-K4  the folded per-component S' (B5, B8) within 2h of each `fix` face.
// G8     images against the full domain (mode F, bspline2): a body symmetric about an
//        interior grid plane on the full grid, against its half on a grid whose face
//        is that plane with a roller (fix + images), at small and large amplitude.
// B11    a uniform particle body force equals the same gravity, bit for bit.
// The material is the interface §9 hypoelastic law (the ElasticIsotropic law).
// Exit 0 + "RESULT: PASS" on success.

#include <LadrunoMPMCore.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace ladruno_mpm;

static const double kPi = 3.14159265358979323846;

static int g_fail = 0, g_checks = 0;

static void check(bool ok, const std::string &what, double got, double lim)
{
  g_checks++;
  if (!ok) {
    g_fail++;
    std::printf("FAIL %s  got %.3e  limit %.3e\n", what.c_str(), got, lim);
  }
}

struct Hypo : public MaterialPoints {
  int nd = 2;
  double lam = 0.0, mu = 0.0;
  int setTrialStrain(int, const double *e, double *s) override
  {
    if (nd == 3) {
      const double tr = e[0] + e[1] + e[2];
      for (int k = 0; k < 3; k++) s[k] = lam * tr + 2.0 * mu * e[k];
      for (int k = 3; k < 6; k++) s[k] = mu * e[k];
    } else {
      const double tr = e[0] + e[1];
      s[0] = lam * tr + 2.0 * mu * e[0];
      s[1] = lam * tr + 2.0 * mu * e[1];
      s[2] = mu * e[2];
    }
    return 0;
  }
  int commitState(int) override { return 0; }
  int revertToLastCommit(int) override { return 0; }
  int pendingRefusals() override { return 0; }
  double normalModulus(int) override { return lam + 2.0 * mu; }
};

static Grid makeGrid(int nd, double h, const int *n, const double *x0)
{
  Grid g;
  g.ndm = nd;
  g.h = h;
  for (int d = 0; d < 3; d++) { g.x0[d] = d < nd ? x0[d] : 0.0; g.n[d] = d < nd ? n[d] : 0; }
  return g;
}

static void nodeCoord(const Grid &g, long ext, double *c)
{
  int ijk[3];
  g.extIJK(ext, ijk);
  for (int d = 0; d < g.ndm; d++) c[d] = g.x0[d] + ijk[d] * g.h;
}

// random admissible positions: uniform, plus every face's 2h band (and the faces)
static std::vector<std::vector<double>> positions(const Grid &g, int nUniform, int nBand,
                                                  std::mt19937_64 &rng)
{
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<std::vector<double>> out;
  const int nd = g.ndm;
  auto rnd = [&](std::vector<double> &x) {
    for (int d = 0; d < nd; d++) x[d] = g.x0[d] + U(rng) * g.n[d] * g.h;
  };
  for (int k = 0; k < nUniform; k++) { std::vector<double> x(nd); rnd(x); out.push_back(x); }
  for (int f = 0; f < 2 * nd; f++) {
    const int d = f / 2;
    const double lo = g.x0[d], hi = g.x0[d] + g.n[d] * g.h;
    for (int k = 0; k < nBand; k++) {
      std::vector<double> x(nd);
      rnd(x);
      const double s = (k < 8) ? 0.0 : U(rng) * 2.0 * g.h;   // the face itself, then the band
      x[d] = (f % 2 == 0) ? lo + s : hi - s;
      out.push_back(x);
    }
  }
  return out;
}

// ---- G1 -------------------------------------------------------------------------
static void g1(Basis b, int nd)
{
  const int n[3] = {6, 5, 5};
  const double x0[3] = {0.3, -0.2, 0.1}, h = 0.25;
  Grid g = makeGrid(nd, h, n, x0);
  std::mt19937_64 rng(232 + 10 * nd + (int)b);
  auto X = positions(g, 10000, 1500, rng);
  const std::string tag = std::string(basisName(b)) + " " + std::to_string(nd) + "D";
  double wS = 0, wPU = 0, wG = 0, wL = 0, wLG = 0;
  int count = 0;
  for (auto &x : X) {
    double lp[3] = {0.1, 0.1, 0.1};
    for (int d = 0; d < nd; d++)
      if (!admissibleAxis(g, b, d, x[d], lp[d])) { check(false, tag + " admissible", 1, 0); }
    std::vector<long> node;
    std::vector<double> S, dS;
    particleSupport(g, b, x.data(), lp, node, S, dS);
    double sS = 0, sG[3] = {0, 0, 0}, sx[3] = {0, 0, 0}, xg[3][3] = {{0}};
    for (size_t k = 0; k < node.size(); k++) {
      double c[3];
      nodeCoord(g, node[k], c);
      wS = std::min(wS, S[k]);
      sS += S[k];
      for (int d = 0; d < nd; d++) {
        sG[d] += dS[nd * k + d];
        sx[d] += S[k] * c[d];
        for (int e = 0; e < nd; e++) xg[d][e] += c[d] * dS[nd * k + e];
      }
    }
    wPU = std::max(wPU, std::fabs(sS - 1.0));
    for (int d = 0; d < nd; d++) {
      wG = std::max(wG, std::fabs(sG[d]));
      wL = std::max(wL, std::fabs(sx[d] - x[d]));
      for (int e = 0; e < nd; e++) wLG = std::max(wLG, std::fabs(xg[d][e] - (d == e ? 1.0 : 0.0)));
    }
    count++;
  }
  check(count >= 10000, tag + " positions >= 1e4", count, 10000);
  check(wS >= 0.0, tag + " S >= 0", wS, 0.0);
  check(wPU <= 1e-13, tag + " |sum S - 1|", wPU, 1e-13);
  check(wG <= 1e-11 / h, tag + " |sum grad S|", wG, 1e-11 / h);
  check(wL <= 1e-12 * h, tag + " |sum S x - x_p|", wL, 1e-12 * h);
  check(wLG <= 1e-11, tag + " |sum x (x) grad S - I|", wLG, 1e-11);
  std::printf("INFO G1 %s: %d positions; min S %.2e, PU %.2e, sum dS %.2e, lin %.2e, lin-grad %.2e\n",
              tag.c_str(), count, wS, wPU, wG, wL, wLG);
}

// ---- G1: cardinal values and C1 continuity ----------------------------------------
static void g1Cardinal()
{
  const double h = 0.37;
  struct CV { int deg; double r; double S; };
  const CV cv[] = {{2, 0.0, 0.75}, {2, 0.5, 0.5}, {2, -0.5, 0.5}, {2, 1.0, 0.125}, {2, -1.0, 0.125},
                   {2, 1.5, 0.0}, {2, -1.5, 0.0}, {3, 0.0, 2.0 / 3.0}, {3, 1.0, 1.0 / 6.0},
                   {3, -1.0, 1.0 / 6.0}, {3, 2.0, 0.0}, {3, -2.0, 0.0}};
  for (const CV &c : cv) {
    double S, dS;
    bsplineFactor(c.deg, c.r * h, h, S, dS);
    char buf[96];
    std::snprintf(buf, sizeof buf, "cardinal bspline%d S(%gh)", c.deg, c.r);
    check(std::fabs(S - c.S) <= 1e-14, buf, std::fabs(S - c.S), 1e-14);
  }
  for (int deg = 2; deg <= 3; deg++) {
    // knots: quadratic at +-h/2, +-3h/2; cubic at 0, +-h, +-2h
    std::vector<double> knots = (deg == 2) ? std::vector<double>{-1.5, -0.5, 0.5, 1.5}
                                           : std::vector<double>{-2.0, -1.0, 0.0, 1.0, 2.0};
    double jump = 0.0, fd = 0.0;
    const double eps = 1e-10 * h;
    for (double k : knots) {
      double Sa, da, Sb, db;
      bsplineFactor(deg, k * h - eps, h, Sa, da);
      bsplineFactor(deg, k * h + eps, h, Sb, db);
      jump = std::max(jump, std::fabs(db - da));
    }
    // dS against a centred difference of S over a sweep of the whole support
    for (int i = -4000; i <= 4000; i++) {
      const double r = i * (deg + 1) * h / 8000.0 + 1e-7 * h, e = 1e-6 * h;
      double S, dS, Sp, dp, Sm, dm;
      bsplineFactor(deg, r, h, S, dS);
      bsplineFactor(deg, r + e, h, Sp, dp);
      bsplineFactor(deg, r - e, h, Sm, dm);
      fd = std::max(fd, std::fabs(dS - (Sp - Sm) / (2.0 * e)));
    }
    char buf[96];
    std::snprintf(buf, sizeof buf, "C1 bspline%d max dS jump at knots", deg);
    check(jump <= 1e-8 / h, buf, jump, 1e-8 / h);
    std::snprintf(buf, sizeof buf, "bspline%d dS vs centred difference of S", deg);
    check(fd <= 1e-6 / h, buf, fd, 1e-6 / h);
  }
}

// ---- G1-K4: the folded S' near each `fix` face -------------------------------------
static void g1k4(Basis b, int nd, bool clamped)
{
  const int n[3] = {6, 5, 5};
  const double x0[3] = {-0.4, 0.15, 0.05}, h = 0.2;
  Grid g = makeGrid(nd, h, n, x0);
  FaceConditions fc;
  for (int d = 0; d < nd; d++) {
    const unsigned m = clamped ? ((1u << nd) - 1u) : (1u << d);   // rollers, or clamped faces
    fc.fixMask[2 * d] = fc.fixMask[2 * d + 1] = m;
  }
  std::mt19937_64 rng(999 + nd + (int)b + (clamped ? 7 : 0));
  auto X = positions(g, 0, 3000, rng);
  const std::string tag = std::string("K4 ") + basisName(b) + " " + std::to_string(nd) + "D " +
                          (clamped ? "clamped" : "rollers");
  double wNeg = 0, wOdd = 0, wPU = 0, wT = 0;
  int nOdd = 0, nEven = 0;
  for (auto &x : X) {
    double lp[3] = {0.05, 0.05, 0.05};
    std::vector<long> node, on;
    std::vector<double> S, dS, W, dW, Sc, dSc;
    particleSupport(g, b, x.data(), lp, node, S, dS);
    foldedSupport(g, b, fc, x.data(), lp, on, W, dW, Sc, dSc);
    // which faces fold this particle's support (axis d: -1 min, +1 max, 0 none)
    int side[3] = {0, 0, 0};
    for (long k : node) {
      int ijk[3];
      g.extIJK(k, ijk);
      for (int d = 0; d < nd; d++) {
        if (ijk[d] < 0) side[d] = -1;
        if (ijk[d] > g.n[d]) side[d] = 1;
      }
    }
    for (int a = 0; a < nd; a++) {
      // faces about which component a is odd (it is fixed there) / even
      int oddAxes = 0, oddAxis = -1;
      for (int d = 0; d < nd; d++) {
        if (side[d] == 0) continue;
        const unsigned m = fc.fixMask[side[d] < 0 ? 2 * d : 2 * d + 1];
        if (m & (1u << a)) { oddAxes++; oddAxis = d; }
      }
      double sS = 0, sOdd = 0, sT[3] = {0, 0, 0};
      for (size_t j = 0; j < on.size(); j++) {
        double c[3];
        nodeCoord(g, on[j], c);
        const double s = Sc[nd * j + a];
        wNeg = std::min(wNeg, s);
        sS += s;
        for (int d = 0; d < nd; d++) sT[d] += s * c[d];
        if (oddAxes == 1) {
          const double xf = side[oddAxis] < 0 ? g.x0[oddAxis] : g.x0[oddAxis] + g.n[oddAxis] * h;
          sOdd += s * (c[oddAxis] - xf);
        }
      }
      if (oddAxes == 1) {
        const double xf = side[oddAxis] < 0 ? g.x0[oddAxis] : g.x0[oddAxis] + g.n[oddAxis] * h;
        wOdd = std::max(wOdd, std::fabs(sOdd - (x[oddAxis] - xf)));
        nOdd++;
      } else if (oddAxes == 0) {
        wPU = std::max(wPU, std::fabs(sS - 1.0));
        for (int d = 0; d < nd; d++)       // tangential: an axis not folded for this particle
          if (side[d] == 0) wT = std::max(wT, std::fabs(sT[d] - x[d]));
        nEven++;
      }
    }
  }
  check(wNeg >= 0.0, tag + " S' >= 0", wNeg, 0.0);
  check(wOdd <= 1e-12 * h, tag + " odd: sum S' (x - x_face)_n", wOdd, 1e-12 * h);
  check(clamped || nEven > 0, tag + " even components sampled", nEven, 1);
  check(wPU <= 1e-13, tag + " even: |sum S' - 1|", wPU, 1e-13);
  check(wT <= 1e-12 * h, tag + " even: sum S' x_t", wT, 1e-12 * h);
  std::printf("INFO G1-K4 %s: %d odd / %d even component samples; min S' %.2e, odd %.2e, PU %.2e, "
              "tang %.2e\n", tag.c_str(), nOdd, nEven, wNeg, wOdd, wPU, wT);
}

// ---- G8: images against the full domain -------------------------------------------
// Plane x = 0. Full grid [-L, L] x [0, H]; half grid [0, L] x [0, H] with `fix xmin x`.
// Body [-a, a] x [y0, y1] (lattice, ppc 2), initial velocity odd in x for v_x and even
// in x for v_y: v_x = V sin(pi x / a) (1 + y), v_y = 0.5 V cos(pi x / (2a)) sin(pi y).
static double g8(double V, int steps, bool verbose, double masscut = 1e-10)
{
  const int nd = 2, ppc = 2;
  const double h = 0.125, L = 1.0, H = 1.0, a = 0.5, y0 = 0.25, y1 = 0.75;
  const double E = 1.0e3, nu = 0.3, rho = 1.0;
  Hypo M;
  M.nd = 2;
  M.lam = E * nu / ((1 + nu) * (1 - 2 * nu));
  M.mu = E / (2 * (1 + nu));
  auto setup = [&](Engine &e, bool half) {
    std::string err;
    const double x0[2] = {half ? 0.0 : -L, 0.0};
    const int n[2] = {half ? (int)std::lround(L / h) : (int)std::lround(2 * L / h), (int)std::lround(H / h)};
    e.setGrid(nd, x0, h, n, err);
    e.par.basis = BASIS_BSPLINE2;
    e.par.masscut = masscut;           // G8-F-cut: 0.05
    if (half) e.fc.fixMask[0] = 1u;            // fix xmin x: the roller, images on
    // particles: the full lattice's x >= 0 half, coordinates mirrored exactly (-x)
    int64_t id = 0;
    const double hs = h / ppc;
    for (int iy = 0; iy < (int)std::lround((y1 - y0) / hs); iy++)
      for (int ix = 0; ix < (int)std::lround(a / hs); ix++) {
        const double xr = ((ix + 0.5) * h) / ppc, y = y0 + ((iy + 0.5) * h) / ppc;
        const double lp[2] = {h / (2.0 * ppc), h / (2.0 * ppc)};
        for (int side = (half ? 1 : 0); side < 2; side++) {
          const double x = side ? xr : -xr;
          const double pos[2] = {x, y};
          const double v[2] = {V * std::sin(kPi * x / a) * (1.0 + y),
                               0.5 * V * std::cos(kPi * x / (2 * a)) * std::sin(kPi * y)};
          e.addParticle(id++, 0, 1, rho, hs * hs, pos, lp, v, nullptr, err);
        }
      }
  };
  Engine full, half;
  setup(full, false);
  setup(half, true);
  const double dt = 0.1 * h / std::sqrt((M.lam + 2 * M.mu) / rho);
  std::string err;
  for (int s = 0; s < steps; s++) {
    if (full.advance(dt, 1.0, M, err) != 0 || half.advance(dt, 1.0, M, err) != 0) {
      check(false, "G8 run: " + err, 1, 0);
      return 1.0;
    }
  }
  // half particle k (ids 0, 1, ...) is full particle 2k+1 (the x > 0 side)
  double dx = 0, xs = 0, ds = 0, ss = 0, maxStrain = 0;
  for (int k = 0; k < half.P.size(); k++) {
    const int q = 2 * k + 1;
    for (int d = 0; d < nd; d++) {
      dx = std::max(dx, std::fabs(half.P.x[nd * k + d] - full.P.x[nd * q + d]));
      xs = std::max(xs, std::fabs(full.P.x[nd * q + d]));
      maxStrain = std::max(maxStrain, std::fabs(full.P.F[4 * q + 3 * d] - 1.0));
    }
    for (int c = 0; c < 3; c++) {
      ds = std::max(ds, std::fabs(half.P.sig[3 * k + c] - full.P.sig[3 * q + c]));
      ss = std::max(ss, std::fabs(full.P.sig[3 * q + c]));
    }
  }
  const double rx = dx / xs, rs = ds / ss;
  if (verbose)
    std::printf("INFO G8 V=%g, %d steps: max|F_ii - 1| %.3e, rel dx %.2e, rel dsig %.2e\n", V, steps,
                maxStrain, rx, rs);
  const std::string tag = "(V=" + std::to_string(V) + ", masscut=" + std::to_string(masscut) + ")";
  check(rx <= 1e-10, "G8 positions " + tag, rx, 1e-10);
  check(rs <= 1e-10, "G8 stresses " + tag, rs, 1e-10);
  return maxStrain;
}

// ---- B11: a uniform particle body force equals the same gravity ----------------------
static void bodyForceEqualsGravity()
{
  Hypo M;
  M.lam = 100; M.mu = 80;
  auto run = [&](bool useBF) {
    Engine e;
    std::string err;
    const double x0[2] = {0, 0};
    const int n[2] = {4, 4};
    e.setGrid(2, x0, 0.25, n, err);
    e.par.basis = BASIS_BSPLINE2;
    e.fc.fixMask[2] = 2u;
    const double geom[4] = {0.0, 0.0, 1.0, 0.5};
    int64_t id = 0;
    e.seedLattice(0, geom, 2, 0, 1, 2.0, id, err);
    const double g[2] = {0.3, -9.81};
    if (useBF) {
      for (int p = 0; p < e.P.size(); p++) { e.bodyForce.push_back(g[0]); e.bodyForce.push_back(g[1]); }
    } else {
      e.par.gravity[0] = g[0]; e.par.gravity[1] = g[1];
    }
    for (int s = 0; s < 50; s++) e.advance(1e-3, 1.0, M, err);
    return e.P.x;
  };
  const auto a = run(false), b = run(true);
  double d = 0;
  for (size_t k = 0; k < a.size(); k++) d = std::max(d, std::fabs(a[k] - b[k]));
  check(d == 0.0, "B11 body force == gravity (bitwise)", d, 0.0);
  // a stale length is refused
  Engine e;
  std::string err;
  const double x0[2] = {0, 0};
  const int n[2] = {2, 2};
  e.setGrid(2, x0, 0.5, n, err);
  const double geom[4] = {0.0, 0.0, 1.0, 1.0};
  int64_t id = 0;
  e.seedLattice(0, geom, 1, 0, 1, 1.0, id, err);
  e.bodyForce.assign(3, 0.0);
  check(e.advance(1e-3, 1.0, M, err) == -1, "B11 stale body-force length refused", 0, 0);
}

// ---- the mode-F MMS mini-ladder (INFO only; the gate runs through the binaries) -----
static double mmsRun(Basis basis, int nCells, double Fmax, double dt, double T)
{
  const double rho = 1.0, E = 1.0, nu = 0.3, c = 2.0, w = c * kPi;
  const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
  const double A = Fmax / (2 * kPi);
  Hypo M;
  M.lam = lam; M.mu = mu;
  Engine e;
  std::string err;
  const double x0[2] = {0, 0};
  const int n[2] = {nCells, nCells};
  const double h = 1.0 / nCells;
  e.setGrid(2, x0, h, n, err);
  e.par.basis = basis;
  e.fc.fixMask[0] = e.fc.fixMask[1] = 1u;
  e.fc.fixMask[2] = e.fc.fixMask[3] = 2u;
  const double geom[4] = {0, 0, 1, 1};
  int64_t id = 0;
  e.seedLattice(0, geom, 2, 0, 1, rho, id, err);
  const int N = e.P.size();
  std::vector<double> X(e.P.x);
  for (int p = 0; p < N; p++)
    for (int i = 0; i < 2; i++) e.P.v[2 * p + i] = w * A * std::sin(2 * kPi * X[2 * p + i]);
  e.bodyForce.assign(2 * N, 0.0);
  const int steps = (int)std::lround(T / dt);
  double t = 0.0;
  for (int s = 0; s < steps; s++) {
    for (int p = 0; p < N; p++) {
      double F[2], g1[2];
      for (int i = 0; i < 2; i++) {
        const double Xi = X[2 * p + i];
        F[i] = 1 + 2 * kPi * A * std::cos(2 * kPi * Xi) * std::sin(w * t);
        g1[i] = -4 * kPi * kPi * A * std::sin(2 * kPi * Xi) * std::sin(w * t);
      }
      for (int i = 0; i < 2; i++) {
        const double ai = -w * w * A * std::sin(2 * kPi * X[2 * p + i]) * std::sin(w * t);
        const double divP = F[1 - i] * (lam + 2 * mu) * g1[i] / F[i];
        e.bodyForce[2 * p + i] = ai - divP / rho;
      }
    }
    if (e.advance(dt, 1.0, M, err) != 0) { std::printf("INFO MMS run failed: %s\n", err.c_str()); return -1; }
    t += dt;
  }
  double s2 = 0;
  for (int p = 0; p < N; p++)
    for (int i = 0; i < 2; i++) {
      const double ue = A * std::sin(2 * kPi * X[2 * p + i]) * std::sin(w * t);
      const double d = e.P.x[2 * p + i] - X[2 * p + i] - ue;
      s2 += d * d;
    }
  return std::sqrt(s2 / N) / A;
}

int main(int argc, char **argv)
{
  for (int nd = 2; nd <= 3; nd++) {
    g1(BASIS_BSPLINE2, nd);
    g1(BASIS_BSPLINE3, nd);
    for (int cl = 0; cl < 2; cl++) {
      g1k4(BASIS_BSPLINE2, nd, cl == 1);
      g1k4(BASIS_BSPLINE3, nd, cl == 1);
    }
  }
  g1Cardinal();
  const double sSmall = g8(0.02, 120, true);
  const double sLarge = g8(3.0, 120, true);
  g8(0.02, 120, true, 0.05);            // G8-F-cut (critic A1): the face-layer cutoff
  g8(3.0, 120, true, 0.05);
  check(sSmall < 0.01, "G8 small amplitude stays small", sSmall, 0.01);
  check(sLarge > 0.05, "G8 large amplitude reaches |F_ii - 1| > 0.05", sLarge, 0.05);
  bodyForceEqualsGravity();
  if (argc > 1 && std::strcmp(argv[1], "mms") == 0) {
    const double Fmax = (argc > 2) ? std::atof(argv[2]) : 1e-3;
    const double c = std::sqrt((1.0 * 0.7 / (1.3 * 0.4)) / 1.0);
    const double dt = 0.1 * (1.0 / 32) / c;
    for (Basis b : {BASIS_BSPLINE2, BASIS_BSPLINE3, BASIS_GIMP, BASIS_LINEAR}) {
      double prev = 0;
      for (int nc : {8, 16, 32}) {
        const double e = mmsRun(b, nc, Fmax, dt, 1.0);
        std::printf("INFO MMS %s Fmax %g 1/h %d e %.4e%s\n", basisName(b), Fmax, nc, e,
                    prev > 0 ? (" rate " + std::to_string(std::log(prev / e) / std::log(2.0))).c_str() : "");
        prev = e;
      }
    }
  }
  std::printf("[checks] %d checks, %d failed\n", g_checks, g_fail);
  std::printf("RESULT: %s\n", g_fail ? "FAIL" : "PASS");
  return g_fail ? 1 : 0;
}
