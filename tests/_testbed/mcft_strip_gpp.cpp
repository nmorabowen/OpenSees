// Standalone g++ gate for Ladruno WP-254 (program #1131): the Model C MCFT strip kernel
// (SRC/material/nD/LadrunoMCFTKernel.h) and its fibre section through the WP-245 shear map
// (SRC/material/section/LadrunoFiberShearKernel.h). No OpenSees.
//   g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/section mcft_strip_gpp.cpp -o mcft_strip
//
// No argument: self-asserting checks (non-zero exit and no PASS line on any failure)
//   K1 lazy scan = full scan: the kernel's lazy, distance-ordered root scan picks the same
//      bracket, and so the same eps_y bit for bit, as evaluating all 81 candidates and taking
//      the first nearest increasing sign change (the reference's rule), on random committed
//      states and strains.
//   K2 steel: backbone closed forms, clip-to-backbone, fracture, and the analytic tangent vs a
//      one-sided difference.
//   K3 path independence: a trial is a function of (committed memory, strain) only; A, B, A
//      gives A's stresses bit for bit.
//   K4 transverse equilibrium: at the solved eps_y, |sy + rho f_st| is at bisection round-off
//      wherever the scan bracketed a root.
//   K5 tangent: the forward-difference condensed tangent vs a central difference with a
//      smaller step, at uncracked and cracked smooth states; uncracked = diag(Ec, Ec/2).
//   C1-C3 (WP-256, -cyclicDegradation on): a monotonic path is bit-identical to the option off;
//      a reversal is bit-identical until the first committed crack re-opening; repeated cycles at
//      one amplitude lose peak shear stress every cycle, while the option off repeats; C4 a committed
//      state re-solves to itself (the damage is implicit).
//   K6 root rule: before the first commit the full range is bisected (bracketed flag 0);
//      after it a root near the committed eps_y is taken (flag 1).
//
// "strip" on stdin: a strip path, printed for tests/test_ladrunoMCFTStrip_cpp.py to compare with
// the numpy reference (tests/_testbed/mcft_strip_ref.py):
//   STRIP fc fp ep e20 smx smy agg kbeta cvci rho fyt fut Es esh esu
//   n  then n lines  ex g commit      -> per line: sx txy ey e1 e2 w fst bracketed
// "stripc": the same with three more strip numbers, ccyc mcyc cdrv (WP-256 cyclic strut degradation)
// "section" on stdin: a Model C section (strips + bars) through the parabolic modelC shear map:
//   SECTION hy hz nstrip nbar
//   nstrip lines  y z A fc fp ep e20 smx smy agg kbeta cvci rho fyt fut Es esh esu
//   nbar lines    y z A fy fu Es esh esu
//   n  then n lines  eps0 kappa gamma commit   -> per line: N M V
#include "LadrunoMCFTKernel.h"
#include "LadrunoFiberShearKernel.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <random>
#include <string>

using namespace ladruno_mcft;
namespace sh = ladruno_shear;

static int nfail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nfail++; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static double rel(double a, double b) { return fabs(a - b) / (fabs(a) + fabs(b) + 1e-300); }

// a lightly tied core strip (rho fyt < fcr) and a cover strip, 3CLH18-like
static Strip coreStrip()
{
  Strip S;
  S.conc.fc = 25.6; S.conc.fp = 26.4; S.conc.ep = 0.00231; S.conc.e20 = 0.0251;
  S.conc.smx = 152.4; S.conc.smy = 457.2; S.conc.derive();
  S.rho = 0.001; S.tie.fy = 400.0; S.tie.fu = 600.0;
  return S;
}
static Strip coverStrip()
{
  Strip S = coreStrip();
  S.conc.fp = 25.6; S.conc.ep = 0.002; S.conc.e20 = 0.0158;
  S.rho = 0.0;
  return S;
}

// the reference's scan, written out: all 81 residuals, first nearest increasing sign change
static double solveFull(const Strip& S, const Hist& H, double ex, double g, int& found)
{
  Out o;
  double lo = kFullLo, hi = kFullHi;
  found = 0;
  if (H.used) {
    double Y[kNScan], F[kNScan];
    for (int i = 0; i < kNScan; ++i) {
      const double cand = (i < 40) ? -kScanOffsets[39 - i] : (i == 40 ? 0.0 : kScanOffsets[i - 41]);
      Y[i] = std::min(std::max(H.ey + cand, kBoxLo), kBoxHi);
      F[i] = resid(S, H, ex, Y[i], g, o);
    }
    double best = INFINITY; int kb = -1;
    for (int j = 0; j < kNScan - 1; ++j) {
      if (!(F[j] <= 0 && F[j + 1] >= 0)) continue;
      const double d = fabs(0.5 * (Y[j] + Y[j + 1]) - H.ey);
      if (d < best) { best = d; kb = j; }
    }
    if (kb >= 0) { lo = Y[kb]; hi = Y[kb + 1]; found = 1; }
  }
  for (int it = 0; it < kNBisect; ++it) {
    const double m = 0.5 * (lo + hi);
    if (resid(S, H, ex, m, g, o) <= 0) lo = m; else hi = m;
  }
  return 0.5 * (lo + hi);
}

static void runChecks()
{
  std::mt19937_64 rng(254);
  std::uniform_real_distribution<double> U(0.0, 1.0);

  // ---------------------------------------------------------------- K1
  {
    printf("K1 lazy scan = full scan\n");
    int n = 0, nfound = 0, bad = 0;
    for (int t = 0; t < 4000; ++t) {
      Strip S = (t % 3 == 0) ? coverStrip() : coreStrip();
      Hist H;
      H.used = true;
      H.ey = -0.002 + 0.02 * U(rng) * U(rng);
      H.em1 = 0.01 * U(rng) * U(rng);
      H.em2 = 0.006 * U(rng) * U(rng);
      H.ept = (U(rng) < 0.3) ? 0.004 * U(rng) : 0.0;
      const double ex = -0.004 + 0.008 * U(rng);
      const double g = (U(rng) - 0.5) * 0.02;
      double eyL; int fl;
      fl = solveEy(S, H, ex, g, eyL);
      int ff; const double eyF = solveFull(S, H, ex, g, ff);
      ++n; nfound += ff;
      if (!(eyL == eyF) || fl != ff) { if (bad < 5) printf("  mismatch t=%d ex=%g g=%g: %.17g vs %.17g (%d/%d)\n", t, ex, g, eyL, eyF, fl, ff); ++bad; }
    }
    CHECK(bad == 0, "K1: %d of %d lazy/full scans differ", bad, n);
    CHECK(nfound > n / 2 && nfound < n, "K1: the sample must exercise both bracketed (%d) and fallback cases of %d", nfound, n);
    printf("  %d states, %d bracketed, %d mismatches\n", n, nfound, bad);
  }

  // ---------------------------------------------------------------- K2
  {
    printf("K2 steel law\n");
    Steel T; T.fy = 400.0; T.fu = 600.0;
    const double ey = T.fy / T.Es;
    CHECK(steelEnv(T, 0.5 * ey) == T.Es * 0.5 * ey, "K2 elastic");
    CHECK(steelEnv(T, 0.005) == T.fy, "K2 plateau");
    const double a = 0.064, r = (a - T.esh) / (T.esu - T.esh);
    CHECK(rel(steelEnv(T, a), T.fy + (T.fu - T.fy) * (1 - (1 - r) * (1 - r))) < 1e-15, "K2 hardening");
    CHECK(rel(steelEnv(T, T.esu), T.fu) < 1e-15, "K2 fu at esu");
    CHECK(steelEnv(T, 0.13) == 1e-3 * T.fy, "K2 fracture");
    // unload from a plastic strain: elastic, then clipped by the compression backbone
    const double ep = 0.01;
    CHECK(steelStress(T, 0.011, ep) == T.Es * (0.011 - ep), "K2 elastic unload");
    CHECK(steelStress(T, 0.0, ep) == -T.fy, "K2 clipped by -fy");
    CHECK(steelStress(T, 0.03, 0.0) == steelEnv(T, 0.03), "K2 tension backbone");
    double maxe = 0.0;
    const double pts[][2] = {{0.001, 0.0}, {0.005, 0.0}, {0.05, 0.0}, {0.011, 0.01}, {-0.03, 0.0}, {0.002, 0.004}};
    for (auto& p : pts) {
      double Et; const double f0 = steelStress(T, p[0], p[1], &Et);
      const double h = 1e-9, fd = (steelStress(T, p[0] + h, p[1]) - f0) / h;
      maxe = std::max(maxe, fabs(Et - fd) / T.Es);
    }
    CHECK(maxe < 1e-5, "K2 tangent vs FD: %g", maxe);
    printf("  closed forms ok, tangent vs FD %.2e (of Es)\n", maxe);
  }

  // ---------------------------------------------------------------- K3
  {
    printf("K3 path independence\n");
    Strip S = coreStrip();
    Hist H; double ey; Out o;
    evaluate(S, H, -2e-4, 0.002, ey, o); commit(S, H, ey, o);
    evaluate(S, H, -2e-4, 0.004, ey, o); commit(S, H, ey, o);
    Out a1, b, a2; double y1, y2, yb;
    evaluate(S, H, -1e-4, 0.006, y1, a1);
    evaluate(S, H, -3e-3, -0.004, yb, b);
    evaluate(S, H, -1e-4, 0.006, y2, a2);
    CHECK(a1.sx == a2.sx && a1.txy == a2.txy && y1 == y2, "K3 A,B,A");
    printf("  A,B,A bit-identical: %s\n", (a1.sx == a2.sx && a1.txy == a2.txy) ? "yes" : "no");
  }

  // ---------------------------------------------------------------- K4
  {
    printf("K4 transverse equilibrium\n");
    Strip S = coreStrip();
    Hist H; double ey; Out o;
    double worst = 0.0;
    for (int k = 1; k <= 60; ++k) {
      const double g = 0.0002 * k, ex = -1e-4 + 2e-5 * k;
      const int f = evaluate(S, H, ex, g, ey, o);
      const double r = fabs(o.sy + S.rho * o.fst);
      if (f || !H.used) worst = std::max(worst, r / (fabs(o.sx) + fabs(o.txy) + 1.0));
      commit(S, H, ey, o);
    }
    CHECK(worst < 1e-9, "K4 relative residual %g", worst);
    printf("  worst |sy + rho fst| / |s| = %.2e\n", worst);
  }

  // ---------------------------------------------------------------- K5
  {
    printf("K5 tangent\n");
    Strip S = coreStrip();
    Hist H; double ey; Out o;
    evaluate(S, H, 1e-5, 0.0, ey, o);
    double D[2][2];
    tangentFD(S, H, 1e-5, 0.0, o, 1e-7, 1e-7, D);
    const double e00 = rel(D[0][0], S.conc.Ec), e11 = rel(D[1][1], 0.5 * S.conc.Ec);
    CHECK(e00 < 1e-6 && e11 < 1e-6 && fabs(D[0][1]) < 1e-6 * S.conc.Ec && fabs(D[1][0]) < 1e-6 * S.conc.Ec,
          "K5 uncracked tangent %g %g %g %g", D[0][0], D[0][1], D[1][0], D[1][1]);
    // cracked smooth state: drive to cracking with commits, then compare with a central difference
    for (int k = 1; k <= 20; ++k) { evaluate(S, H, -2e-4, 0.0004 * k, ey, o); commit(S, H, ey, o); }
    const double ex0 = -2e-4, g0 = 0.0081;
    evaluate(S, H, ex0, g0, ey, o);
    tangentFD(S, H, ex0, g0, o, 1e-7, 1e-7, D);
    const double h = 1e-8;
    double ep_, em_; Out op, om;
    double C[2][2];
    evaluate(S, H, ex0 + h, g0, ep_, op); evaluate(S, H, ex0 - h, g0, em_, om);
    C[0][0] = (op.sx - om.sx) / (2 * h); C[1][0] = (op.txy - om.txy) / (2 * h);
    evaluate(S, H, ex0, g0 + h, ep_, op); evaluate(S, H, ex0, g0 - h, em_, om);
    C[0][1] = (op.sx - om.sx) / (2 * h); C[1][1] = (op.txy - om.txy) / (2 * h);
    double sc = 0.0, er = 0.0;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) sc = std::max(sc, fabs(C[i][j]));
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) er = std::max(er, fabs(D[i][j] - C[i][j]) / sc);
    CHECK(er < 1e-3, "K5 cracked tangent vs central FD: %g", er);
    printf("  uncracked = diag(Ec, Ec/2) to %.1e; cracked forward vs central %.2e\n", std::max(e00, e11), er);
  }

  // ---------------------------------------------------------------- K6
  {
    printf("K6 root rule\n");
    Strip S = coreStrip();
    Hist H; double ey; Out o;
    const int f0 = evaluate(S, H, -1e-4, 0.001, ey, o);
    CHECK(f0 == 0, "K6 before the first commit the full range is bisected");
    commit(S, H, ey, o);
    const int f1 = evaluate(S, H, -1e-4, 0.0012, ey, o);
    CHECK(f1 == 1, "K6 after it the root near the committed eps_y is bracketed");
    printf("  flags %d then %d\n", f0, f1);
  }
}

// ------------------------------------------------------------------ WP-256 cyclic degradation
struct Step { double sx, txy, ey, reopen; };

static std::vector<Step> drive(Strip S, const std::vector<double>& ex, const std::vector<double>& g)
{
  Hist H; std::vector<Step> out;
  for (size_t k = 0; k < g.size(); ++k) {
    double ey; Out o;
    evaluate(S, H, ex[k], g[k], ey, o);
    commit(S, H, ey, o);
    out.push_back({o.sx, o.txy, ey, H.reopen});   // the re-opening this step's trial used (implicit)
  }
  return out;
}

static void runCyclicChecks()
{
  const double CC = 300.0;
  // ---------------------------------------------------------------- C1
  {
    printf("C1 monotonic path unchanged\n");
    std::vector<double> ex, g;
    for (int k = 0; k <= 150; ++k) { ex.push_back(-2e-4); g.push_back(1e-4 * k); }
    int bad = 0;
    for (int which = 0; which < 2; ++which) {
      Strip off = which ? coverStrip() : coreStrip(), on = off;
      on.ccyc = CC;
      const auto a = drive(off, ex, g), b = drive(on, ex, g);
      for (size_t k = 0; k < a.size(); ++k)
        if (!(a[k].sx == b[k].sx && a[k].txy == b[k].txy && a[k].ey == b[k].ey) || b[k].reopen != 0.0) ++bad;
    }
    CHECK(bad == 0, "C1: %d monotonic steps differ with the option on", bad);
    printf("  core and cover, 151 steps each: %d differ, re-opening stays 0\n", bad);
  }
  // ---------------------------------------------------------------- C2
  {
    printf("C2 unchanged until the first re-opening\n");
    // up to gamma 0.006, back to -0.006: identical until the crack re-opens on the reversal
    std::vector<double> ex, g;
    for (int k = 0; k <= 60; ++k) { ex.push_back(-2e-4); g.push_back(1e-4 * k); }
    for (int k = 1; k <= 120; ++k) { ex.push_back(-2e-4); g.push_back(0.006 - 1e-4 * k); }
    Strip off = coreStrip(), on = off; on.ccyc = CC;
    const auto a = drive(off, ex, g), b = drive(on, ex, g);
    int firstDiff = -1, firstReopen = -1;
    for (size_t k = 0; k < a.size(); ++k) {
      if (firstDiff < 0 && !(a[k].sx == b[k].sx && a[k].txy == b[k].txy)) firstDiff = (int)k;
      if (firstReopen < 0 && b[k].reopen > 0.0) firstReopen = (int)k;
    }
    CHECK(firstReopen > 60, "C2: re-opening must start after the reversal (step %d)", firstReopen);
    CHECK(firstDiff >= firstReopen, "C2: a difference at step %d before the first re-opening (%d)", firstDiff, firstReopen);
    CHECK(firstDiff > 0, "C2: the option must act after the re-opening");
    printf("  reversal at step 60; first re-opening committed at %d, first difference at %d\n", firstReopen, firstDiff);
  }
  // ---------------------------------------------------------------- C3
  {
    printf("C3 repeated cycles lose strut strength\n");
    // 8 full cycles at gamma = +-0.006 under axial compression
    std::vector<double> ex, g;
    const int q = 30;
    ex.push_back(-3e-4); g.push_back(0.0);
    for (int c = 0; c < 8; ++c) {
      for (int k = 1; k <= q; ++k) { ex.push_back(-3e-4); g.push_back(0.006 * k / q); }
      for (int k = 1; k <= 2 * q; ++k) { ex.push_back(-3e-4); g.push_back(0.006 - 0.006 * k / q); }
      for (int k = 1; k <= q; ++k) { ex.push_back(-3e-4); g.push_back(-0.006 + 0.006 * k / q); }
    }
    for (int which = 0; which < 2; ++which) {
      Strip off = which ? coverStrip() : coreStrip(), on = off; on.ccyc = CC;
      const auto a = drive(off, ex, g), b = drive(on, ex, g);
      std::vector<double> pa, pb, ra;
      for (int c = 0; c < 8; ++c) {
        const size_t kp = 1 + (size_t)c * 4 * q + q - 1;   // the +peak of cycle c
        pa.push_back(fabs(a[kp + 1].txy)); pb.push_back(fabs(b[kp + 1].txy)); ra.push_back(b[kp + 1].reopen);
      }
      int nonDecr = 0, reopenFlat = 0;
      for (int c = 1; c < 8; ++c) { if (!(pb[c] < pb[c - 1])) ++nonDecr; if (!(ra[c] > ra[c - 1])) ++reopenFlat; }
      const double offDrift = fabs(pa[7] - pa[1]) / pa[1];
      // the strut strength factor 1/(1 + ccyc P) falls every cycle in both strips; the shear stress of
      // the tied strip (the one that carries the shear) falls with it. An untied cover strip
      // (sigma_y = 0) rotates instead: its small residual shear may rise.
      CHECK(reopenFlat == 0, "C3 (%s): re-opening (strut damage) must grow every cycle", which ? "cover" : "core");
      if (!which) {
        CHECK(nonDecr == 0, "C3 (core): peak shear stress must fall every cycle with the option on (%d rises)", nonDecr);
        CHECK(pb[7] < 0.9 * pa[7], "C3 (core): 8 cycles must cost more than 10 %% (%g vs %g)", pb[7], pa[7]);
      }
      CHECK(offDrift < 1e-12, "C3 (%s): the option off must repeat cycle 2 on cycle 8 (%g)", which ? "cover" : "core", offDrift);
      printf("  %s: peak |txy| cycle 1..8 off %.3f -> %.3f (drift %.1e), on %.3f -> %.3f MPa, P %.4f\n",
             which ? "cover" : "core", pa[0], pa[7], offDrift, pb[0], pb[7], ra[7]);
    }
  // ---------------------------------------------------------------- C4
  {
    printf("C4 a committed state re-solves to itself\n");
    // implicit damage: re-evaluating the committed strain under the new memory gives the same stress
    Strip S = coreStrip(); S.ccyc = CC;
    Hist H; int bad = 0, n = 0; double worst = 0.0;
    for (int k = 0; k <= 240; ++k) {
      const double g = 0.006 * sin(k * PI / 30.0), ex = -3e-4;
      double ey, ey2; Out o, o2;
      evaluate(S, H, ex, g, ey, o);
      commit(S, H, ey, o);
      evaluate(S, H, ex, g, ey2, o2);
      ++n;
      const double dd = std::max(fabs(o2.sx - o.sx), fabs(o2.txy - o.txy)) / S.conc.fc;
      worst = std::max(worst, dd);
      if (dd > 1e-12) ++bad;
    }
    CHECK(bad == 0, "C4: %d of %d committed states move when re-solved", bad, n);
    printf("  %d cyclic steps, %d move by more than 1e-12 fc (worst %.1e), re-opening %.4f\n", n, bad, worst, H.reopen);
  }
  }
}

static bool readStrip(FILE* in, Strip& S)
{
  double v[15];
  for (int i = 0; i < 15; ++i) if (fscanf(in, "%lf", &v[i]) != 1) return false;
  S.conc.fc = v[0]; S.conc.fp = v[1]; S.conc.ep = v[2]; S.conc.e20 = v[3]; S.conc.smx = v[4];
  S.conc.smy = v[5]; S.conc.agg = v[6]; S.conc.kbeta = v[7]; S.conc.cvci = v[8]; S.conc.derive();
  S.rho = v[9]; S.tie.fy = v[10]; S.tie.fu = v[11]; S.tie.Es = v[12]; S.tie.esh = v[13]; S.tie.esu = v[14];
  return true;
}

static int runStripPath(FILE* in, bool cyclic = false)
{
  char tok[32];
  if (fscanf(in, "%31s", tok) != 1 || strcmp(tok, "STRIP") != 0) return 2;
  Strip S;
  if (!readStrip(in, S)) return 2;
  if (cyclic && fscanf(in, "%lf %lf %d", &S.ccyc, &S.mcyc, &S.cdrv) != 3) return 2;   // WP-256: ccyc mcyc cdrv
  int n; if (fscanf(in, "%d", &n) != 1) return 2;
  Hist H;
  for (int k = 0; k < n; ++k) {
    double ex, g; int c;
    if (fscanf(in, "%lf %lf %d", &ex, &g, &c) != 3) return 2;
    double ey; Out o;
    const int f = evaluate(S, H, ex, g, ey, o);
    printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %d\n", o.sx, o.txy, ey, o.e1, o.e2, o.w, o.fst, f);
    if (c) commit(S, H, ey, o);
  }
  return 0;
}

struct Fib { double y, z, A; int kind; Strip S; Steel B; };

static int runSectionPath(FILE* in)
{
  char tok[32];
  if (fscanf(in, "%31s", tok) != 1 || strcmp(tok, "SECTION") != 0) return 2;
  double hy, hz; int ns, nb;
  if (fscanf(in, "%lf %lf %d %d", &hy, &hz, &ns, &nb) != 4) return 2;
  std::vector<Fib> F(ns + nb);
  for (int i = 0; i < ns; ++i) {
    if (fscanf(in, "%lf %lf %lf", &F[i].y, &F[i].z, &F[i].A) != 3 || !readStrip(in, F[i].S)) return 2;
    F[i].kind = 0;
  }
  for (int i = ns; i < ns + nb; ++i) {
    Steel& B = F[i].B;
    if (fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf", &F[i].y, &F[i].z, &F[i].A, &B.fy, &B.fu, &B.Es, &B.esh, &B.esu) != 8)
      return 2;
    F[i].kind = 1;
  }
  sh::Spec spec; spec.profile = sh::PARABOLIC; spec.resultant = sh::MODELC; spec.hy = hy; spec.hz = hz;
  std::vector<Hist> H(ns + nb);
  std::vector<double> epb(ns + nb, 0.0), eyT(ns + nb), ebT(ns + nb);
  std::vector<Out> oT(ns + nb);
  int n; if (fscanf(in, "%d", &n) != 1) return 2;
  for (int k = 0; k < n; ++k) {
    double e0, kap, gam; int c;
    if (fscanf(in, "%lf %lf %lf %d", &e0, &kap, &gam, &c) != 4) return 2;
    const double ev[6] = {e0, kap, 0.0, gam, 0.0, 0.0};
    double s6[6] = {0, 0, 0, 0, 0, 0};
    for (size_t f = 0; f < F.size(); ++f) {
      double a[3][6], b[3][6], eps[3], sig[3] = {0, 0, 0};
      sh::maps(spec, F[f].y, F[f].z, a, b);
      sh::strain(a, ev, eps);
      if (F[f].kind == 0) {
        evaluate(F[f].S, H[f], eps[0], eps[1], eyT[f], oT[f]);
        sig[0] = oT[f].sx; sig[1] = oT[f].txy;
      } else {
        ebT[f] = eps[0];
        sig[0] = steelStress(F[f].B, eps[0], epb[f]);
      }
      sh::accumulate(F[f].A, a, b, sig, 0, s6, 0);
    }
    printf("%.17g %.17g %.17g\n", s6[0], s6[1], s6[3]);
    if (c)
      for (size_t f = 0; f < F.size(); ++f) {
        if (F[f].kind == 0) commit(F[f].S, H[f], eyT[f], oT[f]);
        else epb[f] = ebT[f] - steelStress(F[f].B, ebT[f], epb[f]) / F[f].B.Es;
      }
  }
  return 0;
}

int main(int argc, char** argv)
{
  if (argc > 1 && strcmp(argv[1], "strip") == 0) return runStripPath(stdin);
  if (argc > 1 && strcmp(argv[1], "stripc") == 0) return runStripPath(stdin, true);
  if (argc > 1 && strcmp(argv[1], "section") == 0) return runSectionPath(stdin);
  runChecks();
  runCyclicChecks();
  if (nfail == 0) printf("GPP MCFT STRIP GATE: PASS\n");
  else printf("GPP MCFT STRIP GATE: FAIL (%d)\n", nfail);
  return nfail == 0 ? 0 : 1;
}
