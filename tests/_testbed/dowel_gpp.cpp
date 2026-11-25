// Ladruno WP-257: g++ gate for SRC/material/nD/LadrunoDowelKernel.h (bar dowel action).
//
//   g++ -O2 -std=c++17 -I SRC/material/nD tests/_testbed/dowel_gpp.cpp -o dowel
//   ./dowel            -> self-checks D1-D5, "GPP DOWEL GATE: PASS", exit 0
//   ./dowel path       -> reads "fc fy db ls kd su cd" then n rows "gamma sigma commit", prints "tau dg dsig"
//
// D1 capacity closed forms; D2 envelope closed forms; D3 tangent vs central differences (envelope,
// secant, axial coupling); D4 origin-oriented path (unload to zero stress at zero slip, reload on the
// secant, rejoin the envelope, symmetric in sign); D5 axial interaction (yield kills the dowel,
// tension and compression alike); D6 cyclic decay (monotonic path unchanged, fixed-amplitude cycles
// lose strength every cycle, tangent, a committed state re-solves to itself).

#include "LadrunoDowelKernel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace ladruno_dowel;

static int nfail = 0;
static void check(bool ok, const char* what, double a = 0.0, double b = 0.0)
{
  if (!ok) { ++nfail; std::printf("  FAIL %s: %.17g vs %.17g\n", what, a, b); }
}
static bool close(double a, double b, double rel, double abs_ = 0.0)
{
  return std::fabs(a - b) <= rel * std::fabs(b) + abs_;
}

static Params bar()
{
  Params P;
  P.fc = 33.0; P.fy = 331.0; P.db = 25.4; P.ls = 150.0; P.kd = 0.4;
  return P;
}

static int self_checks()
{
  const Params P = bar();
  const double tu = 0.4 * (4.0 / 3.14159265358979323846) * std::sqrt(33.0 * 331.0);

  std::printf("D1 capacity\n");
  double dc;
  check(close(cap(P, 0.0, &dc), tu, 1e-15), "cap(0)", cap(P, 0.0, 0), tu);
  check(dc == 0.0, "dcap(0)", dc, 0.0);
  check(close(cap(P, 0.6 * P.fy, 0), 0.8 * tu, 1e-14), "cap(0.6 fy) = 0.8", cap(P, 0.6 * P.fy, 0), 0.8 * tu);
  check(close(cap(P, -0.6 * P.fy, 0), 0.8 * tu, 1e-14), "cap(-0.6 fy) = 0.8", cap(P, -0.6 * P.fy, 0), 0.8 * tu);
  check(cap(P, P.fy, &dc) == 0.0 && dc == 0.0, "cap(fy) = 0");
  check(cap(P, -1.2 * P.fy, &dc) == 0.0 && dc == 0.0, "cap(-1.2 fy) = 0");

  std::printf("D2 envelope\n");
  const double su = 0.05 * P.db;
  Out o = shear(P, su / P.ls, 0.0, 0.0);
  check(close(o.tau, tu * (1.0 - std::exp(-3.0)), 1e-14), "tau(su) = 0.95 cap", o.tau, tu * (1.0 - std::exp(-3.0)));
  o = shear(P, 0.0, 0.0, 0.0);
  check(o.tau == 0.0, "tau(0) = 0", o.tau, 0.0);
  check(close(o.dg, tu * P.ls / (su / 3.0), 1e-14), "virgin slope", o.dg, tu * P.ls / (su / 3.0));
  o = shear(P, -su / P.ls, 0.0, 0.0);
  check(close(o.tau, -tu * (1.0 - std::exp(-3.0)), 1e-14), "odd in gamma", o.tau, 0.0);
  Params Q = P; Q.su = 2.0;
  o = shear(Q, 2.0 / Q.ls, 0.0, 0.0);
  check(close(o.tau, tu * (1.0 - std::exp(-3.0)), 1e-14), "-dowelSlip", o.tau, 0.0);

  std::printf("D3 tangent vs central differences\n");
  const double gs[] = {1e-5, 3e-4, 2e-3, -7e-4, 1.2e-2};
  const double ss[] = {0.0, 120.0, -250.0, 320.0};
  const double sm[] = {0.0, 0.3, 1.5};
  double worst = 0.0;
  for (double g : gs) for (double s : ss) for (double m : sm) {
    const double sl = std::fabs(g) * P.ls;
    if (std::fabs(sl - m) < 1e-3) continue;              // the kink at s = smax
    const Out c = shear(P, g, s, m);
    const double hg = 1e-7 * std::max(std::fabs(g), 1e-6), hs = 1e-6 * P.fy;
    const double fdg = (shear(P, g + hg, s, m).tau - shear(P, g - hg, s, m).tau) / (2 * hg);
    const double fds = (shear(P, g, s + hs, m).tau - shear(P, g, s - hs, m).tau) / (2 * hs);
    const double eg = std::fabs(fdg - c.dg) / std::max(std::fabs(c.dg), 1.0);
    const double es = std::fabs(fds - c.dsig) / std::max(std::fabs(c.dsig), 1e-3);
    worst = std::max(worst, std::max(eg, es));
    check(eg < 1e-6, "d tau / d gamma", fdg, c.dg);
    check(es < 1e-5, "d tau / d sigma", fds, c.dsig);
  }
  std::printf("   worst relative difference %.2e\n", worst);

  std::printf("D4 origin-oriented path\n");
  double smax = 0.0;
  const double g1 = 4.0e-3;
  const Out a = shear(P, g1, 0.0, smax);
  smax = commitSlip(P, g1, smax);
  check(smax == g1 * P.ls, "commit stores the slip", smax, g1 * P.ls);
  check(shear(P, 0.0, 0.0, smax).tau == 0.0, "zero stress at zero slip");
  check(close(shear(P, 0.5 * g1, 0.0, smax).tau, 0.5 * a.tau, 1e-14), "secant at half", shear(P, 0.5 * g1, 0.0, smax).tau, 0.5 * a.tau);
  check(close(shear(P, -0.5 * g1, 0.0, smax).tau, -0.5 * a.tau, 1e-14), "secant, other sign", 0, 0);
  check(close(shear(P, g1, 0.0, smax).tau, a.tau, 1e-15), "back on the envelope at smax", 0, 0);
  const double g2 = 6.0e-3;
  check(close(shear(P, g2, 0.0, smax).tau, shear(P, g2, 0.0, 0.0).tau, 1e-15), "beyond smax = virgin envelope", 0, 0);
  check(commitSlip(P, 0.5 * g1, smax) == smax, "smaller slip keeps smax");
  check(commitSlip(P, -g2, smax) == g2 * P.ls, "either sign raises smax");

  std::printf("D5 axial interaction\n");
  check(shear(P, g1, P.fy, 0.0).tau == 0.0, "yielded in tension: no dowel");
  check(shear(P, g1, -P.fy, 0.0).tau == 0.0, "yielded in compression: no dowel");
  check(close(shear(P, g1, 0.6 * P.fy, 0.0).tau, 0.8 * shear(P, g1, 0.0, 0.0).tau, 1e-14), "0.6 fy -> 0.8", 0, 0);

  std::printf("D6 cyclic decay\n");
  {
    Params Q = P; Q.cd = 0.5;
    double sm = 0.0, sc = 0.0, D = 0.0;
    for (int i = 1; i <= 20; i++) {                // monotonic: bit-identical to no decay
      const double g = 4e-4 * i;
      const Out a1 = shear(P, g, 0.0, sm), a2 = shear(Q, g, 0.0, sm, sc, D);
      check(a1.tau == a2.tau && a1.dg == a2.dg, "monotonic: decay adds nothing", a2.tau, a1.tau);
      D = commitReslip(Q, g, sm, sc, D); sc = g * Q.ls; sm = commitSlip(Q, g, sm);
    }
    check(D == 0.0, "monotonic: D = 0", D, 0.0);
    sm = sc = D = 0.0;                             // fixed-amplitude cycles: the peak falls every cycle
    double prev = 1e300; int falls = 0;
    const double ga = 3e-3;
    for (int cyc = 0; cyc < 6; cyc++) {
      double peak = 0.0;
      for (int i = 0; i < 40; i++) {
        const double g = ga * std::sin(2.0 * 3.14159265358979323846 * (i + 1) / 40.0);
        const Out o = shear(Q, g, 0.0, sm, sc, D);
        peak = std::max(peak, std::fabs(o.tau));
        const double D1 = commitReslip(Q, g, sm, sc, D), sc1 = std::fabs(g) * Q.ls, sm1 = commitSlip(Q, g, sm);
        const Out r = shear(Q, g, 0.0, sm1, sc1, D1);   // the committed state re-solves to itself
        check(close(r.tau, o.tau, 1e-15, 1e-300), "fixed point", r.tau, o.tau);
        D = D1; sc = sc1; sm = sm1;
      }
      if (cyc > 0 && peak < prev) falls++;
      prev = peak;
    }
    check(falls == 5, "peak falls every cycle", falls, 5);
    const double gs2[] = {1e-3, -2e-3, 2.5e-3};    // tangent inside the re-slip branch
    for (double g : gs2) {
      const Out c = shear(Q, g, 100.0, 0.5, 0.05, 2.0);
      const double hg = 1e-7 * std::fabs(g);
      const double fd = (shear(Q, g + hg, 100.0, 0.5, 0.05, 2.0).tau - shear(Q, g - hg, 100.0, 0.5, 0.05, 2.0).tau) / (2 * hg);
      check(std::fabs(fd - c.dg) < 1e-6 * std::max(std::fabs(c.dg), 1.0), "decay tangent", fd, c.dg);
    }
  }

  if (nfail == 0) { std::printf("GPP DOWEL GATE: PASS\n"); return 0; }
  std::printf("GPP DOWEL GATE: %d FAILURES\n", nfail);
  return 1;
}

static int path_mode()
{
  Params P;
  if (std::scanf("%lf %lf %lf %lf %lf %lf %lf", &P.fc, &P.fy, &P.db, &P.ls, &P.kd, &P.su, &P.cd) != 7) return 2;
  int n;
  if (std::scanf("%d", &n) != 1) return 2;
  double smax = 0.0, sc = 0.0, D = 0.0;
  for (int i = 0; i < n; i++) {
    double g, s; int c;
    if (std::scanf("%lf %lf %d", &g, &s, &c) != 3) return 2;
    const Out o = shear(P, g, s, smax, sc, D);
    std::printf("%.17g %.17g %.17g\n", o.tau, o.dg, o.dsig);
    if (c) { D = commitReslip(P, g, smax, sc, D); sc = std::fabs(g) * P.ls; smax = commitSlip(P, g, smax); }
  }
  return 0;
}

int main(int argc, char** argv)
{
  if (argc > 1 && std::strcmp(argv[1], "path") == 0) return path_mode();
  return self_checks();
}
