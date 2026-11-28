// WP-264: the extracted LadrunoRebarBuckling GA factor (SRC/material/uniaxial/LadrunoRebarBucklingKernel.h) checked
// OpenSees-free:
//
//  BITS   gaFactorT<double> vs a VERBATIM transcription of the pre-WP-264 LadrunoRebarBuckling::gaFactor, bit for bit;
//         and the value of the Dual<1> pass (gaFactorDeriv) == the double factor, bit for bit.
//  DERIV  dfactor/deps by Dual<1> vs a step-plateau central FD (h = 10^-k scaled to the distance to e_cross, the k
//         where successive steps agree best), off the kinks (fs_buck = 1, |fs_buck - 1| = Dft); and vs the shipped FD
//         (h = 1e-8, the epsp clamp at e_cross) it replaces.
//  TAN    the whole GA tangent rho k_bare + dsig/dfactor dfactor/deps with the dual factor vs with the shipped FD: the
//         relative difference distribution (k_bare = E, sBare = the anchor fraction, as the wrapper's GA leg sees it).
//  MUT    break-the-gate: the beta derivative dropped (beta frozen), and the sqrt chain dropped; both must miss.
//
// Prints one line per gate (GATE <name> <ok> <value>) for tests/test_ladrunoRebarBuckling_dual_cpp.py.
// Build: g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/uniaxial tests/ladruno_rebarbuckling_dual_check.cpp

#include <LadrunoRebarBucklingKernel.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using ladruno_dual::Dual;
namespace rb = ladruno_rebar_buckling;

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// --------------------------------------------------------------------------
// Verbatim transcription of LadrunoRebarBuckling::gaFactor and the FD of its tangent as they stood on origin/ladruno
// @ f183b9cd1 (members -> arguments). DO NOT "tidy" it: it is the bit-identity reference.
// --------------------------------------------------------------------------
static double legacyGaFactor(double eps, double e_cross, double lsr, double gaReduction)
{
  const double C = 9.42477796076938;               // 3*pi (upstream constant)
  double denom = e_cross - eps;                     // > 0 in the buckled region
  if (denom < 1.0e-300) denom = 1.0e-300;
  const double fs_buck = sqrt(32.0 / denom) / (C * lsr);
  double betaLoc = 1.0;
  const double gamaLoc = 0.1, Dft = 0.25;
  const double stress_diff = fabs(fs_buck - 1.0);
  if (stress_diff <= Dft) betaLoc = 1.0 - gamaLoc * (Dft - stress_diff) / Dft;
  const double m = (fs_buck < 1.0) ? fs_buck : 1.0;
  return m * betaLoc * (1.0 - gaReduction) + gaReduction;
}
static double legacyDfac(double eps, double e_cross, double lsr, double r)
{
  const double h = 1.0e-8;
  double epsm = eps - h, epsp = eps + h;
  if (epsp >= e_cross) epsp = eps;                 // stay inside the buckled region
  return (legacyGaFactor(epsp, e_cross, lsr, r) - legacyGaFactor(epsm, e_cross, lsr, r)) / (epsp - epsm);
}

// mutants of the derivative (break-the-gate)
template <int MUT>
static double mutantDfac(double eps, double e_cross, double lsr, double gaReduction)
{
  typedef Dual<1> T;
  const T x = T::variable(eps, 0);
  const double C = 9.42477796076938;
  T denom = e_cross - x;
  if (denom < 1.0e-300) denom = T(1.0e-300);
  T fs_buck = sqrt(32.0 / denom) / (C * lsr);
  if (MUT == 2) fs_buck.d[0] = (32.0 / denom.v) * x.d[0];              // M2: the sqrt chain dropped (d sqrt(u) := du)
  T betaLoc = T(1.0);
  const double gamaLoc = 0.1, Dft = 0.25;
  const T stress_diff = fabs(fs_buck - 1.0);
  if (stress_diff <= Dft) betaLoc = 1.0 - gamaLoc * (Dft - stress_diff) / Dft;
  if (MUT == 1) betaLoc.d[0] = 0.0;                                    // M1: beta frozen
  const T m = (fs_buck < 1.0) ? fs_buck : T(1.0);
  return (m * betaLoc * (1.0 - gaReduction) + gaReduction).d[0];
}

struct Dist { std::vector<double> v; void add(double x) { v.push_back(x); }
  double q(double f) const { if (v.empty()) return 0; std::vector<double> s = v; std::sort(s.begin(), s.end());
    return s[(size_t)std::min<double>(s.size() - 1, std::floor(f * (s.size() - 1) + 0.5))]; }
  double mx() const { double m = 0; for (double x : v) m = std::fmax(m, x); return m; } };

int main()
{
  const double lsrs[] = {5.0, 8.0, 10.0, 12.5, 15.0, 20.0};
  const double reds[] = {0.0, 0.2, 0.5, 0.8};
  const double ecs[] = {0.0, 2.0e-3, 1.5e-2};
  long nBits = 0, badBits = 0, badDualVal = 0;
  Dist plateau, shipped, kinkShipped, tanDiff, m1, m2;
  int noPlateau = 0, nKink = 0;
  for (double lsr : lsrs) for (double r : reds) for (double ec : ecs) {
    for (int i = 1; i <= 400; ++i) {
      // eps from just below e_cross down to 6 % below it (log-spaced distance), the buckled region
      const double dist = 1e-7 * std::pow(10.0, 5.8 * i / 400.0);
      const double eps = ec - dist;
      const double fNew = rb::gaFactorT(eps, ec, lsr, r), fOld = legacyGaFactor(eps, ec, lsr, r);
      ++nBits; if (!sameBits(fNew, fOld)) ++badBits;
      double fDual; const double dDual = rb::gaFactorDeriv(eps, ec, lsr, r, &fDual);
      if (!sameBits(fDual, fOld)) ++badDualVal;
      const double dOld = legacyDfac(eps, ec, lsr, r);
      // kinks of the factor: fs_buck == 1 (|.| and min) and |fs_buck - 1| == Dft (beta)
      const double fs = std::sqrt(32.0 / dist) / (9.42477796076938 * lsr);
      const double kinkDist = std::fmin(std::fabs(fs - 1.0), std::fabs(std::fabs(fs - 1.0) - 0.25));
      const double scale = std::fmax(std::fabs(dDual), 1e-6 / dist);   // dfactor/deps ~ factor / dist
      // step-plateau central FD, h = dist 10^-k, k = 2..9
      double g[8];
      for (int k = 2; k <= 9; ++k) {
        const double h = dist * std::pow(10.0, -k);
        g[k - 2] = (legacyGaFactor(eps + h, ec, lsr, r) - legacyGaFactor(eps - h, ec, lsr, r)) / (2.0 * h);
      }
      int best = 0; double qual = 1e300;
      for (int k = 0; k < 7; ++k) { const double q = std::fabs(g[k] - g[k + 1]) / scale; if (q < qual) { qual = q; best = k; } }
      const bool kink = kinkDist < 1e-3;
      if (kink) { ++nKink; kinkShipped.add(std::fabs(dDual - dOld) / scale); continue; }
      if (qual > 1e-6) { ++noPlateau; continue; }
      plateau.add(std::fabs(dDual - g[best + 1]) / scale);
      shipped.add(std::fabs(dDual - dOld) / scale);
      m1.add(std::fabs(mutantDfac<1>(eps, ec, lsr, r) - g[best + 1]) / scale);
      m2.add(std::fabs(mutantDfac<2>(eps, ec, lsr, r) - g[best + 1]) / scale);
      // the whole GA tangent (applyBucklingGA): rho k_bare + dsig/dfactor dfactor/deps
      const double E = 200000.0, fsup = 420.0, ff = 0.1, sBare = -300.0, kBare = E;
      const double rho = (fOld + ff) / (1.0 + ff), dSdF = -(fsup * ff - sBare) / (1.0 + ff);
      const double tNew = rho * kBare + dSdF * dDual, tOld = rho * kBare + dSdF * dOld;
      tanDiff.add(std::fabs(tNew - tOld) / std::fmax(std::fabs(tOld), 1e-12 * E));
    }
  }
  std::printf("BITS factor compares %ld, differ %ld; dual value differs %ld\n", nBits, badBits, badDualVal);
  std::printf("DERIV vs plateau FD: n %zu (no plateau %d, kink-adjacent %d) p50 %.2e p90 %.2e max %.2e\n",
              plateau.v.size(), noPlateau, nKink, plateau.q(0.5), plateau.q(0.9), plateau.mx());
  std::printf("DERIV vs shipped FD (h = 1e-8): p50 %.2e p90 %.2e max %.2e | kink-adjacent p50 %.2e max %.2e\n",
              shipped.q(0.5), shipped.q(0.9), shipped.mx(), kinkShipped.q(0.5), kinkShipped.mx());
  std::printf("TAN dual vs shipped-FD tangent (relative): p50 %.2e p90 %.2e max %.2e\n", tanDiff.q(0.5), tanDiff.q(0.9), tanDiff.mx());
  std::printf("MUT beta frozen: flagged %ld of %zu (max %.2e) | sqrt chain dropped: flagged %ld of %zu (max %.2e)\n",
              (long)std::count_if(m1.v.begin(), m1.v.end(), [](double x) { return x > 1e-6; }), m1.v.size(), m1.mx(),
              (long)std::count_if(m2.v.begin(), m2.v.end(), [](double x) { return x > 1e-6; }), m2.v.size(), m2.mx());
  const bool okBits = nBits > 20000 && badBits == 0 && badDualVal == 0;
  const bool okDeriv = plateau.v.size() > 15000 && plateau.mx() < 1e-6 && noPlateau * 50 <= (int)plateau.v.size();
  const bool okM1 = std::count_if(m1.v.begin(), m1.v.end(), [](double x) { return x > 1e-6; }) > 0;
  const bool okM2 = std::count_if(m2.v.begin(), m2.v.end(), [](double x) { return x > 1e-6; }) > 0;
  std::printf("GATE bits %d %ld\n", (int)okBits, badBits + badDualVal);
  std::printf("GATE deriv %d %.3e\n", (int)okDeriv, plateau.mx());
  std::printf("GATE mut_beta %d %.3e\n", (int)okM1, m1.mx());
  std::printf("GATE mut_sqrt %d %.3e\n", (int)okM2, m2.mx());
  const bool all = okBits && okDeriv && okM1 && okM2;
  std::printf(all ? "REBAR DUAL CHECK: ALL PASS\n" : "REBAR DUAL CHECK: FAIL\n");
  return all ? 0 : 1;
}
