// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno WP-254 (program #1131): the Model C MCFT strip, header-only.
//
// One concrete strip of a shear-critical fibre beam, after the rotating-angle
// Modified Compression Field Theory (Vecchio & Collins 1986) as the "Model C"
// fibre beam uses it (calibrated on Sezen & Moehle 2002 and Lynn et al. 1996).
// The strip is a 2D membrane in the plane of bending: it receives the fibre
// strains (eps_x, gamma_xy) and has ONE internal unknown, the transverse strain
// eps_y, solved from transverse equilibrium against smeared ties:
//
//     r(eps_y) = sigma_y,c(eps_x, eps_y, gamma_xy) + rho_t f_st(eps_y) = 0
//
// Law (tension positive, compression magnitudes c >= 0):
//   - principal strains, coaxial and rotating: e1, e2 and the angle of e1;
//   - compression: Popovics to (ep, fp), then linear to 0.2 fp at e20 and flat
//     (e20 is a strip constant; the caller regularises it with the crushing
//     energy, Coleman & Spacone 2001); f2 = -beta f(c), beta = min(1, 1/(0.8 +
//     C e1)) on the CURRENT e1 (stress-only softening);
//   - tension: Ec x to ecr, then fcr / (1 + sqrt(500 (x - ecr))) (Collins &
//     Mitchell, shifted by ecr);
//   - crack check ("sum"): f1 <= v_ci,max tan(theta) + rho_t max(fyt - f_st, 0),
//     v_ci,max = c_vci sqrt(f'c) / (0.31 + 24 w / (a + 16)), w = e1 s_theta,
//     s_theta = 1 / (|sin th|/smx + |cos th|/smy), th the strut angle from x;
//   - history (committed only): em1 = max principal tension, em2 = max principal
//     compression magnitude; compression un/reloads linearly to the Palermo &
//     Vecchio (2003) plastic offset, tension un/reloads on the secant to the
//     origin; biaxial tension gives f2 = 0, biaxial compression f1 = -f(-e1).
//   - tie steel: elastic / plateau / quadratic hardening / fracture at total
//     strain esu, elastic unloading from a committed plastic strain clipped to
//     the backbone of each sign, no Bauschinger effect.
//
// Cyclic strut degradation (WP-256, opt-in, Strip::ccyc > 0; off = Model C bit for bit):
// Model C's strut remembers only the maximum principal strains, so repeated cycles to
// one drift retrace the same envelope. With the option on, the committed memory adds the
// cumulative crack RE-opening P (the increases of (e1 - ecr)+ that stay below its
// previous maximum, each weighted by the opening reached to the power mcyc, default 0;
// driver 1 also counts the strut's re-compression below its previous maximum the same way;
// a monotonic path adds exactly zero), and the strut stress carries
//     beta_cyc = 1 / (1 + ccyc P),   P = P_committed + the trial's own re-opening increment
// (implicit: the damage of the converged trial is the damage commit stores, so a committed
// state re-solves to itself),
// on top of beta(e1): a damage indicator on the compression response driven by crack
// cycling, in the spirit of the cyclic compression damage of Palermo & Vecchio (2003).
//
// Root rule (the branch rule that matters at the cracking fold, where light
// ties, rho_t fyt < fcr, make the uncracked root vanish):
//   - before the first commit: bisection over [-0.004, 0.15];
//   - after it: scan the committed eps_y plus {0, +-geomspace(1e-7, 0.15, 40)}
//     clipped to [-0.01, 0.15] (81 points), take the INCREASING sign change
//     (r <= 0 then r >= 0) whose midpoint is nearest the committed eps_y, and
//     bisect it; with no such change, bisect [-0.004, 0.15];
//   - 60 halvings, mid-point returned.
// The scan here evaluates r lazily, interval by interval in order of distance;
// it picks the same interval as evaluating all 81 points and taking the first
// minimum, so the result is bit-for-bit the full scan's.
//
// Everything below is written so its floating-point operations follow the
// reference oracle's order (tests/_testbed/mcft_strip_ref.py): the g++ gate
// compares the two to round-off.
//
// No OpenSees dependency: tests/_testbed/mcft_strip_gpp.cpp compiles it alone.

#ifndef LadrunoMCFTKernel_h
#define LadrunoMCFTKernel_h

#include <math.h>
#include <algorithm>

namespace ladruno_mcft {

static const double PI = 3.14159265358979323846;

// ------------------------------------------------------------------- steel
struct Steel {
  double fy, fu, Es, esh, esu;
  Steel() : fy(0.0), fu(0.0), Es(200000.0), esh(0.008), esu(0.12) {}
};

// backbone at strain magnitude a (>= 0)
inline double steelEnv(const Steel& S, double a)
{
  if (a > S.esu) return 1e-3 * S.fy;
  const double ey = S.fy / S.Es;
  if (a <= ey) return S.Es * a;
  if (a <= S.esh) return S.fy;
  double r = (a - S.esh) / (S.esu - S.esh);
  r = std::min(std::max(r, 0.0), 1.0);
  const double q = 1.0 - r;
  return S.fy + (S.fu - S.fy) * (1.0 - q * q);
}

// slope of the backbone at a (the bound's derivative)
inline double steelEnvSlope(const Steel& S, double a)
{
  if (a > S.esu) return 0.0;
  const double ey = S.fy / S.Es;
  if (a <= ey) return S.Es;
  if (a <= S.esh) return 0.0;
  const double r = (a - S.esh) / (S.esu - S.esh);
  return (S.fu - S.fy) * 2.0 * (1.0 - r) / (S.esu - S.esh);
}

// stress at total strain e from committed plastic strain ep; tangent in *Et
inline double steelStress(const Steel& S, double e, double ep, double* Et = 0)
{
  const double ey = S.fy / S.Es;
  const double trial = S.Es * (e - ep);
  const double up = (e > S.esu) ? 1e-3 * S.fy : steelEnv(S, std::max(e, ey));
  const double lo = -((-e > S.esu) ? 1e-3 * S.fy : steelEnv(S, std::max(-e, ey)));
  double f = std::max(trial, lo);
  f = std::min(f, up);
  if (Et != 0) {
    // the forward (increasing-strain) slope: on the upper bound the stress stays
    // clipped; on the lower bound it leaves it elastically
    if (trial >= up)      *Et = (e > S.esu || e <= ey) ? 0.0 : steelEnvSlope(S, e);
    else if (trial < lo)  *Et = (-e > S.esu || -e <= ey) ? 0.0 : steelEnvSlope(S, -e);
    else                  *Et = S.Es;
  }
  return f;
}

// ---------------------------------------------------------------- concrete
struct Concrete {
  double fc;            // f'c (> 0)
  double fp, ep, e20;   // envelope peak stress, peak strain, strain at 0.2 fp (magnitudes)
  double smx, smy;      // crack-spacing parameters [length]
  double agg;           // maximum aggregate size [length]
  double kbeta;         // C in beta = 1/(0.8 + C e1)
  double cvci;          // v_ci,max coefficient
  double Ec, fcr, ecr;  // derived
  Concrete() : fc(0), fp(0), ep(0), e20(0), smx(0), smy(0), agg(19.0), kbeta(170.0), cvci(0.18),
               Ec(0), fcr(0), ecr(0) {}
  void derive()
  {
    Ec = 4700.0 * sqrt(fc);
    fcr = 0.33 * sqrt(fc);
    ecr = fcr / Ec;
  }
};

// compression envelope at magnitude c >= 0
inline double fEnv(const Concrete& C, double c)
{
  if (c <= C.ep) {
    const double x = c / C.ep;
    const double r = C.Ec / (C.Ec - C.fp / C.ep);
    return C.fp * x * r / (r - 1.0 + pow(x, r));
  }
  return std::max(C.fp - 0.8 * C.fp * (c - C.ep) / (C.e20 - C.ep), 0.2 * C.fp);
}

// tension envelope at x
inline double fTen(const Concrete& C, double x)
{
  if (x <= C.ecr) return C.Ec * x;
  return C.fcr / (1.0 + sqrt(500.0 * std::max(x - C.ecr, 0.0)));
}

// ------------------------------------------------------------------- strip
struct Strip {
  Concrete conc;
  Steel tie;
  double rho;           // smeared tie ratio across the strip (0: no ties)
  double ccyc;          // WP-256 cyclic strut degradation, beta_cyc = 1/(1 + ccyc P); 0 = off (Model C)
  double mcyc;          // WP-256 weight exponent: dP = (re-opening) * (opening)^mcyc; 0 = unweighted
  int    cdrv;          // WP-256 driver: 0 = crack re-opening only, 1 = also strut re-compression
  Strip() : rho(0.0), ccyc(0.0), mcyc(0.0), cdrv(0) {}
};

// committed memory of one strip (06 §2.8)
struct Hist {
  double ey;            // committed eps_y (the root reference)
  double em1, em2;      // max principal tension, max principal compression magnitude
  double ept;           // tie plastic strain
  bool   used;          // false until the first commit
  // cyclic degradation memory (WP-256; carried but unused when Strip::ccyc == 0)
  double e1c;           // committed crack opening (e1 - ecr)+ at the last commit
  double reopen;        // cumulative re-opening (and, driver 1, re-compression) below the previous maxima
  double e2c;           // committed strut compression (-e2 - 0.5 ep)+ at the last commit (driver 1)
  Hist() : ey(0.0), em1(0.0), em2(0.0), ept(0.0), used(false), e1c(0.0), reopen(0.0), e2c(0.0) {}
};

// WP-256 cumulative crack re-opening at a trial principal tension e1: the committed P plus
// the part of this step's opening increment that stays below the previous maximum opening.
// A monotonic path adds exactly zero (min(e, emax) - e1c = emax - emax).
inline double reopenTrial(const Strip& S, const Hist& H, double e1, double e2)
{
  double P = H.reopen;
  const double emaxOld = std::max(H.em1 - S.conc.ecr, 0.0);
  const double e1p = std::max(e1 - S.conc.ecr, 0.0);
  const double up = std::min(e1p, emaxOld) - H.e1c;
  // weighted by the opening reached (mcyc > 0): wide cracks grind the strut more than narrow ones
  if (up > 0.0) P = (S.mcyc != 0.0) ? P + up * pow(std::min(e1p, emaxOld), S.mcyc) : P + up;
  if (S.cdrv == 1) {
    // driver 1: also the strut's re-compression below its previous maximum compression
    // measured above 0.5 ep (micro-cracking onset), as the opening is measured above ecr
    const double c0 = 0.5 * S.conc.ep;
    const double cmaxOld = std::max(H.em2 - c0, 0.0);
    const double ccp = std::max(-e2 - c0, 0.0);
    const double uc = std::min(ccp, cmaxOld) - H.e2c;
    if (uc > 0.0) P = (S.mcyc != 0.0) ? P + uc * pow(std::min(ccp, cmaxOld), S.mcyc) : P + uc;
  }
  return P;
}

// WP-256 cyclic strut degradation factor for a given P (1 when the option is off)
inline double betaCyclic(const Strip& S, double P)
{
  return (S.ccyc > 0.0) ? 1.0 / (1.0 + S.ccyc * P) : 1.0;
}

// membrane response at a given eps_y
struct Out {
  double sx, sy, txy;   // stresses
  double e1, e2;        // principal strains
  double w;             // crack width
  double fst;           // tie stress
  double f1, f2;        // principal stresses
  double beta;          // compression softening factor
  double theta;         // strut angle from x [rad]
};

inline void membrane(const Strip& S, const Hist& H, double ex, double ey, double g, Out& o)
{
  const Concrete& C = S.conc;
  o.fst = (S.rho != 0.0) ? steelStress(S.tie, ey, H.ept) : 0.0;
  const double c0 = 0.5 * (ex + ey);
  const double hx = 0.5 * (ex - ey), hg = 0.5 * g;
  const double R = sqrt(hx * hx + hg * hg) + 1e-30;
  const double e1 = c0 + R, e2 = c0 - R;
  const double cos2 = (ex - ey) / (2 * R), sin2 = g / (2 * R);
  const double beta = std::min(1.0, 1.0 / (0.8 + C.kbeta * std::max(e1, 0.0)));
  // Each branch below is computed only where its value is used; the values are the
  // reference's (which evaluates every branch and selects), operation for operation.
  // compression: f2 = -beta f(c), zero in biaxial tension
  double f2 = 0.0;
  if (!(e2 > 0)) {
    const double cc = std::max(-e2, 0.0);
    double fcomp;
    if (H.used && cc < H.em2) {
      const double xx = H.em2 / C.ep;
      const double p2 = std::min(C.ep * (0.166 * (xx * xx) + 0.132 * xx), 0.9 * H.em2);
      double t = (cc - p2) / std::max(H.em2 - p2, 1e-12);
      t = std::min(std::max(t, 0.0), 1.0);
      fcomp = fEnv(C, H.em2) * t;
    }
    else fcomp = fEnv(C, cc);
    f2 = -beta * fcomp;
    if (S.ccyc > 0.0) f2 *= betaCyclic(S, reopenTrial(S, H, e1, e2));   // WP-256, opt-in; off = Model C bit for bit
  }
  // tension, then the crack check 'sum' (06 §2.6); biaxial compression: f1 = -f(-e1)
  const double phi = 0.5 * atan2(sin2, cos2);
  const double th = fabs(0.5 * PI - fabs(phi));
  double f1, w = 0.0;
  if (e1 < 0) f1 = -fEnv(C, std::max(-e1, 0.0));
  else {
    f1 = (H.used && e1 > 0 && e1 < H.em1) ? fTen(C, H.em1) * e1 / std::max(H.em1, 1e-12) : fTen(C, e1);
    if (e1 > 0) {
      const double s_th = 1.0 / (fabs(sin(th)) / C.smx + fabs(cos(th)) / C.smy + 1e-12);
      w = std::max(e1, 0.0) * s_th;
      const double vci = C.cvci * sqrt(C.fc) / (0.31 + 24.0 * w / (C.agg + 16.0));
      const double thc = std::min(std::max(th, 1e-3), 0.5 * PI - 1e-3);
      const double cap = vci * tan(thc) + S.rho * std::max(S.tie.fy - o.fst, 0.0);
      f1 = std::min(f1, cap);
    }
  }
  o.sx = 0.5 * (f1 + f2) + 0.5 * (f1 - f2) * cos2;
  o.sy = 0.5 * (f1 + f2) - 0.5 * (f1 - f2) * cos2;
  o.txy = 0.5 * (f1 - f2) * sin2;
  o.e1 = e1; o.e2 = e2; o.w = w; o.f1 = f1; o.f2 = f2; o.beta = beta; o.theta = th;
}

// transverse residual
inline double resid(const Strip& S, const Hist& H, double ex, double ey, double g, Out& o)
{
  membrane(S, H, ex, ey, g, o);
  return o.sy + S.rho * o.fst;
}

// geomspace(1e-7, 0.15, 40), the scan offsets, exactly as numpy builds them
static const double kScanOffsets[40] = {
  0x1.ad7f29abcaf48p-24, 0x1.353c9f5d0c0ffp-23, 0x1.bd4cbf4d37e0fp-23, 0x1.409d6f79ab10cp-22,
  0x1.cdaf2fc02afe5p-22, 0x1.4c696c6c04d29p-21, 0x1.deabf5227fcb1p-21, 0x1.58a487b98911ep-20,
  0x1.f048bd291c7ebp-20, 0x1.6552d80cb2e73p-19, 0x1.0145b582f1f6cp-18, 0x1.72789a92e9e4cp-18,
  0x1.0abd0cafc45acp-17, 0x1.801a346742d75p-17, 0x1.148d8e2df510dp-16, 0x1.8e3c340a9971bp-16,
  0x1.1eba81ded4afbp-15, 0x1.9ce352e9808f6p-15, 0x1.29474e8ad608cp-14, 0x1.ac1476f08b790p-14,
  0x1.34377b04a3f9cp-13, 0x1.bbd4b42f77403p-13, 0x1.3f8eaf56ed15cp-12, 0x1.cc294e8bc05cfp-12,
  0x1.4b50b5fd4a03fp-11, 0x1.dd17bb8335c79p-11, 0x1.57817d28a733bp-10, 0x1.eea5a3ff20377p-10,
  0x1.6425180f9e4f1p-9, 0x1.006c731c4cd29p-8, 0x1.713fc04b2fe35p-8, 0x1.09dbcbd75b815p-7,
  0x1.7ed5d74051cc4p-7, 0x1.13a4039795ff6p-6, 0x1.8cebe796cb415p-6, 0x1.1dc85f790b704p-5,
  0x1.9b86a6bddbd23p-5, 0x1.284c4364d4762p-4, 0x1.aaaaf67f2f3f4p-4, 0x1.3333333333333p-3,
};

static const int    kNScan = 81;
static const double kBoxLo = -0.01, kBoxHi = 0.15;      // scan box
static const double kFullLo = -0.004, kFullHi = 0.15;   // full-range bracket
static const int    kNBisect = 60;

// Solve eps_y. Returns 1 when the scan bracketed a root near the committed eps_y,
// 0 when it fell back to the full-range bisection (before the first commit, or
// no increasing sign change in the scan). *nEval counts residual evaluations.
inline int solveEy(const Strip& S, const Hist& H, double ex, double g, double& ey, int* nEval = 0)
{
  Out o;
  int ne = 0, found = 0;
  double lo = kFullLo, hi = kFullHi;
  if (H.used) {
    double Y[kNScan], F[kNScan];
    bool have[kNScan];
    for (int i = 0; i < kNScan; ++i) {
      const double cand = (i < 40) ? -kScanOffsets[39 - i] : (i == 40 ? 0.0 : kScanOffsets[i - 41]);
      Y[i] = std::min(std::max(H.ey + cand, kBoxLo), kBoxHi);
      have[i] = false;
    }
    // intervals j = 0..79 ordered by (|mid - ey|, j): the first increasing one is numpy's argmin
    int order[kNScan - 1];
    double dist[kNScan - 1];
    for (int j = 0; j < kNScan - 1; ++j) {
      const double mid = 0.5 * (Y[j] + Y[j + 1]);
      dist[j] = fabs(mid - H.ey);
      order[j] = j;
    }
    std::stable_sort(order, order + kNScan - 1, [&](int a, int b) { return dist[a] < dist[b]; });
    for (int k = 0; k < kNScan - 1; ++k) {
      const int j = order[k];
      for (int i = j; i <= j + 1; ++i)
        if (!have[i]) { F[i] = resid(S, H, ex, Y[i], g, o); have[i] = true; ++ne; }
      if (F[j] <= 0 && F[j + 1] >= 0) { lo = Y[j]; hi = Y[j + 1]; found = 1; break; }
    }
  }
  for (int it = 0; it < kNBisect; ++it) {
    const double m = 0.5 * (lo + hi);
    const double f = resid(S, H, ex, m, g, o);
    ++ne;
    if (f <= 0) lo = m; else hi = m;
  }
  ey = 0.5 * (lo + hi);
  if (nEval) *nEval = ne;
  return found;
}

// Full strip state at (ex, g): solve eps_y, then the membrane there.
inline int evaluate(const Strip& S, const Hist& H, double ex, double g, double& ey, Out& o, int* nEval = 0)
{
  const int found = solveEy(S, H, ex, g, ey, nEval);
  membrane(S, H, ex, ey, g, o);
  return found;
}

// Commit (06 §2.8): o is the converged trial state at eps_y = ey.
inline void commit(const Strip& S, Hist& H, double ey, const Out& o)
{
  const double emPrev = H.em1, em2Prev = H.em2;   // maxima before this commit (WP-256)
  if (H.used) {
    H.em1 = std::max(H.em1, o.e1);
    H.em2 = std::max(H.em2, std::max(-o.e2, 0.0));
  } else {
    H.em1 = std::max(0.0, o.e1);
    H.em2 = std::max(0.0, -o.e2);
  }
  if (S.rho != 0.0) H.ept = ey - steelStress(S.tie, ey, H.ept) / S.tie.Es;
  // WP-256: store the re-opening the converged trial used (reopenTrial with the memory
  // before this commit), and the opening it reached.
  {
    Hist Hprev = H; Hprev.em1 = emPrev; Hprev.em2 = em2Prev;
    H.reopen = reopenTrial(S, Hprev, o.e1, o.e2);
    H.e1c = std::max(o.e1 - S.conc.ecr, 0.0);
    H.e2c = std::max(-o.e2 - 0.5 * S.conc.ep, 0.0);
  }
  H.ey = ey;
  H.used = true;
}

// Condensed tangent d(sx, txy)/d(ex, g) by forward differences (06 §1.5): each
// perturbed state re-solves eps_y from the same committed memory.
inline void tangentFD(const Strip& S, const Hist& H, double ex, double g, const Out& o0,
                      double hEx, double hG, double D[2][2])
{
  double ey; Out o;
  evaluate(S, H, ex + hEx, g, ey, o);
  D[0][0] = (o.sx - o0.sx) / hEx;  D[1][0] = (o.txy - o0.txy) / hEx;
  evaluate(S, H, ex, g + hG, ey, o);
  D[0][1] = (o.sx - o0.sx) / hG;   D[1][1] = (o.txy - o0.txy) / hG;
}

}  // namespace ladruno_mcft

#endif
