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

// Ladruno WP-257 (program #1131, improvement I3): dowel action of a longitudinal bar.
//
// Header-only and dependency-free, so tests/_testbed/dowel_gpp.cpp checks it with plain g++
// (closed forms, finite differences, paths) without an OpenSees build.
//
// A bar fibre of a LadrunoNDFiber section (LadrunoUniaxialFiber -dowel) carries, besides its
// axial stress, a shear stress per shear component k (gxy; and gxz in the 3D view) from the
// slip of the crack it crosses. The fibre's shear strain is turned into a slip over a slip
// length ls (the crack spacing: a smeared crack of spacing ls slips s = gamma ls):
//
//     s_k   = |gamma_k| ls
//     tau_k = sign(gamma_k) cap(sigma) h(s_k, smax_k)
//
// Capacity (Dulacska 1972; Vintzeleou & Tassios 1987): the dowel force of a bar of diameter db
// is V_du = kd db^2 sqrt(fc fy (1 - alpha^2)), alpha = sigma/fy the bar's axial stress ratio
// (kd = 1.3 for a bar in a massive concrete core; smaller where cover splitting governs). Per
// unit bar area, A_b = pi db^2 / 4:
//
//     cap(sigma) = kd (4/pi) sqrt(fc fy) sqrt(1 - alpha^2),   alpha = min(|sigma|/fy, 1)
//
// so a yielded bar carries no dowel force.
//
// Slip law: the envelope h(s) = 1 - exp(-s/s0) with s0 = su/3 (95 % of the capacity at the
// slip su, default 0.05 db, Vintzeleou & Tassios' order of magnitude). Inside the largest
// committed slip smax_k the fibre unloads and reloads on the secant to the origin,
// h = h(smax) s/smax: the pinched, origin-oriented dowel loop. Only smax is history.
//
// Cyclic decay (opt-in, -dowelDecay Cd; Vintzeleou & Tassios: dowels lose strength under slip
// reversals): one committed memory D_k, the cumulative RE-slip (the part of each increase of the
// slip that stays below its previous maximum smax_k, in units of su; a monotonic path adds exactly
// zero), and tau_k carries beta_d = 1/(1 + Cd D_k). Implicit, as LadrunoMCFTStrip -cyclicDegradation:
// the trial adds its own re-slip since the last commit (sc_k), so a committed state re-solves to
// itself. Cd = 0 (default) leaves every value untouched.
//
// Tangent (exact): d tau_k / d gamma_k = cap (h' beta + h beta') ls; d tau_k / d eps = sign h beta
// dcap/dsigma E_t.

#ifndef LadrunoDowelKernel_h
#define LadrunoDowelKernel_h

#include <algorithm>
#include <cmath>

namespace ladruno_dowel {

struct Params {
  double fc;    // concrete strength (> 0)
  double fy;    // bar yield stress (> 0)
  double db;    // bar diameter (> 0)
  double ls;    // slip length: s = |gamma| ls (> 0)
  double kd;    // capacity factor (Dulacska / Vintzeleou-Tassios 1.3)
  double su;    // slip at 95 % of the capacity (default 0.05 db)
  double cd;    // cyclic decay rate per su of re-slip (0 = off)
  Params() : fc(0.0), fy(0.0), db(0.0), ls(0.0), kd(1.3), su(-1.0), cd(0.0) {}
  double slipU() const { return su > 0.0 ? su : 0.05 * db; }
  double s0() const { return slipU() / 3.0; }
  double tauU() const { return kd * (4.0 / 3.14159265358979323846) * std::sqrt(fc * fy); }
  bool valid() const { return fc > 0.0 && fy > 0.0 && db > 0.0 && ls > 0.0 && kd > 0.0; }
};

// capacity per unit bar area at the axial stress sig; dcap = d cap / d sig
inline double cap(const Params& P, double sig, double* dcap)
{
  const double a = std::fabs(sig) / P.fy;
  if (!(a < 1.0)) { if (dcap) *dcap = 0.0; return 0.0; }
  const double r = std::sqrt(1.0 - a * a), t = P.tauU();
  if (dcap) *dcap = -t * a / r * (sig >= 0.0 ? 1.0 : -1.0) / P.fy;
  return t * r;
}

// normalised slip response h in [0, 1) and dh/ds, given the committed maximum slip smax
inline double shape(const Params& P, double s, double smax, double* dh)
{
  const double s0 = P.s0();
  if (s >= smax) {                              // envelope
    const double x = std::exp(-s / s0);
    if (dh) *dh = x / s0;
    return 1.0 - x;
  }
  const double hm = 1.0 - std::exp(-smax / s0); // secant to the origin (smax > s >= 0)
  if (dh) *dh = hm / smax;
  return hm * s / smax;
}

// one shear component: stress, d tau / d gamma, d tau / d sigma
struct Out { double tau, dg, dsig; };

// the cumulative re-slip D (in su) at slip s, from the committed (smax, sc, D); dD = dD/ds
inline double reslip(const Params& P, double s, double smax, double sc, double D, double* dD)
{
  const double up = std::min(s, smax) - sc;
  if (up > 0.0 && s < smax) { if (dD) *dD = 1.0 / P.slipU(); return D + up / P.slipU(); }
  if (dD) *dD = 0.0;
  return (up > 0.0) ? D + up / P.slipU() : D;
}

inline Out shear(const Params& P, double gamma, double sig, double smax, double sc = 0.0, double D = 0.0)
{
  Out o;
  double dc, dh;
  const double c = cap(P, sig, &dc);
  const double s = std::fabs(gamma) * P.ls;
  const double h = shape(P, s, smax, &dh);
  const double sg = (gamma > 0.0) ? 1.0 : (gamma < 0.0 ? -1.0 : 0.0);
  if (P.cd > 0.0) {                             // cyclic decay (opt-in)
    double dD;
    const double Dt = reslip(P, s, smax, sc, D, &dD);
    const double b = 1.0 / (1.0 + P.cd * Dt), db = -P.cd * b * b * dD;
    o.tau = sg * c * h * b;
    o.dg = c * (dh * b + h * db) * P.ls;
    o.dsig = sg * h * b * dc;
    return o;
  }
  o.tau = sg * c * h;
  o.dg = c * dh * P.ls;                         // even in gamma: no sign
  o.dsig = sg * h * dc;
  return o;
}

// the committed re-slip after a converged state at gamma (before smax is updated)
inline double commitReslip(const Params& P, double gamma, double smax, double sc, double D)
{
  return reslip(P, std::fabs(gamma) * P.ls, smax, sc, D, 0);
}

// the committed maximum slip after a converged state at gamma
inline double commitSlip(const Params& P, double gamma, double smax)
{
  const double s = std::fabs(gamma) * P.ls;
  return s > smax ? s : smax;
}

}  // namespace ladruno_dowel

#endif
