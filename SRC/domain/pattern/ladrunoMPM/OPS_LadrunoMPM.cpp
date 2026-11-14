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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 11/2026
//
// OPS_LadrunoMPM — the `mpm` command (WP-232 slice 1, mode F). One row of
// SRC/interpreter/LadrunoCommandTable.h; the shared entry reads its arguments
// through the elementAPI cursor, so the same code serves Python (opensees.pyd)
// and both Tcl engines (classic Tcl bridges here through
// OPS_ResetInputNoBuilder, the pattern of the WP-127 replay command). It needs
// no engine state: mode F has no Analysis object.
//
//   mpm grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp|linear]
//   mpm body <bodyId> box <lo..> <hi..> -ppc <n> -mat <matTag> -rho <rho>
//   mpm body <bodyId> disc <cx> <cy> <R> -ppc ... (2D) | ball <cx> <cy> <cz> <R> ... (3D)
//   mpm particles -file <path>           load a Contract-1 particle set (232a), path as given
//   mpm units <L> <F> <T> [-system <id>]  units Contract 2 carries when no set declared them
//   mpm output <path> [-every N | -dt dt] [-fields F..] [-compress k] [-flush s] | -close
//                                         the .ladruno particle group (232a contract 2)
//   mpm shift <body> <dx..>              translate a seeded body (before the first step)
//   mpm particle <id> <body> <matTag> <rho> <V0> -x <..> -lp <..> [-v <..>] [-s <nv>]
//   mpm velocity <body> uniform <v..> | sine <comp> <axis> <v0> <beta> <xref>
//   mpm stress <body> uniform <nv values> | geostatic <yTop> <rho> <g> <K0>   (2D)
//   mpm fix <face> <components>          face xmin..zmax, components a string over xyz
//   mpm wall <face> <mu>                  stationary Coulomb wall
//   mpm gravity <g..>
//   mpm timeSeries <tsTag>                load factor of the body force (default 1)
//   mpm damping local <alpha>
//   mpm flip <alphaF> | mpm massCutoff <c> | mpm safety <s>
//   mpm criticalTimeStep                  -> dt_c
//   mpm analyze <nSteps> <dt>             -> 0
//   mpm reset
//   mpm info <numParticles|epoch|time|summary>  summary: [mass, momentum.., max |v|, KE]
//   mpm particleData <field>              -> flat list over particles in storage order;
//        field: id body mat m V0 V x v a lp F R sig sigHat epsHat sig0
//   mpm debug shape <x..> -lp <lp..>      -> [node S dS..]* (interface Q4)
//   mpm debug shape1d <h> <lp> <r>        -> S dS
//   mpm debug activeNodes                 -> sorted active node ids at the committed state
//   mpm debug injectRefusal | staleRefusal | aborted   W17 test hooks (a refusal at the next
//                                         particle commit; a stale WP-99 count; the abort latch)
//
// Mode A (WP-232f, tests/mpm_implicit_interface.md §5):
//   mpm mode explicit|implicit|quasiStatic
//   mpm gridvel <name> box <lo..> <hi..> v <value|free ..> [start <t0>] [track <axis>]
//   mpm condForce <name>                  -> ndm doubles at the last accepted step
//   mpm pforce <body> box <lo..> <hi..> <f..>   -> count
//   mpm emptyNodes exclude | renormalize | ghost <gamma>
//   mpm cut <factor> <maxCuts>  |  mpm tagBase <nodeBase> <eleBase>  |  mpm rayleigh <alphaM> <betaK>
//   mpm info g2p|acceptedSteps|retries|domainChanges|connectivityChanges|gridNodes|mode
//   mpm info dtHistory | gridNodeTags     -> lists
//   mpm debug failAttempt <k> [<attempts>] | reseed default|always|skip
// `mpm analyze` in mode A needs the user's analysis: the engine bridges hand it
// in through OPS_LadrunoMPMWithAnalysis (classic Tcl: commands.cpp; Python /
// DL Tcl: OPS_LadrunoMPMDL in OpenSeesCommands.cpp).
//
// Refused in OpenSeesSP / OpenSeesMP (ADR 232 §5.1, v1).

#include "LadrunoMPMPattern.h"

#include <Domain.h>
#include <LoadPattern.h>
#include <TimeSeries.h>
#include <elementAPI.h>
#include <LadrunoThreads.h>
#include <OPS_Globals.h>

#include <algorithm>
#include <cmath>
#include <string.h>
#include <string>
#include <vector>

extern TimeSeries *OPS_getTimeSeries(int tag);

namespace {

// The analysis of the engine whose bridge is running this command (WP-232f).
const LadrunoMPMAnalysisRefs *g_refs = 0;

int fail(const std::string &msg)
{
  opserr << "WARNING mpm " << msg.c_str() << "\n";
  return -1;
}

bool readDoubles(int n, double *out)
{
  if (OPS_GetNumRemainingInputArgs() < n) return false;
  int num = n;
  return OPS_GetDoubleInput(&num, out) == 0;
}

bool readInt(int &out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  int num = 1;
  return OPS_GetIntInput(&num, &out) == 0;
}

int outDoubles(std::vector<double> &v)
{
  int n = (int)v.size();
  if (n == 0) {
    double z = 0.0;
    int zero = 0;
    return OPS_SetDoubleOutput(&zero, &z, false);
  }
  return OPS_SetDoubleOutput(&n, v.data(), false);
}

int outDouble(double x)
{
  int one = 1;
  return OPS_SetDoubleOutput(&one, &x, true);
}

int outInt(int x)
{
  int one = 1;
  return OPS_SetIntOutput(&one, &x, true);
}

LadrunoMPMPattern *host(Domain *dom, const char *sub)
{
  LadrunoMPMPattern *h = ladrunoMPMFindHost(dom);
  if (h == 0)
    opserr << "WARNING mpm " << sub << " -- no MPM grid: run `mpm grid <tag> ...` first\n";
  return h;
}

// mpm grid <tag> -origin .. -h .. -cells .. [-basis gimp|linear]
int cmdGrid(Domain *dom)
{
  int tag;
  if (!readInt(tag)) return fail("grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp|linear]");
  if (ladrunoMPMFindHost(dom) != 0) return fail("grid -- this Domain already has an MPM grid (v1: one per Domain)");
  if (dom->getLoadPattern(tag) != 0) return fail("grid -- a load pattern with tag " + std::to_string(tag) + " exists");
  const int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) return fail("grid -- the model must be ndm 2 (plane strain) or 3");
  double x0[3] = {0, 0, 0}, h = 0.0;
  int n[3] = {0, 0, 0};
  bool haveO = false, haveH = false, haveN = false;
  ladruno_mpm::Basis basis = ladruno_mpm::BASIS_GIMP;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-origin") == 0) {
      if (!readDoubles(ndm, x0)) return fail("grid -origin needs " + std::to_string(ndm) + " values");
      haveO = true;
    } else if (strcmp(opt, "-h") == 0) {
      if (!readDoubles(1, &h)) return fail("grid -h needs a value");
      haveH = true;
    } else if (strcmp(opt, "-cells") == 0) {
      int num = ndm;
      if (OPS_GetNumRemainingInputArgs() < ndm || OPS_GetIntInput(&num, n) != 0)
        return fail("grid -cells needs " + std::to_string(ndm) + " integers");
      haveN = true;
    } else if (strcmp(opt, "-basis") == 0) {
      const char *b = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "";
      if (strcmp(b, "gimp") == 0) basis = ladruno_mpm::BASIS_GIMP;
      else if (strcmp(b, "linear") == 0) basis = ladruno_mpm::BASIS_LINEAR;
      else return fail("grid -basis must be gimp or linear");
    } else {
      return fail(std::string("grid -- unknown option '") + opt + "'");
    }
  }
  if (!haveO || !haveH || !haveN) return fail("grid needs -origin, -h and -cells");
  LadrunoMPMPattern *p = new LadrunoMPMPattern(tag);
  std::string err;
  if (!p->engine().setGrid(ndm, x0, h, n, err)) {
    delete p;
    return fail("grid -- " + err);
  }
  p->engine().par.basis = basis;
  if (!dom->addLoadPattern(p)) {
    delete p;
    return fail("grid -- could not add the MPM host to the Domain");
  }
  return 0;
}

// -ppc n -mat tag -rho rho
bool readBodyOptions(int &ppc, int &mat, double &rho)
{
  ppc = 0; mat = 0; rho = 0.0;
  bool a = false, b = false, c = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-ppc") == 0) a = readInt(ppc);
    else if (strcmp(opt, "-mat") == 0) b = readInt(mat);
    else if (strcmp(opt, "-rho") == 0) c = readDoubles(1, &rho);
    else { opserr << "WARNING mpm body -- unknown option '" << opt << "'\n"; return false; }
  }
  return a && b && c;
}

int cmdBody(LadrunoMPMPattern *h)
{
  const int nd = h->engine().ndm();
  int body;
  if (!readInt(body) || OPS_GetNumRemainingInputArgs() < 1)
    return fail("body <bodyId> box|disc|ball <geometry> -ppc <n> -mat <matTag> -rho <rho>");
  const char *shp = OPS_GetString();
  int shape, ng;
  if (strcmp(shp, "box") == 0) { shape = 0; ng = 2 * nd; }
  else if (strcmp(shp, "disc") == 0) { shape = 1; ng = 3; }
  else if (strcmp(shp, "ball") == 0) { shape = 2; ng = 4; }
  else return fail(std::string("body -- unknown shape '") + shp + "' (box, disc, ball)");
  double geom[6];
  if (!readDoubles(ng, geom)) return fail("body -- " + std::to_string(ng) + " geometry values expected");
  int ppc, mat;
  double rho;
  if (!readBodyOptions(ppc, mat, rho)) return fail("body needs -ppc <n> -mat <matTag> -rho <rho>");
  std::string err;
  const int n = h->seedBody(shape, geom, ppc, body, mat, rho, err);
  if (n < 0) return fail("body -- " + err);
  return outInt(n);
}

// particle <id> <body> <mat> <rho> <V0> -x .. -lp .. [-v ..] [-s ..]
int cmdParticle(LadrunoMPMPattern *h)
{
  const int nd = h->engine().ndm(), nv = h->engine().nv();
  double idd;
  int body, mat;
  double rv[2];
  if (!readDoubles(1, &idd) || !readInt(body) || !readInt(mat) || !readDoubles(2, rv))
    return fail("particle <id> <body> <matTag> <rho> <V0> -x <..> -lp <..> [-v <..>] [-s <..>]");
  if (idd != std::floor(idd)) return fail("particle -- the id must be an integer");
  double x[3], lp[3], v[3] = {0, 0, 0}, s[6] = {0, 0, 0, 0, 0, 0};
  bool hx = false, hl = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-x") == 0) hx = readDoubles(nd, x);
    else if (strcmp(opt, "-lp") == 0) hl = readDoubles(nd, lp);
    else if (strcmp(opt, "-v") == 0) { if (!readDoubles(nd, v)) return fail("particle -v needs ndm values"); }
    else if (strcmp(opt, "-s") == 0) { if (!readDoubles(nv, s)) return fail("particle -s needs " + std::to_string(nv) + " values"); }
    else return fail(std::string("particle -- unknown option '") + opt + "'");
  }
  if (!hx || !hl) return fail("particle needs -x and -lp");
  std::string err;
  if (h->addParticle((int64_t)idd, body, mat, rv[0], rv[1], x, lp, v, s, err) < 0)
    return fail("particle -- " + err);
  return 0;
}

int cmdVelocity(LadrunoMPMPattern *h)
{
  const int nd = h->engine().ndm();
  int body;
  if (!readInt(body) || OPS_GetNumRemainingInputArgs() < 1)
    return fail("velocity <body> uniform <v..> | sine <comp> <axis> <v0> <beta> <xref>");
  const char *kind = OPS_GetString();
  std::string err;
  int n;
  if (strcmp(kind, "uniform") == 0) {
    double v[3];
    if (!readDoubles(nd, v)) return fail("velocity uniform needs ndm values");
    n = h->setBodyVelocity(body, v, err);
  } else if (strcmp(kind, "sine") == 0) {
    if (OPS_GetNumRemainingInputArgs() < 5) return fail("velocity sine <comp> <axis> <v0> <beta> <xref>");
    std::string comp = OPS_GetString(), axis = OPS_GetString();
    double a[3];
    if (!readDoubles(3, a)) return fail("velocity sine <comp> <axis> <v0> <beta> <xref>");
    const int c = comp.empty() ? -1 : comp[0] - 'x', ax = axis.empty() ? -1 : axis[0] - 'x';
    n = h->setBodyVelocitySine(body, c, ax, a[0], a[1], a[2], err);
  } else {
    return fail(std::string("velocity -- unknown kind '") + kind + "'");
  }
  if (n < 0) return fail("velocity -- " + err);
  return outInt(n);
}

int cmdStress(LadrunoMPMPattern *h)
{
  const int nv = h->engine().nv();
  int body;
  if (!readInt(body) || OPS_GetNumRemainingInputArgs() < 1)
    return fail("stress <body> uniform <nv values> | geostatic <yTop> <rho> <g> <K0>");
  const char *kind = OPS_GetString();
  std::string err;
  int n;
  if (strcmp(kind, "uniform") == 0) {
    std::vector<double> s(nv);
    if (!readDoubles(nv, s.data())) return fail("stress uniform needs " + std::to_string(nv) + " values");
    n = h->setBodyStress(body, s, err);
  } else if (strcmp(kind, "geostatic") == 0) {
    double a[4];
    if (!readDoubles(4, a)) return fail("stress geostatic <yTop> <rho> <g> <K0>");
    n = h->setBodyStressGeostatic(body, a[0], a[1], a[2], a[3], err);
  } else {
    return fail(std::string("stress -- unknown kind '") + kind + "'");
  }
  if (n < 0) return fail("stress -- " + err);
  return outInt(n);
}

int cmdFix(LadrunoMPMPattern *h, bool wall)
{
  const int nd = h->engine().ndm();
  if (OPS_GetNumRemainingInputArgs() < 2)
    return fail(wall ? "wall <face> <mu>" : "fix <face> <components>");
  const int f = ladruno_mpm::faceIndex(OPS_GetString());
  if (f < 0 || f >= 2 * nd) return fail("-- the face must be xmin, xmax, ymin, ymax[, zmin, zmax]");
  if (wall) {
    double mu;
    if (!readDoubles(1, &mu) || mu < 0.0) return fail("wall -- mu must be >= 0");
    h->engine().fc.wallMu[f] = mu;
    return 0;
  }
  const std::string comps = OPS_GetString();
  unsigned mask = 0;
  for (char c : comps) {
    const int a = c - 'x';
    if (a < 0 || a >= nd) return fail("fix -- components are a string over x, y[, z]");
    mask |= 1u << a;
  }
  h->engine().fc.fixMask[f] |= mask;
  return 0;
}

int particleField(const ladruno_mpm::Particles &P, const std::string &f, std::vector<double> &out)
{
  out.clear();
  if (f == "id") for (int64_t i : P.id) out.push_back((double)i);
  else if (f == "body") for (int b : P.body) out.push_back(b);
  else if (f == "mat") for (int m : P.mat) out.push_back(m);
  else if (f == "m") out = P.m;
  else if (f == "V0") out = P.V0;
  else if (f == "V") out = P.V;
  else if (f == "x") out = P.x;
  else if (f == "v") out = P.v;
  else if (f == "a") out = P.a;
  else if (f == "lp") out = P.lp;
  else if (f == "F") out = P.F;
  else if (f == "R") out = P.R;
  else if (f == "sig") out = P.sig;
  else if (f == "sigHat") out = P.sigHat;
  else if (f == "epsHat") out = P.eps;
  else if (f == "sig0") out = P.sig0;
  else return fail("particleData -- unknown field '" + f + "'");
  return 0;
}

// [total mass, momentum (ndm), max |v|, kinetic energy] of the committed state
void summary(const ladruno_mpm::Particles &P, int nd, std::vector<double> &out)
{
  double mt = 0.0, mom[3] = {0, 0, 0}, vmax = 0.0, ke = 0.0;
  for (int p = 0; p < P.size(); p++) {
    double v2 = 0.0;
    for (int d = 0; d < nd; d++) {
      mom[d] += P.m[p] * P.v[nd * p + d];
      v2 += P.v[nd * p + d] * P.v[nd * p + d];
    }
    mt += P.m[p];
    ke += 0.5 * P.m[p] * v2;
    vmax = std::max(vmax, std::sqrt(v2));
  }
  out.assign(1, mt);
  for (int d = 0; d < nd; d++) out.push_back(mom[d]);
  out.push_back(vmax);
  out.push_back(ke);
}

// gridvel <name> box <lo..> <hi..> v <value|free ..> [start <t0>] [track <axis>]
int cmdGridVel(LadrunoMPMPattern *h)
{
  const int nd = h->engine().ndm();
  const char *usage = "gridvel <name> box <lo..> <hi..> v <value|free ..> [start <t0>] [track x|y|z]";
  if (OPS_GetNumRemainingInputArgs() < 2) return fail(usage);
  const std::string name = OPS_GetString();
  if (std::string(OPS_GetString()) != "box") return fail(usage);
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0}, v[3] = {0, 0, 0};
  if (!readDoubles(nd, lo) || !readDoubles(nd, hi)) return fail(usage);
  if (OPS_GetNumRemainingInputArgs() < nd + 1 || std::string(OPS_GetString()) != "v") return fail(usage);
  unsigned freeMask = 0;
  for (int d = 0; d < nd; d++) {
    // A value is 'free' or a number. Peek as a string, then step back and read a
    // double: in Python a float token reads as the literal "Invalid String Input!"
    // (element guide, Parser), so parsing the peeked text would refuse it.
    const std::string tok = OPS_GetString();
    if (tok == "free") { freeMask |= 1u << d; continue; }
    OPS_ResetCurrentInputArg(-1);
    if (!readDoubles(1, &v[d]))
      return fail("gridvel -- component " + std::to_string(d) + " is not a number or 'free'");
  }
  double start = -1e300;
  int track = -1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const std::string opt = OPS_GetString();
    if (opt == "start") {
      if (!readDoubles(1, &start)) return fail(usage);
    } else if (opt == "track") {
      const std::string ax = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "";
      track = ax.empty() ? -1 : ax[0] - 'x';
      if (ax.size() != 1 || track < 0 || track >= nd) return fail("gridvel track x|y[|z]");
    } else {
      return fail("gridvel -- unknown option '" + opt + "'");
    }
  }
  std::string err;
  if (h->addGridVel(name, lo, hi, v, freeMask, start, track, err) < 0) return fail("gridvel -- " + err);
  return 0;
}

int cmdDebug(LadrunoMPMPattern *h)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return fail("debug shape|shape1d|activeNodes");
  const std::string what = OPS_GetString();
  ladruno_mpm::Engine &E = h->engine();
  const int nd = E.ndm();
  if (what == "shape1d") {
    double a[3];
    if (!readDoubles(3, a)) return fail("debug shape1d <h> <lp> <r>");
    double S, dS;
    ladruno_mpm::gimpFactor(a[2], a[0], a[1], S, dS);
    std::vector<double> out = {S, dS};
    return outDoubles(out);
  }
  if (what == "shape") {
    double x[3], lp[3] = {0, 0, 0};
    if (!readDoubles(nd, x)) return fail("debug shape <x..> -lp <lp..>");
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char *opt = OPS_GetString();
      if (strcmp(opt, "-lp") != 0 || !readDoubles(nd, lp)) return fail("debug shape <x..> -lp <lp..>");
    }
    std::vector<long> node;
    std::vector<double> S, dS, out;
    ladruno_mpm::particleSupport(E.grid, E.par.basis, x, lp, node, S, dS);
    for (size_t k = 0; k < node.size(); k++) {
      out.push_back((double)E.grid.publicId(node[k]));
      out.push_back(S[k]);
      for (int d = 0; d < nd; d++) out.push_back(dS[nd * k + d]);
    }
    return outDoubles(out);
  }
  if (what == "injectRefusal") { h->injectCommitRefusal(); return 0; }   // W17 test hook
  if (what == "staleRefusal") { h->injectStaleRefusal(); return 0; }     // W17 test hook
  if (what == "aborted") return outInt(h->isAborted() ? 1 : 0);
  if (what == "failAttempt") {                                         // WP-232f retry test hook
    int k, m = 1;
    if (!readInt(k)) return fail("debug failAttempt <step> [<attempts>]");
    if (OPS_GetNumRemainingInputArgs() > 0 && !readInt(m)) return fail("debug failAttempt <step> [<attempts>]");
    h->setDebugFail(k, m);
    return 0;
  }
  if (what == "tangentCheck") {                                        // WP-232f footing diagnosis
    double a[2];
    if (!readDoubles(2, a)) return fail("debug tangentCheck <scale> <eps>");
    std::vector<double> out;
    std::string err;
    if (h->debugTangentCheck(a[0], a[1], out, err) < 0) return fail("debug tangentCheck -- " + err);
    return outDoubles(out);
  }
  if (what == "reseed") {                                              // WP-232f re-seed test hook
    const std::string m = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "";
    if (m == "default") h->setDebugReseed(0);
    else if (m == "always") h->setDebugReseed(1);
    else if (m == "skip") h->setDebugReseed(2);
    else return fail("debug reseed default|always|skip");
    return 0;
  }
  if (what == "activeNodes") {
    std::vector<long> act;
    E.activeNodes(act);
    std::vector<double> out(act.begin(), act.end());
    return outDoubles(out);
  }
  return fail("debug -- unknown query '" + what + "'");
}

}  // namespace

// The double-valued queries (`mpm criticalTimeStep`, `mpm particleData <f>`,
// `mpm info time|summary`), shared by the shared entry and the classic-Tcl
// bridge. The bridge publishes them as Tcl double objects (full precision),
// because classic Tcl's OPS_SetDoubleOutput prints %35.20f, which does not
// round-trip below about 1e-3 (review Q8). Returns 1 (handled, `out` filled),
// 0 (not one of these queries), or -1 (an error, printed).
int ladrunoMPMDoubleQuery(Domain *dom, const std::vector<std::string> &a, std::vector<double> &out)
{
  out.clear();
  if (a.empty()) return 0;
  const std::string &sub = a[0];
  const bool q = sub == "criticalTimeStep" || sub == "particleData" || sub == "condForce" ||
                 (sub == "info" && a.size() > 1 && (a[1] == "time" || a[1] == "summary" ||
                                                     a[1] == "dtHistory" || a[1] == "gridNodeTags"));
  if (!q) return 0;
  if (ladruno_parallelBuild())
    return fail("-- refused in OpenSeesSP / OpenSeesMP: v1 MPM is serial (ADR 232 §5.1)");
  LadrunoMPMPattern *h = host(dom, sub.c_str());
  if (h == 0) return -1;
  ladruno_mpm::Engine &E = h->engine();
  if (sub == "criticalTimeStep") { out.assign(1, h->criticalTimeStep()); return 1; }
  if (sub == "particleData") {
    if (a.size() < 2) return fail("particleData <field>");
    return particleField(E.P, a[1], out) < 0 ? -1 : 1;
  }
  if (sub == "condForce") {
    if (a.size() < 2) return fail("condForce <name>");
    std::string err;
    if (h->condForce(a[1], out, err) < 0) return fail("condForce -- " + err);
    return 1;
  }
  if (a[1] == "time") { out.assign(1, dom->getCurrentTime()); return 1; }
  if (a[1] == "dtHistory") { out = h->dtHistory(); return 1; }
  if (a[1] == "gridNodeTags") { h->gridNodeTags(out); return 1; }
  summary(E.P, E.ndm(), out);
  return 1;
}

int OPS_LadrunoMPM(void);

// The engine bridges (WP-232f, interface §4): classic Tcl (commands.cpp) and
// Python / DL Tcl (OPS_LadrunoMPMDL, OpenSeesCommands.cpp) hand in the user's
// analysis; `mpm analyze` in mode A uses it.
int OPS_LadrunoMPMWithAnalysis(StaticAnalysis *sa, DirectIntegrationAnalysis *ta,
                               StaticIntegrator *si, TransientIntegrator *ti,
                               ConstraintHandler *handler)
{
  LadrunoMPMAnalysisRefs r;
  r.sa = sa;
  r.ta = ta;
  r.si = si;
  r.ti = ti;
  r.handler = handler;
  r.bridged = true;
  const LadrunoMPMAnalysisRefs *saved = g_refs;
  g_refs = &r;
  const int rc = OPS_LadrunoMPM();
  g_refs = saved;
  return rc;
}

int OPS_LadrunoMPM(void)
{
  if (ladruno_parallelBuild())
    return fail("-- refused in OpenSeesSP / OpenSeesMP: v1 MPM is serial (ADR 232 §5.1)");
  Domain *dom = OPS_GetDomain();
  if (dom == 0) return fail("-- no Domain");
  if (OPS_GetNumRemainingInputArgs() < 1)
    return fail("<grid|body|particles|units|output|shift|particle|velocity|stress|fix|wall|gravity|timeSeries|damping|flip|"
                "massCutoff|safety|criticalTimeStep|analyze|reset|info|particleData|debug|mode|gridvel|"
                "condForce|pforce|emptyNodes|cut|tagBase|rayleigh> ...");
  const std::string sub = OPS_GetString();
  if (sub == "grid") return cmdGrid(dom);

  LadrunoMPMPattern *h = host(dom, sub.c_str());
  if (h == 0) return -1;
  ladruno_mpm::Engine &E = h->engine();
  const int nd = E.ndm();

  if (sub == "body") return cmdBody(h);
  if (sub == "particle") return cmdParticle(h);
  if (sub == "particles") {
    // mpm particles -file <path>: a Contract-1 set (232a §1), opened as given
    if (OPS_GetNumRemainingInputArgs() < 2 || std::string(OPS_GetString()) != "-file")
      return fail("particles -file <path>");
    const std::string path = OPS_GetString();
    const int n = h->loadParticleSet(path);
    if (n < 0) return -1;
    return outInt(n);
  }
  if (sub == "units") {
    // mpm units <L> <F> <T> [-system <id>]
    if (OPS_GetNumRemainingInputArgs() < 3) return fail("units <L> <F> <T> [-system <id>]");
    const std::string L = OPS_GetString(), F = OPS_GetString(), T = OPS_GetString();
    std::string system;
    if (OPS_GetNumRemainingInputArgs() >= 2 && std::string(OPS_GetString()) == "-system")
      system = OPS_GetString();
    std::string err;
    if (h->setUnits(L, F, T, system, err) < 0) return fail("units -- " + err);
    return 0;
  }
  if (sub == "output") {
    // mpm output <path> [-every <N> | -dt <dt>] [-fields <F>...] [-compress <k>]
    // mpm output -close
    if (OPS_GetNumRemainingInputArgs() < 1) return fail("output <path> [-every N | -dt dt] [-fields ...] [-compress k] | -close");
    const std::string first = OPS_GetString();
    if (first == "-close") { h->closeOutput(); return 0; }
    ladruno_mpm_io::WriterConfig cfg;
    cfg.path = first;
    std::vector<std::string> fields;
    bool explicitFields = false;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      const std::string opt = OPS_GetString();
      if (opt == "-every") {
        if (!readInt(cfg.every) || cfg.every < 1) return fail("output -every needs N >= 1");
        cfg.dtOut = -1.0;
      } else if (opt == "-dt") {
        if (!readDoubles(1, &cfg.dtOut) || !(cfg.dtOut > 0.0)) return fail("output -dt needs dt > 0");
      } else if (opt == "-flush") {
        if (!readDoubles(1, &cfg.flushSeconds)) return fail("output -flush needs seconds (<= 0, the default: after every output)");
      } else if (opt == "-compress") {
        if (!readInt(cfg.compress) || cfg.compress < 0 || cfg.compress > 9) return fail("output -compress needs 0..9");
      } else if (opt == "-fields") {
        explicitFields = true;
        while (OPS_GetNumRemainingInputArgs() > 0) {
          const std::string f = OPS_GetString();
          if (!f.empty() && f[0] == '-') { OPS_ResetCurrentInputArg(-1); break; }
          if (f.compare(0, 6, "STATE_") == 0)
            return fail("output -fields " + f + " -- v1 admits no STATE_ token (the guide lists the admitted ones)");
          if (ladruno_mpm_io::ParticleWriter::components(f, nd).empty())
            return fail("output -fields " + f + " -- not a v1 particle field (POSITION VELOCITY STRESS "
                        "STRAIN VOLUME GIMP_HALF_LENGTH DEFORMATION_GRADIENT)");
          fields.push_back(f);
        }
      } else {
        return fail("output -- unknown option '" + opt + "'");
      }
    }
    if (!explicitFields) fields = {"VELOCITY", "STRESS", "VOLUME"};
    cfg.fields.push_back("POSITION");          // always written, always first
    for (const std::string &f : fields)
      if (std::find(cfg.fields.begin(), cfg.fields.end(), f) == cfg.fields.end())
        cfg.fields.push_back(f);
    std::string err;
    if (h->setOutput(cfg, err) < 0) return fail("output -- " + err);
    return 0;
  }
  if (sub == "shift") {
    // mpm shift <body> <dx..>: translate a seeded body before the first step (A5)
    int body;
    double dx[3] = {0, 0, 0};
    if (!readInt(body) || !readDoubles(nd, dx)) return fail("shift <body> <dx..>");
    if (h->hasStepped()) return fail("shift -- particles cannot move by hand after the first step");
    return outInt(E.shiftBody(body, dx));
  }
  if (sub == "velocity") return cmdVelocity(h);
  if (sub == "stress") return cmdStress(h);
  if (sub == "fix") return cmdFix(h, false);
  if (sub == "wall") return cmdFix(h, true);
  if (sub == "gravity") {
    double g[3] = {0, 0, 0};
    if (!readDoubles(nd, g)) return fail("gravity needs ndm values");
    for (int d = 0; d < 3; d++) E.par.gravity[d] = g[d];
    return 0;
  }
  if (sub == "timeSeries") {
    int ts;
    if (!readInt(ts)) return fail("timeSeries <tsTag>");
    TimeSeries *s = OPS_getTimeSeries(ts);       // a copy
    if (s == 0) return fail("timeSeries -- no time series " + std::to_string(ts));
    h->setTimeSeries(s);
    return 0;
  }
  if (sub == "damping") {
    if (OPS_GetNumRemainingInputArgs() < 2 || strcmp(OPS_GetString(), "local") != 0)
      return fail("damping local <alpha>");
    double a;
    if (!readDoubles(1, &a) || a < 0.0 || a >= 1.0) return fail("damping local -- alpha in [0, 1)");
    E.par.damping = a;
    return 0;
  }
  if (sub == "flip") {
    double a;
    if (!readDoubles(1, &a) || a < 0.0 || a > 1.0) return fail("flip -- the FLIP fraction is in [0, 1]");
    E.par.flip = a;
    return 0;
  }
  if (sub == "massCutoff") {
    double a;
    if (!readDoubles(1, &a) || a < 0.0) return fail("massCutoff -- c_m >= 0");
    h->setMassCutoffUser(a);      // wins in both modes (mode-A default 1e-6, mode F 1e-10)
    return 0;
  }
  if (sub == "mode") {
    const std::string m = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "";
    LadrunoMPMPattern::Mode md;
    if (m == "explicit") md = LadrunoMPMPattern::MODE_EXPLICIT;
    else if (m == "implicit") md = LadrunoMPMPattern::MODE_IMPLICIT;
    else if (m == "quasiStatic") md = LadrunoMPMPattern::MODE_QUASISTATIC;
    else return fail("mode explicit|implicit|quasiStatic");
    std::string err;
    if (h->setMode(md, err) < 0) return fail("mode -- " + err);
    return 0;
  }
  if (sub == "gridvel") return cmdGridVel(h);
  if (sub == "pforce") {
    int body;
    double lo[3], hi[3], f[3];
    if (!readInt(body) || OPS_GetNumRemainingInputArgs() < 1 || strcmp(OPS_GetString(), "box") != 0 ||
        !readDoubles(nd, lo) || !readDoubles(nd, hi) || !readDoubles(nd, f))
      return fail("pforce <body> box <lo..> <hi..> <f..>");
    std::string err;
    const int n = h->addPointForce(body, lo, hi, f, err);
    if (n < 0) return fail("pforce -- " + err);
    return outInt(n);
  }
  if (sub == "emptyNodes") {
    const std::string m = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "";
    if (m == "exclude") { h->setEmptyNodes(LadrunoMPMPattern::EMPTY_EXCLUDE); return 0; }
    if (m == "renormalize") { h->setEmptyNodes(LadrunoMPMPattern::EMPTY_RENORMALIZE); return 0; }
    double g;
    if (m == "ghost" && readDoubles(1, &g) && g >= 0.0) { h->setEmptyNodes(g); return 0; }
    return fail("emptyNodes exclude | renormalize | ghost <gamma >= 0>");
  }
  if (sub == "cut") {
    double fct;
    int k;
    if (!readDoubles(1, &fct) || !readInt(k) || !(fct > 1.0) || k < 0)
      return fail("cut <factor > 1> <maxCuts >= 0>");
    h->setCut(fct, k);
    return 0;
  }
  if (sub == "tagBase") {
    int nb, eb;
    if (!readInt(nb) || !readInt(eb)) return fail("tagBase <nodeBase> <eleBase>");
    std::string err;
    if (h->setTagBases(nb, eb, err) < 0) return fail("tagBase -- " + err);
    return 0;
  }
  if (sub == "rayleigh") {
    double r[2];
    if (!readDoubles(2, r)) return fail("rayleigh <alphaM> <betaK> (betaK0 / betaKc are refused at particles)");
    h->setRayleighMPM(r[0], r[1]);
    return 0;
  }
  if (sub == "safety") {
    double a;
    if (!readDoubles(1, &a) || !(a > 0.0) || a > 1.0) return fail("safety -- s in (0, 1]");
    E.par.safety = a;
    return 0;
  }
  if (sub == "criticalTimeStep" || sub == "particleData" || sub == "condForce" ||
      (sub == "info" && OPS_GetNumRemainingInputArgs() > 0)) {
    std::vector<std::string> a(1, sub);
    if (sub != "criticalTimeStep" && OPS_GetNumRemainingInputArgs() > 0) a.push_back(OPS_GetString());
    std::vector<double> out;
    const int q = ladrunoMPMDoubleQuery(dom, a, out);
    if (q < 0) return -1;
    // a scalar query answers a scalar (criticalTimeStep, info time); the rest are lists
    const bool scalar = sub == "criticalTimeStep" || (sub == "info" && a.size() > 1 && a[1] == "time");
    if (q > 0) return scalar ? outDouble(out[0]) : outDoubles(out);
    if (sub == "info") OPS_ResetCurrentInputArg(-1);   // not a double query: the int forms below
    if (sub == "particleData") return fail("particleData <field>");
    if (sub == "condForce") return fail("condForce <name>");
  }
  if (sub == "analyze") {
    int n;
    double dt;
    if (!readInt(n) || !readDoubles(1, &dt)) return fail("analyze <nSteps> <dt>");
    if (h->mode() != LadrunoMPMPattern::MODE_EXPLICIT) {
      LadrunoMPMAnalysisRefs none;
      if (h->analyzeImplicit(n, dt, g_refs != 0 ? *g_refs : none) < 0) return -1;
      return outInt(0);
    }
    if (h->analyze(n, dt) < 0) return -1;
    return outInt(0);
  }
  if (sub == "reset") return h->reset();
  if (sub == "info") {
    const std::string what = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetString() : "numParticles";
    if (what == "numParticles") return outInt(E.P.size());
    if (what == "epoch") return outInt(h->particleEpoch());
    bool ok = false;
    const long c = h->infoCounter(what, ok);
    if (ok) return outInt((int)c);
    return fail("info -- numParticles, epoch, time, summary, g2p, acceptedSteps, retries, dtHistory,"
                " domainChanges, connectivityChanges, gridNodes, gridNodeTags or mode");
  }
  if (sub == "debug") return cmdDebug(h);
  return fail("-- unknown sub-command '" + sub + "'");
}
