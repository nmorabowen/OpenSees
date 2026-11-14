// mpm_engine_check — standalone g++ check of the explicit MPM engine core
// (WP-232c; SRC/domain/pattern/ladrunoMPM/LadrunoMPMCore.h) against the frozen
// case file tests/mpm_cases.txt (WP-232b). No OpenSees dependency.
//
//   g++ -O2 -std=c++17 -Wall -Wextra -I SRC/element/solidTransformation
//       -I SRC/domain/pattern/ladrunoMPM tests/mpm_engine_check.cpp -o mpmc
//   ./mpmc tests/mpm_cases.txt
//
// The material is the interface's §9 oracle law, hypoelastic isotropic
// (sig_hat = sig0_hat + D eps_hat), so this checks the engine's numerics — basis,
// P2G / G2P, MUSL, PIC / FLIP, mass cutoff, damping, critical Δt, grid
// conditions, the frictional wall, the kernel entry — independently of the
// OpenSees wiring, which tests/test_ladruno_mpm_engine.py checks through the
// binaries (ElasticIsotropic must reproduce this law). Values marked
// `pending-oracle` are not compared; they are printed as INFO lines so the
// orchestrator can set them against the oracle's when 232b merges.
//
// Exit 0 + "RESULT: PASS" on success.

#include <LadrunoMPMCore.h>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace ladruno_mpm;

static int g_fail = 0, g_checks = 0, g_pending = 0;

// ---- the oracle material (interface §9) ------------------------------------
struct Hypoelastic : public MaterialPoints {
  int nd = 2;
  std::vector<double> lam, mu;    // per particle
  int setTrialStrain(int p, const double *e, double *s) override
  {
    const double l = lam[p], m = mu[p];
    if (nd == 3) {
      const double tr = e[0] + e[1] + e[2];
      for (int k = 0; k < 3; k++) s[k] = l * tr + 2.0 * m * e[k];
      for (int k = 3; k < 6; k++) s[k] = m * e[k];
    } else {
      const double tr = e[0] + e[1];
      s[0] = l * tr + 2.0 * m * e[0];
      s[1] = l * tr + 2.0 * m * e[1];
      s[2] = m * e[2];
    }
    return 0;
  }
  int commitState(int) override { return 0; }
  int revertToLastCommit(int) override { return 0; }
  int pendingRefusals() override { return 0; }
  double normalModulus(int p) override { return lam[p] + 2.0 * mu[p]; }
};

// ---- case model --------------------------------------------------------------
struct Expect {
  std::string obs;
  std::vector<std::string> args;
  std::string op;
  std::vector<double> vals;
  bool pending = false;
  bool fields = false;
  bool hasTol = false;
  double atol = 0.0, rtol = 0.0;
  std::string line;
};

struct Case {
  std::string name, kind;
  int ndm = 2;
  std::string basis = "gimp";
  double origin[3] = {0, 0, 0}, h = 1.0;
  int cells[3] = {1, 1, 1};
  std::map<int, std::pair<double, double>> mats;
  std::vector<std::vector<std::string>> bodies, shifts, velocities, stresses, particles, fixes, walls;
  double point[3] = {0, 0, 0}, plp[3] = {0, 0, 0};
  double gravity[3] = {0, 0, 0};
  double seriesTr = -1.0;
  double damping = 0.0, flip = 1.0, masscut = 1e-10, safety = 0.5, dt = 0.0;
  int steps = 0;
  double atol = 0.0, rtol = 0.0;
  std::vector<Expect> expects;
  std::vector<std::vector<std::string>> compares;
  std::vector<std::vector<std::string>> fields;   // field <id> x .. v .. F .. V .. sig .. sighat .. epshat ..
};

static std::vector<std::string> split(const std::string &s)
{
  std::istringstream ss(s);
  std::vector<std::string> t;
  std::string w;
  while (ss >> w) t.push_back(w);
  return t;
}

static double D(const std::string &s) { return std::stod(s); }

static std::vector<Case> parse(const char *path)
{
  std::vector<Case> out;
  std::ifstream in(path);
  if (!in) { std::printf("FAIL: no cases file %s\n", path); g_fail++; return out; }
  std::string line;
  Case c;
  bool open = false;
  while (std::getline(in, line)) {
    auto t = split(line);
    if (t.empty() || t[0][0] == '#') continue;
    const std::string &k = t[0];
    if (k == "case") { c = Case(); c.name = t[1]; open = true; continue; }
    if (!open) continue;
    if (k == "end") { out.push_back(c); open = false; continue; }
    if (k == "kind") c.kind = t[1];
    else if (k == "desc" || k == "note") {}
    else if (k == "ndm") c.ndm = std::stoi(t[1]);
    else if (k == "basis") c.basis = t[1];
    else if (k == "origin") for (int d = 0; d < c.ndm; d++) c.origin[d] = D(t[1 + d]);
    else if (k == "h") c.h = D(t[1]);
    else if (k == "cells") for (int d = 0; d < c.ndm; d++) c.cells[d] = std::stoi(t[1 + d]);
    else if (k == "material") c.mats[std::stoi(t[1])] = {D(t[3]), D(t[4])};
    else if (k == "body") c.bodies.push_back(t);
    else if (k == "shift") c.shifts.push_back(t);
    else if (k == "velocity") c.velocities.push_back(t);
    else if (k == "stress") c.stresses.push_back(t);
    else if (k == "particle") c.particles.push_back(t);
    else if (k == "point") {
      for (int d = 0; d < c.ndm; d++) c.point[d] = D(t[1 + d]);
      for (int d = 0; d < c.ndm; d++) c.plp[d] = D(t[2 + c.ndm + d]);
    }
    else if (k == "fix") c.fixes.push_back(t);
    else if (k == "wall") c.walls.push_back(t);
    else if (k == "gravity") for (int d = 0; d < c.ndm; d++) c.gravity[d] = D(t[1 + d]);
    else if (k == "series") c.seriesTr = (t[1] == "linear") ? D(t[2]) : -1.0;
    else if (k == "damping") c.damping = D(t[1]);
    else if (k == "flip") c.flip = D(t[1]);
    else if (k == "masscut") c.masscut = D(t[1]);
    else if (k == "safety") c.safety = D(t[1]);
    else if (k == "dt") c.dt = D(t[1]);
    else if (k == "steps") c.steps = std::stoi(t[1]);
    else if (k == "tol") { c.atol = D(t[1]); c.rtol = D(t[2]); }
    else if (k == "compare") c.compares.push_back(t);
    else if (k == "field") c.fields.push_back(t);
    else if (k == "expect") {
      Expect e;
      e.line = line;
      size_t i = 1;
      e.obs = t[i++];
      while (i < t.size() && t[i] != "eq" && t[i] != "le" && t[i] != "ge") e.args.push_back(t[i++]);
      if (i < t.size()) e.op = t[i++];
      while (i < t.size() && t[i] != "tol") {
        if (t[i] == "pending-oracle") e.pending = true;
        else if (t[i] == "fields") e.fields = true;
        else e.vals.push_back(D(t[i]));
        i++;
      }
      if (i < t.size() && t[i] == "tol") { e.hasTol = true; e.atol = D(t[i + 1]); e.rtol = D(t[i + 2]); }
      c.expects.push_back(e);
    }
    else { std::printf("FAIL: case %s unknown key '%s'\n", c.name.c_str(), k.c_str()); g_fail++; }
  }
  return out;
}

// ---- building an engine from a case -----------------------------------------
static bool build(const Case &c, Engine &E, Hypoelastic &M, std::vector<double> &X0)
{
  std::string err;
  if (!E.setGrid(c.ndm, c.origin, c.h, c.cells, err)) { std::printf("FAIL: %s grid: %s\n", c.name.c_str(), err.c_str()); g_fail++; return false; }
  E.par.basis = (c.basis == "linear") ? BASIS_LINEAR : BASIS_GIMP;
  E.par.flip = c.flip; E.par.damping = c.damping; E.par.masscut = c.masscut; E.par.safety = c.safety;
  for (int d = 0; d < 3; d++) E.par.gravity[d] = c.gravity[d];
  const int nd = c.ndm, nv = (nd == 3) ? 6 : 3;
  int64_t next = 0;
  for (const auto &b : c.bodies) {
    // body <id> box lo hi ppc n mat tag rho r | body <id> disc cx cy R ppc ...
    const int bid = std::stoi(b[1]);
    int shape = (b[2] == "box") ? 0 : ((b[2] == "disc") ? 1 : 2);
    const int ng = (shape == 0) ? 2 * nd : nd + 1;
    std::vector<double> geom;
    for (int k = 0; k < ng; k++) geom.push_back(D(b[3 + k]));
    size_t i = 3 + ng;
    int ppc = std::stoi(b[i + 1]), mat = std::stoi(b[i + 3]);
    double rho = D(b[i + 5]);
    if (E.seedLattice(shape, geom.data(), ppc, bid, mat, rho, next, err) < 0) {
      std::printf("FAIL: %s seed: %s\n", c.name.c_str(), err.c_str()); g_fail++; return false;
    }
  }
  for (const auto &sh : c.shifts) {   // shift <body> <ndm floats> (A5): seeded particles move
    double dx[3] = {0, 0, 0};
    for (int d = 0; d < nd; d++) dx[d] = D(sh[2 + d]);
    E.shiftBody(std::stoi(sh[1]), dx);
  }
  for (const auto &pt : c.particles) {
    // particle id body mat rho V0 x .. lp .. v .. s ..
    int64_t id = std::stoll(pt[1]);
    int body = std::stoi(pt[2]), mat = std::stoi(pt[3]);
    double rho = D(pt[4]), V0 = D(pt[5]);
    double x[3], lp[3], v[3], s[6];
    size_t i = 7;
    for (int d = 0; d < nd; d++) x[d] = D(pt[i++]);
    i++;
    for (int d = 0; d < nd; d++) lp[d] = D(pt[i++]);
    i++;
    for (int d = 0; d < nd; d++) v[d] = D(pt[i++]);
    i++;
    for (int k = 0; k < nv; k++) s[k] = D(pt[i++]);
    if (E.addParticle(id, body, mat, rho, V0, x, lp, v, s, err) < 0) {
      std::printf("FAIL: %s particle: %s\n", c.name.c_str(), err.c_str()); g_fail++; return false;
    }
  }
  const int N = E.P.size();
  for (const auto &vl : c.velocities) {
    const int bid = std::stoi(vl[1]);
    for (int p = 0; p < N; p++) {
      if (E.P.body[p] != bid) continue;
      if (vl[2] == "uniform") {
        for (int d = 0; d < nd; d++) E.P.v[nd * p + d] = D(vl[3 + d]);
      } else {   // sine comp axis v0 beta xref
        const int comp = vl[3][0] - 'x', axis = vl[4][0] - 'x';
        const double v0 = D(vl[5]), beta = D(vl[6]), xr = D(vl[7]);
        for (int d = 0; d < nd; d++) E.P.v[nd * p + d] = 0.0;
        E.P.v[nd * p + comp] = v0 * std::sin(beta * (E.P.x[nd * p + axis] - xr));
      }
    }
  }
  for (const auto &st : c.stresses) {   // stress body geostatic ytop rho g K0 (2D)
    const int bid = std::stoi(st[1]);
    const double yt = D(st[3]), rho = D(st[4]), g = D(st[5]), K0 = D(st[6]);
    for (int p = 0; p < N; p++) {
      if (E.P.body[p] != bid) continue;
      const double syy = -rho * g * (yt - E.P.x[nd * p + 1]);
      const double s[3] = {K0 * syy, syy, 0.0};
      for (int k = 0; k < 3; k++)
        E.P.sig0[3 * p + k] = E.P.sigHat[3 * p + k] = E.P.sig[3 * p + k] = s[k];
    }
  }
  for (const auto &f : c.fixes) {
    const int fi = faceIndex(f[1]);
    for (char ch : f[2]) E.fc.fixMask[fi] |= 1u << (ch - 'x');
  }
  for (const auto &w : c.walls) E.fc.wallMu[faceIndex(w[1])] = D(w[2]);
  M.nd = nd;
  M.lam.resize(N); M.mu.resize(N);
  for (int p = 0; p < N; p++) {
    auto em = c.mats.at(E.P.mat[p]);
    const double Ey = em.first, nu = em.second;
    M.lam[p] = Ey * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    M.mu[p] = Ey / (2.0 * (1.0 + nu));
  }
  X0 = E.P.x;
  return true;
}

// ---- observables ---------------------------------------------------------------
struct Obs {
  std::vector<double> comvx;      // com_v_x per recorded state
  std::vector<double> times;
  double maxSpeed = 0.0;
  double momDev = 0.0;            // max |momentum| component deviation from the expectation
  double E0 = 0.0, Eend = 0.0;
  bool refused = false, accepted = true;
};

static double energy(const Engine &E, const Hypoelastic &M)
{
  const int N = E.P.size(), nd = E.ndm(), nv = E.nv();
  double ke = 0.0, se = 0.0;
  for (int p = 0; p < N; p++) {
    double v2 = 0.0;
    for (int d = 0; d < nd; d++) v2 += E.P.v[nd * p + d] * E.P.v[nd * p + d];
    ke += 0.5 * E.P.m[p] * v2;
    double s[6];
    Hypoelastic &MM = const_cast<Hypoelastic &>(M);
    MM.setTrialStrain(p, &E.P.eps[nv * p], s);
    double w = 0.0;
    for (int k = 0; k < nv; k++) w += E.P.eps[nv * p + k] * s[k];
    se += 0.5 * E.P.V[p] * w;
  }
  return ke + se;
}

static void sample(const Engine &E, Obs &o, double t)
{
  const int N = E.P.size(), nd = E.ndm();
  double mt = 0.0, mvx = 0.0;
  for (int p = 0; p < N; p++) {
    mt += E.P.m[p];
    mvx += E.P.m[p] * E.P.v[nd * p];
    double v2 = 0.0;
    for (int d = 0; d < nd; d++) v2 += E.P.v[nd * p + d] * E.P.v[nd * p + d];
    o.maxSpeed = std::max(o.maxSpeed, std::sqrt(v2));
  }
  o.comvx.push_back(mvx / mt);
  o.times.push_back(t);
}

static double period(const Obs &o)
{
  std::vector<double> cr;
  for (size_t k = 0; k + 1 < o.comvx.size(); k++)
    if (o.comvx[k] < 0.0 && 0.0 <= o.comvx[k + 1]) {
      const double a = o.comvx[k], b = o.comvx[k + 1];
      cr.push_back(o.times[k] + (o.times[k + 1] - o.times[k]) * (-a) / (b - a));
    }
  if (cr.size() < 2) return NAN;
  return (cr.back() - cr.front()) / (cr.size() - 1);
}

static std::map<std::string, std::map<std::string, double>> g_scalar;   // case -> obs -> value

static void verdict(const Case &c, const Expect &e, const std::vector<double> &got)
{
  g_checks++;
  if (e.pending) {
    g_pending++;
    std::printf("INFO: %s %s =", c.name.c_str(), e.obs.c_str());
    for (double g : got) std::printf(" %.17e", g);
    std::printf("  (pending-oracle)\n");
    return;
  }
  const double at = e.hasTol ? e.atol : c.atol, rt = e.hasTol ? e.rtol : c.rtol;
  bool ok = !got.empty();
  if (e.op == "eq") {
    ok = ok && got.size() == e.vals.size();
    for (size_t i = 0; ok && i < got.size(); i++)
      ok = std::fabs(got[i] - e.vals[i]) <= at + rt * std::fabs(e.vals[i]);
  } else if (e.op == "le") {
    ok = ok && got[0] <= e.vals[0];
  } else if (e.op == "ge") {
    ok = ok && got[0] >= e.vals[0];
  }
  std::printf("%s: %s %s", ok ? "ok  " : "FAIL", c.name.c_str(), e.obs.c_str());
  for (double g : got) std::printf(" %.10e", g);
  std::printf("  [%s", e.op.c_str());
  for (double v : e.vals) std::printf(" %.10e", v);
  std::printf("]\n");
  if (!ok) g_fail++;
}

static void runCase(const Case &c)
{
  Engine E;
  Hypoelastic M;
  std::vector<double> X0;
  const int nd = c.ndm;
  if (c.kind == "shape1d") {
    for (const auto &e : c.expects) {
      double S, dS;
      gimpFactor(D(e.args[2]), D(e.args[0]), D(e.args[1]), S, dS);
      verdict(c, e, {S, dS});
    }
    return;
  }
  if (c.kind == "shape") {
    E.grid.ndm = nd; E.grid.h = c.h;
    for (int d = 0; d < nd; d++) { E.grid.x0[d] = c.origin[d]; E.grid.n[d] = c.cells[d]; }
    const Basis b = (c.basis == "linear") ? BASIS_LINEAR : BASIS_GIMP;
    std::vector<long> node; std::vector<double> S, dS;
    particleSupport(E.grid, b, c.point, c.plp, node, S, dS);
    std::map<long, std::vector<double>> got;
    for (size_t k = 0; k < node.size(); k++) {
      std::vector<double> v = {S[k]};
      for (int d = 0; d < nd; d++) v.push_back(dS[nd * k + d]);
      got[E.grid.publicId(node[k])] = v;
    }
    std::map<long, bool> listed;
    for (const auto &e : c.expects) {
      const long id = std::stol(e.args[0]);
      listed[id] = true;
      verdict(c, e, got.count(id) ? got[id] : std::vector<double>(nd + 1, 0.0));
    }
    for (auto &kv : got) {   // any other node: zero weight and gradient
      if (listed.count(kv.first)) continue;
      double s = 0.0;
      for (double v : kv.second) s += std::fabs(v);
      g_checks++;
      if (s > 1e-15) { std::printf("FAIL: %s unlisted node %ld non-zero\n", c.name.c_str(), kv.first); g_fail++; }
    }
    return;
  }
  if (!build(c, E, M, X0)) return;
  if (c.kind == "critdt") {
    for (const auto &e : c.expects) verdict(c, e, {E.criticalTimeStep(M)});
    return;
  }
  if (c.kind == "masscut") {
    std::vector<long> act;
    E.activeNodes(act);
    for (const auto &e : c.expects) verdict(c, e, std::vector<double>(act.begin(), act.end()));
    return;
  }
  // ---- run ----
  Obs o;
  double t = 0.0;
  sample(E, o, t);
  o.E0 = energy(E, M);
  auto momentum = [&](std::vector<double> &mv) {
    mv.assign(nd, 0.0);
    for (int p = 0; p < E.P.size(); p++)
      for (int d = 0; d < nd; d++) mv[d] += E.P.m[p] * E.P.v[nd * p + d];
  };
  std::vector<double> mv, momMax(nd, 0.0);
  momentum(mv);
  for (int d = 0; d < nd; d++) momMax[d] = std::max(momMax[d], std::fabs(mv[d]));
  for (int s = 0; s < c.steps; s++) {
    const double dtc = E.criticalTimeStep(M);
    const int bad = E.firstInadmissible();
    if (bad >= 0 || c.dt > E.par.safety * dtc) {
      if (s == 0) o.refused = true;
      o.accepted = false;
      std::printf("note: %s refused at step %d (%s)\n", c.name.c_str(), s,
                  bad >= 0 ? "inadmissible particle" : "dt above s*dt_crit");
      break;
    }
    const double t1 = t + c.dt;
    const double lam = (c.seriesTr > 0.0) ? std::min(t1 / c.seriesTr, 1.0) : 1.0;
    std::string err;
    if (E.advance(c.dt, lam, M, err) != 0) {
      o.accepted = false;
      std::printf("note: %s failed at step %d: %s\n", c.name.c_str(), s, err.c_str());
      break;
    }
    t = t1;
    sample(E, o, t);
    momentum(mv);
    for (int d = 0; d < nd; d++) momMax[d] = std::max(momMax[d], std::fabs(mv[d]));
  }
  o.Eend = energy(E, M);
  const int N = E.P.size(), nv = E.nv(), nn = nd * nd;
  for (const auto &e : c.expects) {
    const std::string &ob = e.obs;
    std::vector<double> got;
    if (ob == "refused") { g_checks++; if (!o.refused) { std::printf("FAIL: %s expected refused\n", c.name.c_str()); g_fail++; } else std::printf("ok  : %s refused\n", c.name.c_str()); continue; }
    if (ob == "accepted") { g_checks++; if (!o.accepted) { std::printf("FAIL: %s expected accepted\n", c.name.c_str()); g_fail++; } else std::printf("ok  : %s accepted\n", c.name.c_str()); continue; }
    if (ob == "particle_fields") {
      if (e.pending || c.fields.empty()) { g_pending++; g_checks++; continue; }
      const double at = e.atol, rt = e.rtol;
      int bad = 0;
      double worst = 0.0;
      for (const auto &fl : c.fields) {
        const int64_t pid = std::stoll(fl[1]);
        int p = -1;
        for (int q = 0; q < N; q++) if (E.P.id[q] == pid) p = q;
        if (p < 0) { bad++; continue; }
        std::vector<double> gotv;
        for (int d = 0; d < nd; d++) gotv.push_back(E.P.x[nd * p + d]);
        for (int d = 0; d < nd; d++) gotv.push_back(E.P.v[nd * p + d]);
        for (int k = 0; k < nn; k++) gotv.push_back(E.P.F[nn * p + k]);
        gotv.push_back(E.P.V[p]);
        for (int k = 0; k < nv; k++) gotv.push_back(E.P.sig[nv * p + k]);
        for (int k = 0; k < nv; k++) gotv.push_back(E.P.sigHat[nv * p + k]);
        for (int k = 0; k < nv; k++) gotv.push_back(E.P.eps[nv * p + k]);
        std::vector<double> expv;
        for (size_t i = 2; i < fl.size(); i++)
          if (std::isdigit((unsigned char)fl[i][0]) || fl[i][0] == '-' || fl[i][0] == '+' || fl[i][0] == '.')
            expv.push_back(D(fl[i]));
        if (expv.size() != gotv.size()) { bad++; continue; }
        for (size_t i = 0; i < gotv.size(); i++) {
          const double dd = std::fabs(gotv[i] - expv[i]);
          worst = std::max(worst, dd / (at + rt * std::fabs(expv[i])));
          if (dd > at + rt * std::fabs(expv[i])) bad++;
        }
      }
      g_checks++;
      std::printf("%s: %s particle_fields (%d particles, worst diff/tol %.3e)\n", bad ? "FAIL" : "ok  ",
                  c.name.c_str(), (int)c.fields.size(), worst);
      if (bad) g_fail++;
      continue;
    }
    if (ob == "translation_error") {
      double worst = 0.0;
      for (int p = 0; p < N; p++) {
        double r2 = 0.0;
        for (int d = 0; d < nd; d++) {
          const double ex = X0[nd * p + d] + D(e.args[d]) * t;
          r2 += (E.P.x[nd * p + d] - ex) * (E.P.x[nd * p + d] - ex);
        }
        worst = std::max(worst, std::sqrt(r2));
      }
      got = {worst};
    } else if (ob == "F_error") {
      double w = 0.0;
      for (int p = 0; p < N; p++)
        for (int a = 0; a < nd; a++)
          for (int b = 0; b < nd; b++)
            w = std::max(w, std::fabs(E.P.F[nn * p + nd * a + b] - (a == b ? 1.0 : 0.0)));
      got = {w};
    } else if (ob == "max_abs_stress") {
      double w = 0.0;
      for (double s : E.P.sig) w = std::max(w, std::fabs(s));
      got = {w};
    } else if (ob == "period_comx") {
      got = {period(o)};
    } else if (ob == "energy_ratio_end") {
      got = {o.Eend / o.E0};
    } else if (ob == "selfweight_error" || ob == "cell_crossing_meter") {
      const double rho0 = D(e.args[0]), g = D(e.args[1]), l0 = D(e.args[2]), yb = D(e.args[3]);
      double w = 0.0, ss = 0.0;
      for (int p = 0; p < N; p++) {
        const double ex = -rho0 * g * (l0 - (X0[nd * p + 1] - yb));
        const double r = E.P.sig[nv * p + 1] - ex;
        w = std::max(w, std::fabs(r));
        ss += r * r;
      }
      got = {(ob == "selfweight_error") ? w / (rho0 * g * l0) : std::sqrt(ss / N) / (rho0 * g * l0)};
    } else if (ob == "max_speed") {
      got = {o.maxSpeed};
    } else if (ob == "momentum") {
      got = momMax;
    } else if (ob == "com_disp" || ob == "com_v") {
      const int bid = std::stoi(e.args[0]);
      double mt = 0.0, a[3] = {0, 0, 0};
      for (int p = 0; p < N; p++) {
        if (E.P.body[p] != bid) continue;
        mt += E.P.m[p];
        for (int d = 0; d < nd; d++)
          a[d] += E.P.m[p] * ((ob == "com_v") ? E.P.v[nd * p + d] : (E.P.x[nd * p + d] - X0[nd * p + d]));
      }
      for (int d = 0; d < nd; d++) got.push_back(a[d] / mt);
    } else if (ob == "com_disp_x" || ob == "com_v_x") {
      const int bid = std::stoi(e.args[0]);
      double mt = 0.0, a = 0.0, a0 = 0.0;
      for (int p = 0; p < N; p++) {
        if (E.P.body[p] != bid) continue;
        mt += E.P.m[p];
        if (ob == "com_v_x") a += E.P.m[p] * E.P.v[nd * p];
        else { a += E.P.m[p] * E.P.x[nd * p]; a0 += E.P.m[p] * X0[nd * p]; }
      }
      got = {(a - a0) / mt};
    } else {
      std::printf("FAIL: %s unknown observable %s\n", c.name.c_str(), ob.c_str());
      g_fail++;
      continue;
    }
    if (!got.empty() && e.args.size() <= 4) g_scalar[c.name][ob] = got[0];
    verdict(c, e, got);
  }
}

int main(int argc, char **argv)
{
  {  // v1 has one cell size for every axis: unequal per-axis sizes are refused
    Engine E;
    std::string err;
    const double x0[2] = {0.0, 0.0}, hs[2] = {0.1, 0.2}, heq[2] = {0.1, 0.1};
    const int n[2] = {4, 4};
    g_checks++;
    if (E.setGrid(2, x0, hs, n, err) || err.empty()) { std::printf("FAIL: unequal cell sizes accepted\n"); g_fail++; }
    g_checks++;
    if (!E.setGrid(2, x0, heq, n, err)) { std::printf("FAIL: equal per-axis cell sizes refused\n"); g_fail++; }
  }
  {  // A5 ghost layers: a block sliding along a roller face, its particle domains
     // overhanging the face (delta/lp = 0.5), keeps its velocity exactly — with
     // truncated support the overhang made sum(grad S) != 0 and the block dragged.
    Engine E;
    Hypoelastic M;
    std::string err;
    const double x0[2] = {0.0, 0.0};
    const int n[2] = {20, 10};
    E.setGrid(2, x0, 0.1, n, err);
    E.fc.fixMask[2] |= 1u << 1;                       // ymin roller
    const double lp[2] = {0.025, 0.025};
    int64_t id = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 8; i++) {
        const double x[2] = {0.5 + 0.05 * i, 0.0125 + 0.05 * j};   // bottom row: y - lp = -0.0125
        const double v[2] = {0.3, 0.0};
        E.addParticle(id++, 1, 1, 1000.0, 0.0025, x, lp, v, 0, err);
      }
    M.nd = 2;
    M.lam.assign(E.P.size(), 1.0e6 * 0.3 / (1.3 * 0.4));
    M.mu.assign(E.P.size(), 1.0e6 / 2.6);
    double worst = 0.0;
    for (int s = 0; s < 200; s++) {
      if (E.advance(1.0e-4, 1.0, M, err) != 0) { worst = 1.0; break; }
      for (int p = 0; p < E.P.size(); p++)
        worst = std::max(worst, std::fabs(E.P.v[2 * p] - 0.3) + std::fabs(E.P.v[2 * p + 1]));
    }
    g_checks++;
    std::printf("%s: A5 roller slide with overhang, max |v - v0| = %.3e\n",
                worst < 1e-12 ? "ok  " : "FAIL", worst);
    if (!(worst < 1e-12)) g_fail++;
  }
  {  // Hand check 1 (A5 F2): friction acts ONCE per step. One unstressed particle
     // moving into a ymin wall. With one particle v_i^n = v_p at every support
     // node; the full operator gives wall nodes v1 = (vx - mu|vy|, 0); FLIP from a
     // single particle gives v_p1 = sum S v1; the second projection gives v* = v_p1
     // at every node, minus ONLY the wall-normal part at wall nodes. So F1 = I + dt L
     // with L = sum v*_i (x) grad S_i, known by hand. Friction applied in both
     // projections changes L_xy.
    Engine E;
    Hypoelastic M;
    std::string err;
    const double x0[2] = {0.0, 0.0};
    const int n[2] = {4, 4};
    E.setGrid(2, x0, 1.0, n, err);
    E.par.masscut = 0.0;
    const double mu = 0.3, dt = 1.0e-3, vp[2] = {0.5, -0.1};
    E.fc.wallMu[2] = mu;
    const double xp[2] = {2.3, 0.2}, lp[2] = {0.25, 0.25};
    E.addParticle(0, 1, 1, 4.0, 0.25, xp, lp, vp, 0, err);
    M.nd = 2; M.lam.assign(1, 100.0); M.mu.assign(1, 100.0);
    std::vector<long> node; std::vector<double> S, dS;
    particleSupport(E.grid, BASIS_GIMP, xp, lp, node, S, dS);
    std::vector<double> v1(2 * node.size());
    double vp1[2] = {0.0, 0.0}, dx[2] = {0.0, 0.0};
    for (size_t k = 0; k < node.size(); k++) {
      int ijk[3];
      E.grid.extIJK(node[k], ijk);
      const bool wall = ijk[1] <= 0;
      v1[2 * k] = wall ? vp[0] - std::min(mu * std::fabs(vp[1]), std::fabs(vp[0])) : vp[0];
      v1[2 * k + 1] = wall ? 0.0 : vp[1];
      for (int a = 0; a < 2; a++) { vp1[a] += S[k] * v1[2 * k + a]; dx[a] += dt * S[k] * v1[2 * k + a]; }
    }
    double L[4] = {0, 0, 0, 0}, Ltwice_xy = 0.0;   // Ltwice: the friction-in-both mutant
    for (size_t k = 0; k < node.size(); k++) {
      int ijk[3];
      E.grid.extIJK(node[k], ijk);
      double vs[2] = {vp1[0], vp1[1]};
      double vsx_twice = vs[0];
      if (ijk[1] <= 0 && vs[1] < 0.0) {
        vsx_twice = vs[0] - std::min(mu * std::fabs(vs[1]), std::fabs(vs[0]));
        vs[1] = 0.0;                                          // no-penetration part only
      }
      for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++) L[2 * a + b] += vs[a] * dS[2 * k + b];
      Ltwice_xy += vsx_twice * dS[2 * k + 1];
    }
    const double Fexp[4] = {1.0 + dt * L[0], dt * L[1], dt * L[2], 1.0 + dt * L[3]};
    E.advance(dt, 1.0, M, err);
    double worst = 0.0;
    for (int k = 0; k < 4; k++) worst = std::max(worst, std::fabs(E.P.F[k] - Fexp[k]));
    for (int a = 0; a < 2; a++) {
      worst = std::max(worst, std::fabs(E.P.x[a] - (xp[a] + dx[a])));
      worst = std::max(worst, std::fabs(E.P.v[a] - vp1[a]));
    }
    g_checks++;
    // non-vacuous: friction in both projections would make L_xy non-zero
    const bool ok = worst < 1e-14 && std::fabs(Ltwice_xy - L[1]) > 1e-3;
    std::printf("%s: hand check, friction once per step (max diff %.3e; L_xy %.3e, twice-mutant %.4e)\n",
                ok ? "ok  " : "FAIL", worst, L[1], Ltwice_xy);
    if (!ok) g_fail++;
  }
  {  // Hand check 2 (A5 F2): the mass cutoff acts inside the step. One stressed
     // particle at rest: v_p1 = -dt (V/m) sigma . sum_{active} grad S_i, with node i
     // active iff S_i > c_m (m_i = S_i m, m_ref = m). Without the cutoff the sum
     // runs over every node and is exactly 0 (sum grad S = 0).
    Engine E;
    Hypoelastic M;
    std::string err;
    const double x0[2] = {0.0, 0.0};
    const int n[2] = {4, 4};
    E.setGrid(2, x0, 1.0, n, err);
    E.par.masscut = 1.0e-4;
    const double dt = 1.0e-3, V = 0.25, m = 250.0, sg[3] = {100.0, -40.0, 25.0};
    const double xp[2] = {3.245, 1.5}, lp[2] = {0.25, 0.25}, v0[2] = {0.0, 0.0};
    E.addParticle(0, 1, 1, m / V, V, xp, lp, v0, sg, err);
    M.nd = 2; M.lam.assign(1, 4.0e3); M.mu.assign(1, 4.0e3);
    std::vector<long> node; std::vector<double> S, dS;
    particleSupport(E.grid, BASIS_GIMP, xp, lp, node, S, dS);
    double gsum[2] = {0, 0};
    int dropped = 0;
    for (size_t k = 0; k < node.size(); k++) {
      if (!(S[k] > E.par.masscut)) { dropped++; continue; }
      gsum[0] += dS[2 * k]; gsum[1] += dS[2 * k + 1];
    }
    const double vexp[2] = {-dt * (V / m) * (sg[0] * gsum[0] + sg[2] * gsum[1]),
                            -dt * (V / m) * (sg[2] * gsum[0] + sg[1] * gsum[1])};
    E.advance(dt, 1.0, M, err);
    const double worst = std::max(std::fabs(E.P.v[0] - vexp[0]), std::fabs(E.P.v[1] - vexp[1]));
    g_checks++;
    const bool ok = worst < 1e-16 && dropped > 0 && std::fabs(vexp[0]) > 1e-9;
    std::printf("%s: hand check, cutoff inside the step (%d node(s) dropped, v = %.6e %.6e, diff %.3e)\n",
                ok ? "ok  " : "FAIL", dropped, E.P.v[0], E.P.v[1], worst);
    if (!ok) g_fail++;
  }
  {  // A5 numbering: id = (i+1) + (nx+3)((j+1) + (ny+3)(k+1)), ghosts included
    Grid g;
    g.ndm = 3; g.h = 1.0; g.n[0] = 3; g.n[1] = 2; g.n[2] = 4;
    bool ok = true;
    for (int k = -1; k <= 5; k++)
      for (int j = -1; j <= 3; j++)
        for (int i = -1; i <= 4; i++) {
          const int ijk[3] = {i, j, k};
          const long want = (i + 1) + 6L * ((j + 1) + 5L * (k + 1));
          int back[3];
          g.extIJK(g.extId(ijk), back);
          ok = ok && g.publicId(g.extId(ijk)) == want && back[0] == i && back[1] == j && back[2] == k;
        }
    g_checks++;
    std::printf("%s: A5 node numbering (3D, ghosts included)\n", ok ? "ok  " : "FAIL");
    if (!ok) g_fail++;
  }
  auto cases = parse(argc > 1 ? argv[1] : "tests/mpm_cases.txt");
  for (const auto &c : cases) runCase(c);
  for (const auto &c : cases)
    for (const auto &cp : c.compares) {
      // compare <obs> [args...] <lt|gt> <other>
      const std::string ob = cp[1], other = cp.back(), op = cp[cp.size() - 2];
      g_checks++;
      if (!g_scalar[c.name].count(ob) || !g_scalar[other].count(ob)) {
        std::printf("FAIL: compare %s %s %s: value missing\n", c.name.c_str(), ob.c_str(), other.c_str());
        g_fail++;
        continue;
      }
      const double a = g_scalar[c.name][ob], b = g_scalar[other][ob];
      const bool ok = (op == "lt") ? (a < b) : (a > b);
      std::printf("%s: compare %s %s %.10e %s %s %.10e\n", ok ? "ok  " : "FAIL", c.name.c_str(),
                  ob.c_str(), a, op.c_str(), other.c_str(), b);
      if (!ok) g_fail++;
    }
  std::printf("[cases] parsed %d cases, %d checks, %d pending-oracle\n", (int)cases.size(), g_checks, g_pending);
  if (g_fail == 0) { std::printf("RESULT: PASS\n"); return 0; }
  std::printf("RESULT: FAIL (%d)\n", g_fail);
  return 1;
}
