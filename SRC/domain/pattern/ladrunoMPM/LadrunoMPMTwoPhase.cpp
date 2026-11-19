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
// LadrunoMPMTwoPhase — two-phase (saturated Biot u-p) particles in the MPM host, route E
// (mode A): WP-232l, frozen interface tests/mpm_twophase_interface.md (P1-P13, §7a).
//
// A two-phase epoch (any particle in a body given `mpm twoPhase`) makes every grid node
// carry ndf = ndm + 1, p last (P4); a node no two-phase particle reaches has its p fixed
// (stray). Each two-phase particle element adds LadrunoUP's blocks through the pure
// kernel (P3). The particle carries p, pd, pdd (P2, A1); the grid takes their
// volume-weighted projection and gives back an incremental update (P5). `mpm pfix`
// drains a face (P8), `mpm plate` is a rigid frictionless plate (P10), `mpm traction` a
// face traction on GIMP / linear grids (P10b). The general support (upSupport, P9) folds
// B-spline ghosts per component. Route S (mode F) is not built here.

#include "LadrunoMPMPattern.h"

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <DirectIntegrationAnalysis.h>
#include <TransientIntegrator.h>
#include <ConstraintHandler.h>
#include <AnalysisModel.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <profiler/ProfilerMacros.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>

namespace {
const char *faceName(int f)
{
  static const char *names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  return (f >= 0 && f < 6) ? names[f] : "?";
}

// int over [lo, hi] of the hat of half-width h centred at xa (P10b)
double hatIntegral(double xa, double h, double lo, double hi)
{
  auto g = [h](double u) { return u - u * u / (2.0 * h); };   // int_0^u (1 - s/h) ds
  double s = 0.0;
  const double a = std::max(lo, xa), b = std::min(hi, xa + h);           // right half
  if (b > a) s += g(b - xa) - g(a - xa);
  const double c = std::max(lo, xa - h), d = std::min(hi, xa);           // left half
  if (d > c) s += g(xa - c) - g(xa - d);
  return s;
}
}  // namespace

// ---------------------------------------------------------------------------
//  Settings (P1, P5, P6, P8, P10, P10b)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::setTwoPhaseBody(int body, const ladruno_mpm::TwoPhaseProps &tp, std::string &err)
{
  if (stepped || !eles.empty()) {
    err = "the two-phase bodies are fixed at the first step of an epoch; `mpm reset` first";
    return -1;
  }
  if (!(tp.Kf > 0.0)) { err = "-Kf must be > 0"; return -1; }
  if (!(tp.n > 0.0) || tp.n > 1.0) { err = "-poro must satisfy 0 < n <= 1"; return -1; }
  if (tp.rhoF < 0.0) { err = "-rhoF must be >= 0"; return -1; }
  if (!(tp.alpha > 0.0) || tp.alpha > 1.0) { err = "-alpha must satisfy 0 < alpha <= 1"; return -1; }
  if (tp.n > tp.alpha) {
    err = "porosity n = " + std::to_string(tp.n) + " exceeds Biot alpha = " + std::to_string(tp.alpha) +
          " (storage 1/Qbar needs 0 < n <= alpha <= 1, as LadrunoUP)";
    return -1;
  }
  for (int d = 0; d < 3; d++)
    if (tp.kbar[d] < 0.0) { err = "-perm values must be >= 0"; return -1; }
  twoPhase[body] = tp;
  nameBody(body);
  return 0;
}

int
LadrunoMPMPattern::setPfix(int face, double p0, std::string &err)
{
  if (face < 0 || face >= 2 * eng.ndm()) { err = "the face must be xmin, xmax, ymin, ymax[, zmin, zmax]"; return -1; }
  if (!std::isfinite(p0)) { err = "p0 must be finite"; return -1; }
  if (ladruno_mpm::isBSpline(eng.par.basis) && p0 != 0.0) {
    err = "a prescribed p0 != 0 on a B-spline face is refused in phase 6 (its image is odd about 0;"
          " interface P8); p0 = 0 (drained) is admitted";
    return -1;
  }
  if (stepped || !eles.empty()) { err = "pressure conditions are fixed at the first step; `mpm reset` first"; return -1; }
  eng.fc.pfixOn[face] = true;
  eng.fc.pfixVal[face] = p0;
  return 0;
}

int
LadrunoMPMPattern::setPlate(int face, double F, std::string &err)
{
  if (face < 0 || face >= 2 * eng.ndm()) { err = "the face must be xmin, xmax, ymin, ymax[, zmin, zmax]"; return -1; }
  if (!std::isfinite(F)) { err = "F must be finite"; return -1; }
  if (stepped || !eles.empty()) {
    // the force may change between calls; the plate itself is fixed at the first step
    if (!eng.fc.plateOn[face]) { err = "a plate is added before the first step; `mpm reset` first"; return -1; }
    eng.fc.plateF[face] = F;
    return 0;
  }
  if (eng.fc.fixMask[face] & (1u << (face / 2))) {
    err = std::string("face ") + faceName(face) + " fixes its normal component: a plate there is refused";
    return -1;
  }
  eng.fc.plateOn[face] = true;
  eng.fc.plateF[face] = F;
  return 0;
}

int
LadrunoMPMPattern::addTraction(const std::string &name, int face, const double *lo, const double *hi,
                               const double *t, std::string &err)
{
  const int nd = eng.ndm();
  if (face < 0 || face >= 2 * nd) { err = "the face must be xmin, xmax, ymin, ymax[, zmin, zmax]"; return -1; }
  if (ladruno_mpm::isBSpline(eng.par.basis)) {
    err = "`mpm traction` is GIMP / linear only in phase 6 (interface P10b); use `mpm plate`";
    return -1;
  }
  for (const Traction &x : tractions)
    if (x.name == name) { err = "a traction named '" + name + "' exists"; return -1; }
  Traction tr;
  tr.name = name;
  tr.face = face;
  for (int d = 0; d < 3; d++) {
    tr.lo[d] = (d < nd) ? lo[d] : 0.0;
    tr.hi[d] = (d < nd) ? hi[d] : 0.0;
    tr.t[d] = (d < nd) ? t[d] : 0.0;
    if (d < nd && !(tr.hi[d] >= tr.lo[d])) { err = "the box needs hi >= lo"; return -1; }
  }
  tractions.push_back(tr);
  return 0;
}

int
LadrunoMPMPattern::setKinematics(int small, std::string &err)
{
  if ((small ? 1 : 0) == eng.par.smallKin) return 0;
  if (stepped || !eles.empty()) {
    err = "the kinematics are fixed at the first step of an epoch; `mpm reset` first";
    return -1;
  }
  eng.par.smallKin = small ? 1 : 0;
  return 0;
}

int
LadrunoMPMPattern::setPressureBlend(double a, std::string &err)
{
  if (!(a >= 0.0) || a > 1.0) { err = "alpha_p must be in [0, 1]"; return -1; }
  eng.par.pBlend = a;
  return 0;
}

int
LadrunoMPMPattern::setParticlePressure(const std::vector<double> &v, std::string &err)
{
  if (stepped) { err = "initial conditions cannot change after the first step (`mpm reset` first)"; return -1; }
  const int N = eng.P.size();
  if ((int)v.size() != N && (int)v.size() != 3 * N) {
    err = "needs N = " + std::to_string(N) + " values (p) or 3N (p, then pd, then pdd), got " +
          std::to_string(v.size());
    return -1;
  }
  for (double x : v) if (!std::isfinite(x)) { err = "every value must be finite"; return -1; }
  for (int p = 0; p < N; p++) {
    eng.P.pw[p] = v[p];
    if ((int)v.size() == 3 * N) { eng.P.pwd[p] = v[N + p]; eng.P.pwdd[p] = v[2 * N + p]; }
  }
  return 0;
}

int
LadrunoMPMPattern::setParticleAcceleration(const std::vector<double> &a, std::string &err)
{
  if (stepped) { err = "initial conditions cannot change after the first step (`mpm reset` first)"; return -1; }
  const int nd = eng.ndm(), N = eng.P.size();
  if ((long)a.size() != (long)nd * N) {
    err = "needs ndm * N = " + std::to_string((long long)nd * N) + " values, got " + std::to_string(a.size());
    return -1;
  }
  for (double x : a) if (!std::isfinite(x)) { err = "every value must be finite"; return -1; }
  eng.P.a = a;
  return 0;
}

int
LadrunoMPMPattern::setBodyPorePressure(int body, double pv, std::string &err)
{
  if (stepped) { err = "initial conditions cannot change after the first step (`mpm reset` first)"; return -1; }
  if (!std::isfinite(pv)) { err = "p must be finite"; return -1; }
  int n = 0;
  for (int p = 0; p < eng.P.size(); p++)
    if (eng.P.body[p] == body) { eng.P.pw[p] = pv; n++; }
  return n;
}

int
LadrunoMPMPattern::twoPhaseInfo(int body, std::vector<double> &out, std::string &err)
{
  const auto it = twoPhase.find(body);
  if (it == twoPhase.end()) { err = "body " + std::to_string(body) + " is not two-phase"; return -1; }
  const ladruno_mpm::TwoPhaseProps &tp = it->second;
  const int nd = eng.ndm();
  out.clear();
  out.push_back(theMode == MODE_EXPLICIT ? 2.0 : 1.0);           // route: 1 E (mode A), 2 S (mode F)
  out.push_back(tp.Kf); out.push_back(tp.n); out.push_back(tp.rhoF);
  for (int d = 0; d < nd; d++) out.push_back(tp.kbar[d]);
  out.push_back(tp.alpha); out.push_back(tp.Ks);
  out.push_back(tp.stabMode); out.push_back(tp.stabValue);
  // alpha_st at the grid spacing, from the first particle of the body (0 if none or off)
  double ast = 0.0;
  if (tp.stabMode == 2) ast = tp.stabValue;
  else if (tp.stabMode == 1)
    for (int p = 0; p < eng.P.size(); p++)
      if (eng.P.body[p] == body && p < (int)mats.size() && mats[p] != 0) {
        const Matrix &D0 = mats[p]->getInitialTangent();
        const double Gd = (nd == 2) ? D0(2, 2) : D0(3, 3);
        const double Kd = D0(0, 1) + 2.0 * Gd / 3.0;
        const double den = Kd + 4.0 * Gd / 3.0;
        if (den > 0.0) ast = tp.stabValue * eng.grid.h * eng.grid.h / den;
        break;
      }
  out.push_back(ast);
  out.push_back(eng.par.pCut);                    // ruling c84f71fb3: route S's p cutoff
  out.push_back(piPerpMax);                       // A' (ruling f69e450dd): max |pi_perp| applied
  return 0;
}

bool
LadrunoMPMPattern::isTwoPhaseParticle(int p) const
{
  return !twoPhase.empty() && p >= 0 && p < eng.P.size() && twoPhase.count(eng.P.body[p]) != 0;
}

bool
LadrunoMPMPattern::anyTwoPhase(void) const
{
  if (twoPhase.empty()) return false;
  for (int p = 0; p < eng.P.size(); p++)
    if (twoPhase.count(eng.P.body[p])) return true;
  return false;
}

bool
LadrunoMPMPattern::anyPlate(void) const
{
  for (int f = 0; f < 6; f++) if (eng.fc.plateOn[f]) return true;
  return false;
}

const double *
LadrunoMPMPattern::mpmFluidDrive(int p)
{
  if ((int)fdBuf.size() < 3) fdBuf.assign(3, 0.0);
  const auto it = twoPhase.find(eng.P.body[p]);
  for (int d = 0; d < 3; d++) {
    double b = eng.par.gravity[d];
    if (it != twoPhase.end() && it->second.fluidBodyGiven) b = it->second.fluidBody[d];
    fdBuf[d] = lambda * b;
  }
  return fdBuf.data();
}

int
LadrunoMPMPattern::tagOfExt(long ext) const
{
  if (ext < 0) {
    const int f = (int)(-ext - 2);
    return (f >= 0 && f < 6) ? plateTag[f] : 0;
  }
  const auto it = liveTag.find(ext);
  return (it == liveTag.end()) ? 0 : it->second;
}

void
LadrunoMPMPattern::retirePlates(void)
{
  Domain *dom = this->getDomain();
  for (int f = 0; f < 6; f++) {
    if (plateTag[f] == 0) continue;
    if (dom != 0) {
      Node *n = dom->removeNode(plateTag[f]);
      delete n;
    }
    liveTags.erase(plateTag[f]);
    plateTag[f] = 0;
  }
}

// P7: "" or the refusal of a two-phase body / plate in this mode-A analysis
std::string
LadrunoMPMPattern::twoPhaseAdmission(const LadrunoMPMAnalysisRefs &A) const
{
  const bool tp = anyTwoPhase();
  for (int f = 0; f < 6; f++)
    if (eng.fc.plateOn[f] && (eng.fc.fixMask[f] & (1u << (f / 2))))
      return std::string("face ") + faceName(f) + " has a plate and fixes its normal component: refused";
  if (anyPlate() || tp) {
    if (eng.par.antiLocking)
      return std::string("`mpm antiLocking fbar` is refused with ") +
             (tp ? "a two-phase body (no tested composition with the u-p coupling; LadrunoUP's answer is"
                   " -formulation bbar; interface P7 (5))"
                 : "a plate (phase 6)");
  }
  if (!tp) return "";
  if (theMode == MODE_QUASISTATIC)
    return "mode quasiStatic is refused with a two-phase body: a static step has no S or Q^T rate terms,"
           " so it cannot consolidate (interface P7 (2)); use mode implicit with Newmark";
  if (A.ta != 0 && A.ta->getIntegrator() != 0 &&
      A.ta->getIntegrator()->getClassTag() != INTEGRATOR_TAGS_Newmark)
    return "a two-phase body runs Newmark only: HHT and GeneralizedAlpha stay refused until a test"
           " admits them (interface U5, P7 (3))";
  const int nv = eng.nv();
  for (int p = 0; p < eng.P.size(); p++) {
    if (!isTwoPhaseParticle(p)) continue;
    for (int k = 0; k < nv; k++)
      if (eng.P.sig0[nv * p + k] != 0.0)
        return "particle " + std::to_string((long long)eng.P.id[p]) + " of two-phase body " +
               std::to_string(eng.P.body[p]) + " carries an initial-stress offset: refused (Q2; interface"
               " P7 (6)); start from zero excess pressure, or stage under gravity";
  }
  bool nonHomogP = false;
  for (int f = 0; f < 6; f++) if (eng.fc.pfixOn[f] && eng.fc.pfixVal[f] != 0.0) nonHomogP = true;
  if (nonHomogP && A.handler != 0 && A.handler->getClassTag() == HANDLER_TAG_PlainHandler)
    return "a prescribed p0 != 0 needs `constraints Transformation`, `Penalty` or `Lagrange`: the Plain"
           " handler assumes every constraint homogeneous (PlainHandler.cpp:92)";
  return "";
}

// ---------------------------------------------------------------------------
//  One mode-A step's preparation in a two-phase epoch or with a plate / traction
//  (interface §1 "Step order, mode A", §2, §3, §7a Q3, Q5). The single-phase path
//  without them is prepareStep, unchanged.
// ---------------------------------------------------------------------------
// The supports, the u-side active set and renormalisation, fix 3's p-activity and S^p, and
// fix 2's affine coefficients and pi, at the particles' CURRENT positions and the host's current
// drive (interface 7b). Used by prepareStepUP and, after every accepted step, by the q + pi
// invariant (ruling on 89718fdad). No side effect on the model.
struct LadrunoMPMPattern::PSup {
  std::vector<long> supStart, node;
  std::vector<double> Sv, dSv, ScG, dScG, aff, SpF, dSpF, piS;
  std::vector<char> keep, inact;
  std::set<long> want, pAct;
  double mmax = 0.0;
};

int
LadrunoMPMPattern::buildPSup(PSup &s, std::string &err)
{
  ladruno_mpm::Particles &P = eng.P;
  const ladruno_mpm::Grid &G = eng.grid;
  const int N = P.size(), nd = eng.ndm(), nc = nd + 1;
  const bool up = anyTwoPhase();
  const bool bspl = ladruno_mpm::isBSpline(eng.par.basis);
  // ---- supports (P9) ----------------------------------------------------------------
  std::vector<long> &supStart = s.supStart, &node = s.node;
  supStart.assign(1, 0);
  node.clear();
  std::vector<double> &Sv = s.Sv, &dSv = s.dSv, &ScG = s.ScG, &dScG = s.dScG;
  Sv.clear(); dSv.clear(); ScG.clear(); dScG.clear();
  const int NA = nd * (nd + 1);              // fix 2: a_d, grad a_d per axis
  std::vector<double> &aff = s.aff;
  aff.assign((size_t)N * NA, 0.0);
  for (int p = 0; p < N; p++) {
    ladruno_mpm::upSupport(G, eng.par.basis, eng.fc, &P.x[nd * p], &P.lp[nd * p], node, Sv, dSv, ScG, dScG,
                           bspl ? &aff[(size_t)NA * p] : 0);
    supStart.push_back((long)node.size());
  }
  std::map<long, double> mass, tpMass;       // fix 3: the two-phase nodal mass, raw W
  double &mmax = s.mmax;
  mmax = 0.0;
  for (int p = 0; p < N; p++) {
    mmax = std::max(mmax, P.m[p]);
    const bool tpp = up && isTwoPhaseParticle(p);
    for (long k = supStart[p]; k < supStart[p + 1]; k++)
      if (node[k] >= 0) {
        mass[node[k]] += Sv[k] * P.m[p];
        if (tpp) tpMass[node[k]] += Sv[k] * P.m[p];
      }
  }
  // fix 2's coefficients follow every rescaling of the weights (B9's rule on a value + gradient)
  auto rescale = [&](double *a, double sum, const double *gsum) {
    for (int d = 0; d < nd; d++) {
      double &v = a[d * (nd + 1)];
      v /= sum;
      for (int q = 0; q < nd; q++) a[d * (nd + 1) + 1 + q] = (a[d * (nd + 1) + 1 + q] - v * gsum[q]) / sum;
    }
  };
  const double mcut = modeAMassCutoff() * mmax;
  const bool ghost = ghostGamma >= 0.0;
  std::vector<char> &keep = s.keep, &inact = s.inact;
  keep.assign(node.size(), 0);
  inact.assign(node.size(), 0);
  std::set<long> &want = s.want;
  want.clear();
  for (int p = 0; p < N; p++) {
    int nk = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (node[k] < 0) { keep[k] = 1; continue; }               // a plate unknown (P10)
      const bool act = ladruno_mpm::imageCutoffMass(G, eng.par.basis, eng.fc, node[k], mass[node[k]]) > mcut;
      if (act || ghost) {
        keep[k] = 1;
        inact[k] = act ? 0 : 1;
        want.insert(node[k]);
        nk++;
      }
    }
    if (nk == 0) {
      err = "particle " + std::to_string((long long)P.id[p]) + " has no active support node (its"
            " mass is below the cutoff everywhere); refused (W15: a zero-DOF element)";
      return -1;
    }
    // renormalize (§7a Q3): sums over the kept lattice nodes' W; every weight of every kept
    // entry -- W, each component, the plate's -- is rescaled by them (B9's rule)
    if (ghostGamma == EMPTY_RENORMALIZE) {
      double sum = 0.0, gsum[3] = {0.0, 0.0, 0.0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++)
        if (keep[k] && node[k] >= 0) {
          sum += Sv[k];
          for (int a = 0; a < nd; a++) gsum[a] += dSv[nd * k + a];
        }
      if (sum > 0.0 && bspl) rescale(&aff[(size_t)NA * p], sum, gsum);
      if (sum > 0.0)
        for (long k = supStart[p]; k < supStart[p + 1]; k++) {
          if (!keep[k]) continue;
          Sv[k] /= sum;
          for (int a = 0; a < nd; a++) dSv[nd * k + a] = (dSv[nd * k + a] - Sv[k] * gsum[a]) / sum;
          for (int c = 0; c < nc; c++) {
            double &sc = ScG[nc * k + c];
            sc /= sum;
            for (int a = 0; a < nd; a++)
              dScG[(nc * k + c) * nd + a] = (dScG[(nc * k + c) * nd + a] - sc * gsum[a]) / sum;
          }
        }
    }
  }

  // ---- fix 3 (P4 amended): p-activity on the two-phase nodal mass; S^p renormalised -------
  std::set<long> &pAct = s.pAct;
  pAct.clear();
  std::vector<double> &SpF = s.SpF, &dSpF = s.dSpF;  // per entry: the final p weight, gradient
  std::vector<double> &piS = s.piS;            // fix 2: pi at the current position, per particle
  piS.assign(N, 0.0);
  if (up) {
    // ruling 3d3699a50 in full: the mass test AND some two-phase particle with S^p > 0 there (a
    // ghost behind a drained-only face keeps its u mass, but its p weight is folded away)
    std::set<long> pReach;
    for (int p = 0; p < N; p++) {
      if (!isTwoPhaseParticle(p)) continue;
      for (long k = supStart[p]; k < supStart[p + 1]; k++)
        if (keep[k] && node[k] >= 0 && ScG[nc * k + nd] > 0.0) pReach.insert(node[k]);
    }
    const double pcut = eng.par.pCut * mmax;
    for (long x : want)
      if (pReach.count(x) && ladruno_mpm::imageCutoffMass(G, eng.par.basis, eng.fc, x, tpMass[x]) > pcut)
        pAct.insert(x);
    SpF.assign(node.size(), 0.0);
    dSpF.assign(nd * node.size(), 0.0);
    for (int p = 0; p < N; p++) {
      if (!isTwoPhaseParticle(p)) continue;
      bool allAct = true, anyAct = false;
      double sum = 0.0, gsum[3] = {0.0, 0.0, 0.0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        if (!keep[k] || node[k] < 0) continue;
        // a node whose p weight is folded away (S^p = 0: a ghost behind a drained-only face) is
        // not a dropped p node -- its W counts as kept, so it never triggers a renormalisation
        const bool pa = pAct.count(node[k]) != 0;
        if (pa) anyAct = true;
        if (pa || !(ScG[nc * k + nd] > 0.0)) {
          sum += Sv[k];
          for (int a = 0; a < nd; a++) gsum[a] += dSv[nd * k + a];
        } else {
          allAct = false;
        }
      }
      if (!anyAct || !(sum > 0.0)) {
        err = "two-phase particle " + std::to_string((long long)P.id[p]) + " has no p-active node (its"
              " two-phase nodal mass is below `mpm pressureCutoff` everywhere); refused (interface 7b,"
              " P4 amended)";
        return -1;
      }
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        if (!keep[k] || node[k] < 0 || !pAct.count(node[k])) continue;
        double sp = ScG[nc * k + nd];
        double gg[3];
        for (int a = 0; a < nd; a++) gg[a] = dScG[(nc * k + nd) * nd + a];
        if (!allAct) {
          sp /= sum;
          for (int a = 0; a < nd; a++) gg[a] = (gg[a] - sp * gsum[a]) / sum;
        }
        SpF[k] = sp;
        for (int a = 0; a < nd; a++) dSpF[nd * k + a] = gg[a];
      }
      if (!allAct && bspl) rescale(&aff[(size_t)NA * p], sum, gsum);
      // fix 2: pi at the projection = rhoF sum_d drive_d a_d (drive = the host's lambda b_f)
      if (bspl) {
        const auto it = twoPhase.find(P.body[p]);
        const double *drv = mpmFluidDrive(p);
        for (int d = 0; d < nd; d++) piS[p] += it->second.rhoF * drv[d] * aff[(size_t)NA * p + d * (nd + 1)];
      }
    }
  }
  return 0;
}

int
LadrunoMPMPattern::prepareStepUP(double hStep, std::string &err)
{
  Domain *dom = this->getDomain();
  ladruno_mpm::Particles &P = eng.P;
  const ladruno_mpm::Grid &G = eng.grid;
  const int N = P.size(), nd = eng.ndm(), nc = nd + 1;
  const bool up = anyTwoPhase();
  const int ndfN = up ? nd + 1 : nd;
  const double t = dom->getCurrentTime();
  const bool bspl = ladruno_mpm::isBSpline(eng.par.basis);
  OPS_PROFILE_SCOPE("mpm.a.prepare");

  PSup s;
  if (buildPSup(s, err) < 0) return -1;
  std::vector<long> &supStart = s.supStart, &node = s.node;
  std::vector<double> &Sv = s.Sv, &dSv = s.dSv, &ScG = s.ScG, &dScG = s.dScG;
  const int NA = nd * (nd + 1);
  std::vector<double> &aff = s.aff, &SpF = s.SpF, &dSpF = s.dSpF, &piS = s.piS;
  std::vector<char> &keep = s.keep, &inact = s.inact;
  std::set<long> &want = s.want, &pAct = s.pAct;
  std::map<long, std::vector<double>> pProj;  // P5a: x^ for p, pd, pdd per p-active node
  if (up) {
    // ---- P5a: the consistent L2 projection M x^ = P^T V x_part over the p-active nodes ----
    std::map<long, int> pidx;
    for (long x : pAct) { const int i = (int)pidx.size(); pidx[x] = i; }
    const int np = (int)pidx.size();
    std::vector<std::map<int, double>> M(np);
    std::vector<double> rhs(3 * (size_t)np, 0.0);
    for (int p = 0; p < N; p++) {
      if (!isTwoPhaseParticle(p)) continue;
      const double V = P.V[p];
      // the q + pi invariant: the projection fits q = p - pi_c (at the epoch's first projection,
      // p_init - pi: A' then makes q representable)
      if ((int)piC.size() != N) piC.assign(N, 0.0);
      const double xv[3] = {P.pw[p] - (piInitApplied ? piC[p] : piS[p]), P.pwd[p], P.pwdd[p]};
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        if (SpF[k] == 0.0 || node[k] < 0) continue;
        const int a = pidx[node[k]];
        for (int q = 0; q < 3; q++) rhs[3 * (size_t)a + q] += V * SpF[k] * xv[q];
        for (long l = supStart[p]; l < supStart[p + 1]; l++) {
          if (SpF[l] == 0.0 || node[l] < 0) continue;
          M[a][pidx[node[l]]] += V * SpF[k] * SpF[l];
        }
      }
    }
    double dmax = 0.0;
    std::vector<double> diag(np, 0.0);
    for (int a = 0; a < np; a++) {
      const auto it = M[a].find(a);
      diag[a] = (it != M[a].end()) ? it->second : 0.0;
      dmax = std::max(dmax, diag[a]);
    }
    const std::string singular =
        "the consistent pressure projection (interface 7b, P5a) is singular or did not reach rel 1e-14:"
        " the particles no longer span the p-active nodes' shape functions (e.g. one particle per cell"
        " on a B-spline grid); step refused, nothing changed";
    for (int a = 0; a < np; a++)
      if (!(diag[a] > 0.0)) { err = singular; return -1; }
    // Jacobi-PCG to rel 1e-14, with the frozen singularity guard (wp-briefs "Freeze of P5a"):
    // false when it does not converge within 4 np + 1000 iterations, or a curvature
    // d^T M d <= 1e-14 max diag(M) d^T d. A zero right side returns 0.
    auto pcg = [&](const std::vector<double> &b, std::vector<double> &xs) -> bool {
      std::vector<double> r(b), z(np), dvec(np), Ad(np);
      xs.assign(np, 0.0);
      double r0 = 0.0;
      for (int a = 0; a < np; a++) r0 += r[a] * r[a];
      r0 = std::sqrt(r0);
      if (!(r0 > 0.0)) return true;
      double rz = 0.0;
      for (int a = 0; a < np; a++) { z[a] = r[a] / diag[a]; dvec[a] = z[a]; rz += r[a] * z[a]; }
      const int maxit = 4 * np + 1000;
      for (int it = 0; it < maxit; it++) {
        double dAd = 0.0, dd = 0.0;
        for (int a = 0; a < np; a++) {
          double s2 = 0.0;
          for (const auto &kv : M[a]) s2 += kv.second * dvec[kv.first];
          Ad[a] = s2;
          dAd += dvec[a] * s2;
          dd += dvec[a] * dvec[a];
        }
        if (!(dAd > 1e-14 * dmax * dd)) return false;     // a (near-)singular direction
        const double al = rz / dAd;
        double rn = 0.0;
        for (int a = 0; a < np; a++) {
          xs[a] += al * dvec[a];
          r[a] -= al * Ad[a];
          rn += r[a] * r[a];
        }
        if (std::sqrt(rn) <= 1e-14 * r0) return true;
        double rzn = 0.0;
        for (int a = 0; a < np; a++) { z[a] = r[a] / diag[a]; rzn += r[a] * z[a]; }
        const double be = rzn / rz;
        rz = rzn;
        for (int a = 0; a < np; a++) dvec[a] = z[a] + be * dvec[a];
      }
      return false;
    };
    // guard (a) amended (ruling on 7b7cdc9cf): the probe, b_k = sin(1.7 k + 0.3) over the
    // p-active rows in storage order. The data's right side P^T V x lies in M's range, so its own
    // PCG converges on a singular M; the probe does not lie in the range of a singular M.
    std::vector<double> xs, bq(np);
    if (np > 0) {
      for (int a = 0; a < np; a++) bq[a] = std::sin(1.7 * a + 0.3);
      if (!pcg(bq, xs)) { err = singular + " (the probe solve)"; return -1; }
    }
    // A' (rulings f69e450dd and 89718fdad, interface 7b): at the epoch's first projection, each
    // two-phase particle's p gains pi's own invisible part, pi_perp = pi - P M^-1 P^T V pi (same
    // p-active set, M, S^p), i.e. q = p_init - P M^-1 P^T V pi(0). pd, pdd unchanged. A drive that
    // ramps from 0 has pi(0) = 0: nothing. P^T V pi_perp = 0: the projection below is unchanged.
    std::vector<double> piPerp;
    bool anyPi = false;
    for (int p = 0; p < N; p++) anyPi = anyPi || piS[p] != 0.0;
    if (!piInitApplied && anyPi && np > 0) {
      std::vector<double> bpi(np, 0.0), y;
      for (int p = 0; p < N; p++) {
        if (!isTwoPhaseParticle(p) || piS[p] == 0.0) continue;
        for (long k = supStart[p]; k < supStart[p + 1]; k++)
          if (SpF[k] != 0.0 && node[k] >= 0) bpi[pidx[node[k]]] += P.V[p] * SpF[k] * piS[p];
      }
      if (!pcg(bpi, y)) { err = singular + " (the pi_perp solve, A')"; return -1; }
      piPerp.assign(N, 0.0);
      for (int p = 0; p < N; p++) {
        if (!isTwoPhaseParticle(p)) continue;
        double s = 0.0;
        for (long k = supStart[p]; k < supStart[p + 1]; k++)
          if (SpF[k] != 0.0 && node[k] >= 0) s += SpF[k] * y[pidx[node[k]]];
        piPerp[p] = piS[p] - s;
      }
    }
    std::vector<double> sol(3 * (size_t)np, 0.0);
    for (int q = 0; q < 3 && np > 0; q++) {
      for (int a = 0; a < np; a++) bq[a] = rhs[3 * (size_t)a + q];
      if (!pcg(bq, xs)) { err = singular; return -1; }
      for (int a = 0; a < np; a++) sol[3 * (size_t)a + q] = xs[a];
    }
    for (const auto &kv : pidx)
      pProj[kv.first] = {sol[3 * (size_t)kv.second], sol[3 * (size_t)kv.second + 1],
                         sol[3 * (size_t)kv.second + 2]};
    // applied once every solve has succeeded (a refused step changes nothing)
    if (!piInitApplied) {          // the epoch's first projection: q = p_init - (P M^-1 P^T V pi)
      if (!piPerp.empty())
        for (int p = 0; p < N; p++) {
          P.pw[p] += piPerp[p];
          piPerpMax = std::max(piPerpMax, std::fabs(piPerp[p]));
        }
      piC = piS;
      piInitApplied = true;
    } else {                       // p <- q + pi(current x, drive): the change of pi since pi_c
      for (int p = 0; p < N; p++) {
        if (!isTwoPhaseParticle(p)) continue;
        P.pw[p] += piS[p] - piC[p];
        piC[p] = piS[p];
      }
    }
  }

  {
    const std::string why = gridvelGhostConflict(t);
    if (!why.empty()) { err = why; return -1; }
  }
  // ---- the grid conditions, in lattice terms (plates: plateId(f)) ----------------------
  std::map<std::pair<long, int>, std::pair<double, bool>> spec;
  auto setSpec = [&](long x, int a, double v, bool constant) -> bool {
    const auto key = std::make_pair(x, a);
    const auto it = spec.find(key);
    if (it != spec.end()) return it->second.first == v;
    spec[key] = std::make_pair(v, constant);
    return true;
  };
  std::vector<std::vector<long>> gvExt(gridVels.size());
  const double tol = 1e-9 * G.h;
  for (long x : want) {
    int ijk[3];
    G.extIJK(x, ijk);
    for (int f = 0; f < 2 * nd; f++) {
      const int d = f / 2;
      const bool on = (f % 2 == 0) ? ijk[d] <= 0 : ijk[d] >= G.n[d];
      if (!on) continue;
      for (int a = 0; a < nd; a++)
        if (eng.fc.fixMask[f] & (1u << a)) setSpec(x, a, 0.0, true);
      if (up && eng.fc.pfixOn[f] && !setSpec(x, nd, eng.fc.pfixVal[f], true)) {
        err = "two `pfix` faces prescribe different p on one node (a corner): refused";
        return -1;
      }
    }
  }
  for (size_t gi = 0; gi < gridVels.size(); gi++) {
    const GridVel &gv = gridVels[gi];
    if (!gv.presentAt(t)) continue;
    double lo[3], hi[3];
    gv.boxAt(t, lo, hi);
    const double fac = gv.stepFactor(t, hStep);
    for (long x : want) {
      int ijk[3];
      G.extIJK(x, ijk);
      bool in = true;
      for (int d = 0; d < nd; d++) {
        const double c = G.x0[d] + ijk[d] * G.h;
        in = in && c >= lo[d] - tol && c <= hi[d] + tol;
      }
      if (!in) continue;
      gvExt[gi].push_back(x);
      for (int a = 0; a < nd; a++) {
        if (gv.freeMask & (1u << a)) continue;
        double cx[3] = {0, 0, 0};
        for (int d = 0; d < nd; d++) cx[d] = G.x0[d] + ijk[d] * G.h;
        const double val = gv.valueAt(a, cx, nd, fac);
        if (!setSpec(x, a, val, val == 0.0)) {
          err = "gridvel '" + gv.name + "' and another condition prescribe different values on one"
                " node component (a malformed case)";
          return -1;
        }
      }
    }
  }
  // P10: a plate face row's normal component is stray (its weight moved to U); the plate
  // unknown's other components (and p) are fixed
  for (int f = 0; f < 2 * nd; f++) {
    if (!eng.fc.plateOn[f]) continue;
    const int d = f / 2;
    for (long x : want) {
      int ijk[3];
      G.extIJK(x, ijk);
      if (ijk[d] != ((f % 2 == 0) ? 0 : G.n[d])) continue;
      if (!setSpec(x, d, 0.0, true)) {
        err = std::string("a gridvel prescribes the normal component of plate face ") + faceName(f) +
              "'s row: refused";
        return -1;
      }
    }
    for (int a = 0; a < nd; a++) if (a != d) setSpec(ladruno_mpm::plateId(f), a, 0.0, true);
    if (up) setSpec(ladruno_mpm::plateId(f), nd, 0.0, true);
  }
  // P4 amended (fix 3): every node that is not p-active has its p fixed (stray)
  if (up)
    for (long x : want)
      if (!pAct.count(x)) setSpec(x, nd, 0.0, true);
  std::vector<SpecRow> newSpec;
  for (const auto &kv : spec)
    newSpec.push_back(SpecRow{kv.first.first, kv.first.second, kv.second.first, kv.second.second});

  // ---- the node diff -------------------------------------------------------------------
  bool nodeSetChanged = (reseedMode == 1);
  if (!nodeSetChanged) {
    if (want.size() != liveTag.size()) nodeSetChanged = true;
    else
      for (long x : want)
        if (liveTag.find(x) == liveTag.end()) { nodeSetChanged = true; break; }
  }
  bool plateMissing = false;
  for (int f = 0; f < 2 * nd; f++) if (eng.fc.plateOn[f] && plateTag[f] == 0) plateMissing = true;
  const bool specChanged = nodeSetChanged || plateMissing || !(newSpec == lastSpec);
  if (specChanged) removeHostSPs();
  auto ensureBases = [&]() {
    if (basesSet) return;
    int maxN = 0, maxE = 0;
    NodeIter &ni = dom->getNodes();
    Node *n;
    while ((n = ni()) != 0) maxN = std::max(maxN, n->getTag());
    ElementIter &ei = dom->getElements();
    Element *e;
    while ((e = ei()) != 0) maxE = std::max(maxE, e->getTag());
    nodeBase = nextNodeTag = maxN + 1;
    eleBase = maxE + 1;
    basesSet = true;
  };
  auto makeNode = [&](double c0, double c1, double c2) -> Node * {
    return (nd == 2) ? new Node(nextNodeTag, ndfN, c0, c1) : new Node(nextNodeTag, ndfN, c0, c1, c2);
  };
  if (nodeSetChanged) {
    std::vector<long> drop;
    for (const auto &kv : liveTag)
      if (reseedMode == 1 || want.count(kv.first) == 0) drop.push_back(kv.first);
    for (long x : drop) retireNode(x);
    ensureBases();
    for (long x : want) {
      if (liveTag.count(x)) continue;
      if (nextNodeTag == INT_MAX) { err = "grid node tags would overflow int (W5: never reused)"; return -2; }
      int ijk[3];
      G.extIJK(x, ijk);
      Node *n = makeNode(G.x0[0] + ijk[0] * G.h, G.x0[1] + ijk[1] * G.h, G.x0[2] + ijk[2] * G.h);
      if (!dom->addNode(n)) {
        delete n;
        err = "could not add grid node " + std::to_string(nextNodeTag) + " (a tag in use? see `mpm tagBase`)";
        return -2;
      }
      liveTags.insert(nextNodeTag);
      liveTag[x] = nextNodeTag++;
    }
  }
  for (int f = 0; f < 2 * nd; f++) {
    if (!eng.fc.plateOn[f] || plateTag[f] != 0) continue;
    ensureBases();
    if (nextNodeTag == INT_MAX) { err = "grid node tags would overflow int (W5: never reused)"; return -2; }
    double c[3] = {0, 0, 0};
    for (int a = 0; a < nd; a++) c[a] = G.x0[a] + 0.5 * G.n[a] * G.h;
    c[f / 2] = G.x0[f / 2] + ((f % 2) ? G.n[f / 2] * G.h : 0.0);
    Node *n = makeNode(c[0], c[1], c[2]);
    if (!dom->addNode(n)) {
      delete n;
      err = "could not add the plate node " + std::to_string(nextNodeTag) + " (a tag in use?)";
      return -2;
    }
    liveTags.insert(nextNodeTag);
    plateTag[f] = nextNodeTag++;
  }

  // ---- P2G (§7a Q5): m, v, a and the two-phase projections ------------------------------
  const bool qs = (theMode == MODE_QUASISTATIC);
  const int R0 = 1 + 2 * nd, RS = R0 + 4;    // r = [m, m v, m a, Vp p, Vp pd, Vp pdd, Vp]
  std::map<long, std::vector<double>> mv;
  for (long x : want) mv[x].assign(RS, 0.0);
  for (int f = 0; f < 2 * nd; f++) if (eng.fc.plateOn[f]) mv[ladruno_mpm::plateId(f)].assign(RS, 0.0);
  for (int p = 0; p < N; p++) {
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (!keep[k]) continue;
      std::vector<double> &r = mv[node[k]];
      r[0] += Sv[k] * P.m[p];
      if (!qs)
        for (int a = 0; a < nd; a++) {
          const double wa = ScG[nc * k + a] * P.m[p];
          r[1 + a] += wa * P.v[nd * p + a];
          r[1 + nd + a] += wa * P.a[nd * p + a];
        }
    }
  }
  // P5a: the projected p, pd, pdd (the consistent solve above); 0 off the p-active set
  for (auto &kv : mv) {
    const auto it = pProj.find(kv.first);
    if (it == pProj.end()) continue;
    for (int q = 0; q < 3; q++) kv.second[R0 + q] = it->second[q];
    kv.second[R0 + 3] = 1.0;
  }
  for (auto &kv : mv) {
    std::vector<double> &r = kv.second;
    for (int a = 0; a < nd; a++) {
      r[1 + a] = (r[0] > 0.0) ? r[1 + a] / r[0] : 0.0;
      r[1 + nd + a] = (r[0] > 0.0) ? r[1 + nd + a] / r[0] : 0.0;
    }
    for (int q = 0; q < 3; q++) r[R0 + q] = (r[R0 + 3] != 0.0) ? r[R0 + q] / r[R0 + 3] : 0.0;
  }
  if (bspl)                                    // WP-232i B6: odd components on the symmetry plane
    for (auto &kv : mv) {
      if (kv.first < 0) continue;
      int ijk[3];
      G.extIJK(kv.first, ijk);
      for (int d = 0; d < nd; d++) {
        unsigned m = 0u;
        if (ijk[d] == 0) m |= eng.fc.fixMask[2 * d];
        if (ijk[d] == G.n[d]) m |= eng.fc.fixMask[2 * d + 1];
        for (int a = 0; a < nd; a++)
          if (m & (1u << a)) { kv.second[1 + a] = 0.0; kv.second[1 + nd + a] = 0.0; }
      }
    }
  Vector dv(ndfN), vel(ndfN), acc(ndfN);
  for (auto &kv : mv) {
    Node *n = dom->getNode(tagOfExt(kv.first));
    if (n == 0) { err = "internal: grid node missing"; return -2; }
    const std::vector<double> &r = kv.second;
    dv.Zero();
    for (int a = 0; a < nd; a++) {
      vel(a) = r[1 + a];
      acc(a) = r[1 + nd + a];
      auto it = spec.find(std::make_pair(kv.first, a));
      if (it != spec.end()) {
        vel(a) = qs ? 0.0 : it->second.first;
        acc(a) = 0.0;
      }
    }
    if (up) {
      dv(nd) = r[R0];
      vel(nd) = r[R0 + 1];
      acc(nd) = r[R0 + 2];
      auto it = spec.find(std::make_pair(kv.first, nd));
      if (it != spec.end()) {
        dv(nd) = it->second.first;
        vel(nd) = 0.0;
        acc(nd) = 0.0;
      }
    }
    n->setTrialDisp(dv);
    n->setTrialVel(vel);
    n->setTrialAccel(acc);
    n->commitState();
  }

  // ---- elements: created, or re-pointed in place, then the u-p blocks ---------------------
  if (eles.empty() && createElements(err) < 0) return -2;
  bool connChanged = false;
  for (int p = 0; p < N; p++) {
    int nk = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) nk += keep[k];
    ID tags(nk);
    std::vector<double> S(nk), dS(nd * nk), vH(nd * nk), aH(nd * nk), Sc(nd * nk), dSc(nd * nd * nk);
    std::vector<double> Sp(nk), dSp(nd * nk), pRef(3 * nk);
    std::vector<char> ina(nk);
    int j = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (!keep[k]) continue;
      tags(j) = tagOfExt(node[k]);
      S[j] = Sv[k];
      ina[j] = inact[k];
      const std::vector<double> &r = mv[node[k]];
      for (int a = 0; a < nd; a++) {
        dS[nd * j + a] = dSv[nd * k + a];
        vH[nd * j + a] = r[1 + a];
        aH[nd * j + a] = r[1 + nd + a];
        Sc[nd * j + a] = ScG[nc * k + a];
        for (int e = 0; e < nd; e++) dSc[(nd * j + a) * nd + e] = dScG[(nc * k + a) * nd + e];
      }
      Sp[j] = SpF.empty() ? 0.0 : SpF[k];                  // fix 3: the p-active, renormalised weight
      for (int e = 0; e < nd; e++) dSp[nd * j + e] = SpF.empty() ? 0.0 : dSpF[nd * k + e];
      for (int q = 0; q < 3; q++) pRef[3 * j + q] = r[R0 + q];
      j++;
    }
    LadrunoMPMPoint *e = eles[p];
    if (e == 0) {
      err = "particle element " + std::to_string(p) + " is gone (removed by hand?)";
      return -2;
    }
    if (!e->connectivityEquals(tags)) connChanged = true;
    if (e->repoint(tags, S, dS, ina, vH, aH, Sc, dSc, 0, 0, G.h) < 0) {
      err = "re-pointing failed";
      return -2;
    }
    if (up) {
      const auto it = twoPhase.find(P.body[p]);
      const ladruno_mpm::TwoPhaseProps *tp = (it != twoPhase.end()) ? &it->second : 0;
      if (e->setTwoPhase(ndfN, tp, G.h, Sp, dSp, pRef, bspl ? &aff[(size_t)NA * p] : 0, piS[p]) < 0) {
        err = "the two-phase blocks of particle " + std::to_string((long long)P.id[p]) + " failed";
        return -2;
      }
    }
  }
  for (int p = 0; p < N; p++) {
    if (eles[p]->getDomain() != 0) continue;
    eles[p]->setRayleigh(rayAlphaM, rayBetaK);
    if (!dom->addElement(eles[p])) {
      err = "could not add particle element " + std::to_string(eles[p]->getTag()) +
            " (a tag in use? see `mpm tagBase`)";
      return -2;
    }
    connChanged = true;
  }
  if (connChanged) nConnChanges++;

  // ---- the grid conditions as this pattern's SPs ------------------------------------------
  if (specChanged) {
    for (const SpecRow &r : newSpec) {
      SP_Constraint *sp = new SP_Constraint(tagOfExt(r.ext), r.dof, r.value, r.constant);
      if (!this->addSP_Constraint(sp)) {
        delete sp;
        err = "could not add a grid-condition SP";
        return -2;
      }
      hostSPTags.push_back(sp->getTag());
    }
    lastSpec = newSpec;
  }
  for (size_t gi = 0; gi < gridVels.size(); gi++) {
    gridVels[gi].nodes.clear();
    for (long x : gvExt[gi]) gridVels[gi].nodes.push_back(liveTag[x]);
  }
  // P10b: the tractions' consistent nodal forces on this step's face-row nodes (x lambda in
  // applyLoad): f_a = t prod_{e != d} int_seg hat_a (thickness 1 in 2D)
  tractionLoads.clear();
  for (const Traction &tr : tractions) {
    const int d = tr.face / 2, row = (tr.face % 2) ? G.n[d] : 0;
    for (long x : want) {
      int ijk[3];
      G.extIJK(x, ijk);
      if (ijk[d] != row) continue;
      double w = 1.0;
      for (int e = 0; e < nd; e++)
        if (e != d) w *= hatIntegral(G.x0[e] + ijk[e] * G.h, G.h, tr.lo[e], tr.hi[e]);
      if (w == 0.0) continue;
      std::vector<double> fv(ndfN, 0.0);
      for (int a = 0; a < nd; a++) fv[a] = w * tr.t[a];
      tractionLoads.push_back(std::make_pair(liveTag[x], fv));
    }
  }

  // ---- Domain::domainChange() from C++, unconditionally (W1, W4) -------------------------
  if (!(reseedMode == 2 && !nodeSetChanged && !connChanged && !specChanged)) {
    dom->domainChange();
    nDomainChanges++;
  }
  return 0;
}

// The q + pi invariant after an accepted step (ruling on 89718fdad): G2P moved the particles
// and added the FLIP increment of q; p is brought to q + pi at the new positions and the current
// drive, so output and queries see p_p = q_p + pi(x_p, drive). Only B-spline grids carry pi.
void
LadrunoMPMPattern::twoPhaseAfterStep(void)
{
  if (!piInitApplied || !ladruno_mpm::isBSpline(eng.par.basis) || !anyTwoPhase()) return;
  PSup s;
  std::string err;
  if (buildPSup(s, err) < 0) {
    opserr << "WARNING mpm analyze -- the pore-pressure term pi could not be re-evaluated after the"
              " step (" << err.c_str() << "); the next projection does it\n";
    return;
  }
  const int N = eng.P.size();
  if ((int)piC.size() != N) piC.assign(N, 0.0);
  for (int p = 0; p < N; p++) {
    if (!isTwoPhaseParticle(p)) continue;
    eng.P.pw[p] += s.piS[p] - piC[p];
    piC[p] = s.piS[p];
  }
}

// ---------------------------------------------------------------------------
//  Debug hooks (§7a)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::debugParticleBlocks(int p, const std::string &path, const double *scales,
                                       std::string &err)
{
  if (theMode == MODE_EXPLICIT) { err = "particleBlocks needs mode implicit (mode A)"; return -1; }
  if (p < 0 || p >= eng.P.size()) { err = "particle index out of range"; return -1; }
  if (elementsRestored && relinkElements(err) < 0) return -1;
  armedWindow = true;
  const int prc = prepareStep(0.0, err);
  if (prc < 0) {
    armedWindow = false;
    if (prc == -2) aborted = true;
    return -1;
  }
  inStep = true;
  stepT0 = this->getDomain()->getCurrentTime();
  int rc = -1;
  std::ofstream os(path.c_str());
  if (!os) err = "cannot open '" + path + "'";
  else if (eles[p] == 0) err = "the particle element is gone";
  else if (eles[p]->writeBlocks(scales, os, lambda, (long long)eng.P.id[p]) < 0)
    err = "the element could not form its trial state (see the message above)";
  else rc = 0;
  inStep = false;
  armedWindow = false;
  return rc;
}

int
LadrunoMPMPattern::debugAssemble(const std::string &path, double dt, const LadrunoMPMAnalysisRefs &A,
                                 std::string &err)
{
  if (theMode != MODE_IMPLICIT || A.ta == 0 || A.ta->getIntegrator() == 0) {
    err = "assemble needs mode implicit and a Transient analysis";
    return -1;
  }
  if (!(dt > 0.0)) { err = "dt must be > 0"; return -1; }
  if (admitAnalysis(A, err) < 0) return -1;
  if (elementsRestored && relinkElements(err) < 0) return -1;
  Domain *dom = this->getDomain();
  armedWindow = true;
  const int prc = prepareStep(dt, err);
  if (prc < 0) {
    armedWindow = false;
    if (prc == -2) aborted = true;
    return -1;
  }
  stepT0 = dom->getCurrentTime();
  inStep = true;
  TransientIntegrator *ti = A.ta->getIntegrator();
  int rc = 0;
  if (A.ta->checkDomainChange() < 0 || ti->newStep(dt) < 0) {
    err = "the analysis could not set up the step (domainChanged / newStep)";
    rc = -1;
  }
  if (rc == 0) {
    AnalysisModel *am = A.ta->getModel();
    const int n = am->getNumEqn();
    std::vector<double> Ad((size_t)n * n, 0.0), bd(n, 0.0);
    std::vector<std::vector<double>> lab(n);
    FE_EleIter &fes = am->getFEs();
    FE_Element *fe;
    while ((fe = fes()) != 0) {
      const ID &id = fe->getID();
      const Matrix &k = fe->getTangent(ti);
      for (int i = 0; i < id.Size(); i++) {
        if (id(i) < 0) continue;
        for (int j = 0; j < id.Size(); j++)
          if (id(j) >= 0) Ad[(size_t)id(i) * n + id(j)] += k(i, j);
      }
      const Vector &r = fe->getResidual(ti);
      for (int i = 0; i < id.Size(); i++) if (id(i) >= 0) bd[id(i)] += r(i);
    }
    DOF_GrpIter &dgs = am->getDOFs();
    DOF_Group *dg;
    while ((dg = dgs()) != 0) {
      const ID &id = dg->getID();
      const Matrix &k = dg->getTangent(ti);
      const Vector &r = dg->getUnbalance(ti);
      Node *nd = dom->getNode(dg->getNodeTag());
      for (int i = 0; i < id.Size(); i++) {
        if (id(i) < 0) continue;
        bd[id(i)] += r(i);
        for (int j = 0; j < id.Size(); j++)
          if (id(j) >= 0) Ad[(size_t)id(i) * n + id(j)] += k(i, j);
        std::vector<double> l(1, (double)dg->getNodeTag());
        if (nd != 0) {
          const Vector &c = nd->getCrds();
          for (int a = 0; a < c.Size(); a++) l.push_back(c(a));
        }
        l.push_back((double)i);
        lab[id(i)] = l;
      }
    }
    std::ofstream os(path.c_str());
    if (!os) {
      err = "cannot open '" + path + "'";
      rc = -1;
    } else {
      os << std::setprecision(17) << "{\"n\": " << n << ",\n \"A\": [";
      for (int i = 0; i < n; i++) {
        os << (i ? ", " : "") << "[";
        for (int j = 0; j < n; j++) os << (j ? ", " : "") << Ad[(size_t)i * n + j];
        os << "]";
      }
      os << "],\n \"b\": [";
      for (int i = 0; i < n; i++) os << (i ? ", " : "") << bd[i];
      os << "],\n \"dofs\": [";
      for (int i = 0; i < n; i++) {
        os << (i ? ", " : "") << "[";
        for (size_t q = 0; q < lab[i].size(); q++) os << (q ? ", " : "") << lab[i][q];
        os << "]";
      }
      os << "]}\n";
    }
  }
  // nothing is committed: the trial state goes back
  ti->revertToLastStep();
  dom->revertToLastCommit();
  inStep = false;
  armedWindow = false;
  return rc;
}
