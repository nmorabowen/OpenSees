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
// LadrunoMPMImplicit — mode A of the MPM host (WP-232f, ADR 232 §5.3 / §5.4,
// frozen interface tests/mpm_implicit_interface.md): particles as
// LadrunoMPMPoint elements, active support nodes as Domain nodes, the user's
// Transient / Static analysis solving each step.
//
// ONE STEP (interface §3): supports at the committed positions; the active set
// (m_a > c_m max m_p; exclude or ghost, D2.3); node diff (retire / create, a
// retired tag never reused, W5); P2G onto the nodes (disp 0, projected v, a,
// committed); elements created or re-pointed in place; the grid conditions as
// this pattern's SP_Constraints; Domain::domainChange() from C++,
// unconditionally (W1, W4); analyzeStep(dt) (W10) or StaticAnalysis::analyze(1,
// false) (W16); the batched G2P in the last particle element's commitState
// (W11, amendment 1). A step that fails before any commit is cut and retried
// with nothing to undo; one that fails after a commit began aborts (W17).

#include "LadrunoMPMPattern.h"

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <LoadPatternIter.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <LadrunoMaterialStatus.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <PFEMAnalysis.h>
#include <TransientIntegrator.h>
#include <StaticIntegrator.h>
#include <LoadControl.h>
#include <ConstraintHandler.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <set>

// ---------------------------------------------------------------------------
//  LadrunoMPMPointHost
// ---------------------------------------------------------------------------
const double *
LadrunoMPMPattern::mpmBodyForce(int p)
{
  const int nd = eng.ndm();
  if ((int)bfBuf.size() < 3) bfBuf.assign(3, 0.0);
  const double mp = eng.P.m[p];
  for (int a = 0; a < nd; a++) {
    double f = mp * eng.par.gravity[a];
    if ((int)pforce.size() >= nd * (p + 1)) f += pforce[nd * p + a];
    bfBuf[a] = lambda * f;
  }
  return bfBuf.data();
}

bool
LadrunoMPMPattern::mpmTrialAdmissible(int p, const double *x)
{
  const int nd = eng.ndm();
  for (int d = 0; d < nd; d++)
    if (!ladruno_mpm::admissibleAxis(eng.grid, eng.par.basis, d, x[d], eng.P.lp[nd * p + d]))
      return false;
  return true;
}

// The batched G2P (interface D1.9, ADR §5.3 amendment 1): the N-th particle
// commit of this Domain::commit runs G2P for every particle — unless a commit
// refusal is pending (a read that does not clear: Domain::commit acts on it
// right after the element loop), in which case nothing moves (W17).
void
LadrunoMPMPattern::mpmParticleCommitted(int)
{
  committedCount++;
  if (committedCount != armedCount) return;
  if (ladrunoPendingCommitRefusals() > 0) {
    g2pSkipped = true;
    return;
  }
  const bool qs = (theMode == MODE_QUASISTATIC);
  for (LadrunoMPMPoint *e : eles)
    if (e != 0) e->g2p(eng.P, eng.par.flip, qs);
  g2pRan = true;
  nG2P++;
}

void
LadrunoMPMPattern::mpmElementGone(int p)
{
  if (p < 0 || p >= (int)eles.size()) return;
  eles[p] = 0;
  if (p < (int)matOwnedByEle.size() && matOwnedByEle[p]) {
    if (p < (int)mats.size()) mats[p] = 0;      // deleted with its element
    matOwnedByEle[p] = 0;
  }
}

double
LadrunoMPMPattern::mpmGhostSpring(int p)
{
  if (ghostGamma < 0.0 || p < 0 || p >= (int)mats.size() || mats[p] == 0) return 0.0;
  const double h = eng.grid.h;
  // M_p0: the largest normal diagonal of the INITIAL tangent (slice 1's normalModulus)
  const Matrix &D = mats[p]->getInitialTangent();
  const int nNormal = (D.noRows() == 6) ? 3 : 2;
  double M = 0.0;
  for (int i = 0; i < nNormal; i++) M = std::max(M, D(i, i));
  return ghostGamma * M * eng.P.V[p] / (h * h);
}

// ---------------------------------------------------------------------------
//  Settings and queries
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::setMode(Mode m, std::string &err)
{
  if (m == theMode) return 0;
  if (stepped || !eles.empty()) {
    err = "the mode is fixed at the first step of an epoch; `mpm reset` first";
    return -1;
  }
  theMode = m;
  return 0;
}

double
LadrunoMPMPattern::modeAMassCutoff(void) const
{
  return masscutUser ? eng.par.masscut : 1.0e-6;   // freeze decision 2: mode A 1e-6, mode F 1e-10
}

int
LadrunoMPMPattern::setTagBases(int nb, int eb, std::string &err)
{
  if (!eles.empty() || !liveTag.empty()) {
    err = "tag bases are fixed once grid nodes or particle elements exist";
    return -1;
  }
  if (nb < 1 || eb < 1) {
    err = "tag bases must be >= 1";
    return -1;
  }
  nodeBase = nextNodeTag = nb;
  eleBase = eb;
  basesSet = true;
  return 0;
}

bool
LadrunoMPMPattern::isGridNode(int tag) const
{
  for (const auto &kv : liveTag)
    if (kv.second == tag) return true;
  return false;
}

long
LadrunoMPMPattern::infoCounter(const std::string &w, bool &ok) const
{
  ok = true;
  if (w == "g2p") return nG2P;
  if (w == "acceptedSteps") return nAccepted;
  if (w == "retries") return nRetries;
  if (w == "domainChanges") return nDomainChanges;
  if (w == "connectivityChanges") return nConnChanges;
  if (w == "gridNodes") return (long)liveTag.size();
  if (w == "mode") return (long)theMode;
  ok = false;
  return 0;
}

void
LadrunoMPMPattern::gridNodeTags(std::vector<double> &out) const
{
  out.clear();
  for (const auto &kv : liveTag) out.push_back((double)kv.second);
  std::sort(out.begin(), out.end());
}

int
LadrunoMPMPattern::addGridVel(const std::string &name, const double *lo, const double *hi,
                              const double *v, unsigned freeMask, double start, int track,
                              std::string &err)
{
  const int nd = eng.ndm();
  for (const GridVel &g : gridVels)
    if (g.name == name) { err = "a gridvel named '" + name + "' exists"; return -1; }
  if (track >= nd) { err = "track axis out of range"; return -1; }
  GridVel g;
  g.name = name;
  for (int d = 0; d < 3; d++) {
    g.lo[d] = (d < nd) ? lo[d] : 0.0;
    g.hi[d] = (d < nd) ? hi[d] : 0.0;
    g.v[d] = (d < nd) ? v[d] : 0.0;
  }
  g.freeMask = freeMask;
  g.start = start;
  g.track = track;
  g.force.assign(nd, 0.0);
  gridVels.push_back(g);
  return 0;
}

int
LadrunoMPMPattern::condForce(const std::string &name, std::vector<double> &out, std::string &err)
{
  for (const GridVel &g : gridVels)
    if (g.name == name) { out = g.force; return 0; }
  err = "no gridvel named '" + name + "'";
  return -1;
}

int
LadrunoMPMPattern::addPointForce(int body, const double *lo, const double *hi, const double *f,
                                 std::string &err)
{
  if (stepped) { err = "pforce is an initial condition: before the first step (`mpm reset` first)"; return -1; }
  const int nd = eng.ndm(), N = eng.P.size();
  if ((int)pforce.size() != nd * N) pforce.resize(nd * N, 0.0);
  int count = 0;
  for (int p = 0; p < N; p++) {
    if (eng.P.body[p] != body) continue;
    bool in = true;
    for (int d = 0; d < nd; d++) in = in && eng.P.x[nd * p + d] >= lo[d] && eng.P.x[nd * p + d] <= hi[d];
    if (!in) continue;
    for (int d = 0; d < nd; d++) pforce[nd * p + d] += f[d];
    count++;
  }
  return count;
}

// ---------------------------------------------------------------------------
//  Admission (interface §4)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::admitAnalysis(const LadrunoMPMAnalysisRefs &A, std::string &err)
{
  Domain *dom = this->getDomain();
  if (!A.bridged) {
    err = "this engine has no analysis bridge for mode A";
    return -1;
  }
  if (theMode == MODE_IMPLICIT) {
    if (A.ta == 0) {
      err = "mode implicit needs a transient analysis: build `analysis Transient` (Newmark, HHT or"
            " GeneralizedAlpha) first";
      return -1;
    }
    if (dynamic_cast<VariableTimeStepDirectIntegrationAnalysis *>(A.ta) != 0) {
      err = "VariableTransient is refused in mode A: the driver owns the step size (ADR 232 §5.1)";
      return -1;
    }
    if (dynamic_cast<PFEMAnalysis *>(A.ta) != 0) {
      err = "a PFEM analysis is refused in mode A";
      return -1;
    }
    TransientIntegrator *ti = A.ta->getIntegrator();
    if (ti == 0) { err = "the transient analysis has no integrator"; return -1; }
    const int tag = ti->getClassTag();
    if (tag == INTEGRATOR_TAGS_LadrunoHHT || tag == INTEGRATOR_TAGS_LadrunoGeneralizedAlpha) {
      err = std::string(tag == INTEGRATOR_TAGS_LadrunoHHT ? "LadrunoHHT (33013)" : "LadrunoGeneralizedAlpha (33014)") +
            " stays refused in mode A until a test admits it (ADR 232 §5.3)";
      return -1;
    }
    if (tag != INTEGRATOR_TAGS_Newmark && tag != INTEGRATOR_TAGS_HHT &&
        tag != INTEGRATOR_TAGS_GeneralizedAlpha) {
      err = "transient integrator class tag " + std::to_string(tag) + " is refused in mode A: only"
            " Newmark, HHT and GeneralizedAlpha re-seed from the committed nodal state in"
            " domainChanged() (ADR 232 §5.3, W1)";
      return -1;
    }
  } else {
    if (A.sa == 0) {
      err = "mode quasiStatic needs a static analysis: build `analysis Static` with `integrator"
            " LoadControl` first";
      return -1;
    }
    const int tag = (A.si != 0) ? A.si->getClassTag() : -1;
    if (tag == INTEGRATOR_TAGS_DisplacementControl) {
      err = "DisplacementControl is refused in mode A: never on a grid node (ADR 232 §5.3), and phase 1"
            " has no FE node to drive (amendment 232f); use `mpm gridvel`";
      return -1;
    }
    if (tag == INTEGRATOR_TAGS_ArcLength || tag == INTEGRATOR_TAGS_ArcLength1 ||
        tag == INTEGRATOR_TAGS_LadrunoArcLength) {
      err = "arc-length is out of scope in mode A (ADR 232 §5.3)";
      return -1;
    }
    if (tag != INTEGRATOR_TAGS_LoadControl) {
      err = "static integrator class tag " + std::to_string(tag) + " is refused in mode A: phase 1"
            " admits LoadControl only (ADR 232 §5.3 amendment 232f)";
      return -1;
    }
  }
  // a prescribed (non-zero) velocity needs a handler that enforces it
  bool nonHomog = false;
  for (const GridVel &g : gridVels)
    for (int d = 0; d < eng.ndm(); d++)
      if (!(g.freeMask & (1u << d)) && g.v[d] != 0.0) nonHomog = true;
  if (nonHomog && A.handler != 0 && A.handler->getClassTag() == HANDLER_TAG_PlainHandler) {
    err = "a prescribed-velocity gridvel needs `constraints Transformation`, `Penalty` or `Lagrange`:"
          " the Plain handler assumes every constraint homogeneous (PlainHandler.cpp:92)";
    return -1;
  }
  for (LadrunoMPMPoint *e : eles)
    if (e != 0 && e->refusedRayleigh()) {
      err = "betaK0 / betaKc Rayleigh damping is refused on particle elements: they name a stiffness"
            " of an earlier support (ADR 232 §6 amendment 232f); use `mpm rayleigh <alphaM> <betaK>`";
      return -1;
    }
  for (int f = 0; f < 6; f++)
    if (eng.fc.wallMu[f] >= 0.0) {
      err = "`mpm wall` is mode F only (multi-body contact is phase 7)";
      return -1;
    }
  if (eng.par.damping != 0.0) {
    err = "`mpm damping local` is mode F only";
    return -1;
  }
  LoadPatternIter &pats = dom->getLoadPatterns();
  LoadPattern *pat;
  while ((pat = pats()) != 0)
    if (pat->getClassTag() == PATTERN_TAG_UniformExcitation ||
        pat->getClassTag() == PATTERN_TAG_MultiSupportPattern) {
      err = "ground-motion patterns are refused in mode A in phase 1";
      return -1;
    }
  std::string which;
  if (this->domainHasFreeDOF(which)) {
    err = "the Domain has a free FE DOF (" + which + "); phase 1 has no FE-MPM seam";
    return -1;
  }
  return 0;
}

// ---------------------------------------------------------------------------
//  Lifecycle of the grid nodes, elements and SPs
// ---------------------------------------------------------------------------
void
LadrunoMPMPattern::removeHostSPs(void)
{
  for (int t : hostSPTags) {
    SP_Constraint *sp = this->removeSP_Constraint(t);
    delete sp;
  }
  hostSPTags.clear();
  lastSpec.clear();
}

void
LadrunoMPMPattern::retireNode(long ext)
{
  auto it = liveTag.find(ext);
  if (it == liveTag.end()) return;
  Domain *dom = this->getDomain();
  if (dom != 0) {
    Node *n = dom->removeNode(it->second);
    delete n;
  }
  liveTag.erase(it);
}

// Mode A torn down: SPs, elements (materials back to the host), grid nodes.
// `mpm reset` and the destructor (`remove loadPattern`). Under `wipe` the
// Domain's containers are already empty: every lookup misses and nothing is
// touched twice.
void
LadrunoMPMPattern::teardownModeA(void)
{
  Domain *dom = this->getDomain();
  removeHostSPs();
  const int nOld = (int)eles.size();
  for (int p = 0; p < (int)eles.size(); p++) {
    LadrunoMPMPoint *e = eles[p];
    if (e == 0) continue;
    if (dom != 0 && dom->getElement(e->getTag()) == e) dom->removeElement(e->getTag());
    if (p < (int)mats.size()) mats[p] = e->releaseMaterial();
    matOwnedByEle[p] = 0;
    e->setLink(std::shared_ptr<LadrunoMPMHostLink>());
    delete e;
  }
  eles.clear();
  std::vector<long> exts;
  for (const auto &kv : liveTag) exts.push_back(kv.first);
  for (long x : exts) retireNode(x);
  if (nOld > 0) eleBase += nOld;          // a removed element's tag is not reused either
  if (dom != 0) dom->domainChange();
}

int
LadrunoMPMPattern::createElements(std::string &err)
{
  Domain *dom = this->getDomain();
  const int N = eng.P.size();
  if (!basesSet) {
    int maxE = 0, maxN = 0;
    ElementIter &ei = dom->getElements();
    Element *e;
    while ((e = ei()) != 0) maxE = std::max(maxE, e->getTag());
    NodeIter &ni = dom->getNodes();
    Node *n;
    while ((n = ni()) != 0) maxN = std::max(maxN, n->getTag());
    eleBase = maxE + 1;
    nodeBase = nextNodeTag = std::max(nextNodeTag, maxN + 1);
    basesSet = true;
  }
  if ((long)eleBase + N > INT_MAX) { err = "particle element tags would overflow int"; return -1; }
  eles.assign(N, (LadrunoMPMPoint *)0);
  matOwnedByEle.assign(N, 0);
  for (int p = 0; p < N; p++) {
    if (mats[p] == 0) { err = "internal: particle " + std::to_string(p) + " has no material"; return -1; }
    LadrunoMPMPoint *e = new LadrunoMPMPoint(eleBase + p, p, this->getTag(), eng.ndm(), mats[p]);
    matOwnedByEle[p] = 1;              // the element owns it now; mats[p] stays a view
    e->setLink(link);
    eles[p] = e;
  }
  return 0;
}

// After a database restore: the elements are in the Domain, the host is new.
int
LadrunoMPMPattern::relinkElements(std::string &err)
{
  Domain *dom = this->getDomain();
  const int N = eng.P.size();
  eles.assign(N, (LadrunoMPMPoint *)0);
  matOwnedByEle.assign(N, 0);
  mats.assign(N, (NDMaterial *)0);
  for (int p = 0; p < N; p++) {
    LadrunoMPMPoint *e = dynamic_cast<LadrunoMPMPoint *>(dom->getElement(eleBase + p));
    if (e == 0 || e->particle() != p || e->hostTag() != this->getTag()) {
      err = "restore: particle element " + std::to_string(eleBase + p) + " is missing or not this"
            " host's (the Domain was changed after the restore)";
      return -1;
    }
    e->setLink(link);
    eles[p] = e;
    mats[p] = e->material();
    matOwnedByEle[p] = 1;
  }
  elementsRestored = false;
  return 0;
}

int
LadrunoMPMPattern::prepareStep(double, std::string &err)
{
  Domain *dom = this->getDomain();
  ladruno_mpm::Particles &P = eng.P;
  const ladruno_mpm::Grid &G = eng.grid;
  const int N = P.size(), nd = eng.ndm();
  const double t = dom->getCurrentTime();

  // ---- supports at the committed positions --------------------------------------
  std::vector<long> supStart(1, 0), node;
  std::vector<double> Sv, dSv;
  for (int p = 0; p < N; p++) {
    ladruno_mpm::particleSupport(G, eng.par.basis, &P.x[nd * p], &P.lp[nd * p], node, Sv, dSv);
    supStart.push_back((long)node.size());
  }
  std::map<long, double> mass;        // ordered: deterministic node creation
  double mmax = 0.0;
  for (int p = 0; p < N; p++) {
    mmax = std::max(mmax, P.m[p]);
    for (long k = supStart[p]; k < supStart[p + 1]; k++) mass[node[k]] += Sv[k] * P.m[p];
  }
  const double mcut = modeAMassCutoff() * mmax;
  const bool ghost = ghostGamma >= 0.0;
  std::vector<char> keep(node.size(), 0), inact(node.size(), 0);
  std::set<long> want;
  for (int p = 0; p < N; p++) {
    int nk = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      const bool act = mass[node[k]] > mcut;
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
  }

  // ---- the grid conditions, in lattice terms ------------------------------------------
  std::map<std::pair<long, int>, std::pair<double, bool>> spec;   // (ext,dof) -> (valueR, constant)
  std::vector<std::vector<long>> gvExt(gridVels.size());
  const double tol = 1e-9 * G.h;
  for (long x : want) {
    int ijk[3];
    G.extIJK(x, ijk);
    for (int f = 0; f < 2 * nd; f++) {
      const int d = f / 2;
      const bool on = (f % 2 == 0) ? ijk[d] <= 0 : ijk[d] >= G.n[d];
      if (!on || eng.fc.fixMask[f] == 0u) continue;
      for (int a = 0; a < nd; a++)
        if (eng.fc.fixMask[f] & (1u << a)) spec[std::make_pair(x, a)] = std::make_pair(0.0, true);
    }
  }
  for (size_t gi = 0; gi < gridVels.size(); gi++) {
    const GridVel &gv = gridVels[gi];
    if (t < gv.start - 1e-12 * std::max(1.0, std::fabs(gv.start))) continue;   // absent before start
    double lo[3] = {gv.lo[0], gv.lo[1], gv.lo[2]};
    if (gv.track >= 0) lo[gv.track] += gv.v[gv.track] * std::max(t - gv.start, 0.0);
    for (long x : want) {
      int ijk[3];
      G.extIJK(x, ijk);
      bool in = true;
      for (int d = 0; d < nd; d++) {
        const double c = G.x0[d] + ijk[d] * G.h;
        in = in && c >= lo[d] - tol && c <= gv.hi[d] + tol;
      }
      if (!in) continue;
      gvExt[gi].push_back(x);
      for (int a = 0; a < nd; a++) {
        if (gv.freeMask & (1u << a)) continue;
        auto key = std::make_pair(x, a);
        auto it = spec.find(key);
        if (it != spec.end() && it->second.first != gv.v[a]) {
          err = "gridvel '" + gv.name + "' and another condition prescribe different values on one"
                " node component (a malformed case)";
          return -1;
        }
        spec[key] = std::make_pair(gv.v[a], gv.v[a] == 0.0);
      }
    }
  }
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
  const bool specChanged = nodeSetChanged || !(newSpec == lastSpec);
  if (specChanged) removeHostSPs();
  if (nodeSetChanged) {
    std::vector<long> drop;
    for (const auto &kv : liveTag)
      if (reseedMode == 1 || want.count(kv.first) == 0) drop.push_back(kv.first);
    for (long x : drop) retireNode(x);
    if (!basesSet) {
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
    }
    for (long x : want) {
      if (liveTag.count(x)) continue;
      if (nextNodeTag == INT_MAX) { err = "grid node tags would overflow int (W5: never reused)"; return -1; }
      int ijk[3];
      G.extIJK(x, ijk);
      const double c0 = G.x0[0] + ijk[0] * G.h, c1 = G.x0[1] + ijk[1] * G.h;
      Node *n = (nd == 2) ? new Node(nextNodeTag, 2, c0, c1)
                          : new Node(nextNodeTag, 3, c0, c1, G.x0[2] + ijk[2] * G.h);
      if (!dom->addNode(n)) {
        delete n;
        err = "could not add grid node " + std::to_string(nextNodeTag) + " (a tag in use? see `mpm tagBase`)";
        return -1;
      }
      liveTag[x] = nextNodeTag++;
    }
  }

  // ---- P2G onto the Domain nodes (interface §3 step 3) -------------------------------
  const bool qs = (theMode == MODE_QUASISTATIC);
  std::map<long, std::vector<double>> mv;     // m, m v (nd), m a (nd)
  for (long x : want) mv[x].assign(1 + 2 * nd, 0.0);
  for (int p = 0; p < N; p++)
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (!keep[k]) continue;
      std::vector<double> &r = mv[node[k]];
      const double w = Sv[k] * P.m[p];
      r[0] += w;
      if (!qs)
        for (int a = 0; a < nd; a++) {
          r[1 + a] += w * P.v[nd * p + a];
          r[1 + nd + a] += w * P.a[nd * p + a];
        }
    }
  Vector zero(nd), vel(nd), acc(nd);
  for (auto &kv : mv) {
    Node *n = dom->getNode(liveTag[kv.first]);
    if (n == 0) { err = "internal: grid node missing"; return -1; }
    const std::vector<double> &r = kv.second;
    for (int a = 0; a < nd; a++) {
      vel(a) = (r[0] > 0.0) ? r[1 + a] / r[0] : 0.0;
      acc(a) = (r[0] > 0.0) ? r[1 + nd + a] / r[0] : 0.0;
    }
    n->setTrialDisp(zero);
    n->setTrialVel(vel);
    n->setTrialAccel(acc);
    n->commitState();
  }

  // ---- elements: created, or re-pointed in place ------------------------------------------
  if (eles.empty() && createElements(err) < 0) return -1;
  bool connChanged = false;
  const bool fresh = !eles.empty() && eles[0] != 0 && eles[0]->getDomain() == 0;
  for (int p = 0; p < N; p++) {
    int nk = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) nk += keep[k];
    ID tags(nk);
    std::vector<double> S(nk), dS(nd * nk), vH(nd * nk), aH(nd * nk);
    std::vector<char> ina(nk);
    int j = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (!keep[k]) continue;
      tags(j) = liveTag[node[k]];
      S[j] = Sv[k];
      ina[j] = inact[k];
      const std::vector<double> &r = mv[node[k]];
      for (int a = 0; a < nd; a++) {
        dS[nd * j + a] = dSv[nd * k + a];
        vH[nd * j + a] = (r[0] > 0.0) ? r[1 + a] / r[0] : 0.0;
        aH[nd * j + a] = (r[0] > 0.0) ? r[1 + nd + a] / r[0] : 0.0;
      }
      j++;
    }
    LadrunoMPMPoint *e = eles[p];
    if (e == 0) { err = "internal: particle element " + std::to_string(p) + " is gone"; return -1; }
    if (!e->connectivityEquals(tags)) connChanged = true;
    if (e->repoint(tags, S, dS, ina, vH, aH) < 0) { err = "re-pointing failed"; return -1; }
  }
  if (fresh) {
    for (int p = 0; p < N; p++) {
      eles[p]->setRayleigh(rayAlphaM, rayBetaK);
      if (!dom->addElement(eles[p])) {
        err = "could not add particle element " + std::to_string(eles[p]->getTag()) +
              " (a tag in use? see `mpm tagBase`)";
        return -1;
      }
    }
    connChanged = true;
  }
  if (connChanged) nConnChanges++;

  // ---- the grid conditions as this pattern's SPs ----------------------------------------
  if (specChanged) {
    for (const SpecRow &r : newSpec) {
      SP_Constraint *sp = new SP_Constraint(liveTag[r.ext], r.dof, r.value, r.constant);
      if (!this->addSP_Constraint(sp)) {
        delete sp;
        err = "could not add a grid-condition SP";
        return -1;
      }
      hostSPTags.push_back(sp->getTag());
    }
    lastSpec = newSpec;
  }
  for (size_t gi = 0; gi < gridVels.size(); gi++) {
    gridVels[gi].nodes.clear();
    for (long x : gvExt[gi]) gridVels[gi].nodes.push_back(liveTag[x]);
  }

  // ---- Domain::domainChange() from C++, unconditionally (W1, W4) -------------------------
  if (!(reseedMode == 2 && !nodeSetChanged && !connChanged && !specChanged)) {
    dom->domainChange();
    nDomainChanges++;
  }
  return 0;
}

void
LadrunoMPMPattern::updateCondForces(void)
{
  const int nd = eng.ndm();
  for (GridVel &gv : gridVels) {
    gv.force.assign(nd, 0.0);
    if (gv.nodes.empty()) continue;
    std::set<int> in(gv.nodes.begin(), gv.nodes.end());
    for (LadrunoMPMPoint *e : eles) {
      if (e == 0) continue;
      const ID &c = e->getExternalNodes();
      for (int k = 0; k < c.Size(); k++)
        if (in.count(c(k))) {
          const double *f = e->internalForceAt(k);
          for (int a = 0; a < nd; a++) gv.force[a] += f[a];
        }
    }
  }
}

// ---------------------------------------------------------------------------
//  `mpm analyze` in mode A (interface §3)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::analyzeImplicit(int nSteps, double dt, const LadrunoMPMAnalysisRefs &A)
{
  Domain *dom = this->getDomain();
  if (dom == 0) { opserr << "ERROR mpm analyze -- the MPM host is not in a Domain\n"; return -1; }
  if (!(dt > 0.0) || nSteps < 0) { opserr << "ERROR mpm analyze -- need nSteps >= 0 and dt > 0\n"; return -1; }
  const int N = eng.P.size(), nd = eng.ndm();
  if (N == 0) { opserr << "ERROR mpm analyze -- no particles\n"; return -1; }
  if (aborted) {
    opserr << "ERROR mpm analyze -- an earlier run was aborted by a material refusal at commit"
              " (ADR 232 W17); `mpm reset` first.\n";
    return -1;
  }
  std::string err;
  if (elementsRestored && relinkElements(err) < 0) {
    opserr << "ERROR mpm analyze -- " << err.c_str() << "\n";
    return -1;
  }
  if ((int)mats.size() != N) {
    opserr << "ERROR mpm analyze -- internal: " << (int)mats.size() << " particle materials for "
           << N << " particles; refused\n";
    return -1;
  }
  if (admitAnalysis(A, err) < 0) {
    opserr << "ERROR mpm analyze -- " << err.c_str() << "\n";
    return -1;
  }
  if ((int)pforce.size() != nd * N) pforce.resize(nd * N, 0.0);
  if (!stepped) {
    snapshot = eng.P;
    haveSnapshot = true;
  }
  dtLog.clear();

  LoadControl *lc = 0;
  double ud = 0.0, umin = 0.0, umax = 0.0;
  if (theMode == MODE_QUASISTATIC) {
    lc = static_cast<LoadControl *>(A.si);           // class tag checked by admitAnalysis
    lc->ladrunoGetIncrement(ud, umin, umax);
    if (ud != dt || umin != dt || umax != dt)
      opserr << "mpm analyze (quasiStatic): the LoadControl increment is dt = " << dt
             << " for this call (the user's " << ud << " [" << umin << ", " << umax
             << "] is put back afterwards; interface §3, freeze decision 5)\n";
  }

  int rc = 0;
  if (writer != 0 && writer->isConfigured() &&
      (!writer->stageOpen() || writer->stageEpoch() != epoch)) {
    if (openOutputStage(err) < 0 || writeOutputNow(err) < 0) {
      opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
      rc = -1;
    } else {
      writer->flush();
    }
  }

  const double tStart = dom->getCurrentTime();
  const double tEnd = tStart + nSteps * dt;
  double cur = dt;
  int done = 0, cuts = 0;
  bool wasCut = false;
  while (rc == 0) {
    const double t = dom->getCurrentTime();
    double h;
    if (!wasCut) {
      if (done >= nSteps) break;
      h = dt;
    } else {
      const double remaining = tEnd - t;
      if (remaining <= 1e-9 * cur) break;
      h = (remaining < cur * (1.0 - 1e-9)) ? remaining : cur;
    }
    int attempt = 0;
    for (;;) {
      attempt++;
      if (prepareStep(h, err) < 0) {
        opserr << "ERROR mpm analyze -- " << err.c_str() << "; step refused, nothing changed\n";
        rc = -1;
        break;
      }
      stepT0 = dom->getCurrentTime();
      failNow = (failStep > 0 && failStep == epochSteps + 1 && attempt <= failAttempts);
      if (failNow)
        opserr << "mpm debug failAttempt: attempt " << attempt << " of step " << failStep
               << " is made to fail before commit\n";
      armedCount = N;
      committedCount = 0;
      g2pRan = g2pSkipped = false;
      inStep = true;
      int res;
      if (theMode == MODE_IMPLICIT) {
        res = A.ta->analyzeStep(h);
      } else {
        lc->ladrunoSetIncrement(h, h, h);
        res = A.sa->analyze(1, false);
      }
      inStep = false;
      failNow = false;
      if (res >= 0 && g2pRan) break;                                // accepted
      if (res >= 0 || committedCount > 0) {
        aborted = true;
        if (g2pSkipped)
          opserr << "ERROR mpm analyze -- a particle material REFUSED this commit (ADR 232 W17): the"
                    " run is aborted with the particles unmoved, until `mpm reset`\n";
        else
          opserr << "ERROR mpm analyze -- the step failed after its commit began (" << res
                 << "); the particle state is not consistent: the run is aborted until `mpm reset`\n";
        rc = -1;
        break;
      }
      // failed before any commit: nothing to undo — cut and retry
      nRetries++;
      if (++cuts > maxCutCount) {
        opserr << "ERROR mpm analyze -- the step at t = " << stepT0 << " failed after " << maxCutCount
               << " cut(s); the particles keep the last accepted step\n";
        rc = -1;
        break;
      }
      cur /= cutFactor;
      wasCut = true;
      const double remaining = tEnd - dom->getCurrentTime();
      h = (remaining < cur * (1.0 - 1e-9)) ? remaining : cur;
      opserr << "mpm analyze -- the step failed before commit at t = " << stepT0
             << "; retried with dt = " << h << " (cut " << cuts << " of " << maxCutCount << ")\n";
    }
    if (rc != 0) break;
    stepped = true;
    done++;
    nAccepted++;
    epochSteps++;
    dtLog.push_back(h);
    updateCondForces();
    if (writeAfterStep(err) < 0) {
      opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
      rc = -1;
    }
  }
  if (writer != 0 && writer->stageOpen() && stepsSinceOutput > 0) {
    if (writeOutputNow(err) < 0) {
      opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
      rc = -1;
    }
  }
  if (writer != 0) writer->flush();
  if (lc != 0) lc->ladrunoSetIncrement(ud, umin, umax);
  dom->flushRecorders();                                          // W16
  return rc;
}

int
LadrunoMPMPattern::writeAfterStep(std::string &err)
{
  if (writer == 0 || !writer->stageOpen()) return 0;
  Domain *dom = this->getDomain();
  stepsSinceOutput++;
  const ladruno_mpm_io::WriterConfig &cfg = writer->config();
  const double t = dom->getCurrentTime();
  bool due;
  if (cfg.dtOut > 0.0) {
    if (nextOutputTime - t > cfg.dtOut * (1.0 + 1e-12)) nextOutputTime = t;
    due = t >= nextOutputTime - 1e-12 * std::max(1.0, std::fabs(nextOutputTime));
  } else {
    due = stepsSinceOutput >= cfg.every;
  }
  if (due && writeOutputNow(err) < 0) return -1;
  writer->flushIfDue();
  return 0;
}

// ---------------------------------------------------------------------------
//  Restart (interface §6): the host's state as one Vector of doubles
// ---------------------------------------------------------------------------
namespace {
struct Packer {
  std::vector<double> d;
  void put(double x) { d.push_back(x); }
  void puts(const std::string &s) { put((double)s.size()); for (char c : s) put((double)(unsigned char)c); }
  template <class T> void putv(const std::vector<T> &v) { put((double)v.size()); for (const T &x : v) put((double)x); }
};
struct Unpacker {
  const Vector &v;
  int i = 0;
  bool bad = false;
  explicit Unpacker(const Vector &x) : v(x) {}
  double get() { if (i >= v.Size()) { bad = true; return 0.0; } return v(i++); }
  std::string gets() { const int n = (int)get(); std::string s; for (int k = 0; k < n && !bad; k++) s.push_back((char)(int)get()); return s; }
  template <class T> void getv(std::vector<T> &out) { const int n = (int)get(); out.resize(n > 0 ? n : 0); for (int k = 0; k < n && !bad; k++) out[k] = (T)get(); }
};

void putParticles(Packer &k, const ladruno_mpm::Particles &P)
{
  k.put(P.ndm); k.put(P.nv);
  k.putv(P.id); k.putv(P.body); k.putv(P.mat); k.putv(P.flags);
  k.putv(P.m); k.putv(P.V0); k.putv(P.lp); k.putv(P.x); k.putv(P.v); k.putv(P.V);
  k.putv(P.F); k.putv(P.R); k.putv(P.eps); k.putv(P.sigHat); k.putv(P.sig); k.putv(P.sig0); k.putv(P.a);
}

void getParticles(Unpacker &u, ladruno_mpm::Particles &P)
{
  P.ndm = (int)u.get(); P.nv = (int)u.get();
  u.getv(P.id); u.getv(P.body); u.getv(P.mat); u.getv(P.flags);
  u.getv(P.m); u.getv(P.V0); u.getv(P.lp); u.getv(P.x); u.getv(P.v); u.getv(P.V);
  u.getv(P.F); u.getv(P.R); u.getv(P.eps); u.getv(P.sigHat); u.getv(P.sig); u.getv(P.sig0); u.getv(P.a);
}
}  // namespace

const int MPM_PACK_VERSION = 1;

int
LadrunoMPMPattern::packState(Vector &out)
{
  Packer k;
  const ladruno_mpm::Grid &G = eng.grid;
  k.put(MPM_PACK_VERSION);
  k.put(G.ndm);
  for (int d = 0; d < 3; d++) k.put(G.x0[d]);
  k.put(G.h);
  for (int d = 0; d < 3; d++) k.put(G.n[d]);
  k.put(eng.par.basis); k.put(eng.par.flip); k.put(eng.par.damping); k.put(eng.par.masscut);
  k.put(eng.par.safety);
  for (int d = 0; d < 3; d++) k.put(eng.par.gravity[d]);
  for (int f = 0; f < 6; f++) { k.put(eng.fc.fixMask[f]); k.put(eng.fc.wallMu[f]); }
  k.put(theMode); k.put(masscutUser); k.put(ghostGamma); k.put(cutFactor); k.put(maxCutCount);
  k.put(rayAlphaM); k.put(rayBetaK);
  k.put(stepped); k.put(haveSnapshot); k.put(epoch); k.put(aborted); k.put((double)nextLatticeId);
  k.put(thickness); k.put(thicknessFixed); k.put(lambda);
  k.put(nodeBase); k.put(nextNodeTag); k.put(eleBase); k.put(basesSet); k.put(eles.empty() ? 0 : 1);
  k.put(nG2P); k.put(nAccepted); k.put(nRetries); k.put(nDomainChanges); k.put(nConnChanges);
  k.put(epochSteps);
  putParticles(k, eng.P);
  putParticles(k, snapshot);
  k.putv(elastic);
  k.putv(pforce);
  k.put((double)liveTag.size());
  for (const auto &kv : liveTag) { k.put((double)kv.first); k.put(kv.second); }
  k.put((double)gridVels.size());
  for (const GridVel &g : gridVels) {
    k.puts(g.name);
    for (int d = 0; d < 3; d++) { k.put(g.lo[d]); k.put(g.hi[d]); k.put(g.v[d]); }
    k.put(g.freeMask); k.put(g.start); k.put(g.track);
    k.putv(g.nodes); k.putv(g.force);
  }
  k.put((double)bodyNames.size());
  for (const auto &kv : bodyNames) { k.put(kv.first); k.puts(kv.second); }
  k.puts(unitL); k.puts(unitF); k.puts(unitT); k.puts(unitSystem); k.put(unitsFromFile);
  out.resize((int)k.d.size());
  for (int i = 0; i < (int)k.d.size(); i++) out(i) = k.d[i];
  return 0;
}

int
LadrunoMPMPattern::unpackState(const Vector &in)
{
  Unpacker u(in);
  if ((int)u.get() != MPM_PACK_VERSION) return -1;
  int ndm = (int)u.get();
  double x0[3], h;
  int n[3];
  for (int d = 0; d < 3; d++) x0[d] = u.get();
  h = u.get();
  for (int d = 0; d < 3; d++) n[d] = (int)u.get();
  std::string err;
  eng.clearParticles();
  if (!eng.setGrid(ndm, x0, h, n, err)) return -1;
  eng.par.basis = (ladruno_mpm::Basis)(int)u.get();
  eng.par.flip = u.get(); eng.par.damping = u.get(); eng.par.masscut = u.get(); eng.par.safety = u.get();
  for (int d = 0; d < 3; d++) eng.par.gravity[d] = u.get();
  for (int f = 0; f < 6; f++) { eng.fc.fixMask[f] = (unsigned)u.get(); eng.fc.wallMu[f] = u.get(); }
  theMode = (Mode)(int)u.get(); masscutUser = u.get() != 0.0; ghostGamma = u.get();
  cutFactor = u.get(); maxCutCount = (int)u.get(); rayAlphaM = u.get(); rayBetaK = u.get();
  stepped = u.get() != 0.0; haveSnapshot = u.get() != 0.0; epoch = (int)u.get();
  aborted = u.get() != 0.0; nextLatticeId = (int64_t)u.get();
  thickness = u.get(); thicknessFixed = u.get() != 0.0; lambda = u.get();
  nodeBase = (int)u.get(); nextNodeTag = (int)u.get(); eleBase = (int)u.get(); basesSet = u.get() != 0.0;
  elementsRestored = u.get() != 0.0;
  nG2P = (long)u.get(); nAccepted = (long)u.get(); nRetries = (long)u.get();
  nDomainChanges = (long)u.get(); nConnChanges = (long)u.get(); epochSteps = (long)u.get();
  ladruno_mpm::Particles P;
  getParticles(u, P);
  eng.setParticles(P);
  getParticles(u, snapshot);
  u.getv(elastic);
  u.getv(pforce);
  liveTag.clear();
  const int nl = (int)u.get();
  for (int i = 0; i < nl && !u.bad; i++) { const long x = (long)u.get(); liveTag[x] = (int)u.get(); }
  gridVels.clear();
  const int ng = (int)u.get();
  for (int i = 0; i < ng && !u.bad; i++) {
    GridVel g;
    g.name = u.gets();
    for (int d = 0; d < 3; d++) { g.lo[d] = u.get(); g.hi[d] = u.get(); g.v[d] = u.get(); }
    g.freeMask = (unsigned)u.get(); g.start = u.get(); g.track = (int)u.get();
    u.getv(g.nodes); u.getv(g.force);
    gridVels.push_back(g);
  }
  bodyNames.clear();
  const int nb = (int)u.get();
  for (int i = 0; i < nb && !u.bad; i++) { const int id = (int)u.get(); bodyNames[id] = u.gets(); }
  unitL = u.gets(); unitF = u.gets(); unitT = u.gets(); unitSystem = u.gets();
  unitsFromFile = u.get() != 0.0;
  if (u.bad) return -1;
  // the per-step SPs are rebuilt by the next step; the materials live in the
  // restored elements (relinked at the next `mpm analyze`)
  hostSPTags.clear();
  lastSpec.clear();
  eles.clear();
  matOwnedByEle.clear();
  mats.assign(eng.P.size(), (NDMaterial *)0);
  return 0;
}
