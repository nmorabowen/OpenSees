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
#include <LinearSOE.h>
#include <profiler/ProfilerMacros.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <optional>
#include <set>

// WP-232g (M0b): sequential profiler phases inside one function whose locals outlive the
// phase, so a braced OPS_PROFILE_SCOPE cannot hold them. MPM_PHASE_BEGIN opens a scope
// held in a std::optional, MPM_PHASE_END closes it; every phase ends before the next
// begins, so the profiler's cursor nests them as siblings. Same gate as OPS_PROFILE_SCOPE
// (ProfilerMacros.h: the name is passed only when the profiler is enabled), and nothing
// at all under OPS_PROFILE_DISABLE. An early return closes an open phase through the
// optional's destructor.
#ifndef OPS_PROFILE_DISABLE
#define MPM_PHASE_BEGIN(var, name)                                                       \
  std::optional<::ops_profiler::ScopedTimerOpt> var;                                     \
  var.emplace(::ops_profiler::theProfiler().enabled() ? (name) : nullptr)
#define MPM_PHASE_NEXT(var, name)                                                        \
  do { var.reset();                                                                      \
       var.emplace(::ops_profiler::theProfiler().enabled() ? (name) : nullptr); } while (0)
#define MPM_PHASE_END(var) var.reset()
#else
#define MPM_PHASE_BEGIN(var, name) ((void)0)
#define MPM_PHASE_NEXT(var, name) ((void)0)
#define MPM_PHASE_END(var) ((void)0)
#endif

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
    // WP-232i B11: the particle body force per unit mass, NOT scaled by the load factor
    if ((int)eng.bodyForce.size() >= nd * (p + 1)) bfBuf[a] += mp * eng.bodyForce[nd * p + a];
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
  OPS_PROFILE_SCOPE("mpm.a.g2p");   // WP-232g (M0b): the batched G2P of every particle
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

int
LadrunoMPMPattern::setAntiLocking(int on, std::string &err)
{
  if ((on ? 1 : 0) == eng.par.antiLocking) return 0;
  // C6: only an accepted step fixes it -- elements a debug hook (tangentCheck, affineCheck)
  // created are re-pointed, corners included, at every step
  if (stepped) {
    err = "the anti-locking setting is fixed at the first step of an epoch; `mpm reset` first";
    return -1;
  }
  eng.par.antiLocking = on ? 1 : 0;
  return 0;
}

double
LadrunoMPMPattern::modeAMassCutoff(void) const
{
  // D2.3 amendment (232f, 2026-11-14): mode A 5e-2 with `renormalize` (the study on #1083);
  // mode F keeps 1e-10. `mpm massCutoff` sets both.
  return masscutUser ? eng.par.masscut : 5.0e-2;
}

// `mpm rayleigh <alphaM> <betaK>`: stored for elements made later AND applied to the live
// ones at once (review item 9: it used to be applied only when the elements were created).
void
LadrunoMPMPattern::setRayleighMPM(double am, double bk)
{
  rayAlphaM = am;
  rayBetaK = bk;
  for (LadrunoMPMPoint *e : eles)
    if (e != 0) e->setRayleigh(am, bk);
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
  return liveTags.count(tag) != 0;      // review item 10: O(1), not a scan per Domain node
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
  // WP-232j (K3, amendment 5): the F-bar tangent is unsymmetric even for elastic
  // materials (J10); a symmetric-storage system would solve a different matrix.
  if (eng.par.antiLocking && A.soe != 0) {
    const int t = A.soe->getClassTag();
    // C3: every symmetric-storage SOE class (classTags.h), and the two whose symmetry is a
    // matType option (PARDISO -matType / -sym, MUMPS -matrixType 1|2), read by the bridges
    const bool solverSym = (t == LinSOE_TAGS_PARDISOGenLinSOE || t == LinSOE_TAGS_MumpsSOE ||
                            t == LinSOE_TAGS_MumpsParallelSOE) && A.pardisoMatType > 0;
    const bool sym = t == LinSOE_TAGS_BandSPDLinSOE || t == LinSOE_TAGS_ProfileSPDLinSOE ||
                     t == LinSOE_TAGS_SymSparseLinSOE || t == LinSOE_TAGS_ProfileSPDLinSOEGather ||
                     t == LinSOE_TAGS_DistributedBandSPDLinSOE ||
                     t == LinSOE_TAGS_DistributedProfileSPDLinSOE ||
                     t == LinSOE_TAGS_SProfileSPDLinSOE || solverSym;
    if (sym) {
      err = "`mpm antiLocking fbar` makes the particle tangent unsymmetric in general: the"
            " symmetric-storage system (class tag " + std::to_string(t) +
            (solverSym ? ", a symmetric matrix type" : "") +
            ") is refused; use an unsymmetric one: UmfPack, SparseGEN, FullGeneral,"
            " `system Pardiso` (matType 0) or `system Mumps -matrixType 0`";
      return -1;
    }
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
  liveTags.erase(it->second);
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

  // WP-232g (M0b): the host's preparation, split into phases (profiler scopes only).
  OPS_PROFILE_SCOPE("mpm.a.prepare");
  MPM_PHASE_BEGIN(ph, "mpm.a.support");      // supports, nodal mass, active set
  // ---- supports at the committed positions --------------------------------------
  std::vector<long> supStart(1, 0), node;
  std::vector<double> Sv, dSv;
  // WP-232i (tests/mpm_bspline_interface.md B8): on a B-spline grid the support is
  // folded (the ghosts behind a `fix` face are images, never Domain nodes): Sv, dSv
  // then hold the unsigned weight W and grad W, Sc, dSc the per-component S' and
  // grad S' (nd, nd x nd per entry); oddP[p] marks a particle with a sign -1, which
  // takes the per-component path in its element. GIMP and linear: unchanged.
  const bool bspl = ladruno_mpm::isBSpline(eng.par.basis);
  std::vector<double> Sc, dSc;
  std::vector<char> oddP(N, 0);
  for (int p = 0; p < N; p++) {
    if (bspl) {
      oddP[p] = ladruno_mpm::foldedSupport(G, eng.par.basis, eng.fc, &P.x[nd * p], &P.lp[nd * p],
                                           node, Sv, dSv, Sc, dSc) ? 1 : 0;
    } else {
      ladruno_mpm::particleSupport(G, eng.par.basis, &P.x[nd * p], &P.lp[nd * p], node, Sv, dSv);
    }
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
      // critic A1: a face-layer node of a `fix` face is tested with the mass of the body
      // plus its reflections (x 2^k); GIMP and linear unchanged
      const bool act = ladruno_mpm::imageCutoffMass(G, eng.par.basis, eng.fc, node[k],
                                                    mass[node[k]]) > mcut;
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
    // `mpm emptyNodes renormalize` (WP-232f round 2, proposed amendment to D2.3): the
    // excluded nodes are dropped AND the kept S, ∇S are renormalised to a partition
    // of unity -- S' = S / ΣS, ∇S' = (∇S − S' Σ∇S) / ΣS over the kept nodes -- so an
    // excluded node is absent, not held at zero increment (no anchoring of the surface).
    if (ghostGamma == EMPTY_RENORMALIZE) {
      double sum = 0.0, gsum[3] = {0.0, 0.0, 0.0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++)
        if (keep[k]) {
          sum += Sv[k];
          for (int a = 0; a < nd; a++) gsum[a] += dSv[nd * k + a];
        }
      if (sum > 0.0)
        for (long k = supStart[p]; k < supStart[p + 1]; k++)
          if (keep[k]) {
            Sv[k] /= sum;
            for (int a = 0; a < nd; a++) dSv[nd * k + a] = (dSv[nd * k + a] - Sv[k] * gsum[a]) / sum;
            if (bspl)            // B9: every component scaled alike, by the unsigned sums
              for (int c = 0; c < nd; c++) {
                double &sc = Sc[nd * k + c];
                sc /= sum;
                for (int a = 0; a < nd; a++)
                  dSc[(nd * k + c) * nd + a] = (dSc[(nd * k + c) * nd + a] - sc * gsum[a]) / sum;
              }
          }
    }
  }

  MPM_PHASE_NEXT(ph, "mpm.a.sp");            // WP-232g: the grid conditions' spec
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

  MPM_PHASE_NEXT(ph, "mpm.a.churn");         // WP-232g: grid-node retire / create
  // ---- the node diff -------------------------------------------------------------------
  bool nodeSetChanged = (reseedMode == 1);
  if (!nodeSetChanged) {
    if (want.size() != liveTag.size()) nodeSetChanged = true;
    else
      for (long x : want)
        if (liveTag.find(x) == liveTag.end()) { nodeSetChanged = true; break; }
  }
  const bool specChanged = nodeSetChanged || !(newSpec == lastSpec);
  // From here on the grid state is being changed: a failure is -2 (the driver aborts and
  // `mpm reset` is required -- review item 8), never "nothing changed".
  MPM_PHASE_NEXT(ph, "mpm.a.sp");            // WP-232g
  if (specChanged) removeHostSPs();
  MPM_PHASE_NEXT(ph, "mpm.a.churn");         // WP-232g
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
      if (nextNodeTag == INT_MAX) { err = "grid node tags would overflow int (W5: never reused)"; return -2; }
      int ijk[3];
      G.extIJK(x, ijk);
      const double c0 = G.x0[0] + ijk[0] * G.h, c1 = G.x0[1] + ijk[1] * G.h;
      Node *n = (nd == 2) ? new Node(nextNodeTag, 2, c0, c1)
                          : new Node(nextNodeTag, 3, c0, c1, G.x0[2] + ijk[2] * G.h);
      if (!dom->addNode(n)) {
        delete n;
        err = "could not add grid node " + std::to_string(nextNodeTag) + " (a tag in use? see `mpm tagBase`)";
        return -2;
      }
      liveTags.insert(nextNodeTag);
      liveTag[x] = nextNodeTag++;
    }
  }

  MPM_PHASE_NEXT(ph, "mpm.a.p2g");           // WP-232g
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
          const double wa = bspl ? Sc[nd * k + a] * P.m[p] : w;   // B8: momentum reflected
          r[1 + a] += wa * P.v[nd * p + a];
          r[1 + nd + a] += wa * P.a[nd * p + a];
        }
    }
  // The NODE's committed v, a take the grid conditions (review item 2): a fixed component
  // v = a = 0, a prescribed-velocity component v = v_d (transient) and a = 0. The analysis
  // never writes a constrained DOF's committed v / a back (DOF_Group::setNodeVel skips
  // loc < 0), so Newmark's start, Rayleigh and Penalty / Lagrange would see the raw
  // projection otherwise. The ELEMENT's vH, aH (the FLIP reference) stay the unconditioned
  // projection -- mode F's v_n = p / m: a velocity a condition imposes at the step start is
  // then part of the FLIP increment, so a body driven from rest picks it up (with the
  // conditioned value as reference, FLIP would never see the jump).
  for (auto &kv : mv) {
    std::vector<double> &r = kv.second;
    for (int a = 0; a < nd; a++) {          // r becomes [m, v, a], unconditioned
      r[1 + a] = (r[0] > 0.0) ? r[1 + a] / r[0] : 0.0;
      r[1 + nd + a] = (r[0] > 0.0) ? r[1 + nd + a] / r[0] : 0.0;
    }
  }
  // WP-232i (B6, found by G8): on a B-spline grid a face-layer node lies on the
  // symmetry plane of its `fix` face, so the components that face fixes project
  // to exactly 0 (body + reflection), not to the one-sided value.
  if (bspl)
    for (auto &kv : mv) {
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
  Vector zero(nd), vel(nd), acc(nd);
  for (auto &kv : mv) {
    Node *n = dom->getNode(liveTag[kv.first]);
    if (n == 0) { err = "internal: grid node missing"; return -2; }
    const std::vector<double> &r = kv.second;
    for (int a = 0; a < nd; a++) {
      vel(a) = r[1 + a];
      acc(a) = r[1 + nd + a];
      auto it = spec.find(std::make_pair(kv.first, a));
      if (it != spec.end()) {
        vel(a) = qs ? 0.0 : it->second.first;
        acc(a) = 0.0;
      }
    }
    n->setTrialDisp(zero);
    n->setTrialVel(vel);
    n->setTrialAccel(acc);
    n->commitState();
  }

  MPM_PHASE_NEXT(ph, "mpm.a.repoint");       // WP-232g
  // ---- elements: created, or re-pointed in place ------------------------------------------
  if (eles.empty() && createElements(err) < 0) return -2;
  bool connChanged = false;
  alPend.assign(eng.par.antiLocking ? N : 0, 0);
  for (int p = 0; p < N; p++) {
    int nk = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) nk += keep[k];
    ID tags(nk);
    std::vector<double> S(nk), dS(nd * nk), vH(nd * nk), aH(nd * nk);
    std::vector<double> Spc, dSpc;          // WP-232i B8: per component, odd folds only
    if (bspl && oddP[p]) { Spc.resize(nd * nk); dSpc.resize(nd * nd * nk); }
    std::vector<char> ina(nk);
    int j = 0;
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      if (!keep[k]) continue;
      tags(j) = liveTag[node[k]];
      S[j] = Sv[k];
      ina[j] = inact[k];
      if (!Spc.empty()) {
        for (int c = 0; c < nd; c++) {
          Spc[nd * j + c] = Sc[nd * k + c];
          for (int a = 0; a < nd; a++) dSpc[(nd * j + c) * nd + a] = dSc[(nd * k + c) * nd + a];
        }
      }
      const std::vector<double> &r = mv[node[k]];
      for (int a = 0; a < nd; a++) {
        dS[nd * j + a] = dSv[nd * k + a];
        vH[nd * j + a] = r[1 + a];          // the unconditioned projection (FLIP reference)
        aH[nd * j + a] = r[1 + nd + a];
      }
      j++;
    }
    LadrunoMPMPoint *e = eles[p];
    if (e == 0) {
      err = "particle element " + std::to_string(p) + " is gone (removed by hand?)";
      return -2;
    }
    if (!e->connectivityEquals(tags)) connChanged = true;
    // WP-232j (J11): the slots of the particle's cell corners in its kept support; a
    // corner dropped by the cutoff (renormalize / exclude) or inactive (ghost) makes
    // the particle a fallback (J5).
    int cs[8];
    int nc = 0;
    if (eng.par.antiLocking) {
      long ce[8];
      double gN[8][3];
      nc = ladruno_mpm::cellCorners(G, &P.x[nd * p], ce, gN);
      bool all = true;
      for (int kc = 0; kc < nc; kc++) {
        cs[kc] = -1;
        int jj = 0;
        for (long k = supStart[p]; k < supStart[p + 1]; k++) {
          if (!keep[k]) continue;
          if (node[k] == ce[kc] && !inact[k]) { cs[kc] = jj; break; }
          jj++;
        }
        if (cs[kc] < 0) all = false;
      }
      alPend[p] = all ? 1 : 2;
    }
    if (e->repoint(tags, S, dS, ina, vH, aH, Spc, dSpc, nc ? cs : 0, nc, G.h) < 0) {
      err = "re-pointing failed";
      return -2;
    }
  }
  // review item 8: whether an element is in the Domain is decided per element, so a
  // partial add earlier never strands the rest
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

  MPM_PHASE_NEXT(ph, "mpm.a.sp");            // WP-232g: SP_Constraints rebuilt
  // ---- the grid conditions as this pattern's SPs ----------------------------------------
  if (specChanged) {
    for (const SpecRow &r : newSpec) {
      SP_Constraint *sp = new SP_Constraint(liveTag[r.ext], r.dof, r.value, r.constant);
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

  MPM_PHASE_NEXT(ph, "mpm.a.domainChange");  // WP-232g: the raise itself
  // ---- Domain::domainChange() from C++, unconditionally (W1, W4) -------------------------
  if (!(reseedMode == 2 && !nodeSetChanged && !connChanged && !specChanged)) {
    dom->domainChange();
    nDomainChanges++;
  }
  MPM_PHASE_END(ph);
  return 0;
}

// `mpm debug tangentCheck <scale> <eps>` (WP-232f footing diagnosis): every particle element
// is armed at the committed state (a prepareStep with no analysis), given the trial nodal
// displacement scale * sin(1.7 j + 0.3) on its DOFs, and its getTangentStiff compared with a
// central difference of getResistingForce. Out: [max rel, mean rel, index of the max, ||K|| of
// it, count]; rel = ||K - K_FD||_F / ||K||_F. Nothing is committed: trial displacements and
// materials are reverted, and the next `mpm analyze` prepares its own step.
int
LadrunoMPMPattern::debugTangentCheck(double scale, double eps, std::vector<double> &out,
                                     std::string &err, int alOverride)
{
  if (alOverride == 1 && !eng.par.antiLocking) {
    err = "-antiLocking fbar needs the deck's `mpm antiLocking fbar` (the corners are set per step)";
    return -1;
  }
  if (theMode == MODE_EXPLICIT) { err = "tangentCheck needs mode implicit or quasiStatic"; return -1; }
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
  double worst = 0.0, sum = 0.0, normW = 0.0;
  int worstP = -1, count = 0;
  for (int p = 0; p < (int)eles.size(); p++) {
    if (eles[p] == 0) continue;
    double rel = 0.0, nk = 0.0;
    if (alOverride == 0) eles[p]->setAntiLockingSuppressed(true);     // WP-232j L2
    const int trc = eles[p]->tangentCheck(scale, eps, rel, nk);
    if (alOverride == 0) eles[p]->setAntiLockingSuppressed(false);
    if (trc < 0) continue;
    sum += rel;
    count++;
    if (rel > worst) { worst = rel; worstP = p; normW = nk; }
  }
  inStep = false;
  armedWindow = false;              // review item 1: no plain `analyze` may run them
  out = {worst, count > 0 ? sum / count : 0.0, (double)worstP, normW, (double)count};
  return 0;
}

// `mpm debug affineCheck <A: ndm x ndm, row-major> <b: ndm>` (WP-232j, L1). Mode F: the
// engine's header-level check. Mode A: every particle element armed at the committed state
// (as tangentCheck), every live grid node given the trial displacement A X + b, every element
// updated once. Out: [max_p |dFbar - dF|_max / |dF - I|_max over averaged particles, max_p
// |deps|_max over them, #averaged, #fallback, #refused]. Nothing is committed.
int
LadrunoMPMPattern::debugAffineCheck(const double *A, const double *b, std::vector<double> &out,
                                    std::string &err)
{
  if (theMode == MODE_EXPLICIT) {               // C2: the production helper of advance()
    std::vector<double> o;
    engineAffineCheck(A, b, o);
    out = {o[0], o[1], o[2], 0.0, o[3], o[4]};
    return 0;
  }
  if (elementsRestored && relinkElements(err) < 0) return -1;
  armedWindow = true;
  const int prc = prepareStep(0.0, err);
  if (prc < 0) {
    armedWindow = false;
    if (prc == -2) aborted = true;
    return -1;
  }
  inStep = true;
  Domain *dom = this->getDomain();
  const int nd = eng.ndm();
  Vector u(nd), zero(nd);
  for (const auto &kv : liveTag) {
    Node *n = dom->getNode(kv.second);
    if (n == 0) continue;
    const Vector &X = n->getCrds();
    for (int a = 0; a < nd; a++) {
      double s = b[a];
      for (int c = 0; c < nd; c++) s += A[nd * a + c] * X(c);
      u(a) = s;
    }
    n->setTrialDisp(u);
  }
  double worst = 0.0, depsMax = 0.0, spur = 0.0;
  int nAvg = 0, nFb = 0, nBad = 0;
  for (LadrunoMPMPoint *e : eles) {
    if (e == 0) continue;
    if (e->update() < 0) { nBad++; continue; }
    const int path = e->antiLockingPath();
    if (path == 2) { nFb++; continue; }
    if (path != 1) continue;
    double rel = 0.0, dm = 0.0;
    e->antiLockingTrial(rel, dm);
    worst = std::max(worst, rel);
    depsMax = std::max(depsMax, dm);
    spur = std::max(spur, e->trialStressChange());
    nAvg++;
  }
  for (const auto &kv : liveTag) {
    Node *n = dom->getNode(kv.second);
    if (n != 0) n->setTrialDisp(zero);
  }
  for (LadrunoMPMPoint *e : eles)
    if (e != 0) e->revertToLastCommit();
  inStep = false;
  armedWindow = false;
  out = {worst, depsMax, (double)nAvg, (double)nFb, (double)nBad, spur};
  return 0;
}

int
LadrunoMPMPattern::setFailDump(const std::string &path, int k)
{
  if (path == "off") {
    failDumpPath.clear();
    LadrunoMPMPoint::diagnose = false;
    return 0;
  }
  if (path.empty() || k < 1) return -1;
  failDumpPath = path;
  failDumpK = k;
  LadrunoMPMPoint::diagnose = true;
  return 0;
}

// The failure diagnosis (the 232f method): the elements still hold the LAST Newton iterate
// (a failed analyze reverts the nodes and materials, not the elements' caches). Nodal residual
// = Sum_e P_e (resisting force minus body force) over free DOFs; a particle's score is the
// largest residual norm among its support nodes. One line per particle:
//   t attempt id x.. score sig.. tangentRatio apex(<1e-2) path(0/1/2) detFc/detF
void
LadrunoMPMPattern::writeFailDump(double t, int attempt)
{
  Domain *dom = this->getDomain();
  const int nd = eng.ndm(), nv = eng.nv();
  std::set<std::pair<int, int>> fixed;
  for (const SpecRow &r : lastSpec) {
    auto it = liveTag.find(r.ext);
    if (it != liveTag.end()) fixed.insert(std::make_pair(it->second, r.dof));
  }
  std::map<int, std::vector<double>> res;
  for (LadrunoMPMPoint *e : eles) {
    if (e == 0 || e->lastResidual() == 0) continue;
    const ID &c = e->getExternalNodes();
    const Vector &Pv = *e->lastResidual();
    if (Pv.Size() != nd * c.Size()) continue;
    for (int k = 0; k < c.Size(); k++) {
      std::vector<double> &r = res[c(k)];
      r.resize(nd, 0.0);
      for (int a = 0; a < nd; a++) r[a] += Pv(nd * k + a);
    }
  }
  std::map<int, double> nodeNorm;
  for (auto &kv : res) {
    double s = 0.0;
    for (int a = 0; a < nd; a++)
      if (!fixed.count(std::make_pair(kv.first, a))) s += kv.second[a] * kv.second[a];
    nodeNorm[kv.first] = std::sqrt(s);
  }
  std::vector<std::pair<double, int>> score;
  for (int p = 0; p < (int)eles.size(); p++) {
    if (eles[p] == 0) continue;
    const ID &c = eles[p]->getExternalNodes();
    double m = 0.0;
    for (int k = 0; k < c.Size(); k++) m = std::max(m, nodeNorm[c(k)]);
    score.push_back(std::make_pair(m, p));
  }
  std::sort(score.begin(), score.end(), [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
    return a.first > b.first || (a.first == b.first && a.second < b.second);
  });
  FILE *f = std::fopen(failDumpPath.c_str(), "a");
  if (f == 0) {
    opserr << "WARNING mpm debug failDump -- cannot open " << failDumpPath.c_str() << "\n";
    return;
  }
  const int n = std::min((int)score.size(), failDumpK);
  for (int r = 0; r < n; r++) {
    const int p = score[r].second;
    LadrunoMPMPoint *e = eles[p];
    std::fprintf(f, "%.12g %d %lld", t, attempt, (long long)eng.P.id[p]);
    for (int a = 0; a < nd; a++) std::fprintf(f, " %.9g", eng.P.x[nd * p + a]);
    std::fprintf(f, " %.6e", score[r].first);
    const double *s = e->trialStress();
    for (int k = 0; k < nv; k++) std::fprintf(f, " %.6e", s[k]);
    const double tr = e->diagTangentRatio();
    std::fprintf(f, " %.6e %d %d %.9g\n", tr, tr < 1e-2 ? 1 : 0, e->antiLockingPath(),
                 e->diagVolumeRatio());
  }
  std::fclose(f);
  (void)dom;
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
  {
    const std::string why = bsplineFaceConflict();   // WP-232i B5, B10
    if (!why.empty()) {
      opserr << "ERROR mpm analyze -- " << why.c_str() << "\n";
      return -1;
    }
  }
  if (!eng.bodyForce.empty() && (int)eng.bodyForce.size() != nd * N) {   // WP-232i B11
    opserr << "ERROR mpm analyze -- the particle body force holds " << (int)eng.bodyForce.size()
           << " values for " << N << " particles (ndm per particle): set it again\n";
    return -1;
  }
  if ((int)pforce.size() != nd * N) pforce.resize(nd * N, 0.0);
  if (!stepped) {
    snapshot = eng.P;
    haveSnapshot = true;
  }
  dtLog.clear();

  LoadControl *lc = 0;
  double ud = 0.0, umin = 0.0, umax = 0.0, ulast = 0.0;
  if (theMode == MODE_QUASISTATIC) {
    lc = static_cast<LoadControl *>(A.si);           // class tag checked by admitAnalysis
    lc->ladrunoGetIncrement(ud, umin, umax, ulast);
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
    // WP-232g (M0b): the wall of one step (prepare, attempts, cuts, output), opened after
    // the loop's exit tests so the last pass of the loop adds no empty call.
    OPS_PROFILE_SCOPE("mpm.a.step");
    int attempt = 0;
    for (;;) {
      attempt++;
      armedWindow = true;
      const int prc = prepareStep(h, err);
      if (prc < 0) {
        if (prc == -2) {
          aborted = true;
          opserr << "ERROR mpm analyze -- " << err.c_str() << ". The grid state was already being"
                    " changed: the run is aborted until `mpm reset` (particles keep the last"
                    " accepted step)\n";
        } else {
          opserr << "ERROR mpm analyze -- " << err.c_str() << "; step refused, nothing changed\n";
        }
        rc = -1;
        break;
      }
      stepT0 = dom->getCurrentTime();
      failNow = (failStep > 0 && failStep == epochSteps + 1 && attempt <= failAttempts)
                ? (failAfterTrial ? 2 : 1) : 0;
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
        lc->ladrunoSetIncrement(h, h, h, -1.0);
        res = A.sa->analyze(1, false);
      }
      inStep = false;
      failNow = 0;
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
      if (!failDumpPath.empty()) writeFailDump(stepT0, attempt);   // WP-232j diagnosis
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
    if (eng.par.antiLocking) eng.countAntiLocking(alPend);     // WP-232j (F1)
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
  if (lc != 0) lc->ladrunoSetIncrement(ud, umin, umax, ulast);   // review item 11
  armedWindow = false;                     // review item 1
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

// WP-232i: 2, + the particle body force. WP-232j (J13): 3, + the anti-locking setting;
// a version-2 pack restores with anti-locking off.
const int MPM_PACK_VERSION = 3;

int
LadrunoMPMPattern::setPackVersion(int v)
{
  if (v != 2 && v != MPM_PACK_VERSION) return -1;
  packVersionOut = (v == MPM_PACK_VERSION) ? 0 : v;
  return 0;
}

int
LadrunoMPMPattern::packState(Vector &out)
{
  Packer k;
  const ladruno_mpm::Grid &G = eng.grid;
  const int ver = packVersionOut ? packVersionOut : MPM_PACK_VERSION;   // WP-232j (L8 hook)
  packVersionOut = 0;                                                   // C6: one save only
  k.put(ver);
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
  k.putv(eng.bodyForce);
  if (ver >= 3) k.put(eng.par.antiLocking);        // WP-232j (J13)
  out.resize((int)k.d.size());
  for (int i = 0; i < (int)k.d.size(); i++) out(i) = k.d[i];
  return 0;
}

int
LadrunoMPMPattern::unpackState(const Vector &in)
{
  Unpacker u(in);
  const int ver = (int)u.get();
  if (ver != MPM_PACK_VERSION && ver != 2) return -1;    // WP-232j (J13): 2 restores as off
  int ndm = (int)u.get();
  double x0[3], h;
  int n[3];
  for (int d = 0; d < 3; d++) x0[d] = u.get();
  h = u.get();
  for (int d = 0; d < 3; d++) n[d] = (int)u.get();
  std::string err;
  eng.clearParticles();
  if (!eng.setGrid(ndm, x0, h, n, err)) return -1;
  {
    const int b = (int)u.get();
    if (b < (int)ladruno_mpm::BASIS_GIMP || b > (int)ladruno_mpm::BASIS_BSPLINE3) return -1;
    eng.par.basis = (ladruno_mpm::Basis)b;
  }
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
  liveTags.clear();
  const int nl = (int)u.get();
  for (int i = 0; i < nl && !u.bad; i++) {
    const long x = (long)u.get();
    liveTag[x] = (int)u.get();
    liveTags.insert(liveTag[x]);
  }
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
  u.getv(eng.bodyForce);
  eng.par.antiLocking = (ver >= 3) ? ((int)u.get() != 0 ? 1 : 0) : 0;
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
