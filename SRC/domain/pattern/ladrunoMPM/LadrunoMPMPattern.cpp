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
// LadrunoMPMPattern — see the header for the design (WP-232, ADR 232 §5.2).

#include "LadrunoMPMPattern.h"

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <LoadPatternIter.h>
#include <NDMaterial.h>
#include <LadrunoJ2.h>
#include <TimeSeries.h>
#include <Vector.h>
#include <Matrix.h>
#include <Parameter.h>
#include <Information.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <Ladruno_mutation.h>
#include <LadrunoMaterialStatus.h>
#include <profiler/ProfilerMacros.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <set>
#include <utility>

// The core mirrors the mutation mode codes so it can stay OpenSees-free.
static_assert(ladruno_mpm::MUT_ZERO == LADRUNO_MUT_ZERO, "MPM mutation code ZERO drifted");
static_assert(ladruno_mpm::MUT_SCALE == LADRUNO_MUT_SCALE, "MPM mutation code SCALE drifted");

// ---------------------------------------------------------------------------
//  The material seam: one owned NDMaterial copy per particle
// ---------------------------------------------------------------------------
class LadrunoMPMPattern::Mats : public ladruno_mpm::MaterialPoints
{
 public:
  Mats(std::vector<NDMaterial *> &m) : mats(m), eps3(3), eps6(6) {}

  int setTrialStrain(int p, const double *e, double *s) override
  {
    NDMaterial *mat = mats[p];
    const int nv = mat->getOrder();
    Vector &eps = (nv == 6) ? eps6 : eps3;
    for (int k = 0; k < nv; k++) eps(k) = e[k];
    const int res = mat->setTrialStrain(eps);
    const Vector &sig = mat->getStress();
    for (int k = 0; k < nv; k++) s[k] = sig(k);
    return (res < 0) ? -1 : 0;
  }
  int commitState(int p) override { return mats[p]->commitState(); }
  int revertToLastCommit(int p) override { return mats[p]->revertToLastCommit(); }
  int pendingRefusals() override
  {
    const int n = ladrunoPendingCommitRefusals();
    if (n > 0) ladrunoClearCommitRefusals();
    return n;
  }
  void clearRefusals() override { ladrunoClearCommitRefusals(); }
  // The INITIAL (elastic) tangent bounds the wave speed (interface Amendment 6,
  // review F6): a softened current tangent (DruckerPrager's apex mCep =
  // 1e-3 mCe) would overestimate the stable Δt by up to 31x.
  double normalModulus(int p) override
  {
    const Matrix &D = mats[p]->getInitialTangent();
    const int nNormal = (D.noRows() == 6) ? 3 : 2;
    double M = 0.0;
    for (int i = 0; i < nNormal; i++) M = std::max(M, D(i, i));
    return M;
  }

 private:
  std::vector<NDMaterial *> &mats;
  Vector eps3, eps6;
};

// ---------------------------------------------------------------------------
LadrunoMPMPattern::LadrunoMPMPattern(int tag)
  : LoadPattern(tag, PATTERN_TAG_LadrunoMPM), adapter(0), series(0), lambda(1.0),
    inStep(false), stepped(false), nextLatticeId(0), haveSnapshot(false),
    epoch(0), writer(0), thickness(1.0), thicknessFixed(false),
    unitL("unspecified"), unitF("unspecified"), unitT("unspecified"), unitsFromFile(false),
    stepsSinceOutput(0), nextOutputTime(0.0), aborted(false)
{
  adapter = new Mats(mats);
}

LadrunoMPMPattern::LadrunoMPMPattern()
  : LoadPattern(0, PATTERN_TAG_LadrunoMPM), adapter(0), series(0), lambda(1.0),
    inStep(false), stepped(false), nextLatticeId(0), haveSnapshot(false),
    epoch(0), writer(0), thickness(1.0), thicknessFixed(false),
    unitL("unspecified"), unitF("unspecified"), unitT("unspecified"), unitsFromFile(false),
    stepsSinceOutput(0), nextOutputTime(0.0), aborted(false)
{
  adapter = new Mats(mats);
}

LadrunoMPMPattern::~LadrunoMPMPattern()
{
  closeOutput();
  detachParameters();    // review F2: no Parameter may keep a pointer into our materials
  for (NDMaterial *m : mats) delete m;
  delete adapter;
}

// Every Parameter that registered a particle material through this host is
// cleaned (Parameter::clean drops its object pointers), so a later
// `updateParameter` after `remove loadPattern` or `mpm reset` is a no-op instead
// of a use-after-free. A parameter the user already removed is skipped.
void
LadrunoMPMPattern::detachParameters(void)
{
  for (const std::pair<Domain *, int> &pt : paramRefs) {
    if (pt.first == 0) continue;
    Parameter *p = pt.first->getParameter(pt.second);
    if (p != 0) p->clean();
  }
  paramRefs.clear();
}

void
LadrunoMPMPattern::nameBody(int body)
{
  if (bodyNames.count(body) == 0)
    bodyNames[body] = "body" + std::to_string(body);
}

// Step 2 of ADR §5.2: Domain::applyLoad(t + Δt) reaches here. Inert outside
// `mpm analyze` (the FE analyses of the same Domain never see the engine).
void
LadrunoMPMPattern::applyLoad(double pseudoTime)
{
  if (!inStep) return;
  if (isConstant != 0)                  // `loadConst` freezes the factor, as for any pattern
    lambda = (series != 0) ? series->getFactor(pseudoTime) : 1.0;
}

void
LadrunoMPMPattern::setTimeSeries(TimeSeries *theSeries)
{
  LoadPattern::setTimeSeries(theSeries);   // the base owns and deletes it
  series = theSeries;
}

double
LadrunoMPMPattern::getLoadFactor(void)
{
  return lambda;
}

int
LadrunoMPMPattern::sendSelf(int, Channel &)
{
  opserr << "ERROR LadrunoMPMPattern " << this->getTag()
         << " -- database save is refused in v1 (ADR 232 §5.2): the particle state"
            " has no restart path yet (phase 1).\n";
  return -1;
}

int
LadrunoMPMPattern::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
  opserr << "ERROR LadrunoMPMPattern -- database restore is refused in v1 (ADR 232 §5.2).\n";
  return -1;
}

void
LadrunoMPMPattern::Print(OPS_Stream &s, int)
{
  const ladruno_mpm::Grid &g = eng.grid;
  s << "LadrunoMPMPattern " << this->getTag() << " (explicit MPM, WP-232): ndm " << g.ndm
    << ", h " << g.h << ", cells";
  for (int d = 0; d < g.ndm; d++) s << " " << g.n[d];
  s << ", " << eng.P.size() << " particles, basis "
    << (eng.par.basis == ladruno_mpm::BASIS_GIMP ? "gimp" : "linear") << ", flip " << eng.par.flip
    << ", damping " << eng.par.damping << endln;
}

LoadPattern *
LadrunoMPMPattern::getCopy(void)
{
  opserr << "WARNING LadrunoMPMPattern::getCopy -- an MPM host cannot be copied\n";
  return 0;
}

// `parameter <t> loadPattern <tag> <args...>`: forwarded to every particle
// material that answers (W12). Each answering material registers itself with
// the Parameter, so `updateParameter` reaches it directly.
int
LadrunoMPMPattern::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  for (NDMaterial *m : mats)
    if (m->setParameter(argv, argc, param) >= 0) res = 0;
  if (res >= 0)
    paramRefs.push_back(std::make_pair(this->getDomain(), param.getTag()));
  if (res < 0)
    opserr << "WARNING LadrunoMPMPattern " << this->getTag()
           << " -- no particle material answers parameter '" << (argc > 0 ? argv[0] : "") << "'\n";
  return res;
}

int
LadrunoMPMPattern::updateParameter(int, Information &)
{
  return 0;
}

int
LadrunoMPMPattern::activateParameter(int)
{
  return 0;
}

// ---------------------------------------------------------------------------
//  Admission (ADR 232 §6, the materials table)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::admitMaterial(int matTag, NDMaterial *&copy, bool &isElastic, std::string &err)
{
  copy = 0;
  NDMaterial *proto = OPS_getNDMaterial(matTag);
  if (proto == 0) {
    err = "no nDMaterial with tag " + std::to_string(matTag);
    return -1;
  }
  const int nd = eng.ndm();
  copy = proto->getCopy(nd == 3 ? "ThreeDimensional" : "PlaneStrain");
  if (copy == 0) {
    err = "nDMaterial " + std::to_string(matTag) + " has no " +
          (nd == 3 ? "ThreeDimensional" : "PlaneStrain") + " form";
    return -1;
  }
  const int ct = copy->getClassTag();
  const bool admitted = ct == ND_TAG_ElasticIsotropicPlaneStrain2d ||
                        ct == ND_TAG_ElasticIsotropicThreeDimensional ||
                        ct == ND_TAG_LadrunoJ2 ||
                        ct == ND_TAG_DruckerPrager3D || ct == ND_TAG_DruckerPragerPlaneStrain;
  if (!admitted) {
    err = "nDMaterial " + std::to_string(matTag) + " (class tag " + std::to_string(ct) +
          ") is not admitted at MPM particles in v1. Admitted: ElasticIsotropic, LadrunoJ2"
          " (without -autoRegularization) and DruckerPrager (ADR 232 §6)";
    delete copy;
    copy = 0;
    return -1;
  }
  if (ct == ND_TAG_LadrunoJ2) {
    LadrunoJ2 *j2 = dynamic_cast<LadrunoJ2 *>(copy);
    if (j2 != 0 && j2->usesAutoRegularization()) {
      err = "LadrunoJ2 " + std::to_string(matTag) + " uses -autoRegularization, which reads"
            " an element length; v1 offers no length at particles (ADR 232 §6, W9)";
      delete copy;
      copy = 0;
      return -1;
    }
  }
  if (copy->getOrder() != eng.nv()) {
    err = "nDMaterial " + std::to_string(matTag) + " has strain order " +
          std::to_string(copy->getOrder()) + ", expected " + std::to_string(eng.nv());
    delete copy;
    copy = 0;
    return -1;
  }
  isElastic = (ct == ND_TAG_ElasticIsotropicPlaneStrain2d ||
               ct == ND_TAG_ElasticIsotropicThreeDimensional);
  return 0;
}

int
LadrunoMPMPattern::addParticle(int64_t id, int body, int matTag, double rho, double V0,
                               const double *x, const double *lp, const double *v,
                               const double *sig0, std::string &err)
{
  if (stepped) {
    err = "particles cannot be added after the first step: v1 fixes the particle set at the"
          " first `mpm analyze` (`mpm reset` returns to it)";
    return -1;
  }
  NDMaterial *copy = 0;
  bool isEl = false;
  if (admitMaterial(matTag, copy, isEl, err) < 0) return -1;
  bool hasS = false;
  if (sig0 != 0)
    for (int k = 0; k < eng.nv(); k++) hasS = hasS || sig0[k] != 0.0;
  if (hasS && !isEl) {
    err = "an initial stress on a non-elastic particle material is refused in v1 (interface Q2):"
          " the offset would bypass the yield check";
    delete copy;
    return -1;
  }
  if (eng.addParticle(id, body, matTag, rho, V0, x, lp, v, sig0, err) < 0) {
    delete copy;
    return -1;
  }
  mats.push_back(copy);
  elastic.push_back(isEl ? 1 : 0);
  if (id >= nextLatticeId) nextLatticeId = id + 1;
  nameBody(body);
  return 0;
}

int
LadrunoMPMPattern::seedBody(int shape, const double *geom, int ppc, int body, int matTag,
                            double rho, std::string &err)
{
  if (stepped) {
    err = "bodies cannot be added after the first step";
    return -1;
  }
  NDMaterial *probe = 0;
  bool isEl = false;
  if (admitMaterial(matTag, probe, isEl, err) < 0) return -1;
  delete probe;
  if (eng.ndm() == 2 && thicknessFixed && thickness != 1.0) {
    err = "lattice seeding uses thickness 1, but a loaded particle set declared " +
          std::to_string(thickness);
    return -1;
  }
  const int before = eng.P.size();
  const int n = eng.seedLattice(shape, geom, ppc, body, matTag, rho, nextLatticeId, err);
  if (n < 0) return -1;
  if (eng.ndm() == 2) thicknessFixed = true;
  nameBody(body);
  for (int p = before; p < eng.P.size(); p++) {
    NDMaterial *copy = 0;
    if (admitMaterial(matTag, copy, isEl, err) < 0) return -1;
    mats.push_back(copy);
    elastic.push_back(isEl ? 1 : 0);
  }
  return n;
}

int
LadrunoMPMPattern::setBodyStress(int body, const std::vector<double> &sig, std::string &err)
{
  if (stepped) {
    err = "initial conditions cannot change after the first step (`mpm reset` first)";
    return -1;
  }
  const int nv = eng.nv();
  if ((int)sig.size() != nv) {
    err = "the stress needs " + std::to_string(nv) + " Voigt components";
    return -1;
  }
  int count = 0;
  for (int p = 0; p < eng.P.size(); p++) {
    if (eng.P.body[p] != body) continue;
    if (!elastic[p]) {
      err = "an initial stress on a non-elastic particle material is refused in v1 (interface Q2)";
      return -1;
    }
  }
  for (int p = 0; p < eng.P.size(); p++) {
    if (eng.P.body[p] != body) continue;
    for (int k = 0; k < nv; k++)
      eng.P.sig0[nv * p + k] = eng.P.sigHat[nv * p + k] = eng.P.sig[nv * p + k] = sig[k];
    count++;
  }
  return count;
}

int
LadrunoMPMPattern::setBodyStressGeostatic(int body, double yTop, double rho, double g,
                                          double K0, std::string &err)
{
  if (stepped) {
    err = "initial conditions cannot change after the first step (`mpm reset` first)";
    return -1;
  }
  if (eng.ndm() != 2) {
    err = "geostatic initial stress is 2D (plane strain) in v1";
    return -1;
  }
  for (int p = 0; p < eng.P.size(); p++)
    if (eng.P.body[p] == body && !elastic[p]) {
      err = "an initial stress on a non-elastic particle material is refused in v1 (interface Q2)";
      return -1;
    }
  int count = 0;
  for (int p = 0; p < eng.P.size(); p++) {
    if (eng.P.body[p] != body) continue;
    const double syy = -rho * g * (yTop - eng.P.x[2 * p + 1]);
    const double s[3] = {K0 * syy, syy, 0.0};
    for (int k = 0; k < 3; k++)
      eng.P.sig0[3 * p + k] = eng.P.sigHat[3 * p + k] = eng.P.sig[3 * p + k] = s[k];
    count++;
  }
  return count;
}

int
LadrunoMPMPattern::setBodyVelocity(int body, const double *v, std::string &err)
{
  if (stepped) {
    err = "initial conditions cannot change after the first step (`mpm reset` first)";
    return -1;
  }
  const int nd = eng.ndm();
  int count = 0;
  for (int p = 0; p < eng.P.size(); p++) {
    if (eng.P.body[p] != body) continue;
    for (int d = 0; d < nd; d++) eng.P.v[nd * p + d] = v[d];
    count++;
  }
  return count;
}

int
LadrunoMPMPattern::setBodyVelocitySine(int body, int comp, int axis, double v0, double beta,
                                       double xref, std::string &err)
{
  if (stepped) {
    err = "initial conditions cannot change after the first step (`mpm reset` first)";
    return -1;
  }
  const int nd = eng.ndm();
  if (comp < 0 || comp >= nd || axis < 0 || axis >= nd) {
    err = "component and axis must be among the model's axes";
    return -1;
  }
  int count = 0;
  for (int p = 0; p < eng.P.size(); p++) {
    if (eng.P.body[p] != body) continue;
    for (int d = 0; d < nd; d++) eng.P.v[nd * p + d] = 0.0;
    eng.P.v[nd * p + comp] = v0 * std::sin(beta * (eng.P.x[nd * p + axis] - xref));
    count++;
  }
  return count;
}

// ---------------------------------------------------------------------------
//  Stepping
// ---------------------------------------------------------------------------
double
LadrunoMPMPattern::criticalTimeStep(void)
{
  return eng.criticalTimeStep(*adapter);
}

bool
LadrunoMPMPattern::domainHasFreeDOF(std::string &which)
{
  Domain *dom = this->getDomain();
  if (dom == 0) return false;
  std::set<std::pair<int, int>> fixed;
  SP_ConstraintIter &sps = dom->getSPs();
  SP_Constraint *sp;
  while ((sp = sps()) != 0) fixed.insert(std::make_pair(sp->getNodeTag(), sp->getDOF_Number()));
  LoadPatternIter &pats = dom->getLoadPatterns();
  LoadPattern *pat;
  while ((pat = pats()) != 0) {
    SP_ConstraintIter &psp = pat->getSPs();
    while ((sp = psp()) != 0) fixed.insert(std::make_pair(sp->getNodeTag(), sp->getDOF_Number()));
  }
  NodeIter &nodes = dom->getNodes();
  Node *nd;
  while ((nd = nodes()) != 0)
    for (int k = 0; k < nd->getNumberDOF(); k++)
      if (fixed.count(std::make_pair(nd->getTag(), k)) == 0) {
        which = "node " + std::to_string(nd->getTag()) + " dof " + std::to_string(k + 1);
        return true;
      }
  return false;
}

int
LadrunoMPMPattern::analyze(int nSteps, double dt)
{
  Domain *dom = this->getDomain();
  if (dom == 0) {
    opserr << "ERROR mpm analyze -- the MPM host is not in a Domain\n";
    return -1;
  }
  if (!(dt > 0.0) || nSteps < 0) {
    opserr << "ERROR mpm analyze -- need nSteps >= 0 and dt > 0\n";
    return -1;
  }
  std::string which;
  if (this->domainHasFreeDOF(which)) {
    opserr << "ERROR mpm analyze -- the Domain has a free FE DOF (" << which.c_str()
           << "). v1 has no FE-MPM seam: an explicit MPM run with free FE DOFs is refused"
              " (ADR 232 §5.2, D2). Fix every FE DOF or remove the FE model.\n";
    return -1;
  }
  if (eng.P.size() == 0) {
    opserr << "ERROR mpm analyze -- no particles\n";
    return -1;
  }
  if (aborted) {
    opserr << "ERROR mpm analyze -- an earlier run was aborted by a material refusal at commit"
              " (ADR 232 W17); the particle state is not consistent. `mpm reset` first.\n";
    return -1;
  }
  if ((int)mats.size() != eng.P.size() || (int)elastic.size() != eng.P.size()) {
    opserr << "ERROR mpm analyze -- internal: " << (int)mats.size() << " particle materials for "
           << eng.P.size() << " particles; refused\n";
    return -1;
  }
  // The `mpm reset` snapshot is the state before the first ACCEPTED step: it is
  // retaken by every analyze until a step is accepted, so a refused or
  // zero-step analyze followed by more bodies is covered (review F1).
  if (!stepped) {
    snapshot = eng.P;
    haveSnapshot = true;
  }
  // .ladruno particle group (232a §2.8): a new stage at the first `mpm analyze`
  // of an epoch, and its output 0 is the state before the first step.
  if (writer != 0 && writer->isOpen() &&
      (!writer->stageOpen() || writer->stageEpoch() != epoch)) {
    std::string err;
    if (openOutputStage(err) < 0 || writeOutputNow(err) < 0) {
      opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
      return -1;
    }
  }

  for (int s = 0; s < nSteps; s++) {
    OPS_PROFILE_SCOPE("step");
    // ---- 1. admissibility and the critical time step ----------------------
    const int bad = eng.firstInadmissible();
    if (bad >= 0) {
      opserr << "ERROR mpm analyze -- particle " << std::to_string((long long)eng.P.id[bad]).c_str()
             << " is not admissible (its centre left the grid, or lp is out of (0, h/2]);"
                " step refused, nothing changed\n";
      return -1;
    }
    const double dtc = eng.criticalTimeStep(*adapter);
    if (dt > eng.par.safety * dtc) {
      opserr << "ERROR mpm analyze -- dt = " << dt << " exceeds s * dt_crit = "
             << eng.par.safety << " * " << dtc << " = " << eng.par.safety * dtc
             << "; step refused, nothing changed (ADR 232 §5.2 step 1)\n";
      return -1;
    }
    // ---- 2. loads and time -------------------------------------------------
    const double t0 = dom->getCurrentTime();
    inStep = true;
    dom->applyLoad(t0 + dt);
    inStep = false;
    // ---- 3-7. the core, with no active element (W9) ----------------------
    Element *saved = ops_TheActiveElement;
    ops_TheActiveElement = 0;
    std::string err;
    int res;
    { OPS_PROFILE_SCOPE("mpm.advance");
      res = eng.advance(dt, lambda, *adapter, err); }
    ops_TheActiveElement = saved;
    if (res == -1) {
      dom->revertToLastCommit();
      opserr << "ERROR mpm analyze -- step failed at t = " << t0 + dt << ": " << err.c_str()
             << "; particles keep their committed state\n";
      return -1;
    }
    if (res == -2) {
      dom->revertToLastCommit();
      aborted = true;      // W17: only `mpm reset` makes the state usable again
      if (writer != 0) writer->flush();
      opserr << "ERROR mpm analyze -- " << err.c_str() << ". Particles were not moved; the run is"
                " aborted until `mpm reset`.\n";
      return -1;
    }
    stepped = true;
    // ---- 8. Domain::commit: committed time, commit tag, recorders ----------
    { OPS_PROFILE_SCOPE("commit");
      if (dom->commit() < 0) {
        opserr << "ERROR mpm analyze -- Domain::commit failed at t = " << t0 + dt << "\n";
        return -1;
      } }
    // ---- 9. the profiler's per-step row --------------------------------------
    OPS_PROFILE_STEP(dom->getCommitTag(), dom->getCurrentTime(), dt, 1);
    // ---- particle output, after Domain::commit: positions and stresses from
    // the same committed step (232a §2.8 synchrony, W11) --------------------------
    if (writer != 0 && writer->stageOpen()) {
      stepsSinceOutput++;
      const ladruno_mpm_io::WriterConfig &cfg = writer->config();
      const double t = dom->getCurrentTime();
      bool due;
      if (cfg.dtOut > 0.0)
        due = t >= nextOutputTime - 1e-12 * std::max(1.0, std::fabs(nextOutputTime));
      else
        due = stepsSinceOutput >= cfg.every;
      if (due && writeOutputNow(err) < 0) {
        opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
        return -1;
      }
    }
  }
  // the last step is never lost (232a §2.8)
  if (writer != 0 && writer->stageOpen() && stepsSinceOutput > 0) {
    std::string err;
    if (writeOutputNow(err) < 0) {
      opserr << "ERROR mpm analyze -- particle output: " << err.c_str() << "\n";
      return -1;
    }
  }
  if (writer != 0) writer->flush();
  dom->flushRecorders();
  return 0;
}

int
LadrunoMPMPattern::reset(void)
{
  epoch++;                          // PARTICLE_EPOCH: the next output opens a new stage
  aborted = false;
  if (!haveSnapshot) return 0;      // nothing has run: already initial
  // Particles added after the snapshot are dropped WITH their materials, flags
  // and ids (review F1): the set returns to exactly the snapshot.
  const size_t keep = (size_t)snapshot.size();
  if (mats.size() > keep) {
    if (!paramRefs.empty()) {
      opserr << "WARNING mpm reset -- particles added after the first step's snapshot are dropped;"
                " the parameters set through `parameter ... loadPattern` are detached\n";
      detachParameters();
    }
    for (size_t q = keep; q < mats.size(); q++) delete mats[q];
    mats.resize(keep);
    elastic.resize(keep);
  }
  eng.setParticles(snapshot);
  nextLatticeId = 0;
  for (int64_t pid : eng.P.id) nextLatticeId = std::max<int64_t>(nextLatticeId, pid + 1);
  for (NDMaterial *m : mats) m->revertToStart();
  stepped = false;
  haveSnapshot = false;     // recaptured at the next analyze (particles may be added again)
  lambda = 1.0;
  return 0;
}

// ---------------------------------------------------------------------------
//  Phase c2: Contract 1 in, Contract 2 out (232a)
// ---------------------------------------------------------------------------
int
LadrunoMPMPattern::loadParticleSet(const std::string &path)
{
  const std::string ctx = "mpm particles -file " + path;
  if (stepped) {
    opserr << "WARNING " << ctx.c_str() << " -- particles cannot be added after the first step\n";
    return -1;
  }
  ladruno_mpm_io::ParticleSet S;
  std::vector<std::string> problems;
  if (ladruno_mpm_io::readParticleSet(path, S, problems) < 0) {
    for (const std::string &p : problems)
      opserr << "WARNING " << ctx.c_str() << " -- " << p.c_str() << "\n";
    opserr << "WARNING " << ctx.c_str() << " -- the particle set is refused ("
           << (int)problems.size() << " problem(s)); nothing was loaded\n";
    return -1;
  }
  // ---- engine refusals: the set must fit this engine (232a §1.2, §0.1) -------
  const int nd = eng.ndm(), nv = eng.nv();
  const ladruno_mpm::Grid &g = eng.grid;
  std::vector<std::string> refuse;
  if (S.ndm != nd)
    refuse.push_back("C1-NDM: the set has NDM " + std::to_string(S.ndm) + ", the model ndm is " +
                     std::to_string(nd));
  else {
    bool equal = true;
    for (int d = 1; d < nd; d++) equal = equal && S.cellSize[d] == S.cellSize[0];
    if (!equal)
      refuse.push_back("the set's CELL_SIZE components differ; v1 has one cell size h for every axis");
    bool same = equal;
    for (int d = 0; d < nd; d++) {
      const double tol = 1e-12 * std::max(1.0, std::fabs(g.x0[d]) + g.n[d] * g.h);
      same = same && std::fabs(S.origin[d] - g.x0[d]) <= tol && S.numCells[d] == g.n[d] &&
             std::fabs(S.cellSize[d] - g.h) <= 1e-12 * g.h;
    }
    if (!same)
      refuse.push_back("the set's GRID (ORIGIN, CELL_SIZE, NUM_CELLS) differs from the deck's `mpm grid`");
    const std::string basis = eng.par.basis == ladruno_mpm::BASIS_GIMP ? "gimp" : "linear";
    if (S.basis != basis)
      refuse.push_back("the set's GRID/BASIS '" + S.basis + "' differs from the deck's '" + basis + "'");
  }
  if (nd == 2 && thicknessFixed && S.thickness != thickness)
    refuse.push_back("PLANE_STRAIN_THICKNESS " + std::to_string(S.thickness) +
                     " differs from the engine's " + std::to_string(thickness));
  if (unitsFromFile && (S.unitLength != unitL || S.unitForce != unitF || S.unitTime != unitT))
    refuse.push_back("the set's units differ from a set already loaded");
  for (size_t b = 0; b < S.bodyId.size(); b++) {
    auto it = bodyNames.find(S.bodyId[b]);
    if (it != bodyNames.end() && it->second != S.bodyName[b] &&
        it->second != "body" + std::to_string(S.bodyId[b]))
      refuse.push_back("body " + std::to_string(S.bodyId[b]) + " is already named '" + it->second +
                       "', the set names it '" + S.bodyName[b] + "'");
  }
  std::map<int, bool> matElastic;
  for (int tag : S.mat) {
    if (matElastic.count(tag)) continue;
    NDMaterial *probe = 0;
    bool isEl = false;
    std::string err;
    if (admitMaterial(tag, probe, isEl, err) < 0) {
      refuse.push_back("C1-MAT: " + err);
      matElastic[tag] = true;
      continue;
    }
    delete probe;
    matElastic[tag] = isEl;
  }
  const int N = (int)S.n;
  for (int p = 0; p < N && refuse.size() < 20; p++) {
    if (eng.hasId(S.id[p]))
      refuse.push_back("particle id " + std::to_string(S.id[p]) + " is already used in this epoch");
    if (!matElastic[S.mat[p]]) {
      bool hasS = S.stress33[p] != 0.0;
      for (int k = 0; k < nv; k++) hasS = hasS || S.stress[nv * p + k] != 0.0;
      if (hasS)
        refuse.push_back("particle " + std::to_string(S.id[p]) + ": an initial stress (STRESS or a "
                         "nonzero STRESS33) on a non-elastic material is refused in v1 (232b Q2)");
    }
  }
  if (!refuse.empty()) {
    for (const std::string &r : refuse)
      opserr << "WARNING " << ctx.c_str() << " -- " << r.c_str() << "\n";
    opserr << "WARNING " << ctx.c_str() << " -- the particle set is refused; nothing was loaded\n";
    return -1;
  }
  // ---- load, all of it --------------------------------------------------------
  for (int p = 0; p < N; p++) {
    std::string err;
    const double rho = S.mass[p] / S.volume[p];
    if (addParticle(S.id[p], S.body[p], S.mat[p], rho, S.volume[p], &S.x[nd * p], &S.lp[nd * p],
                    &S.v[nd * p], &S.stress[nv * p], err) < 0) {
      opserr << "ERROR " << ctx.c_str() << " -- particle " << std::to_string(S.id[p]).c_str()
             << ": " << err.c_str() << " (the set is partly loaded: wipe and fix the file)\n";
      return -1;
    }
    eng.P.m.back() = S.mass[p];     // the set's mass exactly (MASS, or DENSITY x VOLUME once)
  }
  for (size_t b = 0; b < S.bodyId.size(); b++) bodyNames[S.bodyId[b]] = S.bodyName[b];
  if (nd == 2) { thickness = S.thickness; thicknessFixed = true; }
  unitL = S.unitLength; unitF = S.unitForce; unitT = S.unitTime;
  if (!S.unitSystem.empty()) unitSystem = S.unitSystem;
  unitsFromFile = true;
  sourceSetIds.push_back(S.setId);
  std::string abs = path;
  try { abs = std::filesystem::absolute(std::filesystem::path(path)).string(); } catch (...) {}
  sourcePaths.push_back(abs);
  opserr << "mpm particles: " << N << " particles from " << abs.c_str() << " (SET_ID "
         << S.setId.c_str() << ")\n";
  return N;
}

int
LadrunoMPMPattern::setUnits(const std::string &L, const std::string &F, const std::string &T,
                            const std::string &system, std::string &err)
{
  if (L.empty() || F.empty() || T.empty()) { err = "units must be non-empty"; return -1; }
  if (unitsFromFile && (L != unitL || F != unitF || T != unitT)) {
    err = "a loaded particle set declared " + unitL + ", " + unitF + ", " + unitT;
    return -1;
  }
  unitL = L; unitF = F; unitT = T;
  if (!system.empty()) unitSystem = system;
  return 0;
}

int
LadrunoMPMPattern::setOutput(const ladruno_mpm_io::WriterConfig &cfg, std::string &err)
{
  closeOutput();
  writer = new ladruno_mpm_io::ParticleWriter();
  if (writer->open(cfg, eng.ndm(), err) < 0) {
    delete writer;
    writer = 0;
    return -1;
  }
  stepsSinceOutput = 0;
  return 0;
}

void
LadrunoMPMPattern::closeOutput(void)
{
  if (writer != 0) {
    writer->close();
    delete writer;
    writer = 0;
  }
}

int
LadrunoMPMPattern::openOutputStage(std::string &err)
{
  Domain *dom = this->getDomain();
  const int nd = eng.ndm();
  const ladruno_mpm::Particles &P = eng.P;
  idOrder.resize(P.size());
  for (int p = 0; p < P.size(); p++) idOrder[p] = p;
  std::sort(idOrder.begin(), idOrder.end(), [&](int a, int b) { return P.id[a] < P.id[b]; });

  ladruno_mpm_io::StageMeta m;
  m.epoch = epoch;
  m.patternTag = this->getTag();
  m.mode = "explicit";
  m.thickness = thickness;
  m.unitLength = unitL; m.unitForce = unitF; m.unitTime = unitT; m.unitSystem = unitSystem;
  m.sourceSetIds = sourceSetIds;
  for (int d = 0; d < nd; d++) {
    m.origin[d] = eng.grid.x0[d];
    m.cellSize[d] = eng.grid.h;
    m.numCells[d] = eng.grid.n[d];
  }
  m.basis = eng.par.basis == ladruno_mpm::BASIS_GIMP ? "gimp" : "linear";
  for (const auto &kv : bodyNames) { m.bodyId.push_back(kv.first); m.bodyName.push_back(kv.second); }
  for (int q : idOrder) {
    m.id.push_back((long long)P.id[q]);
    m.body.push_back(P.body[q]);
    m.mat.push_back(P.mat[q]);
    m.mass.push_back(P.m[q]);
  }
  if (writer->openStage(m, dom->getCommitTag(), dom->getCurrentTime(), err) < 0) return -1;
  stepsSinceOutput = 0;
  nextOutputTime = dom->getCurrentTime();
  return 0;
}

int
LadrunoMPMPattern::writeOutputNow(std::string &err)
{
  Domain *dom = this->getDomain();
  const int nd = eng.ndm(), nv = eng.nv(), nn = nd * nd;
  const ladruno_mpm::Particles &P = eng.P;
  const ladruno_mpm_io::WriterConfig &cfg = writer->config();
  std::vector<long long> ids;
  ids.reserve(idOrder.size());
  for (int q : idOrder) ids.push_back((long long)P.id[q]);
  std::vector<std::vector<double> > rows(cfg.fields.size());
  for (size_t k = 0; k < cfg.fields.size(); k++) {
    const std::string &f = cfg.fields[k];
    std::vector<double> &r = rows[k];
    for (int q : idOrder) {
      if (f == "POSITION") r.insert(r.end(), &P.x[nd * q], &P.x[nd * q] + nd);
      else if (f == "VELOCITY") r.insert(r.end(), &P.v[nd * q], &P.v[nd * q] + nd);
      else if (f == "STRESS") r.insert(r.end(), &P.sig[nv * q], &P.sig[nv * q] + nv);
      else if (f == "STRAIN") r.insert(r.end(), &P.eps[nv * q], &P.eps[nv * q] + nv);
      else if (f == "VOLUME") r.push_back(P.V[q]);
      else if (f == "GIMP_HALF_LENGTH") r.insert(r.end(), &P.lp[nd * q], &P.lp[nd * q] + nd);
      else if (f == "DEFORMATION_GRADIENT") r.insert(r.end(), &P.F[nn * q], &P.F[nn * q] + nn);
    }
  }
  if (writer->writeOutput(dom->getCommitTag(), dom->getCurrentTime(), ids, rows, err) < 0) return -1;
  stepsSinceOutput = 0;
  if (cfg.dtOut > 0.0) {
    const double t = dom->getCurrentTime();
    while (nextOutputTime <= t + 1e-12 * std::max(1.0, std::fabs(t))) nextOutputTime += cfg.dtOut;
  }
  return 0;
}

// ---------------------------------------------------------------------------
LadrunoMPMPattern *
ladrunoMPMFindHost(Domain *theDomain)
{
  if (theDomain == 0) return 0;
  LoadPatternIter &it = theDomain->getLoadPatterns();
  LoadPattern *p;
  while ((p = it()) != 0)
    if (p->getClassTag() == PATTERN_TAG_LadrunoMPM)
      return static_cast<LadrunoMPMPattern *>(p);
  return 0;
}
