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

// Ladruno WP-188 -- factorization reuse across time steps. See LadrunoFactorReuse.h
// and Ladruno_implementation/188_ladruno_factor_reuse_adr.md.
//
// The step driver below is AcceleratedNewton::solveCurrentStep with three changes:
//   1. the step-start formTangent is skipped while a kept factor is still valid;
//   2. a kept factor is replaced mid-step when it converges too slowly (-reuseIter,
//      -reuseRate), exactly where KrylovAccelerator::updateTangent would replace it;
//   3. a step that fails after starting on a kept factor is undone and retried once
//      with a fresh tangent before the failure is reported.
// The residual, the update and the convergence test are the vanilla calls, in the
// vanilla order, with the vanilla profiler scope names.

#include "LadrunoFactorReuse.h"
#include "LadrunoFactorReuseMPI.h"

#include <AcceleratedNewton.h>
#include <Accelerator.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <TransientIntegrator.h>
#include <LoadControl.h>
#include <LinearSOE.h>
#include <ConvergenceTest.h>
#include <Vector.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <elementAPI.h>
#include <profiler/ProfilerMacros.h>

#include <chrono>   // WP-193
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace {

// Relative tolerance on the step size: a factor formed for dt serves a step of dt'
// only if |dt' - dt| <= kDtRelTol * max(|dt|, |dt'|). Domain times are sums
// t_c + dt, so the recovered dt carries ~1e-16 * t / dt of round-off (3e-12 for a
// 70 s record at dt = 2.5 ms); 1e-8 is far above that and far below any real cut.
const double kDtRelTol = 1.0e-8;

// WP-193: wall clock for the cost model (monotonic; seconds)
double nowSec(void)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// WP-193 test hook: burn `sec` seconds of wall time (a busy wait, not a sleep, so
// the Windows timer granularity does not round it), standing in for a slow factorization.
void burnSeconds(double sec)
{
  if (!(sec > 0.0)) return;
  const double t0 = nowSec();
  volatile double sink = 0.0;
  while (nowSec() - t0 < sec)
    sink = sink + 1.0;
}

// EMA weight of a new sample in the WP-193 running means
const double kCostAlpha = 0.25;
// forgetting factor of the WP-193 n_K(age) regression (one per kept-factor step)
const double kCostForget = 0.98;

void ema(double &m, bool &have, double x)
{
  if (!have) { m = x; have = true; }
  else m += kCostAlpha * (x - m);
}

bool parseOnOff(const char *s, bool &v)
{
  if (s == 0) return false;
  if (strcmp(s, "on") == 0 || strcmp(s, "1") == 0 || strcmp(s, "true") == 0) { v = true; return true; }
  if (strcmp(s, "off") == 0 || strcmp(s, "0") == 0 || strcmp(s, "false") == 0) { v = false; return true; }
  // openseespy stringifies a numeric 1 / 0 through OPS_GetStringFromAll ("1.000...")
  char *end = 0;
  double x = strtod(s, &end);
  if (end != s && *end == 0 && (x == 0.0 || x == 1.0)) { v = (x == 1.0); return true; }
  return false;
}

bool parseIntArg(const char *s, int &v)
{
  if (s == 0 || *s == 0) return false;
  char *end = 0;
  long x = strtol(s, &end, 10);
  if (end == s || *end != 0) return false;
  v = (int)x;
  return true;
}

bool parseDoubleArg(const char *s, double &v)
{
  if (s == 0 || *s == 0) return false;
  char *end = 0;
  double x = strtod(s, &end);
  if (end == s || *end != 0) return false;
  v = x;
  return true;
}

// one value check shared by both parsers; returns 1 ok, -1 bad value
int store(const char *flag, LadrunoFactorReuseOptions &o, bool isInt, int iv, double dv)
{
  if (strcmp(flag, "-reuseFactor") == 0) {
    if (!isInt || iv < 0) { opserr << "WARNING KrylovNewton -reuseFactor needs an integer m >= 0 (0 = off)\n"; return -1; }
    o.maxSteps = iv; o.sawFactor = true; return 1;
  }
  if (strcmp(flag, "-reuseIter") == 0) {
    if (!isInt || iv < 1) { opserr << "WARNING KrylovNewton -reuseIter needs an integer k >= 1\n"; return -1; }
    o.maxIter = iv; o.sawOther = true; return 1;
  }
  if (strcmp(flag, "-reuseRate") == 0) {
    if (isInt) dv = iv;
    if (!(dv >= 0.0)) { opserr << "WARNING KrylovNewton -reuseRate needs a number r >= 0 (0 = off)\n"; return -1; }
    o.maxRate = dv; o.sawOther = true; return 1;
  }
  if (strcmp(flag, "-reuseReport") == 0) {
    if (!isInt || iv < 0) { opserr << "WARNING KrylovNewton -reuseReport needs an integer n >= 0\n"; return -1; }
    o.report = iv; o.sawOther = true; return 1;
  }
  // WP-193 cost-aware trigger tuning
  if (strcmp(flag, "-reuseCostSafety") == 0) {
    if (isInt) dv = iv;
    if (!(dv > 0.0)) { opserr << "WARNING KrylovNewton -reuseCostSafety needs a number s > 0\n"; return -1; }
    o.costSafety = dv; o.sawOther = true; o.sawCostTuning = true; return 1;
  }
  if (strcmp(flag, "-reuseCostWarmup") == 0) {
    if (!isInt || iv < 0) { opserr << "WARNING KrylovNewton -reuseCostWarmup needs an integer w >= 0\n"; return -1; }
    o.costWarmup = iv; o.sawOther = true; o.sawCostTuning = true; return 1;
  }
  if (strcmp(flag, "-reuseCostProbe") == 0) {
    if (!isInt || iv < 0) { opserr << "WARNING KrylovNewton -reuseCostProbe needs an integer p >= 0 (0 = never)\n"; return -1; }
    o.costProbe = iv; o.sawOther = true; o.sawCostTuning = true; return 1;
  }
  return 0;
}

// WP-193: -reuseCost on|off (a word, not a number)
int storeCost(const char *val, LadrunoFactorReuseOptions &o)
{
  bool v = false;
  if (!parseOnOff(val, v)) {
    opserr << "WARNING KrylovNewton -reuseCost needs on|off\n";
    return -1;
  }
  o.cost = v; o.sawOther = true;
  return 1;
}

bool isReuseFlag(const char *flag)
{
  return flag != 0 && (strcmp(flag, "-reuseFactor") == 0 || strcmp(flag, "-reuseIter") == 0 ||
                       strcmp(flag, "-reuseRate") == 0 || strcmp(flag, "-reuseReport") == 0 ||
                       // WP-193
                       strcmp(flag, "-reuseCost") == 0 || strcmp(flag, "-reuseCostSafety") == 0 ||
                       strcmp(flag, "-reuseCostWarmup") == 0 || strcmp(flag, "-reuseCostProbe") == 0);
}

} // namespace

bool
LadrunoFactorReuseOptions::consistent(void) const
{
  if (sawOther && !sawFactor) {
    opserr << "WARNING KrylovNewton: -reuseIter / -reuseRate / -reuseReport need -reuseFactor m "
              "(WP-188); refusing the algorithm\n";
    return false;
  }
  if (sawCostTuning && !cost) {   // WP-193
    opserr << "WARNING KrylovNewton: -reuseCostSafety / -reuseCostWarmup / -reuseCostProbe need "
              "-reuseCost on (WP-193); refusing the algorithm\n";
    return false;
  }
  return true;
}

LadrunoFactorReuse::LadrunoFactorReuse(const LadrunoFactorReuseOptions &o)
  : opts(o), valid(false), fSOE(0), fIntegrator(0), fNumEqn(-1), fDt(0.0), fAge(0),
    haveLastTime(false), lastConvergedTime(0.0), dUstep(0),
    nSteps(0), nStepsReused(0), nFormStart(0), nFormIter(0), nFormRate(0), nFormKrylov(0),
    nFormRetry(0), nRetries(0), nRetryFailed(0), nIterTotal(0), nIterReused(0),
    rank(ladrunoFactorReuseRank()), announced(false), testVoteRank(-1),
    emaTF(0.0), emaTI(0.0), emaNF(0.0), haveTF(false), haveTI(false), haveNF(false),
    rSw(0.0), rSa(0.0), rSn(0.0), rSaa(0.0), rSan(0.0), cycSteps(0), costRefusals(0),
    stepTF(0.0), stepTI(0.0), stepNTF(0), stepNTI(0), testFactorDelay(0.0), testIterDelay(0.0),
    nCostForm(0), nCostKeep(0), nCostIter(0), nCostProbe(0)
{
  // WP-193 test hooks: every factorization formed by this driver also burns
  // LADRUNO_WP193_TEST_FACTOR_DELAY_MS milliseconds of wall time, and every
  // iteration LADRUNO_WP193_TEST_ITER_DELAY_MS (the "slow factorization" and "slow
  // iteration" arms of the test module: a large model's factorization, an expensive
  // state determination). Read once; unset = 0.
  const char *d = getenv("LADRUNO_WP193_TEST_FACTOR_DELAY_MS");
  if (d != 0 && *d != 0)
    testFactorDelay = atof(d) * 1.0e-3;
  d = getenv("LADRUNO_WP193_TEST_ITER_DELAY_MS");
  if (d != 0 && *d != 0)
    testIterDelay = atof(d) * 1.0e-3;

  // Test hook (WP-188 MPI gate): the rank named by LADRUNO_WP188_TEST_VOTE_RANK votes
  // "form" at every step start. With the reduction live every rank then forms; were
  // it compiled out (the OBJECT-library trap, LadrunoFactorReuseMPI.h) the other
  // ranks would keep their factor and the counts would disagree. Read once.
  const char *e = getenv("LADRUNO_WP188_TEST_VOTE_RANK");
  if (e != 0 && *e != 0)
    testVoteRank = atoi(e);
}

LadrunoFactorReuse::~LadrunoFactorReuse()
{
  if (nSteps > 0)
    this->printSummary(opserr, "final");
  if (dUstep != 0)
    delete dUstep;
}

int
LadrunoFactorReuse::parseOPS(const char *flag, LadrunoFactorReuseOptions &o)
{
  if (!isReuseFlag(flag)) return 0;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING KrylovNewton " << flag << " needs a value\n";
    return -1;
  }
  if (strcmp(flag, "-reuseCost") == 0) {   // WP-193: a word (typed under openseespy)
    char buf[64];
    const char *val = OPS_GetStringFromAll(buf, 64);
    return storeCost(val, o);
  }
  const bool wantDouble = (strcmp(flag, "-reuseRate") == 0 || strcmp(flag, "-reuseCostSafety") == 0);
  int numdata = 1;
  if (wantDouble) {
    double v = 0.0;
    if (OPS_GetDoubleInput(&numdata, &v) < 0) {
      opserr << "WARNING KrylovNewton " << flag << ": invalid value\n";
      return -1;
    }
    return store(flag, o, false, 0, v);
  }
  int v = 0;
  if (OPS_GetIntInput(&numdata, &v) < 0) {
    opserr << "WARNING KrylovNewton " << flag << ": invalid value\n";
    return -1;
  }
  return store(flag, o, true, v, 0.0);
}

int
LadrunoFactorReuse::parseArgv(int argc, const char **argv, int &i, LadrunoFactorReuseOptions &o)
{
  const char *flag = argv[i];
  if (!isReuseFlag(flag)) return 0;
  if (i + 1 >= argc) {
    opserr << "WARNING KrylovNewton " << flag << " needs a value\n";
    return -1;
  }
  const char *val = argv[++i];
  if (strcmp(flag, "-reuseCost") == 0)   // WP-193
    return storeCost(val, o);
  if (strcmp(flag, "-reuseRate") == 0 || strcmp(flag, "-reuseCostSafety") == 0) {
    double v = 0.0;
    if (!parseDoubleArg(val, v)) {
      opserr << "WARNING KrylovNewton " << flag << ": invalid value '" << val << "'\n";
      return -1;
    }
    return store(flag, o, false, 0, v);
  }
  int v = 0;
  if (!parseIntArg(val, v)) {
    opserr << "WARNING KrylovNewton " << flag << ": invalid value '" << val << "'\n";
    return -1;
  }
  return store(flag, o, true, v, 0.0);
}

int
LadrunoFactorReuse::agree(int flag)
{
  int v = flag;
  ladrunoFactorReuseAgreeMax(&v, 1);
  return v;
}

void
LadrunoFactorReuse::invalidate(void)
{
  valid = false;
  fSOE = 0;
  fIntegrator = 0;
  fNumEqn = -1;
  fAge = 0;
}

void
LadrunoFactorReuse::factorFormed(LinearSOE *soe, IncrementalIntegrator *integ, int n, double dt)
{
  valid = true;
  fSOE = soe;
  fIntegrator = integ;
  fNumEqn = n;
  fDt = dt;
  fAge = 0;
}

void
LadrunoFactorReuse::printSummary(OPS_Stream &s, const char *when, bool allRanks)
{
  if (rank != 0 && !allRanks) return;
  const long nForm = nFormStart + nFormIter + nFormRate + nFormKrylov + nFormRetry;
  const double pct = (nSteps > 0) ? 100.0 * (double)nStepsReused / (double)nSteps : 0.0;
  s << "KrylovNewton factor reuse (WP-188) " << when;
  if (allRanks) s << " [rank " << rank << "]";
  s << ": steps " << (int)nSteps
    << ", started on a kept factor " << (int)nStepsReused << " (" << pct << " %)"
    << ", factorizations " << (int)nForm << " = start " << (int)nFormStart
    << " + iter " << (int)nFormIter << " + rate " << (int)nFormRate
    << " + krylov " << (int)nFormKrylov << " + retry " << (int)nFormRetry
    << ", retries " << (int)nRetries << " (failed " << (int)nRetryFailed << ")"
    << ", iterations " << (int)nIterTotal << " (on kept factors " << (int)nIterReused << ")"
    << endln;
  if (opts.cost)   // WP-193
    this->printCostSummary(s, when, allRanks);
}

// ---------------------------------------------------------------- WP-193 cost model

void
LadrunoFactorReuse::printCostSummary(OPS_Stream &s, const char *when, bool allRanks)
{
  const double ratio = (haveTF && haveTI && emaTI > 0.0) ? emaTF / emaTI : 0.0;
  s << "KrylovNewton reuse cost model (WP-193) " << when;
  if (allRanks) s << " [rank " << rank << "]";
  s << ": t_factor " << (haveTF ? emaTF : -1.0) << " s, t_iter " << (haveTI ? emaTI : -1.0)
    << " s, ratio " << ratio << ", n_fresh " << (haveNF ? emaNF : -1.0)
    << ", n_kept1 " << (rSw > 0.0 ? this->projectKept(1) : -1.0)
    << ", slope " << (rSw > 0.0 ? this->projectKept(2) - this->projectKept(1) : -1.0)
    << ", starts kept by cost " << (int)nCostKeep << ", starts formed by cost " << (int)nCostForm
    << ", in-step cost refreshes " << (int)nCostIter << ", probes " << (int)nCostProbe
    << endln;
}

bool
LadrunoFactorReuse::costReady(void) const
{
  return opts.cost && nSteps > opts.costWarmup && haveTF && haveTI && haveNF && emaTI > 0.0;
}

// Effective iterations of a step started on a factor of age `age` (>= 1):
// n_K(a) = c0 + c1 * a, a least-squares line through the kept-factor steps seen so
// far, exponentially forgotten (weight kCostForget per new sample, ~50 samples), so
// the slope follows the current phase of the response. c1 is clamped >= 0 (a factor
// does not get better with age); with a single distinct age the line is flat; with
// nothing measured it is n_F (optimistic: the first kept steps get tried).
double
LadrunoFactorReuse::projectKept(int age) const
{
  if (age < 1) age = 1;
  if (!(rSw > 0.0)) return emaNF;
  const double det = rSw * rSaa - rSa * rSa;
  double c1 = 0.0;
  if (det > 1.0e-9 * rSw * rSw)
    c1 = (rSw * rSan - rSa * rSn) / det;
  if (c1 < 0.0) c1 = 0.0;
  const double c0 = (rSn - c1 * rSa) / rSw;
  return c0 + c1 * (double)age;
}

// Mean cost per step, in units of t_I, of a cycle of S steps under the model: the
// factorization plus n_F for the step that formed it plus n_K(1..S-1).
double
LadrunoFactorReuse::cycleMean(int S) const
{
  if (S < 1) S = 1;
  double sum = emaTF / emaTI + emaNF;
  if (S > 1) {
    const double n1 = this->projectKept(1);
    const double c1 = this->projectKept(2) - n1;
    const double m = (double)(S - 1);
    sum += m * (n1 - c1) + c1 * m * (m + 1.0) / 2.0;   // sum_{a=1}^{S-1} (c0 + c1 a)
  }
  return sum / (double)S;
}

// Step end: reduce this step's timing samples (MAX over the ranks: the slowest rank
// sets the wall time), fold them into the running means, then the step's iteration
// count as effective iterations (each factorization beyond the step's own start one
// counts t_F / t_I iterations). Called on every rank at the same two exits, so the
// collective stays in lock-step. `iters < 0` = a failed step (times only).
void
LadrunoFactorReuse::costStepEnd(bool startedFresh, int startAge, int iters, int forms)
{
  if (!opts.cost) return;
  double v[2];
  v[0] = (stepNTF > 0) ? stepTF / (double)stepNTF : -1.0;
  v[1] = (stepNTI > 0) ? stepTI / (double)stepNTI : -1.0;
  ladrunoFactorReuseAgreeMaxDouble(v, 2);
  if (v[0] >= 0.0) ema(emaTF, haveTF, v[0]);
  if (v[1] >= 0.0) ema(emaTI, haveTI, v[1]);
  stepTF = stepTI = 0.0;
  stepNTF = stepNTI = 0;

  if (iters < 0) return;
  const int extra = forms - (startedFresh ? 1 : 0);
  double eff = (double)iters;
  if (extra > 0 && haveTF && haveTI && emaTI > 0.0)
    eff += (double)extra * emaTF / emaTI;
  if (startedFresh)
    ema(emaNF, haveNF, eff);
  else {
    const double a = (double)((startAge < 1) ? 1 : startAge);
    rSw  = kCostForget * rSw  + 1.0;
    rSa  = kCostForget * rSa  + a;
    rSn  = kCostForget * rSn  + eff;
    rSaa = kCostForget * rSaa + a * a;
    rSan = kCostForget * rSan + a * eff;
  }
}

int
LadrunoFactorReuse::solveCurrentStep(AcceleratedNewton &a)
{
  AnalysisModel *theAnaModel = a.getAnalysisModelPtr();
  IncrementalIntegrator *theIntegrator = a.getIncrementalIntegratorPtr();
  LinearSOE *theSOE = a.getLinearSOEptr();
  ConvergenceTest *theTest = a.theTest;
  Accelerator *theAccelerator = a.theAccelerator;

  if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0) || (theTest == 0)) {
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() - setLinks() has";
    opserr << " not been called - or no ConvergenceTest has been set\n";
    return -5;
  }

  if (!announced) {
    announced = true;
    if (rank == 0)
      opserr << "KrylovNewton factor reuse (WP-188) ON: -reuseFactor " << opts.maxSteps
             << " -reuseIter " << opts.maxIter << " -reuseRate " << opts.maxRate
             << ". Steps may start on an earlier step's factor (modified Newton); every step"
                " still passes the same test, but stops at a different point inside its"
                " tolerance than full Newton -- tighten the test ~10x for result-neutral runs"
                " (ADR-188).\n";
    if (rank == 0 && opts.cost)   // WP-193
      opserr << "KrylovNewton reuse cost model (WP-193) ON: -reuseCostSafety " << opts.costSafety
             << " -reuseCostWarmup " << opts.costWarmup << " -reuseCostProbe " << opts.costProbe
             << ". Keep vs refactor is decided from measured wall times; -reuseFactor "
             << opts.maxSteps << " stays the hard age limit (ADR-193).\n";
  }

  const int numEqns = theSOE->getNumEqn();

  if (a.vAccel == 0 || a.vAccel->Size() != numEqns) {
    if (a.vAccel != 0) delete a.vAccel;
    a.vAccel = new Vector(numEqns);
  }
  if (dUstep == 0 || dUstep->Size() != numEqns) {
    if (dUstep != 0) delete dUstep;
    dUstep = new Vector(numEqns);
  }
  Vector &vAccel = *a.vAccel;

  nSteps++;

  // the step size this step was set up with (newStep has already moved the domain to t + dt)
  const double tNow = theAnaModel->getCurrentDomainTime();
  const bool dtKnown = haveLastTime;
  const double dt = dtKnown ? (tNow - lastConvergedTime) : 0.0;

  // ---- decision 1: form at step start, or start on the kept factor --------------
  int form = 1;
  if (valid && dtKnown && fSOE == theSOE && fIntegrator == theIntegrator &&
      fNumEqn == numEqns && fAge < opts.maxSteps) {
    const double scale = (fabs(dt) > fabs(fDt)) ? fabs(dt) : fabs(fDt);
    if (fabs(dt - fDt) <= kDtRelTol * scale)
      form = 0;
  }
  // ---- WP-193: the cost rule, when the hard rules above would keep ------------
  // keep iff s * n_K(a+1) <= A(S): the next step's projected cost on the kept
  // factor against the model's mean cost per step of a cycle of the S steps this
  // factor has served (stop a renewal cycle when the marginal cost reaches the
  // average; optimal for n_K non-decreasing in age, ADR-193 §3). Both sides in
  // units of t_I. Every input is rank-identical (reduced at the step end).
  int costForm = 0;
  bool costProbeStep = false;
  if (opts.cost && form == 0 && costReady()) {
    const double nk = this->projectKept(fAge + 1);
    const double avg = this->cycleMean((cycSteps > 0) ? cycSteps : 1);
    if (opts.costSafety * nk > avg) {
      if (opts.costProbe > 0 && costRefusals >= opts.costProbe)
        costProbeStep = true;     // keep once so n_K(1) is re-measured
      else
        costForm = 1;
    }
  }
  if (costForm) form = 1;

  if (rank == testVoteRank)
    form = 1;
  form = this->agree(form);

  if (opts.cost && costReady()) {   // WP-193 bookkeeping (rank-identical)
    if (costForm && form) { nCostForm++; costRefusals++; }
    else if (form == 0) {
      costRefusals = 0;
      if (costProbeStep) nCostProbe++; else nCostKeep++;
    }
  }
  const int startAge = (form == 0) ? fAge + 1 : 0;   // WP-193
  int stepForms = 0;                                  // WP-193: factorizations formed in this step
  int stepIters = 0;                                  // WP-193: iterations in this step (all attempts)
  double pendingForm = -1.0;                          // WP-193: time of a formation whose first solve is pending

  if (theAccelerator != 0)
    theAccelerator->newStep(*theSOE);

  { OPS_PROFILE_SCOPE("formUnbalance");
  if (theIntegrator->formUnbalance() < 0) {
    opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
    opserr << "the Integrator failed in formUnbalance()\n";
    this->invalidate();
    return -2;
  }
  }

  const bool startedStale = (form == 0);
  if (form) {
    const double tF0 = nowSec();   // WP-193
    { OPS_PROFILE_SCOPE("formTangent");
    if (theIntegrator->formTangent(a.tangent) < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      this->invalidate();
      return -1;
    }
    }
    burnSeconds(testFactorDelay);   // WP-193 test hook (0 unless set)
    pendingForm = nowSec() - tF0;
    stepForms++;
    a.numFactorizations++;
    nFormStart++;
    this->factorFormed(theSOE, theIntegrator, numEqns, dt);
    cycSteps = 1;                   // WP-193: a new factor's cycle
  } else {
    fAge++;
    nStepsReused++;
    cycSteps++;                     // WP-193
  }

  bool onKept = startedStale;     // the factor in use was formed at an earlier step
  int itersOnFactor = 0;

  for (int attempt = 0; attempt < 2; attempt++) {

    theTest->setEquiSolnAlgo(a);
    if (theTest->start() < 0) {
      opserr << "AcceleratedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      this->invalidate();
      return -3;
    }

    dUstep->Zero();
    int k = 1;
    int nTests = 0;
    int result = -1;
    int failCode = 0;

    do {
      const double tIt0 = nowSec();   // WP-193: one iteration, solve through test
      { OPS_PROFILE_SCOPE("linearSolve");
      if (theSOE->solve() < 0) {
        opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
        opserr << "the LinearSysOfEqn failed in solve()\n";
        failCode = -3;
        break;
      }
      }

      vAccel = theSOE->getX();

      if (theAccelerator != 0) {
        if (theAccelerator->accelerate(vAccel, *theSOE, *theIntegrator) < 0) {
          opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
          opserr << "the Accelerator failed in accelerate()\n";
          failCode = -1;
          break;
        }
      }

      // booked before the update: an update that fails part-way has still moved U
      dUstep->addVector(1.0, vAccel, 1.0);

      { OPS_PROFILE_SCOPE("update");
      if (theIntegrator->update(vAccel) < 0) {
        opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
        opserr << "the Integrator failed in update()\n";
        failCode = -4;
        break;
      }
      }
      burnSeconds(testIterDelay);   // WP-193 test hook (0 unless set)

      { OPS_PROFILE_SCOPE("formUnbalance");
      if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";
        failCode = -2;
        break;
      }
      }

      a.numIterations++;
      nIterTotal++;
      if (onKept) nIterReused++;
      itersOnFactor++;

      { OPS_PROFILE_SCOPE("convTest"); result = theTest->test(); }
      nTests++;

      // WP-193: the iteration's wall time. The first solve after a formation
      // carries the factorization: t_F sample = formTangent + that iteration - t_I.
      {
        const double tIt = nowSec() - tIt0;
        if (pendingForm >= 0.0) {
          double tf = pendingForm + tIt - (haveTI ? emaTI : 0.0);
          if (tf < pendingForm) tf = pendingForm;
          stepTF += tf; stepNTF++;
          pendingForm = -1.0;
        } else {
          stepTI += tIt; stepNTI++;
        }
        stepIters++;
      }

      if (result == -1) {
        // ---- decision 2: replace a kept factor that converges too slowly --------
        int refresh = 0;
        if (onKept) {
          if (itersOnFactor >= opts.maxIter)
            refresh = 1;
          else if (opts.maxRate > 0.0 && nTests >= 2) {
            const Vector &norms = theTest->getNorms();
            if (nTests - 1 < norms.Size()) {
              const double n1 = norms(nTests - 1);
              const double n0 = norms(nTests - 2);
              if (n0 > 0.0 && n1 > opts.maxRate * n0)
                refresh = 2;
            }
          }
          // WP-193: in-step cost rule (ski rental): once this kept factor's
          // iterations overrun the projection n_K(a) by t_F / t_I, the overrun
          // has cost as much as a fresh factorization -- buy one.
          if (refresh == 0 && costReady()) {
            const double over = (double)itersOnFactor - this->projectKept(fAge);
            if (over > 0.0 && opts.costSafety * emaTI * over >= emaTF)
              refresh = 3;
          }
          refresh = this->agree(refresh);
        }

        if (refresh != 0) {
          const double tF0 = nowSec();   // WP-193
          { OPS_PROFILE_SCOPE("formTangent");
          if (theIntegrator->formTangent(a.tangent) < 0) {
            opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
            opserr << "the Integrator failed in formTangent()\n";
            failCode = -1;
            break;
          }
          }
          burnSeconds(testFactorDelay);   // WP-193 test hook
          pendingForm = nowSec() - tF0;
          stepForms++;
          a.numFactorizations++;
          if (refresh == 2) nFormRate++; else nFormIter++;   // WP-193: 3 (cost) books as iter
          if (refresh == 3) nCostIter++;
          if (theAccelerator != 0)
            theAccelerator->newStep(*theSOE);     // the subspace belongs to the old factor
          this->factorFormed(theSOE, theIntegrator, numEqns, dt);
          cycSteps = 1;     // WP-193
          onKept = false;
          itersOnFactor = 0;
        } else if (theAccelerator != 0) {
          // the vanilla in-step refresh (KrylovNewton -maxDim)
          const double tF0 = nowSec();   // WP-193
          int ret = theAccelerator->updateTangent(*theIntegrator);
          if (ret < 0) {
            opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
            opserr << "the Accelerator failed in updateTangent()\n";
            failCode = -1;
            break;
          }
          if (ret > 0) {
            burnSeconds(testFactorDelay);   // WP-193 test hook
            pendingForm = nowSec() - tF0;
            stepForms++;
            a.numFactorizations++;
            nFormKrylov++;
            this->factorFormed(theSOE, theIntegrator, numEqns, dt);
            cycSteps = 1;   // WP-193
            onKept = false;
            itersOnFactor = 0;
          }
        }
      }

      a.record(k++);

    } while (result == -1);

    if (failCode == 0 && result >= 0) {
      this->costStepEnd(!startedStale, startAge, stepIters, stepForms);   // WP-193 (no-op when off)
      haveLastTime = true;
      lastConvergedTime = tNow;
      if (opts.report > 0 && (nSteps % opts.report) == 0)
        this->printSummary(opserr, "running");
      return result;
    }

    int code = failCode;
    if (code == 0)
      code = (result == -2) ? -3 : result;

    // ---- decision 3: a step that started on a kept factor is retried once --------
    int retry = (attempt == 0 && startedStale) ? 1 : 0;
    retry = this->agree(retry);
    if (retry == 0) {
      if (attempt == 1) {
        nRetryFailed++;
        if (rank == 0)
          opserr << "KrylovNewton factor reuse (WP-188): the retry with a fresh tangent failed too"
                    " (t = " << tNow << ")\n";
      }
      if (result == -2 && failCode == 0 && attempt == 0) {
        opserr << "AcceleratedNewton::solveCurrentStep() -";
        opserr << "The ConvergenceTest object failed in test()\n";
      }
      this->costStepEnd(!startedStale, startAge, -1, stepForms);   // WP-193: times only
      this->invalidate();
      return code;
    }

    nRetries++;
    if (rank == 0)
      opserr << "KrylovNewton factor reuse (WP-188): step at t = " << tNow
             << " failed on a kept factor (code " << code
             << "); undoing it and retrying with a fresh tangent\n";

    // back to the predictor: every integrator whose update is linear in dU
    // (all transient integrators, LoadControl) is exactly invertible; for the
    // others (path-following controls) the retry continues from the last iterate.
    if (dynamic_cast<TransientIntegrator *>(theIntegrator) != 0 ||
        dynamic_cast<LoadControl *>(theIntegrator) != 0) {
      (*dUstep) *= -1.0;
      { OPS_PROFILE_SCOPE("update");
      if (theIntegrator->update(*dUstep) < 0) {
        opserr << "WARNING KrylovNewton factor reuse: could not undo the failed attempt\n";
        this->invalidate();
        return code;
      }
      }
    }

    { OPS_PROFILE_SCOPE("formUnbalance");
    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";
      this->invalidate();
      return -2;
    }
    }
    const double tF0 = nowSec();   // WP-193
    { OPS_PROFILE_SCOPE("formTangent");
    if (theIntegrator->formTangent(a.tangent) < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      this->invalidate();
      return -1;
    }
    }
    burnSeconds(testFactorDelay);   // WP-193 test hook
    pendingForm = nowSec() - tF0;
    stepForms++;
    a.numFactorizations++;
    nFormRetry++;
    if (theAccelerator != 0)
      theAccelerator->newStep(*theSOE);
    this->factorFormed(theSOE, theIntegrator, numEqns, dt);
    cycSteps = 1;     // WP-193
    onKept = false;
    itersOnFactor = 0;
  }

  // not reached
  this->invalidate();
  return -3;
}
