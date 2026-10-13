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

#include <cmath>
#include <cstdlib>
#include <cstring>

namespace {

// Relative tolerance on the step size: a factor formed for dt serves a step of dt'
// only if |dt' - dt| <= kDtRelTol * max(|dt|, |dt'|). Domain times are sums
// t_c + dt, so the recovered dt carries ~1e-16 * t / dt of round-off (3e-12 for a
// 70 s record at dt = 2.5 ms); 1e-8 is far above that and far below any real cut.
const double kDtRelTol = 1.0e-8;

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
  return 0;
}

bool isReuseFlag(const char *flag)
{
  return flag != 0 && (strcmp(flag, "-reuseFactor") == 0 || strcmp(flag, "-reuseIter") == 0 ||
                       strcmp(flag, "-reuseRate") == 0 || strcmp(flag, "-reuseReport") == 0);
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
  return true;
}

LadrunoFactorReuse::LadrunoFactorReuse(const LadrunoFactorReuseOptions &o)
  : opts(o), valid(false), fSOE(0), fIntegrator(0), fNumEqn(-1), fDt(0.0), fAge(0),
    haveLastTime(false), lastConvergedTime(0.0), dUstep(0),
    nSteps(0), nStepsReused(0), nFormStart(0), nFormIter(0), nFormRate(0), nFormKrylov(0),
    nFormRetry(0), nRetries(0), nRetryFailed(0), nIterTotal(0), nIterReused(0),
    rank(ladrunoFactorReuseRank()), announced(false), testVoteRank(-1)
{
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
  const bool wantDouble = (strcmp(flag, "-reuseRate") == 0);
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
  if (strcmp(flag, "-reuseRate") == 0) {
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
  if (rank == testVoteRank)
    form = 1;
  form = this->agree(form);

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
    { OPS_PROFILE_SCOPE("formTangent");
    if (theIntegrator->formTangent(a.tangent) < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      this->invalidate();
      return -1;
    }
    }
    a.numFactorizations++;
    nFormStart++;
    this->factorFormed(theSOE, theIntegrator, numEqns, dt);
  } else {
    fAge++;
    nStepsReused++;
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
          refresh = this->agree(refresh);
        }

        if (refresh != 0) {
          { OPS_PROFILE_SCOPE("formTangent");
          if (theIntegrator->formTangent(a.tangent) < 0) {
            opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
            opserr << "the Integrator failed in formTangent()\n";
            failCode = -1;
            break;
          }
          }
          a.numFactorizations++;
          if (refresh == 1) nFormIter++; else nFormRate++;
          if (theAccelerator != 0)
            theAccelerator->newStep(*theSOE);     // the subspace belongs to the old factor
          this->factorFormed(theSOE, theIntegrator, numEqns, dt);
          onKept = false;
          itersOnFactor = 0;
        } else if (theAccelerator != 0) {
          // the vanilla in-step refresh (KrylovNewton -maxDim)
          int ret = theAccelerator->updateTangent(*theIntegrator);
          if (ret < 0) {
            opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
            opserr << "the Accelerator failed in updateTangent()\n";
            failCode = -1;
            break;
          }
          if (ret > 0) {
            a.numFactorizations++;
            nFormKrylov++;
            this->factorFormed(theSOE, theIntegrator, numEqns, dt);
            onKept = false;
            itersOnFactor = 0;
          }
        }
      }

      a.record(k++);

    } while (result == -1);

    if (failCode == 0 && result >= 0) {
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
    { OPS_PROFILE_SCOPE("formTangent");
    if (theIntegrator->formTangent(a.tangent) < 0) {
      opserr << "WARNING AcceleratedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      this->invalidate();
      return -1;
    }
    }
    a.numFactorizations++;
    nFormRetry++;
    if (theAccelerator != 0)
      theAccelerator->newStep(*theSOE);
    this->factorFormed(theSOE, theIntegrator, numEqns, dt);
    onKept = false;
    itersOnFactor = 0;
  }

  // not reached
  this->invalidate();
  return -3;
}
