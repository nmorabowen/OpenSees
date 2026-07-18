/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoStaggeredDriver implementation — the shared core of the ADR-73 P2
// iterated fixed-stress driver. Decomposes DirectIntegrationAnalysis::analyzeStep
// so the solid can be re-solved against each fluid iterate, latches the driven
// overlays into external-stepping mode over the whole run (RAII), and drives the
// frozen overlay seam (advanceTrial/commitFluid/revertFluid). See the header +
// Ladruno_implementation/_adr73_implementation_plan.md P2 pins §2.2/§2.3.

#include "LadrunoStaggeredDriver.h"
#include "LadrunoPorousOverlay.h"

#include <DirectIntegrationAnalysis.h>
#include <AnalysisModel.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <classTags.h>
#include <OPS_Globals.h>

#include <vector>
#include <algorithm>
#include <math.h>

namespace ladruno_overlay {

// last-run telemetry (published at every driver exit; `-stats` reads it).
static LadrunoStaggeredStats g_lastStats;

const LadrunoStaggeredStats& LadrunoStaggeredLastStats(void) { return g_lastStats; }

namespace {

// RAII latch over the driven set: sets externalStepping_ (and pScale) ON at
// construction and clears it on EVERY exit path (normal return or early abort).
struct LatchGuard {
  std::vector<LadrunoPorousOverlay*>& driven_;
  LatchGuard(std::vector<LadrunoPorousOverlay*>& driven, double pScale)
   : driven_(driven) {
    for (size_t k = 0; k < driven_.size(); k++) {
      driven_[k]->setResidualPScale(pScale);
      driven_[k]->setExternalStepping(true);
    }
  }
  ~LatchGuard() {
    for (size_t k = 0; k < driven_.size(); k++)
      driven_[k]->setExternalStepping(false);
  }
};

} // anonymous namespace

int LadrunoStaggeredRun(DirectIntegrationAnalysis* dia, Domain* domain,
                        int nSteps, double dt, double tol, int maxIter,
                        double pScale, bool verbose)
{
  if (dia == 0 || domain == 0) {
    opserr << "ERROR LadrunoStaggeredAnalyze -- null transient analysis / domain\n";
    return -1;
  }
  if (nSteps < 0 || dt <= 0.0) {
    opserr << "ERROR LadrunoStaggeredAnalyze -- need nSteps >= 0 and dt > 0 (got "
           << nSteps << ", " << dt << ")\n";
    return -1;
  }
  // panel 2.D robustness-8: a degenerate pScale (0/negative) can turn the
  // residual NaN/negative and be silently accepted as converged; tol<=0 can
  // never converge; maxIter<1 can never iterate. Reject all three loudly.
  if (!(tol > 0.0) || maxIter < 1 || !(pScale > 0.0)) {
    opserr << "ERROR LadrunoStaggeredAnalyze -- need -tol > 0, -maxIter >= 1 "
              "and -pScale > 0 (got " << tol << ", " << maxIter << ", "
           << pScale << ")\n";
    return -1;
  }

  // ---- driven-set discovery (BEFORE any latch) ------------------------------
  // every LadrunoPorousOverlay with staticMode == SM_MARCH AND
  // fluidUpdate == FU_IMPLICIT. HOLD/STEADY/explicit overlays keep their normal
  // hook semantics and are NOT driven.
  std::vector<LadrunoPorousOverlay*> driven;
  {
    LoadPatternIter& lps = domain->getLoadPatterns();
    LoadPattern* lp;
    while ((lp = lps()) != 0) {
      if (lp->getClassTag() != PATTERN_TAG_LadrunoPorousOverlay) continue;
      LadrunoPorousOverlay* ov = static_cast<LadrunoPorousOverlay*>(lp);
      if (ov->staticModeCode()  == LadrunoPorousOverlay::SM_MARCH &&
          ov->fluidUpdateCode() == LadrunoPorousOverlay::FU_IMPLICIT)
        driven.push_back(ov);
    }
  }
  if (driven.empty()) {
    opserr << "ERROR LadrunoStaggeredAnalyze -- no qualifying overlay in the "
              "domain (need a LadrunoPorousOverlay with -staticMode march and "
              "-fluidUpdate implicit). Nothing to drive -- use plain `analyze` "
              "for a domain without a marched implicit overlay.\n";
    return -1;
  }

  AnalysisModel*       model = dia->getModel();
  EquiSolnAlgo*        algo  = dia->getAlgorithm();
  TransientIntegrator* integ = dia->getIntegrator();
  if (model == 0 || algo == 0 || integ == 0) {
    opserr << "ERROR LadrunoStaggeredAnalyze -- the transient analysis is not "
              "fully constructed (model/algorithm/integrator missing)\n";
    return -1;
  }

  // panel 2.D framework-1/2: the repeat-same-step cycle (revertToLastStep +
  // newStep on the SAME step) is only verified state-clean for the one-step
  // implicit family. TRBDF2 keeps a two-step history revertToLastStep does not
  // restore (silently degrades to trapezoidal per iteration); CentralDifference
  // and the explicit family inherit the base no-op revertToLastStep and
  // accumulate state across repeats (and the rate form is CD-unstable, ADR-73
  // §12 P0 pin 3). Whitelist the verified integrators; fail loudly otherwise.
  {
    int itag = integ->getClassTag();
    if (itag != INTEGRATOR_TAGS_Newmark &&
        itag != INTEGRATOR_TAGS_HHT &&
        itag != INTEGRATOR_TAGS_GeneralizedAlpha) {
      opserr << "ERROR LadrunoStaggeredAnalyze -- unsupported transient "
                "integrator (classTag " << itag << "). The iterated staggered "
                "driver re-runs the SAME step against each fluid iterate, which "
                "is verified only for Newmark, HHT and GeneralizedAlpha. "
                "Multi-step (TRBDF2/TRBDF3) and explicit (CentralDifference "
                "family) integrators do not survive the revert/re-newStep "
                "cycle -- use plain `analyze` (fs1) or a supported integrator.\n";
      return -1;
    }
  }

  // telemetry accumulators (published at every exit, success or abort)
  long   nStepsDone = 0, totalSolves = 0, sumIters = 0, maxItersSeen = 0;
  double lastResidual = 0.0, maxResidual = 0.0;

  auto publish = [&]() {
    g_lastStats.nSteps           = nStepsDone;
    g_lastStats.totalFluidSolves = totalSolves;
    g_lastStats.meanIters        = (nStepsDone > 0)
      ? (double)sumIters / (double)nStepsDone : 0.0;
    g_lastStats.maxIters         = maxItersSeen;
    g_lastStats.lastResidual     = lastResidual;
    g_lastStats.maxResidual      = maxResidual;
  };

  // ---- step 0: catch-up any pending -subcycle window (before the latch) ------
  for (size_t k = 0; k < driven.size(); k++)
    if (driven[k]->catchUpPendingWindow() < 0) {
      opserr << "LadrunoStaggeredAnalyze: catch-up fluid sync failed on overlay "
             << driven[k]->getTag() << " (before step 0)\n";
      publish();
      return -6;
    }

  // ---- latch the driven set for the whole run (RAII clears on every exit) ----
  LatchGuard guard(driven, pScale);

  for (int step = 0; step < nSteps; step++) {

    // ---- iteration 1: the base solid solve (overlay forces +Q·pTrial_, and
    //      pTrial_ == pCommitted_ at step start) --------------------------------
    if (model->analysisStep(dt) < 0) {
      opserr << "LadrunoStaggeredAnalyze: analysisStep failed at step " << step
             << " (time " << domain->getCurrentTime() << ")\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -2;
    }
    if (dia->checkDomainChange() < 0) {
      opserr << "LadrunoStaggeredAnalyze: checkDomainChange failed at step "
             << step << "\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -2;
    }
    if (integ->newStep(dt) < 0) {
      opserr << "LadrunoStaggeredAnalyze: newStep failed at step " << step
             << " nit 1 (time " << domain->getCurrentTime() << ")\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -2;
    }
    if (algo->solveCurrentStep() < 0) {
      opserr << "LadrunoStaggeredAnalyze: solveCurrentStep failed at step "
             << step << " nit 1 (time " << domain->getCurrentTime() << ")\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -3;
    }

    int    nit = 1;
    double residual = 0.0;
    {
      bool advOk = true;
      for (size_t k = 0; k < driven.size(); k++) {
        if (driven[k]->advanceTrial(true, dt) < 0) { advOk = false; break; }
        totalSolves++;
      }
      if (!advOk) {
        opserr << "LadrunoStaggeredAnalyze: fluid advance failed at step " << step
               << " nit 1\n";
        domain->revertToLastCommit();
        integ->revertToLastStep();
        for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
        publish();
        return -6;
      }
      for (size_t k = 0; k < driven.size(); k++)
        residual = std::max(residual, driven[k]->lastAdvanceRelChange());
    }
    // panel 2.D robustness-8: a non-finite residual must never satisfy the
    // `residual >= tol` exits below (NaN compares false both ways = silent
    // false-convergence). Abort loudly instead.
    if (residual != residual) {
      opserr << "LadrunoStaggeredAnalyze: non-finite fluid residual at step "
             << step << " nit 1 -- aborting\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -6;
    }

    // ---- inner fixed-point iteration ------------------------------------------
    bool aborted = false;
    int  abortCode = 0;
    while (residual >= tol && nit < maxIter) {
      domain->revertToLastCommit();
      integ->revertToLastStep();
      if (integ->newStep(dt) < 0) {
        opserr << "LadrunoStaggeredAnalyze: newStep failed at step " << step
               << " nit " << (nit + 1) << "\n";
        aborted = true; abortCode = -2; break;
      }
      if (algo->solveCurrentStep() < 0) {
        opserr << "LadrunoStaggeredAnalyze: solveCurrentStep failed at step "
               << step << " nit " << (nit + 1) << " (residual " << residual << ")\n";
        aborted = true; abortCode = -3; break;
      }
      nit++;
      bool advOk = true;
      for (size_t k = 0; k < driven.size(); k++) {
        if (driven[k]->advanceTrial(false, dt) < 0) { advOk = false; break; }
        totalSolves++;
      }
      if (!advOk) {
        opserr << "LadrunoStaggeredAnalyze: fluid advance failed at step " << step
               << " nit " << nit << "\n";
        aborted = true; abortCode = -6; break;
      }
      residual = 0.0;
      for (size_t k = 0; k < driven.size(); k++)
        residual = std::max(residual, driven[k]->lastAdvanceRelChange());
      if (residual != residual) {                 // non-finite guard (see nit 1)
        opserr << "LadrunoStaggeredAnalyze: non-finite fluid residual at step "
               << step << " nit " << nit << " -- aborting\n";
        aborted = true; abortCode = -6; break;
      }
    }
    if (aborted) {
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return abortCode;
    }

    if (residual >= tol) {   // maxIter reached unconverged (pathology per e76)
      opserr << "LadrunoStaggeredAnalyze: step " << step << " reached maxIter "
             << maxIter << " unconverged (residual " << residual << " >= tol "
             << tol << "). classic-L should never hit this on a sane problem.\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -7;
    }

    // ---- converged: final momentum resolve, then commit ----------------------
    // (§3.1 step 3): re-solve the solid once more with the converged p_K forces,
    // then integrator->commit() → Domain::commit → the overlay hook fires LATCHED
    // (commitFluid only, no advance).
    domain->revertToLastCommit();
    integ->revertToLastStep();
    if (integ->newStep(dt) < 0) {
      opserr << "LadrunoStaggeredAnalyze: newStep failed in final resolve at step "
             << step << "\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -2;
    }
    if (algo->solveCurrentStep() < 0) {
      opserr << "LadrunoStaggeredAnalyze: solveCurrentStep failed in final resolve "
                "at step " << step << "\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -3;
    }
    if (integ->commit() < 0) {
      opserr << "LadrunoStaggeredAnalyze: integrator commit failed at step "
             << step << "\n";
      domain->revertToLastCommit();
      integ->revertToLastStep();
      for (size_t k = 0; k < driven.size(); k++) driven[k]->revertFluid();
      publish();
      return -4;
    }

    // telemetry for this committed step
    nStepsDone++;
    sumIters += nit;
    if (nit > maxItersSeen) maxItersSeen = nit;
    lastResidual = residual;
    if (residual > maxResidual) maxResidual = residual;

    if (verbose)
      opserr << "LadrunoStaggeredAnalyze: step " << step << ": iters=" << nit
             << " residual=" << residual << "\n";
  }

  // ---- success: publish telemetry + flush recorders (the analyze convention) --
  publish();
  domain->flushRecorders();
  return 0;
}

} // namespace ladruno_overlay
