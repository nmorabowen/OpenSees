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
// LadrunoStaggeredDriver — the shared core of the ADR-73 P2 iterated fixed-stress
// driver `LadrunoStaggeredAnalyze` (registered in BOTH interpreters). Given the
// active DirectIntegrationAnalysis and its Domain, it drives every qualifying
// LadrunoPorousOverlay (staticMode == SM_MARCH, fluidUpdate == FU_IMPLICIT)
// through the toy `qs_staggered(mode="fs_iter")` recipe: per step it re-solves
// the solid with the current fluid iterate's forces (+Q·pTrial_) and advances the
// fluid (fixed-point reference on repeats — L cancels at convergence = monolithic
// BE), iterating until the max overlay relative-change residual < tol, then does a
// final momentum resolve + integrator commit. The overlay seam
// (advanceTrial/commitFluid/revertFluid + the external-stepping latch) is FROZEN;
// this driver only drives it. See Ladruno_implementation/_adr73_implementation_plan.md
// "P2 — iterated driver + stage transport" and ADR §12 P0/P1.
//
// Scope limits (panel 2.D framework-1/2/3/4): the repeat-same-step cycle is
// whitelisted to the verified one-step implicit integrators (Newmark, HHT,
// GeneralizedAlpha) — TRBDF2's two-step history and the explicit family's no-op
// revertToLastStep do not survive it (loud fatal). The driver's decomposed step
// omits analyzeStep's _RELIABILITY sensitivity pass and the numSubLevels dt
// subdivision — sensitivity users and sub-level configs get plain `analyze`
// semantics only there.

#ifndef LadrunoStaggeredDriver_h
#define LadrunoStaggeredDriver_h

class DirectIntegrationAnalysis;
class Domain;

namespace ladruno_overlay {

// telemetry of the LAST driver run (the `-stats` payload; zeros before any run).
struct LadrunoStaggeredStats {
  long   nSteps;             // committed driver steps
  long   totalFluidSolves;   // total advanceTrial calls (Σ iters over all steps)
  double meanIters;          // mean iterations per committed step
  long   maxIters;           // worst-step iteration count
  double lastResidual;       // residual at the last committed step's convergence
  double maxResidual;        // worst committed-step converged residual
  LadrunoStaggeredStats()
   : nSteps(0), totalFluidSolves(0), meanIters(0.0),
     maxIters(0), lastResidual(0.0), maxResidual(0.0) {}
};

// read the last run's telemetry (file-scope instance in the .cpp).
const LadrunoStaggeredStats& LadrunoStaggeredLastStats(void);

// The iterated fixed-stress driver. Returns 0 on success, negative on abort
// (the `analyze` convention): −2 newStep / domain-change, −3 solve, −4 commit,
// −6 fluid advance, −7 maxIter unconverged, −1 misuse (empty driven set / bad
// args). All-or-nothing per step: on any failure the step is reverted
// (domain->revertToLastCommit + integrator->revertToLastStep + revertFluid) and
// prior committed steps stand. The latch is cleared on EVERY exit path (RAII).
int LadrunoStaggeredRun(DirectIntegrationAnalysis* dia, Domain* domain,
                        int nSteps, double dt, double tol, int maxIter,
                        double pScale, bool verbose);

} // namespace ladruno_overlay

#endif
