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

// Ladruno WP-188 -- factorization reuse across time steps for AcceleratedNewton
// (`algorithm KrylovNewton ... -reuseFactor m ...`). ADR:
// Ladruno_implementation/188_ladruno_factor_reuse_adr.md.
//
// WHAT IT IS. A modified-Newton strategy with lazy refactorization: a step may
// START on the LU of a tangent formed at an earlier step (no formTangent, the SOE
// stays `factored`, MUMPS runs job 3 only), Krylov-accelerated as KrylovNewton
// already is inside a step. The residual is always the exact current residual and
// the step exits only when the user's ConvergenceTest passes, so the equations are
// unchanged; only the iteration path changes, and with it the point inside the
// tolerance ball at which the step stops (ADR §4).
//
// TRIGGERS (all evaluated identically on every rank and then max-reduced, so the
// ranks of a distributed SOE never disagree on job 5 vs job 3):
//   form at step start  no valid factor, or it already served `maxSteps` steps, or
//                       the step size changed, or the SOE / integrator / size /
//                       domain changed since it was formed;
//   -reuseIter k        a kept factor has run k iterations in this step without
//                       convergence -> form the current tangent, restart Krylov;
//   -reuseRate r        a kept factor's test-norm ratio n_i / n_{i-1} > r -> same;
//   failure             a step that ran on a kept factor fails (test -2, solve or
//                       update error) -> undo the step's increments, form the
//                       tangent at the predictor and retry the step ONCE as plain
//                       KrylovNewton before reporting failure.
//
// OFF BY DEFAULT. AcceleratedNewton holds a null pointer unless the option is given,
// and its vanilla solveCurrentStep runs unchanged (byte-identical, WP-188 gate).

#ifndef LadrunoFactorReuse_h
#define LadrunoFactorReuse_h

class AcceleratedNewton;
class LinearSOE;
class IncrementalIntegrator;
class OPS_Stream;
class Vector;

struct LadrunoFactorReuseOptions {
  int    maxSteps = 0;     // -reuseFactor m: steps a factor may serve after the one that formed it (0 = off)
  int    maxIter  = 10;    // -reuseIter k
  double maxRate  = 1.0;   // -reuseRate r (<= 0 disables the rate trigger)
  int    report   = 0;     // -reuseReport n: print the running statistics every n steps (0 = only at the end)
  bool   sawFactor = false; // -reuseFactor given
  bool   sawOther  = false; // any other -reuse* option given
  // fail closed: a -reuse* tuning option without -reuseFactor is a deck error
  bool consistent(void) const;
};

class LadrunoFactorReuse
{
 public:
  explicit LadrunoFactorReuse(const LadrunoFactorReuseOptions &opts);
  ~LadrunoFactorReuse();

  // Parsers shared by the two algorithm-command paths. Both return 1 if `flag` was a
  // reuse option and was consumed with its value, 0 if it is not a reuse option, -1
  // on a missing/invalid value (the caller refuses the command).
  static int parseOPS(const char *flag, LadrunoFactorReuseOptions &opts);             // OPS_Get* stream
  static int parseArgv(int argc, const char **argv, int &i, LadrunoFactorReuseOptions &opts);  // classic Tcl

  // The step driver (replaces AcceleratedNewton::solveCurrentStep when enabled).
  int solveCurrentStep(AcceleratedNewton &algo);

  void invalidate(void);            // domainChanged() and every failure path
  void printSummary(OPS_Stream &s, const char *when, bool allRanks = false);

  const LadrunoFactorReuseOptions &options(void) const { return opts; }

 private:
  void factorFormed(LinearSOE *soe, IncrementalIntegrator *integ, int n, double dt);
  int  agree(int flag);             // MAX over the ranks (identity in serial builds)

  LadrunoFactorReuseOptions opts;

  // the factor currently held by the SOE
  bool   valid;
  LinearSOE *fSOE;
  IncrementalIntegrator *fIntegrator;
  int    fNumEqn;
  double fDt;
  int    fAge;                      // steps started on it since it was formed

  // step bookkeeping
  bool   haveLastTime;
  double lastConvergedTime;
  Vector *dUstep;                   // sum of the increments applied in the current attempt

  // statistics
  long nSteps, nStepsReused, nFormStart, nFormIter, nFormRate, nFormKrylov, nFormRetry;
  long nRetries, nRetryFailed, nIterTotal, nIterReused;
  int  rank;
  bool announced;
  int  testVoteRank;   // test hook, see the constructor
};

#endif
