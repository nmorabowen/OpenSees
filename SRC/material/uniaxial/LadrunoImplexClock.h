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

// Ladruno WP-217: the IMPL-EX extrapolation clock shared by LadrunoCohesiveHinge,
// LadrunoCohesiveHingeBiaxial, ModIMKPeakOriented and Bilin.
//
// IMPL-EX (Oliver et al. 2008; Prazeres et al. 2016 Eq. 17) extrapolates the
// damage driver with the step ratio f = dt_{n+1}/dt_n. The four materials used
// to read dt_{n+1} = ops_Dt on every trial. In a static analysis ops_Dt is the
// load-factor increment (Domain::update), and under DisplacementControl or
// arc-length that increment is part of the solution: it moves with every Newton
// iterate, so the "explicit" step was not explicit (frame-ladruno-validation
// 00a, B3: errors 2.28e-2, 1.53e-2, 1.53e-2, 2.87e-3 Mc for N = 40..320).
//
//   STEP   (default) pseudo-time advances by one per committed step: f = 1.
//          Iterate-independent under every integrator. Exact for uniform
//          steps of the control variable (DisplacementControl with a constant
//          dU, LoadControl with a constant dlambda, Newmark with a constant dt).
//   PSEUDO f = dt_{n+1}/dt_n from ops_Dt, the pre-WP-217 rule, but dt_{n+1}
//          is LATCHED at the first trial of the step, so the step stays
//          explicit. For LoadControl / transient runs with variable steps
//          (cut-backs, VariableTransient). Under DisplacementControl the latch
//          is the predictor's dlambda: explicit, but not a meaningful clock.
//
// f is capped at 2 in both modes (the pre-WP-217 guard).

#ifndef LadrunoImplexClock_h
#define LadrunoImplexClock_h

#include <string.h>

extern double ops_Dt;

class LadrunoImplexClock
{
 public:
  enum Mode { STEP = 0, PSEUDO = 1 };

  LadrunoImplexClock() : mode(STEP), latched(false), dtStep(0.0) {}

  // Parser helper: "step" or "pseudo". Returns false on anything else.
  static bool parse(const char* s, int& m)
  {
    if (s == 0) return false;
    if (strcmp(s, "step") == 0)   { m = STEP;   return true; }
    if (strcmp(s, "pseudo") == 0) { m = PSEUDO; return true; }
    return false;
  }

  // Extrapolation factor for the open step. dtCommit is the dt the material
  // recorded at its last commit (commitDt()). Frozen for the whole step.
  double factor(double dtCommit)
  {
    if (mode == STEP) return 1.0;
    if (!latched) { dtStep = ops_Dt; latched = true; }
    double f = 1.0;
    if (dtStep > 0.0 && dtCommit > 0.0) {
      f = dtStep / dtCommit;
      if (f > 2.0) f = 2.0;
    }
    return f;
  }

  // True when the PSEUDO rule has no ratio to form (first step, a hold, or a
  // non-positive dt). The hinges use it to keep their first-step f = 0.
  bool noRatio(double dtCommit)
  {
    if (mode == STEP) return false;
    if (!latched) { dtStep = ops_Dt; latched = true; }
    return !(dtStep > 0.0 && dtCommit > 0.0);
  }

  // dt to record at commit, and close the step.
  double commitDt(void)
  {
    double dt = latched ? dtStep : ops_Dt;
    latched = false;
    return dt;
  }

  // A reverted step is retried with a new dt: re-latch on its first trial.
  void reopen(void) { latched = false; }

  // Serialization: packs the mode beside the implex flag in one double slot.
  // An old file (slot 0 or 1) reads back as STEP.
  static double pack(bool implex, int m)     { return (implex ? 1.0 : 0.0) + (m == PSEUDO ? 2.0 : 0.0); }
  static bool   unpackImplex(double v)       { return (((int)v) & 1) != 0; }
  static int    unpackMode(double v)         { return (((int)v) & 2) ? PSEUDO : STEP; }

  int    mode;
  bool   latched;
  double dtStep;
};

#endif
