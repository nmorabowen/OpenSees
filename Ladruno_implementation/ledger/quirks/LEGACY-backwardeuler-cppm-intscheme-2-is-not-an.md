---
wp: LEGACY
title: "BackwardEuler_CPPM (IntScheme 2) is NOT an implicit return at low p, and its non-convergence NEVER propagates"
legacy_seq: 367
---
## `BackwardEuler_CPPM` (`IntScheme 2`) is NOT an implicit return at low `p`, and its non-convergence NEVER propagates

**Found 2026-09-05, ADR-92 P0; verified at source; measured on the oracle's corner path.**

Two things, both by design in the shipped code:

1. The low-`p` branch (`tr/3 < m_Pmin`, `:2234` — the one test that omits `m_Presidual`) has
   its Newton **disabled by a literal `errFlag = 0`** (`:2264-2266`, `NewtonIter2_negP`
   commented out, author's note *"tension-cutoff surface ... not working properly. Using
   explicit integrator for the time being"*). Every such step is integrated by
   `explicit_integrator`, i.e. `ModifiedEuler`, and flagged as success. On a Gauss-point path
   onto the `p_min` floor with shear kept on, **58-74 % of all scheme-2 calls went this way.**
2. Away from the floor, the retry ladder on Newton failure is: 50-step forward-Euler warm start
   -> recursive bisection (`mMaxSubStep = 10`) -> `explicit_integrator(...); errFlag = 1;`
   (`:2434-2437`). It gives up and returns an explicit answer **marked converged**. `Check()`'s
   `tr(stress) < 0` test is commented out (`:4907`).

- **Bites:** anyone choosing scheme 2 "because it is implicit" for a free-surface / low-`p`
  problem — the corner of a footing, a slope face, a retaining-wall backfill surface. There it
  costs a 19-unknown Newton per Gauss point and delivers `ModifiedEuler`. ADR-92 D3 was
  written on that assumption and reversed on this measurement.
- **Rule:** at low confinement the only integrator you are ever running is `ModifiedEuler`, so
  cap it (`-maxSubsteps`, #792 T1) rather than route around it. Related: the dead
  `NewtonIter2_negP` entry above.
