---
wp: ADR-94
title: "Narrowing an apex region hands the flank Newton states it can answer WRONG — the cohesion-softening deviator flip (ADR-94 addendum, F8)"
legacy_seq: 444
---
## Narrowing an apex region hands the flank Newton states it can answer WRONG — the cohesion-softening deviator flip (ADR-94 addendum, F8)

A consequence of the fix above, found by adversarial review of #836 and fixed in
the same PR. Once the apex region is the exact (narrower, under dilatancy) one,
near-boundary trials that used to be apex-projected are handed to the flank
scalar Newton. That Newton's `dPhi/dlambda` carries the shipped yield function's
`df/dk = -1` term for a cohesion internal variable **that `f` itself does not
contain** (ADR-97 P0 header finding 2, pinned not fixed). With cohesion
SOFTENING the inconsistency lets it converge — `|Phi| ~ 1e-7`, `rc = 0`, no
exhaustion, no refusal — onto the WRONG ROOT: a state whose deviator points
OPPOSITE the trial deviator, i.e. the map walked through the vertex and out the
other side.

Measured (`ScalarLinearHardeningParameter = -20000`, `etabar = eta`, the ADR-95
cone, one Gauss point):

| (p - p_apex)/q | HS = 0 (correct) | HS = -20000, pre-guard |
|---|---|---|
| 3.0 | p -0.189103, q 0.199762, `s_zz-s_xx` **+0.346** | p -0.478023, q 0.328548, `s_zz-s_xx` **-0.569** |
| 4.3089 | the apex, q 1.0e-15 | p -0.838550, q 0.489254, `s_zz-s_xx` **-0.847** |

Identical at `strict_convergence` 0 and 1.

Three things this teaches:

7. **A yield-function tolerance cannot certify a return map.** Both wrong states
   sit ON the surface to 1e-7. The check that catches them is GEOMETRIC: a
   Drucker-Prager return is a non-negative radial scaling of the trial RELATIVE
   deviator `r = dev(sigma) - alpha` plus a pressure change, so
   `dot(r_ret, r_tr) < 0` is inadmissible for any parameters. That guard is now
   in `Backward_Euler`, reported through the existing `be_flank_failed` channel
   so it inherits the apex fallback and the fail-loud refusal without adding an
   exit.

   **`r`, not `dev(sigma)` — the first version of this guard got that wrong**
   (caught by review round 2), and the two are only the same statement when the
   back stress is zero. With `alpha` antiparallel to the trial deviator, the
   ABSOLUTE deviator can legitimately cross zero while `sqrt(J2(r))` stays
   positive: `alpha = (0.02, 0.02, -0.04, 0, 0, 0)`, `Ht = 0`, associated, trials
   `(q_tr, p_tr) = (0.05, 0.5103)` and `(0.02, 0.3810)`. (On those two trials the
   guard is not what refuses — the apex CLASSIFICATION gets there first; see the
   entry below. The guard would have been the next thing to be wrong, and is
   fixed in the same variable at the same time.) The general
   lesson: **a guard written from a picture of the surface must be written in the
   surface's own variables**, and for any kinematically hardening model that is
   the relative stress, never the raw one.

   Two things the guard cannot see, by construction: a flip whose returned
   relative deviator is smaller than `f_absolute_tol` (that is the scale at which
   the integrator stops distinguishing states at all, so the floor is not a
   tuning knob), and anything at a yield function that does not declare
   `yf_apex_elastic_metric`.
8. **`strict_convergence` is not the safety net people assume.** `be_exhausted`
   — one of layer (b)'s only two triggers — is computed **only** when strict is
   on, and strict is OFF by default. A failure mode that CONVERGES is invisible
   to both.
9. **Narrowing a classification is not automatically conservative.** It moves
   states from a branch that always "works" (a vertex projection cannot fail)
   onto one that can be wrong in a new way. Ask what the newly-routed states do,
   not only whether the classification is now exact.
