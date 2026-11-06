---
wp: LEGACY
title: "ManzariDafalias: m_Pmin silently REBUILDS the whole stress tensor, at sites PR-2's clamp diagnostic does not cover"
legacy_seq: 340
---
## ManzariDafalias: `m_Pmin` silently REBUILDS the whole stress tensor, at sites PR-2's clamp diagnostic does not cover

**Found 2026-08-27, ADR-86 PR-3, by measurement.** ADR-86 PR-2 made the low-`p` clamp
observable by adding a throttled warning at two sites — `ModifiedEuler` (`:1378`) and
`RungeKutta45` (`:1902`). That is not all of them, and the ones it misses are the ones that
fired.

On a confine-first deck at `p = 0.1145 kPa` (IntScheme 1, `m_Presidual` pinned at vanilla's
1.01, `-Pmin` at `LadrunoSANISAND`'s default `1e-3*P_atm = 0.101`), traced per committed step:

```
step 1   p = -0.5647   |dev| = 9.5818     <- mean stress goes NEGATIVE
step 2   p =  0.1010   |dev| = 2.40e-17   <- stress REPLACED by m_Pmin * mI1
```

**Zero warnings were printed.** The same deck at `-Pmin = 0.0101` never resets: its minimum
`|dev|` over the whole 40-step leg is 5.75e-01. End state differs by a factor: `p_end` 4.2254
vs 10.388 kPa.

The uninstrumented resets, all of which zero the deviator AND `alpha`:

| site | guard | writes | deviator |
|---|---|---|---|
| `explicit_integrator` `:1074`/`:1078` | `p_n < m_Presidual` (i.e. true `p` < 0) | `NextStress = m_Pmin * mI1` | **WIPED** |
| `Stress_Correction` `:2555`/`:2557`/`:2624` | true `p` < `m_Pmin` | `p = m_Pmin + m_Presidual`, then `NextStress = p * mI1` | **WIPED** |
| `BackwardEuler_CPPM` `:2213`/`:2229` | `p < m_Pmin`, then true `p` < 0 | `NextStress = m_Pmin * mI1` | **WIPED** |

**And the sharper half of this: the two sites that ARE instrumented are the two that
do the LESS destructive thing.** `ModifiedEuler:1407` and `RungeKutta45:1923` write
`GetDevPart(NextStress) + m_Pmin * mI1` — the deviator is **preserved**. All three
uninstrumented sites write a purely isotropic tensor and call `NextAlpha.Zero()` — the
deviator and the back-stress ratio are **destroyed**. So the warning covers the gentler
rebuild and is silent on the total one.

> **CORRECTED 2026-08-27 by adversarial review, before this entry ever shipped.** A
> first draft of this table carried a fourth row, `Stress_Correction:2608` ("Newton hit
> `maxIter`"), and the summary below said "two of at least **six**". That row is **dead
> code**: line 2608 sits at brace-depth 6 inside the `if (false)` block that opens at
> `:2559` and closes at `:2623`, so it can never execute. The only live effect of that
> whole branch is `p = m_Pmin + m_Presidual` at `:2557` feeding `NextStress = p * mI1`
> at `:2624`, which is *outside* the dead block and is already the row above. One
> logical site had been counted as two. Verified by mechanical brace-depth trace, not by
> eye. The true tally is **five live rebuilds, two of them instrumented**.

Note the second one: PR-2 repaired an analogous `p = m_Pmin` store in `ModifiedEuler` and
established by inspection that it was **dead** (nothing read `p` before the unconditional
recompute). The `Stress_Correction` twin is **live** — `p` is written and then used two
statements later as `NextStress = p * mI1`. Same shape, opposite consequence. Do not
generalise the PR-2 finding to it.

- **If you A/B anything at `p` below a few tenths of a kPa, pin `-Pmin`.** It is not a floor,
  it is a switch that can replace the answer.
- **Do not read "no clamp warning" as "no clamp".** PR-2's diagnostic covers **two of the
  five** live `m_Pmin`-triggered rebuilds — and the three it misses are the three that zero
  the deviator.
- Pinned by `test_pmin_is_behavioural_at_low_confinement`, which asserts the wipe
  categorically (`|dev| < 1e-10` at `p == m_Pmin`) rather than asserting a number.
- **Two more notes on the two diagnostics that DO exist** (adversarial review, 2026-08-27):
  their throttle counters are function-local `static int` incremented non-atomically, so if the
  element/material-level threading planned in ADR-75b ever lands they become a data race (UB, not
  just a garbled count). And their guards are deliberately ASYMMETRIC -- `ModifiedEuler` tests
  `p < m_Pmin + m_Presidual` while `RungeKutta45` tests `p < m_Pmin` -- because the two functions
  define `p` differently (RK45 never adds `m_Presidual` anywhere in its body). Each guard is
  correct for its own function. **Do not "fix" the asymmetry** without checking both `p`
  conventions first.
- **Owed:** instrument the three uninstrumented sites the way PR-2 instrumented its two
  (process-wide budget, not per instance — `getCopy` makes every Gauss point an instance).
  Deliberately out of PR-3's scope, which was three low-risk items.
