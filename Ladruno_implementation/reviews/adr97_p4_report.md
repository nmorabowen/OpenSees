# ADR-97 P4 (wp/97e-numalg-repoint) — report

PR [#829](https://github.com/nmorabowen/OpenSees/pull/829), branch
`wp/97e-numalg-repoint` cut from `wp/97f-explicit-gate` (`c10c8dca7`).
Build `7e93e4381` (final; the measured numbers below were taken on
`c24cda99c`, an ancestor with SRC-identical `ASDPlasticMaterial3D.h`, and
re-verified unchanged after the test-threshold-tightening commit).

## 1. What changed

`tangent_type Numerical_Algorithmic_FirstOrder/SecondOrder` used to
finite-difference `compute_local_stress()` — a simplified single-shot
elastic-predictor / one-step plastic-corrector that evaluates `n`, `m`, `H`
ONCE at the yield-crossing intersection and takes one closed-form `dLambda`
correction, with no resemblance to `Backward_Euler`'s cutting-plane Newton
loop or `Closest_Point`'s coupled implicit solve, and which **neither real
integrator ever calls** (confirmed by `grep -n compute_local_stress` — 7
hits, all inside its own definition or the two functions being rewritten).
This is ADR-94 M3 / ADR-84 P4's "third map" finding: the FD of this function
disagreed with the true consistent tangent by 31% (pre-wp/94c) / 4.6%
(post-wp/94c von Mises correction), not because the FD stencil was bad, but
because it differentiated a map nothing commits.

**The fix:** a new private helper,
`numerical_tangent_of_committed_map(strain_incr, tangent_matrix, second_order, ...)`,
re-enters `this->setTrialStrainIncr()` itself — the SAME dispatch a real host
element drives — for each of the 6 (first-order, forward FD) or 12
(second-order, central FD) perturbed strain increments. Per call:

1. **Snapshot** every member the configured integrator (`Backward_Euler` or
   `Closest_Point`) mutates: `TrialStrain`, `TrialStress`,
   `TrialPlastic_Strain`, the whole `iv_storage` struct, `cp_last_iterations`,
   the four `mutable` scratch buffers (`dsigma`, `depsilon_elpl`,
   `intersection_stress`, `intersection_strain`), and `Stiffness`.
2. Raise `suppress_numerical_tangent` (new per-instance, NSDMI-initialized
   member) — **the recursion guard**: both integrators call
   `ComputeTangentStiffness()` at the tail of every successful commit, so an
   unguarded perturbed sub-call would try to compute ANOTHER numerical
   tangent of its own perturbed state, unbounded.
3. For each perturbation, call `setTrialStrainIncr()`. A nonzero return code
   aborts the loop immediately, restores the snapshot, clears the
   suppression flag, and propagates that SAME code outward — **a
   partial/short tangent is never assembled**. Otherwise read `TrialStress`
   and restore the snapshot before the next perturbation.
4. Assemble the FD into a **local** matrix, never `this->Stiffness`
   column-by-column while perturbing: several `Backward_Euler` early-return
   branches (the Drucker-Prager apex switch, the `special_return` hook
   switch, the "PLASTIC INCONSISTENCY" elastic-fallback) assign `Stiffness`
   UNCONDITIONALLY, bypassing the tangent-type dispatch — a perturbed
   sub-call landing there would clobber already-written columns if the
   accumulator were the real `Stiffness` by reference.
5. Copy the local matrix into `tangent_matrix` in ONE assignment after the
   loop, after a final restore.

`compute_numerical_tangent_firstorder`/`secondorder` become thin wrappers
over this helper (kept for source compatibility). The two
`ComputeTangentStiffness()` `Numerical_Algorithmic_*` dispatch branches are
guarded by `!suppress_numerical_tangent`.

**`compute_local_stress()` fate:** zero remaining callers (grep-confirmed).
Kept, per ADR-97 D2, re-tagged with a "NOT A MAP, dead code after the P4
re-point" banner rather than deleted, for provenance.

**Four non-smooth switch sites left untouched** (out of scope, per the
pre-scope's own recommendation): `Backward_Euler`'s DP-apex switch (~4568 in
the pre-scope's line numbering) and its `special_return` hook switch (~4648),
`Closest_Point`'s own apex switch (~2743), and `cp_apply_tangent_policy`
(~2848). All four already hard-route `Numerical_Algorithmic_*` to the
analytical `Stiffness` before ever reaching `ComputeTangentStiffness()`, so
they are structurally unaffected by anything above — and a real FD across a
non-smooth apex/corner blends two different active sets, which is not
obviously meaningful.

## 2. Measured

Build `c24cda99c` (re-verified unchanged on `7e93e4381`):

| Check | Pre-P4 | Post-P4 (measured) | Gate |
|---|---|---|---|
| `Backward_Euler` + `Numerical_Algorithmic_SecondOrder` vs `fd_check` | (n/a, new gate) | **2.055e-08** | `<= 1e-6` |
| `Backward_Euler` + `Numerical_Algorithmic_FirstOrder` vs `fd_check` | (n/a, new gate) | **3.394e-08** | `<= 1e-6` |
| `Closest_Point` `Algorithmic` vs `Numerical_Algorithmic_SecondOrder` | (n/a, new gate) | **1.209e-08** | `<= 1e-6` |
| Two-cube `testIter` sum, `Backward_Euler`+`Numerical_Algorithmic_SecondOrder` | `Continuum`: **113** (6+41+36+30) | **16** (4+4+4+4) | matches `Closest_Point`+`Algorithmic`'s recorded 16 |
| H6 `errs["Numerical_Algorithmic_SecondOrder"]` | 4.6% (wp/94c) | **2.06e-8** (probe); pinned `< 1e-6` | flipped direction |
| H6 `errs["Numerical_Algorithmic_FirstOrder"]` | 4.6% (wp/94c) | **3.39e-8** (probe); pinned `< 1e-6` | flipped direction |
| Component pin, VonMises `Continuum` vs `Numerical_Algorithmic_FirstOrder` | 0.0197 | **0.2196** | re-pinned hard floor `> 0.1` |

The single-step uniaxial `fd_check` rig's converged strain increment keeps
both forward- and central-difference stencil truncation error far below the
task brief's own stated expectations (`~1e-3` forward / `<=1e-6` central) —
both landed at the same ~1e-8 order. The VonMises component pin moved in the
OPPOSITE direction from a naive "the fix makes things agree more" intuition:
`Continuum` is still only the `dLambda -> 0` limit of the consistent tangent
(ADR-94 M3/H6), and on this load path (`nsteps=20`, non-infinitesimal
`dLambda`) `Numerical_Algorithmic_FirstOrder` now tracks `Backward_Euler`'s
REAL multi-iteration cutting-plane answer instead of the old third map's
coincidentally-smaller disagreement with `Continuum` — the gap got bigger,
not smaller, and that is correct.

## 3. Gate 4 (byte-identity) — ONE deliberate baseline regeneration

`cube/vm/BE/Numerical_Algorithmic_FirstOrder/plastic` moved by up to
5.428e-09 absolute in 40 of 60 committed-stress components. Root-caused, not
assumed: that deck's rig (`test_adr94_hlist_numerics._cube_build`) is
LOAD-CONTROLLED with free DOFs, and `Backward_Euler` calls
`ComputeTangentStiffness()` at the tail of EVERY `setTrialStrainIncr()` —
every outer-Newton trial, not just the converged one. A materially different
returned tangent (now a real FD of `Backward_Euler` itself, vs. the old FD
of an unrelated map) changes the outer Newton's convergence PATH, and a
finite `NormDispIncr` tolerance (`1e-12` here) means a different path lands
at a different point within that tolerance ball — at the same order of
magnitude as the file's own documented cross-platform compiler-noise floor
(5.8e-09, MSVC vs GCC/libm).

Two independent confirmations this is expected, not a regression:

1. The SAME deck's `.../elastic` leg is bit-identical (0 of 60 components
   changed): an elastic `Backward_Euler` trial is exactly linear, so old and
   new FD reproduce the exact same analytical `E` regardless of which map
   is differentiated — no path-dependence to expose.
2. All 22 OTHER decks in the 23-deck baseline — every other `cube/vm/BE/*`
   tangent-type combination, every fully-prescribed tet deck, all four
   explicit-integrator decks — regenerated bit-identical against the P4
   binary.

The baseline was regenerated (verified diff confined to exactly this one
deck's `stress`/`strain` arrays) and the root cause documented at length in
`tests/test_adr97_p4_inertness.py`'s module docstring and
`Ladruno_implementation/LEDGER_quirks.md`.

## 4. Refusal propagation

The pre-scope's ideal reproducer — a perturbation ALONE starving while the
primary call converges — was found impractical to engineer deterministically.
Measured: bisecting `n_max_iterations` to the EXACT boundary a
`Backward_Euler` Mohr-Coulomb tet deck needs (42 fails, 43 succeeds, zero
slack) and running both `Secant` and `Numerical_Algorithmic_SecondOrder` at
niter 43/44/45 — ALL succeeded identically at every value; a `~1e-8`-relative
FD perturbation is far too small to change a quadratically-convergent
Newton's own iteration count in any case actually reachable by search.

Replaced with two reproducers that ARE reliable and test the property that
actually matters (a refusal is never masked into a wrong-but-plausible
tangent):

- `Backward_Euler` + plain `MohrCoulomb_YF` on the ADR-84 P2a exhaustion
  deck (`n_max_iterations 2`, `strict_convergence 1`) refuses IDENTICALLY
  under `Secant` (code `-3`) and under `Numerical_Algorithmic_SecondOrder`
  (code `-3`).
- `Closest_Point` + VM under the ADR-97 P1 gate-6 starved reproducer
  (`hiso=7000.0, niter=1`) still refuses (code `-3`) under
  `Numerical_Algorithmic_SecondOrder`.

Note: `Backward_Euler` needs `strict_convergence 1` to fail loud on a
starved deck (ADR-84 P2a); `Closest_Point` fails loud by default. Naively
reusing a `Closest_Point` starvation reproducer's kwargs for `Backward_Euler`
without adding the flag silently succeeds instead of refusing — this trap is
recorded in `LEDGER_quirks.md`.

## 5. Mutation gate

Reverted `compute_numerical_tangent_firstorder`/`secondorder` to their
pre-P4 `compute_local_stress()`-based bodies on a scratch build (HEAD
`7e93e4381`, never committed). Result: **4 of 6** `test_adr97_p4_numalg.py`
tests plus H6 turned RED, with the measured mutated error (**4.574e-02**)
landing almost exactly on the historically-recorded pre-P4 value (4.6%,
wp/94c) — strong confirmation the suite measures the right thing, not just
"some" regression. The 2 tests that stayed green (refusal propagation) are
correctly testing an orthogonal property this particular mutation does not
touch. Reverted (`git checkout`), rebuilt, re-verified GREEN and
`ladrunoBuild() == HEAD`. Full record: `reviews/adr97_p4_mutation.md`.

## 6. Test totals

Full explicit-list battery (`test_adr84*`, `test_adr94*` minus
`test_adr94_matrix.py` execution, `test_asdplastic_*`, `test_adr97*`):
**236 passed, 2 skipped, 0 failed** (both before and after the mutation
revert).

## 7. Surprises / notes for the owner

- Both `Numerical_Algorithmic_FirstOrder` and `SecondOrder` landed at the
  SAME order of magnitude (~1e-8) on the `fd_check` rig, not the
  forward-vs-central gap (`~1e-3` vs `<=1e-6`) the task brief anticipated —
  a property of this specific single-step rig's small converged strain
  increment, not a general claim about forward-difference accuracy.
- `compute_local_stress()` is confirmed fully dead (zero callers) but was
  NOT deleted, per ADR-97 D2's own instruction — flagged here in case the
  owner wants it gone for real in a follow-up cleanup pass.
- Diagnostic accumulators `GLOBAL_INT_max_iter[ASDP_TAG]` /
  `GLOBAL_DBL_max_error[ASDP_TAG]` are NOT restored across perturbations
  (scope.md's own risk note) — cosmetic only (a perturbed sub-call needing
  more iterations than the real step can bump these per-tag diagnostic
  maxima, printed by `commitState()`'s `cout` line), not read back into any
  committed-state decision.
