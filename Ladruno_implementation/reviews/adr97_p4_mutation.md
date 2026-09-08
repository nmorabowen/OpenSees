# ADR-97 P4 (wp/97e) — mutation gate record

Verifies that the P4 test suite actually depends on the fix, not just on the
binary existing.

## Mutation applied

On a scratch build (HEAD `7e93e4381`, not committed), reverted
`compute_numerical_tangent_firstorder`/`compute_numerical_tangent_secondorder`
back to their PRE-P4 bodies: differentiate `compute_local_stress()` (the
simplified single-shot elastic-predictor / one-step plastic-corrector map
that neither `Backward_Euler` nor `Closest_Point` ever calls) instead of
calling `numerical_tangent_of_committed_map()`. `numerical_tangent_of_
committed_map()` itself, the recursion guard, and the
`ComputeTangentStiffness()` dispatch guards were left untouched — only the
two thin-wrapper bodies were swapped back to the old FD-of-a-third-map
implementation, exactly reversing the P4 re-point at its two call sites.

Built via the same `Ladruno_scripts\build.bat OpenSees OpenSeesPy` /
WMI-launch recipe as the real P4 build. Binary confirmed built from the
mutated source (log: `build_mutation.log`, `dist/bin/opensees.pyd` mtime
past launch).

## Result: RED, as expected

`tests/test_adr97_p4_numalg.py`:

```
Backward_Euler / Numerical_Algorithmic_SecondOrder   rel_err(max) = 4.573508e-02  (gate <= 1e-6)   FAILED
Backward_Euler / Numerical_Algorithmic_FirstOrder    rel_err(max) = 4.573512e-02  (gate <= 1e-6)   FAILED
Closest_Point Algorithmic-vs-NumAlg2nd assembled-K   rel diff     = 4.573506e-02  (gate <= 1e-6)   FAILED
two-cube BE+Numerical_Algorithmic_SecondOrder cost                                                 FAILED (cost no longer collapses to 16)
Closest_Point starved-VM refusal                                                                    passed (untouched by this mutation)
Backward_Euler starved-MC refusal                                                                    passed (untouched by this mutation)
```
4 failed, 2 passed (the two refusal-propagation tests are orthogonal to which
map is differentiated, so they correctly stay green under this mutation).

`tests/test_adr94_hlist_numerics.py::test_H6_no_tangent_option_reproduces_the_consistent_tangent`:
FAILED (both numerical options' `errs` returned to the historical ~4.6%
range instead of ~2-3e-8, tripping the new `< 1e-6` assertions).

The measured mutated error, **4.574e-02**, lands almost exactly on the
pre-P4 measured value recorded in the ADR-97 doc / H6's own docstring table
(**4.6%**, wp/94c) — the mutation reproduces the ORIGINAL defect precisely,
not just "some" regression, which is strong confirmation that the P4 test
suite is measuring the right thing.

## Revert and re-verify GREEN

`git checkout -- SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`
(clean revert to HEAD `7e93e4381`, confirmed via `git status`/`git diff`),
rebuilt (`build_p4_revert.log`), re-ran the full P4 suite: all green again
(see the main P4 report for the exact counts). `ladrunoBuild()` re-verified
== `7e93e4381` after the revert-rebuild.

## Conclusion

The P4 test suite is mutation-adequate for its core claim: undoing the
re-point (reverting to FD-of-`compute_local_stress`) is caught by 5 of the 7
tests that touch `Numerical_Algorithmic_*` behavior (4 in
`test_adr97_p4_numalg.py` + H6), at the exact historical magnitude of the
original defect. The 2 tests that stay green under this mutation
(refusal-propagation) are correctly testing an orthogonal property (does a
refusal from inside the helper propagate) that this particular mutation does
not touch.
