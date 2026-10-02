---
name: ladruno-new-material
description: >
  Checklist for adding or changing a MATERIAL in the Ladruno OpenSees fork (SRC/material/,
  nD or uniaxial: LadrunoSANISAND, LadrunoJ2, LadrunoConcrete3D, Staged*/LogStrain wrappers,
  ASDPlasticMaterial3D changes, ManzariDafalias-family fixes). Use before writing or modifying
  a return map / substepping scheme, getTangent/getInitialTangent, IMPL-EX, getCopy,
  commit/revert, any static or process-wide state, a Tcl/Python material parser, or material
  tests. Every item points to the LEDGER_quirks entry that explains it.
---

# New or changed material — checklist

Read this before adding or changing a material under `SRC/material/`. Each item names the
`Ladruno_implementation/LEDGER_quirks.md` heading to grep for; read that entry when the item
applies. Items marked **[lint]** are enforced by `python ci/check_quirk_patterns.py`.

## Registration (a new material)

- [ ] classTag, manifest row + test, header stamp, ledger row, banner line: same as an element
      (see `AGENTS.md` and the `ladruno-new-element` guide).
- [ ] The Tcl `nDMaterial` command is a hand-written `strcmp` ladder, separate from the Python
      path: register in both. Quirks: "The Tcl `nDMaterial` command is a hand-written strcmp ladder".

## Process-wide state — the recurring one

- [ ] **[lint]** No process-wide mutable state unless it is reset on `wipe`. `wipe()` does not
      recreate the Domain; reset in `Domain::clearAll()`, or mark the declaration
      `// ladruno-lint: wipe-ok <reason>`. Quirks: "`wipe()` does NOT recreate the Domain",
      "HardeningLawStorage is a process-global", "`ManzariDafalias::mElastFlag` is STATIC",
      "`getCommitTag()` is a GLOBAL monotonic counter". Open instance: PR #841 (SANISAND).
- [ ] `getCopy` makes every Gauss point an instance: a "budget" or latch is per instance unless
      you deliberately make it process-wide. A WORK-DONE latch set in place must be copied by `getCopy`,
      and a loud failure must not latch. A REFUSAL latch (bound to one integration point) is NOT
      copied: a clone starts clear, and only sendSelf keeps it. Quirks: "getCopy must PROPAGATE the
      latch", "getCopy of a LATCHED LadrunoNorSand".

- [ ] Vanilla materials that keep per-material data in static `...x[matN]` arrays (PDMY/PIMY
      family): a new field needs the constructor store, the `matCount%20` copy loop AND the
      `recvSelf` reallocation, plus a test that creates >20 materials after the one under test.
      Quirks: "hard-coded its critical-state line".

## Return map, tangent, substepping

- [ ] A non-converged return map must FAIL (return < 0), never commit `f > 0` as success, and
      the refusal must reach `analyze`: several vanilla elements swallow it. Quirks:
      "`Backward_Euler` ACCEPTS a non-converged return map", "swallow the material's refusal".
- [ ] The work of one `setTrialStrain` must be BOUNDED, and hitting the bound must REFUSE
      (return < 0), never grind on and never force-accept. Newton, line-search and Krylov trial
      iterates can be far off the solution path (|Δu| ~ 1e4 has been observed), so any cost that
      scales with |Δε| (substep counts, recursive halvings, local-Newton restarts) needs a cap,
      and exceeding it refuses so the global step is cut. Three incidents, one rule: PDMY
      `setSubStrainRate` asked for ~1e9 substeps per point, a silent hang (WP-135, #874);
      SANISAND `BackwardEuler_CPPM` recursed up to 2^9 halvings, 12–134 s per failing step
      (WP-130, #868); SANISAND `ModifiedEuler` force-accepted failed substeps at `dT_min`
      (WP-127 finding C, SAS-ME fix WP-129, #871). Test it: feed one wild trial increment
      and assert the call returns < 0 within a wall-clock bound.
- [ ] Implement `getInitialTangent()` honestly: the base default returns `getTangent()`, so
      `-initial` silently becomes full Newton. Quirks: "`NDMaterial::getInitialTangent()` DEFAULTS".
- [ ] Substep schemes need error control and yield-drift correction, and must honour the
      tolerance passed in. Quirks: "`IntScheme` 3 (RungeKutta4) and 5 (ForwardEuler) have no
      error control", "IntScheme 1 (ModifiedEuler) IGNORES the `TolR`".
- [ ] The substep error must measure EVERY evolved internal variable (back-stress, fabric,
      ...), not only the stress: a stress-only test accepts an O(1) back-stress jump when both
      stages are elastic in stress. Quirks: "substep error is STRESS-ONLY" (WP-128).
- [ ] A substep scheme must not ACCEPT a substep that failed its error test at the minimum
      step, or return early at `T < 1`, without saying so: count it (WP-127 `substepStats`).
      Quirks: "ACCEPTS a substep that FAILED its error test".
- [ ] A substep count sized from the increment (`|Δε|/h`) must be capped, and past the cap the
      trial refused: a Newton iterate can be ~1e4 and ask for ~1e9 substeps (an apparent hang).
      Quirks: "one wild Newton iterate makes `setSubStrainRate()` ask for".
- [ ] IMPL-EX in a static analysis: `ops_Dt` is pseudo-time and erratic; guard the
      extrapolation factor. Quirks: "IMPL-EX in a STATIC analysis".
- [ ] `revertToStart()` must not reset calibrated constants mid-analysis. Quirks:
      "`ManzariDafalias::revertToStart()` silently restores".
- [ ] Constants that multiply a stress are dimensional: make them unit-consistent or document
      the units. Quirks: "`D_factor` dilatancy sigmoid is DIMENSIONAL".

- [ ] A fork block appended to a base `sendSelf` under the same dbTag and commitTag must NOT have the length of any vector the base sends: FE_Datastore keys vectors by size, and a same-size block overwrites the base state. `static_assert` the size. Quirks: "FE_Datastore keys a sent Vector by its SIZE".

## NaN and silent success

- [ ] Never check divergence with `pNorm(0)` (NaN-blind); use `std::isfinite`. Quirks:
      "`Vector::pNorm(0)` is NaN-BLIND". Eigen `*= 0` keeps NaN garbage:
      "with `*= 0` keep NaN heap garbage".
- [ ] `analyze()` returning 0 does not mean the numbers are finite. Quirks: "`analyze()` returns
      rc=0 on a NaN-poisoned system".

## Finite strain and regularization

- [ ] Don't lift a damage/softening material with the generic `LogStrainNDMaterial` wrapper.
      Quirks: "The generic LogStrainNDMaterial wrapper is UNSOUND".
- [ ] A wrapper around a refusing material (LogStrain, Staged*, InitDefGrad) must forward `LADRUNO_MATERIAL_REFUSED` from the trial WITHOUT touching its staged state, and commit the inner first without advancing its own state on a refusing commit. Quirks: "`LogStrainNDMaterial::setTrialF` DROPPED the inner's refusal".
- [ ] Crack-band `lch` comes from the element through a global. Quirks: "Crack-band materials
      read element size".

## Tests

- [ ] Verifying a tangent needs free equations: a fully prescribed material-point driver cannot
      see a wrong tangent. `setStrain` commits, so a central FD tangent is invalid for
      plasticity. Uniaxial single-element tests leave shear terms at zero; use `NDTest` (3D only).
      Quirks: "has ZERO free equations", "`setStrain` (testUniaxialMaterial) COMMITS",
      "Uniaxial single-element tests leave the shear".
- [ ] Tune test paths against plastic response, not elastic estimates. Quirks:
      "must be tuned against PLASTIC response".
- [ ] Break each new gate on purpose once. Quirks: "A test can be GREEN because of the very bug".
- [ ] A tangent READ as "algorithmic" is not verified: compare it with a finite difference of the
      return map (the replay facility gives one at a real state). Quirks: "`TanType 2` tangent is
      MINUS the derivative of its own return map" (WP-130).
- [ ] Kernel-vs-oracle parity: the tangent at an exact Lode corner or near the vertex is
      round-off-limited (~1e-8) in the oracle itself; gate those bands separately, measure
      per-step increments against the path scale, and gate iteration counts. Quirks: "TANGENT
      parity at an exact Lode corner". Parity is a MODERATE-step statement (<= ~3e-3 strain; at 1e-2 both codes
      decide on round-off), the iteration-count gate applies to the fixed paths only, and near-coalescent
      eigenvalues are a third tangent band. Quirks: "Kernel-vs-O2 PARITY is a MODERATE-STEP statement",
      "THIRD round-off band".
- [ ] A kernel tangent taken w.r.t. the independent TENSOR shear component needs its shear COLUMNS halved by an
      engineering-shear shell (unlike `LadrunoJ2Kernel`). Quirks: "C is d(sigma_tensor)/d(eps_tensor)".
- [ ] Pin a number only from a step a RESIDUAL test converged, and re-measure pins after merging
      `ladruno`. A determinism gate needs no convergence: use `FixedNumIter`. Quirks: "where ONE
      tangent stopped".
- [ ] **[lint]** A `zone_a` test that branches on the platform declares `# ci-coverage:` (L8): PR CI
      is Ubuntu, so a win32-only leg never runs there. Gate only the MKL-specific leg. Quirks:
      "A win32-only `zone_a` test is NEVER run by PR CI".
- [ ] An inner material wrapped by `LogStrain` should provide its own trial elastic strain (mixin `LadrunoElasticStrainProvider`, NDMaterial FIRST base, engineering Voigt) unless it is linear-elastic: otherwise the wrapper's `inv(D0):tau` recovery is wrong. Quirks: "The `LadrunoElasticStrainProvider` mixin".

Found a new trap? Add it to `LEDGER_quirks.md`, then add one line here pointing to it. If the
trap has a greppable pattern, add a rule to `ci/check_quirk_patterns.py` instead.
