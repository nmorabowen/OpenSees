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
      you deliberately make it process-wide. A latch set in place must be copied by `getCopy`,
      and a loud failure must not latch. Quirks: "getCopy must PROPAGATE the latch".

## Return map, tangent, substepping

- [ ] A non-converged return map must FAIL (return < 0), never commit `f > 0` as success, and
      the refusal must reach `analyze`: several vanilla elements swallow it. Quirks:
      "`Backward_Euler` ACCEPTS a non-converged return map", "swallow the material's refusal".
- [ ] Implement `getInitialTangent()` honestly: the base default returns `getTangent()`, so
      `-initial` silently becomes full Newton. Quirks: "`NDMaterial::getInitialTangent()` DEFAULTS".
- [ ] Substep schemes need error control and yield-drift correction, and must honour the
      tolerance passed in. Quirks: "`IntScheme` 3 (RungeKutta4) and 5 (ForwardEuler) have no
      error control", "IntScheme 1 (ModifiedEuler) IGNORES the `TolR`".
- [ ] IMPL-EX in a static analysis: `ops_Dt` is pseudo-time and erratic; guard the
      extrapolation factor. Quirks: "IMPL-EX in a STATIC analysis".
- [ ] `revertToStart()` must not reset calibrated constants mid-analysis. Quirks:
      "`ManzariDafalias::revertToStart()` silently restores".
- [ ] Constants that multiply a stress are dimensional: make them unit-consistent or document
      the units. Quirks: "`D_factor` dilatancy sigmoid is DIMENSIONAL".

## NaN and silent success

- [ ] Never check divergence with `pNorm(0)` (NaN-blind); use `std::isfinite`. Quirks:
      "`Vector::pNorm(0)` is NaN-BLIND". Eigen `*= 0` keeps NaN garbage:
      "with `*= 0` keep NaN heap garbage".
- [ ] `analyze()` returning 0 does not mean the numbers are finite. Quirks: "`analyze()` returns
      rc=0 on a NaN-poisoned system".

## Finite strain and regularization

- [ ] Don't lift a damage/softening material with the generic `LogStrainNDMaterial` wrapper.
      Quirks: "The generic LogStrainNDMaterial wrapper is UNSOUND".
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

Found a new trap? Add it to `LEDGER_quirks.md`, then add one line here pointing to it. If the
trap has a greppable pattern, add a rule to `ci/check_quirk_patterns.py` instead.
