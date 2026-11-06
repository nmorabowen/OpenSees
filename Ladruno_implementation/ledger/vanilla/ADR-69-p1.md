---
wp: ADR-69
title: "ADR-69 P1 -- 3 vanilla row(s)"
files: ["`SRC/element/absorbentBoundaries/LysmerTriangle.{h,cpp}`", "`SRC/analysis/integrator/ExplicitBathe.{h,cpp}`", "`SRC/modelbuilder/tcl/TclModelBuilder.cpp`"]
table: "main"
legacy_seq: [234, 235, 238]
---
| `SRC/element/absorbentBoundaries/LysmerTriangle.{h,cpp}` | `// Ladruno` ADR-69: incident-injection leak publisher. `commitState()` integrates the leak rate `R_injᵀv` (`R_inj = getDamp()·v_gnd` RECOMPUTED, not read from the stage-3-mutating `internalForces` member) trapezoidally per commit and publishes to `Ladruno::EnergyChannelRegistry` (ABSORB_LEAK); `setDomain` declares the channel. 3 diagnostic members, deliberately not serialized. Strictly additive — force path untouched. | ADR-69 P1 |
| `SRC/analysis/integrator/ExplicitBathe.{h,cpp}` | `// Ladruno` ADR-69: LNVD dissipation publisher — `addLocalDamping()` accumulates `Σ α·abs(r_i)·abs(v_i)·dt_substep` into a per-step scratch (sub-step 1 = p·dt, sub-step 2 = (1−p)·dt via the `lnvdInSubStep2` flag set around update()'s internal formUnbalance), `commit()` publishes to LNVD_WORK before `commitDomain()` (recorders read the registry inside), `newStep()` zeroes the scratch so failed/reverted steps never publish; ctor declares the channel when `-lnvd`. | ADR-69 P1 |
| `SRC/modelbuilder/tcl/TclModelBuilder.cpp` | `// Ladruno` ADR-69: `eleLoad -type -lysmerVelocityLoader <dir>` hook (+include). The `LysmerVelocityLoader` class existed since 2011 but NO interpreter path created it — the compliant-base incident-velocity input for `LysmerTriangle` was unreachable, and the `eleLoad` tail *silently returns TCL_OK* for unknown `-type` flags so decks "using" it were no-ops (ADR-69 P0.5 finding). Additive branch modeled on `-BrickW`; dir 1|2|3 validated. | ADR-69 P1 |
