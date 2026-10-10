---
wp: LEGACY
title: "Cyclic softening RC shell wall (ASDShellQ4 + LadrunoRCConcrete) walls on Newton convergence past first crack"
legacy_seq: 80
---
### Cyclic softening RC shell wall (ASDShellQ4 + LadrunoRCConcrete) walls on Newton convergence past first crack
- **Bites:** building a meshed reinforced RC shear wall — `ASDShellQ4` grid on a `LayeredShell`
  (`LadrunoRCConcrete` concrete layers + `PlateRebar(Steel02)` web steel) — and pushing it cyclically
  under `DisplacementControl` to demonstrate panel-scale pinching. It assembles fine and runs the first
  small-drift cycle (real `V`–`δ` hysteresis with dissipation), then `analyze` returns `-3`
  (`NormDispIncr` stalls with a large residual `deltaR`) at the next, larger amplitude — even with
  `-implex`, `KrylovNewton`, and `NewtonLineSearch`/`ModifiedNewton`/`Broyden` fallbacks.
- **Why:** the softening plastic-damage + crack-localization makes the global tangent indefinite/ill-
  conditioned at the load-redistribution events (crack formation, interlock cap engaging across a row of
  elements); a load-/displacement-controlled Newton has no way around the limit/snap-back points. IMPL-EX
  helps the MATERIAL tangent stay SPD-secant but does not fix the STRUCTURAL snap-through. This is the
  textbook reason squat-wall validation is hard, not a bug in the material (the material-point + single-
  shell-element cyclic gates all pass).
- **Fix (the deferred path, not yet done):** drive with an arc-length / indirect-displacement control
  (`LadrunoIndirectControl` / `LadrunoArcLength` are built for exactly this snap-back), or dynamic
  relaxation / quasi-static transient with mass; finer substeps; possibly an IMPL-EX-error step-cut.
  Harness scaffold: `tests/_testbed/rc_wall_harness.py`. Learned 2026-06-17 building the Phase-2b.2c.4
  Tran–Wallace pinching validation for [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]].
