---
wp: LEGACY
title: "-implex looks INERT under ASDShellQ4 in fully-prescribed rigs — it is not: ASDShellQ4 reports the POST-COMMIT state, and IMPL-EX re-integrates implicitly at co…"
legacy_seq: 155
---
### `-implex` looks INERT under ASDShellQ4 in fully-prescribed rigs — it is not: ASDShellQ4 reports the POST-COMMIT state, and IMPL-EX re-integrates implicitly at commit

- **Bites:** A/B-ing `-implex` on/off (LadrunoRCConcrete) with a fully prescribed
  (every-DOF `sp`) single-element bending rig under `ASDShellQ4` + `LayeredShell`:
  the recorded responses are BIT-IDENTICAL (rel 4.7e-16 over 60 softening steps) —
  it looks like the flag is dropped somewhere in the section copy chain. It isn't:
  the same rig under `ShellMITC4` (same section) or `LadrunoSolidShell` (3D view)
  shows the expected ~3.5% IMPL-EX extrapolation difference.
- **Why:** LadrunoRCConcrete's IMPL-EX `commitState` re-integrates implicitly to
  advance the TRUE thresholds (the ASDConcrete3D pattern), so the post-commit
  material state is the implicit one. ASDShellQ4's reported element forces reflect
  the post-commit section state; in a rig with NO free DOFs the extrapolated trial
  stresses are never consulted by any equilibrium iteration, so the recorded curve
  collapses exactly onto the implicit run. ShellMITC4 / LadrunoSolidShell report
  the converged TRIAL state, where the extrapolation lives.
- **Workaround/status (measured 2026-07-07, ADR-66 G7):** not a defect on either
  side, but three consequences: (1) never "verify implex engaged" with a prescribed
  probe under ASDShellQ4 — it is structurally invisible there; (2) cross-element
  parity benchmarks must run implicit-vs-implicit (as G7 does) or accept a
  reporting-path asymmetry masquerading as element deviation (~2% here); (3) on
  free-DOF problems implex under ASDShellQ4 IS active (0.7% path shift measured on
  the G7b rig) — the wall-harness usage is fine. Pinned discriminatingly by
  `tests/test_ladrunoSolidShell_flexure.py::test_implex_reporting_paths`.
