---
wp: LEGACY
title: "LadrunoBrick -formulation eas is a misnomer (it's SSPbrick) → rename to ssp RECOMMENDED"
legacy_seq: 60
---
## LadrunoBrick `-formulation eas` is a misnomer (it's SSPbrick) → rename to `ssp` RECOMMENDED

The `-formulation eas` single-point element is **never** a Simo–Rifai enhanced-assumed-strain
element. It is a **verbatim port of `UWelements/SSPbrick::GetStab`** (McGann/Arduino) — B̄ +
a statically-condensed assumed-strain hourglass `Kstab` baked once at `setDomain` on the
*initial* tangent (no per-step α). Proven by source (`LadrunoBrick.cpp:1887`) and by
byte-identity to `SSPbrick` (6 figs elastic, 4 figs plastic, 0.2% with damage — see the
Lemaitre validation §4.6 + `tests/test_ladrunoBrick_element.py::test_eas_matches_sspbrick_distorted_hex`).

- **RECOMMENDED rename (NOT yet in source):** call it `-formulation ssp`; keep `eas` as a
  deprecated alias (warn → ssp) so the name `eas` is freed for a **true** Simo–Rifai EAS element.
  Suggested impl: enum `Formulation::SSP` with `EAS = SSP` back-compat alias; parser branch in
  `OPS_LadrunoBrick.cpp`; `formulationName→"ssp"`. The owner will land the C++ rename separately
  (docs-only PR keeps the current `eas` name). Tests can use a runtime fallback (`ssp` if built,
  else `eas`) to stay green across the rename.
- **Consequence (why it's not a bug):** being a single-point element, `eas`/`ssp` over-predicts load
  in a plastic/damage stress gradient (Jensen on the concave σ(ε); centroid under-samples) — but it
  shares this *exactly* with the validated `SSPbrick`, and it converges to `bbar` under refinement
  (gap 17.4%→3.7% to h=0.5). Use `bbar` for gradient/fracture fields, `eas`/`ssp` for smooth + cost.
