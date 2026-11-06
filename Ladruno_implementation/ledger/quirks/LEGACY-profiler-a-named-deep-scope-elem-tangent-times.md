---
wp: LEGACY
title: "Profiler: a named deep scope (elem.tangent) times the WHOLE loop including the addA scatter, while its elem_by_type rows time ONLY the kernel call — reading th…"
legacy_seq: 216
---
### Profiler: a named deep scope (`elem.tangent`) times the WHOLE loop **including the `addA` scatter**, while its `elem_by_type` rows time ONLY the kernel call — reading the scope as "element cost" overstates the threadable fraction
- **Bites:** you read `formTangent` (or its `elem.tangent` child) as the element fraction, apply a ">40% element work ⇒ thread it" gate, and over-predict the win. `OPS_PROFILE_SCOPE_DEEP_NAMED(_ops_elemTan, "elem.tangent")` (`IncrementalIntegrator.cpp:109`) wraps the **entire** element loop — `getTangent` *and* `theSOE->addA` — whereas the per-classTag `elem_by_type` bucket, filled by `OPS_PROFILE_FE_ELEM_SCOPE` (`:117`), covers only the `getTangent` call. Threading buys the second, not the first.
- **So:** `scatter = scope_wall − Σ(elem_by_type wall)` is the non-threadable remainder, and it is not small. Measured (ADR-75b L3-0, `lane3/RESULTS_l3a_update_scope.md`): the `addA` share of the tangent loop is **19.1%** on Lane B under UmfPack, **9.5%** under PARDISO — and on Lane A's cheap `forceBeamColumn` tangents the **scatter (44.1 ms) EXCEEDS the kernel (35.4 ms)**, i.e. that assembly loop is scatter-bound, not element-bound.
- **Also:** a loop scope appears at **several places in the tree** — `elem.update` shows up under `newStep`, under `solveCurrentStep/update`, and (Lane A) directly under `solveCurrentStep` from `DisplacementControl`. Summing only one site undercounts by ~2× on Lane A and ~2× on Lane D. Same trap ADR-40b hit with the hidden second `soe.factor`.
- **Workaround/status:** use `Ladruno_files/testbed/perf/lane3/parse_lane3.py`, which sums every site and reports kernel-vs-scatter per loop. *2026-07-25 (ADR-75b L3-0).*
