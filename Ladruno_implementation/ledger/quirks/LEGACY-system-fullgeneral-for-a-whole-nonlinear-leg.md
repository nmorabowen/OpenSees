---
wp: LEGACY
title: "system FullGeneral for a whole nonlinear leg costs ~40x — assemble densely only at the sampling points"
legacy_seq: 284
---
### `system FullGeneral` for a whole nonlinear leg costs ~40x — assemble densely only at the sampling points
- **Bites:** wanting the assembled tangent's spectrum along a collapse path, you switch the analysis to `FullGeneral` and let the leg run. Measured on an H20 leg at 2800 free DOF: **~20 s/step** under `FullGeneral` against a fraction of a second under `system Pardiso`, because every Newton *iteration* pays a dense factorisation. A 250-step leg becomes 80+ minutes, and near a wall — where every ladder rung runs its full iteration count before failing — it effectively stops advancing. A first attempt at this diagnostic burned ~25 min and produced **one** usable sample.
- **Why:** `FullGeneral` is an unblocked dense LU; the cost is per solve, and a Newton leg does 3–125 solves per step.
- **Rule:** keep the sparse solver for the leg and take the dense matrix only when you want it: `wipeAnalysis` → `system FullGeneral` → `integrator LoadControl 0.0` + `algorithm Linear` → `analyze(1)` → `printA -ret` → switch back. The **zero** increment leaves displacements, stresses and pseudo-time untouched (this is `h20_prandtl.py::leg_modes`' idiom applied mid-run), so it costs **one** dense factorisation per sample instead of one per iteration. Two companion traps: a sampler trigger keyed to `DS_BASE` fires only a few halvings from the floor and can yield **zero** converged samples (arm it when the controller first backs off its *maximum* step), and a leg can pass from trigger to termination without converging another step — so take one unconditional sample **after** the loop, at the last converged state.
- **Workaround/status:** implemented as `quad_path_diag.py::sample_tangent` (+ `--cond-at` / `--cond-every` and the terminal sample). *Learned 2026-08-11 (note 82).*
