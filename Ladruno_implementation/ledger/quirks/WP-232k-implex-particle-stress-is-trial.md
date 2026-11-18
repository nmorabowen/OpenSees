---
wp: WP-232k
title: "Under -implex an MPM particle's stress output is the trial (extrapolated) stress, not the committed one"
date: 2026-11-18
---
### `mpm particleData sig` / `sigHat` under SANISAND `-implex` show the extrapolated stress

- **Bites:** anyone reading particle stress from an MPM run whose SANISAND particles use
  `-implex`, in either mode, and comparing it with an FE record or with the material's own
  state. The engine keeps the stress the material returned at the trial call
  (`setTrialStrain` → `getStress`). Under `-implex` that is the extrapolated (companion) stress
  the step was equilibrated with, the stress an FE element's internal force uses. `commitState`
  then corrects the material, but the particle's `sigHat` (and `sig`, its push-forward) keep the
  trial value. An FE recorder reads the material after commit, so it shows the corrected one.
- **Measured** (build `f235a001f`, S1 one-point option legs, all 12 legs and both modes): the
  material's committed stress and its 26-entry state equal the scheme-1 run bit for bit. The
  particle's `sigHat` departs from step 11, the first plastic step, by 0.035 (ps oedo) to 0.233
  (3D shear) relative to scheme 1. Numbers: `tests/data/mpm/soil/s1_options_diagnosis.json`.
- **Read the committed stress with** `mpm particleResponse stress` (the material's own response,
  material frame, per particle). `mpm particleResponse state` is the committed state.
- **Not changed** (ruling d683e0db8): keeping the trial stress is consistent with IMPL-EX. The
  `-implex` boundary-value leg (K0p column, mode A gimp) passes on particle stress at 0.01862
  ρgl0 against 0.02, because after a hold the implex error is small. Test:
  `tests/test_ladruno_mpm_soil_options.py::test_option_implex_material_equals_scheme1`.
