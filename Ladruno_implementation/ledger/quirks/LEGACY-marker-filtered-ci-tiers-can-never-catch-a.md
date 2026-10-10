---
wp: LEGACY
title: "Marker-filtered CI tiers can never catch a cross-test state leak — the tiers must share a process"
legacy_seq: 295
---
### Marker-filtered CI tiers can never catch a cross-test state leak — the tiers must share a process
- **Bites:** `.github/workflows/ladruno.yml` ran `pytest -m "zone_a"` on PRs (ubuntu) and `pytest -m "zone_b"` nightly (self-hosted). The two tiers NEVER shared an interpreter, so a test in one tier that poisons a process-global was invisible to both jobs by construction. Measured on one build: `zone_b` alone is 144 passed; `zone_a`-then-`zone_b` is 1 failed. That is exactly how the ASDConcrete3D tag-cache defect survived — and it would have hidden any other global-state defect just as well.
- **Why it evades the usual guards:** every job is green, every tier is "covered", and the coverage gap is in the SEAM between jobs rather than in any one of them. Running the tiers separately is also the natural thing to do (different runners, different dependencies, different durations), so nothing looks wrong.
- **What catches it:** one job that runs the WHOLE suite unfiltered in a single process. Added as `cross-tier-nightly` (~50 min, self-hosted, nightly + workflow_dispatch).
- **If you must exempt a known failure, pair the exemption with a SENTINEL.** `cross-tier-nightly` deselects the one known-unfixed ASDConcrete3D failure — otherwise the job is red every night and gets ignored, which catches nothing — and then runs the minimal poisoning PAIR and requires it to FAIL. If the defect is ever fixed, the sentinel goes green, the job goes red, and whoever fixed it is told to delete the exemption. A bare exemption rots silently; an exemption with a sentinel cannot.
- **`--deselect` with an unmatched path is SILENTLY IGNORED** — no error, no warning, exit 0. A typo'd deselect leaves the job red forever with nothing explaining why. Verify by collection count, with a deliberately bogus path as the control: real deselect gave `1929/1930 (1 deselected)`, bogus gave `1930 collected` and no complaint.
