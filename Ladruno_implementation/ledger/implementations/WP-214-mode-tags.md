---
wp: WP-214
title: "modalProperties tags residual (R) / Ritz (Z) modes; tighter Penalty state-check margin in residualModes"
date: 2026-10-26
pr: "#TBD"
status: "draft"
banner: "modalProperties -print tags residual / Ritz modes (R / Z)"
---
| **WP-214 — mode-kind tags in `modalProperties`** ([[213_ladruno_modal_basis_completeness_adr]] D6) — when the WP-213 `LadrunoModeKinds` registry holds a valid record for the domain whose eigenvalues match the report's bit for bit, every mode row of `modalProperties -print` / `-file` (sections 2 and 6–10) ends in a trailing comment `  # R` (residual pseudo-mode) or `  # Z` (Ritz vector), eigen rows unmarked, and section 2 adds one `# NOTE:` line per kind (index ranges; Rayleigh periods are not natural; the `-period $T` Sa override when one is recorded). Rows stay numeric for a `#`-aware reader. `-return` gains a `modeKind` key (0/1/2) only then. No record (plain `eigen`, or a new `eigen` after `residualModes`): byte-identical to bac089465 (pinned baseline). Hook `LadrunoModeKinds::report(Domain*, eigenvalues)` on the Ladruno side; the vanilla file only appends its strings. Also: `residualModes` state check now refuses above max(1e-4, 10 f) and warns above max(1e-7, f) (was 100 f / 10 f): a 20 % stiffness loss on the tied-floor Penalty model is refused. No classTag. | Report decoration + guard | — | `SRC/analysis/analysis/{LadrunoRitzKernel.h,LadrunoRitzKernel.cpp,LadrunoResidualModes.cpp}`, `SRC/domain/domain/DomainModalProperties.cpp` (vanilla row), `tests/test_wp214_mode_tags.py`, `tests/data/wp214_modalprops_baseline_v3.txt` | draft | #TBD |
