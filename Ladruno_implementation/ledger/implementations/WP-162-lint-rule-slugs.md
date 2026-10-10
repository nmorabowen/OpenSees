---
wp: WP-162
title: "Quirk-lint rules by slug; per-rule tests; duplicate-definition gate"
pr: "#916"
date: 2026-10-03
status: "draft"
---
| **WP-162 — quirk-lint rules by slug, one self-test file per rule, a duplicate-definition gate** ([[162_lint_rule_slugs]]). `ci/check_quirk_patterns.py` registers its rules once in `RULES` (`rayleigh`, `wipe`, `pointers`, `commit`, `double-load`, `ground-sign`, `sequence`, `ci-coverage`, `dead-decl`, `revert`; import-time uniqueness asserts); findings print the slug; `--only` takes slugs, the L-numbers are deprecated aliases for one release; `--rules-table` generates the `ci/README.md` table. `ci/test_check_quirk_patterns.py` split into `ci/test_quirk_<slug>.py` + `ci/_quirk_testkit.py`. New static-gates step: `ci/check_duplicate_defs.py` (D1, an `ast` function/class-redefined check) + `ruff check --select F811` (`ruff==0.16.10`); F811 alone misses the two-`def _l9` merge of #899/#901. | CI tooling | — (no class tag) | `ci/check_quirk_patterns.py`, `ci/test_quirk_*.py`, `ci/_quirk_testkit.py`, `ci/check_duplicate_defs.py`, `ci/test_duplicate_defs.py`, `ci/README.md`, `.github/workflows/ladruno.yml` | draft | [#916](https://github.com/nmorabowen/OpenSees/pull/916) |
