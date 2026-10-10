---
wp: WP-162
title: "Two PRs taking the next free lint number merge their self-tests silently; F811 misses it"
pr: "#916"
date: 2026-10-03
---
### Two PRs that take "the next free" lint number merge their self-tests SILENTLY into one file -- and ruff F811 does not see it (WP-162, 2026-10-03)
- **Bites:** WP-153 and WP-158 both took quirk rule L9 and both added `def _l9(...)` to `ci/test_check_quirk_patterns.py`. The #899 merge-up of #901 merged the file with no conflict and two `def _l9`; Python keeps the second, so every WP-153 test silently ran WP-158's helper. Caught by eye, renumbered L10 by hand.
- **Why F811 misses it:** pyflakes F811 is "redefinition of UNUSED name", and the first `_l9` is referenced by the tests between the two definitions. Reproduced on the real three-way merge (e2560dd2b / 5300da720): ruff 0.16.10 `--select F811` passes. Ruff has no function-redefined rule (pylint E0102).
- **Status:** rules are named by slug in one `RULES` registry with import-time uniqueness asserts, one `ci/test_quirk_<slug>.py` per rule, and `ci/check_duplicate_defs.py` (D1) fails any `def`/`class` defined twice in one scope; ruff F811 runs beside it. [[162_lint_rule_slugs]]
