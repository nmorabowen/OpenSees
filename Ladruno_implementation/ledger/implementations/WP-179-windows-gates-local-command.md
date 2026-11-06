---
wp: WP-179
title: "Windows-only gates as one local command: an engine-pinning pytest launcher plus a self-maintaining platform test list (no runner, by owner choice)"
pr: "#965"
date: 2026-10-08
status: "draft"
---
| **WP-179 — the Windows-only gates as one local command** ([[BUILD_GOTCHAS]] §4c; issue #934). There is no Windows CI runner, by the owner's choice (a self-hosted runner on a public repo, used on demand, added little over a local command). `Ladruno_scripts/ci_run_pytest.py` refuses to run without `python -S` (exit 91), pins `dist\bin` for itself and its children, asserts the engine file, and optionally requires `ladrunoBuild() == --expect-sha` (exit 90; both codes sit outside pytest's 0–5). `python ci/check_quirk_patterns.py --list-platform-tests` lists every collectable `test_*.py` of ANY tier with a non-portable platform branch, from lint L8's own AST walk, so no list is maintained. That's 12 files today, including 3 non-`zone_a` files that no CI ran. Together: `py -3.12 -S Ladruno_scripts/ci_run_pytest.py -- $(py -3.12 ci/check_quirk_patterns.py --list-platform-tests)`, measured at 175 passed + 2 slow skipped in ~23 min. Also fixes `test_wire_venv_pth_override` (its probe child now starts with `-S`). | test tooling + docs | — (no class tag) | `Ladruno_scripts/ci_run_pytest.py`, `ci/check_quirk_patterns.py`, `ci/test_quirk_ci_coverage.py`, `tests/test_wire_venv_pth_override.py`, `Ladruno_internal/BUILD_GOTCHAS.md` | draft | [#965](https://github.com/nmorabowen/OpenSees/pull/965) |
