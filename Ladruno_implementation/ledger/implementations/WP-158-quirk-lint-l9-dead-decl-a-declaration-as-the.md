---
wp: WP-158
title: "Quirk lint L9 dead-decl — a declaration as the whole unbraced body that shadows (WP-158)"
pr: "#901"
status: "draft"
section: "table"
legacy_seq: 16
---
| **Quirk lint L9 `dead-decl` — a declaration as the whole unbraced body that shadows (WP-158)** — `ci/check_quirk_patterns.py` flags a declaration that is the ENTIRE unbraced body of an `if` / `else` / `for` / `while` and re-declares a parameter or local declared earlier in the same function — the shape of vanilla `ManzariDafalias::ForwardEuler`'s `Vector r(6); if (p > small) Vector r = ...;` (the outer `r` stayed zero). Scans every `SRC` file, vanilla included (like L5–L7); waiver `// ladruno-lint: decl-ok <reason>`. An UNshadowed dead declaration is deliberately not flagged: vanilla `Domain::initialize` declares `Matrix initM(ele->getInitialStiff())` as a loop body on purpose. **Proven on its incident:** flags `ManzariDafalias.cpp:1581` of d63f49750 and is silent after the fix; 0 findings on the rest of `SRC`. 18 self-tests (incident in a vanilla file, the fix, 4 shadowing shapes, 10 must-pass shapes incl. `if constexpr` and the Domain idiom, waiver + short reason + stale); mutants: dropping the shadow test fails 1, the keyword filter 2, the control-prefix test 2; the `{` guard is defensive (no self-test isolates it). Adds ~15 s to a full lint run. Guide: `ladruno-new-material` [lint] item. | lint rule | — | `ci/check_quirk_patterns.py`, `ci/test_check_quirk_patterns.py`, `ladruno-new-material` guide, `ci/README.md`, LEDGER_quirks | draft | [#901](https://github.com/nmorabowen/OpenSees/pull/901) |
