---
wp: WP-140
title: "Quirk lint L7 sequence — unsequenced arithmetic on element accessors (WP-140)"
pr: "#880"
status: "draft"
section: "table"
legacy_seq: 17
---
| **Quirk lint L7 `sequence` — unsequenced arithmetic on element accessors (WP-140)** — `ci/check_quirk_patterns.py` flags a statement that combines TWO or more calls to the reference-returning element accessors (`getResistingForce*`, `getRayleighDampingForces`, `get*Force*`, `getTangentStiff`, `getInitialStiff`, `getMass`, `getDamp`) with `+ - *` — the shape of WP-124 C15 (vanilla `Element::getResponse` `inertialForce`, EXACTLY 0.0 on GCC). Scans every `SRC` file, vanilla included (like L5/L6); calls that are separate function arguments are not flagged; waiver `// ladruno-lint: sequence-ok <reason>`. **Proven on its incident:** flags `Element.cpp:513` of `a2004e0a7^` and is silent on `a2004e0a7`. **Current tree:** one hit, vanilla `IGAKLShell.cpp` (Vector + Matrix·Vector; audited harmless) — waived by owner decision. 14 self-tests (incident in a vanilla file, cross-type, pointer receiver in a non-element file, 8 must-pass shapes, waiver + short reason + stale); breaking the arithmetic test fails 4 of them, dropping the argument-comma guard fails 1. Cross-type is kept on purpose: LadrunoBrick's `getMass()` writes `resid` as a side effect. Guide: `ladruno-new-element` [lint] item. | lint rule | — | `ci/check_quirk_patterns.py`, `ci/test_check_quirk_patterns.py`, `SRC/element/IGA/IGAKLShell.cpp` (comment), `ladruno-new-element` guide, LEDGER_quirks, LEDGER_vanilla_files | draft | [#880](https://github.com/nmorabowen/OpenSees/pull/880) |
