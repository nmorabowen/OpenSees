---
wp: WP-260
title: "TransformationConstraintHandler::enforceSPs() double-updates every constrained element outside elem.update"
date: 2026-11-26
pr: "#1158"
---
### `TransformationConstraintHandler::enforceSPs()` double-updates every constrained element outside `elem.update`

- **Bites:** Any model using `constraints Transformation` with fixed nodes or `equalDOF` ties — notably `DisplacementControl` (which calls `applyLoadDomain()` → `enforceSPs()` and then `updateDomain()` on every Newton/Krylov iteration) and explicit/implicit `newStep()` (which calls `applyLoadDomain()` before `updateDomain()`).
- **Why:** `TransformationConstraintHandler::enforceSPs()` (`SRC/analysis/handler/TransformationConstraintHandler.cpp:502-527`) loops over every `FE_Element` touching a constrained node and calls `theEle->updateElement()`. Immediately afterward, `updateDomain()` calls `Domain::update()`, updating **every** element in the domain again. Materials without a same-strain early return (such as `LadrunoConcrete3D::setTrialStrain`, and `LadrunoCPML::update` on jacket boundary elements) integrate the identical strain twice. In the profiler tree, the `enforceSPs()` call runs outside `elem.update`, so its time appears directly under `update` or `newStep` (e.g., ~29 % of D1b `LadrunoBrick` + `LadrunoConcrete3D` with a fixed base + `equalDOF` cap, and ~18–29 % of `newStep`/`update` on D3a/D3b CPML jackets).
- **Workaround/status:** Accounted in its own named column (`upd_constrained`) in `Ladruno_scripts/wp260_profiles/summarize.py`; proposed as a standalone optimization issue (not modified in `260-w0`). Learned 2026-11-26 (#1158).
