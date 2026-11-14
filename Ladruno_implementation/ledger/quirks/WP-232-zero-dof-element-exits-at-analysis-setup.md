---
wp: WP-232
title: "An element with no DOFs in the Domain kills the process at the first analysis set-up: FE_Element exits on zero DOFs, and the handler wraps every element"
date: 2026-11-13
---
### An element with no DOFs in the Domain kills the process at the first analysis set-up — `FE_Element` calls `exit(-1)`, and the handler wraps every element

- **Bites:** any design that parks a zero-node `Element` in the Domain to borrow the Domain's lifecycle — to receive `commitState`, `revertToLastCommit`, `revertToStart`, `sendSelf` or `setParameter` for an engine that is not a mesh (proposed for WP-232's explicit MPM engine in review, and rejected). It works until something builds an analysis in that Domain; then the process exits, with no loud refusal the user can act on.
- **Why:** the `FE_Element` constructor calls `exit(-1)` when the element has no DOFs (`FE_Element.cpp:126-129`, "element must have 1 dof"), and `PlainHandler::handle()` creates an `FE_Element` for every element the Domain holds, without checking `isActive()` (`PlainHandler.cpp:291-327`). The other handlers were not checked.
- **Workaround / status (2026-11-13, WP-232 — disclosed, not changed):** host such an engine as a `LoadPattern` subclass instead, after ADR-73's `LadrunoPorousOverlay`: patterns get no `FE_Element`, `wipe` deletes them (`Domain::clearAll`), and the `parameter … pattern <tag>` target reaches their `setParameter`. A pattern does not receive `revertToStart` (only a re-`applyLoad`), so the engine needs its own reset command.
