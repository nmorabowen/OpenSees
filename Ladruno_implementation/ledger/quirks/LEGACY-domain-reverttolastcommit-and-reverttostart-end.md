---
wp: LEGACY
title: "Domain::revertToLastCommit() (and revertToStart()) END with return this->update(); — so every element's update() fires ONCE MORE on the just-reverted state"
legacy_seq: 439
---
### `Domain::revertToLastCommit()` (and `revertToStart()`) END with `return this->update();` — so every element's `update()` fires ONCE MORE on the just-reverted state
- **Bites:** any element that MUTATES state inside `update()` — a path-dependent internal
  variable, an augmented-Lagrangian multiplier, a counter — advances it one extra time on the
  state it has just rolled back to. Symptom is tiny and therefore lethal: an exactness test
  that should compare equal is off by one increment evaluated at the *converged* gap. Real
  instance (WP-101, `LadrunoKinematicCoupling`): the AL Uzawa update `λ += D g` was moved into
  `update()` so the constraint converges within a step; a deliberately failed step then left λ
  off its committed value by `1.5e-10` on one row of 27 — `D·g` at the committed state, where
  `g ~ 1e-14` is not exactly zero. On a run with many failed steps that ratchets.
- **Why:** `Domain.cpp` — `revertToLastCommit()` reverts every node and element, restores
  `currentTime`/`dT`, re-applies the load at the committed time, and then does
  `return this->update();` (the trial state has moved, so the elements must be re-informed).
  `revertToStart()` ends the same way. The revert therefore calls `Element::revertToLastCommit()`
  on every element FIRST and `Element::update()` on every element SECOND — a pairing nothing in
  the `Element` interface documents.
- **Also `Domain::recv()`** (the FE_Datastore / parallel restore path) calls `theEle->update()`
  immediately after `recvSelf`, with exactly the same consequence: a database save/restore moved the
  same element's λ by `6.6e-9`. Arm the latch in `recvSelf` too — clearing it there (the obvious
  "reset my transient flags" reflex) is precisely wrong.
- **Workaround/status (2026-09-14):** a **one-shot latch**: the element's
  `revertToLastCommit()`/`revertToStart()`/`recvSelf()` arms a transient `bool`, and the first `update()`
  after it consumes the flag and returns without advancing. Order is safe because the domain
  finishes ALL the reverts before it calls update at all. Do NOT serialize the latch — it only
  ever lives between those two calls. Same shape as the ADR-39 contact-pair revert note above:
  if your element has state in `update()`, it has this problem. See
  [[LEDGER_implementations]] WP-101 row / [PR #839](https://github.com/nmorabowen/OpenSees/pull/839).
