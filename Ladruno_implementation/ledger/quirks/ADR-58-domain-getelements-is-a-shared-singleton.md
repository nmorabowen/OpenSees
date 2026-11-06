---
wp: ADR-58
title: "Domain::getElements() is a SHARED singleton iterator — never iterate the Domain from inside an element callback (ADR-58 P2-S2)"
legacy_seq: 136
---
## `Domain::getElements()` is a SHARED singleton iterator — never iterate the Domain from inside an element callback (ADR-58 P2-S2)

`Domain::getElements()` returns `*theEleIter` after calling `theEleIter->reset()` — **one shared
`SingleDomEleIter`**, not a fresh object (same for `getNodes()` etc.). `Domain::commit()` and
`Domain::update()` walk elements through that single iterator (`while ((e = theEleIter()) != 0)
e->commitState()/update()`). So if an element's `commitState()`/`update()`/`getResistingForce()`
calls **any** Domain method that re-iterates elements — most notably `Domain::calculateNodalReactions()`
(it does `getElements()` + `addResistingForceToNodalReaction` on every element) — the nested
`reset()` **rewinds the iterator the outer loop is using**, so the outer loop terminates early and
**silently skips `commitState()` on every element after the caller**. No crash, no warning — just
some elements never commit (subtly wrong results). This is NOT theoretical: it bit the rigid-body
moment gather, which needed the toe-spring reaction inside `commitState`.
- **Fix pattern:** never iterate the Domain from an element callback. If you need another element's
  force, **cache its tag at `setDomain`** (called from `Domain::addElement`, OUTSIDE any iteration —
  safe) and re-resolve with `Domain::getElement(tag)` + read `getResistingForce()` directly. Cache
  TAGS not `Element*` so a removed element is skipped (`getElement` returns 0), not deref'd after free.
  (`Node`/`Element` do not store incident elements, so there is no per-node shortcut — the setDomain
  scan is the way.)
