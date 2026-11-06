---
wp: LEGACY
title: "Mortar friction gT0/engaged are captured in getResidual and NOT reverted — latent until the C3.2 implicit tangent"
legacy_seq: 126
---
### Mortar friction gT0/engaged are captured in getResidual and NOT reverted — latent until the C3.2 implicit tangent
- **Bites:** ADR-41 C3.2 (NOT C3.1). `revertToLastCommit` drops only `gpTtrial=gpT` for mortar slots
  (`LadrunoContactDomain.cpp`); the engagement origin `gT0`/`engaged` (set once in `addMortarFriction`) are
  never reverted. A rejected Newton step that FIRST-engages a node latches `gT0` from the rejected trial
  config; the retry keeps that stale origin (`engaged` stays true) ⇒ a spurious stick offset. Identical to
  the shipped NTS SEGMENT behavior (which also doesn't revert `engaged`), so NOT a C3.1 regression, and
  **unreachable under C3.1's explicit-only path** (CDL never reverts mid-step). It goes live when the C3.2
  friction tangent lands and an implicit Newton step is rejected. **RESOLVED in C3.2 (#378):** `gT0`/`engaged`
  are double-buffered (`gT0committed`/`engagedCommitted`), promoted in `commit()` and restored in
  `revertToLastCommit()`. Found by the C3.1 gate (MAJOR-2, #377), fixed in C3.2.
- **2026-07-02 (contact-review P2): the SAME fix was finally BACK-PORTED to the NTS lane it was copied
  from.** The paragraph above ("identical to the shipped NTS SEGMENT behavior") documented the shared flaw
  but only excused it as not-a-C3.1-regression — and NTS friction is a first-class IMPLICIT path since
  P3.5 (#361), so the "unreachable under explicit" shield never applied there. NTS `FrictionState` now
  carries the same `gT0committed`/`engagedCommitted` double-buffer (commit promotes / revert restores).
  LESSON: when a gate fixes a defect on a COPIED lane, grep for the source lane the copy came from — a
  ledger sentence acknowledging "the sibling has it too" is a fix obligation, not an absolution. Same PR:
  `LadrunoContactDomain::revertToStart()` (hooked from `Domain::revertToStart` — `ops.reset()`) drops ALL
  contact path state (friction slip/origins, ALM λ_N/λ_T/λ_tie, edge signs, re-emit anchors/fp/trigger,
  NTS-force + nodal-mass caches; NormalField σ+frozen sign KEPT — reference-geometry datum, re-derivable
  identical). Pre-fix, a re-run after `ops.reset()` started from the previous run's committed slip gpT ⇒ a
  large spurious backward stick force at first contact — silently different from a fresh model. Gates:
  `tests/test_contact_review_p2_lifecycle.py` (failed-step retry ≡ never-failed reference bit-tight;
  reset re-run ≡ first run bit-tight — both FAIL pre-fix).
