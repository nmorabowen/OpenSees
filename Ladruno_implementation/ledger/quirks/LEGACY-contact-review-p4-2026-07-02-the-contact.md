---
wp: LEGACY
title: "Contact-review P4 (2026-07-02): the contact handler silently DROPPED every equationConstraint — LadrunoTie + contact = dead ties"
date: 2026-07-02
legacy_seq: 129
---
### Contact-review P4 (2026-07-02): the contact handler silently DROPPED every equationConstraint — LadrunoTie + contact = dead ties
- **Bites:** any model combining `constraints LadrunoContact` with `equationConstraint` — notably the
  fork's own **LadrunoTie (ADR-62)**, which emits ordinary EQ rows. The handler REPLICATES PlainHandler's
  DOF/FE loop (so the contact-FE start tag is knowable), but upstream PlainHandler gained an
  EQ_Constraint block AFTER the replica was written (enforce trivial-identity ⇒ DOF −4; loudly
  warn-and-ignore non-trivial); the replica had NEITHER ⇒ ties ran completely dead with NO diagnostic —
  the exact silent-zero class the handler's non-homogeneous-SP warning guards (review HIGH-3).
- **FIX (upstream parity, no more):** the EQ block is ported — trivial-identity EQs are ENFORCED
  (DOF mark −4, the numberer's EQ code) and non-trivial ones are warned NOT-ENFORCED, pointing at
  `constraints LadrunoProjection`. **RULE: contact + non-trivial ties in ONE analysis remains
  unsupported** (the handlers are mutually exclusive; actual EQ enforcement is the projection
  handler's job) — it now says so instead of silently producing a wrong answer. Gate:
  `tests/test_contact_review_p4_eq_parity.py` (capfd warning assert — FAILS pre-fix on silence).
