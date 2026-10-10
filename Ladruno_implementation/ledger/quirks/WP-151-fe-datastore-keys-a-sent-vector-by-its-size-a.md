---
wp: WP-151
title: "FE_Datastore keys a sent Vector by its SIZE: a fork sendSelf block of the SAME length as the base's vector, under the same dbTag and commitTag, OVERWRITES the…"
legacy_seq: 531
---
### FE_Datastore keys a sent Vector by its SIZE: a fork `sendSelf` block of the SAME length as the base's vector, under the same dbTag and commitTag, OVERWRITES the base state (WP-151)
- **Bites:** `LadrunoSANISAND::sendSelf` sends two vectors, both with `this->getDbTag()` and `commitTag`: the base `ManzariDafalias` state as a `Vector(97)`, then its own Ladruno block. FileDatastore files vectors per `<size>.<commitTag>`, then by dbTag.
  - WP-151 added six entries, and the Ladruno block became exactly 97 long (35 + 17 + 6 + 36 + 3). It landed in the base's slot.
  - Every database round trip then restored the material on another state, silently.
  - Four existing tests caught it: `test_db_roundtrip_carries_presidual`, the two `-implex` round trips, and `test_pre_floor_crosses_the_datastore_wire`.
- **Rule:** A subclass that appends its own send block to a base `sendSelf` under the same dbTag and commitTag must give that Vector a length different from every length the base sends. Sizes are arithmetic, so no grep finds this: `static_assert` it where the size is defined.
- **Workaround/status:** ✅ Fixed (WP-151) with three pieces:
  - a trailing layout tag, `LWIRE_TAG`;
  - ONE `static_assert(LWIRE_SIZE != 97)`, beside the `LWIRE_*` enum in `LadrunoSANISAND.h`;
  - `recvSelf` REFUSES (returns −1) on a tag mismatch, before assigning anything (review #893 M3).

  After the #868 merge the block is 117 long. `ladruno` before WP-151 had 110, and WP-151's first layout had 97.
