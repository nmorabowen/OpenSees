---
wp: LEGACY
title: "FileDatastore silently CLOBBERS same-type objects stored with the same (dbTag, commitTag, SIZE) — pack multi-part sendSelf payloads under distinct dbTags"
legacy_seq: 184
---
### FileDatastore silently CLOBBERS same-type objects stored with the same (dbTag, commitTag, SIZE) — pack multi-part sendSelf payloads under distinct dbTags
- **Bites:** a class whose `sendSelf` sends TWO Vectors (or two IDs) on its one dbTag works fine — until a model size where the second object's length coincides with the first's; then the later write overwrites the earlier entry and `recvSelf` restores garbage (measured: LadrunoPorousOverlay configs where `nPI + 6·nLay + 2·nRN == 38` corrupted the scalar block → huge bogus counts → abort/exit-127 on restore; other coincidences restored silently-wrong state). Non-deterministic-LOOKING because it is config-size-dependent.
- **Why:** `FileDatastore` files entries per (type, size) and keys by (dbTag, commitTag) within that file — same type + same size + same tags = same slot.
- **Workaround/status (2026-07-14, ADR-73 P1):** the upstream matDbTag idiom — grab a second `theChannel.getDbTag()` for the payload object and transmit it inside the first (fixed-size) block. Overlay fixed this way; DB round-trip battery-gated. Audit any future multi-send class for size coincidences.
