---
wp: LEGACY
title: "-noMassCache is per-instance and NOT serialized — the A/B escape hatch is effective on the building rank only"
legacy_seq: 241
---
### `-noMassCache` is per-instance and NOT serialized — the A/B escape hatch is effective on the building rank only
- **Bites:** under SP/MP/DDM, broker-built remote copies of a cache-bearing element default-construct the cache **enabled** (the flag rides no sendSelf vector — deliberate T7/brick policy). A user running `-noMassCache` as a workaround or as the off-arm of a parallel A/B gets the cache silently re-enabled on every non-building rank, so a "cache-off" parallel comparison is not what it claims to be. Results remain bit-identical by the guard construction (G-BYTE), so this is a diagnostics honesty gap, not a correctness one — but it is precisely the parallel context where an escape hatch matters most.
- **Workaround/status:** documented in the `LadrunoMassCache.h` contract (review wave); serializing the flag would cost a stream-format change on six elements for a diagnostic — declined. Treat serial runs as the authoritative A/B arena. *2026-07-27 (ADR-77 review wave).*
