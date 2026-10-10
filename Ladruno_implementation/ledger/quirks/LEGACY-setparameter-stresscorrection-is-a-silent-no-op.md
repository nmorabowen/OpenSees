---
wp: LEGACY
title: "setParameter ... stressCorrection is a silent no-op from every interpreter (vanilla ManzariDafalias, found 2026-09-07 by TIMs)"
date: 2026-09-07
legacy_seq: 387
---
### `setParameter ... stressCorrection` is a silent no-op from every interpreter (vanilla ManzariDafalias, found 2026-09-07 by TIMs)
- **Symptom:** `ops.setParameter('-val', 0, '-ele', ..., 'stressCorrection')` binds (a bad tag moves nothing, `refShearModulus` by the same route moves the stress) but changes nothing: `updateParameter` id 9 reads `info.theInt` (`ManzariDafalias.cpp:897`, same pattern at `:864` for `mElastFlag`) while the interpreters' `setParameter` fill only `theDouble`, and `Information` leaves `theInt` at its default. So the drift correction cannot be switched off from a deck; only the in-model default (ON, four constructors) is ever exercised.
- **Fix (owed, ADR 92 P2 list):** override `updateParameter` in `LadrunoSANISAND` for the ids that read `theInt` and accept `theDouble` too (zero vanilla footprint); until then any "stress-correction-off" arm in a memo was never run.
