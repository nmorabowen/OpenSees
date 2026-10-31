---
wp: LEGACY
title: "ManzariDafalias is ELASTIC until updateMaterialStage ... -stage 1, and the flag is process-wide — a tangent probe that forgets the flip measures Ce and \"passes\""
legacy_seq: 467
---
### `ManzariDafalias` is ELASTIC until `updateMaterialStage ... -stage 1`, and the flag is process-wide — a tangent probe that forgets the flip measures Ce and "passes"
- **Bites:** a material-point probe or FD tangent check on `ManzariDafalias`/`LadrunoSANISAND` finds `dGamma = 0`, TanType 0/1/2 all identical and equal to Ce, and concludes the tangent is fine. It never left the elastic branch.
- **Why:** `mElastFlag` is a `static` class member (one per process, not per instance); the full constructors set it to 0 (ELASTIC, #714), so `integrate()` calls `elastic_integrator` unconditionally and `GetElastoPlasticTangent` is never reached until `updateMaterialStage -material <tag> -stage 1`. Because it is static, constructing ANY new `ManzariDafalias` resets EVERY live instance to elastic — a rebuild-per-probe harness must re-flip after every build.
- **Workaround/status:** by design (staged-gravity idiom). Always flip after construction and assert `dGamma > 0` (state slot 25) before trusting a plastic-state measurement — `tests/test_manzari_ep_tangent_gate.py` does both.
#### Follow-up (apeGmsh #2444, 2026-10-22)
- Constructing a new `ManzariDafalias` no longer resets the stage: only `wipe` does (the default is the file-scope `0` plus the `OPS_clearAllNDMaterial()` reset). A rebuild-per-probe harness that calls `wipe` between builds still starts elastic and must still flip; one that rebuilds WITHOUT `wipe` now keeps the previous flip. The `dGamma > 0` assertion stays the right guard either way. See the `manzaridafalias-melastflag-is-static-one` quirk.
