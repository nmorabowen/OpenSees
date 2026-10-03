---
title: "WP-165 — Ladruno recorder stage semantics: topology-keyed MODEL_STAGE, persistent energy, empty partitions, RUN_ID"
project: Ladruno
type: correctness work package (follows WP-163, WP-164)
status: "in progress — draft PR"
owner: nmora
related:
  - "[[163_ladruno_recorder_hardening_roadmap]] (the review: R3, R6, MP-4, MP-8, MP-9)"
  - "[[164_ladruno_recorder_write_path]] (stage-end envelope finalize this builds on)"
  - "[[ladruno_schema_v1]] (INFO RUN_ID / RUN_ID_SCOPE, MODEL_STAGE EMPTY_PARTITION)"
  - "[[ladruno_apegmsh_contract]] (part-set validation)"
tags: [recorder, multi-stage, energy, partition, wp-165]
updated: 2026-10-03
---

# WP-165 — Ladruno recorder stage semantics

> [!summary] The short version
> Three of the review's findings had one root cause: the recorder treated every move of
> `Domain::hasDomainChanged()` as a new analysis stage. That stamp means "the DOF graph needs
> re-handling" and moves for SP patterns, `eleLoad`, contact `-reemit` re-sorts and more. WP-165
> keys the `MODEL_STAGE` on the actual node/element set, moves the energy integrals out of the
> per-stage source, makes a process with no nodes a valid empty partition, and stamps a run id.

## 1. Changes

| Finding | Before | After |
|---|---|---|
| R6 | every stamp move → new `MODEL_STAGE`: full model copy, all sources rebuilt, envelopes + energy reset. A `-reemit` contact run made one stage per re-sort | a stamp move rebuilds only when an FNV-1a fingerprint of the node / element / pressure-constraint set (tag + object address, so a re-used tag with a new object counts) changed; otherwise the stage, its sources, envelopes and energy continue |
| R3 | energy state lived in the per-stage `EnergyBalanceSource`: every new stage restarted IE/DW/ULW at 0 with a `rate × t_total` jump; integration only over `-T`-recorded samples | `EnergyState` owned by the recorder (survives the rebuild); `advance()` — idempotent per commit tag — called on EVERY commit before the `-T` gate; evaluate() just reads the latest values |
| MP-8 | a process with no nodes of the recorded set → "no nodes to write", error after the stage skeleton was written → a part file without `MODEL/NODES` | zero-length `NODES/ID` / `COORDINATES`, `MODEL_STAGE` attr `EMPTY_PARTITION = 1`, no node/element channels; `ON_DOMAIN` still recorded |
| MP-9 | nothing tied the part files of one run together | `INFO/RUN_ID` + `RUN_ID_SCOPE`: P0 broadcast (PartitionedDomain), `LADRUNO_RUN_ID`, the launcher job id (SLURM, OpenMPI), or a process-unique id marked `process` |
| MP-4 | (stage names diverge across ranks) | mostly removed with R6 — the common rank-local stamp movers no longer create stages; a genuinely rank-local topology change still can |
| docs | 03 §Parallel promised a stamp `MPI_Allreduce` that is compiled out; 03 schema sketch showed `DATA/STEP_k`; ResultIO.h contract said "stateless", "DATA/STEP_<k>", "already partition-reduced"; LadrunoRecorder.h banner "Phase-1 skeleton" | corrected |

## 2. Behaviour notes

- **Energy across a real topology change continues** (staged construction, ADR-51 element removal):
  the dissipated work stays in IE/DW. A user who wants per-stage energy subtracts the first row of the
  stage.
- **`EMPTY_PARTITION`** also covers a serial `-R` run whose region was removed entirely (the WP-163
  R2 scenario now writes an empty stage instead of suspending recording).
- **Cost of R6:** the fingerprint is O(nodes + elements) per stamp move, not per step — a contact run
  that re-sorts every step pays one domain walk per step instead of a full model rewrite.

## 3. Verification

`tests/test_ladruno_recorder_stage_semantics.py`:
- an SP pattern added mid-run keeps one `MODEL_STAGE` and continuous rows (was 2 stages);
- adding a node still starts a new stage;
- energy ULW continues across a topology stage change (one step of work, not a restart) and the closure holds;
- energy at common sample times is identical for `-T nsteps 1` and `-T nsteps 5`;
- `INFO/RUN_ID` honours `LADRUNO_RUN_ID` (scope `user`).

`tests/test_ladruno_recorder_hardening.py` R2 case updated: the removed region is now an `EMPTY_PARTITION`
stage (zero-length `NODES/ID`), no error.
