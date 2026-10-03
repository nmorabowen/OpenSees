---
title: "WP-163 — Ladruno recorder: deep review (red/blue) + hardening; roadmap WP-163/164/165"
project: Ladruno
type: review + work-package plan
status: "draft — WP-163 open; WP-164 (write path) and WP-165 (stage semantics) queued, sequential (same sink code)"
owner: nmora
related:
  - "[[03_ladruno_recorder]] (recorder design; stale spots listed in section 6)"
  - "[[126_partition_reaction_reduction]] (PARTITION_REDUCTION contract)"
  - "[[69_ladruno_energy_recorder_channels_adr]] (energy kernel)"
  - "[[08_analysis_monitor]] (Monitor recorder)"
tags: [recorder, hdf5, robustness, performance, scale, partition, wp-163, wp-164, wp-165]
updated: 2026-10-03
---

# WP-163 — Ladruno recorder: review and hardening roadmap

> [!summary] The short version
> A red-team / blue-team code review (2026-10-03; 5 adversarial reviewers, 5 refuting
> verifiers, read-only) of `SRC/recorder/Ladruno*` + `LadrunoMonitor*` + `EnergyBalance*`
> across architecture, performance, robustness, scale and multi-partition HDF5.
> The design holds (source/sink split, chunked `[T×N×C]` DATA, PARTITION_REDUCTION). The
> defects cluster in three places, and each gets one work package:
>
> | WP | Theme | Size | Changes speed? |
> |---|---|---|---|
> | **WP-163** (this) | Hardening: crashes, silent data loss, parser holes | 1–2 d, all S | no |
> | **WP-164** | Write path: handles, flush cadence, chunk tiling, `-compress`, in-place envelope | 3–5 d | yes (the gains) |
> | **WP-165** | Stage semantics: rebuild on real topology change, energy across stages, empty partitions, Monitor in MP | M | yes for `-reemit` |
>
> Sequential, not parallel: WP-164 rewrites the same sink code WP-163 hardens.

## 1. Findings that survived the blue team

Severities are the corrected (post-refutation) ones. Line numbers are at `fd4ff5983`.

### 1.1 Must fix — wrong results or crash (WP-163 unless noted)

| ID | Defect | Evidence | Fix |
|---|---|---|---|
| R1 | `eigen` → `wipe` → new model with `-N modesOfVibration` → `analyze` **exits the process** (`exit(-1)`): `wipe` never resets `numEigen` (Tcl global, `OpenSeesCommands`), the gate trusts it | `LadrunoRecorder.cpp:1970-1978` → `Domain.cpp:2812`; `Node.cpp:1447` `exit(0)` next | gate on `domain->getNumEigenvalues()`; per node `getNumEigenvectors() > mode` |
| R2 | **Use-after-free** when a stage rebuild fails: stamp stored before `writeModel()`; early returns skip `clearSources()`; stale element `Response*` used next step | `:421-428`, `:660-678` | clear sources before the writers; stage-invalid latch on failure |
| R4 | Duplicate `-N`/`-E` request (or alias `tieForce`/`constraintTieForce`) → **2T rows** in one dataset; 2nd sink's failed `begin()` still marks initialized | parser `:2744-2829`; `Ladruno_Sinks.cpp:126-162` | dedupe at parse; init only on success |
| R5 | **Silent result loss**: no return code checked in `createTimeSeries3d`; failed create → every `accept()` returns quietly; `appendSlab3d` rc ignored → TIME/STEP outrun DATA | `Ladruno_Hdf5.h:331-366`; `Ladruno_Sinks.cpp:150-212` | check rcs, one `opserr`, dead sink; skip TIME/STEP on failed slab |
| M1 | OpenSeesSP + `-G energy` **segfaults P0**: `sweepDomain` walks ShadowSubdomains, `getNodePtrs()==0` (same in standalone `EnergyBalanceRecorder`) | `EnergyBalanceKernel.h:438-445,104`; `EnergyBalanceRecorder.cpp:384,500` | `isSubdomain()` skip |
| R3 | Energy restarts at each MODEL_STAGE with a `rate×t` jump; integrated only on `-T`-gated samples | `Ladruno_DomainResults.cpp:74,125`; `EnergyBalanceKernel.h:221-236` | **WP-165** |
| R6 | `-reemit` contact → full MODEL_STAGE (model copy, envelope + energy reset) per re-sort | `Domain.cpp:2316` → `:419-431` | **WP-165** |

### 1.2 Robustness, smaller (WP-163)

| ID | Defect | Fix |
|---|---|---|
| ROB-8 | translational modes read eigvec rows without `noRows()` / `noCols()` check (`Ladruno_NodeResults.cpp:788-795`) | clamp rows, skip mode ≥ cols |
| ROB-9 | Monitor `-dof 0`/negative → `(*r)(-1)` every frame (`LadrunoMonitorRecorder.cpp:442,319`) | reject `d<1` at parse |
| ROB-11 / MP-10 | failed `initialize()` retried every step, leaks 2 proplists each time (`:396-399,480-485,548-554`) | `init_failed` latch, close proplists |
| ROB-12 | `-G` swallows the next token even if it is an option (`:2605-2619`) | un-consume `-`-prefixed token |
| ROB-13 | `-T dt` no tolerance, no snapping (`:383-390`) | `t - next_t >= -1e-5·dt`, `next_t += dt` |
| ROB-7 | envelope ignores NaN (first NaN sticks; later NaN masked) (`Ladruno_Sinks.cpp:280-301`) | propagate NaN deliberately |
| ROB-10 | size-mismatch rows zero-filled, indistinguishable from real zeros (`Ladruno_ElementResults.cpp:99-116`) | NaN-fill |
| ARCH-4 | `H5Tcopy(H5T_C_S1)` never closed (`Ladruno_Hdf5.h:80,226`) — unbounded with envelope | `H5Tclose` |
| ARCH-10 | destructor `finalizeAllSinks()` outside `H5E_BEGIN_TRY` (`:262-281`) | move inside / `H5Iis_valid` |
| M4 / ROB-3 | `SLURM_NTASKS` under `sbatch` without `srun` → sequential run written as `part-0` of N (`:508-527`; same in `EnergyBalanceRecorder.cpp:207`) | accept SLURM pair only with `SLURM_STEP_ID` |
| MP-3 | SIZE>1 with RANK missing / out of range → every rank `part-0`, TRUNC clobber (`:519-521`) | refuse |
| MP-6 / ARCH-11 | Monitor recorder: every openseesmp rank TRUNCs the same sink file | WP-163: refuse when SIZE>1; WP-165: `.part-N` |

### 1.3 Performance / scale (WP-164)

| ID | Defect | Estimated impact (unmeasured) |
|---|---|---|
| P1 | `-envelope` deletes + recreates every envelope group (uncompressed `nIds×nComp`) + COLUMN_MAP **every recorded step**; doc says "periodically" | ~3.5× the I/O of streaming; file growth |
| P2 | DATA/TIME/STEP reopened + `H5Fflush` every step → partial compressed chunk recompressed ~`ct` times (slab < 128 KiB) | ~1–3 ms/channel/step; long explicit runs |
| P3 | chunk `{ct, nIds, nComp}`: single-entity history inflates the whole dataset; slab-sized RAM spikes; hard cap 2²⁹ f64 per partition | read side 100–1000× |
| P4 | deflate-4 hard-coded on the solver thread | ~30–60 ms per 2.4 MB slab |
| P5 | `Domain::getNode(tag)` per node per channel per step though `Node*` is cached | 10–40 ms/step at 1e5 nodes |
| P7 | reaction channels unsorted; non-reaction channel resets the flag → repeated `calculateNodalReactions` | one residual sweep per extra |
| P8 | `writeSections()` fiber probe on every element every stage, requested or not | minutes for fiber-heavy staged models |
| — | apeGmsh reader `np.asarray(DATA[...])` (`_ladruno.py:424,510`, `_ladruno_element_io.py`) | apeGmsh-side; OOMs first |

### 1.4 Stage semantics / partitions (WP-165)

R3, R6, MP-8 (rank with zero region nodes leaves a broken part file — write a valid empty
partition), MP-4 (stage ids drift across ranks after a rank-local event), MP-9 (`RUN_ID` to
reject stale part files), Monitor `.part-N`.

### 1.5 Refuted (recorded so nobody re-raises them)

FORMAT_VERSION bump (D3 landed before the freeze); SP double commit (DDA commit is a no-op);
boundary-node mass double count (external copies are massless); tie-force reduction (refused
in coupled parallel); partition-counter drift (re-partition unreachable); per-step buffer
allocation; header-compare cost; 64 KiB attribute limit (block merge collapses fibers);
gravity → `loadConst` → nodal pattern making 2 stages (only `eleLoad`/`sp`/`imposedMotion`/
topology edits bump the stamp); envelope delete/recreate "crash window" (order is correct).

## 2. WP-163 scope and acceptance

Each item lands with a regression deck under `Ladruno_scripts/ladruno_recorder_tests/` (or a
zone_a pytest under `tests/`), written to fail on `fd4ff5983` and pass after.

| Item | Regression |
|---|---|
| R1 | `eigen 3; wipe; model; -N modesOfVibration; analyze 1` survives, no MODE datasets (Tcl + py) |
| R2 | `-R` region, `-E stress`; remove all region elements + nodes; analyze → clean error, no crash |
| R4 | `-N displacement displacement` and `-N tieForce constraintTieForce` → DATA shape[0] == T |
| R5 | forced create failure (name clash / bad chunk) → one `opserr`, no silent success |
| M1 | unit-level: `isSubdomain()` skip (SP run not available in CI) |
| ROB-9/12/13, M4, MP-3 | parser / env decks (`rank_env_model.py` pattern) |
| ARCH-4 | envelope run, 500 steps: open-id count flat (`H5Fget_obj_count`) |

Ledgers: `LEDGER_implementations.md` row (WP-163), `LEDGER_quirks.md` entries (eigen gate after
`wipe`; HDF5 type-handle lifetime; SLURM batch env; METIS correction below). No upstream file is
touched except possibly `Domain`/`Node` (not planned — the gates are recorder-side).

## 3. WP-164 step 0 — baseline benchmark (before any write-path change)

Three decks, recorder on vs off, wall time and file size, on `fd4ff5983 + WP-163`:

1. explicit, small channels (energy + a 10-node set), 1e5 steps;
2. medium solid (~1e5 nodes) with `-envelope`;
3. large slab (stress at GPs, ~1e5 bricks), 200 steps.

They become the WP-164 regression gate (speed must improve, outputs bit-identical or
validator-equal).

## 4. Doc / ledger corrections found by the review

- `LEDGER_quirks.md:604-610`: the METIS-4 block is only `OPS_partition()` under
  `_PARALLEL_INTERPRETERS`; **OpenSeesSP.exe (Tcl) partitions automatically via the METIS-4
  legacy API**, so the PartitionedDomain recorder path ships (untested).
- `03_ladruno_recorder.md:331-333` promises a stamp `MPI_Allreduce` the code does not have;
  `:268-277` still shows `DATA/STEP_k`; `:354` says the envelope rewrite is "periodic".
- `Ladruno_ResultIO.h:56-60,101,116`: "stateless", "DATA/STEP_<k>", "already partition-reduced"
  are wrong; `reset()` is dead code.
- `LadrunoRecorder.h:1-8` banner ("Phase-1 skeleton"); `SOLVER_VERSION` hard-coded 3.5.1.
- `Ladruno_DomainResults.cpp:122-124` "takes no spurious increment" is false (R3).
