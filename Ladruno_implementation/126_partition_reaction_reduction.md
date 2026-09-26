# WP-126 — Reactions at partition-boundary nodes are wrong after stitching (OpenSeesMP)

Revision 1. Reproduction only; not adversarially reviewed. From WP-120 (#855) open question 1.

Status: **OpenSees side implemented and verified on the final build; PR #861 marked ready 2026-09-26. The
owner merges.** The owner chose "flag in file +
apeGmsh sums" (2026-09-25). The apeGmsh stitch fix runs in its own session, as a separate apeGmsh PR.

Scoped 2026-09-25. Branch `wp/126-partition-reaction-reduction`, cut from `ladruno` @ `bc5c33453`.
(WP-125 was taken by another session.)

## Problem

WP-120's dead-code scan found `ResultSource::requiresPartitionReduction()` (`SRC/recorder/Ladruno_ResultIO.h:82`,
8 definitions) with **no caller**. The recorder design (`03_ladruno_recorder.md`, D5/D6) explains why it
exists:
- each partition writes its own `part-N.ladruno` and "the reader stitches" (apeGmsh owns stitching);
- "additive" quantities (reactions at partition-boundary nodes, global sums, energy) need a per-step
  reduction instead. That was deferred to "v3b", and the flag marks which sources need it.

The flag was set, but nothing reads it: not the recorder, not the file, not the reader.

## Reproduction (2026-09-25)

`wp126_partition_reactions/repro_model.py`: two trusses into one **fixed support node shared by two
partitions**, statically determinate. Truss 1 (rank 0) is vertical with 10 kN on its free node. Truss 2
(rank 1) is at 45° with 20 kN. Two static load steps (50 %, 100 %). The run used `ladruno`'s binaries
(WP-123's worktree build; its diff touches no recorder or solver code). Serial = the same script on 1 rank.
`wp126_partition_reactions/repro_check.py` reads the results with h5py and apeGmsh (Python 3.11, editable
apeGmsh).

| Layer | Reaction at node 1, final step (Rx, Ry) |
|---|---|
| serial (np = 1): `nodeReaction`, `.ladruno` and apeGmsh `LadrunoReader` | **(20, 30)** |
| np = 2, rank 0 `nodeReaction` / `part-0` stored | (0, 10) |
| np = 2, rank 1 `nodeReaction` / `part-1` stored | (20, 20) |
| np = 2, **apeGmsh `LadrunoMultiPartitionReader` (what a user gets)** | **(0, 10)** — wrong |
| np = 2, `-envelope` MAX: `part-0` / `part-1` | (0, 10) / (20, 20); the true max is (20, 30) and cannot be recovered |

- **The solver is right.** Each rank's reaction is the partial from its own elements
  (`Domain::calculateNodalReactions`), and the partials sum exactly to the serial value.
- **The stitch is wrong.** apeGmsh's `_merge_node_slabs` (`apeGmsh/results/readers/_mpco_multi.py:392`,
  shared by the MPCO and Ladruno multi-partition readers) keeps the **first partition's** value for every
  nodal component. Its docstring: "Boundary nodes appear in multiple partitions with identical kinematics…
  First partition that has the node wins." That is right for displacements and wrong for reactions.
  - A support on a partition interface reports one partition's share.
  - Base shear summed from stitched reactions is too low.
  - A *free* interface node, whose reaction should be ≈ 0, reports a nonzero partial.
- **Envelopes are wrong and cannot be recovered.** `-envelope` gives every node channel an `EnvelopeSink`
  with no check of the flag (`LadrunoRecorder.cpp:1561`). A per-partition max of a partial reaction cannot
  be recombined, because max(a+b) ≠ max a + max b. There is no warning.
- **Scope.** Any OpenSeesMP / openseesmp run that records reactions with the Ladruno recorder and reads them
  through apeGmsh, whenever a supported (or loaded) node lies on a partition interface. The same stitch
  serves MPCO part files. The energy sources (`EnergyBalanceSource`, flag `true`) are already documented as
  a v3b stub.

## Implemented (2026-09-26)

- **`ResultSource::partitionReduction()`** (`Ladruno_ResultIO.h`) returns `"SUM"` when
  `requiresPartitionReduction()`, else `"NONE"`. `EnergyBalanceSource` overrides it to `"UNSUPPORTED"`,
  because a sum would be wrong: shared-node KE is counted in each partition, and RES/ERR are derived. So
  the three-way value, not the old bool, is what reaches the file.
- **Sinks** (`Ladruno_Sinks.cpp`) write `PARTITION_REDUCTION` on every streaming result group
  (`StreamingSink::begin`) and on every envelope group (`EnvelopeSink`).
  - The envelope value is cached in `accept()` as well as `begin()`: the recorder drives an EnvelopeSink
    through `accept()` only. The first build wrote no envelope attribute at all, because the HDF5 string
    writer silently skips an empty string.
- **Recorder** (`LadrunoRecorder.cpp`):
  - `is_partitioned` is kept from `initialize()`.
  - `refusePartitionedEnvelope()` skips, with a warning, any node or domain channel whose reduction is not
    `NONE` when `-envelope` runs partitioned. Element sources are always `NONE`.
  - A one-time warning says partitioned `energyBalance` output is not mergeable.
- **Schema / contract / validator:**
  - schema §7.1 table, plus the stitching note in the naming section;
  - `ladruno_apegmsh_contract.md` bullet;
  - `ladruno_format._validate_partition_reduction`: optional; when present it must be NONE|SUM|UNSUPPORTED;
    a partitioned envelope must be NONE.
- **apeGmsh:** a separate session in the apeGmsh repo (started by the owner from the task chip). It has the
  contract and the real part files below.

### Results

| Check | Result |
|---|---|
| Build (5 targets) | 0 errors; 0 warnings in the touched files |
| `ci`-gated D5 harness (`test_ladruno_harness.py`, incl. new `test_validator_partition_reduction`) | 27/27 |
| `tests/test_ladruno_partition_reduction.py` (zone_a; skips without h5py) | 4/4: DISPLACEMENT=NONE, REACTION_FORCE/UNBALANCED_FORCE=SUM, energy=UNSUPPORTED, serial SUM envelope kept, partitioned reaction envelope refused (single process, `PMI_SIZE=2` subprocess) |
| Two-rank gate `mp_reaction_model.py` + `mp_reaction_check.py` | **ALL PASS**: both parts flag REACTION_FORCE=SUM / DISPLACEMENT=NONE; contract stitch at the shared support = **(20, 30) = serial** (first-wins would give (0, 10)); displacement unchanged; partitioned reaction envelope refused on both ranks with the warning; serial envelope kept |
| Validator on all six real output files | valid |
| **Final full 5-target build**: live recorder regression battery `run_regression.bat` (every gate vs frozen MPCO: nodal/element parity 1e-12, multi-stage, envelopes, local axes, energy, Bezier, truss/zeroLength, frame3D, shell, eigen, f32 precision, rank-env probe, Tcl flag order) | **ALL GATES PASSED**: the attribute changes no recorded value |
| Final build: `mp_parallel` gate (3 ranks, openseesmp) | ALL PASS |
| Final build: new pytest, two-rank gate, D5 harness | 4/4, ALL PASS, 27/27 |

## Proposed fix (owner's decision, 2026-09-25: this option)

1. **Put the flag in the file.** Each result group gets `PARTITION_REDUCTION = "SUM" | "NONE"`, written
   from `requiresPartitionReduction()`. This is a small schema addition; the validator and the apeGmsh
   contract doc `ladruno_apegmsh_contract.md` both live in this repo.
2. **apeGmsh stitch honours it.** Components flagged `SUM` are summed over every partition that holds the
   node; the others keep first-wins. Older files and MPCO have no attribute, so they fall back to the
   component name (`reaction*`). This is an apeGmsh PR under its own conventions.
3. **Envelopes.** In a partitioned run, a source flagged `SUM` cannot be enveloped correctly without a
   per-step `Allreduce`, which D5 keeps out of the recorder. Refuse that channel with a clear warning, or
   write it flagged `PARTIAL` so the reader refuses it. Real reduced envelopes stay v3b.
4. **Tests.** This reproduction becomes a zone gate. Both the openseesmp run and apeGmsh are needed, so it
   likely lives with the MP recorder gate (`Ladruno_scripts/ladruno_recorder_tests/mp_parallel_*`).
   Assertions:
   - stitched == serial for reactions;
   - the `PARTITION_REDUCTION` attributes are present;
   - an enveloped reaction in a partitioned run warns.

## Rejected approaches

- **Reduce in the recorder (per-step `Allreduce`).** Would make stored values right for every reader. But
  D5 keeps the recorder MPI-free, since MS-MPI/libfabric is fragile and it fights SWMR, and it would
  duplicate shared-node data. Kept as the v3b route for envelopes only.
- **Name-based summing only (sum anything called `reaction*` in apeGmsh).** Quick, but it fails for any
  future additive quantity and puts recorder semantics in the reader. Kept only as the fallback for files
  without the attribute.

## Open questions

- ~~Which of fix 3's two envelope behaviours?~~ Refuse, with a warning (owner, 2026-09-25).
- **`CONSTRAINT_TIE_FORCE`** (`ConstraintTieForceSource`) is flagged `NONE`. Whether a constraint force
  at a node shared by two partitions is also a per-partition partial is **unverified**. Settle it with a
  two-rank run like `mp_reaction_*` before anyone relies on stitched tie forces.
- apeGmsh side: nmorabowen/apeGmsh#1179 (draft). It sums SUM rows (vectorised), raises on UNSUPPORTED and
  on partitions that disagree on the kind, raises on partitioned `read_energy`, and falls back by name
  (`REACTION*`/`UNBALANCED*`/`RAYLEIGH*`) for files without the attribute. It was verified end-to-end on
  this WP's real part files: stitched (20, 30) = serial, and the same via the name fallback on a
  pre-WP-126 build.
- Should `nodeReaction`-based recorders (vanilla `Node` recorder, MPCO) get a note in the user guide? They
  have the same partial-per-rank semantics; that is OpenSees behaviour, not a fork bug.
