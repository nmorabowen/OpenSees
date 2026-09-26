# WP-126 — Reactions at partition-boundary nodes are wrong after stitching (OpenSeesMP)

Revision 1. Reproduction only; not adversarially reviewed. From WP-120 (#855) open question 1.

Status: **bug reproduced; draft PR. Fix design awaits the owner's decision (it spans the recorder, the
`.ladruno` schema and apeGmsh).**

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

## Proposed fix (owner's decision)

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

- Which of fix 3's two envelope behaviours (refuse vs flag `PARTIAL`)?
- Should `nodeReaction`-based recorders (vanilla `Node` recorder, MPCO) get a note in the user guide? They
  have the same partial-per-rank semantics; that is OpenSees behaviour, not a fork bug.
