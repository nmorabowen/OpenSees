---
title: ADR-78 parallel contact — session handoff (2026-08-11)
project: Ladruno
status: handoff
owner: nmora
tags:
  - handoff
  - contact
  - parallel
  - mpi
---

# ADR-78 session handoff — 2026-08-11

Cross-library session. The engine half is **ADR-78** (this repo); the emit half
is **apeGmsh ADR 0092**. Neither ships alone, and the two now interlock in a way
worth understanding before touching either.

## What landed

| PR | what |
|---|---|
| fork **#726** | Contact verbs registered in the classic Tcl engine → `OpenSeesMP` can see them. **P0.5** |
| fork **#730** | **P1** — fifteen silent degradations become aborts; MPI job teardown |
| fork **#731** | P1 review fixes — SOFT narrowed, conversion finished, known limitation recorded |
| apeGmsh **#909** | ADR 0092 (the emit-side design) |
| apeGmsh **#910** | **S1** — owner/ghost resolver |
| apeGmsh **#911/#914** | **S2** — INV-4 master-uncuttable |
| apeGmsh **#916** | INV-1 corrected — unique-owner tally |
| apeGmsh **#917** | S2 empty-rank guard + the missing assertions |

**P0, P0.5, P1 done. S1, S2 done.** Next: **P2**, **S3**, then **S4**.

## The one thing to understand first

The deferral in ADR-39 assumed parallel contact needed OpenSees'
`PartitionedDomain`/`Subdomain` actor machinery, dynamic ghost nodes and an
in-Newton allreduce. **That premise does not apply to how this fork is driven.**

apeGmsh partitioned decks are *manual domain decomposition*: one plain `Domain`
per rank inside `if {[getPID] == K}` guards, boundary nodes replicated by tag,
global assembly by the distributed SOE. `DistributedDiagonalSolver::solve()`
sums **both** the diagonal mass `A` and the residual `B` over shared DOFs. So a
contact force assembled into a ghost DOF reaches the owning rank through the
solver, and the contact engine needs no MPI at all. That is why the whole design
is "one owner rank per interaction, ghost the interface" rather than a
distributed contact subsystem.

Measured, not argued: 2-rank vs serial within 1.6e−14 implicit, **bit-identical**
explicit. Harness preserved at `Ladruno_files/testbed/contact_parallel/`.

## Open items, in the order I would take them

1. **`AutoConstraintHandler`'s dead `MPI_Allreduce`.** Confirmed twice from the
   generated `build.ninja` `DEFINES`: the file is in the `OPS_Analysis` object
   library, is not listed per-target, so its `#ifdef _PARALLEL_*` block compiles
   out of **every** target. Under MPI `constraints Auto` sizes each rank's
   penalty from rank-local stiffness. Pre-existing, unrelated to contact, and now
   a one-line fix via the `OPS_CONTACT_PER_TARGET_SOURCES` pattern P1
   established. This is a live correctness bug in shipped code.
2. **P1 known limitation — runtime element removal.** `handle()` re-runs on every
   `domainChanged()`, and `LadrunoContactDomain` has no API to retire a contact
   or prune a surface, while `RemoveRecorder` removes nodes at runtime. A
   collapse run whose removed node belonged to a contact surface now aborts
   mid-run, and the MPI teardown skips `H5Fclose`. **Do not combine
   `recorder Collapse` with a declared contact** until surface pruning or a
   contact-removal command exists. Recorded in ADR-78; this blocks the fork's own
   progressive-collapse roadmap (ADR-51/54) and should probably outrank P2.
3. **P2** — `-soft` refusal under partitioning, `LadrunoContactDomain::sendSelf`
   / `recvSelf`.
4. **S3 / S4** on the apeGmsh side. S4 matters most: it is the layer that can
   supply **element→rank ownership**, which is what makes ADR 0092 INV-1 exact
   instead of a proxy that refuses on an undecidable tie (see below).

## Where the two ADRs interlock

ADR 0092 INV-1 picks the owner rank by master-node majority. That is a *proxy*
for "the rank holding the master surface's backing solid elements", which is what
the fork's `-kn auto` actually needs. The proxy is safe **only because a mis-pick
now fails loudly** — P1 made `-kn auto` failure fatal in all three lanes (NTS,
mortar, edge-edge). If anyone ever softens those aborts back to skips, the emit
side silently regresses. The dependency is not obvious from either file alone.

## Corrections made to this ADR's own claims

Recorded because the pattern matters more than the individual errors:

* **D5 described a stale tree.** Upstream had already added `ladrunoSurfaceNodesOk`
  (loud skip on a missing node); the coordinate-backfill path D5 cited was no
  longer the live one.
* **P1's first `MPI_Abort` was inert.** Written inside `LadrunoContactHandler.cpp`,
  which compiles into the `OPS_Analysis` object library — built once, no parallel
  define. It compiled clean, left both happy paths bit-identical, and did
  nothing. Only the mutation deck caught it.
* **P1's SOFT guard rejected correct models.** It scanned both surfaces and
  claimed the alternative was silence; `LadrunoContactFE.cpp:448` treats a
  massless master as a *fixed* one (supported), and `:516` already warned and
  fell back to the configured kn. Narrowed to slave nodes.
* **P1's verification was misleading.** `grep "; skipped"` returned zero because
  the pattern matched what had been *changed*, not what was *claimed*. Four paths
  survived, one with no message at all.

Every one of those was caught by running something, not by reading. The mutation
decks are the reason this ADR's numbers can be trusted; keep running them.

## Environment notes that cost time

* Build with `Ladruno_scripts\build.bat <targets>` — **not** the root
  `makeWIN.bat` (stale: lp64/static vs the siblings' ilp64/dynamic). `build.bat`
  self-heals a `CMakeCache.txt` written by a foreign cmake (Step 2c), which is
  the failure mode a hand-rolled `cmake .` produces.
* `setup_env.bat` deliberately avoids Intel's `setvars.bat` (it breaks on this
  machine) and handles the oneAPI-2026 Fortran split. Calling `setvars.bat`
  yourself loses MPI and silently drops the `OpenSeesMP` target.
* A hung MPI job leaves orphaned ranks that hold the binary open; the next link
  fails `LNK1104`. Check for stray processes before blaming the build.
* Python side: pin `LADRUNO_OPENSEES_BIN` / `LADRUNO_OPENSEESMP_BIN` and assert
  `mod.__file__`. A bare `import opensees` resolved to another session's
  scratchpad build.
