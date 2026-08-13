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

1. ~~**`AutoConstraintHandler`'s dead `MPI_Allreduce`.**~~ **DONE 2026-08-11.**
   Fixed via the per-target-TU pattern: the reduction now lives in
   `LadrunoAutoPenaltyReduce.cpp`, listed in `OPS_MPI_PER_TARGET_SOURCES` (the
   P1 list, renamed since it now holds two unrelated files).
   `AutoConstraintHandler.cpp` is left with **no `#ifdef _PARALLEL_*` at all** —
   that, not the guard's contents, is what stops the trap reopening. Gated by
   `Ladruno_files/testbed/auto_penalty_mpi/`, whose pre-fix signature was
   **measured by mutation** ( `[1e9, 1e5]` instead of `1e7` twice ), not
   predicted. See ADR-78 §FOLLOW-UP.
2. ~~**P1 known limitation — runtime element removal.**~~ **DONE 2026-08-12
   (#737).** `recorder Collapse` + a declared contact is legal again. The knot was
   that a deck typo and a node the ANALYSIS removed were indistinguishable — both
   first surfaced inside `handle()`, because `LadrunoContactSurface` never saw the
   `Domain` and `addSurface()` did no existence check. Split them by WHEN the tag
   went missing: absent at declaration ⇒ typo ⇒ refused at the deck line; present
   then absent at `handle()` ⇒ removal ⇒ pruned, run continues. An emptied surface
   RETIRES its interaction with a named notice. The `-kn auto` / `-soft` /
   `kn <= 0` aborts are untouched and now have a guard-on-the-guard test, because
   ADR 0092 INV-1 leans on those. See ADR-78 §REMOVAL LANE.
3. **P2** — `-soft` refusal under partitioning, `LadrunoContactDomain::sendSelf`
   / `recvSelf`.
4. ~~**P1 left four tests red on `ladruno`.**~~ **DONE** — re-greened per test by
   asking whether the broken precondition was the SUBJECT (invert) or incidental
   (repair the deck); one of them turned out to be passing *because of* the bug
   it should have caught. All four mutation-verified. See ADR-78 §P1 FALLOUT.
5. **S3 / S4** on the apeGmsh side. S4 matters most: it is the layer that can
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

---

# Session 2 — 2026-08-12

Continued from the 2026-08-11 session above. Four PRs, all merged to `ladruno`
(tip `900501d71`). Three of the four were **not** on the original roadmap: they
were found by measuring, and each one explains why the next was findable.

## What landed

| PR | what |
|---|---|
| **#733** | `constraints Auto` was sizing its MPI penalty from ONE rank's elements — open item 1 |
| **#734** | The four tests P1 orphaned + the collection abort that made the suite unrunnable + the ASDConcrete3D defect report |
| **#736** | `cross-tier-nightly` — the whole suite in ONE process, so cross-test state leaks are catchable at all |
| **#737** | The removal lane — open item 2, which the previous handoff ranked above P2 |

## The chain, because the order is the point

1. **#733** fixed a dead `MPI_Allreduce` in `AutoConstraintHandler`. Same
   OPS_Analysis object-library trap P1 hit; found by taking the previous
   handoff's open item 1 literally.
2. Regression-checking #733 turned up **four tests red on `ladruno`** — P1 had
   changed a contract and touched no test file. Fixing them (#734) meant deciding
   per test whether the broken precondition WAS the subject (invert) or
   incidental (repair the deck).
3. One of those four turned out to be **passing because of the bug it should
   have caught**, and then **my first repair of it was unfalsifiable**. Both
   caught by mutation, neither by review.
4. Trying to run the full suite to check for more revealed it **could not run at
   all** — a missing matplotlib aborted collection, taking ~2000 tests with it.
   That is why the four orphans had gone unnoticed.
5. With the suite finally running end to end: exactly **one** failure, and it was
   ORDER-DEPENDENT. Bisected in six steps to a tag-keyed process-global cache in
   upstream `ASDConcrete3D`. Reported, deliberately not fixed.
6. That defect was invisible to CI **by construction** (`-m zone_a` on PRs,
   `-m zone_b` nightly, never the same process) — hence **#736**.
7. **#737** then took the roadmap's own next item.

## Open items now

**P2 is next**, then S3/S4, then P4 — the list above is current. Plus:

* **ASDConcrete3D tag cache — parked, reported, sentinel armed.** `ladruno`
  carries ONE known-failing test on purpose:
  `test_ladrunoBrick_asdconcrete_bend.py::test_notched_bend_mesh_objectivity`,
  which fails in a full-suite run and passes alone. Do NOT repair, retag, or
  xfail it — it is the only thing detecting a live silent-wrong-answer. Evidence,
  reproducer, the two candidate fixes and the gate they must satisfy:
  `Ladruno_files/testbed/asdconcrete_tag_cache/`. Also ADR-11 §KNOWN DEFECT and
  [[LEDGER_quirks]].
* **`cross-tier-nightly` fires for the first time tonight (06:00 UTC).** Two
  things to know cold: it takes ~50 min, and its **sentinel step REQUIRES the
  ASDConcrete3D defect to still reproduce**. If someone fixes that material the
  sentinel goes green and the job goes RED on purpose, telling you to delete the
  now-unearned `--deselect`. Also: it needs the self-hosted runner online — when
  the runner is offline the job queues and cancels, which is SILENT, so confirm
  the first run actually executed rather than assuming green.

## Traps this session paid for — all recorded in [[LEDGER_quirks]]

Listed here only as an index; the entries carry the detail and the measurements.

* A test can be GREEN because of the very bug it should catch — and the fix for
  it can be unfalsifiable in a new way. Mutate every "X changes nothing"
  assertion.
* Mutation hygiene: anchor on a UNIQUE string (a duplicated anchor edited an
  unrelated test and reported the target as vacuous); restore by writing back the
  bytes you read, NEVER `git checkout -- <path>` (it restores from the INDEX and
  destroyed a finished rewrite after a `git stash pop` de-staged it).
* A behaviour-changing PR that touches no test file leaves the suite asserting
  the old contract, and the failures read as generic breakage.
* A missing OPTIONAL dependency does not skip two tests — it aborts the ENTIRE
  suite at collection.
* Marker-filtered CI tiers can never catch a cross-test state leak.
* `--deselect` with an unmatched path is SILENTLY ignored — verify by collection
  count, with a bogus path as the control.
* `pytest -k` applies to EVERY argument, so it deselects the target you are
  measuring and reports a clean 0.08 s pass.
* Peak load cannot discriminate two softening backbones (`FT*area` either way) —
  integrate dissipated work.
* **Two PRs can merge CLEANLY into a contradictory document.** #733 and #734 both
  edited this file's open-items list; git conflicted on the ADR (both appended a
  section — resolution keeps BOTH) but auto-merged the handoff into two item-5s,
  one of which described work the same file recorded as DONE. Read the merged
  result; a clean merge is not a correct merge.

## Environment notes

* **Run the suite with `opensees_env`**, not base `python3.12`. It has gmsh 4.15,
  matplotlib, h5py, scipy, apeGmsh. Base 3.12 lacks gmsh, so `conftest.py` skips
  the **entire 129-test `zone_b` tier** by design and the run looks green while a
  whole tier never executed. That tier is where the material defect was hiding.
* **`LADRUNO_OPENSEES_BIN` wants the DIRECTORY**, not the `.pyd` path — pointing
  it at the file makes `import opensees` fail with `ModuleNotFoundError`, which
  reads as a broken build. And unpinned in `opensees_env`, `import opensees`
  still resolves to `C:\Program Files\Ladruno\...`, so always assert
  `opensees.__file__`.
* Removing a node that still carries a `NodalLoad` ACCESS-VIOLATES on the next
  `analyze` — measured with `constraints Plain` and no contact, so it is a stock
  lifecycle hazard, not a contact one. Collapse decks must drop the load first.

---

# Session 3 — 2026-08-12

One PR: **P2**, the roadmap item both previous handoffs pointed at. Both halves
shipped and gated by `Ladruno_files/testbed/contact_p2/run.py` (7 cases; the
pre-P2 installed build fails exactly the 3 instruments it should — the best
negative control this lane has had, and free).

## What landed

* **`-soft` refused under partitioning, ALL four lanes** (NTS / rigid plane /
  mortar SOFT=2 / `-edgeSoft`), one named FATAL at the `anySoft` choke point.
  np detection lives in `LadrunoContactAbort.cpp` (`ladrunoContactNumRanks()`),
  NOT the handler — the P1 inert-`#ifdef` trap again, respected this time from
  the start. np=1 (serial and `mpiexec -n 1`) unaffected, mutation-verified.
* **`LadrunoContactDomain::sendSelf`/`recvSelf`** — definitions-only, through
  `Domain::sendSelf/recvSelf` (`domainData` 17→19). `database File`
  save→wipe→restore now carries contact exactly; pre-P2 it restored a
  contact-free model SILENTLY (measured: top block floating, rt=0.0).
* Partitioned `-kn auto` re-measured vs its serial twin on the P2 build:
  1.6e−14 (the P0 number). Serial contact battery 147/147.

## Findings worth knowing cold

1. **P1's mass guard was ONLY ever protecting the fully-ghosted NTS case.**
   Ghost ⇒ zero rank-local mass ⇒ P1 fires. A partition-boundary node has
   PARTIAL mass (nonzero — locally undetectable), and the mortar/edge/plane
   soft lanes were never scanned at all: the pre-P2 build runs a 2-rank mortar
   `-soft` deck to completion, silently. That mortar deck is the honest
   mutation for the refusal; the NTS deck alone can pass via the WRONG abort,
   which is why both decks carry `rho > 0`.
2. **The refusal over-refuses, disclosed:** an MP parametric sweep where every
   rank holds the full model would size SOFT correctly and is refused anyway —
   indistinguishable from manual DD from inside `handle()`. The error message
   says so and names the workaround (explicit `-kn` / serial / a future
   deck-supplied `m_eff`).
3. **The DB stream is now build-lineage-scoped** (17→19): a pre-P2 database
   cannot be restored by a P2+ build and the error names nothing useful. In
   [[LEDGER_quirks]].

## Open items now

**S3/S4 (apeGmsh side) are next**, then **P4** (the validation battery — note
its `-soft` case is now moot-by-refusal on partitioned decks; pounding on 2
ranks should run `-kn auto`/explicit kn). The ASDConcrete3D sentinel and
`cross-tier-nightly` notes from Session 2 stand unchanged.

---

# Session 4 — 2026-08-12 (P4 validation battery)

One PR: **P4**, the fork-side validation battery — the last committed fork phase
of this ADR. All 12 verdicts PASS against a fresh `cebec3f9f` build (hash checked
against the worktree HEAD via the splash before any measurement). Battery at
`Ladruno_files/testbed/contact_p4/` (`python run.py <build-dir>`).

## What landed (measured, headline numbers)

* **Explicit 2-body pounding (NTS `-kn auto`) is BIT-IDENTICAL to serial at
  np=2 AND np=4** — full disp/vel histories at 17 recorded digits, 3000 steps,
  3 impact/separation cycles, ghost == native every step, identical peak base
  reaction. The np=4 deck cuts BOTH bodies across ranks (regular DD boundary +
  contact interface in one model — first time measured together).
* **Contact energy balance closes:** 0.002% of E_peak at contact-open samples
  (serial AND partitioned), per-rank EnergyBalance channels sum to serial at
  5.6e-16 / 6.0e-16 of E_peak. The in-contact RES excursion (up to 66% of
  E_peak) is the penalty-spring storage — contact assembles at the SOE, which
  the recorder sweep never sees — so the closure gate must sit on separated
  samples. Physical, not an error; returned on separation.
* **Mortar mesh-tie across ranks (implicit Mumps + held-load ALM):** tie
  residual 2.2e-19 (np2) / 3.6e-19 (serial) after 15 fixed augmentations; np2
  vs serial 6.5e-16; tip == analytic 2P/E to 10 digits; loaded in TENSION so a
  silently-missing tie cannot pass.
* **np=1/2/4 sanity (NB=8, 40k steps):** 9.1 / 4.9 / 3.2 s, answers
  bit-identical across np. Desktop, startup-inclusive — a sanity check, NOT a
  scaling study (Q-P0GATHER is unmeasurable at this size; it belongs to the
  cluster campaign that would justify P5).
* **Controls:** P2 soft-refusal re-pinned on the pounding deck; INV-6
  ghost-mass mutation runs silently and moves the answer 91% — the battery's
  own serial-comparison instrument catches it.

## The finding that outranks the numbers

**#737 regressed P1's MPI teardown for the missing-ghost case — OPEN DEFECT.**
The removal lane moved missing-node detection from `handle()` (P1's
`ladrunoContactFatal`, measured 1.1 s teardown) to `contactSurface` declaration
time, where the refusal is a rank-local parser `return -1` with no MPI
awareness. The refusing rank aborts its script; the other rank blocks in the
next collective. The preserved P0 gate `contact_parallel/mp_noghost.tcl` now
HANGS on every build since #737 (measured: >30 s, 45 s, and one 17 s
hydra-auto-cleanup reap — nondeterministic). Refusal text is still loud and
named; only the teardown is missing. Fix shape: route the declaration refusal
through `ladrunoContactNumRanks()`/`ladrunoContactFatal()` when np > 1 — the
per-target TU exists for exactly this. Recorded in ADR §P4 finding 1 +
[[LEDGER_quirks]]; NOT fixed in the P4 PR (measurement lane, no SRC changes).

## Traps paid for (details in [[LEDGER_quirks]])

* Orphaned MPI ranks inherit the harness's stdout pipe:
  `subprocess.run(capture_output=True)` blocks forever even AFTER its own
  TimeoutExpired killed mpiexec. File-redirect + `taskkill /F /T` is the only
  robust shape; `run.py`'s `Runner.run` is the reference.
* Recorder files default to `-precision 6`: a serial-vs-parallel parity gate
  reads ~3e-7 from quantization alone and proves nothing tighter. `-precision
  17` on every recorder feeding a gate; under MPI filenames grow `.part-<rank>`.
* ALM across ranks: the held-load augment loop must be FIXED-COUNT on every
  rank (`analyze` is collective, the residual query is rank-local) — a
  residual-keyed break deadlocks. `analyze_augmented`'s while-loop shape does
  not port to partitioned decks as-is.

## Open items now

1. **The #737 teardown gap** (above) — small SRC fix + re-green
   `mp_noghost.tcl` + a P4 `ctl-noghost` gate flip (its verdict line documents
   the expected post-fix behaviour).
2. **apeGmsh S6** (its ADR 0092) if still open, then the ADR-78 table is done
   through P4; **P5 stays deferred** and now has a written measurement
   precondition (cluster-scale interface, §P4 last paragraph).
3. Unchanged from Session 2: ASDConcrete3D sentinel, `cross-tier-nightly`.
