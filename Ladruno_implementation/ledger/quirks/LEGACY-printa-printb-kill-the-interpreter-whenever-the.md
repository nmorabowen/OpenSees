---
wp: LEGACY
title: "printA / printB KILL the interpreter whenever the SOE was never sized — and constraints LadrunoContact is what makes that reachable"
legacy_seq: 163
---
### `printA` / `printB` KILL the interpreter whenever the SOE was never sized — and `constraints LadrunoContact` is what makes that reachable

- **Bites:** any `ops.printA("-ret")` / `ops.printB("-ret")` under `system FullGeneral`
  or `system Diagonal` when the SOE has not been sized -- i.e. before the first
  successful `analyze()`, or after ANY `analyze()` that failed in
  `domainChanged()`. The process dies with exit code -1/255, **no Python
  traceback, no exception, and nothing printed**: `exit(-1)` is a clean exit so
  `faulthandler` is silent, and the `FATAL` text goes to `opserr`, which the pyd
  redirects. It looks exactly like a bug in your own script (reported 2026-08-18
  while debugging ADR-85 T1b).
- **Why it looks like a CONTACT bug, and why it is not:** `printA` ends in
  `LinearSOE::getA()`. Only **two** classes in the entire
  `SRC/system_of_eqn/linearSOE` tree override that virtual -- `FullGenLinSOE` and
  `DiagonalSOE` -- and both used to `exit(-1)` on a null `matA`. Everything else
  inherits `LinearSOE::getA()`'s `return 0` and was always safe (`UmfPack`,
  `BandGeneral`, ... return an empty result). `matA` is allocated in `setSize()`,
  which `StaticAnalysis::domainChanged()` reaches **only after
  `ConstraintHandler::handle()` returns `>= 0`**. Under the upstream handlers
  `handle()` essentially never fails, so nobody hit this. But the contact
  subsystem's whole ADR-78/ADR-85 abort discipline is BUILT on `handle()`
  returning -1 through `ladrunoContactFatal()` for a refused contact -- a failed
  `domainChanged()` is a *designed, routine* outcome there. So: refused contact ->
  `handle()` -1 -> `domainChanged()` fails -> `analyze()` != 0 -> SOE unsized ->
  and the very next thing you reach for to debug the refusal is `printA`, which
  takes the interpreter with it. Hence "printA crashes under LadrunoContact" on a
  deck where contact is not even the thing being measured.
- **Corollary worth internalizing:** after a nonzero `analyze()` under
  `constraints LadrunoContact`, the SOE holds NO tangent for your model. A
  `printA` there was never going to answer your question even when it did not
  crash -- read the FATAL/refusal line above it instead.
- **`printA` had a SECOND, unrelated crash on exactly one system**, found by
  sweeping all ten `system` choices instead of trusting the getA survey:
  `system SparseSYM` died with a real **ACCESS VIOLATION (0xC0000005)**, not a
  clean `exit()`. `SymSparseLinSOE` does NOT override `getA()`, so the fault was
  never in an accessor -- it is in `zeroA()`, reached via `formTangent()`, where
  `penv[size]` dereferences a null `penv` and `blkPtr->beg` a null `first`. The
  `memset(diag, 0, size*sizeof(double))` one line earlier is harmless at size 0,
  which is why the fault appears on the following statement. Lesson: `printA`
  calls `formTangent()` BEFORE `getA()`, so an unsized SOE has two distinct ways
  to die and fixing the accessors alone is not enough. Measure every system.
- **And a THIRD crash in the same class, on the same null-`first` path: the
  DESTRUCTOR.** `~SymSparseLinSOE` walks the row segments with
  `while (1) { if (blkPtr->next == blkPtr) { if (blkPtr != NULL) ...` -- it
  dereferences `blkPtr->next` BEFORE the null check inside the loop. Destroying a
  never-sized SOE therefore access-violates at TEARDOWN. **Its failure shape is
  the reason it hid for so long:** it fires on interpreter exit, after the
  script's final output has already flushed, so the run looks completely normal
  and merely exits nonzero. My own ad-hoc probe grepped stdout for a success line
  and reported SURVIVED for it -- twice -- while the process was dying every time;
  only the pytest gate, which asserts on `returncode`, caught it. **If you write a
  crash probe, assert on the EXIT CODE, not on output.** A process can print
  everything you asked for and still be dead.
- **Status (2026-08-18): FIXED, and wider than first scoped.** Measured
  per-system on an unsized SOE, `printA` died on 2 of 10 and `printB` on 6 of 10.
  Fixed in **12** SOE classes: `getA` (the only two that override it,
  `FullGenLinSOE` + `DiagonalSOE`), `getX`/`getB` in all twelve, and in
  `SymSparseLinSOE` two further null dereferences -- `zeroA()` and the
  destructor. All three accessors in
  both classes now report themselves and return an empty result -- `getA()` returns
  0 (which `OPS_printA` and `LinearSOE::saveSparseA` already branch on), `getX()`/
  `getB()` return a shared size-0 `Vector` (which `OPS_printB` already branches on).
  Nothing in the analysis flow changes: those callers only run post-`domainChanged()`.
  Gate: `tests/test_printa_unsized_soe.py` (zone_a, subprocess-isolated -- on a
  regression the CHILD dies and the parent reports a normal failure instead of the
  whole pytest run being killed). It parametrizes over EVERY serial `system` and
  over both `printA` and `printB`, because the first cut of this gate covered only
  three systems and would have passed while six others still died.
- **Follow-up pass (2026-08-18): the MPI-only half is now fixed too, and the
  remaining three items are closed.** `getX`/`getB` fixed in `MumpsSOE` (which
  also covers `MumpsParallelSOE` -- it derives and overrides neither),
  `MPIDiagonalSOE`, `DistributedDiagonalSOE`, `DistributedSparseGenRowLinSOE` and
  `PetscSOE`; plus `MPIDiagonalSOE::getpartofA` (stale `"FATAL ...::getA"` text
  inside `getpartofA`, and its own `exit(-1)`) and `DiagonalSOE`'s sized
  constructor.
- **The `grep exit(` survey MISSED AN ENTIRE SUBFAMILY, and the MP gate caught it
  on its first run.** Five *parallel* classes override `getB()` with a COLLECTIVE
  merge that had **no null check of any kind** -- straight to `*myVectB` /
  `*vectB`: `MumpsParallelSOE`, `DistributedBandGenLinSOE`,
  `DistributedBandSPDLinSOE`, `DistributedProfileSPDLinSOE`,
  `DistributedSparseGenColLinSOE`. They are invisible to a survey built on
  `grep exit(` **because they have no `FATAL` text to find -- they just crash.**
  Measured on `MumpsParallelSOE` (the only one with a Python door, via MP
  `system Mumps`): `mpiexec -n 2` + `printB` on an unsized SOE => rank 1
  `0xC0000005`, rank 0 down with it at `EXIT STATUS: -1`. A two-step probe
  bisected it cleanly -- `printA` alone survives both ranks (`len=0`), `printB`
  alone kills both. **Method that actually finds this family: enumerate every
  `getX`/`getB`/`getA` DEFINITION in the tree and ask "does it null-guard?",
  rather than grepping for the symptom you already know about.** The same sweep
  cleared four look-alike false positives -- `UmfpackGenLinSOE`, `PFEMLinSOE`,
  `SparsePythonCOOLinSOE`, `SparsePythonCompressedLinSOE` hold their `Vector` BY
  VALUE (`return X;`), so they were always safe. **22 SOE classes total.**
- **A guard inside a COLLECTIVE method must fire on every rank or not at all.**
  These `getB()` overrides have workers send-then-receive while the master gathers
  from each and replies. An early return taken by only some ranks leaves the rest
  blocked forever -- a hang is strictly worse than the crash it replaced. The
  guard is safe here for a specific, checkable reason: the wrappers are allocated
  in `setSize()`, reached only through the collective `domainChanged()`, so in the
  unsized state NO rank has them, every rank returns, and the collective is never
  entered. Whenever you add a bail-out to a collective, write down why the
  predicate is rank-uniform.
- **A `system` your box has and CI does not turns a crash-probe into a false
  alarm.** `system Pardiso` exists only under `#ifdef _PARDISO`, i.e. only where
  MKL is present. On the Linux CI runner `ops.system("Pardiso")` prints
  `WARNING unknown system type Pardiso` and raises `OpenSeesError`, so the child
  exits 1 -- **indistinguishable, to a probe that asserts on the exit code, from
  the very crash the gate exists to detect.** First CI run of this gate: 4 Pardiso
  rows reported "child process DIED", 1965 other cases green, nothing actually
  wrong. Deleting `Pardiso` from the list would be the wrong fix (it is one of the
  six systems measured dying pre-fix, on any box that HAS MKL). The gate now
  probes each `system` once per build and SKIPS the unsupported ones -- and the
  probe concludes "unsupported" **only** when OpenSees itself says
  `unknown system type`; a child that dies any other way is treated as a crash and
  is NOT skipped, so the escape hatch cannot swallow a regression.
- **`stdin=DEVNULL` is not optional in a subprocess-spawning test, and three files
  still lack it.** Under pytest's capture, fd 0 is left in a state `Popen` cannot
  `DuplicateHandle`, so any child that INHERITS stdin dies with
  `OSError: [WinError 6] The handle is invalid` (its sibling is `[WinError 50]`).
  It depends on the handle type the invoking shell supplies, which makes it look
  like a flaky code regression: measured on
  `tests/test_adr74_numberer_1.py` + `tests/test_ladruno_up_mp_smoke.py` with **no
  other file collected**, 20 passed under a piped stdin and 17 FAILED under
  `< NUL`, same binary, same commit. `tests/test_printa_unsized_soe.py` is immune
  in every configuration because it passes `stdin=subprocess.DEVNULL`. If an MPI
  or subprocess test starts failing at `subprocess.py` with an `OSError`, check
  the launcher before you touch the physics -- and add `stdin=DEVNULL`.
- **Deliberately NOT touched:** `ShadowPetscSOE::getX/getB` (same `petsc/` dir, so
  also never compiled) are unguarded too, but their shape is genuinely different
  -- an `MPI_Bcast`-driven shadow/actor protocol where an early return would
  desynchronize the actor, a design question that cannot be settled without a
  build. `badPetscSOE.cpp` is a stale duplicate of `PetscSOE` carrying the same
  defect and is left alone.
- **In MPI, `exit(-1)` in an accessor is a HANG, not a crash — and that is worse.**
  These solvers are collective in every phase. A rank that `exit()`s unilaterally
  leaves its peers blocked in their next MPI call, so the job does not fail with a
  diagnosable error; it stops making progress until the launcher's timeout fires.
  Identical failure shape to the rank-local parser `return -1` fixed in #742. When
  you are tempted to bail out of a rank-local code path, ask what the other ranks
  are waiting for.
- **The build map is NOT what the directory layout implies — check it before you
  promise a fix is compiled.** In `linearSOE/CMakeLists.txt`, both
  `add_subdirectory(petsc)` and `add_subdirectory(mumps)` are **commented out**.
  MUMPS still ships because the top-level `CMakeLists.txt` adds
  `mumps/Mumps*.cpp` straight to the `OpenSeesSP`/`OpenSeesMP` targets, bypassing
  that subdirectory entirely. PETSc has no such bypass, so **`PetscSOE.cpp`
  compiles in no target at all** and neither does `badPetscSOE.cpp` -- correcting
  the assumption that the latter is dead because it is "not in any CMake target":
  it IS listed in `petsc/CMakeLists.txt` (next to a stray `main.cpp`); that file
  is simply never processed. Conversely `DistributedDiagonalSOE.cpp` is compiled
  into **every** target, serial included, because `diagonal/CMakeLists.txt` is
  unconditional -- being named "Distributed" says nothing about where it builds.
- **Compiled, reachable, and exercised are three different things.** Of those five
  MPI classes: `MumpsSOE` and `MPIDiagonalSOE` are genuinely reachable from
  openseespy-MP (`system Mumps`, `system MPIDiagonal`) and are gated by the new MP
  rows. `DistributedDiagonalSOE` is instantiated only under `_PARALLEL_PROCESSING`
  via **Tcl** `system Diagonal` (openseespy's `system Diagonal` has no `#ifdef` and
  always yields plain `DiagonalSOE`), so it is compile-verified but has no Python
  door. `DistributedSparseGenRowLinSOE` is **dead**: it owns `classTag 21` and
  compiles into MP/SP, yet nothing in the tree ever constructs it -- not the Tcl
  `soe_table`, not `OpenSeesCommands.cpp`, not even `FEM_ObjectBrokerAllClasses`,
  which does have a case for its `DistributedDiagonalSOE` sibling. `PetscSOE` is
  not compiled. So 2 of 5 are testable and the ledger says so per row rather than
  implying uniform coverage.
- **`system MPIDiagonal` in a SERIAL build silently gives you a plain
  `DiagonalSOE`** (the `#else` arm of the `_PARALLEL_INTERPRETERS` guard, in both
  the Tcl and Python front-ends). A serial "pass" on that system therefore proves
  nothing about `MPIDiagonalSOE` -- it is a FALSE pass. The MP rows of the gate
  run under `mpiexec` for exactly this reason.
- **Do not let a "no output -> skip" guard swallow the bug you are hunting.** The
  natural way to make an MPI test portable is to skip when no rank output appears
  ("MPI infra?"). But a rank that `exit()`s produces no output either, so that
  guard converts the regression into a green skip. The gate establishes MPI health
  ONCE with a control driver that touches no SOE; after that control passes,
  missing rank output is a hard FAILURE. Same lesson as the destructor bug above,
  one level up: pick a signal that cannot be produced by the failure you are
  looking for.
- **Same family as the entry below**, reached by the other door: that one is
  `setSize()` skipping its wrapper-creation block on a model with ZERO free
  equations (`size == oldSize == 0`); this one is `setSize()` never running at all.
  Both ended in the same `exit(-1)`. If you are adding a new `LinearSOE`
  subclass, do not write `exit()` in an accessor -- return an empty result.
