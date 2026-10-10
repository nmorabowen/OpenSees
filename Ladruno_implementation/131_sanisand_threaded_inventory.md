---
title: WP-131 step 1 — SANISAND shared-state inventory for the threaded update loop (TIMs F19)
project: Ladruno
status: step 1 delivered (inventory, no code); step 2 (code) waits for WP-130
priority: medium
continues: 107_ladruno_openmp_element_loop (section 5.1, "root cause NOT located")
intake: _tims_2d_model_requests_2026-09-25.md, F19
plan: 127_tims_2d_requests_plan.md, WP-131
tags:
  - performance
  - openmp
  - threads
  - state-determination
  - re-entrancy
  - determinism
  - sanisand
---

# WP-131 step 1 — what SANISAND shares, and why IntScheme 1 segfaults

Checkout: `origin/ladruno` at `fb1afe58b`. Every `file:line` below was read on that
checkout. No source file is changed by this step and nothing was built.

Scope, as F19 asks: all shared mutable state reachable from `setTrialStrain`,
`commitState`, `revertToLastCommit`, `getTangent`, `getStress` and the responses of
`ManzariDafalias` (`SRC/material/nD/UWmaterials/ManzariDafalias.{h,cpp}`) and
`LadrunoSANISAND` (`SRC/material/nD/LadrunoSANISAND*.{h,cpp}`), including everything
they call: the tensor helpers, `Vector`/`Matrix`, `NDMaterial`, `opserr`, the
`LadrunoQuad` host and the WP-107 loop itself (`SRC/domain/domain/Domain.cpp:2693-2706`).

## 0. Headline

**The prime suspect for the IntScheme 1 segfault is `opserr`, not a data race.**
Under openseespy, `opserr` is a `PythonStream` whose every `operator<<` calls
`PySys_FormatStderr` (`SRC/interpreter/PythonStream.h:74`). Nothing in the fork
releases the GIL around `analyze` (no `Py_BEGIN_ALLOW_THREADS` / `PyGILState_*`
outside the Python-backed solvers). So the master thread holds the GIL through the
parallel region, and an OpenMP worker that prints calls into CPython with **no
thread state**. In CPython 3.12, `PySys_FormatStderr` reads the calling thread's
state first, and a worker has none. We expect that to be an access violation, and
it happens before `StandardStream` has written a byte.

We did not run a test for this. Every fact it rests on comes from the source or from
WP-107's own logs:

| WP-107 observation (`107_…md` §5.1, `Ladruno_files/testbed/perf/wp107/`) | explained by S1? |
|---|---|
| every bench run is openseespy (`wp107_strip_bench.py:43`, `import opensees`) | the precondition |
| `rc = 3221225477` = `0xC0000005` on every threaded TanType-2 run (`run1/sweep_sanisand_tan2.txt`) | yes, an access violation |
| the **serial** twin `run1/log_t2_sanisand_t1_r0.txt` prints `WARNING ManzariDafalias::ModifiedEuler() … substep cap 1000 reached` at step 4 (the bench passes `-maxSubsteps 1000`, `wp107_strip_bench.py:93`). Each **threaded** log stops dead after the `THREADED` line | yes: the first warning is fired from a worker and kills the process before any text is written |
| TanType 0 at `ds = 0.002`: the serial twin prints **no** warning, and threaded runs are 12/12 clean and bit-identical | yes: nothing prints, so nothing crashes |
| "crashes with **zero** warnings emitted" under `-Pmin 1e-8` | yes. `-Pmin` silences only the low-p clamp (`ManzariDafalias.cpp:1540`). The substep-cap warning (`:1662`) still fires |
| vanilla `--mat manzari` crashes identically | yes. The vanilla arm passes no `-Pmin` (`wp107_strip_bench.py:101`), so the default floor's clamp warning fires on the plastic branch |
| a mutex around `integrate()`, **and** `omp critical` around the whole `update()`, still crash | **yes, and only a thread-affinity cause explains this.** Serialising does not give a worker a Python thread state |
| `KMP_STACKSIZE` / `OMP_STACKSIZE` have no effect ("3/4 still crash") | yes. Stack size is irrelevant, and a run survives only when the first warning happens to be fired by the master |
| elastic material, same element and counts: 6/6 clean | yes: `ElasticIsotropicPlaneStrain2D` never prints |

**Discriminating test, zero code** (for step 2, or sooner). Re-run the crashing
configuration (`--scheme 1 --tan 2 --threads 4`) with `ops.logFile(path, '-noEcho')`
placed before `analyze`. `OPS_logFile` (`SRC/interpreter/OpenSeesOutputCommands.cpp:4767-4772`)
sets `echoApplication = false`, so `PythonStream::operator<<` skips `err_out` and
never touches CPython (`PythonStream.cpp:32-35`). The two outcomes:
- **S1 confirmed:** the crash disappears and the warnings land in the file.
- **S1 refuted, back to the ranking in §4:** it still segfaults.

A second test, also without code: run the same deck under the Tcl `OpenSees.exe`,
where `opserr` is a plain `StandardStream`. It should not crash (it may interleave
lines).

`-noEcho` is a **diagnostic only, not a fix**. `StandardStream` still writes to a
shared `ofstream` (`SRC/handler/StandardStream.cpp:266-268`) from several threads.
The fix is S1's treatment in §1.

## 1. Thread-affinity and output hazards (survive `omp critical`)

| # | State | Declared | Reached from the region by | Race? | Treatment |
|---|---|---|---|---|---|
| **S1** | `opserr` → `PythonStream::err_out` → `PySys_FormatStderr` | `SRC/interpreter/PythonModule.cpp:64-65` (`static PythonStream sserr; opserrPtr = &sserr`); `PythonStream.h:64-75`; all `operator<<` at `PythonStream.cpp:17-95` | any in-region print (list E below) | **Not a race: a thread-affinity fault.** Fatal on the first print from a worker, even when serialised | **Deferred per-thread message buffer.** Every in-region print goes to a `thread_local` buffer tagged with the element's loop index. After the parallel `for`, the master flushes the buffers to `opserr` in index order. The warning ORDER then matches serial at every thread count (F22 needs this). A lock does NOT fix it. `PyGILState_Ensure` from a worker would **deadlock**, because the master holds the GIL while it waits at the loop's implicit barrier. |
| S2 | `PythonStream::msg` (a `std::string` member, `PythonStream.h:60`), written by every `err_out` (`:68`) | shared by every thread | same as S1 | **Yes**: concurrent `std::string` assignment corrupts the heap. Hidden today because S1 kills the process first | Subsumed by S1: nothing prints from the region |
| S3 | `StandardStream::theFile` / `cerr` (`SRC/handler/StandardStream.cpp:257-272` and siblings) | `StandardStream` members | Tcl `OpenSees.exe`, and openseespy after `logFile` | `cerr` interleaves but is not UB. The `ofstream` **is** a data race | Subsumed by S1 |

**Every in-region `opserr` site on the SANISAND path.** These are the ones S1's
buffer must capture (E1–E9). Sites guarded by `debugFlag` are compiled out:
`ManzariDafalias.cpp:63` defines it `static const bool … = false`.

| E# | Site | Guard / budget | Reached when |
|---|---|---|---|
| E1 | `ManzariDafalias.cpp:1540-1549`, ModifiedEuler low-p clamp | `static std::atomic<int>` `:1538`, 10 per process | IntScheme 1 and the default scheme. Also through MaxStrainInc / MaxEnergyInc (`:1261`, `:1337`) when they select ModifiedEuler. Fires when p < m_Pmin + m_Presidual |
| E2 | `ManzariDafalias.cpp:1662-1686`, substep cap | `static std::atomic<int>` `:1660`, 10 per process | Only with a cap set (`mMaxSubstepsInME > 0`, i.e. LadrunoSANISAND `-maxSubsteps`). **This is the site on WP-107's crashing run** |
| E3 | `ManzariDafalias.cpp:2124-2128`, RK45 banner | `static bool do_once` `:2122`, plain read-modify-write | IntScheme 4, the **first RK45 call in the process**. It fires every time, so under threads it crashes whenever a worker reaches RK45 first |
| E4 | `ManzariDafalias.cpp:2180-2188`, RK45 low-p clamp | plain `static int` `:2178` | IntScheme 4 |
| E5 | `LadrunoSANISAND.cpp:2556-2570`, post-latch refusal | plain `static int` `:2554` | After the WP-99 commit latch has fired (an `-implex` companion cap) |
| E6 | `LadrunoSANISAND.cpp:2976-2990`, round-off α_in | plain `static int` `:2974`, once per Gauss point via `mRoundoffAlphaInWarned` | `-flipAlphaIn vanilla`, on plastic increments |
| E7 | `LadrunoSANISAND.cpp:3208-3220`, D2 sign change | plain `static int` `:3207` | `-implex` |
| E8 | `LadrunoSANISAND.cpp:3600` onward, `-implexControl` rungs | `static double` `:3597` plus `static int`s `:3598` and `:3605` | `-implexControl` |
| E9 | `LadrunoSANISANDSasME.cpp:934-938`, SAS-ME refusal (in `ManzariDafalias::ladrunoSasIntegrate`, `:733`). Re-verified on `ladruno` after WP-129 merged (PR #871); the pre-merge draft had it at `:770-776` | **per-instance** warn-once `mLadrunoSas.warned` (`:932-933`); the pre-merge process-wide `static std::atomic<int>` was removed in review (#871 item 2). No shared budget left | IntScheme 129 (SAS-ME), the first refused update at each Gauss point (start-f/α, dT_min, h ≤ 0, low p, drift, α, cap). Harmless today, because `ladrunoThreadSafeUpdate()` refuses the whole family. Once SAS-ME is allowlisted, this is an S1 crash site like E1/E2: route it through the deferred buffer. The per-instance latch fixed the counting, **not** the print from a worker |
| — | `ManzariDafalias.cpp:2984` "Still outside with f" | its `//if (debugFlag)` is commented out, **but the line cannot be reached**: `if (jj == maxIter)` sits inside `for (jj = 0; jj < maxIter; …)` (`:2964`) | never. It is not a suspect |
| — | `ManzariDafalias.cpp:5663` (`Inv`, singular tensor) | none | never: `Inv()` has no call site |
| — | size-check errors in the tensor helpers, `:5401-5627` | none | never: every caller passes 6-vectors and 6×6 matrices |

`LadrunoSANISAND.cpp:3817` (companion) and `:2556`'s commit-side twin run in
`commitState`, which is serial. They are outside the region.

## 2. Process-wide mutable data (true data races, needing concurrency)

| # | State | Declared | Write sites | In the region? which IntScheme/flags | Race under the loop | Treatment | FP-determinism note |
|---|---|---|---|---|---|---|---|
| D1 | `ManzariDafalias::mElastFlag` (**class static**) | `ManzariDafalias.h:299`; defined `ManzariDafalias.cpp:65` | constructors `:232`, `:322`, `:391`; `recvSelf` `:790`; `updateParameter` `:880`, `:886`. All serial | **read only**: `:1044`, `:5017`, `:5059`, `:5101`; `LadrunoSANISAND.cpp:2591`, `:2946`. Every scheme | No: reads only | Leave it for F19. Making it per-instance is PR #841's job (the quirk "`ManzariDafalias::mElastFlag` is STATIC"). **Rule for step 2:** nothing in the region may ever write it | — |
| D2 | the tensor constants `mI1, mIIco, mIIcon, mIImix, mIIvol, mIIdevCon, mIIdevMix, mIIdevCo` (class static) | `ManzariDafalias.h:301-308`; `ManzariDafalias.cpp:67-74`; filled by `initTensorOps` `:75` (`.h:310`) | static initialisation only (grep: no assignment anywhere else) | read everywhere | No | Leave | — |
| D3 | `static const` scalars `one3 … mMaxSubStep` | `ManzariDafalias.h:350-356` | none | read | No | Leave | — |
| **D4** | **`LadrunoImplexGlobals`**, the process-wide ledger: two FP accumulators, `firstCommitter`, eleven `long` counters | `LadrunoSANISAND.cpp:1977-2212`; singleton `instance()` `:1980-1983` (a function-local magic static, so initialisation is thread-safe) | **In the region, and NOT only under `-implex`:** `noteReversalNoiseGuard` at `:2872` and `:2888` (reached from `ladrunoGuardReversalNoise`, which the plain path calls at `:2610`). `noteRefusalLatched` `:2576`. Under `-implex`: `:2709`, `:3206`, `:3282`, `:3385`, `:3554`, `:3585`, `:3699`, `:3730`, plus `:4233`, which runs from the in-region lazy stage flip (`:2599`) under `-implexFlipAbsorb`. Commit phase, serial: `:3776`, `:3814`, `:3945`, `:4021` (`noteCommitRound`), `:4023` (`accumulate`), `:4092`. Reset: `:2219` on `wipe` | **Yes: `long++` loses updates.** WP-107 missed that `nReversalNoise` sits on the **non-implex** path. Its guard `mImplexOpt.enabled` (`LadrunoSANISAND.cpp:1963`) does not cover it, and would not have mattered only because the base class refuses anyway | Integer counters: `thread_local` deltas with an integer reduction at the end of the phase, or `std::atomic<long>`. Both give an exact, order-independent count. FP accumulators: keep them in the serial commit phase. If commit is ever threaded, store the error per instance and replay the sum in element order after the loop | **`sumError` (`:1989`) is the one FP reduction.** Threaded, its last bits would depend on the thread count. `maxError` is a max, which is exact and order-independent. `firstCommitter` needs serial commit order |
| D5 | the warning budgets | atomic: `ManzariDafalias.cpp:1538`, `:1660`. Plain: `:2178`, `LadrunoSANISAND.cpp:2554`, `:2974`, `:3207`, `:3598`, `:3605` (+ `static double :3597`). (SAS-ME's budget is per-instance since #871: see E9, so no D5 entry) | at their `if`s, E1–E9 | as E1–E9 | Atomics: the count is right, but **which** 10 tags get printed depends on the schedule. Plain ints also lose increments | Fold them into S1's deferred buffer, with the budget counted **by the master at flush time**. That makes the printed set deterministic too | — |
| **D6** | **`Matrix::matrixWork` / `intWork`** (class static scratch) | `SRC/matrix/Matrix.cpp:48-52` | `Invert` `:587-617` (free + realloc if too small), `:623-624`; `Solve` `:373` / `:461` (same pattern); **`addMatrixTripleProduct` `:944` / `:1042`**, which writes it at `:969-1020`. That last one is not in 107's list and is common in element `getTangentStiff` | **Live SANISAND use: only `NewtonSol` `:3511`, `:3518`, `:3544`**, via `NewtonIter2` `:3157`, from `BackwardEuler_CPPM` `:2551` / `:2605`. **IntScheme 2 only.** Dead (no call chain): `NewtonIter` `:3061` / `:3072`, `NewtonSol2` `:3933-3965` (only from `NewtonIter3`, which has no caller), `NewtonSol_negP` `:4345-4377` (only from `NewtonIter2_negP`, which has no caller) | **Yes.** A 6×6 fits the fixed area (36 ≤ 400 doubles, 6 ≤ 20 ints), so SANISAND alone never reallocates. But concurrent `DGETRF`/`DGETRI` share one work array and one pivot array, which silently gives a **wrong inverse**. Any thread inverting something larger reallocates under the others: a **use-after-free** | IntScheme 2: replace the three `Invert`s with a local 6×6 LU whose work and pivot arrays are on the stack. This is contained in `ManzariDafalias.cpp`, which WP-130 already edits. A `thread_local` work area inside `Matrix` would be a library-wide vanilla change: out of scope | — |
| D7 | `Matrix` constructors lazily allocate `matrixWork` | `Matrix.cpp:64-66`, `:80-82`, `:125-127`, `:154-156` | first `Matrix` built in the process | not in practice: thousands exist before any analysis | Only if the first-ever `Matrix` were built inside the region | Leave | — |
| D8 | RK45 function-scope statics: `n, d, b, R, dDevStrain, r, nStress, nAlpha, nFabric, ndPStrain, dSigma1-6, dSigma, …, aCep1-6, aCep_thisStep, aD, thisSigma, thisAlpha, thisFabric`, plus `do_once` | `ManzariDafalias.cpp:2122`, `:2137-2144` | throughout `RungeKutta45` (`:2117-2427`) | IntScheme 4 | Yes: hard, every call | Stack locals (for these fixed 6-vectors a stack-backed `Vector(double*, 6)` avoids the heap, as `LadrunoQuad::update` does), or per-instance scratch members. `do_once` goes into S1's buffer, or becomes `std::call_once` | — |
| D9 | `NewtonIter` function-scope statics `sol, R, R2, dX, norms, aux, jaco, jInv` | `ManzariDafalias.cpp:3019-3025` | inside `NewtonIter` | **None: `NewtonIter` has no call site in `SRC/`** (grep; the `SAniSandMS` hits are a different class) | No, because it is dead code | Delete, or leave. **It is not IntScheme 2's hazard.** That is D6 | — |
| D10 | **wrapper class-static return buffers**: "helpers returning const-refs to static buffers" at class scope, the pattern 107's function-scope grep could not see | `LadrunoSANISANDPlaneStrain.h:106-110` / `.cpp:53-57`; `ManzariDafaliasPlaneStrain.h:78-82` / `.cpp:27-31`; `LadrunoSANISAND3D.h:104-105` / `.cpp:50-51`; `ManzariDafalias3D.h:78-79` / `.cpp:27-28` | `LadrunoSANISANDPlaneStrain.cpp`: `getStrain` `:154-156`, `getStress` `:165-167`, `getStressToRecord` `:176-179`, `getTangent` `:196-204`, `getInitialTangent` `:213-221`. `ManzariDafaliasPlaneStrain.cpp`: `:102-104`, `:113-115`, `:124-127`, `:144-152`, `:161-169`. `LadrunoSANISAND3D.cpp`: `:133`, `:141`, `:148`, `:158`. `ManzariDafalias3D.cpp`: `:97`, `:105`, `:112`, `:122` | **Not in loop A with `LadrunoQuad`:** `LadrunoQuad::update` (`LadrunoQuad.cpp:736-791`) calls only `setTrialStrain`, and nothing inside `setTrialStrain` calls a getter (grep). Reached from formUnbalance / formTangent (serial) and the recorders | Not today. **A hard race** the moment any threaded loop (loops B/C), or any host whose `update()` reads `getStress`/`getTangent`, runs it: every instance returns a reference into the same 3/4/9-double buffer | Per-instance `mutable` members (3 + 3 + 4 + 9 + 9 doubles per Gauss point in plane strain) | — |
| D11 | `ManzariDafalias::getPStrain` `static Vector result(6)` | `ManzariDafalias.cpp:5748` | `:5749` | Plane-strain wrappers do not override `getPStrain`, so response 8 (`:664`) reaches the base. It also prints an "error" line on **every** call (`:5747`). Recorder phase, serial | Not today | Return a per-instance member | — |
| D12 | LadrunoSANISAND response statics: `probe*` in `setResponse` (`LadrunoSANISAND.cpp:4499-4600`); `out*` in `getResponse` (`:4610-4697`) | there | there | recorder phase, serial | Not today | Leave. Make them per-instance only if recorders are ever threaded | — |
| D13 | commit-refusal counter | `SRC/material/LadrunoMaterialStatus.h:116-120` | `:123`, from `LadrunoSANISAND::commitState` `:4057` and `ladrunoImplexCommit` `:3872` (serial) | commit phase | Not today (its own note says: atomic if commit is threaded) | Leave | — |
| D14 | `sendSelf`/`recvSelf` statics: `ManzariDafalias.cpp:679`, `:755`; `LadrunoSANISAND.cpp:1726`, `:1805` | there | there | serialisation only | No | Leave | — |
| D15 | object counters and constructor latches: `numManzariDafaliasMaterials` `ManzariDafalias.cpp:77`; `numLadrunoSANISANDMaterials` `LadrunoSANISAND.cpp:150`; `warnedIntSchemeNoErrorControl` `ManzariDafalias.cpp:212`, `:305` | there | parser / constructor | model build, serial | No | Leave | — |
| D16 | `NDMaterial` sensitivity dummies | `NDMaterial.cpp:378`, `:385`, `:398`, `:405`, `:412` | there | not called by SANISAND | No | Leave | — |
| D17 | `ops_Dt` (a global) | e.g. `SRC/database/main.cpp:51` (one definition per target) | integrators, serial | **read** in the region: `LadrunoSANISAND.cpp:2870` (reversal guard); in the commit phase at `:3775`, `:4079` | No | Leave | — |
| D18 | `ops_TheActiveElement`; `LadrunoQuad::shp`/`shpBar` | `Element.cpp:59`; `LadrunoQuad.cpp:72-73` | the loop body; `shapeFunction` | yes | No: already `thread_local` (WP-107) | Leave | — |

**Per-instance state written in the region.** This is fine because `getCopy` gives
every Gauss point its own instance. It is listed so step 2 does not re-audit it:
`mAlpha_in` (`ManzariDafalias.cpp:1039`, `:1041`); every trial member written by
`integrate` (`:1044-1060`); `mSubstepsTakenInME` / `mSubstepCapHitInME` (`:1023-1024`,
`:1643`, `:1645`); `mUseElasticTan` (`:1744`); `mFlipSeen` / `mPrimed` / α
(`LadrunoSANISAND.cpp:2591-2600`, `ladrunoRunStageFlipOnce` `:4164-4247`);
`mRoundoffAlphaInWarned` (`:2973`); `mImplex*`. The plane-strain `getCopy` is
`*clone = *this` (`LadrunoSANISANDPlaneStrain.cpp:95-100`), a deep copy through
`Vector::operator=`.

**Allocation (a speed-up concern, not correctness).** `ModifiedEuler` builds dozens
of heap `Vector`s per substep (the helpers return by value, e.g. `GetDevPart`,
`DoubleDot4_2`, `GetNormalToYield` at `:5191`, `GetStiffness` at `:5110`). The MSVC
CRT heap is thread-safe, so this is not a crash candidate. But at roughly 24 substeps
per point (TIMs §1.2), heap contention may cap step 2's speed-up table. If it does,
stack-backed scratch as in `LadrunoQuad::update` is the remedy.

## 3. Determinism (thread-count identity)

- Each Gauss point's arithmetic is sequential and per-instance. Once D6 and D8 are
  out of the path, a point's result does not depend on the thread count. The
  loop's `reduction(+:sum)` is an integer sum (`Domain.cpp:2693`).
- The only floating-point reduction in scope is D4 `sumError` (`accumulate`,
  `LadrunoSANISAND.cpp:4023`). It must stay serial, or become an ordered replay.
- Output is not bit-level, but F22 cares about it: with atomics (D5), **which**
  warnings print and in **what order** depend on the schedule. S1's ordered flush
  fixes both.
- The D4 integer counters are exact under atomics or a `thread_local` reduction.

## 4. Suspects ranked for the IntScheme 1 segfault

1. **S1: `opserr` → `PySys_FormatStderr` from a worker with no Python thread
   state.** `PythonStream.h:74`; emitted at `ManzariDafalias.cpp:1662` (the site on
   the crashing run) and `:1540`. It matches every row of 107's hypothesis table,
   including the "`omp critical` does not help" row that eliminates every data-race
   explanation. §0 gives the zero-code test.
2. **S2: the `PythonStream::msg` `std::string` race** (`PythonStream.h:60`). This is
   heap corruption from two concurrent printers. It is latent behind S1, would show
   up the moment S1 were papered over with `PyGILState`, and is fixed by the same
   buffer.
3. **D4: `LadrunoImplexGlobals::nReversalNoise++`** on the non-implex path
   (`LadrunoSANISAND.cpp:2872`, `:2888`). A real in-region race that 107 did not
   list. It is benign: lost counts, no crash. It cannot explain the **vanilla**
   crash.
4. **D10: the plane-strain wrapper class-static buffers**
   (`LadrunoSANISANDPlaneStrain.cpp:53-57`, `ManzariDafaliasPlaneStrain.cpp:27-31`).
   This is the pattern the function-scope audit could not see. It is not on loop A's
   path with `LadrunoQuad`, so it cannot be the crash, but it is the first thing that
   breaks when loops B/C are threaded.
5. **D6: `Matrix::matrixWork` (and `addMatrixTripleProduct` using it).** This is not
   on the IntScheme 1 path. It is the real hazard for IntScheme 2, and for any host
   element whose `update()` builds a tangent.

Ruled out by source or by 107's measurements: the function-scope statics on the
IntScheme 1 graph (there are none); stack depth (measured); the element (the
elastic arm); the solver (measured); `debugFlag` prints (compiled out); `:2984`
(unreachable); the heap allocator (thread-safe CRT).

## 5. What WP-130 (F18(c), "de-static `NewtonIter`") already removes

- **Nothing that is live.** `NewtonIter` (`ManzariDafalias.cpp:3005`, statics
  `:3019-3025`) has **no call site**. IntScheme 2 runs `BackwardEuler_CPPM` →
  `NewtonIter2` (`:3122`) → `NewtonSol` (`:3284`), and `NewtonIter2` has no statics.
  Recommendation to WP-130: aim F18(c)'s "remove the static work arrays" at
  **D6**, the three `Matrix::Invert` calls in `NewtonSol` (`:3511`, `:3518`,
  `:3544`). That is the live shared state on IntScheme 2. Delete `NewtonIter` rather
  than de-static it.
- If WP-130 does replace those `Invert`s with local solves, **D6 is closed for
  SANISAND**, and IntScheme 2 can join the allowlist in step 2 (behind the same
  measurement protocol).
- WP-130 touches none of S1, S2, D4, D5, D8 or D10.

## 6. Step 2 (after WP-130): the order of work

1. Run the §0 `-noEcho` test on the WP-107 crashing configuration (t2, 4 and 8
   threads). If S1 is confirmed, it re-states 107 §5.1's root cause.
2. S1 + D5: an in-region deferred-message facility, flushed in index order after the
   loop in `Domain::update`. The facility is element-agnostic, lives next to the
   WP-107 loop, and has one entry point that materials call instead of `opserr`
   while the region is active.
3. D4: integer counters become atomics (or `thread_local` + reduction). `sumError`
   stays serial.
4. D8 (IntScheme 4) and D6 (IntScheme 2, if WP-130 left it) become stack scratch.
   D10 and D11 become per-instance, even though they are not loop A, because they
   are the next trap.
5. Allowlist per scheme: 1, then 2 after WP-130, then 4. Keep refusing under
   `-implex` until D4 is done. Then run F19's test: a few thousand SANISAND points
   at 1/2/4/8 threads with `MKL_NUM_THREADS=1`, checking identical committed curves
   and the speed-up table. The run must include a deck that **prints**, i.e. hits a
   cap or a clamp. WP-107's lesson is that a quiet deck proves nothing.
6. Ledgers: an implementations row; 107's refusal table updated; a quirks row,
   "`opserr` under openseespy calls CPython: never print from an OpenMP worker",
   and, since it is greppable, a `ci/check_quirk_patterns.py` rule for `opserr`
   inside `#pragma omp` regions.
