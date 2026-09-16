---
title: WP-107 — threaded Domain::update element loop (OpenMP), desktop-scoped
project: Ladruno
status: implemented — default OFF, opt-in via `ladrunoThreads`
priority: medium
implements: 75b_ladruno_threaded_assembly_adr (stage L3-1, desktop-scoped re-entry)
amends: 40_ladruno_performance_adr
tags:
  - performance
  - openmp
  - threads
  - state-determination
  - re-entrancy
  - determinism
---

# WP-107 — OpenMP over `Domain::update()` (ADR-75b "loop A")

> **Read [[75b_ladruno_threaded_assembly_adr]] first.** This WP does not re-derive
> its taxonomy, its determinism policy, or its hazard inventory; it implements one
> stage of it and cites the rest.

## 1. Why this exists at all, given ADR-75b closed the lane

ADR-75b §13 **closed Lane 3 for the production/cluster path**. G-L3 measured, on a
3D solid cantilever under `system Mumps` at np=16:

| N (DOF) | loop A | loop B | loop C | solve | verdict |
|---|---|---|---|---|---|
| 28 611 | 1.83% | 8.23% | 1.06% | 87.77% | FAIL |
| 143 811 | 0.58% | 2.07% | 0.31% | 96.70% | FAIL |
| 540 675 | 0.26% | 0.95% | 0.13% | 98.54% | **FAIL by ~42x** |

That verdict stands and nothing here contradicts it. What §13 explicitly left open
is its own last paragraph:

> "**So Lane 3 is at best a desktop-only optimization on models ≲136k DOF** …
> §2's five-loop taxonomy, §3's determinism policy and §4's scatter remedy remain
> correct and are the reusable assets; they should be **cited, not re-derived**, if
> a desktop-scoped case is ever made."

The TIMs plane-strain strip is that case, and it is a different regime in three
ways the cluster deck is not:

1. **It never leaves the desktop.** `system Pardiso` / `BandGeneral`, 10k–80k
   Gauss points, well under the ≲136k-DOF ceiling §13 names.
2. **Its element kernel is ~2 orders of magnitude more expensive per Gauss point
   than G-L3's.** G-L3 ran `stdBrick` + `J2Plasticity` (and re-ran on
   `LadrunoBrick` + `LadrunoJ2`, finding "essentially the same wall"). The strip
   runs **LadrunoSANISAND**, a substepping SANISAND-family model whose
   `ModifiedEuler` integration can take hundreds of substeps at one Gauss point.
   §13's own escape clause — "lifting 0.95% to 40% needs a ~66x more expensive
   element kernel" — is a statement about `LadrunoBrick`/`LadrunoJ2`, not about a
   substepping critical-state model.
3. **The solve is small, not dominant.** G-L3's 98.54% solve was a 40.2 s MUMPS
   factorization at 540k DOF. At strip scale the factorization is milliseconds.

So the gate has to be **re-measured in this regime**, and that measurement — not
this document — is what authorizes using the feature. See §6.

## 2. What is implemented, and what is deliberately not

**Implemented: loop A only** — the `theEle->update()` loop in `Domain::update()`,
i.e. state determination = the material integration.

ADR-75b §2.1 is the reason this is the only loop worth taking on its own terms:
`Domain::update()` contains **no reduction into the SOE** — no `addA`, no `addB`,
no shared accumulator. For an element whose `update()` touches only its own state
and its own materials, the threaded loop performs the identical arithmetic in the
identical order. So:

> a correct threaded loop A is **bit-identical to serial at every thread count**,
> and the fork's existing byte-identical oracles gate it unchanged — no new CI
> mode, no ordered/fast split, no tolerance.

**Not implemented: loops B and C** (`formTangent` / `formUnbalance`). §4 says why.

## 3. The machinery

### 3.1 Build

`option(LADRUNO_OPENMP … OFF)`. `Ladruno_scripts\build.bat` passes
`-DLADRUNO_OPENMP=ON` (override with `set LADRUNO_NO_OPENMP=1`).

**Decision, recorded so it is not re-litigated:** the CMake option defaults OFF so
a bare `cmake` build of this tree stays vanilla-shaped, but build.bat turns it ON,
because a capability that must be recompiled to be tried is a capability nobody
tries — and with the runtime default at 1 thread the serial path is byte-identical
anyway (verified both at 1 thread and with the option compiled out).

**The compile flag is scoped to `OPS_Domain` + `OPS_Utilities`, PRIVATE.** It is
deliberately NOT global. The fork carries 7 pre-existing `#pragma omp` lines in
PFEM (ADR-75b §1) that have always been compiler-ignored no-ops. Putting `/openmp`
on `OPS_Element` would silently activate untested PFEM threading as a side effect
of this WP. Keeping the flag off that target leaves them dead.

**The OpenMP runtime is `libiomp5md`, not `vcomp` — checked, not assumed.** The
serial targets already link `mkl_intel_thread` + `libiomp5md` for desktop PARDISO
(ADR-75 P1b), and MSVC `/openmp` emits `_vcomp_*` calls plus a `vcomp.lib`
default-lib directive, so the expectation was two coexisting OpenMP runtimes. It
is not what happens: `dumpbin /DEPENDENTS distin\opensees.pyd` lists
`libiomp5md.dll`, `mkl_intel_thread.3.dll`, `mkl_core.3.dll` and **no
`vcomp140.dll`** — Intel's runtime exports the Microsoft OpenMP ABI, so the
linker resolved our `_vcomp_*` references against it and there is exactly **one**
runtime in the process.

That is the good outcome, and it is load-bearing enough to re-check after any
change to the MKL link: if `vcomp140.dll` ever appears in that list, the process
has two thread pools, and `KMP_*`/`OMP_*` environment settings will then apply to
only one of them. Note also that this is why `KMP_STACKSIZE` is a valid knob for
our worker threads (used as a diagnostic in §5.1).

The thread count is still applied through an explicit `num_threads(n)` clause and
never `omp_set_num_threads()` — that is about not disturbing MKL's own threading
and not waking PFEM's dormant pragmas, and it holds regardless of which runtime
is linked.

### 3.2 Runtime knob — there is exactly one

ADR-75b **P-5**: "MKL solver threads × OpenMP assembly threads × MPI ranks will
oversubscribe and make every bench lie. Assembly threading must not ship without a
single documented knob."

```
ladrunoThreads            ;# Tcl: query
ladrunoThreads 4          ;# Tcl: request 4
ops.ladrunoThreads()      # Python: query
ops.ladrunoThreads(4)     # Python: request 4
```

plus the `LADRUNO_THREADS` environment variable, which seeds the value once at
first use. **Default is 1** — ADR-40's standing anti-goal is "OpenMP-by-default",
and at 1 thread `Domain::update` takes the serial path, not a one-thread parallel
region. Asking for `n > 1` in a binary built without `LADRUNO_OPENMP` **warns**
rather than silently running serial, because a silently-serial "threaded" run is
how a bench lies.

Nothing in the fork ever calls `omp_set_num_threads()`. The count reaches exactly
one `#pragma omp parallel for … num_threads(n)`.

**Every rejected value is announced** (red-team S7). It used to fail silently on
the env path — the one a bench harness and every job script actually use — which
is precisely how a bench lies:

| `LADRUNO_THREADS` | before | now |
|---|---|---|
| `abc` | ran SERIAL, **no message at all** | `WARNING … "abc" is not an integer -- IGNORED … SERIAL (1 thread).` |
| `2.7`, `4x` | ran SERIAL, no message | `WARNING … has trailing text after the number -- IGNORED …` |
| `0`, `-3` | ran SERIAL, no message | `WARNING … requested 0 threads, which is < 1 -- using 1 (SERIAL).` |
| `99999` | honoured as **1024 threads**, announced as a threaded run, measured **5x slower than serial** | `WARNING … requested 99999 threads but this box reports 24 hardware threads (OpenMP max 24) -- CLAMPED to 24.` |

The ceiling is `omp_get_num_procs()` (falling back to
`std::thread::hardware_concurrency()`, and clamping nothing when neither can
answer) — deliberately NOT `omp_get_max_threads()`, which reports `OMP_NUM_THREADS`
if the environment set it: that is a policy, not the hardware, and silently
clamping a deliberate `ladrunoThreads 8` down to a stray `OMP_NUM_THREADS=2`
would be its own surprise. Both interpreter verbs clamp through the same helper,
so the env path and the verbs say the same thing in the same words.

### 3.3 The loop and the five hazards ADR-75b §5.4 said grep could not see

| ADR-75b hazard | What WP-107 does |
|---|---|
| **H2** element kernel statics | Per-class opt-in `Element::ladrunoThreadSafeUpdate()`, **default `false`**. See §5. |
| **H3** the loop cursor is shared and mutable (`Domain::getElements()` returns the single member `theEleIter`; a naive shared-cursor pull **skips and duplicates** elements) | The iterator is drained into a `std::vector<Element *>` **before** the parallel region; the region is an indexed `for`. |
| **H4** `Node`'s trial-state *getters* lazily heap-allocate (`getTrialDisp` → `createDisp`) | A serial pre-pass calls `getTrialDisp/getTrialVel/getTrialAccel` on every node. Those three getters cover all three `create*` functions, so nothing can allocate inside the region. |
| **H6** the deep profiler's per-element buckets are an unsynchronized shared `std::map` write | **Refuse**: if `theProfiler().deep()` is armed, the element loop runs serial with a one-time warning. Coarse `profiler start -perStep` is unaffected and is what a threaded run is timed with. |
| **H7** `ops_TheActiveElement` is a shared global read per call by ≥7 materials for the lch latch | Made **`thread_local`** (`SRC/element/Element.cpp` + both `extern` declarations). Serial behaviour unchanged. |

Plus, from ADR-75b §2.1b: elements that write **shared node trial state** inside
`update()` — `LadrunoRigidBody` and `ZeroLengthVG_HG` — are excluded automatically,
because neither is on the allowlist (nothing is, by default).

**The allowlist is all-or-nothing.** If *any* element in the domain answers
`false`, the whole loop runs serial. A mixed parallel/serial split was rejected: it
would still let an audited element run concurrently with an un-audited one that
writes shared node state.

**Error handling is deterministic.** `ok` is an **integer** sum, so
`reduction(+:sum)` reproduces the serial `ok += theEle->update()` exactly at every
thread count. The reported failing element is the one with the **lowest index in
serial iteration order**, resolved in a `critical` section entered only on a
failure — so the diagnostic does not depend on which thread got there first.
(MSVC implements OpenMP 2.0, where `reduction(min:)` does not exist; hence the
critical section rather than a min-reduction.)

**Parallel builds.** `PartitionedDomain` and `Subdomain` override
`ladrunoThreadedUpdateAllowed()` to `false`. `OpenSeesSP`/`MP` therefore behave
exactly as before; ADR-75b §11 q6 (hybrid MPI+threads) is untouched and still
deferred.

## 4. Why loops B and C are NOT threaded — the audit, not a shrug

The task this WP came from asked for a Phase 2: parallel tangent/residual
computation with serialised assembly, *conditional on the FE_Element path being
re-entrant*. It is not, and the blocker is structural rather than incidental:

- **`FE_Element::theTangent` is a class-wide pool.** `static Matrix **theMatrices`
  (`FE_Element.cpp:51`, `.h:132`) and `theTangent = theMatrices[numDOF]` hand **one
  Matrix to every same-`numDOF` FE_Element**, and
  `IncrementalIntegrator.cpp:118-120` passes that address straight to `addA`. That
  is a **100% collision** which no *element* allowlist can mitigate (ADR-75b
  §5.4-H1). Under `constraints Transformation` — which the strip uses — the same
  holds one level up (`TransformationFE`: `modMatrices`, `dataBuffer`,
  `localKbuffer`, `dofData`).
- **The element's own tangent buffer is class-static too.** `LadrunoQuad::K` /
  `::P` are `static Matrix`/`static Vector` shared by every instance
  (`LadrunoQuad.cpp:70-71`), and `getTangentStiff()` returns a reference to `K`.
- De-pooling both into per-instance buffers is ADR-40's rank-7 item and carries a
  **serial** memory cost (ADR-75b §11 q9 prices per-element de-statication at
  order 1–4 GB at 325k elements) paid on every run, threaded or not.

So Phase 2 would have meant shipping either (a) a threaded loop over
non-re-entrant code, or (b) a GB-scale serial memory regression, for a loop that is
*not* where this deck's time goes. Neither is acceptable, and the measurement in §6
is what decides whether it is ever worth revisiting.

**If it is revisited, the design is not "atomic scatter".** For a desktop deck the
cheap exact form is a **chunked ordered gather**: compute a chunk of C elements in
parallel into C per-slot buffers, then replay `addA`/`addB` over that chunk in
serial `(FE index, i, j)` order, then advance. Storage is `C × idSize²` (bounded
and small), not `Σ idSize²` (ADR-75b §4.2's 1.5 GB at 1M DOF), and it is
bit-identical because the replay order is the serial order. It still needs
per-instance FE_Element and element buffers first.

## 5. Which elements and materials are thread-safe — the guide paragraph

The rule is **"un-audited is never threaded"**: `Element::ladrunoThreadSafeUpdate()`
and `Material::ladrunoThreadSafeUpdate()` both return `false` in the base class, so
every element and every material in the fork refuses by default and the loop falls
back to serial with a named tag in the warning. Only the classes below opt in.

### Elements

| Class | Threaded? | Why |
|---|---|---|
| **`LadrunoQuad`** — formulations **std / bbar / ssp**, `-geom linear` | **YES**, if every material agrees | `update()`'s three function-scope `static` buffers (`u`, `eps`, `B`) are now **stack-backed** non-owning `Vector`/`Matrix`; the class-static `shp` / `shpBar` are **`thread_local`**; `pts`/`wts` are written only in the constructor; the path only *reads* node state. |
| `LadrunoQuad` — formulation **eas** | **NO** | `formEAStrue()` runs on ~12 shared function-scope statics **and** condenses through matrix inversion, which uses the process-wide `Matrix::matrixWork` scratch. |
| `LadrunoQuad` — `-geom finite` | **NO** | `updateFinite()`'s `static Matrix Fm`. |
| **every other element in the fork** (LadrunoBrick, BezierTri6, BezierTet10, stdBrick, FourNodeQuad, the beam-columns, contact, …) | **NO** | Not audited. ADR-75b §5.4-H2 measured **63 of 215** `Element::update()` bodies carrying non-const function-scope statics (473 declarations), and §5.4-H5 adds `SRC/coordTransformation/` (434 declarations / 18 files) on the path of every beam-column. Each is its own audit. |
| `LadrunoRigidBody`, `ZeroLengthVG_HG` | **NO, and hard** | They write node trial state **shared with other elements** (ADR-75b §2.1b). This is an *ordering* race no reduction policy can fix; only exclusion or redesign works. |

### Materials

| Class | Threaded? | Why |
|---|---|---|
| **`ManzariDafalias`** (and `…3D` / `…PlaneStrain`), **every** `IntScheme` | **NO — measured, see §5.1** | `IntScheme 1` was allowlisted on a complete static audit and **segfaults anyway**. `IntScheme 2` additionally has `NewtonIter()`'s shared static work arrays and calls `Matrix::Invert`/`Solve`; `IntScheme 4` has ~20 static work arrays plus a `static bool do_once`; `3/5` and the MaxStrain/MaxEnergy family are un-audited. Its two `ModifiedEuler` warn budgets are `std::atomic<int>` regardless — their own THREAD-SAFETY note asked for that "if Lane 3 lands", and it is correct independent of this refusal. |
| **`LadrunoSANISAND`** (and its wrappers) | **NO** | Two independent refusals. (1) The base refuses (above). (2) Under `-implex` the diagnostics are a **process-wide ledger** (`LadrunoImplexGlobals`: `maxError`, `sumError`, `count`, four refusal buckets); two accumulators are **floating point**, so this is not fixable with an atomic — a threaded sum would change the reported average's last bits with the thread count, which is precisely the determinism this WP exists to preserve. (2) is kept explicit because it survives any fix to (1). |
| **every other material** | **NO** | Not audited. |

### 5.1 The measurement that removed SANISAND from the allowlist

This is the result that cost WP-107 its intended payoff, and it is the most
transferable thing in the WP.

`ManzariDafalias` under `IntScheme 1` (ModifiedEuler) has **no function-scope
static anywhere on its update call graph** — every `static Vector`/`Matrix` in
that file is in `RungeKutta45`, `NewtonIter`, `getPStrain` or
`sendSelf`/`recvSelf`, all off that path — and no `Matrix::Solve`/`Invert` on it
either. It was allowlisted on exactly that basis, and the first sweep returned
**12/12 runs bit-identical at 1/2/4/8 threads, `maxdiff = 0.000e+00`**.

Then a harder load path was tried, and it **segfaults**: 4/4 at 4 threads on a
6400-element deck, as soon as the PLASTIC branch is exercised in volume. The
elastic branch is clean, which is exactly why the first sweep passed — it had
barely entered plasticity. *A clean threaded run on a mild load path proves
nothing.*

| hypothesis | test | result |
|---|---|---|
| the fork's subclass | vanilla `nDMaterial ManzariDafalias` | crashes identically, 4/4 |
| solver interaction | `BandGeneral` vs `Pardiso` | both crash |
| worker-thread stack overflow | `KMP_STACKSIZE=64M` (libiomp5 is the runtime) | no change |
| `opserr` from inside the region | `-Pmin 1e-8`, so the clamp never warns | crashes with **zero** warnings emitted |
| the element | same element + `ElasticIsotropicPlaneStrain2D`, 10 000 elements, 8 threads | **6/6 clean, bit-identical** |

A bounded second round excluded three more hypotheses, and one of them reframes
the problem — the full table is in [[75b_ladruno_threaded_assembly_adr]] §14.1:

| test | result |
|---|---|
| mutex around the **whole `ManzariDafalias::integrate()`** | 0/4 clean |
| `#pragma omp critical` around the **whole `theEle->update()`** | **0/4 clean** |
| `OMP_STACKSIZE` / `KMP_STACKSIZE` = 256M | 0/4 clean |

Serializing the entire element update means the threads exist and enter the
region but **never run an update concurrently** — and it still faults. **So this
is not a data race between element updates**, which is what every hypothesis
before it had assumed. What is left is something about running this particular
update path on an OpenMP worker thread at all; the elastic path on the same
worker thread, same element, same counts, is clean 6/6. The blocker on getting
further was tooling: no `cdb`/`WinDbg`/`procdump` on the box and **no PDBs** in a
Release build, so no faulting frame could be obtained.

Root cause **not located**, so the family is refused. A located-but-unfixed
hazard is strictly worse than an un-audited one, because the audit manufactures
confidence. The next tool is **ThreadSanitizer**, which ADR-75b §7's correctness
protocol already names as "the only tool that finds the misses `grep` cannot" —
this is the measurement that makes that protocol item load-bearing rather than
belt-and-braces. TSan needs clang/gcc, so it is Esmeralda/Linux work, not MSVC
desktop work.

**Consequence for the WP's own gate.** Loop A was measured at **51.2 % of step**
on the SANISAND deck — it passes ADR-75b's >40 % single-loop gate that G-L3
failed by ~42x at cluster scale, so the *desktop-scoped premise of §1 is
confirmed*. What is not available is the payoff, because the only material that
makes loop A that large is the one now refused. The mechanism is shipped,
default-off, with a correct allowlist; the deck that motivated it cannot use it
until the defect is found.

### What refuses at the loop level (not the class level)

- an **MPI binary** — `OpenSeesSP`, `OpenSeesMP`, `OpenSeesPyMP` — checked first
  and widest (see the box below);
- a `PartitionedDomain` or a `Subdomain`;
- the **deep** profiler gate being armed;
- a thread count of 1 (the default);
- a binary built without `LADRUNO_OPENMP`;
- a domain with fewer than 2 elements.

Every one of these prints and runs the serial loop. None of them is silent.

> **The MPI fence is separate from the `PartitionedDomain` one, and it has to be**
> (red-team B2). The first version of this WP claimed "`PartitionedDomain` /
> `Subdomain` refuse outright, so SP/MP are untouched". Only SP was. Under
> `_PARALLEL_INTERPRETERS` — the `OpenSeesMP.exe` / `openseesmp.pyd` build, and
> the fork's dominant parallel idiom — every rank holds a **plain `Domain`**
> (`SRC/tcl/commands.cpp`'s `#elif _PARALLEL_INTERPRETERS → Domain theDomain;`
> and `SRC/interpreter/OpenSeesCommands.cpp`'s `theDomain = new Domain;`), so the
> overrides never fire and the loop was live on every rank. The count is seeded
> from an environment variable and `mpiexec` propagates the environment, so one
> `set LADRUNO_THREADS=8` in a job script turned an np-8 run into 64-way
> oversubscription, with the one announcement line per process lost in
> rank-interleaved stdout. ADR-75b §11 q6 defers hybrid MPI+threads; this WP had
> shipped it by accident.
>
> The fence keys on the **build**, not on `MPI_Comm_size`: `np == 1` under
> OpenSeesMP is still the parallel code path and is not the desktop case this WP
> measured. It lives in its own per-target translation unit,
> `SRC/utility/LadrunoParallelBuild.cpp`, because `OPS_Domain` is a single OBJECT
> library compiled once with neither parallel define — an `#ifdef _PARALLEL_*`
> written in `Domain.cpp` compiles to **nothing in all five targets**. ADR-78 P1
> measured that trap (`LadrunoContactAbort.h`) and this file copies its answer.

### Announcements are per DOMAIN, not per process

Red-team S3. The audit itself always re-ran on every `Domain::update()`, so
`wipe` + rebuild, a runtime `element`, and `remove element` were handled
correctly from the start — but the **messages** were latched behind three
process-wide `static bool`s, so only the first model in a process ever said
anything. A second model that quietly went serial and one that stayed threaded
looked identical, which is exactly what the announcement exists to prevent; the
perf sweep driver was safe only because it forks one process per run, while
pytest, apeGmsh and any in-process parameter study sat squarely in the blind
spot. The gate is now keyed on `(element generation, outcome, tag, thread
count)`, where the generation is bumped by `addElement` / `removeElement` /
`clearAll`. A changed element set, a different refusing tag, or a changed thread
count speaks again; a steady state stays quiet, so a Newton iteration prints
nothing.

### Two things a future auditor must check, which a `static` grep will not show

1. **`Matrix::matrixWork` / `Matrix::intWork`** (`Matrix.cpp:51-52`) are
   process-wide scratch used by `Matrix::Solve` and `Matrix::Invert` — which
   **free and reallocate** them when the matrix is larger than the fixed work
   area. Any threaded path that calls either is a use-after-free, not merely a
   race. This is why `IntScheme 2/4` and `LadrunoQuad -eas` are refused.
   (`Matrix`'s *constructors* also lazily allocate that buffer, but by the time
   any analysis runs thousands of Matrices have been built on the master thread.)
2. **Re-entrancy is transitive.** A class qualifies only if its whole `update()`
   call graph is re-entrant, across directories — including
   `SRC/coordTransformation/` and any helper that returns a reference to a static
   buffer. ThreadSanitizer over the allowlisted set is what actually closes this;
   the tables above are an audit, not a proof.

## 6. Measurement

Full tables: **`Ladruno_files/testbed/perf/wp107/RESULTS.md`**. Headline:

- **The gate PASSES in this regime.** Loop A is **51.2 % of step** on the
  `LadrunoQuad -bbar` + `LadrunoSANISAND` deck (6 400 elements, `system Pardiso`),
  against **0.26 %** for G-L3's cluster deck at 540 675 DOF. §1's premise — that
  desktop scale plus a substepping critical-state kernel is a different regime —
  is confirmed, not assumed.
- **Bit-identity holds wherever the loop actually threads.** 12/12 runs at
  1/2/4/8 threads with `maxdiff` exactly `0.000e+00` on the allowlisted elastic
  deck (14 400 elements), and the same on the SANISAND deck before it was
  withdrawn. No tolerance was used anywhere.
- **The payoff is not available**, because the material that makes loop A 51 % of
  step is the one §5.1 shows is not thread-safe. On the allowlisted elastic deck
  loop A is 7.8 % of step, so Amdahl caps the win at **1.08x** — and the idle-box
  re-measure (red-team S4) returns **1.03x at 2 and 4 threads and 0.94x at 8, a
  regression**. The mechanism is sound; the deck decides whether it pays, and the
  deck that would pay cannot use it yet.

  The first version of this line reported **1.11x at 8 threads, above the 1.08x
  ceiling this same paragraph computes**, taken on a box running two sibling
  wp/106 jobs and captioned "speed-ups are lower bounds". Both halves were wrong:
  a number above your own Amdahl ceiling is a measurement defect, and contention
  moves this deck's numbers the other way — idle, 8 threads gets *worse*, because
  per-region fork/join and work-queue overhead (paid on every Newton iteration,
  along with the serial snapshot + O(nEle) re-audit + O(nNode) pre-pass in front
  of the loop) is a fixed tax against a prize that is only 8 % to begin with.
- **The serial path is untouched**: `LADRUNO_OPENMP=OFF` vs `ON` at 1 thread is
  byte-identical on both decks, and the affected pytest suites are 128 passed /
  2 skipped / 2 xfailed — plus the 18 new cases in
  `tests/test_wp107_threaded_update.py`, which are the WP's mutation gate
  (red-team B1): with `LadrunoQuad::shp`/`shpBar` reverted to plain `static` the
  pre-existing 128 stay **green** and the new file goes **red, 7 of 18**.
- **The threaded failure path is unreachable** with today's allowlist and
  therefore untested (red-team S5) — see `testbed/perf/wp107/RESULTS.md` §2.2.
- **MPI binaries refuse outright** (red-team B2): `_PARALLEL_INTERPRETERS` gives
  every OpenSeesMP / OpenSeesPyMP rank a plain `Domain`, so the
  `PartitionedDomain`/`Subdomain` overrides never fire there and an inherited
  `LADRUNO_THREADS` would have oversubscribed the whole job silently.

The standing rule from ADR-75b's correctness protocol applies to every number
reported there: **pin `MKL_NUM_THREADS=1`** for identity runs, or the solver's own
~1 ULP jitter (ADR-75b §3 P-6) masquerades as an assembly race.

## 7. Anti-goals reaffirmed

- OpenMP **by default at runtime** stays an anti-goal (ADR-40). The build flag is
  on; the thread count is 1.
- **No threading of any reducing loop** without both modes of ADR-75b P-3.
- **No threading of vanilla elements** before every fork-owned one is done
  (ADR-75b §5.2).
- **No graph colouring** (ADR-75 §4, ADR-75b §4 — settled, not relitigated).
- **No hybrid MPI+threads** (ADR-75b §11 q6 — still deferred).
