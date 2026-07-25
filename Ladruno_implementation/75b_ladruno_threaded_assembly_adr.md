---
title: ADR-75b — Lane 3, threaded element assembly (shared-memory OpenMP)
project: Ladruno
status: proposed — policy SETTLED (§3 determinism, §4 scatter remedy); stage L3-0 MEASURED (gate passes on lanes A/D, and on lane B only under PARDISO); no threading code authorized yet
priority: medium
owner: nmora
amends: 75_ladruno_sparse_direct_strategy_adr
tags:
  - adr
  - performance
  - openmp
  - threads
  - assembly
  - determinism
  - re-entrancy
  - sub-adr
---

# ADR-75b — Lane 3: threaded element assembly

> The **Lane-3 sub-ADR** that [[75_ladruno_sparse_direct_strategy_adr]] §3 / §6-P4 calls for
> ("Likely its own sub-ADR"). Lanes 1–2 (PARDISO desktop, MUMPS cluster) are done; this is the
> remaining lane and the **only** lever for the fork's frame lanes, where solve is 3–31% of runtime.
> It is also **the highest-risk item in the whole ADR-75 family** (§75 §11.1): threaded assembly is a
> whole-codebase re-entrancy invariant, and one missed `static` gives silent, thread-count-dependent
> wrong answers. Family: ADR-40 (perf program) · ADR-40b (dominance) · ADR-75 (solver lanes) ·
> ADR-75a (P0 trade study).

**This document exists to settle two things before any threading code is written**, because both
constrain every stage after them:

1. **§3 — the ordered-reduction / determinism CI policy.** Threaded FP reduction breaks this fork's
   byte-identical oracle discipline *by construction*. **SETTLED.**
2. **§4 — the implicit scatter remedy.** Confirmed as **freeze-sparsity + atomic-scatter** (the Kratos
   pattern), **not** graph coloring. ADR-75 §4 already made this call against two production codes'
   precedent; this section **records it as settled and adds new code evidence**, it does not
   relitigate it. **SETTLED.**

Everything else here is staging and gates. **No threading code is authorized by this ADR.**

---

## 0. Correction to ADR-75 first — the premise of P4(a) is stale

ADR-75 §3-Lane3 and §6-P4(a) both say the first move is to *build* the `elem.update` instrument
("ADR-40b's #1 gap is that the `elem.update` loop … is unscoped. Measure it first"). **That
instrument was built and shipped on 2026-07-06** — `Domain.cpp:2394` carries
`OPS_PROFILE_SCOPE_DEEP_NAMED(_ops_elemUpd, "elem.update")` and `Domain.cpp:2403` carries the
per-classTag `OPS_PROFILE_ELEMPTR_SCOPE`, with the raw-`Element*` macro variant in
`ProfilerMacros.h:114`. ADR-40b's own two addenda of that date already used it to re-measure lanes A
and E (lane A confirmed at 57.6 µs/ele and ~82% of step; **lane E refuted**, killing ADR-68 T6).

So stage (a) is **not** "build the instrument". It is **"use the instrument that exists to produce the
Lane-3 gate table that has never been produced"** — which is a strictly larger and more useful task,
because the >40% gate ADR-75 §3 imposes has never been evaluated *per threadable loop*, only per
profiler *phase*. ADR-75 §3/§6-P4(a) should be amended to say so; this ADR is that amendment.

---

## 1. Why this lane exists (the Amdahl residual)

P1 measured PARDISO at **1.76× against a predicted ~2.2× Amdahl ceiling** on Lane B — so ~34% of even
the *solver-bound* lane is beyond any solver's reach. On the frame lanes the picture is far more
lopsided: solve is 2.9% (lane A) to 30.6% (lane E), and `update` is 64.5% / 35.4%. **No linear solver
helps there at all.** Lane 3 is the only lever.

There is also a second-order effect that *raises* Lane 3's value and is easy to miss: **making the
solver faster mechanically raises the element fraction.** Lane B was 66.4% `linearSolve` / 30.1%
element against UmfPack. PARDISO `-matrixType 2` is measured 1.94–1.96× faster than UmfPack on that
same lane. Holding element time fixed, the element fraction rises toward ~45–48%, i.e. **PARDISO's win
pushes Lane B *across* the >40% gate that ADR-75 §3 sets for Lane 3.** Quantifying that is stage L3-0's
job (§7); it is stated here as the motivation, not as a result.

### Today's shared-memory baseline: there is none (re-verified)

- `find_package(OpenMP)` appears **nowhere** in the fork's own CMake — the only hits are inside bundled
  third-party trees (`mumps-src`, `OTHER/AMGCL`, `OTHER/SuperLU_DIST_4.3`) plus one *comment* at
  `CMakeLists.txt:568` about `libiomp5md` for `mkl_intel_thread`.
- No `/openmp` or `-fopenmp` flag in `CMakeLists.txt`, `Ladruno_scripts/build.bat`, or `setup_env.bat`.
- **7** `#pragma omp` lines exist in `SRC/`, all in PFEM (5 active-looking, 2 commented out). With no
  OpenMP flag they are **compiler-ignored no-ops** — confirming ADR-75 §1's "7 dead lines".

So the whole shared-memory axis is greenfield, and the only multicore parallelism a serial build has
today is MKL's threaded BLAS *inside* a solver.

---

## 2. The three loops — a code-verified taxonomy (the structural contribution)

"Threaded assembly" is not one loop. There are exactly **three** element loops in the implicit path,
and **they have different determinism, different race structure, and different payoff.** Conflating
them is what makes Lane 3 look uniformly terrifying. Separated, one of them is easy.

| # | Loop | Site | Per-element work | Writes to shared state? | FP reduction? |
|---|---|---|---|---|---|
| **A** | `Domain::update()` | `Domain.cpp:2396-2406` | `theEle->update()` | the `ops_TheActiveElement` global (§5.3) + **node trial state, for a named minority of elements** (§2.1b) | **NO — none** |
| **B** | `IncrementalIntegrator::formTangent()` | `IncrementalIntegrator.cpp:112-125` | `getTangent` → `theSOE->addA` | yes — `A[]` | **yes** |
| **C** | `IncrementalIntegrator::formElementResidual()` | `IncrementalIntegrator.cpp:323-337` | `getResidual` → `theSOE->addB` | yes — `B[]` | **yes** |

### 2.1 Loop A has no floating-point reduction — so threading it is bit-identical (for allowlisted elements, §2.1b)

This is the single most consequential fact in this ADR. `Domain::update()`'s body is:

```cpp
while ((theEle = theEles()) != 0) {
  ops_TheActiveElement = theEle;          // <-- a shared GLOBAL write (§5.3)
  OPS_PROFILE_ELEMPTR_SCOPE(_ops_elemUpd, theEle);
  ok += theEle->update();                 // element-private state for MOST elements -- but a
                                          // named minority also write shared NODE state (§2.1b)
}
```

Crucially, **`Domain::update()` contains no reduction into the SOE** — no `addA`, no `addB`, no shared
accumulator. For an element whose `update()` touches only its own state and its own materials, the
threaded loop performs the identical arithmetic in the identical order. Therefore:

> **Threading loop A introduces no change of floating-point summation order, so for an element with no
> cross-element writes a correct threaded loop A is bit-identical to serial at every thread count — and
> the fork's *existing* byte-identical oracles gate it unchanged, with no new CI mode and no
> ordered-mode cost.**

The determinism problem in §3 is confined to loops **B and C**. And loop A is where the frame-lane
money is — **measured** (`lane3/RESULTS_l3a_update_scope.md`): loop A is **80.79% of step on lane A**
and **54.35% on lane D**, with non-kernel overhead of only 0.63%/4.6%. **The highest-payoff loop is
the one with zero determinism cost.** That inverts the natural assumption and it drives §7's staging.

- `ok += theEle->update()` is an **int** accumulation — order-independent, but still a race. Trivially
  fixed with an integer reduction.
- Loop A is also reached from `Domain::update(double,double)` (`Domain.cpp:2415-2422`), which just
  delegates — one implementation point, both entries covered.

### 2.1b The qualification that matters: a minority of elements DO write shared node state

*(Added by the adversarial review, which found the first draft's "writes ONLY this element's own
state" too strong — `RESULTS_l3a_update_scope.md` §7 finding 2.)*

Absence of an SOE reduction is **not** the same as absence of shared writes. At least **13** element
files write **node trial state** (`setTrialDisp` / `setTrialVel` / `setTrialAccel` / `incrTrialDisp`),
and at least one does it squarely on the loop-A path:

- **`LadrunoRigidBody::update()`** (`:490`) calls `imposeSlaveKinematics()` (`:500`), which writes
  `setTrialDisp` (`:518`, `:536`) and `setTrialVel` (`:540`) on **slave nodes shared with other
  elements** — the ADR-58 C3 slaving mechanism.
- Others: `MP_Joint2D`, `MP_Joint3D`, `Adapter`, `ZeroLengthVG_HG`, and the PFEM family
  (`PFEMElement2DCompressible`, `PFEMElement2DQuasi`, `TaylorHood2D`, `BackgroundMesh`, `BCell`,
  `BackgroundFixData`, `PFEMMesher2D/3D`).

This is a **correctness/ordering race, not an FP-order one**: another element reading that node may
see the pre- or post-write value depending on scheduling, so results become nondeterministic in a way
no reduction policy can fix. Consequences:

1. Loop A's bit-identicality is a property of the **allowlisted element set**, not of the loop. The
   per-classTag allowlist (§7) is therefore **load-bearing, not defence-in-depth**.
2. The above classTags are **hard exclusions** for loop-A threading until specifically designed for.
3. Contact is *not* a loop-A exclusion: it is an analysis handler
   (`SRC/analysis/handler/LadrunoContactFE.cpp`, `LadrunoContactHandler.cpp`) over a shared
   `SRC/domain/contact/LadrunoContactDomain.cpp`, so it enters through the **FE_Element** path — an
   exclusion for loops **B/C** instead.
4. Related, and **not** asserted as safe: the draft reasoned that materials are `getCopy()`'d per
   element so cannot alias. True for the standard library, but not exhaustively verified — so it is an
   explicit L3-2 verification item (§11), not a premise.

### 2.2 Loops B and C: the sparsity is **already frozen** — verified, not argued

The Kratos remedy in §4 is "freeze the CSR sparsity, then atomically scatter into it." In OpenSees
**the freeze already exists by construction**, which is why §4's choice is much cheaper here than the
precedent survey implied. Checked in `PARDISOGenLinSOE`:

- `colA[]` / `rowStartA[]` are built **once** in `setSize()` from the DOF graph
  (`PARDISOGenLinSOE.cpp:221-297`) and are not written again anywhere in the assembly path.
- `zeroA()` (`:~305`) zeros **only** `A[]` and clears `factored`. It does **not** touch `colA` /
  `rowStartA`.
- `addA` (`:370-478`) locates each target with a **read-only** linear search over the frozen row
  (`for k in [rowStartA[row]-1, rowStartA[row+1]-1) if colA[k] == col+1`) and then does exactly one
  read-modify-write: **`A[k] += m(i,j)`**.

So the entire race in loop B is **one `+=` on a shared `double`**, at an index computed from
immutable data. Making it `#pragma omp atomic` (or a relaxed `fetch_add` on a `double` via
compare-exchange) is a *localized* change to one line per SOE — no coloring, no thread-private
matrices, no reallocation. `addB` (`:481+`) is the same shape but simpler: `B[pos] += v(i)` with no
search at all.

Two things follow, both material:
1. **§4's decision is not just precedent-backed, it is cheap here.** The expensive half of the Kratos
   pattern (building and maintaining a frozen graph) is pre-paid by `setSize()`.
2. **The `addA` inner linear search is O(rowlen) per entry**, making the scatter O(idSize² × rowlen).
   That is a pre-existing serial inefficiency independent of threading (and exactly the shape the
   fork's own banked lesson warns about: an O(nnz × rowlen) loop where O(nnz) would do).
   **L3-0 measured it: 1 699.2 ms of an ~11.9 s step — ~14% of wall — on the default `-matrixType 0`
   path, and 1.28× SLOWER than UmfPack's scatter.** `-matrixType 2` cuts it 65% (to 591.0 ms) simply
   by scattering half the entries over half-length rows, which is a third measured benefit of
   symmetric storage that P1d never counted. So this is a **separate, serial, zero-determinism-cost
   optimization candidate that on present evidence outranks threading loops B/C** (§11 q5) — not
   Lane-3 work.

---

## 3. DECISION 1 — the determinism / ordered-reduction CI policy · **SETTLED**

### The problem

Floating-point addition is not associative, so any threaded reduction over elements sharing a DOF
changes the summation order and the last bits of the answer. This fork's QA discipline is
**byte-identical / 1e-12 oracles** — Zone-A CI as the C++ no-regression gate, md5-over-the-full-
displacement-field gates (ADR-40b's `-factorOnce` A/B), a 119-test contact battery, and the P1/P1d
solver work reporting **rel err 0.0**. Threaded assembly breaks that by construction (ADR-75 §11.2).

### The precedent, and why this fork must invert its default

LS-DYNA ships this as an explicit consistency switch: `ncpu=-N` selects thread-count-independent
accumulation order at a documented **~10–15%** cost, and the **fast/non-deterministic mode is the
default**. ADR-75 §4/§7 proposed adopting "the same shape."

**Adopt the shape; invert the default.** LS-DYNA can default to fast because its QA is
tolerance-based commercial validation. This fork's QA is *exact*. If threaded assembly defaulted to
fast mode, then every existing oracle test would become thread-count-dependent and the CI signal
would degrade from "bit-identical" to "within tolerance" across the whole suite — losing precisely
the diagnostic power that caught real defects in P1d (a **1-ULP** anomaly on one row was visible and
explicable *because* every other row was 0.0). That loss is global and irreversible; the ~10–15% is
local and recoverable.

There is also a **local precedent already set in this ADR family**: MUMPS **BLR** is an approximate
factorization and was shipped **opt-in, off by default, and explicitly off-limits for the
byte-identical/oracle paths** (ADR-75 §5, §12). Non-deterministic reduction is the same category of
thing — a speed/exactness trade — and gets the same treatment. Deciding otherwise would make the
fork's two "inexactness" features follow opposite rules.

### The policy

**P-1 — Threading stays off by default, per ADR-40's standing anti-goal.** Nothing below weakens
"OpenMP-by-default is an anti-goal." A threaded run is always explicitly requested.

**P-2 — The policy is per loop class, because the determinism cost is not uniform (§2.1).**

| Loop class | Loops | Threaded result vs serial | CI gate | Ordered/fast distinction |
|---|---|---|---|---|
| **Class A** — no FP reduction | `Domain::update` | **bit-identical** (for allowlisted elements — §2.1b) | the **existing** byte-identical oracles, unchanged | **does not exist / not needed** |
| **Class B** — FP reduction into the SOE | `addA`, `addB` | order-dependent | see P-3 | **required** |

**P-3 — For class B: ordered mode is the CI gate AND the default-when-threading-is-on; fast mode is
opt-in.**
- **Ordered mode** must be **bit-identical to serial at every thread count**. It is what CI runs, and
  it is checkable with the existing oracles unmodified — a *stronger* and *cheaper* gate than a
  tolerance.
- **Fast mode** (atomic scatter, §4) is opt-in via an explicit token, carries a documented cost
  delta, is gated at 1e-12 relative, and is **forbidden on oracle / byte-identical / regression
  paths** — the same rule BLR already lives under.
- The flag is a **loud explicit token**, not build-dependent or auto-selected — consistent with
  ADR-75 §5.2's "explicit verbs, not `-auto` magic" decision.

**P-4 — An ordered variant must exist from day one** (ADR-75 §4 lesson 2). §4.2 specifies what it
actually is, because "we'll add ordering later" is how this decision gets silently lost.

**P-5 — One coordinated thread-count policy.** MKL solver threads × OpenMP assembly threads × MPI
ranks will oversubscribe and make every bench lie (ADR-75 §11.4). Assembly threading must not ship
without a single documented knob and a nested-threading rule.

### The one thing that could force a revision — stated up front

Ordered mode being **bit-identical to serial** requires each target scalar to receive its
contributions in *serial element order*, which is a sequential dependency per entry. §4.2's gather
design delivers that, but at a memory cost that scales as Σ(idSize²). **If that memory gate fails at
production solid scale, then "ordered" can only mean *thread-count-independent* (not
*serial-identical*), and the class-B CI gate would have to drop from byte-identical to 1e-12.** That
would be a real weakening of the discipline and it must be taken as an **explicit, separately argued
decision** — not absorbed quietly as an implementation detail. Measuring that memory cost is a gate
in §7 (L3-4), and it is the top open question in §11.

---

## 4. DECISION 2 — the implicit scatter remedy · **SETTLED: freeze-sparsity + atomic-scatter**

**Not relitigated.** ADR-75 §4 evaluated Abaqus, Kratos and LS-DYNA and chose
freeze-sparsity + per-entry atomic scatter (Kratos's proven answer) over graph coloring, explicitly
superseding the perf-skill's coloring-first sketch. That call stands. Recorded here as settled, with
two additions from this session's code reading:

### 4.1 New evidence: the choice is cheaper here than the precedent implied — **fast mode**

§2.2 verified that OpenSees **already freezes the sparsity** (`colA`/`rowStartA` built once in
`setSize`, untouched by `zeroA`), and that the entire race reduces to **one `A[k] += m(i,j)`**. So the
fast-mode implementation is: mark that one statement atomic, per SOE. The expensive half of the
Kratos pattern is pre-paid.

**Hard rule inherited from Kratos, restated because it is easy to violate:** never reallocate or
resize during the threaded phase — only atomic-increment existing scalars. Any path that could grow
`A` mid-assembly (or a `setSize` triggered by a mid-run `domainChanged`) must force the threaded phase
off, not race with it.

**Coloring stays out**, per ADR-75 §4, and is reconsidered only if atomics are *measured* to contend
badly. LS-DYNA-style element-ordering-into-conflict-free-groups remains the stronger-but-larger
option, also demand-gated.

### 4.2 What "ordered mode" concretely is — **gather, not scatter**

P-4 requires an ordered variant from day one, so it needs a real design, not a placeholder. The two
modes turn out to map cleanly onto the two classical assembly algorithms:

| Mode | Algorithm | Parallel over | Race | Result vs serial | Cost |
|---|---|---|---|---|---|
| **fast** | **scatter** + `atomic` (Kratos) | elements | atomics on `A[k]` | order-dependent | low memory; atomic contention |
| **ordered** | **gather** over the same frozen CSR | **matrix rows / nonzeros** | **none — each entry has exactly one writer** | **bit-identical to serial** | Σ(idSize²) storage for element matrices |

Gather is race-free *without atomics* (a thread owns entries, not elements) **and** exact (it sums each
entry's contributions in ascending element order, i.e. serial order). Its price is that all element
matrices must be live simultaneously, since each is read by every row it touches.

**And that price is largely already being paid by the fix for §5.** De-statication done *correctly* —
giving each element **its own** tangent/residual buffer rather than a `thread_local` static — is
simultaneously (i) the re-entrancy fix, and (ii) exactly the per-element storage gather needs. A
`thread_local` static fixes the race but yields only `nthreads` buffers, which does **not** enable
gather. **This is a real design constraint on stage L3-2 and the reason to prefer per-element buffers
over thread-local ones**, even though thread-local is the smaller diff.

Order-of-magnitude for the gather storage, so the gate is concrete: Σ(idSize²) doubles. For Lane B
(3375 `LadrunoBrick`, idSize 24) that is 3375 × 576 × 8 B ≈ **15.5 MB** — trivial. The scaling is
linear in element count at fixed element type, so a ~1M-DOF hex mesh (~300k elements) lands near
**1.4 GB**, which is the same order as the factorization itself and therefore **not** obviously
affordable. Hence the L3-4 memory gate and §11's open question. (These are arithmetic projections from
the element geometry, not measurements.)

---

## 5. The re-entrancy inventory — measured magnitude, and the finding that reorders the plan

ADR-75 §11.1 rates "threaded assembly is a whole-codebase re-entrancy invariant" as the highest risk
in the family. This section puts numbers on it. All counts are `grep` over the current tree
(bdf8adf9f) and are **upper bounds on the hazard**, not all-fatal: many are in `sendSelf`/`recvSelf`
(serialization, off the threaded path). They size the audit, they do not each represent a defect.

### 5.1 The magnitude

- **~5,600** function-/file-scope `static Matrix|Vector|ID` declarations across **587 distinct files**
  in `SRC/element/` + `SRC/material/` (≈1,650 `static Matrix`, ≈3,530 `static Vector`, ≈370
  `static ID`).
- **711** class-level `static Matrix|Vector|ID` members declared in `SRC/element/*.h`.
- `Element` itself owns **static pools shared by every element instance** —
  `Element::theMatrices`, `Element::theVectors1`, `Element::theVectors2`, `Element::numMatrices`
  (`Element.cpp:49-52`). These are the "FE_Element/DOF_Group pools" ADR-40b already flagged as
  blocking rank 7.

**This is not an occasional bug — it is the standard OpenSees idiom.** `getTangentStiff()` returns
`const Matrix &`, so the conventional implementation returns a reference to a function-scope
`static Matrix`. `LadrunoBrick::getTangentStiff` (`:453`) and `getResistingForce` (`:682`) both do
exactly this, over ~12 statics including `stiffJK`, `dd`, `BJ`, `BJtran`, `BK`, `BJtranD`, `shp`,
`Shape` (`LadrunoBrick.cpp:534-544`, `:691`, `:718-721`). Under threading, two elements in
`getTangent` concurrently corrupt each other's buffer — **and the caller then does
`addA(*eleTangent, …)` on a buffer another thread is mid-write on.** Silent, thread-count-dependent
wrong answers: exactly risk §75 §11.1, and the reason ThreadSanitizer + per-classTag gating are
non-negotiable in §7.

### 5.2 The finding that reorders the plan: payoff and threadability are **inverted**

`ForceBeamColumn2d::update()` begins at `ForceBeamColumn2d.cpp:1235`. ADR-40b measured it at
**57.6 µs/ele (~200× its own `getTangent`)**, booking **~82% of lane-A step time**. It is the single
largest measured element cost in the fork. Its interior Newton runs entirely on **function-scope
statics declared inside `update()` itself** — `dv` (:1248), `vin` (:1255), `vr` (:1268), `f` (:1269),
`I` (:1271), `dSe` (:1281), `dvToDo` (:1282), `dvTrial` (:1283), `SeTrial` (:1284), `kvTrial` (:1285),
`Ss`/`dSs`/`dvs`/`fb` (:1342-1345) — 53 statics in the file overall.

And it is a **vanilla upstream file**, which ADR-40b already ruled out of the fork's change budget
("Lane A optimization — out of fork budget for now — vanilla element, vanilla file … upstream-facing
work").

> **So the loop with the highest element fraction is the loop whose element is least threadable, and
> the fix lives in code the fork has decided not to own.** Lane-3 staging therefore cannot be ordered
> by element fraction alone (as ADR-75 §3's single >40% gate implies). It must be ordered by
> **fork ownership × re-entrancy cost × measured fraction**, which puts fork-owned
> `LadrunoBrick`/`LadrunoJ2` (lanes B and D) first and `ForceBeamColumn2d` (lane A) last, or
> upstream.

### 5.3 `ops_TheActiveElement` — a shared global *inside* the loop we most want to thread

Not a `static` in a kernel, so a "de-static the kernels" audit would miss it: `ops_TheActiveElement`
is a **file-scope global** (`SRC/element/Element.cpp:47`, declared `extern` in `SRC/G3Globals.h:45`)
written **per element inside loop A** at `Domain.cpp:2401`, and also in `Element`'s constructor
(`Element.cpp:65`), `Domain.cpp:461`, `OpenSeesCommands.cpp:2865`, and several Ladruno elements
(`LadrunoDispBeamColumn2d.cpp:528`, `3d.cpp:651`).

It is **read by materials** to latch a regularization characteristic length on first
`setTrialStrain`: `LadrunoJ2.cpp:353`, `LadrunoConcrete3D.cpp:347`, `ASDConcrete3DMaterial.cpp:1614`,
`LadrunoRCConcrete.cpp:330`, `LadrunoRCFiniteStrain.cpp`, plus documented reliance in
`BezierTet10.cpp`, `BezierTri6.cpp`, `LadrunoUP.cpp`. It is an **established fork idiom** ("the Phase-3b
lch latch"), not an accident.

Threading loop A races on it, and the failure mode is the worst kind: an element's material latches
**another element's** characteristic length, producing a converged, plausible, **wrong** regularized
softening response. Remedy is `thread_local` (or plumbing the element down the call chain), and it
touches **vanilla files** (`Element.cpp`, `G3Globals.h`) ⇒ `LEDGER_vanilla_files` rows. It is a
prerequisite for L3-1, and it is the concrete proof that stage L3-2 is broader than "grep for
`static`".

---

## 6. Decision summary

1. **§3 determinism policy SETTLED:** per-loop-class. Class A (no reduction) → bit-identical,
   existing oracles gate it, no ordered/fast split. Class B (reduction) → **ordered mode is the CI
   gate and the default when threading is on; fast/atomic mode is opt-in, 1e-12-gated, and forbidden
   on oracle paths** (the rule BLR already lives under). LS-DYNA's *shape*, **inverted default**,
   because this fork's QA is exact rather than tolerance-based.
2. **§4 scatter remedy SETTLED (not relitigated):** freeze-sparsity + **atomic scatter** for fast
   mode — and the freeze is **already there** (`colA`/`rowStartA` survive `zeroA`), so the race is one
   `A[k] +=`. **Coloring stays out**, demand-gated on measured atomic contention. Ordered mode is
   **gather over the same frozen CSR** — race-free without atomics and exact.
3. **Staging is ordered by fork ownership × re-entrancy cost × measured fraction, not by fraction
   alone** (§5.2). Loop A on fork-owned elements first — highest payoff, zero determinism cost.
4. **De-statication must produce per-element buffers, not `thread_local` statics** (§4.2), because
   only per-element storage enables the exact gather mode.
5. **No threading code is authorized yet.** Every stage is gated (§7).

---

## 7. Stages and gates

Each stage lists its **gate** — the measured condition that must hold *before* it starts. ADR-75 §3's
standing gate (>40% element/`update` fraction, because trivially cheap elements can *regress* from
fork-join overhead) applies throughout and is refined per stage.

- **L3-0 — MEASURE. The Lane-3 gate table. ✅ DONE** (`Ladruno_files/testbed/perf/lane3/RESULTS_l3a_update_scope.md`).
  **Headline: the Lane-1 solver win is what makes Lane 3 worth doing.** Same Lane-B model and binary,
  `system Pardiso -matrixType 2` @4T vs UmfPack: element-kernel fraction **35.8% → 74.9%**, solve
  **55.9% → 16.9%** — the >40% gate **fails under UmfPack and passes with huge margin under PARDISO**.
  Lanes A and D are element-bound regardless (**81.9%** and **86.8%** kernel). Loop A alone is
  **80.8%** of lane A and **54.4%** of lane D, with non-kernel overhead of 0.63%/4.6% and no scatter.
  Kernel time is solver-invariant to 2.3% (validity check). `ForceBeamColumn2d::getTangent` measured
  at **0.12 µs/ele**, confirming ADR-75 §3's regression hazard — though the decisive reason not to
  thread lane A's loops B/C is that they are **0.45%/0.79% of step**. Instrumentation tax measured
  (deep vs coarse, 6 rounds): **+3% to +8%**, moving every kernel fraction by ≤1.1 pp. Two defects
  found and fixed by the adversarial review (§2.1b; the F3 re-attribution).
  *(Original scoping:)*
  Use the shipped `elem.update` / `elem.tangent` / `elem.residual` instruments to produce, per lane:
  (i) wall fraction inside **each** of the three loops of §2, separately — never a phase total;
  (ii) the **kernel vs scatter** split within loops B/C (`elem.tangent` vs `formTangent − elem.tangent`),
  which is the Amdahl input for §4's remedy;
  (iii) **mean µs per element per loop** and the element count, against a fork-join/barrier cost — this
  is what decides regression risk, and a fraction alone cannot;
  (iv) the resulting per-loop Amdahl ceiling at 2/4/8 threads;
  (v) Lane B **re-measured under PARDISO**, to test §1's claim that the solver win pushes the element
  fraction across 40%.
  **Deliverable:** `Ladruno_files/testbed/perf/lane3/RESULTS_l3a_update_scope.md`.
  **Gate: none — this is the gate for everything else.**

- **L3-1 — thread loop A (`Domain::update`) for an allowlisted set of fork-owned classTags.**
  Highest payoff, **zero determinism cost** (§2.1). **Target Lane D first** — L3-0 makes it the
  cleanest entry: loop A is 54.4% of step on the fork-owned `LadrunoBrick`, `system Diagonal` removes
  the SOE from the picture entirely (`linearSolve` 0.03%), and de-statication there is already
  authorized (ADR-40a C16 / ADR-68 T3). Needs: `ops_TheActiveElement` → `thread_local` (§5.3), an
  integer reduction on `ok`, a per-classTag **opt-in allowlist** (default empty) that **excludes the
  §2.1b node-writing classTags**, and `find_package(OpenMP)` + a build flag.
  **Gate: ✅ the fraction half is PASSED** — L3-0 measures loop A at 54.4% (lane D) and 80.8%
  (lane A) with per-element update cost 3.20 and 20.53 µs/ele. **Still open:** the barrier cost is
  unmeasured (no OpenMP in the build yet — §11 q3), and the allowlisted classTags must pass §7's
  correctness protocol.
  **Acceptance: bit-identical to serial at 1/2/4/8 threads.** Anything less means a re-entrancy miss,
  not a rounding difference — this is the stage where that gate is free, so it must be enforced
  absolutely.

- **L3-2 — de-static the allowlisted fork-owned kernels into per-element buffers.**
  `LadrunoBrick` first (ADR-40b already authorized "de-static the brick scratch" independent of
  threading, as the ADR-40a C16 carve-out), then `LadrunoJ2`/`LadrunoConcrete3D`. **Per-element, not
  `thread_local`** (§4.2). Explicitly **not** `ForceBeamColumn2d` (§5.2).
  **Gate:** serial performance and answers unchanged (this is a refactor); ThreadSanitizer clean on the
  allowlisted set.

- **L3-3 — explicit / diagonal-SOE residual reduction (Lane D).**
  ADR-75 §3's original "start here" — private-buffer reduction with no factorization interaction.
  Lane D is ~89% element work once fully attributed, so the fraction gate is passed with margin.
  Re-ranked *after* L3-1/L3-2 because loop A is a larger, cheaper win and because Lane D's residual
  path shares the same brick kernels L3-2 fixes.
  **Gate:** L3-2 done for the brick; ordered variant designed per P-4.

- **L3-4 — implicit assembly, loops B and C: fast (atomic scatter) + ordered (gather).**
  The §4 work. **Must ship both modes together** — shipping fast-only would silently make P-3
  undecidable in practice.
  **Gate:** L3-0 shows the B/C kernel fraction justifies it (scatter overhead is *not* threadable
  payoff); **and the §4.2 gather memory measured on a production-scale deck**, since that number
  decides whether the byte-identical class-B CI gate survives (§3's revision trigger, §11).

- **L3-5 — thread-count policy.** One coordinated knob across MKL solver threads × assembly threads ×
  MPI ranks, plus the nested-threading rule (ADR-75 §11.4). **Gate:** must land with or before L3-1;
  benches before it are not trustworthy.

### Correctness protocol (applies to every stage that touches threading)
1. **ThreadSanitizer** on the allowlisted set — the only tool that finds the misses `grep` cannot.
2. **Per-classTag allowlist, default empty.** An un-audited element is never threaded. This converts
   §75 §11.1's "one miss = silent wrong answers" from a whole-codebase invariant into a per-element
   opt-in.
3. **Threaded path behind a default-off flag** so the serial path stays byte-identical (§75 §11.3).
4. **Bit-identical at 1/2/4/8 threads** for class A and for class-B ordered mode. 1e-12 for class-B
   fast mode, never on an oracle path.

---

## 8. Anti-goals
Inherited from ADR-40/75 and still standing: **OpenMP-by-default** · **implicit colored-scatter before
a measured >40% element fraction** · **graph coloring** as the first remedy (§4; demand-gated on
measured atomic contention only) · SIMD before algorithmic work-removal · GPU offload.
Lane-3-specific additions:
- **Threading `ForceBeamColumn2d` (or any vanilla element) before every fork-owned element is done**
  — wrong order on both risk and change-budget grounds (§5.2).
- **`thread_local` statics as the de-statication strategy** — fixes the race, forecloses exact gather
  (§4.2).
- **Shipping class-B fast mode without ordered mode** — makes §3 unenforceable (§7 L3-4).
- **Threading any loop before L3-5's thread-count policy exists** — the benches would mislead
  (§75 §11.4).

## 9. Risk register (Lane-3 specific; refines ADR-75 §11)
1. **Silent re-entrancy miss** — ~5,600 statics / 587 files (§5.1). *Highest.* Mitigated by the
   per-classTag allowlist (default empty) + ThreadSanitizer + bit-identical acceptance, which together
   make a miss *loud* rather than silent for the elements actually threaded.
2. **`ops_TheActiveElement` cross-latch** (§5.3) — a shared global inside loop A, read by 5+ materials
   for regularization length. Failure mode is a converged, plausible, wrong softening response. A
   `static`-only audit misses it entirely.
2b. **Elements that write shared node trial state inside loop A** (§2.1b) — `LadrunoRigidBody` and ≥12
   others. This is an **ordering** race, so unlike everything else in this register it cannot be fixed
   by any reduction policy; only exclusion or redesign works. It is also the finding that makes the
   §7 allowlist load-bearing rather than belt-and-braces. Caught only by the adversarial review, after
   the first draft asserted loop A had no shared writes at all.
3. **The byte-identical class-B CI gate may not survive the gather memory cost** (§3's revision
   trigger, §4.2). This is the decision most likely to need reopening, so it is flagged now rather
   than discovered at L3-4.
4. **Reward remains back-loaded onto the riskiest lane** (§75 §11.7) — the frame-lane payoff lives
   entirely here. §2.1 *reduces* this: loop A is both the biggest frame-lane win and the cheapest,
   determinism-free stage, so the lane has a real early win rather than only a hard finale.
5. **Nested threading / oversubscription** (§75 §11.4) → L3-5 is a hard prerequisite.
6. **Assembly loop is a central chokepoint** shared by static/transient/eigen/sensitivity
   (§75 §11.3) → default-off flag, serial path untouched.
7. **Mid-run `domainChanged`** (contact re-emission ADR-60, element removal ADR-51) can resize the SOE
   and invalidate the frozen sparsity §4 depends on. Must force threading off, not race — same family
   as the `-factorOnce` staleness caveat ADR-40b already banked.

## 10. Ledger / banner plan
**Nothing under `SRC/` is touched by this ADR or by L3-0** (policy, staging, and measurement only):
- `LEDGER_vanilla_files.md` — **no row due.** No upstream file modified. (L3-1 will owe several; see
  below.)
- `LEDGER_implementations.md` — **no row due.** The only new files are perf-testbed tooling
  (`Ladruno_files/testbed/perf/lane3/{parse_lane3.py,sweep_l3a.sh}`) plus env-override knobs in the
  existing phase-0 lane scripts. That follows the standing precedent: the phase-0 and phase-1
  harnesses have no implementation rows either — the ledger tracks engine features, and the profiler
  itself already has its row.
- `LEDGER_quirks.md` — **five rows added in this PR**: the `static`-return-buffer non-re-entrancy
  idiom; `ops_TheActiveElement`; the named-deep-scope-includes-the-scatter trap; the already-frozen CSR
  sparsity; and ADR-40b's stale lane-D `formTangent` number.
- `banner_features.txt` — **not now.** Only once a threading feature is user-visible and `shipped`.

When stages land:
- `LEDGER_vanilla_files.md` — `Domain.cpp` (loop A), `Element.cpp` + `G3Globals.h`
  (`ops_TheActiveElement` → `thread_local`), `IncrementalIntegrator.cpp` (loops B/C), `CMakeLists.txt`
  (`find_package(OpenMP)`), each with a `// Ladruno` marker.
- `LEDGER_implementations.md` — one row per threading feature + its flag token.
- `LEDGER_quirks.md` — this session's entries are due now (§0's stale-premise trap, the
  `static`-return-buffer idiom, `ops_TheActiveElement`, the frozen-sparsity fact).
- `Ladruno_scripts/banner_features.txt` + `patch_banner.py` — only once a threading feature is
  user-visible and `shipped`. **Not now.**

## 11. Open questions
1. **Does the §4.2 gather memory (Σ idSize² doubles) fit at production solid scale?** ~15.5 MB at Lane
   B, ~1.4 GB projected at ~1M DOF. **This single number decides whether the byte-identical class-B CI
   gate of §3 survives** — the highest-leverage unknown in this ADR. Gate at L3-4; measure on a real
   deck.
2. **Does atomic `A[k] +=` actually contend?** If yes, §4's fallback ladder (LS-DYNA element ordering,
   then coloring) activates. Cheap to measure once L3-4 exists; do not pre-optimize.
3. **What is the real fork-join/barrier cost on this toolchain**, in µs? Sets the per-element-cost
   threshold that decides regression risk in every stage's gate. Measure once, reuse.
4. **How much of the ~34% Amdahl residual is loop A vs loops B/C vs neither?** L3-0 answers this; if
   the residual is mostly *neither*, Lane 3's whole ceiling is lower than assumed and the staging
   should shrink.
5. **Replace the `addA` O(idSize² × rowlen) inner search — PROMOTED, now with a number.** L3-0
   measured `PARDISOGenLinSOE::addA` at **1 699.2 ms of an ~11.9 s step (~14% of wall) on the
   default `-matrixType 0` path**, and it is **1.28× slower than UmfPack's** scatter (1 322.8 ms) —
   the *slowest* of the three configurations measured. That is 14% of wall spent locating entries in
   a structure unchanged since `setSize`. A serial, zero-determinism-cost, zero-threading fix on the
   shipped desktop path. **On present evidence this is a better next move than threading loops B/C at
   all** — and it is Lane-1 work, not Lane-3.
6. **Cluster composition:** loop-A threading is the "threads-per-node half" of hybrid MPI+threads
   (ADR-75 §3). Untested interaction with `PartitionedDomain`; defer until L3-1 lands serially.
7. **Does any element alias (rather than `getCopy()`) a material or section across elements?** Not
   verified (§2.1b item 4). If one does, threaded loop A corrupts shared material state. An L3-2
   verification item; ThreadSanitizer over the allowlisted set is the practical check.
8. **Does the F1 gate verdict survive at production size?** Lane B here is 11 520 DOF, and P1c
   measured the solve fraction *growing* with N. At production scale the element fraction is lower
   than the measured 74.85%, so the margin shrinks. Re-check before committing L3-3/L3-4. Not in
   doubt for lanes A and D (81–87%, no solver in the picture).
