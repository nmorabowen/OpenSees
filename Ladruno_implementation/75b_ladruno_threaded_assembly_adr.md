---
title: ADR-75b — Lane 3, threaded element assembly (shared-memory OpenMP)
project: Ladruno
status: PARKED — measurement-gated on G-L3 (§12, answers §11 q8). Policy SETTLED (§3 determinism, §4 scatter remedy); stage L3-0 MEASURED and REVISED after a 64-agent adversarial review (Lane B passes only under PARDISO; Lane D's loop A FAILS the gate once `-commitSolveState` is on ⇒ L3-1 de-authorized). L3-0b work-removal now SHIPPED in full ⇒ every remaining ranked target is a REDUCING loop, so the cheap reduction-free entry point is gone. Neither threading code NOR prerequisite work authorized until G-L3 measures the element fraction at production scale (≥500k DOF)
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

## 2. The FIVE assembly loops — a code-verified taxonomy (the structural contribution)

"Threaded assembly" is not one loop. There are **five** assembly loops, and **they have different
determinism, different race structure, and different payoff.** Conflating them is what makes Lane 3
look uniformly terrifying. Separated, one of them is easy.

> **Corrected by the max-effort review (M1).** The first version of this table listed **three** loops
> and gave `IncrementalIntegrator.cpp` as the sole implementation point. Both were wrong: the
> **transient** path (which is Lane D — this ADR's own recommended first target) assembles in
> `TransientIntegrator.cpp`, and **both** paths carry *DOF_Group* loops that this taxonomy omitted.
> An earlier review pass had explicitly cleared this ("`formTangent` has ONLY the FE_Element loop");
> that clearance read the static path only and is **withdrawn**. Consequences in §4.2 — a gather keyed
> on element matrices alone silently drops nodal mass and nodal loads.

| # | Loop | Site | Per-iteration work | Shared writes? | FP reduction? |
|---|---|---|---|---|---|
| **A** | `Domain::update()` | `Domain.cpp:2396-2406` | `theEle->update()` | `ops_TheActiveElement` (§5.3); node trial state for 2 elements (§2.1b); the shared **iterator cursor** (§5.4) | **NO — none** |
| **B** | element tangent — **two implementations** | static: `IncrementalIntegrator.cpp:112-125`  ·  **transient: `TransientIntegrator.cpp:111-125`** | `getTangent` → `addA` | `A[]`, **and `FE_Element::theTangent`, a class-wide pool (§5.4)** | **yes** |
| **C** | `IncrementalIntegrator::formElementResidual()` | `:323-337` | `getResidual` → `addB` | `B[]`, same pool | **yes** |
| **D** | **DOF_Group tangent** (nodal mass/damping) | `TransientIntegrator.cpp:99-107`, preceded by `addModalDampingMatrix` `:90-95` | `dofPtr->getTangent` → `addA` | `A[]` | **yes** |
| **E** | **DOF_Group unbalance** (nodal loads) | `IncrementalIntegrator::formNodalUnbalance` `:290-309` — runs on **both** paths | `dofPtr->getUnbalance` → `addB` | `B[]` | **yes** |

Loops **D** and **E** are small in wall terms (measured from the committed data: 192.0 ms on lane D
= 0.65% of step; 12.6 ms lane B/UmfPack) so they do not move any §1 fraction — but they are **not
optional for correctness**. `DOF_Group::addMtoTang` (`DOF_Group.cpp:351-358`) adds
`myNode->getMass()`, so on any deck using `ops.mass(...)` the **entire nodal mass reaches `A` only
through loop D**. The fork's own profiler comment at `TransientIntegrator.cpp:108-110` says as much
("the DOF_Group mass/damping loop above is not element-keyed and is left untimed here") — the
information was in the tree and this ADR still missed it.

**Serial composition order**, which any ordered/gather mode must reproduce:
`A` ← (modal damping → all DOF_Groups → all FE_Elements); `B` ← (all FE_Elements → all DOF_Groups).

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

> **Scoped down by a second review pass.** The first version of this section said "at least **13**
> element files write node trial state" and treated all 13 as loop-A exclusions. That was a grep-level
> claim presented as an allowlist: it counted files under `SRC/element/` regardless of whether the
> class is an `Element` at all, and regardless of *which method* does the writing. Resolving both
> shrinks the loop-A hazard from 13 to **2** — which is materially good news for L3-1, and would have
> mis-scoped the allowlist if left as written.

Absence of an SOE reduction is **not** the same as absence of shared writes. Grepping
`setTrialDisp`/`setTrialVel`/`setTrialAccel`/`incrTrialDisp` under `SRC/element/` hits 13 files;
resolving each to its class and its enclosing method gives:

| write site | class | reachable from | loop-A hazard? |
|---|---|---|---|
| `LadrunoRigidBody::update()` `:490` → `imposeSlaveKinematics()` `:500` (writes `:518`, `:536`, `:540`) | `Element` | **`update()`** | **YES** |
| `ZeroLengthVG_HG` `:631` | `Element` | **`update()`** | **YES** |
| `Adapter` `:588`, `:601` | `Element` | `getResistingForce()` | no — **loop C** |
| `PFEMElement2DCompressible` `:222`, `PFEMElement2DQuasi` `:200` | `Element` | `commitState()` | no — commit phase |
| `TaylorHood2D` `:599-601` | `Element` | `setDomain()` | no — setup |
| `MP_Joint2D`, `MP_Joint3D` | **`MP_Constraint`, not `Element`** | constraint handler | no — never in the element loop |
| `BackgroundMesh`, `BCell`, `BackgroundFixData`, `PFEMMesher2D/3D` | **utility/mesher classes, not `Element`** | — | no |

**So the loop-A exclusion list is exactly two elements: `LadrunoRigidBody` and `ZeroLengthVG_HG`**
(plus `Adapter` for loop C, and the two PFEM elements if `commitState` is ever threaded).
`LadrunoRigidBody` is the important one — it writes to **slave nodes shared with other elements** (the
ADR-58 C3 slaving mechanism), so a concurrent reader sees a torn value.

This is a **correctness/ordering race, not an FP-order one**: another element reading that node may
see the pre- or post-write value depending on scheduling, so results become nondeterministic in a way
no reduction policy can fix. Consequences:

1. Loop A's bit-identicality is a property of the **allowlisted element set**, not of the loop. The
   per-classTag allowlist (§7) is therefore **load-bearing, not defence-in-depth**.
2. `LadrunoRigidBody` and `ZeroLengthVG_HG` are **hard exclusions** for loop-A threading until
   specifically designed for. The list is short — but note it was found by grep, and a
   *transitively* reached node write (an element calling a helper that writes) would not show up.
   ThreadSanitizer over the allowlisted set is what actually closes this, not the table above.
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
- `zeroA()` (`PARDISOGenLinSOE.cpp:554-568`) zeros **only** `A[]` and clears `factored`. It does **not** touch `colA` /
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
   path, and 1.28× SLOWER than UmfPack's scatter. ✅ FIXED by ADR-75 P1f
   ([#636](https://github.com/nmorabowen/OpenSees/pull/636)): binary search, 1.09-1.10× end-to-end.
   The search is gone; the `A[k] +=` accumulate — and therefore the race — is unchanged.** `-matrixType 2` cuts it 65% (to 591.0 ms) simply
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

**P-6 — ⚠ THE BASELINE IS NOT BYTE-REPRODUCIBLE WHEN THE SOLVER IS THREADED (P1f, [#636](https://github.com/nmorabowen/OpenSees/pull/636)).**
*Added 2026-07-25 after P1f corrected a claim this policy was resting on.* ADR-75 P1 concluded
"PARDISO is bit-identical at every thread count ⇒ threading introduces no FP drift ⇒ the determinism
concern is **Lane-3-only**". That conclusion came from **one run per thread count**, and it is wrong:
at `MKL_NUM_THREADS=4` a 14³ Lane-B model returns **two distinct tip displacements over 10 runs of
the same binary** (~1 ULP, a 5/5 split); at 1 thread it is 10/10 identical. It is size-dependent —
an 8³ model is reproducible even at 4 threads.

This does **not** weaken §3's policy, but it changes two things about how the policy is *enforced*:

1. **"The existing byte-identical oracles gate it unchanged" is true only with
   `MKL_NUM_THREADS=1`.** Otherwise the *baseline itself* moves ~1 ULP run-to-run, and a threaded-
   assembly bit-identical gate cannot distinguish assembly drift from solver drift — the gate becomes
   unfalsifiable rather than strict. **Every Lane-3 acceptance run pins `MKL_NUM_THREADS=1`.** (The
   real fix, if a threaded reproducible mode is ever wanted, is PARDISO's CNR control `iparm[33]`,
   which *is* available to this fork because we set `iparm[1]=2`, not parallel METIS. Not exposed.)
2. **The determinism problem was never Lane-3-only.** It is a property of every threaded FP
   reduction in the process, and MKL's is already there. Lane 3 adds a *second* source, which is
   exactly why the ordered/fast split of P-3 is worth having — but the ADR should stop implying the
   solver lane is clean.

**And the methodological lesson applies directly to §7's acceptance protocol:** "deterministic" is a
claim about a **distribution**, never about one sample. See §7's protocol item 4, which was written
with the same n=1 flaw and is now fixed.

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

Gather is race-free *without atomics* (a thread owns entries, not elements) **and** exact. Its price
is that all contributions must be live simultaneously, since each is read by every row it touches.

> **Three corrections from the max-effort review — the first draft's exactness argument was
> under-specified and its "storage is free" claim was simply false.**

**(i) The ordering key is a triple, not "element order" (H8).** `addA` issues one independent
`A[k] += m(i,j)` per `(i,j)` pair with **no dedup**, so a single FE_Element can contribute to the
same `A[k]` **more than once** whenever its `ID` repeats an equation number — reachable through
`TransformationFE::setID` + `TransformationDOF_Group::getID` returning the *retained* node's
equations, i.e. any `rigidDiaphragm`/`equalDOF`/`rigidLink` with two element nodes slaved to one
master. So the serial order to replay is **(FE-iteration index, i, j)**, and a gather that
precomputes "element e's total contribution to k" is **not** bit-identical. Also, the loop iterates
**FE_Element** tags, not element tags, and includes constraint FE_Elements (`PenaltySP_FE`,
`PenaltyMP_FE`, `LagrangeMP_FE`) that have **no Element at all** — so the store must be keyed on
FE_Element, not Element. The P-3 acceptance deck must include `rigidDiaphragm` under
`constraints Transformation`, or this case never gets exercised.

**(ii) Element matrices are not the whole of `A` or `B` (M1).** Per §2, loops **D** and **E**
contribute nodal mass, modal damping and nodal loads. A gather that zeroes `A` and reconstructs it
from stored *element* matrices produces a **wrong answer**, not a last-bits difference. Ordered mode
must compose the DOF_Group contributions in their serial position (before the FE contributions for
`A`, after them for `B`).

**(iii) The storage is NOT already paid for by §5's de-statication (H1) — this was the draft's worst
error.** The buffer `addA` actually consumes is **not the element's**: `FE_Element::theMatrices` is a
`static Matrix **` class-wide pool (`FE_Element.cpp:51`, `.h:132`), and
`theTangent = theMatrices[numDOF]` (`:127`, `:137`) hands **one shared 24×24 Matrix to every
same-numDOF FE_Element** — all 3375 bricks of the Lane-B mesh. `IncrementalIntegrator.cpp:118-120`
passes that address straight to `addA`. So de-staticating `LadrunoBrick` gives gather **zero** live
storage, and threading loops B/C is a **100% collision** on that pool that a per-classTag *element*
allowlist cannot mitigate. Under `constraints Transformation` the same holds one level up
(`TransformationFE.cpp`: `modMatrices`, `dataBuffer`, `localKbuffer`, `dofData`).

**Revised conclusion.** Ordered/gather mode requires **FE_Element/DOF_Group de-statication as a hard
prerequisite** — that is ADR-40's rank-7 item, which §5.1 previously mis-cited (see §5.4). The
per-element-vs-`thread_local` preference for *element kernels* still stands for re-entrancy reasons,
but the claim that it comes for free with gather storage is withdrawn; and note that per-element
buffers carry their own **serial** memory cost (§11 q9).

Order-of-magnitude for the gather storage, so the gate is concrete: Σ(idSize²) doubles. For Lane B
(3375 `LadrunoBrick`, idSize 24) that is ≈**15.5 MB** — trivial. Deriving from the deck's actual N³
hex topology, a ~1M-DOF mesh is ~325k elements ⇒ ≈**1.5 GB**, the same order as the factorization
itself and therefore **not** obviously affordable. Hence the L3-4 memory gate and §11's open
question. (Projection from element geometry, not a measurement.)

Order-of-magnitude for the gather storage, so the gate is concrete: Σ(idSize²) doubles. For Lane B
(3375 `LadrunoBrick`, idSize 24) that is 3375 × 576 × 8 B ≈ **15.5 MB** — trivial. The scaling is
linear in element count at fixed element type, so a ~1M-DOF hex mesh (~300k elements) lands near
**~1.5 GB**, which is the same order as the factorization itself and therefore **not** obviously
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
  in `SRC/element/` + `SRC/material/` (≈1,690 `static Matrix`, ≈3,530 `static Vector`, ≈370
  `static ID`).
- **711** class-level `static Matrix|Vector|ID` members declared in `SRC/element/*.h`.
- `Element` owns **static pools shared by every element instance** — `Element::theMatrices`,
  `theVectors1`, `theVectors2`, `numMatrices` (`Element.cpp:49-52`). ⚠ **These are NOT the
  "FE_Element/DOF_Group pools" of ADR-40's rank 7** — the first draft said they were, which is wrong
  and pointed the reader away from the pool that actually matters. `Element::theMatrices` is a
  *third*, distinct pool, reached only from `getDamp`/`getMass`/`getResistingForceIncInertia`, i.e.
  **off the loop-A path**. ADR-40 (`40_ladruno_performance_adr.md:99-102`) names `FE_Element.h:132-133`
  and `DOF_Group.h:150-151` — see **§5.4**, which is the one that races.

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
is a **file-scope global** (`SRC/element/Element.cpp:47`, declared `extern` in **both** `SRC/OPS_Globals.h:71` and `SRC/G3Globals.h:45` — the *former* is the one every reader and `Domain.cpp` actually includes, 329 including files vs 71; best fix is to delete the duplicate declaration)
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

### 5.4 The hazards the inventory's SCOPE structurally could not see (max-effort review)

§5.1 counted statics under `SRC/element/` + `SRC/material/`, and §2.1b grepped for node-trial-state
**setters**. Both scopes are wrong in the same way: the loop-A/B/C call graph leaves those
directories and the write can be a *getter*. Every hazard below produces a **silent,
thread-count-dependent wrong answer with no node write and no FP-order change** — i.e. it is
invisible to §2.1b's allowlist criterion *and* to §3's determinism policy, and it will not reproduce
on every schedule.

| # | Hazard | Where | Why the inventory missed it |
|---|---|---|---|
| **H1** | **`FE_Element::theTangent` is a class-wide pool** — `static Matrix **theMatrices` (`FE_Element.cpp:51`, `.h:132`), `theTangent = theMatrices[numDOF]` (`:127`, `:137`) ⇒ **one 24×24 Matrix shared by all 3375 bricks**, handed straight to `addA` (`IncrementalIntegrator.cpp:118-120`). Same one level up under `constraints Transformation` (`TransformationFE.cpp`: `modMatrices`, `dataBuffer`, `localKbuffer`, `dofData`) | loops B/C | `SRC/analysis/fe_ele/` — outside the scanned dirs. **100% collision**, not mitigated by an element allowlist, by L3-2, or by atomic `A[k] +=` |
| **H2** | **`LadrunoBrick::update()` runs on 16 shared function-scope statics** (`:854 strainE`, `:882 Bbar`, `:883 strainG`, `:909 strainC`, `:910 Bc`, `:911/:965 ulj`, `:930-936 dvol/gaussPoint/strain/shp/Shape/shpBar/BJ`, `:873/:906 shpC`) and `computeLocalDisp()` (`:1202-1206`) returns a static buffer. Corpus scan: **63 of 215** `Element::update()` bodies contain non-const function-scope statics (473 declarations) | loop A | it WAS in scope — but §2.1b's *node-write* list was mis-promoted into a *threadability* list. **This is L3-1's own named first target** |
| **H3** | **The loop cursor is shared and mutable.** `Domain::getElements()` (`Domain.cpp:1564-1569`) resets and returns the single member `theEleIter` (`Domain.h:313`) over `ArrayOfTaggedObjects::myIter`, whose `operator()` does `numDone++; currIndex++`. Same for `AnalysisModel::getFEs()` | A, B, C | the ADR never says how the loop *becomes* a parallel-for. It is **not index-addressable**, so L3-1 needs a snapshot (unbudgeted allocation) or a new random-access accessor on `SRC/tagged/storage/`. A naive shared-cursor pull **skips and duplicates elements** |
| **H4** | **`Node`'s trial-state GETTERS mutate the node** — `getTrialVel/getTrialAccel/getVel/getAccel/getIncrDisp/…` (`Node.cpp:590-706`) lazily `createVel()/createDisp()/createAccel()` on first touch, heap-allocating and assigning members. 41 of 215 `update()` bodies read one | loop A | the grep was over **setters**. Two elements sharing a node both see `trialVel==0`, both allocate ⇒ last-writer-wins, leak, and a live `const Vector&` into a freed buffer |
| **H5** | **`SRC/coordTransformation/` is on loop A's path for every beam-column element** — `LinearCrdTransf2d::getBasicTrialDisp()` (`:311`, `:327`) builds `static double ug[6]` / `static Vector ub(3)` and returns `ub` by reference; called inside `update()` at `LadrunoDispBeamColumn2d.cpp:534`. 434 declarations / 18 files, entirely uncounted | loop A | directory scope. Per-**class** buffer, so the per-element-copy argument does not protect it |
| **H6** | **The profiler instrument inside the loop is itself an unsynchronized shared write.** The named scope is constructed by one thread *outside* the loop (`Domain.cpp:2394`); every `~ElemScope` calls `ProfileNodeLive::addElem` (`Profiler.cpp:115-124`) which lazily builds a `std::map` and does counter RMWs on the master thread's node. `Profiler.h:59-63` states the precondition being violated ("each thread owns its own tree") | A, B, C | it is instrumentation, not physics — so nobody looked. Concurrent `std::map` insertion is UB. **Aggravator: proving L3-1's payoff requires the deep gate armed, i.e. the only configuration where the win is measurable is the UB one** |
| **H7** | **`ops_TheActiveElement` is not a first-touch latch** — `LadrunoJ2.cpp:347-356` (inside `integrate()`, called unconditionally) and `LadrunoConcrete3D.cpp:344-351` read it on **every** call, not just the first. Two readers are missing from §5.3 entirely: `ASDConcrete1DMaterial.cpp:998` and `ASDSteel1DMaterial.cpp:2150` — **uniaxial fiber materials, i.e. the lane-A path**. And the global **persists after the loop**; `thread_local` silently changes that residual from "last element in iteration order" to "last element on the master thread" | loop A | §5.3 characterised it from the guarded readers only. A per-call read is live cross-talk every iteration, not a benign first-touch race |

**Consequences.** (1) **L3-2 becomes a hard prerequisite of L3-1**, not its successor (§7). (2) The
allowlist criterion is **strictly stronger than the node-write list** — an element qualifies only if
its whole `update()` call graph is re-entrant, transitively and across directories. (3) `§10`'s
vanilla-ledger plan grows (see §10). (4) ThreadSanitizer must run **with the deep profiler gate on**,
because that is the configuration the payoff is measured in.

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
  thread lane A's loops B/C is that they are **0.45%/0.79% of step**. Kernel time is solver-invariant to **5.8%**, which is
  this box's stated noise floor. Instrumentation tax measured (deep vs coarse, 6 rounds): **+8.0% on
  lane D**, the only lane where a paired sign test establishes it (6/6, p=0.016; lanes A/B are 4/6 and
  5/6, not separable from noise); it moves every kernel fraction by ≤1.1 pp. **⚠ Lane D's loop-A
  fraction is superseded — see §7 L3-1 and `RESULTS` F9: with the shipped `-commitSolveState` it is
  38.95% and FAILS the gate.** Three adversarial passes found **12** defects in this work
  (§2.1b, §2's missing DOF_Group loops, §4.2's storage claim, §5.4's seven hazards, the F3
  re-attribution, and F9).
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

- **L3-0b — WORK REMOVAL FIRST, before any threading stage. ⚠ NEW, and it outranks everything below.**
  Two shipped, zero-code, zero-risk items that the L3-0 measurement surfaced:
  (i) **`-commitSolveState` on explicit decks** — ADR-67 P-NEW-2, bit-identical on rate-independent
  materials, measured here at **−29.6% wall on Lane D** (24 555 → 17 278 ms, `disp_z` identical to
  all digits in 6 runs); (ii) ~~replace the `addA` O(idSize²×rowlen) search~~ ✅ **SHIPPED as P1f
  ([#636](https://github.com/nmorabowen/OpenSees/pull/636)) — 1.09-1.10× end-to-end**, one day after
  L3-0 named it, and it needed none of Lane 3's machinery. The remaining half of this item is the
  same fix in `UmfpackGenLinSOE::addA`, which has the identical frozen-CSC linear scan (§11 q5). ADR-40's standing order is work-removal before parallelism, and
  both of these are strictly cheaper than any threading stage. Neither is Lane-3 work.

- **L3-1 — thread loop A (`Domain::update`) for an allowlisted set of fork-owned classTags.**
  Still the only **reduction-free** loop, so still the cheapest determinism story and the right place
  to prove the machinery — but **no longer the largest fraction, and its prerequisites are heavier
  than the first draft assumed.**
  > **⚠ GATE RE-EVALUATED — the first draft declared this PASSED on a number that does not survive
  > correct configuration.** It read loop A at 54.4% on Lane D. That figure is the sum of *two*
  > `elem.update` sites — the double constitutive pass — and L3-0b(i) deletes the second one. Measured
  > with `-commitSolveState`: **loop A = 38.95%, which FAILS the >40% gate**, while **loop C rises to
  > 46.33%** and becomes Lane D's dominant loop (`RESULTS` F9). Lane A's 80.79% still passes but is
  > vanilla `ForceBeamColumn2d` (§5.2). **So: no lane currently passes the >40% gate for loop A on a
  > correctly-configured, fork-owned deck.** L3-1 is therefore *not* authorized as the first threading
  > stage; it is authorized only as the machinery-proving stage, and its payoff must be re-argued
  > against the gate rather than assumed.

  **Hard prerequisites (all were previously later stages or absent — §5.4):** de-static
  `FE_Element`/`DOF_Group` pools (P-a); de-static the element kernel — **L3-2 must precede L3-1, not
  follow it** (P-b); an index-addressable element accessor, since the loop is driven by a shared
  mutable cursor and cannot be a `#pragma omp for` (P-c); eager `createDisp/createVel/createAccel`
  at `domainChanged` (P-d); per-thread profiler scopes (P-e). Plus `ops_TheActiveElement` →
  `thread_local` (§5.3, and note §5.4-H7: it is a per-call read for two materials, not a latch), an
  integer reduction on `ok`, the per-classTag opt-in allowlist (default empty), and
  `find_package(OpenMP)` + a build flag.
  **Still open:** barrier cost (§11 q3 — measured ~19.5 µs at 4T against 3.3–11.5 ms of work per
  region, so ~0.2–0.6%, i.e. not the binding constraint).
  **Acceptance: bit-identical to serial at 1/2/4/8 threads.** Anything less means a re-entrancy miss,
  not a rounding difference — this is the stage where that gate is free, so it must be enforced
  absolutely.

- **L3-2 — de-static the allowlisted fork-owned kernels into per-element buffers. ⚠ MOVED AHEAD OF L3-1** (§5.4-H2: `LadrunoBrick::update()` — L3-1's own named target — runs on 16 shared function-scope statics, and 63 of 215 `update()` bodies do).
  Also covers P-a (`FE_Element`/`DOF_Group` pools), which is ADR-40's rank 7 and is a **100% collision** in loops B/C independent of any element allowlist.
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
   > **⚠ Corrected 2026-07-25 (P1f).** As first written this was **one run per thread count** — the
   > exact experimental design that produced ADR-75 P1's wrong "PARDISO is bit-identical at every
   > thread count" claim. With a 50/50 split, a single-sample check **passes half the time**. So:
   > - run **N ≥ 10 repeats per thread count** and require *all* identical — a nondeterminism check
   >   is a claim about a distribution, not a sample;
   > - **also require each configuration to reproduce ITSELF** (10 runs, one binary, one thread
   >   count), not merely to match the serial baseline. P1f found the bug precisely by asking that
   >   question after an A/B reported "ux DIFFERS": the *old* binary showed the identical 5/5 split,
   >   which proved the variation was MKL's and the change was exact;
   > - **pin `MKL_NUM_THREADS=1`** for the whole gate (§3 P-6), or the solver's own ~1 ULP jitter
   >   masquerades as an assembly race;
   > - treat a *failure* as "which layer?" before "what did I break?" — the answer is a two-line
   >   experiment, and getting it backwards costs a day of debugging a correct change.

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
2b. **Elements that write shared node trial state inside loop A** (§2.1b) — **exactly two**,
   `LadrunoRigidBody` and `ZeroLengthVG_HG` (the other 11 grep hits resolve to loop C / `commitState` /
   `setDomain` / non-`Element` classes; "13" is the grep-hit count, not the exclusion list). This is an **ordering** race, so unlike everything else in this register it cannot be fixed
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
- `LEDGER_quirks.md` — **six rows added in this PR**: the `static`-return-buffer non-re-entrancy
  idiom; `ops_TheActiveElement`; the named-deep-scope-includes-the-scatter trap; the already-frozen CSR
  sparsity; ADR-40b's stale lane-D `formTangent` number; and `linearSolve` not being the solver cost
  under `DisplacementControl` (plus PARDISO having no `soe.factor` scope to reveal it).
- `banner_features.txt` — **not now.** Only once a threading feature is user-visible and `shipped`.

When stages land:
- `LEDGER_vanilla_files.md` — `Domain.cpp` (loop A), `Element.cpp` + **`OPS_Globals.h`** (and
  `G3Globals.h`; `ops_TheActiveElement` → `thread_local`), `IncrementalIntegrator.cpp` **and
  `TransientIntegrator.cpp`** (loops B/C/D/E — §2), **`FE_Element.{h,cpp}`, `DOF_Group.{h,cpp}`,
  `TransformationFE.cpp`** (the pools that actually race, §5.4-H1), **`AnalysisModel.cpp` +
  `SRC/tagged/storage/`** (index-addressable iteration, §5.4-H3), **`Node.cpp`** (eager trial-array
  creation, §5.4-H4), **`SRC/coordTransformation/`** (§5.4-H5), **`SRC/utility/profiler/`**
  (per-thread scopes, §5.4-H6), `CMakeLists.txt` (`find_package(OpenMP)`), each with a `// Ladruno`
  marker. **The list roughly tripled after the max-effort review** — that growth is itself the
  strongest single argument that Lane 3 is the highest-risk lane in the family.
- `LEDGER_implementations.md` — one row per threading feature + its flag token.
- `LEDGER_quirks.md` — this session's entries are due now (§0's stale-premise trap, the
  `static`-return-buffer idiom, `ops_TheActiveElement`, the frozen-sparsity fact).
- `Ladruno_scripts/banner_features.txt` + `patch_banner.py` — only once a threading feature is
  user-visible and `shipped`. **Not now.**

## 11. Open questions
1. **Does the §4.2 gather memory (Σ idSize² doubles) fit at production solid scale?** ~15.5 MB at Lane
   B, **~1.5 GB** projected at ~1M DOF. **This single number decides whether the byte-identical class-B CI
   gate of §3 survives** — the highest-leverage unknown in this ADR. Gate at L3-4; measure on a real
   deck.
2. **Does atomic `A[k] +=` actually contend?** If yes, §4's fallback ladder (LS-DYNA element ordering,
   then coloring) activates. Cheap to measure once L3-4 exists; do not pre-optimize.
3. **What is the real fork-join/barrier cost on this toolchain**, in µs? Sets the per-element-cost
   threshold that decides regression risk in every stage's gate. Measure once, reuse.
4. **How much of the ~34% Amdahl residual is loop A vs loops B/C vs neither?** L3-0 answers this; if
   the residual is mostly *neither*, Lane 3's whole ceiling is lower than assumed and the staging
   should shrink.
5. ~~**Replace the `addA` O(idSize² × rowlen) inner search.**~~ ✅ **CLOSED — SHIPPED as ADR-75 P1f
   ([#636](https://github.com/nmorabowen/OpenSees/pull/636), `580cfb5`), one day after L3-0 named it.**
   Linear scan → **binary search** (`ops_pardiso_findCol`), legal because `setSize` *enforces* the
   ascending-CSR invariant rather than assuming it. Measured interleaved A/B at 4 threads:
   **1.098× at 26.5k DOF, 1.091× at 50.7k** (`-matrixType 1`: 1.03×, less because half-storage already
   scatters fewer entries into shorter rows). **This validates L3-0's central argument** — the
   `addA` search was worth taking *before* any threading, and it needed no re-entrancy audit, no
   determinism policy and no OpenMP build. **The `A[k] += m(i,j)` race statement in §2.2/§4.1 is
   unaffected**: P1f changed only how `k` is *found*, not the accumulate. *(Original entry:)* L3-0
   measured it at 1 699.2 ms of an ~11.9 s step (~14% of wall), 1.28× slower than UmfPack's. L3-0 measured `PARDISOGenLinSOE::addA` at **1 699.2 ms of an ~11.9 s step
   (~14% of wall) on the default `-matrixType 0` path**, **1.28× slower than UmfPack's** 1 322.8 ms.
   `UmfpackGenLinSOE::addA` has the *identical* shape on frozen CSC (`for k in Ap[col]..Ap[col+1]: if
   (Ai[k]==row) { Ax[k] += …; break; }`), so **both desktop solvers pay it** — this is not a PARDISO
   patch. (`DiagonalSOE::addA` is `A[pos] += m(i,i)` into a dense diagonal: no search, nothing to
   fix.) That is ~14% of wall spent locating entries in a structure unchanged since `setSize`. Serial,
   zero determinism cost, no threading. **On present evidence a better next move than threading loops
   B/C at all** — and it is Lane-1 work, not Lane-3.
6. **Cluster composition:** loop-A threading is the "threads-per-node half" of hybrid MPI+threads
   (ADR-75 §3). Untested interaction with `PartitionedDomain`; defer until L3-1 lands serially.
7. **Does any element alias (rather than `getCopy()`) a material or section across elements?** Not
   verified (§2.1b item 4). If one does, threaded loop A corrupts shared material state. An L3-2
   verification item; ThreadSanitizer over the allowlisted set is the practical check.
8. ~~**Is Lane 3 still worth doing at all once the cheap work-removal items land?**~~
   ✅ **ANSWERED 2026-07-25 — see §12. Verdict: PARK, do not start; re-entry gated on G-L3.**
   *(Original: L3-0b measures −29.6% on Lane D from one shipped flag, and the `addA` search is ~14%
   of wall on both solvers. Neither needs a re-entrancy audit, a determinism policy, or an OpenMP
   build. **This ADR's own data is the strongest argument for doing them first** — and the honest
   question is whether what remains after them still clears the gate. Re-measure the gate *after*
   L3-0b, not before.)*
9. **What does per-element de-statication cost in SERIAL memory?** L3-2 mandates per-element buffers;
   `LadrunoBrick`'s function-scope statics are ~3.5 KB/ele and its class-level `stiff/resid/mass/xl`
   ~9.6 KB/ele ⇒ order **1–4 GB at the 325k-element scale** where §4.2 already flags 1.5 GB as
   possibly unaffordable — paid in **every** run, including serial. L3-2's gate ("serial performance
   and answers unchanged") does not currently cover memory. Fix the gate.
10. **Does the F1 gate verdict survive at production size?** Lane B here is 11 520 DOF, and P1c
   measured the solve fraction *growing* with N. At production scale the element fraction is lower
   than the measured 74.85%, so the margin shrinks. Re-check before committing L3-3/L3-4. Not in
   doubt for lanes A and D (81–87%, no solver in the picture).

---

## 12. VERDICT on q8 — **PARK Lane 3. Do not start. Re-entry is gated on G-L3.**

*(2026-07-25. Answers §11 q8, which the handoff names as the question to settle before any Lane-3
planning. This is a decision about **sequencing**, not a claim that threaded assembly is unsound.)*

### 12.1 The finding: the reason this lane was attractive no longer describes any ranked target

ADR-75b's structural contribution is §2.1 — *loop A carries no floating-point reduction, so threading
it is bit-identical and the fork's existing byte-identical oracles gate it unchanged, with no new CI
mode and no ordered-mode cost.* That is what made Lane 3 look like it had a cheap entry point.

After L3-0b landed, `RESULTS_l3a_update_scope.md` §3's revised ranking reads:

| rank | target | fraction | reduction? | determinism policy needed? |
|---|---|---|---|---|
| 1 | Lane D **loop C** (`formUnbalance`) | 46.33% | **yes** (`addB`) | **yes — full §3 ordered/fast** |
| 2 | Lane B **loop B** (`formTangent`) | 61.56% | **yes** (`addA`) | **yes**, + the §4.2 gather-memory question |
| 3 | Lane D / B **loop A** | 38.95% / 5.24% | no | no |

**Both targets now ranked above loop A are reducing loops.** The reduction-free case fell to rank 3
and, on Lane D, to *below the ADR's own >40% gate* (38.95%). So the remaining work is precisely the
expensive version of this lane: every ranked target needs the §3 ordered/fast policy, an OpenMP
build, and P-a…P-e in full. **The cheap entry point was real, and it is gone.** Starting now means
paying the whole prerequisite bill before the first measurable win — the opposite of how L3-0b
itself paid off.

### 12.2 Four costs that are now priced, and were not when the lane was opened

1. **Five mandatory prerequisites, none of which pay anything on their own** (P-a…P-e). P-a is not
   negotiable by scoping: `FE_Element::theMatrices` is a **class-wide pool** handed straight to
   `addA` — a 100% collision that no element allowlist can fix.
2. **A serial memory regression is now part of the price.** §11 q9: per-element de-statication is
   order **1–4 GB at 325k elements**, paid in *every* run **including serial**, and L3-2's gate does
   not yet even cover memory. The lane currently proposes spending GB of serial memory to buy a
   **≤2.78× (4T Amdahl)** threaded ceiling on Lane D.
3. **The byte-identical CI gate may not survive its own remedy** — §11 q1 puts the §4.2 gather at
   ~1.5 GB at 1M DOF, i.e. the ordered mode that makes threading exact may be unaffordable exactly
   at production scale.
4. **The loop is not index-addressable** (P-c): a shared mutable cursor, so this is not a
   `#pragma omp for` retrofit but a change to how the domain hands out elements.

### 12.3 The decisive gap: the gate has never been measured where the fork intends to run

Every fraction authorizing Lane 3 comes from **Lane B at 11 520 DOF** or a single explicit deck.
ADR-75's standing correction (§1) is that the production regime is **huge solid nonlinear models**,
that the solve fraction is a **floor that grows with N**, and that *no gate measured at ≤136k DOF
should be treated as production-representative*. §11 q10 states the consequence plainly: at
production scale the element fraction is **lower** than the measured 74.85%.

So the single number that authorizes this lane — "PARDISO flips Lane B's gate, 35.8% → 74.9%" — is
measured at a size **~12× below the 136k DOF that ADR-75 §1 already declares non-representative**,
~90× below 1M DOF, and ~1 600× below the 18.6M-DOF deck ADR-74 treats as the production kill-run —
and it moves in the direction known to shrink it. **Committing prerequisite work against an
unmeasured gate is the confounded-comparison failure mode this program has already paid for three
times.**

### 12.4 What to do instead, and the one measurement that reopens the lane

**G-L3 (the re-entry gate).** On a **production-scale** deck (≥500k DOF, the real solid nonlinear
regime — `system Pardiso -matrixType 2` on the desktop or `system Mumps` on the cluster), profile
with **coarse** instrumentation (`OPS_PROF_FLAGS="-perStep"`, no `-deep`) and report **the fraction
of step held by each individual assembly loop** — loop A (`Domain::update`), loop B (`formTangent`),
loop C (`formUnbalance`) — plus the solve %.

- **The largest SINGLE loop ≥ 40% of step** ⇒ re-authorize Lane 3 *at that loop*, with q1/q9
  answered **before** any code.
- **No single loop reaches 40%** ⇒ **close** Lane 3 as Amdahl-irrelevant in the production regime;
  the residual belongs to Lane 1/2 (solver) work, not assembly.

> **⚠ Read the gate quantity carefully — an earlier draft of this section got it wrong, and the
> error reversed the answer.** The >40% gate is evaluated on a **single loop's** fraction, **not** on
> the aggregate element-kernel fraction. §3's Lane D row is the proof: **kernel 85.30% but
> gate = FAIL**, because the loop being considered (loop A) was 38.95%. A gate written on the
> aggregate would have passed Lane D at 85.30% and re-authorized exactly the stage F9 de-authorized.
> Aggregate kernel % tells you the lane's *ceiling*; a single loop's % tells you what one threading
> stage can actually buy.

**Why coarse, and what coarse cannot give you.** `enabled()` gates coarse phase timing while
`deep()` *additionally* gates the per-element buckets (`Profiler.h:327-334`), so `-perStep` alone
does yield the per-loop phase fractions the gate needs — untaxed. What it **cannot** give is the
kernel-vs-scatter split *within* a loop: that is `scatter = scope_wall − Σ(elem_by_type wall)`, and
`elemByType_` is *"lazy; deep mode only"* (`Profiler.h:149`). That split is a **separate, later**
question — it decides how much of a passing loop is actually threadable (scatter is not), and on
Lane A it was decisive (scatter 44.1 ms **exceeded** kernel 35.4 ms). Take it on a **reduced-size**
deck: at production scale `-deep` is likely prohibitive *and* distortive, since §5 prices it at
0.14–0.39 µs per scope instance and instance count scales with elements × steps (Lane D already
carried 15.25 M instances at only 1 000 elements).

G-L3 costs **one profiled run and no code** — the same shape and cost class as the ADR-75 handoff's
item-1 cluster run, and it should be captured in the *same* run. That is the whole recommendation:
**buy the number before buying the lane.**

**Also worth preferring over threading, on this lane's own evidence:** the −29.6% that de-authorized
L3-1 came from a **shipped flag the deck never set**, found by accident during a review. Nothing in
Lane 3 has yet beaten a deck-configuration audit on return-per-unit-risk. *(Checked, so it is not
offered loosely: ADR-67's sibling P-NEW-1 constant-mass tangent cache is **already default-on**
(`CentralDifferenceLadruno.cpp:89`), so that particular well is dry — the point is the class of win,
not another instance of it.)*

### 12.5 What this verdict does **not** say

- Not "threaded assembly is wrong". §3 and §4 are settled and stay settled; §2's five-loop taxonomy
  is the reusable asset and remains correct.
- Not "the fractions are wrong". Lane D still passes **in aggregate** (85.30% kernel). The objection
  is to starting an expensive lane against a gate measured off-regime.
- Not a claim about the cluster. Hybrid MPI+threads (§11 q6) is untouched and still deferred.

**Status change:** `proposed` → **`parked — measurement-gated (G-L3)`**. No threading code
authorized; no prerequisite work authorized either, which is the part that changed.
