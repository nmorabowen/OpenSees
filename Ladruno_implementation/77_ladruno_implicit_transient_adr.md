---
title: Implicit transient solver study — step anatomy, reuse levers, and the fork-integrator question
project: Ladruno
status: active — C0 + T0 + T0b + T1 + T3 measured; G2 and G3 CLOSED; C0-6 bug FIXED + verified; T2-nodal-mass + rest of C0 patch wave open
priority: medium
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - adr
  - performance
  - transient
  - integrator
  - solver
  - newton
  - sub-adr
---

# ADR-77 — Implicit transient solvers: where does a dynamic step actually spend its time, and do we need our own integrators?

> ADR-77. A **solve-lane** perf sub-ADR under [[40_ladruno_performance_adr]]. It is the
> transient-path sequel to three closed ADRs and inherits their conclusions instead of
> re-deriving them:
>
> - [[49_ladruno_integrator_study_workflow_adr]] / [[49a_integrator_scorecard_2026-06-23]] —
>   graded every integrator; the implicit transient family (Newmark, HHT, HHTGeneralized,
>   GeneralizedAlpha, TRBDF2/3, BackwardEuler, plus legacy) is **mature and complete**, the
>   only named scheme gap being the modern implicit Bathe β1/β2 (Rank 3, "assess then maybe").
> - [[75_ladruno_sparse_direct_strategy_adr]] — PARDISO shipped; the factorization is no
>   longer the dominant per-step cost. [[75b_ladruno_threaded_assembly_adr]] L3-0 measured
>   the consequence: under PARDISO the **element loop rises to ~75%** of the step (Lane B,
>   35.8% → 74.9%).
> - [[76_ladruno_tangent_reuse_adr]] — the "this matrix cannot have changed" premise is
>   treacherous (`getInitialTangent()` defaults to `getTangent()`; invariance predicates must
>   default false). Its §4.5 remainder, `LadrunoModifiedNewton`, is designed and unstarted.

## 1. The two questions this ADR exists to answer

**Q1 (the decision the owner asked for): should the fork commit time to authoring its own
set of implicit transient integrators** (a `LadrunoNewmark`-style family mirroring what
ADR-58/75 did for elements and solvers)? Q1 has **two distinct motives that must not be
conflated**, because they have different answers:

- **Q1a — performance:** would a fork integrator be faster? (§4.1: no — every lever lives
  in a shared layer.)
- **Q1b — ownership for fixes:** upstream's implicit family is *inconsistent* (§4.2 has
  the evidence). Do we need to own the class to apply improvements or fixes, or can we
  patch the vanilla files in place?

**Q2 (the study): where does one implicit dynamic step spend its time under the fork's
current best stack** (`system Pardiso` + `algorithm Newton` + Newmark/HHT), **and which
reuse levers are worth building?**

Q2 produces the evidence; Q1 is decided by gates over that evidence (§6). The position
taken now, subject to those gates: **No — not as a family**, for either motive. §4 argues
it; §6 defines what would overturn it. On Q1b specifically: the fork already patches
vanilla files routinely (that is what `LEDGER_vanilla_files.md` exists for), so ownership
is **not a precondition for fixing** — it is warranted only under the narrow §4.2 ownership rule.

## 2. What is already answered (do not re-study)

| Question | Answered by | Verdict on record |
|---|---|---|
| Which implicit schemes exist / are they any good? | ADR-49a §implicit table | All "keep"; GeneralizedAlpha is the gold standard; consolidations C1/C6 identified, low priority |
| Do we have DDM sensitivity in transient? | ADR-52 W3-I2 | Shipped: `LadrunoHHT`, `LadrunoGeneralizedAlpha` (fork classes exist **for sensitivity**, not performance) |
| Is the factorization the bottleneck? | ADR-75 P1/P1d/P1f | No longer — PARDISO ~O(N^1.45), symmetric −42% memory; `addA` search fixed both solvers |
| Is threading assembly the next move? | ADR-75b L3-0/L3-0b | **Work removal outranks threading**; no loop currently passes the >40% threading gate on a correctly-configured deck |
| Can we skip re-forming an unchanged tangent? | ADR-76 | Engine-side invariance tracking WITHDRAWN; user-side spelling (`ModifiedNewton -initial [-factoronce]`) is the shipped answer; `LadrunoModifiedNewton` designed, unstarted |

What none of them measured: **the transient step itself.** Every ADR-75 number is a
statics/pushover-shaped lane. The P1f correction explicitly warned the transient path is
different — it assembles through **five loops, not three** (`TransientIntegrator.cpp` element
tangent + DOF_Group tangent, plus the `addB` sides), "so a gather keyed on element matrices
would silently drop them" ([[_adr75_session_handoff]] §next-steps 3).

## 3. Anatomy of one implicit transient step (the facts the study starts from)

Per Newton iteration, `TransientIntegrator::formTangent`
(`SRC/analysis/integrator/TransientIntegrator.cpp:70`) does:

1. `zeroA()`;
2. modal damping matrix if present;
3. **DOF_Group loop** — nodal mass/damping into `addA` (`:102`), *not element-keyed, currently
   untimed by the P3 deep profiler* (the hook at `:112` covers only the element loop);
4. **FE_Element loop** — `formEleTangent` per element.

`Newmark::formEleTangent` (`SRC/analysis/integrator/Newmark.cpp:284`) zeroes the element
tangent and re-adds `c1·K + c2·C + c3·M` **every iteration of every step**. And
`FE_Element::addMtoTang` / `addCtoTang` (`SRC/analysis/fe_ele/FE_Element.cpp:383,400`) call
`myEle->getMass()` / `getDamp()` **fresh each time — there is no caching layer anywhere in
this path.** For most elements `getMass()` recomputes a density integral (or at best copies
a stored matrix) and `getDamp()` recombines Rayleigh terms per call.

Known constants being recomputed:

- **M** is constant for essentially every fork model (no adaptive mass, no element
  birth/death mid-analysis in the current decks). Caveats to check, ADR-76-style, before
  claiming invariance: `LadrunoMassScaling*` (explicit-only today), staged constructors
  (`StagedNewmark` diag fill), contact SOFT penalties that inject effective mass.
- **C** under Rayleigh with fixed coefficients is constant **only if `betaK`/`betaKcommit`
  terms are zero or K is elastic** — ADR-76 §5.1 already flagged `betaK≠0` as exactly the
  trap that makes "constant" false. Any C-cache must key on the Rayleigh recipe.
- **K_eff structure**: sparsity is already frozen (`zeroA` never touches `colA`/`rowStartA`,
  ADR-75b §4) and the symbolic factorization is already reused (P1a). The numeric
  recombination `c1·K + c2·C + c3·M` is the part done from scratch.

Plus the L3-0 baseline transplanted from statics: with PARDISO the element loop is ~75% of a
solve-bound lane. If the transient step behaves the same way, **the study's target is the
assembly side, not the scheme and not the factorization** — but that transplant is exactly
what T0 must confirm or kill, because the transient step adds loops statics doesn't have.

> [!warning] **MEASURED 2026-07-26 (T0 + T0b) — true, but only under a threaded solver.**
> At **1 thread** the transplant fails: `elem.tangent` 30.4% → 40.8% under PARDISO while
> `linearSolve` goes 58.5% → **45.2%**, leaving the solve the largest single loop. The
> thread sweep then restored the premise: PARDISO crosses over at **2 threads**, and total
> assembly reaches **59.6% (2 thr) / 64.7% (4 thr) / 73.9% (16 thr)** — the last essentially
> L3-0's 74.9%. **UmfPack never crosses over at any thread count**, so the correct claim is
> *"a fast threaded solver makes assembly dominant"*, not *"transient steps are
> assembly-bound"*. G0 earned its keep twice: the answer is both thread- and
> solver-dependent, which the transplant would have hidden in either direction. Operating
> point: `MKL_NUM_THREADS=4` (wall degrades past 4 on this hybrid P/E part). Full numbers,
> the G2 kill, and the C0 audit: [[77a_c0_t0_results_2026-07-26]].

## 4. Q1 — the fork-integrator question, argued

### 4.1 What a fork implicit integrator family could plausibly buy

| Candidate motive | Where the lever actually lives | Needs a new integrator? |
|---|---|---|
| Skip re-forming constant M/C contributions | `FE_Element`/`DOF_Group` tangent path + integrator `formEleTangent` protocol — **cross-cutting, benefits every integrator at once** | **No** — a `LadrunoNewmark` would fork the wrong layer and orphan HHT/GenAlpha users |
| Reform/refactor less often | Algorithm layer — `LadrunoModifiedNewton` (ADR-76 §4.5, designed) | **No** — it is an `EquiSolnAlgo` |
| Factor once in linear transient | `Linear -factorOnce` exists; its quirks are ADR-76 §5 territory | **No** |
| Numerical damping to buy fewer iterations | Already have HHT/GeneralizedAlpha with full ρ∞ control | **No** — parameter advice, not code |
| DDM sensitivity | Already shipped (`LadrunoHHT`, `LadrunoGeneralizedAlpha`) | Done |
| Modern implicit Bathe β1/β2 (2012) | A genuinely missing *scheme* (ADR-49a Rank 3) — value claim is robustness on contact/sharp nonlinearity, **not speed** | **Yes, but only this one**, and ADR-49a already gated it: benchmark TRBDF2 first |
| Per-element-type timing, stats verbs | P3 profiler already instruments the vanilla classes | No |

The pattern is uniform: **every performance lever lands in a layer that is shared across
integrators.** The integration *schemes* themselves are commodity mathematics that upstream
implements correctly and ADR-49a already graded "keep — mature." Forking them would buy
maintenance surface (every upstream fix to Newmark/HHT now needs porting), a second place
for the sensitivity code to rot, and zero measured speed — the exact anti-pattern ADR-75a's
portfolio-vs-unify trade study warned against on the solver side.

### 4.2 The ownership axis — is upstream inconsistent, and do fixes require owning the class?

**Yes, it is inconsistent — measurably.** A 10-minute source survey (2026-07-26, this
worktree) across the implicit family:

| Class | `INITIAL_TANGENT` | `HALL_TANGENT` | DDM sensitivity | `sendSelf` | notes |
|---|---|---|---|---|---|
| `Newmark` | ✔ | ✔ | ✔ (disp-form only) | sends γ,β,displ | also carries a `converged` flag the others lack |
| `HHT`, `HHTGeneralized`, `GeneralizedAlpha`, `TRBDF2/3`, `WilsonTheta`, `Houbolt` | ✔ | ✔ | ✘ | varies | — |
| `BackwardEuler`, `Collocation`, `Newmark1` | ✔ | **✘** | ✘ | varies | `-hall` silently unsupported |
| `TRBDF2` | ✔ | ✔ | ✘ | **stub `return 0`** | benign today (parameterless) but the *pattern* is the trap the contact handler's P1a stubs already taught us |

Add the structural inconsistencies ADR-49a already booked: `Newmark1` applies Rayleigh
damping at `newStep` (semantics differ from `Newmark`), `TRBDF3` composites with Houbolt
instead of Bathe's third sub-step, `Houbolt`'s brittle step-count bootstrap, the ~24
copy-paste `*HS*`/`_TP` hybrid variants, and DDM sensitivity existing in exactly one class
(now three, counting the fork's `LadrunoHHT`/`LadrunoGeneralizedAlpha`).

**But inconsistency ≠ ownership requirement.** The fork has three delivery vehicles for a
fix, and the precedent record shows the cheap one usually wins:

| Vehicle | Precedents | When |
|---|---|---|
| **Patch vanilla in place** (`// Ladruno` comment + `LEDGER_vanilla_files.md` row + candidate for the upstream-PR campaign) | `BandGenLinSOE::setSize` quicksort ([#593](https://github.com/nmorabowen/OpenSees/pull/593)); `OPS_ModifiedNewton` parser loop (ADR-76 R4); `addA` binary search (P1f/P1g); LAPACK singular-as-success (ADR-76 §6) | Outright bugs and consistency gaps where vanilla behavior is indefensible and the fix preserves results on already-correct inputs |
| **Fork-owned class** (new classTag, ledger row, banner line) | `LadrunoHHT`/`LadrunoGeneralizedAlpha` (DDM); `ExplicitBathe` family; `CentralDifferenceLadruno` | See rule below |
| **Shared-layer feature** (algorithm / FE_Element / SOE) | PARDISO `-krylov`; P3 profiler; `LadrunoModifiedNewton` (planned) | Anything meant to benefit every integrator at once |

**The ownership rule.** Own the class only when at least one holds:

1. **The fix changes numerical results** of an upstream verb on decks that currently run —
   then vanilla must stay intact side-by-side, both as the A/B oracle and because the
   fork's QA is bit-exact (G4). This is why `LadrunoHHT` exists rather than a patched
   `HHT`: adding state to the update path risks perturbing the vanilla trajectory, and
   sensitivity semantics are a behavior change.
2. **The fix is a restructure** (e.g., collapsing the `*HS*`/`_TP` copy-paste explosion,
   or replacing `TRBDF3`'s Houbolt leg with Bathe-3) too invasive to carry as a
   greppable vanilla diff across upstream syncs.
3. **Upstream would reject or indefinitely stall it** *and* the diff is large or
   recurring-conflict-prone. (Small rejected diffs just stay as ledgered vanilla edits —
   that is the fork's normal operating mode.)

Everything else — missing `HALL_TANGENT` branches, `sendSelf` stubs, guard clauses,
parser loops — is a **patch-in-place**, one ledger row each, and feeds the upstream-PR
campaign ([[upstream_pr_campaign]]) as a "consistency wave." Owning eleven integrator
classes to fix twenty small cells of that matrix would convert a page of greppable diffs
into a permanent maintenance surface — and decouple us from every future upstream fix to
those classes.

### 4.3 Position

**Do not author a fork implicit-integrator family — for either motive.** Performance
levers live in shared layers (§4.1); consistency fixes are deliverable as ledgered
vanilla patches (§4.2). Commit the time to, in order:

1. **T0–T2 measurements + the C0 consistency audit** (§5) — cheap, and every later
   decision keys on them;
2. **The C0 patch-in-place wave** — the `HALL_TANGENT`/`sendSelf`/guard-clause cells of
   the §4.2 matrix that C0 confirms, one `LEDGER_vanilla_files.md` row each, batched as
   a consistency wave for the upstream-PR campaign;
3. **`LadrunoModifiedNewton`** (ADR-76 §4.5) — already designed, ~small, algorithm-layer,
   and it is the delivery vehicle for "reform less often" regardless of integrator;
4. **A constant-M/C tangent-contribution cache** *if and only if* T2 shows the headroom
   (gate G2) — designed at the `FE_Element`/integrator-interface layer with an opt-in flag
   and invariance defaulting to **false** (the one surviving design rule of ADR-76 App. A.4);
5. **Implicit Bathe β1/β2** *if and only if* the T3 benchmark shows TRBDF2/HHT actually
   losing on cost-per-simulated-second or robustness on a fork-real contact deck (gate G3).
   That is one class, not a family, and it is a robustness purchase, not a speed one.

## 5. Q2 — the measurement plan

**Environment (fixed by owner, 2026-07-26): this PC** — serial desktop, same sequencing
discipline as ADR-75. Standard stack for all runs unless a stage varies it deliberately:
`system Pardiso` (now available; 1 thread on oracle runs per G4) + `numberer
LadrunoParallel` (ADR-30 — the T0 stall-hunt lesson: never bench on the default RCM at
scale) + P3 profiler + `-stats`. Vehicles: one synthetic lane (Lane-B-style brick block,
scalable N) + one fork-real deck; **complex meshes come from apeGmsh** (via `opensees_env`
py3.12 + the `APEGMSH_OPENSEES_BIN` seam, harness pattern `run_ladruno_integration.ps1`) —
so the fork-real deck need not be limited to what hand-written Tcl can mesh. All timing
runs A/B interleaved, 3 rounds, per the ADR-75b method rules.

- **C0 — consistency audit (no timing, pure source/behavior survey).** Extend the §4.2
  matrix to every cell: `statusFlag` support, sensitivity hooks, `revertToLastStep`
  semantics, `sendSelf`/`recvSelf` completeness, Δt-change handling in `newStep`,
  `domainChanged` state preservation, Rayleigh-damping placement. Each defective cell
  classified **patch-in-place vs needs-ownership** by the §4.2 rule. Deliverable: the
  matrix + the patch-wave PR list. Can run concurrently with T0.
- **T0 — step anatomy.** Newmark + Newton, UmfPack vs PARDISO. Deliver the **five-loop
  split** of one converged transient step: element tangent, DOF_Group tangent, `addA`
  scatter, factor, back-solve, `formUnbalance` (element + DOF_Group `addB`), `update`,
  `commit`. Requires one small instrumentation addition: a deep-profile scope on the
  DOF_Group loop at `TransientIntegrator.cpp:102` (currently untimed, noted in the P3
  comment at `:110`). **Everything downstream keys on this table.**
- **T1 — reuse levers as shipped today.** ✅ **MEASURED — see [[77a_c0_t0_results_2026-07-26]] §6c.**
  Verdict: the winner is **nonlinearity-dependent** (Newton at 1% past yield;
  **KrylovNewton at 2.4x and 7x yield, 1.75x faster than Newton at ×4 and the only
  practical arm at ×10**), every lever loses to plain Newton on a near-elastic deck, and
  the modified/initial-stiffness family is a **robustness** lever (converges at ×10 where
  Newton diverges) rather than a speed one. Same deck: `Newton` vs `ModifiedNewton` vs
  `ModifiedNewton -factoronce` vs `KrylovNewton` (P1c `-krylov` on the PARDISO side), each
  at equal convergence tolerance. Deliver cost-per-converged-step and iterations/step.
  Under PARDISO the textbook ranking may invert (factor is cheap, assembly isn't) —
  measure, don't assume; that is the L3-0b lesson.
- **T2 — constant-M/C headroom (upper bound before any design).** Time spent inside
  `addMtoTang`/`addCtoTang` + `getMass`/`getDamp` + the DOF_Group tangent loop per step.
  That number **is** the ceiling of the cache win of §4.3(4). If it is small, the cache
  dies here at zero design cost.
- **T3 — scheme shootout on cost-per-simulated-second.** ✅ **MEASURED — §6d.**
  **G3 NOT AUTHORIZED** (TRBDF2 1.18-1.25x, best damped 1.00-1.07x, bar is ≥1.3x) ⇒ the
  implicit Bathe β1/β2 lane is **closed with the benchmark on record**. Also **falsified**
  this stage's own hypothesis — numerical damping bought *zero* iteration reduction. And it
  turned up **C0-6**, a confirmed bug in `GeneralizedAlpha::update()` that
  `LadrunoGeneralizedAlpha` inherits. Original scope text: Newmark(γ=½,β=¼) vs HHT(α) vs
  GeneralizedAlpha(ρ∞) vs TRBDF2 on a *nonlinear* fork-real deck: iterations/step ×
  cost/step × achievable Δt at matched accuracy (energy/drift oracle). Feeds gate G3.
  This is the ADR-49a Rank-3 "assess" finally executed.
- **T4 — parallel implicit transient (MUMPS/SP/MP).** Explicitly deferred; not in this
  ADR's first wave. Opened only after T0–T3 close, as its own brief.

## 6. Gates (decision rules, in the fork's numeric style)

- **G0 (study validity):** T0 must be run on both a synthetic and a fork-real deck. No
  conclusion transplants from statics (L3-0) or from the synthetic lane alone.
- **G1 (fork integrator, general):** a new fork integrator class may be authored **only if**
  a measured lever ≥ 15% of step wall on a fork-real deck **and** the lever provably cannot
  be implemented at the algorithm or FE_Element layer. (Both conditions; §4.1's table
  predicts the second is never met for performance motives.)
- **G2 (constant-M/C cache):** design authorized only if T2's headroom ≥ 10% of step wall
  under PARDISO on a fork-real deck. Invariance opt-in, default false, keyed on the
  Rayleigh recipe (`betaK` ⇒ C is not constant, ever).
- **G3 (implicit Bathe β1/β2):** authorized only if T3 shows TRBDF2 *and* HHT/GenAlpha
  each losing ≥ 1.3× on cost-per-simulated-second at matched accuracy, or failing outright,
  on a contact/sharp-nonlinearity fork deck. Otherwise closed with the benchmark on record
  (which is itself the ADR-49a Rank-3 deliverable).
- **G4 (standing order, inherited from ADR-40 via 75b):** work removal before concurrency;
  measurement before design; A/B interleaved with bit-identical oracles where the
  determinism policy applies (PARDISO 1-thread for oracle runs — threaded PARDISO is not
  byte-reproducible, P1f).
- **G5 (ownership for fixes):** a consistency/bug fix ships as a **ledgered vanilla
  patch by default**. Fork-owning the class requires the §4.2 rule — results-changing,
  restructure-scale, or large-and-upstream-rejected — recorded in the C0 matrix cell
  before any `Ladruno*` integrator file is created. "Upstream is inconsistent" alone
  never satisfies G5.

## 7. Non-goals and open items

**Non-goals of ADR-77:** explicit dynamics (owned by ADRs 04/05/52); hybrid-simulation
`*HS*`/`_TP` variants (out of scope per ADR-49a); iterative/Krylov linear solvers for
transient (no credible preconditioner story in-tree; revisit only if T0 shows back-solve
dominating at scales beyond PARDISO's comfort); consolidations C1/C6 (upstream-hygiene,
belongs to the upstream-PR campaign, not here); PFEM.

**Resolved by owner (2026-07-26):** tests run on **this PC** (serial desktop confirmed,
T4 stays deferred); standard stack includes `numberer LadrunoParallel` and `system
Pardiso`; apeGmsh is the mesh source when the deck outgrows hand-written Tcl.

**Open items for the owner:**

1. Pick the fork-real T0/T3 deck: Esmeralda G3, the np8 rung model (has the
   profile-suffix + h5 rollup pattern from the stall hunt), or a purpose-built
   apeGmsh-generated nonlinear deck.
2. `LadrunoModifiedNewton`: fold into this study's wave 2 (recommended — T1 measures the
   exact niche it fills) or keep as a standalone ADR-76 remainder.

## 8. Status log

- 2026-07-26 — ADR drafted; scoping only, no measurements run. Position on Q1 recorded
  (§4.2): no fork integrator family; gated exceptions G2/G3 only. T0 instrumentation gap
  identified (DOF_Group tangent loop untimed).
- 2026-07-26 — **C0-6 FIXED, rebuilt, verified; G3 recomputed** (§6d.3 of the results doc).
  Fix applied to the **vanilla base class**, not `LadrunoGeneralizedAlpha` — the fork class
  inherits `update()`, so one line repairs both, and owning it would have duplicated ~40
  lines while leaving vanilla broken. All 8 acceptance checks pass: αM=1.0 **bit-identical**,
  every other αM's observed order −1.34/−4.41/−6.71 → **+1.30/+1.18/+1.19**. **G3 re-derived
  on the now-valid data: still NOT AUTHORIZED**, and more firmly — the repaired
  GeneralizedAlpha is the second-best scheme (1.05x). Noise floor recorded: wall-clock
  varies ~12% run-to-run, so cost ratios under ~1.2x are not a ranking.
- 2026-07-26 — **T3 measured** (§6d). **G3 NOT AUTHORIZED ⇒ implicit Bathe β1/β2 CLOSED**
  with the benchmark on record (the ADR-49a Rank-3 deliverable). Damping bought zero
  iteration reduction, falsifying §5's T3 hypothesis. **C0-6 found: a real bug in vanilla
  `GeneralizedAlpha::update()`** (αM-weighted acceleration computed then discarded),
  **inherited by `LadrunoGeneralizedAlpha`**; both correct only at αM=1.0. One-line
  patch-in-place, not yet applied. With G2 and G3 both closed, **every gated code item in
  this ADR is now decided, and none authorized a fork integrator** — Q1a and Q1b both land
  on "no" with measurement behind them.
- 2026-07-26 — **T1 measured** (§6c of the results doc). The inference T0b invited —
  "assembly is 65%, so assembly-skipping levers win" — is **false**: `ModifiedNewton` buys
  a 0.33x cheaper iteration for 3.88x more iterations (net 1.27x slower). Winner is
  nonlinearity-dependent; `KrylovNewton` is the standout. The ModifiedNewton family is
  re-cast as a **robustness** lever (converges at 7x yield where Newton diverges), which
  **re-motivates `LadrunoModifiedNewton` on robustness grounds and removes its speed
  justification**. New **C0-5**: the P3 profiler instruments only 3 of 11 solution
  algorithms (patch-in-place). ADR-76 §3's trap does **not** fire on fork materials.
- 2026-07-26 — **C0 + T0 + T0b measured** ([[77a_c0_t0_results_2026-07-26]]). §3's
  assembly premise holds only under a threaded solver (PARDISO crossover at 2 threads,
  assembly 64.7% at 4; UmfPack never crosses). **G2 fails its gate at every thread count
  (4.4-6.9% vs ≥10%) ⇒ constant-M/C cache not authorized.** C0: 4 confirmed findings, all
  patch-in-place, none requiring ownership. T1 is the live stage.
- 2026-07-26 (same day) — **Q1 split into Q1a (performance) / Q1b (ownership for fixes)**
  after owner raised upstream inconsistency. §4.2 added: measured inconsistency matrix
  (`-hall` missing from BackwardEuler/Collocation/Newmark1; sensitivity in Newmark only;
  `TRBDF2::sendSelf` stub; `Newmark1` Rayleigh-at-newStep; TRBDF3/Houbolt structure) +
  the ownership rule + G5. C0 consistency-audit stage added. Environment fixed: this PC,
  `numberer LadrunoParallel`, PARDISO, apeGmsh for meshes.
