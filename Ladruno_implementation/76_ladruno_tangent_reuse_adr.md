---
title: Tangent reuse — stop re-assembling and re-factorizing an unchanged matrix
project: Ladruno
status: R1/R4 shipped + adversarially reviewed; R2 WITHDRAWN as designed (see §8); R3 rejected as specified
priority: medium
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - adr
  - performance
  - solver
  - algorithm
  - newton
  - tangent-reuse
  - sub-adr
---

# ADR-76 — Tangent reuse: `algorithm Newton -initial` re-assembles and re-factorizes an invariant matrix every iteration

> ADR-76. An **assembly/solve-lane** perf sub-ADR under [[40_ladruno_performance_adr]].
> Complements [[75_ladruno_sparse_direct_strategy_adr]]: ADR-75 made the *factorization*
> cheaper (PARDISO) and added the *other* reuse axis (`-krylov`, reuse under a tangent that
> **has** changed); this ADR is about the axis ADR-75's `factored` flag was built for but
> that nothing upstream ever pulls — **not doing the work at all when the tangent has not
> changed**.
>
> Origin: an external issue report from the **TIMs project (SFIM continuum reference model)**
> against `bdf8adf9f`, with measured numbers on a 38 984-DOF `SSPbrickUP`/`ShellMITC4` model.
> Source file: `ISSUE — Newton-initial redundant refactorization.md`.

## 1. Context — the reported defect

`algorithm Newton -initial` performs a full tangent assembly **and** a full numeric
factorization on **every iteration** of every step, of a matrix that in most configurations
cannot have changed since the previous iteration.

`formTangent` sits inside the iteration `do`-loop in `NewtonRaphson::solveCurrentStep`. The
only tangent flag given special treatment is `INITIAL_THEN_CURRENT_TANGENT`, which correctly
forms the initial tangent on iteration 0 and the current tangent thereafter. `INITIAL_TANGENT`
falls through to the generic branch and is re-formed unconditionally. Contrast
`ModifiedNewton::solveCurrentStep`, which forms the tangent **once before** its loop.

> [!note] Citations in this ADR use function-name anchors, not line numbers, for
> `NewtonRaphson.cpp` and `ModifiedNewton.cpp`. The R1/R4 edits shift those files, and an
> earlier draft cited line spans that its own commit invalidated — including
> self-referentially, where the comment block displaced the very `do`-loop it pointed at.
> Line numbers are kept only for files this ADR does not touch.

Each re-formation costs, in order:

1. `theSOE->zeroA()` plus a full loop over every `FE_Element` calling `getTangent()` and
   `addA()` — `IncrementalIntegrator::formTangent` (`:100-125`);
2. the clearing of the SOE `factored` flag — done in each SOE subclass's `zeroA()`
   (`PARDISOGenLinSOE.cpp:560`, `UmfpackGenLinSOE.cpp:299`, …), not in `formTangent` itself —
   which forces the solver to redo its numeric factorization on the next `solve()`.

The engine **already has** the machinery to skip step 2: every sparse solver in the fork gates
its numeric phase on `theSOE->factored` (`UmfpackGenLinSolver`, `PARDISOGenLinSolver.cpp:378`,
`MumpsSolver`) — every *live* one; the unbuilt `PARDISOSymLinSolver` does not, but nothing
references it. The flag is cleared because the matrix is re-assembled. The missing piece is
one level up, at the assembly layer: **nothing asks whether the assembly was necessary.**

### Measured cost (reporter's numbers, `bdf8adf9f`)

SFIM shallow-footing model: 38 984 DOF, 9 248 `SSPbrickUP` (u-p) + `ShellMITC4`,
`Newmark(0.5,0.25)`, MKL PARDISO, 8 threads, i7-13700KF. Unit costs by differencing profiled
`Newton` against `ModifiedNewton` runs (both do 2 triangular solves per step; `Newton` does one
extra assembly and one extra factorization):

| Operation | Cost |
|---|---:|
| One tangent assembly | 325 ms |
| One numeric factorization (PARDISO, 8 threads) | 255 ms |
| One triangular substitution | 32 ms |
| Whole step, `ModifiedNewton` | 934 ms |

Head-to-head, same model / solver / thread count / tolerance, 25 profiled push steps:

| `algorithm` | `formTangent` calls/step | iters/step | s/step | settlement (m) |
|---|---:|---:|---:|---|
| `ModifiedNewton -initial` | 1.00 | 2.00 | **0.936** | −0.000974399 |
| `Newton -initial` | 2.00 | 2.00 | **1.504** | −0.000974399 |

Identical converged settlement, identical iteration count, **1.61x the cost**. On the
reporter's earlier UMFPACK build the same comparison was 6.56 vs 3.88 s/step — the same ratio
the operation counts predict. Note this is a **transient run with `betaK != 0`**, i.e. the case
in which the re-formed matrix genuinely *does* differ between iterations. Even there it buys
nothing: same answer, same iteration count.

### The `-initial`-is-nearly-inert observation (reporter, worth recording)

With `rayleigh $a $b 0 0`, the `-initial` flag is very nearly inert. In the reporter's model
(`Newmark(0.5,0.25)`, `dt=1e-3`, `zeta=0.5`, periods 0.5 s / 0.2 s): `c1 = 1`, `c2 = 2000`,
`b = 0.02274`, so the effective tangent is

```
INITIAL_TANGENT :  1.0 * Ki  +  45.5 * Kt  +  4.02e6 * M
CURRENT_TANGENT :             46.5 * Kt  +  4.02e6 * M
```

`-initial` exchanges about **2%** of the stiffness content *in that model*; the rest re-enters
through `C = ... + betaK*Kt` regardless.

> [!warning] Do not generalise the 2%. It is a two-extreme artifact requiring **both** a large
> damping ratio and a small time step. Same periods at `zeta = 0.05`, `dt = 0.01`:
> `b = 0.00227`, `c2 = 200`, `c2*b = 0.455` — so `-initial` exchanges about **69%** of the
> stiffness content and is doing very nearly what the user expects. The transferable guidance is
> *compute `c2*betaK` for your own model*, not "`-initial` is inert under Rayleigh damping".
> An earlier draft of this ADR, and the ledger entry derived from it, stated the general form.

## 2. What we verified — and the one correction

Every structural claim in the report was checked against the source at `bdf8adf9f`:

| Claim | Verdict | Evidence |
|---|---|---|
| `formTangent` inside the loop; only `INITIAL_THEN_CURRENT` special-cased | ✅ | `NewtonRaphson::solveCurrentStep`, the `do`-loop |
| `ModifiedNewton` forms once before the loop | ✅ | `ModifiedNewton::solveCurrentStep`, before the `do`-loop |
| Static `INITIAL_TANGENT` == `zeroTangent(); addKiToTang();` | ✅ | `StaticIntegrator.cpp:75-76` |
| Newmark `INITIAL_TANGENT` == `c1*Ki + c2*C + c3*M` | ✅ | `Newmark.cpp:295-298` |
| `C` carries `betaK*Kt`, so `betaK != 0` ⇒ genuine change | ✅ | `Element.cpp:222-223` |
| Solvers gate the numeric phase on `factored` | ✅ | `PARDISOGenLinSolver.cpp:378` |
| `OPS_ModifiedNewton` reads exactly one option | ✅ | `OPS_ModifiedNewton` (pre-fix) |
| PARDISO carries no profiler scopes | ✅ | zero `OPS_PROFILE` in `SRC/.../pardiso/` |
| **"Static analyses — redundant, always"** | ❌ **over-claimed** | see below |

### ⚠ Correction: "static ⇒ invariant" is not sound

The report's §3 asserts that under `StaticIntegrator` the assembled initial tangent "is
invariant not merely within a step but for the whole analysis". That is true of the *integrator*
layer — `formEleTangent` really is just `zeroTangent(); addKiToTang();` — but it is **not** true
of the layer underneath. `addKiToTang()` reaches `Element::getInitialStiff()`, and in OpenSees
that name is a convention, not a contract: several implementations return *material* initial
stiffness evaluated on the **current** configuration or the **current** active set.

Confirmed by reading the source:

- **`CorotCrdTransf3d::getInitialGlobalStiffMatrix`** (`:1452-1465`) triple-products the basic
  stiffness through the member `T`, and `T` is recomputed by `compTransfMatrixBasicGlobal()` at
  the tail of `update()` (`:549`). So any `dispBeamColumn`/`forceBeamColumn` on a **3D
  corotational** transformation has a `getInitialStiff()` that changes every iteration.
- **`LadrunoContactFE::addKiToTang`** (`:1590-1650`) deliberately mirrors `addKtToTang`,
  including the augmented **active mask** and the friction cone — both functions of the current
  gaps. Fork contact is configuration-dependent by construction.
- **`updateMaterialStage` is invisible to the topology stamp.** Neither `Domain::updateParameter`
  overload calls `domainChange()` (`Domain.cpp:2426-2457`), and nor does
  `TclUpdateMaterialCommand` / `MaterialStageParameter`. *Caveat on the example:* the obvious
  candidate, `PressureIndependMultiYield::getInitialTangent`, uses `refShearModulus`/
  `refBulkModulus`, which the stage switch does **not** touch — so PDMY/PIMY is a weak
  illustration. The real one is `ManzariDafalias::updateParameter` responseID 1, which flips
  `mElastFlag` and calls `Elastic2Plastic()` with no `domainChange()`. Same pattern in
  `ASDAbsorbingBoundary2D/3D` and `LysmerTriangle`, whose `getInitialStiff` branches on a
  `stage` written via `setParameter`.

### The systemic finding — it is far worse than three special cases

The adversarial review of this ADR (2026-07-25) turned the correction above from "three
awkward exceptions" into "the assumption is false by default". Verified directly:

- **`NDMaterial::getInitialTangent()` IS `getTangent()`.** `SRC/material/nD/NDMaterial.h:64`:
  ```cpp
  virtual const Matrix &getInitialTangent(void) {return this->getTangent();};
  ```
  The base-class default for "initial tangent" is *the current tangent*. Any nD material that
  does not override it hands every solid element a fully state-dependent initial stiffness.
  Contrast `UniaxialMaterial::getInitialTangent` (`UniaxialMaterial.h:68`), which is **pure
  virtual** — uniaxial materials are forced to implement it, which is why fiber/truss models
  are comparatively well behaved and solid models are not.
- **~41 element `getInitialStiff()` implementations are literally `return getTangentStiff();`**
  (count by grep over `SRC/element/`). Among them the `SSPquad`/`SSPbrick`/`SSPquadUP`/
  `SSPbrickUP` family — which is *exactly what the reporter's own SFIM model is built from*.
- Even materials that *do* override it can return a state-updated member: `stressDensity`
  recomputes `initialTangent` from the current trial mean stress inside `getCurrentStress()`;
  the `ManzariDafalias` family returns `mCe`, which the integrator rewrites on every
  `setTrialStrain`; `PM4Sand`/`PM4Silt` rebuild `mCe` in `commitState()`.
- **`Joint2D`/`Joint3D::getInitialStiff` call `theSprings[i]->getTangent()`**, not
  `getInitialTangent()` (`Joint2D.cpp:1135`) — the body is otherwise identical to
  `getTangentStiff()`.
- **`CorotCrdTransf3d`'s `T` is `static`** (`CorotCrdTransf3d.h:135`), i.e. shared by every
  instance, and `DispBeamColumn3d::getInitialStiff` never calls `crdTransf->update()`. So it
  triple-products through whatever `T` the *last element to update* left behind — cross-element
  aliasing, a latent vanilla bug independent of this ADR and worth reporting upstream.

**Consequence, and it is the main design consequence of this ADR:** "the initial tangent is
invariant" is not a mild over-claim, it is false for the default configuration of the very
element/material families the fork's production regime uses. See §4.4 — this is what forces
`isInitialStiffInvariant()` to default **false**.

**Second consequence, user-facing and independent of any code change:** for a solid element
whose nD material does not override `getInitialTangent()`, `algorithm Newton -initial` is not
performing initial-stiffness iteration at all — it silently *is* full Newton, at full cost.
That belongs in the R1 documentation on its own merits. See [[LEDGER_quirks]].

Confirmed **clean**, i.e. genuinely invariant, and worth recording so the audit does not
re-derive them:

- `CorotCrdTransf2d::getInitialGlobalStiffMatrix` — uses only `cosTheta`/`sinTheta`/`L`, all set
  in `initialize()` (`:318-319`); the current-configuration `cosAlpha`/`sinAlpha` are not
  touched.
- `PDeltaCrdTransf3d::getInitialGlobalStiffMatrix` — small-displacement `T_bl` built from `L`
  only; no geometric term, no current displacements.
- `CorotTruss::getInitialStiff` — its `R` is fixed in `setDomain` (`:512-532`), not in `update`.

**Consequence for the fix.** Initial-tangent invariance is **not an integrator-family property**
and cannot be decided by a `StaticIntegrator`-returns-`true` predicate as R3 proposes. It is an
**element/transformation** property that has to be declared and audited. A blanket static skip
would silently change the iteration path on exactly the geometrically-nonlinear and
contact models where `-initial` gets reached for in the first place.

**What is *not* at risk — stated carefully, because the loose version is false.** For a
**force-residual** convergence test (`NormUnbalance`, `NormDispIncr`'s residual siblings), and
*conditional on the step converging*, `K_f` sets only the path to equilibrium and not the
equilibrium, so a skipped re-assembly cannot move a converged answer.

> [!warning] The loose claim "`K_f` can never change the answer" is FALSE, and an earlier draft
> of this ADR shipped it. `NormDispIncr`, `RelativeNormDispIncr` and `EnergyIncr` — among the
> most-used tests in real decks — accept on `dU = K_f^-1 R`, not on `R`. A different `K_f`
> therefore accepts at a *different point*; a stiffer-than-appropriate one yields a small `dU`
> with a non-small `R`, which is the classic false convergence. Changing the path also changes
> whether a step fits inside `maxIter`, and in an adaptive-substepping driver that changes the
> load path and hence the answer. Any acceptance battery for this lane must therefore be run
> under a displacement-based test as well as a force-based one — §6's smoke uses
> `NormUnbalance`, i.e. it dodges precisely the case where the guarantee fails.

## 3. Decision

Ranked by value over risk, mirroring the reporter's own ranking.

| # | Item | Status |
|---|---|---|
| R1 | Document the trap | **shipped** (this PR) |
| R4 | `OPS_ModifiedNewton` parser: loop over all options | **shipped** (this PR) |
| R2 | Tangent-version counter | **designed here, not implemented** |
| R3 | `initialTangentIsInvariant()` shortcut | **rejected as specified**, folded into R2 |

R3 is rejected *as specified* — not as an idea. Its predicate (`true` for `StaticIntegrator`,
`true` for transient iff no `betaK != 0`) is unsound for the reasons in §2, and it would be a
silent behaviour change on corot/contact/PDMY models. What survives of R3 — a cheap early-out
that does not require the full counter — is recoverable as a special case of R2 once the
element-level predicate exists, so there is no reason to build both.

### R1 — what shipped

Documentation only, no behaviour change:

- A block comment at the `NewtonRaphson.cpp` `formTangent` site stating the per-iteration cost,
  pointing at `ModifiedNewton -initial` as the cheap spelling, and recording *why* the obvious
  fix is not applied there (the §2 correction).
- A [[LEDGER_quirks]] entry for the trap and for the `-initial`-is-nearly-inert-under-`betaK`
  observation.
- This ADR, which carries the measured cost table and the user-facing guidance.

The user-facing rule, stated plainly:

> **Want initial-stiffness iteration? Use `algorithm ModifiedNewton -initial`.**
> It is the same iteration `Newton -initial` performs — an initial-stiffness fixed-point — at
> one assembly and one factorization per step instead of one per iteration. Add `-factoronce`
> for one per *analysis*. `algorithm Newton -initial` is the expensive spelling of the same
> algorithm, and users reach for it as *the* robust fallback without knowing that.

### R4 — what shipped

`OPS_ModifiedNewton` (`ModifiedNewton.cpp:53+`) read exactly one option string through an
`if (OPS_GetNumRemainingInputArgs() > 0)` + `else if` chain, so
`algorithm ModifiedNewton -initial -factoronce` silently honoured `-initial` and dropped the
rest. Changed to the `while` form that `OPS_NewtonRaphsonAlgorithm()` already uses. The
combination is the natural way to express "initial stiffness, assembled and factorized once",
which for a static analysis is arguably the cheapest robust algorithm the framework offers.

The first cut of this change was `if` → `while` and nothing else. Adversarial review
(2026-07-25) found four defects in it, all now fixed:

1. **A new null-algorithm path.** `-hall`'s factor read `return 0`s on a parse failure. Under
   the old `if`, that was only reachable with `-hall` as the *first* option; the loop made it
   reachable from any position, so command lines that previously worked could now yield null.
   That is worse than it sounds: openseespy's `OPS_Algorithm` does
   `if (theAlgo != 0) { setAlgorithm(theAlgo); }` and then reports **success**
   (`OpenSeesCommands.cpp:2050`), so a null is silently discarded, the *previous* algorithm
   stays in force, and Python raises nothing. Now warns and keeps the 0.1/0.9 defaults — the
   same rule ADR-75 P1d already records for `-krylov` ("parse failure degrades to the default
   rather than returning 0").
2. **`-factorOnce` (camelCase) was still silently dropped.** `Linear.cpp:67`,
   `ExpressNewton.cpp:78` and both `commands.cpp` sites all accept `-factorOnce`, and it is the
   spelling [[LEDGER_quirks]] itself uses — so the likelier carry-over spelling was being
   dropped by the very fix meant to stop options being dropped. Now accepted, along with
   `-Initial`/`-Secant` for parity with `NewtonRaphson`.
3. **`-secant`/`-initial` did not reset `iFactor`/`cFactor`** the way
   `OPS_NewtonRaphsonAlgorithm()` does. Harmless under an `if` where only one branch could run;
   a loop is precisely the context where a later option must override an earlier one's factors.
   Now reset.
4. **Unknown tokens vanished silently.** Now warned. This also makes the `-hall` arity quirk
   self-diagnosing (see below).

Deliberately **not** fixed: `-hall`'s trailing-factor read still gates on
`OPS_GetNumRemainingInputArgs() == 2`, so in `-hall 0.2 0.8 -factoronce` the factors are
ignored and the defaults used. `NewtonRaphson` carries the identical gate, and fixing it here
alone would leave two parsers for the same flag disagreeing — which is the failure mode this
fork's ledger keeps recording. The unknown-token warning now surfaces the orphaned numerics, so
it is visible rather than silent. Fixing both parsers together is a candidate follow-up.

Honest statement of the behaviour change: this is **not** behaviour-preserving in general. Any
deck passing two or more recognised options now gets the last one to set a given field rather
than the first, which changes the iteration matrix. That is the intent, but it is a silent
results change for such decks. No in-tree deck passes more than one option, so the in-repo blast
radius is zero; single-option decks are byte-identical.

`OPS_ModifiedNewton` is the single choke point for both interpreters: the classic Tcl
`specifyAlgorithm` (`SRC/tcl/commands.cpp:4467`) and the Python/interpreter path
(`SRC/interpreter/OpenSeesCommands.cpp:2017`) both call it. `SRC/runtime/commands/analysis/algorithm.cpp`
parses `algorithm` too and drops `factorOnce` entirely, but that tree appears in **no**
`add_subdirectory()` (`SRC/CMakeLists.txt:8-30`) and is dead in this build — left untouched.

> [!warning] `-initial -factoronce` inherits an existing trap
> `factorOnce` has **no `domainChanged` reset** — `SolutionAlgorithm::domainChanged()` is
> virtual but `ModifiedNewton` does not override it, and the latch only ever resets on a
> convergence *failure* (`ModifiedNewton::solveCurrentStep`, the `result == -2` branch). After a mid-run domain change the SOE is
> re-sized/zeroed but the tangent re-form is skipped ⇒ solve against a stale matrix. Do not
> combine with element removal (ADR-51), contact re-emission (ADR-60), or staged construction.
> This is the same trap [[LEDGER_quirks]] already records for `algorithm Linear -factorOnce`;
> R4 makes it reachable from one more spelling. Arming a `domainChanged()` override on
> `ModifiedNewton` is a candidate follow-up, deliberately out of scope here because it changes
> existing `-factoronce` behaviour.

## 4. R2 design — the tangent version counter

### 4.1 The question to answer

An algorithm wants to ask: *is the matrix currently sitting in the SOE still the matrix that
`formTangent(statFlag, iFact, cFact)` would assemble right now?* That is a function of a tuple:

```
(statFlag, iFact, cFact, integrator coefficients, model state)
```

So the API is not a bare counter but a counter **parameterized by the tangent request**:

```cpp
// IncrementalIntegrator.h  — Ladruno ADR-76
//
// A monotone stamp for the matrix that formTangent(statFlag, iFact, cFact)
// would assemble RIGHT NOW. If two calls return the same non-zero stamp, the
// assembled matrices are bit-identical and the assembly may be skipped.
//
// 0 is reserved and means "unknowable / always re-form" — the SAFE default,
// and what every integrator returns until it has been audited.
virtual unsigned long long tangentVersion(int statFlag, double iFact, double cFact);
```

The algorithm side is then a three-line guard, with no behaviour change when the stamp is 0:

```cpp
unsigned long long v = theIntegrator->tangentVersion(tangent, iFactor, cFactor);
if (v == 0 || v != lastTangentVersion) {
    if (theIntegrator->formTangent(tangent, iFactor, cFactor) < 0) { ... }
    lastTangentVersion = v;
}
```

`lastTangentVersion` must be reset to 0 in `SolutionAlgorithm::domainChanged()` — belt and
braces alongside the topology term below, and the thing `factorOnce` never got.

### 4.2 What the stamp must fold in

| Invalidator | Static (ΣKi) | Transient (c1·Ki + c2·C + c3·M) | Already observable? |
|---|---|---|---|
| Topology / SOE resize (`Domain::domainChange()`) | ✅ | ✅ | `currentGeoTag` — but see §4.3 |
| Δt change (c1,c2,c3) | — | ✅ | integrator-local |
| `iFact`/`cFact` change (HALL) | ✅ | ✅ | argument |
| `statFlag` change | ✅ | ✅ | argument |
| `updateMaterialStage` | ✅ | ✅ | ❌ **no** — see §4.3 |
| `setParameter`/`updateParameter` on a stiffness or density | ✅ | ✅ | ❌ no |
| `rayleigh` re-issued | — | ✅ | ❌ no |
| Element `update()` (new displacements) | ✅ **iff** any element's `getInitialStiff` is configuration-dependent | ✅ iff that, **or** any element has `betaK != 0` | ❌ no |
| `commitState()` | as above | as above | ❌ no |

So the stamp is a fold of three independent terms:

```
tangentVersion = hash( topologyTag,          // Domain::currentGeoTag
                       coefficientTag,       // integrator: dt, statFlag, iFact, cFact, rayleigh
                       stateTag )            // bumped iff the request is state-sensitive
```

and `stateTag` is included **only when the request is state-sensitive**, which is the whole
point:

- `CURRENT_TANGENT` / `CURRENT_SECANT` / `HALL_TANGENT` — always state-sensitive. `stateTag`
  moves on every `update()`, so the stamp always differs and full Newton is untouched. Good:
  the guard must be provably inert for the default algorithm.
- `INITIAL_TANGENT`, static — state-sensitive **iff** any element declares a
  configuration-dependent initial stiffness.
- `INITIAL_TANGENT`, transient — state-sensitive iff that, **or** any element carries
  `betaK != 0`, **or** a `Damping` object with a state-dependent law is attached.

### 4.3 Concrete plumbing constraints found in the source

Three things a naive implementation gets wrong; all verified:

1. **`Domain::hasDomainChanged()` is not a pure query.** It consumes `hasDomainChangedFlag` and
   increments `currentGeoTag` as a side effect (`Domain.cpp`), and the analyses call it exactly
   once per step (`StaticAnalysis.cpp:150`, `DirectIntegrationAnalysis.cpp:206` et al. (~12 call sites, incl. two fork recorders)) to drive
   their own `domainStamp` comparison. A tangent-version scheme **must not call it** — it would
   steal the edge from the analysis. It needs a new pure getter (`int getCurrentGeoTag() const`),
   which is additive and harmless.

2. **`updateMaterialStage` is invisible to the topology stamp.** `MaterialStageParameter` never
   calls `Domain::domainChange()`. The classic PDMY/PIMY gravity-then-plastic workflow therefore
   changes `Ki` with no observable event. Either `MaterialStageParameter::update()` must bump a
   state tag, or the whole scheme must treat every `updateParameter` as an invalidator.

3. **`betaK != 0` must be a model-wide query, and it is not free.** `Element::betaK` is
   per-element and settable at any time via `rayleigh` / `region`. Recomputing "does any element
   carry `betaK != 0`" by looping the model on every iteration would eat the win. It has to be
   cached and invalidated where `setRayleighDampingFactors` is called.

### 4.4 The element predicate

```cpp
// Element.h — Ladruno ADR-76
//
// Does getInitialStiff() depend ONLY on data fixed at setDomain() time
// (initial geometry, initial material tangents)? I.e. is it invariant under
// update() / commitState()?
//
// Default FALSE. "Initial stiffness" is a naming convention in OpenSees, not a
// contract: NDMaterial::getInitialTangent() DEFAULTS to getTangent()
// (NDMaterial.h:64), and ~41 elements implement getInitialStiff() as a straight
// `return getTangentStiff();`. An element must therefore EARN this, per class,
// by audit. Reuse is opt-in per element, never inherited by silence.
virtual bool isInitialStiffInvariant(void) const { return false; }
```

> [!important] Reversed after adversarial review (2026-07-25)
> An earlier draft of this ADR defaulted this to **`true`**, arguing that a `getInitialStiff()`
> returning current-configuration stiffness is "the anomaly", and that default-`false` ships a
> feature that is dead until elements opt in. **That reasoning is dead.** §2 shows the anomaly
> is the *default*: `NDMaterial::getInitialTangent()` is defined as `getTangent()` in the base
> class, so the exception list would have had to enumerate most of `SRC/element/` — and any
> element the audit missed would have been silently, wrongly opted in. Under default-`false`, a
> missed element merely fails to get faster, which is the correct direction for the mistake to
> point. The "ships dead" objection is answered by scope instead: opt in the handful of classes
> that carry the fork's own production models, not the whole tree.

Opt-in candidates, in rough order of value: the fork's own solid elements paired with materials
that genuinely implement `getInitialTangent()`, fiber frame elements on linear/PDelta
transformations (`UniaxialMaterial::getInitialTangent` is pure virtual, so uniaxial-based
elements are the well-behaved population), and `CorotTruss`.

Known **disqualified** from §2: anything on `CorotCrdTransf3d`, `LadrunoContactFE`, the
`SSPquad`/`SSPbrick`/`SSPquadUP`/`SSPbrickUP` family, `Joint2D`/`Joint3D`, the `zeroLength`
contact family, `ASDShellQ4`/`ASDShellT3`, the `UWelements` contact set, and any element whose
`getInitialStiff()` calls `getTangentStiff()`.

**The audit across `SRC/element/` is the gate on P1, not P2.** With default-`false` it is no
longer a safety gate — nothing is unsafe until an element opts in — but it is still the entire
economic case for the lane, because it determines whether anything the fork actually runs can
opt in at all.

### 4.5 Rollout

- **P0 — audit.** Sweep `SRC/element/` for `getInitialStiff()` implementations that read current
  state. Deliverable: the override list, plus a `LEDGER_quirks` entry naming the pattern. No
  behaviour change. *This is the gate; nothing below ships without it.*
- **P1 — opt-in.** `tangentVersion()` + the algorithm guard + `Element::isInitialStiffInvariant`,
  reachable only via an explicit flag (`algorithm Newton -initial -reuseTangent`). Default build
  byte-identical, verified by the §6 battery.
- **P2 — default-on with opt-out**, following the ADR-67 mass-cache precedent (default-ON,
  conservatively invalidated, `-noMassCache` escape). Gated on P0 being complete and the §6
  battery green, including the corot/PDMY/contact regressions.
- **P3 — the static prize.** With the counter in place, a static `Newton -initial` pushover
  factorizes **once for the whole analysis** rather than once per iteration. On the reporter's
  arithmetic (1 000 steps x 5 iterations) that is 5 000 assemblies and 5 000 factorizations
  collapsing to 1. This is the payoff that justifies the counter over any narrower patch.

## 5. Interaction with ADR-75

The two reuse axes are complementary and must not be confused:

| Axis | Mechanism | Pays when |
|---|---|---|
| ADR-75 P1e `-krylov` | retained L/U as a **preconditioner** (PARDISO `iparm[3]`, phase 23) | A **has** changed — the full-Newton case |
| ADR-76 tangent version | **skip the assembly**, leaving `factored` set | A has **not** changed |

ADR-76 also fixes an asymmetry ADR-75 left standing: the `factored` flag can only be *exploited*
by an algorithm that refrains from re-assembling, and upstream ships exactly two such algorithms
(`Linear -factorOnce`, `ModifiedNewton`). ADR-76 is what lets the flag pay on a third.

Note the reporter's §7 instrumentation gap is real and unaddressed here by explicit scope
decision: `PARDISOGenLinSolver.cpp` has no `OPS_PROFILE` scopes where `UmfpackGenLinSolver.cpp`
has `soe.symbolic` / `soe.factor` / `soe.trisolve`, so on a PARDISO run `linearSolve` is one
opaque block and factorization cannot be separated from substitution without differencing two
algorithms. **P1 cannot be measured properly without it** — treat instrumentation parity as a
prerequisite of P1, tracked under ADR-75. The reporter also observed the profiler HDF5 run
attributes recording `threads=1` and `nElem=0` regardless of configuration; also unaddressed,
also worth a separate look.

## 6. Acceptance battery (from the reporter, adopted as-is)

Run with the P1 flag on; each must also be run with the flag off to confirm byte-identity with
the pre-ADR build.

1. **Static pushover, `Newton -initial`** — `formTangent` call count for the analysis equals 1
   (or the number of domain changes), and the displacement history is bit-identical to the
   current build.
2. **Transient, `betaK == 0`, `Newton -initial`** — `formTangent` calls == 1 per step; history
   identical to the current build.
3. **Transient, `betaK != 0`, `Newton -initial`** — behaviour **unchanged**, i.e. `formTangent`
   still called per iteration. *The reporter names this as the test they would most want to see,
   and we agree: it is the guard against an over-eager optimisation.*
4. **`algorithm ModifiedNewton -initial -factoronce`** — both options take effect. ✅ **shipped
   and covered by R4.**

### Test 4 — RUN, 8/8 checks passed (2026-07-25; re-run after the parser hardening)

`Ladruno_implementation/adr76_smoke/` — deck `adr76_factoronce_model.tcl` (classic
`OpenSees.exe`, `profiler` is registered in `SRC/tcl/commands.cpp:1219`) + checker
`adr76_factoronce_check.py` (reads the profiler HDF5 through
`Ladruno_tools/profiler_viewer/profiler_results.py`).

10 `Steel01` trusses in series (`Fy=100`, `E=1000`, `b=0.1`), pushed to `lambda=150` by an
unprofiled full-Newton phase so the chain is **past yield** before the profiled case starts.
That pre-yield is what makes the test possible at all: on a virgin model `Ki == Kt`, so
`-factoronce` and `-initial -factoronce` latch the *same* matrix and no observable separates
them.

| case | options | `formTangent` | iterations | tip ux |
|---|---|---:|---:|---|
| A | `-initial` | 5 | 910 | 6.999999999059916 |
| B | `-initial -factoronce` | **1** | 910 | 6.999999999059916 |
| C | `-factoronce -initial` | **1** | 910 | 6.999999999059916 |
| D | `-factoronce` | 1 | 5 | 7.0 |
| E | `-initial -factorOnce` | **1** | 910 | 6.999999999059916 |
| F | `-Initial -factorOnce` | **1** | 910 | 6.999999999059916 |

- **B vs A** — `-factoronce` survives a preceding `-initial`: 1 assembly instead of 5.
- **C ≡ B**, iteration for iteration — option order is irrelevant, which is the whole fix.
- **D vs B** — `-factoronce` alone latches the *current* (yielded, `0.1*Ki`) tangent and
  converges in 5 iterations; B/C latch the *initial* one and take 910. So `-initial` is not
  being dropped from `-factoronce -initial`.
- **E and F ≡ B** — the camelCase and capitalised spellings, added after the adversarial review
  found `-factorOnce` was still being dropped. Cases E/F would have FAILED the first cut of the
  R4 fix, which is exactly why they are in the battery.
- All converge to the same tip displacement (max deviation 9.4e-10 on 7.0, i.e. the
  convergence tolerance) — `K_f` sets the path, not the answer, exactly as §2 argues.

**Pre-fix expectations need no rebuild to confirm.** Before R4, `-initial -factoronce` *was*
`-initial` alone (= case A) and `-factoronce -initial` *was* `-factoronce` alone (= case D).
Both are measured in the table, and both differ from the post-fix B/C — so the smoke is
discriminating, and it fails loudly on a regression to the old parser.

**Incidental confirmation of the ADR's central claim.** Case A and case B ran the *identical*
910 iterations to the *identical* displacement, on 5 assemblies versus 1. That is the §1
argument reproduced in-repo at toy scale: the extra assemblies buy nothing.

Added by the §2 correction — these are the tests the reporter could not know to ask for:

5. **Static, `dispBeamColumn` on `CorotCrdTransf3d`, `Newton -initial`** — behaviour
   **unchanged** (the element declares itself non-invariant). This is the case that makes R3-as-
   specified wrong; it must fail loudly if the predicate ever regresses to integrator-family
   logic.
6. **Static, PDMY + `updateMaterialStage`, `Newton -initial`** — the stage switch invalidates;
   history identical to the current build.
7. **Fork contact (`LadrunoContactFE`) under `Newton -initial`** — behaviour unchanged.

## 7. Risks / open questions

> [!done] RESOLVED by adversarial review (2026-07-25)
> ~~Is `Element::isInitialStiffInvariant()` defaulting to `true` acceptable?~~ No. It defaults
> **false**, per §4.4. `NDMaterial::getInitialTangent()` is `getTangent()` in the base class, so
> default-`true` would have silently opted in most of `SRC/element/`.

> [!question]
> Given §2, is the lane worth running at all? The reporter's own SFIM model is `SSPbrickUP` +
> `ShellMITC4`, and `SSPbrickUP::getInitialStiff` is one of the ~41 `return getTangentStiff();`
> implementations — so on *that* model the reuse this ADR designs could never fire, and the
> `-initial` flag was never buying an initial-stiffness iteration in the first place. Which
> fork models would actually be able to opt in? That question, not the counter design, is what
> P0 has to answer first, and a negative answer closes the lane at zero further cost.

> [!question]
> Should `ModifiedNewton` get a `domainChanged()` override that re-arms `factorOnce`? It fixes
> a real staleness trap that R4 makes reachable from one more spelling — but it changes existing
> `-factoronce` behaviour, so it is a deliberate follow-up decision, not a drive-by.

> [!question]
> Parallel: `tangentVersion()` must agree across ranks or the SP/MP builds will diverge on which
> ranks skip. The topology term already does (`currentGeoTag` is synchronized); the state and
> `betaK` terms are local. Needs a reduction, or a refusal to reuse under
> `_PARALLEL_PROCESSING`, decided before P1.

- **Determinism.** Skipping assembly cannot perturb FP results (the matrix is bit-identical by
  construction), so this does not interact with the [[75_ladruno_sparse_direct_strategy_adr]]
  Lane-3 determinism policy.
- **Measurement.** Without the ADR-75 PARDISO profiler scopes, P1's win is only observable as
  wall-clock deltas between algorithms — the same differencing the reporter had to resort to.

## 8. Adversarial review outcome (2026-07-25) — R2 is WITHDRAWN as designed

Three reviewers were run against the shipped diff, the §2 audit, and the §4 design. §4 does not
survive. It is kept above as the record of what was tried; **this section supersedes its
rollout plan.**

### 8.1 Why the counter is not safely implementable as specified

- **The §4.2 invalidator set is missing at least ten entries**, and the sharpest is one this ADR
  warns about elsewhere: `Domain::activateElements`/`deactivateElements` (`Domain.cpp:4078-4112`)
  only flip `Element::is_this_element_active` and **never call `domainChange()`** — while every
  contributor gates on it (`FE_Element::addKiToTang` is `if (myEle != 0 && myEle->isActive())`).
  So staged construction removes whole blocks from `A` with `currentGeoTag`, `dt` and every
  element-level predicate unchanged. R2 would re-ship the exact `factorOnce` staleness trap
  under a new name, with default-on planned. Also missing: `modalDamping`/re-run `eigen` (a dense
  block added *before* the element loop, `TransientIntegrator.cpp:86-95`); `region … -rayleigh`
  and `-damping` (per-element writes that never touch `Domain::dmpBetaK`); the `mass` command
  (`Domain::setMass` returns without flagging); `Domain::addParameter`/`removeParameter` (their
  `domainChange()` calls are **commented out**); and `computeSensitivities()`, which re-forms `A`
  *between* `solveCurrentStep()` and `commit()`.
- **The predicate is hung off the wrong class.** `PenaltyMP_FE::getTangent` re-runs
  `determineTangent()` whenever `theMP->isTimeVarying()`, and `TransformationDOF_Group::getT()`
  varies likewise — these are `FE_Element`s with **no `Element` behind them**, so no
  `Element::isInitialStiffInvariant()` override can ever mark them non-invariant.
- **`getInitialStiff()` is state-dependent on the most vanilla elements in the tree the moment
  `-damping` is attached.** `ElasticBeam2d::getInitialStiff` multiplies by
  `theDamping->getStiffnessMultiplier()`, and `SecStifDamping`'s multiplier reads
  `theDomain->getCurrentTime()`, `getDT()`, an activation window, **and the global singleton
  `OPS_GetStaticAnalysis()`**. No per-element boolean can express "depends on wall-clock domain
  time and on whether a StaticAnalysis object exists", and no `stateTag` bumped on `update()`
  can see it.
- **§4.2's `betaK != 0` rule under-specifies `C`.** `Element::getDamp()` has four terms, not two:
  `alphaM*M`, `betaK*Kt`, **`betaK0*getInitialStiff()`** (`Element.cpp:224-225`) and
  **`betaKc*Kc`**. `betaK0` re-injects precisely the corot/contact/damping state-dependence of §2.
- **The stamp contract is self-contradictory.** §4.1 promises a *monotone* stamp where equality
  implies bit-identical matrices; §4.2 specifies `hash(...)`. A hash is not monotone, admits
  collisions, and can produce the reserved `0`. And `currentGeoTag` is not monotone either —
  `Domain::clearAll()` resets it and `setDomainChangeStamp(int)` sets it arbitrarily.
- **The `-krylov` interaction is antagonistic, not complementary** (contra §5). After a CGS win
  `factorsCurrent = false`, so on the next iteration a skipped assembly leaves `factored == true`
  but `factorsCurrent == false` — and `PARDISOGenLinSolver.cpp:314-315` re-enters **phase 23**
  rather than taking the phase-33 shortcut. Correctness is fine; a user arming both flags gets a
  **pessimization**.
- **Parallel is worse than §7 admitted.** The claim that `currentGeoTag` is synchronized is
  unsupported: `Domain::setDomainChangeStamp` has **zero callers in `SRC/`**. And an AND-reduction
  per iteration would be a collective inside a loop with five early `return -1` paths — any rank
  exiting early deadlocks the rest.

### 8.2 Why it is not worth rescuing — the value case inverts

Two findings kill the economics independently of the safety problems:

1. **The flagship measurement is on a model where R2 must never fire.** §1's 1.61x is a transient
   run with `betaK != 0`; §6 acceptance test 3 explicitly demands the guard *not* fire there. And
   the model is `SSPbrickUP` + `ShellMITC4` — and `SSPbrickUP::getInitialStiff` is one of the ~41
   `return getTangentStiff();` implementations, so `-initial` was never buying an
   initial-stiffness iteration on it in the first place. The reporter's problem is solved by R1
   (switch algorithm), which is documentation. R2 would have delivered them **nothing**.
2. **P3's prize already ships.** `ModifiedNewton -initial -factoronce` *already* gives one
   assembly and one factorization per analysis — §6 test 4 caseB measures exactly that
   (`formTangent` = 1 for the whole run). What R2 adds over it is automatic invalidation, and a
   `domainChanged()` override that re-arms `factorOnce` buys that in ~10 lines.

### 8.3 Revised decision

**R2: withdrawn.** Not deferred — withdrawn. Reopening it requires new evidence that some fork
model can actually opt in, which §2 makes unlikely.

**Replacement, and it is small:** carry the `factorOnce` re-arm in a fork algorithm class
(`LadrunoModifiedNewton`, fork's first `EquiSolnAlgo`, class tag in the ≥33000 band) that
overrides `domainChanged()`. That leaves vanilla `ModifiedNewton` byte-identical — no
years-old upstream behaviour changes — while giving the fork a correct
"initial stiffness, assembled and factorized once, invalidated properly" algorithm. Subclassing
needs `ModifiedNewton`'s members moved `private:` → `protected:`, a one-word vanilla edit with
no behavioural effect and a far smaller footprint than duplicating `solveCurrentStep`.

Pair it with an **additive one-time warning** on the vanilla `-factoronce` path: R4 widened the
population that can reach the staleness trap, and a fork class does not help users who do not
know it exists. Note the failure mode is solver-dependent and worst on the fork's own
recommendation: `BandGenLinSOE::setSize` zeroes `A`, so LAPACK fails loudly — but PARDISO
phase 22 on a zeroed matrix takes the perturbed-pivot path and returns `error = 0` with a
garbage solution. Silent-wrong precisely on PARDISO.

**Still worth doing, independent of all the above:** the `SRC/element/` `getInitialStiff()`
audit, not as an R2 gate but as the input to a **documentation** deliverable — users cannot
currently tell whether `-initial` does anything at all on their model. §2 and the
[[LEDGER_quirks]] entry are the beginning of that.

## 9. Provenance

- Issue report: TIMs project (SFIM continuum reference model), against `bdf8adf9f`, binaries of
  2026-07-25. Measurements in §1 are theirs.
- Vanilla files touched by R1/R4: see [[LEDGER_vanilla_files]].
- Quirks recorded: see [[LEDGER_quirks]].
- Family: [[40_ladruno_performance_adr]] (perf program) ·
  [[75_ladruno_sparse_direct_strategy_adr]] (the other reuse axis) ·
  [[40b_phase0_dominance_report]] (lane dominance).
