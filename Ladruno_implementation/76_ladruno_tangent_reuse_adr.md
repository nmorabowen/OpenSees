---
title: Tangent reuse — stop re-assembling and re-factorizing an unchanged matrix
project: Ladruno
status: closed — R1/R4 + LAPACK spin-off shipped; R2 WITHDRAWN; R3 rejected as specified
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

# ADR-76 — `algorithm Newton -initial` re-forms and re-factorizes a tangent nobody asked to change

> [!note] Rewritten 2026-07-26. The first version of this document was amended in
> place five times as adversarial reviews landed, which left ~270 lines of
> confidently-worded superseded plan ahead of the section that withdrew it. This
> version puts the conclusions first. The amended original is in git history
> (`ceff92c83` and earlier); no conclusion changed in the rewrite — only the order
> and the stale connective tissue. Two additions are new to this document but not to
> the record: §6 (the LAPACK spin-off, sourced from the ledger rows) and the §4.5
> implementation detail (sourced from [[_adr76_session_handoff]]).

> ADR-76. An **assembly/solve-lane** perf sub-ADR under [[40_ladruno_performance_adr]].
> Related to [[75_ladruno_sparse_direct_strategy_adr]] — ADR-75 made the *factorization*
> cheaper (PARDISO) and added reuse under a tangent that **has** changed (`-krylov`);
> this ADR examined the axis ADR-75's `factored` flag was built for but that nothing
> upstream ever pulls — not doing the work at all when the tangent has not changed. The
> two axes are **not** complementary in combination — see §4.3.
>
> Origin: an external issue report from the **TIMs project (SFIM continuum reference
> model)** against `bdf8adf9f`, with measured numbers on a 38 984-DOF
> `SSPbrickUP`/`ShellMITC4` model. Provenance in §7.

## 1. Outcome

| # | Item | Status |
|---|---|---|
| R1 | Document the trap | **shipped** — §5.1 |
| R4 | `OPS_ModifiedNewton` parser: loop over all options | **shipped** — §5.2, acceptance §5.3 |
| — | LAPACK singular-as-success fix (audit spin-off) | **shipped** — §6 |
| R2 | Tangent-version counter | **WITHDRAWN** — §4; superseded design in Appendix A |
| R3 | `initialTangentIsInvariant()` shortcut | **rejected as specified** — §4 |

Four conclusions carry the whole ADR:

1. **The premise was false.** The original report framed `Newton -initial`'s
   per-iteration re-assembly as work on a matrix that "cannot have changed". §3 shows
   the assembled initial tangent is state-dependent *by default* for most solid
   element/material pairs — `NDMaterial::getInitialTangent()` **is** `getTangent()` in
   the base class — so on those models `-initial` is not initial-stiffness iteration
   at all: it silently *is* full Newton, at full cost.

2. **The user-facing rule** (R1, and the durable deliverable of this ADR):

   > **Want initial-stiffness iteration? First check you are getting one at all** — §3:
   > on a solid element whose nD material does not override `getInitialTangent()`,
   > `-initial` silently IS full Newton. **If you are, `algorithm ModifiedNewton
   > -initial` is the cheap spelling**: the same iteration at one assembly and one
   > factorization per step instead of one per iteration. They are NOT the same
   > algorithm when the assembled initial tangent is state-dependent (corot3d,
   > contact, `-damping`, or any transient with `betaK`/`betaK0` != 0) — there
   > `Newton` re-linearizes every iteration and `ModifiedNewton` freezes at step
   > start. Add `-factoronce` for one per *analysis*. `algorithm Newton -initial` is
   > the expensive spelling of the same algorithm, and users reach for it as *the*
   > robust fallback without knowing that.

3. **The engine-side fix (R2) is withdrawn, not deferred** — §4. Its invalidator set
   is incomplete in ways that re-ship the known `factorOnce` staleness trap under a
   new name, its predicate cannot be expressed at the class it was hung on, and its
   flagship measurement is on a model where it must never fire. Reopening it requires
   new evidence that some fork model can actually opt in, which §3 makes unlikely. If
   it is ever reopened, the one design decision that survives review is that any
   element-level invariance predicate must default **false** (Appendix A.4 — kept
   because it is the best-argued passage of the withdrawn design).

4. **What replaces it is small** — §4.5: a fork algorithm class
   (`LadrunoModifiedNewton`) that re-arms `factorOnce` on `domainChanged()`, a
   one-time warning on the vanilla `-factoronce` path, and the `SRC/element/`
   `getInitialStiff()` audit repurposed as a *documentation* deliverable. Status:
   designed, **unstarted**.

## 2. Context — the reported defect

`algorithm Newton -initial` performs a full tangent assembly **and** a full numeric
factorization on **every iteration** of every step. Nothing anywhere asks whether the
re-assembly was necessary.

`formTangent` sits inside the iteration `do`-loop in `NewtonRaphson::solveCurrentStep`.
The only tangent flag given special treatment is `INITIAL_THEN_CURRENT_TANGENT`, which
correctly forms the initial tangent on iteration 0 and the current tangent thereafter.
`INITIAL_TANGENT` falls through to the generic branch and is re-formed unconditionally.
Contrast `ModifiedNewton::solveCurrentStep`, which forms the tangent **once before**
its loop.

> [!note] Citations in this ADR use function-name anchors, not line numbers, for
> `NewtonRaphson.cpp`, `ModifiedNewton.cpp`, and the churn-prone solver/SOE files. An
> earlier draft cited line spans that its own commit invalidated — including
> self-referentially, where the comment block displaced the very `do`-loop it pointed
> at. Line numbers are kept only for files this fork rarely touches, and were
> re-verified against the tree on 2026-07-26.

Each re-formation costs, in order:

1. `theSOE->zeroA()` plus a full loop over every `FE_Element` calling `getTangent()`
   and `addA()` — `IncrementalIntegrator::formTangent` (`:100-125`);
2. the clearing of the SOE `factored` flag — done in each SOE subclass's `zeroA()`
   (`PARDISOGenLinSOE::zeroA`, `UmfpackGenLinSOE::zeroA`, …), not in `formTangent`
   itself — which forces the solver to redo its numeric factorization on the next
   `solve()`.

The engine **already has** the machinery to skip step 2: every *live* sparse solver in
the fork gates its numeric phase on `theSOE->factored` (`UmfpackGenLinSolver`,
`PARDISOGenLinSolver::solve`, `MumpsSolver`; the unbuilt `PARDISOSymLinSolver` does
not, but nothing references it). The flag is cleared because the matrix is
re-assembled. The missing piece is one level up, at the assembly layer: **nothing asks
whether the assembly was necessary.**

### Measured cost (reporter's numbers, `bdf8adf9f`)

SFIM shallow-footing model: 38 984 DOF, 9 248 `SSPbrickUP` (u-p) + `ShellMITC4`,
`Newmark(0.5,0.25)`, MKL PARDISO, 8 threads, i7-13700KF. Unit costs by differencing
profiled `Newton` against `ModifiedNewton` runs (both do 2 triangular solves per step;
`Newton` does one extra assembly and one extra factorization):

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
reporter's earlier UMFPACK build the same comparison was 6.56 vs 3.88 s/step — the
same ratio the operation counts predict. Note this is a **transient run with
`betaK != 0`**, i.e. the case in which the re-formed matrix genuinely *does* differ
between iterations. Even there it buys nothing: same answer, same iteration count. (It
is also, per §3, a model on which `-initial` was never buying an initial-stiffness
iteration at all — the fact that ultimately killed R2's economics, §4.2.)

### The `-initial`-is-nearly-inert observation (reporter, worth recording)

With `rayleigh $a $b 0 0`, the `-initial` flag is very nearly inert. In the reporter's
model (`Newmark(0.5,0.25)`, `dt=1e-3`, `zeta=0.5`, periods 0.5 s / 0.2 s): `c1 = 1`,
`c2 = 2000`, `b = 0.02274`, so the effective tangent is

```
INITIAL_TANGENT :  1.0 * Ki  +  45.5 * Kt  +  4.02e6 * M
CURRENT_TANGENT :             46.5 * Kt  +  4.02e6 * M
```

`-initial` exchanges about **2%** of the stiffness content *in that model*; the rest
re-enters through `C = ... + betaK*Kt` regardless.

> [!warning] Do not generalise the 2%. It is a two-extreme artifact requiring **both**
> a large damping ratio and a small time step. Same periods at `zeta = 0.05`,
> `dt = 0.01`: `b = 0.00227`, `c2 = 200`, `c2*b = 0.455` — so `-initial` exchanges
> about **69%** of the stiffness content and is doing very nearly what the user
> expects. The transferable guidance is *compute `c2*betaK` for your own model*, not
> "`-initial` is inert under Rayleigh damping". An earlier draft of this ADR, and the
> ledger entry derived from it, stated the general form.

## 3. What the audit found — "initial tangent" is a convention, not a contract

Every structural claim in the report was checked against the source at `bdf8adf9f`:

| Claim | Verdict | Evidence |
|---|---|---|
| `formTangent` inside the loop; only `INITIAL_THEN_CURRENT` special-cased | ✅ | `NewtonRaphson::solveCurrentStep`, the `do`-loop |
| `ModifiedNewton` forms once before the loop | ✅ | `ModifiedNewton::solveCurrentStep`, before the `do`-loop |
| Static `INITIAL_TANGENT` == `zeroTangent(); addKiToTang();` | ✅ | `StaticIntegrator.cpp:75-76` |
| Newmark `INITIAL_TANGENT` == `c1*Ki + c2*C + c3*M` | ✅ | `Newmark.cpp:295-298` |
| `C` carries `betaK*Kt`, so `betaK != 0` ⇒ genuine change | ✅ | `Element.cpp:222-223` |
| Solvers gate the numeric phase on `factored` | ✅ | `PARDISOGenLinSolver::solve` |
| `OPS_ModifiedNewton` reads exactly one option | ✅ | `OPS_ModifiedNewton` (pre-fix) |
| PARDISO carries no profiler scopes | ✅ | zero `OPS_PROFILE` in `SRC/.../pardiso/` |
| **"Static analyses — redundant, always"** | ❌ **over-claimed** | see below |

### ⚠ Correction: "static ⇒ invariant" is not sound

The report's §3 asserts that under `StaticIntegrator` the assembled initial tangent
"is invariant not merely within a step but for the whole analysis". That is true of
the *integrator* layer — `formEleTangent` really is just
`zeroTangent(); addKiToTang();` — but it is **not** true of the layer underneath.
`addKiToTang()` reaches `Element::getInitialStiff()`, and in OpenSees that name is a
convention, not a contract: several implementations return *material* initial
stiffness evaluated on the **current** configuration or the **current** active set.

Confirmed by reading the source:

- **`CorotCrdTransf3d::getInitialGlobalStiffMatrix`** (`:1452-1465`) triple-products
  the basic stiffness through the member `T`, and `T` is recomputed by
  `compTransfMatrixBasicGlobal()` at the tail of `update()` (`:549`). So any
  `dispBeamColumn`/`forceBeamColumn` on a **3D corotational** transformation has a
  `getInitialStiff()` that changes every iteration.
- **`LadrunoContactFE::addKiToTang`** (`:1590-1650`) deliberately mirrors
  `addKtToTang`, including the augmented **active mask** and the friction cone — both
  functions of the current gaps. Fork contact is configuration-dependent by
  construction.
- **`updateMaterialStage` is invisible to the topology stamp.** Neither
  `Domain::updateParameter` overload calls `domainChange()` (`Domain.cpp:2426-2457`),
  and nor does `TclUpdateMaterialCommand` / `MaterialStageParameter`. *Caveat on the
  example:* the obvious candidate, `PressureIndependMultiYield::getInitialTangent`,
  uses `refShearModulus`/`refBulkModulus`, which the stage switch does **not** touch —
  so PDMY/PIMY is a weak illustration. The real one is
  `ManzariDafalias::updateParameter` responseID 1, which flips `mElastFlag` and calls
  `Elastic2Plastic()` with no `domainChange()`. Same pattern in
  `ASDAbsorbingBoundary2D/3D` and `LysmerTriangle`, whose `getInitialStiff` branches
  on a `stage` written via `setParameter`.

### The systemic finding — it is far worse than three special cases

The adversarial review of this ADR (2026-07-25) turned the correction above from
"three awkward exceptions" into "the assumption is false by default". Verified
directly:

- **`NDMaterial::getInitialTangent()` IS `getTangent()`.** `SRC/material/nD/NDMaterial.h:64`:
  ```cpp
  virtual const Matrix &getInitialTangent(void) {return this->getTangent();};
  ```
  The base-class default for "initial tangent" is *the current tangent*. Any nD
  material that does not override it hands every solid element a fully
  state-dependent initial stiffness. Contrast `UniaxialMaterial::getInitialTangent`
  (`UniaxialMaterial.h:68`), which is **pure virtual** — uniaxial materials are forced
  to implement it, which is why fiber/truss models are comparatively well behaved and
  solid models are not.
- **~41 element `getInitialStiff()` implementations are literally
  `return getTangentStiff();`** (count by grep over `SRC/element/`). Among them the
  `SSPquad`/`SSPbrick`/`SSPquadUP`/`SSPbrickUP` family — which is *exactly what the
  reporter's own SFIM model is built from*.
- Even materials that *do* override it can return a state-updated member:
  `stressDensity` recomputes `initialTangent` from the current trial mean stress
  inside `getCurrentStress()`; the `ManzariDafalias` family returns `mCe`, which the
  integrator rewrites on every `setTrialStrain`; `PM4Sand`/`PM4Silt` rebuild `mCe` in
  `commitState()`.
- **`Joint2D`/`Joint3D::getInitialStiff` call `theSprings[i]->getTangent()`**, not
  `getInitialTangent()` (`Joint2D.cpp:1141`) — the body is otherwise identical to
  `getTangentStiff()`.
- **`CorotCrdTransf3d`'s `T` is `static`** (`CorotCrdTransf3d.h:135`), i.e. shared by
  every instance, and `DispBeamColumn3d::getInitialStiff` never calls
  `crdTransf->update()`. So it triple-products through whatever `T` the *last element
  to update* left behind — cross-element aliasing, a latent vanilla bug independent
  of this ADR, queued for upstream report in [[upstream_pr_campaign]] (Wave 0).

**Consequence, and it is the main design consequence of this ADR:** "the initial
tangent is invariant" is not a mild over-claim, it is false for the default
configuration of the very element/material families the fork's production regime
uses. This is what forces any future invariance predicate to default **false**
(Appendix A.4) — and, with §4.2, what closed the lane.

**Second consequence, user-facing and independent of any code change:** for a solid
element whose nD material does not override `getInitialTangent()`,
`algorithm Newton -initial` is not performing initial-stiffness iteration at all — it
silently *is* full Newton, at full cost. That belongs in the R1 documentation on its
own merits. See [[LEDGER_quirks]].

Confirmed **clean**, i.e. genuinely invariant, and worth recording so a future audit
does not re-derive them:

- `CorotCrdTransf2d::getInitialGlobalStiffMatrix` — uses only
  `cosTheta`/`sinTheta`/`L`, all set in `initialize()` (`:318-319`); the
  current-configuration `cosAlpha`/`sinAlpha` are not touched.
- `PDeltaCrdTransf3d::getInitialGlobalStiffMatrix` — small-displacement `T_bl` built
  from `L` only; no geometric term, no current displacements.
- `CorotTruss::getInitialStiff` — its `R` is fixed in `setDomain` (`:512-532`), not
  in `update`.

**Consequence for any fix.** Initial-tangent invariance is **not an
integrator-family property** and cannot be decided by a
`StaticIntegrator`-returns-`true` predicate as R3 proposed. It is an
**element/transformation** property that has to be declared and audited — and §4.1
shows even that is not enough. A blanket static skip would silently change the
iteration path on exactly the geometrically-nonlinear and contact models where
`-initial` gets reached for in the first place.

**What is *not* at risk — stated carefully, because the loose version is false.** For
a **force-residual** convergence test (`NormUnbalance`, `NormDispIncr`'s residual
siblings), and *conditional on the step converging*, `K_f` sets only the path to
equilibrium and not the equilibrium, so a skipped re-assembly cannot move a converged
answer.

> [!warning] The loose claim "`K_f` can never change the answer" is FALSE, and an
> earlier draft of this ADR shipped it. `NormDispIncr`, `RelativeNormDispIncr` and
> `EnergyIncr` — among the most-used tests in real decks — accept on `dU = K_f^-1 R`,
> not on `R`. A different `K_f` therefore accepts at a *different point*; a
> stiffer-than-appropriate one yields a small `dU` with a non-small `R`, which is the
> classic false convergence. Changing the path also changes whether a step fits
> inside `maxIter`, and in an adaptive-substepping driver that changes the load path
> and hence the answer. Any acceptance battery for this lane must therefore be run
> under a displacement-based test as well as a force-based one — the shipped smoke
> (§5.3) uses `NormUnbalance`, i.e. it dodges precisely the case where the guarantee
> fails.

## 4. Why R2 was withdrawn and R3 rejected

Three reviewers were run against the shipped diff, the §3 audit, and the counter
design now archived as Appendix A. The design did not survive. R3
(`initialTangentIsInvariant()` as an integrator-family shortcut) fell first and
hardest: its predicate (`true` for `StaticIntegrator`, `true` for transient iff no
`betaK != 0`) is unsound for the reasons in §3, and what survives of the idea — a
cheap early-out — would only ever have been a special case of R2, which follows.

### 4.1 Safety — the counter is not implementable as specified

- **The invalidator set (Appendix A.2) is missing at least ten entries**, and the
  sharpest is one this ADR warns about elsewhere:
  `Domain::activateElements`/`deactivateElements` (`Domain.cpp:4080` ff) only flip
  `Element::is_this_element_active` and **never call `domainChange()`** — while every
  contributor gates on it (`FE_Element::addKiToTang` is
  `if (myEle != 0 && myEle->isActive())`). So staged construction removes whole
  blocks from `A` with `currentGeoTag`, `dt` and every element-level predicate
  unchanged. R2 would re-ship the exact `factorOnce` staleness trap under a new name,
  with default-on planned. Also missing: `modalDamping`/re-run `eigen` (a dense block
  added *before* the element loop, `TransientIntegrator.cpp:88` ff);
  `region … -rayleigh` and `-damping` (per-element writes that never touch
  `Domain::dmpBetaK`); the `mass` command (`Domain::setMass` returns without
  flagging); `Domain::addParameter`/`removeParameter` (their `domainChange()` calls
  are **commented out** — re-verified 2026-07-26; the `removeParameter` one is inside
  a dead `/* */` block); and `computeSensitivities()`, which re-forms `A` *between*
  `solveCurrentStep()` and `commit()`.
- **The predicate is hung off the wrong class.** `PenaltyMP_FE::getTangent` re-runs
  `determineTangent()` whenever `theMP->isTimeVarying()`, and
  `TransformationDOF_Group::getT()` varies likewise — these are `FE_Element`s with
  **no `Element` behind them**, so no `Element::isInitialStiffInvariant()` override
  can ever mark them non-invariant.
- **`getInitialStiff()` is state-dependent on the most vanilla elements in the tree
  the moment `-damping` is attached.** `ElasticBeam2d::getInitialStiff` multiplies by
  `theDamping->getStiffnessMultiplier()`, and `SecStifDamping`'s multiplier reads
  `theDomain->getCurrentTime()`, `getDT()`, an activation window, **and the global
  singleton `OPS_GetStaticAnalysis()`**. No per-element boolean can express "depends
  on wall-clock domain time and on whether a StaticAnalysis object exists", and no
  `stateTag` bumped on `update()` can see it.
- **The `betaK != 0` rule under-specifies `C`.** `Element::getDamp()` has four terms,
  not two: `alphaM*M`, `betaK*Kt`, **`betaK0*getInitialStiff()`**
  (`Element.cpp:224-225`) and **`betaKc*Kc`**. `betaK0` re-injects precisely the
  corot/contact/damping state-dependence of §3.
- **The stamp contract is self-contradictory.** Appendix A.1 promises a *monotone*
  stamp where equality implies bit-identical matrices; A.2 specifies `hash(...)`. A
  hash is not monotone, admits collisions, and can produce the reserved `0`. And
  `currentGeoTag` is not monotone either — `Domain::clearAll()` resets it and
  `setDomainChangeStamp(int)` sets it arbitrarily.
- **Parallel is worse than the open-questions list admitted.** The claim that
  `currentGeoTag` is synchronized is unsupported: `Domain::setDomainChangeStamp` has
  **zero callers in `SRC/`**. And an AND-reduction per iteration would be a
  collective inside a loop with five early `return -1` paths — any rank exiting early
  deadlocks the rest.

### 4.2 Economics — the value case inverts

Two findings kill the economics independently of the safety problems:

1. **The flagship measurement is on a model where R2 must never fire.** §2's 1.61x is
   a transient run with `betaK != 0`; acceptance test T3 (§4.4) explicitly demands
   the guard *not* fire there. And the model is `SSPbrickUP` + `ShellMITC4` —
   `SSPbrickUP::getInitialStiff` is one of the ~41 `return getTangentStiff();`
   implementations, so `-initial` was never buying an initial-stiffness iteration on
   it in the first place. The reporter's problem is solved by R1 (switch algorithm),
   which is documentation. R2 would have delivered them **nothing**.
2. **The prize already ships.** `ModifiedNewton -initial -factoronce` *already* gives
   one assembly and one factorization per analysis — §5.3 case B measures exactly
   that (`formTangent` = 1 for the whole run). What R2 adds over it is automatic
   invalidation, and a `domainChanged()` override that re-arms `factorOnce` buys that
   in ~10 lines (§4.5).

### 4.3 The ADR-75 interaction — antagonistic, not complementary

An earlier draft of this ADR presented the two reuse axes as complementary:

| Axis | Mechanism | Pays when |
|---|---|---|
| ADR-75 P1e `-krylov` | retained L/U as a **preconditioner** (PARDISO `iparm[3]`, phase 23) | A **has** changed — the full-Newton case |
| ADR-76 tangent-skip | **skip the assembly**, leaving `factored` set | A has **not** changed |

The mechanisms are distinct, but *arming both at once is a pessimization*, and the
review caught it: after a CGS win `factorsCurrent = false`
(`PARDISOGenLinSolver::solve`, the phase-23 block), so on the next iteration a
skipped assembly leaves `factored == true` but `factorsCurrent == false` — and the
solver re-enters **phase 23** rather than taking the phase-33 shortcut. Correctness
is fine; a user arming both flags gets CGS iterations on every solve of an unchanged
matrix instead of one triangular substitution. The same reasoning applies today to
`ModifiedNewton -factoronce` + `system Pardiso -krylov`: do not combine them.

Two instrumentation gaps the reporter surfaced remain real and are now tracked in
[[75_ladruno_sparse_direct_strategy_adr]] (open items): `PARDISOGenLinSolver` carries
no `OPS_PROFILE` scopes where `UmfpackGenLinSolver` has `soe.symbolic` / `soe.factor`
/ `soe.trisolve`, so on a PARDISO run `linearSolve` is one opaque block and
factorization cannot be separated from substitution without differencing two
algorithms; and the profiler HDF5 run attributes record `threads=1` and `nElem=0`
regardless of configuration.

One non-interaction worth recording: a skipped assembly is bit-identical by
construction (the matrix is simply not touched), so nothing in this lane —
including the §4.5 replacement — interacts with the ADR-75b Lane-3 determinism
policy.

### 4.4 What would have had to be proven — the acceptance bar R2 never met

The reporter proposed, and this ADR adopted, an acceptance battery for the counter.
With R2 withdrawn these are recorded as the bar any reopened design must clear, not
as a plan:

- **T1 — static pushover, `Newton -initial`:** `formTangent` call count for the
  analysis equals 1 (or the number of domain changes), displacement history
  bit-identical.
- **T2 — transient, `betaK == 0`, `Newton -initial`:** `formTangent` == 1 per step;
  history identical.
- **T3 — transient, `betaK != 0`, `Newton -initial`:** behaviour **unchanged**, i.e.
  `formTangent` still called per iteration. *The reporter named this as the test they
  most wanted, and we agree: it is the guard against an over-eager optimisation.*
- **T5 — static, `dispBeamColumn` on `CorotCrdTransf3d`, `Newton -initial`:**
  behaviour unchanged (the element declares itself non-invariant). This is the case
  that makes R3-as-specified wrong; it must fail loudly if a predicate ever regresses
  to integrator-family logic.
- **T6 — static, PDMY + `updateMaterialStage`, `Newton -initial`:** the stage switch
  invalidates; history identical.
- **T7 — fork contact (`LadrunoContactFE`) under `Newton -initial`:** behaviour
  unchanged.

(The battery's T4 — both `ModifiedNewton` options take effect — is not R2's: it gates
R4, ran 11/11, and is the acceptance record in §5.3.) Two requirements apply across
the whole battery: each test must also run with the feature **off** to confirm
byte-identity with the pre-change build, and — per the §3 warning box — each must
additionally run under a displacement-based convergence test, which the original
battery did not require.

### 4.5 What replaces it

**R2: withdrawn.** Not deferred — withdrawn. Reopening it requires new evidence that
some fork model can actually opt in, which §3 makes unlikely.

**Replacement, and it is small** (status: designed, unstarted):

- Carry the `factorOnce` re-arm in a fork algorithm class (`LadrunoModifiedNewton`,
  fork's first `EquiSolnAlgo`, class tag in the ≥33000 band) that overrides
  `domainChanged()`. That leaves vanilla `ModifiedNewton` byte-identical — no
  years-old upstream behaviour changes — while giving the fork a correct "initial
  stiffness, assembled and factorized once, invalidated properly" algorithm.
  Subclassing needs `ModifiedNewton`'s members moved `private:` → `protected:`, a
  one-word vanilla edit with no behavioural effect and a far smaller footprint than
  duplicating `solveCurrentStep`. The correct shape is an `invalidateTangent()`
  helper called from `domainChanged()` **and** from every negative-return path (the
  early-return holes are not accompanied by a domain change), plus a hardened
  `sendSelf` — all three `factorOnce` algorithms serialize the mutable latch.
  **Scope limit:** `domainChanged()` covers only the *topology* subset of
  invalidators — Δt change, staged construction, `updateMaterialStage`,
  `region -rayleigh` and modal damping never reach it (§4.1) — so the
  [[LEDGER_quirks]] warning survives the fix and must keep saying so.
- Pair it with an **additive one-time warning** on the vanilla `-factoronce` path: R4
  widened the population that can reach the staleness trap, and a fork class does not
  help users who do not know it exists. Note the failure mode is solver-dependent and
  worst on the fork's own recommendation: `BandGenLinSOE::setSize` zeroes `A`, so
  LAPACK fails loudly (§6 made sure of that) — but PARDISO phase 22 on a zeroed
  matrix takes the perturbed-pivot path and returns `error = 0` with a garbage
  solution. Silent-wrong precisely on PARDISO.
- **Still worth doing, independent of all the above:** the `SRC/element/`
  `getInitialStiff()` audit, not as an R2 gate but as the input to a
  **documentation** deliverable — users cannot currently tell whether `-initial` does
  anything at all on their model. §3 and the [[LEDGER_quirks]] entry are the
  beginning of that.

## 5. What shipped

### 5.1 R1 — documentation

Documentation only, no behaviour change:

- A block comment at the `NewtonRaphson.cpp` `formTangent` site stating the
  per-iteration cost, pointing at `ModifiedNewton -initial` as the cheap spelling,
  and recording *why* the obvious fix is not applied there (the §3 correction).
- A [[LEDGER_quirks]] entry for the trap and for the
  `-initial`-is-nearly-inert-under-`betaK` observation.
- This ADR, which carries the measured cost table and the user-facing rule (§1).

### 5.2 R4 — the option parser

`OPS_ModifiedNewton` read exactly one option string through an
`if (OPS_GetNumRemainingInputArgs() > 0)` + `else if` chain, so
`algorithm ModifiedNewton -initial -factoronce` silently honoured `-initial` and
dropped the rest. Changed to the `while` form that `OPS_NewtonRaphsonAlgorithm()`
already uses. The combination is the natural way to express "initial stiffness,
assembled and factorized once", which for a static analysis is arguably the cheapest
robust algorithm the framework offers.

The first cut of this change was `if` → `while` and nothing else. Adversarial review
(2026-07-25) found four defects in it, all fixed:

1. **A new null-algorithm path.** `-hall`'s factor read `return 0`s on a parse
   failure. Under the old `if`, that was only reachable with `-hall` as the *first*
   option; the loop made it reachable from any position, so command lines that
   previously worked could now yield null. That was worse than it sounds:
   openseespy's `OPS_Algorithm` did `if (theAlgo != 0) { setAlgorithm(theAlgo); }`
   and then reported **success**, so a null was silently discarded, the *previous*
   algorithm stayed in force, and Python raised nothing. Fixed at both layers: the
   parser now warns and keeps the 0.1/0.9 defaults (the same rule ADR-75 P1d records
   for `-krylov` — "parse failure degrades to the default rather than returning 0"),
   **and** the interpreter choke point now treats a null factory result as an error
   (`OpenSeesCommands.cpp:2060`, commit `074aff56f`) — which closes the
   silent-swallow for **every** algorithm factory at once and restores parity with
   the classic Tcl path, which always errored on null.
2. **`-factorOnce` (camelCase) was still silently dropped.** `Linear.cpp:67`,
   `ExpressNewton.cpp:78` and both `commands.cpp` sites all accept `-factorOnce`,
   and it is the spelling [[LEDGER_quirks]] itself uses — so the likelier carry-over
   spelling was being dropped by the very fix meant to stop options being dropped.
   Now accepted, along with `-Initial`/`-Secant` for parity with `NewtonRaphson`.
3. **`-secant`/`-initial` did not reset `iFactor`/`cFactor`** the way
   `OPS_NewtonRaphsonAlgorithm()` does. Harmless under an `if` where only one branch
   could run; a loop is precisely the context where a later option must override an
   earlier one's factors. Now reset.
4. **Unknown tokens vanished silently.** Now warned. This also makes the `-hall`
   arity quirk self-diagnosing (see below).

Deliberately **not** fixed: `-hall`'s trailing-factor read still gates on
`OPS_GetNumRemainingInputArgs() == 2`, so in `-hall 0.2 0.8 -factoronce` the factors
are ignored and the defaults used. `NewtonRaphson` carries the identical gate, and
fixing it here alone would leave two parsers for the same flag disagreeing — which is
the failure mode this fork's ledger keeps recording. The unknown-token warning now
surfaces the orphaned numerics, so it is visible rather than silent. Fixing both
parsers together is a candidate follow-up.

Honest statement of the behaviour change: this is **not** behaviour-preserving in
general. Any deck passing two or more recognised options now gets the last one to set
a given field rather than the first, which changes the iteration matrix. That is the
intent, but it is a silent results change for such decks. No in-tree deck passes more
than one option, so the in-repo blast radius is zero; single-option decks are
byte-identical.

`OPS_ModifiedNewton` is the single choke point for both interpreters: the classic
Tcl `specifyAlgorithm` (`SRC/tcl/commands.cpp:4467`) and the Python/interpreter path
(`SRC/interpreter/OpenSeesCommands.cpp:2017`) both call it.
`SRC/runtime/commands/analysis/algorithm.cpp` parses `algorithm` too and drops
`factorOnce` entirely, but that tree appears in **no** `add_subdirectory()`
(`SRC/CMakeLists.txt:8-30`) and is dead in this build — left untouched.

> [!warning] `-initial -factoronce` inherits an existing trap
> `factorOnce` has **no `domainChanged` reset** —
> `SolutionAlgorithm::domainChanged()` is virtual but `ModifiedNewton` does not
> override it, and the latch only ever resets on a convergence *failure*
> (`ModifiedNewton::solveCurrentStep`, the `result == -2` branch). After a mid-run
> domain change the SOE is re-sized/zeroed but the tangent re-form is skipped ⇒ solve
> against a stale matrix. Do not combine with element removal (ADR-51), contact
> re-emission (ADR-60), or staged construction. This is the same trap
> [[LEDGER_quirks]] already records for `algorithm Linear -factorOnce`; R4 makes it
> reachable from one more spelling. The §4.5 replacement plan
> (`LadrunoModifiedNewton`) is the designed fix; until it ships, this warning is the
> fix.

### 5.3 Acceptance record — the R4 smoke, 11/11 (2026-07-25, re-run after the parser hardening)

`Ladruno_implementation/adr76_smoke/` — deck `adr76_factoronce_model.tcl` (classic
`OpenSees.exe`, `profiler` is registered in `SRC/tcl/commands.cpp:1219`) + checker
`adr76_factoronce_check.py` (reads the profiler HDF5 through
`Ladruno_tools/profiler_viewer/profiler_results.py`).

10 `Steel01` trusses in series (`Fy=100`, `E=1000`, `b=0.1`), pushed to `lambda=150`
by an unprofiled full-Newton phase so the chain is **past yield** before the profiled
case starts. That pre-yield is what makes the test possible at all: on a virgin model
`Ki == Kt`, so `-factoronce` and `-initial -factoronce` latch the *same* matrix and
no observable separates them.

| case | options | `formTangent` | iterations | tip ux |
|---|---|---:|---:|---|
| A | `-initial` | 5 | 910 | 6.999999999059916 |
| B | `-initial -factoronce` | **1** | 910 | 6.999999999059916 |
| C | `-factoronce -initial` | **1** | 910 | 6.999999999059916 |
| D | `-factoronce` | 1 | 5 | 7.0 |
| E | `-initial -factorOnce` | **1** | 910 | 6.999999999059916 |
| F | `-Initial -factorOnce` | **1** | 910 | 6.999999999059916 |

- **B vs A** — `-factoronce` survives a preceding `-initial`: 1 assembly instead of
  5.
- **C ≡ B**, iteration for iteration — option order is irrelevant, which is the
  whole fix.
- **D vs B** — `-factoronce` alone latches the *current* (yielded, `0.1*Ki`) tangent
  and converges in 5 iterations; B/C latch the *initial* one and take 910. So
  `-initial` is not being dropped from `-factoronce -initial`.
- **E and F ≡ B** — the camelCase and capitalised spellings, added after the
  adversarial review found `-factorOnce` was still being dropped. Cases E/F would
  have FAILED the first cut of the R4 fix, which is exactly why they are in the
  battery.
- All converge to within 9.4e-10 of 7.0 — but note what that spread IS: case D
  (exact tangent) lands on 7.0 exactly, A/B/C/E/F on 6.999999999059916. A different
  `K_f` accepted at a *different point*, which is the mechanism the §3 warning box
  describes, not a refutation of it. It is also not "the convergence tolerance": the
  deck's tolerance is `NormUnbalance 1e-8` on FORCE.

**Pre-fix expectations need no rebuild to confirm.** Before R4,
`-initial -factoronce` *was* `-initial` alone (= case A) and `-factoronce -initial`
*was* `-factoronce` alone (= case D). Both are measured in the table, and both
differ from the post-fix B/C — so the smoke is discriminating, and it fails loudly
on a regression to the old parser. Since 2026-07-26 the checker also anchors every
case's tip displacement against the analytic 7.0 (not merely against the other
cases), so a regression driving all six cases to the same wrong answer fails too —
that is a 12th check added after the recorded 11/11 run; the table above is the
2026-07-25 battery and the anchor has not yet been exercised against a rebuilt
binary.

**Scope limit — this validates the OPTION PARSER and nothing else.** An earlier
draft called case A vs case B "incidental confirmation of the ADR's central claim".
It is not: every case runs `ModifiedNewton`, never `Newton`, so A-vs-B compares
per-step assembly against a per-analysis latch, not §2's per-*iteration* assembly.
And the profiled window is linear (`Steel01`, one hardening branch, monotonic load),
so A's 5 assemblies are provably of the same constant matrix and 910 == 910 is a
tautology. This deck cannot detect tangent staleness of any kind. The checker's
docstring said so; the ADR did not, until this correction.

## 6. Spin-off — the LAPACK band/full solvers reported SUCCESS on a singular matrix

Found while auditing the `factorOnce`/`domainChanged` failure paths for this ADR;
both ledgers record it as an "ADR-76 audit spin-off" and this section is its home in
the ADR family. It is a correctness fix, not a tangent-reuse item — it lives here
rather than in its own ADR because it embodies no open decision (the fix is a
one-token diff), and because its provenance *is* this ADR's audit: a separate ADR
number would add a document whose entire content already sits in the two ledger rows
and the regression deck.

**The defect.** `BandGenLinLapackSolver`, `FullGenLinLapackSolver` and
`BandSPDLinLapackSolver` mapped a positive LAPACK `info` with `return -info+1;`,
which C parses as `(-info)+1`. LAPACK sets `info = i` when `U(i,i)` is exactly zero,
so `info == 1` — a singularity at the **first pivot**, which is exactly what an
all-zero assembled `A` produces — returned **0**, i.e. SUCCESS. Every `info >= 2`
was already negative and detected; position, not singularity, decided whether the
defect fired.

**Why that is silent-wrong, not a return-code nit.** The solvers copy `B` into `X`
*before* calling LAPACK; DGBSV/DGESV/DPBSV do not compute `X` when `info > 0`; and
the algorithm gates on `theSOE->solve() < 0` — so the caller consumed the **load
vector** as a displacement increment. Measured pre-fix (3 free DOF, `Elastic` E=0,
load 7.0, `algorithm Linear`): `analyze` returned **0** with `nodeDisp == 7.0`.
Worse, the fact-check found dozens of call sites that never inspect the solve return
at all (`ArcLength`, `DisplacementControl`, `Newmark.cpp:1047`, the accelerator
family, …), so on those paths the singular solve was invisible by construction.

**The fix.** `return -info;` in all three — always negative in that branch, keeps
the failing pivot index in the magnitude, and is the idiom already used at
`SymBandEigenSolver.cpp:228` / `FullGenEigenSolver.cpp:188`. No caller anywhere in
`SRC/` compares the magnitude. Behaviour change to state plainly: decks that
previously limped past a singular step now fail loudly — which is the point, and it
will surface latent singularities (free nodes under `constraints Plain`,
zero-stiffness materials, fully-released members) that were being masked.

**Two dormant threaded solvers keep singular-as-success holes.**
`ThreadedSuperLU` (`:120`) forwards **pdgstrf**'s positive zero-pivot `info` via
`return info;` — the same shape. `BandSPDLinThreadSolver` never checks the worker
threads' factorization `info` at all after the `threadsDone` wait, so a singular
factor sails straight into the substitution. Neither is compiled by any CMake
target, so both are consequence-free today — recorded in [[LEDGER_quirks]] so the
next person to resurrect a threaded solver inherits the warning, not the bug.

**Regression.**
`Ladruno_implementation/lapack_singular_regression/lapack_singular_model.tcl` — 2
models x 4 solvers (`ProfileSPD` as the always-correct control), 8/8,
dependency-free, self-verifying. Two traps it encodes: the repro needs `n >= 2` (the
solvers' own 1x1 special-case guards intercept a single-DOF singular system before
LAPACK — the first version of the probe reported BandSPD as CLEAN for exactly that
reason), and the free-node model must number the zero row **first** (`info == 1`;
free-node-last gives `info == 3`, which the old code detected). Wired into the
Zone-A CI job (2026-07-26), gated on the deck's terminal marker string — **not** on
the exit code, because `OpenSees.exe` exits 0 after a Tcl parse error (recorded
quirk).

Full analysis, measured pre-fix numbers, and the call-site fact-check:
[[LEDGER_vanilla_files]] (the three solver rows) and [[LEDGER_quirks]].

## 7. Provenance

- Issue report: TIMs project (SFIM continuum reference model), against `bdf8adf9f`,
  binaries of 2026-07-25. Measurements in §2 are theirs. Source file:
  `ISSUE — Newton-initial redundant refactorization.md` — **external; not archived
  in this repo** (searched 2026-07-26: not on this machine). §2 reproduces its
  measured tables in full precisely so that nothing load-bearing lives only in the
  unarchived file.
- Session handoff: [[_adr76_session_handoff]] — includes corrections to earlier
  findings that must not be re-propagated (`LadrunoIMKBeam` is clean;
  `LadrunoDispBeamColumn` clean except `-damping`/corot3d; do NOT freeze `Ki` at the
  reference configuration — implemented and reverted: `betaK0*getInitialStiff()`
  enters the *residual* via `Element::getRayleighDampingForces`, so a frozen `K0`
  changes converged answers and diverges under `-initial` past ~2-8° of chord
  rotation).
- Vanilla files touched by R1/R4 + the LAPACK spin-off: see [[LEDGER_vanilla_files]].
- Quirks recorded: see [[LEDGER_quirks]].
- Family: [[40_ladruno_performance_adr]] (perf program; this ADR amends it — logged
  there 2026-07-26) · [[75_ladruno_sparse_direct_strategy_adr]] (the other reuse
  axis; carries the back-reference and the PARDISO instrumentation open item) ·
  [[40b_phase0_dominance_report]] (lane dominance) · [[upstream_pr_campaign]]
  (CorotCrdTransf3d static-`T` queued in Wave 0).

## Appendix A — the withdrawn R2 design (superseded record)

> [!warning] R2 WAS WITHDRAWN — this appendix is a record, not a plan.
> Everything here was written before the review in §4, which withdrew R2 outright.
> **Do not implement from this appendix.** It is kept because the invalidator
> analysis (A.2-A.3) is the starting inventory for any future attempt, and because
> A.4 is the one piece of the design that *survived* review as a constraint on
> reopening. The original P0-P3 rollout ladder is deleted, not archived — §4.2
> records the only fact from it that matters (the P3 prize already ships as
> `-initial -factoronce`).

In one paragraph: R2 proposed a monotone "tangent version" stamp on
`IncrementalIntegrator`, parameterized by the tangent request
`(statFlag, iFact, cFact)`, letting an algorithm skip `formTangent` when the stamp
had not moved; elements would declare initial-stiff invariance via an audited,
default-false predicate, and the guard would be provably inert for `CURRENT_TANGENT`
requests (whose stamp moves on every `update()`). §4.1 explains why the stamp cannot
be made sound and §4.2 why it would not have paid even if it could.

### A.1 The API sketch

```cpp
// IncrementalIntegrator.h  — Ladruno ADR-76 (WITHDRAWN, never shipped)
//
// A monotone stamp for the matrix that formTangent(statFlag, iFact, cFact)
// would assemble RIGHT NOW. If two calls return the same non-zero stamp, the
// assembled matrices are bit-identical and the assembly may be skipped.
//
// 0 is reserved and means "unknowable / always re-form" — the SAFE default,
// and what every integrator returns until it has been audited.
virtual unsigned long long tangentVersion(int statFlag, double iFact, double cFact);
```

with a three-line guard on the algorithm side and `lastTangentVersion` reset to 0 in
`SolutionAlgorithm::domainChanged()`. (§4.1: the monotone-stamp contract and the
hash-fold specification below contradict each other — one of several reasons this
never left the page.)

### A.2 The invalidator table

| Invalidator | Static (ΣKi) | Transient (c1·Ki + c2·C + c3·M) | Already observable? |
|---|---|---|---|
| Topology / SOE resize (`Domain::domainChange()`) | ✅ | ✅ | `currentGeoTag` — but see A.3 |
| Δt change (c1,c2,c3) | — | ✅ | integrator-local |
| `iFact`/`cFact` change (HALL) | ✅ | ✅ | argument |
| `statFlag` change | ✅ | ✅ | argument |
| `updateMaterialStage` | ✅ | ✅ | ❌ **no** — see A.3 |
| `setParameter`/`updateParameter` on a stiffness or density | ✅ | ✅ | ❌ no |
| `rayleigh` re-issued | — | ✅ | ❌ no |
| Element `update()` (new displacements) | ✅ **iff** any element's `getInitialStiff` is configuration-dependent | ✅ iff that, **or** any element has `betaK != 0` | ❌ no |
| `commitState()` | as above | as above | ❌ no |

The stamp was to fold `hash(topologyTag, coefficientTag, stateTag)` with `stateTag`
included only for state-sensitive requests. §4.1 lists at least ten invalidators
this table is missing — `activateElements`/`deactivateElements`, modal damping,
`region -rayleigh`, `-damping`, `mass`, `addParameter`/`removeParameter`,
`computeSensitivities()`, `betaK0`/`betaKc`, time-varying MP constraints — which is
what withdrew the design.

### A.3 Concrete plumbing constraints found in the source

Verified against source; still true and reusable by any future attempt:

1. **`Domain::hasDomainChanged()` is not a pure query.** It consumes
   `hasDomainChangedFlag` and increments `currentGeoTag` as a side effect, and the
   analyses call it exactly once per step (`StaticAnalysis.cpp:150`,
   `DirectIntegrationAnalysis.cpp:206` et al., ~12 call sites incl. two fork
   recorders) to drive their own `domainStamp` comparison. A tangent-version scheme
   **must not call it** — it would steal the edge from the analysis. It needs a new
   pure getter (`int getCurrentGeoTag() const`), which is additive and harmless.
2. **`updateMaterialStage` is invisible to the topology stamp.**
   `MaterialStageParameter` never calls `Domain::domainChange()`. The classic
   PDMY/PIMY gravity-then-plastic workflow therefore changes `Ki` with no observable
   event. Either `MaterialStageParameter::update()` must bump a state tag, or the
   whole scheme must treat every `updateParameter` as an invalidator.
3. **`betaK != 0` must be a model-wide query, and it is not free.** `Element::betaK`
   is per-element and settable at any time via `rayleigh` / `region`. Recomputing
   "does any element carry `betaK != 0`" by looping the model on every iteration
   would eat the win. It has to be cached and invalidated where
   `setRayleighDampingFactors` is called. (And per §4.1 it under-specifies `C`
   anyway: `betaK0` and `betaKc` are separate terms.)

### A.4 The element predicate — the piece that survived review

```cpp
// Element.h — Ladruno ADR-76 (WITHDRAWN with R2 — but the default is the
// binding precedent for any future invariance predicate)
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
> An earlier draft of this ADR defaulted this to **`true`**, arguing that a
> `getInitialStiff()` returning current-configuration stiffness is "the anomaly",
> and that default-`false` ships a feature that is dead until elements opt in.
> **That reasoning is dead.** §3 shows the anomaly is the *default*:
> `NDMaterial::getInitialTangent()` is defined as `getTangent()` in the base class,
> so the exception list would have had to enumerate most of `SRC/element/` — and any
> element the audit missed would have been silently, wrongly opted in. Under
> default-`false`, a missed element merely fails to get faster, which is the correct
> direction for the mistake to point. The "ships dead" objection is answered by
> scope instead: opt in the handful of classes that carry the fork's own production
> models, not the whole tree.

Known **disqualified** by the §3 audit: anything on `CorotCrdTransf3d`,
`LadrunoContactFE`, the `SSPquad`/`SSPbrick`/`SSPquadUP`/`SSPbrickUP` family,
`Joint2D`/`Joint3D`, the `zeroLength` contact family, `ASDShellQ4`/`ASDShellT3`, the
`UWelements` contact set, and any element whose `getInitialStiff()` calls
`getTangentStiff()`. Known clean: the §3 clean list (`CorotCrdTransf2d`,
`PDeltaCrdTransf3d`, `CorotTruss`), plus in principle the uniaxial-based population
(`UniaxialMaterial::getInitialTangent` is pure virtual). The opt-in inventory the
withdrawn rollout would have started from — and the starting point for the
"evidence some fork model can opt in" test §4.5 sets for reopening — was: the
fork's own solid elements paired with materials that genuinely implement
`getInitialTangent()`, fiber frame elements on linear/PDelta transformations, and
`CorotTruss`.
