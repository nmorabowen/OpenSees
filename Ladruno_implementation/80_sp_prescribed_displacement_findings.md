---
title: Non-homogeneous `sp` + path-dependent material — findings and strengthening candidates
project: Ladruno
status: findings (pre-ADR — evidence only, no decision taken)
priority: medium
owner: nmora
tags:
  - implementation
  - solver
  - constraints
  - integrator
  - cross-pollination
---

# Non-homogeneous `sp` × path-dependent material: what we found, and what could strengthen the fork

> **This is a FINDINGS note, not an ADR.** It records a mechanism traced in the
> source, the measurements that pin it, what Abaqus / Kratos / LS-DYNA do
> instead, and a ranked list of strengthening candidates with their validation
> gates. **No decision is taken and no code is changed.** The ADR is to be
> written on top of this.
>
> Source of the measurements: the Cerro Lindo rung-4 base-fuse campaign
> (`C:\nmb\My Libraries\Cerro Lindo\Informe No3\Models\Fuse FEM\04_solid_fuse\`,
> notes `NOTE_drive_control_RevA.md` §4-§6 and `NOTE_plasticity_RevA.md` §1).
> Date: 4 Aug 2026.

---

## 1. The symptom

A 3D solid model, small-strain `LadrunoJ2` (Gr50, `R_y F_y` = 379.5 MPa,
`H` = 2 000 MPa ⇒ `H/E` = 1 %), static monotonic push. A 180×180×20 cover
plate's outer face — 400 nodes — carries a **non-homogeneous prescribed
displacement** (`sp` inside a `Plain` pattern with a `Linear` series) scaled by
the `LoadControl` factor; opposite face fully fixed; `constraints Transformation`.

Giving the elements **at that driven face** a yielding material instead of an
elastic one costs an order of magnitude in solver work, for an identical answer.
All four runs pinned to `KrylovNewton` so the march is reproducible:

| run | drive on the loaded face | driven-face elements | increments | cutbacks | Newton iterations |
|---|---|---|---|---|---|
| `dcp_el` | prescribed displacement | elastic | 30 | 0 | 149 |
| **`dcp_cov`** | **prescribed displacement** | **plastic** | **324** | **52** | **4 255** |
| `lcp_el` | traction (Neumann) | elastic | 30 | 0 | 152 |
| **`lcp_cov`** | **traction (Neumann)** | **plastic** | **30** | **0** | **145** |

**×28.6 iterations and 52 cutbacks under prescribed displacement; ×0.95 —
nothing — under traction.** Converged answer identical: the device force at a
stated internal rotation is 920.0 vs 920.1 kN, **0.00 %**. Equilibrium residual
9.2e-08 (DC-plastic) against 2.2e-15 (LC-plastic).

⚠ The decisive control: in the **cheap** LC-plastic run the cover **grossly
plastically hinges** (its driven face sags 26.7 mm, 99 % of its own mean
stroke); in the **expensive** DC-plastic run the cover **never yields at all**
(max σ_vM 295 MPa against `f_y` = 379.5, 0.00 % plastic at every step). So
*cover yielding is not the mechanism* — a fully hinging cover converges in 4.6
iterations/increment while a non-yielding one takes 13.4 with 35 cutbacks.

---

## 2. The mechanism, traced in the source

Every reference verified in this tree.

1. **`LoadControl::newStep()` has no predictor.** It calls only
   `theModel->applyLoadDomain(currentLambda)` —
   `SRC/analysis/integrator/LoadControl.cpp:130`. Nothing runs before it.
2. `AnalysisModel::applyLoadDomain` = `myDomain->applyLoad(t)` +
   `myHandler->applyLoad()` — `SRC/analysis/model/AnalysisModel.cpp:556-568`.
   Neither calls `Domain::update()`.
3. `TransformationConstraintHandler::applyLoad()` → `enforceSPs()` —
   `SRC/analysis/handler/TransformationConstraintHandler.cpp:496-524`. It writes
   the full prescribed value into each constrained node
   (`TransformationDOF_Group::enforceSPs(1)` →
   `myNode->setTrialDisp(value, i)`, `TransformationDOF_Group.cpp:1071`;
   `:1069` on the `TRANSF_INCREMENTAL_SP` path) **and then, at lines 518-521,
   calls `theEle->updateElement()` on every element touching a constrained
   node** (`theFEs`, membership from `constrainedNodeSet`).
   `FE_Element::updateElement()` is `myEle->update()` —
   `SRC/analysis/fe_ele/FE_Element.cpp:1116-1125`.
4. ⇒ **The first constitutive evaluation of every increment happens inside
   `newStep`**, with the driven face advanced by the full increment and every
   interior node still at the last converged position. That element layer sees a
   strain increment of **Δδ/h_element** instead of the physical **Δδ/L_model** —
   an overstrain of order `L/h`.
5. Path-independent material: a harmless bad starting point, recovered in 2-3
   iterations. **Path-dependent material:** the layer yields spuriously, its
   consistent radial-return tangent collapses toward
   `2G·H/(3G+H)` ≈ 1 % of elastic (`H/E` = 1 %), the predictor built on it
   overshoots, the layer unloads elastically on the next iteration, and the
   active set oscillates. Nothing survives in committed state — trial states are
   recomputed from the committed state each iteration — so the **answer** is
   untouched while the **path** is wrecked. Conditioning, not accuracy.
6. **Why the pre-update exists at all:** the `sp`'d DOF is genuinely eliminated
   (`TransformationDOF_Group::addSP_Constraint` → `setID(dof,-1)` at `:1034`;
   every numberer assigns equations only to `-2` —
   `SRC/analysis/numberer/DOF_Numberer.cpp:156,298`,
   `PlainNumberer.cpp:107,221`). There is no column for a `K·Δu_prescribed`
   term, and `TransformationFE::getResidual` is only `Tᵗ·R` from element state
   (`SRC/analysis/fe_ele/transformation/TransformationFE.cpp:391-394`). So
   **pre-updating the elements is how the prescribed motion reaches the RHS at
   all.** It is a coherent design, not a bug — it just forces the constitutive
   evaluation into the one state the other codes avoid.

### Corollaries worth keeping

- **`TransformationFE::getTangForce()` is an unimplemented stub** — prints
  `not yet implemented`, zeroes and returns (`TransformationFE.cpp:447-453`).
  That is exactly the hook a `b −= K·Δu_D` route would need (§4 candidate C).
- **`DisplacementControl` already has the correct ordering.**
  `DisplacementControl::newStep` does
  `formTangent → solve → deltaU = dUhat·dλ → theModel->incrDisp(*deltaU) →
  theModel->applyLoadDomain(λ) → theModel->updateDomain()`
  (`SRC/analysis/integrator/DisplacementControl.cpp`, ~lines 35-98). The interior
  is advanced by a tangent-consistent predictor **before** the SPs are enforced
  and the elements updated. **The pattern the fix needs already exists in the
  tree; `LoadControl` simply has nothing in front of `applyLoadDomain`.**
- **`AutoConstraintHandler::applyLoad()` omits the `updateElement()` loop**
  (`SRC/analysis/handler/AutoConstraintHandler.cpp:573-584` runs the two
  `enforceSPs` loops and stops). With a non-homogeneous `sp` the first residual
  would then be stale, `dU ≈ 0`, and `CTestNormDispIncr` — which has **no
  minimum-iteration guard** (`SRC/convergenceTest/CTestNormDispIncr.cpp:160-175`,
  `if (norm <= tol) return currentIter;` with only a `currentIter == 0` check) —
  would report convergence at iteration 1. The increment would commit with only
  the boundary layer moved and equilibrium never re-checked: **a wrong answer,
  not an error.** Source-level inference; **not yet reproduced** (§5 gate 4).
- **Fork ownership context:** `TransformationConstraintHandler.cpp` and
  `TransformationDOF_Group.cpp` are already fork-touched (ADR-74 `handle()`
  O(1) membership rework, PR #595 — searches only, behaviour byte-identical), so
  they carry `// Ladruno` marks and ledger rows. `LoadControl.cpp`,
  `AutoConstraintHandler.cpp`, `CTestNormDispIncr.cpp` and `TransformationFE.cpp`
  are **still vanilla**. See [[LEDGER_vanilla_files]].
- **Adjacent precedent, same code path:** the fork already found and fixed
  `sp ... -subtractInit` being a **no-op on the openseespy path** — a staged deck
  silently yanked the node to the absolute value instead of moving it
  incrementally ("wrong answer, no warning"), gated by
  `tests/test_sp_subtract_init.py`. That row also enumerates the only handlers
  that honour `SP_Constraint::getInitialValue()`: `TransformationDOF_Group:1068`,
  `PenaltySP_FE:125`, `LagrangeSP_FE:143`. See [[LEDGER_vanilla_files]].

---

## 3. What Abaqus, Kratos and LS-DYNA do instead

One principle, three implementations: **the first linear solve of an increment
should be what distributes the boundary motion into the interior, and the
constitutive law should not be evaluated before that distribution has happened.**

| code | how the prescribed increment reaches iteration 1 | law evaluated in the lagging state? | named remedy |
|---|---|---|---|
| **OpenSees** | node written, then constrained elements **pre-updated** so their imbalance becomes the RHS | **yes, every increment** | none |
| **Abaqus/Standard** | not documented (the Theory Guide has no BC-imposition section) | **avoided by the predictor** | **100 % linear extrapolation of the previous increment**, default ON |
| **Kratos** | process writes the value into the nodal database; `ApplyDirichletConditions` zeroes the fixed row **and its RHS entry** — no `K·Δu_D` | **yes, and earlier** | **`use_old_stiffness_in_first_iteration`** |
| **LS-DYNA implicit** | as a **constraint-violation term `g` in the residual**, distributed by the global solve, then scaled by line search | documented equations say no | line search + 4 auto reform triggers + **`MAT_24`/`123` iterate rejection** |

**Abaqus** — `EXTRAPOLATION=LINEAR` is the default: 100 % linear extrapolation of
the previous incremental solution as the starting guess (Analysis User's Guide
Vol 2 §6.1.2 p. 6.1.2-6; §7.2.3 p. 7.2.3-10; 1 % for Riks), **off in the first
increment of a step**, and abandoned when the time increment changes sharply.
Material state is **not** extrapolated — history variables are re-integrated over
the whole increment every iteration (Theory Guide §2.2.1 p. 2.2.1-3), *exactly as
OpenSees does*. **The codes differ in the PREDICTOR, not in the material
update** — which is what localizes the fix. Abaqus's own caution is our mechanism
from the predictor side: extrapolation *"can cause Abaqus/Standard to iterate
excessively"* for *"abrupt changes in the load magnitudes or boundary
conditions"*.

**Kratos** sits in the *same* position as OpenSees and hits it slightly earlier,
but ships the principled remedy as a flag:
`use_old_stiffness_in_first_iteration` selects
`BuildAndSolveLinearizedOnPreviousIteration`
(`residualbased_block_builder_and_solver.h:552-658`) — for iteration 1 it frees
the fixed DOFs, rewinds the database to the **converged** state, assembles `K`
and `b` **there**, then adds `b −= K·Δu_D` and re-fixes. **Kratos's own
regression tests enable it for J2 plasticity, plastic-damage, high-cycle fatigue
and traction-separation, and disable it for elastic trusses and beams** — an
independent statement that the trigger is precisely "path-dependent material +
prescribed displacement". Kratos also has a true displacement-control device
(`DisplacementControlCondition3D1N`) that makes `LOAD_FACTOR` a global unknown
and enters `u − ū = 0` as an extra equation, so the driven DOF stays *free*.
Cross-link: [[REF_ladruno_vs_kratos_assessment]],
[[40a_kratos_crosspollination_amendment]].

**LS-DYNA implicit** treats prescribed motion as a **constraint** whose
linearized violation enters the assembled RHS (Theory R16 §38.2.2, §38.3.1,
eqs. 38.43-38.47), updating as `u + s·Δu` with `s ∈ (0,1]` from line search;
§38.3.6 warns that *"only partial line searches (s < 1) will not satisfy simple
motion constraints… convergence is prevented due to unfulfilled boundary
conditions until all prescribed motion is satisfied to within 1 %"* — a driven
face is a *documented* LS-DYNA difficulty, a different one from ours. Its
defences: BFGS with a Matthies-Strang condition-number test that can **skip** a
bad secant update (§38.3.3); reform on **residual increase** (`DIVERG=1`) and on
an **"energy explosion"** test (§38.3.4); `ILIMIT=1` for full Newton on *"highly
nonlinear problems"*; and **`MAT_24`/`MAT_123` reject the current iterate when
the plastic strain increment is too big** and retry the step, printing
`Material model rejected current iterate` to `d3hsp` (Implicit Analysis Guide
§11.3.9; exposed to user materials as the `reject` flag). ⚠ LS-DYNA has **no
bundled manual** in this environment — those citations come from the public R16
Theory / Keyword manuals and the 2025 R2 Implicit Analysis Guide, not from a
source-verified read.

### Two findings that cut AGAINST the obvious reading

- **Abaqus rates plastic active-set switching as *"usually less severe"* than
  contact** (AUG Vol 2 §7.2.3 p. 7.2.3-1) and gives it no severe-discontinuity
  handling, while it *does* carry a cutback specifically for *"very large plastic
  strain increments"* (§7.2.3 p. 7.2.3-8). So the **overstrain / collapsed-
  tangent** half of the mechanism is the better-supported half; attributing a
  ×28.6 penalty mainly to active-set chatter is a stronger claim than Abaqus
  doctrine predicts. §5 gate 3 discriminates.
- **Abaqus documents the OPPOSITE preference** for submodeling: driving a plastic
  region by **tractions** is the fragile option and by **prescribed
  displacement** the robust one (AUG Vol 2 §10.2.3 p. 10.2.3-6). Do **not**
  generalise our result to "Neumann is intrinsically safer for path-dependent
  materials" — it is safer for this *conditioning* problem and worse for
  uniqueness near a limit load.
- **No code recommends keeping load-introduction regions elastic.** Searched the
  Abaqus manuals specifically: that advice is not there. It is sound practice and
  the current project workaround, but it is a substitute for a missing predictor,
  not an idiom to cite.

---

## 4. Strengthening candidates

**We do NOT need a new `sp` constraint handler.** The SP handling is
self-consistent and gives the right answer; what is missing is a **static
predictor**. OpenSees has predictors in its transient integrators and none in
`LoadControl`.

### Candidate A (primary) — `LadrunoLoadControl` with `-extrapolate <frac>`

A strict superset of stock `LoadControl`, mirroring how `LadrunoArcLength`
relates to stock `ArcLength`. **`-extrapolate 0.0` (default) must be
bit-identical to stock.** In `newStep()`, before `applyLoadDomain`:

1. accumulate `Δu_prev` (the integrator already sees every `dU` in `update()`);
2. `theModel->incrDisp(frac · (d/d_prev) · Δu_prev)`;
3. then `applyLoadDomain(λ)` — `enforceSPs` now pre-updates a layer whose
   neighbours have already moved, so its strain increment is the physical one;
4. then `updateDomain()`.

Guards, from Abaqus's documented behaviour: **off on the first increment** of the
analysis; skip when `d/d_prev` falls outside ~[0.5, 2] (Abaqus abandons
extrapolation when the increment changes sharply — this matters for
caller-driven adaptive marches, which resize often, and it should skip right
after a cutback). It can only move the starting guess, so it **cannot change a
converged answer**.

Cost: one `incrDisp` + one extra `updateDomain` per increment. Risk: low —
the ordering is copied from `DisplacementControl`, already in the tree.
Crosswalk to document: **`EXTRAPOLATION=LINEAR` → `-extrapolate 1.0`**, in the
same style as `*STATIC, RIKS` → `LadrunoArcLength` and `STABILIZE` →
`-stabilize`.

### Candidate B (independent, a real bug) — `AutoConstraintHandler` missing `updateElement()`

Three lines, no design decision: add the loop that
`TransformationConstraintHandler::enforceSPs` has at `:518-521`. Needs a minimal
reproducer first (a single brick with a prescribed face displacement under
`constraints Auto` vs `Transformation`), then the fix, then a regression test in
the style of `tests/test_sp_subtract_init.py`. Worth offering upstream.

### Candidate C (optional) — implement `TransformationFE::getTangForce()`

Enables the Kratos-style `b −= K·Δu_D` route: form the element tangent at the
**committed** state and multiply by the local prescribed increment, so the
prescribed motion reaches the RHS **without** a constitutive evaluation in the
lagging state. Principled and exact. **But** it changes behaviour for every model
with a non-homogeneous `sp` — high blast radius for the same benefit A gets
cheaply. Only worth it if A proves insufficient.

### Candidate D (cheap, diagnostic) — make this class of problem visible

The whole episode cost hours because the pathology is invisible: committed state
is clean, the answer is right, and only the iteration count betrays it. Worth
considering: an opt-in per-**iterate** dump (residual norm, `dU` norm, count of
Gauss points that yielded in this iterate but not in the committed state), and/or
an LS-DYNA-style *"material rejected this iterate"* hook on the plastic
kernels — LS-DYNA's `d3hsp` message would have diagnosed this in one run.

### Explicitly NOT recommended

- **A new SP constraint handler** — the handler is not the defect.
- **`constraints Penalty` / `Lagrange` as a workaround.** Both *do* impose
  through the RHS (`PenaltySP_FE::getResidual` = `alpha·(constraint − (u − u_init))`;
  `LagrangeSP_FE::getResidual` likewise) — but `PenaltyConstraintHandler::handle`
  creates a `PenaltySP_FE` for **every** SP from `getDomainAndLoadPatternSPs()`
  (`PenaltyConstraintHandler.cpp:204-214`), so the fully-fixed face becomes
  penalty-enforced too; and `Lagrange` adds one equation per constrained DOF with
  zero diagonal blocks, making a large 3D solid system indefinite.
- **Adding `updateDomain()` to `LoadControl::newStep`.** *Negative* value: the
  driven layer is already updated by `enforceSPs`, and the interior has not moved
  since the last commit, so re-calling `setTrialStrain` at unchanged strain is
  idempotent. It buys one extra full constitutive sweep per increment (and per
  cutback attempt) for zero behaviour change.
- **`LadrunoArcLength` / `LadrunoIndirectControl` / `LadrunoDynamicRelaxation` /
  `LadrunoStabilizedUnbalance`** — none addresses this. The elastic-cover run has
  **0 cutbacks and no limit point**, so there is nothing for arc-length to do.

---

## 5. Validation gates — run these BEFORE writing the C++

All cheap (hex8-class, ~22.6 k dof), all with the algorithm **pinned** so the
march is deterministic. Save JSON artifacts: the precursor diagnostic
(`diag_convrate.py`) wrote to stdout only, which made its central figures
uncitable — do not repeat that.

1. **Overstrain should scale as Δδ/h.** Vary only the driven cover's
   through-thickness element size, covers plastic. The penalty should fall
   roughly linearly with coarsening and grow with refinement. Report an
   elastic-cover control at each size to separate mesh effects from the penalty.
   This is the mechanism's sharpest prediction.
2. **The penalty should scale with step size.** Covers-plastic at 10 / 20 / 40
   steps, measured per unit of **progress**, not per successful increment.
   ⚠ A prior observation cuts against this ("halving the step made convergence
   worse, 19 → 27 iterations") — but it was taken on an adaptive ladder under
   `Newton -initial` and is a heuristic, not a theorem. Settle it pinned.
   See also [[_adr76_issue_report_newton_initial]].
3. **Which half of the mechanism dominates.** Re-run the covers-plastic
   displacement case pinned to `Newton -initial` instead of `KrylovNewton`: that
   substitutes the elastic tangent for the collapsed plastic one while leaving
   the spurious residual intact. If it recovers most of the ×28.6 the collapsed
   tangent dominates; if not, the overstrained residual does — and Abaqus's
   documentation predicts the latter. **This also bounds how much a predictor can
   buy**, since a predictor attacks the residual half directly.
4. **Reproduce the `AutoConstraintHandler` defect** (candidate B) before fixing
   it. Currently a source-level inference only.

### Acceptance gates for candidate A

- `-extrapolate 0` reproduces the recorded march increment for increment
  (δ = 0.0312 / 0.0625 / 0.0938 / 0.1438 / 0.1938 / 0.2438 / 0.3238 / 0.4038 /
  0.4838 / 0.6118 / 0.7398 / 0.8678 / 1.0726 mm; step growth
  0.0312 → 0.05 → 0.08 → 0.128 → 0.2048; 70 increments, 0 cutbacks).
- `-extrapolate 1.0` collapses the covers-plastic displacement case from
  324 inc / 52 cutbacks / 4 255 iterations toward the elastic 30 / 0 / 149.
  **Anything less than a large improvement means the mechanism is wrong or
  incomplete — report that, do not tune toward the preferred conclusion.**
- The converged answer does not move (920.0 kN at the stated internal rotation
  on the coarse mesh; 1 219.6 kN on the rated mesh).

---

## 6. Open questions

> [!question]
> Does the predictor interact badly with a **caller-driven adaptive march**? The
> `d/d_prev` guard is copied from Abaqus's abandonment criterion, but Abaqus
> owns its own stepper. Ours does not.

> [!question]
> Does the predictor help, hurt or do nothing under **contact** (`LadrunoContact`)
> and under `-geom corot`/`finite`? Both add their own active-set problems on top
> of this one, and both are the next rungs of the campaign that found this. The
> Abaqus caution about extrapolation and "unloading as a result of cracking" is
> the warning shot.

> [!question]
> A precursor measurement recorded an iteration-1 residual of **895 239 N** with
> the covers plastic against **3.9e-9 N** with them elastic. On the corrected
> sequence (§2 step 3-4) a **large** iteration-1 residual is expected for *both*
> materials, because the driven layer is pre-updated either way — so the
> **3.9e-9 figure is the one that now looks wrong**. Both are unreproducible as
> recorded. Regenerate with saved artifacts before either is cited.

> [!question]
> Should the plastic kernels gain an LS-DYNA-style `reject this iterate` hook
> (candidate D)? It is a bigger change than the predictor but it makes the whole
> failure class self-diagnosing, and it composes with the caller's step control.
