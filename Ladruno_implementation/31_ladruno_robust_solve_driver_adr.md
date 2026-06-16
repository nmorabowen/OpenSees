---
title: Robust nonlinear-solution driver — automatic incrementation + stabilization
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - adr
  - solver
  - integrator
  - convergence
  - stabilization
  - observability
  - research
---

# Robust nonlinear-solution driver — automatic incrementation + stabilization

> Siblings: [[20_ladruno_arclength_stabilized_adr]] (the `-stabilize` viscous seam this driver sequences), [[21_ladruno_dynamic_relaxation_adr]] (the matrix-free dynamics fall-through), [[22_ladruno_dissipation_arclength_adr]] (the scoped-but-unbuilt snap-back keystone, classTag 33007), and [[26_ladruno_plane_frontier_adr]] (the regularized materials whose softening this driver is built to push through). This is a **decision document**: it lands on ONE architecture (a Python incrementation spine + named C++ grafts), states what ships in v1, and hard-defers the rest with reasons.

## What

A "robust nonlinear-solution driver" — an Abaqus automatic-incrementation + `*STATIC,STABILIZE` analog, native to OpenSees, that **orchestrates primitives the fork already ships** to solve hard nonlinear FE problems (softening concrete, buckling, snap-through, slack cables, contact) without operator babysitting and without silently corrupting the physical answer.

The **spine** is automatic increment sizing in Python (`ladruno_drive`, an extension of the already-shipped [`ladruno_solve.py`](../Ladruno_scripts/ladruno_solve.py)). The **robustness payload** is the already-shipped `LadrunoArcLength -stabilize` viscous regularization (ADR-20), tested on TRUE equilibrium via the already-shipped `LadrunoStabilizedUnbalance`. The **floor** is the already-shipped matrix-free dynamics (`LadrunoDynamicRelaxation` / `ExplicitBatheLNVD`).

**In scope (v1):** the Python spine with Ramm target-iteration sizing; a JSONL decision log; a STABILIZE reflex; the dynamics fall-through; and **two narrow C++ observability queries** that are on the correctness critical path (without them the dissipation hard-gate is a no-op). **Out of scope / hard-deferred:** the deep-C++ in-loop `RobustNewton`, the b24ac router, a recipe-pack DSL, a Tcl port, and `LadrunoDissipationArcLength` (33007) — each gets its own follow-up ADR gated on a real production model that needs it (§6, §10).

## Why

A multi-year **production research tool** lives or dies on two properties this driver must guarantee simultaneously, and which the naive "algorithm-zoo escalation" straw-man does NOT:

1. **It must not fail to converge** on the hard physics — softening, snap-through, contact chatter, near-singular tangents.
2. **It must not silently change the physical answer** to buy that convergence. A run that "converged" on a heavily-diffused, step-size-dependent path and then feeds a published load-displacement curve is the worst failure mode a research tool has, because it is invisible.

The fork already shipped the hard numerics for (1): the increment controller (`ladruno_solve.py:adaptive_static`, proven against the Gate-2 confined column where naive Newton diverges, [[21_rc3d_validation_gates]]), the viscous STABILIZE seam (`LadrunoArcLength -stabilize`, classTag 33004), the true-equilibrium test (`LadrunoStabilizedUnbalance`), and the matrix-free dynamics fall-through (`LadrunoDynamicRelaxation` 33005, `ExplicitBatheLNVD` 33002). **What is missing is the orchestration that unifies them behind one entry point, the observability that makes (2) auditable, and the honest verification passes that bound artificial diffusion.** This ADR is therefore mostly glue + observability + discipline, NOT new numerics — which is exactly the right risk profile for a tool whose results feed papers.

## Context — what the fork already has, verified against source

The standard OpenSees idiom (swap algorithm/integrator/test/constraints between `analyze()` calls) is legal because the engine reverts on failure. **Verified:**

- **Revert-on-fail is a real contract but a per-class one, NOT a framework guarantee.** `StaticAnalysis::analyze` calls `the_Domain->revertToLastCommit()` + `theIntegrator->revertToLastStep()` on the newStep/solve/commit failure paths (`StaticAnalysis.cpp:171-172, 185-186, 220-221`). `Domain::revertToLastCommit` is a pure **delegator** — it loops `node->revertToLastCommit()` / `ele->revertToLastCommit()` and restores nothing itself. It is only as correct as each material's hand-coded revert. **The headline core (ForceBeamColumn3d, DispBeamColumn3d, ASDConcrete3D) honors it; `PressureDependMultiYield::revertToLastCommit` is a literal `return 0;` no-op (PDMY.cpp:868-872).**
- **`LadrunoArcLength -stabilize`** (33004): viscous-regularized incremental load control; `formTangent` adds `cOverDt·M*`, `formUnbalance` subtracts `f_v=(c/Δt)M*·Δu`; `c` auto-calibrated to a dissipated-energy fraction (default `f=2e-4`, the Abaqus default); **mutually exclusive with the arc-length quadratic** (it IS load control in this mode); records the true pre-`f_v` residual. `M*` is an integrator-owned unity/gershgorin/lumped mass — NOT `Element::getMass()` (zero on zero-density research models, the ADR-20 blocker).
- **`LadrunoStabilizedUnbalance`**: norms `‖λp−f_int‖` WITHOUT `f_v` via the integrator's `getStabilizedTrueResidual()` (LadrunoArcLength.h:135), degrading to SOE `‖B‖` otherwise. A converged stabilized step is a genuine equilibrium point, not merely a regularized one.
- **`LadrunoIndirectControl`** (33006): CMOD / weighted multi-DOF control, monotone through snap-back; needs the control DOF a-priori.
- **`LadrunoDynamicRelaxation`** (33005) + [`relax_to_static.py`](../Ladruno_scripts/relax_to_static.py): matrix-free quasi-static fall-through. **`kineticEnergy` and `residualNorm` are PRIVATE (DR.h:128-129) with no `getResponse` channel** — the rest-detection script today infers rest from `nodeVel`/`nodeDisp` decay.
- **`ExplicitBatheLNVD`** (33002) / `CentralDifferenceLadruno` (33003) + `criticalTimeStep()`: explicit fall-through with queryable `dt_cr`.
- **`ladrunoArcLength <sub>` runtime command** (`OpenSeesCommands.cpp:3129-3188`): exposes `reduceStep|increaseStep|setArcLength|arcLength|deltaLambdaStep|currentLambda|sign|deltaUstepNorm`. **It does NOT expose the dissipation ratio, the true residual, `cVisc`, or a runtime `cVisc` mutator** — these are the missing observability seam (§5, R-OBS).

## Decision — the spine + named grafts

We adopt **Design 1 (Abaqus-auto-increment, faithfully) as the spine**, with four grafts (b24ac-as-a-radius-signal-only, the dynamics micro-burst classifier, the indirect-control polish tail, and stay-resident hysteresis) and one mandatory C++ observability seam. This is the design panel's unanimous ranking (8.7/8.7/8.7) and the correct production risk profile: the hard numerics already exist and are source-verified, so v1 is orchestration + observability.

**Rejected as spines** (kept only as graft sources or deferred): Design 2 (Python state machine — a less-opinionated superset of Design 1, and it misreads the `domainStamp` mechanism); Design 3 (deep-C++ `RobustNewton` — forks the most battle-tested loop in OpenSees, wrong default risk profile, harvested as a deferred enhancement); Design 4 (dynamics-first — inverts the cost economics the wrong way until selective mass scaling ships; its micro-burst and getResponse seam are grafted in).

### The spine: automatic increment sizing, not algorithm-zoo escalation

The driver owns an outer loop over `ops.analyze(1)`. The PRIMARY control is the **step size**, sized to a target iteration count (Ramm / Abaqus `I_d ≈ 6`) from `ops.testIter()`: grow on easy convergence, cut back on hard convergence. The algorithm zoo (`KrylovNewton→Newton→NewtonLineSearch→ModifiedNewton`) is demoted to a single **cheap secondary reflex** tried in-place before halving (it is free — one `ops.algorithm()` call, no rebuild, and `adaptive_static` already ships it). When cutback alone keeps failing AND the failure is a limit point, the driver arms **STABILIZE**. When even stabilized cutback stalls, it hands off to **dynamics**. This is deliberately NOT the straw-man's prominent Newton-variant escalation.

### The escalation ladder (shallow by design)

Each rung is attempted from the last **truly-equilibrium-committed** state (§Risk R-SWITCH explains why "committed" ≠ "equilibrium" and why this matters).

| Rung | Action | Where | Trigger | Cost |
|---|---|---|---|---|
| **0** | Increment cutback / Ramm grow | Python | any fail / clean success | cheap (no rebuild) |
| **0.5** | Algorithm reflex (one step up the ladder, same scale) | Python | fail before halving | cheap (`setAlgorithm` does NOT rebuild) |
| **1** | AUTO-STABILIZE (`-stabilize f -adaptStab` + `LadrunoStabilizedUnbalance`) | Python orchestrating C++ | `scale < 1/32` AND limit-point signature | medium (integrator swap → rebuild) |
| **2** | Path-follower switch (`LadrunoArcLength -adapt` for diffuse limit points; `LadrunoIndirectControl` for snap-back with known control DOF) | Python orchestrating C++ | load factor must reverse (snap-back); STABILIZE cannot follow it | medium |
| **3** | DYNAMICS fall-through (`LadrunoDynamicRelaxation` / `ExplicitBatheLNVD`, `dt` from `criticalTimeStep()`), then **indirect-control polish tail** back onto the branch | Python orchestrating C++ | stabilized Newton stalls; near-zero tangent; contact chatter | expensive (analysis-type swap) |

**Rung 3 is the guaranteed floor** (a dynamics march has no convergence test to fail), but it returns a *regularized rest state*, not a *traced equilibrium branch* — so it is never the publishable answer on its own; the polish tail (graft, from Design 4) re-lands on the true branch with a few CMOD-controlled implicit steps.

### The state machine

```
                 ┌─────────────────────────────────────────────────────────┐
                 │  for each NOMINAL increment:                             │
                 ▼                                                          │
   ┌──────────────────────────┐  analyze(1)==0 (clean)                     │
   │ RUNG0  cut/grow (Ramm)    │ ─────────────────────────► commit, grow ──┤
   └──────────┬───────────────┘                                            │
              │ fail                                                       │
              ▼                                                            │
   ┌──────────────────────────┐  next algo passes                         │
   │ RUNG0.5 algo reflex       │ ─────────────────────────► commit ───────┤
   └──────────┬───────────────┘                                            │
              │ ladder exhausted at scale                                  │
              ▼                                                            │
   ┌──────────────────────────┐  CLASSIFY failure (see below)             │
   │  micro-burst probe        │                                          │
   └───┬──────────┬───────┬───┘                                           │
       │ limit pt │ bad   │ genuine instability                          │
       │          │ step  │ (KE grows & holds)                           │
       ▼          ▼       ▼                                              │
   RUNG1       (back to  RUNG3 dynamics ── rest ── polish tail ──────────┤
   STABILIZE    RUNG0)   │  (or ABORT w/ diagnostic if mechanism)        │
       │                 │                                               │
       │ snap-back       │                                               │
       ▼ detected        │                                               │
   RUNG2 path-follower ──┴───────────────────────────────────────────────┘
       │
       └─► on STABILIZE success: HARD-GATE on dissipation ratio + true residual;
           on phase end: AUTOMATIC c-reduction verification pass (peak-load drift)
```

**Failure classification** (the hardest problem, and the spine's weakest link without help). The driver routes with a layered detector, honestly labelled heuristic until the C++ signals ship:

- **Cheap first pass:** the failure return CODE (`-2` newStep / `-3` solve / `-4` commit) plus the snapshotted `ops.testNorms()` / `ops.testIter()` *of the failed step, captured BEFORE any retry* (the buffer is zeroed by the next `start()`). A `-3`/NaN means the norm buffer is untrustworthy — route on the code, not the trajectory.
- **Tiebreaker (graft, Design 4): the dynamics MICRO-BURST.** From the last committed state, run a short, **dedicated, undamped, fixed-`M*`** DR/explicit burst and read the kinetic-energy *growth rate*: KE-spikes-and-holds ⇒ genuine instability (go dynamic or abort); KE-rings-down-to-rest ⇒ it was just a bad step and the rest point IS the answer. This is the only principled discriminator that needs no eigen-probe, but it must run on a throwaway configuration (NOT the production Cundall-damped DR, which zeroes KE at every peak and would blind the signal — see R-MICROBURST).
- **Hard guard (mandatory):** an **indefinite tangent is NOT, by itself, grounds for dynamics escalation** — every softening point has one. A limit point ⇒ *switch parameterization* (rung 1/2); only growing dynamic energy ⇒ dynamics. Cap **cumulative consecutive STABILIZE-invoked increments without residual recovery**; abort with a "suspected mechanism / missing restraint" diagnostic rather than diffusing a genuine instability into a plausible-but-wrong answer.

## Layered build plan

The honest layering inverts the straw-man (whose "Layer 1 = build C++ auto-stabilization" is already shipped). **v1 is one Python module + two narrow C++ queries.**

### Layer 0 — `ladruno_drive` (Python, the spine, the bulk of v1)

Extends `Ladruno_scripts/ladruno_solve.py`. Adds: (a) Ramm target-iteration sizing of grow/cutback from `ops.testIter()` (today `adaptive_static` uses a fixed `grow=2` and does NOT read `testIter`); (b) the rung-1 STABILIZE reflex and rung-3 dynamics handoff behind one `solve(...)` entry point; (c) the structured **JSONL decision-log writer** (today the only sink is an `on_step(progress)` float — there is no log); (d) the failure classifier + micro-burst probe; (e) **stay-resident hysteresis** (graft, Design 2) on the EXPENSIVE rungs only. No C++. No ledger/banner/classTag rows (it is not a C++ feature).

> **Verified cost-model nuance the grafts must respect:** only `setIntegrator`/`setLinearSOE`/`setNumberer` set `domainStamp=0` and force a full `domainChanged()` rebuild (`StaticAnalysis.cpp:525,558,479`). `setAlgorithm` (line 503) and `setConvergenceTest` (591) do NOT. So the algorithm/test reflex (rung 0.5) is genuinely free; hysteresis must gate the integrator/constraint rungs (1/2/3), not the cheap reflex.

### Layer 1 — STABILIZE (already shipped C++; Python orchestrates)

`LadrunoArcLength -stabilize f -adaptStab` (33004) + `LadrunoStabilizedUnbalance`. **No new integrator.** The ONLY new C++ is the observability seam (Layer 1.5), which is on the correctness critical path.

### Layer 1.5 — the mandatory observability seam (two narrow C++ queries)

These are NOT optional conveniences — the dissipation hard-gate and the DR quasi-staticness certification cannot exist without them (R-OBS, R-DR-ENERGY). Add them as **named getters on the two concrete classes that need them**, each with its own single `OPS_` dispatch (mirroring how `getCriticalTimeStep` was added — there is NO generic `getResponse` channel on the integrator base class, and we do NOT add one in v1):

1. `LadrunoArcLength::getStabilizationDissipationRatio()` returning `dissipVisc / E_strain` (today `dissipVisc` is private, h:177; no accessor), plus surface `getStabilizedTrueResidual()` and a runtime **`cVisc` mutator** (`scaleCVisc(f)`) for the ramp-down (R-RAMPDOWN). Wired into the existing `ladrunoArcLength <sub>` command as new subcommands `dissipRatio | trueResidual | scaleCVisc`.
2. `LadrunoDynamicRelaxation` named getters for `kineticEnergy`, `residualNorm`, and a force-based `restMetric = residualNorm/‖P‖` (the CORRECT, mass-free quasi-staticness metric — see R-DR-ENERGY), via a single `OPS_` dispatch.

### Layer 2 — primitives (mostly shipped; one deferred keystone)

Dynamics fall-through (33005/33002/33003 + `criticalTimeStep()`) is shipped; Layer 0 orchestrates `relax_to_static.py`. Mass-scaling for explicit affordability on production 3D meshes (selective mass scaling) is **roadmapped, not built** — so the dynamics rung is honestly bounded to coarse/uniform meshes until it lands (R-SMS). Regularization-on-by-default is a *modeling-template* decision ([[26_ladruno_plane_frontier_adr]]), NOT a driver concern — explicitly removed from scope.

## API sketch

### Python (openseespy / local build)

```python
from ladruno_drive import solve

# nominal integrator is the user's choice; the driver scales it via a callback
ops.integrator("DisplacementControl", ctrl, dof, du)

result = solve(
    ops,
    set_step   = lambda s: ops.integrator("DisplacementControl", ctrl, dof, du*s),
    n_steps    = 400,
    target_iter= 6,                 # Ramm I_d
    stabilize  = "auto",            # off | auto | on ; "auto" = rung-1 reflex
    f_stab     = 2e-4,              # Abaqus dissipated-energy fraction default
    diss_gate  = 5e-3,             # HARD reject if dissipVisc/E_strain exceeds
    verify     = True,             # MANDATORY c-reduction pass on stabilized phases
    snap_back_control = (node, dof),# arms rung-2 IndirectControl if provided
    dynamics_fallback = "DR",       # DR | LNVD | off
    polish     = True,              # indirect-control tail after a dynamics excursion
    log        = "run.jsonl",       # structured decision log (required for provenance)
)
# result.verdict ∈ {"equilibrium", "regularized", "unverified", "aborted"}
# result.peak_load_drift, result.max_diss_ratio, result.rung_histogram
```

### Tcl

v1 ships **no Tcl port** (R-TCL). The Tcl path retains the already-shipped `ladrunoArcLength reduceStep` Layer-B cut-and-retry idiom. A future C++ lift of the increment-sizing spine into a `StaticAnalysis` flag (ADR-20 §8 follow-up #5) is the tracked route to dual-interpreter parity; it is not built speculatively.

## Alternatives considered

The design panel produced five candidates; three judges independently ranked Design 1 first (8.7).

| Design | Philosophy | Verdict |
|---|---|---|
| **1 — Abaqus-auto-increment** | Increment sizing is the spine; STABILIZE is the payload (on load control, never the follower); algorithm zoo demoted to a reflex. | **CHOSEN SPINE.** Lowest implementation risk (orchestration over source-verified shipped C++); only design that makes step-size-dependence of stabilized results a mandatory automatic check, not a footnote. |
| **2 — Python state machine** | Full state-machine ladder, Python-first, minimal audited C++. | **Rejected as spine, grafted.** A less-opinionated superset of Design 1 without its STABILIZE/verification discipline; **misreads `domainStamp`** ("each setter sets `domainStamp=0`" is false for `setAlgorithm`/`setConvergenceTest`). **Grafted:** stay-resident hysteresis; decision-log-as-recorder framing; the precise location of the observability wall. |
| **3 — Deep-C++ `RobustNewton`** | A new `EquiSolnAlgo` leaf that escalates INSIDE the corrector loop (in-loop tangent ladder, line search, stabilize-arming, `getDivergenceMode()`). | **Rejected as spine, deferred.** Real insight (the per-iteration residual trajectory the engine discards IS only computable in-loop), but forks the 30-line `NewtonRaphson::solveCurrentStep` — a silent drift bug there corrupts EVERY analysis; mid-step K/B mutation breaks the convergence test's monotone-descent assumption; parallel correctness unaddressed. **Harvested:** the divergence-mode *idea*, pulled UP into Python via the micro-burst rather than down into a forked loop. If ever built: a NEW leaf (never edit `NewtonRaphson`) with a gates-off bit-identical parity test. |
| **4 — Dynamics-first** | Dynamics is rung 0 (the guarantee); implicit is an opt-in speed bet that falls through. | **Rejected as spine, grafted.** Inverts the cost economics the wrong way for production until selective mass scaling ships (its own risk section concedes implicit-first-for-speed should be the practical default — i.e. it IS Design 1 with a better fallback). **Grafted:** the dynamics micro-burst classifier; the `getResponse` observability seam; "dynamics-as-guaranteed-fall-through"; the indirect-control polish tail; energy-ratio as a first-class control input. |
| **5 — Continuation-arclength spine** | Riks-style predictor-corrector continuation; b24ac<0 as a rigorous limit-point trigger; load control as degenerate continuation. | **Rejected as spine, partially grafted.** Best mechanics rigor and best path-fidelity, but structurally single-parameter / proportional-loading; its keystone (control-location-free `LadrunoDissipationArcLength` 33007) is scoped-but-unbuilt; honest that `-stabilize` PAUSES path-following (mutual exclusivity). **Grafted with a correction:** b24ac<0 is used ONLY as a *radius-too-large* signal within pure (non-stabilized) arc-length — NOT as a cross-mode instability oracle, because it is dead code in `-stabilize` mode and its own error string proves it conflates step-too-big with limit point. |

## Energy accounting — gate vs. audit (settled; scopes R-OBS / R-DR-ENERGY)

The dissipation signal has **two consumers needing two mechanisms**; conflating them is the root of R-OBS and R-DR-ENERGY.

**Why the source is the integrator, not the recorder.** `f_v = (c/Δt)·M*·Δu` is built inside `LadrunoArcLength::formUnbalance`; its inputs `c`, `M*` (the unity/Gershgorin pseudo-mass) and `Δu` are integrator-private, and `M*` is **never assembled onto nodes/elements**. `EnergyBalanceKernel` only reads `ele/node->getMass/getDamp/getResistingForce`, so it is structurally blind to `f_v`/`M*` and *cannot* reconstruct the artificial work from the domain. Any recorder-side accounting must query the integrator anyway — the getter is unavoidable. (The recorder would own this natively only if stabilization were physical dashpot elements, which [[20_ladruno_arclength_stabilized_adr]] deliberately did not do.)

**Gate (control) → integrator getter.** Layer-1.5 on `LadrunoArcLength` (33004): `getStabilizationDissipatedEnergy()` (cumulative `W_stab=Σ f_vᵀΔu`), `getReferenceStrainEnergy()`, `getStabilizationDissipationRatio()`, and the `scaleCVisc()` actuator, via the `ladrunoArcLength` subcommand. Read synchronously per step, never through recorder IO. DR's settling signal is the analogous mass-free force-residual getter on `LadrunoDynamicRelaxation` (R-DR-ENERGY).

**Audit (trust) → recorder, fed by the integrator.** Key identity: at a stabilized equilibrium `P_ext ≈ F_int + f_v`, so the per-step power balance is `ULW = IE + W_stab + KE (+ DW_phys)`. The recorder's existing residual is therefore **already** `RES = ULW − (KE+IE+DW) = W_stab` — the artificial dissipation sits in `RES`, mislabeled as numerical error. The fix is a **partition, not a new measurement**:

```
SW       := W_stab            (the ONE scalar from the integrator getter — now labeled)
RES_true := RES_old − SW      (genuine numerical non-closure)
```

`SW/E_ref` is then the fidelity/trust metric and `RES_true` the real balance check. Mechanically: `ebkernel::EnergyAccumulator::step()` gains a `stab_work` input; `out[]` widens 6→7 (`KE,IE,DW,SW,ULW,RES,ERR`, `RES` subtracting `SW`); `SW` is the integrator's **cumulative** value (cadence-independent), not re-integrated.

**Three levels (ship L0 in v1):**
- **L0 — no recorder code change:** integrator exposes `W_stab` as a queryable response; record/log it; partition `RES → SW + RES_true` in post-processing. Correct and free given the gate getter.
- **L1 — recorder-native `SW`:** `EnergyBalanceRecorder`/`EnergyBalanceSource` carry `SW` as the 7th component. Wrinkle: the recorder attaches to the `Domain` but the integrator lives in the `Analysis`, so L1 needs a fork hook to hand the recorder the active stabilized-integrator handle (setup wire / `ladruno::activeStabilizer()` registry, serial-only v1). Pure ergonomics (one self-closing balance) — deliberate, not rushed.
- **L2 — per-region `SW`:** needs per-DOF stabilization power to apportion the global `f_v` to `MeshRegion`s. v2.

Net: R-OBS/R-DR-ENERGY are **scoped, not fatal** — the gate is the integrator getter (L0, blocking v1), the audit partition is post-processing (L0, free), and the kernel's physical-mass KE needs **no** change (it is correct for its real job: the true-dynamics audit).

## Risks & must-fixes

Every surviving red-team attack, folded in with severity and mitigation. Severity: **FATAL** = ships a no-op correctness gate or a silently-wrong result unless fixed; **MAJOR** = correctness hazard on the hard paths the driver exists for; **MINOR** = artifact / cost / discipline.

### R-OBS — the dissipation hard-gate is unbuildable today **[FATAL]**

The mandated "on-by-default dissipation hard-gate" reads `dissipVisc/E_strain`, but `dissipVisc` is **private** (LadrunoArcLength.h:177) with no accessor, the `ladrunoArcLength` command exposes no dissipation query (verified `OpenSeesCommands.cpp:3180-3188`), and `getStabilizedTrueResidual` is C++-only. The gate, as described by the panel, silently no-ops.
**Mitigation:** ship the Layer-1.5 seam as a BLOCKING v1 deliverable. **Until it lands, the driver must REFUSE to enter stabilized phases (fail loud, not diffuse silent).** Do NOT describe any hard-gate as "on by default" before the live query exists.

### R-DR-ENERGY — the quasi-staticness gate is computed against the wrong (often zero) mass **[FATAL]**

The `EnergyBalanceRecorder` KE is `½vᵀMv` with the **physical** mass (`EnergyBalanceKernel.h:87,136` use `ele->getMass()`/`node->getMass()`), which is **~0 on exactly the zero-density research models DR exists to rescue.** So an "energy ratio < 1e-4 ⇒ quasi-static" gate built on the recorder ALWAYS passes for the DR phase, certifying an inertia-laden regularized result as static.
**Mitigation:** the DR-phase quasi-staticness criterion is the **force-based** `residualNorm/‖P‖` (mass-free), read from the DR-internal getter (R-OBS Layer-1.5 query #2), NOT any physical-mass KE ratio. Two distinct watchdogs: physical-mass `EnergyBalance` for implicit/stabilized phases; DR-internal force residual for the DR phase.

### R-REVERT-PDMY — `revertToLastCommit` is a per-class contract, not a guarantee **[MAJOR]**

`Domain::revertToLastCommit` delegates and restores nothing itself; `PressureDependMultiYield::revertToLastCommit` is a literal `return 0;` (PDMY.cpp:868-872). The spine's "retry on already-reverted clean state" premise is false at the material level for the PDMY/PIMY soil family (saved today only by the accident that the total-strain `setTrialStrain` path re-derives trial from committed before it is read).
**Mitigation:** (a) scope any soil/SSI recipe to total-strain (displacement-based) elements; (b) ship a **revert-parity + revert-idempotency** regression battery (drive to divergence mid-step, revert, assert `getStress`/`getTangent` bit-match a clean re-run; call revert twice, assert the second is a no-op); (c) record in [[LEDGER_quirks]] that revert is a hand-coded per-class contract; (d) audit every fork-authored material for revert symmetry. **Reframe the ADR's prose away from "reverts every material point's history" — that is false.**

### R-REVERT-RC — a failed revert is invisible to the controller **[MAJOR]**

`StaticAnalysis` discards `revertToLastCommit()`'s return code on every failure path, and the analysisStep path (line 142) calls `revertToLastCommit()` but **not** `theIntegrator->revertToLastStep()` (unlike lines 171-172/185-186/220-221). So a revert that fails to restore, or a `-2`-from-analysisStep that leaves the arc-length per-step snapshot polluted, is invisible to Python.
**Mitigation:** (a) the revert-contract regression test must assert the RETURN CODE of an explicit domain revert, not just that `analyze()` returns negative; (b) after any failed `analyze()`, the driver probes a cheap committed-state invariant (control-DOF `nodeDisp` equals the last value it logged) and **ABORTS** rather than retries if the domain did not roll back; (c) ledgered candidate vanilla edit: add `theIntegrator->revertToLastStep()` to `StaticAnalysis.cpp:142` for parity (marked `// Ladruno`, recorded in [[LEDGER_vanilla_files]]).

### R-SWITCH-PHAT — integrator switch corrupts the reference load unless from TRUE equilibrium **[MAJOR]**

`ArcLength`/`DisplacementControl`/`LadrunoArcLength::domainChanged` all probe `phat` via `currentLambda+=1; formUnbalance(); phat=getB()`, carrying the source comment **"assumes unbalance at last was 0"**. Switching OUT of a `-stabilize` step (the committed residual is `f_v`-polluted) or a loose-tolerance step makes `phat` absorb the residual error — every subsequent arc-length/displacement step is scaled against a corrupted reference load (a SILENT wrong answer on the softening branch). Compounds with R-DR-ENERGY: "committed" ≠ "true equilibrium".
**Mitigation:** only switch integrators from a step that converged on the TRUE static residual (`NormUnbalance`/`LadrunoStabilizedUnbalance`); before any swap, assert `‖λp−f_int‖ < tol`. AL-9-style restored-state regression test: switch mid-softening, assert `phat` equals the pristine reference load to round-off.

### R-SWITCH-SIGN — switching into arc-length predicts the wrong branch direction **[MAJOR]**

`signLastDeltaLambdaStep` is constructor-defaulted to `+1` and is NOT re-derived from committed history in `domainChanged` (LadrunoArcLength.cpp:531-533). A LoadControl→ArcLength switch onto a DESCENDING (softening) branch forces a `+λ` predictor — it climbs back up the branch it just came down, immediately hitting `b24ac<0` or re-diverging. Likewise the arc *radius* is discontinuous across the switch (a nominal radius after a crawling cutback phase overshoots the limit point).
**Mitigation:** seed the predictor sign from the committed path tangent (the sign of the last load-factor increment Python observed before the switch) and seed `arcLength2` from the outgoing `‖deltaUstep‖`, not a recipe constant; take the first post-switch step at reduced scale; cover with a snap-through regression that switches on the descending branch and asserts `λ` decreases on step one. (Needs a small C++ setter to inject the sign, since none exists.)

### R-DIFFUSION — stabilized solutions are step-size / c-dependent **[MAJOR]**

A "converged" STABILIZE run is NOT a step-converged equilibrium path; refining the stepping changes the answer. A single run's log cannot contain this — it is revealed only by comparing two runs at different `c`/step.
**Mitigation:** the automatic **c-reduction verification pass** (re-run the stabilized SEGMENT at halved `c`, bound the **peak-load drift**) is MANDATORY and on-by-default, writing the drift as a provenance stamp. Make it cheap by reusing the committed state and re-running only stabilized segments. The decision log's provenance field is REQUIRED: if the drift was not computed, the run is stamped **`verdict="unverified"`** and that label propagates to result files — so "I skipped the bound" is visible, not invisible. `LadrunoStabilizedUnbalance` proves point-wise equilibrium but is NOT a fidelity certificate; the c-reduction drift is the only trustworthy one.

### R-RAMPDOWN — no ramp-down; STABILIZE biases the whole post-peak branch **[MAJOR]**

`f_v` is injected every `formUnbalance` while stabilize is on, with no mechanism to detect "limit point cleared" and decay `c`. A run that needed one stabilized increment carries the dashpot for the entire remaining softening path, and `E_strain0` is frozen at the first (post-cutback, near-singular, worst-baseline) stabilized commit so the "energy fraction f" guarantee is meaningless on exactly the problems it is invoked for.
**Mitigation:** (a) once N consecutive increments converge with `‖f_v‖` below a tight gate, ramp `c` toward zero via the new `scaleCVisc` runtime mutator (R-OBS) and attempt de-stabilized re-issue; (b) use a windowed/incremental dissipation ratio for `-adaptStab`, not cumulative-over-frozen-baseline; clamp `cVisc` to `[cMin,cMax]`; (c) prefer an explicit `-cVisc` calibrated from a dedicated ELASTIC probe step at the start of the analysis for pre-loaded models.

### R-SNAPBACK — STABILIZE silently mislabels snap-back as a valid curve **[MAJOR]**

In `-stabilize` mode `LadrunoArcLength` is load control (mutually exclusive with the quadratic, h:167-168); when `λ` reverses it produces a regularized-but-WRONG **monotone** path, and `LadrunoStabilizedUnbalance` reports each point as converged true equilibrium. The user gets a smooth load-displacement curve and may publish a wrong branch. The control-location-free fix (`LadrunoDissipationArcLength` 33007) is scoped-but-unbuilt.
**Mitigation:** detect `λ`-reversal in Python (`currentLambda` is queryable) while in load-control/stabilize mode and RAISE a loud "snap-back encountered while path-following is disabled — segment is regularized, not traced" record; route to `LadrunoIndirectControl` with an explicit control DOF, or abort. Never let a stabilized monotone segment inherit the `verdict="equilibrium"` label. Schedule 33007 as a follow-up ADR triggered by a real snap-back model — NOT a v1 commit.

### R-INSTAB — no rigorous instability detector exists in the default solve path **[MAJOR]**

The rigorous signal (tangent inertia / negative-pivot count) is **not computed or exported by any shipped solver the spine runs**: the default BandGen LU (DGBSV) reports only exact singularity; even `ProfileSPDLinDirectSolver` factorizes through a negative pivot silently. The grafted "b24ac<0 rigorous trigger" is dead code in `-stabilize` mode and its own error string ("initial load increment was too large") proves it conflates step-too-big with limit point.
**Mitigation:** (a) DEMOTE every claim of "rigorous instability detection" to **heuristic-with-hard-abort-cap** (the micro-burst tiebreaker + the consecutive-stabilization abort guard) until an inertia query ships; (b) the real future fix is a PARDISO/MUMPS-backed `getNumNegativeEigen()` (PARDISO is in-tree and returns this for free) OR a shifted-inverse smallest-eigenvalue probe at the decision point only — both are NEW C++ on the critical path, deferred with that label; (c) distinguish **limit-point (switch parameterization)** from **instability (dynamics/abort)** as SEPARATE verdicts — an indefinite tangent alone is NOT grounds for dynamics, else the driver routes every softening model into the unaffordable dynamics rung.

### R-MICROBURST — the micro-burst is blinded by the production DR's own damping **[MAJOR]**

The grafted KE micro-burst runs on `LadrunoDynamicRelaxation`, whose default Cundall damping ZEROES all velocities at every KE peak (and `autoRefresh` rebuilds `M*` mid-burst, rescaling the KE normalization). "KE-spikes-and-holds" is exactly the signature Cundall damping is built to suppress — the two classes collapse into one. A true mechanism (missing restraint) is the case the detector must catch and the case kinetic-damped DR most easily fakes a rest for.
**Mitigation:** run the burst on a DEDICATED, undamped (or fixed-light-viscous), `autoRefresh`-OFF, fixed-`M*` throwaway integrator; measure KE GROWTH RATE over the first few steps (exponential-fit) before any peak-zeroing fires; OR use the **pre-damping** `residualNorm` growth rate (computed before Cundall zeroing) as the signal. Before accepting any DR rest as "the answer," require a force-residual check AND (when available) a tangent-singularity probe; abort on a zero-energy mode rather than committing a kinetic-damped mechanism rest.

### R-HANDOFF — implicit↔dynamics state transfer is NOT "free" for v/a and time **[MAJOR]**

"Continuity is free via the shared committed Domain" is true ONLY for displacement. (a) DR `domainChanged` re-seeds `Vhalf` from committed velocity (DR.cpp:281-282), so re-entry after a prior excursion that committed mid-chunk (non-KE-peak, nonzero v) injects a phantom velocity kick. (b) DR advances domain time as pseudo-time and ramps the load through any active TimeSeries unless `loadConst` is issued — nothing enforces it, so a "relaxation" can relax under a moving load. (c) The committed pseudo-time corrupts the next static load factor.
**Mitigation:** make `handoff_to_DR` / `handoff_from_DR` **atomic, tested primitives**: snapshot real `λ=getTime()`; `loadConst('-time', λ)`; zero v/a on entry (DR must start from rest BY CONTRACT — add a `-zeroVel` default-ON option); run DR; on return `setTime(λ)` and `loadConst` before resuming statics. DR-side guard: abort if the active load pattern is not `loadConst` at DR `newStep`. Regression: implicit→DR→return, assert `getTime()` and the next load factor are exact.

### R-THRASH / R-TERMINATION — no global work budget or termination guarantee **[MAJOR]**

Every bound is LOCAL (per-nominal-step `min_scale`; per-integrator `reduceStep` floors at `ellMin²/1e-30` and then **returns 0 unconditionally** — the documented `while ok: reduceStep` C++ idiom spins forever at the floor). Stitching local-terminating loops into a multi-phase ladder does not yield a globally terminating controller; it can livelock between rungs, re-paying O(model) `domainChanged` on every integrator/constraint switch.
**Mitigation:** (a) a single **global budget** (max-total-`analyze(1)`, max wall-clock per nominal increment) + a **monotone phase ladder** (one-way ratchet: once escalated past rung k, do not return below k within the same increment without strict progress in `λ`/control-DOF > tol); (b) route ALL cutback through the bounded Python `scale`/`min_scale`, never the raw C++ `while ok: reduceStep` (or make `reduceStep` return `<0` at the floor); (c) **stay-resident hysteresis** on integrator/constraint rungs (N floor-step failures before paying a swap); (d) gate each expensive correctness pass behind the budget (c-reduction once per SEGMENT, not per increment; micro-burst hard-capped step count; DR relax max-substep cap returning INCOMPLETE). Accept that a multi-year tool must be allowed to return "could not certify within budget" rather than spin. Regression: a pathological no-equilibrium fixture returns INCOMPLETE within a bounded substep count.

### R-LOG-MASK — the decision log transcribes the classifier's verdict, not ground truth **[MAJOR]**

A log records DECISIONS, not GROUND TRUTH: it will read "rung1 STABILIZE armed; converged; ratio under budget" for BOTH a genuine limit point AND a structure with a missing restraint. The global dissipation ratio dilutes a localized diffused mechanism on a large mesh.
**Mitigation:** add INDEPENDENT FALSIFYING signals the log can use to contradict the classifier — the micro-burst (orthogonal physics), the consecutive-stabilize abort cap, and per-REGION dissipation policing (not just global). Emit a top-level `verdict` field (`equilibrium | regularized | unverified | aborted`) that a paper-grade consumer must acknowledge. The log must be allowed to say "I do not trust this result," not merely "here is what I did."

### R-SCOPE — the converged design is over-engineered relative to the MVP **[MAJOR]**

The shipped `adaptive_static` (~70 lines) already clears both RC-3D gates; every differentiating feature is net-new code. The micro-burst classifier in particular is a new stateful cross-paradigm observability-dependent subsystem being added to the "lowest-risk" design to fix its weakest link.
**Mitigation:** ship the MVP and STOP. **v1 = Ramm sizing + JSONL log + STABILIZE reflex (over shipped C++) + the abort-on-runaway-dissipation guard + Layer-1.5 queries.** HARD-DEFER to ranked follow-up ADRs (gated on a real production model that demonstrably needs each): the micro-burst classifier, b24ac routing, the recipe-pack DSL, the phase machine's deeper rungs, snap-back path-following, the deep-C++ `RobustNewton`, the inertia query, and `LadrunoDissipationArcLength`. For the MVP's instability decision, the cheap honest heuristic + the consecutive-stabilization abort guard ALONE prevents the catastrophic failure (stabilizing through a real mechanism) without the dynamics subsystem.

### R-TCL / R-PYTHON-DRIFT — Python-only controller orphans Tcl and drifts from the core **[MINOR]**

A Python-only spine gives Tcl workflows none of the robustness and is invisible to the C++ core; the relied-upon revert contract could silently change under an upstream merge.
**Mitigation:** accept Python-only for v1 explicitly and bound the claim; carry the ADR-20 §8 #5 "lift increment-sizing into a `StaticAnalysis` flag" as the tracked dual-interpreter path (not built speculatively); ship the revert-contract regression test so any drift breaks the build LOUDLY; keep the v1 Python layer thin enough that a future C++ lift is a port, not a rewrite.

### R-SMS / R-PARALLEL — dynamics affordability and parallel signals **[MINOR]**

Selective mass scaling (the affordability enabler for the explicit/DR path on production 3D meshes) is roadmapped, not built — so the dynamics rung's `dt_cr` is throttled by the smallest element and may be 10-100x too slow. Separately, distributed LU solvers have no global inertia and `EnergyBalanceRecorder` records per-partition with double-counted shared-node ULW.
**Mitigation:** honestly bound the dynamics rung to coarse/uniform meshes until SMS ships (cap added-mass fraction ~5%, refuse physical mass-scaling on rate-dependent materials); for distributed runs drop to signal-independent detectors (micro-burst KE via collective reduction; the abort cap) and DISABLE rigorous-inertia claims; compute the dissipation gate from an ALL-RANK reduced sum and stamp global-vs-partition-local in the log. Multi-rank robustness is flagged but not solved here.

### R-REPRO — "reproducible" is overclaimed across builds **[MINOR]**

Rung routing branches on FP/ordering-sensitive `testNorms`/KE traces near decision boundaries; a log records what happened on one run, not a cross-build recipe.
**Mitigation:** distinguish an AUDIT log (what happened — always emitted) from a REPLAY manifest (build/version hash, exact thresholds, numberer/solver/partition, per-step pass/fail bitmap). Bit-identical rung-sequence replay is asserted only on a fixed build; cross-build replay is best-effort and labelled as such.

## Validation plan — the torture-test battery (doubles as a V&V seed)

Each fixture is small, has a reference, and asserts BOTH convergence AND fidelity (the answer did not silently shift). All run under the v1 driver with `verify=True` and a JSONL log; each asserts a `verdict` and a diffusion bound.

1. **Snap-through arch** (`tests/_proto_arch_snapthrough.py`, the shipped prototype oracle: far-branch `uy ≈ -0.217`). Asserts: the driver reaches the far branch; STABILIZE (if armed) carries it with `peak_load_drift` under the c-reduction bound; and — critically — every phase machine path (pure cutback, STABILIZE, dynamics+polish) lands on the SAME far branch. Exercises rung-1 and the limit-point-not-instability discrimination.
2. **Softening RC prism** (single `stdBrick`/`LadrunoQuad` with `ASDConcrete3D`, the Gate-2 confined-column physics). Asserts: completes where naive Newton diverges (already proven for `adaptive_static`); the stabilized peak load matches the unstabilized cutback peak within the diffusion bound; the revert-parity + revert-idempotency battery passes on the softening material. Exercises R-DIFFUSION, R-REVERT.
3. **Buckling / post-buckling column** (Euler column with imperfection). Asserts: the driver passes the limit point WITHOUT routing into dynamics (indefinite tangent ≠ instability — R-INSTAB); the micro-burst classifies the limit point as "bad step / rest IS the answer," NOT "genuine instability." Guards against the over-escalation that defeats production affordability.
4. **Slack-cable snap** (pre-tensioned cable losing tension, a true snap-back). Asserts: snap-back is DETECTED and loudly journalled; STABILIZE is NOT silently accepted as an equilibrium curve (R-SNAPBACK); routing to `LadrunoIndirectControl` with a supplied control DOF traces the branch; absent a control DOF the run is stamped `verdict="regularized"`, not `"equilibrium"`.
5. **Contract regression suite** (not a physics fixture, but gating): the revert-contract return-code test (R-REVERT-RC), the switch-from-true-equilibrium `phat` test (R-SWITCH-PHAT), the descending-branch predictor-sign test (R-SWITCH-SIGN), the handoff time/λ/velocity test (R-HANDOFF), and the no-equilibrium termination test (R-TERMINATION). These break the build LOUDLY if a relied-upon invariant drifts.

A passing battery is the v1 acceptance gate AND the seed for the published V&V appendix (each fixture carries its own diffusion bound and verdict in its log).

## Reserved class tags

**v1 needs NO new Integrator/Algorithm class.** The Layer-1.5 observability seam is named getters + `OPS_` subcommands on the EXISTING `LadrunoArcLength` (33004) and `LadrunoDynamicRelaxation` (33005) — no new class tag, but a `Vector` payload widening on 33004's `sendSelf/recvSelf` for `cVisc`/dissipation state (legacy-read size-guarded, recorded in [[LEDGER_implementations]]).

Deferred class tags (RESERVED in the **Integrator** registry per [[LEDGER_implementations]] convention; bands are per-registry — these do not collide with identical ELE/MAT/RECORDER numbers):

- `INTEGRATOR_TAGS_LadrunoDissipationArcLength = 33007` — already reserved by [[22_ladruno_dissipation_arclength_adr]] (scoped, no code). The snap-back keystone (R-SNAPBACK); the first C++ deliverable of the follow-up roadmap IF a real snap-back model needs it.
- `INTEGRATOR_TAGS_LadrunoRobustNewton = 33008` — **RESERVED, deferred**. Only if the deep-C++ in-loop escalation (rejected Design 3) is ever justified; next free in the integrator band after 33007. Must ship with a gates-off bit-identical parity test against `NewtonRaphson` and must never edit `NewtonRaphson`.

(Current integrator band usage, verified `SRC/classTags.h:1136-1140`: 33000/33002 ExplicitBathe(LNVD), 33003 CentralDifferenceLadruno, 33004 LadrunoArcLength, 33005 LadrunoDynamicRelaxation, 33006 LadrunoIndirectControl; 33007 reserved-no-code.)

## References

- **Belytschko, Liu, Moran, Elkhodary** — *Nonlinear FE for Continua and Structures* (2014), **Ch. 6** (solution methods: Newton/quasi-Newton, line search, arc-length, stability of equilibrium).
- **Crisfield** — *Non-linear FE Analysis of Solids and Structures*, Vol. 1 (1991): arc-length / continuation, limit points, the spherical/cylindrical constraint and its discriminant.
- **Riks, E.** (1979) "An incremental approach to the solution of snapping and buckling problems," *IJSS* 15, 529 — the path-following parent of `*RIKS` and `LadrunoArcLength`.
- **Abaqus** *Analysis User's Guide* — *STATIC, STABILIZE* (automatic stabilization by dissipated-energy fraction, default `2e-4`; attached to `*STATIC`, never to `*RIKS`) and automatic incrementation (`I_d` target-iteration sizing). The faithful spine analog.
- **Ramm, E.** (1981) "Strategies for tracing the nonlinear response near limit points," in *Nonlinear FE Analysis in Structural Mechanics* — the target-iteration step-sizing rule the spine uses.
- **Scott, M.H. & Fenves, G.L.** (2006) "Plastic hinge integration methods for force-based beam-column elements," *J. Struct. Eng.* 132 — and the broader Scott–Fenves regularization / characteristic-length localization-limiter line motivating why softening needs a length scale (the materials this driver pushes through; [[26_ladruno_plane_frontier_adr]]).
- **Kolay, C. & Ricles, J.M.** (2014) "Development of a family of unconditionally stable explicit direct integration algorithms with controllable numerical energy dissipation," *EESD* 43, 1361 — KR-α; the unconditionally-stable explicit option for the dynamics fall-through alongside `ExplicitBatheLNVD` (Noh–Bathe).
- Internal: [[20_ladruno_arclength_stabilized_adr]], [[21_ladruno_dynamic_relaxation_adr]], [[22_ladruno_dissipation_arclength_adr]], [[26_ladruno_plane_frontier_adr]], [[21_rc3d_validation_gates]], [[LEDGER_implementations]], [[LEDGER_quirks]].

## Implementation log

### 2026-06-16 — torture battery executed on a live build; Layer-0 prototype

First runnable validation on the freshly-built `dist` (`ca98b3ccb`). Fixtures under `Ladruno_scripts/robust_solve_tests/`:

- **`torture_softening.py` — PASS.** Single `Concrete02` truss, load target 1.3×peak. LoadControl+Newton dies at the strength peak (strain −0.00185, no state beyond); adaptive halving + the **full** algorithm ladder *also* dies at the peak (strain −0.00200 — there is no far branch to find); DisplacementControl and the rung-3 switch trace the softening branch to strain −0.010. The clean load-control killer and the RC-relevant motivator.
- **`torture_snapthrough.py` — runs; surfaced a finding.** Calibrated snap-through limit `λ_lim ≈ 3.81`. LoadControl dies at the limit, **but** adaptive halving can *dynamically jump* the unstable branch to the far stiffening branch (`λ ≈ 5.72`) — it "succeeds" while **skipping the true equilibrium path**. Folded into R-LOG-MASK: a converged load-control step across a snap-through is NOT proof the traced path is physical; the driver must flag branch-jumps, not only non-convergence. (Snap-through has a far equilibrium; snap-back — fixture 4 — does not, and is the unambiguous case.)

Scenario (d) in each fixture is the Layer-0 rung-3 prototype (load→displacement constraint switch). Next: promote into `robust_drive()` (`adaptive_static` + rung-3 switch + rung-5 dynamics fall-through + JSONL decision log) and wire both fixtures as the pytest acceptance gate.

*(filled in once the plan is being executed; move to `Ladruno_internal/implemented_robust_solve_driver.md` when done.)*