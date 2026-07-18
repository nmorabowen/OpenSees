---
title: "ADR-73 implementation plan — LadrunoPorousOverlay execution (P0 detailed, P1 concrete pins, P2–P4 sketch)"
status: working doc — execution plan, updated as phases land
---

# ADR-73 implementation plan

Companion working doc to [[73_ladruno_porous_overlay_adr]] (merged PR #569;
explicit-lane amendment PR #570). This is the *execution* plan: work packages
(WP), owners, inputs, done-whens. Delete or fold into the ADR's Implementation
log when the runway closes.

## Execution model

- **Main session ("MAIN", Fable)** — interface pinning, adjudication, builds,
  PRs, CI. Never delegated.
- **Opus 4.8 agents ("OPUS")** — substantive production WPs (E7 experiments,
  pattern class, parser, batteries, guide) and the P1 adversarial panel.
  Parallel where files are disjoint.
- **Measurement protocol (P0's analog of ADR-71's independence protocol):**
  the coupled-block math is already twin-verified (ADR-71 P0) and the split
  math spike-verified (PR #567) — P0 here is *measurement*, not re-derivation.
  Every pinned constant (sequencing error, CFL envelope, subcycle curve) must
  come out of E7 with the prediction it confirms/refutes stated FIRST in the
  test header (no post-hoc "expected" values).
- **Branch/PR discipline:** one PR per phase, one branch per PR, base
  `ladruno`, verify PR state==OPEN before every push (auto-merge strands
  follow-ups). Phase branches: `feature/adr73-p0` … `-p4`.
- **Windows gotchas (non-negotiable, from memory):** Write tool for files
  (Bash heredocs eat backslashes); build `cmd /c Ladruno_scripts\build.bat …`;
  worktree Python tests `python -S` + manual paths + `assert
  opensees.__file__`; new sources into `stamp_headers.py` GLOBS + stamped;
  close VS Code before DLL-touching steps.

---

## P0 — toy E7: the measured pins (numpy only, no OpenSees build)

Extends the merged spike (`adr71_meshless_p_spike/`) with a NEW file
`staggered_pins_e7.py` importing `meshless_p_toy.py` as a library — the merged
toy stays frozen as ADR-cited evidence. E7 adds the one thing the toy lacks:
**solid inertia** (lumped M, central-difference solid loop), because the
explicit-lane pins are dynamic questions.

| WP | Owner | Deliverable | Done when |
|---|---|---|---|
| 0.A | MAIN | this plan committed; predictions table pinned below | doc on the P0 branch |
| 0.B | OPUS-1 | `staggered_pins_e7.py` — experiments E7.1–E7.5 | all pins measured, RESULTS.md §E7 written, plots |
| 0.C | MAIN | adjudicate measurements vs predictions; freeze the P1/P3 contract constants in THIS doc; PR | PR merged; pins recorded below as FROZEN |

**E7 experiments (each states its prediction first):**

- **E7.1 — fs1-without-final-resolve** (the cheap v1 element flow: solid
  analyze with p from the last commit → commit → fluid advance; NO second
  solid solve). Prediction: same O(Δt) family as fs1-with-resolve, constant
  factor ≤ 2× worse. Measured over Δt-halving on Terzaghi + the near-undrained
  column. **Decides whether the v1 driver needs the extra solid solve per
  step** (it is the difference between "overlay works with plain `analyze`"
  and "overlay needs a driver from day one").
- **E7.2 — explicit solid + implicit fluid (the P3 lane), stability
  envelope.** CD solid (lumped M) with overlay forces; fluid advanced at
  commit. Prediction (ADR §3.4; hedged ⟨A-7⟩ — the ZPC-1988 proof covers
  their scheme, not our frozen-force variant): stable at the DRAINED-speed
  CFL (the implicit fluid absorbs the undrained stiffening). Sweep Δt across
  [0.5·Δt_drained, 1.2·Δt_undrained⁻¹-scaled] and locate the empirical
  boundary. THE load-bearing P3 pin — if the boundary lands at the undrained
  speed instead, the P3 lane's Δt advisory defaults change, not the
  architecture.
- **E7.3 — multirate curves, BOTH directions.** (a) Subcycle N ∈
  {1, 2, 5, 10, 20, 50} (fluid slower — explicit lane): accuracy vs
  monolithic + stability; prediction: error ~linear in the sync interval
  N·Δt while N·Δt ≪ h²/c_v; pins the `-subcycle auto` formula. (b) Substep
  M ∈ {1, 2, 5, 10} (fluid faster — implicit consolidation lane): coarse-Δt
  Terzaghi with a step load, M BE fluid sub-steps per solid step, Δu
  injected as a linear ramp; prediction: recovers the early-time pressure
  transient a single BE step smears, at reused-factor cost. Pins the
  `-substep` accuracy claim.
- **E7.4 — fully-explicit fluid (P3b).** Lumped S*, forward p-step.
  Predictions: diffusion CFL slack by orders (numbers in ADR §3.4); coupled
  stability at the UNDRAINED-speed CFL, factor √((M_oed+K_f/n)/M_oed) vs the
  drained pencil measured within ±20 %; L not needed (no iteration).
- **E7.5 — removal-step consistency on the dynamic path.** Crack opens
  mid-CD-march (E4 scenario + inertia): prediction — physical stress-wave
  transient only, artifact quantified as the p-jump component that scales
  with 1/Δt (gate: none above solver tolerance); `-onRemoval drain`
  kFactor sweep {1, 10, 100} monotone in drainage time.
- **E7.6 — L-variant decision data ⟨A-2⟩.** classic α²/K_dr vs oedometric on
  BOTH geometries (constrained column AND footing), compressible and
  near-undrained: iters mean/max + any divergence. Prediction: classic never
  diverges (KTJ proof); oedometric converges on both at ~3× fewer iters.
  This is the number pair behind the `-fsL` default/opt-in split.

**P0 exit:** ADR §7 P0 row verbatim + the FROZEN-constants block appended to
this doc (E7.1 default, E7.2 envelope, E7.3 auto-N formula, E7.4 factor).

---

## P1 — `LadrunoPorousOverlay` fs1 implicit, end-to-end

Branch after P0 merges. The ADR fixes the architecture (§3–§4); the pins below
fix the class mechanics. **Full adversarial panel at P1** (ADR §7 policy —
core touch).

| WP | Owner | Deliverable |
|---|---|---|
| 1.A | MAIN | pins below frozen; classTags.h landing plan confirmed |
| 1.B | OPUS | `SRC/domain/pattern/ladrunoPorousOverlay/LadrunoPorousOverlay.{h,cpp}` — pattern class per pins (snapshot, blocks via the ADR-71 kernel headers, CSR+CG fluid solve, advanceTrial/commitFluid/revertFluid, applyLoad via addUnbalancedLoad, removal rescan, `-record` CSV) |
| 1.C | OPUS | `OPS_LadrunoPorousOverlay.cpp` parser (surface per ADR §4.1 incl. `-layer` blocks + v1 `-moduli`; unknown-flag-FATAL; factor-ignored notice) + registration: classTags.h `PATTERN_TAG_LadrunoPorousOverlay 33022` + cross-registry comment, pattern dispatch (Tcl + OPS), FEM_ObjectBrokerAllClasses, CMakeLists, Domain::commit hook (one guarded `// Ladruno` block, contact-engine include precedent) + vanilla-ledger rows + stamps |
| 1.D | OPUS | sendSelf/recvSelf (every ctor arg + region + layers + committed p) + serial DB round-trip test |
| 1.E | OPUS ×2 ∥ | battery, split: **(i) physics** — Terzaghi vs analytic (tol from the E7.1 rate-form curve: late-window relL2 ≤ 2× the toy value at matched Δt, `e7_summary.json` e71 */late keys, + O(Δt) mutual-convergence order ≥ 0.9); two-leg vs monolithic LadrunoUP (mutual Δt-convergence, order ≥ 1); E4-in-OpenSees: real `remove element` crack, both `-onRemoval` policies, curves vs toy reference — **removal gates use the FIXED-WINDOW jump (E7.5: 20·Δt window, spread ≤ 10 %), never the single-commit jump**; sequencing contract test (step-load column: p(0⁺) ≈ q — kills the fluid-first p≡0 class); rate-form-vs-fixed-point contract test (fixed-point single-pass error ≥ 10× rate-form on the same march — pins that the reference rule is wired); L-floor loud-warning test; **(ii) framework** — pattern-timing bit-constancy (forces identical across Newton iterations of one step, Static AND Newmark AND CD); factor-ignored semantics; revertToLastCommit consistency (fluid state rolls with the domain, incl. dpCommitted_); region validation errors (unsupported cell, empty drained set, overlapping overlays on one element → fatal); pInit steady/hydrostatic vs closed form |
| 1.F | OPUS panel ×3 | adversarial gate: (1) **framework-reality critic** — applyLoad/addUnbalancedLoad semantics across integrators, commit ordering, revert paths, hook placement, Rayleigh non-interaction; (2) **math/contract critic** — per-cell Q sign (force = +Q·p, compression-positive p), L assembly, S*/H̃ per-layer, sequencing vs E7 pins, BE window Δu = committed-minus-stored; (3) **robustness critic** — removal edge cases (all cells dead, drained node on dead cell, re-add?), restart from DB, zero-k layer, single-element region, mixed 2D/3D rejection |
| 1.G | MAIN | build, run battery (`python -S`), integrate, ledger amend (33017-row untouched; 33022 row → P1 landed), LEDGER_quirks rows, PR, CI watch |

### WP1 pins (FROZEN at 1.A — agents implement exactly this)

**Command form** (no-braces pattern, H5DRM precedent):

```tcl
pattern LadrunoPorousOverlay $tag -region {$e1 ...} -Kf $Kf -rhoF $rf \
    -perm $k1 $k2 <$k3> -poro $n -moduli $E $nu \
    <-layer {$eles} -perm ... <-poro $n> <-moduli $E $nu>> ... \
    <-alpha $biotAlpha> <-Ks $Ks>            \;# defaults 1.0 / infinite ⟨A-4⟩
    <-thick $t>                              \;# 2D regions, default 1.0 ⟨A-4⟩
    -drained {$n1 ...} <-pInit ...> <-stab ...> \
    <-fsL classic | oedometric | $scale>     \;# default classic α²/K_dr ⟨A-2⟩;
                                              ;#   $scale floors at oedometric
    <-staticMode hold | steady>              \;# static commits: domain "time" is
                                              ;#   the load factor ⟨A-3⟩ — never
                                              ;#   march the fluid on Δλ
    <-onRemoval keep|drain $kF> <-fluidUpdate implicit> \
    <-subcycle auto|$N> \
    <-record $file <$everyN>> <-fluidBody ...> <-dynSeepage on|off>
```

- **Multirate (amended at 1.A per E7.3):** `-subcycle N` = fluid advanced
  every N solid commits with the accumulated Δu (explicit lane; fluid
  slower) — E7.3a confirmed (err ~ N^1.2, θ = 0.089 pins
  `-subcycle auto` → N = max(1, floor(0.09·h²/(c_v·Δt)))); removal events
  force a sync regardless of N. **`-substep M` is REMOVED from the surface**
  (E7.3b refuted the accuracy claim: error flat over M = 1…10; the parser
  rejects it with a pointer to the ADR §12 P0 entry — unknown-flag-FATAL
  covers it naturally, no special case).

- **v1 `-moduli $E $nu` is REQUIRED** (global or per-layer): the overlay does
  not own the solids' materials, so L and α-stab moduli come from explicit
  user input in v1 — a deliberate simplification vs the ADR's
  "auto-from-material" wording (adjudicated here). P2 revisits with the
  updateMaterialStage/parameter transport (PDMY staging needs it anyway).
- `-fluidUpdate explicit` parses but is fatal-NYI until P3b. `-subcycle auto`
  parses at P1 using the E7.3 formula (cheap — needs only c_v from
  moduli+perm).

**Class mechanics:**

- **Snapshot** (first setDomain/applyLoad): for each region element —
  `getNodePtrs()`/`getExternalNodes()`, classify by (ndm, nNodes) into the
  ADR-71 providers (T3/Q4/H8 v1; else fatal with named error). Store per-cell:
  node indices into the region map, per-GP (Np, dNp, dv) from the provider,
  per-cell dense Q_e (via kernel `addQ`). **No global Q matrix** — coupling is
  applied cell-wise both directions.
- **Fluid system:** CSR over region nodes; H, S* = (n/K_f per layer)·M_p +
  α_stab·H̃, **L per `-fsL`: default classic (α²/K_dr per layer)·M_p,
  opt-in oedometric (α²/(K_dr+4G/3))·M_p** via kernel addH/addS/addHtilde.
  Drained BCs by row/col elimination. **Solver pin: CG + Jacobi, rel tol
  1e-10, warm-started from the previous p** — numbering-independent, zero new
  dependencies; direct-factor upgrade only on profile evidence (ADR-40
  discipline).
- **Fluid advance = the RATE form (P0 E7.1 pin — ADR §12 P0 entry
  governs):** `(S* + L + ΔtH)·p_trial = (S* + L)·p₀ + L·Δp_ref − Qᵀ·Δu
  + Δt·f_seep`. **Reference rule:** the FIRST advance of a step uses
  `Δp_ref = dpCommitted_` (previous committed fluid increment → rate form =
  fs1, O(Δt)); any REPEATED advance within the same step uses the last
  trial increment `p_trial_prev − p₀` (→ the ADR §3.1 fixed-point
  iteration; L cancels at convergence = monolithic BE). The fixed-POINT
  form must never run single-pass (O(1) storage error, measured 41–61 %).
- **State discipline (P2-proofed now):** `pCommitted_`, `dpCommitted_`
  (previous committed fluid increment — the rate-form reference),
  `uSnapshot_` (region u at last fluid advance), and THREE methods —
  `advanceTrial()` (solve fluid from current committed/trial u per the
  reference rule above, refresh forces, do NOT commit), `commitFluid()`
  (dpCommitted_ ← pTrial − pCommitted_, then pTrial→pCommitted_, refresh
  uSnapshot_), `revertFluid()` (drop trial, reference rule resets to
  committed). fs1 = advanceTrial+commitFluid back-to-back inside the
  Domain::commit hook. The P2 iterated driver calls advanceTrial repeatedly
  before one commitFluid — no P2 rewrite.
  Serialization pin ⟨A-10, amended 1.A⟩: `pCommitted_` AND `dpCommitted_`
  travel in sendSelf (dpCommitted_ is NOT re-derivable); `uSnapshot_` is NOT
  serialized — recvSelf re-derives it from committed node displacements
  (they travel with the nodes), keeping restart exact.
- **Region-overlap guards ⟨A-13⟩:** at snapshot, fatal if (i) any region
  element is claimed by another LadrunoPorousOverlay (overlays discoverable
  via Domain::getLoadPatterns — classTag 33022 scan), or (ii) any region
  element is a monolithic u-p element (classTag 33017 LadrunoUP or upstream
  UP tags) — double-counted fluid. Battery covers both fatals.
- **Force application** ⟨A-5 verified mechanics⟩: override `applyLoad(time)`
  → per region node `node->addUnbalancedLoad(+Q·p_committed slice)` — the
  exact H5DRM mechanism (H5DRMLoadPattern.cpp:897). Call-path facts (pinned,
  verified): transient integrators apply loads ONCE per step via
  `AnalysisModel::updateDomain(time,dT)` → `Domain::update(newTime,dT)` →
  `applyLoad` (Domain.cpp:2381); iteration-time `updateDomain()` skips load
  application; `formUnbalance` only READS nodal unbalance
  (IncrementalIntegrator.cpp:145–165); `Domain::applyLoad` zeroes all nodal
  unbalance first (Domain.cpp:2060–2074). So commit-refreshed forces reach
  the solid at the NEXT newStep — no re-application assumption anywhere.
  ArcLength-family integrators call applyLoadDomain mid-step → documented
  UNSUPPORTED with the overlay at P1 (parser cannot see the integrator; the
  guide + a battery smoke own it). Sign: force = **+Q·p** on the skeleton
  (compression-positive p; toy convention `u = K⁻¹(f + Q p)`; the
  sequencing contract test pins it physically). Bit-constancy test observes
  `ops.nodeUnbalance` between iterations under Static, Newmark, AND CD.
- **Hook:** `Domain::commit()` gains one guarded block at the ADR-39
  contact-hook slot (Domain.cpp:2233 — AFTER node/element commitState,
  BEFORE `committedTime = currentTime; dT = 0` at :2246–2247, so committed
  disps are final and Δt is still readable there ⟨A-5⟩), following the
  `LadrunoContactDomain` include precedent (Domain.cpp:82). Honors
  `-subcycle`; under a STATIC analysis (detected per `-staticMode`, ⟨A-3⟩)
  it never marches the fluid — `hold` keeps p, `steady` re-solves
  H·p = f_seep. Fallback if the panel kills the hook: applyLoad new-step
  detection (documented, not implemented speculatively).
- **Removal rescan:** at each hook firing, region elements checked against the
  Domain; newly-dead cells → Q_e dropped, H/S per `-onRemoval` (drain →
  k×kFactor in those cells, refactor/re-setup CG operator). Dead-cell GP set
  never shrinks the fluid measure under `keep`.
- **`-record $file <$everyN>`**: CSV `time, p(node1..nodeN)` — the gate-facing
  output until P4's proper channels. Region node order written once as a
  header row; file opened truncate (restart overwrites — documented).
- **dynSeepage default off** (ADR-71 §12 amendment inherited).
- **Test files ⟨A-11⟩:** `tests/test_ladruno_overlay_physics.py` (1.E-i),
  `tests/test_ladruno_overlay_framework.py` (1.E-ii) — Zone-B (fork-local)
  per the testbed convention; `python -S` runner discipline. New C++ files
  into `stamp_headers.py` GLOBS — owned by 1.C.

**LEDGER_quirks rows queued at P1:** (i) overlay pattern factor/TimeSeries
ignored by design; (ii) honest-p LadrunoUP + upstream CentralDifference =
Richardson-unstable p (from ADR §7 P3 — landed early with the P1 battery since
the test needs only monolithic pieces); (iii) one-overlay-per-water-body rule.

---

## P2 — iterated driver + stage transport (PINNED at 2.A, 2026-07-17)

The sketch is superseded by these pins. The overlay seam
(`advanceTrial(firstOfStep, dt)` / `commitFluid` / `revertFluid` with the
unified reference rule) is ALREADY BUILT and battery-pinned at P1 — P2 drives
it; no seam re-derivation.

| WP | Owner | Deliverable |
|---|---|---|
| 2.A | MAIN | these pins committed |
| 2.B | OPUS | `LadrunoStaggeredDriver.{h,cpp}` shared core + overlay latch/param seam + registration (both interpreters) |
| 2.C | OPUS | battery `tests/test_ladruno_overlay_driver.py` (a–f below) |
| 2.D | OPUS ×3 | adversarial panel (framework-reality / math-contract / robustness) |
| 2.E | MAIN | build, battery, panel adjudication, ledgers, ADR §7/§12, PR |

### 2.1 Command surface (FROZEN)

```tcl
LadrunoStaggeredAnalyze $n $dt <-tol $t> <-maxIter $k> <-pScale $s> <-verbose>
LadrunoStaggeredAnalyze -stats
```

- Defaults: `tol = 1e-6`, `maxIter = 500`, `pScale = 1.0` (the e76 protocol
  values; e76 itself used pScale = q — batteries pass `-pScale` to mirror).
- Registered in BOTH interpreters (profiler-command precedent): classic
  `SRC/tcl/commands.cpp` `Tcl_CreateCommand`; interpreter
  `OPS_LadrunoStaggeredAnalyze()` in `OpenSeesCommands.cpp/.h` + wrappers in
  `TclWrapper.cpp` / `PythonWrapper.cpp`. Shared core = free functions in
  `SRC/domain/pattern/ladrunoPorousOverlay/LadrunoStaggeredDriver.{h,cpp}`
  (stamped; CMakeLists of that dir).
- **TRANSIENT-only v1.** Requires the active DirectIntegrationAnalysis
  (classic global `theTransientAnalysis`; interpreter
  `cmds->getTransientAnalysis()`). A static analysis → loud fatal citing
  ⟨A-3⟩ (static domain "time" is the load factor — an iterated Δλ fluid
  march is silently wrong physics). VariableTimeStep analysis drives fine
  through the base API (fixed driver dt).
- Analyze form returns 0 / negative (the `analyze` convention). `-stats`
  returns the LAST run's telemetry as a list
  `{nSteps totalFluidSolves meanIters maxIters lastResidual maxResidual}`
  (Tcl list / Python list of floats; zeros before any run).
- Driven set = every load pattern in the domain with classTag
  `PATTERN_TAG_LadrunoPorousOverlay` AND `staticMode == SM_MARCH` AND
  `fluidUpdate == FU_IMPLICIT`. Empty driven set → loud fatal (misuse).
  HOLD/STEADY overlays are NOT driven and keep their normal hook semantics.

### 2.2 Step recipe (FROZEN — decomposed analyzeStep + toy `fs_iter` verbatim)

Per driver step (toy `qs_staggered(mode="fs_iter")` semantics; iteration
count = number of fluid advances, `nit` in e76):

0. **Catch-up**: before step 1, any driven overlay with a pending
   `-subcycle` window (`commitsSinceSync_ > 0`) gets an early fs1 sync
   (`advanceTrial(true, dtAccum_)` + `commitFluid` + counter reset), one-time
   notice — otherwise the first driver advance would pair a multi-commit Δu
   window with the driver's single dt.
1. Latch ON for the driven set (`setExternalStepping(true)`,
   `setResidualPScale(pScale)`); RAII guard clears it on EVERY exit path.
2. `model->analysisStep(dt)`; `dia->checkDomainChange()`;
   `integrator->newStep(dt)`; `algorithm->solveCurrentStep()` — the
   iteration-1 solid solve (overlay forces `+Q·pTrial`, and
   `pTrial == pCommitted` at step start).
3. `advanceTrial(firstOfStep = (nit==1), dt)` on each driven overlay.
   First advance of a step uses the committed-increment reference
   (`dpCommitted_` — the unified reference rule; a better-than-toy start:
   e76's first iterate used Δp_ref = 0), repeats use the trial increment
   (§3.1 fixed point; L cancels at convergence).
   Residual = max over driven overlays of
   `‖p_k − p_{k−1}‖₂ / (‖p_k‖₂ + pScale)` over the FREE (non-drained) rows —
   the toy `dcon` on the reduced system (amended at 2.E per panel math-6:
   including prescribed drained values in the denominator would under-report
   the residual). Computed inside `advanceTrial`, exposed as
   `lastAdvanceRelChange()`. The first-of-step advance is the rate-form warm
   start, NOT the toy fs_iter's cold Δp_ref = 0 first iterate (panel math-2:
   same fixed point, typically fewer iterations — gate (b)'s band absorbs it).
4. While residual ≥ tol and nit < maxIter:
   `domain->revertToLastCommit()`; `integrator->revertToLastStep()`;
   `newStep(dt)`; `solveCurrentStep()` (solid now sees `+Q·p_k`);
   `advanceTrial(false, dt)`; recompute residual.
5. Converged: **final momentum resolve** (§3.1 step 3 / toy line 219 —
   revert; `newStep`; `solveCurrentStep` with `p_K` forces), then
   `integrator->commit()` → `Domain::commit` → hook fires LATCHED (2.3).
6. `-verbose`: per-step `step i: iters=k residual=r`. `flushRecorders()` on
   the success path (the `analyze` convention).

**Abort semantics (all-or-nothing per step, FROZEN):** on newStep /
solveCurrentStep / advanceTrial failure at ANY iteration, or maxIter reached
unconverged: `domain->revertToLastCommit()`; `integrator->revertToLastStep()`;
`revertFluid()` on all driven overlays; latch cleared; loud print (step, nit,
residual); return negative (−2 newStep, −3 solve — the analyzeStep codes;
−6 fluid solve; −7 maxIter). Prior committed steps stand. e76 says classic-L
never hits maxIter on sane problems — a hit is pathology, never silently
accepted.

Mechanics verified in-tree at 2.A: `checkDomainChange()` public
(DirectIntegrationAnalysis.cpp:619) = the analyzeStep stamp check;
`Domain::revertToLastCommit()` (Domain.cpp:2295) resets
`currentTime = committedTime`, re-applies loads (unbalance zeroed on every
`Domain::applyLoad` — no `+Q·p` accumulation) and does NOT touch overlay
fluid trial state (no overlay hook there — exactly what the iteration
needs); `Newmark::newStep/revertToLastStep` (Newmark.cpp:157/270) form a
consistent repeat-same-step pair (`Ut ← U` after restore keeps the committed
predictor base); `Domain::analysisStep` is a base no-op.

### 2.3 The overlay latch (FROZEN — the double-advance guard)

New overlay state `externalStepping_` (transient bool, false in both ctors,
NOT serialized, `recvSelf` forces false — a driver run is synchronous, no
send can observe it true; pinned defensively) plus `pScale_` (transient,
default 1.0). Under the latch:

- `applyLoad` injects `+Q·pTrial_` (the iterate's forces) instead of
  `+Q·pCommitted_`. Outside iteration the two are equal (`commitFluid`
  syncs), so unlatched P1 behavior is bit-untouched. The latched path also
  covers `Domain::revertToLastCommit`'s internal re-applyLoad (harmless —
  zeroed and re-applied at the next newStep).
- `onDomainCommit` SM_MARCH path = `rescanRemoval` (unchanged) +
  `commitFluid()` + record row + subcycle-counter reset — **NO
  advanceTrial** (the fluid trial is already the converged iterate; an
  advance here would double-march) and NO window accumulation. The ONLY
  commit that fires while latched is the driver's converged commit.
  HOLD/STEADY paths unchanged (those overlays are never latched).
- `-subcycle` interplay: the driver syncs every step by construction;
  configured N > 1 → one-time advisory ("driver overrides -subcycle").
  Counters zeroed at each latched commit so post-driver plain `analyze`
  restarts its window cleanly.
- `revertFluid()` semantics unchanged.

### 2.4 Moduli / stage transport (FROZEN — the parameter route)

- `LadrunoPorousOverlay::setParameter/updateParameter` reached via the
  EXISTING `parameter $ptag loadPattern $overlayTag <arg>` surface — both
  interpreters already dispatch loadPattern targets
  (OpenSeesParameterCommands.cpp:156, TclParameterCommands.cpp:235; no
  vanilla parameter-command touch needed):
  - `E` / `nu` — overlay-global moduli;
  - `layerE $i` / `layerNu $i` — 1-based `-layer` declaration index; setting
    one establishes a layer override where the layer inherited global.
  - `updateParameter` validates with the parser gates (E > 0, 0 ≤ ν < 0.5),
    sets, marks `moduliDirty_`.
- `moduliDirty_` is honored LAZILY at the next fluid use (advanceTrial /
  steady solve): re-resolve per-cell E/ν (global/layer), re-derive per-cell
  L factor (per `-fsL` mode incl. manual-scale floor) and `-stab auto`
  α per cell, REASSEMBLE `aS_` (S + α_stab·H̃) and `aL_` from the retained
  per-GP snapshot data. `aH_`/`fseep_` untouched (perm-only). The CG
  operator is a (cS,cL,cH) combination of those arrays → refreshed for
  free. `maxCv_` recomputed; an already-resolved `-subcycle auto` N stays
  (resolve is once-per-run, documented; the driver bypasses subcycle
  anyway).
- **`updateMaterialStage` does NOT and CANNOT reach the overlay**
  (MaterialStageParameter registers only the first accepting element and
  never scans load patterns — the ADR-71 sibling-broadcast trap,
  family-documented). The transport contract: after a stage flip the USER
  re-sets overlay moduli via the parameter route (PDMY battery + guide pin
  the recipe). A flip without a re-set keeps stage-0 L — stable but
  convergence-degraded; LEDGER_quirks row.
- Serialization: moduli values already travel in config; the dirty flag is
  transient (restore rebuilds the snapshot anyway).

### 2.5 Battery `tests/test_ladruno_overlay_driver.py` (Zone-B, python -S discipline)

- **(a) FIXED-POINT GATE (the gate of the phase):** driver at `-tol 1e-10`
  ≡ monolithic LadrunoUP at the SAME dt on the Terzaghi column (Q4;
  compressible AND near-undrained): relL2(p) and relL2(u) ≤ 1e-6 over the
  march (the E6 exact-equality leg, measured 2.8e-8/4.5e-7 in the toy).
- **(b) e76 telemetry gate:** e76 protocol (tol 1e-6, maxIter 500,
  pScale = q) on the column and a B4-like footing: classic-L no divergence,
  no maxIter hit, mean iters ≤ 1.3× the frozen e76 means (column 11.25 /
  footing 4.35); `-fsL oed` ≤ 1.3× (3.29 / 2.8). The driver's
  rate-form first iterate may only IMPROVE on e76's Δp_ref = 0 start —
  measured means reported in the PR either way (adjudication duty if the
  band breaks).
- **(c) B4 footing CB gate under stab:** driver on the near-undrained
  footing with `-stab auto` — checkerboard metric clean (vs the unstabbed
  control dirty), the ADR-71 B4 methodology.
- **(d) stage/moduli transport:** (i) BIT-TWIN gate: construct with E₁,
  update to E₂ via the parameter route BEFORE marching → p history
  bit-identical (≤ 1e-14) to a control overlay constructed with E₂ — gates
  the full L/α-stab/operator recompute path; (ii) PDMY staged liquefaction
  column (plain-ndf solids + PDMY + overlay + driver; `updateMaterialStage`
  flip mid-analysis + parameter-route moduli re-set) vs the ADR-71 P4
  monolithic reference (`tests/test_ladruno_up_element_pdmy.py` model);
  agreement band adjudicated from measurement, refutations recorded;
  (iii) mid-march L refresh twin-checked vs an oracle recompute.
- **(e) driver/hook interplay:** (i) fresh-model plain `analyze` after a
  driver run in the same process reproduces the P1 anchor (global-state
  leak guard); (ii) same model: plain segment → driver segment → plain
  segment — fluid state continuous, counters clean, record rows == syncs,
  march advisory printed once; (iii) abort path (forced maxIter fail):
  domain time/state rolled back to the last committed step,
  `pTrial == pCommitted`, subsequent plain analyze clean; (iv) `-stats`
  matches a hand-count on a tiny run; (v) driver under a static analysis →
  fatal, and driver with no qualifying overlay → fatal.
- **(f) RIDER:** re-verify the P1 observation "-pInit list overlay crashes
  on FileDatastore restore" post-#577 (may have been the same rho
  corruption). Clean → the list variant joins the DB gate family; still
  crashing → LEDGER_quirks row + report.

### 2.6 Panel (Opus ×3 — core-adjacent control flow, [[feedback_adversarial_gate_when]])

1. **framework-reality critic (load-bearing):** revertToLastCommit inside
   the loop, integrator state across repeated same-step newStep (Newmark
   family verified at 2.A; generality documented), latch coverage of every
   applyLoad path, abort paths incl. inner-analyze failure mid-iteration,
   recorder/commit counts, catch-up sync, sensitivity/subLevel
   non-interaction.
2. **math/contract critic:** driver loop ≡ toy `fs_iter` (residual norm,
   iteration counting, final-resolve placement), reference-rule legality of
   the rate-form first iterate, L cancellation at convergence, moduli
   reassembly vs kernel.
3. **robustness critic:** multi-overlay domains, maxIter edges, `-stats`
   before any run, param validation, latch serialization contract, PDMY
   staging edges, driver re-entry.

### 2.7 Registration / ledgers

- **No new classTag** → no classTags.h edit, no manifest.yaml row
  (check_manifest.py keys on ledger classTag rows — verified by running the
  CI checks locally).
- Vanilla touches (all already-Ladruno-marked upstream files):
  `SRC/tcl/commands.cpp`, `SRC/interpreter/OpenSeesCommands.{cpp,h}`,
  `TclWrapper.cpp`, `PythonWrapper.cpp` → LEDGER_vanilla_files rows.
- New files stamped + in `stamp_headers.py` GLOBS (dir glob may already
  cover); banner ADR-73 line amended to mention the driver
  (`banner_features.txt` → `patch_banner.py`).
- LEDGER_implementations 33022 row → P2; LEDGER_quirks: driver-overrides-
  subcycle, updateMaterialStage-cannot-reach-overlay, rider outcome.

## P3 — explicit lane (PINNED at 3.A, 2026-07-18)

The P3b sketch (below) is untouched. P0's §12 amendments GOVERN the lane:
**L = 0 at Δt ≤ 0.5× the discrete undrained pencil** (the rate form is
CD-unstable at any Δt; the drained-CFL fixed-point-L mode is a documented
accuracy-degraded opt-in NOT built here — refusing speculative surface).

| WP | Owner | Deliverable |
|---|---|---|
| 3.A | MAIN | these pins |
| 3.B | OPUS | `-fsL zero` lane + Δt_cr advisory augmentation + driver refusal + quirks/banner/ledgers |
| 3.C | OPUS | battery `tests/test_ladruno_overlay_explicit.py` (a–g below) |
| 3.D | OPUS ×3 | panel (charter below — advisory-correctness critic is load-bearing) |
| 3.E | MAIN | build, battery, adjudication, ADR §7/§12, PR |

### 3.1 Lane selection (FROZEN): `-fsL zero`

- `-fsL` gains keyword `zero` → `FSL_ZERO = 3` (enum + serialization
  compatible; values already travel). L ≡ 0: `advanceTrial`'s existing RHS
  degenerates to plain BE fluid `(S* + ΔtH)p₁ = S*p₀ − QᵀΔu + Δt·f_seep`
  (both L terms vanish; the reference rule becomes inert). NO new advance
  code path — the shipped formula with aL_ = 0.
- Parser: `zero` bypasses the oedometric floor (the floor guards the
  fs1/iterated lanes) and prints a ONE-TIME loud advisory: this is the
  EXPLICIT-lane setting; a quasi-static implicit fs1 march with L = 0 is the
  naive drained split and diverges in ~4 steps at soil coupling (measured,
  §3.2) — it will fail loudly, not wrongly.
- **P2 driver refusal (FROZEN):** `LadrunoStaggeredAnalyze` excludes
  `FSL_ZERO` overlays from the driven set with a LOUD FATAL (not a silent
  skip): iterating with L = 0 IS the drained split (KTJ-divergent at soil
  coupling). New const accessor `fsLModeCode()`.
- `-subcycle auto|$N` works unchanged under the lane (θ = 0.089 formula was
  MEASURED on the L=0 lane — E7.3a; the hook path already accumulates Δu
  across the window).

### 3.2 Overlay-aware Δt_cr advisory (FROZEN — the discrete undrained pencil)

- P0 pin: the advisory uses the **DISCRETE undrained pencil** (e74: material
  formula √(1+K_f/(n·M_oed)) is ~1.85× conservative, docs-only; discrete
  pencil ≈ 25 % margin vs the measured 1.32× boundary; e72: the L=0
  implicit-fluid boundary = 1.000× the pencil, frozen default 0.5×).
- Mechanism: per-cell **undrained stiffness augmentation**
  ΔK_e = Q_e·S_e⁻¹·Q_eᵀ (the element-local undrained condensation — dense
  nNp-rank block from the overlay's retained Qe + storage data; this is the
  only implementation that yields the DISCRETE pencil rather than the
  material formula). Overlay exposes
  `bool getUndrainedAugmentation(int eleTag, Matrix& Kadd) const`
  (node-major, first-ndm-DOFs layout — the same assumption the HRZ lumping
  makes; returns false for unowned/dead cells).
- Integration point: `CriticalTimeStep.cpp::computeCriticalTimeStep` (the
  SHARED per-element scan — CentralDifferenceLadruno report AND the ADR-36
  SMS pencil both consume it, so SMS scales against the corrected pencil
  exactly as §3.4 item 3 intends). Once per call: scan domain load patterns
  for 33022 overlays (classTag check, the P2 driver idiom); per element,
  if owned, fetch ΔK_e and add into the eigensolve's K. Fork-owned file —
  no vanilla ledger row. Report prints drained vs undrained governing
  pencil + implied factor + the material formula as documentation when any
  overlay contributed.
- Guards: element K size ≠ nNodes·ndm (exotic ndf) → LOUD one-time advisory
  + skip that element's augmentation (never silently optimistic); near-
  singular S_e (Kf → ∞ with stab off) → loud advisory naming the cell, skip
  augmentation, report the material formula as the (infinite) bound. S_e
  factor: dense Cholesky per cell with a rel-pivot floor; computed lazily
  ONCE and cached (moduli-dirty invalidates — S_e depends on α-stab).
- v1 scope: the augmentation reflects the overlay's CURRENT storage
  (post-removal rescan: dead cells return false). `useTangent` semantics
  unchanged (augmentation is state-independent).

### 3.3 Energy advisory (ADR-69, measured not asserted)

- Fact (verified 3.A): EnergyBalanceKernel's ULW = ∫vᵀP_ext dt reads
  `Node::getUnbalancedLoad` — the overlay's `+Q·p` forces are INSIDE the
  external-work channel by construction, so the closure residual accounts
  the coupling work as external load work.
- Deliverable: battery gate (g) MEASURES closure (ERR stays within the
  ADR-69 bound on a CD+overlay run — "documented, not silently absent");
  guide-queued note + LEDGER_quirks row stating the accounting (pore-
  coupling work rides ULW, not a separate channel; P4 may split it out).
  NO recorder code change at P3.

### 3.4 Battery `tests/test_ladruno_overlay_explicit.py` (gates)

- **(a) two-leg vs implicit monolithic on the B2/ZS84 column**: CD
  (CentralDifferenceLadruno, Δt = 0.4× advisory) + quad + overlay(-fsL
  zero) vs implicit Newmark LadrunoUP same mesh (ZS84 config from
  `tests/test_ladruno_up_element_analytic.py` B2, ~line 289): mutual
  Δt-convergence, order ≥ 0.9 (the P1/P2 methodology — exact equality is
  not expected across integrators).
- **(b) measured CFL envelope vs the P0 pin**: Δt sweep across
  [0.6, 1.4]× the DISCRETE undrained pencil on the toy-matched column:
  stable below ~1.0×, unstable above (boundary within ±25 % of 1.0×,
  the e72 spread); the DRAINED pencil demonstrably optimistic (a run at
  0.9× drained pencil ≫ undrained pencil diverges).
- **(c) advisory gate**: `criticalTimeStep()` WITH the overlay returns the
  undrained pencil (ratio to the no-overlay drained value in the direction
  and magnitude of the per-cell material-formula band — loose band, the
  point is direction+magnitude); naive-advisory demonstration = the (b)
  divergence at 0.9× drained. Augmentation twin-check vs a python oracle
  where extractable — else the (b) boundary IS the oracle.
- **(d) incumbent head-to-head**: same ZS84 column, upstream
  `CentralDifference` + `FourNodeQuadUP` vs the overlay lane: accuracy vs
  the implicit reference at matched Δt + per-step wall-clock (report,
  no hard perf gate — numbers to the PR) + **S→0 demo**: Kf ↑ ×1e3
  (near-incompressible): quadUP's coupled route degrades/diverges while
  the overlay lane (stab on) survives — pin the SYMPTOM, whatever it
  measures to be, loudly in the output.
- **(e) Richardson quirk pin**: honest-p LadrunoUP + upstream
  `CentralDifference` on the same column → unstable/garbage p (leapfrogged
  diffusion): loud EXPECTED-BAD gate (the run must NOT track the
  reference — if it ever starts passing, the quirks row is stale) +
  LEDGER_quirks row.
- **(f) `-subcycle` sweep under CD**: N ∈ {1, 2, 5, 10}: all stable, error
  vs N=1 grows ~monotonically (slope band 0.8–1.6 per E7.3a's ~N^1.2),
  `-subcycle auto` resolves N within ±1 of the θ-formula hand-count.
- **(g) energy closure**: ADR-69 EnergyBalanceRecorder on the (a) run:
  |ERR| ≤ the ADR-69 battery bound (read its gates for the number; else
  ≤ 2 %) — proves overlay work rides ULW.
- Runner discipline: python -S bootstrap copied from the P2 battery; e76
  protocol constants where reused; ZS84 config imported/adapted from the
  ADR-71 analytic battery.

### 3.5 Panel charter (Opus ×3)

1. **advisory-correctness critic (load-bearing):** the ΔK_e augmentation
   algebra (is Q S⁻¹ Qᵀ the right element-local undrained condensation for
   the pencil? sign/layout/DOF-map), the shared-kernel touch (SMS pencil,
   -lump variants, damped branch, useTangent), cache invalidation
   (moduli/removal), the S_e-singular guard, cost (per-call augmentation on
   large regions).
2. **lane critic:** `-fsL zero` reachability of every advance path (fs1
   hook, catch-up, SM_STEADY, P2-driver refusal completeness), the parser
   floor bypass, serialization of FSL_ZERO, the divergence advisory
   honesty.
3. **battery critic:** gate design vs the frozen e72/e74 numbers,
   ZS84-config fidelity, the incumbent leg's fairness (is upstream
   CD+quadUP configured at ITS best?), Richardson-pin correctness.

### 3.6 Registration / ledgers

No new classTag. Touched files (all fork-owned): overlay .h/.cpp +
OPS_ parser, LadrunoStaggeredDriver.cpp (refusal), CriticalTimeStep.cpp/.h
(augmentation hook). LEDGER_implementations 33022 → P3;
LEDGER_quirks rows: Richardson (e), energy-ULW accounting (3.3),
`-fsL zero` quasi-static divergence; banner line amended; pipelined-fluid
design note (§3.4 item 2) stays a documented note — no code.

## P3b — `-fluidUpdate explicit` + SMS composability (PINNED at 3b.A, 2026-07-18)

The sketch is superseded by these pins. E7.4 GOVERNS the lane (frozen
constants block): boundary = **1.32×** the discrete undrained pencil on both
benchmark soils; material factor ~1.85× conservative (documentation-only);
diffusion CFL slack **7.2e3×** at realistic k̄; **L not needed** (no
iteration — the reference rule is inert). The advisory pencil (P3 §3.2
augmentation) applies UNCHANGED — the measured explicit-fluid boundary sits
~32 % ABOVE it, so the shipped advisory is conservative for this lane too.

| WP | Owner | Deliverable |
|---|---|---|
| 3b.A | MAIN | these pins |
| 3b.B | OPUS | explicit advance branch + dual-CFL advisory + SMS Kadd wiring + warning retirement |
| 3b.C | OPUS | battery `tests/test_ladruno_overlay_explicit_fluid.py` (gates a–g below) |
| 3b.D | OPUS ×3 | panel (stability/advisory critic load-bearing; lane-reachability; battery) |
| 3b.E | MAIN | build, battery, adjudication, guide/quirks/banner/ledgers, ADR §7/§12, PR |

### 3b.1 Command semantics (FROZEN)

- Parser: `-fluidUpdate explicit` → `FU_EXPLICIT` (the fatal-NYI branch is
  removed; `fluidUpdate_` already serializes — data(20) — so NO wire change).
  The hook's FU_EXPLICIT abort is removed with it.
- **`-fsL` is INERT under FU_EXPLICIT** (no L anywhere in the update; no
  iteration). A non-default `-fsL` alongside `-fluidUpdate explicit` gets a
  one-time notice (not a fatal — model migration must not require flag
  surgery). The oedometric floor / FSL_ZERO advisory logic is untouched (it
  fires only on lanes that use L).
- **`-staticMode hold|steady` stay legal.** SM_STEADY keeps its CG solve of
  H·p = f_seep at static commits — the matrix-free claim is about the
  TRANSIENT march (the gravity-init recipe must keep working; a CG solve per
  static commit is not a march cost). Same for `-pInit steady` (setup-time).
  Documented in the guide.
- **`-subcycle auto|$N` unchanged** (window accumulation already ships; the
  fluid step is dtAccum_). The window is additionally bounded by the
  diffusion CFL — N·Δt ≤ Δt_diff — which the advisory prints; at realistic
  k̄ the slack is ~7e3× so it never binds in practice.
- **`-stab` is march-inert under FU_EXPLICIT** (pinned semantics, toy-exact):
  the lumped diagonal is the ROW-SUM of S* and H̃·1 = 0 exactly, so
  rowsum(S*) = rowsum(S_phys) — the stab matrix cannot enter the explicit
  march. This is the E7.4 oracle's own semantics (`slump = S.sum(axis=1)` on
  the stabilized S). `-stab` still affects the ADVISORY S*_e (measured
  stab-invariant, P3 §12 item 2) and nothing else on this lane. Documented,
  no notice (harmless in both directions).
- **P2 driver refusal already shipped** (FU_EXPLICIT excluded from the driven
  set — LadrunoStaggeredDriver.cpp:118); battery re-pins it as a fatal when
  the driven set comes up empty.
- `-onRemoval keep|drain` unchanged (drain's rebuildHFseep touches H/f_seep
  only; the storage lump is S-side and unaffected — see 3b.2 cache rule).

### 3b.2 The explicit advance (FROZEN — toy `cd_march(fluid="explicit")` verbatim)

- **Storage lump `sLump_`** (size N): per-row sum of the assembled S* CSR
  rows (= lumped physical S since H̃ annihilates constants). Built once when
  assembleStorageAndL completes (snapshot) and rebuilt by
  rebuildModuliCaches (moduli change S). NOT touched by removal rescan
  (fluid measure never shrinks — P1 pin). Guard at build: any row-sum ≤ 0 →
  loud fatal naming the node (the toy's `slump.min() > 0` assert).
- **Update** (inside `advanceTrial` as the FU_EXPLICIT branch — commitFluid /
  revertFluid / record / hook plumbing untouched):
  `pTrial_[i] = pCommitted_[i] + (−dt·(H·p₀)_i − QtDu_i + dt·fseep_i) / sLump_[i]`
  on FREE rows; drained rows hold their prescribed committed value (the
  implicit lane's convention). H·p₀ via the existing applyOp(0,0,1,·);
  QtDu via the existing assembleQtDu (window semantics identical). **NO
  solveFluid call — no CG, no factorization in the march.**
- **The branch always advances from COMMITTED state and ignores
  firstOfStep** (there is no reference increment; dpCommitted_ keeps being
  maintained by commitFluid for serialization compatibility but is never
  read). Consequence, pinned deliberately: a repeated advanceTrial within
  one step is IDEMPOTENT (re-derives the same pTrial_), not a double-march —
  defensive, since no shipped path calls it twice (driver refuses the lane).
- relChange_ still computed (harmless; nothing consumes it on this lane).
- Serialization: zero format changes (fluidUpdate_ travels; sLump_ is
  derived state, rebuilt at snapshot).

### 3b.3 Dual-CFL advisory (FROZEN)

- Overlay accessor `double explicitFluidDiffusionDt() const`: the
  per-fluid-advance forward-Euler diffusion bound
  **Δt_diff = minEdge_² / (2·ndm_·maxCv_)** (the toy's dt_diff formula;
  conservative — min-edge with max-c_v). Returns −1 when not FU_EXPLICIT or
  snapshot not ready. Zero new state (minEdge_/maxCv_ already exist for
  `-subcycle auto`).
- `CTSResult` gains `bool explicitFluid` + `double fluid_diffusion_dt`
  (defaults false / +inf). `computeCriticalTimeStep` (already collecting
  33022 overlays once per call): FU_EXPLICIT overlays with a finite
  diffusion bound set explicitFluid and min-fold fluid_diffusion_dt; the
  bound is then **min-folded into BOTH damped_dt and undamped_dt** (if
  diffusion ever governs — absurd k̄ — the returned advisory is honest;
  governing tags stay the solid element's, the report clarifies). MPI
  reduction: min-reduce fluid_diffusion_dt, OR-reduce explicitFluid
  alongside the existing block.
- Reports (CentralDifferenceLadruno ×2 sites + ExplicitBathe sites — the
  overlayAugmented print blocks): when explicitFluid, append one line:
  explicit-fluid diffusion limit Δt_diff, the slack factor
  Δt_diff/undamped_dt, and the note that `-subcycle N` must keep
  N·Δt ≤ Δt_diff. Governing = min is already reflected in the folded
  numbers.

### 3b.4 SMS composability fix (FROZEN — the Kadd seam wired)

- BOTH builders (`buildMassScaling` AND `buildMassScalingConsistent`)
  collect 33022 overlays once per build (the computeCriticalTimeStep idiom)
  and, per element, sum `getUndrainedAugmentation` blocks into a dense
  Matrix passed to `elementCriticalDt(…, &Kaug)` — SMS then sizes against
  the UNDRAINED per-element pencil. The closed-form scale s = T² + 2Tc is
  UNCHANGED and remains exact: ΔK_e is mass-independent and
  state-independent, so λ_max still scales as 1/s under mass injection —
  the augmented dt_e is the correct sizing input with zero formula changes.
  Size-mismatch (exotic ndf) → skip augmentation for that element with a
  one-time loud note (the CTS behavior, mirrored); self-reported bounds
  keep winning UNCORRECTED (the CTS scan already warns loudly for that
  combo; SMS's self-report skip path is untouched).
- **The blanket `warnIfOverlayPresentSMS` warning is RETIRED** (function
  removed; both call sites replaced by the overlay collection). Replaced
  by: (i) a one-time INFO line when sizing augmented ≥ 1 element (`SMS
  sizing priced the UNDRAINED pencil for N overlay-owned elements
  (ADR-73 P3b)`) — battery-greppable honesty; (ii) the residual
  not-augmented cases keep their existing LOUD per-case advisories, which
  are exactly getUndrainedAugmentation's refusals (singular S_e, size
  mismatch) — never silently optimistic; (iii) a defensive loud warning if
  any 33022 overlay is present with `!snapshotReady()` at sizing time
  (unreachable in practice — the snapshot fires in setDomain at pattern
  creation, long before analysis setup; a FAILED snapshot aborts every
  commit anyway).
- **No LADRUNO/vanilla boundary change**: all touched files are fork-owned
  (LadrunoMassScaling.h, CriticalTimeStep.{h,cpp}, overlay .h/.cpp, OPS_
  parser, CentralDifferenceLadruno.cpp, ExplicitBathe.cpp) — no
  vanilla-ledger rows.

### 3b.5 Battery `tests/test_ladruno_overlay_explicit_fluid.py` (gates)

Runner discipline copied VERBATIM from `tests/test_ladruno_overlay_explicit.py`
(machine python `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\
python.exe`, `-S` + dist\bin pyd eviction + `assert opensees.__file__`,
collection-safe `pytest.skip(allow_module_level=True)`, never module-level
sys.exit). e74 config constants where reused.

- **(a) diffusion-slack sweep**: realistic k̄ sweep (1e-7 … 1e-3 m/s class)
  on the e74-matched column: the advisory's Δt_diff matches the
  h²/(2d·c_v) hand formula, and the slack vs the solid undrained pencil is
  ≥ 1e2× at every realistic point (toy: 7.2e3× at k̄ = 1e-7-class; the gate
  band is loose — the POINT is orders of slack, direction pinned).
- **(b) coupled CFL envelope vs the 1.32× pin**: Δt bisection sweep on the
  e74-matched DRAINING column (NOT the zero-k pathology — fixed-horizon
  honesty per §12 P3 item 3): measured boundary / discrete undrained pencil
  ∈ 1.32× ± 25 % (the e72/e74 spread convention → [0.99, 1.65]); a run at
  0.8× the pencil marches the full horizon stable.
- **(c) equivalence vs the implicit-fluid lane**: same column, FU_EXPLICIT
  vs FU_IMPLICIT + `-fsL zero` at matched Δt (both O(Δt) splits — the
  family's mutual-convergence methodology): inter-lane diff shrinks with
  Δt at observed order ≥ 0.9.
- **(d) SMS composability (THE gate — the inverted P3 quirk)**: overlay
  column where the certified dtTarget under DRAINED pricing provably blows
  up (the P3 quirks scenario): with the fix, (i) the SMS report prices the
  undrained pencil (scaled-element count / added mass reflects the ~26×
  class factor vs a no-overlay control), (ii) a march at the certified
  dtTarget is STABLE over the fixed horizon, (iii) the blanket UNSUPPORTED
  warning no longer prints, the INFO line does. BOTH builders exercised
  (lumped CentralDifferenceSMS + consistent variant).
- **(e) MP — TRANSPOSED at 3b.A (adjudication recorded in §12)**: the
  overlay is per-process serial v1 (P4 §12: a partition without the
  pattern serves no channels); there is NO halo exchange in-tree to smoke.
  The gate becomes: (i) the update is demonstrably partition-LOCAL (axpy on
  overlay-owned CSR state, no Domain-global reduction — code-audited by the
  panel, not run under MPI); (ii) the halo DESIGN for a future
  partitioned overlay is documented in the ADR §8 amendment (QᵀΔu
  contributions summed at shared region nodes + p halo — the explicit
  solid's own nearest-neighbor pattern). NO MP run claimed; the §8
  "dissolution" stays design-demonstrated. DO NOT overclaim in guide/PR.
- **(f) P4 recorder transparency**: the recorder battery's CSV-twin gate
  re-run under FU_EXPLICIT (`recorder ladruno -overlay` HDF5 DATA ==
  `-record` CSV rows, exact; region/drained topology rows intact).
- **(g) zero-drainage secular pumping under fully-explicit fluid —
  MEASURED**: the P3 quirks scenario (k̄ ~ 1e-11, undamped) re-run with
  FU_EXPLICIT at 0.4×/0.5×/1.0× pencil: steps-to-blowup recorded and the
  quirks row amended either way (same pathology class expected — frozen
  forces + no dissipation — but the fluid-side dynamics differ; measure,
  never assume).
- Smokes: driver refusal fatal on an FU_EXPLICIT-only domain; `-fsL
  classic` + explicit → one-time inert notice, march unchanged vs default;
  SM_STEADY gravity-init recipe under FU_EXPLICIT then transient march;
  DB round-trip (send/recv) of an FU_EXPLICIT overlay mid-march.
- No-regression net: P1 physics 6/6, framework 9/9, driver 12/12, explicit
  8/8, recorder 9/9 re-run green.

### 3b.6 Panel charter (Opus ×3)

1. **stability/advisory critic (load-bearing)**: explicit update vs the toy
   oracle (signs, dt placement, window dtAccum, drained-row handling),
   row-sum lumping semantics (H̃ annihilation — is rowsum(S*) really
   stab-free in the C++ assembly given per-layer scatter?), the dual-CFL
   fold (does min-folding damped_dt break any existing consumer?), SMS
   Kadd algebra (scale-invariance of s = T²+2Tc under augmentation;
   betaK-damped form; consistent-builder M_bar interplay), warning
   retirement honesty (enumerate every config that still prices drained —
   is each one loud?).
2. **lane-reachability critic**: every FU_EXPLICIT path — hook (SM_MARCH /
   HOLD / STEADY), subcycle window + removal-forced sync, drain rebuild,
   restart (recvSelf → snapshot → sLump_), getCopy, parser round-trip,
   driver refusal completeness, catchUpPendingWindow, external-stepping
   latch defensiveness, Print output, the P4 channels under the lane.
3. **battery critic**: gate design vs the frozen e74 numbers (is the ±25 %
   band anchored to the right pencil?), fixed-horizon honesty on (b)/(g),
   SMS gate (does (d) actually refute the OLD pathology — would it fail on
   pre-fix code?), recorder-twin fidelity, runner discipline.

### 3b.7 Close-out / ledgers

No new classTag; no vanilla touches expected. Guide §7: SMS caveat LIFTED
(replaced by the composability statement + INFO line), explicit-fluid
recipe + inherited advisories (zero-drainage, incubation, dual-CFL);
LEDGER_quirks: SMS row amended (fixed at P3b, warning retired), secular-
pumping row amended per gate (g), any new rows; banner 33022 line amended
(`-fluidUpdate explicit` + SMS composability) via banner_features.txt →
patch_banner.py; LEDGER_implementations 33022 row appended (P3b, stays
shipped); ADR §7 P3b row closed + §12 P3b entry (record the MP gate
transposition + every refutation — adjudication duty); **pipelined-fluid
design note (§3.4 item 2) stays a documented note — no code** (re-affirmed
from P3 §3.6). apeGmsh emitter/contract row = companion-repo follow-up,
PR-body note only.

## P4 — ecosystem (PINNED at 4.A, 2026-07-18)

The sketch is superseded by these pins. Scope: overlay p-field recorder
channels (LadrunoRecorder + LadrunoMonitorRecorder), the user guide, family
close-out (banner/ledgers/ADR). NO physics/driver/advisory code changes.

| WP | Owner | Deliverable |
|---|---|---|
| 4.A | MAIN | these pins committed |
| 4.B | OPUS | recorder channels: overlay getters + `Ladruno_OverlayResults.{h,cpp}` + LadrunoRecorder `-overlay` + Monitor `-overlay` |
| 4.C | OPUS | battery `tests/test_ladruno_overlay_recorder.py` (gates a–g below) |
| 4.D | OPUS | `Ladruno_implementation/LadrunoPorousOverlay_guide.md` per the outline below |
| 4.E | OPUS panel | adversarial gate: recorder-wiring critic (core-adjacent) + guide documentation-fidelity critic |
| 4.F | MAIN | build, batteries, banner/ledgers, ADR §7 P4 close + §12 entry, PR |

### 4.1 Overlay accessor seam (FROZEN — const getters only, zero physics)

Added to `LadrunoPorousOverlay.h` public block (all trivial const reads of
existing private state; nothing mutates, nothing serializes):

- `const std::vector<int>&    getRegionNodeTags() const;`  (region-node order —
  THE canonical ID order for every channel below)
- `const std::vector<double>& getCommittedP() const;`      (aligned to that order)
- `double pAtNodeTag(int nodeTag) const;`  (via `nodeTagToIndex_`; returns 0.0
  AND sets no error for a non-region tag — callers must pre-validate; the
  Monitor parser fatals at construction on a non-region tag)
- `const std::vector<int>&    getDrainedNodeTags() const;` (as-declared list)
- `bool snapshotReady() const;` (recorder registration must not fire the lazy
  snapshot — see timing pin 4.3)

### 4.2 LadrunoRecorder channels (FROZEN)

- **Command surface:** `recorder ladruno $file ... -overlay <$tag1 $tag2 ...>`
  — bare `-overlay` = every 33022 pattern in the domain; explicit tags =
  exactly those (unknown/non-33022 tag = parser-time fatal, the
  unknown-flag-FATAL house rule). Repeatable/combinable with all existing
  flags (`-N/-E/-G/-kind/-envelope/...`). Duplicate tags deduped with a
  notice (panel F-3). **Amended at 4.F (panel adjudication):** BOTH forms
  fail-fast at PARSE time when unsatisfiable (bare with no 33022 in the
  domain; unknown explicit tag) — the original "fatal at writeModel time"
  wording is unobservable in practice because `Domain::commit` swallows the
  recorder's `record()` return; writeModel re-checks remain as the backstop.
  Consequence (documented, accepted): the overlay pattern must exist BEFORE
  the recorder command — recorder-before-pattern ordering is rejected loudly.
- **Discovery idiom:** scan `domain->getLoadPatterns()` for
  `getClassTag() == PATTERN_TAG_LadrunoPorousOverlay` (the P2 driver /
  Domain.cpp:2248 idiom).
- **Source:** new `OverlayPressureSource : ladruno::ResultSource` in
  `SRC/recorder/Ladruno_OverlayResults.{h,cpp}` (modeled on
  `Ladruno_DomainResults`): `ids()` = `getRegionNodeTags()`; schema
  `name = overlayPressure_<patternTag>`, `components_csv = "p"`,
  `num_components = 1`; `evaluate()` copies `getCommittedP()`;
  `requiresPartitionReduction() = false` (overlay is per-process serial v1).
- **Channel plumbing:** dedicated `overlay_channels` vector +
  `initOverlaySources()` (called from `writeModel()` after the domain-source
  init, LadrunoRecorder.cpp:640 block) + `recordResultsOnOverlays()` (called
  in `record()` after `recordResultsOnDomain()`, :395) — sinks constructed
  with `ResultFamily::OnNodes` so datasets land in
  `RESULTS/ON_NODES/overlayPressure_<tag>/{ID,DATA,STEP,TIME}` with the
  standard self-describing attrs, written by the UNTOUCHED generic
  `StreamingSink`/`EnvelopeSink`. `-envelope`/`-kind` honored exactly like
  other channels (envelope |p|max is the liquefaction use case). Cleanup in
  `clearSources()`/`finalizeAllSinks()` mirrored.
- **Write-once topology rows:** new `writeModelOverlays()` in the
  `writeModel()` block (:618–627), mirroring `writeModelSets()` (:953):
  `MODEL/OVERLAYS/OVERLAY_<tag>/{REGION_NODES, DRAINED_NODES}` int datasets +
  `TAG` attr. Region nodes are ordinary solid nodes already in `MODEL/NODES`
  — no node writing changes. NO `NDIR/NUM_GP/COLUMN_MAP` (element-only
  conventions; this is a nodal scalar).
- **Schema stance:** strictly ADDITIVE under FORMAT_VERSION = 1 (new optional
  groups only; no existing dataset/attr changes; the ladruno_schema_v1
  validator must stay green on non-overlay files unchanged — extend
  `ladruno_format.py` ONLY if it hard-fails on the new groups, warn-tolerant
  preferred, the QUADRATURE-tolerance precedent).
- **Serialization:** the new overlay-tags config field travels in the
  recorder's sendSelf/recvSelf next to the existing request lists (pattern at
  :2085/:2180). MP semantics of the overlay itself remain out of scope (v1
  per-process; documented).

### 4.3 Timing pin (VERIFIED at 4.A — the ordering fact the channels rest on)

`Domain::commit()` fires the overlay hook (`onDomainCommit`, Domain.cpp:2260)
BEFORE the recorder loop (`theRecorders[i]->record`, Domain.cpp:2285) — so
`getCommittedP()` at record time is the freshly committed pⁿ⁺¹ paired with the
step's committed uⁿ⁺¹. Under `-subcycle N>1` the recorder therefore reads the
last SYNCED p between windows (matches the physics — those ARE the forces the
skeleton feels); gate (e) pins the cadence. Under the P2 driver the latched
commit syncs the converged iterate before recorders run — same invariant, no
special case. Registration/writeModel must NOT trigger the overlay's lazy
snapshot (`snapshotReady()` guard: overlays not yet snapshotted at writeModel
time get their topology rows written lazily at first record — or recorder
construction after `analyze` begins is simply documented; agent picks the
simpler honest behavior and the battery pins it).

### 4.4 LadrunoMonitorRecorder channel (FROZEN)

- `recorder Monitor -overlay $tag <-nodes {$n1 ...}> -sink $file <-every $N>
  <-hz $hz>` — overlay mode is EXCLUSIVE of `-node/-region/-dof/-resp`
  (parser fatal on mixing; SWMR columns are frozen at open, keep the modes
  disjoint; `-nodes` is the permitted region-node SUBSET selector — the
  earlier inclusion of `-nodes` in this exclusivity list contradicted the
  grammar line and is corrected here, panel F-3-schema). `-nodes` subset must be region nodes (fatal otherwise); default =
  all region nodes. One overlay per Monitor instance (spawn two recorders for
  two overlays).
- Columns `overlay<tag>.p.node<n>` (node-major, the existing label
  convention), values via `pAtNodeTag()` (committed p — the Monitor records
  post-commit like today's channels). New dataFlag branch in
  `respToDataFlag`/`record()` (LadrunoMonitorRecorder.cpp:43/:204).
  sendSelf/recvSelf stay stubs (sequential-only, as today).

### 4.5 Battery `tests/test_ladruno_overlay_recorder.py` (Zone-B, python -S)

Runner bootstrap copied verbatim from `test_ladruno_overlay_physics.py`
(site-packages re-add, dist/bin pyd eviction + `assert opensees.__file__`,
module-level `pytest.skip(allow_module_level=True)` when unbuilt, h5py reads
with `HDF5_USE_FILE_LOCKING=FALSE`). Gates:

- **(a) CSV-twin gate (THE gate):** one Terzaghi column run with BOTH
  `-record` CSV and `recorder ladruno -overlay`: HDF5 `DATA[k]` ==
  CSV row k to ≤1e-12, `ID` == the CSV header node order ==
  `REGION_NODES`; row count == commit count. Proves the 4.3 ordering AND the
  channel end-to-end.
- **(b) topology rows:** `MODEL/OVERLAYS/OVERLAY_<tag>/{REGION_NODES,
  DRAINED_NODES}` match the model as built; drained rows identically 0 in
  DATA.
- **(c) multi-overlay:** two disjoint overlays → two groups + two topology
  rows, values segregated (cross-check vs per-overlay CSVs).
- **(d) fatals:** unknown tag after `-overlay`; bare `-overlay` with no
  overlay in the domain; Monitor `-overlay` mixed with `-node`; Monitor
  `-nodes` containing a non-region node.
- **(e) subcycle cadence:** `-subcycle 4` run — recorded p piecewise-constant
  between syncs, changes exactly at sync commits.
- **(f) Monitor twin:** Monitor SWMR file (post-run h5py read): columns ==
  `overlay<tag>.p.node<n>`, FRAMES rows == CSV rows at matched commits;
  `-nodes` subset selects the right columns.
- **(g) no-regression:** the existing recorder regression battery
  (`Ladruno_scripts/ladruno_recorder_tests/` harness) green; P1 physics
  battery smoke (gate 1 anchor) green.
- **(h) envelope:** `-kind`/`-envelope` overlay channel == hand-computed
  running max over the CSV (skip-with-note if envelope wiring turns out
  family-gated — adjudicate, never silently drop).

### 4.6 Guide `Ladruno_implementation/LadrunoPorousOverlay_guide.md` (outline FROZEN)

Every claim must trace to shipped behavior, an ADR §12 entry, or a
LEDGER_quirks row — NO invented numbers (the fidelity critic enforces).
Sections:

1. **What / when** — division-of-labor table vs LadrunoUP (ADR §2.6/§6:
   monolithic = implicit statics/TH/B1–B5; overlay = discontinuity/explicit/
   companion lane).
2. **Quick start** — minimal Terzaghi column (Tcl + openseespy).
3. **Command reference** — the full shipped surface (§4.1 as-built incl.
   `-fsL classic|oedometric|$scale|zero`, `-staticMode`, `-subcycle`,
   `-layer`, `-pInit`, `-onRemoval`, `-record`, fatal-NYI `-fluidUpdate
   explicit`), factor/TimeSeries ignored-by-design.
4. **Initialization recipes** — gravity staging with `-staticMode steady|hold`
   (the ⟨A-3⟩ static-time trap: static "time" = load factor, never march on
   Δλ), `-pInit steady|hydrostatic|list`; mixture-ρ gravity convention + the
   staggered-twin BODY-FORCE trap (quad `b` = FORCE/VOLUME = ρ_mix·accel vs
   `LadrunoUP -body` = ACCELERATION — LEDGER_quirks 2026-07-18 row).
5. **Layered profiles** — `-layer` examples; ONE overlay per hydraulically
   connected water body (per-layer k/poro/moduli inside one overlay; separate
   overlays = separate aquifers only).
6. **The iterated driver (P2)** — `LadrunoStaggeredAnalyze` usage + `-stats`;
   integrator WHITELIST Newmark/HHT/GeneralizedAlpha (TRBDF2 two-step-history
   / explicit no-op revertToLastStep = fatal); stage-flip moduli recipe
   (`updateMaterialStage` CANNOT reach the overlay — re-set via
   `parameter ... loadPattern E|nu|layerE|layerNu`); the raise-vs-negative-int
   failure-surface split in openseespy.
7. **The explicit lane (P3)** — `-fsL zero` under CentralDifferenceLadruno at
   Δt ≤ 0.5× the overlay-aware `criticalTimeStep()`; the advisory prints the
   undrained pencil (per-element factor legitimately EXCEEDS the material
   formula — documentation-only, ~1.85× conservative); **zero-drainage
   secular-pumping advisory** (undamped zero-k column has NO asymptotically
   stable Δt — the 0.5× pin is HORIZON-relative there; any physical drainage
   or damping absorbs it — §12 P3 item 3); **SMS + overlay UNSUPPORTED until
   P3b** (sizing prices drained pencils; loud warning ships); quasi-static
   `-fsL zero` diverges by design (~4 steps); `-subcycle` under CD.
8. **Why not CD + u-p elements** — the Richardson quirk incl. the ~1300-step
   INCUBATION warning (tracks the reference then dies; a 400-step validation
   "passes"); upstream CD+quadUP head-to-head numbers (§12 P3 item 4, cited
   not re-measured).
9. **Recording p** — the P4 channels (`-overlay` recorder + Monitor) +
   `-record` CSV; energy accounting note (+Q·p work rides ADR-69 ULW,
   closure holds, attribution merged).
10. **Element removal** — `-onRemoval keep|drain $kF`, fluid survives
    `remove element`, forced sync under `-subcycle`.
11. **Quirks / troubleshooting** — table of the LEDGER_quirks rows with
    one-line symptoms + pointers.

### 4.7 Registration / ledgers / close-out

- No new classTag. New files: `Ladruno_OverlayResults.{h,cpp}` (stamped, into
  `stamp_headers.py` GLOBS + recorder CMakeLists), the battery, the guide.
  Touched fork-owned: overlay .h/.cpp (getters), LadrunoRecorder.{h,cpp},
  LadrunoMonitorRecorder.{h,cpp}. NO vanilla touches expected (recorder
  dispatch already routes `ladruno`/`Monitor`); if one appears, ledger row +
  `// Ladruno` mark.
- Banner 33022 line amended: `+ p-field recorder channels (recorder ladruno
  -overlay / Monitor -overlay)`; `patch_banner.py` re-run.
- LEDGER_implementations 33022 row → **shipped** (P4 files/PR appended).
- ADR §7 P4 row closed + §12 P4 entry; any new quirks → LEDGER_quirks.
- apeGmsh contract row = companion-repo follow-up, noted in the PR body only.

---

## Frozen measured constants (appended by 0.C — P0 landed 2026-07-14)

Source: `Ladruno_implementation/adr71_meshless_p_spike/staggered_pins_e7.py`
(toy imported as a frozen library; predictions printed before measurement;
numbers in `e7_summary.json`, log in `e7_run_full.log`). Refutations are
flagged inline — **three plan predictions changed sign; P1/P3 defaults below
supersede the corresponding sketch text above.**

### E7.1 — v1 sequencing ⚠ PREDICTION REFUTED → P1 PIN CHANGED

- Measured: BOTH plan flavors are O(1)-floored with the fixed-point L form
  (plain 0.61 relL2 at every Δt for either L; resolve+classic 0.41;
  resolve+oed 5.7e-4 = special-case L≈Schur-compliance match, column only).
  Mechanism: single-pass fixed-point L double-counts compliance (storage
  S*+L+C instead of S*+C).
- **FROZEN v1 fluid advance = fixed-stress RATE form**
  `(S*+L)Δp + ΔtH p₁ = −QᵀΔu + L·Δp_prev` (one stored vector: previous
  committed Δp). Plain-`analyze` flow, no driver, NO second solid solve.
  Measured O(Δt): late-window orders +1.00…+1.04, both L, both regimes;
  L=oed near-deadbeat (5.6e-6). Step-load initial layer = one-step O(1)
  (harmless; full-history norm reads O(√Δt) — document in guide).
- Scope: QUASI-STATIC/implicit lane only — see E7.2 carve-out.

### E7.2 — dynamic stability envelope ✔ CONFIRMED + accuracy rider

- Fixed-point-L implicit fluid at commit, CD solid: boundary = **0.987×**
  drained pencil (classic), 0.934× (oed), 0.780× (half-oed), **1.000× the
  undrained pencil for L=0**. ⟨A-7⟩ hedge resolved: drained-CFL stability is
  real.
- ⚠ Rider: that lane carries an O(1) diffusion-rate artifact (0.673 relL2
  vs exact; Δt-independent) — stability true, accuracy false.
- ⚠ The E7.1 rate form is UNSTABLE under CD (every tested Δt ≥ 0.5×
  undrained pencil) — does not transfer to the dynamic lane.
- **FROZEN P3 default: L=0, Δt ≤ 0.5× undrained pencil** (consistent +
  stable; 1.5e-3 vs exact reference). Drained-CFL fixed-point-L mode =
  opt-in, must be documented accuracy-degraded.

### E7.3a — `-subcycle auto` ✔ CONFIRMED (L=0 lane)

- err ~ N^1.2, all N ∈ {1…50} stable; error doubles at N ≈ 1.8 →
  **θ = N·Δt/(h²/c_v) = 0.089**.
- **FROZEN: N_auto = max(1, floor(0.09 · h²/(c_v·Δt)))**, h = min element
  size, c_v = max over layers (conservative).

### E7.3b — `-substep` ⚠ ACCURACY CLAIM REFUTED

- M ∈ {1,2,5,10} at Tv=0.05/step: 0.2452/0.2489/0.2521/0.2534 — flat.
  Coarse-Δt splitting error dominates; fluid sub-resolution buys nothing.
- **FROZEN: `-substep` demoted to plumbing-only (or dropped at P1's
  discretion); the guide must not promise early-time transient recovery.**
  Implementation note: if kept, the rate-form reference is the COARSE
  increment split /M — chaining it per sub-step explodes at M ≥ 5
  (measured).

### E7.4 — explicit-fluid dual CFL ✔ CLASS CONFIRMED / ⚠ ±20 % factor gate REFUTED (safe direction)

- Boundary tracks the discrete undrained PENCIL: dt_cr = **1.32×** pencil,
  both benchmark soils (5.567e-4 / 5.531e-4 s vs pencils 4.206e-4 /
  4.196e-4 s). Material formula √(1+K_f/(n·M_oed)) (21.4 / 12.8) is ~1.85×
  conservative as a Δt_cr predictor (measured effective factors 11.6 /
  6.97). Diffusion CFL slack 7.2e3× at realistic k̄. L=0 confirmed (no
  iteration).
- **FROZEN: `-fluidUpdate explicit` Δt advisory = the discrete undrained
  pencil (≈25 % margin vs measured); material factor formula =
  documentation-only (~2× conservative).**

### E7.5 — removal under inertia ✔ (revised metric) / original commit-jump formulation refuted benignly

- Commit-jump vanishes ~dt^2.6 → **no 1/Δt impulse artifact** (the
  "dt-independent commit jump" expectation was wrong benignly: the
  redistribution spreads over the physical transient). The dt-independent
  invariant is the fixed-window jump (20·dt₀): spread 1.47 % over
  dt-halving ×4 (commit-jump scaling exponent +2.63).
- `-onRemoval drain` kFactor {1,10,100}: Tv90 = 1.004 / 0.938 / 0.933 —
  strictly monotone ✔.
- **FROZEN: P1/P3 removal gates use the fixed-window jump metric, never the
  single-commit jump.**

### E7.6 — `-fsL` default pair ⟨A-2⟩ ✔ CONFIRMED (footing gain smaller)

- classic α²/K_dr: **zero divergences, zero maxit-hits** on column+footing ×
  compressible+near-undrained (means 11.2/4.3/11.3/4.3, max ≤ 18).
- oedometric: also clean everywhere (means 3.3/2.8/3.3/2.8, max ≤ 10);
  speedup **3.4× column / 1.6× footing** (prediction said ~3×; footing is
  1.6×).
- 0.5×oed: maxit-hit (500) — the E6 cliff, reproduced.
- **FROZEN: default `-fsL` classic (proof-backed, divergence-free);
  `-fsL oed` documented opt-in (1.6–3.4×); hard floor: never below oed.**

