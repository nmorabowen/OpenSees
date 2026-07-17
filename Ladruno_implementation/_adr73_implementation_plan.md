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
   `‖p_k − p_{k−1}‖₂ / (‖p_k‖₂ + pScale)` (toy `dcon` verbatim, computed
   inside `advanceTrial`, exposed as `lastAdvanceRelChange()`).
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

## P3 / P3b — explicit lanes (sketch)

- P3: CD + overlay vs implicit monolithic two-leg on B2/ZS84; incumbent
  head-to-head (upstream CentralDifference + FourNodeQuadUP) incl. S→0 demo;
  E7.2 envelope confirmed in C++; overlay-aware Δt_cr advisory — overlay
  exposes per-element undrained factor; integration point into the ADR-65
  advisory machinery pinned at phase start (integrator queries Domain for
  overlay presence); `-subcycle auto` gate; ADR-69 energy-channel advisory.
- P3b: `-fluidUpdate explicit` per ADR §3.4 item 1 — dual-CFL gate, implicit-
  lane equivalence at matched Δt, SMS composability, MP halo-exchange smoke.

## P4 — ecosystem (sketch)

Overlay p-field recorder channels (LadrunoRecorder/Monitor topology rows),
`LadrunoPorousOverlay_guide.md` (division-of-labor table vs LadrunoUP, init
recipes, one-overlay-per-water-body, layered `-layer` examples, mixture-ρ
gravity note), banner line + `patch_banner.py`, ledger → shipped, apeGmsh
contract row (companion repo).

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

