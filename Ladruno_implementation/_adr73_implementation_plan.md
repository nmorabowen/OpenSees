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
| 1.E | OPUS ×2 ∥ | battery, split: **(i) physics** — Terzaghi vs analytic (E7.1-pinned tol); two-leg vs monolithic LadrunoUP (mutual Δt-convergence, order ≥ 1); E4-in-OpenSees: real `remove element` crack, both `-onRemoval` policies, curves vs toy reference; sequencing contract test (step-load column: p(0⁺) ≈ q — kills the fluid-first p≡0 class); L-floor loud-warning test; **(ii) framework** — pattern-timing bit-constancy (forces identical across Newton iterations of one step, Static AND Newmark AND CD); factor-ignored semantics; revertToLastCommit consistency (fluid state rolls with the domain); region validation errors (unsupported cell, empty drained set, overlapping overlays on one element → fatal); pInit steady/hydrostatic vs closed form |
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
    <-subcycle auto|$N | -substep $M> \
    <-record $file <$everyN>> <-fluidBody ...> <-dynSeepage on|off>
```

- **Multirate, both directions** (adjudicated surface addition, 2026-07-13 —
  fold into the ADR §4.1 at the next amendment): `-subcycle N` = fluid
  advanced every N solid commits with the accumulated Δu (explicit lane;
  fluid slower); `-substep M` = M backward-Euler fluid sub-steps per solid
  commit with Δu linearly ramped, factored operator reused (implicit
  consolidation lane; fluid faster — recovers the early-time transient a
  single coarse BE step smears). Mutually exclusive flags; removal events
  force a sync regardless of N; E7.3 pins both curves.

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
  α_stab·H̃, L = (1/(K_dr+4G/3) per layer)·M_p via kernel addH/addS/addHtilde.
  Drained BCs by row/col elimination. **Solver pin: CG + Jacobi, rel tol
  1e-10, warm-started from the previous p** — numbering-independent, zero new
  dependencies; direct-factor upgrade only on profile evidence (ADR-40
  discipline).
- **State discipline (P2-proofed now):** `pCommitted_`, `uSnapshot_`
  (region u at last fluid advance), and THREE methods —
  `advanceTrial()` (solve fluid from current committed/trial u, refresh
  forces, do NOT commit), `commitFluid()` (pTrial→pCommitted, refresh
  uSnapshot_), `revertFluid()` (drop trial). fs1 = advanceTrial+commitFluid
  back-to-back inside the Domain::commit hook. The P2 iterated driver calls
  advanceTrial repeatedly before one commitFluid — no P2 rewrite.
  Serialization pin ⟨A-10⟩: `pCommitted_` travels in sendSelf; `uSnapshot_`
  is NOT serialized — recvSelf re-derives it from committed node
  displacements (they travel with the nodes), keeping restart exact.
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

## P2 — iterated driver + stage transport (sketch; pins at phase start)

- `LadrunoStaggeredAnalyze $n $dt <-tol $t> <-maxIter $k>` command:
  per step — analyze(1,dt) via the existing analysis, then loop
  {revertToLastCommit, overlay advanceTrial, re-analyze} until ‖Δp‖ test,
  then commit (overlay commitFluid rides the hook). Iteration telemetry
  (mean/max) printed and queryable.
- Moduli/stage transport: overlay `setParameter("updateMaterialStage"…)` or
  explicit re-set of `-moduli` via parameter path; dirties L/α-stab/CG
  operator. PDMY staged column gate (vs ADR-71 P4 monolithic reference).
- Fixed-point gate: iterated staggered ≡ monolithic LadrunoUP same-Δt ≤ 1e-6.
- B4 footing staggered CB gate.

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

