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
  commit. Prediction (ADR §3.4): stable at the DRAINED-speed CFL (the implicit
  fluid absorbs the undrained stiffening). Sweep Δt across
  [0.5·Δt_drained, 1.2·Δt_undrained⁻¹-scaled] and locate the empirical
  boundary. THE load-bearing P3 pin.
- **E7.3 — subcycle degradation curve.** N ∈ {1, 2, 5, 10, 20, 50}: accuracy
  vs monolithic + stability. Prediction: error grows ~linearly in N·Δt while
  N·Δt ≪ the diffusion time h²/c_v; pins the `-subcycle auto` formula.
- **E7.4 — fully-explicit fluid (P3b).** Lumped S*, forward p-step.
  Predictions: diffusion CFL slack by orders (numbers in ADR §3.4); coupled
  stability at the UNDRAINED-speed CFL, factor √((M_oed+K_f/n)/M_oed) vs the
  drained pencil measured within ±20 %; L not needed (no iteration).
- **E7.5 — removal-step consistency on the dynamic path.** Crack opens
  mid-CD-march (E4 scenario + inertia): prediction — physical stress-wave
  transient only, no numerical artifact scaling with 1/Δt; `-onRemoval drain`
  kFactor sweep {1, 10, 100} monotone in drainage time.

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
    -drained {$n1 ...} <-pInit ...> <-stab ...> <-fsL auto <$s>> \
    <-onRemoval keep|drain $kF> <-fluidUpdate implicit> <-subcycle $N> \
    <-record $file <$everyN>> <-fluidBody ...> <-dynSeepage on|off>
```

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
- **Force application:** override `applyLoad(time)` → per region node
  `node->addUnbalancedLoad(+Q·p_committed slice)` with factor 1 (pattern
  TimeSeries/factor ignored — one-line notice at creation). Sign: force =
  **+Q·p** on the skeleton (compression-positive p; toy convention
  `u = K⁻¹(f + Q p)`; the sequencing contract test pins it physically).
- **Hook:** `Domain::commit()` gains one guarded block (classTag == 33022 →
  static_cast → `overlay->onDomainCommit()`), following the
  `LadrunoContactDomain` include precedent (Domain.cpp:82). Honors
  `-subcycle`: advance every Nth commit, else only re-snapshot removal state.
  Fallback if the panel kills it: applyLoad new-step detection (documented,
  not implemented speculatively).
- **Removal rescan:** at each hook firing, region elements checked against the
  Domain; newly-dead cells → Q_e dropped, H/S per `-onRemoval` (drain →
  k×kFactor in those cells, refactor/re-setup CG operator). Dead-cell GP set
  never shrinks the fluid measure under `keep`.
- **`-record $file <$everyN>`**: CSV `time, p(node1..nodeN)` — the gate-facing
  output until P4's proper channels. Region node order written once as a
  header row.
- **dynSeepage default off** (ADR-71 §12 amendment inherited).

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

## Frozen measured constants (appended by 0.C — empty until P0 lands)

*(E7.1 sequencing default + tol · E7.2 explicit envelope · E7.3 auto-N
formula · E7.4 undrained factor · E7.5 removal-jump verdict)*
