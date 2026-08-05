# ADR-71 side-spike — meshless pressure field for the Biot u-p problem

**Question** (NMB, 2026-07-13): replace the FE pressure interpolation with a
meshless (MLS) field — particles at FE nodes, element centroids, or Gauss
points — does it make sense, does it work, and is there a performance benefit?

**Method**: `meshless_p_toy.py` — the ADR-70/71 pre-code numpy-oracle idiom.
Same Biot blocks (Q/H/S/H̃), honest-p contract, quasi-static backward Euler as
the shipped LadrunoUP kernel; only the pressure basis is swapped. Credibility
anchors: FEM equal-order and FEM+`-stab auto` reproduce the shipped P0
inf-sup gate results (≥2 spurious modes / exactly 1).

## Verdicts

| pressure space | E1 spurious modes | E2 Terzaghi relL2 (Tv=.1/.5) | E3 CB | E3 sett/ref |
|---|---|---|---|---|
| FEM equal-order (shipped) | 8 | 0.12% / 0.22% | 0.086 | 0.89 |
| FEM + α-stab (shipped) | **1** | 0.25% / 0.25% | 0.015 | 1.00 |
| MLS node cloud (s=2) | 3–4 | 0.16% / 0.23% | 0.068 | 0.89 |
| MLS centroid cloud (s=2) | **1** | 4.5% / 8.9% | 0.010 | 0.95 |
| **MLS GP cloud** (s=1.6–3.2) | **158 = nP−nU floor** | 1.6% / 12.6% | 0.007¹ | **0.09–0.27 locked** |

¹ CB index blind to the GP cloud's failure mode (over-constraint, not
oscillation); its p field peaks at 4× the load with spurious suction.

1. **GP cloud (the original idea): REFUTED, structurally.** The pressure
   space is ~4× richer than the displacement space → nP−nU_free = 158
   guaranteed spurious Schur modes at *every* support radius; skeleton locks
   (settlement 9–27% of reference). Arithmetic, not tuning. Also: no
   particles on the drained boundary (BC by extrapolation; at s=1.6 it
   cannot even evaluate p at a domain corner).
2. **Node cloud: an expensive non-cure for stability.** Trims spurious modes
   (8→3) but does NOT kill the undrained footing checkerboard (CB 0.068 ≈
   unstabilized FEM). Strictly dominated by `-stab auto` for monolithic use.
3. **Centroid cloud ("smoothed P0"): clean inf-sup (1 mode), no checkerboard —
   a rediscovered stabilization — but worst accuracy of the viable options
   (4–9% Terzaghi; no boundary particles).**
4. **E4 — the genuine payoff (fluid life-cycle decoupled from solid
   elements)**: 10×10 column, 80%-width crack removed at Tv=0.05.
   Monolithic removal (upstream `remove element` semantics) kills the flow
   path: p below the crack never reaches 90% dissipation (p/q = 0.35 at
   Tv=1.5 vs 0.03 undamaged). A fluid measure that persists over the crack
   drains fully (Tv90 ≈ 0.99); **MLS node-cloud matches the fluid-persistent
   FEM reference to 0.6%**. Note the reference shows the payoff is the
   *architecture* (separately-discretized persistent fluid), which a
   pressure-only FE overlay would also deliver — meshless is one realization,
   whose specific edge would appear only under large deformation/fragmentation.
5. **E5 — performance: cost, not benefit.** vs FEM-p on the coupled system:
   nnz(H)/row 8.7→41.6 (~4.8×), nnz(A) ×2.5, LU fill ×2.2, factorization
   ×4.9, per-step back-solve ×2.6 (40×40 mesh; ratios grow with size; 3D
   would be steeper). Blocks are geometry-fixed and cacheable (same as the
   element's `buildStaticBlocks`), so assembly is one-time; the recurring
   price is the solver. In a staggered architecture the penalty applies only
   to the (smaller) pressure subproblem.

## E6 — the staggered seam, tried (2026-07-13, same session)

FEM solid + **persistent FEM pressure overlay**, coupled per step only through
the effective-stress seam (solid→fluid: Qᵀ Δu; fluid→solid: Q p), fixed-stress
operator split with the **oedometric modulus L = α²/(K_dr+4G/3)·M_p**
(measured near-optimal: 3.2 iters vs 11.0 for classic α²/K_dr; 0.5× oedometric
diverges — do not undershoot).

- **Accuracy**: iterated fixed-stress ≡ monolithic backward Euler (its fixed
  point): Terzaghi 2.8e-8 rel, near-undrained 4.5e-7 rel. Single-pass (fs1) is
  stable with O(Δt) splitting error: 0.09% (compressible) / 0.4% (undrained)
  vs monolithic at these Δt — usable, but iterate when it matters.
- **The naive drained split DIVERGES in 4 steps** (10 orders of magnitude) in
  BOTH regimes — soil coupling strength τ = (α²/K_dr)/storage ≈ 10³. The
  fixed-stress stabilization is mandatory, not optional.
- **Capability through the seam (E4 rerun)**: staggered fluid-persistent
  matches monolithic fluid+ to **0.00%** (probe curve), same Tv90 = 0.987,
  while element-tied fluid stays trapped (Tv90 = ∞). The E4 payoff needs NO
  monolithic coupling — the seam carries it. Iters ~9 (footing-like), ~3
  (constrained column).
- **Cost (40×40, sparse LU)**: factorization 3.3× CHEAPER than the coupled
  unsymmetric solve (18 vs 61 ms) and both sub-solves are SPD (ProfileSPD /
  MUMPS SYM=2 become legal again — recovers what the honest-p contract's
  unsymmetric tangent costs); per-step cost at 3.3 iters ≈ 1.9× the monolithic
  back-solve (5.4 vs 2.9 ms; splu does not exploit SPD — Cholesky closes most
  of that gap). Iteration count is the whole per-step story and is
  problem-dependent (3–9 measured); it degrades toward the
  incompressible-impermeable limit (known FS property).
- One implementation trap worth recording: a **fluid-first sweep with a stale
  displacement predictor self-"converges" to p≡0 on the first iterate**
  (Δu = 0 → p1 = pk → convergence test passes vacuously). Solid-first
  ordering + final momentum resolve fixes it. Any LadrunoUP staggered
  integrator must gate against exactly this failure (it converges cleanly
  onto the wrong equation — no residual warning).

## Bottom line

- Monolithic meshless-p inside LadrunoUP: **no** — no stability or accuracy
  win over shipped `-stab`/Taylor–Hood, 2.5–5× solver cost, BC and
  connectivity friction (OpenSees has no DOF carrier off-node, static DOF
  graph).
- The idea survives only as the **staggered companion-solver architecture**
  for the liquefaction-with-discontinuities target (ADR-71 §6 / P7 runway):
  pressure field with its own life-cycle, talking to the skeleton through
  the effective-stress seam. E4 is the capability plot that would justify
  that ADR; the honest alternative to name in it is a persistent
  pressure-only FE overlay (same payoff, cheaper) and coupled u-p MPM (the
  field's default for post-failure flow).
- **E6 then tried that architecture and it works end-to-end**: fixed-stress
  staggered + persistent FEM overlay reproduces the monolithic answer
  (0.00% on the crack problem), restores symmetric solvers, factors 3.3×
  cheaper, per-step ~2× at 3 iterations. The eventual ADR's recipe is
  pinned: solid-first sweep, oedometric L, persistent-fluid overlay,
  drained split forbidden.

Plots: `infsup_spectra.png`, `terzaghi_profiles.png`, `footing_maps.png`,
`removal_flowpath.png`, `perf_scaling.png`, `staggered_seam.png`.

## E7 — dynamic pins for ADR-73 P0 (2026-07-14, `staggered_pins_e7.py`)

ADR-73 P0 extension: the toy imported **as a frozen library**, plus the one
thing it lacked — **solid inertia** (row-sum lumped M, central-difference
solid loop). Protocol: every experiment printed its quantitative prediction
before measuring; refutations recorded, not smoothed. Machinery anchor: the
E6 fs1(resolve, oed, 400 steps) number reproduces to 8.88e-4 (E6 recorded
9e-4). Numbers: `e7_summary.json`; log: `e7_run_full.log`.

- **E7.1 — v1 sequencing [plan prediction REFUTED → pin changed].** The plan
  predicted fs1-without-resolve = O(Δt), ≤2× of fs1-with-resolve. Measured
  (relL2 vs same-Δt monolithic, Δt-halving ×5, both regimes): **both plan
  flavors sit on O(1) floors** — plain 0.61 at every Δt (both L), resolve
  +classic 0.41, resolve+oed 5.7e-4 (late-window floor; a special-case match
  where L ≈ true Schur compliance, column only). Root cause: single-pass
  fixed-POINT L double-counts compliance (effective storage S*+L+C instead
  of S*+C; only iteration cancels L). **Repair measured and pinned**: the
  fixed-stress **RATE form** `(S*+L)Δp + ΔtH p₁ = −QᵀΔu + L·Δp_prev` —
  same plain-analyze flow, no driver, no second solid solve, one extra
  stored vector. Late-window orders +1.00…+1.04 (both L, both regimes;
  L=oed near-deadbeat, 5.6e-6 at the finest Δt). The step-load initial layer
  is a one-step O(1) effect (full-history norm reads O(√Δt)) — harmless,
  documented. E6's "fs1 stable with O(Δt) splitting error" survives only as
  "resolve+oed is accurate at one Δt"; it is not a rate statement.
- **E7.2 — CD + implicit fluid envelope [stability CONFIRMED, with an
  accuracy rider and a rate-form carve-out].** Boundaries via Δt-bisection
  (noise IC, 6000 steps, energy blow-up detector), stab OFF, k̄ realistic:
  fixed-point L=classic **0.987×** the drained pencil 2/ω_max(M⁻¹K) —
  the ⟨A-7⟩-hedged drained-CFL claim CONFIRMED; oed 0.934×, half-oed 0.780×,
  **L=0 exactly 1.000× the undrained pencil**. Riders: (1) the drained-CFL
  fixed-point-L lane carries the same O(1) storage artifact on diffusion
  physics (0.673 relL2 vs exact, E7.3a ref) — stability true, accuracy
  false; (2) the E7.1 rate form is **UNSTABLE under CD** at every tested Δt
  (extrapolation feeds oscillatory p back as negative damping). Dynamic
  default pinned: **L=0 at Δt ≤ 0.5× the undrained pencil** (consistent +
  stable, 1.5e-3 vs exact at N=1).
- **E7.3a — subcycle (L=0 lane) [CONFIRMED].** err ~ N^1.2 (slopes
  +1.17…+1.23), all N ∈ {1…50} stable; error doubles at N ≈ 1.8 →
  **θ = N·Δt/(h²/c_v) = 0.089** pins `-subcycle auto`
  N = max(1, floor(0.09·h²/(c_v·Δt))). Artifact refs at N=1: fixed-point
  classic-L 0.673; rate form 66.8 (no blow-up, useless).
- **E7.3b — substep [REFUTED].** Coarse-Δt (Tv=0.05/step) step-load
  Terzaghi, M ∈ {1,2,5,10} ramped-Δu BE sub-steps: 0.2452 / 0.2489 /
  0.2521 / 0.2534 — flat. The binding error is the coarse-Δt splitting/lag,
  not fluid time resolution; `-substep` buys no accuracy here. Demoted from
  the P1 accuracy claim. (First-run bonus: chaining the rate-form reference
  per SUB-step explodes at M ≥ 5 — the reference must be the coarse
  increment split /M.)
- **E7.4 — fully-explicit fluid dual CFL [class CONFIRMED; ±20 % factor
  gate REFUTED, safe direction].** Coupled boundary tracks the undrained
  PENCIL: dt_cr = **1.32×** pencil on BOTH soils (5.567e-4 / 5.531e-4 s vs
  4.206e-4 / 4.196e-4). The ⟨A-1⟩ material factor √(1+K_f/(n·M_oed))
  (21.4 / 12.8) over-predicts the penalty ~1.85× (measured effective 11.6 /
  6.97): the discrete undrained pencil sits below the material formula
  (15.36 vs 21.43 on this mesh) and lagged staggering adds ~1.32×. Diffusion
  CFL slack 7.2e3× at realistic k̄ — the coupled boundary governs, as
  predicted. L=0 throughout (no iteration needed) ✓. Advisory pinned to the
  discrete pencil (≈25 % conservative vs measured); the material formula is
  documentation-only (~2× conservative).
- **E7.5 — removal under inertia [(i) CONFIRMED under the revised metric;
  original formulation refuted benignly; (ii) CONFIRMED].** E4 crack (8
  cells) dropped mid-CD-march, fluid persists. (i) The commit-jump VANISHES
  as ~dt^2.6 — no 1/Δt impulse artifact; the original "dt-independent
  commit jump" prediction was wrong in the benign direction (redistribution
  spreads over the physical transient; one commit samples ever less of it).
  The dt-independent invariant is the fixed-window jump (20·dt₀ after
  removal): spread 1.47 % across dt-halving ×4. Removal gates in
  P1/P3 tests must use the window metric. (ii) `-onRemoval drain` kFactor
  {1,10,100}: probe drainage Tv90 = 1.004 / 0.938 / 0.933 — strictly
  monotone ✓.
- **E7.6 — `-fsL` decision pair ⟨A-2⟩ [CONFIRMED; footing gain smaller than
  predicted].** Iterated fixed-stress, tol 1e-6, maxit 500, both geometries
  × both regimes: **classic never diverges and never hits maxit** (means
  11.2 / 4.3 / 11.3 / 4.3, max ≤ 18) — the KTJ-proof prediction holds; oed
  also converges everywhere (3.3 / 2.8 / 3.3 / 2.8, max ≤ 10). Ratios:
  column **3.4×** (E6's 3.2× reproduced), footing **1.6×** (the open
  measurement — under the ~3× prediction). 0.5×oed cliff: maxit-hit
  (E6 "diverges" reproduced under relative-tol semantics). Default `-fsL`
  classic, oed opt-in, never undershoot oed.

### E7 bottom line (what P1–P3 inherit)

- **Quasi-static v1 (P1):** fluid advance = the **rate-form** fs1 in the
  plain-analyze flow. No driver, no extra solid solve. The fixed-point form
  the toy used must not ship single-pass.
- **Dynamic lane (P3):** default **L=0 at Δt ≤ 0.5× undrained pencil**;
  the drained-CFL fixed-point-L mode is stability-true but
  accuracy-degraded (O(1) diffusion-rate error) — opt-in, documented;
  the rate form is quasi-static-only. `-subcycle auto` θ = 0.09; explicit
  fluid boundary = 1.32× undrained pencil, advisory uses the pencil.
- **P2 driver:** iterated fixed-stress unchanged; classic L default is
  divergence-proof on both geometries, oed 1.6–3.4× faster opt-in.

Plots: `e71_sequencing.png`, `e72_cfl_envelope.png`, `e73_multirate.png`,
`e74_explicit_fluid.png`, `e75_removal_dynamic.png`,
`e76_lvariant_iters.png`.
