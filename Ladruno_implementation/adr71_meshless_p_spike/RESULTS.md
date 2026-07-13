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
