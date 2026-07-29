# hypo_bearing — ADR-79 bearing-backbone campaign (scaffold + scoping findings)

**Status: SCAFFOLD + measured scoping findings. The production campaign needs a
properly sized (graded) mesh or the real SFIM model — see "What blocks the
verdict" below. Do not read the quick-pass numbers as bearing capacities.**

## Goal

Answer the ADR-78 §1 question with the now-complete geometry ladder
(`linear ⊂ corot ⊂ hypo ⊂ hypo -kozenyCarman`, ADR 79): does genuine large
strain (+ porosity/permeability evolution) bend the PDMY footing backbone
toward a bearing limit point, or reduce the reported ~9.3×-Vesic
over-hardening?

`bearing_backbone.py` runs the three viable legs (corot / hypo / hypo+kc) on a
saturated PDMY box, displacement-controlled central footing push, and reports
q/q_Vesic checkpoints + a truncation-honest softening verdict, with CSVs per
leg.

## Scoping findings (2026-07-28, all MEASURED on this script)

1. **Displacement control is mandatory** (inherited from ADR-79 P3): a
   force-controlled footing push on stage-1 PDMY diverges from the first
   increment — near-surface GPs at ~zero confinement make the tangent nearly
   singular and the first Newton iterate unbounded. `-geom linear` grinds
   without converging on this problem class and is excluded from the ladder.
2. **Surface-surcharge regularization is mandatory beyond the P3 smoke size.**
   On the N=4 box the free LATERAL DOFs of zero-confinement surface nodes are
   near-singular under stage-1 PDMY: plain Newton limit-cycles (‖Δu‖ ≈ 25 m
   iterates, residual ~500) and even KrylovNewton stalls at ~1e-6 on the very
   first 2.5 mm push step, for every geometry method. The N=3 P3 smoke sits
   just on the convergent side of the same marginality. A 10 kPa whole-surface
   surcharge (tributary-weighted nodal loads, gravity stage) regularizes it —
   with it the hypo leg marched 19 × 2.5 mm steps where it previously died at
   step 1. The Vesic normalization then needs the `q0·Nq·sq` term (in the
   script).
3. **Box confinement dominates every uniform-1 m-hex config that fits a quick
   run.** With roller sides at 0.5–1.5 B clearance the "footing" is an
   oedometer: q grows self-reinforcingly (PDMY G ∝ √p′ under rising
   confinement) to 27 MPa at s/B = 2.4 % (N=4) and thousands×Vesic (N=3) —
   physically meaningless as bearing. A credible bearing verdict needs ≥ 4–5 B
   side clearance and ~3 B depth, i.e. a GRADED mesh (fine under the footing,
   coarsening outward), not uniform hexes — or the real SFIM mesh, which is
   NOT in this repo (`References/SFIM_model/` is a private path from the
   ADR-78 workspace).
4. Per-step increments must stay ~2.5 mm (the P3-proven size); 25 mm
   first increments hard-fail every method.

## What blocks the verdict

- The real SFIM footing model (mesh + PDMY calibration + staging) — outside
  this repo. Point the script's mesh builder at it, or
- a graded-mesh benchmark builder (apeGmsh is the natural tool) with ≥ 4–5 B
  clearance, plus the surcharge regularization above.

Runtime estimate at the needed size (≈ 1–2k UP-H8, ~120 steps × 3 legs with
the substep-fallback policy): hours, not CI — run it as a campaign, keep the
CSVs.

## Files

- `bearing_backbone.py` — the runner (3 legs, staging idiom, surcharge,
  displacement-controlled push, truncation-honest summary, CSV output).
  `quick` arg = small-box sanity pass (convergence machinery only — its q
  values are the finding-3 artifact, not results).
