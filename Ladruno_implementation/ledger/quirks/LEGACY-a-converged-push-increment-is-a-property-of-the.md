---
wp: LEGACY
title: "A converged push increment is a property of the MESH, not of the problem — recalibrate it after any refinement"
legacy_seq: 252
---
### A converged push increment is a property of the MESH, not of the problem — recalibrate it after any refinement
- **Bites:** reusing a step size proven on one discretization. ADR-79 P3 proved
  2.5 mm push increments on 1 m hexes; on the graded bearing mesh (0.5 m hexes
  at the surface, 4.5 B clearance, 2816 UP-H8) the same 2.5 mm increment does
  not converge **for any** geometry method — `-geom linear` included, which is
  what proves it is the model/mesh and not the geometry lane.
- **What it is NOT** (all measured, so don't re-derive them): not a tolerance
  artifact — Newton diverges (‖Δu‖ ≈ 0.2 m) and KrylovNewton stalls at
  ‖Δu‖ ≈ 3e-5 for *every* test tried, `NormDispIncr` 1e-8…1e-4,
  `RelativeNormDispIncr`, `EnergyIncr`; not a `dt` effect — 2.5 / 25 / 250 s
  all fail identically at the same displacement increment; not the u-p
  stabilization or formulation — `-stab auto 0.10/0.25/0.50`, `-stab off` and
  `-formulation bbar` all fail identically at 0.25 mm.
- **What it IS:** a FIRST-LOADING shock off the `updateMaterialStage 1` PDMY
  state. From the gravity state, 0.05 mm converges 12/12 while 0.10 mm fails —
  but once a few small increments have been taken, the SAME model happily
  accepts 0.4 mm (8x). The threshold is transient, so a fixed ladder is the
  wrong tool.
- **Rule:** drive displacement-controlled pushes with a 2-point linear ramp
  (`timeSeries Path -time 0 T -values u0 u1`) so the increment is carried by
  `dt`, then ADAPT it: halve on failure, grow back after a run of successes,
  and truncate honestly at a floor. This replaces P3's fixed `dt/10` fallback
  and self-tunes across the backbone. Also note KrylovNewton is the *primary*
  algorithm on this problem class, not a fallback — plain Newton diverged at
  every increment tested. *2026-07-28 (ADR-79 bearing campaign).*
