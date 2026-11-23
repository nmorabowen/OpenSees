# WP-232j — MPM anti-locking in rate form: method and results

ADR 232 §9 phase 4; frozen interface `tests/mpm_antilocking_interface.md` (J1–J13) with the
"232j" freeze conditions in `.claude/skills/ladruno-mpm/references/wp-briefs.md`; user guide
[[232_ladruno_mpm_guide]] ("Anti-locking"). Gates were pre-registered at kickoff, before any
phase-4 MPM number; a band is never widened, a miss is reported with its diagnosis.

## Method

- **Scheme (J1).** ΔF̄ = (det ΔF_c / det ΔF)^{1/d} ΔF on the step increment; ΔF_c is the Q1 field of
  the particle's background cell at its centre; d = ndm (plane strain: in-plane, F_zz = 1). The
  carried F is ΔF̄ F_n; the quadrature weight is the carried volume det(F) V₀ (J9).
- **Modes.** Mode F: in `Engine::advance` after the second projection, before the material call.
  Mode A: in `LadrunoMPMPoint::computeTrial`; tangent K = today's assembly at σ̄ + h dᵀ (J10).
- **Fallback (J5, freeze F1).** A particle whose cell has an inactive corner keeps its own
  increment. Every footing run reports averaged / fallback particle-steps and the fallback share in
  the plastic zone (x ≤ 3B, y ≥ −2B).
- **Solution ladder (K6).** Per accepted step: KrylovNewton (1e-6, 40), NewtonLineSearch (1e-6, 60),
  KrylovNewton (1e-5, 80; relaxed), then halve (floor 5e-4 pseudo-time), regrowing ×2.
- **Runs.** `tests/mpm_antilocking_warrant.py` (records in `tests/data/mpm/antilocking/`), and the
  Zone-A files `tests/test_ladruno_mpm_antilocking_{kinematics,tangent,switch,results}.py`.

## Results

Pending the orchestrator's build and warrant run. Each table below is filled from the committed
JSON records, with `ladrunoBuild`.

| Row | Gate | Result | Verdict |
|---|---|---|---|
| L1 | affine → rel ≤ 1e-12; rotation deps ≤ 1e-12; F_zz = 1 | pending | |
| L2 | undeformed ≤ 1e-5; 0.2 % strain ≤ 2× off | pending | |
| L3 | r3 q_lim / 514.16 ∈ [0.97, 1.08] (admission pending, freeze F2) | pending | |
| L4 | q_max within ±10 % of FE B-bar at r1, r2, r3; s = 0.2 m reached | pending | |
| L5 | MMS G2, G4 (bspline2) p ≥ 1.8, e non-increasing | pending | |
| L6 | existing MPM suite unchanged; off bit for bit; memory off ≤ +5 % | pending | |
| L7 | ZERO ≥ 0.46, SCALE ≥ 0.38; switch oracle on ∈ [513.09, 627.11], off outside | pending | |
| L8 | `mpm info antiLocking`, restart v3 / v2, Tcl = Python, refusals | pending | |

## L4 r3 — the miss and its diagnosis (round 4)

pre-registered miss, recorded 2026-11-17 (WP-232j round 4, warrant build 4fc274d51 on esmeralda). Gate row L4 at r3 (16 cells across B): the run STOPPED at t = 57 (s = 0.094 m) on the K6 ladder's halving floor (h 3.66e-4 < 5e-4), with 502 failed attempts, 7 relaxed steps and 109 halvings -- 'a stop is a miss'. Its q_max, 115.64 kPa, is inside [110.28, 134.78]; q fell to about 100 at t = 26 and recovered to 105-108. Diagnosis (b), a model / discretisation limit, not a code defect: (1) the anti-locking tangent is consistent at plastic DruckerPrager states (mpm debug tangentCheck, FD converged at eps 1e-9: mean rel 2.4e-4 on vs 2.0e-4 off at r1 t = 25-60, 2.4e-4 vs 2.3e-4 at r2 t = 30), within L2's 2x; (2) the J5 bookkeeping is inert (0 fallback particle-steps of 3.37e6) and GIMP has no images; (3) the failures are Newton stagnation of non-associated flow, shared by the fork's own FE: LadrunoBrick B-bar on the same domain and ladder needs 887 failed attempts at r3 (8.9 per step; MPM 10.7 per step) and 438 at r2 (MPM 374, which completes), and softens after its peak at s = 0.036 like the MPM; (4) the dumped worst-residual particles of the failed attempts sit mostly below and beside the footing, away from its edge (x 1.3-2.8 m, y -2.3 to -3.8 m) at stagnated residuals of O(1-10) kN, against NormUnbalance 1e-6 / 1e-5, so halving the step does not converge them. The dump's tangent ratio of 1.0 was a diagnostic aliasing defect (quirk WP-232j-ndmaterial-tangent-getters-share-a-buffer), fixed; it never reached the element tangent. Band, ladder and floor unchanged.

> **Re-recorded 2026-11-22 (WP-244, #1112), accepted by the owner.** The paragraph above is the
> 232j record. With the footing box's lower bound at −h/2, the same case on the same apt-MKL
> configuration stops at t = 74 (s = 0.128 m, s/h 1.02) with 725 failed attempts; q_max 115.64 is
> unchanged and L4 r3 is still a miss. The current text is `L4_R3_MISS` in
> `tests/ladruno_mpm_antilocking_harness.py`; the evidence is in `244_mpm_footing_track.md`.

Probes (local build `4fc274d51`, `probe_r4_tangent.py` in the round-4 report): `mpm debug
tangentCheck 1e-6 <eps>` with and without `-antiLocking off` on the same committed states.

| case | t | eps | on: max / mean | off: max / mean |
|---|---|---|---|---|
| r1 | 25 | 1e-9 | 1.08e-3 / 2.43e-4 | 9.28e-4 / 2.00e-4 |
| r1 | 40 | 1e-9 | 1.06e-3 / 2.42e-4 | 9.45e-4 / 2.01e-4 |
| r1 | 60 | 1e-9 | 1.04e-3 / 2.41e-4 | 9.71e-4 / 2.02e-4 |
| r2 | 22 | 1e-9 | 1.04e-1 / 2.81e-4 | 9.60e-2 / 2.17e-4 |
| r2 | 30 | 1e-9 | 1.95e-2 / 2.42e-4 | 1.57e-1 / 2.31e-4 |

(The max over 1152 / 4608 particles is FD branch switching at plastic points, in both legs.)
