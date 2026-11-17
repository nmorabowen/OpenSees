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
