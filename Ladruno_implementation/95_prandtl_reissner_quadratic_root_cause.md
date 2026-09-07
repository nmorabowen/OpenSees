# Note 95 — Why every quadratic element walled on the Prandtl–Reissner collapse, and why none of them does now

**Status:** RESULTS NOTE, 2026-09-07. Branch `wp/95-prandtl-bezier-root-cause`, PR #803.
Plan and pre-registered predictions: `95_prandtl_reissner_quadratic_root_cause_plan.md`.
Phase records: `_adr95_p0_results.md` (instrumentation), `_adr95_p1_results.md` (trajectory),
`_adr95_corner_tangent_review.md` (source review), `_adr95_h4_rank_results.md` (element rank),
`_adr95_p3_results.md` (tet transfer), `_adr95_p4_results.md` (fix + decisive leg).
Predecessors: notes 81 (#721), 82 (#725), 83 (#727), gate #722, TIMs T0–T4.

## 0. The answer

**The quadratic wall was a coding defect in the vanilla UW `DruckerPrager` material, not element
physics. Nothing in any element was changed.**

The two-surface return map (cone f1 plus tension cutoff f2 at I1 = T = √(2/3)·σ_y/ρ) never
assembled the f2 residual: the `Jact(i) == 2` arm is unreachable because `Jact` holds 0/1 flags, so
the second multiplier is identically zero, the cutoff is never enforced, and a Gauss point that
crosses I1 = T keeps its unreturned stress while the material swaps to the corner tangent. The
consistent tangent additionally divided the radial-return term by the *returned* deviatoric norm
instead of the *trial* norm. The result at the first GPs to reach the cutoff, beside the footing
edge, is an operator whose normalised acoustic determinant is 1e5–1e7 times the ambient value. That
is the abrupt 1000× σ_min collapse note 82 §7.3 measured and could not name.

The linear b-bar hex never reaches the cutoff anywhere on its path (I1_min −263 kPa at the plateau):
its element-averaged volumetric strain and coarse resolution keep the heave-zone mean stress
compressive. Quadratic elements resolve the tensile spot at the footing edge, hit the defect, and
die in Newton. Dynamic relaxation (note 83) walked through because it never formed the tangent.
Both facts are now explained by one line of code.

With the return map repaired (commit `31322a47a`, `fix(adr95-p4)`):

| element | pre-fix (allowance named) | fixed material |
|---|---|---|
| LadrunoBrick20 uri, h0 = 1.0 | FLOOR at s/B 0.01114, 0.769 of exact | __H20_FIXED__ |
| TenNodeTetrahedron, tet mesh h0 = 1.0 | FLOOR at s/B 0.00161, 0.391 | TARGET at s/B 0.02, 1.051 (still hardening); long leg **TARGET s/B 0.15, 1.1704, tail 0.11 % — CAPACITY plateau**, 40–181 corner GPs returned to I1 = T, det O(1) |
| BezierTet10 std | FLOOR at s/B 0.00168, 0.410 | **TARGET s/B 0.15, 1.1824, tail 0.11 % — CAPACITY plateau**, 120–203 corner GPs returned, det O(1) |
| BezierTet10 -bbar | TARGET s/B 0.02, 0.713 at matched s/B 0.008 (no corner GP in range) | **TARGET s/B 0.15, 1.0403 of exact, tail 0.02 % — a CAPACITY plateau** |
| LadrunoBrick -bbar (control) | TARGET, 1.0850 | TARGET, 1.0850 (q_max identical to printed digits) |

Coarse-mesh over-strength (1.04–1.09 at h0 = 1.0) is the gate's known from-above convergence
(#722: 1.0842 / 0.9938 / 0.9513 at h0 = 1.0 / 0.5 / 0.25); it is not part of this note's claim.

## 1. What "the Bezier elements fail" turned out to be

Three stacked deficits (plan §0). **A**, the quadratic-class wall, is the defect above and is gone.
**B**, the tet penalty, measured for the first time at matched settlement (P3, s/B = 0.008, pre-fix
material, all legs still on-path): h8bbar 0.7925, BezierTet10 -bbar 0.7128, h20uri 0.6790 — the
Bezier b-bar tet sits *between* the two hex legs, so B is small and not the 1.4× the walled
allowances suggested. **C**, controller allowance, is unchanged and still governs how any walled
number must be quoted.

Bezier as a *basis* is exonerated twice: the pure-Lagrange H20 walls identically (note 81), and the
Bernstein b-bar tet is the first quadratic element in either campaign to reach a genuine plateau.

## 2. How the event was identified (P0–P1)

P0 added a read-only material response `ladrunoBranch` (branch, γ0, γ1, f1/f2 trial, forced-accept,
I1, min-over-directions acoustic determinant normalised by (2G)³) to the UW material, forwarded by
all four elements through `material <gp>`. P1 sampled it at every GP of the note 82 trajectory
(`quad_path_diag.py --branch`, stations every 5e-6 of s/B across the known event):

- 0 → 0.011125: **zero** f2/corner GPs in the quadratic leg, σ_min/scale flat at 2.2e-4.
- Last converged state (0.0111375): **4 corner GPs**, x = ±2.64, z = −1.5 (first row beside the
  footing edge), I1 = 0.80–0.98 ≥ T = 0.816, γ1 = 0 exactly, and the pathological determinant on
  exactly those four. Nothing else changed.
- Linear control over 0 → 0.15: never a single f2/corner GP.
- Repeat on the corrected determinant map (build cf239c9d, `_p1rep`): FLOOR at s/B 0.01124, 0.7716
  (path not bit-identical — note 82 §7.1.1 marginal-decision scatter, wall within 1 %); at the wall
  the same 4 corner GPs, normalised determinant −1.5e7 … −1.0e8 on exactly those four against −0.055
  everywhere else. n = 2, both builds, same identity.

Pre-registered predictions and their fate: H1 (corner branch) **confirmed**; H2 (loss of
ellipticity) **dead** — every plastic GP is non-elliptic in *both* elements (ψ = 0, ν = 0.45) and the
linear element completes the collapse carrying ~400 of them; H3 (Newton algebra) **unneeded**. The
branch-histogram flicker (±40 %) is a sampling artefact of perfectly plastic GPs resting on the
surface (~750 GPs within |f1| < 1e-3 in the linear leg too), not a discriminator. The `count > 3`
forced-accept bailout never fires (no `Jact =` in any log, `n_forced = 0` everywhere).

## 3. The knob (P2)

SY 0.2 → 2 → 20 kPa moves T from 0.82 to 8.2 to 82 kPa. Prediction: the quadratic wall moves out
in s/B; the linear control keeps plateauing. Measured (pre-fix material, correct determinant map):
SY 0.2 → FLOOR at s/B 0.01114 (n = 2: 0.01124); SY 2 → FLOOR at **0.01351** (+21 %); SY 20 → FLOOR at
**0.03791** (+240 %). Monotone in T, as predicted, and every wall is the same event: 4 corner GPs with a
pathological determinant (−3e7 … −5e7 at SY 2, −3e7 at SY 20) at the last converged state, none
before. The linear SY 2 control: TARGET. q is not comparable across SY (it adds cohesion to the oracle); reach and mode are.

## 4. The defect, precisely (source review)

`DruckerPrager.cpp` (upstream, fmckenna 2011): residual/Jacobian assembly switches on the *value*
`Jact(i) == 2` (lines ~521/556) while `Jact` is set to 0/1 (~474–488) → row 1 never written, g(1,1)
stays the dummy 1, dγ1 ≡ 0, f2 never enforced; the KT test then accepts because gTOL = −1e-10 and
the bailout cannot fire. The f2-only branch assembles the f1 residual. The tangent's radial-return
term (~697) divides by the recomputed final ‖η‖, not ‖η_trial‖; with ρ̄ = 0 a cone return cannot
move I1 at all, so an over-cutoff state persists and is promoted to the corner — exactly the
measured signature. The review's claim that lines 688–696 were also wrong was checked by the fix
agent and withdrawn: those expressions are upstream's and correct.

Fix (P4): index-driven assembly of both rows in every active-set combination; trial-norm divide,
guarded at ‖η‖ = 0; no parameter, parsing, tag, or send/recv change. Tests: 20 passed, including a
central-difference check of the corner tangent in all six Voigt directions (1e-3 relative) and the
P0 sentinel flipped to require I1 returned to T. The f1-only branch is **not** bit-identical, because
its tangent denominator was wrong too: gate fastest leg 1.0849417 → 1.0849561 (+1.3e-5), mode
BUDGET → TARGET, 1390 s → 181 s; linear control 1.7e-6 relative, 2628 s → 487 s. Both inside the gate
band (1.0517–1.1167). **Upstream OpenSees master still carries both defects** (dead arms at 495/530,
final-norm divide at 670) — an upstream PR candidate under the campaign's authorship rules.

## 5. Element technology (H4) — measured, secondary, unchanged

Fully plastic single-element tangents (DP perfect plasticity, homogeneous state): near-zero singular
values elastic → plastic: LadrunoBrick -bbar 6 → 7; LadrunoBrick20 uri 12 → 20 (predicted 20);
LadrunoBrick20 std 6 → 10; BezierTet10 6 → 10 (predicted 10); BezierTet10 -bbar 9 → 13. Reduced
integration does lose rank when a patch goes fully plastic, but the wall arrived with a healthy,
flat σ_min and was removed by the material fix alone, so rank loss is not the cause. It remains the
candidate for any residual quadratic-vs-linear difference in plateau quality.

## 6. Rules that held, and one that was added

- No walled number was quoted as a capacity; the identification rested on **states** at matched
  settlement, not on q.
- Controls ran first: the linear leg on the same harness, the same build, the same sampler.
- The n = 2 repeat and the corrected determinant map were run before the tangent numbers were
  quoted.
- New (quirks ledger): **stage every build; never let anything touch `dist/bin` while legs run.**
  A build/test event at 01:38 crashed three unrelated running legs (KERNEL32 trace). A Windows
  reboot at 11:00 killed five more. Every leg here was relaunched via WMI with `ADR95_DIST`
  pointing at a staged copy.

## 7. What is not claimed

- The plateau values at h0 = 1.0 are coarse-mesh numbers; the gate's refinement sequence, not this
  note, says what the exact answer is.
- Deficit B was measured on one mesh pair at one settlement.
- The rate-limited P4 agent left `_adr95_p4_results.md` with placeholders; they are filled here
  and there from the relaunched long legs.
- TIMs' Bezier b-bar legs walled at s/B 0.004–0.006 on a finer graded mesh; their build predates
  the fix. The prediction is that on the fixed material those legs plateau; it is theirs to run.

## 8. Reproducing

```bash
cd Ladruno_files/testbed/hypo_bearing
# staged build: copy dist/bin to a scratch dir, overwrite opensees.pyd from build/build/Release/OpenSeesPy.dll
set ADR95_DIST=<staged bin>
py -3.12 quad_path_diag.py --elem h8bbar --h0 1.0 --ladder linear --cond --cond-at 5e-3 --branch --suffix _p1
py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --cond --cond-at 1.5e-4 --branch --forensics --suffix _p1
py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --sy 2.0 --cond --cond-at 1e-3 --branch --suffix _p2sy2
py -3.12 tet_path_diag.py --elem tet10 --branch --cond-at 5e-4 --cond-every 25 --sfrac 0.02 --budget 200 --suffix _p3
py -3.12 tet_path_diag.py --elem beziertet10bbar --branch --cond-at 5e-3 --cond-every 50 --sfrac 0.15 --budget 200 --tmax 14400 --suffix _p4long
py -3.12 p1_flicker.py qpd_h8bbar_h1.0_p1_branch.npz
py -3.12 p1_plastic_rank.py
python -m pytest tests/test_adr95_dp_branch_response.py tests/test_adr95_dp_corner_fix.py -q
```
Every leg is quoted with its termination mode; `ladrunoBuild()` opens every log.
