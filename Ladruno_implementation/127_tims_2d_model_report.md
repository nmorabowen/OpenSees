---
title: "Why the Ring Walls — the fork's reply to the TIMs 2D-model requests F18–F23"
project: Ladruno
status: living draft (fifth issue, 2026-09-29; updated as the pending results land; see the revision log at the end)
date: 2026-09-29
audience: TIMs project team (2D-model act)
answers: _tims_2d_model_requests_2026-09-25.md (F18–F23)
tags:
  - report
  - evidence
  - sanisand
  - integrator
  - tims
---

# Why the Ring Walls

**Reply to the TIMs 2D-model intake of 25 September 2026 (items F18–F23). Living document, started
28 September 2026.**

This is written for the people who run the strip-footing deck. It answers every item you asked for,
says what shipped and how to use it, and keeps apart four things: what is **shipped** (merged on
`ladruno`), what is **measured** (a number with a committed script behind it), what is **pending**
(an open PR or a run still going), and what is **your decision**. Anything pending is marked as
pending; anything preliminary (a result from work still in progress) is labelled so where it is used.

| | |
|---|---|
| Fork | Ladruno / OpenSees, branch `ladruno` |
| Merged so far | #863 (WP-127), #864 (WP-132), #865 (WP-131 step 1), #866 (WP-133), #869 (WP-128), #870 (WP-136), #871 (WP-129), #872 (WP-134), #874 (WP-135, PDMY hang), #876 (WP-132 guide follow-up), #885 (WP-144/145 plans), #888 (WP-131 re-verify), #846 (WP-109, OpenMP on gcc), #868 (WP-130, CPPM under Newton; `e96f8d77d`), #893 (WP-151, R1; `fd87e396d`) |
| Still open | #878 (WP-138, footing A/B), #892 (WP-150, the wall's mechanism and the regularization memo), #894 (WP-152, the tension cutoff for zero-confinement points; under review) |
| Plan and findings A–D | `Ladruno_implementation/127_tims_2d_requests_plan.md` |

Every line of §2 of your intake was checked against the source before any work started
(`127_tims_2d_requests_plan.md`, preamble: all citations held on `fb1afe58b`; no commit had touched
`SRC/material/nD/UWmaterials/` since your `79e062367`). One reading was refined later, not refuted:
the "abrupt switch at 0.5 kPa" in the error norm is a continuous 1 kPa floor (§1, F18(a)).

---

## 0. The short version

1. **Your deck walls twice, for two different reasons.** `ModifiedEuler` (IntScheme 1) walls first
   because of a set of discrete integration defects: they let the back-stress α leave the bounding
   surface in one accepted substep and then commit states the model cannot reach (§2; WP-128 #869,
   WP-134 #872). On our copy of your deck it floors at **s/B 0.0292**, inside your own band, after
   committing ρ_α up to 13.09 and a spurious +6.4 % stiffening. SAS-ME removes that wall and reaches
   **s/B 0.0508**, where it stops on a **constitutive** singularity of DM04 (the plastic modulus at an
   α_in re-seat), not an integration defect. There is **no peak and no plateau**: q is 966.7 kPa and
   still rising (§4; WP-138 #878, WP-150 #892).
2. **Read the calibration caveat before any footing curve (§4).** The campaign SANISAND set has the
   strength of a very dense sand but dilates ~20–23× (triaxial) / ~8–14× (plane strain) less than
   stress–dilatancy requires, and peaks at 4–16 % strain. Each FE curve is consistent with the classical
   capacity of its own friction angle (§4.10); which one is physical depends on your sand's φ′. We ask
   you to confirm the calibration against your lab data (D7, §5.1).
3. **Your b8 worst point (element 1950, Gauss point 3) cannot be integrated by anything.** Its α is
   6–7× outside the bounding surface. It is not a hard point; it is an inadmissible state that an
   earlier bad increment produced. The right answer is to refuse it, and the new integrator does.
4. **The error floor you asked for in F18(a) is the wrong fix.** Today's norm already has a 1 kPa
   floor, and at a few kPa the substep count is set by stability, not accuracy. The flag exists
   (`-errFloor`, in SAS-ME), but it will not buy what you hoped.
5. **What shipped:** a new integrator, **SAS-ME (`IntScheme 129`)**, that reproduces the paper
   equations to 2e-7 relative and refuses what it cannot integrate (WP-129, #871); per-point
   post-mortem counters and a material-point replay command (WP-127, #863); a deterministic PARDISO
   mode and its guide paragraph (WP-132, #864, #876); PDMY03's critical-state constants as flags, plus a
   verdict on the PDMY "dilation brake" (WP-133, #866); the PDMY substep cap (WP-135, #874).
6. **Our integrator recommendation** (yours and the owner's to decide): `IntScheme 129`, TanType 0,
   TolR 1e-4, `-maxSubsteps 2000`, with the step policy of our runs (§4.11). TolR 1e-3 and TanType 1
   were measured and rejected. It is about the integrator only.
7. **A caution about results you already have.** Every `ModifiedEuler` SANISAND result carries an
   integration error of **6–15 % of the stress increment on 1e-4 strain increments**, measured on
   benign 20–100 kPa states (WP-129). On the footing curve that shows up as 5.5 % on the first step,
   at most 0.64 % from s/B 0.001 to 0.0174, and a spurious upturn near the `ModifiedEuler` wall (§4.4).
8. **Still pending** (updated 2026-09-29; the new results are in §0a): threading SANISAND (F19 step 2);
   the tension cutoff (WP-152, draft #894, under review); the mesh study at the peak (R2/R3); why the
   footing starts too soft (a stiffness ladder and an element check against Tatsuoka et al. 1986 are
   running). CPPM under Newton (#868) and R1 (#893) are **merged**.
9. **Yours to decide** (§5): the p′-floor rule, `D_factor`, the mesh set, what "limit load" means for a
   dense dilatant sand, the calibration (D7), the Lode parameter c ≥ 7/9 (D8) and whether to use R1
   (D9); and the data we need from you (§5.1).


---

## 0a. What changed on 29 September (fifth issue)

Sources for every number here: Esmeralda runs of our copy of your plane-strain strip deck, SAS-ME
(`IntScheme 129`, TolR 1e-4). The R1 legs ran on build `bd93c558d`, in
`~/ladruno_r1/deck{,_toyoura}/runs/<leg>/steps.csv`. The tension-cutoff legs ran on build `8ebde5cbd`, in
`~/ladruno_wp152/deck{,_toyoura}/runs/W_*/steps.csv`. All read on 2026-09-29.

1. **R1 is merged** (#893, `fd87e396d`; opt-in `-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5`).
   - On the campaign set (B/8) it carries the footing past the old wall to s/B 0.093, at q 1459 kPa and
     still rising.
   - Each part alone walls earlier: the floor only at 0.045, the hysteresis only at 0.046, both without the
     cap at 0.0525.
   - **The cap strength κ acts as a guard.** κ 0.25 / 0.5 / 0.75 give identical curves to s/B 0.045
     (902 / 902 / 901 kPa). After that they differ by at most 5 %, and not in κ order.
   - CPPM under a global Newton is merged too (#868, `e96f8d77d`).
2. **The c ≥ 0.78 route is withdrawn.** At footing scale, c = 0.80 without R1 still walls, at s/B 0.048, on
   the compression side (leg `C080_EB_off`). Raising c only delays the wall; R1 is the route out (D8
   revised).
3. **The campaign SANISAND set is a cyclic fit applied to a monotonic problem** (see D7).
   - Its constants are assembled from Gorini's cyclic Messina set, DM04 Toyoura's c and ch, and φ 33° with
     Jaky. No footing curve from it is physical.
   - For drained monotonic use only, we built a **physically bounded monotonic set (PB2)**: nb 1.652,
     A0 0.692, nd 3.5, h0 3.5, D_r 0.47 assumed. It liquefies at N ≈ 1 in cyclic tests, so do not use it for
     cyclic work.
   - With the tension cutoff (item 5), PB2 on B/8 reaches 1866 kPa at s/B 0.15 without a clear peak. That is
     above the classical band of 381–1595 kPa for φ′ 33–42°.
4. **After R1, the limiter is the free surface.**
   - Every realistically dilating sand stopped at s/B 0.01–0.03. The failing Gauss points sit at p′ → 0,
     about 0.27B outside the footing edge, and refuse with `errorAtDTmin`, `tensionAtDTmin(lowP)` or
     `maxSubsteps`. There are zero `loadingNonPosH` refusals.
   - A residual-pressure bracket shows that p_r changes how far a leg gets, not the curve. On the Kimura
     case at s/B 0.045, p_r 2 / 5 / 10 kPa give 820 / 839 / 861 kPa.
5. **A low-confinement separation, historically called the "tension cutoff" (WP-152, draft #894): pending,
   not for reported numbers yet.** *Re-labelled on the night of 2026-09-29.* Every number in items 3, 5, 6 and 7
   below that uses the cutoff comes from the **pre-review build `8ebde5cbd`**. The review fixes (C++ `7f1562c81`,
   PR at `8be9ca280`) gate the E2 trigger on a non-compressing increment, bound E1, refuse the cutoff with
   `-Presidual`, and make re-contact continuous. Re-check legs are running; the numbers will be replaced when
   they land.
   - **What it actually is.** At p′ → 0 the free-surface wall reaches SAS-ME as an accuracy/cost failure (codes 4
     and 9) before any tension. With p_sep = 0, every leg stops at s/B 0.0044 with zero tension refusals. The
     operative rule: a point whose update fails at committed p < p_sep, under a non-compressing increment, is
     treated as separated.
   - The opt-in `-sasTensionCutoff p_sep p_contact` lets a point at zero confinement separate: no tension, no
     shear, weight kept. It re-contacts when the volumetric gap closes. With it, the legs reach their targets.
   - Checks so far:
     - campaign R1 + cutoff is identical to R1 alone up to s/B 0.093;
     - halving p_sep moves the Toyoura peak by 0.2 % (2509 vs 2513 kPa);
     - it agrees with the p_r → 0 extrapolation within 2.3 %.
   - An independent review returned **merge with changes**. The low-pressure trigger is broader than
     intended, and the fixes and extra checks are in progress.
6. **The reference sand against a real footing test** (Gate 1).
   - **Setup.** DM04's own Toyoura calibration (Table 1), against Kimura et al. (1985) Fig. 9: centrifuge at
     30g, B 0.9 m prototype, D_r 85.6 %, load perpendicular to the bedding. Fig. 9 was digitized in full.
   - **The model** (R1 + cutoff, B/8):
     - peak 2015 kPa against 1953 kPa in the test, but at **s/B 0.167 against 0.091**;
     - initial secant stiffness 56 % low; at s/B 0.092, q is 1546 against 1946 kPa.
   - **An adversarial review confirmed units, digitization and parameters.** Its reading:
     - The **peak load is consistent within about ±15 %**. The contributions are mesh orientation (9 %),
       e_max/e_min (±5–8 %) and test scatter (±8 %).
     - The test footing's roughness for that series is unverified. A smooth footing would make the model
       25–35 % high.
     - **The peak is late, and the start is too soft and concave-up**, where the test is concave-down.
   - Kimura's H-bedded test at the same density peaks at s/B 0.163 and 1734 kPa. Bedding alone moves the test
     by about as much as our misfit, and DM04 carries no inherent fabric.
7. **Mesh (R2): the gap opens before the peak.**
   - On Toyoura with the cutoff (B 1.2 m), B/16 is −9 % against B/8 at s/B 0.05, and B/8 sheared 15° is +9 %
     at the peak.
   - This is WP-150 memo §4 case C. If it holds at the peak, it points to Perzyna viscoplasticity inside the
     model (R3b). The B/16 legs are still running.
8. **Running now:**
   - the Kimura leg with e and γ made consistent (e 0.635, γ 15.9 kN/m³);
   - a Kimura-case mesh bracket (B/16, and B/8 sheared);
   - a **stiffness ladder** (footing seating 16 / 40 kPa, K0 0.4, surcharge 5 kPa, G0 × 2, h0 × 2), to decide
     whether the soft start comes from the initial state or from the constitutive law;
   - an element check of DM04 against Tatsuoka et al.'s (1986) low-pressure plane-strain Toyoura tests.

9. **The element check: DM04 against Tatsuoka et al. (1986)** (drained plane strain at 4.9–392 kPa, Toyoura,
   σ1 across the bedding).
   - **Method.** The exact DM04 integrator at each test's e and σ3′. Oracle and driver from the WP-150 testbed;
     Tatsuoka's strains are external; e from V&I's e_max/e_min.
   - **DM04 is not too soft before the peak.**
     - φ′_peak within −1 to +2°.
     - Strain to peak ×0.46–0.73 at 5–10 kPa (DM04 peaks too EARLY) and ×1.04 at 49 kPa (matches).
     - E50 ×0.99 at 49 kPa.
   - **After the peak it is wrong:** too little softening (stress ratio ×1.35 at the peak + 2 %) and 2.3–2.7× too
     much dilation at 8 % strain.
   - **Doubling G0 or h0 moves the element away from the data.**
   - **Reading for the footing.**
     - The soft initial footing response is most likely the missing **small-strain stiffness**: G0 = 125 is a
       working modulus about 3× below Toyoura's G_max, and DM04 has no modulus-decay law.
     - The late footing peak most likely comes from DM04's weak post-peak softening and excess dilation (no
       progressive failure). Neither is a pre-peak defect.
     - A footing ladder (G0 ×2 / ×3 with and without b0 held fixed, seating, K0, surcharge) is running to test
       this. Its first read: the initial state moves the start by 10 % or less; the elastic modulus carries about
       80 % of the stiffening.
   - **If this is confirmed, the remedy is a model option** (a G_max decay law and/or post-peak
     dilatancy/softening). **That is your decision, not a solver setting.**

10. **The mesh question (R2/R3): answered, case C.** This is an acoustic-tensor census of the Toyoura legs
    (build `8ebde5cbd`, R1 + pre-review cutoff; method validated against WP-150's own number).
    - **Timing.** The B/16−B/8 gap passes 2 % at s/B 0.0335. That is exactly where points near the footing start to
      lose ellipticity (1–5 % of them by s/B 0.033–0.047), about 5× before the peak (0.167–0.180). **The peak is
      therefore mesh-dependent.**
    - **Orientation.** Before the peak, the band runs down the element column at the footing edge and follows the
      mesh when the mesh is sheared. The material's own band direction is 45–65° from vertical. Mesh orientation
      alone moves the peak by 8.8 %.
    - **Cause.** The vertical 'punching' bands in our fields are partly imposed by the mesh. The cause is the
      non-associated flow: with associated flow, the non-elliptic count at every peak drops to zero.
    - **What we plan.** A regularizer inside the model: Perzyna viscoplasticity (R3b), opt-in. A nonlocal void ratio
      would miss the points that lose ellipticity while still hardening. The census says the viscosity needed is
      small (H_v/2G ≈ 0.015), so its rate bias and cost should be modest. A plan-only work package is being
      drafted.
    - **The order.** The regularizer is tuned after the constitutive questions of item 9 are settled, because
      changing the model changes when and how bands form.
    - **Your decision (D-b in §5).** The tolerance is ½ the test scatter. Kimura's own N_γ scatter gives about ±4 %,
      and the gap already exceeds that by s/B 0.05.

**What this means for you.**
- The **peak load q_u** is the number we can currently stand behind, within about ±15 % on the reference
  sand.
- The **settlement at peak** and the **pre-peak stiffness** are not yet reliable.
- Your inputs in §5.1 remain the gate to a physical curve for your sand.

---

## 1. Item by item

Status key: **SHIPPED** (merged on `ladruno`), **PENDING** (open PR or run), **THEIR DECISION**,
**NOT DONE** (with the reason).

### F18(a) — an error norm with an absolute floor

**Asked.** `err = ‖dσ₂ − dσ₁‖ / max(2‖σ‖, σ_ref)` as a flag, and substeps and error against a
reference for σ_ref ∈ {0, 0.1, 1, 5} kPa, and where the admitted error is below the Newton tolerance.

**Answer.** The norm you proposed is what `ModifiedEuler` already computes, with σ_ref = 1 kPa. The
"switch at 0.5 kPa" is continuous: `‖dσ₂−dσ₁‖` below ‖σ‖ = 0.5 and `/(2‖σ‖)` above it is exactly
`/max(2‖σ‖, 1)`. A Python port with σ_ref = 1 reproduced today's C++ bit for bit on 640/640 ring
increments. Because ‖σ‖ ≥ √3·p, a floor of 0.1 or 1 kPa only acts below p ≈ 0.29 kPa, and 5 kPa
only below p ≈ 1.4 kPa. To touch your ring at p' ≈ 3.5 kPa, σ_ref would have to exceed about
12 kPa.

More important, **the cost is stability-limited, not accuracy-limited**. At constant p = 2 kPa and
η/M^b = 1.00 the substep count is 60 / 54 / 52 at TolE 1e-6 / 1e-5 / 1e-4. An accuracy-limited Heun
scheme would change 10× over that range; this changes 1.15×. A 20 kPa floor saves 2 of 52 substeps
there and doubles the error. The floor is inert where the error is small and cannot reach the tail
where it is large.

Where the admitted error sits below your Newton tolerance: on every constant-p and active-path state
at p ≥ 2 kPa, today and at σ_ref = 20 (errors 1e-4 to 2e-3 kPa per 1e-5 increment). It fails on the
ring tail (p95 2.9e-2 kPa, max 0.84 kPa), and those tail numbers are understated about 2–3× because
the reference used there was itself α-blind. The conversion from your `NormUnbalance` tolerance to a
stress error at one Gauss point (≈ 2e-3 kPa if the reference load is the footing weight, ≈ 0.12 kPa
if it is the bearing load) is an order-of-magnitude element estimate: **which vector you call "the
reference load" decides it, and only you can say.**

**Status.** SHIPPED as a flag of SAS-ME only; NOT DONE in `ModifiedEuler` (it would change nothing
measurable there).

**How to use it.** `-errFloor $sigRef` on an `IntScheme 129` material; default `P_atm/101` (1 kPa at
`P_atm = 101`), i.e. exactly `ModifiedEuler`'s implicit floor. Leave it at the default.

**Evidence.** `128_sanisand_ring_trace.md` §5 (tables §5.2, §5.3, §5.5); quirks row "The
ModifiedEuler error norm already HAS a 1 kPa floor…" (WP-128, #869); guide
`LadrunoSANISAND_implex_guide.md` §13.1.

### F18(b) — rate-form stages instead of two 6×6 tangents per substep

**Asked.** Compute `dσ = C:dε − Λ·C:m` inside the substep loop, form the tangent once at the end,
match today's results to round-off.

**Answer.** Not done in `ModifiedEuler`, for two reasons. First, it cannot be round-off neutral for
`TanType 2`: the chained "consistent" tangent consumes the per-stage 6×6s
(`127_tims_2d_requests_plan.md` finding D), and that chain has its own defect (it accumulates `T`
where the recurrence needs `dT`; quirks row, WP-129). Second, once the defects in §2 were found,
restructuring `ModifiedEuler` for speed would have made a wrong answer faster. SAS-ME is the
replacement: its stages compute increments at their own state and it forms **one continuum tangent
at the end state** for `TanType 1` and `2`.

**Status.** NOT DONE in `ModifiedEuler`; superseded by SAS-ME (SHIPPED). The allocation-free kernel
that would make the per-substep cost small is on the roadmap (§6).

**Evidence.** Plan finding D; quirks row "`ModifiedEuler`'s `TanType 2` 'consistent' tangent chain
accumulates `T` where the recurrence needs `dT`"; `LEDGER_implementations.md` WP-129 row.

### F18(c) — make `IntScheme 2` (CPPM) usable under a global Newton

**Asked.** Refuse at once instead of 2⁹ recursive halvings, a line search or better start, remove the
static work arrays, rerun F12's bearing deck.

**Answer, and a finding you did not ask for.** **The CPPM's `TanType 2` tangent had the wrong sign in
vanilla `ManzariDafalias`.** `NewtonSol` ends `Cep = -1.0 * CSigma`; the algorithmic tangent is
`+CSigma`. The local return is correct; only the matrix handed to the element is negated, so the
global Newton **diverges from its first iteration** and only the relaxed Krylov rung ever commits a
step. Checked against a finite difference of the return map: the vanilla sign is off by 2.0 relative;
the flipped sign by 1.24e-3 (one local iterate of staleness). That, more than the 2⁹ ladder, is why
F12 found scheme 2 475× shallower. F12 had read the code and called it a genuine algorithmic
tangent; nobody had compared it with a finite difference.

**A sign fix is not a consistent tangent.** Three error sources remain: it is one local iterate
stale (up to 0.27–0.53 relative at the default TolR 1e-7, because the local norm mixes strain and
stress units), the void-ratio dependence is missing from dR/dε (1e-4 to 1e-3), and after a halving
the second half-increment's tangent is handed out. With the fix, the global Newton converges in
about 3 iterations per step but is **superlinear, not quadratic** (median observed order 1.14–1.24).

Measured on F12's bearing deck (same leg, 1200 s budget, driver unchanged), with the recommended
recipe below, run back to back with an `IntScheme 1` control on the same loaded machine:

| arm | s/B at 300 / 600 / 900 / 1200 s | global iterations per committed step | load–settlement vs IntScheme 1 |
|---|---|---|---|
| recommended CPPM recipe | 0.00293 / 0.00421 / 0.00523 / 0.00626 | 3.2 (max 6) | 1.43 / 0.66 / 0.20 / 1.03 % at s/B 0.001 / 0.002 / 0.004 / 0.006 |
| `IntScheme 1`, same load | 0.00138 / 0.00250 / 0.00442 / 0.00698 | 16.8 (max 37) | — |

The recipe is ahead of `IntScheme 1` for the first 900 s (2.1×, 1.7×, 1.2×) and behind at 1200 s,
as refusals cut steps and the driver's 80-subdivision budget ran out (76/80 spent). Vanilla
`IntScheme 2` on the same deck: 3 steps, s/B 0.00002.

**Status.** **PENDING** — PR #868 (WP-130), draft, not merged. The static-array item is done there as
groundwork for F19 (the live shared state was `Matrix::Invert`'s scratch, now a stack-local LU;
`NewtonIter`'s statics are dead code).

**How to use it (once #868 merges).** Recommended recipe for `IntScheme 2` under a global Newton:

```tcl
nDMaterial LadrunoSANISAND $tag <18 params> 2 2 $JacoType $TolF $TolR \
    -cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on
```

`-cppmTangent fixed` is the **default on `LadrunoSANISAND`** (owner decision); `-cppmTangent vanilla`
reproduces the old binary bit for bit. Vanilla `nDMaterial ManzariDafalias` keeps the wrong sign.
`-cppmStart explicit` is **not** in the recipe: on the review's oracle set it doubled the error on 74
of 171 increments. Use a forwarding element (your `LadrunoQuad` is one).

**Evidence.** PR #868 body (tables "F18(c) refusal timing" and "F12 bearing deck, RECOMMENDED
recipe"); `origin/wp/130-sanisand-cppm-under-newton`: guide §9 "IntScheme 2 under a global Newton",
`Ladruno_files/testbed/hypo_bearing/wp130_f18c/tables_recipe.md`.

### F18(d) — per-point fallback from `ModifiedEuler` to CPPM

**Asked.** When `ModifiedEuler` hits `-maxSubsteps`, hand that point's increment to CPPM; refuse only
if CPPM also fails; one-element test.

**Answer.** Built as `-meFallback cppm` (needs `IntScheme 1` and `-maxSubsteps > 0`). One-element
test: a leg that `-maxSubsteps 20` refuses at step 1 runs all 10 steps with the fallback, 20 of 20
capped updates returned by CPPM, stress within 1.3 % of the uncapped integration; where CPPM also
fails the update is refused and nothing is integrated explicitly.

On your footing, though, this lever is small: the cost is **not concentrated in the ring** (§4: the
ring holds about 5–8 % of the substeps), so a per-point fallback would recover under 10 % of the
material time. It also inherits `ModifiedEuler`'s defects for every point it does not rescue. It is
not ported to SAS-ME yet (§6).

**Status.** **PENDING** — PR #868.

**Evidence.** PR #868 body ("F18(d) one-element fallback"); WP-138 census (§4).

### F18(e) — can any integrator take the b8 ring point?

**Asked.** Say plainly if the attached b8 point (p' 0.352 kPa, η 12.87) is one no integrator can
take, and trace how a committed η/M^b ≈ 6 arises.

**Answer. No integrator can take it, and none should.** The dumped state is **inadmissible**: its α
lies 6.26× past the model's bounding surface measured with the Lode angle of n (WP-128), 7.3× with
α's own Lode angle (WP-134), with `b:n = −8.19`. Driven by small probes:

- loading-type probes are "taken" only in the sense that the stress rides the cone around the bad α;
  the state stays inadmissible;
- unloading-type probes either commit **f > 0 as success** (`ModifiedEuler`, at both tolerances), or
  **teleport η from 12.9 to 1.33** in one increment through the force-accept clamp to `Mc` (tight
  `ModifiedEuler`, CPPM's fallback, RK45 even on a 1e-7 increment). A jump of that size in one
  increment is not an integration;
- the exact reference integrates gp 2 and gp 3 honestly under compression (α stays far outside,
  f ≈ 0) and **stops** on gp 3 under shear, where the rate equations are singular (0/0 in the loading
  index) and says so.

How the state arises is §2. In one line: after the low-p clamp sets α = 0, a reversal sets α_in = α = 0;
on the next compression increment h is the 1e10 sentinel, the first Heun stage moves α by about
Δs/p evaluated at the floor pressure while the increment raises p about 20×, and the stress-only
error test passes it at `dT = 1`. **Smallest reproducer:** σ = 0.0101·I, α = α_in = z = 0, one
plane-strain dε_yy = +1e-4. `ModifiedEuler` returns η 10.79 and α at 5.14× the bounding surface in one
substep with rc = 0; the exact equations give η 0.531–0.534 and 0.25×. The ring carries the
signature: all four b8 rows with α outside the bounding surface have α_in ≡ 0 exactly, and no row with
α_in ≠ 0 is outside.

SAS-ME **refuses** both b8 1950 points on entry (`startAlphaOutsideBounding`) and integrates the other
78 rows. That is the behaviour we recommend: refuse with a named code so your step controller cuts the
step, never project α back (a projection would silently rewrite history and hide the upstream defect).

**Status.** Answered (WP-128, WP-134). The refusal is SHIPPED in SAS-ME (WP-129).

**Evidence.** `128_sanisand_ring_trace.md` §0, §2, §4; `134_sanisand_reference_integrator.md` §0.4–0.6,
§6.4–6.6; guide §13.3 (ring row).

### F19 — SANISAND in the threaded state-determination loop

**Asked.** First the inventory of shared mutable state, then per-instance / thread-local / locked
state, then identity and speed-up at 1/2/4/8 threads.

**Answer.** The inventory is done. **The prime suspect for the `IntScheme 1` segfault is not a data
race: it is a print.** Under openseespy, `opserr` goes through `PythonStream` into CPython
(`PySys_FormatStderr`). An OpenMP worker that prints (the `-maxSubsteps` cap warning was the site on
WP-107's crashing run) calls into CPython with no Python thread state. This explains every row of
WP-107's evidence, including "a mutex and `omp critical` still crash", which rules out every
data-race explanation. The fix is a deferred per-thread message buffer flushed in element order after
the loop (which also makes the printed warnings identical at any thread count). The inventory also
found real races for step 2: the process-wide `LadrunoImplexGlobals` counters (on the plain path too,
not only under `-implex`), the warning budgets, `Matrix::Invert`'s shared scratch on the CPPM path, and
class-static return buffers in the plane-strain wrappers.

**Status.** Step 1 (inventory) SHIPPED (#865). Step 2 (code) **PENDING**; it waits on WP-130 (for
`IntScheme 2`) and on the deferred message path. Until then SANISAND stays refused from the threaded
loop.

**Build note for Esmeralda.** PR #846 (WP-109) merged on 2026-09-28 and flipped the CMake default
`LADRUNO_OPENMP` to ON, including gcc builds: until then a bare-cmake Linux build compiled the loop out.
An Esmeralda build from `ladruno` at or after `eeb7847d4` has the threaded loop; SANISAND still runs
serially until step 2.

**Expected gain** (inference, not measured): material update is 85.6 % of your wall time (your §1.2),
so Amdahl bounds 8 threads at about 3.9×; dynamic scheduling and an allocation-free kernel are needed
to approach it.

**Evidence.** `131_sanisand_threaded_inventory.md` §0–§6; PR #846; your intake §1.2.

### F20(a) — a cumulative per-point substep counter

**Answer.** `substepStats`: 17 columns **per integration point** (none process-wide), cumulative since
`revertToStart`, **not** reset by `revertToLastCommit`, carried by `getCopy` and the wire. So a
post-mortem after a failed `analyze` reads the real history instead of the zero you saw. It also counts
what used to be invisible: substeps **force-accepted at `dT_min` after failing the error test**, how
many of those fired the clamp to `Mc`, low-p abandons (the integrator returning at T < 1, silently),
and cap hits. Reading it never changes a number.

```python
s = ops.eleResponse(ele, "material", ip, "substepStats")
substeps, forced, abandoned, caps = s[2], s[5], s[8], s[9]
```

Under SAS-ME the census is `sasStats` (per point, since `revertToStart`, including the last refusal
code). WP-130 (pending) extends `substepStats` to 28 columns for the CPPM.

**Status.** SHIPPED (#863; `sasStats` #871).

**Evidence.** Guide §6.2 (column table); `tests/test_ladruno_sanisand_replay_counters.py`.

### F20(b) — profile scopes inside the integration

**Answer.** Added inside SAS-ME: `sanisand.sasME.predictor`, `.stateDependent`, `.stageArithmetic`,
`.stages`, `.drift`, `.alphaCheck`, `.substeps`, `.update`, `.tangent`. Measured split: about 0.44 ms
per update at a ring state and 0.025 ms at a deep state; at a ring state the stages take roughly half
to 60 % (the state-dependent quantities about 15 % of the total), drift correction about 10 %, the α
check 6–10 %; at a deep state the tangent and drift are about 11 % each.

Not added inside `ModifiedEuler` or the CPPM Newton: `ModifiedEuler` is the integrator we recommend you
leave, and its scopes would not survive the byte-identity constraint cheaply.

**Status.** SHIPPED for SAS-ME (#871); NOT DONE for `ModifiedEuler`/CPPM.

**Evidence.** PR #871 body ("F20(b) profile split"); guide §13.3; `Ladruno_files/testbed/wp129_sasme/out/`.

### F20(c) — a `"tangentEP"` response

**Answer.** `tangentEP` returns the 6×6 continuum elastoplastic tangent at the **committed** state,
whatever `TanType` the deck uses. Checked against a one-sided finite difference at a plastic state
(50 kPa, three strain directions): relative difference below 1e-4.

```python
C = ops.eleResponse(ele, "material", ip, "tangentEP")   # 36 values, row-major
```

**Status.** SHIPPED (#871).

**Evidence.** `tests/test_ladruno_sanisand_sasme.py::test_tangentEP_matches_finite_difference`.

### F21 — replay a dumped material state

**Answer.** `ladrunoSANISANDReplay` puts a private copy of a `LadrunoSANISAND` prototype into a given
(σ, α, α_in, z, e) and drives one strain increment through the same `setTrialStrain` an element uses.
It returns rc, the census of that one update, the returned state (σ, α, α_in, z, e, p, q, f before and
after, the path code) and a per-substep trace (`T, dT, err`, outcome code). Replaying step k+1 from the
committed state of step k reproduces an analysis step to 1e-9.

**Finding A — your CSVs are compression-positive.** The README says "compression negative as
OpenSees stores it", but on all 80 rows `p_kPa = +tr(σ)/3` and every normal stress is ≥ 0: the columns
are the model's internal `mSigma`. So the replay has **no default convention**; you must say which:

```tcl
ladrunoSANISANDReplay $matTag -convention compressionPositive \
    -sigma s11 s22 s33 s12 s23 s31 -alpha ... -alphaIn ... -fabric ... \
    -voidRatio $e -dStrain d11 d22 d33 g12 g23 g31 <-type 3D|PlaneStrain> <-trace 10000>
```

Shear strain is engineering (γ). α, α_in and z are projected to their deviatoric parts with a warning
(b8 row 1859/2 has tr α = 2.3e-3). A Python helper reads your CSVs and runs the standard probes:
`Ladruno_scripts/sanisand_replay.py` (`replay`, `read_ring_csv`, `probes`). For your own dumps from
now on, read `substepStats` in the same dump.

**Status.** SHIPPED (#863).

**Evidence.** Guide §6.3; quirks row "The TIMs ring-point CSVs carry the INTERNAL,
compression-POSITIVE `mSigma`"; `test_replay_reproduces_an_analysis_step`.

### F22 — a deterministic mode

**Asked.** MKL conditional numerical reproducibility for PARDISO, a list of what else is
order-dependent, byte-identical curves twice on 8 threads.

**Answer.**

```tcl
system Pardiso -deterministic              ;# MKL CNR on the AUTO branch + iparm(34)
system Pardiso -cbwr COMPATIBLE            ;# an explicit branch every x86 node can run
```

The first solve prints what MKL actually has in force, e.g.
`PARDISO deterministic mode: MKL CNR branch AUTO, iparm(34)=8 thread(s), CNR ACTIVE`. Measured on a
~22k-DOF push at 8 MKL threads, 5 runs each: mode on, 1 distinct result; mode off, 5 distinct
displacement fields.

Four things to know:

- **Across nodes with different CPUs** (your §1.6 case), AUTO picks a code path per CPU. Pin a branch
  every node can run: **`-cbwr COMPATIBLE`**. The instruction-set branches (`AVX2`, `AVX512`, …) exist
  only on Intel CPUs; on an AMD machine every one of them was refused and only `AUTO` and `COMPATIBLE`
  worked. The thread count must match too.
- **The mode is process-wide and sticky**: it stays on for every later model in the same interpreter.
  MKL refuses to set it once its BLAS/LAPACK dispatch has started (an `eigen` before the `system` line
  triggers this; an earlier PARDISO solve does not). The reliable route is the `MKL_CBWR` environment
  variable set before the process starts.
- **Serial targets only.** `OpenSees.exe` and the sequential `opensees.pyd`; MUMPS and MPI reductions
  are not covered.
- **What else is order-dependent**: the threaded element loop reduces only an integer and is
  bit-identical at 1/2/4/8 threads; SANISAND is refused from it; the `-implex` counters are process-wide
  but serial. See the table in the PARDISO recipe.

**Repeatable is not reliable.** A deterministic mode makes two runs agree; it does not make either run
more correct. Every threaded run is equally correct to machine precision. When a last-bit difference
grows into a 30 % shift in where the wall sits, the model is on a knife edge (a limit point, a yield
state that can flip, a Newton that converges right at its tolerance, an adaptive cut that can go
either way), and a different tolerance, step size or mesh would move it too. Your deck has this
character for a measured reason: with `ModifiedEuler` the stress–strain map is non-smooth (the err = 0
path, §2), so there may be no equilibrium for Newton to converge to. On the fork's own
flip-determinism deck, the first push step has **no reachable equilibrium under any tangent**; a
10⁴× tighter TolR only halves the Newton residual floor (0.20–0.35 kN → 0.12–0.13 kN), and the
"converged" first-step load under `NormDispIncr` moves by about 25 % (9.66 → 7.24) (WP-136). Use
`-deterministic` for regression tests, for reproducing a failure, and for comparing nodes; do not use
it to settle a result.

**Status.** Mode SHIPPED (#864). The "repeatable is not reliable" guide paragraph SHIPPED (#876).
Measured on Windows (AMD); not yet run on Esmeralda.

**Evidence.** `75c_pardiso_solver_recipe.md` Trap 7, "The deterministic mode"; quirks rows WP-132 (CNR
process-wide and sticky; `mkl_cbwr_set` returning -8; `-cbwr AVX2` refused on AMD; `iparm` zeroed at
every symbolic phase); `tests/test_wp132_deterministic_pardiso.py`; `136_flip_test_drift.md`.

### F23(a) — PDMY03's critical-state constants

**Answer.**

```tcl
nDMaterial PressureDependMultiYield03 $tag ... <-ei $e0> <-cs1 $v> <-cs2 $v> <-cs3 $v>
```

Flags, after every positional argument, in any order; defaults 0.6 / 0.9 / 0.02 / 0.7 (the former
hard-coded values), byte-identical when omitted (Python and Tcl baselines captured before any edit).
Found along the way and fixed: the per-material reallocation every 20 materials overwrote every
existing material's constants with the newest one's. Harmless while they were hard-coded; a silent
cross-material leak once they are user-set.

Also found, **not fixed**: `pAtm` is a static member of PDMY01/02/03, so the last material created sets
the atmospheric pressure for every material of that class. Keep one `$pa` per class per process.

**Status.** SHIPPED (#866).

### F23(b) — the PDMY "dilation brake"

**Answer.** Your reading is right that the brake is keyed to void ratio and that a dense sand never
reaches it with the default constants (from e = 0.6 it must dilate 17.5 % volumetrically at 100 kPa,
9.9 % at 1 652 kPa). It is incomplete in a way that matters: **reaching it would not help**.
`isCriticalState()` is a **crossing detector**. It returns 1 only for the increment whose start and
end lie on opposite sides of the line; past the line both are on the same side again and the full
dilatancy rule resumes (measured: the volumetric rate dips at the crossing step and is back within
0.5 % ten steps later). So retuning `ei`/`cs1..3`, now possible on PDMY03 too, moves *when* one
increment loses its dilatancy; **no choice of constants yields a plateau**. That is consistent with
your candidates with a retuned line failing the saturation gate as well. The route to a plateau is a
model whose dilatancy vanishes at critical state by construction (SANISAND's D ∝ M^d(ψ) − η, PM4Sand).
Your ten-candidate and strip numbers were not re-run.

A related defect found in WP-133 and fixed in WP-135 (#874, merged): a wild Newton iterate makes PDMY's
substep count `|Δε|/1e-5` explode to about 1e9 per call, which is why a two-element model "hung" in
`analyze`. With the fix PDMY refuses such a trial in milliseconds. Under `SSPquad` (a host that
discards the refusal) the call is bounded but the step can still be accepted.

**Status.** Note SHIPPED (#866). Hang fix SHIPPED (#874).

**Evidence.** `133_pdmy_notes.md` (b); quirks row "PDMY's 'dilation brake' `isCriticalState()` fires
only on the increment that CROSSES the critical-state line"; PR #874; `tests/test_wp135_pdmy_substep_cap.py`.

---

## 2. What walls the deck under `ModifiedEuler`

> **Scope, since the footing A/B landed (§4).** This section explains the `ModifiedEuler` wall: s/B
> 0.0292 on our copy of your deck, inside your own 0.026–0.041. SAS-ME removes it and stops later, at
> s/B 0.0508, on a different, constitutive cause (§4.6).

### 2.1 The ring, briefly

Just outside the footing edge, the top row of Gauss points sits at p' of a few kPa with η on the
bounding surface. There SANISAND's plastic modulus scales with p and its elastic moduli with √p, so the
rate equations are stiff, and an explicit scheme takes substeps sized by stability, not accuracy. That
part is physics and would cost time under any explicit integrator. It is not what stops a
`ModifiedEuler` run.

### 2.2 What stops it: a chain of discrete defects in `ModifiedEuler`

Each link is a quirks row; the ranking is the one the independent reference integrator (WP-134)
established, which corrected WP-128's first ranking of F.

| role | mechanism | what it does |
|---|---|---|
| **trigger** | **G** — α_in is re-seated once per increment | inside the substeps (α − α_in):n reaches 0, so h is the 1e10 sentinel, and then goes negative, so h < 0 and the α law becomes a repelling relaxation. 37 of 38 crossing substeps have it. The paper resets α_in at the start of each new loading process, which makes h < 0 impossible (0 of 960 runs of the exact reference). |
| **enabler** | **E** — the substep error looks at stress only | a substep that throws α 5–16× outside the bounding surface passes, because both Heun stages have the same stress increment. Adding α to the error alone keeps α inside. |
| **enabler** | **F** — a loading stage with a negative denominator is taken as elastic, and the step factor has no upper cap | both stages then agree exactly, **the error is exactly 0**, and the next substep swallows the rest of the increment. No tolerance can see it. It accounts for all 25 escapes of your `ModifiedEuler` from admissible ring starts, and for 20–65 % stress errors on benign 20–100 kPa states even at TolE 1e-8. |
| adds error | **U9** — K and G are frozen at the committed state for the whole increment | 0.6 / 6 / 24 % of the stress increment at δ = 1e-5 / 1e-4 / 1e-3. Both stages share the same wrong moduli, so the error test cannot see it and a tighter TolR does not shrink it. |
| adds error | **U10** — the loading test uses n:Δσ, not the yield-function gradient | it ignores the −(n:r)dp term, so an isotropic compression that lowers η can be read as plastic. |
| commits it | **C** — at `dT_min` a substep that failed the error test is accepted anyway, with a clamp to `Mc` | the η 12.9 → 1.33 "teleport". Uncounted until WP-127. |
| commits it | `Stress_Correction` gives up silently | when neither correction direction reduces f, it returns the uncorrected state with f > 0 and rc = 0. Worst measured: f = 11.2 kPa at p = 0.58 kPa. |

The error estimate is the reason all of this stayed hidden. It measured only stress, it read zero on
the err = 0 path, and it compared two stages that shared the same frozen moduli. So every failure mode
above produced an increment the estimator called accurate. Tightening TolR, flooring the norm, or
raising `-maxSubsteps` all act on that estimator, which is why none of them moved your wall.

### 2.3 What SAS-ME does instead

SAS-ME (`IntScheme 129`) is a Sloan–Abbo–Sheng-style explicit modified Euler written against the
oracle, not a patch of `ModifiedEuler`:

- exact elastic path (closed form in √p) for the predictor and the intersection;
- every Heun stage evaluates K, G and every state-dependent quantity at its own state (U9);
- stages classified from the true yield gradient; a loading stage with H ≤ 0 is refused or cut, never
  called elastic (F, U10);
- the error covers σ, α **and** z; TolR is always honoured; the step factor is capped at 1.1 with no
  growth after a rejection (E, F);
- the paper's α_in rule inside the increment (G);
- refusal instead of force-accept, with named codes; its own drift correction fails rather than
  returning f > TolF (C);
- a bound check on α after every substep, and a refusal of inadmissible starts.

Measured against the oracle: benign 20–100 kPa states within 2e-7 relative at TolR 1e-7 (5e-5 at
TolR 1e-4), where `ModifiedEuler` is 6–15 % off on 1e-4 increments; the smallest reproducer at
ρ 0.252, η 0.531 in 251 substeps (oracle 0.252 / 0.531; `ModifiedEuler` 5.14 in 1 substep); the ring,
624 of 640 increments integrated, the 16 from b8 1950/2–3 refused, max f at exit 1e-7, no escape. It
costs more per increment: median 17 substeps on the ring against 4 for `ModifiedEuler`, and 4–6 per
1e-5 increment on smooth monotonic chains against 1–3. That is the honest cost the α-blind test was
hiding.

**How to use it.**

```tcl
nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m $h0 $ch $nb \
    $A0 $nd $z_max $cz $Rho  129 $TanType $JacoType $TolF $TolR \
    <-errFloor 1.0> <-alphaBoundTol 0.1> <-alphaEntryTol 2> <-alphaProject 0> \
    <-sasAlphaIn reseat> <-sasErrorVars full> <-maxSubsteps $n> <-Pmin ...> <-Presidual ...>
```

- `TolR` **is** the substep tolerance. Recommended 1e-4 to 1e-7; default 1e-7. Below about 1e-8, large
  low-p increments cannot meet it above `dT_min` and the update refuses. `-honorTolR` is inert (warned).
- The defaults shown are the shipped defaults. `-sasAlphaIn stale` and `-sasErrorVars stress`
  reproduce `ModifiedEuler`'s defects G and E, for attribution only. `-alphaProject 1` projects α
  instead of refusing; it rewrites history and is off by default.
- `-implex` is refused with 129.
- `TanType 1` and `2` both return the continuum tangent at the end state; `0` returns Ce.
- Refusal codes (in `sasStats` and the warning): 1 startOutsideYield, 2 startAlphaOutsideBounding,
  3 startInadmissible, 4 errorAtDTmin, 5 loadingNonPosH, 6 tensionAtDTmin, 7 driftFailed,
  8 alphaOutsideAtDTmin, 9 maxSubsteps.
- Use a **forwarding** element. `LadrunoQuad` (your element) forwards the refusal, so your step
  controller cuts the step.

The configuration we ran on our copy of your deck (WP-138, arm E_B) was
`129 0 1 1e-7 1e-4 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0`, on
`LadrunoQuad -bbar` at B/8 (TolF 1e-7, TolR 1e-4). It is also our integrator recommendation for the
campaign (§4.11).

**Behaviour change you should know about.** Since #871, `LadrunoSANISAND::commitState` refuses to
commit a trial whose last update was refused, and this includes a **`ModifiedEuler` `-maxSubsteps` cap
hit**. Under a forwarding element nothing changes (the step already failed). Under a **discarding**
element (`SSPquad`, `stdBrick`, `BbarBrick`, …) such a deck used to commit the strain without the
stress, silently; it now fails the step and the point latches until `revertToStart`. Also, database
and restart files written by an older build will not load (the wire vector grew).

---

## 3. Cautions for results you already have

1. **Every `ModifiedEuler` SANISAND result carries an integration error of this size.** Per increment:
   6–15 % of the stress increment on 1e-4 strain increments at benign 20–100 kPa states, 20–65 % on the
   err = 0 path, and U9 alone at 0.6 / 6 / 24 % for δ = 1e-5 / 1e-4 / 1e-3, none of it visible to the
   error test or shrinking with TolR. This includes your campaign curves. How much it moves a
   load–settlement curve depends on the deck; on your footing it is 5.5 % on the first step, at most
   0.64 % from s/B 0.001 to 0.0174, and then a **spurious upturn** as the `ModifiedEuler` wall
   approaches (+6.4 % over SAS-ME at s/B 0.0292, from committed states up to ρ_α 13.09; §4.2). Treat any
   `ModifiedEuler` ring-point state, any quantity read from ring points, and any `ModifiedEuler` curve
   near its wall, as unreliable.
2. **The explicit lane's failure dumps contain inadmissible states.** The b8 1950/2–3 rows are not a
   hard point of the material; they are the product of the defects in §2. Do not calibrate or test
   anything against them except a refusal.
3. **Your §1.5 tangent comparison had a cause.** Under `IntScheme 1`, `TanType 2`'s chained tangent
   accumulates `T` where it needs `dT`, and the stress–strain map is non-smooth (the err = 0 path).
   `TanType 0` was the only dependable choice under `ModifiedEuler` for that reason. Under SAS-ME,
   `TanType 1` and `2` are the continuum tangent. Under CPPM, the vanilla `TanType 2` had the wrong sign
   (fixed by default in #868, pending).
4. **Committed steps on the relaxed rung.** Your ladder's last rung (`KrylovNewton` at 10× the
   tolerance) is not a small print item on this deck. On our copy, the share of the settlement accepted
   on that rung is 88.5 % (`ModifiedEuler`), 90.9 % (SAS-ME), 92.9 % (SAS-ME at TolR 1e-3) and 81.7 %
   (SAS-ME at B/16) (§4.4). What that acceptance does to q has not been isolated by any of our arms.
   Report the rung of every committed step, and the share of the settlement committed at the relaxed
   tolerance.
5. **The 30 % run-to-run shift is a signal about the deck, not the solver** (F22). Deterministic mode
   will make the two runs agree; it will not tell you which is right.
6. **The `-Presidual` 1.01 / 5.05 kPa comparison of your §1.4 was made with the defective integrator.**
   Re-measured under SAS-ME (§4.7): no residual pressure from 0.5 to 20 kPa removes the wall,
   0.5 kPa brings it earlier, and larger values reach further only through an apparent cohesion. So it
   does not justify a floor (D1).
7. **PDMY under `SSPquad`**: a refused trial is discarded by the host (WP-135 bounds the time, not the
   acceptance). Prefer a forwarding element (`quad`, `LadrunoQuad`, the u-p family) for any deck that
   depends on a material refusal.

---

## 4. The footing A/B (WP-138)

> **CALIBRATION CAVEAT — read this before any curve in this section.**
>
> The WP-150 element tests (#892, memo §10, commit `e14703ca7`) ran the campaign SANISAND set on the
> exact WP-134 oracle, in drained triaxial and plane-strain compression at p0 = 10 / 50 / 150 / 500 kPa:
>
> - **The strength is that of a very dense sand.** Plane-strain φ′_peak falls from 60.1° to 44.9° as
>   p0 rises from 10 to 500 kPa.
> - **The dilatancy is not.** It dilates **~20–23× less in triaxial** than stress–dilatancy (Bolton 1986)
>   requires for that strength (φ′_cs 33.0° from Mc), and **~8–14× less in plane strain**. The
>   plane-strain figure uses an *estimated* plane-strain critical-state angle (≈ 39.5°), because the set
>   does not reach critical state by 25 % strain.
> - **It peaks late:** at 4–16 % axial strain, where DM04's own lab-calibrated Toyoura set at the same
>   density peaks at 1–5 %. A0 = 0.05 is 14× below Toyoura's value.
> - **The ring dilates even less.** UW's `D_factor` never fires in these tests (p′ ≥ 10 kPa throughout);
>   at the ring (p′ ≈ 3–5 kPa) it cuts the dilatancy further.
>
> So **none of the footing curves here is called physical** until you confirm the calibration against
> your lab data. That covers SANISAND (966.7 kPa at s/B 0.0508, still rising), the fork's
> `DruckerPrager` control (38°, ψ = 0; max 824.2 kPa) and your own `PressureDependMultiYield` control
> (PDMY01, 33° cone; limit point 417.6 kPa at s/B 0.116, your intake §1.1). Our request, in one
> sentence:
>
> *"The campaign SANISAND set reproduces the strength of a very dense sand but dilates ~20–23×
> (triaxial) / ~8–14× (plane strain) less than stress–dilatancy (Bolton 1986) requires, and peaks at
> 4–16 % strain (a lab-calibrated DM04 set at the same density peaks at 1–5 %). Please confirm the
> calibration against your lab data (φ′_peak, strain at peak, dilatancy) before the footing curves are
> used."*
>
> Everything below is about why the integration stops and what it costs. It says nothing about where a
> correctly calibrated footing would peak. §4.10 puts each curve against the classical capacity of its
> own friction angle.

Source for this section, unless stated: the WP-138 report `138_footing_sas_me_ab.md` (#878, draft,
at `762be8332`), with the mechanism, the localization analysis and the capacity bands from the WP-150 memo
`150_sanisand_regularization_memo.md` (#892, draft; §10 at `e14703ca7`, §11 at `00198f278`).

### 4.1 Setup

The fork's own copy of your deck, built from the §1 spec of your intake. Nothing in your Workbench was
run or edited. Plane-strain strip, B = 1.5 m, full width; `LadrunoQuad -bbar`; B/8 (2 430 elements,
9 720 Gauss points); `system Pardiso`; `NormUnbalance` 1e-5 × the applied vertical load (0.0415 kN);
Newton (25) → NewtonLineSearch (40) → KrylovNewton (60, tolerance × 10); ds from 2e-5 m, doubled after
6 good steps up to 1e-3 m, halved on a failed ladder, **FLOOR** when ds < 2e-7 m. Where your spec was
silent we filled the gap and said so (mesh grading reconstructed from your `ring_points_b8.csv`,
γ′ = 9.81 kN/m³, a rough guided footing, F10's step controller): #878 §1, gaps G1–G6.

**The arms that count** ran on Esmeralda, build `ladruno` **`7936ed6e0`**, each alone on its node with
MKL and OpenMP at 1 thread, so their wall clocks compare (#878 §8). All use the material line of your
intake with `-flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0`:

- **E_A**: `ModifiedEuler` (`IntScheme 1`), TanType 0.
- **E_B**: SAS-ME (`IntScheme 129`), TanType 0, TolR 1e-4. The reference SAS-ME arm.
- **E_D**: E_B with TolR 1e-3.
- **E_C2**: E_B with TanType 1, `-maxSubsteps 20000`, KrylovNewton at 1× the tolerance.
- **E_B16**: E_B on a **B/16** mesh (9 720 elements, 38 880 Gauss points).
- **Control**: UW `DruckerPrager` 38°, ψ = 0, on the same deck (a local run).

The earlier local legs (A and B, to s/B 0.0174) are kept for the early-curve comparison and the replay
study (§4.4, §4.5). E_B
reproduced the local B leg to 1e-5 kPa through step 40 (#878 §4).

### 4.2 Where each arm stops

| arm | integrator | s/B at FLOOR | q (kPa) | first `loadingNonPosH` at s/B | refusals (converged-step census) | push wall (h) |
|---|---|---|---|---|---|---|
| E_A | `ModifiedEuler` | **0.0292** | 701.8 | — (no refusal path) | 542 cap hits, 3 943 forced at dT_min | 2.90 |
| E_B | SAS-ME, TolR 1e-4 | **0.0508** | 966.7 | **0.0363** | NonPosH 232, maxSubsteps 209, errorAtDTmin 1 | 4.32 |
| E_D | SAS-ME, TolR 1e-3 | 0.0410 | 808.3 | 0.0334 | maxSubsteps 13 316, errorAtDTmin 190, NonPosH 173 | 4.67 |
| E_C2 | SAS-ME, TanType 1 | 0.0114 | 317.0 | 0.0064 | errorAtDTmin 175 011, maxSubsteps 40 683, NonPosH 409 | 12.35 |
| E_B16 | SAS-ME, B/16 | 0.0135 | 352.8 | 0.0135 | maxSubsteps 115, NonPosH 11 | 3.75 |
| control | `DruckerPrager` 38°, ψ = 0 | 0.15 (target reached) | 752.0 (max 824.2) | — | — | 0.07 |

(#878 §0 and §8. The census sums the per-step refusal lines of the converged steps; the final ladder at
the floor is not in it.)

**Every SANISAND arm stops on the step floor.** The SAS-ME arms reach it on `loadingNonPosH` refusals.
**E_A cannot refuse**: `ModifiedEuler` has no refusal path, so it floors through cap hits and
acceptances forced at dT_min (542 and 3 943 over the run), and it **commits** what it cannot
integrate. Over s/B 0.026–0.0293 it forced 3 786 acceptances and the committed ρ_α reached **13.09**;
SAS-ME over the same window stays at ρ_α ≤ 1.004 with no forced acceptance. The ModifiedEuler curve
bends up there, to **+6.4 %** over E_B at the same s/B, and that stiffening is spurious (#878 §8.4).
E_A's wall sits inside your own `ModifiedEuler` band (s/B 0.026–0.041), so our copy reproduces your
wall. Read your `ModifiedEuler` curves near their wall with this in mind.

### 4.3 The verdict

**SAS-ME moves the wall from s/B 0.0292 to 0.0508, but it does not remove it.** There is **no peak and
no plateau** on this deck: at the wall q = 966.7 kPa and still rising (q_max = q_end; the slope over the
last 0.005 s/B is 0.24× the initial slope). The first `loadingNonPosH` refusal comes at s/B 0.0363, and
from there the refusals accumulate until they end the run (#878 §0, §8.1, §8.3).

**The cause is constitutive, not integration.** No integrator setting lifts it (TolR 1e-3 walls earlier,
TanType 1 much earlier, §4.4), the independent oracle stops at the same kind of state, and the
refusing points are pre-peak (ρ_α < 1). §4.6 gives the mechanism. This revises the first issue of this
report, which called the wall an integration failure: that is true of the `ModifiedEuler` wall
(0.0292; §2), not of the one SAS-ME reaches.

### 4.4 Accuracy and cost

**Per-increment accuracy** (replays of real increments against the oracle, §4.5; #878 §5.3). SAS-ME
sits at **0.6–2e-4·p′** at every checkpoint: that is its TolR 1e-4 error control. `ModifiedEuler`'s
error depends on the state: it is **about 20× worse only at the onset of ring plasticity** (s/B ≈ 0.001),
and **at parity from s/B ≈ 0.01**, apart from isolated low-p outliers (2e-3·p′, a point taken in one
substep). It is not a blanket accuracy factor.

**The load–settlement curves**, `ModifiedEuler` against SAS-ME on the local legs (#878 §5.1):

| s/B range | max \|q_B − q_A\| / q_A |
|---|---|
| 0–0.001 (the worst is the first step, s/B 1.3e-5) | 5.5 % |
| 0.001–0.005 | 0.64 % |
| 0.005–0.0174 | 0.51 % |

Past s/B 0.016 the Esmeralda pair separates, to 6.1 % over 0.016–0.029, as E_A turns up (§4.2).

**Wall clock and substeps per 0.01 s/B** (#878 §8.2; hours / 1e9 substeps, partial intervals prorated):

| arm | 0–0.01 | 0.01–0.02 | 0.02–0.03 | 0.03–0.04 | 0.04–0.05 | whole run, h per 0.01 s/B |
|---|---|---|---|---|---|---|
| E_A (`ModifiedEuler`) | 0.58 / 0.48 | 1.14 / 0.94 | 1.24 / 1.01 | — | — | 0.99 |
| E_B (SAS-ME) | 0.40 / 0.49 | 0.82 / 1.02 | 0.93 / 1.09 | 1.02 / 1.14 | 1.05 / 1.14 | 0.85 |
| E_D (TolR 1e-3) | 0.47 / 0.44 | 0.99 / 0.97 | 0.98 / 0.96 | 1.43 / 1.39 | 7.82 / 7.88 (to 0.041) | 1.14 |
| E_C2 (TanType 1) | 5.12 / 5.05 | 50.0 / 56.6 (to 0.0114) | — | — | — | 10.77 |
| E_B16 (B/16) | 2.46 / 2.71 | 3.65 / 4.16 (to 0.0135) | — | — | — | 2.77 |

SAS-ME takes about the same substeps per unit settlement as `ModifiedEuler` and is **1.3–1.45× cheaper
in wall clock** per unit s/B. E_B's cost is flat at about 1 h per 0.01 s/B from s/B 0.01 to its wall, so
the 0.0508 wall is not a budget stop.

**Where the cost is.** The ring (p′ < 10 kPa) holds 5–8 % of all substeps and the 100 costliest points
7–9 % (#878 §6). The cost is spread over the whole plastic zone, so a per-point fallback for the worst
points (F18(d)) would recover under 10 %. The multiplier is the global iteration count × every point's
update.

**Most of the settlement is accepted on the relaxed rung.** The KrylovNewton rung accepts at 10× the
test tolerance (0.415 kN against 0.0415 kN). The share of the settlement accepted there is **88.5 %**
(E_A), **90.9 %** (E_B), **92.9 %** (E_D) and **81.7 %** (E_B16) (#878 §8.2). The effect of that
acceptance on q is **not isolated by any arm**; E_C2 changed the Krylov tolerance together with two
other settings.

**TolR 1e-3 is not a lever (E_D).** It walls **earlier**, at s/B 0.0410 against 0.0508. It costs the
same as E_B per unit s/B up to s/B 0.038 and then rises to 4.6× E_B's (0.038–0.041); it has 13 316 maxSubsteps
refusals against 209. Its q runs **−3.7 %** (median) below E_B over s/B 0.02–0.041 (range −4.3 % to
−1.9 %), and −1.4 % (median) over 0.001–0.02. No saving, an earlier wall: keep TolR 1e-4 (#878 §8.5).

**TanType 1 is not viable here (E_C2).** Global iterations per step fall (median 11 → 4 → 2), but the
accepted step collapses with them (median ds 1.25e-6 m past s/B 0.01, where E_B runs at 1e-3 m). E_C2
floors at **s/B 0.0114** after 12.35 h, 13× E_B's cost per unit s/B, while its curve stays within
1.45 % of E_B. The consistent tangent buys a step-size collapse, not settlement (#878 §8.6). This
supersedes the stopped local TanType 1 arm of the first issue.

### 4.5 The replay figures, reconciled

The first issue quoted two pairs of figures for the per-increment error on the footing's real
increments and asked which increments each covers. Both are right, for different increments (#878
§5.3). Each row is a Gauss point's committed state plus the strain increment it actually received in
the next converged step, replayed through `ModifiedEuler`, SAS-ME and the WP-134 oracle (Radau, rtol
1e-10); the error is ‖σ − σ_oracle‖ / p′ at the start of the increment.

| committed state (local A run) | increment | points | ME median | SAS-ME median |
|---|---|---|---|---|
| step 25, s/B 0.00141 (onset of ring plasticity) | step 26, ds 3.2e-4 m | 44 | **1.16e-3** | **5.74e-5** |
| step 50, s/B 0.0125 | step 51, ds 1.25e-4 m | 46 | 6.37e-5 | 6.28e-5 |
| step 75, s/B 0.0165 | step 76, ds 5.0e-4 m | 43 | 9.42e-5 | 1.85e-4 |
| step 78, s/B 0.0172 (last converged pair) | step 79, ds 2.5e-4 m | 44 | **6.85e-5** (max 1.98e-3) | **9.99e-5** |

- "ModifiedEuler 1.2e-3·p′ vs SAS-ME 6e-5·p′" is the **step 25 → 26** set.
- "≈ 7e-5 vs ≈ 1e-4" is the **step 78 → 79** set; there the `ModifiedEuler` median is the lower one, its
  maximum is not.
- Both are **medians over ~44 selected worst points** (highest ρ_α, lowest p′, most substeps), not over
  the 9 720 points of the mesh, and **both cover s/B ≤ 0.0174**. No replay exists at the Esmeralda walls.
- No real increment up to s/B 0.017 was refused by either integrator, and the oracle integrated all of
  them.

### 4.6 Why the wall: the mechanism

DM04's plastic modulus has a singularity at an α_in re-seat (#892 memo §1.4; #878 §10):

```
a = (α − α_in):n,   h = b0 / a,   Kp = ⅔·p·h·(b:n)
```

At a re-seat α_in := α, so a = 0 and h is infinite (capped at 1e10 in the code). If b:n > 0 the
modulus is large and positive; if b:n ≤ 0 it goes to −∞ and the increment has no solution. The
singular set is **{a = 0, b:n ≤ 0}**, and `loadingNonPosH` is the SAS-ME refusal that names it.

- **The refusing points are pre-peak** (ρ_α 0.93–0.96, dense, ψ ≈ −0.1, inside the bounding surface) and
  they chatter: E_B makes 10.1 million re-seats (#892 §1.3).
- **It is in the continuum equations, not in SAS-ME.** The WP-134 exact oracle stopped on the same 0/0 at
  your ring point 1950/3 (#872; #892 §1.4).
- **No integrator knob lifts it** (§4.4), and **no boundary-value regularizer** can lift an unbounded
  negative modulus at a point. Any cure is a change to the model (R1, §5).

**How the equations reach it (WP-151 memo §2.2 on #893).** The exact rate
equations reach this set through a **Zeno accumulation of re-seats on the b:n → 0⁺ side**: after a
re-seat, h = ∞ makes α slide along b; with b nearly perpendicular to n that slide rotates n, a turns
negative and the next re-seat follows. At E_B's refuser 1880/1 the re-seat intervals run 1.7e-2,
1.5e-3, 5.6e-5, 2.4e-6, … and accumulate at a finite time where a = 0 and b:n = 1.9e-8 > 0.
**`loadingNonPosH` is only the b:n < 0 exit of that sequence.**

**What makes the wall states singular: the non-convex extension side** (WP-151 memo §2.5; #892 §2.3).
All three wall refusers have n on the extension side (cos 3θ −1.00 / −0.36 / −0.88), while the
non-elliptic band points of §4.8 sit on the compression side (0.09 % / 0.14 % of them extension-side at
the E_B / E_B16 walls): the wall and the bands are separate phenomena. With c = 0.71 < 7/9 the Lode
interpolation is concave on the extension side (§4.9). Driven at c = 0.80, the same five committed wall
states fail **0 of 320** exact trials against **102 of 320** at c = 0.71. These are c = 0.71 states
driven at c = 0.80: a sensitivity test, not a c = 0.80 footing run. That run is under way (§6).

**Two routes out of the wall; the choice is yours (D8, D9).**
1. **R1** (WP-151, #893): an opt-in model-level fix at any c, no recalibration (§5, D9).
2. *(Withdrawn 2026-09-29: at footing scale this only delays the wall; see §0a.)* **A calibration with c ≥ 0.78**: at c = 0.80 the extension strength M_e = c·M_c rises 13 %; it also
   removes the extension ill-conditioning of §4.9.

**Related literature.** Stress overshooting of bounding-surface models at load reversals, where the
reversal memory is re-seated, is a documented source of numerical instability in boundary-value
problems: Chen, Ghorbani, Zhang & Kodikara (2022), "Stress overshooting solution for soil plasticity
models", *Comput. Geotech.* 152, 105008, which finds that the definition of the plastic modulus and
hardening law governs it; and Ghorbani, Chen, Kodikara, Carter & McCartney (2023), "Memory repositioning
in soil plasticity models used in contact problems", *Comput. Mech.* 71, 385–408. We have not verified a
published SANISAND footing that stops on this exact singular set.

### 4.7 Sensitivity ladders — final

Each leg is E_B with one knob changed (the S1 → S4 ablation is cumulative); all legs ran to their end
(#878 §11, records under `Ladruno_files/testbed/footing_sas_me_ab/ladders_final/`). The S and A0/h0
legs shared nodes, so no wall clock is quoted.

| ladder | leg | first NonPosH at s/B | wall (FLOOR) at s/B | q at the wall (kPa) |
|---|---|---|---|---|
| reference | E_B (Presidual 0, e 0.6944, A0 0.05, h0 1.3) | 0.0363 | 0.0508 | 966.7 |
| Presidual | 0.5 / 1 / 2 / 5 / 10 / 20 kPa | **0.0182** / 0.0333 / 0.0346 / 0.0416 / 0.0373 / 0.0535 | 0.0303 / 0.0421 / 0.0434 / 0.0676 / 0.0856 / 0.0964 | 683 / 870 / 901 / 1 286 / 1 602 / 1 979 |
| e_init | 0.65 / 0.75 / 0.80 / 0.85 | 0.0395 / 0.0349 / 0.0525 / 0.0414 | 0.0457 / 0.0431 / 0.0574 / 0.0442 | 1 456 / 500 / 341 / 181 |
| ablation | S1 z_max = 0 → S2 + n_b = 0 → S3 + n_d = 0 → S4 + A0 = 0.001 | 0.0355 / 0.0269 / 0.0347 / **0.0426** | 0.0374 / 0.0623 / 0.0601 / 0.0499 | 790 / 419 / 385 / 353 |
| A0 | 0.02 / 0.10 | 0.0237 / 0.0310 | 0.0315 / 0.0361 | 678 / 835 |
| h0 | × 3 (3.9) | **0.0091** | 0.0188 | 770 |

- **Every leg walls on `loadingNonPosH`: no material switch removes it.** The ablation strips fabric,
  the peak, the critical-state dilatancy surface and finally dilatancy itself; dilatancy off (S4) only
  **delays** the onset (0.0363 → 0.0426). S2–S4 also carry about 40 % of E_B's load at the same s/B (q at
  s/B 0.03: 292 / 277 / 278 kPa against 673). An interim snapshot of 16:20, taken before S4 reached its
  onset, had suggested that only killing the dilatancy clears the refusal; that reading is withdrawn
  (#892 memo §1.4 at `2a82e2046`).
- **Presidual 0.5–20 kPa never clears it**, and the onset is **non-monotonic**: 0.5 kPa brings it
  *earlier* (0.0182) than Presidual 0. A larger Presidual walls later (up to 0.0964) only by stiffening
  the response (1 979 kPa at the wall for 20 kPa): an apparent cohesion, not a cure.
- **e_init 0.65–0.85 never clears it.** A0 is non-monotonic (onset 0.0237 / 0.0363 / 0.0310 at A0
  0.02 / 0.05 / 0.10). **h0 × 3 brings the onset down to 0.0091.**
- **What no leg changed is the Lode ratio c** (§4.6, §4.9). Every leg keeps c = 0.71 < 7/9. The footing
  test of that reading (c = 0.80, R1 off) is running.

**Your §1.4 `-Presidual` comparison (caution 6), re-measured under SAS-ME.** Your conclusion holds in
the sense that matters: no residual pressure from 0.5 to 20 kPa removes the refusal. It does move where
the run stops, in both directions: 0.5 kPa floors earlier (s/B 0.0303 against 0.0508), and the larger
values reach further only by adding strength that is not in the sand. Do not use `-Presidual` to push
past the wall (D1).

### 4.8 Mesh: B/16

E_B16 floors at s/B 0.0135, earlier than B/8 (#878 §9). Up to there, B/16 runs softer than B/8 from
s/B 0.010:

| s/B | 0.002 | 0.005 | 0.008 | 0.010 | 0.012 | 0.013 | 0.0135 |
|---|---|---|---|---|---|---|---|
| q_B16 / q_B8 − 1 | −0.58 % | −0.14 % | −0.91 % | −3.89 % | −4.95 % | −4.84 % | −3.57 % |

**The band is one element wide on both meshes** (full width at half maximum of the incremental shear
strain 1.00–1.11 element sizes), it halves with the element and it follows the mesh lines.

**What it is: non-associated localization that begins while the material is still hardening** (#892
§2.1, Rudnicki & Rice 1975). A plane-strain acoustic-tensor scan of the continuum tangent finds
det ≤ 0 at **16.8 % of the Gauss points at s/B 0.011 on B/8** (9.1 % at s/B 0.0096 on B/16, 16.9 % at
its wall), where H/2G ≈ 1.05–1.07. The same states with associated flow are **elliptic everywhere**; only 17
of 9 720 points are post-peak even at s/B 0.0508. The cause is the flow rule: friction of 45–60° against
dilation of 1–2° (§4 caveat), so the bands are partly a product of the calibration.

**It is not ψ-softening**, so a nonlocal ψ̄ (void-ratio averaging) or a crack band would not treat it.
The wall (§4.6) is a separate matter: the floor refusers are pre-peak on both meshes, and the mesh
dependence of the wall itself is not settled (R2, §6).

### 4.9 A second calibration item: Lode convexity (c ≥ 7/9)

With **c = 0.71 < 7/9**, DM04's Lode interpolation g(θ) is **non-convex at the extension meridian**.

- **What we saw** (WP-151 memo §6.3 on #893, harness and results committed under
  `Ladruno_files/testbed/sanisand_reseat_r1/`). In
  undrained cyclic triaxial (CTXu) at e 0.6944, CSR 0.2, a perturbation of 1e-9 (round-off level) decides
  between **5 % double amplitude at N = 8** and **no 5 % DA by N = 20**, and it does so identically with
  and without the R1 fix: it is a bifurcation of DM04's own axisymmetric extension path. At **c = 0.80**
  the path stays axisymmetric and both give N = 16.
- **Why.** With g = 2c / [(1+c) − (1−c)·cos 3θ], at the extension meridian (θ = 60°): g = c, g′ = 0 and
  g″ = 4.5·c·(1−c). Convexity of the polar curve r(θ) needs r² + 2r′² − r·r″ ≥ 0, which at r′ = 0 is
  r ≥ r″, i.e. c ≥ 4.5·c·(1−c), i.e. **c ≥ 7/9 ≈ 0.778**.
- **Recommendation.** Keep **c ≥ 0.78**, or treat axisymmetric-extension tests, and extension zones in
  a boundary-value problem, as ill-conditioned under the present set (D8).

### 4.10 The three curves against classical bearing capacity

A check on the FE, not on the physics (#892 memo §11, commit `00198f278`). The classical rough-strip
capacity of **this** deck (γ′ 9.81 kN/m³, B 1.5 m, surcharge 7.65 kPa) is q_u = ½·γ′·B·N_γ + q·N_q,
with Martin's (2005) exact N_γ by characteristics (as reproduced by Han et al. 2016, Table 2) at 30°,
35°, 40° and 45°, log-interpolated only between those angles, and the exact N_q:

| φ′ | 30° | 33° | 35° | 38° | 40° | 42° | 45° |
|---|---|---|---|---|---|---|---|
| q_u (kPa) | 250 | 381 | 509 | 812 | 1 121 | 1 595 | 2 753 |
| source | Martin | interp. | Martin | interp. | Martin | interp. | Martin |

Each FE control sits on its own cone:

| control | FE result | classical q_u for its own φ′ | reading |
|---|---|---|---|
| your PDMY01, 33° | 417.6 kPa at s/B 0.116 | 381 kPa | +10 % |
| our `DruckerPrager`, 38°, ψ = 0 | a plateau of ~700–820 kPa over s/B 0.06–0.15 | 812 kPa | at or below the associated value, as expected for ψ < φ |
| SANISAND, campaign set | 967 kPa at s/B 0.05, still rising | its own element φ′_ps,peak at the footing's p′ ≈ 50–150 kPa is 51–55° (T5), so q_u ≥ 2 753 kPa | about ⅓ of its classical capacity mobilized, which is what the late element peak (4–16 %) predicts |

**The three FE curves are each consistent with their own constitutive strength** (PDMY01 33° ≈ 381 kPa
classical, DP 38° ≈ 812 kPa, SANISAND's own 51–55° ⇒ ≥ 2.75 MPa). **The physical capacity follows from
your sand's φ′**, and those inputs are owed (§5.1). With your lab φ′, the physical band is
[q_u(φ′_cs,ps), q_u(φ′_ps,peak at the footing's mean p′)]: the operative angle lies between the
critical-state and the peak angle because of progressive failure and stress level (Lau & Bolton 2011;
Perkins & Madson 2000; Loukidis & Salgado 2011). "Which curve is physical" is the same question as
"what is your sand's operative φ′".

### 4.11 The integrator recommendation (the owner and you decide)

**`IntScheme 129` (SAS-ME), TanType 0, TolR 1e-4, `-maxSubsteps 2000`, on the step policy these runs
used** (#878 §12):

```tcl
nDMaterial LadrunoSANISAND $tag 264.32 0.312885 0.6944 1.3309 0.71 0.027 0.83 0.45 101 0.005 1.3 0.968 3.5 0.05 5.75 12.5 1100 2.0 \
    129 0 1 1e-7 1e-4 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0
```

Step policy: ds0 = 2e-5 m, × 2 after 6 good steps up to 1e-3 m, ÷ 2 on a failed ladder, floor at
ds < 2e-7 m; Newton (25) → NewtonLineSearch (40) → KrylovNewton (60, tolerance × 10); `NormUnbalance`
1e-5 × the applied vertical load; 1 MKL thread with `MKL_CBWR=COMPATIBLE`, or `system Pardiso
-deterministic` on current builds.

Why: it reaches 1.74× the settlement of `ModifiedEuler` at 1.3–1.45× lower cost per unit s/B; its
per-increment error is controlled; and it refuses a state it cannot integrate where `ModifiedEuler`
commits ρ_α 13 and a spurious +6.4 %. Rejected: TolR 1e-3 (earlier wall, no saving) and TanType 1
(step collapse, 13× cost).

**This is a recommendation about the integrator only.** It says nothing about the calibration (the
caveat at the top of this section), and it does not buy a capacity: the wall stays at s/B 0.0508 on this
deck, and about 90 % of the settlement is accepted on the relaxed KrylovNewton rung.

---

## 5. Your decisions

These are recommendations. Nobody outside the calibration and the project can make them.

| # | decision | our recommendation | why |
|---|---|---|---|
| D1 | **The p′-floor rule** | (a) `-Pmin` ≤ 0.5 kPa; (b) `-Presidual 0` (at most 1 kPa if used, and then declared as a regularisation, not physics); (c) report the limit load at floor F and at F/2, and accept the floor if the load moves by less than about 2 %; (d) report the number of Gauss points at the floor at the limit state. **Do not use `-Presidual` to get past the wall** | **re-measured under SAS-ME (§4.7): Presidual 0.5–20 kPa never removes the `loadingNonPosH` wall**; its onset is non-monotonic (0.5 kPa walls *earlier*), and larger values reach further only through an apparent cohesion (q at s/B 0.03: 675 → 835 kPa from 0.5 to 20 kPa). Precedent for a small floor: PM4Sand and numgeo floor p at 0.5 kPa; an apparent cohesion p_r·tan φ times N_c ≈ 30–75 means 1 kPa of p_r can move a ~650 kPa load by 3–8 % |
| D2 | **`D_factor`, the UW low-p dilatancy sigmoid** | decide it explicitly; do not inherit it | it is not in Dafalias & Manzari (2004); it acts below p < 0.05·P_atm = 5.05 kPa, i.e. across most of the ring; with `-Presidual 0` it can suppress dilatancy by up to ~900× (the 86 report §6), on a set that already dilates ~8–23× too little (§4 caveat). **There is no deck flag to switch it off today**; if you want the on/off comparison, ask and we add one |
| D3 | **Mesh** | run and report B/4, B/8, B/16, a B/8 sheared 15° and a B/8 jittered 0.1, with R1 on (R2), before calling any q–s mesh-converged | B/16 runs 4–5 % softer than B/8 from s/B 0.010 and the bands are one element wide and follow the mesh lines (§4.8). The cause is **non-associated localization that starts in the hardening regime**, not ψ-softening, so a nonlocal ψ̄ or a crack band would not treat it |
| D4 | **The definition of "limit load" for a dense dilatant sand** | say which you mean: the peak, a plateau, or q at a fixed s/B (your ADR 65 D6) before the runs, not after | no SANISAND arm shows a peak or a plateau to s/B 0.0508, and the present set peaks at 4–16 % strain in an element test (§4 caveat), so a peak, if one exists, needs large settlement; the fixed-s/B criterion is the one that converged in the 90 report |
| D5 | **Re-running campaign curves** | re-run under SAS-ME the ones that feed a reported number, starting with anything read from ring points or near the `ModifiedEuler` wall | §3.1; §4.2 (+6.4 % spurious stiffening near the ME wall) |
| D6 | **The reference load for `NormUnbalance`** | name the vector | it decides whether the integrator's per-point error sits under your Newton tolerance (F18(a)) |
| D7 | **The calibration (T5)** | confirm the campaign SANISAND set against your lab data before any footing curve is used: φ′_peak, strain at peak, dilatancy | it has the strength of a very dense sand but dilates ~20–23× (triaxial) / ~8–14× (plane strain) less than stress–dilatancy requires and peaks at 4–16 % strain (§4 caveat; #892 §10). If your data agree with the set, the physics check moves to your data; if not, recalibration comes first |
| D8 | **Lode parameter c** | keep c ≥ 7/9 (≈ 0.78), or treat extension paths as ill-conditioned | at c = 0.71 the Lode interpolation is non-convex at the extension meridian, and a round-off perturbation decides a CTXu result (§4.9); it is also what makes the footing's wall states singular. **Revised 2026-09-29:** at footing scale, raising c only delays the wall (c = 0.80 walls at s/B 0.048, on the compression side; §0a), so R1 is the route out. Keep c ≥ 7/9 for the Lode convexity, not as a cure for the wall |
| D9 | **R1, a model-level opt-in fix of the wall** | yours to use or not; built as an opt-in `LadrunoSANISAND` variant, **default OFF**, vanilla `ManzariDafalias` behaviour unchanged (owner approved) | see below |

**D9 in detail: R1** (WP-151 memo on #893, **merged 2026-09-29** at `fd87e396d`; flags `-sasHFloor c_A -sasReseatHyst c_rev [-sasSoftCap κ]`,
IntScheme 129 only, all default OFF and byte-identical). Two **coupled** flags:

- **a floor on h everywhere**: h = b0 / max(a, c_A·√(2/3)·m), with c_A ≈ 1 (the yield cone's α-space
  radius, 4.1e-3 at m = 0.005), used in the α update too;
- **a hysteretic re-seat**: α_in := α only when a < −c_rev·√(2/3)·m (c_rev ∈ {½, 1, 2} all pass).

On the material-point oracle, 320 exact increments from 5 real refuser states (E_B, E_D, E_B16), the
pair fails **0 of 320** where DM04 fails **102 of 320**, and **each flag alone fails** (the floor alone 97;
a weaker floor, c_A = ¼, with the hysteresis 90). In element tests at c_A = 1, monotonic q moves by at
most **2.7e-4·q_max** and drained cycles are unchanged to 4 digits. The CTXu difference that was open is
DM04's own extension bifurcation (§4.9), identical with and without R1. R1 is a constitutive change: it
changes the model you calibrated, which is why the decision is yours. It is implemented as an opt-in
(#893, merged at `fd87e396d`; recommended for the campaign). Its footing results are in §0a.

**Regularization, if it is needed after R1.** Whether it is needed is decided by measurement: the
spread of B/4, B/8 and B/16 (R2) against the scatter of comparable physical footing tests, following
the owner's rule, *"physics and real-world behaviour is the judge"* (#892 §9). If it is needed, the
candidate is **Perzyna-type viscoplasticity inside `LadrunoSANISAND`**, not a Duvaut–Lions wrapper: a
Duvaut–Lions update needs the inviscid solution for the same increment, which is exactly the
computation that refuses (#892 §3, §4 R3). A **shear-band width is likely not a meaningful deliverable**
here: the physical band is ~20·d50 ≈ 4–10 mm for d50 0.2–0.5 mm, against elements of 94–188 mm, so
any width a regularized mesh returns would be a numerical length, not the sand's (#892 §9.1 step 5).

### 5.1 What we need from you

1. **The sand**: d50 and the grading, e_max / e_min, and the target relative density D_r.
2. **Lab data**: drained triaxial and plane-strain tests (φ′_peak, the strain at peak, the dilatancy),
   and any undrained cyclic CSR–N target.
3. **The footing test you treat as the reference.**
4. **For your reference footing test:** footing roughness; how the sand was placed relative to the load direction (pluviation or bedding); and the measured unit weight and e_max/e_min of that batch. On Kimura (1985), bedding alone moves the peak by 11 % and its settlement by 1.8×, and roughness can move the peak by 25–35 % (§0a).
5. **The exact PDMY01 33° parameter set** behind your 417.6 kPa control. The fork holds only its own
   WP-133 PDMY03 stand-in (φ 40°), which has no peak in drained plane strain and is not a physical
   reference (#892 §10).

---

## 6. Roadmap

**Near term.**

| item | what | where |
|---|---|---|
| CPPM under Newton + tangent sign | qualify on the strip | **MERGED**, #868 (WP-130), `e96f8d77d` |
| Footing A/B | the Esmeralda arms and the sensitivity ladders are final | #878 (WP-138), draft |
| R1 and c = 0.80 on the footing | done; results in §0a | **MERGED**, #893 (WP-151), `fd87e396d` |
| Tension cutoff (free surface) | fixes from the independent review; p_sep 0 and p_contact sensitivity; a per-point separation census; an energy balance | #894 (WP-152), draft |
| Why the start is too soft | a footing stiffness ladder (seating, K0, surcharge, G0 × 2, h0 × 2) and an element check against Tatsuoka et al. (1986) | Esmeralda and local, running |
| R2/R3 at the peak | B/16 and sheared legs, on Toyoura and on the Kimura case | Esmeralda, running |
| The decision procedure | **R1 → T5 → R2 → T6 → decisions** (below) | #892 (WP-150 memo §9.1), draft |
| F19 step 2 | the deferred message buffer, the counters, stack scratch; allowlist `IntScheme 1`, then 2 and 129; identity and speed-up at 1/2/4/8 threads on a deck that prints | after #868 |
| PDMY hang | refuse a wild trial instead of ~1e9 substeps | **MERGED**, #874 (WP-135) |
| Guide paragraph | "repeatable is not reliable" | **MERGED**, #876 |
| Refusal-aware line search | treat a material refusal as "backtrack", not "rung failed" (119 LineSearch failures on a refusal in E_B) | follow-up WP (#878 §13) |
| SAS-ME follow-ups | its maxSubsteps refusals on the surface ring (209 in E_B) and the re-seat chatter; from its confirmation review: a 64-sample unload-then-reload locator can miss a very short elastic excursion (up to 5.5e-4 on 5 ring cases), the ψ-driven dead end is relocated to ~30 MPa, not removed, positional TolF/TolR are not range-checked | follow-up WP |

**SAS-ME under IMPL-EX: conditional, not qualified.** SAS-ME under the IMPL-EX companion (opt-in
`-implexAllowScheme129`) is **conditional and not yet qualified**. The tangent identity holds, and IMPL-EX
reduces SAS-ME's substep count; its remaining cost is substeps per companion return at tight TolR. But
SAS-ME refuses more often on adversarial states, the IMPL-EX companion has no global Newton to absorb a
refusal, the qualification gates have not been run for `IntScheme 129`, and there is no softening leg.
It stays opt-in; an earlier "no saving" reading is superseded. The opt-in comes from a study still in
progress and is not on `ladruno`: there, `-implex` is still refused with `IntScheme 129` (§2.3).

**The decision procedure** (WP-150 memo §9.1; each step ends in a number and settles the decisions
after it):

1. **R1**, oracle first: the coupled floor + hysteresis on the WP-134 oracle, then the C++, then E_B and
   E_B16 past their walls.
2. **T5**, the material against the sand's physics: element tests at p′ ≈ 10–500 kPa in plane strain and
   triaxial, checked against Bolton (1986) and **your lab data** (D7). The footing cannot be more right
   than its element.
3. **R2**, the footing mesh study with R1 on: B/4, B/8, B/16, plus a B/8 mesh sheared 15° and a B/8
   jittered 0.1, to s/B 0.15. The deck patch is on #892 (`--mesh b4`, `--mesh-perturb shear:15`,
   `--mesh-perturb jitter:0.1`; the default path is byte-identical). R2 is ready to launch once R1 has a
   commit.
4. **T6**, the footing against real footings, against Perkins & Madson (2000), Loukidis & Salgado
   (2011), Lau & Bolton (2011), Vesić (1973) and the De Beer / Tatsuoka scale effect (#892 §11):
   - capacity: the band [q_u(φ′_cs,ps), q_u(φ′_ps,peak)] from your lab φ′ (§4.10);
   - settlement: dense sand in general shear develops the full mechanism at s/B ≈ 6–8 % (the JGGE
     closure, doi:10.1061/JGGEFK.GTENG-12726); where no clear peak forms, Vesić's rule takes q at
     s/B = 10 %;
   - failure mode: a wedge plus a radial shear zone reaching the surface, against the present straight
     vertical bands, judged first on the sheared-mesh leg;
   - band thickness: 10–20·d50 is sub-grid here, so a shear-band width is expected **not** to be a
     deliverable.
   The tolerance (the scatter of comparable footing tests) has **not** been extracted yet; no number is
   claimed for it.
5. **Decisions**: band width by scale (20·d50 against the element size), the tolerance by the scatter of
   comparable footing tests, R3 (Perzyna in SANISAND) only if the R2 spread exceeds it.

**Performance.** Material cost is (global iterations) × (cost per point per iteration), everywhere in
the plastic zone: every OpenSees iteration re-integrates the whole increment from the committed state.
So the levers, in the order we would pull them:

1. **Fewer global iterations**: a step-size policy that targets a few iterations per step. The
   consistent tangent is **not** that lever on this deck (TanType 1 collapses the step, §4.4).
2. **Threading** once F19 step 2 lands: bounded near 3.9× at 8 threads by your own 85.6 % share
   (inference).
3. **An allocation-free kernel**: SAS-ME builds dozens of heap vectors per substep; a fixed-size stack
   kernel is typically several times faster for 6-vector arithmetic. Unmeasured; it must be
   benchmarked before a number is quoted, and it is also a precondition for threads to scale.
4. **A stiff-point fallback** SAS-ME → CPPM: a smaller lever on this deck than it first looked (§4.4).

**TolR 1e-3 is off the list**: on the footing it gave no saving and an earlier wall (§4.4).

GPU or SIMD batching and surrogate models are **not** levers for this problem: the per-point work is
wildly unbalanced (tens to 10⁵ substeps per point per step), and a validation study must not trade the
full-order answer for speed.

**Longer term: a model that is consistent by construction.** DM04 SANISAND has no energy function for
α, so non-negative dissipation is not guaranteed, and the α_in memory and the h ∝ 1/((α − α_in):n)
singularity are exactly where the bugs, and now the wall, live. Two classes are planned (tags reserved,
no code yet):

- **WP-144 `LadrunoNORSAND`** — NorSand in the Borja & Andrade (2006) form, with DM04's power-law
  critical-state line (so ψ stays bounded as p′ → 0) and a Lode-angle dependence for plane strain;
  hyperelastic, implicit, closed-form tangent, proven non-negative dissipation. About 5–6
  engineer-weeks. **This is the route if T5/T6 show that the model, not the calibration or the mesh, is
  the problem** (#892 §9.1 step 7).
- **WP-145 `LadrunoHySAND`** — the 2026 hyperplastic multisurface sand model: the most rigorous
  thermodynamics and the best cyclic behaviour, the least evidence and no public code; a research
  build of about 10–14 weeks, for the later cyclic SSI work.

One expectation to set now: **no dilatant sand model is truly variational.** Non-associated
dilatancy gives a non-symmetric tangent and no incremental minimum principle. "Thermodynamically
admissible, robust implicit return, consistent tangent" is achievable; "variational" is not. And a
p′-floor or a surcharge is established practice even for the rigorous models, so D1 does not go away.

---

## 7. References

**Intake and plan**
- `Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md` + attachments
  `_tims_2d_model_requests_2026-09-25/` (README, `ring_points_b8.csv`, `ring_points_b16.csv`)
- `Ladruno_implementation/127_tims_2d_requests_plan.md` (findings A–D, work packages)

**Evidence documents**
- `128_sanisand_ring_trace.md` — WP-128, #869 (with its 2026-09-27 correction note after WP-134)
- `134_sanisand_reference_integrator.md` — WP-134, #872 (the oracle; U1–U10)
- `131_sanisand_threaded_inventory.md` — WP-131 step 1, #865; E9 re-verified after WP-129, #888
- `133_pdmy_notes.md` — WP-133, #866
- `136_flip_test_drift.md` — WP-136, #870
- `138_footing_sas_me_ab.md` — WP-138, #878 (draft; on `origin/wp/138-footing-sas-me-ab` at
  `1f22e2bad`): §0 verdict, §5.1 curves, §5.3 replays, §8 Esmeralda arms, §9 B/16, §10 diagnosis,
  §11 ladders (final), §12 default, §13 follow-ups; run records under
  `Ladruno_files/testbed/footing_sas_me_ab/`
- `151_sanisand_reseat_singularity.md` — WP-151, #893 (draft; on
  `origin/wp/151-sanisand-reseat-singularity` at `6e9330a3d`, build `bd93c558d`): §2.2 the Zeno
  re-seat accumulation, §2.5 the non-convex extension side and the c = 0.80 wall-fan test, §5 the fan,
  §6.1–6.2 calibrated behaviour, §6.3 the CTXu gate, §8 recommendation, §9 the flags and C++ gates
- `150_sanisand_regularization_memo.md` — WP-150, #892 (draft; on
  `origin/wp/150-sanisand-regularization-memo` at `2a82e2046`, §1.4 corrected with the final ladders; the T5 figures as of `e14703ca7`): §1.4
  mechanism and the R1 oracle box, §2 localization, §3 options, §4 staged recommendation, §8 test plan
  and the R2 mesh-perturbation patch, §9.1 decision procedure, §10 T5, §11 T6 capacity bands
- `_sand_model_survey_2026-09-27.md`, `_sanisand_external_survey_2026-09-27.md`,
  `144_ladruno_norsand_plan.md`, `145_ladruno_hysand_plan.md` — #885

**Guides**
- `LadrunoSANISAND_implex_guide.md` §6.2 (`substepStats`), §6.3 (replay), §13 (choosing an IntScheme;
  SAS-ME); §9 "IntScheme 2 under a global Newton" on `origin/wp/130-sanisand-cppm-under-newton` (#868)
- `75c_pardiso_solver_recipe.md` Trap 7 and "The deterministic mode" (#864; follow-up paragraph #876)

**Ledger rows** (`LEDGER_quirks.md`): the ring CSV convention (finding A); force-accept at `dT_min`
(finding C); F (negative denominator as elastic); the uncapped step factor; G (α_in once per increment);
`ModifiedEuler` TanType-2 chain (T vs dT); RK45 `dAlpha3/4`; IntScheme 4 non-determinism; U9; U10; the
err = 0 path; stress-only error (E); `Stress_Correction`'s silent give-up; RK45 is IntScheme 45 and not a
reference; the 1 kPa floor and stability-limited cost; the flip-determinism pins; the WP-132 CNR rows;
PDMY03 constants and reallocation; the PDMY crossing detector; static `pAtm`; a refusal under a
discarding element was committed. `LEDGER_implementations.md` rows WP-127, WP-129, WP-132, WP-133,
WP-134.

**PRs.** Merged: #846, #863, #864, #865, #866, #869, #870, #871, #872, #874, #876, #884 (Windows-only
CI gap), #885, #888. Open (drafts): #878 (WP-138), #892 (WP-150), #894 (WP-152), and this report, #887. Merged since the fourth issue: #868 (WP-130), #893 (WP-151).

**Literature cited in §4–§6.** Dafalias & Manzari (2004), *J. Eng. Mech.* 130(6); Bolton (1986),
*Géotechnique* 36(1); Rudnicki & Rice (1975), *JMPS* 23; Vesić (1973), *JSMFD* 99(SM1); Perkins &
Madson (2000), *JGGE* 126(6); Loukidis & Salgado (2011), *Géotechnique* 61(2); Lau & Bolton (2011),
*Géotechnique* 61(8); Tatsuoka et al. (1991), ASCE GSP 27; Martin (2005), *Proc. 11th IACMAG*; Han et
al. (2016), *SpringerPlus* 5, 1482. Full list in the WP-150 memo.

**Earlier replies to the TIMs team.** `86_ladruno_sanisand_tims_report.md` (the hidden cohesion),
`90_ladruno_regularization_tims_report.md` (regularization, `-maxSubsteps`).

---

## Revision log

| date | change |
|---|---|
| 2026-09-28 | First issue. F18(a), (b), (e), F20, F21, F22 (mode), F23(a), (b) answered from merged work; F18(c), (d), F19 step 2, the F22 guide paragraph and the WP-138 footing A/B pending; placeholders in §4. |
| 2026-09-28 | Second issue. §4 final from the WP-138 Esmeralda arms (#878 at `762be8332`): the wall table, the verdict (SAS-ME moves the wall from s/B 0.0292 to 0.0508 and does not remove it; no peak; constitutive), accuracy and cost, the replay figures reconciled, the mechanism (#892), the interim sensitivity ladders (16:20) with caution 6 re-measured, B/16 and non-associated localization, the calibration caveat (#892 §10), Lode convexity c ≥ 7/9, the classical capacity bands (#892 §11 at `00198f278`), the integrator recommendation; placeholders removed. §5: D1 and D3 revised, D7–D9 added (calibration, c, R1), the regularization route and §5.1 "What we need from you". §6: the WP-150 decision procedure, T6 targets, the SAS-ME + IMPL-EX status, #874 and #876 merged, TolR 1e-3 and TanType 1 off the performance list. §0, §2 scope note, §3 cautions 1, 4, 6 and the E_B configuration line (TolR 1e-4, not 1e-7) updated to match. |
| 2026-09-28 | Third issue. §4.7 ladders FINAL (#878 at `1f22e2bad`): every leg walls on `loadingNonPosH`; dilatancy off only delays the onset (0.0363 → 0.0426), so the interim "only killing the dilatancy clears it" is withdrawn. §4.6: the wall states need the non-convex extension side (c = 0.71 < 7/9; c = 0.80 takes the wall fan 102/320 → 0/320, WP-151 §2.5), the wall and the bands are separate phenomena, two routes out (R1 or c ≥ 0.78), and related literature on reversal-memory stress overshooting. R1 and the CTXu finding now cite the WP-151 memo (#893) instead of "preliminary". D8 and the roadmap updated; the R1, c = 0.80 and R2 footing runs are running. |
| 2026-09-28 | Fourth issue. §4.6: the Chen et al. (2022) citation corrected. It is a stress-overshooting study, not a documented SANISAND footing that stops on this singular set; Ghorbani et al. (2023, memory repositioning) added as related literature. |
| 2026-09-29 | Fifth issue. New §0a: R1 merged (#893) with its footing results and κ as a guard; CPPM merged (#868); the c ≥ 0.78 route withdrawn (c = 0.80 walls at s/B 0.048); the campaign set identified as a cyclic fit, and a physically bounded monotonic set (PB2) offered; the free surface as the limiter after R1; the WP-152 tension cutoff (#894, pending review); Gate 1: DM04 Toyoura against Kimura (1985) Fig. 9, digitized, with the peak consistent within ±15 %, a late peak and a soft, concave-up start; mesh case C. §0 item 8, D8, D9, §4.6 route 2, §5.1 and §6 updated. |
| 2026-09-29 (night) | Fifth issue, addendum. §0a item 9: the element check of DM04 against Tatsuoka et al. (1986), with its reading for the footing and the ladder's first read. Item 5: the cutoff re-labelled as a low-confinement separation (E2 is the operative trigger; p_sep = 0 is not viable); every cutoff number marked pre-review (`8ebde5cbd`). |
| 2026-09-29 (late) | §0a item 10: the acoustic census answers R2/R3 as case C (the mesh gap and loss of ellipticity at s/B ≈ 0.033, 5× before the peak; mesh-imposed band orientation; non-associativity is the cause); R3b (Perzyna in-model) is the planned regularizer, tuned after the constitutive questions. |
