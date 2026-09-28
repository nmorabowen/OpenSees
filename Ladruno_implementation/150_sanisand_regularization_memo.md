---
title: "WP-150 — Regularizing SANISAND on the TIMs strip footing: two mechanisms, and a staged design"
project: Ladruno
type: design memo
status: "OWNER DECIDED 2026-09-28: D-a YES (R1 opt-in, oracle-first); D-b/c/d by measurement against physical evidence (§9.1). R1 in progress; R2/R3 not started."
owner: nmora
related:
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
  - "[[134_sanisand_reference_integrator]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[59_ladruno_gradient_concrete_adr]]"
  - "[[_sand_model_survey_2026-09-27]]"
  - "127_tims_2d_model_report (PR #887)"
  - "138_footing_sas_me_ab (PR #878)"
tags: [memo, sanisand, regularization, localization, ellipticity, tims, wp-150]
updated: 2026-09-28
---

# WP-150 — Regularizing SANISAND on the TIMs strip footing

> [!summary] The short version
> The WP-138 Esmeralda legs stop for **two different reasons**. Only one of them is what a regularizer treats.
>
> 1. **The wall (`loadingNonPosH`) is a singularity of the DM04 rate equations, not band softening.**
>    - No committed Gauss point is anywhere near H = 0: min H/2G = 0.92 at E_B's wall.
>    - The refusing points sit **at an α_in re-seat** (a = (α−α_in):n ≈ 0, h → ∞), **near the bounding surface** (b:n small).
>      There, Kp = ⅔·p·(b0/a)·(b:n) is ∞·0, and it goes to −∞ as soon as b:n ≤ 0.
>    - The same 0/0 stopped WP-134's exact Radau oracle at ring point 1950/3.
>    - *Corrected 2026-09-28 by the R1 session's oracle:* the exact rate equations reach this point through a
>      **Zeno accumulation of re-seats**, with b:n → 0 **from above**. H ≤ 0 is only the b:n < 0 exit.
>    - No regularization of the boundary-value problem can lift an unbounded negative modulus at a point.
>      The cure is a model-level fix (**R1: the h floor everywhere plus a hysteretic re-seat, coupled**).
> 3. **The campaign parameter set violates stress–dilatancy** (T5, §10).
>    - Plane-strain φ′_peak is 45–60° with a peak dilation angle of only 0.8–1.9°, and the peak comes at 4–16 %
>      axial strain.
>    - A lab-calibrated DM04 set (Toyoura) at the same implied density peaks at 1–5 % with ψ_max 14–23°.
>    - Before any footing curve is called physical, TIMs' calibration must be checked against their sand's data.
> 2. **The bands are non-associated localization in the HARDENING regime.**
>    - At s/B 0.011, 17 % of the Gauss points already have a singular acoustic tensor (det ≤ 0).
>      At those points H/2G ≈ 1.05 and Kp/2G ≈ 0.03, i.e. still hardening.
>    - The **same state with associated flow is elliptic everywhere**.
>    - Only 17 of 9720 points are post-peak even at s/B 0.0508.
>    - This is Rudnicki–Rice / Sabet & de Borst "structural softening", and it is what makes B/16 run 4–5 % softer than
>      B/8 from s/B 0.010.
>    - A regularizer that acts on ψ-softening (nonlocal void ratio) or on fracture energy (crack band) does not touch it.
>
> **Recommendation.** R0 (confirm, no code) → **R1** (opt-in h floor everywhere + hysteretic re-seat, coupled) →
> **R2** (unregularized B/4–B/8–B/16 re-measure past 0.05) → **R3 only if R2 fails TIMs' tolerance**.
> R3 is Perzyna-type viscoplasticity *inside* `LadrunoSANISAND`. Rejected: the Duvaut–Lions wrapper (it inherits the
> refusal), and nonlocal ψ̄ for this stage (it acts on the wrong mechanism). Cosserat or gradient-in-λ only if a band
> **width** deliverable is ever named.

Reproducers: `Ladruno_files/testbed/hypo_bearing/wp150_regularization/`. The scripts read the WP-138 Esmeralda
checkpoints read-only (`field_*.npz`: committed stress, state, ψ and SAS counters per GP) and mirror the kernel
formulas cited below.

---

## 1. The wall: a singular set of the rate equations

### 1.1 The refusal criterion

The SAS-ME stage (`LadrunoSANISANDSasME.cpp:335-411`) refuses a *loading* stage (N > 0) whose plastic denominator
is not positive:

```
H = Kp + 2G - K·D·qv            (B - C tr n³ ≡ 1, WP-134 §2 derived fact (b); qv = n:α + √(2/3)·m)
Kp = ⅔·p·h·(b:n),   h = b0 / a,   a = (α − α_in):n,   a < 1e-10 → h = 1e10   (ladrunoSasBracketH :162-172)
```

At stage 1 H depends only on the substep's START state. A reversal found there re-seats α_in := α
(`:568-572`), so a = 0 and h = 1e10. The refusal code is RC_NONPOS_H (`:578-582`), and "no cut can change it".

### 1.2 No committed state is close to H = 0

`h_decomp.py` evaluates H term by term at every committed GP:

| checkpoint | s/B | GPs | H ≤ 0 | min H/(2G) | post-peak (b:n < 0) | h at the 1e10 sentinel |
|---|---|---|---|---|---|---|
| E_B step 215 | 0.0407 | 9 720 | **0** | 1.039 | 59 | 7 |
| E_B step 370 | 0.0508 | 9 720 | **0** | 1.038 | 61 | 90 |
| E_B last converged | 0.0508 | 9 720 | **0** | 0.917 | 17 | 99 |
| E_B16 last converged | 0.0135 | 38 880 | **0** | 1.036 | 7 | 1 635 |

Genuine softening is marginal. For the 17 post-peak GPs at E_B's wall, a would have to shrink ≥ 6.7× before
H ≤ 0 (x_crit/x ≤ 0.149).

### 1.3 The refusing points sit at a re-seat, near the bounding surface

The points are from `floor_refusers_in_band.csv`. The state is the last committed one; the counters are cumulative
for that GP (`refuser_stats.py`):

| leg | ele/gp | (x, y) m | p kPa | ψ | a = (α−α_in):n | h | b:n | α_in re-seats | rejected reversals |
|---|---|---|---|---|---|---|---|---|---|
| E_B | 1880/1 | (−0.148, −2.023) | 366 | −0.089 | +2.2e-4 | 2.8e5 | +0.193 | 1 237 | 7 162 |
| E_B | 1879/1 | (−0.148, −2.210) | 334 | −0.091 | −1.6e-3 | **1e10** | +0.030 | 880 | 9 075 |
| E_B16 | 7820/4 | (+0.957, −0.395) | 50 | −0.116 | −3.3e-3 | **1e10** | +0.046 | 602 | 11 109 |

- **All three have |a| below the yield cone's α-space radius √(2/3)·m = 4.1e-3.**
- All three are dense (ψ ≈ −0.1) and inside the bounding surface, with b:n > 0 and ρ_α ≈ 0.95.
- Domain totals: E_B has **10.1 M re-seats and 88.5 M reversal-rejected substeps**; E_B16 has 5.7 M and 63.6 M.
  The band points chatter between "loading" and "reversal" from one Newton iterate to the next.

### 1.4 The mechanism

WP-134 already gives the exact continuous extension (`134_sanisand_reference_integrator.md` §2, derived fact (c)):

```
H_s = a·H = ⅔·p·b0·(b:n) + a·(2G − K·D·n:r)      at a = 0:  H_s = ⅔·p·b0·(b:n)
```

The loading index L = a·N/H_s is 0/0 when b:n → 0, and there is no solution when b:n < 0. So the singular set is
**{a = 0, b:n ≤ 0}**: a new loading process that starts **on or outside** the bounding surface along the new n.

DM04's h = ∞ at a re-seat encodes an elastic-like restart from *inside* the bounding surface. Near the peak, a small
rotation of n (non-coaxial loading in the edge zone, or iterate jitter) re-seats α_in while b:n' ≤ 0. The model then
asks for an infinitely *negative* modulus.

> [!important] Corrected 2026-09-28 — how the set is reached (the R1 session's oracle, exact Radau)
> **Test set:** 5 real refuser states rebuilt from the checkpoints (E_B 1880/1, 1879/1; E_D 1962/1, 2058/1; E_B16
> 7820/4) × 32 directions × 2 magnitudes = 320 increments.
>
> **Mechanism: a Zeno accumulation of re-seats, with b:n → 0 from above.**
> - After a re-seat, h = ∞ makes α slide along b.
> - With b nearly perpendicular to n, that slide rotates n on the m = 0.005 cone, so a < 0 again and the next
>   re-seat follows.
> - At 1880/1 the re-seat intervals are 1.7e-2, 1.5e-3, 5.6e-5, 2.4e-6, … They accumulate at a finite t*, where
>   a = 0, b:n = 1.9e-8 > 0 and |dα/dt| ≈ 2e5.
> - H ≤ 0 is only the b:n < 0 exit of that sequence.
> - The refusers sit with n at an extension-side Lode angle: ρ_b(θ_n) = 1.17–1.31 while ρ_α = 0.93–0.96.
>
> **Consequence: the "b:n ≤ 0 only" floor of the first R1 draft does not work.**
>
> | variant | failures (of 320) |
> |---|---|
> | DM04 | 102 |
> | "b:n ≤ 0 only" floor | 102 |
> | + hysteresis | 102 |
> | floor everywhere alone | 97 |
> | floor everywhere (c_A = 1) **+** re-seat only when a < −c_rev·√(2/3)·m | **0** (min H/X 0.17), for c_rev ∈ {½, 1, 2} |
>
> The dilatancy reading of the ablation below is therefore weakened:
> - the sequence also runs with b:n > 0;
> - S4 carries ~⅓ of E_B's load at the same s/B;
> - ψ-contraction may feed the set, but it is not required.

**Corroboration.**
- WP-134's Radau oracle **stopped** at ring point 1950/3 under shear, where "(α−α_in):n → 1e-10 and b:n → 4e-7
  together" (134 §6.6, row at `:462`). So the singularity is in the DM04 continuum, not in SAS-ME's discretization.
- The WP-138 ablation (provisional, B/8, SAS-ME):

  | leg | NonPosH events | first at s/B |
  |---|---|---|
  | S1 (no fabric) | 8 | 0.0355 |
  | S2 (no peak, nb = 0) | 22 | 0.0269 |
  | S3 (nd = 0) | 1 | 0.0347 |
  | S4 (A0 = 0.001) | **0** | — (run to 0.0374) |

  - Every S1–S3 leg sits at ~⅓ of E_B's load at the same s/B.
  - **Reading:** dilation raises ψ in the band, which contracts the bounding surface onto α. That is how b:n reaches
    0⁻. Without dilatancy, the bounding surface stays an attractor (dα ∝ b) and b:n > 0.
  - h0 × 3 walls much earlier (0.0091): α reaches the bounding surface sooner.
  - `-Presidual` 0.5–20 kPa and e_init 0.65–0.85 all still hit it.
  - The onset is **non-monotonic** in A0 and in Presidual. That is the fingerprint of a singular event, not of a
    smooth limit.

**Owed (R0).** A replay of the *refusing* substep, to show a ≈ 0 and b:n' ≤ 0 at the refusing stage. Use the WP-138
FixedNumIter-1 recipe plus `ladrunoSANISANDReplay -trace` on the SAS-ME snapshot binary.

---

## 2. The bands: non-associated loss of ellipticity at positive hardening

### 2.1 Acoustic tensor of the continuum tangent

`acoustic_vec.py` builds, at every committed GP, D = Cₑ − (Cₑ:R)⊗(Q:Cₑ)/H with Q = n − ⅓·qv·I and the kernel's R.
It then scans the plane-strain acoustic tensor over 721 band orientations. **Associated control:** R → Q, same Kp.
**Blend:** (1−β)·Cₑ + β·D, as in ADR-90 V4.

| checkpoint | s/B | GPs with det ≤ 0 | H/2G there, median (min) | Kp/2G there, median | min det ratio | associated det ≤ 0 | blend det ≤ 0 at β = 0.5 / 0.9 / 0.99 |
|---|---|---|---|---|---|---|---|
| E_B step 45 | 0.0110 | 1 636 (16.8 %) | 1.074 (1.052) | +0.032 | −0.156 | **0** | 0 / 91 / 1 487 |
| E_B step 120 | 0.0294 | 1 924 (19.8 %) | 1.067 (1.047) | +0.022 | −0.172 | **0** | 0 / 199 / 1 691 |
| E_B last | 0.0508 | 2 187 (22.5 %) | 1.060 (0.917) | +0.031 | −0.478 | 1 ¹ | 0 / 588 / 1 949 |
| E_B16 step 55 | 0.0096 | 3 534 (9.1 %) | 1.076 (1.050) | +0.035 | −0.248 | **0** | 0 / 82 / 3 285 |
| E_B16 last | 0.0135 | 6 576 (16.9 %) | 1.072 (1.041) | +0.033 | −0.318 | **0** | 0 / 1 728 / 6 016 |

¹ The footing-edge surface GP at p = 9 kPa, the one genuinely softening point (Kp/2G = −0.16).

**Reading.**
- Ellipticity is lost while every one of these points is still **hardening** (H/2G ≈ 1.05). The associated flow at
  the same Kp is elliptic.
- The cause is the flow rule. R has a volumetric part D/3 ≈ −0.01 against Q's −qv/3 ≈ −0.5.
- At the campaign's h0 = 1.3 the plastic modulus is small against G: Kp/2G ≈ 0.03 once α has travelled a ≈ 1 from α_in.
- This is Rudnicki & Rice (1975), and in FE terms Sabet & de Borst (2019): structural softening from non-associated
  flow with no strain softening at all.
- **Mesh-orientation warning.** The bands in `fields_shear_strain.png` run straight down the element columns at
  x = ±B/2. The late inclined band at E_B's wall runs along element diagonals. An ill-posed problem picks mesh-aligned
  paths, so R2 must include an orientation variant.

### 2.2 What it costs the curve

E_B (B/8) against E_B16 (B/16), q at matched s/B:

| s/B | 0.001 | 0.002 | 0.004 | 0.006 | 0.008 | 0.010 | 0.012 | 0.0135 |
|---|---|---|---|---|---|---|---|---|
| (q16 − q8)/q8 | −0.53 % | −0.58 % | +0.16 % | −0.81 % | −0.91 % | **−3.89 %** | **−4.95 %** | **−3.84 %** |

The two meshes agree to 1 % until the band forms at B/16 (~0.009), then separate by 4–5 %. GATE U (ADR-90 §1.2) saw
the matched-settlement band *contract* on the 3-D deck. Whether this one contracts needs B/4 and B/32-class points (R2).

---

## 3. Options compared

The physical band in sand is 10–20·d50 ≈ 3–6 mm (Mühlhaus & Vardoulakis 1987; Desrues & Viggiani 2004). The
elements here are 94–188 mm. **Any intrinsic length at this scale is a declared numerical ℓ, not the soil's.**
ADR-90 §3.3's honest-framing test applies to ℓ exactly as it does to τ: never tune it against a target load.

| option | regularizes | band WIDTH objective? | lifts the wall (H ≤ 0)? | ellipticity at H > 0 (non-assoc.)? | fork fit | sand-footing evidence | verdict |
|---|---|---|---|---|---|---|---|
| **Duvaut–Lions** (ADR-90 wrapper `LadrunoOverstress` 33022, or WP-F in-model) | rate | no (quasi-static: width is De- and imperfection-set, ADR-90 A0/V3) | **no** — needs the inviscid σ̄ for the same Δε, i.e. the refusing computation | yes, for β ≲ 0.7 (V4; §2.1) | wrapper: V1–V3, stale Cₑ, β per Newton iterate (ADR-90 §4–5); in-model: no closest-point projection exists for a bounding-surface model with a tiny cone | none | **reject** |
| **Perzyna-type viscoplasticity inside LadrunoSANISAND** | rate | no (as above; the length comes only through wave dispersion, de Borst & Duretz 2020) | **yes for bounded Kp** (denominator H + 2Gτ/Δt); not for the ∞·0 set, so it needs R1 first | **yes**: §2.1 needs H_v/2G ≥ 0.17 (s/B 0.011) … 0.44 (0.0508) | intrinsic: the real Cₑ(p) (fixes V1), +2 doubles on the wire, per-instance, latched Δt (ADR-90 D4), refuse with `-implex` (Concrete3D `-eta` precedent) | viscosity regularizes non-associated plasticity (Hageman, Sabet & de Borst 2021) | **R3, conditional** |
| **Nonlocal void ratio / ψ̄** (Gao, Li & Lu 2022; Mallikarachchi & Soga 2020) | ψ-softening | yes, when h < ℓ | no | **no** — it changes how H *evolves* in the band, not the current tangent; the onset here is at Kp > 0 | lagged integral average via a domain component, no new DOF; MP halo needed | **yes**: strip footing (Gao), biaxial (M&S) | **R4, post-peak only** |
| Implicit gradient (u–ē element, ADR-59 architecture) | the averaged variable | yes | no | only if the averaged variable is the plastic multiplier | new ndf-3 element family, mixed-ndf seams, conditioning (ADR-59 R-M4/M5) | clays / concrete | not now |
| Gradient plasticity in λ (de Borst & Mühlhaus 1992) | the plastic multiplier field | yes | partly | yes | mixed u–λ element, research for SANISAND | granular (Vardoulakis & Aifantis 1991) | not now |
| Cosserat / micropolar | rotation gradients | yes (ℓ_c) | no | **yes** (Sabet & de Borst 2019; Hageman et al. 2021) | new element family (ux, uy, ω), no Cosserat SANISAND in the literature, none in `SRC/` | micropolar hypoplastic footings (Tejchman) | only if a WIDTH deliverable is named |
| Element lch / crack band (Pietruszczak & Mróz 1981; Siddiquee et al. 1999) | softening energy | no (1 element) | **worse**: w_phys ≪ h steepens the element's softening (snap-back) | no — there is no fracture energy before the peak | the fork's lch seam exists | yes (scale effect, Tatsuoka group) | not applicable |
| **R1 — model fix of {a = 0, b:n ≤ 0}** | nothing (not a regularizer) | — | **yes** | no | intrinsic, SAS-ME bracket only, default off | — | **do first** |

---

## 4. Recommendation — staged, each stage with its gate

**R0 — confirm the mechanism (no code, about 1 day).**
- Replay the three refusers' failing substep with the WP-138 recipe.
- **Gate:** a ≈ 0 re-seat with b:n' ≤ 0 at the refusing stage in all three.
- If it fails, §1 is wrong and this memo is re-scoped before any code.

**R1 — the h floor everywhere + a hysteretic re-seat, as COUPLED flags (opt-in; owner D-a YES).**

**Owned by the "SANISAND α_in re-seat singularity fix" session, oracle first.** This section records the spec as
corrected by that session's oracle (§1.4 box). The first draft here, a floor only where b:n ≤ 0, fails as often as
DM04.

- **Floor everywhere:** h = b0 / max(a, a_min), with **a_min = c_A·√(2/3)·m** (c_A ≈ 1 → 4.1e-3). The α update uses
  the same h.
- **Hysteretic re-seat:** α_in := α only when a < −c_rev·√(2/3)·m (c_rev ∈ {½, 1, 2} all pass).
  - This is half of the well-posedness fix, not only a cost lever: it plausibly explains the 88 M rejected reversals
    and S4's 26.9 M substeps.
  - It is not WP-129's strain-norm `-reversalTol/-reversalRel`, which IntScheme 129 refuses.
- **Each part alone fails:** floor alone 97/320, c_A = ¼ + hysteresis 90/320. Together: 0/320, min H/X = 0.17.
  An optional softening cap H ≥ ½X also gives 0/320 and never activates in element tests.
- **Element-test cost at c_A = 1:**
  - monotonic |Δq| ≤ 2.7e-4·q_max;
  - drained cycles identical to 4 digits;
  - undrained cyclic N unchanged **except one open outlier, a BLOCKING gate**: CTXu, e0 0.6944, CSR 0.2. DM04 reaches
    5 % DA at N = 8; the variants do not by N = 20. DM04 itself breaks axisymmetry in extension there (c = 0.71 < 7/9).
- **Ablation legs on the footing:** floor alone, hysteresis alone, both, and both + cap.
- **Gate:**
  1. The oracle suite above, including the CTXu outlier.
  2. E_B and E_B16 with both flags pass s/B 0.0508 / 0.0135 with no singular-set refusal, and q–s is unchanged below
     s/B 0.036.
  3. c_A and c_rev ∈ {½, 1, 2} move q at matched s/B by less than the solver floor (0.8–1.4 %, ADR-90 §1.2(iii)).

**R2 — the unregularized re-measure (no code; Esmeralda).**
- Run B/4, B/8, B/16 and one skewed or unstructured B/8 fine zone, with R1, to s/B 0.15.
- Report q(s/B), the §2.1 acoustic census every 5 checkpoints, w₂ (ADR-90 §7.3) and the band paths.
- **Decision rule:** if the matched-settlement band at s/B ∈ {0.05, 0.10} contracts under refinement and sits inside
  TIMs' tolerance (ADR-90 OQ2, still unsupplied), the answer is **disclose** (ADR-90's close-out stands).
  Otherwise go to R3.
- On present evidence (−4 to −5 % B/16 vs B/8 at 0.010–0.013), R3 is likely if the tolerance is below about 5 %.

**R3 — Perzyna-type viscoplastic SANISAND (conditional; C++ behind a default-off flag).**
- Rate form λ̇ = ⟨f⟩ / (2G·τ), with f = SANISAND's own yield function (the stress's distance outside the α-cone).
  The backward-Euler denominator is **H + 2G·τ/Δt**.
- Why Perzyna and not Duvaut–Lions: Perzyna **never needs an inviscid solution**, so it cannot inherit the refusal.
  It uses the real Cₑ(p). It is an ODE in pseudo-time, so SAS-ME's error control integrates it unchanged.
- η = 2G·τ (∝ √p) keeps the Deborah number uniform over the domain.
- **The trade-off, stated up front.**
  - The incremental problem stays elliptic only while τ ≳ 0.2–0.45·Δt_step (§2.1: H_v/2G = τ/Δt).
  - So either the deck caps ds at about 2τ (cost up to ~10×: E_B's mean ds was 2e-4 m against a 2e-5 base), or τ
    grows and so does the rate bias.
  - The bias estimate is f ≈ 2G·τ·λ̇: about 7 kPa in the band at τ = 1e-5 m, and it must be measured.
  - There is **no intrinsic width** (quasi-statics, ADR-90 A0/V3). R3 claims a q–s that converges in h at a declared
    (τ, Δt), never a width.
- **Gate:** C8 (the algorithmic acoustic tensor is elliptic at every committed state); q(s/B) at fixed τ contracts in
  h; Δt-convergence at fixed τ; q(τ)/q(τ → 0) reported per leg at {τ/2, τ, 2τ}; zero steps committed with τ > 0 but
  Δt = 0 (the revert path, ADR-90 §4.2).

**R4 — nonlocal ψ̄ (only if a post-peak branch appears in R2/R3).**
- Gao-type: a nonlocal volumetric-strain increment drives e, averaged with a lag by a domain component, with a
  declared ℓ ≥ 3·h_coarse.
- It is the right tool once ψ-softening dominates, not before.

**Not recommended now:** Cosserat and gradient-in-λ. They are the only routes to an objective *width* for this
mechanism, and at this scale that width is numerical anyway. Open them only on a named width deliverable (ADR-59's gate).

---

## 5. Parameters and how to justify them

| parameter | stage | what it is | calibration | forbidden |
|---|---|---|---|---|
| c_A (a_min = c_A·√(2/3)·m) | R1 | smallest α travel since the last reversal that the memory resolves | tied to the yield cone, so no fit. Default 1, sensitivity {½, 1, 2}. Every refuser had \|a\| < √(2/3)·m; c_A = ¼ fails | fitting c_A to a load |
| c_rev | R1 | re-seat hysteresis, in cone radii | ½, 1 and 2 all pass the oracle set; the owner of R1 picks the default | as above |
| τ (m of settlement; η = 2G·τ) | R3 | relaxation time in pseudo-time | the smallest τ that keeps C8 elliptic at the deck's ds_max (§2.1 gives τ/Δt ≥ 0.45 with margin), then {τ/2, τ, 2τ} | tuning τ to a target q or width (ADR-90 §3.3) |
| ℓ | R4 | nonlocal radius | ≥ 3·h of the coarsest mesh (numerical; the physical ℓ is ~mm), reported at {ℓ/2, ℓ, 2ℓ} | as above |

---

## 6. At the ring (p′ → 0)

- **R1.**
  - b0 ∝ 1/√p, so the capped Kp = ⅔·p·(b0/a_min)·(b:n) ∝ √p. That is the scaling of 2G.
  - The cap's effect is therefore p-independent, and it adds no strength at the free surface.
- **R3.**
  - η must scale with G. With a constant η the ring (G → G(p_min)) is viscosity-dominated, and the overstress acts as
    an apparent cohesion. That is the same failure as `-Presidual` 20 kPa (+30 % q at s/B 0.05, ablation).
  - The absolute overstress 2G·τ·λ̇ → 0 with G.
- **Where genuine softening starts.** The acoustic scan puts it at the footing-edge surface point (p = 9 kPa,
  Kp/2G = −0.16). The ring is where R3's local-uniqueness role would appear first.
- **Unchanged by all of this:** `-Pmin`, `-Presidual` and the UW D_factor sigmoid (< 5.05 kPa). The survey's p′-floor
  rule (§7.1) still governs them.

---

## 7. Tangent, refusal and fork plumbing

- **Refusals.** R1 removes one refusal class and keeps RC_NONPOS_H for genuine softening. Both still leave through
  `LADRUNO_MATERIAL_REFUSED` to the element and to `analyze()`. `LadrunoQuad` is a forwarder (LEDGER_quirks "Element
  refusal roster"). The WP-99 latch is IMPL-EX commit-time only and is not touched.
- **Tangent.**
  - SAS-ME hands out one continuum tangent at the end state. R1 changes only Kp and h inside it.
  - R3's tangent is Cₑ − (Cₑ:R)⊗(Q:Cₑ)/(H + 2G·τ/Δt) and is non-symmetric, so the unsymmetric solver stays mandatory
    (ADR-90 D11).
  - The campaign's TanType 0 (modified Newton) is unaffected.
- **State and wire.**
  - R1 is two option doubles in `LadrunoSasOptions` (c_A, c_rev). That struct is a Ladruno block in the vanilla
    `ManzariDafalias.h`, so it gets a vanilla-ledger row. The R1 session owns this.
  - R3 adds τ, the committed overstress and the latched Δt: getCopy, sendSelf/recvSelf, revertToLastCommit.
  - No statics, so it is thread-safe for WP-131/146.
- **IMPL-EX.** SAS-ME + IMPL-EX is unqualified (concrete session's study), so R1/R3 are SAS-ME-only. R3 hard-refuses
  `-implex`.
- **Checklist.** Implementation follows `.claude/skills/ladruno-new-material/SKILL.md`.

---

## 8. Test plan (the WP-138 footing)

Deck: `~/ladruno_wp138/deck/footing_ab.py` on Esmeralda (E_B settings: IntScheme 129, TolR 1e-4, Pardiso sequential,
MKL_CBWR COMPATIBLE).

**New flags:**
- `--hcap c_A` and `--tau τ`;
- `--dsmax`, to cap ds for R3;
- `--mesh b4`: the graded counts must be integers, so b4 needs its own counts; it cannot be r = ½ of b8;
- a skewed-mesh variant.

| id | stage | legs | pass |
|---|---|---|---|
| T0 | R0 | replay of E_B 1880/1, 1879/1 and E_B16 7820/4 | a ≈ 0 and b:n' ≤ 0 at the refusing stage |
| T1 | R1 | WP-134 oracle suite ± cap, C++ vs oracle | < 0.1 % monotonic; cyclic change reported |
| T2 | R1 | E_B and E_B16 with c_A = 1, plus c_A ∈ {½, 2} on B/8 | past 0.0508 / 0.0135, no singular-set refusal, q within the solver floor across c_A |
| T3 | R2 | B/4, B/8, B/16 and skewed B/8, to s/B 0.15 | the decision rule of §4 R2; acoustic census; w₂; band paths |
| T4 | R3 | τ ∈ {τ/2, τ, 2τ} × {B/4, B/8, B/16}, ds ≤ 2τ, plus one leg at ds_max/2 | C8 elliptic; h-contraction at fixed τ; Δt-convergence; bias reported |
| T5 | physics, element | plane-strain and triaxial element tests of the campaign set at p′ ∈ {10, 50, 150, 500} kPa, R1 on and off | φ′_peak, the strain at peak and the dilatancy against Bolton (1986) and TIMs' lab data |
| T6 | physics, footing | the T3 curves against the dense-sand footing evidence (§9.1 step 4) | q_u / N_γ, s/B at peak and rupture pattern inside the published ranges; mesh spread against test scatter (D-b) |

**The pass criterion TIMs asked for** (mesh-independent q–s past s/B 0.05 at B/8 and B/16, and B/4) is T3 if R2
suffices, and T4 otherwise.

---

## 9. Decisions — owner, 2026-09-28 (relayed by the TIMs orchestrator)

- **D-a: YES.** R1 is an opt-in DM04 variant, default OFF, and it is built **oracle-first**: the WP-134 reference
  integrator with the cap comes before any C++. a_min is tied to m, as proposed (§5).
- **D-b, D-c, D-d: "measure and decide; physics and real-world behaviour is the judge."** No tolerance, ADR filing or
  width deliverable is fixed up front. §9.1 is the procedure that decides them.

### 9.1 The decision procedure

Each step ends in a number, and each open decision is settled by the numbers of the steps before it.

1. **R1, oracle first** (T1, T2; the R1 session). The spec is the coupled floor-everywhere + hysteresis of §4.
   - The WP-134 oracle suite runs with the cap. The 1950/3 shear row must now integrate *through* the former 0/0.
   - Then the C++, then E_B and E_B16 with c_A = 1 past their walls.
2. **The material against the sand's own physics** (T5).
   - Element tests at the footing's stress range (p′ ≈ 10–500 kPa), in plane strain and triaxial: φ′_peak(ψ), the
     strain at peak, and the dilatancy.
   - Checks: the Bolton (1986) relation φ′_peak − φ′_cs ≈ 5·I_R in plane strain, and TIMs' lab data if they have it.
   - The footing cannot be more right than its element. The campaign set already gives q = 967 kPa at s/B 0.05 and
     still rising, against DP's 824 and PDMY's 418. This step says which one is physical.
3. **R2, the footing mesh study with R1 ON** (T3): B/4, B/8, B/16 and a skewed B/8, to s/B 0.15.
4. **The footing against real footings** (T6).
   - The benchmark numbers:
     - q_u or N_γ, compared with the dense-sand evidence for relative density and stress level: Perkins & Madson
       2000; Loukidis & Salgado 2011; Lau & Bolton 2011; De Beer's and Tatsuoka's scale effect;
     - s/B at peak: dense sand in general shear typically peaks at a few percent of B and punches beyond ~10 % in loose
       sand (Vesić 1973). The range is re-read from the sources at this step;
     - the rupture pattern: a wedge plus radial shear zone, versus the present straight vertical bands.
   - Check each band path against the skewed mesh (§2.1) before comparing it with the physical mechanism.
5. **D-d, band width, decided by scale.**
   - Rule: if 20·d50 is below one third of the finest element, a mesh-converged band width is **not a meaningful
     deliverable** and Cosserat or gradient stays closed. For typical sands (d50 ≈ 0.2–0.5 mm), the band is 4–10 mm
     against h = 94–188 mm.
   - TIMs supply d50. The rule does not depend on the answer unless the sand is gravel.
6. **D-b, the tolerance, decided by test scatter.**
   - The tolerance is the scatter of comparable footing tests at matched s/B (repeat tests, or N_γ scatter at the same
     D_r), taken from the step-4 sources.
   - If the step-3 spread (B/4, B/8, B/16, extrapolated) is **below half that scatter**, the answer is disclose; R2 is
     the result.
   - If it is above, **R3** is built, and the same test decides it: its τ-bias must also sit inside the scatter.
7. **D-c, the filing, decided last.**
   - R1 alone suffices → an **ADR-86 follow-up** (model option, no class tag).
   - R3 needed → **ADR-90 revision** (retire the wrapper; WP-F becomes in-model Perzyna; 33022 stays reserved unused).
   - Either way, no new ADR number unless step 4 shows the model itself is the problem. If it does, the answer is the
     survey's NorSand-BA track (WP-144), not a regularizer.

**Inputs owed by TIMs:**
- the sand's grading (d50), e_max/e_min and the target D_r;
- any triaxial or plane-strain data;
- the calibration target for undrained cyclic CSR–N (it decides R1's CTXu outlier);
- the exact PDMY01 33° parameter set behind the 417.6 kPa control;
- any footing test they consider the reference.

Quote each control with its cone: TIMs' limit-point control is **PDMY01 at 33°**; the WP-138 comparison curve is
**UW DruckerPrager, ψ = 0, at 38°**.

---

## 10. T5 results — the campaign set against sand physics (flag OFF baseline)

`t5_element_physics.py` runs drained compression from isotropic states on the exact WP-134 oracle (`uw_model`
options, ε_a to 25 %). Outputs: `out_t5.md` and `out_t5_toyoura_e0.66.md`.

The contrast is DM04's own lab-calibrated Toyoura set (Verdugo & Ishihara data). It runs at e0 = 0.66, i.e.
D_r ≈ 0.83 with e_max 0.977 and e_min 0.597. That matches the D_r the campaign set's strength implies.

| set | test | p0 kPa | φ′_peak ° | ε_a at peak | peak dilatancy | reached critical state by 25 %? | Bolton (1986) check |
|---|---|---|---|---|---|---|---|
| campaign (e0 0.6944) | PS | 10 / 50 / 150 / 500 | 60.1 / 55.0 / 50.7 / 44.9 | 4.4 / 7.2 / 10.2 / 15.5 % | ψ_max 1.9 / 1.6 / 1.3 / 0.8° | **no** (ψ_end −0.09…−0.04) | stress–dilatancy missed by ~3–20× |
| campaign | TX | 10 / 50 / 150 / 500 | 48.6 / 45.8 / 43.1 / 39.1 | 4.1 / 7.3 / 10.7 / 16.6 % | (−dε_v/dε₁)max 0.067 / 0.058 / 0.048 / 0.031 | no | Δφ 15.6 / 12.8 / 10.1 / 6.1° vs 10·(−dε_v/dε₁)max ≤ 0.7°. Strength alone implies D_r 0.80–0.88 |
| Toyoura DM04 (e0 0.66) | PS | same | 51.1 / 49.2 / 47.0 / 43.2 | 1.1 / 2.0 / 3.0 / 4.8 % | ψ_max 23.4 / 21.5 / 19.0 / 14.0° | nearly (ψ_end ≈ −0.03) | Δφ = 0.5–0.7 × (0.8·ψ_max) |
| Toyoura DM04 | TX | same | 40.5 / 39.6 / 38.5 / 36.4 | 1.0 / 1.8 / 2.9 / 4.7 % | 1.10 / 1.00 / 0.88 / 0.64 | nearly | Δφ 9.3 vs 11.0° … 5.3 vs 6.4°: within ~17 % |

**Reading.**
- The campaign set's strength is that of a very dense sand. Its peak comes 3–4× too late in strain against Toyoura at the same D_r.
- It dilates ~15–25× less than such a sand must, by stress–dilatancy (Rowe; Bolton): the triaxial Δφ needs (−dε_v/dε₁)max ≈ 0.6–1.6; the model gives 0.03–0.07.
- A0 = 0.05 is 14× below DM04's Toyoura value. The high strength comes from nb = 3.5 through M_b = M·e^(−n_b·ψ),
  almost decoupled from volume change.
- Three consequences for the footing:
  1. The flow rule is extremely non-associated (friction ~45–60° against dilation ~1–2°). That is the §2 loss of
     ellipticity at positive hardening: the bands are partly a product of the calibration.
  2. A late peak (4–16 %) means the footing needs large settlements to mobilize. That is consistent with "no plateau
     to s/B 0.05" in every WP-138 leg.
  3. None of the three footing curves (SANISAND 967 kPa, DP 38° 824 kPa, PDMY01 33° 418 kPa at their own s/B) can be
     called physical until TIMs' lab data pin φ′_peak, ε at peak and dilatancy for their sand.
- **For the decision procedure:**
  - T6 (the footing benchmarks) is only meaningful after a calibration check.
  - Regularization cannot fix a constitutive mismatch of this size.
  - If TIMs confirm the set is what their data say, the physics check moves to their data. If not, recalibration comes
    first (or the NorSand-BA track, WP-144).

---

## References

- Dafalias, Y. F. & Manzari, M. T. (2004). Simple plasticity sand model accounting for fabric change effects. *J. Eng.
  Mech.* 130(6), 622–634. The h = b0/((α−α_in):n) rule.
- Rudnicki, J. W. & Rice, J. R. (1975). Conditions for the localization of deformation in pressure-sensitive dilatant
  materials. *JMPS* 23, 371–394.
- Sabet, S. A. & de Borst, R. (2019). Structural softening, mesh dependence, and regularisation in non-associated
  plastic flow. *IJNAMG* 43(13), 2170–2183. doi:10.1002/nag.2973.
- Hageman, T., Sabet, S. A. & de Borst, R. (2021). Convergence in non-associated plasticity and fracture propagation
  for standard, rate-dependent, and Cosserat continua. *IJNME* 122, 777–795. doi:10.1002/nme.6561.
- de Borst, R. & Duretz, T. (2020). On viscoplastic regularisation of strain-softening rocks and soils. *IJNAMG* 44(6),
  890–903. doi:10.1002/nag.3046.
- Gao, Z., Li, X. & Lu, D. (2022). Nonlocal regularization of an anisotropic critical state model for sand.
  *Acta Geotech.* 17, 427–439. doi:10.1007/s11440-021-01236-3.
- Mallikarachchi, H. & Soga, K. (2020). Post-localisation analysis of drained and undrained dense sand with a nonlocal
  critical state model. *Comput. Geotech.* 124, 103572.
- Galavi, V. & Schweiger, H. F. (2010). Nonlocal multilaminate model for strain softening analysis. *Int. J. Geomech.*
  10(1), 30–44.
- Liu, H. Y., Abell, J. A., Diambra, A. & Pisanò, F. (2019). Modelling the cyclic ratcheting of sands through
  memory-enhanced bounding surface plasticity. *Géotechnique* 69(9), 783–800. The memory-surface alternative to α_in.
- de Borst, R. & Mühlhaus, H.-B. (1992). Gradient-dependent plasticity: formulation and algorithmic aspects. *IJNME*
  35, 521–539.
- Vardoulakis, I. & Aifantis, E. C. (1991). A gradient flow theory of plasticity for granular materials. *Acta Mech.*
  87, 197–217.
- Mühlhaus, H.-B. & Vardoulakis, I. (1987). The thickness of shear bands in granular materials. *Géotechnique* 37(3),
  271–283.
- Desrues, J. & Viggiani, G. (2004). Strain localization in sand: an overview of the experimental results obtained in
  Grenoble using stereophotogrammetry. *IJNAMG* 28, 279–321.
- Pietruszczak, S. & Mróz, Z. (1981). Finite element analysis of deformation of strain-softening materials. *IJNME* 17,
  327–334.
- Siddiquee, M. S. A., Tanaka, T., Tatsuoka, F., Tani, K. & Morimoto, T. (1999). FEM simulation of scale effect in
  bearing capacity of strip footing on sand. *Soils Found.* 39(4), 91–109.
- Bolton, M. D. (1986). The strength and dilatancy of sands. *Géotechnique* 36(1), 65–78.
- Vesić, A. S. (1973). Analysis of ultimate loads of shallow foundations. *JSMFD* 99(SM1), 45–73.
- Perkins, S. W. & Madson, C. R. (2000). Bearing capacity of shallow foundations on sand: a relative density
  approach. *JGGE* 126(6), 521–530.
- Loukidis, D. & Salgado, R. (2011). Effect of relative density and stress level on the bearing capacity of footings
  on sand. *Géotechnique* 61(2), 107–119.
- Lau, C. K. & Bolton, M. D. (2011). The bearing capacity of footings on granular soils. I: Numerical analysis; II:
  Experimental evidence. *Géotechnique* 61(8), 627–638 and 639–650.
- Tatsuoka, F., Okahara, M., Tanaka, T., Tani, K., Morimoto, T. & Siddiquee, M. S. A. (1991). Progressive failure
  and particle size effect in bearing capacity of a footing on sand. *ASCE GSP* 27, 788–802.
- Wang, W. M., Sluys, L. J. & de Borst, R. (1997). Viscoplasticity for instabilities due to strain softening and
  strain-rate softening. *IJNME* 40, 3839–3864. The consistency-viscoplasticity alternative to Perzyna.
