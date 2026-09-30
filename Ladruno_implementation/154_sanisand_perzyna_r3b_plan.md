---
title: "WP-154 — R3b: Perzyna-type viscoplastic regularization inside LadrunoSANISAND (plan, no code)"
project: Ladruno
type: work-package plan
status: "PLANNED — plan only, no C++ (draft #895). Opened on the 2026-09-29 acoustic census (memo case C). Owner decides τ and adopt / not adopt (§8)."
owner: nmora
related:
  - "[[150_sanisand_regularization_memo]] (PR #892: §2 acoustic tensor, §4 R3 revised 2026-09-29, §9.1 decision procedure)"
  - "[[90_ladruno_viscoplastic_regularization_adr]] (ADR-90: §1.2(iii), §4.2 revert path, §4.5 V4, §7.3 w₂, 33022 reserved)"
  - "[[151_sanisand_reseat_singularity]] (R1, #893, merged)"
  - "[[152_sanisand_tension_cutoff]] (PR #894: the low-confinement separation)"
  - "[[134_sanisand_reference_integrator]] (the oracle)"
  - "[[LadrunoSANISAND_implex_guide]] §13"
tags: [plan, sanisand, sas-me, perzyna, viscoplastic, regularization, ellipticity, tims, wp-154, r3b]
updated: 2026-09-29
---

# WP-154 — R3b: Perzyna viscoplastic regularization inside `LadrunoSANISAND`

> [!summary] The short version
> - **Why now.** On the DM04 Toyoura footing the B/16−B/8 load gap opens at s/B 0.0335, about 5× before the
>   peak (0.167–0.180), exactly where the zone's non-elliptic share passes 1–5 %. That is WP-150 **case C**: a
>   pre-peak, non-associated loss of ellipticity. It is the one case R3b answers and R3a cannot.
> - **What.** An opt-in rate term in SAS-ME (`IntScheme 129`): λ̇ = ⟨f⟩/η, η = 2G·τ, f = SANISAND's own yield
>   function. Per substep the multiplier is Δλ = ⟨f₀ + N⟩ / (H + 2G·τ/δt). Flag `-sasPerzyna τ`, default OFF and
>   then byte-identical. No class tag. R1 and the #894 separation stay ON.
> - **How much.** The census puts the needed viscous hardening at H_v/2G = τ/Δt ≳ 0.015, not the memo's
>   0.17–0.44. At the deck's present ds_max of 4e-4 m, τ ≈ 1.2e-5 – 2.4e-5 m keeps the algorithmic tangent elliptic
>   with **no extra step cap**. So the price is mainly the rate bias, and that is a gate.
> - **Order.** Oracle first (WP-134 + Perzyna), then the census "need vs have" on the saved checkpoints, then
>   C++, then Esmeralda. τ is **re-derived** after any constitutive change (the G_max decay and the post-peak
>   dilatancy are out of scope here, §9).

---

## 1. Why

### 1.1 The evidence: memo case C

**Source.** Acoustic census of 2026-09-29, build `8ebde5cbd`, R1 full set (`-sasHFloor 1 -sasReseatHyst 1
-sasSoftCap 0.5`) plus the pre-review cutoff `-sasTensionCutoff 0.5 1.0`.
- Deck: `deck_toyoura`. Runs: `/mnt/deadmanschest/nmorabowen/ladruno_wp152/deck_toyoura/runs/W_TYR_{b8,b16,b8shear15,fig9_856}`.
- Driver: the WP-150 census tools (`census.py`, `extra.py`, `bands.py`), extended from
  `Ladruno_files/testbed/hypo_bearing/wp150_regularization/acoustic_vec.py` (#892).
- Output: Esmeralda `~/ladruno_wp152/analysis/acoustic_census/out/` (`census_report.md`, `extra.txt`,
  `band_paths.txt`).

The material is DM04 Table 1 (G0 125, M 1.25, c 0.712, m 0.01, h0 7.05, nd 3.5, …), e_init 0.6426 (D_r 0.88),
p_at 100, TanType 0. Every number in this section comes from that census unless another source is named.

| leg | B | mesh | peak q (kPa) at s/B | finished |
|---|---|---|---|---|
| b8 | 1.2 m | B/8 | 2 509 at 0.1798 | yes, to 0.200 |
| b8shear15 | 1.2 m | B/8, 15° sheared | 2 729 at 0.1752 | yes, to 0.200 |
| fig9_856 (Kimura Fig. 9) | 0.9 m | B/8 | 2 015 at 0.1672 | yes, to 0.1995 |
| b16 | 1.2 m | B/16 | (not reached; last s/B 0.0656, q 1 196) | no |

**Timing (case C).**
- The persistent B/16−B/8 gap passes −1 % at s/B 0.0010, −2 % at 0.0335 and −5 % at 0.0395.
- The zone's non-elliptic share (|x| ≤ 1.5B, y ≥ −1.5B) passes 1 % and 5 % at the same place: b16 0.0332 / 0.0336,
  b8 0.0357 / 0.0474, fig9 0.0337 / 0.0458.
- The peaks are at s/B 0.167–0.180, so the gap opens **about 5× before the peak.**
- In the memo's §4 R3 table this is **case C** ("before the peak: the hardening-regime, non-associated onset"). Its
  answer is **R3b**.

**Band orientation (the physics criterion of memo §4 R2 fails).**
- Before the peak the band is one element wide (FWHM/h 1.0–1.2 at s/B 0.035–0.05) and runs down a mesh column
  (b8: x = 0.632 m at every depth, inclination +0.0°).
- On the 15°-sheared mesh it follows the **same column, displaced with the mesh**: RMS(shear15 − b8 column moved
  with the mesh) = **0.001 m**, against 0.039 m for the physical path (h = 0.150 m; s/B 0.035 and 0.05).
- The material's own critical directions at the non-elliptic GPs are **45–65° from vertical** (median 44.4° at
  s/B 0.035, 63.8° at 0.05, 51.8° at the peak).
- Mesh orientation alone moves the peak by **+8.8 %** (2 729 vs 2 509 kPa).

**Mechanism.**
- Ellipticity is lost within **Kp_crit/2G ≈ 0.003–0.015** of each point's own peak. That is the median at the
  non-elliptic GPs, every leg and checkpoint (`extra.txt`), with one two-GP outlier (0.080, b8 s/B 0.011).
- About ½ of the non-elliptic points are locally softening (b:n < 0), and ⅓–½ are still hardening. For example,
  b8 at s/B 0.05 has 0.33 % hardening vs 0.28 % softening; b8 at the peak, 4.43 % vs 9.33 %.
- Kp/2G at those points: p10 down to −0.012 (b8 at s/B 0.10).
- **The associated control (R → Q, same Kp) has 0 non-elliptic points at all three peaks.** The non-associated
  flow is the cause, as on the campaign set (memo §2.1).
- **The ADR-90 V4 viscous blend removes it at β ≤ 0.9:** 0 / 0 / 567 non-elliptic at β = 0.5 / 0.9 / 0.99 (b8
  peak; shear15 0 / 0 / 340; fig9 0 / 0 / 389).

**So the required viscous hardening is H_v/2G ≳ 0.015**, against the memo's 0.17–0.44 for the campaign set (memo
§2.1, WP-138 checkpoints). The Toyoura sand dilates, so its flow rule is much closer to associated. τ is therefore
~10–30× milder than the memo feared.

> [!warning] Pitfall — ν. The census depends strongly on Poisson's ratio.
> - The deck builds the material at ν* = 1/3 (the matdesc, for K0) and switches to **ν = 0.05 for the push**
>   (`footing_ab.py:763-767`, per `census.py`).
> - b8 peak: **13.77 % non-elliptic at ν 0.05, 36.24 % at ν 1/3** (associated 0 in both).
> - The WP-150 scripts hard-code the campaign set. **Every census, C8 check and τ derivation in this WP reads the
>   material from the leg's matdesc and uses the push ν** (§5 O2).

**Context that constrains the plan (relayed by the TIMs orchestrator, 2026-09-29).**
- DM04's pre-peak element response matches Tatsuoka et al. (1986) at 49 kPa.
- After the peak it softens too little and dilates too much.
- A G_max decay test (vanilla `ManzariDafaliasRO`) and a stiffness ladder are running.
- **The model this regularizer is tuned for may change.** τ is therefore derived by a procedure (§2.4), never
  fixed as a number, and that procedure is re-run after any constitutive change (§9).

### 1.2 Why R3b and not the alternatives

| option | verdict | reason (evidence) |
|---|---|---|
| **none (disclose)** | **no** | Memo case A needs the PEAK to converge within ½ the scatter. Here the gap opens 5× before the peak, and orientation alone moves the peak +8.8 %. ADR-90 §1.2(iii) saw q(s/B) contract at τ = 0 on the 3-D deck at s/B ≤ 0.01. This deck contracts only until the bands form (memo §13: ratio 0.23–0.40, then 1.48). |
| **R3a, nonlocal void ratio** (Gao et al. 2022) | **no** | It regularizes ψ-softening: it changes how H evolves in the band, not the current tangent. ⅓–½ of the non-elliptic points are still **hardening** (b:n > 0), at fixed ψ. A nonlocal e cannot reach them (memo §3, §4 R3a). It also needs neighbour plumbing the fork does not have. |
| **Duvaut–Lions** (the ADR-90 wrapper `LadrunoOverstress`, 33022, or WP-F in-model) | **no** | It needs the inviscid σ̄ for the same Δε, so it **inherits every refusal** of the inviscid integrator (memo §3). The generic wrapper is also inexact for Cₑ(p) on every path (ADR-90 §4.5 V1), dissipation goes negative on unloading (V2), and β changes per Newton iterate (§3(3)). |
| **Cosserat / gradient-in-λ** | **no** | These are the only routes to an objective **width**, but no width deliverable is named. Memo §9.1 step 5 (D-d): 20·d50 ≈ 4–10 mm against h = 94–188 mm, so any width at this scale is numerical. New element families. |
| **R3b, Perzyna inside `LadrunoSANISAND`** | **yes** | It acts on the current tangent at any Kp, pre- or post-peak. It never computes an inviscid solution, so it cannot inherit a refusal. It uses the real Cₑ(p) (fixes V1), and η ∝ G keeps the Deborah number uniform. It is local, with no plumbing. The census gives a small H_v. |

R3b does **not** give an intrinsic band width (quasi-statics: ADR-90 A0, V3). Its deliverable is the ADR-90 §3.2
claim, re-based on the peak: *q(s/B) and q_u converge in h and in mesh orientation at a declared (τ, Δt), and their
τ-dependence is measured and reported.*

---

## 2. Design

### 2.1 Rate form

- Yield function: f = ‖s − p·α‖ − √(2/3)·m·p, i.e. SANISAND's own cone, the one SAS-ME already evaluates. Overstress
  is f > 0.
- λ̇ = ⟨f⟩ / η, with **η = 2G(p, e)·τ**. G is the model's own pressure-dependent modulus at the current state, so
  η ∝ √p.
  - The dimensionless overstress f/2G = τ·λ̇ is the same function of the local rate everywhere, which gives a
    uniform Deborah number.
  - At p → 0, η → 0 and the absolute overstress → 0. There is no apparent cohesion at the free surface (memo §6).
    A constant η would behave like `-Presidual` 20 kPa (+30 % q).
- Flow: unchanged. dσ = Cₑ:(dε − dλ·R), dα = dλ·(2/3)h·b, and dz as in DM04. R1's h (floor, hysteresis, soft cap)
  is used as is.
- Everything state-dependent is evaluated at the stage's state (SAS-ME U9).

### 2.2 Discrete form in SAS-ME: denominator H + 2G·τ/δt

Each SAS-ME substep covers pseudo-time δt = dT·Δt, where dT is the substep fraction and Δt the step (§2.3). Each of
its two stages takes the plastic multiplier from **linearized backward Euler in λ**:

```
    N   = df/dσ : Cₑ : Δε_sub             (the existing loading indicator, true yield gradient)
    f₀  = f(stage start)                  (≥ 0 with τ > 0: the stage may start outside the cone)
    H   = Kp + df/dσ : Cₑ : R             (the existing loading denominator; R1's h and cap inside)
    Δλ  = ⟨f₀ + N⟩ / (H + 2G·τ/δt)        (τ = 0:  Δλ = N/H on the surface, SAS-ME as today)
```

- **The τ → 0 limit is today's code.** The early-out keeps it byte-identical (§2.7).
- **The stiff limit is safe.** As δt/τ → 0, Δλ → 0 (elastic). As δt/τ → ∞, Δλ → (f₀ + N)/H, the inviscid return
  including its drift. The per-stage multiplier is L-stable, so no substep is forced by stiffness. SAS-ME's error
  control (σ, α, z; TolR) sizes the substeps as it does now.
- **Stage classification changes.** A stage is viscoplastic when f₀ + N > 0, not only when N > 0: an overstressed
  point relaxes even under a slightly unloading increment. It is elastic when f₀ + N ≤ 0.
- **Refusals.**
  - H ≤ 0 with f₀ + N > 0 is still refused (`loadingNonPosH`) whenever **H + 2G·τ/δt ≤ 0**, and counted separately
    when the viscous term is what kept it positive. With R1's cap, H ≥ κX > 0 already.
  - **RC_START_F** ("start outside yield") is not applied when τ > 0, because a committed f > 0 is the model's
    state. A finiteness check replaces it, and the maximum committed overstress is recorded (§2.8).
  - **The drift correction is not applied when τ > 0.** f > 0 is legitimate, and the relation f = η·Δλ/δt is carried
    by the multiplier itself. The α-bound check (ρ_α ≤ 1 + κ) and every other SAS-ME refusal are unchanged.
- **Tangent: one algorithmic tangent at the end state,**
  D_vp = Cₑ − (Cₑ:R)⊗(Q:Cₑ) / (H + 2G·τ/Δt), with the **step's** Δt. It applies when the last substep was
  viscoplastic, and it is non-symmetric. This is the one-step backward-Euler tangent, which is what C8 is written on.
  - The deck runs TanType 0 (elastic Cₑ: modified Newton), so Newton never sees D_vp. The regularization lives in the
    residual, and C8 is a property of the discrete map, not of the iteration matrix.
  - With TanType 1/2 an unsymmetric solver is mandatory (ADR-90 D11).
- **Ellipticity condition.** D_vp equals the inviscid continuum tangent with Kp replaced by Kp + H_v, where
  **H_v/2G = τ/Δt**. Point by point it is elliptic iff τ/Δt ≥ (Kp_crit − Kp)/2G. From §1.1 the median need is
  0.003–0.015 and the softening tail adds up to 0.012, so **τ/Δt ≥ 0.03 is the provisional design margin**. O2 (§5)
  replaces it with the measured per-GP distribution. The memo's "τ ≳ 0.2–0.45·Δt" was the campaign set's number.

### 2.3 Pseudo-time: τ is in metres of footing settlement

- The push is lane (b) of ADR-90 §3.1: a unit `sp` on the reference node under the 1-argument
  `LoadControl(−ds)`, time = u_y (WP-138 `footing_ab.py`, push setup). So **ops_Dt = −ds**: it is **negative** on
  every push step.
- ADR-90's convention "Δt ≤ 0 ⇒ inviscid" would therefore **switch R3b off silently on the campaign deck.** The
  same trap already bit IMPL-EX (ADR-92 red/blue B1/B2, `LadrunoSANISAND.cpp:3294-3312`), whose fix stores
  **|ops_Dt|**.
- **R3b uses Δt := |ops_Dt|.** τ is in the units of |pseudo-time|, i.e. **metres of footing settlement s** on this
  deck.
- It is reported as **τ̂ = τ/B** (dimensionless), so B 0.9 m and B 1.2 m run at the same τ̂, and as the Deborah
  number De = τ/s_peak.
  - Example: τ = 2.4e-5 m at B 1.2 m gives τ̂ = 2e-5 and De ≈ 1.1e-4 at s_peak ≈ 0.21 m.
- **No latch.** ADR-90 D4 asked for β and Δt to be latched at `newStep`. Under lane (b) Δt is constant within a
  step. A latch keyed on "the first evaluation after a commit or revert" would latch the revert's Δt = 0 (below).
  Instead:
  - Δt is read at every evaluation;
  - a change of |ops_Dt| between two trial evaluations of one step is **counted** (`vpDtChangedInStep`), and must
    be 0 on the campaign lane;
  - `DisplacementControl` and `ArcLength` are not admissible lanes (ADR-90 R11). A material cannot detect the
    integrator, so the counter and the guide carry this.
- **Non-uniform Δt.** The deck's adaptive ds ladder varies Δt 1.2e-5 – 8e-4 m per step (step CSVs of the census
  legs). A smaller ds raises τ/Δt, so ellipticity is safe, but the effective regularization varies. The committed
  |Δt| min/max is reported per leg (ADR-90 M8).

### 2.4 How τ is chosen, and re-chosen

A procedure, not a number (ADR-90 §3.3's honest-framing test: τ is **never** tuned to a target q or width):
1. Fix the leg family's **ds_max** (today 4e-4 m after s/B 0.07 at B 1.2 m, and 8e-4 m at s/B < 0.035).
2. From the O2 census (need per GP), take the **smallest τ̂** such that **τ/2** keeps C8 at ds_max over the
   checkpoints to the peak. The sweep runs {τ/2, τ, 2τ}, so every leg of it must be elliptic.
3. Provisional value from §1.1: τ/2 ≥ 0.03·ds_max. With ds_max = 4e-4·(B/1.2) m this gives **τ̂ = 2e-5**
   (τ = 2.4e-5 m at B 1.2, 1.8e-5 m at B 0.9), and **no step cap below the deck's own**.
4. **Re-run 1–3 after any constitutive change** (G_max decay, stiffness ladder, post-peak dilatancy; §9). Kp_crit
   and Kp both move with the model.

### 2.5 Composition with R1 and with #894

- **R1 (floor, hysteresis, soft cap) stays ON.** R3b cannot lift the α_in singular set: at a re-seat h → ∞ so
  Δλ → 0 (harmless), but where b:n < 0 with Kp → −∞, H + 2G·τ/δt can still go ≤ 0 without the cap (memo §2.3,
  §3).
  - The parser **refuses** `-sasPerzyna` outside `IntScheme 129`.
  - SAS-ME under `-implex` is already refused (`LadrunoSANISAND.cpp:~5003`).
  - It **warns** when the R1 set is incomplete. The ablation stays possible; the guide names the full set.
  - The soft cap acts on the inviscid H. The viscous term is added after it.
- **#894 separation (`-sasTensionCutoff`).**
  - A SEPARATED point (σ = p_min·I, no shear) has no cone. Perzyna is inactive there.
  - At re-contact the point restarts on or inside the cone with zero overstress.
  - ENTRY masks only low-p codes (3/4/6/9). Since η → 0 as p → 0, R3b does not change the low-p regime by
    construction. It can still shift **when** E2 fires, because it changes the substep counts. The separation
    counters are compared with τ = 0 on the same leg (gate G7).
  - Wire and census slots are appended **after** #894's. If #894 has not merged when the C++ starts, WP-154 rebases
    on it.

### 2.6 Flag, parser and options

- `-sasPerzyna τ` (τ ≥ 0, in pseudo-time units; 0 = OFF, the DEFAULT). It is refused with τ < 0, with a non-finite
  τ, and outside `IntScheme 129`.
- It lives in `LadrunoSasOptions` (a Ladruno block in vanilla `ManzariDafalias.h`, so a vanilla-ledger row is updated
  at implementation) as `double vpTau`.
- `sasOptions` (33101) goes 12 → 13 values (after #894's 12).
- Optional: `setParameter sasPerzynaTau` for sweeps without rebuilding (ADR-90 D6). It is not required, because the
  deck builds one material per leg.
- `ladrunoSANISANDReplay` gains a `-dt` argument: the replay loop has no domain clock, so without it the element
  gates would run inviscid.

### 2.7 Off-switch, getCopy, wire, revert

- **τ = 0 is an early-out** around every new line, so the code path is byte-identical (ADR-90 §4.1/D9; `0·NaN`).
  Gate: WP-151's 643-replay baseline plus WP-152's paths, bit-identical.
- **No new state variable.** The overstress is carried by σ itself (f > 0 at commit). The committed state stays
  (σ, α, α_in, z, e). The memo's "+2 doubles" (committed overstress, latched Δt) are not needed.
- **getCopy.** τ travels in the `LadrunoSasOptions` struct copy, the same path as R1/#894 options. Test: a
  `getCopy("PlaneStrain")` point runs the same trajectory as its prototype with `sasOptions` equal.
- **Wire.** `LWIRE_SAS_OPT_N` +1 (τ). The appended sasStats columns ride the existing census block, and the layout
  tag (`kLadrunoSanWireTag`) is bumped. Test: a sendSelf/recvSelf round trip value-checks `sasOptions` and `sasStats`.
- **Revert (ADR-90 §4.2).** `Domain::revertToLastCommit()` sets dT = 0 and re-applies the load, so the first
  evaluation after every cutback runs at Δt = 0.
  - With τ > 0 and Δt = 0 the update takes the **inviscid early-out** (the fork convention, ADR-90 §4.2)
    and is counted as `vpDt0Evals`.
  - The retried step's `newStep` then sets Δt = −ds/2 and SAS-ME restarts from the committed state, so nothing
    leaks.
  - **`vpDt0Commits`, the number of commits whose last trial ran with τ > 0 but Δt = 0, must be 0.** A non-zero
    count means regularized and unregularized steps were mixed. It fails the leg.
- No statics, so it stays thread-safe for the WP-131/146 threaded loop.

### 2.8 Census additions

**`sasStats`,** appended after #894's columns (earlier indices unchanged):

| column | meaning |
|---|---|
| `vpStages` | stages that took the viscoplastic branch |
| `vpRelaxStages` | … of which started with f₀ > TolF (relaxation of an existing overstress) |
| `vpMaxFRel` | max committed f / (√(2/3)·m·p): overstress in cone radii (m = 0.01 makes the cone 0.8 % of p wide, §7 R-4) |
| `vpLastF` | last committed f (kPa) |
| `vpNonPosHHeld` | stages where H ≤ 0 but H + 2G·τ/δt > 0 |
| `vpDt0Evals` / `vpDt0Commits` | §2.7; the second must be 0 |
| `vpDtMin` / `vpDtMax` | committed |Δt| range |
| `vpDtChangedInStep` | §2.3; must be 0 on lane (b) |

**Response `tangentVP`** (proposed id **33102**, the next in the band after 33101; a response id, not a class tag;
confirm it is free at implementation):
- D_vp at the committed state with the last committed Δt, whatever the TanType.
- `tangentEP` (33099) keeps its meaning, the inviscid continuum tangent.
- Together they give the census **"need" (tangentEP) and "have" (tangentVP)** at every checkpoint.
- The per-point acoustic scan (721 directions) stays in the post-processor. It is too costly per commit.

### 2.9 Class tag: none

- R3b is a **model option of `LadrunoSANISAND`'s SAS-ME**, exactly like R1 (#893) and the separation (#894). There
  is no new class, no new wire class and no `classTags.h` edit.
- **33022 stays RESERVED and unused.** It belongs to ADR-90's generic wrapper (`LadrunoOverstress`), which this WP
  does not build and which the memo proposes to retire.
- Taking 33022 would imply a separate `NDMaterial` class. That is the wrapper architecture rejected in §1.2.

---

## 3. What is claimed, and what is not

- **Claimed, if the gates pass:** at a declared τ̂ (and ds_max), the Toyoura footing's q(s/B) and peak q_u contract
  under h-refinement and under mesh orientation. The band leaves the mesh columns. The rate bias is inside ½ of the
  test scatter.
- **Not claimed:** a band width, a τ that is a soil property (ADR-90 S2), or anything about the post-peak physics
  (softening too little, dilating too much: constitutive, §9). **Not claimed either** for the campaign set without
  its own census, since its need is 10–30× larger (memo §2.1).

---

## 4. Gates (from memo §4 R3b, adapted)

| id | gate | measure | pass |
|---|---|---|---|
| **G1 = C8** | the algorithmic acoustic tensor is elliptic at every committed state | min normalised det of n·D_vp·n over 721 directions, every census checkpoint to the peak (s/B 0.005 … peak), `tangentVP`, push ν | 0 non-elliptic GPs at τ/2, τ and 2τ; report the min det and its step |
| **G2** | q(s/B) at fixed τ contracts under h-refinement | B/8 → B/16 (→ B/32 to s/B 0.05 if feasible) at s/B {0.035, 0.05, 0.10, peak} | the gap ratio (B/16−B/32)/(B/8−B/16) < 1, and the B/8–B/16 gap at the peak inside ½ scatter |
| **G3** | … and under mesh orientation | b8 vs b8shear15 at fixed τ | the peak gap below ½ scatter (today +8.8 %) |
| **G4** | Δt-convergence at fixed τ | b8 at ds_max and ds_max/2 | the q(s/B) difference below the solver floor (0.8–1.4 %, ADR-90 §1.2(iii)) |
| **G5** | rate bias | q(τ)/q(τ → 0) at {τ/2, τ, 2τ}, peak and s/B 0.05/0.10; q(τ → 0) from a linear extrapolation in τ (the τ = 0 leg is mesh-dependent, so it is not the reference) | bias at τ inside **½ the test scatter: ±4 %** (from Kimura's N_γ scatter; see §10 Q1) |
| **G6** | the band leaves the mesh | the census band-path tool (`bands.py` / `r2_analysis.py`): path, inclination, FWHM/h, w₂ (ADR-90 §7.3), b8 vs shear15 | RMS(shear15 − b8, physical) < RMS(shear15 − b8, column moved with the mesh) (today 0.039 vs **0.001** m), and the inclination leaves 0 ± 2° toward the critical 45–65° (the proposed threshold, §10 Q5) |
| **G7** | composition | R1 counters (`hFloored`, `hSoftCapped`, `reseatHeld`) and #894 separation entries, τ vs τ = 0 on the same leg | 0 `loadingNonPosH`; separation entries reported, and no new refusal class |
| **G8** | element tests unchanged at rate → 0 | Gate 0 (DM04 Table 1: the 17 Verdugo & Ishihara triaxials) and T5 at e 0.643 (PS/TX, p′ 10–500 kPa), with τ at the footing's De mapped to each test's pseudo-time; Tatsuoka et al. (1986) at 49 kPa | max \|Δq\|/q_max ≤ 2e-3 vs τ = 0 (the C++-vs-oracle scale, 1.9e-3 in memo §12); the Tatsuoka match unchanged |
| **G9** | provenance | per leg: τ, τ̂, De, ds_max, the \|Δt\| min/max, `vpDt0Commits`, `vpDtChangedInStep`, build, `ladrunoBuild()` | `vpDt0Commits` = `vpDtChangedInStep` = 0 |
| **G0** | the off-switch | τ = 0: the WP-151 643-replay baseline plus the WP-152 paths | bit-identical |

**Failure branch (ADR-90 §8 P2, kept):** if G2/G3 fail at the C8-derived τ, the WP does **not** iterate on τ. It
reports, and the owner chooses between a smaller ds_max (cost) at the same τ, or disclosure.

---

## 5. Oracle first, then C++

| phase | content | exit |
|---|---|---|
| **O1 — oracle** | Extend the WP-134 oracle (`Ladruno_scripts/sanisand_reference/`, as a testbed copy with toggles, the WP-151 pattern `Ladruno_files/testbed/sanisand_reseat_r1/sanisand_r1/`) with `viscous=(tau, dt)`. The increment's t ∈ [0, 1] maps to Δt. A new mode has λ̇ = Δt·⟨f⟩/(2Gτ) (no consistency condition); events f → 0 in both directions; R1 options on. Radau is already stiff-safe. | (i) τ → 0 reproduces the inviscid oracle (rtol 1e-8); (ii) τ → ∞ reproduces the elastic path; (iii) the analytic anchor: at critical state (Kp = 0, D = 0) the steady overstress is f = 2G·τ·λ̇ on a drained PS path |
| **O2 — need vs have, no C++** | Run the census on the saved checkpoints with the matdesc and push ν: per GP, (Kp_crit − Kp)/2G, i.e. the **required τ/Δt**, and its distribution to the peak; then C8 on D_vp at the §2.4 τ/2 and ds_max. | the τ̂ of §2.4 fixed by measurement (replaces the provisional 0.03) |
| **O3 — oracle on the footing's states** | The WP-151 fan method (`cxx_fan.py` pattern): committed footing GPs × 32 directions × 2 magnitudes, at τ/Δt ∈ {τ/2, τ, 2τ}/ds_max. Also the Gate 0 paths at 3 De. | the oracle's q bias on the element paths, before any footing run (a G8 preview) |
| **C1 — C++** | §2.2 in `LadrunoSANISANDSasME.cpp`; options, parser, `sasOptions`, `tangentVP`, sasStats, wire, the replay `-dt` (§2.6–2.8) | G0 bit-identity; the C++ matches the oracle **step by step** on O1/O3: benign states ≤ 2e-7 relative at TolR 1e-7 (the WP-129 standard), element paths ≤ 2e-3·q_max; `tangentVP` vs central FD of the step map reported (a one-step tangent over a substepped update is not exact); getCopy, wire and the two-instance interleave |
| **C2 — tests** | `tests/test_ladruno_sanisand_perzyna.py` (fast tier): off-switch, oracle match, revert (`vpDt0Commits`), negative-Δt deck (the §2.3 trap as a regression), `-implex` and non-129 refusals, the separation interplay | green on Zone-A; a mutation gate (drop the viscous term ⇒ red) |
| **E — Esmeralda** | §6 | G1–G7, G9 |

---

## 6. Footing test matrix (Esmeralda)

The deck is `deck_toyoura` (B 1.2 m, γ 15.9, K0 0.5, e_init 0.6426, DM04 Table 1, push ν 0.05). Every leg has R1 full
set + `-sasTensionCutoff` (#894's reviewed values), ds_max 4e-4·(B/1.2) m, ckpt every 5 steps. The Kimura case is
`fig9_856` (B 0.9 m, D_r 0.856, measured q ≈ 1 950 kPa at s/B ≈ 0.092; memo §14).

| # | B | mesh | τ̂ | purpose | gates |
|---|---|---|---|---|---|
| 1 | 1.2 | b8 | 0 (new build, flag OFF) | the build control against `W_TYR_b8` | G0 (on the curve) |
| 2–4 | 1.2 | b8 | τ̂/2, τ̂, 2τ̂ | bias; the C8 margin | G1, G5, G7 |
| 5–7 | 1.2 | b16 | τ̂/2, τ̂, 2τ̂ | h-contraction at each τ | G1, G2 |
| 8 | 1.2 | b8shear15 | τ̂ | orientation, band path | G3, G6 |
| 9 | 1.2 | b8 | τ̂, ds_max/2 | Δt-convergence | G4 |
| 10 | 1.2 | b32 | τ̂, to s/B 0.05 only | the pre-peak onset at a third mesh (if feasible: b16 took 20.7 h to s/B 0.066) | G2 (pre-peak) |
| 11–13 | 0.9 | b8 | τ̂/2, τ̂, 2τ̂ | the Kimura benchmark; bias | G5, T6 |
| 14 | 0.9 | b16 | τ̂ | h-contraction | G2 |
| 15 | 0.9 | b8shear15 | τ̂ | orientation | G3, G6 |

The τ = 0 references are the census legs themselves (`W_TYR_*`, build `8ebde5cbd`). Every leg writes the §2.8
columns and `tangentVP` at checkpoints for G1.

---

## 7. Effort and risks

**Effort, about 2–3 weeks:**
- O1–O3: 3–4 days (O2 needs no engine);
- C1–C2: 4–5 days;
- Esmeralda: 5–8 days wall (b16 legs dominate), plus the analysis.

| # | risk | mitigation |
|---|---|---|
| **R-1 cost** | The step constraint is τ/2 ≥ 0.03·ds_max (§2.4). At τ̂ = 2e-5 the deck's present ds_max already satisfies it. If O2 finds a larger need, either τ grows (bias) or ds_max shrinks (cost ∝ 1/ds_max, i.e. b8 ~8 h → ~16 h at ds/2). | O2 decides before any run; the owner picks between the two (§8) |
| **R-2 rate bias** | The overstress raises q. ADR-90 measured +9.7 / +12.5 / +17 % for the wrapper at De 3e-4 / 1e-3 / 3e-3 (R9). Here De ≈ 1e-4, but the in-band rate is ~B/h larger than the diffuse rate. | G5 with the τ → 0 extrapolation; O3 previews it on the element paths |
| **R-3 low-p interaction** | η → 0 at p → 0 by design, but R3b can shift when #894's E2 fires (substep counts), and so which surface points separate. | G7 compares entries τ vs 0; the separation's own gates are unchanged |
| **R-4 thin cone** | m = 0.01 makes the cone 0.8 % of p wide. An overstress of ~1 kPa at p = 100 kPa is already more than one cone radius, so n follows σ rather than α near the band. The model stays defined (b, d use α), but the peak-strength bias lives here. | `vpMaxFRel` per leg; O3 measures it at the element level |
| **R-5 non-uniform Δt** | The adaptive ladder halves ds, so τ/Δt varies per step. | report the \|Δt\| range; G4 |
| **R-6 ν** | The need scales strongly with ν (§1.1: 13.8 vs 36.2 % non-elliptic). A ν change re-opens τ. | §2.4 step 4 |
| **R-7 the model moves** | G_max decay and post-peak dilatancy work is running (§9). | τ by procedure; re-derive after any change |
| **R-8 revert / clock** | Δt = 0 after a cutback; negative ops_Dt on this deck. | §2.3, §2.7 counters, a regression test |

---

## 8. Owner decisions

1. **τ̂** — the O2-derived value, provisionally 2e-5 (τ = 2.4e-5 m at B 1.2), and whether the fork may trade a
   smaller ds_max (cost) against a larger τ (bias) if O2 finds a larger need.
2. **Adopt or not** — after the §6 matrix: adopt R3b as the TIMs setting (like R1's full set), or disclose (ADR-90
   lane (d)).
3. **Filing (memo §9.1 step 7, D-c)** — if adopted: the ADR-90 revision the memo proposes (retire the wrapper, WP-F
   becomes in-model Perzyna, 33022 reserved and unused). This plan does not edit ADR-90.

---

## 9. Out of scope, and the ordering

- **The constitutive fixes are not part of this WP:** the small-strain stiffness (G_max decay, `ManzariDafaliasRO`
  test, the stiffness ladder), and DM04's post-peak response (it softens too little and dilates too much).
- **Ordering:**
  1. the constitutive choice settles;
  2. then τ is (re-)derived by §2.4 on the settled model;
  3. then the §6 matrix is run or re-run.
- If a constitutive change lands after §6 has run, §2.4 and the G1/G5 legs are repeated before any footing number is
  quoted. O1–C2 carry over unchanged, since the regularizer is independent of the constants.

---

## 10. Left open by the census or ADR-90 — flagged, not guessed

- **Q1 — the scatter.** "±4 %, from Kimura's N_γ scatter" is the orchestrator's number. It is not yet in the fork's
  record (memo §9.1 step 6 still says "to be read from the sources"). Is ±4 % the scatter or half of it? G5 uses it
  as the half.
- **Q2 — the Δt = 0 branch.** The fork convention (ADR-90 §4.2) is **inviscid**, while Perzyna's own Δt → 0 limit is
  **elastic**. The plan follows the convention and makes the choice result-neutral through `vpDt0Commits` = 0.
  Owner to confirm.
- **Q3 — no latch.** This departs from ADR-90 D4, with counters instead (§2.3). Justified on lane (b); confirm.
- **Q4 — the margin.** The census tabulated Kp_crit/2G only at the non-elliptic GPs, and Kp only at p10. The per-GP
  need (Kp_crit − Kp)/2G was not tabulated, so 0.03 is provisional until O2.
- **Q5 — the band-orientation pass line.** The census gives the critical directions (45–65°) but no pass criterion.
  The G6 threshold is proposed.
- **Q6 — ν.** Is the push-stage ν = 0.05 itself part of the model to be settled (§9)? It moves the need ~2.6×.
- **Q7 — B/32** feasibility on Esmeralda (b16 has not reached its peak in 20.7 h).
- **Q8 — R1 enforcement.** The plan warns rather than refuses on an incomplete R1 set.
