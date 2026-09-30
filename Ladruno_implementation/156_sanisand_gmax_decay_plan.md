---
title: "WP-156 — small-strain stiffness (G_max decay, pressure exponent n) in LadrunoSANISAND (plan, no code)"
project: Ladruno
type: work-package plan
status: "PLANNED — plan only, no C++ (draft PR). Opened on the 2026-09-29/30 G_max-decay test (vanilla ManzariDafaliasRO) and the TIMs stiffness requirement. Owner and TIMs decide n, the G_max source and the stiffness tolerance (§8)."
owner: nmora
related:
  - "[[151_sanisand_reseat_singularity]] (R1, #893, merged)"
  - "[[152_sanisand_tension_cutoff]] (PR #894: the low-confinement separation)"
  - "[[154_sanisand_perzyna_r3b_plan]] (PR #895: R3b; its τ is re-derived after any constitutive change)"
  - "[[150_sanisand_regularization_memo]] (PR #892)"
  - "[[134_sanisand_reference_integrator]] (the oracle)"
  - "[[86_ladruno_sanisand_adr]] §4.1 (the RO shadow hazard: GetElasticModuli must stay non-virtual)"
  - "[[LadrunoSANISAND_implex_guide]] §13"
tags: [plan, sanisand, sas-me, small-strain, gmax, ramberg-osgood, stiffness, tims, wp-156]
updated: 2026-09-30
---

# WP-156 — small-strain stiffness in `LadrunoSANISAND`: G_max decay with a pressure exponent n

> [!summary] The short version
> - **Why.** TIMs need the footing's **initial stiffness** to calibrate the TIM macroelement, as well as the limit load.
>   DM04's footing start is 2–3× too soft and concave-up; Kimura's test is concave-down. Of every leg tried, only
>   a G_max decay (vanilla `ManzariDafaliasRO`) gives the right shape, but it reaches only 0.53–0.61 of the test.
>   Vanilla RO also cannot carry the curve to the peak: it has no SAS-ME, no R1 and no cutoff.
> - **What.** An opt-in elastic law inside `LadrunoSANISAND`'s SAS-ME (`IntScheme 129`). G = G_max(p, e)/T. The small-strain modulus is
>   G_max = B_g·p_at·F(e)·(p/p_at)^n, with n configurable. T is an RO/Masing decay on the stress-ratio distance
>   from a reference state (Papadimitriou & Bouckovalas 2002, as in vanilla RO). The decay floors at **DM04's own
>   calibrated G**, so beyond the decay range the model is DM04. b0 stays on G0. The reference state is committed, reverted,
>   copied and sent. The K0 seat is built into the stage flip. Reversal detection uses a relative hysteresis that
>   does not depend on p. Default OFF, and then byte-identical. No class tag.
> - **What it is not.** It does not add an amplitude-dependent plastic modulus. That is a conditional follow-up
>   (§3): the 1e-4 deficit is real at 49 kPa (0.556 against 0.721 on the proxy target), but only 0.440 against
>   0.474 at 5 kPa, where the footing's early soil sits.
> - **Order.** Oracle first (WP-134 plus the decay), with an analytical elastic check of the footing in parallel
>   (a strip on G(z) ∝ z^n, no C++ needed), then C++, then Esmeralda. The WP is independent of WP-154, because stiffness is
>   decided before localization. WP-154's τ is re-derived afterwards.

---

## 1. Why

### 1.1 The owner's framing (2026-09-30, relayed by the TIMs orchestrator; handoff 07 §13)

- **The TIMs objective is unchanged:** the footing's limit load (the plateau), with physical coherence.
- **In addition, TIMs need the footing's INITIAL stiffness** to calibrate the TIM macroelement: the initial
  stiffness, and the stiffness recursion with G/G0. So it is a **deliverable, not a limitation to disclose.**
- DM04's footing start is today **2–3× too soft and concave-up**; the test is **concave-down**.
- This WP must deliver the initial stiffness **without breaking** any of these:
  - the peak;
  - the element response that already matches the lab (Tatsuoka et al. 1986 at 49 kPa);
  - R1 (#893), the separation (#894) or R3b (WP-154).

### 1.2 The evidence

**How numbers are cited.** Every number below was produced on 2026-09-29/30 and names its source.
- "Workbench" means the TIMs Workbench repository (private), branch `work/ape/2d-model` at `9f040d4`, folder
  `Tries/2d-model/`.
- "Esmeralda" means `~/ladruno_wp152/` on the cluster.
- `bin_ro` is build **`8ebde5cbd`** plus only the Python registration of `ManzariDafaliasRO` (+2 lines in
  `OpenSeesNDMaterialCommands.cpp`). It was built on node4, 2026-09-29 23:31 (log `~/ladruno_wp152/build_ro.out`), and
  `G_VAN_ro` on `bin_ro` reproduces `G_VAN` on `bin` step for step (Workbench `model/gmax-ro-2026-09-30/README.md`).

**The footing.** Kimura case: B 0.9 m, e 0.635, D_r 0.90, γ 15.9, B/8, push ν 0.05.
- Deck: Esmeralda `deck_toyoura_ladder/`, driver `footing_ab.py` v3b (`driver_v3.diff` in the Workbench folder).
- Legs: `runs/{G_RO, G_RO_h0x15, G_VAN_ro}` (jobs 148797 / 148803 / 148798, on `bin_ro`), `G_VAN` (148791, `bin`), and
  the SAS-ME ladder legs `L_*` (build `8ebde5cbd`).
- Summary: handoff 07 §13 (Workbench `07 — Handoff, Where the SANISAND Wall Stands.md`, 2026-09-30).
- The test curve is Kimura et al. (1985) Fig. 9, V85.6, from `references/kimura1985_fig9/kimura1985_fig9_digitized.csv`
  (digitized 2026-09-29). The q values below were interpolated from it on 2026-09-30.
- The secant q/(s/B), in MPa per unit s/B, was computed here from the tabulated q.

| leg | q (kPa) at s/B 0.005 / 0.01 / 0.02 | secant (MPa) | shape |
|---|---|---|---|
| **Kimura V85.6 (test)** | 227 / 410 / 743 | 45.4 / 41.0 / 37.2 | **falls: concave-down** |
| DM04, SAS-ME base (vanilla `G_VAN` identical) | 74 / 160 / 349 | 14.8 / 16.0 / 17.5 | rises |
| SAS-ME, G0 ×3, b0 fixed | 183 / 389 / 835 | 36.6 / 38.9 / 41.8 | rises, crosses the test at s/B ≈ 0.012 |
| SAS-ME, G0 ×2, b0 fixed | 134 / 288 / 628 | 26.8 / 28.8 / 31.4 | rises, crosses at ≈ 0.03 |
| **`G_RO`: vanilla RO**, B 750 (G_max 128 MPa at p 100 kPa = 3.08× DM04's G), a1 0.49, γ1 5e-4, κ 2, h0 7.05 | **121 / 230 / 453** | **24.2 / 23.0 / 22.7** | **falls: the right shape** |
| `G_RO_h0x15`: RO, h0 ×1.5 | 130 / 246 / – | 26.0 / 24.6 / – | falls; +7 % |
| seating 16/40 kPa, K0 0.4, surcharge 5 kPa | within 10 % of the base | – | rises |

- **Only G_max decay gives the right shape.** A constant stiffness factor raises the start but keeps the curve
  concave-up, so it crosses the test.
- **The level is still 0.53 / 0.56 / 0.61 of the test** (G_RO against V85.6, at s/B 0.005 / 0.01 / 0.02).
- The test's first 0.002 of s/B is not usable:
  - The CSV header gives a reading uncertainty of ±0.06 mm model, i.e. **±0.002 s/B**, where q < 150 kPa. At s/B 0.002
    the test has q 95 kPa, secant 47.7, so that secant carries about ±100 %.
  - The curve also has a toe: q ≈ 5 kPa at s/B 0.0008.
  - The s/B 0.002 secant is therefore reported, never gated (§5 G5).

**The element** (single element, drained plane strain, e 0.635). Workbench `model/gmax-ro-2026-09-30/element/`:
`engine/*.json` on `bin_ro` via `src/el_test.py`; `oracle/*.json` from `src/ro_oracle.py` on the WP-134 oracle,
paper Options; the table below is from `src/el_table.py`, re-run 2026-09-30.
- The "target" is a **proxy** G/G_max = 1/(1 + (γ/γr)^0.92) with γr = 4e-4·(p/100)^0.5, a Darendeli-form curve. The
  README rates it "medium-low confidence". It is not a Toyoura measurement (§5 G2 replaces it).

| run (engine) | p′0 kPa | G/G_max at γ 1e-6 / 1e-5 / 1e-4 / 1e-3 | φ′_peak | ε_a at peak | E50 (MPa) |
|---|---|---|---|---|---|
| proxy target | 49 | 0.994 / 0.955 / **0.721** / 0.237 | – | – | – |
| `ro_49` (RO, iso start) | 49 | 0.999 / 0.920 / **0.556** / 0.240 | 51.3° | 1.80 % | 36.9 |
| `van_49` (DM04) | 49 | 0.324 / 0.324 / 0.299 / 0.201 | 51.4° | 1.80 % | 34.8 |
| `h15_ro_49` (RO, h0 ×1.5) | 49 | 0.999 / 0.933 / 0.603 / 0.273 | 51.7° | **1.43 %** | **44.7** |
| proxy target | 5 | 0.984 / 0.882 / **0.474** / 0.098 | – | – | – |
| `ro_5` (RO, iso start) | 5 | 0.997 / 0.826 / **0.440** / 0.190 | 56.1° | 0.77 % | 14.0 |

- **DM04 alone matches Tatsuoka et al. (1986) at σ3′ 49 kPa** (Workbench `references/tatsuoka1986_element/compare.md`,
  2026-09-29):
  - Fig. 16(a), e 0.716: E50 26.5 against 26.9 MPa (×0.99), ε_peak 2.14 % against 2.05 % (×1.04);
  - point tests: ε_peak 1.81–2.05 %.
- **RO at 1e-4 reaches only 0.556·G_max** against the proxy's 0.721, because DM04's plastic modulus takes over from
  about 1e-5 (handoff §13 finding 2). At 5 kPa the deficit is small: 0.440 against 0.474.
- **h0 is the wrong lever.** h0 ×1.5 adds only 7 % at the footing start, and it over-stiffens the element:
  - E50 is 44.7 MPa, 1.66× Tatsuoka's 26.9 (the gate was 1.3×);
  - the peak strain is 1.43 % against the lab's 1.81–2.05 %.

  **REJECTED as a lever.** Plain RO keeps ε_peak at 1.80 % (as DM04) and raises E50 by 6 % over DM04 at the same e
  (36.9 against 34.8 MPa).
  - Handoff §13 quotes that E50 as 1.37× Tatsuoka's 26.9 MPa. That test was at e 0.716, not 0.635, so the ratio mixes
    densities. G3 (§5) compares at equal e.

**The contradiction (handoff §13 finding 4).** At the element level G_RO is already stiffer than the lab (whose external
strains under-read), yet the footing is softer than the test (~0.57×). The candidates, in the handoff's order:
- **(a) the pressure exponent.**
  - Toyoura's G_max exponent is about 0.4 (Iwasaki, Tatsuoka & Takagi 1978, S&F 18(1):39–56, at γ ≈ 1e-6).
  - RO and DM04 hard-code 0.5. At a fixed G_max at p_at, n 0.4 gives **+26 % at p′ 10 kPa** and +58 % at 1 kPa.
    Most of the soil that feels the early load sits at p′ 1–20 kPa.
- **(b) fabric**: the V-bedding is Tatsuoka's stiffest direction (δ 90°). DM04 is fabric-free.
- **(c) the test conditions**:
  - the radial g-field at 30g;
  - the unstated seating and roughness;
  - the footing's own weight, about 16 kPa at 30g (red-team A5, `references/kimura1985_fig9/redteam_report.md`);
  - side-wall friction.
- **(d) a small-strain plastic modulus that depends on strain amplitude** (§3).

This WP addresses (a) and delivers the decay. It separates the boundary-value part from the constitutive part with an
analytical elastic check (§5 G4). (b) and (c) are out of scope (§9). (d) is a conditional follow-up (§3).

### 1.3 Why a model option inside `LadrunoSANISAND`, and not vanilla RO

- **Vanilla RO has no SAS-ME, no R1 and no cutoff.** Its legs run on ModifiedEuler and stop at the free surface
  (the WP-152 wall), so they cannot carry the curve to the peak. The peak is the TIMs objective, and R1 plus the cutoff
  are what reach it.
- **A `ManzariDafaliasRO`-style subclass of `LadrunoSANISAND` would be silently INERT under SAS-ME** (source reading,
  2026-09-30):
  - RO's decay lives in two `GetElasticModuli` overloads that **shadow** the base's (`ManzariDafaliasRO.h:93-97`).
  - Inside base member functions those calls bind statically to the **base** body.
  - Vanilla ModifiedEuler (`ManzariDafalias.cpp:1665–2175`) never calls `GetElasticModuli`. It uses the `mG`/`mK`
    that RO's `integrate()`/`commitState()` (`ManzariDafaliasRO.cpp:161-203`) computed at the committed state. That is
    the only way the decay ever reached the G_RO legs, **one step lagged and frozen over the increment**.
  - SAS-ME calls the base `GetElasticModuli` at every stage (U9, `LadrunoSANISANDSasME.cpp:259, 403, 508, 570`), so
    RO's overloads would never run.
  - Making them `virtual` is forbidden: it would start running RO elasticity in every existing RO deck (LEDGER_quirks,
    "Adding `virtual` to `ManzariDafalias::initialize()` or `GetElasticModuli`…").
  - **So the decay must be a flag seam in the base's three `GetElasticModuli` overloads**, like `m_PreElastic` and
    `mUseCurrentVoidRatioInG`, fed by `LadrunoSANISAND`.
- **Vanilla RO has five defects the port must not inherit** (§2.9; recorded in LEDGER_quirks by this PR):
  - no Python registration;
  - an absolute reversal threshold that silently ignores small reversals at low p;
  - trial-time reference writes that revert cannot undo, and an empty `sendSelf`;
  - an isotropic initial reference, so a K0 start is taken as prior shear;
  - the K0 patch drift.

---

## 2. Design

### 2.1 The elastic law

```
    G      = G_max(p, e) / T                       K = 2(1+ν)/(3(1−2ν)) · G      (constant ν, §2.6)
    G_max  = B_g · p_at · F(e) · (p̃/p_at)^n        p̃ = max(p + p_re,e, p_min)    (the existing -pRe / -Pmin floor)
    G_fl   = DM04's own G (G0, frozen e_init, (p̃/p_at)^0.5): the vanilla expression, verbatim
    T_max  = max(G_max / G_fl, 1)                  (decay floor = DM04's calibrated G)
    χ_r    = sqrt(½ (r − r_SR):(r − r_SR)),        r = s/p,  r_SR the reference stress ratio
    η_r    = a1 · G_max(p_SR) · γ1(p_SR) / p_SR,   a1 = κ / (κ + T_max − 1)      (RO's η1, a1 derived)
    γ1(p)  = γ1,ref · (p/p_at)^m_γ                 (m_γ = 0: vanilla RO's p-independent reference strain)
    T      = clamp( 1 + (T_max − 1) · (χ_r / (M·η_r))^(κ−1),  1,  T_max ),   M = 1 first shear, 2 after a reversal (Masing)
```

**The floor is DM04's G, by construction.** This is the design choice that protects the peak.
- Beyond the decay range (T = T_max) the elastic law **is** DM04's calibrated law, so the mid- and large-strain response
  returns to the one that matches Tatsuoka at 49 kPa.
- Vanilla RO reaches the same floor only by tuning a1: its clamp is 1 + κ(1/a1 − 1), so the floor ratio is
  a1/(a1 + κ(1 − a1)), not a1. With a1 0.49 and κ 2 that floor is 0.3245·G_max = 29.26 MPa at 49 kPa, against
  DM04's 29.18 MPa (0.27 %; `el_table.py`).
- Here a1 is **derived** per point from T_max(p, e). With n ≠ 0.5, G_max/G_fl varies with p, and a fixed a1 would put the
  floor off DM04's law everywhere except at one pressure.
- **RO parity** holds for n 0.5 and F = the RO form: the same B_g, γ1, κ as a G_RO leg reproduce RO's T to 0.3 % at
  the floor (gate G1).

**n and the void-ratio function F(e).**
- n defaults to 0.5. It acts on **G_max only**, never on G_fl and never on b0 (§2.7).
- F(e) choices:
  - `ro`: 1/(0.3 + 0.7e²) (Hardin 1978, vanilla RO's form). **The default**, for RO parity.
  - `hardin c_e`: (c_e − e)²/(1 + e). c_e = 2.17 is the Toyoura form of Iwasaki et al. (1978); 2.97 is DM04's.
  - Cross-check: B_g 750 with the `ro` form gives 128 MPa at e 0.635, p 100 kPa. Iwasaki's form,
    900·(2.17 − e)²/(1 + e)·σ0′^0.4 kgf/cm², gives ≈ 128 MPa there too.
- **G_max uses the CURRENT e.** G_fl keeps DM04's frozen `m_e_init`: it is the calibrated law and must stay verbatim
  (the `mUseCurrentVoidRatioInG` note, `ManzariDafalias.h`). The monotonic footing's e moves a few percent, so this
  changes G_max by a few percent. The choice is recorded, not hidden.
- **γr now scales with p.** RO's η1 = a1·G_max·γ1/p makes the reference strain a1·γ1 **independent of p**, whereas
  published curves have γr ∝ p^0.35–0.5 (Darendeli 2001; the proxy uses 0.5). m_γ opens that. Its default is 0, which
  gives RO parity.

**Stage 0 is untouched.**
- In `mElastFlag == 0` the base uses a constant G (no p factor), which keeps the K0 patch at 1e-12.
- The decay acts at stage 1 only.
- This is also why the port does not show RO's reported K0-patch drift of 1e-2 (§2.9).

### 2.2 The reference state and robust reversal detection

**State** (committed, with trial twins):
- r_SR (6) and p_SR: the reference;
- r_tp (6), p_tp and χ_max: the turning-point candidate, i.e. the state at the largest χ_r since the reference;
- M ∈ {1, 2}: first shear or after a reversal;
- a count of reversals.

**The rule.** It runs once per material update, before SAS-ME, from the COMMITTED state and the step's Δε only, so
Newton iterates cannot accumulate anything:
1. **Probe.** Take the elastic trial of the whole Δε with the committed reference and compute χ_r,trial.
2. **Reversal**, when χ_r,trial < χ_max − c_rev·η_r:
   - the reference moves to the committed turning point (r_SR := r_tp, p_SR := p_tp), which is the true Masing turning
     point, not the state one step past it;
   - M := 2 and χ_max := 0;
   - the update then runs with the new reference.
3. **Held.** When χ_max − c_rev·η_r ≤ χ_r,trial < χ_max, nothing moves. It is counted (`gmaxRevHeld`).
4. **Consolidation.** While χ_max < c_rev·η_r (no shear yet), p_SR follows the committed p. This is RO's "first shear"
   η1 refresh, made explicit. Proportional compression at constant r (K0 gravity) therefore never degrades G.
5. **At commit**, when χ_r > χ_max: χ_max := χ_r, r_tp := r and p_tp := p.

**Why this is robust, against RO's `(χ_e − χ_en)·Δχ_n < −1e-14`** (`ManzariDafaliasRO.cpp:173`):
- RO's test is **absolute, in strain², on the product of two increments**. A reversal whose increments are below 1e-7
  is silently ignored.
  - Measured on `bin_ro`, 5 kPa, K0 start: a reference seat with δ 1e-3 keeps G/G_max 0.610 at γ 1e-6, the unseated
    value (`engine/k0_ro_k0d001_5.json` against `k0_ro_iso_5.json`). With δ 0.05 it gives 0.998 (`k0_ro_k0d05_5.json`).
    The same seat works at 49 kPa (`k0_ro_k0_49.json`: 1.001).
- Here the threshold is **relative**, in units of the decay's own scale η_r (dimensionless stress ratio). It is
  therefore the same at every p. It is also a **hysteresis**, so noise at p → 0 (where r = s/p is ill-conditioned) cannot flip the
  reference back and forth.
- c_rev is proposed at 0.05 (§10 Q2). The census counts reversals and held candidates per leg.
- **A reversal inside one large increment** is caught at the next step, as in RO. At the footing's ds this is one
  step of lag on a reversal, which a monotonic push rarely has. It is counted.

**The p floor.** p_SR and p_tp are floored at p_min, as RO floors p_SR. With n < 1, η_r ∝ p^(n−1+m_γ) grows as p → 0,
so the decay at the free surface is gentle rather than singular. `gmaxState` (§2.8) reports η_r per point.

### 2.3 Commit, revert, getCopy and the wire

- **Trial vs committed.** Every reference quantity has a committed twin. `commitState` copies trial → committed.
  **`revertToLastCommit` restores it**, which fixes RO's defect: RO writes `mSigmaSR`/`mDevEpsSR`/`mEta1`/`mIsFirstShear`
  at trial time in `integrate()` (`:173-183`) and has no revert override. `revertToStart` re-seats per §2.5.
- **getCopy.** The options travel in a `LadrunoGmaxOptions` struct copy, the same path as `LadrunoSasOptions` (R1,
  #894). The committed reference is copied too.
  - Test: a `getCopy("PlaneStrain")` and a `getCopy("ThreeDimensional")` point run the same trajectory as the
    prototype, with the reference seated.
- **Wire.** RO's `sendSelf`/`recvSelf` are empty (`ManzariDafaliasRO.cpp:231, 237`), so its restart and MP runs are
  unsafe. Here:
  - the options and the committed reference state are appended to the LWIRE block **after #894's slots** (and after
    WP-154's, if its C++ merges first);
  - the layout tag `kLadrunoSanWireTag` is bumped;
  - `static_assert(LWIRE_SIZE != 97)` is re-checked. FE_Datastore keys a Vector by its size (LEDGER_quirks, WP-151).
  - Test: a database round trip value-checks the reference and the options mid-path, then continues bit-identically.
- No statics, so it stays thread-safe (WP-131/146). `mElastFlag` is the base's static and is only read.

### 2.4 Where it lives in the code (the seam)

- **`ManzariDafalias.h`** (vanilla, Ladruno block) gets a `LadrunoGmax` struct: the options and the trial/committed
  reference, plus `bool enabled` (false in every constructor).
- **The three `GetElasticModuli` overloads** get one branch each, guarded by `mLadrunoGmax.enabled && mElastFlag != 0`,
  which computes §2.1. **Flag false: the existing expression runs verbatim**, with `sqrt`, never `pow(x, 0.5)`, whose
  last bit may differ.
  - This is a vanilla-ledger row update.
  - **No `virtual` anywhere** (the RO shadow hazard, ADR-86 §4.1).
- **The reference rule (§2.2)** lives in `LadrunoSANISAND`'s trial update, before `ladrunoSasIntegrate()`.
- **The K0 seat (§2.5)** goes in `ladrunoRunStageFlipOnce()`, which reaches both flip paths, the `updateParameter`
  fast path and the lazy per-instance trigger (`LadrunoSANISAND.cpp:~4825`).
- **The parser:**
  - it **refuses** `-sasGmax` outside `IntScheme 129`. ModifiedEuler freezes G over the increment and CPPM has its own
    moduli chain, so neither would integrate the decay correctly;
  - it refuses it with `-implex`, as the #894 and WP-154 options are refused;
  - it refuses non-finite or out-of-range values: B_g > 0, 0 < n < 1, γ1 > 0, κ ≥ 1, c_rev ≥ 0.

### 2.5 The K0 seat, built in

- **At the stage flip** (`updateMaterialStage 0 → 1`), `ladrunoRunStageFlipOnce()` seats:
  - r_SR := r_tp := r(σ_committed) and p_SR := p_tp := p(σ_committed);
  - χ_max := 0 and M := 1.

  This happens at the same moment LadrunoSANISAND already seats α on the committed gravity ratio.
- **A deck that starts in stage 1 from `-sigma`** gets the same seat at construction.
- This removes the driver's `--ro-rev-delta` unload/reload trick (`footing_ab.py` v3b). Without it, RO's isotropic
  initial reference (`initialize()`, `:253-258`) reads the K0 ratio as prior shear, and G starts fully degraded:
  `k0_ro_iso_49` gives 0.335·G_max at 1e-6.
- **The #894 separation.**
  - A point that enters SEPARATED re-seats to the isotropic state at p_min: r_SR = 0, M = 1, and T = 1 while separated.
    Its tangent is C_e at p_min with G_max(p_min), i.e. #894's "the model's own moduli floor", evaluated with n.
  - At re-contact the reference is re-seated at the re-contact state, and p_re = p_min + K(p_contact)·g uses the same
    law.
  - Recorded, because it changes #894's numbers when the decay is ON.

### 2.6 Energy, consistency, and how K is tied to G

**The concerns, stated.**
1. DM04's own elasticity (G ∝ √p, constant ν) is **hypoelastic and non-conservative**. Closed stress cycles can
   generate or dissipate energy (Zytynski et al. 1978; Einav & Puzrin 2004). The decay inherits this.
2. The decay adds a **path-dependent reference**, so the "elastic" response is history-dependent. Under the Masing
   rules (a reset at the turning point, M = 2) closed strain cycles dissipate, and that is the small-strain damping,
   ≥ 0. A reference reset that does **not** close a loop can generate energy.
3. At a reversal G jumps from G_max/T to G_max. The stress is continuous (rate form); the stiffness is not.
4. The tangent ignores ∂G/∂σ, so it is not the consistent tangent of the rate form.
5. T has a C⁰ kink at T_max, and for κ < 2 an infinite slope at χ_r = 0.

**The choice.**
- **Hypoelastic rate form**, as in DM04 and Papadimitriou & Bouckovalas (2002). No hyperelastic potential.
- The hyperelastic route (Houlsby, Amorosi & Rojas 2005, which handles any n) couples p and q in the elastic law. It
  would change the default formulation and still would not make the decay conservative.
- Gate G9 **measures** the energy:
  - net work ≥ 0 over closed strain cycles at 3 amplitudes (1e-5, 1e-4, 1e-3) and 3 pressures;
  - the drift over closed stress cycles, reported.
- κ is restricted to κ ≥ 1, and κ = 2 is recommended (T linear in χ_r at the origin).
- The deck runs TanType 0 (the elastic C_e at the committed state, modified Newton), so point 4 costs iterations,
  not correctness. The tangent jump at a reversal is counted.

**K and G.** K = 2(1+ν)/(3(1−2ν))·G with **constant ν**, as in RO and PB2002. K decays with G.
- The alternative, decaying G only and keeping K = K_max(p), would change DM04's volumetric unloading and the K0
  behaviour, and has no support in the RO evidence.
- **ν matters here:** the deck pushes at ν 0.05 (DM04 Table 1), so K = 0.78·G. A rigid strip's vertical stiffness
  scales roughly with G/(1 − ν). §8 D5 asks whether the small-strain ν should differ, and E8 measures the sensitivity.

### 2.7 SAS-ME, R1, the cutoff and R3b

**SAS-ME's elastic predictor** (`ladrunoSasElastic`, exact today because G = g·√x and K = c·G):
- **n ≠ 0.5, decay OFF:** the closed form generalizes, with x^(1−n) = x0^(1−n) + (1−n)·c·g·t·dv on the power branch and
  the floor branch unchanged.
- **Decay ON:** G depends on r through T, so there is no closed form. The predictor becomes an **error-controlled
  sub-integration** of the elastic ODE (embedded Heun, TolR), not the fixed 64 Heun steps of the
  `mUseCurrentVoidRatioInG` fallback. It is counted as `gmaxPredSubsteps`.
- The yield intersection (`ladrunoSasIntersect`) calls the same predictor.

**Stages and error control.**
- Every stage already evaluates K and G at its own state (U9), so the decay enters every stage with **no new stage
  logic**. The reference is fixed within an increment.
- The substep error test (σ, α, z) sizes substeps where G changes fast: first loading from the reference, and the T
  kink. The cost increase is measured (`sasStats` substeps, E2 against E0).

**R1 (floor, hysteresis, soft cap).**
- R1 acts on h and is independent of G, with one exception: **the soft cap's bound H ≥ κX uses X, the elastic part of the
  loading denominator**, and X scales with the stage's current, decayed G.
- That is the consistent reading ("X is the elastic part"), so it needs no code change. G7 reports `hSoftCapped`
  against the OFF leg.

**The cutoff (#894):** §2.5.

**R3b (WP-154), η = 2G·τ.** **Which G: recommended G_max(p, e) (T = 1), not the decayed G.**
- η is a viscosity scale. Tying it to a path-dependent reference would make it jump ×T_max at every reversal and
  make the Deborah number history-dependent.
- G_max keeps η ∝ p^n, a function of the current state only.
- WP-154 re-derives τ after this WP anyway (its §2.4 step 4). The decision belongs to WP-154 (§8 D6).

**b0 (the plastic modulus): stays on G0, as now and as in RO.**
- DM04: h = b0/((α − α_in):n), with b0 = G0·h0·(1 − c_h·e)·(p/p_at)^(−1/2). G0 here is DM04's calibrated constant, not
  the elastic G.
- Keep G0 and the −½ exponent. Do **not** use G_max or the decayed G:
  1. b0 is a calibrated plastic constant. G_RO kept it on G0 and kept ε_peak at 1.80 % and the peak (`ro_49` =
     `van_49`).
  2. b0 on G_max would multiply h by 3.08, i.e. **h0 ×3.08, the lever already rejected at ×1.5** (E50 1.66×, ε_peak
     1.43 %).
  3. b0 on the decayed G would make Kp depend on the elastic reference: two memories acting on one modulus.
- **Stated plainly:** the plastic response still changes, because G enters the loading denominator through the
  2G(B − C·tr n³) − K·D·(n:r) term and the elastic part of every stage. That is the measured E50 +6 % (`ro_49` against
  `van_49`). G3 bounds it.

### 2.8 Flags, options and census

**Flags:**

| flag | meaning | default |
|---|---|---|
| `-sasGmax B_g n γ1 κ` | the decay ON (IntScheme 129 only) | OFF |
| `-sasGmaxVoid ro` \| `-sasGmaxVoid hardin c_e` | F(e) | `ro` |
| `-sasGmaxRefExp m_γ` | γ1(p) ∝ p^m_γ | 0 |
| `-sasGmaxRevTol c_rev` | the reversal hysteresis, in units of η_r | 0.05 |

The echo prints the setting and T_max at p_at.

**`sasOptions`** (33101) gets these values appended after #894's (and WP-154's). The indices are fixed at merge order.

**`sasStats`**, appended:
- `gmaxReversals` and `gmaxRevHeld`, per update and committed;
- `gmaxPredSubsteps`;
- `gmaxLastT` and `gmaxMaxT`, committed;
- `gmaxTangentJumps`.

**A new response, `gmaxState`** (id **33103** proposed, after WP-154's proposed 33102; to be confirmed free at
implementation). At the committed state it returns:
- G_max, G, T, χ_r, η_r;
- r_SR (6), p_SR, M and the reversal count.

It gives TIMs the **G/G_max map of the footing**, which is the input the TIM stiffness recursion with G/G0 needs.

**The replay tool.** `ladrunoSANISANDReplay` needs the reference state in its input, so a replay can start from a
committed reference.

### 2.9 Vanilla RO's defects: what the port fixes, and what is recorded

| # | defect | evidence | in the port |
|---|---|---|---|
| 1 | **No Python registration.** RO is registered for Tcl only (`TclModelBuilderNDMaterialCommand.cpp:600`) | `bin_ro` needed +2 lines (README) | This WP adds the 2-line Python registration (a vanilla-ledger row), so the parity legs run on the release build. RO's numerics are untouched |
| 2 | **An absolute reversal threshold**, `(χ_e − χ_en)·Δχ_n < −1e-14` (`:173`) | a seat with δ 1e-3 does nothing below p′ ≈ 5 kPa (§2.2) | relative hysteresis in units of η_r (§2.2) |
| 3 | **The reference is written at trial time and not restored by `revertToLastCommit`; `sendSelf`/`recvSelf` are empty** (`:173-183`, `:231`, `:237`) | source reading | trial/committed twins, revert, wire (§2.3) |
| 4 | **The K0 start reads as prior shear.** The initial reference is isotropic (`:253-258`) | `k0_ro_iso_49`: 0.335·G_max at 1e-6 | seated at the stage flip (§2.5) |
| 5 | **K0 patch drift 1e-2 instead of 1e-12**, from G ∝ √p in the elastic stage | relayed 2026-09-30, not re-measured here | stage 0 untouched (§2.1). O1 re-measures RO's patch to confirm the mechanism |

Also (§1.3): under ModifiedEuler, RO's decay is applied **one step lagged and frozen** over the increment. The
G_RO numbers therefore carry that integration error, and gate G1 compares the port with the **oracle** RO, not only with the
engine.

LEDGER_quirks gets one entry for these (this PR). Fixing vanilla RO beyond the registration is out of scope, because it would
change existing RO decks.

### 2.10 Class tag: none

This is a **model option of `LadrunoSANISAND`'s SAS-ME**, like R1 (#893), the separation (#894) and R3b (WP-154). There
is no new class, no wire class and no `classTags.h` edit.

---

## 3. The small-strain plastic modulus: a follow-up, conditional, and instrumented here

**The question.** Should this WP also make the plastic modulus depend on strain amplitude near a reversal, as later
SANISAND variants do? SANISAND-MS (Liu, Abell & Pisanò 2019) makes h large inside a memory surface; vanilla
`SAniSandMS` is in `UANDESmaterials/`.

**Recommendation: no. It is a follow-up WP, opened only on the evidence of this one.**
1. **Where the deficit is.** With the decay alone (RO, iso start):
   - at 49 kPa: 0.920 against 0.955 at 1e-5 (−4 %), **0.556 against 0.721 at 1e-4 (−23 %)**, and on target at 1e-3
     (0.240 against 0.237);
   - at 5 kPa: 0.826 against 0.882 at 1e-5, and **0.440 against 0.474 at 1e-4 (−7 %)**.

   The footing's early soil is at p′ 1–20 kPa, where the deficit is ≤ 7 %. It **cannot explain a footing at 0.53–0.61
   of the test.** The exponent n (+26 % G_max at 10 kPa) and the boundary-value error (G4) are larger candidates, and
   both are in this WP.
2. **It touches what already matches the lab.** An amplitude-dependent h stiffens 1e-5–1e-4 only if it decays back
   to DM04's h before mid-strain. Otherwise it is the rejected h0 lever (E50 1.66×, ε_peak 1.43 %).
   - That is a new hardening law with its own memory, its own calibration and its own interaction with R1's floor and
     hysteresis (all of which act on h).
   - It would also move WP-154's Kp_crit.
3. **The proxy target is itself medium-low confidence** at 1e-4 (§1.2). A published Toyoura curve (G2) comes first.
4. **Scope.** The owner asked for the stiffness without breaking the peak or the element match. One modulus at a time
   keeps the attribution clean: the decay (elastic) here, h (plastic) later, if needed.

**What this WP does instead.** The oracle (O2) records the elastic and plastic shares of the secant at
γ 1e-5 / 3e-5 / 1e-4 / 3e-4, at p′ 5 / 20 / 49 kPa. **The trigger for the follow-up:**
- after E3/E4, the footing secant at s/B 0.005–0.02 is still below 0.8× the test with n 0.4;
- the G4 boundary-value error is < 5 %;
- and the element deficit at 1e-4 on the published curve (G2b) is > 0.1 at p′ ≤ 20 kPa.

If all three hold, the owner opens the h-follow-up with that decomposition as its brief.

---

## 4. What is claimed, and what is not

- **Claimed, if the gates pass:**
  - with the declared (B_g, n, γ1, κ, m_γ) calibrated on **element data only**, the Kimura footing's secant at s/B
    0.005–0.02 is concave-down and within the agreed tolerance (§8 D4);
  - the peak, the Tatsuoka match and Gate 0 are unchanged within §5's tolerances;
  - the FE/boundary-value share of any remaining gap is quantified (G4).
- **Not claimed:**
  - a fabric effect;
  - that n = 0.4 is TIMs' sand;
  - the small-strain plastic response (§3);
  - anything about the post-peak (DM04 softens too little and dilates too much; §9).
- **The honest-framing rule** (as WP-154's τ): the decay constants are **never** tuned to the footing. The footing
  is the validation, not the calibration. If G5 fails at the element-calibrated constants, the WP reports and the owner
  decides (§8). It does not iterate on γ1.

---

## 5. Gates

| id | gate | measure | pass |
|---|---|---|---|
| **G0** | the off-switch | flags absent: WP-151's 643-replay baseline plus WP-152's paths | bit-identical |
| **G0b** | continuity | decay ON with B_g, F and n chosen so that G_max ≡ G_fl (T_max = 1): the Gate 0 paths | ≤ 1e-12 relative to OFF (pow vs sqrt) |
| **G1** | RO parity (the port reproduces the evidence) | n 0.5, F `ro`, B_g 750, γ1 5e-4, κ 2, m_γ 0, iso start: the element paths of `oracle/ro_*` (ro_oracle.py) | max \|Δq\|/q_max ≤ 2e-3; G/G_max at 1e-6…1e-3 within 0.005. The engine `bin_ro` differences are reported and explained (one-step lag, §2.9) |
| **G2a** | element G/G_max–γ, decay range | the TOTAL secant G/G_max (elastic + DM04 plastic) on a drained PS/simple-shear path from a K0 and an iso start, p′ 5 / 20 / 49 / 100 kPa, γ ≤ 3e-5, against a **published Toyoura curve**: Iwasaki, Tatsuoka & Takagi (1978), S&F 18(1):39–56, and Kokusho (1980), S&F 20(2):45–60, transcribed at O2 | \|ΔG/G_max\| ≤ 0.05 |
| **G2b** | … the plastic range | same, γ 1e-4 – 1e-3, with the elastic/plastic shares | reported. Pass: no worse than vanilla RO at the same p (0.556 at 1e-4, 49 kPa) and within 0.05 at 1e-3. The deficit is the §3 trigger |
| **G3** | element peak and mid-strain | Gate 0 (DM04 Table 1: the 17 Verdugo & Ishihara triaxials) and the Tatsuoka set (PSL-24, Fig. 6a, Fig. 16a, Fig. 8, the point tests), at the tests' own e | q_peak within 1 %, φ′_peak within 0.3°, ε_peak within ±10 % of plain DM04 at the same e, E50 ≤ 1.10× plain DM04 (RO measured 1.06×), R at ε_peak + 2 % within 2 %. The rejected h0 ×1.5 leg fails it (ε_peak 1.43 against 1.80 %, E50 1.28×), so the gate discriminates |
| **G4** | **the analytical elastic check** | plasticity off: an `ElasticIsotropic` FE with G(z) from the deck's own K0 stress field (G_max(p(z)), the deck's ν, depth and width, B/8 and B/16), and a very wide/deep copy. Against: **Gibson (1967)**, Géotechnique 17:58–67 (G ∝ z, ν ½: exact, Winkler-like); **Booker, Balaam & Davis (1985)**, IJNAMG 9:369–381 (strip on G ∝ z^α); **Gazetas (1991)**, JGE 117(9):1363–1381 (a homogeneous stratum over rigid base, the FE sanity leg). A plane-strain homogeneous half-space has no finite static strip stiffness, so the homogeneous comparison uses the stratum | FE within 5 % of the solution at B/16; the B/8 error and the finite-domain error are reported as E_BVP, which is subtracted from the G5 interpretation |
| **G5** | **the footing initial stiffness** | secant q/(s/B) at s/B 0.002 / 0.005 / 0.01 / 0.02, against Kimura V85.6 (45.4 / 41.0 / 37.2 MPa at 0.005 / 0.01 / 0.02), and the curvature sign | concave-down (the secant falls monotonically from 0.005 to 0.02) **and** within ±20 % at 0.005–0.02 (proposed, §8 D4). The 0.002 secant is reported only (digitizing ±0.002 s/B, §1.2) |
| **G6** | the footing peak | q_peak and s/B at the peak, against the SAS-ME base with R1 + cutoff (the same deck, OFF), and Kimura's 1 952 kPa at s/B 0.091 | q_peak unchanged within 3 %, or closer to the test; no new refusal class; the R1 counters and #894 entries reported |
| **G7** | composition | `hFloored`, `hSoftCapped`, `reseatHeld`, #894 entries, substep counts: ON against OFF on the same leg | 0 `loadingNonPosH`; the changes reported |
| **G8** | n sensitivity | n 0.4 against 0.5 at an equal G_max(p_at) | the Δ secant at s/B 0.005 / 0.02 reported, with its G4 counterpart (the elastic share of Δ) |
| **G9** | energy | the element closed cycles (§2.6) | net work ≥ 0 over closed strain cycles; the stress-cycle drift reported |
| **G10** | mesh | the chosen leg at B/8 against B/16 | secant difference ≤ 5 % at s/B 0.005–0.02 |
| **G11** | provenance | per leg: the flags, B_g, n, F, γ1, κ, m_γ, c_rev, the build (`ladrunoBuild()`), the driver, the output path, the date, and `gmaxReversals` | complete |

**Reliability of the G2 references** (stated, not hidden):
- Both are **cyclic** tests (torsional shear; triaxial), at σ′ ≥ ~20 kPa. The exact ranges are transcribed at O2 from the
  papers, which were not re-read for this plan.
- The footing's early zone at 1–20 kPa is therefore an **extrapolation**.
- Comparing a monotonic backbone with a cyclic secant assumes Masing behaviour.
- A K0 start carries an initial shear stress that the cyclic tests (isotropic consolidation) did not have.
- The fork gets TIMs' own sand data if TIMs supply it (§8 D2).

---

## 6. Oracle first, then C++

| phase | content | exit |
|---|---|---|
| **O1 — oracle** | Extend the WP-134 oracle (`Ladruno_scripts/sanisand_reference/`) as a testbed copy with toggles, the WP-151 pattern (`Ladruno_files/testbed/sanisand_reseat_r1/sanisand_r1/`). Add `gmax=GmaxOptions(B_g, n, void, γ1, κ, m_γ, c_rev)` to `elastic_moduli`/`quantities`, the reference state to `State`, and the §2.2 rule **between** increments (the reference is fixed within one, as in the C++). Workbench `src/ro_oracle.py` is the seed (it patches `M.quantities`). Radau integrates the rate equations exactly | (i) OFF reproduces WP-134 exactly; (ii) n 0.5 / `ro` / B 750 / a1 0.49 / κ 2 reproduces `oracle/ro_*` to 1e-10; (iii) analytic anchors: with plasticity suppressed (h0 → ∞), the first-shear backbone at constant p integrates in closed form, and the Masing unload branch is 2× the backbone; (iv) RO's K0 patch drift re-measured (§2.9 #5) |
| **O2 — calibration and decomposition, no C++** | Transcribe the G2 Toyoura curves and fit (γ1, κ, m_γ) at n 0.4 and 0.5 on the **total** element response, with B_g and F from Iwasaki et al. (1978). Record the elastic/plastic shares (§3). Preview G3 on the Tatsuoka and Gate 0 paths | the constants of E3/E4, the §3 decomposition, and the G3 preview |
| **O3 — the footing's states** | The WP-151 fan method (`cxx_fan.py` pattern): committed footing GPs from the ladder checkpoints × 32 directions × 2 magnitudes, with the O2 constants. No new refusal class | 0 new refusals; a substep-cost estimate |
| **O4 — the analytical elastic check (in parallel, no C++)** | G4 on the current release build: `ElasticIsotropic` per element from G_max(p(z)) of the deck's K0 field, n 0.5 and 0.4, ν 0.05 (and 0.2 for D5), B/8, B/16 and a wide/deep copy; the three references | E_BVP known **before** any footing leg. If E_BVP is large, the footing gap is partly the mesh or the domain, and G5's reading changes |
| **C1 — C++** | §2.4: the seam, the reference rule, the stage-flip seat, the predictor (§2.7), parser, `sasOptions`, `sasStats`, `gmaxState`, commit/revert/getCopy/wire, the replay input, and the RO Python registration | G0/G0b; the C++ matches the oracle **step by step** on O1/O3: benign states ≤ 2e-7 relative at TolR 1e-7 (the WP-129 standard), element paths ≤ 2e-3·q_max; getCopy, the wire round trip and the two-instance interleave; a **revert test** (a reversal in a trial, then `revertToLastCommit`: the reference is restored bit-for-bit) |
| **C2 — tests** | `tests/test_ladruno_sanisand_gmax.py` (fast tier): off-switch, oracle match, the reversal hysteresis at p′ 0.5 / 5 / 49 kPa (the RO δ 1e-3 case as a regression), the K0 seat, revert, wire, the parser refusals (non-129, `-implex`, ranges), and the #894 re-seat | green on Zone-A. Mutation gate: drop the revert restore ⇒ red; replace η_r-relative by RO's absolute threshold ⇒ red; ignore n ⇒ red |
| **E — Esmeralda** | §7 | G5–G8, G10, G11 |

---

## 7. Esmeralda test matrix and effort

The deck is `deck_toyoura_ladder`, the Kimura case: B 0.9 m, e 0.635, γ 15.9, K0 0.5, push ν 0.05, B/8.
- Every leg runs R1 full (`-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5`) plus `-sasTensionCutoff` at #894's reviewed
  values.
- Every leg writes `gmaxState` at checkpoints.
- The references are the ladder legs on build `8ebde5cbd` (L base, `G_RO`).

| # | leg | to s/B | purpose | gates |
|---|---|---|---|---|
| E0 | new build, flags OFF | 0.02 | the build control against L base | G0 (on the curve) |
| E1 | elastic FE, n 0.5 / 0.4, B/8 / B/16 / wide-deep, ν 0.05 (+ 0.2) | small | O4 on the deck geometry | G4 |
| E2 | decay ON, **RO parity** constants (n 0.5, `ro`, B_g 750, γ1 5e-4, κ 2) | 0.02 | the port against `G_RO` (K0 seated by the flip, not the driver) | G1 (footing), G7 |
| E3 | decay ON, n 0.5, O2 constants | 0.02 | calibration effect alone | G5, G8 |
| E4 | decay ON, **n 0.4**, O2 constants (± m_γ) | 0.02 | the exponent | G5, G8 |
| E5 | the chosen leg (E3 or E4) | **the peak** (0.2) | the peak | G6, G7 |
| E6 | the chosen leg, B/16 | 0.02 (0.05 if cheap) | mesh | G10 |
| E7 | the chosen leg, `deck_toyoura` B 1.2 m b8 | the peak | the peak at the second geometry, against `W_TYR_b8` (2 509 kPa at s/B 0.1798) | G6 |

**Cost.**
- The stiffness legs stop at s/B 0.02, a few hours each on B/8.
- The peak legs run about 8 h each (B/8).
- The B/16 legs are the long ones: WP-154 measured b16 at 20.7 h to s/B 0.066.

**Effort, about 3 weeks:**
- O1–O3: 4–5 days;
- O4: 1–2 days, in parallel;
- C1–C2: 5–7 days;
- Esmeralda: 3–5 days wall, plus the analysis.

| # | risk | mitigation |
|---|---|---|
| R-1 | spurious reversals at the free surface (r = s/p ill-conditioned as p → 0) | the η_r-relative hysteresis; the p floor; `gmaxReversals`/`gmaxRevHeld` per leg; C2's low-p case |
| R-2 | substep cost (the predictor is no longer closed-form with the decay ON) | error-controlled, counted; E2 against E0 |
| R-3 | the ellipticity picture moves (G enters Kp_crit/2G) | WP-154 re-derives τ afterwards (its §2.4 step 4) |
| R-4 | the test's own uncertainty at the start: the toe, the ~16 kPa footing weight, roughness, side friction, the g-field | G5 tolerance (§8 D4); red-team items A5, C3, C6, C7 carried into the report |
| R-5 | no published Toyoura curve below ~20 kPa | reliability stated (§5); TIMs' own data (§8 D2) |
| R-6 | ν 0.05 at small strain | E1 at ν 0.2; D5 |
| R-7 | calibration creep (tuning γ1 to the footing) | §4 honest-framing rule; the constants are fixed at O2 |

---

## 8. Owner and TIMs decisions

1. **D1 — n.** 0.4 (Iwasaki et al. 1978, Toyoura at γ ≈ 1e-6), 0.5 (RO, DM04), or TIMs' own. It moves G_max by
   +26 % at 10 kPa.
2. **D2 — the source of G_max for TIMs' sand.** Kimura/Toyoura is the benchmark, but the macroelement is for TIMs'
   target sand. **The TIMs input needed: their target sand's small-strain stiffness data, i.e. G_max or a Vs profile
   (bender elements, seismic CPT or cross-hole), and a G/G_max–γ curve (resonant column or torsional shear), with
   the confining pressures and void ratios of each test.** Without it, the fork calibrates to published Toyoura
   curves and says so.
3. **D3 — which stiffness the macroelement needs** (handoff §13, stiffness track item 3). The choices are the initial
   tangent, a secant at a reference s/B, or the whole G/G0 degradation. It decides whether G5 is read at s/B 0.002,
   0.005 or over a range, and whether `gmaxState` maps are a deliverable.
4. **D4 — the G5 tolerance.** ±20 % is proposed. The red team puts the combined uncertainty at the peak at ±15 %,
   and the start adds the toe, the footing weight and K0 (about 10 % on the early secant, red-team D9).
5. **D5 — the small-strain ν.** Keep DM04's 0.05 for the push, or use a different value? It moves the strip
   stiffness through G/(1 − ν). E1 measures it. The decision is the owner's.
6. **D6 — which G in R3b's η** (with WP-154). G_max is recommended (§2.7).
7. **D7 — adopt or not.** After §7, is the decay part of the TIMs setting, like R1's full set? If yes, WP-154's τ
   is re-derived on it.

---

## 9. Out of scope, and the ordering

- **Out of scope:**
  - fabric and bedding (candidate b);
  - the test conditions (c), beyond reporting them;
  - the small-strain plastic modulus (§3, a conditional follow-up);
  - **R3b** (WP-154, #895);
  - **the explicit campaign** (WP-153);
  - the post-peak dilatancy and softening;
  - fixing vanilla `ManzariDafaliasRO` beyond its Python registration (the defects are recorded in quirks);
  - the stage-0 elastic law.
- **Ordering relative to WP-154.** Stiffness is decided in the first s/B ≈ 0.02. On the Toyoura deck the B/16−B/8 gap
  opens at s/B 0.0335 (WP-154 §1.1), so **stiffness is decided before localization, and the two WPs are
  independent.** Neither waits for the other.
  - **If this WP is adopted (D7), WP-154's τ procedure is re-run on the new model**, because Kp_crit/2G and G move.
    WP-154 §9 already provides for this.
  - The wire and `sasOptions` slots are appended in merge order after #894's (and WP-154's if its C++ merges
    first). The later WP rebases.

---

## 10. Left open — flagged, not guessed

- **Q1 — the published curves' details.** The pressure ranges, void ratios and test types of Iwasaki et al. (1978) and
  Kokusho (1980) were not re-read for this plan (Zotero was unavailable). O2 transcribes them before fitting.
- **Q2 — c_rev.** 0.05·η_r is a proposal. It needs O1 on a small-cycle path at p′ 0.5–5 kPa.
- **Q3 — Booker, Balaam & Davis (1985).** That a strip on G ∝ z^α (0 < α < 1) has a finite static stiffness in plane
  strain is expected but must be confirmed from the source before G4 relies on it. Gibson (1967) covers α = 1.
- **Q4 — RO's K0 patch drift (1e-2)** was relayed, not re-measured here. O1 confirms the mechanism.
- **Q5 — G_max on the current e, G_fl on the frozen e_init.** This mixes the two conventions on purpose (§2.1). The owner
  may prefer both frozen, for a clean DM04 floor at every e.
- **Q6 — the Masing factor.** M = 2 after the first reversal is RO's and PB2002's rule. The monotonic footing hardly
  exercises it, so the cyclic gate (G9) is its only check.
- **Q7 — WP-155.** WP-155 is taken by the pile-contact R0.5 branch (`wp/pile-contact-r05`, `155_pile_contact_r05.md`,
  no PR yet), so this WP is **156**.
