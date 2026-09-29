---
title: "WP-152 — SAS-ME tension cutoff (separation) for zero-confinement points"
project: Ladruno
type: plan + opt-in implementation
status: "PLAN (owner GO, relayed by the TIMs orchestrator 2026-09-29: 'let's do both' — the tension cutoff is the priority fix for the free surface, the p_r bracket stays the cross-check). Building on this branch; the owner reviews the evidence before any merge."
owner: nmora
related:
  - "[[151_sanisand_reseat_singularity]] (R1, #893: this builds on its SAS-ME code)"
  - "[[93_ladruno_sanisand_zero_confinement_adr]] (II.3 consistent low-p cutoff, deferred; III.3 erosion, deferred)"
  - "[[84_ladruno_mc_tension_cutoff_adr]] (the fork's Rankine cutoff pattern, ASDPlasticMaterial3D)"
  - "[[150_sanisand_regularization_memo]] (#892: the p_r -> 0 extrapolation, the acceptance cross-check)"
  - "[[LadrunoSANISAND_implex_guide]] §13"
tags: [plan, sanisand, sas-me, tension-cutoff, separation, free-surface, tims, wp-152]
updated: 2026-09-29
---

# WP-152 — SAS-ME tension cutoff (separation)

## Why

After R1 (WP-151) the binding limiter of a realistically dilating sand under the TIMs footing is the FREE
SURFACE. Two surface Gauss points about 0.23–0.27B outside the footing edge go to p′ → 0 and refuse. The refusal
codes are `errorAtDTmin`, `tensionAtDTmin(lowP)` and `maxSubsteps` (the B/4 leg: 34 × code 4, 310 × code 9 and
2 × code 6 at one point). This stops:
- every Toyoura leg, at s/B 0.010–0.032, with 0 NonPosH and R1 on = off;
- every physically bounded PB/PB2 leg, at s/B 0.009–0.014, even with the campaign deck's 7.65 kPa surcharge.

Deleting elements would remove the passive wedge's weight, which N_γ relies on, and would be mesh-dependent. The
owner's choice is the FLAC / AEM-spring analogue done at the material level. The point loses stiffness and
strength but keeps its weight and its place in the mesh, and the change is reversible.

## Facts it rests on

Citations are to this branch; a research pass was run on 2026-09-29.

- **SAS-ME has ONE low-p test: p + p_r > 0**, with p = tr(σ)/3 + m_Presidual (LadrunoSANISANDSasME.cpp:394, 501,
  845, 885). `-Pmin` is not a SAS-ME threshold: it floors the elastic moduli only (SAS:305-341).
- **The low-p refusals** (SAS:140-161, 668-699, 760-772):
  - code 3 `startInadmissible`, when the committed p0 ≤ 0 (or a non-finite value / a trace check);
  - code 6 `tensionAtDTmin`, when a stage or the predictor reaches p ≤ 0 at dT_min;
  - code 4 `errorAtDTmin` and code 9 `maxSubsteps`, which at low p come from the rate equations' fast α
    evolution (the α and fabric error terms do not scale with p).
- **Vanilla hides this.** `Stress_Correction` sets σ = (p_min + p_r)·I and α = 0 silently when
  tr(σ)/3 < p_min (ManzariDafalias.cpp:3181-3254). ModifiedEuler abandons the increment at dT_min with no refusal
  (MD:1983-1994).
- **ADR-93** deferred a consistent low-p cutoff (II.3) and element erosion (III.3). TIMs' D1 rule: `-Pmin` ≤
  0.5 kPa, `-Presidual 0`, and report the limit load at the floor F and at F/2.
- **Precedents:**
  - Roy et al. 2020, a DM04 footing to a peak in Abaqus: on a negative multiplier, dT/4; at dT ≤ 1e-7, an
    elastic update; "a p′ floor of 0.001 kPa freezes the state";
  - Tejchman & Herle 1999: tensile stress → zero stress and the initial stiffness;
  - FLAC's tension cutoff;
  - the fork's own Rankine cutoff (ADR-84).

  Via the TIMs literature report, [E-sec].

## Design

A per-Gauss-point state machine, SAS-ME only, OPT-IN: `-sasTensionCutoff p_sep p_contact`.

### NORMAL → SEPARATED (entry)

The trigger masks ONLY low-p/tension refusals:
- **E1:** code 6, or code 3 whose cause is p0 ≤ 0. Tension, at any p.
- **E2:** code 4 or 9 with the committed p0 < p_sep. A low-confinement accuracy or cost failure. p_sep = 0
  disables E2 and leaves a pure tension cutoff.

Everything else still refuses loudly, at any p:
- code 5 `loadingNonPosH` (the α_in singularity);
- code 2 (α outside the bounding surface);
- code 3 from a non-finite value or a trace check;
- code 7 (drift) and code 8 (α at dT_min);
- bad input.

Each entry is counted by cause.

### SEPARATED

- Stress: σ = (p_min − p_r)·I, i.e. model p = p_min. No tension, no shear. The weight is untouched (body forces
  act on the nodes).
- α := 0 and α_in := 0 on entry: the only α consistent with an isotropic state. Fabric z is kept, as vanilla does.
- Void ratio: follows the strain, as always.
- Strain increments are absorbed. The volumetric opening since entry is recorded: g = tr ε − tr ε_entry
  (compression positive).
- Tangent: C_e at p_min, the model's own moduli floor (about 1 % of the P_atm moduli at p_min = 0.01 kPa). This is a
  Newton regularisation: the consistent tangent of a constant stress is 0. `tangentEP` returns the same.

### SEPARATED → NORMAL (re-contact)

- Condition: g ≥ g_c = (p_contact − p_min)/K(p_contact), i.e. the gap has closed with some overlap.
- Then σ := p_re·I, with p_re = p_min + K(p_contact)·g ≥ p_contact, and α = α_in = 0. SAS-ME resumes at the next
  increment.
- Hysteresis: p_contact > p_sep, and re-entry needs a qualifying refusal. So a point cannot chatter between the
  two states.

### Energy

- Entry releases the stored elastic energy at p ≤ p_sep: small, dissipated, ≥ 0.
- While separated, the work is p_min·Δε_v, about 0.
- Re-contact is an elastic reload from p_min, so it is conservative.
- Gate: the net work over closed separation cycles is ≥ 0.

### Interfaces

- **Flag:** `-sasTensionCutoff p_sep p_contact` in model stress units (kPa in the campaign).
  - Requires p_sep ≥ 0 and p_contact > max(p_sep, p_min).
  - SAS-ME only. Refused with `-implex`, which is unqualified.
  - OFF by default, and then byte-identical.
  - Starting values: p_sep 0.5 kPa (TIMs D1's p-floor bound), p_contact 1.0 kPa.
- **State:** `sepActive` and tr ε_entry. They are committed, reverted, copied and sent on the wire.
  - The LWIRE block grows; `static_assert(LWIRE_SIZE != 97)` still holds.
  - The layout tag goes 151 → 152.
  - `sasOptions` gains p_sep and p_contact at indices 9–10.
- **Census:** `sasStats` gains `sepEntriesTension`, `sepEntriesLowP`, `sepExits` and `sepActive` (the committed
  0/1). **The length goes 36 → 40**, so consumers must read it from the response.
  - Transitions are counted ONCE, at commit: the update records the trial's transition (`sepEvent`), and
    `LadrunoSANISAND::commitState` counts it.
  - This matters because an element may call the update several times per step, each call starting from the
    committed state. The first build counted a re-contact twice.
- **Echo:** states the setting.

## Evidence gates

1. **Oracle first.** The state machine is a DEFINITION. It is implemented in a driver around the unmodified WP-134
   oracle, whose `p_floor` stop is the exact counterpart of code 6. Element tests:
   - isotropic unloading to separation and back;
   - triaxial extension to separation and back;
   - cyclic open/close: no chatter, net work ≥ 0.

   The C++ must match the driver step by step: the entry step, σ, the re-contact step and p_re. The NORMAL phases
   must agree to SAS-ME's tolerance.
2. **C++ tests:**
   - OFF is byte-identical (the WP-151 643 replays);
   - only E1/E2 are masked (NonPosH, non-finite and so on still refuse);
   - commit/revert and a database round trip, with separation state on the wire;
   - the parser.
3. **Footings** (the orchestrator launches them on Esmeralda):
   - Toyoura B/8 R1 + cutoff vs Kimura: a peak near s/B ≈ 0.09?
   - PB2 B/8 R1 + cutoff vs the classical band;
   - campaign R1 + cutoff identical to R1 wherever no point separates;
   - the R3 peak test (B/8, B/16, shear:15);
   - p_sep vs p_sep/2: TIMs D1's F vs F/2 rule, < 2 %.
4. **Cross-check** against WP-150's p_r → 0 linear extrapolation (`pr_extrapolate.py`, #892 e74cba7ed), on the TYR
   p_r 1/2/5/10 and PB2 p_r 0/2/5 legs. A mismatch is a finding to explain, not a failure.

## Decision point (owner/TIMs)

Separation is a constitutive choice about near-surface sand. Open choices:
- p_sep and p_contact;
- whether fabric is kept or reset;
- the stiffness while separated;
- how separated points are reported in the capacity.

The alternatives stay available as cross-checks: the p_r bracket, `-pRe`, surcharge, and element removal (ADR-51).

## Effort

About 2.5–3 agent-days to the first SHA for the footing gates:

| task | effort |
|---|---|
| oracle driver + element tests | ½ day |
| C++ | 1 day |
| tests | ½ day |
| build and regression | ½ day |
| docs and ledgers | ½ day |

Ownership: the WP-151 (R1) session. The merge gate is as for #893: an adversarial review, Zone-A, and a local
Windows log.
