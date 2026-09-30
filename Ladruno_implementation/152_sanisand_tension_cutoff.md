---
title: "WP-152 — SAS-ME tension cutoff (separation) for zero-confinement points"
project: Ladruno
type: plan + opt-in implementation
status: "IMPLEMENTED, draft #894; advisor review 2026-09-29 answered (a0171df75, 7f1562c81: gated entry, masked-code census, parser refusals, continuous re-contact); footing re-checks on the review build running on Esmeralda. The owner merges after the adversarial re-review and Zone-A."
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

## Implementation status and results (2026-09-29)

The C++ is complete at 8cdc0bdba (draft #894). A Windows build of all 5 targets passes:
- 20 SANISAND/ManzariDafalias test files: 242 passed, 5 skipped, 2 xfailed, 0 failed;
- the static gates.

**Oracle first** (`Ladruno_files/testbed/sanisand_tension_cutoff/tc_oracle.py`, fixture
`tests/data/wp152_oracle_paths.json`). The three element paths start from the C++ post-flip state (σ = 2.0036 kPa
isotropic), with R1 (the full set) and the cutoff at p_sep 0.5, p_contact 1.0:

| path | oracle events (step) | C++ | net work |
|---|---|---|---|
| isotropic out and back | entry 4 (E1), re-contact 96 | the same steps, E1 | ≥ 0 |
| triaxial extension and back | entry 14 (E1), re-contact 110 | the same steps; the entry is **E2** | ≥ 0 |
| 4 isotropic open/close cycles | 4 entries, 4 re-contacts, alternating | the same steps, all E1, no chatter | 0.00147 (oracle 0.00148) |

- Separated steps: σ = p_min·I exactly. After re-contact the C++ agrees with the oracle to ≤ 4e-6.
- On the extension path the entry is E2: SAS-ME's accuracy/cost limit at p0 < p_sep comes first, in the same step
  where the exact trajectory reaches p = 0.
- Before separation on that path, below p ≈ 1 kPa, the deviator differs by ≤ 0.0125 kPa. That is SAS-ME's absolute
  error floor (σ_ref = 1 kPa, `-errFloor`), pre-existing and not the cutoff.

**It touches nothing else.** On the 643 WP-151 replays with the cutoff given:
- the 533 accepted updates are bit-identical;
- all 110 refusals keep their code (101 α_in singularity, 8 α outside the bounding surface, 1 accuracy failure at
  high p);
- none of them is low-p, so none separates.

**Other C++ checks:**
- The separation state survives a database save and restore: a restored point is still separated, and the
  continuation is bit-identical.
- A separated point's `tangentEP` is C_e at p_min.
- The census counts once per committed transition. The first build counted per update call (a re-contact counted as
  2); this is recorded as a quirk.
- The parser refuses: p_sep < 0, p_contact ≤ p_sep, p_contact ≤ p_min, a non-129 scheme, and `-implex`.

**Next: the footing gates** (the orchestrator, on Esmeralda), then the adversarial review and Zone-A.

## Review response (advisor review 2026-09-29: MERGE WITH CHANGES)

Each finding was checked against the code at 8cdc0bdba before it was acted on. Fixes: a0171df75 (gates, census,
parser, ISA) and 7f1562c81 (continuous re-contact, found by the new Newton test). Session: claude-code, 2026-09-29.

| # | finding | verdict | response |
|---|---|---|---|
| 1 | E2 checks only the committed p0; the B/8 top row sits below p_sep in situ, so an accuracy failure under COMPRESSION separates | **confirmed**, with one sub-claim wrong: `TR_REJ_REVERSAL` is set only above dT_min (`:709-710`) and re-seats at it, so it never returns RC_DTMIN | E2 also needs a non-compressing increment (tr Δε ≤ 0, compression positive; this also keeps the elastic predictor's p ≤ p0 < p_sep). Held refusals refuse, counted `sepHeldCompressing`. Tested both ways (code 9 via `-maxSubsteps 1` at p0 ≈ 0.4 kPa). Whether E2 is needed at all: the p_sep = 0 footing legs below. |
| 2 | E1 separates on code 6 at ANY committed p0 | **confirmed** | E1 only at committed p0 ≤ p0max (`-sasSepMaxP0`, default 5·p_contact); above it the code 6 refuses (`sepHeldHighP`). `sepMaxP0` records the largest p0 at a committed entry. Oracle carries the bound (`refused_highp`). |
| 3 | the masked refusal code is lost | **confirmed** | `sepLastCode` (the masked code, set at commit), `sepMaxP0`; tests assert code 6 on an isotropic entry. |
| 4 | with p_r ≠ 0 the separated σ is tensile and K is taken at tr σ/3 + p_Re | **confirmed** (`GetElasticModuli` reads `tr(σ)/3 + m_PreElastic`, ManzariDafalias.cpp:5480) | refused: the cutoff requires `-Presidual 0` (parser and `recvSelf`). The campaign (TIMs D1) uses p_r = 0. |
| 5 | revertToStart under ISA clears the separation but keeps σ | **confirmed**; no footing deck uses ISA (`footing_ab.py` has none) | `ladrunoResetSasSep` split from the census reset, called outside ISA and by the replay command. Found alongside: after "ISA off" the domain update integrates −ε_n for EVERY SANISAND point (pre-existing; a NORMAL point's trial p 1.769 → 1.731 kPa), LEDGER_quirks. |
| 6 | re-contact is volumetric only | **confirmed, by design** | documented (guide §13.5, echo); tested: 50 isochoric steps keep a separated point at p_min, the volumetric closing re-contacts it. |
| 7 | re-contact lands on α = α_in = 0, the 1e10 sentinel without the h floor | **confirmed** (`ladrunoSasBracketH` :228-229) — note it is DM04's own post-reversal state | the cutoff requires `-sasHFloor > 0` (parser and `recvSelf`). |
| 8 | a separated cluster cannot equilibrate; use a force/energy test | **agreed** | the footing driver already uses `NormUnbalance` (a force test); the census now reports separated points in/out of the footprint per step. The `EnergyBalance` recorder is velocity-based and reads all zeros under a static integrator (measured), so it is not an energy test for the push (LEDGER_quirks). |
| 9 | tests never run Newton, plane strain, shear-while-separated, E2 under compression, p_r ≠ 0 | **confirmed** | added all five (+ E1 bound, ISA). The Newton test FOUND A DESIGN FLAW: the first build's jump p_min → p_contact at re-contact leaves a band of top displacement with no equilibrium (a two-brick column failed at step 68). Fixed (7f1562c81): while separated, p = p_min + K(p_contact)·max(g, 0), continuous; the exit state and the oracle events are unchanged. |

**Regression.**
- Windows, the nmora desk (the 8ebde5cbd build was made on this desk as well, in this worktree, so no log is owed by
  another desk): `Ladruno_scriptsuild.bat`, ALL 5 targets, `build 7f1562c818b3` (log
  `build_wp152_review_all_7f1562c81.log`, 2026-09-29 21:02). 20 SANISAND/ManzariDafalias test files:
  **257 passed, 5 skipped, 2 xfailed, 0 failed**. The byte-identity test needs the interpreter's site-packages on
  `PYTHONPATH` under `-S` (LEDGER_quirks); with it, it passes bit-exact against the win32 baseline. WP-152 file 30/30.
- Linux (Esmeralda node4, `~/ladruno_wp152/OpenSees_review`, build dir `build/wp152r_seq`, installed to
  `~/ladruno_wp152/bin_review/opensees.so`, logs `~/ladruno_wp152/logs/build_review2_srun.out`): 28 test files,
  **253 passed, 8 skipped, 2 xfailed, 1 failed** — the failure is the byte-identity child process losing `pytest`
  under `python -S` in the conan venv (an environment artifact, LEDGER_quirks), not a row mismatch.
- The 643 WP-151 replays, R1 vs R1 + cutoff (the cutoff now needs the h floor, so the stored R1-OFF baseline no
  longer applies to it): 635 accepted updates bit-identical, 8 refusals keep their code, 0 separate.

## What the cutoff actually is: E2 is the operative trigger (footing evidence, 2026-09-29)

**At p → 0 the free-surface wall reaches SAS-ME as an accuracy/cost failure (codes 4 and 9) BEFORE any tension.** The
α and fabric error terms do not scale with p, so the substep controller hits dT_min or the substep cap while p is
still positive, and the trajectory never gets to cross p = 0. The design choice is therefore:

> a point whose update fails its accuracy or cost limit at committed p < p_sep, under a non-compressing increment, is
> treated as separated.

E1 (the tension trigger, code 6) is a backstop. The name "tension cutoff" is historical: this is a **low-confinement
separation**.

Evidence (build 7f1562c81, `~/ladruno_wp152/bin_review`, driver `deck_review/footing_ab.py`, Toyoura, R1 + cutoff,
2026-09-29 20:55):
- With p_sep = 0 (E2 off), EVERY leg stops at FLOOR at s/B ≈ 0.0044:
  - B/8 at 0/1.0, 0/0.5 and 0/2.0: q 71.6 kPa, step 38 (`runs/R_TYR_b8_ps0*`);
  - fig9 (B 0.9, Dr 0.856) at 0/1.0: q 61.8 kPa, step 41 (`runs/R_TYR_fig9_ps0`).
- In each case exactly two mirror-symmetric points refuse: B/8 at (±0.87, −0.03) m, fig9 at (±0.65, −0.02) m. Both
  are in the top GP row, just outside the footing edge.
- The codes are 4 and 9 only (B/8: `errorAtDTmin` 56, `maxSubsteps` 50), with **zero code 6**. Nothing separates, so the
  p_contact spread cannot be read at p_sep = 0.
- p_sep = 0 is not viable. The smallest workable p_sep is the defensible default; it is being measured (0.1 / 0.25 /
  0.5, below).

## Footing re-checks on the review build (2026-09-29/30)

Build 7f1562c81 (`~/ladruno_wp152/bin_review/opensees.so`), driver `~/ladruno_wp152/deck_review/footing_ab.py` (the
Toyoura deck with the 44-column census and a per-step `sep_census.csv`), Toyoura, R1 + cutoff, target s/B 0.20, one leg
per node. Output `~/ladruno_wp152/deck_review/runs/R_*`, report `~/ladruno_wp152/review_report.py`, read 2026-09-30
10:27. The pre-fix counterparts are the 8ebde5cbd gates in `deck_toyoura/runs/W_*`.

**Peak and post-peak against the pre-fix build:**

| leg | p_sep / p_contact | peak q (kPa) @ s/B | q @ 0.19 | q @ 0.20 | mode | pre-fix (8ebde5cbd) peak @ s/B, q @ 0.20 |
|---|---|---|---|---|---|---|
| R_TYR_b8_ps05 | 0.5 / 1.0 | 2508.4 @ 0.1795 | 2483.9 | 2446.0 | TARGET | W_TYR_b8: 2508.8 @ 0.1798, 2452.3 |
| R_TYR_b8_ps025_pc05 | 0.25 / 0.5 | 2510.2 @ 0.1800 | — | 2456.0 | TARGET | W_TYR_b8_half: 2513.1 @ 0.184, 2456.8 |
| R_TYR_fig9_ps05 | 0.5 / 1.0 | 2016.6 @ 0.1673 | 1896.1 | 1747.1 | TARGET | W_TYR_fig9_856: 2015.4 @ 0.1672, 1747.8 |
| R_TYR_b8_ps05_pc2 | 0.5 / 2.0 | 2517.2 @ 0.1834 | — | 2458.9 | TARGET | — |
| R_TYR_b8_ps01 | 0.1 / 1.0 | 2507.9 @ 0.1789 | — | 2443.1 | TARGET | — |
| R_TYR_b8_ps0(_pc05, _pc2), R_TYR_fig9_ps0 | 0 / 1.0, 0.5, 2.0 | FLOOR at s/B 0.0044 (71.6; fig9 61.8) | | | FLOOR | — |

- **The review fixes did not move the peak or the post-peak.** The E2 compression gate and the continuous re-contact
  change the peak by −0.02 % (B/8), −0.12 % (the ½ setting) and +0.06 % (fig9); the s/B at the peak by ≤ 0.004; q at
  s/B 0.20 by ≤ 0.26 %. Kimura test V (fig9 conditions): 1953 kPa @ 0.091.
- **p_sep 0.1 / 0.25 / 0.5: 2507.9 / 2510.2 / 2508.4 — spread 0.09 %.** p_sep = 0.1 kPa reaches the target, so the
  smallest workable p_sep tested is 0.1 (the largest committed p0 at an entry in that leg: 0.094 kPa).
- **p_contact spread at the peak:** 0.5 (at p_sep 0.25) 2510.2, 1.0 2508.4, 2.0 2517.2 (+0.35 %).
- **Every leg runs its whole push on E2:** 55–56 E2 entries (masked codes: 9 about 45, 4 about 10) and 0–2 E1 (code 6)
  per leg. The largest committed p0 at an entry stays below p_sep (0.09–0.48 kPa). The compression gate holds back 25–64
  qualifying refusals on B/8 and 224 on fig9, all refused and cut by the driver. The substep cap still refuses
  1450–2170 times per leg (step cuts, none stopping a leg).

**Where the separated points are at the peak** (`sep_census.csv`, `field_step*.npz` positions):
- B/8 (B = 1.2 m, edge at |x| = 0.6): 52–56 GPs separated, **none in the footprint**, |x| 0.63–1.32 m (0.03B–0.6B
  beyond the edge), down to y = −0.42 m (0.35B). **Not the top row only:** 38–40 of them lie below the second GP row —
  the heaving passive wedge, not a surface skin. 0–4 re-contacts by the peak, 4–8 by s/B 0.20.
- fig9 (B = 0.9 m, edge at 0.45): 58 GPs at the peak, none in the footprint, |x| 0.47–0.99 m, down to y = −0.31 m
  (0.35B). Post-peak, at s/B 0.20: 90, of which **2 in the footprint**.

**Reading.** The capacity does not depend on the cutoff's parameters (0.09 % over p_sep, 0.35 % over p_contact) nor on
the review fixes (≤ 0.12 %). But the separated zone is a region of the wedge about 0.35B deep, not a few surface points.
Whether a separated wedge of that size is an acceptable constitutive idealisation is the owner's call. The p_r → 0
cross-check (≤ 2.3 %, pre-fix) bounds its effect on the peak.
