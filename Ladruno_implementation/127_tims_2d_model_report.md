---
title: "Why the Ring Walls: the fork's reply to the TIMs 2D-model requests F18–F23"
project: Ladruno
status: sixth issue, 2026-10-03; ready for the owner's review before sending; not sent (see the revision log at the end)
date: 2026-10-03
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

**Reply to the TIMs 2D-model intake of 25 September 2026, items F18–F23. Sixth issue, 3 October 2026;
first issue 28 September 2026.**

The reply answers every item of the intake. Each statement is one of four kinds: **shipped** (merged on
the fork's `ladruno` branch), **measured** (a value with a committed script, a build and a date),
**pending** (an open change or a run in progress) and **TIMs decision**. Sections 0 and 0a give the state
on 3 October 2026. Where an earlier statement in §§1–6 differs from them, §§0–0a govern; §§1–4 keep the
detailed evidence of the earlier issues.

| | |
|---|---|
| Fork | Ladruno / OpenSees, branch `ladruno` |
| Status of the fork changes | §0a.7, the single status table of this reply |
| Plan and findings A–D | `Ladruno_implementation/127_tims_2d_requests_plan.md` |

Every line of §2 of the intake was checked against the source before any work started
(`127_tims_2d_requests_plan.md`, preamble: all citations held on `fb1afe58b`; no commit had touched
`SRC/material/nD/UWmaterials/` since the intake's `79e062367`). One reading was refined later, not
refuted: the "abrupt switch at 0.5 kPa" in the error norm is a continuous 1 kPa floor (§1, F18(a)).

---

## 0. Summary

1. **The numerical limit is closed.** The footing deck stopped for three successive reasons, each now
   removed by a measured fix.
   - The `ModifiedEuler` wall at s/B 0.0292 came from discrete integration defects (§2). The corrected
     integrator, **SAS-ME** (`IntScheme 129`), removes it.
   - The SAS-ME wall at s/B 0.0508 came from a singularity of the DM04 plastic modulus at an α_in
     re-seat (§4.6). **R1**, an opt-in floor on h with a hysteretic re-seat, removes it.
   - The free-surface limit at s/B 0.01–0.03 came from points at zero confinement in the passive
     wedge. The **low-confinement separation**, an opt-in rule under which such a point carries neither
     tension nor shear and keeps its weight, removes it.

   With the three in place the completed legs reach their target of s/B 0.20, and an independent dynamic-relaxation
   solution reproduces the implicit one to s/B 0.091 within 0.4 % in q (§0a.1).
2. **The footing reaches a peak on Toyoura sand.** With DM04's published Toyoura constants on the
   centrifuge geometry of Kimura et al. (1985), the peak load is 1.03 times the test at D_r 85.6 % and
   1.19 times at D_r 75.1 %. The settlement at peak is 1.6 to 1.8 times the test's. The peak load
   complies with the ±15 % band at D_r 85.6 % and exceeds it at 75.1 %; the timing, the initial stiffness
   and the mechanism do not comply (§0a.2).
3. **The limit load is not mesh-converged.** On Toyoura with B 1.2 m, the B/16 load at s/B 0.20 is
   about 20 % below the B/8 peak, still rising and with no peak. On the Kimura case, at s/B 0.183 the B/16
   leg carries 1757 kPa against the B/8 leg's 2296 kPa peak, still rising (leg running, read
   2026-10-02). The bands lose ellipticity well before the peak and follow the mesh
   lines. Perzyna viscoplasticity inside the material is therefore a prerequisite for a capacity of
   calibration grade, not a refinement (§0a.3).
4. **The initial stiffness is a deliverable, and it is currently about 0.44 of the test.** The secant to
   s/B 0.02 is 16.4 MPa against 37.2 MPa, and the model curve is concave-up where the test is
   concave-down. The remedy is a G_max decay law with a configurable pressure exponent inside the
   material (§0a.4).
5. **Post-peak dilatancy.** DM04 dilates 2.3 to 2.7 times more than the tests of Tatsuoka et al. (1986)
   at 8 % axial strain, and softens less. Under an assumed band volume fraction of 0.3, localisation inside the
   laboratory specimens accounts for about half of that gap; the fraction is not measured, and "about
   half" is a bound under that assumption. Recalibrating the critical-state line comes first; a state-dependent dilatancy
   option follows only if the recalibration is not enough (§0a.5).
6. **The campaign SANISAND parameter set is a cyclic fit applied to a monotonic problem.** No footing
   curve computed with it is called physical. Two sets are offered for drained monotonic loading in its
   place: a physically bounded monotonic set, and DM04's own Toyoura calibration, the set compared with
   Kimura et al. (1985) (§0a.6).
7. **Inputs still owed by TIMs** are listed in §5.1. The small-strain stiffness of the target sand
   (G_max or a Vs profile, and a G/G_max curve) is added in this issue; the initial stiffness cannot be
   closed on the target sand without it.

**The intake items on 3 October 2026.**

| item | request | state | detail |
|---|---|---|---|
| F18(a) | error norm with an absolute floor | shipped as `-errFloor` in SAS-ME; not done in `ModifiedEuler`, where it changes nothing measurable | §1 |
| F18(b) | rate-form stages | superseded by SAS-ME (shipped) | §1 |
| F18(c) | CPPM under a global Newton | shipped, with the corrected tangent sign as the default | §1 |
| F18(d) | per-point fallback `ModifiedEuler` to CPPM | shipped (`-meFallback cppm`); under 10 % of the material time on this deck | §1 |
| F18(e) | the b8 ring point | answered: an inadmissible state; SAS-ME refuses it | §1, §2 |
| F19 | SANISAND in the threaded loop | inventory shipped; code pending; SANISAND stays outside the threaded loop | §1 |
| F20(a)–(c) | counters, profile scopes, `tangentEP` | shipped | §1 |
| F21 | replay of a dumped state | shipped | §1 |
| F22 | deterministic mode | shipped; not yet run on Esmeralda | §1 |
| F23(a)–(b) | PDMY03 constants, the dilation brake | shipped; the brake cannot produce a plateau | §1 |

The integrator recommendation, SAS-ME with R1 and the separation, is in §4.11. It is a decision of the
owner and of TIMs.

---

## 0a. State on 3 October 2026

**Sources.** Esmeralda runs of the fork's copy of the TIMs plane-strain strip deck, SAS-ME at TolR 1e-4,
B/8, unless a row says otherwise. Every value carries its build (fork commit) and the date it was read.
"Toyoura" is the DM04 calibration of Dafalias & Manzari (2004), Table 1. "Kimura case" is the surface
strip footing of Kimura et al. (1985), Fig. 9, B = 0.9 m prototype (30g centrifuge), modelled at
prototype scale.

### 0a.1 The numerical limit

| limit | s/B reached | cause | fix | source |
|---|---|---|---|---|
| `ModifiedEuler` wall | 0.0292 (campaign set) | integration defects G, E, F and C (§2) | SAS-ME | `7936ed6e0`, 2026-09-28 |
| SAS-ME wall | 0.0508 (campaign set) | DM04 plastic-modulus singularity at an α_in re-seat (§4.6) | R1 | `7936ed6e0`, 2026-09-28 |
| free surface | 0.01–0.03 (every realistically dilating set) | updates failing at p′ → 0 in the passive wedge, 0.225B outside the footing edge | low-confinement separation | `8ebde5cbd`, 2026-09-29 |

- **R1** carries the campaign set past its wall to s/B 0.093, at q 1459 kPa and still rising. Each part of
  R1 alone walls earlier, at s/B 0.045 to 0.0525. The cap strength κ 0.25 / 0.5 / 0.75 gives identical
  curves to s/B 0.045 and differs by at most 5 % after, not in κ order; κ = 0.5 is recommended. Source
  `bd93c558d`, 2026-09-29.
- **The c ≥ 0.78 route is withdrawn as a cure.** At footing scale c = 0.80 without R1 walls at s/B 0.048,
  on the compression side (`bd93c558d`, 2026-09-29). The recommendation c ≥ 7/9 stays for Lode convexity
  only (D8).
- **The separation is a recognition, not a cause.** At the peak 52 to 58 Gauss points are separated,
  none under the footing, from 0.03B to 0.6B beyond the edge and down to 0.35B deep: the heaving passive
  wedge. Over the 50 steps before entry p′ falls at 97–100 % of them and all dilate; 79–91 % were already
  below 0.5 kPa when the rule acted. About 40 % of the entries are seeds and only 6–11 % follow a
  neighbour within 10 steps, so the zone is not a numerical cascade. It is a dilatant heave driven by the
  model. Points at p′ ≈ 0 carry no shear, so the general-shear slip surface cannot reach the ground
  there. Source `7f1562c81`, 2026-09-30. Whether that dilation is physical is the subject of §0a.5.

**Verification of the closure.**

| check | value | reference | difference | source |
|---|---|---|---|---|
| reviewed separation, Toyoura B 1.2 m, peak | 2507.9 kPa at s/B 0.1773 | 2508.4 kPa at s/B 0.1795, previous build | −0.02 % | `77c454660` against `7f1562c81`, 2026-10-01 |
| reviewed separation, Kimura case to s/B 0.05 | step and separation records | previous build | byte-identical | `77c454660`, 2026-10-01 |
| review fixes, Kimura case D_r 85.6 %, peak | 2017 kPa at s/B 0.167 | 2015 kPa at s/B 0.167, pre-review | +0.1 % | `7f1562c81` against `8ebde5cbd`, 2026-09-30 |
| build change, Kimura case e 0.635, implicit, peak | 2293.4 kPa at s/B 0.1695 | 2296 kPa at s/B 0.170 | −0.1 % | `2b5cb2f8f` against `8ebde5cbd`, 2026-10-01 |
| p_sep 0.1 / 0.25 / 0.5 kPa | peak spread | | 0.09 % | `7f1562c81`, 2026-09-30 |
| p_contact 0.5 / 1.0 / 2.0 kPa | peak spread | | 0.35 % | `7f1562c81`, 2026-09-30 |
| separation against residual pressure p_r → 0, extrapolated | q | | within 2.3 % | `8ebde5cbd`, 2026-09-29 |
| dynamic relaxation (DR) against implicit, q at s/B 0.05 / 0.091 | 970.73 / 1701.81 kPa | 974.47 / 1705.83 kPa | −0.38 % / −0.24 % | `2b5cb2f8f`, 2026-10-02 |
| DR against implicit, nodal displacement at s/B 0.091 | largest difference 0.45 mm | settlement 85.8 mm | 0.5 % | `2b5cb2f8f`, 2026-10-02 |
| DR at 5e-8 m per step, committed α_in re-seats to s/B 0.005 | 438 | implicit 230 | 1.9 times; q within 0.2 % | `2b5cb2f8f`, 2026-09-30 |

- The review fixes change no footing result. The pre-review label of the fifth issue is lifted, and the
  values from build `8ebde5cbd` stand.
- DR converges to the implicit path: each halving of the rate cuts the re-seat excess about 9 times. The
  remaining excess sits at surface points approaching separation, before they separate. In the loaded
  mechanism zone DR commits fewer re-seats than implicit (22 against at least 193 at s/B 0.091), and the
  points with the excess carry about 1e-5 of the internal work. The DR comparison at the peak is
  pending.

### 0a.2 The peak on Toyoura, against Kimura et al. (1985)

**Setup.** DM04 Toyoura with SAS-ME, R1 and the separation (p_sep 0.5 kPa, p_contact 1.0 kPa), B/8, on the
Kimura geometry. Two tests of Fig. 9 are used, both of series V, loaded across the bedding: V85.6 (D_r
85.6 %, model e 0.652) and V75.1 (D_r 75.1 %). The tests loaded along the bedding are not used, because
DM04 has no fabric and its constants come from specimens loaded across the bedding. Fig. 9 is digitised
in full; two independent reads agree within 30 kPa and give identical peaks. q is the gross pressure
over B.

| quantity | test | model | ratio | verdict | source |
|---|---|---|---|---|---|
| peak, V85.6 | 1953 kPa at s/B 0.091 | 2017 kPa at s/B 0.167 | 1.03 | complies, within ±15 % | `7f1562c81`, 2026-09-30 |
| peak, V75.1 | 1236 kPa at s/B ≈ 0.105 | 1469 kPa at s/B 0.165 | 1.19 | exceeds | `7f1562c81`, 2026-09-30 |
| s/B at peak, V85.6 / V75.1 | 0.091 / ≈ 0.105 (series 0.09 to 0.16) | 0.167 / 0.165 | 1.8 / 1.6 | exceeds | `7f1562c81`, 2026-09-30 |
| secant q/(s/B) to s/B 0.02, V85.6 | 37.2 MPa, concave-down | 16.4 MPa, concave-up | 0.44 | does not comply | `8ebde5cbd`, 2026-09-30 |
| mechanism | general shear | punching captured by the mesh columns | | does not comply | `8ebde5cbd`, 2026-09-29 to 2026-10-01 |
| element: dilation at 8 % axial strain | Tatsuoka et al. (1986) | DM04 | 2.3 to 2.7 | exceeds | exact material-point integration, 2026-09-29 |

- **The ±15 % band** comes from an independent adversarial review of the comparison: mesh orientation
  9 %, e_max/e_min ±5 to 8 %, test scatter ±8 %. The footing roughness of that series is not verified; a
  smooth footing would put the model 25 to 35 % high. The test of the same density loaded along the
  bedding peaks at 1734 kPa at s/B 0.163, so bedding alone moves the test by about as much as the misfit.
- **The element is not too soft before its peak.** Against Tatsuoka et al. (1986), drained plane strain
  at 4.9 to 392 kPa, DM04 matches φ′_peak within −1 to +2° and E50 at 49 kPa with a ratio of 0.99. The
  strain to peak is 0.46 to 0.73 of the test at 5 to 10 kPa and 1.04 at 49 kPa. After the peak DM04
  softens too little (stress ratio 1.35 times the test at peak + 2 %) and dilates too much (table). Exact
  material-point integration, 2026-09-29.
- **q is compared with the test up to the peak only.** Extending the fine mesh zone from 1.5B to 3B deep
  moves the peak by −1.0 % (1996.0 against 2015.4 kPa) but changes the drop from peak to s/B 0.20 from
  −13.3 % to −10.5 %; part of the post-peak softening is made by the mesh grading. Source `8ebde5cbd`,
  read 2026-10-01.
- **On the TIMs footing geometry** (B 1.5 m, q₀ 7.65 kPa) with Toyoura at D_r 0.733 there is no peak by
  s/B 0.20: q is 1938 kPa buoyant and 2303 kPa dry at s/B 0.20, still rising. Source `7f1562c81`,
  2026-09-30. A plateau there needs targets beyond s/B 0.20, or the mechanism fix of §0a.3.

**Consequence.** The peak load is of the right order. Its timing, the initial stiffness and the
mechanism are not. Each discrepancy is assigned to one of three coupled coherence tracks: the
mechanism (§0a.3), the initial stiffness (§0a.4) and the post-peak dilatancy (§0a.5).

### 0a.3 Mesh convergence and the mechanism

| geometry | B/8 | B/16 | difference | source |
|---|---|---|---|---|
| Toyoura, B 1.2 m | peak 2509 kPa at s/B 0.180 | 1966.7 kPa at s/B 0.20, still rising, no peak | about −20 % at s/B 0.20 | `8ebde5cbd`, read 2026-10-02 |
| Kimura case, e 0.635 | peak 2296 kPa at s/B 0.170 | 1757 kPa at s/B 0.183, still rising; leg running | not final: 1757 at s/B 0.183 against the 2296 peak | `8ebde5cbd`, read 2026-10-02 |

- **On the finer mesh no peak forms by s/B 0.20 on Toyoura, and none by s/B 0.183 on the Kimura case,
  where the leg is still running.** The capacity a TIM
  surface would be calibrated from therefore depends on the mesh.
- **Cause: non-associated localisation before the peak.** An acoustic-tensor census of the Toyoura legs
  shows the B/16 to B/8 gap passing 2 % at s/B 0.0335, when 1 to 5 % of the points near the footing lose
  ellipticity, about 5 times before the peak. With associated flow the count of non-elliptic points is
  zero at every peak. Source `8ebde5cbd`, 2026-09-29.
- **The bands are captured by the mesh.** Before the peak they run down the element column at the
  footing edge, where the material's own band direction is 45 to 65° from vertical; after the peak they
  run along a mesh row. A mesh sheared by 15° moves the peak by +8.8 to +9.7 %. A three-fold elastic
  stiffness mobilises the wedge under the footing earlier (s/B 0.05) but leaves the edge bands in the
  mesh columns. Sources `8ebde5cbd`, 2026-09-29 to 2026-10-01.
- **Remedy: Perzyna viscoplasticity inside the material**, opt-in. The census estimates a small viscous
  hardening, H_v/2G ≈ 0.015, so the rate bias and the cost are expected to be modest. A nonlocal void
  ratio would miss the points that lose ellipticity while still hardening. The regularisation is tuned
  after the constitutive changes of §0a.4 and §0a.5, because those change when and how the bands form.
  Status: planned, not built.
- **Proposed acceptance.** Half the test scatter. Kimura's own N_γ scatter gives about ±4 % (D3).

### 0a.4 The initial stiffness

The initial stiffness of the footing is a deliverable: it feeds the calibration of the TIM
macroelement. On the Kimura case (e 0.635) q at s/B 0.005 / 0.01 / 0.02 is 74 / 160 / 349 kPa against the
test's 227 / 410 / 743 kPa (`8ebde5cbd`, 2026-09-30). On V85.6 the secant to s/B 0.02 is 0.44 of the test
(§0a.2).

| variant | effect on the start | effect on the peak | verdict | source |
|---|---|---|---|---|
| initial state: seating 16 / 40 kPa, K0 0.4, surcharge 5 kPa | 10 % or less | | does not comply | `8ebde5cbd`, 2026-09-30 |
| G0 ×3, b0 fixed | matches to s/B 0.01, then above the test; shape concave-up | 2713 kPa at s/B 0.103, 1.39 of the test | exceeds | `8ebde5cbd`, 2026-09-30 |
| G0 ×3, h0 unchanged | concave-up | 3395.1 kPa at s/B 0.20 | exceeds | `8ebde5cbd`, read 2026-10-02 |
| h0 ×1.5 | +7 % at the footing | E50 1.28 of plain DM04 at equal density, against a 1.10 gate | exceeds | element and footing, 2026-09-30 |
| G_max decay (vanilla `ManzariDafaliasRO`, G_max of Iwasaki & Tatsuoka for Toyoura) | secant 22.7 MPa, 0.61 of the test; concave-down | element peak unchanged; the footing run stops at s/B 0.033 (no R1, no separation in the vanilla class) | marginal | vanilla build, 2026-09-30 |

- **The elastic modulus carries about 80 % of the stiffening**; the initial state moves the start by 10 %
  or less. A scaled G0 raises the start but keeps the wrong shape and overshoots the peak. Only G_max
  decay gives the concave-down shape of the test.
- **Most of the remaining gap probably sits at very low confinement.** The G_max of Toyoura scales with
  p′ to a power of about 0.4; the model hard-codes 0.5, which gives 0.79 of the modulus at 10 kPa.
  Fabric, bedding and the test conditions may also contribute.
- **Remedy: an opt-in G_max decay law with a configurable pressure exponent n inside
  `LadrunoSANISAND`.** It keeps SAS-ME, R1, the separation and the Perzyna regularisation, so that the
  stiffness and the limit load come from the same run. It is checked against an elastic strip on a
  modulus increasing with depth (Gibson; Booker et al. 1985; Gazetas 1991) and against Kimura et al.
  (1985). It must not move the peak or the element behaviour that already matches the laboratory.
  Status: planned, not built. The acceptance band for the initial secant is a TIMs decision (D11).

### 0a.5 Post-peak dilatancy

| measure | laboratory (Tatsuoka et al. 1986) | DM04 | ratio | source |
|---|---|---|---|---|
| −ε_v at 8 % axial strain, base set | three tests, 4.9 and 49 kPa | | 2.3 / 2.6 / 2.7 | exact material-point integration, 2026-09-30 |
| drop of stress ratio, peak to 8 % strain | 29 to 34 % | 17 to 27 % | | same |
| A0 ×0.5, n^d 10 to 15, strain to peak held within +10 % | | | 1.8 to 2.4 | same |
| A0 ×0.3, n^d 10, dilation matched | | strain to peak ×1.30 to 1.37; φ′ −0.3 to −0.9°; drop 7 to 12 % | 1.0 | same |
| laboratory shear-band correction under an ASSUMED band volume fraction f ≈ 0.3 (not measured) | | | 1.3 to 1.9 | same date |
| critical-state line recalibrated: e0 0.85, n^b 1.81; φ′, E50, strain to peak and drop held | | | 1.56 / 1.88 / 1.68 | same date |

- **The dilatancy constants A0 and n^d act as one lever and cannot fix both symptoms.** Cutting the
  dilation also cuts the softening, while the laboratory shows less dilation and more softening than
  DM04. One A0 cannot fit both 4.9 and 49 kPa. The fabric constants z_max and c_z are inert on a monotonic
  path.
- **The UW low-pressure dilatancy factor is not the lever.** It acts below p′ = 0.05·P_atm (5.05 kPa),
  only suppresses dilation, and has its half point at 1.05 kPa; the wedge points enter separation where
  it is already about 1e-3, and the laboratory gap sits at 5 to 50 kPa, above its range (D2).
- **Laboratory localisation, under an assumed band volume fraction.** Tatsuoka's specimens form shear
  bands before the peak, and the external strains include them. Under an assumed band volume fraction
  of 0.3, localisation in the specimens accounts for about half of the gap (not measured); "about half"
  is a bound under that assumption, not a measured share.
- **Recalibrating the critical-state line closes part of the rest without code.** The price is a
  critical void ratio of about 0.85 at 4.9 kPa against 0.93 from Verdugo & Ishihara, at the edge of the
  plausible range.
- **Order.** The critical-state recalibration first, element first; then a finite-element biaxial
  specimen; then a footing re-run that checks whether the separated zone shrinks, whether the band
  reaches the surface, and whether the peak moves earlier. An opt-in state-dependent dilatancy option
  only if the recalibration is not enough. Status: planned; the direction is pending the owner's
  decision.

### 0a.6 The campaign parameter set

The campaign SANISAND set is assembled from Gorini's cyclic Messina set, DM04 Toyoura's c and c_h, and
φ 33° with Jaky's K0. In element tests it has the strength of a very dense sand but dilates 20 to 23
times (triaxial) and 8 to 14 times (plane strain) less than stress-dilatancy requires, and it peaks at 4
to 16 % strain (§4, caveat). **It is a cyclic fit applied to a drained monotonic problem**, and no footing
curve computed with it is called physical. Two sets are offered for drained monotonic use:

- **A physically bounded monotonic set (PB2):** n^b 1.652, A0 0.692, n^d 3.5, h0 3.5, D_r 0.47 assumed. It
  liquefies at N ≈ 1 in cyclic tests and is not for cyclic work. On B/8 with the separation it reaches
  1866 kPa at s/B 0.15 without a clear peak, above the classical band of 381 to 1595 kPa for φ′ 33 to 42°
  (`8ebde5cbd`, 2026-09-29).
- **DM04's own Toyoura calibration**, the set behind §0a.2 and the only one so far compared with a
  physical footing test.

Neither replaces a calibration on the target sand (D7, §5.1).

### 0a.7 Status of the fork changes

| change | work package, PR | state on 2026-10-03 |
|---|---|---|
| SAS-ME, counters, replay, `tangentEP`, deterministic PARDISO, PDMY03 items, earlier reports | WP-127 to WP-136; #846, #863 to #866, #869 to #872, #874, #876, #884, #885, #888 | merged |
| CPPM under a global Newton, tangent sign, `-meFallback cppm` | WP-130, #868 | merged (`e96f8d77d`) |
| R1 | WP-151, #893 | merged (`fd87e396d`) |
| DR revert and `commitStats` | WP-153, #899 | merged (`9fcb6cfaa`) |
| `ForwardEuler` (IntScheme 5) defects | WP-158, #901 | merged (`5300da720`) |
| low-confinement separation | WP-152, #894 | merged (`fc2ea4fc8`) |
| sub-step moduli of IntScheme 0, 4, 6, 7, 8, 9 (1, 2, 3, 5, 45 and 129 untouched) | WP-160, #914 | merged (`35dff42ab`) |
| Perzyna viscoplasticity | WP-154, #895 | plan |
| G_max decay with exponent n | WP-156, #898 | plan |
| post-peak dilatancy | none yet | plan, partial |
| footing A/B report; regularisation memo | WP-138, #878; WP-150, #892 | drafts, documents only |
| this reply | WP-147, #887 | draft, for the owner's review |

**For TIMs.**
- The peak load on the reference sand is of the right order, 1.03 to 1.19 of the test, but it is not
  mesh-converged, so no capacity is offered for calibration before the Perzyna regularisation.
- The settlement at peak and the initial stiffness are not yet reliable.
- The inputs of §5.1 remain the gate to a physical curve for the target sand.

---

## 1. Item by item

Status key: **SHIPPED** (merged on `ladruno`), **PENDING** (open PR or run), **THEIR DECISION**,
**NOT DONE** (with the reason).

### F18(a): an error norm with an absolute floor

**Asked.** `err = ‖dσ₂ − dσ₁‖ / max(2‖σ‖, σ_ref)` as a flag, and substeps and error against a
reference for σ_ref ∈ {0, 0.1, 1, 5} kPa, and where the admitted error is below the Newton tolerance.

**Answer.** The proposed norm is what `ModifiedEuler` already computes, with σ_ref = 1 kPa. The
"switch at 0.5 kPa" is continuous: `‖dσ₂−dσ₁‖` below ‖σ‖ = 0.5 and `/(2‖σ‖)` above it is exactly
`/max(2‖σ‖, 1)`. A Python port with σ_ref = 1 reproduced today's C++ bit for bit on 640/640 ring
increments. Because ‖σ‖ ≥ √3·p, a floor of 0.1 or 1 kPa only acts below p ≈ 0.29 kPa, and 5 kPa
only below p ≈ 1.4 kPa. To act on the ring at p' ≈ 3.5 kPa, σ_ref would have to exceed about
12 kPa.

More important, **the cost is stability-limited, not accuracy-limited**. At constant p = 2 kPa and
η/M^b = 1.00 the substep count is 60 / 54 / 52 at TolE 1e-6 / 1e-5 / 1e-4. An accuracy-limited Heun
scheme would change 10× over that range; this changes 1.15×. A 20 kPa floor saves 2 of 52 substeps
there and doubles the error. The floor is inert where the error is small and cannot reach the tail
where it is large.

Where the admitted error sits below the TIMs Newton tolerance: on every constant-p and active-path state
at p ≥ 2 kPa, today and at σ_ref = 20 (errors 1e-4 to 2e-3 kPa per 1e-5 increment). It fails on the
ring tail (p95 2.9e-2 kPa, max 0.84 kPa), and those tail numbers are understated about 2–3× because
the reference used there was itself α-blind. The conversion from the TIMs `NormUnbalance` tolerance to a
stress error at one Gauss point (≈ 2e-3 kPa if the reference load is the footing weight, ≈ 0.12 kPa
if it is the bearing load) is an order-of-magnitude element estimate: **the vector taken as "the
reference load" decides it, and only TIMs can name it (D6).**

**Status.** SHIPPED as a flag of SAS-ME only; NOT DONE in `ModifiedEuler` (it would change nothing
measurable there).

**How to use it.** `-errFloor $sigRef` on an `IntScheme 129` material; default `P_atm/101` (1 kPa at
`P_atm = 101`), i.e. exactly `ModifiedEuler`'s implicit floor. Leave it at the default.

**Evidence.** `128_sanisand_ring_trace.md` §5 (tables §5.2, §5.3, §5.5); quirks row "The
ModifiedEuler error norm already HAS a 1 kPa floor" (the ring trace); guide
`LadrunoSANISAND_implex_guide.md` §13.1.

### F18(b): rate-form stages instead of two 6×6 tangents per substep

**Asked.** Compute `dσ = C:dε − Λ·C:m` inside the substep loop, form the tangent once at the end,
match today's results to round-off.

**Answer.** Not done in `ModifiedEuler`, for two reasons. First, it cannot be round-off neutral for
`TanType 2`: the chained "consistent" tangent consumes the per-stage 6×6s
(`127_tims_2d_requests_plan.md` finding D), and that chain has its own defect (it accumulates `T`
where the recurrence needs `dT`; quirks row). Second, once the defects in §2 were found,
restructuring `ModifiedEuler` for speed would have made a wrong answer faster. SAS-ME is the
replacement: its stages compute increments at their own state and it forms **one continuum tangent
at the end state** for `TanType 1` and `2`.

**Status.** NOT DONE in `ModifiedEuler`; superseded by SAS-ME (SHIPPED). The allocation-free kernel
that would make the per-substep cost small is on the roadmap (§6).

**Evidence.** Plan finding D; quirks row "`ModifiedEuler`'s `TanType 2` 'consistent' tangent chain
accumulates `T` where the recurrence needs `dT`"; the SAS-ME row of `LEDGER_implementations.md`.

### F18(c): make `IntScheme 2` (CPPM) usable under a global Newton

**Asked.** Refuse at once instead of 2⁹ recursive halvings, a line search or better start, remove the
static work arrays, rerun F12's bearing deck.

**Answer, and an additional finding.** **The CPPM's `TanType 2` tangent had the wrong sign in
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

**Status.** SHIPPED (merged 2026-09-29, `e96f8d77d`). The static-array item is done there as
groundwork for F19 (the live shared state was `Matrix::Invert`'s scratch, now a stack-local LU;
`NewtonIter`'s statics are dead code).

**How to use it.** Recommended recipe for `IntScheme 2` under a global Newton:

```tcl
nDMaterial LadrunoSANISAND $tag <18 params> 2 2 $JacoType $TolF $TolR \
    -cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on
```

`-cppmTangent fixed` is the **default on `LadrunoSANISAND`** (owner decision); `-cppmTangent vanilla`
reproduces the old binary bit for bit. Vanilla `nDMaterial ManzariDafalias` keeps the wrong sign.
`-cppmStart explicit` is **not** in the recipe: on the review's oracle set it doubled the error on 74
of 171 increments. Use a forwarding element (the TIMs `LadrunoQuad` is one).

**Evidence.** The CPPM work-package record (tables "F18(c) refusal timing" and "F12 bearing deck, RECOMMENDED
recipe"); `origin/wp/130-sanisand-cppm-under-newton`: guide §9 "IntScheme 2 under a global Newton",
`Ladruno_files/testbed/hypo_bearing/wp130_f18c/tables_recipe.md`.

### F18(d): per-point fallback from `ModifiedEuler` to CPPM

**Asked.** When `ModifiedEuler` hits `-maxSubsteps`, hand that point's increment to CPPM; refuse only
if CPPM also fails; one-element test.

**Answer.** Built as `-meFallback cppm` (needs `IntScheme 1` and `-maxSubsteps > 0`). One-element
test: a leg that `-maxSubsteps 20` refuses at step 1 runs all 10 steps with the fallback, 20 of 20
capped updates returned by CPPM, stress within 1.3 % of the uncapped integration; where CPPM also
fails the update is refused and nothing is integrated explicitly.

On the TIMs footing, though, this lever is small: the cost is **not concentrated in the ring** (§4: the
ring holds about 5–8 % of the substeps), so a per-point fallback would recover under 10 % of the
material time. It also inherits `ModifiedEuler`'s defects for every point it does not rescue. It is
not ported to SAS-ME yet (§6).

**Status.** SHIPPED (merged 2026-09-29, `e96f8d77d`).

**Evidence.** The CPPM work-package record ("F18(d) one-element fallback"); the cost census of the footing A/B (§4).

### F18(e): whether any integrator can take the b8 ring point

**Asked.** A clear statement whether the attached b8 point (p' 0.352 kPa, η 12.87) is one no integrator can
take, and a trace of how a committed η/M^b ≈ 6 arises.

**Answer. No integrator can take it, and none should.** The dumped state is **inadmissible**: its α
lies 6.26× past the model's bounding surface measured with the Lode angle of n (ring trace), 7.3× with
α's own Lode angle (reference integrator), with `b:n = −8.19`. Driven by small probes:

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
78 rows. That is the recommended behaviour: refuse with a named code so the step controller cuts the
step, never project α back (a projection would silently rewrite history and hide the upstream defect).

**Status.** Answered (the ring trace and the reference integrator). The refusal is SHIPPED in SAS-ME.

**Evidence.** `128_sanisand_ring_trace.md` §0, §2, §4; `134_sanisand_reference_integrator.md` §0.4–0.6,
§6.4–6.6; guide §13.3 (ring row).

### F19: SANISAND in the threaded state-determination loop

**Asked.** First the inventory of shared mutable state, then per-instance / thread-local / locked
state, then identity and speed-up at 1/2/4/8 threads.

**Answer.** The inventory is done. **The prime suspect for the `IntScheme 1` segfault is not a data
race: it is a print.** Under openseespy, `opserr` goes through `PythonStream` into CPython
(`PySys_FormatStderr`). An OpenMP worker that prints (the `-maxSubsteps` cap warning was the site on
the earlier crashing threaded run) calls into CPython with no Python thread state. This explains every row of
the evidence of that run, including "a mutex and `omp critical` still crash", which rules out every
data-race explanation. The fix is a deferred per-thread message buffer flushed in element order after
the loop (which also makes the printed warnings identical at any thread count). The inventory also
found real races for step 2: the process-wide `LadrunoImplexGlobals` counters (on the plain path too,
not only under `-implex`), the warning budgets, `Matrix::Invert`'s shared scratch on the CPPM path, and
class-static return buffers in the plane-strain wrappers.

**Status.** Step 1 (inventory) SHIPPED. Step 2 (code) **PENDING**; the CPPM groundwork it
waited on is merged, and it still waits on the deferred message path. Until then SANISAND stays refused from the threaded
loop.

**Build note for Esmeralda.** A fork build change merged on 2026-09-28 flipped the CMake default
`LADRUNO_OPENMP` to ON, including gcc builds: until then a bare-cmake Linux build compiled the loop out.
An Esmeralda build from `ladruno` at or after `eeb7847d4` has the threaded loop; SANISAND still runs
serially until step 2.

**Expected gain** (inference, not measured): material update is 85.6 % of the TIMs wall time (intake §1.2),
so Amdahl bounds 8 threads at about 3.9×; dynamic scheduling and an allocation-free kernel are needed
to approach it.

**Evidence.** `131_sanisand_threaded_inventory.md` §0–§6; the OpenMP build change of 2026-09-28; intake §1.2.

### F20(a): a cumulative per-point substep counter

**Answer.** `substepStats`: 17 columns **per integration point** (none process-wide), cumulative since
`revertToStart`, **not** reset by `revertToLastCommit`, carried by `getCopy` and the wire. So a
post-mortem after a failed `analyze` reads the real history instead of the zero reported in the intake. It also counts
what used to be invisible: substeps **force-accepted at `dT_min` after failing the error test**, how
many of those fired the clamp to `Mc`, low-p abandons (the integrator returning at T < 1, silently),
and cap hits. Reading it never changes a number.

```python
s = ops.eleResponse(ele, "material", ip, "substepStats")
substeps, forced, abandoned, caps = s[2], s[5], s[8], s[9]
```

Under SAS-ME the census is `sasStats` (per point, since `revertToStart`, including the last refusal
code). Since 2026-09-29 `substepStats` has 28 columns, the extra ones for the CPPM.

**Status.** SHIPPED (counters and `sasStats`).

**Evidence.** Guide §6.2 (column table); `tests/test_ladruno_sanisand_replay_counters.py`.

### F20(b): profile scopes inside the integration

**Answer.** Added inside SAS-ME: `sanisand.sasME.predictor`, `.stateDependent`, `.stageArithmetic`,
`.stages`, `.drift`, `.alphaCheck`, `.substeps`, `.update`, `.tangent`. Measured split: about 0.44 ms
per update at a ring state and 0.025 ms at a deep state; at a ring state the stages take roughly half
to 60 % (the state-dependent quantities about 15 % of the total), drift correction about 10 %, the α
check 6–10 %; at a deep state the tangent and drift are about 11 % each.

Not added inside `ModifiedEuler` or the CPPM Newton: `ModifiedEuler` is the integrator recommended for
retirement, and its scopes would not survive the byte-identity constraint cheaply.

**Status.** SHIPPED for SAS-ME; NOT DONE for `ModifiedEuler`/CPPM.

**Evidence.** The SAS-ME work-package record ("F20(b) profile split"); guide §13.3; `Ladruno_files/testbed/wp129_sasme/out/`.

### F20(c): a `"tangentEP"` response

**Answer.** `tangentEP` returns the 6×6 continuum elastoplastic tangent at the **committed** state,
whatever `TanType` the deck uses. Checked against a one-sided finite difference at a plastic state
(50 kPa, three strain directions): relative difference below 1e-4.

```python
C = ops.eleResponse(ele, "material", ip, "tangentEP")   # 36 values, row-major
```

**Status.** SHIPPED.

**Evidence.** `tests/test_ladruno_sanisand_sasme.py::test_tangentEP_matches_finite_difference`.

### F21: replay a dumped material state

**Answer.** `ladrunoSANISANDReplay` puts a private copy of a `LadrunoSANISAND` prototype into a given
(σ, α, α_in, z, e) and drives one strain increment through the same `setTrialStrain` an element uses.
It returns rc, the census of that one update, the returned state (σ, α, α_in, z, e, p, q, f before and
after, the path code) and a per-substep trace (`T, dT, err`, outcome code). Replaying step k+1 from the
committed state of step k reproduces an analysis step to 1e-9.

**Finding A: the TIMs CSVs are compression-positive.** The README says "compression negative as
OpenSees stores it", but on all 80 rows `p_kPa = +tr(σ)/3` and every normal stress is ≥ 0: the columns
are the model's internal `mSigma`. So the replay has **no default convention**; the convention is stated explicitly:

```tcl
ladrunoSANISANDReplay $matTag -convention compressionPositive \
    -sigma s11 s22 s33 s12 s23 s31 -alpha <6 values> -alphaIn <6 values> -fabric <6 values> \
    -voidRatio $e -dStrain d11 d22 d33 g12 g23 g31 <-type 3D|PlaneStrain> <-trace 10000>
```

Shear strain is engineering (γ). α, α_in and z are projected to their deviatoric parts with a warning
(b8 row 1859/2 has tr α = 2.3e-3). A Python helper reads the TIMs CSVs and runs the standard probes:
`Ladruno_scripts/sanisand_replay.py` (`replay`, `read_ring_csv`, `probes`). For later dumps,
`substepStats` is read in the same dump.

**Status.** SHIPPED.

**Evidence.** Guide §6.3; quirks row "The TIMs ring-point CSVs carry the INTERNAL,
compression-POSITIVE `mSigma`"; `test_replay_reproduces_an_analysis_step`.

### F22: a deterministic mode

**Asked.** MKL conditional numerical reproducibility for PARDISO, a list of what else is
order-dependent, byte-identical curves twice on 8 threads.

**Answer.**

```tcl
system Pardiso -deterministic              ;# MKL CNR on the AUTO branch + iparm(34)
system Pardiso -cbwr COMPATIBLE            ;# an explicit branch every x86 node can run
```

The first solve prints what MKL has in force, e.g.
`PARDISO deterministic mode: MKL CNR branch AUTO, iparm(34)=8 thread(s), CNR ACTIVE`. Measured on a
~22k-DOF push at 8 MKL threads, 5 runs each: mode on, 1 distinct result; mode off, 5 distinct
displacement fields.

Four things to know:

- **Across nodes with different CPUs** (intake §1.6), AUTO picks a code path per CPU. Pin a branch
  every node can run: **`-cbwr COMPATIBLE`**. The instruction-set branches (`AVX2`, `AVX512` and others) exist
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
either way), and a different tolerance, step size or mesh would move it too. The TIMs deck has this
character for a measured reason: with `ModifiedEuler` the stress–strain map is non-smooth (the err = 0
path, §2), so there may be no equilibrium for Newton to converge to. On the fork's own
flip-determinism deck, the first push step has **no reachable equilibrium under any tangent**; a
10⁴× tighter TolR only halves the Newton residual floor (0.20–0.35 kN → 0.12–0.13 kN), and the
"converged" first-step load under `NormDispIncr` moves by about 25 % (9.66 → 7.24) (flip-determinism study). Use
`-deterministic` for regression tests, for reproducing a failure, and for comparing nodes; do not use
it to settle a result.

**Status.** Mode SHIPPED. The "repeatable is not reliable" guide paragraph SHIPPED.
Measured on Windows (AMD); not yet run on Esmeralda.

**Evidence.** `75c_pardiso_solver_recipe.md` Trap 7, "The deterministic mode"; quirks rows on the deterministic mode (CNR
process-wide and sticky; `mkl_cbwr_set` returning -8; `-cbwr AVX2` refused on AMD; `iparm` zeroed at
every symbolic phase); `tests/test_wp132_deterministic_pardiso.py`; `136_flip_test_drift.md`.

### F23(a): PDMY03's critical-state constants

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

**Status.** SHIPPED.

### F23(b): the PDMY "dilation brake"

**Answer.** The intake's reading is right that the brake is keyed to void ratio and that a dense sand never
reaches it with the default constants (from e = 0.6 it must dilate 17.5 % volumetrically at 100 kPa,
9.9 % at 1 652 kPa). It is incomplete in a way that matters: **reaching it would not help**.
`isCriticalState()` is a **crossing detector**. It returns 1 only for the increment whose start and
end lie on opposite sides of the line; past the line both are on the same side again and the full
dilatancy rule resumes (measured: the volumetric rate dips at the crossing step and is back within
0.5 % ten steps later). So retuning `ei`/`cs1..3`, now possible on PDMY03 too, moves *when* one
increment loses its dilatancy; **no choice of constants yields a plateau**. That is consistent with
the TIMs candidates with a retuned line failing the saturation gate as well. The route to a plateau is a
model whose dilatancy vanishes at critical state by construction (SANISAND's D ∝ M^d(ψ) − η, PM4Sand).
The TIMs ten-candidate and strip numbers were not re-run.

A related defect, found with the PDMY03 work and fixed (merged): a wild Newton iterate makes PDMY's
substep count `|Δε|/1e-5` explode to about 1e9 per call, which is why a two-element model "hung" in
`analyze`. With the fix PDMY refuses such a trial in milliseconds. Under `SSPquad` (a host that
discards the refusal) the call is bounded but the step can still be accepted.

**Status.** Note SHIPPED. Hang fix SHIPPED.

**Evidence.** `133_pdmy_notes.md` (b); quirks row "PDMY's 'dilation brake' `isCriticalState()` fires
only on the increment that CROSSES the critical-state line"; `tests/test_wp135_pdmy_substep_cap.py`.

---

## 2. What walls the deck under `ModifiedEuler`

> **Scope, since the footing A/B landed (§4).** This section explains the `ModifiedEuler` wall: s/B
> 0.0292 on the fork's copy of the TIMs deck, inside the TIMs band of 0.026–0.041. SAS-ME removes it and stops later, at
> s/B 0.0508, on a different, constitutive cause (§4.6).

### 2.1 The ring, briefly

Just outside the footing edge, the top row of Gauss points sits at p' of a few kPa with η on the
bounding surface. There SANISAND's plastic modulus scales with p and its elastic moduli with √p, so the
rate equations are stiff, and an explicit scheme takes substeps sized by stability, not accuracy. That
part is physics and would cost time under any explicit integrator. It is not what stops a
`ModifiedEuler` run.

### 2.2 What stops it: a chain of discrete defects in `ModifiedEuler`

Each link is a quirks row; the ranking is the one the independent reference integrator
established, which corrected the ring trace's first ranking of F.

| role | mechanism | what it does |
|---|---|---|
| **trigger** | **G**: α_in is re-seated once per increment | inside the substeps (α − α_in):n reaches 0, so h is the 1e10 sentinel, and then goes negative, so h < 0 and the α law becomes a repelling relaxation. 37 of 38 crossing substeps have it. The paper resets α_in at the start of each new loading process, which makes h < 0 impossible (0 of 960 runs of the exact reference). |
| **enabler** | **E**: the substep error looks at stress only | a substep that throws α 5–16× outside the bounding surface passes, because both Heun stages have the same stress increment. Adding α to the error alone keeps α inside. |
| **enabler** | **F**: a loading stage with a negative denominator is taken as elastic, and the step factor has no upper cap | both stages then agree exactly, **the error is exactly 0**, and the next substep swallows the rest of the increment. No tolerance can see it. It accounts for all 25 escapes of the TIMs `ModifiedEuler` from admissible ring starts, and for 20–65 % stress errors on benign 20–100 kPa states even at TolE 1e-8. |
| adds error | **U9**: K and G are frozen at the committed state for the whole increment | 0.6 / 6 / 24 % of the stress increment at δ = 1e-5 / 1e-4 / 1e-3. Both stages share the same wrong moduli, so the error test cannot see it and a tighter TolR does not shrink it. |
| adds error | **U10**: the loading test uses n:Δσ, not the yield-function gradient | it ignores the −(n:r)dp term, so an isotropic compression that lowers η can be read as plastic. |
| commits it | **C**: at `dT_min` a substep that failed the error test is accepted anyway, with a clamp to `Mc` | the η 12.9 → 1.33 "teleport". Uncounted before the per-point counters of F20(a). |
| commits it | `Stress_Correction` gives up silently | when neither correction direction reduces f, it returns the uncorrected state with f > 0 and rc = 0. Worst measured: f = 11.2 kPa at p = 0.58 kPa. |

The error estimate is the reason all of this stayed hidden. It measured only stress, it read zero on
the err = 0 path, and it compared two stages that shared the same frozen moduli. So every failure mode
above produced an increment the estimator called accurate. Tightening TolR, flooring the norm, or
raising `-maxSubsteps` all act on that estimator, which is why none of them moved the TIMs wall.

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
    <-sasAlphaIn reseat> <-sasErrorVars full> <-maxSubsteps $n> <-Pmin $pMin> <-Presidual $pRes>
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
- Use a **forwarding** element. `LadrunoQuad` (the TIMs element) forwards the refusal, so the step
  controller cuts the step.

The configuration run on the fork's copy of the TIMs deck (A/B arm E_B) was
`129 0 1 1e-7 1e-4 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0`, on
`LadrunoQuad -bbar` at B/8 (TolF 1e-7, TolR 1e-4). It is also the base of the integrator recommendation for the
campaign (§4.11).

**Behaviour change.** Since SAS-ME was merged, `LadrunoSANISAND::commitState` refuses to
commit a trial whose last update was refused, and this includes a **`ModifiedEuler` `-maxSubsteps` cap
hit**. Under a forwarding element nothing changes (the step already failed). Under a **discarding**
element (`SSPquad`, `stdBrick`, `BbarBrick` and others) such a deck used to commit the strain without the
stress, silently; it now fails the step and the point latches until `revertToStart`. Also, database
and restart files written by an older build will not load (the wire vector grew).

---

## 3. Cautions for existing results

1. **Every `ModifiedEuler` SANISAND result carries an integration error of this size.** Per increment:
   6–15 % of the stress increment on 1e-4 strain increments at benign 20–100 kPa states, 20–65 % on the
   err = 0 path, and U9 alone at 0.6 / 6 / 24 % for δ = 1e-5 / 1e-4 / 1e-3, none of it visible to the
   error test or shrinking with TolR. This includes the TIMs campaign curves. How much it moves a
   load–settlement curve depends on the deck; on the TIMs footing it is 5.5 % on the first step, at most
   0.64 % from s/B 0.001 to 0.0174, and then a **spurious upturn** as the `ModifiedEuler` wall
   approaches (+6.4 % over SAS-ME at s/B 0.0292, from committed states up to ρ_α 13.09; §4.2). Treat any
   `ModifiedEuler` ring-point state, any quantity read from ring points, and any `ModifiedEuler` curve
   near its wall, as unreliable.
2. **The explicit lane's failure dumps contain inadmissible states.** The b8 1950/2–3 rows are not a
   hard point of the material; they are the product of the defects in §2. Do not calibrate or test
   anything against them except a refusal.
3. **The tangent comparison of intake §1.5 had a cause.** Under `IntScheme 1`, `TanType 2`'s chained tangent
   accumulates `T` where it needs `dT`, and the stress–strain map is non-smooth (the err = 0 path).
   `TanType 0` was the only dependable choice under `ModifiedEuler` for that reason. Under SAS-ME,
   `TanType 1` and `2` are the continuum tangent. Under CPPM, the vanilla `TanType 2` had the wrong sign
   (fixed by default since 2026-09-29, `e96f8d77d`).
4. **Committed steps on the relaxed rung.** The last rung of the TIMs ladder (`KrylovNewton` at 10× the
   tolerance) is not a small print item on this deck. On the fork's copy, the share of the settlement accepted
   on that rung is 88.5 % (`ModifiedEuler`), 90.9 % (SAS-ME), 92.9 % (SAS-ME at TolR 1e-3) and 81.7 %
   (SAS-ME at B/16) (§4.4). What that acceptance does to q has not been isolated by any arm.
   Report the rung of every committed step, and the share of the settlement committed at the relaxed
   tolerance.
5. **The 30 % run-to-run shift is a signal about the deck, not the solver** (F22). Deterministic mode
   will make the two runs agree; it does not tell which is right.
6. **The `-Presidual` 1.01 / 5.05 kPa comparison of intake §1.4 was made with the defective integrator.**
   Re-measured under SAS-ME (§4.7): no residual pressure from 0.5 to 20 kPa removes the wall,
   0.5 kPa brings it earlier, and larger values reach further only through an apparent cohesion. So it
   does not justify a floor (D1).
7. **PDMY under `SSPquad`**: a refused trial is discarded by the host (the hang fix bounds the time, not the
   acceptance). Prefer a forwarding element (`quad`, `LadrunoQuad`, the u-p family) for any deck that
   depends on a material refusal.
8. **Other `ManzariDafalias` integration schemes carried defects of their own**, found and fixed in
   the fork after the fifth issue. Under `IntScheme 5` (`ForwardEuler`) a shadowed variable made every
   n:r term zero. Under `IntScheme 0, 4, 6, 7, 8, 9` a sub-stepped increment used uninitialised elastic
   moduli and a frozen elastic strain; measured before the fix, such an increment did not move the
   stress at all. Schemes 1, 2, 3, 45 and 129 are not affected by the second defect, and the campaign runs
   on 1 and 129. Any earlier result on schemes 0, 4, 5, 6, 7, 8 or 9 is to be re-run on a build that
   includes the fixes (§0a.7). Both defects are present in vanilla OpenSees.

---

## 4. The footing A/B on the campaign set

> **Scope on 3 October 2026.** This section records the A/B of `ModifiedEuler` against SAS-ME on the
> campaign set, before R1 and the low-confinement separation. Its integration findings stand. Its wall at
> s/B 0.0508 is closed (§0a.1), and its campaign-set curves are superseded for any physical reading by
> the Toyoura results of §0a.2 and §0a.3.

> **CALIBRATION CAVEAT, to be read before any curve in this section.**
>
> The element tests of the regularisation memo (§10, commit `e14703ca7`) ran the campaign SANISAND set on
> the exact reference integrator (the oracle), in drained triaxial and plane-strain compression at p0 = 10 / 50 / 150 / 500 kPa:
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
> So **none of the footing curves here is called physical** until TIMs confirm the calibration against
> their laboratory data. That covers SANISAND (966.7 kPa at s/B 0.0508, still rising), the fork's
> `DruckerPrager` control (38°, ψ = 0; max 824.2 kPa) and the TIMs `PressureDependMultiYield` control
> (PDMY01, 33° cone; limit point 417.6 kPa at s/B 0.116, intake §1.1). The request, in one
> sentence:
>
> *"The campaign SANISAND set reproduces the strength of a very dense sand but dilates ~20–23×
> (triaxial) / ~8–14× (plane strain) less than stress–dilatancy (Bolton 1986) requires, and peaks at
> 4–16 % strain (a lab-calibrated DM04 set at the same density peaks at 1–5 %). Confirmation of the
> calibration against the TIMs laboratory data (φ′_peak, strain at peak, dilatancy) is requested before
> the footing curves are used."*
>
> Everything below is about why the integration stops and what it costs. It says nothing about where a
> correctly calibrated footing would peak. §4.10 puts each curve against the classical capacity of its
> own friction angle.

Source for this section, unless stated: the footing A/B report `138_footing_sas_me_ab.md` (the A/B report; draft,
at `762be8332`), with the mechanism, the localization analysis and the capacity bands from
`150_sanisand_regularization_memo.md` (the regularisation memo; draft; §10 at `e14703ca7`, §11 at `00198f278`).

### 4.1 Setup

The fork's own copy of the TIMs deck, built from the §1 spec of the intake. Nothing in the TIMs Workbench was
run or edited. Plane-strain strip, B = 1.5 m, full width; `LadrunoQuad -bbar`; B/8 (2 430 elements,
9 720 Gauss points); `system Pardiso`; `NormUnbalance` 1e-5 × the applied vertical load (0.0415 kN);
Newton (25) → NewtonLineSearch (40) → KrylovNewton (60, tolerance × 10); ds from 2e-5 m, doubled after
6 good steps up to 1e-3 m, halved on a failed ladder, **FLOOR** when ds < 2e-7 m. Where the spec was
silent the gap was filled and declared (mesh grading reconstructed from the TIMs `ring_points_b8.csv`,
γ′ = 9.81 kN/m³, a rough guided footing, F10's step controller): A/B report §1, gaps G1–G6.

**The arms that count** ran on Esmeralda, build `ladruno` **`7936ed6e0`**, each alone on its node with
MKL and OpenMP at 1 thread, so their wall clocks compare (A/B report §8). All use the material line of the
intake with `-flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0`:

- **E_A**: `ModifiedEuler` (`IntScheme 1`), TanType 0.
- **E_B**: SAS-ME (`IntScheme 129`), TanType 0, TolR 1e-4. The reference SAS-ME arm.
- **E_D**: E_B with TolR 1e-3.
- **E_C2**: E_B with TanType 1, `-maxSubsteps 20000`, KrylovNewton at 1× the tolerance.
- **E_B16**: E_B on a **B/16** mesh (9 720 elements, 38 880 Gauss points).
- **Control**: UW `DruckerPrager` 38°, ψ = 0, on the same deck (a local run).

The earlier local legs (A and B, to s/B 0.0174) are kept for the early-curve comparison and the replay
study (§4.4, §4.5). E_B
reproduced the local B leg to 1e-5 kPa through step 40 (A/B report §4).

### 4.2 Where each arm stops

| arm | integrator | s/B at FLOOR | q (kPa) | first `loadingNonPosH` at s/B | refusals (converged-step census) | push wall (h) |
|---|---|---|---|---|---|---|
| E_A | `ModifiedEuler` | **0.0292** | 701.8 | none (no refusal path) | 542 cap hits, 3 943 forced at dT_min | 2.90 |
| E_B | SAS-ME, TolR 1e-4 | **0.0508** | 966.7 | **0.0363** | NonPosH 232, maxSubsteps 209, errorAtDTmin 1 | 4.32 |
| E_D | SAS-ME, TolR 1e-3 | 0.0410 | 808.3 | 0.0334 | maxSubsteps 13 316, errorAtDTmin 190, NonPosH 173 | 4.67 |
| E_C2 | SAS-ME, TanType 1 | 0.0114 | 317.0 | 0.0064 | errorAtDTmin 175 011, maxSubsteps 40 683, NonPosH 409 | 12.35 |
| E_B16 | SAS-ME, B/16 | 0.0135 | 352.8 | 0.0135 | maxSubsteps 115, NonPosH 11 | 3.75 |
| control | `DruckerPrager` 38°, ψ = 0 | 0.15 (target reached) | 752.0 (max 824.2) | — | — | 0.07 |

(A/B report §0 and §8. The census sums the per-step refusal lines of the converged steps; the final ladder at
the floor is not in it.)

**Every SANISAND arm stops on the step floor.** The SAS-ME arms reach it on `loadingNonPosH` refusals.
**E_A cannot refuse**: `ModifiedEuler` has no refusal path, so it floors through cap hits and
acceptances forced at dT_min (542 and 3 943 over the run), and it **commits** what it cannot
integrate. Over s/B 0.026–0.0293 it forced 3 786 acceptances and the committed ρ_α reached **13.09**;
SAS-ME over the same window stays at ρ_α ≤ 1.004 with no forced acceptance. The ModifiedEuler curve
bends up there, to **+6.4 %** over E_B at the same s/B, and that stiffening is spurious (A/B report §8.4).
E_A's wall sits inside the TIMs `ModifiedEuler` band (s/B 0.026–0.041), so the fork's copy reproduces
the TIMs wall. `ModifiedEuler` curves near their wall are read with this in mind.

### 4.3 The verdict

**SAS-ME moves the wall from s/B 0.0292 to 0.0508, but it does not remove it.** There is **no peak and
no plateau** on this deck: at the wall q = 966.7 kPa and still rising (q_max = q_end; the slope over the
last 0.005 s/B is 0.24× the initial slope). The first `loadingNonPosH` refusal comes at s/B 0.0363, and
from there the refusals accumulate until they end the run (A/B report §0, §8.1, §8.3).

**The cause is constitutive, not integration.** No integrator setting lifts it (TolR 1e-3 walls earlier,
TanType 1 much earlier, §4.4), the independent oracle stops at the same kind of state, and the
refusing points are pre-peak (ρ_α < 1). §4.6 gives the mechanism. This revises the first issue of this
report, which called the wall an integration failure: that is true of the `ModifiedEuler` wall
(0.0292; §2), not of the one SAS-ME reaches.

### 4.4 Accuracy and cost

**Per-increment accuracy** (replays of real increments against the oracle, §4.5; A/B report §5.3). SAS-ME
sits at **0.6–2e-4·p′** at every checkpoint: that is its TolR 1e-4 error control. `ModifiedEuler`'s
error depends on the state: it is **about 20× worse only at the onset of ring plasticity** (s/B ≈ 0.001),
and **at parity from s/B ≈ 0.01**, apart from isolated low-p outliers (2e-3·p′, a point taken in one
substep). It is not a blanket accuracy factor.

**The load–settlement curves**, `ModifiedEuler` against SAS-ME on the local legs (A/B report §5.1):

| s/B range | max \|q_B − q_A\| / q_A |
|---|---|
| 0–0.001 (the worst is the first step, s/B 1.3e-5) | 5.5 % |
| 0.001–0.005 | 0.64 % |
| 0.005–0.0174 | 0.51 % |

Past s/B 0.016 the Esmeralda pair separates, to 6.1 % over 0.016–0.029, as E_A turns up (§4.2).

**Wall clock and substeps per 0.01 s/B** (A/B report §8.2; hours / 1e9 substeps, partial intervals prorated):

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
7–9 % (A/B report §6). The cost is spread over the whole plastic zone, so a per-point fallback for the worst
points (F18(d)) would recover under 10 %. The multiplier is the global iteration count × every point's
update.

**Most of the settlement is accepted on the relaxed rung.** The KrylovNewton rung accepts at 10× the
test tolerance (0.415 kN against 0.0415 kN). The share of the settlement accepted there is **88.5 %**
(E_A), **90.9 %** (E_B), **92.9 %** (E_D) and **81.7 %** (E_B16) (A/B report §8.2). The effect of that
acceptance on q is **not isolated by any arm**; E_C2 changed the Krylov tolerance together with two
other settings.

**TolR 1e-3 is not a lever (E_D).** It walls **earlier**, at s/B 0.0410 against 0.0508. It costs the
same as E_B per unit s/B up to s/B 0.038 and then rises to 4.6× E_B's (0.038–0.041); it has 13 316 maxSubsteps
refusals against 209. Its q runs **−3.7 %** (median) below E_B over s/B 0.02–0.041 (range −4.3 % to
−1.9 %), and −1.4 % (median) over 0.001–0.02. No saving, an earlier wall: TolR 1e-4 is kept (A/B report §8.5).

**TanType 1 is not viable here (E_C2).** Global iterations per step fall (median 11 → 4 → 2), but the
accepted step collapses with them (median ds 1.25e-6 m past s/B 0.01, where E_B runs at 1e-3 m). E_C2
floors at **s/B 0.0114** after 12.35 h, 13× E_B's cost per unit s/B, while its curve stays within
1.45 % of E_B. The consistent tangent buys a step-size collapse, not settlement (A/B report §8.6). This
supersedes the stopped local TanType 1 arm of the first issue.

### 4.5 The replay figures, reconciled

The first issue quoted two pairs of figures for the per-increment error on the footing's real
increments and asked which increments each covers. Both are right, for different increments (A/B report
§5.3). Each row is a Gauss point's committed state plus the strain increment it received in
the next converged step, replayed through `ModifiedEuler`, SAS-ME and the reference integrator (the oracle; Radau, rtol
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

DM04's plastic modulus has a singularity at an α_in re-seat (regularisation memo §1.4; A/B report §10):

```
a = (α − α_in):n,   h = b0 / a,   Kp = ⅔·p·h·(b:n)
```

At a re-seat α_in := α, so a = 0 and h is infinite (capped at 1e10 in the code). If b:n > 0 the
modulus is large and positive; if b:n ≤ 0 it goes to −∞ and the increment has no solution. The
singular set is **{a = 0, b:n ≤ 0}**, and `loadingNonPosH` is the SAS-ME refusal that names it.

- **The refusing points are pre-peak** (ρ_α 0.93–0.96, dense, ψ ≈ −0.1, inside the bounding surface) and
  they chatter: E_B makes 10.1 million re-seats (regularisation memo §1.3).
- **It is in the continuum equations, not in SAS-ME.** The exact reference integrator stopped on the same 0/0 at
  the TIMs ring point 1950/3 (reference-integrator report; regularisation memo §1.4).
- **No integrator knob lifts it** (§4.4), and **no boundary-value regularizer** can lift an unbounded
  negative modulus at a point. Any cure is a change to the model (R1, §5).

**How the equations reach it (R1 memo §2.2).** The exact rate
equations reach this set through a **Zeno accumulation of re-seats on the b:n → 0⁺ side**: after a
re-seat, h = ∞ makes α slide along b; with b nearly perpendicular to n that slide rotates n, a turns
negative and the next re-seat follows. At E_B's refuser 1880/1 the re-seat intervals run 1.7e-2,
1.5e-3, 5.6e-5, 2.4e-6, and so on, and accumulate at a finite time where a = 0 and b:n = 1.9e-8 > 0.
**`loadingNonPosH` is only the b:n < 0 exit of that sequence.**

**What makes the wall states singular: the non-convex extension side** (R1 memo §2.5; regularisation memo §2.3).
All three wall refusers have n on the extension side (cos 3θ −1.00 / −0.36 / −0.88), while the
non-elliptic band points of §4.8 sit on the compression side (0.09 % / 0.14 % of them extension-side at
the E_B / E_B16 walls): the wall and the bands are separate phenomena. With c = 0.71 < 7/9 the Lode
interpolation is concave on the extension side (§4.9). Driven at c = 0.80, the same five committed wall
states fail **0 of 320** exact trials against **102 of 320** at c = 0.71. These are c = 0.71 states
driven at c = 0.80: a sensitivity test, not a c = 0.80 footing run. The footing run with c = 0.80 and R1 off
walls at s/B 0.048, on the compression side (`bd93c558d`, 2026-09-29; §0a.1), so raising c only delays the wall.

**Two routes out of the wall; the choice belongs to TIMs (D8, D9).**
1. **R1**: an opt-in model-level fix at any c, no recalibration (§5, D9).
2. *(Withdrawn 2026-09-29: at footing scale this only delays the wall; see §0a.1.)* **A calibration with c ≥ 0.78**: at c = 0.80 the extension strength M_e = c·M_c rises 13 %; it also
   removes the extension ill-conditioning of §4.9.

**Related literature.** Stress overshooting of bounding-surface models at load reversals, where the
reversal memory is re-seated, is a documented source of numerical instability in boundary-value
problems: Chen, Ghorbani, Zhang & Kodikara (2022), "Stress overshooting solution for soil plasticity
models", *Comput. Geotech.* 152, 105008, which finds that the definition of the plastic modulus and
hardening law governs it; and Ghorbani, Chen, Kodikara, Carter & McCartney (2023), "Memory repositioning
in soil plasticity models used in contact problems", *Comput. Mech.* 71, 385–408. No published
SANISAND footing that stops on this exact singular set has been verified.

### 4.7 Sensitivity ladders, final

Each leg is E_B with one knob changed (the S1 → S4 ablation is cumulative); all legs ran to their end
(A/B report §11, records under `Ladruno_files/testbed/footing_sas_me_ab/ladders_final/`). The S and A0/h0
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
  (regularisation memo §1.4 at `2a82e2046`).
- **Presidual 0.5–20 kPa never clears it**, and the onset is **non-monotonic**: 0.5 kPa brings it
  *earlier* (0.0182) than Presidual 0. A larger Presidual walls later (up to 0.0964) only by stiffening
  the response (1 979 kPa at the wall for 20 kPa): an apparent cohesion, not a cure.
- **e_init 0.65–0.85 never clears it.** A0 is non-monotonic (onset 0.0237 / 0.0363 / 0.0310 at A0
  0.02 / 0.05 / 0.10). **h0 × 3 brings the onset down to 0.0091.**
- **What no leg changed is the Lode ratio c** (§4.6, §4.9). Every leg keeps c = 0.71 < 7/9. The footing
  test of that reading (c = 0.80, R1 off) walls at s/B 0.048 (`bd93c558d`, 2026-09-29); R1 is the route
  out (§0a.1).

**The `-Presidual` comparison of intake §1.4 (caution 6), re-measured under SAS-ME.** The intake's conclusion holds in
the sense that matters: no residual pressure from 0.5 to 20 kPa removes the refusal. It does move where
the run stops, in both directions: 0.5 kPa floors earlier (s/B 0.0303 against 0.0508), and the larger
values reach further only by adding strength that is not in the sand. Do not use `-Presidual` to push
past the wall (D1).

### 4.8 Mesh: B/16

E_B16 floors at s/B 0.0135, earlier than B/8 (A/B report §9). Up to there, B/16 runs softer than B/8 from
s/B 0.010:

| s/B | 0.002 | 0.005 | 0.008 | 0.010 | 0.012 | 0.013 | 0.0135 |
|---|---|---|---|---|---|---|---|
| q_B16 / q_B8 − 1 | −0.58 % | −0.14 % | −0.91 % | −3.89 % | −4.95 % | −4.84 % | −3.57 % |

**The band is one element wide on both meshes** (full width at half maximum of the incremental shear
strain 1.00–1.11 element sizes), it halves with the element and it follows the mesh lines.

**What it is: non-associated localization that begins while the material is still hardening** (regularisation memo
§2.1, Rudnicki & Rice 1975). A plane-strain acoustic-tensor scan of the continuum tangent finds
det ≤ 0 at **16.8 % of the Gauss points at s/B 0.011 on B/8** (9.1 % at s/B 0.0096 on B/16, 16.9 % at
its wall), where H/2G ≈ 1.05–1.07. The same states with associated flow are **elliptic everywhere**; only 17
of 9 720 points are post-peak even at s/B 0.0508. The cause is the flow rule: friction of 45–60° against
dilation of 1–2° (§4 caveat), so the bands are partly a product of the calibration.

**It is not ψ-softening**, so a nonlocal ψ̄ (void-ratio averaging) or a crack band would not treat it.
The wall (§4.6) is a separate matter: the floor refusers are pre-peak on both meshes, and the mesh
dependence of the wall itself was not settled; R1 has since closed the wall (§0a.1).

### 4.9 A second calibration item: Lode convexity (c ≥ 7/9)

With **c = 0.71 < 7/9**, DM04's Lode interpolation g(θ) is **non-convex at the extension meridian**.

- **What was seen** (R1 memo §6.3, harness and results committed under
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

A check on the FE, not on the physics (regularisation memo §11, commit `00198f278`). The classical rough-strip
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
| TIMs PDMY01, 33° | 417.6 kPa at s/B 0.116 | 381 kPa | +10 % |
| fork `DruckerPrager`, 38°, ψ = 0 | a plateau of ~700–820 kPa over s/B 0.06–0.15 | 812 kPa | at or below the associated value, as expected for ψ < φ |
| SANISAND, campaign set | 967 kPa at s/B 0.05, still rising | its own element φ′_ps,peak at the footing's p′ ≈ 50–150 kPa is 51–55° (T5), so q_u ≥ 2 753 kPa | about ⅓ of its classical capacity mobilized, which is what the late element peak (4–16 %) predicts |

**The three FE curves are each consistent with their own constitutive strength** (PDMY01 33° ≈ 381 kPa
classical, DP 38° ≈ 812 kPa, SANISAND's own 51–55° ⇒ ≥ 2.75 MPa). **The physical capacity follows from
the target sand's φ′**, and those inputs are owed (§5.1). With the laboratory φ′, the physical band is
[q_u(φ′_cs,ps), q_u(φ′_ps,peak at the footing's mean p′)]: the operative angle lies between the
critical-state and the peak angle because of progressive failure and stress level (Lau & Bolton 2011;
Perkins & Madson 2000; Loukidis & Salgado 2011). "Which curve is physical" is the same question as
"what is the target sand's operative φ′".

### 4.11 The integrator recommendation (a decision of the owner and of TIMs)

**SAS-ME (`IntScheme 129`), TanType 0, TolR 1e-4, `-maxSubsteps 2000`, with R1 and the low-confinement
separation, on the step policy of the runs of §0a.** The material line for the Toyoura legs:

```tcl
nDMaterial LadrunoSANISAND $tag <the 18 constants of the set> \
    129 0 1 1e-7 1e-4 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0 \
    -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 \
    -sasTensionCutoff 0.1 1.0
```

- `-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5` is R1 with κ = 0.5 (merged).
- `-sasTensionCutoff p_sep p_contact` is the low-confinement separation, in model stress units (kPa
  here). p_sep = 0.1 kPa is the smallest value tested and moves the peak by 0.09 % against 0.5 kPa;
  p_contact = 1.0 kPa is the value of the runs of §0a. The flag is on `ladruno` from `fc2ea4fc8`
  (§0a.7) and requires `-Presidual 0`.
- The campaign constants of the first issues are not recommended for monotonic footing work (§0a.6).

Step policy: ds0 = 2e-5 m, × 2 after 6 good steps up to 1e-3 m, ÷ 2 on a failed ladder, floor at
ds < 2e-7 m; Newton (25) → NewtonLineSearch (40) → KrylovNewton (60, tolerance × 10); `NormUnbalance`
1e-5 × the applied vertical load; 1 MKL thread with `MKL_CBWR=COMPATIBLE`, or `system Pardiso
-deterministic` on current builds.

Reasons: on the campaign set SAS-ME reaches 1.74 times the settlement of `ModifiedEuler` at 1.3 to 1.45
times lower cost per unit s/B, its per-increment error is controlled, and it refuses a state it cannot
integrate where `ModifiedEuler` commits ρ_α 13 and a spurious +6.4 % (§4.2, §4.4). R1 and the separation
close the two later limits (§0a.1). Rejected: TolR 1e-3 (earlier wall, no saving) and TanType 1 (step
collapse, 13 times the cost).

**The recommendation concerns the integrator and the numerical limit only.** It does not make a
capacity mesh-converged (§0a.3), it says nothing about the calibration (§0a.6), and most of the
settlement is still accepted on the relaxed KrylovNewton rung (caution 4).

---

## 5. Decisions for TIMs

The table gives recommendations. The decisions belong to the calibration and to the project.

| # | decision | recommendation | reason |
|---|---|---|---|
| D1 | **The p′-floor rule** | (a) `-Pmin` ≤ 0.5 kPa; (b) `-Presidual 0`, which the low-confinement separation requires; (c) the limit load reported at floor F and at F/2, with the floor accepted if the load moves by less than about 2 %; (d) the number of Gauss points at the floor and the number separated, reported at the limit state | under SAS-ME a residual pressure of 0.5 to 20 kPa never removed the wall; 0.5 kPa walls earlier, and larger values reach further only through an apparent cohesion (§4.7). On the Kimura case at s/B 0.045, p_r 2 / 5 / 10 kPa give 820 / 839 / 861 kPa (`8ebde5cbd`, 2026-09-29); the separation agrees with the p_r → 0 extrapolation within 2.3 % |
| D2 | **`D_factor`, the UW low-pressure dilatancy sigmoid** | decided explicitly, not inherited | it is not in Dafalias & Manzari (2004); it acts below p′ = 5.05 kPa and only suppresses dilation. It is not the cause of the over-dilation of §0a.5. No deck switch exists; an opt-in switch is added on request |
| D3 | **Mesh** | B/8, B/16 and a B/8 sheared 15° reported before any curve is called mesh-converged; no capacity of calibration grade before the Perzyna regularisation; acceptance at half the test scatter, about ±4 % from Kimura's N_γ scatter | on Toyoura B/16 is about 20 % below the B/8 peak at s/B 0.20, with no peak; on the Kimura case the B/16 leg carries 1757 kPa at s/B 0.183 against the 2296 kPa B/8 peak, still rising (§0a.3); the bands follow the mesh lines |
| D4 | **The definition of the limit load for a dense dilatant sand** | stated before the runs: the peak, a plateau, or q at a fixed s/B | on Toyoura the Kimura case peaks at s/B 0.165 to 0.18; on the TIMs geometry there is no peak by s/B 0.20 (§0a.2). Where no clear peak forms, Vesić's rule takes q at s/B 0.10 |
| D5 | **Re-running campaign curves** | the curves that feed a reported number re-run under SAS-ME with R1 and the separation; also any curve computed with `IntScheme` 0, 4, 5, 6, 7, 8 or 9 | §3, cautions 1 and 8; §4.2 (+6.4 % spurious stiffening near the `ModifiedEuler` wall) |
| D6 | **The reference load for `NormUnbalance`** | the vector named | it decides whether the integrator's per-point error sits under the Newton tolerance (F18(a)) |
| D7 | **The calibration** | the campaign set replaced for monotonic footing work; PB2 or DM04 Toyoura as interim references; the final set calibrated on the target sand's laboratory data (§5.1) | the campaign set is a cyclic fit applied to a monotonic problem (§0a.6) |
| D8 | **Lode parameter c** | c ≥ 7/9 (≈ 0.78), for Lode convexity only | at c = 0.71 the Lode interpolation is non-convex at the extension meridian, and a round-off perturbation decides a CTXu result (§4.9). At footing scale c = 0.80 without R1 walls at s/B 0.048, so c is not a cure for the wall |
| D9 | **R1** | used for the campaign, with κ = 0.5 | it closes the SAS-ME wall at any c without recalibration (below) |
| D10 | **The separated zone** | the low-confinement separation accepted as the idealisation of the heaving wedge, with p_sep = 0.1 kPa | the peak is insensitive to p_sep (0.09 %) and p_contact (0.35 %), and within 2.3 % of the p_r → 0 extrapolation; the zone is a dilatant heave driven by the model (§0a.1). Whether that dilation is physical is §0a.5 |
| D11 | **Acceptance of the initial stiffness** | a tolerance on the initial secant against the reference test, stated before the stiffness work; ±20 % proposed | the stiffness is a deliverable for the TIM macroelement, and it is 0.44 of the test today (§0a.4) |

**D9 in detail: R1.** Opt-in flags of `LadrunoSANISAND` with `IntScheme 129`, all default off and
byte-identical when off; vanilla `ManzariDafalias` is unchanged. Two coupled parts:

- **a floor on h everywhere**: h = b0 / max(a, c_A·√(2/3)·m), with c_A ≈ 1 (the yield cone's α-space
  radius, 4.1e-3 at m = 0.005), used in the α update too;
- **a hysteretic re-seat**: α_in := α only when a < −c_rev·√(2/3)·m (c_rev ∈ {½, 1, 2} all pass).

On the material-point oracle, 320 exact increments from 5 real refuser states, the pair fails **0 of 320**
where DM04 fails **102 of 320**, and **each part alone fails** (the floor alone 97; a weaker floor,
c_A = ¼, with the hysteresis 90). In element tests at c_A = 1 monotonic q moves by at most
**2.7e-4·q_max**, and drained cycles are unchanged to 4 digits. The CTXu difference that was open is DM04's
own extension bifurcation (§4.9), identical with and without R1. R1 is a constitutive change: it
changes the model being calibrated, which is why the decision belongs to TIMs.

**Regularisation.** The mesh study answers whether it is needed: it is (§0a.3). The candidate is
**Perzyna-type viscoplasticity inside `LadrunoSANISAND`**, not a Duvaut–Lions wrapper: a Duvaut–Lions
update needs the inviscid solution for the same increment, which is exactly the computation that
refuses. **A shear-band width is not a meaningful deliverable** here: the physical band is about
20·d50 ≈ 4 to 10 mm for d50 0.2 to 0.5 mm, against elements of 94 to 188 mm, so any width a regularised
mesh returns is a numerical length, not the sand's.

### 5.1 Inputs owed by TIMs

1. **The sand**: d50 and the grading, e_max / e_min, and the target relative density D_r.
2. **Laboratory data**: drained triaxial and plane-strain tests with the volumetric strain (φ′_peak, the
   strain at peak, the dilatancy and the post-peak drop), at confining pressures down to the
   near-surface range of the footing (5 to 50 kPa) where available; the critical-state void ratio at low
   pressure where available; any undrained cyclic CSR–N target.
3. **The footing test treated as the reference.**
4. **For that reference test:** the footing roughness; how the sand was placed relative to the load
   direction (pluviation or bedding); the measured unit weight and the e_max / e_min of that batch. On
   Kimura et al. (1985) bedding alone moves the peak by 11 % and its settlement by 1.8 times, and
   roughness can move the peak by 25 to 35 % (§0a.2).
5. **Small-strain stiffness of the target sand** (new in this issue; needed for the initial stiffness):
   - G_max at several confining pressures, or a shear-wave velocity profile (bender element, resonant
     column or field Vs), with the confining pressure or depth of each measurement; the pressure
     exponent n follows from these data;
   - a G/G_max–γ curve (resonant column or torsional shear);
   - the stiffness definition the TIM calibration uses: the initial tangent of the footing curve, a
     secant at a stated s/B, or G_max directly.
6. **The exact PDMY01 33° parameter set** behind the 417.6 kPa control. The fork holds only its own
   PDMY03 stand-in (φ 40°), which has no peak in drained plane strain and is not a physical reference.

---

## 6. Roadmap

**Near term.**

| item | what | state on 2026-10-03 |
|---|---|---|
| DR at the peak | DR against implicit at the peak (implicit 2293.4 kPa at s/B 0.1695, `2b5cb2f8f`) | running; expected 3 to 4 October |
| Mechanism | Perzyna viscoplasticity, opt-in, oracle first; the gate matrix includes a sheared mesh and a row-orientation variant | plan; tuned after the two constitutive tracks |
| Initial stiffness | G_max decay with exponent n, opt-in, oracle first; checks against an elastic strip on a modulus increasing with depth and against Kimura et al. (1985) | plan |
| Post-peak dilatancy | critical-state recalibration, element first; then a finite-element biaxial specimen; then the footing; a state-dependent option only if needed | plan, partial; direction pending |
| Post-peak band base | dump pairs over s/B 0.17 to 0.20 to measure the depth of the flat band against the mesh rows | not launched |
| F19 step 2 | the deferred message buffer, the counters, stack scratch; `IntScheme 1`, then 2 and 129; identity and speed-up at 1 / 2 / 4 / 8 threads on a deck that prints | pending |
| Refusal-aware line search | a material refusal treated as "backtrack", not "rung failed" (119 LineSearch failures on a refusal in E_B) | follow-up |
| SAS-ME follow-ups | the maxSubsteps refusals on the surface ring (209 in E_B) and the re-seat chatter; from its confirmation review: a 64-sample unload-then-reload locator can miss a very short elastic excursion (up to 5.5e-4 on 5 ring cases), the ψ-driven dead end is relocated to about 30 MPa, not removed, positional TolF/TolR are not range-checked | follow-up |

**Order of the work.** The two constitutive tracks, the initial stiffness and the post-peak dilatancy,
are decided at the element level first. The Perzyna regularisation is tuned after them, since it depends
on the hardening it regularises. The footing is then run on B/8, B/16 and a sheared B/8, and compared
with real footings:

- capacity: the band [q_u(φ′_cs,ps), q_u(φ′_ps,peak)] from the target sand's laboratory φ′ (§4.10);
- settlement: dense sand in general shear develops the full mechanism at s/B ≈ 6 to 8 % (the JGGE
  closure, doi:10.1061/JGGEFK.GTENG-12726); where no clear peak forms, Vesić's rule takes q at
  s/B = 10 %;
- failure mode: a wedge plus a radial shear zone reaching the surface, against the present bands
  captured by the mesh;
- initial stiffness: the secant against the reference test, within the tolerance of D11.

The tolerance on capacity (the scatter of comparable footing tests) has not been extracted beyond
Kimura's N_γ scatter; no further number is claimed for it.

**SAS-ME under IMPL-EX: conditional, not qualified.** SAS-ME under the IMPL-EX companion (opt-in
`-implexAllowScheme129`) is conditional and not yet qualified. The tangent identity holds, and IMPL-EX
reduces SAS-ME's substep count; its remaining cost is substeps per companion return at tight TolR. SAS-ME
refuses more often on adversarial states, the IMPL-EX companion has no global Newton to absorb a refusal,
the qualification gates have not been run for `IntScheme 129`, and there is no softening leg. The opt-in
comes from a study still in progress and is not on `ladruno`: there, `-implex` is still refused with
`IntScheme 129` (§2.3).

**Performance.** Material cost is (global iterations) × (cost per point per iteration), everywhere in
the plastic zone: every OpenSees iteration re-integrates the whole increment from the committed state.
The levers, in order:

1. **Fewer global iterations**: a step-size policy that targets a few iterations per step. The
   consistent tangent is not that lever on this deck (TanType 1 collapses the step, §4.4).
2. **Threading** once F19 step 2 lands: bounded near 3.9 times at 8 threads by the 85.6 % material share
   of the intake (inference).
3. **An allocation-free kernel**: SAS-ME builds dozens of heap vectors per substep; a fixed-size stack
   kernel is typically several times faster for 6-vector arithmetic. Unmeasured; it is benchmarked
   before a number is quoted, and it is also a precondition for threads to scale.
4. **A stiff-point fallback** SAS-ME → CPPM: a smaller lever on this deck than it first looked (§4.4).

TolR 1e-3 is off the list: on the footing it gave no saving and an earlier wall (§4.4). GPU or SIMD
batching and surrogate models are not levers for this problem: the per-point work is very unbalanced
(tens to 10⁵ substeps per point per step), and a validation study does not trade the full-order answer
for speed.

**Longer term: a model consistent by construction.** DM04 SANISAND has no energy function for α, so
non-negative dissipation is not guaranteed, and the α_in memory and the h ∝ 1/((α − α_in):n) singularity
are where the defects and the wall lived. Two classes are planned, with tags reserved and no code yet:

- **`LadrunoNORSAND`**: NorSand in the Borja & Andrade (2006) form, with DM04's power-law critical-state
  line (so ψ stays bounded as p′ → 0) and a Lode-angle dependence for plane strain; hyperelastic,
  implicit, closed-form tangent, proven non-negative dissipation. About 5 to 6 engineer-weeks. It is the
  route if the comparisons with real footings show that the model, not the calibration or the mesh, is
  the problem.
- **`LadrunoHySAND`**: the 2026 hyperplastic multisurface sand model, with the most rigorous
  thermodynamics and the best cyclic behaviour, the least evidence and no public code; a research build
  of about 10 to 14 weeks, for the later cyclic SSI work.

No dilatant sand model is truly variational. Non-associated dilatancy gives a non-symmetric tangent and
no incremental minimum principle. "Thermodynamically admissible, robust implicit return, consistent
tangent" is achievable; "variational" is not. A p′-floor or a surcharge is established practice even for
the rigorous models, so D1 remains.

---

## 7. References

**Intake and plan**
- `Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md` + attachments
  `_tims_2d_model_requests_2026-09-25/` (README, `ring_points_b8.csv`, `ring_points_b16.csv`)
- `Ladruno_implementation/127_tims_2d_requests_plan.md` (findings A–D, work packages)

**Evidence documents** (fork, `Ladruno_implementation/`)
- `128_sanisand_ring_trace.md`: the ring trace (with its 2026-09-27 correction note after the reference
  integrator)
- `134_sanisand_reference_integrator.md`: the reference integrator, the oracle; U1–U10
- `131_sanisand_threaded_inventory.md`: the threaded-state inventory (F19 step 1), E9 re-verified after
  SAS-ME
- `133_pdmy_notes.md`: the PDMY notes
- `136_flip_test_drift.md`: the flip-determinism study
- `138_footing_sas_me_ab.md`: the A/B report (draft, at `1f22e2bad`): §0 verdict, §5.1 curves, §5.3
  replays, §8 Esmeralda arms, §9 B/16, §10 diagnosis, §11 ladders (final), §12 default, §13 follow-ups;
  run records under `Ladruno_files/testbed/footing_sas_me_ab/`
- `151_sanisand_reseat_singularity.md`: the R1 memo (at `6e9330a3d`, build `bd93c558d`): §2.2 the Zeno
  re-seat accumulation, §2.5 the non-convex extension side and the c = 0.80 wall-fan test, §5 the fan,
  §6.1–6.2 calibrated behaviour, §6.3 the CTXu gate, §8 recommendation, §9 the flags and C++ gates
- `150_sanisand_regularization_memo.md`: the regularisation memo (draft, at `2a82e2046`, §1.4 corrected
  with the final ladders; the T5 figures as of `e14703ca7`): §1.4 mechanism and the R1 oracle box, §2
  localization, §3 options, §4 staged recommendation, §8 test plan and the mesh-perturbation patch, §9.1
  decision procedure, §10 T5, §11 T6 capacity bands
- `_sand_model_survey_2026-09-27.md`, `_sanisand_external_survey_2026-09-27.md`,
  `144_ladruno_norsand_plan.md`, `145_ladruno_hysand_plan.md`: the sand-model surveys and the two
  longer-term model plans

**Guides**
- `LadrunoSANISAND_implex_guide.md` §6.2 (`substepStats`), §6.3 (replay), §13 (choosing an IntScheme;
  SAS-ME); §9 "IntScheme 2 under a global Newton"
- `75c_pardiso_solver_recipe.md` Trap 7 and "The deterministic mode", with the follow-up paragraph

**Ledger rows** (`LEDGER_quirks.md`): the ring CSV convention (finding A); force-accept at `dT_min`
(finding C); F (negative denominator as elastic); the uncapped step factor; G (α_in once per increment);
`ModifiedEuler` TanType-2 chain (T vs dT); RK45 `dAlpha3/4`; IntScheme 4 non-determinism; U9; U10; the
err = 0 path; stress-only error (E); `Stress_Correction`'s silent give-up; RK45 is IntScheme 45 and not a
reference; the 1 kPa floor and stability-limited cost; the flip-determinism pins; the CNR rows of the deterministic mode;
PDMY03 constants and reallocation; the PDMY crossing detector; static `pAtm`; a refusal under a
discarding element was committed. `LEDGER_implementations.md`: the rows of the counters and replay,
SAS-ME, the deterministic mode, the PDMY03 work and the reference integrator.

**Change requests.** The state on 2026-10-03 is in §0a.7, the only list of change requests in this
reply.

**Workbench sources of the sixth issue.** The 2D-model act handoff, §§15–19 (Esmeralda reads of
2026-09-30 to 2026-10-02), and the meeting pages of 1 October 2026 (the Kimura comparison and the three
coherence tracks), from which the wording of §0a.2 to §0a.5 is taken.

**Literature cited in §4–§6.** Dafalias & Manzari (2004), *J. Eng. Mech.* 130(6); Bolton (1986),
*Géotechnique* 36(1); Rudnicki & Rice (1975), *JMPS* 23; Vesić (1973), *JSMFD* 99(SM1); Perkins &
Madson (2000), *JGGE* 126(6); Loukidis & Salgado (2011), *Géotechnique* 61(2); Lau & Bolton (2011),
*Géotechnique* 61(8); Tatsuoka et al. (1986), drained plane-strain tests on Toyoura sand at low
pressure; Tatsuoka et al. (1991), ASCE GSP 27; Kimura, Kusakabe & Saitoh (1985), *Géotechnique* 35(1),
33–45, doi:10.1680/geot.1985.35.1.33; Verdugo & Ishihara (1996), *Soils Found.* 36(2), 81–91; Gibson;
Booker et al. (1985); Gazetas (1991); Martin (2005), *Proc. 11th IACMAG*; Han et
al. (2016), *SpringerPlus* 5, 1482. Full list in the regularisation memo.

**Earlier replies to the TIMs team.** `86_ladruno_sanisand_tims_report.md` (the hidden cohesion),
`90_ladruno_regularization_tims_report.md` (regularization, `-maxSubsteps`).

---

## Revision log

| date | change |
|---|---|
| 2026-09-28 | First issue. F18(a), (b), (e), F20, F21, F22 (mode), F23(a), (b) answered from merged work; F18(c), (d), F19 step 2, the F22 guide paragraph and the footing A/B pending; placeholders in §4. |
| 2026-09-28 | Second issue. §4 final from the Esmeralda arms of the footing A/B: the wall table, the verdict (SAS-ME moves the wall from s/B 0.0292 to 0.0508 and does not remove it; no peak; constitutive), accuracy and cost, the replay figures reconciled, the mechanism, the interim sensitivity ladders with caution 6 re-measured, B/16 and non-associated localization, the calibration caveat, Lode convexity c ≥ 7/9, the classical capacity bands, the integrator recommendation. §5: D1 and D3 revised, D7–D9 added, the regularization route and §5.1. §6: the decision procedure, the SAS-ME + IMPL-EX status. |
| 2026-09-28 | Third issue. §4.7 ladders final: every leg walls on `loadingNonPosH`; dilatancy off only delays the onset. §4.6: the wall states need the non-convex extension side; the wall and the bands are separate phenomena; two routes out (R1 or c ≥ 0.78); related literature. |
| 2026-09-28 | Fourth issue. §4.6: the Chen et al. (2022) citation corrected; Ghorbani et al. (2023) added. |
| 2026-09-29 to 2026-09-30 | Fifth issue and addenda, §0a items 1–12: R1 merged; CPPM merged; the c ≥ 0.78 route withdrawn; the campaign set as a cyclic fit and PB2; the free surface as the limiter after R1; the low-confinement separation; Kimura et al. (1985) as Gate 1; the Tatsuoka element check; the acoustic census (case C); the initial stiffness as a deliverable; the reviewed separation; K75; the separated zone as a dilatant heave; a third coherence track. |
| 2026-10-03 | Sixth issue, for the owner's review before sending. §0 and §0a rewritten as the state on 3 October, replacing items 1–12: the numerical limit closed and its verification; the Kimura comparison; the grading test; the B/16 result with Perzyna as a prerequisite; the initial stiffness; post-peak dilatancy; the campaign set; one status table of the fork changes. §1 statuses, §3 caution 8, §4 scope note and c = 0.80 results, §4.11 with R1 and the separation, §5 (D1–D11), §5.1 small-strain data, §6 rewritten. |
| 2026-10-03 (later) | Voice pass over §§1–7: impersonal throughout; PR and WP numbers only in §0a.7. The lab-localisation share stated as a bound under an assumed band volume fraction of 0.3. The Kimura B/16 value stated at its reading point (1757 kPa at s/B 0.183, leg running), not as a final −24 %. The separation and the sub-step moduli fix merged. |
