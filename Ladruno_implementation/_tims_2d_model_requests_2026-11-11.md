# TIMs 2D model — requests to the fork, 2026-11-11 (F25, F26)

Intake file, written by the TIMs Workbench `2d-model` act (branch `work/ape/2d-model`, commit `7ea0d64`, `Tries/2d-model/`) and committed here the way `_tims_2d_model_requests_2026-09-25.md` was. Two requests, both on `LadrunoSANISAND` under SAS-ME (`IntScheme 129`), both opt-in and byte-identical when off.

**How to work it.** Read `CLAUDE.md` and `AGENTS.md` first and hold them. One `wp/<n>-<slug>` branch and a draft PR per item; F26 belongs in WP-156 (#898) and does not need its own. Oracle first (`Ladruno_scripts/sanisand_reference`), then C++. Run the affected pytest suites before and after and paste the output. `// Ladruno` on every touched upstream line. Ledger fragments, guide and banner in the same PR. **Never run the Workbench's models and never edit anything in the Workbench.** Verify every number cited here against this checkout before relying on it; the act read the oracle at `e96f8d77d` and the C++ at `aa6229bf9`. Answer in a `*_tims_report.md` as before.

---

## 0. Where the act stands (context for both asks)

- **The verified stack is unchanged:** SAS-ME 129, R1 (`-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5`, WP-151), low-confinement separation (`-sasTensionCutoff 0.1 1.0`, WP-152), `-Presidual 0`.
- **The sand is now element-calibrated.** A pre-registered fit of the DM04 constants on the WP-134 oracle (`uw_model` options), against Tatsuoka et al. (1986) plane-strain data (22 φ′_peak and 25 ε_peak points, e 0.65–0.86, σc 4.9–392 kPa, plus three curves) and tied to DM04's own 17 simulations of Verdugo & Ishihara (1996). Result, `TOYOURA_CAL1`: `nb 1.1849, h0 6.5557, ch 0.99080, A0 0.66128, nd 2.6471`, everything else DM04 Table 1. It is DM04 within its own uncertainty; three starts reach the same basin.
- **Two gates fail by a pressure trend no bounded constant set removes.** Means of model − lab at σc 4.9 / 49 / 392 kPa: φ′_peak +1.9 / +0.1 / −2.2°; ln(ε_peak ratio) −0.57 / +0.24 / +0.41. DM04 verbatim shows the same trend (+1.2 / −0.5 / −2.5° and −0.76 / +0.07 / +0.26).
- **The owner retired "Track 3"** (`e0 0.85, nb 1.81`): on the same targets it breaks 12 of the 17 triaxial tests.

---

## F25 — a pressure-clamp low-confinement mode in SAS-ME, for an A/B against separation

### What is to be compared

The project's adviser supplied the source of his own bounding-surface sand model, an NTUA-SAND implementation on the UW `ManzariDafalias` skeleton (private communication; the file is not in the Workbench and is not attached here). Read against `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` at `aa6229bf9`, its low-stress treatment differs as follows:

| | stock `ManzariDafalias` | the adviser's model |
|---|---|---|
| `m_Pmin` | 1e-4·P_atm (0.01 kPa) | 1e-3·P_atm (0.1 kPa) |
| `m_Presidual` | 1e-2·P_atm (1 kPa) | 0 |
| stress correction, p below the floor | floor at `m_Pmin + m_Presidual`; the corrective branch sits under `if (false)` | the three normal stresses are set to `m_Pmin`, **the shear components are kept**, and the routine returns |
| integrators, p below the floor at entry | — | `sigma = dev(sigma) + m_Pmin·I` (the deviator is kept) |
| yield drift after a substep | Sloan correction | Sloan correction, then a bisection on the stress increment, then revert to the committed state |
| low-p dilatancy sigmoid (`D_factor`) | present | present, same constants |

The same branch also updates the back-stress ratio by adding the stress change to it. That mixes a stress with a stress ratio; **do not reproduce it**. The option below must state what it does to α, and why.

### Why an A/B

WP-152's separation uses the same 0.1 kPa and the opposite rule: a separated point carries no shear. On the Kimura et al. (1985) footing (`M1a_sym`-type mesh, B/8) 52–58 Gauss points separate at the peak, between 0.03B and 0.6B beyond the footing edge and down to 0.35B, which is the passive wedge. Whether those points carry shear decides whether the slip surface can reach the ground. The model's peak is late (s/B 0.137 against 0.091 measured) and the mechanism dives under the near-surface strip, so the two rules should differ exactly there.

A clamped point holding its deviator at p = 0.1 kPa lies outside the yield cone. That is the thing being tested, not an oversight; the A/B has to measure what it costs.

### Asked

1. An opt-in SAS-ME low-pressure mode, working name `-sasLowP clamp <p_floor>` (default `p_floor` 0.1 kPa), exclusive with `-sasTensionCutoff`: when a point's update would take p below `p_floor`, set p to `p_floor` and keep the deviator. `-Presidual 0`.
2. The same triggers as WP-152 (E1 and E2, with their compression gate), so the A/B changes one thing.
3. A stated, dimensionally consistent treatment of α and α_in for a clamped point, with the reason. Two candidates: leave α as committed, or project it to the cone at the clamped state. Pick one, say why, and report the other's effect on one element path.
4. Per-point counters in the existing `commitStats` / census style: entries, exits, and the largest f/(√(2/3)·m·p) a clamped point holds, so the inadmissibility is measured.
5. Oracle first: `sanisand_reference` gets the same mode; the C++ matches it on the WP-152 element paths.
6. Byte-identical when off, on an existing WP-152 leg's steps.

### The A/B the act will run

Kimura V85.6 and V75.1, unstructured symmetric mesh at B/8, `TOYOURA_CAL1` and DM04, separation against clamp. Read: peak q and its s/B; count, position and depth of affected points; whether the band reaches the surface; the work done by clamped points against the total.

---

## F26 — a separate pressure exponent for the plastic modulus, in WP-156 (#898)

### Evidence

Oracle only (`sanisand_reference` at `e96f8d77d`, `uw_model`), pre-registered, run 2026-11-11. The probe generalised two lines of `model.py`:

- elastic, `G = G0 p_at (2.97−e)²/(1+e) (p/p_at)^ng` (DM04: `ng` = 0.5);
- plastic, `b0 = G0 h0 (1 − ch e) (p/p_at)^−nh` (DM04: `nh` = 0.5).

Both equal DM04 at p = p_at. At 0.5 the probe takes the original `sqrt` branch and reproduces the unpatched cost to every printed digit.

| Cell | Free | nh | ng | Cost (75 residuals) |
|---|---|---|---|---|
| DM04 Table 1 | — | 0.5 | 0.5 | 56.4 |
| `TOYOURA_CAL1` | 5 constants | 0.5 | 0.5 | 45.2 |
| XC | nh, ng; constants at DM04 | 0.206 ± 0.032 | 0.756 ± 0.050 | 30.7 |
| XA | 5 constants + nh | 0.275 ± 0.033 | 0.5 | 27.3 |
| XB | 5 constants + nh, ng | 0.251 ± 0.044 | 0.619 ± 0.073 | 25.7 |

- **Strain at peak.** RMS of ln(model/lab) falls from 0.48 to 0.21–0.25 on the fit set, and from 0.33 to 0.15–0.20 on a hold-out (σc 9.8 and 98 kPa) that no cell was fitted to.
- **The triaxial tie.** With the constants at DM04 (XC), all 17 of DM04's Verdugo–Ishihara simulations stay within ±10 % in q_peak and in q at 25 % strain. With the constants also free (XA), one leaves (`e 0.833, p0 3000`, undrained, q_peak ×1.18).
- **What it does not fix.** At 392 kPa the model stays about 2° low in φ′_peak at every density (44.2 / 42.2 / 40.4 / 37.6° against 46.1 / 44.3 / 42.3 / 39.6°). Neither exponent moves it; the cause is open.
- **Price.** E50 at 4.9 kPa falls to 4.1–5.2 MPa (DM04 8.3–8.7; lab 13.0 and 6.2 on two near-identical specimens).
- **Caveat.** The lab strains are external; bedding error inflates them most at 4.9 kPa.

### Asked

1. In WP-156, give `b0` its own exponent, working name `-sasPlasticExp nh`, default 0.5, apart from the elastic exponent the plan already makes configurable.
2. Keep the elastic exponent a small-strain quantity, set from G_max data. Left free on these curves it goes to 0.62–0.76, which the fit uses to carry plastic strain; that is not evidence for an elastic law.
3. Gate: with the shipped plastic constants and `nh` in 0.2–0.3, DM04's 17 triaxial tests within ±10 % of the `nh` = 0.5 result, on the oracle and in C++.
4. State how `nh` interacts with R1's floor and soft cap, which act on the same modulus.
5. Byte-identical at `nh` = 0.5.
