# WP-145 — LadrunoHySAND (near-future plan, not started)

**Status: PLANNED.** Recorded 2026-09-27 on the owner's decision after the TIMs F18–F23 program. No code yet. The ND tags are reserved in `LEDGER_implementations.md` (**33026** base, **33027** 3D, **33028** PlaneStrain); they are not in `classTags.h` until the code lands, following the ADR-90 precedent.

## Why

**HySand** (Simonin, Houlsby & Byrne 2026, *Géotechnique* 76(13), doi 10.1680/jgeot.25.00091) is a **hyperplastic multisurface sand model**. It is derived from a Gibbs energy plus convex yield surfaces, so non-negative dissipation holds by construction. In one parameter set it has density lines (B, Γ, Δ), density- and anisotropy-dependent dilation, a consolidation mechanism, and peak and softening. The model survey ([[_sand_model_survey_2026-09-27]]) ranks it **#3** for the TIMs footing, and first for the **long-term, thermodynamically consistent cyclic** role. That is the role DM04 SANISAND fills today without an energy function for α.

It complements WP-144 (LadrunoNORSAND):
- **NorSand-BA** is the consistent **monotonic/footing** model: small state, closed-form tangent, weeks.
- **HySand** is the consistent **cyclic/SSI** candidate: multisurface state, research-grade.

## Known risks (from the survey; tags there)

- **Very new** (2026): no public code found [E, none found]; ABAQUS and PLAXIS UMATs exist but are not public [E]; the integration scheme is in Simonin (2023, DPhil) and Houlsby (2025a), not read [E-sec].
- **No footing evidence**: validation was on monopiles, partly confidential [E].
- **Low confinement:** its own 3D FE needed a **10 kPa surface surcharge "to ensure convergence at low stress"** [E]. The TIMs ring (p' ≈ 0.3–10 kPa) is exactly that regime, so the p'-floor question does not go away.
- **Cost:** a multisurface state of N × (6+) internal variables per point [I].
- **Effort:** the survey estimates about **10–14 weeks** [I]. It is a research build, not a WP-sized port.

## Plan (sketch; to be refined once the integration papers are read)

| Phase | Content |
|---|---|
| P-1 (gate) | Obtain and read Simonin 2023 (DPhil) and Houlsby 2025a (integration); contact the authors about code availability/licensing; decide build vs collaborate. **No code before this.** |
| P0 | Python kernel oracle (multisurface hyperplastic update from the energy + dissipation functions), with D ≥ 0 and energy-balance checks at every step |
| P1 | C++ class + 3D/PlaneStrain wrappers; implicit update; refusal via the WP-99 latch; bounded work |
| P2 | Tests: oracle parity, FD tangent, dissipation census, cyclic drained/undrained elements |
| P3 | Calibration (14 parameters; the CSL maps from DM04, the rest refit from monotonic + cyclic triaxials); the strip deck; then a cyclic SSI benchmark |
| Gate | Adversarial review; ledgers, guide, banner |

**Start condition:** after WP-144 (NorSand) is done, or earlier if TIMs' SSI work needs a consistent cyclic model before then. P-1 (reading and author contact) can start any time and costs nothing.

## Sources

The survey ([[_sand_model_survey_2026-09-27]]); Houlsby & Puzrin (hyperplasticity); Collins & Houlsby 1997; Simonin, Houlsby & Byrne 2026.
