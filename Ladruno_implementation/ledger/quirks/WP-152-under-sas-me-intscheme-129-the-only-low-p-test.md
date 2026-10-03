---
wp: WP-152
title: "Under SAS-ME (IntScheme 129) the ONLY low-p test is p + p_r > 0: -Pmin is NOT an admissibility threshold there (WP-152)"
legacy_seq: 534
---
### Under SAS-ME (`IntScheme 129`) the ONLY low-p test is p + p_r > 0: `-Pmin` is NOT an admissibility threshold there (WP-152)
- **Bites:** reasoning about a footing's surface points from `-Pmin`, as for ModifiedEuler.
  - Under ModifiedEuler, `-Pmin` is the threshold of every vanilla clamp and reset: the entry clamp, `Stress_Correction`'s silent deviator wipe to (p_min + p_r)·I, `Elastic2Plastic`, and CPPM's reset.
  - Under SAS-ME it only floors the elastic moduli (`GetElasticModuli`: `sqrt(max(p + p_Re, p_min)/P_atm)`).
  - The refusal tests are `!(p > 0.0)` on p = tr(σ)/3 + p_r: the start (code 3), stage and predictor tension (code 6), and the drift correction (code 7).
  - With the fork default p_r = 0, a SAS-ME point refuses at tr(σ)/3 ≤ 0, whatever `-Pmin` says. And near p = 0, codes 4 and 9 (the α and fabric error, and the substep count) usually fire first.
- **Rule:** Under IntScheme 129 read low-p behaviour from p_r and from the refusal census (`sasStats` `refLowP`, `refDTmin`, `refCap`), not from `-Pmin`. A floor that SAS-ME should honour must be a declared mechanism.
- **Workaround/status:** WP-152's `-sasTensionCutoff p_sep p_contact` is that mechanism, opt-in: it separates only the low-p/tension refusals, and counts them. [[152_sanisand_tension_cutoff]], [[LadrunoSANISAND_implex_guide]] §13.5.
