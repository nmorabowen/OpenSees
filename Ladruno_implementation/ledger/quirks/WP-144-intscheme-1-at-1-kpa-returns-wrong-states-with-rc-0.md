---
wp: WP-144
title: "LadrunoSANISAND IntScheme 1 (ModifiedEuler) at p' = 1 kPa with 0.05 % strain steps returns wrong states with rc 0: tension, lost symmetry, sigma1/sigma3 32.7"
date: 2026-11-08
---
### LadrunoSANISAND IntScheme 1 (ModifiedEuler) at p' = 1 kPa with 0.05 % strain steps returns wrong states with rc 0: tension, lost symmetry, sigma1/sigma3 32.7

- **Bites:** single material point, Toyoura (DM04 Table 1, P_atm 101), the campaign options of
  `Ladruno_scripts/sanisand_replay.py` (`CAMPAIGN_OPTS`: IntScheme 1, TolF = TolR = 1e-7, -Pmin 0.0101), driven
  by chained `ladrunoSANISANDReplay` calls from an isotropic 1 kPa state in 0.05 % axial steps (K4,
  `Ladruno_files/testbed/norsand_calib/harness/k4_sanisand.py`). Against WP-134's paper integrator:
  drained plane strain e0 0.833 peaks at sigma1/sigma3 32.7 at the first step (paper 4.805; the out-of-plane
  stress falls to 0.13 kPa); drained TC e0 0.833 loses the symmetry of a symmetric (yy = zz) strain increment
  (sigma_zz 0.025 against sigma_yy 1.0) and reaches lateral TENSION (-3.0 kPa), sigma1/sigma3 13,736. Every step
  returns rc 0. Against DM04 `uw_model` (the UW additions LadrunoSANISAND carries; the right reference for the C++),
  5 of the 6 drained cases at 1 kPa depart and one fails outright; the undrained cases agree. None depart at 10 kPa
  or above. (A first reading against DM04 `paper` also flagged the undrained 1 kPa cases, p_min up to 1.5x: that
  was the UW additions, which every LadrunoSANISAND configuration shows, not this integration error.)
- **Why:** at 1 kPa the elastic range is tiny against a 0.05 % step, and ModifiedEuler's substep error control
  (TolE hardcoded 1e-4 unless honorTolR) does not catch it. WP-129 already measured ModifiedEuler 15-100 % off the
  oracle on hard states; this is the silent form of it at the strip surface's pressure.
- **Fix/status:** use SAS-ME (`IntScheme 129`): exact at 0.05 % steps on the same cases (4.804, 3.249, 3.592
  against the paper's 4.805, 3.250, 3.592). IntScheme 1 converges to the paper only at 0.005 % steps. The TIMs
  strip baseline already runs SAS-ME 129 (workbench `APE/attention.md`, 2026-10-29); runs and scripts still on
  IntScheme 1 (the replay wrapper's `CAMPAIGN_OPTS`, the workbench `labs/ape` scripts, which call it
  "mandatory") are suspect wherever p' approaches 1 kPa. K4 runs both and lists every IntScheme 1 departure.
