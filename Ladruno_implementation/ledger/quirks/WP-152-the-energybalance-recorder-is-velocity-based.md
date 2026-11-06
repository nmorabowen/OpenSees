---
wp: WP-152
title: "The EnergyBalance recorder is velocity-based: under a STATIC integrator every column is zero (WP-152)"
legacy_seq: 537
---
### The `EnergyBalance` recorder is velocity-based: under a STATIC integrator every column is zero (WP-152)
- **Bites:** asking for an energy audit of a quasi-static push (`LoadControl`, `DisplacementControl`) with `recorder EnergyBalance`.
  - IE = ∫ F_resᵀ v dt and ULW = ∫ vᵀ P_ext dt integrate NODAL VELOCITIES (`EnergyBalanceKernel.h`); a static integrator never sets them.
  - Measured: the WP-152 two-brick column, 100 static steps through separation and re-contact: 100 rows, every KE/IE/DW/ULW/RES/ERR = 0.
- **Rule:** Use `EnergyBalance` for transient runs only. For a static push, audit work in the driver (the external work ∫ q·B ds against the material work), or at the material point (net work over closed cycles).
- **Workaround/status:** documented; no change to the recorder.
