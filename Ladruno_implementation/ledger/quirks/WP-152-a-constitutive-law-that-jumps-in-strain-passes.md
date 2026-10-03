---
wp: WP-152
title: "A constitutive law that JUMPS in strain passes every material-point test and has NO equilibrium in a structure (WP-152)"
legacy_seq: 536
---
### A constitutive law that JUMPS in strain passes every material-point test and has NO equilibrium in a structure (WP-152)
- **Bites:** a state switch that sets the stress discontinuously at a strain threshold (WP-152's first re-contact: p_min → p_contact at g = g_c).
  - A strain-driven material-point test never sees it: it prescribes the strain, so the stress just jumps.
  - With compliance around the point (a free node, a neighbouring element) there is a band of load with no equilibrium: for the post-jump stress the neighbours must yield, which moves the strain back below the threshold. Newton oscillates between the branches; a step cut cannot help, because the band has a finite width (≈ Δσ / neighbour modulus).
  - Measured: a two-brick oedometric column under gravity, top displacement-controlled, Newton + `NormUnbalance`: step 68 (the re-contact) failed at 30 iterations; band ≈ 1 kPa / 8e4 kPa ≈ 1.25e-5 > the 1e-5 step.
- **Rule:** Make every branch of a state machine continuous in the strain (a stress reached at the switch, not set there), and test a new material state machine under NEWTON with at least one free node, not only strain-driven.
- **Workaround/status:** ✅ WP-152 (continuous closing branch p_min + K(p_contact)·max(g, 0)); `tests/test_ladruno_sanisand_tension_cutoff.py::test_newton_column_under_gravity_separates_and_recontacts`. [[152_sanisand_tension_cutoff]].
