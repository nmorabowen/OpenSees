---
wp: LEGACY
title: "ASDPlasticMaterial3D test paths must be tuned against PLASTIC response, not elastic estimates"
legacy_seq: 299
---
### ASDPlasticMaterial3D test paths must be tuned against PLASTIC response, not elastic estimates
- **Bites:** a lateral-strain drive sized off the ELASTIC trial stiffness (targeting "+1.2e-3 strain reaches +1.15e3 kPa" against a tension cutoff) never actually reaches the cutoff once plastic relaxation kicks in — the measured final state was `s2=s3=-213.6 kPa`, still sitting on the MC ridge, nowhere near the target the elastic estimate promised. An assertion of non-vacuity (e.g. "the cutoff activates somewhere on this path") built on that estimate silently tests nothing.
- **Why:** ASDPlasticMaterial3D's Backward_Euler return map relaxes the stress well below the elastic-predictor trajectory once yielding starts; extrapolating a target strain/stress from `E`/`nu` alone ignores that relaxation entirely.
- **Workaround/status (2026-08-12, PR #741):** sweep the drive numerically (print the actual committed stress path at a few candidate strain magnitudes) and pin the test's non-vacuity assertions to the MEASURED plastic response, not an elastic back-of-envelope number.
- **Layering note recorded with it:** P1's zero-slave-mass abort catches the FULLY-ghosted NTS `-soft` case by accident of completeness (ghost ⇒ no rank-local element ⇒ zero mass). It can never catch a partition-BOUNDARY node (partial mass, nonzero) and never scans the mortar/edge/plane soft lanes at all — the pre-P2 build ran a 2-rank mortar `-soft` deck to completion silently. The P2 refusal (`ladrunoContactNumRanks() > 1 || hostPartitioned` at the `anySoft` choke point) is the actual guard; when writing a mutation deck for a NEW guard, give the model mass so an OLD guard cannot fire first and let the test pass via the wrong abort.
