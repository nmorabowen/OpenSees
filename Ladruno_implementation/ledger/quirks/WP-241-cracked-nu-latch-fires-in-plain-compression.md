---
wp: WP-241
title: "LadrunoRCConcrete's cracked-nu latch (and TS gate) fires in PLAIN COMPRESSION: eps1 is the total principal strain, so Poisson expansion reaches ft/E with no crack"
pr: "#1132"
date: 2026-11-19
---
### LadrunoRCConcrete's cracked-ν latch (and TS gate) fires in PLAIN COMPRESSION: ε1 is the total principal strain, so Poisson expansion reaches ft/E with no crack (WP-241 adds `-crackGate damage`)
- **Bites:** one PlaneStress quad, `-crackedNu 0 -tensStiff vc`, uniaxial compression, laterally free: the latch fires at |εyy| = ε_cr/ν = 5e-4 (≈ 0.42 f'c for the probe backbone) and the next step diverges (WP-241 probe, reproduced by `tests/test_ladrunoRCConcrete_spatial.py`). The TS floor's gate is the same test, so it can inject a lateral "tension" floor in an uncracked strut. In a beam fibre every compressed fibre is laterally condensed, so there it is the common path, not an edge case.
- **Why:** the latch and the TS gate are `ε1 ≥ ε_cr`, and `ε1` is the principal of the **total** strain (MCFT's composite strain); the spine has ν > 0, so a compressed point expands laterally by ν|ε|. MCFT itself has ν = 0 and never sees this.
- **Fix/status (WP-241):** `-crackGate damage` opens both only once tensile damage has started (the tensile threshold past the tension-backbone peak; the Lubliner tensile measure is zero without a tensile principal stress). It is the default with `-addonPlane spatial`; the membrane default stays `strain` (byte-identical) and a shell or panel deck under compression should pass `-crackGate damage`. Gate: `rc_spatial_gpp.cpp` G6 and `test_ladrunoRCConcrete_spatial.py::test_damage_gate_keeps_plane_stress_compression_uncracked`. See [[LadrunoRCConcrete_guide]] §5.2.
