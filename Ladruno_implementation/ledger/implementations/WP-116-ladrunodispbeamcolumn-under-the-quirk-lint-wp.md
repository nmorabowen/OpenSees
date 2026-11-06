---
wp: WP-116
title: "LadrunoDispBeamColumn under the quirk lint (WP-116)"
pr: "#851"
status: "shipped on the branch; mutation-verified"
section: "table"
legacy_seq: 11
---
| **LadrunoDispBeamColumn under the quirk lint (WP-116)** ([[116_dispbeam_stamp]]) — the 2D/3D element files were unstamped and missing from `stamp_headers.py` GLOBS, so the WP-115 quirk lint never scanned them. Stamped; the five L1 Rayleigh sites (getResistingForceIncInertia ×2 per class, 2D getResponse `dampingForces`) traced safe (nothing on the Rayleigh path writes the static P) and converted to a function-local snapshot, bit-identical. First transient Rayleigh regression for the element: differential vs elasticBeamColumn (lumped/nodal/consistent mass × step/ground × alphaM/betaK/betaK0, 2D+3D) plus a closed-form `dampingForces` = betaK·K_e·v_e check. | hygiene + regression gate | 33013 / 33014 (existing) | `Ladruno_scripts/stamp_headers.py`, `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn{2d,3d}.{cpp,h}`, `tests/test_rayleigh_inertia_dispbeam.py` | **shipped on the branch; mutation-verified** | #851 |
