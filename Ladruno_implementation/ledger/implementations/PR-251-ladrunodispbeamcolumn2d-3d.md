---
wp: PR-251
title: "LadrunoDispBeamColumn2d / 3d"
pr: "#251, #254, #255, #258, #260, #267, #269, #271, #274"
status: "shipped"
class_tags: ["ELE_TAG_LadrunoDispBeamColumn2d=33013", "ELE_TAG_LadrunoDispBeamColumn3d=33014"]
banner: "LadrunoDispBeamColumn — regularized disp beam + cohesive hinges"
section: "table"
legacy_seq: 93
---
| **LadrunoDispBeamColumn2d / 3d** ([[32_ladruno_dispbeamcolumn_regularization_adr]], [[33_ladruno_dispbeamcolumn3d_hinge_adr]]) — regularized displacement-based frame element: **Tier-1** per-IP crack-band `lch` channel (mirrors `ForceBeamColumn`, fixes [[LEDGER_quirks]] §59) + Corotational large-disp + `-nl` ½θ² bowing; **Tier-2** embedded strong-discontinuity cohesive rotation-jump hinge (`-hinge`/`-hingeY`/`-hingeBiaxial`, Armero–Ehrlich, guarded static condensation to the basic system before `crdTransf`). No-hinge path bit-identical. **WP-200 parity:** unified `ndm`-dispatch `element LadrunoDispBeamColumn`, fail-closed `LadrunoOptSpec.h` option parser, `.ladruno` HDF5 `endForces` & `hingeState` extraction, OpenMP `ladrunoThreadSafeUpdate/Form/Commit` opt-in (`thread_local` statics). Full build history in the section below (Stage 1 + Tier-2 PR-2a…PR-4b). | Element | **ELE_TAG 33013 (2d) / 33014 (3d)** (Element band; numerically equal to ND_TAG 33013/33014 — per-registry, not a collision) | `SRC/element/ladrunoDispBeamColumn/`, `tests/test_ladrunoDispBeamColumn{2,3}d*.py`, `tests/test_ladruno_beam_consistency.py` | shipped | [#251](https://github.com/nmorabowen/OpenSees/pull/251), [#254](https://github.com/nmorabowen/OpenSees/pull/254), [#255](https://github.com/nmorabowen/OpenSees/pull/255), [#258](https://github.com/nmorabowen/OpenSees/pull/258), [#260](https://github.com/nmorabowen/OpenSees/pull/260), [#267](https://github.com/nmorabowen/OpenSees/pull/267), [#269](https://github.com/nmorabowen/OpenSees/pull/269), [#271](https://github.com/nmorabowen/OpenSees/pull/271), [#274](https://github.com/nmorabowen/OpenSees/pull/274) |
