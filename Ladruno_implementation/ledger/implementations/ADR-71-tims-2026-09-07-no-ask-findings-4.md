---
wp: ADR-71
title: "TIMs 2026-09-07 no-ask findings (4)"
status: "shipped — wp/29c"
section: "table"
legacy_seq: 167
---
| **TIMs 2026-09-07 no-ask findings (4)** — (1) `LadrunoKinematicCoupling` parser REFUSES a slave whose ndf is neither `ndm` nor `ndm + nrot` when `-dof` is omitted (the default component list tied an ndf-4 u-p slave's pressure DOF to θx of the master, silently; `resolveGeometry`'s "lacks DOF" message is behind `!useDefault`); explicit `-dof` unchanged. (2) ADR-71 §3.2 datum rider reworded to what the xfail measures (sealed static is SILENT, solvers return rc = 0). (3) `ladruno_apegmsh_contract.md`: `Results.from_ladruno` marked SHIPPED (`apeGmsh/results/Results.py:491`). (4) static `mElastFlag` sentence in the SANISAND apeGmsh emitter guide. | parser refusal + docs | none | SRC: `SRC/element/ladrunoKinematicCoupling/OPS_LadrunoKinematicCoupling.cpp`. tests: `tests/test_ladrunoKinematicCoupling_element.py` (3 new). docs: `71_ladruno_up_family_adr.md`, `ladruno_apegmsh_contract.md`, `86_ladruno_sanisand_apegmsh_emitter_guide.md`, `LadrunoKinematicCoupling_guide.md`. | **shipped — wp/29c** | PR pending |
