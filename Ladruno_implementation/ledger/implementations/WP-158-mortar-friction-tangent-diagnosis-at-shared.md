---
wp: WP-158
title: "WP-158 — mortar friction tangent diagnosis at shared multi-pair nodes + the -consistanttan recipe"
pr: "#902"
status: "Merged 2026-10-01 (3144e19ba)"
section: "table"
legacy_seq: 183
---
| **WP-158 — mortar friction tangent diagnosis at shared multi-pair nodes + the `-consistanttan` recipe** ([[158_mortar_tangent_diagnosis_consistanttan]]; answers ADR-157 §5 / review #900 findings 1-2). Diagnosis (FD of `printB` vs `printA`): at shared crease nodes the analytic mortar tangent misses the Coulomb pressure coupling Csl (dropped by the symmetric default; `-consistanttan` supplies it) and the cross-crease geometric dD/du, dM/du (no shipped tangent). Ships NO new code path (review #902 rec. A option 2): the opt-in FD pair tangent is parked as `contact_prototypes/adr158_fd_pair_tangent_oracle.patch` (diagnostic oracle; no `-fdTangent`, no wire bump, DB format stays v4). Recipe (rec. B): `-consistanttan` + a non-symmetric solver for mortar Coulomb on curved / faceted / non-matching interfaces, NOT the default (symmetric-SOE contract, ADR-39 P3.5 Q2). Tests: `test_adr158_mortar_consistanttan_multipair.py` (mu>0 multi-pair solid roof: analytic force + <= 6 its with `-consistanttan`, fails on d63f49750; shipped-default linear twin `(not ok) or max(its) > 12`; same converged state). Probes + byte-dump plugin in `contact_prototypes/`. | recipe + tests (no code) | — | `tests/test_adr158_mortar_consistanttan_multipair.py`, `Ladruno_implementation/contact_prototypes/{adr158_fd_pair_tangent_oracle.patch, probe_adr158_*.py, bytedump_plugin.py}` | Merged 2026-10-01 (`3144e19ba`) | [#902](https://github.com/nmorabowen/OpenSees/pull/902) |
