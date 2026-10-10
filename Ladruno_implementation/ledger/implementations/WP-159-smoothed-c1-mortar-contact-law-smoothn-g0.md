---
wp: WP-159
title: "WP-159 — smoothed (C1) mortar contact law: -smoothN <g0> (quadratic normal onset over abs(gbar) < g0 + a C1 friction-onset weight over the same band) and -smoo…"
pr: "#903"
status: "Approved with nits (review 2026-10-02, nits addressed); lands as a fallback law"
section: "table"
legacy_seq: 184
---
| **WP-159 — smoothed (C1) mortar contact law: `-smoothN <g0>` (quadratic normal onset over abs(gbar) < g0 + a C1 friction-onset weight over the same band) and `-smoothT <r>` (the return-map stick/slip corner rounded over abs(rho - cap) < r*cap)** ([[159_mortar_smoothed_contact_law]]; pile-contact R0.8, R3-N1). Opt-in, 3D mortar only, pure penalty (refused unless `-augment never`; refused with `-tie`/`-soft`/`-visc`/NTS; 2D pair = handle-time FATAL). Consistent tangent (the onset coupling eps_N*chi'*tF(x)n only under `-consistanttan`), FD-checked by the oracle (4000 states, 7.6e-9) and on the binary (printA vs FD printB, <= 4.5e-6). Definitions stream v4 -> v5 (+2 mortar slots; v4 still read). Recipe rule: g0 >= c/eps_N for a cohesive interface. R3: bonded lateral 0.83 -> 6.3 mm; alpha lanes still blocked (mortar Tresca slip once engaged from the reference; the shipped law fails the same way) -- R3-N1 not closed. | contact law option | — | `SRC/analysis/handler/LadrunoContactFE.{h,cpp}`, `LadrunoContactHandler.cpp`, `SRC/domain/contact/LadrunoContactDomain.{h,cpp}`, `LadrunoFrictionKernel.h`, `SRC/interpreter/OpenSeesOutputCommands.cpp`, `tests/test_adr159_mortar_smooth_contact.py`, `contact_prototypes/proto_adr159_smooth_normal.py`, `tests/data/adr159_v4_db/` | Approved with nits (review 2026-10-02, nits addressed); lands as a fallback law | [#903](https://github.com/nmorabowen/OpenSees/pull/903) |
