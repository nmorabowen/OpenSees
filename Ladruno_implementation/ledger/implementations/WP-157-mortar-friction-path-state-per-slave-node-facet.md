---
wp: WP-157
title: "WP-157 — mortar friction path state per (slave node, facet pair)"
pr: "#900"
status: "Merged 2026-10-01 (e6275bfa9)"
section: "table"
legacy_seq: 182
---
| **WP-157 — mortar friction path state per (slave node, facet pair)** ([[157_mortar_friction_pair_state]]; resolves LEDGER_quirks MAJOR-1 of the ADR-41 C3.1 gate). `LadrunoContactDomain::MortarFrictionState` keyed (contactTag, slave node, slave-facet ordinal, master-facet ordinal) replaces the friction fields of the per-node `MortarNormalState` (`λ_N` stays per node); `LadrunoContactFE::setMortarMasterFacet` + handler GC marks (3D + 2D); 3D/2D/SOFT=2 mortar friction and both friction tangents read the pair slot; commit/revert/revertToStart iterate it (capstone contract #1; `λ_T` promotion gated by WP-155 `-augment`). No new command, no classTag. One pair per node ⇒ bit-identical. Oracle `proto_adr157_mortar_pair_friction.py` (T1–T4); tests `test_adr157_mortar_pair_friction.py` 12/12 (10/12 on d63f49750). | contact fix | — | `SRC/domain/contact/LadrunoContactDomain.{h,cpp}`, `SRC/analysis/handler/LadrunoContactFE.{h,cpp}`, `SRC/analysis/handler/LadrunoContactHandler.cpp` | Merged 2026-10-01 (`e6275bfa9`) | [#900](https://github.com/nmorabowen/OpenSees/pull/900) |
