---
wp: LEGACY
title: "SRC/system_of_eqn/linearSOE/pardiso/ holds TWO symmetric PARDISO implementations — one of them is dead and unbuilt"
legacy_seq: 211
---
### `SRC/system_of_eqn/linearSOE/pardiso/` holds TWO symmetric PARDISO implementations — one of them is dead and unbuilt
- **Bites:** an agent asked to "add symmetric PARDISO" finds `PARDISOSymLinSOE.{h,cpp}` + `PARDISOSymLinSolver.{h,cpp}` already sitting in the directory and wires *those* up. They are a 2019 contributed prototype (M. Salehi, same author as the `Gen` pair) and are **not listed in `pardiso/CMakeLists.txt`** — only the `Gen` pair is compiled. They carry every defect ADR-75 P1a had to fix in the `Gen` pair and one more of their own: their `setSize` appends adjacency entries **in ID order without sorting**, relying on the adjacency already being ascending, so they would hand PARDISO a CSR whose columns are not guaranteed ascending.
- **Why:** the live symmetric path is `PARDISOGenLinSOE` with `matType != 0` (ADR-75 P1d, [#630](https://github.com/nmorabowen/OpenSees/pull/630)) — one class covering unsym + SPD + symmetric-indefinite, so the hardened factorization-reuse/`mtype`-derivation logic is shared rather than duplicated. The `Sym` pair was left untouched rather than deleted, since it is upstream-contributed and touching it widens the vanilla footprint for no gain.
- **Workaround/status:** use `system Pardiso -matrixType 1|2`. Do **not** wire up `PARDISOSymLin*`; if it ever gets built, its unsorted-column fill is the first thing to fix. Note the ascending-column requirement is now *checked* at the end of `PARDISOGenLinSOE::setSize`, so a future mistake here fails loudly instead of returning a plausible wrong answer. *2026-07-25 (ADR-75 P1d adversarial review).*
