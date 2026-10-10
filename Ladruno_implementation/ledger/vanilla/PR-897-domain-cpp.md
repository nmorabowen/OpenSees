---
wp: PR-897
title: "897 -- 2 vanilla row(s)"
pr: "#897"
files: ["`SRC/domain/domain/Domain.cpp`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [32, 33]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno (ADR-155)`: `Domain::commit()` passes its held-load bracket flag to the contact engine — `theContactDomain->commit(contactAugmenting)` (was `commit()`), so a `-augment request` mortar contact augments only inside `ladrunoBeginAugment`/`ladrunoEndAugment`. Default contacts ignore the argument ⇒ byte-identical. | [#897](https://github.com/nmorabowen/OpenSees/pull/897) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno ADR-155`: `OPS_LadrunoContact` gains the `-mortar`-only options `-augment commit\|request\|never`, `-maxGap d`, `-gapOffset g0`, `-adjust [tol]` (validated after the loop: refused without `-mortar`, gap shift refused with `-tie`), applied through `LadrunoContactDomain::setMortarContactOptions` only when one was given ⇒ byte-identical otherwise. | [#897](https://github.com/nmorabowen/OpenSees/pull/897) |
