---
wp: PR-444
title: "444 -- 2 vanilla row(s)"
pr: "#444"
files: ["`SRC/domain/domain/Domain.cpp`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [35, 39]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno` ADR-60: `Domain::commit()` — after the ADR-39 `theContactDomain->commit()`, raise `domainChange()` when `needsResort(this)` (a slave migrated past the broad-phase band) so the next step re-handles + re-emits the NTS candidate set. Gated `!contactAugmenting` (GA-1) + `theContactDomain!=0`; OFF ⇒ `needsResort`==false ⇒ byte-identical. | [#444](https://github.com/nmorabowen/OpenSees/pull/444) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-60: `OPS_LadrunoContact` gains `-reemit` / `-resortFrac` / `-resortEvery` (finite-sliding NTS re-emit opt-in; passed to `addContact`). Refused with `-mortar`. OFF default ⇒ byte-identical. | [#444](https://github.com/nmorabowen/OpenSees/pull/444) |
