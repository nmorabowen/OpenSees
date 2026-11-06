---
wp: PR-411
title: "411 -- 2 vanilla row(s)"
pr: "#411"
files: ["`SRC/domain/domain/Domain.{h,cpp}`", "`SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`"]
table: "main"
legacy_seq: [64, 65]
---
| `SRC/domain/domain/Domain.{h,cpp}` | `// Ladruno` ADR-41 D1: add a `bool contactAugmenting` flag (in-class `=false` + init `false` in all 4 ctors; declared last → no `-Wreorder`) + inline `set/isContactAugmenting`. In `Domain::commit()` the recorder loop + `commitTag++` are wrapped in `if (!contactAugmenting) { … }` so a held-load within-step augmentation sweep (`analyze_augmented` at a zero load increment) fires NO recorders and bumps NO commitTag — the contact Uzawa λ update (`theContactDomain->commit()`) and `committedTime=currentTime` still run (committedTime a no-op under the zero-increment LoadControl). Flag OFF (default) ⇒ `commit()` byte-identical to stock. Driven by `ladrunoBeginAugment/EndAugment`. | [#411](https://github.com/nmorabowen/OpenSees/pull/411) |
| `SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-41 D1: add the `ladrunoBeginAugment` / `ladrunoEndAugment` commands → `OPS_LadrunoBeginAugment`/`EndAugment` set `Domain::setContactAugmenting(true/false)` (open/close a held-load within-step augmentation sweep; no args, null-domain guarded, idempotent). OPS_ bodies + decls + Py (`Py_ops_LadrunoBeginAugment`/`EndAugment`) + Tcl (`Tcl_ops_…`) wrappers + dual `addCommand`. Additive; no existing command touched. | [#411](https://github.com/nmorabowen/OpenSees/pull/411) |
