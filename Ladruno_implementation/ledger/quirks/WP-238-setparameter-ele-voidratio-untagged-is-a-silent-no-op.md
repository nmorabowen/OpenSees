---
wp: WP-238
title: "`setParameter -val e -ele N voidRatio` on a ManzariDafalias-family material is a SILENT NO-OP without the material tag"
pr: "#1120"
date: 2026-11-19
---
### `setParameter -val e -ele N voidRatio` on a ManzariDafalias-family material is a SILENT NO-OP without the material tag (WP-238)
- **Bites:** setting a per-element void ratio (a random field, a layered deposit) on a mesh that uses one material tag.
  - `ManzariDafalias::setParameter` returns −1 unless `argc >= 2` and `argv[1]` is the material's own tag. The element route hands the material only what follows the element tag, so the untagged form reaches no parameter.
  - The command still returns normally. Measured on one stdBrick: the run after `setParameter -val 0.70 -ele 1 voidRatio` was bit-identical to the run without it.
  - The tagged form works but the tag must be a STRING from Python: `ops.setParameter("-val", e, "-ele", 1, "voidRatio", 1)` fails with "Invalid String Input"; `"1"` works. `ops.parameter(tag, "element", 1, "voidRatio", "1")` + `updateParameter` works.
  - This is the second parameter of this class with the defect; `stressCorrection` was the first (ADR-92 P2-4).
  - Vanilla also stores any value (NaN, negative), allows it mid-analysis in the plastic stage, and leaves `mVoidRatio` and the initial elastic operator stale.
- **Rule:** on `LadrunoSANISAND` use the untagged element form; WP-238 claims it, validates it and refuses it for a point with strain history in the plastic stage. On vanilla `ManzariDafalias` pass the material tag as a string, and set it before the first step. Either way read `psi` at one Gauss point before and after to see that the value landed: a parameter that did nothing looks exactly like one that worked.
