---
wp: LEGACY
title: "OPS_ManzariDafaliasMaterial writes past a 5-element stack array on a deck with >5 trailing optionals"
legacy_seq: 345
---
## `OPS_ManzariDafaliasMaterial` writes past a 5-element stack array on a deck with >5 trailing optionals

**Found 2026-08-27 by adversarial review during ADR-86 PR-3. Vanilla defect, memory-unsafe,
NOT fixed by us — surfaced for a decision per `WORKFLOW_GOTCHAS.md` section 6.**

`SRC/material/nD/UWmaterials/ManzariDafalias.cpp:89-114`:

```cpp
int numArgs = OPS_GetNumRemainingInputArgs();   // :80, counted BEFORE the tag is consumed
...
double oData[5];                                 // :90
...
numData = numArgs - 19;                          // :110
if (numData != 0)
    if (OPS_GetDouble(&numData, oData) != 0) {   // :112  -- writes numData doubles
```

`numArgs` is `19 + k` for `k` trailing optionals, so `numData == k`, **uncapped**. And
`OPS_GetDoubleInput` (`SRC/api/elementAPI_TCL.cpp:316-329`) loops `for (i = 0; i < *numData; i++)
data[i] = ...` with **no bound on the destination array** — it stops only when the input args run
out or a token fails to parse. So a deck supplying more than five trailing numeric arguments
writes past `oData[5]` on the stack:

```tcl
# k = 8 -> writes oData[0..7] into a double[5]: 24 bytes past the array
nDMaterial ManzariDafalias 1 <18 doubles> 1 0 1 1e-7 1e-7 0 0 0
```

This is a THIRD defect in this argument-arithmetic family, and the most serious:

| | defect | consequence |
|---|---|---|
| `SAniSandMS` (ADR 86 sec.7.1) | `numArgs - 19` should be `- 20`; `numData -= 5` should be `-= 3` | `TolF`/`TolR` silently dropped |
| `ManzariDafalias` (this entry) | `numArgs - 19` uncapped into `double[5]` | **stack buffer overflow** |
| `LadrunoSANISAND` | none — see below | hard parse error |

**`LadrunoSANISAND` is immune by construction**, and this is the clearest justification yet for
the design ADR 86 sec.7.1's note argued for. Its parser never computes a count: it consumes every
remaining token one at a time, classifies it, and rejects the sixth positional outright
(`LadrunoSANISAND.cpp`, the `if (nPos >= 5)` guard). `numArgs` is used exactly once, for the
`< 19` minimum check. A count that is never computed cannot be computed wrongly.

- **Do not copy `OPS_ManzariDafaliasMaterial`'s optional-argument block into a new material.**
  It is the pattern to avoid, not the template. Several UW-family parsers share its shape;
  a sweep is owed.
- **FIXED** in the ADR-86 follow-up, with the owner's approval: a `numData > 5` guard that
  **refuses** the deck rather than truncating it (a silent truncation would run the deck with
  arguments the user did not get — the §7.1 defect class again). Gated by
  `tests/test_manzari_safety_pack.py::test_optional_arg_count_is_bounded`, proven by mutation.
- **STILL OWED — the same uncapped pattern is in three more UW parsers**, surveyed in the same
  pass and deliberately left for their own change:
  `ManzariDafaliasRO.cpp:82` (`numArgs - 22` into `double oData[6]`),
  `PM4Sand.cpp:128` (`numArgs - 5` into `oData[24]`),
  `PM4Silt.cpp:129` (`numArgs - 6` into `oData[24]`).
  `SAniSandMS.cpp:134` is **safe from this one** — it bounds its reads with `std::min(numData, 3)`
  and `std::min(numData, 2)` — though it still carries the separate §7.1 silent-drop defect.
