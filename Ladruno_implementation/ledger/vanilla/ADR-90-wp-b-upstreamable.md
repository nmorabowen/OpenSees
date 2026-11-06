---
wp: ADR-90
title: "ADR-90 WP-B -- upstreamable-table row(s)"
files: ["`SRC/material/nD/soil/FluidSolidPorousMaterial.cpp`"]
table: "upstreamable"
legacy_seq: [496]
---
| `SRC/material/nD/soil/FluidSolidPorousMaterial.cpp` | `// Ladruno (ADR-90)` **comment-only, no behaviour change.** Marks the copy constructor's `theSoilMaterial = a.theSoilMaterial->getCopy()` (`:156`, the VOID overload) as load-bearing for `tests/test_ladruno_sanisand.py::test_getcopy_void_carries_the_settings_planestrain`: that test's whole route to `LadrunoSANISANDPlaneStrain::getCopy(void)` depends on this call staying dimension-free rather than ever being changed to forward the soil's own type string. No Python-observable exists to PIN this route (both `getCopy(void)` and a hypothetical type-string forward would clone an untouched, freshly-constructed material identically — there is no window to inject state into the wrapper's inner material between its own construction and the `quad` element's construction-time clone), so a code comment plus this row is what stands in for a red-on-regression test. See `LEDGER_quirks.md`'s InitStrain/StagedStrain entry for why no other wrapper reaches a 2D subclass's `getCopy(void)` at all. | ADR-90 WP-B |
