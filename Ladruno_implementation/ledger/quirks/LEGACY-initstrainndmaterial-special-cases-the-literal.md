---
wp: LEGACY
title: "InitStrainNDMaterial special-cases the literal string \"ThreeDimensional\" when forwarding to getCopy(void); StagedStrainNDMaterial has NO such special case at a…"
legacy_seq: 369
---
## `InitStrainNDMaterial` special-cases the literal string "ThreeDimensional" when forwarding to `getCopy(void)`; `StagedStrainNDMaterial` has NO such special case at all — a `getCopy(void)` override for a 2D subclass has NO reachable caller through either wrapper

**Found 2026-09-04, ADR-90 WP-B, while covering `LadrunoSANISANDPlaneStrain::getCopy(void)`.
CORRECTED 2026-09-04 after code review: the original version of this entry claimed both wrappers
special-case "ThreeDimensional" the same way. Only `InitStrainNDMaterial` does; `StagedStrainNDMaterial`
does not special-case anything, at all — re-read the two side by side below.**

`InitStrainNDMaterial::getCopy(const char *type)` (`InitStrainNDMaterial.cpp:289-327`) special-cases
exactly the literal string `"ThreeDimensional"` (`:292-293`, `if (strcmp(type, "ThreeDimensional")
== 0) return getCopy();`) to forward to its OWN void `getCopy()` — which is what finally calls
`theMaterial->getCopy()` (void) on the inner. That coincides with the "ThreeDimensional" template
its constructor always builds (`material.getCopy("ThreeDimensional")`, `:135,153`), which is why the
3D coverage trick below works. Every OTHER type string (`"PlaneStrain"`, `"AxiSymmetric"`, …)
instead calls `theMaterial->getCopy(type)` (`:310`) — the inner's TYPE-STRING constructor, never its
void clone.

`StagedStrainNDMaterial::getCopy(const char *type)` (`StagedStrainNDMaterial.cpp:293-309`) has NO
such branch: it UNCONDITIONALLY calls `theMaterial->getCopy(type)` (`:303`) for every string,
including `"ThreeDimensional"` itself — there is no `if` on `type` at all. So there is no type
string, not even `"ThreeDimensional"`, that routes a `StagedStrainNDMaterial::getCopy(const char*)`
call to the inner's void `getCopy()`. The ONLY way to reach the inner's void `getCopy()` through
`StagedStrainNDMaterial` is `StagedStrainNDMaterial::getCopy(void)` itself (`:286-291`,
`adoptCopy(theMaterial->getCopy())`) — a bare, zero-argument call, which nothing in this fork's
element/wrapper ecosystem ever issues against an already-live wrapper instance.

**Consequence:** the 3D getCopy(void)-coverage trick (`test_ladruno_sanisand.py`'s
`test_getcopy_void_carries_the_settings` — wrap the prototype in `InitStrain`, drive it with a
3D-only element like `stdBrick` so the element's own type request happens to be
`"ThreeDimensional"`, matching `InitStrain`'s hardcoded template dimension AND its one special
case) **does not port to any 2D subclass, and does not port to `StagedStrainNDMaterial` at ANY
dimension** (it has no special case to exploit in the first place). Wrapping a SANISAND/
ManzariDafalias-family prototype in `InitStrain` or `StagedStrain` and driving it with
`quad ... PlaneStrain` never reaches the PlaneStrain subclass's `getCopy(void)` override at all —
it only re-exercises the base class's `getCopy(const char*)` PlaneStrain branch (typically already
covered directly, unwrapped).

**The route that DOES reach it:** `FluidSolidPorousMaterial` (`SRC/material/nD/soil/
FluidSolidPorousMaterial.cpp`) asks the soil for a DIMENSION-MATCHED copy at its own construction
(`soilMat.getCopy(nd == 3 ? "ThreeDimensional" : "PlaneStrain")`, `:125`), so its `theSoilMaterial`
is already the correctly-dimensioned subclass. Its own `getCopy(const char*)` (`:415-425`) does not
forward the type string to the soil at all — for EVERY type it accepts, it goes through its own
COPY CONSTRUCTOR (`:152-162`), which clones the soil with the dimension-free, no-argument
`theSoilMaterial->getCopy()`. So `nDMaterial FluidSolidPorous $wrap 2 $soilTag $K $atm` +
`element quad ... PlaneStrain $wrap` calls the 2D subclass's `getCopy(void)` once per Gauss point,
regardless of what type string the element passes. Leave the wrapper's own `loadStage` at its
construction default of 0 (never call `updateMaterialStage` on the WRAPPER's tag) and
`FluidSolidPorousMaterial::getStress()` is then an exact pass-through of the soil's own stress, so a
wrapped-vs-unwrapped bit-identical comparison isolates `getCopy(void)` alone. Pattern:
`tests/test_ladruno_sanisand.py::test_getcopy_void_carries_the_settings_planestrain`.

- **If you need to cover a 2D (or AxiSymmetric) subclass's `getCopy(void)`, reach for
  `FluidSolidPorousMaterial`, not `InitStrain`/`StagedStrain`.** Only the 3D case gets the
  convenient coincidence.
- **If you write a new generic ND wrapper and want its `getCopy(const char*)` to reach the inner's
  `getCopy(void)` for a NON-3D view too**, you would need to special-case that view's OWN string the
  same way `"ThreeDimensional"` is special-cased today — nothing currently shipped does this.
