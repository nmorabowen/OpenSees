---
wp: LEGACY
title: "ManzariDafaliasPlaneStrain's null constructor sets the WRONG classTag, and a broker restore never repairs it"
legacy_seq: 344
---
## `ManzariDafaliasPlaneStrain`'s null constructor sets the WRONG classTag, and a broker restore never repairs it

**Found 2026-08-27 by adversarial review during ADR-86 PR-3. Vanilla defect, NOT fixed
by us — surfaced for a decision per `WORKFLOW_GOTCHAS.md` section 6.**

Vanilla is inconsistent between its own two wrappers:

```cpp
// ManzariDafalias3D.cpp:41-42            -- CORRECT
ManzariDafalias3D::ManzariDafalias3D()
  : ManzariDafalias(ND_TAG_ManzariDafalias3D)

// ManzariDafaliasPlaneStrain.cpp:43-44   -- WRONG
ManzariDafaliasPlaneStrain::ManzariDafaliasPlaneStrain()
  : ManzariDafalias()                      // -> NDMaterial(0, ND_TAG_ManzariDafalias)
```

`ManzariDafalias`'s bare null constructor hardcodes `NDMaterial(0, ND_TAG_ManzariDafalias)`
(`ManzariDafalias.cpp:363-364`), so a null-constructed `ManzariDafaliasPlaneStrain` reports
`getClassTag() == ND_TAG_ManzariDafalias` — the DIMENSIONLESS BASE tag.

**And nothing repairs it.** `grep -n setClassTag SRC/material/nD/UWmaterials/ManzariDafalias.cpp`
finds nothing: `recvSelf` restores the 97 data slots and never touches the class tag. So the
object the broker builds at `FEM_ObjectBrokerAllClasses.cpp:2547-2548` carries the wrong tag for
the rest of its life.

Why it is usually invisible: the ordinary `getCopy()` path does `*clone = *this`, and the
compiler-generated member-wise assignment copies `MovableObject::classTag` from a correctly-tagged
source, overwriting the wrong one. The bug only bites on paths that construct the null form and
then READ the tag — i.e. database restore and `OpenSeesMP`. There, if the restored material is
ever re-serialised, `Channel` writes the wrong class tag and the far side brokers a bare
`ManzariDafalias`, whose `getType()` and `getCopy(void)` are "subclass responsibility" + `exit(-1)`.

**We deliberately did NOT replicate it.** `LadrunoSANISANDPlaneStrain`'s null constructor passes
`ND_TAG_LadrunoSANISANDPlaneStrain`, matching what `LadrunoSANISAND3D` and `LadrunoSANISAND`'s own
null constructor already do. For our class the bug would not have been cosmetic: the construction
echo and `Print`'s type guard are BOTH keyed on `getClassTag() == ND_TAG_LadrunoSANISAND`
(`LadrunoSANISAND.cpp`, `echoLadrunoConstants` and `Print`). A restored PlaneStrain carrying the
base tag would therefore have echoed **once per Gauss point** — the exact ~83 MB-of-stderr failure
mode ADR 86 section 4.4's refinement box was written to prevent, arriving through the back door.

- **Do not "align" our wrapper with vanilla's.** The difference is deliberate and load-bearing.
- **If you write a new ND wrapper, pass the explicit class tag to the base constructor.** Copying
  the nearest sibling is how this propagates.
> **CORRECTED 2026-08-28. An earlier version of this entry said the vanilla fix is "one token"
> and "additive". BOTH ARE FALSE, and finding out cost a branch.** The correction is worth more
> than the original entry, so it is kept in full.
>
> **1. The tag cannot be set after construction.** `MovableObject::classTag` is **private**
> (`SRC/actor/actor/MovableObject.h:75`) and there is **no `setClassTag` anywhere in `SRC/`**
> (verified by grep across the tree). So the only route is a base constructor that takes the tag.
>
> **2. The only such base constructor has a SIDE EFFECT.** Diffing
> `ManzariDafalias(int classTag)` (`:302-361`) against `ManzariDafalias()` (`:363-422`) with the
> tag normalised away leaves exactly one substantive difference: the classTag-taking form also
> does **`mElastFlag = 0;`** (`:348`) and the bare form does not. `mElastFlag` is a
> **`static char unsigned`** (`ManzariDafalias.h:202`, defined `= 1` at `ManzariDafalias.cpp:58`)
> — the process-wide stage flag of ADR-86 risk 3, whose whole hazard is that constructing ANY
> Manzari-family material resets the stage for EVERY instance in the process. So repairing the
> tag would newly reset the elastic stage on every broker / database-restore construction of a
> `ManzariDafaliasPlaneStrain`. That is a live behaviour change on the MP path, not a no-op.
>
> **3. It is three wrappers, not one.** `ManzariDafalias3DRO` and `ManzariDafaliasPlaneStrainRO`
> both delegate to `ManzariDafaliasRO()`, which itself delegates to the bare `ManzariDafalias()` —
> and `ManzariDafaliasRO` has **no tag-taking constructor at all** (`ManzariDafaliasRO.h:54-59`).
> Fixing those two additionally requires ADDING a constructor overload to a vanilla class.
> Only `ManzariDafalias3D` is correct today, and it is correct *because* it routes through the
> tag-taking form — which means 3D restores already reset `mElastFlag` and PlaneStrain restores do
> not. That asymmetry exists in vanilla right now.
>
> **DECISION (owner, 2026-08-28): NOT FIXED.** Every available route either carries the
> `mElastFlag` side effect or grows the vanilla API, and the defect is narrow — it bites only when
> an already-broker-restored object is RE-serialised. Recorded here instead. Per D8 the upstream
> call is not ours to force.
>
> **If anyone revisits this:** the thing to look for is a route that sets the tag and nothing else.
> A new base constructor taking only the tag would do it, at the cost of a near-duplicate of a
> 55-line constructor body. Do not simply add `mElastFlag = 0` to the bare form to "make them
> consistent" — that changes every null-constructed Manzari material in the process.
>
> **RECONFIRMED 2026-09-04, ADR-90 WP-B.** The ADR-90 planning brief (F5) listed this as a
> WP-B prerequisite fix, assuming a fix confined to `ManzariDafaliasPlaneStrain.cpp` alone was
> available. It is not: `MovableObject::classTag` is private with no setter anywhere in `SRC/`
> (re-verified), so the only way to change it is a base-constructor call, and every such call
> either carries the `mElastFlag` side effect or requires adding a new constructor overload to
> vanilla `ManzariDafalias` — exactly the two options the 2026-08-28 decision already weighed and
> declined. **Status stays NOT FIXED.** What WP-B added instead is a broker/database round-trip
> test that reproduces the defect end to end —
> `tests/test_manzari_planestrain_classtag_quirk.py::test_manzari_planestrain_classtag_survives_one_roundtrip_but_not_two`
> — so a future change to the null constructor, the broker dispatch table, or the wire format has
> something concrete to check against instead of this prose alone.
