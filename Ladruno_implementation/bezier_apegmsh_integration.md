---
title: BezierTri6 — apeGmsh integration
project: Ladruno
status: working (direct-drive); typed-primitive deferred
owner: nmora
tags:
  - implementation
  - element
  - apegmsh
  - integration
---

# BezierTri6 — apeGmsh integration

How [apeGmsh](https://github.com/nmorabowen/apeGmsh) consumes the
[[04_bezier_elements|BezierTri6]] element (merged to `ladruno`, `ELE_TAG 33000`).
Companion: the **canonical apeGmsh-side implementation spec** lives in the apeGmsh
repo at `docs/plans/bezier-tri6-element.md` (the 6 touch points to make
`ops.element.BezierTri6()` a first-class typed primitive). This file is the
**OpenSees/ladruno-side** view: what works today, the environment split, and the
validated round-trip.

## Two environments (the key constraint)

| | Python | has |
|---|---|---|
| **apeGmsh** | 3.11 | `gmsh` + `apeGmsh` (from `Github/apeGmsh/src`); `import openseespy` |
| **Ladruno fork build** | 3.12 | our `opensees.pyd` (`dist/bin`), where `BezierTri6` lives |

`BezierTri6` is **only** in the fork build — not in PyPI `openseespy`, not upstream.
So apeGmsh's in-process bridge (`apeSees(...).run()`, which imports `openseespy`)
cannot instantiate it, and `apeSees`'s typed-element registry doesn't list it.

## What works today — direct-drive (validated)

Use apeGmsh purely as the **mesh generator**; drive the fork build directly:

1. **py3.11 + apeGmsh** — mesh a straight-sided domain to quadratic (T6) triangles,
   extract `fem`, dump `nodes` + `fem.elements.connectivity` (the T6 rows) to JSON.
2. **py3.12 + `import opensees`** — read the JSON; for each element call
   `ops.element('BezierTri6', eid, *conn6, thick, type, matTag)`.

**Why it just works:** Gmsh's native `tri6` (etype 9) node order — 3 corners then
mid-edges (1-2, 2-3, 3-1) — is **identical** to BezierTri6's. Verified to machine
precision (max mid-node deviation `1.78e-15`) on a real mesh, which simultaneously
confirms ordering **and** straight-sidedness (mid-nodes at edge midpoints ⇒ `P = X`,
the v1 requirement). `fem.elements.connectivity` rows are usable **verbatim** — no
remap.

**Validated round-trip:** cantilever plate (10×2), 465 nodes / 208 straight-sided T6
elements → solve converged, **tip deflection within 0.3 % of the Timoshenko estimate**,
peak σ_xx = Mc/I. Worked two-step scripts at the OpenSees repo root (local helpers,
git-ignored): `bezier_apegmsh_mesh.py` (mesh→JSON) + `bezier_apegmsh_run.py`
(JSON→BezierTri6). The tracked, build-only regression is
`Ladruno_scripts/bezier_tests/test_bezier_tri6.py`.

## The gap — first-class apeGmsh typed primitive (deferred)

To get `ops.element.BezierTri6(pg=...)` (bridge emission + run) instead of
direct-drive, apeGmsh needs a typed-primitive wrapper. Spec written:
`apeGmsh/docs/plans/bezier-tri6-element.md` — 6 touch points modeled on the existing
`SixNodeTri` (element dataclass, `__init__` export, namespace factory, capabilities,
response catalog with `ELE_TAG 33000`), plus the `-bbar`/`-cMass` flags. Caveats:

- **Tcl/Py deck emission** (`ops.tcl/py`) works once the wrapper exists — the deck
  is just `element BezierTri6 …`, run with the fork's `OpenSees.exe` / `import opensees`.
- **Live `run()`** needs the fork build to *be* the `openseespy` apeGmsh imports.
- **Results parsing** keys on `cpp_class_name="BezierTri6"` + `ELE_TAG 33000`, so the
  result round-trip works once the capability/response entries land.

Until then, **direct-drive is the supported path** and needs no apeGmsh change.

## Scope reminder (v1)

Straight-sided meshes only — the element warns if a mid-edge node is off the edge
midpoint. apeGmsh meshes of straight-edged domains satisfy this automatically.
Curved geometry + the Eq.14 Dirichlet mapping is a deferred fork follow-up
([[04_bezier_elements]] *Deferred follow-ups*).
