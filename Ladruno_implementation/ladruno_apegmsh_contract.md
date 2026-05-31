---
title: Ladruno → apeGmsh feature reference
project: Ladruno
status: reference
owner: nmora
tags:
  - apegmsh
  - integration
  - reference
---

# Ladruno → apeGmsh feature reference

apeGmsh ([github.com/nmorabowen/apeGmsh](https://github.com/nmorabowen/apeGmsh))
is the driver/helper for this OpenSees fork: it meshes geometry, emits
Tcl/openseespy, drives the run, and reads results back. This file is the
**apeGmsh-facing record** of the fork-only features — *what exists, the command to
emit, and what to read back* — so apeGmsh work (or an agent on the apeGmsh side)
has one place to look.

It's a plain reference doc, kept current by hand. The companion is
[[LEDGER_implementations]] (which is **authoritative** for each feature's class tag
and shipping PRs); this file adds the apeGmsh-facing *how*. When a fork feature
ships or its command grammar changes, update the row here in the same PR — same
habit as the ledgers, no tooling.

> **Vanilla openseespy is never required.** The fork is opt-in *enhancement*.
> apeGmsh must keep working on stock `openseespy` for everything vanilla OpenSees
> provides; the fork-only features below are simply unavailable there. apeGmsh
> should gate a fork-only feature **at the point of use** (a clear "requires the
> Ladruno fork build" error when you reach for it), never force the fork or fail
> at import.

> **Class tags moved.** Ladruno now allocates fork tags in a **private ≥33000
> band** (BezierTri6, ExplicitBathe, EnergyBalance = 33000 in their respective
> `ELE_/INTEGRATOR_/RECORDER_` namespaces; LNVD 33002; MPCO 33001; CDL 33003). The
> old sub-300 values (272/61/63/26/27/64) are **dead** — don't hardcode them in
> `_response_catalog`. Read the live values from `SRC/classTags.h` /
> [[LEDGER_implementations]], which stay authoritative.

## Feature → apeGmsh reference table

| Feature | OpenSees command (canonical emission) | apeGmsh touch-points | Status |
|---|---|---|---|
| **BezierTri6** | `element BezierTri6 tag n1 n2 n3 n4 n5 n6 thick type matTag [-bbar] [-cMass] [-pressure p] [-rho r] [-bodyForce b1 b2]` | `_response_catalog` `ELE_TAG_BezierTri6` (read from classTags.h, **33000** band); `_ELEM_REGISTRY` `_ElemSpec(gmsh_etype=9, node_reorder=identity` (verified `1.78e-15`)`, slots=(nodes,thick,type,matTag), no transf, has_gauss)`; typed element class `_emit`. `-bbar` valid PlaneStrain/3D only (warns + ignored for PlaneStress). | shipped (direct-drive today; typed primitive deferred — [[bezier_apegmsh_integration]]) |
| **ExplicitBathe** | `integrator ExplicitBathe p [-cfl] [-cflAbort] [-tangent] [-recompute N] [-lump diagonal\|rowsum] [-verbose] [-divergence f]` | integrator emit; `p` **first positional**, `∈(0,1)`, default `0.54`. Exposes `criticalTimeStep()` query. | shipped |
| **ExplicitBatheLNVD** | `integrator ExplicitBatheLNVD p [alpha] [...same opts as ExplicitBathe...]` | as ExplicitBathe; `alpha ∈(0,1)` default `0.8` (FLAC damping). | shipped |
| **CentralDifferenceLadruno** | `integrator CentralDifferenceLadruno [...]` (options mirror the explicit Bathe integrators; correct first-step starter, built-in `dt_cr`, βK guard) | integrator emit. **Coupled mode dropped** → for that case use `NewmarkExplicit 0.5`. | shipped (per-TU verified; full link pending the Ladruno blocker) |
| **EnergyBalance recorder** | `recorder EnergyBalance -file f [-time] [-region tag ...] [-csv\|-xml\|-binary\|-tcp addr port] [-precision N] [-scientific] [-closeOnWrite]` | recorder emit; per-region columns `KE IE DW ULW RES ERR`. **Text sidecar**, *not* MPCO — route to a text reader. The same energy math now also lands inside `.ladruno` (see MPCO row, ADR D8), so prefer that when MPCO is already on. | shipped |
| **Ladruno** | `recorder ladruno file -N <nodal...> -E <element...> [-G energy <regionTag...>] -T dt -R region` | `_response_catalog` `IntRule` enum **must match** `ladruno::detail::ElementIntegrationRuleType`; `Results.from_ladruno` reader keyed on `.ladruno` `FORMAT_VERSION`; partition discovery regex `^(?P<stem>.+?)\.part-(?P<idx>\d+)\.ladruno$`; consumes `basisInfo`/`QUADRATURE` self-declaration (new conforming elements need *no* recorder edit). See the schema notes below. | shipped (schema actively evolving — pin to `FORMAT_VERSION`) |
| **Ladruno brick element(s)** | TBD — higher-order hex, sibling to BezierTri6 on the solid side; will self-declare via the element contract | not yet — heads-up for the apeGmsh element registry. | draft (no plan file yet) |
| **OpenSeesPyMP** | `import openseesmp` (per-rank MPI Python module) | affects which engine `Results.run()` drives (`openseespy` vs `openseesmp`). | shipped |

## Ladruno — reader notes

The `.ladruno` schema is **actively evolving**; treat [[ladruno_schema_v1]]
and the Ladruno-recorder row in [[LEDGER_implementations]] as authoritative, and
**pin the reader to `FORMAT_VERSION`**. Current high points apeGmsh can rely on / watch:

- **File identity:** key the reader on `GENERATOR="Ladruno"` + `FORMAT_VERSION=1`.
  The recorder rename changed `GENERATOR` from `"MPCO_Ladruno"` to `"Ladruno"` — a
  reader matching the old string will reject current files.
- **Time-series is chunked** `[T × nIds × nComp]` (ADR D3) — one growing
  chunked+deflate dataset per result with `STEP[T]`/`TIME[T]` axes, *replacing* the
  old per-step `DATA/STEP_<k>`. The reader must handle chunked **and** legacy.
- **`MODEL/LOCAL_AXES` is now written** (per-class `{ID, FRAME[E×4 quaternion]}`
  from the element `"localAxes"` response) — ElasticBeam3d wired, other beams to
  follow. No silent identity fallback. (This was the old #1 gap; now landing.)
- **Standard-rule `QUADRATURE` by derivation** + `GLOBAL_GP_COORDS` belt-and-
  suspenders, explicit `NDIR`; per-stage `KIND` (`-kind transient|static|eigen`);
  whole-model + per-region energy (`RESULTS/ON_DOMAIN/energyBalance`,
  `ON_REGIONS`, `MODEL/SETS/SET_<tag>`).
- **Still single-stage-safe only where noted** — confirm multi-stage status in the
  recorder's schema/contract docs before relying on more than the first `MODEL_STAGE`.

## Implementation notes & recommended apeGmsh approach (per feature)

For each feature: **what's relevant on the OpenSees/implementation side**, and
**what we think is a good idea** for apeGmsh to do with it. Grounded in apeGmsh's
current seams — the `apeSees` bridge + typed namespaces, `_element_capabilities.py`
(`_ELEM_REGISTRY`/`_ElemSpec`), `_response_catalog.py`, `Results.from_native`/
`from_mpco`, `results.elements.line_stations`, `results.plot.line_force`, the
neutral `_vocabulary`, and per-zone schema versioning (ADR 0023). These are
recommendations; apeGmsh owns its own roadmap.

### BezierTri6 (element)

**OpenSees side.** 6-node quadratic triangle on Gmsh `tri6` (etype 9). The Gmsh
node order is **byte-identical** to the element's control-point order — emit
connectivity verbatim, `node_reorder=identity` (verified `1.78e-15`). It
self-declares geometry to the Ladruno recorder via the `basisInfo`/`integrationPoints`/
`integrationWeights` probes (`FAMILY=bernstein`, `PARAM_DOMAIN=bary`,
`rational=0`), so the recorder needs zero edits. Straight-sided only (v1); `-bbar`
valid PlaneStrain/3D only (warns + ignored for PlaneStress).

**Recommended apeGmsh approach.**
- Add a `_ElemSpec` in `_element_capabilities.py`: `gmsh_etype=9`,
  `node_reorder=identity`, `slots=(nodes, thick, type, matTag)`,
  `mat_family=nDMaterial`, `has_gauss=True`, `needs_transf=False`; plus a typed
  `ops.element.BezierTri6(pg=, thick=, type=, material=, bbar=…)`.
- In `_response_catalog.py`, register the element tag (read **live** from
  `classTags.h` — ≥33000 band, don't hardcode the dead 272) and its barycentric
  Gauss natural coords so result reads map correctly.
- **Direct-drive works today** (mesh in apeGmsh py3.11 → run in the fork build
  py3.12); the typed primitive is the ergonomic upgrade. Keep direct-drive
  documented as the fallback ([[bezier_apegmsh_integration]]).

### ExplicitBathe / ExplicitBatheLNVD (integrators)

**OpenSees side.** `integrator ExplicitBathe p` / `ExplicitBatheLNVD p [alpha]`,
`p` first positional `∈(0,1)` default `0.54`, LNVD `alpha∈(0,1)` default `0.8`.
Explicit, conditionally stable. Exposes `criticalTimeStep()` → per-model `dt_cr`
(valid after first `analyze`/`domainChanged`); damping shrinks `dt_cr`; needs a
lumped mass (`-lump rowsum|diagonal`) for true explicit behavior.

**Recommended apeGmsh approach.**
- apeGmsh's bridge is model-building-centric; analysis config is thin
  (`ops.analyze`, `/opensees/analysis` zone). Add an **integrator/analysis
  declaration surface** so the deck emits e.g. `integrator ExplicitBathe 0.54
  -lump rowsum`.
- **Wire `criticalTimeStep()` into an auto-dt helper** — query `dt_cr` after build,
  run at `dt = safety·dt_cr`. This is the single biggest explicit-run ergonomic
  win; users shouldn't hand-guess `dt`.
- Ensure apeGmsh emits a lumped mass for explicit runs (element `-cMass` off /
  integrator `-lump`), consistently.

### CentralDifferenceLadruno (integrator)

**OpenSees side.** `integrator CentralDifferenceLadruno [...]` — leap-frog CD done
right: correct first-step starter (`v_{n-1/2}=v0−½·dt·a0`), built-in `dt_cr`, clean
full-step velocity, βK guard. **Coupled mode was dropped** → for that case use
`NewmarkExplicit 0.5`. Shipped (PR #22), per-TU verified (full link pending the
Ladruno recorder blocker).

**Recommended apeGmsh approach.**
- Same integrator/analysis surface + `criticalTimeStep()` auto-dt as the Bathe
  integrators. In apeGmsh's explicit helpers, offer CentralDifferenceLadruno for
  pure central-difference and ExplicitBathe for sub-stepped dissipation; **do not
  expose the dropped coupled mode** — route that to `NewmarkExplicit 0.5`.

### EnergyBalance recorder

**OpenSees side.** `recorder EnergyBalance -file f [-region tag…]` → **text
sidecar**, per-region columns `KE IE DW ULW RES ERR`. The energy math is lifted
into the shared `EnergyBalanceKernel.h`, which now **also** feeds the Ladruno
recorder's energy result type: `RESULTS/ON_DOMAIN/energyBalance` + per-region
`ON_REGIONS` (via the `-G energy <regionTag…>` verb), so the same energy lands
**inside** `.ladruno`.

**Recommended apeGmsh approach.**
- **Prefer the `.ladruno`-embedded energy channel** over the text sidecar when the
  Ladruno recorder is already on — one lineage-tracked file, viewer-readable. Add a
  `Results` energy accessor (e.g. `r.energy(region=…)` reading `ON_DOMAIN/energyBalance`
  + `ON_REGIONS`).
- For the standalone recorder, a small **text → DataFrame** reader (cols
  `KE/IE/DW/ULW/RES/ERR`) is enough; route it via a `text` file-format
  discriminator, not the HDF5 reader.
- Surface **`ERR` (energy-balance error %)** as a solution-quality indicator in
  history plots / the viewer — a cheap, high-value diagnostic for explicit runs.

### Ladruno recorder (`.ladruno`)

**OpenSees side.** `recorder ladruno file -N <nodal…> -E <element…> [-G energy
<regions…>] -T dt -R region`. **Self-describing** — elements declare basis/
quadrature, so new conforming elements record with zero recorder edits. Current
schema: chunked `[T×nIds×nComp]` time-series (replaces per-step `DATA/STEP_<k>`);
**`MODEL/LOCAL_AXES` now written** (per-class quaternion `FRAME` from the element
`localAxes` response — ElasticBeam3d wired, more beams to follow); per-stage `KIND`
(transient/static/eigen); standard-rule `QUADRATURE` by derivation +
`GLOBAL_GP_COORDS` + explicit `NDIR`; energy `ON_DOMAIN`/`ON_REGIONS` +
`MODEL/SETS`. `IntRule` enum must match `ladruno::detail::ElementIntegrationRuleType`.
Partitions: `stem.part-N.ladruno` (0-based).

**Recommended apeGmsh approach.**
- Add a typed `ops.recorder.Ladruno(…)` (apeGmsh-side name; sibling of the existing
  `ops.recorder.MPCO`) and a `Results.from_ladruno(path, *, model=…)` reader
  (sibling of `from_mpco`), keyed on `GENERATOR="Ladruno"` + `FORMAT_VERSION=1` with
  a two-version window (mirror ADR 0023) — note the `GENERATOR` string is
  `"Ladruno"`, not the old `"MPCO_Ladruno"`. Reuse the partition-merge logic (regex
  `^(?P<stem>.+?)\.part-(?P<idx>\d+)\.ladruno$`); handle chunked **and** legacy
  time-series.
- **The big one — beam orientation is now solvable from `.ladruno` directly.** The
  apegmsh-helper skill's §7.2 (*"MPCO carries no beam vecxz… don't read MPCO
  LOCAL_AXES"*) is **out of date for the Ladruno recorder**: it writes `MODEL/LOCAL_AXES`
  as a per-class quaternion `FRAME` for wired beams. So apeGmsh's `line_force` /
  section diagrams can orient straight from `.ladruno` (for wired classes) instead
  of forcing the native `/opensees/transforms` path or the deferred injection seam.
  **Recommend:** teach the `.ladruno` reader to populate
  `results.elements.line_stations` orientation from the quaternion `FRAME`; fall
  back to native `vecxz` only when `LOCAL_AXES` is absent (unwired classes). Update
  the skill note once the reader lands.
- Map the Ladruno recorder's component names onto apeGmsh's neutral `_vocabulary`
  (`axial_force`, `shear_y/z`, `torsion`, `bending_moment_y/z` + conjugate strains)
  so `line_force(…)` works without per-feature shims.
- Read the `IntRule` enum values **from the fork** (the contract point with
  `ladruno::detail::ElementIntegrationRuleType`) — don't hardcode. Pin the reader to
  `FORMAT_VERSION` and treat the schema as actively evolving.

## Related docs (the normative detail)

This is a quick reference; the deep specs live next door:

- [[ladruno_element_contract]] — the element-side `setResponse` grammar +
  normative ordering conventions.
- [[ladruno_schema_v1]] — the on-disk `.ladruno` HDF5 layout.
- [[bezier_apegmsh_integration]] — BezierTri6 direct-drive round-trip + the
  two-environment (py3.11 apeGmsh / py3.12 fork build) split.
- [[LEDGER_implementations]] — authoritative class tags + shipping PRs.

## Maintenance log

- 2026-05-30 — Created as a lean apeGmsh-facing reference. (An earlier draft
  proposed a machine-readable manifest + CI gate + a compiled `ladrunoCapabilities()`
  runtime probe + a cross-repo sync bot; dropped as over-engineered for a
  solo-owned pair of repos — the plain doc is the deliverable. Revisit a
  structured/automated contract only if a second maintainer joins, apeGmsh ships to
  external users, or apeGmsh consumes a machine-readable feature list in code.)
- 2026-05-31 — Rebased onto current `ladruno` and refreshed: class tags moved to
  the **private ≥33000 band** (old sub-300 values dead), CentralDifferenceLadruno
  shipped (PR #22), Ladruno gained chunked time-series, `MODEL/LOCAL_AXES`,
  per-stage `KIND`, and energy ON_DOMAIN/ON_REGIONS. Dropped the hardcoded numeric
  tag column in favor of pointing at the ledger — exactly the field that went stale.
- 2026-05-31 — Added the per-feature **Implementation notes & recommended apeGmsh
  approach** section (EnergyBalance, ExplicitBathe/LNVD, CentralDifferenceLadruno,
  BezierTri6, Ladruno recorder) — what's relevant on the OpenSees side + what we
  recommend apeGmsh do, grounded in the apegmsh-helper skill's seams. Added the
  `GENERATOR="Ladruno"` reader note. (These rode in after PR #27 merged + the
  recorder rename swept the doc; re-landed here against the renamed base.)
