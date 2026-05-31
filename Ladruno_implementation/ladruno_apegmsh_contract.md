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
and the MPCO row in [[LEDGER_implementations]] as authoritative, and **pin the
reader to `FORMAT_VERSION`**. Current high points apeGmsh can rely on / watch:

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
  MPCO docs before relying on more than the first `MODEL_STAGE`.

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
