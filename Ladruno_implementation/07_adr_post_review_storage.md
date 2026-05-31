---
title: ADR — MPCO_Ladruno storage architecture, post adversarial review
project: Ladruno
status: accepted
date: 2026-05-30
supersedes-parts-of: 03_mpco_ladruno.md
tags:
  - adr
  - recorder
  - schema
  - decision
---

# ADR — MPCO_Ladruno storage architecture (post adversarial review)

## Context

The `mpcoLadruno` recorder (a bespoke HDF5 sibling of the frozen `MPCORecorder`,
ADR [[03_mpco_ladruno]]) reached working node+element value parity (1e-12, single
and multi-stage). Before investing further we ran a **6-lens adversarial review**
(schema/geometry, HDF5 perf, parallel/MPI, reader ergonomics, maintainability,
steelman-the-alternative) + a synthesis judge.

**Verdict: proceed-with-changes — do NOT pivot the format.** All six lenses high
confidence; the steelman lens (job: find a better off-the-shelf format) reached
*proceed* — no standard container (VTKHDF / XDMF / CGNS / exodus / ADIOS2 / Zarr)
models the FEM-result semantics we need: fiber **and** resultant sections, per-step
co-rotational element frames, and natural-coordinate Gauss data with weights. The
review also confirmed the `ResultSource`/`ResultSink` seam is genuinely modular
(fixes are localized, not schema-wide). The two most damning findings (fabricated
self-description; bucket-as-dataset) were already fixed on `ladruno` by the
BezierTri6 work ([#18](https://github.com/nmorabowen/OpenSees/pull/18)).

This ADR records the decisions taken in response. Full review artifact:
session task `wov6bf1q1.output`. The *implementation* lives in
[[06_quadrature_global_gp_plan.md]].

## Decisions

### D1 — Keep the bespoke HDF5 format (do not adopt an off-the-shelf container)
No standard format represents fiber+resultant sections, per-step co-rotational
frames, and natural-coord Gauss-with-weights together. VTKHDF higher-order is
still maturing (2025) and models cells+fields, not sections. **Status: confirmed.**
Consequence: we own the schema and (at least) the apeGmsh + STKO_to_python readers;
mitigate drift with a shared validator (D5).

### D2 — Belt-and-suspenders geometry: write static `GLOBAL_GP_COORDS`
In addition to the parametric `BASIS` + `QUADRATURE` descriptor, the recorder writes
a **static, C++-computed `GLOBAL_GP_COORDS[nElem × nGP × ndim]`** (reference config)
per element group. **Why:** it is a write-time oracle that catches basis/ordering
bugs at build time, eliminates two-reader basis-math drift, and lets legacy elements
be consumed with zero reader-side basis functions. **Consequence:** one extra static
dataset per group (written once at model time, no per-step cost); the parametric
descriptor is retained for deformed-config reconstruction (co-rotational/large-strain).
The reader treats `GLOBAL_GP_COORDS` as authoritative and falls back to `x(ξ)=ΣRᵢXᵢ`
only when it is absent.

### D3 — Chunked extensible time-series instead of one dataset per step  ✅ IMPLEMENTED
Replace per-step `RESULTS/.../DATA/STEP_<k>` datasets with **one chunked, extensible
`[T × nIds × nComp]` dataset per result** + `TIME[T]`/`STEP[T]` vectors, with
chunk+shuffle+deflate. **Why:** per-step datasets explode the dataset count (≈10⁵–10⁶
tiny datasets on long transients), carry per-dataset metadata that can exceed the
payload for small-result/high-step runs, defeat compression along the (compressible)
time axis, and make a single time-history an O(T) dataset-open. **Consequence:**
smaller files + O(1) slice reads; localized to `StreamingSink::begin/accept` behind
the sink seam; the apeGmsh reader co-versions its `STEP_<k>` loop to a hyperslab.
**Do it now**, while only the in-repo parity harness + co-developed reader depend on
`STEP_<k>` — the migration only gets more expensive once a consumer ships against it.

> **Implemented (PR #36):** `MPCOL_Hdf5.h` gained `createTimeSeries3d`/`appendSlab3d`
> + extensible 1-D `TIME`/`STEP` helpers; `StreamingSink::begin` creates the
> `[0×nIds×nComp]` dataset + axes, `accept` appends one slab/step. All families flow
> through the one sink seam. Reader (`ladruno_format.iter_step_slices`) reads chunked
> *or* legacy per-step transparently, so a chunked `.ladruno` still diffs against the
> per-step frozen `.mpco`. Verified: values byte-identical (parity 1e-12), full
> regression green (80/80·96/96·144/144·72/72·108/108·standard_quad·pytest 10/10).
> Deferred within D3: the opt-in per-step `LOCAL_AXES` frame series adopts the same
> chunked shape when LOCAL_AXES lands.

### D4 — Explicit `NDIR` attribute, decoupled from `len(ORDER)`
Each element group carries an `NDIR` int = number of `GP_PARAM` columns = parametric
dimension (line 1, quad/tri 2, hex/tet 3). **Why:** simplices (`PARAM_DOMAIN="bary"`)
have a parametric dimension that differs from the polynomial `ORDER`; the prior
validator equated `ndir = len(ORDER)`, which is a tensor-product assumption that
breaks for tri/tet. **Consequence:** validator and readers read `NDIR`, not
`len(ORDER)`; `ORDER` keeps one entry per parametric direction (replicated total-degree
for simplices), so no format churn on the landed BezierTri6 group.

### D5 — Shared validator + write-time round-trip oracle; phased `FORMAT_VERSION=1`
`validate_ladruno.py` becomes a **shared importable module both readers gate on in
CI**, plus a mandatory write→read→reconstruct-`x(ξ)`→compare-to-`GLOBAL_GP_COORDS`
oracle (≤1e-12). **Why:** with two readers and no external validation surface, a writer
bug both reader and writer share would be invisible; the oracle is the only mechanism
that catches silent basis/ordering errors. `FORMAT_VERSION=1` is **held** until the
remaining Phase-0 items land (standard-rule QUADRATURE, `KIND`, `LOCAL_AXES`, `NDIR`,
shared validator) — nothing external consumes v1 yet, so there is no break to manage.

## Remaining (tracked, not decided away)
Done since this ADR: standard-rule QUADRATURE + GLOBAL_GP_COORDS (incl. higher-order),
chunked layout (D3 ✅), **`KIND` ✅** (`-kind` + eigen auto, PR #38), **`LOCAL_AXES` ✅**
(per-class quaternion FRAME from the element `"localAxes"` response; ElasticBeam3d wired,
PR #38). Still open: extend `"localAxes"` to the remaining beams (dispBeam/forceBeam/
ElasticBeam2d — identical `getLocalAxes` pattern); **shared validator + CI round-trip
oracle (D5)** then **freeze `FORMAT_VERSION=1`**; per-class result column-naming still in
the reader catalog (make `COLUMN_MAP/COMP_NAMES` authoritative); **parallel path**
unimplemented (`sendSelf`/`recvSelf` stubs → N ranks would race one file under MP);
envelopes (`EnvelopeSink` built, unused). Sequencing: → remaining-beam localAxes →
shared validator (D5) + freeze v1 → parallel → envelopes.
