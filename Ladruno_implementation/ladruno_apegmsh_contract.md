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
| **Analysis monitor** (live) | `recorder Monitor -node <n…> -dof <d…> -resp disp\|vel\|accel\|reaction -sink file.h5 [-every K] [-hz H] [-region tag]` → SWMR-HDF5 stream tailable *while the run is live*; stop via `remove recorder $tag` | **apeGmsh must implement a consumer** — a typed `ops.recorder.Monitor(…)` emitter + a live tailing reader (open `swmr=True`, `ds.id.refresh()`) feeding the viewer's live plots. File is self-describing (`COLUMNS` var-strings, `STEP`/`TIME`/`FRAMES`), `FORMAT="ladruno-monitor"` `FORMAT_VERSION=1`. Also a valid at-rest file (read post-run like any results h5). class tag `RECORDER_TAGS_LadrunoMonitorRecorder`=33002. | **shipped (v1, sequential, nodal scalars)** — see [[08_analysis_monitor]] | 
| **Embedded reinforcement** (`LadrunoEmbeddedRebar` + `LadrunoBondSlip`) | `element LadrunoEmbeddedRebar tag rebarNode {nHost h… \| -host eleTag} {-shape N… \| -xi ξ…} -dir dx dy [dz] (-perfect kAxial \| -bond matTag [-bondScale πd·Ltrib]) [-kt {k\|auto} -ktAlpha a] [-corot -xiB ξ…] [-enforce penalty\|al] [-bipenalty {-dtcr dt \| -wcap β}]`; axial law `uniaxialMaterial LadrunoBondSlip tag τmax s1 s2 s3 τf α [-Gf Gf] [-s0 s0]` | **apeGmsh owns the inverse map** (global bar point → host ξ) and should ship a `g.reinforce(host=<set>, bars=<set>, bond=…, enforce=…, explicit=…)` generator: lay out `corotTruss`/beam rebar along bar paths, locate the host + inverse-map each rebar node to ξ, emit `-host -xi` primitives (host owns the weights via `getInterpolationWeights`). **`-xi` is 3D-only** (LadrunoBrick/BezierTet10 override it); 2D hosts need apeGmsh-computed `-shape`. ELE **33005** / MAT **33002** (read live). | **shipped (element + §10 roadmap); apeGmsh `g.reinforce` generator TO IMPLEMENT** — full grammar/theory in [[LadrunoEmbeddedRebar_guide]] |
| **Absorbing boundaries** (`ASDAbsorbingBoundary2D/3D` + `LysmerTriangle`) | `element ASDAbsorbingBoundary2D tag n1..n4 G v rho thick btype [-fx tsx] [-fy tsy]`; `element ASDAbsorbingBoundary3D tag n1..n8 G v rho btype [-fx tsx] [-fy tsy] [-fz tsz]`; `element LysmerTriangle tag i j k rho Vp Vs [len] [stage]` | **apeGmsh owns the generator** — a `g.absorbing_boundary(faces=…, soil=…, rayleigh=…, base_input=…)` that extrudes a **conforming one-element-thick ghost layer** outward from the domain faces (inner face = existing soil nodes, outer face = new ghost nodes), emits **one** element per ghost brick with `btype` ORed from face membership (corners/edges combined → single element), sets `G,v,rho` (+`thick` 2D) from the adjacent soil, assigns **Rayleigh on the absorbing elements** (required — free-field column damping), and emits the staging hook `setParameter -val 1 -ele <tags> stage` between gravity and transient. Base input is a **velocity** TimeSeries (within/÷2 motion) on **bottom** elements only. Axis-aligned faces only. Upstream tags **185/219/220**. **No new results reader** (only scalar `eleResponse stage\|G\|v\|rho\|E`); read physics from adjacent soil nodes — just tag absorbers to exclude from contour plots. | **upstream shipped; apeGmsh `g.absorbing_boundary` generator + staging hook TO IMPLEMENT** — full theory/arch/contract in [[lysmer_asd_absorbing_boundaries_guide]] |
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
  from the element `"localAxes"` response) — Elastic/Force/Disp **and**
  Mixed/GradientInelastic beam-columns (2D+3D) are wired. No silent identity
  fallback. (This was the old #1 gap; now landing.)
- **Standard-rule `QUADRATURE` by derivation** + `GLOBAL_GP_COORDS` belt-and-
  suspenders, explicit `NDIR`; per-stage `KIND` (`-kind transient|static|eigen`);
  whole-model + per-region energy (`RESULTS/ON_DOMAIN/energyBalance`,
  `ON_REGIONS`, `MODEL/SETS/SET_<tag>`).
- **Still single-stage-safe only where noted** — confirm multi-stage status in the
  recorder's schema/contract docs before relying on more than the first `MODEL_STAGE`.

### 📣 2026-06 recorder hardening — what apeGmsh can rely on now

A round of `LadrunoRecorder` fixes landed on `ladruno` (PRs **#200** + **#201**,
both merged, CI-green). The apeGmsh-facing effects — **no command-grammar changes**,
all additive/robustness:

- **Shell per-layer stresses now record.** For a layered shell (`LayeredShell` /
  `MembranePlateFiberSection`), `-E material.fiber.<resp>` (e.g.
  `material.fiber.stress`) now emits a real per-layer `ON_ELEMENTS` bucket —
  it previously produced **no bucket at all**. Output is byte-identical to
  `section.fiber.<resp>` (the recorder swaps `section`→`material` for shells
  internally; both verbs are equivalent for shells now). Reminder: the `-E` verb is
  a **single dot-joined token**, never space-separated args.
  - **Caveat (pre-existing, unchanged):** the shell-layer `COLUMN_MAP` encodes each
    `(gauss-point, layer)` pair as a distinct running `GAUSS_ID`, with
    `FIBER_ID = -1`, `SECTION_TAG = -1`, and component names that may read
    `UnknownStress`. The **data is complete and correct** (`nGP × nLayer × nComp`
    columns); only the per-column *layer label* is implicit. Decode layers as
    `gp = GAUSS_ID // nLayer`, `layer = GAUSS_ID % nLayer`, or read `nLayer` from
    `MODEL/SECTION_ASSIGNMENTS`. (A richer per-column locator is deferred — it would
    perturb the frozen-MPCO parity oracle for marginal gain.)
- **Partition schema is uniform again under OpenSeesMP.** With `-envelope` and/or
  `-precision f32`, worker ranks now write the **same** schema/precision as rank 0
  (the recorder serializes `envelope_mode` + the f32 flag across `sendSelf`/
  `recvSelf`). Before this, `.part-0` could be ENVELOPES/f32 while `.part-N` were
  time-series/f64 — breaking the stitch-on-read. The stitcher can again assume one
  family + precision across the whole `.part-*` set.
- **`STORED_PRECISION` is honest in envelope mode.** An `-envelope` file now reports
  `f64` (envelope `MIN/MAX/ABSMAX` are always f64); `-precision f32` narrows only
  the **time-series** `DATA` datasets. The reader can trust the attribute per file
  to choose its diff tolerance.
- **The recorder never aborts the analysis.** An element whose `setResponse` emits
  an output-tag nesting the recorder doesn't recognize (or two same-classTag
  elements with differing node counts) used to `exit(-1)` the whole run; it now
  warns and **drops just that bucket**. Directly relevant to apeGmsh-driven runs
  with Bézier / custom elements — an unmappable result degrades gracefully instead
  of killing the job. A malformed bare `-E section` / `-E material` (keyword, no
  sub-verb) is likewise a no-op now, not a segfault.

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

### Analysis monitor (live `recorder Monitor`)  — **TO IMPLEMENT on apeGmsh side**

**OpenSees side.** A live streaming recorder (`08_analysis_monitor.md`). You pick
nodes/dofs/response, call `analyze(N)` once, and it streams selected nodal scalars
(disp/vel/accel/reaction) to a **SWMR-HDF5** file a separate process can tail *while
the analysis is still running* (Shape B — the engine pushes from `record()` at each
commit). Two rate gates bound the cost so it never perturbs the solve: `-every K`
(step decimation) + `-hz H` (wall-clock throttle). The append is fire-and-forget —
a slow/absent reader can never stall `analyze()`. Self-describing layout: root attrs
`FORMAT="ladruno-monitor"`/`FORMAT_VERSION=1`/`GENERATOR="Ladruno"`; `COLUMNS`
(variable-length strings, e.g. `node5.disp.dof1`); chunked/unlimited `STEP[T]`,
`TIME[T]`, `FRAMES[T×nCols]`. **One file, two lifetimes** — the live sink is also a
valid at-rest results file. v1 = nodal scalar time-histories, sequential; element/
region channels, a profiler-telemetry tee, and the parallel path are follow-ups.

> **`STEP` is a global monotonic id, not a 0-based per-analysis step.** OpenSees'
> `getCommitTag()` does not reset on `wipe()` — treat `STEP` as a monotonically
> increasing id, key plots off `TIME`/array index, and never assume `STEP[0]==0`.

**Recommended apeGmsh approach (this is the deliverable — the monitor is shipped,
apeGmsh's consumer is not).**
- **Emitter:** a typed `ops.recorder.Monitor(nodes=…, dof=…, resp=…, sink=…,
  every=…, hz=…)` (sibling of `ops.recorder.MPCO`/`.Ladruno`) that emits the
  `recorder Monitor …` line; gate at point-of-use with a "requires the Ladruno fork
  build" error on stock openseespy.
- **Live reader / viewer:** this is the actual UX — a tailing consumer that opens the
  sink `h5py.File(path, 'r', swmr=True)`, calls `ds.id.refresh()` on a timer (set
  `HDF5_USE_FILE_LOCKING=FALSE`), and pushes new `FRAMES` rows into the live plot for
  the picked channels. The eventual OpenSees-side path is a FastAPI `/live` SSE
  sidecar (Ladruno P1) → apeGmsh can either consume that SSE or tail the file
  directly; tailing the file directly is the simplest first cut.
- **At-rest reuse:** because the sink is a normal results h5, a post-run
  `Results.from_monitor(path)` (or fold into the existing reader) gives the same
  channels as a plain time-history — free, no live machinery.
- **Selection ↔ picker:** the `COLUMNS` labels are authoritative and self-describing;
  drive the viewer's channel picker straight off them (no per-feature shim).

### Embedded reinforcement (`LadrunoEmbeddedRebar` + `LadrunoBondSlip`) — **`g.reinforce` generator TO IMPLEMENT on apeGmsh side**

**OpenSees side.** A penalty **coupling** element (ELE **33005**) that ties one
discrete rebar node to a solid host element's nodes via shape-function weights, so
a rebar mesh embeds in a **non-matching** concrete mesh. **Translations-only** tie
(the bar's own axial stiffness lives on a separate `corotTruss`/beam). Axial law:
perfect bond (`-perfect kAxial`) or a `LadrunoBondSlip` τ–s law (MAT **33002**,
CEB-FIP/MC2010 backbone; `-bond matTag -bondScale πd·Ltrib`). Host-agnostic:
`-host eleTag -xi ξ η ζ` lets the host own the weights via `getInterpolationWeights`
(**3D hosts only today** — LadrunoBrick trilinear, BezierTet10 Bernstein), while
`-shape N…` accepts apeGmsh-supplied weights for any host. Plus `-kt auto`
(host-stiffness-scaled transverse penalty), `-corot` (objective under large host
rotation), `-enforce penalty|al` (augmented Lagrangian → near-exact bond at
moderate stiffness), and `-bipenalty -dtcr|-wcap` (explicit critical-time-step
control). **Full grammar, theory, responses, and use cases: [[LadrunoEmbeddedRebar_guide]].**

**Recommended apeGmsh approach (the deliverable — element shipped, generator not).**
- **Ship a `g.reinforce(host=<set>, bars=<set>, bond=…, enforce=…, explicit=…)`
  generator.** This is naturally apeGmsh-side: the irreducible step is **point
  location / inverse map** (global bar point → host natural coord ξ), which apeGmsh
  owns (it has the geometry + mesh). It should:
  1. lay out longitudinal/transverse rebar as `corotTruss` (default) or beam
     elements along the bar paths, carrying the steel material
     ([[LadrunoUniaxialJ2_guide]] / [[LadrunoRebarBuckling_guide]]);
  2. for each rebar node, find the containing host element and **inverse-map to ξ**
     with a *guarded* Newton (relative tolerance, det/condition guard, explicit
     out-of-bounds policy — ADR 20 D3);
  3. emit `element LadrunoEmbeddedRebar … -host <hostEle> -xi ξ η ζ -dir … (-perfect …
     | -bond <LadrunoBondSlip> -bondScale <πd_b·Ltrib>)` — prefer **`-xi`** (host owns
     the weights, single source of truth) over emitting `-shape`;
  4. for the `-bond` path, emit the `LadrunoBondSlip` material and compute
     `bondScale = πd_b · L_trib` from the bar geometry.
- **Explicit runs:** add `-bipenalty -dtcr <dt>` (or `-wcap β`) so the stiff penalty
  tie doesn't collapse the explicit step, then wire `ops.criticalTimeStep()` (now
  accounts for the ties via the §10.6.1 seam) into the same auto-dt helper as the
  explicit integrators above.
- **2D hosts:** no element implements `getInterpolationWeights` in 2D yet, so for a
  2D quad/tri host apeGmsh must compute the **`-shape`** weights itself (bilinear /
  area coords); 3D hosts (LadrunoBrick/BezierTet10) take `-xi`.
- **General node-on-continuum ties (mesh stitch, point/part embed, frame-into-solid):
  use the fork `LadrunoEmbeddedNode` (ELE 33006), NOT the bar-in-solid rebar element.**
  `LadrunoEmbeddedRebar` is the *anisotropic bar-in-solid* tool (axial bond + transverse
  penalty); the **isotropic** node tie is `LadrunoEmbeddedNode` — see its own section below.
  (`ASDEmbeddedNodeElement` stays a valid **2D-native fallback** where no fork solid host is
  in play, but it is implicit-only and tri/tet-only.)

### General node-to-host embedment (`LadrunoEmbeddedNode`) — **`g.embed` generator TO IMPLEMENT on apeGmsh side**

**OpenSees side.** A penalty **coupling** element (ELE **33006**) — the **isotropic** sibling
of `LadrunoEmbeddedRebar` over the same kernel — that ties one constrained node to a host
element's nodes via shape-function weights, so an arbitrary node embeds in a **non-matching**
host mesh (mesh stitch, point/part embed, beam/shell-into-solid, SSI). **v1 SUPPORTED CORE =
the U translational tie + `g0` stress-free birth (no jolt on staged addition) + penalty/AL/
bipenalty enforcement + `-kt auto` conditioning** — host-agnostic (hex/tet/quad), implicit +
explicit. **UR/UP/D9/`-corot` are EXPERIMENTAL (not validated) — keep them off in the
generator.** This is the **drop-in fork upgrade** of the `ASDEmbeddedNodeElement` ties apeGmsh
already emits for non-matching meshes. **Full grammar, theory, responses, use cases (incl. the
generator contract): [[LadrunoEmbeddedNode_guide]].**

**Recommended apeGmsh approach (element shipped, generator not).** Ship a
`g.embed(nodes=<set>, host=<set>, k="auto", enforce="penalty", explicit=None, staged=True)`
generator. apeGmsh owns the irreducible **point location / inverse map** (global point → host
ξ). It should:
1. for each constrained node, find the containing host element and **inverse-map to ξ** (the
   same guarded Newton as `g.reinforce`, ADR 20 D3);
2. emit `element LadrunoEmbeddedNode <tag> <cNode> -host <hostEle> -xi ξ η ζ -k auto -kAlpha 1e3
   [-enforce penalty|al] [-bipenalty -wcap β|-dtcr dt]` — prefer **`-xi`** (host owns the
   weights) on **3D** hosts; for **2D** hosts (no `getInterpolationWeights` override) compute
   the bilinear/area weights and emit **`-shape N…`**;
3. default to **U-only** (rotations free — Abaqus `*EMBEDDED REGION`-consistent) and **`g0` on**
   (stress-free birth on a staged add). Keep `-rot`/`-pressure`/`-matN`/`-corot` **off**
   (experimental, opt-in only with a "not validated" note);
4. **frame→solid moment connections:** emit a **string of ties along the embedded stub**
   (≥2 nodes spanning a lever arm) so the moment transfers as a force couple — do **not** use
   the single-point `-rot` spin tie (see [[LadrunoEmbeddedNode_guide]] §10.5).
- **Explicit runs:** add `-bipenalty -wcap β` (host-aware) so the stiff penalty tie doesn't
  collapse the step; `ops.criticalTimeStep()` accounts for the ties (host-reduced `μ`) via the
  same `getExplicitCriticalTimeStep` seam as the rebar.

### Absorbing boundaries (`ASDAbsorbingBoundary2D/3D` + `LysmerTriangle`) — **`g.absorbing_boundary` generator TO IMPLEMENT on apeGmsh side**

**OpenSees side.** Three **upstream** elements (tags 185/219/220). The two ASD
elements are a FLAC/PLAXIS-style **free-field boundary + compliant base**, *not* bare
dashpots: each element straddles the truncation interface carrying a soil-side (SS)
node column shared with the mesh and a free-field (FF) "ghost" column outside it.
In the absorbing stage it superposes a 1-D **free-field column** (K/M/Rayleigh-C on
the FF nodes), a **FF→soil traction** `σ·n` injecting the free-field motion, **Lysmer
dashpots on the relative velocity `V_ff−V_ss`** (so only the scattered field is
absorbed), and on bottom faces a **compliant base** `2ρV·v_incident`. It runs a
**one-way stage machine** (`Stage_StaticConstraint=0` penalty roller/tie for gravity
→ `setParameter stage=1` → `Stage_Absorbing`), saving the gravity reactions so the
static stress survives the switch. `LysmerTriangle` is the bare dashpot (clean
diagonal `C` → the explicit-friendly option). **Full theory/architecture/impl with
`file:line` citations: [[lysmer_asd_absorbing_boundaries_guide]].**

**Recommended apeGmsh approach (the deliverable — elements shipped, generator not).**
- **Ship `g.absorbing_boundary(faces=…, soil=…, thickness=…, rayleigh=(aM,bK),
  base_input={…})`.** The irreducible step is geometric and apeGmsh-native:
  **extrude a conforming one-element-thick ghost layer** outward from the chosen
  domain faces — inner face reuses the existing soil boundary nodes, outer face is
  new ghost nodes (connected only to the absorbing elements), ghost width matched to
  the adjacent soil element size (it sets the free-field column stiffness `∝ lx/ly`).
- **One element per ghost brick**, node order free (the element self-sorts by
  coordinate), `btype` ORed from which **domain** faces the brick abuts — **edges and
  corners get the combined string and a single element** (`"BL"`, `"LF"`, `"BLF"`…);
  the C++ does the tributary weighting, so never overlap two.
- **Axis-aligned faces only** (3D `computeNmatrix` asserts vertical normals are ±X/±Y).
  Set `G,v,rho` from the adjacent soil (`G=E/(2(1+v))`); 2D also needs `thickness`.
- **Assign Rayleigh damping to the absorbing elements** (emit `region … -ele <absorb>
  -rayleigh aM 0 bK 0`) — **required, not optional**: `addCff` reads the element's own
  `alphaM/betaK`; without it the free-field column is undamped and resonates.
- **Input is a *velocity* TimeSeries** (the within/÷2 motion) on `-fx/-fy/-fz` of the
  **bottom** elements only; sides free-field automatically. Alternative: a DRM ring
  just inside the layer (full 3-D field) — then no `-fx/-fy/-fz`. Not both.
- **Layered / nonlinear soil — resolve props PER GHOST CELL.** The element is
  **linear-elastic** (no constitutive model; the free-field column never yields), so:
  give each cell the **small-strain `G₀`** of its local layer (`V_s=√(G₀/ρ)`), **not**
  a global or secant value; for pressure-dependent soil (`PDMY`/`PM4Sand`) `G₀` is
  **depth-graded** → `soil=` should accept a scalar, a per-layer table, **or a
  callable `f(x,y,z)→(G0,v,rho)`** and resolve per cell. Keep nonlinearity in the
  **interior** (place the boundary far enough out that boundary soil stays low-strain
  — warn if a nonlinear/structure region is within ~1–2 wavelengths). The
  gravity-before-`stage`-switch ordering is **mandatory** for pressure-dependent soil
  (builds `K₀` confining stress under the stage-0 roller/tie). Tune the absorber
  Rayleigh to the **small-strain** damping. Full rationale: §9.6 of the guide.
- **Staging hook:** add an analysis-side `g.opensees.stage_absorbing()` that emits
  `setParameter -val 1 -ele <absorb_tags> stage` **between** the gravity `Static`
  solve (`loadConst`) and the `Transient` block. One-way; any other transition aborts.
- **Results:** no new reader needed — these expose only scalar `eleResponse
  stage|G|v|rho|E`. **Tag the absorbing elements** so they're excluded from
  stress/contour plots (meaningless there); read boundary reactions/energy from the
  adjacent soil nodes. A `stage` probe is a cheap "did the switch fire?" diagnostic.
- **Solver hint to surface:** the ASD tangent is **non-symmetric** (one-way FF→soil)
  and **penalty**-based; default the deck to `system UmfPack`/`FullGeneral` (not
  `BandSPD`) when an absorbing layer is present, and avoid extra `fix`/`equalDOF` on
  boundary DOFs.

## Related docs (the normative detail)

This is a quick reference; the deep specs live next door:

- [[ladruno_element_contract]] — the element-side `setResponse` grammar +
  normative ordering conventions.
- [[ladruno_schema_v1]] — the on-disk `.ladruno` HDF5 layout.
- [[bezier_apegmsh_integration]] — BezierTri6 direct-drive round-trip + the
  two-environment (py3.11 apeGmsh / py3.12 fork build) split.
- [[LadrunoEmbeddedRebar_guide]] — embedded-reinforcement coupling element: full
  grammar, theory (tie, anisotropic traction, penalty/AL, co-rot, bipenalty), and
  the `g.reinforce` use cases.
- [[LadrunoEmbeddedNode_guide]] — general node-to-host embedment (ELE 33006): the
  isotropic U tie + `g0` stress-free birth + penalty/AL/bipenalty, the responses, and
  the `g.embed` generator contract (incl. frame→solid moment ties).
- [[20_ladruno_embedded_reinforcement_adr]] — embedded reinforcement ADR (D1–D6 +
  the §10 Mode-P roadmap; the guarded inverse-map contract is D3).
- [[lysmer_asd_absorbing_boundaries_guide]] — Lysmer/ASD free-field absorbing
  boundaries: full theory, OpenSees architecture/implementation, and the §9
  `g.absorbing_boundary` ghost-layer generator contract.
- [[LEDGER_implementations]] — authoritative class tags + shipping PRs.

## Maintenance log

- 2026-06-07 — Added **Absorbing boundaries** (`ASDAbsorbingBoundary2D/3D` +
  `LysmerTriangle`, upstream tags 219/220/185). Elements are shipped upstream; the
  apeGmsh **`g.absorbing_boundary` ghost-layer generator** + the `stage_absorbing()`
  staging hook are **TO IMPLEMENT**. Full deep dive (free-field boundary + compliant
  base, stage machine, 3D static condensation, ghost-layer contract) in the new
  [[lysmer_asd_absorbing_boundaries_guide]].

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
- 2026-05-31 — Added the **Analysis monitor** (live `recorder Monitor`, SWMR-HDF5,
  tag 33002) — shipped v1 on the OpenSees side, but its apeGmsh consumer (typed
  emitter + live tailing reader/viewer + at-rest `from_monitor`) is **not yet
  implemented and must be**. See [[08_analysis_monitor]].
- 2026-06-04 — Added **Embedded reinforcement** (`LadrunoEmbeddedRebar` ELE 33005 +
  `LadrunoBondSlip` MAT 33002): the element + bond-slip + full §10 roadmap (`-host`/
  `-xi`, `-kt auto`, `-corot`, `-enforce penalty|al`, `-bipenalty` + the `-cfl`
  seam) are **shipped**; the apeGmsh **`g.reinforce` generator** (inverse map → emit
  `-host -xi` primitives, `LadrunoBondSlip` + `bondScale`, `-bipenalty` for explicit)
  is **TO IMPLEMENT**. Full grammar/theory/use in the new [[LadrunoEmbeddedRebar_guide]].
  Recorded the routing rule: node-on-continuum ties (frame on a 2D boundary) go to
  upstream `ASDEmbeddedNodeElement`, not this bar-in-solid element.
