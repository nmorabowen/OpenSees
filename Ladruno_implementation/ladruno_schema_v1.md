---
title: Ladruno HDF5 schema v1
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - recorder
  - schema
  - spec
---

# Ladruno HDF5 schema v1

The on-disk contract for the `ladruno` recorder. apeGmsh is the **canonical
consumer** and co-owns this spec; `STKO_to_python` gets a parallel reader branch.
Companion: the ADR [[03_ladruno_recorder]]. This file is the layout `ladruno::` writes and
apeGmsh reads.

> **Why diverge from STKO's `.mpco`?** Two planned element families — *second-order /
> Bézier* elements and *Belytschko co-rotational beams* — break the STKO schema's
> core assumption that the reader can reconstruct geometry by pushing parent
> coordinates through *hardcoded* low-order Lagrange shape functions. v1 replaces the
> two closed enums (`ElementGeometryType`, `ElementIntegrationRuleType`) + implicit
> reader basis with an **explicit, self-describing `BASIS` + `QUADRATURE` descriptor**,
> and **inverts the dispatch** so elements declare their own basis/quadrature.

## File identity

- **Extension: `.ladruno`** — the canonical extension for *all* Ladruno recorders
  (suffixable, e.g. `model.ladruno`, `model.something.ladruno`). A distinct extension
  keeps these files from being accidentally opened in STKO desktop (they are not
  STKO-readable by design, D2). apeGmsh's discovery glob targets `.ladruno`.
- **Partitioned runs (D5):** `<stem>.part-<N>.ladruno`, **N 0-indexed and contiguous
  from 0**. This is a hard interop contract — apeGmsh discovers siblings via the regex
  `^(?P<stem>.+?)\.part-(?P<idx>\d+)\.ladruno$` and *errors* on a gap. When stitching,
  a node present in several parts is combined per its result's `PARTITION_REDUCTION`
  (§7.1): kept once for `NONE`, **summed** for `SUM` (e.g. the reaction at a support on
  a partition interface — each part holds only its own elements' share; WP-126).

## Design principles

0. **Self-sufficient (canonical).** A `.ladruno` file carries **everything** apeGmsh
   needs — geometry, regions, beam local axes, section index, GP coords+weights,
   provenance — so `Results.from_ladruno` equals the native `model.h5` path with **no
   sibling file**. There is no privileged native path; this *is* the native path.
1. **Self-describing.** Every element group carries enough to reconstruct its
   geometric map and quadrature with **no element-class-specific reader logic**:
   `x(ξ) = Σ Rᵢ(ξ)·Xᵢ`, `Rᵢ` = the declared basis (rational if NURBS).
2. **Inverted dispatch (D-schema-1).** Elements declare `BASIS`/`QUADRATURE` via
   `setResponse` (generalizing the existing `integrationPoints` / IGA hooks). The
   recorder consults a legacy class-tag table **only as a fallback** for elements that
   predate the protocol. New elements never edit the recorder.
3. **Versioned.** `INFO/GENERATOR="Ladruno"` + `INFO/FORMAT_VERSION=1`. Readers
   branch on these; the schema may evolve without silent misreads.
4. **Parity by derivation.** For the ~60 legacy OpenSees elements, the recorder
   *derives* `BASIS`/`QUADRATURE` from the old `(GEOMETRY, INTEGRATION_RULE)` pair, so
   v1 needs zero element edits. Self-declaration is what *new* elements opt into.
5. **Keep what's good.** `MODEL_STAGE` staging, natural-coordinate Gauss points,
   fiber sections in section-local axes, per-result self-description — all retained.
   Fix what's clunky: store quadrature **weights**; replace the flattened `META`
   `;`-string with a structured column map.

---

## 1. Top-level tree

```
/
├── INFO
│   ├── GENERATOR        (string)  "Ladruno"
│   ├── FORMAT_VERSION   (int)     1
│   ├── FORMAT_MINOR_VERSION (int) 1          ← absent = 0 (§2a)
│   ├── SOLVER_NAME      (string)  "OpenSees"
│   ├── SOLVER_VERSION   (int[3])
│   ├── SPATIAL_DIM      (int)     1 | 2 | 3
│   └── PROVENANCE       (group)   empty from the recorder; consumer attributes, §7b.1
│
└── MODEL_STAGE[<id>]                            (one per analysis stage, §2a)
    ├── attrs STEP(int), TIME(double), KIND(string "transient"|"static"|"eigen")  ← per stage, §2a
    ├── MODEL
    │   ├── NODES         { ID[nN], COORDINATES[nN×ndim] }
    │   ├── ELEMENTS      { <group> … }           ← see §3 (BASIS/QUADRATURE)
    │   ├── LOCAL_AXES    { <classTag>-<ClassName>[<g>]/{ FRAME[nE×4], ID[nE] } }  ← §5
    │   ├── SECTION_ASSIGNMENTS  { SECTION_<tag>[<Class>] … }   ← fiber + resultant (§6)
    │   └── SETS          { SET_<regionTag> { NODES, ELEMENTS } }   ← self-contained regions
    └── RESULTS
        ├── ON_NODES      { <RESULT> … }          ← §7.1
        ├── ON_ELEMENTS   { <result> … }          ← §7.2 (structured column map)
        ├── ON_DOMAIN     { <scalar> … }          ← §7.3  (DEFERRED — energy work)
        └── ENVELOPES     { ON_{NODES,ELEMENTS,DOMAIN}/<name> … }   ← §7.4 (v3)
```

---

## 2. `INFO`

| Attribute | Type | Value | Notes |
|---|---|---|---|
| `GENERATOR` | string | `"Ladruno"` | **reader branch key** — distinguishes from STKO `.mpco` |
| `FORMAT_VERSION` | int | `1` | bump on any breaking layout change |
| `FORMAT_MINOR_VERSION` | int | `1` | apeGmsh #2034: a semantic, layout-compatible revision. Absent means 0. Minor 1 = per-stage `KIND` from the running analysis, and the stage boundaries of §2a |
| `SOLVER_NAME` | string | `"OpenSees"` | |
| `SOLVER_VERSION` | int[3] | e.g. `[3,5,1]` | `OPS_VERSION` split on `.` |
| `SPATIAL_DIM` | int | 1/2/3 | from first node, enforced uniform |
| `PARTITIONED` | int | 0/1 | 1 if this is one `.part-N` of a partitioned run |
| `PARTITION_ID` | int | `N` | this file's 0-based partition index (0 if not partitioned) |
| `NUM_PARTITIONS` | int | total | partition count of the set (1 if not partitioned) |
| `RUN_ID` | string | e.g. `slurm-81234.0` | WP-165 (MP-9): identity of the run that wrote this file; identical across the part files of one run unless `RUN_ID_SCOPE == "process"` |
| `RUN_ID_SCOPE` | string | `broadcast` / `user` / `launcher` / `process` | where `RUN_ID` came from: P0 broadcast (PartitionedDomain), `LADRUNO_RUN_ID`, the launcher's job id (SLURM / OpenMPI), or process-unique (no shared identity available — readers must not cross-check it) |

`MODEL_STAGE[<stamp>]` may carry `EMPTY_PARTITION = 1` (WP-165, MP-8): this process had
no nodes of the recorded set (e.g. a `-R` region that lives on other ranks); its
`MODEL/NODES/ID` and `COORDINATES` are zero-length and it holds no node/element results.

### 2a. `MODEL_STAGE` boundaries and `KIND` (`FORMAT_MINOR_VERSION = 1`, apeGmsh #2034)

A recorder opens a new `MODEL_STAGE` at a recorded step when any of these holds since
the open stage began:

1. the node / element / pressure-constraint set changed (WP-165 topology fingerprint);
2. the user ran `domainChange` (Tcl or openseespy): a **stage barrier**. Internal
   re-stamps of the domain-change stamp (an SP pattern, `eleLoad`, contact `-reemit`)
   are NOT barriers and stay folded into the open stage;
3. the running analysis flipped between static and transient, whatever the `KIND`
   label (`-kind` and eigen promotion set the label only), so no stage mixes static
   and transient steps.

The run's first recorded step opens the first stage; a barrier before it opens nothing.
A repeated kind is never a boundary: static → `wipeAnalysis` → static is one stage
unless a `domainChange` separates them.

`KIND` is decided per stage, at the stage's first record:

| Precedence | Source | `KIND` |
|---|---|---|
| 1 | a modal request (`-N modesOfVibration[Rotational]`) | `"eigen"` |
| 2 | `-kind <transient\|static\|eigen>` | that value for every stage (a label only: barriers and flips still split) |
| 3 | the running analysis: `analysis Static` | `"static"` |
|   | `analysis Transient` / `VariableTransient` / `PFEM` | `"transient"` |
| 4 | no analysis visible yet (an explicit `record` first, or an `eigen …; record` snapshot) | provisional `"static"`, rewritten by the first recorded step that sees an analysis; a warning if the stage closes still unresolved |

The transient analysis the classic Tcl `eigen` command builds when none is defined is
not a running analysis (openseespy deletes it after the solve; classic Tcl keeps it)
until an `analysis` or `analyze` command, so `eigen …; record; … analysis Static` is
one static stage in both interpreters.

The `<id>` in `MODEL_STAGE[<id>]` is the domain-change stamp when the stamp moved past the
open stage's id, else open id + 1 (a flip or barrier can arrive without a stamp move).
Ids are strictly increasing and never reused; readers sort them numerically and must
not assume they equal `Domain::hasDomainChanged()`.

Minor 0 files (no `FORMAT_MINOR_VERSION`): `KIND` was the recorder-wide `-kind`
(default `"static"`, so a transient stage read `static`), and from WP-165 until minor 1
a stage opened only on a topology change, so a same-topology static → transient run
was one stage.

Parallel: on openseesmp every rank runs the same commands, so barriers and flips agree
across `part-N` files. An OpenSeesSP (PartitionedDomain) worker runs no interpreter: it
sees neither the barrier nor the analysis type, so its part file follows topology-only
stages with a provisional `KIND` and warns (WP-165 doc §4).
A new `MODEL_STAGE` starts only when the node/element/pressure-constraint SET changes
(WP-165, R6) — not on every domain-change stamp move (SP patterns, `eleLoad`, contact
re-emit continue the current stage).

A reader that doesn't find `GENERATOR=="Ladruno"` must treat the file as a
legacy STKO `.mpco` (or refuse it).

**Partition manifest.** When the recorder runs under a parallel launcher it writes
`<stem>.part-<N>.ladruno` (one file per process) and stamps each with the three
`PARTITION*` attributes so a reader can validate the contiguous-from-0 set it globs.
The partition index `N` is resolved without MPI (the recorder is in the shared
sequentially-compiled `OPS_Recorder` library): from the `sendSelf`/`recvSelf`
broadcast counter in the PartitionedDomain/OpenSeesMP path, or from the launcher's
`PMI_RANK`/`PMI_SIZE` env in the interpreter-per-rank `openseesmp` path. Single
process → `PARTITIONED=0`, filename verbatim.

---

## 3. `MODEL/ELEMENTS` — the centerpiece

Elements are grouped by **identical `(classTag, BASIS, QUADRATURE-in-parametric-space)`**.
One dataset per group; the connectivity rows differ per element, the parametric basis
and quadrature are shared (that's the grouping criterion).

```
MODEL/ELEMENTS/<classTag>-<ClassName>[<g>]
   dataset:  CONNECTIVITY  [nElem × (1 + NUM_CTRL)]  int   rows: (elemTag, c1..cK)
   attrs (BASIS descriptor):
     TOPOLOGY        string   "line"|"tri"|"quad"|"tet"|"hex"|"wedge"|"pyramid"|"custom"
     FAMILY          string   "lagrange"|"serendipity"|"bernstein"|"nurbs"|"custom"
     ORDER           int[ndir]   polynomial order per parametric direction (e.g. [2,2])
     PARAM_DOMAIN    string   "[-1,1]" | "[0,1]" | "bary"
     RATIONAL        int      0|1   (1 ⇒ CTRL_WEIGHT present)
     NUM_CTRL        int      control points / nodes per element (== K above)
     NUM_GP          int      Gauss points per element
   child group QUADRATURE/
     GP_PARAM        [NUM_GP × ndir]  double   parametric coords, in PARAM_DOMAIN
     GP_WEIGHT       [NUM_GP × 1]     double   integration weights        ← NEW vs STKO
   sibling dataset (rational only):
     CTRL_WEIGHT     [nElem × NUM_CTRL] double  rational weights wᵢ, row-aligned to CONNECTIVITY
```

`<g>` is a 0-based disambiguator across basis/quadrature variants of the same class.
The `BASIS`/`QUADRATURE` block is the **complete** geometric self-description — there
is no reader-side per-class shape-function table.

### 3.1 The self-declaration protocol (D-schema-1)

> The exact element-side `setResponse` signatures (attr names, response types, C++
> skeletons) are specified in [[ladruno_element_contract]]. This section is the
> recorder-side view of the same protocol.

The recorder builds the descriptor for each element with this precedence:

1. **`setResponse("basisInfo", …)`** — element returns FAMILY / TOPOLOGY / ORDER /
   PARAM_DOMAIN / RATIONAL. *(new canonical keyword; generalizes IGA `IGAOrder`.)*
2. **`setResponse("integrationPoints", …)`** → `GP_PARAM`  *(already implemented by
   every `ForceBeamColumn*`; IGA exposes `IGAKnot1P/2P/3P`).*
3. **`setResponse("integrationWeights", …)`** → `GP_WEIGHT`  *(new; IGA exposes
   `IGAWeight`).*
4. **`setResponse("controlPointWeights", …)`** → `CTRL_WEIGHT`  *(rational only).*
5. **Legacy fallback** — if `basisInfo` is absent, derive `BASIS`/`QUADRATURE` from the
   old `getGeometryAndIntRuleByClassTag` table (Lagrange assumption, standard-rule
   abscissae/weights from the rule enum). This covers all ~60 existing elements for v1
   parity with **no element edits**.

**Consequence:** adding a Bézier or Belytschko element = implement `basisInfo`
(+ weights for Bézier). The recorder is never touched. This is the maintainability win.

### 3.2 Legacy → BASIS derivation (v1 parity table, abbreviated)

| Old GEOMETRY | Old INTEGRATION_RULE | Derived BASIS | Derived QUADRATURE |
|---|---|---|---|
| `Quadrilateral_4N` | `Quad_GL_2` (2×2) | quad / lagrange / [1,1] / [-1,1] | GP_PARAM = ±1/√3 tensor; GP_WEIGHT = 1 |
| `Hexahedron_8N` | `Hex_GL_2` | hex / lagrange / [1,1,1] / [-1,1] | 8-pt tensor GL |
| `Line_2N` | `Line_GL_1` | line / lagrange / [1] / [-1,1] | 1-pt, w=2 |
| `Line_2N` | `Custom` (FB beam) | line / lagrange / [1] / [-1,1] | GP_PARAM from `integrationPoints`, GP_WEIGHT from rule |
| `Triangle_3N` | `Tri_GL_1` | tri / lagrange / [1] / bary | centroid, w=1/2 |

(The full table mirrors `integration-rules-and-gauss.md` §1–5; the recorder fills
weights for standard rules from the known abscissae.)

### 3.3 Reader reconstruction contract (universal)

```
# global coordinate of Gauss point j of element e:
ξ   = QUADRATURE/GP_PARAM[j]                       # parametric
X   = COORDINATES[ CONNECTIVITY[e, 1:] ]           # control points / nodes
if RATIONAL:
    w = CTRL_WEIGHT[e]
    R = (w * B(ξ; FAMILY, ORDER)) ; R /= R.sum()   # rational basis
else:
    R = B(ξ; FAMILY, ORDER)                         # Bernstein/Lagrange/serendipity
x_global = R @ X
```

apeGmsh implements `B(ξ; FAMILY, ORDER)` **once per family**, not once per element
class. That is the whole point of the divergence.

### 3.4 Variable node count within a class — CSR connectivity (WP-226)

Some classes carry a **different node count per element**: `LadrunoDistributingCoupling`
(RBE3) and `LadrunoKinematicCoupling` (RBE2) have one footprint per element. A
rectangular `CONNECTIVITY` cannot hold them. Before WP-226 the recorder logged one error
per odd-sized element and left it out of the file. Such a group is now written as
CSR:

```
MODEL/ELEMENTS/<classTag>-<ClassName>[<g>]
   attr     CONNECTIVITY_LAYOUT = "csr"
   attr     NUM_CTRL            = -1          (varies per element)
   dataset  ID                    [nElem]          int   element tags
   dataset  CONNECTIVITY_OFFSETS  [nElem + 1]      int   0 = OFFSETS[0] <= … <= OFFSETS[nElem] = len(NODES)
   dataset  CONNECTIVITY_NODES    [OFFSETS[nElem]] int   nodes of element e = NODES[OFFSETS[e]:OFFSETS[e+1]]
   (no CONNECTIVITY, no GLOBAL_GP_COORDS, no CTRL_WEIGHT; BASIS attrs as §3)
```

- **Data-driven.** The writer switches to CSR only when the class's node counts differ
  in that file. A class whose elements all share one count (every RBE3 with 4
  independents, say) keeps the rectangular `CONNECTIVITY`, byte-identical to before.
  A reader must handle both.
- **Old readers skip, never misread.** The group has no `CONNECTIVITY`, and every
  known reader (`ladruno_format`, apeGmsh `FEMData` / results / tag map) guards on
  its presence. So a reader written before §3.4 drops these elements, as it did
  before, and never misreads a padded table. This is why there is no
  `FORMAT_MINOR_VERSION` bump: the change is additive (§9).
- Node order inside a row is the element's `getExternalNodes()` order. For RBE3/RBE2
  that is the reference node, then the independents or slaves.
- `RESULTS/ON_ELEMENTS` buckets are unchanged. They key by element `ID`, so a coupling
  that answers a request still pairs with its CSR row by tag.

---

## 4. `MODEL/NODES`

`ID [nN×1] int` (node tags, write order) + `COORDINATES [nN×ndim] double`. For
higher-order/IGA elements these are **control points** (geometry may not interpolate
them) — the `BASIS` descriptor tells the reader which. Pressure nodes filtered out as
in MPCO; pressure surfaces via the `PRESSURE` nodal result.

---

## 5. `MODEL/LOCAL_AXES` + per-step frame (D-schema-3)

> **The #1 apeGmsh gap.** MPCO writes **no** beam local axes (`MODEL/LOCAL_AXES` is
> empty for disp/force/elastic beam-columns even with `-E … localAxes`), so apeGmsh
> cannot orient line/section force diagrams from a `.mpco` file at all — it falls back
> to an **identity quaternion that masks "rotation lost in export" as "no rotation"**,
> or forces the native `model.h5` path. Ladruno **must** write real per-element
> beam frames. The recorder derives the frame from the element's `CrdTransf`/`vecxz`
> when the `localAxes` response is absent — never silently defaults to identity.

- **Reference frame (static):** one group per element class (mirrors `ELEMENTS`
  grouping, matching apeGmsh's reader): `MODEL/LOCAL_AXES/<classTag>-<ClassName>[<g>]/{
  FRAME[nE×4], ID[nE] }`, unit quaternion `(w,x,y,z)` taking global → element-local at
  the reference config. **Covers beams** (the gap above), shells, and any element with
  a local frame. A missing frame is an error, not an identity default.
- **Per-step frame (opt-in, NEW):** for co-rotational / large-rotation elements
  (Belytschko beams) that report a current frame, the recorder writes a time series:

```
RESULTS/ON_ELEMENTS/LOCAL_AXES/
   ID            [nE×1] int
   DATA/STEP_<k> [nE×4] double     current-config quaternion
```

Enabled by the `localAxes` element request when the element reports a step-varying
frame. Lets apeGmsh orient deformed-state section resultants (N, V, M) correctly
without recomputing the frame from node motion. Elements with a fixed frame write only
the static dataset.

---

## 6. `MODEL/SECTION_ASSIGNMENTS` — fiber **and** resultant sections

```
SECTION_<secTag>[<SectionClass>]
   attrs:  ID(int), NAME(string), KIND(string "fiber"|"resultant")
   ASSIGNMENT       [nPairs × 2] int     (elemTag, gpIdx0based)
   # fiber sections only:
   FIBER_DATA       [nFibers × 3] double (y, z, area|thickness)   section-local axes
   FIBER_MATERIALS  [nFibers × 1] int
   # resultant sections (Belytschko): no fiber arrays — RESULTANTS names the components
   RESULTANTS       (string)             e.g. "N,Vy,Vz,T,My,Mz"
```

`KIND="resultant"` makes the non-fiber case explicit (Belytschko beams, elastic
sections): the element↔gp↔section mapping is still recorded, but there are no fibers —
the section response is force/deformation resultants (`DATA_TYPE=BeamForceDeformation`).
Shell layered sections keep the 1-based→0-based fiber shift handling from MPCO.

---

## 7. `RESULTS`

### 7.0 Streaming time-series layout (D3 — replaces per-step `DATA/STEP_<k>`)

**Implemented** (ADR [[07_adr_post_review_storage]] D3). Every streaming result group
(node, element, domain, region) stores its history as **one chunked, extensible
dataset** rather than a `DATA` *group* of per-step datasets:

```
<result-or-bucket>/
   ID    [nIds × 1] int
   DATA  [T × nIds × nComp] double   ← extensible on T, chunked (layout: §7b.2)
   STEP  [T] int                     ← commit/step id per slab (was the STEP attr)
   TIME  [T] double                  ← pseudo-time per slab (was the TIME attr)
```

One `[1 × nIds × nComp]` slab is appended per committed step (`H5Dset_extent` +
hyperslab). **Why:** the old per-step layout exploded to 10⁵–10⁶ tiny datasets on long
transients, defeated time-axis compression, and made a single time-history an O(T)
dataset-open; the chunked form gives O(1) slice reads and compresses the (highly
compressible) time axis. **Values are identical** to the per-step form — only the on-disk
organization changed (parity harness verifies 1e-12 against the frozen per-step `.mpco`).
The reader treats a `DATA` *dataset* (chunked) and a `DATA` *group* (legacy/`.mpco`)
interchangeably (`ladruno_format.iter_step_slices`). The notation `DATA/STEP_<k>` below is
the **superseded** per-step form, kept for cross-referencing the frozen recorder.

### 7.1 `ON_NODES/<RESULT>`

As MPCO, retained: group attrs `DISPLAY_NAME, COMPONENTS, DIMENSION, DESCRIPTION,
TYPE, DATA_TYPE`; `ID[nN×1]`; results stored per the **§7.0 chunked layout** (`DATA
[T×nN×nComp]` + `STEP`/`TIME` axes). Modes of vibration keep the `STEP_<k>` *group* →
`MODE_<i>` with `LAMBDA/OMEGA/FREQUENCY/PERIOD`.

**`PARTITION_REDUCTION`** (string attr, added by WP-126, written on every result group —
`ON_NODES`, `ON_ELEMENTS`, `ON_DOMAIN`, `ON_REGIONS`, and the `ENVELOPES` groups).
It tells a reader how to combine one id's rows when the id appears in several
`.part-N` files (a node on a partition interface):

| Value | Meaning | Reader obligation when stitching |
|---|---|---|
| `NONE` | consistent: every partition holds the same value (kinematics, element results) | keep any one copy |
| `SUM` | additive: each partition holds the partial from its own elements (`reactionForce`, `reactionMoment`, `…IncludingInertia`, the Rayleigh variants `RAYLEIGH_FORCE`/`RAYLEIGH_MOMENT`, `unbalancedForce`, `unbalancedMoment`) | **sum** the copies, component-wise |
| `UNSUPPORTED` | partial per partition and not recoverable by a sum (`energyBalance`: shared-node KE counted in each partition; `RES`/`ERR` derived) | do not present a merged value; warn |

Optional for readers: files written before WP-126 lack it. A reader should then fall back
to the result name (the `REACTION*` / `RAYLEIGH*` / `UNBALANCED*` groups are `SUM`). In a
**partitioned** file an `ENVELOPES` group is always `NONE`: the recorder refuses to
envelope a `SUM`/`UNSUPPORTED` source there, with a warning, because a per-partition
extreme of a partial cannot be recombined (max(a+b) ≠ max a + max b). Serial files may
carry a `SUM` envelope (there is only one partition). Validator:
`ladruno_format._validate_partition_reduction`.

#### 7.1.1 Kinematic node results: names and `COMPONENTS` labels (apeGmsh #2128)

Translation and rotation are **two separate result groups**. There is no combined
6-component displacement group. A reader that draws a deformed shape reads
`DISPLACEMENT`, and adds `ROTATION` only if it needs nodal rotations. Mode shapes follow
the same split. Values below are copied from the writer source. `COMPONENTS` is one
comma-separated string attr (`h5::group::createResultGroup`, `SRC/recorder/Ladruno_Hdf5.h`).

| Recorder token (`-N`) | Group `ON_NODES/<name>` | `DISPLAY_NAME` | ndm = 1 | ndm = 2 | ndm = 3 | `DIMENSION` | Writer |
|---|---|---|---|---|---|---|---|
| `displacement` | `DISPLACEMENT` | `Displacement` | `Ux` | `Ux,Uy` | `Ux,Uy,Uz` | `L` | `DisplacementSource::build()` |
| `rotation` | `ROTATION` | `Rotation` | `Rx,Ry,Rz` | `Rz` | `Rx,Ry,Rz` | `A` | `RotationSource::build()` |
| `modesOfVibration` | `MODES_OF_VIBRATION(U)` | `Modes of Vibration (U)` | `Ux` | `Ux,Uy` | `Ux,Uy,Uz` | `""` | `ModesOfVibrationSource::build()` |
| `modesOfVibrationRotational` | `MODES_OF_VIBRATION(R)` | `Modes of Vibration (R)` | `Rx,Ry,Rz` | `Rz` | `Rx,Ry,Rz` | `""` | `ModesOfVibrationRotationalSource::build()` |

All four writers are in `SRC/recorder/Ladruno_NodeResults.cpp`. The token-to-type mapping
is in the `-N` option parser in `SRC/recorder/LadrunoRecorder.cpp`.

**The labels depend on ndm only.** They do not depend on ndf, the element family, or
the analysis. `ndm` is `INFO/SPATIAL_DIM`: the coordinate count of the **first node** in
the domain (`LadrunoRecorder.cpp`, model-write setup, rejected outside 1..3). Each
`build()` fixes the labels once, at construction. A mixed-ndf model therefore has one
label set for every node.

**ndf affects only the values** (`bufferNodeResponseVec3u` / `bufferNodeResponseVec3r`,
`SRC/recorder/Ladruno_Types.h`):

- `DISPLACEMENT` holds the node's first `ndm` DOFs of `getTrialDisp()`. Missing DOFs are
  written as `0.0`.
- `ROTATION` with ndm = 2 holds DOF index 2, or `0.0` if the node has ndf ≤ 2.
  With ndm = 3 it holds DOF indices 3, 4, 5 when ndf > 5, and all zeros otherwise.
  With ndm = 1 it writes three columns of zeros.
- The modal groups use the same rules on the eigenvector rows. A node with no
  eigenvector for the current mode stays `0.0`.

**Reader caveat (not a label ambiguity).** `ROTATION` (and `MODES_OF_VIBRATION(R)`)
copies a DOF index. It does not check what the DOF means. With ndm = 2 and ndf = 3 on a
u-p node (the `*UP` / pressure-DOF elements), `Rz` holds the **pore pressure**, not a
rotation. A reader that needs a true rotation must check the element family itself.
`DISPLACEMENT` is not affected.

#### 7.1.2 Every other `ON_NODES` group: names and `COMPONENTS` labels (apeGmsh #2214)

This completes §7.1.1. The tables list every group the writer can put in
`ON_NODES`. Unless a row says otherwise, the writer is in
`SRC/recorder/Ladruno_NodeResults.cpp`. The `-N` / `-NS` tokens are parsed in
`SRC/recorder/LadrunoRecorder.cpp`. §7.1.1 rules apply here too: labels depend on ndm
(`INFO/SPATIAL_DIM`) only, are fixed at construction, and missing DOFs are written as `0.0`.

Each group belongs to one of two families, by the buffer helper it uses
(`SRC/recorder/Ladruno_Types.h`):

- **Translational** (`bufferNodeResponseVec3u`): one label per ndm, from the node's
  first `ndm` DOFs.
- **Rotational** (`bufferNodeResponseVec3r`): ndm = 2 gives one `z` label from DOF index
  2; ndm = 3 gives `x,y,z` from DOF indices 3, 4, 5 (when ndf > 5); ndm = 1 gives
  `x,y,z` written as zeros. The §7.1.1 u-p caveat applies to every rotational group:
  with ndm = 2 and ndf = 3, the `z` column holds whatever DOF 2 is.

**Translational groups**

| `-N` token | Group `ON_NODES/<name>` | `DISPLAY_NAME` | ndm = 1 | ndm = 2 | ndm = 3 | `DIMENSION` | Writer | Node vector read |
|---|---|---|---|---|---|---|---|---|
| `velocity` | `VELOCITY` | `Velocity` | `Vx` | `Vx,Vy` | `Vx,Vy,Vz` | `L/t` | `VelocitySource::build()` | `getTrialVel()` |
| `acceleration` | `ACCELERATION` | `Acceleration` | `Ax` | `Ax,Ay` | `Ax,Ay,Az` | `L/t^2` | `AccelerationSource::build()` | `getTrialAccel()` |
| `reactionForce` | `REACTION_FORCE` | `Reaction Force` | `RFx` | `RFx,RFy` | `RFx,RFy,RFz` | `F` | `ReactionForceSource::build()` | `getReaction()` |
| `reactionForceIncludingInertia` | `REACTION_FORCE_INCLUDING_INERTIA` | `Reaction Force Including Inertia` | `RFx` | `RFx,RFy` | `RFx,RFy,RFz` | `F` | `ReactionForceSource::build()`, called by `ReactionForceIncInertiaSource` | `getReaction()` |
| `rayleighForce` | `RAYLEIGH_FORCE` | `Rayleigh Force` | `RFx` | `RFx,RFy` | `RFx,RFy,RFz` | `F` | `ReactionForceSource::build()`, called by `RayleighForceSource` | `getReaction()` |
| `unbalancedForce` | `UNBALANCED_FORCE` | `Unbalanced Force` | `Fx` | `Fx,Fy` | `Fx,Fy,Fz` | `F` | `UnbalancedForceSource::build()` | `getUnbalancedLoad()` |
| `unbalancedForceIncludingInertia` | `UNBALANCED_FORCE_INCLUDING_INERTIA` | `Unbalanced Force Including Inertia` | `Fx` | `Fx,Fy` | `Fx,Fy,Fz` | `F` | `UnbalancedForceSource::build()`, called by `UnbalancedForceIncInertiaSource` | `getUnbalancedLoadIncInertia()` |
| `constraintTieForce` (alias `tieForce`) | `CONSTRAINT_TIE_FORCE` | `Constraint Tie Force` | `TFx` | `TFx,TFy` | `TFx,TFy,TFz` | `F` | `ConstraintTieForceSource::build()` | `getProjectionTieForce()` |

**Rotational groups**

| `-N` token | Group `ON_NODES/<name>` | `DISPLAY_NAME` | ndm = 1 | ndm = 2 | ndm = 3 | `DIMENSION` | Writer | Node vector read |
|---|---|---|---|---|---|---|---|---|
| `angularVelocity` | `ANGULAR_VELOCITY` | `Angular Velocity` | `RVx,RVy,RVz` | `RVz` | `RVx,RVy,RVz` | `L/t*A` | `AngularVelocitySource::build()` | `getTrialVel()` |
| `angularAcceleration` | `ANGULAR_ACCELERATION` | `Angular Acceleration` | `RAx,RAy,RAz` | `RAz` | `RAx,RAy,RAz` | `L/t^2*A` | `AngularAccelerationSource::build()` | `getTrialAccel()` |
| `reactionMoment` | `REACTION_MOMENT` | `Reaction Moment` | `RMx,RMy,RMz` | `RMz` | `RMx,RMy,RMz` | `F*L` | `ReactionMomentSource::build()` | `getReaction()` |
| `reactionMomentIncludingInertia` | `REACTION_MOMENT_INCLUDING_INERTIA` | `Reaction Moment Including Inertia` | `RMx,RMy,RMz` | `RMz` | `RMx,RMy,RMz` | `F*L` | `ReactionMomentSource::build()`, called by `ReactionMomentIncInertiaSource` | `getReaction()` |
| `rayleighMoment` | `RAYLEIGH_MOMENT` | `Rayleigh Moment` | `RMx,RMy,RMz` | `RMz` | `RMx,RMy,RMz` | `F*L` | `ReactionMomentSource::build()`, called by `RayleighMomentSource` | `getReaction()` |
| `unbalancedMoment` | `UNBALANCED_MOMENT` | `Unbalanced Moment` | `Mx,My,Mz` | `Mz` | `Mx,My,Mz` | `F*L` | `UnbalancedMomentSource::build()` | `getUnbalancedLoad()` |
| `unbalancedMomentIncludingInertia` | `UNBALANCED_MOMENT_INCLUDING_INERTIA` | `Unbalanced Moment Including Inertia` | `Mx,My,Mz` | `Mz` | `Mx,My,Mz` | `F*L` | `UnbalancedMomentSource::build()`, called by `UnbalancedMomentIncInertiaSource` | `getUnbalancedLoadIncInertia()` |

The plain, `…IncludingInertia` and Rayleigh variants share labels. Only the group name
tells them apart. Each reaction source passes its own `reactionFlag()` to the reaction
computation, so the three read different reactions through the same `getReaction()`.

**Scalar groups (one component, label independent of ndm)**

| Token | Group `ON_NODES/<name>` | `DISPLAY_NAME` | `COMPONENTS` | `DIMENSION` | Writer |
|---|---|---|---|---|---|
| `-N pressure` | `PRESSURE` | `Pressure` | `p` | `F/L^2` | `PressureSource::build()` |
| none (one group per ready `LadrunoPorousOverlay`) | `overlayPressure_<patternTag>` | `Overlay pore pressure <patternTag>` | `p` | `F/L^2` | `OverlayPressureSource::buildSchema()`, `SRC/recorder/Ladruno_OverlayResults.cpp` |

The `PRESSURE` **label** is fixed. The **value** depends on the element family attached
to the node (`PressureSource::evaluate()` / `scanUPContract()`): it reads DOF `ndm` (index
2 in 2D, 3 in 3D) only when the node has exactly `ndm + 1` DOFs. A node attached to a
`LadrunoUP` element is read from the displacement vector. A node attached only to an
upstream `*UP` element (FourNodeQuadUP, BrickUP, Nine_Four_Node_QuadUP,
Twenty_Eight_Node_BrickUP, BBarFourNodeQuadUP, BBarBrickUP, SSPquadUP, SSPbrickUP) is read
from the velocity vector, because their DOF is ∫p dt. A node attached to both reads the
displacement vector and the recorder warns once. Every other node reads the velocity
vector, and is `0.0` unless it has exactly `ndm + 1` DOFs. The overlay group name is
lowercase-camel and has no `-N` token. The recorder adds it from `initOverlaySources()`
in `LadrunoRecorder.cpp`.

**Sensitivity groups (`-NS <token> <gradIndex>`).** The gradient index `<g>` is part of
the group name and of every label. The families and ndm rules are as above:

| `-NS` token | Group `ON_NODES/<name>` | ndm = 1 | ndm = 2 | ndm = 3 | Writer |
|---|---|---|---|---|---|
| `displacementSensitivity` | `DISPLACEMENT_SENSITIVITY_<g>` | `d_Ux/d_q<g>` | `d_Ux/d_q<g>,d_Uy/d_q<g>` | `d_Ux/d_q<g>,d_Uy/d_q<g>,d_Uz/d_q<g>` | `DisplacementSensitivitySource::build()` |
| `rotationSensitivity` | `ROTATION_SENSITIVITY_<g>` | `d_Rx/d_q<g>,d_Ry/d_q<g>,d_Rz/d_q<g>` | `d_Rz/d_q<g>` | `d_Rx/d_q<g>,d_Ry/d_q<g>,d_Rz/d_q<g>` | `RotationSensitivitySource::build()` |
| `velocitySensitivity` | `VELOCITY_SENSITIVITY_<g>` | `d_Vx/d_q<g>` | `d_Vx/d_q<g>,d_Vy/d_q<g>` | `d_Vx/d_q<g>,d_Vy/d_q<g>,d_Vz/d_q<g>` | `VelocitySensitivitySource::build()` |
| `angularVelocitySensitivity` | `ANGULAR_VELOCITY_SENSITIVITY_<g>` | `d_RVx/d_q<g>,d_RVy/d_q<g>,d_RVz/d_q<g>` | `d_RVz/d_q<g>` | `d_RVx/d_q<g>,d_RVy/d_q<g>,d_RVz/d_q<g>` | `AngularVelocitySensitivitySource::build()` |
| `accelerationSensitivity` | `ACCELERATION_SENSITIVITY_<g>` | `d_Ax/d_q<g>` | `d_Ax/d_q<g>,d_Ay/d_q<g>` | `d_Ax/d_q<g>,d_Ay/d_q<g>,d_Az/d_q<g>` | `AccelerationSensitivitySource::build()` |
| `angularAccelerationSensitivity` | `ANGULAR_ACCELERATION_SENSITIVITY_<g>` | `d_RAx/d_q<g>,d_RAy/d_q<g>,d_RAz/d_q<g>` | `d_RAz/d_q<g>` | `d_RAx/d_q<g>,d_RAy/d_q<g>,d_RAz/d_q<g>` | `AngularAccelerationSensitivitySource::build()` |

Their `DISPLAY_NAME` is `<Quantity> Sensitivity d_<sym>/d_q<g>`, and `DIMENSION` is `""`.
`<g>` is the 1-based sensitivity parameter index. If `<g>` is out of range
(`< 1` or `> getNumParameters()`), the group is written as all zeros and the recorder warns.

### 7.2 `ON_ELEMENTS/<result>` — structured column map (replaces the `;`-string)

```
ON_ELEMENTS/<result>/<classTag>-<ClassName>[<g>:<h>]
   attr NUM_COLUMNS (int)
   COLUMN_MAP/                          ← structured, replaces flattened META
     LEVELS        [k]  int   block path depth code (0=Elem 1=Gauss 2=Section 3=Fiber 4=Material)
     GAUSS_ID      [k]  int   0-based gp, -1 if N/A
     SECTION_TAG   [k]  int   -1 if N/A
     FIBER_ID      [k]  int   -1 if N/A
     NUM_COMP      [k]  int
     MULTIPLICITY  [k]  int
     COMP_NAMES    (string ATTR on COLUMN_MAP)  newline-separated, k lines, each a CSV
                                                of component names (one line per block)
   SECTION_MAP   [nElem × NUM_GP] int   section tag per (elem, gp); -1 if none   ← NEW
   ID            [nElem × 1] int
   DATA          [T × nElem × NUM_COLUMNS] double + STEP[T] + TIME[T]   ← §7.0 chunked
```

`<g>` = element/basis group; `<h>` = response-shape (header) bucket — two same-class
elements with different fiber counts land in different `<h>`. Invariant:
`NUM_COLUMNS == Σ MULTIPLICITY[i]·NUM_COMP[i]`. The structured `COLUMN_MAP` removes the
bespoke `;`/int-code decoder STKO required.

- **`SECTION_MAP` (perf, NEW):** a pre-built `[nElem × NUM_GP]` section-tag index so the
  reader resolves `(elem, gp) → SECTION_<tag>` in O(1). Today apeGmsh re-walks
  `SECTION_ASSIGNMENTS` `O(n_sections × n_elems)` on **every** read; writing this once
  on the recorder side deletes that walk.
- **Homogeneity is enforced at write.** Every element in a `[<g>:<h>]` bucket has an
  identical column schema by construction. The recorder **fails loud** if an element
  would break bucket homogeneity, rather than emitting a ragged bucket the reader
  silently skips (apeGmsh's current failure mode).
- **Beam internal forces** are recorded as **section/station resultants at the
  `GP_PARAM` positions** with a single documented sign convention (internal-section
  force, consistent across stations). This replaces MPCO's `localForce` end-force path,
  whose hardcoded station-2 sign flip and synthetic `[-1,+1]` stations only work for the
  2-station case. apeGmsh's `line_force` diagrams read these directly.

### 7.3 `ON_DOMAIN/<scalar>` — global results **(DEFERRED)**

Single `[1 × nComp]` row per step, stored per the §7.0 chunked layout (`DATA[T×1×nComp]`).
Energy balance, base shear, total
reaction. **Built only after the separate energy work lands** — Ladruno will
*consume* `EnergyBalanceRecorder` ([[project_energy_balance_feature]]), not reinvent
it. Spec stub here for forward compatibility.

### 7.4 `ENVELOPES/ON_{NODES,ELEMENTS,DOMAIN,REGIONS}/<name>` — (v3a, IMPLEMENTED)

Per the ADR (D6/D7): componentwise `MIN`/`MAX`/`ABSMAX` + `ARG_STEP`, per-`MODEL_STAGE`,
periodic in-place rewrite. **Opt in with the `-envelope` recorder flag** — every
non-modal node/element/domain channel then uses an `EnvelopeSink` (running extremes)
instead of the time-series `StreamingSink`; `ModesOfVibration` (eigen) is never
enveloped. The accumulators are written only on `finalize()`, so the recorder
rewrites them in place every recorded step (crash safety) and once at teardown.
Each `ENVELOPES/<family>/<name>` group carries the **same self-describing attrs**
as the time-series result group (`COMPONENTS`/`DISPLAY_NAME`/`DIMENSION`/…, §7.1) so
the extremes are named in-file. v3a covers consistent quantities (zero-MPI); v3b
reduced (reaction/global) envelopes need the per-step Allreduce, which is compiled
out of this shared sequential recorder TU (see `LEDGER_quirks.md`) — deferred.
**Element envelopes** nest one level (`ENVELOPES/ON_ELEMENTS/<display>/<bucket>`, like
the time-series tree) and carry the **same structured per-column `COLUMN_MAP`** as the
time-series `ON_ELEMENTS` path (§7.2): `LEVELS`/`GAUSS_ID`/`SECTION_TAG`/`FIBER_ID`/
`NUM_COMP`/`MULTIPLICITY` + `COMP_NAMES` + the `NUM_COLUMNS` attr. This is the
authoritative gauss/section/fiber column structure — element results have an *empty*
result-level `COMPONENTS` attr, so the `COLUMN_MAP` is how an envelope column is named
and located. The `EnvelopeSink` delete-recreates its leaf group on every in-place flush,
so the recorder re-stamps the `COLUMN_MAP` after each `finalize()`.

```
ENVELOPES/ON_NODES/<name>/                          (node/domain: flat)
   ID      [nN×1] int
   MIN     [nN×nComp]   MAX [nN×nComp]   ABSMAX [nN×nComp]   double
   ARG_STEP[nN×nComp]   int     step index of each extreme (commitTag, global across ranks)

ENVELOPES/ON_ELEMENTS/<display>/<bucket>/           (element: nested + COLUMN_MAP)
   ID      [nE×1] int
   MIN/MAX/ABSMAX [nE×NUM_COLUMNS] double   ARG_STEP [nE×NUM_COLUMNS] int
   COLUMN_MAP/  { LEVELS, GAUSS_ID, SECTION_TAG, FIBER_ID, NUM_COMP, MULTIPLICITY,
                  COMP_NAMES (ATTR) }                ← §7.2, identical structure
```

---

## 7.5 `ON_PARTICLES` — MPM particle results (`FORMAT_MINOR_VERSION = 2`; WP-232a)

`MODEL_STAGE[*]/MODEL/PARTICLES` and `RESULTS/ON_PARTICLES` add MPM particle output.
The change is additive. `FORMAT_VERSION` stays 1, and a file that holds a particle group writes
`FORMAT_MINOR_VERSION = 2`.

**The full contract is
[[232a_mpm_contracts]] §2. It is the only copy until the move into this document, which is
deferred.** It covers the field-by-field layout, the rule codes, the crash rule for a torn
`OFFSETS` entry, and the fixtures under `tests/data/mpm/`. Its executable form is
`ladruno_format.validate()` (`_validate_particles`).

---

## 7b. File robustness & provenance (OpenSees-side)

The apeGmsh consultation surfaced two robustness failures the recorder owns at the
source (the user's stated priority: performance + robustness on OpenSees first):

- **Deterministic close + periodic flush.** Writing an MPCO file onto a synced FS
  (OneDrive/SeaDrive/Dropbox) native-crashes the kernel at `wipe()` — the deferred
  `H5Fclose` races the sync client and leaves the HDF5 write-flag stuck;
  `HDF5_USE_FILE_LOCKING=FALSE` does not fix it. Combined with the `exit(-1)`-leaves-
  files-unflushed history ([[project_mpco_exit_crash]]), Ladruno must: flush every
  recorded step (as MPCO does), **rewrite envelope/aggregate datasets periodically**
  (§7.4), and provide an **explicit flush/close recorder control** so a user can force a
  clean close in the same cell before copying off local disk. (WP-164 replaced the
  per-step flush with the `-flush` cadence. What a killed run actually keeps, and how to
  read it, is §7b.2.)
- **Provenance / lineage slot.** The recorder creates `INFO/PROVENANCE` as an
  **empty group** at open (`LadrunoRecorder.cpp`, beside the `RUN_ID` stamps) and
  writes nothing into it. It is the one place in a `.ladruno` file that a consumer may
  write into. The rules are in §7b.1.

### 7b.1 `INFO/PROVENANCE`: consumer-written pairing attributes

*Ratified in apeGmsh #2079 (2026-10-15, open points 2 and 3); the apeGmsh side is
apeGmsh ADR 0122.* `.ladruno` is apeGmsh's canonical results format. After a run,
apeGmsh stamps the attributes below so its viewer can pair the results file with the
model that produced it (`<stem>.h5`), without a sidecar file.

**Ownership.** The fork owns the `.ladruno` format: the layout, the buckets, the token
names, `FORMAT_VERSION`/`FORMAT_MINOR_VERSION` and `ladruno_format.py`. A consumer
(apeGmsh) owns only the `INFO/PROVENANCE` attributes it stamps, as listed here.

**Write rules.** A consumer that stamps provenance:
1. writes **attributes only**: no datasets, no subgroups;
2. writes **only in `INFO/PROVENANCE`**: never in `INFO` itself, `MODEL_STAGE[*]`, or
   any other group, and never modifies an attribute the recorder wrote;
3. writes **only on a closed file**: after the recorder has closed it (`remove
   recorder`, `wipe`, or interpreter exit) and with no OpenSees process holding it
   open, and with **no other process** (OpenSees, a viewer, an `h5py` session) holding
   it open. If the file will not open for writing, the consumer skips stamping and
   warns; it never disables HDF5 file locking to force the write. A file that was never
   closed cleanly (a crash or a killed process) is not stamped: the consumer SHOULD
   require `ladruno_format.validate()` to pass first. For a partitioned run, the consumer
   stamps each `<stem>.part-<N>.ladruno` it pairs;
4. writes **lowercase** attribute names. The case marks the owner: lowercase in
   `PROVENANCE` belongs to the consumer, and UPPERCASE stays reserved for the fork (an
   optional, later recorder-side `-provenance` flag would write UPPERCASE names here).
   The `MODEL_HASH`/`RESULTS_HASH` slots this section used to reserve were never written
   by the recorder; they are withdrawn as consumer slots and remain fork-reserved names.
   The lowercase names are a deliberate exception to this document's UPPERCASE naming
   convention: they are the only names a consumer writes.

| Attribute | Type | Meaning |
|---|---|---|
| `model_hash` | string | apeGmsh model lineage hash, copied byte-equal from `<stem>.h5` `/meta/lineage` |
| `fem_hash` | string | apeGmsh FEM-snapshot lineage hash, copied byte-equal from `<stem>.h5` `/meta/lineage` |
| `stage_names` | string (JSON array of strings) | program stage names, one per `MODEL_STAGE[*]` in stage order; written only when the stage count and kinds match the program's stages |
| `element_index_space` | string | `"fem"` when every element ID in the file is a paired FEM element id; absent otherwise (the reader treats absence as "unverified") |
| `unsupported` | string (JSON array of strings) | categories present in the source that the consumer did not read or convert (e.g. layered shells, energy, MPCO springs); `"[]"` when none |
| `writer` | string | the apeGmsh component and version that stamped (or, for a converted file, produced) the file, e.g. `apeGmsh <version>` or `apeGmsh mpco->ladruno <version>`. A converted file keeps `INFO/GENERATOR="Ladruno"`; `writer` is what marks it as converted |
| `opensees_backend` | string | *optional; added 2026-10-19 (apeGmsh G10).* Which OpenSees ran the analysis: exactly `"ladruno"` (this fork), `"stock"` (upstream OpenSees), or `"unknown"`. The `.ladruno` recorder exists only in the fork, so a file the recorder wrote is always `"ladruno"`. `"stock"` can occur only on a converted file (`writer` = `apeGmsh mpco->ladruno …`) whose source `.mpco` came from upstream OpenSees |
| `opensees_build_sha` | string | *optional; added 2026-10-19 (apeGmsh G10).* The git commit the OpenSees binary was built from: 40 lowercase hex characters, exactly as the fork's `ladrunoBuild` command returns it (Tcl `ladrunoBuild`, Python `ops.ladrunoBuild()`; CMake stamps `git log -1 --format=%H` into `OPENSEES_VERSION`). `"unknown"` when the build reported no hash (`ladrunoBuild` itself returns `"unknown"` then) or the backend could not be queried. The consumer writes it as returned and does not shorten it |

All of these are scalar attributes stored as variable-length UTF-8 strings (the apeGmsh viewer's stamp check expects vlen UTF-8, not fixed-length). JSON values are compact, UTF-8,
and never NaN. A reader must tolerate any subset of these attributes being absent (a
file that was never stamped has an empty `PROVENANCE`). `opensees_backend` and
`opensees_build_sha` are newer than the other six, so files stamped before 2026-10-19
lack them and stay valid. A reader treats their absence as unknown.

**No fork-side equivalent today.** The recorder writes no build hash into a `.ladruno`
file. `INFO/GENERATOR="Ladruno"` marks the writer, and `INFO/SOLVER_VERSION` is a fixed
`[3,5,1]`, not a commit. The build hash is available only at run time, through
`ladrunoBuild`. That is why the consumer stamps it here. If the recorder later writes its
own build hash (an UPPERCASE name, e.g. via the `-provenance` flag), the fork's attribute
takes precedence and `opensees_build_sha` becomes redundant.

**Partitioned sets.** The attributes are per file. `stage_names` describes that file's
own `MODEL_STAGE[*]` (part files may differ, WP-165). A reader treats a partitioned set as
paired only when every part carries the same `model_hash` and `fem_hash`; a set that was
stamped only partly (a consumer that stopped part-way) is unpaired.

**Versioning.** This is additive. `FORMAT_VERSION` stays `1` and
`FORMAT_MINOR_VERSION` does not move: the recorder writes nothing new, and a v1 reader
that ignores `PROVENANCE` still reads every file. `ladruno_format.validate()` does not
inspect `PROVENANCE`, so a stamped file stays conformant.

### 7b.2 A killed run: the writer's guarantees and the reader rule (WP-233, #1082)

*Applies to every `.ladruno`: `recorder ladruno` since WP-233 (#1082), and the MPM particle
writer since WP-232c (#1081), whose particle-specific crash rule stays in
[[232a_mpm_contracts]] §2.5.*

**What a kill used to do.** Before WP-233 the recorder wrote the earliest HDF5 format (v1
B-tree chunk indexes), without SWMR. It deflated a partially filled `DATA` chunk at every
flush, so the chunk moved and its old space was reused, while the chunk index still pointed
there. A kill (TerminateProcess, a scheduler kill, a hard Ctrl-C) could therefore leave
`DATA` that would not read ("filter returned failure during read"), losing the whole
unfinished chunk (up to 1,024 steps, often the whole run). It could also leave `DATA` that
read silently wrong: zeros, or another result's bytes. That happened at 72–79 % of the kill
points inside a flush. The old promise here, "a hard crash loses at most the last `-flush`
seconds", did not hold (quirk `WP-233-recorder-killed-without-swmr`).

**The writer's guarantees.**

- **SWMR except during a layout step.** The file uses the latest format
  (`H5Pset_libver_bounds(LATEST, LATEST)`; extensible-array chunk indexes). It is in
  SWMR-write mode (`H5Fstart_swmr_write`) except during a **layout step**: a recorded step
  that creates objects. SWMR forbids creating objects or attributes. The layout steps are:
  - a stage's first recorded step (`writeModel`, every result's group with its `ID` and
    empty `DATA`/`TIME`/`STEP`, the element `COLUMN_MAP`s, the envelope groups, and the
    groups of a modes result);
  - the rewrite of a provisional `KIND` (§2a), at most once per stage;
  - a modes snapshot (`DATA/STEP_<k>/MODE_<j>`, each time the eigenvalues change);
  - an envelope whose first sample arrives after the stage's first step.
- **How a layout step runs.** If SWMR is armed, the recorder flushes, closes the file (a
  clean close clears the SWMR flags) and reopens it read-write. It then creates the objects
  and holds the step's rows. Next it flushes and arms SWMR with no object open. Finally it
  appends the held rows and flushes. No streaming row is ever written outside SWMR. If SWMR
  cannot be armed, the recorder disables itself with an error rather than writing on
  without it.
- **What a layout step may touch.** It writes new objects, the superblock, and one
  in-place header write:
  - the root's header gains the new stage's `MODEL_STAGE` link;
  - a modes result's `DATA` gains a `STEP_<k>` link;
  - the stage's header gets the new `KIND`.

  That is one write because the root and each modes `DATA` keep their links compact in an
  object header presized for 1,024 stages (1,024 snapshots; about 44 KB each, and HDF5 caps
  one header chunk near 1,500 such links), and because the provisional `KIND` is a
  fixed-length 9-byte string rewritten in place. Without that, the replay showed a kill
  inside the layout flush leaving the root unreadable, with every stage lost.
  **Past 1,024 stages** (or 1,024 snapshots in one modes result), adding a link needs a new
  header chunk, and a kill in a later layout step can again make earlier stages (snapshots)
  unreadable. This was measured with the root presized for 8 and 21 stages: 416 of 2,607
  kill points left `MODEL_STAGE[1]` unreadable. The recorder warns when a file reaches
  1,024 stages.
- **Streaming layout.** Under SWMR no freed file space is reused. A `DATA` chunk is
  therefore written so that it never moves:
  - with `-compress N > 0` (the default, 1) and a slab (one step of the result) of at
    least 4 KiB, it holds one step (ids tiled to 1 MiB) and is compressed and written once;
  - otherwise (`-compress 0`, or a slab under 4 KiB, where a per-step chunk compresses worse
    than raw and costs time), it is unfiltered (about 64 KiB, stacking up to 1,024 steps),
    keeps its size and is rewritten in place.

  `TIME`/`STEP` stay unfiltered, chunk 1,024.
- **What survives a kill:**
  - every step up to the last completed flush (`-flush` cadence; the first recorded step of
    a stage always flushes), exact;
  - every stage that ended, in full;
  - the modes snapshots and envelopes as below.

  What is in flight may be partly present. `DATA`, `TIME` and `STEP` are separate datasets,
  and SWMR orders each dataset's own writes, not one dataset against another. So a killed
  file can hold them at different lengths ("ragged"); in the replay this happened at about
  40 % of the in-flush kill points.

**The reader rule (opening any `.ladruno`).**

1. Open read-only. A clean close leaves the superblock's consistency flags at 0x00, and a
   plain open works.
2. If that fails with "already open for write" (or "truncated file", which is checked
   first), open in SWMR-read mode (h5py: `h5py.File(path, "r", swmr=True)`). This is a
   kill while SWMR was armed (flags 0x05), which is almost every kill. Never `h5clear` first:
   the stored end-of-allocation lags the data, and only an SWMR-read open skips that check.
   No `HDF5_USE_FILE_LOCKING=FALSE` is needed after a kill (measured, Windows).
3. If that also refuses ("not already open for SWMR writing") and the flags are 0x01, the
   run was killed during a layout step. Recover a **copy**, never the original:
   - set the flags to 0x00;
   - set the stored end-of-allocation (superblock offset 28) to max(stored EOA, file size);
   - pad the copy with zeros up to it if the file is shorter;
   - re-checksum the superblock (lookup3 over bytes 0–43, stored at 44);
   - open the copy plainly.

After any of the three steps, **the last stage may fail to open or read**: the kill fell
while it was being laid out, and it holds no committed step, so skip it. An earlier stage
that will not open is corruption. This can follow a plain open too. sec2 writes
WRITE_ACCESS into the superblock only when it next flushes it, so a file killed inside a
layout step can still carry flags 0x00. The real-binary kill stress found that at 7 of 16
kills with a stage every 20 steps.

A reader that cannot open in SWMR-read mode (e.g. h5wasm) may apply step 3's repair to a
copy of any file that step 1 refuses, flags 0x05 included.

**The crash rule (what to read).**

- A file is written by one recorder, so every streaming result of a stage shares one step
  axis. Read **T' = the minimum, over the stage's streaming results, of min(len `DATA`,
  len `TIME`, len `STEP`)**, and ignore rows past T'. A closed file has all these lengths
  equal.
- A modes snapshot (`DATA/STEP_<k>`) that does not open, or lacks a `MODE_<j>`, is the
  snapshot being written: skip it, along with any later snapshot.
- An envelope's `MIN`/`MAX`/`ABSMAX`/`ARG_STEP` entries are each the extreme over some
  prefix of the stage's recorded steps. An in-place rewrite torn between two steps can mix
  two neighbouring prefixes. Each value is a true extreme of the run so far, but entries
  may come from consecutive steps.
- `KIND` is `static`, `transient` or `eigen`. A provisional `static` that was never
  rewritten means no analysis was visible at any recorded step (§2a).

**Reference implementation:** `ladruno_format.py`:
- `open_ladruno(path)` (steps 1–2; step 3 raises `StageLayoutKilledError`);
- `repair_copy(src, dst)` (step 3);
- `stage_report(f)` (every stage's T', or the torn last stage);
- `recover_ladruno(src, dst)` (step 3, then `stage_report`);
- `committed_steps(stage)` (T');
- `LadrunoReader`, which applies both rules.

`validate()` checks **closed** files and refuses a killed one (`FILE-KILLED`). Read a killed
file by these rules; do not validate it. The particle file's helpers in
`tests/mpm_contracts.py` delegate to the same functions.

**Evidence:**
- `Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py` is the exhaustive replay:
  every byte-level kill point of an h5py mirror of the writer, read by these rules;
- `recorder_binary_kill_stress.py` is the kill stress of the real binary, per lifecycle
  window;
- the numbers are on #1082 and in the quirk.

## 7c. apeGmsh consumption contract (interop summary)

The `.ladruno` reader (apeGmsh `Results.from_ladruno`, `STKO_to_python` parallel
branch) is the canonical consumer. Hard contract points:

| Reader needs | Schema provides | Replaces apeGmsh workaround |
|---|---|---|
| Structured column schema | `COLUMN_MAP` (§7.2) | ~200 lines of `;`-string parse |
| GP coords + weights, per group | `QUADRATURE/{GP_PARAM,GP_WEIGHT}` (§3) | hardcoded `RESPONSE_CATALOG` lookup |
| Beam local frame, per element | `LOCAL_AXES/<class>/FRAME` (§5) | identity-quaternion fallback |
| `(elem,gp)→section` O(1) | `SECTION_MAP` (§7.2) | per-read `SECTION_ASSIGNMENTS` walk |
| Stage kind | `MODEL_STAGE/@KIND` | "assume transient" |
| Format/version gate | `INFO/{GENERATOR,FORMAT_VERSION}` | none → silent multi-partition corruption |
| Contiguous partition discovery | `<stem>.part-<N>.ladruno` (0-based) | — (already matched, now pinned) |
| Model ↔ results pairing | `INFO/PROVENANCE` consumer attributes (§7b.1) | a sidecar file or a native `results.h5` copy |

Self-contained ⇒ no `model_h5=` sibling: regions (`SETS`), local axes, section index
all live in the one file.

## 8. Worked examples

**(a) Legacy 4-node quad — parity, derived descriptor.**
`MODEL/ELEMENTS/156-FourNodeQuad[0]`: `TOPOLOGY=quad FAMILY=lagrange ORDER=[1,1]
PARAM_DOMAIN=[-1,1] RATIONAL=0 NUM_CTRL=4 NUM_GP=4`; `GP_PARAM` = 2×2 tensor of ±1/√3;
`GP_WEIGHT` = [1,1,1,1]. No element edit — derived from the legacy table.

**(b) Quadratic rational Bézier quad — self-declared.**
`…-BezierQuad[0]`: `TOPOLOGY=quad FAMILY=bernstein ORDER=[2,2] PARAM_DOMAIN=[0,1]
RATIONAL=1 NUM_CTRL=9 NUM_GP=9`; `GP_PARAM`/`GP_WEIGHT` from `integrationPoints` /
`integrationWeights`; `CTRL_WEIGHT[nElem×9]` from `controlPointWeights`. Reader maps
via the rational Bernstein path in §3.3 — no new reader class.

**(c) Belytschko co-rotational resultant beam.**
`…-BelytschkoBeam[0]`: `TOPOLOGY=line FAMILY=lagrange ORDER=[1] PARAM_DOMAIN=[-1,1]
NUM_GP=1` (or as declared). Section under `SECTION_ASSIGNMENTS` with `KIND=resultant`,
`RESULTANTS="N,Vy,Vz,T,My,Mz"`, no fiber arrays. Per-step frame under
`RESULTS/ON_ELEMENTS/LOCAL_AXES/DATA/STEP_<k>` so resultants orient in the deformed
state.

**(d) Force-based fiber beam (Lobatto) — parity.**
`…-ForceBeamColumn3d[0]`: line/lagrange/[1], `GP_PARAM` from `integrationPoints`
(preserved if weights sum to 2 ⇒ Lobatto on [-1,1]), `GP_WEIGHT` populated (was missing
in STKO). Fiber section as today.

---

## 9. Versioning & compatibility

- **`FORMAT_VERSION=1` is FROZEN (2026-05-30).** This document is the stable v1
  contract. It is now machine-enforced: `ladruno_format.validate()` is the executable
  spec (L2) and `ladruno_format.round_trip_oracle()` is a mandatory write-time geometry
  check (D5) — both gate every PR in CI (`static-gates`, pure-Python, no build) against a
  synthetic fixture, and run on real recorder output in the per-element gates. A
  conformant writer/reader is whatever passes both.
- Additive changes (new optional groups) keep the version; any change that breaks a v1
  reader bumps it. A semantic change a v1 reader still parses bumps
  `INFO/FORMAT_MINOR_VERSION` instead (absent = 0); minor 1 = §2a (apeGmsh #2034). New element families (Bézier/IGA/Belytschko) opt in via self-declaration
  (`basisInfo`) without a version bump — that is the whole point of the inverted dispatch.
- apeGmsh and `STKO_to_python` gate on `GENERATOR`/`FORMAT_VERSION` and SHOULD import
  `ladruno_format` (validator + oracle) so both readers gate on one shared spec.
- These files are **not** STKO-desktop readable by design (D2). The parity harness
  compares *values* against a sibling legacy `.mpco`, not bytes.

## 10. Open items

> [!question]
> Exact `basisInfo` response packing — one packed response vs. several tagged responses
> (`basisFamily`, `basisOrder`, …). Lean: a small structured response to minimize
> `setResponse` round-trips. Settle when the first Bézier element is implemented.

> [!question]
> `PARAM_DOMAIN` convention per family — Bernstein/Bézier natural on `[0,1]`, Lagrange
> on `[-1,1]`, simplex on barycentric. Store per-group (done) rather than normalizing,
> so the element's native convention is preserved and the reader adapts.

> [!question]
> Control-point connectivity for IGA patches that span multiple Bézier elements
> (Bézier extraction operator `C`). v1 assumes element-local control points
> (extracted). If un-extracted NURBS patches are recorded directly, add an optional
> `EXTRACTION` dataset per group. Defer until the IGA/Bézier element design is fixed.

## Implementation log

- 2026-05-30 — Schema v1 drafted. Three decisions settled: inverted self-declaration
  dispatch (legacy table as fallback), flat `GP_PARAM`+`GP_WEIGHT` quadrature, static +
  opt-in per-step local axes. Centerpiece = self-describing `BASIS`/`QUADRATURE` so
  Bézier + Belytschko need no reader-side per-class logic. Grounded in the unused IGA
  `setResponse` hook (`IGAOrder`/`IGAKnot*`/`IGAWeight`) as the existing precedent.
- 2026-05-30 — apeGmsh consultation (skill + read of `src/apeGmsh/results/readers/
  _mpco_*.py`). apeGmsh's reader **independently validated** the v1 design (its top-5
  wishlist ≈ COLUMN_MAP + GP coords/weights + beam local axes + format version). Folded
  in: extension `.ladruno` (suffixable family); **self-contained** file (no sibling
  model.h5); per-class `LOCAL_AXES` covering beams (the #1 gap — fixes the identity-
  quaternion mask); `SECTION_MAP` O(1) index; per-stage `KIND`; beam internal forces as
  section resultants at `GP_PARAM` (kills the `localForce` station-2 sign-flip hack);
  write-time homogeneity enforcement; `INFO/PROVENANCE` lineage slot; deterministic
  close + periodic flush + explicit flush control (synced-FS / crash robustness);
  pinned `<stem>.part-<N>.ladruno` 0-based-contiguous partition naming.
- 2026-05-30 — **D5 shared validator + round-trip oracle + v1 FREEZE.** `ladruno_basis`
  gained an importable `shape_functions(topology, num_nodes, xi)` (the independent
  reimplementation of C++ `computeGlobalGP`, moved out of `standard_quad_check`'s inline
  copy — Step E). `ladruno_format` gained `round_trip_oracle(path)`: for every element
  group with GP_PARAM + GLOBAL_GP_COORDS it reconstructs `x(GP_PARAM[k])` from
  MODEL/NODES + CONNECTIVITY and compares to the recorder's static `GLOBAL_GP_COORDS`
  (≤1e-12) — catching GP-ordering/basis bugs at write time. `make_synthetic` now emits
  `NDIR` + `GLOBAL_GP_COORDS` on every group and a `bary`-domain simplex (tri, NDIR=2 vs
  ORDER=(1,)) exercising finding (f). Both validator + oracle gate every PR in CI
  (`ladruno.yml` `static-gates`, pure-Python, no build) and run on real recorder output
  in `standard_quad_check`. Verified: harness 22/22; oracle clean on real
  quad/tri/brick/quad9/tet10/hex20. `FORMAT_VERSION=1` declared **frozen** (§9).
- 2026-05-31 — **Parallel per-partition output (review finding (c)).** The recorder
  now writes `<stem>.part-<N>.ladruno` (one file per process) + an INFO partition
  manifest (`PARTITIONED`/`PARTITION_ID`/`NUM_PARTITIONS`, §2), so N ranks no longer
  race to write one file (the corruption hazard). Implemented `sendSelf`/`recvSelf`
  (were bare `return 0`) to broadcast the recorder *config* P0->workers (filename,
  requests, frequency, region/sets, kind) for the PartitionedDomain/OpenSeesMP path —
  p_id carried as P0's send counter. The recorder lives in the shared sequentially-
  compiled `OPS_Recorder` lib (no `_PARALLEL_PROCESSING`, no MPI link — the frozen
  recorder's lone Allreduce is likewise compiled out here), so the partition index is
  resolved WITHOUT MPI: from the broadcast counter, else from the launcher's
  `PMI_RANK`/`PMI_SIZE` env (interpreter-per-rank `openseesmp`). **Verified on a real
  4-rank `mpiexec` openseesmp run:** 4 contiguous `mp_out.part-0..3.ladruno`, zero
  HDF5 lock errors (vs the pre-fix single-file collision), each schema-valid with
  PARTITION_ID==rank / NUM_PARTITIONS==4. Sequential output unchanged (no PMI env ->
  verbatim; 80/80,96/96,144/144,108/108, bezier, standard_quad, localaxes, pytest 23/23).
  Documented limits: the PartitionedDomain `partition` command is METIS-4-blocked in
  this build so that path is correct-by-construction (mirrors frozen) but not runtime-
  tested here; cross-rank MODEL_STAGE stamp sync (frozen's Allreduce) is absent in this
  shared-TU build (single-stage unaffected).
- 2026-05-31 — **Envelopes wired (v3a, ADR D7) — `-envelope` flag.** The built-but-
  unused `EnvelopeSink` is now reachable: the channel `sink` fields became
  `ResultSink*` and, under `-envelope`, every non-modal node/element/domain channel
  is built as an `EnvelopeSink` (componentwise MIN/MAX/ABSMAX + ARG_STEP under
  `RESULTS/ENVELOPES/<family>/<name>`, §7.4) instead of a `StreamingSink`. A
  `finalizeAllSinks()` runs each recorded step (periodic in-place rewrite for crash
  safety) and at teardown (EnvelopeSink defers all output to `finalize()`).
  `ladruno_format` gained a `envelopes()` reader + `_validate_envelopes` (shape +
  `MIN<=MAX` + `ABSMAX==max(|MIN|,|MAX|)`). **Gate (`envelope_{model,check}`):** one
  transient with BOTH a time-series and an `-envelope` recorder on the same model
  (identical commitTags); the checker reduces the time series to componentwise
  extremes and matches the ENVELOPES datasets — **ALL PASS (MIN/MAX/ABSMAX/ARG_STEP
  to 1e-9)**. Pytest 24/24. No regression in time-series mode
  (80/80,96/96,144/144,108/108, bezier, standard_quad, localaxes). Deferred: v3b
  reduced (reaction/global) envelopes need the Allreduce (dead code in this TU);
  element-envelope COLUMN_MAP component naming (a v1 follow-up, ties into finding (h)).
- 2026-05-31 — **Finding (h) resolved — authoritative names in-file + self-describing
  envelopes.** The adversarial review's finding (h) ("apeGmsh RESPONSE_CATALOG still
  supplies canonical names/counts/GP-walk-order") is closed by the cumulative work: the
  `.ladruno` format now carries everything the catalog supplied for MPCO — natural GP
  coords (`GP_PARAM`+`GLOBAL_GP_COORDS`, Steps B/D), per-result component names
  (node results' `COMPONENTS` attr e.g. `Ux,Uy,Uz`; element results' `COLUMN_MAP/
  COMP_NAMES` e.g. `sigma11,sigma22,sigma12` from the element's own `ResponseType`
  tags), component counts, and GP-walk-order (`COLUMN_MAP/GAUSS_ID` + GP_PARAM order).
  The one gap — envelope output had bare arrays with no names — is fixed here:
  `EnvelopeSink` now caches the source schema metadata and writes the **same
  self-describing result-group attrs** (`COMPONENTS`/`DISPLAY_NAME`/`DIMENSION`/…) as
  the time-series `StreamingSink`. The validator checks `COMPONENTS` count == nComp;
  `envelope_check` verifies the names round-trip (envelope `COMPONENTS` == time-series
  `COMPONENTS`). The catalog (≈2500 lines) is now **optional for `.ladruno`** (its MPCO
  role was chiefly the GP coords MPCO omits). Remaining naming follow-up: a per-column
  `COLUMN_MAP` on element ENVELOPES (gauss/section/fiber structure, not just result-
  level COMPONENTS).
- 2026-10-15 — **apeGmsh #2034: per-stage `KIND` and stage barriers (`FORMAT_MINOR_VERSION=1`).**
  `KIND` is resolved per stage from the running analysis (`OPS_GetStaticAnalysis` /
  `OPS_GetTransientAnalysis`), with eigen promotion first and `-kind` as an override.
  A `MODEL_STAGE` now also opens on a user `domainChange` (an interpreter-side barrier
  counter, `ladrunoNoteUserStageBarrier()`) and on a static↔transient flip; internal
  re-stamps stay folded (WP-165). Stage ids stay strictly increasing when no stamp moved.
  Gate: `ladruno_recorder_tests/stage_kind_check.py` (13 Tcl + 7 openseespy cases,
  Tcl and openseespy stage lists asserted identical). Review follow-up: the flip
  splits under `-kind` too (a label only), and the classic-Tcl `eigen` analysis does
  not count as running until `analysis`/`analyze`.
- 2026-10-15 — **apeGmsh #2079: consumer-written `INFO/PROVENANCE` (§7b.1).** `.ladruno`
  is apeGmsh's canonical results format. apeGmsh may stamp six lowercase pairing
  attributes into the recorder's empty `PROVENANCE` group, on a closed file only. The fork
  owns the format; the consumer owns only those attributes. Doc only: no recorder change,
  no version bump.
- 2026-10-18 — **apeGmsh #2128: kinematic node result names and labels (§7.1.1).** The
  schema documents the exact `DISPLACEMENT` / `ROTATION` / `MODES_OF_VIBRATION(U|R)` group
  names and their `COMPONENTS` per ndm, from `Ladruno_NodeResults.cpp`. This replaces the
  "e.g. `Ux,Uy,Uz`" example in the 2026-05-31 line as the reference. Doc only: no recorder
  change, no version bump.
- 2026-10-19 — **apeGmsh G10: two more consumer `PROVENANCE` attributes (§7b.1).**
  `opensees_backend` (`ladruno`/`stock`/`unknown`) and `opensees_build_sha` (the
  `ladrunoBuild` hash). Both are optional, so older files stay valid. Doc only: no
  recorder change, no version bump.
- 2026-11-14 — **WP-233 (#1082): a killed run (§7b.2).** The recorder writes under HDF5
  SWMR (latest format), the time-series `DATA` layout changes so that SWMR does not bloat
  the file, and every reader applies the reader rule and the crash rule of §7b.2. The
  contents of a closed file are unchanged: same groups, datasets, attributes and dtypes.
  Readers need HDF5 ≥ 1.10. No version bump (owner, #1082): a killed file must be opened
  before `INFO` can be read, so the rule cannot be gated on an attribute.
- 2026-10-19 — **apeGmsh #2214: every `ON_NODES` group's name and labels (§7.1.2).**
  Adds velocity, acceleration, reaction, Rayleigh, unbalanced and tie-force groups, the
  angular groups, `PRESSURE`, `overlayPressure_<tag>` and the `-NS` sensitivity groups,
  each with its writer from `Ladruno_NodeResults.cpp` / `Ladruno_OverlayResults.cpp`.
  Doc only: no recorder change, no version bump.
