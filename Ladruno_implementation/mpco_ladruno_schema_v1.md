---
title: MPCO_Ladruno HDF5 schema v1
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

# MPCO_Ladruno HDF5 schema v1

The on-disk contract for the `mpcoLadruno` recorder. apeGmsh is the **canonical
consumer** and co-owns this spec; `STKO_to_python` gets a parallel reader branch.
Companion: the ADR [[03_mpco_ladruno]]. This file is the layout `mpcol::` writes and
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
  `^(?P<stem>.+?)\.part-(?P<idx>\d+)\.ladruno$` and *errors* on a gap.

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
3. **Versioned.** `INFO/GENERATOR="MPCO_Ladruno"` + `INFO/FORMAT_VERSION=1`. Readers
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
│   ├── GENERATOR        (string)  "MPCO_Ladruno"
│   ├── FORMAT_VERSION   (int)     1
│   ├── SOLVER_NAME      (string)  "OpenSees"
│   ├── SOLVER_VERSION   (int[3])
│   └── SPATIAL_DIM      (int)     1 | 2 | 3
│
└── MODEL_STAGE[<stamp>]                         (one per Domain::hasDomainChanged())
    ├── attrs STEP(int), TIME(double), KIND(string "transient"|"static"|"eigen")
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
| `GENERATOR` | string | `"MPCO_Ladruno"` | **reader branch key** — distinguishes from STKO `.mpco` |
| `FORMAT_VERSION` | int | `1` | bump on any breaking layout change |
| `SOLVER_NAME` | string | `"OpenSees"` | |
| `SOLVER_VERSION` | int[3] | e.g. `[3,5,1]` | `OPS_VERSION` split on `.` |
| `SPATIAL_DIM` | int | 1/2/3 | from first node, enforced uniform |

A reader that doesn't find `GENERATOR=="MPCO_Ladruno"` must treat the file as a
legacy STKO `.mpco` (or refuse it).

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
> skeletons) are specified in [[mpco_ladruno_element_contract]]. This section is the
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
> or forces the native `model.h5` path. MPCO_Ladruno **must** write real per-element
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
   DATA  [T × nIds × nComp] double   ← extensible on T, chunked + shuffle + deflate(4)
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
reaction. **Built only after the separate energy work lands** — MPCO_Ladruno will
*consume* `EnergyBalanceRecorder` ([[project_energy_balance_feature]]), not reinvent
it. Spec stub here for forward compatibility.

### 7.4 `ENVELOPES/ON_{NODES,ELEMENTS,DOMAIN}/<name>` — (v3)

Per the ADR (D6/D7): componentwise `MIN`/`MAX`/`ABSMAX` + `ARG_STEP`, per-`MODEL_STAGE`,
periodic in-place rewrite. v3a covers consistent quantities (zero-MPI); v3b adds
reduced (reaction/global) envelopes alongside the energy work.

```
ENVELOPES/ON_NODES/<name>/
   ID      [nN×1] int
   MIN     [nN×nComp]   MAX [nN×nComp]   ABSMAX [nN×nComp]   double
   ARG_STEP[nN×nComp]   int     step index of each extreme (commitTag, global across ranks)
```

---

## 7b. File robustness & provenance (OpenSees-side)

The apeGmsh consultation surfaced two robustness failures the recorder owns at the
source (the user's stated priority: performance + robustness on OpenSees first):

- **Deterministic close + periodic flush.** Writing an MPCO file onto a synced FS
  (OneDrive/SeaDrive/Dropbox) native-crashes the kernel at `wipe()` — the deferred
  `H5Fclose` races the sync client and leaves the HDF5 write-flag stuck;
  `HDF5_USE_FILE_LOCKING=FALSE` does not fix it. Combined with the `exit(-1)`-leaves-
  files-unflushed history ([[project_mpco_exit_crash]]), MPCO_Ladruno must: flush every
  recorded step (as MPCO does), **rewrite envelope/aggregate datasets periodically**
  (§7.4), and provide an **explicit flush/close recorder control** so a user can force a
  clean close in the same cell before copying off local disk.
- **Provenance / lineage slot.** `INFO/PROVENANCE` (group) reserves
  `MODEL_HASH`/`RESULTS_HASH` string slots so apeGmsh can chain its blake2b lineage
  (`fem_hash → model_hash → results_hash`). OpenSees stamps a stable results identity;
  apeGmsh stamps the model side. Lets the reader detect model/results drift.

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

- `FORMAT_VERSION=1` is this document. Additive changes (new optional groups) keep the
  version; any change that breaks a v1 reader bumps it.
- apeGmsh and `STKO_to_python` gate on `GENERATOR`/`FORMAT_VERSION`.
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
