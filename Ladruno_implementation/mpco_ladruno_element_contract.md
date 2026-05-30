---
title: MPCO_Ladruno element contract
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - recorder
  - element
  - contract
  - spec
---

# MPCO_Ladruno element contract

The **element-side API** that `MPCORecorderLadruno` consumes. An element that honors
this contract is recorded — including its geometry, quadrature, local frame, and
results — with **zero edits to the recorder**. This is the bridge between the planned
element work (second-order/Bézier, Belytschko beams) and the recorder.

Companions: the schema [[mpco_ladruno_schema_v1]] (what lands on disk) and the ADR
[[03_mpco_ladruno]] (why). Where the schema says "the element declares X via
`setResponse`", **this** file is the exact signature.

> **Design rule (from the ADR's inverted dispatch):** the recorder asks; the element
> answers. New element classes are never added to a recorder-side dispatch table — they
> implement the `setResponse` keywords below. The recorder falls back to a legacy
> class-tag table only for elements that predate this contract.

All signatures use the standard OpenSees element response pair:

```cpp
Response* Element::setResponse(const char **argv, int argc, OPS_Stream &output);
int       Element::getResponse(int responseID, Information &eleInfo);
```

and the `OPS_Stream` tag protocol (verified `SRC/handler/OPS_Stream.h:48-54`):

```cpp
int tag(const char *name);                 int tag(const char *name, const char *val);
int endTag();
int attr(const char *name, int value);     int attr(const char *name, double value);
int attr(const char *name, const char *value);
```

Two transport mechanisms are used, by data kind:
- **Metadata** (strings, scalars) → emitted as `tag`/`attr` calls on the passed
  `output` stream during the `setResponse` call (the recorder's stream harvests them).
- **Bulk arrays** (coords, weights, frames) → returned as an `ElementResponse` whose
  `getResponse` fills a `Vector`/`Matrix` (the `ForceBeamColumn::"integrationPoints"`
  pattern).

---

## 0. Probe lifecycle — when the recorder calls what

| Phase | Recorder calls | Element answers with |
|---|---|---|
| **Setup** (`initElementRecorders`) | `"basisInfo"` | geometry metadata (Part A.1) |
| Setup | `"integrationPoints"` | `GP_PARAM` array (A.2) |
| Setup | `"integrationWeights"` | `GP_WEIGHT` array (A.3) |
| Setup | `"controlPointWeights"` | `CTRL_WEIGHT` array — rational only (A.4) |
| Setup | `"localAxes"` | reference frame (Part B.1) |
| Setup | the user's `-E` result tokens | result descriptor tree (Part C) |
| **Per step** (`recordResultsOnElements`) | the cached result `Response*` | `getResponse` fills values |
| Per step | `"localAxes"` (if per-step frame enabled) | current frame (B.2) |

**Precedence / fallback** (matches schema §3.1): if `"basisInfo"` returns null, the
recorder derives `BASIS`/`QUADRATURE` from the legacy class-tag table assuming a
Lagrange basis + standard-rule abscissae/weights. So **standard Lagrange elements need
implement nothing new**; only non-standard geometry (Bézier/IGA) and frame-bearing
elements (beams/shells) must opt in.

---

## Part A — Geometry self-declaration

Required for any element whose geometry/quadrature is **not** a standard low-order
Lagrange rule already in the legacy table (i.e. all Bézier/Bernstein/NURBS/IGA, and any
higher-order or custom-quadrature element).

### A.1 `"basisInfo"` — the geometry descriptor (metadata via stream)

The element emits one `ElementBasis` tag carrying the full `BASIS` block, then returns a
**non-null sentinel** response so the recorder knows the probe succeeded:

```cpp
} else if (strcmp(argv[0], "basisInfo") == 0) {
    output.tag("ElementBasis");
    output.attr("topology",    "quad");        // line|tri|quad|tet|hex|wedge|pyramid|custom
    output.attr("family",      "bernstein");   // lagrange|serendipity|bernstein|nurbs|custom
    output.attr("paramDomain", "[0,1]");       // "[-1,1]" | "[0,1]" | "bary"
    output.attr("rational",    1);             // 0|1  (1 ⇒ controlPointWeights present)
    output.attr("numCtrl",     9);             // control points / nodes per element
    output.attr("numGP",       numGP);         // Gauss points per element
    output.attr("orderU",      2);             // polynomial order, parametric dir 1
    output.attr("orderV",      2);             // dir 2 (omit for line)
    output.attr("orderW",      0);             // dir 3 (omit for line/surface)
    output.endTag();
    theResponse = new ElementResponse(this, RESP_BASIS_INFO, ID(1)); // sentinel; payload ignored
}
```

The recorder reads the attrs off its capturing stream and discards the sentinel. A
`return 0` (null) here means "I have no special basis" → recorder uses the legacy
fallback. **Attr names are the contract** — spell them exactly.

> **Reserved keywords / pure-emit rule.** `"basisInfo"`, `"integrationPoints"`,
> `"integrationWeights"`, and `"controlPointWeights"` are **reserved for the recorder** —
> it is the sole caller. The `"basisInfo"` branch must be **pure-emit**: tag/attr calls
> only, no element-state mutation, so it is safe to invoke at setup and idempotent. The
> stream is the **only** metadata channel; an element that mutates state here (or relies
> on a side effect of the call) is non-conformant, because the recorder may probe with a
> capturing stream that performs no downstream work.

### A.2 `"integrationPoints"` — `GP_PARAM` (bulk array)

Parametric coordinates of every Gauss point, in the `paramDomain` declared above. Return
a `Matrix(numGP, ndir)` (a `Vector(numGP)` is accepted for `ndir == 1`):

```cpp
} else if (strcmp(argv[0], "integrationPoints") == 0) {
    theResponse = new ElementResponse(this, RESP_INT_POINTS, Matrix(numGP, ndir));
}
// getResponse:
if (responseID == RESP_INT_POINTS) {
    Matrix xi(numGP, ndir);
    for (int g = 0; g < numGP; ++g)
        for (int d = 0; d < ndir; ++d)
            xi(g, d) = gpParam_[g][d];     // in paramDomain, NOT global coords
    return eleInfo.setMatrix(xi);
}
```

`ForceBeamColumn*` already implement this (1-D `Vector`). Reuse that path.

### A.3 `"integrationWeights"` — `GP_WEIGHT` (bulk array)

```cpp
} else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, RESP_INT_WEIGHTS, Vector(numGP));
// getResponse → eleInfo.setVector( w )   // length numGP, the quadrature weights
```

Strongly recommended (enables volume/energy integrals; STKO MPCO dropped weights). If
absent for a standard rule, the recorder fills them from the known abscissae.

### A.4 `"controlPointWeights"` — `CTRL_WEIGHT` (rational only)

Required **iff** `basisInfo` declared `rational=1` (NURBS / rational Bézier). One weight
per control point, in the element's local control-point order:

```cpp
} else if (strcmp(argv[0], "controlPointWeights") == 0)
    theResponse = new ElementResponse(this, RESP_CTRL_WEIGHTS, Vector(numCtrl));
// getResponse → eleInfo.setVector( w )   // length numCtrl, the rational weights wᵢ
```

The reader's geometric map then uses `Rᵢ(ξ) = wᵢ Bᵢ(ξ) / Σⱼ wⱼ Bⱼ(ξ)` (schema §3.3).

---

## Part B — Local frame (fixes the #1 apeGmsh gap)

MPCO writes **no** beam local axes, so apeGmsh can't orient line/section force diagrams.
Any element whose results are direction-dependent (beams, shells) **must** answer
`"localAxes"`. Belytschko / co-rotational beams must also enable the per-step frame.

### B.1 `"localAxes"` — reference frame (static)

Return the 9 direction cosines `[vx(3) | vy(3) | vz(3)]` (existing OpenSees convention;
the recorder converts to the on-disk quaternion `FRAME`):

```cpp
} else if (strcmp(argv[0], "localAxes") == 0)
    theResponse = new ElementResponse(this, RESP_LOCAL_AXES, Vector(9));
// getResponse → eleInfo.setVector( {vx[0..2], vy[0..2], vz[0..2]} )
```

For frame elements this is exactly the `CrdTransf` basis — delegate to it. **Do not
return identity unless the frame truly is global**: the recorder treats a missing frame
as an error, not a silent default (apeGmsh's identity fallback masks lost data).

### B.2 Per-step frame (co-rotational / large rotation)

For co-rotational formulations the frame evolves. Declare it once so the recorder enables
a per-step `LOCAL_AXES` time series, then answer `"localAxes"` with the **current**
(deformed-config) basis on every step:

```cpp
} else if (strcmp(argv[0], "basisInfo") == 0) {
    ...
    output.attr("frameTimeVarying", 1);   // ⇒ recorder records LOCAL_AXES per step
    ...
}
```

The recorder re-queries `"localAxes"` each `record()` and writes
`RESULTS/ON_ELEMENTS/LOCAL_AXES/DATA/STEP_<k>` (schema §5). Elements with a fixed frame
omit the attr and are recorded once in `MODEL/LOCAL_AXES`.

---

## Part C — Result output protocol (the descriptor tree)

This is the **existing** MPCO `setResponse` descriptor protocol, documented here as the
contract new elements follow for stress/force/fiber results. The recorder builds the
`COLUMN_MAP` (schema §7.2) directly from these tag/attr calls. Emit a nested tree and
close every level with `endTag()`:

```cpp
output.tag("ElementOutput");
output.attr("eleTag", this->getTag());
// (node attrs node1..nodeK optional)

for (int g = 0; g < numGP; ++g) {
    output.tag("GaussPoint");
    output.attr("number", g + 1);          // 1-BASED; recorder stores g (0-based)
    output.attr("eta",    gpParam_[g][0]); // optional: parametric coord (custom rules)
    output.attr("weight", gpWeight_[g]);   // optional

    output.tag("SectionOutput");           // a section at this GP
    output.attr("secTag", secTag);

      // --- fiber section: one FiberOutput per fiber ---
      output.tag("FiberOutput");
      output.attr("yLoc", y); output.attr("zLoc", z); output.attr("area", A);
      output.attr("matTag", matTag);
        output.tag("ResponseType", "sigma");   // component name(s)
        output.tag("ResponseType", "epsilon");
      output.endTag();  // FiberOutput

    output.endTag();    // SectionOutput
    output.endTag();    // GaussPoint
}
output.endTag();        // ElementOutput
```

Level tags the recorder understands (→ `COLUMN_MAP/LEVELS` codes):
`ElementOutput`(0) · `GaussPoint`(1) · `SectionOutput`(2) · `FiberOutput`(3) ·
`NdMaterialOutput`/`UniaxialMaterialOutput`(4). Component names come from
`tag("ResponseType", <name>)`. **Use the canonical component names** that map to
apeGmsh's vocabulary (e.g. `sigma_xx`, `axial_force`, `bending_moment_z`) — see §C.3.

### C.1 Homogeneity (hard requirement)

Every element the recorder groups into one bucket must emit an **identical** descriptor
shape (same GP count, same section/fiber structure, same components). Elements that vary
(e.g. different fiber counts) are auto-split into separate buckets by the recorder — but
within a `setResponse` call the tree must be internally consistent. The recorder
**fails loud** on a ragged tree rather than emitting a bucket apeGmsh silently skips.

### C.2 Section `KIND` — fiber vs resultant (Belytschko)

A **resultant** (non-fiber) section — Belytschko beams, elastic sections — declares
itself so the recorder writes `SECTION_ASSIGNMENTS/.../KIND="resultant"` with no fiber
arrays (schema §6):

```cpp
output.tag("SectionOutput");
output.attr("secTag", secTag);
output.attr("kind",   "resultant");        // vs "fiber" (default if omitted)
output.attr("resultants", "N,Vy,Vz,T,My,Mz");
  output.tag("ResponseType", "N");
  output.tag("ResponseType", "Mz"); ...
output.endTag();
```

### C.3 Beam internal-force sign convention

Beam section/station forces are recorded as **internal section forces** at the
`GP_PARAM` stations, with **one consistent sign convention across all stations** (tension
+, sagging + per the section's local axes). This replaces MPCO's `localForce` end-force
path and its hardcoded station-2 sign flip (valid only for 2 stations). Canonical
component names: `axial_force`, `shear_y`, `shear_z`, `torsion`, `bending_moment_y`,
`bending_moment_z` (and conjugate strains `axial_strain`, `curvature_y/z`).

---

## Response-ID allocation

Use a private contiguous block of `responseID`s per element (the IGA hook used 2–6).
Keep them stable once shipped (they're not serialized, but downstream code probes by
keyword, so the integer is element-internal). Suggested local constants:

```cpp
enum { RESP_BASIS_INFO = 1001, RESP_INT_POINTS, RESP_INT_WEIGHTS,
       RESP_CTRL_WEIGHTS, RESP_LOCAL_AXES /* , element results 1.. */ };
```

Element **result** IDs (stress, force, …) stay in the element's existing scheme; the
geometry/frame probes above are additive and independent.

---

## Ordering convention (normative)

The reader evaluates the basis **once per family** (`R = B(ξ; FAMILY, ORDER)`, schema
§3.3) and forms `x = R @ X` with `X = COORDINATES[CONNECTIVITY[e, 1:]]`. For that dot
product to be correct, basis value `R[i]` must line up with control point `i` in
connectivity order. The contract pins that alignment so **no permutation is transported**.
The rule splits on the *shape of the parametric domain*, not just the family — a Bézier
**tet/tri** is a simplex and is **not** indexed like a tensor-product quad/hex:

- **Tensor-product topologies (`line`, `quad`, `hex`) with `family ∈ {bernstein, nurbs}`**
  — control points in `CONNECTIVITY`, and the rows of `GP_PARAM` / `CTRL_WEIGHT`, are in
  **lexicographic tensor order, fastest in U, then V, then W**:

  ```
  index(i,j,k) = i + (orderU+1)·j + (orderU+1)(orderV+1)·k     # i,j,k = 0..order_dir
  ```

- **Simplex topologies (`tri`, `tet`) with `family = bernstein`** — control points are a
  barycentric multi-index set `α`, `|α| = n` (count `C(n+d,d)`, **not** `(n+1)^d`), but
  the **emission order is the standard FE hierarchical order: all vertices first, then
  mid-edge nodes (then face, then interior for higher order)** — *not* multi-index
  lexicographic. This matches the reference Bézier elements (Kadapa 2018, IJNME 117:543)
  and is **byte-identical to the existing quadratic-Lagrange T6/T10 node order**, so the
  reader reuses the Lagrange simplex node order and only swaps the basis to Bernstein.

  For the quadratic triangle (`BezierTri6`, the v1 reference element) the pinned order is:

  ```
  P1 = ξ₃² (corner)   P2 = ξ₁² (corner)   P3 = ξ₂² (corner)      # vertices
  P4 = 2ξ₁ξ₃ (edge 1–2)  P5 = 2ξ₁ξ₂ (edge 2–3)  P6 = 2ξ₂ξ₃ (edge 3–1)   # mid-edges
  ```

  with edge order (1-2, 2-3, 3-1) matching the standard T6. The Bézier **tet** (Tet10,
  Kadapa §5) follows the same rule: 4 vertices, then 6 mid-edges in the conventional
  quadratic-tet edge order. `GP_PARAM` rows are the rule's points in the element's GP
  loop order (`PARAM_DOMAIN="bary"`, stored as the `d` free area coords); `CTRL_WEIGHT`
  is per control point in the same vertices-then-edges order (`=1` for non-rational
  Bézier — `BezierTri6` is `rational=0`, so `CTRL_WEIGHT` is omitted entirely).

- **`wedge` / `pyramid`** — mixed tensor⊗simplex; defer the pinned order until such an
  element is on the table (no Bézier wedge/pyramid is planned for v1).

- **`family = lagrange`** — retains the legacy OpenSees node numbering; the reader's
  `lagrange` evaluator already matches it via the legacy class-tag table. New elements
  should prefer `bernstein` if they want the pinned order.

In all cases the element is responsible for emitting `CONNECTIVITY` (via the recorder's
geometry pass) and `CTRL_WEIGHT`/`GP_PARAM` (via the responses below) **in the same
order** the reader's `B(ξ; family, order, topology)` enumerates. A mismatch is silent
(wrong `x_global`, no exception) — this is the single most important invariant for a
custom-geometry element, and the simplex/tensor split above is exactly where it bites.

---

## Worked skeletons

### (a) Quadratic rational Bézier quad

```cpp
Response* BezierQuad::setResponse(const char **argv, int argc, OPS_Stream &out) {
  if (strcmp(argv[0], "basisInfo") == 0) {
    out.tag("ElementBasis");
    out.attr("topology","quad"); out.attr("family","bernstein");
    out.attr("paramDomain","[0,1]"); out.attr("rational",1);
    out.attr("numCtrl",9); out.attr("numGP",numGP_);
    out.attr("orderU",2); out.attr("orderV",2);
    out.endTag();
    return new ElementResponse(this, RESP_BASIS_INFO, ID(1));
  }
  if (strcmp(argv[0],"integrationPoints")==0)
    return new ElementResponse(this, RESP_INT_POINTS, Matrix(numGP_,2));
  if (strcmp(argv[0],"integrationWeights")==0)
    return new ElementResponse(this, RESP_INT_WEIGHTS, Vector(numGP_));
  if (strcmp(argv[0],"controlPointWeights")==0)
    return new ElementResponse(this, RESP_CTRL_WEIGHTS, Vector(9));
  // ... element results (stress/strain) via the Part C tree ...
}

int BezierQuad::getResponse(int responseID, Information &eleInfo) {
  if (responseID == RESP_INT_POINTS) {        // GP_PARAM, lexicographic in (u,v)
    Matrix xi(numGP_, 2);
    for (int g = 0; g < numGP_; ++g) { xi(g,0) = gpU_[g]; xi(g,1) = gpV_[g]; }
    return eleInfo.setMatrix(xi);
  }
  if (responseID == RESP_INT_WEIGHTS)         // GP_WEIGHT, same GP order as GP_PARAM
    return eleInfo.setVector(gpW_);
  if (responseID == RESP_CTRL_WEIGHTS)        // CTRL_WEIGHT, lexicographic ctrl order
    return eleInfo.setVector(ctrlW_);         // MUST match CONNECTIVITY emission order
  // ... element results ...
  return -1;
}
```

The `MUST match CONNECTIVITY` comment is load-bearing: `ctrlW_` and the connectivity the
element reports to the recorder are both bound to the lexicographic order pinned above.

### (b) Belytschko co-rotational resultant beam

```cpp
Response* BelytschkoBeam::setResponse(const char **argv, int argc, OPS_Stream &out) {
  if (strcmp(argv[0],"basisInfo")==0) {
    out.tag("ElementBasis");
    out.attr("topology","line"); out.attr("family","lagrange");
    out.attr("paramDomain","[-1,1]"); out.attr("rational",0);
    out.attr("numCtrl",2); out.attr("numGP",numGP_);
    out.attr("orderU",1);
    out.attr("frameTimeVarying",1);          // co-rotational ⇒ per-step frame
    out.endTag();
    return new ElementResponse(this, RESP_BASIS_INFO, ID(1));
  }
  if (strcmp(argv[0],"localAxes")==0)         // current basis from the co-rot transform
    return new ElementResponse(this, RESP_LOCAL_AXES, Vector(9));
  if (strcmp(argv[0],"force")==0 || strcmp(argv[0],"section")==0) {
    out.tag("ElementOutput"); out.attr("eleTag", getTag());
    for (int g=0; g<numGP_; ++g) {
      out.tag("GaussPoint"); out.attr("number", g+1); out.attr("eta", gpParam_[g]);
      out.tag("SectionOutput"); out.attr("secTag", secTag_);
      out.attr("kind","resultant"); out.attr("resultants","N,Vy,Vz,T,My,Mz");
      out.tag("ResponseType","axial_force"); /* ... */
      out.endTag(); out.endTag();
    }
    out.endTag();
    return new ElementResponse(this, /*result id*/ 1, Vector(6*numGP_));
  }
}
```

---

## Conformance checklist

What an element must implement, by category:

| Element kind | basisInfo | int. points/weights | ctrl weights | localAxes | per-step frame | result tree |
|---|---|---|---|---|---|---|
| Standard Lagrange solid/shell | optional¹ | optional¹ | — | if framed | — | required (Part C) |
| Higher-order / Bézier / IGA | **required** | **required** | if rational | if framed | — | required |
| NURBS / rational Bézier | **required** | **required** | **required** | if framed | — | required |
| Belytschko / co-rotational beam | recommended² | recommended² | — | **required** | **required** | required (resultant) |

¹ Covered by the legacy class-tag fallback; implement only to override the derived basis.
² Line/Lagrange is derivable from the legacy table; declare explicitly to also carry
`frameTimeVarying`.

---

## Settled rules

- **`integrationPoints` is authoritative (was an open question).** When an element
  implements `basisInfo` *and* falls under a standard rule in the legacy table, the
  element's `integrationPoints`/`integrationWeights` **win**; the recorder validates
  against the derived abscissae and **warns** (does not fail) on mismatch. The element is
  the source of truth for its own quadrature.
- **Ordering (tensor)** — pinned lexicographic U→V→W for `bernstein`/`nurbs` on
  `line`/`quad`/`hex` (see *Ordering convention* above). The former GP-flattening open
  item is closed for tensor-product topologies by that rule.
- **`basisInfo` sentinel** — keep the non-null `ID(1)` sentinel; it is cheap and
  unambiguous, and the *Reserved keywords / pure-emit rule* (A.1) makes the metadata
  channel well-defined regardless of how the recorder's capturing stream is implemented.

## Open items

> [!note]
> **Simplex ordering — RESOLVED (2026-05-30).** Confirmed against the reference element
> (`bezierFEM/opensees implementation/BezierTri6.cpp`, Kadapa 2018): simplex Bézier
> elements use the **standard FE hierarchical order (vertices then mid-edges)**, *not*
> multi-index lexicographic. The contract is pinned to that order above; it coincides with
> the legacy Lagrange T6/T10 numbering, so no `CTRL_ORDER` hint is needed.

## Implementation log

- 2026-05-30 — Contract drafted. Three parts: (A) geometry self-declaration
  (`basisInfo` metadata stream + `integrationPoints`/`integrationWeights`/
  `controlPointWeights` numeric responses), (B) local frame (`localAxes` static +
  `frameTimeVarying` per-step), (C) the existing descriptor-tree result protocol +
  section `KIND` + beam sign convention. Grounded in `OPS_Stream.h:48-54` and the IGA /
  `ForceBeamColumn::integrationPoints` precedents. Worked skeletons for a rational Bézier
  quad and a Belytschko co-rotational beam.
- 2026-05-30 — Contract review pass (recorder-side). Closed the geometry-ordering hole:
  added the normative **Ordering convention** (lexicographic U→V→W for bernstein/nurbs,
  mandated so no permutation is transported and `R[i]`↔ctrl-point `i` is guaranteed).
  Added a complete `BezierQuad::getResponse` showing where `GP_PARAM`/`GP_WEIGHT`/
  `CTRL_WEIGHT` are filled in that order. Added the **Reserved keywords / pure-emit rule**
  to A.1 (recorder is sole caller; `basisInfo` branch mutates no state). Promoted the
  three open questions to **Settled rules**: `integrationPoints` authoritative
  (warn-on-mismatch), ordering pinned, sentinel kept. Schema §3 datasets unchanged.
- 2026-05-30 — **Simplex ordering fix + resolution.** The first ordering rule was
  tensor-product only (`index = i + (orderU+1)j + …`) and silently wrong for `tet`/`tri`
  (simplices, count `C(n+d,d)` not `(n+1)^d`). Split *Ordering convention* into tensor
  (line/quad/hex) vs simplex (tri/tet) vs deferred wedge/pyramid. **Verified against the
  actual reference element** `bezierFEM/.../BezierTri6.cpp` + `bezier_tri6_element.py`
  (Kadapa 2018, IJNME 117:543): the element uses the **standard FE hierarchical order
  (vertices then mid-edges)** — N₁=ξ₃², N₂=ξ₁², N₃=ξ₂² corners; N₄/N₅/N₆ mid-edges —
  **not** multi-index lexicographic (my initial assumption, now corrected). This matches
  legacy Lagrange T6/T10 numbering ⇒ reader reuses it, swaps only the basis. `BezierTri6`
  is `rational=0` so `CTRL_WEIGHT` is omitted. Open item closed.
