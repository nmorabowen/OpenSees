---
title: MPCO_Ladruno — modular recorder fork
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - recorder
  - adr
---

# MPCO_Ladruno — modular recorder fork

> **ADR** (Architecture Decision Record). Captures the four settled decisions of
> 2026-05-29 plus the module/architecture plan they imply. Companion docs:
> [[../Ladruno_internal/01_compilation_journal]] (build), and a future
> `mpco_ladruno_schema_v1.md` (the on-disk HDF5 spec, co-designed with apeGmsh).

## What

A **sibling** HDF5 recorder forked from `SRC/recorder/MPCORecorder.cpp`, exposed as
`recorder mpcoLadruno …`. It reorganizes the 6,901-line monolith into ~10 focused
translation units (all under `namespace mpcol`), retains the full nodal + element
result functionality, and adds two genuinely new capabilities the original cannot
express:

1. **Global / domain-scalar results** (`ON_DOMAIN`) — e.g. energy balance, base
   shear, total reactions: one `[1 × nComp]` row per step.
2. **Envelope results** (`ENVELOPES`) — min / max / abs-max over time plus the step
   at which the extreme occurred, for nodes, elements, *and* domain scalars.

**In scope (this plan, v1):** the modular split at *semantic* parity with the
current recorder, in a new apeGmsh-native HDF5 schema, with a value-level parity
harness. **Out of scope (deferred):** the global/`ON_DOMAIN` channel (energy balance
is being implemented separately — MPCO_Ladruno will *consume* it, not build it, so the
whole `ON_DOMAIN` channel waits on that landing) and the envelope channel. The
envelope work itself splits: **v3a** = zero-MPI envelopes for consistent quantities
(no dependency on the energy work), **v3b** = reduced (reaction/global) envelopes,
which share the per-step reduction plumbing with the energy/`ON_DOMAIN` channel and
co-land with it. The v1 abstraction (`ResultSource`/`ResultSink`) is designed up front
so all of this is additive, not a refactor.

## Why

Four drivers, all confirmed:

- **apeGmsh canonical recorder** — the end goal is apeGmsh emitting *and* consuming
  this as its native output format (replacing the STKO-shaped `.mpco` for apeGmsh
  workflows). This is what licenses the schema divergence: apeGmsh, not STKO,
  becomes the source of truth for the layout.
- **Global / energy outputs** — energy balance ([[project_energy_balance_feature]],
  `EnergyBalanceRecorder` classTag 26) and other whole-model scalars have no home in
  MPCO today; both engines are strictly per-node / per-element.
- **Envelope outputs** — peak drift / peak stress without post-processing every step;
  requires stateful cross-step accumulation the current stateless writers can't do.
- **Maintainability** — a 6,901-line single file is hard to extend safely; the
  element engine in particular concentrates years of element-specific tribal
  knowledge (workarounds for SSPbrick, ForceBeamColumn3d, shell keyword swap).

**Hard constraint:** `MPCORecorder.cpp` is *frozen* — byte-identical to upstream
master, deliberately untouched for ~4.5 years so STKO read-compatibility never
regresses (see [[project_mpco_exit_crash]]). MPCO_Ladruno must not edit it. The
price we accept: Layers 0/1/2/4 are **copied**, not shared.

## Decisions (settled 2026-05-29)

| # | Decision | Rationale | Consequence |
|---|----------|-----------|-------------|
| **D1** | **New sibling recorder**, not in-place refactor or shared lib. | Keeps the frozen file byte-identical; STKO + upstream diff untouched. | Layers 0/1/2/4 duplicated. Acceptable. |
| **D2** | **apeGmsh-native schema** — diverge from STKO. | apeGmsh is the canonical consumer; freedom to fix MPCO's clunky parts (flattened `META` columns). | STKO desktop + current `STKO_to_python`/`Results.from_mpco` will **not** read v1 files until their readers are extended. Mitigated by a `GENERATOR`/`FORMAT_VERSION` tag so readers branch cleanly. |
| **D3** | **v1 = refactor-to-parity** (semantic, not byte). Global=v2, Envelope=v3. | De-risk: prove the modular base reproduces every existing result value before adding features. | "Parity" is verified by a value-level harness (read both files, assert equality to 1e-12), since divergent schema rules out a byte diff. |
| **D4** | **Wrap the whole module in `namespace mpcol`.** | Avoid ODR/duplicate-symbol collisions with the frozen TU (`LibraryLoader`, `libload::`, `h5::`, `H5*` macros all have external linkage). | A second HDF5 `LibraryLoader` singleton exists at runtime → a benign second `dlopen`/refcount of the HDF5 lib. Noted, accepted. |
| **D5** | **Keep MPCO's per-partition parallel model** (`part-N.mpco`, reader stitches). No parallel HDF5. | Verified MPCO's record loop only communicates one int/step (stage-stamp `Allreduce`, line 4594); result data is never communicated. Parallel HDF5 (MPI-IO collective writes) is fragile on the MS-MPI/libfabric stack ([[project_openseespymp_plan]]) and fights SWMR. | Stitching is **apeGmsh-owned** (canonical reader). `sendSelf` keeps broadcasting only the spec. |
| **D6** | **Envelope phasing: zero-MPI consistent quantities first (v3a), reduced quantities later (v3b).** | min/max commute with partition-reduction *only* for consistent quantities; kinematics + element results need **zero** new comm. Reaction-at-boundary + global envelopes need per-step reduction — same plumbing as the deferred energy/`ON_DOMAIN` channel. | v3a is shippable with no MPI risk. v3b lands with/after the energy work it shares machinery with. |
| **D7** | **Envelope semantics: componentwise min/max/abs-max + `ARG_STEP`.** | Matches design-check needs (independent peak per component); `ARG_STEP` (step of each extreme) enables concurrent-state reconstruction without extra storage. | Per-`MODEL_STAGE` reset; periodic in-place rewrite for crash safety. |

## Where

New code (all under `SRC/recorder/mpco_ladruno/`, all in `namespace mpcol`):

| File | From monolith lines | Layer |
|------|--------------------|-------|
| `MPCOL_Hdf5Library.{h,cpp}` | 200–555 (`libload`, `LibraryLoader`) | 0 — platform / runtime HDF5 load |
| `MPCOL_Hdf5.{h,cpp}` | 928–1243 (`h5`) | 1 — HDF5 C++ wrapper |
| `MPCOL_Types.h` | 556–927, 1244–1477 (utils, structs, **enums**, `ProcessInfo`, `OutputFrequency`, `Timer`) | 2 — shared vocabulary |
| `MPCOL_ResultIO.h` | *new* | 2 — `ResultSchema` / `ResultSource` / `ResultSink` interfaces |
| `MPCOL_NodeResults.{h,cpp}` | 1478–2594 (`mpco::node`) | 3a — nodal sources |
| `MPCOL_ElementResults.{h,cpp}` | 2595–4279 (`mpco::element`) | 3b — element source (the hard port) |
| `MPCOL_DomainResults.{h,cpp}` | *new* (consumes EnergyBalance logic) | 3c — global/scalar sources (**v2**) |
| `MPCOL_Sinks.{h,cpp}` | *new* (`StreamingSink`, `EnvelopeSink`) | 3d — persistence strategies (**v3** for envelope) |
| `MPCOL_Serializer.h` | 4280–4401 | 4 — parallel send/recv |
| `MPCORecorderLadruno.{h,cpp}` | 4410–6901 (`private_data`, lifecycle, `OPS_…`) | 5 — orchestration + command API |

Modify (3 external touch-points — the entire upstream-diff surface, all additive,
all under `#ifdef _HDF5`; record in [[../Ladruno_internal/01_compilation_journal]]):

- `SRC/interpreter/OpenSeesOutputCommands.cpp:132` — add
  `recordersMap.insert({"mpcoLadruno", &OPS_MPCOLadrunoRecorder});`
- `SRC/classTags.h:1207` — add `RECORDER_TAGS_MPCOLadrunoRecorder` (next free value).
- `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp:2799` — add broker case.
- `SRC/recorder/CMakeLists.txt` — add the new sources under the `_HDF5` guard.

Reference (copy patterns / verify against): `SRC/recorder/MPCORecorder.cpp` (frozen,
read-only), `EnergyBalanceRecorder` (classTag 26).

## How

### Core abstraction — split *what to compute* from *how to persist*

Today `ResultRecorder::record()` fuses two responsibilities: produce a flat buffer
for the step, and write a `STEP_<k>` dataset. That fusion is exactly why envelopes
and globals don't fit. v1 breaks them apart:

```cpp
namespace mpcol {

struct ResultSchema {                 // self-describing column metadata
    std::string name, display_name, components_csv, dimension, description;
    int num_components;
    ResultType::Enum result_type;     // Generic / Modal / ...
    ResultDataType::Enum data_type;   // Scalar / Vectorial / Tensorial / Beam / Shell
};

// WHAT to compute. node / element / domain all implement this.
class ResultSource {
public:
    virtual const ResultSchema& schema() const = 0;
    virtual const std::vector<int>& ids() const = 0;   // node tags / elem-gp ids / {0} for domain
    virtual void evaluate(const ProcessInfo&, std::vector<double>& buffer) = 0; // this step
};

// HOW to persist. swappable per source.
class ResultSink {
public:
    virtual void begin(const ProcessInfo&, const ResultSource&) = 0;            // groups + ID dataset
    virtual void accept(const ProcessInfo&, const ResultSource&,
                        const std::vector<double>& buffer) = 0;                 // one step's buffer
    virtual void finalize(const ProcessInfo&) = 0;
};

} // namespace mpcol
```

- `StreamingSink::accept` → writes `STEP_<k>` (today's behavior, exact parity).
- `EnvelopeSink` → keeps running `min`/`max`/`absmax` buffers + an `arg-step`
  buffer; `accept` updates them in place; `finalize` writes them once. **Because it
  wraps any `ResultSource`, envelopes come for free for nodes, elements, and domain
  scalars — zero per-family reimplementation.** This is the scalability win.

The three families become three `ResultSource` flavors:

- **Node sources** — the existing 26-class hierarchy, lightly adapted: the current
  `bufferResponse()` *is* `evaluate()`; the metadata block in the constructor *is*
  `schema()`. Cleanest port (the seam already exists).
- **Element source** — an *adapter* (not a rewrite) over the existing
  `OutputDescriptor` / `OutputDescriptorStream` / collection machinery:
  `schema()` derives from `OutputDescriptorHeader`; `evaluate()` runs the response
  collection into the flat buffer. **Port faithfully — keep every workaround.** This
  is the highest-risk file; treat it as a structural move, not a redesign.
- **Domain source** (v2) — `ids() == {0}`, `evaluate()` returns a handful of
  scalars. Energy balance source mirrors / consumes `EnergyBalanceRecorder`'s
  per-region computation rather than reinventing it.

### Lifecycle (unchanged shape, new wiring)

`setDomain → record → (lazy initialize: open file, write INFO, writeModel,
build sources+sinks) → on each qualifying step: for each (source,sink) pair
sink.accept(source.evaluate()) → flush`. Model-stage staging
(`hasDomainChanged()` stamp → new `MODEL_STAGE[<stamp>]` group) is preserved
verbatim. At analysis end, `finalize()` lets envelope sinks flush.

### Element-side contract

The element `setResponse` API the recorder consumes — `basisInfo` (geometry metadata),
`integrationPoints`/`integrationWeights`/`controlPointWeights` (quadrature), `localAxes`
(static + `frameTimeVarying` per-step), and the descriptor-tree result protocol — is
specified in [[mpco_ladruno_element_contract]]. It is the bridge between the Bézier /
Belytschko element roadmaps and this recorder: an element honoring it is recorded with
no recorder edits.

### Public API (command syntax)

v1 mirrors `recorder mpco` exactly (parity), with the new verb:

```
recorder mpcoLadruno <file.mpco> -N displacement reactionForce \
                                 -E force section.fiber.stress \
                                 -T dt 0.0  -R <region?>
```

Forward-looking additions (reserved now, implemented v2/v3):

```
  -G energy baseShear            # v2: global/domain scalars  -> ON_DOMAIN
  -ENV -N displacement           # v3: wrap a request in an EnvelopeSink -> ENVELOPES
```

### Schema (apeGmsh-native, versioned) — sketch

Keep MPCO's genuinely good ideas; fix its worst ergonomics.

```
/
├── INFO            attrs: GENERATOR="MPCO_Ladruno", FORMAT_VERSION=1, NDIM, UNITS, SOLVER
└── MODEL_STAGE[<stamp>]
    ├── MODEL       NODES, ELEMENTS, LOCAL_AXES, SECTION_ASSIGNMENTS, SETS   (≈ as today)
    └── RESULTS
        ├── ON_NODES/<name>/{ID, DATA/STEP_k}
        ├── ON_ELEMENTS/<name>/{ID, META, DATA/STEP_k}
        ├── ON_DOMAIN/<name>/{DATA/STEP_k}                 # NEW (v2), shape [1 × nComp]
        └── ENVELOPES/ON_{NODES,ELEMENTS,DOMAIN}/<name>/   # NEW (v3)
                {ID, MIN, MAX, ABSMAX, ARG_STEP}
```

**Keep:** self-describing groups; natural-coord `GP_X` in [-1,+1]; fiber sections in
section-local axes; per-`MODEL_STAGE` staging.
**Improve:** replace the flattened `META` columns (the `;`-separated `COMPONENTS`
string + `GAUSS_IDS`/`MULTIPLICITY`/`NUM_COMPONENTS` parallel arrays — MPCO's worst
reader ergonomics) with structured, attribute-rich per-result metadata that apeGmsh
and `STKO_to_python` can parse without the bespoke decoder.

> The concrete layout is specified in [[mpco_ladruno_schema_v1]] (drafted
> 2026-05-30), **co-designed with apeGmsh's reader**. Its centerpiece is a
> self-describing `BASIS`/`QUADRATURE` descriptor per element group (replacing the two
> closed enums) so second-order/Bézier elements and Belytschko co-rotational beams need
> no reader-side per-class logic. `FORMAT_VERSION` lets readers branch; apeGmsh owns
> the spec going forward.

### apeGmsh integration

apeGmsh's `Results.from_mpco`-equivalent gains a branch on
`INFO/GENERATOR == "MPCO_Ladruno"`. Because apeGmsh co-designs the schema, this is a
clean read path, not a reverse-engineering exercise. Longer term, apeGmsh's
exporters emit `recorder mpcoLadruno` as the default, making it the canonical
output. `STKO_to_python` gets a parallel reader branch (same `FORMAT_VERSION` gate)
so existing post-processing keeps working against the new files.

### Parallel model & envelopes

**Parallel (kept from MPCO, D5).** `sendSelf`/`recvSelf` broadcast only the recorder
*spec* from p0 to all ranks (extend the `Serializer` payload for new request types;
keep enum values stable). Each rank records its own subdomain to `part-<p_id>.mpco`.
The **only** collective op in the record loop is the existing single-int
`MPI_Allreduce` (line 4594) that synchronizes the `hasDomainChanged()` stamp so every
rank names its `MODEL_STAGE[<stamp>]` group identically — result data is never
communicated. apeGmsh stitches part files on read (it owns the contract now).

**Envelopes via a sink, not a new hierarchy.** `EnvelopeSink` wraps any
`ResultSource`; the only parallel-aware knob is a per-source
`requiresPartitionReduction()` flag. Quantities split into three classes by whether
their per-step value is consistent across the partitions that hold the entity
(min/max commute with partition-reduction *only* when it is):

| Class | Examples | Strategy | Comm |
|-------|----------|----------|------|
| **Consistent** | disp/vel/accel; element stress/strain/force (single-owner) | local accumulate → reader stitches | **none** (v3a) |
| **Additive** | reactions at partition-boundary nodes | per-step `Allreduce` → accumulate on the summed value (per-partition max is wrong: argmax steps differ, `max_t(a+b) ≠ max_t a + max_t b`) | per-step (v3b) |
| **Global** | base shear, total reaction, energy | per-step `Allreduce` of scalars → accumulate | per-step (v3b) |

The additive and global paths are the **same** reduction plumbing the deferred
energy/`ON_DOMAIN` channel needs — so v3b co-lands with that work.

**Envelope semantics (D7):** componentwise `MIN`/`MAX`/`ABSMAX` + `ARG_STEP` (step
index of each extreme; `commitTag` is global across ranks so it's a valid local id).
Envelopes are **per `MODEL_STAGE`** — flush-and-reset on `hasDomainChanged()`.
Because `finalize()` has no clean OpenSees callback (only destructor / stage change),
the sink **periodically rewrites** envelope datasets in place alongside the existing
per-step `H5Fflush`, so a kernel crash ([[project_mpco_exit_crash]]) loses at most the
last N steps of envelope.

### Testing — layered, frozen-MPCO-as-oracle

Divergent schema ⇒ no byte diff. The frozen `recorder mpco` is the **value oracle**:
both recorders run on the same domain in one process, and we compare *information, not
layout*. Five layers, fastest first; all under `Ladruno_scripts/`, run on
`opensees_venv` (h5py + `STKO_to_python`). Build (oneAPI/Conan) is the bottleneck — loop
is *build once → run all model variants → fast pytest*.

- **L0 — C++ unit tests (no run).** Pure-logic pieces the Source/Sink split isolates:
  BASIS-derivation from legacy enums, `COLUMN_MAP` build, `Serializer` round-trip,
  `EnvelopeSink` min/max, basis-function reconstruction. Fake `ResultSource` → assert
  `ResultSink` output.
- **L1 — parity gate (v1 acceptance).** One script, both recorders. A pytest normalizes
  both to `dict[(node,component,step)]→v` and
  `dict[(elem,gp,section,fiber,component,step)]→v` (legacy via `STKO_to_python`, new via
  the h5py reader), asserts equality to 1e-12. A component-name **alias map** handles
  renames (`localForce`→`axial_force`). Parity covers only the **intersection** (what
  MPCO records); new fields (weights, `GP_PARAM`, `SECTION_MAP`) are checked by L2/L3.
  Canonical model = one element per geometry/rule family (truss, dispBeamColumn fiber,
  forceBeamColumn Lobatto, ASDShellQ4, stdBrick).
- **L2 — schema validator (`validate_ladruno.py`).** Pure h5py; the **executable form of
  the spec**. Asserts every schema-v1 invariant (BASIS+QUADRATURE present;
  `NUM_COLUMNS == Σ MULTIPLICITY·NUM_COMP`; `GP_*` shapes vs `numGP`; `LOCAL_AXES` for
  framed elements; `SECTION_MAP` vs `SECTION_ASSIGNMENTS`; part-file contiguity).
  Written first; every produced file passes through it.
- **L3 — reconstruction.** Read BASIS+QUADRATURE+`GP_PARAM`, reconstruct global Gauss
  coords via `x(ξ)=ΣRᵢ(ξ)Xᵢ`, compare to independently-computed coords. Validates the
  reader contract; reused verbatim for Bézier (Bernstein) later.
- **L4 — parallel parity + robustness.** serial vs 2-rank MP (merge `part-0/1.ladruno`,
  assert == serial) validates D5 + merge + naming; plus crash-safety (kill mid-run, file
  readable to last flush) and synced-FS close.

v3 adds envelope-correctness (independent Python max vs `ENVELOPES`, `ARG_STEP` check);
the deferred energy/`ON_DOMAIN` work adds energy-conservation checks.

### Execution plan — the diamond (agent orchestration)

Clean seams make Phase 2 genuinely parallel, but the dependency graph forces a diamond:
narrow contract → wide fan-out → narrow integration. Phased **discrete agents with
build-gates run between phases** (not one autonomous workflow — compile-and-link is
unforgiving and the gates depend on the slow build).

1. **Phase 1 — narrow, sequential (linchpin).** `MPCOL_Types.h`, `MPCOL_ResultIO.h` (the
   Source/Sink contract), Layers 0/1 in `namespace mpcol`, a trivial registering
   recorder. **Gate: builds + links next to the frozen recorder** — proves the
   `namespace mpcol` ODR fix (the biggest unknown). Nothing fans out before this.
2. **Phase 2 — wide, parallel (worktrees, against frozen Phase-1 headers).** Agent B
   element-source adapter + `basisInfo` capture (long pole, start first, "evaluate don't
   transcribe" the workarounds; codex-rescue second opinion on the OutputDescriptor
   port); Agent A node-source port; Agent C sinks; Agent D (Python, **test-first**) the
   L1/L2/L3 harness against the spec.
3. **Phase 3 — narrow, sequential (integration).** Wire Layer 5 (lifecycle, `OPS_`
   parser, 3 external touch-points, CMake), build, run the harness. **Gate: parity 1e-12.**

Highest-leverage first move: **Agent D's validator + parity harness** — pure Python, no
build, turns the spec executable so Phase 1's first file is validated the moment it
exists.

## Risks / open questions

> [!question]
> **Energy / `ON_DOMAIN` source of truth — DEFERRED.** Energy balance is being
> implemented separately; MPCO_Ladruno will *consume* it. When that lands, confirm
> whether `EnergyBalanceRecorder` (classTag 26) exposes its per-region computation in
> a form the domain source can call directly, or whether we lift it. The `ON_DOMAIN`
> channel and v3b reduced-envelopes wait on this.

> [!question]
> **Element engine port fidelity.** The 2595–4279 region holds load-bearing
> workarounds (SSPbrick auto-GaussPoint insert, ForceBeamColumn3d unclosed-tag
> auto-close, shell `material` vs `section` keyword swap). Port as a structural move
> with the parity harness as the regression net — do **not** "clean up" these.

> [!question]
> **Parallel (`sendSelf`/`recvSelf`).** The new `ON_DOMAIN` / `ENVELOPES` requests
> must be added to the `Serializer` payload. Keep enum numeric values stable once
> assigned. Envelope state is per-rank; decide whether p0 reduces or each part writes
> its own envelope (lean: per-part, reduce in the reader).

- **ODR / linkage (D4):** resolved by `namespace mpcol`. Verify no leaked global
  symbols when both recorders are linked (`nm`/`dumpbin` spot-check post-build).
- **Double HDF5 load:** two `LibraryLoader` singletons → benign second `dlopen`;
  confirm no double-free on shutdown.
- **classTag allocation:** pick a free `RECORDER_TAGS_*` value; record it so it
  doesn't collide with other Ladruno additions (energy=26, ExplicitBatheLNVD=63).
- **CMake:** new sources under the existing `_HDF5` guard; mirror MPCORecorder's
  entry. Cross-check against [[project_opensees_cmake_gotchas]].

## Implementation log

*(filled in during execution; move to `Ladruno_internal/implemented_mpco_ladruno.md`
when v1 lands)*

- 2026-05-29 — ADR drafted. Monolith parsed (5-layer map). Four decisions settled:
  sibling fork, apeGmsh-native schema, refactor-to-parity v1, `namespace mpcol`.
- 2026-05-30 — Parallel + envelope design settled (D5–D7): keep per-partition files
  (verified record-loop comm = one int/step, line 4594); envelopes via `EnvelopeSink`
  wrapping any source, phased zero-MPI-consistent (v3a) → reduced (v3b, co-lands with
  energy); componentwise absmax + `ARG_STEP`, per-stage reset, periodic rewrite.
  Energy/`ON_DOMAIN` deferred (built elsewhere; we consume).
