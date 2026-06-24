# apeGmsh ↔ Ladruno fork — feature scoping report

**Date:** 2026-06-24
**Method:** multi-agent scoping sweep (11 agents: 3 inventory → 7 per-domain → 1 synthesis).
Inventory drawn from `LEDGER_implementations.md`, `banner_features.txt`, the
`Ladruno_implementation/*` ADRs/guides, `ladruno_apegmsh_contract.md`, and the
`apegmsh-helper` skill references (`api-cheatsheet.md`, `ladruno.md`).
**Scope:** which fork-only features apeGmsh should expose, the natural Gmsh-driven use
cases, and the gap/effort/priority to wire them. Companion to
[`ladruno_apegmsh_contract.md`](ladruno_apegmsh_contract.md) (the emit/read handshake).

> This is a forward-looking plan, not a commitment. Effort is S/M/L (small/medium/large).
> Priority is P0 (foundational) → P4 (fork-side defer). "Fork-side defer" = do **not**
> build an apeGmsh generator until the fork ships the underlying feature.

> **Correction (2026-06-24) — read apeGmsh `main` source directly.** This report was
> built from the (stale) apegmsh-helper skill + fork contract. Verified against apeGmsh
> v2.0.0 source: the **inverse-map machinery** and **`g.reinforce`** are **already
> shipped** — so backlog items **#1–#3 are largely done**, not P0 build work. The
> remaining *embedded* gap is **`g.embed`** (`LadrunoEmbeddedNode`) only. Also: RBE3
> **`weighting="area"` tributary weights ARE computed apeGmsh-side** (the resolver),
> contradicting the "#6 silent `-area` cliff" framing here — the cliff is only in the
> *fork's* `-area` parser, which apeGmsh correctly sidesteps by emitting explicit `-w`.
> The pipeline is now **online**: apeGmsh's live runner can target the fork build via a
> backend-selection seam, with an end-to-end integration harness (see
> [[ladruno_apegmsh_contract]] → "Backend selection & integration harness").
>
> **Update (2026-06-24) — the P0 generators are now SHIPPED.** `g.embed`
> (`LadrunoEmbeddedNode`, backlog #2) shipped in [apeGmsh#721], and the **`contactSurface`
> + `contact` pair** (backlog #4/#5) shipped in [apeGmsh#722] — both stacked on the
> backend/harness PR [apeGmsh#720], each covered by a `ladruno_fork` live test. Contact
> v1 is **core-first**: NTS penalty + mortar frictionless/friction + auto `-outward`;
> `-soft`/`-visc`/`-consistanttan`/`-geomtan` + self-contact remain the deferred
> follow-up. So the three genuinely-missing P0 generators (#2/#4/#5) are all done; what
> remains is the deferred contact extensions + the lower-priority bulk/robust/cohesive
> backlog below.

---

## 1. Executive summary

apeGmsh↔fork integration is **uneven but well-mapped, and the gaps are almost entirely
emit-side**. Two domains are effectively done; four are partial-to-absent and they all
converge on a single fact:

> **apeGmsh owns the geometry the fork elements deliberately do not carry** — named faces,
> point→host-xi inverse maps, mesh size fields, tributary areas/lengths, CAD face normals.
> So the *generator* is the natural — often the *only correct* — home for these primitives.

The single highest-leverage missing piece is the **straight-sided inverse-map machinery**
(global point → host element xi), which unlocks `g.reinforce` + `g.embed` across the
embedded, cohesive/bond-slip, and bulk-constitutive RC domains at once.

Three headline opportunities:
1. The **`contactSurface` + `contact`** generator pair — named-face-pair contact is the
   fork's pounding/impact/collapse north star and has **zero** apeGmsh surface today.
2. The shared **inverse-map → `g.embed` → `g.reinforce`** ladder — RC reinforcement,
   bond-slip pull-out, non-conformal submodeling.
3. Closing the **correctness cliffs already latent in "shipped" paths** (see §8).

---

## 2. Domain status

| Domain | Status | Headline |
|---|---|---|
| Coupling (RBE2 33012 / RBE3 33011) | ✅ EXPOSED (mature) | Typed emitters + readers ship. Only real gap: RBE3 `-area` silent equal-weights cliff. |
| Output & recorders (.ladruno / Monitor / EnergyBalance) | ✅ EXPOSED (most mature) | Read-side co-designed with apeGmsh. Remaining work is doc/contract hygiene + additive read paths. |
| Bulk constitutive (plasticity/concrete) | 🟡 PARTIAL | Materials reach the deck via stringly-typed passthrough; missing the constitutive-*aware* wiring. |
| Embedded elements (rebar / node ties) | ✅ SHIPPED | `g.reinforce` (LadrunoEmbeddedRebar) + `g.embed` (LadrunoEmbeddedNode, [apeGmsh#721]) both ship; live-tested on the fork. `g.constraints.embedded` (vanilla ASDEmbeddedNode) kept as fallback. |
| Robust solve & integrators | 🟡 PARTIAL | Explicit lane well-served (auto-dt). Static-robust lane absent (no ArcLength/IndirectControl/DR emitter). |
| Contact & progressive collapse | ✅ SHIPPED (core) | `g.constraints.contact(formulation='nts'\|'mortar')` emits `contactSurface`+`contact`+`LadrunoContact` ([apeGmsh#722]); live-tested. Deferred: `-soft`/`-visc`/`-consistanttan`/`-geomtan` + self-contact. Element-removal stays upstream. |

---

## 3. Top cross-domain use cases (ranked)

1. **RC member: discrete rebar embedded in a non-conforming concrete solid** (`g.reinforce`)
   — bar EDGE PG threaded through a meshed concrete VOLUME; apeGmsh inverse-maps each bar
   point into its host element and emits `LadrunoEmbeddedRebar -host -xi -dir`, owning the
   bar axis, `bondScale = π·d_b·L_trib`, and host weights. The v2 RC-3D deliverable and the
   convergence point of three domains.
2. **Face-pair pounding / impact / progressive-collapse contact** (NTS + SOFT explicit) —
   two named surface PGs → `contactSurface(master)` + `contactSurface(slave/-slave-segments)`
   + `contact(-kn/-mu | -soft SOFSCL)` under `CentralDifferenceLadruno`; master CAD-face
   normal auto-fills `-outward`. The meshed FACE *is* `LadrunoContactSurface`.
3. **Beam/shell RP into a solid face — moment transfer (RBE3) and rigid platen (RBE2)** —
   a 6-DOF point PG into a meshed solid FACE. *Already exposed and mature* — a reference
   point, not a build.
4. **Non-matching mortar/ALM interface between two independently-meshed parts** — `g.parts`
   assembly, shared interface as master face + `-slave-segments` faceted face (no boolean
   weld), `contact(-mortar -epsN -augTol -ngp)`. The accuracy lane.
5. **Crack-band-objective concrete with `lch` pulled from the mesh** — the meshed element
   size *is* the crack-band length `LadrunoConcrete3D`/`ASDConcrete3D` need; apeGmsh owns
   the size field → emits `-autoRegularization`/per-region `-lch` and auto-suggests the
   unsymmetric solver.
6. **Non-matching solid stitch / submodel + staged stress-free birth** (`g.embed`) — refined
   Part embedded in a coarse host volume; `LadrunoEmbeddedNode -host -xi -k auto` with `g0`
   birth. Drop-in upgrade over the unconditioned 1e18-penalty `ASDEmbeddedNodeElement` path.
7. **Robust softening/buckling static solve with a Gmsh-named control node** (`robust_drive`)
   — maps a named control point → (node, dof, du) and hands off to
   `Ladruno_scripts/robust_drive.py`.

---

## 4. Prioritized backlog

| # | Item | Domain | Kind | Effort | Pri |
|---|------|--------|------|--------|-----|
| 1 | Shared straight-sided inverse-map machinery (global point → host xi, guarded Newton, ADR-20 D3). Foundation for `g.embed`/`g.reinforce`. | Embedded (shared) | api-wiring | M | **P0** |
| 2 | `g.embed(nodes, host, k='auto', enforce, explicit, staged)` → `LadrunoEmbeddedNode -host -xi -k auto`; default U-only + `g0` on; gate UR/UP/D9/`-corot` off. | Embedded | api-wiring | M | **P0** |
| 3 | `g.reinforce(host, bars, bond, enforce, explicit)` → `LadrunoEmbeddedRebar -host -xi -dir (-perfect | -bond LadrunoBondSlip -bondScale=π·d_b·L_trib) -kt auto`. The headline RC generator. | Embedded / Bulk / Cohesive | api-wiring | L | **P0** |
| 4 | `contactSurface` generator: walk a named surface PG's boundary facets into `LadrunoContactSurface` flat-connectivity/stride (`-master`/`-slave`/`-slave-segments`). | Contact | api-wiring | M | **P0** |
| 5 | `contact()` verb: typed emitter for the verified `-kn/-mortar/-epsN/-mu/-epsT/-cohesion/-tauMax/-augTol/-maxAug/-ngp/-tie/-visc/-soft/-outward/-consistentNormal` grammar + auto `constraints('LadrunoContact')`; auto-fill `-outward` from master CAD-face normal. | Contact | api-wiring | M | **P0** |
| 6 | Make RBE3 `weighting='area'` actually compute per-node tributary areas apeGmsh-side and emit explicit `-w`. Fork `-area` is a silent stub (equal weights, no error). | Coupling | api-wiring | M | **P0** |
| 7 | Wire crack-band `lch` from the mesh: emit `-autoRegularization`/per-region `-lch` for `LadrunoConcrete3D`/`ASDConcrete3D`. Confirm `LadrunoBrick getCharacteristicLength` override is lifted. | Bulk | api-wiring | M | P1 |
| 8 | Auto-suggest/emit the unsymmetric solver (UmfPack/FullGeneral) when `LadrunoConcrete3D` is present (non-symmetric tangent → silent mis-solve under symmetric solver). | Bulk | modeling-guidance | S | P1 |
| 9 | Disambiguation doc: `g.constraints.tie/mortar/embedded` today emit upstream `ASDEmbeddedNodeElement` (per-node 1e18 tie), NOT fork `LadrunoContact -tie/-mortar` nor `LadrunoEmbeddedNode`. Name the fork generators unambiguously; route hex/tet/quad hosts to the fork path. | Contact / Embedded | doc | S | P1 |
| 10 | Three typed cohesive/bond material emitters (uniaxial `LadrunoBondSlip`, uniaxial `LadrunoCohesiveHinge`, nDMaterial `LadrunoCohesiveHingeBiaxial`). | Cohesive | api-wiring | S | P1 |
| 11 | Add `LadrunoDispBeamColumn2d/3d` to `_ELEM_REGISTRY` (Line2, needs_transf, beamIntegration slot) with typed emitter exposing `-hinge/-hingeY/-hingeBiaxial`. Without it the cohesive-hinge family is unreachable. | Cohesive | api-wiring | M | P1 |
| 12 | Typed `ops.integrator.*` emitters for `LadrunoArcLength` (incl. `-adapt`/`-stabilize`), `LadrunoIndirectControl`, `LadrunoDynamicRelaxation`. Static-robust lane has no emitter today. | Robust-solve | api-wiring | M | P1 |
| 13 | Typed `_MatSpec` registry/schemas for fork nDMaterials (`LadrunoJ2`, `LadrunoConcrete3D`, `LogStrain`, `StagedStrain`/`InitDefGrad`, `ASDConcrete3D`); document wrapper-nesting order rules. | Bulk | api-wiring | M | P2 |
| 14 | `ops.robust_drive(...)` bridge mapping a named Gmsh control point/PG → (node, dof, du); wire `LadrunoStabilizedUnbalance` test + `ladrunoArcLength` query → `robust_drive.py`. | Robust-solve | api-wiring | M | P2 |
| 15 | RBE2/RBE3 modeling docs: auto-derive `-dof` from slave-PG ndf (3-DOF face → `-dof 1 2 3`) + ship the RBE2-vs-RBE3 decision table (rigid vs flexible patch, master/dependent role-flip). | Coupling | doc | S | P2 |
| 16 | 2D-host weight path: `-xi` via `getInterpolationWeights` is 3D-only; for `LadrunoQuad/CST` hosts the generator must compute bilinear/area `-shape` weights itself. | Embedded / Cohesive | api-wiring | M | P2 |
| 17 | Reconcile stale contract rows (BezierTri6 + Monitor consumer marked deferred but shown shipped in `ladruno.md`) and document the contact-augment recorder-skip contract (STEP is monotonic, not 0-based/contiguous). | Output | doc/contract | S | P2 |
| 18 | `ops.analyze_augmented(...)` helper wiring the right query (`ladrunoMortarPenetration`/`ladrunoMortarTieResidual`) for the contact/tie FACES apeGmsh just emitted. | Robust / Contact | api-wiring | M | P2 |
| 19 | Mass-scaling (SMS lumped + consistent/Olovsson) + HRZ lumping flags on the explicit integrator emitters, plus a Gmsh-region selector for which elements to scale. | Robust-solve | api-wiring | M | P2 |
| 20 | `results.energy`/`results.regions` read path for ON_DOMAIN/ON_REGIONS energy (energy_balance stackplot with RES/ERR%); pin `from_ladruno` to FORMAT_VERSION with a two-version window; shell per-layer stress decode. | Output | api-wiring | M | P3 |

---

## 5. Quick wins (do first — all S effort, several are correctness fixes)

- **Unsymmetric-solver auto-suggest** whenever `LadrunoConcrete3D` is present — a symmetric
  solver *silently mis-solves* its non-associated tangent. (#8)
- **tie/mortar/embedded vocabulary disambiguation doc** — today those verbs silently emit
  upstream `ASDEmbeddedNodeElement`, not gap-aware fork contact. A user asking for "mortar
  contact" gets the wrong thing with no error. (#9)
- **Three typed cohesive/bond material emitters** — prerequisite that unblocks two delivery
  paths. (#10)
- **Auto-fill contact `-outward`** from the master CAD-face OCC normal — a correctness win
  the raw command can't offer on curved/closed surfaces. (#5, sub-item)
- **RBE2 `-dof` auto-derive** from slave-PG ndf + the **RBE2-vs-RBE3 decision table**. (#15)
- **Reconcile stale Output contract rows** + the contact-augment STEP-monotonicity note. (#17)
- **Force RP ndf=6** (via `g.decouple_node`/`ops.ndf`) and wire `-k auto -host` auto-resolution
  from the slave surface PG for RBE2/RBE3 — otherwise the run is refused or needs a hand-picked
  numeric `-k`.

---

## 6. Biggest bets

1. **The shared inverse-map → `g.embed` → `g.reinforce` ladder** (M+M+L). One piece of
   machinery (global point → host xi) unlocks RC discrete reinforcement, CEB-FIP bond-slip
   pull-out, non-conformal submodeling, and `LadrunoConcrete3D`-in-cages across three domains.
   Highest-leverage investment in the whole roadmap; fully pre-specced in the contract.
2. **The `contactSurface` + `contact` generator pair** (M+M). Named-face-pair contact is the
   fork's headline pounding/impact/progressive-collapse target with zero apeGmsh surface
   today; also unlocks the mortar accuracy lane, the `-tie` surface bond (strictly stronger
   than the shipped `ASDEmbeddedNode` tie), and the `analyze_augmented` companion.
3. **Constitutive-aware wiring for the bulk concrete lane** (`lch`-from-mesh +
   unsymmetric-solver + typed `_MatSpec` + fc/ft/Gf calibration helper). Removes the two
   most error-prone manual steps in RC-solid modeling.
4. **The `LadrunoDispBeamColumn2d/3d` registry entry + the cohesive-hinge family.** Unlocks
   mesh-objective fracture-energy plastic hinges (Gf per hinge, biaxial Mz-My) for frame
   pushover/collapse — a whole capability currently unreachable from apeGmsh.

---

## 7. Per-domain detail

### 7.1 Contact & progressive collapse — 🔴 ABSENT
Fork ships the whole `LadrunoContact` subsystem (NTS penalty, mortar/ALM, Coulomb/Tresca
friction, viscous stabilization, SOFT=1/2 explicit penalty, mesh-tying, rigid plane,
within-step augmentation) exposed only via raw Tcl/Py. apeGmsh has **no** surface for any of
it. **Build the `contactSurface` + `contact` generator pair** keyed on "a contact interaction
is exactly two named meshed faces." Auto-fill `-outward` from the master CAD-face normal.
Disambiguate the existing `g.constraints.tie/mortar` (silently emits upstream
`ASDEmbeddedNodeElement`). **Defer:** element-removal collapse is *upstream OpenSees, not a
fork feature* — out of scope. Self-contact is fork-deferred (ADR-47 #3).

### 7.2 Coupling constraints (RBE2/RBE3) — ✅ EXPOSED (mature)
Both `LadrunoKinematicCoupling` (33012) and `LadrunoDistributingCoupling` (33011) ship typed
emitters (`g.constraints.kinematic_coupling`/`distributing_coupling`) + readers; the
meshed-FACE-as-named-handle workflow gives every RP-to-surface coupling a clean entry point.
Remaining work is three high-value low-effort items already pre-specced in the `*_guide.md`
§10.6/§10.7: (1) make `weighting='area'` actually compute tributary areas (the fork `-area`
is a **deliberate stub** → silent equal weights, a real correctness cliff); (2) auto-derive
`-dof` from slave-PG ndf; (3) the RBE2-vs-RBE3 decision doc (the #1 documented wiring
mistake). Medium follow-ups: `-bipenalty -dtcr` integrator-aware emission, `-k auto -host`
auto-resolution, forcing RP ndf=6. **Do not** invent finite-rotation support (fork-side
defer) — just document the small-relative-rotation validity boundary.

### 7.3 Embedded elements — 🟡 PARTIAL
`LadrunoEmbeddedRebar` (33005) + `LadrunoEmbeddedNode` (33006) + `LadrunoBondSlip` (33002)
ship with full grammars, a shared kernel, and per-element generator contracts already in the
guides + contract §10.6. apeGmsh emits only vanilla `ASDEmbeddedNodeElement`. The enabling
fact: apeGmsh owns the geometry (bar points, host volumes, meshed edge lengths, tributary
lengths) and therefore the inverse map the elements deliberately don't carry. Build order:
(1) shared straight-sided inverse map once; (2) `g.embed` (simpler, isotropic — immediately
upgrades the `ASDEmbeddedNodeElement` path with conditioned `-k auto` + explicit + staged
`g0`); (3) `g.reinforce` (adds `-dir`, `bondScale`, bond law); (4) the 2D-host `-shape`
path. **Defer:** curved hosts, dowel action, experimental UR/UP/D9/`-corot` modes.

### 7.4 Cohesive & bond-slip interface materials — 🔴 ABSENT
None of the three materials (`LadrunoCohesiveHinge` 33003, `LadrunoCohesiveHingeBiaxial`
33004, `LadrunoBondSlip` 33002) nor their consuming elements have an apeGmsh emitter.
Delivered through two shipped fork mechanisms: (a) bond-slip as the axial slot of an embedded
bar-in-solid (needs `g.reinforce`); (b) cohesive hinges as embedded strong-discontinuity
moment-rotation laws *inside* frame elements (needs `LadrunoDispBeamColumn2d/3d` in the
registry). Prioritize: g.reinforce (L), DispBeamColumn registry entry (M), the three typed
material emitters (S). **Modeling caveat:** the cohesive hinge is NOT a surface/zeroLength
interface — it lives inside the frame element; document the delivery mechanism + the
softening solution-control caveat (DisplacementControl/ArcLength/IMPL-EX, never load control).
**Defer:** the surface cohesive-zone element (ADR-56) is a fork-side sketch — the cleanest
Gmsh fit (masonry joints, RC cold joints, debonding) but blocked on a fork element that
doesn't exist yet.

### 7.5 Bulk constitutive (plasticity/concrete) — 🟡 PARTIAL
Every material ships fork-side and reaches the deck via generic `ops.nDMaterial.<Type>` string
passthrough; apeGmsh already meshes the solids and selects faces/volumes/Parts. Basic
confined-column / large-strain / staged-construction use cases work today. Missing the
*constitutive-aware* wiring: (1) crack-band `lch` from the mesh (high value, most error-prone
manual step); (2) unsymmetric-solver auto-suggest + softening guidance for `LadrunoConcrete3D`
(S, prevents silent mis-solves); (3) typed `_MatSpec` schemas (discoverability); (4)
wrapper-nesting order docs (StagedStrain outermost, isotropic under LogStrain,
`LadrunoJ2Finite` for combined+large-rotation). The rebar-cage `g.reinforce` generator and
the `.ladruno` RC damage post-processor are the next tier (overlap embedded/output). **Defer:**
confined-fiber 1D view, `-eta`+`-implex`, finite prestrain.

### 7.6 Robust solve & integrators — 🟡 PARTIAL
Split sharply along analysis-vs-geometry. The **explicit lane is well-served**: typed
`ExplicitBathe`/`ExplicitBatheLNVD`/`CentralDifferenceLadruno` emitters + the high-value
`ops.critical_time_step()` + `ops.analyze_explicit(duration=, safety=)` auto-dt driver, with
mesh-derived prerequisites already wired. The **static-robust lane is absent**:
`LadrunoArcLength`/`IndirectControl`/`DynamicRelaxation` have no emitter; `robust_drive` and
`analyze_augmented` have no entry point. Prioritize: typed emitters for the static
path-followers + DR (reuse the explicit grammar); an `ops.robust_drive(...)` bridge whose one
real geometric anchor is mapping a named control point/PG → (node, dof, du);
`ops.analyze_augmented(...)` tied to the contact/mortar FACES apeGmsh emits. Add mass-scaling
(SMS/HRZ) flags + a region selector. **Defer:** the half-increment-residual gate (ADR-52
W1-I1b) has no geometric anchor → doc/OpenSeesPy guidance only; track the ExplicitBathe* 6→1
flag collapse until the fork refactor lands.

### 7.7 Output & recorders — ✅ EXPOSED (most mature)
The `.ladruno` recorder was co-designed with apeGmsh's reader; `Results.from_ladruno` + the
Monitor consumer are shipped. Remaining work is mostly doc/contract hygiene: (1) reconcile
the stale contract rows so read-side status stops being double-sourced; (2) document the
contact-augment recorder-skip contract (`ladrunoBeginAugment`/`EndAugment` →
`Domain::contactAugmenting` keeps the time axis clean during held-load Uzawa sweeps; a reader
assuming 0-based/contiguous STEP, or a helper that forgets to wrap the held-load `analyze`,
corrupts results). Medium additive read-side features to confirm-or-finish: FORMAT_VERSION
pinning + chunked/legacy handling, shell per-layer decode (`gp = GAUSS_ID // nLayer`), the
ON_DOMAIN/ON_REGIONS energy ResultLevel + `energy_balance` plot, multi-stage read coverage.
**Defer:** envelopes (ENVELOPES v3), Monitor element/region/parallel channels.

---

## 8. Correctness cliffs surfaced (latent bugs, not just gaps)

These already bite users on supposedly-working paths — worth fixing regardless of the broader roadmap:

1. **RBE3 `weighting='area'` silently degrades to equal weights.** The fork `-area` parser is
   a stub (`OPS_LadrunoDistributingCoupling.cpp` warns and falls back). If apeGmsh forwards
   `-area` instead of computing `-w`, the user gets equal weights with no error. → #6.
2. **`g.constraints.tie/mortar/embedded` silently emit the wrong element.** They produce
   upstream `ASDEmbeddedNodeElement` (per-node 1e18 stiffness tie), not the gap-aware fork
   `LadrunoContact -tie/-mortar` or the conditioned `LadrunoEmbeddedNode`. A "mortar contact"
   request gets a rigid embed tie with no gap/KKT/friction and no warning. → #9.
3. **`LadrunoConcrete3D` under a symmetric solver mis-solves.** Its non-associated tangent is
   non-symmetric; nothing currently forces UmfPack/FullGeneral. → #8.
4. **Contact-augment time-axis corruption.** Any apeGmsh augmentation helper that forgets to
   wrap the held-load `analyze` in `ladrunoBeginAugment`/`EndAugment` writes spurious
   frames/steps into `.ladruno`. → #17.

---

## 9. Fork-side defers — do NOT build apeGmsh generators against these yet

Surface cohesive-zone element (ADR-56 sketch); self-contact (ADR-47 #3); curved-host inverse
maps; quantitative dowel action; UR/UP/D9/`-corot` embedded modes; finite-rotation geometric
stiffness for RBE2/RBE3; confined-fiber 1D `LadrunoConcrete3D` view + `-eta`+`-implex`;
half-increment-residual gate (no geometric anchor → doc only); ExplicitBathe* 6→1 flag
collapse (track until refactor lands). **Upstream-not-fork:** progressive-collapse element
removal (`recorder Collapse` / Talaat-Mosalam / Elwood).
