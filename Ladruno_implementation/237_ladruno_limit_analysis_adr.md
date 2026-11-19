---
title: "ADR 237 — Certified finite-element limit analysis in the Ladruno fork (WP-237)"
project: Ladruno
type: ADR
status: "Proposed (2026-11-18) — design only, no code"
priority: high
owner: fork maintainer
wp: "WP-237 (provisional — confirmed when the limit-analysis program issue claims it)"
related:
  - "[[_wp237_orchestration_plan]]: work packages, roles, gates"
  - "[[87_ladruno_depth_with_width_adr]]: warrant package, merge rules"
  - "[[232_ladruno_mpm_general_driver_adr]]: the engine-with-its-own-file precedent"
  - "Fork ADR 236 (draft PR #1117): material history transfer; not a dependency of Level A"
  - "apeGmsh ADR 0136 (draft PR apeGmsh#3052): adaptive mesh convergence; needed from P4"
tags: [adr, limit-analysis, conic-optimization, bounds, certificate, apegmsh, wp-237]
updated: 2026-11-18
---

# ADR 237 — Certified finite-element limit analysis

> **Status.** Proposed. Nothing here is built, and no row below is accepted until the owner accepts
> this ADR. Each decision carries one of two marks: **[owner]** was decided by the owner on
> 2026-11-18; **[proposed]** is a working default or a recommendation of this document.
>
> **Number.** 237 is the next number above every WP found on 2026-11-18: no PR title contains
> `WP-237`, and the highest `wp/*` branch on the remote is `wp/236-state-transfer-adaptivity`. It
> stays provisional until the program issue claims it (plan §8). That issue does not exist yet.
>
> **Checked against.** Fork claims: `cfe478254` (`origin/ladruno` on 2026-11-18), cited as
> `path::symbol`. Fork ADR 236: PR #1117 head `f7f89a4c6`. apeGmsh ADR 0136: PR apeGmsh#3052 head
> `668b99a44`. apeGmsh ADRs 0121 and 0134: `main` at `60d0517bc`. Only apeGmsh ADR text was read, not
> its source.
>
> **Discipline.** Under ADR-87 an unmeasured claim is a claim. §19.3 lists everything this document
> could not verify. §10 lists where checking the working design against mathematics or code changed
> it; those items are proposed amendments, not settled changes.

---

## 0. Decision summary

| Question | Answer | Mark |
|---|---|---|
| What is built | 3D finite-element Lower and Upper bounds for perfectly plastic, associated solids, with an independent certificate | [owner] 3D; [proposed] the rest |
| Where | In the fork, as a one-shot engine behind one `limitAnalysis` command row | [owner] fork; [proposed] the shape |
| Input | A separate versioned problem file `<stem>.la.h5`, not the Domain | [proposed] O3 |
| Output | The engine's own `.ladruno` writer, not `recorder ladruno` | [proposed] |
| Materials | A `LimitMaterial` registry, separate from `NDMaterial`. Von Mises and Drucker-Prager first; Tresca and Mohr-Coulomb in P5 | [proposed] O4 |
| Lower | Discontinuous linear stress on straight tetrahedra (§8) | [proposed] |
| Upper | P2 velocity on straight tetrahedra (§9), **with a locking problem in true 3D that this ADR measured and has not resolved (§10.1)** | [proposed], open |
| Claim language | Certified bound, bracket, estimate, inadmissible (§11) | [proposed] |
| Domain | Every result states a `domain_claim`; box-labelled bounds first (§12) | [proposed] O5 |
| Solvers | Open source only; no solver linked by default; no component that needs a licence request (§14) | [owner] the rule; [proposed] the mechanism |
| FEniCSx | An external test reference, outside the fork (§15) | [owner] |
| Remeshing | Level A through apeGmsh ADR 0136, amended (§17.1) | [owner] central; [proposed] the amendments |
| Restart remeshing | Not a dependency of Level A; coordinated under ADR 236 (§17.3) | [proposed] |

---

## 1. Problem and scope

### 1.1 Problem

An elastoplastic analysis gives a load path, but says little about how far its peak is from the true
collapse load of the idealised problem. Limit analysis answers that question directly. For a perfectly
plastic solid with associated flow, a statically admissible stress field gives a lower bound on the
collapse multiplier, and a kinematically admissible velocity field gives an upper bound
[DPG52]. A bracket from both sides is an error statement that does not rely on a reference solution.

The fork has no such capability. A search of `SRC/`, `tests/`, `ci/` and `Ladruno_scripts/` at
`cfe478254` finds no occurrence of `limitAnalysis` or `LadrunoLimit`. Files with `limit` in their
names are unrelated (`SRC/material/uniaxial/limitState/`, the `*HSIncrLimit` integrators, and
`Ladruno_implementation/81_quadratic_hex_limit_load_measurement.md`).

### 1.2 Goal

A general oracle for fork users: certified 3D bounds, a query workflow for combined loading of a rigid
footing, an apeGmsh front and back end, and adaptive meshing that makes the bounds tight enough to use.

### 1.3 Non-goals for version 1

- Non-associated flow, softening, hardening, state-dependent strength, and history. The bound theorems
  do not hold for non-associated materials in general [DD93]. A result here is a statement about the
  associated idealisation, never about a calibrated soil model.
- Replacing elastoplastic analysis. The engine is a check on such models.
- Large deformation, seepage coupling, or dynamics.
- Reading mesh, loads, supports or footing from the Domain.
- Curved elements. The rigour arguments in §8 and §9 need straight tetrahedra.
- OpenSeesSP and OpenSeesMP. Version 1 refuses them, as `mpm` does.
- Commercial solvers or commercial reference codes, as dependencies or as gates.

---

## 2. Formulation in brief

Let `K(x)` be the convex set of admissible stresses, with support function
`d(x, D) = sup { s : D, s in K(x) }`, where `D` is a strain rate. Stress is tension positive. Loads are `t0 + lambda t1` on
surfaces and `b0` in the body. `lambda` is the multiplier.

- **Lower.** If `sigma` satisfies equilibrium with `lambda`, traction continuity, traction boundary
  conditions and `sigma(x) in K(x)` everywhere, then `lambda <= lambda*`.
- **Upper.** If `v` satisfies the velocity boundary conditions, has `W1(v) > 0`, and has finite
  dissipation, then `lambda* <= (Diss(v) - W0(v)) / W1(v)`, where `Diss` is total dissipation and `W0`, `W1`
  are the work rates of the fixed and scaled loads.

Both problems are convex. With the cones of §7 they are conic programs
[Mak06, Mak07, Kra07, BH21]. Earlier nonlinear-programming and review treatments are [LS02a, LS02b, Slo13].

The two criteria of the first phase, with `s` and `e` the deviators of stress and strain rate, `norm` the
Frobenius norm and `I1 = tr(sigma)`:

| Criterion | Yield | Admissible strain rate | Dissipation |
|---|---|---|---|
| Von Mises, `k` = shear strength | `norm(s) / sqrt(2) <= k` | `tr(D) = 0` | `sqrt(2) k norm(e)` |
| Drucker-Prager | `norm(s) / sqrt(2) + alpha I1 <= k`, `alpha > 0` | `tr(D) >= 3 sqrt(2) alpha norm(e)` | `k tr(D) / (3 alpha)` |

Outside the admissible cone the dissipation is `+infinity`. §10.4 depends on this. Both rows follow from
the definition of `d`.

---

## 3. The fork: what exists, with its limits

Every row was read at `cfe478254`.

| Asset | Use here | Its limit |
|---|---|---|
| `SRC/interpreter/LadrunoCommandTable.h::LADRUNO_COMMAND`, row `"mpm"` | one row registers the command in all engines | `ci/check_ladruno_commands.py` fails a command registered elsewhere (H1). `tests/test_ladruno_command_registry.py` holds both shipped engines to the table |
| `SRC/domain/pattern/ladrunoMPM/LadrunoMPMCore.h` | precedent for a pure kernel: its header says "NO OpenSees dependency", and it is unit-tested with `g++` alone | it depends on the fork's hypo-kernel header; the limit kernel would not need that |
| `SRC/domain/pattern/ladrunoMPM/LadrunoMPMIO.h::readParticleSet`, `::ParticleWriter` | precedent for a contract-checked sidecar with rule codes, and for an engine-owned `.ladruno` writer | both need `_HDF5` and refuse without it |
| `tests/mpm_interface.md`, `tests/mpm_cases.txt`, `tests/mpm_reference.py`, `tests/mpm_contracts.py` | precedent for a frozen interface, cases, an oracle and an executable contract | the C++ author must not read the oracle |
| `SRC/utility/LadrunoThreads.h::ladruno_parallelBuild`, used in `SRC/domain/pattern/ladrunoMPM/OPS_LadrunoMPM.cpp` | refusal in OpenSeesSP and OpenSeesMP | — |
| `SRC/api/packages.cpp::getLibraryFunction` | existing run-time loader for a library and symbol, on Windows and POSIX | it appends the platform suffix and was written for element and material packages; it is a precedent, not yet an adapter loader |
| `CMakeLists.txt`, `option(LADRUNO_CMS ...)` | precedent for an opt-in library; its comment requires that OFF leaves stock and serial targets unchanged at source-list level | it is OFF by default |
| `Ladruno_implementation/ladruno_schema_v1.md` §2a, §3, §7.5, §9 | `.ladruno` layout, element self-description, the additive `ON_PARTICLES` precedent, version rules | §7.3 still labels `ON_DOMAIN` as deferred, although `SRC/recorder/LadrunoRecorder.cpp` creates the group; the text needs reconciling before this engine relies on it |
| `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py::validate` | executable schema | see W1 |
| `SRC/system_of_eqn/linearSOE/mumps/` | an open sparse solver already in the fork's build, for a possible tailored method | only under the conditions of §14.4 |

| # | Wall | Consequence |
|---|---|---|
| W1 | `ladruno_format.py::validate` rejects any stage `KIND` other than `transient`, `static` or `eigen`. apeGmsh ADR 0134 records that its TypeScript reader throws a whole-file layout error on an unknown `KIND` | a `lower` or `upper` stage kind would break both. §6 keeps `KIND = "static"` |
| W2 | The engine has no Domain host, so `Domain::clearAll` will not reach it. All three wipe paths (`SRC/tcl/commands.cpp::wipeModel`, `SRC/interpreter/OpenSeesCommands.cpp::OpenSeesCommands::wipe`, `SRC/interpreter/PythonAnalysisBuilder.cpp::ops_wipeModel`) call `OPS_clearAllNDMaterial` | `mpm` is reached through its `LoadPattern` host; this engine needs an explicit hook (§4.3) |
| W3 | `ci/check_manifest.py::main` is driven by new class tags. A feature with no class tag is not forced to have a manifest row | ADR-87 D1 would pass silently here unless the row is added deliberately (§16) |
| W4 | Zone-A installs `pytest pyyaml numpy sympy scipy h5py conan` (`.github/workflows/ladruno.yml`) and no conic solver | operator and checker tests can run in Zone-A; solver-dependent ones cannot without changing that line (§16) |
| W5 | The repository's `COPYRIGHT` is the University of California text, limiting use and distribution by purpose | any bundling of third-party code needs a licence review; see O10. This is an engineer's reading, not legal advice |

---

## 4. Architecture

### 4.1 Layers

```text
apeGmsh or any writer
  geometry, labels, straight tet10 mesh, adaptive driver
      |  <stem>.la.h5
      v
Fork: limitAnalysis
  LadrunoLimitKernel.h      operators; plain arrays; no OpenSees types
  LadrunoLimitAssembler     builds a solver-neutral ConicProgram
  ConicSolver adapter       run-time, versioned C ABI; no solver linked by default
  writer                    .ladruno stages for Lower and Upper; program export
      |  <stem>.ladruno
      v
Independent Python checker
  rebuilds operators from the problem file, recomputes the result, issues the certificate
      |
      v
apeGmsh Results and viewer; adaptive indicators; comparisons with OpenSees models
```

**D1 [proposed].** The engine follows the shape of the MPM engine: a pure kernel, an assembler, its own
sidecar and writer, one command row. It does not extend `Element`, `NDMaterial` or an integrator. Limit
analysis has no load path, no tangent and no commit.

**D2 [proposed].** The engine does not certify itself. Only the independent checker issues a certificate
(§11).

### 4.2 Data model

```text
LimitProblem     mesh; regions (element -> material, constant unit weight);
                 labelled boundary facets; footings; dead loads
Discretisation   once per mesh: Lower and Upper operators
Query            many per mesh: support(d) | ray(F0, F1)
ConicProgram     sparse A, b, c; cones {zero, nonnegative, soc(n), psd(n)}; map back to fields
Solution         fields, multipliers, status, timings, inputs for the certificate
```

### 4.3 Command surface [proposed]

One `LADRUNO_COMMAND("limitAnalysis", ...)` row, with subcommands in the style of `mpm`: `problem -file`,
`material`, `solver`, `query`, `solve`, `output`, `export`, `info`, `reset`.

- `export` always writes the exact conic program, so independent solvers can run identical data.
- The feature sits behind a CMake option. A disabled build leaves existing targets unchanged.
- `wipe` must clear any process-wide state. The proposed route is a local-extern reset called from
  `SRC/material/nD/NDMaterial.cpp::OPS_clearAllNDMaterial`, the hook used by WP-104 and already reached
  by all three wipe paths. That is one upstream edit with a `vanilla/` ledger row. If the engine holds
  no process-wide state, no edit is needed; 237f decides.
- Classic Tcl is either wired or recorded as absent (O9).

---

## 5. Problem file [proposed, O3]

`<stem>.la.h5` is a versioned contract owned by the fork. 237a freezes it field by field, with rule
codes, as `Ladruno_implementation/232a_mpm_contracts.md` did for MPM.

It must carry: format identity and version; units; physical node coordinates; tet10 connectivity with a
stated node order; element regions; material records; labelled boundary facets with a facet kind;
footing reference points, footprints and interface laws; dead loads; and a content hash.

Why not the Domain, as the first route:

- A Lower bound needs tractions that are exact pointwise on each facet. Consistent nodal loads do not
  carry that information.
- An Upper bound needs supports and rigid footing motion imposed exactly. Penalty coupling elements do
  not do that.
- A strength record is not an `NDMaterial` state.

The fork builds its facet table from connectivity and checks straightness itself. Physical coordinates
cross the file boundary because parametric conventions differ between codes.

A later `-fromDomain` adapter may help Tcl users. Its output is an estimate unless it satisfies every
condition of §8, §9 and §11.

---

## 6. Output [proposed]

The engine writes `.ladruno` itself, following `ParticleWriter` and the schema process.

| Data | Location |
|---|---|
| Mesh | `MODEL/NODES`, `MODEL/ELEMENTS` |
| Upper velocity | `RESULTS/ON_NODES/velocity` |
| Lower stress at tet vertices | `RESULTS/ON_ELEMENTS`, vertices as points |
| Upper strain rate and dissipation | `RESULTS/ON_ELEMENTS` |
| Lower multipliers and recovered velocity | `RESULTS/ON_ELEMENTS` |
| Multiplier, resultants, status, timings | `RESULTS/ON_DOMAIN/limitAnalysis` |
| Raw solution, program hash, cone list | `MODEL/LIMIT_ANALYSIS` |
| Provenance | `INFO` attributes |

**D3 [proposed].** Lower and Upper are separate `MODEL_STAGE` groups with `KIND = "static"` (W1) and a new
attribute naming the bound. The additions are optional groups and attributes, so the schema's §9 rules
apply; the version bump is fixed by 237f.

Three items stay open until tested on fixtures, and this ADR does not claim them: which element group
name and class tag the writer uses; how today's readers treat `MODEL/LIMIT_ANALYSIS`; and how `ON_DOMAIN`
is specified, given the note in §3.

Certificate results are written on the consumer side, not by the engine.

---

## 7. Limit-material registry [proposed]

**D4.** A separate `LimitMaterial` registry. No new virtuals on `NDMaterial`, and no mix-in. A strength
criterion supplies a cone representation and a support function. It has no strain, state or tangent.

| Criterion | 3D cone | Phase |
|---|---|---|
| Von Mises | second-order | P0 |
| Drucker-Prager | second-order | P0 |
| Tresca | semidefinite | P5 |
| Mohr-Coulomb | semidefinite, two 3x3 blocks in the support form [MM08, Kra08] | P5 |

In 3D, Drucker-Prager is not Mohr-Coulomb. A result names its criterion, and a matched parameter set is
valid only for the state it was matched to.

A `-fromND` route may come later as an explicit table that refuses unknown classes. It never infers a
strength from a state-dependent model.

---

## 8. Rigorous Lower bound [proposed]

Field: stress linear in each straight tetrahedron, discontinuous between elements (24 unknowns per
element) [BH21 §4.2; Sloan88 and Mak06 for the lineage].

| Condition | Imposed as | Why it is exact |
|---|---|---|
| Equilibrium | `div sigma + b = 0` once per element, `b` constant per element | linear stress has constant divergence |
| Interior facet | `sigma n` continuous at the facet's three vertices | traction is linear on a flat facet |
| Traction facet | `sigma n = t0 + lambda t1` at three vertices | exact only if `t0`, `t1` are linear on the facet; the contract requires that |
| Free facet | `sigma n = 0` at three vertices | same |
| Symmetry or roller | zero tangential traction at three vertices | same |
| Fixed facet | no traction condition | a support may supply any traction |
| Yield | `sigma in K` at four vertices | `K` is convex and the element is the convex hull of its vertices, with constant strength per element |
| Footing resultant | closed-form integral of the linear interface traction | force and moment integrands are polynomial on a flat facet |
| Interface law | conditions on interface traction at three vertices | the law is a convex set and the traction is linear |

If strength varies inside an element, the argument survives only if the admissible set stays convex in
`(x, sigma)`. Version 1 uses constant strength per element.

---

## 9. Rigorous Upper bound [proposed, open]

Field: P2 velocity on straight tetrahedra [BH21 §4.1; YSK94 and Mak07 for the lineage].

| Condition | Imposed as | Why it is exact or safe |
|---|---|---|
| Straightness | checked; curved tetrahedra refused or explicitly straightened | a curved map makes strain rate non-linear |
| Flow rule | strain rate in the flow cone at four vertices | strain rate is linear; the cone is convex |
| Dissipation | `(V / 4) sum_i d(D_i)` | safe: see below |
| Velocity supports | eliminated or imposed exactly, no penalty | a penalty leaves an inadmissible field |
| Bonded rigid footing | all footprint nodes follow `v0 + omega x (x - x_ref)` | a rigid field is linear, so P2 represents it exactly on a flat facet |
| Work rates | exact integrals | see below |

**Dissipation is over-estimated, never under-estimated.** On a tetrahedron, `D(x) = sum_i L_i(x) D_i`,
with barycentric `L_i`. `d` is convex and positively homogeneous, so `d(D(x)) <= sum_i L_i(x) d(D_i)`.
Integrating with `int L_i dV = V / 4` gives the rule above.

**Work rates must not use the vertex rule.** `int N dV` for P2 is `-V/20` at a vertex and `V/5` at an edge
midpoint. A vertex rule would be wrong for gravity. Traction work of linear `t` against quadratic `v` is a
cubic on the facet and is integrated exactly.

P1 constant-strain tetrahedra are not the baseline. They are the classic locking case for isochoric flow
[NPR74].

---

## 10. Where checking changed the working design (proposed amendments)

The design handed to this ADR is treated as frozen. These items are therefore amendments for the owner and
orchestrator to accept or reject.

### 10.1 Continuous P2 velocity appears to lock for isochoric flow in true 3D

For von Mises and Tresca the flow rule forces `tr(D) = 0`. On a P2 tetrahedron `tr(D)` is linear, so that
is four equality constraints per element. By Euler's relation `E ~ N + T` on a large mesh, so unknowns
minus constraints is about `3(N + E) - 4T ~ 6N - T`. Tetrahedral meshes commonly have `T` near `6N`
(recalled; the meshes below had not reached it). The count leaves almost no margin before boundary
conditions, and the theory of exactly divergence-free elements points the same way: stable quadratic
pairs on tetrahedra rely on special splits [Zha11].

**Scoping probe, 2026-11-18.** Linear algebra only, with no optimiser and no bound computed. It uses no
OpenSees build. It asks whether a continuous P2 field with `tr(D) = 0` at every vertex of every
tetrahedron can move a rigid footprint. Method in Appendix A. numpy 2.3.5, scipy 1.16.3, gmsh 4.15.0.
No output file was kept; the table is the record.

| Problem and mesh | Tets | Free DOF | Constraints | Rank | Nullity | Relative residual of the footing motion |
|---|---:|---:|---:|---:|---:|---:|
| Thin slab, gmsh unstructured, `h = 0.14` | 612 | 2509 | 2448 | 2326 | 183 | 9.3e-15 |
| Thin slab, gmsh unstructured, `h = 0.10` | 909 | 3760 | 3636 | 3478 | 282 | 1.2e-8 |
| Thin slab, structured, `n = 8` | 384 | 1645 | 1536 | 1258 | 387 | 6.6e-15 |
| Thin slab, structured, `n = 16` | 1536 | 6881 | 6144 | 5078 | 1803 | 1.5e-14 |
| 3D square footing, structured, `n = 4` | 384 | 1101 | 1536 | 1024 | 77 | 1.1e-1 |
| 3D square footing, structured, `n = 6` | 1296 | 4209 | 5184 | 3657 | 552 | 1.2e-14 |
| 3D square footing, gmsh unstructured, `h = 0.28` | 384 | 975 | 1536 | 975 | 0 | 1.5e-1 |
| 3D square footing, gmsh unstructured, `h = 0.20` | 734 | 2070 | 2936 | 2070 | 0 | 3.8e-2 |
| 3D square footing, gmsh unstructured, `h = 0.16` | 1563 | 4716 | 6252 | 4715 | 1 | 1.3e-2 |

Reading, within the probe's limits:

- On the unstructured 3D meshes tested, there is essentially no isochoric field, and the footing motion
  is **infeasible**. The residual falls with `h` but is not zero. An interior-point solver with a finite
  tolerance could still report near-feasibility, so the checker must not accept it.
- The structured mesh has redundant constraints and behaves differently, so a structured mesh can hide
  the problem.
- **The thin slab passes.** The planned formulation gate, the thin-slab strip, cannot detect this.

Limits of the probe: small meshes; walls and base fully fixed, with no roller case; a footprint selected by
node position, not conformal to the mesh; `phi = 0` only. Drucker-Prager with `alpha > 0` has an inequality
cone and was not probed. It should be expected to worsen as `alpha` falls, but that is not measured.

**Proposed amendment A1.**

1. Add a P0 gate: true-3D Upper feasibility for `phi = 0` on unstructured tetrahedra, at several sizes,
   next to the thin slab. 237a freezes it and 237d measures it.
2. Treat facet velocity jumps (§10.2) [SK95, Kra05] as a likely requirement for isochoric criteria, not as an option.
   A count for fully discontinuous P2 leaves about `14T` freedom (`30T` unknowns, `4T` isochoric
   constraints, six normal-jump constraints on each of about `2T` facets), but freedom is not accuracy.
3. Evaluate alternatives in 237d: Powell-Sabin-type splits [Zha11] and higher-order velocity. Neither is
   adopted here.
4. If none works, the Upper scope is reduced or reformulated before any C++ (K1, K2).

### 10.2 A P2 velocity jump is quadratic, so facet-vertex checks are not enough

The working design imposes jump and interface cones at facet vertices. That is right for the Lower bound,
where interface traction is linear. It is not right for an Upper jump. The jump of two P2 fields on a flat
facet is quadratic, and so is the difference between P2 soil velocity and a rigid footing. Even all six
nodal values are not enough: with vertex values `1` and midside values `0`, a P2 function on a triangle is
`-1/3` at the centroid.

**Proposed amendment A2.** Impose jump conditions on the six Bernstein coefficients of the quadratic jump:

- `b_i = j_i` at the vertices, and `b_ij = 2 j_ij - (j_i + j_j) / 2` for the edge midpoints.
- The jump is a convex combination of the `b`, because the quadratic Bernstein polynomials are non-negative
  and sum to one. If every `b` lies in the convex jump cone, the whole jump does.
- Each polynomial integrates to `A / 6`, so `(A / 6) sum_b d_j(b)` over-estimates jump dissipation by the
  argument of §9.

This is sufficient, not necessary, and it costs some accuracy. It applies to velocity jumps and to
unbonded interfaces: no tension and Coulomb. Bonded footings need no jump.

### 10.3 `domain_claim` belongs to each bound first

The two bounds extend differently (§12). A box with fixed walls has a collapse load no lower than the
half-space: its mechanisms extend to the half-space by zero. So its Upper bound holds for both, and its
Lower bound holds only for the box.

**Proposed amendment A3.** Each certified bound carries its own `domain_claim`. A bracket's claim is the
strongest domain for which both bounds hold.

### 10.4 Certificate repair needs conditions the working design does not state

- **Lower.** Scaling toward an interior point needs a *statically admissible* point strictly inside `K`
  everywhere. Zero stress with `lambda = 0` works for a cohesive, weightless problem without dead load. With
  self-weight or dead load, it does not satisfy equilibrium, and for `c = 0` zero stress is on the boundary.
- **Upper.** For von Mises and Drucker-Prager, `d` is infinite outside the flow cone. A finite charge
  cannot pay for a cone violation. The field must be corrected into the cone.

**Proposed amendment A4.** The Lower repair uses a separately computed and separately checked interior
field. The Upper repair corrects the velocity field and recomputes `Diss`. Where neither is available the
result stays an estimate, with that reason named.

### 10.5 Stage kinds and the wipe hook

W1 and W2 change the earlier plan for `LB` and `UB` stages and for `wipe`. See D3 and §4.3.

---

## 11. Certificate and terminology [proposed]

### 11.1 Terms

| Term | Meaning |
|---|---|
| **Certified bound** | passed the independent checker |
| **Estimate** | everything else, including a solver's claimed optimum and an uncertified `AlmostSolved`; its reason is stated |
| **Inadmissible** | an estimate that cannot satisfy a gate |
| **Bracket** | a certified Lower and a certified Upper for the same stated problem |
| **Gap** | `(Upper - Lower) / ((Upper + Lower) / 2)` |

Write `[Lower, Upper]` only for a certified bracket. Otherwise write `estimate: L=..., U=...`. "Rigorous" means
rigorous for the stated idealisation, up to the reported floating-point floor. Every result reports its
`domain_claim`, criterion, associated-flow status and interface law.

A gate cannot pass on an estimate. A certified Lower above a certified Upper for one problem stops the
program (K3).

### 11.2 What the checker does

The checker is independent of the reference and of the engine. From the problem file and saved fields it:

1. verifies the problem hash, the labels and that every tetrahedron is straight;
2. rebuilds the operators through its own code;
3. **Lower:** recomputes equilibrium, facet jumps, boundary tractions and resultants, and evaluates yield from
   principal stresses or another route not used by the solver's cone map;
4. **Upper:** recomputes boundary and footing constraints, work rates, flow-cone membership and dissipation from
   the support function, including §10.2 for any jump;
5. checks that `domain_claim` is supported;
6. applies the repairs of §10.4;
7. reports the repaired value and the remaining floating-point floor;
8. checks `Lower <= Upper` for a pair.

Thresholds are named constants, justified by scale and calibrated on exact cases before use. They may not be
changed to rescue a run. A solver status is an input to this process, never a verdict.

### 11.3 False-green defences

The checker must reject fixtures with: one perturbed stress component; stress scaled by 1.01; a facet traction
jump; a wrong load normalisation; swapped Lower and Upper labels; a curved tetrahedron; a mismatched problem hash;
an unsupported `domain_claim`. A test reads the run record and compares it with the requested case, so a plausible
number from the wrong case fails.

---

## 12. Domain claims [proposed, O5]

| Value | Meaning |
|---|---|
| `box_stated` | a bound for the stated box and its exact wall conditions |
| `unverified` | extension to a wider domain is not established; for that wider domain the value is an estimate |
| `halfspace_by_zero_extension` | an Upper field with zero velocity on all truncation faces extends by zero |
| `halfspace_by_extension_elements` | a Lower field has a verified admissible extension |

Rules:

- Roller walls can permit tangential motion, so an Upper bound with rollers is not automatically a half-space bound.
- A Lower bound in a box with fixed walls is a bound for that box.
- With self-weight, a traction-free truncated Lower field does not automatically extend.
- A reference solution is compared only with a claim that covers its domain. A box used against a half-space value
  must also contain the classical mechanism.

Version 1 produces box-labelled bounds, then an extent study. Extension elements come later.

---

## 13. Query design [proposed]

Let `F` be the six-component footing resultant about a declared reference point, with moments divided by a declared
length. Let `K_F` be the set of resultants the soil can carry with dead loads fixed. `K_F` is convex, since it is the
linear image of a convex set.

| Query | Lower | Upper |
|---|---|---|
| `support(d)` | maximise `d . F`; gives a point `F_L` inside `K_F` | prescribe footing velocity `q = d`; `Diss - W0` gives the half-space `d . F <= Diss - W0` containing `K_F` |
| `ray(F0, F1)` | maximise `lambda` with `F = F0 + lambda F1` | `q` free with `F1 . q = 1`; `lambda <= Diss - W0 - F0 . q` |
| fixed components | a ray, with the fixed components in `F0` | same |

Report `conv(lower points) subset K_F subset intersection(upper half-spaces)`. Every Upper solution also yields a
half-space with its own `q`, which can be added to the outer set.

Assembled operators are reused across queries on one mesh. A direction-specific adaptive mesh is allowed, because any
valid mesh gives a valid bound. Bounds combine only for exactly the same stated problem: geometry, material, interface,
reference point and `domain_claim`.

Interface laws:

| Law | Lower | Upper | Note |
|---|---|---|---|
| Bonded | no extra law | footprint follows the rigid motion | — |
| No tension | linear inequality at three vertices | jump condition, §10.2 | — |
| Coulomb | cone at three vertices | associated jump cone, §10.2 | the theorem treats it as associated, so sliding dilates |

For real non-dilatant friction, an associated interface result is not a bound on the physical problem.

---

## 14. Open-solver policy

### 14.1 Rule [owner]

No component that needs a licence request enters the fork for this feature, and no required gate needs one. No MOSEK
code, adapter, licence check or test requirement enters the fork, nor any equivalent commercial solver. The fork's
existing host toolchain, including Intel oneAPI and MKL where its build already uses them, is outside this rule. New
limit-analysis code does not depend on proprietary solvers.

### 14.2 Mechanism [proposed]

- **D5.** The default build links no solver and needs no Rust.
- **D6.** A run-time adapter with a small versioned C ABI connects a solver. `getLibraryFunction` is the precedent for
  loading it.
- **D7.** `limitAnalysis export` always works, so a run can be solved outside the process and checked.
- How users obtain an adapter is open (O15).

### 14.3 Candidates

Licences come from the GitHub licence API on 2026-11-18 unless stated. Re-check before adding or shipping anything.
This table was built from licence metadata and project documentation; no GPL or LGPL source was opened.

| Tool | Licence | Role |
|---|---|---|
| Clarabel.rs, Clarabel.cpp [GC24] | Apache-2.0 | primary interior-point candidate, with second-order and semidefinite cones |
| faer | MIT | Clarabel's optional supernodal direct solver |
| OpenBLAS | BSD-3-Clause | BLAS/LAPACK option for Clarabel's semidefinite support (`sdp-openblas`) |
| QOCO [CA25] | BSD-3-Clause | candidate independent second-order-cone solver; not yet tried here |
| SCS [OCPB16] | MIT | first-order cross-check and large-problem experiment; not a certificate by itself |
| HiGHS | MIT | linear-program tasks only |
| CSDP | EPL-2.0 | possible small semidefinite cross-check in the harness |
| ECOS | GPL-3.0 | optional out-of-process comparison on exported programs; never linked, never required, source not read |
| MUMPS | CeCILL-C (`OTHER/MUMPS/include/dmumps_c.h`) | sparse linear solver for a possible tailored method |

Clarabel build features `pardiso-panua` (its README says a separate licence is required) and `sdp-mkl` are out of policy
for adapters. The Python wheel in the development environment, version 0.11.1, solved a trivial problem with `qdldl` and
`faer` and reported `mkl` and `panua` as not available.

No claim is made here about which solver reaches certification at scale. That is 237d's question.

### 14.4 A tailored solver

Consider a tailored interior-point method on MUMPS only when all three hold:

1. after P0's workstation test, open general solvers cannot produce a certified result on the reference 3D footing within
   one cluster node's time and memory;
2. at least 80% of solve time is in the KKT factorisation and solve;
3. a standalone test on the same reduced system is at least three times faster.

---

## 15. FEniCSx stays outside the fork [owner]

FEniCSx may serve as an independent test reference with open solvers. It is never a fork dependency.

- The lane, 237e, lives outside the fork and outside apeGmsh; its location is O11.
- It assembles the same frozen cases independently and solves them with open solvers. Agreement is judged on the same
  exported problem class and exact cases.
- DOLFINx and UFL report LGPL-3.0 and Basix reports MIT. FFCx reports no SPDX match, so its licence file must be read first.
- That lane may read FEniCS code. Fork workers may not. No code passes from that lane into the fork.
- No fork gate requires it.

---

## 16. Verification and warrant [proposed]

The harness follows `Ladruno_implementation/testbed/00_canonical_testbed.md` and ADR-87.

| Level | Purpose |
|---|---|
| L0 | operators: constant and linear stress, rigid motion, footing resultants, cone membership |
| L0b | from P1: Python and C++ programs match entry by entry |
| L1 | checker: good fixture passes, §11.3 fixtures fail |
| L2 | identical exported programs through independent open solver paths |
| L2e | external FEniCSx assembly, when active, outside the fork |
| L3 | exact reductions and property gates |
| L4 | true 3D: property gates, **the §10.1 gate**, later Tresca references |
| L5 | associated OpenSees models compared with certified brackets |
| L6 | adaptive: gap reduction, re-resolved declarations, budget and size-field fidelity |

Property gates, where a case supports them: certified Lower not above certified Upper; nested uniform refinement does not
lower the Lower or raise the Upper; rigid rotation of the problem leaves both unchanged; doubling cohesion on a weightless
problem doubles both; a valid half model matches the full model. Reference gates are one-sided, plus a maximum gap, so a
bound cannot pass from the wrong side.

Reference cases enter from a checked source, with date, units, load definition and tolerance. First case: the plane-strain
strip with `(2 + pi) cu` [Pra21], as a thin slab with a matched von Mises strength.

Placement (O1): the frozen contract, Python reference and checker in `tests/`, with solver-free tests in Zone-A; solver
studies in `Ladruno_files/testbed/limit_analysis/`.

Warrant for the C++ family at 237f: non-self-referential oracles; a `LIMIT` mutation family through
`CMakeLists.txt::LADRUNO_MUTATE_FAMILY`, `Ladruno_scripts/mutation_gate.py::FAMILIES` and `SRC/Ladruno_mutation.h`; a guide;
both shipped interpreters or a recorded absence; ledger fragments; a banner line; and a manifest row added deliberately, with
`tag_symbol: null` if there is no class tag (W3). Mutants should fail when facet continuity is dropped, gravity is zeroed, the
resultant map changes sign, or a cone map is weakened.

---

## 17. Dependencies

### 17.1 apeGmsh ADR 0136 (Level A) [proposed]

Limit analysis is stateless, so each new mesh reruns the whole analysis. That is Level A of ADR 0136. The draft needs these
changes before it can serve this program. They are requests to that ADR's process, through T; this ADR does not decide them.

| Item | Draft text at `668b99a44` | Needed |
|---|---|---|
| Estimator decoupling | D1's driver "estimates"; S3's gates are built on Kirsch and ZZ-based targets | the driver takes an estimator callable, so its core does not need ZZ |
| A1 indicator input | D3 `results.error.indicator(field, kind=...)` ranks by a recorded field; INV-6 refuses "a quadratic element in version 1" for the estimator | an element array can drive size, so tet10 limit fields never need ZZ |
| A2 stop predicate | the loop targets an error level | a study-specific stop: certified gap, budget, stagnation, iteration limit |
| A3 coarsening opt-in | D2: "It refines only" | add when budget studies need it |
| A4 3D size-field gate | S2 gates name a prescribed size, the model-data refusal and a tensor path | a scalar 3D list-data gate with budget and fidelity |
| A5 portable size field | `SizeField` is built from the old mesh | one study can seed another |

Also required on the limit-analysis side: labels, supports, loads and footprints are verified after each rebuild; second-order
curving is handled explicitly, because the bounds need straight tetrahedra.

Kept out of ADR 0136: bound, certificate, pair and solver concepts; a fork dependency in S1 to S3; two meshes per variant row;
in-session remeshing; footing-specific rules; 3D anisotropy in the first slice.

Adaptive loop (237h): build a fresh session, write the problem file, run both bounds, certify, compute element values from
certified fields, test the stop, build the next size field. A failed certification is logged and not counted. The best certified
Lower and lowest certified Upper are kept only across meshes of exactly the same stated problem. Every indicator is unverified
until tested for this formulation; published adaptive schemes [Lya05, Cir08, Mun08] are starting points only.

### 17.2 apeGmsh bridge [proposed]

A separate apeGmsh ADR, ratified by its maintainer, covers the bridge. ADR 0134 is the precedent: emission only, a sidecar written
together with the deck or not at all, new census families by exact entry, primitives only after the fork feature is merged and in a
release named by `FORK_PIN`, and a member Warranted only when its own `live-fork` test passes (ADR 0121 D2 and D3). All apeGmsh work
goes through its T coordinator. Whether to use new families or a token under `analysis` is O6.

### 17.3 Fork ADR 236 [proposed]

**Level A does not depend on ADR 236.** The engine uses no material history, no `getHistory` or `setHistory`, and no commit primitive.
No limit-analysis phase waits for it.

The restart lane is coordinated separately: ADR 236 S1 to S5 and ADR 0136 S4 to S6. ADR 232's FE-to-MPM handoff waits on the same
primitive, so ADR 236 S1 must be scheduled with the MPM program, without duplicate work. The open questions in ADR 236 §6 stay with
that ADR.

---

## 18. Phases and kill criteria [proposed]

| Phase | WP | Delivers | Exit |
|---|---|---|---|
| P0 | 237a to 237e | contract, Python reference, checker, solver study, external reference | L0 to L3 pass; certified thin-slab bracket contains `(2 + pi) cu`; **§10.1 gate decided**; property gates pass where they apply; study states which open path certifies at what size; exported programs replay through a second solver path; any FEniCSx mismatch resolved or recorded before I1 |
| P1 | 237f | C++ kernel, assembler, command, export, writer, adapter | L0b identity; checker certifies engine output; mutation lane fails on operator faults; disabled build unchanged; shipped interpreters handled or absence recorded |
| P2 | L1 to L4, apeGmsh | bridge and reader | emitter parity; stock refusal; sidecar atomicity; `Results.limit` fixtures; `live-fork` proof on the pin |
| P3 | 237g | V-H-M and 3D queries | strip V-H and V-M certified brackets; sliding correct for the stated interface; query reuse measured; one six-component 3D demonstration labelled by criterion and domain |
| P4 | 237h | adaptive 3D | certified gap falls at fixed budget on references; every counted iteration passes the checker and re-resolution; best bounds tracked correctly; 3D size-field fidelity proven |
| P5 | 237i | semidefinite cones, jumps, facet output | cone maps agree with principal-stress evaluation; rough circular Tresca bracketed against its checked reference [ES60]; Mohr-Coulomb reduces correctly to plane-strain references; open-solver cost measured |
| P6 | 237j | performance and advanced meshing | each item shows measured gain over P4 |

If A1 and A2 are accepted, facet jumps move from P5 toward P0 or P1 for isochoric criteria.

| Kill | Trigger | Action |
|---|---|---|
| K1 | no open path certifies 3D Drucker-Prager Lower and Upper at about 20,000 tetrahedra within one hour and one workstation's memory. The threshold is provisional | stop before C++; reformulate, rescale, prototype a tailored open solver, or reduce scope |
| K2 | the thin-slab strip fails after operator tests pass; **proposed addition:** the §10.1 gate fails and no accepted remedy exists | stop and rethink the formulation |
| K3 | certified Lower above certified Upper for one problem | stop the program |
| K4 | solver integration disturbs the default build | keep or repair the run-time adapter before going on |
| K5 | after adaptivity, the certified gap cannot separate the result it must judge | label the 3D result limited evidence and ask the owner |

P0 records the fixed-mesh gap but does not kill on it. Thresholds are not tuned to reach a reference, and neither is an element,
tolerance or mesh rule.

---

## 19. Risks, decisions, sources

### 19.1 Risks

| Risk | Consequence | Mitigation |
|---|---|---|
| Upper locking in 3D (§10.1) | no usable Upper for low-friction criteria on general meshes | P0 gate, jumps, alternatives, K2 |
| Interior-point accuracy | solutions may not repair into certificates | scaling, formulation work, checker as the only judge, K1 |
| 3D program size | memory and time limit the mesh, so gaps stay wide | 237d measures; adaptivity; §14.4 |
| Wide fixed-mesh gaps | a bracket is true but not useful | P4; K5 |
| Thin-slab false confidence | P0 passes and true 3D fails | §10.1 gate |
| Reference and engine share an error | agreement proves nothing | separate authors, frozen contract, independent checker, external lane |
| Misread as a soil-model result | an associated bound is applied to a non-associated problem | required labels (§11.1), non-goals |
| Domain misstatement | a box result is cited as a half-space result | per-bound claims (§10.3, §12) |
| Reader breakage | `.ladruno` additions break existing readers | D3, fixtures before relying on readers |
| Licence drift | a dependency changes terms | re-check before adding or shipping (§14.3) |
| ADR 0136 changes shape | the adaptive phase slips | P0 to P3 do not depend on it |

### 19.2 Open owner decisions

None of these is decided by this ADR.

| ID | Question | Recommendation | Needed by |
|---|---|---|---|
| O1 | Confirm P0 placement (§16) | as stated | first push |
| O2 | Confirm WP-237 onward, fork ADR 237, an apeGmsh limit ADR through T, and amending ADR 0136 | yes | program issue |
| O3 | Separate problem file for the engine | yes | P1 |
| O4 | Von Mises and Drucker-Prager first; Tresca and Mohr-Coulomb in P5 | yes, subject to A1 | P1 |
| O5 | Box-labelled bounds first, extent study, extension elements later | yes | any half-space claim |
| O6 | New apeGmsh families or a token under `analysis` | new families | L1 |
| O7 | Amend apeGmsh#3052 before acceptance, or follow up | before, if T agrees | C0 |
| O8 | Is an ADR-87 D4 out-of-family verdict required | decide with this ADR | 237f merge |
| O9 | Classic Tcl wired or recorded absent | wire if inexpensive | 237f |
| O10 | A separate licence for the OpenSees-free kernel | keep the option open; seek advice first | external release |
| O11 | Where the FEniCSx lane lives | separate private repository or folder | 237e |
| O12 | Do fork merges stay with the owner, or are they delegated | owner unless the program issue says otherwise | first ready-flip |
| O13 | ADR 236 §6 questions | decide in D0; coordinate S1 with MPM | D1 |
| O14 | ADR 0136 questions: Level A first, driver location, Bezier gate | Level A first; driver package; measure Bezier | C0 |
| O15 | How users obtain a solver adapter: separate release, optional build, or user-built | decide after 237d | 237f |
| O16 | Accept amendments A1 to A4 (§10) | accept; A1 first, because it can change P0 | 237a freeze |

### 19.3 Unverified or partly verified

- **Literature.** Bibliographic data below were checked against Crossref or arXiv on 2026-11-18. Of the papers, only [BH21]
  was read for content, in its arXiv version. Attributions to the others are from memory and are marked `(m)` in the list.
- `T` near `6N` for production tetrahedral meshes (§10.1): recalled.
- Applicability of [Zha11] to this problem: title and metadata only.
- Solver reach at scale, `AlmostSolved` behaviour, and any earlier review probe of Upper or P1 elements: no reproducible record in
  this repository, so no numbers are given. 237d must measure them.
- The rough circular reference value is not stated here. It is entered in P5 after checking [ES60].
- Square-footing references [Sal04, Lya07] are candidates only; values and conditions are unchecked.
- A self-weight strip reference (Martin's `N_gamma`, IACMAG 2005) was not found in Crossref and is unverified.
- Licences of FFCx and Mmg (GitHub reports no SPDX match) and gmsh, CVXOPT and Hypatia: not verified here.
- Reader behaviour for new `.ladruno` element groups and for `MODEL/LIMIT_ANALYSIS`: untested.
- Drucker-Prager locking: not probed.

### 19.4 References

| Key | Reference |
|---|---|
| DPG52 | Drucker, Prager, Greenberg (1952). Extended limit design theorems for continuous media. *Q. Appl. Math.* 9(4), 381-389. doi:10.1090/qam/45573 (m) |
| Pra21 | Prandtl (1921). Über die Eindringungsfestigkeit (Härte) plastischer Baustoffe und die Festigkeit von Schneiden. *ZAMM* 1. doi:10.1002/zamm.19210010102 (m) |
| ES60 | Eason, Shield (1960). The plastic indentation of a semi-infinite solid by a perfectly rough circular punch. *ZAMP* 11(1), 33-43. doi:10.1007/BF01591800 (m) |
| NPR74 | Nagtegaal, Parks, Rice (1974). On numerically accurate finite element solutions in the fully plastic range. *CMAME* 4(2), 153-177. doi:10.1016/0045-7825(74)90032-2 (m) |
| Sloan88 | Sloan (1988). Lower bound limit analysis using finite elements and linear programming. *IJNAMG* 12(1), 61-77. doi:10.1002/nag.1610120105 (m) |
| DD93 | Drescher, Detournay (1993). Limit load in translational failure mechanisms for associative and non-associative materials. *Géotechnique* 43(3), 443-456. doi:10.1680/geot.1993.43.3.443 (m) |
| YSK94 | Yu, Sloan, Kleeman (1994). A quadratic element for upper bound limit analysis. *Eng. Comput.* 11(3), 195-212. doi:10.1108/02644409410799281 (m) |
| SK95 | Sloan, Kleeman (1995). Upper bound limit analysis using discontinuous velocity fields. *CMAME* 127(1-4), 293-314. doi:10.1016/0045-7825(95)00868-1 (m) |
| LS02a | Lyamin, Sloan (2002). Lower bound limit analysis using non-linear programming. *IJNME* 55(5), 573-611. doi:10.1002/nme.511 (m) |
| LS02b | Lyamin, Sloan (2002). Upper bound limit analysis using linear finite elements and non-linear programming. *IJNAMG* 26(2), 181-216. doi:10.1002/nag.198 (m) |
| Kra05 | Krabbenhøft, Lyamin, Hjiaj, Sloan (2005). A new discontinuous upper bound limit analysis formulation. *IJNME* 63(7), 1069-1088. doi:10.1002/nme.1314 (m) |
| Lya05 | Lyamin, Sloan, Krabbenhøft, Hjiaj (2005). Lower bound limit analysis with adaptive remeshing. *IJNME* 63(14), 1961-1974. doi:10.1002/nme.1352 (m) |
| Mak06 | Makrodimopoulos, Martin (2006). Lower bound limit analysis of cohesive-frictional materials using second-order cone programming. *IJNME* 66(4), 604-634. doi:10.1002/nme.1567 (m) |
| Mak07 | Makrodimopoulos, Martin (2007). Upper bound limit analysis using simplex strain elements and second-order cone programming. *IJNAMG* 31(6), 835-865. doi:10.1002/nag.567 (m) |
| Kra07 | Krabbenhøft, Lyamin, Sloan (2007). Formulation and solution of some plasticity problems as conic programs. *IJSS* 44(5), 1533-1549. doi:10.1016/j.ijsolstr.2006.06.036 (m) |
| Kra08 | Krabbenhøft, Lyamin, Sloan (2008). Three-dimensional Mohr-Coulomb limit analysis using semidefinite programming. *Commun. Numer. Methods Eng.* 24(11), 1107-1119. doi:10.1002/cnm.1018 (m) |
| MM08 | Martin, Makrodimopoulos (2008). Finite-element limit analysis of Mohr-Coulomb materials in 3D using semidefinite programming. *J. Eng. Mech.* 134(4), 339-347. doi:10.1061/(ASCE)0733-9399(2008)134:4(339) (m) |
| Cir08 | Ciria, Peraire, Bonet (2008). Mesh adaptive computation of upper and lower bounds in limit analysis. *IJNME* 75(8), 899-944. doi:10.1002/nme.2275 (m) |
| Mun08 | Muñoz, Bonet, Huerta, Peraire (2008). Upper and lower bounds in limit analysis: adaptive meshing strategies and discontinuous loading. *IJNME* 77(4), 471-501. doi:10.1002/nme.2421 (m) |
| Sal04 | Salgado, Lyamin, Sloan, Yu (2004). Two- and three-dimensional bearing capacity of foundations in clay. *Géotechnique* 54(5), 297-306. doi:10.1680/geot.2004.54.5.297 (m) |
| Lya07 | Lyamin, Salgado, Sloan, Prezzi (2007). Two- and three-dimensional bearing capacity of footings in sand. *Géotechnique* 57(8), 647-662. doi:10.1680/geot.2007.57.8.647 (m) |
| Zha11 | Zhang (2011). Quadratic divergence-free finite elements on Powell-Sabin tetrahedral grids. *Calcolo* 48(3), 211-244. doi:10.1007/s10092-010-0035-4 (m) |
| Slo13 | Sloan (2013). Geotechnical stability analysis. *Géotechnique* 63(7), 531-571. doi:10.1680/geot.12.RL.001 (m) |
| OCPB16 | O'Donoghue, Chu, Parikh, Boyd (2016). Conic optimization via operator splitting and homogeneous self-dual embedding. *JOTA* 169(3), 1042-1068. doi:10.1007/s10957-016-0892-3 (m) |
| BH21 | Bleyer, Hassen (2021). Automated formulation and resolution of limit analysis problems. *Comput. Struct.* 243, 106341. doi:10.1016/j.compstruc.2020.106341; arXiv:2005.04779 |
| GC24 | Goulart, Chen (2024). Clarabel: an interior-point solver for conic programs with quadratic objectives. arXiv:2405.12762 (m) |
| CA25 | Chari, Açıkmeşe (2025). QOCO: a quadratic objective conic optimizer with custom solver generation. arXiv:2503.12658 (m) |

---

## Appendix A. Method of the §10.1 probe

The probe is throwaway and is not tracked; this recipe is the record. It is not a reference implementation.

1. **Unknowns.** Three velocity components at every vertex and edge midpoint of straight tetrahedra.
2. **Constraint rows.** For each tetrahedron and each of its four vertices, one row for `tr(D) = 0`, from the gradients of the P2 shape
   functions at that vertex.
3. **Unit-cube footing.** Base and four sides fully fixed. Footprint: nodes on `z = 1` with `0.25 <= x, y <= 0.75`, prescribed to
   `v = (0, 0, -1)`. The rest of the top is free.
4. **Slab.** `0 <= x, z <= 1`, thickness `h` in `y` (`1/n` for the structured mesh). `v_y = 0` on both `y` faces; `x = 0`, `x = 1` and
   `z = 0` fully fixed; footprint nodes on `z = 1` with `0.375 <= x <= 0.625` prescribed to `v = (0, 0, -1)`.
5. **Meshes.** Structured: each cube cell split into six tetrahedra around its main diagonal. Unstructured: gmsh 4.15.0 default 3D
   meshing of an OpenCASCADE box, with minimum and maximum size both `h`.
6. **Measures.** With `A` the free-column matrix and `b` the prescribed-motion right-hand side: rank from singular values above
   `1e-10` of the largest; nullity = free DOF minus rank; residual = `|A x - b| / |b|`, with `x` the least-squares solution at
   `rcond = 1e-10`.

---

## Review log

- **2026-11-18, drafting.** Written from the owner's decisions and the working design, then checked against the fork at `cfe478254`.
  The check produced W1 to W5 and amendments A1 to A4. No independent review has been run on this draft.

## Implementation log

Nothing built.
