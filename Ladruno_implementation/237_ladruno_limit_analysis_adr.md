---
title: "ADR 237 — Certified finite-element limit analysis in the Ladruno fork (WP-237)"
project: Ladruno
type: ADR
status: "Proposed (2026-11-18) — design only, no engine code. Amendments A1 to A4 accepted by the owner"
priority: high
owner: fork maintainer
wp: "WP-237 (claimed by program issue #1118)"
related:
  - "Program issue #1118: claims WP-237 onward"
  - "[[_wp237_orchestration_plan]]: work packages, roles, gates"
  - "[[87_ladruno_depth_with_width_adr]]: warrant package, merge rules"
  - "[[232_ladruno_mpm_general_driver_adr]]: the engine-with-its-own-file precedent"
  - "Fork ADR 236 (draft PR #1117): material history transfer; not a dependency of Level A"
  - "apeGmsh ADR 0136 (draft PR apeGmsh#3052): adaptive mesh convergence; needed from P4"
tags: [adr, limit-analysis, conic-optimization, bounds, certificate, apegmsh, wp-237]
updated: 2026-11-18
---

# ADR 237 — Certified finite-element limit analysis

> **Status.** Proposed. No engine exists. Each decision carries one of two marks: **[owner]** was decided
> by the owner on 2026-11-18; **[proposed]** is a working default or a recommendation of this document
> and stays open until accepted.
>
> **Number.** Program issue #1118 claims WP-237 onward, including the lettered packages.
>
> **Checked against.** Fork claims: `cfe478254` (`origin/ladruno` on 2026-11-18), cited as `path::symbol`.
> Fork ADR 236: PR #1117 head `f7f89a4c6`. apeGmsh ADR 0136: PR apeGmsh#3052 head `668b99a44`. apeGmsh
> ADRs 0121 and 0134: `main` at `60d0517bc`. Only apeGmsh ADR text was read, not its source. The draft
> P0 contract is `tests/limit_interface.md` on `wp/237a-limit-contract` at `7909baa4e`. It is not on this
> branch and is cited as "draft contract"; it may change before the freeze.
>
> **Discipline.** Under ADR-87 an unmeasured claim is a claim. §19.3 lists what could not be verified.

---

## 0. Decision summary

| Question | Answer | Mark |
|---|---|---|
| What is built | 3D finite-element Lower and Upper bounds for perfectly plastic, associated solids, with an independent certificate | [owner] 3D; [proposed] the rest |
| Where | In the fork, as a one-shot engine behind one `limitAnalysis` command row, present in every build | [owner] fork; [proposed] the shape |
| Input | A separate versioned problem file `<stem>.la.h5`, not the Domain | [proposed] O3 |
| Output | The engine's own `.ladruno` writer, not `recorder ladruno` | [proposed] |
| Materials | A `LimitMaterial` registry, separate from `NDMaterial`. Von Mises and Drucker-Prager first; Tresca and Mohr-Coulomb in P5 | [proposed] O4 |
| Lower | Discontinuous linear stress on straight tetrahedra (§8) | [proposed] |
| Upper | Continuous P2 is a **control**. P0 must carry a jump-based isochoric formulation or a proven alternative (§9) | [owner] the gate and the whole-facet rule; [proposed] the formulation |
| True-3D gate | Non-structured meshes, at least three sizes, fixed and roller walls, conformal footprint, rank and residual, then certified quality (§10.1) | [owner] the gate; [proposed] its detail |
| Certified | An exact affine projection, then admissible inequality repair, then every check again (§11) | [owner] admissible-only repair; [proposed] the procedure |
| Domain | Each bound carries its own `domain_claim`; box-labelled bounds first (§12) | [owner] per bound; [proposed] O5 |
| Solvers | No component that needs a licence request. No solver linked by default (§14) | [owner] the rule; [proposed] the mechanism and the toolchain reading |
| FEniCSx | An external test reference, outside the fork, in a separate private repository (§15) | [owner]; [proposed] its licence and reading rules |
| Remeshing | Level A through apeGmsh ADR 0136, with requests R1 to R6 (§17.1) | [owner] central; [proposed] the requests |
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
| `SRC/interpreter/LadrunoCommandTable.h::LADRUNO_COMMAND`, row `"mpm"` | one row registers the command in all engines | see W6. No workflow runs `ci/check_ladruno_commands.py` directly; the Zone-A test `tests/test_ladruno_command_registry.py` calls it and is the CI enforcement |
| `SRC/domain/pattern/ladrunoMPM/LadrunoMPMCore.h` | precedent for a pure kernel: its header says "NO OpenSees dependency", and it is unit-tested with `g++` alone | it depends on the fork's hypo-kernel header; the limit kernel would not need that |
| `SRC/domain/pattern/ladrunoMPM/LadrunoMPMIO.h::readParticleSet`, `::ParticleWriter` | precedent for a contract-checked sidecar with rule codes, and for an engine-owned `.ladruno` writer | both need `_HDF5` and refuse without it |
| `tests/mpm_interface.md`, `tests/mpm_cases.txt`, `tests/mpm_reference.py`, `tests/mpm_contracts.py` | precedent for a frozen interface, cases, an oracle and an executable contract | the C++ author must not read the oracle |
| `SRC/utility/LadrunoThreads.h::ladruno_parallelBuild`, used in `SRC/domain/pattern/ladrunoMPM/OPS_LadrunoMPM.cpp` | refusal in OpenSeesSP and OpenSeesMP | — |
| `SRC/api/packages.cpp::getLibraryFunction` | existing run-time loader for a library and symbol, on Windows and POSIX | it appends the platform suffix and was written for element and material packages; it is a precedent, not yet an adapter loader |
| `CMakeLists.txt`, `option(LADRUNO_CMS ...)` | precedent for an opt-in library; its comment requires that OFF leaves stock and serial targets unchanged at source-list level | it is OFF by default. An option cannot remove a command row (W6) |
| `Ladruno_implementation/ladruno_schema_v1.md` §2a, §3, §7.5, §9 | `.ladruno` layout, element self-description, the additive `ON_PARTICLES` precedent, version rules | §7.3 still labels `ON_DOMAIN` as deferred, although `SRC/recorder/LadrunoRecorder.cpp` creates the group; the text needs reconciling before this engine relies on it |
| `Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py::validate` | executable schema | see W1 |
| `SRC/system_of_eqn/linearSOE/mumps/` | an open sparse solver in some fork builds, for a possible tailored method | parallel targets only, when `MUMPS_DIR` is defined; see §14.4 |

| # | Wall | Consequence |
|---|---|---|
| W1 | `ladruno_format.py::validate` rejects any stage `KIND` other than `transient`, `static` or `eigen`. apeGmsh ADR 0134 records that its TypeScript reader throws a whole-file layout error on an unknown `KIND` | a `lower` or `upper` stage kind would break both. §6 keeps `KIND = "static"` |
| W2 | The engine has no Domain host, so `Domain::clearAll` will not reach it. All three wipe paths (`SRC/tcl/commands.cpp::wipeModel`, `SRC/interpreter/OpenSeesCommands.cpp::OpenSeesCommands::wipe`, `SRC/interpreter/PythonAnalysisBuilder.cpp::ops_wipeModel`) call `OPS_clearAllNDMaterial` | `mpm` is reached through its `LoadPattern` host; this engine needs an explicit hook (§4.3) |
| W3 | `ci/check_manifest.py::main` is driven by new class tags. A feature with no class tag is not forced to have a manifest row | ADR-87 D1 would pass silently here unless the row is added deliberately (§16.3) |
| W4 | Zone-A installs `pytest pyyaml numpy sympy scipy h5py conan` (`.github/workflows/ladruno.yml`) and no conic solver | operator and checker tests can run in Zone-A; solver-dependent tests use a separate required lane (§16.2) |
| W5 | The repository's `COPYRIGHT` is the University of California text, limiting use and distribution by purpose | any bundling of third-party code needs a licence review; see O10. This is an engineer's reading, not legal advice |
| W6 | `LadrunoCommandTable.h` has no conditional rows. Each row is expanded in every build, and `tests/test_ladruno_command_registry.py` requires the shipped engines to register exactly what the columns say | `limitAnalysis` cannot appear only in some builds. It stays one row and refuses at run time when it cannot work (§4.3) |

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

- **The row is in every build.** The table has no conditional rows, so the entry points named by the row are always compiled and linked.
- **Unavailable means a run-time refusal.** The engine's implementation may sit behind a CMake option. If it is compiled out, or HDF5 or a
  required adapter is missing, `limitAnalysis` refuses with a message naming what is missing, as `mpm particles` does without `_HDF5`. It also
  refuses in OpenSeesSP and OpenSeesMP.
- **What an engine-off build changes:** the row and its refusing entry point. Existing commands do not change.
- `export` writes the exact conic program wherever the engine is compiled, so independent solvers can run identical data.
- `wipe` must clear any process-wide state. The proposed route is a local-extern reset called from
  `SRC/material/nD/NDMaterial.cpp::OPS_clearAllNDMaterial`, the hook used by WP-104 and already reached
  by all three wipe paths. That is one upstream edit with a `vanilla/` ledger row. If the engine holds
  no process-wide state, no edit is needed; 237f decides.
- Classic Tcl is either wired or recorded as absent through the row's classic column (O9).

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
element) [BH21 §4.2; Sloan88 and Mak06 for the lineage]. The draft contract names it
`lower.p1d-vertex.v0`.

| Condition | Imposed as | Why it is exact |
|---|---|---|
| Equilibrium | `div sigma + b = 0` once per element, `b` constant per element | linear stress has constant divergence |
| Interior facet | `sigma n` continuous at the facet's three vertices | traction is linear on a flat facet |
| Traction facet | `sigma n = t0 + lambda t1` at three vertices | exact only if `t0`, `t1` are linear on the facet; the contract requires that |
| Free facet | `sigma n = 0` at three vertices | same |
| Symmetry or roller | zero tangential traction at three vertices | same |
| Fixed facet | no traction condition | a support may supply any traction |
| Yield | `sigma in K_e` at four vertices | `K_e` is convex and the element is the convex hull of its vertices |
| Footing resultant | closed-form integral of the linear interface traction | force and moment integrands are polynomial on a flat facet |
| Interface law | conditions on interface traction at three vertices | the law is a convex set and the traction is linear |

**Constant-strength conditions, for both bounds.** `k`, `alpha` and the body force are constant in each
element, and each element belongs to one region. No element straddles a material boundary. Spatially
varying strength is represented by piecewise-constant element values, and the bound is then a bound for
that piecewise-constant problem. A facet between two regions is covered in §9.2.

---

## 9. Rigorous Upper bound

### 9.1 Conditions every Upper formulation must meet [proposed]

The velocity is P2 on straight tetrahedra [BH21 §4.1; YSK94 and Mak07 for the lineage].

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
midpoint. Traction work of linear `t` against quadratic `v` is cubic on the facet and is integrated exactly.

### 9.2 Formulations in P0

Formulation ids follow the draft contract. A new formulation gets a new id, so existing cases keep their
meaning.

| Id | Field | Role in P0 |
|---|---|---|
| `upper.p2c-vertex.v0` | continuous P2 | **control**. It is run and reported, and it is not expected to pass the true-3D gate for `phi = 0` |
| `upper.p2d-jump-bernstein.v0` | P2 per element, with charged facet jumps | **candidate** for isochoric criteria. The draft contract still refuses it, so 237a must admit it or replace it before the freeze |
| another element | must be proven by the same gates | admissible alternative, for example a split-element or higher-order space [Zha11]. Nothing is adopted here |

**Ownership.** 237a freezes the jump-based formulation or a proven alternative. 237b assembles it, 237c
certifies it, and 237d measures it through the true-3D gate. It is not deferred to P5 or 237i.

**Jump formulation [proposed]**, after the discontinuous-velocity line of work [SK95, Kra05]. On a facet between elements `e-` and `e+`, `n` is the outward normal of
`e-` and `j = v+ - v-`.

- The jump acts like a thin layer with strain rate `sym(j (x) n)`, so the jump dissipation is
  `d(sym(j (x) n))`.
- Von Mises gives zero normal jump and dissipation `k norm(j_t)`. The normal jump is a quadratic that
  vanishes when its six nodal values vanish, so that condition is an exact linear equality.
- Drucker-Prager gives a convex jump cone. It is imposed on the Bernstein coefficients of §10.2.
- On a facet between regions, either side's strength is valid, because any admissible stress field satisfies
  yield on both sides. The contract names one side; the draft contract uses `e-`.
- Facets on fixed walls and on bonded footprints keep exact nodal motion unless 237a admits a boundary jump
  with a stated proof.

A count gives about `14T` free unknowns for von Mises (`30T` unknowns, `4T` isochoric constraints, six
normal-jump constraints on each of about `2T` facets). Freedom is not accuracy, and the jump cones add
solver cost. Neither is established here.

P1 constant-strain tetrahedra are not a baseline. They are the classic locking case for isochoric flow
[NPR74].

---

## 10. Accepted amendments

The owner accepted A1 to A4 on 2026-11-18. The acceptance covers the four rules below. The detailed
formulation, thresholds and case values remain proposed until 237a is frozen.

### 10.1 A1: a true-3D Upper feasibility gate [owner]

For von Mises and Tresca the flow rule forces `tr(D) = 0`. On a P2 tetrahedron `tr(D)` is linear, so that is
four equality constraints per element. By Euler's relation `E ~ N + T` on a large mesh, so unknowns minus
constraints is about `3(N + E) - 4T ~ 6N - T`. Tetrahedral meshes commonly have `T` near `6N` (recalled), which
leaves almost no margin before boundary conditions.

**Scoping probe, 2026-11-18.** Linear algebra only, with no optimiser and no bound computed. It uses no OpenSees
build. numpy 2.3.5, scipy 1.16.3, gmsh 4.15.0. Appendix A gives the method. No output file was kept, and the
table is the record.

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

The probe is small, uses fixed walls only, and selects its footprint by node position. It motivates the gate and
does not replace it. **The thin slab passes, so a thin-slab pass says nothing about true 3D.** A structured mesh can
also hide the problem.

The tracked home for the gate's driver and records is `Ladruno_files/testbed/limit_analysis/`. This branch contains
no driver or result there, and this ADR reports none.

**Gate definition [proposed detail].**

| Item | Requirement |
|---|---|
| Meshes | non-structured tetrahedra, committed or reproducibly generated, at **three or more sizes** |
| Walls | every size and angle with **fixed** walls and with **roller** walls |
| Footing | a **conformal** footprint: its boundary lies on mesh edges and its facets are labelled |
| Criteria | `phi = 0` (von Mises) and Drucker-Prager matched in plane strain at `phi = 10, 20, 30` degrees. These are the draft contract's placeholders until frozen |
| Formulations | the control and each candidate, reported separately |
| Step 1, algebra | rank and nullity of the active linear constraints, and the scaled feasibility residual, with the statuses of the draft contract's section 11 |
| Step 2, quality | a **certified** Upper at every size, and a certified gap to the certified Lower of the same problem no larger than the frozen maximum at the finest size. The gap trend is recorded |

Feasibility alone does not pass. A tolerance may not be relaxed to turn `formulation_infeasible` into a pass, and
`AlmostSolved` is not success. If the control fails and a candidate passes, the gate passes and the control's failure
is recorded. If no candidate passes, K2 applies.

### 10.2 A2: a jump must be proven over the whole facet [owner]

A P2 jump on a flat facet is quadratic, and so is the difference between P2 soil velocity and a rigid footing. Even all
six nodal values are not enough: with vertex values `1` and midside values `0`, a P2 function on a triangle is `-1/3` at
the centroid.

The proposed proof uses the six Bernstein coefficients:

- `c_a = j_a` at the vertices, and `c_ab = 2 j_ab - (j_a + j_b) / 2` for the edge midpoints.
- The jump is a convex combination of the `c`, because the quadratic Bernstein polynomials are non-negative and sum to
  one. If every `c` lies in the convex jump cone, the whole jump does.
- Each polynomial integrates to `A / 6`, so `(A / 6) sum_c d_j(c)` over-estimates the jump dissipation.

The condition is sufficient, not necessary. Another whole-facet proof is acceptable.

### 10.3 A3: `domain_claim` belongs to each bound [owner]

See §12. A bracket claims only the domain both of its bounds support.

### 10.4 A4: certificate repair must be admissible [owner]

See §11. Outside the flow cone `d` is infinite, so a finite dissipation charge cannot pay for a violated cone. A Lower
interior field must itself satisfy equilibrium and the fixed loads. If no admissible construction exists, the result stays
an estimate.

---

## 11. Certificate and terminology

### 11.1 What "certified" means [proposed]

A value is a **certified Lower** or **certified Upper** for one problem hash, formulation id, query and `domain_claim`
when the independent checker, using only the problem file and saved fields, has built a final field and shown all of:

1. the identity, geometry and straightness checks pass;
2. every linear equality of the formulation holds on the final field within `tau_eq` (§11.3);
3. every inequality and cone condition required by §8, §9 and §13 holds on the final field, evaluated by a route that
   does not use the solver's cone map;
4. the value is computed from that final field by the formulas of §2, with exact work and resultant integrals;
5. the `domain_claim` is supported by the rules of §12.

This is a double-precision certificate with a reported floor. It is not interval arithmetic and not a machine proof. It
says that the field exhibited meets the stated conditions to rounding level. It says nothing about whether the
idealisation represents a real material.

### 11.2 Procedure [proposed; repair rule owner-accepted]

1. **Identity.** Check hashes, query, formulation, labels, volumes and straightness.
2. **Rebuild.** Build the operators through the checker's own code. Record the incoming residuals. A solver status is an
   input, never a verdict.
3. **Exact affine projection, first.** Project onto the formulation's linear constraints, using a stable least-squares
   solve that handles dependent rows.
   - Lower: equilibrium, facet traction continuity, boundary tractions, and resultants with the load relation.
   - Upper: supports, rigid footing motion, query normalisation, isochory where it applies, and normal jumps where they
     are zero.
   - An inconsistent constraint set produces no certificate. That is the expected outcome for an isochoric control that
     has locked.
4. **Inequality repair, second.** Use only a stated admissible construction that stays in the affine set.
   - Lower: a convex combination with a strictly admissible field that satisfies the same equalities at its own multiplier.
     Zero stress qualifies only with no fixed load and positive strength everywhere.
   - Upper: a convex combination with a strictly admissible velocity field that satisfies the same linear conditions. For
     von Mises the flow condition is an equality, so step 3 handles it.
   - A dissipation charge is never a repair.
5. **Rerun every check** of §11.1 on the repaired field, from the start. Any failure leaves an estimate.
6. **Value and floor.** Compute the value from the final field, and report the remaining floating-point floor.
7. **Pair check.** For a pair, test `Lower <= Upper`, with the floors reported.

### 11.3 Thresholds [proposed]

Equality acceptance is at rounding level and aware of conditioning. It is never a solver tolerance.

- With `a_i` a constraint row, the scaled residual is
  `r_i = abs(a_i z - b_i) / (abs(a_i) . abs(z) + abs(b_i))`.
- With unit roundoff `u` and `kappa` an estimate of the condition of the row-scaled constraint operator on its row space,
  from the factorisation used in step 3, the threshold is `tau_eq = C_eq u max(1, kappa)`.
- If `tau_eq` exceeds the named cap `TOL_EQ_CAP`, the constraints are too ill-conditioned to certify. The result stays an
  estimate, and the cap is not raised.
- `C_eq`, `TOL_EQ_CAP` and the inequality margins are named constants, justified by scale and calibrated on exact cases
  before first use. They may not be changed to rescue a run.

An intake limit may decide whether projection is worth attempting. It does not certify anything. The draft contract's fixed
`TOL_EQ_FINAL` must be replaced by this rule, or justified against it, before the freeze.

### 11.4 Terms

| Term | Meaning |
|---|---|
| **Certified bound** | met §11.1 |
| **Estimate** | everything else, including a solver's claimed optimum and an uncertified `AlmostSolved`; its reason is stated |
| **Inadmissible** | an estimate that cannot satisfy a gate |
| **Bracket** | a certified Lower and a certified Upper for the same stated problem |
| **Gap** | `(Upper - Lower) / ((Upper + Lower) / 2)` |

Write `[Lower, Upper]` only for a certified bracket. Otherwise write `estimate: L=..., U=...`. "Rigorous" means rigorous for
the stated idealisation, up to the reported floating-point floor. Every result reports its `domain_claim`, criterion,
associated-flow status and interface law.

A gate cannot pass on an estimate. A certified Lower above a certified Upper for one problem stops the program (K3).

### 11.5 Negative fixtures

Each must be rejected, with its reason code.

| Applies to | Fixture |
|---|---|
| Both | wrong load normalisation; swapped Lower and Upper labels; curved tetrahedron; mismatched problem hash; unsupported `domain_claim`; equality residual just above the threshold |
| Lower | one perturbed stress component; stress scaled by 1.01; facet traction jump; wrong gravity in equilibrium |
| Upper | `tr(D) != 0` at one vertex of an isochoric field |
| Upper | Drucker-Prager flow-cone violation at one vertex |
| Upper | footprint node off the rigid motion |
| Upper | non-zero support velocity |
| Upper | non-zero wall velocity under a zero-extension claim |
| Upper | non-zero normal jump at one facet node under isochoric flow |
| Upper | jump admissible at all six facet nodes but not across the facet, using the §10.2 example |
| Upper | wrong gravity work |
| Upper | wrong dissipation quadrature |
| Upper | non-positive reference work rate |
| Upper | a finite dissipation charge offered in place of restoring a violated cone |

A test reads the run record and compares it with the requested case, so a plausible number from the wrong case fails.

---

## 12. Domain claims

**[owner]** Each certified bound carries its own `domain_claim`. A bracket claims only the domain both bounds support.
**[proposed]** The values and rules below; O5 for the order of work.

| Value | Meaning |
|---|---|
| `box_stated` | a bound for the stated box and its exact wall conditions |
| `unverified` | extension to a wider domain is not established; for that wider domain the value is an estimate |
| `halfspace_by_zero_extension` | the admissible field extends by zero, under the conditions below |
| `halfspace_by_extension_elements` | a Lower field has a verified admissible extension; reserved until such data exists |

- **Upper by zero extension:** every truncation face is fixed and no load acts outside the mesh.
- **Lower by zero extension:** every truncation face is traction-free, all body forces are zero and no load acts outside
  the mesh.
- Roller walls can permit tangential motion, so they do not support the Upper claim.
- A fixed wall supports no Lower extension.
- With self-weight, a traction-free truncated Lower field does not extend by zero.

### 12.1 Fixed-wall box lemma

Let box `B` lie inside half-space or half-plane `H`, with the same footing, strength and loads, and no load outside `B`.
Let `lambda_B` and `lambda_H` be their collapse multipliers.

1. **Every fixed-wall mechanism of `B` extends by zero to `H`** with the same dissipation and work. So
   `lambda_H <= lambda_B`, and a certified Upper for the fixed-wall box is also an Upper for `H`.
2. **If `B` contains the support of a kinematically admissible mechanism of `H` with value `lambda_ref`**, that mechanism is
   admissible in `B`. So `lambda_B <= lambda_ref`.
3. **If `lambda_ref = lambda_H` is exact**, then `lambda_B = lambda_H`. The box bracket must contain the reference value.

Result 3 is a comparison for that box. It does not change the Lower bound's claim, which stays `box_stated`.

### 12.2 Extent conditions for the strip comparison

For footing half-width `b`, with `beta = pi/4 + phi/2`, `r0 = b / cos(beta)` and `r1 = r0 exp((pi/2) tan(phi))`, the Prandtl
mechanism reaches `x = b + 2 r1 cos(pi/4 - phi/2)` from the centre and depth at most `r1`. For `phi = 0` these are `3b` and
`sqrt(2) b`.

The `(2 + pi) cu` comparison requires:

- fixed walls and base;
- a box that strictly contains that extent, with the margin frozen in the case;
- rollers on both slab faces, so that the slab reproduces plane strain;
- a weightless problem, with the load definition recorded.

With roller walls, result 2 survives and result 1 does not, so only the Lower inequality is available.

Version 1 produces box-labelled bounds, then an extent study. Extension elements come later.

---

## 13. Query design [proposed]

Let `F` be the six-component footing resultant about a declared reference point. Let `K_F` be the set of resultants the soil
can carry with dead loads fixed. `K_F` is convex, since it is the linear image of a convex set.

**Rotation scaling.** With a declared length `ell > 0`, reported coordinates are `F_hat = (Fxyz, Mxyz / ell)` and
`q_hat = (u, ell w)`, where `u` is footing velocity and `w` angular velocity. Then `F_hat . q_hat = F . u + M . w`, so power is
preserved. A support direction is given in these coordinates, and changing `ell` changes what "direction" means but not the
set `K_F`.

| Query | Lower | Upper |
|---|---|---|
| `support(d_hat)` | maximise `d_hat . F_hat`; gives a point inside `K_F` | prescribe `q_hat = d_hat`; `Diss - W0` gives the half-space `d_hat . F_hat <= Diss - W0` containing `K_F` |
| `ray(F0, F1)` | maximise `lambda` with `F = F0 + lambda F1` | `q` free with `F1 . q = 1`; `lambda <= Diss - W0 - F0 . q` |
| fixed components | a ray, with the fixed components in `F0` | same |

Report `conv(lower points) subset K_F subset intersection(upper half-spaces)`. Every Upper solution also yields a half-space
with its own `q`.

Assembled operators are reused across queries on one mesh. A direction-specific adaptive mesh is allowed, because any valid
mesh gives a valid bound. Bounds combine only for exactly the same stated problem: geometry, material, interface, reference
point, `ell` and `domain_claim`.

**Interface laws.** Take `n` out of the soil, `t = sigma n` on the soil, `t_n = t . n` (tension positive), `t_t` the tangential
part, and `j = v_footing - v_soil`, with normal part `j_n` (opening positive) and tangential part `j_t`.

| Law | Traction set | Upper condition | Interface dissipation |
|---|---|---|---|
| `bonded` | no interface limit | `j = 0` | none |
| `no_tension` | `t_n <= 0`; **no interface limit on shear**, so shear is limited only by soil yield | `j_t = 0`, `j_n >= 0`: opening without slip | zero |
| `no_tension_smooth` | `t_n <= 0`, `t_t = 0` | `j_n >= 0`, `j_t` free | zero |
| `coulomb` (`mu > 0`, adhesion `c_i >= 0`) | `norm(t_t) + mu t_n <= c_i` | `j_n >= mu norm(j_t)` | `c_i j_n / mu` |

The Lower side imposes the traction set at three facet vertices. The Upper side uses the whole-facet proof of §10.2. The draft
contract admits only `bonded` in version 0.

The Coulomb law is associated, so sliding dilates. For real non-dilatant friction, an associated interface result is not a bound
on the physical problem.

---

## 14. Open-solver policy

### 14.1 Rule

**[owner]** No licence requests. Nothing that needs a licence goes into the fork. FEniCSx may be used outside the fork, with open
solvers.

**[proposed] reading of that rule for this feature:**

- No MOSEK code, adapter, licence check or test requirement enters the fork, nor any equivalent commercial solver. No required
  gate needs one.
- The fork's existing host toolchain, including Intel oneAPI and MKL where its build already uses them, is treated as outside the
  rule. This carve-out is an interpretation that the owner has not separately ruled on. New limit-analysis code does not depend
  on it, and no limit-analysis gate requires MKL.

### 14.2 Mechanism [proposed]

- **D5.** The default build links no solver and needs no Rust.
- **D6.** A run-time adapter with a small versioned C ABI connects a solver. `getLibraryFunction` is the precedent for loading it.
- **D7.** `limitAnalysis export` always works where the engine is compiled, so a run can be solved outside the process and checked.
- How users obtain an adapter is O15.

### 14.3 Candidates

Licences come from the GitHub licence API on 2026-11-18 unless stated. Re-check before adding or shipping anything. This table was
built from licence metadata and project documentation; no GPL or LGPL source was opened.

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
| MUMPS | CeCILL-C (`OTHER/MUMPS/include/dmumps_c.h`) | sparse linear solver for a possible tailored method; see §14.4 for availability |

Clarabel build features `pardiso-panua` (its README says a separate licence is required) and `sdp-mkl` are out of policy for adapters.
The Python wheel in the development environment, version 0.11.1, solved a trivial problem with `qdldl` and `faer` and reported `mkl`
and `panua` as not available.

No claim is made here about which solver reaches certification at scale. That is 237d's question.

### 14.4 A tailored solver

Consider a tailored interior-point method on MUMPS only when all three hold:

1. after P0's workstation test, open general solvers cannot produce a certified result on the reference 3D footing within one cluster
   node's time and memory;
2. at least 80% of solve time is in the KKT factorisation and solve;
3. a standalone test on the same reduced system is at least three times faster.

**MUMPS is not available in every build.** In `CMakeLists.txt` it is configured inside the `MPI_FOUND` branch and only when `MUMPS_DIR`
is defined; otherwise the flag is `_NOMUMPS`. `MUMPS_LIBRARIES` is linked into `OpenSeesSP`, `OpenSeesMP` and `OpenSeesPyMP`, not into
the serial `OpenSees` and `OpenSeesPy` targets. Version 1 of this engine refuses the parallel targets. A tailored method would therefore
need its own link decision, outside the default build.

---

## 15. FEniCSx stays outside the fork

**[owner]** FEniCSx may serve as an independent test reference with open solvers, outside the fork. The lane lives in a separate private
repository.

**[proposed] working rules for the lane:**

- 237e assembles the same frozen cases independently and solves them with open solvers. Agreement is judged on the same exported problem
  class and exact cases.
- Licences are checked before use. DOLFINx and UFL report LGPL-3.0 and Basix reports MIT. FFCx reports no SPDX match, so its licence file
  must be read first.
- The lane's author may read FEniCS code. Fork workers do not read it for this feature, and no code passes from that lane into the fork.
- No fork gate requires the lane. A disagreement is resolved or recorded before I1 if the lane has started.

---

## 16. Verification and warrant [proposed]

The harness follows `Ladruno_implementation/testbed/00_canonical_testbed.md` and ADR-87.

| Level | Purpose |
|---|---|
| L0 | operators: constant and linear stress, rigid motion, footing resultants, cone membership, exact work integrals |
| L0b | from P1: Python and C++ programs match entry by entry |
| L1 | checker: good fixture passes, §11.5 fixtures fail |
| L2 | identical exported programs through independent open solver paths |
| L2e | external FEniCSx assembly, when active, outside the fork |
| L3 | exact reductions, the gravity cases of §16.1 and property gates |
| L4 | true 3D: the §10.1 gate, property gates, later Tresca references |
| L5 | associated OpenSees models compared with certified brackets |
| L6 | adaptive: gap reduction, re-resolved declarations, budget and size-field fidelity |

Property gates, where a case supports them: certified Lower not above certified Upper; nested uniform refinement does not lower the Lower or
raise the Upper; rigid rotation of the problem leaves both unchanged; doubling cohesion on a weightless problem doubles both; a valid half model
matches the full model. Reference gates are one-sided, plus a frozen maximum gap, so a bound cannot pass from the wrong side.

Reference cases enter from a checked source, with date, units, load definition and tolerance. The first reference is the plane-strain strip with
`(2 + pi) cu` [Pra21], as a thin slab with matched von Mises strength, under §12.2.

### 16.1 Gravity

**Level-ground von Mises footing cases do not detect gravity errors.** Take a level free surface at `z = 0`, a footing on it, gravity along `z`,
and zero normal velocity on the other boundaries. Then gravity work is `gamma int v_z dV = gamma int grad(z) . v dV`. Integrating by parts, with
`div v = 0` and zero normal jumps, leaves only boundary terms, and each is zero because `z = 0` or `v . n = 0` there. On the Lower side, a hydrostatic
stress balances gravity with zero deviator. Zeroing gravity therefore leaves both bounds unchanged, and a mutant survives such a case.

P0 therefore requires cases where gravity changes the answer:

| Case | Setup | Closed-form statement | Assertions on certified values |
|---|---|---|---|
| **Shear layer (exact)** | layer of thickness `H` on a fixed base or stronger substrate, von Mises `k`, dead body force with shear component `b_x`, top shear load `lambda`, side tractions taken from the exact field | Lower field `sigma_xz = lambda + b_x (H - z)` and rigid translation with a base jump give `lambda* = k - b_x H` | with a base jump: bracket contains `k - b_x H`, within the frozen gap. Without one: `L >= k - b_x H` and `U <= k - b_x (H - h / 2)`, using a first layer of height `h` |
| **Uniaxial column (sandwich)** | prism, free sides, roller base, top compressive load `lambda`, unit weight `gamma`, uniaxial strength `S` | the uniaxial stress field gives `S - gamma H`, and the homogeneous isochoric or dilating field gives `S - gamma H / 2`. So `S - gamma H <= lambda* <= S - gamma H / 2`. Neither value is claimed to equal `lambda*` | `L >= S - gamma H`, `L <= S - gamma H / 2`, `U >= S - gamma H`, `U <= S - gamma H / 2` |

**These are the mutant-killing cases.** Zeroing gravity moves the certified values to `k` and `S`, above the upper limits in the last column, so both a Lower and an Upper gravity mutant fail. A sign flip fails the same way. The fields are in the P2 space, so the column and the unjumped shear layer also work for the continuous control. The column statement for Drucker-Prager uses its compressive uniaxial strength.

The draft contract's `vm_body_force_shear_block` uses a body-force multiplier. It is a useful exact Lower case and does not replace these dead-load cases. 237a freezes the final cases.

Self-weight references for frictional soils enter only after their primary sources are checked.

### 16.2 Where tests run

**O1 [proposed].** The frozen contract, Python reference and checker live in `tests/`. Solver studies live in `Ladruno_files/testbed/limit_analysis/`.

| Lane | Runs | Solver |
|---|---|---|
| Zone-A or `static-gates` | contract, operators, checker and negative fixtures; no OpenSees build is needed for the P0 Python tests | none |
| **`limit-open-solver` (new, proposed required check)** | L2 and the reduced L3 set on exported programs | a pinned open solver, Clarabel first, launched as a **child process** that reads the exported program and writes a solution file |
| Testbed studies | 237d's size and scale runs | open solvers; long runs launched by the orchestrator |

In `limit-open-solver` a missing or failing solver **fails the job**. It is never a skip, and `pytest.importorskip` is not used. The job is pure Python, like the
MPM contract step in `static-gates`. Making it a required check is a branch-protection change for the owner.

### 16.3 Warrant for the C++ family

At 237f: non-self-referential oracles; a `LIMIT` mutation family through `CMakeLists.txt::LADRUNO_MUTATE_FAMILY`, `Ladruno_scripts/mutation_gate.py::FAMILIES` and
`SRC/Ladruno_mutation.h`; a guide; both shipped interpreters or a recorded absence; ledger fragments; a banner line; and a manifest row added deliberately, with
`tag_symbol: null` if there is no class tag (W3). Mutants should fail when facet continuity is dropped, gravity is zeroed, the resultant map changes sign, or a cone
map is weakened.

---

## 17. Dependencies

### 17.1 apeGmsh ADR 0136, Level A [proposed]

Limit analysis is stateless, so each new mesh reruns the whole analysis. That is Level A of ADR 0136. **ADR 0136 currently says 2D first**: its review log records "2D
first" as confirmed, and its named S3 benchmarks are 2D. This program requires 3D, so the draft needs the requests below. They go through that ADR's process and T; this ADR
decides none of them.

| Request | Draft text at `668b99a44` | Needed | Earlier label |
|---|---|---|---|
| **R1, estimator decoupling** | D1's driver "estimates"; S1 and the Kirsch and L-shape gates of S3 rest on ZZ. INV-6 refuses "a quadratic element in version 1" | the driver takes an estimator callable that may return a plain element array, so tet10 limit fields never need ZZ | decoupling, A1 |
| **R2, non-`Results` run product** | D1: `run(fem, k)` "returns a `Results`" | `run` may return an opaque study product, here the certified bounds and fields, which the driver passes to the estimator and stop rule | new |
| **R3, stop predicate** | the loop targets an error level | a study-specific stop: certified gap, budget, stagnation, iteration limit | A2 |
| **R4, 3D size-field gate** | S2 gates name a prescribed size, the model-data refusal and a tensor path | a scalar 3D list-data gate with budget and fidelity | A4 |
| **R5, straight tet10 output** | no statement about straight second-order elements | `build` can produce straight tet10 meshes, or refuse, with straightness verified after each rebuild | new |
| **R6, later options** | D2: "It refines only"; `SizeField` is built from the old mesh | coarsening opt-in and a portable size field, when measurements call for them | A3, A5 |

**Adaptive gate I3:** 237g, C2, C3 and R1 to R5. R6 joins only if measurements show it is needed.

Kept out of ADR 0136: bound, certificate, pair and solver concepts; a fork dependency in S1 to S3; two meshes per variant row; in-session remeshing; footing-specific rules;
3D anisotropy in the first slice.

Adaptive loop (237h): build a fresh session, write the problem file, run both bounds, certify, compute element values from certified fields, test the stop, build the next
size field. Labels, supports, loads and footprints are verified after each rebuild. A failed certification is logged and not counted. The best certified Lower and lowest
certified Upper are kept only across meshes of exactly the same stated problem. Every indicator is unverified until tested for this formulation; published adaptive schemes
[Lya05, Cir08, Mun08] are starting points only.

### 17.2 apeGmsh bridge [proposed]

A separate apeGmsh ADR, ratified by its maintainer, covers the bridge. ADR 0134 is the precedent: emission only, a sidecar written together with the deck or not at all, new
census families by exact entry, primitives only after the fork feature is merged and in a release named by `FORK_PIN`, and a member Warranted only when its own `live-fork` test
passes (ADR 0121 D2 and D3). All apeGmsh work goes through its T coordinator. Whether to use new families or a token under `analysis` is O6.

### 17.3 Fork ADR 236 [proposed]

**Level A does not depend on ADR 236.** The engine uses no material history, no `getHistory` or `setHistory`, and no commit primitive. No limit-analysis phase waits for it.

The restart lane is coordinated separately: ADR 236 S1 to S5 and ADR 0136 S4 to S6. ADR 232's FE-to-MPM handoff waits on the same primitive, so ADR 236 S1 must be scheduled with
the MPM program, without duplicate work. The open questions in ADR 236 §6 stay with that ADR.

### 17.4 Shared resources [proposed]

These are touched by more than one program, so changes to them merge one at a time, with affected tests rerun on the merged result.

| Resource | State at `cfe478254` | Who needs it | Rule |
|---|---|---|---|
| `.ladruno` `FORMAT_MINOR_VERSION` and `ladruno_format.py::validate` | `LadrunoRecorder::initialize` writes minor 1; the MPM `createFile` writes minor 2; `validate` is the executable schema | 237f writer; ADR 236 S4 (D7), which plans a minor bump and a `validate()` arm | the next minor number is assigned when the PR is ready to merge, and recorded on #1118. Each change updates the schema text, `validate()` and fixtures together |
| `SRC/material/nD/NDMaterial.cpp::OPS_clearAllNDMaterial` | carries the fork's process-wide resets | 237f, if the engine keeps process-wide state; ADR 236 INV-5 | one reset call per feature, order-independent, with a `vanilla/` ledger row |
| `SRC/interpreter/LadrunoCommandTable.h` | 36 rows, no conditional rows | `limitAnalysis`; ADR 236 D2 (`ladrunoGetHistory`, `ladrunoSetHistory`, `ladrunoCommitState`, plus a flag on the existing `ladrunoSetNodeTrial`) | one row per command, with both columns matching what every build registers (§4.3) |
| apeGmsh `FORK_PIN` | names one fork release; ADR 0134 requires its keys to move together | L4; D6 (ADR 0136 S4); other fork-backed bridge work | T sequences every move. One move may carry several features, and each member still needs its own `live-fork` pass |

---

## 18. Phases and kill criteria [proposed]

### 18.1 Phases

| Phase | WP | Delivers | Exit |
|---|---|---|---|
| P0 | 237a to 237e | contract, Python reference, checker, solver study, external reference | §18.2 |
| P1 | 237f | C++ kernel, assembler, command, export, writer, adapter | L0b identity for every P0-admitted formulation; checker certifies engine output; mutation lane fails on operator faults; a build without the engine differs only by the command's refusing entry point, and existing commands are unchanged; shipped interpreters handled or absence recorded |
| P2 | L1 to L4, apeGmsh | bridge and reader | emitter parity; stock refusal; sidecar atomicity; `Results.limit` fixtures; `live-fork` proof on the pin |
| P3 | 237g | V-H-M and 3D queries | strip V-H and V-M certified brackets; sliding correct for the stated interface; query reuse measured; one six-component 3D demonstration labelled by criterion and domain |
| P4 | 237h | adaptive 3D | I3; certified gap falls at fixed budget on references; every counted iteration passes the checker and re-resolution; best bounds tracked correctly; 3D size-field fidelity proven |
| P5 | 237i | semidefinite criteria and richer interfaces | cone maps agree with principal-stress evaluation; rough circular Tresca bracketed against its checked reference [ES60]; Mohr-Coulomb reduces correctly to plane-strain references; open-solver cost measured |
| P6 | 237j | performance and advanced meshing | each item shows measured gain over P4 |

### 18.2 P0 exit

P0 is complete only when all of these are demonstrated:

1. The orchestrator has frozen 237a, including the admitted Upper formulation or proven alternative.
2. L0, L1, L2 and L3 pass. Python-to-C++ identity is not part of P0.
3. The certified fixed-wall thin-slab bracket contains `(2 + pi) cu`, under §12.1 and §12.2.
4. Each reference bracket meets its **frozen maximum gap**. The value is recorded in the case before the first solve of that case and is not changed after seeing a result.
5. The true-3D gate of §10.1 passes, at both algebra and certified-quality steps.
6. The gravity cases of §16.1 pass, and the gravity mutants fail.
7. Every §11.5 fixture is rejected.
8. Property gates pass for each P0 case where they apply.
9. The solver study states which open path certifies at what size, and K1 is decided.
10. Exported programs replay through a second solver path, and the `limit-open-solver` lane is green.
11. If 237e has started, any mismatch with FEniCSx is resolved or recorded before I1.

### 18.3 Kill criteria

| Kill | Trigger | Action |
|---|---|---|
| K1 | no open path certifies 3D Lower and Upper at about 20,000 tetrahedra, within one hour and one workstation's memory, for **each of `phi = 0, 10, 20, 30` degrees** (`phi = 0` is von Mises; the others are plane-strain-matched Drucker-Prager). The size, time and angles are provisional until frozen | stop before C++; reformulate, rescale, prototype a tailored open solver, or reduce scope. A reduced scope that drops `phi = 0` is an owner decision |
| K2 | the thin-slab strip fails after operator tests pass, or no candidate passes the true-3D gate | stop and rethink the formulation |
| K3 | certified Lower above certified Upper for one problem | stop the program |
| K4 | solver integration disturbs the default build | keep or repair the run-time adapter before going on |
| K5 | after adaptivity, the certified gap cannot separate the result it must judge | label the 3D result limited evidence and ask the owner |

A failed gap gate keeps P0 open but is not a kill rule by itself. Thresholds are not tuned to reach a reference, and neither is an element, tolerance or mesh rule.

---

## 19. Risks, decisions, sources

### 19.1 Risks

| Risk | Consequence | Mitigation |
|---|---|---|
| No Upper formulation passes the true-3D gate | no usable Upper for low-friction criteria on general meshes | control plus candidates in P0, K2 |
| Jump cones are expensive | the program grows and solvers reach their limits sooner | 237d measures; selective jumps; §14.4 |
| Interior-point accuracy | solutions may not repair into certificates | scaling, formulation work, checker as the only judge, K1 |
| No strictly admissible repair field | valid-looking results stay estimates | a stated construction per case; the result stays an estimate otherwise |
| 3D program size | memory and time limit the mesh, so gaps stay wide | 237d; adaptivity; §14.4 |
| Wide fixed-mesh gaps | a bracket is true but not useful | frozen maximum gaps; P4; K5 |
| A gravity error passes level-ground cases | wrong self-weight results look validated | §16.1 |
| Reference and engine share an error | agreement proves nothing | separate authors, frozen contract, independent checker, external lane |
| Misread as a soil-model result | an associated bound is applied to a non-associated problem | required labels (§11.4), non-goals |
| Domain misstatement | a box result is cited as a half-space result | per-bound claims (§12) |
| Shared-resource collision | format, wipe, command or pin changes conflict | §17.4 |
| Reader breakage | `.ladruno` additions break existing readers | D3, fixtures before relying on readers |
| Licence drift | a dependency changes terms | re-check before adding or shipping (§14.3) |
| ADR 0136 stays 2D-first | the adaptive phase slips | R1 to R5; P0 to P3 do not depend on it |

### 19.2 Decisions

**Closed by the owner, 2026-11-18.** Former O16: amendments A1 to A4 accepted. Former O11: the FEniCSx lane lives in a separate private repository.

**Working defaults, proposed, not decided.**

| ID | Default | Confirm before |
|---|---|---|
| O1 | P0 contract, reference and checker in `tests/`; solver studies in `Ladruno_files/testbed/limit_analysis/` | first push |
| O2 | Fork ADR 237; apeGmsh limit ADR through T; amending ADR 0136 | the dependent step |
| O3 | Separate problem file for the engine | P1 |
| O4 | Von Mises and Drucker-Prager first; Tresca and Mohr-Coulomb in P5 | P1 |
| O5 | Box-labelled bounds first, extent study, extension elements later | any half-space claim |

**Open.**

| ID | Question | Recommendation | Needed by |
|---|---|---|---|
| O6 | New apeGmsh families or a token under `analysis` | new families | L1 |
| O7 | Amend apeGmsh#3052 before acceptance, or follow up | before, if T agrees | C0 |
| O8 | Is an ADR-87 D4 out-of-family verdict required | decide with this ADR | 237f merge |
| O9 | Classic Tcl wired or recorded absent | wire if inexpensive | 237f |
| O10 | A separate licence for the OpenSees-free kernel | keep the option open; seek advice first | external release |
| O12 | Do fork merges stay with the owner, or are they delegated | owner unless the program issue says otherwise | first ready-flip |
| O13 | ADR 236 §6 questions | decide in D0; coordinate S1 with MPM | D1 |
| O14 | ADR 0136 questions: Level A first, driver location, Bezier gate | Level A first; driver package; measure Bezier | C0 |
| O15 | How run-time solver adapters are distributed | separate from the default fork build; decide before release | 237f release |
| O17 | Make `limit-open-solver` a required check | yes, once it is stable | I1 |
| O18 | Confirm the reading of §14.1 on the existing MKL toolchain, and the FEniCS reading rules of §15 | as written | 237f; 237e |

### 19.3 Unverified or partly verified

- **Literature.** Bibliographic data below were checked against Crossref or arXiv on 2026-11-18. Of the papers, only [BH21] was read for content, in its arXiv version. Attributions to the
  others are from memory and are marked `(m)` in the list.
- `T` near `6N` for production tetrahedral meshes (§10.1): recalled.
- Applicability of [Zha11]: title and metadata only.
- The jump formulation's accuracy and cost: not measured.
- Solver reach at scale and `AlmostSolved` behaviour: no reproducible record in this repository, so no numbers are given.
- The rough circular reference value is not stated here. It is entered in P5 after checking [ES60].
- Square-footing references [Sal04, Lya07] are candidates only; values and conditions are unchecked.
- A self-weight strip reference (Martin's `N_gamma`, IACMAG 2005) was not found in Crossref and is unverified.
- Licences of FFCx and Mmg (GitHub reports no SPDX match) and gmsh, CVXOPT and Hypatia: not verified here.
- Reader behaviour for new `.ladruno` element groups and for `MODEL/LIMIT_ANALYSIS`: untested.

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

## Appendix A. Method of the §10.1 scoping probe

This records the method of the table in §10.1. It is not a reference implementation and it is not the gate. The gate's tracked driver belongs in
`Ladruno_files/testbed/limit_analysis/`.

1. **Unknowns.** Three velocity components at every vertex and edge midpoint of straight tetrahedra.
2. **Constraint rows.** For each tetrahedron and each of its four vertices, one row for `tr(D) = 0`, from the gradients of the P2 shape functions at that vertex.
3. **Unit-cube footing.** Base and four sides fully fixed. Footprint: nodes on `z = 1` with `0.25 <= x, y <= 0.75`, prescribed to `v = (0, 0, -1)`. The rest of the top is free. This footprint is
   selected by node position and is not conformal, which is one reason the probe does not satisfy the gate.
4. **Slab.** `0 <= x, z <= 1`, thickness `h` in `y` (`1/n` for the structured mesh). `v_y = 0` on both `y` faces; `x = 0`, `x = 1` and `z = 0` fully fixed; footprint nodes on `z = 1` with
   `0.375 <= x <= 0.625` prescribed to `v = (0, 0, -1)`.
5. **Meshes.** Structured: each cube cell split into six tetrahedra around its main diagonal. Unstructured: gmsh 4.15.0 default 3D meshing of an OpenCASCADE box, with minimum and maximum size both `h`.
6. **Measures.** With `A` the free-column matrix and `b` the prescribed-motion right-hand side: rank from singular values above `1e-10` of the largest; nullity = free DOF minus rank; residual =
   `norm(A x - b) / norm(b)`, with `x` the least-squares solution at `rcond = 1e-10`.

---

## Review log

- **2026-11-18, drafting.** Written from the owner's decisions and the working design, then checked against the fork at `cfe478254`. The check produced walls W1 to W6 and amendments A1 to A4.
- **2026-11-18, review fix round.** Program issue #1118 recorded; the owner's acceptance of A1 to A4 and the FEniCSx location recorded. Continuous P2 became a control, and P0 now owns a jump-based
  formulation or a proven alternative. The certificate, true-3D gate, gravity cases, P0 exit, build rule, required solver lane, shared-resource table and ADR 0136 requests R1 to R6 were rewritten.

## Implementation log

Nothing built.
