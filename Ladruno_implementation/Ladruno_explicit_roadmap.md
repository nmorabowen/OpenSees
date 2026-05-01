---
title: Ladruno — Explicit & Multi-Body Research Roadmap
project: Ladruno
status: living
priority: high
owner: nmora
tags:
  - index
  - roadmap
  - planning
  - explicit
  - multi-body
  - ls-dyna
created: 2026-05-01
updated: 2026-05-01
---

# Ladruno — Explicit & Multi-Body Research Roadmap

> [!summary]
> Roadmap for extending OpenSees with capabilities currently missing or weak: **explicit dynamics performance**, **rigid bodies + joints**, **contact**, **meshfree (SPH/EFG)**, **modular damage**. Direction is informed by an architectural comparison with LS-DYNA. Each item below will graduate to its own per-feature plan in this folder using [[_template]] when work begins.

> [!info] Related skills and references
> - LS-DYNA learning companion: `~/.claude/skills/ls-dyna/` — indexed against the four LS-DYNA manuals in `Dropbox/UANDES EC/Lit Review/Explicit/LS-DYNA/`. Use this for any LS-DYNA Theory/Vol I/Vol II citation.
> - OpenSees source / patterns: `anthropic-skills:opensees-expert`.
> - Continuum mechanics & FEM theory: `anthropic-skills:fem-mechanics-expert`.
> - Compilation context: [[../Ladruno_internal/01_compilation_journal|compilation journal]].

---

## 1. Why this roadmap exists

### 1.1 Why explicit — the class of problems we want to reach

This roadmap is not "let's add explicit because LS-DYNA has it." It's **"let's open the class of seismic / geotechnical / nonlinear FEM problems that implicit fundamentally can't reach in our domain."** The capability gap matters; LS-DYNA is the best architectural reference for how to close it.

#### Where implicit hits walls in our problem space

OpenSees is excellent at implicit nonlinear dynamics within the regime where Newton iteration converges. The walls show up where the regime breaks down:

- **Path-dependent constitutive models in their softening / state-change regime.** PM4Sand, SANISAND, PDMY (sand liquefaction); multi-surface plasticity. Once the tangent goes singular or non-positive-definite, Newton fails. Line search, arc-length, step cutting are workarounds — but they *bound the strain regime we can actually study*. Liquefaction triggering is at the edge of this regime; post-triggering is past it.

- **Concrete fracture and quasi-brittle damage** — the canonical implicit-killer in our domain. Smeared damage models (CDPM, CSCM, KCC, MAT_72R3, Winfrith) live in their softening regime where the tangent goes non-positive-definite and Newton diverges; arc-length helps for monotonic but breaks under cyclic seismic loading. Discrete cracking via cohesive / traction-separation surfaces (LS-DYNA `*MAT_138/184/185/186/240`) introduces **snap-back** — load-displacement curve with negative slope *and* negative load increment. Element erosion (concrete spalling, RC fragment ejection under blast / impact) is a discrete state change. Crack branching produces multiple competing localization paths and Newton picks one and gets stuck. **Explicit handles all of this natively**: no tangent needed, snap-back is just dynamics, erosion is a switch, branching is local. This is the path to KCC/RHT/CSCM-class concrete simulation in OpenSees — currently the practical ceiling is whatever stays in the hardening regime.

- **Buckling and postbuckling.** Bifurcation makes the tangent singular by definition. Limit-point snap-through is genuinely dynamic — static analysis misses the kinetic-energy phase. Imperfection-sensitive postbuckling is multi-modal and unstable. Specific cases that bite us: **steel braces in CBF / SCBF systems** undergoing repeated cyclic buckling-yielding (the implicit recovery from each buckling event is fragile under realistic ground motion); **RC column folding** under axial+lateral combined loading; **slender shear walls** with out-of-plane buckling under in-plane shear; **plate / shell wrinkling and local buckling** in cold-formed sections; **liquid storage tank shell buckling** under seismic loads (sloshing + shell buckling coupling, currently essentially impossible in OpenSees). **Explicit makes these direct**: no tangent operator needed, dynamic snap-through is just integration, multi-modal interaction emerges naturally from the dynamics.

- **Contact and gap formation.** Soil-foundation interface uplift, pile-soil gap opening, inter-story pounding between adjacent buildings, foundation rocking, soil-pile-cap separation. Implicit contact is notoriously hard to converge. Explicit penalty contact is **robust by construction** — it just steps through.

- **Wave propagation through nonlinear media.** DRM with realistic high-frequency content; SSI with rocking, sliding, base separation; blast and impulse loads. Implicit wastes effort iterating to convergence at every small step, but the *accuracy-required* step is already small (CFL-scale). The implicit Newton overhead is pure cost in this regime.

- **Discrete state changes.** Element erosion, fracture initiation, gap opening/closing, friction stick-slip, dewatering. Each transition kills implicit convergence; explicit doesn't notice.

- **Massive 3D continuum models.** SSI domains with $10^6$+ DOF where the implicit linear solve dominates cost. Explicit per-step cost is $O(N)$ with no global solve, and parallelizes cleanly — including to GPU.

#### Concrete capabilities this opens

In the seismic / geotechnical / nonlinear-FEM context this group works in (UANDES + collaborators):

1. **Realistic collapse and damage simulation** past the hardening regime — into softening, localization, fracture. Currently the practical limit of OpenSees nonlinear simulation. Direct relevance: collapse-margin assessment for performance-based earthquake engineering.
2. **Concrete fracture from initiation through propagation, branching, and erosion.** Smeared damage models (CDPM/CSCM/KCC-class) deep into their softening regime; cohesive surfaces with snap-back; element erosion for spalling and fragment ejection under blast / impact; multi-crack interaction and branching. Opens the KCC, RHT, CSCM concrete-modeling space — currently inaccessible in OpenSees.
3. **Buckling and postbuckling with full dynamic response.** Cyclic buckling of CBF braces under realistic ground motion; RC column folding; slender pier P-Delta past the limit point; plate / shell wrinkling and local-distortional-global interaction; liquid storage tank shell buckling under seismic + sloshing. Many of these problems today need special-purpose models with careful arc-length tuning; explicit makes them direct.
4. **Soil liquefaction *with* the discontinuities** — gap formation, sand boils, lateral spreading patterns. Tracks state changes the implicit PM4Sand pipeline can't follow.
5. **Pounding studies for dense urban retrofit** — Latin American context particularly relevant; multiple buildings + slab + soil contact in one model.
6. **Foundation uplift and rocking** with real contact mechanics, not modeled-as-stiffness-loss.
7. **Pile-soil interaction with realistic gap mechanics** — beyond the p-y spring abstraction.
8. **Tsunami debris impact, slope failure, post-earthquake landslide / tunnel response** — once SPH / EFG come online ([[#5.8]] / [[#5.9]]), problems with severe deformation and free-surface flow that no OpenSees mesh can handle today.
9. **Direct cross-validation against industry codes** — the same OpenSees constitutive models running in an explicit framework can be benchmarked head-to-head against LS-DYNA / Abaqus/Explicit using identical material parameters. This is a verification door that closed commercial codes don't open.

#### Strategic doors

- **An open-source explicit-dynamics platform serious enough for academic research.** None exists at OpenSees-level domain depth. LS-DYNA / Abaqus are closed. Code_Aster, FEAP, FEniCS lack the seismic-engineering vocabulary. OpenSees + explicit fills this gap.
- **Bridge to meshfree / DEM / particle methods.** Explicit is the natural place for SPH-FEM, DEM-FEM, and particle-based fluid coupling. Implicit can't do these cleanly. This is the long-term meshfree research direction ([[#5.8]], [[#5.9]]).
- **Industrial relevance** — pounding, liquefaction, blast, impact studies that today require commercial licenses.
- **Research compounding** — every explicit element / material / contact algorithm becomes a contribute-able, citable piece. The platform compounds.

#### Honest tradeoffs

> [!warning] Explicit is not a panacea
> The argument is "broaden the class of problems we can solve," not "replace implicit." Explicit has real costs:
> - **CFL-bounded timestep** scales with smallest element size and stiffest wave speed — fine meshes are expensive.
> - **Quasi-static / slow-loading problems** are inefficient in explicit unless you use mass scaling, which introduces inertia artifacts.
> - **Frequency-domain accuracy at low frequencies** can suffer when mass scaling is aggressive.
> - **Energy balance** must be checked carefully; spurious hourglass / contact / mass-scaling energy can grow without warning unless monitored.
> Implicit OpenSees stays the right tool for static pushover, modal analysis, low-frequency seismic, and well-behaved dynamic problems. The roadmap *adds* a capability, doesn't replace one.

### 1.2 What LS-DYNA gives us that's missing in OpenSees

LS-DYNA and OpenSees occupy nearly disjoint capability spaces:

- **OpenSees strengths**: composable analysis architecture (six orthogonal axes — `ConstraintHandler` × `DOF_Numberer` × `SystemOfEqn` × `ConvergenceTest` × `EquiSolnAlgo` × `Integrator`), open source, polymorphic Element/Material framework, deterministic parallel reproducibility, Tcl/Python scripted modeling, fiber sections.
- **LS-DYNA strengths**: explicit-first solver throughput (vectorized blocks of 128 elements — Theory §29.1, p. 575), full contact subsystem with bucket sort (Theory §26), rigid-body dynamics with joint catalog (Theory §25), meshfree (SPH §38, EFG §39), modular damage (`*MAT_ADD_DAMAGE_GISSMO`), mass scaling (`*CONTROL_TIMESTEP DT2MS`).

This roadmap captures **what we want to bring from the LS-DYNA side into OpenSees, and how to do it without forking the framework**. Each port respects OpenSees' modular analysis architecture (don't merge things that should stay swappable).

---

## 2. Architectural insights driving the plan

Distilled from the LS-DYNA ↔ OpenSees comparison conversation. These are the *load-bearing* facts that drive item priorities and design choices below.

### 2.1 Where the explicit performance gap actually lives

OpenSees dispatches per element via virtual table; LS-DYNA dispatches per **block of 128 homogeneous elements** sorted by material model (Theory §29.1, pp. 575–577 — confirmed: "groups of 128 elements or possibly some larger integer multiple of 64 are utilized"). Disjoint blocks (Fig. 29.1, p. 577) make the scatter race-free.

> [!important]
> The performance cost of OpenSees' polymorphism is **concentrated in explicit analyses**. In implicit Newton, the global linear solve dominates and per-element virtual dispatch is invisible. Therefore: any "let's go fast" design effort should target the explicit element loop, not the implicit one.

### 2.2 Rigid bodies are not elements

In LS-DYNA, a rigid body is a **parallel object kind** (Theory §25, pp. 511–521): owns 6 DOFs (CoM translation + rotation), has its own SO(3)-aware integrator, kinematically slaves the nodes of its part. Joints (`*CONSTRAINED_JOINT_*`) are constraints **between** rigid bodies — they remove relative DOFs, they don't compute strain.

OpenSees has **no analogue**. `RigidLink`, `RigidDiaphragm` are `MP_Constraint`s — they relate DOFs of existing nodes, but there is no "object that owns 6 DOFs and integrates rotational dynamics." Adding one is a **new DomainComponent kind**, not a new Element. This is the prerequisite for any joint family.

### 2.3 Geometric nonlinearity — data layout determines the natural fit

Two storage conventions, with sharply different consequences:

| | Stores | Natural fit | Awkward fit |
|---|---|---|---|
| **OpenSees** | Reference config $X$ + displacement $u$ | **TL** — element computes $F = \nabla(X+u)$ on demand | UL — would require "rebasing" the reference each step |
| **LS-DYNA** | Current (deformed) config; explicit time-stepping updates nodes directly | **UL with Jaumann rate** (Theory §18.1, p. 243) | TL — element must explicitly remember reference |

Implication: OpenSees solids that need GNL today encode TL inside the element itself (e.g., Abell's `FourNodeTetrahedron`). LS-DYNA defaults to UL+Jaumann at the material interface, with TL reserved for hyperelastic and SPH-FORM=7/8.

Corotational is a third path used by *both* frameworks for slender members:
- OpenSees: `CrdTransf` (`Linear` / `PDelta` / `Corotational`) for beams/columns.
- LS-DYNA: built into the element formulation — Belytschko-Tsay shell (Theory §7.1, p. 111), Belytschko beam (§4.1, p. 73).

> [!decision] Open architectural decision
> For new explicit OpenSees solids, support **both UL+Jaumann and TL** as per-element formulation choices. UL+Jaumann is what makes block-dispatch SIMD viable (matches LS-DYNA's pattern); TL keeps consistency with the existing OpenSees nonlinear-solid convention and stays accurate for hyperelasticity. The framework should not pick one. See [[#5.10 UL+Jaumann vs TL kinematic convention]].

### 2.4 Wrapper at the integrator level avoids forcing element rewrites

The non-invasive performance path: extend `Element` with an **opt-in batch entry point** that the integrator calls for groups of same-class-tag elements. Default fallback is the existing per-element loop, so nothing breaks. New elements that opt in fill SoA buffers and run a SIMD-friendly inner loop. Heterogeneous models still work — they just batch one type at a time.

This is the architectural move that lets OpenSees adopt LS-DYNA-style throughput **without forking the Element hierarchy**. See [[#5.2 Batch element dispatch interface]].

### 2.5 Contact is its own subsystem, not "just elements"

LS-DYNA contact (Theory §26, pp. 523–568) runs every step, **before** element forces, with its own data structures (bucket sort §26.11, slave search §26.6, segment-based penalty §26.7.3). OpenSees contact today (`ZeroLength`, `Contact2D/3D`) is per-element with no broad-phase search and no segment logic — it's a different design.

A faithful port wants a **`ContactDomain` parallel to `Domain`** with its own broad-phase + narrow-phase machinery. See [[#5.7 ContactDomain subsystem]].

---

## 3. Roadmap at a glance

Items ordered for execution feasibility. Each links to a future per-feature plan that will be written from [[_template]] when the item activates.

| # | Item | Effort | Prereq | Note |
|---|------|:-----:|:------:|------|
| 5.1 | [[01_selective_mass_scaling\|Selective mass scaling integrator]] | S | — | First implementation. Bounded scope. |
| 5.2 | [[02_batch_element_dispatch\|Batch element dispatch interface]] | M | — | Non-invasive opt-in; defines the contract for fast-path elements. |
| 5.3 | [[03_soa_native_explicit_hex\|Reference SoA-native explicit hex]] | M | 5.2 | Validates 5.2 with real benchmarks. |
| 5.4 | [[04_modular_damage_decorator\|Modular damage decorator (GISSMO-style)]] | S | — | Independent track; useful for everyone. |
| 5.5 | [[05_rigid_body_domain_component\|RigidBody DomainComponent + SO(3) integrator]] | L | — | Heaviest item; gates joints. |
| 5.6 | [[06_joint_family\|Joint family as MP_Constraint subclasses]] | M | 5.5 | Capability gap closed. |
| 5.7 | [[07_contact_domain\|ContactDomain subsystem (bucket sort + node-to-segment + segment-to-segment)]] | L | 5.5 helpful | Largest capability gap. |
| 5.8 | [[08_sph_element_family\|SPH element family (FORM 1/7/8, MLS variant)]] | L | 5.7 (bucket sort) | Major research direction; meshfree story part 1. |
| 5.9 | [[09_efg_meshfree\|EFG / meshfree (MLS shape functions, integration cells)]] | L | 5.8 helpful | Meshfree story part 2. |
| 5.10 | [[10_kinematic_convention\|UL+Jaumann vs TL — supporting both]] | — | — | Architectural decision, not a single feature. Lives across 5.3 / 5.7 / 5.8. |
| 5.11 | [[11_parallel_scatter_coloring\|(Optional) Graph coloring for parallel scatter]] | M | 5.3 benchmarks | Only if single-threaded batch isn't enough. |
| 5.12 | [[12_gpu_offload\|GPU offload of explicit element kernels]] | L | 5.3 benchmarks | Single largest performance lever after batch dispatch. |
| 5.13 | [[13_matrix_free_implicit\|Matrix-free implicit (Krylov + element-local K·v)]] | L | 5.2 batch interface | Modern HPC FEM trick; scales to $10^7$+ DOF. Independent of explicit work. |
| 5.14 | [[14_xpbd_relaxation\|XPBD-style quasi-static relaxation mode]] | M | 5.2 | Research direction. Quasi-static problems with severe softening where Newton fails. |

**Effort key**: S ≈ a few hundred lines / weeks. M ≈ thousands of lines / months. L ≈ multi-paper or thesis-chapter scale.

---

## 4. Background — what each framework gives us

Quick reference table; full discussion lives in the comparison conversation and in `~/.claude/skills/ls-dyna/references/manual_index.md` §15.

| Topic | LS-DYNA | OpenSees today | Gap → action |
|---|---|---|---|
| Explicit time integration | Central difference only (Theory §24) | `CentralDifference` exists | Add **mass scaling** ([[#5.1]]) |
| Element dispatch | Vectorized blocks of 128 (Theory §29.1) | Per-element virtual dispatch | Add **batch interface** ([[#5.2]]) |
| Hex hourglass control | 6 IHQ options + Puso (Theory §3.2–3.3) | `SSPbrick` only | Future: port Puso assumed-strain |
| Shell formulations | Belytschko-Tsay, Hughes-Liu, FI #16, mesh-free (Theory §7–11) | `ShellMITC4`, `ShellNLDKGQ`, etc. | Future: explicit shell options |
| Contact | Penalty / segment-based / mortar / tied / eroding + bucket sort (Theory §26) | `ZeroLength`, `Contact2D/3D` | Add **`ContactDomain`** ([[#5.7]]) |
| SPH | `*ELEMENT_SPH` with FORM 0–16 (Theory §38) | None | Add **SPH** ([[#5.8]]) |
| EFG / meshfree | `*SECTION_SOLID_EFG` (Theory §39) | None | Add **EFG** ([[#5.9]]) |
| Rigid bodies + joints | Full joint catalog (Theory §25) | `RigidLink`/`RigidDiaphragm` (constraints only — no rigid-body objects) | Add **`RigidBody`** ([[#5.5]]) + **joints** ([[#5.6]]) |
| PML / non-reflecting | `*MAT_PML_*` material wrappers (Vol II p. 1565+) | `PML*` elements (Petracca/Abell) | Compare designs (separate research note) |
| Damage add-ons | `*MAT_ADD_DAMAGE_GISSMO` (Vol II p. 129) — modular | Per-material, not modular | Add **decorator** ([[#5.4]]) |
| Cohesive | `*MAT_138/184/185/186/240` + interface elements | `ZeroLengthContactNTS2D` and similar | Out of scope here; revisit after 5.7 |

---

## 5. Per-item plans (sketch level)

Each subsection captures the design intent; the full per-feature plan will be created in this folder under the linked filename when work begins.

### 5.1 Selective mass scaling integrator

**Why**: without this, the smallest element in the mesh bounds the global timestep — and in any realistic 3D SSI / pile-soil / contact-zone model there will be a few tiny elements (interface zones, refinement near loads) that throttle the entire run. Selective mass scaling lets us keep the fine local mesh while running at a global timestep set by the *bulk* of the model. **This is the gating prerequisite for any practical 3D explicit SSI work in OpenSees.**

**Goal**: match LS-DYNA `*CONTROL_TIMESTEP DT2MS` / `IMSCL` semantics — add fictitious mass to nodes whose elements have $\Delta t_e < \Delta t_{\text{target}}$, optionally restricted to a node set.

**Where it lives**: `SRC/analysis/integrator/CentralDifferenceSMS.cpp` (or extend existing `CentralDifference` with a flag — TBD; standalone class is cleaner).

**Algorithm**:
1. At analysis start: compute $\Delta t_e$ per element from $\ell_e / c_e$ where $c_e$ is the wave speed (every element/material can compute these).
2. For each element with $\Delta t_e < \Delta t_{\text{target}}$, compute $\Delta m_e$ such that the scaled $\Delta t_e \geq \Delta t_{\text{target}}$.
3. Distribute $\Delta m_e$ to the element's nodes.
4. Diagnostics: track total added mass; warn at >5% of model mass; warn when affected nodes carry significant inertia in the global modes of interest.

**Reference reading**:
- LS-DYNA Theory §22.1 (PDF p. 489) — solid timestep formulas.
- LS-DYNA Vol I `*CONTROL_TIMESTEP` (PDF p. 1890) — `DT2MS`, `IMSCL`, `LCTM`, `MS1ST` semantics.
- Olovsson & Simonsson 2006 (Comput. Struct.) — selective mass scaling theory; consistent vs lumped.

> [!question] Open
> Lumped or consistent mass scaling? Lumped is what LS-DYNA does by default; consistent (Olovsson) preserves global frequencies better but requires solving an extra small system per element block. Probably start lumped, study consistent later.

→ Future plan: [[01_selective_mass_scaling]].

---

### 5.2 Batch element dispatch interface

**Why**: closes the explicit-throughput gap. Without this, OpenSees explicit will be 5-10× slower than LS-DYNA on the same mesh, and the platform argument fails on performance grounds even if the capabilities are there. The opt-in design means **we don't fork the framework** — every existing element keeps working exactly as it does today; only new elements that opt in get the fast path.

**Goal**: opt-in fast path that lets the explicit integrator process homogeneous element groups without per-element virtual dispatch. Existing elements unchanged.

**Where it lives**:
- Modify: `SRC/element/Element.h` (add the new virtual with default impl in `Element.cpp`).
- Modify: `SRC/analysis/integrator/CentralDifference.cpp` (use the batch entry).

**Interface sketch**:

```cpp
class Element : public DomainComponent {
public:
    // ... existing virtuals unchanged ...

    // Optional fast path. Default fallback dispatches to per-element calls,
    // so existing elements keep working without modification.
    virtual int batchAssembleResisting(
        const std::vector<Element*>& batch,    // homogeneous: same class tag
        Vector& globalForce,
        const AnalysisModel& model);
};
```

**Integrator side**:

```cpp
std::map<int, std::vector<Element*>> byType;
for (Element* e : *domainPtr)
    byType[e->getClassTag()].push_back(e);

for (auto& [tag, batch] : byType)
    batch[0]->batchAssembleResisting(batch, globalForce, *theModel);
```

**Why this works**: in explicit, only the internal force vector is needed (no global tangent assembly). The batch path can fill SoA buffers (positions, velocities, history vars), run a SIMD-friendly inner loop, scatter once.

**Why this only helps explicit**: implicit Newton is dominated by the linear solve; per-element virtual cost is invisible there.

> [!warning] Race-free scatter
> First version: single-threaded scatter is fine and already a big win. For shared-memory parallel scatter, need graph coloring of element connectivity (LS-DYNA Theory §29.1, p. 577, Fig. 29.1). Defer to [[#5.11]].

> [!note] Leave room for matrix-free implicit
> The batch interface should be designed so that an analogous `batchApplyTangent(batch, v_in, Kv_out)` can be added later **without redesigning the contract**. Matrix-free Krylov methods (CG, GMRES with element-local $K \cdot v$ application instead of forming $K$ globally) are the modern HPC FEM trick — see [[#11.3]]. We're not building it now, but if we name and shape `batchAssembleResisting` such that a sibling method fits the same data layout, future-us has the option. Mirrors what MFEM and deal.II do.

→ Future plan: [[02_batch_element_dispatch]].

---

### 5.3 Reference SoA-native explicit hex

**Why**: the canonical 3D continuum element for explicit SSI / soil / concrete work. Validates the [[#5.2]] interface with a real production element so the contract is right *before* downstream items lock against it. Also: the benchmark target — head-to-head against an existing OpenSees `stdBrick`/`bbarBrick` and against LS-DYNA `ELFORM=1` with the same material — that tells us whether the batch interface actually pays off.

**Goal**: one production-quality 8-node hex element built around the [[#5.2]] interface from day one. Validates the contract and produces benchmark numbers.

**Decisions**:
- **Kinematic convention**: probably **UL with Jaumann** (matches LS-DYNA ELFORM=1; SoA-friendly). Could add a TL variant later. See [[#5.10]].
- **Hourglass control**: start with Flanagan-Belytschko viscous (IHQ=2); Puso assumed-strain (IHQ=9 in LS-DYNA, Theory §3.3 pp. 46–52) is a follow-up.
- **Material interface**: `NDMaterial::setTrialStrainIncr` already exists; needs to handle the rate-of-deformation input correctly.

**Reference**:
- LS-DYNA Theory §3.1 (PDF p. 39) — volume integration / 1-pt vs full integration.
- LS-DYNA Theory §3.2 (PDF pp. 40–46) — Flanagan-Belytschko hourglass.

→ Future plan: [[03_soa_native_explicit_hex]].

---

### 5.4 Modular damage decorator (GISSMO-style)

**Why**: turns *any* plasticity model in OpenSees (`J2`, `PressureDependentMultiYield`, `PM4Sand`, custom `NDMaterial`s) into a damage-capable model — without rewriting the constitutive integration. Today, every concrete / damage / failure model in OpenSees is a one-off (`Concrete02` has its own damage; `ConcreteCM` has its own; `ASDConcrete3DMaterial` has its own). A decorator pattern means a single damage implementation is reusable across the catalog. **Also: works in implicit *and* explicit** — independent of the rest of the explicit roadmap, immediately useful.

**Goal**: a decorator `NDMaterial` that wraps any plasticity-capable `NDMaterial` and adds damage with stress-state-dependent thresholds. Mirrors `*MAT_ADD_DAMAGE_GISSMO` (LS-DYNA Vol II p. 129).

**Why it's worth doing now**: independent of explicit work; useful for any seismic or impact problem where plasticity + damage matters. Should be ~500 lines.

**Interface sketch**:

```cpp
class GISSMODamageMaterial : public NDMaterial {
    NDMaterial* base;          // wrapped plasticity
    // damage state: D, equivalent plastic strain, stress triaxiality history
public:
    int setTrialStrain(const Vector& strain) override {
        base->setTrialStrain(strain);
        // update damage from base's plastic strain increment + current triaxiality
        // degrade stress / tangent by (1-D)
    }
};
```

**Reference**:
- LS-DYNA Vol II `*MAT_ADD_DAMAGE_GISSMO` (PDF pp. 129–141).
- Neukamm, Feucht, Haufe (DYNAmore papers, 2008+) — GISSMO theory.

> [!warning] Regularization is non-negotiable
> This decorator MUST accept fracture energy $G_f$ and apply crack-band scaling on day one. Without that, it's mesh-dependent garbage. See [[#10.5.1]] for the requirement spec and [[#10]] for the full regularization context.

→ Future plan: [[04_modular_damage_decorator]].

---

### 5.5 RigidBody DomainComponent + SO(3) integrator

**Why**: enables the entire class of problems where parts of the model are stiff enough that treating them as deformable is just wasted compute — and where their *rotational* motion is what matters. **Foundation rocking** (currently modeled as zero-length rotational springs, losing the actual rigid-body kinematics); **multi-body machine-foundation systems**; **rigid pile caps** under deformable pile groups; **base-isolated structures** where the isolator separates a rigid superstructure from a deformable foundation; **anchor / weight blocks** in retaining systems. Also: prerequisite for the joint family ([[#5.6]]) — joints are constraints *between* rigid bodies, so the rigid-body abstraction has to exist first.

**Goal**: a new `DomainComponent` kind (parallel to `Element` and `MP_Constraint`) representing a rigid body with 6 DOFs (translation + rotation), its own mass + inertia tensor, and an SO(3)-aware integrator for the rotational state.

**Why it's heavy**:
- Rotation can't be integrated with plain central difference because Euler angles don't compose linearly. Need quaternion + exponential map, or Lie-group integrator.
- Inertia tensor must be transformed each step.
- Constraint handler interaction: when a `RigidBody` slaves a node, the framework needs to know that node's DOFs are derived, not free.

**Key references**:
- LS-DYNA Theory §25 (PDF pp. 511–521) — LS-DYNA's rigid-body dynamics.
- Simo & Vu-Quoc 1986 (Comput. Methods Appl. Mech. Engrg.) — the canonical reference for rigid-body dynamics with large rotations in finite-element context.
- Betsch & Steinmann 2001 — energy-momentum conserving schemes for rigid bodies.
- Holzapfel "Nonlinear Solid Mechanics" Ch. 2 for rotation parameterization.

**Where it lives**: `SRC/domain/rigid/RigidBody.{h,cpp}` plus an integrator extension or a side-channel handler in the explicit step.

> [!question] Open
> Should `RigidBody` be a standalone `DomainComponent` (cleanest) or a special `Element` (existing infrastructure handles it more easily)? Standalone is more correct architecturally but requires more framework changes.

→ Future plan: [[05_rigid_body_domain_component]].

---

### 5.6 Joint family as MP_Constraint subclasses

**Why**: opens dynamic multi-body modeling — base isolators and friction pendulums (revolute / spherical), rocking foundations with constrained tipping (planar with limits), mechanical hinges in seismic-protective systems, vehicle suspensions for bridge-vehicle interaction, dampers connecting rigid blocks. Today these are built awkwardly with `zeroLength` elements + nonlinear materials; joints make the relative-motion constraint *the* primitive, which is what it actually is mechanically. Architectural advantage over LS-DYNA: because joints are `MP_Constraint`s in OpenSees, the user picks penalty / Lagrange / transformation through `ConstraintHandler` — LS-DYNA hardcodes penalty.

**Goal**: revolute, spherical, cylindrical, planar, universal, translational, locking — the LS-DYNA `*CONSTRAINED_JOINT_*` catalog. Built on top of [[#5.5]].

**Architecture**: each joint type is an `MP_Constraint` subclass that eliminates specific relative DOFs between two `RigidBody` objects. Preserves the OpenSees `ConstraintHandler` swap (penalty / Lagrange / transformation) — that's the architectural advantage over LS-DYNA, where penalty is hardcoded.

**Reference**:
- LS-DYNA Vol I `*CONSTRAINED_JOINT_*` (PDF pp. 939–989).
- Theory §25.1 (PDF p. 513).

→ Future plan: [[06_joint_family]].

---

### 5.7 ContactDomain subsystem

**Why**: **the single largest capability gap between OpenSees and the rest of the nonlinear-FEM world.** Without it, none of the urban / SSI contact-driven research is possible: pounding between adjacent buildings, foundation uplift, pile-soil gap formation, soil-retaining-wall sliding, slab-on-soil separation, sliding bearings, friction pendulums (post-uplift), rocking foundations with toe contact, post-buckling member-to-member contact in collapsing frames. The current OpenSees contact (`ZeroLength`, `Contact2D/3D`) requires the user to pre-define every contact pair, which is fine for simple bench problems but breaks for any realistic geometry. A real `ContactDomain` with a broad-phase search makes contact a *first-class capability* the way it is in LS-DYNA / Abaqus / Code_Aster.

**Goal**: a contact engine parallel to `Domain`, with broad-phase (bucket sort) + narrow-phase (penalty / segment-based / tied) algorithms running every step.

**Components**:
1. **Bucket sort** (LS-DYNA Theory §26.11, PDF p. 545–550). Spatial hash with cell size ≈ max element size. Reused by SPH ([[#5.8]]).
2. **Slave search** + closest-point projection (Theory §26.6).
3. **Penalty algorithms**:
   - Standard node-to-segment penalty (Theory §26.7.1) — start here.
   - Soft-constraint penalty `SOFT=1` (§26.7.2) — Courant-stable stiffness.
   - Segment-based penalty `SOFT=2` (§26.7.3) — robust for corner / edge / T-intersection cases.
4. **Tied interfaces** (§26.9) — for bonded contact / mesh-tied submodeling.
5. **Friction**: Coulomb radial-return.

**Why heavy**: ~hundreds of thousands of lines in LS-DYNA. For OpenSees, target a focused subset (node-to-segment penalty + bucket sort + Coulomb friction) as v1.

**Reference**:
- LS-DYNA Theory §26 (PDF pp. 523–568) in full — the canonical chapter.
- Wriggers "Computational Contact Mechanics" — textbook treatment.
- Puso & Laursen 2004 — mortar contact (defer to v3).

→ Future plan: [[07_contact_domain]].

---

### 5.8 SPH element family

**Why**: opens the class of problems with **severe deformation, fragmentation, free-surface flow, and topology change** — situations where a mesh-based formulation simply can't continue: tsunami debris impact on structures, coseismic landslides and slope failure, post-failure granular flow (sand boils, lateral spreading runout), tunnel collapse, fragment ejection from spalling concrete under blast / impact. These problems exist in the user's research domain (Chile / Andean geomorphology, post-earthquake reconnaissance, blast-resistant design) and currently have **no path** in OpenSees. SPH is the well-established meshfree choice; getting it into OpenSees alongside the existing FEM elements (with FEM-SPH coupling via contact) is the bridge to a much larger problem class.

**Goal**: SPH for deformable solids, with multiple formulation variants:
- FORM=1 — standard Eulerian-kernel continuity-density form.
- FORM=7/8 — total Lagrangian SPH (kernel anchored to reference, kills tensile instability).
- FORM=12 — MLS-based.

**Components**:
- `SPHParticle` element (or a particle subclass).
- `NeighborList` helper (reuses bucket sort from [[#5.7]]).
- Cubic B-spline kernel (Theory Eqs. 38.1–38.3, PDF p. 637).
- Kernel gradient with renormalization (Eq. 38.6).
- Variable smoothing length evolution (DERIV ODE — `*CONTROL_SPH` Vol I p. 1844).
- Monaghan artificial viscosity for shocks (Theory Eqs. 38.22–23, p. 641).
- CFL timestep (Eq. 38.24, p. 642).

**Constitutive interface**: SPH particles use `NDMaterial` for stress update — this matches OpenSees' existing material interface. SPH-specific damage / failure handled at particle level (e.g., particle deletion).

**Reference**: LS-DYNA Theory §38 (PDF pp. 637–643) and Vol I `*CONTROL_SPH` / `*SECTION_SPH`. See also iteration-1 eval-2 answer in `~/.claude/skills/ls-dyna-workspace/iteration-1/eval-2-sph-deformable-solid/with_skill/outputs/answer.md`.

→ Future plan: [[08_sph_element_family]].

---

### 5.9 EFG / meshfree

**Why**: meshfree alternative for the regime where standard FEM fails from mesh distortion but you don't need full SPH-style fragmentation — large-strain forming, large-strain plasticity in geomaterials before they actually fluidize, **fracture initiation in continuum without pre-inserted cohesive surfaces**, and the dynamic crack-propagation problems where XFEM was historically the answer (the user has previously flagged XFEM as a research interest). EFG retains a continuum stress field and convergence proofs that SPH lacks; the cost is heavier per evaluation. Together, [[#5.8]] and [[#5.9]] cover the full meshfree spectrum from continuum-with-flexibility (EFG) to fully-discontinuous (SPH).

**Goal**: Element-Free Galerkin for solids. Continuation of the meshfree push.

**Components**:
- MLS shape functions (Theory §39.1, PDF p. 645).
- Integration constraint + strain smoothing (§39.2 — what makes EFG converge).
- Lagrangian strain smoothing for path-dependent problems (§39.3).
- Imposition of essential BCs (§39.5 — known weak spot for meshfree).
- Optional mesh-free shell variant (§39.6).

**Reference**: LS-DYNA Theory §39 (PDF pp. 645–657).

→ Future plan: [[09_efg_meshfree]].

---

### 5.10 UL+Jaumann vs TL — supporting both

**Why**: the kinematic convention determines what hyperelasticity, plasticity, and explicit performance look like at the element level. Different problem classes need different conventions: hyperelastic rubber wants TL, J2 plasticity at high rate wants UL+Jaumann, SPH-FORM=7/8 wants TL. **Forcing one convention bakes a research limit into the framework.** The decision here is to *not* decide framework-wide — make space for both, document the tradeoffs, let each per-element implementation choose.

**Not a single feature** — an architectural decision that affects 5.3, 5.7, 5.8 and any future explicit element work.

The framework should support:
- **TL formulations** (existing OpenSees convention; natural given OpenSees stores reference config + displacement).
- **UL formulations with Jaumann objective rate** (LS-DYNA convention; required to make block-dispatch SIMD viable for plasticity-dominated explicit problems).

Per-element choice. Document the tradeoff:

| Aspect | TL | UL+Jaumann |
|---|---|---|
| Storage layout | Reference-anchored — natural in OpenSees | Each element must remember reference internally |
| Performance | Heavier per element (recompute $F$ each step) | Lighter, SoA-friendly inner loop |
| Hyperelasticity | Native — $W(F)$ is reference-frame | Requires extra rotation handling, Jaumann pathology under simple shear |
| Plasticity | Standard | Standard, matches LS-DYNA's default |
| Use case | High-strain hyperelastic, careful nonlinear statics | High-rate, large-rotation, plasticity-dominated explicit |

**Recommendation**: new explicit elements default to UL+Jaumann; TL stays available for hyperelastic and SPH-FORM=7/8.

→ Future plan: [[10_kinematic_convention]] (will collect the decision rationale and the per-element choices made downstream).

---

### 5.11 (Optional, later) Graph coloring for parallel scatter

**Why**: shared-memory parallel speedup on top of [[#5.2]] / [[#5.3]]. The block-dispatch interface alone removes virtual-call overhead; parallel scatter removes the assembly-side bottleneck that remains. Only pursue if benchmarks demand it, because the implementation complexity (graph coloring, race-free coloring updates after element birth/death) is non-trivial and OpenSees' existing single-threaded determinism is something we don't want to give up casually.

Only pursue if [[#5.3]] benchmarks show single-threaded batch isn't enough.

**Reference**: LS-DYNA Theory §29.1 (PDF pp. 575–577, especially Fig. 29.1 — "Group of 48 elements broken into 4 disjoint blocks"). The GM/Katnik/Benson 1988–89 approach.

→ Future plan: [[11_parallel_scatter_coloring]].

---

### 5.12 GPU offload of explicit element kernels

**Why**: explicit dynamics is the GPU-natural workload in scientific FEM. There's no global linear solve to fight (the hard part of GPU FEM); the inner loop is **gather → element kernel → scatter**, exactly the pattern GPUs eat for breakfast. LS-DYNA has had GPU offload for a decade; OpenSees has nothing. Once [[#5.2]] / [[#5.3]] are stable, this is the single largest performance lever still on the table — credibly **another 2–3× on top of CPU batch dispatch** for typical SSI / continuum problems.

**Where it lives**:
- New: `SRC/element/gpu/` for GPU kernels of fast-path elements.
- New: a portability layer (Kokkos / RAJA / SYCL / a Taichi-like DSL — Q below).
- Modify: `CMakeLists.txt` for CUDA / HIP / oneAPI build options.

**Design constraints inherited from [[#5.2]]**:
- Element batch is already SoA — that's the GPU-friendly layout.
- The batch is already homogeneous in element type — perfect for one CUDA kernel per type.
- The scatter race problem ([[#5.11]]) becomes the GPU race problem; both want graph coloring.

**Likely target stack**:
- **CUDA + HIP via portability layer**. Don't write CUDA directly; that locks us to NVIDIA. Kokkos (used by MFEM, Trilinos) and RAJA (used by LLNL codes) are the proven options. SYCL (Intel oneAPI) is the open-standard alternative.
- **Material kernels are the hard part**, not element kernels. Element kernels are mostly numerical integration of $\int B^T \sigma$ — straightforward parallel patterns. Material kernels (PM4Sand, ASDPlasticMaterial3D, CDPM) have branchy state-machine code that fights GPU SIMT execution. Need to refactor materials into "pure-functional" form: `(strain_n, state_n) → (stress, state_{n+1})` with no hidden mutation, so they map to GPU threads.

**Realistic scope for v1**:
- One element type (the SoA hex from [[#5.3]]).
- One material model (something simple — `ElasticIsotropic3D` or `J2Plasticity`).
- Single-GPU only (no multi-GPU until benchmarks justify it).
- Gather/kernel/scatter on GPU; everything else (integrator update, recorders, output) on CPU with PCIe transfer per step.

**Hard limits**:
- PCIe bandwidth caps the speedup ceiling. If the analysis is dominated by data transfer (small mesh, frequent CPU sync), GPU might be slower than CPU. Profile before celebrating.
- Branchy materials may give little GPU speedup. The classic "GPU wins by 100×" stories use simple elastic kernels; nonlinear plasticity is much closer to 2–3×.
- Determinism: order of summation in scatter is non-deterministic on GPU without atomics; atomics are slow. Same tradeoff discussed in [[#8.2]].

**Reference**:
- MFEM GPU implementation (LLNL) — `mfem/fem/bilinearform_kernels.cpp`-style files.
- NVIDIA Warp (https://github.com/NVIDIA/warp) — Python-fronted GPU physics, useful pattern.
- LS-DYNA's documented GPU element types (LS-DYNA R10+ release notes; Vol I `*CONTROL_GPU`).

> [!question] Open
> Portability stack — Kokkos (heaviest, most proven), RAJA (LLNL favorite), SYCL (open standard), or a Taichi-style DSL with custom code generation? Pick after profiling [[#5.3]]; the answer depends on whether materials can stay in C++ or need a DSL.

→ Future plan: [[12_gpu_offload]] (lower priority than [[#5.7]] / [[#5.8]]; revisit after benchmarks from [[#5.3]] confirm CPU batch dispatch isn't enough for the target problem sizes).

---

### 5.13 Matrix-free implicit (Krylov with element-local $K \cdot v$)

**Why**: closes the implicit-FEM scaling ceiling. OpenSees' `system SparseGEN` (and friends) form the global stiffness $K$ and store it sparsely. For 3D problems above ~$10^6$ DOF this becomes memory-bound — the matrix itself is gigabytes, the sparse direct factorization is super-linear, and you hit a wall well below the problem sizes that modern HPC FEM codes (MFEM, deal.II) handle routinely. **Matrix-free Krylov methods change the cost model**: never form $K$, never store it; instead, compute $K \cdot v$ on demand by looping over elements and applying $k_e \cdot v_e$. Memory-bandwidth-limited rather than memory-capacity-limited; scales to $10^7$+ DOF; GPU-friendly.

This is the single biggest implicit-FEM scaling improvement available, and it's **independent of the explicit work** — useful for the static / pushover / modal side of OpenSees that this roadmap otherwise leaves alone.

**Where it lives**:
- Modify: `SRC/element/Element.h` to add an opt-in `batchApplyTangent(batch, v_in, Kv_out)` mirror of `batchAssembleResisting` ([[#5.2]]).
- New: `SRC/system_of_eqn/matrixFree/MatrixFreeSOE.{h,cpp}` — implements the `LinearSOE` interface without storing $K$; defers to elements for the apply.
- New: `SRC/analysis/algorithm/` extensions — Krylov solvers (CG for SPD, GMRES for general) that work against `MatrixFreeSOE`.
- New: preconditioners under `SRC/system_of_eqn/precond/` — Jacobi (cheap), block-Jacobi (per-element block inverse), AMG (via hypre or AMGX wrapper for the heavy version).

**Algorithm sketch**:

```cpp
// Element opt-in: sibling to batchAssembleResisting
class Element {
public:
    virtual int batchApplyTangent(
        const std::vector<Element*>& batch,
        const Vector& v,         // global trial vector
        Vector& Kv,              // output: K * v
        const AnalysisModel& model);
};

// MatrixFreeSOE just remembers the domain; "solve" runs Krylov iteration
class MatrixFreeSOE : public LinearSOE {
    int solve() override {
        return cg_or_gmres(*this, /*matvec=*/&MatrixFreeSOE::applyK,
                           rhs, solution, preconditioner, tolerance);
    }
    void applyK(const Vector& v, Vector& Kv) {
        // bucket elements by class tag (same as #5.2)
        // call batchApplyTangent on each bucket
        // assemble Kv (with constraint handler contribution)
    }
};
```

**Why this works**: matrix assembly is the memory bottleneck in implicit FEM. Element-local $k_e \cdot v_e$ is cache-friendly (gather → small matvec → scatter, the same pattern as explicit). For a 20-iteration GMRES solve, you do 20 element-loop matvecs vs one assembly + one factorization; the matvecs are *cheap per iteration* and memory-bandwidth-limited, while the factorization is super-linear in DOF count.

**The catch — preconditioning is everything**:
- Without a preconditioner, Krylov convergence depends on $\kappa(K)$ (condition number). For mechanics problems it can be in the millions; convergence stalls.
- With a good preconditioner (block-Jacobi, AMG, or domain-decomposition with local solves), convergence is back to a few dozen iterations.
- For ill-conditioned problems (very stiff/soft contrasts in SSI, near-singular tangent during plasticity localization), Krylov can stall hard. **Sparse direct stays the more robust fallback** for those cases — the architectural decision is "matrix-free *available*," not "matrix-free *only*."

**What it unlocks**:
- 3D continuum problems with $10^7$ DOF (today: limited by sparse factorization memory).
- Implicit-explicit hybrid analyses (subcycling).
- Implicit GPU FEM (matrix-free is the only practical route there).
- Differentiable implicit solves (auto-diff through Krylov is much cleaner than through sparse factorization — relevant for [[#12]]).

**Reference**:
- deal.II matrix-free guide: https://www.dealii.org/current/doxygen/deal.II/group__matrixfree.html — the canonical reference.
- MFEM `mfem/fem/bilinearform.cpp` (matrix-free bilinear form integration) — production HPC code.
- Trottenberg, Oosterlee, Schüller "Multigrid" (2001) — for AMG preconditioning theory.
- Saad "Iterative Methods for Sparse Linear Systems" (2nd ed., SIAM 2003) — Krylov reference.

> [!question] Open
> **Initial preconditioner choice**: ship Jacobi only (cheap, sometimes enough), Jacobi + block-Jacobi (covers most of mechanics), or go straight for AMG via a hypre wrapper (best convergence but heavyweight dependency)? Probably block-Jacobi first; AMG when block-Jacobi convergence isn't enough.

> [!question] Open
> **Matrix-free for nonlinear problems**: line search and trust-region methods need an inexact $K$ to make iterate decisions. Can the matrix-free $K$ be queried with a coarse approximation cheaply, or do we always pay the full element loop? Newton-Krylov literature has answers (Knoll & Keyes 2004 review); needs a per-feature plan to commit.

→ Future plan: [[13_matrix_free_implicit]].

---

### 5.14 XPBD-style quasi-static relaxation mode

**Why**: opens an under-explored regime — **quasi-static problems with severe softening where Newton fails**. Today's options for these are arc-length (fragile under cyclic loading), explicit dynamic with damping (energy-balance polluted, mass-scaling artifacts), or just "we couldn't run that case." **Position-Based Dynamics** (Müller et al. 2007) and especially **XPBD** (Macklin et al. 2016) are the modern game-physics descendants of dynamic relaxation. They iterate constraint corrections at the position level until the residual relaxes. Bring the engineering accuracy bar back, and you get a third quasi-static option: robust through softening, no dynamic artifacts, and (importantly) cheap per iteration.

This is **research territory** — the marriage of XPBD's mathematics with engineering constitutive models hasn't been done at production quality. But there's a clear research thread, and the architectural fit with [[#5.2]] is good.

**The basic idea**: instead of solving $K \cdot du = R$ (Newton) or integrating $\ddot{u} = M^{-1}(F_{\text{ext}} - F_{\text{int}})$ (explicit dynamics), iterate position corrections directly:

$$\Delta u_i = -\frac{C_i + \alpha_i \, \lambda_i}{\nabla C_i^T M^{-1} \nabla C_i + \alpha_i}, \qquad \lambda_i \mathrel{+}= \Delta \lambda_i$$

where $C_i$ is a constraint (in mechanics, typically the local stress-strain residual mapped to a position constraint), $\nabla C_i$ is its gradient, $\alpha_i = 1/k_i$ is the **compliance** (XPBD's key innovation — physical stiffness recovery), and $\lambda_i$ are accumulating Lagrange multipliers.

You loop element-by-element (Gauss-Seidel) applying constraint corrections, sweep until residuals drop below tolerance. **It's matrix-free by construction** — only local element data is needed at each step, exactly the [[#5.2]] pattern.

**Where it lives**:
- New: `SRC/analysis/integrator/XPBDRelaxation.{h,cpp}` — the relaxation driver.
- New: `Element` extension for constraint-form output — each element exposes its residual as a position-constraint with stiffness/compliance info.
- Or: reformulate constitutive update as constraint gradient — reuses the existing `setTrialStrain` / `getStress` / `getTangent` interface, the relaxation loop just calls them differently.

**What this unlocks**:
- **Pushover with severe softening** — RC frame collapse where every Newton step fights damage; XPBD just relaxes through.
- **Form-finding** for cable nets, tensile structures, soft inflatables — XPBD's home turf in graphics, directly applicable.
- **Pseudo-static SSI** with foundation gap-opening (currently quasi-static contact is a nightmare in OpenSees).
- **Soft-soil consolidation** at the slow-loading end where dynamic explicit is overkill.

**The honest caveat**: XPBD with a real constitutive model is a research problem. The Macklin paper handles elastic and elastoplastic; nonlinear hardening, kinematic plasticity, multi-surface plasticity (PM4Sand, SANISAND), damage — these need work. **The first per-feature plan should target a narrow validated case** (elastic-perfectly-plastic 3D continuum) to establish the technique, then extend.

**Reference**:
- Müller, Heidelberger, Hennix, Ratcliff (2007) "Position Based Dynamics," J. Visual Comm. Image Repr. — the original.
- Macklin, Müller, Chentanez (2016) "XPBD: Position-Based Simulation of Compliant Constrained Dynamics," Proc. ACM MIG — the modern compliant version.
- LS-DYNA Theory §31 (PDF pp. 593–598) — dynamic relaxation, the conventional engineering analogue. XPBD generalizes this with proper compliance treatment.
- NVIDIA FleX papers — production-quality XPBD with elastoplasticity.

> [!question] Open
> Is XPBD a **third "analysis mode"** alongside static and transient (modify the OpenSees `analysis` command to take `XPBDRelaxation`), or a special **integrator within static** (similar to how arc-length is)? Probably its own mode — the iteration structure is too different from Newton.

→ Future plan: [[14_xpbd_relaxation]].

---

### 5.15 Known gaps not yet prioritized

Items the fracture / buckling motivation in [[#1.1]] surfaced but that haven't been promoted to full roadmap entries yet. Listed here so the gaps are visible and don't get forgotten.

- **Cohesive / interface elements** for discrete cracking. LS-DYNA `*MAT_138/184/185/186/240` family (Vol II pp. 1022, 1283–1289, 1605). Needed for: discrete crack propagation in concrete, mode-I/II/III fracture energy-based debonding, delamination of composite layers. Architecturally fits as a new `Element` family + `UniaxialMaterial` (or small `NDMaterial`) for the traction-separation law. **Modular damage ([[#5.4]]) handles smeared damage; cohesive elements handle discrete cracks** — both are needed for the full concrete-fracture story.
- **Element erosion machinery** for failed-element removal. Currently OpenSees has no mechanism to remove an element from the global assembly mid-run. Needed for: concrete spalling, fragment ejection, post-failure of softened soil zones. Architecturally a small `Domain` extension (an "alive" flag on each element + assembly skipping) plus a recorder for erosion events. Smaller than it sounds; could be done alongside [[#5.4]].
- **Explicit shells.** Buckling and large-deformation collapse need shell elements compatible with the batch interface ([[#5.2]]). Path forward: adapt Petracca's `ASDShellQ4` to opt into the batch interface. Rich subject — Belytschko-Tsay vs Hughes-Liu vs MITC4 vs the Hu-Washizu fully-integrated #16 (LS-DYNA Theory §§7, 9, 10). Probably a thesis chapter on its own.
- **XFEM / phase-field fracture.** The user has flagged XFEM as a research interest. EFG ([[#5.9]]) is a step in that direction; XFEM would add enrichment functions to the standard FE basis at known crack locations. Different design from EFG — defer until [[#5.9]] is in.
- **DEM (discrete element method) coupling.** Granular media + structure problems. Could share infrastructure with SPH ([[#5.8]]) and contact ([[#5.7]]) but is its own particle type. Long-term meshfree direction.

These will graduate to numbered items once the foundational work (5.1-5.7) is far enough along that prioritizing them is meaningful.

---

## 6. Open architectural questions

These are the cross-cutting decisions that affect multiple roadmap items. They're collected here (rather than buried in per-feature plans) so the design tradeoffs can be debated as a whole. Each has a current leaning ("decision direction") but that direction is **not final** — these stay open until the relevant per-feature work is far enough along to make the call with real data. Each question links back to the per-item plans where it first arose.

### Q1 — RigidBody as DomainComponent or Element?

**Where it arises**: [[#5.5]] (RigidBody implementation), [[#8.7]] (architectural tradeoff).

**The choice**:
- **(a) Standalone `DomainComponent`** parallel to `Element` and `MP_Constraint`. New top-level object kind. Rigid bodies have their own iteration in the analysis step.
- **(b) Special `Element`** subclass that overrides `getResistingForce` to return zero, `getTangentStiff` to return zero, and uses internal SO(3) state.

**Pro (a) — DomainComponent**:
- Rigid bodies *aren't* elements mechanically — they don't compute strain, don't have a constitutive model, don't fit the per-Gauss-point pattern. The abstraction matches reality.
- Recorders, parallel decomposition, energy-balance computation, output formatters never need to "remember to skip rigid elements" — those code paths just iterate the components they care about.
- Rotational DOFs live naturally on the rigid body itself; nodes constrained to it are pure followers.
- Matches LS-DYNA's mental model (Theory §25 treats them as a parallel subsystem).

**Pro (b) — special Element**:
- Smaller initial framework change. The element-loop infrastructure (assembly, parallel send/recv, recorder hookups) is already wired up; piggyback on it.
- Faster to a working prototype.
- Easier to merge upstream (`git diff` is smaller; reviewers see fewer touched files).

**Hidden cost analysis**:
- Approach (b) leaks the rigid-body-ness into every place that iterates elements. Today maybe 5 such places (recorders, energy, parallelism, output, broker). In two years, 10. Each code change has to remember "is this element rigid? skip if so." Compounds as debt.
- Approach (a) costs more upfront — touching the analysis loop, recorder framework, and `sendSelf/recvSelf` infrastructure to make them aware of a new component kind — but is a *one-time* cost.

**Decision direction**: **(a) DomainComponent**. Pay the architectural cost upfront. The "cleaner abstraction" argument compounds in our favor over the lifetime of the codebase.

**Decision criterion to revisit**: after the prototype implementation in [[#5.5]], if the framework changes balloon (more than ~15 files touched in `SRC/` outside `domain/rigid/`), reconsider whether (b) is the pragmatic move.

---

### Q2 — Lumped or consistent mass scaling?

**Where it arises**: [[#5.1]] (selective mass scaling implementation).

**The choice**:
- **Lumped scaling** (LS-DYNA `DT2MS`): added inertia $\Delta m_e$ at a node = whatever lumped fraction makes $\Delta t_e = \Delta t_{\text{target}}$. Diagonal mass matrix preserved.
- **Consistent scaling** (Olovsson & Simonsson 2005, *Comput. Methods Appl. Mech. Engrg.*): scales mass *anisotropically* along the high-frequency element modes only, leaving low-frequency response untouched. Requires per-element block-mass treatment.

**Why lumped is the LS-DYNA default**: it's simple, robust, and the diagonal-mass property is what makes explicit explicit — you invert mass element-wise, no global solve. Lumped scaling preserves this.

**Why consistent matters for our domain**: Olovsson's papers show that **for the same accuracy target, consistent scaling allows ~10× more aggressive mass scaling** than lumped. For seismic SSI specifically — where preserving the structure's fundamental period is non-negotiable — that's a very large deal. Lumped scaling shifts modes; consistent doesn't (or does so far less).

**Cost of consistent**:
- Mass becomes block-diagonal (per-element) rather than fully diagonal. Inverting it requires a per-element block solve every step.
- Breaks the "lumped diagonal mass" assumption that pervades explicit OpenSees machinery. Probably needs to be its own analysis mode rather than a flag on the existing one.
- Implementation is meaningfully heavier than lumped — Olovsson reports several hundred lines of element-level code per element type.

**Decision direction**: **Ship lumped first** ([[#5.1]] v1). Document its frequency artifacts honestly; provide modal-shift diagnostics. **Park consistent as a separate research item** ([[#5.1]] v2 or its own future plan).

**Decision criterion to revisit**: if benchmarks on representative SSI models show lumped mass scaling shifts the fundamental period by >1% at the mass-scaling levels needed to make the analysis tractable, consistent moves up the priority list.

---

### Q3 — Integrator-side or Element-side batch dispatch?

**Where it arises**: [[#5.2]] (batch dispatch interface), [[#9]] (deep dive on dispatch).

**The choice**:
- **(a) Integrator-side** (current sketch in [[#5.2]]): the existing `Domain` holds individual `Element*`s; the integrator calls `batchAssembleResisting(batch, ...)` after grouping by `getClassTag()` each step.
- **(b) Element-group-side**: a new `ExplicitElementGroup` `DomainComponent` owns N homogeneous elements; the user (or a builder) explicitly creates groups; the integrator iterates groups directly.

**Pro (a) — integrator-side**:
- Non-invasive. Existing Tcl/Python `element` commands work unchanged. Existing models look identical.
- Heterogeneous models get partial speedup automatically (whichever element types opt into the batch interface get the fast path).
- No changes to `Domain` ownership semantics or serialization.

**Con (a)**:
- Groups are recomputed every step (cheap — just a `std::map<int, vector<Element*>>` build over $N$ elements — but real overhead). Could cache and invalidate on element birth/death.
- The integrator itself has to be aware of the batch interface. Adds coupling between analysis and element infrastructure.

**Pro (b) — element-group-side**:
- More explicit: the user knows their elements are batched. Easier to reason about performance.
- SoA buffers can be allocated once at group creation, not rebuilt each step.
- Natural unit for GPU offload later — a group is exactly what you'd ship to a CUDA kernel.
- Group-level recorders become possible (energy per group, hourglass per group, etc.).

**Con (b)**:
- Invasive. Existing models don't get the speedup unless rewritten as groups (or unless a "auto-grouping builder" is added that constructs groups from existing element lists).
- Serialization (`sendSelf`/`recvSelf`) becomes complex — a group has to send all its sub-elements, and parallel decomposition has to decide whether groups are atomic units or split-able.
- Tcl/Python API needs new commands.

**Decision direction**: **(a) integrator-side for v1** ([[#5.2]] / [[#5.3]]). Cheapest path to benchmarks; non-invasive; existing models benefit. **(b) revisit as v2** if grouping overhead at every step is measurably costly, or when GPU offload becomes a priority.

**Decision criterion to revisit**: profiling [[#5.3]] benchmarks. If the group-rebuild step shows up as >5% of total time, or if group-level pre-allocation of SoA buffers gives a major speedup, (b) earns its place.

---

### Q4 — How does `ContactDomain` fit into the analysis flow?

**Where it arises**: [[#5.7]] (contact subsystem).

**The choice**:
- **(a) Side-effect on global force vector**: ContactDomain is queried by the integrator each step; computes contact forces from current positions/velocities; *adds* them to the global $F$ before the integrator does the velocity/displacement update. No tangent stiffness contribution; pure force injection. This is LS-DYNA's pattern (Theory §26 doesn't touch the SOE).
- **(b) Element-equivalent assembly**: each active contact pair contributes a stiffness $K_c$ + force $F_c$ to the global linear system, the same way an Element does. Assembled into the SOE; participates in Newton iteration if any. This is Abaqus/Standard's pattern.

**Pro (a) — side-effect**:
- Trivially works in explicit. No SOE entanglement, no constraint-handler interaction.
- Can use rich, intricate contact algorithms (segment-based, mortar) without modifying the analysis pipeline.
- Failure modes are local — a bad contact stiffness explodes one contact, not the whole step.
- Matches the LS-DYNA reference architecture, so porting algorithms is simpler.

**Con (a)**:
- Doesn't generalize to implicit. If implicit contact ever becomes a target, you need (b) on top.
- Recorders for contact forces need their own subsystem (can't reuse the existing element-force recorder).
- Parallelism: contact-force assembly has to participate in the same domain-decomposition + send/recv protocol as elements; needs care.

**Pro (b) — assembled**:
- Generalizes to implicit. Newton iteration works.
- Reuses existing OpenSees assembly machinery — recorders, parallel decomposition, energy balance hook in for free.
- Cleaner from a "everything is a contribution to the SOE" standpoint.

**Con (b)**:
- Heavier — contact tangent stiffness is non-trivial (and often non-symmetric). LS-DYNA-style segment-based and mortar algorithms don't naturally produce a clean $K_c$.
- Constraint handler interaction is messy. Penalty contact wants to add to the SOE; Lagrange contact wants to add DOFs.
- Implicit contact convergence is a research topic in itself — opening that can has its own cost.

**Hybrid**: a `ContactDomain` that supports both modes, switchable per analysis. Probably overkill for v1; defer.

**Decision direction**: **(a) side-effect for v1** ([[#5.7]]). Matches the explicit-first roadmap. Revisit (b) if and when implicit contact becomes a target — currently it's not.

**Decision criterion to revisit**: if a research project demands implicit contact (e.g., quasi-static foundation rocking with arc-length under monotonic loading), the hybrid mode becomes necessary.

---

### Q5 — Should the modular damage decorator accept any `NDMaterial`?

**Where it arises**: [[#5.4]] (modular damage decorator).

**The technical issue**: GISSMO and similar damage decorators drive damage evolution off the *equivalent plastic strain increment* $\Delta \bar{\varepsilon}^p$ from the wrapped material. Most plasticity models track this internally, but not all expose it via a stable interface.

**The choices**:
- **(i) Strict interface**: only accept `NDMaterial`s that implement a new mixin like `IPlasticityCapable::getPlasticStrainIncrement()`. Compile-time enforced via inheritance or runtime via `dynamic_cast`.
- **(ii) Permissive with runtime check**: accept any `NDMaterial`; query plastic strain via the existing `getResponse(string)` interface (every `NDMaterial` has this); fail with a clear error message at first call if the response isn't there.
- **(iii) Type-tagged registry**: introduce a small enum/set of "damage-capable" tags; new materials declare which they support.

**Pro (i) — strict**:
- Compile-time safety. Can't pass an incompatible material.
- Clean documentation — readers see the interface and know what's required.

**Con (i)**:
- All existing plasticity `NDMaterial`s need to be retrofitted to implement the new interface. Touches dozens of files in `SRC/material/nD/`. Big PR, big merge headache.
- Adding a new damage-driving variable later (e.g., dissipated energy density) means a new interface, repeat the retrofit.

**Pro (ii) — permissive**:
- Zero changes to existing materials. Decorator works with any `NDMaterial` that already exposes plastic strain via `getResponse`.
- Late-binding errors *do* appear at model setup (the decorator queries the response on its first `setTrialStrain` call), so the error is findable, not buried 1000 steps in.

**Con (ii)**:
- "Stringly-typed" API; typos in response names give cryptic errors.
- No compile-time safety.

**Pro (iii) — registry**:
- Middle ground: explicit declaration without modifying every existing material.
- Discoverable: the user can ask "what damage-capable materials exist?"

**Con (iii)**:
- Introduces a registry; one more global table to maintain.
- Materials declaring capabilities they don't actually support (silently broken) is a real risk.

**Decision direction**: **(ii) permissive with runtime check, plus a clear setup-time validation pass**. The decorator queries the wrapped material's response immediately on construction (not on first use); fails loudly with a list of expected response names if the material doesn't provide them.

**Decision criterion to revisit**: if the user-facing error messages turn out to be confusing in practice, switch to (iii) and document a curated list of compatible materials.

---

### Q6 — `getCharacteristicLength()` on Element or Domain?

**Where it arises**: [[#10.5.2]] (regularization framework gap).

**The need**: crack-band regularization ([[#10.3.1]]) requires knowing each element's characteristic length $h_e$ to scale the softening modulus. OpenSees has no systematic API for this today; each element computes (or doesn't compute) its own size for its own purposes.

**The choices**:
- **(a) Element-level**: `virtual double Element::getCharacteristicLength() const`, default implementation returns $\sqrt[3]{V_e}$ for solids, $\sqrt{A_e}$ for shells, length for beams/trusses.
- **(b) Domain-level**: a query like `Domain::getCharacteristicLengthAtNode(nodeTag, radius)` that returns a neighborhood-scale length — needed for nonlocal models that average over a radius.
- **(c) Both**: element-level for the cheap default; domain-level for nonlocal-specific queries.
- **(d) Gauss-point-level**: extend down to per-Gauss-point characteristic length, since elements with anisotropic shape or higher integration order need finer granularity.

**Pro (a) — element-level**:
- Locally computable; the element knows its own geometry.
- Matches LS-DYNA: each element has an implicit characteristic length used by all crack-band-capable materials.
- Default implementation in `Element.cpp` covers 95% of cases; subclasses override only when needed.
- Adoption is incremental — adding the virtual is a 1-file change to `Element.h/cpp`.

**Con (a)**:
- Doesn't help nonlocal models that need a *radius* in physical units, not a per-element length.
- "Characteristic length" is ambiguous for highly elongated or distorted elements. $\sqrt[3]{V_e}$ over-estimates the relevant length for a sliver hex; $V_e / A_{\max}$ (volume divided by max face area) is better but not standard.

**Pro (d) — Gauss-point-level**:
- LS-DYNA's `*MAT_ADD_DAMAGE_GISSMO` actually uses a per-Gauss-point length for elements with non-uniform integration. More correct.

**Con (d)**:
- Requires elements to expose Gauss-point geometry, which most don't today. Bigger framework change.

**Pro (b) — domain-level**:
- Required for nonlocal-integral regularization ([[#10.3.2]]) — needs a radius and a list of points within it.

**Con (b)**:
- More expensive per query (spatial search via bucket sort).
- Not what crack-band actually needs.

**Decision direction**: **(c) start with (a), add (b) when nonlocal models actually require it**.
- Add `Element::getCharacteristicLength()` with a sensible default ($\sqrt[3]{V_e}$ for solids, etc.) as part of [[#5.4]].
- Document the formula used, and let elements override.
- Expose a Gauss-point-level extension as an optional second virtual (`getCharacteristicLengthAtGaussPoint(int)`) that defaults to the element-level value. Materials that need finer granularity (likely just GISSMO and future damage models) call the more specific version.
- Defer the domain-level neighborhood query until nonlocal-integral regularization becomes a priority (likely v3+, after [[#5.7]] builds the bucket-sort infrastructure).

**Decision criterion to revisit**: if the elongated-element issue (sliver hex, high-aspect tet) shows up in benchmark fracture problems, refine the default formula to $V_e / A_{\max}$ or expose multiple length measures.

---

### Q7 — Phase-field as its own roadmap item, merged with XFEM, or unified under cohesive?

**Where it arises**: [[#10.5.6]] (regularization), [[#5.15]] (known gaps).

**The technical situation**: discrete-fracture modeling has three families that share the goal but use completely different machinery:
- **Cohesive elements** (already in [[#5.15]] as a future item): element-level traction-separation laws on pre-defined interfaces. Most established; most engineering-friendly.
- **XFEM / GFEM**: extends the FE basis with enrichment functions at known crack locations. Cracks are discrete; mesh-cutting and crack-tip enrichment needed. Implementation lives at the element level (modified shape functions).
- **Phase-field fracture** (Francfort-Marigo, Miehe et al.): a smeared damage field $d(\mathbf{x})$ as an extra global PDE. Cracks emerge automatically; branching and merging come for free. Implementation is domain-level (a new analysis mode with a coupled damage-mechanics solver).

**The choice**:
- **(α) One roadmap item: "Discrete fracture methods"** combining XFEM and phase-field as a single research direction, with cohesive as a foundation.
- **(β) Two separate roadmap items** (XFEM and phase-field), with cohesive as a separate prerequisite.
- **(γ) Phase-field never on the roadmap** — defer indefinitely; XFEM stays in [[#5.15]] as a "known gap."

**Pro (α) — unified**:
- Looks coherent on paper. "We support discrete fracture" is a single capability claim.
- Allows shared infrastructure decisions (e.g., the cohesive-element foundation).

**Con (α)**:
- The implementations are nothing alike. XFEM = element-level enrichment, modified shape functions, level-set tracking. Phase-field = global gradient PDE, staggered solver, fine meshes. Bundling them produces a roadmap item that would actually have to be split during implementation anyway.
- One "discrete fracture" plan would either be vague (helpful to nobody) or too long (unread).

**Pro (β) — separate items**:
- Each method gets its own focused per-feature plan.
- Different research interests can pursue them independently.
- Prerequisites stay distinct: XFEM needs level-set / mesh-cutting infrastructure; phase-field needs a coupled-PDE solver.

**Con (β)**:
- Two extra roadmap items inflates the plan further.

**Pro (γ) — phase-field deferred**:
- Phase-field is research-grade, expensive, and only really clean for monotonic loading. Cyclic seismic problems are not its strength.
- Cohesive elements (already [[#5.15]] item 1) cover most engineering fracture needs.

**Con (γ)**:
- Phase-field is the modern academic mainstream for fracture mechanics. Deferring it limits research-paper alignment with current literature.
- The user explicitly flagged phase-field-adjacent interest (XFEM mention earlier).

**Decision direction**: **(β) — separate items, when each is ready to be promoted from [[#5.15]]**. Specifically:
- **Cohesive elements stay where they are** in [[#5.15]] as the highest-priority fracture follow-up. They're the foundation for both XFEM and phase-field eventually (XFEM uses cohesive surfaces internally for the enrichment; phase-field reduces to cohesive in the limit $\ell_0 \to 0$).
- **XFEM and phase-field each become their own future roadmap items** when promoted, with their distinct prerequisites called out.
- For the phase-field item specifically, the prerequisites are: (i) cohesive elements working ([[#5.15]] item 1), (ii) a coupled-PDE solver mode in OpenSees (which doesn't exist today and would be its own infrastructure work — possibly the gradient-enhanced damage path from [[#10.3.3]] is the lighter precursor).

**Decision criterion to revisit**: if a specific research project (a thesis, a paper, an industry collaboration) wants phase-field specifically, promote it. Otherwise let it sit until cohesive ships and the gradient-enhanced damage question becomes concrete.

---

### Question lifecycle

When a question gets resolved (the per-feature work makes the call clear, or a benchmark settles it), it moves out of this section and into the relevant per-feature plan's *Implementation log*, with the date and the resolving evidence. The point of this section is *open* questions — closed ones belong with the work that closed them.

> [!note]
> If a new architectural question surfaces during per-feature work, add it here with a backlink to where it arose. Don't bury it in an issue tracker that nobody will re-read.

---

## 7. Reference: LS-DYNA Theory Manual sections this roadmap leans on

PDF page numbers from `LS-DYNA Theory.pdf` (R-version embedded in file). Cross-check via `~/.claude/skills/ls-dyna/references/manual_index.md`.

| Topic | Section | PDF pages |
|---|---|---|
| Volume integration / 1-pt vs full int | §3.1 | 39 |
| Flanagan-Belytschko hourglass | §3.2 | 40–46 |
| Puso assumed-strain hourglass | §3.3 | 46–52 |
| Belytschko-Tsay shell, co-rotational | §7 | 111–122 |
| Hughes-Liu shell | §10 | 143–158 |
| Stress update overview / Jaumann | §18 | 243–258 |
| Timestep formulas | §22 | 489–493 |
| Central difference + stability + subcycling | §24 | 499–509 |
| Rigid body dynamics + joints | §25 | 511–521 |
| Contact algorithm (full) | §26 | 523–568 |
| Bucket sort | §26.11 | 545–550 |
| Vectorization (block-of-128 architecture) | §29.1 | 575–577 |
| Parallelization (shared memory) | §29.2 | 578–580 |
| SPH | §38 | 637–643 |
| EFG | §39 | 645–657 |

LS-DYNA Vol II (`*MAT_*` and `*EOS_*`) and Vol I (keyword cards) are referenced inline in each per-item plan.

---

## 8. Tradeoffs and design decisions

This section captures the architectural and engineering tradeoffs we'll keep running into as the roadmap unfolds. Most are not single-decision moments — they're tensions that resurface in every per-feature plan, and being explicit about them now saves rehashing each time.

### 8.1 Polymorphism vs raw throughput

**Tension**: OpenSees' per-element virtual dispatch is what makes it a *framework* — anyone can drop in a new element and the analysis just works. But that same dispatch is the thing that costs performance in the explicit element loop. LS-DYNA gives up per-element polymorphism (block-of-128 SoA dispatch — Theory §29.1) to win throughput.

**Decision**: keep polymorphism at the user-facing API (Tcl / Python / `Element*` is still a virtual base) and recover performance through the opt-in batch interface ([[#5.2]]). New elements that opt in get LS-DYNA-class throughput; existing elements pay zero cost and keep working unchanged.

**What we accept**: a two-tier element ecosystem. Some elements are "fast path" (SoA-native, batch-aware), others are "slow path" (per-element, the existing convention). Documenting which is which becomes part of the per-element README.

**What we don't compromise on**: the user shouldn't have to know which tier their element is in. The integrator picks automatically based on `getClassTag()`.

### 8.2 Determinism vs parallel scaling

**Tension**: OpenSees insists on bit-exact reproducibility — serial *and* parallel. Every component implements `sendSelf/recvSelf`; the assembly is deterministic; rerunning the same model on the same hardware gives identical results. LS-DYNA explicitly gives this up (Theory §29.2, p. 579: "the order of operations will vary from run to run … variations in nodal accelerations and sometimes even velocities are observable") for SMP/MPP performance.

**Decision**: don't give up determinism. Graph-colored scatter ([[#5.11]]) preserves it within shared-memory; domain-decomposed `sendSelf/recvSelf` preserves it across MPI ranks. Pay the constant-factor cost; the value for research and regression testing is irreplaceable.

**Future GPU consideration**: deterministic atomics on GPU are slow. When/if Ladruno targets GPU, this tradeoff comes back — at minimum, a per-run "deterministic / fast" flag like LS-DYNA's "ordered summation" option. Defer the decision; flag it as a future research call.

### 8.3 UL+Jaumann vs TL kinematic convention

**Tension**: discussed at length in [[#5.10]]. UL+Jaumann is faster and SoA-friendly, matching LS-DYNA's default for solids; but it has known pathologies (Jaumann oscillations under simple shear of hyperelastic rubber — Theory §18.1, p. 243) and stores per-element reference state implicitly. TL is what OpenSees does today, natively handles hyperelasticity, and matches the framework's "store reference, compute current" data layout — but is heavier per element and harder to SoA-batch.

**Decision**: support both per-element. New explicit solids default to UL+Jaumann ([[#5.3]]); hyperelastic / SPH-FORM=7/8 elements use TL. The framework doesn't pick a convention — each formulation declares its own.

**Risk**: a model can mix elements with different conventions, and care must be taken at coupling interfaces (e.g., contact between a UL hex and a TL hyperelastic block). Document this; add diagnostics that warn at convention boundaries.

### 8.4 Mass scaling artifacts vs explicit tractability

**Tension**: without mass scaling, the smallest element bounds the global timestep. With mass scaling, low-frequency response shifts because added inertia changes effective frequencies. Selective mass scaling ([[#5.1]]) is the practical compromise.

**Decision**: implement selective mass scaling with strong diagnostics:
- Track and report total added mass as % of model mass.
- Hard-warn at >5%, error at >10% by default (configurable).
- Report which nodes received added mass and how much.
- Compare modal frequencies pre- and post-scaling (cheap eigenproblem on the relevant DOFs) when the user asks.

**What we accept**: the user must validate that their problem is rate-insensitive enough for mass scaling to be physical. We give them tools to check; we don't pretend the artifacts don't exist.

**Specific to seismic**: aggressive mass scaling can shift the structural fundamental period — fatal for seismic response. Document a workflow: run the modal first, set $\Delta t_{\text{target}}$ to keep mode 1 frequency error below 1%.

### 8.5 Penalty vs Lagrange constraint enforcement

**Tension**: contact and joints both need a constraint enforcement strategy. Penalty is fast and simple but has tunable stiffness with no physical meaning — and penetration depends on stiffness. Lagrange is exact but adds DOFs and is hard for explicit. Augmented Lagrangian is the middle ground, popular in implicit but rarely used in explicit.

**Decision**:
- **Contact** ([[#5.7]]) defaults to penalty (LS-DYNA's standard `SOFT=0/1/2` family). Mortar / augmented Lagrangian as a v3 if implicit contact ever becomes a target.
- **Joints** ([[#5.6]]) leverage OpenSees' existing `ConstraintHandler` swap — the user picks penalty / Lagrange / transformation per analysis. **This is an architectural win over LS-DYNA**, where joint enforcement is hardcoded.

**What we accept**: penalty stiffness becomes a user knob. Document the LS-DYNA defaults (`SLSFAC = 0.10`, bulk-modulus·area²/volume scaling — Theory §26.7.1, Eq. 26.14) as the starting point.

### 8.6 Implicit vs explicit at the user level

**Tension**: which solver should a user pick for which problem? The answer is "it depends," and we want to make the depends-on explicit so users don't pick wrong.

**Heuristic table** (will be documented for users):

| Problem class | Recommended | Reason |
|---|---|---|
| Static pushover, modal, eigenvalue | Implicit | Equilibrium-finding is the goal; explicit doesn't help |
| Low-frequency seismic, well-behaved nonlinear | Implicit | Newton converges; large step sizes possible |
| High-rate dynamic (impact, blast) | Explicit | CFL step ≈ accuracy step; no convergence concerns |
| Severe softening / damage / fracture | Explicit | Newton fights tangent singularities |
| Contact-dominated (pounding, rocking) | Explicit | Explicit penalty contact is robust by construction |
| Liquefaction post-triggering | Explicit | State changes kill implicit |
| Multi-body with joints | Either | Pick per problem; rigid-body cost is similar |
| Quasi-static large deformation (forming) | Either | Explicit with mass scaling is common; implicit with arc-length also works |

**Decision**: provide both, document when each is appropriate, don't deprecate implicit. The default suggestion in tutorials: "implicit unless you have a reason."

**Bridge case — implicit-explicit subcycling** (LS-DYNA Theory §24.4, p. 504): out of scope for v1. Revisit if specific problems demand it.

### 8.7 RigidBody as DomainComponent vs special Element

**Tension**: covered in [[#5.5]] and Section 6 Q1. Architecturally, rigid bodies are *not* elements — they don't compute strain, they don't have a constitutive model, they have rotational DOFs that don't fit the standard nodal-translation pattern. But the OpenSees element-loop infrastructure is heavily wired up; piggybacking on it is the path of least resistance.

**Decision**: pay the architectural cost upfront. Make `RigidBody` a parallel `DomainComponent` kind. Reasoning: every code path that iterates elements (recorders, parallel decomposition, output) would have to learn to skip rigid bodies if we made them special elements — that's debt that compounds, vs paying once now for a clean abstraction.

**Risk**: more framework changes. The element-loop infrastructure has to gain awareness of a new component kind. This is acknowledged as the heaviest item on the roadmap.

### 8.8 Build complexity vs capability

**Tension**: every new feature adds compile time, dependencies, source-file count, and maintenance burden. Ladruno already has a non-trivial compilation story (see [[../Ladruno_internal/01_compilation_journal|the journal]]).

**Decision**: every Ladruno feature lives behind a CMake flag, off-by-default in `Conf.cmake`. Users who don't need contact don't pay for compiling it. Mass scaling and the batch interface are universal (no flag); rigid bodies, joints, contact, SPH, EFG each get their own flag.

**Maintenance commitment**: each merged feature gets a smoke test in `EXAMPLES/Ladruno/`. A feature without a test ages into a feature that nobody trusts.

### 8.9 Research-paper churn vs maintainable contribution

**Tension**: every item on this roadmap could be shipped as a one-off paper artifact (private fork, dies on next OpenSees merge) or as a maintainable upstream contribution (takes 2-3× longer, but compounds). The opportunity cost runs in both directions.

**Decision**: maintainable, by default. Each feature graduates through:
1. Per-feature plan in this folder.
2. Working prototype in a Ladruno branch.
3. Tests + docs.
4. Upstream PR to OpenSees/OpenSees when stable, with the Ladruno fork keeping in-flight work.

**What we accept**: slower paper output per researcher, but the platform compounds — features stay alive, can be combined, get used by others. Over a 3-5 year horizon this is the better trade for the group.

**Exception**: throwaway research prototypes (probing whether an approach works at all) live outside this folder, in `Ladruno_internal/spikes/` or similar. Once an approach is validated, it graduates to a per-feature plan here.

### 8.10 The "explicit is easier" misconception

**Tension**: there's a common belief that explicit is "easier" than implicit because there's no Newton iteration to converge. **This is false** and it's worth stating clearly so future-us doesn't fall into it.

Explicit has its own discipline:
- **Energy balance** — hourglass energy, contact energy, mass-scaled inertia work, sliding-interface energy must all stay small. Without monitoring, a run can be silently wrong.
- **Hourglass control** — under-integrated elements need stabilization that affects accuracy. The Puso vs Flanagan-Belytschko vs Belytschko-Bindeman question is research in itself (LS-DYNA Theory §3.2-3.3).
- **Mass scaling judgment** — see [[#8.4]].
- **Contact stiffness tuning** — penalty stiffness is a user knob.
- **CFL discipline** — material wave speed has to be tracked carefully; rate-dependent softening can drop wave speed mid-run, requiring step recomputation.
- **Convergence is silent** — there's no Newton residual to alert you that something's wrong. You only find out from physical sanity checks (energy, momentum, displacement bounds).

**Decision**: build diagnostics first-class. Energy-balance recorder, hourglass-energy monitor, added-mass tracker. The Ladruno explicit story is *explicit dynamics with discipline*, not *explicit dynamics as the easy option*.

---

## 9. Deep dive: block dispatch explained

This is the spine idea of the entire performance argument in [[#5.2]]. Worth understanding properly because every downstream item assumes it. We'll build it up from first principles, ground it in actual code, and finish with order-of-magnitude numbers.

### 9.1 The two element loops, side by side

The simplest framing is to look at what each framework does inside the inner loop of an explicit dynamic step.

**OpenSees today** (paraphrased from `SRC/analysis/integrator/CentralDifference.cpp`):

```cpp
// One explicit step: compute internal forces from every element
for (Element* e : domain.getElements()) {
    e->update();                          // [virtual call]
    const Vector& F_internal = e->getResistingForce();  // [virtual call]
    assembleIntoGlobalForce(e, F_internal);
}
```

**LS-DYNA** (paraphrased from Theory §29.1, p. 575):

```fortran
! Process LLT-LFT+1 elements (e.g. 128) in one shot
DO I = LFT, LLT
   X1(I) = X(1, IX1(I))            ! gather node 1 x-coord into local vec
   Y1(I) = X(2, IX1(I))
   Z1(I) = X(3, IX1(I))
   VX1(I) = V(1, IX1(I))
   ! ... gather all 8 nodes' positions and velocities
ENDDO

! ... compute strain rate from velocities, all 128 at once ...
! ... call material kernel with 128-wide stress arrays ...
! ... compute forces, all 128 at once ...

DO I = LFT, LLT                    ! scatter forces back
   RHS(:, IX1(I)) = RHS(:, IX1(I)) + FORCE(:, 1, I)
   ! ... 8 nodes per element ...
ENDDO
```

> [!important] The key visible difference
> OpenSees does `e->update()` *N* times — one virtual call per element. LS-DYNA does it *once* for all 128 elements together. **Same physics, completely different inner-loop shape.**

### 9.2 Why per-element dispatch is expensive — three layers

Many engineers' first reaction is "it's just a function call, how slow can it be?" The answer is "the function call itself is fast; everything *around* the function call is slow." Three separate effects compound:

#### 9.2.1 The virtual call itself

```cpp
e->update();
```

In C++, when `update()` is virtual, this becomes:

```
mov  rax, [rdi]         ; load vtable pointer from the object
mov  rax, [rax + 0x40]  ; load function pointer from vtable slot
call rax                ; indirect call
```

Three memory loads + one indirect call. On modern x86, **~5–10 ns** in the best case (everything in L1 cache). Not catastrophic per call, but if you have $10^6$ elements and $10^5$ time steps, that's $10^{11}$ calls — **~500–1000 s of pure dispatch overhead** before any real work.

#### 9.2.2 The indirect branch mispredict — bigger than it looks

Modern CPUs are speculative: they guess where each branch is going and start executing ahead. Direct calls (compile-time known target) are predicted near-perfectly. **Indirect calls** through a vtable are the hardest case — the CPU's branch predictor has to guess based on history.

In a model with one element type, the predictor gets it. In a heterogeneous model with 5 element types interleaved, the predictor misses constantly. Each mispredict is **~15–25 cycles wasted** flushing the pipeline. At 4 GHz, that's ~5 ns *per mispredict*, on top of the call cost.

But here's the worse case: even with one element type, the *function called by* `update()` may itself contain virtual calls (to the material). So you've got `Element::update` → `NDMaterial::setTrialStrain` → `NDMaterial::getStress` — **three indirect calls per element**, each a potential mispredict.

#### 9.2.3 Cache and prefetch — usually the dominant cost

This is the one most people don't think about, and it's often **10× larger than the function-call cost**.

When you do `e->update()` per element, you're chasing pointers: `e` → its node pointers → node coordinates → its material → material state. Each pointer chase is potentially a cache miss. A cache miss to L3 is ~30 ns; a miss all the way to main memory is ~100 ns.

The CPU's hardware prefetcher cannot help you here. Prefetchers detect *strided* access patterns ("the program is reading addresses A, A+64, A+128, …") and pre-load the next line. They cannot predict where `e` points to next, because the pointer chain is data-dependent.

**The result**: the dominant cost in a per-element loop is often *waiting for memory*, not actual arithmetic.

### 9.3 Memory layout: AoS vs SoA

To understand why LS-DYNA is fast, you have to understand how it lays out data.

**Array-of-Structs (AoS)** — what OpenSees naturally does:

```
Memory:  [Elem0: x,y,z,vx,vy,vz,...] [Elem1: x,y,z,vx,vy,vz,...] [Elem2: ...]
                                     ^                          ^
                                     cache line boundary        cache line boundary
```

To process all the x-coordinates of element 0's first node, you load a cache line that *also* contains element 0's velocity, material state, history variables, etc. — most of which you don't need yet. Cache utilization: poor.

**Struct-of-Arrays (SoA)** — what LS-DYNA does inside its element kernels:

```
Memory:  X1: [x_e0, x_e1, x_e2, x_e3, x_e4, ..., x_e127]   ← contiguous, all hot
         Y1: [y_e0, y_e1, ..., y_e127]
         Z1: [z_e0, z_e1, ..., z_e127]
         VX1: [vx_e0, ...]
         ...
```

Now, when you process the x-coords of node 1 across all 128 elements, you read a perfectly contiguous chunk. The hardware prefetcher locks on instantly. SIMD instructions (AVX-512: 8 doubles in one register) eat through this layout natively — the compiler vectorizes the loop with no help.

That `DO I = LFT, LLT` Fortran loop in §9.1 isn't slow Fortran — it's the *fastest possible* memory access pattern on modern hardware. Each iteration of the inner loop touches one new cache line.

> [!info] The "blocks of 128" choice
> LS-DYNA picks 128 because:
> - It's a multiple of 64 (the typical vector-register-line size on Cray-1, where this design crystallized; today it works because 128 doubles = 8 AVX-512 registers fully loaded).
> - It's small enough that the local SoA buffers fit comfortably in L1 cache (~32 KB on modern x86).
> - It's large enough to amortize the gather cost.
> Theory §29.1 (p. 575) is candid: "Larger groups give a marginally faster code, but can reduce computer time sharing efficiency because of increased core requirements."

### 9.4 What batch dispatch actually does — the gather/kernel/scatter pattern

Now we can describe what `batchAssembleResisting` ([[#5.2]]) does mechanically. The pattern is **gather → kernel → scatter**:

```cpp
// What an opt-in element override looks like — pseudocode
int MyExplicitHex::batchAssembleResisting(
    const std::vector<Element*>& batch,
    Vector& globalForce,
    const AnalysisModel& model)
{
    const int N = batch.size();   // e.g. 128

    // ---- GATHER phase ----
    // Allocate local SoA buffers (or reuse static ones)
    static thread_local std::vector<double> X1(N), Y1(N), Z1(N), VX1(N), VY1(N), VZ1(N);
    // ... 8 nodes worth of these ...

    for (int i = 0; i < N; ++i) {
        MyExplicitHex* e = static_cast<MyExplicitHex*>(batch[i]);
        // Pull node 1 of element i into the SoA buffers
        const double* coords = e->getNodeCoords(0);   // direct, non-virtual
        X1[i] = coords[0];  Y1[i] = coords[1];  Z1[i] = coords[2];
        // velocities, history vars, ...
    }

    // ---- KERNEL phase ----
    // Now we run the element formulation 128-wide, no virtual calls,
    // perfect cache locality, compiler auto-vectorizes
    for (int i = 0; i < N; ++i) {
        // strain rate from velocities and shape function derivs
        // call material at gauss point — also batched
        // accumulate internal force into a local FORCE[24][N] buffer
    }

    // ---- SCATTER phase ----
    // Push forces back into the global force vector
    for (int i = 0; i < N; ++i) {
        // ... add FORCE[:, i] to globalForce at the element's DOF tags
    }

    return 0;
}
```

The integrator side is dead simple:

```cpp
// CentralDifference::computeInternalForces() — paraphrased
std::map<int, std::vector<Element*>> byType;
for (Element* e : *domainPtr) {
    byType[e->getClassTag()].push_back(e);  // bucket by type tag (cheap)
}

// Each bucket gets ONE call (which dispatches internally to all its members)
for (auto& [tag, batch] : byType) {
    batch[0]->batchAssembleResisting(batch, globalForce, *theModel);
    // ↑ ONE virtual call for the whole batch, not one per element
}
```

> [!note] The crucial observation
> The virtual call still happens — once per element type, not once per element. With 5 element types and $10^6$ elements, you go from $10^6$ virtual calls per step to **5 virtual calls per step**. Five orders of magnitude fewer.

### 9.5 The opt-in default fallback

The non-invasive part is the base-class default:

```cpp
// In Element.cpp — default fallback for elements that don't override
int Element::batchAssembleResisting(
    const std::vector<Element*>& batch,
    Vector& globalForce,
    const AnalysisModel& model)
{
    // Old per-element loop. No SoA, no batching, no speedup.
    // But also: no breakage of any existing element.
    for (Element* e : batch) {
        e->update();
        const Vector& F = e->getResistingForce();
        assembleIntoGlobalForce(e, F, globalForce);
    }
    return 0;
}
```

This is what makes the design non-invasive: **every existing OpenSees element keeps working with zero modification**. They go through the slow path. Only elements written for the explicit roadmap (starting with [[#5.3]] — the SoA-native hex) override the default to opt into the fast path.

A user model that mixes one fast-path element (5000 SoA hexes) with one slow-path element (200 legacy beams) gets the fast path on the hexes and the slow path on the beams — automatic, no user intervention.

### 9.6 The scatter problem — why parallel scatter is hard

The gather is trivially parallel — each element reads from independent positions. The kernel is trivially parallel — each element computes independently. **The scatter is not.**

The reason: **shared nodes**. If element 5 and element 17 both touch node 42, and two threads run their scatter steps simultaneously:

```
Thread A (elem 5):   tmp = globalForce[42];  tmp += F5_at_42;  globalForce[42] = tmp;
Thread B (elem 17):  tmp = globalForce[42];  tmp += F17_at_42; globalForce[42] = tmp;
```

If A and B interleave, you lose one of the contributions — the classic read-modify-write race. Three solutions, in order of complexity:

1. **Serial scatter.** Single-threaded loop after the parallel kernel. Cheap and correct, but caps your speedup. *This is what [[#5.2]] does in v1.*
2. **Atomic accumulation.** Use `atomic_add` for each scatter. Correct but slow (atomics serialize through cache coherency).
3. **Disjoint coloring.** Pre-process the elements into groups (colors) such that no two elements in the same color share a node. Each color is then scattered in parallel; colors process serially.

LS-DYNA Theory §29.1 (p. 577, Fig. 29.1) shows option 3: "Group of 48 elements broken into 4 disjoint blocks." Each block can scatter in parallel because the blocks share no nodes.

This is what [[#5.11]] does — graph-color the elements into disjoint sets, then dispatch each color in parallel. Worth it only if benchmarks show the serial scatter is the bottleneck. **Probably worth it for any model with $>10^5$ elements running on $>4$ cores.**

### 9.7 Why this only matters in explicit, not implicit

Implicit Newton step structure:

```
1. Element loop: assemble K (tangent stiffness) and R (residual force)    ← O(N) work
2. Solve K·du = R for du                                                  ← O(N^{1.5}) to O(N^3)
3. Convergence check
4. Update u
5. Repeat until converged (multiple Newton iterations per step)
```

For implicit, step 2 — the linear solve — typically dominates the runtime. Even on a sparse direct solver, it's super-linear in problem size, while the element assembly is linear. By the time you have $>10^5$ DOFs, the solve is 80%+ of the time. The element loop is in the noise; per-element virtual dispatch overhead is invisible.

Explicit step structure:

```
1. Element loop: assemble F (internal force vector only — no K!)          ← O(N) work
2. Update u, v, a using lumped mass:  u_{n+1} = ... cheap arithmetic ...   ← O(N) work
3. Boundary conditions
4. (no convergence test, no iteration)
```

There's no global linear solve. The whole step is the element loop plus some cheap vector arithmetic. **The element loop is the *entire* workload.** Every nanosecond per element matters because there's nothing else competing for the time budget.

This is why the batch-dispatch idea pays off in explicit and not in implicit. Architecturally, there's no reason you couldn't add a batch path for implicit too — but the gains would be tiny because you'd still be dominated by the linear solve.

### 9.8 Order-of-magnitude — what we're actually buying

Concrete back-of-envelope for a representative explicit problem:

**Setup**: $10^6$ hex elements, $10^5$ time steps, single thread, single element type.

**OpenSees today** (per element, per step):
- Virtual call overhead: ~5 ns
- Indirect branch mispredict (occasionally): ~5 ns
- Pointer chase to material + node data: ~30–100 ns (depends on cache state)
- Actual element work (small — a hex with elastic material is maybe 2 μs): ~2000 ns
- **Total per element-step**: ~2050 ns, of which ~50 ns is dispatch overhead (~2.5%)

For one element type and good cache behavior, dispatch is small. **But:**

**Realistic case** (multiple element types, larger material state, cold caches between elements):
- Dispatch overhead inflates to ~200–500 ns per element-step
- Element work doesn't change
- **Total per element-step**: ~2500 ns, of which ~500 ns is dispatch (~20%)

At $10^6 \times 10^5 = 10^{11}$ element-steps, that 20% dispatch overhead is **~5–10 hours of pure pipeline waste** in a 24-hour analysis.

**OpenSees with batch dispatch** (one virtual call per type per step, SoA-native kernels):
- Dispatch becomes negligible (5 calls per step instead of $10^6$)
- Memory access is contiguous SoA, prefetcher works, no cache misses inside the kernel
- SIMD vectorization kicks in — element kernel work *itself* can be 2-4× faster
- **Estimated total per element-step**: ~700–1000 ns

Net realistic speedup: **~2.5–3.5×** end-to-end on a representative model. Not 10×; not 100×. But this is on the *element loop*, which IS the explicit workload — so it translates directly to wall-clock.

LS-DYNA gets a slightly higher absolute number because their element kernels were hand-vectorized for Cray-1 in 1978 and have been polished for 47 years. We don't need to hit that ceiling; we need to close the gap from "10× slower than LS-DYNA" to "1.5–2× slower than LS-DYNA," which is plenty for the platform argument to work.

### 9.9 What this gives you mentally

The way to think about the design: **every level of the call stack should know about the same number of elements at once**.

- The user's model defines $10^6$ elements.
- The integrator sees $10^6$ elements organized into ~5 type-buckets.
- The element class's `batchAssembleResisting` sees a batch of ~$2 \times 10^5$ elements of one type.
- That method internally processes them in chunks of ~128 (to fit local SoA buffers in L1 cache).
- The CPU's SIMD units process 8 doubles at a time inside the chunk.

When the levels are *aligned* like this, the data flows smoothly through the cache hierarchy. When the levels are *mis-aligned* (per-element dispatch in a $10^6$-element problem), the data thrashes through the cache hierarchy on every call.

OpenSees today is mis-aligned; LS-DYNA is aligned. The batch interface is what realigns OpenSees, **without breaking anything that's already working**.

> [!summary] The one-sentence version
> Per-element virtual dispatch isn't slow because of the function call — it's slow because *the cache hierarchy can't predict what comes next*. Batch dispatch fixes that by making the next thing predictable: more of the same.

---

## 10. Deep dive: localization and regularization

The capability claims in [[#1.1]] about concrete fracture, liquefaction post-triggering, and buckling-into-collapse all rest on softening constitutive models. **Softening BVPs are mathematically ill-posed without regularization**, and the consequences of ignoring this are not "small numerical errors" — they're "the answer changes when you refine the mesh, in ways that look superficially reasonable but are physical garbage." This section makes the regularization requirement a first-class part of the roadmap, not a footnote.

### 10.1 The localization problem in plain language

Take the simplest possible case: a 1D bar in tension with a stress-strain law that softens past peak.

```
σ ↑
  |   ___
  |  /   \___
  | /        \___ ← softening branch (negative slope)
  |/             \___
  +-----------------→  ε
```

You apply prescribed displacement at the ends and pull. As long as everything is hardening, every element strains uniformly and life is good. Past peak — when one element happens to be slightly weaker, or the mesh seeds a tiny imperfection — that element keeps straining while its neighbors *unload elastically*. Strain runs away in one element; the rest of the bar springs back.

This is **localization**: the deformation concentrates into a thin band. Now run the same problem with a finer mesh:

| Mesh | Localization band width | Energy dissipated to fracture |
|---|---|---|
| 4 elements | 1/4 of the bar | 25% of original |
| 40 elements | 1/40 of the bar | 2.5% |
| 400 elements | 1/400 of the bar | 0.25% |
| 4000 elements | 1/4000 | 0.025% |

> [!warning] The pathology
> The localization band is always **one element wide**. Refining the mesh shrinks the band proportionally, and the dissipated fracture energy goes to zero as $h \to 0$. The mesh-refined model isn't more accurate — **it's a fundamentally different problem with less energy dissipation, less ductility, more brittle response, sharper load drop**. This is not a numerical bug; the FEM is solving exactly what you asked. The continuum BVP is what's wrong.

Mathematically: when the constitutive tangent has negative eigenvalues, the governing PDE loses ellipticity (Hadamard 1903). The solution becomes non-unique; the localization band can be anywhere; in the continuum limit its width is zero. **Real materials don't have zero-width fracture bands** — concrete cracks are a few aggregate diameters wide, soil shear bands are 10–20 grain diameters wide, metal necking is set by void distribution. Real materials have a length scale; the local-continuum theory does not.

**Regularization** = injecting a length scale back into the model so the BVP becomes well-posed and the localization band has finite width regardless of mesh.

### 10.2 Why this matters for our problem space

Every capability claim in [[#1.1]] that involves softening hits this:

- **Concrete fracture (CDPM, CSCM, KCC)** — the canonical case. Without regularization, the load-displacement curve is mesh-dependent, the dissipated energy goes to zero, and the failure pattern is mesh-aligned (cracks follow element boundaries instead of physical directions).
- **Liquefaction post-triggering** — once PM4Sand softens past the residual ratio, strain localizes into a shear band. Lateral-spreading runout depends on the band width, which depends on regularization.
- **Steel necking / fracture** — Johnson-Cook with damage softens past peak; tube buckling and folding patterns depend on the localization band.
- **RC column folding under axial+lateral** — concrete crushes (softens), longitudinal bars buckle (softens through P-Δ + yielding), the failure mode depends on the order in which softening initiates and how it localizes.
- **Slope failure / landslide initiation** (relevant for the SPH push, [[#5.8]]) — strain-softening clay or sand develops a slip surface whose width is set by regularization.

A roadmap that promises these capabilities but punts on regularization is selling capabilities it doesn't have. **Every softening-capable item must specify its regularization strategy in its per-feature plan.**

### 10.3 The regularization toolbox

Six families, in order of how often they appear in practice for our domain.

#### 10.3.1 Crack-band / characteristic-length scaling (Bažant 1976, 1983)

The pragmatist's choice. Used in **every major concrete model in LS-DYNA** (CSCM, CDPM, KCC, MAT_72R3, Winfrith) and most concrete models everywhere.

**Idea**: parameterize the softening branch by the *fracture energy* $G_f$ (energy per unit crack area, a real material property) instead of by the ultimate strain. Then scale the softening modulus per element so that the energy dissipated by an element going from peak to zero stress equals $G_f / h_e$, where $h_e$ is a characteristic length of the element (typically $\sqrt[3]{V_e}$ for hexes, $\sqrt{A_e}$ for shells / 2D).

Concretely, for a linear-softening damage law $\sigma = f_t (1 - D)$ with $D$ evolving from 0 to 1 over a strain range $[\varepsilon_p, \varepsilon_u]$:

$$\varepsilon_u = \varepsilon_p + \frac{2 G_f}{f_t \, h_e}$$

Element gets brittle at small $h_e$, ductile at large $h_e$, but **dissipated energy per unit fracture area is constant**.

**What this fixes**: mesh-independence of total dissipated energy; load-displacement curves converge as the mesh refines.

**What this does NOT fix**: the localization band is still one element wide. Below an element size of ~$G_f / f_t$ the fracture energy needed exceeds the element's capacity and the model gives up. Doesn't fix mesh-orientation bias (cracks still want to align with element edges).

**Why it's the default anyway**: cheap (just a per-element scaling), no global coupling, easy to add to any softening model. Good enough for engineering when the mesh is reasonable (element size ≈ aggregate size for concrete).

#### 10.3.2 Nonlocal integral (Pijaudier-Cabot & Bažant 1987)

Replace the local damage-driving variable (plastic strain, equivalent strain, etc.) with a **weighted spatial average** over a neighborhood:

$$\bar{\varepsilon}^p(\mathbf{x}) = \frac{\int_\Omega w(\|\mathbf{x} - \boldsymbol{\xi}\|) \, \varepsilon^p(\boldsymbol{\xi}) \, dV}{\int_\Omega w(\|\mathbf{x} - \boldsymbol{\xi}\|) \, dV}$$

where $w$ is a bell-shaped kernel (Gaussian, polynomial bump) with characteristic radius $\ell$. The damage / yield evolution then uses $\bar{\varepsilon}^p$ instead of the local $\varepsilon^p$.

**Effect**: a true material length scale $\ell$ enters the model. Localization bands have a finite width $\sim 2$–$3\ell$ regardless of mesh. Hadamard ellipticity is restored.

**Cost**:
- Each Gauss point needs a list of neighbors within radius $\ell$ — a spatial query (bucket sort comes back here, same machinery as [[#5.7]]).
- The integral averaging is a non-local operation — hard to parallelize across MPI ranks; even shared-memory is tricky because of the read-write dependency between elements.
- LS-DYNA exposes this via `*MAT_NONLOCAL` (Vol II p. 202), which wraps any underlying material. The manual is candid that it's expensive — "use only where needed."

**When it's worth it**: when crack-band isn't enough (very fine meshes; cases where you need the band *width* to be physical, not just the energy).

#### 10.3.3 Gradient-enhanced damage (de Borst, Peerlings, Aifantis)

Add a Laplacian term to the evolution of the nonlocal variable:

$$\bar{\varepsilon}^p - \ell^2 \nabla^2 \bar{\varepsilon}^p = \varepsilon^p$$

Mathematically equivalent to a Helmholtz-kernel nonlocal model (the Green's function of the Helmholtz operator is exactly the bell kernel of 10.3.2), but solved as a **coupled PDE** rather than a spatial integral.

**Effect**: same as 10.3.2 — proper length scale, finite-width band.

**Cost**: adds a DOF per node (the nonlocal variable $\bar{\varepsilon}^p$). Requires solving a coupled system; in implicit, this is just an extra Newton block. In explicit it's harder — typically you do an inner iterative solve on the gradient PDE every step, or you treat it explicitly with its own CFL.

**Status in LS-DYNA**: not directly supported. (Phase-field, which is gradient-enhanced damage in disguise, is also not in LS-DYNA.)

#### 10.3.4 Viscous / rate-dependent regularization — the explicit-dynamics gift

Any positive rate dependence in the constitutive law introduces a length scale through

$$\ell \sim c \cdot \tau$$

where $c$ is the elastic wave speed and $\tau$ is the relaxation time of the rate dependence (e.g., the Cowper-Symonds parameter $1/p$, or the Johnson-Cook reference strain rate inverse).

**Effect**: the localization band cannot shrink below $\ell$ because faster-than-$c$ propagation is inadmissible. Even a small amount of rate dependence stabilizes the BVP under dynamic loading.

> [!important] This is a major argument for explicit + rate-dependent materials
> If you run softening concrete or softening soil under dynamic loading with a rate-dependent material, **the dynamics regularizes the problem for free**. CSCM, CDPM-Rel3, KCC, and Johnson-Cook all have rate dependence built in for this reason. Static analyses with the same materials need crack-band on top because the dynamic regularization isn't active.

**Caveat for explicit**: mass scaling can interact with this — if you scale mass aggressively, the effective wave speed in the scaled region drops, and the dynamic regularization length shrinks. Document this in [[#5.1]].

#### 10.3.5 Cosserat / micropolar continua

Add rotational DOFs at each point, with a coupled moment stress and a length scale built into the constitutive relations. Theoretically clean (Cosserat brothers 1909, Eringen 1960s), introduces a real microstructural length. **Rarely used in production** — DOF count nearly doubles, calibration is hard, most materials don't have measured Cosserat parameters.

Not on the Ladruno radar; mentioned for completeness.

#### 10.3.6 Phase-field fracture (Francfort-Marigo 1998; Miehe, Hofacker, Welschinger 2010)

Recasts cracks as a *smeared damage field* $d(\mathbf{x}) \in [0, 1]$ that evolves by minimizing a regularized Griffith energy:

$$\mathcal{E}[\mathbf{u}, d] = \int_\Omega (1-d)^2 \psi^+(\boldsymbol{\varepsilon}) \, dV + \int_\Omega G_c \left( \frac{d^2}{2\ell_0} + \frac{\ell_0}{2} |\nabla d|^2 \right) dV$$

where $\ell_0$ is an explicit length-scale parameter (the "regularization length"; smaller = sharper crack representation).

**Effect**: cracks emerge as smeared damage automatically; branching, merging, and arbitrary propagation paths happen with no special tracking. State-of-the-art for monotonic concrete fracture, especially in research.

**Cost**: the damage field is a global PDE (gradient-enhanced family), so each step is a coupled solve. Length scale $\ell_0$ has to be small (~$h_e$ to a few $h_e$) to represent cracks crisply, which means **fine meshes are required** — phase-field is expensive in absolute terms. But it gives the cleanest crack mechanics of any regularization scheme.

**Status in LS-DYNA**: not present. **Status in OpenSees**: not present. A potential future research direction; aligned with the user's XFEM interest.

### 10.4 What LS-DYNA actually does — survey

| Material / framework | Default regularization | Notes |
|---|---|---|
| `*MAT_159 / CSCM` | Crack-band (Theory + Vol II p. 1139) | Plus optional rate dependence |
| `*MAT_273 / CDPM` (Concrete Damage Plasticity Model) | Crack-band + viscous (rate) | Both active by default |
| `*MAT_072R3 / KCC (Concrete Damage Rel3)` | Crack-band | Vol II p. 582 |
| `*MAT_84 / Winfrith Concrete` | Crack-band, smeared crack approach | Theory p. 384 |
| `*MAT_15 / Johnson-Cook` | Rate dependence (inherent) | Damage softening regularized by rate |
| `*MAT_120 / Gurson` (porous plasticity) | None by default; `*MAT_NONLOCAL` wrapper available | |
| `*MAT_ADD_DAMAGE_GISSMO` | Crack-band via per-element fracture-energy input | The user supplies $G_f$ and characteristic length; GISSMO scales internally |
| `*MAT_NONLOCAL` | Generic nonlocal-integral wrapper | Wraps any history variable; expensive |

**Pattern**: LS-DYNA's default story is **crack-band + (free, when applicable) viscous**. Heavy options exist (`*MAT_NONLOCAL`) but are opt-in. No phase-field, no gradient-enhanced.

### 10.5 What this means for Ladruno — decisions

#### 10.5.1 Modular damage decorator ([[#5.4]]) — non-negotiable

The decorator **must accept characteristic length and fracture energy as inputs and apply crack-band scaling internally**, on day one. Without this, the decorator is mesh-dependent garbage and the "concrete fracture" capability claim is hollow.

Per-feature plan for [[#5.4]] should specify:
- `Gf` (fracture energy, $J/m^2$) as input.
- Either auto-derive characteristic length $h_e$ from the wrapped element (preferred), or accept it as input.
- Hard-warn if the user-supplied $h_e$ is so small that $G_f / (f_t \cdot h_e)$ implies unrealistic ductility.

#### 10.5.2 Framework gap — `Element::getCharacteristicLength()`

OpenSees today has no systematic way to ask an element "what's your size?" Some elements (`bbarBrick`, the new tetrahedra) compute it internally for their own purposes; most don't expose it.

**Proposal**: add a `virtual double Element::getCharacteristicLength() const` to the base class, default implementation returns $\sqrt[3]{V_e}$ for solids, $\sqrt{A_e}$ for shells, length for beams/trusses. This is small, low-risk, and unblocks crack-band-correct damage for every regularization-aware material/decorator going forward.

> [!question] Q1 (open)
> Should `getCharacteristicLength()` be on `Element` or on `Domain` (queryable from a node-set)? Element-level is more local and natural for damage models. Domain-level might be needed for nonlocal integral models that need the search radius to be set in physical units.

#### 10.5.3 Cohesive elements — naturally regularized

Future cohesive-element work ([[#5.15]]) doesn't need crack-band scaling — the cohesive law is parameterized directly by $G_f$ and traction strength, the band width is the cohesive surface itself, and the length scale is implicit. **Cohesive is therefore the path of least regularization-headache** for discrete cracks.

#### 10.5.4 Explicit + rate-dependent default story

For dynamic problems with softening, encourage rate-dependent material variants (Cowper-Symonds, Johnson-Cook with rate, CSCM, CDPM). Document in tutorials: "if your material softens and your problem is dynamic, prefer rate-dependent variants — you get viscous regularization for free." This is a **documentation / examples** deliverable, not a code one, but it should accompany [[#5.1]] and [[#5.3]].

> [!warning]
> Mass scaling ([[#5.1]]) interacts with viscous regularization. Effective wave speed drops in mass-scaled regions, the regularization length $\ell \sim c \cdot \tau$ shrinks. The diagnostics in [[#5.1]] should warn when mass scaling is high *and* the model uses softening rate-dependent materials.

#### 10.5.5 Nonlocal integral — keep on the radar, not v1

A nonlocal-decorator analogue to GISSMO is appealing (wrap any plasticity, add nonlocal averaging) but the parallel-scatter / MPI complications are significant. **Park as a future item**; revisit after [[#5.7]] (which builds the bucket-sort infrastructure that nonlocal also needs).

#### 10.5.6 Phase-field — research-grade future direction

Aligned with the user's XFEM interest and the long-term meshfree push ([[#5.8]] / [[#5.9]]). Distinct enough from EFG that it deserves its own item if pursued. Cost-benefit vs cohesive elements is the open question — phase-field is more general (automatic branching, no pre-defined surfaces) but much more expensive.

> [!question] Q2 (open)
> Phase-field fracture as a future numbered roadmap item, or merge into the existing 5.12 "known gaps" entry on XFEM? My instinct: separate item because the implementation patterns are completely different from EFG/XFEM (gradient PDE vs enrichment functions).

### 10.6 The regularization checklist for every per-feature plan

Going forward, every per-feature plan in this folder that involves softening must answer:

1. **Does the model soften?** (If no, skip.)
2. **What's the regularization strategy?** (Crack-band / rate / nonlocal / none.)
3. **Where does the length scale come from?** (`getCharacteristicLength()`, user input, material property.)
4. **What happens at the limits?** (Very fine mesh: does $G_f / (f_t \cdot h_e)$ go non-physical? Hard error or warn?)
5. **Diagnostics**: track the regularization quantities (band width, dissipated energy, $h_e$ used per element). The same energy-balance discipline mentioned in [[#8.10]] applies to fracture energy specifically.

Items currently affected: [[#5.4]] (modular damage), [[#5.7]] (cohesive contact, future), and any entry in [[#5.15]] that touches fracture (cohesive elements, element erosion). Items NOT affected: [[#5.1]], [[#5.2]], [[#5.5]], [[#5.6]], [[#5.10]], [[#5.11]] — these don't introduce softening on their own.

> [!summary] One-sentence framing
> Localization is what happens; regularization is what you do about it. The roadmap's fracture / damage / liquefaction promises are cashable only if every softening-capable item in this plan ships with a documented regularization strategy from day one.

---

## 11. Performance landscape and realistic ceiling

When you watch a 1-million-particle physics demo run at 60 fps in a game engine and compare it to a 100,000-element OpenSees model that takes hours, the gap looks indefensible. It isn't — but the explanation isn't "FEM people are bad at code." This section lays out the actual landscape so the roadmap's performance promises stay honest.

### 11.1 Why games and AI feel orders of magnitude faster

Four separate reasons compound:

#### 11.1.1 They're not solving the same problem

| | What's actually happening | Where the difficulty hides |
|---|---|---|
| **AI inference** | Dense matrix-multiplication chains on uniformly-shaped tensors | The hard part is *training* (weeks on GPU clusters); inference is just a forward pass |
| **Game physics** | Position-based dynamics, mass-spring, simple rigid contact, iterative Gauss-Seidel constraint projection | "Looks right" is the bar; energy can drift a few percent; penetrations of 1 mm are masked by graphics |
| **Scientific FEM** | Sparse irregular nonlinear constitutive integration with branchy state-machine material laws | Stress to 4 sig figs at every Gauss point matters; crack initiation, liquefaction triggering, collapse mechanisms depend on it |

Position-Based Dynamics (Müller et al. 2007) literally trades physical correctness for stability — it modifies positions directly to satisfy constraints without going through forces. You can't use that for earthquake engineering and pretend the result is real (though XPBD ([[#5.14]]) recovers some of the physics back, which is exactly why it's interesting).

#### 11.1.2 Accuracy bar separated by orders of magnitude

Games: "looks plausible at 60 fps" → 16 ms/frame budget; whatever physics gets done is shipped. Convergence isn't a concept — there's no test for it.

FEM: predict response well enough to inform design or research. PM4Sand benchmarked against centrifuge tests to a few percent. Energy balance checked to 1%. Stress at the crack tip to single-digit percent.

If you relax FEM's accuracy bar to game-physics levels, **you get game-physics speeds**. PBD-on-FEM-with-real-materials is a known research direction (NVIDIA FleX, Macklin et al.) — the speed is real, but the physics has to be very deliberately chosen.

#### 11.1.3 Memory and dispatch discipline

Games and AI grew up after the cache-cliff era and were written *for* it. Modern game engines use Entity-Component-System (ECS) architectures (Unity DOTS, Unreal Mass Entity, Bevy) — bit-packed components, SoA layout, hand-written SIMD inner loops, lock-free message passing. AI frameworks sit on cuBLAS/cuDNN, hand-tuned by NVIDIA to within a few percent of theoretical peak.

Scientific FEM grew up before this discipline became cheap. OpenSees inherits an OO-first late-90s C++ architecture (virtual dispatch on every element). LS-DYNA inherits Cray-1-vectorized Fortran from 1978 (faster than OpenSees on the inner loop, slower than a modern cuBLAS kernel by an order of magnitude).

[[#5.2]] is the move that closes the OO-dispatch gap *partially* — but FEM will never hit Unity's per-frame budget because we're solving harder problems.

#### 11.1.4 Funding asymmetry

The unstated reason. AI is the largest investment in human history. Games are a $200B+/year industry with full-time graphics and physics performance engineers per studio. NVIDIA tunes cuDNN with an engineer-decade per release. Scientific FEM runs on NSF/DOE grants; the OpenSees user community is several thousand people; the developer community is dozens at most.

A meaningful chunk of the speed gap is just the resource gap.

### 11.2 The actual solver landscape

For our purposes — large nonlinear FEM with explicit ambitions and the kind of problems Section 1.1 promises:

#### 11.2.1 Linear solvers

- **Sparse direct** (MUMPS, PARDISO, SuperLU, KLU): what OpenSees uses today via `system SparseGEN` / `Mumps`. Robust, no parameter tuning. Cost is super-linear in 3D — $O(N^{1.5})$ to $O(N^2)$ depending on bandwidth. Hits a wall above ~$10^6$ DOF on a workstation.
- **Iterative Krylov** (CG, GMRES, BiCGStab, MINRES): $O(N)$ per iteration. Convergence depends on conditioning; needs preconditioning to be useful.
- **Algebraic multigrid (AMG)** preconditioning (hypre, ML/MueLu, AMGX): the modern fast option for elliptic-like problems. Scales to $10^9$ DOF on clusters.
- **GPU sparse** (cuSPARSE, AMGX, Ginkgo): sparse iterative on GPU is now production-grade. Sparse *direct* on GPU remains an open problem — sparse triangular solve has branch divergence and irregular memory access that fights GPU SIMT execution.

OpenSees today: sparse direct only (effectively). Adding matrix-free Krylov ([[#5.13]]) is the path forward.

#### 11.2.2 Domain decomposition

OpenSees' parallel story is built on this:

- **METIS / SCOTCH / Zoltan**: graph partitioners that split a mesh into balanced subdomains while minimizing cut edges. METIS is the workhorse.
- **Additive Schwarz** (overlapping subdomains): simple, scales OK. Used as a preconditioner.
- **FETI** (Finite Element Tearing and Interconnecting, Farhat & Roux 1991): condense interior DOFs per subdomain, solve interface problem with Lagrange multipliers. Scales to $10^4$+ subdomains. The HPC standard for elasticity-class problems.
- **BDD / BDDC** (Mandel 1993, Dohrmann 2003): better scalability than FETI for nearly-incompressible problems. Newer; less mature in community codes.
- **OpenSees today**: domain decomposition via `Channel`-based `sendSelf/recvSelf`. Works for the scale of problems people actually run; not at FETI/BDDC sophistication.

Adopting FETI/BDDC inside OpenSees would be a separate research direction. Not currently on the roadmap; flagged here for completeness.

#### 11.2.3 GPU computing for FEM

This is where it gets interesting for the roadmap:

- **Explicit dynamics is GPU-natural**. No global linear solve. The `gather → kernel → scatter` pattern is what GPUs are built for. Game physics engines (PhysX, FleX, Bullet on GPU) have done this for years; LS-DYNA and Altair RADIOSS have GPU-offload modes. [[#5.12]] is OpenSees catching up.
- **Implicit FEM on GPU is hard but feasible**. You need GPU-friendly sparse linear algebra (now fine), GPU-friendly preconditioners (AMG via AMGX is mature), and matrix-free element kernels ([[#5.13]] gives us the contract). MFEM and deal.II demonstrate it.
- **Materials are the hard part**, not elements. Element kernels are mostly $\int B^T \sigma \, dV$ — straightforward parallel patterns. Material kernels (PM4Sand, ASDPlasticMaterial3D, CDPM) have branchy radial-return / multi-surface return-map code that fights GPU SIMT. Refactoring materials into pure-functional `(strain_n, state_n) → (stress, state_{n+1})` form (no hidden mutation) is the design discipline that unlocks GPU.

### 11.3 Matrix-free as the modern fast pattern

This deserves its own subsection because it's the pattern that connects [[#5.2]] (explicit batch dispatch), [[#5.13]] (matrix-free implicit), and [[#5.12]] (GPU offload) — all three want the same data layout.

The standard FEM pattern (what OpenSees does today):

```
1. Assemble K globally — sparse matrix, gigabytes for 10^6 DOF
2. Solve K · du = R via sparse direct factorization
3. Update u
```

The bottleneck is **memory**: storing $K$, factoring $K$. Not flops — memory.

The matrix-free pattern (what MFEM, deal.II do):

```
1. (no assembly)
2. Krylov iteration: at each step, compute K · v on the fly by 
   looping over elements and applying k_e · v_e — never form K
3. Update u
```

Cost moves from memory-capacity-bound to memory-bandwidth-bound. The element-local matvec is the same cache-friendly gather/kernel/scatter pattern as explicit. **CPU SIMD, GPU SIMT, and explicit batch dispatch all want the same shape of code.**

Architecturally, this means [[#5.2]]'s `batchAssembleResisting` and [[#5.13]]'s `batchApplyTangent` are siblings — same data layout, same gather/scatter, just a different inner kernel. Designing them coherently from the start is what makes the full performance story (explicit + matrix-free implicit + GPU + SIMD) compose cleanly.

### 11.4 The realistic ceiling for Ladruno

Stacked speedups, in execution order:

| Source | Speedup over OpenSees today | Cumulative | Status |
|---|---|---|---|
| **Batch dispatch** ([[#5.2]] / [[#5.3]]) | 2–3× on element loop | 2–3× | Roadmap v1 |
| **GPU offload** ([[#5.12]]) on supported elements | 2–4× over CPU batch on GPU-friendly elements | 4–12× | Roadmap medium-term |
| **Matrix-free implicit** ([[#5.13]]) for $>10^6$ DOF problems | 5–10× on memory-bound implicit (and unlocks problem sizes that were previously infeasible) | n/a — separate workload | Roadmap medium-term |
| **XPBD relaxation** ([[#5.14]]) on softening quasi-static | 5–10× over implicit Newton with arc-length, on the problems where Newton was failing anyway | n/a — separate workload | Roadmap research |
| **AMG / FETI / BDDC parallel solvers** (not on roadmap) | 2–10× on large parallel runs | additional | Future |

Best case stacking on a representative workload (large 3D explicit SSI, GPU available): **roughly 50–100× over OpenSees today**, on workloads where everything aligns. **Not 1000×.** The 1000× wins in AI come from specialized hardware (TPUs) on specialized workloads (dense matmul) that scientific FEM doesn't have.

This is honestly framed in [[#8.10]] already, but the realistic ceiling deserves its own statement: **make it 50–100×, not 1000×, and the platform argument still works**. Don't oversell.

### 11.5 Implications for grant proposals and PhD scoping

- **For funding agencies**: the realistic-ceiling story is what we should write in proposals. "Open-source explicit dynamics with matrix-free implicit, scaling to $10^7$ DOF, with cross-platform GPU support" — that's a fundable claim. "GPU-accelerated FEM rivaling AI inference speeds" — that's not, and reviewers will see through it.
- **For PhD students**: each of [[#5.12]] / [[#5.13]] / [[#5.14]] is dissertation-scale work. Don't treat them as side quests on top of an applications thesis. Either own the methods development, or use the methods someone else owns.
- **For users**: the documentation should explain the 50–100× ceiling, what enables it (problem class + element type + hardware), and what *doesn't* get there (general implicit nonlinear quasi-static with sparse direct on a single workstation — that's still bound by sparse solver performance).

> [!summary] One paragraph
> Game physics and AI are faster than scientific FEM because they solve different problems with different accuracy bars on different hardware with different funding. The architectural patterns from those worlds (SoA, batch dispatch, matrix-free, GPU-friendly kernels, position-based relaxation) are portable and worth porting. The realistic ceiling is 50–100× over OpenSees today on the workloads where it stacks, which is plenty for the platform argument. Don't promise more.

---

## 12. Beyond OpenSees — what a clean-slate rewrite would look like

This roadmap is incremental: extend OpenSees, don't replace it. That's the right call for the next few years — there's a community, a constitutive-model library, a workflow ecosystem, a pile of validated benchmark cases. But it's worth writing down what we'd build *if the constraints were different* — both as a horizon to steer toward and as a way to test whether any of [[#5]]'s decisions are short-sighted compared to where the field is going.

This section is **not a plan**. It's a thought experiment that sharpens the plan.

### 12.1 The framing question

If we threw out the OpenSees codebase tomorrow and built a scientific FEM framework from scratch in 2026, with full knowledge of what AI infrastructure, modern HPC, and game physics taught us — **what would we build?**

The answer is not "OpenSees but in Rust." That misses the point. The actual lessons are about **data layout, kernel design, differentiability, portability, and the relationship between the user-facing model and the runtime**.

### 12.2 Core architectural choices

#### 12.2.1 Storage: columnar, SoA from day one

Mesh data lives in a columnar store (Apache Arrow-style). Nodes table, elements table, materials table, history-variables table. Each is a contiguous array of one field, not an array of structs.

Why: every kernel becomes a stream over contiguous memory. The CPU prefetcher loves it. SIMD just works. GPU offload is a memcpy of a column, not a deep object graph. This is **the** lesson from AI tensor frameworks and game ECS architectures. Everything else follows from this choice.

Concrete: think Arrow + Parquet for on-disk persistence; in-memory layout that maps directly to GPU buffers.

#### 12.2.2 Element and material kernels as pure functions

```
element_kernel: (X, u, state_in)  → (F_internal, state_out)
material_kernel: (strain, state_in) → (stress, state_out)
```

No hidden mutation, no virtual dispatch in the inner loop, no framework callbacks. Each kernel is a pure function from input arrays to output arrays. Same kernel runs on CPU SIMD, on GPU, on a TPU-class accelerator, in the browser via WebGPU.

This is the LS-DYNA design discipline ([[#9]]) taken to its logical conclusion: elements and materials are *kernels*, not objects.

#### 12.2.3 Matrix-free, differentiable, and architecture-portable

- **Matrix-free by default**: never form $K$ unless asked. $K \cdot v$ is the primitive operation. Sparse direct is available as a backend, not the architectural assumption.
- **Differentiable from day one**: every kernel has a backward (adjoint) pass. Forward and reverse-mode automatic differentiation through the entire solver. Enables design optimization, parameter inversion, neural surrogate training, physics-informed neural networks (PINNs). This is what JAX-FEM and DiffTaichi prove out at small scale — bring it to engineering scale.
- **Architecture-portable kernels**: write the math once in a DSL or via a portability layer (Kokkos, RAJA, SYCL, MLIR-based). Compile to CPU SIMD, CUDA, ROCm, Metal, WebGPU, oneAPI. Don't lock to one vendor.

#### 12.2.4 DSL for constitutive models

The hardest part of OpenSees today is writing a new material. It's hundreds of lines of C++ with virtual methods, `commitState`/`revertToLastCommit` discipline, `sendSelf/recvSelf`, `getResponse`-string-tags. Most of that is framework boilerplate, not constitutive math.

A modern framework gives the researcher a DSL where they write the **math**:

```python
@material
def pm4sand(eps_rate, state):
    # ... clean Pythonic constitutive math ...
    return stress, new_state
```

The DSL compiler (think Taichi, JAX, MLX) emits efficient CPU + GPU code, generates auto-diff for free, handles commit/revert state semantics automatically. **Researchers write physics; the runtime handles performance.**

This is a productivity win of probably 10× for new constitutive model development. The current OpenSees model says "researchers should learn modern C++ to contribute"; the modern framework says "the runtime adapts to the researcher."

#### 12.2.5 Modeling frontend: declarative + scripted

The user-facing API stays familiar — Python for scripted models, plus a declarative "describe the analysis you want" layer for routine runs:

```yaml
analysis:
  type: pushover
  domain: my_building.hdf5
  load_pattern: lateral_increasing
  output:
    drift: roof
    base_shear: foundation
  parametric_sweep:
    fy: [40, 50, 60, 70]
```

The runtime decides how to parallelize, where to run (workstation, cluster, cloud), what to cache, what to checkpoint.

Existing OpenSees Tcl/Python scripts are first-class — preserved as a frontend layer over the new core for backward compatibility. This is what protects the existing user community during transition.

#### 12.2.6 Cloud-native, reproducible execution

Each analysis run is a function of `(model, parameters, ground motion, …) → (results)`. Trivially deployable to:
- A workstation
- A campus cluster (SLURM)
- AWS Lambda / GCP Cloud Run for short jobs
- Kubernetes for parametric studies

Every run produces a manifest with input hashes, RNG seeds, software version, hardware. Provenance is embedded, not bolted on. **Reproducibility is the default; you have to opt out of it.**

#### 12.2.7 AI integration first-class

- **Surrogate models** (NN-based response prediction trained on past simulations) are a configuration option: `material: surrogate(neural_pm4sand_v3)`.
- **Differentiable solver** enables PINN training, model calibration via gradient descent on observation residuals.
- **Active learning loops**: the framework knows which input regions have high model uncertainty and steers parametric studies toward them.

None of this is sci-fi — it exists in fragments today (JAX-FEM, NeuralOperators, DeepXDE). What's missing is a serious engineering-domain framework that integrates them.

### 12.3 What this would unlock

- **10–100× performance** on compute-bound workloads, *consistently*, not just on cherry-picked benchmarks.
- **Real-time interactive analysis** for moderate problems — adjust a load, see the stress redistribute live (WebGPU in a browser).
- **Differentiable structural design** — gradient-descend through the simulation to find optimal section sizes for a target performance metric.
- **Inverse problems at scale** — invert from observed seismic response to soil parameters; calibrate constitutive parameters automatically against test data.
- **Surrogate-augmented simulation** — replace the 5% of element types that are hot but not interesting with neural surrogates; keep the 95% as physical models. 100× speedup on the right problems.
- **Education on modern hardware** — students writing new constitutive models in a DSL get GPU performance on the first try, instead of a 6-month detour through C++ optimization.

### 12.4 What's preserved from OpenSees — the actual value

OpenSees' value is not the framework — it's:
- The **constitutive model library** (PM4Sand, PDMY, SAniSandMS, PressureDependentMultiYield, Concrete02, ConcreteCM, Steel02, ASDConcrete3DMaterial, ASDPlasticMaterial3D, ...). Decades of validated implementations.
- The **earthquake-engineering vocabulary** — DRM, fiber sections, ground motion patterns, recorders set up the way seismic researchers think.
- The **community** — thousands of users, validated benchmarks, published comparisons.
- The **academic legitimacy** — papers cite OpenSees results; new methods are validated against it.

A clean-slate rewrite **must port these forward**, not invent new ones. The constitutive models are the irreplaceable artifact. The framework is the disposable scaffold around them.

This is also what protects the rewrite from boil-the-ocean syndrome: the path forward is "**port models, replace framework**," not "rewrite everything."

### 12.5 Risks of starting from scratch

- **Boil-the-ocean syndrome.** Years of work to reach feature parity. By the time the rewrite is usable, the field has moved.
- **Loss of community during transition.** Users follow working code, not promising rewrites. If the new framework isn't drop-in usable, people stay on OpenSees.
- **Constitutive porting is months per model.** PM4Sand alone took years to validate the first time; re-porting is faster but not trivial. ASDPlasticMaterial3D's template-metaprogramming framework (Abell-Petracca) doesn't translate to a DSL trivially.
- **Validation re-do.** Every benchmark in every published comparison needs to be re-run. Significant labor; no immediate research output.
- **Funding mismatch.** "Rewrite from scratch" doesn't get NSF money. "Add capability X to OpenSees" does. The funding landscape favors the incremental path even when the rewrite would be technically better.

### 12.6 The pragmatic path: Ladruno as the bridge

The whole point of [[#1]] through [[#10]] is that **Ladruno is the incremental path**. But the architectural decisions made in Ladruno can either be informed by the rewrite vision or oblivious to it.

Concrete bridges:
- **[[#5.2]] batch dispatch** is the inner-loop pattern of the rewrite, retrofitted to OpenSees. The contract should be *stable enough* that a future rewrite can keep it.
- **[[#5.13]] matrix-free implicit** is the modern HPC FEM pattern. Its `batchApplyTangent` interface is what a clean-slate framework would have from day one.
- **[[#5.4]] modular damage decorator** is composable-constitutive thinking — the rewrite's DSL would support this natively.
- **Pure-functional materials** (refactoring constitutive models to have no hidden mutation, with `state_in / state_out` explicit) is the hardest discipline. Doing it incrementally inside OpenSees prepares the models for any future port.
- **Reproducibility instrumentation** — adopting input-hash + RNG-seed + version-manifest practices in Ladruno's per-feature plans is small effort now and infrastructure later.

These are the design choices that **make Ladruno's work portable to a future rewrite**, instead of being throwaway scaffolding. Worth being explicit about — every per-feature plan should ask "is this decision compatible with the clean-slate vision in Section 12, or are we baking in legacy?"

> [!important] The honest synthesis
> We're not going to rewrite OpenSees. The community, the constitutive models, and the funding model all argue for incremental work — Ladruno. **But every architectural decision in Ladruno should be made as if a rewrite is coming in 5–10 years.** That changes how we design interfaces, lay out data, structure constitutive models, and document our choices. The rewrite may never happen — but designing as if it might is what keeps Ladruno's work from being legacy on day one.

### 12.7 A name for the horizon

This vision needs a label so we can refer to it in future per-feature plans without retyping the whole thing. **Tentatively, *Ladruno Core*** — a name that signals "the rewrite-ready substrate," distinguishing it from *Ladruno* (the OpenSees fork) without committing us to actually building it.

When per-feature plans say "Ladruno Core compatible," they mean: this design choice could survive a future rewrite, because it follows the principles in this section.

---

## 13. Living document conventions

- This file is the index. Per-feature plans live as siblings in this folder, named `NN_short_slug.md`, copied from [[_template]].
- When a plan starts being executed, the per-feature file gains an `## Implementation log` section. When complete, the file moves to `Ladruno_internal/implemented_<name>.md` and this index updates the link.
- Open architectural questions live as `> [!question]` callouts here until they get resolved (then move to the per-feature plan or to a closed-decisions appendix).
- LS-DYNA citations: always use the manual + section + PDF page form (e.g. *Theory §26.11, PDF pp. 545–550*) so future-you can verify.
- OpenSees source paths: relative to `OpenSees/SRC/`.
- Cross-link freely with `[[wikilinks]]`. Renaming files is fine — Obsidian keeps links alive.

---

*Next concrete step: when you're ready to start, create `01_selective_mass_scaling.md` from [[_template]] and start the design sketch.*
