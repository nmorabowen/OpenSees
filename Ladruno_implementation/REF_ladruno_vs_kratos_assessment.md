---
title: Ladruno (OpenSees fork) vs Kratos Multiphysics — Comparative Assessment
status: reference
date: 2026-06-22
note: >-
  Ladruno facts are source-anchored to C:\Users\nmb\Documents\Github\OpenSees
  (branch ladruno). Kratos facts are from the kratos skill references built from
  C:\Users\nmb\Documents\Github\Kratos. Performance is assessed architecturally
  (parallel model, solvers, data structures), NOT from head-to-head benchmarks.
---

# Ladruno (OpenSees fork) vs Kratos Multiphysics

## 0. The fundamental asymmetry (read this first)

These are not the same *kind* of thing, and every comparison below is colored by it:

- **Kratos Multiphysics** is a large, multi-organization (CIMNE/TUM/Siemens/…) **general multiphysics framework**: a C++ core + ~48 plug-in *applications* spanning structural, fluid, geomechanics, DEM, particle methods, contact, IGA, optimization, FSI/co-simulation. HPC-first (threaded assembly + MPI), BSD-licensed, pip-installable, large community.
- **Ladruno** is a **single-maintainer, governance-heavy research fork of OpenSees**, additive-only (no core edits), focused on **nonlinear solid/RC structural mechanics and explicit/robust dynamics**, delivered as a curated Windows product with an ADR-driven design process.

So the honest framing is **breadth-and-scale (Kratos) vs depth-and-curation-in-a-niche (Ladruno)**. Each "wins" exactly where its design center is.

---

## 1. Architecture

### Kratos
- **Plugin/application model**: `KratosApplication::Register()` + `KRATOS_REGISTER_ELEMENT/CONDITION/CONSTITUTIVE_LAW`, `KratosComponents` name→prototype maps, a dotted-path `Registry` tree, and a `Kernel` that imports apps and merges variable lists. Each app is its own shared lib + pybind module (`kratos/includes/kernel.h`, `kratos_application.h`, `registry.h`).
- **Data model**: `Model` → `ModelPart` (sub-model-parts reference, not copy, entity subsets), `Variable<T>`-keyed data with an explicit **historical (time-buffered) vs non-historical** split, and **DOFs keyed by `Variable`** (`DISPLACEMENT_X`, `WATER_PRESSURE`, …). Adding a physics field = adding a Variable DOF — multiphysics is native to the data model (`model_part.h`, `node.h`, `variable.h`).
- **Solving layer**: generic `SolvingStrategy`/`Scheme`/`BuilderAndSolver`/`ConvergenceCriteria` templated over `<TSparseSpace, TDenseSpace, TLinearSolver>`, assembled declaratively from JSON. Python `AnalysisStage` orchestrates; `Parameters` (JSON) with `ValidateAndAssignDefaults`.

### Ladruno / OpenSees
- **Monolithic SRC tree**, classic OpenSees OOP: `TaggedObject → DomainComponent → Element/NDMaterial/...`, the **trial/commit/revert** state cycle, `classTags.h` + `FEM_ObjectBroker` for serialization, `OPS_*` elementAPI factories wired into both Tcl and Python interpreters. **Positional DOFs** (`ndf` per node), not variable-keyed.
- **Ladruno's architectural contribution is process + additivity, not a new core**:
  - Purely additive **leaf classes in a private per-registry classTag band ≥ 33000** (one above each upstream frontier → rebase-safe), every vanilla touch `// Ladruno`-tagged and ledgered. No application/plugin seam — the deliberate inverse of Kratos.
  - **ADR-driven governance**: 37 numbered ADRs with adversarial multi-lens review, P0 falsification tests *before* feature code, severity-rated risk registers, and three always-current ledgers (implementations / vanilla-files / quirks).
  - **Constraints as penalty/AL elements** (RBE2/RBE3, embedded node/rebar) rather than handler-managed `MP_Constraint`s, plus one new `LadrunoProjectionHandler` only where explicit dynamics demands it.
  - A **Python robust-solve orchestration spine** above the standard analysis objects (escalation ladder), with only two narrow C++ observability getters.

### Verdict — Architecture
| | Kratos | Ladruno/OpenSees |
|---|---|---|
| Module boundary | **True plugin apps** (clean to ship/compose) | Deliberately none (additive in-tree, upstreamable) |
| Extensibility seam | Registry + Variables (multiphysics-native) | classTag + trial/commit (structural-native) |
| DOF model | Variable-keyed (any field) | Positional (`ndf`) |
| Solver composition | Generic/templated, JSON-assembled | Component commands, interpreter-assembled |
| Governance | Institutional (CI, releases, many contributors) | **Unusually rigorous for a fork** (ADR + falsification + ledgers), single maintainer |

**Takeaway:** Kratos's architecture is built for *composing physics modules at scale*; Ladruno's is built for *disciplined, rebase-safe accretion onto a structural code*. Ladruno's per-feature governance is arguably tighter than typical Kratos app development; Kratos's framework-level modularity and registry are far more capable for multiphysics.

---

## 2. FEM models

### Kratos
Enormous breadth (file-anchored in the kratos skill):
- **Structural**: solids (small / total-Lagrangian / updated-Lagrangian, B-bar, mixed u/p), shells (thin DKQ/DKT, thick MITC, corotational), membranes, beams (corotational, Timoshenko), trusses, cables, nodal masses/springs.
- **Constitutive**: `ConstitutiveLawsApplication` — generic *templated* small/finite-strain plasticity, damage, plastic-damage, viscoplasticity, composites; the standard `ConstitutiveLaw` CalculateMaterialResponse* contract, one CL clone per integration point.
- **Geomechanics**: U-Pw coupled poromechanics (small + updated-Lagrangian, diff-order/Taylor-Hood, FIC-stabilized), Mohr-Coulomb/Coulomb return mapping with tension cutoff, retention laws (Van Genuchten), interface elements, UMAT/UDSM external models, staged construction.
- **Beyond structures**: fluids (CFD/FSI/biomedical/hydraulics), DEM + coupling, MPM/PFEM particle methods, contact, IGA, convection-diffusion, thermal, optimization, ROM.
- **Gap**: **no first-class fiber-section beam abstraction** (section response comes from resultant CLs + `Properties`, not fiber integration).

### Ladruno
Deep but narrow, with a coherent **"one core, many views"** doctrine (header-only OpenSees-free kernels, selector commands):
- **Solids**: `LadrunoBrick` (one class = std / B-bar / URI+hourglass / **EAS** Simo–Rifai, general-ν locking cure), Bézier (Bernstein) Tri6/Tet10 with **positive lumped mass for explicit**, a `solidTransformation` layer (`-geom linear|corot|finite`) where **corotational reuses the whole small-strain material library** and finite strain is multiplicative **Hencky log-strain**.
- **Plane**: `LadrunoQuad`/`LadrunoCST` selector elements (collapse upstream FourNodeQuad/SSPquad/EnhancedQuad/ConstantPressureVolumeQuad).
- **Beams**: `LadrunoDispBeamColumn` (disp-based, **per-IP tributary-length regularization** + embedded strong-discontinuity cohesive hinge, 2D + 3D biaxial), `LadrunoIMKBeam` (concentrated-plasticity IMK macro — no upstream equivalent).
- **RC & coupling**: host-agnostic **embedded rebar** (bond-slip CEB-FIP, Dhakal–Maekawa buckling), **RBE2/RBE3** general 3D couplings, embedded-node ties.
- **Materials**: `LadrunoConcrete3D` (Menétrey–Willam 3-invariant plasticity + dual scalar damage, CDPM2-grade, all five reduced views, non-symmetric tangent, IMPL-EX/Duvaut–Lions tiers), `LadrunoJ2` (Chaboche combined hardening — Abaqus `*PLASTIC,COMBINED`), `LadrunoJ2Finite` (co-rotating backstress), Lemaitre coupled ductile damage, biaxial cohesive hinge (B–K mode-mix), RC shell (MCFT + aggregate interlock + tension stiffening on the existing `ASDShellQ4`).
- Inherits all of base OpenSees (fiber-section **force-based** beams, the broad uniaxial/nD library, ASDConcrete3D/ASDSteel1D, ASDPlasticMaterial3D template framework).

### Verdict — FEM models
- **Multiphysics / soil / fluid / particle / contact / IGA → Kratos, decisively.** Ladruno is structural-only.
- **Nonlinear RC & seismic structural depth → Ladruno (on top of OpenSees).** Production-grade concrete (triaxial confinement + cracked-membrane MCFT), combined-hardening steel, IMK hinges, embedded rebar with bond-slip/buckling, RBE2/RBE3, and selector elements with explicit-ready lumped mass are a tighter, more curated kit than Kratos's structural app for that specific work.
- **Fiber-section beam-column** (force-based distributed plasticity) is an OpenSees strength Kratos structurally lacks — relevant for frame/seismic modeling and a clean cross-pollination candidate.
- Both have **generic plasticity frameworks** (Kratos `GenericSmallStrainIsotropicPlasticity`; OpenSees `ASDPlasticMaterial3D`; Ladruno's J2/concrete kernels) — convergent design, different ergonomics.

---

## 3. Performance & numerics

### Kratos
- **Genuine HPC framework**: shared-memory parallelism via `block_for_each` (**native threaded element/assembly loops**) *plus* MPI (Trilinos + Metis partitioning), scalable to thousands of cores. Solver backends include AMGCL multigrid, Eigen, Trilinos. Generic block/elimination builder-and-solvers; implicit and explicit strategies.
- This is the structural advantage: **threaded assembly + distributed-memory solvers are first-class**.

### Ladruno / OpenSees
- **Base OpenSees is process-parallel (MPI: SP partitioned-domain, MP independent-interpreter) with NO native OpenMP element-loop threading** (confirmed: the only `#pragma omp` in SRC is PFEM meshing). **Ladruno does not change this.**
- Ladruno's performance/numerics maturity is on a **different axis — single-rank robustness + explicit dynamics + instrumentation**, all *shipped* (source-verified, not just ADR'd):
  - Explicit: correct-first-step central difference with a built-in CFL/`dt_cr` guard, **HRZ mass-conserving lumping**, **selective mass scaling** (`CentralDifferenceSMS`, LS-DYNA DT2MS semantics) — the prerequisite for practical 3D explicit work.
  - Robust implicit: viscous-**stabilized arc-length** (Abaqus `*STATIC,STABILIZE`-like), **indirect/CMOD control** for snap-back, matrix-free **dynamic relaxation**, a true-unbalance convergence test, and a Python escalation-ladder driver with provenance verdicts.
  - A real **profiler** (state-determination ns/call, per-element-class timing, memory census — base OpenSees `Timer` is a Windows no-op) and an HDF5/SWMR recorder.
- **Weakest axis: parallel scaling.** Python MPI (`openseesmp`) is build-system packaging of the existing MPI surface; the one true scaling fix (a distributed `LadrunoParallelNumberer` replacing the rank-0-serialized numbering) is **designed but not built**.

### Verdict — Performance & numerics
- **Large 3D / multiphysics / many-core scaling → Kratos.** Threaded assembly + MPI + scalable solvers are architectural wins Ladruno cannot match today.
- **Hard-to-converge nonlinear quasi-statics (softening/snap-through/snap-back) and instrumented explicit dynamics → Ladruno.** Its robust-solve toolkit, SMS, and profiler are more turnkey for research-grade structural runs than what Kratos exposes out of the box.
- Net: **Kratos scales; Ladruno survives the hard nonlinear path and tells you why.** Different definitions of "performance."

---

## 4. Usability

### Kratos
- **Declarative Python + JSON**: `ProjectParameters.json` / `MaterialParameters.json` / `.mdpa`, orchestrated by an `AnalysisStage`. `pip install KratosMultiphysics-all`.
- **Ecosystem**: GiD / Salome / Gmsh front-ends, broad docs/wiki/tutorials, large international community, BSD license, cross-platform (Win/Linux/macOS).
- **Cost**: a conceptually heavy model (ModelPart/sub-model-parts, processes, solver wrappers, the variable database) to learn; GUI pre/post is partly via commercial GiD.

### Ladruno
- **Turnkey Windows product**: a real Inno Setup **GUI installer** with venv auto-wiring and `.pth` injection — non-developers go from `setup.exe` to `import opensees` without a build toolchain (vs OpenSees's source-build / `pip openseespy`).
- **`import openseesmp`** MPI Python module with bundled Intel MPI runtime; openseespy alias for drop-in compatibility (openseessp deliberately not shipped, with documented rationale).
- **Results UX**: a richer `.ladruno` HDF5 recorder (energy/envelopes/global results) **plus a live SWMR streaming monitor with a FastAPI/React web sidecar** — no stock-OpenSees analogue.
- **Opinionated, validated modeling**: RC-3D recipes with empirical validation gates, a "blessed" adaptive-stepping Python helper layer (`Ladruno_scripts/`), selector-style commands (`-formulation`), and an ADR/ledger/guide documentation corpus with a feature banner at import.
- **Toolchain**: `apeGmsh` as first-class mesher/driver/post-processor; STKO compatibility preserved via the frozen MPCO recorder.
- **Cost**: single-maintainer, **Windows-centric**, Python-3.12-ABI-locked, no pip wheels, niche audience; `.ladruno` not yet readable by STKO.

### Verdict — Usability
- **Newcomer wanting multiphysics + broad pre/post + community + cross-platform → Kratos.**
- **Researcher wanting a curated, installer-delivered, instrumented nonlinear-structural toolkit on Windows with live results and opinionated recipes → Ladruno.**
- For the **structural/earthquake niche specifically**, the OpenSees lineage (OpenSeesPy) has the larger global user base than Kratos's structural app; Ladruno layers a markedly better install + results experience on top of that lineage.

---

## 5. Cross-pollination opportunities

This is the actionable payoff (and extends the kratos skill's `opensees-bridge.md`).

### Kratos → Ladruno/OpenSees
1. **Threaded assembly.** Kratos's `block_for_each` threaded element loop is the model for the OpenSees performance gap (no OpenMP element threading). Hard — OpenSees elements return shared static matrices (not thread-safe) — but it's the highest-leverage perf idea.
2. **Variable-keyed DOFs / multiphysics data model.** If Ladruno ever needs coupled fields (thermal-mechanical, U-Pw), Kratos's `Variable`-DOF + historical/non-historical database is the cleaner design than positional DOFs.
3. **U-Pw poromechanics & retention laws.** Kratos GeoMechanics (`U_Pw_small_strain_element`, Mohr-Coulomb `CoulombImpl`, Van Genuchten) is well-factored reference if soil coupling enters Ladruno's scope.

### Ladruno → Kratos
1. **Robust-solve orchestration**: the stabilized arc-length + indirect control + dynamic relaxation + escalation-ladder pattern is more turnkey than Kratos's strategy stack for softening/snap-back.
2. **Fiber-section / IMK beam-columns**: Kratos lacks a first-class fiber-section beam; OpenSees's force-based fiber beam and Ladruno's IMK macro are strong donors.
3. **Selective mass scaling + HRZ lumping** for practical explicit dynamics.
4. **Governance discipline**: the ADR + P0-falsification + ledger process is transferable to any framework (including individual Kratos apps).
5. **Live SWMR monitoring** as a results-UX pattern.

### Convention bridge (already in the skill)
3D Voigt order matches: Kratos `[xx,yy,zz,xy,yz,xz]` (`constitutive_law.cpp:21`) ≡ OpenSees nDMaterial `[xx,yy,zz,xy,yz,zx]` (zx≡xz); engineering shear on both — so 3D continuum stress/strain ports 1:1. State-commit maps `setTrialStrain/getStress/getTangent → CalculateMaterialResponse*` and `commitState → FinalizeMaterialResponse*`.

---

## 6. Bottom-line matrix

| Dimension | Winner | Why (one line) |
|---|---|---|
| **Architecture — modularity/multiphysics** | Kratos | True plugin apps + registry + Variable DOFs |
| **Architecture — disciplined accretion** | Ladruno | ADR/falsification/ledger governance, rebase-safe additivity |
| **FEM — breadth (fluid/soil/DEM/contact/IGA)** | Kratos | 48 apps, genuine multiphysics |
| **FEM — nonlinear RC/seismic depth** | Ladruno | CDPM2 concrete, MCFT RC shell, IMK, embedded rebar, fiber beams (OpenSees) |
| **Performance — parallel scaling** | Kratos | Threaded assembly + MPI + scalable solvers |
| **Performance — nonlinear robustness & explicit** | Ladruno | Stabilized arc-length, indirect control, SMS, profiler |
| **Usability — ecosystem/community/cross-platform** | Kratos | GiD/Salome, pip, BSD, large community |
| **Usability — turnkey install + live results (Windows niche)** | Ladruno | GUI installer, SWMR monitor, opinionated recipes |

**One sentence:** *Kratos is the better general-purpose, parallel, multiphysics framework; Ladruno is the better-curated, more robust, better-instrumented tool for nonlinear RC/seismic structural research on Windows — and the two have a clean, specific set of ideas worth trading in both directions.*
