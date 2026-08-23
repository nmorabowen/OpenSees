---
title: "ADR 86 — Geomechanics-library cross-pollination, a poroelastic oracle corpus, and the MPM question"
status: DECIDED 2026-08-23, then REVISED after adversarial review — see §11. D1 = unsaturated IN (§3.3); D2 = A+B (§5.5, **affordability premise weakened — re-confirm**); D3 = all three skills (§7). **X1/F2 WITHDRAWN — the fixed-stress split is already shipped in ADR-73.** No code yet; §6 is the committed plan.
---

# ADR 86 — Three libraries, an oracle corpus, and whether MPM belongs in this fork

**Status:** scoping report. No code, no ledger rows. Reviews
[Anura3D](https://github.com/Anura3D/Anura3D_OpenSource),
[CODE_BRIGHT](https://deca.upc.edu/en/projects/code_bright) and
[OpenGeoSys 6](https://www.opengeosys.org/) against the fork's shipped
poromechanics lane, proposes a verification corpus, and assesses the Material
Point Method. Follows the precedent of
[[50_aem_opensees_scoping_report]] (assess-a-method → skills as deliverables)
and [[40a_kratos_crosspollination_amendment]] (borrow-from-a-sibling-code).
Companions: [[71_ladruno_up_family_adr]], [[73_ladruno_porous_overlay_adr]],
[[78_ladruno_up_corot_adr]], [[79_ladruno_geom_hypo_adr]].

---

## 0. Executive summary

| Question | Verdict |
|---|---|
| **Cross-pollination** | **OGS only** for code (BSD-3). CODE_BRIGHT = theory + BBM. Anura3D = formulation reference, **license-blocked** for code. Six harvests named in §3 — but **X1 (fixed-stress) is WITHDRAWN**: ADR-73 already ships it, measured. The surviving cheap item is the benchmark corpus. |
| **Oracle examples** | **Yes, but a narrower gap than first claimed.** ADR-71 §7.1 has no Mandel, Cryer, de Leeuw, Booker–Small or Gibson — and its closed-form gates are **all 1D columns**, so there is *no multi-dimensional closed-form coupled pressure gate*. (The original "all monotonic / a sign-flipped Q would survive" claim was refuted — §4.3, §11.) Corpus in §4.3; Tier-1 is pure numpy, no build. |
| **MPM in OpenSees** | The tree ships particle + background-grid infrastructure (`SRC/element/PFEMElement/`, Zhu & Scott) — but source reading (§11) cut this from "far cheaper" to **"integration risk removed; the reusable part is the easy ~15%"** (§5.1). **DECIDED: A + B**, on a premise that has since weakened — §5.5 flags it for re-confirmation. Own ADR; phases in §6-M. |
| **Way forward** | §6 — three tracks: **F** (harvests + oracles), **U** (unsaturated program), **M** (MPM program). F1 is unblocked today and costs no build. |
| **Unlocks (why)** | §12 — U removes the water-table ceiling (rainfall landslides, tailings, wetting collapse); M removes the post-failure ceiling (runout, breach extent, debris); together, one continuous analysis from ground motion to deposition that no surveyed code does. T bridges to design practice via macroelement calibration. |
| **Skills** | **All three approved** — `opensees-poromechanics`, `soil-constitutive-models`, `material-point-method`. §7. |

The strategic caveat, stated once and then assumed: **none of the three
libraries targets what this fork targets.** OGS and CODE_BRIGHT are nuclear-waste
/ geothermal / CO₂ codes — slow, unsaturated, thermally driven, geological
timescales. Anura3D is landslide runout. The fork's roadmap names **seismic
liquefaction with discontinuities**, and its PM4Sand / PM4Silt / ManzariDafalias
/ PDMY stack plus DRM, absorbing boundaries, contact and element removal has no
counterpart in any of them. The harvest below is therefore *targeted mechanisms*,
never a change of direction.

---

## 1. Context — what the fork already holds

Shipped, from the seam review of this tree:

- **`LadrunoUP`** (33017) — Biot u-p saturated-porous continuum, one class over
  T3/Q4/H8 (equal-order + McGann stabilization) and Bézier T6/Tet10
  (Taylor–Hood), `-formulation std|bbar`, `-geom linear|corot`. The load-bearing
  divergence from upstream is the **honest pressure DOF** (p, not ∫p·dt), which
  buys well-posed statics at the price of an unsymmetric tangent.
  ([[71_ladruno_up_family_adr]], [[LadrunoUP_guide]])
- **`LadrunoPorousOverlay`** — a pore-fluid field with its own life-cycle, so
  element removal no longer seals the flow path. **Its coupling is the
  fixed-stress staggered split, solid-first** (ADR-73 § "Decision" item 3), with
  the Kim–Tchelepi–Juanes modulus, and measured ≡ monolithic backward Euler
  (crack probe 0.00 %, Terzaghi 2.8e-8, factorization ×3.3 cheaper). ADR-71 §12
  records **P7 as OPENED by ADR-73**, not as unstarted research.
  ([[73_ladruno_porous_overlay_adr]])
- **Soil constitutive stack** — `SRC/material/nD/UWmaterials/` (PM4Sand,
  PM4Silt, ManzariDafalias ±RO, BoundingCamClay, J2CyclicBoundingSurface,
  DruckerPrager) and `SRC/material/nD/soil/` (PDMY{,02,03}, PIMY,
  FluidSolidPorousMaterial), plus CycLiqCP/SP, stressDensityModel,
  ASDPlasticMaterial3D.
- **Open phases** — ADR-71 **P5** (Q9/H20 providers, demand-driven), **P6**
  (hybrid u/p incompressibility, reserved). **P7 is NOT open** — it was opened
  and carried by [[73_ladruno_porous_overlay_adr]] (ADR-71 §12, 2026-07-13).
  ADR-71 §6's three-route list predates that and was superseded inside its own
  document; reading §6 alone gives a false picture of the fork (§11, finding 1).

Not held, verified by grep over `SRC/`:

- **No unsaturated capability of any kind.** No suction, no retention curve, no
  van Genuchten, no relative permeability. The only hits for those terms are
  RC-concrete tension retention and the Qz pile materials. ADR-71 §3.6 fences
  the element to saturated Biot deliberately; ADR-71 §3.6 Zone III (full Biot
  u-p-U) is likewise out of scope by decision.

That absence is the single largest white space these libraries speak to, and it
is the subject of decision **D1** (§3.3).

---

## 2. The license map — the hard constraint

**OpenSees is not BSD.** The grant lives in **`COPYRIGHT`** at the repo root —
there is no `LICENSE` file — and it is a UC Regents grant:
use/modify/distribute for *"educational, research, and non-profit entities for
noncommercial purposes only"*, others *"for internal purposes only"*, with
commercial distribution requiring a separate agreement with UC's Office of
Technology Licensing.

| Code | License | Can we take source? |
|---|---|---|
| **OpenGeoSys 6** | BSD-3-Clause | ✅ yes, with attribution |
| **CB-Geo `mpm`** (§5) | MIT | ✅ yes, with attribution |
| **`cb-geo/mpm-benchmarks`** | CC-BY-4.0 | ✅ inputs/data, with attribution |
| **Anura3D** | LGPL-3.0-or-later | ❌ **no** |
| **Karamelo** (§5) | GPL | ❌ no |
| **CODE_BRIGHT** | none published; **binaries only** | ❌ n/a — no source exists to take |

The Anura3D verdict is the one that bites. A copyleft licence forbids adding
restrictions to the licensed material; the OpenSees grant adds exactly such a
restriction (noncommercial / internal-use-only). Copying, porting, or
hand-translating Anura3D source into this fork and shipping the result puts the
two in conflict — and a careful Fortran→C++ transcription of an LGPL file is
still a derivative work of that file. The only clean arrangement would be a
separately distributed, dynamically linked component with relinking preserved,
which for Fortran MPM kernels inside an OpenSees element is not a real option.

> [!warning] **Anura3D is read-only.** Cite it, learn from it, reproduce its
> *published equations* from the papers. Do not open its `.FOR` files with an
> editor while writing fork code. If any of §3.3 proceeds, this gets a
> `LEDGER_quirks.md` row.

This is a reading of licence text by an engineer, not legal advice. The
mechanism above is concrete enough to plan around; if anything in §3.3 becomes
real, it deserves a proper check.

---

## 3. Track A — cross-pollination

### 3.1 What each library actually is

**OpenGeoSys 6** — C++, CMake, BSD-3, ~30k commits, 32 releases, PETSc/MPI
domain decomposition, VTK output, a `ogs` wheel on PyPI. Processes:
`LiquidFlow`, `RichardsFlow`, `TwoPhaseFlow`, `HeatConduction`,
`HeatTransportBHE`, `ComponentTransport`, `SmallDeformation`,
`LargeDeformation`, `PhaseField`, `HydroMechanics`, `RichardsMechanics`,
`ThermoHydroMechanics`, `ThermoRichardsMechanics`, `TH2M`,
`SteadyStateDiffusion`. Monolithic for all coupled processes; **staggered
available for TH, HM and phase-field**.

**CODE_BRIGHT** — FEM THMC for geological media, UPC/CIMNE, funded by the
nuclear-waste consortium (SKB, Posiva, GRS, ANDRA). Distributed as **compiled
64-bit binaries installed as a GiD problemtype**; Intel Fortran runtime
required. No source. User's Guide dated April 2026, so actively maintained.
Home of the **BBM (Barcelona Basic Model)**, including BBM interface elements.

**Anura3D** — MPM (not FEM), Fortran, LGPL-3, ~174 stars, ~173 commits, **no
tagged releases**, GiD + ParaView. Its multiphase formulations are real and
production-grade: `MPMDYN2PhaseSP.FOR`, `MPMDYN3PhaseSP.FOR`,
`TwoLayerFormulation.FOR`, `MPMDYNConsolidation.FOR`,
`MPMDYNUnsatConsolidation.FOR`, `MPMDynamicExplicit.FOR`. Its **constitutive
library is nearly empty** — `src/Soilmodels/` is `A3DLinearElasticity.f`,
`A3DMohrCoulombStandard.f`, `A3DBingham.f` plus UMAT twins and `aba_param.inc`;
real models arrive through `ExternalSoilModel.for` as user DLLs on an Abaqus
UMAT ABI.

### 3.2 The six harvests, ranked

| # | Harvest | From | Cost | Target |
|---|---|---|---|---|
| ~~**X1**~~ | ~~Fixed-stress split for staggered HM~~ — **WITHDRAWN, already shipped** | — | — | see below |
| **X2** | **Poroelastic benchmark corpus** | OGS + classics | low | §4, own PR |
| **X3** | **Unsaturated u-p** (retention + rel-perm + BBM) | OGS `RichardsMechanics` + CODE_BRIGHT theory | **high** | new ADR — **gated on D1** |
| **X4** | **Media-property abstraction** (MPL pattern) | OGS `MaterialLib/MPL` | medium | `LadrunoUP` parser evolution; prerequisite of X3 |
| **X5** | **Constitutive-plugin ABI** (MFront/MGIS) | OGS + Anura3D UMAT hook | medium | own architecture ADR |
| **X6** | **Explicit 2/3-phase + double-point** formulations | Anura3D (papers only) | reading | ADR-71 **P7** option set |

**X1 — fixed-stress split.** OGS uses the fixed-stress split for staggered HM
"since it turned out advantageous", iterating H → M → H → M. ADR-71 §6 names
three routes for explicit/staggered u-p and does not name fixed-stress — from
which this report originally inferred a gap. **That inference was wrong, and the
harvest is withdrawn.** ADR-73 *is* the fixed-stress ADR: "Coupling =
fixed-stress staggered split, solid-first", carrying the Kim–Tchelepi–Juanes
modulus, measured ≡ monolithic backward Euler (crack probe 0.00 %, Terzaghi
2.8e-8, both sub-solves SPD, factorization ×3.3 cheaper), and already citing
Kim–Tchelepi–Juanes 2011, Mikelić–Wheeler 2013 (the contraction proof) and
Turska–Schrefler 1993 (the Δt/h² bound). ADR-71 §6's list was superseded by
ADR-71 §12 in the same document. Reading §6 in isolation produced a
confidently-stated hole in a lane the fork had already measured — the exact
failure mode [[LEDGER_quirks]] exists to prevent.

**What genuinely remains open** (much smaller, and not a harvest from OGS):
whether a fixed-stress mode belongs *inside* the monolithic `LadrunoUP` element,
as distinct from the overlay that implements it today. That is a fork question,
not a cross-pollination item, and it is demand-driven.

**X4 — the property abstraction.** `LadrunoUP` today takes `-Kf -poro -perm
-alpha -Ks` as **scalars**. That is correct and sufficient for saturated Biot.
The moment saturation enters, permeability and storage become *functions of
state* and scalars stop working. OGS's MPL (medium / phase / component ×
property, with declared variable dependence) is the design that survives that
transition. Worth adopting **before** X3, not during it.

**X5 — the plugin ABI.** Two independent codes converged on the same answer:
constitutive models as external shared libraries behind a generic interface
(OGS→MFront/MGIS; Anura3D→Abaqus UMAT DLL). This is the same disease ADR-71 §1.2
diagnosed (N geometries × M formulations = hand-copied classes that drift), one
level up. It must be weighed against what `ASDPlasticMaterial3D` already
provides in-tree and against the fork's build/packaging story (a DLL-loading
material would interact with the installer and the lazy `.pth` alias fixed in
#735). **Own ADR; do not decide
here.**

**X6 — Anura3D as reading.** Its two- and three-phase MPM runs **fully
explicit**, in production, which is precisely the regime ADR-71 §6 filed as
"research". And `TwoLayerFormulation.FOR` (double-point: separate solid and
liquid material points) is convergent evolution with ADR-73's overlay — both
answer "the fluid must outlive the skeleton's continuum". Read the papers, cite
the equations, take nothing from the tree (§2).

### 3.3 Decision D1 — is unsaturated in scope?

X3 is the largest item and the only one that changes the fork's stated
boundaries. ADR-71 §3.6 fenced `LadrunoUP` to saturated Biot **on purpose**, and
the roadmap's target (seismic liquefaction) lives below the water table where
that fence is correct. Unsaturated buys: partially saturated slopes, rainfall
infiltration, capillary barriers, embankment dams above phreatic surface,
compacted-fill collapse, and the BBM class of models.

**DECIDED 2026-08-23 — D1 = YES.** The unsaturated lane is open. Consequences,
binding on everything downstream:

- **X4 before X3 — the chosen order, on grounds of avoiding a double refactor**
  (softened after review, §11-6). `LadrunoUP`'s current surface (`-Kf -poro
  -alpha -Ks` scalars, `-perm` an ndm-vector) cannot *express*
  saturation-dependent properties, but nothing prevents prototyping unsaturated
  behind flags with kernel closures S_r(p), k_rel(S_r) first — that is how most
  codes grew Richards support. The argument for X4-first is that doing so means
  refactoring twice; it is prudence, not necessity, and it should not be
  defended as though it were necessity.
- **ADR-71 §3.6 is amended, not discarded.** The saturated envelope
  (Zienkiewicz–Chang–Bettess zones I/II, u-p valid, zone III out) still governs
  the *dynamic* validity of the formulation. Unsaturated adds a second axis
  (degree of saturation) to a fence that was drawn on a rate axis; the two are
  independent and both stay.
- **The effective-stress seam is where the physics enters.** Today the material
  sees σ′ and p never enters the constitutive update. Unsaturated breaks that
  cleanly-drawn line — Bishop-type effective stress makes χ(S_r) a constitutive
  quantity, and BBM takes *suction as an independent state variable* rather than
  folding it into σ′. **Which of the two the fork adopts is the first real
  decision of track U** (U2), and it determines whether the existing
  `NDMaterial` contract survives untouched.
- **CODE_BRIGHT is promoted** from "reference only" to the primary theory source
  for this track (§3.1), alongside OGS `RichardsMechanics` as the worked C++
  precedent.

### 3.4 What we deliberately do NOT take

- **TH2M / thermal / chemical processes.** Out of scope. The fork has no thermal
  field ambition and importing one would be a strategy change, not a harvest.
- **Anura3D's constitutive models.** Three models, all weaker than what the fork
  ships. Going there would be a downgrade.
- **CODE_BRIGHT's validation emphasis.** Bentonite barriers over geological time
  is a different problem class; its benchmark set does not transfer to seismic.
- **OGS's process/solver architecture wholesale.** It is a well-built code, but
  the fork is an OpenSees fork and its Element/Material/Domain contracts are
  fixed. Borrow patterns, not layering.

---

## 4. Track B — the oracle corpus

### 4.1 The fork already has the idiom

`tests/` holds 335 files and a mature, honest verification pattern, exemplified
by `tests/ladruno_up_reference.py`: an **independent numpy re-derivation** of
every kernel block from first principles — own shape functions, own Jacobian,
own quadrature, own assembly — written under an explicit independence protocol
("reads NOTHING under `SRC/`"), gated at ≤1e-9 against the C++ kernel. Paired
with `tests/ladruno_up_kernel_check.cpp` and `tests/ladruno_up_cases.txt`.

So the question is not *can we build oracles* — the machinery exists and is
good. The question is **which physics is currently unguarded**.

### 4.2 Three tiers, and the trap

- **Tier 1 — closed form / series.** An analytic solution to the *same* boundary
  value problem. This is verification proper. Cheap: pure numpy, no build, no
  OpenSees.
- **Tier 2 — cross-code comparison.** Another code's answer to the same input.
  This is **not verification** and must never be gated as if it were.
- **Tier 3 — physical.** Centrifuge and shake-table data (VELACS, LEAP;
  Dewoolkar 1996 is already cited in ADR-71 §3.6). Validation, not verification.

> [!warning] **The cross-code trap.** OGS's `HydroMechanics` is not our u-p. Its
> pressure sign convention, its storage definition, its Biot-α treatment, its
> element technology and its time integrator all differ. A number that disagrees
> tells you nothing about which code is wrong until every one of those has been
> reconciled — which usually costs more than deriving the analytic solution. Use
> Tier 2 as a **plausibility screen and a bug-smell detector**, never as a gate.
> Where OGS ships an analytic reference *alongside* its benchmark, that reference
> is Tier 1 and the OGS run is decoration.

### 4.3 The proposed corpus

ADR-71 §7.1 pins **B1** ZCB80 periodic-load layer, **B2** ZS84 Example-1 column,
**B3** Boone–Ingraffea poroelastic column, **B4** McGann checkerboard footing,
**B5** Simon–Zienkiewicz–Paul dynamic column, plus Terzaghi 1D at P1. Grepped
and confirmed: **Mandel, Cryer, de Leeuw, Booker–Small and Gibson appear
nowhere in the ADR.**

**The first version of this section overstated the gap and has been corrected**
(§11, finding 3). The claim that every pinned benchmark is monotonic is false —
B1 is a periodic-load layer gated on phasor amplitude and phase, and B5 pins
step/sine/spike legs with p-oscillation explicitly gated. The claim that a
sign-flipped or transposed Q could survive the battery is also false, twice
over: ADR-71 **P0** already cross-checks the Q block against an independent
numpy oracle at ≤1e-9, and **B3** (Boone–Ingraffea) gates the undrained
Skempton response p(0⁺) = 0.410 MPa — a hard number generated purely by coupling
and storage, and the only α≠1 gate. ADR-73's own P0 additionally exercised the
Mandel–Cryer dip in the numpy-toy idiom (PR #575).

**The gap that actually survives is narrower and still worth closing:** every
closed-form-gated benchmark in the fork is a **1D column** (B4 is 2D but gated
only on a checkerboard index). There is **no multi-dimensional closed-form
coupled pressure-history gate anywhere in the fork.** Mandel and Cryer are the
cheapest such gates in existence, and the bug class they guard is
dimension-coupling and 3D-provider physics — O2's stated target — not a sign
flip. O1–O3 remain a worthwhile single PR on that basis, and on no other.

| ID | Benchmark | Physics probed | Reference | Bug class it catches | Phase |
|---|---|---|---|---|---|
| **O1** | **Mandel** (2D strip, Abousleiman generalization) | Mandel–Cryer non-monotonic pressure | closed form | **coupling-block sign/transpose; storage scaling** | Tier 1, now |
| **O2** | **Cryer** (sphere under hydrostatic step) | same effect, 3D, no boundary-layer excuse | closed form | 3D block placement; H8/Tet10 providers | Tier 1, now |
| **O3** | **de Leeuw** (2D consolidation) | Mandel family cross-check | series | independence check on O1 | Tier 1, now |
| **O4** | **Terzaghi 1D** *(already at P1)* | drainage/rate race | series | baseline | shipped path |
| **O5** | **Booker & Small** strip-footing consolidation | 2D consolidation with a real footing | series | mixed drainage BCs | Tier 1 |
| **O6** | **Gibson–England–Hussey** finite-strain consolidation | large-strain consolidation | series | `-geom corot` u-p ([[78_ladruno_up_corot_adr]]) | Tier 1, gated on need |
| **O7** | OGS `hydro-mechanics/consolidationbenchmark` **staggered** | staggered HM convergence | OGS (Tier 2) + analytic where given | X1 fixed-stress work | with X1 |
| **O8** | OGS `richards-mechanics/*` | unsaturated poromechanics | OGS (**Tier 2 — screen only**) | active under D1 = yes | U3 |
| **O8a** | **Srivastava–Yeh (1991)** exponential-K transient infiltration; **Tracy** exact 2D/3D Richards solutions; **Philip** series | unsaturated flow, **closed form** | analytic (**Tier 1**) | the §4.2 doctrine applied to track U — O8 alone would gate unsaturated work on a cross-code number | U2/U3 |
| **O9** | CB-Geo `mpm-benchmarks` hydrostatic column / sliding block | MPM baseline | CC-BY-4.0 data | **only if §5 proceeds** | gated |

**O1–O3 are the recommendation.** They are pure numpy in the established
`*_reference.py` idiom, need no build, no C++ and no new element, and they close
a real hole in the strongest lane the fork owns. Estimated one focused PR.

---

## 5. Track C — MPM in OpenSees

### 5.1 The finding: the tree already ships particles and a background grid

`SRC/element/PFEMElement/` (Minjie Zhu; Zhu & Scott, *Computers & Structures*
2014 and *J. Struct. Eng.* 2022) contains, beyond the PFEM elements themselves:

- **`Particle.h/.cpp`** — a particle carrying `coord`, `velocity`, `accel`,
  `pressure`, `pdot`, a group tag and a sub-step `dt`, with `moveTo`/`move`
  convection.
- **`BackgroundMesh.h/.cpp`, `BackgroundGrid`, `BCell`, `BNode`** — a **fixed
  structured background grid** with `setRange`/`setBasicSize`, cell indexing
  (`getIndex`/`lowerIndex`/`upperIndex`/`nearIndex`/`getCorners`), particle
  bookkeeping (`addParticles`, `gatherParticles`, `moveParticles`,
  `convectParticle`, `moveFixedParticles`), grid construction (`gridNodes`,
  `gridFluid`, `gridFSI`, `gridEles`), shape-function evaluation on the grid
  (`getNForTri`, `getNForTet`), particle↔grid transfer (`interpolate`), teardown
  (`clearGrid`, `clearBackground`, `clearAll`) and `remesh()`.
- It **creates real `Node` objects and adds them to the `Domain` every step**
  (`new Node(...)` → `domain->addNode(...)`), then tears them down.
- It is **exposed to both interpreters** — `OPS_getBgMesh()` is reached from
  `OpenSeesMiscCommands.cpp` and `OpenSeesOutputCommands.cpp`, with its own
  recorder path.

**Corrected after review (§11, finding 2).** The API inventory above is accurate
— every named method exists — but the *reuse* claim it was used to support was
overstated, and two items must be struck or qualified:

- **`BackgroundGrid` is effectively dead code.** Nothing outside itself,
  `BackgroundStructure.h` and `CMakeLists.txt` references it, and its
  `GridIndex` is 2D-only. The live grid is a `std::map<VInt, BCell/BNode>` inside
  `BackgroundMesh`. It was listed here from a file listing, not from source.
- **The transfers that exist are the wrong operator pair for MPM.** P2G is an
  SPH-style quintic-kernel weighted average of velocity/pressure onto new grid
  nodes — *non-conservative, not mass-weighted*; G2P (`interpolate`) is pure PIC
  with a hidden dt-halving convection clamp, on a different basis. MPM needs
  same-basis, mass/momentum-conserving transfers in both directions.
- **Grid nodes are not ndf-flexible** — `gridNodes` hard-codes ndm-velocity nodes
  plus paired ndf=1 pressure nodes bound by `Pressure_Constraint`.
- **There is no parallel story at all.** `BackgroundMesh` is a static singleton,
  has no `sendSelf`/`recvSelf`, and carries an in-source admission:
  `// TODO: setCenterNode will add nodes to domain. Can't be in parallel.`
- **The per-step lifecycle has its own analysis machinery** —
  `SRC/analysis/analysis/PFEMAnalysis.{h,cpp}` and
  `SRC/analysis/integrator/PFEMIntegrator.h`. Per-step `domainChanged` / SOE
  resize is a real cost and a real deliverable, not a free ride.

**Honest restatement.** The genuinely generic, reusable code is roughly the
container and index bookkeeping — `Particle`, `BCell`, `BNode`, the map-based
grid — against a fluid-specific `BackgroundMesh.cpp` of ~4,100 lines. What §5.1
actually buys is **integration precedent**: proof that the `Domain` tolerates
per-step node churn, and a worked interpreter/recorder pattern for a
particle method inside OpenSees. That is a **risk reduction, not a cost-category
change** — perhaps the easy 15 %. The §5.2 delta below is the hard core of MPM
and remains untouched by this finding.

### 5.2 The delta — what MPM needs that PFEM-BackgroundMesh lacks

Not small, but bounded and *listable*:

1. **Particles must carry constitutive state.** PFEM particles carry (x, v, a,
   p, ṗ). MPM particles must carry **σ, F, volume, mass, and the full material
   history** — i.e. an `NDMaterial` instance per particle, with commit/revert
   semantics. This is the largest single item and it touches the fork's
   state-commit cycle.
2. **The grid solve changes character.** PFEM-BG solves incompressible
   Navier–Stokes by fractional step. MPM solves **momentum with internal force
   assembled from particle stress** (∫Bᵀσ over the particle's domain, mapped to
   grid nodes). Different operator, same grid.
3. **Shape functions must be upgraded.** Linear grid functions produce the
   **cell-crossing error**: a particle crossing a cell boundary meets a
   discontinuous shape-function gradient and emits spurious stress oscillation.
   Production MPM answers this with GIMP, CPDI/CPDI2 or B-spline bases. Not
   optional for a soil code.
4. **Stress-update scheme.** USF / USL / MUSL variants, and the particle
   volume/domain update rule that goes with the chosen basis.
5. **Δt machinery.** Reuses the fork's existing explicit stack conceptually but
   the critical time step is grid-based, not element-based.
6. **Conservative, same-basis transfers.** Replacing the SPH-average P2G and PIC
   G2P (§5.1) with mass/momentum-conserving operators, and choosing a PIC/FLIP
   or APIC blend.
7. **Grid-node disentanglement.** Freeing grid nodes from the ndm-velocity +
   `Pressure_Constraint` pairing so they can carry MPM's DOF layout.
8. **An analysis driver.** MPM's own analysis/integrator classes alongside
   `PFEMAnalysis`/`PFEMIntegrator`, and an answer to per-step `domainChanged`
   and SOE-resize cost.
9. **Parallel — in or out.** The substrate is provably serial (§5.1). The fork
   ships SP/MP, and §5.5's standing justification is framed against an MPI code.
   This must be decided at M0, not discovered at M5.

### 5.3 The honest pathology list

MPM is not a strictly better FEM. Known and documented:

- **Cell crossing** (§5.2-3) — the signature MPM artifact.
- **Quadrature error** — particles are quadrature points that drift into bad
  configurations; integration accuracy degrades as they cluster or void.
- **Volumetric locking** — near-incompressibility (i.e. *undrained soil*, i.e.
  the fork's actual target) locks, needing F-bar / TLMPM-class remedies. The
  fork already fought this battle in FEM ([[81_quadratic_hex_limit_load_measurement]]).
- **Null-space / ringing instability** — grid modes the particles cannot see.
- **Energy conservation** — PIC damps, FLIP rings; every code picks a blend.

And one strategic collision worth naming plainly:

> [!important] **MPM's contact is automatic, and that undoes an investment.**
> In standard MPM, bodies sharing grid nodes interact through the grid with
> **automatic no-slip contact** — free, robust, and *not what this fork has spent
> ADRs 39/41/47/48/55/56/57/85 building*. The fork's differentiator in the
> collapse space is an explicit, frictional, augmented-Lagrange contact engine
> with real interface physics. Adopting MPM for the discontinuity problem does
> not extend that work; it **routes around it**. The same tension exists with the
> element-removal / FDEM-lite line ([[50_aem_opensees_scoping_report]],
> [[51_ladruno_element_removal_adr]], [[54_ladruno_fdem_lite_adr]]), which
> targets the same use case by a different road.

### 5.4 Three architectures

**A — a full MPM lane in the fork.** Build particles-with-material, the
momentum-form grid solve, GIMP/CPDI, the explicit driver, two-phase MPM. Cheaper
than it looks (§5.1) but still a multi-ADR program. It would compete directly
with **CB-Geo `mpm`** — MIT-licensed, C++, MPI+OpenMP, explicit/implicit/
semi-implicit solvers, linear/quadratic/GIMP/CPDI/CPDI2 bases, tested to ~15,000
cores, with its own benchmark repo — a code we may legally read *and reuse*, and
which will remain better at runout than a fork lane ever will.

**B — sequential FEM→MPM handoff.** The published hybrid (arXiv 2412.08040,
liquefaction-induced tailings-dam failure): **FEM captures the initiation
mechanism, MPM simulates the runout**, with a single transfer at an "optimal
window" — after liquefaction reaches critical depth, before mesh distortion
becomes fatal. This maps onto the fork's actual strengths with unusual
precision: `LadrunoUP` + PM4Sand/PDMY + DRM + absorbing boundaries is a strong
*initiation* engine, and initiation is the half that MPM codes are weakest at.
Under this architecture the fork's deliverable is **a state-export contract**,
not an MPM solver.

**C — interop only.** Ship nothing MPM-shaped; define an export of
stress/state/geometry at a chosen instant, consumable by CB-Geo or Anura3D. This
is architecture B minus the ambition to own the receiving end — and it is
architecture B's first phase regardless.

### 5.5 Verdict and decision D2

**Recommendation as filed — C now, B as the horizon, A explicitly not.**
*(Superseded by the decision at the end of this section. Retained because its
reasoning names the two obligations that option A now inherits.)*

Reasoning, in order of weight:

1. The fork's differentiator is **initiation physics**, not runout kinematics.
   Runout is a solved, crowded, well-licensed space.
2. **MPM's grid contact routes around the fork's contact engine** (§5.3), which
   is one of its strongest assets. Building A would put two answers to the same
   question in one tree, competing for the same maintenance.
3. CB-Geo is MIT and can be read and reused — so the "we need MPM in-house to
   learn from it" argument does not hold.
4. The §5.1 finding means A stays *cheap to revisit* if the driver ever changes.
   Recording it is most of its value; building it is not.

The counter-argument, stated fairly: if the goal is a **single deck, one
analysis, initiation-through-runout, with the fork's own materials all the way
down**, then B's handoff seam is a real physical approximation (state transfer
loses history and the transfer window is a modelling choice), and only A
avoids it. That is a legitimate ambition — it is just a much larger program than
anything currently on the roadmap, and it should be chosen deliberately.

**DECIDED 2026-08-23 — D2 = A **and** B.** The MPM lane is built in-tree, and it
is driven from the FEM initiation lane through the sequential handoff.

**Stated plainly, after review (§11, finding 5):** the owner chose the larger
program, and reasons 1–3 above were never refuted — they were fenced. The
earlier gloss here ("B is the workflow, A implements its second half") quietly
redefined B, whose §5.4 deliverable was explicitly *a state-export contract, not
an MPM solver*. What A+B actually means is: **the fork owns
initiation-through-runout end to end, accepting the maintenance and the
duplication that reasons 1–3 warned about, in exchange for a single continuous
analysis with the fork's own materials throughout.** That is a legitimate choice;
it is not a refutation of the objections, and the objections are therefore
converted into standing gates rather than closed.

> [!warning] **The affordability premise has weakened since D2 was taken.**
> D2 was decided partly on §5.1's "far cheaper than assumed". Source reading
> (§11, finding 2) reduced that to "integration risk removed; the reusable part
> is the easy ~15 %", and added four structural deltas including *no parallel
> story at all*. The decision is not thereby wrong — but it was taken on a
> stronger premise than now survives, and **should be re-confirmed at M0 rather
> than inherited.**

Two obligations follow directly, and neither is optional:

1. **Contact reconciliation is now a deliverable, not a tension.** §5.3's
   collision — MPM's automatic grid contact vs the fork's explicit frictional
   ALM engine — must be *resolved in writing before M1 code*, because the answer
   determines whether MPM bodies in this fork share a grid (automatic, no-slip)
   or route through `LadrunoContactDomain` (explicit, frictional, and much
   harder). Phase **M0** gates on it.
2. **A maintenance answer versus CB-Geo.** CB-Geo is MIT, tested to ~15,000
   cores, and will stay ahead on runout throughput. The fork's MPM lane must be
   justified by what CB-Geo *cannot* do — run the fork's own materials
   (PM4Sand/PM4Silt/PDMY), inside the fork's own analysis stack, continuous with
   the FEM half. That claim should be written down at M0 and re-tested at every
   phase gate; if it ever stops being true, the lane should stop.

The §5.1 finding is what makes this affordable: the particle store, the
structured background grid, per-step Domain node creation/teardown, particle
convection, particle↔grid interpolation, and the interpreter/recorder plumbing
already exist and already build.

---

## 6. Way forward

Three tracks. **F** is near-term and cheap; **U** and **M** are each multi-ADR,
multi-month programs opened by D1 and D2. They are largely independent and can
run concurrently in separate worktrees.

### 6.1 Track F — harvests and verification (start now)

| Phase | Scope | Depends on | Cost | Gate |
|---|---|---|---|---|
| **F1** | **Oracle corpus O1–O3** (Mandel, Cryer, de Leeuw) as numpy references in the `tests/*_reference.py` idiom; extend ADR-71 §7.1 | none — **unblocked today** | 1 PR, no build | each reproduces its published closed form; `LadrunoUP` matches within a pinned tolerance, or the disagreement is a finding |
| **F2** | ~~X1 fixed-stress split~~ — **WITHDRAWN** (§3.2, §11-1: ADR-73 ships it, measured). Replaced by a one-line docs fix: make ADR-71 §6's three-route list point forward to §12/ADR-73 so the next reader is not misled | F1 | trivial | ADR-71 §6 no longer readable as "fixed-stress absent" |
| **F3** | **X5 constitutive-plugin ABI** study — MFront/MGIS vs `ASDPlasticMaterial3D` vs status quo | independent | medium | own ADR, no code before it |
| **F4** | Housekeeping: Anura3D licence check (§2); `LEDGER_quirks.md` rows for the §5.1 PFEM-BackgroundMesh finding and the Anura3D read-only rule | none | small | rows land |

> [!note] **D1 and D2 raised F3's priority.** The constitutive-interface question
> was previously a nice-to-have. It now has two new consumers arriving at once:
> track U wants **BBM** (suction as an independent state variable — a shape the
> current `NDMaterial` contract does not express), and track M wants **materials
> living on particles** with per-particle commit/revert. Deciding the interface
> *after* both have hard-coded around it is the expensive order. F3 should land
> before U4 and before M1.

### 6.2 Track U — unsaturated poromechanics (proposed **ADR 87**)

| Phase | Scope | Gate |
|---|---|---|
| **U0** | Scoping ADR. **The first real decision: Bishop-type effective stress with χ(S_r), versus net stress + suction as an independent state variable (the BBM route).** Also: which OGS process is the model (`RichardsMechanics`), what the parser surface becomes, air-phase treatment (passive vs active) | adversarial gate — this choice determines whether the `NDMaterial` contract survives untouched |
| **U1** | **X4 property abstraction** — permeability, storage and retention as declared functions of state; `LadrunoUP`'s scalar surface retro-fitted behind it. **Sub-decision: full MPL-style framework vs a minimal `LadrunoPorousMedium`/retention object** — OGS needs MPL for ~15 processes × media × phases; this fork has *one* element family, so the big-code pattern must be priced, not imported | **bit-identical recorder output** on the pinned ADR-71 B1–B5 decks (FP operation order preserved), per the fork's existing flag-off/determinism precedent. Tolerance-based gates apply to O1–O3 separately; do not conflate the two |
| **U2** | Retention + relative-permeability library as first-class objects (van Genuchten, Brooks–Corey), each with its own numpy oracle | oracle agreement; hysteresis explicitly scoped in or out |
| **U3** | `LadrunoUP` unsaturated mode per the U0 decision | OGS `richards-mechanics` benchmarks as a Tier-2 screen (the §4.2 trap applies); analytic infiltration where one exists |
| **U4** | **BBM** as an `nDMaterial` | CODE_BRIGHT worked examples; LC/SI yield-curve behaviour reproduced; wetting-collapse demonstrated |
| **U4b** | **Composition with the existing soil stack** — how unsaturated mode interacts with PM4Sand / PM4Silt / PDMY / `updateMaterialStage`, or an explicit fence ("v1 unsaturated is elastic + BBM only") | written and gated either way; not discovered during U3 |
| **U5** | Validation + guide + ledgers + banner | published guide in the `LadrunoUP_guide` idiom |

### 6.3 Track M — MPM (proposed **ADR 88**)

| Phase | Scope | Gate |
|---|---|---|
| **M0** | **Seam study, the two obligations, and three decisions.** Read `BackgroundMesh.cpp` end to end; map reuse-vs-fork per component against the corrected §5.1/§5.2 lists. **Resolve obligation 1 (contact); write obligation 2 (the CB-Geo justification); decide parallel IN or OUT** — and if OUT, state what obligation 2 can still mean against an MPI code. **Re-confirm D2 on the corrected affordability premise (§5.5 warning).** | written seam map + contact decision + parallel decision + justification + D2 re-confirmation, all adversarially reviewed. **No M1 code before this lands.** |
| **M1** | Single-phase explicit MPM kernel + **numpy oracle, no OpenSees build** (the ADR-71 P0 pattern): particle state (σ, F, V, m + an `NDMaterial` per particle), P2G/G2P transfer, internal force from particle stress, momentum solve, USF/USL/MUSL choice | oracle ≤ 1e-9; then CB-Geo `mpm-benchmarks` (uniaxial stress, hydrostatic column, sliding block on incline — CC-BY-4.0, attributed) |
| **M2** | Shape-function upgrade: GIMP and/or B-spline | **measured** cell-crossing test — stress oscillation as a particle crosses a cell, linear vs upgraded basis, with a number |
| **M3** | OpenSees integration on top of the M0 seam decision: Domain/Element/Recorder wiring, interpreter commands, Δt machinery | end-to-end deck runs from Python and Tcl; recorder output verified |
| **M4** | Locking under near-incompressibility — i.e. **undrained soil, the actual target**. F-bar / TLMPM-class remedy | the ADR-81 limit-load discipline applied to MPM: plateau or no ship |
| **M5** | **Two-phase MPM** — where the lane meets `LadrunoUP`'s physics. Anura3D's published 2-phase/3-phase and the double-point literature as *reading* (§2) | consolidation benchmark reproduced on particles |
| **M1b** | **The state-export contract**, pulled forward from M6 per §5.4-C's own logic ("C is architecture B's first phase regardless"): enumerate what transfers at the handoff and what is lost. Deliverable is a written contract + emitter; it plays to the fork's stated strength and every architecture needs it first | contract published; round-trips into CB-Geo on one deck |
| **M6** | **The B handoff consumed in-tree** — FEM state export → MPM particle initialization; transfer-window study | one deck runs initiation in FEM and runout in MPM continuously |
| **M7** | Contact: implement whatever M0 decided | per M0 |

### 6.4 Track T — TIM macroelement calibration on shipped capability (proposed **ADR 89**)

Added 2026-08-23 after cross-checking this plan against the UANDES/Ladruño **TIMs**
campaign (`tim-macroelement` skill; the vault is authority). It is the only track
here that authors **no new physics** — it applies capability the fork already
ships to a live external campaign, and it closes two of that campaign's named
open questions.

The finding that shapes it: **TIM and TIM-G are drained**, and the SFIM reference
model's measured pore pressure is hydrostatic — so it sits on the *drained
branch*, and P1–P4 need no u-p element at all. Dropping u-p makes every node
ndf = 3, which is exactly what the fork's contact engine requires (§8, the
`ndf==3` guard).

| Phase | Scope | Gate |
|---|---|---|
| **T0** | **Seam check.** `LadrunoBrick -formulation bbar` + `ManzariDafalias` + `LadrunoContactDomain` on the SFIM geometry. All three ship; the combination has never been run | the three compose, or the incompatibility is named and rowed |
| **T1** | **Correctness gate.** Reproduce SFIM's verified vertical response in the drained limit — 5.29 mm at 200 kPa via 1.72 / 2.99 / 4.34 mm, 0.97 mm under gravity alone | matched to the campaign's own tolerance, or the discrepancy is the finding. **Nothing downstream means anything until this passes** |
| **T2** | **Replace the interface layer with real contact.** Gorini's 50-element reduced-strength layer sits at **93.4 % of its own strength under gravity alone**, voiding any leg that must initiate a mechanism in that band. Frictional NTS/mortar contact admitting separation has no layer strength to pre-mobilise | closes campaign **OQ11**; uplift and combined capacities stop being over-estimated |
| **T3** | **Volumetric relief for a real limit load.** [[81_quadratic_hex_limit_load_measurement]] measured linear hex + B-bar at **0.9977 with a plateau**, against the TIMs campaign's own locked linear hex at 0.40 (ADR-81's table cites TIMs directly) | a pushover that plateaus, so P1's "ultimate load" is a limit load and not a displacement threshold on an asymptote |
| **T4** | **The failure campaign** — uniaxial legs then two-component sweeps, per `calibration_workflow.md` steps 4–5; three-mesh bands per the campaign's ADR 65 | ultimate surface fitted, convex, C¹ |
| **T5** *(only if pore pressure builds)* | P5 / TIMg-UP. **Requires lifting the `ndf==3` contact guard** — the u-p + contact seam, its own ADR | gated on U-track and on measured pore-pressure build-up, not assumed |

Two corrections owed back to the TIMs vault, both measured in this fork:

- **The capacity anchor is off by ≈ 2.4×, not ~52×.** [[79_ladruno_geom_hypo_adr]] §9
  measured the collapse load: q = 1108 kPa at s/B = 10 % and 1152 kPa at the
  alternative criterion, against a 1525 kPa Davis (ψ = 0) anchor — 0.73–0.76 of
  it. The over-strength survives at ~2.2×. The mode is **punching**, not general
  shear.
- **`ndf==3`** (§8) is a hard prerequisite for any contact-on-u-p plan.

### 6.5 Sequencing

**F1 starts today** — no build, no decision, and it closes a named verification
hole in the fork's strongest lane. **M0 can start in parallel**: it is reading and
adjudication, not code, and it is the phase most likely to change the shape of
track M. **U0** likewise. F3 should precede U4 and M1.

**T0–T1 are the highest value-per-week in this document** and are independent of
every decision above: they author no physics, need no ADR to land first, and pay
into a campaign that is running now. The author's recommendation is that **T
leads**, with U and M proceeding behind it — a *sequencing* preference, not a
challenge to D1 or D2, both of which the owner reaffirmed on 2026-08-23 with the
long-horizon value in §12 explicitly in view.

Everything after those first phases is gated on their outputs. Nothing in U or M
should be scheduled beyond its first phase until that phase has landed — both
tracks are large enough that a plan drawn now past U1/M1 would be fiction.
---

## 7. Skills — decision D3 (ALL THREE APPROVED)

Per the house pattern (own private repo + Seafile manuals, never binaries in
git), and following [[50_aem_opensees_scoping_report]], whose deliverables were
two skills rather than code.

**DECIDED 2026-08-23 — all three approved.** Sequencing: S1 and S2 are wanted
now (they serve tracks F and U, which start immediately); S3 is best written
*alongside* M0 so the seam study feeds it rather than the other way round.

**S1 — `opensees-poromechanics` (APPROVED — highest value, write first).**
The fork has `LadrunoUP`, `LadrunoPorousOverlay`, ADRs 71/73/78/79 and an
open P7, and **no skill covers any of it.** `quake-research` covers
liquefaction framing, `explicit-dynamics` the integrator side, `kratos` the
U-Pw cross-check, `opensees-expert` the code — nothing covers Biot theory,
consolidation, staggered coupling, storage/permeability conventions, the
drained/undrained envelope, or unsaturated extension. Scope: Biot u-p and u-p-U,
consolidation and its analytic solutions (the §4.3 corpus becomes its reference
set), stabilization and inf-sup, monolithic vs staggered vs fixed-stress,
explicit poromechanics, retention curves and BBM. Manuals: Zienkiewicz et al.
*Computational Geomechanics*, Coussy *Poromechanics*, Lewis & Schrefler, the
**CODE_BRIGHT theoretical manual**, OGS process docs.
Repo `Github\Poromechanics Skill` → `opensees-poromechanics-skills`; Seafile
`opensees-poromechanics-manuals`.

**S2 — `soil-constitutive-models` (APPROVED).**
The symmetric counterpart to `opensees-concrete`, and equally missing. The fork
ships PM4Sand, PM4Silt, ManzariDafalias ±RO, PDMY{,02,03}, PIMY, BoundingCamClay,
CycLiqCP/SP, stressDensityModel, DruckerPrager, and ADR-84 is live work on MC
tension cutoff — with no skill covering selection, calibration, staging
(`updateMaterialStage`), or implementation of any of them. Scope: critical-state
soil mechanics, bounding-surface and multi-yield-surface plasticity, cyclic
mobility, calibration from CPT/SPT/lab, BBM and unsaturated models, and the
fork's own implementation seams. Does not overlap `tim-macroelement` (which
lumps the soil away) or `opensees-concrete`.

**S3 — `material-point-method` (APPROVED under D2 = A+B).**
No longer a keep-the-option-alive skill — with A committed this becomes working
reference material for an implementation track, and should be scoped accordingly:
deeper on formulation and pathology, lighter on survey. Scope: MPM /
GIMP / CPDI / CPDI2 / B-spline bases, TLMPM, USF/USL/MUSL stress updates, the
§5.3 pathology list *with remedies*, two-phase and double-point MPM, the
FEM→MPM handoff literature, and CB-Geo / Karamelo / Anura3D / NairnMPM as
reference implementations — **with §2's licence boundaries stated inside the
skill**, so a future agent reading it does not paste LGPL Fortran into the fork.
Repo `Github\MPM Skill` → `mpm-skills`; Seafile `mpm-manuals`.

Not recommended: a standalone `opensees-pfem` skill — fold the §5.1 findings
into S3 or `opensees-expert`.

---

## 8. Risks and open questions

- **D1 said yes, so U1 is a refactor of shipped, validated code.** The
  byte-identical gate on the saturated battery is the only thing standing between
  a property abstraction and a silent regression in `LadrunoUP`. Do not soften it.
- **The Anura3D licence** (§2) needs a real check before any of X6 turns into
  code, and a `LEDGER_quirks.md` row regardless.
- **Cross-code comparison is not verification** (§4.2). The main way this report
  could do harm is by an OGS number being gated as truth.
- **MPM vs the contact/element-removal line** (§5.3) is now a live design
  problem rather than a strategic tension: D2 chose to build both. **M0 must
  resolve it in writing before any M1 code**, so the fork does not end up with two
  unreconciled answers to "what happens when two bodies touch".
- **Track M's justification is perishable** (§5.5 obligation 2). CB-Geo will stay
  ahead on runout throughput. Re-test the "only the fork can run fork materials
  continuously from initiation through runout" claim at every M gate; if it stops
  being true, stop the lane rather than finishing it out of momentum.
- **Nothing here is measured.** Every claim about the fork's state is grepped
  from this tree; every claim about the three libraries is from their published
  docs and source listings. No benchmark in §4.3 has been run. The first honest
  number arrives at F1.
- **This report's own failure mode, now demonstrated (§11-1):** grepping one
  section of a long ADR and inferring absence. ADR-71 §6's route list was
  superseded by ADR-71 §12 *in the same file*, and the superseding work
  (fixed-stress) was shipped and measured in ADR-73. A confident "genuine hole"
  survived into the executive summary. **Any future claim in this fork of the
  form "X is missing" must cite the ADR's §12 implementation log and the
  companion ADRs, not a keyword grep.**
- **D2's premise moved after the decision** (§5.5 warning). Re-confirm at M0.
- **Track M's serial substrate** (§5.1) versus the fork's SP/MP story and the
  CB-Geo comparison is an unresolved contradiction until M0 decides parallel
  in or out.

---

## 9. References

**Libraries**
- Anura3D — <https://github.com/Anura3D/Anura3D_OpenSource> (LGPL-3.0)
- CODE_BRIGHT — <https://deca.upc.edu/en/projects/code_bright> (binaries; UPC/CIMNE)
- OpenGeoSys 6 — <https://www.opengeosys.org/>, <https://gitlab.opengeosys.org/ogs/ogs> (BSD-3)
- CB-Geo MPM — <https://github.com/cb-geo/mpm> (MIT); benchmarks <https://github.com/cb-geo/mpm-benchmarks> (CC-BY-4.0); Berkeley fork <https://github.com/geomechanics/mpm>
- Karamelo — <https://github.com/adevaucorbeil/karamelo> (GPL); de Vaucorbeil et al., *Comput. Particle Mech.* 8:767 (2021)
- NairnMPM — <https://osupdocs.forestry.oregonstate.edu/index.php/NairnMPM>

**OGS specifics**
- Staggered consolidation benchmark — <https://www.opengeosys.org/6.5.7/docs/benchmarks/hydro-mechanics/consolidationbenchmark/>
- MFront/MGIS interface — <https://www.opengeosys.org/docs/userguide/features/mfront/>

**PFEM in OpenSees (§5.1)**
- Zhu & Scott, "Modeling fluid–structure interaction by the particle finite
  element method in OpenSees", *Computers & Structures* 132 (2014)
- Zhu & Scott, "A PFEM background mesh for simulating fluid and frame structure
  interaction", *J. Struct. Eng.* 148(6) (2022)
- In-tree: `SRC/element/PFEMElement/{Particle,BackgroundMesh,BackgroundGrid,BCell,BNode}.{h,cpp}`

**MPM method**
- Wilson et al., "Distillation of the MPM cell crossing error…", *IJNME* (2021)
- Hybrid FEM→MPM for tailings-dam failure — arXiv:2412.08040
- Liquefaction-induced dam failure, a case for MPM — arXiv:2111.13584
- Stabilised semi-implicit double-point MPM for soil–water coupling — arXiv:2401.11951

**Fork**
- [[71_ladruno_up_family_adr]], [[73_ladruno_porous_overlay_adr]],
  [[78_ladruno_up_corot_adr]], [[79_ladruno_geom_hypo_adr]],
  [[50_aem_opensees_scoping_report]], [[40a_kratos_crosspollination_amendment]],
  [[LadrunoUP_guide]]

---

## 10. Implementation log

*(empty — no code authored under this ADR)*

---

## 11. Adversarial review log

### 2026-08-23 — Fable panel, post-decision gate. **NOT CLEAN — one blocker, three majors.**

Reviewer was pointed at seven load-bearing claims and instructed to verify
against the tree rather than the document. Findings, and disposition:

| # | Finding | Verdict | Disposition |
|---|---|---|---|
| 1 | **X1/F2's "fixed-stress hole" does not exist** — ADR-73 *is* the fixed-stress ADR (shipped, measured ≡ monolithic BE, ×3.3 cheaper factorization; Kim–Tchelepi–Juanes, Mikelić–Wheeler, Turska–Schrefler already cited). ADR-71 §12 records P7 as opened by ADR-73 | **BLOCKER — REFUTED** | X1 withdrawn (§3.2); F2 rewritten to a docs fix (§6.1); §0, §1 corrected; failure mode recorded in §8 |
| 2 | **PFEM reuse overstated** — `BackgroundGrid` is dead code; P2G/G2P are non-conservative SPH-average / PIC on different bases; grid nodes ndf-locked to `Pressure_Constraint`; **no parallel story** (static singleton, no `sendSelf`, in-source "Can't be in parallel" TODO); analysis-side lifecycle omitted; "far cheaper" never quantified | **MAJOR — WEAKENED** | §5.1 requantified ("integration risk removed, easy ~15 %"); four deltas added to §5.2; **D2 flagged for re-confirmation at M0** (§5.5) |
| 3 | **Oracle inference refuted** — B1 is periodic and B5 pins oscillation, so "all monotonic" is false; ADR-71 P0 already gates Q vs numpy at 1e-9 and B3 gates the Skempton p(0⁺) = 0.410 MPa; ADR-73 P0 already exercised the Mandel–Cryer dip | **MAJOR — REFUTED as stated** | §4.3 rewritten; the surviving gap restated as *no multi-dimensional closed-form coupled pressure gate*; O1–O3 retained on that basis only |
| 4 | **Track M internally inconsistent** — the state-export contract (which §5.4-C calls "B's first phase regardless") was scheduled at M6, last; and no phase delivers parallel although obligation 2 benchmarks against a 15k-core MPI code | **MAJOR — REFUTED** | export pulled forward to **M1b**; parallel in/out added to the M0 gate |
| 5 | A+B coherence gloss was partly rationalisation; obligations are real but structurally hard to honor given finding 4 | MINOR — WEAKENED | §5.5 paragraph rewritten to say plainly that the objections stand and are converted to gates |
| 6 | "X4 before X3 **without exception**" asserts necessity where there is prudence; MPL may be over-scaled for a one-element-family fork; U1's gate conflates byte-identity with tolerance gates | MINOR — WEAKENED | §3.3 softened; U1 gains an MPL-vs-minimal sub-decision and a restated bit-identical-recorder-output gate |
| 7 | The licence file is **`COPYRIGHT`**, not `LICENSE` — a tell that §2 was written from memory of upstream rather than this tree | MINOR | corrected §2 |
| 8 | Track U silent on the existing soil constitutive stack; U3 gated on a Tier-2 number while Tier-1 unsaturated closed forms exist and go unnamed | MINOR | **O8a** added (Srivastava–Yeh, Tracy, Philip); **U4b** added (composition or explicit fence) |

**Reviewer's overall verdict:** *"Not safe to act on as written."* F1 safe once
finding 3's justification is corrected; **F2 must not start**; D2 rests on a
premise source-reading does not support at the stated strength.

**Claims the reviewer checked and cleared:** §2's grant terms and the LGPL
incompatibility mechanism (including the derivative-work claim for hand-
translation and the adequacy of the hedging); §1's no-unsaturated grep and the
soil-stack inventory; ELE tag 33017; §4.1's 335 test files and the
`ladruno_up_reference.py` independence protocol at ≤1e-9; §4.2's Tier doctrine;
§4.3's narrow grep (Mandel/Cryer/de Leeuw/Booker–Small/Gibson genuinely absent
from ADR-71); §5.1's raw API inventory and interpreter exposure; §5.3's pathology
list; §3.1's external facts on all three libraries plus CB-Geo; §7's skill-gap
inventory.

**Author's verification.** Findings 1, 2 (`BackgroundGrid`), 3 (ADR-73 dip) and
7 were independently re-checked in-tree before acceptance, not taken on the
reviewer's word. All four confirmed.

---

## 12. What this program unlocks

The near-term phases justify themselves on their own gates. This section answers
the larger question — *what becomes possible that is not possible today* — because
D1 and D2 were reaffirmed on it and it should be written down rather than assumed.

### 12.1 Today's two ceilings

The fork computes **initiation** well: when and where soil fails, under real
seismic input, with pore-pressure generation, frictional contact and element
removal. Two things it cannot compute at all:

1. **Anything above the water table.** `LadrunoUP` is saturated Biot by decision.
2. **Anything after the mesh distorts.** Every FE lane ends where the elements
   tangle — which is at, or just after, the failure the analysis was run to find.

### 12.2 Track U removes the water-table ceiling

- **Rainfall-induced landslides** — the largest geohazard class the fork cannot
  touch today. The trigger *is* suction loss on wetting; without retention
  behaviour there is no mechanism to model.
- **Tailings dams** — the beach and the crest are partially saturated. There is
  no honest tailings model without unsaturated behaviour, and tailings is where
  the liquefaction-plus-runout problem actually bites.
- **Compacted fill and wetting collapse** — BBM's home ground.
- **Embankment dams above the phreatic surface**; rapid drawdown with a real
  unsaturated zone rather than an assumed one.
- **Seasonal and climate-driven serviceability** — expansive soils, shrink–swell,
  pavement subgrades.
- The **infiltration coupling** generally: the hydrology side of geotechnics,
  which the fork currently cannot reach.

### 12.3 Track M removes the post-failure ceiling

Risk is probability × consequence. **The fork computes the first factor and
cannot compute the second at all.** Track M is the consequence half:

- **Landslide and debris-flow runout** — where the material goes, how far, how fast.
- **Tailings-dam breach inundation extent** — the number that governs downstream
  risk, and the reason the hybrid FEM→MPM literature exists.
- **Liquefaction-induced lateral spreading to large displacement**, past where a
  Lagrangian mesh survives.
- **Pile driving, penetration, and installation effects.**
- **Progressive collapse with debris** — combined with the shipped contact and
  element-removal stack, this is the [[50_aem_opensees_scoping_report]] / AEM
  target reached by a different road.

### 12.4 Together: the thing nothing else does

U and M are individually useful and jointly differentiating. The combination
unlocks **one continuous analysis from ground motion → pore-pressure generation →
triggering → failure-surface formation → runout → deposition, with the same
materials and the same code throughout.**

Nothing in the surveyed field does both ends:

| code | initiation | runout |
|---|---|---|
| Anura3D | weak — three constitutive models total (§3.1) | yes |
| CB-Geo MPM | weak — no seismic input stack | yes, at scale |
| OpenGeoSys | coupled flow, but no seismic and no large deformation | no |
| CODE_BRIGHT | THM, geological timescales | no |
| FLAC / PLAXIS | strong | stops where the mesh distorts |
| **this fork, after U + M** | **strong — PM4Sand/PDMY + DRM + absorbing boundaries + contact** | **yes** |

That is the strategic claim, and it is the one §5.5 obligation 2 requires be
re-tested at every gate: the fork's MPM lane is justified by *continuity with its
own initiation physics*, never by out-running CB-Geo.

### 12.5 Track T unlocks the near-term one

A calibrated macroelement turns an entire continuum campaign into a handful of
parameters. That is the only way soil–structure interaction survives into a
performance-based assessment — nobody runs a soil continuum under thirty ground
motions. **Track T is the bridge from the fork's physics to design practice**,
and it uses only what already ships.

### 12.6 What this does NOT unlock — stated so the ledger is honest

- **Nothing here makes the fork faster.** That is ADR-40/75's lane.
- **Nothing here helps the concrete or steel side.**
- **Every track is permanent maintenance.** Three new lanes, two new physics
  domains and three new skills all compete for the same reviewer, and the fork's
  own history (§11) shows what happens when a lane is planned faster than it is
  read.
