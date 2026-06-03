# RC 3D Modeling — End-to-End Recipe & Gap Map

**Status:** scoping doc (no code). Purpose: walk the *complete* pipeline for a
reinforced-concrete 3D **solid** model using only what exists on `ladruno`
today, marking every step ✅ works / ⚠️ works-but-manual / ❌ missing. The point
is to expose exactly where the gaps bite before we write any new C++.

**Scope target:** cyclic RC solid models — confined columns, shear walls,
beam-column joints — under **both implicit and explicit** analysis.

**Companion docs:** [[09_ladruno_brick]] (concrete host), [[11_brick_asdconcrete_integration]]
(crack-band handshake), [[14_ladruno_rebar_buckling_adr]] + [[13_ladruno_uniaxial_j2_adr]]
(rebar), [[10_ladruno_j2_plasticity]] (continuum steel/J2). Scoping verdict +
embedment landscape recorded in agent memory `project_rc3d_embedment_scoping`.

---

## TL;DR — readiness

| Pillar | Status | Note |
|---|---|---|
| 3D solid elements + large-disp geometry | ✅ production | `LadrunoBrick`, `BezierTet10` |
| Reinforcing steel (uniaxial) | ✅ production | `LadrunoUniaxialJ2` + `LadrunoRebarBuckling` |
| 3D concrete continuum | ⚠️ usable-with-care | `ASDConcrete3D`: manual hardening laws + lch; **dilatancy is the hidden knob** |
| **Rebar ↔ concrete coupling** | ❌ **critical gap** | perfect-bond only; no 1D bond-slip; embedding is penalty (implicit-only) |
| **Rebar layout / embedding tooling** | ❌ **critical gap** | no generator; `generateInterfacePoints` is Tcl-only |
| Solver (cyclic + softening) | ⚠️ usable-with-care | everything present; adaptive stepping is DIY (Python) |
| Post-processing | ⚠️ manual | recorder ships; apeGmsh `.ladruno` reader not built |

Two ❌ gaps stand between today and a credible cyclic joint/column:
**(1) a 1D bond-slip mechanism that survives explicit**, and **(2) tooling to
lay out and embed the rebar cage**.

---

## The pipeline, step by step

### Step 1 — Mesh the concrete  ✅ / ⚠️
- **Hex:** `LadrunoBrick` via apeGmsh direct-drive meshing. ✅
- **Tet:** `BezierTet10` (quadratic Bézier tet, has crack-band `getCharacteristicLength()`
  + F-bar). ✅ — easier to fill an arbitrary rebar cage.
- ⚠️ `LadrunoBrick` does *not* yet expose `getCharacteristicLength()` (BezierTet10 does);
  with `ASDConcrete3D` on hex you currently lean on the element's inter-node default
  or set the material length-scale by hand. *(Quick win: lift the override onto LadrunoBrick.)*

### Step 2 — Concrete material  ⚠️
```python
# ASDConcrete3D needs explicit tension/compression hardening laws + lch_ref.
ops.nDMaterial('ASDConcrete3D', matTag, E, nu,
               '-Te', *Te, '-Ts', *Ts,        # tension strain/stress points
               '-Ce', *Ce, '-Cs', *Cs,        # compression strain/stress points
               '-lch_ref', lch_ref, ...)       # crack-band reference length
```
- ✅ The model itself is capable (Lubliner triaxial, plastic-damage, crack planes, crack-band).
- ⚠️ **No plug-and-play `fc/ft/Gf` parametrization** — you build the hardening curves yourself.
- ⚠️ **Confinement: do NOT pre-inflate `fc` (Mander) in 3D** — the hoops develop lateral
  pressure and the triaxial surface responds. BUT whether ASDConcrete3D actually delivers
  emergent confined-strength gain is **uncertain and must be tested** (ADR 20 Gate 2): it has
  a Lubliner plastic-damage surface with a Kc triaxial parameter, but **NO tunable dilation
  angle**, and the backbone peak is set by the user's compression curve — so confined `fcc`
  may not rise without calibrating the backbone per confinement level. The relevant levers are
  **Kc + the compression backbone**, NOT a dilation angle (there is no dilation input).
  Validate via a standalone triaxial-vs-Mander unit test before relying on it.
- ✅ Damage-scaled hourglass (Tier-A) + crack-band handshake already wired for LadrunoBrick
  (`-formulation eas|uri-stiffness`, NORMAL analysis only — see [[11_brick_asdconcrete_integration]]).

### Step 3 — Lay out the rebar cage  ❌
- ❌ **No layout generator.** Longitudinal bars, ties/hoops, cover, spacing, lap positions
  must all be placed by hand as line geometry. This is the workflow wall.
- Target tooling: parametric apeGmsh helper — given (member type, dims, bar layout, tie
  spacing, cover) → emit bar/tie line elements + embedding + bond automatically.

### Step 4 — Rebar elements + material  ✅
```python
# Bar material: combined-hardening J2 + buckling wrapper
ops.uniaxialMaterial('LadrunoUniaxialJ2', barCore, E, fy, ...)         # Chaboche ratcheting
ops.uniaxialMaterial('LadrunoRebarBuckling', barTag, barCore, '-lsr', s_over_d, ...)
# Bar element: corotational truss (large-disp) carrying the bar material
ops.element('corotTruss', eleTag, ni, nj, A, barTag)
```
- ✅ Steel side is best-in-class (cyclic ratcheting, Dhakal–Maekawa/Gomes–Appleton buckling,
  re-straightening, low-cycle fatigue via composition).
- Note: the bar is a **1D axial** law — dowel action / multi-axial bar stress is out of scope v1.

### Step 5 — Embed the rebar in the concrete  ❌ (the crux)
Today's options, all deficient for the explicit+implicit requirement:

| Option | How | Implicit | Explicit | Bond |
|---|---|---|---|---|
| **Matching mesh + `equalDOF`** | rebar nodes coincide with concrete nodes | ✅ | ✅ | perfect only; needs conformal mesh (cage ⇒ infeasible) |
| **`ASDEmbeddedNodeElement`** | tie rebar node to projected triangle (3 nodes) / tet volume (4 nodes) | ✅ | ❌ penalty K≈1e18 collapses `dt_cr` | perfect only |
| **`generateInterfacePoints` → `EmbeddedBeamInterfaceP/EP`** | auto point-location + beam-in-solid | ✅ | ❌ penalty 1e12 / Lagrange | EP uses a *3D soil-pile* interface law, not a 1D rebar bond |
| | | | | **`generateInterfacePoints` is Tcl-only — unreachable from openseespy** |

```python
# What works TODAY (implicit, perfect bond, non-matching, rebar node into a concrete tet):
ops.element('ASDEmbeddedNodeElement', tag, rebarNode, c1, c2, c3, c4)   # +R4 ⇒ tet volume
#   c1..c3            ⇒ projected onto that triangular surface
#   '-K', K           ⇒ penalty (default 1e18) — fine implicit, fatal for explicit dt
```
- ❌ **No bond-slip** anywhere — every path is perfect bond. Joints/splices/hinge zones are
  dominated by bond degradation, so perfect bond is *qualitatively* wrong.
- ❌ **No explicit-safe embedding** — all are penalty/Lagrange.

### Step 6 — Analysis  ✅ / ⚠️
```python
# Implicit cyclic pushover
ops.constraints('Transformation')          # needed if we move to transformation-MPC embedding
ops.integrator('DisplacementControl', ctrlNode, dof, d_incr)
ops.algorithm('NewtonLineSearch')          # or KrylovNewton / BFGS for softening
# Explicit collapse / dynamic
ops.integrator('CentralDifferenceLadruno') # or 'ExplicitBathe'
```
- ✅ Implicit (arc-length/DisplacementControl + line-search/Krylov/BFGS) and explicit
  (CentralDifferenceLadruno / ExplicitBathe + energy balance) all present.
- ⚠️ **Adaptive stepping is DIY** — wrap `analyze()` in a halve-on-fail / grow-on-success
  Python loop. This is the *correct* place for it (not a new integrator). A blessed,
  tested wrapper proc should ship as fork tooling.

### Step 7 — Record & post-process  ✅ / ⚠️
- ✅ `LadrunoRecorder` (.ladruno HDF5: GP stress/strain, energy, local axes, envelopes);
  fiber/element recorders for rebar stress.
- ⚠️ **apeGmsh `.ladruno` reader / RC post-processor not built** → custom h5py for damage
  fields + per-bar stress histories. Rebar (fiber) results are not a first-class recorder family.

---

## Where the gaps bite — ranked

1. **❌ No 1D bond-slip that survives explicit** — the single missing physical mechanism.
2. **❌ No rebar-cage layout + embedding tooling** — the workflow wall (and the auto-embed
   helper that does exist, `generateInterfacePoints`, is stranded in Tcl).
3. **⚠️ Concrete triaxial/dilatancy not validated** — confinement is emergent, so the model
   is only as good as ASDConcrete3D's `fcc(p)` response and dilation angle. Needs a calibration
   battery, not a confinement "feature".
4. **⚠️ Adaptive stepping + post-processing** are DIY — tooling, not source.

---

## Design direction (decided 2026-06-03)

The reference codes (Abaqus `*EMBEDDED ELEMENT`, LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID`,
DIANA bond-slip reinforcement) converge on: reinforcement keeps its own DOFs, located in
host elements; perfect bond = host-interpolation tie; bond-slip = a **1D τ–s interface**
along the bar axis. None pre-inflate `fc` (confinement emergent; dilatancy governs).

Target architecture for the fork:

- **One embedded-reinforcement element, two modes — and the mode split is FORCED by the
  solver, not just a preference:**
  - *perfect bond, implicit* via **transformation / DOF-elimination** (Abaqus-style) — adds no
    stiffness, no penalty parameter. The clean implicit path.
  - *explicit* MUST use the **penalty / bond-slip mode**. `CentralDifferenceLadruno` does a
    trivial **diagonal** M⁻¹ solve (`system Diagonal`, [CentralDifferenceLadruno.cpp:186]).
    A transformation MPC condenses rebar DOFs into host DOFs and **breaks mass diagonality** →
    kills the O(N) explicit solve. Penalty coupling keeps mass diagonal (it only adds stiffness),
    which is exactly why LS-DYNA CBIS uses penalty for explicit. Cost: the penalty stiffness +
    stiff/short rebar throttle `dt_cr` → mitigate with **rebar mass scaling** (standard practice).
  - So: transformation = implicit perfect-bond; penalty(+bond-slip, mass-scaled) = explicit.
    Both modes are *required*, not optional.
- **Host-agnostic** via a generic isoparametric inverse-map → works in **both** BezierTet10 (tet)
  and LadrunoBrick (hex).
- **New 1D bond-slip τ–s uniaxial material** (CEB-FIP backbone, cyclic degradation) dropped
  into the element's bar-axial slot; transverse compatibility tied kinematically.

Two genuinely-new pieces: **(a)** the transformation-based embedded-reinforcement element
(the missing Abaqus primitive, explicit-safe), **(b)** the 1D bond-slip material. Everything
else — host elements, rebar materials, concrete, solvers — already exists.

### Suggested order
1. **This recipe** validated end-to-end on a perfect-bond confined column (implicit), to
   surface any remaining plumbing surprises. *(quick win: lift `getCharacteristicLength()`
   onto LadrunoBrick first.)*
2. **Triaxial/dilatancy validation** of ASDConcrete3D (`fcc(p)` vs test data; dilation sweep).
3. **ADR + build:** embedded-reinforcement element (both modes, host-agnostic) + 1D bond-slip material.
4. **Tooling:** apeGmsh rebar-layout generator (emits bars + embed + bond) + `.ladruno` RC reader.
5. **Blessed Python adaptive-analyze wrapper.**
