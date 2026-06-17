---
title: LadrunoDispBeamColumn — regularized displacement-based frame (per-IP lch + embedded softening hinge)
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - adr
  - element
  - frame
  - regularization
  - softening
  - embedded-discontinuity
  - eas
  - large-displacement
  - research
---

# LadrunoDispBeamColumn — regularized displacement-based frame

> Siblings: [[14_ladruno_imk_beam]] (the **lumped**-plasticity end-hinge beam; this is its **distributed**-plasticity, fiber-section, regularized cousin — they share the cohesive-law idea but NOT the condensation algebra, see §"False reuse" risk), [[19_ladruno_brick_eas_simo_rifai]] + [[20_ladruno_brick_eas_stabilization]] (the EAS / static-condensation machinery this element ports), [[26_ladruno_plane_frontier_adr]] (the regularized crack-band materials whose softening this element is built to carry on a frame), [[31_ladruno_robust_solve_driver_adr]] (`ladruno_drive`, the snap-back solve orchestration this element's Stage-2 gates run under), and the FB analog `RegularizedHingeIntegration` (Scott & Hamutçuoğlu 2008). This is a **decision document**: it lands on a two-tier element, sequences what ships in v1, and pins the corotational-composition contract that one scoping pass got wrong.

## What

A new fork-authored, distributed-plasticity, fiber-section **displacement-based** frame element, `LadrunoDispBeamColumn` (2D + 3D), that handles material softening with proper regularization and supports large displacement. It ships in **two tiers on one element**:

- **Tier-1 — `lch` channel (the cheap floor, must-ship).** Mirror `ForceBeamColumn` exactly: store the per-IP tributary length `current_section_lch = wt[i]·L` in the `update()` integration loop and override `getCharacteristicLength()` to return it. This is the fix stock `DispBeamColumn` never got, and it is the *only* piece needed to make existing crack-band / auto-regularizing materials (`ASDConcrete1D -autoRegularization`, `ASDSteel1D`, `LadrunoUniaxialJ2`+`LadrunoLemaitreDamage`) mesh-objective on a displacement-based fiber frame.
- **Tier-2 — embedded strong-discontinuity hinge (the robust endpoint).** An enhanced-kinematics rotation jump `κ = B·d + G·α` with `α` (the rotation jump) statically condensed at the element level, governed by a **discrete cohesive moment–rotation law `M([[θ]])` carrying fracture energy per hinge** (units of energy — **no `lch` to calibrate**). This is the Armero–Ehrlich (2006) / Jukić–Brank–Ibrahimbegović (2013) construction. It is mesh- *and* integration-objective, and removes the residual integration-rule sensitivity Tier-1 still carries.

**In scope (v1):** Tier-1 in 2D + 3D with large displacement via the existing Corotational coordinate transform; Tier-2 in **2D only**. **Out of scope / deferred:** the 3D embedded hinge (biaxial/torsional jump on the `CorotCrdTransf3d` quaternion-triad — its own later ADR), DDM sensitivity for the condensed-`α` element, and any nonlocal/gradient (internal-length) regularization (a different and heavier formulation, not pursued for a line element).

## Why

Stock `DispBeamColumn` **silently** produces mesh-dependent softening response. It does **not** override `getCharacteristicLength()`, so it inherits `Element::getCharacteristicLength()` ≈ the full element length (`SRC/element/Element.cpp:682`, returns the min inter-node distance). Auto-regularizing materials therefore smear the softening branch over the *whole element* instead of the *localizing integration point*, so `Gf` is regularized against the wrong length and the global response drifts with both mesh size and integration-point count. `ForceBeamColumn` solved this years ago (`current_section_lch = wt[i]*L`, `SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp:1357`, returned at `:3464`); the displacement-based sibling was never given the same fix. The pathology is documented in [[LEDGER_quirks]] (§"Crack-band materials read element size via `ops_TheActiveElement->getCharacteristicLength()`", quirks ledger lines 59–90).

Beyond the cheap fix, the *theoretically correct* endpoint for a softening displacement-based frame is not "scale the constitutive law by a length" (which stays objective only inside a size window and can trigger material-level snap-back below a minimum element size) but to **add back the kinematic mode the displacement element is missing** — a rotation jump (a kink) — and let a discrete cohesive law carry the fracture energy directly. The displacement-based formulation is the *natural* home for this enrichment (the defect is a missing kinematic freedom; the cure adds exactly that freedom), whereas the force-based element's analog is the integration-level `RegularizedHingeIntegration`. The fork already owns the EAS + static-condensation machinery to build it, in `LadrunoBrick` and `ASDShellQ4`.

This element is the distributed-plasticity workhorse meant to pair with the fork's RC material stack (ASDConcrete1D + rebar buckling + bond-slip); the lumped IMK beam ([[14_ladruno_imk_beam]]) remains the concentrated-plasticity option.

## Where

- **New code:**
  - `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn2d.{h,cpp}` — Tier-1 + Tier-2.
  - `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn3d.{h,cpp}` — Tier-1 only in v1.
  - `SRC/element/ladrunoDispBeamColumn/OPS_LadrunoDispBeamColumn.cpp` — one `ndm`-dispatching factory (the `LadrunoIMKBeam` pattern), parsing nodes, transf, beam integration, section(s), `-lch` mode, hinge options.
  - `SRC/element/ladrunoDispBeamColumn/CMakeLists.txt` (+ parent `add_subdirectory`).
  - `SRC/material/uniaxial/LadrunoCohesiveHinge.{h,cpp}` (Tier-2) — energy-based discrete `M([[θ]])` cohesive law, `Gf_hinge` in, near-rigid penalty with a guarded floor, `getEnergy()`. The element accepts **any** `UniaxialMaterial` in the hinge slot; this is the default.
- **Reference (copy patterns from, read-only):**
  - `SRC/element/dispBeamColumn/DispBeamColumn2d.cpp` (skeleton, `update()` IP loop ~`:655-679`, B-operator ~`:670`), `DispBeamColumn3d.cpp`, `DispBeamColumnNL2d.cpp` (the `½θ²` membrane/bowing basic strain).
  - `SRC/element/dispBeamColumn/TimoshenkoBeamColumn2d.cpp:610` (interdependent-interpolation anti-locking B-operator, `phi = 12EI/(GA_s L²)`).
  - `SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp:1354-1357,3464-3473` (the `lch` channel to mirror **verbatim**).
  - `SRC/element/Element.cpp:682` (the default to override away from).
  - `SRC/element/ladrunoBrick/LadrunoBrick.cpp:2657-2735` (`formEAStrue`: inner-Newton-on-`α` + static condensation **algebra**) and `:362-412,2823-2930` (committed enhanced state + fixed-layout gated serialization).
  - `SRC/element/shell/ASDShellQ4.cpp:2104-2139` (EAS condensation by explicit `Invert`; the `EASData` persistent-state container + conditional-send discipline).
  - `SRC/coordTransformation/CorotCrdTransf2d.cpp:546,591,935` (`getGlobalResistingForce(pb,p0)`, `getGlobalStiffMatrix(kb,pb)`, `getGeomStiffMatrix(pb)` — the **3-DOF-only** basic-system contract, see the pinned invariant in §How).
  - `SRC/material/uniaxial/ASDConcrete1DMaterial.cpp:997-1000` (the once-only `lch` latch on first `setTrialStrain`).
- **Modify (fork plumbing, same PR per CLAUDE.md):**
  - `SRC/classTags.h` (tags added **only on merge**, see Reserved tags below).
  - `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` (null-ctor cases — missing this crashes parallel/DB `recvSelf`).
  - interpreter registration at **both** sites: `SRC/interpreter/OpenSeesElementCommands.cpp` (Python) and the Tcl element command table.
  - `Ladruno_scripts/banner_features.txt` → `python Ladruno_scripts/patch_banner.py` (never hand-edit the C strings).
  - `Ladruno_implementation/LEDGER_implementations.md` row; `LEDGER_quirks.md` cross-ref resolving the `:59` pathology; `LEDGER_vanilla_files.md` row **only if** `MPCORecorder.cpp` is touched (see open question Q1).

### Reserved tags

Per-registry bands are independent (Element / uniaxial / Integrator each have their own 33000-space), so a number reused across registries is not a collision.

- `ELE_TAG_LadrunoDispBeamColumn2d = 33013` — RESERVED, not yet built. Free in the **Element** registry (used: 33000,02–08,11,12; 33009/33010 reserved by [[26_ladruno_plane_frontier_adr]] for VEM/SBFEM; 33016 reserved for LogStrain2D). Numerically equals `ND_TAG_InitDefGradNDMaterial 33013` — **not** a collision (different registry, per the ledger rule).
- `ELE_TAG_LadrunoDispBeamColumn3d = 33014` — RESERVED, not yet built (Element registry; numerically equals `ND_TAG_StagedStrainNDMaterial 33014`, again not a collision).
- `MAT_TAG_LadrunoCohesiveHinge` — RESERVED, number TBD against the uniaxial band at implementation.

## How

### Kinematics (default Timoshenko)

Use shear-flexible **Timoshenko** kinematics with the interdependent-interpolation anti-locking B-operator copied from `TimoshenkoBeamColumn2d.cpp:610` (`phi` **frozen at the initial elastic value** and documented as a fixed interpolation parameter, not a tracked physical quantity). Euler–Bernoulli is recoverable as the `phi → 0` (rigid-shear) limit, so it is not a separate element.

### Tier-1 — the `lch` channel

In the `update()` IP loop, set `current_section_lch = wt[i] * crdTransf->getInitialLength()` (the **reference** length — frame-invariant; using the deformed length silently de-calibrates `Gf` as the element rotates/stretches) **immediately before** each `theSections[i]->setTrialSectionDeformation(e)`. Override:

```cpp
double LadrunoDispBeamColumn2d::getCharacteristicLength(void) {
    return (current_section_lch > 0.0) ? current_section_lch
                                       : Element::getCharacteristicLength();
}
```

`-lch {ip|element|<value>}`: `ip` (per-IP tributary) is the **only safe default**; `element` (full `L`) is a guarded A/B-debug option that **must emit a loud `opserr` warning** that it re-enables multi-IP energy double-counting; `<value>` pins a user crack-band width.

### Tier-2 — embedded strong-discontinuity hinge

Enhanced curvature `κ = B·d + G·α`, where `G` is a (regularized-)Dirac enhancement concentrated at **one** hinge section per element. `α` is element-internal, converged inside `update()` and statically condensed:

```
K_basic = K_dd − K_dα · K_αα⁻¹ · K_αd      (3×3, basic system)
pb_basic = condensed 3-vector basic force
```

Borrow the inner-Newton + condensation **algebra** from `LadrunoBrick::formEAStrue` and the persistent-state/serialization discipline from `ASDShellQ4` — but **not** LadrunoBrick's `globalizeStiff` seam (see pinned invariant). `G` **must** satisfy the EAS/incompatible-modes patch-test orthogonality (`∫G dx = 0` against constant section forces) so the no-hinge state byte-reduces to Tier-1. The hinge is the **softening branch of the same fiber section at the localizing IP** (Armero–Ehrlich consistent), not a parallel spring — so the fiber section must be frozen / unload elastically once the hinge opens (see landmine #1).

### PINNED INVARIANT — corotational composition

> One scoping pass framed corotational composition as a low-risk drop-in ("operates in the basic system exactly like LadrunoBrick routes through `globalizeStiff`"). **That is a false analogy and must not be repeated.** `CorotCrdTransf2d::getGlobalStiffMatrix(kb,pb)` (`:591`) and `getGlobalResistingForce(pb,p0)` (`:546`) accept **only a 3×3 basic stiffness and a 3-vector basic force**, and the transform **owns and builds its own geometric stiffness from `pb`** (`getGeomStiffMatrix(pb)`, `:935`). There is **no seam for element-internal DOFs.** LadrunoBrick's `globalizeStiff` globalizes a *full condensed 24-DOF* stiffness in a co-rotated CORE frame the element owns — a fundamentally different mechanism.
>
> Therefore: **condense `α` to the 3-DOF basic `K_basic`/`pb` BEFORE calling `crdTransf`.** Frame-invariance is inherited **only because `α` is a basic (corotated, rigid-body-free) rotation jump whose cohesive `M([[θ]])` law is work-conjugate to the basic moment** — this is a load-bearing assumption, not an incidental one.

### Large-displacement decision

In scope for v1 via the existing **Corotational** transform (`CorotCrdTransf2d` for Tier-2; `CorotCrdTransf2d/3d`, `PDelta`, `Linear` all free through the `CrdTransf` API for Tier-1). **Do not** add element-internal geometric stiffness. Default basic strain: the `DispBeamColumnNL2d`-style `½θ²` bowing term, because a softening column under large displacement is exactly where P-Δ + bowing drive snap-back and the purely linear strain does not recover that coupling.

### Four explicit design rules (the correctness landmines)

1. **No energy double-counting.** Tier-1 smeared regularization and the Tier-2 discrete hinge are *mutually exclusive* localization mechanisms, not additive layers. When the hinge opens, the fiber section at that IP **must** freeze / unload elastically, or `Gf` is counted twice. Mandatory mutual-exclusion switch + an energy-balance regression (`total dissipation == prescribed Gf`). Corollary: a localizing section whose fibers are themselves softening (e.g. ASDConcrete1D) **cannot** also carry a Tier-2 hinge — validate non-softening fibers at the hinge IP in `setDomain`.
2. **`K_αα` passes through zero at the cohesive peak.** This is the activation event *every* hinge undergoes, not an edge case. Do **not** port LadrunoBrick's residual-Newton `Kaa.Solve` verbatim — use a **closed-form / return-mapping jump update** (jump as primary unknown). ADR-20 already found scalar `β·Kaa` stabilization is not the cure.
3. **The `lch` in-loop assignment is load-bearing.** Auto-regularizing materials latch `lch` **once** on first `setTrialStrain`; all `N_fiber × N_IP` copies latch during the single `setTrialSectionDeformation` while `current_section_lch` holds that IP's value. Move it out of the IP loop (e.g. to `setDomain`) and every fiber silently gets the wrong band. Pin with an in-loop-assignment invariant + a white-box probe proving section-`i` material received `wt[i]·L` (not min-node-distance) on its first read. Defensively set `ops_TheActiveElement = this` at the top of `update()`.
4. **The condensation reuse is algebra-only.** IMK ([[14_ladruno_imk_beam]]) condenses against a *constant analytic* `4EI/L` interior; this condenses `α` against the *nonlinear, numerically-integrated, possibly-indefinite* fiber-section tangent of the localizing IP. There is no closed-form `F_el + F_h` inverse — re-cost the cohesive law + homogenization-at-activation rule to ~5–6 days. The *brick's* algebra ports; the IMK "reuse verbatim" premise does not.

### Staged plan

| Stage | Goal | Exit gate |
|---|---|---|
| **0 — Skeleton (2D)** | Timoshenko fiber beam that byte-reduces to stock; lock the verification harness before any new physics | FD-tangent (1e-6 rel) ✓; bit-identical to stock `DispBeamColumn2d` with Linear transf + no hinge ✓; rigid-body-rotation objectivity of the elastic corotational path ✓ |
| **1 — Tier-1 `lch` + large-disp (2D+3D)** | Cheap must-ship mesh-objectivity fix + Corotational, independently mergeable PR | Cantilever 1/2/4/8-elem objectivity (peak ≤2%, dissipated energy ≤3% via EnergyBalance, post-peak slope ≤5%); **negative control on stock `DispBeamColumn` must FAIL** the energy gate; per-IP `lch` white-box probe ✓; elastic large-rotation benchmarks (Bathe–Bolourchi cantilever, end-moment-rolls-to-circle, Lee's frame / von Mises arch) within published tolerances; 2D+3D build + serialize across SP/MP; ledgers/banner/ADR landed in the same PR |
| **2 — Tier-2 embedded hinge (2D)** | `lch`-free softening hinge composing with corotational; mesh- AND integration-objective | Constant-moment patch test to machine precision; reduce-to-Tier-1 when no hinge active; **pre-cracked finite-rotation invariance** (rotate an open-hinge element 90/180° → zero force/dissipation increment); single-element displacement-controlled `∫M d[[θ]] == Gf` to machine precision (solver-independent, run BEFORE the multi-mesh arc-length test); integration-objectivity sweep (split fixed-rule-vary-nIP and fixed-nIP-swap-rule) shows the hinge removes the residual nIP drift Tier-1 keeps; no energy double-count (section frozen at activation); state cycles correctly across commit/revert and SP/MP; combined softening+large-disp collapse test under [[31_ladruno_robust_solve_driver_adr]] |
| **3 — 3D hinge + DDM (deferred)** | Biaxial/torsional jump on `CorotCrdTransf3d`; optional sensitivity | own ADR; 2D gate suite lifted to 3D incl. a 3D finite-rotation invariance test; out of scope until prioritized |

### State & serialization

`α / αCommit` + per-IP **irreversible** localization flags are committed state (`commitState` ships `αCommit`; forgetting that flags are committed can resurrect/lose a hinge on revert and corrupt irreversibility). Use a **fixed layout sized to `numSections`** (an `α` vector + a flag `ID`), not a variable count-prefixed payload (the LadrunoBrick precedent). `K_αα` is rebuilt in `setDomain`, never serialized.

## Risks / open questions

> [!question]
> **Q1 — MPCO scope?** Is stock `.mpco`/STKO output a required deliverable, or is the (generic) Ladruno recorder sufficient? `MPCORecorder.cpp::getGeometryAndIntRuleByClassTag()` is a hard-coded tag switch — a new tag falls through to a point-cloud default with wrong per-IP coordinates, so MPCO support means **two new case entries + a `LEDGER_vanilla_files` row** (touches a vanilla file). The Ladruno recorder is generic and needs nothing.

> [!question]
> **Q2 — Tiers coexist or globally exclusive?** Intended design is coexist (fiber away from the hinge may still use Tier-1 `lch` for distributed softening; the hinge uses `Gf`), with the mandatory mutual-exclusion switch **only at the active hinge IP**. Confirm.

> [!question]
> **Q3 — Basic strain template.** Default to the `DispBeamColumnNL2d`-style `½θ²` bowing term (recommended; raises 2D core effort) over the cheaper linear-strain `DispBeamColumn2d` template?

> [!question]
> **Q4 — Cyclic validation oracle.** Which named PEER Structural Performance Database specimen (axial ratio, reinforcement, cyclic protocol) is the experimental capstone? Pick one consistent with the existing rc3d recipes.

> [!question]
> **Q5 — 3D timing.** Confirm v1 ships the 3D element as **Tier-1 only**, deferring the 3D biaxial/torsional embedded hinge to its own ADR (given the quaternion-triad finite-rotation complexity).

- **Solver dependency (largely mitigated).** Stage-2 mesh-objectivity / collapse gates need a snap-back-capable follower. **`LadrunoIndirectControl` (Integrator tag 33006) is already built** (`SRC/analysis/integrator/LadrunoIndirectControl.cpp`, CMOD/indirect control, monotone through snap-back) — it is the primary follower. `LadrunoArcLength -stabilize` and `LadrunoDynamicRelaxation` also exist (DR is a rest-state corroborator only, never a descending-branch tracer). The only genuinely unbuilt piece is the *dissipation* arc-length variant for multi-crack ([[22_ladruno_dissipation_arclength_adr]], RESERVED). Add a cheap solver-independent gate (single-element `∫M d[[θ]] == Gf`) FIRST so constitutive bugs are caught without the path-follower.
- **Tangent inconsistency under mid-step hinge activation / line search** — the corotational geometric stiffness is built from `pb`, so composition is correct only if the inner `α` update is tightly converged at every globalize call. Mitigation: converge-`α`-inside-`update()`, cache, and have both `getResistingForce` and `getTangentStiff` read the single cached converged state — an explicit tested contract.
- **PDelta + softening hinge** cannot be assumed a free rider for post-peak geometric response; either validate once or explicitly mark `PDelta`+softening unsupported in v1.
- **Mixed strain-reference subtlety** — axial uses deformed chord length while curvature/shear integrate over initial `L0` (corotated small-strain approximation); state the assumption + validity envelope.
- **Backwards compatibility:** new element, no change to existing model behavior. Stock `DispBeamColumn` is untouched (its mesh-dependence remains; the [[LEDGER_quirks]] entry stays as the documented reason to prefer this element).

## Implementation log

### 2026-06-16 — Stage 0 + Tier-1 (2D) landed on branch `guppi/ladruno-dispbeamcolumn`

- **Files:** `SRC/element/ladrunoDispBeamColumn/{LadrunoDispBeamColumn2d.{h,cpp}, OPS_LadrunoDispBeamColumn.cpp, CMakeLists.txt}`. Cloned from `DispBeamColumn2d` via symbol-rename (faithful clone), then Tier-1 edits applied.
- **Tier-1 implemented:** `current_section_lch` member set to `wt[i]*crdTransf->getInitialLength()` (reference length) inside the `update()` IP loop immediately before each `setTrialSectionDeformation`; `getCharacteristicLength()` override returning it; `ops_TheActiveElement = this` defensively at the top of `update()`; `-lch {ip|element|<value>}` flag (default `ip`), with `element` emitting the energy-double-count warning. `lchMode`/`userLch` serialized (data Vector 16→18); `current_section_lch` is transient (not serialized).
- **Registration (the gotcha):** a Ladruno element needs registration in **THREE** places, not two — `classTags.h` (`ELE_TAG_LadrunoDispBeamColumn2d=33013`), `FEM_ObjectBrokerAllClasses.cpp` (include + case), `OpenSeesElementCommands.cpp` `functionMap` (OpenSeesPy path), **and `SRC/element/TclElementCommands.cpp` `ladrunoElementTable` (the standalone Tcl `OpenSees.exe` path)**. Missing the last one builds & links clean but yields `element ... not known` only in the Tcl exe. Recorded in [[LEDGER_quirks]].
- **Build note:** editing `classTags.h` forces a wide recompile; a CMake reconfigure that adds a new `add_subdirectory` needs the build re-run once to actually compile the new TU (first run regenerates `build.ninja` but may not compile the new sources). Pre-existing MUMPS bootstrap in `build.bat` can spuriously re-enter on a stale `mumps-src/`.
- **Verified (Tcl `OpenSees.exe`):** (1) element creates + analyzes; (2) **reduce-to-stock is bit-identical** — elastic 2-element cantilever tip displacement equals vanilla `dispBeamColumn` to round-off (RELDIFF 0.0); (3) all `-lch` variants parse and run, `-lch element` warns as designed.
- **Tier-1 lch-delivery PROVEN (2026-06-16):** a single 2-point-Lobatto element with an `ASDConcrete1D -autoRegularization 100` tension-softening fiber section, pulled past peak under displacement control, gives the SAME peak (~280 N ≈ ft·A) but a DIFFERENT post-peak tail for `-lch ip` (lch = wt·L = 50 → residual 11.5 N) vs `-lch element` (lch = L = 100 = lch_ref → residual 0.3 N, the raw backbone). The 2× fracture-energy scaling (`lch_ref/lch`) is exactly the auto-regularization response, confirming the per-IP `lch` reaches the material and the flag controls it.
- **Mesh-objectivity VERIFIED (2026-06-17):** concrete fiber-section cantilever, `ASDConcrete1D -autoRegularization`, pushed past peak, dissipated energy = ∫V·dδ:
  `-lch ip` → 33.49 / 31.14 / 30.47 (N=2/4/8, N4→N8 = 2.1%); `-lch element` → 20.63 / 19.25 / 18.86 (N4→N8 = 2.0%). BOTH modes converge (regularization works, no mesh-dependence pathology — auto-reg cancels lch for objectivity), and they converge to DIFFERENT values (ip ≈ 1.6× element) because the per-IP band feeds the localizing Gauss–Lobatto IP its correct tributary length `wt·Lₑ` whereas `element`/stock over-estimates the band and under-dissipates — the §59 failure, quantified.
- **Regression test landed:** `tests/test_ladrunoDispBeamColumn2d_element.py` (12 tests, pass via OpenSeesPy — also the first exercise of the `functionMap`/OpenSeesPy path): reduce-to-stock, `-lch` accept/reject (incl. inf/nan/<=0), lch-delivery (ip > 1.2× element), mesh-objectivity convergence.
- **Honest boundary:** Tier-1 only helps lch-CONSUMING materials (ASDConcrete1D/ASDSteel1D/LadrunoUniaxialJ2+Lemaitre). A non-regularizing material (e.g. `Concrete02`) ignores `getCharacteristicLength` → stays mesh-dependent regardless of `-lch`; that general case is what Tier-2 (embedded hinge) addresses.
- **Stage-1 large-disp (Corotational) VERIFIED (2026-06-17):** `tests/test_ladrunoDispBeamColumn2d_element.py::test_corotational_large_displacement_matches_stock` — elastic cantilever driven into large deflection under a Corotational transform matches stock `dispBeamColumn` bit-identically (large-disp is in scope via the existing `CrdTransf`, no element-side geometric code). Shipped #255.
- **Stage-1 `½θ²` NL-strain toggle SHIPPED (2026-06-17):** `-nl` flag adds the `DispBeamColumnNL2d` bowing strain `ε₀ = v(0)/L + ½θ²` (θ from the Hermitian slope interpolation), with the matching geometric tangent in `getBasicStiff` (axial-force + B'ksC coupling) and the bowing term in the force/tangent `q`-loops; `getInitialStiff` stays linear; `nlGeom` serialized (data Vector 18→19). Default `nlGeom=0` keeps the linear basic strain (reduce-to-stock unchanged). Verified: reduces to linear at small deformation; under an axially-restrained cantilever at finite rotation the `-nl` bowing builds restraint tension that stiffens the member (|tip| drops >2% vs linear). 15/15 tests pass.
- **NOT yet (next):** 3D sibling (Tier-1 + `-nl`), then Stage-2 (embedded hinge).

*(move to `Ladruno_internal/implemented_<name>.md` when Stage 1 merges to `ladruno`)*
