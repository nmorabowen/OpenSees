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

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

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
- `MAT_TAG_LadrunoCohesiveHinge = 33003` — BUILT (Stage 2, 2026-06-17). Uniaxial band (used: 33000 J2, 33001 RebarBuckling, 33002 BondSlip; 33003 was the next free slot).

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
- **3D sibling SHIPPED (2026-06-17):** `LadrunoDispBeamColumn3d` (`ELE_TAG 33014`) — clone of `DispBeamColumn3d` + the Tier-1 per-IP `lch` channel + the `-nl` biaxial bowing `ε₀ = v(0)/L + ½(θz²+θy²)` (ported verbatim from `DispBeamColumnNL3d`: biaxial geometric tangent in `getBasicStiff`, bowing in the force/tangent `q`-loops). Reached through the same `LadrunoDispBeamColumn` command (ndm-dispatch: ndm2/ndf3 → 2D, ndm3/ndf6 → 3D). `lchMode`/`userLch`/`nlGeom` serialized (data 16→19). Verified (`tests/test_ladrunoDispBeamColumn3d_element.py`, 13 tests): reduce-to-stock bit-identical, `-lch`/`-nl` accept + inf/nan/≤0 reject, Corotational large-disp == stock, `-nl` reduces-to-linear / stiffens under finite rotation. 28/28 (2D+3D) pass.
- **Adversarial review of -nl + 3D (2026-06-17):** 10-agent workflow, verdict **merge-ready, 0 must-fix**; geometric tangent verified term-by-term against NL2d/NL3d. Two low-sev should-fixes applied: (1) 3D `getTangentStiff` applied the damping stiffness-multiplier only on the `-nl` branch (via `getBasicStiff`) — moved the multiplier to `getInitialStiff` (mirroring 2D) so both tangent branches are consistent under `-damp`; (2) corrected the stale dispatcher comment (3D is built) + documented `-nl`. Nice-to-have left: 3D `getTangentStiff` recomputes-then-discards the linear kb when `-nl` (perf only). 28/28 still green.
- **NOT yet (next):** Stage-2 (embedded hinge). Handoff in [[ladruno_handoff]] (Track 3).

### 2026-06-17 — Stage 2 started: `LadrunoCohesiveHinge` cohesive material (merged [#264](https://github.com/nmorabowen/OpenSees/pull/264))

The Tier-2 hinge's **discrete cohesive law** ships first, standalone, before any
element wiring — it is the independently-mergeable, lowest-risk piece and it
unblocks the *cheap solver-independent energy gate* (`∫M d[[θ]] == Gf`) the ADR
says to land FIRST (before the path-follower / multi-mesh tests).

- **Files:** `SRC/material/uniaxial/LadrunoCohesiveHinge.{h,cpp}` (`MAT_TAG 33003`),
  `tests/test_ladrunoCohesiveHinge_material.py` (10 tests).
- **Law (rigid–softening cohesive + irreversible secant-to-origin damage):**
  near-rigid penalty `Kpen` pre-peak (hinge CLOSED, `M = Kpen·κ`) up to `M = Mc`,
  then softening — **exponential** (default) `M = Mc·exp(−a(κ−κ0))` or **linear**
  `M = Mc(1 − (κ−κ0)/κf)` — calibrated so the **monotonic envelope integrates to
  EXACTLY `Gf`**: `κ0 = Mc/Kpen`, `Esoft = Gf − Mc²/(2Kpen)`, `a = Mc/Esoft`,
  `κf = 2Esoft/Mc`. Strain = rotation jump `[[θ]]`, stress = cohesive moment `M`.
  **No `lch` to calibrate** — `Gf` is consumed directly (the whole point of the
  discrete hinge vs. Tier-1 smearing).
- **Guarded penalty floor (ADR "near-rigid penalty with a guarded floor"):** both
  calibrations require `Kpen > Mc²/(2Gf)` (else pre-peak energy alone exceeds `Gf`).
  Default `Kpen = penaltyRatio·floor` (`penaltyRatio` default 1000 → pre-peak energy
  ≈0.1% of `Gf`, near-rigid). A user `-penalty K` below the floor is **clamped with
  a loud `opserr` warning**, not accepted.
- **Irreversibility:** `κmax = max|κ|` is the monotone damage driver; unload/reload
  is the secant `Ksec = Menv(κmax)/κmax` to the origin (isotropic damage) — no
  strength recovery, no residual moment at zero jump.
- **`getEnergy()`** overridden: trapezoidal path work `∫M d[[θ]]` (exact for the
  piecewise-linear LINEAR envelope).
- **Parser:** `uniaxialMaterial LadrunoCohesiveHinge tag Mc Gf <-penalty K | -penaltyRatio r> <-exp|-linear>`;
  rejects `Mc≤0`/`Gf≤0`. Registered at the 5 uniaxial sites (`classTags.h`,
  `FEM_ObjectBrokerAllClasses.cpp`, `OpenSeesUniaxialMaterialCommands.cpp` (Py),
  `TclModelBuilderUniaxialMaterialCommand.cpp` (Tcl), `CMakeLists.txt`).
- **Verified (OpenSeesPy, 10/10):** energy gate `∫M d[[θ]] == Gf` — LINEAR exact to
  `1e-9` (finite-jump break), EXP converges to `2e-3` (tail truncation); `peak == Mc`
  & near-rigid activation at `κ0`; initial tangent `== Kpen`; irreversible secant
  unload (returns to ~0 moment) + no reload recovery; sign symmetry; guarded floor
  clamps & still dissipates `Gf`; DB sendSelf/recvSelf round-trips the committed
  `κmax/work`. Existing 28/28 DispBeamColumn tests unchanged (only `classTags.h` +
  uniaxial plumbing touched).
- **NEXT:** wire the cohesive material into the element as the Tier-2 enhanced
  rotation-jump `κ = B·d + G·α` with `α` statically condensed to the 3-DOF basic
  system BEFORE `crdTransf` (the PINNED INVARIANT) — the genuinely hard EAS piece.

### 2026-06-17 — Adversarial pre-implementation review of Stage-2 element wiring → CORRECTIONS (4-agent review)

Before coding the element-side hinge, a 4-agent adversarial review (corotational
composition / EAS-condensation+activation / energy-double-count+freezing /
kinematics+serialization+sequencing) was run against this ADR's Tier-2 plan. It
found **several load-bearing errors in the plan as written**. The plan is amended
as follows; the original §Tier-2 / landmines text is superseded where they conflict.

1. **Kinematics: the shipped base is EULER–BERNOULLI, not Timoshenko.** Stage-1
   shipped a faithful clone of `DispBeamColumn2d` (`LadrunoDispBeamColumn2d.cpp:491`
   == `DispBeamColumn2d.cpp:671`, cubic-Hermite linear curvature; the §"Kinematics
   (default Timoshenko)" text and the `TimoshenkoBeamColumn2d.cpp:610` reference are
   WRONG and were never built). The hinge `G` operator is therefore built against
   the E–B curvature field `κ_d = (1/L)[(6ξ−4)θ_i + (6ξ−2)θ_j]`. Do NOT pull in the
   Timoshenko `phi`/shear-row B-operator — that would break reduce-to-Tier-1.

2. **"Freeze the fiber section" is unimplementable AND unnecessary — use the real
   Armero–Ehrlich split.** There is no `SectionForceDeformation` API to freeze/elastic-
   unload a section; `setTrialSectionDeformation` always re-evaluates fibers from the
   strain handed in. The correct mechanism: the **bulk section sees the BOUNDED
   enhanced curvature** `κ_bulk = B·v + Ḡ·α`, where `Ḡ = −1/L` (constant, bounded —
   the smooth part of the incompatible mode). As `α` grows, `Ḡ·α` *unloads* the bulk
   (its curvature decreases), so the bulk dissipates nothing further and unloads on
   its own constitutive law — no freeze flag, no double count. ALL post-peak
   dissipation is carried by the cohesive `M([[θ]])` on `α`. Delete the contradictory
   "softening branch of the same section / not a parallel spring" language: post-
   activation this is operationally a discrete hinge in **series** with an elastically-
   unloading bulk (Jukić–Brank–Ibrahimbegović 2013 embedded-discontinuity beam).

3. **The Dirac never enters quadrature; orthogonality is machine-exact.** The jump
   (Dirac) part is handled DIRECTLY by the cohesive law (`M_coh(α)` added to the
   enhancement residual), not integrated. Only the bounded `Ḡ = −1/L` is integrated:
   `Σ wt_k Ḡ = −(1/L)·Σ wt_k = −1` exactly (Lobatto/Legendre weights sum to L), and
   the jump contributes `+1`, so `∫G dx = 0` holds to machine precision under ANY
   rule — resolving the "Dirac-at-one-IP vs ∫G=0 is unsatisfiable" objection.

4. **The condensation algebra (per element, scalar α, ONE hinge):**
   - enhancement residual `h(α) = Σ_k wt_k · Ḡ · M_sec,k(κ_bulk) + M_coh(α) = 0`
   - `K_αα = Σ_k wt_k · Ḡ · EI_sec,k · Ḡ + dM_coh/dα`  (bulk term ≥0; cohesive term <0
     post-peak) — **`K_αα` is sign-discontinuous at activation and INDEFINITE across
     the whole softening branch, not merely zero "at the peak"** (the landmine-#2
     premise was mis-stated). The condensed tangent `K_basic = K_vv − K_vα K_αα⁻¹ K_αv`
     contains `K_αα⁻¹` regardless of how α is found — **"closed-form jump update" does
     NOT escape the singularity.** Mitigation that ACTUALLY works: a **guarded
     reciprocal** `1/K_αα` with a magnitude floor (`|K_αα| ≥ ε·K_αα0`), the bulk
     `Ḡᵀ EI Ḡ` providing a positive stabiliser. Inner Newton on the scalar α with the
     same floor.
   - `K_vα = Σ_k wt_k · B_kᵀ · EI_sec,k · Ḡ` (3-vector), `K_αv = K_vαᵀ`. Axial DOF
     untouched in v1 (linear basic strain only).
   - **The basic FORCE needs NO explicit condensation correction**: at converged α
     (`h=0`) the sections already hold `κ_bulk`, so `q = Σ wt_k B_kᵀ M_sec,k` IS the
     condensed basic force. Only the TANGENT gets the `−K_vα K_αα⁻¹ K_αv` correction.

5. **`-nl` bowing × α couples the jump into the axial/geometric channel** (the `½θ²`
   membrane strain makes the axial force depend on the bending rotations, which the
   hinge curvature perturbs). v1 **forbids `-nl` + `-hinge` together** (parser error);
   the hinge uses the linear basic strain. (Revisit the cross-term algebra later.)

6. **Reduce-to-Tier-1 must be byte-identical → GATE the α machinery** on a per-element
   `hingeOn` flag, exactly like the `-nl` flag. With no `-hinge` option the force/
   stiffness/update paths take the IDENTICAL code path (and FP summation order) they
   take today — never "condensation with zero operands" (which perturbs at O(ulp) and
   fails the `RELDIFF==0.0` gate).

7. **State/serialization & cached-tangent (the non-Newton hole):** `getTangentStiff()`
   can be called WITHOUT a preceding `update()` under ModifiedNewton / KrylovNewton /
   line-search / arc-length — so a naive "converge α in update(), cache, readers
   trust cache" contract reads a STALE α under exactly the solvers Stage-2's later
   gates use. **PR-2a restricts to monotonic full-Newton / DisplacementControl** (the
   solver-independent gates need nothing else); the idempotent re-converge-α-from-
   commit hardening is deferred to PR-2b. The localization/`αCommit` flag flips ONLY
   in `commitState`; `revertToLastCommit` restores α to `αCommit` and rebuilds the
   hinge tangent (else a rejected step resurrects/loses a hinge). Serialize the hinge
   material like the sections (classTag+dbTag+sendSelf) + scalar `αCommit`/`hingeOn`/
   `hingeLoc`; gate so old (no-hinge) streams are unaffected.

8. **`setDomain` "non-softening fibers at hinge IP" check is unbuildable** — there is
   NO `isSoftening()` query anywhere in `SRC/`. Downgrade to a **documented user
   contract** for v1 (the hinge carries `Gf`; pairing it with a softening fiber
   section at the same IP double-counts and is the user's responsibility). A real
   `isSoftening()` interface is a separate vanilla-file change, deferred.

9. **Energy gate must measure ELEMENT total dissipation** (`∫F·v` via EnergyBalance)
   on an OTHERWISE-ELASTIC section `== Gf` — the cohesive material's own `Twork==Gf`
   gate (already green) cannot detect the section re-integrating the same energy.

#### PR-2a scope (the minimal correct first increment — what ships next)

- `-hinge <matTag> [loc]` (or `-hinge -Mc Mc -Gf Gf ...` auto-building a
  `LadrunoCohesiveHinge`): ONE rotation-jump `α` (scalar) condensed per element with
  the guarded reciprocal; default hinge location mid-element. Linear basic strain
  only (mutually exclusive with `-nl`). Linear + Corotational transf, monotonic.
- Everything GATED on `hingeOn` so the no-hinge path is byte-identical to Stage-1.
- **Gates (solver-independent, no path-follower):** (a) reduce-to-Tier-1 `RELDIFF==0.0`
  with `-hinge` absent; (b) constant-moment patch test to machine precision (a single
  element under pure moment: the hinge carries exactly the applied moment, zero
  spurious stress); (c) single-element `∫M d[[θ]] == Gf` under DisplacementControl;
  (d) FD-tangent of the condensed `K_basic` (incl. an open-hinge state); (e) element
  total-dissipation `== Gf` on an elastic section; (f) commit/revert state cycle.
- **Deferred to PR-2b+:** integration-objectivity sweep, pre-cracked finite-rotation
  invariance, `-nl`+hinge cross-terms, non-Newton/line-search α re-convergence, the
  `ladruno_drive` collapse test (needs the still-RESERVED dissipation arc-length), 3D.

### 2026-06-17 — PR-2a SHIPPED: element-side embedded hinge (`-hinge $matTag`, 2D)

The corrected formulation above, implemented and validated. Reuses `ELE_TAG 33013`
(no new tag — same element, new gated option).

- **Files:** `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn2d.{h,cpp}` +
  `OPS_LadrunoDispBeamColumn.cpp` (usage), `tests/test_ladrunoDispBeamColumn2d_hinge.py` (8 tests).
- **What landed:** `-hinge $matTag` adds a single scalar rotation jump `α` carried by any
  `UniaxialMaterial` (default `LadrunoCohesiveHinge`). `update()` runs an inner Newton on
  `α` (sections see `κ_bulk = B·v − α/L`; residual `h = −Σwt·M_sec + M_coh(α)`; **guarded**
  `1/K_αα` with a magnitude floor against the bulk term). `getTangentStiff()` subtracts the
  cached condensation `K_vα K_αα⁻¹ K_αv` from the 3×3 basic stiffness **before** `crdTransf`
  (pinned invariant); `getResistingForce()` is UNCHANGED (sections hold the converged
  `κ_bulk`, so `q` is already the condensed basic force). Inner-Newton convergence uses a
  running moment scale `hingeMscale` (~Mc) so the tol does not collapse when a fully-broken
  LINEAR hinge carries `M→0`. ALL gated on `hingeOn`; `-hinge`+`-nl` rejected at parse.
  commit/revert/revertToStart + sendSelf/recvSelf (separate hinge-material send like the
  sections, + `hingeOn`/`hingeJumpCommit`, data Vector 19→21).
- **Verified (OpenSeesPy, 8/8; full 46/46 element+material+hinge):** reduce-to-Tier-1
  **bit-identical** with `-hinge` absent; constant-moment **patch test** — the hinge carries
  exactly the applied moment to `1e-9` (enhancement equilibrium `h=0`); **energy gate** — a
  cantilever past the hinge peak dissipates `Gf` (LINEAR `19.999992` vs `20`, ~4e-7 rel) and
  peaks at `Mc`; **element total dissipation == Gf** on an elastic section (closed
  load/unload cycle, no bulk double-count); **tangent consistency** — every step through
  EXP softening converges under NormDispIncr-1e-12 Newton ≤10 iters; commit/revert + DB
  sendSelf/recvSelf round-trip of an open hinge.
- **Honest scope:** `setDomain` non-softening-fiber check is a documented user contract (no
  `isSoftening()` exists). FD-tangent is via tight-Newton convergence (the testbed has no
  element-DOF FD harness). PR-2b carries the deferred list above.

### 2026-06-17 — PR-2b: 2D Stage-2 objectivity / invariance / robustness gates (test-only)

Empirical probing showed the PR-2a hinge **already passes** the deferred Stage-2 validation
gates — the implementation is robust; PR-2b locks them in as regressions (no source change).

- **Probed first (not assumed):** under a transverse-load cantilever with the **Corotational**
  transform driven to a **74° tip rotation**, the hinge dissipates `Gf = 20.0000` exactly with
  zero step failures → the pinned invariant (condensed basic K/q through `crdTransf`) composes
  correctly with the rotating frame. And under **Newton / ModifiedNewton / NewtonLineSearch /
  KrylovNewton** the dissipation is identically `Gf` (the residual is always evaluated
  post-`update()`, so tangent reuse / line search never corrupt the converged dissipation) —
  the non-Newton "stale-α" hole the review flagged does **not** bite in practice, so the
  deferred idempotent-re-converge hardening is **unnecessary**.
- **Gates added** (`tests/test_ladrunoDispBeamColumn2d_hinge.py`, 8→15): corotational
  large-rotation `Gf`-dissipation (tip rot > 0.5 rad); orientation invariance (member at 0° vs
  90° under Corotational → identical M–θ path + `Gf` to round-off); integration-objectivity
  sweep (Lobatto nIP ∈ {2..6} → invariant M–θ + `Gf` to 1e-6, the discrete hinge has NO
  residual nIP drift); solver-robustness (4 algorithms all dissipate `Gf`). Full 53/53 green.
- **2D Stage-2 is now gate-complete** except the `ladruno_drive` snap-back collapse test (still
  blocked on the RESERVED dissipation arc-length, [[22_ladruno_dissipation_arclength_adr]]) and
  the `-nl`+hinge cross-terms. **3D** is the next real frontier (own ADR — quaternion-triad
  finite-rotation biaxial/torsional jump).

*(move to `Ladruno_internal/implemented_<name>.md` when Stage 1 merges to `ladruno`)*
