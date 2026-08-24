---
title: Selective mass scaling for explicit dynamics (CentralDifferenceSMS)
project: Ladruno
status: ready-to-implement
priority: high
owner: nmora
tags:
  - implementation
  - adr
  - integrator
  - explicit
  - mass
  - mass-scaling
---

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

# Selective mass scaling for explicit dynamics (`CentralDifferenceSMS`)

> **Design / ADR (pre-implementation).** Implements roadmap
> [[Ladruno_explicit_roadmap]] §5.1 (the future plan it points at,
> `[[01_selective_mass_scaling]]`, is superseded by this ADR). A new
> `TransientIntegrator` subclass — **classTag `INTEGRATOR_TAGS_CentralDifferenceSMS =
> 33007`** (next free integrator slot) — that adds fictitious nodal mass to the
> elements throttling the timestep, so a fine local mesh can run at the global
> timestep set by the *bulk* of the model. LS-DYNA `*CONTROL_TIMESTEP DT2MS`
> semantics. **Recommended-after [[35_ladruno_hrz_lumped_mass_adr]]** — SMS *functions* on
> the shipped `Diagonal` lump (per-node shares already positive), but HRZ makes "the
> element's mass" well-defined and conserving across the full element zoo, so 35 should
> land first. Lumped scaling ships first; consistent (Olovsson) scaling is parked as v2.

## What

`CentralDifferenceSMS` is a thin **subclass of `CentralDifferenceLadruno`**: it inherits
the correct first-step starter, the leap-frog state, the dual-velocity output, the βK
guard and the `dt_cr` reporting verbatim, and adds exactly one thing — a
**mass-scaling pass at `domainChanged()`**. The pass:

1. computes each element's stable step `Δtₑ` from the per-element `K v = λ M v`
   eigensolve already factored out of `CriticalTimeStep`;
2. for every element with `Δtₑ < Δt_target`, computes the mass scale
   `sₑ = (Δt_target / Δtₑ)²` and the increment `Δmₑ = (sₑ − 1)·mₑ`;
3. **distributes `Δmₑ` to the element's nodes** and adds it to those nodes' mass
   (additive `Node::setMass`, re-baselined on every `domainChanged`);
4. reports total added mass as a fraction of model mass, warns past a cap.

The scaling logic lives in a shared header-only `Ladruno::buildMassScaling(...)`
util so a future `ExplicitBatheSMS` is a 20-line clone. `Δt_target` is user-supplied
(direct, like `DT2MS`), with an optional `-targetFactor f` mode that derives it from a
percentile of the element-`Δtₑ` distribution so the bulk governs.

**In scope:** lumped (diagonal) mass scaling on `CentralDifferenceSMS`; optional
restriction to an element/node set; added-mass diagnostics; Tcl + Python command. **Out
of scope:** consistent/Olovsson anisotropic scaling (v2, `> [!question]` below);
`ExplicitBatheSMS` (trivial follow-up once the util is proven); wave-speed `ℓ/c` Δtₑ
(we reuse the eigensolve — see M3); automatic modal-mass validation (diagnostic hook
only).

## Why

Verbatim from roadmap §5.1: *"without this, the smallest element in the mesh bounds the
global timestep — and in any realistic 3D SSI / pile-soil / contact-zone model there
will be a few tiny elements (interface zones, refinement near loads) that throttle the
entire run … **This is the gating prerequisite for any practical 3D explicit SSI work
in OpenSees.**"* Today OpenSees has **no** mass-scaling facility of any kind
(`LadrunoFictitiousMass` is for quasi-static DR/arc-length, not timestep
acceleration). LS-DYNA has had `DT2MS` since forever; this closes the single biggest
practical-scalability gap on the explicit axis.

## Where

- New code:
  - `OpenSees/SRC/analysis/integrator/CentralDifferenceSMS.{h,cpp}` — subclass of
    `CentralDifferenceLadruno`.
  - `OpenSees/SRC/analysis/integrator/LadrunoMassScaling.h` — header-only
    `namespace Ladruno`. **Note:** `LadrunoFictitiousMass.h` is a lifecycle/idiom
    precedent only — it iterates `getFEs()`/`getID()` (equation numbers, for SOE
    injection), which is the **wrong granularity** here. For per-node `setMass` the
    correct mapping is `Element::getNodePtrs()` (as used in `CriticalTimeStep` /
    `EnergyBalanceKernel`).
- Modify:
  - `OpenSees/SRC/classTags.h` — `#define INTEGRATOR_TAGS_CentralDifferenceSMS 33007`
    (integrator band is 33000–33006 today; 33007 is next-free). Add the inline
    "per-registry, no collision" annotation: `33007` is *also* `ELE_TAG_LadrunoQuad` in
    the **element** registry, but the G2 classtag gate keys collisions within-family
    (`ci/check_classtags.py`), so this is safe and already precedented (e.g. 33004 =
    `LadrunoArcLength` + `LadrunoIMKBeam2d`).
  - `OpenSees/SRC/analysis/integrator/CriticalTimeStep.{h,cpp}` — **small refactor**:
    extract the per-element `(Mₑ,Kₑ) → ωₑ,max` kernel (currently inline in the loop at
    lines 222-282) into a reusable helper so both `computeCriticalTimeStep` and the
    mass-scaling pass call it (avoids a second copy of the D8-safe fetch+lump+eigensolve).
  - `OpenSees/SRC/domain/node/Node` — no source change; the pass uses the existing
    `getMass()`/`setMass()` (read-modify-write, since `setMass` overwrites — `Node.cpp:1288`).
  - `OpenSees/SRC/api/...` + `OpenSees/SRC/tcl/...` integrator command parsing
    (wherever `CentralDifferenceLadruno`/`ExplicitBathe` are dispatched) — register the
    new keyword.
  - `OpenSees/SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` +
    `TclPackageClassBroker` — register classTag 33007 (the D7 pattern from ADR 04).
- Reference: `LadrunoDynamicRelaxation` (`buildFictitiousMass()` called from
  `domainChanged()` at line 292; `addDiagonalToSOE` from `formTangent()` at line 179) —
  the exact lifecycle slot and injection idiom to mirror, except we inject into **node
  mass** not the SOE (see M2).
- Build: header-only util; the new `.cpp/.h` join the explicit-integrator build group.
  Stamp headers + `stamp_headers.py` GLOB. Banner: add a `banner_features.txt` line and
  run `patch_banner.py` (per CLAUDE.md banner workflow). Ledger:
  `LEDGER_implementations.md` row (classTag 33007, integrator, files, status) — the
  class-tag allocation is reserved *here* and added to `classTags.h` only when code
  merges.

## How

### Formulation

For an element with lumped mass `mₑ` and stiffness `kₑ`, the undamped central-difference
stable step scales as

- ωₑ,max² = λₑ,max  (largest generalized eigenvalue of `Kₑ vₑ = λ Mₑ vₑ`),
- Δtₑ = 2 / ωₑ,max.

Scaling **all** of the element's mass by a factor `s` leaves `Kₑ` unchanged and divides
the eigenvalue by `s`, so ωₑ,max → ωₑ,max/√s and Δtₑ → Δtₑ·√s. To raise Δtₑ to a target
Δt_target (only when Δtₑ < Δt_target):

```
sₑ  = (Δt_target / Δtₑ)²          ≥ 1
Δmₑ = (sₑ − 1) · mₑ               (the added "fictitious" mass for element e)
```

**betaK-damped sizing (v1.1, 2026-06-20).** Stiffness-proportional (betaK) Rayleigh
damping shrinks the explicit stable step (ξ = betaK·ω/2 grows with ω), so the undamped
`sₑ = (Δt_target/Δtₑ)²` UNDER-scales and the element is still unstable at Δt_target. With
`c = betaK/Δtₑ` (= ½·betaK·ωₑ,max) the betaK-damped step at scale `s` is
`Δt_d(s) = (2/ωₑ,max)(√(s + c²) − c)`, which inverts in **closed form** to
`sₑ = T² + 2·T·c`, `T = Δt_target/Δtₑ`. This reduces to the undamped `T²` when betaK = 0
(no-damping models byte-identical) and injects more mass when betaK > 0 so the *damped*
step reaches Δt_target. Mass-proportional (alphaM) damping is intentionally **excluded** —
it does not reduce the high-frequency step that governs explicit stability, and folding it
in across scales is non-monotonic. The skip test also uses the damped step (an element
whose undamped Δtₑ exceeds Δt_target but whose *damped* step does not is now correctly
scaled). Tested: `test_massScaling_validation.py::test_betaK_damped_sizing`.

`Δmₑ` is distributed to the element's nodes in proportion to their existing lumped share
(HRZ lump from [[35_ladruno_hrz_lumped_mass_adr]]), then **added** to each node's mass
matrix diagonal. The scaling factor `sₑ` is a *per-element estimate* — it assumes all of
element e's mass is multiplied by `sₑ`, whereas the pass injects into shared **nodal**
mass, which the per-element pencil (`ele->getMass()`) cannot see. The realized global
stable step is therefore **conservatively bounded** (added diagonal mass only lowers
eigenvalues, so the step never drops below the un-scaled value and the bound runs in the
design's favor) but is **not** exactly Δt_target "by construction." The achieved step is an
*optimistic estimate* reported for diagnostics and **validated empirically** by tests 1–2,
not a guarantee. The fidelity cost is the frequency shift the inertia induces, monitored by
`ΔM / M_model` (cheap) and, optionally, a modal check (advanced).

Total added mass `ΔM = Σₑ Δmₑ` is tracked and capped (see `-maxAddedMass`).

### API

```
# Tcl / Python
integrator CentralDifferenceSMS $dtTarget <-maxAddedMass $frac> \
           <-set $nodeOrEleSetTag> <-targetFactor $f> <-verbose> \
           <-lump rowsum|diagonal|hrz> <-tangent>
# NOTE: -cfl / -cflAbort / -recompute are rejected with SMS (see warning below)
```

- `$dtTarget` — target stable step (required in direct mode; `DT2MS`-style).
- `-targetFactor $f` — alt: derive `dtTarget` from the f-percentile of element Δtₑ
  (e.g. `f=0.1` → bulk governs, smallest 10% get scaled).
- `-maxAddedMass $frac` — hard cap; abort (or warn-and-clamp) if ΔM/M_model exceeds it.
  Default warn at 0.05 (roadmap §5.1 diagnostic).
- `-set` — restrict scaling to a node/element set (LS-DYNA part-restricted scaling).
- inherited `-lump / -tangent` flags from `CentralDifferenceLadruno`.

> [!warning]
> **`-cflAbort` and `-recompute` are a hard parse error when combined with SMS.** The
> inherited `newStep()` re-runs `computeCriticalTimeStep()`, which reads `ele->getMass()`
> (the *un-augmented* element mass SMS never touches) and, under `-cflAbort`, hard-aborts
> (`return -2`) precisely when `Δt > dt_cr` — i.e. exactly the configuration SMS creates by
> design. (`-cflAbort` also auto-enables the step-1 check, so it fires even without
> `-recompute`.) Either reject these flags at parse time, or feed the eigensolve a
> nodal-augmented mass source so it sees the scaled masses. v1 rejects them; the augmented
> eigensolve is a v2 extension point.

### Internals / data flow

```cpp
// CentralDifferenceSMS::domainChanged()  (after base-class state alloc, before first step)
int CentralDifferenceSMS::domainChanged() {
  if (CentralDifferenceLadruno::domainChanged() < 0) return -1;
  removePriorScaling();                       // re-baseline: undo last pass's ΔM
  Ladruno::MassScalingReport rep =
      Ladruno::buildMassScaling(theModel, dtTarget, lumping, setRestriction,
                                /*out*/ addedNodalMass);
  applyScaling(addedNodalMass);               // node->setMass(getMass() + Δ) per node
  reportAddedMass(rep);                        // total, %, governing elements; warn>cap
  return 0;
}
```

`buildMassScaling` reuses the refactored `CriticalTimeStep` per-element kernel for Δtₑ.
`removePriorScaling`/`applyScaling` keep an integrator-owned record of the ΔM injected
per node (a `std::map<nodeTag, Vector>` or parallel arrays) so re-entrant
`domainChanged` calls re-baseline rather than compound, and the destructor restores the
original node masses (no permanent mutation of shared `Domain` state).

### Testing

Zone-A pytest (`tests/test_centralDifferenceSMS_integrator.py`):

1. **Single-DOF sanity** — two-mass chain, one tiny stiff element; assert ΔM applied
   only to the offending element's nodes, and post-scaling `criticalTimeStep()` ≥
   dtTarget.
2. **Stable-run gate** — a mesh with one element 100× smaller than the bulk: without
   SMS the run needs ~100× more steps / diverges at the bulk Δt; with
   `CentralDifferenceSMS dtTarget=bulk_dt` it runs stably and the response matches the
   small-Δt unscaled reference within tolerance (the central claim of mass scaling).
3. **Added-mass accounting** — assert reported ΔM/M_model matches the analytic Σ Δmₑ;
   `-maxAddedMass` aborts when exceeded.
4. **Frequency-shift bound (on an SSI-representative model, not a cantilever)** — the
   v2-escalation decision (M4) is about SSI rocking/soil modes, so run modal analysis
   before/after scaling on a soil-column + structure model; assert Δf₁ of the modes that
   matter < the documented target (roadmap: <1% at acceptable scaling levels) for a modest
   dtTarget, and *demonstrate* the shift growing as dtTarget is pushed (the honest-tradeoff
   test that feeds the lumped-vs-Olovsson decision).
5. **Re-baseline** — call `domainChanged` twice (e.g. after a mesh change); assert ΔM
   does not compound and the destructor restores original node masses.
6. **Reduce-to-base** — `dtTarget ≤ min Δtₑ` ⇒ no element scaled ⇒ results bit-identical
   to `CentralDifferenceLadruno`.

Reference reading (carry into the test/validation doc): LS-DYNA Theory §22 (solid
timestep), Vol I `*CONTROL_TIMESTEP` (`DT2MS`/`IMSCL`/`MS1ST`), Olovsson & Simonsson
2006 (selective scaling theory).

## Decisions

| # | Decision | Rationale | Consequence / extension point |
|---|----------|-----------|-------------------------------|
| M1 | **Subclass `CentralDifferenceLadruno`**, add only a `domainChanged` scaling pass | Inherits the correct starter, leap-frog, dual-velocity, βK guard, dt_cr — all the robustness already reviewed in [[05_robust_central_difference]] | Promote a few base members to `protected` + a `virtual` hook; ExplicitBatheSMS later reuses the util, not the class |
| M2 | Inject via **additive nodal mass** (`Node::setMass` read-modify-write), not SOE-level | Energy still *closes* because the `EnergyBalanceRecorder` KE sums **both** element (`ele->getMass()`) and nodal (`node->getMass()`) mass into `g_ke`, and the leap-frog `M` does the same via element+nodal `addMtoTang` — so injecting into the nodal term keeps KE and the inverted `M` consistent. Matches roadmap §5.1 step 3 and LS-DYNA semantics. SOE injection (the `LadrunoFictitiousMass` route) only changes the LHS and would desync the recorder's KE | Must re-baseline on `domainChanged` and restore on destruct (M5); see parallel risk |
| M3 | Δtₑ from the per-element eigensolve, not wave-speed `ℓ/c` | The eigensolve (`CriticalTimeStep`, D8-safe, DSYGVX) is element/material-agnostic and handles beams/shells/anisotropy `ℓ/c` can't | **Requires a refactor** — `maxGeneralizedEigenvalue` is file-static and `computeCriticalTimeStep` returns only the global min, so the per-element `(Mₑ,Kₑ)→ωₑ` kernel must be extracted to a reusable entry point (it is *not* "already factored out"); `ℓ/c` stays an optional cheap-estimate extension point |
| M4 | **Lumped scaling v1; consistent (Olovsson) v2** | Roadmap §6 Q2 decision: lumped is `DT2MS`-default and keeps `M` diagonal (the whole point of explicit). Consistent allows ~10× more aggressive scaling but needs per-element block solves — research-grade | Revisit if M4-validation (test 4) shows lumped shifts f₁ > 1% at needed scaling |
| M5 | **Re-baseline + restore** node masses (integrator owns the ΔM record), with an explicit **parallel contract** | `domainChanged` is re-entrant and node mass is shared `Domain` state; *and* under OpenSeesMP a partition-boundary node receives Δmₑ only from rank-local elements while the only MPI reduction is the scalar dt — so shared-node masses would desync across ranks and the global `M⁻¹` would be rank-inconsistent | Restore eagerly at the **start** of every `domainChanged` (don't rely on residency at teardown); **v1 = sequential or partition-interior nodes only, error on a boundary node**; v2 = `MPI_Allreduce` the per-shared-node injected ΔM. `std::map<nodeTag,Vector>` keyed on a stable identity |
| M6 | New standalone class + **classTag 33007**, shared util header | User decision (vs. flag-on-existing); isolates the mutation-of-node-state behavior to an opt-in integrator | Broker registration (D7 pattern); banner + ledger rows part of the merging PR |

## Risks / open questions

> [!question]
> **Lumped vs consistent (Olovsson) scaling.** Lumped shifts global frequencies; for
> seismic SSI the fundamental period is non-negotiable. Ship lumped (M4) and gate on
> test 4. If f₁ drift exceeds 1% at the scaling needed for real 3D SSI meshes, elevate
> consistent scaling to a v2 ADR.

> [!done] IMPLEMENTED v1.1 (2026-06-20)
> **Constrained nodes (equalDOF / rigidDiaphragm / rigidLink / generic MP_Constraint).** Adding mass to a
> constrained (slave) node interacts with the constraint handler and with the
> dt_cr-ignores-constraints caveat inherited from ADR 04/05. The shipped v1 only *warned*;
> **v1.1 EXCLUDES** any sub-target element touching an MP-constrained **slave** node
> (`Domain::getMPs()->getNodeConstrained()`): the injected mass would be redistributed
> through the constraint's `Tᵀ M T` and the dt boost would not land. Such elements are
> skipped, counted (`nConstrained`), and **reported as still governing** at their un-scaled
> `dt_e` (so the user knows dtTarget was not delivered for them — lower dt or remove the
> constraint — rather than silently under-delivering / mis-distributing mass). **SP/fixed
> nodes are deliberately NOT excluded** (mass on a removed DOF is inert, and the motivating
> SSI/pile case puts refinement right on fixed supports). The guard is keyed on MP *slave*
> nodes only. Tested: `tests/test_massScaling_validation.py::test_constrained_element_excluded_free_scaled`
> (constrained truss excluded + slave stays massless, free truss still scaled `1/2`). The
> deeper constraint-aware scaling (project ΔM through the constraint) is still future work,
> tied to the explicit-constraint-projection handler ([[project_explicit_constraint_projection]]).

> [!question]
> **dt_cr recompute is blind to nodal augmentation.** The per-element pencil uses element
> `Mₑ`, not nodal mass, so any post-scaling recompute reports the *un-improved* element
> Δtₑ — which is why `-cflAbort`/`-recompute` are rejected with SMS (see API warning;
> they would hard-abort a run that is in fact stable). Resolution: report the achieved
> stable step as an **optimistic analytic estimate** from the scale factors (not a
> guarantee — the realized global step is only conservatively bounded, see §Formulation);
> an honest post-hoc check requires a global `M⁻¹K` eigensolve (expensive, not done in
> v1).

- **Energy balance**: `EnergyBalanceKernel` already sums **both** element and nodal mass
  into KE (verified: `ele->getMass()` and `node->getMass()` both feed `g_ke`), so closure
  holds under nodal scaling with no recorder change — just **assert closure in test 2**.
  (The earlier "recorder bug to fix" contingency is moot.) Watch the half-step vs
  full-step velocity convention separately when interpreting KE.
- **Numerical**: guard `Δtₑ > 0` (massless/constraint elements → skip, not divide-by-zero);
  clamp pathological `sₑ`.
- **Lump self-consistency (SF-2)**: the Δtₑ estimate (eigensolve lump) and the Δmₑ
  distribution (HRZ lump) **must use the same lump** — otherwise they disagree (e.g.
  Diagonal ρAL/3 vs HRZ ρAL/2 ⇒ ~√(3/2)≈1.22× mis-sizing). Harmless in the default
  lumped-mass workflow (diagonal≡rowsum there), real under `-cMass 1`. v1: mandate HRZ for
  both, or document the estimate is non-conservative under `-lump diagonal`.
- **Backwards compat**: new integrator, no impact on existing models; `dtTarget ≤ min Δtₑ`
  is a no-op (test 6).
- **ABI/build**: header-only util, standard `.cpp/.h`; classTag/broker registration is
  the only cross-cutting touch (follow D7 exactly to avoid a parallel-stream break).
- **Ledger/banner debt**: `LEDGER_implementations.md` row + `classTags.h` define +
  `banner_features.txt` line + `patch_banner.py` rebuild are part of the merging PR
  (CLAUDE.md REQUIRED). Manifest/G9 gate will fail without the ledger row.

## Implementation log

- **2026-06-19 — adversarial review (45-agent workflow, 6 dimensions × refute-by-default verify).**
  Verdict: **sound-with-fixes**. Confirmed-correct (do not touch): the eigenvalue scaling
  math (`s` divides λ ⇒ Δtₑ→Δtₑ·√s; added diagonal mass keeps M SPD and only lowers
  eigenvalues); energy closure (recorder sums element+nodal KE); classTag 33007 free and
  G2-safe; the subclass + re-baseline concept. Folded in three must-fixes: (MF-1)
  inherited `-cflAbort`/`-recompute` would spuriously hard-abort a stable SMS run (reads
  un-augmented `ele->getMass()`) → **rejected at parse with SMS** (API warning); (MF-2)
  `map<nodeTag,Vector>` re-baseline is **unsafe under OpenSeesMP** shared/ghost nodes →
  added explicit parallel contract to M5 (v1 = interior-only/error on boundary, restore
  eagerly); (MF-3) "≥ Δt_target by construction" **overclaimed** → reworded to a
  conservatively-bounded *optimistic estimate*, validated empirically (Formulation,
  open-Q3). Should-fixes: Δtₑ/Δmₑ lump self-consistency; constrained-all-nodes no-op
  reporting; test-4 moved to an SSI model; corrected M2/M3 wording, the FictitiousMass
  precedent granularity (`getNodePtrs` not `getFEs`), and the moot recorder-bug
  contingency.

*(filled in once executing; move to `Ladruno_internal/implemented_selective_mass_scaling.md` when done)*
