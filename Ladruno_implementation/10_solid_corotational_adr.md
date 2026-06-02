---
title: ADR — Corotational geometry method for solid elements (SolidTransformationCorot)
project: Ladruno
status: proposed
priority: high
owner: nmora
tags:
  - implementation
  - element
  - solid
  - kinematics
  - corotational
  - adr
---

# ADR — Corotational geometry method for solid elements

**Status:** proposed (2026-06-01) · **Supersedes:** the "v2 corot" placeholder in
[[solid_transformation_wrapper]] · **Consumer:** [[09_ladruno_brick]] · **Seam
interface:** `SRC/element/solidTransformation/SolidTransformation.h`

---

## 1. Context

The solid geometry-method layer ships two of three methods:

| method | strain regime | rotation regime | material consumed | status |
|---|---|---|---|---|
| `linear` (v1) | small | small | small-strain lib (`setTrialStrain`) | **merged** #71 |
| `corot`  (v2) | **small** | **large** | small-strain lib (`setTrialStrain`) | **this ADR** |
| `finite` (v3) | large | large | `FiniteStrainNDMaterial` (`setTrialF`) | **merged** #76 |

`corot` is the missing capability. Its slot in the factory is reserved
(`METHOD_COROT = 1`) and commented out
(`SolidTransformationLinear.cpp:94`); no `SolidTransformationCorot.{cpp,h}`
exists in any branch.

**The defining purpose (settled with the user 2026-06-01).** Corotational exists
to deliver geometric nonlinearity *while reusing the existing small-strain
material library* — `ElasticIsotropic3D`, `J2Plasticity`, Drucker-Prager, etc.,
all driven by `setTrialStrain(ε)`. It strips the (large) rigid-body rotation so
the constitutive core only ever sees *small* strain in a co-rotating frame.

> **Why not TL.** A Total-Lagrangian finite-strain element (Green-Lagrange `E`,
> 2nd Piola-Kirchhoff `S`, `∂S/∂E`) was considered and rejected for this slot.
> The user's framing settles it: *"What would be the sense of having corot for a
> TL material, if the element consumes the tensor directly?"* — i.e. if a
> material already hands back a finite stress/tangent tensor, the spatial
> `finite` path (v3, already shipped) consumes it directly and corot adds
> nothing. Corot's entire value is material reuse without a finite material. A
> TL material track, if ever wanted, is a *separate* ADR and does **not** use
> this wrapper. Note also: the merged `FiniteStrainNDMaterial` is **spatial**
> (`getStress`→Cauchy σ, `getSpatialTangentTensor`→spatial `c`); it exposes no
> PK2 / `∂S/∂E` hook, so TL is not even reachable from today's material side.

**Seam readiness (from the 2026-06-01 codebase assessment).**
- `globalizeForce` / `globalizeStiff` have **no call site yet** — `linear` and
  `finite` are both identity pass-throughs, so `LadrunoBrick::formResidAndTangent`
  never calls them. **Corot is the first method that actually uses seam 3.** The
  wiring must be added (and must stay a no-op for `linear` → byte-identical).
- `localizeDisp` *is* wired, via `LadrunoBrick::computeLocalDisp()`
  (`LadrunoBrick.cpp:904-922`).
- The identity aliasing guards (`if (&out != &in)`) in Linear/Finite are safe
  only for identity. Corot's `R̃ Pᵀ k Pᵀ R̃ᵀ` reads all of `kCore` while writing
  `kGlobal`, so the **wrapper owns its own scratch buffers** and reads the core
  matrix fully before writing the result (the element may keep passing one
  buffer).
- Method id is serialized as `… + 1000*getMethodID()` into `idData(28)`
  (`LadrunoBrick.cpp:2217-2220`); `recvSelf` rebuilds via
  `SolidTransformation::create(id)`. Corot = `1000`.

**Reusable assets in-tree** (no need to roll our own SO(3) algebra):
- `SRC/element/shell/ASDEICR.h` — element-independent corotational projector /
  spin matrices (Felippa CU-CAS-00-03; Felippa & Haugen CMAME 2005). **The
  template.** `ASDShellQ4CorotationalTransformation`'s
  `calculateLocalDisplacements` / `transformToGlobal` map 1:1 onto our
  `localizeDisp` / `globalizeForce`+`globalizeStiff`.
- `SRC/matrix/Versor.h`, `GroupSO3.h` — quaternion + SO(3) exp/log, `Spin`,
  consistent-linearization tangent operators. Use for the algebra once `R` is in
  hand.
- **Gap:** there is *no* polar-decomposition-of-`F` helper anywhere in the tree
  (`Versor::from_matrix` extracts `R` from an *already orthogonal* matrix). We
  write the polar primitive.

---

## 2. Decision

Implement **`SolidTransformationCorot`** as an **element-independent
corotational (EICR)** geometry method:

1. `getStrainMeasure()` → `SmallStrain` (the element runs its normal
   `std`/`bbar`/`uri`/`eas` B-matrix path unchanged on the de-rotated displacement
   and feeds `ε` to ordinary small-strain materials via `setTrialStrain`).
2. `update(refCrds, curCrds)` → extract a single **element-level** rigid rotation
   `R` by **iterative Higham polar decomposition** of a best-fit deformation
   measure (see §3 for the best-fit choice). Store `R` and the centroid.
3. `localizeDisp(uGlobal, uCore)` → produce the **deformational** displacement
   `u_d` by removing rigid translation (centroid) and rigid rotation `R`:
   `u_d,a = Rᵀ(x_a − x_c) − (X_a − X_c)`. *Real work, real buffers.*
4. `globalizeForce` / `globalizeStiff` → EICR transform back to global:
   `f = R̃ Pᵀ f_d` and `k = R̃ (Pᵀ k_d P + K_geo) R̃ᵀ`, where `P` is the EICR
   rigid-body projector and `K_geo` is the rotational/equilibrium geometric
   stiffness built from the corotated internal force `f_d` (hence `fCore` is
   passed to `globalizeStiff`). Reuse `ASDEICR.h` for `P`, `Spin`, and the
   geometric blocks.
5. `frameTimeVarying()` → `true`; `getCurrentFrame(R)` → the current `R` (for the
   recorder local frame).

### 2a. Sub-decision — element-level single `R` (NOT per-Gauss-point)

This is the "best-fit-F scope" question flagged in memory. **Decision:
element-level single `R` per element**, matching EICR for shells/frames.

- *Element-level (chosen).* One rigid rotation for the whole element; cheaper
  (one polar per element per iteration); is the textbook corotational element
  (Crisfield, Felippa); the residual within-element deformation (incl. bending
  modes) is carried by the de-rotated small-strain core. Valid precisely in the
  corot regime (small strain, large rotation).
- *Per-GP `R` (rejected for v2).* Polar-decompose `F` at each Gauss point. More
  accurate when there are large rotation *gradients* inside one element, but it
  is essentially finite-strain-lite, costs a polar per GP per iteration, and
  blurs the line with the `finite` path we already have. If per-GP rotation
  matters, the answer is mesh refinement or `-geom finite`, not per-GP corot.

### 2b. Solid simplification — no rotational DOFs

Brick nodes carry **only translational DOFs** (3/node). The heavy machinery in
`CorotCrdTransf3d` (incremental nodal triads, quaternion composition, Spurrier
extraction per node) exists to advance *rotational* DOFs — solids have none.
Consequences:
- `R` is a pure function of current nodal *positions* ⇒ **stateless**. It is
  recomputed from `(refCrds, curCrds)` every evaluation; **no committed triad**,
  so `commitState`/`revertToLastCommit`/`revertToStart` stay the base no-ops and
  `sendSelf`/`recvSelf` need no extra corot state.
- This is materially simpler than the shell/frame EICR and removes the
  rotation-update bookkeeping that dominates `CorotCrdTransf3d`.

---

## 3. Open design choice inside the decision — the best-fit measure for `R`

Two equivalent-quality routes to the element rotation; pick one in
implementation (both polar-decompose via Higham; neither uses eigen-of-`U`, so
repeated principal stretches are a non-issue, `R` unique for `det>0`):

- **(A) Procrustes best-fit of nodal positions.** `H = Σ_a (x_a−x_c)(X_a−X_c)ᵀ`,
  then `R = polar(H)`. Pure EICR; independent of shape functions; standard.
- **(B) Volume-averaged deformation gradient.** `F̄ = (1/V)∫ F dV` (reuses the
  brick's existing `F` shape-gradient machinery from the `finite` path),
  `R = polar(F̄)`. Ties into code we already have.

**Recommendation: (A) Procrustes** — it is the canonical EICR best-fit, needs no
quadrature, and is the cleanest match to `ASDEICR.h`. Revisit only if a patch
test exposes a frame-invariance defect.

---

## 4. Consumer wiring (LadrunoBrick) — the integration points

1. **Add the seam-3 calls** to `formResidAndTangent` (the `std`/`bbar`/`uri`/`eas`
   path, `LadrunoBrick.cpp:747+`): after the core `f`/`k` are assembled in the
   corotated frame, call `theGeom->globalizeForce(f, f)` and
   `theGeom->globalizeStiff(k, f, k)`. For `linear`/`finite` these are identity ⇒
   **byte-identical preserved**; verify the `-geom linear` 1e-12 reproduction
   gate still passes after adding the calls.
2. **`localizeDisp`** is already wired (`computeLocalDisp`); corot returns the
   deformational `u_d` there, so the existing B-matrix loops need no change.
3. **Factory:** uncomment `case METHOD_COROT: return new SolidTransformationCorot();`
   in `SolidTransformationLinear.cpp` and `#include` the corot header.
4. **CMake/Makefile:** add `SolidTransformationCorot.cpp` to
   `SRC/element/solidTransformation/{CMakeLists.txt,Makefile}`.
5. **Tcl/Python parse:** map `-geom corot` → `geomMethodID = METHOD_COROT` in the
   brick's input parsing (wherever `-geom finite/linear` is parsed).

---

## 5. Acceptance gates (the test battery — Zone-A)

The defining corotational gates, in priority order:

1. **Rigid-body rotation ⇒ zero stress (THE gate).** Rotate the element through an
   arbitrary large `R` (e.g. 90°, 179°) with zero deformation: internal force and
   every Gauss-point stress must be ~0 (≤1e-10). The EICR projector `P` must
   guarantee self-equilibrium.
2. **Consistent tangent / quadratic Newton.** Finite-difference the global tangent
   against the residual at a deformed, rotated state; must match to ~1e-6. (This
   is what forces the full EICR `K_geo`, not just `R̃ Pᵀ k_d P R̃ᵀ`.)
3. **`-geom linear` byte-identical** after wiring the globalize calls (1e-12 vs
   pre-wiring).
4. **Small-rotation limit ⇒ linear.** Corot under tiny rotation reproduces `linear`.
5. **Large-rotation cantilever benchmark.** A slender elastic cantilever / Euler
   elastica vs. analytical or a fine `finite` reference — the headline physics.
6. **Plasticity in the rotated frame.** `J2Plasticity` core under large rotation +
   small plastic strain: yields correctly in the corotated frame (proves material
   reuse).

---

## 6. Consequences

**Positive.** Fills the last gap on the geometry-method axis; unlocks
large-displacement analysis (buckling, snap-through, large-rotation frames-of-
solids) for the *entire* small-strain material library with zero new material
code; reuses `ASDEICR.h`/`Versor`/`GroupSO3`; stateless ⇒ no serialization or
commit complexity.

**Negative / risks.**
- First real exercise of seam 3 — the globalize wiring is genuinely new code in
  the brick; the `linear` byte-identical gate is the safety net.
- **EAS × corot** interaction is the fragile corner (EAS condensation happens in
  the corotated frame; the condensed enhanced modes must ride along correctly).
  *Scope decision:* v2 ships **`std` + `bbar` first**; `uri` and `eas` under corot
  are a follow-up once the std gate battery is green.
- Writing the polar-of-`F` primitive (no in-tree helper). Higham Newton iteration
  is short and robust; unit-test it standalone against known `R`.
- Accuracy degrades as strain grows (corot assumption) — documented limit, not a
  bug; `-geom finite` is the escape hatch.

---

## 7. Implementation plan (step order)

1. **Polar primitive** — `Higham` Newton polar decomposition `F = R·U` as a small
   static helper (its own unit test vs. synthetic `R`). No eigen.
2. **`SolidTransformationCorot.{h,cpp}`** — implement the five seams (§2), EICR
   `P`/`K_geo` via `ASDEICR.h`, internal scratch buffers. `getStrainMeasure` =
   `SmallStrain`, stateless commit/revert/serialize.
3. **Factory + build** — uncomment `create()` case, add to CMake/Makefile.
4. **Brick wiring** — add `globalizeForce`/`globalizeStiff` calls to
   `formResidAndTangent`; `-geom corot` parse; confirm `linear` 1e-12.
5. **Test battery** — gates §5, Zone-A (`std`+`bbar`).
6. **Ledgers + banner** — `LEDGER_implementations.md` row (class-agnostic: no new
   classTag, wrapper is not a `MovableObject`); `LEDGER_quirks.md` for the
   seam-3-first-consumer wiring note; banner line if user-facing.
7. **PR on a fresh branch off `ladruno`** (per CLAUDE.md `--base ladruno`).

**Deferred to follow-ups:** `uri`/`eas` under corot; per-GP rotation (only if a
benchmark demands it); incremental-frame tracking (only if rotation >π ever
misbehaves — unlikely with total-geometry polar); the consistent-tangent
completion (§8).

> **STATUS UPDATE (2026-06-02): v2.0 SHIPPED — PR #88, merged to `ladruno`
> (squash `13a15ea9`).** `SolidTransformationCorot` with the §2 design; 14
> Zone-A tests green incl. ADR Gate 5 (load-driven cantilever, matches `-geom
> finite` <1% at ~26% tip) and Gate 6 (J2-plasticity-in-rotated-frame). An
> 8-lens adversarial-review sweep (33 confirmed findings; core math re-derived
> & confirmed) hardened it. The v2.0 tangent is the reduced form described in
> §2 — the completion is §8 below.

---

## 8. v2.1 — consistent-tangent completion (numerical-tangent entry point)

**Problem.** The v2.0 tangent is `K = R̄ k_d R̄ᵀ + (G1 + G1ᵀ)`; it omits two
O(strain) geometric terms — the **polar-Hessian** `∂Γ/∂u · m` and the
**material-frame ∂R coupling** `R k_d (∂Rᵀ/∂u) x⁰`. Because the *residual* is the
exact energy gradient, Newton still converges to the correct equilibrium, but at
a **linear/super-linear rate** — large-rotation runs need load sub-stepping (the
Gate-5 cantilever needs ~20 increments where `-geom finite`'s full tangent takes
far fewer). v2.1 closes the gap.

### 8.1 Decision — ship the *numerical* tangent first

Rather than first derive the algebraically heavy, degeneracy-prone analytic
polar-Hessian, the v2.1 entry point is an **element-level numerical (perturbation)
tangent**, exposed as a debug/robustness flag. This is uniquely clean here because
**the corot residual is the exact gradient**, so a finite difference of it *is* the
full consistent tangent — including both deferred terms — with no new derivation.

It serves three purposes at once:
1. **Troubleshooting isolator.** Swap to it when an analysis converges badly:
   converges well ⇒ the reduced analytic tangent was the cause (expected);
   still struggles ⇒ the tangent is **not** the culprit — look at the residual,
   conditioning, or the model. That separation is the headline debugging value.
2. **Robustness escape hatch.** Delivers full-tangent (≈quadratic) convergence
   *today*, at a compute cost, for stiff / strongly geometrically-nonlinear runs.
3. **Validation oracle.** When the analytic polar-Hessian + material-∂R terms
   eventually land (§8.4), they MUST match this numerical tangent — it is the
   arbiter, reusing the existing FD-tangent gate machinery.

### 8.2 Design / behaviour

- **New element option `-tangent {analytic|numerical}`**, default `analytic`
  (v2.0 behaviour, byte-unchanged). Parsed in `OPS_LadrunoBrick.cpp` alongside
  `-geom`/`-formulation`. Primarily meaningful for `-geom corot`; for
  `linear`/`finite` (already-consistent tangents) it is a debug cross-check only.
- **Where it lives: `LadrunoBrick`, NOT the wrapper.** The seam `globalizeStiff`
  only receives `k_d`/`f_d` and cannot re-evaluate the residual at perturbed
  DOFs; the element can. Add `LadrunoBrick::formNumericalTangent()`:
  ```
  for each dof d in [0, 3·numNodes):
      perturb node trial disp:  u_d ± h
      K[:,d] = ( R(u+h e_d) − R(u−h e_d) ) / (2h)     // R(·) = global resisting force
      restore u_d
  optionally symmetrise
  ```
  where `R(·)` runs the element's **own** residual path (`computeLocalDisp` →
  `theGeom->update` refreshes the corot frame → material `setTrialStrain` →
  `globalizeForce`). Central difference. This captures **all** dependencies —
  the rotation `R`, the spin map, and the material response — i.e. the true
  consistent tangent.
- **Perturbation size:** `h = h_rel · (1 + ‖u‖∞)`, `h_rel ≈ 1e-7…1e-6` (central
  diff ⇒ O(h²) truncation balanced against round-off; the FD-tangent test already
  uses `h=1e-6`).
- **State hygiene (critical):** perturbation touches **trial** state only —
  `setTrialStrain` is trial, `commitState` is never called inside the loop, and
  each DOF is restored after its ±h pair. Committed material state (and the
  stateless corot frame, recomputed each call) is left intact. For a
  path-dependent inner material (J2) this yields the algorithmic/consistent
  tangent, exactly as the material's own `getTangent` would.
- **Serialization:** pack the tangent-mode int next to `formulation`/`geom`/
  `mass`/`hourglass` in `idData(28)` (analytic = 0 ⇒ byte-identical stream), and
  rebuild on `recvSelf`.
- **Cost:** `2·(3·numNodes)` residual evals per tangent — **48 for the 8-node
  hex**. Acceptable for debugging or a hard model; **not** a default.

### 8.3 Diagnostics (complementary)

- **`eleResponse(tag, "tangentGap")`** → `‖K_analytic − K_numerical‖ / ‖K‖` at
  the current state — a developer probe to *measure* how far the reduced tangent
  is from consistent (a regression watch as the analytic terms land).
- **Help/doc note** on `-geom corot`: convergence trouble ⇒ try smaller load
  steps, `-tangent numerical`, or `-geom finite`. (An auto-warning is impractical
  — the element cannot observe Newton non-convergence; keep it in the docs +
  the `tangentGap` probe.)

### 8.4 Acceptance gates (Zone-A)

1. **Correctness:** `-tangent numerical` gives the **same** converged Gate-5
   cantilever tip displacement as `-tangent analytic` and as `-geom finite`
   (the tangent must not change the answer — only the path there).
2. **Robustness payoff:** the cantilever converges in **strictly fewer load
   steps / total iterations** with `numerical` than with `analytic` (e.g.
   converges at `nsteps` where `analytic` diverges, or far fewer Newton iters).
3. **Consistency:** at small strain (where the v2.0 analytic tangent is already
   exact) `numerical ≈ analytic` to ~1e-5 (sanity that both code paths agree).
4. **State safety:** a committed multi-step run with `-tangent numerical` matches
   `-tangent analytic` step-for-step (the perturbation loop leaves committed
   material/plastic state uncorrupted) — exercised with the J2 material.

### 8.5 Follow-up — the analytic terms (true v2.1)

Once the numerical tangent is the validated oracle, derive and add the two
deferred analytic terms in `SolidTransformationCorot::globalizeStiff`:
- material-frame coupling `R k_d (∂Rᵀ/∂u) x⁰` (reuses the spin map Γ already
  computed in `update()`), and
- the polar-Hessian `∂Γ/∂u · m` (second derivative of the Procrustes polar; the
  hard, degeneracy-prone piece — branch on stretch multiplicity).
Gate: tighten the FD-tangent test from the documented bound to ~1e-6 at finite
strain, and confirm the cantilever matches `numerical`'s step count. This makes
quadratic Newton the default `analytic` behaviour at no per-tangent cost, and the
`-tangent numerical` flag then survives purely as the debugging isolator (purpose
1) and the oracle.

### 8.6 Implementation plan

1. `-tangent {analytic|numerical}` parse + serialize (default analytic, byte-safe).
2. `LadrunoBrick::formNumericalTangent()` (perturb-restore central FD over the
   residual path); route `getTangentStiff`/`formResidAndTangent` to it when the
   mode is `numerical`.
3. `eleResponse "tangentGap"` + doc note.
4. Zone-A tests (§8.4): correctness, robustness payoff (step count), small-strain
   consistency, committed-state safety with J2.
5. Ledger row update + PR off `ladruno`.
6. (Later) the analytic terms (§8.5) + tightened FD gate.
