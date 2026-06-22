# ADR 39 — P2 design: NTS penalty narrow phase (frictionless)

> **Revised after the full 4-lens adversarial design gate** (Workflow wp3cr60mf →
> SALVAGEABLE-WITH-CHANGES). The gate caught a self-contradiction (the main-term-
> only tangent cannot pass the design's own FD-on-rotated gate) + two silent-wrong
> BLOCKERs. All folded in below, flagged `[GATE]`. Verdict + fixes archived in the
> loop log. Parent ADR `39_..._adr.md`; P1 done; loop `_adr39_loop_state.md`.

## P2 goal + split `[GATE: split mandatory]`

Real NTS penalty contact (frictionless). **Split into two shippable rungs** — the
rigid plane isolates penalty+assembly+active-set from the four deformable-path
risks (auto-kₙ reachability, projection robustness, normal orientation, faceted
tangent):
- **P2a** — slave nodes vs a **rigid analytical plane** (+ inclined plane),
  explicit + implicit, **`-kn $val`** (no auto). Connectivity = {slave} only,
  B = nᵀ, no projection kernel, no master mass. Ships clean.
- **P2b** — faceted-master projection rung → two deformable `LadrunoBrick` blocks
  + Hertz + `-kn auto` + SOFT floor + the **∂n/∂u** tangent (FD-on-rotated gate) +
  the `ZeroLengthContactASDimplex` cross-check.

Friction is P3; bucket sort P2.5.

## Kinematics — closest-point projection (P2b) `[GATE: bounded + oriented]`

Slave `x_s`, master segment nodes `x_i`, bilinear/linear `φ_i(ξ,η)`,
`x̄=Σφ_i x_i`, tangents `g_α=∂x̄/∂ξ_α`, `n=(g_1×g_2)/‖g_1×g_2‖`.

- **Projection** = Newton on `(x_s−x̄)·g_α=0` (26.10) — but **bounded** `[GATE
  BLOCKER-2]`: cap 10 iters; reject segment if `|detK| < ε‖g_1‖‖g_2‖`; return a
  **"no valid projection" sentinel** on non-convergence so the scan skips (never
  assemble garbage). Do NOT copy SimpleContact3D's unbounded `while` (`:600-635`).
  Flat tri-3 → closed-form (the de-risk path); quad-4 → bounded Newton.
- **Outward normal is DERIVED, not trusted from winding** `[GATE BLOCKER-1 — the
  #1 silent-wrong bug]`: at setDomain orient each segment normal from the master
  **solid-element centroid**: flip `g_1×g_2` so `n·(x̄ − x_elem_centroid) > 0`.
  (A flipped/inconsistent winding makes `gₙ<0` never fire → contact silently
  passes through, or some pairs attract.) Gate: flip a block's winding → identical
  force.
- **Concave-corner / multi-facet tie-break** `[GATE MAJOR-4]`: among segments with
  a valid IN-BOUNDS projection pick **max penetration** (most-negative gₙ); if none
  in-bounds but the node is geometrically inside, project to the nearest
  edge/vertex with the **averaged adjacent-facet normal**. "Node penetrates but no
  segment claims it" = a gate **assertion**, never a silent zero.
- Gap `gₙ = n·(x_s − x̄)`; penetration ⇔ `gₙ < 0` AND in-bounds.

## Force (penalty) + B-operator `[verified by the gate]`

Macaulay form `[GATE MINOR-15]`: `tₙ = kₙ⟨−gₙ⟩₊` (penetrating branch only; hard
zero otherwise — never write bare `−kₙgₙ`, the adhesion bug). Active pair:
- slave `f_s = −tₙ n`, master node `f_i = +φ_i(ξ̄) tₙ n`.
Gap operator (1×nDOF) over `[u_s ; u_1..u_nseg]`: `B = [ nᵀ | −φ_1 nᵀ ... ]`,
residual `r = Bᵀ tₙ`. **Self-equilibrated for any kₙ** (ΣF=0, ΣM=0 at the
converged projection) — gate-verified; add an explicit `|Σf|<1e-10·|F_c|` gate.
First-order-exact-with-frozen-ξ̄ holds for **interior** projections only `[GATE
MAJOR-10]`; at a clamped edge/vertex use the reduced edge normal or document the
O(error).

## Tangent `K_c` — the corrected, honest NTS tangent `[GATE BLOCKER-3, Q-P2-tan]`

`K_c = kₙ Bᵀ B  +  (∂n/∂u block)`, dropping ONLY the `O(gₙ·κ)` curvature term.
- **Formula fix `[T1]`:** the main term is `kₙ BᵀB`, NOT `kₙ Bᵀ(n⊗n)B` — `B`
  already contains `n`; the latter double-applies the projector.
- **Include ∂n/∂u `[T2 — resolves the self-contradiction]`:** for a FLAT segment
  under rigid rotation `∂n/∂u ≠ 0` (n rotates with the nodes), magnitude `O(tₙ)∝kₙ`
  — finite, NOT negligible. Main-term-only would FAIL FD-on-rotated by construction.
  The ∂n/∂u normal-variation block is the standard Laursen/Wriggers NTS tangent and
  is **symmetric** for frictionless (a δ²gₙ Hessian) → symmetric solver stays valid.
- FD-on-**flat-rotated** PASSES; FD-on-**curved** carries a documented `O(gₙ·κ)`
  residual that shrinks with refinement (assert convergence, not exact pass).
- `getTangent` returns **bare** `K_c`; Newmark scales by `c1` via `addKtToTang(c1)`
  (no internal pre-multiply — double-scale risk `[T7]`). Explicit never asks.
- **Implicit re-projects per Newton iteration; explicit freezes per step** `[GATE
  MAJOR-9]`: the reference re-projects in `update()` every iter, commits geometry
  at commitState. Recomputing `n` fresh while freezing ξ̄ violates the exactness
  assumption + injects a spurious residual moment — keep `(ξ̄,x̄,n,B)` a consistent
  set per evaluation.

## Active set `[GATE Q-P2-active-set]`

Explicit steps through (fine). **Implicit needs an explicit anti-chatter rule** —
freeze the active set per step (status from the predictor, held through the
iteration; the semismooth fix, recommended) OR mandate `algorithm NewtonLineSearch`.
Bake the choice into the implicit gate scripts or they won't converge.

## kₙ `[GATE MAJOR-6]`

- **P2a: `-kn $val` mandatory** (no auto). The rigid-plane/inclined/2-block gates
  set kₙ explicitly.
- **P2b: `-kn auto`** = `kₙ = f_si·K·A²/V` (26.14a, f_si=0.10) + SOFT floor
  `max(½·SOFSCL·m*/Δt_c², …)` (26.15, deformable-only — undefined for the rigid
  plane). **K/A/V has NO reachable Element API** (`Element.h:62` = only
  `getCharacteristicLength`): reconstruct at setDomain from `material->getInitialTangent()(0,0)`
  of a master GP + V from nodal coords, **cache it** (not reachable at getResidual),
  retain the master element pointer.

## State / g_n0 `[GATE MAJOR-5]`

Default = **no offset + ABORT on initial penetration** (matches the reference, keeps
the negative-E_contact diagnostic meaningful). Stress-free-relief (StagedStrain
analogue) is a **separate opt-in** that also offsets the release condition (release
when `gₙ ≥ gₙ0`) and makes the gate measure ABSOLUTE penetration. Do NOT ship the
offset as silent default. Pair state on the Domain-owned `LadrunoContactDomain`
(committed via the P1b hooks); adapter is a stateless view.

## Connectivity `[GATE Q-P2-conn]`

One adapter per slave; `getID = {slave} ∪ {all master-surface nodes}` is a valid
superset for a SINGLE STATIC master surface. **Enforce, don't assume:** surface
**immutable after setDomain** (hard error on post-analyze edit); cheap assertion in
getResidual that the projected segment's nodes ⊆ getID (`opserr`+skip, not silent
drop). One adapter ↔ one master surface.

## Location `[GATE Q-P2-loc]`

`LadrunoContactFE` STAYS in OPS_Analysis (handler dir); cross-lib direction
Analysis→Domain (the established direction). Put the reusable projection/penalty/
tangent math in a header-only OpenSees-free **`LadrunoContactKernel.h`** (mirrors
`LadrunoJ2Kernel.h`) so both libs include it with no link inversion.

## Reuse correction `[GATE MINOR-14]`

`ContactMaterial3D` has NO penalty-normal tangent to lift (its `getTangent` is a
Lagrange form: `(0,0)=0`, `(0,3)=1`). P2's `kₙBᵀB`+∂n/∂u is authored fresh; reuse
= the projection algebra only (re-derived as pure fns). Friction `C_ss`/`C_sl` = P3.

## Independent REFERENCE for the gate `[GATE: replaces the bogus "1e-12 abs"]`

Penalty has an `O(load/kₙ)` penetration baseline + Newton converges to `tolGap`, so
absolute 1e-12 is wrong. Three tiers at RELATIVE tol (matches the fork's ~1e-7 J2
precedent):
1. **rigid plane → analytic `g = P/kₙ`** (+ inclined `g = P(n·d̂)/kₙ`), rel 1e-8. [P2a]
2. **single slave vs single fixed master segment → cross-check a hand-placed
   `ZeroLengthContactASDimplex` pre-defined pair**, rel 1e-6 — **THE numerical
   oracle for the deformable kernel**. [P2b]
3. **Hertz** = the refinement-convergence benchmark; 2-spring series = order-of-
   magnitude sanity only (ignores master bending/NTS distribution). [P2b]

## Acceptance battery

**P2a:** penetration `=P/kₙ` (explicit+implicit, rel 1e-8); inclined `=P(n·d̂)/kₙ`;
release/reopen → exact F=0; **sign** (penetration→restoring, never attract);
explicit 500-step stable; implicit converges under the anti-chatter rule;
E_contact load–unload returns to 0 `[MINOR-17]`.
**P2b:** faceted-master projection rung first; 2 deformable blocks; **`K_c`-vs-FD
on flat-ROTATED** (PASS) + on CURVED (O(gₙ·κ) shrinks w/ refinement); Hertz
`p(r)=p₀√(1−(r/a)²)` converges; ASDimplex cross-check (rel 1e-6); **self-equilibrium
`|Σf|<1e-10·|F_c|`**; **oblique 30°/45° plane** (off-axis n⊗n); patch-test
oscillation characterization + exact resultant; mesh-ratio sweep (master coarse vs
fine); winding-flip → identical force; **convergence-RATE** gate (pre-declared).

## Files

- `SRC/domain/contact/LadrunoContactKernel.h` — NEW header-only pure-fn math
  (project, gap, penalty, B, K_c incl. ∂n/∂u, normal-orientation), OpenSees-free.
- `LadrunoContactDomain.{h,cpp}` — pair state (gₙ0, active flag, rigid-plane def),
  surface→element link + cached K/A/V (P2b), immutability enforcement.
- `LadrunoContactSurface.{h,cpp}` — Node*/coords at setDomain, segment iteration,
  derived outward normal, master-element pointer.
- `LadrunoContactFE.{cpp,h}` — real connectivity + real getResidual/getTangent via
  the kernel; per-step (explicit) / per-iteration (implicit) projection.
- `tests/test_adr39_contact_p2a.py`, `tests/test_adr39_contact_p2b.py`.

## Resolved open questions (all closed by the gate)
- Q-P2-proj → bounded Newton + sentinel + derived normal + tie-break; both tri/quad.
- Q-P2-tan → **add ∂n/∂u** (drop only curvature); main term `kₙBᵀB`. The decisive fix.
- Q-P2-conn → superset OK for one static surface; ENFORCE immutability + assert.
- Q-P2-active-set → freeze-per-step (or NewtonLineSearch) for the implicit gate.
- Q-P2-loc → stay OPS_Analysis; shared header-only kernel.
