# ADR 39 — P2 design: NTS penalty narrow phase (frictionless)

> Pre-code design for P2 — the first REAL mechanics. Goes through the FULL
> multi-agent adversarial gate BEFORE coding. Grounded in the ADR's Theory
> grounding (LS-DYNA §26.7 + de Souza Neto / Laursen / Wriggers). Parent:
> `39_ladruno_contact_domain_adr.md`; P1 done (`_adr39_p1_design.md`,
> `_adr39_p1b_design.md`); loop `_adr39_loop_state.md`.

## P2 goal

Turn the zero-force P1 adapters into real **node-to-segment (NTS) penalty
contact** (frictionless): closest-point projection → gap → penalty normal force
`F_c` (+ consistent tangent `K_c` for the implicit leg), assembled over **real
per-pair connectivity**. Still **brute-force** candidate pairing (bucket sort is
P2.5). Friction is P3.

## Kinematics — closest-point projection (LS-DYNA §26.6/26.10)

A slave node at `x_s = X_s + u_s`. A master segment with nodes `x_i = X_i + u_i`
and bilinear (quad-4) / linear (tri-3) shape functions `φ_i(ξ,η)`. The master
point `x̄(ξ,η) = Σ φ_i(ξ,η) x_i`, tangents `g_α = ∂x̄/∂ξ_α`, normal
`n = (g_1 × g_2)/‖g_1 × g_2‖`.

**Projection** = find `(ξ̄,η̄)` s.t. `(x_s − x̄)·g_α = 0` (orthogonality). Newton
on the 2×2 system (26.10), ≤10 iters; flat tri-3 → closed form (no Newton).
**Gap** `g_n = n·(x_s − x̄)`; penetration ⇔ `g_n < 0` AND `(ξ̄,η̄)` inside the
facet (with the §26.6 small overlap tolerance ~1.02 for edges).

## Force (penalty) + B-operator

Pair active (`g_n < 0`): normal traction magnitude `t_n = -k_n g_n` (>0).
- slave force `f_s = -t_n n` (pushes slave out);
- master node force `f_i = +φ_i(ξ̄) t_n n`.

By the orthogonality result (ADR Q-TAN), `δg_n = n·(δu_s − Σ φ_i δu_i)` — the
`∂ξ̄/∂u` term drops from the **residual**, so the force is FIRST-ORDER EXACT with
a once-per-step-frozen `(ξ̄,η̄)`. Define the gap operator `B` (1 × nDOF) over the
pair's DOFs `[u_s ; u_1..u_nseg]`: `g_n = B·U + g_n0`, `B = [ nᵀ | -φ_1 nᵀ ... ]`.
Residual `r = Bᵀ t_n` (assembled into the adapter's resid over its `getID`).

## Tangent `K_c` (implicit leg, scaled by c1)

Consistent tangent (frictionless, de Souza Neto / Laursen):
`K_c = k_n Bᵀ (n⊗n) B`  — the MAIN penalty term — **plus** the dropped
`O(g_n·κ)` curvature + `∂n/∂u`/`∂ξ̄/∂u` terms. **P2 ships only the main term**
(lagged geometric tangent, Q-TAN): exact for flat segments + small rotation,
degrades on curved/rocking/large-slide. Newmark assembles `addKtToTang(c1)` ⇒ the
adapter's `getTangent` returns `c1·K_c` (P1 gate MAJOR-3). Explicit never asks.

## `k_n auto` (LS-DYNA 26.14a + SOFT floor)

`k_n = f_si·K·A²/V` for a solid master segment (K = bulk modulus, A = segment
face area, V = element volume), `f_si` default **0.10**; with the SOFT=1 Courant
floor `k = max(½·SOFSCL·m*/Δt_c², f_si·K·A²/V)` (26.15). `-kn $val` overrides.
(Reading the master element's K/A/V via the surface→element link — Q-P2-kn.)

## Connectivity (P1 empty → P2 real) — the numbering shift

P2 adapters declare **real connectivity** = the pair's DOFs. **Granularity (P2):
one adapter per SLAVE NODE**, connectivity = `{slave node} ∪ {ALL master-surface
nodes}` (conservative static superset; bucket sort P2.5 prunes to neighbour
segments). The adapter, in `getResidual`, brute-force-scans its master surface for
the penetrated segment, projects, and assembles. Because connectivity is now
non-empty, the gate is **answer-to-1e-12** vs a reference (NOT bitwise — the
graph changes; P1 gate Q-P1-3). The frozen `getID` MUST be a superset of every
DOF any active pair can touch, else `addB` silently drops force (P1 gate / Q-WIRE).

## State (P2 introduces the first pair state)

`g_n0` (initial-gap capture for stress-free start) per pair, on the Domain-owned
`LadrunoContactDomain` (committed via the P1b `Domain::commit`/`revert` hooks). The
adapter is still a stateless VIEW reading pair state by key. Frictionless ⇒ the
only state is `g_n0` + the active/inactive flag; friction slip is P3.

## First sub-rung (de-risk before deformable–deformable)

1. **slave nodes vs a RIGID analytical plane** (normal + point, no master
   compliance) — penetration `= P/k_n` exact, unambiguous.
2. then **two deformable LadrunoBrick blocks**.

## Files

- `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — add pair state (`g_n0`,
  active flag), the projection + penalty kernel (or a `LadrunoContactKernel.h`),
  rigid-plane support.
- `SRC/domain/contact/LadrunoContactSurface.{h,cpp}` — resolve Node* + coords at
  setDomain; segment iteration; `getSegmentNodes(s)`.
- `SRC/analysis/handler/LadrunoContactFE.{cpp,h}` — real connectivity (per slave
  node) + real `getResidual`/`getTangent` (project, gap, penalty, B, K_c·c1).
  Move to `SRC/domain/contact/`? (Q-P2-loc — it needs Node trial disp + the
  surfaces; decide.)
- reuse: SimpleContact3D projection algebra (`:575-684`, RE-DERIVE as pure fns);
  `LadrunoEmbeddedKernel` scalar helpers.
- `tests/test_adr39_contact_p2.py`.

## P2 acceptance battery (the meaty one)

- rigid-plane: penetration `= P/k_n` exact (explicit + implicit).
- 2-block crush: penetration vs analytic 2-spring series.
- **`K_c` vs finite-difference on a ROTATED/curved config** (a flat crush passes
  even with `∂n/∂u` missing — false confidence; Q-TAN).
- frictionless **release/reopen → exact return to F=0** (active-set on/off).
- **Hertz** sphere/cylinder-on-halfspace → analytic `p(r)=p₀√(1−(r/a)²)`, peak
  `p₀` + radius `a` converge under refinement.
- **contact patch test** — report the NTS oscillation magnitude vs mesh ratio +
  confirm resultant equilibrium EXACT (NTS fails it — characterize, not pass/fail).
- `E_contact` bounded ≥ 0 frictionless; NEGATIVE ⇒ flag initial penetration.
- 1e-12 (not bitwise) vs an independent reference; explicit stable.

## Open questions for the gate

> [!question] Q-P2-proj
> Newton 2×2 projection robustness: non-convergence on near-degenerate segments,
> the concave-corner / edge cases (§26.6 fails at sharp concave corners). Fallback?
> Flat tri-3 closed form vs general Newton — support both or quad-only?

> [!question] Q-P2-tan
> Is shipping ONLY `k_n Bᵀ(n⊗n)B` (dropping `O(g_n·κ)`) acceptable for the P2
> implicit gate, given the gate runs `K_c`-vs-FD on a ROTATED config? Or must P2
> include the `∂n/∂u` term to pass FD on the rotated case? (This decides whether
> the "lagged tangent" claim survives its own test.)

> [!question] Q-P2-conn
> One-adapter-per-slave with connectivity = {slave ∪ ALL master nodes}: is the
> all-master-nodes superset too wide (dense K_c block, fill) even for P2's small
> tests, or fine until P2.5 prunes it? Risk: a slave whose active segment's nodes
> are (trivially) all in the superset — OK; but verify no active DOF escapes `getID`.

> [!question] Q-P2-active-set
> Active-set chatter: a pair flickering active/inactive across steps/iterations —
> in implicit Newton this can prevent convergence (non-smooth). Frictionless
> penalty is C0 at g_n=0; is that enough, or do we need a smoothing/persistence
> rule for the P2 implicit gate? (Explicit just steps through.)

> [!question] Q-P2-loc
> Does `LadrunoContactFE` move to `SRC/domain/contact/` (OPS_Domain) now that it
> needs the surfaces + projection kernel? Or stay in OPS_Analysis and call into
> OPS_Domain? (Cross-lib direction; P1 kept it in handler dir.)
