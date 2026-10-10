---
wp: LEGACY
title: "Mortar friction committed slip is last-writer-wins (order-dependent) at SHARED slave nodes — fenced to matched/explicit C3.1, must be guarded before non-matche…"
legacy_seq: 125
---
### Mortar friction committed slip is last-writer-wins (order-dependent) at SHARED slave nodes — fenced to matched/explicit C3.1, must be guarded before non-matched friction
- **Bites:** ADR-41 C3.1 mortar friction at a slave node shared by ≥2 (slave-facet, master-facet) pairs.
  The per-global-node committed slip `st.gpTtrial` (`LadrunoContactFE::addMortarFriction`) is a plain
  OVERWRITE: each facet visiting the node computes its OWN LOCAL `gbarT` (its own clip/projection) and the
  LAST facet evaluated in the residual sweep wins the committed slip. The *force* is still deterministic
  (every facet reads the same read-only committed `gpT`, so `R(u)` is clean — no singular solve), but the
  committed plastic slip carried to the next step depends on FE-tag ordering. The normal gap dodged this
  with an idempotent delta-accumulator keyed `(c,node,feTag)` (`accumulateMortarGap`); the friction slip has
  no equivalent because the slip is a return-map OUTPUT, not a linear accumulation.
- **Why it's fenced (for now):** C3.1 ships matched-facet + explicit (CDL) only — one facet per node, so the
  race never fires (the battery is matched). It is within the design's accepted "standard-basis LOCAL
  approximation at shared nodes" ([[_adr41_c3_design]]). But it is UNGUARDED and untested for non-matched
  friction. **Before C4 / non-matched frictional meshes:** add a shared-node friction regression + either a
  per-(node,feTag) slip reconciliation or an explicit area-weighted blend. Found by the C3.1 adversarial gate
  (MAJOR-1, #377). **C3.3 update (#379):** the tangential multiplier `λ_T` (committed from `lambdaTtrial`,
  written per-facet last-writer-wins) and `gpT` BOTH inherit this — unlike the normal `λ_N` Uzawa, which
  augments from the order-INDEPENDENT global accumulator `gtGlobal/aGlobal`. So the per-node `λ_T`/`gpT`
  reconciliation is the same single fix for the whole friction state. Still fenced to matched/explicit; the
  C3.3 gate (MINOR-1) re-confirmed it is inherited, not introduced.
  **WP-155 update (2026-09-30, review #897):** the fence is now crossed in practice. The apeGmsh pile
  ladder runs μ = 1 and cohesion-only mortar on a NON-matching skin/hole, and its cohesion-only
  multi-step failure persists under `-augment never` (so not the Uzawa) and is plausibly THIS defect.
  `-augment never` removes only the λ_T half. See [[155_pile_contact_r05]] §5–§6. Recommended as the
  next contact slice.
  **WP-157 (2026-09-30) — RESOLVED for friction.** The state is now per (slave node, slave facet, master
  facet): `MortarFrictionState` keyed (contactTag, node, sf, mf), so each pair re-reads what its own
  return map wrote. The second failure mode surfaced along the way: on a CURVED/creased interface a
  pair also read a `λ_T`/`gpT` lying in a NEIGHBOUR's tangent plane (a normal leak |t·n| ≈ cap from step
  2 on; oracle T3). Pinned by `test_adr157_mortar_pair_friction` (creased roof: analytic force and
  multi-step force-control convergence; both fail on d63f49750). See [[157_mortar_friction_pair_state]].
  **Lifecycle change (review #900, finding 3).** A pair refused by `-maxGap`
  for an epoch loses its friction state and re-engages fresh (D4). A pair that
  stays paired but inert keeps `λ_T`/`gpT` frozen and re-applies them when it
  re-enters. Before WP-157, siblings sharing the node refreshed that state.
- **C4 update (#381) — RESOLVED for the TIE path; STILL FENCED for FRICTION.** C4 mesh-tying hits shared
  slave nodes immediately (non-matching meshes are the whole point), so the pre-req had to be discharged
  before relying on it. The tie state (`λ_tie`, the full 3-vec relative displacement `r_I`) does NOT inherit
  the bug, because `r_I = Σ D u_s − Σ M u_m` is a **LINEAR accumulation** — not a return-map OUTPUT — so it
  uses the SAME order-independent global accumulator `λ_N` already uses (`accumulateMortarTie` delta-update
  keyed `(c,node,feTag)` → `rtGlobal`, Uzawa'd in `commit()` no-clamp). The FORCE reads each facet's LOCAL
  `r` (deterministic R(u), the C2.2 rule); the GLOBAL `r` feeds only the commit Uzawa + the `‖r‖` query.
  Pinned by `test_adr41_mortar_c4_1::test_c4_1_shared_node_order_independent` (a slave node shared by 2 tie
  facets ⇒ either facet order gives a BIT-identical converged solution) + oracle T6. **The FRICTION
  last-writer-wins (gpT/λ_T) is STILL fenced** — C4 is mutually exclusive with friction, so the tie never
  touches the friction slip; non-matching FRICTIONAL meshes still need the per-(node,feTag) slip
  reconciliation (an area-weighted blend, since the slip is a return-map output, not linear).
