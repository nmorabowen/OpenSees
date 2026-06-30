# LadrunoTie — kinematic mortar mesh-tie via constraint emission (the projection handler)

> ADR-62. Status: **PLAN / pre-code — the constructive successor to the shelved ADR-61.**
> Family: ADR-30 (LadrunoProjectionHandler — the enforcement, SHIPPED) · ADR-41 (mortar
> D/M + `-tie` C4 — the pairing) · ADR-39 (ContactDomain bucket-sort + projection) ·
> ADR-61 (contact bipenalty — SHELVED; the penalty route this replaces). Next free ADR
> slot is 62 (60 held for finite-sliding re-emission; 61 = bipenalty). Author: N.
> Mora-Bowen (Ladruno).

---

## What

Enforce a **non-conforming explicit mesh-tie** as **kinematic constraints routed through
the shipped `LadrunoProjectionHandler` (ADR-30)** — *not* a penalty (SOFT / ALM /
bipenalty). A mesh-tie is **static** (a permanent bond that never forms, breaks, or
slides), so we pair the slave surface to the master surface **once at setup**, compute the
node-to-facet shape-function weights, and **emit one ordinary OpenSees `MP_Constraint` per
slave node** (`u_s = Σ_i N_i(ξ̄) u_{m,i}`). The user runs `constraints LadrunoProjection`;
the projection handler then enforces the tie by **M-orthogonal acceleration projection** —
**exact, momentum-conserving, Δt-neutral, keeping the `system Diagonal` recipe**.

This is the enforcement strategy every mature explicit code uses for ties, and the one the
penalty-based ADR-61 (shelved) was the wrong tool for:
- **LS-DYNA** — tied interfaces use the **kinematic constraint method** ("the first
  approach is now used for tying interfaces", Theory §26.1–26.2, §26.9), *not* penalty.
- **Abaqus** — surface `*TIE` is a **kinematic MPC** (exact slave-DOF handling, TG §6.6).
- **Kratos** — ties are `MasterSlaveConstraint`s (DOF-level master-slave relations).

The fork already owns **both halves**: the mortar/NTS pairing + shape weights (ADR-39/41)
and the kinematic enforcement (ADR-30). **The only new code is the glue** — a setup-time
constraint *generator*. There is **no new enforcement code and no runtime contact machinery
for a tie.**

---

## Why

ADR-61 established that penalty-based explicit ties are a bad trade:
- **SOFT** (soften `k`): leaves residual penetration `δ/h = ε_iface·CFL²/(4·α_m·SOFSCL)`.
- **Bipenalty** (add mass): the P0 oracle showed the reduced-mass bound `μ ≤ min(m_s,m_m)`
  forces **~100× two-sided interface-mass inflation** for a deformable–deformable tie —
  shelved.
- **ALM**: exact, but **implicit only**.

All three are workarounds for using a *penalty* where the physics wants a *constraint*. A
tie is a **bilateral, permanent, non-sliding** bond — i.e. a linear MP-constraint
`u_s − Σ N_i u_{m,i} = 0`. Enforced kinematically it is **exact** (no penetration, unlike
SOFT), **Δt-neutral** (no penalty spring → cannot raise `ω_max`, unlike a stiff penalty),
**momentum-conserving** (the projection is, by construction), and adds **no fictitious
mass** (unlike bipenalty). It dominates every penalty variant for the explicit tie — and
reuses shipped, validated code.

The key realization ADR-61 missed: **for a tie you should not penalize at all.** The
reduced-mass wall that killed bipenalty only exists for a *penalty spring*; kinematic
projection has no spring, so the wall does not exist.

---

## Where

### New code — a setup-time constraint generator (no enforcement code)

- **`SRC/.../LadrunoTie` generator** (command + a small builder) — at setup: (1) reuse the
  ADR-39 bucket-sort + closest-point projection to pair each **slave node** to **one master
  facet** (tri-3/quad-4); (2) evaluate the master shape functions `N_i(ξ̄)` at the slave's
  projection point; (3) emit, per slave node, an OpenSees **`MP_Constraint`** (one per
  translational DOF) with constraint matrix `Ccr = [N_1 … N_{nps}]` tying the slave to the
  master facet nodes. Solid–solid (ndf 3) first.
- **Reuse, unchanged:** `LadrunoProjectionHandler` (ADR-30, HANDLER 33001) for enforcement;
  the mortar/NTS projection + shape-weight kernels (ADR-41); the bucket-sort (ADR-39).
- **Parser:** `LadrunoTie -slaveSet … -masterSurface …` (Tcl + Python), emitting standard
  `MP_Constraint`s the existing `constraints LadrunoProjection` consumes. Optionally surface
  it as an enforcement mode on the existing mortar tie: `contact … -tie -kinematic`.

### Modify vanilla — NONE expected

The generator emits standard `MP_Constraint` objects via the public Domain API (the same
objects `equalDOF`/`rigidLink` create). No upstream edit; any touch → `LEDGER_vanilla_files`.

### classTags

**None** — the generator creates ordinary `MP_Constraint`s (no new tag); enforcement reuses
`HANDLER_TAG_LadrunoProjectionHandler 33001`. (If the generator is later realized as a
distinct constraint-set object needing a tag, the free contact ELE slot is **33016**.)

---

## How

### The constraint (per slave node, captured once)

Slave node `s` projects onto master facet `f` at parametric `ξ̄` (the ADR-39 closest-point
projection). The tie is the **collocation (node-to-segment) relation**

```
u_s = Σ_i N_i(ξ̄) u_{m,i}          (per translational DOF; N_i = master facet shape fns)
```

emitted as an `MP_Constraint` with `Ccr = [N_1 … N_{nps}]`, retained DOFs = the master facet
nodes, constrained DOF = the slave node. The projection handler then solves, per
connected-component group, `a_proj = L(LᵀML)⁻¹LᵀM a_raw` with `L=[I;C]` — exact enforcement
of the tie at the acceleration level every step, momentum-conserving, on the untouched
diagonal mass.

**Why node-to-segment (collocation), not integral-mortar weights:** a node-wise constraint
ties each slave node to **one** master facet, so the slave sets are disjoint and **no slave
node appears in two constraints** — which is exactly what the projection handler v1 requires
(it **refuses MP-chains and doubly-constrained DOFs**, §BLOCKER-1). True integral-mortar
weighting (`∫N dΓ`, two-sided, dual basis) couples slave nodes to each other → MP-chains →
needs the handler's deferred chain support; deferred to P2. Collocation is the standard,
robust v1 tie (it is exactly how Abaqus surface-to-surface `*TIE` and LS-DYNA
`*CONTACT_TIED_NODES_TO_SURFACE` collocate) and is **exact for matching meshes, optimal for
non-matching**.

### Decision 1 — Static pairing (setup-time emission), not a runtime contact handler

A mesh-tie never forms/breaks/slides, so the pairing + weights are computed **once at
setup** and frozen into `MP_Constraint`s. Consequences:
- **No custom contact constraint handler** for the tie — the generator runs at model-build,
  the shipped `LadrunoProjectionHandler` does enforcement. (This sidesteps the
  "only-one-constraint-handler" composition problem entirely: the tie *is* just constraints.)
- **No `LadrunoContactFE`, no per-step narrow phase, no `ladrunoBuildNodalMass`** for a tie.
- Large-deformation drift: the frozen weights are valid while relative rotation at the
  interface stays small (the projection handler's frozen-`Ccr` limit, ~0.1 rad). A tie does
  not slide, so this is benign; genuine finite-sliding re-emission is the ADR-60 hook, out
  of scope.

### Decision 2 — Reuse the projection handler verbatim; the tie inherits its guarantees AND its limits

Inherited **for free** (ADR-30, shipped + validated): exactness, momentum conservation,
Δt-neutrality, `system Diagonal` compatibility, the tie-force query
(`ladrunoProjectionTieForce`) and HDF5 recorder, prescribed-motion compatibility.

Inherited **limits** (→ the BLOCKERs): `system Diagonal` only; every tied (slave) DOF needs
nonzero lumped mass; no MP-chains / double-constraints; IC must sit on the constraint
manifold (`-projectICs` to snap); partition-interior only (no MPI cross-partition in v1);
frozen small-rotation `Ccr`.

### Decision 3 — Relationship to the penalty tie family (SOFT / ALM / bipenalty)

`LadrunoTie` (kinematic) is the **default and preferred** explicit mesh-tie. The penalty
family remains for the cases kinematic projection can't serve:
- **SOFT tie** — when the user wants a *compliant* bond (some give) or cannot satisfy the
  projection handler's requirements (e.g. a massless interface node, an MP-chain topology).
- **ALM tie** — **implicit** runs (the projection handler is explicit-only).
- **Bipenalty** — **removed from consideration** (ADR-61 shelved).

Selection: `LadrunoTie` for explicit + disjoint surfaces + massed nodes (the common case);
fall back to SOFT/ALM where its requirements aren't met. The generator **must detect and
report** (named errors) the cases it must hand off — massless slave DOF, chain topology,
non-Diagonal system, implicit integrator — pointing at the SOFT/ALM alternative.

### Decision 4 — Energy / momentum (trivially clean vs the penalty schemes)

A kinematic tie does **no work** (an exact bilateral constraint) → it creates no energy; the
projection is momentum-conserving by construction (ADR-30). So the energy oracle is trivial
compared to the penalty schemes: there is **no penalty strain-energy offset** (unlike the
SOFT/bipenalty RES≈0 problem of ADR-61) and **no fictitious-mass KE to account** (unlike
bipenalty Route B). Total mechanical energy closes directly.

---

## Design-gate BLOCKERs

**BLOCKER-1 (LEAD) — chain-free, single-constraint topology.** The projection handler v1
**refuses** MP-chains (a slave retained elsewhere) and doubly-constrained DOFs. The generator
**must** guarantee: (a) slave and master surfaces are **node-disjoint**; (b) each slave node
is tied to **exactly one** master facet (collocation). Gate: the generator detects a node
that would be both slave and master, or a slave paired to ≥2 facets, and **refuses with a
named error** (hand off to SOFT, or require the user to designate disjoint surfaces). Two
node-disjoint surfaces with one-facet-per-slave-node is the supported, validated topology.

**BLOCKER-2 — every tied slave DOF must carry lumped mass.** The projection keeps slave DOFs
in the equation set (it does NOT eliminate them), so a massless tied DOF is refused by the
handler. Real solid/shell interface nodes carry mass, so this is normally satisfied — but the
generator must **check at emission** and give the recipe-bearing message (this is the
*opposite* of bipenalty's massless-rescue: here mass is a precondition, and it's present).

**BLOCKER-3 — IC on the constraint manifold.** A non-conforming tie generally starts with the
slave surface *not* exactly on the master surface (a small initial gap). The committed initial
`u,v` must satisfy `u_s = Σ N_i u_{m,i}` or the handler aborts. Gate: emit with `-projectICs`
semantics (snap the slave IC onto the facet) OR document that the as-built geometry must be
conforming-at-the-interface. Decide whether the snap perturbs the initial stress state
(it moves the slave node) — for a *built-in* tie the snap is part of defining the bond.

**BLOCKER-4 — shells / rotational ties (ndf 6).** Solid–solid (ndf 3, translational) is clean.
A shell tie must also tie rotational DOFs, which then need a small rotational lumped mass (the
handler's existing `rigidLink -beam` / 3D-diaphragm rule). Deferred to P3; v1 is solid–solid.

---

## Phased implementation plan (oracle-first, mirroring the fork's method)

- **P0 — single-node tie + projection (build-free + tiny regression).** Emit one
  `u_s = Σ N_i u_{m,i}` constraint, route through the projection handler. *Oracles:*
  (a) **exactness** — `u_s − Σ N_i u_{m,i} → 0` to machine precision (numpy + OpenSees);
  (b) **`dt_cr` unchanged** vs the untied model (no penalty spring); (c) **momentum
  conserved** through a transient; (d) **energy closes** (no offset). Reuses the ADR-30 P0
  falsification harness.
- **P1 — full disjoint-surface non-conforming tie (the shipped feature).** Many slave
  nodes, one master facet each, solid–solid, non-matching meshes. *Oracles:* patch test
  (constant stress transmits exactly across a non-conforming tie), a cantilever/bar split by
  a non-conforming tie matches the monolithic mesh, `dt_cr` = the untied bulk, momentum +
  energy clean. Compare penetration vs the SOFT tie (should be **zero**, vs SOFT's `δ/h`).
- **P2 (successor) — integral-mortar (two-sided, dual-basis) ties** needing the handler's
  **deferred MP-chain support**, and **shell rotational ties** (BLOCKER-4).
- **P3 (successor) — finite-sliding re-emission** (the ADR-60 hook) if a tie must survive
  large interface rotation.

---

## Adversarial gate decision

- **P0/P1 — lighter review.** Both halves are already adversarially gated and shipped (ADR-30
  projection handler; ADR-41 mortar pairing). The new code is a thin, setup-time constraint
  generator; the BLOCKERs are topology checks with named refusals. The real risk (chains /
  double-constraints) is caught by BLOCKER-1's gate, not novel math.
- **P2 (integral-mortar + chains) — full gate**, since it requires extending the projection
  handler's MP-chain handling (genuinely new enforcement math).

---

## Ledger / classTag bookkeeping

- `LEDGER_implementations.md` — one row: *LadrunoTie kinematic mesh-tie generator* (emits
  `MP_Constraint`s for the projection handler), **no new class tag**, status per phase, PR.
- `LEDGER_vanilla_files.md` — none expected (uses the public Domain `MP_Constraint` API).
- `LEDGER_quirks.md` — one entry: *a kinematic mesh-tie is emitted as MP-constraints and
  enforced by `LadrunoProjectionHandler`; it inherits that handler's requirements (Diagonal
  SOE, massed tied DOFs, disjoint chain-free surfaces, on-manifold ICs) — the generator must
  refuse-and-hand-off to SOFT/ALM where they aren't met.*
- `classTags.h` — no change. Banner — a new `banner_features.txt` line if shipped distinctly.

---

## Open questions — need sign-off before coding

- **OQ-1 (surface API).** Expose as a standalone `LadrunoTie -slaveSet -masterSurface` command,
  or as `-kinematic` enforcement mode on the existing mortar `-tie`? (I lean standalone — the
  tie is a constraint generator, not contact runtime.)
- **OQ-2 (collocation vs integral-mortar for v1).** Confirm **node-to-segment collocation**
  for P1 (chain-free, projection-handler-compatible, standard `*TIE` behavior), deferring
  integral-mortar to P2. (I lean yes — it's the only chain-free option for the v1 handler.)
- **OQ-3 (IC handling).** For a non-conforming as-built interface, snap the slave IC onto the
  master facet (`-projectICs`), or require the user to build the interface conforming and
  refuse off-manifold ICs? (BLOCKER-3.)
- **OQ-4 (scope).** Ship **P1 (solid–solid, collocation) only**, deferring shells (P3) and
  integral-mortar/chains (P2)? (I lean yes.)
- **OQ-5 (is this the right successor to ADR-61?).** Confirm we pursue the kinematic route
  (this ADR) and leave bipenalty shelved. This is the constructive replacement the LS-DYNA /
  Abaqus / Kratos comparison points to.
