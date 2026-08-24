# LadrunoTie — kinematic mortar mesh-tie via constraint emission (the projection handler)

> ADR-62. Status: **SHIPPED — P1 (collocation, #454) + P2 (integral-mortar, #455) + P3
> (shell/rotational ndf-6, #459) + P2.1 (`-dual` biorthogonal/sparse mortar, #462) + P3.1
> (Hermite w–θ shell edge transfer `-hermite`, #467) all merged; `LadrunoTie` is the
> auto-generator.** Backlog: shell-to-solid (the P4 rung — its own decision record,
> **ADR-64** `64_ladruno_shell_to_solid_tie_adr.md`, planned/awaiting sign-off), and
> mortar-Hermite (P3.1b, needs a kernel per-GP hook). Quadratic (quad8/tri6)
> mortar facets are **ADR-78** (`78_ladruno_quadratic_mortar_tie_adr.md`, accepted
> 2026-08-04 — the apeGmsh ADR 0086 D2 "both" follow-up). See `kinematic_tie_handoff.md`.
> The constructive successor to the shelved ADR-61.
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
distinct constraint-set object needing a tag, the next free contact ELE slot is **33022** — the old 33016 reservation was superseded by `ELE_TAG_LadrunoLST` (ADR-70 P3).)

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

**BLOCKER-4 — shells / rotational ties (ndf 6). SHIPPED in P3.** A shell tie ties the rotational
DOFs (4,5,6) with the SAME per-slave transfer P (`θ_s = Σ P_sk θ_{m,k}`) — collocation weights or
the mortar D⁻¹M row, either mode. The default `-dof` becomes 1..6 for a 3D ndf-6 node
(`-dof`-selectable). Rotational lumped mass: `ShellMITC4`/`ASDShellQ4` `getMass()` NEGLECT rotary
inertia (verified in source), so BLOCKER-2 was made **per-DOF** — it refuses (named) a tied
rotational DOF that carries no mass, and the user adds nodal rotary mass (`mass $n 0 0 0 mrx mry mrz`)
or drops the rotations with `-dof 3 1 2 3`. The handler needs NO change (every (node,dof) is an
independent union-find vertex ⇒ master-only rotational rows). Scope = shell-to-shell (same ndf);
shell-to-solid (ndf mismatch ⇒ needs a `θ=½∇×u` rotation-from-translation coupling, not a straight P)
is guarded with a named refusal and deferred. **Honest limit:** P is linearly complete, so rotations
and CYLINDRICAL bending with the curvature axis ∥ the interface cross EXACTLY (the FE patch test);
a curvature varying ALONG the interface leaves an O(h²) residual on the quadratic transverse
displacement `w` (never on θ) — the Hermite w–θ transfer is deferred as P3.1.

---

## Phased implementation plan (oracle-first, mirroring the fork's method)

- **P0 — single-node tie + projection (build-free) — DONE, all-green.**
  `kinematic_tie_validation/proto_p0_kinematic_tie.py` (a real non-conforming axial bar-tie)
  validates the projection math `a_proj = L(LᵀML)⁻¹LᵀM a_raw`: (a) **exactness** —
  `G·a_proj = −5.6e-17`, projector idempotent to `1e-16`; (b) **`dt_cr` neutral** —
  projection `dt_cr ≥` the unconstrained value with **zero penetration**, while a PENALTY
  enforcing the SAME tie collapses `dt_cr` ~2300× to reach `1e-7` penetration (the thesis,
  quantified); (c) **no-work** — `|f·v_admissible| = 2.5e-16` (tie force ⟂ admissible
  motion → energy-clean); (d) **no fictitious mass** added. *Next:* the OpenSees regression
  (needs a build) reusing the ADR-30 P0 falsification harness.
- **P1 — non-conforming tie on the REAL solver — concept DONE, all-green (no build).**
  `kinematic_tie_validation/proto_p1_kinematic_tie_opensees.py` runs a weighted
  multi-master non-conforming tie (`u_4 = 0.7 u_2 + 0.3 u_3`) as an `equationConstraint`
  on the **shipped** `CentralDifferenceLadruno` + `LadrunoProjection` + `system Diagonal`:
  (a) **EXACT** — tie error `6.7e-18` over 400 steps; (b) **Δt-NEUTRAL** — `dt_cr` tied =
  untied = `0.0316228`, ratio **1.000000** (the kinematic tie adds no stiffness; a stiff
  penalty would collapse it); (c) **load-carrying** — the slave bar moves (load crossed the
  non-conforming interface). The whole ADR thesis validated end-to-end on the real solver
  **with no new code** — because `equationConstraint` already carries the weighted row and
  the projection handler enforces it. *Remaining for the shipped feature:* the `LadrunoTie`
  auto-generator (geometry pairing → emit these constraints) — an ergonomic layer (needs a
  build), not a correctness question. Then the solid–solid patch test + SOFT-penetration
  comparison (should be **zero** vs SOFT's `δ/h`).
- **P2 — integral-mortar ties (`-mortar`) — SHIPPED.** Assemble the global mortar operators
  `D_IJ=∫N_I^s N_J^s dΓ`, `M_IK=∫N_I^s φ_K^m dΓ` over the clipped overlap (reusing
  `LadrunoMortarKernel::integratePair` verbatim), condense ONCE at setup `u_s = P u_m`,
  `P = D⁻¹M`, and emit per-slave `EQ_Constraint`s. **The "needs MP-chain support" premise was
  WRONG/avoidable:** pre-inverting D globally makes every P row tie to MASTER nodes only ⇒ no
  chains ⇒ NO handler change and NO kernel change (the shipped projection handler accepts dense,
  master-only rows). Standard basis ⇒ P dense (one large interface group); the `-dual` biorthogonal
  basis (diagonal D ⇒ sparse P) is the shipped large-interface optimization (P2.1, below). Guards:
  coverage-ratio (self-clip `fullCov`, refuse a slave protruding past the master), conforming-gap,
  post-solve partition-of-unity. Oracle `proto_p2_mortar_tie.py` (13/13) + `tests/test_ladrunoTie_mortar.py` (6/6).
- **P2.1 — dual / biorthogonal basis (`-mortar -dual`). SHIPPED.** ψ=Aᵉ·N per slave facet,
  `Aᵉ=diag(∫N)·(Dᵉ)⁻¹` ⇒ `D_dual` DIAGONAL ⇒ `P=D_dual⁻¹M` LOCAL/sparse (each slave ties only to
  masters under its own facet support ⇒ small local handler groups at large interfaces). Built from
  the SAME `integratePair` Dᵉ/Mᵉ (a per-facet npsS×npsS solve) ⇒ NO kernel change, NO handler change;
  default OFF = standard dense P byte-identical. Preserves linear completeness (row-sum LUMPING D
  would break the patch — the oracle measures 0.28 err lumped vs 2e-15 dual). Oracle
  `proto_p2_1_dual_mortar.py` (12/12: D diagonal, sparse, dual≠lumped) + 4 FE tests in
  `tests/test_ladrunoTie_mortar.py` (dual patch/split/explicit + `-dual`-requires-`-mortar` refusal).
- **P3 — shell / rotational ties (ndf 6). SHIPPED.** The rotational DOFs tie with the SAME P
  (`θ_s = Σ P_sk θ_{m,k}`), both modes; default `-dof` → 1..6 for a 3D ndf-6 node; per-DOF BLOCKER-2
  refuses a massless rotational DOF (shells neglect rotary inertia — add nodal rotary mass or drop
  the rotations with `-dof 3 1 2 3`); handler unchanged; shell-to-solid guarded + deferred (needs
  `θ=½∇×u`). Oracle `proto_p3_rotational_tie.py` (8/8) + `tests/test_ladrunoTie_shell.py` (6/6).
  Linearly complete ⇒ rotations + aligned-cylindrical bending EXACT; along-interface quadratic w
  is O(h²) (Hermite w–θ transfer = P3.1).
- **Successors** — **P3.1** (Hermite w–θ shell transfer for along-interface bending),
  **shell-to-solid** ties, and **finite-sliding re-emission** (the ADR-60 hook) if a tie must
  survive large interface rotation.

---

## Adversarial gate decision

- **P0/P1 — lighter review.** Both halves are already adversarially gated and shipped (ADR-30
  projection handler; ADR-41 mortar pairing). The new code is a thin, setup-time constraint
  generator; the BLOCKERs are topology checks with named refusals. The real risk (chains /
  double-constraints) is caught by BLOCKER-1's gate, not novel math.
- **P2 (integral-mortar) — full gate DONE.** Two adversarial lenses (mortar math + robustness)
  + the handler-topology investigation + the numpy oracle. Finding: the condensation math is
  correct (no bug); the one fix was the coverage guard (key off the per-node full tributary area
  `cover/fullCov`, not the interface average, + a post-solve partition-of-unity backstop) since
  DGESV flags only an exact-zero pivot, not near-singular D. The feared "new enforcement math"
  (handler chain support) turned out unnecessary — global D⁻¹ pre-inversion keeps the topology
  master-only, so the gate's real risk reduced to the setup-time integral's conditioning.
- **P3 (shell / rotational) — targeted verification (not a full loop).** Per the gate-when rule,
  P3 is a gated, well-tested mechanical extension mirroring a proven sibling (the same P, one extra
  per-DOF loop) — so a single focused adversarial pass on the C++ delta sufficed. Finding: the
  per-DOF mass-diagonal index map is exact (cross-checked against `ASDShellQ4::getMass()`, `index=j*6`),
  no false-accept of a massless rotational DOF, no off-by-one/null-deref/iterator bug. Two flags:
  (a) an ndf-6 *beam* slave now defaults to tying rotations too ⇒ false-refuses a pre-P3 lumped-mass
  beam tie (INTENDED default per OQ-P3; escape hatch `-dof 3 1 2 3`, documented); (b) the master-DOF
  existence was unchecked ⇒ added a named shell-to-solid refusal guard (its happy path is covered by
  all shell-to-shell/solid-solid tests; the refusal branch is unreachable via openseespy since ndf is
  model-wide, so it is defense-in-depth). The numpy oracle (`proto_p3`, 8/8) proves the rotational
  transfer + the exact/limit boundary; the FE tests (6/6) confirm the constant-moment patch crosses
  a non-conforming ShellMITC4 tie exactly.
- **P2.1 (dual basis) — oracle + one focused review (novel math, but strong coverage).** The dual
  condensation is new math, so it earned a numpy oracle (`proto_p2_1`, 12/12: D diagonal, linear
  completeness, sparsity, dual≠lumped) AND a targeted adversarial pass on the C++ delta. Finding: NO
  bug — the per-facet accumulation for a multi-facet node is correct, partition of unity holds exactly
  (`Dᵉ·1 = Mᵉ·1 = cᵉ` over the same clipped overlap, since `Σφ=1` per Gauss point), `De.Solve` is safe
  for large Nm (work buffers sized by npsS≤4, not the RHS), the coverage guard + `Ddual≤0` refuse leave
  no partial-facet leak, and the standard path is byte-identical (`P.Zero()` before `D.Solve` is a
  no-op since Solve does `x=b`). Row-sum lumping was explicitly rejected (breaks the patch — oracle T6).

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
