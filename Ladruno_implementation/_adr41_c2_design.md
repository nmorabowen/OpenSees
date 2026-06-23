---
title: ADR-41 C2 design / handoff — frictionless commit-cycle ALM (mortar enforcement)
project: Ladruno
status: in progress — C2.0 #373 + C2.1 #374 shipped; C2.2 (Uzawa-on-commit) is next (see the C2.2 handoff section)
owner: nmora
tags:
  - implementation
  - contact
  - mortar
  - alm
  - handoff
---

# ADR-41 Track C2 — frictionless commit-cycle ALM: wire the mortar kernel into enforcement

> **START HERE (next session).** C1 shipped (#369 + #370 hardening): the header-only
> `LadrunoMortarKernel.h` produces `D_IJ`, `M_IK`, `g̃_I` per slave/master facet pair, oracle
> 30/30, gates structurally robust. **C2 turns that geometry into a SOLVED frictionless contact**:
> the per-slave-node normal multiplier `λ_N`, the penalty→Uzawa enforcement, and the FE wiring that
> assembles `F_c`/`K_c` into the real OpenSees analysis. Friction is C3 (consumes
> `LadrunoFrictionKernel.h`); mesh-tying is C4. C2 is **enforcement only** — no new mechanics kernel.
>
> Read first: this doc, then ADR-41 §How "Mortar integration scheme" + the integration-points list
> (getTangent/c1, Uzawa-rides-commit, active-set, revertToLastCommit invariant), the capstone
> [[48_ladruno_contact_capstone_adr]] C2 row + D1, and the SHIPPED seams cited below.

## Deliverable (C2 scope only)

| File | Role | C2? |
|---|---|---|
| `SRC/domain/contact/LadrunoContactSurface.{h,cpp}` | **extend**: a `SLAVE_SEGMENTS` kind (faceted slave; mortar needs slave–slave `D`, not just a node list). | **✅ C2.0** |
| `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` | **extend**: a `MORTAR` adapter mode (mirror the shipped `SEGMENT` mode) — inline narrow phase calling `LadrunoMortarKernel::integratePair`, assembling `F_c`/`K_c`. | **✅ C2.1/C2.2** |
| `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` | **extend**: per-slave-node `MortarNormalState{λ_N, λ_N_trial}` keyed `(contactTag, slaveSegIndex, localNode)` + `commit`/`revertToLastCommit` Uzawa update + GC (mirror `FrictionState`). | **✅ C2.2** |
| `SRC/analysis/handler/LadrunoContactHandler.cpp` | **extend**: a mortar pairing loop in `handle()` (bucket-sort broad phase → one adapter per slave-facet × candidate master-facet). | **✅ C2.1** |
| parser (`OpenSeesOutputCommands.cpp` + Py/Tcl wrappers) | **extend**: `-slave-segments` surface kind + `mortar` contact selector with `-epsN -augTol -maxAug -ngp`. | **✅ C2.0/C2.2** |
| `contact_prototypes/proto_c2_alm.py` | numpy oracle: assemble `F_c`/`K_c` from C1's `D,M,g̃`; FD-check `K_c` on a rotated config; Uzawa drive `‖g̃‖→tol` at finite `epsN`; release→F=0; eqn-count invariant. | **✅ build first** |
| `LadrunoMortarPair.{h,cpp}`, `LadrunoMortarSegment.{h,cpp}` | **NOT BUILT** — the shipped NTS path keeps narrow phase INLINE in the adapter + state on the Domain (there is no `LadrunoContactPair`). Mortar mirrors that; a separate pair/segment object is ADR-47 refactor scope, not C2. | ❌ |

**No new DOFs, no saddle-point solver, no custom `EquiSolnAlgo`.** This is the ADR-41 Q-DOF/Q-DRIVER
resolution: `λ_N` is **Domain-side per-node state** (Uzawa), augmented **once per `Domain::commit()`**
(the `LadrunoEmbeddedRebar::commitState` precedent), so C2 needs **zero analysis-layer code**.

## The mechanics (frictionless mortar ALM, made concrete)

C1 gives, per slave facet `s` paired with master facet `m`: `D_IJ` (slave–slave), `M_IK`
(slave–master), `g̃_I` (weighted gap), and the per-GP normal `n` already folded into `g̃`. C2 adds
the **normal multiplier** and the **enforcement**:

1. **Nodal multiplier (standard basis).** One `λ_I` per **slave facet node** `I`, work-conjugate to
   the weighted gap `g̃_I` (NOT per-GP — per-GP is only the C3 friction slip `s_e`). The contact
   pressure field is `p(ξ) = Σ_I N_I^s(ξ) λ_I` (standard, not dual ⇒ `D` non-diagonal — ADR-47 owns
   biorthogonal). The **effective nodal pressure** under penalty+Uzawa, compression-only KKT:
   `t_I = min(0, λ_I + epsN · ḡ_I)`,  with `ḡ_I = g̃_I / a_I`, `a_I = Σ_J D_IJ` (the C1 row sum =
   `∫N_I^s dΓ`, so `ḡ_I` is the area-normalized weighted gap, a true length). `epsN` has pressure/length
   units; `-epsN auto` ← the ADR-39 `ladrunoResolveAutoKn` (owning-solid stiffness), reused.
2. **Contact forces** (normal `n` per GP already in C1; for the assembly use the per-facet
   force operator built from `D`,`M`):
   - slave node `K`:  `f_K^s = +n · Σ_I D_KI t_I`
   - master node `L`: `f_L^m = −n · Σ_I M_IL t_I`
   (Equivalently `F_c = Bᵀ t` with `B` the C1 mortar B-operator over slave∪master facet nodes — ADR-41
   §How step 4–5. The **single `n`** form is exact only per-GP; assemble inside the GP loop, normal
   inside, exactly as C1 does for `g̃`.)
3. **Tangent** `K_c = ∂F_c/∂u`. MVP = the **material/penalty** block `epsN · (active rows of D,M)ᵀ(…)`
   — the `∂(D,M,n)/∂u` geometric terms are the C1 linearization STUB (`dD,dM,dn,dξ`), deferred to a
   convergence-refinement sub-phase (analog of ADR-39's deferred `∂n/∂u`). Gate `K_c` vs FD on a
   ROTATED config to 1e-6 with the geometric terms OFF first (small-rotation correctness), then fold
   them when a named gate needs them.
4. **Uzawa augmentation (the ALM win).** Penalty alone leaves residual penetration `g̃ = O(1/epsN)`.
   Uzawa drives `g̃ → 0` at FINITE `epsN` by updating, **once per `Domain::commit()`**:
   `λ_I ← min(0, λ_I + epsN · ḡ_I)` (compression clamp). Across load steps this augments for free on
   stock `Newton`; within a step (when a gate needs it) a held-load **`analyzeAugmented` Tcl/Py proc**
   (zero-increment re-commits reading `‖ḡ‖_∞`) loops to `augTol`.
5. **Active set / KKT.** Node `I` active iff `λ_I + epsN·ḡ_I < 0` (penetration/compression). The set is
   **frozen within an augmentation sweep** (re-solving between augmentations must NOT call
   `domainChanged()` ⇒ no re-number ⇒ eqn count constant — a named gate). Membership changes only
   between physical steps (ADR-39 epoch model).

## The seams it plugs into (SHIPPED — grounded, file:line)

The Explore map (verified against source) — mirror these, do **not** invent a new architecture:

- **Adapter narrow phase is INLINE** in `LadrunoContactFE::getResidual()` (`SEGMENT` branch ~`.cpp`
  L197–241) and `addKtToTang()` (~L265–300); there is **no `LadrunoContactPair`**. Add a `MORTAR`
  mode beside `RIGID_PLANE`/`SEGMENT` (enum ~`.h` L139): a ctor taking slave-facet nodes ∪ master-facet
  nodes + `epsN`/ALM params, a `getResidual` branch calling `LadrunoMortarKernel::integratePair` then
  the `F_c = Bᵀt` loop, and an `addKtToTang` branch for `K_c`.
- **`c1` contract (already solved for NTS).** `getTangent` overrides the base and returns its own
  matrix; fetch `c1` via `IncrementalIntegrator::getCFactor()` and return `c1·K_c`, override
  `addKtToTang` to its own scatter (no double-`c1`); under explicit CDL return a **zero** tangent
  (LHS mass-only). Reuse verbatim.
- **State on the Domain**, keyed for rebuild-survival. `FrictionState` pattern:
  `LadrunoContactDomain::getOrCreateFrictionState(contactTag, slaveTag, segIndex)` (~`.cpp` L136),
  keyed `(c,s,g)`; `commit()` promotes trial→committed (~L174), `revertToLastCommit()` drops trial
  (~L186), GC `frictionGCBegin/Mark/End` (~L145–171) prunes dead slots each `handle()`. Add a parallel
  `MortarNormalState` map + `mortarGC*` + extend the **Domain hooks** (`Domain::commit` ~L2195,
  `revertToLastCommit` ~L2233 already call the contact domain — the integrator-agnostic choke point).
- **Handler factory.** `LadrunoContactHandler::handle()` builds the bucket-sort grid (~`.cpp` L306–332)
  and emits one adapter per (slave,segment) candidate (~L335–396), `numFe++`, `addFE_Element`. Add a
  mortar loop after the NTS loop: same grid, but slave-FACET × candidate master-FACET, emit the MORTAR
  adapter, `mortarGCMark`. **MP_Constraint composition is refused** (handler warns) — make "a mortar
  node may not be an equalDOF/rigidDiaphragm slave" an explicit C2 error (ADR-41 Q-CONSTR).
- **Surface rep gap.** `LadrunoContactSurface` has only `SLAVE_NODES` / `MASTER_SEGMENTS` (enum ~`.h`
  L37). Mortar needs a **faceted slave** for `D` ⇒ add `SLAVE_SEGMENTS` (flat connectivity +
  `nodesPerSeg`, exactly like `MASTER_SEGMENTS`) + the `-slave-segments nps tag…` parser branch.
- **classTags.** `HANDLER 33002` covers the handler; adapters ride `numFe` (no ELE_TAG). Kernels are
  header-only (no tag). C2 needs **no new classTag** (per ADR-41 §classTags; `33001` in the contact
  band stays reserved-but-unbuilt).

## Phasing (explicit-discipline, mirrors the P-series)

- **C2.0 — slave-segment surface + parser.** `SLAVE_SEGMENTS` kind + `-slave-segments`/`mortar`
  parser, Py+Tcl. Gate: a mortar surface round-trips; null-mortar build byte-identical.
- **C2.1 — mortar PENALTY adapter (λ≡0).** The load-bearing first phase (analog of P2a/P2b-1). Inline
  narrow phase + `F_c = Bᵀt`, `t_I = epsN·⟨−ḡ_I⟩`, the penalty `K_c`. Gates: **`K_c` vs FD on a rotated
  config 1e-6** (geometric terms off); static penetration `ḡ = O(1/epsN)`; null path byte-identical;
  explicit CDL zero-tangent. **Oracle `proto_c2_alm.py` first** (assemble F/K from C1 D,M,g̃; FD-check).
- **C2.2 — Uzawa ALM on the commit cycle.** `MortarNormalState` + Domain commit/revert Uzawa update +
  active-set freeze + `analyzeAugmented` proc. Gates (ADR-41 P2): penetration `g̃→` an `epsN`-independent
  tol within `maxAug` (the ALM win penalty can't reach); **release/reopen → exact F=0**; **Hertz**
  sphere/cylinder `p(r)=p₀√(1−(r/a)²)` peak+radius converge under refinement; `E_contact ≥ 0`; **SOE eqn
  count CONSTANT across augmentations**; Kratos `ALMFrictionlessMortarContactCondition` trend cross-check.

## Scope discipline (don't pull C3/C4/ADR-47 work into C2)

- **Frictionless only.** No `λ_T`, no `LadrunoFrictionKernel` call — that is C3 (the kernel is already
  shipped/A1; C3 wires it on the mortar `λ_T` form). No mesh-tying (`-tie`) — C4.
- **Standard basis** (`D` non-diagonal, per-facet `D`-solve for any nodal `λ` recovery). Dual/biorthogonal,
  true-LM saddle point + inf-sup, self-contact, nodal-normal smoothing = **ADR-47**.
- **No new DOFs / no custom algorithm** (Q-DOF/Q-DRIVER resolved = Uzawa-over-penalty on `Domain::commit`).
  A bespoke `LadrunoAugmentedNewton` (tag `33001`, reserved) is built ONLY if a named gate fails
  *specifically because within-step augmentation is required and the held-load proc is insufficient*.
- **Geometric tangent `∂{D,M,n}/∂u` deferred** to a convergence sub-phase (C1 left the stub) — the
  material/penalty `K_c` carries the MVP, exactly as ADR-39 shipped NTS without `∂n/∂u`.

## Risks / hard parts (flagged honestly)

- **The `epsN`/gap units + the `ḡ_I = g̃_I/a_I` normalization** is the silent-error spot (C1's `g̃` is
  weighted, units length·area). Pin it in the oracle FIRST: a constant-pressure block must give
  `t_I == −p̄` and `F_c` == the consistent nodal load; the **patch test from C1 is the guard**.
- **Active-set chatter** across augmentations would re-number (eqn count changes) — freeze the set within
  a sweep (gate it). The ADR-39 epoch model is the precedent.
- **`D` is not diagonal** ⇒ a per-facet `D`-solve if you ever need nodal `λ` from a GP pressure; the
  Uzawa update works directly on `g̃_I`/`λ_I` so the MVP avoids the solve — keep it that way.
- **Adapter rebuild vs state**: `λ_N` MUST live on the Domain keyed by a STABLE id (global slave-facet
  ordinal + local node), like `FrictionState`; an adapter-local `λ` would vanish each `handle()`.
- **`revertToLastCommit` invariant**: `λ_N` is committed-only (mutated solely in `commit`), never on a
  trial — a rejected step's revert is then automatically safe (EmbeddedRebar precedent). Assert it.

## References

- ADR-41 (`41_ladruno_mortar_alm_contact_adr.md`) §How (mortar scheme, integration points, Uzawa-rides-
  commit, active-set, revert invariant, Q-DOF/Q-DRIVER/Q-CONSTR), §Testing (P2 gates + Hertz + Kratos
  cross-check); capstone [[48_ladruno_contact_capstone_adr]] C2 + D1 rows.
- Shipped seams: `LadrunoContactFE.{h,cpp}` (SEGMENT mode + c1 contract), `LadrunoContactDomain.{h,cpp}`
  (FrictionState + commit/revert/GC), `LadrunoContactHandler.cpp` (handle() factory + bucket sort),
  `LadrunoContactSurface.{h,cpp}` (kinds), `LadrunoMortarKernel.h` (C1 `integratePair`), the ALM
  precedent `LadrunoEmbeddedRebar.cpp::commitState`.
- Literature: Popp/Gee/Wall (2010) std/dual mortar ALM; Wriggers *Computational Contact Mechanics*
  (Uzawa/augmented Lagrange); Hertz contact (analytic gate). Skills: `fem-mechanics-expert`,
  `opensees-expert`.
- Workflow: `gh pr create/merge --admin`; oracle runs with **`python3.12`**; watch **Zone-A** before
  merging C++ (C2 touches real TUs, unlike C1) — see [[ladruno-oracle-python]],
  `Ladruno_internal/WORKFLOW_GOTCHAS.md`.

## Oracle-first checklist (do these in `proto_c2_alm.py`, in order)

1. From C1's `D,M,g̃` on a matched + a non-matched pair, assemble `F_c`,`K_c` for a CONSTANT pressure;
   confirm `t_I == −p̄` and the nodal load is the consistent `∫N_I p̄`.
2. **`K_c` vs FD** w.r.t. nodal coords on a ROTATED config to 1e-6 (geometric terms off).
3. **Uzawa**: at finite `epsN`, loop `λ_I ← min(0, λ_I+epsN·ḡ_I)` → `‖ḡ‖_∞ → augTol` in `≤ maxAug`;
   show the converged penetration is `epsN`-INDEPENDENT (the ALM win penalty cannot reach).
4. **Release/reopen → exact F=0** (active set drops, `λ` resets).
5. **Hertz** (1D/axisymmetric reduction) `p(r)` peak + contact radius converge under refinement.
6. Eqn-count invariant across augmentations (active set frozen within a sweep).

Only after 1–6 pass in numpy, wire C2.1 then C2.2 into the C++ seams above, re-confirm via a standalone
`g++` assembly check + Zone-A (C2 is NOT header-only — it links into the build).

---

# C2.2 handoff — Uzawa on the commit cycle (START HERE for C2.2)

> **Status (2026-06-23).** C2.1 shipped (#374): the `MORTAR` `LadrunoContactFE` mode assembles the
> **penalty** force/tangent (`t_I = min(0, epsN·ḡ_I)`, λ≡0) over brute-force facet pairs; the C2.1
> ctor already stores `slaveFacetIndex` + `contactTag` for the λ_N state lookup C2.2 adds. The Uzawa
> mechanics are validated in `proto_c2_alm.py` (18/18: T3 epsN-independent convergence, T4 release).
> C2.2 turns the penalty into true ALM: penetration → tol at FINITE epsN (the headline accuracy win).

## The crux — per-node `λ_N` but per-facet adapters (resolve THIS first)

The augmented nodal pressure is `p_I = λ_I + epsN·ḡ_I` for each **global slave node** `I`, with
`ḡ_I = g̃_I^global / a_I^global` where `g̃_I^global = Σ_facets ∫N_I g_N dΓ` and `a_I^global = Σ_facets ∫N_I dΓ`
(a slave node is shared by several facets). The variational force is a GLOBAL assembly `f = (D^global) p`
— but C2.1 evaluates one **per-facet** adapter at a time, so when an adapter assembles a shared node's
force it does **not** yet know the global gap. This is the one non-obvious decision; the others (commit-
cycle update, active set) follow the shipped `FrictionState` pattern.

**Why summing per-facet forces is correct IF every facet uses the SAME global `p_I`:**
`f_K = Σ_facets Σ_I D_KI^facet p_I n = Σ_I (Σ_facets D_KI^facet) p_I n = Σ_I D_KI^global p_I n` ✓.
So the ONLY requirement is that each adapter reads the **global** `p_I` (same value for a shared node
across its facets). `λ_I` is global and committed (no problem). The trouble is the **global penalty gap**
`ḡ_I^global` inside `p_I` — it needs the cross-facet sum.

**Recommended resolution — lagged global gap on the Domain (one extra Domain map, no new pass):**
1. `LadrunoContactDomain` gains a per-global-slave-node map keyed `(contactTag, slaveNodeTag)`:
   `struct MortarNormalState { double lambdaN; double lambdaNtrial; double gtAccum; double aAccum; double gtGlobal; double aGlobal; }`
   (`lambdaN` ≤ 0, committed-only — mutated solely in `commit()`, so `revertToLastCommit` is automatically
   safe, the EmbeddedRebar invariant).
2. **Adapter `getResidual` (per facet):** for each of its slave nodes `I`, ADD this facet's `g̃_I^facet`
   and `a_I^facet` into `gtAccum`/`aAccum`. Form the per-node pressure with the **previous** iteration's
   global gap: `p_I = min(0, λ_I + epsN · gtGlobal_I / aGlobal_I)` (lagged by one residual eval). Assemble
   `f^s = −(D·p)n`, `f^m = +(Mᵀ·p)n` and the tangent `epsN·B̃ᵀdiag(act/a)B̃⊗(n⊗n)` exactly as C2.1 but with
   the active mask from `p_I` and `a_I` the **global** `aGlobal_I`.
3. **End-of-iteration promotion:** when does `gtAccum → gtGlobal`? Each *residual sweep* over all adapters
   fills `gtAccum`; promote `gtGlobal = gtAccum; aGlobal = aAccum; gtAccum = aAccum = 0` at the START of the
   next sweep. Cleanest trigger: the FIRST adapter to touch a node in a sweep promotes+zeros it (a per-node
   `epoch`/`touched` flag reset by `commit`/`revertToLastCommit`), OR (simpler) promote in `commit()` and
   `revertToLastCommit()` and accept that within-iteration the gap lags — the lag → 0 at convergence (the
   contact gap is stationary at equilibrium; standard Uzawa-penalty practice). **Start with the simplest:
   accumulate during getResidual, read `gtGlobal` from the last committed/promoted sweep.** Validate against
   the oracle's T3 (converged penetration must be epsN-independent — if the lag corrupts it, escalate to a
   per-node `touched`-flag promotion).
4. **Uzawa update in `LadrunoContactDomain::commit()`** (the existing `Domain::commit` hook, beside the
   friction promotion): for every node slot, `λ_I ← min(0, λ_I + epsN · gtGlobal_I / aGlobal_I)`, then
   promote the accumulators. ONE Uzawa step per commit = the EmbeddedRebar precedent; across load steps this
   augments for free on stock `Newton`. Within-step augmentation (when a gate needs it) = the held-load
   `analyzeAugmented` proc (zero-increment re-commits reading `‖ḡ‖_∞` until `< augTol`), a tested Tcl/Py
   recipe — NOT a custom `EquiSolnAlgo` (Q-DRIVER).
5. **Active set frozen within an augmentation sweep** (membership from `λ_I + epsN·ḡ_I < 0`, fixed across
   re-solves between augmentations ⇒ no `domainChanged()` ⇒ **eqn count constant** — a named gate). It
   changes only between physical steps.
6. **GC** (mirror `frictionGC*`): `mortarNormalGCBegin/Mark/End` each `handle()` prunes node slots no live
   mortar pair references.

**Alternative if the lag misbehaves:** a Domain **pre-pass** that, before the residual loop, walks the
mortar pairs once to fill `gtGlobal`/`aGlobal` at the current trial (no lag). Cleaner numerically but needs
a hook to run before FE assembly — heavier. Try the lagged version first.

**Do NOT** use a per-facet-local λ (one λ per facet-node copy): it is variationally inconsistent at shared
nodes and will fail the patch test. λ is per GLOBAL slave node.

## C2.2 deliverables

| File | Change |
|---|---|
| `LadrunoContactDomain.{h,cpp}` | `MortarNormalState` map keyed `(contactTag, slaveNodeTag)` + `getOrCreateMortarNormalState` + the Uzawa update in `commit()` + accumulator promotion + `mortarNormalGC*`. |
| `LadrunoContactFE.{cpp,h}` | MORTAR `getResidual`/`addMortarTang`: read `λ_I` + global `ḡ_I` from the Domain (lazy `theDomain` re-fetch, like friction), accumulate `g̃_I^facet`/`a_I^facet`, use `p_I = min(0, λ_I + epsN·ḡ_I)`. The ctor already carries `theDomain`? — C2.1 passes `contactTag`/`slaveFacetIndex` but NOT `theDomain`; ADD a `Domain*` arg (mirror the SEGMENT ctor) so the adapter can reach the engine. |
| handler | pass `theDomain` to the MORTAR ctor; `mortarNormalGCMark` per built pair. |
| `analyzeAugmented` proc | a tested Tcl/Py recipe (held-load zero-increment re-commits to `augTol`). |
| `proto_c2_alm.py` | already covers the mechanics; add an assembled-multi-facet shared-node check if the lag resolution needs pinning. |

## C2.2 gates (capstone C2 / ADR-41 P2)

- penetration `g̃ → an epsN-INDEPENDENT tol` within `maxAug` (the ALM win penalty can't reach) — the headline;
- **release/reopen → exact F=0** and `λ→0`;
- **Hertz** sphere/cylinder `p(r)=p₀√(1−(r/a)²)` peak+radius converge under refinement (a real elastic mesh —
  the half-space coupling the C2 oracle could not model; this is where the deferred Hertz gate lands);
- **SOE eqn count CONSTANT across augmentations** (no new DOFs; active set frozen within a sweep);
- Kratos `ALMFrictionlessMortarContactCondition` trend cross-check (concept-level).

After C2.2: **C3** = friction (adopt the shipped `LadrunoFrictionKernel` on the mortar `λ_T` form); **C4** =
mesh-tying (`-tie`, zero-gap = active set frozen ON).
