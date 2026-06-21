---
title: "ADR-30 P4c — prescribed MASTER (Tier 2) — design"
project: Ladruno
type: design
status: IMPLEMENTED — KINEMATIC imposition (the known-RHS projection below was tried first and ABANDONED)
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"   # §2.4 the spec
  - "[[_adr30_p4b_design]]"                                 # Tier-1 (shipped #323); §2 derived this math
  - "[[_adr30_loop_state]]"
tags: [adr, constraints, explicit-dynamics, projection, prescribed-motion, tier2, kinematic]
---

# ⚠️ DECISION (2026-06-20): KINEMATIC imposition, NOT the known-RHS acceleration projection

The known-RHS *acceleration* projection in §1–§5 below (ADR §2.4's "overwrite a before
projecting") was implemented and then **abandoned after empirical measurement**: a prescribed
master's displacement is externally imposed (GM-integrated for `imposedMotion`, or held constant
for a plain SP) while the slave is leap-frog integrated — they do NOT co-integrate, so projecting
only the *acceleration* leaves the displacement tie `u_c = C u_master` with a **non-converging
O(1e-3) drift** (measured: the error GREW as `dt` diverged from the GroundMotion's internal
`dtInt`, because the two use different integration schemes/steps). Unlike every other constraint
in this handler (machine-precision because both DOFs leap-frog the same projected accel), an
accel-only tie to an imposed master is fundamentally approximate. A constant prescribed master is
worse still: `a_p=0`, so the slave never moves from its IC.

**Shipped approach — KINEMATIC imposition.** A slave driven PURELY by prescribed master(s)
(no free retained master) is itself fully determined: `u_c = Σ_k C_k u_{m_k} + delta`, and
likewise `v_c`/`a_c` (no delta). It is therefore **excluded from the equation set** (`eqn = -1`,
like its masters) and its `u/v/a` are **imposed directly** from the masters each step in
`applyLoad()` (after the masters' own disp is imposed) — exact tracking, zero drift, exactly like
Transformation's kinematic elimination. This is P4b's Tier-1 mechanism extended to slaves. A slave
tied to BOTH a free and a prescribed master (mixed) is refused (the displacement tie cannot be
held by the free-DOF projection). A group left with zero free equations is refused (clean error,
ADR §5.2 — was a segfault). Implementation is **handler-only** (`LadrunoProjectionHandler`):
`classifyDerivedSlaves()` detects + excludes derived slaves and refuses mixed; `applyLoad()`
imposes them; `doneNumberingDOF()` drops prescribed masters + skips derived slaves. The projector
is UNTOUCHED. Empirically: slave-vs-master drift = 0, free-DOF vs Transformation = 0 (machine eps).

Tests `tests/test_adr30_projection_p4c.py` (TC1 imposedMotion→equalDOF master, TC2 constant SP
master, TC3 rigidLink derived-ux + projected-uy, TC4 mixed refused, TC5 prescribed-slave refused,
TC6 zero-free-DOF refused). Quirk found: Transformation+CDL drops a rigidLink uy tie when the
master also carries an imposedMotion SP — projection gets it right, so TC3 self-verifies the
exact ties rather than comparing to Transformation.

---

# (HISTORICAL — abandoned) P4c — prescribed MASTER via known-RHS acceleration projection

**Goal.** Let a prescribed-motion DOF (non-homogeneous SP / `imposedMotion`) be a **master** in
an MP/EQ group and DRIVE its slaves: `a_c = C_f·a_f + C_p·a_p`, with `a_p` the prescribed
acceleration. This is the literal ADR §2.4 "overwrite `a` on those DOFs before projecting the
rest". P4b (shipped #323) handles prescribed DOFs NOT in a group and **refuses** prescribed
masters with a named error; P4c removes that refusal and implements the known-RHS projection.

## 1. The known-RHS projection (math, from P4b §2 — gate-validated SOUND)

A group's DOFs partition into free-retained `f` (eqn≥0, adjustable, has SOE mass), prescribed
masters `p` (eqn=−1, accel known each step = `a_p`, NO equation/SOE mass), and slaves `c`
(eqn≥0) with `a_c = C_f a_f + C_p a_p`. Prescribed masters are NOT rows of `L` (no eqn, no mass)
— they are external knowns contributing `C_p a_p` to the slave rows. Minimizing
`Σ_{i∈f∪c} m_i (a_i − a_raw,i)²` s.t. `a_c = C_f a_f + C_p a_p`:

```
ã_f = a_raw,f                          (free-retained rows)
ã_c = a_raw,c − C_p a_p                (slave rows: subtract the known prescribed contribution)
a_f = (L_f^T M̃ L_f)^{-1} L_f^T M̃ ã    L_f = [I_f ; C_f] (nrows×nRetFree),  M̃ = diag(m over f∪c rows)
a_c = C_f a_f + C_p a_p
a_p : unchanged (prescribed); it carries the support reaction m_p(a_raw,p − a_p)
```

The shipped `project()` is exactly the `nPres=0` (`C_p`=∅) special case. **Momentum** is
conserved on the free directions; the prescribed DOF is an external driver (momentum is not
conserved across a prescribed-motion boundary — its tie force IS the support reaction).

## 2. Where `a_p` comes from — node-read (Q1 RESOLVED)

`a_p` = the prescribed acceleration at the current step. **Source = `node->getTrialAccel(dof)`**
read inside `project()`:
- `imposedMotion` (`ImposedMotionSP`): `applyConstraint(time)` sets the node trial accel from
  the GroundMotion. This runs in `Domain::applyLoad` inside `updateDomain(time,dt)`.
- plain non-homogeneous SP: supplies disp only → node accel stays 0 → `a_p = 0` (the same
  documented limitation as Tier-1; use `imposedMotion` when accel matters).

**Timing (the design-gate Q1, now traced and RESOLVED):** every `project()` call site in CDL is
preceded by an `updateDomain(...)→applyLoad` that sets the prescribed node accel:
- starter `project(*Aprev=a0)`: the `firstStep` block runs `updateDomain(t0,dt)` at
  `CentralDifferenceLadruno.cpp:472` **before** `project(a0)` at :501 → node holds `a_p(t0)`.
- per-step `project(*Aproj)` in `update()`: `newStep` already ran `updateDomain(t_{n+1})` →
  node holds `a_p(t_{n+1})`; `update()` does not re-`setResponse` before `project()`.
So reading the node accel in `project()` is correct at BOTH sites — no SP/GM-direct-read needed
(the gate's Q1 mitigation is unnecessary once the starter `updateDomain` order is traced).

**Projector ↔ node coupling:** the projector currently is Domain-agnostic. P4c gives it a
`Domain*` (set once, e.g. at `buildMass`) + per-group prescribed `(nodeTag,dof)` and resolves
`theDomain->getNode(nodeTag)->getTrialAccel()(dof)` per `project()` (tags not raw Node*, per the
P4b code-gate convention; nPres is tiny so the per-step lookup is negligible).

## 3. Group structure changes (projector)

Add to `Group`:
- `nRet` stays = # FREE retained columns (prescribed masters are NOT counted/rows).
- `L` (nrows×nRet), `m`, `LtML` — over free-retained + slave rows only (as shipped).
- NEW `Matrix Cp` (nCon × nPres) — slave-row coefficients on each prescribed master.
- NEW `ID presNode, presDof` (nPres) — the prescribed masters to read `a_p` from.

`buildMass`: unchanged for `L_f^T M̃ L_f` (already over free cols + slave rows). `project()`:
read `a_p` (nPres), `ã_c = a_raw,c − Cp·a_p`, solve `a_f`, `a_c = (L a_f)_c + Cp·a_p`, scatter
`a_f`/`a_c`, leave prescribed nodes' eqn untouched (they have none). Tie force for free+slave as
shipped; prescribed-master reaction deferred (eqn=−1 ⇒ no Diagonal-SOE mass to form `m_p`).

## 4. Handler changes

- `doneNumberingDOF`: a master vtx with eqn<0 that IS in `prescribedKey` (P4b set) becomes a
  **prescribed column** instead of the P4b refusal: assign it a prescribed-index, accumulate its
  slave-row coefficients into `Cp`, record `(nodeTag,dof)`. Free masters (eqn≥0) unchanged.
- Pass `Domain*` + `Cp` + `presNode/presDof` to `projector->addGroup(...)` (extend the signature).
- Keep the prescribed-**slave** refusal (overconstraint) and the SP-on-slave path.
- `getValue()`/Tier-1 `applyLoad` displacement imposition is UNCHANGED (prescribed master node
  still gets its disp imposed each step).

## 5. IC (checkIC / snapICs)

The slave IC relation now has a prescribed-master disp term: `u_c = C_f u_f + C_p u_p + delta`.
`checkIC` must add `Σ C_p · u_p(committed)` (read the prescribed master's committed node disp);
`snapICs` likewise. Without this, a slave tied to a prescribed master at a nonzero offset would
false-trigger an IC violation. (For `a_p` we read trial accel mid-step; for IC we read committed
disp at `domainChanged`.)

## 6. Scope / tests (`tests/test_adr30_projection_p4c.py`)
- **TC1** `imposedMotion` ground accel on a `rigidLink -bar` MASTER → the slave follows
  `a_c = a_master`; vs Transformation+FullGeneral reference (rel<1e-6).
- **TC2** `equalDOF` master prescribed (constant SP) → slave tracks the master's prescribed disp.
- **TC3** `rigidDiaphragm` master driven by ground motion → corner slaves follow rigid-body
  transport of the prescribed master motion; vs dense reference.
- **TC4** prescribed master + a FREE slave-side spring (mixed) → momentum on free directions,
  support reaction sign.
- **TC5** still-refused: prescribed SLAVE (overconstraint).
- Regression: P4b TP1–TP4 unchanged; full ADR-30 battery + Zone-A.

## 7. Open questions for the gate
- **Q1 (node-read timing)** — RESOLVED above (starter `updateDomain(t0)` precedes `project(a0)`);
  the gate should CONFIRM by tracing, not re-flag.
- **Q2** IC with prescribed-master disp: is reading committed node disp in `checkIC` correct, and
  is the `delta` (already `Uc0 − Σc·Ur0`) double-counting the prescribed master's `Ur0`?
- **Q3** Tie-force/reaction on the prescribed master (eqn=−1): report via node mass, or document
  as unreported? (P4b reports free+slave tie forces only.)
- **Q4** `rigidDiaphragm`/`rigidLink -beam` prescribed master: the transport `C_p` columns must
  carry the rotational/offset terms correctly (general-C, already handled for free masters).
- **Q5** a group with ONLY prescribed masters + slaves (no free retained, nRet=0): is `a_c = C_p
  a_p` directly (no solve)? Handle the `nRet==0 && nPres>0` path.
