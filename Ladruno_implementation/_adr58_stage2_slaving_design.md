---
title: "ADR-58 P2 Stage 2 — slave-following design + the slaving-mechanism finding"
project: Ladruno
type: design / handoff (pre-implementation)
related:
  - "[[58_ladruno_rigid_body_adr]]"
  - "[[_rigidbody_session_handoff]]"
updated: 2026-06-24
---

# ADR-58 P2 Stage 2 — slave finite-rotation following: design + the crux finding

Pre-implementation design for the last chunk of the rocking MVP. Read after the
[[_rigidbody_session_handoff|handoff]]. **Stage 1 (free-spin SO(3) rotation) is merged**
(#428); the body rotates correctly in its own state, but the **slaves don't track the
rotation yet** and **applied moments are ignored**. Stage 2 closes both + the Housner gate.

## The goal

Slaves must follow the body's FINITE rotation: `x_i(t) = x_c(t) + R(t)·d_i⁰`, i.e.
`u_i = u_c + (R−I)·d_i⁰`, with `R` from the side-channel orientation `q` (element state).
Today the slave MPs use a frozen transport block (R=I), so slaves follow translation only.

## What the framework actually does (verified — the load-bearing constraints)

A source-grounded design pass + direct verification established the mechanics that any
solution must respect:

1. **The Transformation handler reconstructs a slave as the HOMOGENEOUS map
   `δu_slave = T·δu_retained`** — no constant term. `Uc0`/`Ur0` are never read by the
   handler (grepped). So `(R−I)d⁰` can reach a slave only through the **internal
   (retained) node's DOFs** that `T` multiplies — specifically its rotation DOFs (the
   transport/cross-product columns). (`TransformationDOF_Group.cpp:543-563, 890-922`.)

2. **`TRANSF_INCREMENTAL_MP` is defined** (`TransformationDOF_Group.h:45`), so the
   handler does **`incrTrialDisp(T·δu_retained)`** with `δu_retained = u −
   modTrialDispOld` (the retained node's per-step INCREMENT). This is an **incremental
   corotational** slaving: if `T` is rebuilt each step from the current lever arm
   `r=R·d⁰`, the slave increment is `δu_c + δθ × r`, which accumulates to track the
   finite rotation *incrementally-exactly* — provided `δθ` (the internal node's rotation
   DOF increment) equals the body's per-step incremental rotation.

3. **⚠ THE CRUX FINDING (corrects the design pass).** `getT()` calls `getConstraint()`
   FIRST (`:543`), then the handler does `incrTrialDisp(T·δu_retained)` (`:560`). So an
   **exact `setTrialDisp` inside `getConstraint()` does NOT survive** — the handler
   overwrites it. The MP_Joint3D "length correction" is likewise overwritten within the
   step (it works only as a slow cross-step nudge). **Therefore the slave follows
   whatever `T·δu_retained` gives — the incremental transport result — NOT an exact
   placement.** The design pass's "exact-placement escape hatch survives" is WRONG.

4. **The real unsolved knot.** For `δu_slave = δu_c + δθ × r` to track the body, the
   internal node's rotation-DOF **increment** `δθ` must equal the body's per-step
   incremental rotation. But the internal node's rotation DOFs are integrated by the
   **global central-difference solve** (free, unforced ⇒ `δθ ≈ 0`), while the **truth is
   the side-channel `q`** (Stage 1, off the global solve — required by D3). The two are
   disconnected. Driving the global-solve DOFs from `q` each step means (a) fighting the
   solve's own update of those DOFs, and (b) handling the `LogSO3` **π-wrap** (a freely
   spinning body's rotation vector wraps; a difference-of-totals jumps there, corrupting
   one step's `δθ`).

## The reconciliation question (resolve this FIRST, with a run-it probe)

How is the internal node's incremental rotation `δθ` driven from `q` so that
`δu_slave = δu_c + δθ × r` tracks the body, robustly across the π-wrap? Candidates:

- **(C1) Kinematic-mirror via `Node::incrTrialDisp` on the internal node.** OpenSees
  Nodes carry a `rotation` Versor auto-composed from rotational-DOF increments
  (`Node.cpp:842`). Each commit, the element pushes the body's per-step **increment**
  `δθ_body = LogSO3(q_n · q_{n-1}^{-1})` (always small ⇒ no wrap) into the internal
  node's rotation DOFs such that the handler's `δu_retained` rotation part = `δθ_body`.
  Question to verify: can the element write a per-step rotation INCREMENT on the
  retained node that survives into `modTrialDispOld`/`u` differencing, given the global
  solve also writes those DOFs? Probe: drive `δθ` on the internal node, record a slave's
  position vs `x_c + R·d⁰` over a pure spin.
- **(C2) Restrain the internal rotation DOFs + a second "carrier" node** the element
  drives, with the slave MP retained on the carrier. Avoids fighting the global solve,
  at the cost of a second private node per body.
- **(C3) Element-imposed slave increments** — give up the MP for the rotational part and
  have the element `incrTrialDisp` each slave with `δθ_body × r` after the handler runs
  (an `update()`/`commitState` hook). Slaves stay MP-eliminated for translation; the
  rotational increment is added on top. Watch ordering vs. the handler and the toe
  reaction path.

Decide with a **20-line run-it probe per candidate** (pure-spin body, check a slave
tracks `x_c + R·d⁰` in position AND velocity `v_c + ω×r`) before the full build — the
"run it, don't just review it" lesson has caught a silent bug twice this session, and
this is the highest-risk seam.

## Then: moment gather + Housner (the rest, mechanical once slaving works)

- **Moment gather** (`gatherCoMTorque`, currently `≡0`): `m = Σ (R·d_i⁰) × f_i`, with
  `f_i` the slave's current external force read via `Node::getUnbalancedLoad()`
  (translational part) at commit. P3 gate: a static slave point load → exact CoM moment.
- **Housner rocking gate**: rigid block on a compression-only toe — a `zeroLength`
  ENT/gap spring at the bottom-corner slaves (a TEST-HARNESS contact; the contact engine
  is fenced out per D5). Check the analytic free-rocking half-period + post-impact
  ω-ratio `(1 − (3/2)sin²α)`.

## Rejected (with why) — from the design pass, still valid

- **Write a Versor onto the Node so the MP reads it** — nothing reads `Node::rotation`
  except `getTrialRotation()`, which the MP/handler never call (`Node.cpp:709-718`). Dead.
- **Time-varying MP alone, internal rotation DOFs ~0** — `T·0 = 0`; slaves don't rotate.
  Necessary but insufficient; needs the driven `δθ`.
- **Free slaves + element-prescribed disp (no MP)** — slaves become zero-mass free DOFs
  ⇒ singular SOE; SP-restraining them kills the toe contact reaction. Worse on 3 axes.
- **Move rotation into the global-solve internal DOFs (drop the side channel)** — D3:
  the dense/time-varying spatial inertia can't survive `DiagonalSOE` and central
  difference gets no gyroscopic term ⇒ wrong rotational dynamics. The side channel stays.

## State of the worktree

`feat/adr58-rigidbody-p2-stage2` (off latest `ladruno`). A **started, NOT-yet-wired**
`SRC/element/ladrunoRigidBody/LadrunoMP_RigidOffset.h` exists (the time-varying MP skeleton,
base-storage design) — keep or discard depending on the reconciliation choice (C1/C2 keep
a custom MP; C3 may not need one). No Stage-2 code is built or committed.
