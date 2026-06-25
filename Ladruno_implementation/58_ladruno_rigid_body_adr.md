---
title: "ADR 58 — RigidBody DomainComponent + SO(3) integrator (LS-DYNA / Abaqus aligned)"
project: Ladruno
type: ADR / scoping (no code)
status: P1 in progress — implemented as a zero-stiffness Element (D1 fallback fired at P1 mapping); see Implementation log
related:
  - "[[Ladruno_explicit_roadmap]]"          # §5.5 (this ADR promotes it) + §5.6 joints + Q1 architecture
  - "[[24_ladruno_coupling_constraints_adr]]" # CNRB-as-constraint; the deferred condensation integrator
  - "[[29_ladruno_kinematic_coupling_rbe2_adr]]" # RBE2 penalty rigid tie (the constraint cousin)
  - "[[51_ladruno_element_removal_adr]]"      # the debris "no home for a free rigid body" gap (H6 FALSIFIED seam)
  - "[[48_ladruno_contact_capstone_adr]]"
  - "[[LEDGER_implementations]]"
tags:
  - adr
  - rigid-body
  - domain-component
  - so3-integrator
  - joints
  - explicit-dynamics
  - debris
  - ls-dyna
  - abaqus
updated: 2026-06-24
review: "v2 — strengthened after a 6-lens adversarial review (architecture / SO(3) theory / explicit-mass / citation-fidelity / contact-debris / red-team completeness). Source-verified findings folded into D1–D3, D5; added D7–D9; corrected the slave-map, citations, and P2 gate."
---

# ADR 58 — RigidBody DomainComponent + SO(3) integrator

**Status:** scoped, **no code**. This ADR promotes roadmap item **§5.5** of
[[Ladruno_explicit_roadmap]] to a numbered decision-capture. It fixes the
**architecture** (what kind of object a rigid body is), the **rotational-integration
family**, the **node-slaving / mass topology**, and the **user-facing I/O surfaces
(loads, ICs, recorders)** — so that *when* a driver appears (rocking foundations,
base isolation, multi-body, or AEM debris) the build is a focused extension rather
than a fresh design. It does **not** commit to a build, it scopes **explicit-only
v1**, and it does **not** block the edge-edge contact branch in flight
(`feat/adr57-e1-router`).

> [!info] Four things called "rigid" — keep them apart
> 1. **`rigidLink` / `rigidDiaphragm`** — `MP_Constraint` *factories* that relate
>    DOFs of **existing** nodes. No object owns rotational dynamics.
> 2. **CNRB-as-constraint** ([[24_ladruno_coupling_constraints_adr|ADR 24]] §2.1,
>    [[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29 / RBE2]]) — a reference node
>    *rigidly drives* a slave set; shipped in the fork as a **penalty** tie. Rides
>    on existing nodal DOFs.
> 3. **`RigidBody` object (THIS ADR)** — a DOF-owning object that integrates its
>    **own SO(3) rotational dynamics** and kinematically slaves the nodes of its
>    part. **Closest in-tree precedent: `Joint2D`/`Joint3D`** — an `Element` that
>    creates a **private internal `Node`** + internal `MP_Constraint`s to slave
>    external nodes (`SRC/element/joint/Joint3D.cpp:320-327,466`). RigidBody is
>    *that pattern + finite-rotation SO(3) DOFs on the master + LS-DYNA-style mass
>    condensation* — **not** a thing with "no analogue."
> 4. **rigid SURFACE** — a *faceted* boundary rigidly attached to a body frame and
>    swept to the current config each step. This is the only thing the contact
>    narrow phase can consume (it reads node coordinates via `getCrds()`); a
>    point-mass `RigidBody` has **no faces**. The engine's existing `RigidPlane` is
>    a **static** analytic plane, not a moving rigid master. **A rigid surface is a
>    distinct, unscoped concept** handed off to the AEM/debris ADR (see D5).

---

## 1. Context & problem

OpenSees has no object that *owns 6 DOFs and integrates rotational dynamics*.
`rigidLink` / `rigidDiaphragm` are `MP_Constraint` factories; the fork's RBE2/RBE3
(ADR 28/29) are still **constraints among nodal DOFs**, not a rigid-body *object*.

> [!warning] Honest framing (corrected in v2)
> An earlier draft claimed OpenSees has "no analogue." That is **false** and was
> removed. `Joint2D`/`Joint3D` already create a private internal `Node` and slave a
> node set through internal `MP_Constraint`s; `Subdomain : public Element, public
> Domain` is an Element that owns internal DOFs; `Pressure_Constraint : public
> DomainComponent` is already a *new top-level kind* with a bounded Domain quartet.
> The genuine differentiator for a rigid body is **(i) finite-rotation SO(3) DOFs
> on the master** and **(ii) LS-DYNA-style mass condensation that keeps the
> *translational* system diagonal** — not the existence of a node-owning master.

That leaves a class of problems either un-modelable or modeled by proxy:

- **Foundation rocking / uplift** — today a zero-length rotational spring at a base
  node, which throws away the true tipping kinematics and the *moving* contact
  point. A rigid body + toe contact recovers the actual rocking response (the
  headline use case — and the reason **Housner's analytic rocking block** is a
  named acceptance benchmark in §7).
- **Base isolation** — a rigid superstructure on a deformable foundation separated
  by the isolator; the superstructure should not be thousands of wasted DOFs.
- **Rigid pile caps / footings / anchor & weight blocks.**
- **Multi-body machine-foundation systems.**

And two structural reasons it is *foundational*, not just one more feature:

- **It gates the joint family (§5.6).** Joints are constraints **between** rigid
  bodies — the abstraction has to exist first.
- **It is the missing "home/owner" for contact debris** ([[51_ladruno_element_removal_adr|ADR 51]]).
  But see D5 / §6 — ADR 51's H6 gate *already falsified* the full debris-handoff
  seam; a RigidBody supplies the home/owner, **not** the contact "consumer."

---

## 2. Reference-code alignment

| Capability | LS-DYNA | Abaqus | Nastran | OpenSees today | Gap |
|---|---|---|---|---|---|
| **Rigid-body object** (owns 6 DOFs, own inertia, own rotational integrator) | `*MAT_RIGID` part + CNRB (Theory, Rigid-Body chapter) | `*RIGID BODY` (ref node + `R3D*`/`RB*` elements; pin vs tie nodes; analytical rigid surfaces) | **RBE2/RBAR** (rigid MPC *elements*; **no** dynamic rigid-body object) | `rigidLink`/`rigidDiaphragm`; **`Joint2D/3D`** (internal-node master) | SO(3) master + condensation |
| **Node-set rigidified** | `*CONSTRAINED_NODAL_RIGID_BODY` (CNRB) | `*COUPLING, kinematic` | RBE2 | RBE2 penalty (ADR 29) | general 3D 6-DOF object |
| **Joints between bodies** | `*CONSTRAINED_JOINT_*` | connectors / `*MPC` | RJOINT/CBUSH | `zeroLength` + piecemeal | §5.6 (gated on this) |
| **Free debris body** | eroding contact → free rigid fragment | — | — | **— (none)** | ADR 51 (seam falsified) |

> [!caution] LS-DYNA citation precision (corrected in v2)
> The LS-DYNA Theory Manual is **not** in any in-repo / Seafile reference library,
> so equation/page numbers cannot be verified here. v1 cited `§25, pp.511–521`,
> `eqs 25.1–25.4 / 25.5–25.8 / 25.17–25.18 / 25.30` with false precision (and
> ADR-24 disagrees on the joint range, `25.19–25.30`). **Treat these as
> release-pinned TODOs:** cite a specific LS-DYNA Theory Manual **release**
> (R7/R11/R13 renumber chapters) and re-verify each equation against that PDF at
> P1, or attach page images to the implementation log. Do **not** let the joint
> ADR (§5.6) inherit the penalty-cap formula `k ≤ (2·TSSF/Δt·Ω_max)²` without an
> independent units/derivation check.

### 2.1 The load-bearing reference fact (verified against ADR-24, not invented)

LS-DYNA integrates a rigid body as a **parallel subsystem**, condensing the
slaved nodes' masses/inertias **to the body center of mass** so the global mass
stays manageable — the same insight ADR-24 §2.2 reached: *"transformation-densifies-
mass is an OpenSees artifact, not a law; LS-DYNA keeps it diagonal by dedicated,
constraint-specific lumping."* **Crucial v2 caveat (see D3):** that "diagonal"
lesson is about the **translational** slave masses. It does **not** license a
diagonal *rotational* inertia — a general 3×3 inertia tensor is dense and cannot
pass through `DiagonalSOE`. Joints (§5.6) are penalty constraints *between* bodies,
with the stiffness **capped** so they never control the step.

---

## 3. Goals / non-goals

**Goals.**
- A DOF-owning **`RigidBody`** (CoM translation + finite-rotation orientation), a
  mass + 3×3 inertia tensor, and an **explicit, momentum-conserving SO(3)
  integrator** (default; see D2 for the open scheme choice).
- **Explicit-only v1.** (v1's "works implicitly" goal from the draft is **removed**
  — it contradicted the parked implicit-tangent question; implicit is a recorded
  follow-up. See D2/§6.)
- **Explicit-safe mass topology**: the **translational** CoM DOFs keep a diagonal
  mass (DiagonalSOE-valid); the **rotational** state is integrated in the
  body/principal frame in a **side channel** (D3). Slaved nodes are kinematic
  followers contributing **no free DOFs**, with a **force/moment gather** back to
  the CoM (D3).
- The finite-rotation slave map is the **exact** form
  **`u_i = u_R + (R − I)·r_i`** (the linearized `u_R + θ_R×r_i` was a v1 error — see
  the correction note in D3 / §6).
- A clean seam for the **joint family (§5.6)** and a **non-precluding** hook for
  contact debris (D5).
- The **user-facing surfaces a rigid body actually needs to run**: node-to-body
  assignment (D7), loads/gravity/applied-moment + initial conditions (D8), and
  recorder/output + restart (D9). Without these, P1/P2 cannot even be set up.

**Non-goals (this ADR).**
- The **joint catalog** (`*CONSTRAINED_JOINT_*`) — roadmap §5.6, a separate ADR.
- **Implicit operation** (finite-rotation tangent) — recorded follow-up, not v1.
- **Contact wiring / a moving rigid-SURFACE / rigid-vs-rigid contact / the
  broad-phase rewrite** — explicitly fenced out; belongs to the AEM/debris ADR
  (D5, and ADR-51 `Q-BUCKETSORT-REWRITE`).
- **Deformable→rigid switching mid-analysis** (`*DEFORMABLE_TO_RIGID`); **flexible
  multibody / superelement reduction**; **follower (configuration-dependent)
  loads**; **sensitivity/DDM**.

---

## 4. Decisions

### D1 — `RigidBody` is a standalone `DomainComponent` — argued honestly against the in-tree precedents
The roadmap chose to pay the architectural cost upfront (Q1, lines 648–669). v2
keeps the decision but **grounds it against what actually ships**, because the v1
"no analogue / biggest-risk" framing was false:

- **Precedents that DO exist** (so this is bounded, not open-ended surgery):
  - `Joint2D/Joint3D` — `Element` + **private internal `Node`** + internal
    `MP_Constraint`s (`Joint3D.cpp:320-327,466`). The node-owning-master pattern.
  - `Subdomain : public Element, public Domain` (`Subdomain.h:58`) — an Element
    that owns internal DOFs; element-iterating paths tolerate it today.
  - `Pressure_Constraint : public DomainComponent` (`Pressure_Constraint.h:58`) —
    a **new top-level kind already in `Domain`**, whose entire footprint is a
    bounded quartet (`add/remove/get/getPCs` + a `TaggedObjectStorage` + an iter
    class, `Domain.{h,cpp}`). This is the copy-paste template for "a new kind."
- **The real differentiator** (the honest reason for a dedicated object, not
  "skip-debt"): finite-rotation **SO(3) DOFs on the master** + the **mass
  condensation** that keeps the translational system diagonal. The skip-debt
  argument is a wash (a new kind also needs its *own* recorder/energy/parallel
  paths), so it is **not** the load-bearing reason.
- **"DOF ownership without a Node" is resolved, not a headline risk.** DOFs enter
  the system only through `DOF_Group`s; a node-less `DOF_Group(int tag, int ndof)`
  already exists (Lagrange multipliers), and the **Joint3D private-internal-`Node`
  route** makes the numberer/`DOF_Group`/recorder/IC machinery work **unmodified**.
  **Decision: v1 carries a private internal `Node` for the 6 CoM DOFs.** The
  residual question is narrower — how SO(3) orientation is carried *alongside* that
  node's 3 rotational DOFs (the D2 problem), not whether DOFs can exist.
- **The genuinely new cost** (over `Pressure_Constraint`, which owns no equations):
  DOF/equation participation + a **`FEM_ObjectBroker` construction arm** for
  parallel/database `recvSelf` (serial-only v1 — see §6).

> [!decision] Recorded alternative + falsifiable kill-criterion (carried forward from roadmap Q1)
> Alternative **(b): a zero-stiffness `Element`** (`getResistingForce`/
> `getTangentStiff` → 0, internal SO(3) state) — the Joint pattern taken further.
> **P1 gate (roadmap line 669):** if the DomainComponent scaffolding touches **> ~15
> files in `SRC/` outside `domain/rigid/`**, fall back to (b). D1 is a
> decision-with-fallback, not an assertion.

### D2 — Rotational integration: **explicit, momentum-conserving exp-map** (open scheme), reuse the in-tree SO(3) kernel
v2 demotes this from "strategy fixed" to an **open decision with a default**,
because the two scheme families the draft named together are mutually exclusive:

- **The body-frame EOM** (write it down — the gyroscopic term is the whole
  difficulty): **`I_body·ω̇ + ω × (I_body·ω) = m_body`**. The `ω×Iω` nonlinearity
  is exactly what a naive central-difference split gets wrong.
- **Impossibility note:** no scheme is simultaneously *explicit*,
  *unconditionally stable*, and *energy-AND-momentum conserving*. The D3
  diagonal-mass / explicit mandate therefore **forces the momentum-conserving
  family** and **accepts bounded energy drift** (→0 as Δt→0).
  - **Default (v1): explicit exponential-map / Verlet on SO(3)** (Krysl & Endres
    2005) — conserves linear + angular momentum exactly, energy not conserved,
    conditionally stable.
  - **Energy-momentum (Simo–Wong 1991 / Betsch–Steinmann 2001) is implicit** — a
    per-step nonlinear solve. **Out of scope** for the explicit path; an
    implicit-track follow-up only.
- **Reuse the fork's own SO(3) toolkit** (the single most effort-relevant fact,
  omitted in v1): `SRC/matrix/GroupSO3.h` (`ExpSO3`, `LogSO3`, `TExpSO3`,
  `dExpSO3`, Rodrigues coefficients with small-angle fallbacks) and
  `SRC/matrix/Versor.h` (unit quaternion, `conj_mult`/`normalize`), already used by
  the corotational frame elements (Perez–Filippou 2024). Carry orientation as a
  `Versor`; advance via `ExpSO3`/`TExpSO3`; relative rotations via
  `Versor::conj_mult`. **This lowers the SO(3)-algebra portion of the "L" effort to
  reuse, not new code.**

### D3 — Mass topology: translational diagonal (DiagonalSOE-valid) **vs** rotational dense (body-frame side channel)
This is the **blocker correction** in v2 — the draft's "global mass stays diagonal /
DiagonalSOE stays valid" is **false for a general 3D body**, verified against the
solver:

- A rigid body's rotational sub-block **is** the inertia tensor, **dense** in the
  global frame unless aligned to principal axes — and D2 rotates it every step.
  `DiagonalSOE::addA` reads only `m(i,i)` and either **row-sums** off-diagonals onto
  the diagonal (`lumpDiagonal`) or **drops** them (`DiagonalSOE.cpp:194-205`). So
  `Ixy/Ixz/Iyz` — the gyroscopic/Euler coupling the Dzhanibekov gate exists to test
  — **cannot survive** a `DiagonalSOE`.
- **Decision — split the topology:**
  - The **3 translational CoM DOFs** keep a **diagonal** lumped mass and go through
    the global (Diagonal)SOE. *This* is where the ADR-24 "keep it diagonal" lesson
    applies, and it is true.
  - The **3 rotational DOFs** carry a dense, configuration-dependent inertia and are
    integrated in the **body/principal frame** (constant **diagonal `I_body`**, the
    3 principal moments) via the **exp-map side channel** — *not* fed to the global
    `DiagonalSOE`. `R·I_body·Rᵀ` is used only for spatial-frame output/coupling.
    This is the **"side-channel handler in the explicit step"** the roadmap already
    floated (line 342), now committed.
- **`dt_cr` (reworded — it was crediting condensation for the wrong reason):** the
  body contributes **no element stiffness**, so it **never appears in the CFL
  eigensolve** (`CriticalTimeStep` iterates `getElements()` only). The real `dt_cr`
  coupling enters through **joint/contact penalties** — carry the LS-DYNA penalty
  **cap** into the joint ADR and the D5 hook, else a stiff joint silently collapses
  `dt_cr`.
- **Force/moment gather (the missing dual of slaving):** every external load,
  gravity, damping, and contact/joint reaction on a slaved node must map to the 6
  CoM equations via the transpose of the slaving map — **`Σ f_i`** (force),
  **`Σ r_i × f_i`** (moment). **P3 acceptance gate:** a slaved-node point load
  produces the correct CoM moment.
- **SMS interaction (was absent):** the body is invisible to the element-based SMS
  integrators (`CentralDifferenceSMS`, `LadrunoMassScaling`) and its slaved nodes
  are the existing **slave-node hazard** the mass-scaler already guards. State that
  the body neither benefits from nor interferes with SMS, and its CoM DOFs must not
  gate `dt_target`.

### D4 — The object is canonical; **CNRB / RBE2 stay as the cheap penalty front-end**
`RigidBody` is the source of truth. A "make these nodes rigid" spelling (CNRB)
can be **sugar** that builds a `RigidBody` from a node selector (D7) + computes
CoM/inertia by condensation. The **penalty RBE2 tie (ADR 29)** stays as the
*no-core-surgery* alternative for users who only need a node set to *move* rigidly;
the two coexist (penalty tie = cheap, rides existing nodes; `RigidBody` = exact,
own dynamics, joint-ready).

### D5 — Contact debris: a **non-binding API-shape constraint only** (reconciled with ADR-51)
v2 reduces this from a forward design hook to an API non-preclusion clause, because
**ADR-51 H6 already adversarially *falsified* the debris-handoff seam**
(`51_…:295-308`: *"No home: `LadrunoContactDomain` has no struct/API for a free
rigid body"* → de-scoped to an inert stub or cut):

- A `RigidBody` legitimately supplies the debris **home/owner** (it owns +
  integrates the motion) but **not** the **consumer**: the contact engine consumes
  pre-declared **meshed node-segment surfaces** (it reads `getCrds()` off real
  Domain nodes); a 6-DOF point-mass body has **no faces**, and its slaved nodes
  contribute **no free DOFs** for a contact adapter to bind to. `RigidPlane` is a
  **static** analytic plane, not a moving rigid master.
- **Decision:** the `RigidBody` API must expose only its **body-frame transform
  (CoM + orientation + inertia)** so a *future* rigid-**SURFACE** / debris layer can
  be built on top **without re-opening this object**. **No** contact wiring, **no**
  surface, **no** broad-phase visibility is designed or promised here.
- **Fenced out explicitly:** a moving-rigid-body master narrow phase, rigid-vs-rigid
  contact, and the broad-phase BVH/octree rewrite — all belong to the AEM/pillar-2
  ADR (`ADR-51 Q-BUCKETSORT-REWRITE`).

### D6 — Behind its own CMake flag, off-by-default
Per fork convention (roadmap line ~1022): `RigidBody` compiles behind a dedicated
`Conf.cmake` flag, off by default; flag-off byte-identical.

### D7 — Node-to-body assignment API + over-constrained-node rule (NEW)
The first thing a user touches, absent in the draft. **Decision:** v1 takes an
**explicit node-tag list** (region/part-tag as later sugar). **Rule:** a node already
in another `RigidBody`, an `equalDOF`/RBE2, or a contact `MP_Constraint` is
**rejected** at `setDomain` with a clear error (no silent precedence) — this is the
conflict the constraint handler would otherwise resolve wrongly.

### D8 — Loads + initial conditions on a body (NEW)
No existing load class targets a non-Node/non-Element DOF owner; the P1 "ballistic
under gravity" and P2 "free-spin" gates are otherwise untestable. **Decision (via the
D1 private internal `Node`, which makes these mostly free):**
- **Gravity/self-weight** derived from the condensed CoM mass; **applied
  force/moment** at the CoM (reuse the `NodalLoad` path on the internal node).
- **Initial conditions:** initial CoM velocity **and initial angular velocity**
  (required to set up the Dzhanibekov gate) **and** initial orientation, injected
  through the internal node's `DOF_Group` IC channel.
- **Follower (configuration-dependent) loads: deferred** (non-goal).

### D9 — Recorder / output + restart (NEW)
`Node`/`Element` recorders bind to fixed kinds; the P2 verification *requires*
recording orientation and angular momentum. **Decision:** expose a response menu —
CoM `{disp, vel, accel}`, **orientation as a quaternion** (Euler as a derived
convenience), and **angular velocity / momentum** — readable by `NodeRecorder`
pointed at the private internal node (free via D1), with a thin `RigidBody` response
alias. **Restart:** serialize the orientation quaternion and **renormalize on
`recvSelf`**.

---

## 5. What it unlocks (payoff summary)

- **Modeling capability not correctly doable today**: foundation rocking/uplift with
  the true moving contact point (Housner regime), base isolation (rigid
  superstructure), rigid pile caps / footings / anchor blocks, multi-body machine
  foundations.
- **The joint family (§5.6)** — revolute / spherical / cylindrical / planar /
  universal / translational, with the OpenSees constraint-handler swap.
- **The debris home/owner for AEM** — a separated element can become a free rigid
  body that *owns and integrates* its motion (the *consumer*/contact side remains
  the AEM ADR's problem, D5).
- **Performance** (secondary): stiff parts collapse from a full mesh to 6 DOFs.

Not on the critical path for the edge-edge contact branch — this is the *next
foundation*, not a current blocker.

---

## 6. Resolved decisions & remaining open questions

**Resolved in v2** (were open in the draft): DOF ownership → **private internal
`Node`** (D1); implicit-vs-explicit → **explicit-only v1** (D2/§3); mass topology →
**translational-diagonal + rotational side-channel** (D3); D5 → **API-shape
constraint only**.

> [!question] SO(3) scheme tuning (D2)
> Explicit exp-map/Verlet (Krysl–Endres) is the default, but the exact variant,
> the orientation-storage cadence (`Versor::normalize` frequency), and `LogSO3`
> degeneracy-band handling near π need a prototype + the §7 P2 gate before lock-in.

> [!question] Integration-loop hook
> The explicit integrators apply **one** central-difference formula over the flat
> DOF vector; the body's rotation DOFs must **not** receive that linear update.
> Decide the hook: either the rotation DOFs are **not** global DOFs and the body
> self-integrates outside the `DOF_Group` sweep (an analysis-loop call into bodies
> each step), or the integrator is special-cased. Default leans self-integrating
> (matches D3's side channel). Cost goes into the effort estimate.

> [!question] Parallel (serial-only v1)
> The whole contact subsystem is serial-only (`sendSelf/recvSelf` stubs); ADR-51
> makes OpenSeesMP a *requirement* for the collapse line. So the **debris** use case
> the body motivates is blocked on parallelizing contact first — an independent
> effort. v1 `RigidBody` is **serial-only**; state it, don't hide it.

- **Class-tag/broker:** add a **`RIGIDBODY_TAGS`** block in `classTags.h` (ladruno
  band ≥ 33000) + a `FEM_ObjectBroker` arm; the band is cheap, the broker/parallel
  arm is the real cost. (If D1 reverses to a special Element, the next free
  `ELE_TAG` is **~33015** — the band is taken to 33014, not 33002 as the draft said.)
- **Energy/KE accounting:** the `EnergyBalanceRecorder`/KE path reads mass off
  Nodes/Elements; the body's translational+rotational KE must be added as a touched
  surface (or it reads free via the internal node).
- **Constraint-handler interaction, damping, staged add/remove, sensitivity** —
  enumerated as touched surfaces; damping (Rayleigh on the CoM) matters for rocking
  response and must be in the MVP.

---

## 7. Phasing sketch (explicit-only; MVP = rocking foundation, joints/debris excluded)

> [!note] Effort honesty
> §5.5's single "L" flattens a multi-year program. **P1 alone** is a major
> framework effort (the new kind + internal-node DOF ownership + condensation).
> **The headline rocking-foundation MVP = P1+P2+P3 + D7/D8/D9 + toe-contact** —
> i.e. everything *except* joints (P4) and debris (P5). **P4 and P5 are NOT part of
> this build** (separate ADRs).

| Phase | Scope | Gate |
|---|---|---|
| **P1** (M–L) | `DomainComponent` + private internal `Node` (D1) + node selector (D7) + translational diagonal-mass condensation (D3) + loads/IC/recorder plumbing (D8/D9) | free body under gravity = ballistic; node followers exact; **≤ ~15 files outside `domain/rigid/`** (else fall back to D1(b)); flag-off byte-identical |
| **P2** (M) | body-frame **SO(3)** exp-map integrator (D2, reuse `GroupSO3`/`Versor`) + side-channel rotational solve | **angular momentum** `‖L(t)−L(0)‖/‖L(0)‖ < 1e-10` over N intermediate-axis flips; **energy drift bounded, →0 as Δt→0** (do NOT require energy conservation); torque-free symmetric-top precession rate vs analytic; **Housner** free-rocking half-period + post-impact ω-ratio vs analytic |
| **P3** (M) | finite-rotation slaving `u_i = u_R + (R−I)r_i` + **force/moment gather** to CoM + constraint-handler "derived DOF" marking | meshed-part-as-rigid vs stiff-deformable reference within tol; slaved-node point load → correct CoM moment |
| **P4** (separate ADR) | hand-off seam to the **joint family** (§5.6) | — (not built here) |
| **P5** (separate ADR) | contact-debris: a moving rigid-**SURFACE** + the AEM consumer | — (not built here; D5 only guarantees the object API does not preclude it) |

---

## 8. References

- **LS-DYNA Theory Manual, Rigid-Body Dynamics chapter** — CNRB mass condensation,
  6-DOF integration, kinematic slave update, joint penalty cap. **⚠ release-pin and
  re-verify all equation/page numbers at P1** (manual not in-repo; R7/R11/R13
  renumber). Reconcile the joint-eq range with ADR-24.
- **Krysl & Endres 2005**, *Explicit Newmark/Verlet algorithm for the rotational
  dynamics of rigid bodies*, IJNME — the **explicit exp-map default** (D2).
- **Simo & Wong 1991**, *Unconditionally stable algorithms for rigid-body dynamics
  that exactly preserve energy and momentum*, IJNME 31:19–52; **Simo, Tarnow & Wong
  1992**, CMAME 100:63–116 — the energy-momentum (implicit) family.
- **Betsch & Steinmann 2001** — energy-momentum conserving schemes (implicit;
  out-of-scope track). *Pin the exact paper title/journal at P1.*
- **Housner 1963**, *The behavior of inverted pendulum structures during
  earthquakes* — analytic rocking-block period + restitution (P2 gate).
- **Holzapfel**, *Nonlinear Solid Mechanics*, Ch. 2 — rotation parameterization.
- ~~Simo & Vu-Quoc 1986~~ — **removed**: that is the geometrically-exact *rod* paper
  (CMAME 58), **not** rigid-body dynamics. (Also fix the same mis-citation at
  `Ladruno_explicit_roadmap.md:338`.)
- **In-tree substrate (D2):** `SRC/matrix/GroupSO3.h` (Exp/Log/TExpSO3),
  `SRC/matrix/Versor.h` (unit quaternion). **Precedent (D1):**
  `SRC/element/joint/Joint3D.cpp` (internal-node master), `SRC/domain/subdomain/
  Subdomain.h`, `SRC/domain/constraints/Pressure_Constraint.h` (new-kind quartet).
- Fork siblings: [[24_ladruno_coupling_constraints_adr|ADR 24]] (CNRB-as-constraint,
  the diagonal-mass lesson), [[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]] (RBE2
  penalty tie), [[51_ladruno_element_removal_adr|ADR 51]] (debris home gap, H6
  falsified seam).

---

## Implementation log

### 2026-06-24 — P1 mapping ⇒ **D1 kill-criterion FIRED ⇒ pivot to the Element route (D1b)**
A source-grounded P1 plan (Plan-agent pass, verified against `SRC/`) flipped the
headline D1 decision via its own recorded fallback:

- **File budget:** the D1(a) *DomainComponent* route lands at **~15–16 files**
  outside a new `domain/rigid/` dir — the `Domain` add/remove/get/iter quartet +
  `TaggedObjectStorage` + iter classes in **both** brokers
  (`FEM_ObjectBrokerAllClasses` **and** `TclPackageClassBroker`) + the 4-ctor
  allocation. At/over the **>~15 tripwire**. The D1(b) *Element* route is **6 files**
  (classTags + 3 interpreter + 2 CMake); all object state rides inside the element
  (the `Joint3D` precedent touches **zero** Domain/broker quartet files).
- **The decider (a surprise the scoping ADR did not anticipate):**
  `Domain::commit()` / `revertToLastCommit()` iterate **Nodes and Elements only**
  (`Domain.cpp:2183-2237`). A `DomainComponent` gets **no per-step commit/revert/
  update callback** — so D1(a) would *still* need a custom analysis-loop hook (the
  open §6 "integration-loop hook" question). An **`Element` inherits all three for
  free**, which also makes the P2 SO(3) side-channel integration land naturally in
  `commitState()`/`update()`.

**Decision:** P1 (and the build) is a **zero-stiffness `Element`**:
`getTangentStiff`/`getResistingForce`/`getMass` → 0; a **private internal 6-DOF CoM
`Node`** (Joint3D pattern) carries the condensed diagonal translational mass; slaves
are tied by internal `MP_Constraint`s; SO(3) state self-integrates in `commitState`
(P2). **This supersedes the D1(a) framing above** — D1's recorded alternative (b)
+ kill-criterion is now the chosen path. Other v2 decisions (D3 mass split, D7–D9,
the §7 gates) stand unchanged.

- `classTag` reserved: **`ELE_TAG_LadrunoRigidBody = 33015`** (next free ELE slot
  after 33014; per-registry band).
- Home: **`SRC/element/ladrunoRigidBody/`** (Element tree — moved from the scoped
  `SRC/domain/rigid/`, since it is now an Element, matching the
  `SRC/element/ladrunoKinematicCoupling/` sibling).
- D6 flag note: the fork has **no compile-time-flag precedent** (ADR-39 contact is
  compiled unconditionally, runtime-gated). Honoring "flag-off byte-identical" is
  trivial on the Element route (no edits to shared Domain/broker TUs); a new
  `option(OPS_Use_RigidBody … OFF)` pattern is introduced in `Conf.cmake`.

**P1 increment 1:** ADR pivot recorded; `classTags.h` tag; the `LadrunoRigidBody.h`
Element interface (internal-node + MP-slaving lifecycle + condensed mass).

### 2026-06-24 — P1 COMPLETE (gate PASS)
Full P1 built and verified. `LadrunoRigidBody.{h,cpp}` + `OPS_LadrunoRigidBody.cpp`
+ `CMakeLists.txt` under `SRC/element/ladrunoRigidBody/`; 5 wiring edits (classTags,
element dispatch Tcl+Py, broker, element CMake); test
`Ladruno_scripts/rigidbody_tests/test_p1_ballistic.py`.

- **Design:** zero-stiffness Element; private internal 6-DOF CoM `Node` carries the
  condensed mass (m_body diagonal + body-frame inertia, floored); rigid-link
  `MP_Constraint` per slave (`u_i=u_R+(R−I)d_i`), rectangular 3×6 (3-DOF solid
  slaves) or square 6×6 (6-DOF); slaves zeroed to massless followers (D3); cached
  `Domain*` for teardown (Joint3D pattern).
- **Pre-build adversarial review (3-lens):** cleared the two real worries (Ccr
  transport-block signs, mass double-count) and the latent teardown leak; fixes
  folded in (cached domain + tag-based cleanup, explicit `<string.h>`, doc accuracy).
- **The bug static review can't catch:** `resolveCentroidAndMass()` set `valid=false`
  only on failure and left it untouched on success, so the `if(!valid) return` guard
  fired on the success path too → internal node never created → null `M0` deref →
  segfault on the first analyze step. Found by **running** the gate (exactly the
  review's #2 "run it" finding). Fixed: helper returns a status code, not `valid`.
- **Gate PASS:** ballistic free-fall `z=−½gt²` rel-err **5e-4** (central-difference
  discretization), node-follower abs-err **0.0** (exact, R=I); both Ccr paths.
- **Build gotcha (recorded for next session):** `Ladruno_scripts\build.bat` must be
  invoked via a real cmd/PowerShell (`cmd /c ... | tail` from Git Bash silently
  no-ops — only the cmd banner prints and nothing recompiles). Verify artifact
  mtimes (`dist/bin/opensees.pyd`) after every build before trusting a re-run.

### 2026-06-24 — P2 Stage 1 COMPLETE (Dzhanibekov gate PASS)
The body-frame SO(3) rotational integrator, **free-spin** (no applied moment yet).

- **Topology (resolves the §6 integration-loop question):** the body integrates its
  own rotation in `commitState`, OFF the global solve (D3) — the dense inertia can't
  survive a `DiagonalSOE`, so rotation must be a side channel. The internal node's
  rotation DOFs stay free-but-unforced; the body's orientation lives in element state.
- **Scheme (the key move):** carry **spatial angular momentum `L`** as the integrator
  state, not `ω`. `ω_body = Ibody⁻¹·(qᵀ·L)`; orientation by a **2nd-order midpoint**
  exp-map `q ← q·exp(ω_mid·dt)` (reuses `Versor::from_vector`); `L ← L + m·dt`. For
  free spin `L` is never touched ⇒ `‖L‖` conserved to machine precision *by
  construction*; the Dzhanibekov flip emerges from `ω_body` re-derived as `q` rotates.
  No eigendecomposition (just the constant `Ibody⁻¹`).
- **Surprise found by running it:** a *forward* (step-start) `ω` sample is dissipative
  — energy bled 7% and the free body spiralled to its major axis instead of flipping.
  The **midpoint** sample fixed it (energy drift 7.2e-2 → 7.6e-7).
- **`Versor::rotate` gotcha:** it uses left-scalar mult (`2.0*vec`) which `VectorND<3>`
  lacks; replaced with a local `rotateVec` using only `vec*double` + one-arg `.cross`.
- **D9 recorder:** `orientation` (quaternion), `omega` (spatial), `omegaBody` (the flip
  metric), `angularMom` (spatial `L`). IC via `-omega wx wy wz` (body frame).
- **Gate PASS** (`test_p2_dzhanibekov.py`): `‖L‖` rel-dev 0.0; energy drift 7.6e-7 @
  dt=1e-3, 1.6e-7 @ dt=5e-4 (2nd-order); **7 intermediate-axis flips**.
- **Adversarial review (focused):** core math confirmed correct; fixes folded in —
  `revert`/`revertToStart` now restore the integration clock (B1/N1), a `started` flag
  guards a giant first-`dt` on a time jump (B3), anisotropic slave mass is averaged
  (N4). Stage-1 free-spin means **applied moments are ignored** (gathered in Stage 2);
  made explicit in `commitState`.

**Next (P2 Stage 2):** make the *slaves* track finite rotation (the §6 tension — MP
constraints are linear/homogeneous; the side-channel `q` must drive the slave
kinematics). Then the **force/moment gather** (D3, `gatherCoMTorque`) so an applied
moment / toe-contact reaction drives the spin, and the **Housner rocking** gate.
Reuse `GroupSO3.h`/`Versor.h`; this is the deferred architectural design pass.

### 2026-06-24 — P2 Stage 2 COMPLETE (slave-following + moment gather + Housner gate PASS)
The rocking-foundation MVP closes. Three pieces, each validated by **running an
isolated probe** before the integrated gate (the "run it, don't just review it"
lesson held a third time — see the Housner penalty/lag finding below).

- **Slave-following — mechanism C3 wins (NOT the time-varying MP).** The §6 crux:
  the exact map `u_i = u_R + (R−I)d_i⁰` has a finite, nonlinear `(R−I)d_i⁰` term that
  the Transformation handler's homogeneous `u_c = T·u_r` (and its `TRANSF_INCREMENTAL_MP`
  `slave += T·δu_retained`) cannot carry — and the side-channel `q` is not a retained
  DOF. **Resolution: the element imposes each slave's exact trial displacement via
  `setTrialDisp` in `update()`.** A pre-implementation design pass feared this is
  overwritten by the handler — but that finding was about `setTrialDisp` *inside a
  custom MP's `getConstraint()`* (a different path). The element's `update()` hook runs
  **after** the handler's incremental slave update and **before** residual formation
  (`CentralDifference::newStep` → `applyLoad`/enforceSPs → `Domain::update`→element
  `update()` → `formUnbalance`), so the element wins the last write. **A probe
  (`test_p2_slavetrack.py`) confirmed slaves track `(R−I)d⁰` to 5e-16**, with only a
  one-step explicit half-step lag (the committed slave reflects the orientation it was
  imposed at; `q` advances later in `commitState`). The MPs are KEPT (slaves stay
  derived/non-singular for the solve + the translational force gather). The design
  pass's time-varying-MP skeleton (`LadrunoMP_RigidOffset.h`) was **discarded** — C3
  needs no custom MP.

- **Moment gather** `gatherCoMTorque() = Σ (R·d_i⁰) × f_i` with the **current** spatial
  lever arm (a frozen `d_i⁰` keeps a constant restoring moment and destroys the rocking
  equilibrium at θ=α — verified catastrophic). `f_i = −P` from `getResistingForce()` of
  the external elements incident on the slaves (Newton's third law vs the
  `addResistingForceToNodalReaction` reaction convention). Those elements are **cached
  by tag at `setDomain`** and re-resolved via `Domain::getElement()` in the gather —
  reading forces this way avoids `calculateNodalReactions`, which re-walks the **shared**
  `SingleDomEleIter` that `Domain::commit()` is already iterating (a reentrant
  `reset()` there silently aborts the commit loop — a real trap, found by source
  reading). Tag-resolve (not raw pointers) makes a *removed* incident element safe
  (this fork has runtime element removal). `commitState` reordered so the gather uses
  the **pre-advance** orientation (the config the toe force was evaluated at): gather →
  midpoint exp-map advance → `L += m·dt`. Probe `test_p2_gather.py` (a torsional
  oscillator) confirms the restoring sign, no off-axis leak, clean oscillation.

- **Housner rocking gate (`test_p2_housner.py`) — the headline benchmark.** A slender
  block released from rest at θ0 = α/2, rocking about a toe **pinned by a stiff elastic
  `zeroLength`** (an *element*, so the gather reads its reaction; ADR D5 keeps the real
  contact engine fenced out). The free-rocking quarter-period (release → first upright)
  matches Housner's `(1/p)cosh⁻¹(1/(1−θ0/α))`, `p=√(W·Rdist/I₀)`, to **1.22%**. This
  exercises the whole chain AND its consistency: the parallel-axis `M·Rdist²` part of
  `I₀` emerges from the CoM-translation channel coupling with the `I_cm` side channel.
  **Run-it finding:** an over-stiff pin injects a parasitic restoring (period biased
  fast) through the one-step rotation/translation lag; the bias scales `~k·dt` and
  **converges to the analytic as dt→0** (sweep: k=3e4 gives −5.1%→−2.5%→−1.2% as dt
  halves). Operating point k=3e4 (toe penetration `W/k≈3e-4` — rigid for the physics),
  dt=1e-5. **The penalty-impact restitution `(1−1.5 sin²α)` is intentionally NOT gated**
  — a finite-stiffness toe does not reproduce Housner's instantaneous angular-momentum
  transfer; the quarter-period (impact-free) is the rigorous, achievable gate.

- **No-regression:** with no incident elements `nIncident==0 ⇒ gatherCoMTorque≡0`, so
  the P1 ballistic and P2-S1 Dzhanibekov gates are byte-path unchanged (re-verified).
- **Focused adversarial review:** reentrancy / `f_i=−P` sign / force-vector layout all
  verified correct against the framework; the one block-worthy finding (dangling
  incident `Element*` under runtime element removal) was fixed by caching tags +
  `getElement()` re-resolve. M2 (velocity-dependent incident elements) + 6-DOF-slave
  rotation follow are documented v1 preconditions.

### 2026-06-24 — SHIPPED: Zone-A pytest + banner + M2 (consistent slave velocity)
The two "shipped" follow-ups + the cheapest deferred item landed.
- **Zone-A pytest** `tests/test_ladrunoRigidBody_element.py` (6 cases:
  ballistic / Dzhanibekov / slave-tracking / slave-velocity / gather-sign / Housner;
  Housner at dt=2e-5/±4% for CI speed, the standalone keeps dt=1e-5/±3%). Manifest
  flipped `planned → shipped`; `check_manifest` warns 4 → 3. **Splash-banner** line
  added via `banner_features.txt` + `patch_banner.py`. ([#435](https://github.com/nmorabowen/OpenSees/pull/435))
- **M2 — consistent slave velocity** (the recommended deferred item): `update()` now
  re-imposes `v_i = v_R + ω_spatial × (R·d_i⁰)` alongside the disp, so a
  velocity-dependent incident element (dashpot / viscous contact) and slave velocity
  recorders are coherent with the finite-rotation disp (was the M2 review finding).
  `test_slave_velocity` confirms `v_i = ω×r` to machine precision. **Only slave
  ACCELERATION is now un-imposed** (accel-dependent incident elements remain
  inconsistent — a documented, narrow v1 limit).

**Status: SHIPPED** (manifest `shipped` + banner line). Remaining v1 limits:
serial-only `sendSelf/recvSelf` (stubs), explicit-only (no implicit tangent),
accel-dependent incident elements. **Deferred to follow-up ADRs:** parallel
serialization (gated on parallelizing the contact subsystem anyway) and the
implicit finite-rotation tangent — both large-effort, demand-free; pull when a
concrete multi-body / quasi-static-rigid / debris use case lands.
