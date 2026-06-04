---
title: "ADR 23 — Coupling constraints (RBE2 / RBE3 / linear equation), Abaqus + LS-DYNA aligned"
project: Ladruno
type: ADR / scoping (no code)
status: scoped — no code; decision-capture for a coupling-constraint family
related:
  - "[[20_ladruno_embedded_reinforcement_adr]]"
  - "[[LadrunoEmbeddedRebar_guide]]"
  - "[[ladruno_integrators_guide]]"
  - "[[LEDGER_implementations]]"
tags:
  - adr
  - constraints
  - coupling
  - rbe2
  - rbe3
  - distributing-coupling
  - explicit-dynamics
  - abaqus
  - ls-dyna
updated: 2026-06-04
---

# ADR 23 — Coupling constraints: RBE2 / RBE3 / linear equation (Abaqus + LS-DYNA aligned)

**Status:** scoped, **no code**. This is a decision-capture for a *family* of
coupling constraints OpenSees lacks, grounded in how **Abaqus** and **LS-DYNA**
realize them — and in particular in LS-DYNA's explicit-solver treatment, since the
fork's center of gravity is explicit dynamics. It does **not** commit to a build;
it fixes the design space and the enforcement strategy so that *when* a driver
appears (explicit SSI, load redistribution, moment-transfer embedded members) the
build is a focused extension of the [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]]
kernel rather than a fresh design.

---

## 1. Context & problem

OpenSees exposes only a thin slice of the coupling-constraint space:
`equalDOF` (equality on selected DOFs), `equalDOF_Mixed`, `rigidLink` (a single
slave tied to one master, beam/bar), and `rigidDiaphragm` (planar, in-plane only).
Missing, and routinely wanted:

- **Distributing coupling / RBE3** — a reference node whose motion is the
  *weighted average* of a node set, adding **no stiffness** (load/BC application to
  a control point, spread over a surface). **No OpenSees equivalent at all.**
- **General kinematic coupling / RBE2** — a reference node rigidly driving an
  arbitrary 3D slave set with full 6-DOF kinematics (`rigidLink` does one slave;
  `rigidDiaphragm` is planar).
- **General linear constraint equation** `Σ_k c_k u_{n_k,d_k} = 0` with arbitrary
  coefficients — `equalDOF` only does identity.

The fork already hit the wall that blocks the "obvious" route: a **transformation
(DOF-elimination) MPC densifies the mass** (`TᵀMT` via the Transformation handler)
→ wrong inertia in a `DiagonalSOE` → **explicit-hostile**. That is exactly why
[[20_ladruno_embedded_reinforcement_adr|ADR 20]]'s **Mode T** (perfect-bond by
DOF elimination) was deferred indefinitely, and why the embedded element shipped as
**penalty (Mode P)** with bipenalty for the explicit `dt`.

The strategic question this ADR answers: **is "transformation densifies mass" a law
of nature, or an OpenSees implementation artifact — and what does LS-DYNA do?**

---

## 2. Reference-code alignment

| Capability | Abaqus | LS-DYNA | Nastran | OpenSees today | Gap |
|------------|--------|---------|---------|----------------|-----|
| **Distributing coupling** (ref ↔ weighted avg, *no* stiffening) | `*COUPLING, distributing` / `*DISTRIBUTING` | **`*CONSTRAINED_INTERPOLATION`** | **RBE3** | **— (none)** | **the big one** |
| **Kinematic coupling** (ref → rigid slave set, 6-DOF) | `*COUPLING, kinematic` | **`*CONSTRAINED_NODAL_RIGID_BODY`** (CNRB) | **RBE2** | `rigidLink` (1 slave), `rigidDiaphragm` (planar) | general 3D |
| **General linear equation** `Σ cₖuₖ=0` | `*EQUATION` | `*CONSTRAINED_LINEAR_{GLOBAL,LOCAL}` | MPC | `equalDOF` (identity only) | coefficients + N-node |
| **Embedded member** | `*EMBEDDED ELEMENT` | `*CONSTRAINED_BEAM_IN_SOLID` | — | **`LadrunoEmbeddedRebar` ✓** / `ASDEmbeddedNodeElement` | **DONE** (ADR 20) |
| **Surface tie** (non-matching) | `*TIE` | `*CONTACT_TIED_*` | — | penalty contact (partial) | mortar tie |
| **Joints** | connectors / `*MPC` | `*CONSTRAINED_JOINT_*` | RJOINT/CBUSH | `zeroLength` + joint elems (piecemeal) | catalog |
| **Welds** | — | `*CONSTRAINED_SPOTWELD` / weld | — | — | niche |

### 2.1 How LS-DYNA enforces these in **explicit** (the load-bearing detail)

LS-DYNA implements three generic constraint strategies and chose deliberately
among them — this is the most useful evidence for the fork:

- **Kinematic constraint method** (Theory §26.2): impose the constraint by a
  *transformation of the slave nodal DOFs* (DOF elimination). It **does** preserve
  explicit efficiency — *"the mass is lumped to the extent that only the global
  DOFs of each master node are coupled"* — i.e. LS-DYNA takes deliberate care to
  keep the eliminated system diagonal. **But** it has zoning/penetration problems
  (a finer master than slave kinks the interface), so LS-DYNA **moved away from it
  for tying**.
- **Penalty method** (Theory §26.3): *"the first approach [kinematic] is now used
  for tying interfaces"* was superseded — penalty is *preferred today* because it
  excites little hourglassing, conserves momentum exactly, needs no special
  intersection logic, and — crucially — *"the computed time step size is unaffected
  by the existence of the interfaces"* when the penalty stiffness is chosen
  **~ the interface element stiffness order**.
- **Rigid bodies (CNRB)** (Theory §25): a genuinely different, explicit-native
  treatment — the constrained nodes' masses/inertias are **condensed to the body
  center of mass** (eq 25.5–25.8), the body is integrated as a **6-DOF rigid body**
  (eq 25.1–25.4), and the slave nodes are updated **kinematically** (eq 25.17–25.18).
  The mass stays **diagonal**; there is *no* global `TᵀMT` densification.
- **Joints** (Theory §25.1): penalty, with the penalty stiffness **capped** so it
  does not control the step — `k ≤ (2·TSSF / Δt·Ω_max)²` (eq 25.30). This is the
  LS-DYNA analog of the fork's `-kt auto` (stiffness-cap) and `-bipenalty` (mass-cap).
- **RBE3 (`*CONSTRAINED_INTERPOLATION`)** (Vol I): a single **dependent** node is
  interpolated from a weighted set of **independent** nodes; the dependent node's
  **mass and rotary inertia are redistributed** to the independents; per-independent
  **translational (`TWGHT*`) and rotational (`RWGHT*`) weights**; the dependent
  DOFs (`DDOF`, up to `123456`) **include rotations** — it is LS-DYNA's tool for
  **shell-brick and beam-brick interfaces** (the DOF-mismatch case). In explicit it
  is a **sequential constraint elimination** (order-dependent; a dependent/independent
  node cannot be dependent in another constraint).

### 2.2 The two takeaways

1. **"Transformation densifies mass" is an OpenSees artifact, not a law.** LS-DYNA
   eliminates DOFs (kinematic method, CNRB) *and* keeps the mass diagonal — by
   **dedicated, constraint-specific mass lumping** (condense-to-center for rigid
   bodies; couple-only-master-DOFs for the kinematic method). The generic OpenSees
   Transformation handler does the naive `TᵀMT` and densifies. So a *dedicated*
   fork constraint treatment **could** be explicit-clean — but it is a custom
   element/integrator path, not the stock handler.
2. **For *tying*, LS-DYNA picked penalty** — and the fork already owns a
   best-in-class penalty kernel (anisotropic coupling, augmented-Lagrangian for
   near-exact at moderate stiffness, `-kt auto` stiffness-cap = LS-DYNA's eq 25.30
   idea, `-bipenalty` mass-cap, energy accounting). This is the path of least
   resistance and most reuse.

---

## 3. Goals / non-goals

**Goals.**
- A **distributing coupling / RBE3** capability (the genuine OpenSees gap), usable
  in **both implicit and explicit**, that also serves the **DOF-mismatch interface**
  (frame/shell on a continuum) the way `*CONSTRAINED_INTERPOLATION` does.
- A **general kinematic coupling / RBE2** generalizing `rigidLink` to a 3D slave set.
- A **general linear equation** primitive (`Σ cₖuₖ=0`).
- **Explicit-safe by construction** — reuse the fork's penalty + `-kt auto` +
  bipenalty machinery; no analysis-core surgery.

**Non-goals (this ADR).**
- Surface-to-surface **mortar tie** (`*TIE`) — needs surface pairing/segment
  projection/integration; a separate, larger effort.
- A **joint catalog** (`*CONSTRAINED_JOINT_*`) — OpenSees covers much piecemeal.
- **Welds / tie-break** — niche.
- The **LS-DYNA-style rigid-body condensation integrator** (CNRB done by mass
  condensation) — recorded as the explicit-native alternative, but it is a custom
  integrator path (deferred; penalty is the v1).

---

## 4. Decisions

### D1 — Penalty-first enforcement (explicit-safe; reuse the ADR-20 kernel)
Build the family as **penalty elements** on the [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]]
pattern: a constraint gap `g`, traction `t`, `P = Bᵀt`, `K = BᵀDB`, with the
`-enforce {penalty | al}` family, `-kt auto` (stiffness-cap = LS-DYNA eq 25.30), and
`-bipenalty` (mass-cap) carried over verbatim. This is LS-DYNA's own verdict for
tying (§2.1) and gives explicit capability with **zero analysis-core change**. The
**transformation / DOF-elimination** route stays deferred (the OpenSees handler
densifies mass; the explicit-clean version needs LS-DYNA-style dedicated lumping —
D3 alternative).

### D2 — v1 = distributing coupling / RBE3 (`LadrunoDistributingCoupling`)
The highest-value gap (no OpenSees equivalent) and the DOF-mismatch tool. A
**reference node R** (translations + rotations) tied to a weighted set of
**independent nodes** `{I_i, w_i}`:
- **Translation:** `u_R = (Σ w_i u_i)/(Σ w_i)` (weighted centroid).
- **Rotation:** derived from a **moment-balanced least-squares fit** of the
  independent translations about R (the classic RBE3 derivation) — *not* a free DOF.
  Needs the independent-node **positions** relative to R (geometry, not shape
  functions).
- **No stiffening of the independent set** beyond the penalty reaction, which *is*
  the load-distribution mechanism (a force at R is spread to the `I_i` by `w_i` and
  geometry — work-conjugate to the kinematics). This is the defining RBE3 property.
- **Mass handling:** mirror `*CONSTRAINED_INTERPOLATION` — **redistribute R's mass/
  inertia to the independents** (so there is no spurious massless reference node in
  explicit), *or* fall back to `-bipenalty` on R. Decide per the explicit driver.
- Reuses the embedded element's B-operator (now with **user/area weights** instead
  of shape functions, and a **rotational reference block**), AL, bipenalty, `-host`/
  set-based authoring, energy diagnostics, serialization.

### D3 — Kinematic coupling / RBE2 (`LadrunoKinematicCoupling`)
A reference node R rigidly driving a slave set: `u_i = u_R + θ_R × (x_i − x_R)`.
Two routes:
- **(a) Penalty rigid tie (v1, fork-style):** penalize `g_i = u_i − (u_R + θ_R×r_i)`
  per slave — explicit-safe, reuses the kernel, `-kt auto`/bipenalty apply. Adds a
  (penalty) stiffness to the slaves — the RBE2/rigid behavior.
- **(b) LS-DYNA-style condensation (deferred):** condense slave masses to R, integrate
  R as a 6-DOF body, update slaves kinematically (Theory §25). Explicit-native, no
  added stiffness, but a **custom element + integrator hook** — bigger lift. Record
  as the principled alternative if penalty conditioning proves limiting.

### D4 — General linear equation (two tiers)
- **Tier 1 (quick win):** a `constraint`/`equation` command exposing the **existing
  single-retained `MP_Constraint`'s full `Ccr` matrix** — arbitrary
  `c₁ u_{a,i} + c₂ u_{b,j} = 0` between **two** nodes. The C++ object already
  supports this; only a user-facing command is missing. Covers a large slice of
  `*EQUATION` / `*CONSTRAINED_LINEAR` / `*MPC` for pairs.
- **Tier 2 (deferred):** **N-node** general equation — needs the **multi-retained
  constraint object** the stock `MP_Constraint` can't express (the same wall as
  Mode T). Implicit-domain; the penalty alternative is a degenerate
  `LadrunoDistributingCoupling` with all-equal weights.

### D5 — `dt` preservation reuses both fork knobs (LS-DYNA-validated)
- **`-kt auto`** = the LS-DYNA stiffness-cap (penalty `k ~ host stiffness`, eq 25.30
  / §26.3 "time step unaffected"). Default conditioning.
- **`-bipenalty`** = the complementary mass-cap (Hetherington–Askes), for when a
  light reference/slave node would still collapse the explicit step.
Both already shipped on the embedded element ([[LadrunoEmbeddedRebar_guide]] §7 (bipenalty)); they
drop straight onto the coupling elements.

---

## 5. Use cases — when to reach for each

### 5.1 Distributing coupling / RBE3 (`LadrunoDistributingCoupling`)
The **"apply / transmit without over-stiffening"** tool: a control point drives or
is driven by a node set through a *weighted average*, leaving the set free to deform.
- **Control-point load/BC application** — a concentrated force or moment at a
  reference node spread over a bearing surface or node ring (column base on a
  footing, actuator/jack load on a specimen face, prestress-anchor head) **without**
  the artificial stiffness a rigid tie would inject.
- **Frame / shell on a continuum (the DOF-mismatch interface)** — a beam or shell
  node connected to a solid face, transmitting force but free to rotate relative to
  the bulk: LS-DYNA's `*CONSTRAINED_INTERPOLATION` shell-brick / beam-brick case
  (tunnel lining on soil, embedded plate, base-plate on concrete). *(For a single
  translation-only tie where ASD suffices, see [[20_ladruno_embedded_reinforcement_adr|ADR 20]]
  / the guide; RBE3 is for the weighted-set, moment-carrying version.)*
- **Pile head ↔ cap / reaction distribution** — spread pile-head reactions across a
  cap node set; introduce a column reaction into a mat.
- **Sub-model load introduction (St. Venant)** — apply section resultants at a cut
  face, smoothly distributed, in a zoom/sub-model.
- **Lumped equipment mass/inertia** — attach a point mass (equipment, tank, node
  mass) to a flexible panel via weighted distribution — RBE3 redistributes the
  reference mass to the set (D2), so explicit runs see no spurious massless node.

### 5.2 Kinematic coupling / RBE2 (`LadrunoKinematicCoupling`)
The **"rigid region / rigid driver"** tool: the set moves rigidly with the control
node (adds rigid stiffness — use when that is physically intended).
- **Rigid end zones / panel zones** — a beam–column joint core or rigid offset
  driven rigidly by a control node (a 3D, arbitrary-set generalization of
  `rigidDiaphragm`).
- **Loading platen / rigid boundary** — drive a specimen face uniformly through one
  control node (displacement-controlled test, rigid footing, rigid loading head).
- **Rigid foundation block / pile cap** — tie a node set into one rigid body.
- **Fastener / bolt head** — rigidly tie a bolt-hole node ring to a fastener
  reference node.

### 5.3 General linear equation (`equation` command, D4 Tier 1)
Arbitrary `Σ cₖ u_{nₖ,dₖ} = 0` — the primitive most ad-hoc constraints reduce to.
- **Periodic / cyclic-symmetry BCs** — tie a face's DOFs to the opposite face
  (`u_left = u_right`, or with a coefficient for a phase/offset): RVE / unit-cell
  homogenization, repeating tunnel rings, periodic soil columns.
- **Inclined roller / skew support** — `u_x = tan θ · u_y` on a non-axis-aligned
  boundary.
- **Scaled / lever / gear relations** — `u_a = r · u_b` (mechanism ratio, rigid
  lever arm, displacement amplification).
- **Symmetry on a non-coordinate plane** — enforce the normal-displacement relation
  directly instead of rotating the whole model to an axis.

> [!tip] Picking the right one
> **Load/BC at a point, surface stays flexible → RBE3.** **Region must move as a
> rigid body → RBE2.** **A precise algebraic relation between a few DOFs →
> `equation`.** When unsure between RBE3 and RBE2: RBE3 *distributes* (soft, no
> stiffening); RBE2 *constrains* (rigid, stiffening).

## 6. Implementation sketch (reuse map — no code yet)

~70% is the [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]] kernel:

| Concern | Reused from EmbeddedRebar | New for coupling |
|---------|---------------------------|------------------|
| `g = u_R − Σ w_i u_i`, `P=Bᵀt`, `K=BᵀDB` | ✓ B-operator | weights `w_i` (user/area), not shape fns |
| `-enforce {penalty\|al}` (Uzawa) | ✓ | — |
| `-kt auto`, `-bipenalty`, `getExplicitCriticalTimeStep` | ✓ | — |
| energy/diagnostic responses, serialization, broker | ✓ pattern | new classTags |
| reference node **rotational DOFs** | — | RBE3 moment-balance / RBE2 `θ×r` |
| **mass redistribution** to independents | — | `*CONSTRAINED_INTERPOLATION`-style |
| set-based authoring | `-host`/`-xi` pattern | apeGmsh `g.couple(ref, set, mode, weights)` |

New classTags (Ladruno band ≥33000), broker cases, `stamp_headers.py`, ledger rows,
banner line — per the usual Definition-of-Done.

---

## 7. Risks / open questions

- **RBE3 rotational derivation** — the moment-balanced least-squares fit of the
  reference rotation from independent translations is the subtle part; validate
  against Nastran RBE3 / Abaqus distributing on a canonical load-redistribution
  patch before trusting it.
- **Mass redistribution vs bipenalty** for the explicit reference node — pick one
  per driver; redistribution matches LS-DYNA but changes the independents' inertia.
- **Penalty RBE2 over-stiffening** — a stiff penalty rigid tie adds real stiffness
  to the slaves; the LS-DYNA condensation (D3b) avoids it but is a bigger lift.
- **Sequential vs simultaneous enforcement** — LS-DYNA enforces explicit constraints
  *sequentially* (order-dependent); a penalty family is simultaneous (additive
  stiffness) and avoids that pitfall — a genuine fork advantage to note.
- **No driver yet** — like ADR 20 §10.7 (Nitsche), this stays scoped-only until a
  concrete need (explicit SSI, control-point load application, moment-transfer
  embedded shell/beam) justifies the build.

---

## 8. Build order (when a driver lands)

1. **Tier-1 linear-equation command** (expose `MP_Constraint` `Ccr`) — the cheap
   warm-up, days not weeks, no new element.
2. **`LadrunoDistributingCoupling` / RBE3** — the flagship gap; penalty kernel +
   weighted B-operator + reference rotation + mass redistribution; `-kt auto` /
   `-bipenalty` for explicit. Validate vs Nastran RBE3 / Abaqus distributing.
3. **`LadrunoKinematicCoupling` / RBE2 (penalty)** — sibling, reuses the kernel.
4. *(deferred)* LS-DYNA-style rigid-body **condensation** path (D3b) and the
   **N-node** multi-retained constraint (D4 Tier 2) — only if penalty conditioning
   or implicit exactness demands it.

---

## 9. References

- **LS-DYNA Theory Manual** §25 *Rigid Body Dynamics* (PDF pp. 511–521 — CNRB mass
  condensation eq 25.1–25.18; joints by capped penalty eq 25.19–25.30); §26.1–26.4
  *Contact-Impact* (PDF pp. 521–525 — kinematic constraint vs penalty vs distributed;
  penalty preferred for tying).
- **LS-DYNA Keyword Manual Vol I** `*CONSTRAINED_INTERPOLATION` (PDF pp. 919–924 —
  RBE3, mass redistribution, shell-brick/beam-brick interface, `DDOF`/`TWGHT`/`RWGHT`),
  `*CONSTRAINED_NODAL_RIGID_BODY` (p. 1023), `*CONSTRAINED_LINEAR_{GLOBAL,LOCAL}`
  (pp. 1012–1018), `*CONSTRAINED_BEAM_IN_SOLID` (p. 883).
- **Abaqus**: `*COUPLING` (kinematic / distributing), `*EQUATION`, `*MPC`,
  `*EMBEDDED ELEMENT`, `*TIE`.
- **Nastran**: RBE2 (kinematic), RBE3 (distributing/interpolation), MPC.
- **Fork:** [[20_ladruno_embedded_reinforcement_adr|ADR 20]] (the kernel + Mode-T
  deferral), [[LadrunoEmbeddedRebar_guide]] §7 (bipenalty) (the dt knobs), [[LadrunoEmbeddedRebar_guide]].
