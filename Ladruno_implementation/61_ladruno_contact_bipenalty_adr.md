# LadrunoContact — BIPENALTY enforcement (mass+stiffness penalty pair)

> ADR-61. Status: **SHELVED — designed, gated, NOT built (decision taken 2026-06-30).** The
> design is complete and verified sound through 4 adversarial rounds; the P0 kill-gates +
> the P0 sizing oracle then showed the feature is served badly even in its one "unique"
> niche (deformable–deformable tie ⇒ ~100× two-sided interface-mass inflation), so it was
> shelved. The §How design is implementation-ready if a hard rigid-master high-stress case
> ever arrives. See §P0 kill-gate results.
> **Revision v6** — records the P0 sizing-oracle finding (`bipenalty_validation/`) + the
> shelve decision. (v1 sizing unsound → v2; v2 mis-classified rigid-plane → v3 mortar-tie
> only; v3 energy/SMS claims unsound → v4 conservation oracle + refuse SMS; v5 P0 closed
> form; v6 P0 oracle → ~100× wall → SHELVED.)
> Family: ADR-39 (ContactDomain / SOFT) · ADR-41 (mortar/ALM / `-tie`) · ADR-57
> (edge-edge) · ADR-29 (RBE2 bipenalty, the reuse source) · ADR-36/38 (SMS / consistent
> mass + energy registry). Next free ADR slot is 60 (held for finite-sliding
> re-emission) → this is **61**.

---

## P0 kill-gate results (run before any code — the build/shelve decision)

The design below is complete and adversarially verified, but P0 gates whether it is worth
building. Both gates were measured. **Verdict: recommend SHELVE** — sound design, but it
clears its value bar only in a narrow corner.

### Gate 2 — RBE2 redundancy: **PASS, but narrow.**
`LadrunoKinematicCoupling` (RBE2, shipped, with `-bipenalty`) is a **rigid driver**:
`u_i = u_R + θ_R×d_i`, with *explicitly* **no weights, no shape functions, no fit**
(ADR-29 L53/77/152/177). It binds a slave *set* rigidly to one control node and **cannot**
represent a facet-weighted, surface-to-surface, **non-conforming** mortar bond (each slave
node tied to a shape-function-weighted combination of the master facet it projects onto).
RBE3 has weights but is reference→set load distribution, not a surface tie. ⇒ The
facet-weighted explicit mortar-tie **is genuinely unique** — *but only for NON-CONFORMING
interfaces*; for conforming meshes `equalDOF`/RBE2 already serve.

### Gate 1 — SOFT penetration: **SOFT is adequate for the common case.**
SOFT pins the contact mode at `ω·Δt = 2√SOFSCL` by softening `k → k_soft =
SOFSCL·4·m_eff/Δt²`. The residual tie gap under interface force `P`, non-dimensionalized for
a mortar-tie node in a 3D explicit solid run (`Δt = CFL·h/c`, `c² ≈ E/ρ`,
`m_eff ≈ α_m·ρ·h³`, `P ≈ σ·h²`), is a **closed form** — penetration is just the interface
strain × element size × a tuning factor:

```
δ_soft / h = ε_iface · CFL² / (4·α_m·SOFSCL) ,    ε_iface = σ/E
```

| interface stress | SOFSCL=0.1 (default) | SOFSCL→1 (stiffest stable) |
|---|---|---|
| ε=10⁻³ (moderate)          | 0.4 % of h | 0.04 % of h |
| ε=10⁻² (stress-concentrated)| 4 % of h   | 0.4 % of h  |

So SOFT's penetration is **sub-percent of element size** for low–moderate interface stress
or whenever SOFSCL can be pushed toward 1 (linear/benign runs) — almost always within tie
tolerance. Gate 1 *fails* (SOFT is good enough) outside the high-stress corner.

### The deciding tradeoff (why SHELVE)
Bipenalty and SOFT-at-max sit at the **same Courant margin** (`ω·Δt = 2·safety`); they
differ only in which variable absorbs the constraint. Choosing `k_p = β·k_soft,max`:
- **penetration** `δ_bip = δ_soft,max / β`  (β× better), **but**
- **interface mass** `m_p ≈ (β/safety² − 1)·m_s`  (≈ β× inflated).

**Bipenalty buys a β× penetration reduction at the cost of a β× tied-node-inertia inflation,
at equal stability.** To beat SOFT's floor 10×, you inflate the interface inertia ~10× —
which corrupts local wave-transmission / impact dynamics across the tie, arguably a *worse*
trade than SOFT's penetration wherever interface inertia matters.

### Oracle sub-finding — the deformable–deformable niche is even worse (proto_p0)

The build-free P0 sizing oracle (`bipenalty_validation/proto_p0_bipenalty_sizing.py`, run
green) sharpened the tradeoff into a hard wall. The contact mode's reduced mass obeys the
**physical bound `μ ≤ min(m_s, m_m)`**, so a stiff `k_p` between two *light* bodies cannot
be slowed by mass on one side — the lighter side caps `μ`. Consequences measured:

- **Tie to a (near-)rigid/fixed master** (flexible part on a foundation/platen): one-sided
  `m_p` works cleanly — `ω·Δt` lands exactly at `2·safety`; only the flexible side is
  perturbed. **The only clean case.**
- **Deformable–deformable non-conforming tie** (the Gate-2 "unique niche"): slave-only
  sizing returns **`unfixable` (rhs ≤ 0)** the moment `k_p` exceeds SOFT's stable floor.
  The fix needs **two-sided `m_p` on BOTH interfaces, measured at ~80–120× the physical
  nodal mass** (`k_p = 50× k_soft,max`) — a ~100× inertia concentration on the very
  interface explicit dynamics is meant to resolve. **Far worse than SOFT's sub-% penetration.**

So bipenalty's *only* sound application is the rigid-master tie — and even there SOFT's
penetration is already sub-% (the master being heavy/rigid makes SOFT's `m_eff` the slave's).

**Recommendation: SHELVE** (DECISION TAKEN — user shelved after this oracle, 2026-06-30).
The oracle finding removed the last plausible niche: the deformable–deformable tie that
Gate 2 flagged as unique is served terribly (~100× two-sided inflation), and the
rigid-master tie is closely covered by SOFT. The §How design remains complete and
implementation-ready should a hard rigid-master high-stress case ever arrive.

---

## What

Add a **bipenalty** enforcement mode to the Ladruno contact **mesh-tie** (the persistent,
bilateral zero-gap bond — ADR-41 C4 `isTie`): hold a stiff, *accurate* stiffness penalty
`k_p`, and where the real assembled inertia at the tied DOF is **insufficient** to keep the
penalty oscillator under the Courant limit, add a co-located **mass penalty** `m_p` that
raises it to exactly the threshold. The explicit critical step `dt_cr` is then **not
throttled by the penalty** — the price is a bounded *local inertia* perturbation.

**Scope is the tie ONLY** (Round-2 finding: every other contact mode — RIGID_PLANE, NTS
SEGMENT, MORTAR-contact, EDGE_EDGE — is **unilateral**; its `getResidual` gates the force on
`gap<0` and the active set changes every step with no `domainChanged`, which is the hard
problem this ADR defers, §Decision 1). The tie is the one mode with **no gap clamp**
(`LadrunoContactFE.h:114-117`, the full 3-vector `r→0` force, unconditional), so its active
set is constant and `m_p` is constant — the only case where bipenalty is clean. Concretely
the shipped feature is **a facet-weighted (mortar) permanent bond between non-conforming
meshes, made stiff-and-`dt_cr`-neutral under the explicit `CentralDifferenceLadruno`.**

It is the algebraic **dual** of the already-shipped SOFT scheme (ADR-39 B1/B2 + ADR-57
E5), but with a correction the v1 draft got wrong (Round-1 finding A): both SOFT and
bipenalty are governed by the **assembled** gap-mode frequency
`ω² = k·(B M⁻¹ Bᵀ)`, not the isolated-oscillator `√(k/m)`.

| | knob | holds fixed | solves for | sizing quantity |
|---|---|---|---|---|
| **SOFT** (shipped) | `-soft SOFSCL` | the mass | a softer **stiffness** `k_soft` | `k_soft = SOFSCL·4·m_eff/Δt²` |
| **BIPENALTY** (this ADR) | `-bipenalty` | the stiffness `k_p` | added **mass** `m_p` | `m_p = max(0, k_p·(Δt/2safety)² − m_eff,real)` |

where `m_eff = 1/(B M⁻¹ Bᵀ)` is the **same assembled gap-mode mass** `gapModeInvMass()`
already computes for SOFT. **Default OFF ⇒ byte-identical when unused**; explicit-only
(inert under implicit, like SOFT).

The scalar machinery is **not new**: ADR-29's `LadrunoKinematicCoupling` (RBE2) and
ADR-28's RBE3 ship `-bipenalty` with the critical-ratio sizing and a **massless-node
scan** — but note (finding A) that the RBE2 code lumps `m_p` **only on DOFs that are
actually massless** (`LadrunoKinematicCoupling.cpp:422-423`, `if (nodeMassDd>0) continue`),
which is exactly the regime where `√(k/m_p)` *is* the assembled frequency. This ADR
**generalizes** that guard correctly: add `m_p` only by the deficit
`k_p·(Δt/2safety)² − m_eff,real`, reusing `gapModeInvMass` for `m_eff,real`.

---

## Why

A stiff contact penalty `k_p` is needed for accuracy (penetration ≈ load/`k_p`). Under
explicit central difference it collapses `dt_cr ~ 2/√(k_p·B M⁻¹ Bᵀ)`. SOFT fixes this by
*softening* `k_p` to whatever the live `Δt` allows — correct for dynamics, but it **caps
accuracy** (the penetration floor is set by the soft stiffness, not the user's tolerance).
Bipenalty instead keeps `k_p` stiff and adds the minimum mass that makes the penalty mode
Courant-stable — accurate + `dt_cr`-neutral, paying a bounded local mass instead of
penetration.

**But the value is unproven and the niche is narrow** (Round-1/2 findings D). Three
existing facilities already cover most of the tie space:
- **SOFT tie** — `k_soft` is already large for small `Δt`, so a SOFT mortar-tie may
  *already* be tight enough; the accuracy delta is **unmeasured**.
- **ALM `-tie`** (ADR-41) — gives the **exact** bond (`r→0`) but on the **implicit** lane.
- **RBE2 `-bipenalty`** (ADR-29, shipped) — already a stiff, `dt_cr`-neutral, explicit-safe
  permanent tie, but **node-to-node** (a rigid point driver `u_i=u_R+θ_R×d_i`): no facet
  weighting, no `∫N_I dΓ` mortar projection.

So bipenalty's **only surviving unique niche** is the intersection of all three gaps: a
**facet-weighted (mortar) tie between NON-CONFORMING meshes, under an EXPLICIT integrator,
held stiff**. If the user's tie is node-conforming, the answer is **use RBE2 `-bipenalty`;
do not build this.** ⇒ This ADR is gated by **two P0 kill-gates** (below): (1) SOFT
mortar-tie penetration exceeds tolerance at a realistic `Δt`; (2) the tie is genuinely
facet-weighted / non-conforming (RBE2 cannot express it). If either fails, **stop.**

---

## Where

### New code — NONE new beyond existing files (the SOFT pattern)

Bipenalty is a **mode/flag on existing classes**, **no new class tag** (verified against
`SRC/classTags.h`; SOFT added none either). Touched:

- `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — a `bpScale`/`bpDt` member + an
  `mpNormal()/mpTangent()` sizing pair + a **new decomposed `gapModeMassTerms()` helper**
  (returns `{invMproj_slave, s_master}` — the fused `gapModeInvMass` does NOT expose the
  two terms the deficit formula needs, Round-2 finding A2).
- `SRC/analysis/handler/LadrunoContactHandler.{h,cpp}` — parse `-bipenalty -dtcr Δt
  [-safety s]`; **build the Route-B nodal-mass machinery the handler does NOT currently
  have** (Round-2 finding): a `bipenaltyInjected` map, a `removeBipenaltyMass()` called at
  the **top of every `handle()`** (re-baseline — without it `m_p` COMPOUNDS each
  domainChanged), and a **non-empty destructor** that subtracts against the committed
  Domain (without it `m_p` LEAKS to later stages). `ladrunoBuildNodalMass` is only a *read*
  cache rebuild — none of this is inherited from it or from SMS.
- Parser: `constraints LadrunoContact -bipenalty …` (Tcl + Python), mirroring `-soft`.
- **Reuse** `LadrunoEmbeddedKernel` (`criticalTimeStep`, `maxAbsDiagonal`) only as helper
  arithmetic; the **sizing inequality uses the new decomposed gap-mode helper**, not the
  bare `massPenaltyDtcr` (which is the massless-limit form — Round-1 finding A).

### Modify vanilla

- **NONE expected.** Route B writes `Node` mass via the public `Node::setMass`/`addMass`
  API that `CentralDifferenceSMS` already uses; the commit/restore lives in the handler.
  Any upstream touch → `LEDGER_vanilla_files.md` + `// Ladruno` comment.

### classTags

**None.** Verified: bipenalty is a flag on the shipped `LadrunoContactFE` +
`LadrunoContactHandler (33002)`. *Note (finding D6):* if a contact element ever needs a
broker tag, the next free ELE slot is **33022** — the old 33015/33016 reservations are SUPERSEDED (33015 = LadrunoRigidBody ADR-58, 33016 = LadrunoLST ADR-70 P3, 33017 = LadrunoUP ADR-71, 33018 = LadrunoBrick20, 33019 pencilled H27, 33020 SolidShell, 33021 CSTPair);
the historical ADR-39 "33015" note is dead. This ADR allocates nothing.

---

## How

### What sets `dt_cr` (the corrected scheme)

Central difference is stable when `Δt ≤ 2/ω_max(M⁻¹K)`, stability factor **1.0**, no
Noh-Bathe bonus (`CentralDifferenceLadruno.h:90-93`). A contact penalty `k` contributes a
mode with eigenvalue **`ω²_contact = k · (B M⁻¹ Bᵀ)`**, with the **real assembled masses**
in `M⁻¹`. Worked 2-DOF (slave mass `m_s` + lumped `m_p`, master `m_m`, `B=[n|−n]`):

```
ω²_contact = k_p·( 1/(m_s+m_p) + 1/m_m )       (rigid/fixed master ⇒ 1/m_m → 0)
```

This equals `k_p/m_p` **only** when `m_s=0` and `m_m→∞` — the massless/rigid limit. For
finite masses, `√(k_p/m_p)` *under-estimates* the real frequency, so the v1 sizing
under-stabilized (18–58% in worked cases). **Corrected sizing** — pick the smallest `m_p`
that makes the assembled mode Courant-stable, `ω_contact·Δt ≤ 2·safety`:

```
let s_target = (2·safety/Δt)² / k_p          # upper bound on B M⁻¹ Bᵀ  [1/mass]
                                             #   1/s_target = k_p·(Δt/2safety)²  [mass]
                                             #   (the §What table writes this MASS form;
                                             #    reciprocal of s_target — same quantity)
let m_s_proj = 1 / invMproj_slave(n)         # gap-NORMAL-projected slave mass
let s_master = Σ N_i²·invMproj_master_i(n)   # master-only term of B M⁻¹ Bᵀ
# (BOTH from a NEW decomposed helper — the fused gapModeInvMass returns only
#  s0 = 1/m_s_proj + s_master, which does NOT expose the two terms — see note)
let rhs = s_target − s_master                # capacity left after the master's own mass
if rhs ≤ 0:   m_p = +∞ unattainable  →  the MASTER body alone throttles dt;
              FALL BACK to SOFT (warn) — bipenalty cannot fix a too-light master.
              (Cannot occur for the tie when one side is fixed; only the deferred
               two-deformable-body case.)
else:         m_p = max(0,  1/rhs − m_s_proj)
```

**Decomposed-helper note (Round-2 finding A2).** `gapModeInvMass(n,N)`
(`LadrunoContactFE.cpp:426`) returns the **single fused scalar**
`1/m_s_proj + Σ N_i²·invMproj_i` and does NOT expose the slave/master terms separately,
which the deficit formula needs. The implementation adds a small **decomposed sibling**
(e.g. `gapModeMassTerms() → {invMproj_slave, s_master}`) refactored from the
already-separable `ladrunoNodeMass`/`ladrunoInvMassProj` blocks. So §Where/§How say
"decomposed sibling," NOT "reuse `gapModeInvMass`." All terms are `ladrunoInvMassProj(·,n)`
projections on the active gap normal (raw `Node::getMass` would be wrong off-axis).

Two consequences the v1 draft missed, both load-bearing:
1. **`m_p` depends on the real `m_eff`** — it is the *deficit* `1/rhs − m_s`, **not**
   independent of it. When the node is already heavy enough (`s0 ≤ s_target`), `m_p = 0`
   — bipenalty adds nothing. This is the RBE2 massless-guard generalized correctly: add
   mass only where the real inertia is *insufficient*, not merely zero.
2. **Two light deformable bodies can be unfixable** (`rhs ≤ 0`): the master's own
   `Σ N_i²/m_m,i` already exceeds `s_target`, so no slave-side `m_p` neutralizes `dt`.
   This is a hard reason the deformable-deformable case is **out of P1 scope** — SOFT,
   which sizes `k` from the full `m_eff`, is the robust choice there.

The tangential (stick) mode is sized the same way with `B_t` and `m_eff,t` (the worst-case
of the two basis tangents to `n`), reusing the `softKt` worst-case-tangent logic — but to
solve for `m_p,t` from the **assembled** `B_t M⁻¹ B_tᵀ`, not from `k_t` alone (finding A4).

### Cadence — `m_p` is set once per `domainChanged`, NOT per step (finding B)

The explicit `M` is assembled and factored **once** at the first step after each
`domainChanged` (`TransientIntegrator::formTangent` inside `CentralDifferenceLadruno::
newStep`'s `firstStep` guard); `LadrunoContactHandler::handle()` (which builds the FE
adapters + `ladrunoBuildNodalMass`) runs **only** at `domainChanged` too. The contact
narrow phase (active/inactive) is re-evaluated every `getResidual` (each step), but **the
inverted mass is frozen between domainChanges**. Therefore:

- A **mass** penalty cannot be sized "per step from the live `Δt`" the way SOFT's
  *stiffness* is (SOFT is applied in the residual, which *is* per-step). `m_p` is
  **domainChanged-cadence**. Size it once, from the `Δt` that will be used, and **guard**:
  if `newStep` later sees a `Δt` different from the sizing `Δt` (mid-run change without a
  domainChanged), warn/refuse — `ω_p·Δt` would silently exceed `2`. This guard is a P0
  requirement.
- **Neither Route A nor Route B gives finer-than-domainChanged cadence** (finding B4) —
  Route B's hoped-for "inject more often" is refuted. The route choice is purely about
  **energy bookkeeping + Node-state hygiene**, not active-set responsiveness.

### Decision 1 — Scope: the MESH-TIE only (Round-2 redraw)

**Ship `P0 (two kill-gates) → P1 (mortar mesh-tie) only.** Defer **everything unilateral**
— RIGID_PLANE, NTS SEGMENT, MORTAR-contact, EDGE_EDGE — to a successor ADR.

The dividing line is the **gap clamp, not "one fixed side"** (Round-2 finding 1/4 — v2's
"fixed master ⇒ constant active set" was wrong). Verified in code: the tie path
(`LadrunoContactFE.h:114-117`) assembles the full 3-vector `r→0` force **unconditionally**;
every other mode gates the force on `gap<0` (`getResidual` RIGID_PLANE branch
`if (g < 0.0)`, `LadrunoContactFE.cpp:584`; SEGMENT/MORTAR/edge use the same `⟨−gap⟩₊` /
`min(0,·)` clamp). **A rigid plane is therefore UNILATERAL** — a slave dropped onto a rigid
floor closes the pair at impact (BLOCKER-1) and the active set flips inside `getResidual`
with no `domainChanged` (BLOCKER-2), exactly like a two-body contact. The fixed master only
removes the `rhs≤0` failure; it does NOT make the active set constant. So **only the tie**
(no clamp ⇒ constant active set ⇒ constant `m_p` ⇒ no activation event) is clean.

This makes the **build-at-all decision ride on P0's RBE2-redundancy gate** (Round-2
finding 2): with scope = mortar-tie, if RBE2 `-bipenalty` (node-to-node) already covers the
user's tie, there is nothing unique to ship. The surviving niche is exactly the
**facet-weighted (mortar) tie on non-conforming meshes under explicit CD**.

**The unilateral problem (why everything else is deferred — but NOT "intractable").** A
co-located `m_p` that appears/disappears with the active set causes spurious impact energy
(`½ m_p v²`, a discontinuous inertia change inconsistent with the half-step velocity
history — *not* a clean `Δp=m_p·v` momentum kick, finding A2) and a stale inverted `M`
(active set flips with no `domainChanged`; forcing one re-seeds the whole leap-frog
starter). **There is a known route in, recorded for the successor ADR (Round-2 finding 3):
a PERMANENT, slave-surface-only mass floor** — present whether or not the pair is closed,
so there is no activation event, `M⁻¹` never goes stale, and the `½m_p v²` discontinuity
vanishes (this is standard LS-DYNA contact-added-mass practice, and is energy-accounted
exactly like lumped SMS). The deferral is a **sequencing** choice (the floor still needs
converged-penetration + friction-drift validation), **not** "appearing mass is unsolvable."
v2's blanket "research problem / reject the floor" framing was too strong; the *event-free
slave-surface* floor (not v2's rejected "floor on every potential-contact node") is the
design.

### Decision 2 — Where the mass penalty lives: **Route B (nodal injection) everywhere** (resolves OQ-1)

Three candidate routes were considered; **Route B wins** on the energy analysis (findings
C1/C2/C5):

- **Route A — contact-FE `addMtoTang`** (currently a deliberate no-op so contact adds zero
  mass). *Mechanically assemblable* (the FE is in the `formTangent` element loop), but the
  injected `m_p` is **invisible to the entire EnergyBalanceRecorder** (which sweeps Domain
  *elements* and *nodes* only — the contact FE is an `FE_Element` *adapter*, not a Domain
  Element, so the recorder never calls anything on it). Route A would therefore need its
  own Olovsson-style energy-registry **publisher built from scratch**
  (mirroring `LadrunoMassScalingEnergy`), and even then the single-`owner` registry can't
  host consistent-SMS + Route-A bipenalty simultaneously (finding C5). Rejected.
- **Route B — nodal-mass injection like SMS.** At `domainChanged`, add `m_p` to `Node`
  mass for each tied node, record it in a `bipenaltyInjected` map, and **re-baseline**
  (subtract before recompute) every `handle()`. *Energy closes for free*: the
  EnergyBalanceRecorder KE term reads `Node::getMass()` (`EnergyBalanceKernel.h:153`), so an
  injected `m_p` is auto-counted — **add the missing KE, never subtract** (finding C5: added
  explicit mass is *physical* inertia, like lumped SMS); no double-count (element vs nodal
  mass are disjoint KE buckets). **Chosen.**
- **Route C — real Element `getMass`.** Not available (contact is FE adapters). Rejected.

**Handler-lifecycle requirement (Round-2 NEW finding — load-bearing).** The handler does
**NOT** inherit SMS's re-baseline/restore. Verified: `LadrunoContactHandler::~…()` is empty
(`LadrunoContactHandler.cpp:73`), `clearAll()` never touches `Node` mass (`:935`), and
`handle()` only rebuilds a *read* cache (`ladrunoBuildNodalMass`). A naïve Route-B write
would therefore **compound** `m_p` every `domainChanged` and **leak** it permanently to a
later stage (e.g. a modal stage after the explicit one → wrong frequencies). The handler IS
persistent (held by `DirectIntegrationAnalysis`, not rebuilt per step), so the fix is
viable but must be **built**: the `bipenaltyInjected` map + a `removeBipenaltyMass()` at the
top of every `handle()` + a non-empty destructor that restores against the committed Domain.
A **bipenalty-only** teardown/no-compounding gate is required, separate from the
SMS-coexistence gate.

**Invariant note:** under `-bipenalty`, "contact adds zero mass" is intentionally broken —
via `Node` mass (Route B), not the FE `addMtoTang` (which stays a no-op). LEDGER_quirk.

### Decision 3 — Sizing details, safety, normal & tangential

- **Fix `k_p`; size `m_p` once per `domainChanged`** from the corrected deficit formula
  (§What-sets-`dt_cr`), via the **new decomposed gap-mode helper** (not the fused
  `gapModeInvMass` — finding A2). `Δt`-change guard mandatory (finding B3).
- **`safety` default 0.9** (RBE uses 1.0; we want the *smallest* safe `m_p` since bigger
  `m_p` = more perturbation, so `safety` close to 1). Monotonicity verified: smaller
  `safety` ⇒ larger `m_p` (more conservative). **OQ-3.**
- **No separated-pair floor in P1** — the tie has no separated state, so the question
  doesn't arise for the shipped scope. (The slave-surface floor is the *deferred*
  unilateral design, §Decision 1, not a P1 element.)
- **Tangential**: `m_p,t` from the assembled `B_t M⁻¹ B_tᵀ` worst-case tangent; inherits
  the SOFT isotropy caveat (use isotropic lumped mass under explicit).

### Decision 4 — Relationship to SOFT and ALM

- **Bipenalty vs SOFT: mutually exclusive on the *normal* spring** (opposite sizings of
  the same mode). Not claimed to be a deep principle (finding D5) — per-mode mixing
  (SOFT-normal + bipenalty-tangential) is *algebraically* possible but **out of scope**;
  it's a narrow corner of an already-narrow niche.
- **Explicit-only**, gated by `dynamic_cast<CentralDifferenceLadruno*>`, like SOFT. Under
  any implicit integrator (incl. ALM) `-bipenalty` is **inert** ⇒ byte-identical. ALM is
  the implicit-accuracy story; bipenalty is the explicit-accuracy story; they don't
  compose on one explicit run.
- **Not with SMS in P1** — the `dynamic_cast` catches `CentralDifferenceSMS` subclasses too;
  P1 **refuses** the SMS+bipenalty combination (Decision 5 — sizing-staleness deferred).

### Decision 5 — Interaction with SMS: **DEFERRED — refuse the combination in P1** (Round-3 correction)

v2/v3 claimed bipenalty automatically sizes from pre-SMS mass because `handle()` runs
before the integrator's `domainChanged()`. **Round-3 found this is FALSE after the first
domainChanged.** Verified: `CentralDifferenceSMS::domainChanged()` runs its `removeScaling()`
at its own top (`CentralDifferenceSMS.cpp:171`) — but that whole method is invoked at
`DirectIntegrationAnalysis` **line ~467**, *after* `handle()` (~423). So at `handle()` time
the **previous cycle's SMS injection is still on the nodes** (SMS removes it later in the
same cycle). Bipenalty reading `ladrunoBuildNodalMass` at `handle()` would therefore size
`m_eff` off **SMS-inflated** mass on every cycle past the first — over/under-sizing `m_p`
depending on whether SMS scaling drifts (stable only if SMS scaling is monotone, not
guaranteed in a nonlinear run).

**Resolution: SMS + bipenalty coexistence is DEFERRED to the successor ADR. P1 REFUSES the
combination** — if `-bipenalty` is active under a `CentralDifferenceSMS`/`…SMSConsistent`
integrator, warn and disable bipenalty (or error). This keeps the shipped scope clean
(bipenalty sizes off genuine physical mass, no SMS interference) and removes the fragile
ordering entirely. A proper fix (size from a physical-mass representation that excludes
*both* injection maps, e.g. element-density mass + original nodal `mass`) is successor-ADR
work, alongside the unilateral floor.
- **Energy (Route B alone, SMS deferred).** Route B mass is nodal ⇒ auto-counted in KE
  (like lumped SMS); **no registry conduit, no subtraction** (finding C5 corrects v1's
  "extends the Olovsson accounting" — nothing to extend).
- **Closure oracle — corrected AGAIN (Round-3 finding 1).** The penalty's **internal**
  energy `½k_p·gap²(t)` lives on the FE adapter the recorder never visits, so `RES` carries
  a per-step (non-constant) offset for *any* penalty contact ⇒ "expect RES≈0" is wrong. v3's
  replacement (a bipenalty-on vs -off **difference** test, asserting the KE Δ == `½ Σ m_p v²`)
  is **also unsound**: a stiff tie's `m_p` is *large*, so it changes the dynamics — the two
  runs do NOT share `gap(t)`/`v(t)`, the `½k_p·gap²` histories don't cancel, and `v` is
  ambiguous. **Correct oracle = energy CONSERVATION of the bipenalty run itself:** in a
  free-vibration / no-external-work tie test, the recorder total `KE+IE+DW` (now counting
  `m_p` in KE) must stay bounded — drift ≤ integrator tolerance — because the only
  un-recorded term, `½k_p·gap²`, is **small and bounded for a stiff tie** (`gap ≈ load/k_p
  → 0`). Pair it with the **accuracy** check (the tie holds: relative displacement
  `r ≈ load/k_p`, unchanged by `m_p` — mass doesn't soften the static bond) and the
  **stability** check (BLOCKER-A). The heavier publish-`½k_p·gap²`-into-IE conduit (→ RES→0
  honestly) is successor-ADR work.

### Decision 6 — Massless-node case

The corrected deficit formula **is** the massless fix: when `m_eff,real → ∞` (`s0`
large / node massless on the gap DOF), `m_p = 1/rhs − 0 = k_p·(Δt/2safety)²` — a real mass
that makes `M⁻¹` well-posed where SOFT today falls back + warns (`m_eff ≤ 0`). Reuses the
RBE "node already carries mass → don't double-count" guard, generalized to "carries
*enough* mass."

---

## Design-gate BLOCKERs

Re-scoped: the BLOCKERs that were unilateral-only now **define the deferral boundary**
(they are the reason the unilateral set is a successor ADR, not gates inside this one). What
remains as live gates for the shipped P0+P1:

**BLOCKER-A (LEAD) — corrected sizing proven on finite masses.** The headline P1 oracle
(`dt_cr` unchanged with stiff `k_p`) must pass on a **mortar tie with finite slave mass**,
using the **decomposed** gap-mode deficit formula — *not* the v1 `k_p(Δt/2)²` form (the
massless limit, which masks the bug). Gate: 2-DOF + multi-node numeric eigenvalue check
that `ω_contact·Δt ≤ 2·safety` with the **injected** `m_p` in `M`.

**BLOCKER-B — `m_p` reaches the SAME diagonal `system Diagonal` inverts, in the SAME
representation it was sized in.** Route B writes `Node` mass; sizing reads the assembled
cache. Gate: a unit test reads the factored diagonal, confirms the injected `m_p` changed
`M⁻¹`, and that sizing-mass == injection-mass representation (else `ω_p·Δt` is mis-targeted).

**BLOCKER-C — `Δt`-change guard.** Mid-run `Δt` change without a domainChanged invalidates a
domainChanged-cadence `m_p`. Gate: the guard warns/refuses (finding B3).

**BLOCKER-D — Route-B handler-lifecycle (bipenalty-only; Round-2 finding 4).** The handler
builds the `bipenaltyInjected` map + per-`handle()` re-baseline (`removeBipenaltyMass` runs
at the TOP of `handle()`, using the stored map + cached Domain, before adapters/cache are
rebuilt) + a non-empty destructor (restores against the cached committed Domain, as SMS's
`appliedDomain` does). Gate: a bipenalty-only run **does not compound** `m_p` across multiple
domainChanges **and** restores original `Node` masses exactly after teardown (a following
modal stage sees unperturbed frequencies). *(SMS coexistence is deferred — Decision 5 —
so P1 REFUSES `-bipenalty` under an SMS integrator; no coexistence gate ships.)*

**BLOCKER-E — energy CONSERVATION oracle (not RES≈0, not a difference test).** Round-3:
the `½k_p·gap²` offset is per-step (kills RES≈0) AND the on/off difference test is unsound
(large `m_p` changes the dynamics ⇒ runs don't share `gap(t)`). Gate: the **bipenalty run's
own** recorder total `KE+IE+DW` drifts ≤ integrator tolerance in a free-vibration tie test
(the un-recorded `½k_p·gap²` is small & bounded for a stiff tie), plus the tie-holds
accuracy check. Do **not** assert RES≈0 or a cross-run KE difference.

**Deferred (successor ADR):** the unilateral set (spurious-impact-energy, stale-inverted-
mass / no per-step write path, friction-stick orientation drift, the `rhs ≤ 0` too-light-
master fallback) + the **slave-surface permanent-floor** design that resolves them; **and
SMS+bipenalty coexistence** (size from a physical-mass representation excluding both
injection maps).

---

## Phased implementation plan

Each phase: numpy oracle → build-free C++ self-check → OpenSees regression.

- **P0 — TWO KILL-GATES + plumbing (no injection yet).**
  1. **Quantify (finding D1):** measure SOFT **mortar-tie** penetration vs a stated
     tolerance at a representative explicit `Δt`. Within tolerance ⇒ **stop** (pointer to
     SOFT/ALM).
  2. **RBE2-redundancy (finding D3, now decisive):** confirm the target tie is genuinely
     **facet-weighted / non-conforming** — something `LadrunoKinematicCoupling -bipenalty`
     (node-to-node) cannot express. Node-conforming ⇒ **stop, use RBE2.**
  3. If both pass: wire `-bipenalty -dtcr Δt [-safety s]` parse + the **decomposed gap-mode
     sizing** (sized, not yet injected). *Oracle:* deficit formula vs numpy 2-DOF
     eigenvalue; `-bipenalty` off ⇒ byte-identical.
- **P1 — mortar mesh-tie bipenalty (Route B).** Inject the deficit `m_p` via nodal mass,
  with the handler-built re-baseline + restore (BLOCKER-D); **refuse `-bipenalty` under an
  SMS integrator** (Decision 5). *Oracles:* **(a)** `dt_cr` **unchanged** with stiff `k_p`
  on a **finite-mass** tie (BLOCKER-A — headline); **(b)** tie-holds accuracy = the stiff
  penalty (no penetration cap, `m_p` doesn't soften the bond); **(c)** energy
  **conservation** — the bipenalty run's `KE+IE+DW` drift ≤ tol (BLOCKER-E), NOT RES≈0 nor a
  difference test; **(d)** massless-node stability (zero-mass tied node: SOFT falls
  back+warns, bipenalty stable); **(e)** bipenalty-only no-compound + restore-after-teardown
  (BLOCKER-D). **This is the whole shipped feature; trust gate.**
- **Deferred successor ADR** — unilateral RIGID_PLANE / NTS / mortar-contact / edge-edge,
  friction, the appearing/disappearing-mass problem; **the slave-surface permanent-floor**
  is its design starting point. Full multi-agent adversarial gate.

Oracle infra mirrors the SOFT validation (`proto_b1_soft_penalty.py`) and the RBE
bipenalty legs (`tests/test_ladrunoKinematicCoupling_element.py`).

---

## Adversarial gate decision

- **P0, P1** — **lighter review.** With scope cut to the mortar tie (constant active set),
  P1 is a near-direct port of the already-gated RBE2 bipenalty + the already-gated SOFT
  infrastructure; the corrected sizing + the handler-lifecycle build are the real work, and
  BLOCKER-A/D1/E carry them. Rounds 1–2's full adversarial passes already did the hard
  attacking (this ADR is their distillate).
- **Deferred successor ADR** — **full multi-agent gate**, as the unilateral findings (incl.
  the slave-surface floor) are the "novel math + core invariant break" profile.

---

## Ledger / classTag bookkeeping

- `LEDGER_implementations.md` — one row: *contact bipenalty mode* (flag on
  `LadrunoContactFE` + `LadrunoContactHandler`), **no new class tag**, status per phase, PR.
- `LEDGER_vanilla_files.md` — only if a vanilla touch slips in (expected NONE).
- `LEDGER_quirks.md` — four entries: (1) *`-bipenalty` sizes `m_p` from the **decomposed
  assembled** gap-mode mass (slave & master terms separately), not `k_p(Δt/2)²` — the bare
  form is the massless-limit and under-stabilizes finite masses; the fused `gapModeInvMass`
  is insufficient*; (2) *`-bipenalty` breaks "contact adds zero mass" via `Node` mass
  (Route B), explicit-only, default off — the FE `addMtoTang` stays a no-op*; (3) *the
  contact HANDLER does NOT inherit SMS's re-baseline/restore — a Route-B nodal write needs a
  `bipenaltyInjected` map + per-`handle()` re-baseline + non-empty destructor or `m_p`
  compounds per domainChanged and leaks to later stages*; (4) *`-bipenalty` is REFUSED under
  a `CentralDifferenceSMS` integrator in v1 — at `handle()` time the previous cycle's SMS
  injection still sits on the nodes (SMS's `removeScaling` runs later, inside the integrator
  `domainChanged`), so bipenalty would size off SMS-inflated mass; coexistence is deferred*.
- `classTags.h` — **no change** (next free contact ELE slot is **33022**; 33015-33021 are all taken — the earlier 33015/33016 reservations were superseded by ADR-58/70/71/72).
- Banner — no new feature line unless we surface "contact bipenalty" distinctly.

---

## Open questions — need sign-off before coding

- **OQ-A (THE decision — build at all?) — P0 NOW RUN; recommendation is SHELVE.** Both
  kill-gates were measured (§P0 kill-gate results): Gate 2 passes but is narrow (unique only
  for non-conforming ties); Gate 1 shows SOFT is adequate (sub-% of `h`) outside a
  high-interface-stress corner, and bipenalty's accuracy edge costs a ~β× interface-mass
  inflation at equal stability. **Your call:** (i) accept SHELVE — leave the ADR as a
  complete, gated design; or (ii) green-light P1 because you have the specific corner case
  (non-conforming + explicit + high interface stress + tight tolerance). This is the only
  question that needs you; everything else is resolved.
- **OQ-5 (scope) — RESOLVED to "mortar-tie only"** (Round-2: rigid-plane is unilateral, not
  constant-active-set). Defer all unilateral to a successor ADR. Confirm.
- **OQ-1 (route) — RESOLVED to Route B** (nodal injection; the only route whose KE
  auto-counts). Confirm.
- **OQ-2 (floor) — RESOLVED:** no floor in P1 (the tie has no separated state); the
  **slave-surface** permanent floor is the *deferred* unilateral design (v2's blanket
  rejection was of the wrong, all-nodes floor). Confirm.
- **OQ-3 (safety default).** `safety = 0.9`?
- **OQ-6 (sizing mode) — RESOLVED to `-dtcr` only.** Confirm.

---

## Adversarial review log

**Round 1 (4 independent reviewers, code-grounded).** Converged on:
- **A (mechanics).** The v1 identity `ω_p=√(k_p/m_p)` is the assembled max frequency only
  in the massless/rigid-master limit (the RBE2-guarded regime); for finite masses
  `ω²=k_p·(B M⁻¹ Bᵀ)` under-stabilizes 18–58%. Fix: deficit sizing on `gapModeInvMass`.
  The "independent of `m_eff` — that's the whole point" sentence was the bug; deleted.
- **B (architecture).** `m_p` is domainChanged-cadence, not per-step (mass ≠ residual);
  neither route beats that cadence; active set changes with no domainChanged ⇒ unilateral
  has no per-step write path ⇒ deferred. `Δt`-change guard added.
- **C (energy).** Contact FE adapter is invisible to EnergyBalanceRecorder ⇒ RES can't →0
  for any penalty contact (oracle redefined). Route B nodal mass auto-counts KE (add, not
  subtract); Route A would need a registry conduit it can't share with SMS ⇒ Route B
  chosen. "(c) route through nodal mass" makes ΔKE *accounted*, not *physical* ⇒ only a
  constant `m_p` (P1) is clean.
- **D (scope).** Value unquantified vs SOFT; every BLOCKER is unilateral-only; P1 may
  duplicate RBE2 `-bipenalty`. ⇒ Added P0 kill-gate, cut to P1-only, rejected the floor,
  `-dtcr`-only, fixed stale 33015→33016.

**Round 2 (3 reviewers: verify-math, attack-route/energy, steelman-the-cut).** Verified the
v2 corrections HOLD (deficit algebra sound; Route-B KE auto-count + no-double-count; SMS
insensitivity; safety² monotonicity) — strong convergence on fundamentals. New load-bearing
fixes folded into v3:
- **A2 (math).** The fused `gapModeInvMass` does not expose slave/master terms separately →
  a new decomposed helper is required (the "reuse `gapModeInvMass`" wording was wrong).
  Reconciled the reciprocal `s_target` ↔ mass forms.
- **Route/energy.** SMS-before-bipenalty ordering is **forced by snapshot timing**, not a
  choice (corrected rationale). Energy oracle "known constant offset" was wrong — `½k_p·gap²`
  is per-step → switched to a **bipenalty-on/off energy-difference** test. **NEW:** the
  handler has **no** SMS-style re-baseline/restore (empty dtor, no per-`handle()` subtract) →
  Route-B `m_p` would compound + leak → must build the lifecycle machinery; split BLOCKER-D.
- **Scope (the big one).** **RIGID_PLANE is unilateral** (verified `getResidual` gates on
  `g<0`) → v2 mis-classified it. P1 cut to **mortar-tie only**; build-at-all now rides on the
  RBE2-redundancy gate (surviving niche = facet-weighted non-conforming explicit tie). The
  permanent **slave-surface** floor recorded as the *known route* into the deferred
  unilateral case ("intractable" was too strong).

**Round 3 (code-grounded self-review — Agent classifier was down; done via Read/Grep on
source).** Verified the v3 fundamentals against code: the DIA call order (`handle()` @423
before `integrator->domainChanged()` @467), the decomposed helper is a trivial split of
`gapModeInvMass` (`LadrunoContactFE.cpp:429-433`, no vanilla edit), and the handler-lifecycle
prescription is sound (map survives `clearAll`; destructor can cache the Domain like SMS).
Two NEW load-bearing findings folded into v4:
- **R3-1 (energy oracle).** v3's bipenalty-on/off **difference** test is unsound — a stiff
  tie's `m_p` is large and changes the dynamics, so the two runs don't share `gap(t)`/`v(t)`
  and the `½k_p·gap²` histories don't cancel. Replaced with an **energy-CONSERVATION**
  (bounded-drift) check on the bipenalty run itself (the un-recorded `½k_p·gap²` is small &
  bounded for a stiff tie) + the tie-holds accuracy check.
- **R3-2 (SMS staleness).** v3's "snapshot timing ⇒ bipenalty sizes off pre-SMS physical
  mass" is FALSE after the first domainChanged: SMS's `removeScaling` runs at
  `CentralDifferenceSMS.cpp:171`, *inside* the integrator `domainChanged` (@467) which is
  AFTER `handle()` (@423), so the prior cycle's SMS mass lingers on the nodes at sizing time.
  **Resolution: defer SMS+bipenalty coexistence; P1 REFUSES the combination.** BLOCKER-D2
  dropped; Decision 4 + the quirk updated.

**Round 4 — CONVERGED.** Re-attack on v4 surfaced only minor items (the energy-conservation
oracle should note "no *secular* drift; bounded oscillation = the un-recorded `½k_p·gap²`" —
already implied; and the SMS-refusal narrows utility, which feeds OQ-A, not a design flaw).
No new load-bearing findings ⇒ the loop's stop condition is met. The core (deficit sizing,
Route B, mortar-tie-only scope, RBE2-redundancy build gate) has been stable since v2–v3 and
is verified against source. The only unresolved item is OQ-A, a user sign-off, not a
technical gap.

**P0 kill-gates — RUN (v5).** Gate 2 PASS-but-narrow; Gate 1 shows SOFT adequate outside a
high-stress corner; the β-reciprocal penetration↔interface-mass tradeoff ⇒ **recommend
SHELVE**. See §P0 kill-gate results.
