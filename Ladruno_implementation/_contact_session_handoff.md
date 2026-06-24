---
title: Ladruno contact — session handoff (post C4 + D2)
project: Ladruno
status: handoff — mortar lane DONE (C1→C4) + viscous stabilization DONE (D2.1 NTS + D2.2 mortar)
owner: nmora
tags:
  - contact
  - handoff
  - roadmap
---

# Ladruno contact — session handoff

> Pick-up doc for the next session. The **authoritative** status/roadmap is the capstone
> [[48_ladruno_contact_capstone_adr]] (status-of-record table + unified roadmap); this is a
> point-in-time snapshot + recommended next step. Deferred set: [[47_ladruno_contact_deferrals_adr]].

## State (as of this handoff)

- **B3 SHIPPED** ([#389](https://github.com/nmorabowen/OpenSees/pull/389)) — the consistent ∂n/∂u
  geometric NORMAL tangent (`contact … -geomtan`) + the `ladrunoContactForce` query + a Hertz
  benchmark. Closes the geometric-tangent deferral the NTS lane carried since P2b (C2/C3/C4 all
  deferred it). `K_geom = kn·gN·∂²gN/∂u²` ⇒ quadratic Newton on curved/large-sliding; symmetric;
  gated off-default (indefinite → opt-in like `-consistanttan`). Battery **83→91/91**.
  NOTE: quantitative elliptic-Hertz `p(r)` is resolution-bound by NTS-penalty pressure recovery —
  4 new [[LEDGER_quirks]] entries (handler imposed-disp; contact force not in nodeReaction;
  geom tangent indefinite; curved-indenter penalty bootstrap). Finer-mesh + displacement-control /
  D1 augmentation is the documented follow-up.
- **Imposed-displacement under `LadrunoContact` now WARNS** ([#392](https://github.com/nmorabowen/OpenSees/pull/392)) —
  the contact handler is Plain-style (must inject the FE adapters in `handle()`), so like stock
  `PlainHandler` it does NOT enforce a non-homogeneous SP (imposed displacement). It was doing so
  SILENTLY; now it emits the `PlainHandler`-style warning pointing at `DisplacementControl` (which
  DOES drive a DOF via the load factor with this handler). Not a regression; full imposed-SP support
  (a Transformation-style contact handler) is deferred — no consumer.
- **B1 SHIPPED** ([#402](https://github.com/nmorabowen/OpenSees/pull/402)) — the LS-DYNA §26.15 **SOFT=1 Courant-stable explicit penalty**
  (`contact … -soft <SOFSCL>` / `contactPlane … -soft <SOFSCL>`, NTS lane). Under explicit
  `CentralDifferenceLadruno` the contact kₙ is replaced each step by `k_soft = SOFSCL·4·m_eff/dt²`
  (SOFSCL default 0.10) so the contact's own `ω·dt = 2√SOFSCL ≤ 2` ⇒ explicit dt_cr is NEVER
  penalty-throttled — impact/pounding/recontact runs at the STRUCTURAL dt (the 50_/51_ progressive-
  collapse / element-removal enabler). `m_eff` is the gap-mode generalized mass from the ASSEMBLED
  nodal mass (nodal `mass` + element-density mass — the handler pre-builds a per-node cache; the B1
  adversarial gate caught + fixed sourcing it from `Node::getMass()` alone). Explicit-only (dt via
  `dynamic_cast`→CDL); implicit + `-soft`-absent byte-identical. Battery **91→98/98**.
- **B2 SHIPPED** ([#406](https://github.com/nmorabowen/OpenSees/pull/406)) — the LS-DYNA §26.15 **SOFT=2 segment-based explicit penalty**
  (`contact … -mortar -soft <SOFSCL>`). The shipped B1 SOFT=1 is NODE-to-segment ⇒ it MISSES (1) partial
  facet overlap with no slave node in-bounds, (2) a node sliding off a master segment edge, (3) a
  T-intersection. B2 REUSES the shipped `LadrunoMortarKernel` clip→Gauss (D,M,g̃,n over the slave∩master
  overlap) and distributes a Courant-stable penalty `p_I = min(0, k_soft,I·ḡ_I)`, `k_soft,I =
  SOFSCL·4·m_eff,I/dt²` with the SEGMENT-BASED gap-mode mass `m_eff,I = 1/(B_I M⁻¹ B_Iᵀ)`,
  `B_I = [D,−M]/a_I` (reuses the B1 assembled nodal-mass cache) ⇒ each contact mode `ω·dt = 2√SOFSCL ≤ 2`,
  dt_cr un-throttled, while CATCHING the corner/edge/T cases NTS misses. New `LadrunoContactFE::
  addSoft2Penalty` + an explicit-only fast path at the top of the MORTAR `getResidual`
  (`softScale>0 && dynamic_cast<CentralDifferenceLadruno*>`; strictly additive ⇒ `softScale=0`
  byte-identical); implicit falls through to the shipped penalty/ALM (explicit-only, like B1).
  FRICTIONLESS MVP (⊥ `-tie/-mu/-cohesion/-tauMax/-visc`). 3-reviewer adversarial gate PASS (1 MINOR
  folded — the per-node Courant bound is necessary-not-sufficient; coupled `K_c` ω·dt ~2× ⇒ a
  `SOFSCL>0.25` warning + [[LEDGER_quirks]] note; default 0.10 safe). Battery **98→104/104**.
- **D1 SHIPPED** ([#411](https://github.com/nmorabowen/OpenSees/pull/411)) — within-step augmentation sign-off (see Track D1 below): the
  `analyze_augmented` proc is now a first-class no-corruption recipe (`ladrunoBeginAugment`/
  `ladrunoEndAugment` + the `Domain::contactAugmenting` commit guard). Battery **106→112/112**.
- **Contact battery: 112/112** — `tests/test_adr39_contact_p*.py` (+ `_p2b2c_hertz.py`, `_p4_soft.py`,
  `_p5_soft2.py`) + `test_adr41_mortar_c2_{0,1,2}.py` + `_c3_{1,2,3}.py` + `_c4_{1,2}.py` +
  `test_adr41_viscous_d2.py` + **`test_adr41_d1_augmentation.py`**.

### Shipped this session
| PR | What |
|---|---|
| #381 / #382 / #383 | **C4 mortar mesh-tying** — permanent bond (zero-gap limit of contact): full 3-vec weighted relative-displacement `r_I → 0`, no-clamp ALM `λ_tie`. Resolved the shared-node MAJOR-1 quirk for the tie (linear `r` ⇒ the order-independent global accumulator `λ_N` uses). |
| #384 | Capstone status-of-record sync (C4). |
| #385 / #386 | **D2.1 viscous stabilization (NTS)** — `-visc μ_c`, RIGID_PLANE + SEGMENT. |
| #387 / #388 | **D2.2 mortar viscous** — extends `-visc` to the mortar contact ⇒ **D2 COMPLETE**. |
| #389 / #390 | **B3 consistent ∂n/∂u normal tangent** (`-geomtan`) + `ladrunoContactForce` query + Hertz benchmark (#390 = ledger ref backfill). Closes the NTS geometric-tangent deferral; battery 83→91/91. |
| #392 | **Imposed-disp warning** — `LadrunoContact` now warns (was silent) on a non-homogeneous SP, points at `DisplacementControl` (B3 follow-up). |
| [#402](https://github.com/nmorabowen/OpenSees/pull/402) | **B1 SOFT=1 Courant-stable explicit penalty** — `-soft <SOFSCL>` (NTS); explicit dt_cr no longer penalty-throttled. Adversarial gate fixed a MAJOR (m_eff must use the assembled, not nodal-only, mass). |

Workflow fix: recorded the `gh pr checks --watch` premature-exit trap in
[[../Ladruno_internal/WORKFLOW_GOTCHAS]] §3 (verify the actual Zone-A run by id, not the watch exit).

### Build / test recipe (this Windows box)
- Build: `Ladruno_scripts\build.bat OpenSeesPy` (exit 0).
- Run: `python3.12` with `PYTHONPATH=dist/bin LADRUNO_OPENSEES_QUIET=1 PYTHONIOENCODING=utf-8`.
- Oracles (numpy-only): `Ladruno_implementation/contact_prototypes/proto_*.py`.
- Zone-A CI is the merge gate; **verify the actual run id**, not `gh pr checks --watch` (it can exit
  before Zone-A registers): `RID=$(gh run list --branch <br> --limit 1 --json databaseId -q '.[0].databaseId'); gh run watch $RID; gh run view $RID --json conclusion`.
- Merge: `gh pr merge <n> --admin --squash`, base `ladruno`.

## Roadmap — what's left

### ✅ Done
- **Track A** (shared kernels): A1 friction extract #367, A2 projection extract #365.
- **Track C** (mortar lane): C1 kernel → C2 commit-cycle ALM → C3 friction → C4 mesh-tying. The whole
  implicit/accuracy lane.
- **Track D2** (viscous stabilization): NTS #385 + mortar #387.

### 📋 Remaining — COMMITTED (scoped, oracle-first ready)

**Track B — NTS explicit-stability hardening** (the explicit-first lane; independent of mortar):
- ~~**B1** — P4 `SOFT=1` Courant-stable penalty~~ **DONE [[#402](https://github.com/nmorabowen/OpenSees/pull/402)]** (see State above). Follow-ups:
  source m_eff from the assembled M for the parallel/distributed (row-sum lumped) SOE — serial-only
  today. **SOFT on the tangential `kt` is now DONE** ([#409](https://github.com/nmorabowen/OpenSees/pull/409) — `softKt` sizes the friction stick
  penalty Courant-stably `k_soft_t = SOFSCL·4·m_eff_t/dt²` so a stiff `kt` no longer throttles dt_cr via
  the stick mode; battery 104→106). (The mortar SOFT penalty B1 listed is delivered by **B2** above.)
- ~~**B2** — P5 `SOFT=2` segment-based penalty (corner / edge / T-intersection robustness)~~ **DONE
  [#406](https://github.com/nmorabowen/OpenSees/pull/406)** (see State above). Follow-ups: frictional / viscous SOFT=2; the perpendicular
  edge-edge case (cos_t→0 ⇒ the mortar clip degenerates) needs a dedicated edge-edge treatment.
- ~~**B3** — P2b-2c consistent `∂n/∂u` normal tangent + Hertz~~ **DONE [#389]** (see State above).
  Follow-up spun off: a quantitative FE Hertz harness (contact-force recorder ✅ via
  `ladrunoContactForce`; a robust curved-indentation driver — displacement-control or D1
  within-step augmentation — for the elliptic-`p(r)` / compliant-radius match).

**Track D1 — within-step augmentation refinement**: ✅ **DONE** ([#411](https://github.com/nmorabowen/OpenSees/pull/411)) — the held-load
`analyze_augmented` proc is now a FIRST-CLASS no-corruption recipe. It brackets its zero-increment
`LoadControl` re-commit loop with new `ladrunoBeginAugment`/`ladrunoEndAugment` commands that set a
`Domain::contactAugmenting` flag making `Domain::commit()` skip the recorder loop + `commitTag++`
(the Uzawa λ update still runs). So the within-step sweep drives `‖ḡ‖→augTol` (contact) /
`‖r̄‖→augTol` (tie) inside ONE physical step without recorder / load / time corruption (the capstone
contract #3 sign-off). `try/finally` ⇒ the flag can never stick ON. Oracle `proto_d1_augmentation.py`
11/11; `test_adr41_d1_augmentation.py` 6/6; battery 106→**112/112**; 3-reviewer gate PASS (2 MINOR
folded). This ALSO unlocks the quantitative Hertz follow-up via displacement-control-free curved
indentation. **The committed contact roadmap (Tracks A→D) is now COMPLETE.**

### 🚫 Deferred — ADR-47 ([[47_ladruno_contact_deferrals_adr]]; each gated behind a re-open trigger that has NOT fired)
1. Dual / biorthogonal mortar basis (diagonal `D`, cheap nodal λ, LBB-optimal).
2. True-LM / saddle-point enforcement + inf-sup stabilization.
3. Self-contact (a surface vs itself).
4. Full slide-line Hermite smoothing (4) / averaged nodal-normal smoothing (4a) — faceted-normal chatter.
5. Anisotropic / elliptic friction (two principal μ).
6. Pressure/velocity-dependent `μ(N,v)` wired into contact (machinery exists in `frictionModel/`, unwired).
7. Custom global-Uzawa `LadrunoAugmentedNewton` `EquiSolnAlgo`.
8. NTN / NTS-via-mortar-weights unification.
9. Coupled-multiphysics surfaces (pore / thermal / acoustic) + pressure penetration.

### Disclosed limitations (not deferrals — won't change)
- Full ALM is **implicit-only** (explicit = single-pass penalty under the mass-only CDL).
- Both lanes are **penalty / Uzawa**, not exact-Lagrange stick.

## Recommended next step

**The committed contact roadmap (Tracks A→D) is COMPLETE** (D1 shipped this session). What remains is
all OPTIONAL / deferred:
- **Small follow-ups** (each its own small PR, oracle-first): ~~viscous SOFT=2~~ **DONE** ([#412](https://github.com/nmorabowen/OpenSees/pull/412) —
  `-mortar -soft -visc`; the D2.2 normal damper on the SOFT=2 active set in `addSoft2Penalty`,
  `μ_c=0` byte-identical; battery 112→114/114). ~~frictional SOFT=2~~ **DONE** ([#414](https://github.com/nmorabowen/OpenSees/pull/414) —
  `-mortar -soft -mu/-cohesion/-tauMax`; Courant-stable Coulomb/Tresca via a per-node segment `softKt`
  [the B1-kt rule, n→t] + the shipped return map over the soft pressure; penalty friction λ_T≡0;
  reuses the committed `MortarNormalState` slot; only `-tie` still refused; μ=0 byte-identical; oracle
  `proto_soft2_friction.py` 12/12, `test_adr41_soft2_friction.py` 4/4; battery 114→119/119).
  **⇒ SOFT=2 now supports viscous + Coulomb/Tresca friction.** assembled/row-sum `m_eff` for the
  parallel/distributed SOE — **BLOCKED, not small**: the contact handler `sendSelf`/`recvSelf` are P1a
  single-process stubs (the whole contact subsystem is serial-only), so parallel `m_eff` requires
  parallelizing the handler first (its own ADR-scale effort, not a follow-up).
- **NOT small** — the perpendicular edge-edge treatment (a dedicated edge-to-edge contact for the
  `cos_t→0` case the mortar face-overlap clip degenerates on; a new algorithm — scope it as its own
  ADR-47-adjacent item).
- **Validation** — the quantitative elliptic-Hertz `p(r)` study, now unblocked by D1's
  displacement-control-free curved indentation.
- Everything else is **ADR-47** ([[47_ladruno_contact_deferrals_adr]]), deferred behind re-open triggers.

Process per fork standard: oracle-first → C++ → build → adversarial gate → PR on `ladruno`,
keep the capstone + ledgers current in the same PR.

### ⏭️ Queued for a FRESH session (not this one) — two INDEPENDENT efforts

Both are from-scratch DESIGN work (a new ADR each, not a follow-up PR), so start them in a clean
context with the relevant skills loaded, not by extending an implementation session. They do NOT
depend on each other — scope them as two separate ADRs/sessions.

**1. Perpendicular edge-edge contact ADR** (the `cos_t→0` mortar-clip degeneration — a NEW algorithm).
✅ **DELIVERED — [[57_ladruno_edge_edge_contact_adr]]** (2026-06-24, design-only): segment-segment
closest point + common-perpendicular normal `n=(ê_s×ê_m)/‖ê_s×ê_m‖` + the two-stage cos_t routing
detector (facet-normal alignment routes to edge-edge as →0; edge-direction `sinθ_edge`→0 refuses the
parallel degeneracy; clip-area is the partition arbiter). New header `LadrunoEdgeKernel.h` + a new
`LadrunoContactFE` `EDGE_EDGE` mode + `LadrunoContactDomain::EdgeEdgeState`; reuses the bucket sort
(edge pairs derived inside the near facet neighborhood — zero new spatial structure), the friction
kernel, the SOFT mass cache, and a one-scalar ALM. Phase gates E0→E7 (oracle-first). Fenced:
predefined-surface only — self-contact/runtime-discovery edge-edge → ADR-55, smoothing → ADR-47 #4a.
Capstone B2 row + component table + ADR-47 ledger row #10 updated. **Implementation is a later session.**
Original kickoff prompt (for reference):
> Scope the perpendicular edge-edge contact ADR for the Ladruno fork (the `cos_t→0` case where the
> mortar face-overlap clip degenerates — two facets meeting edge-on, no in-plane overlap to integrate).
> START by reading `Ladruno_implementation/48_ladruno_contact_capstone_adr.md` (the status-of-record
> table + the B2 "perpendicular edge-edge deferred" note) + `47_ladruno_contact_deferrals_adr.md` +
> this handoff (`_contact_session_handoff.md`). Use the `opensees-contact` skill (its
> `broad_phase_collision` edge-edge / self-contact module + `contact_formulations`) and
> `abaqus-theory-contact-loading` for the closest-point-between-two-segments / edge-edge kinematics.
> DELIVER a design ADR (NOT code): the governing kinematics (segment-segment closest point, the
> edge normal + gap, the degenerate-`cos_t` detector that routes a pair to edge-edge vs face-mortar),
> how it slots into `LadrunoContactBucketSort` + `LadrunoContactDomain` + the `LadrunoContactFE`
> adapter, the per-phase gates (oracle-first), and an explicit fence vs ADR-47. Number it in the ADR
> sequence; keep it design-only — implementation is a later session.

**2. Contact handler parallelization** (unblocks parallel/row-sum `m_eff` for OpenSeesMP).
The `LadrunoContactHandler::sendSelf`/`recvSelf` are P1a single-process stubs ⇒ the whole contact
subsystem is serial-only. Kickoff: scope an ADR for parallelizing the contact handler under
`PartitionedDomain`/OpenSeesMP — `sendSelf`/`recvSelf` serialization of the `LadrunoContactDomain`
state, partition-boundary pairing, and the cross-partition mass reduction the `m_eff` SOFT sizing
needs. Read the same capstone + handoff; use the `opensees-performance` skill (parallelism /
`PartitionedDomain` / `ParMETIS`) + `opensees-expert`. Design-only ADR first.

**3. Validation (smaller, could ride an implementation session):** the quantitative elliptic-Hertz
`p(r)` study, now unblocked by D1's displacement-control-free curved indentation — a finer-mesh harness
+ the `ladrunoContactForce` recorder. Not a new algorithm; lower priority than the two ADRs above.
