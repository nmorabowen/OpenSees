---
title: Validation Plan — LadrunoEmbeddedNode v1 (world-class node-to-host embedment)
project: Ladruno
status: plan
priority: high
owner: nmora
tags:
  - validation
  - plan
  - element
  - embedded-node
  - mesh-tie
  - penalty-method
  - explicit-dynamics
created: 2026-06-07
related:
  - "[[23_ladruno_embedded_node_adr]]"
  - "[[20_ladruno_embedded_reinforcement_adr]]"
  - "[[LadrunoEmbeddedRebar_guide]]"
  - "[[LadrunoBondSlip_guide]]"
  - "[[17_finite_strain_validation_plan]]"
---

# Validation Plan — `LadrunoEmbeddedNode` v1

> [!abstract] What this document is
> The **verification & validation (V&V)** battery that earns `LadrunoEmbeddedNode` (ELE
> 33006) its **"world-class node-to-host embedment"** claim — benchmarked against **Abaqus
> `*EMBEDDED REGION`** and **LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID`**. It is a **plan**: each
> row is a contract a Zone-A / Zone-B test fulfils, with an **oracle** and a **tolerance**.
>
> **Scope is deliberately narrow** (the [[23_ladruno_embedded_node_adr#14 v1 scoping status supported core vs experimental credibility re-scope 2026-06-07|ADR §14]] re-scope):
> this battery validates the **v1 supported core only** — the **U** translational tie, the
> **`g0` stress-free birth**, and **penalty / AL / bipenalty** enforcement. **UR, UP, D9 and
> `-corot` are EXPERIMENTAL and explicitly OUT of this battery** (§5).

Companion docs: [[23_ladruno_embedded_node_adr]] (the element + the v1 scoping §14),
[[20_ladruno_embedded_reinforcement_adr]] (the sibling rebar element + the shared
`-enforce`/bipenalty/`-cfl` machinery), `tests/test_ladrunoEmbeddedNode_element.py` (where
these contracts live today).

---

## 1. The bar — what "world-class" means here

A node-to-host embedment is world-class when it matches the **capability** of the commercial
references **and exceeds** the one OpenSees precedent on conditioning and explicit dynamics.

| Reference | What it does | Our v1 core parity |
|---|---|---|
| **Abaqus `*EMBEDDED REGION`** | translational DOFs of embedded nodes **interpolated** from the host; embedded **rotations left FREE**; embedded nodes **born consistent** at activation | **U tie** (host-agnostic) + **rotations free by default** + **`g0` stress-free birth** |
| **LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID`** | **penalty**, **non-matching**, normal+tangential; **slip/bond optional**; rigorous friction = separate `*CONTACT` | **penalty/AL/bipenalty** core; D9 bond/cohesive is the *optional* (experimental) add-on; rigorous friction → `LadrunoContact` |
| **`ASDEmbeddedNodeElement`** (the OpenSees precedent) | tri/tet-only, raw `1e18` penalty, **implicit-only** | **exceed:** hex/tet/quad hosts, **`-kt auto`** conditioning, **AL**, and **clean EXPLICIT** (`getExplicitCriticalTimeStep` + `CriticalTimeStep` fold-in) |

**The differentiator to foreground:** correctness (patch / equilibrium / objectivity) is table
stakes; the world-class claim rests on **stress-free staged birth** + **clean implicit *and*
explicit** conditioning. Tests T3 and T5 are the headline.

## 2. Two-tier strategy (license-free rigorous gate)

- **Tier 1 — rigorous, in-tree, analytical oracles (the PASS/FAIL gate).** Every v1-core claim
  is provable with a closed-form or exact-interpolation oracle and **no commercial license**.
  This tier alone earns the credibility claim.
- **Tier 2 — cross-code credibility (Abaqus / LS-DYNA).** Documented decks + a comparison
  metric, run **if a license is available**; otherwise a **documented analytical surrogate**
  (shear-lag / springs-in-series) stands in. **Not a gate** — a credibility cross-check.

Zone tags follow the testbed (`testbed/00_canonical_testbed`): **Zone-A** = upstreamable,
OpenSeesPy/numpy-only, travels with the PR, gates CI; **Zone-B** = fork-local, needs gmsh.

## 3. The battery

### T1 — Patch / consistency (translational-tie correctness) · Tier 1 · Zone-A · 2D + 3D
**Claim:** a node embedded at an arbitrary **interior** point ξ of a single host element
**follows the interpolated host motion** under an affine field.
**Setup:** one host element (`LadrunoBrick` hex / `BezierTet10` / `FourNodeQuad`); prescribe an
affine displacement `u(x) = a + G·x` on the host nodes; embed a free node at an off-centroid ξ
(weights via `-xi`); solve. **Oracle (exact):** the embedded node's converged displacement
equals `Σ N_i(ξ) u_host,i = a + G·x(ξ)` (partition of unity reproduces affine fields exactly).
**Tolerance:** `1e-8` (relative), i.e. machine-exact up to the penalty residual `F/K_u` → use
`-enforce al` (or large `K_u`) so the residual is below tol.
**Status:** *partially covered* — `test_xi_offcentroid_matches_trilinear` checks the host
**weights**; **gap to add:** the affine-field *follow* test (free embedded node tracks the
interpolated motion), 2D + 3D.

### T2 — Pull-out vs Abaqus Embedded Region · Tier 2 · Zone-B (Zone-A surrogate) · rotations free
**Claim:** the coupling force / host reaction of a bar (or beam) node embedded in a solid and
pulled matches the Abaqus Embedded Region response, **with embedded rotations left free**
(Abaqus-consistent).
**Setup:** a solid block (meshed, non-matching), a single embedded node pulled by a prescribed
displacement; measure the coupling force vs imposed slip and the host reaction.
**Oracle:** *Tier-2 primary* = **analytical surrogate** — a **springs-in-series / shear-lag**
model (embedment penalty `K_u` in series with the host's local stiffness `k_host ≈ E·lch`),
coupling force `F = K_eff·δ`, `1/K_eff = 1/K_u + 1/k_host`; *Tier-2 optional* = an Abaqus
`*EMBEDDED REGION` deck (documented) when a license is available.
**Tolerance:** surrogate `≤ 5%` on `K_eff` (mesh-dependent `k_host`); Abaqus cross-check `≤ 5%`
on the force–slip curve. **Rotations free** asserted (no `-rot`).
**Status:** *new* — Zone-B meshed case + the Zone-A analytical-surrogate leg.

### T3 — Stress-free staged birth (`g0`) · Tier 1 · Zone-A · **headline**
**Claim:** an embedment added to an **already-deformed** host (staged construction) is born
**force- and stress-free** — no jolt, zero relative gap at activation — and tracks only
**subsequent** host motion. The decisive "world-class vs toy" test.
**Setup:** settle a `LadrunoBrick` host under load and commit (stage 1); **then** add a fresh
slave node tied into the deformed host (stage 2); read the element at activation.
**Oracle (exact):** initial resisting force `≈ 0`; relative gap `≈ 0`; captured offset
`= −Σ N_i u_host,i`. **Contrast (negative control):** with `-absolute` (no `g0`) the initial
force `= K_u·|Σ N_i u_host,i|` (the spurious spike) and the slave is **yanked** to the host
point.
**Tolerance:** force `< 1e-3`, gap `< 1e-9`; contrast force matches `K_u·|gap|` to `1e-3`.
**Status:** **DONE** — `test_staged_activation_born_stress_free`,
`test_staged_activation_absolute_yanks`, `test_staged_capture_noop_when_undeformed` (PR #214).
This plan row is **already fulfilled**; listed because it is the headline claim.

### T4 — LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID` coupling benchmark · Tier 2 · Zone-B (Zone-A surrogate)
**Claim:** the coupling-force response of a **non-matching beam-in-solid** reproduces the
LS-DYNA penalty-coupling behaviour (normal + tangential), with bond/slip as the optional add-on.
**Setup:** a beam line embedded through a solid block via a string of `LadrunoEmbeddedNode`
ties (one per beam node), loaded transversely and axially.
**Oracle:** *Tier-2 primary* = **analytical surrogate** — distributed penalty bedding
(Winkler-like) `q = K_u·w` per tie, closed-form tip response of a beam on an elastic bed;
*Tier-2 optional* = an LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID` deck (documented: `CDIR`, coupling
stiffness, penalty) when a license is available. **Comparison metric:** tip displacement and
total coupling force vs imposed load.
**Tolerance:** surrogate `≤ 5%`; LS-DYNA cross-check `≤ 10%` (penalty-stiffness dependent).
**Status:** *new* — Zone-B meshed case + the Zone-A surrogate leg.

### T5 — Explicit stability (the bipenalty differentiator) · Tier 1 · Zone-A · **headline**
**Claim:** with `-bipenalty` the global explicit run is **stable at the self-reported
`getExplicitCriticalTimeStep`** and the bound is **conservative** (it is *reported*, not left to
the per-element eigensolve that sees the massless coupling as `λ_max = 0`).
**Setup:** a slave node tied to a fixed host, mass penalty `m_p` from `-dtcr`/`-wcap`; integrate
with `CentralDifferenceLadruno` at `0.9×` and `1.1×` the reported `dt_cr`.
**Oracle (closed form):** `dt_cr = 2√(m_p / k_eff)` with `m_p = k_eff·(dt/2)²` (`-dtcr`) or
`m_p = k_eff/(β·ω_host)²` (`-wcap`); stable at `0.9×`, divergent at `1.1×`.
**Tolerance:** `dt_cr` self-report matches the closed form to `1e-9`; the stable/unstable split
brackets the bound. **Plus:** the `CriticalTimeStep` fold-in carries the embedded bound into
`ops.criticalTimeStep` (the embedded `1e-3` governs over a stiff host's `~0.04`).
**Status:** *partially covered* — `test_bipenalty_dtcr_self_report` checks the `m_p`/`dt_cr`
formula; **gap to add:** the actual `0.9×`/`1.1×` explicit-SDOF stable/unstable run + the
`criticalTimeStep` governance check (mirror the rebar's `test_bipenalty_governs_cfl_critical_step`).

### T6 — Cheap correctness adds (fold into the gate) · Tier 1 · Zone-A
Three near-free oracles that harden the core claim:
- **T6a — Partition-of-unity force split + equilibrium:** the force the tie pushes onto the
  hosts splits as `N_i·(slave force)` and balances the slave (`P = Bᵀt`). Exact, `1e-6`.
  *Covered:* `test_force_split_by_N`.
- **T6b — Objectivity under rigid host rotation:** the isotropic tie is frame-objective —
  a prescribed rigid rotation of the host produces **no** coupling force (no `-corot` needed).
  Exact. *Gap to add* (a clean rigid-rotation leg; the ADR §8 item 1 contract).
- **T6c — `-kt auto` conditioning:** `K_u = α·max|K_host(i,i)|` scales linearly with α and keeps
  a bounded condition number across host stiffness/mesh. *Covered:*
  `test_k_auto_scales_linearly_with_alpha`.

## 4. Pass/fail gate & status map

| Test | Tier | Zone | Oracle | Tol | Status |
|---|---|---|---|---|---|
| **T1** patch/consistency | 1 | A | exact affine interpolation | 1e-8 | weights covered; **affine-follow leg TODO** |
| **T2** pull-out (Abaqus) | 2 | B (+A surrogate) | springs-in-series / Abaqus deck | 5% | **TODO** |
| **T3** stress-free `g0` birth | 1 | A | zero force/gap at activation | 1e-3/1e-9 | **DONE (#214)** |
| **T4** beam-in-solid (LS-DYNA) | 2 | B (+A surrogate) | Winkler bed / LS-DYNA deck | 5–10% | **TODO** |
| **T5** explicit stability | 1 | A | `2√(m_p/k_eff)`; 0.9×/1.1× | 1e-9 | formula covered; **0.9/1.1 run + CFL governance TODO** |
| **T6a** force split | 1 | A | `N_i` split + equilibrium | 1e-6 | DONE |
| **T6b** objectivity | 1 | A | zero force under rigid rotation | 1e-6 | **TODO** |
| **T6c** `-kt auto` | 1 | A | linear-in-α, bounded κ | — | DONE |

**Gate to ship "validated v1":** all **Tier-1** rows (T1, T3, T5, T6a–c) green. Tier-2 (T2, T4)
is a credibility cross-check — surrogate legs in CI, commercial decks documented for later.

## 5. Explicitly OUT of v1 validation (experimental — do not claim)

> [!warning] Not validated, not in this battery (ADR §14.3)
> - **UR (`-rot`)** — `½ curl(u)` **spin** coupling, **not** lever-based moment transfer;
>   degenerates to a rigid element-wide spin on CST/TET4 (UR-4). A documented spin-restraint,
>   not a validated moment-transfer mechanism. Real moment transfer is a **future ADR**.
> - **UP (`-pressure`)** — niche poromechanics; experimental.
> - **D9 interface materials (`-matN/-matT*`) + `-corot`** — an **interface/contact** capability
>   adjacent to the `LadrunoContact` / [[LadrunoBondSlip_guide|LadrunoBondSlip]] lineage; kept as
>   an experimental convenience. Uncoupled per-direction uniaxials only **approximate** friction
>   (fixed `ElasticPP` slip, not `μ·N`).
>
> These have their own Zone-A *mechanics* tests today (they work as coded), but a **mechanics
> test is not a validation pass** — promotion to "supported" requires its own oracle-backed
> battery and the must-fix items in ADR §14.4 (notably the D9 `getInitialStiff` initial-tangent
> fix).

## 6. References
- **Abaqus** Analysis User's Guide — *Embedded elements* (`*EMBEDDED REGION`): translational
  constraint to the host, rotations free, geometric tolerance / activation.
- **LS-DYNA** Keyword Manual — `*CONSTRAINED_BEAM_IN_SOLID`, `*CONSTRAINED_LAGRANGE_IN_SOLID`
  (penalty coupling, `CDIR`, coupling stiffness); `*CONTACT_*` for rigorous friction.
- **Precedent:** `ASDEmbeddedNodeElement` (Petracca/ASDEA) — the implicit-only tri/tet penalty
  this element exceeds on hosts/conditioning/explicit.
- **Machinery reused:** [[20_ladruno_embedded_reinforcement_adr|ADR 20 §10]] (`-enforce` /
  bipenalty / `-cfl`), [[23_ladruno_embedded_node_adr|ADR 23]] (the element + the §14 scope).
- **Testbed:** `testbed/00_canonical_testbed` (Zone-A/B contract, tiered checks, tolerances).
