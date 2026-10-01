---
title: "WP-158 — Mortar friction tangent at shared multi-pair nodes: FD diagnosis and an opt-in FD pair tangent"
project: Ladruno
type: ADR (amends the ADR-41 C2/C3 mortar tangent; status row in the ADR-48 capstone)
status: "PR open (#902, not merged) — stacks on #900 (WP-157)"
owner: nmora
related:
  - "[[41_ladruno_mortar_alm_contact_adr]] (C2.2 normal ALM, C3.2 friction tangent, C3.3 -consistanttan)"
  - "[[48_ladruno_contact_capstone_adr]] (status-of-record row)"
  - "[[157_mortar_friction_pair_state]] (§5: the two follow-ups this slice answers)"
  - "apeGmsh piles-validation ladder/R0_6_fork_slice/REVIEW_fable.md findings 1-2, ladder/R2c_prime/VERDICT.md (C-mu), ladder/R0_7_fork_slice/BUILDER.md"
tags: [adr, contact, mortar, friction, tangent, shared-node, pile, wp-158]
updated: 2026-10-01
---

# WP-158 — Mortar friction tangent at shared multi-pair nodes

> [!summary] The short version
> An FD check of the assembled mortar tangent against the residual finds two missing terms at
> slave nodes shared by several facet pairs on a creased interface. (1) The Coulomb pressure
> coupling `Csl = −μ εN t̂⊗n`. The default symmetric tangent drops it; `-consistanttan` already
> supplies it. With μ = 1 it is as large as the normal stiffness, which makes Newton linear
> (ratio 0.24). (2) The geometric terms `∂D/∂u`, `∂M/∂u` of the thin cross-crease pairs. They
> dominate a cohesion-only crease. This slice adds `-fdTangent [hRel]`, an opt-in tangent that
> takes a central finite difference of each pair's own residual. It carries both terms, and on
> the roof tests Newton takes 3–5 iterations instead of 10–15. The residual is unchanged and the
> option is off by default: the 60-file battery is hash-identical on 145 846 analyze snapshots.
> The R2-C′ pile does **not** need it. The fix there is the existing `-consistanttan` on a
> non-symmetric solver (Pardiso is mtype 11). With it, C-mu at εT = εN/10 converges in every load
> case, with K within −1.2…−3.3 % of B. On the faceted pile the geometric FD terms make Newton
> worse: the clip topology is non-smooth there (section 3).

## 1. Question

Review #900 finding 1: on the force-controlled shared-ridge roof, Newton is linear (10 and 14
iterations, ratio ≈ 0.27), while the split-ridge twin is quadratic. Finding 2: no binary test
covers μ > 0 at multi-pair nodes. On the R2-C′ pile, Coulomb μ = 1 diverges at εT = εN, fails at
εT = εN/10, and converges only at εT = εN/100, and then it is 5–13 % too soft. Is the assembled
tangent inconsistent with the residual at shared nodes, and which term is responsible?

## 2. FD diagnosis

**Probe** (`contact_prototypes/probe_adr158_mortar_tangent_fd.py`). The tangent `K` comes from
`printA -ret` and the residual `R(u)` from `printB -ret`. The trial state is set with
`setNodeDisp … -commit`, a node-level commit that leaves the contact Domain's path state
untouched. Two probe traps cost a day; both are now in LEDGER_quirks:

- `setNodeDisp` without `-commit` starts from the committed vector, so setting the DOFs one by
  one resets the earlier ones.
- `printA -ret` returns the column-major `Matrix` buffer, so `reshape(n, n)` gives `Kᵀ`.

Springs and bricks are not `update()`d by `setNodeDisp`, so their exact linear part is added
analytically. A Python Newton loop drives the same residual with `K`, with the FD `K`, or with
a numpy re-assembly. It replicates the OpenSees iterates digit for digit (`|K_os − K_np| ≤ 2e-14`).
The numpy pair replica reuses the C1 oracle's `mortar_pair`.

**Element-less roof, cohesion only** (the review model, shared ridge). The per-pair FD
(numpy, at the step-3 slip iterate) gives these errors:

| Pair class | Relative tangent error |
|---|---|
| Main flank pairs | ≤ 0.9 % (the frozen-D,M penalty term is accurate) |
| Cross-crease pairs (a slave facet clipped against the other flank's master facet; overlap ≈ 1e-4 of a facet) | ≈ 100 % |

The cross-crease overlap width grows linearly with the ridge penetration, so `p·∂D/∂u` is the
same order as the material term. Assembled, the spectral radius of `I − K⁻¹K_FD` is 0.033
with a shared ridge and 0.0016 with a split ridge. Newton on the FD tangent converges in
4 iterations; the shipped tangent needs 8. `-consistanttan` changes nothing, because with
cohesion only `∂cap/∂N = 0`.

**Solid-backed roof, μ = 1** (a brick layer under the slave, the realistic case: a pile skin is
a solid). The shipped tangent converges linearly with ratio 0.244 (15 iterations). Two other
tangents both converge in 4 iterations:

- the numpy assembly of the analytic tangent **with** `Csl` and **without** any geometric term;
- the FD tangent.

So `Csl` is the dominant missing term. The per-pair numbers agree: slipping main pairs have a
73 % error without `Csl` and 1e-3 with it. In OpenSees, `-consistanttan` gives 3–4 iterations
per step against 8–15 without it.

**Element-less roof, μ > 0.** Newton two-cycles even with the FD tangent. The cause is the
per-facet normal active set at nodes that lift off (`p ≈ 0`), a residual kink and not a tangent
error. An element-less slave is not a model of a pile, so the tests use the solid-backed roof.

**Answer to review finding 1.** Two terms are missing: the Coulomb `Csl` (when μ > 0) and the
cross-crease geometric terms. The suspected "per-pair `N_I` against nodal `λ_N`" mismatch is not
present. `λ_N` is a committed constant within a step, and the per-pair pressure enters the
residual and the tangent identically.

## 3. Decisions

### D1 — `-fdTangent [hRel]`, opt-in, 3D mortar (non-tie) only

`LadrunoContactFE::addMortarTangFD` replaces the analytic pair tangent (normal ALM + C3.2/C3.3
friction) with a central difference of this pair's own static residual. That residual comes from
`mortarPairForceAt`, a side-effect-free mirror of the `getResidual` MORTAR branch. It reads
`λ_N` and the WP-157 pair friction state and writes nothing. A pair that has not engaged takes
`gT0` at the base configuration, so every column shares the origin the residual captures there.

- **Step.** `h = hRel · (longest current slave facet edge)`, with `hRel` defaulting to 1e-7.
- **Columns.** Every DOF of the pair, slave and master.
- **Refusals.** When the kernel refuses a perturbed evaluation (the clip's convexity, sliver or
  back-map guards flipping within ±h), the column falls back to the one-sided difference on the
  surviving side. Differencing across a refusal would put `f/h` into the tangent.
- **Unchanged paths.** The initial-stiffness path (`addKiToTang`) keeps the analytic SPD stick
  tangent. The viscous term stays in `C`. 2D and `-tie` are out of scope; the command surface
  refuses `-tie` and non-mortar contacts.
- **Solver.** The tangent is non-symmetric and needs a non-symmetric solver. The command warns.

**Why opt-in, not default.**

- It is non-symmetric, so a symmetric SOE would corrupt it (the `-consistanttan` rule).
- It changes the iteration path, so converged results differ at round-off. Single-pair byte
  identity of results would not hold.
- It costs `2·3·(npsS+npsM)` = 48 pair integrations per tangent per active pair.
- On the faceted pile it is worse than the analytic tangent (D3).

Off, nothing in the shipped path runs. `mortarActive` was split into `mortarActiveX` with the same
arithmetic, and the battery is hash-identical.

**Rejected alternatives.**

- *Analytic geometric linearization* (`∂D, ∂M, ∂n, ∂ξ` of the clipped overlap). This is the
  ADR-41 C2 deferral, still large. The FD tangent measures how much it would buy (D3) before
  anyone writes it.
- *Freeze the base branch in the FD columns* (normal active set, stick/slip/τmax per node).
  Implemented and measured. On the roof it was neutral to worse (cohesion: 5, 3, 7, 7 against
  5, 3, 4, 5 for the central difference) and it did not help the pile. Dropped.
- *A minimum-overlap weight for the cross-crease pairs.* Dropping a pair below a threshold puts a
  jump into the residual and changes converged answers. Not done.
- *`-consistanttan` by default.* It is non-symmetric, which the default symmetric-solver contract
  forbids (ADR-39 P3.5 Q2).

### D2 — `-consistanttan` is the pile recipe for μ > 0

Pardiso in this fork is mtype 11 (real unsymmetric), and UmfPack, Mumps and FullGeneral are also
non-symmetric, so `-consistanttan` is safe on every solver the pile ladder uses. The R2-C′ C-mu
runs used neither flag.

### D3 — measured limits

**Roof tests** (`tests/test_adr158_mortar_fd_tangent.py`). Iterations per step:

| Model | Shipped | `-consistanttan` | `-fdTangent` |
|---|---|---|---|
| Solid roof, μ = 0.2 | 8, 15, 15, 12 | 3, 4, 4, 4 | 3, 3, 3, 3 |
| Element-less shared-ridge cohesion roof | 6, 4, 10, 15 | 6, 4, 10, 15 | 5, 3, 4, 5 |

**R2-C′ pile, P5-M1, μ = 1.** The flags are `-adjust -gapOffset -0.001 -augment never -maxGap 0.1`
with Pardiso. K is relative to B at the same level.

| εT/εN | Tangent | P0 / H / M / V (iterations; final norm) | K vs B |
|---|---|---|---|
| 0.1 | shipped (R2-C′) | 23 ok / H fails / M round-off / V diverges | V −89 % |
| 0.1 | `-consistanttan` | 16 / 30 (1.3e-8) / 25 / 17, all converged | HH −2.05, HM −1.20, MM −1.70, free −3.06, V −3.34 % |
| 1 | `-consistanttan` | 16 / diverges / diverges / 12 | — |
| 0.1 | `-fdTangent` | P0 and V stall at 0.07–0.2 | — |
| 1 | `-fdTangent` | all four stall | — |

The pile is unusable at εT = εN with either tangent. For reference, the shipped εT = εN/100 runs
were 5–13 % soft against B.

**Why the FD tangent hurts on the pile.** The hole (n + 4 segments) and the skin (n) share
vertices at the meshed configuration. The clipped overlap of a pair therefore changes topology
inside ±h, and the one-sided slopes differ by about 100 % (measured: `|∂f/∂u|` ≈ 1e6 against a
penalty scale `εN·a` ≈ 3e5, under the 1 mm prestress `p` = 1e4 kPa). The FD tangent's first
iterate norm grows like `h^-1/2` (3.8e2, 3.8e3, 5.9e4 for `hRel` = 1e-5, 1e-7, 1e-9), and Newton
stalls at 0.1–5 for every `hRel`. The faceted mortar residual is genuinely non-smooth along that
path, and the analytic tangent, which drops the geometric terms, is the better-conditioned
choice there. The remaining εT = εN divergence is not a tangent defect this slice can fix. It
belongs with the R2-C′ finding C′-2 (polygon mismatch) and the mesh rule (≥ M2 around the pile).

## 4. Verification

**New tests** (`tests/test_adr158_mortar_fd_tangent.py`, 7 cases):

- (a) The solid-backed creased roof with μ = 0.2 in full slip, the first binary μ > 0 multi-pair
  test. The vertical contact force equals the analytic rigid descent
  `2A·εN(δ+w)cos α·(cos α + μ sin α)` to 2e-4 (measured 3e-5), and `Fx = 0`. Newton takes ≤ 6
  iterations per step with `-fdTangent` and with `-consistanttan`. A twin pins that the shipped
  symmetric default needs more than 12.
- (b) The review's shared-ridge cohesion roof: ≤ 6 iterations per step with `-fdTangent`.
- (c) The FD tangent converges to the same state (1e-6 of the largest displacement, for μ and for
  cohesion), and the refusals.

Result: 7/7 on the new binary. On fc75db7f3 (the #900 head) the four `-fdTangent` cases fail,
because the flag does not exist there. ADR-157 stays 12/12.

**Contact battery.** The 60 files of the WP-157 battery give 379 passed on fc75db7f3 and 379 on
this head. The byte-dump plugin hashes every node's displacement and the time after each
`analyze`. All 145 846 snapshots in 279 tests are identical. The plugin is checked in as
`contact_prototypes/bytedump_plugin.py`, which answers #900 review finding 4.

## 5. Follow-ups

- **εT = εN on the faceted pile.** Measure it on M2/M3 with `-consistanttan`; finding C′-2
  shrinks about 400× from M1 to M3.
- **An analytic geometric linearization** only if a smooth-crease production case needs it. The
  FD tangent is its oracle.
- **The element-less μ > 0 roof two-cycle.** It is a per-facet active-set kink at lift-off
  nodes, not a tangent issue.
