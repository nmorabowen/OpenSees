---
title: "WP-158 — Mortar friction tangent diagnosis at shared multi-pair nodes + the -consistanttan recipe"
project: Ladruno
type: ADR (amends the ADR-41 C2/C3 mortar tangent; status row in the ADR-48 capstone)
status: "Merged (#902, 3144e19ba, 2026-10-01) — stacks on #900 (WP-157); reshaped per review #902 (no new command, no wire bump)"
owner: nmora
related:
  - "[[41_ladruno_mortar_alm_contact_adr]] (C2.2 normal ALM, C3.2 friction tangent, C3.3 -consistanttan)"
  - "[[48_ladruno_contact_capstone_adr]] (status-of-record row)"
  - "[[157_mortar_friction_pair_state]] (§5: the two follow-ups this slice answers)"
  - "apeGmsh piles-validation ladder/R0_6_fork_slice/REVIEW_fable.md findings 1-2, ladder/R2c_prime/VERDICT.md (C-mu), ladder/R0_7_fork_slice/BUILDER.md"
tags: [adr, contact, mortar, friction, tangent, shared-node, pile, consistanttan, wp-158]
updated: 2026-10-01
---

# WP-158 — Mortar friction tangent diagnosis + `-consistanttan` recipe

> [!summary] The short version
> An FD check of the assembled mortar tangent against the residual finds two missing terms at
> slave nodes shared by several facet pairs on a creased interface. (1) The Coulomb pressure
> coupling `Csl = −μ εN t̂⊗n`. The default symmetric tangent drops it; the shipped opt-in
> `-consistanttan` already supplies it. With μ = 1 it is as large as the normal stiffness, which
> makes Newton linear (ratio 0.24). (2) The geometric terms `∂D/∂u`, `∂M/∂u` of the thin
> cross-crease pairs. They dominate a cohesion-only crease, and no shipped tangent has them.
> **What ships is the diagnosis and a recipe, no new code path:** for mortar Coulomb friction on a
> curved, faceted or non-matching interface, add `-consistanttan` and use a non-symmetric solver
> (D2). It stays opt-in; it is not made the default (D2). On the R2-C′ pile it makes C-mu at
> εT = εN/10 converge in every load case, with K within −1.2…−3.3 % of B. The finite-difference
> pair tangent that measured all this is **parked as a diagnostic oracle patch**
> (`contact_prototypes/adr158_fd_pair_tangent_oracle.patch`), not shipped (D1): it has no
> application win, it duplicates the residual with no sync guard, and on the faceted pile it makes
> Newton worse (D3). The binary is unchanged from #900 (WP-157), so the battery is hash-identical.

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

### D1 — no shipped FD tangent; the FD pair tangent is a parked oracle

R0.7 first shipped an opt-in `contact … -mortar -fdTangent [hRel]` (central FD of each pair's own
residual). Review #902 (recommendation A, option 2) asked to drop it from the binary, and this
ADR does: no `-fdTangent` token, no `MortarContact::fdTangent` slot, and the database wire format
stays at **v4** (the v4→v5 bump is gone, so every v4 database written by #900 builds still reads).

**Why not ship it.**

- No application benefit, measured: it stalls on the faceted pile at every `hRel` (D3), and the
  pile recipe is the existing `-consistanttan` (D2). Its only win is the cohesion-only crease
  (5, 3, 4, 5 iterations against 6, 4, 10, 15), a test geometry.
- A second copy of the residual with no sync guard: the evaluator (`mortarPairForceAt`, about 90
  lines) mirrors `getResidual` + `addMortarFriction`, and nothing pins that the two agree. The
  next change to the friction residual would silently de-sync it.
- The wire bump would be paid by every database for a slot that is 0.0 in every production deck.
- It is fragile as written (review #902 findings 2–4): `hRel` = 1e-5 stalls the cohesion ridge and
  1e-3 fails both roofs; it was silently inert on a 2D pair; a never-engaged node entering contact
  inside ±h takes `gT0 = 0`.

**The parked oracle.** `contact_prototypes/adr158_fd_pair_tangent_oracle.patch` is the full R0.7
C++ (`mortarActiveX` split, `mortarPairForceAt`, `addMortarTangFD`, the `-fdTangent` parser, the
v5 slot). It applies with `git apply` on the WP-157 tree and builds a **diagnostic** binary: the
tangent oracle any future analytic geometric linearization must match. Apply it only to a scratch
build; it bumps the database wire format, and findings 2–4 above are open in it. What it does:

- **Step.** `h = hRel · (longest current slave facet edge)`, `hRel` default 1e-7 (1e-7 and 1e-9
  work on both roofs; 1e-5 and 1e-3 do not).
- **Columns.** Every DOF of the pair, slave and master: `tang(:, c) −= ∂f/∂u_c` by central
  differences of a side-effect-free mirror of the residual (it reads `λ_N` and the WP-157 pair
  state and writes nothing).
- **Refusals.** When the kernel refuses a perturbed evaluation, the column falls back to the
  one-sided difference on the surviving side.
- The initial-stiffness path keeps the analytic SPD stick tangent; the tangent is non-symmetric.

**Rejected alternatives** (measured with the oracle).

- *Analytic geometric linearization* (`∂D, ∂M, ∂n, ∂ξ` of the clipped overlap). The ADR-41 C2
  deferral, still large. The oracle measures what it would buy (D3) before anyone writes it.
- *Freeze the base branch in the FD columns* (normal active set, stick/slip/τmax per node).
  Neutral to worse on the roof (cohesion: 5, 3, 7, 7 against 5, 3, 4, 5) and no help on the pile.
- *A minimum-overlap weight for the cross-crease pairs.* It puts a jump into the residual and
  changes converged answers. Not done.

### D2 — `-consistanttan` is the recipe for mortar Coulomb friction on curved / non-matching interfaces; it is not the default

**Recipe.** For a 3D mortar contact with μ > 0 on a curved, faceted or non-matching interface
(slave nodes shared by facet pairs whose tangent planes differ: a pile skin in a hole, a crease, a
faceted cylinder), add `-consistanttan` and use a non-symmetric solver: FullGeneral, UmfPack,
Mumps (unsymmetric) or Pardiso (this fork's default mtype 11, real unsymmetric). It restores the
Coulomb pressure coupling `Csl`, the dominant missing term (section 2): the μ = 0.2 solid roof
takes 3, 4, 4, 4 iterations against 8, 15, 15, 12. On the pile, use εT = εN/10 with it (D3).
It does not help a cohesion-only crease (there `∂cap/∂N = 0`, and the missing terms are the
geometric ones); there Newton is linear but converges, and a globalised algorithm
(`NewtonLineSearch`) keeps it robust across platforms (ADR-157 §4 (c)).

**Not the default** (review #902 recommendation B):

- Symmetric SOEs corrupt it: ProfileSPD, BandSPD, SparseSYM, `Pardiso -sym 1|2`, Mumps
  symmetric. The symmetric default friction tangent is a tested contract
  (`test_adr41_mortar_c3_2.py`, and the 2D twin in `test_adr85_contact2d_t2_friction.py`),
  recorded in ADR-39 P3.5 Q2.
- It changes the iteration path. The converged state is the same fixed point only to the solver
  tolerance, so every frictional battery hash would move, and "byte-identical" could not be
  claimed. The honest statement, "same converged state to `tol`", is what test (c) pins.
- The handler cannot auto-select: a `ConstraintHandler` sees the `AnalysisModel`, not the
  `LinearSOE`, and `LinearSOE` has no symmetry query. If a default ever flips, flip it per SOE
  after adding an `isSymmetric()` query, with the C3.2 ProfileSPD test retargeted.

The R2-C′ C-mu runs used neither flag; the pile ladder and apeGmsh's pile emit should carry
`-consistanttan` for Coulomb mortar contact.

### D3 — measured limits

**Roof tests** (`tests/test_adr158_mortar_consistanttan_multipair.py`; the last column was
measured with the oracle patch, D1). Iterations per step:

| Model | Shipped | `-consistanttan` | FD oracle |
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
| 0.1 | FD oracle | P0 and V stall at 0.07–0.2 | — |
| 1 | FD oracle | all four stall | — |

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

**New tests** (`tests/test_adr158_mortar_consistanttan_multipair.py`, 3 cases). The model is the
solid-backed creased roof with μ = 0.2 in full slip, under force control: the first binary μ > 0
multi-pair test (#900 review finding 2).

- (a) With `-consistanttan`, the vertical contact force equals the analytic rigid descent
  `2A·εN(δ+w)cos α·(cos α + μ sin α)` to 2e-4 (measured 3e-5), `Fx = 0`, and Newton takes ≤ 6
  iterations per step (3, 4, 4, 4). On d63f49750 (per-node path state, before WP-157) it fails
  with 3, 11, 19, 23 iterations.
- (a′) The shipped symmetric default on the same model is linear: more than 12 iterations in some
  step (8, 15, 15, 12 on Windows), **or** no convergence. The gate is `(not ok) or max(its) > 12`,
  because that baseline is platform-sensitive (review #902 finding 1).
- (c) When the default also converges, it reaches the same state as `-consistanttan` to 1e-6 of
  the largest displacement (skipped, not failed, if the default does not converge).

There is no refusal test: this slice adds no command surface.

**Platforms.** Windows/MSVC and Linux/gcc (Esmeralda, a CI-like build with the bundled reference
BLAS/LAPACK); see the PR body and piles-validation `ladder/R0_7_fork_slice/RESHAPE.md`.

**Contact battery.** This slice changes no C++ file, so the binary is the WP-157 (#900) code. The
60-file WP-157 battery is run on the #900 binary and on this head with the byte-dump plugin
(`contact_prototypes/bytedump_plugin.py`, which hashes every node's displacement and the time
after each `analyze`); the counts and every snapshot hash are identical (numbers in the PR body).

## 5. Follow-ups

- **εT = εN on the faceted pile.** Measure it on M2/M3 with `-consistanttan`; finding C′-2
  shrinks about 400× from M1 to M3.
- **An analytic geometric linearization** only if a smooth-crease production case needs it. The
  parked FD patch (D1) is its oracle; before relying on it, fix its `gT0` fallback (review #902
  finding 4) and add one FD-vs-analytic pin on a single-pair flat stick case.
- **apeGmsh pile emit:** carry `-consistanttan` for Coulomb mortar contact (D2).
- **The element-less μ > 0 roof two-cycle.** It is a per-facet active-set kink at lift-off
  nodes, not a tangent issue.
