---
title: "LadrunoKinematicCoupling — RBE2 / Kinematic-Coupling Element"
project: Ladruno
type: reference / user guide
status: v1 shipped (rigid driver, U + rotation, -dof component select, penalty/AL, derived K_r, bipenalty default-off, g0 birth, 2D+3D); built + 16/16 Zone-A + 4-lens code review; MERGED to ladruno (PR #221, 2026-06-09). v1.1 (WP-101, PR #839): the within-step route is the ADR-41 held-load augmentation sweep with the DEFAULT `-alUpdate commit` (measured to close the gate for every algorithm × integrator); `-alUpdate iter` is opt-in and guarded to full Newton + LoadControl; committed-λ revert; K_t selection rule + conditioning guard; battery 52/52
element: element LadrunoKinematicCoupling
classTag: 33012 (element)
related:
  - "[[29_ladruno_kinematic_coupling_rbe2_adr]]"
  - "[[24_ladruno_coupling_constraints_adr]]"
  - "[[LadrunoDistributingCoupling_guide]]"
  - "[[LadrunoEmbeddedNode_guide]]"
  - "[[ndf_and_mixed_models_guide]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - guide
  - element
  - rbe2
  - kinematic-coupling
  - rigid-body
  - rigid-link
  - ndf-mismatch
  - moment-transfer
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
updated: 2026-09-14
---

# LadrunoKinematicCoupling — RBE2 / kinematic coupling

`element LadrunoKinematicCoupling` (classTag **33012**) makes a **reference node R**
**rigidly drive** a set of **slave nodes** `{S_i}`: every tied slave DOF follows R's
rigid-body motion,

```
u_i = u_R + θ_R × d_i          (slave translations follow R, with the moment arm d_i)
θ_i = θ_R                       (slave rotations follow R, where the slave carries them and the DOF is tied)
```

with `d_i = x_i − x_R`. A load or prescribed displacement at R is transmitted to the set
**with rigid stiffness** — the region moves as one rigid body. It is the OpenSees
realization of **Nastran RBE2** / **Abaqus `*COUPLING, kinematic`** / **LS-DYNA
`*CONSTRAINED_NODAL_RIGID_BODY`**, and it generalizes the classic `rigidLink`/`rigidDiaphragm`
to an arbitrary node set, an offset reference, selectable tied components, finite penalty
(or augmented-Lagrangian) enforcement, and explicit-dynamics safety.

This is the single user/developer reference. It covers the **theory** (the rigid tie, the
gaps, the transport sign, the force transfer), the **capabilities** (the option grammar),
the **explicit-safety** machinery, the **implementation** (where the math lives,
serialization, responses), the **validation** battery, and **use cases**. The design record
is [[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]]; the family scoping is
[[24_ladruno_coupling_constraints_adr|ADR 24]].

> [!important] RBE2 vs RBE3 — the single most common wiring mistake
> In **RBE2 (this element)** the *single* node R is the **master**: it rigidly drives the
> *many* slaves and **adds stiffness** to them (the region becomes rigid). In **RBE3
> ([[LadrunoDistributingCoupling_guide|distributing coupling]])** the *single* node is the
> **dependent**, interpolated from the *many* independents, and the set is **not** stiffened.
> The "one vs many" roles **flip**. Reach for **RBE2** when the region must genuinely move as
> a **rigid body** (a loading platen, a rigid offset, a rigid connection block); reach for
> **RBE3** to **introduce/transmit a load or BC at a point while the region stays flexible**.
> Using RBE2 where RBE3 is wanted **over-stiffens** the patch (artificial local rigidity);
> using RBE3 where RBE2 is wanted leaves the region too soft.

---

## 1. Theory — the rigid coupling

### 1.1 The kinematic constraint (rigid-body driver)

R is the **master**. Each tied slave DOF is constrained to R's rigid-body field. With
`d_i = x_i − x_R` the moment arm from R to slave `i`, the per-slave constraint gaps are

```
g_i^t = u_i − ( u_R + θ_R × d_i )        (translation rows — the transport block)
g_i^r = θ_i − θ_R                         (rotation rows — identity blocks, where tied)
```

`g_i^t = 0` says the slave translation is R's translation **plus** the rigid lever
`θ_R × d_i`; `g_i^r = 0` says the slave rotation equals R's rotation. Unlike RBE3 there is
**no weighted centroid, no position-inertia `I_c`, no eigensolve, no least-squares fit** —
the tie is a *direct* per-slave rigid link. This makes RBE2 substantially simpler than its
distributing sibling (it reuses RBE3's DOF-scatter / AL / bipenalty / `g0` scaffold but
replaces the geometry resolve and the `B` build).

> [!important] Transport sign flip vs RBE3 (don't copy the line verbatim)
> RBE2 carries `−θ_R × d_i` in `g_i^t`, so the transport block is
> `∂g_i^t/∂θ_R = +[d_i]_×`, the **negative** of the operator RBE3 uses (RBE3 carries
> `+θ_R × (x_c − x_R)`). The `B` assembly is `B += −transOp`. Copying RBE3's `buildB` line
> as-is reverses every offset moment; the validation suite asserts the **sign** of the
> transport couple (test 2), not just its magnitude.

### 1.2 Selectable tied components (`-dof`) and the ragged gap layout

`-dof c1 … cK` chooses which slave components are dependent (Nastran's RBE2 dependent-DOF
list): components `1..ndm` are translations, `ndm+1..ndm+nrot` are rotations. Default = every
DOF the slave possesses. One element therefore subsumes:

- **u–p slaves (ndf 4 in 3D, ndf 3 in 2D) MUST carry an explicit `-dof 1 2 3`** (or `1 2`): the
  default list is refused for them since 2026-09-07, because without the refusal the pressure
  DOF (node index `ndm`) was tied to θ_R as if it were a rotation — silently, the "lacks DOF"
  skip message being reachable only with an explicit list. apeGmsh's `kinematic_coupling`
  guards the same case on its side (its PR #1100).
- **translation-only** tie (`-dof 1 2 3` in 3D) — slaves follow R's translation + transport, but spin freely;
- **translation + transport** of an offset reference;
- **full rigid** tie (default on 6-DOF slaves) — translations *and* rotations driven.

Because slaves can be a **mix of 3-DOF and 6-DOF** nodes, the number of tied DOFs varies
per slave: the gap vector is **ragged**. The element resolves the layout **once** at
`setDomain` into three parallel index arrays — `gapNode` (which slave), `gapDof`
(node-local component), `gapIsRot` (translation vs rotation row) — so a flat uniform stride
can never mis-index the moment reaction silently.

### 1.3 Penalty formulation, residual, and consistent tangent

Penalize the gaps with a per-row penalty `D_i = diag(K_t I, K_r I)` (translation rows get
`K_t`, rotation rows get `K_r`):

```
Π_p = ½ Σ ( K_t |g_i^t|² + K_r |g_i^r|² )         (penalty potential)
t    = D g                                          (coupling traction, per row)
P    = Bᵀ t                                         (resisting force)
K    = Σ Bᵢᵀ D_i Bᵢ = BᵀDB                          (consistent tangent)
```

`g = B u − g0` is **linear** in the DOFs and `B` is **constant** (geometry only), so
`K = BᵀDB` is the **exact** penalty Hessian — symmetric, PSD, state-independent (hence
`getInitialStiff ≡ getTangentStiff`). `B` (size `nGap × nDOF`) is built once at `setDomain`.
Its non-zero blocks for tied slave `i`:

| ∂/∂ | `g_i^t` (translation rows) | `g_i^r` (rotation rows) |
|---|---|---|
| `u_R`  | `−I`              | `0` |
| `θ_R`  | `+[d_i]_×` (transport, sign-flipped vs RBE3) | `−I` |
| `u_i`  | `+I`              | `0` |
| `θ_i`  | `0`               | `+I` |

The force transfer is the work-conjugate transpose: a force at R balances against the
slave reactions (`Σ reaction = −F`), and a **moment** at R enters even a **3-DOF slave
face** as a **self-equilibrated force couple** (`Σ f_i = 0`, `Σ d_i × f_i = −M`) — the same
ndf-mismatch moment-transfer driver as RBE3, but here the patch is held **rigid** rather
than left free to deform.

> [!note] Reference configuration only (finite-rotation boundary)
> In the reference configuration `B` is constant and `K` is exact. Under a **finite** rigid
> rotation the lever `θ_R × d_i` is the first-order (small-rotation) form; a finite-rotation
> update of `d_i` (with its geometric-stiffness term) is deferred (ADR 29), shared with
> RBE3. Fine for the small-relative-rotation regime of rigid offsets, platens, and
> connection blocks.

### 1.4 Degeneracy & refusals

RBE2 has no `I_c` to go rank-deficient, so its degeneracy handling is about **ill-posed
ties**, not rotation-axis dropping:

| Situation | Handling |
|---|---|
| **Self-tie** (a slave tag `==` the reference node) | **refused** — but *inertly*: `valid=false`, `nGap=0`, matrices allocated zero (it does **not** early-return before `allocate`, which would null-`K` crash) |
| **Duplicate slave** | refused inertly (same path) |
| **All-coincident slaves** (`d_i = 0`) with a rotation tied | the default `K_r = K_t·ℓ²` is **floored** strictly positive (a zero `ℓ²` can't leave a tied rotation unpenalized); translation-only + bipenalty reports a **finite** `dtcr` (no `2√(I_p/0) = +Inf` trap) |
| Reference node with too few DOFs for the tied rotations | refused at `setDomain` |

Health check: `eleResponse $tag tiedDOFs` returns the total tied-DOF count (`Σ` over slaves
of their tied components). `0` means the element went inert (a refusal fired); a full rigid
3D tie of `N` 6-DOF slaves should read `6N`.

### 1.5 No damping; mass via bipenalty (default OFF)

A pure coupling carries **no physical damping** (`getDamp ≡ 0`; Rayleigh factors are refused
so a `βK` can't shrink the explicit step). Unlike RBE3 the reference node is **often a real
massed node** (a platen, a footing, an equipment block), so the mass-penalty machinery
**defaults OFF** — see §5.

---

## 2. Capabilities — the option grammar

```
element LadrunoKinematicCoupling $tag $refNode $N $s1 ... $sN
        [-dof $c1 ... $cK]               # dependent components on each slave (default: all the slave has)
        [-k {$Kt | auto}] [-kAlpha $a] [-host $eleTag]   # auto needs a representative -host
        [-kr $Kr]                        # rotational-tie penalty (default DERIVED K_t·ℓ²)
        [-enforce {penalty | al}]        # default penalty
        [-alUpdate {commit | iter}]      # AL Uzawa cadence; default commit (see §4.2)
        [-bipenalty {-dtcr $dt | -wcap $beta}]   # default OFF (R is often massed)
        [-absolute]                      # opt out of g0 initial-gap (offset) capture
```

| Token | Meaning | Default |
|---|---|---|
| `$refNode` | the **master** reference node (ndf ≥ ndm+nrot: **6** in 3D, **3** in 2D) | — |
| `$N $s1..sN` | count + tags of the **slave** nodes (ndf ≥ ndm; may mix 3-/6-DOF) | — |
| `-dof $c1..cK` | dependent components per slave: `1..ndm` trans, `ndm+1..ndm+nrot` rot | all DOFs the slave has — **only for slaves of ndf `ndm` or `ndm+nrot`**; any other ndf (an ndf-4 u-p node, say) is REFUSED at the parser without an explicit `-dof`, because the default would tie node DOF `ndm+1` (its pressure) to a master rotation (TIMs 2026-09-07 no-ask 1) |
| `-k $Kt` / `-k auto` | translational penalty; `auto` = `kAlpha·max\|K_host(i,i)\|` (needs `-host`) | `1e12` |
| `-kAlpha $a` | multiplier for `-k auto` | `1e3` |
| `-host $eleTag` | one **representative** slave-side element, used only to scale `-k auto` / `-wcap` | none |
| `-kr $Kr` | rotational-tie penalty; omit to **derive** `K_r = K_t · ℓ²` (floored) | derived |
| `-enforce` | `penalty` or `al` (augmented Lagrangian, implicit) | `penalty` |
| `-alUpdate` | where the AL Uzawa recursion advances: `commit` = once per **committed step** (λ frozen inside a step ⇒ safe with **every** algorithm and integrator, and the cadence the §4.3 augmentation sweep turns into a true outer Uzawa loop); `iter` = once per **equilibrium iteration** — **opt-in, expert, refused outside full Newton + LoadControl** (§4.4). Ignored (with a warning) without `-enforce al` | `commit` |
| `-bipenalty` | explicit critical-step control (see §5); **needs** `-dtcr`/`-wcap` budget | **off** |
| `-absolute` | keep the **absolute** tie (skip stress-free `g0` capture); `-noInitGap` alias | off (capture on) |

Build the reference node with the rotational DOFs it needs — `ndf 6` in 3D, `ndf 3` in 2D
(per-node `-ndf`, see [[ndf_and_mixed_models_guide]] §1.3); the element **refuses** a
reference node with too few DOFs. Slave nodes need at least translations (`ndf ≥ ndm`); a
slave needs `ndf ≥ ndm+nrot` for its rotations to be tieable.

---

## 3. Penalty & auto-scaling (`-k`, `-kr`)

`-k $Kt` sets the translational penalty directly (default `1e12`). `-k auto` scales it from a
representative slave-side host element's initial-stiffness diagonal,
`K_t = kAlpha · max|K_host(i,i)|` — mesh/material-independent conditioning — but it requires
**`-host`** because an RBE2 set has **no single host element**. Without `-host`, use a numeric
`-k`.

> [!important] What `-k auto -host` is FOR — and what it is NOT
> `-k auto` is a **conditioning** control, inherited from
> [[LadrunoEmbeddedNode_guide]] §3 where it was introduced: `K_u = a·max_i|K_host(i,i)|`
> off `host->getInitialStiff()`, the max-absolute diagonal (which already carries `~E·l_ch`
> units), default `a = 1e3`. Its stated job there is to make the penalty **mesh- and
> material-independent in conditioning** — "a coarse stiff host and a fine soft host both get
> a well-scaled tie" — and the same guide warns explicitly against passing ASD's `1e18` as
> `a` ("not `E`-scaled — condition-number blow-up"). That is the property you want in an
> **explicit / dynamic** run, where an over-stiff tie collapses `dt_cr`, and in an
> **embedded** use where the tie must not dominate the host it lives inside. [[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]] §6
> carries it over to RBE2 with one restriction: "**`-k auto` is undefined for a node set**
> (no host element). Default to numeric `-k`; `-k auto` only with an explicit representative
> `-host` element among the slaves' parents" — hence the parse-time refusal of bare `-k auto`.
>
> **It is therefore NOT a rigidity setting, and it cannot hold a rigid footing.** By
> construction it pins `K_t` to the *host's* order of magnitude, which is exactly the order at
> which the penalty tie is *comparable* to the structure it is tying — the regime where the
> residual gap is largest. Measured: on the TIMs strip's elastic rigidity gate (master-node
> push vs a direct footprint push, tol `1e-6`) `-k auto -host <soil element>` left
> **`1.4e-4`** of the push un-transmitted, versus `5.7e-8` for the flat default `1e12`. The
> Zone-A gate (§9, a 2×2×2 block of the same soil) reproduces the mechanism: `-k auto`
> resolves to `K_t ≈ 7.9e6` and leaves `2.1e-3`. If the footing must genuinely be rigid,
> either raise `-k` (§3.1) or — better — keep the conditioning-friendly `K_t` and add
> `-enforce al` (§4.2), which converges the constraint *at* that moderate stiffness.

### 3.1 Choosing `K_t` — the rigidity/conditioning trade-off

Two effects run in opposite directions, and both are measurable:

**Rigidity error falls as `1/K_t`.** With the master's motion prescribed, the only thing
resisting the tie is the host, so the residual gap is the tie force divided by the penalty:

```
err  ≡  max_i |g_i| / |push|  ≈  c / K_t
```

`c` is a property of the model (host stiffness × footprint), not a tolerance. On the Zone-A
gate (2×2×2 `stdBrick`, `E = 45 000 kPa`, `ν = 0.3`, `B = 1.5 m`) `c = 1.66e4`, constant to
five digits over `K_t = 1e6 … 1e12`. The TIMs strip shows the same law with its own `c`:
`K_t = 5e9 → 6.7e-7`, and `K_t = 1e12 → 5.7e-8` — the last point is no longer on the line
because it has hit the round-off floor of the assembled solve, which is the first sign you
have gone too far.

**Conditioning degrades linearly in `K_t`.** The tie contributes `K_t·BᵀB` to the global
matrix; once that block is many orders above everything else, the matrix is numerically
singular even though it is formally non-singular. Measured on the TIMs strip
(**101 583 DOF**): at `K_t = 1e12` **Pardiso limped on perturbed pivots and then died**, and
**SuperLU failed its first factorisation**; at `K_t = 5e9` the *same leg* ran to a clean
plateau. The flat default `1e12` is a small-fixture default, not a production one.

**The rule.** Express the band relative to the host element's diagonal stiffness scale
`k_host = max_i |K_host(i,i)|` (what `-k auto` reads):

| `K_t / k_host` | What you get |
|---|---|
| `≲ 1e1` | the "rigid" patch is visibly soft — the tie is a spring, not a constraint |
| **`1e2 … 1e4`** | **recommended.** Rigidity error `1e-3 … 1e-5` of the push, conditioning untouched. Pair with `-enforce al` (§4.2) to take the error to round-off *without* leaving the band |
| `1e4 … 1e6` | works, tightening; watch the solver's pivot warnings on large models |
| `> 1e6` | **warned.** Rigidity stops improving (round-off floor) while the factorisation degrades. This is the regime that killed the 101 583-DOF leg |

**Worked example — the TIMs strip.** Soil `E = 45 000 kPa`, footing `B = 1.5 m`, elements
~`0.2 m`, so `k_host ~ E·l_ch ~ 1e4`. The band is then `K_t ≈ 1e6 … 1e8`; `5e9` (`~5e5·k_host`)
is at the tight end and was the value that carried the leg to a clean plateau at `6.7e-7`;
`1e12` (`~1e8·k_host`) is two decades past the warning and is the value that made the system
near-singular. With `-enforce al` the *same* strip would be run at `~1e6–1e7` and reach the
tolerance exactly (§4.2) — that is the recommended production setting.

> [!note] Why the warning is not at parse time
> With a numeric `-k` **and** a `-host` named, the element warns once when
> `K_t > 1e6 · k_host`. The check lives in `resolveAutoKt()` (first use, `ktResolved`-guarded)
> rather than in the parser because it needs `host->getInitialStiff()`: at parse time the host
> element may not be in the domain yet (element order in a deck is free) and, if it is, may not
> have had `setDomain()` called — calling `getInitialStiff()` there is a null-node dereference
> waiting to happen. `resolveAutoKt()` is the first point where the value genuinely exists, and
> it still runs before the first factorisation. **No `-host`, no warning** — there is nothing
> to compare against, which is itself a reason to name a representative `-host` even when you
> are passing a numeric `-k`.

The **rotational** penalty `K_r` is, by default, **derived** from `K_t` and the geometry:

```
K_r = K_t · ℓ²,   ℓ² = (lever-weighted length scale, floored to max_i |d_i|²)
```

This is mandatory, not cosmetic: a translation gap is a length and a rotation gap is a
rotation (radians), so a single penalty would be off by `ℓ²` and wreck the conditioning (and
the explicit step). For the flat-face fixture (`a = 1`, all `|d_i|² = 2`), `ℓ² = 2` so
`K_r = 2·K_t` (validation test 14). Override with `-kr $Kr` only with a specific reason.

---

## 4. Enforcement strategies (`-enforce`)

### 4.1 `penalty` (default)
The gaps are driven toward zero by `K_t`/`K_r`; the **force/moment transfer is exact for any
penalty** (§1.3), and the residual *kinematic* gap is `O(1/K)` — precisely `c/K_t`, see §3.1.
On a small fixture the default `1e12` is stiff enough that the region is rigid to round-off;
on a production-size model it is **not** a safe default (§3.1, the 101 583-DOF case).

### 4.2 `al` (augmented Lagrangian)

Adds per-gap-row multipliers `λ` (size `nGap`) carried on the **same** tangent, with the
traction `t = D g + λ`. The Uzawa recursion is `λ ← λ + D g`; **where it advances** is
`-alUpdate`, and that choice is not cosmetic — it decides whether the residual is still a
function of the displacements.

**`-alUpdate commit` (DEFAULT).** One update per **committed step**. Inside a step `λ` is
**frozen**, so the residual `r(u) = f − S u − Bᵀ(λ + D B u)` is a genuine function of `u` and
the tangent `S + BᵀDB` is its exact Jacobian. Every algorithm and every integrator works,
because that is the contract each of them assumes. Across steps it is a first-order Uzawa: a
multi-step push tightens the tie step by step (measured on the §9 gate, 5 DisplacementControl
steps: `1.63e-2 → 3.33e-3`), but a **single** step is exactly the penalty answer.

**To close the constraint inside one step, do not change the cadence — add an outer loop.**
See §4.3.

### 4.3 Closing the constraint **within** a step — the held-load augmentation sweep

The consistent way to drive `g → 0` inside a step is an **outer Uzawa loop whose inner solve
sees a fixed residual**: converge at frozen `λ`, update `λ`, repeat. The fork already ships
that mechanism — [[41_ladruno_mortar_contact_adr|ADR-41]] D1's held-load augmentation sweep —
and `LadrunoKinematicCoupling` participates in it **with no new hook**: `Domain::commit()`
still runs every element's `commitState()` during the sweep (`Domain::contactAugmenting` only
suppresses the recorder loop and the `commitTag` bump), and with `-alUpdate commit` that
`commitState` **is** the outer update.

```python
ops.element('LadrunoKinematicCoupling', 1, ref, n, *skin,
            '-dof', 1, 2, 3, '-k', 1.0e6, '-enforce', 'al')      # -alUpdate commit = default

ops.integrator('LoadControl', 1.0)
ops.analyze(1)                                   # the real step (any algorithm/integrator)

ops.ladrunoBeginAugment()                        # recorders + commitTag frozen
ops.integrator('LoadControl', 0.0)               # HOLD everything, whatever drove the step
for _ in range(10):
    ops.analyze(1)                               # inner solve at FIXED lambda ...
    if ops.eleResponse(1, 'constraintViolation')[0] < tol:
        break                                    # ... then commitState does lambda += D g
ops.ladrunoEndAugment()
```

Use `LoadControl 0.0` for the held passes **whatever drove the real step** — a zero-increment
`DisplacementControl` is degenerate. Measured on the §9 gate at a moderate `K_t = 1e6`
(`≈1.3e2 · k_host`, where the penalty alone leaves `1.63e-2`), **20/20 cells converge**:

| driving integrator | algorithms | passes | final `max|g|/push` |
|---|---|---|---|
| `LoadControl` | Newton, ModifiedNewton, KrylovNewton, BFGS, Broyden | 5 | `6.14e-11` |
| `DisplacementControl` | Newton, ModifiedNewton, KrylovNewton, BFGS, Broyden | 4 | `9.53e-10` |

(both convergence tests, `NormUnbalance 1e-6` and `NormDispIncr 1e-12`, give the same cells).
This is the **supported** within-step route: it costs a handful of extra linear solves, it is
algorithm- and integrator-agnostic, and each inner solve is an ordinary penalty problem.

### 4.4 `-alUpdate iter` — opt-in, expert, and narrowly valid

`iter` advances the recursion inside `update()`, once per equilibrium iteration:

```
λ_{k+1} = λ_k + D g(u_k)                                   # in update()
r_k     = f − S u_k − Bᵀ( λ_{k+1} + D g(u_k) )
T       = S + BᵀDB                                          # tangent unchanged
```

The fixed point is still exact (`Δu = 0 ⇒ Δλ = 0 ⇒ D g = 0`) and under **full Newton +
LoadControl** it reaches it in a single step. But `update()` advances `λ` **before** the force
is formed, so the tie force carries `λ_k + 2·D·g(u_k)` against a tangent that linearises a
single `D·g`, and `λ_k` is **path-dependent**. The residual is therefore **not a function of
`u`**, and every secant / accelerated / re-solving method is fed `(δu, δr)` pairs that describe
no Jacobian. Full Newton survives only because its contraction on this gate happens to be
~0.008.

Measured on the §9 gate (linear elastic, `K_t = 1e6`) — each cell is the number of **failed
steps**:

| algorithm | `LoadControl` (1 step) | `DisplacementControl` (5 steps) |
|---|---|---|
| **Newton** | **0 — `err = 1.3e-12`** | 5 (all) |
| ModifiedNewton | refused | 5 (all) |
| KrylovNewton | refused | 5 (all) |
| BFGS | refused | 5 (all) |
| Broyden | refused | 5 (all) |

Before the guard, those cells did not refuse — they *failed*: `DisplacementControl` stagnated
geometrically (residual `78.85 → 70.62 → 62.17`, ratio 0.978, no `maxIter` rescues it, because
it re-solves `dLambda` every iterate against a residual that moves independently of `u`);
KrylovNewton diverged (`87.99 → 864.8 → 7267.9`) or, at a loose `1e-6`, "converged" to a gap
`1e3×` worse than Newton's; Broyden reached `2.5e275`. ModifiedNewton survives *this linear*
model but fails 10/10 on a `LadrunoBrick` bbar + `LadrunoJ2` host.

So `iter` is **refused at the first `update()`** unless the active algorithm is full Newton
**and** the active static integrator is `LoadControl` (read via `OPS_GetAlgorithm` /
`OPS_GetStaticIntegrator` / `OPS_GetTransientIntegrator`; silent while no analysis exists, so
the `Domain::addElement`-time `update()` does not trip it). The refusal is loud — the analysis
fails and the message names the offending class tags and points at §4.3. The parser echoes the
same window when `-alUpdate iter` is parsed.

> [!warning] `iter` buys exactness at the price of iterations, and it tracks the *displacement*
> tolerance — not the gap
> "`g = 0` exactly" means "to the algorithm's convergence tolerance". Measured (Newton,
> LoadControl, `K_t = 1e6`; penalty is `1.63e-2` at every row and converges in **2** iterations):
>
> | `NormDispIncr` tol | `iter` gap/push | `iter` iterations |
> |---|---|---|
> | `1e-14` | `2.98e-14` | 9 |
> | `1e-10` | `5.87e-11` | 7 |
> | `1e-8` | `2.68e-09` | 6 |
> | `1e-6` | `1.25e-07` | 5 |
>
> and it needs the iteration budget: with `maxIter = 3` the `iter` leg **fails** while the
> penalty leg converges in 2. Budget 5–9 iterations where penalty needs 1–2 (on a nonlinear
> `LadrunoJ2` host the same ratio showed up as 553 vs 110 total iterations under
> ModifiedNewton). §4.3's sweep costs extra *solves* instead, but keeps every inner solve a
> well-posed penalty problem.

### 4.5 AL is implicit-only — and under an explicit integrator it is refused *by consequence*

Under an explicit integrator there is no equilibrium iteration for the recursion to converge
against. Worse, the combination is not merely useless but **unusable**: `-bipenalty` is
refused together with `-enforce al` (the flag is dropped at parse with a warning), so a
massless tied DOF has **no mass source at all** and the explicit step is singular. Measured on
the §9 gate under `CentralDifferenceLadruno` (`analyze` return code):

| | result |
|---|---|
| `-enforce penalty -bipenalty -dtcr 1e-6` | **0** (runs) |
| `-enforce penalty`, no mass source | `−2` |
| `-enforce al -bipenalty -dtcr 1e-6` | `−2` — `-bipenalty` was dropped at parse |

`CentralDifference` and `ExplicitBathe` fail at step 0 the same way. **A parse-time refusal is
impossible** — the integrator does not exist when the element is declared — so the element
emits a one-time warning at the first `update()` when it sees a transient integrator active
with `-enforce al`. For explicit runs use `-enforce penalty` with `-bipenalty` (§5).

### 4.6 Bookkeeping (both cadences)

- `λ` is snapshotted at `commitState` (`lambdaCommitted`) and **rolled back by
  `revertToLastCommit`** — a failed/retried step must not inherit the multipliers of a
  discarded trial state. (Before WP-101 `revertToLastCommit` was a bare `return 0`, correct
  only while `λ` never moved inside a step.)
- `Domain::revertToLastCommit()`, `Domain::revertToStart()` **and `Domain::recv()`** all call
  `update()` right after they change the state, so a one-shot latch keeps `λ` from advancing on
  the just-reverted / just-restored state (see [[LEDGER_quirks]]; the un-latched database
  round-trip moved `λ` by `6.6e-9`).
- `sendSelf`/`recvSelf` carry `alUpdate` and `lambdaCommitted` (header version 2 — a newer
  payload is now **refused** rather than mis-read), payload `3·nGap`.

---

## 5. Explicit stability — bipenalty (`-bipenalty`, default OFF)

A penalty tie to a **massless** tied DOF has an unbounded frequency → zero stable step in
explicit central difference. RBE2's twist vs RBE3: **the master R is frequently a real,
massed node** (platen / footing / equipment), so `-bipenalty` **defaults OFF** and only fires
where you ask for it.

When on, the element scans **every tied DOF of R *and* every slave** and lumps a penalty mass
only on those that are **actually massless** — a **massless slave** is the RBE2-specific
hazard that RBE3's R-centric lumping would miss (validation test 9). Per-DOF stiffness for
sizing is the **Gershgorin row-sum** of the assembled penalty tangent (`≥ λ_max` ⇒
conservative; no eigensolve — a deliberate simplification vs the ADR's `jacobi3`):

```
-dtcr $dt :  m_p(DOF) = k_dof · (dt/2)²   on each massless tied DOF
             ⇒ self-reported critical step = min over lumped DOFs of 2√(m_p/k_dof) = dt
-wcap $β  :  m_p(DOF) = k_dof / (β·ω_host)²                 (needs -host for ω_host)
```

A `k_dof ≤ 0` guard skips untied/zero-stiffness DOFs, which is what prevents the
all-coincident `2√(I_p/0) = +Inf` trap. Penalty mass is lumped **diagonally** (`getMass`
stays diagonal → `DiagonalSOE`-safe). If **every** tied DOF already carries mass, the element
lumps **nothing** (R's own mass is not double-counted) and the self-report returns `0` —
"no opinion" (validation test 10). The bound is exposed via
`Element::getExplicitCriticalTimeStep` (folded into `ops.criticalTimeStep` / `-cflAbort`), and
Rayleigh factors are refused so a spurious `βK` can't shrink the step.

> [!note] `-dtcr` is a user-asserted target
> For a host-less node set there is no single host to derive an exact reduced-mass bound from;
> `-dtcr` sets the step you want the coupling to permit. Prefer `-dtcr` over `-wcap` for RBE2.

---

## 6. Deferred / not-yet-built

| Feature | Status |
|---|---|
| Finite-rotation geometric stiffness (large relative rotation) | deferred (§1.3) — shared with RBE3 |
| LS-DYNA CNRB mass-condensation integrator (true rigid-body inertia on R) | deferred (ADR 24 D3b) |
| General N-node linear-equation primitive `Σ cₖuₖ = 0` | separate element (ADR 24 D4) |
| RBE3 / distributing coupling (the flexible sibling) | shipped — [[LadrunoDistributingCoupling_guide]] |

---

### 6.x Corotational transport (scoping only — TIMs request 2026-09-07, F6)

**Status: scoped, not built; TIMs' G4 measures whether it matters first.**

What exists. `buildB()` (`LadrunoKinematicCoupling.cpp:335-350`) fills the gap operator
once, in `setDomain`, from the **reference** lever arms `d_i = X_i − X_R`
(`resolveGeometry`, `:200-275`): the translation row is
`u_i − u_R − [d_i]_× θ_R`, i.e. the small-rotation form `u_i = u_R + θ_R × d_i`.
`B` is constant, so `K = k·BᵀB` and `r = k·Bᵀg` are geometrically linear. The
reference-point transfer `Q_R2(z) = Q_R2(z_M) + Q_1·(z − z_M)` that TIMs post-process is
exact under that linearisation only; at a rocking angle `θ` the neglected term in the
slave position is `O(θ²)·|d_i|` (for `θ = 2°`, `6e-4·|d_i|`), and the moment carried through
the platen is off by the same order.

What would change (three parts, ADR 29 §2.4 conventions kept).

1. *Kinematics.* Replace `[d_i]_× θ_R` by the finite-rotation map: with
   `R(θ_R)` (Rodrigues on the reference node's current rotation vector, or an
   incremental rotation composed at each `update()`), the gap becomes
   `g_i = (x_R + R d_i) − x_i` with `x = X + u`. `B` is then rebuilt every
   `update()` from the **current** `R d_i` (the corotated arm), which is exactly the
   `-geom corot` idiom the fork uses on `LadrunoUP`/`LadrunoBrick` (ADR 78/79): the
   element does not become finite-strain, only the rigid arm rotates with the
   master.
2. *Tangent.* `K = k·BᵀB` is no longer the whole consistent tangent: the
   derivative of `R d_i` with respect to `θ_R` adds a geometric block
   `k·Σ_i gᵢ·∂²(R d_i)/∂θ_R²` (the "rotation-of-the-arm" term, symmetric for a
   penalty energy, of the same shape as the contact `∂n/∂u` block of ADR-39 B3).
   Dropping it keeps Newton convergent but linear-rate at large `θ`; keeping it
   restores quadratic convergence. Under `-enforce al` the augmented multiplier
   rides on the same `B`, unchanged.
3. *Rotation-row bookkeeping.* `θ_i − θ_R` for slaves that carry rotations stays
   additive only for small increments; a consistent version composes rotations
   (quaternion or rotation-vector update). For a translation-only slave skin (the
   TIMs footing) this part is not exercised.

Cost. One `R(θ_R)` per master per iteration and an `nGap × nDOF` rebuild of `B`
each `update()` (it is already allocated); the geometric block is a rank-`nGap`
correction assembled in the same loop. No new class tag, no new command: it
would be `-geom corot` on the existing element, default `linear` so every
existing deck is byte-identical. Explicit lane: `-bipenalty` mass is unaffected
(it is diagonal in the node DOFs), but the critical time step would pick up the
geometric block's contribution to the stiffness bound.

Trigger. Build it only if TIMs' G4 (rocking-angle sensitivity of the reference
transfer) exceeds their tolerance; the estimate above says a few degrees is
sub-percent. Until then this section is the scope, and the deferral is
recorded, not silent.

## 7. Diagnostics & responses

`eleResponse $tag <name>` / `recorder Element -ele $tag -<name>`:

| Response | Aliases | Size | Meaning |
|---|---|---|---|
| `force` | `couplingForce` | nGap | the coupling traction `t = D g` (incl. `λ` under AL) |
| `gap` | — | nGap | the constraint gap `g` (relative to `g0` if captured) |
| `kt` | `k` | 1 | resolved translational penalty `K_t` |
| `kr` | — | 1 | resolved (derived or user) rotational penalty `K_r` |
| `lambda` | `augLambda` | nGap | AL multiplier vector |
| `dtcr` | `dtCritical` | 1 | self-reported explicit critical step (`0` = no opinion; `−1` = bipenalty off) |
| `tiedDOFs` | `nGap` | 1 | total tied DOF count (`0` ⇒ the element went inert — a refusal fired) |

`tiedDOFs` is the quick health check: it should equal `Σ` over slaves of their tied
components (e.g. `6N` for a full rigid 3D tie of `N` 6-DOF slaves, `3N` for translation-only).

---

## 8. Implementation map

| Concern | Where |
|---|---|
| Element | `SRC/element/ladrunoKinematicCoupling/LadrunoKinematicCoupling.{h,cpp}` |
| Parser | `…/OPS_LadrunoKinematicCoupling.cpp` |
| Geometry + ragged layout | `resolveGeometry()` (`d_i = x_i − x_R`, `gapNode`/`gapDof`/`gapIsRot`, floored `ℓ²`, self-tie/duplicate refusals) — resolved **once** at `setDomain` from `getCrds()` |
| Constant gap operator | `buildB()` → `Matrix* B` (`nGap × nDOF`); `g = B·u − g0`; transport block `+[d_i]_×` (**sign-flipped** vs RBE3) |
| Residual / tangent | `getResistingForce` (`Bᵀt`), `getTangentStiff` (`BᵀDB`); `getInitialStiff ≡` tangent; per-row penalty via `rowPenalty(row) = gapIsRot(row) ? Kr : Kt` |
| Auto / derived penalties | `resolveAutoKt()` (`-k auto` off `-host`, derive `K_r = K_t·ℓ²`, **and the `K_t > 1e6·k_host` conditioning warning** — §3.1) |
| AL multipliers | `commitState()` (the DEFAULT `commit` cadence + the `lambdaCommitted` snapshot; this is also the outer update of the §4.3 augmentation sweep), `update()` (the opt-in `iter` cadence **and** its `refuseIterCadence()` guard + the AL-under-transient note), `revertToLastCommit()` / `recvSelf()` (roll back / restore + arm the one-shot latch) — §4.2–§4.6 |
| Damping bypass | `getDamp` / `getRayleighDampingForces` return element-owned **zeroed** `C0`/`dampF` (a no-op `setRayleighDampingFactors` alone would crash transient — see [[LEDGER_quirks]]) |
| Explicit | `getMass` (`M0` diagonal, massless-scan lumps), `resolveBipenalty()` (Gershgorin row-sum over R **and** slaves), `getResistingForceIncInertia`, `getExplicitCriticalTimeStep` |
| Serialization | `sendSelf`/`recvSelf` carry `dofSel`/weights flags/`λ`/`g0` + options; geometry, ragged layout & `B` **recomputed** on recv from coords |
| Plumbing | `classTags.h` (`ELE_TAG_LadrunoKinematicCoupling = 33012`), `FEM_ObjectBrokerAllClasses.cpp`, `OpenSeesElementCommands.cpp`, CMake — see [[LEDGER_vanilla_files]] |

### 8.1 Build-bug lessons folded in (don't relearn)
- **`getDamp` override is mandatory**, not optional: the no-op `setRayleighDampingFactors`
  means the base `Element` never allocates its lazy damping slot, so the base `getDamp` would
  deref `theMatrices[−1]` the first time a **transient** integrator forms the C-tangent → hard
  crash in the element's *primary* (dynamic) use. The 16 quasi-static-leaning tests can't catch
  it without the transient-Newmark smoke (test 11). Any element with a no-op
  `setRayleighDampingFactors` **must** override `getDamp()` + `getRayleighDampingForces()`.
- **Refuse inertly, after `allocate`** — self-tie / duplicate-slave set `valid=false` but still
  build zeroed matrices; an early-return-before-`allocate` would leave a null `K` → crash.
- **`getMass` double-zero guard** — `resolveBipenalty` is `bpResolved`-guarded so a second call
  doesn't wipe `M0`.
- **`lambdaAL` / `g0` first-call detection before resize** — resizing destroys the size-0
  fresh-vs-recv signal; check the size *before* resizing.
- **Test gotcha:** `ops.mass(node, …)` sizes the mass vector by the **model** `ndf`, not the
  node's — to mass a 6-DOF reference node, use `model -ndf 6` (massing a `-ndf 6` node in a
  `-ndf 3` model → "Node::setMass incompatible matrices").

---

## 9. Validation

Zone-A battery `tests/test_ladrunoKinematicCoupling_element.py` — **52/52** (16 at v1, +3 at
the 2026-09-07 u-p refusal, +7 at WP-101, +26 at the WP-101 review round 1: the 12-cell
algorithm × integrator × test sweep, the 10-cell augmentation sweep, and 4 guard cases),
plus full
Zone-A 633-pass no-regression and a 4-lens adversarial **code** review (2 CRITICAL + 6 MAJOR
folded in; ADR 29). The kinematic tests exploit a clean fact: with the slaves otherwise free,
their only stiffness is the penalty tie, so equilibrium drives each slave **exactly** onto R's
rigid prediction (gap → 0, independent of `K`).

- **Kinematics:** rigid translation (slaves follow, zero gap); **rotation + transport with the
  asserted SIGN** on an offset reference (the transOp sign-flip guard); full-rigid drives slave
  rotation (`θ_i = θ_R`, `tiedDOFs = 24`); translation-only `-dof 1 2 3` leaves slave rotation
  free (`tiedDOFs = 12`); single rigid link `N = 1` (generalized `rigidLink` with a moment arm).
- **Statics:** force balance (`Σ reaction = −F`); **moment into a 3-DOF face** as a
  self-equilibrated couple (`Σ f = 0`, `Σ d_i × f_i = −M`).
- **Robustness:** AL solves at finite `K`; derived `K_r = K_t·ℓ²` (`= 2·K_t` on the flat face);
  all-coincident → finite `dtcr` + floored `K_r > 0`; self-tie refused inertly.
- **Explicit / dynamics:** bipenalty massless-reference self-report (`dtcr == dt`);
  **massless-slave scanned** (the RBE2-specific hazard); **massed-R not double-counted**
  (`dtcr == 0`); **transient-Newmark smoke** (the `getDamp` regression guard).
- **Serialization:** FE_Datastore `sendSelf`/`recvSelf` round-trip (probes geometry-derived
  `kr` + `tiedDOFs` → fails if recv didn't reconstruct the element from coords).
- **Algorithm × integrator × test sweep (WP-101 r1)** — `Newton / ModifiedNewton /
  KrylovNewton` × `LoadControl / DisplacementControl` × `NormUnbalance / NormDispIncr`, 12
  cells, asserting per cell that **penalty** converges, that **`commit`** converges and equals
  penalty under `LoadControl`, and that **`iter`** either meets `1e-8` (full Newton +
  LoadControl only) or is **refused** — never a silent wrong answer. This sweep is nine lines
  and it is what the first cut of WP-101 was missing.
- **Augmentation sweep (WP-101 r1)** — §4.3 driven for 5 algorithms × 2 integrators, asserting
  `max|g|/push ≤ 1e-8` within one step and a bounded pass count.
- **Rigidity gate (WP-101)** — a 2×2×2 `stdBrick` block (`E = 45 000 kPa`, `ν = 0.3`,
  `B = 1.5 m`, base fixed) whose 9 top-face nodes are a footing skin driven by a fully
  prescribed 6-DOF master; the metric is `max|g| / |push|`, cross-checked against a **direct
  push of the same footprint** (leg B). It pins: the penalty law `err = c/K_t` over three
  decades (`c = 1.66e4`, constant to 5 %); `-enforce al` closing to `2.98e-14` in **one** step
  at `K_t = 1e6` where penalty leaves `1.63e-2`, with the base reaction matching leg B to
  `1e-9`; `-alUpdate commit` reproducing the penalty gap exactly (the legacy pin); `-k auto
  -host` measured **not** rigid (`K_t ≈ 7.9e6`, `err = 2.1e-3`); `λ` restored after a
  deliberately failed step; and both new warnings (over-stiff `K_t` vs `-host`, `-alUpdate`
  without `-enforce al`).

---

## 10. Use cases & recipes

### 10.1 Loading platen / rigid bearing on a deformable face (the driver)
A rigid platen pressing on a solid specimen face: drive the face nodes rigidly from one
control node, prescribe its motion or load it, and the patch translates+rotates as one rigid
body (contrast RBE3, which would let the face deform under the load):

```python
ops.model('basic', '-ndm', 3, '-ndf', 3)
# ... solid face nodes 101..104 (ndf 3) ...
ops.node(1, xc, yc, zc, '-ndf', 6)                         # platen control node (master, 6-DOF)
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-dof', 1, 2, 3,                               # tie translations (3-DOF slaves)
            '-k', 'auto', '-host', 5001)                   # scale penalty off solid element 5001
ops.sp(1, 3, -0.01)                                        # press the platen down 10 mm
```

### 10.2 Rigid offset / rigid link (generalized rigidLink)
Tie a node to a master through a rigid arm — a beam end to an offset working point, a
sensor node rigidly carried by a member, a rigid connection block. `N = 1` reduces exactly to
`rigidLink` but with full moment-arm transport and an offset reference:

```python
ops.element('LadrunoKinematicCoupling', 7, masterNode, 1, slaveNode, '-k', 1.0e10)
```

### 10.3 Rigid body / rigid diaphragm over an arbitrary node set
Make a group of nodes move as one rigid body driven by a single reference (a rigid footing, a
rigid cap, an equipment skid, an arbitrary-shape diaphragm that `rigidDiaphragm` can't
express). Use `-dof` to tie only the in-plane DOFs for a diaphragm, or all 6 for a full rigid
body.

### 10.4 Beam/shell-to-solid rigid moment transfer (ndf-mismatch, rigid variant)
The rigid counterpart of RBE3's signature transfer: a 6-DOF beam/shell node framing into a
3-DOF solid face, delivering the member's moment into the continuum as a force couple **while
holding the patch rigid**. Choose RBE2 here when the joint block is genuinely stiff (a thick
gusset, an embedded plate); choose [[LadrunoDistributingCoupling_guide|RBE3]] when the face
must stay flexible.

```python
ops.node(1, xb, yb, zb, '-ndf', 6)                         # beam node (6-DOF master)
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104, '-dof', 1, 2, 3)
ops.load(1, 0, 0, -P, 0, Mx, 0)                            # force + moment at the master
```

### 10.5 Explicit dynamics (impact / SSI / blast)
R is **often massed** here, so bipenalty stays off unless a tied DOF is massless. If R or a
slave is massless, add `-bipenalty -dtcr <dt>` at/below your explicit step and query the bound:

```python
ops.mass(1, m, m, m, J, J, J)                              # massed master -> no penalty mass lumped
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-k', 1.0e8, '-bipenalty', '-dtcr', 2.0e-6)    # fires only on massless tied DOFs
dt_cr = ops.eleResponse(1, 'dtcr')[0]                      # honored by ops.criticalTimeStep / -cflAbort
```

### 10.6 Staged construction (stress-free birth)
`g0` capture is **on by default**: a coupling added to an already-deformed model is born
force-free (the gap is measured relative to the activation state, including any pre-existing
offset between R's predicted position and the slaves). Add `-absolute` only for the legacy
absolute snap-to-master tie.

### 10.7 apeGmsh integration — generator contract (recommended)
For a rigid platen/connection-on-solid interface, apeGmsh should emit one
`LadrunoKinematicCoupling` per master: resolve the host face the master drives, pass those
nodes as the slave set, choose `-dof` by the slave ndf (translations for a 3-DOF face), name a
representative face element via `-host` for `-k auto`, and add `-bipenalty -dtcr` only when a
tied node is massless. A `g.couple(ref, set, mode='kinematic', dof=...)` wrapper is the natural
API surface (cf. ADR 24 §6 and the RBE3 `mode='distributing'` sibling).

### 10.8 Pitfalls
- **Expecting RBE2 to leave the patch flexible** → it doesn't; it makes the set **rigid**.
  Use [[LadrunoDistributingCoupling_guide|RBE3]] for a load-introduction-without-stiffening.
- **`tiedDOFs == 0`** → the element went inert: a **self-tie** (slave tag equals the reference)
  or a **duplicate slave**. Check the slave list.
- **Reference node built with too few DOFs** (e.g. `ndf 3` in 3D) → refused at `setDomain`.
  Build it `ndf 6` (3D) / `ndf 3` (2D).
- **`-k auto` without `-host`** → refused (no single host for a node set); use a numeric `-k`.
- **Expecting `-k auto` to make the patch rigid** → it won't; it is a *conditioning* control
  and pins `K_t` to the host's own order (§3, measured `1.4e-4` of the push left on the TIMs
  gate). Raise `-k` into the §3.1 band, or add `-enforce al`.
- **Reaching for `-k 1e12` on a production model** → the rigidity error stops improving
  (round-off floor) while the factorisation degrades; a 101 583-DOF strip went near-singular
  (Pardiso perturbed pivots → failure; SuperLU failed its first factorisation) at `1e12` and
  ran clean at `5e9`. Stay in `1e2…1e4 × k_host` and use `-enforce al` for tightness (§3.1).
- **Expecting `-enforce al` to converge the tie in ONE step on its own** → it does not: the
  default `-alUpdate commit` updates `λ` once per *committed* step, so a single push returns
  the plain penalty gap. Wrap the step in the §4.3 held-load augmentation sweep — that is the
  supported within-step route and it works with every algorithm and integrator.
- **Reaching for `-alUpdate iter` to avoid the sweep** → it is refused outside full Newton +
  LoadControl, and for good reason: it makes the residual path-dependent, so
  `DisplacementControl` fails 5/5 steps with every algorithm and KrylovNewton / BFGS / Broyden
  diverge (§4.4).
- **`-enforce al` in an explicit run** → refused by consequence: `-bipenalty` is dropped when
  `-enforce al` is given, leaving a massless tied DOF with no mass source (§4.5).
- **Single penalty across translation and rotation** → don't override `-kr` with `K_t`; the
  derived `K_t·ℓ²` is there for a reason (§3).
- **Bipenalty default surprise** → unlike RBE3, RBE2 bipenalty is **off** by default (R is
  usually massed). Turn it on explicitly if any tied node is massless in an explicit run.

---

## 11. References & related
- **[[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]]** — the design record (formulation, the
  transport sign flip, the ragged layout, explicit safety, the 4-lens review must-fixes).
- **[[24_ladruno_coupling_constraints_adr|ADR 24]]** — the coupling-constraint family (RBE2 /
  RBE3 / linear equation) and the Abaqus/LS-DYNA/Nastran alignment.
- **[[LadrunoDistributingCoupling_guide]]** — the **flexible sibling (RBE3)**: R is the
  dependent weighted-average of the set, adds no stiffness. Reach for it when RBE2 would
  over-stiffen.
- **[[LadrunoEmbeddedNode_guide]]** — the node-into-host **embedment** tie (shares the
  penalty/AL/bipenalty kernel).
- **[[ndf_and_mixed_models_guide]]** — mixed-ndf models; §4 the rotational-compatibility gap;
  §7 the explicit zero-mass-rotational trap.
- **Nastran** RBE2 · **Abaqus** `*COUPLING, kinematic` · **LS-DYNA**
  `*CONSTRAINED_NODAL_RIGID_BODY`.
