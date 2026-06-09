---
title: "ADR 29 — LadrunoKinematicCoupling (RBE2 / kinematic coupling): implementation spec"
project: Ladruno
type: ADR / implementation spec
status: BUILT — formulation locked, 4-lens adversarially reviewed, element shipped; Zone-A 16/16 + full Zone-A 633 no-regression
related:
  - "[[24_ladruno_coupling_constraints_adr]]"        # parent: the coupling-constraint family (D3)
  - "[[28_ladruno_distributing_coupling_rbe3_adr]]"  # sibling: RBE3 (roles flip — R is dependent there)
  - "[[23_ladruno_embedded_node_adr]]"               # the reused DOF-scatter / bipenalty scaffold
  - "[[20_ladruno_embedded_reinforcement_adr]]"      # the shared penalty/AL/bipenalty kernel lineage
  - "[[LEDGER_implementations]]"
tags: [adr, constraints, coupling, rbe2, kinematic-coupling, rigid-body, explicit-dynamics, abaqus, ls-dyna, nastran]
updated: 2026-06-09
---

# ADR 29 — `LadrunoKinematicCoupling` (RBE2 / kinematic coupling)

**Status:** planned. Formulation locked and **adversarially reviewed** (4-lens sweep:
mechanics / OpenSees-integration / explicit-dynamics / degeneracy). This is the
implementation spec; code has not landed. It is the **D3** build under
[[24_ladruno_coupling_constraints_adr|ADR 24]] (the coupling-constraint family) and the
rigid sibling of the shipped RBE3 [[28_ladruno_distributing_coupling_rbe3_adr|`LadrunoDistributingCoupling`]]
(ELE 33011).

> [!info] What RBE2 is, in one line
> A **reference node R** that **rigidly drives** a set of **slave nodes** — each slave
> follows R's rigid-body motion `u_i = u_R + θ_R × (x_i − x_R)` (and `θ_i = θ_R` where the
> slave carries rotations). A load/displacement at R is transmitted to the set *with* rigid
> stiffness. = Abaqus `*COUPLING, kinematic` = LS-DYNA `*CONSTRAINED_NODAL_RIGID_BODY` =
> Nastran **RBE2**. Contrast RBE3 (distributing coupling): R is the *weighted average* of
> the set, adding no stiffness. **In RBE2 the single node R is the master/independent; in
> RBE3 it is the dependent/slave.** The roles flip — the classic wiring mistake.

---

## 1. Driver & goal

OpenSees has no general 3D rigid kinematic coupling: `rigidLink` ties **one** slave to one
master (beam/bar offset), and `rigidDiaphragm` is **planar / in-plane only**. The routinely-
wanted capability — a control node rigidly driving an *arbitrary 3D slave set* with
*selectable* dependent DOFs — has no equivalent. RBE2 is the **"rigid region / rigid driver"**
tool: rigid end zones, panel zones, loading platens, rigid footings/pile caps, fastener
heads, eccentric rigid offsets, lumped rigid masses with rotary inertia (see §5).

**v1 scope (locked):** a reference node R rigidly driving N slaves, 2D + 3D, with a
**`-dof` selection** of which slave DOFs are dependent — one element subsuming the three
constraint *modes* (translation-only, translation+transport, full rigid body). Penalty +
augmented-Lagrangian, explicit-safe via bipenalty. Kinematically-complete transport term
`θ_R × d_i` for the offset between R and each slave. classTag **33012**.

**Non-goals (v1):** finite-rotation configuration update (geometric stiffness — §4.2);
the LS-DYNA-style rigid-body **mass-condensation integrator** (ADR 24 D3b — penalty is v1);
weighted averaging (that is RBE3, not RBE2 — RBE2 has **no weights**); the general N-node
linear-equation primitive (ADR 24 D4).

---

## 2. Formulation (locked, adversarially verified)

### 2.1 Notation
Reference node **R**: translations `u_R`, rotations `θ_R`, position `x_R`. Slave nodes
`S_i`, `i=1..N`: translations `u_i`, rotations `θ_i` (if 6-DOF), position `x_i`. The
**lever** (reference configuration) `d_i = x_i − x_R`. `ndm ∈ {2,3}`; `nrot = 3` (3D) or
`1` (2D drilling). `[d_i]_×` = the skew/cross matrix (`[d_i]_× v = d_i × v`).
**2D convention (read carefully — see review M1/M4):** the *prediction* `θ × d` for scalar
`θ` is the column `θ·[ −d_{i,y}, d_{i,x} ]ᵀ`; but the `g_i^t` **B-block** is
`∂g_i^t/∂θ_R = −∂(θ×d_i)/∂θ_R`, i.e. the **negative** of that — the 2D transport B-column is
`[ +d_{i,y}, −d_{i,x} ]ᵀ`. (In 3D the B-block is `+[d_i]_×`; see §2.4.)

**There is no weighted centroid and no position-inertia `I_c`.** Unlike RBE3, RBE2 fits
nothing — it is a *direct, per-slave* rigid tie. The Jacobi eigensolver, `I_c`, its
pseudo-inverse and the degeneracy machinery of ADR 28 §3 are **not needed**.

### 2.2 Kinematic constraint (the relation RBE2 enforces)
Each slave follows R's rigid-body field (small rotation):
```
u_i = u_R + θ_R × d_i          (translation)
θ_i = θ_R                       (rotation, where tied)
```
**Constraint gaps** (per slave):
```
g_i^t = u_i − ( u_R + θ_R × d_i )        (size ndm)   — for tied translation DOFs
g_i^r = θ_i − θ_R                         (size nrot)  — for tied rotation DOFs (6-DOF slave)
```
Using `θ_R × d_i = [θ_R]_× d_i = −[d_i]_× θ_R`, the translation gap is linear in the DOFs:
`g_i^t = u_i − u_R + [d_i]_× θ_R`. The total gap vector stacks the *tied* rows of each
slave; its length is `Σ_i (n_tied_trans,i + n_tied_rot,i)`.

### 2.3 The `-dof` selection — three modes, one element
`-dof {c1 c2 …}` lists the dependent **component numbers** applied to **every** slave
(Nastran-RBE2 / Abaqus-kinematic semantics): `1..ndm` = translations, `ndm+1..ndm+nrot`
= rotations. Default (no `-dof`) = **all DOFs each slave possesses**. The three "constraint
types" are then just selections:

| Mode | `-dof` (3D) | Gaps active | Behavior | Typical use |
|---|---|---|---|---|
| **Translation + transport** | `1 2 3` | `g_i^t` only | slave translations rigid to R incl. moment-arm `θ_R×d_i`; slave rotations (if any) float | 6-DOF node → 3-DOF solid face, rigid variant |
| **Full rigid body** | `1 2 3 4 5 6` (default for 6-DOF slaves) | `g_i^t` + `g_i^r` | slave moves as a rigid extension of R | rigid block, lumped rigid mass+inertia, panel zone |
| **Translation-only** | `1 2 3` on 6-DOF slaves | `g_i^t` only | identical operator to mode 1 — never ties slave rotation even when present | tie position, leave slave spin free |

**2D numbering:** `1,2` = translations, `3` = drilling (there is no third translation). A
`-dof 3` against a 2-DOF (no-drilling) node is skipped+warned, exactly like a 3D rotation
component on a 3-DOF slave.

A listed rotation component on a slave that lacks it (3-DOF slave under `-dof 4 5 6`; 2-DOF
node under `-dof 3`) is **skipped for that slave + warned** (not an error — mixed slave sets
are legal). Two refusal/validation rules (review D4, D5):
- **Parse-time** (`ndm` known): every `-dof` component must satisfy `1 ≤ c ≤ ndm+nrot`;
  reject out-of-range (`0`, `7`, negative, `>6`) and a syntactically empty `-dof` list there.
- **setDomain** (ndf known): if the *effective* selection (after ∩ each slave's ndf) is
  empty, **refuse inertly** (`valid=false`, zeroed B — never null matrices / no hard exit).

### 2.4 Penalty enforcement & the work-conjugate reaction
Penalty potential `Π_p = Σ_i ( ½ K_t |g_i^t|² + ½ K_r |g_i^r|² )`; tractions
`t_i^t = K_t g_i^t` (+`λ_i^t` for AL), `t_i^r = K_r g_i^r` (+`λ_i^r`). Resisting force
`P = Bᵀ t`, tangent `K = Σ_i Bᵢᵀ D_i Bᵢ`, `D_i = diag(K_t I, K_r I)` over the tied DOFs.

**The B-operator** (constant in the reference configuration):

| ∂/∂ | `g_i^t` (trans rows) | `g_i^r` (rot rows) |
|---|---|---|
| `u_i`  | `+I`                | — |
| `θ_i`  | —                   | `+I` |
| `u_R`  | `−I`                | — |
| `θ_R`  | `+[d_i]_×`  (transport) | `−I` |

The `θ_R` transport block reuses the **same `transOp` cross-product machinery** as RBE3's
`buildB` — but with the per-slave lever `d_i = x_i − x_R` and **`B += −transOp(·;d_i)`, a
sign FLIP versus RBE3** (review M1, MAJOR): RBE3 carried `+θ_R×(x_c−x_R)` and assembled
`B += +transOp`; RBE2's prediction is `−θ_R×d_i`, and since `transOp = ∂(θ×d)/∂θ = −[d_i]_×`,
the correct block is `−transOp = −(−[d_i]_×) = +[d_i]_×`. **Do not copy RBE3's `buildB`
`+= transOp` line verbatim — flip the sign**, else every offset-moment reaction reverses
(validation §8 case 2 must assert the reaction *sign*, not just magnitude). RBE3's `crossOp`
(the `I_c⁺ C_i` lever block) is **not used** — RBE2's rotation rows are the trivial identity
blocks `+I`/`−I`.

**Reactions delivered** (= `−` the relevant row of `Bᵀt`):
```
slave i translation:  +t_i^t            slave i rotation:  +t_i^r
reference  force:      −Σ_i t_i^t
reference  moment:     −Σ_i d_i × t_i^t  −  Σ_i t_i^r
```
The element is **self-equilibrated** for *any* K: `Σ force = 0` and `Σ moment_R = 0`
(the rigid-body B does no work on rigid-body motion). Load conservation is exact
regardless of penalty magnitude — the RBE2 analog of ADR 28 §2.3 Findings 5 & 7.

### 2.5 Reuse boundary — the transport block is in-element
Identical to ADR 28 §2.4: the `g_i^t` row couples `θ_R` to slave *translations* through the
lever `[d_i]_×` — a **cross-DOF-class block** the per-node `LadrunoEmbeddedKernel`
structurally cannot assemble (`block(p,q)=c_p c_q D`). RBE2's per-slave translation tie is
itself so simple (direct `+I`/`−I`) that the kernel buys nothing and there are **no
shape-function weights** to consume. Therefore **RBE2 does not depend on the shared
kernel** — it assembles its own `B`/`K` (copying the EmbeddedNode/RBE3 *scaffold*, not the
kernel). "~70% scaffold reuse" holds for DOF-scatter / bipenalty / AL / serialization;
`resolveGeometry`+`buildB` are new (and simpler than RBE3's).

---

## 3. Ill-posed cases (no `I_c`, so far fewer than RBE3)

Resolved **once** at setDomain (coords known; freeze for the analysis). RBE2 sheds the
`I_c`/fit degeneracies entirely (no eigensolve, no pseudo-inverse), **but introduces three
new classes RBE3 never had** (review: ragged gap layout, the all-coincident `K_r,eff=0`
corner, `-dof` range). The table grew accordingly:

| Case | Action |
|---|---|
| **R has no rotation DOFs** (`ndf_R < ndm+nrot`) | if `-dof` was defaulted/translation-only: drop transport + all `g_i^r`, degrade to `g_i^t = u_i − u_R`, **WARN** ("reference carries no rotations; multi-slave equalDOF, no moment transfer"). **If `-dof` explicitly named a rotation component: REFUSE** (review D3 — silently dropping the user's moment-transfer intent is a silent-wrong-answer). |
| **Slave lacks a `-dof`-listed component** (3-DOF slave under `-dof 4 5 6`; 2-DOF under `-dof 3`) | skip that slave's gap row + **WARN** (mixed sets legal). |
| **Empty effective `-dof`** (nothing tied after ∩ ndf) | **REFUSE inertly** at setDomain (`valid=false`, zeroed B — not a hard exit). Syntactic-empty `-dof` is caught earlier at parse. |
| **`-dof` component out of range** (`0`,`7`,neg,`>ndm+nrot`) | **REJECT at parse time** (review D4) — `ndm` is known in the parser; gives the error at the offending line, avoids an OOB tied-DOF index build. |
| **Duplicate slave tag**, or **slave tag == refNode** | **REFUSE** (review D6): self-tie is a zero gap row tying R to itself (always a user error); a duplicate doubles that node's penalty + over-counts the moment reaction (silent wrong answer). RBE3's `setDomain` lacks this check too — fix both. |
| **N = 1** | valid — a general single rigid link with moment arm (= `rigidLink` + transport). |
| **All slaves coincident with R** (`d_i = 0 ∀i`) | transport block = 0. **If a rotation is tied:** fine *iff* the default `K_r` is floored off the lever (see §6 — `ℓ²=0` must NOT zero `K_r`). **If translation-only:** R's rotation is genuinely unconstrained → drop it, **WARN**, and **do NOT lump `I_p` on it** (guard the `2√(I_p/K_r,eff)` divide against `K_r,eff=0` → +Inf). Review D1, **CRITICAL**. |
| **Single coincident slave, `d_i = 0`** | fine — transport zero, pure position tie (only the all-coincident + translation-only combination above is degenerate). |

There is **no `W≤0` / negative-weight / collinear / `I_c`-rank failure** — RBE2 has no
weights and no fit.

**Ragged per-slave gap layout (review D2, MAJOR).** Unlike RBE3 (fixed `(ndm+nrot)×nDOF`
B), RBE2's gap length `Σ_i(n_tied_trans,i + n_tied_rot,i)` **varies per slave** (a 6-DOF
slave under `-dof 1..6` adds 6 rows; a 3-DOF slave adds 3). setDomain MUST build an explicit
per-slave layout table — `gapRowStart[i]`, `nTiedTrans[i]`, `nTiedRot[i]` (the `-dof`∩ndf
result) — and `buildB` + the reaction assembly (`−Σ d_i×t_i^t − Σ t_i^r`) MUST index through
it. A flat uniform-stride assumption mis-indexes a mixed 3-DOF/6-DOF set → a self-equilibrated
but *wrong* moment reaction (silent). §8 case 4 must include a **mixed-ndf** slave set.

---

## 4. Consistent tangent & objectivity

### 4.1 Reference-configuration tangent (v1)
`B` is constant (reference geometry only), `t = D g` is linear in the DOFs, so
`K = Σ_i Bᵢᵀ D_i Bᵢ` is the **exact penalty Hessian** — symmetric PSD, constant.
`getInitialStiff = getTangentStiff`. Resolve `d_i` eagerly at setDomain.

### 4.2 Objectivity boundary (documented, not a bug)
Same boundary as ADR 28 §4.2. Under rigid translation of the whole assembly `g_i^t → 0`
exactly. Under a **finite** rigid rotation, both `g_i^t` (through the linearized
`θ_R × d_i`) and `g_i^r` vanish only to **first order** in the relative rotation, because
`θ_R` is a finite nodal-rotation DOF while the tie is linear in `u_i`. Acceptable for the
small-relative-rotation rigid-region driver. **Finite-rotation support is a future
geometric-stiffness term** (`∂Bᵀ/∂x · t` from co-rotating `d_i`), *not* a parser tweak —
the "exact & constant tangent" claim is explicitly scoped to the reference configuration.

---

## 5. Use cases — when to reach for `LadrunoKinematicCoupling`

The **"rigid region / rigid driver"** tool: a slave set moves rigidly with one control
node — reach for it whenever that rigid stiffness is *physically intended*.

- **Rigid end zones / panel zones** — a beam–column joint core or rigid offset driven
  rigidly by a control node; a 3D, arbitrary-set generalization of `rigidDiaphragm`.
- **Loading platen / displacement-controlled boundary** — drive a specimen face uniformly
  through one control node (actuator, rigid loading head, compression platen); prescribe
  `u`/`θ` at R, the whole face follows.
- **Rigid footing / foundation block / pile cap** — tie a node set into one rigid body that
  translates and rotates as a unit.
- **Fastener / bolt head / rivet** — rigidly tie a bolt-hole node ring to one fastener
  reference node (the connector's rigid shank).
- **Rigid offset / lever arm / moment arm** — connect a load or mass at an eccentric control
  point through the kinematically-exact transport `θ_R × d_i` (what `rigidLink` does for
  *one* slave, generalized to a set).
- **Lumped rigid mass with rotary inertia** — attach a rigid equipment block so its mass
  *and* rotary inertia drive a flexible support patch as one body (full-rigid `-dof 1..6`);
  contrast RBE3's *soft* mass distribution.
- **Frame/shell ↔ solid moment transfer, rigid variant** — a 6-DOF beam/shell node into a
  3-DOF solid face where the interface is meant to be *rigid* (embedded anchor plate);
  translation+transport carries the moment as a force couple and holds the patch rigid.
  *(When the face must stay flexible, use RBE3 — the explicit RBE2-vs-RBE3 choice.)*
- **Rigid tie of repeated nodes** — force a stiffener line / 3D diaphragm edge to share one
  master's motion exactly.

> [!tip] RBE2 vs RBE3 (the wiring choice)
> **Region must move as a rigid body, stiffness intended → RBE2** (this element).
> **Apply/transmit a load at a point while the set stays flexible → RBE3**
> ([[28_ladruno_distributing_coupling_rbe3_adr|`LadrunoDistributingCoupling`]]). In RBE2 R is
> the **master** (drives the set); in RBE3 R is the **dependent** (averaged from the set).

---

## 6. Explicit-dynamics safety (the fork's center of gravity)

**The RBE2 role-inversion breaks RBE3's bipenalty story — do NOT mirror RBE3 §5 naively
(review E1/E2, CRITICAL+MAJOR).** In RBE3 the reference R is *structurally* the massless
interpolation point and the slaves are the always-massed host nodes, so lumping `m_p`/`I_p`
on slot-0 R is always safe. **RBE2 inverts both premises:** R is the *master*, frequently a
**real massed node** (a loading platen, footing, or the lumped equipment block of §5's own
use cases — these attach a `mass` to R); and the **slaves are the dependents**, which can
themselves be **massless** (a geometric attachment node, a fastener ring, a mid-face node
with no rotary mass). Both inversions are correctness hazards:

- **`-bipenalty` defaults OFF for RBE2** (vs RBE3's default-on). Unconditionally lumping
  `m_p` on slot-0 R **double-counts mass** when R already carries nodal mass `M_R`
  (`getMass` contributions are *summed* into global M) — silently inflating R's inertia and
  the inertial residual `m_p·a_R`, corrupting the dynamics. This is a *correctness* bug, not
  only a CFL bug. At setDomain, query `theNodes[i]->getMass()` on the tied DOFs; **lump a
  penalty mass only where the node is actually massless**, and warn if `-bipenalty` is given
  while R carries mass.
- **Scan EVERY tied DOF of R *and every slave*, not just slot 0** (review E2). A massless
  slave DOF has tie stiffness (`K_t` translation / `K_r` rotation) but no inertia → its own
  CFL singularity that R-only bipenalty misses → silent instability. Generalize RBE3's
  slot-0 lump to a per-DOF massless scan over all `1+N` nodes.
- **Two stiffness paths into `θ_R`** (the RBE2-specific subtlety vs RBE3's single path):
  1. *transport* — each `g_i^t` depends on `θ_R` via `[d_i]_×`, contributing
     `K_t Σ_i (|d_i|²I − d_i⊗d_i)` to the 3×3 `K_{θθ}` block;
  2. *direct rotation tie* — each `g_i^r = θ_i − θ_R` contributes `+K_r·I` (the `−I` block),
     i.e. `N_rot·K_r·I`.
  So `K_{θθ} = K_t Σ_i(|d_i|²I − d_i⊗d_i) + N_rot·K_r·I`. **Use `K_r,eff = λ_max(K_{θθ})`,
  NOT a scalar `Σ|d_i⊥|²` (review E3, MAJOR):** the cross terms `d_i⊗d_i` split the
  eigenvalues, so for a planar slave cloud (a platen!) the out-of-plane axis carries the
  full `Σ|d_i|²` while in-plane axes carry less — a scalar average under-estimates `λ_max`
  and over-reports `dt_cr` → unstable. Form the 3×3 block at setDomain and run the
  already-present `jacobi3` (2D: scalar `K_θθ = K_t·Σ|d_i|² + N_rot·K_r`). **Without a stored
  `K_r,eff` the rotation mode is unbounded/unreported → silent instability** (ADR 28 §5 trap)
  — the #1 explicit must-fix.
- **Penalty-mass sizing must match the *actual* DOF stiffness (review E4, MAJOR).** The
  inherited RBE3 `resolveBipenalty` sizes `iPenalty` against `Kr` and `mPenalty` against
  `K_t`, making the `2√(m/k)` self-report a tautology that returns the user's `-dtcr` while
  the true stiffness is larger. For RBE2: size the **rotation** penalty mass and the
  `criticalTimeStep` denominator against **`K_r,eff`** (not `Kr`); size R's **translation**
  penalty against **`N·K_t`** (R's single translation DOF is pulled by N unweighted ties in
  parallel — `K_{u_R u_R}=N·K_t`, not `K_t`), else the translational `dt_cr` is over-reported
  by `√N`. Self-report `min( 2√(m_p/N·K_t), 2√(I_p/K_r,eff) )`, skipping any class whose
  stiffness floors to 0 (the all-coincident, no-rotation-tied case — §3).
- **Default `K_r` must be decoupled from the lever (review D1, CRITICAL).** Units: `K_t` is
  force/length, `K_r` is moment/rad = force·length, so `K_r = K_t·ℓ²`. But `ℓ² = Σ|d_i|²/N`
  **= 0 when all slaves are coincident with R**, which would zero the *direct* rotation-tie
  penalty even though that tie has nothing to do with any lever — the rigid body would
  silently fail to constrain slave spin. Floor it: `ℓ² = max(Σ|d_i|²/N, L_char²)` with
  `L_char` a model length scale (largest `|d_i|`, or the max slave-pair distance, or a tiny
  unit-length ε as last resort). User override `-kr` bypasses the default.
- **bipenalty default OFF; auto-enable only on a detected massless DOF under an explicit
  integrator.** Mass redistribution is **not applicable** (RBE2 adds stiffness, it does not
  average mass).
- **AL is implicit-only**: Uzawa updates `λ` in `commitState`; under explicit central
  difference it fires every sub-step with no equilibrium iteration → `λ` ratchets on the
  un-converged gap. **Warn if `-enforce al` meets an explicit integrator** (reuse the
  EmbeddedNode/RBE3 interlock).
- **Refuse Rayleigh factors** — a pure penalty coupling carries no physical damping; a `βK`
  would spuriously shrink `dt_cr`. **AND override `getDamp` / `getRayleighDampingForces`**
  to return element-owned zero matrices — the no-op `setRayleighDampingFactors` never
  allocates the base `Element`'s lazy damping slot (`index` stays −1), so the base path
  dereferences `theMatrices[-1]` the first time a **transient** integrator forms the
  C-tangent (Newmark nonzero c-factor) → **hard crash**. This exact bug was confirmed
  empirically on RBE3 (ADR 28 §8b, exit 5) and must be carried over from day one (D ≡ 0 is
  physically correct; mass/inertia come from `getMass` + bipenalty).
- **`-k auto` is undefined for a node set** (no host element). Default to numeric `-k`;
  `-k auto` only with an explicit representative `-host` element among the slaves' parents.

---

## 7. Implementation map

### 7.1 Files (new element dir, mirrors RBE3)
```
SRC/element/ladrunoKinematicCoupling/
  LadrunoKinematicCoupling.{h,cpp}     # element: DOF scatter, d_i resolve, per-slave B, K_r,eff
  OPS_LadrunoKinematicCoupling.cpp     # parser
  CMakeLists.txt
tests/test_ladrunoKinematicCoupling_element.py
```

### 7.2 Reuse (copy the scaffold from RBE3 / EmbeddedNode)
- **Copy** `SRC/element/ladrunoDistributingCoupling/LadrunoDistributingCoupling.{h,cpp}`:
  per-node DOF scatter (`nodeNdf`/`dofOffset`/`nDOF`/`allocate`), penalty/AL (`lambda`,
  Uzawa in `commitState`), bipenalty (`mPenalty`/`iPenalty`, `getMass` diagonal-on-slot-0,
  `setRayleighDampingFactors` refusal, **`getDamp`/`getRayleighDampingForces` overrides**,
  `getExplicitCriticalTimeStep`), g0 initial-gap capture (`-absolute` opt-out), the
  `transOp` machinery, `sendSelf`/`recvSelf` header-vector + recompute-geometry-on-recv,
  `setResponse`/`getResponse`, `Print`.
- **Replace** `resolveGeometry` + `buildB`: **drop** `jacobi3`, the weighted centroid,
  `I_c`/`Icplus`/`Pproj`/`nKept`, `rvec`-from-centroid, `weights`, `ell2`-from-centroid.
  **Add** per-slave `d_i = x_i − x_R`, the `-dof` tied-DOF bookkeeping (`tiedTrans`,
  `tiedRot` per slave from the component list ∩ each slave's ndf), the per-slave B (§2.4),
  and the two-path `K_r,eff` (§6).
- **Do NOT inherit these RBE3 behaviors verbatim** (the review's must-fixes — the copy will
  carry the bug otherwise): (1) the `buildB` transport line — **flip the sign** to `−transOp`
  (M1); (2) `resolveBipenalty` sizing — use `K_r,eff`/`N·K_t`, not `Kr`/`K_t` (E4); (3)
  bipenalty **default-on** — make it **opt-in + massless-scan over all nodes** (E1/E2);
  (4) `Kr = Kt·ell2` default — **floor `ell2`** (D1); (5) `setDomain` missing dup/self-tie
  checks (D6).

### 7.3 Plumbing (Definition-of-Done — mirror ADR 28 §6.3)
1. `SRC/classTags.h` — `#define ELE_TAG_LadrunoKinematicCoupling 33012` (next free ELE tag
   after `LadrunoDistributingCoupling=33011`; ELE band is per-registry — no collision with
   `ND_TAG_LadrunoJ2Finite=33012`).
2. `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — `#include` + `case … return new`
   (the real switch, not the `FEM_ObjectBroker.cpp` stub).
3. `SRC/interpreter/OpenSeesElementCommands.cpp` — forward decl + 2 `functionMap.insert`:
   `"LadrunoKinematicCoupling"` + the **leading-lowercase camelCase** alias
   `"ladrunoKinematicCoupling"` (review I2 — NOT all-lowercase; keep interior CamelCase, per
   the RBE3 precedent at lines 679-680). Vanilla edit.
4. `SRC/element/CMakeLists.txt` — `add_subdirectory(ladrunoKinematicCoupling)`. Vanilla edit.
   (No kernel include needed — RBE3's CMakeLists never explicitly pulled in `LadrunoEmbeddedKernel`;
   it rides every ladruno dir's `target_include_directories` PUBLIC export. RBE2 simply omits the
   `#include <LadrunoEmbeddedKernel.h>`. Confirmed by integration review.)
5. **`tests/test_ladrunoKinematicCoupling_element.py` must EXIST** (review I1) — it is the
   gated artifact the manifest `pytest:` field points at; without it the Axis-1 manifest gate
   fails. The §8 battery is its content.
6. `Ladruno_scripts/stamp_headers.py` — add a new **dir-glob line pair**
   `"SRC/element/ladrunoKinematicCoupling/*.cpp"`, `"…/*.h"` to GLOBS (review I3 — GLOBS is
   per-directory, not individual files) + rerun (4-author header).
7. `Ladruno_scripts/banner_features.txt` — one line + `python patch_banner.py` + rebuild
   (auto-regenerates both `PythonModule.cpp` and `tclMain.cpp` FEATURES blocks — no hand-edit).
8. `Ladruno_implementation/LEDGER_implementations.md` — one feature row.
9. `Ladruno_implementation/LEDGER_vanilla_files.md` — rows for items 1–4 (`// Ladruno …` marks).
10. `Ladruno_implementation/testbed/manifest.yaml` — full row (element feature_class), with
    `pytest: tests/test_ladrunoKinematicCoupling_element.py` and `tcl: WAIVED(...)` (review I5
    — the shared `functionMap` serves Tcl+Python; verified via the Zone-A pytest battery, per
    the RBE3 row's WAIVE).
11. `Ladruno_implementation/24_ladruno_coupling_constraints_adr.md` §D3 — mark RBE2 **built**;
    optional `LadrunoKinematicCoupling_guide.md`.

> [!note] classTag 33012 is next-*usable* (review I4): 33009/33010 are **reserved** (VEM/SBFEM),
> not free — do not backfill them. The ELE registry jumps `LadrunoCST=33008 → RBE3=33011 → 33012`.

### 7.4 Proposed command
```
element LadrunoKinematicCoupling $tag $refNode $N $s1 ... $sN
        [-dof $c1 ...]                           # dependent components on each slave (default: all available)
        [-k {$Kt | auto}] [-kAlpha $a] [-host $eleTag]   # auto only with -host
        [-kr $Kr]                                # rotational-tie penalty (default K_t*ell^2)
        [-enforce penalty|al]
        [-bipenalty {-dtcr $dt | -wcap $beta}]
        [-absolute]                              # opt out of g0 initial-gap capture
```
**No `-w`/`-area`** (RBE2 has no weights). **`-dof` is the new flag** vs the RBE3 parser; the
parser **range-validates** each component (`1 ≤ c ≤ ndm+nrot`) and rejects an empty list
(review D4). **`-bipenalty` defaults OFF** for RBE2 (review E1) — opt-in, and auto-enabled
only when a massless tied DOF is detected under an explicit integrator.

### 7.5 Serialization
Serialize the `-dof` component list, ref-dof selection, `enforce`/bipenalty/g0 state
(EmbeddedNode/RBE3 hdr+payload pattern) + a header version field. **Recompute `d_i`** on
`recvSelf` from `getCrds()` (setDomain runs on the recv path) — fewer floats, no drift.

---

## 8. Validation battery (Zone-A unless noted)
1. **Rigid-body**: rigid translation + (small) rotation of R+slaves ⇒ all `g_i=0`, zero force.
2. **Force/moment transfer, OFFSET R + SIGN check** (review M1): unit `F` & `M` at R ⇒
   `Σ_i t_i^t = F`, `Σ_i d_i × t_i^t (+Σ t_i^r) = M`, no spurious moment-from-force AND the
   **sense** of the transport couple is correct (apply pure `M` about one axis at an offset
   R, assert the slave force-couple direction — catches a sign-flipped `transOp`, not just
   magnitude).
3. **`-dof` modes**: `-dof 1 2 3` (translation+transport) vs `1 2 3 4 5 6` (full rigid) ⇒
   slave rotation floats vs locked; assert the operator difference.
4. **Mixed-ndf slave set** (review D2): one 6-DOF + one 3-DOF slave under `-dof 1..6` ⇒
   rotation skipped+warned on the 3-DOF slave, ragged gap layout assembles a *correct* moment
   reaction (the indexing-bug regression). Plus the all-3-DOF `LadrunoBrick` face: moment as a
   force couple, no rotational stiffness leaked into the solid.
5. **Reduce to `rigidLink`**: N=1 + translation DOFs ⇒ matches OpenSees `rigidLink -bar`.
6. **AL** (implicit): gap → 0 as Uzawa iterates.
7. **bipenalty + massless-DOF scan** (review E1/E2/E3/E4): (a) massless control-node R stable
   below `-dtcr`, unstable above, self-report bounds both classes via `K_r,eff=λ_max(K_{θθ})`;
   (b) **massed R + `-bipenalty` warns and does NOT double-count** R's nodal mass; (c) a
   **massless slave** DOF is detected and bounded (not just R).
8. **Transient Newmark smoke** — the `getDamp` index-landmine regression (RBE3 crashed here).
9. **Serialization**: `sendSelf`/`recvSelf` roundtrip incl. recomputed `d_i` + `-dof` layout.
10. **Degeneracy negatives** (review D1/D4/D6): `-dof 7` rejected at parse; slave==refNode and
    duplicate slave refused; **all-coincident slaves** translation-only ⇒ R rotation left free,
    no `+Inf` critical step, no NaN; all-coincident with a tied rotation ⇒ floored `K_r > 0`.

---

## 9. Adversarial-review log (2026-06-09, 4-lens sweep)

Mechanics / OpenSees-integration / explicit-dynamics / degeneracy. The mechanics **core was
confirmed correct** by independent re-derivation (self-equilibrium, force/moment transfer,
constant PSD tangent, objectivity scoping). The load-bearing findings were on the
**RBE3→RBE2 role inversion**: RBE3's bipenalty/mass story does not survive when R becomes the
(often massed) master and the slaves become the (possibly massless) dependents.

| # | Lens | Severity | Finding | Resolution in this ADR |
|---|---|---|---|---|
| E1 | explicit | **CRITICAL** | bipenalty unconditionally lumps `m_p`/`I_p` on R, but RBE2's R is often a **real massed node** (platen/footing/equipment) → additive mass double-count = correctness bug | §6: `-bipenalty` **default OFF**; query `Node::getMass()`, lump only on massless DOFs, warn if R massed |
| D1 | degeneracy | **CRITICAL** | all `d_i=0` ⇒ default `K_r=K_t·ℓ²=0` (zeros a tied rotation) and `2√(I_p/K_r,eff=0)=+Inf` | §3 + §6: floor `ℓ²=max(Σ|d_i|²/N, L_char²)`; guard the divide; don't lump `I_p` on an unconstrained rotation |
| E2 | explicit | **MAJOR** | massless **slave** DOF missed by R-only bipenalty → silent CFL filtering | §6: scan every tied DOF of R **and all slaves**; lump where `getMass`=0 |
| E3 | explicit | **MAJOR** | scalar `K_r,eff≈K_t·Σ|d_⊥|²` non-conservative; CFL needs `λ_max(K_{θθ})` | §6: form 3×3 `K_{θθ}`, run `jacobi3`, take `λ_max` |
| E4 | explicit | **MAJOR** | inherited sizing (`iPenalty`↦`Kr`, `mPenalty`↦`K_t`) makes self-report a tautology | §6: size rotation vs `K_r,eff`, R-translation vs `N·K_t` |
| M1 | mechanics | **MAJOR** | "reuse `transOp` with a `+` sign" mislabels the needed **sign FLIP** (`B += −transOp`) → flipped transport block if copied verbatim | §2.1/§2.4: explicit `B += −transOp`, 2D B-column `[d_y,−d_x]ᵀ`; §8 case 2 asserts the sign |
| D2 | degeneracy | **MAJOR** | ragged per-slave gap layout (mixed 3/6-DOF) → mis-indexed moment reaction (silent) | §3: mandate `gapRowStart/nTiedTrans/nTiedRot` layout table; §8 case 4 mixed-ndf test |
| D4 | degeneracy | **MAJOR** | `-dof` component range unvalidated (`0`,`7`,neg) → OOB/typo masking | §2.3/§7.4: **parse-time** range check `1≤c≤ndm+nrot` |
| D6 | degeneracy | **MAJOR** | duplicate slave / slave==refNode unchecked → over-stiffen / self-tie (silent) | §3: REFUSE both (inherited RBE3 gap — fix both) |
| I1 | integration | MAJOR(bk) | manifest `pytest:` path must point at a real test file | §7.3 item 5: test file is a gated DoD artifact |
| D3 | degeneracy | MINOR | R-without-rotation silently drops a user's explicit rotation `-dof` | §3: degrade only if defaulted; **REFUSE** if rotation explicitly requested |
| D5 | degeneracy | MINOR | empty-`-dof` refusal point | §2.3: syntactic-empty at parse, effective-empty inert-refuse at setDomain |
| D7 | degeneracy | MINOR | 2D drilling component numbering undocumented | §2.3: explicit `1,2`=trans, `3`=drilling + 2D test |
| I2 | integration | MINOR | alias is **leading-lowercase** camelCase, not all-lowercase | §7.3 item 3 |
| I3 | integration | MINOR | `stamp_headers.py` GLOBS is a per-dir glob, add a line pair | §7.3 item 6 |
| I4/I5 | integration | MINOR | 33009/33010 reserved-not-free; record `tcl: WAIVED` | §7.3 note + item 10 |
| — | mechanics | OK | self-equilibrium, F/M transfer, `K=ΣBᵀDB` PSD/constant, objectivity scoping: **CONFIRMED** | §2.4, §4 (verified, keep) |
| — | explicit | OK | `getDamp`+`getRayleighDampingForces` override mandate **complete**; AL-under-explicit warn correct | §6 (keep) |
| — | integration | OK | broker file, command-dispatch pattern, CMake (no kernel needed), parser signature: **CONFIRMED** | §7.3 (verified) |

**Net:** the spec is mechanically sound; the explicit-dynamics section was substantially
hardened (the bipenalty default flipped to OFF, the massless-scan generalized to all nodes,
`K_r,eff` upgraded to `λ_max`, sizing corrected) and the degeneracy table grew from 5 to ~9
rows. Locked for build.

---

## 10. References
- **ADR 24** §2 (Abaqus/LS-DYNA/Nastran alignment), §2.1 (LS-DYNA explicit constraint
  treatment), D3 (kinematic coupling / RBE2), D5 (`dt` knobs). **ADR 28 / RBE3** (the
  sibling + the reused scaffold, `transOp`, the `getDamp` landmine §8b). **ADR 23 /
  EmbeddedNode** + **ADR 20 / EmbeddedRebar** (the scaffold lineage).
- **Nastran** RBE2 (rigid body, dependent components). **Abaqus** `*COUPLING, kinematic` /
  `*KINEMATIC COUPLING`. **LS-DYNA** `*CONSTRAINED_NODAL_RIGID_BODY` (CNRB; Theory §25 mass
  condensation — the deferred D3b alternative). **Felippa IFEM** (MPC/constraint chapters).
  de Souza Neto/Perić/Owen Ch. 2 (linearisation).
