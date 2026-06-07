---
title: "ADR 28 — LadrunoDistributingCoupling (RBE3 / distributing coupling): implementation spec"
project: Ladruno
type: ADR / implementation spec (pre-code)
status: planned — formulation locked, adversarially reviewed; ready to build
related:
  - "[[24_ladruno_coupling_constraints_adr]]"   # parent: the coupling-constraint family
  - "[[23_ladruno_embedded_node_adr]]"          # sibling: node-into-host (UR ≠ RBE3 rotation)
  - "[[20_ladruno_embedded_reinforcement_adr]]" # the shared penalty/AL/bipenalty kernel
  - "[[ndf_and_mixed_models_guide]]"            # the ndf-mismatch driver
  - "[[LEDGER_implementations]]"
tags: [adr, constraints, coupling, rbe3, distributing-coupling, ndf-mismatch, moment-transfer, explicit-dynamics, abaqus, ls-dyna, nastran]
updated: 2026-06-07
---

# ADR 28 — `LadrunoDistributingCoupling` (RBE3 / distributing coupling)

**Status:** planned. Formulation locked and **adversarially reviewed** (4-lens sweep:
mechanics / OpenSees-integration / explicit-dynamics / degeneracy, 2026-06-07).
This is the implementation spec; code has not landed. It is the first concrete build
under [[24_ladruno_coupling_constraints_adr|ADR 24]] (the coupling-constraint family),
triggered by the **ndf-mismatch moment-transfer driver** ([[ndf_and_mixed_models_guide]] §4).

> [!info] What RBE3 is, in one line
> A **reference node R** whose motion is the *weighted least-squares fit* of a set of
> **independent nodes** — equivalently, a load applied at R is **distributed** to the
> independents as a statically-equivalent force set, **adding no stiffness** to them.
> = Abaqus `*COUPLING, distributing` = LS-DYNA `*CONSTRAINED_INTERPOLATION` = Nastran RBE3.
> Contrast RBE2 (kinematic coupling): one reference rigidly *drives* a slave set (adds
> rigid stiffness). **In RBE3 the single node is the dependent/slave; in RBE2 it is the
> master.** The roles flip — the classic wiring mistake.

---

## 1. Driver & goal

The concrete need is the **ndf-mismatch interface**: a 6-DOF beam/shell reference node
framing into a face/ring of 3-DOF solid nodes, transmitting **moment into the continuum**.
A translation-only tie (`equalDOF_Mixed`, the embedded U-tie) lets the beam node spin
free → the moment leaks ([[ndf_and_mixed_models_guide]] §4). RBE3 transmits the moment
**correctly**: it becomes a self-equilibrated *force couple* over the solid face (how a
continuum physically carries a moment — as distributed traction), with **no rotational
DOF or rotational stiffness injected into the solid** and the **face free to deform**
(unlike RBE2, which would rigidify the patch).

**v1 scope (locked):** full **U + rotation** distributing coupling in one element,
2D + 3D, penalty + augmented-Lagrangian, explicit-safe via bipenalty. Arbitrary
reference-node location (kinematically-complete transport term). classTag **33011**.

**Non-goals (v1):** finite-rotation configuration update (geometric stiffness term —
see §4.6); mass redistribution to independents (D2 — deferred, see §6.4); RBE2 kinematic
coupling and the general N-node linear equation (separate ADR-24 build items).

---

## 2. Formulation (locked, adversarially verified)

### 2.1 Notation

Reference node **R**: translations `u_R`, rotations `θ_R`. Independent nodes `I_i`,
`i=1..N`: translations `u_i`, user weights `w_i > 0`. Position vectors `x_R`, `x_i`.

- `W = Σ w_i`
- **weighted centroid** `x_c = (1/W) Σ w_i x_i`  — *always the origin for the fit, never `x_R`*
- `r_i = x_i − x_c`  ⇒  `Σ w_i r_i = 0` (the property that decouples translation/rotation)
- **weighted position-inertia tensor**
  - direct:  **I_c** = Σ w_i [ (r_i·r_i) **1** − r_i ⊗ r_i ]
  - index:   (I_c)_{jk} = Σ w_i [ (r_i·r_i) δ_{jk} − r_{i,j} r_{i,k} ]
  - this is a *geometric* tensor (a property of the point cloud + weights), **not** a
    penalty stiffness and **not** a thin-plate inertia (see §3 — flat faces are full-rank).

`nrot = 3` (3D) or `1` (2D, drilling). `[r_i]_×` = the skew/cross-product matrix
(`[r_i]_× v = r_i × v`); in 2D the `1×2` row `[ −r_{i,y},  r_{i,x} ]`.

### 2.2 Kinematic constraint (the relation RBE3 enforces)

R's motion is the **weighted least-squares rigid-body fit** of the independent
translation field. Minimising `Φ = Σ w_i | u_i − (a + ω × r_i) |²` over `(a, ω)` gives
(verified by hand, mechanics review Finding 4 — both decouplings hinge on `Σ w_i r_i = 0`):

- translation:  `a = (1/W) Σ w_i u_i`           (weighted centroid translation)
- rotation:     `ω = I_c⁻¹ Σ w_i (r_i × u_i)`   (moment-balanced LSQ rotation)

**Constraint gaps** (kinematically complete for an *offset* reference node — decision
locked; mechanics Finding 1 + degeneracy Task 7):

```
g_t = ( u_R + θ_R × (x_c − x_R) ) − (1/W) Σ w_i u_i          (size ndm)
g_r =   θ_R − I_c⁻¹ Σ w_i (r_i × u_i)                         (size nrot)
```

The `θ_R × (x_c − x_R)` term is the rigid link from R to the centroid; dropping it
(the naïve `g_t = u_R − Σ(w_i/W)u_i`) is silently wrong whenever `x_R ≠ x_c` and force
and moment are transmitted together — the most likely failure vs Nastran/Abaqus. Keeping
`r_i` measured from `x_c` preserves the clean decoupled `I_c`; the offset enters *only*
through this one extra block in `g_t`.

### 2.3 Penalty enforcement & the work-conjugate force split

Penalty potential `Π_p = ½ K_t |g_t|² + ½ K_r |g_r|²`; tractions `t_t = K_t g_t`
(+`λ_t` for AL), `t_r = K_r g_r` (+`λ_r`). Resisting force `P = Bᵀ t`, tangent
`K = Bᵀ D B` with `D = diag(K_t I, K_r I)`.

**The B-operator** (constant in the reference configuration):

| ∂/∂ | `g_t` (trans) | `g_r` (rot) |
|---|---|---|
| `u_R`  | `I`                       | `0` |
| `θ_R`  | `[ (x_c − x_R) ]_×ᵀ`      | `I` |
| `u_i`  | `−(w_i/W) I`              | `−I_c⁻¹ (w_i [r_i]_×)` |

The force delivered to independent node `i` (= `−` its row of `Bᵀt`) is the **textbook
RBE3 split** (verified: full force *and* moment balance for **any** K — load conservation
is exact regardless of penalty magnitude, mechanics Findings 5 & 7):

```
f_i = (w_i/W) t_t          ← force part:  Σ f_i^F = t_t,  Σ r_i×f_i^F = 0
    + w_i (α × r_i),  α = I_c⁻¹ t_r   ← moment part: Σ r_i×f_i^M = t_r, Σ f_i^M = 0
```

### 2.4 The rotation block is NEW — the shared kernel CANNOT assemble it

**Critical reuse boundary** (explicit + integration reviews). The `g_r` row couples
`θ_R` to the independents' *translations* through lever arms `I_c⁻¹ w_i [r_i]_×` — a
**cross-DOF-class block**. `LadrunoEmbeddedKernel::assembleTangent`'s per-node scalar
form `block(p,q) = c_p c_q D` **structurally cannot** express it. Therefore:

- **Translation block** → reuse `LadrunoEmbeddedKernel` *unmodified* (it is already
  weight-agnostic — never assumes `Σw=1`; integration review verified, no risk to the
  rebar/node callers). Assemble via the **EmbeddedNode scatter** (R is 6-DOF), not the
  rebar's direct full-matrix write.
- **Rotation block + the `θ_R↔u_i` and `θ_R↔u_R` transport blocks** → a **new
  RBE3-specific assembler** in the element. Needs `getCrds()`, `I_c`, its pseudo-inverse,
  and the lever-arm B. This is the genuinely new code; "~70% kernel reuse" (ADR 24 §6)
  holds for *everything except* this load-bearing piece.

> [!warning] Do NOT reuse the EmbeddedNode UR path for `g_r`
> EmbeddedNode UR is `θ = ½ curl(u) = skew(∇u)` from host shape-function **gradients**
> (a continuum spin). RBE3's rotation is a least-squares fit from nodal **positions**
> `r_i` and `I_c`. Different operators — must not be conflated (integration Finding 3,
> confirmed: EmbeddedNode has no `getCrds`/`I_c`/pseudo-inverse).

---

## 3. Degeneracy of `I_c` — flat face is FINE; collinear/N=1 are not

The plan's original worry ("flat face is degenerate") was a **false premise** — debunked
numerically (degeneracy review Task 2). `I_c` is a *position*-inertia: for a square face
(4 nodes, side 2, unit weights) `I_c = diag(4,4,8)`, eigenvalues `(4,4,8)` — the
normal/drilling axis is the **largest** (`= Σw_i|r_i|²`), all three reference rotations
resolvable, `cond = 2`. **The main use case (beam node on a 2D-spread flat face) is
perfectly conditioned.** Regularization is needed only for:

| Case | `I_c` | Action |
|---|---|---|
| **Coplanar, 2D spread** (the flat face) | full rank | nothing special |
| **Collinear** (edge of nodes, N=2) | exactly **one** zero eig (the line axis) | drop that rot axis, leave free, warn |
| **N = 1** | `I_c = 0` (all zero) | tie translation only; drop all rotation; warn |
| **Near-collinear strip** | one tiny eig | drop below floor (not Tikhonov) |
| **W ≤ 0 / negative weights** | `x_c` = NaN / `I_c` indefinite | **refuse** at setDomain |

### 3.1 Locked regularization algorithm (degeneracy review, recommendation (b))

At `setDomain` (resolve **once**, freeze for the analysis — coords are known here):

```
1. W = Σ w_i ;  REFUSE if any w_i < 0 or W ≤ 1e-12·max w_i ;  WARN any w_i == 0 (inactive).
2. x_c = Σ w_i x_i / W ;  r_i = x_i − x_c .            # origin = x_c ALWAYS, never x_R.
3. if N == 1:  translation-only; drop ALL rotation rows; WARN; return.
4. I_c = Σ w_i [ (r_i·r_i)1 − r_i⊗r_i ] .
5. L² = max_i |r_i|² .                                  # characteristic length²
   symmetric eigensolve  I_c = Σ λ_k v_k v_kᵀ .         # 3×3 SPD/PSD (Jacobi/dsyevd)
6. λ_floor = max( 1e-8·λ_max ,  1e-10·W·L² ) .          # unit-invariant (I_c scales as L²)
7. keep_k = (λ_k > λ_floor) ;
   I_c⁺ = Σ_{keep} (1/λ_k) v_k v_kᵀ ;  P = Σ_{keep} v_k v_kᵀ .   # spectral pseudo-inverse
   for each dropped k: WARN "reference rotation about axis v_k unconstrained"; leave free.
8. θ_R^t = I_c⁺ Σ w_i (r_i × u_i) ;  enforce only P·g_r .
9. capture g0 / gr0 in the SAME basis (P·gr0).          # consistency with the drop (§5)
2D: replace 5–8 with scalar I_c = Σ w_i|r_i|² ; drop drilling iff I_c < 1e-10·W·L².
```

**Spectral drop, NOT Tikhonov** — verified numerically: dropping projects noise out
(spurious θ → 0); Tikhonov leaves residual noise *and* biases the well-conditioned axes
(pollutes good DOFs to fix a bad one). The drop scheme also **matches Nastran/Abaqus**,
which leave an unresolvable reference rotation free + warn rather than inverting `I_c`.
Freeze the eigenbasis + keep-mask at setDomain so a borderline eigenvalue can't flip
across `λ_floor` mid-run (→ discontinuous stiffness / non-convergence).

---

## 4. Consistent tangent & objectivity

### 4.1 Reference-configuration tangent (v1)
`B` is constant (weights + reference geometry + frozen `I_c⁺`), `t = D g` is linear in
the DOFs, so `K = Bᵀ D B` is the **exact penalty Hessian**, symmetric PSD, and constant.
`getInitialStiff = getTangentStiff` (resolve `I_c` eagerly at setDomain so it is valid).

### 4.2 Objectivity boundary (documented, not a bug)
Under rigid translation of the whole assembly, `g_t → 0` exactly. Under a **finite**
rigid rotation, `g_r` vanishes only to **first order** (`O(φ_rel²)`, where `φ_rel` is R's
rotation *relative to the set's frame*) because `θ_R` is a finite nodal-rotation DOF while
the LSQ fit is linear in `u_i`. Acceptable for the small-relative-rotation beam/shell-on-
solid driver. **Finite-rotation support is a future geometric-stiffness term** (`∂Bᵀ/∂x · t`
from co-rotating `r_i`), *not* a parser tweak — the "exact & constant tangent" claim is
explicitly scoped to the reference configuration (mechanics Findings 2 & 5).

---

## 5. Explicit-dynamics safety (the fork's center of gravity)

RBE3's topology is **inverted** vs the embedded family: the reference node R is the
**massless interpolation point** (not the host). Treatment (explicit review):

- **R maps to slot 0** (the `+I` block) so `getMass` lumps the bipenalty `m_p` (trans)
  and `I_p` (rot) on R, and `getExplicitCriticalTimeStep` self-reports the bound.
- **Derived rotational stiffness** `K_r,eff = K_t · Σ (w_i/W) |r_i⊥|²` (the lever-weighted
  translational penalty; same `ℓ²` that fixes the K_t/K_r units mismatch, ℓ²≈Σw|r|²/W).
  `I_p` is lumped against **`K_r,eff`**, and the self-report takes
  `min( 2√(m_p/K_t), 2√(I_p/K_r,eff) )`. **Without a stored `K_r,eff` the rotation mode
  is unbounded and unreported → `CriticalTimeStep` silently filters R's massless rotation
  ([[ndf_and_mixed_models_guide]] §7) → unstable, no diagnostic.** This is the #1 explicit
  must-fix.
- **bipenalty is the explicit default**; mass redistribution (D2) is **deferred** (more
  CFL-hazardous — per-element estimate becomes inconsistent with the assembled M, plus
  double-count risk). If both are ever enabled, they must be **mutually exclusive**.
- **AL is implicit-only**: Uzawa updates `λ` in `commitState`, which under explicit
  central-difference fires every sub-step with no equilibrium iteration → `λ` ratchets on
  the un-converged gap, giving no accuracy benefit. **Warn if `-enforce al` meets an
  explicit integrator.** (The EmbeddedNode has *no* AL/bipenalty interlock to inherit — a
  draft premise was wrong.)
- Refuse Rayleigh factors (pure coupling carries no physical damping; a `βK` would
  spuriously shrink `dt_cr`) — reuse the EmbeddedNode refusal.
- **`-k auto` is undefined for a node set** (no host element). Default to numeric `-k`;
  `-k auto` only with an explicit representative `-host` element among the independents'
  parents. Flagged design gap; do not silently pick "one independent's parent".

---

## 6. Implementation map

### 6.1 Files (new element dir, mirrors EmbeddedNode)
```
SRC/element/ladrunoDistributingCoupling/
  LadrunoDistributingCoupling.{h,cpp}     # element: DOF scatter, I_c resolve, new rot assembler
  OPS_LadrunoDistributingCoupling.cpp     # parser
  CMakeLists.txt                          # + include path ../ladrunoEmbeddedRebar (shared kernel)
tests/test_ladrunoDistributingCoupling_element.py
```

### 6.2 Reuse (copy/consume from the embedded family)
- **Consume** `SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.{h,cpp}` for the
  **translation** block (unmodified — weight-agnostic).
- **Copy the scaffold** from `LadrunoEmbeddedNode.{h,cpp}`: per-node DOF scatter
  (`nodeNdf`/`dofOffset`/`nDOF`/`allocate`), penalty/AL (`lambda`, Uzawa in `commitState`),
  bipenalty (`mPenalty`/`iPenalty`, `getMass` diagonal-on-slot-0, `setRayleighDampingFactors`
  refusal, `getExplicitCriticalTimeStep`), g0 initial-gap capture, `sendSelf`/`recvSelf`
  header-vector pattern, `setResponse`/`getResponse`.
- **NEW**: `getCrds`-based `I_c` assembly + spectral pseudo-inverse (§3.1); the rotation
  + transport B blocks and their tangent (§2.3–2.4); the derived `K_r,eff` (§5).

### 6.3 Plumbing (Definition-of-Done — verified file:line by integration review)
1. `SRC/classTags.h` — `#define ELE_TAG_LadrunoDistributingCoupling 33011` (after CST=33008;
   33009/33010 left to the VEM/SBFEM paper reservations).
2. `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` — `#include` + `case … return new …()`
   (the real switch; **not** the `FEM_ObjectBroker.cpp` stub).
3. `SRC/interpreter/OpenSeesElementCommands.cpp` — forward decl + 2 `functionMap.insert`
   (CamelCase + lowercase alias). Vanilla edit.
4. `SRC/element/CMakeLists.txt` — `add_subdirectory(ladrunoDistributingCoupling)`. Vanilla edit.
5. `Ladruno_scripts/stamp_headers.py` — add new files to GLOBS + rerun (4-author header).
6. `Ladruno_scripts/banner_features.txt` — one line + `python patch_banner.py` + rebuild.
7. `Ladruno_implementation/LEDGER_implementations.md` — one feature row.
8. `Ladruno_implementation/LEDGER_vanilla_files.md` — rows for items 1–4 (`// Ladruno …` marks).
9. `Ladruno_implementation/testbed/manifest.yaml` — full row (element feature_class).

### 6.4 Serialization
Serialize weights `w_i`, ref-dof selection, `enforce`/bipenalty/g0 state (EmbeddedNode
hdr+payload pattern) + a **header version field**. **Recompute `I_c` (and `r_i`, the
eigenbasis/keep-mask) on `recvSelf`** from `getCrds()` (setDomain runs on the recv path) —
fewer floats, no drift.

### 6.5 Proposed command
```
element LadrunoDistributingCoupling $tag $refNode $N $i1 ... $iN
        [-w $w1 ... $wN | -area]                # weights (default equal); -area = tributary
        [-refdof $d1 ...]                        # which ref DOFs are dependent (default all 6/3)
        [-k {$Kt | auto}] [-kAlpha $a] [-host $eleTag]   # auto only with -host
        [-enforce penalty|al]
        [-bipenalty {-dtcr $dt | -wcap $beta}]
        [-absolute]                              # opt out of g0 initial-gap capture
```

---

## 7. Validation battery (Zone-A unless noted)

1. **Rigid-body**: rigid translation + (small) rotation of ref+set ⇒ `g_t=g_r=0`, zero force.
2. **Force/moment balance with an OFFSET ref** (`x_R ≠ x_c`): unit `F` & `M` at R ⇒
   `Σ f_i = F`, `Σ (x_i−x_c)×f_i = M_c = M + (x_R−x_c)×F`, **no spurious moment-from-force**
   (the transport-term regression — §2.2).
3. **vs Nastran RBE3 / Abaqus distributing** on a canonical load-redistribution patch.
4. **Degeneracy**: collinear set + N=1 ⇒ warn, free axis, **no NaN**; flat face ⇒ full rank
   (assert all 3 rotations resolved).
5. **ndf-mismatch end-to-end**: 6-DOF beam node on a 3-DOF solid face (LadrunoBrick) ⇒
   moment transmitted, face deforms, no rotational stiffness in the solid.
6. **AL** (implicit): gap → 0 as Uzawa iterates.
7. **bipenalty**: `getExplicitCriticalTimeStep` self-report bounds **both** classes; a
   massless ref runs stable below the `-dtcr` target, unstable above.
8. **Serialization**: `sendSelf`/`recvSelf` roundtrip (incl. recomputed `I_c`).

---

## 8. Adversarial-review log (2026-06-07, 4-lens sweep)

| # | Lens | Severity | Finding | Resolution in this ADR |
|---|---|---|---|---|
| 1 | mechanics | CRITICAL | R≠centroid drops the transport term | §2.2 kinematically-complete gap (locked) |
| 2 | mechanics | MAJOR | "exact & constant tangent" overclaimed | §4 scoped to reference config |
| 3 | mechanics | MAJOR | K_t/K_r units mismatch → conditioning | §5 `K_r,eff = K_t·ℓ²` |
| 4 | explicit | CRITICAL | massless ref node topology inverted | §5 R→slot 0, m_p/I_p on R |
| 5 | explicit | CRITICAL | rotation block uncoupleable by shared kernel | §2.4 new RBE3 assembler |
| 6 | explicit | CRITICAL | rotation mode unbounded (no `K_r`) → §7 trap | §5 derived `K_r,eff`, min-report |
| 7 | explicit | HIGH | `-k auto` undefined for a node set | §5 numeric default, `-host` opt-in |
| 8 | explicit | HIGH | mass redistribution (D2) CFL-hazardous | §6.4 deferred; bipenalty v1 |
| 9 | explicit | HIGH | AL meaningless under explicit | §5 implicit-only + warn |
| 10 | integration | HIGH | classTag 33009 soft-reserved (VEM) | §6.3 take 33011 |
| 11 | integration | HIGH | wrong broker file + missing element-cmd dispatch | §6.3 corrected |
| 12 | integration | MED | resolve I_c eager / recompute on recv | §6.4 |
| 13 | degeneracy | CRITICAL | **flat face is NOT degenerate** (false premise) | §3 (debunked numerically) |
| 14 | degeneracy | HIGH | drop (pseudo-inverse), not Tikhonov | §3.1 |
| 15 | degeneracy | HIGH | W≤0 / negative weights refuse | §3.1 step 1 |
| 16 | mechanics/deg | — | LSQ fit, duality, 2D, load-conservation: **CORRECT** | §2 (verified, keep) |

---

## 9. References
- **ADR 24** §2 (Abaqus/LS-DYNA/Nastran alignment), §2.1 (LS-DYNA explicit constraint
  treatment), D2 (mass redistribution), D5 (`dt` knobs). **ADR 23 / EmbeddedNode** (the
  reused scaffold; UR ≠ RBE3 rotation). **ADR 20 / EmbeddedRebar** (the shared kernel).
- **Nastran** RBE3 (interpolation/distributing). **Abaqus** `*COUPLING, distributing`.
  **LS-DYNA** `*CONSTRAINED_INTERPOLATION` (Keyword Vol I; `DDOF`/`TWGHT`/`RWGHT`, mass
  redistribution). **Felippa IFEM** (MPC/constraint chapters). de Souza Neto/Perić/Owen
  Ch. 2 (linearisation) for the tangent.
