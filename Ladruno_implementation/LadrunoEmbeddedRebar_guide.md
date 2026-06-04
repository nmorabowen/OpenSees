---
title: "LadrunoEmbeddedRebar — Embedded-Reinforcement Coupling Element"
project: Ladruno
type: reference / user guide
status: shipped (Mode P + bond-slip + §10 roadmap: -host/-xi, -kt auto, energy, -corot, -enforce/AL, bipenalty + -cfl seam)
element: element LadrunoEmbeddedRebar
classTag: 33005 (element) · LadrunoBondSlip 33002 (axial law)
related:
  - "[[20_ladruno_embedded_reinforcement_adr]]"
  - "[[19_rc3d_modeling_recipe]]"
  - "[[22_rc3d_conformal_recipe]]"
  - "[[ladruno_integrators_guide]]"
  - "[[central_difference_ladruno_guide]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - guide
  - element
  - embedded-reinforcement
  - bond-slip
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
  - rc3d
updated: 2026-06-04
---

# LadrunoEmbeddedRebar — embedded-reinforcement coupling element

`element LadrunoEmbeddedRebar` (classTag **33005**) ties **one discrete rebar
node** to the nodes of a **solid host element** through precomputed shape-function
weights, so a reinforcement mesh can be embedded in a **non-matching** concrete
mesh. It is a pure **coupling** (bond) element — the bar's own axial stiffness
lives on a separate `corotTruss`/beam along the bar; this element only enforces
the bar↔concrete attachment, with either **perfect bond** or a **1D bond-slip**
law in the axial direction.

This is the single user/developer reference for the element. It covers the
**theory** (the tie, the anisotropic traction, the penalty formulation), the
**capabilities** (the full option grammar and what each switch does), the
**implementation** (where the math lives, serialization, the vanilla seams), and
the **use cases** (host queries via apeGmsh, perfect vs bond-slip, explicit RC
workflows, pitfalls). The design record is
[[20_ladruno_embedded_reinforcement_adr|ADR 20]].

> [!info] Where it fits
> `LadrunoEmbeddedRebar` is the **v2** deliverable of the fork's RC-3D plan: the
> tool for cases the conformal (shared-node) path can't mesh — non-grid bar
> positions, hoops/ties, beam–column joints. For regular ≤8-bar members the
> conformal recipe ([[22_rc3d_conformal_recipe]]) is cheaper and needs no new
> element. Concrete hosts: [[20_ladruno_embedded_reinforcement_adr|LadrunoBrick / BezierTet10]].
> Rebar materials: [[LadrunoUniaxialJ2_guide]], [[LadrunoRebarBuckling_guide]].

---

## 1. Theory — the coupling

### 1.1 The tie and the gap

The element holds external nodes $[\,\text{rebar},\,\text{host}_1,\dots,\text{host}_M\,]$
and enforces the kinematic tie that the rebar node follows the host
interpolation:

$$g(u) = u_r - \sum_{i=1}^{M} N_i(\xi)\,u_i^{host} = 0 ,$$

where $N_i(\xi)$ are the host element's nodal shape-function weights at the rebar
point's natural coordinate $\xi$ (fixed at build time), and $g$ is the
**current-configuration relative displacement** (the *constraint gap*). Because
$g$ and the $N_i$ at a fixed material point are both objective, the tie is
**exactly objective** under rigid host rotation.

### 1.2 Anisotropic traction — axial vs transverse

The gap is split along the **bar axis** $\hat{d}$ (unit vector):

$$s = g\cdot\hat{d}\quad(\text{axial slip}),\qquad
  g_t = g - s\,\hat{d}\quad(\text{transverse gap}),$$

and the coupling traction is

$$t = F_\text{axial}(s)\,\hat{d} + k_t\,g_t ,$$

with the **axial law** $F_\text{axial}$ either perfect-bond ($k_\text{axial}\,s$)
or a bond-slip stress, and $k_t$ the **transverse penalty**. The split is
deliberate (ADR D6): the bar resists bond/slip along its axis but is held
transversely by penalty — a bare `corotTruss` is rank-1 ($\propto n\otimes n$),
so a bar parallel to a host face would leave a near-floating transverse DOF
without $k_t$.

### 1.3 Resisting force and tangent

With the discrete operator $B = [\,I \mid -N_1 I \mid \dots \mid -N_M I\,]$
(so $g = B\,u$),

$$P = B^\mathsf{T} t,\qquad K = B^\mathsf{T} D B,\qquad
  D = \frac{\mathrm{d}F_\text{axial}}{\mathrm{d}s}\,\hat{d}\otimes\hat{d}
      + k_t\,(I - \hat{d}\otimes\hat{d}).$$

Block-wise, the rebar block of $P$ is $+t$ and host block $i$ is $-N_i t$; the
stiffness blocks are $c_p c_q D$ with $c = [1, -N_1, \dots, -N_M]$. This is the
anisotropic operator of ADR D2 (**not** an isotropic $iK\,B^\mathsf{T}B$, which
would violate virtual-work symmetry).

### 1.4 Mode P (penalty) — why penalty, and what's deferred

The element is **Mode P**: it enforces the tie by **penalty** — it adds
**stiffness only**, so it works in **both implicit and explicit** (it never
densifies the mass matrix). Perfect bond is Mode P with a large axial penalty;
bond-slip plugs a τ–s law into the axial slot.

> [!note] Mode T (transformation / DOF-elimination) is deferred indefinitely
> A perfect-bond tie could also be imposed by **eliminating** the rebar DOFs via a
> constraint $u_r = \sum_i N_i u_i^{host}$ (no penalty). That needs a
> **multi-retained** `MP_Constraint` the stock single-retained object cannot
> express, and it is implicit-only (it densifies $T^\mathsf{T} M T$ → wrong inertia
> in a `DiagonalSOE`). All forward development is on the Mode P
> `-enforce` family. `-enforce transformation` is a hard parse-time error.

### 1.5 Host-agnostic interpolation

The element stores the weights $N_i$ and never needs to know the host element
type. The weights come from one of:

- the **apeGmsh generator** (precomputed from the host type + reference-config
  natural coords), or
- a **host query**: `-xi ξ η ζ` calls the host's
  `Element::getInterpolationWeights` (a fork-added vanilla virtual, [[LEDGER_vanilla_files]]),
  overridden by `LadrunoBrick` (trilinear) and `BezierTet10` (quadratic
  Bernstein). v1 supports **straight-sided** hosts; curved-host inverse maps are
  deferred.

### 1.6 The rebar element itself (truss vs beam)

The embedding couples **translations only**, so the *rebar* element type is the
user's choice on the rebar mesh (ADR D6, the Abaqus `*EMBEDDED ELEMENT` model):

- **`truss` (default, `corotTruss`)** — axial only, cheapest, explicit-clean. Its
  transverse stiffness comes from the Mode-P $k_t$ tie, not the element.
- **`beam` (opt-in)** — adds the bar's distributed flexural/transverse stiffness
  (an *approximate*, mesh-limited dowel contribution). Requires a bar-axis twist
  guard (spurious zero-energy mode) and no double-count with bond-slip.

---

## 2. Capabilities — the option grammar

```
element LadrunoEmbeddedRebar tag rebarNode
        {nHost h1..hN | -host eleTag}          # host nodes: explicit or by element
        {-shape N1..NN | -xi x1..x_ndm}        # weights: explicit or queried from host
        -dir dx dy [dz]                        # bar axis (reference)
        ( -bond matTag [-bondScale bs] | -perfect kAxial )   # axial law
        [-kt {kt | auto}] [-ktAlpha a]         # transverse penalty (numeric or auto-scaled)
        [-corot {-xiB b.. | -shapeB N..}]      # co-rotated bar axis (large rotation)
        [-enforce {penalty | al}]              # constraint treatment
        [-bipenalty {-dtcr dt | -wcap β}]      # explicit critical-time-step control
```

| Switch | Meaning | Default | PR |
|--------|---------|---------|----|
| `nHost h1..hN` / `-host eleTag` | host nodes explicit, or pulled from a host element | — | #169 / #175 |
| `-shape N1..NN` / `-xi x..` | weights explicit, or queried via `getInterpolationWeights` | — | #169 / #175 |
| `-dir` | reference bar axis $\hat{d}$ | — | #169 |
| `-perfect kAxial` | perfect-bond axial penalty | — | #169 |
| `-bond matTag [-bondScale bs]` | `LadrunoBondSlip` τ–s law in the axial slot; `bs = perimeter·L_trib` | bs=1 | #169 |
| `-kt {kt\|auto}` | transverse penalty, numeric or host-stiffness-scaled | 1e12 | #169 / #177 |
| `-ktAlpha a` | multiplier for `-kt auto` ($k_t=a\cdot\max\lvert K_\text{host}\rvert$) | 1e3 | #177 |
| `-corot -xiB../-shapeB..` | co-rotate the bar axis each step | off (frozen) | #180 |
| `-enforce {penalty\|al}` | penalty (Mode P) or augmented Lagrangian | penalty | #181 |
| `-bipenalty {-dtcr dt\|-wcap β}` | explicit critical-time-step control | off | #186 / #190 |

---

## 3. Axial law — perfect bond vs bond-slip

### 3.1 `-perfect kAxial`

A linear axial penalty $F_\text{axial} = k_\text{axial}\,s$. Use a large
`kAxial` to approximate perfect bond. The artificial axial energy
$\tfrac12 k_\text{axial} s^2$ is reported in `penaltyEnergy` (it is not physical
dissipation).

### 3.2 `-bond matTag` (LadrunoBondSlip) + `-bondScale`

A physical bond-slip law: $F_\text{axial} = \tau(s)\cdot(\text{bondScale})$, where
`matTag` is a `LadrunoBondSlip` uniaxial (classTag **33002**, CEB-FIP / MC2010
monotonic backbone) and **`bondScale = perimeter · L_trib`** converts the τ–s
*stress* law into a nodal *force* (perimeter $\pi d_b$ × tributary length).

`LadrunoBondSlip` backbone: a linear segment for $s<s_0$ (kills the power-law
$\mathrm{d}\tau/\mathrm{d}s\to\infty$ singularity), an ascending power law
$\tau_\text{max}(s/s_1)^\alpha$, a plateau, softening to $\tau_f$, residual
$\tau_f$; `-Gf` regularizes the softening triangle by a bond fracture energy. See
ADR D4. The physical bond work is reported in `bondEnergy` (single-sourced from
the material's own `energy` response).

> [!warning] Softening needs displacement/arc/IMPLEX control
> Past $\tau_\text{max}$ the backbone has a **negative tangent** — load control
> diverges. Use DisplacementControl / ArcLength / IMPLEX on the softening branch.

---

## 4. Transverse penalty & auto-scaling (`-kt`, `-kt auto`)

The transverse hold $k_t$ can be a raw number (`-kt 1e8`) or **auto-scaled** to
the host's own stiffness (`-kt auto`):

$$k_t = a \cdot \max_i \lvert K_\text{host}(i,i)\rvert \;\approx\; a\cdot E_\text{host}\cdot\ell_\text{ch},
\qquad a = \texttt{-ktAlpha}\ (\text{default } 10^3).$$

This makes the constraint conditioning **mesh/material-independent**: the
constraint residual $\lVert g_t\rVert \sim O(\sigma/(aE))$ shrinks while
$\kappa(K)\sim a\,\kappa(K_\text{host})$ stays bounded. `-kt auto` requires the
`-host` form (it reads the host's `getInitialStiff`) and resolves lazily on first
assembly. Floor $a\gtrsim 10$ (below it the bar floats transversely); ceiling
$a\lesssim 10^4$ for implicit Krylov solvers.

> [!tip] Porting from ASDEmbeddedNodeElement
> ASD's `m_K = 1e18` is **not** $E$-scaled — do **not** pass `1e18` as `-ktAlpha`
> (condition-number blow-up). `-kt auto -ktAlpha 1e3` is the conditioned default.

---

## 5. Enforcement strategies (`-enforce`)

### 5.1 `penalty` (Mode P, default)

The shipped kernel: traction $t = F_\text{axial}\hat{d} + k_t g_t$, tangent
$K=B^\mathsf{T}DB$. The constraint gap is $O(\sigma/k)$ — exact only as
$k\to\infty$. Perfect-bond conditioning at very large $k$ is the price.

### 5.2 `al` (augmented Lagrangian)

Adds a per-element multiplier $\lambda\in\mathbb{R}^{ndm}$ (committed state):

$$t = F_\text{axial}\hat{d} + k_t g_t + \lambda ,$$

with the **tangent unchanged** ($\lambda$ is constant within an inner solve). A
per-step **Uzawa update** fires once per converged step in `commitState`:

$$\lambda \leftarrow \lambda + \Delta\lambda,\quad
  \Delta\lambda = k_t g_t\ (\text{transverse, always}) + k_\text{axial}^\text{perfect} s\,\hat{d}\ (\text{axial, perfect-bond only}).$$

The converged gap contracts geometrically, so $g\to 0$ at **moderate** $k_t$
(no $k\to\infty$). This is the headline win: AL reaches near-exact constraint
satisfaction at $k_t\approx k_\text{host}$, dodging the ill-conditioning of a
stiff perfect-bond penalty and *enlarging* the explicit $\Delta t_\text{cr}$. With
a bond law, AL augments **only** the transverse constraint — the physical τ–s
axial force is left untouched. No analysis-core change (`commitState` is
driver-called once per converged step). Under `-corot`, the transverse multiplier
is re-projected each step so a rotating axis can't leak it into a spurious axial
force.

### 5.3 Not built / deferred

`-enforce nitsche` (research-grade — needs a host-stress/∂N query that re-couples
to host types) and `-enforce transformation` (deferred indefinitely, §1.4) are
**rejected at parse time** with a guiding message.

---

## 6. Large rotation — co-rotated bar axis (`-corot`)

The tie ($g$, $N_i$) is already objective; the only large-rotation defect is the
**constitutive split** against a *frozen* $\hat{d}$. With $\hat{d}$ fixed, a pure
rigid host rotation $Q$ registers spurious axial slip ($s=(Qg)\cdot\hat{d}_0 \ne g\cdot\hat{d}_0$)
and a non-objective traction.

`-corot` cures it by recomputing the bar axis each step as the **secant** between
the embed point (weights $N_i$) and a second point **B** along the bar
(weights $N_i^B$ from `-xiB`/`-shapeB`), using **current** host node positions:

$$\hat{d}(t) = \operatorname{normalize}\!\Big(\textstyle\sum_i N_i^B x_i - \sum_i N_i x_i\Big),
\quad x_i = X_i + u_i .$$

This captures host rigid rotation exactly, restoring $t\to Qt$. Default **off**
(frozen `-dir`, bit-identical). v1 omits the $\partial\hat{d}/\partial u$ geometric
tangent term (standard EICR practice — exact for explicit, converges for moderate
per-step rotation). It **lifts the small-rotation restriction** for Mode P; the
small-*stretch* $L_\text{trib}\cdot\lambda$ correction stays deferred.

---

## 7. Explicit stability — bipenalty (`-bipenalty`)

A stiff penalty tie on a near-massless rebar node injects a spurious
high-frequency mode that collapses the explicit critical step
$\Delta t_\text{cr}=2/\omega_\text{max}$. **Bipenalty** (Hetherington & Askes 2009)
pairs the stiffness penalty $k_\text{eff}=\max(k_\text{axial},k_t)$ with a **mass
penalty** $m_p$, lumped on the rebar (slave) node only, so the spurious frequency
$\omega_p=\sqrt{k_\text{eff}/m_p}$ is bounded **independently of stiffness** — and
the global explicit step is stabilized at a value you choose.

| mode | mass penalty | reported $\Delta t_\text{cr}$ |
|------|--------------|-------------------------------|
| **`-dtcr` $\Delta t^\*$** | $m_p = k_\text{eff}(\Delta t^\*/2)^2$ | $\Delta t^\*$ exactly (no host query) |
| **`-wcap` $\beta$** | $m_p = k_\text{eff}/(\beta\,\omega_\text{host})^2$ | $2/(\beta\,\omega_\text{host})$ |

with $\omega_\text{host}\approx\sqrt{\lVert K_\text{host}\rVert_\infty/\lVert M_\text{host}\rVert_\infty}$
and the self-report $\Delta t_\text{cr}^{\,\text{elem}}=2\sqrt{m_p/k_\text{eff}}$.

**Lumped-on-slave is a hard constraint**: the faithful $m_p B^\mathsf{T}B$ has host
off-diagonal blocks that would densify M and corrupt a `DiagonalSOE`; lumping the
full $m_p$ on the rebar node keeps M diagonal (DiagonalSOE-safe) and preserves the
frequency bound.

> [!warning] Honest limits (verified by adversarial review)
> - **The per-element `CriticalTimeStep` eigensolve cannot see the tie.** For the
>   coupling element's own pencil $K v=\lambda\,\operatorname{diag}(m_p I,0)\,v$ the
>   massless free host DOFs slave out the constraint (Schur complement $=0$ ⇒
>   $\lambda_\text{max}=0$). The bound is delivered by **global** stabilization +
>   the **self-report**, and carried into the scan by the §10.6.1 seam below.
> - **CDL does not auto-replace `dt`** — `-cfl` *reports*, `-cflAbort` *guards*; you
>   still pass `dt`. Read `eleResponse "dtcr"` to choose it.
> - **`-wcap` reports a heavy-host approximation** (assumes host heavier than
>   $m_p$; true for concrete/steel). A lighter host node's own per-element
>   $\Delta t_\text{cr}$ enters the same global min independently.
> - **Explicit only** — $m_p$ is real added mass; do not enable in implicit
>   transient. **Gated on `-enforce penalty`** (auto-disabled under `al`).

**The `-cfl` seam (§10.6.1).** A vanilla `virtual double Element::getExplicitCriticalTimeStep()`
(default −1 = "no opinion"; [[LEDGER_vanilla_files]]) is folded into
`CriticalTimeStep`'s per-element loop — a non-negative return replaces the
eigensolve for that element. `LadrunoEmbeddedRebar` returns $2\sqrt{m_p/k_\text{eff}}$,
so `ops.criticalTimeStep()` / `-cflAbort` now honor the embedded bound.

See the deep treatment (theory + derivation of the $\lambda_\text{max}=0$
degeneracy) in [[20_ladruno_embedded_reinforcement_adr|ADR 20]] §10.6 / §10.6.1.

---

## 8. Diagnostics & responses

`eleResponse $tag <name>`:

| name | returns |
|------|---------|
| `slip` | axial slip $s$ |
| `bondStress` | τ (bond) or $k_\text{axial}s$ (perfect) |
| `axialForce` | $F_\text{axial}$ |
| `force` / `localForce` | traction vector $t$ |
| `kt` / `penalty` | resolved transverse penalty (useful with `-kt auto`) |
| `penaltyEnergy` | **artificial** energy $\tfrac12 k_t\lVert g_t\rVert^2$ (+ $\tfrac12 k_\text{axial}s^2$ perfect-bond) |
| `constraintViolation` | transverse gap $\lVert g_t\rVert$ |
| `bondEnergy` / `bondDissipation` | **physical** bond work $\propto\!\int\!\tau\,\mathrm{d}s$ (0 for perfect bond) |
| `dir` / `barAxis` | current working bar axis (co-rotated under `-corot`) |
| `augLambda` / `lambda` | AL multiplier |
| `gap` | raw constraint gap $g$ |
| `mpenalty` | bipenalty mass penalty $m_p$ |
| `dtcr` | self-reported critical step $2\sqrt{m_p/k_\text{eff}}$ |

> [!tip] Netting artificial energy
> For an energy balance, bin embedded elements in their own `MeshRegion` and
> subtract $\sum$ `penaltyEnergy` to recover the net physical balance; `bondEnergy`
> is the genuine τ–s dissipation.

---

## 9. Implementation map

| Concern | Location |
|---------|----------|
| Element + parser | `SRC/element/ladrunoEmbeddedRebar/{LadrunoEmbeddedRebar.{cpp,h},OPS_LadrunoEmbeddedRebar.cpp}` |
| Axial law | `SRC/material/uniaxial/LadrunoBondSlip.{cpp,h}` (MAT 33002) |
| Tie / traction / tangent | `formBandTraction`, `getResistingForce`, `getTangentStiff` |
| AL Uzawa | `commitState` (per converged step) |
| Co-rotated axis | `currentBarAxis` (`-corot` secant) |
| Auto-kt | `resolveAutoKt` (lazy, from host `getInitialStiff`) |
| Bipenalty $m_p$ | `resolveBipenalty`, `getMass` (lumped), `getResistingForceIncInertia` (+$m_p a$), `setRayleighDampingFactors` (no-op) |
| Self-reported $\Delta t_\text{cr}$ | `getExplicitCriticalTimeStep` → folded in `CriticalTimeStep::computeCriticalTimeStep` |
| Host shape query | `Element::getInterpolationWeights` (vanilla virtual; `LadrunoBrick`/`BezierTet10` overrides) |
| Serialization | 18-field header; $m_p$/$k_t$ re-resolved after `recvSelf`; broker case 33005/33002 |

**Vanilla (upstream) edits** (logged in [[LEDGER_vanilla_files]]): `Element.{h,cpp}`
(two added virtuals — `getInterpolationWeights`, `getExplicitCriticalTimeStep`),
`CriticalTimeStep.cpp` (the self-report fold-in), plus the command/broker
registrations. New fork code is in [[LEDGER_implementations]].

---

## 10. Validation

Zone-A batteries (no host continuum needed for the mechanics — prescribed host
displacements play the host role):

- **Coupling mechanics** — perfect axial spring $k\cdot u$; anisotropic transverse
  $k_t\cdot u$; bond-slip axial follows $\tau(s)$; partition-of-unity force split
  to hosts + equilibrium; serialization round-trip.
- **`-host`/`-xi`** — real `LadrunoBrick` host: node auto-resolve, centroid
  weights all $0.125$, off-centroid weights match the trilinear formula.
- **`-kt auto`** — $k_t$ scales linearly with `-ktAlpha`; drives the transverse
  force; `penaltyEnergy`/`constraintViolation`/`bondEnergy` channels.
- **`-corot`** — axis follows $Q\hat{d}_0$ under rigid host rotation (frozen
  without); the force path uses the rotated axis; `-shapeB` ≡ `-xiB`.
- **`-enforce al`** — Uzawa drives the perfect-bond gap → 0 at moderate $k_t$
  (penalty leaves $P/k_t$); the multiplier carries the load; bond-slip axial
  untouched; `penalty` bit-identical to default.
- **`-bipenalty`** — off bit-identical; `-dtcr` formula + `dtcr=dt`; `-wcap`
  inverse-β scaling; host/massless guards; AL disables it; explicit SDOF stable at
  $0.9\times$ / unstable at $1.1\times$ the target; `-cfl` governance (the embedded
  bound governs `ops.criticalTimeStep`).

File: `tests/test_ladrunoEmbeddedRebar_element.py` (+ `tests/test_ladrunoBondSlip_material.py`).

---

## 11. Use cases & recipes

### 11.1 Building a tie (host query)

```python
# rebar node 1 embedded in LadrunoBrick host 100 at its centroid,
# perfect bond axially, auto-scaled transverse hold
ops.element("LadrunoEmbeddedRebar", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
            "-dir", 1.0, 0.0, 0.0, "-perfect", 1.0e8,
            "-kt", "auto", "-ktAlpha", 1.0e3, "-enforce", "penalty")
```

The apeGmsh generator (planned `g.reinforce(host, bars, bond, mode)`) runs the
guarded inverse map (global bar point → ξ) at build time and emits the `-host -xi`
primitives, so users never hand-compute weights.

### 11.2 Choosing the axial law

| Goal | Use |
|------|-----|
| Stiff bond, no slip modeling | `-perfect <large k>` (+ `-enforce al` for near-exact at moderate k) |
| Realistic bond-slip / pull-out | `-bond <LadrunoBondSlip> -bondScale <πd·L_trib>` + disp/IMPLEX control |

### 11.3 Implicit static pushover

Use `-enforce al` for near-exact perfect bond at well-conditioned stiffness; omit
`-bipenalty` (no mass term in statics). Both pre-build gates needed
KrylovNewton + step-halving — wrap pushovers in the adaptive-stepping proc.

### 11.4 Explicit RC dynamics

1. Build ties with `-enforce penalty -bipenalty -dtcr <dt>` (or `-wcap β`).
2. `integrator CentralDifferenceLadruno -cflAbort`; `analyze(1, dt)` to compute
   the limit; read `ops.criticalTimeStep()` (now accounts for the ties).
3. Run with `dt ≤ ops.criticalTimeStep()`; `-cflAbort` aborts rather than going
   unstable. Audit per-tie via `eleResponse "dtcr"` / `"mpenalty"`.

### 11.5 Pitfalls

- `-xi` / `-kt auto` / `-wcap` / `-xiB` all require the **`-host`** form.
- Softening bond-slip under load control diverges — use disp/arc/IMPLEX.
- `-bipenalty` is **explicit-only**; auto-disabled under `-enforce al`.
- Corot needs a point **B** along the bar (`-xiB`/`-shapeB`).
- Quantitative **dowel** action is mesh-limited — the default `truss` is
  axial-only; `beam` gives only an approximate transverse contribution (v2 dowel
  interface deferred).

---

## 12. References & related

- **Design record:** [[20_ladruno_embedded_reinforcement_adr|ADR 20 — Embedded Reinforcement + Bond-Slip for 3D RC]] (D1–D6, §10 roadmap, §10.6 bipenalty).
- **RC-3D context:** [[19_rc3d_modeling_recipe]], [[22_rc3d_conformal_recipe]], [[21_rc3d_validation_gates]].
- **Reference codes:** Abaqus `*EMBEDDED ELEMENT` (perfect-bond translations-only),
  LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID` (penalty + bond + mass scaling), DIANA
  bond-slip reinforcement (own-DOF rebar + CEB-FIP τ–s).
- **Bond-slip:** CEB-FIP Model Code 2010; ADR D4.
- **Augmented Lagrangian:** Simo & Laursen (1992); Wriggers, *Computational Contact
  Mechanics* §6.4; Bertsekas (1982).
- **Bipenalty:** Hetherington & Askes (2009) *Comput. Struct.* 87:1474; Askes &
  Hetherington (2011); Caleyron et al. (2013). Explicit stability: Belytschko et
  al. (2014) §6.7; Hughes (1987) Ch. 9.
- **Companions:** [[ladruno_integrators_guide]], [[central_difference_ladruno_guide]],
  [[LadrunoUniaxialJ2_guide]], [[LadrunoRebarBuckling_guide]].
