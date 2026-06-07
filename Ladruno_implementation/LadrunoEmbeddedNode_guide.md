---
title: "LadrunoEmbeddedNode — General Node-to-Host Coupling (Embedment) Element"
project: Ladruno
type: reference / user guide
status: v1 SUPPORTED CORE shipped (U tie + g0 stress-free birth + penalty/AL/bipenalty); UR/UP/D9/-corot EXPERIMENTAL
element: element LadrunoEmbeddedNode
classTag: 33006 (element)
related:
  - "[[23_ladruno_embedded_node_adr]]"
  - "[[27_ladruno_embedded_node_validation_plan]]"
  - "[[LadrunoEmbeddedRebar_guide]]"
  - "[[LadrunoBondSlip_guide]]"
  - "[[ladruno_apegmsh_contract]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - guide
  - element
  - embedded-node
  - mesh-tie
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
  - staged-construction
  - ssi
updated: 2026-06-07
---

# LadrunoEmbeddedNode — general node-to-host coupling element

`element LadrunoEmbeddedNode` (classTag **33006**) ties **one constrained node** to
the nodes of a **host element** through precomputed shape-function weights `N_i(ξ)`,
so an arbitrary node embeds in a **non-matching** host mesh. It is the **isotropic
sibling** of [[LadrunoEmbeddedRebar_guide|LadrunoEmbeddedRebar]] over the same shared
coupling kernel: same penalty / augmented-Lagrangian / `-kt auto` / bipenalty
machinery, but with an **isotropic** translational tie (no bar axis, no axial/
transverse split, no bond-slip). It is the fork's general **embedment** primitive —
mesh stitching, point embedding, beam/part-into-solid, soil-structure interaction.

This is the single user/developer reference. It covers the **theory** (the tie, the
gap, stress-free birth, the penalty formulation), the **capabilities** (the option
grammar and what each switch does), the **implementation** (where the math lives,
serialization, responses), the **validation** (the battery that earns the claim), and
the **use cases** — including the **contract for the apeGmsh team** to emit ties
automatically (§10.6). The design record is [[23_ladruno_embedded_node_adr|ADR 23]].

> [!important] v1 scope — supported core vs experimental (read this)
> Only the **U translational tie + `g0` stress-free birth + penalty/AL/bipenalty
> enforcement** is the **validated, world-class core** (the part to rely on and to have
> apeGmsh emit by default). **UR (rotation), UP (pressure), D9 (material interface) and
> `-corot` are EXPERIMENTAL** — they work as coded but are **not validated**; treat them
> as research switches, not production features. See [[23_ladruno_embedded_node_adr|ADR §14]]
> for the rationale and [[27_ladruno_embedded_node_validation_plan|the validation plan]]
> for the battery. This guide flags every experimental switch with **⚗️**.

> [!info] Where it fits
> Use `LadrunoEmbeddedNode` when meshes **don't conform** — an independently-meshed
> sub-block in a coarse host, a beam/shell node embedded in a solid, a pile/foundation
> node in a soil continuum, a control/recorder/mass node tied into a solid interior.
> For **reinforcement** (a bar with bond-slip) use [[LadrunoEmbeddedRebar_guide|LadrunoEmbeddedRebar]];
> for **rigorous frictional contact** use the scoped `LadrunoContact`. It exceeds the
> one OpenSees precedent (`ASDEmbeddedNodeElement`, implicit-only, tri/tet-only, raw
> `1e18` penalty) on **host coverage** (hex/tet/quad), **conditioning** (`-kt auto`/AL),
> and **explicit dynamics** (bipenalty + self-reported `dt_cr`).

---

## 1. Theory — the coupling

### 1.1 The tie and the gap (U)

The element ties the constrained ("slave") node `c` to a point **interpolated inside a
host element** at natural coordinate `ξ`. External nodes are `[c, host_1, …, host_M]`.
The translational **gap** is the current relative displacement of the constrained point:

$$\mathbf{g} \;=\; \mathbf{u}_c \;-\; \sum_{i=1}^{M} N_i(\xi)\,\mathbf{u}_{\text{host},i}
\qquad (\text{the first } n_{dm} \text{ DOFs of each node}).$$

The weights `N_i(ξ)` are the host's shape functions evaluated at the embedded point;
they form a **partition of unity** (`Σ N_i = 1`), so a rigid host translation gives
`g`-invariance and an affine host field is reproduced exactly. The tie is **isotropic**
and therefore **already frame-objective** under rigid host rotation — unlike the
anisotropic rebar, it needs **no co-rotation** for the U core.

### 1.2 Stress-free birth — the `g0` offset capture (the headline)

A constraint that enforces the **absolute** tie `u_c = Σ N_i u_host` would, when added
to an **already-deformed** host (staged construction), activate with `g = −N·u_host ≠ 0`
and **yank** the slave by the full accumulated host displacement. To avoid that, at
`setDomain` the element captures the gap **once** (`g0`, the parent
`ASDEmbeddedNodeElement` `m_U0` pattern) and thereafter drives **all** traction from the
**relative** gap:

$$\mathbf{t} \;=\; K_u\,(\mathbf{g} - \mathbf{g}_0),\qquad
  \mathbf{g}_0 = \mathbf{g}\big|_{\text{activation}}.$$

The element is **born force- and stress-free**: at activation `(g − g0) = 0`, so the
coupling force, the penalty energy, and (in the experimental D9 mode) the interface
material's strain are all zero. It then resists only **subsequent** relative motion. The
capture is a **no-op when the element is added at the undeformed state** (`g0 = 0`), so
non-staged models are byte-identical to the absolute tie. Default **on**; `-absolute`
opts out (a deliberate snap-to-host, or legacy behaviour).

### 1.3 Isotropic traction, resisting force and tangent

The constitutive law is the isotropic degenerate of the rebar's anisotropic `D`:

$$\mathbf{t} = K_u\,(\mathbf{g}-\mathbf{g}_0),\qquad \mathbf{D}_u = K_u\,\mathbf{I}_{n_{dm}}.$$

With the discrete operator `B = [ I | −N_1 I … −N_M I ]`, the resisting force and tangent are

$$\mathbf{P} = \mathbf{B}^{\mathsf T}\mathbf{t},\qquad \mathbf{K} = \mathbf{B}^{\mathsf T}\mathbf{D}_u\,\mathbf{B}.$$

The kernel assembles these in a compact `(1+M)·n_dm`-packed space, then the element
**scatters** each node's `n_dm`-block into the full per-node-`ndf` layout (`nDOF = Σ ndf_i`),
so the constrained / host nodes may carry **more** DOFs than `n_dm` — only the first
`n_dm` translational DOFs are coupled by the U core. Because `D_u = K_u·I` is
state-independent, `getInitialStiff ≡ getTangentStiff` is exact for the core.

### 1.4 Penalty formulation — implicit *and* explicit

The tie is a **penalty** coupling: it adds **stiffness only** (`K = BᵀD_uB`), never a
constraint equation that eliminates DOFs. Two consequences:

- It composes with **any** solver/integrator — implicit static/transient **and**
  explicit (`CentralDifferenceLadruno`, `ExplicitBathe`).
- The constraint is satisfied to `≈ P/K_u` (a small residual gap). Tighten it with
  `-k auto` conditioning, or make it near-exact at **moderate** `K_u` with augmented
  Lagrangian (§4). Raw huge penalties (the ASD `1e18` anti-pattern) wreck conditioning
  and the explicit `dt_cr` — **do not** port them.

### 1.5 Host-agnostic interpolation

The weights `N_i(ξ)` come from the host, not the element:

- **`-xi ξ…`** — query the host via the vanilla `Element::getInterpolationWeights(ξ, N)`
  (overridden on the 3D fork solids `LadrunoBrick` trilinear / `BezierTet10` Bernstein).
  Needs the `-host` form. **3D hosts only** today.
- **`-shape N_1…N_M`** — supply the weights explicitly (any host, any dimension). This is
  how **2D** hosts and any non-overriding host are served: the generator (apeGmsh)
  computes the weights and passes them. The **inverse map** (global point → ξ) always
  lives in the generator, never the element.

### 1.6 Rotations are free by default (Abaqus-consistent)

By default the element couples **only translations**; the constrained node's rotation
(or pressure) DOFs are **left free**. This matches Abaqus `*EMBEDDED REGION` (embedded
rotations free — the embedded beam/shell carries its own bending) and is the correct
default for a **pinned** embedment. Moment transfer through a *single* point is a pin and
cannot carry moment — embed a **length** (≥2 ties with a lever arm) so the moment is
carried as a force couple through the translational ties (§10.5). The experimental `-rot`
tie couples rotations to the host **continuum spin** `½ curl(u)`, which is **not** a
validated moment-transfer mechanism (§6).

---

## 2. Capabilities — the option grammar

```
element LadrunoEmbeddedNode tag cNode {nHost h1..hN | -host eleTag}
        {-shape N1..NN | -xi x1..x_ndm}            # host weights (supply or query)
        [-k {Ku | auto}] [-kAlpha a]               # isotropic penalty (numeric or host-scaled)
        [-enforce {penalty | al}]                  # constraint enforcement
        [-bipenalty {-dtcr dt | -wcap beta}]       # explicit critical-step control
        [-absolute]                                # opt out of g0 stress-free birth
        # ---- EXPERIMENTAL (⚗️ not validated) ----
        [-pressure [-kp Kp]]                       # ⚗️ UP: tie the u-p pressure DOF
        [-rot [-kr {Kr|auto}] [-krAlpha a]         # ⚗️ UR: tie rotations to host continuum spin
              {-dNdx N1x N1y [N1z].. | -gradXi ξ..}]
        [-normal nx ny [nz] [-orient ox oy oz]     # ⚗️ D9: per-direction interface materials
         -matN tag | -matT1 tag | -matT2 tag ..]   #       (cohesive / gap / bond / approx-friction)
        [-corot]                                   # ⚗️ D9: co-rotate the material frame
```

| Group | Switch | Meaning |
|-------|--------|---------|
| **Host** | `nHost h…` / `-host eleTag` | host nodes explicit, or pulled off a host element |
| **Weights** | `-shape N…` / `-xi ξ…` | supply weights, or query the host (3D, needs `-host`) |
| **Penalty** | `-k Ku` / `-k auto` `-kAlpha a` | numeric, or `K_u = a·max\|K_host(i,i)\|` (needs `-host`) |
| **Enforce** | `-enforce penalty\|al` | penalty (default) or augmented Lagrangian |
| **Explicit** | `-bipenalty -dtcr dt` / `-wcap β` | mass penalty for an explicit `dt_cr` (needs `-host` for `-wcap`) |
| **Staged** | `-absolute` | keep the absolute tie (no `g0` capture) |
| ⚗️ **UP** | `-pressure -kp Kp` | also tie the pressure DOF (all nodes u-p) |
| ⚗️ **UR** | `-rot -kr {Kr\|auto} -krAlpha a` + gradients | also tie rotations to `½curl(u)` |
| ⚗️ **D9** | `-normal`/`-orient` + `-matN/-matT*` + `-corot` | material-driven interface |

`-rot` and `-pressure` are **mutually exclusive** (the extra DOF is rotation *or*
pressure — ambiguous in 2D `ndf=3`); passing both is a parse error.

---

## 3. Penalty & auto-scaling (`-k`, `-k auto`)

`-k Ku` sets a numeric penalty (default `1e12`). `-k auto -kAlpha a` resolves it from the
host's own initial stiffness:

$$K_u = a \cdot \max_i |K_{\text{host}}(i,i)|\qquad (\text{via } \texttt{host->getInitialStiff()}),$$

the max-absolute diagonal (which already carries `~E·lch` units). Default `a = 1e3`. This
makes the penalty **mesh- and material-independent** in conditioning — a coarse stiff host
and a fine soft host both get a well-scaled tie. Requires the `-host` form. **Do not** pass
ASD's `1e18` as `a` (it is not `E`-scaled — condition-number blow-up). A non-numeric `-k`
token (a forgotten value, e.g. `-k -enforce`) is **rejected at parse time**, not silently
turned into `K_u=0`.

---

## 4. Enforcement strategies (`-enforce`)

### 4.1 `penalty` (default)
Adds `K = BᵀD_uB`; the constraint holds to a residual gap `≈ P/K_u`. Implicit + explicit.

### 4.2 `al` (augmented Lagrangian)
Adds a per-element multiplier `λ` (size `n_dm`) carried with the **same** tangent, plus a
per-step **Uzawa** update in `commitState`: `λ ← λ + K_u·(g − g0)`. This drives the gap to
**near zero at MODERATE `K_u`** (no `K→∞`), fixing both conditioning and the explicit
`dt_cr` blow-up. The isotropic tie has **no preferred axis**, so — unlike the rebar — there
is **no directional re-projection** of `λ`. `nitsche` / `transformation` are rejected at
parse time. `-bipenalty` is auto-disabled under `-enforce al` (AL needs no mass penalty).

---

## 5. Explicit stability — bipenalty (`-bipenalty`)

A pure penalty coupling with zero mass is an explicit-`dt_cr` false-safe: the per-element
eigensolve sees the massless coupling as `λ_max = 0` and reports no bound, yet the global
penalty mode is stiff. **Bipenalty** lumps a small **mass penalty** `m_p` on the
constrained node so the coupling mode has a finite, bounded frequency, and the element
**self-reports** its critical step.

Two budgets:

| Budget | `m_p` | `dt_cr` |
|--------|-------|---------|
| `-dtcr dt` | `m_p = k_eff·(dt/2)²` | the requested `dt` (slave-side) |
| `-wcap β` | `m_p = k_eff/(β·ω_host)²`, `ω_host² = ‖K_host‖/‖M_host‖` | `2/(β·ω_host)` (host-aware; needs `-host`) |

**Host-side reduced mass.** The mass penalty is lumped on the **slave**, but the coupling
stiffness `N_i²·k_eff` also loads the (penalty-massless) host DOFs, so the true
relative-mode step uses the **reduced mass** `μ = m_p·M_h/(m_p+M_h)`:

$$\Delta t_{\text{cr}} = 2\sqrt{\mu/k_\text{eff}},\qquad \mu = \frac{m_p\,M_h}{m_p+M_h}.$$

When the host element is queryable, the report **tightens with `μ`** (`M_h` = smallest
positive host mass diagonal); a massless host emits a warning. Reporting `2√(m_p/k_eff)`
(i.e. `M_h→∞`) is **unconservative for a light/coarse host** — the `-dtcr`-without-host form
is therefore a *user-asserts-`dt`* contract. Prefer `-wcap` (host-aware) for an
automatically host-bounded step.

The self-reported `dt_cr` is the **min over active DOF classes** and is folded into
`ops.criticalTimeStep` via the vanilla `Element::getExplicitCriticalTimeStep` →
`CriticalTimeStep` seam, so `-cflAbort` honours the embedded bound. Off ⇒ `m≡0` ⇒
bit-identical. The element **refuses Rayleigh factors** (a pure penalty coupling carries no
physical damping, so a `βK` can't spuriously shrink `dt_cr`).

---

## 6. ⚗️ Experimental modes (not validated)

These are flag-gated and **off by default**; the core never sets them. Documented for
completeness — do not rely on them in production.

- **UP — pressure tie (`-pressure -kp Kp`).** Ties the constrained node's pressure DOF
  (index `n_dm`, the u-p convention) to the host's interpolated pressure,
  `g_p = p_c − Σ N_i p_{\text{host},i}`, `t_p = K_p g_p`. Activated only if all coupled nodes
  are u-p (`ndf ≥ n_dm+1`), else degrades to U-only. Niche poromechanics / SSI.
- **UR — rotation tie (`-rot -kr {Kr|auto}` + gradients).** Ties the constrained node's
  rotations to the host **continuum spin** `θ = ½ curl(u) = skew(∇u)` at the embedded point
  (3D: 3 rotations; 2D: the drilling `R_z`), read from the host cartesian gradients `∂N/∂x`
  (`-dNdx` explicit, `-gradXi`/`-xi` queried). **This is SPIN coupling, not lever-based
  moment transfer:** on a CST/TET4 host `∂N/∂x` is element-constant so it degenerates to a
  rigid element-wide spin. Use it as a rotational restraint, **not** a validated
  shell-to-solid moment connection (for that, embed a length — §10.5).
- **D9 — material-driven interface (`-normal` + `-matN/-matT*` + `-corot`).** Each local
  direction (normal `e_0`, tangents) carries a **uniaxial material** instead of the bare
  penalty: cohesive separation, unilateral contact/gap (`ENT`/`ElasticPPGap`), elastic
  bedding, bond. Coulomb friction is only **approximate** (uncoupled per-direction — a fixed
  `ElasticPP` slip, not `μ·N`; rigorous → `LadrunoContact`). `-corot` co-rotates the frame
  with the host spin. **Known must-fix before promotion:** `getInitialStiff` aliases the
  state-dependent D9 tangent (ADR §14.4(b) / §15 #1).

---

## 7. Diagnostics & responses

`eleResponse(tag, …)`:

| Response | Returns |
|----------|---------|
| `force` | coupling traction `t` (size `n_dm`) |
| `gap` | the **relative** gap `g − g0` (size `n_dm`) |
| `constraintViolation` | `‖g − g0‖` (scalar) |
| `k` / `penalty` | resolved `K_u` |
| `penaltyEnergy` | `½ K_u ‖g−g0‖²` |
| `augLambda` / `lambda` | AL multiplier (size `n_dm`) |
| `mpenalty` | resolved mass penalty `m_p` |
| `dtcr` | self-reported critical step (host-reduced when available) |
| `initGap` / `offset` | the captured activation offset `g0` (zero if never captured / `-absolute`) |
| ⚗️ `pgap`/`pforce`/`plambda` | UP diagnostics |
| ⚗️ `rgap`/`rforce`/`rlambda`/`kr` | UR diagnostics |
| ⚗️ `localGap`/`localForce`/`normal` | D9 local-frame diagnostics |

The `gap`/`initGap` pair is the staged-construction audit: after activation, `gap ≈ 0` and
`initGap = −Σ N_i u_host` confirms a stress-free birth.

---

## 8. Implementation map

| Concern | Location |
|---------|----------|
| Element + parser | `SRC/element/ladrunoEmbeddedNode/{LadrunoEmbeddedNode.{cpp,h},OPS_LadrunoEmbeddedNode.cpp}` |
| Shared coupling kernel | `SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.{cpp,h}` (gap, `B`-assembly, auto-`k`, bipenalty, `dt_cr`) |
| Gap (relative) | `computeGap` / `computeGapP` / `computeGapR` (subtract `g0`/`gp0`/`gr0`) |
| Stress-free birth | `captureInitialGap` (once, in `setDomain`) |
| Traction / tangent | `formTransTraction`, `getResistingForce`, `getTangentStiff` (`getInitialStiff ≡` for the core) |
| AL Uzawa | `commitState` (per converged step) |
| Auto-`k` | `resolveAutoKt` (lazy, from host `getInitialStiff`) |
| Bipenalty `m_p` / `dt_cr` | `resolveBipenalty`, `getMass` (lumped), `getResistingForceIncInertia` (+`m_p a`), `getExplicitCriticalTimeStep` (host-reduced `μ`) |
| Host weight query | `Element::getInterpolationWeights` (vanilla; `LadrunoBrick`/`BezierTet10` overrides) |
| ⚗️ UR gradients | `Element::getInterpolationGradients` (vanilla; same hosts) |
| Serialization | 29-field header + payload; `g0`/`gp0`/`gr0` + `g0Computed` round-tripped; broker case 33006 |

**Vanilla (upstream) edits** ([[LEDGER_vanilla_files]]): `Element.{h,cpp}`
(`getInterpolationGradients` for the UR slice; `getInterpolationWeights` /
`getExplicitCriticalTimeStep` shared with the rebar), `CriticalTimeStep.cpp` (the
self-report fold-in), command/broker registrations. Fork code in [[LEDGER_implementations]].

---

## 9. Validation

The **supported-core** battery (Zone-A; prescribed host displacements play the host role
where no continuum is needed) — full plan in
[[27_ladruno_embedded_node_validation_plan|the validation plan]]:

- **U mechanics** — isotropic spring `K·g`; partition-of-unity force split to hosts +
  equilibrium; objectivity under rigid host rotation; serialization round-trip.
- **`-host`/`-xi`** — real `LadrunoBrick`/`BezierTet10`: node auto-resolve, centroid weights
  `0.125`, off-centroid weights match the trilinear/Bernstein formula.
- **`-k auto`** — `K_u` scales linearly with `α`; bounded condition number.
- **Stress-free staged birth** — settle the host, then embed onto the deformed host: zero
  force + zero relative gap at activation; offset `= −N·u_host`; slave not yanked; the
  `-absolute` contrast shows the jolt; undeformed add is byte-identical.
- **`-enforce al`** — Uzawa drives the gap → 0 at moderate `K`; `penalty` bit-identical.
- **`-bipenalty`** — off bit-identical; `-dtcr` formula + `dtcr=dt`; **host reduced-mass
  `μ`** tightens the report on a queryable host; `-wcap` host-aware; explicit SDOF stable at
  `0.9×` / unstable at `1.1×`.

File: `tests/test_ladrunoEmbeddedNode_element.py`. The element passed a full adversarial
review (ADR §15): the supported core is sound (zero blockers, zero core mechanics bugs).

---

## 10. Use cases & recipes

### 10.1 Embed a node into a solid host (host query)

```python
# constrained node 1 embedded in LadrunoBrick host 100 at its centroid,
# host-stiffness-scaled penalty, born stress-free (g0 capture on by default)
ops.element("LadrunoEmbeddedNode", 1, 1, "-host", 100, "-xi", 0.0, 0.0, 0.0,
            "-k", "auto", "-kAlpha", 1.0e3, "-enforce", "penalty")
```

### 10.2 Non-matching solid↔solid mesh tie

Tie each node of a refined sub-block to the coarse host element it falls inside. The
generator (apeGmsh, §10.6) locates the host and inverse-maps each node to `ξ`, emitting
`-host -xi` per tie. Use `-k auto` for conditioning; `-enforce al` if you need the gap
near-exact under implicit pushover.

### 10.3 Staged construction (the `g0` workflow)

```python
# stage 1: build + settle the host, commit (host is now deformed)
# ... analyze the host under gravity, then:
ops.node(500, *xyz)                                   # the new (fresh) node, u = 0
ops.element("LadrunoEmbeddedNode", 7, 500, "-host", 42, "-xi", *xi, "-k", "auto")
assert abs(ops.eleResponse(7, "force")[0]) < 1e-6     # born stress-free — no jolt
assert ops.eleResponse(7, "gap")[0] == 0.0            # relative gap zero at activation
# subsequent stage-2 loading is resisted normally
```

Pass `-absolute` only to deliberately **snap** a misplaced node onto the host point.

### 10.4 Explicit dynamics (SSI / impact)

1. Build ties with `-enforce penalty -bipenalty -wcap β` (host-aware) or `-dtcr <dt>`.
2. `integrator CentralDifferenceLadruno -cflAbort`; read `ops.criticalTimeStep()` (now
   accounts for the ties, host-reduced).
3. Run with `dt ≤ ops.criticalTimeStep()`. Audit per-tie via `eleResponse "dtcr"` /
   `"mpenalty"`. Pairs naturally with the absorbing-boundary / PML stack.

### 10.5 Frame → 2D/3D solid (DOF-mismatch interface)

A frame node (`ndf=3` in 2D: `u_x,u_y,R_z`) into a plane solid (`ndf=2`): the U tie couples
the translations and **leaves `R_z` free** → a clean **pin** (Abaqus-consistent). For a
**moment** connection, **do not** use the single-point `-rot` spin tie — instead **embed a
length** of the frame (≥2 nodes spanning a lever arm `d`), each with a U tie: a moment `M`
at the frame end is then carried as a **force couple** `F = M/d` through the translational
ties + the frame's own bending over the embedded stub (the embedded-beam-in-solid
mechanism). The plane mesh must be fine enough to **resolve** the couple (the lever arm
should span ≥1 element).

### 10.6 apeGmsh integration — generator contract (TO IMPLEMENT)

apeGmsh already emits `ASDEmbeddedNodeElement` for non-matching ties; `LadrunoEmbeddedNode`
is the **drop-in fork upgrade** (hex/tet/quad hosts + explicit + AL + stress-free birth).
**apeGmsh owns the inverse map** (global point → host `ξ`); the element only consumes
weights. Ship a generator, e.g.:

```python
g.embed(nodes=<set>, host=<set>,            # constrained nodes + the host region
        k="auto", kAlpha=1e3,               # conditioning
        enforce="penalty",                  # or "al"
        explicit=None,                      # or {"wcap": 2.0} / {"dtcr": 1e-4}
        staged=True)                        # g0 stress-free birth (default)
```

For each constrained node the generator:

1. **locates** the host element containing the node and **inverse-maps** to `ξ` (the
   guarded straight-sided inverse map — same as the rebar `g.reinforce`);
2. emits, per node,
   `element LadrunoEmbeddedNode <tag> <cNode> -host <hostEle> -xi ξ η ζ -k auto -kAlpha 1e3 [-enforce …] [-bipenalty …]`;
3. **3D hosts** (`LadrunoBrick`/`BezierTet10`) own the weights via `getInterpolationWeights`
   → use `-xi`. **2D hosts** have no `getInterpolationWeights` override yet → apeGmsh
   computes the bilinear/area weights and emits **`-shape N1..NM`** instead (the element is
   `ndm`-parametric and accepts either form);
4. defaults to **U-only** (rotations free — Abaqus-consistent) and **`g0` on** (stress-free
   on staged add). Expose `enforce`/`explicit` pass-throughs; keep `-rot`/`-pressure`/`-matN`
   **off** (experimental — opt-in only, with a "not validated" note).

> [!info] Routing note for apeGmsh
> Use `LadrunoEmbeddedNode` for the **fork** path (hex/tet/quad, explicit, staged). For a
> **2D-native** projection where the upstream element already works, `ASDEmbeddedNodeElement`
> remains a valid fallback — but it is implicit-only and tri/tet-only, so prefer the fork
> element wherever a fork solid host is in play.

### 10.7 Pitfalls

- `-xi` / `-k auto` / `-wcap` require the **`-host`** form. `-xi` is **3D-only** (use
  `-shape` for 2D hosts).
- The single-point tie is a **pin** — for moment transfer embed a length (§10.5), not `-rot`.
- `-bipenalty` is **explicit-only**; auto-disabled under `-enforce al`.
- `-dtcr` without `-host` bounds only the **slave** side — unconservative on a light host;
  prefer `-wcap`.
- UR/UP/D9/`-corot` are **experimental / not validated** — don't ship them as production.
- A non-numeric `-k`/`-kr` value (a forgotten number) is now a **parse error**, not a silent
  zero penalty.

---

## 11. References & related

- **Design record:** [[23_ladruno_embedded_node_adr|ADR 23 — General LadrunoEmbeddedNode]]
  (decisions D1–D10, §14 v1 scope, §15 adversarial review).
- **Validation:** [[27_ladruno_embedded_node_validation_plan|the v1 validation plan]] (the
  two-tier battery: in-tree rigorous + Abaqus/LS-DYNA credibility).
- **Sibling:** [[LadrunoEmbeddedRebar_guide|LadrunoEmbeddedRebar]] (the anisotropic
  bar-in-solid coupling over the same kernel) + [[LadrunoBondSlip_guide|LadrunoBondSlip]].
- **apeGmsh contract:** [[ladruno_apegmsh_contract]] (the generator touch-points).
- **Reference codes:** Abaqus `*EMBEDDED REGION` (translations interpolated, rotations free,
  born consistent); LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID` (penalty, non-matching, optional
  slip; rigorous friction on a separate `*CONTACT` card); `ASDEmbeddedNodeElement`
  (Petracca/ASDEA — the tri/tet implicit-only precedent this exceeds).
- **Continuum rotation (UR):** de Souza Neto et al. (2008) §3.
- **Augmented Lagrangian:** Simo & Laursen (1992); Wriggers, *Computational Contact
  Mechanics* §6.4; Bertsekas (1982).
- **Bipenalty / explicit:** Hetherington & Askes (2009); Askes & Hetherington (2011);
  Belytschko et al. (2014) §6.7.
