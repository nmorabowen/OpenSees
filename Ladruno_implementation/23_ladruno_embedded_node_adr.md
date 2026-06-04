---
title: "ADR 23 — General LadrunoEmbeddedNode (node-to-host coupling)"
project: Ladruno
type: ADR + implementation plan
status: proposed (scope locked: U + UP + UR, shared-kernel sibling of LadrunoEmbeddedRebar)
element: element LadrunoEmbeddedNode
classTag: 33006 (element) — reserved here, free in the ELE_TAG ladruno band
related:
  - "[[20_ladruno_embedded_reinforcement_adr]]"
  - "[[LadrunoEmbeddedRebar_guide]]"
  - "[[22_rc3d_conformal_recipe]]"
  - "[[ladruno_element_contract]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - adr
  - element
  - embedded-node
  - mesh-tie
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
  - ssi
updated: 2026-06-03
---

# ADR 23 — General `LadrunoEmbeddedNode` (node-to-host coupling)

**Status:** proposed. **Scope locked** (session 2026-06-03): full ASD-parity DOF
coupling — **translations U + pressure UP + rotation UR**, in **both 2D and 3D
domains**, plus an opt-in **material-driven interface mode** (D9: per-direction uniaxial
materials → cohesive / unilateral-contact / elastic / approximate-friction interfaces) —
built as a **sibling** of [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]]
over a **shared header-only coupling kernel**. Deliverable = this ADR + the phased
implementation plan in §9.

**Companions:** [[20_ladruno_embedded_reinforcement_adr|ADR 20]] (the rebar element
this generalizes), [[LadrunoEmbeddedRebar_guide]] (shipped machinery reused here),
[[09_ladruno_brick]] / [[06_bezier_tet10]] (fork solid hosts), [[22_rc3d_conformal_recipe]]
(the non-element alternative). Memory: `project_rc3d_embedment_scoping`.

> [!info] Sibling: [[24_ladruno_coupling_constraints_adr|ADR 24 — coupling constraints (RBE2/RBE3/linear equation)]]
> This ADR embeds a node **into a host element** (continuum); ADR 24 constrains
> **node sets** (RBE3 weighted-average, RBE2 rigid, `Σcₖuₖ=0`) — no host element.
> They share this coupling kernel and were scoped in parallel (both first numbered
> "23"; the coupling-constraints one is now **24**). The **frame-on-continuum /
> DOF-mismatch** interface is served by both: here via host embedding **with rotation
> coupling (UR)**, there via **weighted node-set distribution** (RBE3 /
> `*CONSTRAINED_INTERPOLATION`). `LadrunoEmbeddedNode` (33006) is the concrete
> node-embedding element; ADR 24 references it rather than re-specifying it.

---

## 1. Context & problem

The fork ships `LadrunoEmbeddedRebar` (ELE 33005): a **penalty coupling element** that
ties **one node** to a **host element's nodes** through precomputed shape-function
weights `N_i(ξ)`, with an **anisotropic** traction (a bar axis `dir` + transverse
penalty + an axial bond-slip law). Stripped of the directional split, that element is
*already* a general node-to-host tie — the isotropic case is `k_axial = k_t = K`,
`D = K·I`, `dir` irrelevant.

What is **missing** is a first-class, fork-native **general node-embedding element**:
an isotropic position tie of an arbitrary node into a host continuum, optionally also
tying the node's **rotation** (to carry moment from beam/shell members into a solid)
and its **pressure** DOF (for u-p saturated-soil hosts). OpenSees has exactly one such
element today — upstream `ASDEmbeddedNodeElement` (Petracca/ASDEA) — and it has gaps
that block the fork's workflows:

- **Host-locked.** It **hardcodes 3-node-triangle / 4-node-tet** shape functions
  ([ASDEmbeddedNodeElement.cpp:64-80], `TRI_*`/`TET_*` paths). It **cannot** embed into
  an 8-node `LadrunoBrick` hex or a 10-node `BezierTet10` — the fork's primary solids.
- **Implicit-only in practice.** Penalty `m_K = 1.0e18` ([ASDEmbeddedNodeElement.h:124])
  with zero mass — the same explicit-`dt_cr` false-safe ADR 20 §10.6 dissected. No
  bipenalty, no `-cfl` self-report.
- **Ill-conditioned.** A raw `1e18·√V` penalty ([:1166]); no host-stiffness auto-scale,
  no augmented Lagrangian.

`LadrunoEmbeddedRebar` already solved every one of those (host-agnostic
`getInterpolationWeights`, `-kt auto`, `-enforce al`, bipenalty + `getExplicitCriticalTimeStep`
+ the `CriticalTimeStep` fold-in). The opportunity is to **reuse that validated
machinery** for the general node case the rebar element's anisotropy doesn't serve.

## 2. Goals / non-goals

**Goals:**
- A general **node-to-host** coupling element, host-agnostic, that works in **both 2D
  (`ndm=2`) and 3D (`ndm=3`)** domains (the dimensionality is a first-class requirement,
  see **D8**) — with the fork's `LadrunoBrick` / `BezierTet10` (3D) and `BezierTri6` /
  standard quad (2D) hosts, and any host implementing the shape-fn virtuals.
- **Isotropic translational tie (U)** — the mesh-stitch / point-embed primitive.
- **Pressure tie (UP)** — embed into u-p coupled hosts (all nodes u-p).
- **Rotational tie (UR)** — tie the constrained node's rotations to the host's continuum
  rotation field, so a beam/shell node embedded in a solid transmits moment.
- **Pluggable per-direction interface materials (D9)** — any local direction can carry a
  **uniaxial material** instead of a bare penalty, turning the element into a node↔host
  *interface*: cohesive separation, unilateral contact/gap, elastic bedding, bond — the
  bond-slip slot generalized to all directions.
- Inherit, for free, the rebar element's **penalty + AL** enforcement, **`-kt auto`**
  host-scaling, **bipenalty + explicit `dt_cr`** governance, serialization, responses,
  and energy accounting — via a **shared kernel** (no copy-paste).
- openseespy-native; set-to-set generation through apeGmsh's existing assembly path
  (it already emits `ASDEmbeddedNodeElement` for non-matching ties — this is the
  drop-in fork upgrade).

**Non-goals (v1):**
- **No *coupled* frictional contact.** The D9 interface materials are **uncoupled
  per-direction uniaxials** — faithful for cohesive / unilateral-normal / elastic / bond,
  but Coulomb friction (tangential capacity = `μ·N`, coupled to the live normal force) is
  only an *approximation* here (a fixed `ElasticPP` slip force). Rigorous coupled friction
  + projection kinematics + contact search is the **separate, scoped `LadrunoContact`
  element** (memory `project_ladruno_contact`); this element is its lightweight,
  search-free, non-matching-host *complement*, not its replacement.
- **No curved-host inverse map** — straight-sided hosts in v1 (same restriction as ADR 20
  D3); the global→ξ inverse map stays in the apeGmsh generator.
- **No `-enforce nitsche` / `transformation`** — same deferrals as ADR 20 §10.
- **UR is an approximate, mesh-limited moment tie** (continuum rotation = skew(∇u),
  resolution-limited like dowel action) — documented, not silently sold as exact.
- **No multi-node constrained set** — one constrained node per element (loop for many).
- Finite-strain `L_trib`/Jacobian stretch corrections deferred (small-strain host
  assumption, as ADR 20 D6).

## 3. The DOF-scope decomposition (drives the whole build)

The three coupling levels are **not** equal cost. The decisive fact, read from the
upstream UR path ([ASDEmbeddedNodeElement.cpp:1048-1173]):

| Level | Constraint enforced | Host query needed | Cost |
|-------|--------------------|--------------------|------|
| **U** (translations) | `u_c = Σ N_i(ξ) u_i^host` | **weights `N_i`** (have it) | trivial — degenerate isotropic case of the rebar kernel |
| **UP** (pressure) | `p_c = Σ N_i(ξ) p_i^host` (all nodes u-p) | **weights `N_i`** (have it) | cheap — one extra row, an ndf-compatibility check |
| **UR** (rotation) | `θ_c = skew(∇u_host)\|_ξ` (continuum rotation) | **gradients `∂N_i/∂x`** — **NOT available today** | moderate — needs a new vanilla `Element` virtual |

**UR is the only one needing new infrastructure.** The host continuum rotation is the
skew part of the displacement gradient — in the local triangle frame ASD builds it from
the **cartesian shape-function derivatives** `dNdX` ([:1104-1106], [:1153-1162]):

$$\theta^{loc}_x = \sum_i \partial_Y N_i\,u^z_i,\quad
  \theta^{loc}_y = -\sum_i \partial_X N_i\,u^z_i,\quad
  \theta^{loc}_z = \tfrac12\sum_i(\partial_X N_i\,u^y_i - \partial_Y N_i\,u^x_i),$$

then rotates to global with the host triangle's orientation `R`. **This requires the
host's `∂N_i/∂x` at ξ, which `getInterpolationWeights` does not return** — it is exactly
the `getInterpolationGradients` / `getDeformationGradient(ξ)` surface ADR 20 §10.5 flagged
as "research-grade, re-opens the base-class surface." For straight-sided simplex hosts
`∂N/∂x` is constant; for `LadrunoBrick`/`BezierTet10` it varies → must be evaluated at ξ.

**2D collapses the rotation to a single drilling DOF** (`Rz`): the continuum rotation is
the one in-plane curl `θ_z = ½(∂_x u_y − ∂_y u_x)` ([ASDEmbeddedNodeElement.cpp:858]
`TRI_2D_UR`), tied to the node's `Rz`; UP is the same scalar tie. Everything else (U gap,
weights, penalties, AL, bipenalty) is `ndm`-parametric and identical between 2D and 3D —
the rebar kernel is already written this way (`int ndm; // 2 or 3`).

**Consequence for the plan:** ship U + UP first (zero new vanilla surface beyond what's
already in `Element`), then add UR behind a **new vanilla virtual**
`Element::getInterpolationGradients(ξ, dNdx)` (default `return -1`; overridden on
`LadrunoBrick` + `BezierTet10`). This staging keeps the high-value, low-risk capability
(mesh tie, pressure tie, SSI) unblocked while the rotation tie — which re-opens the
base-class — lands as its own reviewed slice.

## 4. Decisions

### D1 — One new element, shared kernel with the rebar (no refactor of 33005)
Ship **`LadrunoEmbeddedNode`** (classTag **33006**, ELE band). Extract the
constitutive-agnostic numerics shared with `LadrunoEmbeddedRebar` into a **header-only
`LadrunoEmbeddedKernel.h`** (the proven [[project_ladruno_j2|LadrunoJ2Kernel]] pattern):
the gap/`B`-matrix assembly, the penalty/AL traction-and-tangent update, the `-kt auto`
formula, the bipenalty `m_p` and self-reported `dt_cr` formulas. Each element keeps its
**own** DOF bookkeeping, `setDomain`, and serialization (they genuinely differ — the
node element carries rotation/pressure DOFs on the constrained node; the rebar is
translations-only with a `dir`). The shipped 33005 contract and serialization version are
**untouched** — the kernel is additive; the rebar's call sites are refactored to route
through it **bit-identically** (a no-op refactor, gated by its existing Zone-A battery).

> **Why not** add an `-isotropic` mode to 33005 (rejected): bloats a shipped element,
> mis-names a general mesh tie "EmbeddedRebar", and UR/UP do not belong on a rebar.
> **Why not** a full strategy-merge into one `LadrunoEmbedded` (rejected for now): it
> touches the 33005 serialization version and contract — higher risk for no v1 gain. The
> kernel extraction captures ~90% of the reuse at ~10% of that risk; a later merge stays
> open if a third consumer appears.

### D2 — Constitutive: isotropic, per-DOF-class penalties
Traction is the **isotropic** degenerate of the rebar `D`:
- **U:** `t_u = K_u·g_u`, `g_u = u_c − Σ N_i u_i^host`, `D_u = K_u·I` (no `dir`, no
  axial/transverse split — the isotropic tie is **already objective** under rigid host
  rotation, so **no `-corot` is needed**, unlike the anisotropic rebar — record as a
  LEDGER_quirk).
- **UP:** `t_p = K_p·g_p`, `g_p = p_c − Σ N_i p_i^host` (scalar).
- **UR:** `t_r = K_r·g_r`, `g_r = θ_c − θ^{host}(ξ)` with `θ^{host}` the continuum
  rotation of §3; `B_r` carries the `R·dNdx` rotation block.

Each class gets its **own** penalty so they can be conditioned independently (ASD shares
one `m_K·√V` across translation and rotation — dimensionally crude). All three default to
**`auto`** host-scaling (D3).

**The penalty is the degenerate case of a per-direction constitutive (see D9).** The kernel
is written so each coupling direction's traction is `t_d = mat_d.stress(g·e_d)` with tangent
`mat_d.tangent` — a bare penalty is `mat_d = Elastic(K_d)`. This is the same abstraction the
rebar already uses for its axial slot (`LadrunoBondSlip` or perfect-bond penalty); D9 just
exposes it on the general element's directions.

### D3 — Auto-scaled penalties (host-stiffness matched, the rebar `-kt auto` generalized)
Reuse the rebar's lazy host-stiffness auto-scale ([[LadrunoEmbeddedRebar_guide]] §4):
`K_u = α_u · ‖K_host‖_∞ / lch · lch = α_u·‖K_host‖_∞`-class scaling via
`host->getInitialStiff()` + `host->getCharacteristicLength()`. **Rotation needs its own
dimensional scale** — `K_r ~ K_u·lch²` (so a unit rotation gap costs a comparable energy
to a unit translation gap over the element length); expose `-krAlpha`. Pressure `K_p`
scales to the host's pressure-block diagonal. Defaults: `α_u = 1e3` (the rebar default),
`α_r`, `α_p` tuned in the §8 sweep. The ASD `1e18` raw value is the **anti-pattern** —
documented loudly (porting users must not pass `1e18` as an `α`).

### D4 — Enforcement: penalty + augmented Lagrangian (`-enforce {penalty|al}`)
Identical to the rebar's `-enforce` family (kernel-shared). AL adds a per-element
multiplier `λ` (now spanning the active DOF classes — `ndm` + rot + p) with the **same
tangent** and a per-step **Uzawa** update in `commitState`. Headline win unchanged:
near-exact constraint at **moderate** `K` (no `K→∞`), fixing both conditioning and the
explicit `dt_cr` blow-up. `nitsche`/`transformation` rejected at parse time.

### D5 — Explicit capability: bipenalty + `dt_cr` self-report (kernel-shared)
Reuse the rebar's bipenalty (`-bipenalty {-dtcr dt | -wcap β}`): a mass penalty `m_p`
lumped **on the constrained (slave) node only** (DiagonalSOE-safe), with
`k_eff = max(K_u, K_r-scaled, K_p)`, the closed-form `eleResponse "dtcr" = 2√(m_p/k_eff)`,
and the already-shipped `Element::getExplicitCriticalTimeStep` → `CriticalTimeStep`
fold-in. Gated on `-enforce penalty`; off ⇒ `m_p≡0` ⇒ bit-identical, explicit-safe.
**This is the single biggest gain over `ASDEmbeddedNodeElement`** — it makes a general
node tie usable under `CentralDifferenceLadruno` / `ExplicitBathe`.

### D6 — Host-agnostic interpolation; UR adds a gradient virtual
- **U / UP:** `getInterpolationWeights(ξ, N)` — already on `Element` (vanilla virtual,
  [[LEDGER_vanilla_files]]), overridden today on the **3D** solids `LadrunoBrick`
  (trilinear) / `BezierTet10` (Bernstein). For **2D** host-query convenience, add the same
  override to `BezierTri6` and a standard quad (`FourNodeQuad`); **the explicit `-shape`
  form already works for any 2D host today** (the kernel is `ndm`-parametric), so 2D is not
  blocked on these overrides — they are the `-xi` ergonomic. `-host eleTag -xi ξ..` reuses
  the rebar parser path verbatim.
- **UR:** add **`virtual int Element::getInterpolationGradients(const Vector& xi,
  Matrix& dNdx)`** (default `return -1`; vanilla, ledgered). Override on `LadrunoBrick`
  (∂ of `0.125∏(1+ξ_Iξ)`) and `BezierTet10` (∂ of the quadratic Bernstein, via the
  existing private `shapeFunctions`/Jacobian) for 3D, and on `BezierTri6`/quad for 2D.
  Parser errors with a "host lacks gradient support → use U/UP only or supply `-dNdx`" hint
  when the host returns −1. An explicit `-dNdx` escape hatch (like `-shape`) feeds Zone-A
  unit tests a fake host.
- The **inverse map** (global point → ξ) stays in the apeGmsh generator regardless
  (ADR 20 §9) — the element only ever consumes `ξ`/weights.

### D7 — DOF-class activation is auto-detected from the constrained node's ndf (ASD model)
Mirror [ASDEmbeddedNodeElement.cpp:370-413]: read the constrained node's `ndf` at
`setDomain` and activate UR only if the node has rotational DOFs **and** the user set
`-rot`; activate UP only if all nodes are u-p **and** the user set `-pressure`. Flags are
**opt-in** (default = U only) so the common mesh-tie case stays minimal and bit-stable.
Allowed ndf: **2D = {2, 3}**, **3D = {3, 4, 6}**, matching ASD's table. **The 2D `ndf=3`
case is ambiguous — it is *either* `(u_x,u_y,R_z)` *or* `(u_x,u_y,p)`** — and is
disambiguated **only** by which flag the user set (`-rot` ⇒ drilling, `-pressure` ⇒ u-p);
this is a hard requirement (silent mis-binding otherwise), tested in §8.

### D8 — Dimensionality (2D + 3D) is first-class and `ndm`-parametric
The element runs in both `ndm=2` and `ndm=3`, resolved from the host/constrained-node
coordinates at `setDomain` (as the rebar element already does — `int ndm; // 2 or 3`).
Everything in the shared kernel (the U gap, weights, isotropic `D=K·I_{ndm}`, penalties,
AL multiplier sized to the active DOF classes, bipenalty `m_p`, `dt_cr`) is written once
over `ndm` and needs **no per-dimension code**. The only dimension-specific code is the
**UR rotation block**: 3D ties three rotations to `skew(∇u)` (a 3-vector), 2D ties the
single drilling `R_z` to the one in-plane curl `θ_z = ½(∂_x u_y − ∂_y u_x)` (§3). Both
read `∂N/∂x` from `getInterpolationGradients`. v1 supports straight-sided hosts in both
dimensions; curved-host inverse maps stay in the apeGmsh generator (D6).

### D9 — Material-driven interface mode (node↔interpolated-host-point `zeroLength`)
Generalize the rebar's axial-material slot to **all coupling directions**: any local
direction may carry a **uniaxial material** instead of a penalty. Mechanically the element
becomes a **`zeroLength` whose far end is a *point interpolated inside a non-matching host*
** rather than a real node — the novel capability (no existing OpenSees element gives a
material-driven interface against an interpolated host point).

- **Constitutive (per-direction, uncoupled).** With a local frame `{e_d}` (normal + tangents,
  D9.1), `t = Σ_d f_d(g·e_d)\,e_d`, `D = Σ_d k_d\,e_d⊗e_d`, where `(f_d, k_d)` is the
  uniaxial `stress`/`tangent` of `mat_d` (penalty ⇒ `mat_d = Elastic`). This is exactly the
  rebar's anisotropic `D` with each slot now a material — the rebar (axial=`LadrunoBondSlip`,
  transverse=penalty) and the isotropic node tie are both special cases.
- **What it models (faithfully):** **cohesive interface** (softening uniaxial on the normal
  ± shear, `-Gf`-regularized like `LadrunoBondSlip`); **unilateral contact / gap** (`ENT`
  elastic-no-tension or `ElasticPPGap` on the normal → stiff in compression, separates in
  tension); **elastic bedding / interface spring** (`Elastic`); **bond** (the existing slot).
- **Uncoupled limit (honest):** per-direction uniaxials cannot couple tangential capacity to
  the live normal force, so **Coulomb friction is only an approximation** (a fixed `ElasticPP`
  slip force, not `μ·N`). Rigorous coupled friction → the scoped `LadrunoContact` element
  (non-goal here). State this loudly in the guide.
- **D9.1 Local frame (opt-in).** Penalty-only ties (U/UP/UR) need **no** frame — isotropic
  `D=K·I` is frame-free and objective. A material on a *specific* direction needs one:
  user-supplied `-normal nx ny [nz]` (+ `-orient` tangents) **or** auto from the host face
  normal at ξ. The frame **co-rotates** through the existing `-corot` machinery (the normal
  follows the deformed host), which is *required* for a meaningful unilateral-contact normal.
- **D9.2 Enforcement/explicit interplay.** AL (`-enforce al`) augments only penalty
  directions — material directions carry physical force, untouched (mirrors the rebar's
  "AL augments transverse only, bond axial untouched"). Bipenalty `k_eff` takes the max over
  the **initial tangents** of all active directions (`mat_d.getInitialTangent()`), so a stiff
  contact-normal penalty is bounded for explicit just like `k_t`. Softening interface
  materials need DisplacementControl / IMPLEX (same caveat as bond, D4 of ADR 20).
- **Energy.** Physical interface dissipation (cohesive/bond ∮f·dg) reports in
  `bondEnergy`/`interfaceEnergy` (single-sourced from the material's own `energy` response);
  penalty directions stay in `penaltyEnergy` — the artificial/physical split already built
  for the rebar carries over verbatim.

### D9.x — Grammar (representative; firm at build time)
```
element LadrunoEmbeddedNode tag cNode {nHost h..|-host eleTag} {-shape N..|-xi ξ..}
        [-rot] [-pressure]                                  # extra DOF classes (D7)
        [-k Ku | -kt auto -ktAlpha a]                       # isotropic/penalty default
        [-normal nx ny [nz] [-orient ...]]                  # local frame (D9.1, material mode)
        [-matN tag] [-matT1 tag] [-matT2 tag]               # per-direction interface materials
        [-enforce {penalty|al}] [-corot ...] [-bipenalty ...]
```
`-matN/-matT*` select the material-driven interface; any direction without a `-mat*` falls
back to its penalty. `-matN` requires `-normal`.

## 5. Reference-code alignment
- **Abaqus** `*EMBEDDED ELEMENT` (translations-only DOF-elim, rotations of an embedded
  beam left free) → our **U** (and the deliberate choice to make **UR opt-in**, since
  Abaqus leaves it off by default).
- **LS-DYNA** `*CONSTRAINED_*_IN_SOLID` (penalty + mass scaling for explicit) → our **U +
  bipenalty**.
- **ASDEmbeddedNodeElement** (tri/tet penalty, optional UR/UP) → the **direct precedent**;
  we match its U/UR/UP capability and **exceed** it on host coverage (hex/tet),
  conditioning (`-kt auto` + AL), and explicit (bipenalty + `dt_cr`).
- **Continuum rotation tie (UR):** the host rotation `θ = ½ curl(u) = skew(∇u)` is the
  standard small-strain continuum rotation (de Souza Neto §3); tying a structural node's
  drilling/bending rotation to it is the established way to inject moment into a
  rotation-free continuum — approximate and mesh-resolution-limited (note honestly).
- **Interface materials (D9):** the `zeroLength`/`zeroLengthContact` family (uniaxial
  materials along local axes) — D9 is that idea with the far end an *interpolated host
  point*. `ENT`/`ElasticPPGap` for unilateral normal contact, a softening uniaxial for
  cohesive separation (cf. cohesive-zone models, [[15_lemaitre_ductile_damage_adr]] context),
  `ElasticPP` for approximate friction. Coupled friction → scoped `LadrunoContact`
  (ASDimplex IMPL-EX, the rigorous node-to-surface element).

## 6. Use cases (ranked, with the existing consumer)
1. **Non-matching solid↔solid mesh tie** — stitch independently-meshed parts (refined
   sub-block in a coarse host, part-to-part transitions). **apeGmsh already emits
   `ASDEmbeddedNodeElement` for exactly this** in its Part-based assembly — this element
   is the drop-in fork upgrade that adds hex/tet hosts + explicit + AL. *Strongest signal.*
2. **Beam/shell node embedded in a solid** — piles-in-soil, anchors, headed studs,
   embedded stiffeners/connectors carrying moment (**UR**).
3. **Embedded point entities** — lumped mass / applied load / control / recorder node tied
   into a solid interior without conforming mesh (**U**).
4. **Soil-structure interaction** — foundation/pile nodes into a soil continuum; **UP**
   for saturated u-p soil; explicit capability pairs with the absorbing-boundary/PML stack.
5. **Global-local submodeling** — drive a local model's boundary nodes from a global
   solid's displacement field.
6. **Embedded interface / contact (D9)** — cohesive debonding of an inclusion, a unilateral
   contact/gap between a node and a non-matching host face, an elastic bedding layer, or an
   approximate-friction sliding interface — all without conforming the meshes. The
   material-driven cousin of the dedicated `LadrunoContact` element.

## 7. Implementation sketch (Definition-of-Done registration checklist)
- `SRC/element/ladrunoEmbeddedNode/{LadrunoEmbeddedNode.{cpp,h},OPS_LadrunoEmbeddedNode.cpp,CMakeLists.txt}`.
- `SRC/element/ladrunoEmbeddedNode/LadrunoEmbeddedKernel.h` (header-only; **also** wired
  into `LadrunoEmbeddedRebar` as a no-op refactor — or sited under a shared path if
  preferred; decide at build time).
- `SRC/classTags.h`: add `ELE_TAG_LadrunoEmbeddedNode 33006` (next free in the ladruno
  ELE band; 33005 is the rebar).
- `SRC/element/Element.{h,cpp}`: add `virtual int getInterpolationGradients(const Vector&,
  Matrix&)` (default −1) — **UR slice only**; `// Ladruno` markers + [[LEDGER_vanilla_files]]
  row (vtable change, additive, recompile-all).
- Overrides: `getInterpolationGradients` on `LadrunoBrick` + `BezierTet10` (3D) and
  `BezierTri6` + `FourNodeQuad` (2D); the matching `getInterpolationWeights` overrides on
  the 2D hosts too (3D ones already exist) for the `-xi` convenience.
- Register: forward-decl + `functionMap` in `OpenSeesElementCommands.cpp`; **broker**
  `case 33006` in `FEM_ObjectBrokerAllClasses.cpp` + a serialization round-trip test.
- Ship-time obligations (CLAUDE.md): [[LEDGER_implementations]] row; `stamp_headers.py` on
  the new files (add the new dir to its GLOBS); `banner_features.txt` line +
  `patch_banner.py`; `testbed/manifest.yaml` row (the classTag fast-gate fails without it).

## 8. Validation plan (Zone-A; prescribed host displacements play the host role)
Mirrors the rebar battery ([[LadrunoEmbeddedRebar_guide]] §10), adding the rotation/pressure
legs. **Every mechanics leg runs in both `ndm=2` and `ndm=3`** (parametrized fixture); the
2D-specific cases are called out in items 5 and 9:
1. **U mechanics** — isotropic spring `K·g`; partition-of-unity force split to hosts +
   equilibrium; objectivity under prescribed rigid host rotation (passes **without**
   `-corot`, the isotropic-tie objectivity check); serialization round-trip. (2D + 3D.)
2. **`-host`/`-xi`** — real `LadrunoBrick` host: node auto-resolve, centroid weights all
   `0.125`, off-centroid weights match trilinear; real `BezierTet10` host.
3. **`-kt auto`** (and `-krAlpha`/`-kpAlpha`) — penalties scale with α; bounded κ.
4. **UP mechanics** — pressure gap → `p_c = Σ N_i p_i`; ndf-incompatibility rejected.
5. **UR mechanics** — `getInterpolationGradients` recovers `∂N/∂x` on hex/tet (3D) **and
   quad/tri (2D)** (FD-checked); continuum-rotation tie reproduces `θ = skew(∇u)` (3D
   3-vector) and `θ_z = ½(∂_x u_y − ∂_y u_x)` (2D drilling) for a prescribed host
   shear/bending field; **moment transmitted** from a constrained rotation into the host
   force split (2D + 3D); UR auto-off when the node has no rotational DOFs.
6. **`-enforce al`** — Uzawa drives the U/UR/UP gaps → 0 at moderate `K` (penalty leaves
   `P/K`); multiplier carries the load; `penalty` bit-identical.
7. **`-bipenalty`** — off bit-identical; `-dtcr` formula + `dtcr=dt`; `-wcap` inverse-β;
   explicit SDOF stable at `0.9×` / unstable at `1.1×` the target; `-cfl` governance.
8. **Zone-B (apeGmsh, follow-up)** — a real non-matching two-block tie (coarse hex host +
   refined sub-block) under implicit pushover and explicit dynamics; head-to-head vs
   `ASDEmbeddedNodeElement` on a tet/tri case where both apply (parity check).
9. **2D disambiguation** — `ndf=3` constrained node binds to drilling `R_z` under `-rot`
   and to pressure `p` under `-pressure` (never silently swapped); a 2D non-matching
   quad-block tie as the planar Zone-B analog.
10. **D9 interface materials** — `-matN Elastic(K)` reproduces the penalty bit-identically;
    `ENT`/`ElasticPPGap` on the normal opens a gap in tension and holds in compression;
    a softening uniaxial dissipates the right `∮f·dg` (reported in `interfaceEnergy`);
    `ElasticPP` tangential gives a fixed slip force (documented uncoupled-friction check);
    local frame co-rotates under `-corot`; bipenalty `k_eff` picks up the material initial
    tangents (2D + 3D).

## 9. Phased build order (the implementation plan)

**Phase 0 — kernel extraction (no behavior change).** Pull the rebar element's
constitutive-agnostic numerics into `LadrunoEmbeddedKernel.h`; refactor
`LadrunoEmbeddedRebar` to call it; prove **bit-identical** via its existing Zone-A battery
(`tests/test_ladrunoEmbeddedRebar_element.py`). *Risk: low — pure refactor under a green
test. Gate: rebar suite unchanged.*

**Phase 1 — `LadrunoEmbeddedNode` U + UP, 2D + 3D (no new vanilla surface).** New element
(33006) over the kernel: isotropic `K_u` (sized to `ndm`), optional `K_p`, `-enforce
{penalty|al}`, `-kt auto`, `-host`/`-xi`, full serialization + responses + bipenalty.
Reuses `getInterpolationWeights` only (works in 2D via the explicit `-shape` form even
before the 2D host overrides land). Ship the apeGmsh-tie use case. *Risk: low — degenerate
config of proven math + the trivial pressure row; `ndm`-parametric. Gate: Zone-A items
1–4, 6–7, 9 (all 2D + 3D).*

**Phase 2 — UR (rotation) behind the new gradient virtual, 2D + 3D.** Add
`Element::getInterpolationGradients` (vanilla, ledgered) + the host overrides
(`LadrunoBrick`/`BezierTet10` for 3D, `BezierTri6`/`FourNodeQuad` for 2D); add the `-rot`
path, the `R·dNdx` rotation `B`-block (3-rotation 3D / single-drilling 2D), and `-krAlpha`.
*Risk: moderate — re-opens the base-class surface (vtable/recompile-all) and the
continuum-rotation formulation. Gate: Zone-A item 5 (FD gradient check + moment-transmission,
2D + 3D).*

**Phase 2b — D9 material-driven interface mode.** Add the local-frame plumbing
(`-normal`/`-orient`, host-face-normal auto, co-rotation reuse) and the per-direction
`-matN`/`-matT*` slots over the kernel's `(e_d, mat_d)` abstraction; gap/cohesive/elastic
interfaces + approximate friction. *Risk: low–moderate — the kernel already drives a
uniaxial in the rebar axial slot; the new work is the general local frame + multi-slot
parsing + the contact-friendly defaults. No new vanilla. Gate: Zone-A item 10.*

**Phase 3 — apeGmsh generator + Zone-B + docs.** Teach the apeGmsh assembly path to emit
`LadrunoEmbeddedNode` (replacing/optioning `ASDEmbeddedNodeElement`); Zone-B items 8–9; the
`LadrunoEmbeddedNode_guide.md` (mirror the rebar guide); banner + ledger + manifest rows.

**Scope locked:** Phases 0–2b = the C++ (kernel + element + UR virtual + interface mode);
Phase 3 = tooling + docs. UR (Phase 2) is the only phase touching vanilla beyond the
already-shipped virtuals; Phase 2b adds no vanilla surface.

## 10. Risks
- **R1 (kernel extraction regresses 33005).** Mitigated: Phase 0 is gated bit-identical by
  the shipped rebar battery; the kernel is additive (no contract/serialization change).
- **R2 (UR continuum rotation is approximate).** Real but bounded: documented as
  mesh-limited (like dowel action, ADR 20 D6); opt-in; not the default. Quantify error in
  Zone-A item 5 vs an analytic shear/bending host field.
- **R3 (gradient virtual re-opens the base class).** Same surface ADR 20 §10.5 already
  weighed; additive default-`−1` virtual, recompile-all but no behavior change for
  non-overriding elements. Confined to Phase 2.
- **R4 (rotation/pressure penalty conditioning).** Distinct `α_r`/`α_p` (D3) instead of
  ASD's single shared `m_K` — swept in Zone-A item 3; AL (D4) gives the moderate-`K`
  escape if penalty conditioning bites.
- **R5 (explicit `dt_cr` for the rotation mode).** `k_eff` must include the
  `lch²`-scaled `K_r` so bipenalty bounds the rotation mode too — covered by D5's
  `max(...)` over active classes; tested in item 7 with `-rot` on.
- **R6 (D9 friction is uncoupled — mis-use risk).** Users may reach for `-matT` expecting
  true Coulomb friction. Mitigated by (a) loud docs ("approximate, fixed slip force, not
  `μ·N`"), (b) pointing to the scoped `LadrunoContact` for coupled friction, (c) the §8
  item-10 test that pins the uncoupled behavior as *intended*. Cohesive/gap/elastic/bond
  uses are exact — the caveat is friction-specific.

## 11. Ship-time obligations (CLAUDE.md, do in-PR)
- [[LEDGER_implementations]]: `LadrunoEmbeddedNode` row (ELE 33006) + the
  `getInterpolationGradients` note.
- [[LEDGER_vanilla_files]]: `Element.{h,cpp}` (the new gradient virtual),
  `LadrunoBrick`/`BezierTet10` (overrides), command/broker registrations — each with
  `// Ladruno` source markers.
- [[LEDGER_quirks]]: "isotropic node tie needs NO co-rotation (already objective); only
  the anisotropic rebar split does" + "UR continuum rotation needs host ∂N/∂x, not weights".
- `stamp_headers.py` GLOBS += the new element dir; `banner_features.txt` += a line +
  `patch_banner.py`; `testbed/manifest.yaml` += the 33006 row.
- Base PRs on `ladruno`; one logical PR per phase (avoid the stranded-commit trap,
  `feedback_stranded_commits_after_automerge`).

## 12. References & related
- **Generalizes:** [[20_ladruno_embedded_reinforcement_adr|ADR 20]] (the rebar element +
  the §10 `-enforce`/bipenalty/`-cfl` machinery reused here).
- **Precedent:** `ASDEmbeddedNodeElement` (Petracca/ASDEA) — tri/tet penalty with optional
  UR/UP; this ADR matches its DOF coverage and exceeds it on hosts/conditioning/explicit.
- **Kernel pattern:** [[project_ladruno_j2|LadrunoJ2Kernel]] (header-only shared numerics
  consumed by two element/material classes).
- **Interface materials (D9):** `zeroLength`/`zeroLengthContact` (uniaxials along local
  axes); `ENT`, `ElasticPPGap`, `ElasticPP` (OpenSees uniaxials); the scoped `LadrunoContact`
  element (memory `project_ladruno_contact`) for coupled frictional node-to-surface contact.
- **Continuum rotation:** de Souza Neto et al. (2008) §3; Wriggers (2006) §6 (penalty/AL).
- **Bipenalty / explicit:** Hetherington & Askes (2009); Belytschko et al. (2014) §6.7.
- **RC-3D context:** [[19_rc3d_modeling_recipe]], [[22_rc3d_conformal_recipe]].
