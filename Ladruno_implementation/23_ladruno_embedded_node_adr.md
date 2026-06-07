---
title: "ADR 23 — General LadrunoEmbeddedNode (node-to-host coupling)"
project: Ladruno
type: ADR + implementation plan
status: proposed, revised after adversarial review (wf_bbe77ee8 — 1 blocker fixed, Phase 0 GO); scope U + UP + UR + D9 interface, shared compiled-helper-kernel sibling of LadrunoEmbeddedRebar
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

> [!important] v1 scope (2026-06-07 credibility re-scope) — read **§14** first
> The build below (§1–§13) scoped the *full* U·UP·UR·D9 surface. After a credibility
> review we **re-scoped v1** to a focused, validated core: **U translational tie + `g0`
> stress-free birth + penalty/AL/bipenalty enforcement** (the only "validated, world-class"
> path), with **UP / UR / D9 / `-corot` kept in the tree but flagged EXPERIMENTAL (not
> validated)**. The benchmark framing (Abaqus `*EMBEDDED REGION`, LS-DYNA
> `*CONSTRAINED_BEAM_IN_SOLID`), the supported/experimental split, and the corrected
> must-fix list are in **§14**; the battery that earns the claim is
> [[27_ladruno_embedded_node_validation_plan]].

**Status:** proposed, **revised after adversarial review** (workflow `wf_bbe77ee8`, 33
agents: 25/26 findings survived verification, 1 confirmed blocker; verdict
**proceed-with-changes**, **Phase 0 = GO**). The must-fix list is folded in below; the
review record is §13. **Scope locked** (session 2026-06-03): full ASD-parity DOF
coupling — **translations U + pressure UP + rotation UR**, in **both 2D and 3D
domains**, plus an opt-in **material-driven interface mode** (D9: per-direction uniaxial
materials → cohesive / unilateral-contact / elastic / approximate-friction interfaces) —
built as a **sibling** of [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]]
over a **shared compiled helper kernel** (`LadrunoEmbeddedKernel.{h,cpp}` taking framework
handles as parameters — **not** a header-only OpenSees-free leaf; M3/KX-1). Deliverable =
this ADR + the phased implementation plan in §9.

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
- Inherit, via the **shared compiled helper kernel** (no copy-paste), the rebar's
  **penalty + AL** enforcement, **`-kt auto`** host-scaling (for the **translational `K_u`**
  — UR/UP scaling needs *additional* host queries the rebar does not supply, see D3),
  **bipenalty + explicit `dt_cr`** governance, responses, and energy accounting.
  (Serialization is *not* inherited — each element keeps its own, D1.)
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
| **UR** (rotation) | `θ_c = skew(∇u_host)\|_ξ` (continuum rotation) | **gradients `∂N_i/∂x`** *in addition to* weights `N_i` — gradients **NOT available today** | moderate — needs a new vanilla `Element` virtual |

**UR is the only one needing new infrastructure.** Note it needs **both** queries: the
UR `B`-matrix's *translational* rows still use the weights `N_i`
([ASDEmbeddedNodeElement.cpp:1145-1148]), the *rotational* rows use `∂N_i/∂x`
([:1155-1158]). Only the gradients are new. The host continuum rotation is the
skew part of the displacement gradient — in the local triangle frame ASD builds it from
the **cartesian shape-function derivatives** `dNdX` ([:1104-1106], [:1153-1162]):

$$\theta^{loc}_x = \sum_i \partial_Y N_i\,u^z_i,\quad
  \theta^{loc}_y = -\sum_i \partial_X N_i\,u^z_i,\quad
  \theta^{loc}_z = \tfrac12\sum_i(\partial_X N_i\,u^y_i - \partial_Y N_i\,u^x_i),$$

then rotates to global with the host triangle's orientation `R`. **The rotational rows
require the host's `∂N_i/∂x` at ξ, which `getInterpolationWeights` does not return** (the
translational rows still use `N_i` via `getInterpolationWeights`) — the gradients are exactly
the `getInterpolationGradients` / `getDeformationGradient(ξ)` surface ADR 20 §10.5 flagged
as "research-grade, re-opens the base-class surface." For straight-sided simplex hosts
`∂N/∂x` is constant; for `LadrunoBrick`/`BezierTet10` it varies → must be evaluated at ξ.

**UR sharpened caveat (UR-4):** on a **CST (3-node tri) / TET4 (4-node tet)** host `∂N/∂x`
is element-wise *constant*, so the UR constraint reduces to a single **element-wide
rigid-spin tie** — there is no intra-element rotational gradient and the moment couple's
lever arm is the element size, unresolvable by refining the surrounding mesh. **Moment-critical
embedments** (use case 2 — anchors, headed studs) therefore need a **higher-order host**
(e.g. `BezierTet10`, where `∂N/∂x` varies with position) or conforming refinement at the
connection. §8 item 5 reports the transmitted-moment error vs host *type and size*, not just
the displacement gap.

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

### D1 — One new element, shared compiled-helper kernel with the rebar (no refactor of 33005)
Ship **`LadrunoEmbeddedNode`** (classTag **33006**, ELE band). Extract the
constitutive-agnostic numerics shared with `LadrunoEmbeddedRebar` into a **shared compiled
helper translation unit `LadrunoEmbeddedKernel.{h,cpp}`**: the gap/`B`-matrix assembly, the
penalty/AL traction-and-tangent update, the `-kt auto` formula, the bipenalty `m_p` and
self-reported `dt_cr` formulas.

> **M3/KX-1 — NOT a header-only OpenSees-free leaf.** The earlier framing called this "the
> [[project_ladruno_j2|LadrunoJ2Kernel]] header-only pattern." The adversarial review
> refuted that analogy: unlike `LadrunoJ2Kernel` (pure scalar/tensor math, no framework
> types), this kernel's hot parts — `resolveAutoKt`, `resolveBipenalty`, `formBandTraction`
> (with `-corot`), the bipenalty `dt_cr` path — **require live `Domain*` / `Element*` /
> `Node**` handles** (it queries the host's `getInitialStiff`/`getCharacteristicLength` and
> the current node positions). So the kernel is a **normal compiled `.h`+`.cpp`** whose
> functions take those framework handles as **explicit parameters** — architecturally clean,
> captures the genuine reuse, avoids copy-paste, but it is a compiled helper, **not** a
> header-only leaf. (Only the truly handle-free scalar formulas — `k_eff`, `m_p`, `dt_cr` —
> could live header-only; that is a small fraction, so the compiled-TU framing is the right
> one.) The "~90%" below is a *reuse-vs-risk ratio* vs a full merge, not a code-fraction claim.

Each element keeps its **own** DOF bookkeeping, `setDomain`, and serialization (they
genuinely differ — the node element carries rotation/pressure DOFs on the constrained node;
the rebar is translations-only with a `dir`). The shipped 33005 contract and serialization
version are **untouched** — the kernel is additive; the rebar's call sites are refactored to
route through it **bit-identically** (a no-op refactor, gated by its existing Zone-A battery).

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

### D3 — Auto-scaled penalties (host-stiffness matched)
> **M2/HON-1/ES-2 — corrected against the shipped rebar.** The earlier
> `K_u = α_u·‖K_host‖_∞/lch·lch` was algebraically vacuous (the `lch` cancel) and falsely
> attributed a `getCharacteristicLength()` call to the rebar's `resolveAutoKt`, which in fact
> reads `host->getInitialStiff()` **only** (max absolute diagonal; no `lch`).

- **Translational `K_u` (inherited from the rebar `-kt auto`, [[LadrunoEmbeddedRebar_guide]] §4):**
  `K_u = α_u · max_i |K_host(i,i)|` via `host->getInitialStiff()` (max-abs-diagonal — the
  diagonal already carries `~E·lch` units). **No `getCharacteristicLength()` in this path.**
- **Rotational `K_r` (NEW work, not inherited):** rotation needs a *moment/rotation* scale,
  so `K_r = α_r · K_u · lch²` with `lch = host->getCharacteristicLength()` — this is the
  **first** `getCharacteristicLength()` call in the embedded-element family (the rebar kernel
  never makes it; it exists and is geometry-true on `LadrunoBrick`/`BezierTet10`/`BezierTri6`).
  Dimensionally required so a unit rotation gap costs energy comparable to a unit translation
  gap. Expose `-krAlpha`.
- **Pressure `K_p` (NEW work):** scales to the host's pressure-block diagonal (a new
  pressure-diagonal accessor / the u-p block of `getInitialStiff`); `-kpAlpha`.

Defaults: `α_u = 1e3` (the rebar default), `α_r`, `α_p` tuned in the §8 sweep. **UR/UP
auto-scaling is therefore NOT free** — it needs host queries beyond the rebar's `-kt auto`
(budgeted into Phase 1/2). The ASD `1e18` raw value is the **anti-pattern** — documented
loudly (porting users must not pass `1e18` as an `α`).

### D4 — Enforcement: penalty + augmented Lagrangian (`-enforce {penalty|al}`)
Same `-enforce` family as the rebar (kernel-shared *structure*). AL adds a per-element
multiplier `λ` spanning the active DOF classes (`ndm` + rot + p) with the **same tangent**
and a per-step **Uzawa** update in `commitState`. **Per-class increments (ES-3, spell out —
do not "match the rebar's pattern" blindly):** `Δλ_u = K_u·g_u`, `Δλ_r = K_r·g_r`,
`Δλ_p = K_p·g_p`. Headline win unchanged: near-exact constraint at **moderate** `K`
(no `K→∞`), fixing both conditioning and the explicit `dt_cr` blow-up.
`nitsche`/`transformation` rejected at parse time.

> **ES-3 — the rebar's corot λ-reprojection does NOT generalize to the isotropic tie.**
> The rebar re-projects `λ` onto the transverse plane each commit
> ([LadrunoEmbeddedRebar.cpp:304-311]) precisely because a rotating bond axis could leak the
> multiplier into the *axial material* slot. The **isotropic U/UP/UR penalty tie has no
> preferred axis**, so this reprojection has **no analog** and must **not** be imported here
> — `λ` accumulates without directional projection. (The one place it *does* return is the
> **D9 material-direction** case — see D9.2/M4, where a co-rotating frame can again leak a
> penalty multiplier onto a material axis.)

### D5 — Explicit capability: bipenalty with PER-DOF-CLASS (stiffness, inertia) pairs
> **M1/ES-1 (BLOCKER, the one the review caught) — a single scalar `k_eff`/`m_p` CANNOT
> bound the rotation or pressure mode.** The rebar's bipenalty is dimensionally homogeneous
> (everything translational): `getMass` lumps `m_p` only on DOFs `[0, ndm)`
> ([LadrunoEmbeddedRebar.cpp:498-521]) and `effectiveCouplingStiffness` returns a bare
> `max(k_axial, k_t)`. The node element carries **rotational/pressure DOFs whose penalty has
> different units** (`K_r` is moment/rotation), so a translational-only mass leaves those
> modes **unbounded** in explicit. Reusing the rebar's scalar pattern verbatim is wrong.

`-bipenalty {-dtcr dt | -wcap β}` with **matched (stiffness, inertia) pairs per active DOF
class**, each lumped on the constrained (slave) node only (DiagonalSOE-safe):
- **Translational:** `(K_u, m_p)` on DOFs `[0, ndm)` → `dt_u = 2√(m_p/K_u)` (the rebar form).
- **Rotational:** `(K_r, I_p)` with a **rotational inertia** `I_p = m_p·lch²` lumped on the
  slave's rotation DOF indices `[ndm, ndm+n_rot)` → `dt_r = 2√(I_p/K_r)`. **Self-consistency:**
  under D3's `K_r = K_u·lch²`, `dt_r = 2√(m_p·lch²/(K_u·lch²)) = 2√(m_p/K_u) = dt_u` — the
  translational budget self-bounds the rotation mode **iff** `getMass` actually lumps
  `I_p = m_p·lch²` on the rotation DOFs (not zero).
- **Pressure:** `(K_p, m_pp)` on the pressure DOF, **or** omit pressure bipenalty and require
  implicit enforcement for u-p hosts (decide in Phase 1).

The `getMass` loop must therefore extend past `ndm` to the active rotation/pressure DOF
indices. Report `dt_cr = min` over active classes of `2√(m_class/k_class)`; the
`eleResponse "dtcr"` self-report and the already-shipped `Element::getExplicitCriticalTimeStep`
→ `CriticalTimeStep` fold-in carry that min. Gated on `-enforce penalty`; off ⇒ all
`m≡0` ⇒ bit-identical, explicit-safe. **This (corrected) is the single biggest gain over
`ASDEmbeddedNodeElement`** — it makes a general node tie usable under
`CentralDifferenceLadruno` / `ExplicitBathe`.

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
Allowed ndf: **2D = {2, 3}**, **3D = {3, 4, 6}**, matching ASD's table.

**M5/D7-1 — the disambiguation is a PARSE-TIME guard, not setDomain precedence.** The 2D
`ndf=3` case is ambiguous — *either* `(u_x,u_y,R_z)` *or* `(u_x,u_y,p)`. ASD resolves this by
**rejecting `-rot` and `-pressure` together at parse time** ([ASDEmbeddedNodeElement.cpp:277],
"Cannot use both -rot and -p"). We do the same: **mutual exclusion is enforced at parse time
with a clear error** (primary mechanism); the `setDomain` `if/else if` ndf precedence is only
a defensive backstop. Tested in §8 item 9 (passing both flags ⇒ parse error, not silent
precedence). NB ASD's flag is `-p`; this element uses `-pressure` (the guard wording tracks
the new name).

**D8-1 — 3D is unambiguous (guard against a copy-paste bug):** in 3D `ndf=4` ⇒ UP and
`ndf=6` ⇒ UR with no flag disambiguation needed, so the 2D flag-priority logic must **not**
be applied in the 3D branch.

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
  **D9-5 — inherited tangent approximation:** the co-rotating frame inherits the rebar's
  **dropped `∂e_d/∂u` consistent-tangent term** ([LadrunoEmbeddedRebar.h:152-154]): for a
  stiff contact normal under large *per-step* host rotation the tangent is inexact in the same
  way the rebar's is — **frame-objective for explicit** (EICR is exact there), **converges
  under step-halving for implicit**, but may **slow Newton** for large-rotation contact steps.
  Port that caveat verbatim into the guide (also R7).
- **D9.2 Enforcement/explicit interplay.**
  - **M4/D9-3 — AL must re-project the multiplier off material directions.** AL augments only
    the penalty directions; the material directions carry physical force. But mirroring the
    rebar ([LadrunoEmbeddedRebar.cpp:304-311], `la = λ·e_mat; λ −= la·e_mat`), at each commit
    the accumulated `λ` must be **re-projected to subtract any component along each
    material-driven direction `e_d`**. Under `-corot`, **frame rotation alone** can shift a
    previously-transverse multiplier onto a material axis — so the purge is needed **even when
    no Uzawa step fires on the material slot**. (This is the one place ES-3's "no reprojection"
    does *not* apply — it returns precisely for the D9 material-direction case.) New Zone-A
    item 11.
  - **M6/ES-4/HON-2 — bipenalty `k_eff = max(initial tangents)` only for NON-STIFFENING
    materials.** `mat_d.getInitialTangent()` is a safe `k_eff` upper bound **only** when the
    material is monotone-non-stiffening (`ENT`, `ElasticPPGap` compression side, softening
    cohesive, `ElasticPP`). A **stiffening** uniaxial (e.g. a gap material whose re-contact
    tangent exceeds its initial/open tangent) silently breaks the closed-form `dt_cr` (the
    once-and-latch resolve is safe for `LadrunoBondSlip` but not in general). **Prohibit
    stiffening materials on any bipenalty direction** (parser error / guide), or require
    `-dtcr` set from the closed-contact stiffness, or use `-enforce al` (no `dt_cr` concern).
    NB the D9 grammar confines `-mat*` to *translational* normal/tangent directions, so the
    ES-1 rotation-unit problem does **not** arise in the material slots.
- **Energy (D9-4 — unit contract, not "verbatim").** Physical interface dissipation
  (cohesive/bond `∮f·dg`) reports in `interfaceEnergy`/`bondEnergy` (single-sourced from the
  material's own `energy` response); penalty directions stay in `penaltyEnergy`. The
  artificial/physical split carries over **in structure**, but the **unit contract differs**:
  a D9 `mat_d` is driven by the displacement gap `g·e_d` (metres) and **returns force**
  (`stress()` in N, not N/m²), so `energy()` is already in Joules and needs **no**
  `bondScale = perimeter·L_trib` converter — unlike the rebar, whose `bondMat` works in stress
  units and the element applies `bondScale`.

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
- `SRC/element/ladrunoEmbeddedNode/LadrunoEmbeddedKernel.{h,cpp}` (a **compiled helper TU**,
  M3/KX-1 — its functions take `Domain*`/`Element*`/`Node**` handles as parameters; **also**
  wired into `LadrunoEmbeddedRebar` as a no-op refactor — or sited under a shared path if
  preferred; decide at build time).
- `SRC/classTags.h`: add `ELE_TAG_LadrunoEmbeddedNode 33006` (next free in the ladruno
  ELE band — `classTags.h:924` is the last occupied ELE tag, `LadrunoEmbeddedRebar=33005`;
  REG-2. Per-registry bands ⇒ no collision with INTEGRATOR `33006`).
- `SRC/element/Element.{h,cpp}`: add `virtual int getInterpolationGradients(const Vector&,
  Matrix&)` (default −1) — **UR slice only**; `// Ladruno` markers + [[LEDGER_vanilla_files]]
  row (vtable change, additive, recompile-all).
- Overrides: `getInterpolationGradients` on `LadrunoBrick` + `BezierTet10` (3D) and
  `BezierTri6` + `FourNodeQuad` (2D); the matching `getInterpolationWeights` overrides on
  the 2D hosts too (3D ones already exist) for the `-xi` convenience.
- Register (REG-1, both steps, per the rebar precedent): forward-decl + **both** PascalCase
  **and** camelCase aliases in `functionMap` (`"LadrunoEmbeddedNode"` *and*
  `"ladrunoEmbeddedNode"`, cf. `OpenSeesElementCommands.cpp:667-668`); **`add_subdirectory(ladrunoEmbeddedNode)`
  in `SRC/element/CMakeLists.txt`** (the per-dir `CMakeLists.txt` is inert without it, cf.
  line 37); **broker** `case 33006` in `FEM_ObjectBrokerAllClasses.cpp` + a serialization
  round-trip test.
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
   shear/bending field; **transmitted-moment error reported vs host TYPE and SIZE** (UR-4 —
   CST/TET4 give a single rigid-spin tie; quantify vs `BezierTet10`), not just the
   displacement gap; UR auto-off when the node has no rotational DOFs.
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
    `ElasticPP` tangential gives a fixed slip force (documented uncoupled-friction check);
    local frame co-rotates under `-corot`. **Split legs (ES-5):** the *softening* leg
    (`∮f·dg` in `interfaceEnergy`) runs under **DisplacementControl**; the *bipenalty-`k_eff`*
    leg uses an **elastic/pre-peak** material only. **M6 gap sub-check:** `ElasticPPGap` on the
    normal with `-bipenalty -dtcr` set from the *open* tangent — confirm the explicit SDOF
    stays stable as the gap closes (or that the bound must be set from the closed-contact
    stiffness). (2D + 3D.)
11. **AL × material × corot (M4/D9-3)** — `-enforce al` + one penalty direction + one material
    direction + `-corot`: after N commit steps under large frame rotation, the projection of
    `λ` onto the material axis stays below tolerance (the reprojection purge works); without
    the purge it leaks (negative control).

## 9. Phased build order (the implementation plan)

**Phase 0 — kernel extraction (no behavior change). GO (review-cleared).** Pull the rebar
element's constitutive-agnostic numerics into a **compiled helper TU
`LadrunoEmbeddedKernel.{h,cpp}` whose functions take `Domain*`/`Element*`/`Node**` handles as
explicit parameters** (M3/KX-1 — bake this in from the outset; do **not** attempt a
header-only OpenSees-free leaf, or `resolveAutoKt`/`resolveBipenalty`/`formBandTraction` will
hit the wall that they cannot be made handle-free and force a mid-phase re-architecture).
Refactor `LadrunoEmbeddedRebar` to call it; prove **bit-identical** via its existing Zone-A
battery (`tests/test_ladrunoEmbeddedRebar_element.py`). *Risk: low — pure refactor under a
green test; none of the M1–M6 blockers touch the rebar's behavior. Gate: rebar suite
unchanged.*

**Phase 1 — `LadrunoEmbeddedNode` U + UP, 2D + 3D (no new vanilla surface).** New element
(33006) over the kernel: isotropic `K_u` (sized to `ndm`), optional `K_p`, `-enforce
{penalty|al}`, `-kt auto`, `-host`/`-xi`, full serialization + responses + bipenalty.
Reuses `getInterpolationWeights` only (works in 2D via the explicit `-shape` form even
before the 2D host overrides land). Ship the apeGmsh-tie use case. *Risk: low — degenerate
config of proven math + the trivial pressure row; `ndm`-parametric. Gate: Zone-A items
1–4, 6–7, 9 (all 2D + 3D).*

**Phase 2 — UR (rotation) behind the new gradient virtual, 2D + 3D. BUILT (2026-06-04).** Added
`Element::getInterpolationGradients` (vanilla, ledgered) + the host overrides on the two 3D fork
solids `LadrunoBrick` (`shp3d`) / `BezierTet10` (`computeJacobian`); the `-rot` path, the
continuum-rotation `B`-block (3-rotation 3D `θ=½Σ(∇N_i×u_i)` / single-drilling 2D), `-kr {val|auto}` +
`-krAlpha`, the M5 `-rot`/`-pressure` parse guard, and the M1/ES-1 per-DOF-class bipenalty
(`I_p=K_r·(dt/2)²` on the rotation DOFs ⇒ `dt_r=dt_u`). Gradients reach the element via `-dNdx`
(explicit, any host), `-gradXi ξ` (host query), or `-xi` auto-query on a gradient-capable host.
**Decision vs the ADR §3 sketch:** for a 3D VOLUME host all 9 ∇u components are available, so the
element uses the **pure continuum rotation `skew(∇u)`** (½ on all three, frame-objective, no local
`R`) rather than ASD's planar-surface mixed convention (full-slope bending + ½ drilling). **2D host
gradient overrides (`BezierTri6`/`FourNodeQuad`) DEFERRED** — 2D UR runs today via the explicit
`-dNdx` form (the kernel is `ndm`-parametric; this is exactly how `-shape` already serves 2D weights),
so 2D is not blocked. *Risk realized: moderate (vtable/recompile-all). Gate: Zone-A item 5 — FD/analytic
continuum-rotation recovery (2D drilling + 3D, `-dNdx`), real LadrunoBrick + BezierTet10 shear-field
moment transmission (`-xi`), M5 guard, UR auto-off, per-class bipenalty dt_cr, AL rotation gap→0.*

**Phase 2b — D9 material-driven interface mode. BUILT v1 (2026-06-04).** Added the
local-frame plumbing (`-normal`/`-orient` → orthonormal frame) and the per-direction
`-matN`/`-matT1`/`-matT2` slots: each translational local direction carries a uniaxial
material (penalty fallback K_u) — `t=Σ_d f_d(g·e_d) e_d`, `D=Σ_d k_d e_d⊗e_d` in FORCE
units (D9-4, no bondScale). Full material lifecycle (getCopy/commit/revert/setTrialStrain)
+ broker serialization; AL re-projects λ off material directions (M4/D9-3); bipenalty
`k_eff=max(K_u, init tangents)` with a stiffening-material warning (M6/ES-4). Models
cohesive / unilateral-gap (ENT/ElasticPPGap) / elastic-bedding / bond; approximate
friction (fixed ElasticPP slip). v1 used the REFERENCE frame. No new vanilla. *Gate:
Zone-A item 10 (fixed-frame subset) — Elastic == penalty, ENT unilateral gap, ElasticPP
fixed slip, AL reprojection, bipenalty k_eff, 2D, -matN-requires-normal, local responses.*

**Phase 2b v2 — D9 frame co-rotation (`-corot`). BUILT (2026-06-07).** The local material
frame now CO-ROTATES with the host continuum rotation `θ_host = skew(∇u)|_ξ` at the embedded
point: `frameCur = R(θ_host)·frame` (3D Rodrigues exponential map of the axial vector / 2D
drilling planar rotation), so a directional contact normal follows the deformed host. Reuses
the **Phase-2 continuum-rotation machinery** verbatim — `θ_host` from the host gradients
`∂N/∂x` via the shared `rotOper` (factored into `hostContinuumRotation`), supplied by
`-dNdx`/`-gradXi`/`-xi` exactly as UR. `-corot` is material-mode-only (parse-time reject
otherwise — the isotropic/penalty tie is already frame-objective) and resolves `corotActive`
in `setDomain` (needs `gradN`). The co-rotated `frameCur` drives the traction/tangent, the AL
λ-reprojection (so **frame rotation alone** can no longer leak the multiplier onto a material
axis — M4/D9-3 now exercised under rotation), and the `localGap`/`localForce`/`normal`
responses. The **dropped `∂e_d/∂u` consistent-tangent term (R7/D9-5)** is inherited from the
rebar: residual exact, tangent inexact (frame-objective for explicit, step-halving for
implicit, may slow Newton on stiff-normal large-per-step-rotation contact); NB when the host
DOFs are prescribed the cNode tangent is in fact exact (`θ_host` is host-only). Serialization
extends to `corot` + `gradN` (hdr→27; gradN now sent whenever UR **or** `-corot`). No new
vanilla. *Gate: Zone-A — frame normal = `R(θ_z)·n` under a prescribed-host drilling rotation
(2D) with a `-corot`-off negative control; `-corot`-requires-material parse guard; **item 11**
(AL × material × corot — λ·e_0^cur stays ~0 under frame rotation).* **DEFERRED still:** the
auto host-face-normal (ill-defined for an interior volume embedding — belongs to
`LadrunoContact`'s surface projection); large-rotation rigorous contact remains `LadrunoContact`.

**Phase 2c — initial-gap (offset) capture for stress-free staged activation. BUILT
(2026-06-07).** v1 computed every gap as a pure **trial-displacement difference**
(`g = u_c − Σ N_i u_host`, and likewise `g_p`/`g_r`), so the penalty enforced an **absolute**
tie `u_c = Σ N_i u_host`. An element added **mid-analysis to an already-deformed host** (staged
construction) therefore activated with `g = −N·u_host ≠ 0` and **yanked the slave** by the full
accumulated host displacement (a spurious force spike) — the parent `ASDEmbeddedNodeElement`
avoids this with its `m_U0` capture, which v1 dropped. **Decision D10 (capture ON by default;
`-absolute` opts out):** at `setDomain` — after the DOF layout + `upActive`/`rActive` + host
gradients are resolved — the element captures each **active** gap **once** (`g0`/`gp0`/`gr0`,
guarded by `g0Computed`) and thereafter drives **all** traction from the **relative** gap
`(g − g0)`. Implemented by subtracting the offset **inside** `computeGap`/`computeGapP`/
`computeGapR`, so **every** consumer sees the relative gap: penalty force, tangent, AL Uzawa,
the recorder `gap`/`pgap`/`rgap` probes, **and — the stress-free crux — the D9 material
`setTrialStrain`**, so interface materials are born **unstrained** (a cohesive law at its origin,
a gap material open, bond un-slipped), not merely the penalty zeroed. Composes with `-corot`:
`g0` is a **global** vector subtracted **before** the frame projection `g·e_d`, so at activation
`(g−g0)=0` ⇒ all material strains zero regardless of the (already-rotated) frame. **Force-free
AND stress-free**, not just force-free. The capture is a **no-op at the undeformed state**
(`g0=0`) ⇒ the entire v1 battery stays byte-identical. `-absolute` (alias `-noInitGap`) keeps the
legacy absolute tie (e.g. a deliberate snap-to-host). UR is linearized, so `gr0` subtraction is
**exact for small inter-stage rotation, approximate for large** (consistent with UR's
mesh-limited scope, R2). Serialization adds `initGapCapture`/`g0Computed` + `g0`/`gp0`/`gr0`
(hdr→29; `g0Computed` restored so `recvSelf` does **not** re-capture). New `initGap`/`offset`
diagnostic response. **No new vanilla.** *Gate: Zone-A — staged born-stress-free on a deformed
`LadrunoBrick` (zero force + zero relative gap at activation; offset = −host-centroid-disp; slave
not yanked on the next step), `-absolute` negative control (full-gap jolt + slave yanked to the
host point), undeformed no-op (byte-identical), FE_Datastore send/recv round-trip of the captured
state.*

**Phase 3 — apeGmsh generator + Zone-B + docs.** Teach the apeGmsh assembly path to emit
`LadrunoEmbeddedNode` (replacing/optioning `ASDEmbeddedNodeElement`); Zone-B items 8–9; the
`LadrunoEmbeddedNode_guide.md` (mirror the rebar guide); banner + ledger + manifest rows.

**Scope locked:** Phases 0–2c = the C++ (kernel + element + UR virtual + interface mode +
the v2 `-corot` frame co-rotation + the 2c initial-gap capture); Phase 3 = tooling + docs. UR
(Phase 2) is the only phase touching vanilla beyond the already-shipped virtuals; Phases 2b
(`-corot`) and 2c (initial-gap capture) add **no vanilla surface**.

> **Review gate (post-`wf_bbe77ee8`):** **Phase 0 = GO now.** **Hold Phase 1+ until M1–M6
> (§13) are reflected in the code design** — in particular M1/ES-1 (per-class bipenalty
> inertia) gates the Phase-1 `K_u`/bipenalty `getMass` and the Phase-2 `K_r`; M2/HON-1
> (the corrected `K_u = α_u·max|K_host(i,i)|`, no `lch`) gates the Phase-1 auto-scale.
> Building from the *pre-revision* D3/D5 text would ship a no-op `lch` cancel and an
> unbounded rotation mode. (This ADR revision folds them in, so the text is now correct.)

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
- **R5 (explicit `dt_cr` for the rotation/pressure modes) — the confirmed blocker (M1/ES-1).**
  A single scalar `k_eff`/`m_p` lumped on translational DOFs **cannot** bound the
  rotation/pressure modes (different units; `getMass` lumps only `[0,ndm)`). Resolved by D5's
  **per-class (stiffness, inertia) pairs** — rotational `I_p = m_p·lch²` lumped on the rotation
  DOFs (then `dt_r = dt_u` under `K_r = K_u·lch²`), pressure `(K_p, m_pp)` or implicit-only.
  Tested item 7 with `-rot`/`-pressure` on. **Must be in the Phase-1 `getMass` from day one.**
- **R6 (D9 friction is uncoupled — mis-use risk).** Users may reach for `-matT` expecting
  true Coulomb friction. Mitigated by (a) loud docs ("approximate, fixed slip force, not
  `μ·N`"), (b) pointing to the scoped `LadrunoContact` for coupled friction, (c) the §8
  item-10 test that pins the uncoupled behavior as *intended*. Cohesive/gap/elastic/bond
  uses are exact — the caveat is friction-specific.
- **R7 (D9 co-rotated-frame tangent is inexact — D9-5).** The `-corot` interface normal drops
  the `∂e_d/∂u` consistent-tangent term (inherited from the rebar): frame-objective for
  explicit, converges under step-halving for implicit, but may **slow Newton** for large
  per-step rotation in stiff-normal contact. Documented in D9.1 + the guide; not a
  correctness defect (residual stays exact).

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
- **Kernel pattern (with a caveat, M3/KX-1):** [[project_ladruno_j2|LadrunoJ2Kernel]] is the
  *spirit* (one numeric core, two consumers) but **not** the literal model — J2Kernel is
  header-only/OpenSees-free, whereas `LadrunoEmbeddedKernel` is a **compiled helper TU** taking
  framework handles as parameters (its auto-kt/bipenalty/corot paths need `Domain*`/`Node**`).
- **Interface materials (D9):** `zeroLength`/`zeroLengthContact` (uniaxials along local
  axes); `ENT`, `ElasticPPGap`, `ElasticPP` (OpenSees uniaxials); the scoped `LadrunoContact`
  element (memory `project_ladruno_contact`) for coupled frictional node-to-surface contact.
- **Continuum rotation:** de Souza Neto et al. (2008) §3; Wriggers (2006) §6 (penalty/AL).
- **Bipenalty / explicit:** Hetherington & Askes (2009); Belytschko et al. (2014) §6.7.
- **RC-3D context:** [[19_rc3d_modeling_recipe]], [[22_rc3d_conformal_recipe]].

## 13. Adversarial review record (workflow `wf_bbe77ee8`, 2026-06-03)
6 reviewers (one per load-bearing claim) → each finding adversarially verified against source
→ lead synthesis. **33 agents, 26 findings, 25 survived verification, 1 confirmed blocker.
Verdict: proceed-with-changes; Phase 0 = GO.** Must-fixes folded into the text above:

| ID | Sev | Status | Fix | Where |
|----|-----|--------|-----|-------|
| **ES-1 (M1)** | blocker | confirmed | per-DOF-class `(K, inertia)` pairs; rotational `I_p=m_p·lch²` lumped on rotation DOFs (scalar `k_eff`/`m_p` can't bound rotation/pressure) | D5, R5 |
| **HON-1/ES-2/KX-3 (M2)** | major | confirmed/partial | `K_u = α_u·max\|K_host(i,i)\|` (getInitialStiff diagonal, **no `lch`** — the old formula cancels); `K_r = α_r·K_u·lch²` is NEW work (first `getCharacteristicLength()` call) | D3, Goals |
| **KX-1 (M3)** | major | confirmed | kernel is a **compiled helper TU** taking framework handles, **not** a header-only J2Kernel leaf | Status, D1, §7, §9 P0, §12 |
| **D9-3/ES-3 (M4)** | major | confirmed | AL must re-project `λ` off material directions (needed under `-corot` even with no Uzawa step); the *isotropic* tie keeps no such reprojection | D4, D9.2, item 11 |
| **D7-1 (M5)** | major | confirmed | parse-time `-rot`/`-pressure` mutual-exclusion guard (primary, not setDomain precedence) | D7, item 9 |
| **ES-4/HON-2 (M6)** | major | partial | bipenalty `k_eff=max(initial tangents)` only for non-stiffening materials; gap-material guard | D9.2, item 10 |
| UR-1 | minor | confirmed | UR needs `∂N/∂x` **in addition to** `N_i` (translational rows still use weights) | §3 |
| UR-4 | minor | partial | CST/TET4 `∂N/∂x` constant ⇒ single rigid-spin tie; moment-critical needs higher-order host | §3, item 5 |
| D9-5 | minor | partial | port the rebar's dropped-`∂e_d/∂u` tangent caveat into D9.1 | D9.1, R7 |
| D9-4 | minor | partial | D9 material returns **force** (no `bondScale`); energy split carries over "in structure" not verbatim | D9.2 |
| KX-4 | minor | partial | serialization is per-element, not "inherited for free" | Goals |
| D8-1 | minor | partial | 3D ndf mapping is unambiguous — don't apply 2D flag-priority in the 3D branch | D7 |
| REG-1 | minor | confirmed | register PascalCase **+** camelCase aliases; `add_subdirectory` in parent CMake | §7 |
| REG-2 | minor | confirmed | 33006 free (`classTags.h:924` = last ELE tag 33005) | §7 |
| UR-2, UR-3, D9-1, D9-2, D9-6, D6-1 | minor | confirmed | **positive verifications** — core formulas/architecture sound, no change | — |

Source cross-checks (worktree `great-chebyshev-643d1a`): `LadrunoEmbeddedRebar.cpp:150-185`
(resolveAutoKt — getInitialStiff diagonal only), `:190-195` (scalar k_eff), `:304-311` (corot
transverse λ-reprojection), `:498-521` (getMass lumps `[0,ndm)` only); `Element.h:62/71/87`
(getCharacteristicLength/getInterpolationWeights/getExplicitCriticalTimeStep present,
getInterpolationGradients absent); `ASDEmbeddedNodeElement.cpp:277` (parse-time `-rot`/`-p`
guard); `classTags.h:924` (33005 last ELE tag, 33006 free).

---

## 14. v1 scoping & status — supported core vs experimental (credibility re-scope, 2026-06-07)

> [!important] Read this first if you are deciding what to rely on
> `LadrunoEmbeddedNode` grew to **five** flag-gated capabilities (U · UP · UR · D9 ·
> enforcement). Even gated, the wide surface **dilutes the validation story** and was the
> wrong way to earn the original goal: a **world-class embedment** that starts closing the
> gap with Abaqus / LS-DYNA. This section re-scopes v1 to a **focused, validated core**; the
> rest stays in the tree but is **explicitly experimental (not validated)**. *Focused-and-
> validated beats wide-and-half-tested.* The validation battery that earns the "world-class"
> claim is [[27_ladruno_embedded_node_validation_plan]].

### 14.1 The benchmark — what a world-class embedment actually is
- **Abaqus `*EMBEDDED REGION` / `*EMBEDDED ELEMENT`:** the **translational** DOFs of embedded
  nodes are constrained to (interpolated from) the host element; embedded **rotations are left
  FREE** (the embedded beam/shell carries its own bending); embedded nodes are **born
  consistent** with the host at activation (no jolt). → our **U tie + g0 stress-free birth**,
  with rotations free by default.
- **LS-DYNA `*CONSTRAINED_BEAM_IN_SOLID` / `*CONSTRAINED_..._IN_SOLID`:** a **penalty**,
  **non-matching** coupling with normal + tangential terms; **slip/bond is an OPTIONAL add-on**
  and **rigorous friction is a SEPARATE `*CONTACT` card**. → our **penalty/AL/bipenalty**
  enforcement is the core; D9 cohesive/contact/friction/bond is the *optional* (experimental)
  add-on, and rigorous frictional contact belongs to the separate `LadrunoContact` lineage.

The decisive differentiator over the one OpenSees precedent (`ASDEmbeddedNodeElement`,
implicit-only, raw `1e18` penalty, tri/tet-only) is: **host-agnostic** (hex/tet/quad),
**`-kt auto` conditioning** (no `1e18`), **AL** (near-exact at moderate `K`), and **clean
EXPLICIT** capability (bipenalty + self-reported `getExplicitCriticalTimeStep` + the
`CriticalTimeStep` fold-in). That implicit-**and**-explicit conditioning story is the headline
the validation plan foregrounds.

### 14.2 v1 SUPPORTED CORE (ships as "validated")
The **only** path that carries the validated, world-class claim:
1. **U — the isotropic translational shape-function tie** into a non-matching host
   (`-shape`/`-xi`/`-host`), host-agnostic, 2D + 3D. Frame-objective (no `-corot` needed).
2. **Stress-free birth — the `g0` offset capture** (the parent `ASDEmbeddedNodeElement`
   `m_U0` pattern), so an element added to an already-deformed host (staged construction) is
   born force- and stress-free. **DONE — shipped in PR #214** (default ON; `-absolute` opts
   out). This is the test that separates a world-class embedment from a toy.
3. **Enforcement for the U tie:** `penalty | al` (augmented Lagrangian, near-exact at moderate
   `K`) + **bipenalty** (explicit `dt_cr`, the genuine differentiator) + **`-kt auto`**
   host-stiffness conditioning. (NB the UR-only `-kr auto` and the D9-only `k_eff` bipenalty
   path ride with their *experimental* features, not the core.)

### 14.3 EXPERIMENTAL — flag-gated, NOT validated (kept in place, no migration commitment)
These remain available behind their flags but are **not part of the v1 validation battery**
and **must not be sold as validated**. They stay on the element (no deprecation/migration
plan); promotion to "supported" requires their own validation pass.
- **UP (pressure, `-pressure`)** — niche poromechanics (u-p saturated hosts). Experimental.
- **UR (rotation, `-rot`)** — ties the constrained node's rotations to the host **continuum
  SPIN** `θ = ½ curl(u) = skew(∇u)`. **This is NOT the lever-based moment transfer** Abaqus /
  LS-DYNA use for shell-to-solid coupling; on a **CST / TET4** host `∂N/∂x` is element-constant
  so it degenerates to a single **rigid element-wide spin** (UR-4). It is a documented
  *spin-restraint / regularization*, not a validated moment-transfer mechanism. **Real
  continuum moment transfer is a separate future track (its own ADR), not this curl tie.**
- **D9 (per-direction cohesive / unilateral-contact / approximate-friction / bond materials),
  incl. the v2 `-corot` frame co-rotation** — this is conceptually an **interface / contact**
  capability, **not an embedment feature**: it is adjacent to the `LadrunoContact` (scoped,
  `project_ladruno_contact`) and [[LadrunoBondSlip_guide|LadrunoBondSlip]] lineage. It is kept
  here as an **experimental convenience**; rigorous coupled friction is `LadrunoContact`'s job
  (uncoupled per-direction uniaxials only approximate Coulomb friction — a fixed `ElasticPP`
  slip, not `μ·N`). Experimental.

### 14.4 v1 must-fix correctness list (code tasks — NOT this docs PR)
Status-corrected against what has already shipped:
- **(a) `g0` stress-free capture — DONE (PR #214).** No jolt on staged addition; relative gap
  drives all traction; the D9 `setTrialStrain` sees the relative gap (materials born
  unstrained). Validated by `test_staged_*` in `tests/test_ladrunoEmbeddedNode_element.py`.
- **(b) `getInitialStiff` must use the per-direction INITIAL tangent in D9 mode.** Today
  `getInitialStiff()` aliases `getTangentStiff()` → `formTransTraction()` → `setTrialStrain()`,
  so in **D9 mode** the "initial" stiffness is **state-dependent** and **mutates material
  state during a query**. This is a real bug **for D9 only** (experimental); for the U core
  (`matMode 0`) `getInitialStiff = K_u·I` is exact and state-independent, so it does **not**
  affect the v1-core claim. Fix gates D9 promotion.
- **(c) Add a VERSION field to `sendSelf`/`recvSelf`.** The serialization has **no version
  field** (`hdr(0)=tag`, `hdr(1)=ndm`, …) yet the format has changed every phase (hdr grew to
  **29** in #214). Originally this was "add the version *before* the g0 format change" — but
  #214 already changed it, so this is now **retroactive hygiene**: add a version int now and
  accept that pre-#214 DBs are already incompatible (acceptable for this fast-moving internal
  fork — DB format is not externally stable, per `LEDGER_quirks`). Lower urgency than a
  correctness bug.

- **(d) `bipenalty` `dt_cr` host-side reduced mass — FIXED (review #2, §15).** The
  self-reported critical step lumped the penalty mass on the slave only and reported
  `2√(m_p/k_eff)`, ignoring that the penalty stiffness also loads the host DOFs; it was
  **unconservative** in `-dtcr` mode on a light/coarse host. Now tightened with the reduced
  mass `μ = m_p·M_h/(m_p+M_h)` when the host element is queryable (warn on a massless host;
  `-dtcr`-without-host is the documented user-asserts-dt contract). This is a **core** item
  and the only review finding with real analysis consequence for the supported v1.

The adversarial review **`wf_7e5f8152-dfe` (§15)** confirmed the core is sound (zero blockers,
zero core mechanics bugs) and fixed the 4 **core** hardening items (#2/#3/#4/#15). Item (b)
above (`getInitialStiff` in D9) was re-confirmed and remains the key **D9-promotion** gate;
the other experimental findings are batched per §15. **This re-scope was docs-only** (ADR +
validation plan + status markers); the core fixes landed as a follow-up code PR. See
[[27_ladruno_embedded_node_validation_plan]] for the battery, §9 for the build history, and
§15 for the full review record.

---

## 15. Adversarial review record (workflow `wf_7e5f8152-dfe`, 2026-06-07)

8 review dimensions read the element + kernel source; **one adversarial verifier per finding
tried to refute it** against the code; a lead synthesised. **27 agents, 18 findings, 17
confirmed, 1 refuted, 0 uncertain.** Verdict: **the v1 SUPPORTED CORE (U tie + `g0` birth +
penalty/AL/bipenalty) is SOUND and shippable — zero blockers, zero confirmed core *mechanics*
bugs.** The experimental modes (UR/UP/D9/`-corot`) carry real issues but every one is gated
behind a flag the core never sets (the clean `matMode==0` vs `matMode==1`/UR/UP separation is
the element's saving grace). The review independently re-confirmed the §14.4(b) `getInitialStiff`
must-fix.

| # | Sev | Scope | Finding | Status |
|---|-----|-------|---------|--------|
| **#2** | major | **core** | `bipenalty` self-reported `dt_cr` ignored the host-side reduced mass (`μ = m_p·M_h/(m_p+M_h)`) — **unconservative** in `-dtcr` mode on a light/coarse host (penalty stiffness loads the penalty-massless host DOFs) | **FIXED** — tighten with `μ` when the host is queryable (`minPosDiagonal(host->getMass())`); warn on a massless host; `-dtcr`-without-host is the documented user-asserts-dt contract |
| #3 | minor | core | `getExplicitCriticalTimeStep` could substitute the rotational `dt` when the translational class was active-but-unbounded (`m_p==0`, `k>0`), advertising a false bound | **FIXED** — return −1 (no certificate) when any active class is unbounded. *(Defensive: unreachable today — `resolveBipenalty` always zeroes `m_p` and `I_p` together via early return.)* |
| #4 | minor | core | `-k`/`-kr` followed by a non-numeric token (a forgotten value, e.g. `-k -enforce`) → `atof`→0 → silently disabled tie | **FIXED** — `strtod` end-pointer guard, reject non-numeric (≠ `auto`) |
| #15 | nit | core | `setDomain(0)` wrote through `theNodes` when null (default-constructed element, before `recvSelf` allocates) | **FIXED** — `if (theNodes != 0)` guard. *(Latent: not reached by the normal broker→recvSelf→setDomain lifecycle.)* |
| #1 | major | exp | `getInitialStiff` aliases `getTangentStiff` in D9 → returns the current (softened/zero) interface tangent **and** mutates material trial state during a const K0 query (= §14.4(b)). #7/#8 share this root cause | open — gates D9 promotion (§14.4(b)) |
| #6 | minor | exp | `-corot` frame has no activation baseline (`θ0`) — staged birth into a rotated host references absolute host spin | open — `-corot`/D9 only |
| #7 | minor | exp | `-corot` drops the `∂e_d/∂u` tangent term (documented R7) — costs quadratic convergence only | open (documented) |
| #8 | minor | exp | UP guard uses `ndf>=ndm+1` not `==`: a 3D `ndf=6` node binds a rotation DOF as "pressure" | open — UP only |
| #9 | minor | exp | `-rot` on a 2D u-p node (`ndf=3`) binds the pressure DOF as drilling rotation (count-only guard, no semantics backstop) | open — UR only |
| #5 | minor | exp | D9 `force`/`localForce` responses call `setTrialStrain` (mutate material trial state during a read) | open — D9 only |
| #10 | minor | exp | `localForce` response string shadowed by the global-traction case (D9 local case unreachable under that name) | open — D9 only |
| #11 | minor | exp | `-normal`/`-orient` silently ignored (consumed but dropped) when no `-mat*` given | open — D9 only |
| #12 | nit | exp | UP scatter assumes ndf homogeneity across coupled nodes (no ASD-style check) | open — UP only |
| #13 | nit | exp | `-corot`-without-material reports the gradient error before the clearer matMode-required error (diagnostic ordering) | open — cosmetic |
| #14 | nit | tooling | `<4 args` usage string omits the Phase-1b/2/2b options | open — cosmetic |
| — | — | — | (refuted ×1; one collateral inaccuracy in #15's "parent has this guard" rationale — the parent uses a `std::vector`, no such branch — noted, the underlying null write was still real) | — |

**Action taken (this pass):** the 4 **core** items (#2, #3, #4, #15) were FIXED + a host-reduced-mass
`dt_cr` regression test (`test_bipenalty_dtcr_host_reduced_mass`) and a `-k`-guard test added.
The **experimental** items are batched for their modes' promotion (the D9 set gated on §14.4(b)).
