---
title: "ADR 42 — LadrunoBuckle (prestressed modal + linear buckling): design spec"
project: Ladruno
type: ADR / design spec
status: draft
priority: high
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"          # cross-code synthesis — this ADR = its row "B"
  - "[[modal_gap_study/03_kratos_source]]"      # PRIMARY template: Kratos prebuckling_strategy (Jia-Mang)
  - "[[modal_gap_study/02_abaqus_theory]]"      # Abaqus *BUCKLE / linear-perturbation-about-base-state
  - "[[modal_gap_study/01_opensees_current_state]]" # ground-truth file:line audit (dead buckling path)
  - "[[46_ladruno_complex_modal_adr]]"          # SIBLING ADR 46 — complex/state-space modal (QZ on Φᵀ{M,C,K}Φ)
  - "[[43_ladruno_feast_eigensolver_adr]]"      # SIBLING ADR 43 — FEAST/Sturm robust+parallel eigensolver
  - "[[44_ladruno_frequency_domain_adr]]"       # SIBLING ADR 44 — modal FRF / SSD / random
  - "[[20_ladruno_arclength_stabilized_adr]]"              # NON-goal boundary: nonlinear post-buckling (arc-length)
  - "[[21_ladruno_dynamic_relaxation_adr]]"     # NON-goal boundary: quasi-static collapse path
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
tags: [adr, analysis, integrator, eigen, buckling, prestressed-modal, stress-stiffening, geometric-stiffness, jia-mang, arpack]
updated: 2026-06-22
---

# ADR 42 — `LadrunoBuckle` (prestressed modal + linear buckling)

> **Strategic role (load-bearing assessment — see [[modal_gap_study/00_SYNTHESIS]] §6).**
> **Standalone analysis type — modest unlock multiplier.** Valuable in its own right (stability /
> buckling factors are a missing analysis type, and prestressed modal gives correct frequencies
> under load), but little downstream builds on it. Its one genuine cross-link: the
> geometric-stiffness eigenpath could feed limit-point / bifurcation detection into the arc-length
> solvers (ADRs 20/22). **Opportunistic — build when a specific project asks.**

**Status:** draft. Design only — **no code has landed.** classTag **33021 RESERVED, not yet
built** (HANDLER/ANALYSIS band; this is a buckling *analysis/integrator*, not an element or
material). This is the design spec for two complementary capabilities that share one piece of
machinery: **(a) prestressed (stress-stiffened) modal** and **(b) linearized buckling** via the
Kratos/Jia–Mang two-tangent difference.

> [!info] What `LadrunoBuckle` is, in one line
> A small buckling **analysis/integrator** that (a) lets `eigen` run on a *preloaded* tangent so
> you get frequencies-under-load / stress-stiffened modes, and (b) snapshots **two static tangents**
> at consecutive load levels, differences them into `ΔK = K_ref − K`, and feeds the generalized
> pair **`K_ref φ = λ ΔK φ`** to the *existing* ARPACK `eigen` to extract buckling load factors and
> mode shapes — with an eigenvalue-tracking continuation loop. **No new `Kg` assembler:** each
> element's own corotational / P-Δ geometric nonlinearity already builds `K0+Kg` into its tangent.

---

## 1. Driver & goal

OpenSees today has **no first-class buckling procedure** and **no prestressed-modal recipe** that is
documented and verified. The ground-truth audit is unambiguous:

- The `eigen` parser carries a comment listing `// 2 - buckling` as a `generalizedAlgo` mode
  (`SRC/runtime/commands/analysis/analysis.cpp:265–268`), but `generalizedAlgo` is a **bool** — the
  buckling branch is **never wired**. It is **dead code**. Same in `OpenSeesCommands.cpp`.
- The element API *has* geometric-stiffness hooks — `FE_Element::addKgToTang`
  (`SRC/analysis/fe_ele/FE_Element.cpp:423–436`) → `Element::getGeometricTangentStiff` — but the
  base implementation **returns a zeroed matrix** (`SRC/element/Element.cpp:801–811`); only
  `PFEMElement2DCompressible` overrides it, and `addKgToTang` is used *only* in `PFEMIntegrator`
  (`SRC/analysis/integrator/PFEMIntegrator.cpp:403`), **never in any eigen path**.
- The eigen assembly forms **only** `Kt` (`StaticAnalysis::eigen`,
  `SRC/analysis/analysis/StaticAnalysis.cpp:285–293`, `elePtr->addKtToTang(1.0)`) and mass `M`
  (`:300–323`, `addMtoTang`). There is **no `Kg` assembly anywhere in the eigen path**.

So the only route to "buckling" today is to run a nonlinear incremental analysis and *watch the
tangent go singular* — a fragile, non-self-documenting procedure that gives no mode shape and no
load factor. And "modal under load" requires the user to know, undocumented, that corotational /
P-Δ transforms fold `Kg` into `Kt` at the current state.

**Goal.** Make two queries first-class, mirroring Abaqus `*BUCKLE` and a post-load `FREQUENCY` step:

1. **Buckling load factors + mode shapes** — `buckling N` → returns the lowest `N` critical load
   factors `λᵢ` (multipliers on the live load) and the buckling mode shapes (on the nodes, so
   existing recorders/visualizers work).
2. **Frequencies-under-load / stress-stiffened modes** — a documented, verified
   *static-step-then-`eigen`* recipe that confirms `Kg` flows into the tangent (tension stiffens,
   compression softens), exposed with a clarifying `eigen -prestressed` alias so intent is explicit.

These directly serve: slender steel members and frames (Euler / frame buckling), shell/plate
stability, guyed/cable structures and tensioned membranes (frequency rises under tension), and any
model where preload reshapes the modes (the seismic "soft-storey under gravity" check).

---

## 2. Decision summary

**Adopt the Kratos `prebuckling_strategy` (Jia–Mang consistent linearization): a two-tangent
*difference*, not the naive single-snapshot `(K0 + λ Kg) φ = 0`.**

The classical textbook linear-buckling eigenproblem is

$$(K_0 + \lambda\, K_g)\,\phi = 0,$$

where `K0` is the small-displacement stiffness and `Kg` is a **separately assembled** geometric
stiffness driven by the perturbation stress field. Implementing that form in OpenSees would require
**a real `Kg` assembler** — i.e. every element would have to override `getGeometricTangentStiff()`
to return a non-zero initial-stress matrix, and the eigen assembly would have to form `Kg` into the
"M" slot. That is a large, element-by-element effort and is exactly the dead path that returns zero
today.

The Jia–Mang form sidesteps the assembler entirely. Run **genuine nonlinear static steps**; at each
converged load level the element's *own* corotational / P-Δ geometric nonlinearity has already
folded `K0+Kg(σ)` into its tangent (verified: `PDeltaCrdTransf{2,3}d.cpp`,
`CorotCrdTransf{2,3}d.cpp` produce the current-state geometric tangent). Snapshot two tangents at
consecutive load levels and **difference** them:

$$\Delta K = K_{\text{ref}} - K,$$

then solve the generalized pair with the eigensolver OpenSees already ships:

$$\boxed{\,K_{\text{ref}}\,\phi = \lambda\,\Delta K\,\phi\,}$$

`ΔK` *is* the discrete derivative of the tangent with respect to the load multiplier — exactly the
incremental geometric stiffness, recovered for free from two real tangents. The eigenvalue `λ` is a
multiplier on the incremental stiffness change; mapped through the load increment it gives the
critical load factor.

**Why this reuses everything (the elegance):**

| Need | Naive `(K0+λKg)φ=0` | **Jia–Mang two-tangent (adopted)** |
|---|---|---|
| `Kg` assembler | **new**, per element | **none** — element tangent already has `K0+Kg` |
| Eigen solve | new "Kg-in-M" assembly | **existing** `eigen` — feed `(K_ref, ΔK)` as `(A, M)` |
| Prestress / preload | manual | **free** — `K_ref` is captured at a preloaded state |
| Accuracy near nonlinear prebuckling | poor (single far `Kg`) | **better** — consistent linearization at the load level + continuation |
| New C++ | element + eigen + assembly | **one small analysis/integrator** (`LadrunoBuckle`) |

Prestressed modal (a) is then *almost free*: it is the same "capture a preloaded tangent" act,
followed by the ordinary `K φ = ω² M φ` solve — no `ΔK` at all. The two capabilities share the
**tangent-capture** machinery; they differ only in what goes in the second matrix slot (`M` for
modal, `ΔK` for buckling).

This is ADR row **"B — Prestressed modal + buckling"** in
[[modal_gap_study/00_SYNTHESIS]]; the template is verbatim Kratos
`prebuckling_strategy.hpp` (see [[modal_gap_study/03_kratos_source]] §2).

---

## 3. Scope-fence & classTag

`#define ANALYSIS_TAG_LadrunoBuckle 33021` (HANDLER/ANALYSIS band; **RESERVED, not yet built**).
Record in `SRC/classTags.h` + `LEDGER_implementations.md` at reservation time so no sibling collides
(ADR 46 complex-modal and ADR 43 FEAST will want neighboring tags).

**In scope**

- **(a) Prestressed / stress-stiffened modal** — documented static-step-then-`eigen` recipe + an
  `eigen -prestressed` intent alias; verify `Kg` flows from corot / P-Δ.
- **(b) Linearized buckling** via the two-tangent difference `K_ref φ = λ ΔK φ`, fed to ARPACK `eigen`.
- **Eigenvalue-tracking continuation** (initial / small-probe / path-following increments;
  convergence on `|λ/λ_prev| < ratio`) — the Jia–Mang accuracy upgrade.
- Reporting buckling load factors `λᵢ`, critical loads `P_cr = P_base + λᵢ·Q`, and mode shapes on
  the nodes (reuse `setEigenvector`).

**NOT in scope (explicit boundaries)**

- **Nonlinear / imperfection-seeded post-buckling collapse** — that is the **arc-length** family
  ([[20_ladruno_arclength_stabilized_adr|ADR 20 (stabilized arc-length)]]) and quasi-static collapse via dynamic relaxation
  ([[21_ladruno_dynamic_relaxation_adr|ADR 21]]). `LadrunoBuckle` estimates the *bifurcation*; it
  does **not** trace the secondary path. (Mirrors Abaqus's `*BUCKLE` → seed → `STATIC, RIKS` split.)
- **A general `Kg` element assembler** — deliberately avoided (§2).
- **Plasticity / rate effects in the buckling estimate** — like Abaqus `*BUCKLE`, the linearized
  factor is an *elastic* critical-load estimate; material nonlinearity in the prebuckling steps
  shapes `K_ref`/`ΔK` but the eigen estimate itself ignores plastic-rate effects.
- **A new eigensolver** — reuse ARPACK. Robustness/parallel/band-targeting is ADR 43's job.

---

## 4. Formulation

### 4.1 Base-state linear perturbation (shared frame)

Both capabilities ride the linear-perturbation-about-a-base-state idea (Abaqus TG §2.1.1). The
tangent assembled at a base state splits into

$$K = \underbrace{K_{\text{mat}}}_{\text{material}}
    + \underbrace{K_{\text{geom}}(\sigma_{\text{base}})}_{\text{initial-stress / geometric}}
    + \underbrace{K_{\text{load}}}_{\text{load stiffness (follower)}} .$$

The key OpenSees fact: with a **corotational** or **P-Δ** `CrdTransf`, the element returns
`K_mat + K_geom(σ_base)` *as its ordinary tangent* at the current state — stress stiffening is
already "in the tangent." (`PDeltaCrdTransf{2,3}d.cpp`, `CorotCrdTransf{2,3}d.cpp`.) No extra
assembly is needed to *use* it; the only question is *capturing* it at the right state (§6).

### 4.2 Prestressed modal (capability a)

After a (possibly nonlinear) static preload step to state `σ_base`, solve the symmetric generalized
eigenproblem on the **preloaded tangent**:

$$\big(K(\sigma_{\text{base}}) - \omega^2 M\big)\,\phi = 0,
  \qquad K = K_{\text{mat}} + K_{\text{geom}}(\sigma_{\text{base}}).$$

Because `K_geom` is present, `K` **need not be positive definite** — a compressed member softens
(`ω` drops, `ω²` can even go negative near instability), a tensioned cable stiffens (`ω` rises).
This is the guitar-string / tensioned-cable physics and the seismic "modes under gravity" check.
**No new math** — it is the existing `eigen` on a tangent captured *after* a load step instead of at
the undeformed state.

### 4.3 Linearized buckling (capability b) — the two-tangent difference

**Classical form (for contrast).** With `K0` the base-state stiffness and `Kg` the initial-stress
stiffness of the linear stress response to a live load `Q`, the critical multipliers solve

$$(K_0 + \lambda\, K_\Delta)\,\phi = 0,
  \qquad P_{cr} = P_{\text{base}} + \lambda\, Q. \tag{classical}$$

This needs a *separately assembled* `K_Δ = K_g(\Delta\sigma(Q))` — the assembler we refuse to build.

**Adopted Jia–Mang form.** Drive two genuine nonlinear static steps and capture the converged
tangents `K_ref` (reference, lower load) and `K` (advanced load). Difference them:

$$\Delta K = K_{\text{ref}} - K \;\approx\; -\frac{dK}{d\lambda}\,\Delta\lambda . \tag{difference}$$

Solve the generalized pair

$$K_{\text{ref}}\,\phi = \lambda^\star\,\Delta K\,\phi, \tag{eigenproblem}$$

and map the eigenvalue to the critical load factor through the load increment used between the two
snapshots:

$$\lambda_{\text{crit}} = \lambda_{\text{prev}} + \lambda^\star \cdot \Delta\lambda_{\text{load}},
  \qquad P_{cr} = \lambda_{\text{crit}}\, P_{\text{ref live}} . \tag{critical}$$

The smallest `λ★` (the mode whose incremental stiffness vanishes first) gives the governing buckling
load; `φ` is the buckling mode shape.

**Why differencing two real tangents is equivalent and assembler-free.** Write the load-dependent
tangent as `K(λ) = K_mat + λ Kg₁ + O(λ²)` near the base state (the geometric term is, to leading
order, linear in the stress and hence in the load multiplier). Then between two levels separated by
`Δλ_load`,

$$\Delta K = K(\lambda_{\text{ref}}) - K(\lambda_{\text{ref}}+\Delta\lambda_{\text{load}})
          \;=\; -\,\Delta\lambda_{\text{load}}\; K_{g1} + O(\Delta\lambda^2),$$

so `ΔK` *is* (minus, scaled) the incremental geometric stiffness `Kg` — recovered numerically from
two tangents the elements already produce, rather than assembled analytically. Substituting,
`K_ref φ = λ★ ΔK φ` is the discretized `(K0 + λ Kg) φ = 0` with `λ = −1/(λ★ Δλ_load)` to leading
order; the continuation loop (§4.4) removes the `O(Δλ²)` error by re-linearizing at the converging
load. The sign of `ΔK` and the branch chosen for `λ★` are handled in the continuation/extraction
(see Risk R1, §9): take the orientation `K_ref − K` so that a *destabilizing* (compressive) live
load yields a positive critical factor.

### 4.4 Eigenvalue-tracking continuation loop

A single difference at one arbitrary pair of levels is the classic linear-buckling estimate (v1,
§7). Jia–Mang's accuracy upgrade (v2) is an **eigenvalue-tracking continuation** that re-linearizes
at successively better load levels until the smallest eigenvalue stabilizes (verbatim structure from
Kratos `prebuckling_strategy.hpp`):

1. **Iteration 1 — initial increment.** Apply `initial_load_increment · (λ★ + λ_prev)`; solve the
   pair; set `λ = λ★ · Δλ_load`.
2. **Odd iterations — probe.** Advance a *tiny* `small_load_increment · λ_prev`, re-solve the eig,
   and test convergence `|λ / λ_prev| < convergence_ratio`.
3. **Even iterations — path-following.** Advance `path_following_step · λ` toward the critical load.
4. Loop until the smallest eigenvalue stabilizes ⇒ that is the critical buckling factor. Defaults
   (Kratos-proven): `initial_load_increment = 1.0`, `small_load_increment = 5e-4`,
   `path_following_step = 0.5`, `convergence_ratio = 0.05`.

Each iteration is: one nonlinear static solve (advance load), capture tangent, difference against the
stored reference, one ARPACK generalized solve. The outer loop is small (handful of iterations), not
hundreds of restarts.

### 4.5 Contrast summary

| | classical `(K0+λKg)φ=0` | **adopted `K_ref φ = λ ΔK φ`** |
|---|---|---|
| `Kg` source | analytic per-element assembler (must build) | numerical difference of two real tangents |
| Prebuckling nonlinearity | ignored (single far `Kg`) | captured (real nonlinear steps + re-linearization) |
| Implementation | element + assembly + eigen | tangent capture + existing eigen + small loop |
| Equivalence | — | `ΔK ≈ −Δλ·Kg` ⇒ same eigenproblem to leading order; continuation removes higher order |

---

## 5. Public API

### 5.1 Buckling — Tcl

```tcl
# ... build model, apply DEAD/base load in one pattern, apply LIVE/perturbation load in another ...
constraints Transformation       ;# elimination-style (Kratos forces this; see R-notes)
numberer RCM
system BandGeneral               ;# any LinearSOE; ArpackSOE wraps it for the shift-invert
integrator LoadControl 1.0
algorithm Newton
analysis Static

buckling 5 \
    -liveLoadPattern   2 \
    -refLoadFactor     1.0 \
    -loadIncrement     0.10 \
    -continuation      1 \
    -initialInc        1.0 -smallInc 5.0e-4 -pathStep 0.5 -ratio 0.05 \
    -solver            genBandArpack
# returns: TCL list of the N lowest buckling load factors {λ1 λ2 ... λN}
```

Outputs:
- **Return value** — the `N` lowest buckling load factors `λᵢ` (smallest = governing).
- `set lam [lindex [buckling 5 ...] 0]` — governing factor.
- **Mode shapes** — written to the nodes via `setEigenvector`, so
  `nodeEigenvector $node $mode` and any eigen recorder / MPCO visualizer works unchanged.
- **Critical loads** — `bucklingResponse -criticalLoads` → `{P_cr,1 ... P_cr,N}` with
  `P_cr,i = λ_i · (live-load resultant)`; or compute in script from `λᵢ` and the applied live load.

### 5.2 Buckling — openseespy

```python
import openseespy.opensees as ops
# ... model, dead-load pattern 1, live-load pattern 2, Static analysis set up ...
factors = ops.buckling(5,
                       '-liveLoadPattern', 2,
                       '-refLoadFactor', 1.0,
                       '-loadIncrement', 0.10,
                       '-continuation', 1,
                       '-initialInc', 1.0, '-smallInc', 5.0e-4,
                       '-pathStep', 0.5, '-ratio', 0.05,
                       '-solver', 'genBandArpack')   # -> [λ1, λ2, ..., λ5]
phi1 = [ops.nodeEigenvector(n, 1) for n in ops.getNodeTags()]   # mode-1 shape
Pcr  = ops.bucklingResponse('-criticalLoads')                   # -> [Pcr1, ..., Pcr5]
```

### 5.3 Prestressed modal — recipe + intent alias

No new command strictly required — it is **static-step-then-`eigen`** — but ship an explicit alias so
intent is documented and the `Kg`-in-tangent contract is asserted:

```tcl
# 1. preload to the base state (gravity / tension / axial load)
pattern Plain 1 Linear { load ... }
system BandGeneral; numberer RCM; constraints Transformation
integrator LoadControl 1.0; algorithm Newton; analysis Static
analyze 1                       ;# tangent now carries Kg(σ_base)

# 2. modal on the preloaded tangent
set lambdas [eigen -prestressed 6]    ;# alias for "eigen 6 on the CURRENT tangent, do NOT rebuild K0"
#   plain `eigen 6` also works IF the active integrator left the preloaded tangent in place;
#   `-prestressed` makes the no-rebuild contract explicit and verifies a corot/PDelta transform is present.
foreach lam $lambdas { puts "f = [expr sqrt($lam)/(2.0*acos(-1.0))] Hz" }
```

```python
ops.analyze(1)                                  # preload
lambdas = ops.eigen('-prestressed', 6)          # frequencies under load
import math
freqs = [math.sqrt(l)/(2*math.pi) for l in lambdas]
```

`eigen -prestressed N`: assemble the eigen `A` matrix from the **current** element tangents (which
already contain `Kg`) instead of zeroing and reforming `K0` at the undeformed state; warn if no
element uses a geometric `CrdTransf` (then "prestressed" is a no-op and the user should be told).

---

## 6. OpenSees integration points (file:line)

1. **Tangent capture around a static step — the crux.** `StaticAnalysis::eigen`
   (`SRC/analysis/analysis/StaticAnalysis.cpp:244–348`) currently *always* rebuilds `A` from
   `addKtToTang(1.0)` at the current state (`:285–293`) and `M` from `addMtoTang` (`:300–323`).
   `LadrunoBuckle` needs to **snapshot the assembled tangent at two load levels** rather than
   reform it once. Two clean options:
   - capture the `LinearSOE`'s assembled `A` immediately after each converged static step's
     `formTangent` (the analysis already holds `K_ref` and `K` as SOE matrices), or
   - re-run the eigen-style assembly loop (`:285–293`) at each level into a stored `Matrix`.
   The captured matrices are then handed to the EigenSOE as the generalized pair.

2. **Feeding `(K_ref, ΔK)` to the EigenSOE.** The eigen pair is consumed via `theEigenSOE->addA(...)`
   (the `K`/`A` slot, `StaticAnalysis.cpp:288`) and `theEigenSOE->addM(...)` (the `M` slot, `:306` /
   `:318`). For buckling we put **`K_ref` into `addA`** and **`ΔK` into `addM`**, then call
   `theEigenSOE->solve(numMode, /*generalized=*/true, /*findSmallest=*/true)` (`:330`) — ARPACK's
   shift-invert generalized solver (`ArpackSOE`/`ArpackSolver`, `bmat='G'`, mode 3) handles
   `A φ = λ M φ` with an indefinite `M=ΔK` (it does not require `M` SPD, only the shifted
   `(A−σM)` factorable). Mode shapes go back to the nodes through `setEigenvector` (`:343`),
   exactly as modal does — so recorders are free.

3. **The dead "buckling" flag.** `analysis.cpp:265–268` (`// 2 - buckling` on a `bool
   generalizedAlgo`) and its `OpenSeesCommands.cpp` twin: either **wire** this to `LadrunoBuckle`
   or **remove** the misleading comment. Decision: remove the dead comment/branch and register a
   real `buckling` command (cleaner than overloading `eigen`). If the removal touches the upstream
   parser, log it in `LEDGER_vanilla_files.md` with a `// Ladruno` marker.

4. **Geometric stiffness origin (verify, don't build).** `Element::getGeometricTangentStiff`
   returns zero (`SRC/element/Element.cpp:801–811`) and `FE_Element::addKgToTang`
   (`SRC/analysis/fe_ele/FE_Element.cpp:423–436`) is unused by eigen — confirming we must **not**
   rely on it. The geometric term lives in the corot / P-Δ transforms
   (`SRC/coordTransformation/{PDelta,Corot}CrdTransf{2,3}d.cpp`), folded into the ordinary element
   tangent. P1 must *verify* this flow (test that a tensioned member's `eigen` frequency rises).

5. **Command registration.** New `buckling` / `bucklingResponse` and the `eigen -prestressed` alias:
   register in the xara runtime (`SRC/runtime/commands/analysis/analysis.cpp`, near the `eigen`
   registration at `:250`) and the legacy/parallel + Python paths
   (`SRC/interpreter/OpenSeesCommands.cpp:270` region), mirroring how `eigen` is dual-registered.

6. **Analysis/integrator class.** New `SRC/analysis/integrator/LadrunoBuckle.{h,cpp}` (or a thin
   `SRC/analysis/analysis/LadrunoBucklingAnalysis.*`) that orchestrates: static load step → tangent
   capture → difference → EigenSOE pair → continuation loop. Model it on `StaticAnalysis::eigen`
   plus `EigenIntegrator` (`SRC/analysis/integrator/EigenIntegrator.cpp:191–196`).

---

## 7. Phased roadmap + gates

| Phase | Deliverable | Gate (must pass to proceed) |
|---|---|---|
| **P1 — Prestressed modal (recipe + verify Kg flows)** | Document static-step-then-`eigen`; add `eigen -prestressed` alias + "no geometric transform" warning. **No buckling code yet.** | **G1:** tensioned cable/string `eigen` frequency *rises* with tension and matches `f = (n/2L)√(T/μ)` (corot transform) within tol; compressed column frequency *drops* toward zero as `P → P_cr`. Confirms `Kg` is in the captured tangent. |
| **P2 — Two-tangent buckling + continuation** | `LadrunoBuckle` analysis/integrator: capture two static tangents, difference, feed `(K_ref, ΔK)` to ARPACK `eigen`, extract `λᵢ`/modes; v1 single-difference, then v2 eigenvalue-tracking continuation loop; `buckling` + `bucklingResponse` commands; remove dead buckling path. | **G2:** simply-supported Euler column → `λ·P_applied ≈ π²EI/L²` within tol; mode 1 = half-sine. v1 vs v2 continuation both converge; v2 tighter. Sign/branch of `ΔK` correct (positive critical factor for compressive live load). |
| **P3 — Validation battery + docs** | Full oracle battery (§8); user-guide section; ledger/banner/header obligations. | **G3:** all §8 cases within tol vs analytic; buckling `λ·P` cross-checks a fine arc-length collapse load within a few %; prestressed frequencies match closed form. Adversarial gate per fork policy (novel-ish analysis core ⇒ run it). |

Loop discipline: phased, each gate adversarially reviewed where the fork policy calls for it (P2/P3
touch the analysis core and are novel-ish ⇒ gate; P1 is a documented recipe + thin alias ⇒ lighter).

---

## 8. Validation / oracle battery

**Linear buckling — Euler column, multiple boundary conditions.** Critical load

$$P_{cr} = \frac{\pi^2 E I}{(K L)^2},$$

with effective-length factor `K`: pinned–pinned `K=1`, fixed–free (cantilever) `K=2`,
fixed–fixed `K=0.5`, fixed–pinned `K≈0.699`. Run `buckling` with a unit axial live load; assert
`λ·P_applied = P_cr` for each BC; assert mode-1 shape (half-sine / quarter-cosine etc.).

**Plate buckling — simply-supported rectangular plate, uniaxial compression.**

$$N_{cr} = k\,\frac{\pi^2 D}{b^2}, \qquad D = \frac{E t^3}{12(1-\nu^2)},$$

with buckling coefficient `k` (e.g. `k=4` for a long SS plate under uniform uniaxial compression).
Mesh a shell plate; assert `λ·N_applied = N_cr` and the checkerboard/half-wave mode.

**Frame buckling.** A portal frame or a fixed-base column with a beam — compare `λ·P` to a known
frame-buckling reference (alignment-chart / effective-length result, or a published value). Sanity
that the procedure handles assembled multi-member `Kg`.

**Prestressed modal — string / cable frequency shift under tension.** Taut string transverse modes

$$f_n = \frac{n}{2L}\sqrt{\frac{T}{\mu}},$$

`μ` = mass per length, `T` = tension, `n = 1,2,...` (the guitar-string analytic). Preload a corot
cable to tension `T`, run `eigen -prestressed`, assert each `f_n`. Vary `T` and check the
`f ∝ √T` scaling — the cleanest demonstration that `Kg` enters the modal tangent. Add the
compressed-column companion: frequency drops to zero as `P → P_cr` (ties buckling and prestressed
modal to the *same* `Kg`).

**Cross-checks.**
- **vs analytic** — every case above has a closed form (primary acceptance).
- **vs fine arc-length collapse** — run an imperfection-seeded arc-length ([[20_ladruno_arclength_stabilized_adr|ADR 20 (stabilized arc-length)]])
  to the limit/bifurcation load on the Euler column and the plate; the linear-buckling
  `λ·P` should bound/approximate the nonlinear collapse load (within a few % for near-linear
  prebuckling; the *gap* itself documents the linear-buckling validity envelope).

---

## 9. Risk register

> [!question] **R1 — Sign / branch of `ΔK`.**
> `ΔK = K_ref − K` vs `K − K_ref` flips the sign of `λ★`, and the *physically meaningful* branch is
> the smallest-magnitude positive critical factor. **Mitigation:** fix the orientation so a
> destabilizing (compressive) live load gives a positive `λ_crit`; in P2/G2 assert the sign against
> the Euler analytic; in the continuation, track the branch by eigenvalue continuity, not by raw
> index.

> [!question] **R2 — Choosing the two load levels.**
> Too-close levels ⇒ `ΔK ≈ 0`, ill-conditioned pair, noisy `λ`. Too-far levels ⇒ the `O(Δλ²)`
> linearization error pollutes the estimate. **Mitigation:** sensible `-loadIncrement` default
> (~10% of an estimated critical), the continuation loop's `small_load_increment` probe for the
> final tighten, and a conditioning warning when `‖ΔK‖/‖K_ref‖` is below a floor.

> [!question] **R3 — Conditioning when `ΔK` is near-singular.**
> `ΔK` (an incremental stiffness) is generally **indefinite and rank-deficient** (only the
> load-affected DOFs change), so the generalized pair `K_ref φ = λ ΔK φ` has many infinite/garbage
> eigenvalues. **Mitigation:** ARPACK shift-invert with `findSmallest` targets the few finite small
> `λ`; the shift `σ` must be chosen so `(K_ref − σ ΔK)` is factorable. Document that buckling
> *needs* a non-trivial shift; expose `-shift` (the `shift` plumbing exists in `ArpackSOE` but the
> `eigen` command hard-zeros it at `analysis.cpp:271`). Surface spurious-eigenvalue filtering
> (drop `λ` whose residual is large or whose magnitude is implausible).

> [!question] **R4 — Dead-code path cleanup.**
> The `// 2 - buckling` comment on a `bool` (`analysis.cpp:265–268`) is misleading and the
> `addKgToTang`/`getGeometricTangentStiff` hooks return zero. **Mitigation:** remove the dead
> comment/branch (log in `LEDGER_vanilla_files.md`); decide *not* to revive `addKgToTang` (we use the
> tangent-difference, so the zero-returning hook stays unused — note this in `LEDGER_quirks.md` so a
> future agent doesn't "fix" it expecting buckling to start working).

> [!question] **R5 — corot vs P-Δ `Kg` fidelity.**
> P-Δ transforms carry only the *leading* geometric term (good for frames/columns, exact for the
> classical P-Δ string-stiffness); corotational carries the full large-rotation geometric tangent
> (needed for shells/cables/snap-through). The buckling factor's accuracy inherits the transform's
> fidelity. **Mitigation:** validate Euler with **both** transforms (G2); document that
> shell/plate/cable buckling should use **corotational**; warn if a buckling run finds no geometric
> transform at all (then `ΔK` is ~0 and the result is meaningless).

> [!question] **R6 — Constraint handler interaction.**
> Kratos forces an **elimination** builder for prebuckling (penalty/Lagrange perturb the spectrum
> with huge/zero diagonal entries that pollute the small buckling eigenvalues). **Mitigation:**
> recommend `constraints Transformation` (elimination-style) for `buckling`; warn (don't hard-fail)
> on `Penalty`/`Lagrange`; mirror the `mass/stiffness_matrix_diagonal_value` regularization trick
> on constrained DOFs if needed.

> [!question] **R7 — Tangent staleness / capture point.**
> The captured tangent must be the **converged** tangent at each load level (a mid-iteration tangent
> is wrong). **Mitigation:** capture only after the static step reports convergence;
> assert `hasDomainChanged` consistency between the two snapshots (same DOF numbering / sparsity, as
> the generalized pair requires identical patterns for `K_ref` and `ΔK`).

---

## 10. Ledger / header / PR obligations

- **`LEDGER_implementations.md`** — add a row at *reservation*: feature `LadrunoBuckle (prestressed
  modal + linear buckling)`, kind `analysis/integrator`, classTag **33021**, files
  `SRC/analysis/integrator/LadrunoBuckle.*` (+ command registrations), status
  `RESERVED → in-progress` per phase, PR `#tbd`. Update status as P1→P3 land.
- **`SRC/classTags.h`** — `#define ANALYSIS_TAG_LadrunoBuckle 33021` with a `// Ladruno ADR 42`
  marker; reserve it in the same PR that creates the class so no sibling (ADR 46/43) collides.
- **Banner** — add a line to `Ladruno_scripts/banner_features.txt` (e.g.
  `Linear buckling + prestressed modal (LadrunoBuckle)`), then run
  `python Ladruno_scripts/patch_banner.py` and rebuild — **do not hand-edit** the C strings. Every
  `shipped` ledger row needs a matching banner line.
- **Header stamp** — run `Ladruno_scripts/stamp_headers.py` on the new `LadrunoBuckle.{h,cpp}` (add
  them to its GLOBS) so they carry the four-author LADRUNO header. `--check` in CI.
- **`LEDGER_vanilla_files.md`** — if the dead `// 2 - buckling` comment/branch is removed from
  `analysis.cpp` / `OpenSeesCommands.cpp`, add a row (file, why = "removed dead/unwired buckling
  flag, replaced by `buckling` command", PR) and mark the edit with a `// Ladruno` comment.
- **`LEDGER_quirks.md`** — record: "buckling uses the **tangent-difference** (`K_ref − K`), *not*
  the zero-returning `addKgToTang`/`getGeometricTangentStiff` hooks — leave those unused"; and
  "ARPACK `eigen` shift is hard-zeroed at `analysis.cpp:271`; buckling needs a real `-shift`."
- **PR base** — base on **`ladruno`** (default branch), one logical PR per phase; verify the prior
  PR shows `state == MERGED` before stacking the next (fork auto-merges fast). Use `--base ladruno`
  on every PR (never hand-stack with `--base <prev-branch>`).

---

## 11. Cross-references

- [[modal_gap_study/00_SYNTHESIS]] — this ADR is row **"B — Prestressed modal + buckling"** (build
  order A→B→C→D; B is the cheapest high-value add).
- [[modal_gap_study/03_kratos_source]] §2 — the `prebuckling_strategy.hpp` template (two-tangent
  difference + continuation defaults) ported here near-verbatim.
- [[modal_gap_study/02_abaqus_theory]] GAP 2 — `*BUCKLE` `(K0+λK_Δ)φ=0` and the
  linear-perturbation-about-base-state frame (the conceptual parent of §4).
- [[46_ladruno_complex_modal_adr|ADR 46]] (complex/state-space modal) and
  [[44_ladruno_frequency_domain_adr|ADR 44]] (modal FRF / SSD / random) — sibling modal ADRs that
  also build on `eigen`; ADR 46's `Φ` basis can be the *preloaded* basis from this ADR (prestressed
  complex modes).
- [[43_ladruno_feast_eigensolver_adr|ADR 43]] (FEAST / Sturm robust + parallel eigensolver) —
  buckling benefits directly: FEAST contour targeting + the Sturm-sequence negative-pivot count make
  the indefinite `ΔK` buckling solve **robust** (certify no critical mode is missed) and parallel.
  When ADR 43 lands, `buckling` should be re-pointable at the FEAST backend with no formulation
  change (it is still `A φ = λ M φ` with `(A,M) = (K_ref, ΔK)`).
- [[20_ladruno_arclength_stabilized_adr|ADR 20 (stabilized arc-length)]] / [[21_ladruno_dynamic_relaxation_adr|ADR 21]] — the
  **non-goal** boundary: nonlinear/imperfection post-buckling collapse. `LadrunoBuckle` estimates
  the bifurcation; arc-length traces the secondary path. The two compose (buckling mode → arc-length
  imperfection seed), exactly the Abaqus `*BUCKLE` → `STATIC, RIKS` workflow.
