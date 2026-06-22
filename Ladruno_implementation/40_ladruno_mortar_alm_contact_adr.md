---
title: LadrunoMortar — variationally-consistent mortar contact with Augmented-Lagrangian (Uzawa) enforcement + a shared FrictionalLaw kernel
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - contact
  - mortar
  - augmented-lagrangian
  - friction
  - kernel
---

# LadrunoMortar — mortar/ALM contact + FrictionalLaw

## What

ADR-40 is the **implicit, accuracy-first** complement to ADR-39's explicit, robust-first
node-to-segment (NTS) penalty `ContactDomain`. It ships **two coupled deliverables** as
**additive leaf code inside the same ADR-39 ContactDomain subsystem** (an *extension*, not a
sibling subsystem):

1. **A segment-to-segment (mortar) contact formulation** that is variationally consistent —
   the contact virtual work is *integrated* over the slave surface rather than collocated at a
   single node, so it transmits a **constant interface pressure** across a non-matching mesh
   (the contact patch test ADR-39 NTS provably fails, Risk **Q-PATCH** in ADR-39). Enforcement
   is **Augmented Lagrangian via per-Gauss-point Uzawa over a penalty kernel** — **ZERO new
   global DOFs** — recovering the exact constraint at *finite* penalty (resolves ADR-39's
   deferred Risk **Q-AL**). Frictionless first, then frictional (Coulomb + Tresca).

2. **A pluggable, header-only, OpenSees-free `FrictionalLaw` kernel** (Coulomb + Tresca, with
   consistent — including the non-symmetric pressure-coupling — derivatives) that serves **both**
   the ADR-39 NTS-penalty path *and* this mortar path from one numpy-oracle-verified return map.

**Committed scope (this ADR):** `P0`→`P4` = the two kernels (friction + projection), the
mortar narrow phase with **overlap clipping** (so the patch-test claim is honest), frictionless
Uzawa ALM (the MVP), and frictional mortar. **Hard-deferred to a successor ADR-42:** dual
(biorthogonal) Lagrange shape functions, true-LM/saddle-point enforcement, self-contact, and
NTN/NTS-via-mortar-weights. See *Risks* for the rejection rationale on each.

**Critical sequencing reality (every reviewer flagged this):** ADR-39 has shipped **P1 only** —
`LadrunoContactFE` is an *empty-connectivity zero adapter* (`FE_Element(tag,0,0)`,
`getResidual`/`getTangent` return `Zero()`), `LadrunoContactDomain::commit()` is a no-op
counter, and there is **no broad phase, no narrow phase, no projection kernel, and no friction
class** in the tree. ADR-40 therefore does **not** "dock onto a working penalty engine"; it
**co-owns and builds** the shared narrow-phase machinery (projection kernel, real adapter
connectivity, friction kernel) that ADR-39 P2/P3 will also consume. This is stated as a hard
dependency, not assumed away.

## Why

- **The accuracy bar.** ADR-39 itself documents (`39_..._adr.md` §What, Risk Q-PATCH) that NTS
  "does not pass the contact patch test … interface-stress fidelity awaits mortar (v3)." The
  Kratos `ContactStructuralMechanicsApplication` reference sets exactly this bar: mortar
  conditions with ALM and penalty, frictionless + frictional, passing the patch test. ADR-40
  delivers that bar earlier and at lower cost than full dual-mortar by using a **clipped
  Gauss-point mortar** + **Uzawa ALM**.
- **Exact constraint at finite penalty.** Pure penalty (ADR-39) leaves residual penetration
  `g = P/k_n`. ALM/Uzawa drives `g → 0` at *finite* `k_n` over an outer augmentation loop,
  resolving ADR-39 Risk Q-AL.
- **One friction law, two consumers.** A return map is the single most reused, most
  bug-prone, most oracle-testable piece of any contact code. Factoring it header-only (mirroring
  the proven `LadrunoJ2Kernel.h` discipline) and shipping it **first** de-risks both the ADR-39
  NTS path (which will *adopt* it for its planned P3 friction) and this mortar path.
- **Fork-convention fit.** Additive leaf classes, header-only OpenSees-free kernels, a
  selector-style command (`-formulation mortar -enforce alm`), and a phased loop with
  adversarial gates — exactly how every prior Ladruno element/material landed.

## Where

### New code

| Path | Role |
|---|---|
| `SRC/material/nD/LadrunoFrictionKernel.h` | **Header-only, OpenSees-free** Coulomb+Tresca tangential radial-return + normal KKT law; consistent `C_ss` and the **non-symmetric** `C_sl` block. Sits beside `LadrunoJ2Kernel.h` / `LadrunoConcrete3DKernel.h` / `LadrunoRCKernel.h` (one-line CMake add). **Deliverable (2).** |
| `SRC/domain/contact/LadrunoContactProjection.h` | **Header-only, plain-`double` only** closest-point projection (re-derived as pure functions from `SimpleContact3D::project/GetPoint/UpdateBase`) → `{xi_bar[2], gap_N, n[3], t1[3], t2[3], surface-metric g[2][2], phi_master[4]}`. **Bounded** Newton (≤10 iters + `|detK| < eps·‖g1‖‖g2‖` guard + non-convergence sentinel) — carries forward ADR-39 P2b's `contact_p2_handoff.md` BLOCKER. |
| `SRC/domain/contact/LadrunoMortarKernel.h` | **Header-only, plain-`double` only** mortar segment integration: slave/master **overlap polygon clip** on the projected master plane → sub-triangle Gauss rule → `D`, `M`, weighted gap `g̃`, and the consistent linearization `dD,dM,dn,dξ` terms. numpy-oracle tested build-free. |
| `SRC/domain/contact/LadrunoMortarPair.{h,cpp}` | **Domain-side path state** for one slave-facet ↔ master pairing: per slave-GP `{committed elastic slip s_e_n, λ_N, λ_T, slipFlag, last master-segment id + frame}`; `commit()`/`revertToLastCommit()`. The stateless-view target for the FE adapter (the ADR-40 analogue of ADR-39's planned `LadrunoContactPair`). |
| `SRC/domain/contact/LadrunoMortarSegment.{h,cpp}` | Narrow-phase formulation object owned by `LadrunoContactDomain`: drives the GP loop, calls `LadrunoContactProjection` + `LadrunoMortarKernel` + `LadrunoFrictionKernel`, returns `F_c` and `K_c`. Mortar sibling of the (planned ADR-39) NTS narrow phase. |
| `SRC/analysis/algorithm/equiSolnAlgo/LadrunoAugmentedNewton.{h,cpp}` | **The Uzawa outer-loop driver.** A `NewtonRaphson` subclass whose `solveCurrentStep()` wraps the inner Newton, reads `‖g̃‖_∞` from the `LadrunoContactDomain`, performs the multiplier update on a **frozen active set**, and re-solves until `‖g̃‖_∞ < augTol` or `maxAug` — committing **once** per physical step. (See *Integration points* for why this, not a Tcl `analyze`-loop, is mandatory.) |
| `tests/contact/ladruno_friction_reference.py` | numpy oracle for `LadrunoFrictionKernel.h` (stick/slip/cone return + FD-check of `C_ss`,`C_sl`). |
| `tests/contact/ladruno_mortar_reference.py` | numpy oracle for projection + clipped `D`/`M` (partition-of-unity + constant-pressure patch). |

### Modify vanilla — **NONE new beyond ADR-39's already-ledgered footprint.**

- `SRC/Domain.{h,cpp}` — **no new edit.** ADR-39 already added the `LadrunoContactDomain*`
  member + `commit()`/`revertToLastCommit()` hooks; ALM rides those exact seams.
- `SRC/classTags.h` — **one-line addition only if needed (see below).** `LadrunoContactFE`
  and the mortar adapter are `FE_Element`s injected by the handler with a **runtime** tag
  (`new LadrunoContactFE(numFe++)`); they are *not* `FEM_ObjectBroker`-registered, so **no new
  `ELE_TAG` is required for the MVP.** A new `equiSolnAlgo` *does* need a broker tag for
  Tcl construction → reserve **`EquiALGORITHM_TAGS_LadrunoAugmentedNewton = 33001`** in the
  algorithm registry (next free in the Ladruno band; per-registry, no collision).
- The new files are **Ladruno-authored** (ledgered in `LEDGER_implementations.md`, not
  `LEDGER_vanilla_files.md`), plus their `CMakeLists.txt` lines.
- `SRC/analysis/handler/LadrunoContactHandler.cpp` — **extend** (already a Ladruno file, not
  vanilla): emit the mortar adapter when a contact declares `-formulation mortar`.

### classTags (verified against the live `SRC/classTags.h`)

| Tag | Decision |
|---|---|
| **ELE band** | Used 33000–33014; **33009/33010 reserved** (VEM/SBFEM, ADR-26). The brief's claim *"LadrunoContactFE took ELE_TAG 33015"* is **false** — no `ELE_TAG_LadrunoContactFE` exists; the adapter carries no broker tag. If a mortar adapter is *ever* broker-serialized (parallel/DB), the next genuinely-free ELE tag is **33015**; reserve it then, not now. |
| **HANDLER band** | `HANDLER_TAG_LadrunoContactHandler = 33002` is **reused** (formulation is a per-contact selector, not a new handler). Next free = 33003 (unused). |
| **Algorithm band** | New `LadrunoAugmentedNewton` → `EquiALGORITHM_TAGS_LadrunoAugmentedNewton = 33001` (Ladruno band; the projection handler took the analogous 33001 in the *handler* registry — different registry, not a collision). |
| **FRN band** | `LadrunoFrictionKernel` is a header-only kernel, **NOT** a `FrictionModel` subclass — it does **not** consume an `FRN_TAG` (those are for bearing/isolator `frictionModel`s). Explicitly disjoint. |

### Reference (re-derive, do **not** lift)

- `SRC/element/UWelements/SimpleContact3D.{h,cpp}` — `project()` (Newton on `R=(d·g1,d·g2)`
  with the curvature term `K2=-(d·(x1-x2+x3-x4))·0.25`), `GetPoint()` (bilinear `N`),
  `UpdateBase()` (`g1,g2`), normal `n=g1×g2/‖·‖`. ~110/1270 lines re-derived as pure functions.
  *Note its unbounded `while` loop is a defect — we bound it.*
- `SRC/material/nD/UWmaterials/ContactMaterial3D.{h,cpp}` — **verified** consistent tangent:
  slip `C_ss = stiffness·scale·(g − R⊗R)` with `R = g·r`, `scale = 1 − γ/‖s_e‖`; the
  **non-symmetric** `C_sl = R·frictionCoeff`; stick `C_ss = stiffness·g`, `C_sl = 0`
  (`ContactMaterial3D.cpp` getTangent, lines ~303–339). Re-derive into the kernel, dropping its
  defects: `static int mFrictFlag` (process-global) and the stubbed `getInitialTangent`.
- `SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedRebar.cpp` — **the in-fork ALM precedent**:
  `commitState()` (lines ~297–326) stores `lambda` on the object, does the Uzawa update
  `lambda(k) += kt·gt(k)` at the converged state (frozen-within-solve so the tangent stays
  exact), and **re-projects `lambda` onto the current frame** under corotation. ADR-40 clones
  this pattern *per slave-GP* on the `LadrunoMortarPair`.
- `SRC/analysis/handler/Lagrange*` — `LagrangeDOF_Group`/`LagrangeMP_FE` — read to **confirm
  the rejected** true-LM path: `LagrangeDOF_Group` has **only** constraint-bound ctors
  (`(int,SP_Constraint&)`, `(int,MP_Constraint&)`, `(int,EQ_Constraint&)`); there is **no**
  node-less `(tag,count)` ctor, so it cannot be reused verbatim for free traction DOFs. DOF
  count freezes at `handle()` → active-set churn forces a full re-number. (See Risk Q-DOF.)
- ADR-39 (`39_ladruno_contact_domain_adr.md`) — the parent subsystem; ADR-20
  (`20_ladruno_embedded_reinforcement_adr.md`) — the Uzawa/AL §10.4 precedent.

## How

### Public API (Tcl + Python)

Extends the (ADR-39-owned) `contact`/`contactSurface` command family with a `-formulation`
selector and ALM options. `nts` (ADR-39) and `mortar` (ADR-40) are siblings on one command.

```tcl
# constraint handler (ADR-39, reused unchanged)
constraints LadrunoContact

# define meshed surfaces (ADR-39 LadrunoContactSurface; mortar needs the slave meshed
# as SEGMENTS, not just nodes -> -kind slaveSegments)
contactSurface 1 -kind masterSegments -faces $masterFaceSet
contactSurface 2 -kind slaveSegments  -faces $slaveFaceSet

# NTS penalty (ADR-39)
contact 10 -master 1 -slave 2 -formulation nts   -kn auto -kt auto -mu 0.3

# MORTAR + ALM (ADR-40) — frictionless
contact 11 -master 1 -slave 2 -formulation mortar -enforce alm \
        -kn auto -epsN auto -augTol 1e-8 -maxAug 20 -ngp 2

# MORTAR + ALM, frictional (Coulomb); requires an unsymmetric solver
contact 12 -master 1 -slave 2 -formulation mortar -enforce alm \
        -kn auto -kt auto -epsN auto -epsT auto -augTol 1e-8 -maxAug 20 \
        -friction coulomb -mu 0.3 -cohesion 0.0
system UmfPackGen           ;# unsymmetric (Coulomb C_sl); Tresca stays symmetric
algorithm LadrunoAugmentedNewton -augTol 1e-8 -maxAug 20   ;# the Uzawa OUTER loop driver

# mesh-tying (zero-gap permanent mortar = the active set frozen ON; no inequality)
contact 13 -master 1 -slave 2 -formulation mortar -enforce penalty -tie
```

```python
# openseespy mirror (selector style)
ops.constraints('LadrunoContact')
ops.contactSurface(1, '-kind','masterSegments','-faces', masterFaces)
ops.contactSurface(2, '-kind','slaveSegments',  '-faces', slaveFaces)
ops.contact(11, '-master',1,'-slave',2,'-formulation','mortar','-enforce','alm',
                '-kn','auto','-epsN','auto','-augTol',1e-8,'-maxAug',20,'-ngp',2)
ops.algorithm('LadrunoAugmentedNewton','-augTol',1e-8,'-maxAug',20)
```

### Class hierarchy / data flow

```
Domain
 └─ LadrunoContactDomain*            (ADR-39; owns surfaces + contact defs + PATH STATE)
     ├─ LadrunoContactSurface        (ADR-39; +slaveSegments kind in ADR-40)
     ├─ Contact{ formulation, enforce, epsN, epsT, mu, augTol, maxAug, ... }
     └─ for each mortar contact:
         └─ LadrunoMortarSegment      (narrow phase; one per slave facet over an epoch)
             └─ LadrunoMortarPair     (Domain-side per-GP state: s_e_n, λ_N, λ_T, frame)

handle()  [LadrunoContactHandler, 33002]
 └─ injects one LadrunoContactFE adapter per contact def (runtime tag, NO classTag)
     getResidual(I) -> Σ_GP w·J·B_gpᵀ·t(g̃, slip, λ)        [both integrator families]
     getTangent(I)  -> Σ_GP w·J·B_gpᵀ·D·B_gp  scaled per the integrator (see below)
        - NTS  branch -> (ADR-39 planned) single-node B
        - MORT branch -> LadrunoMortarSegment clipped-GP B
     state lives on LadrunoMortarPair (Domain) -> adapter is a STATELESS VIEW (ADR-39 rule)

LadrunoAugmentedNewton::solveCurrentStep()   [the OUTER Uzawa loop]
 repeat:  inner NewtonRaphson to equilibrium at FROZEN λ
          read ‖g̃‖_∞ from LadrunoContactDomain
          λ_N += epsN·⟨g̃_N⟩₋  (clamp λ_N ≤ 0, compression-only KKT) ; λ_T outer update
 until ‖g̃‖_∞ < augTol or maxAug ;  then return -> ONE Domain::commit() promotes state
```

### Mortar integration scheme (clipped Gauss-point mortar)

For each slave facet paired (via the ADR-39 broad phase) with candidate master facet(s):

1. **Project** the slave facet onto the master plane; **clip** the slave/master overlap into a
   polygon on the master parametric plane; **triangulate** it; place a low-order Gauss rule
   (default 3-pt per sub-triangle; `-ngp` for the slave-facet rule order). *Clipping is in the
   MVP, not deferred* — without it the partition-of-unity `Σ_K M_IK = Σ_J D_IJ` is violated at
   master inter-element boundaries and the patch test fails (every mechanics reviewer flagged
   this; D1/D4's un-clipped "mortar-lite" headline was over-stated).
2. At each integration point: `LadrunoContactProjection` → `{gap_N, n, t1, t2, g[2][2], phi_master}`.
   Accumulate `D_IJ += w·J·N_I^s N_J^s`, `M_IK += w·J·N_I^s phi_K^m`.
3. **Weighted (variationally consistent) gap**, normal **inside** the integral:
   `g̃_I = Σ_GP w·J·N_I^s · ( (x_s − x_m(ξ̄))·n )` — *not* `n_I·(D x_s − M x_m)` with `n` pulled
   out (that factoring is exact only for a flat single-normal facet; D1's form was flagged).
4. **B-operator** `B_gp` spans the slave-facet nodes ∪ reachable master-segment nodes
   (the ADR-39 **conservative-static-superset** connectivity over a topology epoch).
5. `F_c += w·J·B_gpᵀ·t`, `K_c += w·J·B_gpᵀ·D·B_gp`, with `t`,`D` from the FrictionLaw + normal law.

What is **sacrificed vs full dual-mortar** (stated honestly): standard (not biorthogonal)
multiplier basis → `D` is **not** diagonal, so the per-node `λ` update is local only after the
per-facet `D`-solve; and it is **not** LBB/inf-sup optimal (deferred to ADR-42). ALM at finite
`epsN` is the mitigation; the patch-test gate **reports** any residual pressure oscillation.

### FrictionalLaw interface (C++ signature sketch)

Header-only, `namespace ladruno_contact_kernel`, plain `double`s + `<math.h>`, no OpenSees
types — numpy-oracle-testable without building OpenSees. The **surface metric `g` is an
input** (slip tangent uses `g`, not `I₂` — D4's `I₂` form was flagged wrong off
flat/orthonormal facets; the verified source builds `R = g·r`).

```cpp
namespace ladruno_contact_kernel {

enum FrictionType { COULOMB = 0, TRESCA = 1 };

struct FrictionParams {
  int    type;        // COULOMB | TRESCA
  double mu;          // Coulomb friction coefficient (0 for pure Tresca)
  double cohesion;    // c, adhesive intercept
  double tauMax;      // Tresca shear cap (pressure-independent)
  double epsT;        // tangential penalty / regularization
};

// STATELESS pure return map (one slave Gauss point), tangential plane.
//   tN        : current normal traction (compression < 0; |tN| used in the cone)
//   gT[2]     : incremental tangential gap (slip increment) in the {t1,t2} frame
//   gMetric   : 2x2 covariant surface metric g_ab = g_a · g_b
//   s_e_n[2]  : COMMITTED elastic slip (in), s_e_out[2] : trial elastic slip (out, to commit)
//   lambdaT[2]: ALM tangential multiplier (in/out; pass null for pure penalty)
//   tT[2]     : returned tangential traction
//   Css[2][2] : consistent tangential tangent  dtT/dgT
//   Csl[2]    : NON-SYMMETRIC pressure-coupling dtT/dtN  (Coulomb: mu·sign(tN)·r ; Tresca: 0)
//   slip      : 1 if slipping, 0 if sticking
void integrate(double tN, const double gT[2], const double gMetric[2][2],
               const FrictionParams& p,
               const double s_e_n[2], double s_e_out[2],
               double* lambdaT /*nullable*/,
               double tT[2], double Css[2][2], double Csl[2], int* slip);

// One-sided normal KKT law (compression only). gN < 0 => contact.
//   ALM:    tN = -(lambdaN + epsN*<-gN>)         (Uzawa, lambdaN stored on the pair)
//   penalty: lambdaN = 0
void normalTraction(double gN, double lambdaN, double epsN, int alm,
                    double* tN, double* dtN_dgN, int* active);

} // namespace
```

**Coulomb vs Tresca = one code path**: cone `f = ‖tT‖_g − (mu·|tN| + cohesion)` (Coulomb) vs
`f = ‖tT‖_g − tauMax` (Tresca); slip is a radial return `tT = scale·tT*`, `scale = cap/‖tT*‖_g`,
consistent `Css = epsT·scale·(g − r⊗r)` with `r = tT*/‖tT*‖_g`, and `Csl = mu·sign(tN)·R`
(Coulomb, **non-symmetric**) vs `Csl = 0` (Tresca → **symmetric**, the safe first bring-up).

**Sharing with ADR-39 (the central question resolved):** ADR-40 **writes**
`LadrunoFrictionKernel.h`; ADR-39's planned P3 NTS friction will **adopt** it (its explicit
IMPL-EX branch is a thin wrapper that extrapolates the committed state and calls the same
`integrate()`). There is **no existing `LadrunoContactFriction` to refactor** — D2/D3/D4's
"refactor ADR-39's friction with bit-for-bit regression" premise is false (verified: no such
class in the tree). The correct ordering is **forward adoption**, and it only holds if ADR-40's
kernel lands before ADR-39 reaches P3.

### Integration points

- **getTangent / c1 contract (resolved concretely, not "as ADR-39").** The shipped
  `LadrunoContactFE::getTangent` **overrides** the base and returns its own matrix, so the base
  `formEleTangent → addKtToTang(c1)` chain is bypassed. Returning an *already-c1-scaled* tangent
  **and** inheriting the base `addKtToTang` would apply `c1` **twice**. Decision: the adapter
  **fetches `c1` itself** inside the overridden `getTangent(Integrator*)` and returns
  `c1·K_c`, **and overrides `addKtToTang` to a no-op-beyond-its-own-scatter** so no double
  scaling occurs. `c1` is obtained via `IncrementalIntegrator::getCFactor()` (the portable
  virtual); under explicit `CentralDifferenceLadruno` the adapter returns a **zero** tangent so
  the LHS stays mass-only. This is a P1 gate (FD-vs-`K_c` on a rotated config catches a wrong
  factor immediately).
- **The Uzawa outer loop is an `EquiSolnAlgo`, NOT a Tcl `analyze`-loop.** A Tcl
  `while {…} {analyze 1}` re-solve advances `committedTime`, bumps `commitTag`, **fires recorders
  per augmentation**, and (under `LoadControl`) advances the load each call — all wrong. The
  augmentation must loop `solveCurrentStep` **pre-commit** at frozen load, then commit once.
  `LadrunoAugmentedNewton` (a `NewtonRaphson` subclass) owns that loop; the `λ` update happens
  **in the algorithm between inner solves**, and `Domain::commit → LadrunoContactDomain::commit`
  only **freezes** the converged `λ` and promotes `s_e`. (This is the single most important
  correction across all 16 reviews — D1/D2/D3/D4 all under-specified the driver.)
- **revertToLastCommit invariant.** `λ_N`,`λ_T` and `s_e` are **committed-only** (mutated solely
  in `commit`), never on a trial, so a rejected step's `revertToLastCommit` is automatically
  safe — identical to the EmbeddedRebar precedent (whose `revertToLastCommit` does *not* touch
  `lambda` because `lambda` only changes at `commitState`). `LadrunoMortarPair::revertToLastCommit`
  asserts this for the per-GP arrays.
- **Active-set / KKT (frictionless normal).** Active iff `λ_N + epsN·gN < 0`; `λ_N` clamped
  `≤ 0` (compression). The active set is **frozen within an augmentation sweep** (so re-solving
  between augmentations never triggers `domainChanged()` → no re-number); membership changes only
  between physical steps (ADR-39 epoch model, Q-EPOCH).
- **Segment-switch frame consistency.** When a GP re-projects onto a *different* master segment
  across an epoch, `λ_T` (expressed in the old tangent frame) is **re-projected** onto the new
  frame (the EmbeddedRebar `lambda -= (λ·dirCur)dirCur` move), never carried stale.
- **MP-constraint composition (Q-CONSTR, resolved not inherited).** `LadrunoContactHandler`
  (verified) does **not** enforce `MP_Constraint`s — it warns. ADR-40 makes this an **explicit
  P1 restriction**: a mortar slave/master node may **not** simultaneously be an `equalDOF` /
  `rigidDiaphragm` slave under this handler. The gate combines a mortar interface with a
  `rigidDiaphragm` and asserts the documented error (delegation to a base Transformation handler
  is an ADR-42 item, not in this scope).
- **Unsymmetric solver (Coulomb).** `Csl ≠ 0` makes `K_c` unsymmetric → frictional Coulomb
  **requires** an unsymmetric SOE (`UmfPackGen`/`FullGen`), exercised at P3; Tresca (`Csl=0`)
  stays symmetric and is the first frictional bring-up. Precedent: `LadrunoConcrete3D` already
  ships non-symmetric.

### Testing / oracle battery (numeric thresholds, named cross-checks)

| Phase | Gate (PASS threshold) | Oracle |
|---|---|---|
| **P0** | FrictionLaw matches reference to **1e-10**; `Css`,`Csl` vs FD to **1e-6**; `Csl≠0` Coulomb / `=0` Tresca. | `tests/contact/ladruno_friction_reference.py` (1D-plasticity radial return) |
| **P0.5** | Projection `gap/n/ξ̄` vs closed form on a **tilted plane** + a **known-curvature** facet; bounded-Newton sentinel never assembles garbage. | analytic |
| **P1** | Clipped `D`/`M`: row-sums `Σ_K M_IK = Σ_J D_IJ` (partition of unity) to **1e-12**; **constant-pressure patch test on a deliberately non-matched 2-block mesh** → interface pressure constant to **≤1e-6 relative**; `K_c` vs FD on a **rotated** config to **1e-6**. | numpy oracle + analytic patch |
| **P2 (MVP)** | Frictionless Uzawa: penetration `g̃ → O(epsN-independent tol)` within `maxAug` (the ALM win penalty cannot reach); **release/reopen → exact F=0**; **Hertz** sphere/cylinder `p(r)=p₀√(1−(r/a)²)` peak+radius converge under refinement; `E_contact ≥ 0`; SOE eqn count **constant across augmentations**. | analytic Hertz; cross-check vs **Kratos CSMA `ALMFrictionlessMortarContactCondition`** on an identical mesh (concept-level magnitude/trend, not bit-compare) |
| **P3** | Frictional: **sliding-block-on-incline** `a = g(sinθ − μcosθ)`, static for `tanθ<μ`; ironing/sliding-patch shear smooth + resultant; **Tresca cap**; Coulomb `Csl` FD-checked; unsymmetric path exercised; converged Coulomb answer **Δt-independent** (the implicit virtue IMPL-EX cannot match). | analytic + Kratos CSMA frictional cross-check |
| **P4** | Mesh-tying (`-tie`, zero-gap) transmits a uniform/linear traction across a non-matched interface = single-block stress (exact patch); large-sliding curved-interface convergence. | analytic |

`null-mortar` gate (all phases): with no mortar contact defined, the build is **byte-identical
to stock** (the empty adapter is graph-neutral). Note honestly: an **active** mortar adapter
declares real connectivity, so the numberer permutation differs from stock — the bitwise claim
holds **only** for the null case (D1 review fix).

## Risks / open questions

> [!question] Q-DOF (true Lagrange/traction DOFs vs Uzawa) — **RESOLVED: Uzawa.**
> OpenSees DOFs are positional (ndf/node) and the equation set freezes at
> `domainChanged() = clearAll → handle → numberDOF → setSize`. True LM DOFs are *feasible*
> (a node-less `DOF_Group` minted in `handle()`, as `LagrangeConstraintHandler` does) but the
> multiplier **count freezes at `handle()`**, so a churning contact active set forces a full
> re-number per change. `LagrangeDOF_Group` is **constraint-bound** (no `(tag,count)` ctor) and
> its `LagrangeMP_FE` assembles a `−α·C` coupling, **not** a mortar `Bᵀn` block — so it is *not*
> reusable verbatim (D2's claim corrected). **Decision:** ship Uzawa-over-penalty (`λ` as
> Domain-side per-GP state, zero new DOFs), grounded in the **verified in-fork EmbeddedRebar
> `commitState` AL precedent**. True-LM / saddle-point (with the inf-sup stabilization it
> genuinely needs) is **hard-deferred to ADR-42**.

> [!question] Q-DRIVER (the Uzawa outer loop) — **RESOLVED: custom EquiSolnAlgo.**
> The per-step `Domain::commit` update (the EmbeddedRebar pattern) gives augmentation that
> converges *across* load steps for free, but a true *within-step* augmentation loop (capped,
> gap-gated) is **not** free and **must not** be a Tcl `analyze`-loop (double-commit / recorder /
> load-advance corruption). `LadrunoAugmentedNewton` (a `NewtonRaphson` subclass, broker tag
> 33001) is the additive-leaf driver; it commits once per physical step.

> [!question] Q-MORTARLITE (full dual-mortar D/M vs mortar-lite) — **RESOLVED: clipped GP mortar, dual deferred.**
> Un-clipped slave-GP "mortar-lite" (D1/D4) does **not** pass the non-matched patch test
> (partition-of-unity broken at master element boundaries). ADR-40 ships **overlap clipping in
> the MVP** with a **standard** multiplier basis (`D` non-diagonal, per-facet solved). **Dual
> (biorthogonal) basis** (diagonal `D`, cheap nodal `λ`) and **LBB-optimal** treatment are
> **deferred to ADR-42** — finite-`epsN` ALM is the interim mitigation; the patch gate reports
> any residual oscillation rather than hiding it.

> [!question] Q-DEP (ADR-39 maturity) — **OPEN dependency, stated not assumed.**
> ADR-39 shipped **P1 only** (zero adapter; no broad phase, narrow phase, projection kernel, or
> friction). ADR-40 must either be **sequenced after ADR-39 P2 (broad phase + projection)**, or
> **co-own** the projection kernel and adapter-connectivity work as its own P0.5/P1 deliverables.
> This ADR takes the **co-own** stance for the projection + friction kernels (they are on the
> critical path for both ADRs) and **depends on** ADR-39 P2.5 for the broad-phase pair candidates.

> [!question] Q-CONSTR (rigidDiaphragm/equalDOF composition) — **RESOLVED for this scope: restricted.**
> The contact handler does not enforce MP constraints. A mortar contact node may not also be an
> MP slave; gated at P1 with a `rigidDiaphragm`+mortar model asserting the documented error.
> Base-handler delegation is ADR-42.

> [!question] Q-EPOCH / Q-GRAN (inherited from ADR-39, sharpened for mortar)
> Mortar adapter connectivity (slave facet ∪ reachable master nodes) is **wider** than NTS →
> denser frozen `getID` supersets, larger fill, more frequent re-number under sustained sliding.
> Add a mortar-specific epoch-cost gate (re-number frequency + fill growth vs sliding distance).

> [!question] Q-IMPLFILL (implicit SOE fill) — **OPEN.**
> ADR-39's "system Diagonal costs nothing" holds only for the explicit LHS. The implicit mortar
> tangent `c1·epsN·BᵀB` is genuine off-diagonal coupling (slave facet ↔ master nodes) that a
> `Diagonal` SOE cannot store — mortar **mandates** a profile/sparse/`UmfPack` SOE; document and
> profile fill/bandwidth growth.

> [!question] Q-EXPLICIT (explicit ALM) — **DISCLOSED limitation.**
> Under `CentralDifferenceLadruno` the tangent is mass-only, so Uzawa degenerates to a single-pass
> penalty force. Full ALM mortar is **implicit-only** (this is the accuracy-first lane,
> complementary to ADR-39's explicit-first NTS). Do not promise explicit ALM.

- **Dependencies:** header-only kernels (no external deps); numpy for oracles; an unsymmetric SOE
  (`UmfPackGen`, already in-tree) for frictional Coulomb.
- **Numerical:** Uzawa is linearly convergent — provide `-epsN auto` (reuse ADR-39's `-kn auto`
  rule) + `maxAug` + an `epsN`-ramp heuristic; gate the `epsN`-independence claim explicitly.
- **Backwards compatibility:** purely additive; `null-mortar` is byte-identical to stock.

## Implementation log

*(filled in once the plan is being executed; move to `Ladruno_internal/implemented_mortar_alm_contact.md` when done)*

## Adversarial review log

Four independent designs (D1 mesh-tying-first, D2 full-ALM-mortar, D3 extend-ContactDomain,
D4 frictionlaw-first), each scored by four lenses: **A** = Assembly/DOF/enforcement,
**M** = Mortar mechanics/correctness, **R** = Reuse/footprint/fork-fit, **S** = Scope/phasing/risk.

| Design | A | M | R | S | Worst flaws found |
|---|---|---|---|---|---|
| **D1** mesh-tying-first | 8 | 6 | 8 | 7 | Patch-test claim contradicts un-clipped slave-GP integration; `g̃` factors `n` out of the integral (wrong off flat); "docks onto shipped state" — state machinery is P1b stubs; classTag premise false. |
| **D2** full-ALM-mortar | 8 | 7 | 8 | **4** | **Fatal:** assumes ADR-39 machinery exists "verbatim" — it does not; `LagrangeDOF_Group(tag,nLM)` node-less ctor **does not exist**; P1 "refactor ADR-39 friction" refactors nonexistent code; implicit SOE fill cost (Q-IMPLFILL) unaddressed. |
| **D3** extend-ContactDomain | 7 | 5 | 7 | 7 | Uzawa driver as Tcl `analyze`-loop double-commits/fires recorders (must be an `EquiSolnAlgo`); collocated mortar-lite likely **fails** the patch test; "reuse ADR-39 P2b kernel/bucket sort" — unbuilt; classTags confused across registries. |
| **D4** frictionlaw-first | 8 | 6 | **9** | 6 | "Three nested loops"/capped outer Uzawa is **not** free in the analysis driver; GPTS without overlap clipping fails the patch test; slip tangent uses `I₂` not metric `g`; "refactor ADR-39 `LadrunoContactFriction`" — class doesn't exist; classTag 33009/33010 collide with ADR-26 reservations. |

**Synthesis decision.** Spine = **D4 (frictionlaw-first + per-GP Uzawa on the verified
EmbeddedRebar `commitState` precedent)** — highest reuse/footprint score and the only design that
correctly grounds zero-DOF ALM in a *real in-fork* precedent. Grafted:
- **D3's "extension not sibling"** framing + the `-formulation/-enforce` selector command
  (one subsystem, two formulations) — cleaner than D1/D2's parallel-plumbing language.
- **D1's mesh-tying** insight folded in as the **degenerate `-tie` (active-set-frozen-ON)** case
  of the same mortar kernel (P4), not a separate MVP — it exercises the kernel but mesh-tying
  alone needs no `λ`, so it is not the headline.
- **The metric-aware FrictionalLaw interface** (D4/D2) with `g` as an input and the verified
  `C_ss/C_sl` from `ContactMaterial3D` — correcting D4's `I₂` error.
- **D2/D3's honest implicit-SOE-fill and unsymmetric-solver** risk accounting (Q-IMPLFILL).

**Required-fix dispositions (folded in):**
1. **Patch test ⇄ mortar-lite contradiction** (all M lenses): **overlap clipping moved into the
   P1 MVP**; the headline is "passes the non-matched patch test" with clipping, and the gate is a
   hard ≤1e-6 numeric falsifier — *not* an oscillation-magnitude metric.
2. **Uzawa driver** (all A/S lenses): **custom `LadrunoAugmentedNewton` `EquiSolnAlgo`**, never a
   Tcl `analyze`-loop; commits once per step; `λ` updated between inner solves on a frozen set.
3. **c1 double-scaling** (D2/D4 A): adapter fetches `c1` via `getCFactor()`, returns `c1·K_c`,
   overrides `addKtToTang` to avoid re-scaling; zero tangent under explicit.
4. **classTags** (all R lenses): the brief's "33015 taken by LadrunoContactFE" is **false**; MVP
   needs **no `ELE_TAG`**; only a new **algorithm tag 33001** is reserved; D4's 33009/33010
   suggestion (collides with ADR-26) **dropped**; `FrictionKernel` explicitly **not** an `FRN_TAG`.
5. **ADR-39 maturity** (D2/D3/D4 S — fatal/major): premise that ADR-39 ships a working penalty
   engine + a `LadrunoContactFriction` class is **false**; reframed as **co-own the shared
   kernels + depend on ADR-39 P2.5 broad phase** (Q-DEP). No "bit-for-bit regression" of code
   that doesn't exist — friction sharing is **forward adoption**.
6. **`g̃` consistency** (D1 M): normal kept **inside** the integral; flat-facet factoring rejected.
7. **Metric in slip tangent** (D4 M): `g` is a kernel input; `R = g·r`.
8. **Q-CONSTR / Q-EPOCH / Q-IMPLFILL / Q-EXPLICIT**: resolved-with-restriction or
   disclosed-with-gate, not "inherited" (ADR-39 shipped scaffolding, so there is nothing to inherit).
9. **Single-maintainer realism** (all S): committed scope **cut to P0→P4**; dual basis, true-LM,
   self-contact, NTN/NTS-via-mortar-weights **hard-deferred to ADR-42**.

**Hard-defer ledger:** full dual/biorthogonal mortar (Q-MORTARLITE), true-LM saddle-point +
inf-sup stabilization (Q-DOF), self-contact, simplified MPC/NTN/NTS-via-mortar-weights,
base-handler MP delegation (Q-CONSTR) → all **ADR-42**, each with the rejection reason above.
