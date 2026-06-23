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

> **Part of the definitive contact plan — see the [[48_ladruno_contact_capstone_adr]] capstone**
> (architecture, contracts, status-of-record, unified roadmap). ADR-41 is the **detailed design of
> record for the mortar / ALM accuracy-first lane**; global status & sequencing live in the capstone.
> The "Sequencing reality" and Q-DEP notes below were corrected 2026-06-22 against the live tree
> (ADR-39 has shipped P1→P3.5; this lane docks onto a working penalty NTS engine).

## What

ADR-41 is the **implicit, accuracy-first** complement to ADR-39's explicit, robust-first
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
Uzawa ALM (the MVP), and frictional mortar. **Hard-deferred to a successor ADR-47:** dual
(biorthogonal) Lagrange shape functions, true-LM/saddle-point enforcement, self-contact, and
NTN/NTS-via-mortar-weights. See *Risks* for the rejection rationale on each.

**Sequencing reality (UPDATED 2026-06-22 against the live tree — the original "P1 only" framing
below is STALE and was rewritten).** ADR-39 has since shipped **P1 → P3.5** (PRs #345–#361): a
working penalty NTS engine. Verified in the tree today:
- **Broad phase shipped** — `SRC/domain/contact/LadrunoContactBucketSort.h` (P2.5, #358).
- **Narrow phase + projection kernel shipped** — `SRC/domain/contact/LadrunoContactKernel.h`:
  bounded closest-point projection Newton (`project`, ~L105, `maxIt` cap + `|detK|` degenerate
  reject), winding-immune normal, gap, penalty traction (P2b, #354).
- **Friction return map + consistent tangent shipped** — same header: `frictionReturnMap` (Coulomb,
  P3 #360) and `frictionTangentBlock` (the consistent friction tangent, symmetric default +
  `-consistanttan` non-symmetric opt-in, P3.5 #361).
- **Real adapter connectivity shipped** — `LadrunoContactFE` SEGMENT ctor builds
  `FE_Element(tag, 1+nps, 3·(1+nps))` and assembles real `Bᵀ·tn` + friction; **not** a zero adapter.
- **Real path-state + commit shipped** — `LadrunoContactDomain::commit()` promotes committed friction
  slip (`gpT = gpTtrial`) per per-pair `FrictionState`; **not** a no-op counter.

So ADR-41 **does dock onto a working penalty NTS engine.** Of its five planned new files, **two
already have shipped equivalents** inside `LadrunoContactKernel.h` (projection + friction kernel),
and **three are genuinely mortar-specific and absent** (`LadrunoMortarKernel`, `LadrunoMortarPair`,
`LadrunoMortarSegment`). The remaining ADR-41-specific work is therefore: the **mortar narrow phase
with overlap clipping + D/M**, **Uzawa ALM** (now commit-cycle — Q-DRIVER), the **per-GP mortar
state**, and **optionally extracting** ADR-39's in-`LadrunoContactKernel.h` friction math into the
shared header-only `LadrunoFrictionKernel.h` (consolidation, not a greenfield build). See Q-DEP
(re-resolved) for the corrected dependency status.

> *Original (now-stale) text, preserved for the record:* "ADR-39 has shipped **P1 only** —
> `LadrunoContactFE` is an empty-connectivity zero adapter … no broad phase, no narrow phase, no
> projection kernel, and no friction class … ADR-41 **co-owns and builds** the shared narrow-phase
> machinery." This was true at drafting; ADR-39 has since shipped that machinery (P2→P3.5).

## Why

- **The accuracy bar.** ADR-39 itself documents (`39_..._adr.md` §What, Risk Q-PATCH) that NTS
  "does not pass the contact patch test … interface-stress fidelity awaits mortar (v3)." The
  Kratos `ContactStructuralMechanicsApplication` reference sets exactly this bar: mortar
  conditions with ALM and penalty, frictionless + frictional, passing the patch test. ADR-41
  delivers that bar earlier and at lower cost than full dual-mortar by using a **clipped
  Gauss-point mortar** + **Uzawa ALM**.
- **Exact constraint at finite penalty.** Pure penalty (ADR-39) leaves residual penetration
  `g = P/k_n`. ALM/Uzawa drives `g → 0` at *finite* `k_n` over an outer augmentation loop,
  resolving ADR-39 Risk Q-AL.
- **One friction law, two consumers.** A return map is the single most reused, most
  bug-prone, most oracle-testable piece of any contact code. Factoring it header-only (mirroring
  the proven `LadrunoJ2Kernel.h` discipline) and **extracting it from ADR-39's already-shipped
  P3/P3.5 friction math in `LadrunoContactKernel.h`** (PRs #360/#361) gives both the ADR-39 NTS path
  and this mortar path one oracle-tested return map. (Build direction reversed — see §How / Q-DEP.)
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
| `SRC/domain/contact/LadrunoMortarPair.{h,cpp}` | **Domain-side path state** for one slave-facet ↔ master pairing: per slave-GP `{committed elastic slip s_e_n, λ_N, λ_T, slipFlag, last master-segment id + frame}`; `commit()`/`revertToLastCommit()`. The stateless-view target for the FE adapter (the ADR-41 analogue of ADR-39's planned `LadrunoContactPair`). |
| `SRC/domain/contact/LadrunoMortarSegment.{h,cpp}` | Narrow-phase formulation object owned by `LadrunoContactDomain`: drives the GP loop, calls `LadrunoContactProjection` + `LadrunoMortarKernel` + `LadrunoFrictionKernel`, returns `F_c` and `K_c`. Mortar sibling of the (planned ADR-39) NTS narrow phase. |
| ~~`SRC/analysis/algorithm/equiSolnAlgo/LadrunoAugmentedNewton.{h,cpp}`~~ **(DEFERRED — see Q-DRIVER revision + the Abaqus-scope addendum)** | **No longer a P2 deliverable.** The MVP augments **per `Domain::commit()`** (the verified `LadrunoEmbeddedRebar::commitState()` precedent), so the multiplier update needs **zero new analysis-layer code** — it rides the `LadrunoContactDomain::commit()` seam ADR-39 already owns. A bespoke `NewtonRaphson` subclass is an **optional later upgrade**, built only if a named P2/P3 gate fails *specifically because within-step augmentation is required* and the held-load driver loop is insufficient. Tag `33001` stays **reserved-but-unbuilt** (same reserve-when-needed posture as the `ELE_TAG`). |
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
  exact), and **re-projects `lambda` onto the current frame** under corotation. ADR-41 clones
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
selector and ALM options. `nts` (ADR-39) and `mortar` (ADR-41) are siblings on one command.

```tcl
# constraint handler (ADR-39, reused unchanged)
constraints LadrunoContact

# define meshed surfaces (ADR-39 LadrunoContactSurface; mortar needs the slave meshed
# as SEGMENTS, not just nodes -> -kind slaveSegments)
contactSurface 1 -kind masterSegments -faces $masterFaceSet
contactSurface 2 -kind slaveSegments  -faces $slaveFaceSet

# NTS penalty (ADR-39)
contact 10 -master 1 -slave 2 -formulation nts   -kn auto -kt auto -mu 0.3

# MORTAR + ALM (ADR-41) — frictionless
contact 11 -master 1 -slave 2 -formulation mortar -enforce alm \
        -kn auto -epsN auto -augTol 1e-8 -maxAug 20 -ngp 2

# MORTAR + ALM, frictional (Coulomb); requires an unsymmetric solver
contact 12 -master 1 -slave 2 -formulation mortar -enforce alm \
        -kn auto -kt auto -epsN auto -epsT auto -augTol 1e-8 -maxAug 20 \
        -friction coulomb -mu 0.3 -cohesion 0.0
system UmfPackGen           ;# unsymmetric (Coulomb C_sl); Tresca stays symmetric
algorithm Newton            ;# MVP: stock Newton — λ augments per Domain::commit (EmbeddedRebar pattern)
# within-step augmentation, when a gate proves it needed, is a documented held-load proc:
#   analyzeAugmented $augTol $maxAug    ;# zero-increment re-commits reading ‖g̃‖ (NOT a custom algorithm)
# the bespoke `algorithm LadrunoAugmentedNewton ...` is a DEFERRED upgrade (see Q-DRIVER), not the MVP.

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
ops.algorithm('Newton')   # MVP: λ augments per commit (EmbeddedRebar pattern); see Q-DRIVER.
# within-step augmentation = a held-load analyzeAugmented(augTol, maxAug) recipe, not a custom algorithm.
```

### Class hierarchy / data flow

```
Domain
 └─ LadrunoContactDomain*            (ADR-39; owns surfaces + contact defs + PATH STATE)
     ├─ LadrunoContactSurface        (ADR-39; +slaveSegments kind in ADR-41)
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

MVP augmentation = commit-cycle Uzawa on STOCK NewtonRaphson (EmbeddedRebar precedent):
 stock analyze step: newStep -> NewtonRaphson::solveCurrentStep (FROZEN λ) -> Domain::commit
   Domain::commit -> LadrunoContactDomain::commit():
       λ_N += epsN·⟨g̃_N⟩₋  (clamp λ_N ≤ 0, compression-only KKT) ; λ_T outer update
       promote s_e ; one Uzawa update per commit  (exactly LadrunoEmbeddedRebar::commitState)
 within-step convergence (if needed): driver holds load fixed (zero-increment integrator)
   and re-commits until ‖g̃‖_∞ < augTol — a documented Tcl/Py `analyzeAugmented` proc,
   NOT a custom EquiSolnAlgo.   [bespoke LadrunoAugmentedNewton = deferred upgrade, see Q-DRIVER]
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
per-facet `D`-solve; and it is **not** LBB/inf-sup optimal (deferred to ADR-47). ALM at finite
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

**Coulomb vs Tresca = one cone, not two code paths** (sharpened against Abaqus TG §5.2.3,
which unifies them as `tau_crit = min(mu*p, tau_max)`): the single critical stress is
`cap = min(mu·|tN| + cohesion, tauMax)` — so a **Coulomb law with a pressure-independent cap**
(`mu>0` *and* finite `tauMax`) is just the general case, pure Tresca is its `mu=0` corner, and
pure Coulomb is `tauMax = ∞`. Slip is a radial return `tT = scale·tT*`,
`scale = cap/‖tT*‖_g`, consistent `Css = epsT·scale·(g − r⊗r)` with `r = tT*/‖tT*‖_g`. The
pressure-coupling block is **active only on the branch where the cap actually depends on `tN`**:
`Csl = mu·sign(tN)·R` when `mu·|tN|+cohesion < tauMax` (Coulomb-controlled, **non-symmetric**),
and `Csl = 0` when the `tauMax` cap is binding *or* `mu=0` (Tresca-controlled → **symmetric**,
the safe first bring-up). This `min()`-selected `Csl` is the FD-checked discriminator at the P3
gate (oracle must exercise both sides of the `min`).

**Sharing with ADR-39 (UPDATED — the build DIRECTION reversed).** The original plan was
"ADR-41 **writes** `LadrunoFrictionKernel.h` first; ADR-39 P3 **adopts** it (forward adoption)."
That premise is now **inverted by reality**: ADR-39 P3/P3.5 already **shipped** its friction return
map + consistent tangent **inside** `LadrunoContactKernel.h` (`frictionReturnMap` ~L211,
`frictionTangentBlock` ~L241) *before* ADR-41 wrote anything. So the corrected options are:
**(a) adopt/extend in place** — ADR-41's mortar narrow phase calls the existing
`LadrunoContactKernel.h` friction math directly (fastest; couples mortar to the NTS header); or
**(b) extract-then-share** — lift ADR-39's in-header friction math out into the standalone,
OpenSees-free, separately-oracle'd `LadrunoFrictionKernel.h` ADR-41 specified, and have *both* the
NTS path and the mortar path consume it (cleaner; a **refactor of shipped, validated code**, not a
greenfield write — and it must keep ADR-39's P3/P3.5 gates green bit-for-bit). Recommended: **(b)**,
because the mortar path needs the `λ_T`/ALM tangential form and the `min(μ|tN|+c, τ_max)` cone
unification anyway, and a shared header is the right home for both. Note: ADR-39's shipped kernel is
**penalty-`kt` NTS-flavored** and **dropped IMPL-EX from the v1 ship** (it survives only in the
Python oracle) — the extraction must generalize it to the mortar `λ_T` form, not just copy it.

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
- **The Uzawa augmentation rides the commit cycle (MVP), not a custom `EquiSolnAlgo`** —
  *revised after a source-grounded deep-dive of the analysis loop; supersedes the earlier
  "must be an `EquiSolnAlgo`" stance.* The `λ` update lives in `LadrunoContactDomain::commit()`
  (`Domain::commit` already calls it), updated **once per commit** exactly as
  `LadrunoEmbeddedRebar::commitState()` does (`lambda += kt·g`, frame re-projection, committed-only).
  This augments **across load steps** for free on **stock `NewtonRaphson`**. The earlier claim
  that a Tcl loop is categorically wrong was overstated: of its three objections only **`LoadControl`
  auto-advance** (`LoadControl::newStep`, `currentLambda += deltaLambda`) is a real correctness
  issue, and it is defused by a **zero-increment integrator** during the augmentation sweep;
  **recorder rows per augmentation** (`Domain::commit` fires recorders) are *cosmetic*; and
  **"double-commit"** is **not** a correctness hazard (commit is idempotent on a converged state —
  EmbeddedRebar relies on exactly this). Within-step augmentation, *when a gate proves it
  necessary*, is first attempted via a documented **Tcl/Py `analyzeAugmented` proc** (held-load
  zero-increment re-commits reading `‖g̃‖`), shipped as a tested recipe. A bespoke
  `LadrunoAugmentedNewton` is the **last-resort** upgrade — subclassing `NewtonRaphson` (whose
  header explicitly says "not expected … to be subclassed") is the fragile path, not the mandatory
  one. (See Q-DRIVER and the Abaqus-scope addendum for the full cost/risk delta and the promotion
  trigger.)
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
  (verified) does **not** enforce `MP_Constraint`s — it warns. ADR-41 makes this an **explicit
  P1 restriction**: a mortar slave/master node may **not** simultaneously be an `equalDOF` /
  `rigidDiaphragm` slave under this handler. The gate combines a mortar interface with a
  `rigidDiaphragm` and asserts the documented error (delegation to a base Transformation handler
  is an ADR-47 item, not in this scope).
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
> genuinely needs) is **hard-deferred to ADR-47**.

> [!question] Q-DRIVER (the Uzawa loop) — **RE-RESOLVED: commit-cycle Uzawa primary; custom EquiSolnAlgo deferred.**
> *Revised after a source-grounded deep-dive of `StaticAnalysis::analyze` / `IncrementalIntegrator::commit`
> / `Domain::commit` / `LoadControl::newStep` / `LadrunoEmbeddedRebar::commitState` (see the
> Abaqus-scope addendum for the full trace). Supersedes the prior "RESOLVED: custom EquiSolnAlgo".*
> **Primary (P2 MVP):** update `λ` **once per `Domain::commit`** inside `LadrunoContactDomain::commit()`
> — the verified in-fork EmbeddedRebar pattern, **zero new analysis-layer code**, on stock
> `NewtonRaphson`. This augments *across* load steps for free; the P2 gates test the *converged*
> answer, which it reaches.
> **The earlier "Tcl loop is wrong" argument was overstated.** Of its three objections: (1) `LoadControl`
> **load auto-advance** (`LoadControl::newStep`) is the only real correctness issue → defused by a
> **zero-increment integrator** during the sweep; (2) **recorder rows per augmentation** are *cosmetic*;
> (3) **"double-commit"** is **not** a correctness hazard (`Domain::commit` is idempotent on a converged
> state — EmbeddedRebar already depends on this). So *within-step* augmentation, when needed, is a
> **documented Tcl/Py `analyzeAugmented` proc** (held-load re-commits), not a class.
> **Deferred upgrade:** the bespoke `LadrunoAugmentedNewton` (`NewtonRaphson` subclass, tag 33001
> reserved-but-unbuilt) is built **only if** a named P2/P3 gate fails *specifically because within-step
> augmentation is required* **and** the `analyzeAugmented` proc proves insufficient — e.g. a single
> monotonic step that must land at `augTol` with no step structure to amortize across. Subclassing
> `NewtonRaphson` (whose header says it is "not expected … to be subclassed") is the fragile
> last resort, not the MVP mechanism.

> [!question] Q-MORTARLITE (full dual-mortar D/M vs mortar-lite) — **RESOLVED: clipped GP mortar, dual deferred.**
> Un-clipped slave-GP "mortar-lite" (D1/D4) does **not** pass the non-matched patch test
> (partition-of-unity broken at master element boundaries). ADR-41 ships **overlap clipping in
> the MVP** with a **standard** multiplier basis (`D` non-diagonal, per-facet solved). **Dual
> (biorthogonal) basis** (diagonal `D`, cheap nodal `λ`) and **LBB-optimal** treatment are
> **deferred to ADR-47** — finite-`epsN` ALM is the interim mitigation; the patch gate reports
> any residual oscillation rather than hiding it.

> [!question] Q-DEP (ADR-39 maturity) — **RE-RESOLVED 2026-06-22: dependency SATISFIED, not open.**
> *The original "ADR-39 shipped P1 only → ADR-41 must co-own and build the shared kernels" is STALE.*
> ADR-39 has shipped **P1 → P3.5** (#345–#361). Every machinery ADR-41 depended on now exists in the
> tree: **broad phase** (`LadrunoContactBucketSort.h`, P2.5), **projection kernel** + **normal/penalty
> law** (`LadrunoContactKernel.h`, P2b), **friction return map** (P3) + **consistent tangent** (P3.5),
> **real adapter connectivity** (`LadrunoContactFE` SEGMENT ctor), and **committed per-pair path state**
> (`LadrunoContactDomain::FrictionState` + real `commit()`). ADR-41 therefore **docks onto a working
> penalty NTS engine** and consumes these directly — it does **not** co-own/build them.
> **Residual dependency (the only one left):** if ADR-41 wants the shared header-only
> `LadrunoFrictionKernel.h`, that is a **refactor-extract** of ADR-39's shipped in-`LadrunoContactKernel.h`
> friction math (must keep P3/P3.5 gates green bit-for-bit), generalized to the mortar `λ_T` form — not
> a greenfield write. **ADR-39 PENDING work does NOT block ADR-41:** ADR-39 P4 (SOFT=1 Courant-stable
> penalty), P5 (SOFT=2 segment penalty), and P6 (mesh-tying) are **explicit-stability tiers**; ADR-41 is
> the implicit/accuracy-first lane and is orthogonal to them.

> [!question] Q-CONSTR (rigidDiaphragm/equalDOF composition) — **RESOLVED for this scope: restricted.**
> The contact handler does not enforce MP constraints. A mortar contact node may not also be an
> MP slave; gated at P1 with a `rigidDiaphragm`+mortar model asserting the documented error.
> Base-handler delegation is ADR-47.

> [!question] Q-NORMAL (faceted-master normal discontinuity / contact-point chatter) — **OPEN, raised by Abaqus TG §5.1.2.**
> The mortar weighted gap `g̃` integrates over the slave facet, but the **normal `n` at each GP is
> still taken from a single faceted master segment**, so `n` *jumps* across master inter-element
> boundaries. Abaqus documents (§5.1.2) that exactly this slope discontinuity makes the contact
> point **oscillate between segments** under sliding, and mitigates it with **slide-line smoothing**
> (C1-continuous rounded junctions) — or, in small-sliding (§5.1.1), with an **averaged nodal-normal
> field `N(X)`** evaluated as a smooth combination of adjacent segment normals. Overlap clipping
> fixes the *partition-of-unity / pressure* problem (Q-MORTARLITE) but **not** this normal-direction
> chatter; the `λ_T` frame re-projection on segment switch (the EmbeddedRebar move) reduces *stale*
> tangents but does not smooth the normal itself. **Decision:** ship the faceted per-GP normal in
> P1–P4 (honest), add a **smoothing/chatter gate** to the P3 sliding-patch test (monitor active-GP
> normal flips per increment under sustained sliding), and **defer averaged-nodal-normal smoothing
> to ADR-47** alongside the dual basis — both are surface-representation upgrades, not enforcement
> changes. Optional interim mitigation: the Abaqus **viscous normal pressure**
> `p_visc = mu_c·v_rel` (§5.2.1) damps status flip-flop near the threshold (disallowed under any
> future arc-length lane, where velocity is undefined).

> [!question] Q-SLIDING (small- vs finite-sliding cost lane) — **OPEN optimization, framed by Abaqus TG §5.1.1/§5.1.2.**
> Abaqus splits contact into **small-sliding** (the master tangent plane each slave GP sees is
> frozen affine in the master node coordinates at step start — no re-search, **symmetric**, no
> re-number, correct under arbitrary rotation but only small tangential motion) and **finite-sliding**
> (re-project every increment). The ADR's "epoch model" (frozen connectivity superset + re-number
> between physical steps) is effectively *finite-sliding with a lagged active set*. For the large
> class of problems with small interface slip (seated/bolted joints, tie-like contact, the `-tie`
> P4 case), a `-sliding small` mode that freezes the affine projection plane per step would be
> **cheaper and symmetric** and would sidestep most of Q-EPOCH's re-number/fill cost. Flag as a
> post-P4 optimization lane, not MVP scope.

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

> [!question] Q-VISCOUS (viscous contact stabilization) — **OPEN, residue from the Abaqus adversarial scope. Best ROI of the residue.**
> Abaqus offers a velocity-proportional normal pressure `p_visc = mu_c·v_rel` (TG §5.2.1) to damp
> contact **chattering / snap-through** near status flips — exactly the fork's flagship
> **pounding / rocking / uplift** regime. **Verified gap:** no velocity-proportional contact term
> exists anywhere (`ZeroLengthContactASDimplex` even zeros its `getDamp()`), and ADR-41 currently
> carries it only as an *unfunded mitigation note* under Q-NORMAL. It is **cheap** — a
> `v_rel`-proportional residual term + a `getDamp()` contribution on the adapter (ASDimplex already
> exposes the `getDamp()` seam). **Recommendation:** promote to a funded, gated **P3.5/P4 option**
> (`-visc <mu_c>` on the contact def), off by default, auto-disabled under any future arc-length
> lane (velocity undefined). Not in the committed P0→P4 spine yet — decide before P3 freezes.

> [!question] Q-MUDEP (pressure-/velocity-dependent friction coefficient) — **OPEN residue, adopt-later.**
> Abaqus's `mu` may depend on pressure / slip-rate / temperature (TG §5.2.3). The fork **already
> has the machinery** in the *wrong module*: `SRC/element/frictionBearing/frictionModel/`
> (`VelDependent`, `VelPressureDep` (Constantinou 1996), `VelNormalFrcDep`, `FrictionModel::setTrial`)
> computes `mu(N,v)` — but for **isolator bearings**, never wired to a contact surface.
> Velocity-weakening `mu(v)` governs sliding-interface seismic dissipation, so this is a real
> structural payoff. **Recommendation:** after the constant-`mu` kernel ships (P3), let
> `LadrunoFrictionKernel` take `mu` from a `FrictionModel`-style callback instead of a constant.
> **Adopt later**, not MVP.

> [!question] Q-ANISO (anisotropic / elliptic friction) — **OPEN residue, defer to ADR-47.**
> Two principal `mu` (friction ellipse, TG §5.2.3) via the scaled-shear transform — the existing
> radial return is reused after scaling, so it is a clean extension, not a rewrite. But it is
> **niche for structural-seismic** (orthotropic rock joints, laminated/fabric interfaces) and
> doubles the return-map state + oracle surface. **Recommendation:** defer to **ADR-47**. If ever
> pursued, read TG §5.2.3 Eq.5.2.3-10/11 directly — the skill flags its ellipse RHS as schematic,
> not the verbatim manual normalization.

- **Dependencies:** header-only kernels (no external deps); numpy for oracles; an unsymmetric SOE
  (`UmfPackGen`, already in-tree) for frictional Coulomb.
- **Numerical:** Uzawa is linearly convergent — provide `-epsN auto` (reuse ADR-39's `-kn auto`
  rule) + `maxAug` + an `epsN`-ramp heuristic; gate the `epsN`-independence claim explicitly.
  - **`-epsT auto` from an elastic-slip bound (Abaqus TG §5.2.3).** Do **not** ship a bare
    tangential-penalty number. Abaqus sizes its stick stiffness `k_stick` so the reversible
    "elastic slip" stays under `gamma_crit = 0.5%` of the average contact-element length; mirror
    that: `epsT_auto = c_T · cap / (gamma_crit · L_facet)` with `gamma_crit ≈ 5e-3·L_facet`, so a
    sticking GP slips at most ~0.5 % of a facet before the cone engages. This makes the tangential
    penalty **mesh-length-aware** (the dimensional sibling of ADR-39's `-kn auto`) rather than an
    opaque constant, and gives the P3 stick gate a physical pass tolerance (`|s_e| ≤ gamma_crit`)
    instead of an arbitrary one.
- **Backwards compatibility:** purely additive; `null-mortar` is byte-identical to stock.

## Abaqus Theory Guide cross-check (independent second source)

The ADR's mechanics were derived from Kratos `ContactStructuralMechanicsApplication`. The Abaqus
Theory Guide (Parts 5–6) is an **independent second source**; where it agrees it raises confidence
in a committed decision, where it diverges it surfaces a real option. This section is secondary to
the Kratos derivation and the source-verified OpenSees re-derivations above — it is *grounding*, not
a new requirement.

### Decisions Abaqus independently confirms (confidence ↑, no change)

| ADR decision | Abaqus confirmation (TG cite) |
|---|---|
| **ALM/Uzawa over a penalty kernel** (resolves Q-AL) | Abaqus/Standard "hard" contact *is* a mixed formulation `p = lambda + k*(h − h̄)` with a small **reference stiffness** `k*` and a **local Newton** driving `(h − h̄) → 0` — textbook augmented Lagrange (§5.2.1). Uzawa-over-penalty is the same family Abaqus ships in production, not a budget substitute. |
| **Reject true-LM / saddle-point** (Q-DOF) | A pure Lagrange multiplier `lambda` has no self-stiffness → a **zero on the Jacobian diagonal** that must be guarded against rigid-body modes; the reference stiffness `k*` exists precisely to remove it (§5.2.1). The inf-sup/LBB concern the ADR cites for deferring dual-mortar is the same pathology. |
| **Coulomb ⇒ unsymmetric solver** (`Csl≠0`, `UmfPackGen`) | The consistent friction Jacobian is "**non-symmetric whenever sliding occurs**" — Abaqus *recommends the unsymmetric solver* for frictional sliding to keep Newton quadratic (§5.2.3). Exactly the ADR's P3 requirement. |
| **Friction = rate-independent plasticity return map** | Abaqus states the friction constitutive structure "mirrors rate-independent plasticity": `tau_crit` = yield surface, stick = elastic region, slip = associated flow along `tau_i/tau_eq`, integrated by **backward-Euler radial return** (§5.2.3). Independent confirmation of the `ContactMaterial3D`-derived kernel. |
| **`g̃` (weighted gap) with the normal inside the integral** | Abaqus builds overclosure `h = (x_s − x_proj)·n` from the *projected* normal and needs `dh`, `d²h` per point (§5.1, §5.2.1); it never factors a single `n` out of a surface integral. Matches the ADR's rejection of D1's flat-facet factoring. |

### The augmentation driver — Abaqus's local augmentation tipped the verdict (Q-DRIVER re-resolved)

- **Local (per-point) augmentation vs the global Uzawa outer loop.** Abaqus drives the
  augmentation with a **local Newton *inside the contact constitutive law*** at each point
  (`(h − h̄) → 0` to tolerance), so there is **no global outer augmentation loop** — the
  augmentation is invisible to the analysis driver. This (plus a source-grounded deep-dive of the
  OpenSees analysis loop) **flipped the original ADR choice**: the MVP now uses **commit-cycle
  Uzawa** (the `LadrunoEmbeddedRebar::commitState` pattern — `λ` updated once per `Domain::commit`,
  stock `NewtonRaphson`, **zero new analysis code**), and the bespoke global `LadrunoAugmentedNewton`
  is **deferred** behind a within-step-convergence trigger. See **Q-DRIVER** (re-resolved) and the
  adversarial-scope addendum below for the full cost/risk trace. Abaqus's local-augmentation
  posture is the production-grade confirmation that no global driver is required for correctness.

### Abaqus contact scope — adversarial review result (what OpenSees can learn, triaged)

A four-agent adversarial pass (three lenses — kinematics / enforcement / friction — plus a
ROI skeptic, each reading the live `SRC` contact tree and the TG) scoped the *whole* "what can
OpenSees learn from Abaqus, contact-wise" question. Headline: **the premise mostly inverts** —
11 of 12 candidate learnings are **already shipped or already committed** in ADR-39/41, and OpenSees
is arguably *ahead* in one place (`ZeroLengthContactASDimplex` has an IMPL-EX contact option
Abaqus/Standard lacks).

**Already shipped (do not reinvent):** friction-as-plasticity radial return + non-symmetric
consistent tangent (`ContactMaterial3D.cpp:333-334`); closest-point projection (`SimpleContact3D`,
`BeamContact*`); cohesion + tension cutoff; the unsymmetric-solver path; pressure/velocity-dependent
`mu` *machinery* (`frictionModel/`, just unwired); a second IMPL-EX friction implementation
(`ZeroLengthContactASDimplex`).

**Already planned here (folded last revision):** unified `min(mu|tN|+c, tauMax)` cone; `-epsT auto`
from `gamma_crit = 0.5%·L_facet`; AL/Uzawa at finite penalty; nodal-normal smoothing (Q-NORMAL,
→ ADR-47); small-vs-finite-sliding lane (Q-SLIDING).

**Genuine net-new residue (now risk-noted above):** **Q-VISCOUS** viscous stabilization
(best ROI — promote to a funded P3.5/P4 option), **Q-MUDEP** wire `mu(N,v)` into contact
(adopt-later), **Q-ANISO** elliptic friction (→ ADR-47). Plus the `SimpleContact3D::project()`
**unbounded Newton loop** defect (`SimpleContact3D.cpp:600`) — already on the ADR-41 P0.5 critical
path (bounded re-derivation).

**Unanimous skip / ADR-47 (scope-creep guard):** coupled multiphysics surfaces
(pore-pressure / thermal / Joule / frictional-heat / acoustic), pressure penetration, exact-Lagrange
stick, full slide-line Hermite smoothing, rigid *analytical*-surface subsystem, true-LM / dual-mortar.
All low-ROI for a single-maintainer structural-seismic fork or already on the ADR-47 ledger — keep
them out of the committed P0→P4.

### Net effect on the plan

Folded into the body: the unified `min(mu|tN|+c, tauMax)` cone (§5.2.3), the
`gamma_crit = 0.5%·L_facet` elastic-slip basis for `-epsT auto` (§5.2.3), the **Q-DRIVER
re-resolution** (commit-cycle primary, custom algorithm deferred), and five new risk notes
(**Q-NORMAL**, **Q-SLIDING**, **Q-VISCOUS**, **Q-MUDEP**, **Q-ANISO**). The only committed-scope
*change* is Q-DRIVER demoting `LadrunoAugmentedNewton` out of the P2 deliverable list (now a
gated upgrade); the P0→P4 spine is otherwise unchanged, and no new deferral leaves the ADR-47 ledger.

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
2. **Uzawa driver** (all A/S lenses): ~~custom `LadrunoAugmentedNewton` `EquiSolnAlgo`~~ —
   **SUPERSEDED by the Q-DRIVER re-resolution** (post-deep-dive). MVP is **commit-cycle Uzawa on
   stock `NewtonRaphson`** (λ per `Domain::commit`, EmbeddedRebar pattern); within-step augmentation
   is a held-load `analyzeAugmented` proc; the custom `EquiSolnAlgo` is a gated, deferred upgrade.
   The original disposition correctly rejected a *naïve* `analyze`-loop, but over-generalized that
   into "a custom algorithm is mandatory" — only `LoadControl` auto-advance is a real hazard, and a
   zero-increment integrator defuses it.
3. **c1 double-scaling** (D2/D4 A): adapter fetches `c1` via `getCFactor()`, returns `c1·K_c`,
   overrides `addKtToTang` to avoid re-scaling; zero tangent under explicit.
4. **classTags** (all R lenses): the brief's "33015 taken by LadrunoContactFE" is **false**; MVP
   needs **no `ELE_TAG`**; only a new **algorithm tag 33001** is reserved; D4's 33009/33010
   suggestion (collides with ADR-26) **dropped**; `FrictionKernel` explicitly **not** an `FRN_TAG`.
5. **ADR-39 maturity** (D2/D3/D4 S — fatal/major): **SUPERSEDED (drafting-time record; re-resolved
   — see Sequencing-reality §What / Q-DEP).** At drafting, the premise that ADR-39 ships a working
   penalty engine was treated as false and friction sharing as "forward adoption." Both are now
   **stale**: ADR-39 has since shipped P1→P3.5 (a working penalty NTS engine *with* the friction
   kernel), so the corrected framing is **dock onto the shipped engine + extract/share its in-tree
   friction math** (reversed direction), not co-own-and-build.
6. **`g̃` consistency** (D1 M): normal kept **inside** the integral; flat-facet factoring rejected.
7. **Metric in slip tangent** (D4 M): `g` is a kernel input; `R = g·r`.
8. **Q-CONSTR / Q-EPOCH / Q-IMPLFILL / Q-EXPLICIT**: resolved-with-restriction or
   disclosed-with-gate, not "inherited" (ADR-39 shipped scaffolding, so there is nothing to inherit).
9. **Single-maintainer realism** (all S): committed scope **cut to P0→P4**; dual basis, true-LM,
   self-contact, NTN/NTS-via-mortar-weights **hard-deferred to ADR-47**.

**Hard-defer ledger:** full dual/biorthogonal mortar (Q-MORTARLITE), true-LM saddle-point +
inf-sup stabilization (Q-DOF), self-contact, simplified MPC/NTN/NTS-via-mortar-weights,
base-handler MP delegation (Q-CONSTR) → all **ADR-47**, each with the rejection reason above.
