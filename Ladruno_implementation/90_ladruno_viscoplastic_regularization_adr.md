---
title: "ADR 90 — Duvaut–Lions viscoplastic regularization wrapper (LadrunoDuvautLions)"
project: Ladruno
type: ADR
status: "P0 complete (oracle + A0); P1 wrapper not started"
priority: high
owner: nmora
related:
  - "[[_adr90_regularization_planning_brief]]"
  - "[[_adr90_orchestration_plan]]"
  - "[[_adr90_a0_results]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[31_ladruno_concrete3d_adr]]"
  - "[[59_ladruno_gradient_concrete_adr]]"
  - "[[32_ladruno_dispbeamcolumn_regularization_adr]]"
  - "[[87_ladruno_depth_with_width_adr]]"
  - "[[testbed/00_canonical_testbed]]"
tags:
  - adr
  - material
  - wrapper
  - regularization
  - viscoplastic
  - duvaut-lions
  - localization
  - sanisand
  - tims
aliases:
  - ADR-90
  - LadrunoDuvautLions
  - duvaut-lions wrapper
updated: 2026-09-04
---

# ADR 90 — Duvaut–Lions viscoplastic regularization wrapper

> [!abstract] Status — P0 complete, P1 not started
> **Number 90 on purpose.** 88 is cited throughout `SRC/` and the ledgers by PR #778 (H5DRM
> higher-order elements); 89 is proposed for Track T by
> `86_geomechanics_libraries_and_mpm_scoping_report` §6.4. ND class tag **33022 is RESERVED here
> and not yet in `SRC/classTags.h`** — the ND registry's highest is 33021
> (`SRC/classTags.h:594`); the 33022 in the EigenSOE (`:57`) and PATTERN (`:646`) registries are
> per-registry and deliberately not collisions.
>
> **P0 is done and is what this ADR is written against.** The numpy oracle and the A0 bar
> (commits `cc7c7f7a5`, `8863a468c`, PR #783) answered the one question that had to be answered
> before any C++: *is a generic `NDMaterial` wrapper — which can only see `inner.getStress()`,
> hence is two-track — an adequate stand-in for true Duvaut–Lions?* The answer turned out to be a
> **theorem** rather than a tolerance (§4.3), and the owner decided **D2** on it at checkpoint CP1
> on 2026-09-04: build the generic two-track wrapper.
>
> **Nothing in `SRC/` has changed.** No build exists for this ADR. Every number below is numpy,
> regenerable by `python3.12 tests/_testbed/run_a0_sweep.py`.

---

## 1. Driver and problem — stated as a measurement, not an intuition

The TIMs / APE macroelement calibration campaign fits ultimate-surface loci from radial pushovers
on `LadrunoSANISAND`. Those loci live on the **post-peak / collapse** branch of a **non-associated
softening** material, which is where the rate-independent boundary-value problem loses ellipticity
and the answer stops depending on the material and starts depending on the mesh.

The fork has no gate for this. Every objectivity gate it ships regularizes **dissipated energy**
through a characteristic length `lch` — `test_ladrunoRCConcrete_meshobj.py`,
`test_ladrunoSolidShell_softening.py` G5, `test_lemaitre_notched_bar.py`,
`test_ladrunoBrick_asdconcrete*.py` — and one of them says so in its own docstring:
*"lch regularizes the dissipated ENERGY, not the localization WIDTH"*
(`tests/test_lemaitre_notched_bar.py:21-23`). Greps for `rudnicki|bifurcat|shear band|band width`
in `tests/` return nothing (brief F6).

**The collapse load of a non-associated softening strip depends on the thickness of the mechanism,
not only on the energy it dissipates.** Energy objectivity is therefore not sufficient for the
named consumer.

### 1.1 The size of the problem, measured

A0 (`Ladruno_implementation/_adr90_a0_results.md` §5.1) — 1-D softening bar, `L = 100 mm`,
`E = 20000 MPa`, `σ_Y = 20 MPa`, `H = −E/400`, 10 % imperfection, `N ∈ {20, 40, 80, 160}`,
rate-independent (`τ = 0`):

| N | h [mm] | w1 | w2 | w3 | w2/h | W₅₀ | W₅₀ ratio |
|---|---|---|---|---|---|---|---|
| 20 | 5.000 | 5.000 | 5.000 | 5.000 | 1.00 | 12.150 | — |
| 40 | 2.500 | 2.500 | 2.500 | 2.500 | 1.00 | 6.075 | 2.000 |
| 80 | 1.250 | 1.250 | 1.250 | 1.250 | 1.00 | 3.038 | 2.000 |
| 160 | 0.625 | 0.625 | 0.625 | 0.625 | 1.00 | 1.519 | 2.000 |

The band width **is** the element size, by all three definitions, exactly; the dissipated work
halves with `h`. And the ill-posedness has a second, free signature: the number of uniform load
steps the rate-independent problem needs in order to converge at all grows 4000 → 16000 → 32000 →
64000 across the same four meshes, while **every** regularized run completes at 250 steps on
**every** mesh (§5.8 of the results note).

## 2. Named consumer, and what is settled

**Named consumer (the ADR-59 un-descope gate (a)):** vault 65 **P1** — the ultimate-surface loci
fitted from radial pushovers on SANISAND, whose post-peak values are mesh-sensitive by
ill-posedness. This ADR exists to serve that consumer and no other; if P1 is withdrawn, this ADR
is withdrawn with it.

Settled on both sides; **not re-opened here** (brief §1):

| # | Settled | Source |
|---|---|---|
| S1 | The method is a Duvaut–Lions viscoplastic **wrapper at the `NDMaterial` level**, not a modified `ManzariDafalias`. | vault 64 §6; vault 65 **D5** |
| S2 | **τ is a declared numerical parameter, not a soil property.** The deliverable is *"mesh-independent given a declared τ"*, characterised on a Deborah number. Uniqueness of the non-associated collapse load is **not** restored, and that is disclosed. | vault 10 line 75; vault 65 D3 + D5 |
| S3 | Element frozen: `LadrunoBrick -formulation bbar`; **tetrahedra prohibited** on failure legs. | vault 65 D4; vault 71; fork gate `tests/test_r3_prandtl_collapse_gate.py` (#722) |
| S4 | **Provenance is output, not documentation:** every regularized number ships with engine hash, element/formulation, τ, rate, three-mesh band, ultimate criterion. | vault 65 D3; `ops.ladrunoBuild()` (#718) |
| S5 | The characterisation protocol: vary the relaxation parameter and the ramp duration **independently**; matched pairs must collapse onto the ratio. | vault 14 |
| S6 | Primary calibration instrument = the displacement-controlled radial monotonic curve; probes radial, swipe as cross-check. | vault 65 D2, D7 |
| S7 | The wrapper is the smaller half; the validation campaign is the larger half. | vault 64 §6 |
| S8 | The papers are in hand; no acquisition phase. | skill `tim-macroelement/references/library_map.md` |

## 3. The quasi-static contradiction, and the lane decision

Vault 65 makes the **displacement-controlled** radial probe the primary instrument (D2, D7) and
Duvaut–Lions the regularizer (D5). On this engine those are in tension three ways:

1. **The regularizing effect decays toward the rate-independent limit** and a slow static push
   *is* that limit (de Borst et al. 1993; brief F4). Their internal length `ℓ = 2mc_e/E` is a
   **wave-speed** quantity — under quasi-statics there is no intrinsic length at all. A0 §5.2
   makes this quantitative: the converged band width is a function of `De = τ/T` alone, and it
   moves a **factor 12 for a factor 3 in De** (3.53 mm at De = 3e‑4 → 42.0 mm at De = 1e‑3).
2. **Under `DisplacementControl` the pseudo-time increment is erratic.** `ops_Dt` is
   `Domain::dT` = current − committed pseudo-time (`Domain.cpp:2080-2082`, `:2125`, `:2392`);
   under `LoadControl(dλ)` it is uniform and positive, under `DisplacementControl` the λ-increment
   swings by ~1e5 and goes negative after `loadConst` (`LEDGER_quirks.md:1440-1443`, `:1612`).
   `β = Δt/(τ+Δt)` is then not even well defined step to step, so **De is set by the integrator,
   not by the deck**.
3. A0 §5.5 shows the transient is Δt-dependent even though PV3's steady overstress is not: `w2`
   moves 2.06 % over a 32× step refinement at fixed De, converging at first order.

### 3.1 Alternatives table (mandatory, ADR-59 kill-list item)

| Option | What it means | Cost | Verdict |
|---|---|---|---|
| **(a) Transient lane, τ replaces β_K** | Run the radial probes as slow damped transients — the reference model's *native* lane — with Rayleigh **β_K = 0** and the wrapper's τ as the one declared relaxation time. Δt is physical, `De = τ/T` is a deck quantity, and vault 14's machinery applies unchanged. Vault 14's β_K artifact was *already* a Kelvin–Voigt relaxation of 23 ms; this replaces an accidental regularizer with a declared one. | Wall clock (vault 14: ~17 min/run on the twin); the coupled lane is still blocked by the contact `ndf` guard. | **RECOMMENDED PRIMARY** |
| (b) Static lane under uniform `LoadControl` pseudo-time | Keep statics but drive the push with uniform pseudo-time so `Δt = 1/nsteps` exactly and `De = τ` at `T = 1`. Displacement targets via a displacement-controlled *load* series, not the `DisplacementControl` integrator. | Loses arc-length path-following past a limit point; vault 60 D16 stepping guards need re-verification. | **ACCEPTABLE FALLBACK** — must be gated identical to (a) at matched De on the P2 deck |
| (c) Strain-driven internal clock | Replace Δt with an accumulated strain measure so the "rate" is path-length based and integrator-independent. | A different constitutive model (Perzyna-in-arclength), no literature anchor, no closed-form gate. | **REJECTED** unless (a) and (b) both fail |
| (d) Do nothing to the lane; report the three-mesh band | Run `DisplacementControl` as-is and disclose the spread (vault 10's stance). | Zero engine work; the residue stays as large as §1.1 measures. | **This is the NEGATIVE CONTROL of the ADR, not an option** |

### 3.2 The headline claim, worded to survive review

> *The regularized quantity `q(De; h)` converges in `h` at fixed, declared De, and its
> De-dependence is measured and disclosed.*

Not "mesh-independent collapse load". A0 §5.2 shows why the weaker wording is the only honest one:
at De = 1e‑4 the load–displacement **curve** converges (`curveL2` 0.217 → 0.040) while the
**width** does not (ratios still tracking 0.5–0.86). *Width convergence and curve convergence are
different gates and must be reported separately.*

### 3.3 Reconciliation with the fork's own shipped alternatives (ADR-59 kill-list item)

| Shipped alternative | What it regularizes | Why it does not serve this consumer |
|---|---|---|
| **Crack-band `lch`** (`LadrunoConcrete3D`, `LadrunoRCConcrete`, `ASDConcrete3D`, Lemaitre) | The **dissipated energy** per unit area of band, by scaling the softening branch with the element's characteristic length. | It leaves the band exactly one element wide. Its own gate says so (`tests/test_lemaitre_notched_bar.py:21-23`), and A0 §1.1 is that statement in numbers. A collapse load that depends on mechanism *thickness* is not fixed by fixing the energy. |
| **`LadrunoConcrete3D -eta`** (the fork's *shipped* Duvaut–Lions, `SRC/material/nD/LadrunoConcrete3DKernel.h:1492-1503`) | The same blend as this ADR, but applied **inside** the CDPM2 kernel on the *effective* stress and on `κ_p`, before damage is chained. | It is twelve lines that need the committed internal effective stress and `κ_p` — quantities no generic `NDMaterial` seam exposes. It cannot be lifted; it can only be re-written (brief F1). It is also **Tier-1 only** (`!mp.implex`) and inert in the BeamFiber view. And it is a *damage* model: ADR-31 §4.4 refused the nominal-stress blend for exactly the reason this ADR inherits as **D1** — relaxing a nominal stress relaxes damage as well as plasticity. |
| **ADR-32 embedded discontinuity** (`LadrunoDispBeamColumn*`) | Localization in a **frame element**, by a kinematic enhancement (an embedded softening hinge) with per-IP `lch`. | It is a 1-D element-kinematics fix. There is no 3-D continuum version, and building one is ADR-59's descoped territory. |
| **ADR-59 gradient/non-local concrete** | Localization width, correctly, by a true internal length. | **Descoped.** Its kill list (unverified seam claim, non-byte-identical off switch, wrong regularized variable, "free" infrastructure reuse, no named consumer, no reconciliation) is the standard this ADR is written against; §12 lists the artifacts that discharge each item. |

**Conclusion:** nothing shipped in the fork regularizes localization *width* for a plasticity-type
3-D continuum material. This ADR is not "mostly re-wiring" of anything.

## 4. Formulation

### 4.1 The two-track (TT) update

Per material point, per step, with `Δε` the total strain increment and `C_e` the inner material's
**initial** tangent read *before* `setTrialStrain` (decision **D3**):

```
    inner.setTrialStrain(eps_{n+1})                       # the inner runs INVISCIDLY on the
    sigma_bar = inner.getStress()                         # total strain path — it never sees tau
    C_bar     = inner.getTangent()

    sigma_tr  = sigma_vp,n + C_e : delta_eps               # viscoplastic elastic predictor
    beta      = dt / (tau + dt)
    sigma_vp  = (1 - beta) * sigma_tr + beta * sigma_bar
    D_vp      = (1 - beta) * C_e      + beta * C_bar
```

`β ∈ [0,1]`. `β → 1` (τ → 0) is the inviscid return; `β → 0` (τ → ∞ at finite Δt) is the frozen
elastic trial. The tangent blend is the *exact* algorithmic derivative of the stress blend, not an
approximation: `∂σ_vp/∂ε = (1−β)C_e + β ∂σ_bar/∂ε`, verified against a central finite difference
to `3.31e-10` on both the 1-D and the 6×6 J2 model (oracle PV5).

### 4.2 The `τ == 0 or Δt <= 0 ⇒ β = 1` fork convention

Mathematically `Δt → 0` is the *elastic* limit (`β → 0`, frozen trial). Operationally a missing or
zero pseudo-time increment must **never** silently turn a material elastic. The fork therefore
gates to `β = 1`, i.e. **inviscid**, whenever `τ ≤ 0` or `Δt ≤ 0` — the same convention the shipped
`-eta` already uses (`LEDGER_quirks.md:1612`, `SRC/material/nD/LadrunoConcrete3DKernel.h:1486-1488`).
The elastic limit is reached only through a large τ at finite Δt, which is the correct knob.

Consequence, inherited and disclosed: like `-implex`, the wrapper is **inert without a positive
`ops_Dt`**. That is precisely why §3's lane decision exists.

### 4.3 The A0 theorem — TT **is** true Duvaut–Lions on proportional monotonic paths

True Duvaut–Lions (Simo & Hughes §2.7; Simo, Kennedy & Govindjee 1988) projects the
*viscoplastic* trial state `(σ_vp,n + C_e Δε, q_vp,n)` onto the yield surface and relaxes **both**
the stress and the internal variables toward that projection:

```
    (sigma_bar, q_bar) = P( sigma_vp,n + C_e : delta_eps , q_vp,n )
    sigma_vp = (sigma_tr + (dt/tau) sigma_bar) / (1 + dt/tau) = (1-beta) sigma_tr + beta sigma_bar
    q_vp     = (q_n      + (dt/tau) q_bar    ) / (1 + dt/tau) = (1-beta) q_n      + beta q_bar
```

A generic wrapper cannot re-seed `q`, so it cannot do this — or so the planning brief assumed
(F2). It is wrong.

> **Theorem (1-D, from rest, monotonic loading).** For the 1-D associated model with **any**
> hardening function `K(α)` — linear, softening, or nonlinear — the two-track wrapper and true
> Duvaut–Lions produce the **identical** stress path.
>
> *Proof.* The TDL update gives `σ_{n+1} = σ_tr − βE(ᾱ − α_n)` and `α_{n+1} = α_n + β(ᾱ − α_n)`,
> hence `σ_{n+1} + Eα_{n+1} = σ_tr + Eα_n`. So `ψ_n = σ_n + Eα_n` advances by exactly `EΔε` every
> step — elastic steps included, where `ᾱ = α_n`. From rest `ψ_0 = 0`, therefore
> **`α_n = ε_n − σ_n/E`**: the relaxed internal variable is exactly the plastic strain of the
> viscoplastic stress itself. Substituting into the projection equation
> `σ_tr − E(ᾱ − α_n) = K(ᾱ)` gives `Eε_{n+1} − Eᾱ = K(ᾱ)`, i.e. `σ̄ = K(ᾱ)` with
> `ᾱ = ε_{n+1} − σ̄/E` — which is the definition of the **inviscid** stress at the total strain
> `ε_{n+1}`. TDL's projection target therefore *is* `inner.getStress()` on the total strain path,
> which is exactly what the two-track wrapper blends toward. Both return
> `(1−β)σ_tr + β σ_inviscid(ε_{n+1})`. ∎

**Where it holds, measured** (oracle `run_tt_vs_tdl_point`, max relative `|σ_TT − σ_TDL|`):

| case | measured |
|---|---|
| perfect plasticity | `3.20e-14` |
| 1-D **linear**, `H/E ∈ {+0.10, +0.02, 0, −0.02, −0.10}` × `De ∈ {0 … 0.3}` | `9.23e-14` |
| 1-D **exponential** (nonlinear) hardening AND softening, `α_end/α_f ≈ 4.75` | `2.80e-13` |
| J2 (Voigt-6), **proportional** path, `H/E = ±0.05` | `4.72e-15` |
| **the whole A0 bar**, `De ∈ {0, 3e-4, 1e-3, 3e-3}` × `N ∈ {40, 80, 160}` × 2 imperfection conventions — w2 / peak load / work / curve L2 | `1.31e-5` / `1.8e-16` / `1.25e-6` / `4.82e-6` |

**Where it fails, measured** — the wrapper's declared approximation:

| case | De = 0.01 | De = 0.10 |
|---|---|---|
| J2, **non-proportional** (axial 120 steps, then shear 80), `H/E = +0.05` | `4.02e-3` | **`4.42e-2`** |
| J2, **non-proportional**, `H/E = −0.05` | `8.99e-4` | **`2.62e-2`** |
| 1-D linear, load / **unload** / reload, `H/E = +0.10` | `4.63e-2` | `3.14e-1` |
| 1-D linear, load / **unload** / reload, `H/E = −0.02` | `9.23e-2` | **`3.33e-1`** |

The error grows monotonically with De, which is the signature of a genuine constitutive difference
rather than round-off.

**The boundary is over PATHS, not over material classes.** Brief F2 framed the split as
*perfect plasticity vs hardening*; that framing is refuted — hardening, softening and strongly
nonlinear hardening are all in the *exact* regime. The real boundary is
**proportional-and-monotonic** vs **non-proportional-or-unloading**. SANISAND sits on the wrong
side of it (α-tensor, fabric `z`, ψ-driven `M^b`), so the wrapper over SANISAND is an
approximation whose size on the real material is **not yet measured** — the J2 surrogate above is
the only quantification we have. That is a WP-D measurement, not a WP-C blocker (§10 OQ9, §11 R8).

### 4.4 The closed-form 1-D anchor

The steady overstress of the discrete backward-Euler fixed point is `σ* − σ_Y = E·ε̇·τ`, **exactly
and independently of Δt** (the Δt cancels: `(1−β)/β = τ/Δt`). Measured to `4.26e-15` relative over
`Δt ∈ {1, 0.25, 0.05}` for both TT and TDL. This is the ADR's non-self-referential oracle for the
verification manifest (ADR-87), and it is what the C++ byte check will be pinned against.

## 5. Architecture

`SRC/material/nD/LadrunoDuvautLions.{h,cpp}` — a generic `NDMaterial` wrapper, ND class tag
**33022**, modelled on `StagedStrainNDMaterial` (33014, `SRC/classTags.h:587`), which is the
fork's own adopting-wrapper precedent.

**Command (Tcl and Python; flags after positionals, ADR-86 parser rule):**

```
    nDMaterial LadrunoDuvautLions $tag $innerTag -tau $tau
```

| Item | Decision | Why, with the source that forces it |
|---|---|---|
| **Adopting ctor** | Take the inner by tag, `getCopy("ThreeDimensional")` it once at construction, own the copy, delete it in the dtor. | `StagedStrainNDMaterial` pattern. |
| **`setTrialStrain(v)` / `setTrialStrainIncr(v)`** | Keep an **absolute** committed strain; `setTrialStrainIncr` reconstructs `ε_{n+1} = ε_committed + Δε` and calls the same core. Never forward an *increment* to the inner. | The blend needs the total strain path for the inner AND the increment for the predictor; only an absolute state gives both. The strain-rate overloads discard the rate upstream (`ManzariDafalias3D.cpp:88-91`) so no rate can be recovered from the inner. |
| **`getCopy(const char *type)`** | Resolve **every supported type explicitly** — `"ThreeDimensional"`, `"3D"`, `"PlaneStrain"`, `"PlaneStrain2D"` — build the inner's view with `inner->getCopy(type)`, wrap it, and **never call the inner's `getCopy(void)`**. Return `0` gracefully on an unsupported type (no `exit(-1)`). | `StagedStrainNDMaterial::getCopy(const char*)` (`SRC/material/nD/StagedStrainNDMaterial.cpp:293-309`) has **no** `"ThreeDimensional"` special case — it always forwards the string; only `InitStrainNDMaterial.cpp:290-292` special-cases it to `getCopy(void)`. `ManzariDafalias::getCopy(void)` is `exit(-1)` (`ManzariDafalias.cpp:447-465`), and only the four strings above exist (`:534-555`). An implicit route is a process abort. |
| *(do not lean on this)* | `FluidSolidPorousMaterial::getCopy(const char *code)` (`SRC/material/nD/soil/FluidSolidPorousMaterial.cpp:415-425`) **ignores** its `code` argument and routes through the copy ctor, whose `:156` calls the soil's `getCopy(void)`. It is the only plane-strain route in the tree that reaches a wrapped material's void `getCopy` — a **coincidence**, recorded as a quirk (§12), never a dependency of this wrapper or its tests. |
| **Copy, never alias** | Deep-copy every `Vector`/`Matrix` returned by the inner into the wrapper's own members before using or returning it. | `ManzariDafalias3D::getStress/getStrain/getTangent` return **static class buffers** (`ManzariDafalias3D.h:78-79`). Aliasing them makes a second material instance silently overwrite the first. |
| **`setParameter`** | Handle `"tau"` locally; on a **tag miss**, forward `(argv, argc, param)` to the inner and return its result. | `updateMaterialStage` resolves by string-tag match *inside* `ManzariDafalias` (`ManzariDafalias.cpp:791-828`); without forwarding, the static `mElastFlag` never flips and a staged geostatic run silently stays elastic. Precedent: `FluidSolidPorousMaterial.cpp:320-333`, `InitStrainNDMaterial.cpp:436-473`, and `StagedStrainNDMaterial.cpp:397-399` (which forwards unconditionally — this wrapper must **not**, because it owns `"tau"`). |
| **τ as a `Parameter`** | `setParameter("tau")` / `updateParameter`, so a De sweep or a staged activation does not need reconstruction. | Decision **D6**; the shipped `-eta` cannot do this. |
| **`revertToLastCommit`** | Roll the wrapper's own `(σ_vp, ε)` back to the committed values and forward to the inner — **documenting that the inner's is an empty stub** (`ManzariDafalias.cpp:513-517`). | Pre-existing limitation, not a new one. It must be in the guide, not discovered at a step cutback. |
| **`getResponse` tokens** | `overstress` (‖σ_vp − σ_bar‖), `beta`, `dt`. | These are what makes De reportable per run (S4). |
| **Non-uniform-Δt diagnostic** | Track the committed Δt; warn **once** when it varies by more than a declared threshold between commits, naming the integrator. | §3 item 2 / brief F3's converse: a rate artifact when τ is comparable to a static pseudo-step. This is a new quirk the ADR owes (§12). |
| **`sendSelf` / `recvSelf`** | `ID(4)` = `{tag, inner classTag, inner dbTag, nstate}` → `Vector` of the wrapper state (`τ`, committed `σ_vp`, committed `ε`, last Δt) → **delegate** to `inner->sendSelf`. On receive, `theBroker.getNewNDMaterial(ID(1))` if the inner is null, then `setDbTag(ID(2))`. | Verbatim the `StagedStrainNDMaterial.cpp:316-…` protocol, which is the tested one. |
| **Registration** | ×3 (Tcl `TclNDMaterialCommands`, `OpenSeesNDMaterialCommands.cpp` functionMap, `FEM_ObjectBroker`) + CMake. | The Tcl/Python double-dispatch trap is the first entry in `LEDGER_quirks.md`. |
| **Scope guard** | Refuse (graceful, with a message) an inner whose `getCopy("ThreeDimensional")` fails; **do not** attempt to detect damage models — document them as out of scope instead. | Decision **D1**; damage cannot be detected at the seam. |

## 6. Claims, each with its gate

| # | Claim | Gate | Class |
|---|---|---|---|
| **C1** | τ = 0 is **byte-identical** to the inner material, including the instruction path | PV1 clone at the material point and at the element; g++ byte check vs the oracle fixture. Oracle already at `max|Δσ| = 0.0` **exactly** for both variants, both point models, `H < 0 / = 0 / > 0`, loading and unloading. | correctness |
| **C2** | 1-D steady overstress `= E·ε̇·τ`, independent of Δt | PV3 clone. Oracle: `4.26e-15` over `Δt ∈ {1, 0.25, 0.05}`. | correctness |
| **C3** | Tangent `= (1−β)C_e + β C_ep`, matches a central FD | PV4 / PV5 clones. Oracle: blend identity `0.0`; FD `3.31e-10`. | correctness |
| **C4** | **Band width converges under h-refinement at fixed De** (positive) **and does not at τ = 0** (negative control) | The new width gate, §7. Oracle A0: w2 ratios `0.972 / 1.010 / 1.004 / 1.001` at the finest pair for `De ≥ 3e-4`, against exactly `0.500` at τ = 0. | **the deliverable** |
| **C5** | Collapse load at fixed De converges in h (the three-mesh band collapses) | Slow-tier gate, capacity three-clause rule. | **the deliverable** |
| **C6** | De-dependence of width and `q_u` is **measured and disclosed** | vault-14-style τ × T sweep with matched pairs at **different step counts** (§7.4). | disclosure |
| **C7** | The wrapper is transparent to `updateMaterialStage`, `getCopy(type)`, database round-trip and the MP wire | Zone-A + roundtrip + a `test_fspm_over_manzari_family.py` twin; mutation gate `-DLADRUNO_MUTATE_DUVAUTLIONS` must turn the suite red. | wiring |

**Claims this ADR explicitly does NOT make.** It does not restore uniqueness or the bound theorems
(S2). It does not introduce an intrinsic length in statics (§3, brief F4 — corroborated by A0's
factor-12-per-factor-3 De sensitivity). It does not claim mesh independence without a declared De.
It does not apply to damage or plastic-damage inner materials (D1). It does not change the
tetrahedron prohibition (S3). And after §4.3 it does not claim to be true Duvaut–Lions on
**non-proportional or unloading** local strain paths — there it is a declared approximation of
measured size on a J2 surrogate and unmeasured size on SANISAND.

## 7. Acceptance case

The shape is fixed by de Borst 1993 Fig. 18 and vault 65 D5; the *numbers* are TIMs' to set (OQ2).

| Leg | Deck | Material | Element | Status |
|---|---|---|---|---|
| **A0** | 1-D softening bar, `N = 20/40/80/160` | numpy oracle, TT vs TDL | — | **DONE** — `_adr90_a0_results.md`, commits `cc7c7f7a5` / `8863a468c` |
| **A1** | Plane-strain biaxial, symmetric half, smooth ends, imperfection at mid-height; 3 meshes + 2 orientations | `DruckerPragerPlaneStrain` (non-associated, no state dependence — a clean bifurcation) | `quad PlaneStrain` | P2 |
| **A2** | Same deck | `LadrunoSANISANDPlaneStrain` (33021) | `quad PlaneStrain` | P2, with the OQ7 constraint below |
| **A3** | Same specimen extruded one element thick | `LadrunoSANISAND` (33019/33020) | **`LadrunoBrick -formulation bbar`** | P2 — **the only admissible leg** (S3, D7) |
| **A4** | R3 Prandtl–Reissner strip (existing gate) | SANISAND / DP | B-bar hex | P2 regression at the declared De |

### 7.1 A0's four corrections to the acceptance protocol

**(i) The imperfection must be a mesh-convergent FIELD, not de Borst's one element.** A0 §5.2 ran
both. With a one-element defect (whose physical size shrinks with N) the width **never** converges
at any De — the peak stays one element wide (`w3 = h` at every mesh) and `w2` drifts toward the
specimen length — although the load–displacement *curve* does converge for `De ≥ 1e-3`
(`curveL2` 0.116 → 0.027 → 0.001). With a graded notch over a **fixed physical length** the width
converges: `w2` ratios 0.972 / 1.010 / 1.004 / 1.001. The reason is structural, not numerical: an
imperfection field that does not converge under refinement cannot produce a converging width, and
a quasi-static Duvaut–Lions bar has no intrinsic length of its own to supply one.

> **And the notch must be graded, not flat.** A *flat* fixed-length weak zone (the brief's
> convention (b) read literally: 2/4/8/16 elements all at 0.9 σ_Y) makes the continuum solution
> piecewise constant with the zone boundary on a mesh line, so **every mesh represents it exactly**
> and the result is bit-identical in N at every De > 0. That is a gate that cannot fail. Use a
> parabolic notch with a tie-breaking offset.

**(ii) The De window.** The brief proposed `De ∈ {0.01, 0.1, 1}`. On the A0 bar deck all three are
far past the point where the specimen still localizes:

| De | N | P_peak | ε_p max | inelastic ratio | w2 | w3 |
|---|---|---|---|---|---|---|
| 0.01 | 20 / 160 | 23.79 / 23.80 | 0.0263 / 0.0289 | 31.1 / 31.1 | 97.9 / 98.1 | 100 / 100 |
| 0.10 | 20 / 160 | 59.25 / 59.26 | 0.0177 / 0.0180 | 11.1 / 11.1 | 99.8 / 99.8 | 100 / 100 |
| 1.00 | 20 / 160 | 264.8 / 264.8 | 0.0068 / 0.0069 | 0.76 / 0.76 | 0 | 0 |

At 0.01 and 0.1 the "band" is the **whole 100 mm bar**. At **De = 1 the load never peaks at all** —
it rises monotonically to 265 N, 14.7× the rate-independent peak of 18 N — so there is no post-peak
branch and "the width converges" is a statement about an essentially elastic bar. **On the bar deck
the regularized-and-still-localizing window is `De ≈ 3e-4 … 1e-3`**, two orders of magnitude below
the brief's guess. The window is deck-dependent (it scales with the softening ductility and the
loading rate), so **the SANISAND-deck window is to be MEASURED in P2, not inherited from A0.** A0
fixes the *method* for finding it: sweep De down until the response both localizes (a post-peak
branch exists, `inelastic ratio > 1`) and converges.

**(iii) The τ = 0 negative control must FAIL, and must use a unique weakest point.** A0 §5.1: with
either a one-element defect or a graded notch, `w1 = w2 = w3 = h` **exactly** at all four meshes
and `W₅₀` halves with `h` (ratios 2.000 / 2.000 / 2.000). With a flat weak zone the tie lets
round-off pick the sub-band and `w2/h` wanders 2.0 → 14.0 with no pattern — the control becomes
uninterpretable rather than failing cleanly.

**(iv) The De collapse must be run at DIFFERENT step counts or it is tautological.** Because
`β = Δt/(τ+Δt) = 1/(1 + nsteps·De)` and the strain increment is `u_max/nsteps`, two runs with the
same `(De, nsteps)` are **bit-identical** whatever `(τ, T)` produced that De — measured to 12
significant digits. Run the matched pairs with different step counts; the collapse then holds to
**0.47 % in width, 0.009 % in peak load, 0.025 % in work**, and that residue is the Δt-transient
effect of §3 item 3, not a De violation.

### 7.2 Measurements per leg

Band width by the **threshold-free** metric (below), peak load, post-peak curve, dissipated work,
`De = τ/T`, **the step count**, and `ladrunoBuild()`.

### 7.3 The width metric (OQ5, settled)

`w2 = √(12·Var)` with `Var = Σ p_e[(x_e − x̄)² + h²/12] / Σ p_e` over the **post-peak plastic-strain
increment** profile. The `h²/12` term is the within-element variance of the piecewise-constant
profile; it makes a one-element band read **exactly h** and a k-element top hat read **exactly
k·h**, so the number is comparable across meshes. Without it a one-element band reads 0. Pinned by
`tests/test_duvaut_lions_oracle.py::test_band_width_metric_is_calibrated`.

Threshold-based metrics are **not** admissible as the headline: on the same A0 runs the
threshold-count `w1` and the FWHM `w3` disagree by up to a factor 40 (risk R4 confirmed).

### 7.4 Gates

τ = 0 → width ∝ h (the control must fail objectivity). τ > 0 at fixed De → `w(h/4)/w(h/2)` within
the band TIMs declares, **and** the load–displacement curves converge, **reported separately**.
De × {½, 1, 2} at different step counts → width monotone in De and matched pairs collapse (C6).

### 7.5 Parameters TIMs must supply (OQ2)

Target band width relative to B; ramp duration / strain rate; the De they will run at (informed by
the P2 measurement of the SANISAND-deck window, not by A0's bar numbers); tolerance bands; and the
ultimate criterion (OQ1, Prof. Gorini).

## 8. Phases and exit gates

| Phase | Content | Exit gate | Warrant items |
|---|---|---|---|
| **P0 — ADR + numpy oracle** ✅ **COMPLETE 2026-09-04** | This ADR; `tests/_testbed/duvaut_lions_ref.py`; `tests/_testbed/run_a0_sweep.py`; `tests/test_duvaut_lions_oracle.py` (9 `zone_a` cases, **22.6 s measured**, numpy only — imports no OpenSees); `_adr90_a0_results.md`. | PV1–PV6 green for **both** architectures over **both** point models; A0 width ∝ h at τ = 0 and convergent at fixed De; TT-vs-TDL quantified. **All met.** Commits `cc7c7f7a5`, `8863a468c` (PR #783). | none (docs + `tests/_testbed`) |
| **P1 — the C++ wrapper** | `LadrunoDuvautLions.{h,cpp}` per §5; Tcl + Python; `classTags.h` 33022 with the `// N. Mora-Bowen (Ladruno) — …` comment. | g++ byte check against the oracle fixture; Zone-A: byte gate (material point **and** `LadrunoBrick`), PV3 overstress, non-tautology guard (viscous ≠ inviscid by > 1e-3), stage-forwarding over `LadrunoSANISAND`, database round-trip, MP wire, **mutation gate red**. | manifest `ledger/ladruno_duvaut_lions.yml`, `-DLADRUNO_MUTATE_DUVAUTLIONS`, guide stub, ledger rows, the three new quirks entries (§12) |
| **P2 — acceptance case** | Legs A1 → A3 in `tests/` (slow tier) with the G5-style positive + negative pairing; the De sweep; **the SANISAND-deck De window measured**; the §3 lane decision (a)/(b) made on measurement; **the non-proportionality error leg** (§10 OQ9). | C4, C5, C6 green; A4 regression unchanged within the De family; the non-proportionality error reported against the campaign's own three-mesh band. | slow-tier docstring wall times; capacity three-clause rule |
| **P3 — TIMs integration** | Provenance fields (τ, De, rate, step count) in the campaign's output schema; one SFIM/APE radial probe with declared τ on the chosen lane; out-of-family verdict. | Joint sign-off; the vault note records the acceptance case as **verified**, not claimed. | `reviews/handoff_adr90.md` → verdict |
| **P4 — optional** | Explicit/transient tier (`CentralDifferenceLadruno`); `-implex` interplay; mesh-aware τ via the `getCharacteristicLength()` latch (`LadrunoConcrete3D.cpp:347-350`). | Only on a named consumer. | — |
| **WP-F — PARKED** | Material-specific `-tau` **inside** `LadrunoSANISAND` (true DL on σ, α, z using the `protected` base members). | **Trigger:** fires only if A2/A3 measure the non-proportional wrapper error **above the campaign's own three-mesh band**. A0 gives no reason to fire it now. | same warrant as P1 |

## 9. Decisions

| # | Decision | Status | Rationale |
|---|---|---|---|
| **D1** | **Scope: plasticity-type inner materials only.** Damage / plastic-damage models are documented as out of scope; not guarded, because damage cannot be detected at the seam. | **DECIDED** | Brief F2; ADR-31 §4.4 precedent (a nominal-stress blend relaxes *damage* as well as plasticity). |
| **D2** | **Relaxation form: the generic two-track blend.** The wrapper's **declared validity domain is proportional-and-monotonic local strain paths**, where §4.3 proves and measures it to be *exactly* true Duvaut–Lions. Outside it — non-proportional or unloading — it is a declared approximation, measured at 9.0e-4 … 4.4e-2 (J2, De 0.01 → 0.10) and 3.3e-1 (1-D unloading). | **DECIDED 2026-09-04 (CP1, owner)** | **Correction to brief F2:** the boundary is *proportional-and-monotonic vs non-proportional-or-unloading*, **not** *perfect-vs-hardening*. Linear, softening and strongly nonlinear hardening are all in the exact regime (`≤ 2.8e-13`). |
| **D3** | `C_e` for the trial predictor = `inner.getInitialTangent()` read **before** `setTrialStrain`. | **DECIDED** | SANISAND rewrites its initial tangent every step (`ManzariDafalias3D.cpp:146-149`); the predictor must use the committed-state `C_e`. |
| **D4** | Δt source `ops_Dt`; `Δt ≤ 0 ⇒ β = 1` (inviscid); **plus a non-uniform-Δt diagnostic** that warns once when Δt varies by more than a declared fraction between commits. | **DECIDED** | Brief F3 and its converse; `LEDGER_quirks.md:1612`. A0 §5.5 shows the transient is Δt-sensitive at fixed De (2.06 % in w2 over 32×), so the warning is not cosmetic. |
| **D5** | Calibration lane: §3 **(a)** transient with β_K = 0 primary; **(b)** uniform-pseudo-time static fallback; the two gated equal at matched De. | **DECIDED (to be confirmed by P2 measurement)** | Brief F3/F4; vault 14. |
| **D6** | τ is a `Parameter` (`setParameter("tau")`), so De sweeps and staged activation need no reconstruction. | **DECIDED** | S5; the shipped `-eta` cannot do this. |
| **D7** | The element for the final acceptance leg is the **B-bar hex** (A3). No result on another element is admissible for the campaign. | **DECIDED** | S3. |
| **D8** | Command shape `nDMaterial LadrunoDuvautLions $tag $innerTag -tau $tau` — flags after positionals. | **DECIDED** | ADR-86 parser rule; apeGmsh emitter guide. |
| **D9** | Off-switch: τ = 0 byte-identical **including the instruction path** (ADR-31's "same instructions" rule). | **DECIDED** | ADR-59 kill-list item B2. Oracle already at exactly `0.0`. |
| **D10** | Class/command name `LadrunoDuvautLions` (was OQ8). A Perzyna variant would be a different class, not a flag. | **DECIDED** | Says what it is; leaves the Perzyna name free. |

## 10. Open questions

- **OQ1** [Prof. Gorini] — the ultimate criterion. The acceptance case can run without it; its `q_u`
  gate cannot be *cited* until it lands.
- **OQ2** [TIMs] — the acceptance-case numbers (§7.5). The fork drafts; TIMs sets.
- **OQ3** [fork] — **CLOSED by A0, 2026-09-04.** Two-track is adequate; §4.3 says exactly where and
  why, and the answer is a theorem rather than a tolerance.
- **OQ4** [both, P2] — lane (a) vs (b) on the actual radial probe; and whether vault 14's β_K family
  and the wrapper's τ family collapse onto one De when both are present. They should (same
  Kelvin–Voigt structure) — but that is a claim to test, not to assume.
- **OQ5** [fork] — **CLOSED by A0.** The threshold-free width metric is §7.3, unit-pinned.
- **OQ6** [fork, P2] — interaction with SANISAND's substepping / `-Pmin` whole-tensor resets. The
  inner is called once per trial and is unaware of the wrapper, so nominally none — verify on A2.
- **OQ7** [fork] — **REWRITTEN after WP-B (PR #785).** The `ManzariDafaliasPlaneStrain` null-ctor
  classTag anomaly (`LEDGER_quirks.md:4826`) is a **recorded owner decision NOT to fix**, not an
  outstanding bug. Consequence for this ADR: **leg A2 must not depend on a broker / database
  restore of the plane-strain material.** Build A2's plane-strain models directly in the deck; keep
  the round-trip coverage (C7) on the 3-D view, where the classTag is correct. `getCopy(void)`
  coverage for `LadrunoSANISANDPlaneStrain` landed in WP-B.
- **OQ8** — **CLOSED**, see D10.
- **OQ9** [fork, P2 — NEW] — **how large is the two-track approximation on SANISAND itself?** §4.3
  measured a J2 surrogate (4.4 % at De = 0.10 on a two-leg non-proportional path). The real
  material's radial pushover under a footing rotates its stress path continuously. This is the
  non-proportionality error leg added to WP-D at CP1, and it is WP-F's trigger.

## 11. Risks

| # | Sev | Risk | Mitigation / status |
|---|---|---|---|
| **R1** | ~~BLOCKING~~ → **RETIRED 2026-09-04** | *"Two-track ≠ Duvaut–Lions for hardening materials; the oracle passes on perfect plasticity and the claim silently over-reaches."* | **Retired by the §4.3 theorem**, which shows hardening is not the discriminator at all, and by the measurements that bracket it (`≤ 2.8e-13` inside the domain, `4.4e-2 … 3.3e-1` outside). The residual concern is re-issued as **R8**, which is a different risk with a different mitigation. |
| **R2** | BLOCKING | Regularization inert or erratic on the campaign's static lane (F3/F4); "mesh-converged" claimed from a lane where De is undefined. | §3 lane decision + the D4 Δt-uniformity diagnostic + De **and step count** reported per run. |
| **R3** | MAJOR | Band width is De-dependent; a reviewer reads "converges in h" as "mesh-independent". | §3.2 claim wording; C6 disclosure is mandatory. A0 quantifies the exposure: **factor 12 in width for factor 3 in De**. |
| **R4** | MAJOR | Width metric threshold-dependent and reverses across meshes. | **CONFIRMED as real** — `w1` and `w3` disagree by up to 40× on the same A0 run. Mitigated by quoting the threshold-free `w2` (§7.3), unit-pinned. |
| **R5** | MAJOR | A 2-D result on `quad` does not transfer to the B-bar hex. | A3 is the admissible leg; A1/A2 are ladders, not deliverables (D7). |
| **R6** | MAJOR | Static-buffer aliasing / implicit `getCopy(void)` / stage-flag forwarding — each a silent wrong answer. | §5 makes each an explicit design rule with the source that forces it; each becomes a Zone-A test in P1, plus the `test_fspm_over_manzari_family` twin. |
| **R7** | MINOR | Slow-tier cost (the R3 Prandtl gate alone is 61 min; A3 × 3 meshes × {τ=0, τ>0} × 3 De ≈ a day of compute). | Wall times in docstrings; nightly Zone-B for the sweep, Zone-A slow for one leg. |
| **R8** | **MAJOR — NEW** | **The wrapper's error on rotating stress paths is unmeasured until WP-D.** A footing pushover rotates the principal directions continuously, which is exactly the regime §4.3's theorem excludes. The J2 surrogate says 0.09 % at De = 0.01 and 4.4 % at De = 0.10, growing with De — SANISAND, with an α-tensor and fabric `z` that the wrapper cannot re-seed, could be worse. Shipping the wrapper with a *proved* correctness statement makes it easy to forget that the proof does not cover the consumer's own load path. | The declared validity domain is written into D2, §6's NOT-claimed list and the guide. **OQ9** is a named P2 leg with a number attached. **WP-F is parked with an explicit trigger** rather than cancelled: it fires if A2/A3 measure the error above the campaign's own three-mesh band. |

## 12. Ledger obligations

**Discharged in this PR (#783):**

- `LEDGER_implementations.md` — the ADR 90 row: `ND_TAG_LadrunoDuvautLions` **33022 — RESERVED, not
  yet built**, status *P0 complete — oracle + A0; P1 not started*.
- `Ladruno_implementation/README.md` — the index line for 90.
- `LEDGER_quirks.md` — the three A0 entries: the two tautological gates (fixed-`nsteps` De collapse;
  a flat fixed-length imperfection); the residual-hardening default that silently made a "perfectly
  plastic" bar harden; and Newton's need for a consistent-tangent predictor post-peak at N ≥ 40.
- A dated **CORRECTION** callout in `_adr90_regularization_planning_brief.md` at §2 F2, §7 D2 and
  §8 OQ7 — appended, not rewritten.

**Owed by P1:**

- `SRC/classTags.h` 33022 with the `// N. Mora-Bowen (Ladruno) — …` comment and the per-registry
  non-collision note (EigenSOE `:57`, PATTERN `:646`).
- `Ladruno_implementation/ledger/ladruno_duvaut_lions.yml` (verification manifest, ≥ 1
  non-self-referential oracle — §4.4's closed form is it).
- `-DLADRUNO_MUTATE_DUVAUTLIONS` mutation build flag; the family suite must turn red.
- `LadrunoDuvautLions_guide.md` stub, both interpreters, `banner_features.txt` line **on ship**.
- Two further quirks entries: (i) the **rate artifact** when τ is comparable to a static pseudo-step
  (brief F3's converse — the ledger has no entry for it); (ii) **`FluidSolidPorousMaterial::getCopy(const char*)`
  ignores its `code` argument** and routes through the copy ctor, making it the only plane-strain
  path in the tree that reaches a wrapped material's `getCopy(void)` — a coincidence, not a
  contract.
- The **F8 stale-doc corrections** (`LadrunoConcrete3D_guide.md:396,480,593` still list `-eta` as
  "not yet (P3)"/"deferred"; `31_ladruno_concrete3d_adr.md:91` claims `-eta` works "under any tier"
  while the kernel gates it to Tier-1) — landed by WP-B in PR #785; verify on merge.

## 13. Implementation log

| Date | Event |
|---|---|
| 2026-09-04 | **Planning brief** written (`_adr90_regularization_planning_brief.md`) from three read-only source sweeps: the shipped `-eta` kernel, the SANISAND attachment surface and its five traps, and the fork's ADR acceptance bar. Number 90 allocated (88 taken by PR #778, 89 proposed for Track T); ND tag 33022 reserved. Commit `3feec0fc9`. |
| 2026-09-04 | **WP-A / A0 executed.** `tests/_testbed/duvaut_lions_ref.py` + `run_a0_sweep.py` + `tests/test_duvaut_lions_oracle.py` (9 `zone_a` cases, 22.6 s). PV1–PV6 green for both architectures over both point models. **The TT ≡ TDL theorem found and proved**, with the boundary measured. A0 bar: H1 confirmed (width ≡ h exactly, work halves), H2 confirmed for `De ≥ 3e-4` with a mesh-convergent imperfection field and **refuted** for a one-element defect, H3 confirmed but tautological in its naive form, H4 → TT adequate. Commits `cc7c7f7a5`, `8863a468c`. |
| 2026-09-04 | **Checkpoint CP1 — owner decision.** **D2 approved as recommended**: build the generic two-track wrapper; declared validity domain = proportional-and-monotonic local strain paths; a non-proportionality error leg (OQ9) added to WP-D on the SANISAND deck; **WP-F parked** with the trigger *"fires only if A2/A3 measure the non-proportional error above the campaign's own three-mesh band"*. This ADR written against that decision. |
| 2026-09-04 | WP-B (PR #785) landed the SANISAND plane-strain prerequisites and the `-eta` doc corrections, and established that the plane-strain classTag anomaly is a recorded owner decision **not** to fix — which rewrote **OQ7**. |

## References

- **de Borst, R., Sluys, L. J., Mühlhaus, H.-B. & Pamin, J. (1993).** "Fundamental issues in finite
  element analyses of localization of deformation." *Engineering Computations* **10**, 99–121.
  Seafile `05_numerics_localization/deborst1993.pdf`. The comparison of regularization strategies
  (rate dependence, gradient, Cosserat, non-local), the canonical plane-strain biaxial acceptance
  deck (smooth ends, symmetric half, an imperfect element near mid-height, non-associated
  Drucker–Prager — their Fig. 18 with Duvaut–Lions), and the internal length `ℓ = 2 m c_e / E`.
  **The caveat this ADR is built around**, paraphrased: rate dependence works for both failure
  mechanisms, but its applicability is limited to transient loading and the regularizing effect
  falls away rapidly for slow loading or near the rate-independent limit. A0 §5.2 is that sentence
  in numbers.
- **Simo, J. C. & Hughes, T. J. R. (1998).** *Computational Inelasticity.* Springer. §2.7 —
  viscoplastic regularization; the Duvaut–Lions backward-Euler closed form and the relaxation of
  the internal variables that §4.3 reformulates. Box 3.2 — the J2 algorithmic tangent used by the
  oracle's 3-D point model.
- **Simo, J. C., Kennedy, J. G. & Govindjee, S. (1988).** "Non-smooth multisurface plasticity and
  viscoplasticity. Loading/unloading conditions and numerical algorithms." *IJNME* **26**,
  2161–2185. The closest-point-projection formulation of Duvaut–Lions this ADR implements.
- **Perzyna, P. (1966).** "Fundamental problems in viscoplasticity." *Advances in Applied
  Mechanics* **9**, 243–377. The alternative overstress family; a Perzyna variant would be a
  different class (D10), not a flag on this one.
- **Rudnicki, J. W. & Rice, J. R. (1975).** "Conditions for the localization of deformation in
  pressure-sensitive dilatant materials." *JMPS* **23**, 371–394. Why a **non-associated**
  pressure-sensitive material localizes while still on a hardening branch — the reason the named
  consumer's collapse loads are mesh-sensitive in the first place.
- **Needleman, A. (1988).** "Material rate dependence and mesh sensitivity in localization
  problems." *CMAME* **67**, 69–85. The original demonstration that rate dependence restores
  well-posedness and sets a length through the wave speed — and, read carefully, the origin of the
  quasi-static caveat that §3 turns into a lane decision.
