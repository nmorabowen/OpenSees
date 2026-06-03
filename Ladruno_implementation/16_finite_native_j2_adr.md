---
title: ADR — Finite-strain-NATIVE combined-hardening J2 with a co-rotating backstress (§14.11 "v2")
project: Ladruno
status: accepted
priority: medium
owner: nmora
tags:
  - adr
  - material
  - finite-strain
  - plasticity
---

# ADR — Finite-strain-native combined-hardening J2 (co-rotating backstress)

> **Primary source.** de Souza Neto, Perić & Owen, *Computational Methods for
> Plasticity* (2008), **§14.11** (finite-strain kinematic hardening) on top of the
> Ch. 14 multiplicative / logarithmic-strain framework already used by
> `LogStrainNDMaterial`. Companion docs: [[10_ladruno_j2_plasticity]] (the 3D
> material + kernel), [[09_finite_strain_material_wrapper]] (the LogStrain v1 path).

## Status

**Accepted.** The numerical algorithm is **prototyped and verified in numpy**
(`tests/ladrunoj2_finite_native_reference.py`, 5/5 green) — no C++ yet.

**Locked (2026-06-02 scoping review):**
- **Class name `LadrunoJ2Finite`**, **classTag `ND_TAG_LadrunoJ2Finite = 33012`**
  (next free ND tag after `LadrunoJ2 = 33011`).
- **Minimal scope** — one focused PR that lands the native material, co-rotation,
  registration, the **consistent tangent (incl. the co-rotation term, see below)**,
  and **flips the strict xfail**. No backlog bundling (IMPL-EX, `LadrunoJ21D`,
  tabulated iso stay out).
- Review confirmed all C++ reuse dependencies exist on disk: `LadrunoJ2Kernel.h`,
  `LogStrainKernel.h`, `LadrunoHardening.h`, `FiniteStrainNDMaterial`,
  `LadrunoBrick -geom finite`, and a polar primitive in `SolidTransformationCorot`.

## Context — the §14.11 boundary v1 cannot cross

`LadrunoBrick -geom finite → LogStrain → LadrunoJ2` (shipped, PR #97) gives
large-strain combined-hardening J2 by wrapping the **unchanged** small-strain
material. It is **exact for the isotropic spine** but **NOT frame-indifferent for
the kinematic part under large rotation** — pinned as a strict xfail
(`tests/test_ladrunoJ2_finite.py::test_finite_combined_hardening_large_rotation_objectivity_is_v2`).

**Why (the mechanism).** The log-strain (MATISU) algorithm co-rotates the elastic
state automatically, `bᵉᵗʳ = f_Δ bᵉ_n f_Δᵀ`, so the deviatoric stress `s` emerges
in the **current** frame.
- **Isotropic** yield uses only `‖s‖` and the scalar `ε̄ᵖ` — both rotation-invariant
  (`‖Q s Qᵀ‖ = ‖s‖`) — so it is exact.
- **Kinematic** yield uses the relative stress `ξ = s − α`. The backstress `α` is a
  *second-order tensor* (the centre of the yield surface). In v1 it lives inside the
  inner material in a **fixed frame**, and the wrapper has **no channel** to rotate
  it (`LadrunoJ2` exposes `α` read-only via `getResponse`; there is no setter). So
  `‖s − α‖` subtracts two tensors in different frames ⇒ wrong yield, wrong flow.

This is a **framework** limit, not a wiring bug: the direct Box-14.4 chain shows the
identical failure (verified). It is exactly the kinematic-hardening-at-finite-strain
case dSNPO defers to §14.11. The 33-agent PR-#97 review independently confirmed it is
fundamental to the "unchanged inner" v1 architecture.

**This is an nD / continuum concern only.** For the 1D `LadrunoUniaxialJ2` the
backstress is a *scalar* `X` along the bar axis — there is no frame to rotate, so 1D
is immune. (This is the precise meaning of "kernel-sharing is nD↔finite-strain only,
not 1D.")

## Decision

Build a **finite-strain-native** J2: a `FiniteStrainNDMaterial` subclass that **owns
all the finite-strain state itself** — `bᵉ_n`, `ε̄ᵖ_n`, **and the backstresses
`α_{n,k}` as SPATIAL tensors** — and **co-rotates the backstress** by the incremental
material rotation each step before running the return map:

```
α̃_{n,k} = R_Δ · α_{n,k} · R_Δᵀ ,   R_Δ = polar(f_Δ)
```

It **reuses the verified return map (`LadrunoJ2Kernel.h`) VERBATIM** — the kernel is
frame-agnostic, so the only genuinely new code is the α push-forward plus carrying /
serializing `α` as spatial state. (This is the payoff of having extracted the kernel
in PR #97.) The 6×6 form is `LogStrainNDMaterial`-style; the element side is the
already-shipped `LadrunoBrick -geom finite`.

### Why a native material and not a wrapper patch

`LogStrain`'s contract is "wrap an UNCHANGED small-strain `NDMaterial` as a black
box." It structurally cannot reach the inner's `α`. A native material that owns the
state is the only place the co-rotation can live.

## Verified algorithm (numpy oracle — `ladrunoj2_finite_native_reference.py`)

```
state (per GP): F_n, Be_n, ebarP_n, alpha_n[N]   # α carried as SPATIAL tensors

setTrialF(F):
    f_Δ    = F · inv(F_n)
    Be_tr  = f_Δ · Be_n · f_Δᵀ                     # elastic predictor  (reuse LogStrainKernel)
    eps_tr = ½ ln(Be_tr)                           #   (degeneracy branches already handled)
    R_Δ    = polar(f_Δ)                            # incremental material rotation
    α̃_k    = R_Δ · alpha_n[k] · R_Δᵀ               # ← THE NEW STEP (§14.11 push-forward)
    (τ, Δεᵖ, alpha_{n+1}, ebarP_{n+1}, Δγ, D)
           = LadrunoJ2Kernel::returnMap(eps_tr as elastic-trial, α̃, ebarP_n)
    σ      = τ / J ;   spatial tangent via LogStrainKernel::spatial_tangent_full
    stage  Be_{n+1} = exp[2(eps_tr − Δεᵖ)] , alpha_{n+1} (CURRENT frame), ebarP_{n+1}
commit: promote the staged state.
```

`R_Δ = polar(f_Δ)` is the incremental, **incrementally-objective** material rotation
(Hughes–Winget style). For a pure rigid increment `f_Δ = Q` it gives `α̃ = Q α Qᵀ`,
exactly cancelling the rotation of `s`; as `f_Δ → I`, `R_Δ → I`, so it reduces to the
small-strain / non-rotating case exactly.

**Prototype verification (numpy, vs the v1 fixed-frame control):**

| Property | co-rotating α | fixed-frame α (v1) |
|---|---|---|
| superposed rigid `Q` on plastic state → `σ_rot = Q σ Qᵀ` | **3.2e-13** ✓ | 6.17 ✗ |
| continuously rotating+stretching plastic path → `σ_rot = R σ_body Rᵀ` | **7.9e-13** ✓ | 2.37 ✗ |
| no-rotation (symmetric `F`) path → reduces to v1 | 2.4e-12 ✓ | — |

So the polar push gives **full path objectivity** (not merely superposed-rotation
objectivity) and the correct small-strain reduction, to ~1e-12.

## Reuse — what already exists

- **Return map:** `SRC/material/nD/LadrunoJ2Kernel.h` (`ladruno_j2_kernel::returnMap`).
- **Hencky kinematics + spatial tangent + degeneracy branches:**
  `SRC/material/nD/LogStrainKernel.h` (`trial_Be`, `hencky_voigt`,
  `spatial_tangent_full`, `be_from_hencky_voigt`).
- **Polar decomposition (Higham/SVD):** already implemented in
  `SRC/element/solidTransformation/SolidTransformationCorot.{cpp,h}` — extract the
  3×3 polar-rotation primitive (or duplicate the few lines) for `R_Δ`.
- **Base + element seam:** `FiniteStrainNDMaterial` (`setTrialF`,
  `getSpatialTangentTensor`) and `LadrunoBrick -geom finite` — both shipped.

## Where (files, when built)

- `SRC/material/nD/LadrunoJ2Finite.{h,cpp}` — the `FiniteStrainNDMaterial`
  subclass; `#define ND_TAG_LadrunoJ2Finite 33012` in `classTags.h` (next free ND
  tag after `LadrunoJ2 = 33011`).
- Reuse `LadrunoJ2Kernel.h`, `LogStrainKernel.h`, `LadrunoHardening.h` (the same
  shared σ_y), the polar primitive.
- `FEM_ObjectBroker` registration, the `nDMaterial` command, `CMakeLists.txt`,
  banner line, stamp-headers glob.
- Tests: the numpy oracle (`ladrunoj2_finite_native_reference.py` — written, now
  also carries the consistent-tangent helpers) is the target; the tangent recipe is
  pinned by `tests/test_ladrunoJ2_finite_native_tangent.py` (3 tests: A+B recipe
  exact-to-FD, channel-B small-but-present, channel-B objective). Add a g++
  kernel-level check and a built `LadrunoBrick -geom finite` acceptance
  (necking-with-reversal), and **flip the v1 strict-xfail's positive twin**.

## Validation plan / definition of done

1. **The gate:** the C++ material passes `test_native_objective_under_superposed_rotation`
   (σ_rot = Q σ Qᵀ on a plastic combined-hardening state) — the property v1 fails.
2. Matches the numpy oracle step-for-step along the rotating+stretching plastic path
   (~1e-9), and reduces to the v1 / small-strain result with no rotation.
3. Element-level: `LadrunoBrick -geom finite` over the native material — consistent
   tangent vs FD (with channel B wired, this must hold at **tight** tolerance even in
   a saturated-backstress state, not just near-symmetric `F`), a finite-rotation
   **cyclic** load (buckling-brace-like) showing a correct Bauschinger loop, plus the
   existing finite battery. Material-level tangent (incl. the co-rotation term) is
   pre-pinned in numpy per the section above.
4. The `…_is_v2` strict xfail in `test_ladrunoJ2_finite.py` gets a passing native twin.

## Scope / non-goals

- **In:** isotropic + Chaboche AF kinematic (the full `LadrunoJ2` law), 3D, the
  co-rotated backstress, spatial tangent, serialization of `α` as spatial state.
- **Out (v2.x):** plane-stress finite (the §14.7 nested route), the IMPL-EX hook,
  tabulated iso (gated on the small-strain tabulated mode), thermomechanical coupling.

## Consistent tangent — the co-rotation term (RESOLVED with data, 2026-06-02)

The Cauchy stress depends on the trial `F` through **two** channels:
- **(A)** the trial elastic strain `εᵉᵗʳ = ½ ln(f_Δ bᵉ_n f_Δᵀ)` — the standard
  log-strain dependence, captured by the **already-shipped** spatial tangent
  (`LogStrainKernel::spatial_tangent_full` / numpy `spatial_tangent_a`).
- **(B)** the **co-rotated backstress** `α̃ = R(f_Δ) α_n R(f_Δ)ᵀ` — `α̃` depends on
  `F` through `R = polar(f_Δ)`. This channel is **new** to the native material and
  is what a naive "reuse the log-strain tangent" would omit.

**Is it needed? — measured in numpy** (`tests/test_ladrunoJ2_finite_native_tangent.py`,
helpers on the oracle: `NativeFiniteJ2.tangent_dsigma_dF`, `.channelB_dsigma_dF`):

| property | finding |
|---|---|
| channel-B size `‖∂σ/∂F_full − ∂σ/∂F_(A only)‖ / ‖∂σ/∂F‖` | **~1e-4** (stiff metals, `G≫σ_y`) up to **~2e-3** (soft elasticity / strong kinematic) |
| vs. the existing tangent-FD gate (rtol **2e-4**) | **straddles it** ⇒ exact tangent needs B |
| effect on Newton robustness | none (≪1% off ⇒ still ~quadratic) |
| frame-objectivity of the tangent | channel-B fraction is **rotation-invariant** (objective) |
| channel B is non-zero only when | the step is **plastic** (elastic predictor is α-independent) |

**Verdict: include channel B** (it crosses the tangent-test tolerance), but it does
**not** justify the heavy machinery of a fully-analytic `∂R/∂F` (Sylvester) plus a
kernel extension for `∂τ/∂α̃`.

**The recipe (validated to ~1e-14 vs full FD):** keep channel A **analytic** (the
shipped spatial tangent) and add channel B **numerically by perturbing only `R`** —
hold `εᵉᵗʳ` and `J` at the base point, perturb `F → R(f_Δ)` in the 6 (sym) / 9
directions, and finite-difference the extra Cauchy-stress response through the
**unchanged** return map (≈9 cheap scalar-Newton solves, only when plastic). This
adds ~15 lines, needs **no** kernel re-derivation, and is exact-to-FD. Proven clean
and additive: `∂σ/∂F_full = ∂σ/∂F_(A, R-frozen) + ∂σ/∂F_(B, R-only)` to machine
precision across stiff/soft × first-yield/saturated × no-rot/large-rot.

(Residual open item — **RESOLVED 2026-06-02 by step-refinement, no swap needed.** See
the P1 entry under "Follow-ups" below: a numpy step-refinement study along a severe
non-coaxial path — `F(s)=exp(s·(D+W))`, 114° rotation + ~19% isochoric stretch with
stretch-rate `D` and spin `W` deliberately not coaxial — shows (i) the shipped polar
scheme converges at order ≈1 to a well-defined AF magnitude, and (ii) the polar(f_Δ)
vs exact continuous-spin (exp-map) transport difference is *second order* in the step
(≈O(Δs²)), 3–4 orders of magnitude below the first-order discretization error. The two
transports share the same limit, so the exact exp-map transport cannot meaningfully
improve accuracy. Confirms the expectation — the AF evolution runs in the correctly-
rotated frame via the unchanged return map, exactly as the isotropic part does.
Tests: `tests/test_ladrunoJ2_finite_native_steprefine.py` (5) +
`tests/ladrunoj2_finite_native_steprefine_reference.py`.)

## Status — MERGED (PR #127, 2026-06-02) + adversarial review

Shipped as `SRC/material/nD/LadrunoJ2Finite.{h,cpp}` (classTag 33012) via
**[PR #127](https://github.com/nmorabowen/OpenSees/pull/127)** (rebased onto
`ladruno` past the Lemaitre-damage merge; re-verified green post-rebase). Built
(OpenSeesPy, exit 0); full J2 battery green (34 passed / 1 xfailed — the v1
LogStrain-wrapper non-objectivity test correctly *stays* xfail, with its passing
native twin `test_native_objective_under_superposed_rotation`).

**34-agent adversarial review (7 dimensions, each finding refuted-or-confirmed by an
independent skeptic): NO confirmed correctness bug.** Re-derived correct and cleared:
the polar decomposition `R = f_Δ U⁻¹`, the `α̃ = R α Rᵀ` push-forward (tensor-Voigt,
loss-free symmetrisation), the elastic-trial feed (`epsP_n=0`) + `εᵉ`/`bᵉ` recovery,
the engineering↔tensor shear bridging, the channel-B `Ḟ=LF` map
`cmatB_ijkl = Σ_m (∂σ_ij/∂F_km) F_lm`, and frame-indifference. The one real (medium)
item was a **test-coverage gap**: the C++ channel-B tangent path was executed but not
numerically gated. Fixes applied this PR:
- **Closed the gap** — `tests/ladrunoj2_finite_native_tangentB_check.cpp` +
  `tests/test_ladrunoJ2_finite_native_tangentB_cpp.py`: the C++ channel-B tensor is
  compared (rel-Frobenius <1e-3) to the numpy oracle's `channelB_dsigma_dF` mapped
  through the identical `Σ_m·F_lm` transform — channel B is purely constitutive (no
  geometric term), so this is a clean, commit-semantics-free gate that catches any
  F→L index/factor/sign error.
- **Serialization** now tested — `test_native_finite_database_roundtrip` (FE_Datastore
  round-trip of a committed plastic finite state; skips without DB) + the false
  "database restore" docstring removed.
- **Hardening:** per-instance `K0init` (was a process-shared function-static `Matrix`
  in `getInitialTangent`); `polarRotation` eigenvalue clamp (no rank-deficient/NaN `R`
  on a near-singular channel-B perturbation); `override` on all contract methods
  (signature-drift → compile error, not a silent base-default fallback); backstress
  slot-5 `(0,2)==(2,0)` convention comment.

Element-level FD-of-K is intentionally **not** gated: channel B scales with the
committed backstress, but `eleResponse("stiff")` re-forms post-commit (α≠0) while a
from-virgin force-FD carries α=0 — different committed bases, and an FE_Datastore
replay fights OpenSees commit/sp semantics. Channel B is therefore gated off-element
(material-level exact-to-FD + the g++ tensor check above); the element side is covered
by the v1 assembly gate + Newton convergence on the rotated-plastic objectivity path.

## Follow-ups (next session)

v1 (PR #127) is complete for its scope (3D, full Voce+Chaboche, co-rotated backstress,
channel-A+B consistent tangent, serialization). Remaining work, prioritized — a future
session can start straight from here:

**P1 — physics caveat to close (mild correctness flavor, cheap): ✅ RESOLVED 2026-06-02 — no swap.**
- **AF magnitude under *simultaneous* large rotation AND large stretch.** Objectivity was
  proven; this was the one un-pinned *accuracy* question (`R_Δ=polar(f_Δ)` transport vs the
  exact §14.11 exponential-map transport differ at higher order). **Verified by step-refinement**
  (numpy oracle, no build needed): driver `tests/ladrunoj2_finite_native_steprefine_reference.py`
  drives the SAME return map + Hencky kinematics the C++ uses along `F(s)=exp(s·(D+W))` with a
  **non-coaxial** stretch-rate `D` and spin `W` (114° rotation + ~19% isochoric stretch) — the
  regime where the transport choice actually matters (a fixed-axis stretch is degenerate: all
  `U(s)` commute ⇒ every transport coincides). Test `tests/test_ladrunoJ2_finite_native_steprefine.py`
  (5 cases, all green) establishes:
  - the `polar` driver is **bit-identical** to the shipped `NativeFiniteJ2` along the path
    (it is the same logic, not a re-model);
  - the shipped polar scheme **converges at observed order ≈1** to a well-defined AF magnitude
    (no O(1) inconsistency);
  - polar(f_Δ) and the **exact** continuous-spin transport (composed micro-polars → exp-map of
    the vorticity, M→∞) **share the same limit**, their gap shrinking at **order ≈1.9** (second
    order in the step);
  - at every coarse resolution the polar-vs-exact *transport-choice* difference is **3–4 orders
    of magnitude below** the step-discretization error (e.g. N=12: ~3e-4 of it) — i.e. dominated
    by the unavoidable first-order return-map discretization.

  **Verdict: keep the incremental polar transport. The exact §14.11 exp-map transport would not
  measurably improve accuracy.** The caveat is closed; no code change.

**P2 — deferred features (each its own PR; the "Out (v2.x)" non-goals):**
- **IMPL-EX** `-implex` — ✅ **SHIPPED (PR [#134](https://github.com/nmorabowen/OpenSees/pull/134), 2026-06-02).** See the dedicated
  section below.
- **Plane-stress / dimensional finite views** (§14.7 nested route). `LadrunoJ2Finite` is
  **3D-only** (`getType`=="ThreeDimensional"); `LadrunoJ2` already has 5 views. Finite
  plane-stress needs the nested out-of-plane iteration. **DEFERRED — no consumer:** the only
  element driving a `FiniteStrainNDMaterial` (via `setTrialF`) is the 3D `LadrunoBrick`; there
  is no 2D finite element in the fork, and standard 2D continuum elements drive materials via
  `setTrialStrain` (disabled here). So dimensional finite views are orphaned until a 2D
  finite element exists (a plane-strain study can already be done with a 3D brick + fixed
  out-of-plane DOFs). Confirmed 2026-06-02.
- **Tabulated / Bézier isotropic curve** — gated on the small-strain `LadrunoJ2` tabulated
  mode landing first (shared `LadrunoHardening.h`).
- **Thermomechanical coupling.**

### IMPL-EX (`-implex`) — Δγ-extrapolation, SHIPPED PR [#134](https://github.com/nmorabowen/OpenSees/pull/134) (2026-06-02)

The classic Oliver–Huespe–Cante implicit/explicit split, on the **plastic multiplier**
(not the damage variable — that is the small-strain `LadrunoJ2 -implex` Lemaitre path).

**Scheme.** Every step the implicit co-rotated return map runs and is **committed** — so the
committed finite history (`bᵉ_n, ε̄ᵖ_n, α_{n,k}`) is byte-identical to a fully-implicit run.
Only the stress/tangent *reported to the solver* are replaced, freezing two history quantities:
the extrapolated multiplier `Δγ̃ = Δγ_n·(Δt_{n+1}/Δt_n)` (uniform step ⇒ `Δγ_n`) and the
committed log-strain flow direction `N_n`, **co-rotated by the SAME `R_Δ = polar(f_Δ)`** as the
backstress:

    εᵉ̃ = εᵉᵗʳ − Δγ̃ · Ñ_n,   Ñ_n = R_Δ N_n R_Δᵀ;   τ̃ = ℂᵉ : εᵉ̃;   σ̃ = τ̃/J.

Two consequences make it clean: (i) `Ñ_n` co-rotates ⇒ `σ̃` is **frame-indifferent** (objective
to ~7e-13); (ii) the backstress **drops out** of `σ̃` (it appears only through the implicit
`N_n`, already frozen) ⇒ the reported material tangent is the **constant SPD elastic operator**
— no plastic `h`, and the co-rotation "channel B" vanishes. The price is an O(Δt) consistency
error (the explicit stress → the implicit stress at **observed order ≈2** under step refinement
for smooth flow) and a one-step lag at first yield (`Δγ_n = 0` ⇒ `σ̃` = elastic trial).

**Why for the *finite* J2:** plain hardening J2 has an SPD implicit tangent already, so IMPL-EX
is not a robustness cure on its own; its payoff is the factor-once constant SPD tangent for
**explicit / quasi-static** use and as the SPD-tangent enabler when finite J2 is later paired
with a **softening** law (finite Lemaitre). It mirrors the small-strain damage IMPL-EX (frozen
extrapolated `D~ ⇒ (1−D~)Dtan_eff`); here the frozen quantity is `Δγ`.

**Implementation:** entirely in `LadrunoJ2Finite.{h,cpp}` — **no kernel change**. The flow
direction is recovered from the implicit return as `N = Δεᵖ/Δγ`; committed `Δγ_n, N_n` added to
the state + sendSelf/recvSelf; default off ⇒ bit-identical to the implicit material.

**Validation (oracle-first, no build needed for the constitutive proof):**
- `tests/ladrunoj2_implex_reference.py` + `..._finite_implex_reference.py` (numpy oracles) +
  `test_ladrunoJ2_implex_reference.py` (4) + `test_ladrunoJ2_finite_implex_reference.py` (5):
  committed history == implicit (byte), explicit→implicit order≈2, objective 7e-13, constant
  SPD elastic tangent, finite→small-strain reduction, documented yield-onset lag.
- `tests/ladrunoj2_finite_implex_check.cpp` + `test_ladrunoJ2_finite_implex_cpp.py`: standalone
  g++ mirror of the C++ `useImplex` math (Voigt conventions, flow-direction recovery,
  co-rotation) vs the oracle to 1e-9 — the bridge the 3×3 oracle cannot catch.
- Element acceptance: `test_ladrunoJ2Finite_element.py` (`LadrunoBrick -geom finite` with
  `-implex`) — built-binary gate.

**P3 — validation / perf (not code gaps; nice-to-haves):**
- **Finite cyclic Bauschinger / buckling-brace element demonstrator** — v1 ships objectivity
  + tangent but no explicit cyclic-loop element test (the validation-plan item §3 floated it).
- **`-formulation bbar/uri/eas` coverage with the native material** — only `std` is tested;
  the element drives `setTrialF` regardless of formulation, so these should work but aren't
  gated (note: `bbar`+finite = F-bar, already shipped on the element side).
- **Analytic channel B** — ✅ **SHIPPED (PR [#139](https://github.com/nmorabowen/OpenSees/pull/139), 2026-06-02).** The numeric
  R-perturbation (~18 return-map calls per GP-tangent) is replaced by the analytic chain
  `∂R/∂F` (polar-rotation derivative: `Ω` skew solving `(tr U·I − U) ω = axial(A−Aᵀ)`,
  `A = Rᵀ df`) ∘ `∂α̃/∂R` ∘ `∂τ/∂α̃` (return-map backstress sensitivity:
  `∂Δγ/∂α̃_m = −n/(h D_m)`, `dM = −dα̃_m/D_m + Mp·dΔγ`, `dn = (dM − n(n:dM))/‖M‖`,
  `ds = −2G(dΔγ·n + Δγ·dn)`) — no return-map calls, no FD step. Default; the numeric path
  is retained behind `-DLADRUNO_J2FINITE_CHANNELB_NUMERIC` for validation. Oracle
  `tests/ladrunoj2_finite_channelB_reference.py` + `test_ladrunoJ2_finite_channelB.py` (5:
  polar-deriv vs FD ~1e-9, analytic vs numeric ~3e-8 stiff & soft, A+B = full-tangent,
  elastic⇒0) + g++ `test_ladrunoJ2_finite_channelB_cpp.py` (C++ analytic vs oracle ~1e-7).

**Also shipped in PR [#139](https://github.com/nmorabowen/OpenSees/pull/139) (material robustness/docs):**
- **Softening-parameter guard** (the LEDGER_quirks `h>0` caveat): the OPS factory now warns
  when the minimum isotropic hardening slope `σ_y'_min = Hiso + min(0,Qinf)·b < 0` (the
  consistent tangent may become indefinite; below `−3G` it loses positive-definiteness),
  recommending `-implex` (constant SPD tangent) or crack-band regularization. Warn, not
  reject — `-implex` makes softening well-posed.
- **IMPL-EX uniform-Δt assumption documented** — `Δγ̃ = Δγ_n` bakes in `Δt_{n+1}/Δt_n = 1`;
  flagged in the factory usage string + `Print` (accuracy degrades under variable Δt; wiring
  the dt-ratio is future work, blocked on materials not cleanly receiving Δt).

NB the kernel-sharing is **nD↔finite-strain only**, not 1D — there is no `LadrunoJ2Finite`
analogue for the uniaxial twin (a 1D backstress is a scalar with no frame to co-rotate).
