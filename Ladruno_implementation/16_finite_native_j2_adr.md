---
title: ADR — Finite-strain-NATIVE combined-hardening J2 with a co-rotating backstress (§14.11 "v2")
project: Ladruno
status: proposed
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

**Proposed.** The numerical algorithm is **prototyped and verified in numpy**
(`tests/ladrunoj2_finite_native_reference.py`, 5/5 green) — no C++ yet.

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

- `SRC/material/nD/LadrunoFiniteJ2.{h,cpp}` (or `LogStrainCombinedJ2`) — the
  `FiniteStrainNDMaterial` subclass; classTag in the 3301x band (next free after
  `LadrunoJ2=33011` — reserve in `classTags.h`).
- Reuse `LadrunoJ2Kernel.h`, `LogStrainKernel.h`, `LadrunoHardening.h` (the same
  shared σ_y), the polar primitive.
- `FEM_ObjectBroker` registration, the `nDMaterial` command, `CMakeLists.txt`,
  banner line, stamp-headers glob.
- Tests: the numpy oracle (`ladrunoj2_finite_native_reference.py` — already written)
  is the target; add a g++ kernel-level check and a built `LadrunoBrick -geom finite`
  acceptance (necking-with-reversal), and **flip the v1 strict-xfail's positive twin**.

## Validation plan / definition of done

1. **The gate:** the C++ material passes `test_native_objective_under_superposed_rotation`
   (σ_rot = Q σ Qᵀ on a plastic combined-hardening state) — the property v1 fails.
2. Matches the numpy oracle step-for-step along the rotating+stretching plastic path
   (~1e-9), and reduces to the v1 / small-strain result with no rotation.
3. Element-level: `LadrunoBrick -geom finite` over the native material — consistent
   tangent vs FD, a finite-rotation **cyclic** load (buckling-brace-like) showing a
   correct Bauschinger loop, plus the existing finite battery.
4. The `…_is_v2` strict xfail in `test_ladrunoJ2_finite.py` gets a passing native twin.

## Scope / non-goals

- **In:** isotropic + Chaboche AF kinematic (the full `LadrunoJ2` law), 3D, the
  co-rotated backstress, spatial tangent, serialization of `α` as spatial state.
- **Out (v2.x):** plane-stress finite (the §14.7 nested route), the IMPL-EX hook,
  tabulated iso (gated on the small-strain tabulated mode), thermomechanical coupling.

## Open question (low risk)

The polar-incremental push `R_Δ = polar(f_Δ)` is verified **objective + small-strain
exact** to ~1e-12. dSNPO §14.11 presents the transport tied to the exponential-map
integrator; the two agree for objectivity and the small-strain limit. The one thing
the objectivity tests do **not** pin is the *magnitude* convergence of the AF
evolution under simultaneous large rotation **and** large stretch — verify by
step-refinement against a fine-increment reference when building, and if a
discrepancy appears adopt the exact §14.11 exp-map transport. (Expected to be a
non-issue: the AF evolution itself runs in the correctly-rotated frame via the
unchanged return map, exactly as the isotropic part does.)
