---
title: ADR — Rebar-buckling wrapper (LadrunoRebarBuckling)
project: Ladruno
status: proposed
priority: medium
owner: nmora
tags:
  - implementation
  - material
  - uniaxial
  - buckling
  - steel
  - adr
---

# ADR — Rebar-buckling wrapper (`LadrunoRebarBuckling`)

**Status:** proposed (2026-06-02) · **Wraps:** any tension-compression
`UniaxialMaterial` (designed for [[13_ladruno_uniaxial_j2_adr]] `LadrunoUniaxialJ2`,
but composes over `Steel02`/`Steel4` too) · **Registry:** `UniaxialMaterial` ·
**Oracle:** `ReinforcingSteel -DMBuck` (verbatim DM formula port).

> **Supersedes** the "core flag, NOT a wrapper" stance in
> [[13_ladruno_uniaxial_j2_adr]] §6.2 **for the Dhakal–Maekawa case.** That note
> was right for a buckling model that alters the *plastic flow*; DM does not —
> it is a post-hoc **modification of the bare-bar stress** (`σ_buckled =
> f(σ_bare, ε, λ)`), which a stress-modifying wrapper expresses exactly while
> leaving the core a pure oracle. Verified against the OpenSees implementation:
> `ReinforcingSteel::Buckled_stress_Dhakal(ess, fss)` is literally a function of
> (strain, bare stress).

---

## 1. Context

### 1.1 The phenomenon, and what is / isn't in scope

A reinforcing bar in compression between two ties behaves as a **strut buckling
over the tie spacing**: a 1-D column-on-restraints problem with a single
slenderness `λ = s/d` (unsupported length / bar diameter). Past a
slenderness-dependent strain the bar bows, and the **section-average**
compressive stress drops well below the bare-bar stress. This is genuinely a
*fiber-material* phenomenon (unlike plate or global buckling) and is the clean,
canonical case for a constitutive buckling model.

**In scope:** average compressive stress–strain degradation of a single bar fiber
due to bar buckling between restraints, with cyclic re-straightening. Rebar /
RC fiber sections.

**Explicitly NOT in scope** (different scale, different tool — recorded so the
boundary stays sharp):
- **Steel-profile local (plate) buckling** — 2-D plate buckling with post-buckling
  reserve + boundary-condition (`k`) dependence; no clean fiber model. Use
  element-level IMK (`Bilin`/`ModIMKPeakOriented02`, already in-tree) or the
  shell/continuum-corotational route (`LadrunoJ2` + `-geom corot` + Lemaitre).
- **Global member buckling** (`KL/r`) — geometric; the corotational element's job
  (`LadrunoBrick -geom corot` / beam corotational transforms + an imperfection),
  never a material. **Do not double-count** with a buckling fiber (see §8).

### 1.2 Why a wrapper (the clean-separation argument)

DM computes the buckled *average* stress as a modification of the bare-bar curve
at the same strain — the bar's material still follows its own σ–ε; the
cross-section average drops because of the bowing. So the natural decomposition
is **core = the bar material (a pure oracle), wrapper = the geometric buckling
average on top.** Benefits:
- `LadrunoUniaxialJ2` stays byte-untouched — the V7 oracle property is preserved.
- **Stackable, opt-in layers:** `Fatigue ∘ Buckling ∘ J2` (or Lemaitre instead of
  Fatigue) — each independently testable.
- **Reusable** over any bare-bar `UniaxialMaterial` (`Steel02`, `Steel4`, ours).
- Strictly more modular than `ReinforcingSteel`, which bakes buckling *into* one
  monolithic material.

---

## 2. Decision

1. **New `UniaxialMaterial` wrapper `LadrunoRebarBuckling`** that holds a pointer
   to a wrapped bare-bar material, forwards `setTrialStrain`, reads its stress +
   tangent, and applies the **Dhakal–Maekawa** buckled-average modification in
   the compression branch; tension and pre-onset pass through unchanged.
2. **Law selectable** `-model dm|ga` — Dhakal–Maekawa (default) and
   Gomes–Appleton (the two rebar buckling laws OpenSees `ReinforcingSteel`
   offers). Both are `σ_buckled = f(σ_bare, ε, λ)` ⇒ both wrapper-shaped.
3. **classTag `MAT_TAG_LadrunoRebarBuckling = 33001`** (next free in the Ladruno
   *uniaxial* band after `LadrunoUniaxialJ2` = 33000).
4. **Analytic consistent tangent** `dσ/dε = r·k_bare + σ_bare·(∂r/∂ε)` (product
   rule) — improving on `ReinforcingSteel`, which finite-differences it
   (`Buckled_mod_Dhakal`). FD fallback available for the GA law if its `∂r/∂ε` is
   awkward.
5. **Port the DM formula verbatim** from `ReinforcingSteel::Buckled_stress_Dhakal`
   but **re-implement the cyclic bookkeeping cleanly** — the upstream formula is
   entangled with `ReinforcingSteel`'s branch-state machine (`TBranchNum`,
   `BackStress`, `Backbone_f`); the wrapper carries its own minimal state instead
   (§4).

---

## 3. The model — Dhakal–Maekawa (structural/constitutive section)

Reference: Dhakal & Maekawa (2002), *J. Struct. Eng.* 128(9) — "Reinforcing bar
stress–strain relations including buckling" + the cyclic companion. Mirrors
`ReinforcingSteel::Buckled_stress_Dhakal` (verified in-tree).

**Slenderness gate.** `λ = s/d`. Buckling negligible for `λ ≲ 5`; the wrapper is
a pure pass-through when `λ ≤ λ_min` (so `λ=0`/unset ⇒ exactly the bare material).

**Onset strain** (slenderness-dependent, measured from the unloading point at max
tensile strain `ε_max` — buckling is referenced to the tension excursion):

```
ε*  = −max(7, 55 − 2.3·√(f_y/E·2000)·λ) · ε_y          (compressive, magnitude grows with λ)
```

Larger `λ` ⇒ `ε*` reached sooner (earlier, more severe buckling). `ε_y = f_y/E`.

**Buckled average stress** (compression, relative to the bare stress `σ_bare` at
the same strain; `e` = strain measured from the tensile-unload anchor):

```
f*L = backbone(ε*)                          # bare stress at onset
f*  = f*L · β · (1.1 − 0.016·√(f_y/E·2000)·λ),   floored at −0.2 f_y   # residual level
                                            # intermediate branch  −ε_y > e ≥ ε* :
σ_buckled = σ_bare · [ 1 − (1 − f*/f*L)·(e+ε_y)/(ε*+ε_y) ]   = r(e,λ)·σ_bare
                                            # deep branch  e < ε* : ramp to the −0.2 f_y floor
```

So in the operative branch `σ_buckled = r(e,λ)·σ_bare` — a **multiplicative
reduction `r ∈ [resid, 1]` on the live bare stress**, gated by `e` and `λ`. This
is the exact shape that makes the wrapper clean.

**Consistent tangent** (analytic, product rule):

```
dσ_buckled/dε = r·(dσ_bare/dε) + σ_bare·(∂r/∂ε) = r·k_bare + σ_bare·r'
```

`r' = ∂r/∂ε` is closed-form on each branch (linear in `e`), so no FD needed —
the upstream `Buckled_mod_Dhakal` central-differences this; we do it analytically
and the FD becomes a *test oracle* (V-tangent), not the shipping path.

**Cyclic re-straightening.** On tension reload a bowed bar straightens and the
response transitions back toward the bare curve. The wrapper tracks `ε_max` (the
max tensile strain anchor `e_cross = ε_max − f_sup/E`) and the
buckled/unbuckled branch, reproducing DM's path dependence — but with its **own**
small state, not `ReinforcingSteel`'s `TBranchNum` machine.

**Gomes–Appleton (`-model ga`).** Alternative rebar law,
`σ_buckled = f_buck(λ, e_cross−ε)·…` (mirrors `Buckled_stress_Gomes`); same
wrapper plumbing, different `r`.

---

## 4. Architecture — the wrapper

```cpp
class LadrunoRebarBuckling : public UniaxialMaterial {
  UniaxialMaterial* theBar;     // wrapped bare-bar (owns a getCopy)
  int    model;                 // DM | GA
  double lsr;                   // slenderness s/d   (lsr<=lsrMin => pass-through)
  double alpha, beta, ...;      // DM/GA shape params (defaults per ReinforcingSteel)
  double fy, E;                 // queried from the bar (or supplied) for eStar/eY
  // committed cyclic state (the ONLY state the wrapper owns):
  double CmaxTenStrain;         // ε_max anchor for e_cross / re-straightening
  int    Cbranch;               // buckled / pre-onset / straightening flag
  ...
};

int setTrialStrain(double eps, double rate) {
  theBar->setTrialStrain(eps, rate);            // core evolves on the TRUE strain
  double sBare = theBar->getStress();
  double kBare = theBar->getTangent();
  if (lsr <= LSR_MIN || in_tension_or_pre_onset(eps)) {
     Tstress = sBare;  Ttangent = kBare;        // pass-through (oracle preserved)
  } else {
     double r  = dm_reduction(eps, sBare, lsr, state);   // ∈ [resid,1]
     double rp = dm_reduction_slope(eps, sBare, lsr, state);
     Tstress  = r * sBare;
     Ttangent = r * kBare + sBare * rp;          // analytic consistent tangent
  }
  update_buckling_state(eps);                    // ε_max anchor, branch, re-straighten
  return 0;
}
```

- **State hygiene:** the wrapper's `commitState`/`revertToLastCommit`/`revertToStart`
  manage *only* `CmaxTenStrain`/`Cbranch`; the bar's plastic state is committed by
  forwarding to `theBar->commitState()`. The bar is the real material state; the
  wrapper is the geometric overlay.
- **`sendSelf`/`recvSelf`:** serialize wrapper params + `(CmaxTenStrain, Cbranch)`
  **and** the wrapped material (tag + class via the broker, the standard nested-
  material idiom — see `FatigueMaterial`/`MinMaxMaterial`).
- **`getInitialTangent`** = `theBar->getInitialTangent()` (buckling does not change
  the elastic stiffness).
- **Files / wiring** (mirrors `LadrunoUniaxialJ2`): new
  `SRC/material/uniaxial/LadrunoRebarBuckling.{h,cpp}` + `OPS_*` parser;
  `classTags.h` (33001), `FEM_ObjectBrokerAllClasses`, `OpenSeesUniaxialMaterial
  Commands` (Python registry), `TclModelBuilderUniaxialMaterialCommand` (Tcl),
  `material/uniaxial/CMakeLists.txt`. **Reference:** `FatigueMaterial`/
  `MinMaxMaterial` (the wrap-a-material idiom incl. nested send/recv);
  `ReinforcingSteel.cpp` `Buckled_stress_Dhakal`/`_Gomes` (the formulae to port).

---

## 5. Public API

```tcl
# bare bar (pure core) -> buckling layer -> optional fracture
uniaxialMaterial LadrunoUniaxialJ2    10  $E -iso voce ... -kin 3 ...
uniaxialMaterial LadrunoRebarBuckling 11  10  -lsr $s_over_d  -model dm  <-alpha α> <-fy $fy>
uniaxialMaterial Fatigue              12  11                 ;# Chaboche ∘ buckling ∘ LCF
# fiber section consumes tag 12
```

- `-lsr` (slenderness `s/d`) is the only required buckling input; `-model dm`
  (default) / `-model ga`; `-fy`/`-E` optional (queried from the bar if omitted).
- `-lsr 0` (or ≤ `LSR_MIN`) ⇒ exact pass-through (the no-buckling identity gate).

---

## 6. Testing / oracle matrix (Zone-A)

| # | Test | Oracle |
|---|---|---|
| B0 | `-lsr 0` pass-through ≡ bare material | identical σ/tangent to the wrapped `LadrunoUniaxialJ2` (the identity gate) |
| B1 | tension + pre-onset compression ≡ bare | no modification before `ε*` |
| **B2** | monotonic compression past onset, several `λ` | **bit-for-bit vs `ReinforcingSteel -DMBuck`** (same `Buckled_stress_Dhakal`) — the porting oracle |
| B3 | larger `λ` ⇒ earlier onset + lower residual | monotone trend in `ε*`, residual → `−0.2 f_y` floor |
| B4 | cyclic re-straightening | tension reload recovers toward bare; loop shape vs `ReinforcingSteel` |
| B5 | consistent-tangent FD | analytic `r·k+σ·r'` vs one-step central FD (same fresh-material idiom as J2 V6) |
| B6 | composition `Fatigue ∘ Buckling ∘ J2` builds + runs | rupture still triggers on the buckled response |
| B7 | `sendSelf`/`recvSelf` round-trip (nested material) | state + wrapped bar reconstructed |

Smoke (L0): single truss / 1-fiber RC section push-pull-cycle, Zone-A pytest.

---

## 7. Risks / open questions

> [!question]
> **Cyclic branch logic is the fiddly part.** DM's re-straightening transitions
> are where `ReinforcingSteel` spends its `TBranchNum` complexity. The wrapper
> must reproduce them with minimal state (`ε_max` anchor + a branch flag). Build
> monotonic (B2) first, then layer the cyclic rules (B4); keep B2 bit-for-bit as
> the anchor.

> [!question]
> **Don't double-count with global buckling.** If the element is corotational
> (global buckling geometric) AND the fiber buckles (local), calibrate so they
> don't overlap. Document: rebar buckling fiber is for *bar* buckling between
> ties, not member buckling.

- **Softening ⇒ localization.** The post-onset negative tangent localizes; for
  fiber-section elements this is usually acceptable (it's the section average),
  but flag the same regularization caveat as Lemaitre / ASDConcrete
  ([[project_bezier_charlen]]) and note IMPL-EX as the robustness hook if implicit
  Newton stalls.
- **Tangent we improve on the precedent** — analytic vs upstream FD; B5 guards it.
- **No `exit()`**, real `revert*`/`Print`/`setParameter`; nested-material
  serialization correctness (the classic wrapper footgun).
- **Backwards compat:** new material, additive; composes, doesn't replace.

---

## 8. Implementation plan (step order)

1. **`LadrunoRebarBuckling.{h,cpp}`** wrapper skeleton: hold + forward to the
   wrapped material, pass-through everything (`-lsr 0` identity). Gate B0/B1 green
   first — proves the wrapper plumbing before any buckling math.
2. **DM monotonic compression** — port `Buckled_stress_Dhakal` (`ε*`, `f*`,
   `r(e,λ)`), analytic `r'` for the tangent. Gate **B2 bit-for-bit vs
   `ReinforcingSteel -DMBuck`** + B3 + B5 (tangent FD).
3. **Cyclic re-straightening** — `ε_max` anchor + branch flag; gate B4.
4. **GA law** (`-model ga`) — port `Buckled_stress_Gomes`; reuse the plumbing.
5. **Parser + registration + build** — classTag 33001, broker, Python/Tcl
   registries, CMake (mirror `LadrunoUniaxialJ2`); nested-material `sendSelf`/
   `recvSelf` (B7) + composition test (B6).
6. **Build-control:** `LEDGER_implementations` row (classTag 33001, uniaxial,
   wrapper); `LEDGER_vanilla_files` rows for the registry/broker/classTags edits;
   `manifest.yaml` row (G9 gate — **do not forget**, it blocked #99); banner line
   → `patch_banner.py`.
7. **PR off `ladruno`** (`--base ladruno`); verify OPEN before any follow-up push;
   add the `manifest.yaml` row in the same PR so the classTag gate stays green.

**Deferred:** steel-profile plate-buckling law (route c — a *different* material,
not this one); interaction-calibration guidance with corotational global buckling.
