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

**Status:** implemented — v1 = monotonic Dhakal–Maekawa (2026-06-02) · **Wraps:**
any tension-compression `UniaxialMaterial` (designed for
[[13_ladruno_uniaxial_j2_adr]] `LadrunoUniaxialJ2`, but composes over
`Steel02`/`Steel4` too) · **Registry:** `UniaxialMaterial` · **Oracle:** the
`ReinforcingSteel::Buckled_stress_Dhakal` formula (ported; see the B2 caveat in
§6 — a *tight-tolerance*, not literally bit-for-bit, match because
`ReinforcingSteel` works in **log** strain on its own monolithic backbone while
this wrapper works in **engineering** strain on a generic bar).

> [!note] As-built scope (v1 + GA)
> Implements the **monotonic** degradation for **both** laws — the clean
> `σ_buckled = r(e,λ)·σ_bare` form (DM's `TBranchNum%4 > 1` path; GA's
> stress-only `Buckled_stress_Gomes`). **`-model dm`** (default) and **`-model
> ga`** are both shipped. The DM **cyclic re-straightening** branch (upstream's
> `TBranchNum%4 ≤ 1`, backstress-anchored Menegotto-Pinto blend off
> `Tfa`/`BackStress`) is **NOT** a stress-only modification of `σ_bare` and is
> **deferred to v2 (test B4)**; the as-built state is the minimal `ε_max`/`σ_max`
> anchor + cached onset backbone `fStarL` (no `Cbranch` flag — see §4). Where
> §3/§4 below describe richer cyclic behavior, read it as the v2 target.
>
> **GA specifics:** `-reduction r` (RS GABuck `r`, clamped [0,1]) and `-fsufrac g`
> (RS GABuck `gama` = `fsu_fraction`) are the GA knobs; GA does **not** use `-fy`.
> The wrapper defaults `-reduction 0.0` (full GA buckling when `lsr>0`) — this
> **deviates from RS's `r=1` default** (which means *no* buckling), chosen so
> `-model ga -lsr X` actually buckles, matching the wrapper's "lsr>0 ⇒ buckle"
> contract. NB the upstream `Buckled_stress_Gomes` **shadows** its own
> `beta`/`gama` locals (hardcoded 1.0 / 0.1), so the GABuck `beta` argument is
> dead — faithfully reproduced. GA tangent: hybrid analytic — `ρ·k_bare +
> (∂σ/∂factor)·(∂factor/∂ε)` with `∂factor/∂ε` from a cheap central FD of the
> closed-form scalar `factor` (no material probe), strictly better than RS's
> `Buckled_mod_Gomes` (which FDs the whole stress with `σ_bare` frozen, then adds
> `k_bare` instead of `ρ·k_bare`).

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
   **Both shipped:** `dm` (default) and `ga`. GA adds `-reduction`/`-fsufrac`
   knobs and does not use `-fy` (see the as-built note above).
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
f*L = backbone(ε*)                          # bare stress at the onset strain ε*
                                            # (ε* treated as an ABSOLUTE backbone
                                            #  strain, mirroring upstream
                                            #  Backbone_f(eStar); == the onset's
                                            #  true absolute strain when e_cross=0)
f*  = f*L · β · (1.1 − 0.016·√(f_y/E·2000)·λ),   floored: if f* > −0.2 f_y ⇒ f* = −0.2 f_y
                                            # intermediate branch  −ε_y > e ≥ ε* :
σ_buckled = σ_bare · [ 1 − (1 − f*/f*L)·(e+ε_y)/(ε*+ε_y) ]   = r(e,λ)·σ_bare
                                            # deep branch  e < ε* :
σ_buckled = σ_bare · ( f* − 0.02·E·(e − ε*) ) / f*L,   then floored: if σ_buckled > −0.2 f_y
                                            #                              ⇒ σ_buckled = −0.2 f_y
                                            # (C0-continuous at e=ε*: both branches give σ_bare·f*/f*L there)
```

`β` here is the single DM residual-shape parameter (the `-alpha` flag /
`ReinforcingSteel`'s `beta`, default 1.0) — there is no second independent
shape constant.

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

**Cyclic re-straightening (v2 — NOT in v1).** On tension reload a bowed bar
straightens and the response transitions back toward the bare curve. In
`ReinforcingSteel` this is the `TBranchNum%4 ≤ 1` path, which is **not** a
stress-only modification of `σ_bare`: it computes the buckled stress off the
branch-start anchor stress `Tfa` and a Menegotto-Pinto blend through
`BackStress = MP_f(e_cross−ε_y)`. Because it depends on more than `(σ_bare, e,
λ)`, the clean wrapper form does not express it; **v1 omits it** and applies the
monotonic `r·σ_bare` form whenever `e < −ε_y`. v1 *does* track the `ε_max`/
`σ_max` anchor (so `e_cross` and the onset shift with the last tensile
excursion), but the re-straightening branch shape itself is **deferred to v2
(test B4)**. NB the anchor `e_cross = ε_max − σ_max/E` uses the stress recorded
at the max-strain sample as a proxy for upstream's reversal-point `f_sup`; these
coincide for clean monotonic-tension-then-reverse histories and diverge under
sub-step overshoot — a documented v1 approximation.

**Gomes–Appleton (`-model ga`) — SHIPPED.** Alternative rebar law (mirrors
`Buckled_stress_Gomes`), same wrapper plumbing. Gate: pass-through unless
`ε < e_cross` (buckles for any compression past the anchor, unlike DM's
`e < −ε_y`). With `fs_buck = √(32/(e_cross−ε)) / (3π·λ)`, `m = min(1, fs_buck)`,
a local `β(stress_diff)` knee at `Dft=0.25`, and `factor = m·β·(1−r) + r`:

```
σ_buckled = f_sup·g − (factor + g)·(f_sup·g − σ_bare)/(1 + g)
```

where `f_sup = σ_max` (anchor stress), `g = fsu_fraction` (`-fsufrac`), `r =
reduction` (`-reduction`). For monotonic-from-virgin (`f_sup=0`) this collapses
to a clean `σ_bare·(factor+g)/(1+g)`. `r=1 ⇒ factor=1 ⇒` no buckling; `r=0 ⇒`
full GA. **Faithful quirk:** upstream shadows its `beta`/`gama` with hardcoded
locals (1.0 / 0.1), so the GABuck `beta` argument is dead — reproduced.

---

## 4. Architecture — the wrapper

```cpp
// As-built v1 (monotonic DM). The cyclic `Cbranch` flag / update_buckling_state
// machinery sketched in earlier drafts is a v2 concern; v1 needs only the
// ε_max/σ_max anchor + a cached onset backbone.
class LadrunoRebarBuckling : public UniaxialMaterial {
  UniaxialMaterial* theBar;     // wrapped bare-bar (owns a getCopy)
  int    model;                 // MODEL_DM (MODEL_GA reserved, parse-guarded off)
  double lsr;                   // slenderness s/d   (lsr<=0 => exact pass-through)
  double alpha;                 // the single DM residual-shape factor (== RS beta)
  double fy, E;                 // for eStar/eY and the deep-branch slope (E from -E or the bar)
  // committed cyclic state (the ONLY state the wrapper owns):
  double CmaxTenStrain, CmaxTenStress;   // ε_max / σ_max anchor -> e_cross
  double CfStarL, CfStarLcross;          // cached onset backbone f*L + the e_cross it was computed for
  // (+ trial twins T..., trial output Tstrain/Tstress/Ttangent/Tr)
};

int setTrialStrain(double eps, double rate) {
  theBar->setTrialStrain(eps, rate);            // core evolves on the TRUE strain
  double sBare = theBar->getStress();
  double kBare = theBar->getTangent();
  // update the ε_max/σ_max anchor (invalidate the f*L cache if it moves)
  if (eps > TmaxTenStrain) { TmaxTenStrain = eps; TmaxTenStress = sBare; invalidate_cache(); }
  if (lsr <= 0.0 || fy <= 0.0 || E <= 0.0) {
     Tstress = sBare;  Ttangent = kBare;        // identity gate (oracle preserved)
     return 0;
  }
  applyBuckling(eps, sBare, kBare);             // e=eps−e_cross; if e<−ε_y: Tstress=r·sBare,
                                                // Ttangent=r·kBare+sBare·r' (analytic); else pass-through
  return 0;
}
```

- **State hygiene:** the wrapper's `commitState`/`revertToLastCommit`/`revertToStart`
  manage *only* `CmaxTenStrain`/`CmaxTenStress`/`CfStarL`/`CfStarLcross`; the
  bar's plastic state is committed by forwarding to `theBar->commitState()`. The
  bar is the real material state; the wrapper is the geometric overlay. **`getCopy`
  must NOT route through a ctor that calls `revertToStart()`** — that would wipe
  the state `bar.getCopy()` just preserved (the classic wrapper footgun;
  `FatigueMaterial`/`MinMaxMaterial` avoid it the same way) — so the value ctor
  initializes the wrapper doubles inline and `getCopy` restores `C*`/`T*` after.
- **`sendSelf`/`recvSelf`:** serialize wrapper params + committed anchor/cache
  `(CmaxTenStrain, CmaxTenStress, CfStarL, CfStarLcross)` **and** the wrapped
  material (tag + class via the broker, the standard nested-material idiom — see
  `FatigueMaterial`/`MinMaxMaterial`). Trial state is recomputed on the next
  `setTrialStrain` (the `fStarL` cache rebuilds via the probe), so it is not sent.
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

- `-lsr` (slenderness `s/d`) is the primary buckling input; `-model dm`
  (default; `ga` reserved). When `-lsr > 0`, **`-fy` is required** — a generic
  `UniaxialMaterial` exposes no yield-stress getter, so `f_y` cannot be inferred
  (only `E` can, from `theBar->getInitialTangent()` when `-E` is omitted). `-E`
  **should equal the wrapped bar's elastic modulus**: it sets `ε_y = f_y/E` and
  the deep-branch slope `0.02·E`, so a mismatched `-E` desyncs the onset/floor
  geometry from the bare curve.
- `-lsr 0` (or `≤ 0`) ⇒ exact pass-through (the no-buckling identity gate).

---

## 6. Testing / oracle matrix (Zone-A)

| # | Test | Oracle |
|---|---|---|
| B0 | `-lsr 0` pass-through ≡ bare material | identical σ/tangent to the wrapped `LadrunoUniaxialJ2` (the identity gate) |
| B1 | tension + pre-onset compression ≡ bare | no modification before `ε*` |
| **B2** | monotonic compression past onset, several `λ` | the ported `Buckled_stress_Dhakal` formula, fed the bare bar's own `σ_bare(ε)` + its onset backbone `f*L = σ_bare(ε*)` — a **tight-tolerance** match (rel ~1e-6), NOT literally bit-for-bit vs `ReinforcingSteel -DMBuck`: RS works in **log** strain on its **own** monolithic backbone, this wrapper in **engineering** strain on a **generic** bar, so the two materials cannot coincide to machine precision. Oracle = the formula port, verified self-consistently against the wrapped bar. |
| B3 | larger `λ` ⇒ earlier onset + lower residual | monotone trend in `ε*`, residual → `−0.2 f_y` floor |
| B4 | cyclic re-straightening | tension reload recovers toward bare; loop shape vs `ReinforcingSteel` |
| B5 | consistent-tangent FD | analytic `r·k+σ·r'` vs one-step central FD (same fresh-material idiom as J2 V6) |
| B6 | composition `Fatigue ∘ Buckling ∘ J2` builds + runs | rupture still triggers on the buckled response |
| B7 | `sendSelf`/`recvSelf` round-trip (nested material) | state + wrapped bar reconstructed |
| GA0 | `-model ga -lsr 0` pass-through ≡ bare (full push-pull) | identity gate for GA |
| GA1 | tension pass-through (`ε ≥ e_cross`) | GA gate before the anchor crossing |
| **GA2** | monotonic compression, several `λ` and `reduction` | the ported `Buckled_stress_Gomes` formula fed `σ_bare(ε)` + `f_sup` anchor (rel ~1e-6) |
| GA3 | `reduction`=1 ⇒ no buckling; `reduction`=0 ⇒ knock-down; larger `λ` ⇒ more knock-down | the GA blend + slenderness trend |
| GA5 | consistent-tangent FD (smooth region, off the `fs_buck≈1` knee) | hybrid analytic `ρ·k + (∂σ/∂factor)·factor'` vs central FD |

Smoke (L0): single truss / 1-fiber RC section push-pull-cycle, Zone-A pytest.
**Implemented:** B0/B1/B2/B2b/B3/B5 (DM) + GA0/GA1/GA2/GA3/GA5 (GA) + **B6
composition** (`Fatigue ∘ RebarBuckling ∘ Steel02`: buckling propagates through
the outer Fatigue layer, and Fatigue's `-min` rupture still triggers on the
buckled response) + **B7 `sendSelf`/`recvSelf` round-trip** (via `FE_Datastore`
`database_roundtrip`, both DM and GA incl. the GA `gaReduction`/`gaFsuFrac`
fields) — all green; `LadrunoUniaxialJ2` regression unchanged. (Only **B4
cyclic** remains for a follow-up.)

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
3. **Cyclic re-straightening** — `ε_max` anchor + branch flag; gate B4. *(deferred)*
4. ✅ **GA law** (`-model ga`) — ported `Buckled_stress_Gomes`; reused the
   plumbing (PR #119 = DM; GA shipped in the follow-up PR). `-reduction`/`-fsufrac`
   knobs; gates GA0/GA1/GA2/GA3/GA5.
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
