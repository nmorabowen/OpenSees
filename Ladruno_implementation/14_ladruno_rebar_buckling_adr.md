---
title: ADR — Rebar-buckling wrapper (LadrunoRebarBuckling)
project: Ladruno
status: v1 implemented; v2 cyclic design hardened
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
> §3/§4 below describe richer cyclic behavior, read it as the v2 target — now
> fully specified in **§9 (v2 cyclic re-straightening design)**.
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
3. **Cyclic re-straightening** — additive residual-gap memory (`σ = σ_bare − g`,
   `g ∈ [0, g_law]`) + a `{PASS, BUCKLING, RESTRAIGHTEN}` branch machine; gate B4.
   *(design hardened in §9 — multi-agent review + design panel; code deferred)*
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

---

## 9. v2 — Cyclic re-straightening — DESIGN (hardened; decisions locked, no code yet)

> [!note] v2 status (2026-06-02)
> This section was hardened by a multi-agent adversarial review (45 agents, 39
> findings → 31 verified) and a 3-candidate design panel for the residual-gap
> memory. **Four owner decisions are locked** and folded in below:
> 1. **Smoothstep** closing factor (not linear) — C1-clean at both Phase-2 seams.
> 2. **`k_rev = theBar->getInitialTangent()`** — the frozen reversal tangent.
> 3. **Both `L_rs` forms** ship, selected by `-restraighten` (§9.4).
> 4. **Residual-gap memory pulled forward** — the RESTRAIGHTEN→BUCKLING re-entry
>    is **C0-continuous**; the earlier "v2-A break / v2-B fixes it later" split is
>    **superseded** (there is no C0 break to defer). Read any remaining "v2-A/v2-B"
>    phrasing elsewhere in this doc as historical.

### 9.1 The problem, and why RS can't be ported here

When a bar that has buckled in compression is unloaded and reloaded in tension,
the bowed bar **straightens** before it can carry full tension: the
section-average stress climbs back toward the bare-bar curve over a finite
re-straightening strain span. v1 ignores this — past `−ε_y` it always applies the
monotonic `r·σ_bare` knock-down, so on tension reload it tracks the bare curve
immediately (no straightening lag) and on re-compression it re-buckles from the
current anchor. That is correct for monotonic compression but wrong for cyclic.

`ReinforcingSteel` captures the lag in its `Buckled_stress_Dhakal`
`TBranchNum%4 ≤ 1` branch:

```
BackStress = MP_f(e_cross − e_y);                       // RS Menegotto–Pinto curve eval
aveStress  = Tfa·(1 − (1−f*/f*L)·(e+e_y)/(e*+e_y));     // buckled value off the BRANCH-START stress Tfa
return BackStress − (BackStress − f_ss)·(BackStress − aveStress)/(BackStress − Tfa);
```

This is **not portable** to a material-agnostic wrapper: `Tfa` (branch-start
stress), `BackStress = MP_f(…)`, and `Cfa[]` are all RS-internal
Menegotto–Pinto hysteresis state. The wrapper only ever sees the wrapped bar's
`σ_bare`/`k_bare` and the buckled value *it itself* produced. So v2 **reformulates
re-straightening generically** rather than porting RS's branch.

### 9.2 Decision — generic re-straightening via an additive residual-gap memory

> **Locked:** v2 does **not** reproduce RS's MP-anchored blend. It carries the
> straightening lag in **one persistent committed scalar** — an additive stress
> **deficit `g ≥ 0` measured straight down from the bare curve** — using only
> quantities the wrapper already owns. Core stays a pure oracle; wrapper stays a
> geometry overlay. Physically motivated, **not** bit-equal to RS (RS works in
> log strain on its own MP backbone; §6 B2 already abandoned bit-for-bit).

**The invariant (the whole design in one line).** The wrapper output is *always*

```
σ = σ_bare(ε) − g ,        with   g ∈ [0, g_law(ε)] ,
g_law(ε) = σ_bare(ε) − σ_v1(ε)  ≥ 0          (the v1 monotonic buckling-law gap)
```

where `σ_v1` is the v1 DM/GA result (`applyBucklingDM/GA`'s `Tstress`). **v1 is
exactly the case `g == g_law`; full straightening is `g == 0`.** This single
bracket is what makes every hard property fall out: reduce-to-v1 is bit-identical
(no reload ⇒ `g = g_law`), energy cannot be injected (`g ≥ 0` ⇒ `σ ≤ σ_bare`),
and the gap is trapped in a compact interval that does not grow with cycle count
(convergence — see below).

**Per-call setup** (all pure functions of committed state + the trial strain;
computed once in `setTrialStrain` after `theBar->setTrialStrain(eps)`):

```
sBare = theBar->getStress();   kBare = theBar->getTangent();   eY = fy/E;
k_rev = theBar->getInitialTangent();              // DECISION 2: frozen, always-elastic unload
e_cross_live = TmaxTenStrain - TmaxTenStress/E;   // live v1 onset anchor (.cpp:260); only ADVANCES
```

**Stress law per branch** (`σ = sBare − g`; `tangent = kBare − dg/dε`):

```
PASS  (onset not reached):
    g = 0;   σ = sBare;   tangent = kBare.                     # exact v1 pass-through (.cpp:264-266)

BUCKLING:
    call applyBucklingDM/GA VERBATIM into scratch (σ_v1, k_v1);  g_law = sBare − σ_v1  (≥ 0)
    g = min( g_law,  Cg + max(0, g_law − Cg_law_reentry) )      # residual-gap reconciliation
    σ = sBare − g
    # virgin / never-straightened:  Cg = 0, Cg_law_reentry = 0  ⇒  g = g_law  (BIT-IDENTICAL v1)
    # carried offset (after an interrupted straighten): g re-grows from Cg by the v1 increment,
    #   clamped to g_law so it re-merges onto the v1 envelope by the deepest committed buckle.

RESTRAIGHTEN  (episode latches Ce_rev, Cs_rev, Ce_cross_rs, Cf_bare_cross, CL_rs all frozen):
    # Phase 1 — elastic unload (Ce_rev ≤ eps < Ce_cross_rs): no straightening yet.
    σ = Cs_rev + k_rev*(eps − Ce_rev);                         # tangent = k_rev
    # Phase 2 — tension-side straightening (eps ≥ Ce_cross_rs): close the gap.
    q  = clamp((eps − Ce_cross_rs)/CL_rs, 0, 1)
    B(q) = 1 − (3q² − 2q³)                                     # SMOOTHSTEP closing factor (DECISION 1)
    D0 = Cf_bare_cross − [Cs_rev + k_rev*(Ce_cross_rs − Ce_rev)]   # gap at the crossing (committed)
    g  = D0·B(q);   σ = sBare − D0·B(q)
```

**`B(q)` is smoothstep, `B(q) = 1 − S(q)` with `S(q) = 3q²−2q³`** (DECISION 1):
`B(0)=1`, `B(1)=0`, and `B'(0)=B'(1)=0`, so **both** Phase-2 seams (the crossing
`q=0` and the completion `q=1`) are **C1-clean** — superseding the v2-A linear
`(1−q)·shape(q)` form, which was C1-discontinuous at both seams by `D0/L_rs`. A
concave polynomial would clean only one seam at a time (it trades the two ends
off); smoothstep zeroes the slope at both ends at zero cost.

**`D0` uses the cached `Cf_bare_cross`, NOT a live probe** (review fix). The
committed-clone `backboneStressAt` (`cpp:238-250`) reverts to the *current*
committed reload step before probing, so probing `Ce_cross_rs` after the bar has
plastified during reload does **not** return the bare stress that was at the
crossing. Instead `Cf_bare_cross = σ_bare(Ce_cross_rs)` is **latched once at the
Phase-1→Phase-2 crossing commit** (`= theBar->getStress()` when `eps` first
reaches `Ce_cross_rs`) and frozen for the episode — making `D0` an exact pure
function of committed state.

**Consistent tangent (signs verified by differentiating the stress law):**

```
PASS:            dσ/dε = kBare                                   (> 0)
BUCKLING:        dσ/dε = k_v1   — because within a step Cg, Cg_law_reentry are frozen so
                 dg/dε = d(g_law)/dε = kBare − k_v1  ⇒  kBare − (kBare − k_v1) = k_v1.
                 (= the exact v1 analytic tangent: intermediate kBare·ratio+sBare·dratio .cpp:296;
                  deep (kBare·num+sBare·(−0.02E))/fStarL .cpp:309; floored 0 .cpp:305. A carried
                  offset does NOT change the tangent ⇒ re-buckle stiffness == virgin-buckle stiffness.)
RESTRAIGHTEN P1: dσ/dε = k_rev                                   (> 0, elastic unload)
RESTRAIGHTEN P2: dσ/dε = kBare − D0·B'(q)/CL_rs = kBare + D0·6q(1−q)/CL_rs.
                 D0 ≥ 0 and 6q(1−q) ≥ 0 on [0,1] ⇒ term ≥ 0 (bar stiffens pulling straight);
                 the term vanishes at q=0 and q=1 ⇒ dσ/dε = kBare at BOTH seams (C1-clean).
```

> **Tangent-formula correction (review must-fix):** the earlier v2-A line
> `dσ = kBare + D0·d[(1−q)·shape(q)]/dε` carried the wrong sign in the general
> case — differentiating `σ = σ_bare − (closing)·D0` gives `dσ = kBare −
> D0·d(closing)/dε`. With the smoothstep closing factor `B(q)` this is
> `kBare − D0·B'(q)/CL_rs = kBare + D0·6q(1−q)/CL_rs` (the sign now comes out
> right because `B'(q) = −6q(1−q) ≤ 0`).

**C0 at all five transitions** (stress matched across each):

| # | transition | check |
|---|---|---|
| 1 | PASS→BUCKLING (onset) | at onset `g_law = 0` (v1 is C0 at onset, `.cpp:264`); `Cg=Cg_law_reentry=0` ⇒ `g=0` ⇒ `σ=sBare` both sides. |
| 2 | BUCKLING→RESTRAIGHTEN (reversal) | Phase-1 at `eps=Ce_rev`: `σ = Cs_rev` = committed buckled stress. (Tangent kink `k_v1→k_rev` is an accepted reversal-class kink, like Steel02.) |
| 3 | Phase 1→Phase 2 (crossing) | `q=0 ⇒ B(0)=1 ⇒ σ = Cf_bare_cross − D0 = Cs_rev + k_rev·(Ce_cross_rs−Ce_rev)` = Phase-1 endpoint. |
| 4 | RESTRAIGHTEN→PASS (completion) | `q=1 ⇒ B(1)=0 ⇒ σ=sBare`; **C1-clean** (`B'(1)=0`). |
| 5 | RESTRAIGHTEN→BUCKLING (re-compress, **the crux**) | latch `Cg`=live deficit, `Ce_reentry=Cstrain`, `Cg_law_reentry=g_law(Ce_reentry)`. Next BUCKLING at `eps=Ce_reentry`: `g = min(g_law, Cg + 0) = Cg` ⇒ `σ = sBare − Cg` = the RESTRAIGHTEN value. **No stress jump — the v2-A break is eliminated.** |

**Energy.** `g ∈ [0, g_law]` ⇒ the output is bracketed `[σ_v1, σ_bare]`: the
wrapper can only *remove* stress relative to bare, never add it, so no spurious
injection is possible. The clean compress→elastic-unload→full-reload loop
encloses positive area (the only recoverable segment is the Phase-1 elastic line;
the reload path lies below bare and rejoins it from below) ⇒ dissipative.

**Convergence.** `g` is clamped `0 ≤ g ≤ g_law(ε)` *every* step, and `g_law` is
itself bounded (DM floors `σ` at `−0.2 f_y`, `.cpp:289/302-305`). `g` is **reset**
to a law value or to 0 each episode, never integrated — there is no per-cycle
accumulator, hence no ratchet/runaway. The deepest-buckle envelope `(Ce_bow,
Cg_bow)` only ever *deepens monotonically*; re-buckling to the same depth gives
the same gap (idempotent envelope).

### 9.3 State machine — committed-latch / commit-time transition (path-independent)

Add (mirroring every v1 state as a **committed+trial pair**) the following
**new committed fields** (each gets a `T…` trial twin; `int` enums stored as
`double`, exact cast on read like `model` at `.cpp:540`):

| # | field | meaning | new? |
|---|---|---|---|
| 1 | `Cbranch` | `{PASS=0, BUCKLING=1, RESTRAIGHTEN=2}` regime of the converged step | new |
| 2 | `Cstrain` | committed engineering strain (v1 has only `Tstrain`) — needed for `dEps` direction | new (review [BUG]) |
| 3 | `Cstress` | committed wrapper output stress (the `s_rev` continuity-anchor source) | new (review [BUG]) |
| 4 | `Cg` | **the residual additive deficit `g ≥ 0` below bare** — the single crux carrier | new |
| 5 | `Ce_rev` | reversal-latched strain (`= Cstrain` at BUCKLING→RESTRAIGHTEN) | new |
| 6 | `Cs_rev` | reversal-latched stress (`= Cstress` at that transition) | new |
| 7 | `Ce_cross_rs` | frozen episode anchor crossing `= CmaxTenStrain − CmaxTenStress/E` | new |
| 8 | `Cf_bare_cross` | `σ_bare(Ce_cross_rs)` latched at the Phase-1→2 crossing commit (exact `D0`) | new |
| 9 | `CL_rs` | re-straightening span (§9.4) | new |
| 10 | `Ce_reentry` | strain at the RESTRAIGHTEN→BUCKLING re-compression (offset anchor) | new |
| 11 | `Cg_law_reentry` | `g_law(Ce_reentry)` cached at re-entry (the offset reference) | new |
| 12 | `Ce_bow` | deepest committed buckle strain (envelope bound) | new (optional) |
| 13 | `Cg_bow` | `g_law(Ce_bow)` (envelope amplitude the offset re-merges to) | new (optional) |

Kept from v1 unchanged (`.cpp:116-119`): `CmaxTenStrain`, `CmaxTenStress`,
`CfStarL`, `CfStarLcross`. Fields 12-13 are **bounds, not stress-law carriers**:
the `min(g, g_law)` clamp re-merges the offset onto the v1 envelope without them,
so an implementer may drop them (10 new instead of 12) — but keeping them makes a
mid-re-buckle serialized state reconstruct exactly, so **keeping them is
recommended**. `k_rev` is **not** stored (recomputed `= getInitialTangent()`).

**Idempotence rule (the path-dependence fix, review [BUG]):** `setTrialStrain`
decides the *effective* branch for the current trial from the **committed**
branch + committed latches + the current trial strain, and computes σ/tangent —
**recomputed from scratch every call, never incrementally latched on the trial
path**. Every quantity feeding the law is committed or a stateless re-call:
`g_law` re-calls `applyBucklingDM/GA` (its `fStarL` cache keys on `e_cross`, a
memoized pure function — not a path latch); `D0` reads only committed
`{Cf_bare_cross, Cs_rev, Ce_rev}`; the only trial-path writes are
`Tstrain/Tstress/Ttangent/Tr`. The **authoritative branch+latch transition is
persisted only in `commitState`** from the converged `Tstrain` vs `Cstrain`.

> **Scope of the idempotence claim (review precision):** it covers the new
> RESTRAIGHTEN state and its committed latches. The **retained v1 anchor update**
> (`TmaxTenStrain`/`TmaxTenStress`, `.cpp:384-389`) is a monotone trial-path
> latch that is *intentionally* **not** re-seeded from committed on entry — the
> PASS→BUCKLING onset gate and the v1 BUCKLING read the live trial anchor and are
> path-dependent across Newton sub-iterates **exactly as in v1**. This is a
> tolerated pre-existing property *required* to stay so by the v1-bit-identical
> mandate (B2/B4a). The RESTRAIGHTEN blend is immune because `Ce_cross_rs` is
> frozen from the committed anchor.

Transitions (`tol = 1e-12 + 1e-8·eY` dead-band; `dEps = eps − Cstrain`;
`e_cross_live = CmaxTenStrain − CmaxTenStress/E`; **rows within a `from` group
are first-match-wins, top to bottom**):

| from | condition (per law) | to | latched at commit (all from committed state) |
|---|---|---|---|
| PASS | DM: `eps − e_cross_live < −eY`; GA: `eps < e_cross_live` | BUCKLING | — (`Cg` starts at `g_law` since `Cg=Cg_law_reentry=0`) |
| PASS | else | PASS | — |
| BUCKLING | `dEps ≤ −tol` or `|dEps|<tol` (compress/hold) | BUCKLING | if `eps<Ce_bow`: `Ce_bow=eps, Cg_bow=g_law(eps)` (monotone deepen) |
| BUCKLING | `dEps > +tol` (turns tensile) | RESTRAIGHTEN | `Ce_rev=Cstrain`, `Cs_rev=Cstress`, `Ce_cross_rs=e_cross_live` (frozen), `CL_rs` per §9.4 |
| RESTRAIGHTEN | (entering Phase 2 first time, `eps≥Ce_cross_rs`) | RESTRAIGHTEN | `Cf_bare_cross = theBar->getStress()` (the crossing bare stress) |
| RESTRAIGHTEN | `q ≥ 1` (fully straightened) | PASS | `Cbranch=PASS`, `Cg=0`, clear episode latches |
| RESTRAIGHTEN | `dEps < −tol` before `q=1` (re-compress) | BUCKLING | `Cg`=live deficit, `Ce_reentry=Cstrain`, `Cg_law_reentry=g_law(Ce_reentry)` — **C0-continuous** |
| RESTRAIGHTEN | else | RESTRAIGHTEN | — (Phase 1/2 split on frozen `Ce_cross_rs`) |

- **v1 BUCKLING math byte-untouched:** the new code only *dispatches* on the
  branch; `BUCKLING` calls the existing `applyBucklingDM/GA` verbatim into scratch
  (`σ_v1, k_v1`) and the reconciliation collapses to `g = g_law` whenever there
  was no prior straightening. The v1 anchor update + `fStarL` cache invalidation
  stay first in `setTrialStrain`, exactly as now ⇒ **B4a is bit-identical by
  literal code reuse, not numerical coincidence** (do **not** add a separate PASS
  early-return — route PASS/pre-onset through `applyBucklingDM/GA`, which already
  passes through at `.cpp:264/341`).
- **RESTRAIGHTEN is law-agnostic** (`g = D0·B(q)` over the wrapper's own
  `Cs_rev`/`D0`), so one branch serves DM and GA; only the PASS→BUCKLING *onset*
  is law-specific. A GA-buckled C0 assertion at the crossing (B4-GA) is the only
  thing that *proves* this rather than asserts it.
- **Two crossings coexist by design:** the *live* `e_cross_live` drives the
  BUCKLING onset; the *frozen* `Ce_cross_rs` drives the active RESTRAIGHTEN blend.
  Correct, but a subtle invariant — document inline in the code.
- **Implementation checklist — field-wiring table (every new field × every call
  site).** A new committed field omitted from `getCopy` is the wrapper-state-loss
  footgun one level deeper (silent, no compiler error; v1 `getCopy` hand-copies
  12 members, `.cpp:451-462`). Each of the 13 fields above (+ trial twin) must be
  wired through **all** of: value-ctor inline-init · no-arg-ctor inline-init ·
  `getCopy` manual copy · `commitState` T→C · `revertToLastCommit` C→T ·
  `revertToStart` zero (`Cbranch=PASS`, `Cg=0`) · `sendSelf`/`recvSelf` index
  (§9.6). Additional call sites the v1 checklist omitted:
  - **`Print`** (`.cpp:568-596`) — dump `branch` + latches in both text and JSON
    so the path-dependent machine is inspectable.
  - **`-restraighten` parser** (§9.4) — new flag in `OPS_LadrunoRebarBuckling` +
    a `restraightenMode`/`restraightenC` member + optional `setParameter` id 3.
  - **`Tr` on RESTRAIGHTEN** — define `Tr = (sBare≠0)? σ/sBare : 1.0` (documented
    as *effective* reduction, since the branch is additive not multiplicative),
    set every RESTRAIGHTEN trial; `Tr=1.0` on the Phase-1 sub-branch — so
    `getResponse("reduction")` stays meaningful instead of reporting a stale `r`.

### 9.4 Calibration of `L_rs` — both forms ship (DECISION 3)

`L_rs` scales with the plastic bow accumulated in compression, modulated by
slenderness. The **`-restraighten` flag selects the form**; both feed the *same*
smoothstep Phase 2, only `CL_rs` differs. `CL_rs` is computed **once** at the
BUCKLING→RESTRAIGHTEN latch and frozen for the episode.

**Option (1) — `-restraighten c <value>` (default, `c=1.0`):**
```
CL_rs = max( c·|Ce_rev − Ce_cross_rs| , eY )
```
The `eY` floor prevents `D0/CL_rs` blow-up / instant `q`-saturation when the
reversal is just past onset (`|Ce_rev−Ce_cross_rs| → 0`).

**Option (2) — `-restraighten lambda` (slenderness-set span):**
```
CL_rs = f(λ)·eY ,    f(λ) = 0.5 + 0.10·max(0, λ − 5) ,   clamped to [0.5, 6.0]
```
- **Monotone in λ** — bow depth (hence the straightening span) grows with
  slenderness `λ = s/d`.
- **Sane threshold:** at `λ ≈ 5` (where the onset magnitude saturates,
  `eStarMag` floored at 7, `.cpp:271`) `f(5)=0.5` ⇒ `CL_rs = 0.5·eY`: a
  barely-buckled bar straightens within half a yield strain.
- **Slope tied to the DM onset family:** onset deepens by `~2.3·kfac·eY` per unit
  `λ` (`kfac = √(fy/E·2000)`, `.cpp:269-270`); `0.10·eY` per unit `λ` is a
  conservative `~2%` fraction of that excursion, `O(eY)`-scaled.
- **Upper clamp 6.0** caps `CL_rs` at `6·eY` so a very slender bar does not give
  an unphysically long, near-flat reload (which would under-predict dissipation).

> **Floor magnitude (review precision):** the `eY` floor caps the Phase-2 linear
> term at `D0·E/fy`. In a deep-buckled reversal `s_rev` sits at the `−0.2 f_y`
> residual while `σ_bare` at the crossing is `O(f_y)`, so `D0` can reach `~1.2 f_y`
> and the reload tangent `~ kBare + 1.2·E ≈ 2.2 E` — the floor prevents blow-up
> but does **not** cap near `E`; the reload stays markedly stiffer than bare (the
> smoothstep distributes this over the span, and IMPL-EX mitigates implicit
> thrash). For DM the floor only binds near the shallow `−eY` onset; a deep
> `eStar` reversal has `|Ce_rev−Ce_cross_rs| ≳ 7·eY` and never approaches it — so
> the degenerate near-onset case is primarily a GA concern.

The DM-2002 cyclic-companion fit of `f(λ)` (and of the smoothstep itself) is the
deferred calibration target (§9.6); option (2) exists precisely so a user can fit
it.

### 9.5 Tests (B4 battery) — loop-shape + limits + brackets + structural Newton

**Harness (locked):** all B4 cyclic rows use the **committing `_drive`** harness
(`tests/...:49`) — the fresh-material `_oneshot` (`:62`) cannot express a reversal
and never exercises `commitState`/branch dispatch, so it is **forbidden** for any
B4 row. The existing `_dm_ref`/`_ga_ref` ports (`:93/:246`) are **virgin-only**
(`e_cross=0`, `fsup=0`) and **cannot** serve as the Phase-2 oracle; add a
closed-form `_restraighten_ref(eps, σ_bare(eps), s_rev, e_rev, e_cross_rs, L_rs,
k_rev)` mirroring §9.2 (both phases, smoothstep) so B4b/B4d/B4e check a falsifiable
formula, not just a trend.

| # | Test | Oracle |
|---|---|---|
| B4a | **reduce-to-v1 (frozen golden)** | capture v1 numeric output of B2 (`λ∈{6,8,12}`), B2b (the 80-step ramp) and GA2 as **golden array literals committed into the test** (no live v1 binary exists post-v2); assert v2 reproduces them to `rel≤1e-12` through **both** `_oneshot` *and* committing `_drive`, DM and GA; and that a no-reversal `_drive` leaves `branch ∈ {PASS,BUCKLING}` every step (the new commit path is exercised but proven inert). |
| B4b | compression→full tension reload rejoins bare | at `q≥1`, σ/tangent == bare (PASS restored); C0 at reversal, crossing, completion; matches `_restraighten_ref`. |
| B4c | elastic-unload segment recovers at `k_rev` | Phase-1 tangent `≈ k_rev = getInitialTangent` (guards: no flat/zero-stiffness unload). |
| B4d | re-straightening lag is **bracketed** | strictly inside `(Ce_cross_rs, Ce_cross_rs+L_rs)`: v2 σ is **strictly below `σ_bare(eps)`** AND **strictly above the extended Phase-1 elastic line** `s_rev + k_rev·(eps−e_rev)` (the "above v1 immediate-rejoin" bound is **dropped** — v1 already == bare past onset, so that bracket is empty); enclosed loop area `∮σ dε > 0` for the clean loop only. |
| B4d' | RS overlay (informational, NOT asserted) | overlay a `ReinforcingSteel -DMBuck` cyclic loop; eyeball shape only. |
| B4e | consistent-tangent FD on RESTRAIGHTEN | analytic vs central FD, probed strictly **inside** Phase 2 (away from `q=0,1` and the Phase-1/2 seam). |
| B4f | partial reload → re-compress (RESTRAIGHTEN→BUCKLING) | **C0-CONTINUITY** at re-entry: assert `|Δσ|` at the re-compression instant is `≈ 0` (machine, pinned as a regression — the v2 residual-gap closes the v2-A break); plus (1) `|σ| < |σ_bare|` once past the live onset (knock-down resumed), (2) post-re-compression stress == residual-gap formula at the live anchor. |
| B4g | **structural Newton across the seams** | drive a `Truss`/`DisplacementControl` path through real `analyze()` that commits a BUCKLING step then reverses into RESTRAIGHTEN, plus DC steps straddling `q=0` and `q=1`; assert `analyze()==0` (Newton converges) or document the step-cut requirement. |
| B4h | **v2 serialization round-trip** | commit one model in BUCKLING and one in RESTRAIGHTEN (`0<q<1`), round-trip via `FE_Datastore` (mirror B7); assert `branch` + all latches survive `recvSelf` and the next `setTrialStrain` reproduces σ/tangent bit-for-bit; assert a v2 `recvSelf` of v1-width data hard-errors per the §9.6 version check. |
| B4-GA | **GA-cyclic** (B4b/B4d/B4e with `-model ga`) | a prior tension excursion (`e_cross<0`), compress past `e_cross` into BUCKLING, reverse, reload; same asserts **plus an explicit C0 assertion at `q=0` for a GA-buckled reversal** — the one line proving RESTRAIGHTEN is law-agnostic. |
| B4i | near-onset floor | reverse one increment past onset (GA: past `e_cross`; DM: past `−eY`); assert `CL_rs` clamped to `eY`, Phase-2 tangent finite/bounded (`≈ kBare` since `D0→0`), `q` monotone in `[0,1]` with no overshoot. |

### 9.6 Compatibility / risk

- **Default-off identity preserved**: v2 changes only the post-reversal branch;
  `-lsr 0` and monotonic paths are bit-identical to v1 (B4a guards it).
- **Serialization version in the ID, NOT `data(0)` (review must-fix):** `data(0)`
  already holds `model` (`.cpp:489`) — a version there would silently corrupt the
  model selector. **Locked:** widen `dataID` from `Vector(3)` to `Vector(4)` and
  write `dataID(3) = schemaVersion` (v1 implicitly 1; v2 writes 2). `recvID` is
  read *before* `recvVector`, so `recvSelf` branches on `dataID(3)` to choose the
  `Vector` width (11 for v1; `11 + N_new` for v2 — `Vector(24)` keeping all 13
  new fields, or `Vector(22)` dropping `Ce_bow/Cg_bow`) and zero-inits the v2
  fields when reading a v1 stream (`Cbranch=PASS`, `Cg=0`, …). A v2 `recvSelf`
  that reads an unsupported version **hard-errors** (`opserr << "...unsupported
  serialization version N (need 2)"; return -3`) rather than relying on a
  channel-dependent width-mismatch. The `Vector` sits strictly between the ID send
  and the trailing `theBar->sendSelf`, so growing it leaves the tag/classTag/dbTag
  handshake and the nested `theBar` send/recv ordering untouched. `int branch`
  stored as `double` is exact (cast on read like `model`).
- **Both laws**: one law-agnostic RESTRAIGHTEN branch (DM + GA).
- **Honest divergence**: physically-motivated reformulation, not an RS port; the
  RS overlay (B4d') stays informational, never asserted.
- **Softening/localization (carried forward from §7):** the RESTRAIGHTEN branch
  adds a second negative-then-stiff regime on top of the v1 post-buckling
  softening — the same implicit-Newton stall and fiber-section localization
  hazards apply; mitigation = **IMPL-EX** (structure-only hook) + step-cut, and
  the crack-band/regularization caveat ([[project_bezier_charlen]]) applies
  unchanged.
- **Known calibration debt** (documented, not silent): the smoothstep `B(q)` and
  both `L_rs` forms are physically-motivated and scale-consistent with the DM
  onset coefficients but **not yet fit** to the DM-2002 cyclic companion;
  `-restraighten` exposes both forms for that fit. The re-buckle re-merge shape
  (`g` re-grows by the v1 increment, clamped onto `g_law`) is C0 and bounded but
  uncalibrated — same shape-debt class. Deep multi-interrupt histories are
  individually C0/bounded (proven) but the cumulative loop shape is validated only
  by B4f/B4g — an acknowledged validation gap.

---

**Deferred:** steel-profile plate-buckling law (route c — a *different* material,
not this one); interaction-calibration guidance with corotational global buckling;
DM-cyclic `L_rs`/smoothstep fit against the Dhakal–Maekawa 2002 cyclic companion.
