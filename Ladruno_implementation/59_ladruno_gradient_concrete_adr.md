---
title: "ADR 59 — Gradient/nonlocal regularization for the Ladruno concrete models (LadrunoGradientConcrete): scoping"
project: Ladruno
type: ADR / scoping
status: "DESCOPED → P0 research probe. PROPOSED architecture adversarially reviewed (5 red-team critics); 3 BLOCKING + ~9 MAJOR findings folded in. v1 reduced to a numpy/1-element probe, gated on (a) a named band-WIDTH-objectivity consumer and (b) reconciliation with ADR 32's shipped embedded-discontinuity alternative. NO code."
related:
  - "[[31_ladruno_concrete3d_adr]]"      # the local crack-band model this regularizes (LadrunoConcrete3D 33017)
  - "[[19_ladruno_rc_shell_adr]]"        # LadrunoRCConcrete (33015) — RC shell; P3 STRUCK (its lch seam is NOT the same)
  - "[[32_ladruno_dispbeamcolumn_regularization_adr]]" # SIBLING that REJECTED gradient + shipped an embedded discontinuity — must reconcile
  - "[[09_ladruno_brick]]"               # the host element the coupled u–ē element would fork from (fixed 24-DOF today)
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_quirks]]"
tags: [adr, regularization, gradient-damage, nonlocal, concrete, plastic-damage, mesh-objectivity, coupled-field, descoped]
updated: 2026-06-24
---

# ADR 59 — Gradient/nonlocal regularization for the Ladruno concrete models

**Status: DESCOPED to a P0 research probe.** The implicit-gradient architecture was scoped from the
regularization-seam review, then **adversarially reviewed by 5 red-team critics** (constitutive, FEM/
element, tangent/claims, strategy, numerics). The sweep found the architecture *implementable* (BrickUP
proves a 4-DOF/node coupled element assembles) but the proposal **over-sold on three load-bearing
claims** and **strategically premature**. This revision folds the findings in: v1 is reduced to a
**numpy + single-element probe**, and a full build is **gated** on a named consumer and an explicit
comparison against the fork's *own* shipped alternative (ADR 32). **No code.**

> [!warning] What the sweep changed
> The original one-liner — *"the constitutive change is a seam swap; the work is the element + one
> extra field"* — is **false**. The kernel has **no single `ε̃ → ω` seam** (B1); the off-switch is
> **not** byte-identical at `c=0` (B2); the nonlocal variable is **wrong** (B3). The element is also
> **not** the only real cost. See §9.

---

## 1. Driver & problem (tempered)

Both `LadrunoConcrete3D` (ADR 31) and `LadrunoRCConcrete` (ADR 19) regularize softening with
**crack-band** (Bažant-Oh): the kernel scales its softening thresholds by `lch`
(`eps_f=Gf/(ft·lch)`, `eps_fc=Gc/(fc·lch)`, kernel:989-990). This is **local** and robust. Its
documented limits: energy- not band-width-objective in compression (ADR 31 R5/R15); a snap-back
element-size cap `h ≤ 2EGf/ft²`; no length scale in the BVP.

A gradient/nonlocal scheme adds a real internal length `ℓ=√c` and a well-posed BVP — **but the sweep
established this is a *research* capability, not a structural-analysis need, and it inverts the fork's
mesh economy** (M4). Crack-band's coarse-mesh, energy-objective behaviour is the right *engineering*
default for the fork's seismic columns/walls; this ADR does **not** propose replacing it.

> [!note] Honest framing (B3 + M9)
> `lch` is a transparent numerical regularizer. `ℓ=√c` is **also** calibrated empirically (R3 itself
> says "calibrate ℓ to the FPZ ~few `d_a`" — a proxy), is **not independently identifiable**, and a
> *constant* ℓ spuriously broadens the damage band as damage grows. Gradient trades a bounded, honest
> fudge for an unbounded identifiability problem **plus** cost. It only earns its keep where
> band-WIDTH (not energy) objectivity is genuinely required — a consumer this ADR must name.

## 2. Where on the stack — decision (retained, claim corrected)

Implicit-gradient **over** nonlocal-integral stands (the sweep confirmed it: element-local coupling →
**parallel-friendly**, no Gauss-point ghost exchange across partitions; the strongest surviving
argument). But the **"rides the standard Domain for free"** framing is **struck (M5)**: OpenSees
`Node`s carry a **fixed ndf**; a 4-DOF (`u`+`ē`) node **cannot share a mesh** with a 3-DOF brick. The
whole gradient region must be ndf=4 with an `equalDOF` transition at its boundary — the known u-p/
`BrickUP` modeling cost. The architecture is *possible* (BrickUP proves it), not *free*.

## 3. The constitutive seam — CORRECTED (B1, B2, B3)

The original "swap local `ε̃` for nodal `ē`" is **mis-specified**. Verified against the kernel:

- **B1 (BLOCKING) — there is no single `ε̃ → ω` map.** ω is **dual-channel**: ω_t solved against
  `sigtMax = max(·, E·ε̃)` (a *stress*), ω_c solved against `sigcMax` = a **principal effective
  stress, not ε̃ at all** (kernel L1034/1036/1040-1043; the kernel's own comment notes `E·ε̃` is
  ft-scaled and "could never onset ω_c"). A single nonlocal `ε̃` regularizes the **tensile** channel
  at best and leaves **compressive crush — the headline motivation — fully local.** ⇒ either average a
  driver that genuinely feeds both channels, or carry **TWO** nonlocal fields `ē_t`, `ē_c` (**+2 DOF/
  node, not +1**, ~doubling the element cost).
- **B3 (BLOCKING) — wrong variable.** `ε̃ = equivStrainGeneral(σ̄)` is a **post-return-map effective-
  stress functional**; averaging it creates a plasticity↔damage feedback and the Jirásek spurious-
  damage-growth pathology. **Grassl's own nonlocal CDPM averages a strain-rate driver `κ_d`** (loading-
  gated, monotone, non-negative), **not** `ε̃(σ̄)`. The implicit-gradient field must source from `κ_d`,
  and the loading/monotone-max gates (currently local, kernel L1014/L1038) must move onto the
  **nonlocal** variable, or damage leaks into elastic-unloading/compressive neighbours (B3-adjacent).
- **B2 (BLOCKING) — `c=0` is NOT the off-switch.** At `c=0`, `R_ē` reduces to `M·ē = ∫N ε̃ dV` — an
  **L2 projection** of `ε̃` onto the (deliberately) ē interpolation space, **not** pointwise `ε̃` at the
  Gauss points. So `c=0` drives ω with a *smoothed* field ≠ the local kernel. **Byte-identity holds
  ONLY on a separate local-bypass host** (plain `LadrunoBrick`, `R_ē` not assembled, ε̃ fed pointwise).
  The P0/P1 "`c=0` byte-for-byte" gates as originally written are **unattainable** — they must be
  restated as a *regularizer-disabled* bypass, not `c=0`.

## 4. Architecture (revised)

### 4.1 Coupled element `LadrunoGradientBrick` (fork from `LadrunoBrick`, today fixed 24-DOF)
Per node `u`(3) **+** the nonlocal field(s). With B1, v1 must decide **one** common-driver field vs
**two** channel fields. `R_ē` is a **coercive Helmholtz** (`NᵀN + c∇Nᵀ∇N`) — **SPD, no inf-sup/LBB**:
the original "inf-sup / lower-order ē" language is a **category error (struck)**; **equal-order Q1/Q1
is the standard, accepted choice** for implicit-gradient damage. BC: homogeneous Neumann is the
conventional choice but is **known-imperfect** — it over-predicts damage in a boundary layer ~ℓ at a
notch tip (M8), exactly where the SENB gate cracks initiate → the gate must include a boundary-layer
diagnostic.

### 4.2 / 4.3 Tangent reuse — CORRECTED (M6, and one true positive)
- `K_ēu = −∫Nᵀ(dε̃_local/dε)B` **IS** a clean reuse of `det_deps` — *valid precisely because in the
  damage-plastic split σ̄(ε) is independent of ē* (damage reapplied downstream). `ε̃` never sees ω, so
  there is no damage-feedback term to add. (Confirmed by the sweep — the one genuine "half-built" piece.)
- `K_uē = ∂σ/∂ē` reuses only the **outer** IFT shell (the ω-solve coefficients + the `σ_t⊗dω/σ_c⊗dω`
  rank-1 assembly, ~30%). The **driver gradients** (`dkdt*`, `dnorm_deps`, `dxs_deps`, `dbc_deps`,
  `dwtw_deps`, kernel L1297-1381) are `d/dε` chains that must be **re-derived w.r.t. ē** (~70%), gated
  on the unmade B1/B3 choice of *what ē drives*. "Mostly re-wiring" is **withdrawn**.

### 4.4 Field scaling (NEW — M3)
`K_uu ~ O(E)` (1e4-1e5) vs `K_ēē ~ O(1)` — a **4-10 order block mismatch** ⇒ ill-conditioned monolithic
Jacobian that degrades the *mandated* unsymmetric direct solver (ADR 31 R13) and lets the global Newton
test **silently under-resolve `ē`**. v1 MUST choose `ē`-field non-dimensionalization (or row/col
equilibration) in the ADR, and gate on `cond(K)` + an independent ē-residual tolerance.

## 5. Reconciliation with ADR 32 (NEW — M1, the decisive strategic finding)
The **frame track of this same fork already faced this exact question and chose the opposite**: ADR 32
**explicitly rejected** nonlocal/gradient ("a different and heavier formulation… not pursued") and
**shipped** a two-tier answer — per-IP `lch` (Tier-1) + an **embedded strong-discontinuity hinge that
carries `Gf` directly, with no `lch`/ℓ to calibrate** (Tier-2, gate-complete 2D). The **3-D analog of
that** — an embedded crack/E-FEM band in `LadrunoBrick`, reusing the fork's **own EAS + static-
condensation machinery** (`LadrunoBrick::formEAStrue`) — is `lch`-free, band-width-objective, fork-
consistent, needs **no extra global DOF, no fine mesh, and composes with the explicit/IMPL-EX tiers**.
**This ADR cannot advance until it explains why gradient beats the embedded-discontinuity path the fork
already shipped for frames.** This is the single biggest open question (R8).

## 6. Phased roadmap — DESCOPED

- **P0 (the only committed phase) — research probe, no C++.** numpy 1-D screened-Poisson `ē` oracle +
  the **corrected** off-switch (a *local-bypass* byte gate, NOT `c=0`) + a **width-objectivity-vs-ℓ +
  ℓ-identifiability** study + a **spurious-initiation** check (no ω>0 in the far elastic field) on the
  `κ_d`-driven (not `ε̃(σ̄)`) field. Answers B1/B3/M9 cheaply and reversibly.
- **GATE TO P1 (all must hold):** (a) a **named deliverable that requires band-WIDTH (not energy)
  objectivity** — none exists today; (b) the §5 **ADR-32 reconciliation** shows embedded-discontinuity
  is insufficient for it; (c) an **alternatives table** (element-local driver averaging — which rides
  the existing seam with *no* new DOF; projected/eigen-`lch` per ADR 19 D5 Option B; 3-D embedded
  discontinuity) shows each is insufficient.
- **P1+ (NOT committed):** `LadrunoGradientBrick`, two-field-vs-one decision, the scaling fix, the
  mixed-ndf transition, the re-derived `K_uē`. **Tier-1 only, with a HARD refuse of `-implex`** (M2).
- **P3 (RC shell) — STRUCK.** ADR 19's `lch` is a single in-plane scalar with no per-layer/per-direction
  path and a through-thickness xfail; gradient-on-shell is a **separate hard problem**, not a corollary.

## 7. Alternatives (expanded — M1/M5/the §6 gate)
The original 3-sentence dismissal of nonlocal-integral is retained for *that* scheme, but the sweep
surfaced **cheaper alternatives the ADR must evaluate first**: (1) **element-local averaging of the
damage driver** over a Gauss-point stencil — rides the existing material seam, **no new DOF, no coupled
solve**; (2) **projected/eigen `lch`** (ADR 19 D5 Option B already scoped) — kills the directional-bias
defect cheaply; (3) **3-D embedded strong discontinuity** (the ADR-32-consistent, `lch`-free path).
Gradient is the heavyweight cure for *one* defect (true band-width objectivity) with **no named
consumer** — defects (heuristic size, snap-back cap, directional bias) are addressable by (1)+(2)
without a coupled field.

## 8. Risk register (rebuilt from the 5-critic sweep)

| # | Sev | Risk | Resolution / status |
|---|---|---|---|
| **R-B1** | **BLOCKING** | no single `ε̃→ω` seam; ω is dual-channel and ω_c is **not** ε̃-driven ⇒ one field leaves compression local | §3: average a true common driver, or **2 fields (+2 DOF)**; P0 decides |
| **R-B2** | **BLOCKING** | `c=0 ≠ byte-identical` — it's an L2 projection of ε̃ | §3/§6: off-switch = a **local-bypass host**, not `c=0`; restate the gate |
| **R-B3** | **BLOCKING** | averaging `ε̃(σ̄)` is the wrong variable (feedback + spurious growth) | §3: nonlocalize the **strain-rate driver `κ_d`**; move loading/max gates onto the nonlocal field |
| R8 | MAJOR | the fork's frame track **rejected gradient** for an embedded discontinuity | §5: reconcile before P1; the 3-D embedded-discontinuity may **supersede** this ADR |
| R-M2 | MAJOR | gradient is **Tier-1-only** — incompatible with Tier-2 IMPL-EX (returns secant≠Ceff; σ has no ē-dep) and Tier-3 explicit (inertia-less `ē`) ⇒ forfeits the robustness stack | re-rated from MINOR; **hard-refuse `-implex`** in v1; document Tier-3-incompatible or scope a staggered implicit-`ē` solve |
| R-M3 | MAJOR | u–ē block scaling 4-10 orders ⇒ ill-conditioned, under-resolves `ē` | §4.4: non-dimensionalize `ē` / equilibrate; `cond(K)` + ē-residual gate |
| R-M4 | MAJOR | **cost inversion** — needs `h≲ℓ/3` in the FPZ + 4 DOF/node; expensive where crack-band was coarse-cheap | gradient is research/high-fidelity, **not a drop-in**; keep crack-band default; P2 resolution-floor gate |
| R-M5 | MAJOR | mixed-ndf — "rides the Domain for free" false; whole region ndf=4 + equalDOF transition | §2: documented user modeling; not "free" |
| R-M6 | MAJOR | `K_uē` ~70% new derivation, gated on the B1/B3 driver choice | §4.3: withdrawn "mostly re-wiring" |
| R-M8 | MAJOR | Neumann BC unjustified ⇒ boundary-damage layer ~ℓ at the notch | justify as known-imperfect; boundary diagnostic in the gate |
| R-M9 | MAJOR | ℓ trades the `lch` fudge for an unidentifiable ℓ fudge + constant-ℓ broadening | P0 identifiability study before any C++ |
| (fix) | — | "inf-sup / lower-order ē" is a category error (ē block is coercive SPD) | §4.1: **equal-order Q1/Q1**; drop inf-sup language |

## 9. What the sweep retained (so the de-risking isn't lost)
- implicit-gradient **over** nonlocal-integral (parallel locality) — sound.
- `K_ēu = det_deps` reuse — genuinely valid (σ̄ is ē-independent).
- the **P0 probe** — cheap, reversible, and the right place to settle B1/B3/M9 before any C++.
- the architecture is *implementable* (BrickUP) — this is descoped on **strategy + correctness of the
  spec**, not feasibility.

## 10. Ledger / obligations
- No `LEDGER_implementations.md` row — no code, descoped.
- `LEDGER_quirks.md` — record: (a) crack-band `lch` is energy- not width-objective (the honest, bounded
  limit gradient was meant to cure); (b) CDPM2 damage is **dual-channel** and ω_c is a *stress*-driven,
  not ε̃-driven, variable (the fact that breaks a single-field nonlocal); (c) at `c=0` an implicit-
  gradient field is the **L2 projection** of the driver, not the pointwise driver.
