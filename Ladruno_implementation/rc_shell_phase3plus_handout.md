---
title: "RC shell stack — forward handout: Phase 3 and beyond"
project: Ladruno
type: developer handoff / forward plan
status: planning (cyclic material physics COMPLETE; this is the roadmap for what remains)
owner: nmora
related:
  - "[[19_ladruno_rc_shell_adr]]"
  - "[[LadrunoRCConcrete_guide]]"
  - "[[LEDGER_quirks]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[LogStrain_guide]]"
  - "[[LadrunoBrick_reference]]"
tags:
  - rc-shell
  - roadmap
  - handoff
  - tension-stiffening
  - finite-strain
  - solid-shell
  - punching
---

# RC shell stack — forward handout: Phase 3 and beyond

This is the **developer handoff** for everything left in the RC shell stack after the cyclic
constitutive core landed. It turns the terse Phase-3/4/5 bullets of the source ADR
([[19_ladruno_rc_shell_adr]]) into actionable plans: for each remaining item — *goal, what to
build, what to reuse, the decisions already locked vs still open, the acceptance gates, the
gotchas to watch, and the rough effort*. Read the ADR for the full reasoning; read this to
**pick up and execute**.

---

## Where we are (the baseline this builds on)

**Shipped and merged to `ladruno` (the cyclic material physics is COMPLETE):**

| Capability | Flag | Phase / PR |
|---|---|---|
| MCFT compression softening (strength axis) | `-beta` (+`-lublinerReduced`) | 1 (#155/#192) |
| Aggregate-interlock shear bound (`v_ci,max`) | `-interlock` | 2a (#239) |
| Cyclic friction-slip | `-cyclic` | 2b.1 (#245) |
| Shell-element integration + serialization | — | 2b.2a (#246) |
| Second crack (X-cracking) + Archard wear | `-xcrack` | 2b.2b (#253) |
| Crack-shear retention curves | `-shearRetention {mcft\|const\|dsfm\|rots}` | 2b.2c.1 |
| Rigid-rotation objectivity gate | — (test) | 2b.2c.2 |
| Crack-closure-on-normal (verified already correct) | — (test) | 2b.2c.3 |
| Meshed-wall **quasi-static explicit** cyclic validation | — (Zone-B) | 2b.2c.4a (#266) |
| **IMPL-EX** robustness | `-implex` | 4a (#263) |

**The solver story is settled:** cyclic softening RC walls run as **quasi-static explicit**
(`CentralDifferenceLadruno`) — see [[LadrunoRCConcrete_guide]] §7.4 for the recipe and the four
traps (element mass via `-rho`; `ASDShellQ4` has no `dt_cr` → manual wave-speed bound; no
`equalDOF`; mass-proportional damping only). `-implex` is the implicit-path alternative.

**What remains** is one validation item + three forward phases. None of them is blocked by the
others except where noted.

---

## The remaining-work map

| Item | Kind | Blocks on | Effort | Recommended order |
|---|---|---|---|---|
| **V — Tran–Wallace experiment calibration** | validation (Zone-B) | *experimental data from the user* | M | when data is available |
| **Phase 3 — tension stiffening + directional `lch`** | material (small) + maybe 1 vanilla edit | nothing | S–M | **do first (code-only)** |
| **Phase 4b — finite-strain view** (`LadrunoRCFiniteStrain`) | material view | reuses shipped `LogStrain`/`FiniteStrainNDMaterial` | M | after Phase 3 |
| **Phase 5 — `LadrunoSolidShell`** (33020) | NEW element | wants Phase 4b (the 3D material view) | **L (the big one)** | last |

`S`≈a focused slice, `M`≈a multi-day phase, `L`≈a multi-PR element with net-new EAS code.

> [!tip] Recommended sequencing
> **Phase 3 → Phase 4b → Phase 5**, with **V (Tran–Wallace)** slotted in whenever the
> experimental loops become available (it needs your input, not more code infrastructure). Phase
> 3 is the cleanest self-contained next slice; Phase 5 is the heavy lift and genuinely wants the
> Phase-4b 3D material view first.

---

## V — Tran–Wallace squat-wall experiment calibration (deferred 2b.2c.4 item 2)

**Goal.** The ADR's *primary* squat-wall gate: match a **named physical experiment**'s cyclic
loops, asserting **pinching shape + cumulative hysteretic energy**, not just "it runs." 2b.2c.4a
proved the solver + the panel mechanism; this proves the *numbers*.

**What's needed (the gating inputs — from you, not code):**
- A specimen: **Tran–Wallace RW-A20-P10** (the ADR's named target) or a PEER squat-wall with
  published cyclic data.
- Specimen geometry + axial-load ratio + the **reinforcement layout** (smeared web ratios/angles
  + boundary-element bars).
- **Digitized measured hysteresis loops** (lateral force vs top drift) to compare against.

**What to build (once the data exists):**
- A graded mesh of the actual wall (`gmsh`/`apeGmsh`; boundary elements as discrete
  `PlateRebar(LadrunoRebarBuckling)` layers, smeared web steel inside the kernel via `-rho`-style
  homogenization — see the composite-`ε1` note below).
- Drive the measured loading protocol under the quasi-static explicit recipe ([[LadrunoRCConcrete_guide]] §7.4).
- Assert: peak-capacity envelope within ~10–15%; **cumulative dissipated energy** within a stated
  band; **pinching shape** (a reload-stiffness / waist metric) tracks the measured trend.

**Reuse.** The existing harness `tests/_testbed/rc_wall_harness.py` (the structured grid + layer
build) and the Zone-B explicit recipe in `tests/test_ladrunoRCConcrete_wall.py`.

**Gotchas.** Calibration is **iterative** — `β`/`Kc`, the retention mode (`-shearRetention`), and
the wear knobs (`-degKappa/-degSlipRef`) are the dials; the boundary-rebar buckling
([[LadrunoRebarBuckling_guide]]) drives the strength-degradation tail. Don't over-fit; report the
calibrated parameter set explicitly. **Composite `ε1`:** MCFT softening + tension stiffening are
defined on the *reinforced* average principal tensile strain — the smeared web steel must be
homogenized **inside** the kernel (or its strain shared), else the concrete over-softens (see
Phase 3 + ADR risk "MCFT needs composite `ε1`").

---

## Phase 3 — Tension stiffening + crack-band / `lch` hardening

**Goal.** (a) Average tensile stress *between* cracks from bonded reinforcement (tension
stiffening), for slabs and distributed-reinforcement walls; (b) make in-plane softening
**size-objective on an inclined-crack (rotated) mesh**, not just an axis-aligned one.

**What to build.**
1. **Tension-stiffening plateau** as an opt-in mode. Two candidate laws (ADR §"Tension
   stiffening"): MCFT/Bentz `σ_t_avg = ft/(1+√(c·ε1))` and Collins–Mitchell
   `α1·α2·ft/(1+√(500·ε1))`. Currently bakeable into the `-Te/-Ts` backbone (the documented
   default); Phase 3 promotes it to a real `-tensStiff {vc|cm}` knob (the ADR reserves the name)
   driven by the **composite `ε1`** (reinforced, not bare-concrete) so it doesn't double-count
   with a separate steel layer.
2. **`lch` resolution — pick and own ONE (ADR decision D5):**
   - **Option A (default, recommended):** accept the element's scalar in-plane `lch` via
     `getCharacteristicLength()` (already encodes `ASDShellQ4`'s EAS `/2`); treat through-thickness
     crush as the layer thickness; **document that out-of-plane bending-crack energy is not
     through-thickness size-objective on the director shell** (that's Phase 5's solid-shell).
   - **Option B (only if inclined-crack objectivity must be exact):** a small **Ladruno-tagged,
     `LEDGER_vanilla_files`-recorded** edit to `LayeredShellFiberSection`/`ASDShellQ4` plumbing a
     per-direction `lch` to the layer — and **drop the "zero vanilla edit" claim** for that case.

**Reuse.** The spine's `buildBackbone`/`adjust()` E-consistency machinery; the existing
`getCharacteristicLength()` channel; `tests/_testbed/rc_shell_ref.py` oracle.

**Decisions — locked vs open.** Locked: tension stiffening rides the **composite `ε1`**; flags-off
reduces to the bare fracture-energy backbone (baseline-identical). **Open:** which `lch` option
ships (A vs B), and the crush-band floor `lch_t ≥ max(t_i, k·d_agg)`.

**Acceptance gates (ADR Phase 3).**
- In-plane **mesh-objectivity on a ROTATED (inclined-crack) notched panel**: peak ~1–3%, energy
  ~3–5% across a 2× refine — **calibrated on the elastic shell patch first**, NOT imported from the
  brick numbers.
- A **loud failure** if the scalar `lch` fallback is silently used in a softening run (no silent
  default — ADR D5).
- Tension-stiffening: average tension above the bare softening curve, reducing to it as `ε1` grows;
  reduce-to-baseline with the flag off.

**Gotchas / risks.**
- **`lch` is a single in-plane scalar on this seam** and `ASDShellQ4` halves it under EAS. A single
  scalar mis-regularizes ~45° struts by up to √2 (ADR risk). Option A documents this away; Option B
  fixes it at the cost of a vanilla edit.
- **`lch_t = t_i` is NOT a safe compression default** — it can trigger snapback that the material's
  `gmin` clamp silently discards; floor it at a physical crush-band width.
- **Composite `ε1` consistency** must be reconciled everywhere (softening *and* stiffening use it).

**Effort:** S–M. Tension stiffening is a small additive mode; the `lch` work is S for Option A, M
for Option B (the vanilla plumb + ledger).

---

## Phase 4b — Finite-strain view (`LadrunoRCFiniteStrain`)

**Goal.** A large-strain RC material for the solid-shell host and for `-geom finite` runs;
symmetric-solver, explicit-friendly cyclic. This is the third `getType()` view of the existing
class (the small-strain `3D`/`PlaneStress`/`PlateFiber` views already ship).

**What to build.** `LadrunoRCFiniteStrain` as a `FiniteStrainNDMaterial` subclass on **classTag
33015** (the ADR-reserved finite view — confirm no collision at build time):
`setTrialF(F)` → Hencky log-strain → the existing kernel `returnMap3D` → Cauchy + spatial tangent.
Carry IMPL-EX through (clamped `ε1`/`β`, the implex-vs-implicit error monitor + step-cut).

**Reuse (all verified present on `ladruno`):** `FiniteStrainNDMaterial` base,
`LogStrainNDMaterial` ([[LogStrain_guide]]), and the **`LadrunoJ2Finite` pattern**
([[09_finite_strain_material_wrapper]]) — this is the exact same lift that gave J2 its finite view.
The small-strain kernel is reused verbatim; only the `F`-seam adaptor is new.

**Decisions — locked vs open.** Locked: the seam is `setTrialF` → log-strain → kernel (NOT routing
`F` through the section — the section owns no Gauss-point `F`; ADR D2). **Open:** whether to reuse
the shipped `LogStrainNDMaterial` directly as a wrapper vs a dedicated subclass that calls the
kernel inside (the J2 family did the latter for the co-rotating-state case).

**Acceptance gates (ADR Phase 4b).**
- IMPL-EX `O(dt)` on a **smooth proportional** path (excluding damage-activation/closure steps).
- A cyclic **energy-balance band** on the wall.
- **Rigid-rotation-of-a-cracked-point xfail (dSNPO §14.11):** the directional crack/interlock state
  is objective only for the corotational *element* route, NOT the `setTrialF` *material* view under
  large rotation. Pin this as a strict xfail (same boundary as `LadrunoJ2Finite`'s kinematic
  hardening). The *isotropic-damage spine* IS objective under the wrapper; the *directional* state
  is the xfail.

**Gotchas / risks.**
- **`det F ≤ 0` guard** in the log-strain setTrialF (the J2-finite review added exactly this).
- The §14.11 boundary is a **framework limit, not a bug** — document it; don't try to "fix" it
  without a co-rotating-crack finite-native formulation (out of scope).

**Effort:** M. It is a well-trodden lift (the J2 family is the template), but the
finite-strain + directional-state + IMPL-EX interactions need careful testing.

---

## Phase 5 — `LadrunoSolidShell` (33020) — through-thickness host (the big one)

**Goal.** The one genuine **elemental** blind spot a director shell cannot represent:
through-thickness `σ33` for **punching / bearing / 3D-stress crush**. A narrow specialist, **not**
a co-equal flexural host.

**What to build.** An 8-node, 3-translational-DOF solid-shell, classTag **`ELE_TAG` 33020**:
- **Genuine state-dependent EAS-on-`E33`** — persistent committed `alpha[nEAS]`, per-Newton
  condensation with the consistent **damaged** tangent, serialized send/recvSelf. This is **net-new
  Simo–Rifai code**, NOT a port.
- ANS/MITC transverse shear + ANS membrane/trapezoidal cures (anti-locking).
- Selectable **multi-layer / Gauss–Lobatto `n_z`** (a single-layer ~2-GP solid-shell cannot resolve
  cracked-RC `σ(z)` or the migrating neutral axis).
- Directional / projected `lch`; a `cond(S)` corotation guard for thin shells.

**Reuse.** `SolidTransformation` linear/corot/finite ([[solid_transformation_wrapper]]); the
**Phase-4b `LadrunoRCFiniteStrain`** material view (the reason 4b comes first); `LadrunoBrick`'s
damage-scaled `Kstab` **as stabilization only, re-tuned for thin shells** — NOT as the EAS template.

> [!warning] The two FATAL-avoided traps (ADR risk register — do not re-cost these)
> 1. **`LadrunoBrick`'s "EAS" is NOT a reusable template.** It is the SSPbrick stabilized
>    single-point scheme that condenses enhanced modes ONCE at setup with the *initial elastic*
>    tangent into a CONSTANT `Kstab` — it cannot soften. Genuine EAS-on-`E33` for softening RC needs
>    persistent `alpha`, per-Newton condensation with the consistent damaged tangent, and serialized
>    state. **Budget it as net-new code.**
> 2. **A single nodal-cloud corotation `R` is ill-conditioned for thin shells** (`cond(S) ~ (L/t)²`),
>    so a single `R` cannot remove the through-thickness-varying bending rotation. Ship as a
>    moderate-rotation / punching specialist with a `cond(S)` guard; thin large-rotation routes to
>    `-geom finite` or a shell-aware corotation.

**Decisions — locked vs open.** Locked: solid-shell is a **punching/bearing specialist**, not a
co-equal flexural host (benchmark its moment-curvature vs `LayeredShellFiberSection` before any
flexural claim). **Open:** shell-aware corotation vs `-geom finite` for thin large-rotation;
whether `SolidTransformationCorot` needs its two deferred geometric-tangent terms (validate on a
**softening snap-back**, not the elastic pinched-cylinder, before softening-RC use).

**Acceptance gates (ADR Phase 5).**
- Element correctness: elastic **patch** + **pinched-cylinder** + **Scordelis–Lo**.
- A **softening snap-back** (mesh-objective dissipation; Newton/arc-length convergence with the
  incomplete-geometric + secant tangent).
- A slab **punching** benchmark (the headline — the thing the director shell structurally cannot do;
  the ADR documents the director-shell punching blind spot as an xfail that this discharges).
- An **EAS internal-mode growth-stability** check under post-peak softening.
- A documented, **validated** `ASDShellQ4`(6-DOF) ↔ solid-shell(3-DOF) rigid-link / `equalDOF`
  connection recipe (moment continuity across the DOF mismatch).

**Gotchas / risks.** The EAS internal modes can grow unstably under post-peak softening (gate it).
The 6-DOF↔3-DOF edge loses moment continuity (connection recipe is a deliverable). `dt_cr` for
explicit on this element — decide whether to implement a real per-element estimate (cf. the
`ASDShellQ4`-has-none gotcha in [[LEDGER_quirks]]) or document the manual bound.

**Effort:** **L — the big one.** Multi-PR: (P5.1) the element skeleton + ANS/MITC + patch tests;
(P5.2) the genuine EAS-on-`E33` + softening; (P5.3) punching benchmark + the connection recipe.

---

## Cross-cutting risks (apply across Phase 3–5)

These are the ADR open questions that recur; keep them in view:

- **Consistent tangent vs eigenprojector degeneracy.** `dP/dε ∝ (p_i⊗p_j + p_j⊗p_i)/(λ_i − λ_j)` is
  **singular at equibiaxial / hydrostatic** states (wall corners, biaxial slab bending). The current
  material uses the **fixed-projector secant** (omits `dP/dε`) — a deliberate v1 choice. If a future
  phase wants the full consistent tangent, it must add the **Miehe (1993) / dSNPO §A** coalescence
  regularization (blend to the symmetric average when `|λ_i − λ_j| < tol`), or it produces
  indefinite/NaN tangents in states walls/slabs visit routinely.
- **`lch` is a single in-plane scalar** on the director-shell seam (Phase 3 decision; Phase 5
  needs directional/projected `lch` for through-thickness).
- **MCFT needs composite `ε1`** — smeared web steel homogenized inside the kernel; discrete boundary
  bars are separate `PlateRebar` layers (recurs in V and Phase 3).
- **The shell path is `PlateFiber`-only** — `ASDShellQ4` consumes a *section*, not an `nDMaterial`;
  all wall shear flows through the order-5 condensed `PlateFiber`. Any new analysis is done in that
  5-component condensed setting, not a free-standing 3×3.

---

## Quick-start pointers for the next session

- **Source:** `SRC/material/nD/{LadrunoRCKernel.h, LadrunoRCConcrete.{cpp,h}}` (kernel is
  header-only, OpenSees-free, numpy-oracle-testable).
- **Oracle:** `tests/_testbed/rc_shell_ref.py`. **Batteries:**
  `tests/test_ladrunoRCConcrete_{material,shell,implex,objectivity,wall}.py`.
- **Build:** edit in the worktree → copy to the main checkout → **bump the file mtime** (Copy-Item
  preserves the old mtime → ninja skips the rebuild — see [[LEDGER_quirks]]) → `build.bat OpenSeesPy`
  via the PowerShell tool → confirm the `.pyd` timestamp advanced.
- **Schema:** currently **v3** (hard-checked). Any new committed state bumps it; keep send/recv
  counts balanced (the RC_DATA arithmetic).
- **Build-control:** every phase updates the three ledgers (`LEDGER_implementations`,
  `LEDGER_quirks`, `LEDGER_vanilla_files` *if* Phase-3 Option B touches a vanilla file) **in the same
  PR**, and the ADR phase subsection. New fork-authored files get the LADRUNO header stamp.
- **Reviews:** the project convention is implement → adversarial review → fix → PR. Each phase here
  warrants it (the §14.11 xfail, the EAS softening stability, and the tangent degeneracy are the
  likely review hot-spots).
