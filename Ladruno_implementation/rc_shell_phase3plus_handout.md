---
title: "RC shell stack — forward handout: Phase 4b and beyond (Phase 3 shipped)"
project: Ladruno
type: developer handoff / forward plan
status: planning — cyclic material physics + Phase 3 (tension stiffening + crack-band lch) COMPLETE; remaining = V (Tran–Wallace), 3b structural rotated-mesh gate (staged), Phase 4b (finite view), Phase 5 (solid-shell)
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
  - crack-band
  - regularization
  - finite-strain
  - solid-shell
  - punching
---

# RC shell stack — forward handout: Phase 4b and beyond

This is the **developer handoff** for everything left in the RC shell stack after the cyclic
constitutive core **and Phase 3 (tension stiffening + crack-band `lch` regularization)** landed.
It turns the terse Phase-4/5 bullets of the source ADR ([[19_ladruno_rc_shell_adr]]) into
actionable plans: for each remaining item — *goal, what to build, what to reuse, the decisions
already locked vs still open, the acceptance gates, the gotchas to watch, and the rough effort*.
Read the ADR for the full reasoning; read this to **pick up and execute**.

> [!note] What changed since this handout was first written (2026-06-18)
> **Phase 3 is SHIPPED.** `-tensStiff {vc|cm}` tension stiffening (3a, PR #273) and
> `-autoRegularization $lch_ref` crack-band regularization (3b, PR #277, **Option A** — scalar
> in-plane `lch`, zero vanilla edit) are merged. The Phase-3 section below is kept as a *shipped*
> record; the one remaining 3b item is the **structural** rotated-mesh objectivity gate (staged).
> The live forward work is **V → Phase 4b → Phase 5**.

---

## Where we are (the baseline this builds on)

**Shipped and merged to `ladruno` (cyclic material physics + Phase 3 are COMPLETE):**

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
| **Tension stiffening** (stress floor on the live principal axis) | `-tensStiff {vc\|cm}` (+`-tensStiffC/-tensStiffAlpha`) | **3a (#273)** |
| **Crack-band (Bažant–Oh) regularization** (mesh-objective softening energy) | `-autoRegularization $lch_ref` | **3b (#277)** |

Schema is now **v5** (hard-checked); all flags **default OFF ⇒ baseline-identical**.

**The solver story is settled:** cyclic softening RC walls run as **quasi-static explicit**
(`CentralDifferenceLadruno`) — see [[LadrunoRCConcrete_guide]] §7.4 for the recipe and the four
traps (element mass via `-rho`; `ASDShellQ4` has no `dt_cr` → manual wave-speed bound; no
`equalDOF`; mass-proportional damping only). `-implex` is the implicit-path alternative.

**What remains** is one validation item, one staged Phase-3 gate, and two forward phases. None is
blocked by the others except where noted.

---

## The remaining-work map

| Item | Kind | Blocks on | Effort | Recommended order |
|---|---|---|---|---|
| **V — Tran–Wallace experiment calibration** | validation (Zone-B) | *experimental data from the user* | M | when data is available |
| **3b-struct — structural rotated-mesh objectivity gate** | validation (Zone-B) | a softening-localization solver setup | S–M | optional, anytime |
| **Phase 4b — finite-strain view** (`LadrunoRCFiniteStrain`) | material view | reuses shipped `LogStrain`/`FiniteStrainNDMaterial` | M | **do first (code-only)** |
| **Phase 5 — `LadrunoSolidShell`** (33020) | NEW element | wants Phase 4b (the 3D material view) | **L (the big one)** | last |

`S`≈a focused slice, `M`≈a multi-day phase, `L`≈a multi-PR element with net-new EAS code.

> [!tip] Recommended sequencing
> **Phase 4b → Phase 5**, with **V (Tran–Wallace)** slotted in whenever the experimental loops
> become available (it needs your input, not more code infrastructure) and **3b-struct** whenever
> a localization-solver harness is convenient (the regularization *mechanism* is already proven at
> the material point). Phase 4b is the cleanest self-contained next code slice; Phase 5 is the
> heavy lift and genuinely wants the Phase-4b 3D material view first.

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

## Phase 3 — Tension stiffening + crack-band / `lch` hardening (SHIPPED)

Both halves landed; this is the *shipped record* (full detail in [[LadrunoRCConcrete_guide]] §4.7–4.8,
ADR Phase 3, and the `LEDGER_implementations` row).

**3a — Tension stiffening (`-tensStiff {vc|cm}`, PR #273).** A stress **floor** on the *live*
in-plane principal tensile axis `p1`: inject `Δ = σ_ts(ε1) − nᵀσn` along `p1` (only when `Δ>0`),
active **only post-crack** (`ε1 ≥ ε_cr`). `σ_ts = ft/(1+√(c·ε1))` (`vc`/Bentz) or
`α·ft/(1+√(500·ε1))` (`cm`/Collins–Mitchell), driven by the **composite `ε1`** (= the membrane
principal tensile strain the MCFT `β` already uses — perfect-bond compatibility, no separate steel
layer needed). Equibiaxial floors both in-plane normals; consistent tangent + IMPL-EX freeze.
**Monotonic-scope (v1):** re-inflates on unload (no `ε1max` memory) — use for pushover.

**3b — Crack-band (Bažant–Oh) regularization (`-autoRegularization $lch_ref`, PR #277).**
**Option A chosen** (the ADR D5 decision is now resolved): a faithful clone of `ASDConcrete3D`'s
auto-regularization — the post-peak softening is rescaled so `g_reg = G_f0·(lch_ref/lch)`, hence
the physical dissipated energy `g_reg·lch` is **mesh-objective**. `lch` is the element's scalar
`getCharacteristicLength()` (EAS-aware) or an explicit `-lch`, **latched once**; **loud failure**
if unresolved (no silent default). **Zero vanilla edit** — Option B (per-direction vanilla plumb)
was *not* taken; the single-scalar √2 inclined-strut mis-regularization is the documented residual.

> [!warning] The Phase-3 bug worth remembering (3b review)
> `regularize()` must re-apply `adjust()` (E-cap + monotone plastic strain + non-decreasing damage)
> after each scaling step — the reference does, and omitting it lets the post-peak plastic strain go
> non-monotone and corrupts `q`/damage **~6% on steep-softening backbones, invisible to every energy
> gate** (the fracture energy uses only `x,y`). Caught by a 3-agent review; gated now with a
> steep-damage plastic-strain-monotonicity check in `rc_reg_gpp.cpp`. See [[LEDGER_quirks]].

**Still open — 3b-struct (staged Zone-B gate).** The regularization *mechanism* is proven at the
material point (`g_reg·lch` constant across `lch` — oracle R1 + g++ + Zone-A). The **structural**
proof is not yet wired:

- **Goal.** In-plane **mesh-objectivity on a ROTATED (inclined-crack) / notched panel**: peak ~1–3%,
  total dissipated energy ~3–5% across a 2× refine — calibrated on the elastic shell patch first,
  NOT imported from the brick numbers.
- **What to build.** A localizing panel (a notch or a weak element row so the crack picks one band)
  at two mesh refinements, `-autoRegularization` on, driven to full softening; assert peak load +
  total dissipated energy mesh-objective. A straight **Bažant-bar** (one row of `ASDShellQ4`, weak
  first element — gmsh-free, hence Zone-A) is the cleanest first cut; the inclined notched panel is
  the harder Zone-B version.
- **Gotcha (why it's staged).** Softening localization in series **snap-backs** (the elastic
  unloading of the stiff part can outrun the band's dissipation) — displacement control may not
  converge; heavy regularization (long ductile tail) helps, else use arc-length / quasi-static
  explicit. This is the standard difficulty, not a material bug.
- **Effort:** S–M.

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
- **`lch` is a single in-plane scalar** on the director-shell seam — **resolved as Option A** in
  Phase 3b (`getCharacteristicLength()` or `-lch`, latched once, loud-fail if unresolved). The
  residual is a √2 mis-regularization of ~45° struts and no through-thickness objectivity; Phase 5's
  solid-shell needs directional/projected `lch`.
- **MCFT needs composite `ε1`** — the membrane principal tensile strain (shared by concrete + smeared
  steel under perfect-bond compatibility) the kernel already computes; tension stiffening (3a) and
  `β` both ride it. Discrete boundary bars are separate `PlateRebar` layers (recurs in V).
- **The shell path is `PlateFiber`-only** — `ASDShellQ4` consumes a *section*, not an `nDMaterial`;
  all wall shear flows through the order-5 condensed `PlateFiber`. Any new analysis is done in that
  5-component condensed setting, not a free-standing 3×3.

---

## Quick-start pointers for the next session

- **Source:** `SRC/material/nD/{LadrunoRCKernel.h, LadrunoRCConcrete.{cpp,h}}` (kernel is
  header-only, OpenSees-free, numpy-oracle-testable).
- **Oracle:** `tests/_testbed/rc_shell_ref.py` (now also `fracture_energy`/`regularize`/`tens_stiff`
  + the A1/T1/R1 gates). **Batteries:**
  `tests/test_ladrunoRCConcrete_{material,shell,implex,objectivity,wall,tensstiff,reg}.py` plus the
  CI-wired g++ gates `tests/test_ladrunoRCConcrete_{tensstiff,reg}_cpp.py` (compile+run
  `tests/_testbed/rc_{tensstiff,reg}_gpp.cpp` — the *only* tangent / fidelity gates, since a converged
  static stress is tangent-independent).
- **Build:** edit in the worktree → copy to the main checkout → **bump the file mtime** (Copy-Item
  preserves the old mtime → ninja skips the rebuild — see [[LEDGER_quirks]]) → `build.bat OpenSeesPy`
  via the PowerShell tool → confirm the `.pyd` timestamp advanced. **Test interpreter:**
  `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe` with the `dist\bin`
  `os.add_dll_directory` + `sys.path.insert` bootstrap (set `LADRUNO_OPENSEES_QUIET=1`).
- **Schema:** currently **v5** (hard-checked; v2=IMPL-EX, v3=shearRetFactor, v4=tension stiffening,
  v5=regularization). Any new committed state bumps it; keep send/recv counts balanced (the
  `RC_NSCALAR`/`RC_DATA` arithmetic) — and add a NON-default-value round-trip test (a dropped wire
  slot reverts to the ctor default and a default-valued round-trip can't catch it).
- **Build-control:** every phase updates the three ledgers (`LEDGER_implementations`,
  `LEDGER_quirks`, `LEDGER_vanilla_files` *if* a vanilla file is touched) **in the same PR**, and the
  ADR phase subsection. New fork-authored files get the LADRUNO header stamp. **CI/merge:** this fork
  auto-merges `ladruno` fast — a PR that sits goes **DIRTY**; merge `ladruno` in, and **backfill any
  manifest row a sibling PR left behind** (the G9 gate fails on the inherited gap — see
  [[feedback_stale_pr_ledger_ci]]; hit again in 3b for `ND_TAG_LadrunoCohesiveHingeBiaxial`).
- **Reviews:** the project convention is implement → adversarial review → fix → PR. Each phase here
  warrants it (the §14.11 xfail, the EAS softening stability, and the tangent degeneracy are the
  likely review hot-spots).
