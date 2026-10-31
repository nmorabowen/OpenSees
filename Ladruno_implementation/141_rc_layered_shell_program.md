---
title: "WP-141 — RC layered-shell program: keep the host vanilla, strengthen the materials"
project: Ladruno
type: ADR / work-package program
status: "PROPOSED 2026-09-27; revised 2026-10-31 after the #882 review (H1–H3). Owner decisions D3a (additive vanilla lch seam, re-scoped after P2e) and D4 (phase order) TAKEN. Design record, no code. P1a, P3c, P3d blocked on #877."
owner: nmora
related:
  - "[[19_ladruno_rc_shell_adr]]"          # LadrunoRCConcrete (33015) — the RC shell material this program strengthens
  - "[[31_ladruno_concrete3d_adr]]"        # LadrunoConcrete3D (33017) — CDPM2; P4 adds its confined plate view
  - "[[66_ladruno_solidshell_adr]]"        # LadrunoSolidShell — the punching/σ33 zoom, O4 directional-lch item
  - "[[64_ladruno_shell_to_solid_tie_adr]]" # LadrunoTie -shellSolid seam
  - "[[59_ladruno_gradient_concrete_adr]]" # gradient regularization — descoped, stays descoped here
  - "[[LadrunoRCConcrete_guide]]"
  - "[[LEDGER_quirks]]"
external:
  - "PR #877 (wp/concrete3d-integration @ f19d55bd1, OPEN, CONFLICTING vs ladruno on 2026-10-31): CDPM2 B1/B2 (the content of #867, which closes as superseded) + Phase C (#873: -crackedNu/-betaC, VC-1986 TS default c = 200) + material failures return the ADR-86b/94 sentinel LADRUNO_MATERIAL_REFUSED (−33086), which cuts the step"
  - "Validation study C:/Users/nmora/Documents/gitAPE/ladruno-concrete-validation (fix_plan.md phases A–E; family 03 = layered shells)"
tags: [adr, program, rc, shell, layered-shell, crack-band, regularization, confinement, validation]
updated: 2026-10-31
---

# WP-141 — RC layered-shell program

**Question the owner asked (2026-09-27):** is the fork ready to model reinforced concrete with
nonlinear layered shells — smeared cracking, rebar layers, element-size (crack-band) regularization —
using our materials, and are there better implementations or technologies?

**Answer:** keep the layered shell, keep the host classes vanilla, and put the work into our two
concrete materials. The architecture (fixed crack + interlock + smeared rebar layers on a MITC-type
shell) is how ATENA, VecTor and PARC_CL do it. The measurable gaps are the **regularization length**
and **transverse shear**, not the crack model.

Tasks, models, effort and oracles: [[_wp141_implementation_plan]].

This document records the assessment (§2), the decisions (§3), and four phases (§4). It was assembled
from four read-only code/doc audits of `ladruno` @ `64a0341a6`, a literature scan, PR #877, and the
external validation study. Each phase gets its own `wp/<n>-<slug>` branch when it is opened.

---

## 1. The two concrete materials are different models

A recurring confusion, stated once:

| | `LadrunoRCConcrete` (ND 33015, ADR-19) | `LadrunoConcrete3D` (ND 33017, ADR-31) |
|---|---|---|
| Spine | verbatim clone of vanilla `ASDConcrete3D` (spectral T/C split, scalar dt/dc, Lubliner envelope, backbone laws) — `LadrunoRCKernel.h` | CDPM2 (Grassl 2013): Menétrey–Willam surface, non-associated flow, confinement-dependent ductility, dual ωt/ωc — `LadrunoConcrete3DKernel.h` |
| Built for | cracked RC membranes and shells | triaxial / confined solids |
| Extra physics | MCFT β, fixed-crack interlock, cyclic slip, X-crack + wear, tension stiffening, crack band, IMPL-EX | crack band on Gf and Gc, IMPL-EX, Duvaut–Lions `-eta`, confined BeamFiber view (`-hoop`) |
| Flags off | reproduces `ASDConcrete3D` | — |

Both expose a native PlateFiber view (their own σ33 = 0 condensation). No code is shared.

## 2. Assessment — state on `ladruno` @ `64a0341a6` (+ PR #877 where noted)

### 2.1 Host (all vanilla, no `Ladruno` edits)

| Piece | Fact | Where |
|---|---|---|
| `ASDShellQ4` | AGQ6-I membrane + EAS (default on), MITC4 shear, `-corotational`, 2×2 GPs | `SRC/element/shell/ASDShellQ4.cpp` |
| lch | min node-pair distance, **halved under EAS** ("localizes in a row of Gauss points") | `ASDShellQ4.cpp:1858-1867`, `Element.cpp:714-741` |
| `LayeredShellFiberSection` | midpoint rule per layer (≥ 3 layers); `getCopy("PlateFiber")` per layer, `exit(-1)` on null | `LayeredShellFiberSection.cpp:175-191` |
| transverse shear | same γxz, γyz in every layer (uniform through thickness), **no √(5/6)** (commented out) | `LayeredShellFiberSection.cpp:433-447`, `:528-559` |
| `PlateRebarMaterial` | bar strain ε11c²+ε22s²+γ12cs, stress σ·{c²,s²,cs}, zero transverse shear; **no `setResponse` forwarding** | `PlateRebarMaterial.cpp:199-264` |
| `PlateFiberMaterial` | σ33 Newton, abs tol 1e-8, 20 its, **returns 0 when not converged** — the path `ASDConcrete3D` takes in a shell | `PlateFiberMaterial.cpp:195-258` |
| `PlateFromPlaneStress` | elastic `gmod` transverse shear (ADR-19 rejected it for RC) | `PlateFromPlaneStressMaterial.cpp:225-251` |
| explicit | `criticalTimeStep()` = −1 on `ASDShellQ4`; manual dt | `LEDGER_quirks.md:1490-1495` |
| recording | per-layer via `material.fiber.<resp>` in MPCO and `.ladruno` recorders | `LEDGER_quirks.md:735-766` |

### 2.2 Materials

- **`LadrunoRCConcrete`** — the tested shell material (`tests/test_ladrunoRCConcrete_shell.py`, `_meshobj`,
  `_tensstiff`, `_wall`, `_objectivity`). Guarded, damped σ33 Newton. It still returns a bare −1 on
  failure (legacy code, `LadrunoRCConcrete.cpp:371, 623-640` @ f19d55bd1); `ladrunoMaterialMustCut`
  accepts both −1 and the sentinel, so the step is cut either way. lch latched
  once (`-autoRegularization lch_ref`), loud failure if none. Tension stiffening is **monotonic only**
  and is **not** lch-scaled. No smeared steel inside (the ADR-19 `nWebRebar`/`rho[8]` design was never
  built). PR #877 adds `-crackedNu`, `-betaC`, and `-tensStiff vc` default `c = 200`; the study's
  PV20/PV19/PV27 pure-shear panels are within +6.4 / +5.0 / −0.9 % of test (τ_max 4.53 / 4.15 /
  6.30 MPa) with `-crackedNu 0` and `C = 0.34/ε'c` (189 for PV20, not the fixed 170 — the gate
  numbers are taken with it). Nothing gates these numbers yet; P2 G-PV does.
- **`LadrunoConcrete3D`** — PlateFiber view exists but is **untested in a shell**. On #877
  (@ f19d55bd1) a kernel return-map failure already refuses with the sentinel: the fibre returns it,
  `LayeredShellFiberSection::setTrialSectionDeformation` sums the fibre codes, `ASDShellQ4` adds the
  section code and returns it, and the step is cut. What is left is the σ33 = 0 condensation itself
  (`LadrunoConcrete3D.cpp:~590-620`): an undamped 30-iteration Newton, tolerance
  1e-9·max(|σ|, 1), whose **miss warns and accepts** (the code comment defers the honest refusal to
  P1a). Study fix-plan item **B4**; WP-141 P1a takes the condensation part. `-tcTemper proj` is the
  default now (`:125`). lch re-read every call (not latched, `:~460-468`), including the Gc → ε_fc
  table, so it inherits the EAS halving (G10); silent fallback to `-lch`. Plastic dissipation not
  regularized (~30 % lch-dependent). `-hoop` is inert outside the BeamFiber view.
- **`ASDConcrete3D`** (vanilla) — only via the silent `PlateFiberMaterial` wrapper; reads the scalar
  `getCharacteristicLength()` (`ASDConcrete3DMaterial.cpp:1615`), so G10 reaches it too.
- **Steel** — any uniaxial inside `PlateRebar`. `ASDSteel1D -auto_regularization` reads
  `getCharacteristicLength()/2` (`ASDSteel1DMaterial.cpp:2150`); on `ASDShellQ4` + EAS that is
  (min node distance)/4, unrelated to the bar direction or the tie spacing. The D3a seam is a new
  virtual with the scalar as default, so it changes none of these scalar readers
  (`LadrunoConcrete3D`, `ASDConcrete3D`, `ASDSteel1D`); for them `-noeas` (R1) is the only fix. `LadrunoUniaxialJ2` and
  `LadrunoRebarBuckling` (`-lsr s/d`) are element-independent. `LadrunoJ2` as a layer is an isotropic
  **plate**, right for steel plates/liners, wrong for rebar.

### 2.3 Known gaps, ranked by accuracy per effort

| # | Gap | Size of the error | Evidence |
|---|---|---|---|
| G1 | one scalar in-plane lch for every crack direction | 45° crack in a square linear element: band √2·h, dissipation ≈ √2·Gf (**~41 % too ductile**); elongated elements follow the short side | ADR-19:236-237; RTD 1016 §2.4.1.7; DIANA default for linear elements is √(2A) |
| G2 | element-size band in reinforced layers | too brittle when h > mean crack spacing s_rm, especially in cyclic runs where TS is off | RTD 1016 §2.4.3.1 (G_F^RC = max(1, h/s_rm)·G_F); ATENA §2.2.8 |
| G3 | tension stiffening + bare steel with no crack check | overestimates post-yield tension capacity by ≈ A_c·σ_ts; TS re-inflates on unload | Vecchio–Collins 1986; `LadrunoRCConcrete_guide.md:281-287` |
| G4 | transverse shear uniform, not degrading by depth | no out-of-plane shear/punching failure in the director shell; solid-shell zoom under-predicts PG-1 ~2.5–3× | `LadrunoSolidShell_guide.md:130-143`; Hrynyk & Vecchio 2015 (VecTor4 mean 1.01) |
| G5 | silent failures | on #877 only the `LadrunoConcrete3D` σ33-condensation miss (warn and accept) and vanilla `PlateFiberMaterial` (returns 0 unconverged) remain; kernel failures refuse | §2.1, §2.2 |
| G6 | through-thickness confinement impossible (σ33 ≡ 0) | boundary-element cores lose CDPM2's confined strength/ductility | §4 P4 |
| G7 | observability | bar stress/strain/buckling state unreachable inside `PlateRebar` | `PlateRebarMaterial.cpp` (no `setResponse`) |
| G8 | doc drift | "lch = edge length" is half that under EAS (`LEDGER_quirks.md:1523`, `tests/test_ladrunoRCConcrete_meshobj.py:149` labels); ADR-66:134-136 says in-kernel smeared steel shipped (it did not); the midpoint-rule quirk header (`LEDGER_quirks.md:2616`) says "≈2% at 5 uniform layers" — for n equal layers of one material the bending loss is exactly 1/n² (4 % at 5, 1 % at 10; per layer ∫z²dz = t·z_i² + t³/12), the measured 2.07 % was the G7 mixed stack | §2.2 |
| **G10** | **EAS halves lch, but the band does not localize into one Gauss-point row** | `ASDShellQ4` returns (min node distance)/2 under EAS (default), assuming a band half an element wide. In the G-A periodic cell and the Bažant bar all 4 GPs of the band element soften together, so crack-band energy is **doubled**: W/(Gf·t·ℓ_c) = 1.94–1.97 with EAS vs 0.99–1.00 with `-noeas` (θ = 0, n = 1/4/8). `test_ladrunoRCConcrete_meshobj.py` gates only mesh *ratios* (1.000) and measures 1.94 absolute. Affects vanilla `ASDConcrete3D` users of `ASDShellQ4` too. **P2e verdict:** single-row localization is a bifurcation. It happens only when the band's bounding node lines can stretch differently (free edges): there h/2 is right, in constrained kinematics h is right. So no fixed factor works under EAS; `-noeas` + lch = h is deterministic (0.99–1.04). Details are in the plan's P2e outcome. | G-A report 2026-09-28 (study try `work/claude/wp141-ga-inclined-crack` @ 24fca01); P2e (`work/claude/wp141-p2e-eas-drilling` @ 352dc4d) |
| **G11** | **elastic drilling springs bridge inclined cracks** | 45°, n = 8, `-noeas`: drilling springs store 0.84·Gf·t·ℓ_c at the end of loading and hold 3.6 % stress at 4·U_max (the crack never fully releases). `-drillingNL` cuts both ~10× but made damage diffuse at θ = 0, n = 4 (ratio 3.25). Relevant for diagonal cracking in squat walls. | same |
| **G12** | **elongated elements: lch is the short side, the band follows the crack normal** | min node-pair distance ignores direction. P2e notched beam, 2:1 elements (h_x × h_y = 20 × 10 mm), crack normal along x: the band is h_x (or h_x/2) wide but lch = h_y, so the local ratio is 2.06 under EAS and 1.58 under `-noeas` (rows still opening at the stop; the length error alone predicts 2.0). `-noeas` does not fix it; this, not the 45° √2, is the measured job of the D3a seam (P3a). | P2e finding 3 |
| G9 | in-plane frame of rebar layers | without `-local`, stock openseespy 3.7.1.x and the fork's build give `ASDShellQ4` in-plane frames 90° apart, so a `PlateRebar` at angle 0 can run the wrong way (apeGmsh live test, 2026-09-27) | §4 P1c |

## 3. Decisions

- **D1 — Keep the layered shell as the RC workhorse** (walls, slabs, cores: flexure, membrane,
  in-plane shear). Out-of-plane shear and punching are outside its envelope: post-check (RTD 1016
  §2.5.1) and, where it governs, a `LadrunoSolidShell` patch through `LadrunoTie -shellSolid`.
  *Rejected:* a VecTor4-style 9-node degenerate layered element (new element + validation for one
  blind spot); solid-only RC.
- **D2 — Two concrete materials with roles, not one.** `LadrunoRCConcrete` = cracked membrane/slab
  layers; `LadrunoConcrete3D` = confined crushing cores. `LayeredShell` assigns materials per layer.
  `ASDConcrete3D` is retired as a shell layer (the RC spine with flags off reproduces it, through a
  guarded view). No in-kernel smeared steel: steel stays in `PlateRebar` layers; the concrete only
  learns, as constants, the direction, ratio and yield stress of the bars that stiffen it (P3d crack
  check). It holds no steel state. *Deferred:* moving the RC shell physics onto the CDPM2 spine to
  get one concrete — large port, and CDPM2's unconfined plane-stress post-peak is snap-back prone.
- **D3 — Host stays vanilla, with two strictly additive vanilla edits** (each `// Ladruno`, each a
  `LEDGER_vanilla_files` row, default behaviour bit-identical):
  - **D3a (owner decision 2026-09-27, re-scoped 2026-10-31 after P2e):** a directional
    characteristic-length seam — a new virtual on `Element` whose default returns today's scalar,
    overridden by `ASDShellQ4` with the **full** directional length h(n), with no EAS factor. This
    **reverses ADR-19's "Option A" (zero vanilla edit)**, which accepted the √2 residual. P2e moved
    its job: the EAS factor (G10) is solved by `-noeas` with no seam, and on structured meshes the
    inclined-crack error is dominated by locking, which h(n) cannot fix. What the seam still fixes is
    the elongated-element error (G12, 1.58–2.06 on 2:1 elements), and it is correct only paired with
    `-noeas`: under EAS the band is h or h/2 by a bifurcation the element cannot see, so the seam
    reports EAS and the regularizing material warns or refuses (§6 Q2). Decided, not yet built;
    sequenced after P3c/P3d (§4).
  - **D3b:** `PlateRebarMaterial::setResponse/getResponse` forward unknown requests to the wrapped
    uniaxial material.
  A fork section class (parabolic transverse shear, Gauss/Lobatto through thickness) is **not** built
  until a punching or one-way-shear consumer is named.
- **D4 — Order (owner decision 2026-09-27):** P1 cleanup → P2 validation baseline → P3
  `LadrunoRCConcrete` regularization → P4 `LadrunoConcrete3D` confined plate view.
- **D5 — Not doing:** gradient / nonlocal / phase-field (ADR-59 descoped it; the literature finds no
  structural RC-shell result that beats a well-done crack band); sequentially linear analysis (not
  established for shells or cyclic loading); a plane-stress-projected return map (condensation is
  ≈ 2 % of wall time on a J2 plate shell; Newton iteration count is the cost).

## 4. Phases

### P1 — Cleanup (small; coordinate with #877 and study B4)

**Prerequisites for P1a (and P4):** #877 merged, and its two open kernel decisions settled: the
`-subIncr` default (`deterministic` today; `refuse` landed as an option in 0ba6a269a) and the
`-tangent` policy (chain rule vs last piece across the deterministic ladder). Ownership: #877 owns
`LadrunoConcrete3D.cpp` until it merges; P1a then branches from `ladruno` and owns the condensation
block of `setTrialStrain` only. The old agreement "the study will not touch `LadrunoConcrete3D.cpp`
until P1a" was overtaken: the study's `-subIncr refuse` work landed on #877. Merge order: #877 → P1a.

| Item | Gate |
|---|---|
| P1a `LadrunoConcrete3D` σ33 condensation (the only silent path left on #877) — the PlateFiber part of study fix_plan B4, **handed to WP-141 by the study session (2026-09-27)**; branch after #877. **Decided:** (1) a miss returns `LADRUNO_MATERIAL_REFUSED`, never −1 (ADR-94; `test_R5_*` expects the sentinel only); (2) tolerance relative to `ft`; (3) a safeguarded scalar solve, not only damping: σ33(ε33) is piecewise discontinuous at the deterministic-ladder boundaries (study register 31 M1: jumps up to 1.7 % σ, 7 % κp) and can change sign without a root, so Newton runs inside a bracket with a bisection/secant fallback, and "sign change across a jump, no root within tol" is a refusal; (4) an inner iterate whose `integrate()` failed is not used — today the next ε33 step takes that call's elastic-fallback stress and tangent, and `lastStatus` sees only the last call; a failed iterate shrinks the bracket toward the last good one; (5) dead and crushed points (ω ≥ `omegaDead`) are a defined case of the pivot (σ33 flat in ε33), solved by the bracket, not divided by a zero pivot. **Blocked on #877, and what decides it:** the condensed-tangent gate. Under a chain-rule `-tangent`, FD 1e-6 stands; under last-piece (FD error 6–121 %) it cannot pass, and the gate becomes the global Newton rate on PV20 (the condensed tangent feeds the Newton that stalled in study B1). Also `-subIncr refuse` as default would move inner failures into the bracket (20 absorbed in one PV20 run). | (a) PV20 with `LadrunoConcrete3D` layers: 0 accepted misses (every miss is a cut step, counted by `returnFailures`); (b) τ(γ) within 1 % of the #877 baseline wherever that baseline converged (4.67 MPa still rising at γ = 0.02, +9.6 % vs test: P1a is a robustness change, and accuracy against test is not its gate — CDPM2 has no MCFT compression softening); (c) the condensed `LadrunoSolidShell` softening tests (`test_ladrunoSolidShell_explicit.py::test_softening_band_sms_completes`, `test_ladrunoSolidShell_softening.py::test_g5_*`), which became hard failures when a first cut refused the miss, pass by converging; (d) a `LadrunoConcrete3D`-in-`ASDShellQ4` Zone-A test (there is none today) |
| P1b D3b `PlateRebar` response forwarding (vanilla, additive) | `eleResponse(e,'material',gp,'fiber',k,'stress')` still returns the 5-comp plate stress; a new key reaches the bar (`LadrunoRebarBuckling` state, `ASDSteel1D` damage) |
| P1c Guidance: in shells, `ASDSteel1D` without `-auto_regularization`, `-buckling` with the **tie spacing**; always pin `ASDShellQ4 -local` when a section has `PlateRebar` layers (G9 — confirm the frame difference in the fork and upstream sources first). Separate from study register 16: apeGmsh used to emit `-localCS`, which the parser ignores (fixed on apeGmsh main); the PV panels have no pin and no `-noeas` and match only because the fork's default frame happens to be right — G-PV reruns them pinned | RC guide + `LEDGER_quirks` rows |
| P1d Doc drift G8 — after #877 merges (it edits ADR-19 and the RC guide) | grep clean |
| P1f **R1** (plan §3, P2e): quirk rows + guide — crack-band RC on `ASDShellQ4` ⇒ `-noeas` with lch = h, `-drillingNL` where localized inclined cracks form; `test_ladrunoRCConcrete_meshobj.py` gates the **absolute** W/(Gf·A) under `-noeas` (today ratios only, `:155`, and its u_y = 0 chain reads 1.94 under EAS); apeGmsh warns when a crack-band layer meets EAS. Owns study defect 21 (G10/G11). Becomes guidance only after **R2** | R2: cantilever wall, n = 1/2/4/8 across, EAS vs `-noeas` vs Timoshenko (apeGmsh: 0.74× tip displacement for one element), plus a reinforced strip with bars **parallel** to the crack under EAS and `-noeas` (P2e did not test it; steel may suppress the one-column bifurcation) |
| P1e apeGmsh: `PlateFiber`, `ShellLayer` guards (uniaxial → points to `PlateRebar`), `RCLayeredShell`/`RebarMesh` builder, `ASDShellQ4(no_eas=)` | built 2026-09-27 on apeGmsh branch `feat/layered-shell-rc-primitives` (uncommitted), live 9/9 on the fork and on stock openseespy 3.7.1.2; `PlateRebar`/`PlateFromPlaneStress` (#1182) and `-local`/drilling flags (#1183) were already merged |

### P2 — Validation baseline (reuse study family 03; add only what measures P3)

The study already runs PV19/PV20/PV27; #877 has no gate on them (its only 5 % is the `_meshobj`
ratio). P2 gates are run **before** P3 so every P3 change has a before/after number, and each one
must be able to fail: a gate the baseline already passes is labelled a no-regress check, not a
target. All P2 runs use the #877 build, `LadrunoLoadControl`, `ASDShellQ4 -noeas` with lch = h (P2e)
and `-local` pinned (P1c), so P2 baseline runs wait for #877.

- **G-PV no-regress (P3c, P3d, and P3b if it is ever built):** PV19/PV20/PV27, the study's MCFT flag
  set with `-crackedNu 0` and C = 0.34/ε'c, `-local` pinned, EAS and `-noeas`. With `-noeas` lch
  becomes h, so the backbone authoring length (`lch_ref = h/2` today, a no-op under EAS) and the
  study's l_max check (2071–2545 mm vs h = 890) are redone first. Gate: τ_max within ±7 % of test on every panel. Baseline
  +5.0 / +6.4 / −0.9 % (EAS) and PV20 4×4 `-noeas` 4.276 MPa (+0.4 %, P2e B1): passes today. It
  exists because G-C is uniaxial: a rule that is right in a tie can still remove tension
  stiffening in a biaxial panel before the peak. `-betaC` alone gives 3.61 MPa on PV20 (−15 %), so
  any change that zeroes tension stiffening before the peak fails it.
- **G-A inclined crack (plan O2c):** rig = the plan's periodic unit cell (the rotated-PV20 rig of the
  2026-09-27 draft is dropped); oracle W_diss/(t·ℓ_c) = Gf with ℓ_c = L·√(p²+q²) for tanθ = p/q;
  **no tolerance — recorded, not gated**, a locking detector. The partial baseline did not show the
  expected √2: 45° `-noeas` gave 2.13 (drilling 0.84, a 3-element staircase) and tanθ = ½ locked to
  (1,0) (0.45). Neither h(n) nor `-noeas` brings 45° near 1: that needs drilling bridging (P3f) and
  locking fixed. The directional-length gate is the elongated-element strip of P3a; P3b's revival
  test is a seeded strip at θ (plan D-P2).
- **G-B size vs crack spacing (P3c's gate):** reinforced membrane tie of length L = 12·s_rm, mesh
  h/s_rm ∈ {½, 1, 2, 4}, **elastic** steel (no yield bifurcation), ρ ≈ 1 %. Two arms:
  - *TS off — energy (T1).* With perfect bond every element cracks once, so plain crack band
    dissipates n·G_F·A_c = (L/h)·G_F·A_c, while the tie forms L/s_rm cracks. Run until every
    concrete Gauss point is fully softened and gate `W_conc / (G_F·A_c·L/s_rm) = 1.00 ± 0.05` at
    h/s_rm ∈ {1, 2, 4}. Expected baseline: 1.00 / 0.50 / 0.25 (= s_rm/h), so it fails at 2 and 4
    and P3c must fix it. h/s_rm = ½ is recorded, not gated: expected 2.0 at baseline **and** after
    P3c, because min(h, s_rm) leaves it alone (§6 Q6).
  - *TS on — EC2 eq. 7.9 mean strain*, ±15 % in stabilized cracking, mesh spread ≤ 10 pp. The
    2026-09-28 run already passes (+7.0–7.1 %, 0.1 pp), so this arm is no-regress for P3c and P3d.
  - Why the 2026-09-28 G-B run (`work/claude/wp141-gb-gc-tension`) cannot gate P3c: its TS-off arm
    was compared against EC2 7.9, which contains tension stiffening, so +55–62 % is the expected
    answer, not an error; no TS-off run finished (h/s_rm = 1 stopped at ε 1.44e-3, 2 at 2.35e-3,
    4 in a solver-artefact region); the "7.0 pp spread" is two partial meshes; it ran on the
    pre-#877 June build (stock `LoadControl`, c = 500, lch_ref = h/2 under EAS, P2e's 2× case);
    and its spacing was EC2 7.11's s_r,max = 646 mm, not a mean spacing (§6 Q5).
- **G-C post-yield capacity (P3d's gate):** one element, ρ = 0.3 / 0.5 %, elastic-perfectly-plastic
  steel, uniaxial tension along the bars. TS off is the control (1.0024 / 1.0014: passes, tests
  nothing). TS on today gives `A_s·f_y + A_c·σ_ts`: +117 / +70 % at ε_y on the pre-#877 build
  (c = 500); re-measure on #877 (c = 200: smaller, still far outside). Gate after P3d: `N/(A_s·f_y)
  ≤ 1.01` at every strain after cracking and `1.00 ± 0.01` on [ε_y, 10·ε_y].

### P3 — `LadrunoRCConcrete` regularization (the accuracy phase)

Order after P2e (plan §3, re-ranking confirmed by the owner 2026-10-31): **P3d → P3c → P3a**; P3b parked;
P3e optional. **P3c and P3d are blocked on #877**: they edit the tension-stiffening and
regularization paths #877 changes (VC-1986 TS default, c = 200, `-crackedNu`), and their G-B/G-C/G-PV
baselines must be re-measured on its build. The rules below are fixed here; their gate numbers are
set on the #877 build.

- **P3d Tension-stiffening crack check** (replaces the 2026-09-27 cutoff "σ_ts → 0 once ε1 reaches
  the steel yield strain"). New flag, default off ⇒ byte-identical: one entry per bar set i —
  angle α_i (same frame as the `PlateRebar` angle, i.e. the element local frame, so P1c's `-local`
  pin applies), ρ_i = A_s,i / A_c,ef (relative to the concrete that carries tension stiffening),
  f_y,i, E_s,i. With n1 the crack normal the TS law already uses:

  `σ_ts ≤ σ_cap = Σ_i ρ_i · max(0, f_y,i − f_s,i) · cos²θ_i`,  `f_s,i = clip(E_s,i · n_iᵀ ε n_i, ±f_y,i)`,

  θ_i the angle between bar i and n1 (the MCFT crack check, Vecchio & Collins 1986; Bentz, Vecchio &
  Collins 2006). Three properties make it the right rule:
  - in a tie along the bars it is exactly "steel stress at the crack ≤ f_y", so N ≤ A_s·f_y at every
    strain and the post-yield double count (G3) is gone;
  - it is continuous: σ_cap reaches 0 only when the last bar crossing the crack yields, so no fade
    band is needed (the tangent has a kink, checked one-sided);
  - it is inactive wherever the crack steel is below yield, and an inactive cap leaves σ and the
    tangent untouched.

  The bar strain is read from the concrete layer's own strain: exact in membrane action; in bending
  it is off by the curvature times the layer-to-bar distance (documented). The MCFT shear-on-crack
  check (v_ci) stays deferred; the interlock covers it.

  *Why the cutoff on ε1 was wrong.* ε1 is a principal strain, not a bar strain. In PV20 (study
  `runs/03_layered_shell/pv_panels/PV*_rc_mcft*.csv`, `_beta`; provenance below) ε1 passes the y-bar yield strain 1.5e-3 at γ ≈ 1.7–2.0e-3, a fifth of
  γ_peak (8.1e-3 at ν = 0, 10.6e-3 at ν = 0.2). The literal rule removes tension stiffening for
  most of the pre-peak curve, which is the `-betaC`-alone response, 3.61 vs 4.25 MPa (−15 % vs
  test). A cutoff on "any bar yielded" fails the same way: the y bars yield at γ = 3.4e-3 / 4.75e-3,
  before the peak, while the x bars never yield and still carry tension across the crack. The
  crack check sees that reserve: at τ_max, σ_cap/σ_ts = 1.93 / 2.13 / 3.67 (PV20 / PV19 / PV27,
  ν = 0.2) and 1.23 / 1.35 / 2.95 (ν = 0), so the cap is inactive up to the peak on all three
  panels and their pre-peak curves are unchanged.

  *Provenance of these numbers.* The CSVs were committed in study 3a137e9 (2026-09-26) and ran on the
  **local 2026-06-25 build**, pre-#877 (pv_panels README, "Results (local build 2026-06-25)"): stock
  `LoadControl`, `-tensStiffC 200` passed explicitly, β with the hard-wired 170, **no
  `-crackedNu`**. They are the "before" row of the study's C2 table, not the #877 configuration
  G-PV gates (`-crackedNu 0`, C = 0.34/ε'c), whose CSVs are not committed. `-crackedNu 0` behaves
  like the ν = 0 runs, so the relevant PV20 margin is the thin one, ≈ 1.2, not 1.9. The 3.61 MPa
  `-betaC`-alone figure is from the same June build (ν = 0.2). All of it is re-measured in P2f on
  the #877 build before P3d's design is frozen (§6 Q7).
- **P3c Crack-spacing band.** For layers the user marks as reinforced, the tension softening uses
  `h_eff = min(h, s_θ)`, `s_θ = 1/(|cosθ|/s_rx + |sinθ|/s_ry)` (RTD 1016 §2.4.3.1), θ the crack
  normal at onset measured from x, h the element length in use (the scalar today, h(n) after P3a).
  Plain layers keep h. Tension only: crushing keeps the existing `regularize()` rescale. The onset
  rescale happens inside a Newton iteration, so the latched value is a trial variable committed only
  in `commitState` (a cut step must not keep it). s_rx, s_ry are **mean** spacings (§6 Q5). The
  apeGmsh RC builder marks layers inside h_c,ef of a rebar layer. The P3c session also checks the
  TS-on interaction (possible double count of stabilized cracking) against G-B's TS-on arm.
- **P3a Seam (D3a), after P3c/P3d and the plan's R2 (`-noeas` bending cost).**
  `virtual double Element::getCharacteristicLength(const Vector& n)` — `n` is a unit direction in the
  frame of the strains the element passes to its materials; the default ignores `n` and returns
  `getCharacteristicLength()`. `ASDShellQ4` overrides with the projection formula at the element
  centre, `h(n) = 2 / Σ_a |∇N_a · n|` (Oliver 1989; Govindjee, Kay & Simo 1995), in its local frame,
  and returns the **full** h(n) with EAS on or off: no fixed factor is right under EAS (G10). Under
  EAS the element also reports it (the interface is chosen in the P3a session; the scalar path stays
  untouched), and a material with `-autoRegularization` warns or refuses (§6 Q2). Check on a square
  Q4 of side a: `n = x` gives a, `n` at 45° gives √2·a; on a 2:1 element, `n` along the long side
  gives the long side. The consumer is the P3c onset rescale; P3b would be the second.
- **P3b Directional latch — parked.** The 45° error on structured meshes is locking and drilling
  bridging (G11), which h(n) does not fix. Revisit only when a seeded strip at θ shows a band one
  element wide along θ (plan D-P2). The design stays as written: at crack onset call the seam with
  the crack normal and rescale once, before any softening; same trial/commit rule as P3c.
  Crack onset is one definition shared by the P3b/P3c rescale, the interlock freeze and the
  `-crackedNu` switch (today two separate tests of "cracked").
- **P3e (optional) cyclic tension stiffening:** ε1,max envelope + secant unload (ADR-19 deferred item).
- **Gates.** Each gate names what fails it at baseline or under the wrong rule:
  - *P3d:* the O4a material point (σ_cap closed form at θ ∈ {0°, 30°, 60°}, one and two bar sets;
    a dropped cos² fails it); G-C (no rule fails it: +117 / +70 % pre-#877); G-PV (the ε1 or
    any-bar cutoff fails it at ≈ −15 % on PV20), with the pre-peak τ–γ curves of PV19/PV20/PV27
    bit-identical to the pre-P3d binary; G-B TS-on arm unchanged (its σ_s ≤ 0.8 f_y keeps the cap
    inactive).
  - *P3c:* G-B TS-off energy at h/s_rm ∈ {1, 2, 4} (baseline 1.00 / 0.50 / 0.25 → 1.00 ± 0.05;
    `max` instead of `min` fails it); G-B TS-on and G-PV no-regress.
  - *P3a:* seeded strip, `-noeas`, 2:1 and 1:2 elements, crack normal along each side:
    W/(Gf·t·ℓ_c) = 1.00 ± 0.05 (baseline on 2:1 with the normal along the long side ≈ 2.0; P2e beam
    1.58 partial). G-A recorded only.
  - All: flags-off byte identity on the RC suite; numpy oracle + g++ check for each new path;
    mutation gate per ADR-87.

### P4 — `LadrunoConcrete3D` confined plate view

Confinement from cross-ties through the wall thickness is impossible while σ33 ≡ 0. Make it passive,
reusing the BeamFiber machinery (`LadrunoConcrete3DKernel.h`: `hoopStress` `:~2735`,
`driveConfinedFiber` `:~2775-2900` @ f19d55bd1) with one confined direction. #877 added dead-point
and crushed-freeze handling to that path (the dead-aware effective-stress return, `:~2487`, review
M2), and the deterministic sub-increment ladder runs inside its nested lateral solve; P4 inherits
both, and has the same prerequisites as P1a (the two #877 kernel decisions):

- residual `σ33 + σ_h(ε33) = 0` (today `σ33 = 0`); hoop stiffness added to the condensation pivot
  and to `condenseTangent()`;
- `K = k_e·ρ_t·E_s`, `f_y,cap = k_e·ρ_t·f_yt` — ρ_t the volumetric ratio of the through-thickness ties,
  k_e Mander's effectiveness; single direction ⇒ rectangular ties are fine (the BeamFiber view is
  circular-only);
- **hoop memory:** a committed hoop plastic strain (today `hoopStress = min(K·ε, fy)` is path-
  independent and unloads along its loading curve — wrong for cyclic walls); serialized;
- **nominal** stress balance (the BeamFiber view balances the *effective* stress; acceptable at the
  peak where ω ≈ 0, not in a cyclically damaged boundary element) — decide and document the divergence;
- only core layers of boundary-element sections get `-hoop`; `hoopK = 0` ⇒ byte-identical to today,
  **dead and crushed points included** (plan O5b, O5e);
- **not** a constant f_l: it would pre-compress the layer at zero load and put compression through the
  thickness during the boundary element's tension excursions.

**Depends on** study B1/B2 (CDPM2 plastic potential and compressive-damage drive, in #877; #867 closes as superseded):
confined accuracy is meaningless before those land. **Gates:** material-point oracle for the confined
plate path; `hoopK = 0` bit-identity; a boundary-element prism as shell-with-confined-view vs
`LadrunoBrick` + embedded ties (study E3, Sheikh & Uzumeri); Thomsen & Wallace RW2 boundary elements
when study family 03 reaches it.

## 5. Vanilla footprint (planned)

| File | Edit | Phase |
|---|---|---|
| `SRC/element/Element.h/.cpp` | new virtual `getCharacteristicLength(const Vector&)`, default = scalar | P3a (after P3c/P3d) |
| `SRC/element/shell/ASDShellQ4.h/.cpp` | override with the full projection length h(n), no EAS factor; report EAS to the caller. The scalar `getCharacteristicLength()` (and its `/2`) is untouched | P3a (after P3c/P3d) |
| `SRC/material/nD/PlateRebarMaterial.h/.cpp` | forward `setResponse`/`getResponse` to the uniaxial | P1b |

## 6. Open questions

1. Under `-corotational`, confirm the frame of the material's crack normal equals the frame the
   override projects in (both are `ASDShellQ4`'s local system — verify, don't assume).
2. ~~Whether the EAS factor ½ should apply to the directional length (G-A decides).~~ Settled by P2e:
   no fixed factor is right under EAS, so the seam returns the full h(n). **Decided (owner,
   2026-10-31):** warn once now. A material with `-autoRegularization` on an EAS `ASDShellQ4` prints
   one warning that names the up-to-2× energy error. It refuses with the material sentinel only once
   P3a lands, and that change goes in the P3a release note. Refusing today would break every
   existing EAS crack-band model before a fix exists.
3. How the apeGmsh builder should expose "reinforced layer" marking for P3c without leaking OpenSees
   flags into geometry code.
4. Whether `ASDShellQ4` should also override `getExplicitCriticalTimeStep` (explicit cyclic walls use a
   manual dt today) — out of scope unless P2/P4 runs need it.
5. **Which spacing P3c uses — decided (owner, 2026-10-31).** RTD's rule wants the mean crack spacing; EC2 eq.
   7.11 gives s_r,max. Adopted: s_rm = s_r,max/1.7 (the characteristic-to-mean ratio behind EC2
   7.11; MC 1990 β = 1.7), so the G-B tie of the 2026-09-28 run has s_rm ≈ 380 mm, not 646. The
   G-B energy gate checks the rule against whatever s_rm is chosen; whether that s_rm is physical is
   a separate check against a tension-tie test with counted cracks.
6. **h < s_rm in reinforced layers (P3c design session).** With perfect bond (`PlateRebar` shares
   the strain) nothing stops every element of a uniformly stressed reinforced zone from cracking,
   so plain crack band over-dissipates by s_rm/h (G-B h/s_rm = ½: expected 2.0). `min(h, s_θ)`
   leaves that case alone, as RTD does. Using h_eff = s_θ for every reinforced layer would fix it
   and would change today's fine-mesh results. **Owner (2026-10-31): decide from the G-B ½ arm on
   the #877 build (P2f), not before.**
7. P3c and P3d wait for #877. If #877 changes the TS law or its default before merging, the P3d
   margins quoted in §4 (σ_cap/σ_ts at τ_max) are re-measured on the merged build before P3d's
   design is frozen.
8. **Upstream report to the `ASDShellQ4` authors** (P2e evidence, no patch): needs owner approval
   and is not tracked anywhere yet; open a fork issue when approved.
9. **P3f (research):** a drilling law that follows the section's current secant stiffness (vanilla
   `ASDShellQ4`, owner decision; an upstream suggestion is the alternative).

## 7. Sources

Code audits (file:line above). Literature (from the scan; items marked ⚠ were read only as abstracts):
Oliver 1989, IJNME 28(2); Govindjee, Kay & Simo 1995, IJNME 38 (doi:10.1002/nme.1620382105); Jirásek &
Bauer 2012, Comput. Struct. 110–111:60–78; Slobbe, Hendriks & Rots 2013, Eng. Fract. Mech. 109 ⚠;
Hendriks & Roosen (eds.) RTD 1016-1:2022 (read in full); ATENA Theory manual §2.2.7–2.2.8 (read in
full); Vecchio & Collins 1986, ACI J. 83(2); Bentz, Vecchio & Collins 2006, ACI Struct. J. 103(4)
(SMCFT; crack check) ⚠; EN 1992-1-1 §7.3.4 (eqs. 7.9, 7.11); Belarbi & Hsu 1994, ACI Struct. J. 91(4) ⚠;
Hrynyk & Vecchio 2015, J. Struct. Eng. 141(12) ⚠; Belletti, Scolari & Vecchi 2017 (PARC_CL 2.0) ⚠;
Dhakal & Maekawa 2002, J. Struct. Eng. 128(9).
