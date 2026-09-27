---
title: "WP-141 — RC layered-shell program: keep the host vanilla, strengthen the materials"
project: Ladruno
type: ADR / work-package program
status: "PROPOSED 2026-09-27. Owner decisions D3a (additive vanilla lch seam) and D4 (phase order) TAKEN. No code yet."
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
  - "PR #877 (wp/concrete3d-integration, OPEN): CDPM2 fixes #867 + Phase C (#873: -crackedNu/-betaC, VC-1986 TS default, bare -1 cuts the step)"
  - "Validation study C:/Users/nmora/Documents/gitAPE/ladruno-concrete-validation (fix_plan.md phases A–E; family 03 = layered shells)"
tags: [adr, program, rc, shell, layered-shell, crack-band, regularization, confinement, validation]
updated: 2026-09-27
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
  `_tensstiff`, `_wall`, `_objectivity`). Guarded, damped σ33 Newton, returns −1 on failure. lch latched
  once (`-autoRegularization lch_ref`), loud failure if none. Tension stiffening is **monotonic only**
  and is **not** lch-scaled. No smeared steel inside (the ADR-19 `nWebRebar`/`rho[8]` design was never
  built). PR #877 adds `-crackedNu`, `-betaC`, and `-tensStiff vc` default `c = 200`; the study's
  PV19/PV20/PV27 pure-shear panels are within +6.4 / +5.0 / −0.9 % of test with `C = 0.34/ε'c`.
- **`LadrunoConcrete3D`** — PlateFiber view exists but is **untested in a shell** and its σ33 Newton is
  undamped with a fixed tolerance and **returns 0 on non-convergence** (`LadrunoConcrete3D.cpp:440-459`;
  unchanged on the #877 branch). Study fix-plan item **B4**; WP-141 P1a takes the PlateFiber part. lch re-read every call
  (not latched), silent fallback to `-lch`. Plastic dissipation not regularized (~30 % lch-dependent).
  `-hoop` is inert outside the BeamFiber view.
- **`ASDConcrete3D`** (vanilla) — only via the silent `PlateFiberMaterial` wrapper.
- **Steel** — any uniaxial inside `PlateRebar`. `ASDSteel1D -auto_regularization` reads
  `getCharacteristicLength()/2` (`ASDSteel1DMaterial.cpp:2150`); on `ASDShellQ4` + EAS that is
  (min node distance)/4, unrelated to the bar direction or the tie spacing. `LadrunoUniaxialJ2` and
  `LadrunoRebarBuckling` (`-lsr s/d`) are element-independent. `LadrunoJ2` as a layer is an isotropic
  **plate**, right for steel plates/liners, wrong for rebar.

### 2.3 Known gaps, ranked by accuracy per effort

| # | Gap | Size of the error | Evidence |
|---|---|---|---|
| G1 | one scalar in-plane lch for every crack direction | 45° crack in a square linear element: band √2·h, dissipation ≈ √2·Gf (**~41 % too ductile**); elongated elements follow the short side | ADR-19:236-237; RTD 1016 §2.4.1.7; DIANA default for linear elements is √(2A) |
| G2 | element-size band in reinforced layers | too brittle when h > mean crack spacing s_rm, especially in cyclic runs where TS is off | RTD 1016 §2.4.3.1 (G_F^RC = max(1, h/s_rm)·G_F); ATENA §2.2.8 |
| G3 | tension stiffening + bare steel with no crack check | overestimates post-yield tension capacity by ≈ A_c·σ_ts; TS re-inflates on unload | Vecchio–Collins 1986; `LadrunoRCConcrete_guide.md:281-287` |
| G4 | transverse shear uniform, not degrading by depth | no out-of-plane shear/punching failure in the director shell; solid-shell zoom under-predicts PG-1 ~2.5–3× | `LadrunoSolidShell_guide.md:130-143`; Hrynyk & Vecchio 2015 (VecTor4 mean 1.01) |
| G5 | silent failures | `LadrunoConcrete3D` plate view and vanilla `PlateFiberMaterial` return 0 unconverged | §2.1, §2.2 |
| G6 | through-thickness confinement impossible (σ33 ≡ 0) | boundary-element cores lose CDPM2's confined strength/ductility | §4 P4 |
| G7 | observability | bar stress/strain/buckling state unreachable inside `PlateRebar` | `PlateRebarMaterial.cpp` (no `setResponse`) |
| G8 | doc drift | "lch = edge length" is half that under EAS (`LEDGER_quirks.md:1523`, `tests/test_ladrunoRCConcrete_meshobj.py:149` labels); ADR-66:134-136 says in-kernel smeared steel shipped (it did not); the midpoint-rule quirk header (`LEDGER_quirks.md:2616`) says "≈2% at 5 uniform layers" — for n equal layers of one material the bending loss is exactly 1/n² (4 % at 5, 1 % at 10; per layer ∫z²dz = t·z_i² + t³/12), the measured 2.07 % was the G7 mixed stack | §2.2 |
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
  learns the steel yield strain (P3). *Deferred:* moving the RC shell physics onto the CDPM2 spine to
  get one concrete — large port, and CDPM2's unconfined plane-stress post-peak is snap-back prone.
- **D3 — Host stays vanilla, with two strictly additive vanilla edits** (each `// Ladruno`, each a
  `LEDGER_vanilla_files` row, default behaviour bit-identical):
  - **D3a (owner decision 2026-09-27):** a directional characteristic-length seam — a new virtual on
    `Element` whose default returns today's scalar, overridden by `ASDShellQ4`. This **reverses
    ADR-19's "Option A" (zero vanilla edit)**, which accepted the √2 residual.
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

| Item | Gate |
|---|---|
| P1a `LadrunoConcrete3D` plate view: tolerance relative to `ft`, damped Newton, return −1 on failure — the PlateFiber part of study fix_plan B4, **handed to WP-141 by the study session (2026-09-27)**; branch after #877 | PV20 with `LadrunoConcrete3D` layers: 0 silent failures (every failure is a cut step); a `LadrunoConcrete3D`-in-`ASDShellQ4` Zone-A test (there is none today) |
| P1b D3b `PlateRebar` response forwarding (vanilla, additive) | `eleResponse(e,'material',gp,'fiber',k,'stress')` still returns the 5-comp plate stress; a new key reaches the bar (`LadrunoRebarBuckling` state, `ASDSteel1D` damage) |
| P1c Guidance: in shells, `ASDSteel1D` without `-auto_regularization`, `-buckling` with the **tie spacing**; always pin `ASDShellQ4 -local` when a section has `PlateRebar` layers (G9 — confirm the frame difference in the fork and upstream sources first) | RC guide + `LEDGER_quirks` rows |
| P1d Doc drift G8 — after #877 merges (it edits ADR-19 and the RC guide) | grep clean |
| P1e apeGmsh: `PlateFiber`, `ShellLayer` guards (uniaxial → points to `PlateRebar`), `RCLayeredShell`/`RebarMesh` builder, `ASDShellQ4(no_eas=)` | built 2026-09-27 on apeGmsh branch `feat/layered-shell-rc-primitives` (uncommitted), live 9/9 on the fork and on stock openseespy 3.7.1.2; `PlateRebar`/`PlateFromPlaneStress` (#1182) and `-local`/drilling flags (#1183) were already merged |

### P2 — Validation baseline (reuse study family 03; add only what measures P3)

The study already runs PV19/PV20/PV27 (±5 % gate in #877). P2 adds three gates that isolate G1–G3,
run **before** P3 so every P3 change has a before/after number:

- **G-A mesh bias:** PV20 on a 4×4 mesh rotated 0°, 22.5°, 45° to the load axes. Metric: spread of
  peak τ and of dissipated energy. Expected today: energy grows toward ≈ √2 at 45°.
- **G-B size vs crack spacing:** PV20 at 1×1, 2×2, 4×4, plus a reinforced tension tie with a known
  crack spacing (element size above and below s_rm). Metric: post-cracking stiffness and energy.
- **G-C post-yield capacity:** a lightly reinforced panel in uniaxial tension with `-tensStiff`.
  Metric: yield plateau vs A_s·f_y.

### P3 — `LadrunoRCConcrete` regularization (the accuracy phase)

- **P3a Seam (D3a).** `virtual double Element::getCharacteristicLength(const Vector& n)` — `n` is a
  unit direction in the frame of the strains the element passes to its materials; the default ignores
  `n` and returns `getCharacteristicLength()`. `ASDShellQ4` overrides with the projection formula at the
  element centre, `h(n) = 2 / Σ_a |∇N_a · n|` (Oliver 1989; Govindjee, Kay & Simo 1995), in its local
  frame, keeping the upstream EAS factor ½. Check on a square Q4 of side a: `n = x` gives a,
  `n` at 45° gives √2·a.
- **P3b Directional latch.** New flag (default off ⇒ byte-identical): at crack onset — the event that
  already freezes the interlock crack normal — call the seam with the crack normal and rescale the
  softening once, before any softening has happened. Freezing the normal must no longer require
  `-interlock`. Second crack (`-xcrack`): v1 regularizes with the first crack only (documented).
- **P3c Crack-spacing band.** For layers the user marks as reinforced: `h_eff = min(h(n), s_rm,θ)`,
  with `s_θ = 1/(cosθ/s_x + sinθ/s_y)` (RTD 1016 §2.4.3.1). Plain layers keep `h(n)`. The apeGmsh RC
  builder marks layers inside h_c,ef of a rebar layer.
- **P3d Tension-stiffening cutoff.** `-tensStiffEy ε_y`: σ_ts → 0 once ε1 reaches the steel yield
  strain (RTD's cap in h_c,ef; removes the post-yield double count without coupling the concrete to a
  steel model). Full MCFT crack check stays deferred.
- **P3e (optional) cyclic tension stiffening:** ε1,max envelope + secant unload (ADR-19 deferred item).
- **Gates:** P2 G-A/B/C re-run (target: G-A energy spread ≤ 10 %); flags-off byte-identity on the RC
  suite; numpy oracle + g++ check for each new path; mutation gate per ADR-87.

### P4 — `LadrunoConcrete3D` confined plate view

Confinement from cross-ties through the wall thickness is impossible while σ33 ≡ 0. Make it passive,
reusing the BeamFiber machinery (`LadrunoConcrete3DKernel.h:1585-1700`) with one confined direction:

- residual `σ33 + σ_h(ε33) = 0` (today `σ33 = 0`); hoop stiffness added to the condensation pivot
  and to `condenseTangent()`;
- `K = k_e·ρ_t·E_s`, `f_y,cap = k_e·ρ_t·f_yt` — ρ_t the volumetric ratio of the through-thickness ties,
  k_e Mander's effectiveness; single direction ⇒ rectangular ties are fine (the BeamFiber view is
  circular-only);
- **hoop memory:** a committed hoop plastic strain (today `hoopStress = min(K·ε, fy)` is path-
  independent and unloads along its loading curve — wrong for cyclic walls); serialized;
- **nominal** stress balance (the BeamFiber view balances the *effective* stress; acceptable at the
  peak where ω ≈ 0, not in a cyclically damaged boundary element) — decide and document the divergence;
- only core layers of boundary-element sections get `-hoop`; `hoopK = 0` ⇒ byte-identical to today;
- **not** a constant f_l: it would pre-compress the layer at zero load and put compression through the
  thickness during the boundary element's tension excursions.

**Depends on** study B1/B2 (CDPM2 plastic potential and compressive-damage drive, in #867/#877):
confined accuracy is meaningless before those land. **Gates:** material-point oracle for the confined
plate path; `hoopK = 0` bit-identity; a boundary-element prism as shell-with-confined-view vs
`LadrunoBrick` + embedded ties (study E3, Sheikh & Uzumeri); Thomsen & Wallace RW2 boundary elements
when study family 03 reaches it.

## 5. Vanilla footprint (planned)

| File | Edit | Phase |
|---|---|---|
| `SRC/element/Element.h/.cpp` | new virtual `getCharacteristicLength(const Vector&)`, default = scalar | P3a |
| `SRC/element/shell/ASDShellQ4.h/.cpp` | override with the projection formula | P3a |
| `SRC/material/nD/PlateRebarMaterial.h/.cpp` | forward `setResponse`/`getResponse` to the uniaxial | P1b |

## 6. Open questions

1. Under `-corotational`, confirm the frame of the material's crack normal equals the frame the
   override projects in (both are `ASDShellQ4`'s local system — verify, don't assume).
2. Whether the EAS factor ½ should apply to the directional length unchanged (G-A decides).
3. How the apeGmsh builder should expose "reinforced layer" marking for P3c without leaking OpenSees
   flags into geometry code.
4. Whether `ASDShellQ4` should also override `getExplicitCriticalTimeStep` (explicit cyclic walls use a
   manual dt today) — out of scope unless P2/P4 runs need it.

## 7. Sources

Code audits (file:line above). Literature (from the scan; items marked ⚠ were read only as abstracts):
Oliver 1989, IJNME 28(2); Govindjee, Kay & Simo 1995, IJNME 38 (doi:10.1002/nme.1620382105); Jirásek &
Bauer 2012, Comput. Struct. 110–111:60–78; Slobbe, Hendriks & Rots 2013, Eng. Fract. Mech. 109 ⚠;
Hendriks & Roosen (eds.) RTD 1016-1:2022 (read in full); ATENA Theory manual §2.2.7–2.2.8 (read in
full); Vecchio & Collins 1986, ACI J. 83(2); Belarbi & Hsu 1994, ACI Struct. J. 91(4) ⚠;
Hrynyk & Vecchio 2015, J. Struct. Eng. 141(12) ⚠; Belletti, Scolari & Vecchi 2017 (PARC_CL 2.0) ⚠;
Dhakal & Maekawa 2002, J. Struct. Eng. 128(9).
