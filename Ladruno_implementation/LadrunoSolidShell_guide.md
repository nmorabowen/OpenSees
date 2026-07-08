# LadrunoSolidShell — user & validation guide (ADR 66, `ELE_TAG` 33020)

`LadrunoSolidShell` is an 8-node solid-shell: a hexahedral element that carries a
**full 3D constitutive state at its integration points** (so `σ33`, transverse shear,
and a real through-thickness stress profile are all present) while remaining locking-free
in the thin-plate/shell limit through **assumed natural strain (ANS)** shear+thickness
tying and a **state-dependent enhanced-assumed-strain (EAS)** ζ-mode. It is the ADR-19
"Phase 5" host: the element that can represent behaviours a director (Reissner–Mindlin)
shell structurally cannot — most importantly, **punching shear**.

This guide publishes the G1–G10 validation record and the modelling recipes. It is the
companion to the design record in `66_ladruno_solidshell_adr.md`.

---

## 1. Element command

```
element LadrunoSolidShell $tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $matTag \
        <-nz $nz> <-quadz gauss|lobatto> <-formulation ans|std> \
        <-geom linear> <-lumped>
```

- `$n1..$n8` — hex connectivity; nodes 1–4 are the bottom face traversed so the
  bottom-face normal points to the top face (right-hand rule). Curved/radial stacks
  must keep `detJ > 0` (see G3 connectivity note below).
- `$matTag` — any 3D `nDMaterial` (elastic, `LadrunoConcrete3D`, `LadrunoRCConcrete`,
  `LadrunoJ2`, …). The element is **transparent** to it (G5 hosting-parity gate).
- `-nz` — through-thickness integration points (default 2). `gauss` 2–5 or `lobatto`
  3–7. **`lobatto` n=2 is rejected** (trapezoid rule cannot integrate the ζ² bending
  term). Use `lobatto` when you need a sampling point exactly on each face.
- `-formulation` — `ans` (default: ANS shear+thickness + EAS) or `std` (plain
  displacement, A/B locking control). `std` **locks** in bending/shear and is for
  locking-demonstration only.
- `-geom linear` — v1 assembles directly in the global frame (the `SolidTransformation`
  corot/finite seam is P5.4 backlog; the parser rejects `corot`/`finite` and
  `FiniteStrainNDMaterial`).
- `-lumped` — row-sum/HRZ diagonal mass for explicit dynamics (see G9). Default is
  consistent mass.
- Density rides the material (`getRho`); there is no `-rho` flag. The crack-band
  characteristic length supplied to regularizing materials is `√(midsurface area)`
  (verified against `LadrunoConcrete3D`'s `-autoRegularization` latch to 1e-9).

Responses: `stresses`, `strains`, `alpha`/`eas`, `stiff`/`stiffInitial`,
`material $gp $arg` (per-GP; `$gp` in `1..4·nz`).

---

## 2. Formulation (one paragraph)

Covariant small-strain `E_ij = sym(g_i · u_,j)` mapped per-GP to a 6×6 engineering-strain
operator `T` from the contravariant basis, feeding any 3D material. **ANS**: Dvorkin–Bathe
transverse shear (`γ13` tied at `(0,∓1,ζ)`, `γ23` at `(∓1,0,ζ)`) + Betsch–Stein bilinear
`E33` from the mid-surface corners kills transverse-shear and curvature-thickness locking.
**EAS**: one `ζ·α₁` enhanced mode, **state-dependent** (live-σ inner Newton, scalar `Kaa`
condensation, committed `alphaCommit`), relaxes the Poisson thickness lock. 2×2 in-plane ×
`-nz` through thickness.

---

## 3. Validation record (G1–G10)

| Gate | What | Result |
|---|---|---|
| **G1** | membrane / bending / **thickness** patch tests | exact to solver tol; the σ33-through-a-distorted-mesh thickness patch catches a broken `∫M dV = 0` |
| **G2** | slenderness locking sweep `t/L ∈ {1e-1,1e-2,1e-3}`, 1 element thick | center deflection flat (spread < 0.05); `-formulation std` locks 8–11× |
| **G3** | Scordelis-Lo + pinched cylinder (`test_ladrunoSolidShell_benchmarks.py`) | SL 0.948 / 0.974 / **0.989** / 0.992 (4/8/16/24); PC 0.378 / 0.745 / **0.927** / 0.970 — **16×16 == canonical MITC4 ≈ 0.93** |
| **G4** | curved / trapezoidal-mesh bending (Betsch-Stein) | no curvature-thickness locking; 27° fibre-tilt trapezoid ans 0.854 vs std collapse |
| **G5** | crack-band energy objectivity (`test_ladrunoSolidShell_softening.py`) | **0.9 % work spread** across a 4× in-plane-square mesh range with `-autoRegularization`+`-implex`; fixed-`lch` control at the 0.25× band-volume ratio (the disease it cures) |
| **G6** | EAS α stability deep post-peak | α bounded/convergent on the steep `-lch 50` law; the P5.1 freeze-α guard suffices (fallback ladder dormant) |
| **G7** | moment-curvature vs `LayeredShell` (`test_ladrunoSolidShell_flexure.py`) | see §4 |
| **G8** | **PG-1 punching** (`test_ladrunoSolidShell_punching.py`) | see §5 — **the ADR-19 discharge** |
| **G9** | explicit path (`test_ladrunoSolidShell_explicit.py`) | `-lumped` mass = thickness CFL `dt_cr = t/c`, halves with thickness, bisection-validated (stable 0.9×, unstable 1.5×); SMS ≥ 15× fewer steps on a sliver-bound mesh |
| **G10** | 6-DOF↔3-DOF seam (`test_ladrunoSolidShell_seam.py`) | see §6 |

### 4. G7 — moment-curvature parity (D1, the flexural-claim gate)

Same RC wall section (`t=200`, `ρ_tot=0.6 %` at `z=±80`) built two ways: the director-shell
reference (`ASDShellQ4` + `LayeredShell` of `LadrunoRCConcrete` PlateFiber layers +
`PlateRebar`) vs a `LadrunoSolidShell` stack (3D `LadrunoRCConcrete` at the GPs, `Steel01`
interface trusses). Classic M–κ (dummy-node `equationConstraint` section driver, N held
exactly), implicit both arms, κ → 6e-5 (deep cracking + steel yield):

| Case | dM / Mmax | notes |
|---|---|---|
| elastic anchor | exact both | layered midpoint-rule `Σh³/12` deficit ≈ 2 % (predicted, not an error) |
| M–κ, N = 0 | **2.2 %** | ε₀-migration 1.9 %, N held to machine zero |
| M–κ, N = −0.1·fc·Ag | **1.6 %** | ε₀-migration 5.1 % |
| restrained (ε₀=0), 6-element stack | dM 2.1 %, dN 16 % of ft·t·b | |

**O2 decided — element stacking, not per-GP-layer assignment.** One core element through
`[−80,80]` leaves dN ≈ 62 % of ft·t·b (a single EAS ζ-mode + nodal-uz mean cannot relax the
kinked cracked-state σ33 profile); a 6-element stack converges it to ≈ 4 % relative.
Stacking is also where interface rebar naturally attaches.

*IMPL-EX reporting quirk:* under a fully-prescribed rig, an `-implex` material inside
`ASDShellQ4` reports **bit-identical** to the implicit run (ASDShellQ4 reports the
post-commit state; IMPL-EX re-integrates implicitly at commit), whereas trial-state hosts
(`LadrunoSolidShell`, `ShellMITC4`) show the extrapolated stresses. Parity gates therefore
run implicit-vs-implicit.

### 5. G8 — PG-1 punching headline (the ADR-19 discharge)

Specimen: **Guandalini & Muttoni PG-1** (EPFL 2004; ACI Struct. J. 106(1) 2009). 3.0×3.0 m
slab, `h=250`, `d=210`, 260 mm square column, `ρ=1.5 %` (Ø20@100, `fy=573`), `fc=27.7`,
**measured `Ec=25.7 GPa`**. Published **`V_R = 1023 kN`** (includes 73 kN self-weight + rig
inside the critical perimeter; net applied 950 kN), `ψ_R ≈ 7.4 mrad` (10.2 mm at r=1.38 m).
Quarter model, symmetry planes, edge roller support, column loaded through an elastic stub,
tension steel as `Steel01` trusses on the `z=+40` interface plane (⇒ `d=210`).

**The scientific result — the ADR-19 blind spot, discharged:**

| Model | Behaviour at ~1.55 mm column deflection |
|---|---|
| `LadrunoSolidShell` stack | **forms a through-thickness shear CONE**: the near-column concrete is fully damaged in the upper core (`ω → 1.0` at the z=125..225 layer, not just the tension face) while the tension steel is still only ~15 % of yield — shear-damage localizing through the thickness before the flexural steel engages IS punching |
| `ASDShellQ4` + `LayeredShell` (director shell) | **no shear cone, no limit**: it hardens MONOTONICALLY and is still gaining capacity at the same deflection — no `σ33`, transverse shear is a resultant with no failure surface, so it cannot form the cone; it can only harden or eventually fail in flexure |

*What is gated, and why not the tangent/drop.* The gate reads the **stable pair**: the
through-thickness cone (`ω → 1.0` in the upper core) with **steel below yield** on the
solid, versus **monotonic hardening** (a strictly non-decreasing, still-climbing load) on
the director shell. The load-deflection tangent *does* collapse toward the limit, but that
signature lives in the deep-softening region where `LadrunoConcrete3D` triggers **thousands
of material return-map step-cuts** — the run crawls (the last ~0.3 mm to the limit took
~440 s in probing) *and* the tangent is step-schedule sensitive (measured 0.09–0.22 across
step counts). The full post-peak drop is an HPC/offline exercise (every BC idealization
tried — roller, edge-clamped, rigid-punch — stalls on an EAS wall or grinds). The cone +
unyielded-steel + monotonic-shell trio is stable (ω=1.00, steel 0.14–0.15 across step
schedules and BCs) and runs in ~27 s. See `test_ladrunoSolidShell_punching.py`.

**Capacity — an honest, documented under-prediction (NOT a ±20 % match).** With a physical
roller edge support the model reaches its punching limit at **≈ 330 kN (total) vs 1023 kN**
— a ~2.5–3× conservative under-prediction. This is physical, not a bug: isotropic tension
damage degrades the shear/compression stiffness across the flexural crack, so the
diagonal-tension resistance is shed at ~half the expected shear stress (~0.9 vs ~1.8 MPa).
It is **mesh-insensitive downward** — a finer mesh gives ~329 kN (lower) and stalls on the
EAS wall near 1.65 mm — and `ft` is a weak lever (+16 % capacity for +37 % `ft`). The
support idealization matters: an over-restrained edge-clamped ring inflates the limit to
~710 kN (−23 %), so the honest number carries its BC. Genikomsou & Polak (2015) reach −16 %
only with a *calibrated* CDP on a ~20 mm mesh (~250 k-element model) — HPC/offline, out of
reach of a fork pytest. The G8 capacity gate is therefore a **documented
uncalibrated-conservative regression band**; the validated science is the punching-limit
MODE and the director-shell discharge. A calibrated quantitative capacity study is future
offline work.

### 6. G10 — 6-DOF ↔ 3-DOF seam recipe

Butt an `ASDShellQ4` panel (ndf 6) to a `LadrunoSolidShell` patch (ndf 3) with
**`LadrunoTie -shellSolid`** (Lagrange static / `LadrunoProjection` + `-lumped` + master-θ
rotary mass for explicit). A constant-moment state crosses the seam **exactly**: interface
rotation `= θ/2` to 1e-9, station moment `= θEI/L` to machine precision on both sides, GP
`σxx = M·z/I` at every solid GP (no seam boundary layer); membrane at 1e-14; explicit
`dt_cr`-neutral with 1.7e-17 tie violation over 300 CDL steps. **Rejected candidates
(failure mode quantified):** translations-only `equalDOF` rows = moment hinge (M=0);
`rigidLink 'beam'` across ndf = a **silent-disconnect trap** (warns, adds nothing, solves
disconnected); `-hermite` = shell-edge↔shell-edge vocabulary, refused.

---

## 7. Modelling recipes

- **Through thickness:** stack 2–3+ elements (O2). For an RC section put a stack interface
  on each reinforcement plane and attach `Steel01`/`Steel02` trusses there
  (area = `As_per_width × tributary`, both directions). Cover elements can be thin
  `-nz 2`; resolve the core with more points if you need the σ33 profile.
- **Softening concrete:** `LadrunoConcrete3D -autoRegularization -implex` + constant-dλ
  displacement control. IMPL-EX gives an SPD-ish secant so plain Newton (or KrylovNewton)
  tracks localization without step-cutting; keep in-plane elements ~square so
  `lch = √A` matches the band width (directional-`lch` is O4 backlog).
- **Large 3D models:** use a sparse unsymmetric solver (`system UmfPack`) — the damaged
  tangent is non-symmetric and `BandGeneral` is slow on a 3D block.
- **Explicit dynamics:** `-lumped` + `system Diagonal` + `CentralDifferenceLadruno`
  (`-lump hrz` pencil) or `CentralDifferenceSMS $dtTarget -lump hrz`. `dt_cr` is
  thickness-bound (`t/c`). Under SMS drive support motion with a **smoothstep** protocol,
  never a Linear ramp (a velocity step × √s mass-scaling shock cracks the loaded face).

## 8. Known limitations

- `-formulation std` locks (bending/shear/Poisson) — demonstration only.
- **Punching capacity is under-predicted** by the uncalibrated coarse model (§5); the mode
  is right, the load number is conservative. A quantitative capacity match needs a
  calibrated CDP at fine (~20 mm) mesh.
- The state-dependent EAS inner Newton makes the per-step cost high; explicit
  quasi-static runs of large concrete models are expensive, and very fine meshes can hit
  an EAS convergence wall in the near-column stress-concentration region.
- v1 is `-geom linear` only (corot/finite = P5.4).

## 9. Files

`SRC/element/ladrunoSolidShell/` (element) · `tests/test_ladrunoSolidShell_*.py`
(element, benchmarks, softening, flexure, seam, explicit, punching) ·
`Ladruno_implementation/66_ladruno_solidshell_adr.md` (design record).
