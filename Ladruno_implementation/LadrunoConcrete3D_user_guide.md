---
title: "LadrunoConcrete3D — user guide (API + modeling guidance)"
project: Ladruno
type: user guide
status: SHIPPED — `nDMaterial LadrunoConcrete3D` (classTag 33017) is callable from Tcl/Python on `ladruno` (kernel g++-verified through P3b analytic damaged tangent; wrapper verified on Zone-A CI). ALL dimensional views ship (3D + PlaneStrain/AxiSymmetric/PlateFiber/PlaneStress, #299); cyclic temper + robustness tiers deferred. For the full theory/architecture/usage/quirks reference see LadrunoConcrete3D_guide.
related:
  - "[[LadrunoConcrete3D_guide]]"          # the COMPLETE four-part reference
  - "[[LadrunoConcrete3D_dev_handoff]]"   # the implementer's brief
  - "[[31_ladruno_concrete3d_adr]]"
  - "[[19_ladruno_rc_shell_adr]]"          # LadrunoRCConcrete (shell/MCFT sibling)
updated: 2026-06-19
---

# LadrunoConcrete3D — user guide (short)

> [!info] This is the short API spec. The **complete reference** (theory · architecture · usage ·
> quirks) is **[[LadrunoConcrete3D_guide]]**.

> [!success] Build status (2026-06-19) — SHIPPED
> `nDMaterial LadrunoConcrete3D` (classTag 33017) is **callable from Tcl and Python** on `ladruno`. The
> C++ kernel is g++ byte-verified against the numpy oracle through the analytic damaged tangent (P3b),
> and the wrapper is verified end-to-end on Zone-A CI (full build + an openseespy stdBrick battery).
> **ALL dimensional views ship** — 3D + the Phase-2 reduced PlaneStrain / AxiSymmetric / PlateFiber /
> PlaneStress (one `dim`-mode class, reached via the element's `getCopy(type)`; finite-strain via
> `nDMaterial LogStrain`). The cyclic temper and the robustness tiers (`-eta`/`-implex`) are deferred.

## 1. What it is, and when to use it

`LadrunoConcrete3D` is a **CDPM2-grade** (Grassl et al. 2013) 3D **solid** concrete model:
effective-stress **Menétrey–Willam plasticity** (smooth 3-invariant surface with Lode dependence and
a compression cap) + **dual scalar damage** (tension/compression). It is the model to reach for when
**confinement, triaxiality, and loading-path generality** matter.

**Choose it over the alternatives when:**

| Need | Use |
|---|---|
| Solid continuum, confined columns/joints/deep beams, varying triaxiality, high confinement, smeared confinement | **`LadrunoConcrete3D`** (this) |
| In-plane RC walls / shells, MCFT membrane shear | `LadrunoRCConcrete` (33015) — the shell sibling |
| Robust general 3D concrete, calibrated-band confinement (`p/fc ≲ 0.2`), already good | `ASDConcrete3D` (upstream) |

vs **ASDConcrete3D**: ASD captures confinement *emergently* (calibrated band, ≈Mander ±5% to
`p/fc=0.2`) but is single-meridian, capless, with no dilation knob. `LadrunoConcrete3D` makes
confinement **constitutive** — the **ductility measure** is the mechanism generalization of Mander's
pressure→ductility law, so it extrapolates across triaxiality and high confinement.

## 2. Command

```tcl
nDMaterial LadrunoConcrete3D $tag $E $nu $fc $ft $Gf $Gc  \
    <-e $e | -kupfer $fccRatio>                           \
    <-Df $Df>  <-As $As>  <-rho $rho>                     \
    <-hardening $qh0 $Hp>                                 \
    <-ductility $Ah $Bh $Ch $Dh>                          \
    <-lch $lch>  <-autoRegularization>  <-implex>
```
`E ν fc ft Gf Gc` are **positional and required**; `fc`, `ft` are **positive magnitudes** (the model
uses **compression-negative** internally). `-autoRegularization` pulls the crack-band length from the
parent element (mesh-objective); otherwise `-lch` sets it (default 1.0).

## 3. Parameters & defaults

| Param | Meaning | Default |
|---|---|---|
| `E, nu` | Young's modulus, Poisson's ratio (positional) | required (`E>0`, `0≤ν<0.5`) |
| `fc, ft` | uniaxial compressive / tensile strength, **positive** (positional) | required (`0<ft<fc`) |
| `Gf, Gc` | tensile / compressive fracture energy, crack-band (positional) | required (`>0`) |
| `-e` | deviatoric eccentricity, **∈ (0.5, 1] (hard convexity bound)** | derived from `-kupfer` |
| `-kupfer` | equibiaxial/uniaxial strength ratio `fcc/fc` (sets `e`) | 1.16 (Kupfer → `e≈0.52`) |
| `-Df` | dilatancy (non-associated flow); `<1` realistic, `=1` max | 1.0 (set `<1`, e.g. 0.85) |
| `-As` | softening ductility (`x_s` in compression); `≥1` | 2.0 |
| `-rho` | mass density | 0 |
| `-hardening qh0 Hp` | initial yield fraction `qh0`; hardening modulus `Hp` | 0.3, 0.5 |
| `-ductility Ah Bh Ch Dh` | confinement-ductility (Eq.33) — **calibrate from peak strains** | 0.08, 0.003, 2.0, 1e-6 |
| `-lch` / `-autoRegularization` | crack-band length: fixed / from the element | 1.0 / off |
| ~~`-eta`/`-implex`~~ | Duvaut–Lions viscosity / Tier-2 IMPL-EX | **not in v1 (P3)** |

**Calibration notes:**
- `e` is a *validation* target, not a fit knob — leave it derived from `-kupfer` (1.16) unless you have
  measured biaxial data; it lands at the canonical `e≈0.52`.
- `-ductility` parameters set the **absolute peak strain and its growth with confinement**; the
  literature defaults give the right *trend* but recalibrate `Ah/Bh/Ch/Dh` to your concrete's measured
  peak strains. Keep `Hp` small (large `Hp` can make the hardening non-monotone).
- `-Gf`/`-Gc` (P2) are crack-band-regularized → results are **mesh-objective in dissipated energy**.

## 4. Solver requirement — READ THIS

The consistent tangent is **non-symmetric** (non-associated flow + the spectral return). You **must
use an unsymmetric solver**:
```tcl
system UmfPack           ;# or FullGeneral
```
A symmetric solver (`ProfileSPD`, `BandSPD`, `Mumps` symmetric) will **silently mis-solve** — there is
no runtime *enforcement*, but the parser **prints a warning at material creation**. This holds even at
`Df=1` (associated). (The future IMPL-EX tier is SPD ⇒ a symmetric solver will be fine there.)

## 5. What you get — shipped vs deferred

- **Shipped (v1):** the full **monotonic** response — pre-peak plasticity, the triaxial **strength
  envelope** (Kupfer biaxial, confined `fcc(σ₃)`), **confinement-dependent ductility**, **the peak +
  tension/compression softening + automatic unilateral recovery**, crack-band-regularized — with the
  **analytic damaged tangent** (Tier-1 implicit). Effective-stress plasticity is monotonic; the peak and
  softening are the **damage** part. **All dimensional views** — 3D + PlaneStrain / AxiSymmetric /
  PlateFiber / PlaneStress (the element picks one via `getCopy(type)`) + the finite-strain view via
  `nDMaterial LogStrain`. (PlaneStress/PlateFiber enforce σ_zz=0 by a nested ε_zz Newton + static
  condensation; unconfined plane-STRESS post-peak softening is snap-backy → robust pre-peak / confined.)
- **Deferred:** cyclic (`β_c` + the compression→tension temper, P2f), robustness tiers (`-eta`/`-implex`,
  explicit — P3), the confined-fiber 1D view ("Mander by mechanism").

## 6. Worked example skeleton

```python
import openseespy.opensees as ops
ops.model("basic", "-ndm", 3, "-ndf", 3)
# ... 8 nodes of a brick + restraints ...
ops.nDMaterial("LadrunoConcrete3D", 1, 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0,  # E nu fc ft Gf Gc
               "-Df", 0.85, "-autoRegularization")
ops.element("stdBrick", 1, 1,2,3,4,5,6,7,8, 1)   # or LadrunoBrick / SSPbrick
ops.system("FullGeneral")          # REQUIRED — unsymmetric tangent (or UmfPack)
ops.constraints("Transformation"); ops.numberer("RCM")
ops.test("NormDispIncr", 1.0e-8, 100); ops.algorithm("Newton")
ops.integrator("DisplacementControl", 2, 1, d)   # drive a face; step-CUT past the softening peak
ops.analysis("Static")
```
Notes: **unconfined softening is snap-backy** on a single implicit element — drive with **adaptive
step-cutting** (halve the increment on a failed `analyze`), use a confined cell, or an arc-length /
indirect-control integrator. Compression hardens through its pre-peak range so plain fixed-step
DisplacementControl reaches `−fc` fine; tension reaches its limit point at once (onset `ε0=ft/E`).

## 7. Recorders

Per-Gauss-point via the element's `material` response — `stress`, `strain`, `tangent`,
`effectiveStress` (the undamaged `σ̄`, record alongside `stress` to see the damage knock-down),
`damage` (`[ω_t, ω_c]`), `kappaP`, `plasticStrain`:
```python
ops.eleResponse(eleTag, "material", gp, "damage")     # -> [omega_t, omega_c]
```

---

**Provenance:** constitutive law pinned to Grassl et al. 2013 (CDPM2, arXiv:1307.6998), verified by a
numpy oracle (`tests/_testbed/concrete3d_ref.py`) and two adversarial review rounds. See
[[LadrunoConcrete3D_dev_handoff]] for the implementation brief and [[31_ladruno_concrete3d_adr]] for
the decision record.
