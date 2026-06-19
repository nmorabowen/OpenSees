---
title: "LadrunoConcrete3D — user guide (intended API + modeling guidance)"
project: Ladruno
type: user guide
status: API SPEC (short) — the constitutive law is verified (oracle, through P2 damage + analytic tangent); the C++ KERNEL has the return map, the analytic effective tangent, and the damage stress update (g++-verified); the nDMaterial wrapper is still pending. For the full theory/architecture/usage/quirks reference see LadrunoConcrete3D_guide.
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

> [!warning] Build status (2026-06-19)
> The **constitutive law is fully verified** (numpy oracle through P2 dual damage + the analytic
> damaged tangent, 2 adversarial-review rounds). The **C++ kernel** has the effective return map, the
> analytic effective tangent, **and the dual-damage nominal-stress update** (all g++ byte-verified).
> **Not yet shipped:** the analytic *damaged* tangent in the kernel + the `nDMaterial LadrunoConcrete3D`
> wrapper (classTag 33017). So the model is **not callable from Tcl/Python yet** — this is the intended
> interface; commands below will work once the wrapper lands.

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

## 2. Command (intended)

```tcl
nDMaterial LadrunoConcrete3D $tag $E $nu $fc $ft  \
    <-rho $rho>                                   \
    <-e $e | -fcc $fccRatio>                      \
    <-Df $Df>                                     \
    <-hardening $qh0 $Hp>                         \
    <-ductility $Ah $Bh $Ch $Dh>                  \
    <-Gf $Gf -Gc $Gc>                             \
    <-eta $eta> <-implex>
```
`fc`, `ft` are **positive magnitudes**; the model uses **compression-negative** internally.

## 3. Parameters & defaults

| Param | Meaning | Default |
|---|---|---|
| `E, nu` | Young's modulus, Poisson's ratio | required |
| `fc, ft` | uniaxial compressive / tensile strength (positive) | required |
| `-rho` | mass density | 0 |
| `-e` | deviatoric eccentricity, **∈ (0.5, 1] (hard convexity bound)** | derived from `-fcc` |
| `-fcc` | equibiaxial/uniaxial strength ratio `fcc/fc` (sets `e`) | 1.16 (Kupfer → `e≈0.52`) |
| `-Df` | dilatancy (non-associated flow); `<1` realistic, `=1` max | calibrated (set <1) |
| `-hardening qh0 Hp` | initial yield fraction `qh0`; hardening modulus `Hp` | 0.3, 0.5 |
| `-ductility Ah Bh Ch Dh` | confinement-ductility (Eq.33) — **calibrate from peak strains** | CDPM2 literature defaults |
| `-Gf -Gc` | tensile / compressive fracture energy (**P2 damage**, crack-band) | from `fc` (P2) |
| `-eta` | Duvaut–Lions viscosity (0 = inviscid) | 0 |
| `-implex` | engage Tier-2 IMPL-EX (robust, error-controlled) | off |

**Calibration notes:**
- `e` is a *validation* target, not a fit knob — leave it derived from `-fcc` (1.16) unless you have
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
A symmetric solver (`ProfileSPD`, `BandSPD`, `Mumps` symmetric) will **silently mis-solve** — there
is no runtime guard. This holds even at `Df=1` (associated). For the IMPL-EX tier (`-implex`) the
tangent is SPD and a symmetric solver is fine.

## 5. What you get at each phase

- **Now (P0/P1, oracle-verified):** the **plasticity** — pre-peak nonlinearity, the correct
  triaxial **strength envelope** (Kupfer biaxial, confined `fcc(σ₃)`), and **confinement-dependent
  ductility**. Effective-stress plasticity is **monotonic** — it reaches the failure surface and
  hardens; **there is no stress peak or softening yet.**
- **P2 (damage):** the **peak and post-peak softening** (tension cracking + compression crushing),
  crack-band-regularized (`-Gf`, `-Gc`), with crack/crush pattern output. *This is what turns the
  monotonic plasticity curve into a realistic stress–strain response.*
- **P3+:** robustness tiers (`-implex`, explicit), finite strain (`-geom finite` via `LogStrain`),
  confined-fiber view (1D fibers, "Mander by mechanism").

## 6. Worked example skeleton (once built)

```tcl
# Confined triaxial element test (single brick), unsymmetric solver
model BasicBuilder -ndm 3 -ndf 3
nDMaterial LadrunoConcrete3D 1  30000.0 0.2  30.0 3.0  -fcc 1.16 -Df 0.7
# ... nodes, an 8-node brick (e.g. LadrunoBrick or stdBrick) with material 1 ...
system UmfPack                 ;# REQUIRED — unsymmetric tangent
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-8 50
algorithm KrylovNewton          ;# softening/confined paths want robust iteration
integrator LoadControl 0.01     ;# + adaptive step-cut on non-convergence
analysis Static
analyze 100
```
Notes: confined/softening paths typically need **adaptive step-cutting** (or `-implex`), exactly as
the confinement validation found for ASDConcrete3D. Use core/cover material assignment for columns
(confined view on core fibers once the P5 confined-fiber view lands).

## 7. Recorders (intended, P2+)

Stress/strain/tangent plus the concrete-specific fields once damage lands: `damage` (tension/
compression), `crackWidth`/`crackPattern`, `crushPattern`, `equivalentPlasticStrain`, `kappa_p`,
`implexError`. (P0/P1 expose stress/strain/tangent + `kappa_p`.)

---

**Provenance:** constitutive law pinned to Grassl et al. 2013 (CDPM2, arXiv:1307.6998), verified by a
numpy oracle (`tests/_testbed/concrete3d_ref.py`) and two adversarial review rounds. See
[[LadrunoConcrete3D_dev_handoff]] for the implementation brief and [[31_ladruno_concrete3d_adr]] for
the decision record.
