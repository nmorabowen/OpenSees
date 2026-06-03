---
title: "Ladruno nonlinear RC shell stack — a header-only RC kernel on the ASDShellQ4 + LayeredShellFiberSection frontier"
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - rc-shell
  - constitutive
  - mcft
  - plane-stress
  - compression-softening
  - shear-retention
  - tension-stiffening
  - crack-band
  - solid-shell
  - adr
---

# Ladruno nonlinear RC shell stack — a header-only RC kernel on the ASDShellQ4 + LayeredShellFiberSection frontier

**What.** A fork-authored, header-only, OpenSees-free reinforced-concrete constitutive kernel
(`LadrunoRCKernel.h`) cloned from the proven `LadrunoJ2Kernel.h` "one core, many views" pattern,
consumed through thin `nDMaterial` views that ride the **existing** `ASDShellQ4` (tag 203) +
`LayeredShellFiberSection` seam with at most one adaptor. The kernel keeps `ASDConcrete3DMaterial`'s
verified-correct plastic-damage spine (spectral effective-stress split, dual `dt`/`dc` damage,
Lubliner–Lee–Fenves biaxial envelope, crack-band `lch` regularization, IMPL-EX, crack-closure) and
adds the four physics layers that production RC codes have and `ASDConcrete3D` lacks: compression
(transverse-tension) softening, a degrading aggregate-interlock shear-transfer law, tension
stiffening, and a consistent algorithmic tangent that honors the new strength couplings. An
optional `LadrunoSolidShell` host carries the through-thickness state for punching/bearing — the one
genuine *elemental* blind spot a director shell cannot represent.

**Why.** The deficiency that makes OpenSees over-predict squat-wall diagonal-strut capacity and
miss cyclic pinching is **constitutive, not elemental.** For a wall loaded in its plane the
structural shear `V` lives in the **membrane** block (`gxy`/`tau_xy`), so it is governed by the
plane-stress constitutive law, not by transverse shear and not by element technology. `ASDShellQ4`
(AGQ6-I membrane + 4-DOF EAS + MITC4 transverse shear + Hughes–Brezzi drilling + optional
corotational) already equals the LS-DYNA `ELFORM=16` fully-integrated assumed-strain shell, and
`LayeredShellFiberSection` already does the standard director-stack integration. Rebuilding either
re-treads solved element technology for zero new physics. The **keystone seam decision** is therefore:
keep the element and section at the frontier, deliver the missing physics as a material that drops
into the **5-component `PlateFiber`** layer view the section already requests, and reserve the
finite-strain / `sigma_33`-carrying path for a separate solid-shell host. This document also folds in
the hard findings the adversarial review surfaced — several "free reuse" claims do not survive
contact with the source and are re-scoped below rather than buried.

---

## Background theory — the mechanics

### Generalized shell strain and the director-stack split

The shell elements compute an 8-component generalized strain

```
E = [ e_xx  e_yy  g_xy | k_xx  k_yy  k_xy | g_xz  g_yz ]
      \___ membrane ___/  \___ bending ___/  \_ shear _/
```

and call `Section::setTrialSectionDeformation(E)` (verified `ASDShellQ4.cpp:2009`). The section
integrates through the thickness. `LayeredShellFiberSection` places one midpoint per layer at
`sg[i] = (2i+1)/n − 1`, `wg[i] = 2*t_i/h`, and at signed distance `z_i = 0.5*h*sg[i]` forms the
per-layer plane state `e(z) = e_m − z*kappa`, handing each layer a **5-component `PlateFiber` strain**
`[e_xx e_yy g_xy g_yz g_xz]` (verified `LayeredShellFiberSection.cpp:420-452`). Resultants are
`N = sum w_i sigma_i`, `M = sum w_i z_i sigma_i`, `V = sum w_i tau_i`; the 8×8 tangent is the
transposed assembly `A_sig^T (dd_5x5) A_eps` (`getSectionTangent`, lines 508-670). The constructor
hard-codes `fibers[i]->getCopy("PlateFiber")` and `exit(-1)`s if it returns null
(`LayeredShellFiberSection.cpp:175-180`) — a load-bearing constraint, below.

### Plane-stress condensation (`sigma_33 = 0`)

A 3D law enters a shell layer by condensing `sigma_33 = 0` via a nested scalar Newton on
`eps_33`: `PlateFiberMaterial.cpp:213-258` iterates `Tstrain22 -= condensedStress/dd22` with
`dd22 = threeDtangent(2,2)`, capped at `maxCount = 20`, tolerance `1e-8`. **Critical caveat**
(de Souza Neto §9.4): for a *plastic-damage* law `dd22 -> (1−d)E -> 0` on the softening branch,
so the bare `1/dd22` step explodes or sign-flips, and the loop **returns 0 even when unconverged**.
This is a genuine hazard, addressed in the Decision (the RC `PlateFiber` view must guard and report
non-convergence, not inherit the silent `return 0`).

### MCFT / FSAM compression softening

The verified gap: `ASDConcrete3D`'s `equivalentCompressiveStrainMeasure` uses only the negative
*effective-stress* principals; there is **no** reduction of compressive strength by transverse
tensile cracking strain. The Modified Compression Field Theory (Vecchio & Collins 1986) closes this
with a softening factor on the compressive **strength**:

```
fc_eff = beta * fc' ,    beta = 1 / (0.8 + 170 * eps_1) <= 1 ,    d(beta)/d(eps_1) = -170 * beta^2
```

where `eps_1` is the average principal **tensile** (cracking) strain transverse to the strut. The
DSFM variant (Vecchio 2000) uses `beta = 1/(1 + Cs*Cd)`, `Cd = 0.35*(-eps_1/eps_2 − 0.28)^0.8`.

> **Insertion-point correctness (folded-in fatal finding, D4).** `beta` must scale the **strength /
> stress axis** (`fc'`, equivalently the effective-compressive backbone value `q` fed into `dc`, or a
> contraction of the Lubliner activation surface to `beta*fc'`). It must **not** scale the *strain
> abscissa* fed to `equivalentCompressiveStrainMeasure()`. The abscissa is built from effective-stress
> principals and indexes a non-linear backbone (`ASDConcrete3DMaterial.cpp` ~2403, 2512-2522); scaling
> it merely slides the lookup point along a fixed curve, so the realized peak is `hc(beta*xc)` — an
> uncontrolled value that equals `beta*fc'` only for a linear-through-origin backbone, which concrete
> is not. The acceptance test is a closed-form check: hold confined compression, sweep `eps_1`, assert
> realized peak `|sigma_c| = beta(eps_1)*fc'`.

> **Double-counting against the existing biaxial envelope (folded-in high finding, D4).** The Lubliner
> envelope (`fb = 1.16 fc`, `Kc`, `gamma = 3(1−Kc)/(2Kc−1)`, verified `ASDConcrete3DMaterial.cpp`
> ~2483-2497) already raises the effective-stress measure in the tension–compression quadrant where
> struts live. Stacking `beta` on top double-penalizes that quadrant and can flip the verified
> *over*-prediction into *under*-prediction. With `beta` active, the tension–compression `k1`/`gamma`
> interaction must be reduced or disabled so transverse tension is counted once; `beta` and `Kc` must
> be calibrated **jointly** against a tension–compression panel battery (Kupfer points + Vecchio–Collins).

**OpenSees already ships MCFT prior art** (verified present): `SRC/material/nD/reinforcedConcretePlaneStress/`
(`RAReinforcedConcretePlaneStress` rotating-angle, `FAReinforcedConcretePlaneStress` fixed-angle,
`ConcreteL01/Z01`, `SteelZ01`), `ConcreteMcftNonLinear5/7`, and the `FSAM` material. The gap is
specific to `ASDConcrete3D`, *not* to OpenSees as a whole. These materials are therefore validation
oracles and a "what physics to include" reference, not code to ignore (see Prior art, below).

### Shear retention, rotating vs fixed crack, and the objectivity boundary

> **The shear-retention / crack-model coupling (folded-in high finding, D4/D5).** A pure
> **rotating-crack** model is coaxial — stress and strain share eigenvectors, so there is **no
> independent crack-plane shear stiffness**; cracked shear is fully determined by the principal normal
> stiffnesses (this is exactly why `ASDConcrete3D`'s cracked shear is an emergent spectral byproduct).
> The moment you impose an independent retention `beta_sr * G` on a **stored crack normal `n`**, you
> have re-introduced a fixed crack frame whose response is orientation-history dependent and **fails a
> rigid-rotation objectivity test** — the FSAM/fixed-crack problem the rotating choice was meant to
> avoid. You cannot have both "rotating and coaxial" and "independent `beta_sr*G` on a stored normal."

This forces a clean choice. Cyclic **pinching** physically arises from slip and aggregate-interlock
degradation on a **formed (fixed) crack** across load reversals — a rotating-coaxial model
structurally cannot accumulate it. Since cyclic squat-wall pinching is an explicit goal, the v1
default cannot be plain rotating MCFT. The Decision adopts **fixed-crack-with-degrading-interlock
(FSAM-style) / DSFM (rotating + slip)** as the cyclic core, accepting that the stored crack frame is
a directional internal variable subject to the dSNPO §14.11 large-rotation objectivity boundary
(below), and uses the **corotational element frame** as the supported large-rotation route.

The interlock magnitude can be bounded by the MCFT crack-shear limit
`v_ci,max = 0.18*sqrt(fc') / (0.31 + 24w/(a_g+16))`, crack width `w = eps_1 * s_theta`, with crack
spacing `s_theta` from `lch` or reinforcement spacing.

### Tension stiffening

Between cracks, bonded reinforcement carries average tensile stress, raising the *average* concrete
tension above bare fracture-energy softening:

```
sigma_t_avg = ft / (1 + sqrt(c * eps_1))           (MCFT / Bentz)
sigma_t_avg = alpha1*alpha2*ft / (1 + sqrt(500*eps_1))   (Collins–Mitchell)
```

> **MCFT requires composite `eps_1` (folded-in high finding, D5).** MCFT softening and tension
> stiffening are defined on the **reinforced** average principal tensile strain. If smeared web steel
> lives *outside* the kernel as a separate section layer, the concrete kernel computes `beta` from a
> bare-concrete `eps_1` (no tension stiffening) and over-softens, while the steel layer carries
> tension independently — the two never share the MCFT compatibility equation. Therefore the kernel
> **homogenizes smeared web steel** (or at least receives the rebar strain) for the membrane core;
> external `PlateRebar(LadrunoRebarBuckling)` layers are reserved for **discrete boundary-element
> bars** where buckling matters, not the smeared web that controls strut softening.

### Consistent algorithmic tangent (incl. d(beta)/d(eps_1) and the projector derivative)

The continuum damage-coupled tangent is `D = (I − D_dmg) C0 − sigma~ (x) dD/deps`. The
`sigma~ (x) dD/deps` term is precisely the strength-coupling the secant drops. With `beta = beta(eps_1)`:

```
D_alg = (I − dt_bar PT − dc_bar PC) C0
        − sigma~^- (x) (d dc/d beta)(d beta/d eps_1)(d eps_1/d eps)
        − sigma~^+ (x) (d dt/d eps_eq)(d eps_eq/d eps)
        − (shear/interlock cross term)
        + (dP/deps : sigma~)              [eigenprojector-rotation term]
```

with `d beta/d eps_1 = -170 beta^2` and `d eps_1/d eps = p_1 (x) p_1`.

> **Re-derive against the ACTUAL update, not the secant (folded-in high finding, D4).**
> `ASDConcrete3D`'s operative tangent (`ASDConcrete3DMaterial.cpp:2462-2469`) already uses the
> **nominal** `dt_bar`/`dc_bar` with **fixed** projectors and is itself only an approximate secant —
> it omits `d dt_bar/d eps`, the plastic-strain derivative, the `R`/`mix_dam` coupling
> `d dc/d dt`, and `dP/deps`. The RC tangent must be derived against the real update sequence, not
> bolted onto the misquoted secant, or it still will not be quadratically convergent.

> **Eigenprojector degeneracy (folded-in high finding, D4).** `dP_i/deps ∝ (p_i (x) p_j + p_j (x) p_i)/(lambda_i − lambda_j)`
> is **singular** at equibiaxial / hydrostatic states (`lambda_i -> lambda_j`), which walls (corners)
> and slabs (biaxial bending) visit routinely. Use the Miehe (1993) / de Souza Neto §A
> perturbation–limit regularization with a coalescence tolerance (blend to the symmetric average when
> `|lambda_i − lambda_j| < tol`); a v1 fallback is to keep the **fixed-projector secant** for the
> `dP` term and add only the scalar `beta`/`beta_sr` cross-terms — a defensible, explicitly-stated
> choice rather than an open footnote.

### Crack-band regularization and direction

Bazant–Oh: scale the softening modulus so dissipated energy `g_f = G_f / lch_band` is
mesh-objective. The correct `lch` is the element size **projected onto the crack normal**.

> **`lch` is delivered out-of-band as a single in-plane scalar (folded-in high finding, D2/D3/D4/D5).**
> Verified: the only `lch` channel is `lch = ops_TheActiveElement->getCharacteristicLength()`
> (`ASDConcrete3DMaterial.cpp:1614-1618`), latched **once** (`if (!regularization_done)`), applying the
> **same scalar to tension and compression**; `ASDShellQ4::getCharacteristicLength()`
> (`ASDShellQ4.cpp:1858-1868`) returns the **min nodal distance** (not `sqrt(A)`), **halved when EAS is
> active**, and `setParameter` has **no `lch` case**. There is no per-direction, per-layer, or
> per-crack-normal `lch` channel and an `nDMaterial` cannot override an *element* method. Consequences
> (all real): (a) "directional lch rides the unchanged seam" is **false**; (b) `lch = sqrt(A)` computed
> in the section diverges from the element's EAS-aware value; (c) a single getCharacteristicLength scalar
> mis-regularizes inclined (~45°) wall struts by up to sqrt(2); (d) through-thickness crush-band width is
> not deliverable on this seam. The Decision resolves this explicitly rather than assuming it away.

### Objectivity via the corotational element frame

For pure rigid motion `x_a = Q X_a + c`, the EICR/solid polar fit gives `R = Q`, so deformational
`u_d = R^T(x_a−x_c) − (X_a−X_c) = 0` and the material sees zero strain — frame-indifference holds for
the de-rotated path. The trifecta review established the corotational **element** path is objective
for kinematic hardening because `u_d` is in the reference frame, whereas the `setTrialF` log-strain
**material** wrapper is not (dSNPO §14.11). For RC:

> **Scalar damage is frame-indifferent; the crack frame is not (folded-in high finding, D1/D5).** The
> `dt`/`dc` scalars are objective, but spectral projectors, a stored crack normal, and a directional
> interlock state are directional internal variables that inherit the §14.11 boundary under large
> rotation. The corotational **element** frame (`ASDShellQ4 -corotational`, `LadrunoBrick -geom corot`)
> is the supported objective route for moderate-to-large wall rotation. The `LadrunoRCFiniteStrain`
> (`setTrialF`) view is **objective only for the isotropic-damage spine**; its fixed-crack/interlock
> directional state under large rotation is pinned **xfail** until a co-rotating-crack finite-native
> view exists.

> **Single solid corotation is not a shell corotation (folded-in high finding, D1).** For a thin
> solid-shell (`t << L`), `SolidTransformationCorot` extracts **one** element-average `R` from a polar
> fit of the nodal cloud `H = sum (x_a−x_c)(X_a−X_c)^T`; the eigenvalue ratio of `S=(H^T H)^{1/2}`
> scales as `(L/t)^2`, so the bending rotation is ill-conditioned and a single `R` cannot remove the
> through-thickness-varying bending rotation that is the dominant shell strain. The solid corot is only
> reliable in moderate-rotation / thicker regimes; thin large-rotation cases route to `-geom finite` or
> a shell-aware corotation, with a `cond(S)` guard in `update()`.

### LS-DYNA comparison (manual citations and portable lessons)

LS-DYNA independently arrived at the **"thin element, thick material"** division this ADR adopts. The
portable lessons (Theory Manual, R-versions; sections verified against the Theory Manual, per-`MAT`
Vol II page numbers **flagged unverified** below):

- **Belytschko–Lin–Tsay `ELFORM=2`** (Theory §7): one-point co-rotational Mindlin shell with viscous
  perturbation **hourglass control** (shape vector `tau_I`, eq 7.15; hourglass stress rates eq
  7.18a-c scaled by `r in [0.01,0.05]`). Resultant integration `f^R = ∫ sigma dz`, `m^R = −∫ z sigma dz`
  (eqs 7.13a/b) is **identical** to `LayeredShellFiberSection`. Lesson: any reduced-integration shell
  host must degrade hourglass stiffness with damage — exactly the `LadrunoBrick`
  `Kstab <- max(1%, 1−max(dt,dc))` pattern.
- **Hughes–Liu `ELFORM=1`** (Theory §10): degenerated-brick director shell; lamina frame enforces
  `sigma_33 = 0` (§10.2.2), Hughes–Carnoy iterative thickness-thinning (eqs 10.37-10.39). This is the
  director shell the panel explicitly declines to re-tread.
- **Fully-integrated `ELFORM=16`** (Theory §9): Hu–Washizu three-field, Pian–Sumihara assumed in-plane
  strain (anti-locking), Bathe–Dvorkin MITC transverse shear (§9.4), Belytschko–Leviathan drill
  projection (eq 9.27). **Direct analog of `ASDShellQ4`** (AGQ6-I/EAS + MITC4 + drilling) — confirming
  `ASDShellQ4` is already at the LS-DYNA frontier.
- **Layered-shell transverse shear** (Theory §11): FOSDT gives a constant `gamma_xz` that violates
  zero-traction faces and "could yield very stiff behavior in sandwich and laminated shells"; LS-DYNA
  reconstructs a parabolic `tau_xz` from `d tau_xz/dz = −d sigma_x/dx` (eqs 11.10-11.14). **Portable
  lesson with a sharp caveat (folded-in finding, D6):** this is a **1D-bending postprocessor** needing
  an in-plane *stress gradient* the per-Gauss-point material API does not own. It is valid as an
  **output recovery** (`setResponse('tauxz_profile')`), **not** as a tangent-bearing constitutive law;
  it cannot ride the existing seam "with one adaptor." v1 transverse shear is a **local degrading shear
  law**, not the gradient-coupled equilibrium profile.
- **Concrete cards.** MAT_084/085 Winfrith (smeared crack, explicit cracked-plane **shear-retention**
  factor), MAT_172 Concrete_EC2 (layered RC shell, explicit tension stiffening + EC2 compression
  curve), MAT_159 CSCM (cap + damage), MAT_072R3 K&C (three-surface, eta-interpolated residual =
  degrading shear), MAT_273 CDPM2 (two damage variables — directly comparable to `ASDConcrete3D`'s
  `dt`/`dc`). **Verdict:** every production RC card carries an explicit degrading-shear and/or
  compression-softening term that `ASDConcrete3D` lacks — independent confirmation that the four missing
  layers are exactly what a "great" RC material needs.

> **Honesty on LS-DYNA framing (folded-in finding, D6).** LS-DYNA corroborates the *architectural
> division* (prior art for "thin element, thick material"); it does **not** "prove the kernel correct."
> Kernel correctness is established only by the experiment + oracle battery. The Vol II per-`MAT` page
> numbers (MAT_159 ~p.1139, MAT_172 ~p.1208, MAT_072R3 ~p.582, MAT_273 ~p.1890) are carried from the
> brief and **must be spot-checked against the on-disk R15/R16 PDFs**, or cited by `MAT` number +
> manual version only, before appearing as durable citations.

---

## Decision — the architecture across six dimensions

**D1 — Hosts.** Two hosts, one seam. (a) **Reuse `ASDShellQ4` (203) unmodified** as the primary RC
shell host for walls and slabs; reuse `ASDShellT3` (204) and `ShellNLDKGQ` (157) as-is. (b) **Build
`LadrunoSolidShell` (classTag 33020)** — 8-node, 3 translational DOF/node, full 3D state carrying
`sigma_33`, with mandatory EAS-on-`E33` + ANS/MITC transverse shear + ANS membrane/trapezoidal cures —
as a **narrow specialist for punching/bearing/3D-stress and large strain**, not a co-equal flexural
host. *Rejected:* a new director shell (duplicates `ASDShellQ4`/Felippa AFEM Ch.31-36 for zero physics);
a 6-DOF drilling solid-shell (defeats brick/contact parity, cannot use `SolidTransformationCorot`);
bolting `sigma_33` onto `ASDShellQ4` (corrupts the proven Reissner–Mindlin element); full-integration
solid-shell (locks catastrophically — Belytschko Ch.8).

> Re-scoped from the proposal: (1) **`LadrunoBrick`'s "EAS" is NOT a reusable template.** It is the
> SSPbrick stabilized-single-point scheme that condenses enhanced modes **once at setup with the initial
> elastic tangent** into a **constant** `Kstab` (`LadrunoBrick.h:270-281`). Genuine EAS-on-`E33` for
> softening RC needs **persistent committed `alpha[nEAS]`, per-Newton condensation with the consistent
> *damaged* tangent**, and serialized send/recvSelf — net-new Simo–Rifai code, not a port. A constant
> `Kstab` would leave the thickness mode elastically stiff exactly where the element must soften.
> (2) The single-layer solid-shell has ~2 z-Gauss points; it cannot resolve cracked-RC `sigma(z)` or the
> migrating neutral axis, so it ships with **selectable multi-layer / Gauss–Lobatto `n_z`** and is
> benchmarked vs `LayeredShellFiberSection` before any flexural claim. (3) `ASDShellQ4`(6-DOF) ↔
> solid-shell(3-DOF) edges lose moment continuity; a documented, **validated** rigid-link/`equalDOF`
> connection recipe is a D1 deliverable.

**D2 — Deformation seam.** The **8-component generalized-strain → `Section::setTrialSectionDeformation`
→ per-layer 5-component `PlateFiber`** path is the default and only mandatory seam. **Do not route `F`
through the section** (the Section has no Gauss-point `F`, no current geometry, and cannot assemble the
element-owned geometric stiffness — first-principles Belytschko split, corroborated by
`LayeredShellFiberSection.cpp:420-504` only ever seeing the 8-vec and assembling an 8×8). The `F`-seam
(`FiniteStrainNDMaterial::setTrialF`, `LogStrainNDMaterial`) is a **separate optional host** for the
solid-shell and `LadrunoBrick -geom finite`. *Rejected:* section-chooses-reduce-vs-passthrough; pass-`F`
everywhere (pays the eigen-log spatial-tangent cost on every layer for zero in-plane-shear gain, inherits
§14.11); the generic `PlateFiberMaterial` wrapper (redundant condensation on top of the kernel's own).
**Verified:** `FiniteStrainNDMaterial.h` and `LogStrainNDMaterial.{h,cpp}` are present on this branch —
the solid-shell `F`-view is reuse, not a dependency on unshipped assets.

**D3 — Section.** Use `LayeredShellFiberSection` essentially as-is via the `PlateFiber` view (below);
keep **midpoint-per-layer** integration as the default with **8–12 layers** for nonlinear RC bending
(warn below 6); offer per-layer `-quad lobatto2` only for homogeneous/elastic layers. The section seam
stays **strictly order-5 `PlateFiber`** — any `PlaneStress` order-3 view lives **inside** the RC material,
never at the section boundary, and never via `PlateFromPlaneStress`'s elastic-`gmod` transverse shear.

> Re-scoped from the proposal (D3): (1) **A subclass that does not override `getCopy()` is SLICED to
> the base type.** `LayeredShellFiberSection::getCopy()` hard-constructs `new LayeredShellFiberSection(...)`
> and `ASDShellQ4.cpp:740` clones per-Gauss-point — a subclass overriding "only ~3 methods" silently runs
> vanilla with the new classTag. If a fork section subclass is authored, it **must** override
> `getCopy()`/`sendSelf`/`recvSelf` + broker registration (the full `MovableObject` contract), or compose
> instead of inherit. **v1 avoids this entirely:** the RC `PlateFiber` view drops into the *unmodified*
> `LayeredShellFiberSection`; no fork section is needed for Phase 1. (2) **Smeared rebar must not be a
> rho-weighted overlay at a concrete layer's `z_i`** (double-counts concrete `(1)c + rho*s` instead of
> `(1−rho)c + rho*s`); discrete boundary rebar is its own thin layer at the true bar depth with the
> overlapping concrete reduced. (3) `setResponse('damage')` on a shell section is a **recorder/diagnostic**
> output, **not** a `Kstab` coupling — `ASDShellQ4` has no damage-scalable hourglass stiffness analogous
> to `LadrunoBrick`. (4) `lch_t = t_i` is **not** a safe compression default — it triggers snapback and the
> material's `gmin` clamp silently discards it; floor it at a physical crush-band width.

**D4 — RC kernel physics.** Author `LadrunoRCKernel.h` keeping the `ASDConcrete3D` spine verbatim and
adding: compression softening `beta = 1/(0.8+170*eps_1)` **applied to the strength/stress axis** (with the
Lubliner tension–compression interaction reduced to avoid double-counting); a **fixed-crack /
DSFM-with-slip** degrading aggregate-interlock shear law (chosen over plain rotating-coaxial **because
cyclic pinching needs a fixed plane**); opt-in tension stiffening (default = fracture-energy, so flags-off
reduces to baseline); the **full consistent tangent** re-derived against the real update with a
**regularized eigenprojector derivative** (or fixed-projector-secant v1 fallback); IMPL-EX with
**clamped extrapolated `eps_1` and `beta in [floor,1]`**; and crack-band regularization fed by the `lch`
resolution chosen in D5. *Rejected:* secant-only tangent (drops the real `beta` Jacobian → squat-wall
divergence); editing `ASDConcrete3DMaterial.cpp` (answers only `ThreeDimensional`, forfeits the
multi-view + ledger story); replicating MAT_172/MAT_159 verbatim (closed, knob-heavy).

**D5 — Software architecture.** Header-only OpenSees-free `LadrunoRCKernel.h` (namespace
`ladruno_rc_kernel`), numpy-oracle-testable before any OpenSees link, cloning `LadrunoJ2Kernel.h`.
**Views and classTags:**

| Class | classTag | `getType()` | order | host / status |
|---|---|---|---|---|
| `LadrunoRCPlateFiber` | 33014 | `PlateFiber` | 5 | **the shell path** — rides `LayeredShellFiberSection` `getCopy("PlateFiber")` |
| `LadrunoRCPlaneStress` | 33013 | `PlaneStress` | 3 | 2D continuum / membrane quad **only**; oracle exerciser — **NOT a shell host** |
| `LadrunoRCFiniteStrain` | 33015 | `ThreeDimensional` | — | optional solid-shell / `-geom finite`; xfail under large-rotation directional state |

> Re-scoped from the proposal (D5): (1) **`ASDShellQ4` consumes a `SectionForceDeformation`, not an
> `nDMaterial`** (`ASDShellQ4.h:128`, `.cpp:740`). The "`PlaneStress` drops directly into `ASDShellQ4`"
> claim is false — the **shell path is the `PlateFiber` view only.** `PlaneStress` (33013) hosts a 2D
> continuum element and is the clean oracle target; it is shipped for that, not as a wall solution. (2)
> **"Byte-identical across all three views" is downgraded.** A native plane-stress map and a 3D map of a
> *directional* (crack-normal-bearing) law do not commute with `sigma_33=0` condensation or large-rotation
> log-strain. The accurate guarantee: `returnMapPS` is byte-shared between the `PlaneStress` view and the
> `PlateFiber` in-plane block; the FiniteStrain view uses a **distinct 3D map** validated independently.
> Identity-gate regressions are **per-view reduce-to-`ASDConcrete3D`**, not cross-view. (3) **`getCopy(type)`
> self-routing mints sibling classTags** — centralize history pack/unpack in a shared header so all
> `sendSelf`/`recvSelf` use one order; per-view db-roundtrip + MP parity tests (the `LadrunoJ2` family had
> to harden exactly this). Defer the speculative `33016` reservation until a real brick consumer exists.

**D5 — `lch` resolution (the binding cross-cut).** Pick and own **one**:
- **Option A (v1 default, recommended):** accept the element's **scalar in-plane** `lch` via
  `element->getCharacteristicLength()` (which already encodes `ASDShellQ4`'s EAS `/2` correction);
  treat through-thickness crush regularization as the **layer thickness handled in the material from a
  value the layer already knows**; **document that out-of-plane bending-crack energy is not
  through-thickness size-objective on the director-shell host** (that physics belongs to the solid-shell).
- **Option B (when inclined-crack mesh-objectivity is required):** add a small, **Ladruno-tagged,
  `LEDGER_vanilla_files`-recorded** edit to `LayeredShellFiberSection`/`ASDShellQ4` plumbing a
  per-direction `lch` to the layer — and **drop the "zero vanilla edit" claim** for that case.
There is **no silent `lch` default**: a mesh-objectivity test fails loudly if the scalar fallback is used
in a softening run.

**D6 — External + validation.** LS-DYNA is prior art for the architecture and a physics reference, not a
correctness proof. The validation battery (below) maps to the two-zone testbed, **corrects the SFI-MVLEM
oracle misidentification** (it is **FSAM/fixed-strut**, verified `FSAM.cpp:10-22` — *not* MCFT; it
brackets rather than confirms a rotating-angle kernel; `ConcreteMcftNonLinear`/`RAReinforcedConcretePlaneStress`
are the genuine in-code MCFT oracles), and **documents the punching blind spot** as an xfail.

---

## Architecture & interfaces

### Kernel/views layout

```
SRC/material/nD/
  LadrunoRCKernel.h          # header-only, OpenSees-free; namespace ladruno_rc_kernel
  LadrunoRCConcreteLaws.h    # shared backbones (beta softening, tension stiffening, interlock) + RCHist serializer
  LadrunoRCSteel.h           # smeared web-steel homogenization (composite eps_1 for MCFT)
  LadrunoRCPlateFiber.{h,cpp}    # ND_TAG 33014  (the shell path)
  LadrunoRCPlaneStress.{h,cpp}   # ND_TAG 33013  (2D continuum / oracle)
  LadrunoRCFiniteStrain.{h,cpp}  # ND_TAG 33015  (solid-shell / finite; FiniteStrainNDMaterial subclass)
SRC/element/ladrunoSolidShell/
  LadrunoSolidShell.{h,cpp}      # ELE_TAG 33020  (optional through-thickness host)
```

### Kernel contract (clone of `LadrunoJ2Kernel.h`)

```cpp
namespace ladruno_rc_kernel {
  enum { STATUS_OK = 0, STATUS_SINGULAR = 1, STATUS_NO_CONVERGE = 2 };
  struct Params {
    double E, nu, fc, ft, eps_c0, Gf_t, Gf_c, Kc, ag;
    int  crackModel;          // 0 = fixed/DSFM-slip (default for cyclic), 1 = rotating (monotonic-only)
    int  shearRet;            // const | dsfm | rots
    int  tensStiff;           // fracture(default) | vc | cm
    int  tangentMode;         // algorithmic(default) | secant | numerical
    bool implex;
    double betaSrMin;         // aggregate-interlock residual floor
    double lch_in, lch_thick; // resolved per D5 Option A/B
    int    nWebRebar;         // smeared web steel homogenized into the membrane core
    double rho[8], theta[8];  // web ratios/angles (MCFT composite eps_1)
  };
  struct RCHist { /* xt_max, xc_max, xt_pl, xc_pl; crack_frame; interlock state; eps1_commit; implex commit set */ };

  int returnMapPS (const Params&, const double eps3[3], const RCHist&,
                   double sig3[3], double Dtan3[3][3], RCHist&, double* outResidual = 0);  // {exx,eyy,gxy}
  int returnMap3D (const Params&, const double eps6[6], const RCHist&,
                   double sig6[6], double Dtan6[6][6], RCHist&, double* outResidual = 0);

  inline double betaCompr (double e1) { return 1.0 / (0.8 + 170.0 * e1); }   // applied to STRENGTH axis
  inline double dBetaCompr(double e1) { double b = betaCompr(e1); return -170.0 * b * b; }
}
```

### View / NDMaterial contract

- `getType()` returns `"PlateFiber"` (33014, order 5 `{exx,eyy,gxy,gyz,gxz}`), `"PlaneStress"` (33013,
  order 3), `"ThreeDimensional"` (33015).
- `getCopy(const char* type)` self-routes to the sibling view (centralized `RCHist` (de)serializer so all
  three pack identically). The shell path requires `getCopy("PlateFiber")` to return a real order-5 clone —
  `LayeredShellFiberSection` `exit(-1)`s otherwise (the hard requirement, `:175-180`).
- `LadrunoRCPlateFiber`: in-plane block via `returnMapPS`; the `sigma_33=0` condensation uses a
  **guarded, line-searched / bisection-fallback** inner Newton that **propagates a non-convergence code**
  (must not inherit `PlateFiberMaterial`'s silent `return 0`); transverse shear is a **real degrading
  local law**, not elastic `gmod`.
- `LadrunoRCFiniteStrain : public FiniteStrainNDMaterial` — `setTrialF(F)` → Hencky → kernel → Cauchy +
  spatial tangent; objective for the isotropic spine only.
- `setResponse`: `"damage"->(dt,dc)`, `"beta"`, `"betaShear"`, `"eps1"`, `"crackNormal"`,
  `"stress"`, `"strain"`, `"tangent"`, `"tauxz_profile"` (output-recovery only).

### How it rides the existing seam

`ASDShellQ4 -corotational` → `LayeredShellFiberSection` (unmodified) → per layer
`getCopy("PlateFiber")` → `LadrunoRCPlateFiber` (33014) → `returnMapPS` + transverse-shear block.
Boundary-element rebar = discrete `PlateRebar(LadrunoRebarBuckling)` layers at true bar depth; smeared web
steel is **inside** the kernel for MCFT composite `eps_1`.

---

## Implementation path

Each phase is independently shippable and reduces to baseline when its flags are off.

### Phase 1 — Compression softening on the EXISTING element/section (closes the headline gap)
- **Build:** `LadrunoRCKernel.h` + `LadrunoRCPlateFiber` (33014) with the spine cloned from
  `ASDConcrete3D`, plus `beta = 1/(0.8+170*eps_1)` **on the strength axis** (Lubliner t–c interaction
  reduced to avoid double-counting) and the consistent tangent's scalar `d beta/d eps_1` cross-term
  (fixed-projector-secant v1 for `dP`). `lch` per **D5 Option A** (scalar in-plane). Smeared web steel
  homogenized for composite `eps_1`.
- **Reuses:** `LadrunoJ2Kernel.h` pattern; `ASDConcrete3D` spine; `LayeredShellFiberSection` +
  `ASDShellQ4` unmodified (`getCopy("PlateFiber")` seam); MCFT prior-art materials as oracles;
  `PlateFiberMaterial` condensation pattern (hardened).
- **New:** the kernel, the order-5 view, the guarded `sigma_33=0` inner Newton.
- **Deliverable:** `nDMaterial LadrunoRCPlateFiber` usable in a `LayeredShellFiberSection` under
  `ASDShellQ4` with **zero element/section edit**.
- **Acceptance:** (i) closed-form `|sigma_c| = beta(eps_1)*fc'` under swept transverse tension;
  (ii) reduce-to-`ASDConcrete3D` on the raw `ThreeDimensional` path to ~1e-6/1e-7 *trajectory* tolerance
  (not byte-identity — `ASDConcrete3D` default tangent is damaged-secant, `ASDConcrete3DMaterial.cpp` ~1728);
  (iii) Kupfer biaxial-envelope regression (no double-counting overshoot); (iv) forward-difference tangent
  check at pre-peak, post-tensile-peak, deep-compression-softening, and a rotating-axis state;
  (v) hardened condensation: a strain history crossing the compressive peak (negative `dd22`) asserts inner
  convergence + `|sigma_33| < tol` + a propagated non-convergence code.

### Phase 2 — Degrading shear / cyclic pinching
- **Build:** fixed-crack/DSFM-slip degrading aggregate-interlock law (`-shearRetention`) + crack-closure;
  bounded by `v_ci,max`.
- **Reuses:** Phase-1 kernel; `ASDConcrete3D` crack-closure spectral reassembly.
- **New:** the interlock/slip state and its tangent cross-term.
- **Deliverable:** cyclic membrane shear law with pinching.
- **Acceptance:** cyclic shear-panel hysteresis vs an MCFT oracle and experiment — assert **pinching shape +
  cumulative hysteretic energy**, not just peak; **ablation** softening-ON/retention-OFF vs both-ON proving
  the retention term is load-bearing; **rigid-rotation objectivity** test on a cracked state (the
  stored-normal form must pass the supported corotational-element route).

### Phase 3 — Tension stiffening + crack-band/`lch` hardening
- **Build:** VC/CM tension-stiffening plateaus (opt-in); resolve `lch` per **D5 Option A or B**.
- **New:** if Option B, the ledgered vanilla `lch` plumb.
- **Deliverable:** slab/distributed-reinforcement tension stiffening; size-objective in-plane softening.
- **Acceptance:** in-plane mesh-objectivity on an **inclined-crack (rotated) mesh** (peak ~1–3%, energy
  ~3–5% across 2× refine, calibrated on the elastic shell patch first — **not** imported from the brick
  numbers); a loud failure if a scalar `lch` fallback is used in softening.

### Phase 4 — Finite-strain view + IMPL-EX
- **Build:** `LadrunoRCFiniteStrain` (33015) via `setTrialF`/`LogStrainNDMaterial`; IMPL-EX with clamped
  `eps_1`/`beta` and an implex-vs-implicit error monitor + step-cut.
- **Reuses:** `FiniteStrainNDMaterial`, `LogStrainNDMaterial` (verified present); `LadrunoJ2Finite` pattern.
- **Deliverable:** large-strain RC for the solid-shell; symmetric-solver explicit-friendly cyclic runs.
- **Acceptance:** IMPL-EX `O(dt)` on a **smooth** proportional path (excluding damage-activation/closure
  steps); cyclic energy-balance band on the wall; rigid-rotation-of-a-cracked-point **xfail** (§14.11).

### Phase 5 — `LadrunoSolidShell` (33020) — optional through-thickness host
- **Build:** 8-node, 3-DOF, **genuine state-dependent EAS-on-`E33`** (persistent `alpha`, per-Newton
  consistent-damaged-tangent condensation, serialized) + ANS/MITC transverse shear + ANS membrane;
  selectable multi-layer/Gauss–Lobatto `n_z`; directional/projected `lch`; `cond(S)` corot guard.
- **Reuses:** `SolidTransformation` linear/corot/finite; `LadrunoRCFiniteStrain`; `LadrunoBrick` damage-
  scaled `Kstab` (as stabilization, re-tuned for thin shells, **not** as the EAS template).
- **Deliverable:** punching/bearing/3D-stress RC host with `sigma_33`.
- **Acceptance:** elastic patch + pinched-cylinder + Scordelis-Lo (element correctness) **and** a
  **softening snap-back** (mesh-objective dissipation, Newton/arc-length convergence with the
  incomplete-geometric + secant tangent), a slab **punching** benchmark, and an EAS internal-mode
  growth-stability check under post-peak softening.

---

## Validation plan (mapped to the two-zone testbed)

**Zone-A (upstreamable pytest):**
1. Membrane/bending/shear constant-stress **patch** tests through `ASDShellQ4` + view (elastic) → ~1e-8.
2. **Objectivity, split:** (a) pure rigid rotation via the corotational **element** path → must PASS;
   (b) superposed finite-rotation + deviatoric strain via the FiniteStrain **material** view → **xfail**
   (§14.11 directional-state mechanism).
3. **reduce-to-`ASDConcrete3D`** identity gate (flags off) on the **raw `ThreeDimensional`** view, stress
   *trajectory* ~1e-6/1e-7 (bypasses condensation; floor set by damaged-secant tangent + condensation
   residual).
4. **`beta` strength oracle** — closed-form `|sigma_c| = beta(eps_1)*fc'`.
5. **Forward-difference tangent** at pre-peak / post-tensile-peak / deep-compression / rotating-axis /
   **equibiaxial** (finite, symmetric — guards the eigenprojector regularization).
6. **Condensation robustness** — `PlateFiber`/`PlaneStress` `sigma_33=0` Newton across the compressive
   peak (negative `dd22`): inner convergence + residual + non-convergence code (the most likely hidden
   source of false global non-convergence).
7. **IMPL-EX vs implicit** `O(dt)` on a smooth proportional path; **Kupfer** biaxial-envelope no-overshoot.

**Zone-B (gmsh, fork-local):**
8. **Squat shear wall** (aspect <~1.5), cyclic — vs **FSAM-backed SFI-MVLEM as a fixed-angle bracket** and
   an **MCFT in-code oracle** (`ConcreteMcftNonLinear`/`RAReinforcedConcretePlaneStress`), with a named
   **physical experiment** (e.g. Tran–Wallace RW-A20-P10 or a PEER specimen) as the **primary** gate —
   asserting diagonal-strut capacity **and** pinching shape + hysteretic energy.
9. **Slender flexural wall** (aspect >~2) — fiber flexural backbone + boundary-element rebar buckling via
   `PlateRebar(LadrunoRebarBuckling)`.
10. **Two-way slab nonlinear bending** vs yield-line/experiment.
11. **In-plane mesh-objectivity** on a **rotated (inclined-crack)** notched panel + the chosen `lch`
    resolution; loud failure on scalar fallback in softening.
12. **Punching-shear blind-spot — xfail with a written mechanism note + `LEDGER_quirks` entry:** a
    director/condensed-`PlaneStress` shell carries constant transverse shear and no `sigma_33` and
    **cannot form a punching cone**; refinement cannot manufacture the missing kinematics. The documented
    resolution is `LadrunoSolidShell` (Phase 5).
13. *(Optional)* `tauxz_profile` **output recovery** vs the LS-DYNA §11 reconstruction (recorder only).

---

## Risks & open questions

> [!question] **(FATAL→mitigated) `beta` insertion point — strength axis, not strain abscissa.**
> Scaling the `equivalentCompressiveStrainMeasure` abscissa does not realize `fc_eff = beta*fc'` on a
> non-linear backbone. Mitigated by acting on the strength/stress axis with the closed-form acceptance
> test (Phase 1-i). *If that test cannot pass, Phase 1 is blocked* — this is the single most important
> correctness gate.

> [!question] **(FATAL→avoided) `LadrunoBrick` EAS is not a solid-shell EAS template.** Its constant,
> initial-tangent `Kstab` cannot soften. `LadrunoSolidShell` (Phase 5) budgets genuine state-dependent
> Simo–Rifai EAS-on-`E33` as net-new code; do not re-cost it as a port.

> [!question] **(HIGH) Rotating vs fixed crack is forced by the cyclic goal.** Plain rotating-coaxial
> has no fixed plane to accumulate slip → cannot pinch; an independent `beta_sr*G` on a stored normal is
> fixed-crack masquerading as rotating and fails rigid-rotation objectivity. v1 = fixed-crack/DSFM-slip,
> accepting the §14.11 directional-state boundary and using the corotational **element** as the objective
> route. *Open:* MCFT vs DSFM as the calibration default — pick against the squat-wall battery.

> [!question] **(HIGH) `lch` is a single in-plane scalar on this seam.** No per-direction/per-layer/per-
> crack-normal channel exists, and `ASDShellQ4` halves it under EAS. v1 = **Option A** (scalar in-plane,
> through-thickness bending-crack objectivity explicitly out of scope on the director host). Option B
> (a ledgered vanilla plumb) only if inclined-crack objectivity must be exact. *Open:* which option ships
> in Phase 3, and the crash-band floor `lch_t >= max(t_i, k*d_agg)`.

> [!question] **(HIGH) Consistent tangent vs eigenprojector degeneracy.** The default-algorithmic tangent
> must be re-derived against the real `ASDConcrete3D` update (nominal `dt_bar`/`dc_bar`, `R`/`mix_dam`
> coupling, plastic-strain derivatives) and must regularize `dP/deps` at equibiaxial states (Miehe limit)
> — else the default path produces indefinite/NaN tangents in states walls/slabs routinely visit. v1
> fallback = fixed-projector secant for `dP` + scalar `beta` cross-terms, explicitly stated. *Open:*
> whether to ship the full regularized projector derivative or the fallback in Phase 1.

> [!question] **(HIGH) MCFT needs composite `eps_1`.** Smeared web steel must be homogenized **inside**
> the kernel; external `PlateRebar` layers are discrete boundary bars only. Reconcile this everywhere
> (the proposal text was internally inconsistent on rebar location).

> [!question] **(HIGH) The shell path is `PlateFiber`-only.** `ASDShellQ4` consumes a section, not an
> `nDMaterial`; the `PlaneStress` view is a 2D-continuum/oracle material, not a wall solution. All wall
> shear flows through the order-5 `PlateFiber` condensation — the `beta`/shear/objectivity analysis is
> done in the 5-component condensed setting, not a free-standing 3×3.

> [!question] **(HIGH→narrowed) Solid-shell corotation + single-layer flexure.** The single nodal-cloud
> `R` is ill-conditioned for thin shells (`cond(S) ~ (L/t)^2`) and 2 z-Gauss points cannot resolve
> cracked-RC flexure. `LadrunoSolidShell` ships as a punching/bearing/3D-stress specialist with selectable
> `n_z`, a `cond(S)` guard, and a moment-curvature benchmark vs `LayeredShellFiberSection` before any
> flexural claim — **not** a co-equal flexural host. *Open:* shell-aware corotation vs `-geom finite` for
> thin large-rotation cases; whether `SolidTransformationCorot` needs its two deferred geometric tangent
> terms (validated on a **softening** snap-back, not the elastic pinched-cylinder) before softening-RC use.

> [!question] **(MEDIUM) "Byte-identical across views" is downgraded.** A directional law's native PS map
> and 3D map do not commute with condensation/log-strain. The guarantee is `returnMapPS` shared between
> the `PlaneStress` view and the `PlateFiber` in-plane block; the FiniteStrain view is validated
> independently with a reduce-to-plane-stress consistency regression.

> [!question] **(MEDIUM) IMPL-EX is `O(1)` at damage-activation/closure.** Certify `O(dt)` only on smooth
> paths; on the cyclic wall assert energy/peak bands and exclude activation/closure steps. Clamp
> extrapolated `eps_1` and `beta`; monitor implex-vs-implicit error with step-cut.

> [!question] **(MEDIUM) SFI-MVLEM is FSAM, not MCFT** (verified `FSAM.cpp:10-22`). It brackets a
> rotating-angle kernel rather than confirming it; the physical experiment is the primary squat-wall gate
> and `ConcreteMcftNonLinear`/`RAReinforcedConcretePlaneStress` are the in-code MCFT oracles.

**Nonlocal / E-FEM boundary.** A **local** crack-band RC kernel composes for free into the existing seam.
**Nonlocal / gradient-damage** regularization needs neighbor-state averaging or an extra nodal Helmholtz
field — OpenSees has no general nonlocal section/element machinery, no integration-point neighbor map at the
section level, and no extra-DOF gradient field on shells. The **E-FEM** embedded-discontinuity path (ADR #18)
needs condensed internal DOF + an enriched section/element. Both are **out of scope** for this material-only
stack and belong to a future section/element ADR; the RC kernel does not pretend to deliver them.

**Objectivity caveat (restated).** Scalar `dt`/`dc` are frame-indifferent; the stored crack frame and
interlock direction are directional internal variables bounded by dSNPO §14.11 under large rotation. The
corotational **element** path is the supported objective route; the FiniteStrain material view is xfail for
the combined large-rotation + directional-state case until a co-rotating-crack finite-native view is built.

**Provenance note.** `FiniteStrainNDMaterial.h`, `LogStrainNDMaterial.{h,cpp}`, `LadrunoJ2Kernel.h`, the
`reinforcedConcretePlaneStress/` family, `ConcreteMcftNonLinear5/7`, and `FSAM.{h,cpp}` are all verified
present on this branch. The LS-DYNA Vol II per-`MAT` page numbers are **not yet spot-checked** and must be
verified or reduced to `MAT`-number citations before the ADR is finalized.

---

## Design discussion log

| Dim | Decision | Strongest adversarial finding (severity) | Resolution |
|---|---|---|---|
| **D1** Hosts | Reuse `ASDShellQ4`; build `LadrunoSolidShell` (33020) as a punching/bearing specialist | `LadrunoBrick` EAS is constant initial-tangent SSP, not a softening-EAS template (**fatal**) | Strike the reuse claim; spec genuine state-dependent Simo–Rifai EAS-on-`E33` as net-new; demote solid-shell from co-equal flexural host to 3D-stress specialist with selectable `n_z` |
| **D2** Seam | 8-vec generalized strain → section → 5-comp `PlateFiber`; `F` is a separate optional host | Directional `lch` cannot ride the unchanged seam; `F`-asset files unverified (**high/medium**) | Files verified present; lead with the first-principles Belytschko split; resolve `lch` via D5 Option A/B (drop "zero edit" if Option B) |
| **D3** Section | `LayeredShellFiberSection` as-is via `PlateFiber`; midpoint, 8–12 layers | A subclass not overriding `getCopy()` is **sliced** to vanilla on per-GP clone (**fatal**) | v1 needs no fork section (view drops into the unmodified section); if subclassed, full `MovableObject` contract; rebar as own thin layer; `setResponse('damage')` = recorder only; floor `lch_t` |
| **D4** Kernel | Spine + strength-axis `beta` + fixed-crack interlock + re-derived tangent + clamped IMPL-EX | `beta` wired to the strain **abscissa**, not the strength axis — does not realize `fc_eff=beta*fc'` (**fatal**) | Apply `beta` to the strength/stress axis with a closed-form acceptance test; reduce Lubliner t–c interaction to avoid double-counting; regularize the eigenprojector derivative |
| **D5** Software | Header-only kernel + 3 views (33013/33014/33015); shell path = `PlateFiber` | "`PlaneStress` drops into `ASDShellQ4`" is false (section, not `nDMaterial`); RC physics already exists in OpenSees (**high**) | Shell path is `PlateFiber`-only; `PlaneStress` = 2D/oracle; add a prior-art section + use those materials as oracles; downgrade "byte-identical across all views"; centralize `RCHist` serialization |
| **D6** External + validation | LS-DYNA = architecture prior art; 10-class battery; punching xfail | **SFI-MVLEM is FSAM (fixed-strut), not MCFT** — breaks the squat-wall isolation logic (**high**) | Correct the oracle: FSAM brackets, MCFT in-code materials confirm, experiment is primary; demote reduce-to-baseline from bit-identity to trajectory tolerance; add cyclic-energy/condensation/equibiaxial/objectivity gates; rescope §11 reconstruction to output-recovery |

---

## References

**Textbooks**
- T. Belytschko, W. K. Liu, B. Moran, *Nonlinear Finite Elements for Continua and Structures* — Ch. 3
  (strain/stress measures), Ch. 8 (locking, ANS/EAS, B-bar), Ch. 9 (plates/shells, plane stress).
- C. A. Felippa, *Advanced Finite Element Methods (AFEM)* — Ch. 22–23 (ANS transverse shear,
  Hughes–Brezzi drilling), Ch. 31–36 (director / 4–6 DOF shell); EICR (Felippa–Haugen).
- E. A. de Souza Neto, D. Perić, D. R. J. Owen, *Computational Methods for Plasticity* (2008) — Ch. 6
  (softening/regularization), Ch. 9 (plane-stress-projected return mapping, §9.4), Ch. 12 (damage),
  Ch. 14 (multiplicative finite strain, Box 14.3 MATISU, §14.11 objectivity boundary), App. A
  (eigenprojector derivatives / Miehe limit).
- F. J. Vecchio, M. P. Collins, "The Modified Compression-Field Theory for Reinforced Concrete Elements
  Subjected to Shear," *ACI Journal* 83(2), 1986; F. J. Vecchio, "Disturbed Stress Field Model (DSFM),"
  *J. Struct. Eng.* 126(9), 2000.
- T. T. C. Hsu, Y.-L. Mo, *Unified Theory of Concrete Structures* (CSMM / fixed-angle softened truss);
  Z. P. Bažant, B. H. Oh, "Crack band theory," *Mat. & Struct.* 16, 1983.

**LS-DYNA manuals** (Theory sections verified; per-`MAT` page numbers flagged for spot-check)
- LS-DYNA *Theory Manual* (R15/R16): §7 Belytschko–Lin–Tsay shell + hourglass control; §9 fully-integrated
  Hu–Washizu/MITC shell (`ELFORM=16`); §10 Hughes–Liu degenerated-brick shell (`ELFORM=1`); §11 layered-shell
  equilibrium transverse-shear reconstruction (eqs 11.10–11.14).
- LS-DYNA *Keyword User's Manual* Vol II: MAT_084/085 Winfrith, MAT_159 CSCM, MAT_172 Concrete_EC2,
  MAT_072R3 K&C, MAT_273 CDPM2 (cite by `MAT` number + manual version pending page verification).

**OpenSees source (verified this session, paths relative to `SRC/`)**
- `element/shell/ASDShellQ4.{h,cpp}` — `:128`/`:740` consumes a `SectionForceDeformation`;
  `:1858-1868` `getCharacteristicLength` (min-distance, EAS `/2`); `:2009` `setTrialSectionDeformation`.
- `material/section/LayeredShellFiberSection.cpp` — `:175-180` hard `getCopy("PlateFiber")` + `exit(-1)`;
  `:222` base-type `getCopy()`; `:420-504` 8→5 map + resultant integration; `:508-670` 8×8 tangent.
- `material/nD/PlateFiberMaterial.cpp` — `:213-258` `sigma_33=0` nested Newton (`return 0` even unconverged).
- `material/nD/ASDConcrete3DMaterial.cpp` — spectral split + Lubliner envelope (~2452-2497); secant tangent
  (2462-2469); IMPL-EX (2383-2398, default off ~402); `lch` latch (1614-1618); `regularize`/`gmin` (947-990).
- `material/nD/FiniteStrainNDMaterial.h`, `material/nD/LogStrainNDMaterial.{h,cpp}` — `setTrialF` / Hencky
  seam (present on branch).
- `material/nD/LadrunoJ2Kernel.h` — header-only "one core, many views" template.
- `material/nD/FSAM.{h,cpp}` — **Fixed-Strut-Angle-Model** backing SFI-MVLEM (`:10-22`).
- `material/nD/reinforcedConcretePlaneStress/` (`RA`/`FA`ReinforcedConcretePlaneStress, `ConcreteL01/Z01`,
  `SteelZ01`) and `material/nD/ConcreteMcftNonLinear5/7.{h,cpp}` — in-code MCFT prior art / oracles.
- `element/LadrunoBrick.{h,cpp}` — `:270-281` constant initial-tangent SSP `Kstab` (the non-template);
  damage-scaled `Kstab <- max(1%,1−max(dt,dc))`.
