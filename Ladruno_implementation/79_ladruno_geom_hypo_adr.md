# ADR 79 — `-geom hypo`: hypoelastic rate-form updated-Lagrangian geometry method

**Status:** COMPLETE — P0+P1 merged (#679), P2 merged (#680), P3 PDMY smoke shipped (#681). The staged program of this ADR is done, and the bearing-limit campaign it was built to answer has been **run and reported in §8**: the geometry ladder produces no bearing limit point and hardens the backbone slightly, refuting the ADR-78 §1 hypothesis.
**Element:** `LadrunoBrick` (ELE 33002) first; `LadrunoUP` (ELE 33017) in P2
**Seam extended:** `SolidTransformation` gains `METHOD_HYPO = 3` (marker only)
**Kernel:** `SRC/element/solidTransformation/LadrunoHypoKernel.h` (header-only, OpenSees-free)

---

## 1. Why

ADR 78 opened `-geom corot` on the u–p family: large ROTATION + geometric
stiffness, small deformational strain. The remaining gap is genuine LARGE
STRAIN for soils. `-geom finite` (the `setTrialF(F)` path) is blocked on the
material side: the `LogStrainNDMaterial` lifting adaptor states its own limit —
*"v1 assumes a LINEAR-elastic inner law (so Cᵉ = inv(initial tangent) is
constant and the εᵉ recovery is exact)"* — and the target soil model
`PressureDependMultiYield` is pressure-dependent hypoelastic (G ∝ p′^0.5),
which violates that assumption directly. No finite-strain soil model exists in
the fork or upstream, and writing one is a multi-year constitutive program.

The industry answer — what Abaqus NLGEOM does for every rate-form
(hypoelastic) material (Abaqus Theory Guide §1.4.3 strain-increment
integration, §1.5.3 objective rates, §1.5.4 state storage; via the
`abaqus-theory-foundations` skill) — is an **objective-rate updated-Lagrangian
integration**: compute an objective (rigid-rotation-free) strain INCREMENT per
step from midpoint-configuration gradients (Hughes & Winget 1980), hand it to
the SAME small-strain constitutive routine, keep the material's tensorial
state consistent across rotation, and assemble on the CURRENT configuration
with the initial-stress geometric term. PDMY-class soil laws are already
rate-form — they consume strain increments and produce stress increments —
so this route gives large strain WITH the existing soil library, unchanged.

---

## 2. The design decision: wrapper vs seam vs element lane

Three candidate homes for the hypo path were on the table. The deciding
constraint fell out of a fact about OpenSees materials, not about elements:

**The classic Hughes–Winget prescription is unimplementable against a stock
`NDMaterial`.** HW (and Abaqus's UMAT contract) rotate the stored stress AND
every stored tensorial internal variable (backstress, surface centers) forward
by the incremental rotation ΔR at the top of each step. That works when the
host code owns the state arrays (Abaqus stores STRESS/STATEV and hands them to
the routine). An OpenSees `NDMaterial` owns its state privately —
`PressureDependMultiYield`'s nested-surface centers, a J2's backstress — and
exposes no "rotate your tensors by ΔR" entry point. Adding one per material
would mean touching the whole library, which is exactly what this ADR exists
to avoid.

**Resolution: use the unrotated-configuration (Green–Naghdi) variant of the
same objective integration** (Flanagan–Taylor lineage; the PRONTO/Sierra
approach; Abaqus itself uses the Green–Naghdi rate for all Explicit VUMAT and
structural elements — Theory Guide Table 1.5.3-1). Instead of co-rotating the
material state forward each step, DE-ROTATE the objective strain increment
into a fixed material frame:

    Δε_mat = Rᵀ_mid · Δε_spatial · R_mid ,     R = polar(F)

and let the material live its entire life in that fixed unrotated frame —
its stress, its backstress, its surface centers never need rotating, because
its input never rotates. The stress is pushed forward `σ = R σ_mat Rᵀ` only at
assembly time. Objectivity class is identical (Jaumann and Green–Naghdi differ
only under combined finite rotation + finite shear; GN's simple-shear response
is monotone where Jaumann's famously oscillates — a strictly better ceiling).
And the fixed-material-frame convention is exactly the one `-geom corot`
already established (`SolidTransformationCorot.h` header notes: the material
and its backstress live in a fixed reference frame), so the kinematic-hardening
objectivity test pattern carries over unchanged.

With that settled, the three candidates:

1. **A fourth `SolidTransformation` method (seam-resident logic).** Rejected.
   The seam's LOCKED division of labour gives the transformation exactly four
   jobs (strain measure, localizeDisp, globalizeForce/Stiff, recorder frame) on
   whole-element vectors. The hypo path needs PER-GAUSS-POINT machinery —
   midpoint gradients, per-GP polar rotation, per-GP accumulated strain — that
   the seam interface cannot see (its `update()` receives only nodal
   coordinates). Cramming per-GP state into the wrapper would break the
   "formulation owns B" contract and make the wrapper element-shape-aware.
   Precedent: `-geom finite` hit the same wall and became an element lane with
   `SolidTransformationFinite` as a pure MARKER (all seam methods identity).

2. **A per-GP material wrapper (`HypoIncrementalNDMaterial`).** Rejected. The
   wrapper would need Δε and R per step — data outside the `NDMaterial`
   interface — so the element must `dynamic_cast` and call custom methods
   anyway; the wrapper buys no decoupling. It costs a user-facing extra object
   (`nDMaterial HypoWrap 2 1` around every soil material — the LogStrain
   friction, again), a full recorder/response passthrough layer (the SEAM-1
   lesson), and an ownership split of the per-GP state. With the unrotated-
   frame resolution above, there is nothing material-side left to wrap: the
   inner material is driven by plain `setTrialStrain()` with an accumulated
   feed — the element can do that directly.

3. **An element-internal lane (the `-geom finite` pattern).** **Chosen.**
   `SolidTransformationHypo` is a marker (METHOD_HYPO = 3, all seam methods
   identity, `getStrainMeasure() == SmallStrain` — honest: the material IS
   driven via `setTrialStrain`), parse/serialize rides the existing method-id
   plumbing, and the element implements `updateHypo()` +
   `formResidAndTangentHypo()` beside `updateFinite()` +
   `formResidAndTangentFinite()`, reusing the UL assembly skeleton (spatial
   gradients via F⁻¹, current-config dv = J·dV, initial-stress geometric term).
   The per-GP math lives in a header-only, OpenSees-free kernel
   (`LadrunoHypoKernel.h`, the `LadrunoUPKernel.h` discipline) so P2 reuses it
   from `LadrunoUP` without touching the brick.

**User experience result (the point of the whole exercise):**

    nDMaterial PressureDependMultiYield 1 ...
    element LadrunoBrick 1  n1..n8  1 -geom hypo

— no adaptor object, no material change, large strain.

---

## 3. The algorithm (per Gauss point, per evaluation)

State: committed accumulated feed strain `ε_feed_n` (6 Voigt components/GP —
the ONLY persistent hypo state) plus whatever the material commits internally.
Committed nodal displacement `u_n` comes from `Node::getDisp()` (committed);
trial `u_{n+1}` from `getTrialDisp()`. Every Newton iteration recomputes from
committed state — no accumulation-within-iteration hazard.

1. **Deformation gradients** from reference shape gradients `∂N/∂X`:
   `F_n = I + Σ u_n,a ⊗ ∂N_a/∂X`, `F_{n+1}` likewise, `F_mid = ½(F_n + F_{n+1})`
   (Hughes–Winget's midpoint configuration `x̄ = ½(x_n + x_{n+1})`).
2. **Midpoint spatial gradients** `∂N_a/∂x̄ = (∂N_a/∂X)·F_mid⁻¹`; incremental
   displacement gradient `G = Σ Δu_a ⊗ ∂N_a/∂x̄`, `Δu = u_{n+1} − u_n`.
3. **Objective strain increment** `Δε_spatial = sym(G)`. HW's midpoint makes
   this EXACTLY zero for an arbitrarily large rigid-rotation increment (the
   sym part of `2·tan(Δθ/2)·ŵ` vanishes identically — verified in the kernel
   oracle), not merely O(Δθ²); the O(Δθ³/step) error quoted for HW applies to
   combined stretch+rotation paths.
4. **Midpoint rotation** `R_mid = polar(F_mid)` (Higham iteration, the corot
   wrapper's algorithm). De-rotate: `Δε_mat = R_midᵀ · Δε_spatial · R_mid`.
   Using the TOTAL polar rotation (not a ΔR product chain) keeps R stateless,
   drift-free, and restart-safe.
5. **Accumulate and feed:** `ε_feed = ε_feed_n + Δε_mat` (Voigt, engineering
   shear), then `setTrialStrain(ε_feed)`. The material perceives exactly the
   increment `Δε_mat` relative to its committed state — the Green–Naghdi-rate
   constitutive integration, by construction, on an UNCHANGED material.
6. **Assembly** (current configuration, `R_{n+1} = polar(F_{n+1})`):
   - push forward `σ = R σ_mat Rᵀ`;
   - `f_a = ∫ σ·∇ₓN_a dv`, `dv = J·dV`, spatial gradients via `F_{n+1}⁻¹`;
   - `K = ∫ BᵀcB dv + ∫ GᵀΣG dv` with `c` = the material 6×6 tangent
     bond-rotated by `R_{n+1}` and the standard symmetric initial-stress term
     `(K_geo)_{(a,i)(b,k)} = δ_ik ∫ ∂Nₐ/∂x_j σ_jl ∂N_b/∂x_l dv`;
   - body force stays a dead load per REFERENCE volume (the finite-lane
     convention, consistent with reference-config mass).

`commitState()` copies the trial feed to `ε_feed_n`; revert restores it;
`sendSelf` ships it (§8).

### 3.1 Conventions, stated honestly

- **Rate law on Cauchy stress** (`σ̇_GN = C:d`), the FLAC/geotech convention.
  No `J`/work-conjugacy (Kirchhoff) correction: at 20 % volumetric strain the
  difference from a Kirchhoff-based law is O(10 %) — that is the accepted
  hypoelastic trade, and it is what PDMY's calibration assumes anyway (soil
  moduli are fitted to true-stress/engineering-ish data at exactly these
  strains). The P1 gate vs `-geom finite`+LogStrain MEASURES this 1/J factor
  and asserts it, turning the discrepancy into a verified property.
- **Hypoelastic elasticity is not hyperelastic**: closed strain cycles through
  large strain dissipate/generate O(ε³) energy. Documented ceiling; irrelevant
  for the driving application (soils are dissipative from the first increment).
- **Simple shear at very large γ**: the Green–Naghdi rate is monotone (no
  Jaumann oscillation) but stiffens vs the analytic hyperelastic reference for
  γ ≳ 1. Gate 5 pins the element against the numpy oracle integrating the SAME
  algorithm — element-vs-oracle is exact; oracle-vs-physics is the documented
  hypoelastic ceiling.
- **Tangent policy**: `∫BᵀcB + ∫GᵀΣG` omits the rate-convection terms of the
  exact GN consistent tangent (O(σ) vs O(C), unsymmetric). Same "residual
  exact, tangent dominant-term" policy as corot v2.0 / ADR 78 §3.3; K stays
  SYMMETRIC (a real solver advantage over F-bar finite). Newton converges to
  the correct equilibrium; only the asymptotic rate is affected.
- **Recorder frame**: `stress`/`strain` material responses report in the
  UNROTATED material frame (σ_mat, ε_feed) — the corot convention, and what
  the objectivity gates want to see invariant. `frameTimeVarying()` stays
  false in v1 (per-GP frames don't fit the element-level recorder contract);
  a global-frame push-forward channel is a demand-driven follow-up.

---

## 4. Scope and guards (v1)

- `-formulation std` only. bbar/F-bar-style mean dilatation in rate form is a
  known extension (apply the mean-dilatation split to `Δε` and the assembly B)
  but is NOT needed by the driving application: drained PDMY is not
  incompressibility-locked, and undrained saturation arrives via the u–p
  element in P2, where incompressibility is carried by the pressure FIELD, not
  by the skeleton operator. uri/ssp/eas: reserved (same status as corot).
- Requires a SMALL-STRAIN material — a `FiniteStrainNDMaterial` is rejected at
  parse time (the existing GEOM-1 converse guard already does this for any
  non-finite method).
- `-damp` reserved (mirror of the finite guard); `-bulkViscosity` stripped
  (existing non-linear-geom guard).
- 3D H8 lane only in v1 (the kernel is dimension-generic; 2D and Tet lanes are
  demand-driven).
- `-geom linear` stays bit-identical: every hypo branch is behind
  `isHypo()`, which is false unless METHOD_HYPO was parsed.

---

## 5. Staging

- **P0 — kernel + oracle (this PR):** `LadrunoHypoKernel.h` (midpoint
  gradients, Higham polar, Δε build/de-rotation, Voigt bond rotation of σ and
  C) + numpy oracle `tests/hypo_reference.py` + g++ cross-check
  `tests/test_hypo_kernel_cpp.py` (the `test_ladruno_up_kernel_cpp.py`
  pattern: oracle emits cases, C++ harness diffs to ≤1e-9, plus standalone
  invariants — rigid-increment Δε ≡ 0, polar orthogonality, rotation-of-
  identity-tangent invariance).
- **P1 — dry solid (this PR):** the hypo lane on `LadrunoBrick` (H8, std).
  Gates (tests/test_ladrunoBrick_hypo.py) — ALL GREEN; measured:
  cantilever tip at ~26 % deflection hypo=1.5783 vs corot=1.5738 (0.29 %) vs
  finite=1.5806 (0.15 %); simultaneous-rotation J2-kin objectivity error
  6.6e-3 (the predicted O(Δθ·Δε) ≈ 5e-3 scale); rigid 0.5 rad/step to 2.5 rad
  stress-free at 1e-7·E; uniaxial-strain 25 % exact vs (λ+2μ)lnλ at 5e-4 and
  the finite+LogStrain ratio = J to 1e-3; shear γ=2 element==oracle at 1e-8:
  1. rigid-rotation zero stress through large multi-step rotation (tight —
     §3 step 3 makes it near-identity, asserted at solver tolerance);
  2. `-geom linear` unchanged + `-geom hypo` ≡ linear at infinitesimal strain;
  3. agreement with `-geom corot` at moderate strain under rotation;
  4. uniaxial-strain elastic compression to 25 %: exact vs the analytic
     hypoelastic target `σ = (λ+2μ)·ln λ`, and the measured gap to
     `-geom finite`+LogStrain asserted ≈ the 1/J factor (§3.1);
  5. elastic simple shear to γ = 2 vs the numpy oracle (exact) — documents
     the GN signature (monotone, no Jaumann oscillation);
  6. J2 kinematic-hardening cyclic objectivity in the hypo path's OWN
     convention (fixed material frame — the corot test pattern): (a)
     deform-then-rigid-rotate freezes the material-frame stress EXACTLY (the
     HW rigid-increment property), (b) simultaneous rotation+stretch agrees
     to the documented incremental-objectivity error scale;
  7. database round-trip: the guarded Vector(48) feed-strain payload survives
     save/restore — the post-restore step matches the uninterrupted run to
     1e-9.
- **P2 — the u–p element (second PR; #677 merged 2026-07-28, unblocking it):**
  `-geom hypo` on `LadrunoUP`, reusing the kernel. Per-block decisions (the
  ADR-78 §3 discipline):
  - **u-rows** — the P1 brick lane verbatim: spatial total-stress internal
    force `∫Bᵀ(σ′_pushed − αp·m) dv` on the current configuration (the −αp leg
    against the SAME spatial gradients IS −Q_cur·p, so K/Q frame-AND-
    configuration consistency is structural — ADR 78 §3.2 carried over);
    `K = ∫BᵀcB dv +` the total-stress symmetric geometric term.
  - **p-row coupling** — the ADR-78 GN22/GN11 incremental pairing verbatim,
    with the coupling increment now the TRUE volume change: per GP
    `α·tr(Δε)·dv/Δt`. `tr(Δε)` sums to `ln J` in the midpoint limit (the
    Δln J continuity) and is EXACTLY zero on a rigid increment — the ADR-78
    chord defect cannot arise by construction. S (evolved) + αH̃ go
    incremental with it (mixing schemes on the p-row is the measured ADR-78
    pump). Tangent: Q_cur / Q_curᵀ dominant-term blocks.
  - **H/S** — rebuilt per evaluation on the CURRENT configuration:
    `H = ∫(∇ₓN_p)ᵀ k_sp (∇ₓN_p) dv` with the skeleton-attached permeability
    pushed spatially, `k_sp = kc(J)·R·diag(k̄)·Rᵀ`;
    `S = ∫N_pᵀ (1/Q̄(n(J))) N_p dv`. The KINEMATIC porosity
    `n(J) = 1 − (1−n₀)/J` (clamped) is always on under hypo — it is
    kinematics, not a model. **Kozeny–Carman `k(J)` scaling is a constitutive
    CHOICE and OPT-IN via `-kozenyCarman`** (normalized to 1 at n₀, floored at
    1e-6 so a crushed GP degrades to near-impermeable without singularizing
    H). The equal-order H̃ stays the reference unit-α Laplacian (a
    stabilization term; O(ε) effect).
  - **Seepage drive** — spatial gradients are already global-frame: no
    pull-back (the corot Rᵀ exists because corot assembles in the material
    frame). Full-tensor k_sp loop replaces the diagonal-k kernel call.
  - **Mass / body / initial stiffness / Rayleigh** — reference-configuration
    (mass and dead load per reference volume keep the summed-reaction weight
    identity exact; `getInitialStiff` keeps R=I/J=1 semantics per ADR 78
    §3.6; Rayleigh C stays the reference build — per-GP rotations define no
    single element frame, the βK parts carry a documented O(rotation)
    approximation in the damping forces only).
  - **Scope v1**: 3D H8 equal-order std lane only (parser + ctor sanitizer
    both police; 2D plane-strain embedding, TH vertex-p, bbar-in-rate-form
    reserved). Serialization: payload widened 30→31 (`-kozenyCarman` flag) +
    the committed feed strain as a guarded extra `Vector(nGP·6)` (the P1
    pattern). New response `hypoState` = per-GP `[J, n(J), kcScale]`.
  - **Gates** (tests/test_ladruno_up_element_hypo.py): parser guards; static
    + transient-undrained rigid rotation (the tr(Δε) rigid-exactness);
    consolidation reduces to linear (refinement-aware); drained equivalence
    vs the P1 dry hypo brick (roller sides, `-b` conversion, summed-reaction
    weight check); porosity/J + KC-formula identities via `hypoState` +
    KC-slows-consolidation direction; undrained large-strain storage
    `p ≈ Q̄·α·|ln λ|`. A Gibson-family closed-form large-strain consolidation
    oracle is deferred (demand-driven) — the physics identities above pin the
    same machinery without the PDE oracle.
- **P3 — PDMY01 smoke (SHIPPED, tests/test_ladruno_up_hypo_pdmy.py):** a
  TIMs-style footing-in-miniature (3×3×3 saturated PDMY box) through the full
  staging idiom (elastic gravity → `updateMaterialStage 1` → plastic settle →
  DISPLACEMENT-controlled footing push to 5 % penetration) under `-geom hypo`
  with `-geom corot` as the anchor twin. Two scoping findings, both recorded
  in the test docstring:
  1. the push MUST be displacement-controlled — a force-controlled footing
     push on stage-1 PDMY diverges from the first increment (near-surface
     GPs at ~zero confinement, G ∝ p′^d ⇒ near-singular tangent ⇒ unbounded
     first iterate; hypo correctly refuses the inverted iterate, linear
     grinds without converging);
  2. `-geom linear` is not a usable twin on this problem (it grinds even
     displacement-controlled) — corot anchors, and the hypo-vs-corot gap IS
     the measured geometric-nonlinearity content: 0.3 % at 0.25 %
     penetration, 1.8 % at 1 %, 15.4 % at 5 % (hypo STIFFER — UL assembly on
     the compacted configuration + n(J) storage), J ∈ [0.984, 1.022].
  Every push step runs through the PDMY-battery substep-fallback policy
  (KrylovNewton dt/10) — the designed mechanism, not an anomaly. The full
  bearing-limit-point campaign on the real SFIM mesh is analysis work, not
  CI.

---

## 6. Serialization

Per-GP persistent hypo state = committed `ε_feed_n` only (8 GP × 6 = 48
doubles). Shipped as an EXTRA guarded `Vector(48)` after the standard `dData`,
gated on the geometry method id in `idData(28)` (the EAS-alpha pattern —
non-hypo streams are byte-identical; eas+hypo cannot co-occur, parser-refused,
so the guarded-send order is deterministic). `R` is stateless (recomputed from
nodal positions), so nothing else ships. The idData(28) packing
(`+ 1000·methodID`) already accommodates id 3 without widening.

---

## 7. Files

| file | change |
|---|---|
| `SRC/element/solidTransformation/SolidTransformation.h` | `METHOD_HYPO = 3` |
| `SRC/element/solidTransformation/SolidTransformationHypo.h/.cpp` | marker method (identity seams) |
| `SRC/element/solidTransformation/SolidTransformationLinear.cpp` | `create()` case 3 |
| `SRC/element/solidTransformation/LadrunoHypoKernel.h` | P0 kernel (header-only) |
| `SRC/element/ladrunoBrick/LadrunoBrick.h/.cpp` | hypo lane: state, updateHypo, formResidAndTangentHypo, commit/revert/send/recv |
| `SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp` | `-geom hypo` parse + guards |
| `tests/hypo_reference.py` | numpy oracle (emits `tests/hypo_cases.txt`) |
| `tests/hypo_kernel_check.cpp` + `tests/test_hypo_kernel_cpp.py` | P0 g++ cross-check |
| `tests/test_ladrunoBrick_hypo.py` | P1 gates |
| `SRC/element/ladrunoUP/LadrunoUP.h/.cpp` | P2 hypo lane: state, updateHypo, current-config blocks, incremental p-row, hypoState response, serialization 30→31 + guarded feed |
| `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp` | P2: `-geom hypo` parse + guards, `-kozenyCarman` |
| `tests/test_ladruno_up_element_hypo.py` | P2 gates |

No upstream file is touched; no new class tag (the seam is tag-free by
design).

---

## 8. Bearing-backbone campaign — result (2026-07-29)

The ADR-78 §1 question this ADR was built to answer — *does genuine large
strain bend the PDMY footing backbone toward a bearing limit point, or reduce
the ~9.3×-Vesic over-hardening?* — has been measured. Full write-up, CSVs and
plot: `Ladruno_files/testbed/hypo_bearing/RESULTS.md`.

**Answer: no, and no — the ladder hardens the backbone slightly instead.**

The real SFIM model is not in this repo (and is not on the machine), so the
campaign runs a self-contained graded benchmark built with apeGmsh: B = 2 m
smooth square footing on a 20 × 20 × 8 m box (4.5 B side clearance, 4 B depth,
which is what the scoping work showed is required to escape the oedometer
regime), 2816 `LadrunoUP` H8 at 0.5 m under the footing grading outward,
13 872 u-p DOF, saturated PDMY (φ = 33°), 10 kPa surcharge, displacement-
controlled drained push. Three legs in parallel, 3.2–4.1 h each.

1. **No limit point on any rung.** Every leg is still hardening where it ends:
   `dq/ds` decays to 15 % (corot) / 25 % (hypo) of its initial value and never
   approaches zero; `q` is monotone with its maximum at the last converged
   point. hypo reaches **5.31 × Vesic at the full s/B = 15 % target** and is
   still climbing (q/q_Vesic = 1.00 / 1.87 / 2.85 / 3.59 / 4.21 / 4.77 / 5.31
   at s/B = 1 / 2.5 / 5 / 7.5 / 10 / 12.5 / 15 %), with the tangent flat at
   6847 kPa/m over the final 200 steps.
2. **`hypo` is STIFFER than `corot`, and the gap grows** — +2.1 % / +6.2 % /
   +14.8 % of q at s/B = 1 / 2.5 / 5 %. *(These, and every other figure in this
   section, are the geometric content of a volumetrically LOCKED mesh — the
   whole ladder runs `-formulation std`. See the 2026-07-30 follow-up below,
   which measures a locking lever several times larger.)* Same sign and
   comparable magnitude to
   the P3 smoke (§5 P3: 0.3 / 1.8 / 15.4 % at 0.25 / 1 / 5 %), so the P3 result
   was not a small-model artifact. Large strain moves this model *away* from a
   limit point — the ADR-78 §1 hypothesis that missing large-deformation
   kinematics was suppressing the limit point is **refuted**.
3. **`-kozenyCarman` is inert on a drained push** (hypo vs hypo+kc agree to
   ~1e-4 over the whole backbone), which is the expected physics: k(J) only
   sets dissipation rate, and the drained backbone is k-independent. Its lever
   is undrained/rapid loading, which this campaign does not exercise.

No leg reached the s/B = 0.15 target, for two different reasons that
`RESULTS.md` keeps distinct: `corot` ended on its own at s/B = 0.060 when the
adaptive increment fell below its floor on a heavily distorted near-footing
mesh, while the two hypo legs were **stopped by the operator** at s/B ≈ 0.076
after 4 h, still converging comfortably. Neither ending is a limit load — the
two are distinguished by `dq/ds`, still clearly positive. So the campaign
cannot exclude a limit point *past* the last converged step; it does establish
that none appears in the span the three legs share (s/B ≤ 0.060).

Two things worth carrying forward. `corot` hit its convergence floor at
s/B = 0.060 while `hypo` was still stepping at 0.076, i.e. the rate-form UL
lane is better conditioned than rotate-only kinematics once elements are
genuinely distorted — a small independent point for
`-geom hypo`. And the over-hardening itself is now pointed at boundary
confinement rather than kinematics: the same footing at 0.5–1.5 B clearance ran
to tens of × Vesic as a pure oedometer, versus 2.7–3.6 × here at 4.5 B.
(Confinement stays a *measured* cause, but it is no longer the only one: the
2026-07-30 follow-up adds volumetric locking as a second and larger lever, and
the element probe below raises a third candidate in the Vesić anchoring.)

**Follow-up (same day):** a `-geom linear` BASE rung was added and the hypo
legs re-run to the full s/B = 0.15 target. Two things changed. Against base,
`corot` is 5.4 % SOFTER and `hypo` 7.5 % STIFFER — opposite signs — so the
14.8 % hypo-vs-corot gap is not all large-strain hardening; about half is corot
dropping below small strain. And at full depth `dq/ds` flattens to 6847 kPa/m
(−2.1 % across s/B 0.125→0.150, a 0.3 % band over the final 200 steps), which
is a plateau in the TANGENT and therefore the stronger claim that no limit
point is coming further along either.

A single-element deviatoric probe (`element_probe.py`) then settled the
material-level question and raised a bigger one. PDMY *does* reach perfect
plasticity — τ flat to 1.0000 of peak out to γ = 0.6 at every confinement,
τ_f/σ_v constant at 0.770–0.782 — so the absent limit point is a
boundary-value artifact, not a constitutive impossibility. But the mobilised
friction angle at failure is **50.4°, not the 33° supplied**, because the
failure surface is a Lode-independent cone calibrated in triaxial compression.
Vesić's N_γ is a Mohr–Coulomb formula: 35.2 at 33° against ~834 at 50.4°, a
factor of 23.7, which EXCEEDS the 16.3× over-strength measured against the
γ-term anchor. A large part of the reported over-strength may therefore be an
anchoring artifact. This is a hypothesis — 50.4° was measured in simple shear
while a bearing mechanism is nearer plane strain — and the rigorous check is a
limit analysis with the actual cone.

**Follow-up (2026-07-30) — the whole ladder was measured on a LOCKED mesh, and
locking is the larger lever.** Every rung above runs `-formulation std`, so the
ladder measures geometry on a volumetrically locked discretisation. `-geom hypo`
refuses bbar (`OPS_LadrunoUP.cpp`, P2: *"bbar in rate form is reserved"*), but
that guard sits inside the `METHOD_HYPO` branch and is the parser's only
`-formulation` guard — **corot composes with bbar unchanged** (ADR 78 §3.1; `Q`
is barred through the same `dNuBar_`). So corot is where the locking lever is
measurable today, with no new C++. Legs `corot_std` / `corot_bbar`, drained,
same mesh / material / push:

| s/B | corot `std` | corot `bbar` | **bbar/std** | hypo/corot | hypo/base | `std` rerun / archived |
|---|---|---|---|---|---|---|
| 1.0 % | 0.98 | 0.75 | **0.769** | 1.021 | 1.010 | **1.000** |
| 2.5 % | 1.76 | 1.16 | **0.657** | 1.062 | 1.033 | **1.000** |
| 5.0 % | 2.48 | — | — | 1.148 | — | **1.000** |

1. **Locking is ~5.5× the geometry lever.** At s/B = 2.5 % relieving it removes
   **34.3 %** of q, where the hypo-vs-corot rung adds 6.2 % and the ladder's
   *entire* content over small strain (hypo/base) is 3.3 %. The largest
   identified effect on this backbone is the element formulation, not the
   kinematics — and the lever was still growing (23 % → 34 %) where the data
   ends.
2. **The baseline is engine-matched, and that was not free.** `corot_std`
   re-runs `std` on the current `dist/bin` because the committed
   `backbone_corot.csv` came from a different build (a fresh worktree's `.pyd`
   predated ADR-78 and refused `-geom corot` outright — `LEDGER_quirks`). It
   reproduces the archive to **1.000 at all three checkpoints**, so the locking
   ratio is clean of engine drift rather than assumed to be.
3. **Relieving locking COSTS reach.** `corot_bbar` truncated at s/B = 0.0364 —
   earliest of every leg in the campaign (linear 0.0469, corot_std 0.0550, corot
   0.0602, hypo complete). Its early retry count was much better than `std`'s
   (3 at step 100 vs 23 at step 200) and then it walled. So the conditioning
   argument this section credits to the richer-kinematics lane does **not**
   generalise to the unlocked formulation, and the comparison span is capped at
   s/B ≤ 3.64 %.
4. **The interesting cell remains unmeasurable.** Locking and geometry are not
   additive — a locked mesh is over-stiff, so it under-develops the deformation
   geometry feeds on, and the honest reading of item 2 is that the geometric
   content of an *unlocked* mesh is unknown. It cannot be obtained in the
   rate-form lane at all while `hypo` refuses bbar. Lifting that is real work,
   not a flag: the barred volumetric operator must be rebuilt per step in the
   CURRENT configuration (`dNuBar_` is a precompute-once reference-config array),
   and the p-row's incremental Δln J with kinematic n(J) must be barred the same
   way or the discrete mass balance stops matching the volume change the solid
   sees, breaking the drained/undrained consistency P2 pinned. `Q` already
   follows `formulation_`; Δln J and n(J) do not. "Reserved" reads as principled,
   not lazy.

None of this disturbs item 2's refutation of the ADR-78 §1 hypothesis — missing
large-deformation kinematics still does not explain the over-hardening. It
demotes the *size* of that result: geometry is a few percent, locking is tens.
An undrained pair (`*_und`, scoping finding 9) is in flight; undrained loading
is where bbar's near-incompressibility lever should be largest, and those
numbers are not in this table.

Validation of the campaign model (all in `RESULTS.md`): global vertical
equilibrium exact to 0.0000 %, mesh volume exact, the closed-form surcharge
reaction correction exact, commanded ≡ measured settlement to machine
precision at every step, and Pardiso ≡ UmfPack to 5.2e-13.

---

## 9. Numerical collapse load with the actual cone (2026-07-30)

§8 ended by naming its own next step: *"This is a hypothesis — 50.4° was
measured in simple shear while a bearing mechanism is nearer plane strain — and
the rigorous check is a limit analysis with the actual cone."* That check has
been run. Full write-up, CSVs, figure:
`Ladruno_files/testbed/hypo_bearing/COLLAPSE.md`.

**Answer: the anchoring error is real but is ≈ 2.4×, not 52×, and the
over-strength is re-scaled rather than dissolved.**

### The instrument

PDMY cannot produce a collapse load — 20 nested surfaces approach the envelope
asymptotically, which is why 11 h of pushing ended still hardening. So the
material is replaced by an **elastic–perfectly-plastic `DruckerPrager` sharing
the measured failure surface** (ρ = √2·α = 0.344531), on the same mesh, the
same footing, the same surcharge and the same reaction correction; drained, so
`LadrunoUP` gives way to a 3-DOF `LadrunoBrick` carrying γ' as a body force.
A collapse load is a property of the failure surface and the flow rule, not of
the hardening path to them, so the substitution is legitimate — and it is the
only thing changed.

The surrogate is verified rather than assumed: a strain-controlled plateau
probe measures α = 0.24407 against the target 0.24362, and the +0.18 % is
exactly the σ_y = 0.2 kPa apex offset (predicted +0.000444, measured 0.00044).
Running `cone_probe.py`'s own stress-controlled probe on the surrogate returns
−0.01 %, so **that instrument is unbiased** — the 1.5° gap between PDMY's
measured txc-equivalent and its supplied φ = 33° is not a probe artifact.

### Three things wrong with the Chen–Han route

| basis | φ | N_q | N_γ | q_u (square) |
|---|---|---|---|---|
| nominal φ = 33° — what is quoted | 33.00 | 26.1 | 35.2 | **637.5** |
| cone, txc equivalent | 31.51 | 21.9 | 28.1 | 518.3 |
| cone, plane-strain equivalent (Chen–Han) | 53.72 | 673.1 | 1836.6 | **26 711** |
| **cone, plane-strain THEN Davis ψ = 0** | **38.87** | **55.0** | **90.3** | **1525.0** |

1. **The matching is unstable here.** tan φ = √(9α²/(1−12α²)) has a pole at
   α = 1/√12 = 0.28868 — above which no plane-strain Mohr–Coulomb equivalent
   exists — and the measured cone sits at **84.4 %** of it. A ±2 % move in α,
   inside `cone_probe`'s own 3.9 % Lode spread, swings q_u from 17.3 to
   43.8 MPa.
2. **It assumes associated flow and this set is non-dilatant** (`dil1 = dil2 =
   0`). Davis's ψ = 0 reduction (tan φ* = sin φ) takes 53.72° → 38.87° and N_q
   from 673 to 55, a factor of **12.2** by itself.
3. **Its coefficients are far outside their calibration range**, and the model
   says so: the measured sensitivity of the collapse load to cohesion is
   53.7 kPa per kPa of σ_y, where Vesić's `c·N_c·s_c` at 53.72° predicts 565 —
   **11× too large**.

### What was measured

Square footing, actual cone, non-associated, `-formulation bbar`, campaign mesh:
**q = 1140 kPa at the s/B = 10 % ultimate criterion** (q_max 1179.1 at the
full s/B = 0.15, where it plateaus at 0.6 % of its initial tangent; tail
extrapolation 1362 kPa; σ_y → 0 costs 0.92 %). That is
0.75–0.92 of the Davis anchor and **1/23 of the Chen–Han prediction**.

So the benchmark's PDMY backbone — 3384 kPa at s/B = 15 % — is **2.2× the
correct anchor and ~3.0× the measured collapse load of its own cone.**
Re-anchoring removes a factor of 2.4 from the 5.31× §8 reported. It does not
remove the anomaly; locking, large-strain embedment, boundary confinement and
mesh coarseness still have to carry the remaining ~2.2×.

### Validation, and two things it exposed

The plane-strain Prandtl–Reissner problem (weightless soil, surcharge q₀) has
`q_u = q₀·N_q` **exactly** — upper and lower bounds coincide. On a
tensor-graded plane-strain strip mesh at a mild cone (φ_ps = 27.47°), the
non-associated leg holds **139.2 kPa against the exact 138.9 (1.0020)** with
dq/ds = 0, out to the full s/B = 0.25 target. The machinery is correct.

- **Associated legs never plateau.** Every one of them, at every cone, hardens
  without reaching a limit point — `nq20_assoc` walls out 10–12 % *above* its
  own exact answer, and the associated SQUARE leg ends at s/B = 0.0139 still
  at 66 % of its initial tangent with 64.6 % of the mesh yielded and the
  plastic zone on both the sides and the base. Dilatancy at ψ = φ = 53.7° has
  to displace the surrounding soil and a fixed box resists it. The flow-rule
  bracket is therefore one-sided — but the non-associated rung is the one that
  matches PDMY's `dil1 = dil2 = 0`, so it is the rung the verdict needs.
- **A verification cone must be checked against its own initial state.** The
  first attempt (φ_txc = 20° at PDMY's moduli) sits at **m = 0.950 of yield
  under gravity alone**, because the elastic K₀ = ν/(1−ν) = 0.507 mobilises
  19.1° against a 20° cone. That leg measured its initial condition and is
  reported void; raising ν to 0.45 fixes it, legitimately, since a collapse
  load does not depend on elastic constants. At the real cone m = 0.579.

### Limits of this result

- **No refined collapse load, and the way refinement fails is diagnostic.** The
  non-associated series at the real cone terminates *earlier* at every
  refinement, in BOTH load cases independently — weightless Prandtl at
  s/B = 0.0124 / 0.0030 / 0.0014 with the termination tangent rising
  26 → 35 → 63 % of initial, and full (γ' + q₀) at s/B = 0.0520 / 0.0198 /
  0.0060 rising 19 → 26 → 69 %, at 4 / 8 / 16 elements across B. That
  is not a sequence converging on a collapse load. The CONTROL does the
  opposite — the mild cone plateaus at both 4 and 8 elements across B (1.0020
  and 0.9514 of the exact answer), so neither refinement nor the solver is the
  problem. A perfectly plastic
  **non-associated** frictional solid is past the Rudnicki–Rice localization
  threshold (the critical hardening modulus for banding is positive under
  non-associated flow), so the continuum problem has lost ellipticity, band
  width is set by the element size, and refinement resolves ever-thinner bands
  in a progressively worse-conditioned discrete problem. A mesh-converged
  number at this cone needs regularization (viscoplastic, gradient, Cosserat)
  or an explicit dynamic solve — not a finer mesh. The exact-oracle control
  shows 4 elements across B suffices at φ_ps = 27.5°; nothing here shows it
  suffices at 53.7°. This is the honest limit of the study.
- **Not a sharp limit point.** dq/ds decays as a power law in settlement
  (∝ s^−1.48, 41 432 → 302 kPa/m over s/B = 0.5 → 15 %), which is a punching
  signature — now reproduced by a *perfectly plastic* material, which rules out
  hardening as its cause.
- **Locking cannot be measured at collapse**, and that inverts §8's
  locking-leg item 3: on a perfectly plastic cone `std` is the formulation that
  cannot be pushed at all (walls out at s/B = 0.0031 even with an algorithm
  ladder), where on hardening PDMY it was `bbar` that cost reach. Over the
  shared span std/bbar runs 1.063 → 1.150 across s = 1 → 6 mm, still growing.
- **The plastic zone REACHES the lateral boundary at collapse.** At the full
  s/B = 0.15 the fully-mobilised zone (880 of 2816 elements) puts 16 of the 352
  elements in the outermost column at yield; the base stays clear. A
  14.5 B-clearance control mesh (`build_mesh_big.py`, 8064 hexes) puts the cost
  at about **3.4 %**: big/small runs 0.974 → 0.9686 → 0.9647 → 0.9643 → 0.9663
  over s = 5 → 110 mm, i.e. it BOTTOMS OUT and turns back up rather than
  diverging. It has not reached the full plateau (s/B = 5.6 % of 15 %), so this
  is a strong indication rather than a closed measurement — but the boundary is
  worth a few percent, not a factor, and it moves the quoted collapse load DOWN,
  widening the gap to PDMY's 3384 kPa rather than closing it. (The runner's first verdict said
  "contained"; that test compared the centroid to 90 % of the domain extent,
  and on a graded mesh the outermost hex is 3.1 m wide with its centroid at
  8.45 m of 10. Now tested by element-column membership.)

---

## 10. References

- Hughes, T.J.R. & Winget, J. (1980), "Finite rotation effects in numerical
  integration of rate constitutive equations arising in large-deformation
  analysis", IJNME 15.
- Flanagan, D.P. & Taylor, L.M. (1987), "An accurate numerical algorithm for
  stress integration with finite rotations", CMAME 62 (the unrotated-
  configuration formulation).
- Abaqus Theory Guide §1.4.3 (strain-increment integration), §1.5.3 (Jaumann
  vs Green–Naghdi, Table 1.5.3-1), §1.5.4 (state storage / ΔR co-rotation) —
  via the `abaqus-theory-foundations` skill.
- de Souza Neto, Perić & Owen (2008) §14 — the contrast case: the
  multiplicative/log-strain path (`-geom finite`), and why it needs a
  linear-elastic inner law in the fork's adaptor.
- ADR 09/10 (SolidTransformation seam, corot), ADR 78 (u–p corot, the
  incremental p-row findings that P2 inherits).
- Chen, W.F. & Han, D.J. (1988), *Plasticity for Structural Engineers* — the
  Drucker–Prager ↔ Mohr–Coulomb matching relations §9 uses and refutes.
- Davis, E.H. (1968), "Theories of plasticity and the failure of soil masses",
  in *Soil Mechanics: Selected Topics* — the non-associated reduction
  tan φ* = cos ψ sin φ / (1 − sin φ sin ψ), which is what the §9 anchor uses.
- Rudnicki, J.W. & Rice, J.R. (1975), "Conditions for the localization of
  deformation in pressure-sensitive dilatant materials", JMPS 23 — why a
  perfectly plastic NON-associated solid is already past the banding threshold,
  which is the diagnosis §9 gives for the failed refinement series.
- Prandtl (1920) / Reissner (1924) — the weightless surcharge-bearing problem
  `q_u = q₀·N_q`, whose coincident upper and lower bounds make it the exact
  oracle §9's validation leg is checked against.
