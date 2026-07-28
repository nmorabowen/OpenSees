# ADR 78 — `LadrunoUP -geom corot`: geometric nonlinearity on the Biot u–p family

**Status:** shipped (v1)
**Element:** `LadrunoUP` (ELE 33017, ADR 71)
**Seam reused:** `SolidTransformation` (`SRC/element/solidTransformation/`, ADR 09/10)
**Scope:** `-geom corot`, 3D lanes only (H8, Bézier Tet10). `-geom finite` stays reserved.

---

## 1. Why

The TIMs macroelement-calibration project (shallow footing on saturated sand,
`References/SFIM_model/`) cannot produce a bearing limit point: pushed to 5 % of
the footing width the model hardens monotonically to 9.3× the Vesic capacity.
One surviving cause is that no u–p element in the fork or upstream carries any
large-deformation support — the wedge **rotation** a bearing mechanism requires
is not representable. The drained branch has a validated workaround (a
total-stress `LadrunoBrick -geom corot` twin with buoyant weight, §6 gate ★);
the undrained/partially-drained branch — the entire point of the TIMg-UP
macroelement — has none. It needs geometric nonlinearity on a genuine u–p
element.

### 1.1 Why corot and not finite

`-geom finite` needs a `FiniteStrainNDMaterial` driven by `setTrialF(F)`. The
lifting adaptor, `LogStrainNDMaterial`, states its own limit: *"v1 assumes a
LINEAR-elastic inner law (so Cᵉ = inv(initial tangent) is constant and the εᵉ
recovery is exact)"* (`SRC/material/nD/LogStrainNDMaterial.h`). The target soil
law in every driving application is `PressureDependMultiYield` — pressure-
dependent hypoelastic, G ∝ p′^0.5 — which violates that assumption directly.
So finite-strain u–p would ship with **no usable soil constitutive law**. Corot
drives the ordinary `setTrialStrain()` path (ADR 10 §1: corot exists to
"deliver geometric nonlinearity while reusing the existing small-strain
material library") and works with PDMY01 unchanged.

### 1.2 The honest ceiling

Corot buys **large rotation and geometric (load) stiffness, not large
deformational strain**. A fully developed shear band with large local strains is
outside it; so is any effect driven by the change of the pore volume's
*magnitude* (porosity/permeability evolution, §4). What it does buy is the
difference between a mechanism that can rotate and one that cannot — which is
what a bearing failure is.

---

## 2. What already existed

- `SolidTransformation` seam (`METHOD_LINEAR|COROT|FINITE`) with
  `SolidTransformationCorot`: Procrustes/polar R from the nodal cloud,
  deformational displacement `u_d = Rᵀ(x−x_c) − (X−X_c)`, exact-gradient
  `globalizeForce`, `globalizeStiff = R̄ k R̄ᵀ + (G1+G1ᵀ)` with the v2.0
  higher-order deferrals (ADR 10 §5/§7). Stateless; serialized as a method id.
- `LadrunoBrick` as the consuming precedent, including the **COROT-1** idiom
  (global-frame dead loads stay OUT of the core force fed to `globalizeForce`
  and are added back unrotated) and the initial-stiffness R=I refresh (sweep #7).
- `LadrunoUP` (ADR 71) with the honest-p block contract:
  tangent `[K, −Q; 0, H]`, damp `[C_Ray, 0; Qᵀ, S+H̃]`, mass `[M, 0; 0, 0]`,
  residual `[∫Bᵀσ′dV − Q·p − body ; H·p − f_seep]`, solid-only Rayleigh.

The corot wrapper is 3D-only by construction: a planar node cloud gives a
rank-2 Procrustes cross-covariance, `det(H) = 0`, and the polar rejects it.
**`-geom corot` is therefore gated to ndm == 3 (H8, BTET10) at parse time.**
2D corot would need a dedicated planar-rotation wrapper — out of scope.

---

## 3. The u–p-specific design: per-block decisions

The element assembles four coupled blocks. Each gets an explicit frame
decision. Notation: `R` is the current polar rotation, `R̄ = blockdiag(R)` on
the solid DOFs, "core frame" = the de-rotated corotational frame in which the
small-strain material lives (fixed reference axes — ADR 10's kinematic-
hardening argument carries over unchanged).

### 3.1 K (solid stiffness) — the straightforward case

`update()` feeds the material `ε = B·u_d` with the **reference-configuration
B** and `u_d = localizeDisp(...)`, exactly as `LadrunoBrick`. `getTangentStiff`
builds the core-frame `K_s = ∫BᵀDB dV` and globalizes it through
`globalizeStiff(K_s, f_core, ·) = R̄ K_s R̄ᵀ + K_geo(f_core)`. B-bar composes
unchanged (it is a formulation choice inside the core frame; the F-bar
distinction belongs to `finite`).

### 3.2 Q (coupling) — the main correctness risk, resolved by construction

The risk named in the brief: Q carries a B-operator, and if the coupling and
the stiffness are evaluated in different frames the element is wrong in the
exact regime corot exists for. **Resolution: the −Q·p leg of the u-residual is
assembled in the CORE frame and routed through the same `globalizeForce` call
as the internal force.** I.e. the core force handed to seam 3 is

    f_core = ∫ Bᵀσ′ dV − Q·p        (both terms reference-B, core frame)

which is exactly the total-stress internal force `∫Bᵀ(σ′ − αp·m) dV` in the
corotated frame (p is a scalar and does not rotate). Consequences, all wanted:

1. **Frame consistency of K and Q is structural**, not a matching exercise —
   one force vector, one globalize.
2. **Pore pressure feeds the geometric stiffness**: `globalizeStiff` receives
   the same `f_core`, so K_geo is built from the *total* internal force — the
   physically correct prestress for a saturated medium (the effective-stress
   split does not apply to the load-stiffness).
3. The centroid back-reaction in `globalizeForce` stays zero: the Q·p force set
   is self-equilibrated per component (`Σ_a ∇N_a = 0` by partition of unity),
   same as `ΣBᵀσ′ = 0`.

The **tangent** u-p block becomes `−Q_rot`, `Q_rot = R̄·Q` (rows rotated,
columns untouched — pressure carriers do not rotate). This is the dominant
term of `∂[globalizeForce(f_core)]/∂p`; the Γᵀ·∂m/∂p sensitivity is dropped, in
line with the corot v2.0 tangent policy (residual exact, tangent
dominant-term; Newton converges to the right equilibrium regardless — ADR 10
§5). `Q_rot` is refreshed once per `update()` alongside R.

### 3.3 Qᵀ (storage-rate coupling) — the trap the design review caught

The mass-balance rate leg is `Qᵀ·u̇_d` (rate of skeleton volumetric strain).
The obvious implementation — damp p-row block `(R̄Q)ᵀ` contracted against the
integrator's nodal velocities — is **systematically wrong under rotation**,
and the failure analysis is the load-bearing result of this ADR:

*The chord defect.* The discrete velocity of a body rotating by Δθ per step is
dominated by the chord `(R_{n+1} − R_n)X/Δt`. Contracted in any frame F, its
apparent volumetric rate is `tr(Fᵀ(R_{n+1}−R_n))/Δt`; for F = R_{n+1} that is
`+2(1−cos Δθ)/Δt` — a **one-signed systematic apparent dilation** (an
expansion rate in the tension-positive convention, driving spurious undrained
*suction*; contracting in the committed frame F = R_n flips it to an apparent
compression — either way the bias never averages out), O(Δθ²) per step, and
*not removable by any velocity-linear operator*: the chord genuinely contains
a planar stretch defect that no instantaneous spin field can represent (spin
fields are traceless), so subtracting the best-fit rigid velocity (via the
seam's Γ) does not help, and a geodesic-midpoint frame kills the chord term
but picks up an equal-order defect from the integrator's velocity-history
terms. Cumulative over a rotation θ_tot in N steps: `Δε_v ≈ θ_tot²/N`.
Amplified by the storage modulus `Q̄ ≈ K_f/n` in undrained problems this is
**order-1 spurious pore pressure at bearing-mechanism rotations** (θ_tot ~
0.3, N ~ 300, Q̄ ~ 5·10⁶ kPa ⇒ ~10³ kPa of fictitious p). Unacceptable for the
driving application — and invisible in drained or small-rotation runs.

*The cure is incremental — and it must take the WHOLE p-row with it.* First
attempt: swap only the coupling to `Qᵀ·(u_d − u_d,committed)/Δt` (exactly zero
under rigid motion — `u_d ≡ 0` at every rigidly rotated configuration) while
leaving `S·ṗ` on the integrator's Newmark velocity. **Measured: unstable.** On
the consolidation column the corot run grew a ±100·q pore-pressure oscillation
at Δt = 0.02 where the linear run decays — mixing a backward-difference
coupling with Newmark-ṗ storage breaks the skew-symmetry of the discrete
coupling pair, and the broken pair pumps the (numerically damped) structural
ringing. The stable shape is the Book's canonical u–p pairing — GN22 on the
solid, GN11/backward difference on the **entire** p-row rate block:

    p-row (corot, transient):  H·p₊ − f_seep + Qᵀ·Δu_d/Δt + (S + αH̃)·Δp/Δt

Measured after the fix: the corot−linear consolidation gap converges cleanly
first-order (2.3e-2 → 1.3e-2 → 6.4e-3 → 4.2e-3 of the p scale over
Δt = 0.16/0.08/0.04/0.02), and rigid-motion exactness holds (`Δu_d = Δp = 0`).
This is also what production finite-deformation poromechanics codes do
(continuity discretized with the *incremental* volume change, cf. Abaqus
NLGEOM pore-fluid; Kratos U-Pw pairs Newmark displacements with a
backward-Euler pressure row). Cost: one committed nU-vector (`uDn_`, advanced
in `commitState`, rebuilt by `setDomain` from committed nodal displacements
after a restore — deliberately not serialized; committed p comes from the
nodes) and `Domain::getDT()` (falls back to the velocity form when no step
context exists, i.e. post-commit reporting, where `Δu_d = Δp = 0` anyway).

*Tangent policy.* `getDamp` keeps `(R̄Q)ᵀ` as the p-row coupling block — the
dominant Newton sensitivity of the incremental term (exact up to the
integrator's γ/β collocation factor on this one block, absorbed as an
inexact-tangent term under the corot v2.0 "residual exact, tangent
dominant-term" policy; the tangent's −R̄Q u-p block keeps its transpose
pairing). Under `-geom linear` none of this is active and the P1 velocity
form runs bit-identically.

*Reporting caveat (review finding 3).* `Domain::commit` zeroes `dT` before
recorders run, so a **post-commit** `getResistingForceIncInertia` (element
force recorders, `reactions -dynamic`) takes the `dt ≤ 0` fallback: the
reported p-row dynamic force uses the velocity form — converged velocities
are NOT zero, so at finite rotation rates the *reported* (never the solved)
p-row force carries the chord term. Solve accuracy is unaffected; recorded
element dynamic forces on rotating saturated elements should be read with
this in mind.

*Symmetrization caveat (review finding 6, tangent side of §3.5's damping
note).* The seam's `globalizeStiff` symmetrizes, so a non-associated soil
tangent (unsymmetric BᵀDB) assembles faithfully under `-geom linear` but as
`sym(R̄KR̄ᵀ + K_geo)` under corot — an inexact-tangent (convergence-rate)
difference between the two geom modes, inherited from the LadrunoBrick corot
lane and consistent with the v2.0 policy; recorded here for the u-p lane.

The rigid-rotation gate therefore asserts excess p at solver-tolerance level
(note the Q̄ amplification applies to Newton tolerance too: p_noise ~
Q̄·α·tol_disp/h — the test converges tightly and asserts orders below the
un-fixed failure mode, orders above roundoff).

### 3.4 H (Darcy) and S (storage) — left in the reference configuration, verified

Decision: **H, S, and the equal-order stabilization H̃ stay exactly as
assembled by `buildStaticBlocks()` (reference frame), unrotated.** This is not
just "second-order under small strain" — under **pure rotation it is exact**:

    H_spatial = ∫ (∇ₓN)ᵀ k_spatial (∇ₓN) dV,   ∇ₓN = R∇_X N,
    k_spatial = R k_mat Rᵀ  (permeability is a skeleton-attached material
                             property; its principal axes rotate with the body)
    ⇒ H_spatial = ∫ (∇_X N)ᵀ k_mat (∇_X N) dV = H_reference.

Same cancellation for H̃ (scalar α) and trivially for S (no gradients). The
anisotropic-k case is *why* the reference frame is the right one: keeping k on
the material axes is the physics, and rotating H while keeping k global-axis
would be wrong. Under deformational strain the current-config integrals differ
from reference by O(ε) — the corot premise, accepted and documented.

**The one place the frame does bite: the seepage drive.** `f_seep` and the
Darcy-flux response contract the **material-frame** gradient `∇_X N` against
the **global-frame** gravity/acceleration vector. The drive is therefore pulled
back:

    drive_mat = Rᵀ (b_f − ü_global)         (ü interpolated at the GP)

and the reported `darcyFlux` (response 4) is computed in the material frame and
pushed forward `q_global = R·q_mat` so recorders stay global-frame. Under
linear (R = I) both are bit-identical to the P1 code.

### 3.5 Rayleigh damping (solid-only contract)

`C_Ray = αM·M_s + βK·K_s + βK0·K_s0 + βKc·K_sc` is built in the core frame
(the element's own solid-K copies, per the ADR 71 contract) and then rotated as
one block: `C_glob = R̄ C R̄ᵀ` via `globalizeStiff(C, 0, ·)` (zero force ⇒ pure
rotation, the sweep-#7 idiom). The αM·M term is exactly invariant under this
(M's 3×3 nodal blocks are `m·I`), so the single call rotates precisely the
stiffness-proportional parts — the corotated K feeding βK·K_s **is** the same
object the contract intends, expressed in the frame the global velocities live
in. `K_geo` is deliberately NOT in the damping (damping is material, not
load-stiffness). `getRayleighDampingForces` uses the same rotated C.

*Caveat (review #5):* the seam's `globalizeStiff` symmetrizes its output, so
for a NON-associated material tangent (unsymmetric BᵀDB — the PDMY-class
regime) the corot damping is `sym(R̄CR̄ᵀ)`, i.e. the skew part of the βK term
is dropped, where the linear path keeps it. This enters the residual through
`damp·v`, so it is a (small, βK-scaled) physical damping difference, accepted
and recorded here rather than silently implied by the equality above.

### 3.6 M, initial stiffness, and the untouched rest

- `getMass`: `∫ρNᵀN dV` contains no gradients and no direction — rotation-
  invariant. Untouched.
- `getInitialStiff`: **stays the reference-configuration operator** `[K₀, −Q;
  0, H]` with R = I semantics — consistent with "initial" and with
  `LadrunoBrick`'s sweep-#7 decision (which refreshes the frame to cur == ref
  before globalizing K₀). We do not even need the refresh: the u–p element
  never lets `getInitialStiff` touch the live frame, so the corot state cannot
  go stale through it.
- `addInertiaLoadToUnbalance`: −M·(R·üg) with invariant M. Untouched.
- Body force (u-rows): a **global-frame dead load** — kept out of `f_core` and
  added after `globalizeForce`, the COROT-1 idiom verbatim.
- `stresses`/`stressesTotal` responses report the material's σ′ (and σ′−αp) in
  the **corotated material frame** — which is precisely what the rigid-rotation
  gate wants to see unchanged. Documented in the guide; the current-frame R is
  available for post-rotation if a global-frame recorder is ever needed
  (`frameTimeVarying` machinery, deferred until asked for).

### 3.7 Deferred by decision (not omission)

- **Porosity / permeability evolution with J** (n = 1 − (1−n₀)/J,
  Kozeny–Carman): requires tracking volumetric stretch — *outside the corot
  premise of small deformational strain*. It belongs to the `finite` phase and
  is deferred WITH that phase, blocked on the material-side adaptor (§1.1).
- **Higher-order corot tangent terms** (polar-Hessian, lever-arm/material-frame
  coupling, Γᵀ∂m/∂p): the ADR-10 v2.1 list plus the u–p coupling member; same
  rationale (exact residual, dominant-term tangent).
- **2D corot** (planar-rotation wrapper): no driving application; the polar
  rejects planar clouds today with a clear message.
- **Explicit/`-dynSeepage on` interaction**: the drive pull-back covers the
  residual; no additional corot terms are added to the (already FD-gated)
  dynamic-seepage path.

---

## 4. Sibling item, explicitly NOT this ADR: contact on ndf > 3 nodes

The contact handler guards `ndf == 3` at six sites
(`LadrunoContactHandler.cpp:808, 835, 973, 992, 1146, 1158`) and
`LadrunoContactFE` hard-sizes `3·(1+nps)`. What u–p models need is only the
**mechanical subset**: relax to `ndf >= 3`, fill the adapter's ID map from each
node's first three DOFs, pressure rides as a passenger. **This is NOT ADR 47
deferral 9** (gap flow, pressure penetration, thermal coupling — genuinely
coupled contact physics), which stays deferred; the ndf-relaxation is a
plumbing item to be sequenced independently. Recorded here so the deferral is
not misread as covering it.

---

## 5. Interface

    element LadrunoUP $tag $n1..$nk $matTag ... <-geom linear|corot>

- `-geom corot` **requires ndm == 3** (H8 or Bézier Tet10 lanes; the TH BTET10
  composes — corot acts on the 3 solid DOFs of every node, carriers and
  mid-edge alike). 2D is rejected with the planar-cloud reason.
- `-geom finite` is rejected with the material-blocker message (§1.1), so the
  error teaches the user *why*.
- Composes with `-formulation std|bbar`, `-pOrder`, `-stab`, `-lumped`,
  `-dynSeepage` unchanged.
- Serialization ships the method id (payload widened 29 → 30); `recvSelf`
  rebuilds via `SolidTransformation::create` (stateless, nothing else to ship).

**`-b` convention reminder (the silent factor-of-ρ trap):** `LadrunoUP -b` is
an **acceleration** (multiplied by the material's saturated ρ, like
`SSPbrickUP`); `LadrunoBrick -b` is a **force per unit volume**. Any
UP-vs-brick twin model must convert: `b_brick = ρ_sat · b_UP` (or buoyant
γ′ for the drained twin). The equivalence gate (§6 ★) carries a summed-reaction
weight check that catches a mix-up instantly; nothing else does.

---

## 6. Validation gates (tests/test_ladruno_up_element_corot.py)

1. **No-regression:** `-geom linear` is **bit-identical** to the pre-change
   element — the linear code path is untouched (every corot branch is behind
   `corotActive_`), and the existing ADR-71 battery (equiv/analytic/H8/T3/TH/
   PDMY/B5/init-recorders) runs on the default path unchanged. A direct test
   asserts `-geom linear` ≡ flag-omitted trajectories.
2. **Rigid-body rotation (the frame-mismatch catcher):** a saturated H8 block
   spun through large finite rotation with no deformation → σ′ unchanged in
   the material frame and excess p zero — an identity by §3.3's incremental
   coupling, asserted at solver-tolerance level (Q̄ amplifies Newton noise
   too), many orders below the un-fixed chord-defect failure mode. Both a
   static leg and a transient undrained leg (per-DOF Path series so every
   *converged* state is exactly rigid).
3. **Patch/objectivity under finite rotation:** the ADR-10 objectivity pattern
   on the u–p element — rotate + deform vs deform: material-frame σ′ agrees;
   drained and undrained (low-k, fast-load) limits.
4. **Consolidation under corot at small strain reproduces linear** — with the
   caveat that §3.3's incremental coupling is a (deliberately) different
   discretization of the same term, so corot vs linear trajectories differ at
   O(Δt) even at zero rotation. The gate is therefore refinement-aware:
   the corot−linear gap must shrink with Δt-halving and sit under a small
   absolute cap at the run Δt.
5. **★ Drained-equivalence gate:** with k large (fully drained) and hydrostatic
   p, `LadrunoUP -geom corot` must reproduce a **dry `LadrunoBrick -geom
   corot`** twin of the same mesh carrying buoyant weight
   (`∇·σ + γ_sat·ĝ = 0  ⇔  ∇·σ′ + γ′·ĝ = 0`). The equivalence is *measured*
   on the TIMs footing model (−0.33 % … −0.53 % over the whole backbone, base
   reaction = γ′V exact) — so the CI-sized twin gates at sub-percent, and a
   disagreement means the coupling blocks are wrong. Includes the summed-
   reaction weight check (§5).
6. **Seepage-drive frame:** a rotated block under gravity drive reaches the
   hydrostatic p field measured along the *rotated* depth axis (catches a
   missing Rᵀ pull-back, which gates 2–5 cannot see with their zero-gravity or
   unrotated-flow setups).

---

## 7. Files

| file | change |
|---|---|
| `SRC/element/ladrunoUP/LadrunoUP.h` | geom member + corot scratch + ctor arg |
| `SRC/element/ladrunoUP/LadrunoUP.cpp` | corot branches in update/tangent/resist/damp/flux; serialization +1 |
| `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp` | `-geom linear|corot` parse + ndm/finite guards |
| `Ladruno_implementation/71_ladruno_up_family_adr.md` | reserved-axis notes point here |
| `Ladruno_implementation/LadrunoUP_guide.md` | user-facing `-geom` section |
| `tests/test_ladruno_up_element_corot.py` | gates §6 |

No upstream file is touched; no new class tag (the seam is tag-free by design).
