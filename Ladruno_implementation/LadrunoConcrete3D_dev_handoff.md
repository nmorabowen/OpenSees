---
title: "LadrunoConcrete3D — developer / C++-implementer handoff guide"
project: Ladruno
type: handoff guide
status: P0/P1 oracle + C++ KERNEL return map/tangent (g++-verified) + P2 dual damage DONE in oracle — P2a tensile ω_t, P2b compressive ω_c + α_c split, P2c UNIFIED TENSOR split + automatic unilateral crack-closure, P2d single-step tensor update + NUMERICAL damaged consistent tangent (the reference for the C++ analytic dual-projector tangent). NEXT = P2e (β_c cyclic + the ANALYTIC dual-projector tangent + multiaxial apportioning) → C++ returnMapDamaged port + the nDMaterial wrapper (carries the 33017 define + foot-gun guards).
related:
  - "[[31_ladruno_concrete3d_adr]]"          # the ADR (decision record)
  - "[[project_ladruno_concrete3d]]"          # the agent-memory pointer
  - "[[10_ladruno_j2_plasticity]]"            # the kernel pattern + return-map IMPL-EX donor
  - "[[19_ladruno_rc_shell_adr]]"             # the shell/MCFT sibling (33015)
updated: 2026-06-17
---

# LadrunoConcrete3D — handoff to the C++ implementer

This is the working brief for whoever writes the **C++ return map**. The **numpy oracle**
(`tests/_testbed/concrete3d_ref.py`) is the *verified specification* — when in doubt, the oracle is
right and this doc explains it. Everything here has been **pinned to Grassl et al. 2013** (CDPM2,
IJSS 50:3805 / arXiv:1307.6998) by equation number and **adversarially reviewed twice** (a 15-agent
ADR scoping panel and a 5-agent final review). The physics is correct; what remains is a faithful
C++ port plus the items in §7.

> **One-line summary:** CDPM2-grade solid concrete = effective-stress **Menétrey–Willam plasticity**
> (3-invariant surface, non-associated flow, confinement-aware ductility) + **dual scalar damage**.
> classTag **33017** (ND band). Solid/triaxial sibling of `LadrunoRCConcrete` (33015, shell/MCFT).
> **P0/P1 = plasticity only and is MONOTONIC (no peak); the peak/softening is the DAMAGE part (P2).**

---

## 1. Files & how they relate

| File | Role |
|---|---|
| `tests/_testbed/concrete3d_ref.py` | **The spec.** numpy oracle: surface, return maps, hardening, ductility, spectral tensor return, consistent tangent. Run it: `python tests/_testbed/concrete3d_ref.py` → P0/P1/HARDENING/TANGENT all PASS. |
| `tests/test_ladrunoConcrete3D_material.py` | pytest gates (11/11). Mirrors the oracle gates. |
| `tests/_testbed/concrete3d_kernel_check.cpp` | standalone **g++ self-check** of the kernel surface+hardening identities (`g++ -std=c++17 -I. … && ./c3dchk`). |
| `SRC/material/nD/LadrunoConcrete3DKernel.h` | **Header-only, OpenSees-free C++ kernel — return map + analytic tangent DONE.** Surface (`yieldF`), hardening (`qh1Of/qh2Of/ductilityXh`), `eccentricityFromKupfer`, `returnMapPrincipal`/`returnMapHardening` (analytic 3×3/4×4 Jacobians + apex + honest-f), spectral `returnMapTensor`, non-symmetric analytic `consistentTangent`, and the public `returnMap(strain, State)` all implemented + g++ oracle-verified. P2 damage / wrapper still stubbed. |
| `tests/_testbed/gen_concrete3d_fixture.py` + `concrete3d_oracle_fixture.txt` | oracle numeric-dump generator + committed fixture (physical driven paths + tangent cases) that `concrete3d_kernel_check.cpp` diffs against. |
| `Ladruno_implementation/31_ladruno_concrete3d_adr.md` | the decision record (read §4 for the formulation, §8 risk register). |

The kernel follows the `LadrunoJ2Kernel` "one core, many views" doctrine (plain doubles + `<cmath>`,
no OpenSees), so the same verified core serves the small-strain `nDMaterial` AND finite strain via
`LogStrain` (P4). **Tensor convention** (LadrunoJ2 lineage): symmetric tensors as `{00,11,22,01,12,02}`
with **true tensor** off-diagonals (J2 squares them directly); the engineering↔tensor conversion is
the wrapper's job at the OpenSees boundary. **Sign: compression NEGATIVE; `fc`, `ft` positive.**

---

## 2. The surface — VERIFIED CORRECT (Grassl Eq.12–21)

Invariants (`invariants()`): `σ̄_V = I₁/3` (Eq.12), `ρ̄ = √(2J₂)` (Eq.13), `θ̄ = ⅓ arccos((3√3/2)J₃/J₂^{3/2})`
(Eq.14). **Code stores `ξ = I₁/√3`, and `ξ/(√3·fc) = σ̄_V/fc` exactly** — that identity is why the
code looks slightly different from the paper but is byte-identical to it. `J₃ = det(deviator)`.

**Yield (Eq.18)** — `yieldF(sig, mp, qh1, qh2)`:
```
f_p = { [1−qh1]·(ρ̄/(√6 fc) + σ̄_V/fc)² + √(3/2)·ρ̄/fc }²
    + m0·qh1²·qh2·[ ρ̄/(√6 fc)·r(θ̄,e) + σ̄_V/fc ] − qh1²·qh2²
```
- The `[1−qh1]` term is the **ellipsoidal hardening cap** (closes the surface while `qh1<1`, vanishes
  at peak). **Lode `r` is on the friction term ONLY**, never the cap/quad.
- At `qh1=qh2=1` → **Eq.21** = `(3/2)ρ̄²/fc² + m0[ρ̄·r/(√6 fc)+σ̄_V/fc] − 1 = 0` = Menétrey–Willam.
- `m0 = 3(fc²−ft²)/(fc·ft)·e/(e+1)` (Eq.20). `r(θ̄,e)` = Willam–Warnke (Eq.19), `r(0)=1/e`, `r(π/3)=1`.
- **Eccentricity `e ∈ (0.5,1] EXCLUSIVE** — hard convexity precondition. Derive `e` from the Kupfer
  equibiaxial ratio `fcc/fc` (default 1.16 → `e≈0.52`). `eccentricityFromKupfer` is declared in the
  kernel but **undefined** — define it (bisection mirror of the oracle) in the build PR.

Verified machine-exact: uniaxial compression `(−fc,0,0)` and tension `(+ft,0,0)` on `f=0` to ~1e-16
**for any `e`/`m0`** (proves the fc-normalization), `m0` falls out of the tensile-meridian condition,
Kupfer recovers `e≈0.52`, `r` matches the canonical published form to 0.0, degenerates to von Mises
at `e=1`. C++ `yieldF`/`lodeR`/`m0Of`/`invariants` byte-match the oracle.

---

## 3. The return map — semi-implicit, two maps (Grassl Eq.22–36)

**Key idea (semi-implicit):** freeze `θ̄` at the trial value → the deviatoric return is **radial**
(direction = trial deviator), so the whole thing is a small Newton in **invariant `(ξ,ρ)` space** —
this **sidesteps the `∂θ̄/∂σ` Lode-corner singularity** entirely. After convergence, reconstruct:
`p_new = ξ/√3`, `s_new = (ρ/ρ_tr)·s_tr` (radial scaling preserves eigenvectors), `σ = p_new·I + s_new`.

There are **two maps** (load-bearing — decides which Jacobian you write analytically):

### 3a. Perfect-plastic map (`return_map_principal`) — 3 unknowns `(ξ,ρ,Δλ)`, ANALYTIC Jacobian
```
R1: ξ  = ξ_tr − 3K·Δλ·m_v
R2: ρ  = ρ_tr − 2G·Δλ·m_s
R3: f(ξ,ρ; θ̄_tr) = 0                       (failure surface Eq.21)
```
- `K = E/(3(1−2ν))`, `G = E/(2(1+ν))`. Factors `3K` (volumetric) and `2G` (deviatoric) are derived
  from `dε^p = Δλ·∂g/∂σ` with `Cᵉ = K I⊗I + 2G I_dev`.
- **Flow gradients (Eqs.22–25) — the potential is LODE-INDEPENDENT (no `r`!):**
  `m_v = ∂g/∂ξ = Df·m0/(√3 fc)` (volumetric, dilatancy-scaled), `m_s = ∂g/∂ρ = 3ρ/fc² + m0/(√6 fc)`.
  ⚠️ **Do NOT put `r` in the flow.** The *yield* gradient (Jacobian row 3) keeps `r`; the *flow*
  drops it. This asymmetry is correct CDPM2 — don't "fix" it.
- **Apex:** if the deviatoric return drives `ρ≤0`, return to the hydrostatic-tension vertex
  `ξ_apex = √3·fc/m0`, `ρ=0`. (3a's apex branch is at `return_map_principal`; the trigger is a
  transient `ρ<0` heuristic — see §6 known gaps.)
- The analytic 3×3 Jacobian is in the oracle (`return_map_principal`); FD-check your C++ version.

### 3b. Hardening map (`return_map_hardening`) — 4 unknowns `(ξ,ρ,Δλ,κp)`, oracle uses NUMERICAL Jacobian
Adds the hardening variable `κp` and uses the full Eq.18 surface (`qh1(κp)`, `qh2(κp)`):
```
R3: f_p(ξ,ρ,θ̄_tr; κp) = 0                               (Eq.18 hardened)
R4: κp = κp_n + Δλ·(‖m‖/xh(σ̄_V))·(2cosθ̄_tr)²            (Eq.32)
```
- `‖m‖ = √(m_v² + m_s²)` — **this is the correct Frobenius norm of `∂g/∂σ` precisely because the
  code uses `ξ=I₁/√3`** (HW orthonormal axial coord). If you ever switch to `σ̄_V=I₁/3` you must
  rescale by √3. The `(2cosθ̄)²` factor (Eq.32) = 4 on the tensile meridian, 1 on the compressive.
- **Hardening (Eqs.30–31):** `qh1Of(κp,qh0,Hp)`, `qh2Of(κp,Hp)` — already in the kernel.
- **Ductility `xh` (Eqs.33–36):** `ductilityXh(σ̄_V, fc, Ah,Bh,Ch,Dh)`, `Rh = −σ̄_V/fc − 1/3` (Eq.34).
  **`Rh` is the confinement term** — compression (`σ̄_V<0`) → `Rh>0` → larger `xh` → slower `κp` →
  more plastic strain to reach the failure surface = **confinement-dependent ductility by mechanism**
  (gated HD: strain at `κp=1` grows 0.0012→0.0032 from unconfined to `p/fc=0.1`). `xh` is C0 **and**
  C1 continuous at `Rh=0` (Eh, Fh enforce it).
- **The oracle uses a NUMERICAL 4×4 Jacobian** (FD, `du=1e-8`). **You owe the ANALYTIC 4×4 Jacobian**,
  including the piecewise `∂xh/∂σ̄_V` across `Rh=0` (Eq.33) and cubic `∂qh1/∂κp` (Eq.30). Add an oracle
  helper that returns the residual vector so you can FD-check your analytic Jacobian column-by-column.

**Reduce-to-P1 safety net (HA):** at `qh0=1, Hp=0` → `qh1=qh2=1` → the hardening map == the
perfect-plastic map to **2.6e-13** over compressive/deviatoric states. Use this as your first C++ gate.

**Why "monotonic, peak = damage":** with `Hp>0`, `qh2>1` for `κp>1`, so the effective-stress surface
keeps *growing* past `κp=1` — there is **no stress peak in plasticity alone**. The failure surface is
reached exactly at `κp=1` (`σ11 = fc` for uniaxial, gated HB). The peak and softening come from the
**P2 damage** part. Do not expect a peak from P0/P1.

---

## 4. The consistent tangent (P1-tangent gate)

The numerical consistent tangent (`consistent_tangent`) is the oracle reference; **you write the
analytic tangent and FD-check it against the oracle.** Verified facts you must reproduce:
- **Non-symmetric for non-associated flow** (`‖C−Cᵀ‖/‖C‖ ≈ 0.46` at `Df=0.3`). ⇒ **Tier-1 REQUIRES
  an unsymmetric solver** (`UmfPack`/`FullGeneral`) — unconditionally (see next).
- **Non-symmetric even in the associated limit** (`e=1,Df=1` → ≈0.024). The cause is the
  **frozen-eigenvector spectral recompose** (`σ=V·diag(σₚ)·Vᵀ` with `V` held fixed drops the
  `dV/dε` spin terms), **NOT** the Lode `θ`-freeze. Falsified (gate T3c): principal-space associated
  state is machine-symmetric (~2e-10); the asymmetry appears **only with shear** and scales
  **linearly** with it. *Implication:* don't try to symmetrize Tier-1 by removing the θ-freeze — the
  asymmetry lives in the eigenprojection terms.
- Elastic step → `C == C_elastic` (5.5e-13); reduces to principal for diagonal trials (0); frame-
  objective (2.6e-16); **quadratic-Taylor convergent** (ratio ~4 as δ halves — gate T4).
- The tangent gate runs on the **perfect-plastic map** (it has the analytic inner Jacobian). The
  hardening map's analytic tangent is yours to write + FD-check.

---

## 5. The C++ build-PR checklist — status

**KERNEL increment (DONE — this PR, `LadrunoConcrete3DKernel.h`, g++ oracle-numeric-dump verified):**

- [x] **Defined `eccentricityFromKupfer`** (inline bisection mirror of the oracle, + the
  `equibiaxialStrength` helper); the routine searches strictly inside the convexity band `(0.5, 1)`.
  *(Kept inline in the header — the LadrunoJ2Kernel "all-inline" doctrine — so the standalone g++
  self-check stays a single TU; the earlier "declaration-only in a `.cpp`" plan was superseded.)*
- [x] **`returnMap` implemented** (as `returnMapPrincipal` / `returnMapHardening` / `returnMapTensor`
  / public `returnMap(strain, State)`):
  - perfect-plastic 3-unknown Newton with **analytic 3×3 Jacobian** (matches the oracle to ~1e-13);
  - hardening 4-unknown Newton with **analytic 4×4 Jacobian** (incl. `∂xh/∂σ̄_V` across `Rh=0`,
    `∂qh1/∂κp`), self-FD-checked to ~1e-9 and oracle-matched to the ~1e-8 reference floor;
  - **hydrostatic-tension apex fallback** for both maps;
  - **honest convergence** — `f` recomputed at the returned stress with its *own* Lode angle;
  - **non-symmetric analytic `Dtan6`** — spectral lift (de Souza Neto Box A.6) of the principal
    Jacobian via IFT on the inner Newton residual; FD-matches the oracle numerical tangent to
    ~2e-10 (perfect-plastic) / ~7e-7 (hardening). The Lode directional gradient `dr/dw` uses an
    isolated scalar central-difference (corner-singular in closed form; the rest is closed-form).
  - `sigEffImplicit[6]` exposed (== `sigma` at P1; the LogStrain `bᵉ` fix, ADR R3).
- [x] **Committed g++ oracle-numeric-dump diff** — `concrete3d_kernel_check.cpp` (+ `gen_concrete3d_fixture.py`
  → `concrete3d_oracle_fixture.txt`); wired into `tests/test_ladrunoConcrete3D_material.py` (CI compiles
  + runs it where g++ is available). Floors: perfect-plastic stress 1e-9 / kp 1e-10 (analytic oracle
  Jacobian ⇒ machine-precision); hardening stress 1e-6 / kp 1e-7 (the oracle hardening reference is
  itself ~1e-8, numerical Jacobian, §6.3); tangent 1e-6.

**WRAPPER increment (DEFERRED — lands with P2 damage, where the stress peak first appears):**

- [ ] `#define ND_TAG_LadrunoConcrete3D 33017` in `SRC/classTags.h` (still deferred — **no orphan
  tag**; the LEDGER row reserves it, mirroring the LogStrain2D 33016 convention).
- [ ] The `nDMaterial LadrunoConcrete3D` class (`.cpp/.h` + FEM_ObjectBroker + Tcl/Py parser +
  commit cycle), and there:
  - **Enforce `mp.m0 == m0Of(fc,ft,e)`** in one place; never let the user set `m0` independently
    (`yieldF` trusts `mp.m0`); reject/clamp user-supplied `e≤0.5`.
  - **Require an unsymmetric solver unconditionally** (Tier-1) — document; warn on a symmetric solver.

---

## 6. Known gaps & caveats — state these honestly

> **PR #249 adversarial-review correction (2026-06-17):** a 5-lens panel found that in the C++ kernel
> the hardening apex branch could falsely report `converged=True` while teleporting a deep-COMPRESSION
> trial to the hydrostatic-TENSION vertex with `κp<κp_n` (even `κp<0`) — ~2% of converged hardening
> returns were inadmissible — directly contradicting the old item-1 claim below. **Fixed in the
> kernel:** every converged plastic hardening return is now gated on admissibility (`κp≥κp_n`,
> `Δλ≥0`, finite, on-surface); an inadmissible/non-converged state honestly reports
> `converged=false` and falls back to the elastic predictor (never poisons `κp`-history), so the
> caller cuts the step. The C++ deliberately diverges from the oracle's (equally-arbitrary) apex
> teleport here — it is the safe reference. A fuzz regression (`run_robustness` in
> `concrete3d_kernel_check.cpp`) asserts **0 inadmissible** converged returns over 180k trials. The
> *accurate* apex tangent (rank-deficient, ~0) is still owed — the kernel returns the elastic
> operator at a converged apex (≈K too stiff; flagged, slows but doesn't corrupt Newton). The panel
> also caught: the whole test file lacked the `zone_a` CI marker (so nothing ran in CI — fixed); the
> g++ check now adds per-entry + asymmetry-norm tangent checks and a NaN-guard regression.

1. **Hardening apex** is only a hydrostatic-vertex projection; the dedicated apex/Koiter sub-algorithm
   is owed. The honest-`f` flag + admissibility gate now genuinely report `converged=False` near the
   apex (see the PR #249 note above — pre-fix this *claim was false*, the apex lied `converged=True`).
2. **Off-meridian first-yield drift** (`κp~0`, small surface + steep `qh1` ramp): the semi-implicit
   return leaves O(0.1) off-surface drift that does **not** vanish under refinement — needs
   **sub-stepping near first yield off-meridian**. Diagnostic `HE` records it; only the compressive
   meridian is gated (the axisymmetric driver can't make off-meridian states).
3. **Hardening map = numerical Jacobian** in the oracle → its stress is ~1e-8; the analytic 4×4 is
   yours, with no oracle analytic reference (FD-check against the residual helper).
4. **Tangent is exercised on the perfect-plastic map only** (it has the analytic inner Jacobian) —
   *not* because the hardening FD tangent is noisy (it's equally step-stable).
5. **Flow is Lode-independent (simplified `mg`, no `r`)** per Eq.22–29 — CDPM2 design, not a shortcut.
   v1 uses a constant-`Df` volumetric flow; the full `mg` potential (Eq.23) is a follow-on (doesn't
   change peak strength — flow-independent).
6. **`qh1` can go non-monotone/negative for large `Hp`** (`−Hp·κp(κp−1)(κp−2)` term): Hp=2 → qh1 dips
   to 0.045; Hp=5 → −1.075. Default Hp=0.5 is safe; **warn if `Hp` would make `qh1` non-monotone**.
7. **Effective-stress plasticity is MONOTONIC; peak/softening = P2 damage** (not yet built).
8. **HB gate's 2% tolerance** is a step-discretization band (linspace resolution), not a physics
   tolerance — interpolate `σ11` at exactly `κp=1` to tighten if desired.

---

## 6b. P2 damage — status (started 2026-06-17, `guppi/concrete3d-p2-damage`, WIP)

CDPM2 damage pinned to **Grassl 2013 §2.3** by equation (fetched from the source):
`σ = (1−ω_t)σ̄_t + (1−ω_c)σ̄_c` (Eq.1); equivalent strain `ε̃` (Eq.37), uniaxial `ε̃ = σ̄_t/E`
(Eq.38); onset `ε₀ = f_t/E` ⟺ `q_h2 = ε̃/ε₀ = 1` ⟺ **κ_p = 1** (damage starts exactly at the
P1 failure surface — pre-peak pure plasticity, post-peak damage, no double-count); bilinear/
exponential softening (Eq.51-59), crack-band `ε_f → w_f/h`, `G_Ft = f_t·w_f/2`; `α_c` T/C split
(Eq.46); compressive exponential + confinement ductility (Eq.55-57).

**DONE + VERIFIED (oracle `run_p2_gate`, pytest `test_p2_damage_gate`, 13/13 Zone-A):**
- **D1 — the peak is damage:** the nominal uniaxial-tension stress **peaks exactly at `f_t`** then
  softens to ~0 (`ω_t→1`), while the P1 effective stress stays **monotonic** — proving the
  `(1−ω_t)σ̄` architecture and the onset-at-`κ_p=1` coupling. *(P1 alone has no peak.)*
- **D2 — the ADR §4.3 BLOCKING crack-band `G_f` energy gate PASSES:** dissipation·lch = `G_f`
  (0.1000 for lch ∈ {50,100,200}, rel err 2e-5) — **size-objective**.

The fix that unlocked D2 was the **exact CDPM2 inelastic-strain split** (fetched from the paper):
`ε_i = κ_dt1 + ω_t·κ_dt2` (Eq.52), with `κ̇_dt1 = ‖ε̇_p‖/x_s` (Eq.44, scaled accumulated PLASTIC
strain) and `κ_dt2 = (κ_dt−ε₀)/x_s` (Eq.45); `x_s = 1+(A_s−1)R_s`, `R_s = −√6 σ̄_V/ρ̄` for `σ̄_V≤0`
else 0 (Eq.56-57, `x_s=1` in uniaxial tension). `ω_t` is **implicit** (`ε_i` carries `ω_t`) → a 1-D
Newton; exponential softening `σ_t,nom = f_t·exp(−ε_i/ε_f)` (Eq.55), `ε_f = w_f/lch = G_f/(f_t·lch)`,
so `∫σ d ε_i = f_t·ε_f = G_f/lch` by construction. The lumped `ε_i ≈ ε_tot−σ̄/E` failed (twice)
because it dropped the `ω_t·κ_dt2` term and starved under tension's tiny ductility — **the precise
equations were essential, not a refinement** (the P0 `qh1·qh2`-guess lesson, reconfirmed).

**P2b — compression `ω_c` + the `α_c` T/C split — DONE + VERIFIED** (PR after #259; oracle
`run_p2b_gate`, pytest `test_p2b_compression_damage_gate`, Zone-A 14/14):
- **C0:** `α_c` (Eq.46) = 0 (uniaxial tension) / 1 (uniaxial compression); the **general equivalent
  strain Eq.37** `ε̃` = `ε₀` on the failure surface for *any* state (verified comp + tension) and
  reduces to Eq.38 in uniaxial tension.
- **C1:** nominal uniaxial-compression **peaks at `f_c`** (1.4% err) then softens, P1 effective
  stress monotonic.
- **C2:** crack-band `G_c` softening-law WIRING (`∫|σ|dε_i·lch = G_c/lch` BY CONSTRUCTION via
  `ε_i = κ_dc1 + ω_c·κ_dc2`, Eq.47-49,52,55, `ε_fc=G_c/(f_c·lch)`) — catches `ε_fc` errors, NOT an
  independent objectivity proof (see the PR #261 review note below).
- Ductility `x_s = 1+(A_s−1)R_s`, `R_s=−√6 σ̄_V/ρ̄` (Eq.56-57) wired (`x_s=A_s` in uniaxial
  compression — the confinement-ductility hook). **`β_c` (Eq.50) DROPPED for the MONOTONIC slice**
  (it is the smooth damage↔plasticity transition for CYCLIC loading) → restore as a P2c refinement.

> **PR #261 adversarial-review (workflow, 2026-06-17) — folded in:** (1) **CRITICAL** the implicit
> `ω_t`/`ω_c` Newton clamp-stalled to `ω=0` on steep-softening paths (small `ε_f`, e.g. large lch) →
> the cracked material spuriously **HEALED** (stress jumped back to the full effective stress); the
> residual `F(ω)` is non-monotone. **Fixed:** a bracketed (bisection-safeguarded) root solve
> `_solve_omega_bracketed` (the root is always bracketed on `[0,1]` during damage). Regression
> `test_p2_no_spurious_healing` drives the failing high-lch regimes. (2) **The `G_f`/`G_c` "objectivity
> gate" was TAUTOLOGICAL** — integrating `σ_nom` over `ε_i` (the softening law's own abscissa) returns
> `G_f/lch` by construction (invariant to E, ν, D_f, A_s). Relabelled D2/C2 honestly as the `ε_f`
> WIRING check. (3) **The FE-visible TOTAL fracture energy is NOT lch-objective (~33% spread)** because
> the effective-PLASTICITY dissipation is per-volume and un-regularized — a **CDPM2 damage-only-
> regularization characteristic** (CDPM2 regularizes the damage softening only). Now measured + REPORTED
> (D3/C3), not hidden; regularizing the plastic dissipation (ADR §4.3 [MAJOR]) is a documented P2c task.
> (4) onset-at-`κ_p=1` now directly gated (D0/C0b). (5) coverage gap: the C2 leg holds `x_s` constant
> (uniaxial `R_s=1`), so it does not exercise the confinement-ductility effect — a multi-confinement
> `G_c` leg is a P2c gap. **No correctness defect found in Eq.37/46/52/56-57 or the α_c split (the
> reviewers re-derived them clean).**

`κ_dt1`/`κ_dt2`/`κ_dc1`/`κ_dc2`, `α_c`, and Eq.37/38 are in the oracle.

**P2c — the UNIFIED TENSOR dual-damage update + unilateral crack-closure — DONE + VERIFIED**
(oracle `run_p2c_gate`, pytest `test_p2c_tensor_damage_gate`, Zone-A 16/16). P2a/P2b validated the two
scalar softening laws on *separate* uniaxial-stress drivers; P2c fuses them into the **single
constitutive update the C++ `nDMaterial` wrapper will call** — `spectral_split_principal`,
`apply_damage_principal`, `damaged_stress_tensor`, `drive_damaged_unified`:
- **One equivalent strain `ε̃` (Eq.37) + `α_c` (Eq.46) drive BOTH histories** — tension `κ_dt ← ε̃`
  (full, Eq.43); compression `κ_dc ← α_c·ε̃` (Eq.47); the plastic parts `κ_dt1=∫‖ε̇_p‖/x_s`,
  `κ_dc1=∫α_c·‖ε̇_p‖/x_s` (Eq.44/48, `β_c=1` monotonic); the damage-scaled parts `κ_d2=∫ dκ_d/x_s`
  (Eq.45/49). The **equations were re-pinned from the arXiv full text** (`ε̇_t=ε̇` Eq.43 vs
  `ε̇_c=α_c·ε̇` Eq.47, `κ̇_dt1=(1/x_s)‖ε̇_p‖` Eq.44 — tension has **no** `α_c` factor).
- **Nominal `σ = (1−ω_t)σ̄_t + (1−ω_c)σ̄_c` (Eq.1)** — `ω_t`/`ω_c` solved by the bracketed root,
  each driven by the **extreme effective principal of its sign** (`max⟨σ̄_i⟩₊` / `max⟨−σ̄_i⟩₊`), which
  collapses to the single axial stress P2a/P2b used → **DT1 reduces to P2a EXACTLY** (4.9e-13;
  confirming Eq.37 ≡ `σ̄/E` for *all* uniaxial-tension `σ̄`, not just on the surface), **DT2 to P2b**
  to the one onset-crossing step (2e-3 MPa). The genuinely multiaxial apportioning (biaxial/triaxial
  peak: extreme-principal vs `‖σ̄_t‖` norm; `/x_s` onset harmonization across channels) is a P2d gate.
- **UNILATERAL recovery is AUTOMATIC + tier-independent (ADR §4.3 BLOCKING)** — because the split is
  recomputed from the **converged effective stress every step**, a principal flipping negative (crack
  closing) is routed into the `(1−ω_c)` channel and is **no longer** multiplied by `(1−ω_t)`. **DT3:** a
  tension→compression reversal (`ω_t→1`) recovers `nominal/effective → 1.000` in early compression with
  **zero extra state** (no `s_rec·g_close` knob needed for the standard full-recovery assumption). The
  partial-recovery `s_rec` knob + the dual-projector analytic tangent + `β_c` cyclic stay **P2d**.
- **DT4** the damaged stress is frame-objective (6.4e-16) — `damaged_stress_tensor` eigendecomposes
  `σ̄`, splits the eigenvalues, recomposes on the SAME eigenvectors, so the spectral recompose carries
  the rotation and `ω_t`/`ω_c` (invariant) are unchanged.

**P2d — the SINGLE-STEP tensor update + the DAMAGED CONSISTENT TANGENT (oracle reference) — DONE +
VERIFIED** (oracle `run_p2d_gate`, pytest `test_p2d_damaged_tangent_gate`, Zone-A 17/17). The path
drivers carry state in locals; P2d factors the per-step logic into an explicit **`damaged_step_tensor`
(state, deps) → (σ_nom, new state, diagnostics)** — the exact contract the C++ `setTrialStrain` mirrors —
so one update can be **perturbed for the tangent** (fixed committed state, vary the increment) and
**chained step-by-step** (== the path driver). The committed `state` = `{σ̄, κ_p, ε, et_max, κ_dt1/2,
κ_dc/1/2}`; the plastic-strain increment uses the **tensor Frobenius norm** `‖Δε_p‖`, `ε_p = ε − C⁻¹:σ̄`
(generalizes the principal compliance the path drivers used; off-diagonals via the tensor elastic `C`).
- The **damaged consistent tangent** `dσ/dε = (1−ω_t)dσ̄_t/dε + (1−ω_c)dσ̄_c/dε − σ̄_t⊗∂ω_t/∂ε −
  σ̄_c⊗∂ω_c/∂ε` (the full dual-projector operator, ADR §4.3 [MAJOR]) is provided **NUMERICALLY** here
  (`damaged_consistent_tangent` = central difference of `damaged_step_tensor`) — the **honest reference
  the C++ ANALYTIC tangent is FD-checked against**, exactly the P1-tangent pattern (oracle numerical
  #247 → C++ analytic #249). The analytic spectral-projector derivatives `∂P_T/∂σ` + the implicit-ω
  chain rule are the **C++ deliverable**.
- **Gates:** **TD0** the single-step update reproduces the P2c path driver (tension 1e-14; compression
  ~5e-9, the eigendecompose-route floor over the long path); **TD1** with no damage the tangent reduces
  to the elastic `C` (7.8e-14) and, on a plastic-but-pre-peak step (`ω=0`), to the **P1 effective
  tangent** (3.3e-11) — both the `(1−ω)` factor and the `−σ̄⊗∂ω` rank-update vanish pre-onset; **TD2**
  on the softening branch the tangent is **DEGRADED + INDEFINITE** (`C[0,0]≈−3.2e3 < 0`,
  `λ_min(symC)≈−3.4e3 < 0`) — **the concrete Tier-2 IMPL-EX motivation**; **TD3** non-symmetric
  (`‖C−Cᵀ‖/‖C‖≈0.31`) for non-associated + damaged flow ⇒ unsymmetric solver (ADR §4.4); **TD4** the
  update is frame-objective (7.9e-16); **TD5** the tangent stays **finite** across a tension→compression
  **reversal** (`ω_t≈0.9`) and near an **eigenvalue crossing** (`σ̄≈` hydrostatic, where the analytic
  `∂P_T/∂σ` is singular — the C++ analytic tangent must regularize the crossing or accept Tier-2 drops).

**NEXT (P2e):** `β_c` cyclic (Eq.50) + the **ANALYTIC** dual-projector tangent (spectral `∂P_T/∂σ` +
implicit-ω chain rule, FD-checked against `damaged_consistent_tangent`) + multiaxial-damage apportioning
(extreme-principal vs `‖σ̄_t‖` norm; `/x_s` onset harmonization across channels) + plastic-dissipation
regularization (the D3/C3 caveat) → then the **C++ port** (`returnMapDamaged` over the P1 kernel:
spectral split + `_solve_omega_bracketed` per channel + nominal recompose + the analytic damaged
tangent) + the `nDMaterial` wrapper (lands classTag 33017 + the foot-gun guards).

## 7. Roadmap context

P0 surface ✓ → P1 return-map/hardening/tangent ✓ (oracle) → **C++ kernel return map + analytic
tangent ✓ (this PR, g++ oracle-dump verified)** → **P2 dual damage `ωt`/`ωc` + crack-band (where the
peak/softening comes from) + the nDMaterial wrapper (33017 define + foot-gun guards, §5)** → P3
robustness tiers (Tier-2 IMPL-EX freezes plastic state + damage; Tier-3 explicit) → P4 finite-strain
(`LogStrain`, clean — isotropic, no co-rotating backstress) → P5 confined-fiber view (§4.6 hoop-spring
condensation, "Mander by mechanism" in 1D fibers) → P6 auto-hybrid switch.

PRs so far (all → `ladruno`): **#240** P0 surface + P1 return map · **#244** hardening · **#247**
consistent tangent · **#248** review fixes · **#249** C++ kernel return map + analytic tangent +
g++ oracle-numeric-dump gate · **#259** P2a tensile damage + crack-band `G_f` · **#261** P2b
compressive damage + `α_c` split + crack-band `G_c` · **#284** P2c unified tensor dual-damage
split + automatic unilateral crack-closure · **(this PR)** P2d single-step tensor update + numerical
damaged consistent tangent.
