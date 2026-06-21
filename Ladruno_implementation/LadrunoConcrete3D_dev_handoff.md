---
title: "LadrunoConcrete3D — developer / C++-implementer handoff guide"
project: Ladruno
type: handoff guide
status: SHIPPED to `ladruno` — kernel (return map + analytic damaged tangent, g++-verified) + nDMaterial wrapper (classTag 33017) + ALL dimensional views (3D + PlaneStrain/AxiSymmetric/PlateFiber/PlaneStress, #299) + P3 Tier-2 IMPL-EX (oracle #301 → review-hardened #304 → C++ kernel port + `-implex` wrapper #309) + P3 Duvaut–Lions `-eta` (oracle #316 → C++ kernel port + `-eta` wrapper #318) + P2f cyclic `β_c` ORACLE + C++ kernel port (faithful CDPM2 compressive ductility, g++-verified, #321) + P2g monotone-`ω` no-heal cyclic damage (oracle + C++ kernel + wrapper, secant unload + SPD unload tangent, #325) + P2h `-ctTemper {none|alphat|proj}` compression→tension damage temper (oracle + C++ kernel + wrapper, #327) + P3 Tier-3 explicit (do_tangent-independence gate + CentralDifference softening demo, #328 — robustness trilogy complete). NEXT = multiaxial apportioning → P4 finite-strain. See §0 / §6b for the current-state handoff.
related:
  - "[[31_ladruno_concrete3d_adr]]"          # the ADR (decision record)
  - "[[project_ladruno_concrete3d]]"          # the agent-memory pointer
  - "[[10_ladruno_j2_plasticity]]"            # the kernel pattern + return-map IMPL-EX donor
  - "[[19_ladruno_rc_shell_adr]]"             # the shell/MCFT sibling (33015)
updated: 2026-06-20
---

# LadrunoConcrete3D — handoff to the C++ implementer

This is the working brief for whoever writes the **C++ return map**. The **numpy oracle**
(`tests/_testbed/concrete3d_ref.py`) is the *verified specification* — when in doubt, the oracle is
right and this doc explains it. Everything here has been **pinned to Grassl et al. 2013** (CDPM2,
IJSS 50:3805 / arXiv:1307.6998) by equation number and **adversarially reviewed twice** (a 15-agent
ADR scoping panel and a 5-agent final review). The physics is correct, and the C++ port + wrapper +
all dimensional views + Tier-2 IMPL-EX are **shipped** — **§0 is the current-state handoff** for the
next session; §1–§6 remain the verified technical reference; §7 lists what is left (Duvaut–Lions,
cyclic, Tier-3).

> **One-line summary:** CDPM2-grade solid concrete = effective-stress **Menétrey–Willam plasticity**
> (3-invariant surface, non-associated flow, confinement-aware ductility) + **dual scalar damage**.
> classTag **33017** (ND band). Solid/triaxial sibling of `LadrunoRCConcrete` (33015, shell/MCFT).
> **P0/P1 = plasticity only and is MONOTONIC (no peak); the peak/softening is the DAMAGE part (P2).**

---

## 0. Current state (2026-06-20) — handoff for the next session

Everything below is **shipped to `ladruno`** unless marked. Methodology throughout: **oracle-first**
(numpy `concrete3d_ref.py` = the verified spec) → **g++ kernel byte-check** (`concrete3d_kernel_check.cpp`
vs a regenerated `concrete3d_oracle_fixture.txt`, run by the pytest gate; g++ is available locally so the
check runs without the 30-min OpenSees build) → **wrapper** → **element battery** → one PR per increment
off `ladruno` (fast auto-merge ⇒ fresh branch each time; predict the next PR number for ledger refs).

**SHIPPED:**
- **Kernel + `nDMaterial LadrunoConcrete3D` (classTag 33017)** — MW surface, non-assoc flow, confinement
  ductility, semi-implicit return, dual damage ωt/ωc, crack-band, analytic damaged tangent. g++-verified
  (#240–#294).
- **ALL dimensional views (#299)** — 3D + PlaneStrain/AxiSymmetric/PlateFiber/PlaneStress via one
  `dim`-mode class (`setupDim`/`vmap`/`condense`); `getCopy(type)` selects the view (parser builds the 3D
  prototype). PlaneStress/PlateFiber = nested ε22-Newton + static condensation. Verified self-referentially
  vs the 3D material (no CDPM2 peer). `tests/test_ladrunoConcrete3D_element.py`.
- **P3 Tier-2 IMPL-EX (#301 oracle → #304 review → #309 C++/`-implex`)** — `-implex` reports an explicit
  extrapolated stress (plastic incr + dual damage frozen, `r=Δt/Δt_n` clamped [0,2]) + the degraded-elastic
  secant `D_dam(ω~):C0`; commits the exact implicit state. Kernel `State` += `wt/wc/dwt/dwc/depl[6]/dt_n`;
  `returnMap` has a `dt` param + Tier-2 branch; wrapper reads `dt=ops_Dt`. g++ check **B5** (`NIMPLEX`
  block) pins the reported stress to the oracle ~1e-16. Oracle gate `run_p3_implex_gate` / pytest
  `test_p3_implex_gate` (PI1–PI6).
- **P3 Duvaut–Lions `-eta` ORACLE (#316)** — viscoplastic relaxation at the PLASTIC level (ADR §4.4): the
  inviscid effective return + `κ_p` are relaxed toward the elastic trial by the Simo–Hughes closed form
  `σ̄ = (1−β)σ̄_tr + β σ̄_inviscid`, `β = Δt/(η+Δt)`; damage then follows from the **relaxed** effective
  stress, and the effective consistent tangent blends `C_eff ← (1−β)C_elastic + β C_inviscid`. `make_material`
  carries `eta`; `damaged_step_tensor`/`damaged_tangent_analytic`/`damaged_consistent_tangent` take an
  optional `dt`. **`η→0` (or `dt≤0`) ⇒ `β=1` ⇒ byte-identical to Tier-1** (a missing time increment falls
  back to inviscid, NOT to the elastic limit `β→0`). Oracle gate `run_p3_eta_gate` / pytest
  `test_p3_eta_gate` (PV1–PV6): PV1 `η=0` byte-exact; PV2/PV3 the **closed-form 1-D overstress oracle**
  (`duvaut_lions_1d_discrete/_analytic` — exact discrete steady overstress `= E·ε̇·η`, dt-independent, +
  an order-1 transient); PV4 the tensor kernel == the `(1−β)`-blend; PV5 the viscous damaged tangent (FD +
  the pre-onset blend identity on a genuinely-plastic confined-compression state); PV6 overstress-norm
  monotone in `η`.
- **P3 Duvaut–Lions `-eta` C++ KERNEL PORT + `-eta` WRAPPER (#318)** — the build PR for #316's oracle.
  Kernel `returnMap`: after the inviscid `returnMapTensor`, when `!implex && eta>0 && dt>0` form `sig_tr`
  via `elasticPredTensor` and blend `sig_eff ← (1−β)sig_tr + β sig_eff`, `kp ← (1−β)kp_n + β kp`, then
  blend `Dtan6 ← (1−β)C0 + β Dtan6` (so `damagedTangent` chains its damage linearization through the
  blended `C_eff`); `damagedUpdate` runs on the relaxed `sig_eff`; `sigEffImplicit` = the relaxed
  effective stress. Wrapper: parse `-eta η` (guard `η≥0`), member `eta`, set `p.eta`, pass `dt=ops_Dt`,
  serialize `eta` (`LC3D_NDATA`+1); warn `-eta`+`-implex` runs the IMPL-EX path inviscid. g++ byte-check
  **B6** (`NETA`): viscous nominal == oracle `damaged_step_tensor(...,dt)` (~1e-14), the `dt≤0` inviscid
  fallback == oracle byte-for-byte, + a non-tautology guard (max viscous−inviscid gap >1e-3). Element
  `test_eta_*` (LoadControl inviscid-limit + genuine-effect, db roundtrip). v1 = Tier-1 + `eta` only.

**TWO LIVE QUIRKS (also in [[LEDGER_quirks]]) — carry these forward:**
1. **Dual-damage IMPL-EX secant is SPD only on SINGLE-SIGN principal states.** `D_dam(ω~):C0` does not
   commute (two branch slopes 1−ωt≠1−ωc), so on a mixed-sign high-ω state (tension crack + lateral
   compression, ωt>~0.97) the symmetric part is INDEFINITE. Genuine dual-damage consistency and
   unconditional SPD are mutually exclusive — the contract is conditional (gate PI5 pins it).
2. **`-implex` + `DisplacementControl` past a limit point DIVERGES.** `dt=ops_Dt` is the load-factor (λ)
   increment under DisplacementControl — non-uniform near limit points → over-extrapolates ω. `-implex`
   is for implicit TRANSIENT dynamics + uniform LoadControl; for softening+DisplacementControl use Tier-1
   + an unsymmetric solver (or Tier-3 explicit). **A clean future fix: source the IMPL-EX `dt` from a
   monotone control parameter, not λ** — worth doing before promoting `-implex` for quasi-static softening.

**SHIPPED (P2f β_c — oracle + C++ kernel port, #321):**
- **P2f cyclic `β_c` (#321)** — the full CDPM2 `β_c` (Eq.50) restored into the compressive-damage plastic
  driver `κ_dc1` in BOTH the oracle AND the C++ kernel (`damagedUpdate` + `damagedTangent`); makes
  compression markedly more ductile (faithful CDPM2, user decision). g++-verified (the `dmg_compression`
  byte-check stress ~3.5e-15, its analytic damaged tangent with the `∂β_c/∂ε` term matches numerical
  ~8e-6). No wrapper/serialization change (β_c is computed from existing params, always-on). See §6b for
  the C1/C2/P2e re-gating, the `∂β_c/∂ε` analytic-tangent term (composite micro-FD through the return map
  — the kernel mirrors it), and the KEY finding that the cyclic story needs a separate monotone-`ω_c` fix.

**SHIPPED (P2g monotone-`ω` no-heal cyclic damage — oracle + C++ kernel + wrapper, #325):**
- **P2g monotone-`ω` (#325)** — `ω_t`/`ω_c` were re-solved each step against the LIVE effective drive
  stress; the inelastic histories `κd*` are already monotone, so the live drive was the ONLY non-monotone
  input — on elastic unload it drops with the histories frozen and the bracketed solve relaxed `ω` back
  (the crack HEALED, the #321 F4 diagnostic). Fix: drive each `ω` with the running MAX of its channel's
  drive stress (`sigt_max`/`sigc_max`, two new committed scalars in the kernel `State` AND the wrapper,
  serialized `LC3D_NDATA`+2), so `ω = ω(κ_d)` is monotone-nondecreasing, the unload follows the degraded
  secant `(1-ω)σ̄`, and the analytic damaged tangent drops the `−σ⊗∂ω` rank-update on an unloading channel
  (per-channel "advancing the max" flag) ⇒ the **SPD secant** `D_dam:C_eff` (vs the INDEFINITE loading
  tangent TD2). On any monotonic path `max == live` ⇒ **byte-identical** to pre-P2g (DT1/DT2/P2e/P2f and
  every monotonic g++ case unchanged). Oracle `run_p2g_gate` G1-G6 (tension no-heal + secant unload; reload
  retraces; compression no-heal; reversal monotone; reduce-to-monotonic; SPD unload tangent) + pytest
  `test_p2g_monotone_damage_gate` + g++ **discriminating** `dmg_cyclic_unload` byte-check (`nom_sig_err=0`;
  a live-drive kernel would heal and diverge) + element `test_cyclic_no_heal_unload`. The F4 diagnostic is
  now a real gate (`F4_ok`, no heal).

**SHIPPED (P2h compression→tension damage temper `-ctTemper` — oracle + C++ kernel + wrapper, #327):**
- **P2h `-ctTemper {none|alphat|proj}` (#327)** — the DT5 fix, behind a flag (user decision: offer BOTH
  temper modes). Literal CDPM2 (Eq.43, `κ̇_dt=ε̃̇`, no `(1−α_c)`) accumulates the TENSILE damage history in
  compression, so a compression excursion pre-damages a tension reload to ~0. A tensile weight `w_t` scales
  the `κ_dt1/κ_dt2` accumulation: `none` (default) `w_t=1` (literal CDPM2, byte-identical); `alphat`
  `w_t=1−α_c` (restores tension after compression, both monotonic backbones EXACT — the clean mode);
  `proj` `w_t = ‖P+ Δε_p‖/‖Δε_p‖`, `P+` = projection onto the POSITIVE effective-STRESS principal
  directions (tension restored; lightly softens the monotonic tension backbone). **KEY:** the plastic
  strain's OWN `<Δε_p>+` is NOT a valid shield — compression's dilatant flow makes the lateral plastic
  strains positive; project onto the tensile-STRESS frame. Analytic `∂w_t/∂ε` = `−∂α_c/∂ε` (alphat) /
  micro-FD (proj). Oracle `run_p2h_gate` H0-H4; g++ `dmg_cttemper_alphat` (biaxial, discriminating) +
  `dmg_cttemper_proj` (pure tension); wrapper parses `-ctTemper`, serializes the int (`LC3D_NDATA`+1).
  No new classTag/banner.

**SHIPPED (P3 Tier-3 explicit — validation/demo, #328):** the robustness trilogy is now complete
(Tier-1 implicit · Tier-2 IMPL-EX `-implex` · Duvaut–Lions `-eta` · **Tier-3 explicit**). No source
change was needed — the kernel already runs with `do_tangent=false` and the wrapper already exposes
`getRho()` for the explicit mass. Deliverables: (a) g++ **B7** gate — `returnMap(do_tangent=false)`
committed stress + state is **byte-identical** to `do_tangent=true` (`t3=0`; ADR §398 "Tier-3 committed
== Tier-1"), certifying the material is exact under an explicit solver that never factorizes the
indefinite softening tangent; (b) element `test_tier3_explicit_softening` — drives unconfined tension
into deep softening under `CentralDifference` (prescribed-displacement ramp, `system Diagonal`,
`-rho` mass), which **completes every step** where Tier-1 implicit `DisplacementControl` stalls at the
limit point, with the nominal stress peaking ~`ft` then degrading and `ω_t` developing.

**NEXT INCREMENTS (each its own oracle-first PR):**
- **Multiaxial-damage apportioning + plastic-dissipation regularization** — the remaining P2 refinements
  (extreme-principal vs `‖σ̄_t‖` norm; the D3/C3 un-regularized plastic dissipation caveat).
- **P4 finite-strain** — `nDMaterial LogStrain` wrapping the 3D material (clean: isotropic, no
  co-rotating backstress); already free via the wrapper, needs a validation gate.

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

**KERNEL DAMAGE increment (DONE — P3a, this PR, g++ byte-verified):** the P2 monotonic damage
subsystem ported into `LadrunoConcrete3DKernel.h`, mirroring the oracle verbatim:
- [x] **Damage kinematics** — `equivStrainGeneral` (Eq.37), `alphaCompression` (Eq.46),
  `damageDrivers` (`ε̃`, `α_c`, `x_s` Eq.56-57), `solveOmegaBracketed` (the bisection-safeguarded
  implicit-`ω` root, PR #261 no-spurious-healing); `As` added to `Params`.
- [x] **`damagedUpdate`** — the dual-damage NOMINAL stress `σ = (1−ω_t)⟨σ̄⟩₊ + (1−ω_c)⟨σ̄⟩₋` (Eq.1):
  re-eig the converged effective stress (⇒ **automatic unilateral**), accumulate the histories
  (clamped at `ε0`, full `‖Δε_p‖` via `plasticStrain6` isotropic compliance + tensor Frobenius), the
  **physical FLOOR** `> 1e-6·strength` on each `ω`-drive (review-fix), recompose. `State` extended with
  `sigEff` + the 6 damage-history fields; `returnMap` now returns the nominal `σ` and keeps
  `sigEffImplicit` = effective.
- [x] **g++ byte-check** — 4 `DMG` fixture cases (tension / reversal / confined-compression / shear
  committed damage states + a probe) reproduce the oracle `damaged_step_tensor` nominal stress to
  **~1e-14** (machine precision — single-step committed state, no multi-step amplification). Wired into
  the existing pytest g++ gate.
- [x] **DONE (P3b #291):** the ANALYTIC dual-projector DAMAGED tangent in the kernel — `returnMap`
  upgrades `Dtan6` to `D_dam:C_eff − σ̄_t⊗∂ω_t/∂ε − σ̄_c⊗∂ω_c/∂ε` (oracle `damaged_tangent_analytic`)
  on a converged effective return; self-verified vs a numerical central-diff of the same C++ stress.

**WRAPPER increment — DONE (#292 wrapper, #299 reduced views, #309 IMPL-EX; see §0):**

- [x] `#define ND_TAG_LadrunoConcrete3D 33017` in `SRC/classTags.h` (DEFINED, #292).
- [x] The `nDMaterial LadrunoConcrete3D` class (`.cpp/.h` + FEM_ObjectBroker + Tcl/Py parser + commit
  cycle + flat-`Vector` serialization), with:
  - `m0` enforced as `m0Of(fc,ft,e)` in the ctor (never user-set); `e∈(0.5,1]` guarded; full foot-gun
    guards (E>0, 0≤ν<0.5, 0<ft<fc, Gf,Gc>0, As≥1).
  - the unsymmetric-solver warning (and, under `-implex`, the conditional-SPD note).
  - **all dimensional views** (`dim`-mode, #299) and **`-implex` Tier-2** (#309) — see §0.

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

## 6b. P2 damage — status (SHIPPED through P2e; cyclic `β_c` = P2f, still owed — see §0/§7)

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

**P2e — the ANALYTIC dual-projector DAMAGED CONSISTENT TANGENT (the C++ deliverable) — DONE +
VERIFIED** (oracle `run_p2e_gate`, pytest `test_p2e_analytic_damaged_tangent_gate`, Zone-A 18/18),
FD-checked against the P2d numerical reference. Structure (ADR §4.3 [MAJOR]):
**`C = D_dam : C_eff − σ̄_t⊗∂ω_t/∂ε − σ̄_c⊗∂ω_c/∂ε`**:
- **`C_eff`** = the P1 EFFECTIVE consistent tangent `∂σ̄/∂ε` (numerical in the oracle; **ANALYTIC in
  the C++ kernel, #249** — so the C++ assembles its own `C_eff` and adds the damage linearization).
- **`D_dam`** = `isotropic_tangent(...)` — the spectral derivative of the per-principal damaged stress
  with `ω` FROZEN (**de Souza Neto Box A.6** operational form `dY:S = Σ_a y'_a(E_a:S)E_a +
  Σ_{a≠b}G_ab E_a S E_b`, `G_ab=(y_a−y_b)/(λ_a−λ_b)` → `y'_a` via l'Hôpital at coalescence — the SAME
  machinery as the P1 spectral tangent). This IS the `(1−ω_t)∂σ̄_t/∂ε + (1−ω_c)∂σ̄_c/∂ε` dual-projector
  secant.
- **`∂ω/∂ε`** via the **implicit-function theorem** on the bracketed ω-solve
  (`H = ∂F/∂ω = D[(1−ω)κ_d2/ε_f − 1]`), chained through the histories. The scalar sub-gradients
  `∂ε̃/∂σ̄`, `∂x_s/∂σ̄`, `∂α_c/∂σ̄` are **isolated scalar micro-FDs** (the LadrunoJ2/P1 "Lode directional
  gradient by micro-FD" pattern); `∂λ_extreme/∂σ̄` is the analytic **eigenprojection**
  (`E_max`/`E_min`, with the **tensor double-contraction weight `[1,1,1,2,2,2]`** on the
  Voigt off-diagonals — the one subtle bug found+fixed: the micro-FD sub-gradients are already
  per-Voigt-component, but the analytic eigenprojection needs the ×2 on shear); `∂‖Δε_p‖/∂ε` closed form.
- **Gates** (`damaged_tangent_analytic` == `damaged_consistent_tangent`, rel ~1e-10): **PE0** the
  spectral `dY/dX` (X², the damage function, near-degenerate l'Hôpital); **PE1** tension-damaged;
  **PE2** confined/triaxial compression-damaged (smooth); **PE3** shear + non-associated (`Df=0.3`,
  off-diagonal spectral terms); **PE4** load reversal; **PE5** reduces to elastic `C` / the P1 effective
  tangent pre-onset. **PE6 — the `σ̄_lat=0` Macaulay kink** (uniaxial-STRESS compression, the lateral
  eigenvalues at the tension/compression boundary): a **valid-subgradient** point, NOT a bug — the
  analytic tangent agrees with the central difference on the **loaded axial component** (rel ~3e-8) and
  differs only on the **~zero-stress lateral directions** by `O(ω_c−ω_t)` (the central diff crosses the
  kink). The C++ picks the committed-sign subgradient (same family as the eigenvalue-crossing the ADR
  flags).

**P2e ADVERSARIAL REVIEW (4-dim workflow, folded in — review-fixes PR):** the math core held (every
shipped uniaxial `σ11`/`ω_c` correct, no return-map/tangent error), three findings fixed: (1) **[MAJOR]
spurious `ω_t` in pure compression** — the `ω`-solve was gated only by `sig_t_drive>0`, so the ~1e-10 MPa
lateral-Newton residual (sign-dependent) drove `ω_t→1` in pure uniaxial-STRESS compression (mask-hidden
behind the compressive channel, but it poisons the tensile history). FIXED with a **physical floor**:
solve `ω_t` only when `sig_t_drive > 1e-6·ft` (mirror `ω_c`/`fc`), at all 3 sites (driver, single-step,
analytic-tangent). Gated: `DT0_pure_compression_wt < 1e-6`. (2) **[MAJOR] DT3 was tautological** — in its
`comp_early` window `ω_t≡0` (the driver floors it once all principals are compressive), so it never fed
`apply_damage_principal` a *(compressive principal, live `ω_t`)* pair; a broken unilateral routing still
passed. FIXED by testing the routing DIRECTLY + unconditionally (`DT0_unilateral` +
`DT0_compr_wt_invariant`: a compressive principal is carried by `(1−ω_c)` ONLY, byte-invariant to a live
`ω_t=0.95`); DT3 kept as the end-to-end check with an honest docstring. (3) **[MINOR] PE6** docstring
named the wrong kink directions (it's the lateral-normal AND coupled in-plane shear, not just lateral) —
reworded; the gate already only asserts the axial column.

**KNOWN LIMITATION (DT5 diagnostic, reported): compression→tension damage coupling.** Per literal CDPM2
**Eq.43** `κ̇_dt = ε̃̇` (FULL equivalent strain, **no `(1−α_c)` factor** — re-confirmed from the arXiv
source), so a compression excursion accumulates `κ_dt1/κ_dt2` and **pre-damages a subsequent tension
reload — today to ZERO tensile strength** (DT5 reports `tension-after-compression peak ≈ 0` vs fresh
`ft`). The monotonic tension/compression responses (DT1/DT2) are correct; this is the **cyclic T/C
coupling** (the dropped `β_c` Eq.50 + the open `α_t`-weighting question: literal-CDPM2 full-`ε̃` vs a
tensile-plastic-strain projection) and is **P2f** scope. Tracked, not gated.

**P2f `β_c` cyclic — ORACLE DONE (#321, `guppi/concrete3d-p2f-betac`).** The full CDPM2 `β_c` (Eq.50)
`= f_t·q_h2·√(2/3)/(ρ̄·√(1+2D_f²))` is restored into the compressive-damage plastic driver `κ_dc1`
(Eq.48) at all four sites (`drive_uniaxial_compression_damaged`, `drive_damaged_unified`,
`damaged_step_tensor`, `damaged_tangent_analytic`); `beta_c(sig_eff,kp,mp)` helper (ρ̄>0 guard, clamp
[0,1] — inactive in the damaging regime). In MONOTONIC compression `β_c ≈ f_t/(f_c√(1+2D_f²)) ≈ 0.058`,
so it makes compression **markedly MORE DUCTILE** than the `β_c=1` simplification (the chosen "faithful
CDPM2" direction — user decision 2026-06-20; post-peak stress shifts by **~23 MPa** vs `β_c=1`). The
analytic damaged tangent gains the `∂β_c/∂ε` term (a composite micro-FD THROUGH the return map, since
`β_c` depends on both `ρ̄(σ̄)` and `q_h2(κ_p)` — the C++ kernel will get `∂κ_p/∂ε` analytically from the
4-unknown return-map IFT). **Re-gates:** C1 holds (peak=`f_c`, softens, eff-monotone); C2 needed an
analytic exponential TAIL (`ε_fc·|σ_c,last|`) because the more-ductile response truncates the
by-construction `G_c` integral over a fixed strain (5.000 at all lch after the tail); P2e analytic==
numerical STILL holds with the new `∂β_c` term. **C++ port (same PR):** the kernel `damagedUpdate` +
`damagedTangent` mirror the `β_c` blend + the `∂β_c/∂ε` composite micro-FD (the kernel re-runs
`returnMapTensor` per component — `β_c` depends on both `ρ̄` and `q_h2(κ_p)`, and `∂κ_p/∂ε` isn't
otherwise exposed to `damagedTangent`); g++ `dmg_compression` byte-check stress ~3.5e-15, tangent ~8e-6.
Gate `run_p2f_gate`/pytest `test_p2f_gate` (F1 closed-form `β_c`; F2 monotonic backbone; F3 non-tautology
stress-gap ~23 MPa + `β_c∈(0,1)`; F4 a REPORTED diagnostic). **KEY FINDING — the cyclic story is bigger
than `β_c`:** the oracle (and kernel) solve `ω_c` IMPLICITLY against the CURRENT effective stress every
step, so an elastic UNLOAD lets `ω_c` relax back (HEALS); a cyclic-correct response needs `ω_c` driven by
the MONOTONE history (`ω_c ← max` over the path), a separate fix touching every driver + the committed
state + the C++ kernel.

**NEXT (P2f remaining, each its own slice):** (1) **monotone-`ω_c`** cyclic damage (no healing on
unload) — the F4 diagnostic; (2) the **compression→tension temper** (DT5; the `α_t`-weighting question) +
multiaxial-damage apportioning + plastic-dissipation regularization (D3/C3).

## 7. Roadmap context

P0 surface ✓ → P1 return-map/hardening/tangent ✓ → **C++ kernel return map + analytic tangent ✓** →
**P2 dual damage `ωt`/`ωc` + crack-band ✓** → **nDMaterial wrapper (33017) ✓** → **ALL dimensional
views ✓ (#299)** → **P3 robustness: Tier-2 IMPL-EX ✓ (oracle #301 → review #304 → C++/`-implex` #309;
freezes plastic state + damage)** → **Duvaut–Lions `-eta` ✓ (oracle #316 → C++ kernel + `-eta` wrapper
#318, the rate term, §0)** → **P2f cyclic `β_c` ✓ (oracle + C++ kernel port #321, faithful CDPM2 compressive ductility, §6b)** →
**P2g monotone-`ω` no-heal ✓ (oracle + C++ kernel + wrapper #325, secant unload + SPD unload tangent, §0)** →
**P2h `-ctTemper` compression→tension damage temper ✓ (oracle + C++ kernel + wrapper #327, none/alphat/proj, §0)** →
**P3 Tier-3 explicit ✓ (g++ B7 do_tangent-independence gate + element CentralDifference softening demo #328 — robustness trilogy complete)** →
**NEXT: multiaxial apportioning → P4 finite-strain (`LogStrain`, clean — already free via the
wrapper)** → P5 confined-fiber view (§4.6 hoop-spring condensation, "Mander by mechanism") → P6
auto-hybrid switch.

PRs (all → `ladruno`): **#240** P0+P1 · **#244** hardening · **#247** tangent · **#248** review ·
**#249** C++ kernel + g++ gate · **#259** P2a `G_f` · **#261** P2b `α_c`+`G_c` · **#284** P2c unified
tensor split + auto-unilateral · **#285** P2d single-step + numerical tangent · **#286** P2e analytic
damaged tangent · **#287** PE2 cross-platform · **#288** P2e review (ω floor) · **#289** P3a C++ damage
stress · **#290** guide · **#291** P3b C++ damaged tangent · **#292** nDMaterial wrapper (33017) ·
**#293** handout · **#294** wrapper convention tests · **#299** Phase-2 reduced views · **#301** P3
IMPL-EX oracle · **#304** IMPL-EX review fixes · **#309** IMPL-EX C++ port + `-implex` · **#316**
Duvaut–Lions `-eta` oracle · **#318** Duvaut–Lions `-eta` C++ kernel port + `-eta` wrapper · **#321**
P2f cyclic `β_c` oracle + C++ kernel port (faithful CDPM2 compressive ductility) · **#325** P2g
monotone-`ω` no-heal cyclic damage (oracle + C++ kernel + wrapper; secant unload + SPD unload tangent) ·
**#327** P2h `-ctTemper {none|alphat|proj}` compression→tension damage temper (oracle + C++ kernel + wrapper) ·
**#328** P3 Tier-3 explicit (g++ B7 do_tangent-independence gate + element CentralDifference softening demo).
