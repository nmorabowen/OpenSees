---
title: "WP-144a — LadrunoNORSAND equation sheet (AB06/BA06 re-derived, curved CSL, WW ζ, Q-cap)"
project: Ladruno
type: equation sheet
status: "G2 owner decision applied 2026-10-01 (§16.6 item 14: exponential specific-volume update v = v0 exp(tr eps), dv/deps = v; every v-term of (S.26), (S.31)-(S.32), (S.40), (S.45)-(S.46), §1.4, §13, §14 updated and re-verified; scripts in Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_vexp/). G1 revision applied 2026-10-01 (12 items from the G1 oracle/test round, §16.6; owner decision: WW rho = 1/2 refused). G0 fix round 2026-09-30 (Adversary PASS-WITH-FIXES, 7 items; independent Opus numeric check PASSED; owner approved the refusal rule). Every derivative sympy/FD-checked; G0 scripts in the session scratchpad p0a/, G1 scripts in Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_r2/."
owner: nmora (Deriver: Fable, P0a)
related:
  - "[[144_ladruno_norsand_plan]] (design §2, oracles §5, roster §6)"
  - "[[134_sanisand_reference_integrator]] (the O1 oracle template)"
tags: [equation-sheet, norsand, critical-state, hyperelasticity, sand, wp-144]
updated: 2026-10-01
---

# WP-144a — LadrunoNORSAND equation sheet

This is the only document the O1/O2 oracle authors, the kernel author and the test author read. It is complete and
self-contained for the **small-strain** model in **principal space**. The finite-strain version (AB06 §3) is the same
algebra with the substitutions in §1.4; it is obtained in the fork by the LogStrain wrapper, not by new kernel code.

Sources: **[AB06]** Andrade & Borja 2006, IJNME 67:1531 (PDF pages given as AB06 p.N, N = 1531+page−1);
**[BA06]** Borja & Andrade 2006, CMAME 195:5115 (BA06 p.N). Tags: **[E]** read in the source (equation + page);
**[I]** my derivation or inference (one-line note); **sympy-checked** = verified symbolically or by finite differences
in `scratchpad/p0a/{invariants,zeta,energy,hardening,returnmap,cap,dissipation}.py` (all runs reported in §16.4).

Equation numbers `(S.n)` are this sheet's. `(AB06 nn)`, `(BA06 n.nn)` are the papers'.

---

## 1. Notation and conventions

### 1.1 Signs and invariants
- **Compression negative** (both papers and OpenSees). p < 0 throughout, q ≥ 0, π_i < 0 (image pressure), π_i* < 0.
- Principal stresses σ_a, a = 1,2,3 (AB06 write τ_a, the Kirchhoff stress; in small strain σ ≡ τ). δ_a := 1 for every a
  (AB06's "δ_a = 1"). Sums over repeated principal indices are written explicitly (`Σ_c`); no implicit summation.
- Invariants (AB06 8, p.1535; BA06 2.5, p.5117) [E]:

      p = (σ₁+σ₂+σ₃)/3,   ξ_a = σ_a − p,   R := ‖ξ‖ = (Σ_c ξ_c²)^{1/2}  (AB06's χ),
      q = √(3/2) R,       n̂_a = ξ_a/R,
      y := Σ_c ξ_c³ / R³ = cos3θ / √6,      θ = (1/3) arccos(√6 y) ∈ [0, π/3].                      (S.1)

  **Trap:** AB06 eq 8 defines y = tr ξ³/χ³ = (1/√6) cos3θ. y is *not* cos3θ. Eqs 19–20 are consistent with this y
  (sympy-checked). Everything below uses y as in (S.1).
- **θ origin** (AB06 p.1535) [E]: θ = 0 is the **tension corner** (triaxial extension, σ₁ = σ₂ < σ₃, one principal
  stress less compressive than the other two); θ = π/3 is the **compression corner** (triaxial compression, one
  principal stress more compressive). Check (sympy-checked): σ = (−100,−100,−400) → θ = π/3; σ = (−100,−100,+50) → θ = 0.
- Compression-corner form of the strength: on the yield surface q = −p η/ζ(θ); ζ(π/3) = 1 gives q = −pM at the image
  point in compression, ζ(0) = 1/ρ gives q = −pMρ in extension (so ρ = M_e/M_c, §15).

### 1.2 Strains
- Principal elastic strains ε^e_a; ε^e_v = Σ_c ε^e_c; e^e_a = ε^e_a − ε^e_v/3; ε^e_s = √(2/3) ‖e^e‖ (AB06 6, BA06 2.4) [E].
  These are **tensor** strains. Engineering shear γ = 2ε_ij enters only at the OpenSees shell boundary (Voigt
  conversion), never inside the kernel; in principal space there is no shear component.
- Total strain ε, plastic strain ε^p, ε = ε^e + ε^p (additive, small strain; BA06 Box 1) [E].
- Plastic strain-rate invariants: ε̇^p_v = tr ε̇^p, ε̇^p_s = √(2/3) ‖ε̇^p − (1/3)ε̇^p_v 1‖ (AB06 35) [E].
- Specific volume v = 1 + e. **Exponential update [I; G2 owner decision 2026-10-01, decision 1 = option b]:**

      v = v₀ exp(tr ε)   ⇔   v_{n+1} = v_n exp(tr Δε),        ∂v/∂ε_a = v  (every principal component; shears do not enter),

  with ε the **total** strain, v₀ the initial specific volume, and inside the local return (total strain fixed)
  ∂v/∂ε^e_a = 0 (sympy-checked, vexp_sympy.py (a)–(b)). The two forms are the same map (exp of a sum); the kernel
  uses the incremental one, so only tr Δε of the increment it receives matters. This **supersedes** BA06 Box 2
  step 6b, v = v₀(1 + tr ε) with ∂v/∂ε_a = v₀ [E], which was the G0/G1 rule: the two agree to first order
  (v₀eˣ − v₀(1 + x) = v₀x²/2 + O(x³), x = tr ε; d/dx of the difference is 0 at x = 0) and the exponential form is
  **exact under the LogStrain wrapper** (x = ln J, so v = v₀J, the AB06/BA06 finite-strain v = v₀ det F; §1.4). Why
  (G2 measurement, test_g2_logstrain.py): with the linear update the wrapper read v = v₀(1 + x) against the finite
  oracle's v₀eˣ, a gap v₀(1 + x − eˣ) (closed form to 1e−10), amplified ~25–30× into τ and π_i once the path dilates
  toward the critical state: 2.2e−3 / 2.3e−3 relative (τ / π_i) at 20 % drained TXC, 5.6e−3 / 6.0e−3 in TXE fork mode.
  On the TXC_paper drained path of §13 the gap at the path end is −4.8e−4 (tr ε = 0.0236, −v₀x²/2 = −4.7e−4;
  vexp_fd.py (3)).
- **Which v multiplies Π_v in the tangent [I; G2; FD-checked].** The explicit trial-strain derivative of the step is
  ∂v_{n+1}/∂ε̃_b = v_{n+1} (ε̃ = ε^e_n + Δε at fixed ε^e_n, vexp_sympy.py (b)): the **converged** specific volume of
  the step, not v_n and not v₀. Measured (vexp_fd.py (1), K2 paper set, drained TXC state 0.01 rad off the corner,
  central FD of the return map over ε̃ with v = v_n exp(tr Δε)): a^{ep} with v_{n+1} misses the FD by 2.4e−8 (h = 10⁻⁶;
  the O(h²) floor), with v_n by 1.7e−6, with v₀ by 1.8e−5 (v_n = 1.691, v₀ = 1.701); on a state built mid-path with
  v_n = 1.45 against v₀ = 1.701: 7.0e−9 / 3.8e−7 / 1.1e−4. The kernel's `return_map(…, v, vfac)` must therefore be
  called with vfac = v = v_{n+1}. DONE (G2, 2026-10-01): the kernel, O2 and O1 implement the exponential
  v-update and pass vfac = v_{n+1}; the kernel-vs-oracle parity was re-run.
- **v₀ stays a separate committed datum [I; G1, inverted at G2].** It no longer enters any derivative (the G1 note
  "dv/dε uses v₀, not v" is withdrawn), but it is still needed for the identity v = v₀ exp(tr ε) (a Gauss point whose
  strain history starts at ε = 0), for `revertToStart` (v := v₀) and for `initialState` (v = v₀). A state built
  mid-path (restart, a test fixture, a Gauss point initialised from a stress state) carries v and v₀ separately; the
  tangent is now insensitive to a wrong v₀, so the G1 FD trap is gone, but revertToStart would restart from the wrong
  specific volume. The G1 measurement (r2_neutral_v0.py: v-form 1.1e−4 vs v₀-form 7e−7 under the linear update) is
  kept in §16.4 as the record of the superseded rule.

### 1.3 Parameters (one table for the whole sheet)

| symbol | meaning | source |
|---|---|---|
| p₀, κ̂, ε^e_{v0}, μ₀, α₀ | BA06 energy: reference pressure (<0), elastic compressibility, reference volumetric strain, shear modulus, coupling | BA06 2.3 |
| M | critical stress ratio in compression (θ = π/3) | AB06 10 |
| N, N̄ | curvature of F and of Q on the meridian plane (0 ≤ N̄ ≤ N < 1) | AB06 10, 14 |
| ρ, ρ̄ | ellipticity of F and of Q (WW: ½ < ρ ≤ 1, **ρ = ½ refused**, §4.2; GA: 7/9 ≤ ρ ≤ 1). The same range applies to ρ̄ [I; G1] | AB06 11–12 |
| β := (1−N)/(1−N̄) ≤ 1 | volumetric non-associativity | AB06 p.1537, BA06 2.29 |
| χ (< 0, ≈ −3.5) | maximum-dilatancy coefficient, D* = χ ψ_i. **This is BA06's α (AB06's α); renamed to avoid AB06's χ = ‖ξ‖.** χ̄ := χ/β (AB06's ᾱ). | BA06 2.26, 2.29 |
| h | hardening constant (280 in both papers) | AB06 43 |
| λ̃, v_{c0} | "paper" CSL: v_c = v_{c0} − λ̃ ln(−p) | AB06 41 |
| e₀, λ_c, ξ, p_a | "fork" CSL: e_c = e₀ − λ_c (−p/p_a)^ξ (DM04 form) | plan §2.3 |
| c₁, c₂ | Q-cap blend bounds, η₁ = c₁M ≤ η₂ = c₂M (planar cap: c₁ = c₂; BA06's χ_cap = 0.10 → c₁ = c₂ = 0.10) | BA06 2.76, §10 |
| v₀ | initial specific volume (state input; committed datum for v = v₀ exp(tr ε), `initialState`, `revertToStart`; enters **no** derivative since G2, §1.2) | BA06 Box 2; §1.2 |

### 1.4 Finite-strain mapping (stated once; AB06 §3, BA06 §3.3) [E]
Replace ε^e_a by the principal elastic logarithmic stretches ε^e_a = ln λ^e_a, σ_a by the principal Kirchhoff stresses
τ_a, the trial elastic strain ε^{e,tr} by ε̃_a = ln λ̃_a from b^{e,tr} = f_{n+1} b^e_n f^t_{n+1}, and v = v₀ det F with
∂v/∂ε̃_a = v (BA06 p.5132). **Since G2 the v-mapping is identical to small strain, not a substitution [I; G2]:** the
LogStrain wrapper feeds the inner kernel ε_feed = ε_feed,n + (ε^{e,tr} − ε^e_n), whose increment has
tr Δε_feed = ½ ln det b^{e,tr} − ½ ln det b^e_n = ln det f_{n+1} = ln(J_{n+1}/J_n), so the kernel's own
v_{n+1} = v_n exp(tr Δε) is v₀ J exactly (sympy-checked, vexp_sympy.py (e); measured to 1e−10 at G2) and its
∂v/∂ε̃_a = v is the BA06 derivative. The local residual, Jacobian, nested π_i loop and ã^{ep}_ab are then
**identical** with no v₀ → v replacement anywhere (AB06 46–69; BA06 p.5130 "identical"). Only the assembly of the
spatial tangent differs (§9.5).

---

## 2. Elastic energy module

### 2.1 Interface (what the plastic part needs from the energy)
Ψ(ε^e_v, ε^e_s) → p = ∂Ψ/∂ε^e_v, q = ∂Ψ/∂ε^e_s, the 2×2 Hessian D = [[D₁₁, D₁₂],[D₂₁, D₂₂]] with
D₁₁ = ∂p/∂ε_v, D₁₂ = ∂p/∂ε_s = D₂₁ = ∂q/∂ε_v, D₂₂ = ∂q/∂ε_s (BA06 2.50; symmetric because Ψ exists) [E], and
the principal-space elastic tangent a^e_ab = ∂σ_a/∂ε^e_b. Principal stresses (AB06 63, BA06 3.41) [E]:

      σ_a = p δ_a + √(2/3) q n̂^e_a,     n̂^e_a = e^e_a/‖e^e‖   (= n̂_a of (S.1), co-axial, when q > 0).      (S.2)

General a^e_ab for **any** Ψ(ε_v, ε_s) (BA06 3.42 written in principal directions; sympy-checked against autodiff):

      a^e_ab = D₁₁ δ_aδ_b + √(2/3) D₁₂ (δ_a n̂_b + n̂_a δ_b) + (2/3) D₂₂ n̂_a n̂_b
             + (2q/(3ε^e_s)) (δ_ab − (1/3)δ_aδ_b − n̂_a n̂_b).                                              (S.3)

AB06 eq 64 writes a^e_ab = K δ_aδ_b + 2μ(δ_ab − δ_aδ_b/3) + √(2/3) d (δ_a n̂_b + n̂_a δ_b) with K = D₁₁, 3μ = D₂₂,
d = D₁₂. **That form drops the (2/3)(D₂₂ − q/ε_s) n̂_a n̂_b term and is exact only when q is linear in ε_s at fixed
ε_v** (true for BA06's energy, not for a general one) [I; sympy-checked: both forms agree for BA06's energy].
Limit ε_s → 0: replace q/ε_s by D₂₂ (valid for the BA06 energy; a general energy must supply its own limit) [I].

### 2.2 BA06 energy (BA06 2.2–2.3 p.5117; AB06 4–5 p.1534) [E]

      Ψ = Ψ̃(ε_v) + (3/2) μ^e(ε_v) ε_s²,   Ψ̃ = −p₀ κ̂ exp ω,   ω = −(ε_v − ε_{v0})/κ̂,   μ^e = μ₀ + (α₀/κ̂) Ψ̃ = μ₀ − α₀ p₀ e^ω.   (S.4)

Stresses and Hessian (BA06 2.51–2.52 p.5124; sympy-checked):

      p = p₀ e^ω [1 + (3α₀/(2κ̂)) ε_s²],          q = 3 (μ₀ − α₀ p₀ e^ω) ε_s,
      D₁₁ = −(p₀/κ̂) e^ω [1 + (3α₀/(2κ̂)) ε_s²] = −p/κ̂,   D₂₂ = 3μ₀ − 3α₀ p₀ e^ω,   D₁₂ = D₂₁ = (3 p₀ α₀ ε_s/κ̂) e^ω.   (S.5)

- Bulk modulus K = D₁₁ = −p/κ̂ ∝ p (exact, also for α₀ ≠ 0). Shear modulus μ^e = μ₀ − α₀p₀e^ω; for α₀ = 0 it is the
  constant μ₀ and the response decouples (D₁₂ = 0). Both papers' runs use α₀ = 0 [E, BA06 Table 1, AB06 §6.1].
- Convexity: det D = D₁₁D₂₂ − D₁₂². For α₀ = 0, det D = −3μ₀p₀e^ω/κ̂ > 0 always (sympy-checked); with α₀ ≠ 0 it
  must be tabulated (the plan's §2.5 gate).
- Isotropic compression (ε_s = 0): p(ε_v) = p₀ exp(−(ε_v − ε_{v0})/κ̂), i.e. ε_v = ε_{v0} − κ̂ ln(p/p₀) (K1.1).
- Conservative: W = ∮ σ:dε = ∮ dΨ = 0 on any closed elastic loop (K1.2).

### 2.3 Houlsby–Amorosi–Rojas slot (paper not yet in the library) [I]
An alternative Ψ_HAR(ε_v, ε_s) plugs in through §2.1 only. It must supply: p, q, the symmetric Hessian D, the
ε_s → 0 limit of q/ε_s, its region of positive-definiteness of D (and of the full 6×6 Hessian, which additionally needs
q/ε_s > 0), and a reference state where p = p₀. Nothing in §3–§12 depends on the energy except through (S.3) and D.
The dissipation argument (§11) does not use the energy at all.

---

## 3. Invariant derivatives in principal space

All entries sympy-checked against direct differentiation of (S.1).

First derivatives (AB06 19–20 p.1537):

      ∂p/∂σ_a = (1/3) δ_a,      ∂q/∂σ_a = √(3/2) n̂_a,
      y_a := ∂y/∂σ_a = 3ξ_a²/R³ − 3 (Σ_c ξ_c³) ξ_a/R⁵ − δ_a/R,                                          (S.6)
      θ_a := ∂θ/∂σ_a = −(2/√6) csc3θ · y_a.

Identities: Σ_c y_c = 0, Σ_c ξ_c y_c = 0, Σ_c n̂_c = 0, Σ_c n̂_c θ_c = 0, Σ_c σ_c θ_c = 0 (sympy-checked; used in §5, §11).

Second derivatives (AB06 27–28 p.1538):

      n̂_ab := ∂n̂_a/∂σ_b = (1/R) (δ_ab − (1/3) δ_aδ_b − n̂_a n̂_b),
      y_ab := ∂y_a/∂σ_b = 6 ξ_a δ_ab/R³ − 3 (Σ_c ξ_c³)/R⁵ (δ_ab − (1/3)δ_aδ_b − 5 ξ_aξ_b/R²)
                         + (δ_a ξ_b + ξ_a δ_b)/R³ − 9 (ξ_a ξ_b² + ξ_a² ξ_b)/R⁵        (no sum on a, b),      (S.7)
      θ_ab := ∂θ_a/∂σ_b = −(2/√6) csc3θ · y_ab − 3 cot3θ · θ_a θ_b.

### 3.1 The corners θ = 0, π/3 (csc3θ singular) and the safe evaluation [I; sympy-checked]
θ(σ) has a conical (|φ|-type) singularity at the corners: θ_a is bounded but discontinuous and θ_ab is unbounded there,
while every quantity the model needs (ζ, ∂ζ/∂σ, ∂²ζ/∂σ²) is regular. **Never form θ_a or θ_ab.** Use y as the third
variable and ζ = ζ̂(y):

      ∂ζ/∂σ_a = ζ_y y_a,    ∂²ζ/∂σ_a∂σ_b = ζ_yy y_a y_b + ζ_y y_ab,
      ζ_y := dζ/dy = −(2/√6) csc3θ · ζ'(θ),     ζ_yy := d²ζ/dy² = (2/3) csc²3θ [ζ''(θ) − 3 ζ'(θ) cot3θ].      (S.8)

(S.8) is an exact rewriting of AB06's θ-form (sympy-checked: q_a, q_ab identical). ζ_y is finite at both corners for
both shape functions; ζ_yy is finite at θ = 0, and at θ = π/3 it is finite for GA but **diverges like 1/|θ−π/3| for
WW** (§4.3) — harmless, because it multiplies y_a y_b = O(|θ−π/3|²). Corner branch (|sin3θ| < 10⁻⁸):

      ζ_y → −√6 ζ''(0)/9  at θ = 0,      ζ_y → +√6 ζ''(π/3)/9  at θ = π/3,      ζ_yy y_a y_b → 0 (set ζ_yy := 0).   (S.9)

The corner values are O(|θ−θ_c|) accurate at the switch, i.e. 10⁻⁸ relative. Also clamp the arccos argument √6 y to
[−1, 1] before taking θ (round-off can push it outside at the corners).

**Precision loss in (S.8) near the corners [I; measured].** In double precision the bracket ζ'' − 3ζ' cot3θ cancels
catastrophically: for ρ = 0.7, θ = 10⁻⁶ gives ζ_yy = −0.19392 (exact −0.19390), θ = 10⁻⁷ gives −0.1974, θ = 10⁻⁸ gives
−0.49 (the Adversary's evaluation gave +0.164): **all digits are lost for |θ − θ_c| ≲ 10⁻⁷.** This is harmless only
because ζ_yy enters solely through ζ_yy·q·y_a·y_b = O(|θ−θ_c|²) (measured product ≤ 2e−15 for θ ≤ 10⁻⁷). Two rules:
(i) sin3θ, cos3θ (and cot3θ) in (S.8) must be computed from the **same** θ used in ζ'(θ), ζ''(θ) — never from y and θ
independently; (ii) oracles and tests must **not** unit-test ζ_y or ζ_yy alone within |θ−θ_c| < 10⁻⁴ of a corner; test
q_a, q_ab, Ω, Ω_a instead (which are regular), or ζ_y at |θ−θ_c| ≥ 10⁻⁴ against the corner limits (S.9).

### 3.2 q = 0 (hydrostatic): the vertex rule [I; decided at G0, restated at G1]
R = 0 makes n̂, y, n̂_ab, y_ab undefined. On the axis the yield function loses its q-term: **F := pη(p, π_i)** (q dropped,
ζ does not enter), so its gradient is purely volumetric, **f_a = F_p/3** (every a), and F = 0 there means η = 0, i.e.
p = π_c. Q has a vertex on the axis: the deviatoric flow √(3/2) ζ̄ n̂_a + ζ̄_y q y_a has magnitude √(2/3)Ω →
√(2/3)·√(3/2)ζ̄ ≠ 0 along **every** non-axial approach, but no direction on the axis itself. **Rule (all modes, cap or
no cap): if R < R_tol at the trial state or at a local iterate, set n̂ := 0, y_a := 0, n̂_ab := 0, y_ab := 0, Ω := 0,
Ω_a := 0, so q_a = (1/3)βF_p δ_a (purely volumetric), ε^p_s = 0 and π_i stays at π_{i,n}.** Justification: (i) it is
the limit w → 0 of the cap (S.36), so the capped and uncapped models agree on the axis; (ii) it keeps (S.20) exact
(both sides zero) instead of the contradictory "no deviatoric flow but Ω = √(3/2)ζ̄"; (iii) the return from an axial
trial state beyond π_c is then the 1-D axial return to π_c with frozen π_i, which is what BA06's planar cap does. Ω is
therefore discontinuous at the axis (a vertex, not a smoothness bug).

**What the two integrators do there [I; G1, measured].** Without a cap, a near-isotropic compression path (deviatoric
to volumetric strain-increment ratio 0.2/3, the G1 "AMP_STOP" path) does not merely approach the axis: the plastic
flow at small η is compactive and its deviatoric part exceeds the imposed deviatoric strain rate, so the stress
reaches **q → 0 in finite time**. After that instant the continuum rate problem has no consistent solution — the
rate form needs a flow direction and the vertex has none, and with f purely volumetric the consistency condition
cannot hold for a non-isotropic ε̇: the deviatoric part of a^e:ε̇ is invisible to f = (F_p/3)δ yet moves the stress
off the axis, where F regains its ζq term. The rate oracle therefore **stops** there (status `vertex_reached`; measured at increment 13 of 40
on the AMP_STOP path, last state p = −236.7 kPa, q = 0). Backward Euler has no such problem: an axial (or
R < R_tol) trial state is returned **along the axis to p = π_c with π_i frozen** by the rule above, i.e. hydrostatic
compression beyond π_c is **perfectly plastic** (no volumetric hardening in this model, §10.1). Measured (O2, apex state
π_c = p = −100 kPa, three hydrostatic steps of tr Δε = −3e−3): p = −100.000000 after every step, π_i unchanged,
Δε^p_v = tr Δε exactly, ε^p_s = 0, D = −p·|tr Δε| = 0.300 per step. A non-axial trial state whose return crosses R_tol
is a corner problem for the local Newton and falls under the bounded local-work refusal (the O2 oracle refuses the
AMP_STOP step 12 (`refused_at` index) with `local_noconv` after 2⁸ substeps, re-measured at G2 under the exponential v-update of §1.2 (before it: step 13, `local_linesearch`); the kernel refuses at the same step with the same finest reason `LOCAL_NOCONV`, §9.1). K1.5 (flow-rule identity)
is tested only at plastic steps with q > 0, where (S.20) has no 0/0. K2 (no cap) starts isotropic at −100 kPa
**inside** the surface (elastic), so the rule never fires there. R_tol: 10⁻⁸·|p| (relative) [I; confirmed by the
oracle census: both oracles use it].

---

## 4. ζ(θ, ρ): shape functions

Both satisfy ζ(0) = 1/ρ (tension corner), ζ(π/3) = 1 (compression corner) and ζ'(0) = 0 for every ρ, and
ζ'(π/3) = 0 **on the admissible ranges** (GA: 7/9 ≤ ρ ≤ 1; WW: ½ < ρ ≤ 1 — not at ρ = ½, §4.2) (sympy-checked).
The deviatoric section is the polar curve r(θ) = 1/ζ(θ); it is convex iff r² + 2r'² − r r'' ≥ 0 on [0, π/3]
(checked numerically on a 4001-point grid; sign flips confirmed at the ranges quoted by AB06). **Every ζ range in
this section applies to ρ̄ (the ζ̄ of Q, §5.2) exactly as to ρ** [I; G1]: the parser checks both with the same rule.

### 4.1 Gudehus–Argyris (AB06 11 p.1535) [E] — option, refused for ρ < 7/9

      ζ = [(1+ρ) + (1−ρ) cos3θ]/(2ρ) = [(1+ρ) + √6 (1−ρ) y]/(2ρ),
      ζ' = −3(1−ρ) sin3θ/(2ρ),   ζ'' = −9(1−ρ) cos3θ/(2ρ),   ζ''(0) = −9(1−ρ)/(2ρ),   ζ''(π/3) = +9(1−ρ)/(2ρ),
      ζ_y = √6 (1−ρ)/(2ρ) = const,   ζ_yy = 0.                                                             (S.10)

Convex for 7/9 ≤ ρ ≤ 1 (numerically: min r²+2r'²−rr'' = −2.7e−3 at ρ = 7/9 − 10⁻³, +2.7e−3 at 7/9 + 10⁻³).
Linear in y ⇒ C^∞ in stress, no corner issue at all.

### 4.2 Willam–Warnke (AB06 12 p.1536) [E] — default

With c := cosθ, A := 4(1−ρ²), B := 2ρ − 1:

      ζ = [A c² + B²] / [2(1−ρ²) c + B (A c² + 5ρ² − 4ρ)^{1/2}].                                           (S.11)

This is 1/r_WW with r_WW the 1975 Willam–Warnke elliptic trace with r_t/r_c = ρ, θ from the tensile meridian [I].
Convex for ½ ≤ ρ ≤ 1 (numerically: min = −8.9e−15 at ρ = 0.5, +0.011 at 0.55, −214 at 0.45), but the **admissible
range is the open-ended (½, 1]: ρ = ½ exactly is refused, for ρ and for ρ̄** (owner decision 2026-10-01) [I; G1;
sympy-checked, r2_zeta_half.py]. At ρ = ½ the constant B = 2ρ − 1 vanishes and (S.11) collapses to ζ = 2cosθ: the
elliptic trace degenerates to a straight line (r cosθ = ½, the Rankine triangle), the convexity measure is identically
0, and **ζ'(π/3) = −√3**, not 0 — the compression corner is a **vertex** of the deviatoric section, so the
mirror-continued ζ has a kink in θ, ζ∘σ is only C⁰ there, and the corner branch (S.9), the C² argument of §4.3 and the
regularity of q_a all fail. For every ρ > ½ the square root in (S.11) equals 2ρ − 1 at θ = π/3 and ζ'(π/3) = 0
identically (sympy, with ρ = ½ + s, s > 0), but the corner curvature blows up as ρ → ½⁺:

      ζ''(π/3) = 3(1 − ρ²)/(2ρ − 1)²     (= 9.5625 at 0.7, 3 at 0.8; 209 at 0.55, 5.5e3 at 0.51, → ∞ at ½⁺),

so the corner constants of (S.9) and the FD floor of §4.3 degrade continuously towards ½; in double precision the
raw formula already loses ζ'(π/3) at ρ = ½ + 10⁻⁶ (evaluates to −5e−5 instead of 0). The refusal is at ½ exactly; a
value like 0.51 is admissible but a poor choice. Derivatives in θ:
ζ' = (dζ/dc)(−sinθ), ζ'' = (d²ζ/dc²) sin²θ − (dζ/dc) cosθ, with dζ/dc, d²ζ/dc² by the quotient rule (O2: generate
with sympy; the closed forms are long and add nothing). Corner constants (needed for (S.9); numeric, sympy-checked):

| ρ | ζ''(0) | ζ''(π/3) | ζ'''(π/3) | ζ_y(0) | ζ_yy(0) | ζ_y(π/3) |
|---|---|---|---|---|---|---|
| 0.7 | −1.1208791209 | 9.5625 | 174.94390 | 0.3050646566 | −0.1938959572 | 2.6025828517 |
| 0.71 | −1.0828693089 | 8.4336734694 | 137.80287 | 0.2947196961 | −0.1851525834 | 2.2953551841 |
| 0.8 | −0.75 | 3.0 | 20.784610 | 0.2041241452 | −0.1111111111 | 0.8164965809 |
| 1.0 | 0 | 0 | 0 | 0 | 0 | 0 |

ζ_yy(0) = (2/81) ζ''''(0) + (2/9) ζ''(0) (valid where ζ is even about the corner; sympy-checked at θ = 0).

### 4.3 Smoothness of WW at the compression corner [I; sympy-checked]
(S.11) is an analytic function of cosθ, but the physical locus on [0, π/3] is continued by mirror symmetry about
θ = π/3, and (S.11) is **not** mirror-symmetric there (ζ(π/3−d) ≠ ζ(π/3+d) from the raw formula). Consequently
ζ'''(π/3) ≠ 0 (table), the mirrored ζ is C² but not C³ in θ, ζ_yy ~ −ζ'''(π/3)/(27 |θ−π/3|) (checked: ζ_yy·δ →
−6.479 for ρ = 0.7, matches −ζ'''/27), and ζ∘σ is **C² in stress**. Effects:
- q_a is C¹, q_ab is continuous, the local Jacobian is continuous: Newton is unaffected.
- The consistent tangent is continuous with a |φ|-kink across the compression meridian. Central-difference tangent
  tests **at the exact TXC corner converge only O(h)** for symmetry-breaking perturbations (measured: rel. err 1.9e−2,
  1.6e−3, 1.6e−4, 1.6e−5 for h = 10⁻³…10⁻⁶) and O(h²) for axisymmetric ones (1.9e−11 at h = 10⁻⁶). The test author
  must run FD-tangent tests off the corners (or accept O(h) at the WW compression corner). GA has no such kink.
- At θ = 0 (tension corner) WW is even in θ and fully regular.

---

## 5. Yield function F and plastic potential Q

### 5.1 F (AB06 9–10 p.1535) [E]

      F(σ, π_i) = ζ(θ, ρ) q + p η(p, π_i),
      η = M [1 + ln(π_i/p)]                         (N = 0),
      η = (M/N) [1 − (1−N) (p/π_i)^{N/(1−N)}]       (N > 0).                                               (S.12)

Image point p = π_i: η = M. Compression-side apex on the axis (η = 0): π_c = π_i /(1−N)^{(1−N)/N} (AB06 p.1538), π_c = e·π_i
for N = 0. Tension apex p → 0: η → M/N (N > 0), → ∞ (N = 0). On the surface η ∈ [0, M/N] for p ∈ [π_c, 0].
Inverse (BA06 2.8): π_i/p = exp(η/M − 1) (N = 0); π_i/p = [(1−N)/(1 − ηN/M)]^{(1−N)/N} (N > 0).

Derivatives (AB06 24, 26, 51, 60; BA06 2.62–2.63; all sympy-checked, both branches):

      F_p := ∂F/∂p = η + p ∂η/∂p = (η − M)/(1−N)          [N = 0: M ln(π_i/p)]
      F_pp := ∂²F/∂p² = −(M/(1−N)) (1/p) (p/π_i)^{N/(1−N)}   [N = 0: −M/p]
      F_π := ∂F/∂π_i = M (p/π_i)^{1/(1−N)}                  [N = 0: M p/π_i]
      F_pπ := ∂²F/∂p∂π_i = (M/((1−N) p)) (p/π_i)^{1/(1−N)}  [N = 0: M/π_i]
      ∂F/∂q = ζ,   ∂F/∂θ = ζ' q   (θ-form),   ∂F/∂y = ζ_y q  (y-form).                                     (S.13)

Note ∂η/∂p = (F_p − η)/p and ∂η/∂π_i = F_π/p (used by the cap, §10). Gradient in principal space (AB06 60):

      f_a := ∂F/∂σ_a = (1/3) F_p δ_a + √(3/2) ζ n̂_a + ζ_y q y_a.                                          (S.14)

### 5.2 Q and the two readings of AB06 — what the kernel implements

AB06 eq 13–14 (p.1536) define Q = ζ(θ, ρ̄) q + p η̄(p, π̄_i) with η̄ of the form (S.12) with N̄, and BA06 2.11 fix the
free parameter π̄_i by Q = 0 at the current stress, i.e. η̄ = −ζ̄q/p. On the yield surface that gives **η̄ = (ζ̄/ζ) η**
(AB06 p.1538 line before eq 31) — call this **reading B**. But AB06 eq 22 states ∂Q/∂p = β ∂F/∂p, and every downstream
result (eqs 23, 25, 37–42, 51, 54, Box 2) uses it — call this **reading A**. The two agree iff ζ̄ = ζ (ρ̄ = ρ, in
particular the 2-invariant BA06 where BA06 p.5119 says η̄ = η directly). For ρ̄ ≠ ρ:

      reading A:  ∂Q/∂p = β (η − M)/(1−N) = (η − M)/(1−N̄),
      reading B:  ∂Q/∂p = (η̄ − M)/(1−N̄) = ((ζ̄/ζ)η − M)/(1−N̄).                                            (S.15)

**Suspected paper inconsistency** (§16.2): the flow rule the paper implements (23) is reading A while its dissipation
proof (30–32) is reading B. **This sheet adopts reading A** (the paper's algorithm, its tangent, and the K2 benchmark
were produced with it). Reading A has a genuine potential [I; sympy-checked that its gradient is (S.16)]:

      Q_A(σ, π_i) := ζ(θ, ρ̄) q + β p η(p, π_i),      i.e. "π̄_i is fixed by η̄(p, π̄_i) = η(p, π_i)"  (BA06 p.5119).

Equivalently, the plastic potential is (S.12) with N → N̄, ρ → ρ̄ and π̄_i/p = [(1−N̄)/(1 − ηN̄/M)]^{(1−N̄)/N̄}, so that
F and Q_A share the same stress ratio η but not the same Lode scaling. §11 gives the dissipation condition for reading A.

Derivatives of Q_A (AB06 22, 26; sympy-checked):

      ∂Q/∂p = β F_p,   ∂²Q/∂p² = β F_pp,   ∂²Q/∂p∂π_i = β F_pπ,   ∂Q/∂q = ζ̄,   ∂Q/∂y = ζ̄_y q,
      ζ̄ := ζ(θ, ρ̄), ζ̄_y, ζ̄_yy as in §4 with ρ̄.                                                          (S.16)

### 5.3 Flow vector and its Hessian (AB06 16–17, 21, 23, 25; sympy-checked)

      ε̇^p = λ̇ Σ_a q_a m^a,   m^a = n^a ⊗ n^a (eigenprojections of σ, co-axial with ε^e by isotropy),
      q_a := ∂Q/∂σ_a = (1/3) β F_p δ_a + √(3/2) ζ̄ n̂_a + ζ̄_y q y_a,                                        (S.17)
      q_ab := ∂q_a/∂σ_b = (1/9) β F_pp δ_aδ_b + √(3/2) ζ̄ n̂_ab + ζ̄_y q y_ab + ζ̄_yy q y_a y_b
                        + √(3/2) ζ̄_y (n̂_a y_b + y_a n̂_b)          (symmetric),                            (S.18)
      q_{a,π} := ∂q_a/∂π_i = (1/3) β F_pπ δ_a.                                                              (S.19)

Plastic strain-rate invariants (AB06 36–38; sympy-checked):

      ε̇^p_v = λ̇ Σ_c q_c = λ̇ β F_p,      ε̇^p_s = λ̇ √(2/3) Ω,
      Ω := [ (3/2) ζ̄² + (ζ̄_y q)² S_y ]^{1/2},   S_y := Σ_c y_c²     (θ-form: (ζ̄'q)² Σ_c θ_c²; 2-invariant: Ω = √(3/2)),   (S.20)
      Ω_a := ∂Ω/∂σ_a = (1/Ω) [ (3/2) ζ̄ ζ̄_y y_a + ζ̄_y q (ζ̄_yy q y_a + √(3/2) ζ̄_y n̂_a) S_y + (ζ̄_y q)² Σ_c y_c y_ca ].   (S.21)

(S.21) is AB06 eq 55 in the y-form (sympy-checked against direct differentiation; the paper's eq 55 has implicit sums on c).
Ω depends on σ only (not on π_i) — until the cap of §10 is switched on.

Dilatancy (AB06 39) [E]: D := ε̇^p_v/ε̇^p_s = √(3/2) β F_p/Ω = √(3/2) (η − M)/((1−N̄) Ω). 2-invariant: D = (η−M)/(1−N̄)
(BA06 2.25). D < 0 (compaction) for η < M, D > 0 (dilation) for η > M, D = 0 at the image point.

---

## 6. Critical state line: two modes

The CSL enters the model **only** through ψ_i (hence π_i*). Define Λ(π_i) := π_i ∂ψ_i/∂π_i.

**Paper mode** (AB06 41 p.1540; BA06 2.23) [E]:  ψ_i = v − v_{c0} + λ̃ ln(−π_i),  ∂ψ_i/∂π_i = λ̃/π_i,  Λ = λ̃.
(v_{c0} is the critical specific volume at unit pressure; ψ_i < 0 dense.)

**Fork mode** (plan §2.3; sympy-checked):

      e_c(π) = e₀ − λ_c (−π/p_a)^ξ,    ψ_i = e − e_c(π_i) = (v − 1) − e₀ + λ_c (−π_i/p_a)^ξ,
      ∂ψ_i/∂π_i = λ_c ξ (−π_i/p_a)^ξ / π_i,     Λ(π_i) = λ_c ξ (−π_i/p_a)^ξ.                               (S.22)

In both modes ∂ψ_i/∂v = 1. Every formula below is written with Λ(π_i); AB06's λ̃ is the paper-mode value of Λ.
The fork Λ → 0 as π_i → 0 (the reason for the change: no log divergence), and ψ_i is finite at π_i = 0.

---

## 7. D* = χ ψ_i and the limit image pressure π_i*

D* = χ ψ_i (BA06 2.26; AB06 40 writes D* = α ψ_i) [E]. Setting D = D* in §5.3 gives the limit stress ratio and, through
the inverse of (S.12), π_i* (AB06 40, 42; BA06 2.27–2.28; sympy-checked: η(p, π_i*) = η* exactly and D(η*) = χψ_i):

      η* = M + √(2/3) χ (1−N̄) ψ_i Ω,
      π_i*/p = exp( √(2/3) χ̄ ψ_i Ω / M )                                  (N = 0),
      π_i*/p = [ 1 − √(2/3) χ̄ ψ_i Ω N/M ]^{(N−1)/N} =: B^{(N−1)/N}           (N > 0),        χ̄ = χ/β = χ(1−N̄)/(1−N).   (S.23)

The branch is selected by N (AB06 label the cases "N̄ = N = 0" and "0 ≤ N̄ ≤ N ≠ 0"; N̄ enters only through χ̄). The
N → 0 limit of the power branch is the exp branch (checked numerically). **Guard:** the base B must stay positive:
B > 0 ⇔ ψ_i > −M(1−N)/(√(2/3) |χ| (1−N̄) Ω N) for χ < 0 (K2 parameters: ψ_i > −0.643). Refuse/flag otherwise.

Derivatives of π_i* = Π(p, Ω, ψ_i) (AB06 54, 56 re-derived; sympy-checked, both branches):

      Π_ψ := ∂π_i*/∂ψ_i = √(2/3) χ̄ Ω (1−N) π_i* / (M − √(2/3) χ̄ ψ_i Ω N)        [N = 0: π_i* √(2/3) χ̄ Ω/M]
      Π_Ω := ∂π_i*/∂Ω   = √(2/3) χ̄ ψ_i (1−N) π_i* / (M − √(2/3) χ̄ ψ_i Ω N)      [N = 0: π_i* √(2/3) χ̄ ψ_i/M]
      Π_p := ∂π_i*/∂p   = π_i*/p                                                                            (S.24)
      π*_a := ∂π_i*/∂σ_a |_{ψ_i} = (π_i*/(3p)) δ_a + Π_Ω Ω_a                       (= AB06 54)
      ∂π_i*/∂π_i |_σ = Π_ψ Λ(π_i)/π_i      (paper: Π_ψ λ̃/π_i; fork: Π_ψ λ_c ξ (−π_i/p_a)^ξ/π_i).

---

## 8. Hardening

Rate law (AB06 43; BA06 2.35–2.37 with f(x) = hx) [E]:

      π̇_i = h (π_i* − π_i) ε̇^p_s = √(2/3) h λ̇ (π_i* − π_i) Ω.                                                (S.25)

Hardening modulus (AB06 44–45; BA06 2.31) [E]: H := −(1/λ̇) F_π π̇_i = −M (p/π_i)^{1/(1−N)} √(2/3) h (π_i* − π_i) Ω.
Since p/π_i > 0: H > 0 (hardening) iff π̇_i < 0 iff π_i* < π_i (|π_i*| > |π_i|), H = 0 iff π_i = π_i*, H < 0 (softening)
iff |π_i| > |π_i*|. sgn H = sgn(D* − D) = sgn(η* − η) (BA06 2.32).

Backward-Euler form (AB06 52, Box 2 step 7f) [E]:

      π_i = π_{i,n} + √(2/3) h Δλ (π_i* − π_i) Ω,     π_i* = Π(p, Ω, ψ_i(v, π_i)),   v = v_n exp(tr Δε) = v₀ exp(tr ε_{n+1}).   (S.26)

(v-update: §1.2, G2 owner decision; BA06 Box 2 step 6b has v = v₀(1 + tr ε_{n+1}), superseded.)

Nested scalar residual and its derivative (AB06 61–62) [E, Λ generalised; FD-checked]:

      r(π_i) = π_i − π_{i,n} − √(2/3) h Δλ (π_i* − π_i) Ω,
      r'(π_i) = 1 + √(2/3) h Δλ Ω [ 1 − (Λ(π_i)/π_i) Π_ψ ],     c := r' at the converged π_i.                (S.27)

Solve by Newton from π_{i,n} at every local iterate (p, Ω, v fixed). **Root-selection contract [I; G1]:** without a
cap r(π_i) is monotone: r' = 1 + √(2/3)hΔλΩ[1 − (Λ/π_i)Π_ψ] > 1 because Λ/π_i < 0 and Π_ψ > 0 under the B > 0 guard
of §7, so the Newton root is unique; with the smooth cap Ω depends on π_i and r(π_i) can **fold** (r' changes sign, up
to three roots, §10.2). The nested solve must return the root **continuous with π_{i,n}**: the first sign change of r
found by a scan from π_{i,n} in the direction of −r(π_{i,n}) with a step much smaller than the dip width of §10.2
(≤ 10⁻⁴|π_{i,n}|), refined by safeguarded Newton/bisection inside that bracket; never a factor-2 geometric bracket,
which can enclose the far w = 1 root. A nested solve that fails, or whose root jumps by more than the scan width
between local iterates, rejects the local step (Δλ backtracked, then the increment substepped, §9.1). At the selected
root c = r'(π_i) > 0, so the closed-form sensitivities below (and the CTO of §9.3 built from them) remain valid.
Implicit derivatives of the converged π_i (AB06 57–58, 69; FD-checked to 10⁻⁹ in both CSL modes):

      P_a := ∂π_i/∂σ_a |_{Δλ, v} = (√(2/3) h Δλ / c) [ Ω π*_a + (π_i* − π_i) Ω_a ],
      Π_b := ∂π_i/∂ε^e_b = Σ_a P_a a^e_ab,
      Π_λ := ∂π_i/∂Δλ = (√(2/3) h Ω / c) (π_i* − π_i),
      Π_v := ∂π_i/∂v  = (√(2/3) h Δλ Ω / c) Π_ψ.                                                            (S.28)

(AB06 eq 53 has ∂π_i/∂ε^e_b on both sides; (S.28) is its solved form, = AB06 57.)

---

## 9. Return map (small strain, principal space) and consistent tangent

### 9.1 Algorithm (AB06 Box 2, BA06 Box 2) [E]
Given ε^e_n (tensor), π_{i,n}, v_n, total strain ε_{n+1} (so Δε); v₀ is carried but not used by the step (§1.2):
0. Specific volume: v_{n+1} = v_n exp(tr Δε) (§1.2; the step's v in (S.26)–(S.28) and vfac = v_{n+1} in (S.31)).
1. Trial: ε^{e,tr} = ε^e_n + Δε. Spectral: ε^{e,tr} = Σ_a ε̃_a m^a. (The converged ε^e, σ and ε̇^p share these m^a.)
2. σ^tr_a from §2 at ε̃. **Trial contract [I; G1]: the step is plastic iff F(σ^tr, π_{i,n}) > F_tol := 10⁻¹⁰·|p₀|**;
   otherwise elastic: ε^e = ε^{e,tr}, π_i = π_{i,n}, tangent = a^e (§9.4). The threshold is in stress units (F has
   them) and relative to the reference pressure, the same scale as the r₄ normalisation below; both oracles implement
   this value and the kernel must reproduce it (an elastic/plastic decision that differs between kernel and oracle is
   a contract failure, not a tolerance issue). Neutral increments: a pure shear increment of size h on a coaxial
   yielded state gives F^tr − F_n = O(h²) (measured 1.06e7·h² kPa on the K2 drained TXC state at |p| ≈ 158 kPa, r2_neutral_v0.py),
   so the threshold is crossed at h* ≈ 3e−8; above it the increment is a plastic step with Δλ = O(h²) ≥ 0, below it
   an elastic step with a drift bounded by F_tol. Backward Euler has no elastic/plastic chatter (contrast §12). Else:
3. Unknowns x = (ε^e₁, ε^e₂, ε^e₃, Δλ), start x = (ε̃, 0). Residual (AB06 48):

      r_a(x) = ε^e_a − ε̃_a + Δλ q_a(σ(ε^e), π_i),   a = 1,2,3;      r₄(x) = F(σ(ε^e), π_i),                (S.29)

   where π_i = π_i(ε^e, Δλ) is the converged root of (S.27) at the current iterate (nested Newton; Ω, p from σ(ε^e)).
4. Newton: x ← x − J⁻¹ r until ‖r‖ small (AB06 report 4–5 iterations, quadratic). KKT: Δλ ≥ 0. Local work is
   bounded (iteration cap, backtracking line search on the scaled residual, nested-solve failure or root jump = step
   rejected, §8); **on a refusal the kernel substeps the increment** (halving, to a stated depth, e.g. 2⁸) before it
   reports the refusal upward [I; G1]. A substepped increment returns the **chained** consistent tangent of §9.6
   (owner decision 2026-10-01), never the tangent of its last sub-increment.
5. Update: ε^p_{n+1} = ε^p_n + Δλ Σ_a q_a m^a, σ = Σ_a σ_a m^a, state (σ, e or v, π_i).

Scaling note [I]: r₄ is in stress units, r₁₋₃ in strain; normalise (e.g. r₄/|p₀|) for the convergence test.

### 9.2 The 4×4 Jacobian (AB06 49–51, 57–60), every entry expanded

      J_ab = δ_ab + Δλ [ Σ_c q_ac a^e_cb + q_{a,π} Π_b ],                a, b = 1..3
      J_a4 = q_a + Δλ q_{a,π} Π_λ,
      J_4b = Σ_c f_c a^e_cb + F_π Π_b,
      J_44 = F_π Π_λ,                                                                                       (S.30)

with a^e_cb (S.3), q_ac (S.18), q_{a,π} (S.19), f_c (S.14), F_π (S.13), Π_b, Π_λ (S.28). Checked against a central FD
of (S.29) at 5e−9 relative in fork mode and paper mode (returnmap.py). J is not symmetric.

### 9.3 Consistent tangent in principal directions (AB06 65–69) [E; FD-checked]
b := J⁻¹. The explicit dependence of r on the trial strain ε̃_b at fixed x is through −ε̃_a and through v (§1.2,
∂v_{n+1}/∂ε̃_b = v_{n+1}):

      ∂r_k/∂ε̃_b |_x = −δ_kb [k ≤ 3] + s_k δ_b,    s_k := Δλ q_{k,π} Π_v v (k ≤ 3),   s₄ := F_π Π_v v,   v = v_{n+1},   (S.31)
      ∂x_i/∂ε̃_b = −Σ_k b_ik ∂r_k/∂ε̃_b = b_ib − (Σ_k b_ik s_k) δ_b,
      ã^{ep}_ab := ∂σ_a/∂ε̃_b = Σ_{c≤3} a^e_ac ∂x_c/∂ε̃_b.                                                  (S.32)

(AB06 67 + 69 as written, with v the converged specific volume of the step; BA06 2.71 carries v₀ for small strain,
superseded by the G2 decision of §1.2 [I; G2]. G0/G1 had s_k with v₀; FD-checked with the exponential update in
vexp_fd.py (1): 2.4e−8 and 7.0e−9 at h = 10⁻⁶ against 1.8e−5 / 1.1e−4 for the v₀ form and 1.7e−6 / 3.8e−7 for v_n.)
**ã^{ep} is non-symmetric in general, even for associative
flow (N̄ = N, ρ̄ = ρ)**: the π_i-sensitivities Π_b, Π_λ (through π_i*(p, Ω, ψ_i)) and the v-term Π_v v δ_b in (S.31)
enter the rows and columns differently. Non-associativity only adds to this. **The kernel must never select, and the
shell must never advertise, a symmetric solver/tangent for this material** (checked: major asymmetry 0.106 relative
in returnmap.py for the K2 set).

### 9.4 Closed-form spectral tangent, small strain [I; FD-checked with non-coaxial perturbations, 2e−8]
With m^a = n^a⊗n^a, m^{ab} = n^a⊗n^b from the **trial** strain, ε̃_a the trial elastic principal strains, and the
convention (A⊗B)_{ijkl} = A_ij B_kl, (C:dε)_ij = C_ijkl dε_kl:

      C = dσ/dε = Σ_a Σ_b ã^{ep}_ab m^a ⊗ m^b + (1/2) Σ_{a≠b} g_ab ( m^{ab} ⊗ m^{ab} + m^{ab} ⊗ m^{ba} ),
      g_ab = (σ_a − σ_b)/(ε̃_a − ε̃_b),        repeated eigenvalues (|ε̃_a − ε̃_b| < tol):  g_ab = ã^{ep}_aa − ã^{ep}_ab.   (S.33)

Derivation: σ = Σ σ_a(ε̃) m^a(ε̃); dm^a = Σ_{b≠a} (n^b·dε·n^a)/(ε̃_a − ε̃_b)(n^b⊗n^a + n^a⊗n^b); pairing (a,b) with (b,a)
gives g_ab. The **½** is correct for the small-strain derivative of an isotropic tensor function; the finite-strain
Lie-derivative form (AB06 81) has **no ½** and different γ_ab (§9.5). The repeated-eigenvalue limit follows from
σ_a(ε̃_a, ε̃_b, ·) = σ_b(ε̃_b, ε̃_a, ·) (isotropy). Elastic step: ã^{ep} = a^e, same (S.33) (checked, 4e−9). C has minor
symmetries, no major symmetry: use an unsymmetric solver.

**Repeated-eigenvalue row convention [I; P1; documented contract].** The limit g_ab = ã^{ep}_aa − ã^{ep}_ab is not
symmetric in a ↔ b when ã^{ep} is non-symmetric: the (b,a) row gives g_ba = ã^{ep}_bb − ã^{ep}_ba, and the two coincide
exactly only for *exactly* repeated eigenvalues (swap symmetry of the isotropic map σ_b(ε̃) = σ_a(swap_ab ε̃);
sympy-checked, chain_sympy.py (c)). Inside the switch tolerance |ε̃_a − ε̃_b| < 10⁻¹⁰ they differ by O(|ε̃_a − ε̃_b|),
so the spectral sum has C4_abab ≠ C4_baab at that level (the P1 parity drills measured ~1e−8 relative on C;
kernel_parity/scratch_indep/spin_drill.py). **Contract:** O2 and the kernel compress C4 to 6×6 reading only the rows
with i ≤ j ({00,11,22,01,12,02}), i.e. the shear row (0,1) is C4_01kl, built from g_01 = ã_00 − ã_01, never C4_10kl;
both implement exactly this (O2 `_spectral` + the test/parity `c6_of`, kernel `compress_c4`). The chained tangent of
§9.6 is compressed by the same rows after the a^e assembly of (S.47); inside the band its own limit is the Φ rows
averaged by a^e's input symmetrisation (§9.6, "m = 1"), which equals the (S.33) row only under exact swap symmetry.
Symmetrising g_ab or averaging the two rows in (S.33) would move the kernel–oracle parity at the 1e−8 level and is
not adopted. The limit itself stays the right one: for ε̃_a → ε̃_b the
quotient (σ_a − σ_b)/(ε̃_a − ε̃_b) → ã_aa − ã_ab along the row that is kept.

### 9.5 Finite-strain assembly (for the LogStrain wrapper and for K2 in the oracles) (AB06 81–82) [E]

      c̃ = Σ_a Σ_b c̃_ab m^a⊗m^b + Σ_{a≠b} γ̃_ab (m^{ab}⊗m^{ab} + m^{ab}⊗m^{ba}),
      c̃_ab = ã^{ep}_ab − 2 τ_a δ_ab,     γ̃_ab = (τ_b λ̃_a² − τ_a λ̃_b²)/(λ̃_b² − λ̃_a²),   ε̃_a = ln λ̃_a,            (S.34)

with ã^{ep}_ab = ∂τ_a/∂ε̃_b from (S.32) **unchanged** (since G2 the small-strain (S.31) already carries v = v_{n+1};
before G2 this line read "with v₀ → v"; §1.4). Repeated stretches (|λ̃_a − λ̃_b| < tol, e.g. isotropic states,
which the LogStrain wrapper meets at every start): γ̃_ab → (ã^{ep}_bb − ã^{ep}_ba)/2 − τ_a [I; from τ_b λ̃_a² − τ_a λ̃_b² =
τ_a(λ̃_a² − λ̃_b²) + (τ_b − τ_a)λ̃_a² and dε̃/d(λ̃²) = 1/(2λ̃²); checked numerically, O(Δε̃) convergence].
The total spatial tangent is a^{ep} = c̃ + τ⊕1,
(τ⊕1)_ijkl = τ_jl δ_ik (AB06 p.1534 definition of ⊕). BA06 3.48 warns that earlier papers carry a spurious ½ on the
spin sum in this finite-strain form.

**Re-derived and FD-checked at G1 [I; G1; r2_finite_tangent_fd.py].** a^{ep} is the derivative, at the current
configuration, of the nominal stress referred to that configuration: with f → (1 + hE) f (any E, not necessarily
symmetric), P(E) := τ (1 + hE)^{−T} and dP/dh|₀ = a^{ep}:E. Derivation: b^{e,tr} → (1+hE) b^{e,tr} (1+hE)ᵀ gives, in the
trial eigenbasis with μ_a = λ̃_a² and E_ab := n^a·E·n^b, dμ_b = 2hμ_b E_bb (so dε̃_b = hE_bb and dτ_a = hΣ_b ã^{ep}_ab E_bb)
and dm^a = hΣ_{b≠a} (E_ba μ_a + μ_b E_ab)/(μ_a − μ_b) (n^b⊗n^a + n^a⊗n^b); subtracting hτEᵀ from P and collecting the
(a,b) pair gives γ̃_ab(E_ab + E_ba) + τ_b E_ab with exactly the γ̃_ab of (S.34) — **no ½** — and the diagonal
Σ_b ã_ab E_bb − τ_a E_aa = Σ_b c̃_ab E_bb + (Eτ)_aa; the τ⊕1 term is (Eτ)_ij = E_il τ_lj. Numerically (K2 paper set,
ρ = 0.7/ρ̄ = 0.8, a plastic state 5 steps into the (S.43) protocol with distinct stretches, off the corner): central FD of
P over the nine unit E_kl vs (S.34): 6.2e−7, 6.2e−9, 8.7e−10 for h = 10⁻⁵, 10⁻⁶, 10⁻⁷ (O(h²)); the ½-spin variant
misses by 0.137 and dropping τ⊕1 by 8.6e−3, so the check discriminates both. (S.44) against the direct contraction
n_j a_ijkl n_l: 3.9e−16 over 50 random n. a^{ep} has no minor symmetry in (k,l) (the τ⊕1 term): the acoustic tensor must
be built from the full a^{ep}, never from a symmetrised one. **This FD check is a G1 gate test** (§16.3), not only an
author script.

### 9.6 Chained consistent tangent across substeps [I; P1; owner decision 2026-10-01; sympy- and FD-checked]
Not in AB06/BA06 (their Box 2 has no substepping). A substepped increment (§9.1 step 4) returns the exact derivative
of its **final** stress with respect to the **total** strain increment Δε, propagated through every sub-increment.
The tangent of the last sub-increment alone is **not** that derivative: measured 0.51 (m = 2) to 0.81–0.90 (m = 4)
relative error against the FD of the whole increment on the AMP_STOP cap path, 3.6e−2 on a generic m = 8 step
(chain_fd.py (D)); plan §2.8 records why that cost is not acceptable.

**State map of one sub-increment.** Fractions α_k > 0, k = 0..m−1, Σ_k α_k = 1 (uniform α_k = 1/m on the O2/kernel
ladder; any recursive-halving sequence is covered by the same algebra). State z_k = (ε^e_k [full tensor], π_{i,k},
v_k); v₀ does not enter the chain (§1.2); z_0 = the committed state at n. Sub-increment k: trial ε̃_k = ε^e_k + α_k Δε,
with eigen-pairs (ε̃_a, n^a) and m^a = n^a⊗n^a, m^{ab} = n^a⊗n^b; v_{k+1} = v_k exp(α_k tr Δε) (§1.2, G2; was
v_k + v₀ α_k tr Δε); then (S.29) with
(S.27) nested at (ε̃, π_{i,k}, v_{k+1}) gives x = (ε^e_a, Δλ) and π_{i,k+1}, and ε^e_{k+1} = Σ_a ε^e_a m^a (the m^a of
the trial, §9.1). Write z_{k+1} = Φ(z_k, α_k Δε). The ladder level m and every sub-increment's branch (elastic,
plastic, vertex, cap-active) are *decisions*, not differentiated: Φ is differentiated inside its branch, and the
chain is one-sided across a branch switch exactly as (S.33) is across the elastic/plastic switch.

**Partial derivatives of a plastic sub-increment (implicit function theorem on the converged (S.29) + (S.27)).**
Split (S.30) as J = A + t Π_xᵀ with A := ∂r/∂x|_{π_i} (the terms of (S.30) without Π), t := ∂r/∂π_i|_x =
(Δλ q_{1,π}, Δλ q_{2,π}, Δλ q_{3,π}, F_π), Π_x := ∂π_i/∂x = (Π_1, Π_2, Π_3, Π_λ) of (S.28); b := J⁻¹, u := b t,
κ := Π_x·u, c := r'(π_i) of (S.27), Π_v of (S.28). The 4-vector r has no explicit dependence on π_{i,n} or v (both
enter only through the nested root of the scalar residual of (S.27), written ρ(π_i; x, π_{i,n}, v) here to keep it
apart from r: ∂π_i/∂π_{i,n}|_x = −ρ_{π_{i,n}}/c = 1/c, ∂π_i/∂v|_x = Π_v, since ρ_{π_{i,n}} = −1 and ∂ψ_i/∂v = 1), and at
fixed v ∂r_k/∂ε̃_b|_{x,π_i} = −δ_kb (k ≤ 3), ∂π_i/∂ε̃_b|_x = 0. Solving the 5×5 system row by row (chain_sympy.py (a)):

      ∂x/∂ε̃_b |_{π_{i,n}, v} = b_{·b}  (column b of b),      ∂π_{i,n+1}/∂ε̃_b = w_b := Σ_{c≤3} Π_c b_cb + Π_λ b_4b,
      ∂x/∂π_{i,n} = −u/c,                                  ∂π_{i,n+1}/∂π_{i,n} = (1 − κ)/c,
      ∂x/∂v = −u Π_v,                                      ∂π_{i,n+1}/∂v = (1 − κ) Π_v.                     (S.45)

Relation to (S.31)–(S.32): there v is tied to ε̃ by ∂v_{n+1}/∂ε̃_b = v_{n+1} (§1.2), so the (S.32) column is the ε̃
column of (S.45) plus v_{n+1} times its v column, b_{·b} − u Π_v v_{n+1}, because Σ_k b_ik s_k = Σ_k b_ik t_k Π_v v =
u_i Π_v v (s_k = t_k Π_v v by (S.31); checked symbolically, chain_sympy.py (a) with v₀, vexp_sympy.py (f) with a
general v'(ε̃)). The spin terms are those of §9.4 applied to the map ε̃ ↦ ε^e_{k+1} (an isotropic
tensor function of ε̃ at fixed π_{i,n}, v, with eigenvalue map ∂ε^e_a/∂ε̃_b = b_ab): Φ^ε_ε̃ has the (S.33) form with
diagonal block b_ab (a, b ≤ 3) and spin g^Φ_ab = (ε^e_a − ε^e_b)/(ε̃_a − ε̃_b), limit b_aa − b_ab for |ε̃_a − ε̃_b| < tol.
The π_{i,n} and v columns move only eigenvalues (m^a depends on ε̃ alone): ∂ε^e_{k+1}/∂π_{i,n} = −Σ_a (u_a/c) m^a,
∂ε^e_{k+1}/∂v = −Σ_a u_a Π_v m^a.

**Recursion.** Columns are indexed by the kernel's six Δε components J with input tensors E_J (E_J = e_k⊗e_k for a
normal slot, e_k⊗e_l + e_l⊗e_k for a shear slot: the independent tensor shear component, ε_kl and ε_lk moved
together, so tr E_J = 1 or 0 and 𝕀:E_J = E_J). S^ε_k := ∂ε^e_k/∂Δε_J (a 3×3 per column), S^π_k := ∂π_{i,k}/∂Δε_J,
S^v_{k+1} := ∂v_{k+1}/∂Δε_J = v_{k+1} (Σ_{j≤k} α_j) tr E_J (closed form [I; G2]: v_{k+1} = v_n exp((Σ_{j≤k} α_j) tr Δε)
differentiated; sympy-checked for m = 4 general fractions, vexp_sympy.py (c); the cumulative fraction is the only
extra state, v_{k+1} being the sub-increment's own converged v; before G2 the factor was v₀). With
T_k := S^ε_k + α_k E_J = ∂ε̃_k/∂Δε_J, T̂_bb := n^b·T_k·n^b, S^ε_0 = 0, S^π_0 = 0:

      plastic:  S^ε_{k+1} = Φ^ε_ε̃ : T_k − Σ_a m^a [ (u_a/c) S^π_k + u_a Π_v S^v_{k+1} ],
                S^π_{k+1} = Σ_b w_b T̂_bb + ((1 − κ)/c) S^π_k + (1 − κ) Π_v S^v_{k+1};
      elastic:  S^ε_{k+1} = T_k,   S^π_{k+1} = S^π_k   (Φ^ε_ε̃ = 𝕀, u = 0, w = 0; v still advances).            (S.46)

Returned tangent, with a^e(ε^e_m) in the (S.33) form with ã → a^e and the eigen-data of the **final** ε^e_m (spin
(σ_a − σ_b)/(ε^e_a − ε^e_b), limit a^e_aa − a^e_ab):

      C = dσ_{n+1}/dΔε = a^e(ε^e_m) : S^ε_m,    C[I][J] = (a^e : S^ε_m[J])_ij for I = (i,j), i ≤ j.            (S.47)

Minimal carried state per Δε column: S^ε (6 symmetric components; the reference keeps the full 3×3, below), S^π (1),
the cumulative fraction (1, shared): 6 + 1 + 1 rows × 6 columns. Cost per plastic sub-increment, beyond the return
map itself: b = J⁻¹ (already needed for the CTO), the 4-vectors u = b t and w = Π_x b, two scalars, and per column one
rotation of T_k into the trial basis, the nine-entry spectral product and the rotation back (≈ 10³ flops for the six
columns, far below one local Newton iterate with its nested π_i solve); at the end one a^e assembly and six
contractions. Each sub-increment's chain uses its own (J, c, Π, V); a refused level discards its sensitivities and
the finer level restarts from S_0 = 0.

**Branches.** Elastic sub-increment: (S.46) elastic line. Vertex (§3.2, R < R_tol): Ω = 0 ⇒ Π_x = 0, Π_λ = 0, Π_v = 0,
c = 1, κ = 0; (S.45) collapses to ∂π_{i,n+1}/∂π_{i,n} = 1 (π_i frozen), ∂x/∂π_{i,n} = −u with t = (Δλ β F_{ppi}/3 δ_a,
F_π) and the vertex q_a, and no v column; the deviatoric columns of the FD leave the branch (Ω is discontinuous on
the axis), so only the derivative along 1 is two-sided there (chain_fd.py (F)). Smooth cap (§10.2): nothing changes
in form; the cap enters only through q_a, q_ab, q_{a,π}, Ω, Ω_a, Ω_π of (S.36) in A, t, and through (S.37) in c and
Π_x; at the selected root c > 0 (§8), so (S.45) is well posed. Planar cap (§10.1): the corner is a branch, one-sided.

**m = 1 reduces to (S.33).** S^ε_1 = Φ^ε_ε̃ : E_J − Σ_a m^a u_a Π_v v_1 tr E_J (v_1 = v_{n+1}) has diagonal block
b_ab − u_a Π_v v_{n+1} = the (S.32) ∂x/∂ε̃, so a^e : S^ε_1 has diagonal block ã^{ep}_ab; the spin of the composition of two (S.33)-form operators is
g^A_ab (g^B_ab + g^B_ba)/2 (sympy, chain_sympy.py (b)), here (σ_a − σ_b)/(ε^e_a − ε^e_b) · (ε^e_a − ε^e_b)/(ε̃_a − ε̃_b)
= g_ab of (S.33). Exact for distinct trial eigenvalues (measured 1e−15 and 8e−17 relative, chain_fd.py (C)); an
elastic increment gives a^e exactly. Inside the repeated-eigenvalue band |ε̃_a − ε̃_b| < 10⁻¹⁰ the two limits differ:
(S.33) keeps the row ã_aa − ã_ab, the chain gives (a^e_aa − a^e_ab)·½[(b_aa − b_ab) + (b_bb − b_ba)] (the Φ limit rows
are averaged when a^e symmetrises its input); the two coincide under the exact swap symmetry of isotropy (sympy,
chain_sympy.py (c)) and differ by O(|ε̃_a − ε̃_b|) otherwise, the ~1e−8 level of the §9.4 row-convention note. **Contract:**
the O2 chained reference keeps the full 3×3 column tensors (the (S.33) row convention per operator) and lets each
operator symmetrise its input; the kernel does the same, and compresses (S.47) by the i ≤ j rows of `compress_c4`.
The regression check "C_chain(m = 1) = C_(S.33)" is therefore exact (round-off) at distinct trial eigenvalues and
≤ ~1e−8 inside the band; the FD check does not see the band (h ≫ tol).

**Verification (chain_sympy.py; chain_fd.py on Esmeralda with O2's `_step_once`, fractions held fixed across the FD
points, same branch pattern at every point).** The chain_fd.py numbers below were measured at P1 under the linear
v-update (S^v with v₀); the G2 re-check with the exponential update, S^v = v_{k+1} cum tr E_J and O2's finite-mode
`_step_once` (which already implements v_{k+1} = v_k exp(α_k tr Δε), vfac = v_{k+1}) as the integrator, is
vexp_fd.py (2) on the (B) increment: m = 8 (EEEEPPPP) 1.0e−7, 1.1e−9, 1.5e−9; m = 2 (EP) 7.5e−8, 7.4e−10, 1.4e−9;
α = (½, ¼, ⅛, ⅛) 9.8e−8, 1.5e−9, 1.8e−9; m = 1 7.5e−8, 7.8e−10, 9.7e−10 and = (S.33) to 2.8e−15 (h = 10⁻⁶, 10⁻⁷, 10⁻⁸);
the same with the state's v forced to 1.45 against v₀ = 1.70: 2.7e−8 … 4.2e−10 and 1.9e−15, while the v₀ variant of
S^v sits at 2.2e−6 (real v) and 2.2e−5 … 2.6e−5 (v = 1.45) for every h. Symbolic: the 5×5 block identities behind (S.45) (t parametrised as
A u/(1 − κ) so no inverse is formed), the composition rule of two spectral operators, the (S.32) reduction, the
swap-symmetry limits — all pass. Numerical, kernel 6×6 convention, max over columns of ‖C_J − FD_J‖/‖C_J‖:
(A) AMP_STOP smooth cap (§10.2 defaults, K2 set ρ = 0.7/ρ̄ = 0.8), n = 40: O2's ladder substeps all 30 plastic
increments 11–40 (step 11 at m = 2, 12–40 at m = 4, patterns PP/PPPP, w = 0.02–0.064); chained vs FD 1.1e−5 → 8.6e−8
→ 8.7e−9 (step 11) and 3.7e−6…7.9e−6 → 3.9e−8…7.9e−8 → 7.6e−10…1.2e−8 (steps 12–40) for h = 10⁻⁶, 10⁻⁷, 10⁻⁸ (O(h²)
to the round-off floor); last-sub-increment CTO 0.51 (m = 2), 0.81–0.90 (m = 4). (B) generic plastic increment with
all three shears (no cap, θ 0.57 → 0.71): forced m = 8 (pattern EEEEPPPP) 1.0e−7, 1.1e−9, 2.0e−9; m = 2 (EP) 7.5e−8,
6.9e−10, 1.4e−9; last-sub-increment CTO 3.6e−2 at m = 8 (and equal to the chain at EP, as it must be: an elastic first
half makes T_1 = 𝕀). (C) m = 1: 1.1e−15, 7.8e−17 (plastic), 0 (elastic). (E) non-uniform fractions (a recursive-halving
shape, Σα_k = 1): α = (½, ¼, ⅛, ⅛) on the (B) increment (EPPP) 9.8e−8, 1.6e−9, 2.1e−9; α = (¼, ¼, ¼, ⅛, ⅛) on AMP
step 20 (PPPPP) 4.8e−6, 5.2e−8, 2.7e−9; α = (½, ¼, ¼) on AMP step 11 (PPP) 8.8e−6, 6.6e−8, 8.9e−9. (F) vertex branch,
hydrostatic plastic step of tr Δε = −3e−3 from the apex (no cap, π_i = −46.4758 kPa, exactly repeated trial
eigenvalues, p stays at −100 kPa, π_i frozen): chained m = 1 vs (S.33) 8.4e−17; C:1 = 0 to 2.3e−16 of max|C| for
α = (½, ¼, ¼) and m = 1, FD along 1 agrees to 2.3e−16 … 2.3e−11; the six-column FD switches branch on the deviatoric
columns (one-sided, as stated above). The expected agreement is set by the FD truncation (O(h²), ~1e−6·h/10⁻⁶) down
to the round-off floor ~1e−9; every case meets it.

**Re-measured under the exponential v-update (G2 owner decision 2026-10-01; v_{n+1} = v_n exp(tr Δε), vfac = v_{n+1},
S^v = v_{k+1} Σα tr E_J).** The numbers in the paragraph above were taken at P1 under the linear update and are kept
as that record. The O2 post-decision self-check (`o2_algo/README.md`, "Self-check re-run after the exponential
v-update", Esmeralda, `selfcheck chain`; same ‖C_J − FD_J‖/‖C_J‖ measure, h = 10⁻⁶ / 10⁻⁷ / 10⁻⁸, fractions held fixed)
gives, on the README's own increments (not the (B) increment of the paragraph above, so the generic-increment values
differ from it by the increment, not by the update): (A) AMP_STOP smooth cap, n = 40, all 30 substepped increments
3.4e−6…1.1e−5 / 3.5e−8…8.6e−8 / 2.0e−9…1.5e−8 (last-sub-increment CTO 0.51–0.90, unchanged); (B) generic non-coaxial
plastic increment (fork WW, three shears, θ = 0.82) forced m = 8 2.74e−8 / 1.57e−10 / 9.8e−10, m = 2 2.56e−8 / 3.5e−10 /
1.3e−9; (C) m = 1 vs (S.33) 1.0e−15, ladder m = 1 bit-identical to (S.33): True; (E) non-uniform α (½,¼,⅛,⅛) on (B)
2.71e−8 / 3.4e−10 / 1.6e−9, (¼,¼,¼,⅛,⅛) on AMP step 20 4.76e−6 / 5.1e−8 / 3.9e−9, (½,¼,¼) on AMP step 11
8.77e−6 / 6.5e−8 / 9.7e−9; (F) vertex: chain m = 1 vs (S.33) 4.6e−17, max|C:1|/max|C| 3.8e−16 (m = 1) and 3.2e−16
(α = ½,¼,¼), FD along 1 at 0 and 2.0e−11. Every FD agreement stays at its O(h²)/round-off floor, which is what shows the
v_{k+1} factors are the consistent derivatives (the v₀ form would sit at ~1e−5…1e−4). The symbolic and vexp_fd.py
checks are quoted in the paragraph above.

---

## 10. Q-cap

### 10.1 BA06 planar cap (BA06 2.76–2.79 p.5126) [E]
For η(p, π_i) < χ_cap M (χ_cap user parameter, e.g. 0.10): Q := −p, so q_a = −(1/3) δ_a, q_ab = 0, q_{a,π} = 0,
Ω = 0 (no deviatoric plastic flow ⇒ ε̇^p_s = 0 ⇒ **π_i does not evolve**), ε̇^p_v = −λ̇ (compaction). Residual:
ε^e_a − ε̃_a − Δλ/3 = 0, F = 0. Jacobian: J_ab = δ_ab, J_a4 = −1/3, J_4b = Σ_c f_c a^e_cb, J_44 = 0. ∂r/∂ε̃|_x = (−I, 0).
The switch at η = χ_cap M is a corner of Q (discontinuous q_a and tangent). BA06 note a smooth cap is possible.
Consequence, both caps [I]: hydrostatic compression beyond π_c is perfectly plastic (no volumetric hardening in this
model); bounded p on the compression side. Open item §16.

**Corner sliding [I; G1, measured].** The planar cap replaces the vertex of §3.2 by a corner at η = c₁M, and on a
near-isotropic path the state can be *attracted* to it: on the uncapped side (w = 1, η ≥ c₁M) the flow drives η down,
on the capped side (w = 0, purely volumetric flow) η goes up, so the state **slides along the corner** η = c₁M
(Filippov sliding) with no classical rate solution — neither branch is consistent for a finite time. The rate oracle
stops there (status `cap_sliding`; measured at increment 17 of 40 on the AMP_STOP path with c₁ = 0.10, last state
p = −161.9 kPa, q/|p| = 0.094, i.e. η = ζ(θ)q/|p| at the corner η₁ = 0.12, π_i frozen at −80 kPa since the capped
side has Ω = 0). **Convention at η = c₁M exactly: w = 1** (the
uncapped branch), in both oracles and the kernel. Backward Euler returns to one side or the other per step and
chatters across the corner; the O2 oracle refuses the AMP_STOP step 16 (`local_linesearch` after 2⁸ substeps). This is
the reason for the smooth cap of §10.2 (plan §2.7), which has no corner to slide on.

### 10.2 Smooth cap [I] — fork extension, unified with 10.1 and with "no cap"
Blend weight on the stress ratio of (S.12), η = η(p, π_i) (well defined at q = 0, unlike −q/p):

      t = clamp((η − η₁)/(η₂ − η₁), 0, 1),  η₁ = c₁M ≤ η₂ = c₂M,
      w = S(t) = t³(10 − 15t + 6t²)  (C² quintic; C¹ cubic 3t²−2t³ also admissible),  w_η := dw/dη = S'(t)/(η₂ − η₁),
      S'(t) = 30 t²(1−t)².                                                                                  (S.35)
      q_a = −(1/3) δ_a + w ( q^u_a + (1/3) δ_a ),        q^u_a = (S.17)   [so q_a = q^u_a for w = 1, = −δ_a/3 for w = 0].

Planar cap = c₁ = c₂ = χ_cap (w a step, w_η ≡ 0); no cap = w ≡ 1. Derivatives, with g_a := q^u_a + δ_a/3,
η_p := ∂η/∂p = (F_p − η)/p, η_π := ∂η/∂π_i = F_π/p (FD-checked at 1e−8 to 1e−9 for w = 0.19 and 0.51, row-wise, and a
w_η-dropping mutant is caught):

      q_ab = w q^u_ab + w_η η_p (1/3) g_a δ_b,
      q_{a,π} = w q^u_{a,π} + w_η η_π g_a,
      Ω = w Ω^u,   Ω_a = w Ω^u_a + w_η η_p (1/3) δ_a Ω^u,   Ω_π := ∂Ω/∂π_i = w_η η_π Ω^u,                    (S.36)

and, because Ω now depends on π_i, the nested loop and the sensitivities change to

      r'(π_i) = 1 − √(2/3) h Δλ [ (Π_Ω Ω_π + Π_ψ Λ/π_i − 1) Ω + (π_i* − π_i) Ω_π ]  (= (S.27) when Ω_π = 0),
      P_a = (√(2/3) h Δλ / c) [ Ω ( (π_i*/(3p)) δ_a + Π_Ω Ω_a ) + (π_i* − π_i) Ω_a ],   Π_λ, Π_v unchanged in form.   (S.37)

Here π_i* is (S.23) with the capped Ω; the D*-identity (K1.6) holds exactly wherever w = 1 (η ≥ η₂ < M, so at every
drained peak). Dissipation with the cap: D^p = λ̇[−p(1−w) + w D^p_u/λ̇] ≥ 0 whenever the uncapped D^p_u ≥ 0 (§11) [I].
Recommended defaults: quintic S, c₁ = 0.05, c₂ = 0.15 (to be set by the oracle census; BA06's 0.10 lies between).

**The fold of the nested π_i solve [I; G1, measured; r2_fold.py].** With the cap on, (S.37) reads r' = 1 − G̃ with the
Ω_π terms, and the loop gain

      G := √(2/3) h Δλ |π_i* − π_i| Ω^u w_η |η_π|          (the magnitude of the (π_i* − π_i)Ω_π term of (S.37))

reaches **G ≈ 1 where the state enters the ramp** (η ≈ η₁ = c₁M, w small, w_η up to S'(½)/(η₂ − η₁) = 15.6 at the
defaults): r' changes sign and r(π_i) **folds**, with up to three roots. This is not an absurd-overshoot artefact: on
the mild near-isotropic AMP_STOP path (K2 paper set, ρ = 0.7/ρ̄ = 0.8) at a **volumetric strain step of 7.5e−4** (n = 40)
the converged plastic steps reach G = 0.959, and at the step-8 iterate with Δλ = 7.877e−4 the residual has roots at
π_i = −80.017 (η = 0.065, w = 0.001), −80.483 (η = 0.077, w = 0.021) and −104.7 (η = 0.55, w = 1): the near root is
0.02 % of |π_{i,n}| away, the dip between the first two roots is 0.6 % of |π_{i,n}| (0.07 % at Δλ = 8.0e−4), and a
factor-2 bracket from π_{i,n} (first point −160) encloses the far w = 1 root at −104.7. G falls with the step (0.885,
0.701, 0.351 at n = 80, 160, 320). The backward-Euler solution exists at every step: a monolithic 5-unknown Newton on
(ε^e, Δλ, π_i) converges to a scaled residual of 1e−13 (≤ 6 iterations, n = 40) — but it may land on a different root
(its converged steps show G up to 2.8, i.e. r' < 0 there), so **the backward-Euler problem is non-unique at this step
size** and the root-selection contract of §8 is what makes the kernel's answer the one continuous with π_{i,n}: scan
step ≤ 10⁻⁴|π_{i,n}| (≪ the dip width), no factor-2 bracket, Δλ bounded with backtracking when the nested solve fails or
jumps, substep on refusal (§9.1). The closed-form CTO (§9.3) is kept: at the selected root c = r' > 0 and (S.37) is
exact. With the shipped near-root selection the O2 oracle completes the AMP_STOP path at every n from 20 to 320.

---

## 11. Dissipation

### 11.1 At yield, reading A (the implemented flow rule) [I; numerically swept]
D^p := λ̇ Σ_a σ_a q_a = λ̇ [ p β F_p + q ζ̄ ] (using Σ_a σ_a y_a = 0, Σ_a σ_a n̂_a = R). At yield q = −pη/ζ:

      D^p = −λ̇ p [ (M − η)/(1−N̄) + (ζ̄/ζ) η ] ,     −λ̇ p ≥ 0.                                             (S.38)

The bracket is linear in η; on the surface η ∈ [0, M/N]; at η = 0 it is M/(1−N̄) > 0; at η = M/N it is
M − (M/N)[1 − (ζ̄/ζ)(1−N̄)] ≥ 0 ⇔ **ζ̄(θ)/ζ(θ) ≥ β = (1−N)/(1−N̄) for every θ** (sympy: solve gives r ≥ (1−N)/(1−N̄)).
For WW and GA, ζ̄/ζ is monotone in θ with minimum ρ/ρ̄ at θ = 0 when ρ < ρ̄, and ≥ 1 when ρ ≥ ρ̄ (numeric sweep), so

      **condition A:   N̄ ≤ N   and   ρ/ρ̄ ≥ (1−N)/(1−N̄)   (automatically true if ρ ≥ ρ̄).**                   (S.39)

(N̄ ≤ N **is** forced by (S.38): at θ = π/3, ζ̄/ζ = 1 for every ρ, ρ̄, so the condition at the compression corner is
1 ≥ β ⇔ N̄ ≤ N — the 2-invariant condition of BA06 2.20. Condition A is exactly min_θ ζ̄/ζ ≥ β, of which N̄ ≤ N is the
θ = π/3 member and ρ/ρ̄ ≥ β the θ = 0 member.)

### 11.2 At yield, reading B (AB06 29–34 p.1538–1539) [E]
With η̄ = (ζ̄/ζ)η: D^p = −λ̇ p [ (M − η̄)/(1−N̄) + η̄ ] ≥ 0 ⇔ 1 − (N̄ ζ̄/(N ζ))[1 − (1−N)(p/π_i)^{N/(1−N)}] ≥ 0, the bracket
being in (0, 1] for p ∈ [π_c, 0], hence N̄ζ̄/(Nζ) ≤ 1 for all θ, satisfied if **N̄ ≤ N and ζ̄ ≤ ζ ⇔ ρ ≤ ρ̄** (AB06 32).
AB06 33–34: ρ = (3 − sin φ_c)/(3 + sin φ_c), ρ̄ = (3 − sin ψ_c)/(3 + sin ψ_c), so ρ ≤ ρ̄ ⇔ ψ_c ≤ φ_c (dilatancy angle at
critical ≤ friction angle) — direction confirmed from the paper.

Table (dissipation.py, min over θ ∈ [0, π/3] and η ∈ [0, M/N] of the bracket):

| M, N, N̄, ρ, ρ̄ | min A | min B | cond A | cond B (ρ ≤ ρ̄) |
|---|---|---|---|---|
| 1.2, .4, .2, .7, .8 (K2) | +0.375 | +0.750 | yes | yes |
| 1.2, .4, .4, .7, .8 | **−0.375** | 0.000 | no | yes |
| 1.2, .4, .2, .8, .7 | +0.750 | +0.643 | yes | no |
| 1.33, .3, .2, .71, .9 | **−0.382** | +0.554 | no | yes |

So the paper's "ρ ≤ ρ̄" does **not** guarantee D^p ≥ 0 for the flow rule the paper (and we) implement. **Decided
(owner, G0, 2026-09-30; plan §2.2, §5.2 K1.9): the parser hard-refuses unless N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄)
(condition A), and only WARNS on ρ > ρ̄** (which violates the paper's ψ_c ≤ φ_c reading but not dissipation under
reading A). An independent numeric check (Opus, 1116 parameter cases) found condition A sharp.

**Condition A has two independent parts [I; G1].** (i) **N̄ ≤ N** is the θ = π/3 member of min_θ ζ̄/ζ ≥ β (where
ζ̄/ζ = 1 whatever ρ, ρ̄) and is required **on its own, also when ρ ≥ ρ̄**: with ρ = ρ̄ (any value, e.g. the 2-invariant
ρ = ρ̄ = 1) and N̄ > N the η = M/N endpoint of (S.38) is M(1 − N̄/N) < 0 at every θ, so D^p < 0 near the tension apex
although ρ/ρ̄ = 1 ≥ β trivially holds. (ii) **ρ/ρ̄ ≥ (1−N)/(1−N̄)** is the θ = 0 member; it is implied by ρ ≥ ρ̄ (then
ρ/ρ̄ ≥ 1 ≥ β) but adds a genuine restriction when ρ < ρ̄: the second table row (ρ = 0.7 < ρ̄ = 0.8, N̄ = N so β = 1)
satisfies (i) and fails (ii). The parser must test both; passing one never excuses the other.

### 11.3 Hardening-conjugate term and the dropped coupling (BA06 §2.5 p.5121–5122) [E]
With π_i := ∂Ψ/∂ε^p_s the reduced inequality is σ:ε̇^p − π_i ε̇^p_s ≥ 0 (BA06 2.41). Since ε̇^p_s = √(2/3)λ̇Ω ≥ 0 and
π_i < 0, the second term is ≥ 0 by itself; the inequality holds whenever §11.1 does. Ψ^p is built by integrating
−π̇_i = h(π_i − π_i*)ε̇^p_s twice (BA06 2.43–2.45); because π_i* depends on ε^e (through p and Ω) the stress acquires
the term ∫∫ h ∂π_i*/∂ε^e dε^p_s dε^p_s = O(Δε_s^{p2}) (BA06 2.46–2.47), dropped in the model. BA06 justify it as small
before critical state and exactly zero at it (h(π_i − π_i*) → 0).

Sign note [I; paper-level oddity, not a sheet error]: with π_i = ∂Ψ/∂ε^p_s < 0, Ψ^p is **decreasing** in ε^p_s, and
the conjugate term −π_i ε̇^p_s in BA06 2.41 is positive, i.e. hardening **adds** dissipation rather than storing energy.
In the usual hardening thermodynamics the stored part is positive and subtracts. This follows from BA06's choice of
sign convention for the internal variable (a pressure, negative) and does not affect any equation here; recorded so
that nobody "fixes" the sign in the kernel's dissipation output. The kernel reports D^p = σ:ε̇^p (§11.1) as the
dissipation, not σ:ε̇^p − π_i ε̇^p_s.

### 11.4 Does any of this depend on the CSL form? [I]
No. The proof of §11.1/11.2 uses only: F = 0, the flow direction (S.17), p ≤ 0, λ̇ ≥ 0 and the range of η on the
surface. None of ψ_i, π_i*, λ̃, e₀, λ_c, ξ appears. The hardening term needs only π_i < 0, which holds in both CSL
modes provided π_i* < 0, i.e. the base B > 0 guard of §7 (a ψ_i bound whose *value* depends on the CSL but whose
*form* does not). The dropped coupling term keeps its structure; only ∂π_i*/∂ε^e = Π_p ∂p/∂ε^e + Π_Ω ∂Ω/∂ε^e + Π_ψ ∂ψ_i/∂ε^e
changes numerically (∂ψ_i/∂ε^e = 0 at fixed π_i and v in either mode). **The dissipation argument survives the
power-law CSL unchanged.** The plan's expectation (§2.2) is confirmed; the oracle census remains the evidence.

---

## 12. Continuum rate form (for the O1 rate oracle)

State: ε^e (symmetric tensor), π_i, v (algebraic from total strain). Rate equations (AB06 Box 1; BA06 Box 1) [E]:

      σ = σ(ε^e)  (§2, spectral),   ε̇^e = ε̇ − λ̇ q,   q = Σ_a q_a m^a,   f = Σ_a f_a m^a,
      π̇_i = √(2/3) h λ̇ (π_i* − π_i) Ω,   ψ_i = ψ_i(v, π_i) algebraic,   v̇ = v tr ε̇  (⇔ v = v₀ exp(tr ε)).        (S.40)

(v-rate: §1.2, G2 owner decision [I; G2]; AB06/BA06 Box 1 small strain have v̇ = v₀ tr ε̇, superseded. O1 keeps v
algebraic from the total strain, v = v₀ exp(tr ε); if it integrates v̇ instead, the ODE is v̇ = v tr ε̇. The
multiplier (S.41) and the continuum tangent (S.42) do **not** contain dv/dε: π_i* enters the consistency condition
only through its current value in π̇_i, so they are unchanged.)
Consistency Ḟ = f:σ̇ + F_π π̇_i = 0 with σ̇ = a^e:(ε̇ − λ̇q) gives the closed-form multiplier (AB06 44–45) [E]:

      λ̇ = ⟨ f : a^e : ε̇ ⟩ / ( f : a^e : q + H ),    H = −M (p/π_i)^{1/(1−N)} √(2/3) h (π_i* − π_i) Ω,       (S.41)

with λ̇ = 0 when F < 0, and on the surface by the tie-break below. a^e is the 4th-order elastic tangent in spectral
form: (S.33) with ã^{ep} → a^e (S.3) and ε̃ → ε^e (this includes the spin terms; O1 needs them because ε̇ is a general
tensor).

**Neutral loading and the elastic/plastic tie-break [I; G1; r2_neutral_v0.py].** The textbook rule "λ̇ = 0 if the
numerator N := f:a^e:ε̇ ≤ 0" makes the mode decision **chatter** on a coaxial yielded state under a shear increment:
with f and σ diagonal in the same basis and ε̇ = ε̇₁₂(e₁⊗e₂ + e₂⊗e₁), a^e:ε̇ is off-diagonal and N = 0 **exactly**
(measured 0.0 for both shear directions against N/scale = −0.21 for an axial increment). On the elastic branch F then
grows at second order (F = O(ε̇₁₂²t²)), the yield-crossing event fires, the plastic branch is entered with N ≈ 0,
the rule sends it back to elastic, and so on. Tie-break (contract for the rate oracle):

      scale := ‖f‖ ‖a^e‖ ‖ε̇‖  (Frobenius norms; a^e as the 6×6 Mandel matrix),   tol := 10⁻¹²,
      enter the plastic mode  iff  F ≥ −F_tol and N > +tol·scale, or the elastic branch crosses F = 0 (then N ≥ 0 by
                              construction and λ̇ = ⟨N⟩/den ≥ 0 is consistent even for N = 0);
      leave the plastic mode  iff  N < −tol·scale (genuine unloading);
      |N| ≤ tol·scale on the surface with no mode history → elastic.

The asymmetric thresholds are the **hysteresis**: a neutral increment keeps the mode it is in, so the elastic → plastic
switch happens once (at the crossing) and the plastic branch, whose consistency condition Ḟ = 0 is enforced in rate
form, carries the rotating-axes case with λ̇ ≈ 0 and no drift. F_tol is the surface tolerance of the oracle (10⁻⁸ M|p|).
The denominator den = f:a^e:q + H ≤ 0 is never mapped to elastic: it is the loss-of-uniqueness stop below. Backward
Euler needs no tie-break (the trial threshold of §9.1 decides; F^tr = O(h²) for a neutral increment).
Continuum elastoplastic tangent:

      a^{ep} = a^e − (a^e:q) ⊗ (f:a^e) / ( f:a^e:q + H ).                                                   (S.42)

O1 integrates (S.40)–(S.41) as an ODE in (ε^e, π_i) under strain control; at H = 0 the denominator stays positive
(f:a^e:q > 0 for the parameter ranges of interest — the oracle should assert it). The denominator vanishing is the
loss-of-uniqueness limit (snap-back), to be reported not hidden.

---

## 13. K1 closed forms (plan §5.2), expected values

1. **Isotropic hyperelastic compression** (BA06 energy, ε_s = 0): p(ε_v) = p₀ exp(−(ε_v − ε_{v0})/κ̂). K2 parameters:
   ε_v = −0.01 → p = −100 e¹ = −271.828 kPa.
2. **Closed elastic loop**: W = ∮σ:dε = 0 to ≤ 10⁻¹² (relative to ∮|σ:dε|) and state return to round-off, any loop
   (including non-coaxial), because σ = ∂Ψ/∂ε^e.
3. **ζ corners**: ζ(0) = 1/ρ, ζ(π/3) = 1, ζ'(0) = ζ'(π/3) = 0 for WW and GA on the admissible ranges: WW ρ ∈ (½, 1]
   (ρ = ½ exactly is refused: ζ'(π/3) = −√3, a vertex, §4.2), GA ρ ∈ [7/9, 1] (refused for ρ < 7/9). The same ranges
   are enforced on ρ̄ [G1].
4. **Image point**: at p = π_i, F = 0 ⇔ q = M|π_i|/ζ(θ): q = M|π_i| in TXC (θ = π/3), q = ρM|π_i| in TXE (θ = 0).
5. **Flow rule**: at every plastic step with q > 0 (the vertex rule of §3.2 makes both sides zero on the axis),
   Δε^p_v/Δε^p_s = √(3/2) β F_p/Ω evaluated at the converged state (exact for the backward-Euler map since
   Δε^p = Δλ q(σ_{n+1}, π_{i,n+1})); measured in returnmap.py: −1.2274809384 both ways; O2 oracle over a full drained
   TXC path: max deviation 2.5e−15 [G1]. **The identity holds per increment for backward Euler only; for the rate
   form it holds POINTWISE** [I; G1; r2_k1_protocols.py]: over a finite increment the ratio Δε^p_v/Δε^p_s =
   ∫D dε^p_s / ∫dε^p_s is the ε^p_s-weighted mean of D(t), which lies between neither endpoint in general (measured,
   O1 drained TXC, K2 paper set: a 2e−3 increment near the start gives chord −0.4081 against D(start) = −0.4618 and
   D(end) = −0.3583, 5e−2 from each; 1e−2 increments later 4e−4 to 1.4e−3 from each). A rate-oracle check must
   therefore be pointwise: probe from the state with two small increments d, d/2 and Richardson-extrapolate,
   D(start) ≈ 2r(d/2) − r(d); the chord error is O(d) (2.2e−4, 1.1e−4, 5.5e−5, 2.8e−5 at d = 8e−6 … 1e−6) and the
   extrapolation recovers D(start) to 2.7e−9 (d = 4e−6) and 6.8e−10 (d = 2e−6), i.e. the ODE tolerance.
6. **Peak identity** (stated exactly): whenever π_i = π_i*, D = χψ_i and H = 0 (by construction of π_i*; sympy-checked
   D(η*) = χψ_i). H = 0 is the **stationary point of the yield-surface size**. It coincides with the peak of η only on a
   constant-p path: dη/dt = (∂η/∂p)ṗ + (∂η/∂π_i)π̇_i, so on a conventional drained triaxial path (ṗ ≠ 0) the η-peak
   precedes H = 0 by the term (∂η/∂p)ṗ [I]. **Location protocol [I; G1; r2_k1_protocols.py]:** H = 0 is not a step
   end point, so it is located by the sign change of a := π_i* − π_i between consecutive plastic states k−1, k of the
   coarse path (sgn H = −sgn a, every other factor of H being positive), and the zero is then refined by **bisection
   on the sub-step fraction f ∈ (0, 1) with a single integrator sub-step of size f·Δ from state k−1** (one
   backward-Euler step for O2, one Radau increment for O1 — the same constitutive call, so the located point is a
   genuine state of the integrator, not an interpolant) until |a| ≤ 10⁻¹⁰|π_i|; at that state |D − χψ_i| ≤ tol
   (D from (S.20)/(S.23) closed forms at the state's stress, χψ_i from (S.22)). Measured (O2, drained TXC, K2 paper
   set): bracket at states 20/21 (a = −0.129, +0.431 kPa), 25 bisection sub-steps, |a|/|π_i| = 2.8e−11,
   |D − χψ_i| = 4.0e−11 (D = χψ_i = +0.155730). Secondary coarse check: the linear interpolant of b = D − χψ_i to
   a = 0 is 7.9e−6 against a bracket variation of 3.9e−3 (b(a) ∝ a + O(a²), so the interpolant is quadratically
   small). An end-point check "at the step where H changes sign, |D − χψ_i| ≤ tol" is wrong by O(step).
7. **Undrained critical state**: isochoric ⇒ v (hence e) constant (tr ε = 0 ⇒ v = v₀ exp(0) = v₀ exactly, unchanged
   by the G2 update; the 1e−12 constancy check of the G1 suite stands); critical state ⇔ ψ_i = 0, π_i = π_i* = p, D = 0,
   H = 0, η = M (θ = π/3): p_cs = −p_a ((e₀ − e)/λ_c)^{1/ξ} (fork; requires e < e₀), p_cs = −exp((v_{c0} − v)/λ̃)
   (paper); q_cs = M|p_cs|/ζ(θ) = M|p_cs| in TXC. At the CS, ε̇^p_v = 0 ⇒ ε̇^e_v = 0 ⇒ p stationary: a fixed point.
8. **Drained CS asymptote**: at large ε_s, ψ_i → 0, π_i → p, η → M (i.e. −ζ(θ)q/p → M, so q/|p| → M/ζ(θ)), D → 0, H → 0.
9. **Dissipation**: D^p = Δλ Σ_a σ_a q_a ≥ 0 every plastic step (0 on elastic ones) under (S.39); a set violating
   (S.39) or N̄ > N is refused; ρ > ρ̄ warned (§11.2).
10. **Specific volume on a drained path [I; G2]**: at every committed state v = v₀ exp(tr ε) with ε = ε^e + ε^p the
   total strain (round-off, ≤ 1e−12 relative), in both oracles and the kernel, for any increment sequence and any
   substepping (exp of a sum: the composed sub-increment updates equal the single-increment one exactly,
   vexp_sympy.py (b)–(c)); under LogStrain v = v₀J with J = det F to the same tolerance (§1.4). The closed forms of
   items 5–8 use ψ_i(v, π_i) at that v. The superseded linear update differs by v₀(1 + x − eˣ) ≈ −v₀x²/2, x = tr ε
   (−4.8e−4 at the end of the TXC_paper path, x = 0.0236; vexp_fd.py (3)): a check at 1e−12 discriminates the two
   updates on any path with |tr ε| ≳ 1e−6.

---

## 14. K2 benchmark: AB06 §6.1 stress-point localization (AB06 p.1551–1553)

**Parameters** [E, p.1551]: κ̂ = 0.01, ε^e_{v0} = 0 at p₀ = −100 kPa, μ₀ = 5400 kPa, α₀ = 0; λ̃ = 0.0135, M = 1.2,
N = 0.4, N̄ = 0.2, h = 280; v₀ = 1.59 (initial specific volume), v_{c0} ≈ 1.81 ("≈" in the paper); Willam–Warnke;
case 1: ρ = ρ̄ = 1; case 2: ρ = 0.7, ρ̄ = 0.8. χ is **not stated** in §6.1; α ≈ −3.5 is the paper's only value
(AB06 p.1540, BA06 2.26) — assume χ = −3.5 [I]. Paper CSL, BA06 energy, no cap mentioned (assume none) [I].

**Loading** [E, eq 98–99, read from the page]: relative deformation gradients

      f₁ = diag(1 + λ₂, 1 − λ₁, 1),   f₂ = diag(1, 1 − λ₂, 1 + λ₁),   λ₁ = 1×10⁻³, λ₂ = 4×10⁻⁴,
      F = f₂^{n₂} · f₁^{n₁},  n₁ = 10 (f₁, "mostly compressed"), then f₂ until localization; n = n₁ + n₂.        (S.43)

Diagonal ⇒ principal directions fixed; log strains ε = diag(n₁ ln(1+λ₂), n₁ ln(1−λ₁) + n₂ ln(1−λ₂), n₂ ln(1+λ₁)).

**Known result** [E, p.1551]: ρ = 0.7/ρ̄ = 0.8 localizes at **n = 22**; ρ = ρ̄ = 1 at **n = 26**. Fig 5 plots the
normalised minimum determinant from step 10, crossing zero near those steps (normalisation not stated; presumably by
the step-10 value [I]). Fig 6: the minimum for ρ = 0.7 is at φ = π/2 [E]. **Reading of φ = π/2 [I; G1]:** the
localization normal is perpendicular to the **intermediate** principal direction; under (S.43) that direction is e₁
at the crossing (stretched by f₁, untouched by f₂; by then e₃, stretched by n₂ f₂ steps, is the least compressive axis
and e₂ the most), so n lies in the e₂–e₃ plane. In AB06's
eigenvalue-ordered basis (a = 1 most compressive, 2 intermediate, 3 least) with α = (sinθ sinφ, cosφ, cosθ sinφ) this is
φ = π/2 exactly; both oracles give φ = 90.0°, θ ≈ 35° (two mirror wells ±θ, one normal reported). The earlier wording
"(n in the n¹–n³ plane)" referred to the paper's axis labels and is withdrawn.

**Localization criterion** [E, AB06 83–86, Remark 4]: F(A) = inf_n det A(n) = 0, A_ik = n_j a^{ep}_ijkl n_l with
a^{ep} = c̃ + τ⊕1 (finite strain, §9.5) built from the **consistent** tangent (Remark 4). In the principal basis, with
n = Σ_a α_a n^a, α = (sinθ sinφ, cosφ, cosθ sinφ) (AB06 88, Fig 2):

      Â_aa = α_a² c̃_aa + σ_n + Σ_{c≠a} α_c² γ̃_ca,    Â_ab = α_a (c̃_ab + γ̃_ab) α_b  (a ≠ b),   σ_n := Σ_c τ_c α_c²,   (S.44)

(AB06 85; transcribed, consistent with (τ⊕1)_ijkl = τ_jl δ_ik [I]). Search: coarse sweep of (θ, φ) ∈ [0,π]², then
Newton (AB06 90–97). For a small-strain acoustic tensor drop σ_n and use (S.33) (AB06 Remark 3).

**Reproducibility** [I]: specified well enough — energy, plasticity, CSL parameters, loading, ζ. **Not specified
(the swept unknowns):** (a) the initial π_i (or π_c): the surface "A" in Fig 3 crosses the axis near τ₃'' ≈ −225 kPa
⇒ p_c ≈ −130 kPa (= BA06's plane-strain p_c), giving π_{i,0} = p_c (0.6)^{1.5} ≈ −60.4 kPa; (b) χ (assume −3.5);
(c) v_{c0} ("≈ 1.81"); (d) the crossing criterion: first step with min det ≤ 0, or the linear-interpolated zero
crossing of Fig 5's normalised curve; (e) the (θ, φ) grid. Inferable: the initial stress (point O in Figs 3–4 sits on
the hydrostatic axis at τ'' ≈ −173 = √3·(−100) ⇒ σ₀ = −100 kPa isotropic, ε^e = 0, consistent with ε_{v0} = 0 at p₀),
and v₀ = 1.59 with v = v₀J (which the small-strain kernel now reproduces exactly under LogStrain, §1.2/§1.4).

**K2 gate (plan §5.2, decided at G0)**: the gate is the **ordering** (ρ = 0.7 localizes before ρ = 1) and the **gap**
(n₁.₀ − n₀.₇ ≈ 4 steps) inside a stated band, with a reported sensitivity table over the swept unknowns; exact n is a
sanity check, not a gate. Band [I]: with π_{i,0} ∈ {−60, −80, −100} kPa, χ ∈ {−3.0, −3.5, −4.0}, v_{c0} ∈ {1.80, 1.81,
1.82} and both crossing criteria, require for every combination: n₀.₇ < n₁.₀, gap ∈ [2, 6]; and for the nominal
combination (π_{i,0} = −60.4, χ = −3.5, v_{c0} = 1.81, first-step criterion) n₀.₇ ∈ [19, 25], n₁.₀ ∈ [23, 29]. Report
the full table; a miss outside the band is a finding against the sheet or the oracle, not something to tune away.
The small-strain kernel reproduces this only through the LogStrain wrapper with (S.34); the oracles can run it directly
in finite strain (ε^e_a = ln λ^e_a, v = v₀J, τ). Since G2 the small-strain kernel's v under LogStrain is v₀ exp(ln J) =
v₀J exactly (§1.2); the superseded linear update gave v₀(1 + ln J), off by O(10⁻⁴) here (and the ~25–30× amplified
state error that motivated the decision, §16.6 item 14).

**Why the two oracles do not report the same n, and why that is not a disagreement [I; G1, measured by the K2 lag
diagnosis].** AB06 Remark 4 allows the consistent tangent for the localization check only "for a small enough load
step". At the (S.43) step size backward Euler **lags the continuum by about one step**: a first-order time error of
the state (and of its CTO, which is the derivative of the discrete map). Richardson extrapolation of O2 at 4–32
substeps per nominal step reproduces the O1 crossing to 0.003 step. The reported n therefore depends on the integrator
*and* on the crossing criterion: at the nominal set, O1 gives **22/26** by the first step with min det ≤ 0
(interpolated crossings 21.55/25.11); O2 at one step per nominal increment gives **23/27** first-step, or 22/26 by
rounding its interpolated crossings (22.40/26.47). **Cross-oracle agreement is tested by convergence under
substepping, not by equal n**; the gate stays ordering + gap (which both oracles pass over the whole sensitivity
table), and the paper's 22/26 is the sanity check it was declared to be.

**K2b (paper-mode regression, graphical)** [E, BA06 §4.2, Figs 11–13]: BA06 have no single-point curves; the closest
is the *homogeneous* 3D cube (v = 1.63, 100 kPa lateral pressure, vertical compression) which deforms homogeneously
until localization: peak nominal axial stress ≈ 375 kPa at ≈ 8 % nominal axial strain, volume change ≈ −0.040 m³
(2 m³ specimen, dilative), localization at ≈ 8.7 %. Parameters: Table 1 (κ̂ 0.03, α₀ 0, μ₀ 2000 kPa, p₀ −100 kPa,
ε_{v0} 0), Table 2 (λ̃ 0.04, N 0.4, N̄ 0.2, h 280), v_{c0} = 1.915, p_c = −130 kPa (plane-strain case; assumed the same
for the cube [I]); **M and α are not tabulated**: M = 1.2 from Fig 2 (CSL through q ≈ 120 at p = −100) and AB06 [I],
α = −3.5 [I]. A soft, graphical gate only.

---

## 15. Fork-mode parameter mapping from TIMs' DM04 calibration

| DM04 (TIMs) | LadrunoNORSAND | transfer |
|---|---|---|
| e₀ = 0.83, λ_c = 0.027, ξ = 0.45, p_a | e₀, λ_c, ξ, p_a of (S.22) | **direct** (same e_c(p) form; p_a must be the same numeric value TIMs used, e.g. 101.325 kPa) |
| M_c = 1.3309 | M | **direct**: η = M at the image point in compression (ζ(π/3) = 1) |
| c = M_e/M_c = 0.71 | ρ | **direct**: ζ(0) = 1/ρ ⇒ M_e = ρ M_c ⇒ ρ = c = 0.71. WW admissible (> ½, §4.2); GA refused (< 7/9) |
| G₀ (G = G₀ p_a (2.97−e)²/(1+e) √(p/p_a)), ν | μ₀, κ̂ (BA06) or the HAR slot | **refit**: BA06 gives constant μ₀ (α₀ = 0) and K = −p/κ̂; match at a representative p: μ₀ = G(p_rep, e), κ̂ = −p_rep/K(p_rep) with K = 2(1+ν)G/(3(1−2ν)). √p shear stiffness needs the HAR energy (§2.3) |
| ψ = e − e_c(p) | ψ_i = e − e_c(π_i) | different argument (image pressure): DM04's dilatancy A_d and ψ enter differently; **refit χ** from peak dilatancy vs ψ_i |
| h₀, c_h, n_b, n_d, A_d | h, N, N̄, ρ̄, χ, c₁, c₂ | **refit** (P3): h from pre-peak stiffness/strain to peak; N, N̄ from volumetric curves; ρ̄ from (S.39) and the dilatancy angle (§11.2: ρ̄ = (3 − sin ψ_c)/(3 + sin ψ_c)); default ρ̄ = ρ (deviatoric associativity, Lade & Duncan per AB06 Remark 1) |
| fabric, α_in, z | — | none (monotonic model) |
| initial e | v₀ = 1 + e | direct |

Constraint check for TIMs: with ρ = ρ̄ = 0.71, (S.39) holds for any N̄ ≤ N. If ρ̄ > ρ is chosen, need
0.71/ρ̄ ≥ (1−N)/(1−N̄), e.g. N = 0.3, N̄ = 0.2 ⇒ ρ̄ ≤ 0.811.

---

## 16. Open items, could-not-read, suspected paper typos

### 16.1 G0 decisions (all decided 2026-09-30)
1. **Reading A of Q** (§5.2, AB06 23) is the model. Confirmed by the Adversary's independent re-derivation.
2. **Dissipation bound** (§11): hard refuse unless N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄); WARN on ρ > ρ̄. Owner-approved;
   plan §2.2, §5.2 K1.9 and K2 already say so.
3. **Hydrostatic vertex rule** (§3.2): Ω := 0 at q = 0 in every mode; π_i frozen there.
4. **K2 gate** is ordering + gap in the band of §14, with the sensitivity table; exact n is a sanity check.
5. Still open (not blocking): cap defaults (c₁, c₂, quintic) — set by the oracle census; HAR energy is a slot (§2.3),
   BA06's energy with α₀ = 0 is P0's energy.
6. **G1 (owner, 2026-10-01): Willam–Warnke ρ = ½ exactly is refused**, for ρ and for ρ̄; the admissible WW range is
   (½, 1] (§4.2: at ½ the compression corner is a vertex, ζ'(π/3) = −√3). GA stays [7/9, 1].

### 16.2 Suspected paper typos / inconsistencies
- **AB06 eq 22₁** ∂Q/∂p = β ∂F/∂p contradicts eqs 13–14 + "Q = 0 on the surface" (BA06 2.11) unless ζ̄ = ζ; eqs 30–31
  use the other reading (§5.2). Not a typo in one symbol; an inconsistency between the flow rule and the dissipation
  proof for ρ̄ ≠ ρ. Consequence: AB06 32 "ρ ≤ ρ̄" is not the dissipation condition of the implemented model (§11.1).
- **AB06 eq 64** omits the (2/3)(D₂₂ − q/ε_s) n̂n̂ term; exact for the paper's energy only (§2.1). BA06 3.42 is general.
- **AB06 eq 53** has ∂π_i/∂ε^e_b on both sides (implicit); eq 57 is the solved form. Not a typo, a presentation trap.
- **AB06 eq 8** y = cos3θ/√6, not cos3θ (§1.1). Eqs 19–20 are consistent with it.
- **AB06 eq 12 "smooth"**: WW is C¹ in θ at both meridians but only C² in stress at the compression meridian
  (ζ'''(π/3) ≠ 0, §4.3). No formula error; affects FD-tangent tests only.
- **AB06 eq 55**: the sums over c in θ_c² and θ_c θ_ca are implicit. Written out in (S.21).
- **AB06 eq 42 / Box 2 7(e)** case labels "N̄ = N = 0" / "0 ≤ N̄ ≤ N ≠ 0": the branch is decided by N alone; N̄ enters
  via ᾱ = α/β. Not an error.
- **BA06 Table 2** lacks M and α (§14 K2b). BA06 Table 1 (κ̂ 0.03, μ₀ 2000) differs from AB06 §6.1 (0.01, 5400): two
  parameter sets, not a typo.
- **AB06 §6.1** v_{c0} "≈ 1.81"; χ, π_{i,0}, initial stress and the Fig 5 normalisation not stated (§14).
- **BA06 3.48** flags a spurious ½ in earlier papers' finite-strain spin sum; the small-strain form (S.33) legitimately
  has the ½ (§9.4). Do not "fix" one with the other.

### 16.3 Open items for O1/O2/kernel
- The tension apex p → 0 (q → 0 with η → M/N): the vertex rule of §3.2 applies to the direction; the p-floor rule of
  the plan (§2.7) handles the magnitude and is outside this sheet. R_tol value to be confirmed by the census.
- No volumetric hardening in the cap region (§10.1): isotropic compression beyond π_c is perfectly plastic. Report, do
  not fix in P1.
- The B > 0 guard of §7 (very dense states): decide refuse vs clamp.
- **[G1, superseded]** The nested π_i Newton with the cap on does not merely "have multiple roots at absurd
  overshoots": it **folds** (G ≈ 1) where the state enters the smooth-cap ramp, on a mild near-isotropic path at a
  volumetric strain step of 7.5e−4, and the backward-Euler problem is non-unique there (§10.2). The contract is the
  root-selection rule of §8 (root continuous with π_{i,n}, scan step ≤ 10⁻⁴|π_{i,n}|, no factor-2 bracket), bounded
  Δλ with backtracking on a nested failure or jump, and substepping on refusal (§9.1). The closed-form CTO is kept.
- Undamped local Newton diverged once in the cap checks at an extreme trial state; the kernel needs the usual step
  control (the checklist's bounded local work). Not a formula issue.
- **[G1, closed]** (S.34) is re-derived (§9.5) and (S.34)/(S.44) are FD-checked (r2_finite_tangent_fd.py, 8.7e−10 and
  3.9e−16). **The FD check of the finite-strain tangent (S.34) and of the acoustic tensor (S.44) against a central
  difference of the nominal stress P = τ(1+hE)^{−T} over the nine unit E_kl is a G1 gate test** (tests/, not an author
  script): it must discriminate the ½-spin variant (BA06 3.48) and the missing τ⊕1 term, which it does at 0.137 and
  8.6e−3 against ≤ 1e−8 for the correct form.
- **[G1]** The rate oracle stops at `vertex_reached` (no cap) and `cap_sliding` (planar cap) on near-isotropic paths
  (§3.2, §10.1); these are properties of the continuum model, not oracle defects. The kernel does not see them (it
  substeps and, on the axis, returns with π_i frozen); the smooth cap removes the corner. The O2 oracle's refusal on
  those paths after 2⁸ substeps is the expected report.
- **[G1]** Round-off-negative Δλ on a neutral increment (F^tr = O(h²) just above F_tol, §9.1): the KKT check Δλ ≥ 0
  is strict in both oracles; whether the kernel tolerates Δλ ≥ −tol_λ·(scale) is for the kernel census. Not a formula
  issue.

### 16.4 Verification record (scratchpad/p0a, all run 2026-09-30)
- invariants.py: (S.6)–(S.7), (S.13)–(S.14), (S.16)–(S.21) exact (sympy) or ≤ 1e−15 (random points).
- zeta.py: §4 corner values, θ- and y-derivatives, corner limits, convexity sweeps, ζ̄/ζ minima.
- energy.py: (S.3)–(S.5) exact; (S.3) vs autodiff 2e−12.
- hardening.py: (S.23)–(S.24), (S.27)–(S.28) vs FD ≤ 1e−9, both CSL modes; D*-identity 2e−16.
- returnmap.py: (S.30) vs FD 4.6e−9 (fork), 4.7e−9 (paper); (S.33) vs FD 2.3e−8 (plastic, non-coaxial), 3.7e−9
  (elastic); TXC corner branch O(h²) axisymmetric, O(h) symmetry-breaking (§4.3); D^p ≥ 0; flow-rule identity exact;
  2-invariant Ω = √(3/2).
- cap.py: (S.36)–(S.37) vs FD 1e−8/3e−9 at w = 0.19, 0.51 (rows 1–3 and row 4 separately); w_η-mutant caught.
- dissipation.py: §11 table and the symbolic endpoint conditions.
- G0 fix round (2026-09-30): ζ_yy precision loss re-measured (§3.1: −0.19392 / −0.1974 / −0.49 at θ = 10⁻⁶/10⁻⁷/10⁻⁸,
  ρ = 0.7); γ̃ repeated-stretch limit checked numerically (5823.6 vs 5825.2 at Δε̃ = 3e−4, O(Δε̃)); returnmap.py re-run
  unchanged (Jacobian 4.6e−9, CTO 2.3e−8, major asymmetry 0.106).

G1 revision checks (2026-10-01, run on Esmeralda with the WP-144 venv against the shipped oracles; scripts in
`Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_r2/`, logs in the session scratchpad `r2/`):
- r2_zeta_half.py (sympy): ζ(θ; ½) = 2cosθ, ζ'(π/3; ½) = −√3; ζ'(π/3; ½+s) = 0 for s > 0; ζ''(π/3) = 3(1−ρ²)/(2ρ−1)²
  (→ ∞ at ½⁺; 9.5625 / 3 at 0.7 / 0.8 match the §4.2 table); convexity measure ≡ 0 at ½, +1.7e−4 at 0.501, −1.0e3 at 0.499.
- r2_finite_tangent_fd.py (O2 kernel): (S.34) vs central FD of P = τ(1+hE)^{−T}: 6.2e−7 / 6.2e−9 / 8.7e−10 at
  h = 10⁻⁵/10⁻⁶/10⁻⁷; ½-spin variant 0.137, no-τ⊕1 variant 8.6e−3; (S.44) vs direct contraction 3.9e−16.
- r2_fold.py (O2 kernel): loop gain G at converged plastic steps 0.959 / 0.959 / 0.885 / 0.701 / 0.351 at n = 20 … 320
  on the AMP_STOP smooth-cap path; three roots of r(π_i) at the n = 40 step-8 iterate (−80.017, −80.483, −104.7 kPa),
  dip width 0.47 kPa (Δλ = 7.877e−4) and 0.052 kPa (8.0e−4); monolithic 5-unknown Newton completes n = 40/80/160
  (scaled residual ≤ 1e−12, ≤ 6 iterations) with G up to 2.8 on its own branch.
- r2_o1_stops_vertex.py (O1 + O2): O1 `vertex_reached` at 13/40 (no cap), `cap_sliding` at 17/40 (planar, c₁ = 0.10);
  O2 axial return from the apex: p = −100.000000 held, π_i frozen, Δε^p_v = tr Δε, ε^p_s = 0, D = 0.300/step; O2 refuses
  the same paths at steps 13 (none) / 16 (planar) with `local_linesearch` after 2⁸ substeps (2026-09-30, linear v-update). Re-measured at G2 (2026-10-01) under the exponential v-update v = v₀ exp(tr ε) (§1.2, decision 1 = option b): planar unchanged (step 16, `local_linesearch`); no cap now refuses at step 12 with `local_noconv` after 2⁸ substeps, O2 and the kernel alike (kernel: `SUBSTEPS_EXHAUSTED`, finest `LOCAL_NOCONV`, 256 substeps, same step); the finest reason (and the step, 13 -> 12) of the no-cap stop changed with the v-update; the cause inside the near-vertex local solve was not investigated, both oracles are a refusal at the same corner either way.
- r2_neutral_v0.py (O2 kernel): N = f:a^e:ε̇ = 0.0 exactly for both shear directions on a coaxial TXC yielded state
  (−0.21·scale for axial); F^tr − F_n = 1.06e7·h² kPa (h = 10⁻⁴…10⁻⁶), threshold crossed at h* = 3.2e−8; v₀-vs-v tangent
  FD: 7.0e−7 (v₀) vs 1.1e−4 (v) at v = 1.45, 2.4e−6 vs 2.0e−5 at v = 1.691 (state 0.01 rad off the WW corner).
- r2_k1_protocols.py (O1 + O2): flow-rule chord vs point value (O1): 5e−2 for a 2e−3 increment, O(d) probes,
  Richardson 2.7e−9 / 6.8e−10; O2 increment ratio = point value to 2.5e−15; H = 0 bisection (O2): 25 sub-steps,
  |a|/|π_i| = 2.8e−11, |D − χψ_i| = 4.0e−11; coarse interpolant 7.9e−6 vs bracket variation 3.9e−3.
- No algebra of the G0 sheet failed a re-check; the only formula-level change is the status of (S.34) (transcribed →
  re-derived + FD-checked) and the explicit ζ''(π/3) closed form in §4.2.

G2 revision checks (2026-10-01, run locally with the py3.12 G2 venv against the shipped O2 oracle; scripts in
`Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_vexp/`):
- vexp_sympy.py (sympy, all exact): ∂(v₀e^{tr ε})/∂ε_a = v, no shear dependence, ∂v/∂ε^e_a = 0 at fixed total strain;
  ∂v_{n+1}/∂ε̃_b = v_{n+1}; group property v₀e^{tr ε_n}e^{tr Δε} = v₀e^{tr ε_{n+1}}; S^v_{k+1} = v_{k+1}(Σ_{j≤k}α_j) tr E_J
  for m = 4 general fractions and all six columns; v₀(1 + x − eˣ) = −v₀x²/2 − v₀x³/6 + O(x⁴), zero slope at 0;
  v₀e^{ln J} = v₀J; IFT column b_{·b} − uΠ_v v' for a general v'(ε̃).
- vexp_fd.py (O2, finite-mode `_step_once` as the integrator of the new rule): (S.31) with vfac = v_{n+1} vs central FD of
  the return map 2.4e−6 / 2.4e−8 (h = 10⁻⁵ / 10⁻⁶) on the real state (v_n = 1.691, v₀ = 1.701) and 7.0e−7 / 7.0e−9 on a
  mid-path state (v_n = 1.45); vfac = v₀ 1.8e−5 and 1.1e−4, vfac = v_n 1.7e−6 and 3.8e−7 (h = 10⁻⁶). Chain (S.46) with
  S^v = v_{k+1} cum tr E_J vs FD of the whole increment, (B) increment: 1.0e−7 / 1.1e−9 / 1.5e−9 (m = 8), 7.5e−8 / 7.4e−10 /
  1.4e−9 (m = 2), 9.8e−8 / 1.5e−9 / 1.8e−9 (α = ½,¼,⅛,⅛), 7.5e−8 / 7.8e−10 / 9.7e−10 (m = 1; = (S.33) to 2.8e−15), h = 10⁻⁶
  / 10⁻⁷ / 10⁻⁸; mid-path v = 1.45: 2.7e−8 … 4.2e−10, (S.33) 1.9e−15; the v₀ variant of S^v 2.2e−6 (real v), 2.2e−5 …
  2.6e−5 (v = 1.45). Gap on the TXC_paper path end: tr ε = 0.0236, v₀(1 + x − eˣ) = −4.76e−4 (−v₀x²/2 = −4.73e−4).

### 16.5 G0 review record (2026-09-30)
- **Adversary (Fable, independent re-derivation): PASS-WITH-FIXES**, seven items, all applied in this revision:
  (1) §3.2 hydrostatic guard was self-contradictory → vertex rule Ω := 0; (2) §14 K2 "±1 step" → ordering/gap band with
  sensitivity table; (3) §11.1 parenthetical on N̄ ≤ N corrected (forced by θ = π/3); (4) §3.1 double-precision loss of
  ζ_yy near corners documented, same-θ rule, no isolated unit tests; (5) §9.3 CTO non-symmetric even when associative,
  no symmetric solver; (6) (S.34) repeated-stretch limit of γ̃_ab added; (7) §11.3 BA06 conjugate-term sign oddity noted.
- **Independent numeric check (Opus): PASSED** — dissipation reading A sharp over 1116 cases; 4×4 Jacobian 1.9e−9 vs
  FD; CTO 9.8e−11 (5e−10 on the shear columns).
- **Owner**: approved the refusal rule (hard refuse unless N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄); WARN on ρ > ρ̄).

### 16.6 G1 revision record (2026-10-01)
Twelve items raised by the G1 oracle/test round (O1 rate oracle, O2 return map, the K1/K2 gate tests and the cap and
K2-lag triages), applied in one pass; item 13 added at P1 (owner decision 2026-10-01); item 14 added at G2 (owner
decision 2026-10-01, exponential v-update). Tags: [E] read in the source, [I; G1] / [I; P1] / [I; G2] this revision's
derivation or decision.

| # | section | change | tag | evidence |
|---|---|---|---|---|
| 1 | §1.3, §4, §4.2, §13.3, §15, §16.1 | WW admissible range (½, 1] for ρ and ρ̄; ρ = ½ refused (ζ = 2cosθ, ζ'(π/3) = −√3, vertex); "ζ'(π/3) = 0" restricted to ρ > ½; ζ''(π/3) = 3(1−ρ²)/(2ρ−1)² closed form; all ζ ranges apply to ρ̄ | [I; G1], owner decision | r2_zeta_half.py |
| 2 | §3.2 | vertex rule restated: F := pη at q = 0, f_a = F_p/3; q → 0 in finite time on a near-isotropic path, no consistent rate solution after it (O1 `vertex_reached`); backward Euler returns along the axis to π_c with π_i frozen; hydrostatic compression beyond π_c perfectly plastic | [I; G1] | r2_o1_stops_vertex.py |
| 3 | §10.1 | planar cap: Filippov sliding on the corner η = c₁M (O1 `cap_sliding`), convention w = 1 at η = c₁M exactly; motivates the smooth cap (plan §2.7) | [I; G1] | r2_o1_stops_vertex.py |
| 4 | §8, §10.2, §16.3 | fold of the nested π_i solve (loop gain G ≈ 1 entering the ramp) at a volumetric strain step of 7.5e−4; BE solution exists (monolithic 5-unknown Newton) but is non-unique; root-selection contract (root continuous with π_{i,n}, scan step ≤ 10⁻⁴|π_{i,n}|, no factor-2 bracket), bounded Δλ with backtracking, substep on refusal; closed-form CTO kept | [I; G1] | r2_fold.py |
| 5 | §12 | neutral-loading tie-break for the rate form: scale = ‖f‖‖a^e‖‖ε̇‖, tol = 10⁻¹², asymmetric enter/leave thresholds (hysteresis); N = 0 exactly for a shear increment on a coaxial yielded state | [I; G1] | r2_neutral_v0.py |
| 6 | §9.1 | trial contract F^tr > 10⁻¹⁰|p₀| made explicit; neutral increments F^tr = O(h²), h* ≈ 3e−8; substepping on refusal | [I; G1] | r2_neutral_v0.py |
| 7 | §1.2 | v₀ (initial specific volume) is a separate state datum from v; a mid-path state must carry it; measured tangent error when v is used instead | [I; G1] | r2_neutral_v0.py |
| 8 | §13.5 | flow-rule identity per increment for BE only, pointwise for the rate form (ε^p_s-weighted mean over an increment); rate-oracle check must be Richardson-probed | [I; G1] | r2_k1_protocols.py |
| 9 | §13.6 | H = 0 location protocol: sign change of π_i* − π_i, bisection with single sub-steps to |a| ≤ 10⁻¹⁰|π_i| | [I; G1] | r2_k1_protocols.py |
| 10 | §11.2 | condition A has two independent parts; N̄ ≤ N required on its own, also when ρ ≥ ρ̄ | [I; G1] | §11 table (G0) |
| 11 | §14 | AB06 Remark 4 "small enough load step": BE lags the continuum by ≈ 1 step; O1 22/26, O2 23/27 first-step (22/26 by rounding the interpolated crossing); agreement tested by convergence under substepping; Fig 6 φ = π/2 reread as n ⊥ the intermediate principal direction (e₁ under (S.43), n in the e₂–e₃ plane; both oracles φ = 90.0°, θ ≈ 35°, two mirror wells) | [I; G1] | K2 lag diagnosis (tests/scratch_k2lag), tests/out/k2_sensitivity.md |
| 12 | §9.5, §16.3 | (S.34) re-derived (nominal-stress derivative, no ½) and FD-checked with (S.44); the FD check is now a G1 gate test | [I; G1] | r2_finite_tangent_fd.py |
| 13 | §9.1, §9.4, §9.6 | **P1 owner decision: chained tangent.** A substepped increment returns the exact derivative of its final stress w.r.t. the total Δε, chained through every sub-increment (any fractions α_k, Σα_k = 1): state map z_{k+1} = Φ(z_k, α_kΔε), IFT columns (S.45) (−J⁻¹∂r/∂(π_{i,n}, v) through the nested solve: −u/c, −uΠ_v, (1−κ)/c, (1−κ)Π_v), recursion (S.46), assembly C = a^e(ε^e_m):S^ε_m (S.47); m = 1 reduces to (S.33) exactly at distinct trial eigenvalues; the repeated-eigenvalue limit of (S.33) keeps the i ≤ j row (C4_01kl from g_01 = ã_00 − ã_01; O2 and kernel; ~1e−8 effect inside the 10⁻¹⁰ band) as the documented contract; supersedes "tangent of the last sub-increment" (plan §2.8, O2 `step`, kernel `step_ex`, 0.5–0.9 off the FD) | [I; P1], owner decision 2026-10-01 | chain_sympy.py (all OK), chain_fd.py at P1 (linear v-update): AMP_STOP n = 40 m = 2–4 → 4e−8…9e−8 at h = 10⁻⁷, generic m = 8 → 1.1e−9, non-uniform α → 1.6e−9…6.6e−8, vertex 8e−17; **re-measured under the exponential v-update** (G2, 2026-10-01; o2_algo/README "Self-check re-run after the exponential v-update", at h = 10⁻⁷): AMP_STOP n = 40 3.5e−8…8.6e−8, generic m = 8 1.6e−10, non-uniform α 3.4e−10…6.5e−8, m = 1 vs (S.33) 1.0e−15, vertex 4.6e−17 |
| 14 | §1.2, §1.3, §1.4, §8 (S.26), §9.1, §9.3 (S.31), §9.5, §9.6 (state map, (S.45) relation, S^v, m = 1), §12 (S.40), §13.7, §13.10, §14 | **G2 owner decision (2026-10-01, decision 1 = option b): exponential v-update.** v = v₀ exp(tr ε) ⇔ v_{n+1} = v_n exp(tr Δε), dv/dε = v·1 (not v₀·1) everywhere: s_k = t_k Π_v v_{n+1} in (S.31) (the converged v of the step; FD discriminates v_{n+1} from v_n and v₀), S^v_{k+1} = v_{k+1} cum tr E_J in (S.46), v̇ = v tr ε̇ in (S.40) ((S.41)–(S.42) unchanged: no dv/dε in them), §1.4 v-mapping identical to finite strain (no v₀ → v substitution left; §9.5 "with v₀ → v" withdrawn), §13.7 isochoric endpoint unchanged (v = v₀ exactly), new §13.10 drained-path identity. v₀ stays a committed datum (v = v₀ exp(tr ε), initialState, revertToStart) but enters no derivative; the G1 item 7 warning is inverted. Motivation (G2 measurement, test_g2_logstrain.py): under LogStrain the linear update read v = v₀(1 + x) against the finite oracle's v₀eˣ, x = ln J, a gap v₀(1 + x − eˣ) (closed form to 1e−10) amplified ~25–30× into τ and π_i: 2.2e−3 / 2.3e−3 at 20 % drained TXC, 5.6e−3 / 6.0e−3 in TXE fork mode; the exponential update makes LogStrain(LadrunoNorSand) exact, v = v₀J, and leaves small strain unchanged to first order. Supersedes BA06 Box 2 step 6b / Box 1 v̇ = v₀ tr ε̇ [E] and BA06 2.71's v₀. Implementation DONE (G2, 2026-10-01; kernel `LadrunoNorSandKernel.h`, O2 `_step_once` small-strain branch and `chain_propagate`, O1 v̇ all moved to the exponential update with vfac = v_{n+1}; parity re-run) | [I; G2], owner decision 2026-10-01 | vexp_sympy.py (all exact), vexp_fd.py: (S.31) 2.4e−8 / 7.0e−9 (v₀ form 1.8e−5 / 1.1e−4), chain m = 8/2/non-uniform/1 → 1.1e−9 … 1.8e−9 at h = 10⁻⁷…10⁻⁸ (v₀ variant 2.2e−6 … 2.6e−5), m = 1 = (S.33) to 2.8e−15 |

Items 1–13: no FD-checked algebra of the G0 sheet changed; status of (S.34): transcribed → re-derived and FD-checked
(8.7e−10). Item 14 is the first formula-level change to FD-checked algebra: the v-factor in (S.31) (v₀ → v_{n+1}) and
in the S^v closed form of (S.46) (v₀ → v_{k+1}), both re-verified (§16.4, G2 checks).
