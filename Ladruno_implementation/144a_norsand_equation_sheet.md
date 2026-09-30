---
title: "WP-144a — LadrunoNORSAND equation sheet (AB06/BA06 re-derived, curved CSL, WW ζ, Q-cap)"
project: Ladruno
type: equation sheet
status: "G0 fix round applied 2026-09-30 (Adversary PASS-WITH-FIXES, 7 items; independent Opus numeric check PASSED; owner approved the refusal rule). Every derivative sympy/FD-checked; scripts in the session scratchpad p0a/."
owner: nmora (Deriver: Fable, P0a)
related:
  - "[[144_ladruno_norsand_plan]] (design §2, oracles §5, roster §6)"
  - "[[134_sanisand_reference_integrator]] (the O1 oracle template)"
tags: [equation-sheet, norsand, critical-state, hyperelasticity, sand, wp-144]
updated: 2026-09-30
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
- Specific volume v = 1 + e. Small strain: v = v₀ (1 + tr ε) with ε the **total** strain and v₀ the initial specific
  volume (BA06 Box 2 step 6b) [E]. Hence ∂v/∂ε_a = v₀ for every principal component, and inside the local return
  (total strain fixed) ∂v/∂ε^e_a = 0.

### 1.3 Parameters (one table for the whole sheet)

| symbol | meaning | source |
|---|---|---|
| p₀, κ̂, ε^e_{v0}, μ₀, α₀ | BA06 energy: reference pressure (<0), elastic compressibility, reference volumetric strain, shear modulus, coupling | BA06 2.3 |
| M | critical stress ratio in compression (θ = π/3) | AB06 10 |
| N, N̄ | curvature of F and of Q on the meridian plane (0 ≤ N̄ ≤ N < 1) | AB06 10, 14 |
| ρ, ρ̄ | ellipticity of F and of Q (WW: ½ ≤ ρ ≤ 1; GA: 7/9 ≤ ρ ≤ 1) | AB06 11–12 |
| β := (1−N)/(1−N̄) ≤ 1 | volumetric non-associativity | AB06 p.1537, BA06 2.29 |
| χ (< 0, ≈ −3.5) | maximum-dilatancy coefficient, D* = χ ψ_i. **This is BA06's α (AB06's α); renamed to avoid AB06's χ = ‖ξ‖.** χ̄ := χ/β (AB06's ᾱ). | BA06 2.26, 2.29 |
| h | hardening constant (280 in both papers) | AB06 43 |
| λ̃, v_{c0} | "paper" CSL: v_c = v_{c0} − λ̃ ln(−p) | AB06 41 |
| e₀, λ_c, ξ, p_a | "fork" CSL: e_c = e₀ − λ_c (−p/p_a)^ξ (DM04 form) | plan §2.3 |
| c₁, c₂ | Q-cap blend bounds, η₁ = c₁M ≤ η₂ = c₂M (planar cap: c₁ = c₂; BA06's χ_cap = 0.10 → c₁ = c₂ = 0.10) | BA06 2.76, §10 |
| v₀ | initial specific volume (state input) | BA06 Box 2 |

### 1.4 Finite-strain mapping (stated once; AB06 §3, BA06 §3.3) [E]
Replace ε^e_a by the principal elastic logarithmic stretches ε^e_a = ln λ^e_a, σ_a by the principal Kirchhoff stresses
τ_a, the trial elastic strain ε^{e,tr} by ε̃_a = ln λ̃_a from b^{e,tr} = f_{n+1} b^e_n f^t_{n+1}, and v = v₀ det F (so
∂v/∂ε̃_a = v instead of v₀; BA06 p.5132). The local residual, Jacobian, nested π_i loop and ã^{ep}_ab are then
**identical** (AB06 46–69; BA06 p.5130 "identical"). Only the assembly of the spatial tangent differs (§9.5).

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

### 3.2 q = 0 (hydrostatic): the vertex rule [I; decided at G0]
R = 0 makes n̂, y, n̂_ab, y_ab undefined. F = pη there and ζ does not enter. Q has a vertex on the axis: the
deviatoric flow √(3/2) ζ̄ n̂_a + ζ̄_y q y_a has magnitude √(2/3)Ω → √(2/3)·√(3/2)ζ̄ ≠ 0 along **every** non-axial
approach, but no direction on the axis itself. **Rule (all modes, cap or no cap): if R < R_tol at the trial state or
at a local iterate, set n̂ := 0, y_a := 0, n̂_ab := 0, y_ab := 0, Ω := 0, Ω_a := 0, so q_a = (1/3)βF_p δ_a (purely
volumetric), ε^p_s = 0 and π_i stays at π_{i,n}.** Justification: (i) it is the limit w → 0 of the cap (S.36), so the
capped and uncapped models agree on the axis; (ii) it keeps (S.20) exact (both sides zero) instead of the
contradictory "no deviatoric flow but Ω = √(3/2)ζ̄"; (iii) the return from an axial trial state beyond π_c is then the
1-D axial return to π_c with frozen π_i, which is what BA06's planar cap does. Ω is therefore discontinuous at the
axis (a vertex, not a smoothness bug); Newton iterates that cross R_tol are a corner problem and fall under the bounded
local-work refusal. K1.5 (flow-rule identity) is tested only at plastic steps with q > 0, where (S.20) has no 0/0.
K2 (no cap) starts isotropic at −100 kPa **inside** the surface (elastic), so the rule never fires there.
R_tol: 10⁻⁸·|p| (relative) [I; to be confirmed by the oracle census].

---

## 4. ζ(θ, ρ): shape functions

Both satisfy ζ(0) = 1/ρ (tension corner), ζ(π/3) = 1 (compression corner), ζ'(0) = ζ'(π/3) = 0 (sympy-checked).
The deviatoric section is the polar curve r(θ) = 1/ζ(θ); it is convex iff r² + 2r'² − r r'' ≥ 0 on [0, π/3]
(checked numerically on a 4001-point grid; sign flips confirmed at the ranges quoted by AB06).

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
Convex for ½ ≤ ρ ≤ 1 (numerically: min = −8.9e−15 at ρ = 0.5, +0.011 at 0.55, −214 at 0.45). Derivatives in θ:
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

      π_i = π_{i,n} + √(2/3) h Δλ (π_i* − π_i) Ω,     π_i* = Π(p, Ω, ψ_i(v, π_i)),   v = v₀(1 + tr ε_{n+1}).   (S.26)

Nested scalar residual and its derivative (AB06 61–62) [E, Λ generalised; FD-checked]:

      r(π_i) = π_i − π_{i,n} − √(2/3) h Δλ (π_i* − π_i) Ω,
      r'(π_i) = 1 + √(2/3) h Δλ Ω [ 1 − (Λ(π_i)/π_i) Π_ψ ],     c := r' at the converged π_i.                (S.27)

Solve by Newton from π_{i,n} at every local iterate (p, Ω, v fixed). Implicit derivatives of the converged π_i
(AB06 57–58, 69; FD-checked to 10⁻⁹ in both CSL modes):

      P_a := ∂π_i/∂σ_a |_{Δλ, v} = (√(2/3) h Δλ / c) [ Ω π*_a + (π_i* − π_i) Ω_a ],
      Π_b := ∂π_i/∂ε^e_b = Σ_a P_a a^e_ab,
      Π_λ := ∂π_i/∂Δλ = (√(2/3) h Ω / c) (π_i* − π_i),
      Π_v := ∂π_i/∂v  = (√(2/3) h Δλ Ω / c) Π_ψ.                                                            (S.28)

(AB06 eq 53 has ∂π_i/∂ε^e_b on both sides; (S.28) is its solved form, = AB06 57.)

---

## 9. Return map (small strain, principal space) and consistent tangent

### 9.1 Algorithm (AB06 Box 2, BA06 Box 2) [E]
Given ε^e_n (tensor), π_{i,n}, total strain ε_{n+1} (so Δε), v₀:
1. Trial: ε^{e,tr} = ε^e_n + Δε. Spectral: ε^{e,tr} = Σ_a ε̃_a m^a. (The converged ε^e, σ and ε̇^p share these m^a.)
2. σ^tr_a from §2 at ε̃; F(σ^tr, π_{i,n}) < 0 → elastic: ε^e = ε^{e,tr}, π_i = π_{i,n}, tangent = a^e (§9.4). Else:
3. Unknowns x = (ε^e₁, ε^e₂, ε^e₃, Δλ), start x = (ε̃, 0). Residual (AB06 48):

      r_a(x) = ε^e_a − ε̃_a + Δλ q_a(σ(ε^e), π_i),   a = 1,2,3;      r₄(x) = F(σ(ε^e), π_i),                (S.29)

   where π_i = π_i(ε^e, Δλ) is the converged root of (S.27) at the current iterate (nested Newton; Ω, p from σ(ε^e)).
4. Newton: x ← x − J⁻¹ r until ‖r‖ small (AB06 report 4–5 iterations, quadratic). KKT: Δλ ≥ 0.
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
b := J⁻¹. The explicit dependence of r on the trial strain ε̃_b at fixed x is through −ε̃_a and through v (§1.2):

      ∂r_k/∂ε̃_b |_x = −δ_kb [k ≤ 3] + s_k δ_b,    s_k := Δλ q_{k,π} Π_v v₀ (k ≤ 3),   s₄ := F_π Π_v v₀,        (S.31)
      ∂x_i/∂ε̃_b = −Σ_k b_ik ∂r_k/∂ε̃_b = b_ib − (Σ_k b_ik s_k) δ_b,
      ã^{ep}_ab := ∂σ_a/∂ε̃_b = Σ_{c≤3} a^e_ac ∂x_c/∂ε̃_b.                                                  (S.32)

(AB06 67 + 69 with v → v₀ for small strain, BA06 2.71.) **ã^{ep} is non-symmetric in general, even for associative
flow (N̄ = N, ρ̄ = ρ)**: the π_i-sensitivities Π_b, Π_λ (through π_i*(p, Ω, ψ_i)) and the v-term Π_v v₀ δ_b in (S.31)
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

### 9.5 Finite-strain assembly (for the LogStrain wrapper and for K2 in the oracles) (AB06 81–82) [E]

      c̃ = Σ_a Σ_b c̃_ab m^a⊗m^b + Σ_{a≠b} γ̃_ab (m^{ab}⊗m^{ab} + m^{ab}⊗m^{ba}),
      c̃_ab = ã^{ep}_ab − 2 τ_a δ_ab,     γ̃_ab = (τ_b λ̃_a² − τ_a λ̃_b²)/(λ̃_b² − λ̃_a²),   ε̃_a = ln λ̃_a,            (S.34)

with ã^{ep}_ab = ∂τ_a/∂ε̃_b from (S.32) with v₀ → v. Repeated stretches (|λ̃_a − λ̃_b| < tol, e.g. isotropic states,
which the LogStrain wrapper meets at every start): γ̃_ab → (ã^{ep}_bb − ã^{ep}_ba)/2 − τ_a [I; from τ_b λ̃_a² − τ_a λ̃_b² =
τ_a(λ̃_a² − λ̃_b²) + (τ_b − τ_a)λ̃_a² and dε̃/d(λ̃²) = 1/(2λ̃²); checked numerically, O(Δε̃) convergence].
The total spatial tangent is a^{ep} = c̃ + τ⊕1,
(τ⊕1)_ijkl = τ_jl δ_ik (AB06 p.1534 definition of ⊕). BA06 3.48 warns that earlier papers carry a spurious ½ on the
spin sum in this finite-strain form. I have not re-derived (S.34); it is transcribed.

---

## 10. Q-cap

### 10.1 BA06 planar cap (BA06 2.76–2.79 p.5126) [E]
For η(p, π_i) < χ_cap M (χ_cap user parameter, e.g. 0.10): Q := −p, so q_a = −(1/3) δ_a, q_ab = 0, q_{a,π} = 0,
Ω = 0 (no deviatoric plastic flow ⇒ ε̇^p_s = 0 ⇒ **π_i does not evolve**), ε̇^p_v = −λ̇ (compaction). Residual:
ε^e_a − ε̃_a − Δλ/3 = 0, F = 0. Jacobian: J_ab = δ_ab, J_a4 = −1/3, J_4b = Σ_c f_c a^e_cb, J_44 = 0. ∂r/∂ε̃|_x = (−I, 0).
The switch at η = χ_cap M is a corner of Q (discontinuous q_a and tangent). BA06 note a smooth cap is possible.
Consequence, both caps [I]: hydrostatic compression beyond π_c is perfectly plastic (no volumetric hardening in this
model); bounded p on the compression side. Open item §16.

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
      π̇_i = √(2/3) h λ̇ (π_i* − π_i) Ω,   ψ_i = ψ_i(v, π_i) algebraic,   v̇ = v₀ tr ε̇.                        (S.40)

Consistency Ḟ = f:σ̇ + F_π π̇_i = 0 with σ̇ = a^e:(ε̇ − λ̇q) gives the closed-form multiplier (AB06 44–45) [E]:

      λ̇ = ⟨ f : a^e : ε̇ ⟩ / ( f : a^e : q + H ),    H = −M (p/π_i)^{1/(1−N)} √(2/3) h (π_i* − π_i) Ω,       (S.41)

with λ̇ = 0 when F < 0 or the numerator is ≤ 0. a^e is the 4th-order elastic tangent in spectral form: (S.33) with
ã^{ep} → a^e (S.3) and ε̃ → ε^e (this includes the spin terms; O1 needs them because ε̇ is a general tensor).
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
3. **ζ corners**: ζ(0) = 1/ρ, ζ(π/3) = 1, ζ'(0) = ζ'(π/3) = 0 for WW and GA; WW convex for ρ ∈ [½, 1], GA for [7/9, 1]
   (parser refuses outside; GA refused for ρ < 7/9).
4. **Image point**: at p = π_i, F = 0 ⇔ q = M|π_i|/ζ(θ): q = M|π_i| in TXC (θ = π/3), q = ρM|π_i| in TXE (θ = 0).
5. **Flow rule**: at every plastic step with q > 0 (the vertex rule of §3.2 makes both sides zero on the axis),
   Δε^p_v/Δε^p_s = √(3/2) β F_p/Ω evaluated at the converged state (exact for the backward-Euler map since
   Δε^p = Δλ q(σ_{n+1}, π_{i,n+1})); measured in returnmap.py: −1.2274809384 both ways.
6. **Peak identity** (stated exactly): whenever π_i = π_i*, D = χψ_i and H = 0 (by construction of π_i*; sympy-checked
   D(η*) = χψ_i). H = 0 is the **stationary point of the yield-surface size**. It coincides with the peak of η only on a
   constant-p path: dη/dt = (∂η/∂p)ṗ + (∂η/∂π_i)π̇_i, so on a conventional drained triaxial path (ṗ ≠ 0) the η-peak
   precedes H = 0 by the term (∂η/∂p)ṗ [I]. Test it as: at the step where H changes sign, |D − χψ_i| ≤ tol.
7. **Undrained critical state**: isochoric ⇒ v (hence e) constant; critical state ⇔ ψ_i = 0, π_i = π_i* = p, D = 0,
   H = 0, η = M (θ = π/3): p_cs = −p_a ((e₀ − e)/λ_c)^{1/ξ} (fork; requires e < e₀), p_cs = −exp((v_{c0} − v)/λ̃)
   (paper); q_cs = M|p_cs|/ζ(θ) = M|p_cs| in TXC. At the CS, ε̇^p_v = 0 ⇒ ε̇^e_v = 0 ⇒ p stationary: a fixed point.
8. **Drained CS asymptote**: at large ε_s, ψ_i → 0, π_i → p, η → M (i.e. −ζ(θ)q/p → M, so q/|p| → M/ζ(θ)), D → 0, H → 0.
9. **Dissipation**: D^p = Δλ Σ_a σ_a q_a ≥ 0 every plastic step (0 on elastic ones) under (S.39); a set violating
   (S.39) or N̄ > N is refused; ρ > ρ̄ warned (§11.2).

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
the step-10 value [I]). Fig 6: the minimum for ρ = 0.7 is at φ = π/2 (n in the n¹–n³ plane) [E].

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
and v₀ = 1.59 with v = v₀J.

**K2 gate (plan §5.2, decided at G0)**: the gate is the **ordering** (ρ = 0.7 localizes before ρ = 1) and the **gap**
(n₁.₀ − n₀.₇ ≈ 4 steps) inside a stated band, with a reported sensitivity table over the swept unknowns; exact n is a
sanity check, not a gate. Band [I]: with π_{i,0} ∈ {−60, −80, −100} kPa, χ ∈ {−3.0, −3.5, −4.0}, v_{c0} ∈ {1.80, 1.81,
1.82} and both crossing criteria, require for every combination: n₀.₇ < n₁.₀, gap ∈ [2, 6]; and for the nominal
combination (π_{i,0} = −60.4, χ = −3.5, v_{c0} = 1.81, first-step criterion) n₀.₇ ∈ [19, 25], n₁.₀ ∈ [23, 29]. Report
the full table; a miss outside the band is a finding against the sheet or the oracle, not something to tune away.
The small-strain kernel reproduces this only through the LogStrain wrapper with (S.34); the oracles can run it directly
in finite strain (ε^e_a = ln λ^e_a, v = v₀J, τ). Small-strain v = v₀(1 + ln J) vs v₀J differs by O(10⁻⁴) here.

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
| c = M_e/M_c = 0.71 | ρ | **direct**: ζ(0) = 1/ρ ⇒ M_e = ρ M_c ⇒ ρ = c = 0.71. WW admissible (≥ ½); GA refused (< 7/9) |
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
- The nested π_i Newton can have multiple roots for large Δλ with the cap on (Ω depends on π_i); start from π_{i,n}
  and keep local steps bounded (plan's bounded-work rule). Observed only at absurd trial overshoots in the checks.
- Undamped local Newton diverged once in the cap checks at an extreme trial state; the kernel needs the usual step
  control (the checklist's bounded local work). Not a formula issue.
- (S.34) and (S.44) are transcribed, not re-derived; O1/O2 should FD-check the finite-strain acoustic tensor before K2.

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

### 16.5 G0 review record (2026-09-30)
- **Adversary (Fable, independent re-derivation): PASS-WITH-FIXES**, seven items, all applied in this revision:
  (1) §3.2 hydrostatic guard was self-contradictory → vertex rule Ω := 0; (2) §14 K2 "±1 step" → ordering/gap band with
  sensitivity table; (3) §11.1 parenthetical on N̄ ≤ N corrected (forced by θ = π/3); (4) §3.1 double-precision loss of
  ζ_yy near corners documented, same-θ rule, no isolated unit tests; (5) §9.3 CTO non-symmetric even when associative,
  no symmetric solver; (6) (S.34) repeated-stretch limit of γ̃_ab added; (7) §11.3 BA06 conjugate-term sign oddity noted.
- **Independent numeric check (Opus): PASSED** — dissipation reading A sharp over 1116 cases; 4×4 Jacobian 1.9e−9 vs
  FD; CTO 9.8e−11 (5e−10 on the shear columns).
- **Owner**: approved the refusal rule (hard refuse unless N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄); WARN on ρ > ρ̄).
