---
title: "WP-144 — LadrunoNORSAND: NorSand in the Andrade & Borja (2006) form, with a curved CSL and a hyperelastic energy"
project: Ladruno
type: implementation plan
status: "PLANNED, rev 2. Owner go 2026-09-30 (consistency by construction + full ownership of the model). Plan re-checked against the source papers; P0 (Python oracle) is next. No code yet."
owner: nmora
related:
  - "[[_sand_model_survey_2026-09-27]] (evidence, ranking)"
  - "[[145_ladruno_hysand_plan]] (the cyclic companion)"
  - "[[151_sanisand_reseat_singularity]] (what DM04 + SAS-ME + R1 reaches)"
  - "[[134_sanisand_reference_integrator]] (the oracle template)"
  - "[[128_sanisand_ring_trace]]"
tags: [plan, norsand, critical-state, hyperelasticity, sand, tims, wp-144]
updated: 2026-09-30
---

# WP-144 — LadrunoNORSAND

**Status: PLANNED (rev 2).** No code yet. The ND tags are reserved in `LEDGER_implementations.md`
(**33023** base, **33024** 3D, **33025** PlaneStrain). They go into `classTags.h` only when the code lands
(the ADR-90 precedent). Note for that commit: 33023 and 33025 already exist in the **EigenSOLVER / EigenSOE**
families (FeastEigenSolver, LadrunoCMS). `ci/check_classtags.py` checks collisions per family, so this is not a
collision, but the `#define` comments must say so ("different registry, deliberately NOT a collision"), as for
33020/33021.

Rev 1 (2026-09-27) was written from the survey. Rev 2 (2026-09-30) re-checks every design item against the
two source papers, read in full text:
- **[AB06]** Andrade & Borja (2006), IJNME 67:1531, the 3-invariant version.
- **[BA06]** Borja & Andrade (2006), CMAME 195:5115, the 2-invariant base.

Tags: **[E]** read in the source; **[I]** our inference; **[R]** recollection, still to be read.

## 1. Why, and the start condition

**Owner decision, 2026-09-30:** build it. The reason is *a model that is consistent by construction for design
use, and one we fully control*. It is not "DM04 cannot reach a limit load". WP-151 R1 has since taken DM04 +
SAS-ME past the old wall. With the full set on the B/8 deck it reaches s/B 0.054 with no `loadingNonPosH`; B/4
reaches 0.127 and stops on low confinement. So rev 1's "start after WP-138" gate is superseded.

What NorSand buys, against DM04 + R1:
- A **hyperelastic** energy, so an elastic closed loop does no work **[E, BA06 §2.1]**.
- **Proven non-negative plastic dissipation** under stated parameter bounds **[E, AB06 §2.2.1; BA06 §2.2]**.
- A small implicit local problem with a **closed-form consistent tangent** **[E, AB06 §3]**.
- State = (stress, e, π_i): **no back-stress, no α_in memory, no h ∝ 1/((α−α_in):n) singularity**. Those are the
  roots of every SANISAND integration defect (128, 134, 151) **[I]**.
- A code base we write from the equations, with an oracle, so every term is ours to audit and change.

**Honest limits (keep these in every report):**
- **Admissible, not variational.** Flow is non-associated, so the tangent is non-symmetric (unsymmetric solver).
- **Admissible up to a stated approximation.** BA06 show that the plastic free energy depends on εᵉ through π_i
  and drop the O(Δε_sᵖ²) coupling term in the stress. They justify it as small before the critical state
  **[E, BA06 §2.5]**. We inherit that approximation and say so.
- **Monotonic model.** It has no fabric and weak cyclic behaviour. DM04 stays the cyclic/SSI model; HySand
  (WP-145) is the long-term consistent cyclic candidate.
- **p′ → 0 is not solved by the model.** π_i* ∝ p, so the free-surface ring still needs the floor rule (§2.7).
  The B/4 stop in WP-151 is exactly this failure class.

## 2. Design

### 2.1 Source formulation: AB06, not BA06

BA06 is 2-invariant; the third invariant is explicitly left for later **[E, BA06 closure]**. **AB06 is the
3-invariant version** and is the primary source **[E]**:
- Yield: F = ζ(θ, ρ)·q + p·η(p, π_i), with η = (M/N)[1 − (1−N)(p/π_i)^(N/(1−N))] for N > 0 (eqs. 9–10).
- Plastic potential: Q = ζ(θ, ρ̄)·q + p·η̄(p, π̄_i), with N̄ and ρ̄ (eqs. 13–14). Non-associativity comes from
  both the volume (N̄) and the deviatoric plane (ρ̄).
- Flow in principal directions; closed-form first and second derivatives of Q (eqs. 17–28).
- BA06 is kept as the 2-invariant special case (ζ ≡ 1) and as the oracle's first regression target
  (AB06 Remark 2).

### 2.2 Dissipation bounds: now parameter checks

AB06 §2.2.1 proves D = σ:ε̇ᵖ ≥ 0 at yield if
- **N̄ ≤ N** (from BA06), and
- **ρ ≤ ρ̄**, which is equivalent to ψ_c ≤ φ_c: the dilation angle at critical is at most the friction angle at
  critical (eqs. 32–34) **[E]**.

BA06 §2.2 adds the hardening part: with π_i conjugate to ε_sᵖ, σ:ε̇ᵖ − π_i ε̇_sᵖ ≥ 0, since ε̇_sᵖ = λ̇ ≥ 0 and
π_i < 0 **[E]**.

The parser **refuses** parameter sets outside these bounds (a hard error, not a warning).

What is still open (P0): **do the bounds survive the curved CSL (§2.3)?** The CSL enters only through π_i*
(the limit image pressure) and ψ_i, not through F or Q at fixed π_i. So I expect the proof to carry over
unchanged **[I]**. The oracle's D ≥ 0 census at every step is the evidence, not this expectation.

### 2.3 Critical state line: power law, DM04's form

Both papers use the log CSL, with ψ_i = v − v_c0 + λ̃ ln(−π_i) evaluated **at π_i** (AB06 eq. 41) **[E]**.
Because π_i* ∝ p (eq. 42), the log still diverges as p → 0 **[I]**. Adopt
e_c(π) = e0 − λc (−π/pa)^ξ, with **ψ_i = e − e_c(π_i)**. e0, λc and ξ transfer from TIMs' DM04 calibration
(0.83, 0.027, 0.45).

Consequence: every ∂ψ_i/∂π_i term in AB06's tangent changes (eqs. 54–57: the λ̃/π_i factor becomes
λc ξ (−π_i/pa)^ξ / π_i). The P0 oracle derives these terms symbolically and checks them against finite
differences before any C++ is written.

### 2.4 Lode dependence: Willam–Warnke, because TIMs' c = 0.71

AB06 offers two ζ(θ, ρ) **[E]**:
- **Gudehus–Argyris** (eq. 11): convex **only for 7/9 ≤ ρ ≤ 1**.
- **Willam–Warnke** (eq. 12): smooth and convex for **½ ≤ ρ ≤ 1**.

TIMs' DM04 calibration uses **c = 0.71 < 7/9** (WP-151 §2.5). Gudehus–Argyris would be non-convex there.
**Default to Willam–Warnke.** Offer Gudehus–Argyris as an option, and refuse it for ρ < 7/9.

A C² alternative with Mohr–Coulomb-like sections is Mortara (2026, Acta Mech. 237:1399). It ships FORTRAN and
Octave convexity checks. Keep it as an option for later, not P1.

### 2.5 Elastic energy: Houlsby–Amorosi–Rojas, not BA06's

BA06's energy (Houlsby 1985 type, eqs. 2.2–2.4) gives a bulk modulus **∝ p** (exponent 1). A shear modulus
μ = μ0 + (α0/κ̂)·p0·exp(ω) is also ∝ p. BA06's own runs used **α0 = 0, μ0 = 2000 kPa: a constant shear
modulus** **[E, BA06 Table 1]**, which also sidesteps the energy's coupling term.

TIMs' DM04 has G ∝ √p. To transfer it, use the **Houlsby, Amorosi & Rojas (2005)** energy with n ≈ 0.5
**[R: read the paper before P0]**.
- The dissipation proof (§2.2) does not depend on the choice of elastic energy **[I]**: it uses only F, Q and
  the flow rule at yield.
- HAR-type energies lose convexity above a stress ratio that depends on n and the stiffness ratio
  **[R, survey "could not verify"]**. The deck's ring reaches η ≈ 2.1.
- **P0 gate:** tabulate the Hessian's minimum eigenvalue over the (p, η) range actually visited by the WP-138/151
  footing's Gauss points. If it goes non-positive in that range, go back to BA06's energy with α0 = 0 (constant
  shear modulus) and document the loss of √p shear stiffness.

### 2.6 Hardening

AB06/BA06 use π̇_i = h (π_i* − π_i) ε̇_sᵖ with a **constant h** (h = 280 in both papers) **[E, AB06 eq. 43]**.
Jefferies' original H = H0 − Hy·ψ is a different law.
- **P1 implements the constant-h law.** It is the one whose tangent AB06 derive.
- Hy is a later, separately gated option (the tangent terms of eqs. 53–58 change).
- The P3 calibration refits **h** (not H0/Hy).

The limit image pressure π_i* comes from the maximum dilatancy D* = χ ψ_i through eq. 42 (N̄ > 0 branch).
A&B use χ ≈ −3.5 for sands **[E]**.

### 2.7 Potential cap and the low-confinement rule

- **Q-cap.** The yield function and Q have a corner at the compression end of the p axis, and η can go
  negative. BA06 add a cap on Q, placed by a parameter (e.g. 0.10), and allow a smooth cap instead **[E, BA06
  §2.6 Remark]**. **Use the smooth cap:** a planar cap puts a corner in the consistent tangent.
- **p′ floor.** A documented floor (Pmin ≤ 0.5 kPa), projected and **counted, not refused**. Report the limit
  load at floor F and F/2, and accept the floor if the load changes by less than about 2 % (survey §7.1).
  Report the number of Gauss points at the floor at the limit state.

### 2.8 Return map and tangent: AB06 Box 2, in principal elastic strains

**Local unknowns** x = (ε₁ᵉ, ε₂ᵉ, ε₃ᵉ, Δλ): a **4×4** Jacobian (AB06 eqs. 48–49). The residual is
εₐᵉ − εₐᵉ,tr + Δλ·qₐ = 0 for a = 1, 2, 3, plus F = 0.

A **nested scalar Newton** solves for π_i (eq. 61) inside each local iterate, because π_i depends on ψ_i, which
depends on π_i. BA06 note that folding it in as a 5th unknown costs "about the same" **[E]**. Start nested, as
the paper does, so the oracle and the paper match term by term.

The **consistent tangent** is closed-form in spectral form (AB06 §3). It is **non-symmetric** when N̄ ≠ N or
ρ̄ ≠ ρ.

This supersedes rev 1's "3 local unknowns": that is the 2-invariant BA06 algorithm.

**Finite strain.** AB06 are multiplicative and return in principal elastic **log** stretches; the algebra is
identical to the small-strain return in principal elastic strains **[E, AB06 §3; BA06 §3]**. So: implement
the **small-strain** kernel, and get finite strain by wrapping it in the fork's existing **LogStrain** ND
wrapper. That is a test we can run (§3 P2), not new finite-strain code.

### 2.9 Class shape

- A standalone class in the `LadrunoJ2` style, not the ASDPlasticMaterial3D kit (the kit has no hyperelastic
  energy or ψ-driven image hardening).
- 3D and PlaneStrain wrappers, like the SANISAND family.
- Header-only, OpenSees-free **kernel** + a thin NDMaterial shell (the kernel-oracle doctrine).
- Refusal through the WP-99 commit latch from day one, with refusal codes the element can forward (F7 roster),
  and bounded local work (the material checklist's rules).
- B-bar: BA06 compared standard integration with B-bar near the critical state and saw little difference in
  their plane-strain runs **[E]**. Record the element pairing used in the strip deck.

### 2.10 Optional P4: nonlocal ψ regularisation

ψ-softening localises, so post-peak results depend on the mesh. Mallikarachchi & Soga (2020, Comput. Geotech.
124:103572) regularise NorSand by **nonlocal averaging of the void ratio** (a Galavi–Schweiger weighting plus
softening scaling). They report mesh-consistent force–displacement curves and band widths in drained biaxial
compression **[E, review abstract and conclusions]**. It fits NorSand naturally: ψ is its only softening state
variable.

Scope decision: P1–P3 report the **peak** as the design quantity. B/8 and B/16 show the post-peak dependence;
P4 is run only if TIMs needs the post-peak branch.

## 3. Plan

| Phase | Content | Effort |
|---|---|---|
| **P0** | Python oracle of AB06 Box 2 with the power-law CSL, Willam–Warnke ζ, the HAR energy and the smooth Q-cap. Symbolic derivatives; FD check of every derivative and of the tangent. D ≥ 0 census at every step. The HAR convexity table over the footing's (p, η) cloud (§2.5 gate). Regression: ζ ≡ 1 + log CSL + BA06 energy reproduces BA06's single-point curves. | ~1 week |
| **P1** | Header-only kernel + `LadrunoNorSand` NDMaterial + 3D/PlaneStrain wrappers; 4×4 spectral return + nested π_i; closed-form tangent; parser with the §2.2 bound refusals; `sendSelf`/`recvSelf`; responses (ψ, ψ_i, π_i, D, dissipation increment, floor and refusal counts); commit-latch refusal | ~2 weeks |
| **P2** | Tests: kernel parity vs the oracle; FD tangent; D ≥ 0 census; LogStrain-wrapped finite-strain parity vs the oracle at large strain; byte-identity of everything else; refusal on discarding elements; bounded work; mutation gate | ~1 week |
| **P3** | Calibration from TIMs' data: transfer the CSL (e0, λc, ξ), Mc = M_tc 1.3309, ρ from c = 0.71, the HAR energy from G0/√p and ν, and the initial e. Refit χ from peak dilatancy vs ψ, h from pre-peak stiffness and strain to peak, and N, N̄, ρ̄ from the volumetric curves. Single-point gates vs DM04 (oracle `uw_model`) and PM4Sand. Strip deck at B/8 and B/16 with the floor-sensitivity report. | ~1 week |
| **P4** (optional) | Nonlocal ψ (§2.10) | ~1–1.5 weeks |
| **Gate** | Adversarial review (new maths); banner line; ledgers; the material guide | — |

**Total: about 5–6 engineer-weeks without P4.** P0 is shorter than rev 1 estimated, because the 3-invariant
proof exists (AB06 §2.2.1). P1 is longer, because the return is the 4-unknown spectral one.

## 4. Sources

Primary (read in full text, 2026-09-30):
- **[AB06]** Andrade, J. E. & Borja, R. I. (2006). Capturing strain localization in dense sands with random
  density. *IJNME* 67:1531–1564. doi:10.1002/nme.1673. §2 (constitutive), §2.2.1 (dissipation), §3 + Boxes 1–2
  (return map, CTO).
- **[BA06]** Borja, R. I. & Andrade, J. E. (2006). Critical state plasticity, Part VI. *CMAME* 195:5115–5140.
  Author PDF: https://geomechanics.civil.northwestern.edu/Papers_files/ccvi.pdf. §2.1 (energy), §2.2
  (dissipation), §2.5 (free-energy coupling), §2.6 (return map, Q-cap), Tables 1–2.

Supporting (the owner's local library `SOILS_rev`, reviews read):
- Mallikarachchi, H. & Soga, K. (2020). *Comput. Geotech.* 124:103572 (nonlocal NorSand, §2.10).
- Jefferies, M. (2022). On the fundamental nature of the state parameter. *Géotechnique* 72(12):1082.
- Jefferies, M., Shuttle, D. & Been, K. (2015). Principal stress rotation as cause of cyclic mobility.
  *Geotech. Res.* 2(2):66 (the NorSand-PSR route, out of scope).
- Mortara, G. (2026). Lode dependence incorporating the Mohr–Coulomb deviatoric section. *Acta Mech.*
  237:1399 (C² ζ option).
- Sheng, Sloan & Yu (2000). *Comput. Mech.* 26:185 (critical-state implementation practice).

To read before P0 **[R]**:
- Houlsby, G. T., Amorosi, A. & Rojas, E. (2005). *Géotechnique* 55(5):383 (the HAR energy and its convexity
  limit).
- Borja, R. I., Tamagnini, C. & Amorosi, A. (1997). *JGGE* 123(10) (convexity of pressure-dependent
  hyperelasticity).
- Jefferies, M. (1993). *Géotechnique* 43(1); Jefferies & Shuttle (2002). *Géotechnique* 52(9) (the original
  NorSand and its 3-invariant form, to map parameter names).
- The survey ([[_sand_model_survey_2026-09-27]]) holds the rest of the citations.
