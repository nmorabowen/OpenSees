---
title: "WP-144 — LadrunoNORSAND: NorSand in the Andrade & Borja (2006) form, with a curved CSL and a hyperelastic energy"
project: Ladruno
type: implementation plan
status: "PLANNED, rev 2.1. Owner go 2026-09-30 (consistency by construction + full ownership). Rev 2.1 adds the oracle / known-result gates (§5) and the orchestration roster (§6). P0a (equation sheet) is next. No code yet."
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

AB06 §2.2.1 states D = σ:ε̇ᵖ ≥ 0 at yield if N̄ ≤ N and ρ ≤ ρ̄ (eqs. 32–34) **[E]**.

**That condition is wrong for the flow rule AB06 actually integrate** (G0, 2026-09-30) **[I]**.
- The paper's flow rule, Box 2 and tangent use ∂Q/∂p = β ∂F/∂p with β = (1−N)/(1−N̄) (eq. 22, "reading A").
- The proof (eqs. 30–31) uses a different potential with η̄ = (ζ̄/ζ)η ("reading B").
- With reading A, D/λ̇ on the yield surface is linear in r = (p/π_i)^{N/(1−N)} and binds at r = 0. That
  gives β ≤ ζ̄/ζ for all θ.
- For Willam–Warnke, ζ̄/ζ is monotone in θ, so the bound is set at the corners. The extension corner gives
  ζ̄/ζ = ρ/ρ̄.

**Adopted condition** (owner decision 2026-09-30):
- **N̄ ≤ N and ρ/ρ̄ ≥ (1−N)/(1−N̄).**
- **Hard refusal** in the parser if it fails.
- **Warning only** for ρ > ρ̄ (the paper's ψ_c ≤ φ_c reading).
- Counterexample to the paper's test: N̄ = N, ρ = 0.7, ρ̄ = 0.8 passes ρ ≤ ρ̄ yet dissipates negatively at the
  extension corner at small |p|.
- Evidence: the equation sheet [[144a_norsand_equation_sheet]] §5.2/§11. It was confirmed by the Adversary's
  independent re-derivation (26×26 (ρ, ρ̄) grid × 2001 θ, no interior excess) and by the orchestrator's own
  check.

BA06 §2.2 adds the hardening part: with π_i conjugate to ε_sᵖ, σ:ε̇ᵖ − π_i ε̇_sᵖ ≥ 0, since ε̇_sᵖ = λ̇ ≥ 0 and
π_i < 0 **[E]**.

**Curved CSL: settled at G0.** The dissipation argument does not involve the CSL. Only Λ(π_i) = π_i ∂ψ_i/∂π_i
changes, in the tangent (sheet §11.4) **[I]**. The oracle's D ≥ 0 census at every step remains the evidence.

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
- **Willam–Warnke** (eq. 12): smooth and convex for **½ < ρ ≤ 1**. At ρ = ½ exactly the compression corner becomes a vertex (ζ′(π/3) = −√3, found by O1 at G1). The parser **refuses ρ = ½ and ρ̄ = ½** (owner decision 2026-10-01). The range applies to ρ̄ as well.

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

## 5. Oracles and known results (what "correct" means, fixed before any code)

**Principle.** Every test compares against something that was **not produced by the code under test**: a closed
form, a published number, a lab curve, or an oracle written independently of it. The expected value and its
tolerance are written in the test file before the implementation exists.

### 5.1 Two independent oracles (P0)

| Oracle | What it is | Written from | Its job |
|---|---|---|---|
| **O1, rate oracle** | The continuum rate equations (AB06 Box 1 + the curved CSL + the chosen energy), integrated as an ODE by SciPy Radau at rtol 1e-10. It uses the closed-form plastic multiplier; no return map. This is the WP-134 pattern. | The equation sheet (§6, P0a) | The **truth** for any strain path. It is independent of all the return-map algebra. |
| **O2, algorithmic oracle** | AB06 Box 2 in Python: backward Euler, the 4×4 spectral return, nested π_i, and the closed-form tangent. Derivatives are generated with sympy, not typed by hand. | The equation sheet, **by a different agent than O1** | The **reference for the C++**, which must match it to ~1e-10 relative, because it is the same algorithm in the same arithmetic. |

Gates between them:
- O2 → O1 with first-order convergence as Δε → 0, on every K-path below.
- O2's tangent matches finite differences to ~1e-7 (central differences, several step sizes).
- A mismatch between O1 and O2 is a finding against the equation sheet or one of the two codes. It is never
  resolved by editing one oracle to agree with the other.

### 5.2 Known results (determinate expectations)

**K1, closed-form identities** (exact up to tolerance; T0m in the manifest):
1. Hyperelastic isotropic compression: p(ε_v) matches the energy's closed form.
2. A closed elastic strain loop does **zero net work** (≤ 1e-12) and returns the state exactly (conservative
   energy, BA06 §2.1).
3. ζ(θ): ζ = 1 at the compression corner and 1/ρ at the extension corner. The Willam–Warnke section is convex
   for ρ ∈ (½, 1] (½ refused, §2.4), and the parser refuses outside that range, for ρ and ρ̄.
4. The yield function gives η = M·ζ-scaled at p = π_i (the image-stress definition, AB06 eq. 10).
5. Flow rule: at every plastic step the measured ε̇_vᵖ/ε̇_sᵖ equals AB06 eq. 39 evaluated at the state.
6. **Peak identity:** at a drained peak, π_i = π_i* and D = χ ψ_i.
7. **Undrained critical state endpoint:** an isochoric triaxial ends on the CSL at
   p_cs = −pa·((e0 − e)/λc)^(1/ξ) with q = M_tc·|p_cs| (closed form from the power-law CSL).
8. **Drained critical state:** at large shear strain, ψ → 0, η → M(θ) and D → 0.
9. **Dissipation:** D ≥ 0 at every plastic step, and D = 0 at every elastic step. A parameter set with N̄ > N or
   ρ/ρ̄ < (1−N)/(1−N̄) is refused (§2.2). The counterexample set N̄ = N, ρ = 0.7, ρ̄ = 0.8 is a test: refused
   by the parser, and negative D in the oracle when forced.

**K2, published benchmark (T1): AB06 §6.1, the single-point localization test.**
- Parameters, all given in the paper **[E]**: κ̂ 0.01, ε_v0 0 at p0 −100 kPa, μ0 5400 kPa, α0 0; λ̃ 0.0135,
  M 1.2, N 0.4, N̄ 0.2, h 280; v 1.59, v_c0 1.81; Willam–Warnke.
- Loading: f₁ for n₁ = 10 steps, then f₂ until localization.
- **Published result:** ρ = 0.7 / ρ̄ = 0.8 localizes at **n = 22**; ρ = ρ̄ = 1 at **n = 26**.
- Loading read from the PDF (P0a): λ₁ = 1e-3, λ₂ = 4e-4, f₁ = diag(1+λ₂, 1−λ₁, 1), f₂ = diag(1, 1−λ₂, 1+λ₁).
- **Gate form (set at G0):** the paper leaves the initial π_i, χ, the exact v_c0 and the crossing criterion
  unspecified, and each can move n by more than 1.
  - **The gate is the ordering:** ρ = 0.7 localizes before ρ = 1, with a gap of about 4 steps, inside a stated
    band.
  - A sensitivity table over π_i0 and χ is reported.
  - Exact n = 22 / 26 is a sanity check, not a pass/fail.
- It needs finite strain (LogStrain) and a min det(n·A·n) sweep over directions (done in the oracle).
- This case also exercises the log CSL and BA06's energy: it is the "paper mode" regression before the fork's
  extensions are switched on.

**K3, laboratory data (T1, fit quality reported, sanity-gated):**
- Ottawa F65 monotonic drained triaxials (Vasko 2014 / LEAP-2017, in the owner's library).
- TIMs' own drained triaxials (P3).
- Gate: the calibrated curves stay inside the specimen scatter to the peak. The residual is reported, not hidden.

**K4, cross-model (soft gate):** single-point drained/undrained triaxial and plane-strain compression against
the DM04 oracle (`uw_model`, WP-134) and PM4Sand, with the same CSL. Peak q and ε_v must agree within the
triaxial scatter; differences are explained, not tuned away.

**K5, boundary-value problem (P3):**
- The strip deck at B/8 and B/16.
- Floor sensitivity: F vs F/2, accepted if the load changes by less than 2 %.
- **Loukidis & Salgado (2011), Géotechnique 61(2):107**, in the owner's library: Nγ as a function of relative
  density and stress level, at the deck's density and stress. Also cross-check against the Lau (2011) and
  Lyamin (2007) bearing-capacity results.
- These are the design-level known results. A miss here is reported, with the mechanism, not calibrated away.

### 5.3 How the gates map to the manifest

- **T0m** = K1 plus the kernel-vs-O2 parity.
- **T1** = K2, K3 and FD-tangent.
- The mutation gate targets semantic mutants: drop a tangent term, flip N̄/N, freeze π_i, remove the Q-cap,
  bypass the floor count, swap WW for Gudehus–Argyris. Every mutant must be killed by a named K-test. A
  survivor is a coverage gap to be recorded, as in the ADR-92 gate.

## 6. Execution and orchestration

The orchestrator (the lead session) does not write the maths or the kernel. It owns:
- the equation sheet's sign-off;
- the task specs;
- the independence rule between oracles;
- builds and CI;
- the ledgers;
- the decision at every gate.

Agents get narrow specs and return short reports: at most ~300 words plus file paths. The orchestrator reads
diffs of the load-bearing parts, not whole files.

### 6.1 Roster: model and effort per role

| Role | Model | Effort | Why this tier |
|---|---|---|---|
| **Deriver**: equation sheet, symbolic derivatives, the curved-CSL proof check | Fable | high | Dense maths where a sign error survives everything downstream. No downgrade. |
| **Oracle author A** (O1) | Opus | high | Numerics, written from the sheet only |
| **Oracle author B** (O2) | Fable | high | Deliberately a different model from O1, to decorrelate errors |
| **Kernel author**: header-only C++ kernel + tangent | Opus | high | Must match O2 to 1e-10 |
| **Shell/wiring**: NDMaterial shell, wrappers, parser, dispatch, sendSelf/recvSelf, responses | Sonnet | medium | Pattern work against the `ladruno-new-material` checklist and the SANISAND wrappers |
| **Test author**: K1/K2 pytest, mutation mutants | Sonnet | medium | The expected values come from the spec, not from running the code |
| **Adversary**: derivation and C++ review at the P0 and final gates | Fable (whole artifact) + an independent **Opus** numeric check (narrow, own code) | max / high | Two independent reviewers on new maths (the plan's adversarial-gate requirement), from different models than the author. No Codex/ChatGPT (owner policy, 2026-09-30). |
| **Third-family review** (optional, G0 and final gate) | Grok, run by the owner in Cursor | — | A model outside the Claude family. The orchestrator writes a self-contained prompt packet; the owner pastes it in Cursor and relays the report. There is no headless Cursor CLI on this machine. |
| **Calibrator**: K3/K4/K5, P3 | Opus | high | Judgment on data fit and mechanisms |
| **Clerk**: ledgers, banner, manifest row, grep sweeps | Haiku | low | Mechanical, verified by CI gates |

Model and effort are set per role by agent definitions in `.claude/agents/` (session-local, created at kickoff).

### 6.2 Sequence and gates

1. **P0a, equation sheet** (Deriver).
   - Writes `Ladruno_implementation/144a_norsand_equation_sheet.md`: AB06/BA06 equations re-derived in our
     notation, with the curved CSL, the WW ζ and its derivatives, the Q-cap, the energy module, and the AB06 §6.1
     eq. 98 values read from the PDF page.
   - The sheet is the **only** document later agents read. Nobody re-reads the PDFs; that saves tokens.
   - **Gate G0:** Adversary sign-off plus orchestrator review.
2. **P0b ∥ P0c**: O1 (author A) and O2 (author B) are written in parallel, from the sheet only. **P0d** (Test
   author) writes the K1/K2 tests against **expected values from the sheet and the paper**. It runs them against O1
   and O2.
   - **Gate G1:** O2→O1 convergence, FD tangent, K1 all green, K2 ordering and gap reproduced (with the π_i0/χ sensitivity table), D ≥ 0 census,
     convexity table (§2.5).
3. **P1**: the kernel (Kernel author), then the shell and wiring (Shell/wiring). Builds use `build.bat` in this
   worktree, run in the background.
   - **Gate G2:** kernel-vs-O2 parity at ≤ 1e-10, K1/K2 through OpenSees, LogStrain large-strain parity.
4. **P2**: the mutation gate, refusal and bounded-work tests, byte-identity of everything else.
   - **Gate G3:** Zone-A green and a mutation score at or above the floor.
5. **P3**: K3/K4/K5 (Calibrator). The strip runs need Esmeralda (owner or remote orchestrator).
6. **Final gate**: Adversary pass on the whole diff; Clerk does the ledgers, banner and manifest. Flip to ready.
   The owner merges.

### 6.3 Token discipline (without trading quality)

- **One read of the sources.** The PDFs are read once, into the equation sheet (P0a). Every later agent reads
  the sheet (a few kB), not the papers.
- **Templates are reused, not rediscovered.** The agents' specs point at WP-134's oracle, `LadrunoJ2`'s kernel
  shape, the SANISAND wrappers and the `ladruno-new-material` skill by path.
- **Bounded loops.** At most two fix-and-review rounds per gate before the orchestrator escalates to the owner.
- **Cheap tiers only where CI verifies the output** (the Clerk). Never on maths, oracles or the kernel.
- **Builds and long runs go to the background**, with no polling.
- **Findings, not transcripts.** Agent reports list findings and paths. Evidence stays in the files.

### 6.4 Blockers

- **Houlsby, Amorosi & Rojas (2005)** is not in the library.
  - Mitigation: the energy is a swappable module. P0 starts with BA06's energy (which the K2 benchmark needs
    anyway). HAR is added when the paper arrives, behind its convexity gate.
- **Esmeralda access** for K5 strip runs: owner or remote orchestrator, as in WP-151.

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
