# -*- coding: utf-8 -*-
"""Round-3 edits (p-floor, pi_i0 rule, energy bookkeeping, scan-step correction) of the WP-144a equation sheet.

Each edit is an exact-anchor replacement that must match ONCE; a drifted anchor aborts before anything is written.
    python apply_sheet_round3.py <sheet.md> [<out.md>]        (no <out.md>: edit in place)
Written this way because the Deriver session was bound to another worktree and could not edit the sheet directly.
"""
import sys

SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else SRC
text = open(SRC, encoding="utf-8").read()
EDITS = []


def rep(old, new):
    EDITS.append((old, new))


# ---------------------------------------------------------------------------------------------- frontmatter
rep('status: "HAR energy option derived and gated 2026-10-02 (',
    'status: "p-floor round 2026-10-03 (§16.6 item 16: §9.7 floor operator Π_f — a strain-space projection at fixed '
    'deviatoric elastic strain onto p = −p_min, applied to the trial and to the converged state, counted, never '
    'refuses; (S.48)–(S.55) with (S.32f)/(S.46f) FD-checked to ≤ 5e−8 and every named mutant ≥ 2.5e−3 apart; §5.4 '
    'unified π_i0 rule (S.53); §2.4 energy-option bookkeeping; §8/§10.2 nested-scan step corrected from 10⁻⁴ to the '
    'shipped 10⁻³ with the ramp-width criterion (S.56), O2 and kernel unchanged; scripts in '
    'Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_floor/). HAR energy option derived and gated 2026-10-02 (')
rep("updated: 2026-10-02\n", "updated: 2026-10-03\n")

# ---------------------------------------------------------------------------------------------- §1.3 row
rep("enters **no** derivative since G2, §1.2) | BA06 Box 2; §1.2 |\n",
    "enters **no** derivative since G2, §1.2) | BA06 Box 2; §1.2 |\n"
    "| p_min | p′ floor (§9.7): every trial and committed state has p ≤ −p_min, by the projection Π_f (S.48); ≥ 0, "
    "0 = off (the pre-round-3 refusals). Default 5·10⁻³ p_ref with p_ref := |p₀| (BA06) or p_a (HAR): 0.5 / 0.505 kPa | "
    "plan §2.7; §9.7 |\n")

# ---------------------------------------------------------------------------------------------- §2.3 wording
rep("a tensile mirror image that the kernel must refuse (the p-floor rule of plan §2.7 handles the approach).",
    "a tensile mirror image that is never evaluated: a trial strain with ε* ≤ 0 is a floor event (§9.7: Π_f of (S.48) "
    "acts in strain space before any stress is formed; with p_min = 0 it is the refusal `elastic_domain`, §2.4), and a "
    "local iterate with ε* ≤ 0 is an evaluation failure like p ≥ 0 (line-search backtrack).")
rep("(4) §9.1: F_tol = 10⁻¹⁰·|p₀| and the r₄/|p₀| scaling (kernel.h:936, 1067) use **p₀ := −p_a** under HAR.",
    "(4) §9.1: F_tol = 10⁻¹⁰·|p₀|, the r₄/|p₀| scaling (kernel.h:936, 1067) and the default p_min = 5·10⁻³|p₀| (§9.7) "
    "use **p₀ := −p_a** under HAR.")
rep("parser validation (k, g > 0, 0 ≤ n < 1, p_a > 0; refuse ε* ≤ 0 as a p-floor hit), O2 `elastic`/`energy_psi`/`invert`,\n"
    "O1 `model.py` energy.\n",
    "parser validation (§2.4; ε* ≤ 0 at the trial is a floor event, §9.7), O2 `elastic`/`energy_psi`/`invert`,\n"
    "O1 `model.py` energy. The complete branch list is §2.4.\n"
    "\n"
    "### 2.4 Energy-option bookkeeping (`-energy BA06|HAR`; owner decision (a) 2026-10-02) [I]\n"
    "\n"
    "One flag, two parameter sets, and the table below is the complete list of places where a code branches on it.\n"
    "Anything not in it is energy-blind by construction (§2.1: the plastic part sees the energy only through p, q, D,\n"
    "a^e of (S.3), the inverse map and the reference pressure).\n"
    "\n"
    "| quantity | BA06 (default; paper mode; K2/K2b) | HAR | who branches |\n"
    "|---|---|---|---|\n"
    "| parameters | p₀ < 0, κ̂ > 0, ε_{v0}, μ₀ > 0, α₀ (0 in both papers) | k > 0, g > 0, 0 ≤ n < 1, p_a > 0 — and **none** of the five BA06 values | parser; kernel `Params` (+ `energy`, k, g, n_e); O1 `params.py`; O2 `params.py`; `sendSelf`/`recvSelf`; `Print` echo |\n"
    "| reference pressure p_ref | \\|p₀\\| | p_a (p₀ := −p_a) | F_tol = 10⁻¹⁰ p_ref (§9.1); r₄/p_ref (kernel.h:936, 1067; O2 `scaled_norm`); the p_min default 5·10⁻³ p_ref (§9.7); O1's surface tolerance |\n"
    "| p, q, D₁₁, D₁₂, D₂₂ | (S.5); D₁₂ = 0 for α₀ = 0 | (S.5h)/(S.5h'); D₁₂ ≠ 0 always | kernel `elastic()`, O2 `elastic`, O1 `model.py` |\n"
    "| q/ε_s in (S.3) and its ε_s → 0 limit | q/ε_s = D₂₂ identically (α₀ = 0): t3/t4 of (S.3) collapse | q/ε_s = 3g p_a (ϖ/p_a)^n ≠ D₂₂; limit D₂₂\\|_{q=0} | the same three; (S.3) in full |\n"
    "| domain of ε^e | all of ℝ³ (p = p₀e^ω < 0 always) | ε* > 0; a trial with ε* ≤ 0 is a **floor event** (§9.7; with p_min = 0 the refusal `elastic_domain`); a local iterate with ε* ≤ 0 is an evaluation failure (backtrack), like p ≥ 0 | kernel `elastic()` (a failure code, never NaN), O2 `EvalError`, O1 (event stop) |\n"
    "| inverse map (`initialState`) | closed form for α₀ = 0, 2×2 Newton otherwise | (S.5h''), closed form | kernel `invert_elastic`, O2 `invert_elastic`, O1 `initial_state` |\n"
    "| Ψ (K1.2 loop work; the floor energy E_f of (S.52)) | (S.4) | (S.4h) | O2 `energy_psi`, O1; the kernel only if it reports E_f |\n"
    "| floor closed form ε_{v,f}(ε_s), ε'_f | (S.49) | (S.50) | kernel, O2, O1 (the operator Π_f) |\n"
    "| K1 expected values | 13.1, 13.12 | 13.1h, 13.11, 13.13, 13.14 | test author |\n"
    "| parameter transfer (§15) | μ₀, κ̂ refit | g, k from G₀, ν, e_ref; p_a = p_atm | calibrator; echo |\n"
    "\n"
    "Parser refusals (hard; both oracles and the shell). HAR: k ≤ 0, g ≤ 0, n < 0, n ≥ 1 (n = 1 is HAR05 eq 47–48, another\n"
    "closed form, not shipped), p_a ≤ 0; any of `-p0 -kappa_hat -eps_v0 -mu0 -alpha0` given together with `-energy HAR`\n"
    "(refused, never ignored: HAR replaces α₀ and the other four are not read). BA06: any of `-k -g -n` given. Both:\n"
    "p_min < 0; an initial state with p ≥ 0 after the floor (only possible with p_min = 0); ε* ≤ 0 at `initialState` under\n"
    "HAR (the inverse map (S.5h'') needs p < 0). No new warnings.\n"
    "\n"
    "FD checks that must be re-run under HAR before the option ships, because D₁₂ ≠ 0 and D₂₂ ≠ q/ε_s make the t2 and\n"
    "t3/t4 terms of (S.3) live for the first time (under BA06 with α₀ = 0 both were inert): (S.30) Jacobian vs central FD\n"
    "(returnmap.py-type, fork and paper CSL); (S.33) non-coaxial CTO off the corners; (S.34) with (S.44) against the\n"
    "nominal-stress FD (r2_finite_tangent_fd-type); the §9.6 chain cases (A) (B) (C) (E) (F); the floor tangents\n"
    "(S.51)/(S.54) (floor_fd-type, under HAR: the round-3 scratch ran them under BA06 with α₀ = 5 as the D₁₂ ≠ 0 stand-in\n"
    "and the elastic floor operator under HAR symbolically, §16.4); kernel-vs-O2 parity on a HAR path; K1.2 (closed loop)\n"
    "under HAR. The mutant \"HAR silently replaced by BA06\" is killed by 13.11 (η = 3gε_s exactly on a constant-volume\n"
    "elastic shear, against p constant under BA06), by 13.1h and 13.13 (the domain edge and the floor closed form differ\n"
    "from the BA06 values by orders of magnitude), and by any FD-tangent test with q > 0 (the t4 term of (S.3) is zero\n"
    "under BA06 α₀ = 0 and not under HAR).\n")

# ---------------------------------------------------------------------------------------------- §5.4 pi_i0 rule
rep("D < 0 (compaction) for η < M, D > 0 (dilation) for η > M, D = 0 at the image point.\n\n---\n\n## 6. Critical state line: two modes",
    "D < 0 (compaction) for η < M, D > 0 (dilation) for η > M, D = 0 at the image point.\n"
    "\n"
    "### 5.4 Initial image pressure: the unified π_{i0} rule [I; owner decision (c) 2026-10-02; sympy-checked, floor_sympy.py (f)]\n"
    "\n"
    "For both CSL modes and every N: the rule is pure yield-surface geometry, the CSL enters only ψ_{i0} afterwards. Inputs:\n"
    "the initial stress after the floor of §9.7 (p_init < 0, q_init, θ_init), ρ, and the cap's upper bound c₂ (smooth cap:\n"
    "c₂; planar: c₂ = c₁ = χ_cap; no cap: c₂ := 0). With\n"
    "\n"
    "      η_init := ζ(θ_init, ρ) q_init/|p_init|   (0 on the axis, §3.2),       η* := max(η_init, c₂ M),\n"
    "\n"
    "π_{i0} is the inverse of (S.12) at (p_init, η*):\n"
    "\n"
    "      π_{i0} = p_init exp(η*/M − 1)                                        (N = 0),\n"
    "      π_{i0} = p_init [ (1−N)/(1 − η* N/M) ]^{(1−N)/N}                       (N > 0; requires η* < M/N).          (S.53)\n"
    "\n"
    "(Numbered with the round-3 equations of §9.7.) Refuse at `initialState` if η* ≥ M/N: no surface passes through the\n"
    "state. Checked: η(p_init, π_{i0}) = η* exactly in both branches (10 random sets, 30 digits); the N → 0 limit of the\n"
    "power branch is the exp branch (gap O(N): 2.4e−10 at N = 10⁻⁹); η* = M gives π_{i0} = p_init (the image point); η* = 0\n"
    "gives p_init (1−N)^{(1−N)/N}, the apex through p_init (the pre-round-3 O2 default for an isotropic start);\n"
    "d ln|π_{i0}|/dη* = (1−N)/(M − η*N) > 0, a larger η* is a larger surface. Then ψ_{i0} = e − e_c(π_{i0}) by (S.22) (fork)\n"
    "or v − v_{c0} + λ̃ ln(−π_{i0}) (paper), and the B > 0 guard of §7 is checked at (π_{i0}, ψ_{i0}).\n"
    "\n"
    "What the rule does. If η_init ≥ c₂M the surface passes through the initial stress (F = 0 exactly: the K0 deck states,\n"
    "η_init ≈ 0.75 > 0.2). If η_init < c₂M the state is inside, F = |p_init|(η_init − c₂M) < 0, and at p = p_init the\n"
    "surface sits at the ramp's upper end η = c₂M, so the whole cap ramp [c₁M, c₂M] lies inside the initial surface at that\n"
    "pressure. The first yield point depends on the path: a constant-p path yields at η = c₂M exactly (w = 1, cap\n"
    "inactive); a drained TXC path from an isotropic state raises |p| and meets the surface lower — K2 set, p_init = −100\n"
    "kPa, c₂ = 0.15: π_{i0} = −50.9959 kPa, first yield at q = 11.35 kPa, p = −103.78 kPa, η = 0.109, where w = 0.34 (inside\n"
    "the ramp, cap partly active) — against the apex start (η* = 0, π_{i0} = −46.4758 kPa), which is plastic from the\n"
    "first increment at w = 0. That is the `ramp_end` start of the P3 ruling: elastic up to a finite q, no apex\n"
    "substepping. K2 values (K1.15): η* = c₂M → −50.995881; η* = 0 → −46.475800; η_init = 0.75 → −71.554175; η* = M → −100.\n"
    "Limits: p_init → 0⁻ gives π_{i0} ∝ p_init (the floor of §9.7 is applied first, so |p_init| ≥ p_min); c₂M → M/N is the\n"
    "refusal; the no-cap case (c₂ = 0) with an isotropic start is the apex rule, unchanged from G1. Code: O2\n"
    "`pi_of_eta(P, p_init, max(eta_init, c2*M))`; the kernel's `initialState` and O1 `initial_state` must take the same\n"
    "c₂ (0 for cap = none); a `-pi0` given on the deck overrides the rule (the paper-mode K2 benchmark needs\n"
    "π_{i0} = p_c (0.6)^{1.5}, §14).\n"
    "\n---\n\n## 6. Critical state line: two modes")

# ---------------------------------------------------------------------------------------------- §8 scan step
rep("found by a scan from π_{i,n} in the direction of −r(π_{i,n}) with a step much smaller than the dip width of §10.2\n"
    "(≤ 10⁻⁴|π_{i,n}|), refined by safeguarded Newton/bisection inside that bracket; never a factor-2 geometric bracket,\n"
    "which can enclose the far w = 1 root.",
    "found by a scan from π_{i,n} in the direction of −r(π_{i,n}) with the fixed step PI_SCAN_REL·|π_{i,n}|,\n"
    "PI_SCAN_REL = 10⁻³ (range PI_SCAN_MAX·PI_SCAN_REL = |π_{i,n}| with PI_SCAN_MAX = 1000), refined by safeguarded\n"
    "Newton/bisection inside that bracket; never a factor-2 geometric bracket, which can enclose the far w = 1 root. The\n"
    "step is set by the cap-ramp width (S.56), not by the dip width, which vanishes at the fold: the G1 text\n"
    "\"≤ 10⁻⁴|π_{i,n}| (≪ the dip width)\" was the wrong criterion and is withdrawn (round 3, §10.2).")

# ---------------------------------------------------------------------------------------------- §9.1 steps 1f, 5f
rep("1. Trial: ε^{e,tr} = ε^e_n + Δε. Spectral: ε^{e,tr} = Σ_a ε̃_a m^a. (The converged ε^e, σ and ε̇^p share these m^a.)\n",
    "1. Trial: ε^{e,tr} = ε^e_n + Δε. Spectral: ε^{e,tr} = Σ_a ε̃_a m^a. (The converged ε^e, σ and ε̇^p share these m^a.)\n"
    "   1f. **Floor at the trial (§9.7; p_min > 0):** if p(ε̃) > −p_min, or ε̃ is outside the energy's domain (HAR: ε̃* ≤ 0),\n"
    "   replace ε̃ := Π_f(ε̃) of (S.48) (same m^a; counted as a trial floor event). Steps 2–4 see the floored trial.\n")
rep("5. Update: ε^p_{n+1} = ε^p_n + Δλ Σ_a q_a m^a, σ = Σ_a σ_a m^a, state (σ, e or v, π_i).\n",
    "5. Update: ε^p_{n+1} = ε^p_n + Δλ Σ_a q_a m^a, σ = Σ_a σ_a m^a, state (σ, e or v, π_i).\n"
    "   5f. **Floor after convergence (§9.7):** if p(ε^e_{n+1}) > −p_min, replace ε^e_{n+1} := Π_f(ε^e_{n+1}) (counted as a\n"
    "   post floor event); π_i, v, ε^p and D^p of the step are those already formed. The tangent is then (S.32f); a\n"
    "   substepped increment chains with (S.54).\n")

# ---------------------------------------------------------------------------------------------- §9.7 (new)
S97 = r"""
### 9.7 The p′ floor (plan §2.7; owner decision (b) 2026-10-02; `-pmin`) [I; sympy- and FD-checked: floor_sympy.py, floor_fd.py, floor_rate.py, floor_fold.py]

Why. π_i* ∝ p (S.23), K and G → 0 as p → 0 under HAR (K → 0 under BA06), F_pp ∝ 1/p (S.13): at the free-surface ring
(p′ down to 0.19 kPa in the TIMs dumps, 0.57 kPa at s/B 0.002 on the Kimura deck) the model is not wrong but empty, and a
HAR trial strain can leave the energy's domain. The plan's rule: a floor p_min, projected and counted, never refused,
with the F vs F/2 report as the acceptance.

**Options weighed.** (i) A projection of the state onto p ≤ −p_min, after the return and of the trial before it; (ii) a
regularised energy with constant moduli below |p| = p_min; (iii) a clamp of p inside the pressure-dependent terms of the
return map only. (ii) keeps a potential but fixes only the stiffness: F, π_i* and F_pp still collapse at p → 0, the
plastic machinery still needs a guard, and K1.1/K1.1h change below p_min. (iii) is a partial substitution with no
consistent derivative (clamped terms under an unclamped F) and would change the model term by term in silence.
**Adopted: (i), as the operator Π_f below**, applied at two points of every (sub-)increment: the trial (§9.1 step 1f)
and the converged state (step 5f). It never acts inside the local Newton (an iterate with p ≥ 0 or outside the domain
stays an evaluation failure → line-search backtrack, bounded work) and never inside the nested π_i solve.

**The operator (S.48): strain space, deviatoric elastic strain fixed.** With ε_v = tr ε^e, e = dev ε^e, ε_s = √(2/3)‖e‖,
n̂^e = e/‖e‖ (0 at ε_s = 0), let ε_{v,f}(ε_s) be the unique solution of p(ε_v, ε_s) = −p_min at fixed ε_s (unique because
∂p/∂ε_v = D₁₁ > 0, §2):

      Π_f(ε^e) := ε^e − ((ε_v − ε_{v,f}(ε_s))/3) 1    if p(ε^e) > −p_min or ε^e ∉ dom Ψ,      := ε^e otherwise;
      principal: ε_{f,a} = ε_a − Δε^f_v/3,    Δε^f_v := ε_v − ε_{v,f}(ε_s) > 0 whenever active.                    (S.48)

Π_f is co-axial with its argument, preserves every eigenvalue difference (n̂^e, θ and the spectral basis are unchanged),
is idempotent, and is defined for an out-of-domain HAR trial (it needs only ε_s). The alternative "keep q, move p" (the
hydrostatic projection in stress space) coincides with (S.48) under BA06 with α₀ = 0 (q = 3μ₀ε_s depends on ε_s alone)
but needs a stress, which an out-of-domain trial does not have; under HAR the two differ by the D₁₂ coupling (q_f =
3g p_a ε_s (ϖ_f/p_a)^n grows with the floored pressure — the stiffness-floor reading of SANISAND's `-Pmin`), and (S.48)
is adopted so that one operator serves the trial and the converged state in both energies.

Closed forms (sympy-exact, floor_sympy.py (a)–(c)):

      BA06:  ε_{v,f} = ε_{v0} − κ̂ ln[ p_min / ( |p₀| (1 + 3α₀ε_s²/(2κ̂)) ) ],
             ε'_f := dε_{v,f}/dε_s = 3α₀ε_s / (1 + 3α₀ε_s²/(2κ̂))      (= 0 for α₀ = 0).                              (S.49)
      HAR:   x := ϖ_f/p_a solves  x² − a x^{2n} − b = 0,   a := 3k(1−n)g ε_s²,   b := (p_min/p_a)²;
             n = ½:  x = [a + (a² + 4b)^{1/2}]/2  (closed form);
             general n: f(x) := x² − a x^{2n} − b has f(0) = −b < 0 and one stationary point x_s = (na)^{1/(2−2n)} with
             f(x_s) = −(1−n) a x_s^{2n} − b < 0, hence exactly one root, in [x_s, x_hi], x_hi := max((2a)^{1/(2−2n)}, 2√b)
             (f(x_hi) ≥ x_hi²/2 − b > 0): safeguarded Newton/bisection there, |f| ≤ 10⁻¹⁴ (b + a x^{2n}), ≤ 100 iterations
             (contract constants; the bracket is exact, so this never refuses);
             ε*_f = (p_min/p_a) / (k(1−n) x^n),   ε_{v,f} = 1/(k(1−n)) − ε*_f,   q_f = 3g p_a ε_s x^n,
             ε'_f = n p_min q_f / ( (1−n) p_a² x² + n p_min² )     (= −D₁₂/D₁₁ at p = −p_min; > 0).                   (S.50)

The x-equation is p = −p_min and q = 3g p_a ε_s (ϖ/p_a)^n put into ϖ² = p² + k(1−n)q²/(3g), exact for every n (symbolic;
6 random states × n ∈ {½, 0.3, 0.7} to 1.4e−28; the n = ½ root exact). In both energies ε'_f = −D₁₂/D₁₁|_{p = −p_min},
the implicit derivative of p(ε_{v,f}(ε_s), ε_s) ≡ −p_min (symbolic for BA06; 2.4e−15 against a 30-digit FD for HAR at
the TIMs constants).

**Tangent (S.51), also written (S.32f).** When Π_f is active its principal block and spin are

      Φ_ab := ∂ε_{f,a}/∂ε_b = δ_ab − 1/3 + (1/3) ε'_f √(2/3) n̂^e_b,     spin g_ab = (ε_{f,a} − ε_{f,b})/(ε_a − ε_b) = 1 exactly,
      repeated-eigenvalue limit Φ_aa − Φ_ab = 1 + (1/3) ε'_f √(2/3) (n̂_a − n̂_b) = 1 at exactly repeated eigenvalues;   (S.51a)

Φ := I when inactive. Write Φ^tr := Φ at the trial (step 1f) and Φ^post := Φ at the converged state (step 5f). The m = 1
tangent of a floored step replaces (S.32) by

      elastic:  ã^{ep}_f = a^e(ε̃_f) Φ^tr;
      plastic:  ã^{ep}_f = a^e(ε^e_f) Φ^post [ b_{·b} Φ^tr − u Π_v v_{n+1} 1ᵀ ]_{rows ≤ 3},
                i.e. (S.32) with b_{·b} → b Φ^tr and the v-column Σ_k b_ik s_k = u_i Π_v v_{n+1} (S.45) **unchanged**: it
                multiplies δ_b of the RAW trial strain, because v = v_n exp(tr Δε) is built from the total strain and the
                floor does not touch it;                                                                              (S.32f)

and the 4th-order C is (S.33) with ã → ã^{ep}_f and the spin (σ_a − σ_b)/(ε̃_a − ε̃_b) on the raw trial ε̃ (unit spin of
Π_f, composition rule of §9.6 "m = 1"). Properties. **δ : C_f = 0 whenever the last operator applied is an active Π_f**
(p is pinned, no strain increment changes it; 2.5e−31 HAR, 1e−31 BA06 symbolic, ≤ 1.1e−15 in the O2 steps): a floored
point has no volumetric stress response and the shear stiffness of the energy at |p| = p_min — 2μ₀ under BA06 α₀ = 0;
under HAR of the order 2G(p_min), G(p_min) = g p_a (p_min/p_a)^n = 5769 kPa at the TIMs set — nonzero, which is what the
floor buys against G → 0 as p → 0. This is the exact linearisation of the floored map and is one-sided like the
elastic/plastic switch: a compressive next increment leaves the floor and gets the full stiffness. No bulk stiffness is
faked; if a global solver needs one on a floored patch, that is a solver-side regularisation to be declared, not a
kernel change. FD record (floor_fd.py: the O2 kernel with Π_f wrapped around `return_map`, BA06 K2 set with α₀ = 0 and
α₀ = 5 as the D₁₂ ≠ 0 stand-in, p_min scaled to 50 kPa, states off the WW corners, h = 10⁻⁶/10⁻⁷): (A) elastic + trial
floor 8.1e−13 / 1.6e−11; (B) plastic + trial floor 5.2e−9 / 1.7e−8 (α₀ = 0), 5.2e−9 / 5.1e−9 (α₀ = 5); (C) plastic + post
floor 4.6e−8 / 5.0e−10 (α₀ = 0), 4.2e−8 / 4.4e−10 (α₀ = 5), O(h²). Mutants: unprojected a^e 0.67–0.69; ε'_f dropped
(α₀ = 5) 3.5e−3 with δ:C = 1.1e−2 (HAR TIMs, symbolic: 1.45e−2); v-column tied to the floored trial (ã^{ep} Φ^tr) 2.6e−3 /
2.5e−3; Φ^tr dropped 0.93; Φ^post dropped 2.0–2.1. The HAR elastic floor operator itself: (S.51a) against a 30-digit FD
of σ(Π_f(ε)), 5e−19 (TIMs set, non-degenerate state); BA06 α₀ = 5: 1e−22.

**Chained substep tangent (S.54), also written (S.46f).** Each floored sub-increment contributes its two operators to the
recursion (S.46): with T_k := S^ε_k + α_k E_J (raw trial sensitivity) and T^f_k := Φ^tr_k : T_k,

      plastic:  S^ε_{k+1} = Φ^post_k : [ Φ^ε_ε̃ : T^f_k − Σ_a m^a ( (u_a/c) S^π_k + u_a Π_v S^v_{k+1} ) ],
                S^π_{k+1} = Σ_b w_b T̂^f_bb + ((1−κ)/c) S^π_k + (1−κ) Π_v S^v_{k+1};
      elastic:  S^ε_{k+1} = T^f_k,   S^π_{k+1} = S^π_k;      S^v_{k+1} = v_{k+1} (Σ_{j≤k} α_j) tr E_J unchanged (raw trace),   (S.54)

Φ^tr_k, Φ^post_k in the (S.33) form with the block (S.51a) and unit spin (identity when inactive), all in the
sub-increment's trial basis (Π_f is co-axial, so the basis of §9.6 is unchanged); the assembly (S.47) uses a^e at the
final, floored, ε^e_m. Floor events in different sub-increments compose like any other branch decision (one-sided
across an activation). FD (floor_fd.py (D), m = 2, pattern FP-/FP-, fractions held fixed): 2.9e−9 / 1.8e−9 (α₀ = 0),
2.9e−9 / 5.5e−10 (α₀ = 5); the mutant that keeps the floored states but omits the Φ operators: 0.98–0.99.

**What is given up — the stated regularisation.**
1. Energy. A projection at fixed total strain does no external work and raises the stored energy by

      E_f := Ψ(ε^e_f) − Ψ(ε^e) = ∫_{ε_{v,f}}^{ε_v} |p(ε'_v, ε_s)| dε'_v ∈ [0, p_min Δε^f_v]   (pre-floor state in dom Ψ),   (S.52)

   because |p| < p_min on the projected interval. BA06 α₀ = 0 exactly: E_f = κ̂ (p_min − |p|), Δε^f_v = κ̂ ln(p_min/|p|),
   and the bound is 1 − t ≤ −ln t, t = |p|/p_min (symbolic). Measured E_f/(p_min Δε^f_v) ∈ [0.21, 1) (BA06, α₀ = 0 and 5)
   and [0.37, 1) (HAR TIMs) over a grid of 4 ε_s × 5 |p|/p_min, never above the bound (margin 5e−4 at t = 0.999). For an
   out-of-domain trial only the bound W_f = p_min Δε^f_v is stated (the mirror-branch Ψ is unphysical). The floor is
   therefore an **energy source** of at most W_f := p_min ε^f_v per Gauss point, ε^f_v := Σ Δε^f_v over the history (trial,
   post and init events). The NorSand dissipation D^p = Δλ σ:q_a (S.38) is formed at the pre-floor converged state and
   keeps D^p ≥ 0 under (S.39) unchanged; the second-law bookkeeping of a floored point reads σ:dε = dΨ + D^p − dE_f with
   dE_f ≥ 0: **D ≥ 0 holds for the plastic flow, and the floor's violation is exactly E_f, bounded by W_f, counted and
   reported**, vanishing with p_min — which is what the F vs F/2 report tests.
2. Yield consistency. After a post floor F(σ_f, π_i) ≠ 0 in general: F_p = (η − M)/(1−N) with dp = −(p_min − |p_c|) < 0
   gives F(σ_f) < 0 (inside) on the dry side η > M, the usual approach to the tension apex, and F(σ_f) > 0 by at most
   |F_p| (p_min − |p_c|) + ζ |q_f − q_c| (O(p_min); q_f = q_c under BA06 α₀ = 0) on the wet side η < M. No F ≤ F_tol
   invariant is claimed at a floored committed state; the next increment's trial test (§9.1 step 2) sees F^tr > F_tol and
   returns. π_i is never projected: with |p| ≥ p_min the hardening target |π_i*| = |p| B^{(N−1)/N} stays bounded away
   from 0 on any bounded ψ_i, which is what keeps π_i from collapsing with p.
3. The volumetric stiffness at the floor is zero (above), and a floored committed state can sit slightly outside F
   (item 2). Nothing else changes: ε^p, π_i, v, D^p, the vertex rule (Π_f at ε_s = 0 is the axis closed form), the cap.

**Determinism, counting, interfaces.**
- Activation: p(ε^e) > −p_min (1 − 10⁻¹²) or ε^e ∉ dom Ψ. After a projection p = −p_min to round-off and the test is
  false (idempotent, no re-projection chatter). Closed forms (S.49), (S.50) for n = ½; the general-n scalar solve has the
  exact bracket above. The floor introduces **no refusal**: with p_min > 0 the trial-side `p_or_pi_nonneg` (its p part)
  and the HAR domain failure at the trial become floor events; `local_*`, `pi_*`, `B_nonpos` and π_i ≥ 0 are unchanged.
  p_min = 0 switches the floor off and restores the pre-round-3 refusals (`-pmin 0`); p_min < 0 is refused.
- Counted. Per (sub-)increment: floor_tr, floor_post ∈ {0, 1}. Per Gauss point, committed: n_f,tr, n_f,post, n_f,init ∈
  {0, 1}, ε^f_v = Σ Δε^f_v ≥ 0, W_f = p_min ε^f_v, and at_floor := (p_committed > −p_min (1 + 10⁻¹⁰)). A substepped
  increment sums its sub-increments; a refused increment counts nothing (state frozen). Shell response `floor` =
  (at_floor, n_f,tr, n_f,post, ε^f_v, W_f); `stepInfo` gains floor_tr/floor_post of the last step. Deck report (plan §2.7,
  K5): the number of Gauss points with at_floor at the limit state, Σ W_f over the mesh, and the limit load at p_min and
  p_min/2 (accepted if it moves by < 2 %).
- `initialState`: ε^e := Π_f(invert(σ₀)) (counted, n_f,init = 1); the π_{i0} rule (S.53) then uses the floored p_init. At a
  floored point the deck's σ₀ is replaced by σ(ε^e_f) (p = −p_min; under HAR q also scaled by (ϖ_f/ϖ)^n) and the first
  equilibrium iteration absorbs the O(p_min) imbalance. Refusing instead would kill every deck whose first row of Gauss
  points sits above p_min (geostatic p ≈ 0.64 σ_v at K₀ = 0.4554).
- LogStrain provider (plan §2.8, option c): `getElasticStrain` returns the committed ε^e_f, which reproduces σ exactly
  and lies in dom Ψ by construction, so the wrapper never receives an out-of-domain strain; it rebuilds b^e_n from it,
  so Δε^f joins the inelastic part of the multiplicative split; v = v₀J is untouched (v from the total deformation);
  (S.34) takes ã^{ep}_f with the raw trial stretches.
- Default: **p_min = 5·10⁻³ p_ref**, p_ref = |p₀| (BA06) or p_a (HAR, §2.4): 0.5 kPa on the K2 set, 0.505 kPa on the TIMs
  set — the plan's "≤ 0.5 kPa", in the deck's stress units through p_ref (no unit-blind constant). The SANISAND precedent
  10⁻³ p_atm (0.1 kPa) is a *stiffness* floor under the moduli, not a stress projection, and is not transferred. The ring
  states (0.19–10 kPa) and the Kimura deck (0.57 kPa at s/B 0.002) make 0.5 kPa engage at the free-surface ring only.
- O1 (rate form, S.55). On the floor the continuum statement is a second mechanism ε̇^f = λ̇_f P, P := 1/3, with ṗ = 0
  added to Ḟ = 0 (Koiter):

      [ f:a^e:q + H    f:a^e:P ] [ λ̇   ]   [ f:a^e:ε̇ ]
      [ P:a^e:q        P:a^e:P ] [ λ̇_f ] = [ P:a^e:ε̇ ],      λ̇, λ̇_f ≥ 0 (a negative multiplier drops its mechanism),
      ε̇^e = ε̇ − λ̇ q − λ̇_f P;    floor only: λ̇_f = (P:a^e:ε̇)/(P:a^e:P)    (BA06 α₀ = 0: P:a^e:P = K, λ̇_f = tr ε̇ exactly).   (S.55)

  The floor is active while p = −p_min and the unconstrained ṗ would be > 0, released when λ̇_f would turn negative.
  The O2/kernel split (1f → return → 5f) is the projection splitting of this constrained flow and converges to it at
  **first order** (floor_rate.py: BA06 K2 set, α₀ = 0 and 5, p_min = 50 kPa, a loading path along the stress deviator
  with slight expansion that keeps both mechanisms active — min λ̇ = 0.38, min λ̇_f = 0.20 over every RHS evaluation —
  Radau at rtol 10⁻¹⁰; the split with m = 1 … 64 sub-increments misses the rate solution by 7.2e−3 … 1.4e−4 in σ and
  2.2e−2 … 4.1e−4 in π_i with observed orders 0.87, 0.92, 0.96, 0.98, 0.99, 0.99; the floor-only branch under BA06
  α₀ = 0 gives λ̇_f = tr ε̇ to round-off with P:a^e:P = K = 5000 kPa). O2 → O1 first-order convergence on a floored path is a gate test; O1 may also run the
  split at its output increments, but then it is not an independent oracle for the floor.

**Tests and mutants (K1.12–K1.15, §13; named for the mutation gate, plan §5.3).** M-F1 floor not applied (the
pre-round-3 refusal, or the state passed through with |p| < p_min): K1.13 expects a committed p = −p_min and no refusal.
M-F2 applied but not counted: K1.12/13 assert n_f, ε^f_v, W_f. M-F3a unprojected tangent (plain a^e or ã^{ep}): δ:C ≠ 0,
0.67–2.1 off the FD. M-F3b ε'_f dropped: δ:C ≠ 0 under HAR or α₀ ≠ 0 only (invisible under BA06 α₀ = 0, so the test runs
under HAR). M-F3c v-column tied to the floored trial: 2.5e−3 off the FD of a plastic trial-floored step. M-F3d Φ
operators omitted from the chain: 0.98 off the FD of a substepped floored increment. M-F4 the HAR floor computed with
the BA06 inverse, or HAR replaced by BA06 at the floor: ε_{v,f} differs by orders of magnitude (0.0530 vs 9.84e−4;
K1.12 vs K1.13). M-F5 trial floor skipped (post only): a HAR out-of-domain trial refuses instead of flooring. M-F6
projection direction wrong (q kept): under HAR ε_s^e must be unchanged and q_f = 3g p_a ε_s x^n (K1.14). M-F7 π_i or v
altered by the projection: K1.12 asserts both unchanged. M-F8 default reference wrong (|p₀| vs p_a) or units: parser
echo. M-F9 floor applied inside the local Newton: the K1.12 (C) values differ (the post-floored q is not the in-Newton
one).
"""
rep("checks are quoted in the paragraph above.\n\n---\n\n## 10. Q-cap",
    "checks are quoted in the paragraph above.\n" + S97 + "\n---\n\n## 10. Q-cap")

# ---------------------------------------------------------------------------------------------- §10.2 scan step
rep("the root-selection contract of §8 is what makes the kernel's answer the one continuous with π_{i,n}: scan\n"
    "step ≤ 10⁻⁴|π_{i,n}| (≪ the dip width), no factor-2 bracket, Δλ bounded with backtracking when the nested solve fails or\n"
    "jumps, substep on refusal (§9.1).",
    "the root-selection contract of §8 is what makes the kernel's answer the one continuous with π_{i,n}: scan\n"
    "step 10⁻³|π_{i,n}| (the ramp-width criterion (S.56) below; the G1 value 10⁻⁴ is withdrawn), no factor-2 bracket, Δλ\n"
    "bounded with backtracking when the nested solve fails or jumps, substep on refusal (§9.1).")
S56 = r"""

**Scan step: 10⁻³, not 10⁻⁴ — the fold geometry (round 3, 2026-10-03; floor_fold.py) [I; measured].** The G1 text set the
scan step by "≪ the dip width" (0.47 kPa at the step-8 iterate, 0.052 kPa at Δλ = 8.0e−4) and wrote ≤ 10⁻⁴|π_{i,n}|,
while O2 and the kernel ship PI_SCAN_REL = 10⁻³ (plan §2.8). The dip criterion is the wrong one: the dip width vanishes
at the fold (the near and middle roots annihilate), so no fixed step is "≪ the dip", and none is needed. What the scan
must guarantee is never to select the far root in silence, which needs a scan point between the middle root π₂ and the
first local extremum of |r| beyond it (at distance d_x from π_{i,n}): a step that skips both near roots still lands
before that extremum, sees |r| grow, and raises `pi_fold` — the designed backtrack — as long as PI_SCAN_REL ≪ d_x. d_x
scales with the width of the cap ramp in π_i at fixed p, a pure parameter constant because η(p, π_i) depends on p/π_i
only:

      W_ramp := 1 − π_i(η₂)/π_i(η₁) = 1 − [ (1 − c₂N)/(1 − c₁N) ]^{(1−N)/N}   (N > 0),    = 1 − e^{−(c₂ − c₁)}   (N = 0);
      contract:  PI_SCAN_REL ≤ W_ramp/10,  enforced in validate() (refuse a narrower smooth cap). Defaults: W_ramp = 0.0606;
      c₁ = 0.05 with c₂ = 0.07 gives 0.0122 (admissible), with c₂ = 0.06 gives 0.0061 (refused).                      (S.56)

Measured on the AMP_STOP path (K2 paper set, c₁ = 0.05, c₂ = 0.15, shipped O2, exponential v-update): the ramp width is
6.060e−2 |π_i(η₂)| at p = −100, −160, −172 kPa (3.09, 4.94, 5.32 kPa; the estimate (c₂ − c₁)|π_i|(π_i/p)^{N/(1−N)}
overstates it by 5 %). At the converged iterate of every non-substepped plastic step the distance to the middle root is
d₂ ≥ 2.19e−3, 3.94e−3, 9.69e−3 |π_{i,n}| (n = 40, 80, 160: 2–10× the step) and d_x ≥ 5.49e−2, 5.05e−2, 4.22e−2 |π_{i,n}|
(42–55× the step, ≈ 0.7–0.9 W_ramp). At the step-8 iterate (n = 40, π_{i,n} = −80.000 kPa) the roots sit at 0.0172 /
0.4828 / 24.68 kPa from π_{i,n} for Δλ = 7.877e−4 (the G1 numbers), at 0.1348 / 0.1868 / 24.92 kPa for 8.0e−4 (dip
0.052 kPa = 6.5e−4|π_{i,n}|, narrower than the 10⁻³ step — and the near root is still the one bracketed, 7 evaluations),
and the near pair is gone at 8.05e−4: both scan steps (10⁻³ with PI_SCAN_MAX = 1000 and 10⁻⁴ with 10 000, same range)
select the same root at every Δλ ≤ 8.0e−4 and both raise `pi_fold` at 8.05e−4, 8.1e−4 and 8.3e−4. The path runs to
completion at n = 20, 40, 80, 160 with 2341, 3069, 3072, 349 `pi_fold` backtracks inside the line searches (none a
refusal; substep levels up to 2⁸ at n = 20, 2² at n = 40, 2¹ at n = 80, none at n = 160). A 10⁻⁴ step would resolve dips
down to 10⁻⁴ at 10× the nested cost and a 10× shorter range unless PI_SCAN_MAX rose to 10⁴ (a hardening step moving
π_i by more than 10 % would otherwise refuse with `pi_nobracket`). **Ruling: 10⁻³ is the correct contract value; the
sheet is corrected (§8, §16.3, §16.6 item 4); O2 and the kernel need no change**; the P2 census's 0 violations in
12 476 replays is consistent with the geometry. New requirement for O2/kernel `validate()`: the (S.56) refusal.
"""
rep("With the shipped near-root selection the O2 oracle completes the AMP_STOP path at every n from 20 to 320.\n",
    "With the shipped near-root selection the O2 oracle completes the AMP_STOP path at every n from 20 to 320." + S56)

# ---------------------------------------------------------------------------------------------- §13 items 12-15
S13 = r"""12. **Floor projection, BA06 (K1.12; §9.7, (S.48)–(S.49), (S.52)) [I; floor_sympy.py (g), floor_fd.py]:** K2 set,
   p_min = 5·10⁻³|p₀| = 0.5 kPa, α₀ = 0: ε_{v,f} = −κ̂ ln(p_min/|p₀|) = 0.0529831737; a trial at p^tr = −p_min/2 has
   ε_{v,tr} = 0.0599146455, so Δε^f_v = κ̂ ln 2 = 6.931471806e−3, W_f = 3.465735903e−3 kPa, E_f = κ̂ (p_min − |p^tr|) =
   2.5e−3 kPa (≤ W_f); committed p = −0.5 exactly, with q, n̂, π_i, v unchanged; tangent block ã_f = 2μ₀ (δ_ab − 1/3) =
   10 800 (δ_ab − 1/3) kPa, δ:C = 0. The same with the floor scaled to 50 kPa is the FD record of §9.7: (A) p^tr −33.3 →
   −50.00 elastic, Δε^f_v = 4.066e−3; (B) plastic from a floored trial, p_c = −50.70; (C) wet side, p_c = −48.76 → −50.00
   by the post floor, η = 0.953.
13. **HAR floor and the out-of-domain trial (K1.13; (S.50)) [I; floor_sympy.py (g)]:** TIMs set (§2.3), p_min = 5·10⁻³ p_a =
   0.505 kPa: ε_{v,f}(ε_s = 0) = 9.83645033e−4 (domain edge 1/(k(1−n)) = 1.05849170e−3); G(p_min) = g p_a (p_min/p_a)^{1/2} =
   5769.156 kPa, floored shear stiffness 2G = 11 538.31 kPa, volumetric response 0. An isotropic state at p = −1 kPa
   (ε_v = 9.53167838e−4) with Δε_v = +1e−4 is still in the domain (p^tr = −2.56e−3 kPa) and floors to ε_{v,f} with
   Δε^f_v = 6.9522805e−5, W_f = 3.5109017e−5 kPa; with Δε_v = +1.1e−4 the trial is **out of the domain** (ε_v = 1.0632e−3 >
   1.0585e−3) and floors to the same ε_{v,f} with Δε^f_v = 7.9522805e−5 — no refusal (M-F1, M-F5). Under BA06 the same two
   increments are ordinary elastic steps (M-F4).
14. **HAR floored state under shear (K1.14; (S.50)–(S.51)) [I; floor_sympy.py (d), (g)]:** TIMs set, ε_s^e = 2e−4 at the
   floor: x = ϖ_f/p_a = 0.09185198375, ε_{v,f} = 1.04102893e−3, q_f = 14.8362051 kPa (η_f = 29.38; the ring's high-η states
   are of this kind), ε'_f = 0.08679793; ε_s^e is unchanged by Π_f (M-F6); δ:C_f = 0 with the ε'_f term and δ:C/max|C| =
   1.45e−2 without it (M-F3b).
15. **π_{i0} rule (K1.15; (S.53)) [I; floor_sympy.py (f)]:** K2 set, p_init = −100 kPa: η* = c₂M = 0.18 → π_{i0} = −50.995881;
   η* = 0 → −46.475800 (apex; the §9.6 (F) value); η_init = 0.75 → −71.554175; η* = M → −100. First yield on a drained TXC
   path from the isotropic state with π_{i0} = −50.9959: q = 11.3469, p = −103.7823 kPa, η = 0.1093, w = 0.337; on a
   constant-p path η = 0.18 exactly.
"""
rep("which is the discriminating check between the two energies.\n",
    "which is the discriminating check between the two energies.\n" + S13)

# ---------------------------------------------------------------------------------------------- §16.1 item 5
rep("   gated 2026-10-02 (§2.3, §16.6 item 15); BA06's energy with α₀ = 0 stays the default and the paper-mode energy.\n",
    "   gated 2026-10-02 (§2.3, §16.6 item 15); BA06's energy with α₀ = 0 stays the default and the paper-mode energy.\n"
    "   p′ floor: designed 2026-10-03 (§9.7, owner decision (b)); π_{i0} rule: §5.4 (owner decision (c)); the nested-scan\n"
    "   step is 10⁻³ (§10.2 (S.56)).\n")

# ---------------------------------------------------------------------------------------------- §16.3
rep("- The tension apex p → 0 (q → 0 with η → M/N): the vertex rule of §3.2 applies to the direction; the p-floor rule of\n"
    "  the plan (§2.7) handles the magnitude and is outside this sheet. R_tol value to be confirmed by the census.\n",
    "- The tension apex p → 0 (q → 0 with η → M/N): the vertex rule of §3.2 applies to the direction; the magnitude is the\n"
    "  p′ floor of §9.7 (round 3; no longer outside this sheet). R_tol value to be confirmed by the census.\n")
rep("root-selection rule of §8 (root continuous with π_{i,n}, scan step ≤ 10⁻⁴|π_{i,n}|, no factor-2 bracket), bounded",
    "root-selection rule of §8 (root continuous with π_{i,n}, scan step 10⁻³|π_{i,n}| — corrected at round 3 from the\n"
    "  10⁻⁴ first written here, §10.2 (S.56) — no factor-2 bracket), bounded")
rep("  the parser refusals (k, g > 0, 0 ≤ n < 1, p_a > 0) and the ε* ≤ 0 domain refusal (p-floor); confirm that the §3.2\n"
    "  vertex rule and the §10 cap need nothing new (they use a^e through (S.3) only). The ring dumps' README says\n"
    "  \"compression negative\" but the CSVs are compression-positive (§2.3) — a note for the TIMs side, not for this sheet.\n",
    "  the parser refusals of §2.4 and the ε* ≤ 0 floor event of §9.7; confirm that the §3.2\n"
    "  vertex rule and the §10 cap need nothing new (they use a^e through (S.3) only). The ring dumps' README says\n"
    "  \"compression negative\" but the CSVs are compression-positive (§2.3) — a note for the TIMs side, not for this sheet.\n"
    "- **[round 3, 2026-10-03] p′ floor gate tests** (test author; expected values from §9.7 and §13.12–13.14): the FD\n"
    "  checks of (S.32f) (cases (A) elastic + trial floor, (B) plastic + trial floor, (C) plastic + post floor) and of (S.54)\n"
    "  (a substepped floored increment) under BA06 α₀ = 0 and under HAR, off the WW corners, with the mutants M-F1…M-F9\n"
    "  named in §9.7; O2 → O1 first-order convergence on a floored path ((S.55) in O1); kernel-vs-O2 parity on a floored\n"
    "  path including the counters; the HAR out-of-domain trial (13.13) through the OpenSees shell; the `floor` response\n"
    "  and the F vs F/2 deck report (K5). Design decisions still open for the owner: whether the kernel reports E_f (needs\n"
    "  Ψ) or only W_f; the global-solver behaviour on a floored patch (zero bulk stiffness is the consistent tangent, §9.7).\n"
    "- **[round 3] validate() gains the (S.56) refusal** (PI_SCAN_REL ≤ W_ramp/10) in O2 and the kernel — the only code change\n"
    "  the scan-step ruling asks for; PI_SCAN_REL = 10⁻³ and PI_SCAN_MAX = 1000 stay.\n"
    "- **[round 3] π_{i0} rule (S.53)**: O2's `initial_state` default (apex through p_init) becomes the c₂-rule; the kernel's\n"
    "  `initialState` and O1 take the same rule with c₂ = 0 for cap = none; an explicit `-pi0` overrides.\n")

# ---------------------------------------------------------------------------------------------- §16.4
rep("  off-axis modulus ratios. har_k1.py: the §13 items 1h and 11h values at 20 digits.\n",
    "  off-axis modulus ratios. har_k1.py: the §13 items 1h and 11h values at 20 digits.\n"
    "\n"
    "p-floor round checks (2026-10-03, local py312 venv with the O2 oracle imported read-only; scripts and logs in\n"
    "`Ladruno_files/testbed/norsand_oracle/tests/scratch_sheet_floor/`):\n"
    "- floor_sympy.py (sympy/mpmath, 30 digits): (a) the (S.50) x-equation ⇔ p = −p_min for every n (symbolic), the n = ½\n"
    "  root exact, 6 random states × n ∈ {½, 0.3, 0.7} to 1.4e−28, the f(x_s) identity 1.5e−31; (b) ε'_f = −D₁₂/D₁₁ symbolic\n"
    "  (HAR stress form, BA06 with α₀) and 2.4e−15 vs a 30-digit FD (HAR TIMs); (c) (S.49) p = −p_min symbolic; (d) (S.51a)\n"
    "  vs a 30-digit FD of σ(Π_f(ε)): 5e−19 (HAR TIMs), 1e−22 (BA06 α₀ = 5); δ:C_f = 0 to 1e−31; unit spin 4e−31; the\n"
    "  ε'-dropped mutant at δ:C/max|C| = 1.45e−2 / 9.1e−5; (e) (S.52) symbolic for BA06 α₀ = 0, a 4 ε_s × 5 |p|/p_min grid\n"
    "  for HAR TIMs and BA06 α₀ = 5 and 0: 0 ≤ E_f ≤ p_min Δε^f_v with min E_f/bound 0.37 / 0.21 / 0.21; (f) (S.53)\n"
    "  identities 4e−31 / 0, N → 0 limit 2.4e−10 at N = 10⁻⁹, special values exact, monotonicity 9e−20; (g) the §13.12–13.14\n"
    "  numbers.\n"
    "- floor_fd.py (O2 kernel, BA06 K2 set, α₀ = 0 and 5, p_min = 50 kPa, off-corner states, h = 10⁻⁶/10⁻⁷): the §9.7 FD\n"
    "  record — (A) 8e−13 / 1.6e−11, (B) 5.2e−9 / 1.7e−8 and 5.2e−9 / 5.1e−9, (C) 4.6e−8 / 5.0e−10 and 4.2e−8 / 4.4e−10, (D)\n"
    "  2.9e−9 / 1.8e−9 and 2.9e−9 / 5.5e−10; mutants 0.67–2.1 (projection, operators) and 2.5e−3–3.5e−3 (v-column, ε'_f).\n"
    "  The first draft put (B) and (C) at the exact TXC corner and measured the O(h) of §4.3 (1.0e−4 → 1.0e−5), not a floor\n"
    "  error; a second draft tied the v-column to the floored trial and measured the 5e−3 that is now the M-F3c mutant.\n"
    "- floor_rate.py (O2 primitives, BA06 K2 set, α₀ = 0 and 5, p_min = 50 kPa, Radau rtol 10⁻¹⁰, a loading path on\n"
    "  which both mechanisms of (S.55) stay active, min λ̇ = 0.38 / 0.39, min λ̇_f = 0.20 / 0.20): the split (1f → return →\n"
    "  5f) vs the rate solution at m = 1, 2, 4, 8, 16, 32, 64: σ 7.2e−3 → 1.4e−4 (α₀ = 0), 7.5e−3 → 1.4e−4 (α₀ = 5), π_i\n"
    "  2.2e−2 → 4.1e−4, observed orders 0.87, 0.92, 0.96, 0.98, 0.99, 0.99 (first order); p = −p_min to 6e−9 along the\n"
    "  rate path; floor-only branch λ̇_f = tr ε̇ exactly (P:a^e:P = K = 5000 kPa). A first draft on a NorSand-unloading\n"
    "  path had λ̇ < 0 (one mechanism only) and is recorded as the reason the path is a loading one.\n"
    "- floor_fold.py (shipped O2 on the AMP_STOP smooth-cap path): the §10.2 (S.56) numbers (d₂ and d_x per n, the step-8\n"
    "  roots through the fold under both scan steps, the ramp widths and the (S.56) constant).\n")

# ---------------------------------------------------------------------------------------------- §16.6
rep("item 14 added at G2 (owner\ndecision 2026-10-01, exponential v-update).",
    "item 14 added at G2 (owner\ndecision 2026-10-01, exponential v-update); item 15 at the HAR round (2026-10-02); item 16 at the p-floor round\n(2026-10-03).")
rep("root-selection contract (root continuous with π_{i,n}, scan step ≤ 10⁻⁴|π_{i,n}|, no factor-2 bracket)",
    "root-selection contract (root continuous with π_{i,n}, scan step ≤ 10⁻⁴|π_{i,n}| [corrected to 10⁻³ at item 16], no factor-2 bracket)")
ROW16 = ("| 16 | §1.3, §2.3, §2.4 (new), §5.4 (new), §8, §9.1, §9.7 (new), §10.2, §13.12–15, §16.1, §16.3, §16.4 | "
         "**p-floor round (2026-10-03; owner decisions (a)–(e) of 2026-10-02).** (1) The p′ floor as the operator Π_f (S.48): "
         "a strain-space projection at fixed deviatoric elastic strain onto p = −p_min, applied to the trial (step 1f) and "
         "to the converged state (step 5f), never inside the local Newton; closed forms (S.49) BA06 (any α₀) and (S.50) HAR "
         "(n = ½ closed, general n bracketed); tangent (S.51)/(S.32f) with the v-column on the raw trial, chain (S.54)/(S.46f); "
         "δ:C_f = 0 at a floored point (the exact linearisation; the energy's deviatoric stiffness at |p| = p_min); an energy "
         "source E_f ∈ [0, p_min Δε^f_v] (S.52) counted as W_f, with D^p ≥ 0 of the flow untouched; F ≤ F_tol not claimed at a "
         "floored committed state; counters and the `floor` response; initialState projects; the LogStrain provider returns "
         "the floored ε^e; default p_min = 5·10⁻³ p_ref (0.5 / 0.505 kPa); O1 rate form (S.55) with the split converging at "
         "first order; mutants M-F1…M-F9; K1.12–14. (2) The unified π_{i0} rule (S.53) with its limits and the K2 values (K1.15). "
         "(3) §2.4 energy-option bookkeeping: the branch table, the parser refusals (BA06 flags with HAR refused, not ignored), "
         "the FD re-run list under HAR. (4) The nested-scan step: 10⁻⁴ (G1 text) → 10⁻³ (O2/kernel), with the ramp-width "
         "criterion (S.56) and a new validate() refusal; O2 and the kernel unchanged. | [I; round 3], owner decisions "
         "2026-10-02 | floor_sympy.py (all OK), floor_fd.py (all OK), floor_rate.py, floor_fold.py |\n")
rep("| har_sympy3.py (all ≤ 3e−30 / exact), har_gate.py, har_k1.py |\n",
    "| har_sympy3.py (all ≤ 3e−30 / exact), har_gate.py, har_k1.py |\n" + ROW16)
rep("in the S^v closed form of (S.46) (v₀ → v_{k+1}), both re-verified (§16.4, G2 checks).\n",
    "in the S^v closed form of (S.46) (v₀ → v_{k+1}), both re-verified (§16.4, G2 checks). Item 16 adds new FD-checked\n"
    "algebra ((S.48)–(S.55)) and changes no existing formula; its one correction to existing text is the scan-step value of\n"
    "§8/§10.2 (10⁻⁴ → 10⁻³), a contract constant, not a derivative.\n")

# ---------------------------------------------------------------------------------------------- apply
for old, new in EDITS:
    n = text.count(old)
    if n != 1:
        sys.exit(f"ANCHOR MATCHED {n} TIMES (need 1):\n---\n{old[:200]}\n---")
    text = text.replace(old, new)
if "RATE_RESULT" in text:
    sys.exit("RATE_RESULT placeholder still present: fill it from floor_rate.log before applying")
open(OUT, "w", encoding="utf-8", newline="\n").write(text)
print(f"applied {len(EDITS)} edits -> {OUT}")
