# -*- coding: utf-8 -*-
"""Round-3b edits (2026-10-03) of the WP-144a equation sheet: the five Adversary defects of round 3, the owner's
approval of option (i) and the orchestrator recommendation beside the two open §16.3 questions.

Applies to the CURRENT working-copy sheet (the round-3 edits already applied, uncommitted). Each edit is an
exact-anchor replacement that must match EXACTLY ONCE; a drifted or duplicated anchor aborts before anything is
written. Line endings and encoding are preserved byte for byte (the sheet is LF, UTF-8).
    python apply_sheet_round3b.py <sheet.md> [<out.md>]        (no <out.md>: edit in place)
Written this way because the Deriver session cannot write into the WP-144 worktree.
"""
import sys

SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else SRC
with open(SRC, encoding="utf-8", newline="") as fh:
    text = fh.read()
EDITS = []


def rep(old, new):
    EDITS.append((old, new))


# ---------------------------------------------------------------------------------------------- frontmatter
rep('status: "p-floor round 2026-10-03 (',
    'status: "round 3b 2026-10-03 (the Adversary\'s five round-3 defects fixed on the round-3 sheet: §9.7 items 2–3 — '
    'the post floor is a wet-side event under BA06 only; under HAR the D₁₂ coupling of the plastic shear strain makes '
    'the dry-side pattern FPf, reproduced on the TIMs set with (S.32f)/(S.54) FD-checked on the HAR law (K1.14b, '
    '2.2e−7 / 2.2e−9 at h = 10⁻⁷/10⁻⁸, O(h²)); the wet-side F(σ_f) = O(p_min) stated with its size (+0.27 … +0.82 '
    'p_min) and the self-limiting contraction (|Δε^p_v| ≤ 2.7e−5 < ε*_f = 7.5e−5); §10.2 (S.56) the W_ramp words '
    'corrected to 1 − π_i(η₁)/π_i(η₂) and the validate() refusal gated to cap = smooth; §2.4/§15 p_a is one flag for '
    'HAR and the fork CSL, TIMs value 101 kPa, and the HAR→BA06 mutant is killed by expected-value tests and O2 parity, '
    'not by self-FD; owner approval of option (i) Π_f recorded (§9.7, §16.1, §16.3) with the orchestrator '
    'recommendation beside the two still-open §16.3 questions; scripts in '
    'gitAPE/ladrunoNORDSAD/handoff/round3b/scratch/). p-floor round 2026-10-03 (')

# ---------------------------------------------------------------------------------------------- §1.3 rows
rep("reference pressure (> 0; p = −p_a at ε^e = 0; also the p₀ of the §9.1 scalings). Replaces the five BA06 entries "
    "above (α₀ included) | HAR05 40–41 |",
    "reference pressure (> 0; p = −p_a at ε^e = 0; also the p₀ of the §9.1 scalings; **the same single parameter as "
    "the fork-CSL p_a below — one `-p_a` flag, one value, §2.4**). Replaces the five BA06 entries above (α₀ "
    "included) | HAR05 40–41 |")
rep('| e₀, λ_c, ξ, p_a | "fork" CSL: e_c = e₀ − λ_c (−p/p_a)^ξ (DM04 form) | plan §2.3 |',
    '| e₀, λ_c, ξ, p_a | "fork" CSL: e_c = e₀ − λ_c (−p/p_a)^ξ (DM04 form); p_a shared with the HAR energy (one flag, '
    '§2.4; TIMs: 101 kPa) | plan §2.3 |')

# ---------------------------------------------------------------------------------------------- §2.4 p_a, mutant
rep("| parameters | p₀ < 0, κ̂ > 0, ε_{v0}, μ₀ > 0, α₀ (0 in both papers) | k > 0, g > 0, 0 ≤ n < 1, p_a > 0 — and "
    "**none** of the five BA06 values |",
    "| parameters | p₀ < 0, κ̂ > 0, ε_{v0}, μ₀ > 0, α₀ (0 in both papers) | k > 0, g > 0, 0 ≤ n < 1, p_a > 0 (p_a is "
    "the one `-p_a` flag, shared with the fork CSL (S.22)) — and **none** of the five BA06 values |")
rep("(refused, never ignored: HAR replaces α₀ and the other four are not read). BA06: any of `-k -g -n` given. Both:\n",
    "(refused, never ignored: HAR replaces α₀ and the other four are not read). BA06: any of `-k -g -n` given — `-p_a`\n"
    "is **not** in this list: it is the fork-CSL parameter as well, read under `-csl fork` in either energy, and under\n"
    "BA06 with the paper CSL it is unused and only echoed [round 3b]. **One flag, one value:** with `-energy HAR -csl\n"
    "fork` the energy's p_a and the CSL's p_a are the same number, the p_atm of the DM04 calibration the set came from —\n"
    "TIMs: **101 kPa** (the campaign set's `Patm 101`, ring-dump README); 101.325 kPa is only O2's inactive default and\n"
    "the stale example of the pre-round-3b §15 row. The p_min default 5·10⁻³ p_a = 0.505 kPa (§9.7) and K1.13/K1.14\n"
    "are at 101 kPa and move with p_a. Both:\n")
rep("(S.51)/(S.54) (floor_fd-type, under HAR: the round-3 scratch ran them under BA06 with α₀ = 5 as the D₁₂ ≠ 0 stand-in\n"
    "and the elastic floor operator under HAR symbolically, §16.4); kernel-vs-O2 parity on a HAR path; K1.2 (closed loop)\n"
    "under HAR. The mutant \"HAR silently replaced by BA06\" is killed by 13.11 (η = 3gε_s exactly on a constant-volume\n"
    "elastic shear, against p constant under BA06), by 13.1h and 13.13 (the domain edge and the floor closed form differ\n"
    "from the BA06 values by orders of magnitude), and by any FD-tangent test with q > 0 (the t4 term of (S.3) is zero\n"
    "under BA06 α₀ = 0 and not under HAR).\n",
    "(S.51)/(S.54) (floor_fd-type, under HAR: the round-3 scratch ran them under BA06 with α₀ = 5 as the D₁₂ ≠ 0 stand-in\n"
    "and the elastic floor operator under HAR symbolically; round 3b ran the plastic dry-side FPf case and the chains\n"
    "FPf,FPf / -P-,-Pf under HAR with the O2 return map on the HAR law, K1.14b, §16.4); kernel-vs-O2 parity on a HAR\n"
    "path; K1.2 (closed loop) under HAR. The mutant \"HAR silently replaced by BA06\" is killed by 13.1h (p on the\n"
    "isotropic axis), 13.11 (η = 3gε_s exactly on a constant-volume elastic shear, against p constant under BA06), 13.13\n"
    "and 13.14 (the domain edge and the floor closed forms differ from the BA06 values by orders of magnitude), and by\n"
    "kernel-vs-O2 parity on a HAR path with O2 running the HAR law — **not by any FD-tangent test**: a consistent swap\n"
    "(BA06 stress with the BA06 tangent) passes its own FD check, because an FD test compares a kernel with itself; the\n"
    "t4 term of (S.3) being live under HAR is seen only against an expected value or an independent oracle [round 3b].\n")

# ---------------------------------------------------------------------------------------------- §9.7
rep("**Adopted: (i), as the operator Π_f below**, applied at two points of every (sub-)increment:",
    "**Adopted: (i), as the operator Π_f below** (owner-approved 2026-10-03; §16.1 item 5, §16.3), applied at two "
    "points of every (sub-)increment:")
rep("2.5e−3; Φ^tr dropped 0.93; Φ^post dropped 2.0–2.1. The HAR elastic floor operator itself: (S.51a) against a 30-digit FD\n"
    "of σ(Π_f(ε)), 5e−19 (TIMs set, non-degenerate state); BA06 α₀ = 5: 1e−22.\n",
    "2.5e−3; Φ^tr dropped 0.93; Φ^post dropped 2.0–2.1. The HAR elastic floor operator itself: (S.51a) against a 30-digit FD\n"
    "of σ(Π_f(ε)), 5e−19 (TIMs set, non-degenerate state); BA06 α₀ = 5: 1e−22. **Under HAR with the plastic part live**\n"
    "(round 3b, floor_fd_har.py: the O2 return map on the HAR law, TIMs set, p_min = 0.505 kPa, off-corner): the dry-side\n"
    "case with both floors active, FPf (K1.14b), 2.2e−5 / 2.2e−7 / 2.2e−9 at h = 10⁻⁶/10⁻⁷/10⁻⁸ (ratios 100 and 101: O(h²);\n"
    "the increment is 2e−5, so h = 10⁻⁶ is 5 % of it); mutants Φ^post dropped 0.14, Φ^tr dropped 0.12, both dropped 2.2.\n")
rep("2.9e−9 / 5.5e−10 (α₀ = 5); the mutant that keeps the floored states but omits the Φ operators: 0.98–0.99.\n",
    "2.9e−9 / 5.5e−10 (α₀ = 5); the mutant that keeps the floored states but omits the Φ operators: 0.98–0.99. Under HAR\n"
    "(round 3b, floor_fd_har.py (D1)–(D2), m = 2, fractions held fixed, K1.14b): pattern FPf,FPf 2.3e−7 / 2.3e−9 and\n"
    "pattern -P-,-Pf 1.2e−6 / 1.2e−8 (h = 10⁻⁷/10⁻⁸, ratio 100 both); Φ operators omitted 2.1 for both.\n")
rep("2. Yield consistency. After a post floor F(σ_f, π_i) ≠ 0 in general: F_p = (η − M)/(1−N) with dp = −(p_min − |p_c|) < 0\n"
    "   gives F(σ_f) < 0 (inside) on the dry side η > M, the usual approach to the tension apex, and F(σ_f) > 0 by at most\n"
    "   |F_p| (p_min − |p_c|) + ζ |q_f − q_c| (O(p_min); q_f = q_c under BA06 α₀ = 0) on the wet side η < M. No F ≤ F_tol\n"
    "   invariant is claimed at a floored committed state; the next increment's trial test (§9.1 step 2) sees F^tr > F_tol and\n"
    "   returns.",
    "2. Yield consistency. After a post floor F(σ_f, π_i) ≠ 0 in general. **Which side of M triggers a post floor depends\n"
    "   on the energy** [round 3b; floor_fd_har.py]: the return from a (floored) trial changes p by D₁₁Δε^e_v + D₁₂Δε^e_s with\n"
    "   Δε^e = −Δλ q_a. Under BA06 α₀ = 0 (D₁₂ = 0, q independent of ε_v) a dry-side return (η > M, Δε^p_v > 0, dilative)\n"
    "   compresses p and a wet-side one (Δε^p_v < 0) relaxes it, so there the post floor is a **wet-side** event (K1.12 (C))\n"
    "   and FP- is the only dry-side pattern (K1.12 (B)). Under HAR (D₁₂ < 0) the plastic shear strain adds\n"
    "   D₁₂·(−Δε^p_s) > 0, which can beat the volumetric term: TIMs set, p_min = 0.505 kPa, surface state p = −0.6 kPa,\n"
    "   η = 1.2M, θ = 0.271 (off-corner), ψ_i = −0.10, expansion 2e−5 + shear 2e−5 along n̂ — trial p = −0.461, floored to\n"
    "   −0.505, return D₁₁Δε^e_v = −0.084 kPa against D₁₂Δε^e_s = +0.100 kPa, p_c = −0.488 > −p_min, post floor active:\n"
    "   the pattern **FPf on the dry side** (K1.14b). At that state the floored trial returns above the floor (FPf) for shear\n"
    "   amplitudes up to 8e−5 at any expansion ≥ 2e−5 and below it (FP-) from 1.6e−4 on, where the larger Δλ lets the\n"
    "   dilative volumetric term win; a pure expansion floors elastically (FE-) — the single-step map of fpf_region_har.py.\n"
    "   The sign of F(σ_f) after a post floor follows dF ≈ F_p dp + ζ dq, F_p = (η − M)/(1−N), dp = −(p_min − |p_c|) < 0,\n"
    "   dq = q_f − q_c ≥ 0 (= 0 under BA06 α₀ = 0; > 0 under HAR, where q grows with the floored |p| through ϖ): dry side\n"
    "   (F_p > 0) F(σ_f) < 0, inside, unless the ζ dq term wins (the FPf case: F(σ_f)/p_min = −0.004); wet side (F_p < 0)\n"
    "   F(σ_f) > 0 by O(p_min) — the first-order size |F_p| (p_min − |p_c|) + ζ (q_f − q_c), measured in item 3. No F ≤ F_tol\n"
    "   invariant is claimed at a floored committed state; the next increment's trial test (§9.1 step 2) sees F^tr > F_tol and\n"
    "   returns.")
rep("3. The volumetric stiffness at the floor is zero (above), and a floored committed state can sit slightly outside F\n"
    "   (item 2). Nothing else changes:",
    "3. The volumetric stiffness at the floor is zero (above), and a floored committed state can sit **outside F by\n"
    "   O(p_min)** after a wet-side post floor [round 3b; wet_floor_har.py]: HAR TIMs set, a state at the floor on the surface\n"
    "   at η = 0.5M (θ = 0.454, off-corner), ψ_i ∈ {0, +0.04, +0.08, +0.12}, one backward-Euler pure-shear step along n̂ of\n"
    "   engineering shear strain Δγ := √2 ‖dev Δε‖ ∈ {10⁻⁴, 10⁻³, 10⁻²} (pattern -Pf throughout): F(σ_f)/p_min = +0.27\n"
    "   (Δγ 10⁻⁴, any ψ_i) … +0.82 (ψ_i +0.12, Δγ 10⁻²: p_c = −0.239 → −0.505, q 0.304 → 0.414, η_c = 1.327), against the\n"
    "   first-order size of item 2, 1.37 p_min, there; the Adversary's state gives +0.5 … +1.35. **The wet-side contraction\n"
    "   is self-limiting against the HAR domain edge**: the converged return moves ε^e_v up from the floor by |Δε^p_v| ≤\n"
    "   2.7e−5 over the table (1.2e−5 at Δγ 10⁻⁴, 2.0e−5 at 10⁻³, 2.7e−5 at 10⁻² and ψ_i +0.12), below the margin\n"
    "   ε*_f(ε_s) = (p_min/p_a)/(k(1−n) x^n) = 7.49e−5 at ε_s = 0 and 7.3e−5 at the state (it shrinks with shear: 1.75e−5\n"
    "   at the K1.14 ε_s = 2e−4), so the converged state never leaves dom Ψ (the Adversary: |Δε^p_v| ≤ 3.5e−5 < 7.5e−5,\n"
    "   no refusal up to Δγ = 10⁻² at ψ_i ≤ +0.12). The full Newton steps do overshoot the edge (6–37 rejected iterates per\n"
    "   step at Δγ ≥ 10⁻³, min ε* down to −1.2e−3) and the line-search backtrack of §2.4 absorbs them: no refusal in the\n"
    "   table. Nothing else changes:")
rep("M-F9 floor applied inside the local Newton: the K1.12 (C) values differ (the post-floored q is not the in-Newton\n"
    "one).\n",
    "M-F9 floor applied inside the local Newton: the K1.12 (C) values (BA06, wet side) and the K1.14b FPf values (HAR, dry\n"
    "side) differ — Δλ, π_i and q_c are those of the unconstrained return; a floor inside the Newton converges to a\n"
    "different iterate (p_c pinned at −p_min) and, under HAR, would also hide the domain overshoots that the line search\n"
    "is meant to backtrack (item 3).\n")

# ---------------------------------------------------------------------------------------------- §10.2 (S.56)
rep("      W_ramp := 1 − π_i(η₂)/π_i(η₁) = 1 − [ (1 − c₂N)/(1 − c₁N) ]^{(1−N)/N}   (N > 0),    = 1 − e^{−(c₂ − c₁)}   (N = 0);\n"
    "      contract:  PI_SCAN_REL ≤ W_ramp/10,  enforced in validate() (refuse a narrower smooth cap). Defaults: W_ramp = 0.0606;\n"
    "      c₁ = 0.05 with c₂ = 0.07 gives 0.0122 (admissible), with c₂ = 0.06 gives 0.0061 (refused).                      (S.56)\n",
    "      W_ramp := 1 − π_i(η₁)/π_i(η₂) = |π_i(η₁) − π_i(η₂)| / |π_i(η₂)|   (at fixed p; |π_i(η₁)| < |π_i(η₂)|)\n"
    "             = 1 − [ (1 − c₂N)/(1 − c₁N) ]^{(1−N)/N}   (N > 0),    = 1 − e^{−(c₂ − c₁)}   (N = 0);\n"
    "      contract (cap = smooth only, c₁ < c₂):  PI_SCAN_REL ≤ W_ramp/10,  enforced in validate() (refuse a narrower\n"
    "      smooth cap); planar (c₁ = c₂, W_ramp ≡ 0) and no cap have no ramp and are not subject to it. Defaults: W_ramp =\n"
    "      0.0606; c₁ = 0.05 with c₂ = 0.07 gives 0.0122 (admissible), with c₂ = 0.06 gives 0.0061 (refused).            (S.56)\n"
    "\n"
    "(Round 3b: the round-3 text defined W_ramp as 1 − π_i(η₂)/π_i(η₁), which is −0.0645 on the K2 set; the displayed\n"
    "formula, the 0.0606 and the measured 6.060e−2 were the ratio written above all along — round3b_sympy.py (b)–(c). An\n"
    "ungated refusal would have rejected every planar and no-cap model, since W_ramp = 0 there.)\n")
rep("12 476 replays is consistent with the geometry. New requirement for O2/kernel `validate()`: the (S.56) refusal.\n",
    "12 476 replays is consistent with the geometry. New requirement for O2/kernel `validate()`: the (S.56) refusal, gated\n"
    "to cap = smooth.\n")

# ---------------------------------------------------------------------------------------------- §13 K1.12, K1.14b
rep("−50.00 elastic, Δε^f_v = 4.066e−3; (B) plastic from a floored trial, p_c = −50.70; (C) wet side, p_c = −48.76 → −50.00\n"
    "   by the post floor, η = 0.953.\n",
    "−50.00 elastic, Δε^f_v = 4.066e−3; (B) plastic from a floored trial, p_c = −50.70 (dry side, FP-: under BA06 α₀ = 0\n"
    "   the return from a floored trial always compresses p, §9.7 item 2); (C) wet side, p_c = −48.76 → −50.00 by the post\n"
    "   floor, η = 0.953 (-Pf). The dry-side pattern FPf exists under HAR only: K1.14b.\n")
rep("   1.45e−2 without it (M-F3b).\n"
    "15. **π_{i0} rule (K1.15; (S.53))",
    "   1.45e−2 without it (M-F3b).\n"
    "   **14b. HAR dry-side FPf and the chained patterns (K1.14b; (S.32f), (S.54) on the HAR law) [I; round 3b;\n"
    "   floor_fd_har.py]:** TIMs elastic set with M = 1.3309, ρ = ρ̄ = 0.71 (WW), fork CSL (e₀ 0.83, λ_c 0.027, ξ 0.45,\n"
    "   p_a 101) and the K2 plastic constants N 0.4, N̄ 0.2, χ −3.5, h 280 (TIMs' own are a P3 refit), p_min = 0.505 kPa,\n"
    "   p₀ := −p_a in the §9.1 scalings. Surface state p = −0.6, q = 0.6999 kPa, η = 1.2M = 1.5971, θ = 0.271, π_i =\n"
    "   −0.74366 kPa, v = 1.72704 (ψ_i = −0.10), ε_v = 9.851e−4, ε_s = 3.335e−5 (D₁₁ 13 525, D₁₂ −6 235, D₂₂ 28 255 kPa).\n"
    "   Increment Δε = (2e−5/3) 1 + 2e−5 n̂ (expansion + shear): p^tr = −0.4607 → floored, Δε^f_v,tr = 3.82e−6 → return\n"
    "   Δλ = 1.157e−5, η_c = 1.8200, π_i = −0.74335, Δε^p_v = +7.07e−6 (dilative), Δε^p_s = +1.57e−5, p_c = −0.4877\n"
    "   (D₁₁Δε^e_v = −0.0843 kPa, D₁₂Δε^e_s = +0.0997 kPa) → post floor p = −0.505, q 0.6619 → 0.6705, F(σ_c)/p_min =\n"
    "   6e−13, F(σ_f)/p_min = −0.0041 (inside). Tangent (S.32f) vs central FD over the three principal strains: 2.2e−5 /\n"
    "   2.2e−7 / 2.2e−9 at h = 10⁻⁶/10⁻⁷/10⁻⁸ (O(h²); the increment is 2e−5); δ:C_f = 2e−16; committed p = −p_min to\n"
    "   2e−15. Mutants: Φ^post dropped 0.14, Φ^tr dropped 0.12, both dropped (plain ã^{ep}) 2.2; v-column tied to the\n"
    "   floored trial 1.2e−6 (tr Δε = 2e−5 makes the v-column negligible here — M-F3c stays a K2 BA06 test). Chains\n"
    "   (S.54), m = 2, fractions ½, ½: twice the increment gives FPf,FPf (the increment itself halves to -P-,FPf), 2.3e−7 /\n"
    "   2.3e−9 at h = 10⁻⁷/10⁻⁸; tr Δε = 1.4e−5 with shear 2e−5 gives -P-,-Pf (unsplit: -P- with p_c = −0.5072; the\n"
    "   halved one crosses the floor in its second return), 1.2e−6 / 1.2e−8; Φ operators omitted from the chain 2.1 for\n"
    "   both. (The Adversary's independent reproduction of the same case: p_c = −0.490, +0.10 vs −0.085 kPa, 8e−9 / 5e−8\n"
    "   and 9e−9 / 6e−8 at its own h.) A test author runs it with O2's `elastic` swapped for the HAR law until O2 ships\n"
    "   `-energy HAR`.\n"
    "15. **π_{i0} rule (K1.15; (S.53))")

# ---------------------------------------------------------------------------------------------- §15
rep("| e₀ = 0.83, λ_c = 0.027, ξ = 0.45, p_a | e₀, λ_c, ξ, p_a of (S.22) | **direct** (same e_c(p) form; p_a must be the "
    "same numeric value TIMs used, e.g. 101.325 kPa) |",
    "| e₀ = 0.83, λ_c = 0.027, ξ = 0.45, p_a | e₀, λ_c, ξ, p_a of (S.22) | **direct** (same e_c(p) form; p_a is the one "
    "`-p_a` flag shared with the HAR energy, §2.4, and takes the numeric value TIMs used: **101 kPa** (campaign set "
    "`Patm 101`), not 101.325) |")

# ---------------------------------------------------------------------------------------------- §16.1
rep("   p′ floor: designed 2026-10-03 (§9.7, owner decision (b)); π_{i0} rule: §5.4 (owner decision (c)); the nested-scan\n"
    "   step is 10⁻³ (§10.2 (S.56)).\n",
    "   p′ floor: designed 2026-10-03 (§9.7, owner decision (b)); **option (i), the strain-space projection Π_f (S.48),\n"
    "   approved by the owner 2026-10-03** (the two solver/reporting questions of §16.3 stay open); π_{i0} rule: §5.4\n"
    "   (owner decision (c)); the nested-scan step is 10⁻³ (§10.2 (S.56), refusal gated to cap = smooth).\n")

# ---------------------------------------------------------------------------------------------- §16.3
rep("  and the F vs F/2 deck report (K5). Design decisions still open for the owner: whether the kernel reports E_f (needs\n"
    "  Ψ) or only W_f; the global-solver behaviour on a floored patch (zero bulk stiffness is the consistent tangent, §9.7).\n",
    "  and the F vs F/2 deck report (K5); under HAR the K1.14b FPf case and the chains FPf,FPf / -P-,-Pf (round 3b).\n"
    "  **Owner decision 2026-10-03: option (i) of §9.7, the strain-space projection Π_f (S.48), is APPROVED.** Two design\n"
    "  questions stay **OPEN for the owner**; the orchestrator's recommendation is recorded beside each and decides nothing:\n"
    "  (a) the global-solver behaviour on a fully floored patch, where the consistent tangent has zero bulk stiffness\n"
    "  (δ:C_f = 0, §9.7) — recommendation: the kernel keeps the exact consistent tangent; any added stiffness would be\n"
    "  TANGENT-ONLY, opt-in, default off, so it changes the convergence path but never the converged answer and leaves the\n"
    "  FD and parity gates valid (they test the exact tangent with the option off); (b) whether the kernel reports E_f\n"
    "  (needs Ψ) or only W_f — recommendation: W_f is always counted in the kernel; E_f is an on-demand response evaluated\n"
    "  from the closed-form Ψ ((S.4)/(S.4h)) at the committed state, never on the step path.\n")
rep("- **[round 3] validate() gains the (S.56) refusal** (PI_SCAN_REL ≤ W_ramp/10) in O2 and the kernel — the only code change\n"
    "  the scan-step ruling asks for; PI_SCAN_REL = 10⁻³ and PI_SCAN_MAX = 1000 stay.\n",
    "- **[round 3] validate() gains the (S.56) refusal** (PI_SCAN_REL ≤ W_ramp/10, **cap = smooth only**: with c₁ = c₂ or\n"
    "  no cap W_ramp = 0 and an ungated check would refuse every planar/none model — round 3b) in O2 and the kernel — the\n"
    "  only code change the scan-step ruling asks for; PI_SCAN_REL = 10⁻³ and PI_SCAN_MAX = 1000 stay.\n")

# ---------------------------------------------------------------------------------------------- §16.4
rep("- floor_fold.py (shipped O2 on the AMP_STOP smooth-cap path): the §10.2 (S.56) numbers (d₂ and d_x per n, the step-8\n"
    "  roots through the fold under both scan steps, the ramp widths and the (S.56) constant).\n",
    "- floor_fold.py (shipped O2 on the AMP_STOP smooth-cap path): the §10.2 (S.56) numbers (d₂ and d_x per n, the step-8\n"
    "  roots through the fold under both scan steps, the ramp widths and the (S.56) constant).\n"
    "\n"
    "Round 3b checks (2026-10-03, local py312 venv with the O2 oracle imported read-only; scripts and logs in\n"
    "`gitAPE/ladrunoNORDSAD/handoff/round3b/scratch/`, to be copied into `tests/scratch_sheet_floor/` with the sheet):\n"
    "- har_patch.py: the HAR law (S.5h)–(S.5h'') as a drop-in for O2's `elastic`/`energy_psi`/`invert_elastic` (the\n"
    "  plastic part is energy-blind, §2.1), the (S.50) floor closed form and Π_f under HAR; round3b_sympy.py (a): its\n"
    "  stress forms D₁₁, D₁₂ = D₂₁, D₂₂, q/ε_s vs direct differentiation of (S.4h) at three TIMs states 2.2e−31, the\n"
    "  (S.5h'') round trip 2e−35, the K1.13/K1.14 ε_{v,f} values to 4e−18.\n"
    "- round3b_sympy.py (b)–(d): (S.56) 1 − π_i(η₁)/π_i(η₂) = [(1 − c₂N)/(1 − c₁N)]^{(1−N)/N} on 8 random sets 7e−31,\n"
    "  the round-3 words 1 − π_i(η₂)/π_i(η₁) = −0.0645 on K2 against +0.0606, the N → 0 limit exact, the floor_fold (C)\n"
    "  ratio identical, W_ramp ≡ 0 at c₁ = c₂; the HAR floor margin ε*_f(0) = (p_min/p_a)^{1−n}/(k(1−n)) = 7.485e−5.\n"
    "- floor_fd_har.py (O2 return map on the HAR law, TIMs set, p_min 0.505 kPa, off-corner): the K1.14b record — FPf\n"
    "  2.2e−5 / 2.2e−7 / 2.2e−9 at h = 10⁻⁶/10⁻⁷/10⁻⁸, δ:C_f 2e−16, mutants 0.12–2.2 and 1.2e−6 (v-column); chains FPf,FPf\n"
    "  2.3e−7 / 2.3e−9 and -P-,-Pf 1.2e−6 / 1.2e−8, Φ omitted 2.1; fpf_region_har.py: the single-step pattern map at that\n"
    "  state (FE- / FPf / FP- / -P- over expansion 1e−5 … 3.2e−4 and shear 0 … 3.2e−4).\n"
    "- wet_floor_har.py (the same O2 on HAR, a state at the floor at η = 0.5M, ψ_i 0 … +0.12, Δγ 10⁻⁴ … 10⁻², one BE\n"
    "  step with the floor at the trial and after convergence): F(σ_f)/p_min +0.27 … +0.82, converged contraction\n"
    "  |Δε^p_v| ≤ 2.7e−5 < ε*_f, 0–37 domain overshoots per step backtracked by the line search, no refusal (§9.7 item 3).\n")

# ---------------------------------------------------------------------------------------------- §16.6
rep("decision 2026-10-01, exponential v-update); item 15 at the HAR round (2026-10-02); item 16 at the p-floor round\n"
    "(2026-10-03). Tags:",
    "decision 2026-10-01, exponential v-update); item 15 at the HAR round (2026-10-02); item 16 at the p-floor round\n"
    "(2026-10-03); item 17 at round 3b (2026-10-03, the Adversary's round-3 defects and the owner's approval of option\n"
    "(i)). Tags:")
rep("| [I; round 3], owner decisions 2026-10-02 | floor_sympy.py (all OK), floor_fd.py (all OK), floor_rate.py, floor_fold.py |\n"
    "\n"
    "Items 1–13:",
    "| [I; round 3], owner decisions 2026-10-02 | floor_sympy.py (all OK), floor_fd.py (all OK), floor_rate.py, floor_fold.py |\n"
    "| 17 | §1.3, §2.4, §9.7 (items 2–3, FD record, M-F9), §10.2 (S.56), §13.12, §13.14b (new), §15, §16.1, §16.3, §16.4 | "
    "**Round 3b (2026-10-03): the Adversary's five round-3 defects; owner approval of option (i).** (1) MAJOR — the "
    "'dry side: the return compresses p, F(σ_f) < 0; post floor = wet-side event' statement of §9.7 item 2, K1.12 (B) and "
    "M-F9 was BA06-only: under HAR the D₁₂ coupling of the plastic shear strain (+0.100 kPa) beats the volumetric "
    "compression of the dilative return (−0.084 kPa) and the dry-side pattern is FPf (TIMs set, p = −0.6 kPa, η = 1.2M, "
    "p_c = −0.488); reproduced independently, (S.32f) FD-checked on the HAR law 2.2e−7 / 2.2e−9 (O(h²)), the chains "
    "FPf,FPf 2.3e−7 / 2.3e−9 and -P-,-Pf 1.2e−6 / 1.2e−8 added (K1.14b); the dF ≈ F_p dp + ζ dq sign rule written for "
    "both energies. (2) MINOR — (S.56) words: W_ramp is 1 − π_i(η₁)/π_i(η₂) (the round-3 '1 − π_i(η₂)/π_i(η₁)' is −0.0645 "
    "on K2); formula and numbers unchanged. (3) MINOR — the (S.56) validate() refusal gated to cap = smooth (W_ramp ≡ 0 "
    "for planar/none). (4) MINOR — p_a is one `-p_a` flag for HAR and the fork CSL (not in the BA06 refusal list); the "
    "TIMs value is 101 kPa (campaign `Patm 101`; §15's 101.325 withdrawn; p_min 0.505 kPa and K1.13/K1.14 unchanged); "
    "the HAR→BA06 mutant is killed by 13.1h/13.11/13.13/13.14 and O2 parity, not by self-FD. (5) MINOR — §9.7 item 3 "
    "states the wet-side F(σ_f) magnitude, +0.27 … +0.82 p_min (first-order size 1.37 p_min; Adversary +0.5 … +1.35), "
    "and the self-limiting contraction |Δε^p_v| ≤ 2.7e−5 < ε*_f = 7.5e−5 (Adversary 3.5e−5), Newton overshoots of the HAR "
    "domain edge backtracked, no refusal. Owner: option (i) Π_f approved 2026-10-03 (§9.7, §16.1); the two §16.3 "
    "questions (E_f vs W_f; floored-patch solver behaviour) stay open with the orchestrator recommendation recorded "
    "(exact tangent kept, any added stiffness tangent-only/opt-in/default off; W_f always counted, E_f on demand from the "
    "closed-form Ψ). No formula changed. | [I; round 3b], owner approval 2026-10-03 | har_patch.py, round3b_sympy.py "
    "(all OK), floor_fd_har.py (all OK), wet_floor_har.py |\n"
    "\n"
    "Items 1–13:")
rep("algebra ((S.48)–(S.55)) and changes no existing formula; its one correction to existing text is the scan-step value of\n"
    "§8/§10.2 (10⁻⁴ → 10⁻³), a contract constant, not a derivative.",
    "algebra ((S.48)–(S.55)) and changes no existing formula; its one correction to existing text is the scan-step value of\n"
    "§8/§10.2 (10⁻⁴ → 10⁻³), a contract constant, not a derivative. Item 17 changes no formula: it corrects the words of\n"
    "(S.56), the energy-dependence of the §9.7 dry/wet statements, the p_a bookkeeping and the mutant-kill list, and adds\n"
    "the HAR FD record (K1.14b) and the wet-side magnitudes.")

# ---------------------------------------------------------------------------------------------- apply
missing = [old for old, _ in EDITS if text.count(old) != 1]
if missing:
    for old in missing:
        n = text.count(old)
        print(f"ANCHOR {'MISSING' if n == 0 else 'DUPLICATED (%d)' % n}: {old[:90]!r}", file=sys.stderr)
    sys.exit(f"refusing: {len(missing)} of {len(EDITS)} anchors do not match exactly once; nothing written")
for old, new in EDITS:
    text = text.replace(old, new, 1)
with open(OUT, "w", encoding="utf-8", newline="") as fh:
    fh.write(text)
print(f"applied {len(EDITS)} edits: {SRC} -> {OUT}")
