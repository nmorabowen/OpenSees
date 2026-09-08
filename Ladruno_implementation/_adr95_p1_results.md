# ADR-95 P1 — VERDICT (interim, 2026-09-07, build db19a30e)

**H1 CONFIRMED in its precise form; H2 (strong) and the flicker reading are DEAD; H3 not needed.**

- `h20uri h0=1.0 --branch`: MODE = FLOOR at s/B 0.01114, q 0.7689 of exact (ALLOWANCE). Up to the
  last healthy station (s/B 0.011125) there are **zero** f2 or corner GPs, `n_forced = 0`,
  σ_min/scale 2.2e-4 flat, cond 9e4. At the wall station: **4 GPs on the cone+cutoff CORNER
  (branch 3) for the first time**, at x = ±2.638, z = −1.5 (elements 79/129, the first row beside
  the footing edge), with I1 = 0.80–0.98 kPa ≥ T = √(2/3)·SY/ρ = 0.816 kPa (f2 = +0.03…+0.21),
  and the normalised acoustic determinant at EXACTLY those 4 GPs is −3e5 … −6.6e7 (elsewhere
  −0.06). Every other GP is unchanged. ⇒ the abrupt 1000× σ_min collapse of note 82 §7.3 is the
  UW DruckerPrager **corner-branch consistent tangent**, which is pathological (entries ~1e2–1e3
  × 2G), entering the global operator through the first GPs to reach the tension cutoff.
- `h8bbar h0=1.0` control: MODE = TARGET, CAPACITY, 1.0850 of exact at this coarse h0 (matches the
  gate's coarse sequence). **Never** reaches f2/corner anywhere over s/B 0→0.15; I1_min −263.
- **detAmin < 0 at every plastic GP in BOTH elements** (non-associated ψ=0, ν=0.45): the linear
  element completes the collapse with ~400 non-elliptic GPs ⇒ loss of ellipticity per se is NOT
  the event (H2 strong form dead). Branch-histogram flicker (±40 % quadratic, but ~280 switches
  per station in the linear too) is a sampling artefact of perfectly plastic GPs resting on the
  surface (|f1| < 1e-3 at ~750 GPs) — NOT a discriminator.
- Forensics `top_unbalance` snapshot is the surcharge reference pattern (identical for both
  events) — the harness's unbalance capture is not usable; ignore.

**Why the linear bbar hex never hits it:** its element-averaged volumetric strain and coarse
resolution keep the heave-zone mean stress below the cutoff; quadratic elements resolve the tensile
spot beside the footing edge. Consistent with TIMs' dilatant (ρ̄=ρ) leg walling EARLIER (dilatancy
drives I1 up) and with DR walking through (no tangent formed).

**P2 launched (knob = SY 0.2 → 2 → 20 kPa ⇒ T 0.82 → 8.2 → 82 kPa):** prediction — the quadratic
wall moves OUT in s/B (possibly to a plateau); the linear leg keeps plateauing. q is NOT comparable
across SY (SY adds cohesion to the oracle); the observable is reach and termination mode.
Logs `p2_*.log`, CSVs `qpd_*_p2sy*`.

**CAVEATS from the P0 report (received after P1 ran):** (1) build db19a30e carried a wrong 1/2
factor on the shear columns of the Voigt→4th-order map, so **P1's `detAmin` MAGNITUDES are
invalid** (slots 0–6 — branch, gammas, f-trials, forced, I1 — are sound); the 4-GP outlier is
re-measured by the repeat leg `_p1rep` on the corrected pyd (build cf239c9d, staged copy, log
`p1_h20uri_rep.log`), which also serves as the rule-mandated n=2 repeat. (2) **`gamma(1)` is
structurally zero — the `Jact(i)==2` residual arm is dead code**: crossing the cutoff swaps the
tangent to the corner operator but returns NO stress, so I1 stays above T (P0 sentinel test:
15.0 vs T=1.099). That is the defect in one sentence: a tangent that does not correspond to the
return actually performed.

Open: cause-vs-effect is settled by the 5e-6 station spacing (corner absent at 0.011125, present at
the last converged state 0.0111375); the definitive falsifier is a `-noTensionCutoff` option (P4).
