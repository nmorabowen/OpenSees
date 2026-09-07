# ADR-95 ASD-DP cross-check — VERDICT (2026-09-07, merged build feb358fda, staged dist_p5)

**The trigger is implementation-independent; the defect is not.** Same deck, `ASDPlasticMaterial3D`
DruckerPrager (η = 3α, ξ_c = SY/√3, ψ = 0; NO tension cutoff exists in `DruckerPrager_YF`, and its
`CHECK_APEX_REGION`/`APEX_STRESS` are stubs — ADR-94 open item):
- `h8bbar` linear control: **TARGET s/B 0.15, q 147.40 = 1.0611 of exact** (UW-DP on the same mesh:
  1.0850 — a 2.2 % cone-fit difference between the two implementations, not a wall question); 0 GPs
  with mean stress ≥ 0 anywhere on the path, 0 NaN.
- `h20uri`: **FLOOR at s/B 0.01804, q 122.45 = 0.8815 (ALLOWANCE)**, with the global tangent HEALTHY
  at the wall (σ_min/scale 2.38e-4, cond 8.5e4) and **exactly 4 GPs at mean stress ≥ 0 (I1 max
  +0.77 kPa) appearing at the last two stations and never before** — the same four footing-edge GPs
  that hit UW's cutoff, now hitting ASD's cone APEX, where the return map is a stub and the ladder
  collapses without a tangent event. No NaN (the ADR-94 B4 fix is in this build).
⇒ Both fork DP implementations fail at the first tensile Gauss points beside the footing edge, each
through its own defect: UW via the dead corner/cutoff branch (FIXED, #803), ASD via the unimplemented
apex return (OPEN, ADR-94). ASD walls later (0.018 vs 0.011) because without a cutoff I1 rises
further before the apex is reached. Legs: `asd_h8bbar.log`, `asd_h20uri.log`, `asd_path_diag.py`.

---
(previous blocked-state record follows)

# ADR-95 ASD Drucker-Prager cross-check -- BLOCKED before either leg ran

**VERDICT: neither leg (h8bbar control, h20uri) produced a usable path.** Both
die to NaN within ~10 converged steps (s/B <= 0.00064, ~0.4% of the s/B=0.15
target) as soon as the material goes plastic -- including the h8bbar CONTROL
leg, which under the vanilla UW material stays elastic/mildly plastic far
longer and plateaus cleanly. MODE for both = FLOOR (step collapses hunting a
converged state that no longer exists). This is NOT the quadratic-wall
question being reproduced or refuted -- neither leg got anywhere near the
corner/tension-cutoff regime note 95 studied. **Question left open.**

## Parameter mapping (verified correct)
UW cone: `||dev sigma|| + rho*I1 = sqrt(2/3)*SY` (I1 = tr(sigma), tension+).
ASD cone (`DruckerPrager_YF.h`): `sqrt(J2) + eta*p = xi_c`, p = tr/3, also
tension-positive (R4 in `reviews/adr94_verdict.md`: the header comment saying
compression-positive is wrong). Matching coefficients: **eta = 3*alpha**
(rho=sqrt(2)*alpha cancels the sqrt(2)), **xi_c = SY/sqrt(3)**; both cones'
apex I1 = sqrt(2/3)*SY/rho match to fp precision (asserted in-script).
`DP_etabar = 0` (psi=0 non-associated). `E = 9*K*G/(3*K+G)`, `nu = 0.45` at
the same K,G as the UW leg. Verified with SY forced to 1e6 (forces elastic):
1-D patch test max error 1.5e-14 -- mapping is exact, not the blocker.

## Root cause of the blocker
`DruckerPrager_YF.h:64-65` and `DruckerPrager_PF.h:68-69` (as they exist in
this worktree, `wp/95-prandtl-bezier-root-cause`):
```
VoigtVector pressure_part;
pressure_part *= 0.0;   // multiplies UNINITIALISED Eigen storage by zero
```
This is **ADR-94 finding B4**, already root-caused and fixed on `origin/ladruno`
(`wp/94a-fail-loud`, merged `a493b15c9` / PR #809, `dev_part.setZero()` /
`VoigtVector::Zero()`). The wp/95 branch was cut from `ladruno` at
`76668332f` -- before that merge -- so the P4-staged binary (`ADR95_DIST`,
build `0388889713bcdb7f6bb1a4bc6fac09f5b4e1d67a`) was compiled from the
**unfixed** source. Every Backward_Euler plastic-corrector iteration calls
this derivative; on this Windows heap it reliably decodes as NaN, silently
committed (`analyze()==0`) -- exactly B2+B4's documented signature.

## Does ASD-DP have a tension cutoff?
**No.** One yield surface only: `f = sqrt(J2) + eta*p - xi_c`. No f2. Its
`CHECK_APEX_REGION`/`APEX_STRESS` (`DruckerPrager_YF.h:105-121`) are stubs
("Implement!!!", return `false` / zero vector) -- so even the framework's
apex special-case never activates; nothing stands between the cone and B4.

## Why this is not a rebuild
Per instructions, no C++ edit/rebuild was performed; no other staged dist
(`dist_fixed`, `dist_bin_p0` = build `cf239c9d...`) carries the wp/94a fix
either. Re-running this script against a binary built past `a493b15c9` is
the natural next step and should reach the corner/tension regime cleanly.

## PR #815 (wp/94c) read against this deck — PREDICTION written before the re-run (2026-09-07 16:35)

#815 makes the apex site live: `check_apex_region` = `(p − p_apex) ≥ η·q` (Euclidean normal-cone
test in the (p, √J2) plane), `apex_stress` = hydrostatic p_apex = ξ_c/η, projection
dε^p = C⁻¹(σ_trial − σ_apex) guarded by |f(σ_apex)| ≤ tol. Its own comment states the caveat: the
exact test lives in the ELASTIC metric, `(p − p_apex) ≥ (K·η̄/G)·q`, and the Euclidean one coincides
with it only when K·η̄/G = η.

**On this deck η̄ = 0 (ψ = 0) and ν = 0.45 (K/G ≈ 9.7).** The exact condition degenerates to
`p ≥ p_apex`: with zero dilatancy the flank return cannot move p at all (the same fact that made
UW's over-cutoff state persist), so EVERY trial state beyond the apex belongs to the apex return.
The Euclidean test instead sends any such state with q > η·(p − p_apex) to the flank map, which
then cannot close f. The #815 guard only catches an apex whose own f is off-surface; a
misclassified flank return is not caught. The footing-edge GPs arrive at p ≈ p_apex with small but
nonzero q, i.e. exactly in the misclassified wedge.

**Prediction:** the ASD-DP `h20uri` leg on a #815 build still walls at s/B ≈ 0.018 with a healthy
tangent and ~4 GPs at mean stress ≥ 0; the linear control still plateaus (never reaches the apex).
If instead it plateaus, the Euclidean wedge is narrower than argued here and the prediction is
withdrawn. **Fix if the prediction holds:** classify in the elastic metric using K, G and the
potential's dilatancy (the PF is available to the integrator even if not to the YF signature), or,
signature-free, order the returns "flank first; if the flank map fails and p > p_apex, apex
projection" — the apex projection already adds the volumetric plastic strain a ψ = 0 cone lacks.
