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
