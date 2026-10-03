---
wp: ADR-66
title: "getResistingForceIncInertia MUST snapshot the shared static resid before calling getRayleighDampingForces() — else stiffness-proportional Rayleigh silently dro…"
legacy_seq: 148
---
## `getResistingForceIncInertia` MUST snapshot the shared static `resid` before calling `getRayleighDampingForces()` — else stiffness-proportional Rayleigh silently drops element inertia (ADR-66 P5.1)

An element that builds its residual in a **static class-member** `Vector resid` (the OpenSees
norm) and does `formInternal(); formInertia(); resid += getRayleighDampingForces();` has a hidden
re-entrancy bug: `Element::getRayleighDampingForces()` (Element.cpp:347/349) calls
`this->getTangentStiff()` (when `betaK != 0`) or `this->getInitialStiff()` (when `betaK0 != 0`),
which re-enter the element's own form routine and **`resid.Zero()`** it. The `resid +=` then adds
the damping force to a resid that has been wiped back to `f_int` (or, for the first `betaK0` call
before `Ki` is cached, to **zero**), so the returned unbalance is missing `M·a` — and Newton still
CONVERGES (the Newmark tangent keeps its `c3·M` term) to a wrong dynamics solution, or (as observed
for this element) fails to converge outright. **Fix = the LadrunoBrick donor pattern
(LadrunoBrick.cpp:688-700): `static Vector res(24); res = resid;` BEFORE the Rayleigh call, then
accumulate into `res`.** SYMPTOM: a transient with `rayleigh 0 <betaK> 0 0` gives quantitatively
wrong periods/accelerations (or diverges) while a `betaK=0` run is fine. GATE: compare a `betaK=0`
transient to a tiny-`betaK` one — light damping must change the peak <5%; the bug collapses the
tiny-`betaK` run to quasi-static (or non-convergence). This is INVISIBLE to every static/patch gate
— a dynamic Rayleigh regression is mandatory for any new element with mass. (Verified by
reverting the fix + rebuilding: the gate fails `analyze -3` on the buggy binary, passes on the fixed.)

**2026-07-11 recurrence (ADR-70 P4a adversarial gate): the WHOLE plane family had it.**
`LadrunoQuad`/`LadrunoCST`/`LadrunoLST` under `-geom finite` and `LadrunoCSTPair` all used the
unsnapshotted pattern; their `getTangentStiff()` → `formFinite`/`formPair(1)` refills the shared
static `P` (small-strain paths are safe — those fill only `K`). Extra sting on the pair: the
re-entry also wiped `−Q`, so UniformExcitation ground-motion loads were dropped too. All four fixed
with the donor snapshot in the same PR; parametrized regression gate
`tests/test_ladrunoplane_dynamics.py::test_dynamic_rayleigh_preserves_inertia[pair|cst|quad|lst]`.
LESSON: the quirk was already in this ledger and the P4a author even reasoned about it in a code
comment — and still missed that the clobber happens INSIDE `getRayleighDampingForces()`. When a
quirk names a pattern, grep for the pattern (`P += this->getRayleighDampingForces`), don't reason
about the instance.
**`Vector::pNorm(0)` is NaN-BLIND — never use it for a divergence/NaN check (2026-07-02).** `pNorm(0)`
implements max via `value = (fabs(data) > value) ? fabs(data) : value`; every comparison against NaN is
FALSE, so NaN entries are silently SKIPPED and the returned max is never NaN. The explicit integrators'
circuit breakers (`A_max = U.pNorm(0); if (A_max != A_max || A_max == inf)`) therefore only ever fired on
±Inf: an all-NaN acceleration (NaN material state, poisoned IC, `inf − inf` residual) sailed through and
was COMMITTED into the node state. Fixed by `vectorIsFinite()` (`CriticalTimeStep.{h,cpp}`, std::isfinite
scan) in `CentralDifferenceLadruno`/`ExplicitBathe`/`LadrunoDynamicRelaxation::update()`. Related trap when
WRITING the test: a *seeded* Inf displacement degenerates to NaN before reaching the breaker
(`inf + dt·(−inf) = NaN`), so the honest Inf-path fixture is a genuinely divergent run driven to its first
overflow, and the honest NaN fixture is a NaN committed displacement (`Fint = k·NaN`). See
`tests/test_explicit_nan_breaker.py`.

**betaKinit/betaKcomm MUST be summed with betaK for any explicit stable-step bound (2026-07-02).**
`C = αM·M + βK·K + βK0·K_init + βKc·K_commit` — all three β slots are stiffness-proportional and shrink the
explicit stable step identically (ξ = β·ω/2 grows with ω), and `rayleigh 0 0 βKinit 0` is the MOST COMMON
form in practice (chosen precisely to avoid the current tangent). Reading `getRayleighDampingFactors()(1)`
alone made the SMS damped sizing AND the damped dt_cr estimate blind to it → under-scaled → unstable at
dtTarget for exactly the users the betaK feature targets. Use `stiffnessRayleighBeta(ele)`
(`CriticalTimeStep.{h,cpp}`): per-slot clamp at 0, then sum. Exact at the initial state (K == K0 == Kc);
conservative under softening — the right side to err on for a stability bound.

**The `-divergence` KE proxy FALSE-TRIPS on free vibration at velocity troughs (2026-07-02; FIXED review-P3 #475 -- baseline is now the running MAX of KE).**
The breaker compares per-step `ke/prevKE` against the factor, with `prevKE` updated every step it is
positive. In plain free vibration the velocity passes through ~0 every half period; the step nearest the
zero leaves `prevKE ≈ ε`, and the next steps' quadratic KE regrowth off that near-zero floor produces an
unbounded ratio — phase luck decides whether a given trough exceeds the factor (observed: an SDOF at
dt=0.005, ω=10 tripping factor 10 at its second trough). Workaround in tests/models: excite a
CIRCULAR-motion state (equal springs x+y, quadrature seed) whose total KE is constant, or keep
`-divergence` for monotonic-divergence detection only. Real fix (PR-3 diagnostics batch candidate):
compare against a running MAX of KE, not the previous step.
**ADR-63 P2.5 — the AUTO outward sign is a per-slave MAJORITY vote; a LOCAL closest-point vote beats the
aggregate-normal·global-chord coin-flip (2026-07-01).** The frozen global sign on the auto (no-`-outward`)
path used to be `sign(Σ_a σ_a n_a · (slaveCentroid − masterCentroid))` — an AGGREGATE normal (which nearly
cancels on a curved/domed master) dotted against a GLOBAL chord (which goes ~tangent to the field when the
slave cloud grazes the master edge-on) ⇒ `vote·seed ≈ 0` is a coin-flip and a tiny wrong-signed component
flipped the WHOLE field inward → silent pass-through even with `-smoothNormal` (F2/F3/F5). FIX =
`LadrunoContactProjection::voteSignRobust`: each slave projects onto its NEAREST facet (closest-point,
clamped) to get that slave's LOCAL coherent unit normal `n̂ = σ_{s*}·newell̂(s*)`, then votes
`w = n̂·(slave − surfaceCentroid)` (surfaceCentroid = slot-average of the facet nodes, an INTERIOR
reference for an open convex patch); the surface takes the DISTANCE-WEIGHTED majority `sign(Σ w)`. TWO
adversarial-forced choices: (i) the LOCAL normal (not the aggregate `Σσn`) supplies the lateral
component the aggregate lacked at an edge-grazing slave — the actual F2/F3/F5 fix; (ii) the INTERIOR
CENTROID reference (not the local footpoint separation) keeps a single slave seeded slightly
PENETRATING voting outward — a footpoint separation points inward for such a slave and with one slave
there's no majority to protect it (the P1 sign gate does exactly this; adversarial F2). `|w|`-weighting
then lets a clearly-separated majority dominate. The vote runs on the REFERENCE coords of BOTH master
and slaves (adversarial F1 — the DEFORMED master vs reference slave mix mis-signs on restart / mid-run
recapture; `setNormalField` takes `refSegCoords`, the DEFORMED `segCoords` still drives the per-handle
field). RESIDUAL (LOW): a non-convex open patch whose centroid is not interior ⇒ pass `-outward`.
KEY POINTS: (1) it decides only the
ONE frozen sign (D2/F1) — still captured once, still frozen; (2) `-outward` given ⇒ the aggregate
`sign(vote·outward)` path is UNCHANGED (byte-behavior preserved) — the robust vote runs only on the auto
path; (3) `nVoted==0` (no slave projected) ⇒ fall back to the aggregate seed vote; (4) a genuinely two-sided
cloud yields margin≈0 ⇒ the SAME `conf<0.1` handler warning fires (the ambiguity is DETECTED, recommend
`-outward`); (5) a disconnected multi-shell master is still refused at `propagateOrientation` — a
per-connected-component vote (run `voteSignRobust` per component vs its own nearest slaves) is the
Q-MULTISHELL follow-up; (6) the slave coords fed to the vote are the REFERENCE coords (config-independent —
captured once); (7) RESIDUAL: the degenerate-BLEND fallback still orients by the aggregate seed (review
GAP-2), so a degenerate blend AND an edge-grazing cloud together can still drop a pair (fails safe) — pass
`-outward` for that compound corner. `-smoothNormal` OFF stays byte-identical; no classTag; no vanilla touch.

**An ABSOLUTE tolerance on a DIMENSIONAL residual silently killed all contact away from the origin
(2026-07-02, contact-review fix PR-1).** `LadrunoContactProjection::project()` converged on
`|R| < 1e-12` where `R = d·g_α` has units **length²**: its floating-point noise floor is `~eps·|X|·|g|`,
so for coordinates far from the origin the test was UNREACHABLE — e.g. a plain mm-unit building
(h~500 mm facets at x~5e4 mm, noise ~3e-10) failed **200/200** projections; every pair evaluated
inactive and contact SILENTLY vanished (slave free-falls, no warning). Never caught because every
contact gate ran at unit scale near the origin, where the products happen to be exact in binary.
FIX = a scale-free **parametric-step escape** `|dξ|+|dη| < 1e-8` (parent coords are O(1)) checked AFTER
the Newton update. HONEST behavior contract (the adversarial gate REFUTED a stronger "bit-preserving"
claim): NO previously-converging input is ever LOST (0/5M trials), and on a FLAT facet the escape exits
one iteration early with the IDENTICAL (ξ,η); on a WARPED facet Gauss-Newton contracts only linearly, so
~19% of converging warped cases exit with a footpoint within ~tolP parametric of the residual-converged
answer (measured drift ≤ ~1e-9; gap error second-order ⇒ physically nil; the full 199-test contact+tie
battery, incl. its exact-`==` byte-identity gates, is unchanged) — in-bounds classification can flip ONLY
for a footpoint within the 1e-9 parent-boundary slack (~1e-10·h of slave positions). Oracle:
`contact_prototypes/proto_projection_offset.cpp` (13/21 checks fail on the pre-fix header); in-solver:
`tests/test_contact_projection_offset.py`. THREE general lessons: (1) any tolerance compared against a
quantity with length units must be RELATIVE to a local metric (the `detK` degeneracy guard had the twin
disease — length⁴ vs a length² floor — now the pure angle test `detK < 1e-14·K00·K11`, i.e. sin²θ);
(2) the SAME absolute test is too LOOSE at micro scale: for h≲1e-5 it passes at the INITIAL GUESS
(R ~ |d||g| < 1e-12 before any iteration), so micro-facets get centre-footpoint "convergence" — gap on a
flat facet is footpoint-independent so contact stays functional, just parametrically sloppy (documented,
unchanged; tightening it would break byte-identity for nothing); (3) test meshes at unit scale near the
origin CANNOT catch dimensional-tolerance bugs — put one offset/scaled case in every geometric oracle.
Same review pass: the bucket-grid cap arithmetic used 32-bit `long` (LLP64 Windows) — a diverging
`clipPct=0` re-emit feed with per-axis cell counts in the [~5e4, 2e9] window wrapped the product, exited
the cap loop early, and `grid_.assign()` could request a ruinous allocation (`bad_alloc` kills the
process mid-run). FIX = per-axis pre-clamp to 5000.0 IN DOUBLE before the int cast (also removes the
double→int UB; the total is capped at min(nSeg,5000) cells anyway) + `long long` product arithmetic.
**OPEN FOLLOW-UP (found by the PR-1 adversarial gate, pre-existing, NOT fixed here):** the mortar
back-map `LadrunoMortarKernel::inverseIsomap2D` has the SAME disease — `tolR = 1e-13` ABSOLUTE on a
LENGTH-unit residual (aux-plane UV ~ facet size h, noise floor ~eps·h): measured 0/400 GPs dead at
h≤500 but **241/400 dead at h=5e3 and 378/400 at h=5e4**, and the caller silently SKIPs a failed GP ⇒
`-mortar` contact and `LadrunoTie -mortar` quietly lose most of their integration for facet edges
≳ 2000 length units. Same fix pattern (parametric-step escape or eps-relative tolR) + its own oracle —
file with the review-fix PR-3 hardening batch. Same-theme, likely benign: `isConvex2` tol=1e-12
(length²) and `dedupe` tol=1e-12 (length) in the same header.
**Collapsing a command family into flags? DIFF THE DEFAULTS, not just the grammar (2026-07-02).**
The ADR-52 W1-E2 collapse kept every deprecated alias parsing byte-compatibly, but the UNIFIED command
carried its own `-lump` default (RowSum, upstream-compatible for the bare dt_cr estimate) while the
alias impl defaulted Diagonal ("matches the system Diagonal run") — so `ExplicitBathe p -sms dt` and
`ExplicitBatheSMS p dt` sized the SAME model's scaling with DIFFERENT lumping: 3.73 vs 31.88 injected
mass (8.5x) on a consistent-mass beam, silently. Per-combo byte-identity tests passed because the test
elements had rowsum == diagonal (Truss — the [[project_zonea_link_blocker]] CDL battery caught that
equivalence once before). Fixed (review-P2): `-sms` without an explicit `-lump` flips the sizing
default to Diagonal. LESSON: when merging commands, enumerate every DEFAULT each retired parser had and
prove the merged parser reproduces them per mode — and put a rotational-DOF (rowsum≠diagonal) model in
the byte-identity battery.
