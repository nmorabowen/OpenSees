# ADR-94 wp/94f — ASD-DP apex fallback for zero-dilatancy decks — VERDICT

**FIXED, and the ADR-95 ASD leg now reaches the target.** Build
`ff47275fd4a463d13a721a2138b48024a72463f2` (branch base; staged
`scratchpad/dist_94f/bin`, `ladrunoBuild()` verified). PR
[#832](https://github.com/nmorabowen/OpenSees/pull/832), branch `wp/94f-asd-apex-fallback`.

| gate | result |
|---|---|
| `tests/test_adr94f_asd_apex_fallback.py` | **4/4** — and **1 failed / 3 passed on the #815 build** (`c945f9a8b`), so it gates the fix |
| `tests/test_adr97_p4_inertness.py` | **10/10** — `Backward_Euler` still BYTE-IDENTICAL on all 23 baseline decks (baseline `3622d6214`, pre-94f) |
| `tests/test_adr94c_numerics.py` + `test_adr94_redblue_numerics.py` | **11/11** |
| ADR-95 `asd_path_diag.py --elem h8bbar` (control) | **MODE=TARGET**, s/B 0.15000, q_max **150.7104 = 1.0850** of exact 138.907 — unchanged from #815 |
| ADR-95 `asd_path_diag.py --elem h20uri --cond` | **MODE=TARGET**, s/B **0.15000** (pre-fix FLOOR **0.01122**), q_max **135.5529 = 0.9758** vs the repaired UW-DP **0.9757** |

`scalar Newton exhausted` lines on the `h20uri` leg: **0** (pre-fix **435**). No
refusal of any kind, no NaN, no apex-guard message. Final `[tension]` census:
`GP seen 1600, tension(mean>=0) 6, NaN 0, I1 [-276, +0.777]` — the footing-edge
Gauss points still go into tension, they are simply returned to the apex now
instead of walling the leg. Legs: `Ladruno_files/testbed/hypo_bearing/asd_{h8bbar,h20uri}_94f.log`,
`asd_*_h1.0_94f{,_tension}.csv`.

## The defect

`DruckerPrager_YF::CHECK_APEX_REGION` is EUCLIDEAN, `(p − p_apex) ≥ η·q`. The exact
test is in the ELASTIC metric, `(p − p_apex) ≥ (K·η̄/G)·q`, which at η̄ = 0 (zero
dilatancy — the ADR-95 deck) degenerates to `p ≥ p_apex`, because a non-dilatant
flank return has a traceless flow direction and cannot move `p` at all. Trials in
the wedge `η·q > (p − p_apex) > 0` were therefore handed to a flank map with no
solution (`f = q + η·p − ξ_c ≥ η·p − ξ_c > 0` for every `q ≥ 0`); its scalar Newton
exhausted and, under `strict_convergence`, the step was refused. The #815
`|f(σ_apex)| ≤ tol` guard cannot catch this — a *misclassified* state never reaches it.

## The fix (both layers in the integrator; the YF signature cannot see K, G, η̄)

**(a) Elastic-metric pre-check.** `Backward_Euler` unions the YF's Euclidean answer
with `cp_apex_region` — ADR-97 wp/97b's family-agnostic elastic-metric classification
(linearised cone step; APEX iff the returned deviator flips sign), reused verbatim, so
no dilatancy accessor was needed. Gated by the new `yf_apex_elastic_metric` trait,
specialized **only for `DruckerPrager_YF`**; every other yield function keeps its own
Euclidean answer and is byte-identical here.

**(b) Flank-first apex fallback.** At exactly the sites that would otherwise return
`LADRUNO_MATERIAL_REFUSED` — Newton exhaustion under `strict_convergence`, singular
local tangent, NaN, the strict plastic-inconsistency branch (the first three converted
from early `return`s to a `be_flank_failed`/`be_flank_reason` break) — a trial whose
mean stress exceeds `apex_stress().meanStress()` rewinds stress / plastic strain /
internal variables to the elastic predictor and takes the apex projection under the
same `|f(σ_apex)| ≤ tol_yf` guard. Otherwise the identical refusal is issued.
**Which YFs the fallback can touch:** those declaring `yf_has_apex` —
`DruckerPrager_YF`, `MohrCoulomb_YF<NO_HARDENING>`, `HoekBrown_YF<NO_HARDENING>`,
`TensionCutoff_YF<NO_HARDENING>`. For the latter three it can only convert a refusal
into an admissible vertex state; non-strict decks and every YF without an apex are
untouched (proved by the 23-deck inertness gate).

The wp/94c projection body is factored into `be_apex_project` and shared by both sites,
unchanged. On the ADR-95 deck layer (a) does all the work — the fallback message never
printed — but (b) is what makes the ordering robust to a wrong classification.

## Files

`SRC/material/nD/ASDPlasticMaterial3D/{ASDPlasticMaterial3D.h, YieldFunctionBase.h,
YieldFunctions/DruckerPrager_YF.h}` (three `LEDGER_vanilla_files` rows, #832) +
a `LEDGER_quirks` entry. All edits marked `// Ladruno (ADR-94 wp/94f)`.
