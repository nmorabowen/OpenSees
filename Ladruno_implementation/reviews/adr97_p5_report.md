---
title: ADR-97 P5 — experimental_integrator gate (D5) + StiffSoilShear NaN root cause (M9) (report)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P5 (`wp/97f-explicit-gate`) — implementation report

**PR** [#827](https://github.com/nmorabowen/OpenSees/pull/827) (draft, based on
`wp/97d-cp-hoekbrown` = P3 [#825](https://github.com/nmorabowen/OpenSees/pull/825))
· **plan** [[97_ladruno_asdp_closest_point_adr]] · **pre-investigation**
`p5/findings.md` (session scratchpad; the root causes, the fix scripts, and the
draft test below were all produced there before this session opened a fresh
build).

## 1. What shipped

Two independent pieces on the same build, per the plan's P5 row.

### A. D5 — the `experimental_integrator` gate

A new `Begin_Integration_Options` token, `experimental_integrator` (int,
default 0), gates the four explicit integrators —`Forward_Euler`,
`Forward_Euler_Subincrement`, `Modified_Euler_Error_Control`,
`Runge_Kutta_45_Error_Control` — behind an opt-in. ADR-94 M2/M8 found real,
unfixed defects in them (empty drift-correction `if` bodies, an error
controller that accepts unconditionally at `rk45_dT_min`); rewriting four
integrators is a different ADR, so this is the *gate*, not the fix — the same
treatment ADR-94 gave `Backward_Euler_LineSearch` and
`Runge_Kutta_45_Error_Control_old`, one notch softer because these four have
real users (opt-in, not outright refusal).

The check is resolved **after** the whole `Begin_Integration_Options` loop,
in the same place as the existing `Algorithmic`/`Closest_Point` cross-check,
so `experimental_integrator` may appear before or after `integration_method`
in a deck. `Backward_Euler` and `Closest_Point` are untouched.

### B. M9 — StiffSoilShear step-1 NaN, two independent bugs, both fixed

**Bug A — `StiffSoilShear_YF.h`, unguarded `cot(phi)` at `phi == 0`.**
`qf = (c*cot(phi) + sigma3) * 2*sin(phi) / (1 - sin(phi))`. `cot(0) = Inf`;
`Inf * 2*sin(phi)` (`== Inf*0` at `phi==0`) is the IEEE-754 indeterminate
NaN. `phi == 0` is a normal, legitimate cohesive-only / undrained-clay
choice, not an edge case. Fixed by an exact algebraic rearrangement
(multiply through by `sin(phi)`):

```
qf = 2*(c*cos(phi) + sigma3*sin(phi)) / (1 - sin(phi))
```

numerically identical to the old formula for any `phi != 0` (matched to
~1e-14 relative in a standalone g++ probe) and gives the correct Tresca limit
`qf -> 2c` as `phi -> 0` instead of NaN.

**Bug B — `StiffSoilShear_PF.h`, unguarded normalization at an exactly
hydrostatic trial stress.** `PLASTIC_FLOW_DIRECTION` central-differences a
Mohr-Coulomb-type potential and does `vv_out /= norm` with no zero-guard. At
an exactly hydrostatic trial stress `computeMobilizedDilatancy()` returns
`psi_m == 0` exactly (its `sigma1 - sigma3` numerator is exactly zero), which
kills the potential's `I1*sin(psi)/3` term; the remainder,
`cos(lodeAngle)*sqrt(J2)`, is an EVEN function of the sign of a `+-ds`
perturbation, so the central-difference numerator is the exact zero vector on
all six Voigt axes — `norm == 0.0` and `vv_out /= norm` was the indeterminate
0/0. Fixed with a zero-guard:

```cpp
if (norm > 1e-12) {
    vv_out /= norm;
} else {
    vv_out.setZero();
}
```

the same `setZero()` precedent the ADR-94 DruckerPrager `NaN*0` fix used.

**This IS the ADR-94 R3a "6/194 non-finite cloud points" finding.** The
pre-investigation re-ran that exact probe (`p5/pf_cloud_probe.cpp`, same
seed/cloud/params) and found the 6 failing points are precisely the six
hydrostatic-axis points (`sigma = (p,p,p,0,0,0)` for `p` in
`{0,1,5,20,50,100}`) — not "random general (non-diagonal) points" as
`Ladruno_implementation/_adr94_components.md` finding 3 originally said. That
sentence is corrected in this PR.

`Closest_Point` continues to refuse the whole StiffSoil family (D3,
unchanged) — the fix is Backward_Euler-only. As a consequence, D6's
`ELASTICITY_STRESS_DERIVATIVE` for `StiffSoil_EL` is **not required**: the
block only matters under `Closest_Point`, which StiffSoil never reaches.

## 2. Measured

**Gate 4 (byte identity).** The 23-deck / 282-row `Backward_Euler` baseline
(`Ladruno_implementation/adr97_oracle/baselines/be_secant_baseline_3622d6214.json`)
is still byte-identical, *including* the four decks that build an explicit
integrator (`cube/vm/{Forward_Euler,Forward_Euler_Subincrement,
Modified_Euler_Error_Control,Runge_Kutta_45_Error_Control}/Continuum`), which
now pass `experimental_integrator 1` from the gate-4 dumper
(`dump_hist.py`) — D5 is a parse-time gate, not a behavior change, so the
committed history is untouched. 10/10 (fast slice + full sweep).

**D5 refusal coverage** (`tests/test_adr97_p6_failloud.py`, appended
section): for each of the four explicit methods, checked in **both**
directions — refused with the flag unset, refused with it explicitly `0`,
accepted with it `1`. A typo'd option name (`experimental_integratr`) is
still rejected (the ADR-94 wp/94a unknown-token contract). `Backward_Euler`
and `Closest_Point` construct identically across all three flag states
(unset / 0 / 1). One dedicated test re-confirms `tangent_type Algorithmic` is
still refused on an *opted-in* explicit integrator — this closes a coverage
gap D5 would otherwise open silently in the pre-existing
`test_algorithmic_is_refused_with_any_other_integrator`: three of its four
parametrized methods (`Forward_Euler`, `Modified_Euler_Error_Control`,
`Runge_Kutta_45_Error_Control`) would now be refused by D5 *first*, for an
unrelated reason, the exact "refused, but for the wrong reason" class of
defect the P3 report's `HB_sigma_ci` finding warned about. Two child-process
tests assert the refusal message text on the process's real OS-level
stdout/stderr — pytest's `capfd` cannot see this `.pyd`'s `cout`/`opserr`
output on this build (the ADR-94 `_run_child` idiom, reused here as a local
helper).

**StiffSoilShear/StiffSoilCap** (`tests/test_adr97_p5_stiffsoil.py`, new):

- Bug A regression: `MC_phi = 0.0`, one isotropic-compression step —
  `codes[0] == 0`, every committed stress finite.
- Bug B regression: a purely isotropic (hydrostatic) leg — every committed
  stress finite at every step, including the first (the exact hydrostatic
  trial that used to hit the 0/0). **Measured on this build:** the deeper
  path still cuts off with a clean refusal (`rc == -3`) a couple of steps in
  — a separate, legitimate Newton/step-size limit, not a NaN and not Bug B
  (confirmed identical in the pre-investigation's own post-fix probe,
  `p5/run_iso.log`). The assertion checks `codes[0] == 0` and that every code
  is in `{0, -3}` (ADR-94 B2's invariant: no committed state is ever
  non-finite, whatever the step codes are).
- Hardening path: a single-leg triaxial-compression ramp (16 steps, well past
  first yield) never commits `q` past the hyperbolic law's own asymptote
  `q_a = qf/Rf`, computed from the corrected `qf`. This invariant holds
  regardless of the accumulated `EpsQpShear` internal variable (hardening
  only slides `q` up *toward* `q_a`; it can never push `q` past it), which
  avoids needing to query the internal accumulator directly — a
  hand-computed upper bound on `eps_qp_shear` was rejected during drafting
  because it can under- or over-estimate depending on the sigma3 range (see
  §3).
- `StiffSoilCap` smoke test: `phi` in `{0, 30}` x cap IV in `{0, 100}`, all
  finite (StiffSoilCap shares neither Bug A nor Bug B — no defect found, no
  regression to guard against beyond "still builds and stays finite").
- `Closest_Point` refusal for both `StiffSoilShear` and `StiffSoilCap`,
  checked in both directions (constructible under `Backward_Euler` with the
  identical parameters, refused under `Closest_Point`).

**Full battery** (`test_adr84*`, `test_adr94*` except `test_adr94_matrix.py`
execution, `test_asdplastic_*`, `test_adr97*`): **230 passed, 2 skipped
(pre-existing)**.

`tests/test_adr94_matrix.py` was edited (`_int_opts()` now sends
`experimental_integrator 1` unconditionally — harmless for `Backward_Euler`,
required for the four explicit integrators it measures) but never *executed*
in this session; verified only by `pytest --collect-only -q
test_adr94_matrix.py` (1 test collected). See §3 for why.

## 3. Found while implementing

1. **`test_adr94_matrix.py`'s module-scoped fixture rewrites its own tracked
   doc on ANY execution, and a bare `--ignore=<filename>` inside a
   glob-expanded pytest invocation does not reliably exclude it.** Running
   `pytest test_adr94*.py --ignore=test_adr94_matrix.py -q` (the natural way
   to run "everything ADR-94 except the matrix file") still collected and
   *executed* `test_adr94_matrix.py`'s one test, silently overwriting
   `Ladruno_implementation/_adr94_matrix.md` — a hand-annotated, git-tracked
   doc with real attribution analysis — with the sweep's own regenerated,
   unannotated version. Caught by `git status` before committing;
   `git checkout -- Ladruno_implementation/_adr94_matrix.md` recovered it.
   The fix going forward: build the file list explicitly in the shell
   (`ls test_adr94*.py | grep -v matrix`) rather than trusting `--ignore`
   against a glob-expanded argument list. `pytest --collect-only` alone
   (never executing the file) is safe, since fixtures are lazy. Recorded in
   [[LEDGER_quirks]].
2. **A hand-computed upper bound on `EpsQpShear` is not a safe substitute
   for the real internal variable in an admissibility oracle.** The first
   draft of the hardening-path test (from the pre-investigation) used a
   fixed guess (`eps_qp_shear=0.05`) for the hyperbolic law's hardening
   term. A quick offline check showed the TRUE `eps_qp_failure` for this
   deck's parameters ranges roughly 0.029–0.058 depending on `sigma3` along
   the path — meaning the fixed guess is sometimes an over-estimate (which
   would silently WEAKEN the check by giving undue hardening credit) and
   sometimes an under-estimate (which would produce false failures). Neither
   direction is safe without actually querying the accumulated internal
   variable. Replaced with the asymptote invariant `q <= q_a = qf/Rf`, which
   is true independent of the hardening state and needs no internal-variable
   query.
3. **`drive()`'s multi-leg "owner" bookkeeping (copied from
   `test_adr97_p1_smooth.py`) requires each dof's target to change in at
   most one leg-to-leg transition** — legs are meant to progressively engage
   previously-idle dofs, not re-target one that is already moving. A
   two-leg "isotropic consolidation, then deepen the same axial strain"
   path is therefore not expressible with this helper (`AssertionError:
   path not component-disjoint`); a single proportional ramp to a
   non-isotropic target already mixes isotropic and deviatoric loading and
   was used instead.
4. **`_child_script`'s hardcoded `iv_type` string must match the
   registered specialization exactly, including the second concatenated
   hardening-law clause.** `IV_VM_LIN` in `test_adr97_p1_smooth.py` is
   `"BackStress(TensorLinearHardeningFunction):YieldStress(ScalarLinear
   HardeningFunction):"`, not just the `BackStress(...)` half; a first draft
   of the child-process test script used only the first half and every
   subprocess construction failed with "Material not found for input
   specification" — not a D5 refusal, a wrong-deck bug in the test itself,
   caught immediately by the child's own stdout.
5. **The `cot(phi)`/`0/0` bug class generalizes**: any `cot(x)`/`tan(x)`/
   `1/sin(x)` term that is later multiplied by `sin(x)` (or a factor
   containing it) in the same expression is a candidate for the exact same
   `Inf*0 -> NaN` failure at `x == 0`. Recorded in [[LEDGER_quirks]] as a
   general pattern to check for.

All five are in [[LEDGER_quirks]] (1, 4, 5 as new entries) or noted here
where they are test-authoring lessons rather than repo-wide quirks (2, 3).

## 4. Not done / left for the owner

- No mutation-testing gate was run for this WP (not requested in scope; D5
  is a parser-level gate and M9 is two narrowly-scoped arithmetic fixes, not
  new algorithmic machinery of the kind P1–P3's mutation gates targeted).
- `StiffSoilCap`'s `pc0 == 0` (uninitialised cap) gives an anomalous first
  response (wrong physics, not NaN) per the pre-investigation's
  `stiffsoilcap_drive.py` probe — out of scope here, recorded for a future
  pass.
- D6's `ELASTICITY_STRESS_DERIVATIVE` for `StiffSoil_EL`/`DuncanChang_EL`
  remains unimplemented; not required while `Closest_Point` refuses the
  whole family (D3).

## See also

[[97_ladruno_asdp_closest_point_adr]] · [[reviews/adr97_p3_report]] ·
[[LEDGER_vanilla_files]] · [[LEDGER_quirks]] ·
[[_adr94_components]] (the corrected finding 3)
