---
title: ADR 94 — R1-A verdicts (numerics lane: H1, H6, H7, H8)
project: Ladruno
status: draft
owner: nmora
tags: [implementation, material, review]
---

# ADR 94 R1-A — numerics lane verdicts

Measured on build `52314165a14b2a8cd6f846347b19c5bdd73528e7` (`ops.ladrunoBuild()`
verified against HEAD). Tests: `tests/test_adr94_hlist_numerics.py`, **11 passed
in 0.25 s** (order-independent; 26 passed together with
`test_adr84_p2a_strict_convergence.py` + `test_asdplastic_mctc.py`). Oracle:
`Ladruno_implementation/adr94_oracle/hex8_tangent.py` (validated against the
OpenSees `stdBrick` elastic tangent to **4.4e-16**).

No C++ was edited. Line numbers are `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`.

---

## H1 — class-static `Stiffness` — **CONFIRMED · blocker**

**The number that decides it.** In a two-element model whose elements are in
genuinely different states (their stand-alone tangents differ by **13.56 %**),
the two assembled diagonal blocks are **bit-identical (0.000e+00)**, and the
plastic element's block is **13.19 %** away from its own correct tangent while
the elastic element's block reproduces its own to **4.0e-11**.

**Mechanism.** `Stiffness` (and `dsigma`, `depsilon_elpl`,
`intersection_stress/strain`) are `static` members of
`ASDPlasticMaterial3D<E,Y,P,tag>` (declared 4116–4120, defined 4150+). `tag` is
the *template* class tag, so the sharing spans every Gauss point, every element
**and every material tag** of one YF×PF×EL combination. The integrators write
`Stiffness` (2088, 2217–2233, 2303, 2354 → `ComputeTangentStiffness`, and 1370,
1425, 1601, 2407, 2560, 2669, 3090, 3437); `getTangent()` (698) only copies it.
Hosts `setTrialStrain` every GP in `update()` and call `getTangent()` in a later
loop (`Brick.cpp` 1069 vs 1201) ⇒ **the whole model is assembled with the
tangent of the last Gauss point integrated.** Swapping which element is plastic
swaps which block is wrong — the winner is decided by domain iteration order,
never by the element.

**Blast radius.** All 46 registered specializations; all hosts that separate
`update()` from `getTangent()` — `Brick`, `LadrunoBrick`, `TenNodeTetrahedron`,
`BbarBrick`, `SSPbrick`, `BrickUP`, the Bezier solids. Invisible in every
existing test because all drivers are single-element homogeneous strain.
Also the permanent ADR-75b blocker: a threaded element loop cannot be
deterministic while a class-static holds per-GP state.

**Pinned by** `test_H1_one_static_tangent_is_shared_by_every_element`,
`test_H1_shared_tangent_follows_the_last_element_integrated`
(+ `test_H1_oracle_matches_opensees_elastic_stiffness` validating the oracle).

---

## H6 — cutting plane, not closest point — **PARTIAL** (accuracy half REFUTED on VM; tangent half CONFIRMED · major)

**Classification: CONFIRMED by reading.** The loop at 2244–2320 re-evaluates
`n`, `m`, `H` at the current iterate and corrects incrementally
(`σ ← σ − δλ·E·m`, IVs advanced with `h` at the updated stress) — Ortiz–Simo
cutting plane, not a closest-point projection.

**Accuracy claim REFUTED on VonMises.** For `f = ‖dev s‖ − √(2/3)σ_y` with
associated flow and linear isotropic hardening the return direction does not
rotate, `Φ(λ)` is exactly linear, and the cutting plane lands on the
closest-point answer in one Newton step. Measured against the numpy closed form
at `Δλ = 2.3e-3`: **relative stress error 3.7e-14**. ⇒ **the H6 row's proposed
"convergence-order study on VM + linear hardening" is degenerate** and cannot
distinguish the two maps. A rotating-normal YF (Lode-dependent MC, kinematic
hardening, stress-dependent elasticity) is required; that is R4's matrix, not
this lane's.

**Tangent claim CONFIRMED, and it is the operational finding.** On a
homogeneous plastic state (so H1 is invisible), *no* `tangent_type` reproduces
the consistent tangent of the committed map. Measured vs
`hex8_K(vm_consistent_tangent)` at `Δε_zz = −3.6e-3`:

| `tangent_type` | error vs consistent | global Newton iters |
|---|---|---|
| `Continuum` | **57.3 %** | 3 |
| `Secant` (**the shipped default**) | 79.9 % | 16 |
| `Elastic` | 102.5 % | 23 |
| `Numerical_Algorithmic_FirstOrder` | 31.0 % | 4 |
| `Numerical_Algorithmic_SecondOrder` | 31.0 % | 4 |

The numerical pair are closest yet still 31 % out because they differentiate
`compute_local_stress()` — a third map that is never the one committed (ADR-84
P4). The `Continuum` error is first order in the step (57.3 / 46.2 / 31.3 /
19.1 / 10.7 % at 1/2/4/8/16 steps) — it is the `Δλ→0` limit, so it is only ever
right in the limit the step vanishes.

**Actionable:** the DEFAULT `Secant` costs **5.3×** the global iterations of
`Continuum` on this step and is the *worst* of the two cheap options. Changing
the default is upstream-facing (ADR-84 §3), but documenting `Continuum` as the
recommended setting costs nothing.

**Blast radius.** Every specialization, every host. **Severity major**
(convergence cost + the impossibility of a consistent tangent), *not* a wrong
answer: stresses are exact.

**Pinned by** `test_H6_backward_euler_is_exact_for_von_mises` (REFUTED),
`test_H6_no_tangent_option_reproduces_the_consistent_tangent`,
`test_H6_continuum_tangent_error_is_first_order_in_the_step`.

**Owner decision D2:** recommend **document + rename, no rewrite in this ADR** —
but the review should record that a closest-point BE would also buy the
consistent tangent that no option currently offers.

---

## H7 — the `dLambda + deltaLambda < 0` fallback — **CONFIRMED · major** (and worse than the row states)

**The number that decides it.** With `H_iso = −120000` (softening steeper than
`n:E:m = 2G = 53846`) the branch at 2298 fires on **iteration 0**, so nothing
has been corrected: the commit is the **elastic predictor to 0.000e+00**, the
step returns **0 (success)**, and `f_VM` at the four commits is **+3.50, +31.49,
+59.49, +87.48 kPa** against `f_absolute_tol = 1e-6`. The material behaves
purely elastically for the whole run while drifting monotonically further
outside the surface, and the global Newton converges every step.

**Row correction.** The H7 row predicts "a partially corrected iterate". For the
reachable case it is the *elastic predictor exactly* — which is worse, because
it is a plausible-looking state. (A partial iterate is possible only if the
branch fires at `iter > 0`, i.e. on a negative overshoot larger than the
accumulated `dLambda`; not observed here.)

**New finding not in the plan.** `strict_convergence 1` does **not** gate this
branch — it returns 0 before the `be_strict` check at 2338. So this is an
**eighth** silent-accept site, inside the DEFAULT integrator, in addition to the
seven H5 lists. Measured: flag-on is bit-identical to flag-off, all four codes
`0`, `f_VM = +87.5`.

**Blast radius.** Any specialization that can reach `H > n:E:m` — every YF with
a negative `ScalarLinearHardeningParameter` (VonMises, DruckerPrager,
MohrCoulomb, the StiffSoil pair), and any softening/strain-localizing model
added later. Hosts: all — the material reports success, so no host can catch it.

**H7 second half (`depsilon_elpl` provenance) — CONFIRMED by reading ·
doc-only today.** `depsilon_elpl` is assigned only at 516
(`compute_local_stress`), 1430/1441 (`Forward_Euler`) and 1606–1623
(`Forward_Euler_Subincrement`). `Backward_Euler` **never** assigns it, yet
`ComputeTangentStiffness()` (458/460, 472/474) feeds it to `pf(...)` and
`yf.hardening(...)` for `Continuum` **and `Secant` (the default)**. It is
therefore zero, or stale from another integrator / another GP's numerical
tangent. **Blast radius today is zero**: no registered PF or YF reads its
`depsilon` argument (verified by grep over `PlasticFlowDirections/*.h` and
`YieldFunctions/*.h`), so it is unobservable from Python. It becomes a live
wrong-tangent bug the moment a stress-dilatancy PF lands — which is exactly what
ADR-92/93 (SANISAND) points at. Fix it with H1, since both are the same
class-static defect.

**Pinned by** `test_H7_inconsistency_branch_commits_the_elastic_predictor`,
`test_H7_strict_convergence_does_not_gate_the_inconsistency_branch`.

---

## H8 — `Backward_Euler_LineSearch` — **CONFIRMED · major**

**The number that decides it.** On the ADR-84 MC tet leg, `Backward_Euler`
completes **20/20** steps; `Backward_Euler_LineSearch` completes **2/20**. The
integrator whose name promises robustness is strictly less robust than the plain
one.

Three independently measured defects:

1. **`n_max_iterations` is inert.** 2385 hardcodes `max_iter = 30` (and
   `tol_rel = 1e-8`). Measured: `n_max_iterations` 2 and 100 give bit-identical
   histories. On the same rig `Backward_Euler` at `n_max_iterations 2` is the
   ADR-84 exhaustion reproducer (worst `f_MC` **77.6** vs **6.3e-4** at 100), so
   the option is not inert in general — only here.
2. **`strict_convergence` never reaches it.** `be_strict` appears only at 2081,
   2086, 2184, 2338, all inside `Backward_Euler`. Measured: flag on/off
   bit-identical. The fork's one loud-failure switch is a no-op the moment a
   user selects this integrator.
3. **The "line search" cannot cut α, and the "substepping" truncates the
   step.** The acceptance test (2477–2486) is on a *linear prediction*
   `Φ + (dΦ/dλ)·dl`; for an unclipped Newton direction `dl = −α·Φ/(dΦ/dλ)` it
   reduces to `|1−α| ≤ 1 − 1e-4·α`, true for every `α ∈ (0,1]` ⇒ `α = 1` is
   accepted on the first try, always. The split loop (2578–2596) does **not**
   chain substeps: on failure it halves `dEps` and solves **one** reduced
   increment, overwriting `TrialStrain` with `CommitStrain + dEps/2^k` and
   returning **success for a strain the element never asked for**. (Confirmed by
   reading; not reproduced as a *successful* truncation — on the MC rig the
   inner Newton failed at all six split levels and returned −1.)
   `Eelastic = et(CommitStress)` is also hoisted outside the split loop (2379),
   so a stress-dependent `EL` is frozen at commit for every substep (H15's
   concern, same site).

**Blast radius.** Only decks that select `integration_method
Backward_Euler_LineSearch` — but such a deck silently loses `n_max_iterations`
*and* `strict_convergence`, which is precisely the H13 "silent
misconfiguration" failure mode. Recommend either fixing all three or refusing
the option.

**Pinned by** `test_H8_line_search_ignores_n_max_iterations`,
`test_H8_line_search_ignores_strict_convergence`,
`test_H8_line_search_is_less_robust_than_plain_backward_euler`.

---

## Method notes for the next lane

* `printA -ret` (dense) returns **empty** for every SOE but `FullGeneral` —
  `getA()` is null elsewhere (`OpenSeesCommands.cpp` 2718). Use
  `printA('-sparse','-ret')`, which works with `UmfPack` and returns
  `{rowIndices, colIndices, values}`. `printA` calls `formTangent()` itself, so
  it reads the tangent **after** the last `update()` — exactly the H1 window.
* A homogeneous single-element rig (x,y fixed everywhere, z fixed at the base,
  z **loaded** at the top) has 4 free DOFs, constant strain, and therefore
  measures the material tangent with H1 switched off. Use it for anything that
  is not H1.
* `Brick::setResponse("strains")` reports the **material's** `getStrain()`
  (line 1880), not `B·u` — so a strain-truncating integrator (H8 defect 3) is
  directly observable by comparing it to the imposed field.
* `cout` from the `.pyd` does not survive pytest's `capfd` reliably; assert on
  physics, not on the 161 `cout` lines (H12).

## See also

[[94_asdplastic_review_plan]] · [[84_ladruno_mc_tension_cutoff_adr]] ·
[[40_ladruno_performance_adr]] (ADR-75b threading, blocked by H1) ·
`tests/test_adr94_hlist_numerics.py` · `Ladruno_implementation/adr94_oracle/hex8_tangent.py`
