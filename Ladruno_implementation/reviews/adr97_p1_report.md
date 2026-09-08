---
title: ADR-97 P1 — Closest_Point + Algorithmic for the smooth families (report)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P1 (`wp/97b-cp-smooth`) — implementation report

**PR** [#819](https://github.com/nmorabowen/OpenSees/pull/819) (draft, based on
`wp/97a-plan-oracles`) · **plan** [[97_ladruno_asdp_closest_point_adr]] ·
**oracles** `Ladruno_implementation/adr97_oracle/`

## 1. What shipped

`integration_method Closest_Point` — the fully implicit closest-point return map
solved as ONE coupled Newton in `x = (σ_{n+1}, q_{n+1}, dλ)`:

```
R_sigma = s - s_tr + dl * E(s) * m(s, q)      (6 rows)
R_q     = q - q_n  - dl * h(s, q, m(s,q))     (n_q rows)
R_f     = f(s, q)                             (1 row)
```

and `tangent_type Algorithmic` — the exact consistent tangent of *that* map,
obtained by solving `J Z = [E; 0; 0]` with the ONE factorization already in hand
at convergence.

Scope: **VonMises** and **Drucker–Prager (flank + apex)** with **Null /
LinearScalar / LinearTensor / ArmstrongFrederick** hardening solved implicitly
at n+1 — **20 of the 46 registered specializations**. Everything else is refused
at parse time naming the ADR phase that will deliver it (MohrCoulomb + MCTC =
P2, HoekBrown = P3, StiffSoil = P5).

Six new interface members, each with an **inert base default plus an opt-in
trait**, so an unconverted family compiles unchanged and is refused rather than
silently approximated:

| block | macro | base default | trait |
|---|---|---|---|
| `df/dq` (uncontracted) | `YIELD_FUNCTION_IV_DERIVATIVE` | zero | `yf_has_cp_derivatives` |
| `dm/dσ` (6×6) | `PLASTIC_FLOW_STRESS_DERIVATIVE` | one-time-warned central difference | `pf_has_cp_derivatives` |
| `dm/dq` | `PLASTIC_FLOW_IV_DERIVATIVE` | zero | (same) |
| `dh/dq`, `dh/dm` | `HARDENING_FUNCTION_{IV,M}_DERIVATIVE` | zero | `hardening_policy_has_cp_derivatives` |
| `dE/dσ : m` | `ELASTICITY_STRESS_DERIVATIVE` | zero (D6) | `el_is_stress_dependent` |

`Backward_Euler` is untouched (D1). Build: `ec6091c4f`.

## 2. Gate-by-gate results

### Gate 1 — material-level correctness (`tests/test_adr97_p1_smooth.py`)

Committed stress against the P0 numpy CPPM oracles, 10 steps per leg, same
conventions and the same truncated `SQRT_2_over_3 = 0.816496580928`:

| family | path | rel. error |
|---|---|---|
| VM perfect | triaxial / simple shear / rotating normal | 7.4e-14 / 3.1e-13 / 2.4e-12 |
| VM `H = 7000` | triaxial / simple shear / rotating normal | 4.9e-13 / 5.9e-13 / 1.7e-12 |
| VM + AF (`ha=15000, cr=300`) | triaxial / simple shear / rotating normal | 1.3e-13 / 9.9e-13 / **9.2e-9** (α to 8.0e-9) |
| DP cone, associated (`η=η̄=0.4`) | compression | 4.9e-13 |
| DP cone, non-associated (`η̄=0.2`) | compression + shear | 5.1e-13 |
| DP **apex**, associated | hydrostatic tension | exact (σ = 50·**1**) |

The 9.2e-9 on the AF rotating-normal path is the **yield tolerance accumulated
over 20 steps**, not a defect: `f_absolute_tol` defaults to 1e-6 on a stress of
~45, i.e. 2e-8 relative per step. Every other path is at the 1e-13 floor.

* worst committed `|f|` over a whole path: **7.1e-15** (gate: ≤ tol);
* local Newton iterations, read from the new `cp_iterations` material response:
  **1** (perfect), **1** (linear hardening), **3** (Armstrong–Frederick) —
  gate ≤ 5, and exactly what the P0 oracle takes.

### Gate 2 — the tangent (same file)

`Algorithmic` against a **central difference of the binary's own assembled
internal force** on a free-DOF rig (`adr97_oracle/fd_tangent_driver.py`), no
numpy reference anywhere:

| rig | tangent | rel. error |
|---|---|---|
| 4 free DOFs (uniaxial strain) | `Backward_Euler` / `Continuum` — **negative control** | **0.573447** (pinned, reproduced to the digit) |
| 4 free DOFs | `Closest_Point` / `Algorithmic` | **1.51e-11** |
| 12 free DOFs (sheared) | `Closest_Point` / `Algorithmic` | **2.15e-10** |

ADR-94 two-cube heterogeneous model, per-step `ops.testIter()`:

| | step 1 | 2 | 3 | 4 | total |
|---|---|---|---|---|---|
| `Backward_Euler` / `Continuum` | 6 | 41 | 36 | 30 | 113 |
| `Closest_Point` / `Algorithmic` | 4 | 4 | 4 | 4 | **16 (7.1×)** |

The ADR's gate was written `testIter() <= 3`; **measured it is 4** — one
iteration to leave the previous step's converged tangent, two quadratic, one to
satisfy `NormDispIncr` at 1e-9. The test now asserts `<= 4` **and** at least a
2× reduction against `Continuum`, and says why in its body.

### Gate 3 — path error (same file)

Step-refinement error at N = 10 against each map's **own** N = 160 limit, on the
VM + Armstrong–Frederick rotating-normal path:

| map | error |
|---|---|
| `Closest_Point` | 5.68e-3 |
| `Backward_Euler` (cutting plane) | 2.31e-2 |
| **ratio** | **4.07×** (P0 oracle predicted 4.3×) |

With **linear** hardening the two maps agree to **1.45e-16** — which is exactly
why ADR-94 H6 could not see this defect: the shipped cutting plane is only path
dependent when `h` is not constant along the iterates, i.e. for
Armstrong–Frederick (22 of 46 specializations carry it).

The ADR states gate 3 as "step-halving reproduces the state to ≤ 1e-9"; the P0
oracle's 1.1e-12 is the per-step invariance to the **Newton start guess**, which
cannot be reached from Python (the start is not an option). The refinement-error
contrast above is the same statement made measurable from the binary.

### Gate 4 — `Backward_Euler` inertness (`tests/test_adr97_p4_inertness.py`)

* **23 decks / 282 committed-stress rows BYTE-IDENTICAL** against the
  pre-change binary `3622d6214`, each re-run in a FRESH SUBPROCESS
  (`Ladruno_implementation/adr97_oracle/baselines/`). Comparison is `==` on the
  doubles, not `allclose`. The whole file runs in **6.0 s** — a child costs
  0.2 s — so both the representative slice and the full sweep run on every push
  (it was written as `@pytest.mark.slow` on the assumption that 23 interpreter
  starts would cost minutes; measuring was cheaper than assuming).
* CP ≡ BE on **non-rotating-normal** perfectly plastic decks:
  0.0 (VM simple shear), 1.6e-33 (VM triaxial), 3.8e-17 (DP compression) —
  the two maps are the same point when the flow direction does not rotate.
* CP ≠ BE on Armstrong–Frederick: **2.3e-2**. An integrator that silently *was*
  `Backward_Euler` would pass the agreement half alone; this is the other half.

### Gate 5 — mutation

Dropping the `dl · dm/dσ` term of `Xi` from the Jacobian on a scratch build
turns `Algorithmic` into **exactly** `Continuum`: the free-DOF finite-difference
error goes 1.51e-11 → **0.573447** (4 DOF) and 2.15e-10 → **0.670133** (12 DOF),
both matching the P0 oracle's `Backward_Euler`/`Continuum` values to every digit
printed. 4 tests killed, 39 survivors — every survivor a case the mutation
provably cannot reach (the residual is untouched, so the committed stress stays
exact). Full record: [[reviews/adr97_p1_mutation]].

### Gate 6 — fail-loud (`tests/test_adr97_p6_failloud.py`)

13/13. Starved `n_max_iterations` refused on `LadrunoBrick` and **swallowed on
`stdBrick`** (pinned negative control, ADR-94 B2 — the reason every refusal test
here uses `LadrunoBrick`); `Algorithmic` refused with all four other
integrators and accepted with `Closest_Point`; `Closest_Point` refused for
MohrCoulomb; unknown `integration_method` / `tangent_type` tokens still
rejected; `strict_convergence` byte-inert on a converging deck and still
refusing a starved one; the ADR-94 B4 hydrostatic-tension reproducer commits a
finite, admissible apex stress.

### Gate 7 — portability

Zone-A dispatched on `wp/97b-cp-smooth`; the whole ASDP translation unit also
passes `g++ -std=c++17 -fsyntax-only` locally with the build's own include set,
which is how the GCC-only defects below were caught before the 20-minute build.

## 3. Two P0 findings honoured, not fixed

Both would change `Backward_Euler`, which **D1** forbids in this PR.

1. **`DruckerPrager_YF`'s cohesion IV is commented out of its own `f`** while
   `yf_hardening` still contributes `df/dk = -1` times its rate. `Closest_Point`'s
   `df/dq` must be the TRUE derivative of the `f` it solves, so a Drucker–Prager
   with cohesion hardening is **perfectly plastic under `Closest_Point`** and
   hardening under `Backward_Euler`. Pinned by
   `test_gate1_dp_cohesion_hardening_is_perfectly_plastic_under_cp`, which asserts
   BOTH halves — it turns red the day the yield function is fixed, which is the
   point.
2. **Armstrong–Frederick's saturation branch is kept.** The plan suggested
   `Closest_Point` drop it (it is non-differentiable, and the implicit AF update
   is already a contraction toward `‖α‖ = ha/cr`). It is kept because `f` is
   SHARED with `Backward_Euler`; `Closest_Point` differentiates the branch `f`
   actually takes, which keeps its Jacobian exactly consistent with its own
   residual, and the P0 oracle mirrors the same branch.

## 4. Found while implementing

* **The elastic-metric apex test degenerates on the hydrostatic axis.** The test
  is a sign test on `dot(dev_ret, dev_tr)` over the linearised cone step; with
  `dev_tr == 0` it reads `0 < 0` — CONE — and the cone Newton then has no flow
  direction at all. That state is ADR-94 B4's own hydrostatic-tension reproducer.
  Fixed with a `‖dev σ_tr‖ <= tol_f` short circuit (in the YIELD tolerance, so it
  stays unit consistent per ADR-94 M5).
* **`dm_dsigma_buffer` needs `this->`.** It lives in the dependent base
  `PlasticFlowBase<T>`; unqualified it is ill-formed on GCC and accepted by MSVC —
  the trap that cost the ADR-94 wave a red CI, caught here by the `-fsyntax-only`
  pre-flight instead.
* **`utuple_concat_unique_type` de-duplicates internal variables by TYPE.** A
  `VonMises_YF<BackStress<TensorLinear>, …>` paired with
  `VonMises_PF<BackStress<NullTensor>>` gets **two** back stresses that evolve
  independently. This is why `df_dq` / `dm_dq` select on
  `std::is_same<IVType, AlphaHardeningType>` inside each functor: each
  differentiates with respect to the variable IT reads and returns zero for the
  other, which is the correct derivative.
* **A `Path` time series returns zero at its last time point.** A prescribed
  strain driver ending exactly on the final `-time` entry unloads the whole path
  in one step, and a plasticity material then reports a perfectly plausible
  ON-SURFACE stress that is the wrong point on the surface.

All four are in [[LEDGER_quirks]].

## 5. Mutation gate

See [[reviews/adr97_p1_mutation]].

## 6. Open

* `Numerical_Algorithmic_*` still differentiate `compute_local_stress()`, a third
  map nothing commits (ADR-97 **P4**).
* `Closest_Point` uses `Eigen::FullPivLU` on the (14×14 for the P1 families)
  Newton system — chosen for its singularity detection, not for speed. P6's
  measurement decides whether a partial-pivot factorization with an explicit
  conditioning guard is worth it.
* The apex consistent tangent is exact only while `d σ_apex/dq == 0` (true for
  every yield function this ADR ships, and checked at run time by a finite
  difference); a yield function whose apex moves with an internal variable gets a
  one-time warning and the same zero. P5.
* `Closest_Point` is **not** the default and `Algorithmic` is **not** the default
  tangent. The flip is a separate PR gated on P6 (D1).
