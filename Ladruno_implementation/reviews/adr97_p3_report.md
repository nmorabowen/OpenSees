---
title: ADR-97 P3 — closest point on the CURVED Hoek–Brown surface (report)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P3 (`wp/97d-cp-hoekbrown`) — implementation report

**PR** [#825](https://github.com/nmorabowen/OpenSees/pull/825) (draft, based on
`wp/97c-cp-principal` = P2 [#824](https://github.com/nmorabowen/OpenSees/pull/824)
→ P1 [#819](https://github.com/nmorabowen/OpenSees/pull/819) →
[#817](https://github.com/nmorabowen/OpenSees/pull/817)) · **plan**
[[97_ladruno_asdp_closest_point_adr]] · **P2** [[reviews/adr97_p2_report]]
· **oracle** `Ladruno_implementation/adr97_oracle/cppm_hb.py`

## 1. What shipped

`integration_method Closest_Point` + `tangent_type Algorithmic` for the
**Hoek–Brown family**, as a principal-stress-space return to a **curved**
surface (Clausen & Damkilde 2008).

P2's principal machinery is *extended*, not duplicated: the spectral
decomposition, the eigenprojection back-transform (rotation term
`(y_i−y_j)/(x_i−x_j)`, its l'Hôpital limit on a degenerate trial eigenvalue,
`Rs^-1` built as the Voigt image of `Q^T · Q`), `cp_apply_tangent_policy`, the
`d eps^p = E^-1 (sigma_tr − sigma_ret)` convention and the strict-convergence
contract are all P2's, verbatim. What is new is the return itself, because the
Hoek–Brown meridian is a curve, not a plane.

| piece | how |
|---|---|
| the Newton unknown | `w`, with `arg = (w²)^(1/a)` ⇒ `y1 = T − (σci/mb)(w²)^(1/a)` and yield row `y_i − y_j − σci w²`. This is a **requirement**, measured by the P0 oracle: with `y1` as the unknown the map diverges near the apex (the surface exists only for `arg ≥ 0`; the first step from the elastic predictor overshoots to `arg = −2.2e-03` and the next Jacobian is singular). `(w*w)^(1/a)` rather than `w^(2/a)`, because the latter is NaN for a negative iterate; the two agree for `w > 0`, the residual is then even in `w`, and a Newton that crosses zero converges to `−w*` and returns the same stress. |
| flow direction | **normalized** in the residual, the multiplier rescaled by `|m|`. `|m| ~ arg^(a−1)` blows up exactly where near-apex returns land (460.6 there against 3.9 on an ordinary face point); the oracle measured the worst-case Newton count at 6 un-normalized and 5 normalized. |
| the potential | a **frame-consistent** Hoek–Brown potential `g_ij = y_i − y_j − σci (s − mb_psi·y_i/σci)^a`, with analytic principal-space derivatives, in `Closest_Point`'s own path. `dm/dy` has exactly one non-zero entry, at `(i,i)` — that entry **is** the curvature term. See §3. |
| region selection | apex by an **exact elastic-metric cone test** — the domain's tangent cone at the vertex is the negative octant (the meridian meets the hydrostatic axis vertically), so the apex region is exactly `apex + cone{D3·m_ij(apex)}`, decided by the dual/facet form of that cone; then the face; then the two edges ordered by the **face return's own margins** `y1−y2`, `y2−y3`, which for a curved surface *are* the exact signed boundary functions (P2's boundary planes need a piecewise-linear surface and do not generalise). |
| the tension regime | **no Rankine face return exists, and none is missing.** `f_shear ≤ 0` already forces `y1 ≤ T`, so the composite's zero set is the shear surface plus its own natural vertex; the oracle sampled 4000 surface points and `f_tension` won on **none** of them (`min(f_shear − f_tension) = 1.004e-04`). Outside the domain the tension branch wins exactly on `{y1 > T and y3 > T}`, which **is** the header's `CHECK_APEX_REGION` set. A Rankine return onto `y1 = T` from the oracle's `[645,265,255]` trial leaves `f_shear = +123.34` kPa (2.31 % of the strength scale) — inadmissible. **The shear/tension corner IS the apex.** |
| tangent | per region from the converged Jacobian (`dy/dx = (dy/dz)·J^-1·[I;0]`), so the `dm/dy` curvature term is carried exactly; **rank 0 at the apex**, as for the DP (P1) and MC (P2) apices. |
| tolerance | admissibility carries a **gradient-scaled floor**, `64·eps·max(\|y1\|, scale)·(1 + a·mb·arg^(a−1))`. `\|df/dy1\|` diverges at the apex, so an absolute 1e-10 gate is unattainable within ~1e-2 kPa of the vertex (the oracle's own max `\|f\|` there is 1.08e-10 against a floor of 1.16e-08). The Hoek–Brown instance of ADR-94's `f_relative_tol` lesson. |

`Backward_Euler` and every YF/PF path it executes are untouched (ADR-97 D1).

**Build.** `d0de2d7c8` is the last commit that changes anything under `SRC/`;
everything after it is tests and docs. The gates were first measured on that
binary and re-measured, identically, on the post-mutation-revert rebuild, whose
`ops.ladrunoBuild()` stamps `6e937fabb6bd4183a994bb8227ff69c328084992` —
`git diff d0de2d7c8 6e937fabb -- SRC` is empty. Note that `ladrunoBuild()`
alone could NOT have told the mutated binary from the restored one (the mutation
was never committed, so both stamp the same hash); the restore was verified by
BEHAVIOUR — 32/32 on the restored build against 11 failures on the mutated one.

## 2. Support count — 23 of 46

`HoekBrown_YF` appears in **7** registered specializations and `HoekBrown_PF` in
**6**, and in exactly **one** of them are both of the Hoek–Brown family:
`HoekBrown_YF<BackStress<Null>> × HoekBrown_PF<BackStress<Null>>`.

**20 (P1) + 2 (P2) + 1 (P3) = 23 of 46.** The other **eleven** Hoek–Brown
pairings stay refused at parse time under the same family-marker equality rule
P2 introduced (`ladruno_cp_principal_family` non-zero only when the YF's and the
PF's markers are equal and non-zero *and* every internal variable is inert). All
eleven are pinned by `static_assert` in the g++ pre-flight; the seven that are
reachable from a deck are also pinned at run time, **in both directions**.

## 3. The potential decision — implemented in P3's path, recorded for the owner

`HoekBrown_PF::g` is evaluated in the **un-negated (tension-positive) frame**
while `HoekBrown_YF` negates first, and then destructures the ascending tuple as
`[sigma3, sigma2, sigma1]` and feeds the tree's most **compressive** principal
into `arg = mb_psi·sigma3/sigma_ci + s`. On any compressive state that `arg` is
negative, so `g` always takes its `else` branch and collapses to a **Tresca**
potential. The P0 oracle quantified all four consequences:

* the header's own central-difference `dg/dσ` at `[-2000,-6000,-25000]` is
  `[1,0,-1,0,0,0]` for `mb_psi = mb`, `mb/2` **and** `0` — `HB_mb_psi` has **no
  effect** and the flow is exactly non-dilatant (trace 0);
* it is **32.86°** off the frame-consistent normal, at `mb_psi == mb`, i.e.
  exactly where the deck is asking for *associated* flow;
* at the apex all six of its flow directions have negative trace
  (`trace(D3·m) = −7.78e+08`), so the hydrostatic direction is not in its return
  cone and a trial pushed past the tensile corner has **no return to the apex at
  all** — the mechanism behind the residual recorded in
  `tests/test_adr94_hlist_hb.py`;
* same trial, same surface: `y = [-3642.77, -6441.62, -25123.72]` (plastic
  volumetric strain +2.208e-05) under the intended potential, against
  `y = [-3278.27, -6000.00, -23721.73]` (−6.8e-21, zero dilation) under the
  header's.

Fixing `g` changes `Backward_Euler`, which **D1** keeps byte-identical, so it is
the **owner's separate PR** (recorded in [[LEDGER_quirks]] and the ADR
implementation log). `Closest_Point` uses the frame-consistent potential in its
own code path and the difference is **pinned by a test in both directions**
(gate 4 below), which turns red the day `g` is fixed — that is its purpose. When
`HB_mb_psi == HB_mb` the P3 potential **is** the yield function's own gradient,
so associated flow is exact under `Closest_Point`.

## 4. Gate-by-gate results

### Gate 1a — the regions, against the P0 oracle (`tests/test_adr97_p3_hoekbrown.py`)

One step of prescribed total strain `eps = E^-1 sigma_tr`, so the elastic
predictor **is** the oracle's trial state. `σci = 50000 kPa, mi = 10, GSI = 60,
D = 0, E = 5e7 kPa, ν = 0.25, HB_mb_psi = mb` — the same material as
`tests/test_adr94_hlist_hb.py`.

| trial | region | rel. error vs the oracle | committed `\|f\|` | `cp_iterations` |
|---|---|---|---|---|
| face, no shear | face | **1.448e-15** | 1.09e-11 | 4 |
| face, sheared (rotation term live) | face | **1.735e-15** | 1.82e-11 | 4 |
| edge `y1 == y2` | line1 | **7.259e-16** | 3.64e-12 | 3 |
| edge `y2 == y3` | line2 | **1.268e-15** | 3.64e-12 | 3 |
| apex, hydrostatic tension | apex | **0.0** | 0.0 | 1 |
| apex, slightly deviatoric | apex | **0.0** | 0.0 | 1 |
| `[645,265,255]` — the header says **apex** | **face** | **1.851e-13** | 1.73e-11 | 5 |

Pinned at 1e-10 relative, measured at 1e-13 or better; the same justification as
P2's tighter-than-1e-6 pin (the only platform-variable step is a 3×3 symmetric
eigen-decomposition, and every pinned trial is either non-degenerate or returns
a state invariant under the eigenvector ambiguity).

The last row is the one that matters for the header: `CHECK_APEX_REGION` is the
**Euclidean** octant `x_tr ≥ T`, and `D3·(octant)` is a strict subset of the
octant, so it always **over**-claims. On this trial `APEX_STRESS` would commit
`T·[1,1,1]` instead of the correct face return — 122.6 kPa, **2.29 % of the
strength scale**, of silently lost strength (admissible, but wrong).
`Closest_Point` classifies where `E` is in scope and never calls it.

**Curved-face Newton counts: 4, 4, 5 — the ADR's `≤ 5` gate, met.** That number
is the whole reason for the natural variable and the normalized flow direction.

### Gate 1b — admissible at every commit

| path | steps | plastic | worst committed `f` | `cp_iterations` seen |
|---|---|---|---|---|
| triaxial compression | 10 | 8 | **0.0** | {0, 3} |
| simple shear | 10 | 10 | **+7.3e-12** | {4} |
| rotating principal directions (2 legs) | 20 | 18 | **+7.3e-11** | {0, 3, 4} |

Gate `1e-8 · strength_scale = 5.35e-05`. The rotating path is checked to
actually rotate: the first eigenvector moves **16.80°** between the legs.

### Gate 1c — the ADR-94 tension plateau, sharpened

Uniaxial **stress** tension (lateral faces free), so the stress path is exactly
`σ = (σxx, 0, 0)` and the limit is the root of
`σxx = σci·(s − mb·σxx/σci)^a`.

| quantity | value |
|---|---|
| exact closed-form limit | **244.4854419069** kPa |
| `Closest_Point` plateau | **244.4854419** kPa — **1.16e-16 relative** |
| the apex `T` (ADR-94's own pin) | 245.0151819 kPa |
| steps committed | 24 of 60, then `analyze -> -3` |

ADR-94's test pins BE's plateau at `T` with `rel = 5e-2`, which does not
separate the two: `T` is **0.216 % above** the exact limit, because at `y1 = T`
the clamp leaves `f_shear = y1 > 0` and the surface is reached strictly *before*
the vertex. `Closest_Point` lands on the exact root. Both integrators still stop
at the limit state — correctly: on that edge the consistent tangent is rank 1 in
principal space and the free lateral DOFs have no stiffness left.

And a trial pushed **2× / 20× past the tensile corner** now commits the finite
apex `T·[1,1,1]` (`f` = 0.0 / 5.7e-14) rather than stalling, because the
frame-consistent potential has a return to the vertex where the shipped Tresca
`g` has none.

### Gate 2 — the consistent tangent

`Algorithmic` against a **central difference of the binary's own assembled
internal force**. No numpy reference anywhere.

| rig | region | rel. error |
|---|---|---|
| free-node, face axis aligned (separation 1.33e-02) | face — curvature **and** rotation terms live | **3.342e-09** |
| free-node, face sheared (separation 1.20e-01) | face | **3.869e-09** |
| `fd_tangent_driver` oedometric, load `-3.0e4` | edge, **degenerate eigenvalue** (l'Hôpital branch) | **1.485e-10** (`rel_fro` 1.915e-10) |

The **apex tangent is rank 0** by construction (the vertex does not move), so
there is no finite-difference gate on it and none is missing — the same
statement P1 makes for the Drucker–Prager apex and P2 for the Mohr–Coulomb
vertex. The free-node rig (P2's) is still the only one that reaches a genuine
face state with three separated principal stresses.

**Iteration contrast — an outcome, not a ratio.** On the Hoek–Brown oedometric
deck (6 steps, `-3.0e4`/node):

| deck | result |
|---|---|
| `Closest_Point` / `Algorithmic` | 2, 2, 2, 2, 4, 4 = **16**, `analyze -> 0` |
| `Backward_Euler` / `Continuum` | **did not converge** (`analyze -> -3`) |
| `Backward_Euler` / `Secant` | **did not converge** (`analyze -> -3`) |

### Gate 3 — hardening: there is none to gate

`HoekBrown_YF` and `HoekBrown_PF` are registered with
`BackStress<NullHardeningTensorFunction>` and nothing else, and
`YIELD_FUNCTION_HARDENING` returns 0.0. The principal map has no `q`-row, and
the all-IVs-inert requirement is folded into `ladruno_cp_principal_family` at
**compile** time, so a hypothetical `HoekBrown_YF<ArmstrongFrederick..>` is
refused at parse time rather than integrated with its hardening ignored. The
test asserts the committed back stress is exactly zero after 10 plastic steps.

### Gate 4 — `Backward_Euler` inertness, and the potential gap pinned

* `tests/test_adr97_p4_inertness.py` re-run on this binary: **23 decks / 282
  committed-stress rows byte-identical** against the pre-ADR-97 baseline
  `3622d6214`, each in a fresh subprocess. The Hoek–Brown decks are among them,
  so the two HB headers P3 edits are covered directly.
* **The potential gap, pinned in both directions:**

  | | principal return | plastic `eps_vol` |
  |---|---|---|
  | `Closest_Point` | `[-3642.773683, -6441.622240, -25123.715278]` | **+2.208111e-05** |
  | `Backward_Euler` | `[-3278.269706, -6000.000004, -23721.730311]` | **+2.09e-13** |

  Relative stress gap **5.910e-02**, i.e. **26.20 % of the strength scale**.
  CP's plastic volumetric strain reproduces the oracle's `+2.208e-05` for the
  frame-consistent potential; BE's is the oracle's `-6.8e-21` (zero dilation) up
  to the cutting plane's own residue. The test asserts BOTH halves, so it turns
  red the day the shipped `g` is fixed.

* **Step refinement replaces P2's exact step independence.** A closest point
  onto a *plane* along a proportional path is step independent (P2: 1.5e-16 to
  5.2e-16 over N = 1…40). A *curved* surface cannot be: after step k the state
  sits on the surface, so step k+1's trial starts from a curved point.

  | N | rel. vs the oracle (the one-step map) | rel. vs N = 160 |
  |---|---|---|
  | 1 | **1.448e-15** | 2.594e-03 |
  | 4 | 1.184e-03 | 1.413e-03 |
  | 10 | 2.011e-03 | 5.885e-04 |
  | 40 | 2.474e-03 | 1.260e-04 |
  | 160 | 2.601e-03 | 0 |

  N = 1 **is** the oracle, and the refinement sequence converges monotonically
  to an incremental limit **0.26 %** away. The test pins both, because a silent
  fall-back to an incremental cutting plane would look like a much larger
  version of exactly this.

### Gate 6 — fail loud

* **All seven runtime-reachable mixed Hoek–Brown pairings refused**, each
  checked in BOTH directions — it must *construct* under `Backward_Euler` and be
  *refused* under `Closest_Point`:
  `MohrCoulomb_YF × HoekBrown_PF`, `HoekBrown_YF × VonMises_PF` (Null and
  TensorLinear), `HoekBrown_YF × MohrCoulomb_PF`,
  `HoekBrown_YF × DruckerPrager_PF`, `VonMises_YF × HoekBrown_PF`,
  `DruckerPrager_YF × HoekBrown_PF`. The remaining four (ArmstrongFrederick
  variants) are pinned by `static_assert` in the pre-flight.
* the matched `HoekBrown_YF × HoekBrown_PF` pairing is **accepted** (the
  positive half — without it the refusals above pass for a build in which
  `Closest_Point` is refused for everything);
* `tangent_type Algorithmic` still refused with `Backward_Euler` (D2 survives P3);
* `strict_convergence 1` **byte-inert** on a converging HB deck (gap exactly 0.0);
* exactly hydrostatic tension, hydrostatic + 1e-15 and + 1e-9 deviators all
  commit the finite admissible vertex — the class of state ADR-94 B4 used to
  commit NaN through — and hydrostatic **compression** stays elastic
  (`f = -2.695e+03`).

**Full battery:** `test_adr84*`, `test_adr94*` (except `test_adr94_matrix.py`),
`test_asdplastic_*`, `test_adr97*` — **204 passed, 2 skipped, 0 failed**. Both
skips are the pre-existing ADR-94 R3a time-box skips.

## 5. Two vacuous tests found and fixed

`HB_sigma_ci` is **not a parameter** — the name is `HB_sigci`. Both
`tests/test_adr97_p6_failloud.py::_mat_hb` (P1) and the Hoek–Brown rows of
`tests/test_adr97_p2_principal.py` (P2) used the wrong spelling, so those decks
were rejected for a **missing parameter** and the refusal assertions never
exercised the family gate at all. Both are corrected here, and the assertions
P3 inverts (`HoekBrown_YF × HoekBrown_PF` is now accepted) are inverted in place
with the reason recorded — the same treatment P2 gave P1's MohrCoulomb
assertion.

This is why every gate-6 row in the P3 file is checked in **both** directions: a
refusal test built on a wrong parameter list refuses for the wrong reason and
looks exactly like a passing gate.

## 6. Found while implementing

1. **`HB_sigma_ci` vs `HB_sigci`** — §5.
2. **A `template class` explicit instantiation is still the wrong pre-flight
   shape** (P2's finding 4), and `#define private public` now needs care of its
   own: GCC 15's libstdc++ declares `basic_stringbuf::__xfer_bufptrs` private
   and re-declares it later, so flipping the keyword before `<sstream>` is a
   hard `-Wtemplate-body` error. Include the standard library and Eigen *first*,
   then flip the keyword for the project's own headers.
3. **`AllASDInternalVariableTypes.h` and `AllASDHardeningFunctions.h` have no
   include guard.** Including either directly in a translation unit that also
   includes `ASDPlasticMaterial3D.h` is a redefinition storm.
4. **A numpy transcription of the C++ is worth writing before the build.** The
   region layout, the `Ydz` maps and the analytic Jacobian were validated
   against the oracle (all seven regions ≤ 1.9e-13, same iteration counts) and
   against a finite difference of the transcription's own return, in seconds,
   before spending twenty minutes on a build.

All of 1–3 are in [[LEDGER_quirks]].

## 7. Mutation gate

See [[reviews/adr97_p3_mutation]].

## 8. Open

* **`HoekBrown_PF::g` is still wrong** (§3). The owner's separate PR; this WP
  supplies the warrant and the pin.
* **The apex cone test costs a facet enumeration** (≤ 15 cross products over ≤ 6
  de-duplicated generators) on every plastic step. It is exact and cheap in
  absolute terms, but it is not free; P6's measurement decides whether the
  associated case (`mb_psi == mb`, where the cone collapses to
  `C3 (x − apex) ≥ 0`) deserves a short circuit.
* **Region classification can solve the face Newton and then discard it** when
  the trial belongs to an edge. Worst case is three 4×4 Newtons per plastic
  Gauss point instead of one. Also P6.
* `HB_mb_psi > HB_mb` is **refused** rather than supported: the potential's own
  `arg` goes negative on part of the yield surface, where the flow direction does
  not exist. No registered deck does this.
* Mixed HB pairings (eleven specializations) and `StiffSoil` (P5) remain refused.
* `Closest_Point` is still **not** the default and `Algorithmic` is **not** the
  default tangent (D1; the flip is gated on P6).
