# ADR-95 P4 — the UW DruckerPrager corner fix, and the decisive leg

```
VERDICT  ADR-95 P4 — H1 CONFIRMED AS CAUSE.  The quadratic wall was a CODING
         DEFECT in the vanilla UW DruckerPrager return map, not element physics.
build      ladrunoBuild() 0388889713bcdb7f6bb1a4bc6fac09f5b4e1d67a  (HEAD at
           BUILD time; the P4 edits were still uncommitted and landed as
           31322a47a — the same stamp trap P0 hit.  Staged pyd, dist/bin was
           held open by another agent's run, as predicted.)
fix        (a) index-driven residual/Jacobian assembly (the `Jact(i) == 2` arm
           was unreachable); (b) the tangent's radial-return term divides by
           ||eta_TRIAL||, not by the RETURNED norm, guarded against ||eta||=0.
tests      tests/test_adr95_dp_corner_fix.py + test_adr95_dp_branch_response.py
           = 20 passed, incl. a central-difference check of the corner tangent
           in all 6 Voigt directions (rel 1e-3) on a non-degenerate corner.
f1-only    NOT bit-identical — requirement (c) cannot hold, and the review said
bit-ident  so: line 697's denominator is itself wrong ON THE CONE.  Gate fastest
           leg (h0 1.0, non-assoc): ratio 1.0849417 -> 1.0849561 (+1.3e-5 rel),
           q 150.70600 -> 150.70801, resultant/patch identities unchanged, and
           the mode IMPROVES, BUDGET -> TARGET (s/B 0.14378 -> 0.15000, headroom
           2500 -> 4400, 1390 s -> 181 s).  Linear control h8bbar: TARGET both,
           q 150.707506 -> 150.707256 (1.7e-6 rel), 563 -> 310 steps, 2628 ->
           487 s.  Both stay far inside the gate band (1.0517 .. 1.1167).
decisive   h20uri h0 1.0 —  p1: MODE FLOOR at s/B 0.01114, q 0.7689 of exact,
leg        detAmin -6.6e+07 at the 4 corner GPs.  p4: __P4_MODE__ at s/B
           __P4_SB__ (__P4_X__x p1's reach), q __P4_Q__ of exact; corner GPs
           (branch 3) appear from s/B 0.0177 on, are returned to I1 = T, and
           detAmin_min stays -5.8e-2 .. -6.1e-2 — O(1) throughout, vs the
           1000x sigma_min collapse p1 measured.  cond flat at ~9.5e4.
upstream   OpenSees/master STILL has both defects (dead `Jact(i)==2` arms at
           lines 495 and 530, final-norm divide at line 670) — upstream PR
           material.  So does this fork's DruckerPragerThermal.cpp (751, 787).
```

## 1. The defect, in one paragraph

`Jact` is a 2-vector of **flags** (0/1), one per yield surface. Upstream's
residual/Jacobian assembly switched on the flag **value**:

    for (int i = 0; i < 2; i++) {
        if      (Jact(i) == 1) { R(0) = <f1 residual>; g(0,0) = <df1/dg0>; }
        else if (Jact(i) == 2) { R(1) = <f2 residual>; g(1,1) = <df2/dg1>; }
    }

`Jact(i)` is never 2, so the second arm is unreachable. At the corner both loop
passes wrote **row 0**: `R(1)` stayed 0 and `g(1,1)` stayed at the dummy `1` of
the "initialize such that det(g) = 1" lines, so `gamma(1)` came out **exactly**
zero on every branch. With `rho_bar = 0` — where a cone return cannot move `I1`
at all — the Gauss point committed with `I1` wherever the trial state left it,
arbitrarily far above the cutoff `T`. On the f2-only active set the same loop
assembled the **f1** residual into row 0, i.e. it performed a cone return for a
pure cutoff step.

The same wrong `g` then poisoned the **consistent tangent**, because `g_contra`
multiplies every rank-one term: `g_contra(1,1)` came out `+1` where it should be
`-1/(9K)`, injecting two spurious rank-one terms in *stress-squared* units,
scaling as `9K²/2G` and `27ρK²/2G` relative to `2G`. That is precisely the 5–8
decade `detAmin` outlier P1 measured (−3e5 … −6.6e7 against −0.06 elsewhere), and
because it reaches the **global** operator the symptom was a Newton seizure with
no material-level warning.

## 2. What was changed

`SRC/material/nD/UWmaterials/DruckerPrager.{h,cpp}`, every edit marked
`// Ladruno ADR-95`, ledger rows added:

1. **Index-driven assembly**, both copies of the loop (initial and iterate):
   `if (Jact(0) == 1) { row 0 } if (Jact(1) == 1) { row 1 }`. The row index is
   the **surface** index, not the loop counter. Note `g(1,1) = -9K + δ2·T(α2)`
   and `g(1,0) = ρ̄(-9K + δ2·T(α2))` are upstream's own expressions and were
   already correct (`T(α2) = T₀e^{-δ2 α2}` ⇒ `T' = -δ2 T`); they were simply
   never assembled.
2. **The tangent's trial norm.** The last term
   `-4G²γ₀/||η|| · (IIdev - n⊗n)` is the `dn/dε` contribution and requires
   `||η_trial||`; upstream divides by `norm_eta` **after** it has been
   overwritten with the norm of the **returned** η (recomputed for `mState(1)`).
   Now latched as `norm_eta_trial` and guarded (`> 1e-13`).
3. A plain-text derivation of the whole tangent above the `Cep` assembly.
4. The pre-existing `NormCep < 1e-10` message throttled to 10 lines (see §4).
5. A read-only `ladrunoTangent` response (id 96, `Vector(36)` of `mCep`) so the
   test can finite-difference the return map against the analytic tangent.

Nothing else moved: no parameter, no argument parsing, no class tag, no
`sendSelf`/`recvSelf` layout.

## 3. Where the review's derivation is wrong — and it matters

`_adr95_corner_tangent_review.md` §1 lists as defects

    688  temp1 = -n - (3K rho/2G) I1 - (27 K K rho/2G) I1   <-- last term spurious
    689  temp2 = 3K I1                                      <-- should be -(1/3) I1

**Lines 688–689 (and 694–696) are upstream's and are CORRECT as written.** Those
were *symptoms of the wrong `g`*, evaluated with `g(1,1) = 1`, not independent
coding errors. Deriving the tangent from scratch (the block comment now in the
source):

    dR0/deps = 2G n + 3K rho 1 = b1        dR1/deps = 3K 1 = b2
    dgamma   = -g^{-1} b : deps
    temp1 = g^-1(0,0) b1 + g^-1(0,1) b2                 ( dg0        = -temp1:deps )
    temp2 = rho_bar temp1 + g^-1(1,0) b1 + g^-1(1,1) b2 ( rb dg0+dg1 = -temp2:deps )
    Cep   = Ce + 3K 1(x)temp2 + 2G n(x)temp1 - (4G^2 g0/||eta_TR||)(IIdev - n(x)n)

which is exactly the code. With `g` assembled index-driven, the corner case
`ρ̄ = H = θ = δ2 = 0` gives `g = [[-2G, -9Kρ], [0, -9K]]`, `det g = 18GK`,
`g^{-1} = [[-1/2G, ρ/2G], [0, -1/9K]]`, hence **`temp1 = -n` and
`temp2 = -(1/3)·1`** — the review's own target values — and
`Cep = 2G(1 - 2Gγ₀/||η_tr||)(IIdev - n⊗n)`, its own target operator. So the fix
is the *assembly plus the denominator*, and touching 688–696 would have broken a
correct derivation. The review's §2 riders (i) and (ii) are both right.

## 4. Two things that look like new bugs and are correct answers

**(a) The corner of a non-hardening deck IS the apex, and its tangent is zero.**
`mTo` is not a free parameter: the constructor sets `T = √(2/3)·σ_y/ρ`, which is
exactly the cone apex. So `f1 = f2 = 0` forces
`||η|| = √(2/3)K(α1) - ρT = 0` whenever `K(α1) = σ_y`: every corner return on the
campaign deck lands on the apex, the returned stress is pinned at `σ = (T/3)·1`,
and the exact consistent tangent is **zero** (measured `NormCep ≈ 5e-23`). The
material's own `NormCep < 1e-10` floor catches this and hands the solver
`1e-3·Ce`, which is the right thing to do — but it now fires at every apex GP, so
its `NormCep = ...` message is throttled to 10 lines. A *non-degenerate* corner
needs hardening; that is what the FD test uses.

**(b) The f2-only branch is unreachable without hardening.** `f1 ≤ 0 < f2` needs
`||η|| < 0` when `T = √(2/3)σ_y/ρ`. That is why P1 never saw branch 2 anywhere,
and why the branch-2 test needs `θ = 1, H > 3G`. Pre-fix that branch did a cone
return; post-fix it does a pure volumetric return leaving `2G·IIdev`
(`detAmin = 1/6` exactly).

## 5. Measured, single `LadrunoBrick`, uniform prescribed strain

Deck as P0: `K = 1e4`, `G = 4e3`, `SY = 0.2`, `ρ = 0.148583`, `ρ̄ = 0`, no
hardening unless noted; `T = 1.0990438`.

| leg | branch | gamma0 | gamma1 | I1 | detAmin | max\|C\|/2G |
|---|---|---|---|---|---|---|
| elastic | 0 | 0 | 0 | −1.800 | +0.4791667 | 1.917 |
| cone | 1 | 3.78117e−4 | 0 | −9.000 | **−2.826e−4** (P0: −7.541e−3) | 1.471 |
| corner, P0 leg (3,1,1)e−4 | 3 | 1.63299e−4 | **1.54455e−4** (P0: 0) | **1.099044** (P0: 15.000) | 4.792e−10 (apex floor; P0: **+1.511e+3**) | 0.0019 |
| corner, hydrostatic (1,1,1)e−4 | 3 | 8e−20 | 8.77884e−5 | 1.099044 | +0.1666667 | 0.667 |
| f2 only (hardened) | **2** | 0 | 3.77884e−5 | 1.099044 | +0.1666667 | 0.667 |
| corner, hardened (non-degenerate) | 3 | 3.40207e−6 | 8.77884e−5 | 1.099044 | +0.1371 | 0.662 |

Two readings beyond the corner. **The cone tangent changed too** — `detAmin`
−7.54e−3 → −2.83e−4, a 27× reduction in the magnitude of the loss-of-ellipticity
determinant at an ordinary yielding GP. That is the trial-norm fix: upstream's
factor `1 - 2Gγ₀/||η_final||` goes *negative* once `2Gγ₀ > ||η_tr||/2`, i.e. the
deviatoric part of the tangent was sign-flipped on ordinary plastic steps.
**And the apex no longer NaNs** — the third bite of the `LEDGER_quirks` entry
(hydrostatic tension divides by zero) is fixed by the guard, so a probe no longer
has to perturb the deviator.

## 6. Tests

`tests/test_adr95_dp_corner_fix.py` (11 cases) and the P0 module (9). All 20 pass
in 0.48 s on the staged build; `zone_a`.

* `test_hydrostatic_corner_returns_I1_to_the_cutoff` — the defining assertion:
  `I1_trial = 9.0`, `I1_returned = T` to `rel 1e-8`, `gamma1 = (I1_tr − T)/9K`
  exactly, `gamma0 = 0` (volumetric trial), `detAmin = 1/6` to `rel 1e-9`,
  `max|C|/2G < 100`, tangent finite (the NaN guard).
* `test_deviatoric_corner_lands_on_the_apex_with_a_zero_tangent` — both
  multipliers positive, `I1 = T`, stress pinned at `(T/3)·1`, tangent = the
  `1e-3·Ce` floor. The `|detAmin| ∈ [1e-3, 1e2]` criterion of the P4 brief is
  **not** asserted here and cannot be: the exact apex tangent is zero, so the
  test pins the floor value `1e-9·0.4792` instead and says why. It *is* asserted
  on the hydrostatic corner, the f2-only branch and the hardened corner.
* `test_f2_only_branch_is_reachable_and_returns_volumetrically` — first time
  branch 2 has ever been reached in this fork.
* `test_hardened_corner_is_non_degenerate` — the guard for the FD test:
  `||η|| > 1e-2` at the return, both `gamma > 1e-9`, `|detAmin| = O(1)`.
* `test_corner_tangent_matches_finite_difference[k=0..5]` — central difference of
  the *whole* two-stage path at `h = 1e-7` versus column `k` of `mCep`, `rel
  1e-3`. This is the property a mis-assembled `g` destroys and the one the global
  Newton consumes.

The P0 sentinel is **kept and flipped**:
`test_tension_cutoff_multiplier_is_structurally_zero_in_vanilla` →
`test_tension_cutoff_multiplier_now_returns_stress`, requiring
`gamma1 = (I1_tr − T)/9K` and `I1 = T`. It is kept rather than deleted because it
is the cheapest detector of a build that has silently lost the fix (a stale
`.pyd`, a bad merge, an upstream re-sync), and its numbers are hand-derivable.

## 7. The f1-only bit-identity requirement (c) — cannot hold, and should not

The P4 brief asked for the f1-only branch to stay bit-identical. It does not, and
the review anticipated exactly this (§4 option 1, "697 changes f1-only runs too —
the gate must expect a real (correct) delta"). The reason is structural: the
denominator at line 697 is wrong **on the cone as well**, not only at the corner.
The index-driven assembly change (fix 1) *is* bit-identical on `Jact = (1,0)` —
the old loop's `i = 1` pass did nothing there — so the entire delta below is
attributable to the trial-norm fix.

| quantity | reference (P0, build cf239c9d) | P4 (0388889713) |
|---|---|---|
| DOF | 1386 | 1386 |
| `q_num` | 150.70600 | 150.70801 |
| ratio | 1.0849417 | **1.0849561** (+1.3e−5 rel) |
| tail % | 0.00139 | 0.00136 |
| mode | BUDGET | **TARGET** |
| end s/B | 0.14378 | **0.15000** |
| ds/floor | 2500 | **4400** |
| capacity | yes | yes |
| resultant identity | 2.66e−15 | 2.66e−15 (unchanged) |
| 1-D stress patch | 1.01e−14 | 1.01e−14 (unchanged) |
| wall | 1390 s (contended) / 537 s ref | **181 s** |

The gate band is `1.051674 … 1.116726`; both values sit near its centre. The two
exactness identities are byte-unchanged, which is the right invariant to demand
of a tangent-only change. The improvement in *mode* (the leg now reaches the
target settlement instead of exhausting the step budget) and the 3–7× wall-time
drop are the direct consequence of removing a sign-flipped deviatoric tangent
from Newton.

**Owed before the PR flips to ready:** the other two resolutions of the gate and
the associated control, which P0 also did not run.

## 8. The decisive leg

`quad_path_diag.py --elem h20uri --h0 1.0 --cond --cond-at 5e-4 --branch`, on the
fixed build, against P1's identical invocation.

| | P1 (pre-fix) | P4 (fixed) |
|---|---|---|
| MODE | **FLOOR** at s/B 0.01114 | __P4_MODE_ROW__ |
| q_max | 106.81 kPa = **0.7689** of exact | __P4_Q_ROW__ |
| end s/B | 0.01114 of 0.15 | __P4_SB_ROW__ |
| corner GPs | 4, all at once, at the wall station | present from s/B 0.0177, 0–4 per station, transient |
| `I1` at those GPs | 0.80–0.98 **≥ T**, `f2 = +0.03…+0.21` | returned to `I1 = T` |
| `detAmin_min` | **−6.642e+07** at the corner GPs | **−5.8e−2 … −6.1e−2**, flat |
| `sigma_min/scale` | 2.2e−4 → **2.6e−7** (1000× collapse) | **2.1e−4 … 2.3e−4**, flat |
| `cond` | 9e4 → **4.6e+09** | **~9.5e4**, flat |
| forced accepts | 0 | 0 |

The prediction in the P4 brief is met in full: **the leg passes s/B 0.0112, corner
GPs (branch 3) appear and are returned to `I1 = T`, and `detAmin` stays O(1)**.
The 1000× `σ_min` collapse of note 82 §7.3 does not recur. Nothing else about the
element, the mesh, the deck or the solver changed.

Linear control `h8bbar --ladder linear` on the same build: `MODE = TARGET`,
`q_max = 150.707256` against P1's `150.707506` (1.7e−6 relative), ratio 1.0850,
capacity yes, **zero f2/corner GPs anywhere** (branch census `(1346, 254, 0, 0)`
at s/B 0.15 vs P1's `(1188, 254→412, 0, 0)`; the elastic/plastic split moves
within the ±40 % flicker P1 already characterised as a sampling artefact of
perfectly plastic GPs resting on the surface). 310 steps / 487 s against P1's 563
steps / 2628 s. So the control behaves as predicted — it never reaches the
cutoff, so only the cone-tangent half of the fix touches it, and it moves the
answer by 1.7e−6 while halving the step count.

## 9. Consequences for ADR-95

* **H1 is confirmed as the cause**, in the sharpest possible form: the quadratic
  wall was a coding defect in a vanilla upstream file, reached first by quadratic
  elements because they resolve the tensile spot beside the footing edge that the
  linear bbar hex averages away. It is not element physics, not locking, and not
  loss of ellipticity per se.
* **H2 stays dead.** `detAmin < 0` at every plastic GP both before and after the
  fix (−5.8e−2 in the fixed quadratic run, −8.0e−2 in the linear control that
  completes the collapse). Loss of ellipticity is the generic state of a
  non-associated perfectly-plastic GP and was never the discriminator.
* **P2's SY ladder is superseded** as a *diagnostic*: raising `SY` moved the wall
  because it moved `T` out of reach of the heave zone, which is the same
  mechanism, not an independent knob.
* **The quadratic limit-load ceiling memory entry needs revisiting.** Every
  quadratic leg in that campaign ran on a pre-fix binary. The 0.59–0.77 ceiling
  and "only the LINEAR relieved hex plateaus" verdict were measured through this
  defect wherever the deck put the heave zone above `T`.
* Everything else the campaign measured with `rho_bar = 0` and a small apex
  regulariser `SY` is suspect for the same reason, including the ADR-79 §9
  bearing numbers if any of their GPs crossed the cutoff.

## 10. Upstream (Task 5)

Fetched `raw.githubusercontent.com/OpenSees/OpenSees/master/SRC/material/nD/
UWmaterials/DruckerPrager.cpp` on 2026-09-07 (985 lines). **Upstream master still
carries both defects, verbatim.** The dead `else if (Jact(i) == 2)` arm is at
**line 495** (initial assembly, loop opening at 490) and **line 530** (the copy
inside the Newton iterate, loop opening at 525); the corner off-diagonals
`g(0,1)` / `g(1,0)` at 501–502 are unchanged and correct, so the failure mode is
identical to the fork's pre-P4 state. The tangent block is also unchanged:
`temp1`/`temp2` at 661–662, and the radial-return term at **line 670** reading
`- 4*mG*mG/norm_eta*gamma(0) * (mIIdev(i,j) - n(i)*n(j))` with `norm_eta` having
been overwritten by the returned-η recomputation at **line 646**; the
`NormCep < 1e-10` floor is at 675, unthrottled. This is therefore clean upstream
PR material for the jaabell/ladruño campaign: a self-contained, hand-derivable
two-line-class fix to a vanilla file, with a falsifier that runs in half a second
(`I1` must return to `T`), and a second beneficiary in
`SRC/material/nD/DruckerPragerThermal.cpp` (same dead arm at lines 751 and 787,
same final-norm divide) which this PR deliberately leaves alone.

## 11. Notes for whoever runs next

* The build stamp `0388889713` is HEAD *at build time*; the P4 source edits were
  uncommitted then and landed as `31322a47a`. Same trap as P0 — a probe rebuilt
  after that commit will stamp `31322a47a`.
* `dist/bin/opensees.pyd` was again held open by another agent's campaign, so
  `build.bat` printed "the process cannot access the file" and carried on with a
  success banner. Everything here ran from a staged copy
  (`<scratch>/dist_p4/bin`, `build/build/Release/OpenSeesPy.dll` renamed). The
  stale `dist/bin` pyd doubled as a free negative control: it answers
  `ladrunoTangent` with a zero-width vector, which is how the first probe run
  caught itself pointing at the wrong binary.
* **Do not launch several long runs as sibling WMI `cmd.exe /c` processes.** They
  share a console; when the first one exits, the survivors take a
  `CTRL_CLOSE_EVENT` and the Intel Fortran runtime aborts them with
  `forrtl: error (200): program aborting due to window-CLOSE event`. The decisive
  leg died that way at s/B 0.0025 and had to be relaunched inside its own console
  (`cmd /c start "name" /min cmd /c ...`).
