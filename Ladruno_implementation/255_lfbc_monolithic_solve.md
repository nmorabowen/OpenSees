---
title: WP-255 — LadrunoForceBeamColumn3d -monolithic (Model C's element solve)
project: Ladruno
tags:
  - element
  - force-based
  - mcft
  - wp-255
---

# WP-255 — monolithic element state determination for LadrunoForceBeamColumn3d

Program #1131 (MCFT shear-critical fibre beam). Follows WP-254 (`LadrunoMCFTStrip`, PR #1149), whose
material reproduces the reference model ("Model C") to round-off but whose column stops on 4 of 8 lab
tests after the shear failure. Lab tests: Sezen & Moehle 2002 (2CLD12, 2CHD12, 2CLD12M, 2CLH18, 2CMH18),
Lynn et al. 1996 (3CLH18, 3CMH18, 3CMD12).

## Status

Prototype, 2026-11-24 (draft PR). 14 new unit tests and the 130 existing LadrunoForceBeamColumn tests pass on
Linux (esmeralda build). Lab gate met with the defaults: 8/8 complete, J 0.149 (§5).

## 1. Problem

`LadrunoForceBeamColumn3d::evalBulkState` is the classical nested force-based state determination
(Spacone et al. 1996; Neuenhofer & Filippou 1997): for a target `v` it iterates

```
section:  e_i += f_s,i (b_i q - s_i(e_i)),   f_s = K_s^-1
element:  q   += F^-1 (v - v_r),             F = sum w_i b_i^T f_s,i b_i
```

After a shear failure every MCFT strip of a section can be on a flat branch and every bar past fracture.
`K_s` is then exactly singular: `f_s` does not exist, `F` cannot be inverted, and the iterate runs away
(section strains ~1e14 seen in WP-254). Every WP-254 stop is
`evalBulkState() -- could not invert flexibility`. Making that branch fall through (issue filed separately)
does not complete the runs: the iteration itself is the limit.

Model C never forms `f_s`: its element solve is one Newton on all section deformations and the end forces,
with the section tangents and a backtracking line search.

## 2. Decision

An **opt-in mode of the existing element**, not a new element: `-monolithic [maxIter [tol]]`. Every other
feature (integration rules, sections, `-shear`, member loads, recorders, `revertToLastCommit`, the
`update()` early returns) is shared; only the state determination changes. Default behaviour is untouched.

Unknowns and residual (`n` sections of orders `o_i`, `x` of size `sum o_i + 6`):

```
x   = (e_1, ..., e_n, q)
R_i = (s_i(e_i) - b_i q - s_p,i) / sc_i                        section equilibrium
R_c = (sum_i w_i L b_i^T e_i + v_sh(q) + v_T(q) - v) / sc_c    compatibility
J   = [ K_s,i / sc_i        -b_i / sc_i         ]
      [ w_i L b_i^T / sc_c   dv_sh/dq / sc_c    ]
```

- `s_p,i`: member loads (`computeSectionForces`), `v_sh`: `-shear GAvy GAvz` flexibility, `v_T`: the default
  torsion flexibility `DefaultLoverGJ` when no section carries T — the same terms as `evalBulkState`.
- Scales (Model C): force rows (P, Vy, Vz) `1e5`, moment rows (Mz, My, T) `1e8`, compatibility `1e-5`; the
  axial row `du` is scaled as a mean strain (`1e-5 L`) since Model C prescribes N and has no such row.
  `-monoScale sForce sMoment sCompat` overrides (units of the model).
- Newton from the last **trial** solution (`vs`, `Se`), `dx = J^-1 (-R)` by LAPACK `dgesv`; backtracking
  line search, `lambda` halved up to `-monoHalve n` times (12), accepted on `||R_t|| < (1 - 1e-4 lambda) ||R||`;
  converged on `||R||_2 < tol` (1e-5), at most `maxIter` (40) iterations. Exactly `mcft_oracle.Element.solve`.
- `-monoAccept aForce aMoment [aCompat]` (off by default): Model C's acceptance of a stagnated solve — every
  section force residual `< aForce`, moment residual `< aMoment`, scaled compatibility `< aCompat` (1e-2).
  Model C's calibrated runs use 10 kN / 10 kN m.
- **Start points.** Newton starts from the last trial solution. With `-monoRestart 1` (default) a failed solve
  is retried from the committed state (Model C restarts every solve from its last converged, committed
  state; a trial state left by a poor global iterate can be far off). Near a fold it decides which branch the
  solve lands on: on 2CLD12M only the restart drops onto the failed branch at the shear failure, as Model C
  does (§5). It is the lab gate setting. Also a parameter (`parameter ... element <tag> monoRestart`,
  `updateParameter` 0/1) so a driver can switch it per attempt.
- `-monoStall nIter [ratio]` (default 0 = off; ratio 0.99): stop a Newton whose `||R||` fell by less than 1 %
  over the last `nIter` iterations. Meant to make failing solves cheap (the global Newton calls the element once per global
  iteration, up to ~200 times per failed attempt with fallback algorithms; Model C once per step). It also
  cuts off slow solves that would have converged and changed the path on 3CLH18 and 2CLD12M (§5), hence off
  by default. Also a parameter (`monoStall`).
- **Element tangent** `kv = dq/dv` from the same matrix with six right-hand sides,
  `J X = B`, `B(c,k) = 1/sc_c,k` on the compatibility rows; `kv = X(q rows, :)`. This is the condensation
  `(sum w b^T K_s^-1 b)^-1` without forming any `K_s^-1`. Only the whole block matrix must be non-singular:
  a section with `K_s = 0` is fine as long as the element has a stiff path through the others.
- A section is never inverted: in this mode `revertToLastCommit` no longer calls `getSectionFlexibility`.
- `eleResponse monoStats` = [on, last iterations, last scaled residual, accepted solves (calls), largest
  accepted force residual].
- Refused with embedded hinges (`-hinge*`, which have their own Schur solve) and `-axialRelease`.

### Considered

1. **New element** (`LadrunoMonoBeamColumn3d`): cleaner file, but duplicates ~6,000 lines of shared
   behaviour (loads, recorders, sensitivity, commit/revert bookkeeping). Rejected.
2. **Harden the nested loop only** (fall through on a singular `F`, bound the inner iterate): needed anyway,
   filed as its own issue; WP-254's experiment shows it does not complete the runs.
3. **Displacement-based column** (several DB elements): no inner loop, but changes Model C's integration, so
   the calibration does not carry over.
4. **Monolithic with N prescribed inside the element** (Model C's exact unknown set, `du` free): impossible
   in an element that receives all six `v` from the global Newton; the global Newton solves the axial DOF
   instead. Same converged state; the outer loop sees a consistent `dq/dv`.

## 3. Equivalence and differences with Model C

- Converged states solve the same equations as the nested iteration (gate: §4, a cyclic biaxial fibre column,
  `q` within 1e-6) and as Model C's `solve_mono`.
- Two-level Newton: Model C solves `(e, q2, q3)` with N fixed in one system; here the element solves for the
  global iterate's `v` and the global Newton moves `du`. A bad global iterate can ask the element for a `v`
  with no nearby solution; the driver then sub-steps.
- No `_rebranch` (Model C's own runs did not use it either); the strip's branch rule makes `s(e)` jump at
  folds, so the element Newton can stagnate there. Model C crosses those with acceptance and its ladder; so
  does the lab driver (§5).
- Section tangent: per-fibre forward differences inside `LadrunoMCFTStrip` (WP-254), not Model C's
  section-level differences. Newton paths differ in round-off; converged states agree to the tolerance.

## 4. Tests

`tests/test_ladrunoForceBeamColumn3d_monolithic.py` (zone A, 14):
- monolithic vs nested on a Concrete02/Steel02 fibre column with aggregated Vy, Vz, axial load and a cyclic
  biaxial drift path: basic forces within 1e-6 relative at every step;
- elastic section with `-shear GAvy`: tip displacement = `PL^3/3EI + PL/GAvy` to 1e-9;
- the returned tangent is finite and symmetric for a symmetric section;
- parser: refuses bad `-monolithic`, `-monoScale`, `-monoAccept`, `-monoHalve`, `-monoRestart`, `-monoStall` values and
  `-axialRelease -monolithic`; accepts the full option set; `monoStats` on and off.

## 5. Lab gate (program #1131)

Model C's driver (2e-4 rad sub-steps, halving to 1/64, fallback algorithms, 10 kN global acceptance, ladder
with the overshoot predictor) and the element defaults plus Model C's acceptance:
`-monolithic -monoAccept 10e3 10e6 1e-2` (restart from committed on, stall exit off). Scored against the
Model C reference and the tests:

| Specimen | steps | wall s | e fork | e Model C | d80 fork / Model C / test |
|---|---|---|---|---|---|
| 2CLD12 | 2385/2385 | 341 | 0.083 | 0.083 | 1.42 / 1.42 / 0.93 |
| 2CHD12 | 1639/1639 | 81 | 0.168 | 0.168 | 0.92 / 0.92 / 0.59 |
| 2CLD12M | 791/791 | 157 | 0.155 | 0.154 | 2.83 / 2.83 / 2.83 |
| 3CLH18 | 374/374 | 167 | 0.134 | 0.171 | 0.52 / 0.52 / 0.52 |
| 2CLH18 | 737/737 | 243 | 0.303 | 0.307 | 0.78 / 0.78 / 1.42 |
| 2CMH18 | 461/461 | 54 | 0.052 | 0.059 | — |
| 3CMH18 | 987/987 | 94 | 0.166 | 0.165 | 0.52 / 0.52 / 0.52 |
| 3CMD12 | 1147/1147 | 106 | 0.131 | 0.131 | 1.10 / 1.10 / 0.82 |

**8/8 complete, J 0.149 (Model C 0.155; gate 0.17), d80 within one level 6/7, every specimen < 6 min**
(Model C's own code: 1,000–5,400 s). No run stops on a singular flexibility.

Branch selection at the shear failure: on 2CLD12M at 56 mm Model C snaps from 254 to 66 kN. Without the
restart, the fork crosses that point on a 10 kN-accepted step, stays on the unfailed branch, which folds near
60 mm, and stops at 771/791. With it, the fork drops at the same record point and follows Model C within
5–7 kN. An earlier "grind" of 2CLD12 with the restart on was a lab-driver bug (a half-successful overshoot
predictor left the node past its target and the step direction kept the old sign), not the element.

## 6. Not done before a PR

- `sendSelf`/`recvSelf` do not carry the mode (an OpenSeesMP copy runs nested). Add the 12 values to the
  data vector and a database round-trip test.
- 2D element (`LadrunoForceBeamColumn2d`): same change, not needed by the program.
- Guide/manual entry, ledger row, banner line, stamp headers.
- Windows build and the Windows test run (this prototype was built and tested on Linux only).
