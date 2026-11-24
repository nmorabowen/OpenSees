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

Prototype, 2026-11-24. Not committed. 14 new unit tests and the 130 existing LadrunoForceBeamColumn tests pass on
Linux (esmeralda build). Lab results in §5.

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
  state; a trial state left by a poor global iterate can be far off). The retry costs a second solve exactly
  where solves fail, and on 2CLD12 it leads onto a branch that grinds (§5), so the lab driver keeps it off
  and switches it on only as a last resort, through `parameter ... element <tag> monoRestart` /
  `updateParameter` (0/1).
- `-monoStall nIter [ratio]` (default 6, 0.99): stop a Newton whose `||R||` fell by less than 1 % over the last
  6 iterations. Meant to make failing solves cheap (the global Newton calls the element once per global
  iteration, up to ~200 times per failed attempt with fallback algorithms; Model C once per step). It also
  cuts off slow solves that would have converged and changed the path on 3CLH18 and 2CLD12M (§5): the lab
  runs use `-monoStall 0`. Also a parameter (`monoStall`).
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

Model C's driver (2e-4 rad sub-steps, halving to 1/64, fallback algorithms, 10 kN acceptance, ladder with the
overshoot predictor), element `-monolithic -monoAccept 10e3 10e6 1e-2 -monoStall 0`, restart switched on only to
retry a record point about to be given up:

- **7 of 8 complete**, including 3CLH18 and the four WP-254 stoppers (2CLD12, 2CLH18, 3CMD12, 3CMH18). No run
  stops on a singular flexibility. 2CLD12M stops at 771/791 in its final collapse; with the restart always on
  it completes (e 0.155, Model C 0.154), but that setting sends 2CLD12 onto a branch that grinds.
- Lab score e equals Model C's on 2CLD12, 2CHD12, 2CMH18, 3CMD12 and within 0.005 on 3CMH18, 2CLH18; 3CLH18
  0.140 (Model C 0.171). J 0.190 over 8 with 2CLD12M incomplete, ~0.150 with its restart run (Model C 0.155).
- 57–910 s per specimen (Model C's own code: 1,000–5,400 s).

What remains is the strip law's folds (its branch rule makes `s(e)` jump), crossed with acceptance and the
ladder, and the two-level cost (the element is solved per global iteration). Details and every configuration
tried: juanito-frames `research/11_wp255_report.md`.

## 6. Not done before a PR

- `sendSelf`/`recvSelf` do not carry the mode (an OpenSeesMP copy runs nested). Add the 12 values to the
  data vector and a database round-trip test.
- 2D element (`LadrunoForceBeamColumn2d`): same change, not needed by the program.
- Guide/manual entry, ledger row, banner line, stamp headers.
- Windows build and the Windows test run (this prototype was built and tested on Linux only).
