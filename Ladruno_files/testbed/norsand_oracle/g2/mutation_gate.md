# WP-144 G2 shell mutation gate

Rounds: 2026-10-01 (M1 to M7, branch `wp/144-ladruno-norsand` at ab6503a0c) and 2026-10-02 (MX1 to MX3, section "Round 3" below). Nothing committed.
Kernel mutants are covered separately (`kernel_parity/mutate_kernel.sh`, 15/15, re-run 2026-10-02 on Esmeralda). This gate mutates the SHELL
(`SRC/material/nD/LadrunoNorSand.cpp`) one mutant at a time.

Procedure per mutant: apply the edit, full `build.bat` from a fresh cmd.exe
(`call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat`, all 5 targets, 60-120 s incremental),
confirm `dist\bin\OpenSees.exe` and `opensees.pyd` mtimes moved, run the G2 files, record failures,
`git checkout` the file, confirm `git diff -- SRC` is empty.

Test sets (all with `--runslow`):
- Zone A: `tests/test_ladruno_norsand.py` + `tests/test_ladruno_norsand_element.py` (40 tests, ~95 s).
- Zone B: `Ladruno_files/testbed/norsand_oracle/g2/test_g2_shell_parity.py` + `test_g2_logstrain.py`
  (91 passed + 6 strict xfail, ~50-70 s). RUN IT FROM THE g2 DIRECTORY. Run from the repo root with
  `tests/conftest.py` in play it skips every case ("gmsh not installed"), which would make Zone B look green
  while testing nothing.

Baseline (unmutated): A 40 passed, B 91 passed + 6 xfailed.

## Result: 7 of 7 killed, no survivors

| # | Mutant (edit) | Zone A failures | Zone B failures | Named killer(s) |
|---|---|---|---|---|
| M1 | `getTangent`: shear-column factor 0.5 -> 1.0 | 15 | 66 | `test_elastic_response_and_initial_tangent_match_the_closed_form`; `test_brick_assembled_tangent_vs_fd[*]` (4) + `_perturbed_nonuniform`; `test_planestrain_quad_assembled_tangent_vs_fd[*]` (4); `test_newton_is_quadratic_*`; Zone B `test_shell_matches_o2_on_path` (54 of 54 plastic/tangent paths), `test_gate_discriminates_mapping_mutants`, `test_wrapper_spatial_tangent_equals_o2_s34_assembly` |
| M2 | `unpackState`: `v0 := v` on recv | 2 | 0 | `test_database_roundtrip_carries_the_committed_state_with_v0_distinct`, `test_database_roundtrip_with_committed_differing_from_trial` (Zone A only) |
| M3 | `revertToStart`: v0 not restored (takes the pre-revert committed v) | 1 | 1 | `test_revert_to_start_restores_the_initial_state_including_v0`; Zone B `test_reset_replay_is_bitwise` |
| M4 | getCopy clones share the committed state with the source (static map keyed by tag; written on commit, read before integrate) | 3 | 2 | `test_getcopy_clones_are_history_isolated`; also `test_discarding_element_aborts_the_commit_and_latches_until_revert`, `test_revert_to_start_restores_the_initial_state_including_v0`; Zone B `test_reset_replay_is_bitwise`, `test_specific_volume_is_v0_times_detF_in_the_wrapper_and_both_oracles` (incidental) |
| M5 | refusal latch never set (`latched = true` removed in `commitState`) | 3 | 2 | `test_discarding_element_aborts_the_commit_and_latches_until_revert`, `test_latch_and_warnings_carry_the_finest_refusal_reason`, `test_latch_survives_a_database_round_trip`; Zone B `test_refusing_paths_refuse_through_the_shell` (2) |
| M6 | `getTangent` returns the hyperelastic tangent at the trial state | 13 | 66 | same set as M1: FD brick/quad gates, Newton order tests, `test_database_roundtrip_carries_...`, Zone B shell parity (54 paths), LogStrain (S.34 assembly, min-det curve, K2 localization) |
| M7 | `getTangent` symmetrised, 0.5(T+T^T) | 13 | 59 | FD brick (4) + perturbed, FD quad (4), Newton order tests (3), `test_database_roundtrip_carries_...`; Zone B shell parity (48 paths + summary + smooth-cap AMP chained tangent), LogStrain (S.34 assembly, min-det curve, K2 localization step) |

## Final state

Source reverted after each mutant (`git diff -- SRC` empty, checked after every one and at the end).
Full `build.bat` (5 targets, no errors) rebuilt from the clean source; all five binaries have fresh mtimes and the
OpenSees.exe / opensees.pyd sizes equal the baseline. Re-run: Zone A 40 passed, Zone B 91 passed + 6 xfailed.
"Clean-state" here means clean SOURCE state with a full 5-target build; I did not run `build.bat clean`
(it wipes build/ and dist/ and re-resolves Conan; the incremental route above is the fork's documented route
for source edits).

## Coverage notes (not survivors, but thin spots)

1. M2 is killed by Zone A only. No Zone B test goes through a database round trip.
2. M4 is killed by exactly one test that targets it directly (`test_getcopy_clones_are_history_isolated`).
   The Zone B shell-parity paths drive a single brick with uniform strain, so all 8 Gauss points hold the same
   history and a shared committed state is invisible to them; the Zone B hits are incidental. The brick/quad FD
   and Newton-order tests also pass under M4, because a shared state keeps residual and tangent mutually
   consistent. A non-uniform-strain element path compared to O2 would add an independent kill.
3. M3 as literally worded ("leave v0 alone") is an EQUIVALENT mutant: v0 is set only in `initialState` and the
   kernel never changes it, so the committed v0 always equals the initial v0. The mutant used makes v0 follow the
   pre-revert committed v, which is what a v0-not-restored bug looks like through the observable state.
4. M1 mutates `getTangent` only. `getInitialTangent` carries the same 0.5 factor and was not mutated; the
   closed-form test `test_elastic_response_and_initial_tangent_match_the_closed_form` checks it directly, but that
   specific edit was not run.
5. M7 makes the slow Newton tests expensive (Zone A 530 s instead of 95 s: the lost quadratic convergence runs
   long refusal ladders). Cost, not a correctness issue.

## Round 3 (2026-10-02): the provider route, the fallback and the commit counter

Branch `wp/144-ladruno-norsand` at 7089b394c plus the uncommitted round-3 changes (the `LadrunoElasticStrainProvider`
mixin and its use in `LogStrainNDMaterial::setTrialF`; `Domain::commit` zeroing the commit-refusal counter before the
element loop). Nothing committed. Three new mutants, MX1 to MX3, each with a full `build.bat` (all 5 targets), plus a
control, MX0.

THE GATES UNDER TEST (all in `g2/test_g2_logstrain.py`; expected values from closed forms written into the docstrings
before the code was run, tolerances from the task):
- `test_rigid_rotation_is_objective_with_pressure_dependent_shear_modulus[2.0, 50.0]`: the two former
  `xfail(strict)` "owner decision 2 pending" tests, now REAL. A 0.2 rad rotation after 10 plastic steps: stress
  `R sigma R^T` to 1e-10 relative to max|sigma|, `pi_i` to 1e-10 |pi_i|, `v` to 1e-12; a second oblique rotation
  (0.35 rad) and a held step again at 1e-10. Measured on the final build: the stress at most 5e-14 at every stage, `pi_i` unchanged exactly, `v` 4e-16.
- `test_rigid_rotation_is_objective_for_constant_shear_modulus` (alpha0 = 0 control, same extra stages).
- `test_provider_identity_committed_be_is_exp_two_eps_e[0, 2, 50]`: the committed b^e of `LogStrain(LadrunoNorSand)`
  equals `exp(2 eps^e)` of the material's own `elasticStrain`, read through a held step after an oblique rotation (all
  three shear components non-zero), `expm` on 3x3 tensors, 1e-12 relative. Measured 4e-16 to 1e-15.
- `test_non_provider_inner_keeps_the_d0_inversion_fallback[ElasticIsotropic, LadrunoJ2]`: Hencky closed form
  (1e-12), objectivity (1e-10), committed b^e = `exp(2 C tau)` with C the isotropic compliance (1e-12). Measured
  1e-14 to 1e-13. The J2 path is asserted plastic.
- Zone A `tests/test_ladruno_norsand.py::test_getcopy_of_a_latched_source_is_not_latched` (MX3).

Run recipe: py312g2 venv, `PYTHONPATH` = a one-line `sitecustomize.py` adding `dist\bin` to the DLL search path, and the
stale `tests\opensees.pyd` copy removed (BUILD_GOTCHAS 4b: it was a copy of the baseline build and would have shadowed
every mutant build). Sets: A = `tests/test_ladruno_norsand.py` + `_element.py` + `_k1.py` (59 tests, about 4 min);
B = `g2/test_g2_shell_parity.py` + `test_g2_logstrain.py`, run FROM `g2/` (113 tests, about 100 s); L = the LogStrain /
finite-strain regression (`test_logstrain_plastic_protocol`, `_reference`, `_tangent_and_j2`, `test_logstrain2d`,
`_2d_reference`, `test_finite_strain_L1_analytical`, `test_ladrunoJ2_finite_element`, `test_ladrunoJ2Finite_element`;
59 tests). Revert: the tree holds uncommitted changes, so `git checkout` is NOT usable; the three mutated files
(`LogStrainNDMaterial.cpp`, `LadrunoNorSand.cpp`, `Domain.cpp`) were copied aside before the first mutant and copied back
after each one, and the SHA-256 of all ten changed SRC files was compared with the pre-mutation list after each mutant.

Baseline (unmutated, the build before any mutant): A 59 passed, B 113 passed (0 xfail, 0 xpass), L 59 passed.

### Result: 3 of 3 killed, no survivors (and the control behaves as predicted)

| # | Mutant (edit) | Zone A | Zone B | L | Named killer(s) and measured violation |
|---|---|---|---|---|---|
| MX1 | `LogStrainNDMaterial.cpp`: `if (prov != 0)` -> `if (false && prov != 0)` (route disabled, always the `inv(D0)` fallback) | 0 | 5 | 0 | `test_rigid_rotation_is_objective_with_pressure_dependent_shear_modulus[2.0]` stress 2.3e-3 (gate 1e-10) and `[50.0]` 3.6e-2; `test_provider_identity_committed_be_is_exp_two_eps_e[0.0, 2.0, 50.0]` b^e error 7.6e-3, 7.9e-3, 1.2e-2 (gate 1e-12) |
| MX2 | `LadrunoNorSand.cpp` `ladrunoGetElasticStrain`: shear returned as tensor, not doubled (`epsE(i) = sT.eps_e[i]`) | 0 | 6 | 0 | `test_rigid_rotation_is_objective_for_constant_shear_modulus` (second rotation, 1.6e-2), `..._with_pressure_dependent_shear_modulus[2.0]` 1.6e-2 and `[50.0]` 2.0e-2 (second rotation); `test_provider_identity_...[0.0, 2.0, 50.0]` b^e error 1.8e-3, 1.7e-3, 1.2e-3 |
| MX3 | `Domain.cpp`: the `ladrunoClearCommitRefusals()` before the element loop removed | 1 | 0 | 0 | `tests/test_ladruno_norsand.py::test_getcopy_of_a_latched_source_is_not_latched` ("step 1: the clone of a latched material started refused"); the only failure of 59 |
| MX0 (control) | `LogStrainNDMaterial.cpp` := the pre-G2 file from `git HEAD` (no provider route at all) | n/a | 5 | 0 | the same five as MX1, to the digit (it is the same behaviour); the fallback outputs are BIT-IDENTICAL to the current source, below |

Why Zone A does not see MX1 and MX2: the Zone A files drive the material through bricks and quads without the
LogStrain wrapper; the provider contract is a LogStrain-level fact and is gated in Zone B. Why Zone B does not see MX3:
the stale counter needs a `commitState()` outside any `Domain::commit()` (a getCopy source), which is what the Zone A
test does.

### The fallback is byte-identical to the pre-G2 LogStrain (MX0)

A script (`bitdump.py`, session scratch, not committed) drives `LogStrain` over `ElasticIsotropic` and over `LadrunoJ2`
(plastic: Mises(tau) 10.2 against the elastic line 17.5) through 25 coaxial steps, a z rotation or an oblique rotation, two
held steps and one general non-coaxial step, 30 steps per path, and records every step's Cauchy stress, Hencky strain and
the 24 x 24 element stiffness: 4 paths x 30 steps x 588 doubles = 70,560 doubles per build. Compared as raw 64-bit
patterns, the build with the pre-G2 `LogStrainNDMaterial.cpp` (MX0), the build with the route disabled (MX1), the one with
MX2 and the baseline / final builds are bit-identical on all four paths. (A recorded reference file is not committed: libm
`exp`/`log` have CPU-dispatched variants, so a bit pattern is not portable across machines; the committed gate is the
closed-form one above, and this measurement is the byte-identity evidence for this machine and toolchain.)

### Test-design finding from the first MX2 run

The first version of the rotation gates held F for one step after the rotation and compared the stress. Under MX2 only the
provider-identity test failed. Reason (a property of the plastic-inner protocol, not a bug): a held step feeds the inner
`eps_tr - eps_n`, which is zero whatever the committed b^e is, so it cannot see an error in b^e. A second, non-zero
rotation starts from the committed b^e and does; it was added (oblique axis (0.3, -1, 0.5), 0.35 rad) and MX2 is now also
killed by all three rotation tests. The v1 logs are kept in the session scratch (`g2mx/logs_v1`).

### Final state

After MX0 the source was restored and the pre-mutation SHA-256 list matched (all ten files, checked after each of the four
mutants and again at the end). A full `build.bat` (5 targets, no errors, all binaries re-stamped) was run from the
restored source. Re-run on it: B 113 passed; L 59 passed; A 59 passed; the sibling regressions (SaniSand
`test_ladruno_sanisand.py`, commit-refusal `test_ladrunoQuad_sanisand_implex_commit_refusal.py` 6,
`test_wp104_implex_refusals_wipe_reset.py` 3, `test_cdl_commit_solve_state.py` 7) pass; the Esmeralda oracle set
(`norsand_oracle/tests` + `kernel_parity`) 460 passed in 14:41; `mutate_kernel.sh` 15 of 15 mutants fail the parity gate.
`ci/check_quirk_patterns.py` 0 findings, `stamp_headers.py --check` current, `ci/check_classtags.py` OK,
`ci/check_manifest.py` OK. `build.bat clean` was not run (it wipes `build/` and `dist/`).

## Round 4 (2026-10-02): the G2-close refusal propagation and the commit-depth guard

Branch `wp/144-ladruno-norsand` at 2ce4e6d74 plus the uncommitted G2-close fix (`LogStrainNDMaterial.cpp/.h` and `LogStrain2D.cpp`
propagate `LADRUNO_MATERIAL_REFUSED`, inner commits first; `Domain.cpp` `LadrunoCommitDepthGuard`; StagedStrain warning) and the
G2 CLOSE section of `test_g2_logstrain.py` (GR0 to GR5). Nothing committed. Full `build.bat` (all 5 targets, 0 `error C`) per
build, from a fresh cmd.exe; mtimes of all five artifacts checked fresh each time. The stale `tests/opensees.pyd` (2026-10-02 02:14,
a copy of the pre-fix build; BUILD_GOTCHAS 4b) was moved aside for the whole campaign, runs used `PYTHONPATH` = a boot dir adding
`dist\bin` to the DLL search path, and after the final build a fresh copy was put back (byte-identical to `dist\bin\opensees.pyd`).
`ops.ladrunoBuild()` read inside pytest returns the HEAD hash, as it must; it cannot show the uncommitted fix, so freshness is by mtime.

Sets (all `--runslow`): A = `tests/test_ladruno_norsand.py` + `_element.py` + `_k1.py`; B = the whole `g2/` directory (run FROM
`g2/`); L = every LogStrain / finite-strain file (`test_logstrain*`, `test_finite_strain_*`, `test_finitestrain2d*`,
`test_ladrunoJ2_finite*`, `test_ladrunoJ2Finite_element`, `test_ladruno{Brick,Concrete3D,cst,lst,quad}_finite`,
`test_bezierTet10_finite`, `test_ladrunoRCFiniteStrain`); W = the commit-refusal files (`test_ladrunoQuad_sanisand_implex_commit_refusal`,
`test_wp104_implex_refusals_wipe_reset` and `_tcl`, `test_cdl_commit_solve_state`, `test_adr85_contact2d_t0_refusals`); S =
`tests/test_ladruno_sanisand.py`. Revert: the tree holds uncommitted changes, so each mutated file was saved aside first and copied
back; after EVERY mutant the SHA-256 list of all changed SRC files and `git diff -- SRC` equalled the pre-mutation copies.

Baseline (clean fixed source, full build): A 59 passed; B 120 passed; L 211 passed + 1 xfailed (the pre-existing strict v1-wrapper
kinematic-objectivity xfail in `test_ladrunoJ2_finite.py`); W 26 passed; S 17 passed + 1 skipped (the documented 2-rank placeholder).
GR0 to GR5: 7 passed (GR1 and the adversary script `refusal_under_logstrain.py` flipped: finite route now -3 at the trial, refusal
latched = 0, retry 0, the same as the -geom linear control).

### Result: 2 of 3 required mutants killed; MR2 SURVIVES (unreachable); two extra variants survive for the same reason class

| # | Mutant (edit) | A | B | L | W | Killer / verdict |
|---|---|---|---|---|---|---|
| MR0 (control) | `LogStrainNDMaterial.cpp`, `LogStrain2D.cpp`, `Domain.cpp` := git HEAD (the pre-fix source) | n/r | 2 fail | 0 | 0 | `test_GR1_finite_refused_trial_cuts_the_step_and_the_smaller_retry_converges` and `test_GR5_nested_commit_guard_is_in_place_...`: the gates detect the original defect |
| MR1 | `setTrialF` ignores the inner return code (`setTrialStrain(...)` result dropped, `if (false)`) | 0 | 1 fail | 0 | 0 | KILLED by `test_GR1_finite_...` only (behavioural: -4 + latch instead of -3) |
| MR2 | wrapper commits despite an inner refusal (both guards in `commitState`: the early `return rc` and the `trialRefused` branch, removed) | 0 | 0 | 0 | 0 | SURVIVOR: unreachable, see below |
| MR3 | depth guard removed (`Domain.cpp` := HEAD: unguarded clears, no `LadrunoCommitDepthGuard`) | 0 | 1 fail | 0 | 0 | KILLED by `test_GR5_...`, a STATIC source-text gate (not behavioural) |
| MR2A (extra) | only the early `if (rc == LADRUNO_MATERIAL_REFUSED) return rc;` removed (the `trialRefused` branch stays) | n/r | 0 | 0 | 0 | SURVIVOR; equivalent in practice: the second branch is reached under the same condition and returns the sentinel without advancing (only the refusal is declared twice) |
| MR3B (extra) | the guard exists but `outermost()` always returns true (the original bug, behaviourally) | n/r | 0 | 0 | 0 | SURVIVOR: GR5 greps for the `if (ladrunoDepth.outermost())` text, it does not execute the predicate |

(n/r = not run: the sets A and the mutated file are independent; MR0, MR2A and MR3B were run on B, L and W.)

WHY MR2, MR2A, MR3B CANNOT BE KILLED FROM THE TEST BED (the premise is gated by GR4 and GR5, which are static):
- MR2/MR2A: the wrapper's `commitState` refusal branches run only when a host COMMITS after a refused trial. Every in-tree caller of
  `setTrialF` (LadrunoBrick, LadrunoQuad, LadrunoCST, LadrunoCSTPair, LadrunoLST, BezierTet10, InitDefGrad) tests `< 0` and aborts
  `update()`, an aborted update never reaches `commitState()` (static, explicit and the CentralDifference family all test
  `updateDomain() < 0`), `FiniteStrainNDMaterial::setTrialStrain` is a hard error so `NDTest SetStrain` cannot drive a LogStrain, and
  `NDTest` has no `setTrialF` verb. GR4 therefore gates only the premise (every call site tests the return). The branch is defensive
  code with no test that executes it. To make it testable the product needs an `NDTest SetF` (or equivalent) verb; then the dynamic gate
  is: refuse a trial, commit anyway, assert -4 + latch, then a held step must still see `b^e_n` unchanged.
- MR3B: nesting needs an in-process `Subdomain` (a `_PARALLEL`-only class); the module has no subdomain command and `getNP() == 1`.
  GR5 checks the clears sit under `if (ladrunoDepth.outermost())` and that the guard is constructed first; it does not evaluate the
  predicate. A semantic regression of `outermost()` itself is not caught.

### Bit-identity of the non-refusing inners (MX0 repeated for this fix)

`bitdump.py` (the Round-3 script, copied; session scratch, not committed) drives `LogStrain` over `ElasticIsotropic` and over plastic
`LadrunoJ2` (25 coaxial steps, z or oblique rotation, two held steps, one general non-coaxial step): 4 paths x 30 steps x 588 doubles =
70,560 doubles (Cauchy stress, Hencky strain, 24 x 24 element stiffness). Compared as raw 64-bit patterns:
- MR0 (pre-fix wrapper + `Domain.cpp` from git HEAD) versus the fixed build: BIT-IDENTICAL on all four paths (70,560 of 70,560).
- the Round-3 final dump versus the fixed build: BIT-IDENTICAL on all four paths.
- the final clean rebuild versus MR0: BIT-IDENTICAL on all four paths.
A non-refusing inner is unchanged by the fix to the last bit on this machine and toolchain (not a portable claim: libm `exp`/`log`).

### Esmeralda (independent of the Windows build; the fix touches no kernel)

`norsand_oracle/tests` (G1) + `kernel_parity`, synced from the worktree with `SRC/material/nD/LadrunoNorSand{Kernel.h,.h,.cpp}`,
`nohup` + poll: 460 passed in 14:20 (the same count as Round 3). `mutate_kernel.sh` was not re-run (no kernel source changed).

### Final state

Source restored after each mutant (SHA-256 of the changed SRC files and `git diff -- SRC` identical to the pre-mutation copies, checked after
every one). Final full `build.bat` from the clean fixed source: 5 targets, exit 0, 0 `error C`, all five mtimes fresh; A 59, B 120,
L 211 + 1 xfail, W 26, S 17 + 1 skip, all green. `build.bat clean` was not run.

## Round 3b (2026-10-03): the shell mutants of the HAR energy, the p' floor and the unified pi_i0 (NOT RUN: written for the Mutate step)

New gates (all need the round-3b build; Zone A = `tests/`, Zone B = this folder; `--runslow` for the slow one):
- A: `tests/test_ladruno_norsand_har_floor.py` (52 tests): echo, every parser refusal with its code, the HAR closed forms (K1.1h, K1.11, DM04, the
  gate-table eigenvalues), the floor (K1.12, K1.13, K1.14, responses, counters, tangent, initial projection, K1.15), the cube past the domain edge,
  the element FD tangent under HAR, the database round trip of the floor block.
- B: `test_g2_shell_parity.py` (round 3b: floor counters and responses compared per step with O2 on the kernel_parity FLOOR_* / HAR_* / PI0_* paths,
  delta : C = 0 at floored steps, the p_min = 0 HAR refusal, the pi_i0 / initial projection test), `test_g2_logstrain.py` (HAR paths through the
  wrapper at 1e-10, HAR rigid-rotation objectivity, HAR provider identity), `test_g2_floor_sensitivity.py` (slow: F vs F/2 on the strip BVP).

| # | Mutant (edit of `SRC/material/nD/LadrunoNorSand.cpp` / `LadrunoNorSandKernel.h`) | Killers (expected) |
|---|---|---|
| S1 | `-energy HAR` parsed, the BA06 law built (HAR -> BA06) | shell parity HAR_* paths and `test_har_paths_run_the_har_law_not_ba06`; Zone A K1.1h / K1.11 / gate table; NOT any FD test (round 3b, A4) |
| S2 | BA06 constants accepted under HAR (code 201 removed) | `test_har_refuses_each_ba06_constant_code_201` (5 flags) |
| S3 | HAR constants accepted under BA06 (code 202 removed), or `-p_a` added to that list | `test_ba06_refuses_each_har_constant_code_202` (6 flags + the p_a positive controls) |
| S4 | k / g / n / p_a range checks removed or off by one (codes 21-24), `-pmin < 0` accepted (25) | `test_har_kernel_range_refusals_codes_21_to_25`, `test_pmin_negative_refused_code_25_under_ba06_too` |
| S5 | (S.56) refusal absent / ungated (planar, none refused) / with the inverted W_ramp / factor 10 changed | `test_smooth_cap_scan_gate_code_26_only_for_smooth` |
| S6 | trial floor skipped (post only): the out-of-domain HAR trial refuses (M-F5) | `test_k1_13_har_floor_in_and_out_of_the_domain_is_counted_not_refused[1.1e-4]`, parity FLOOR_HAR_K113_out_of_domain |
| S7 | floor not applied (pass-through) or applied at p_min = 0 | `test_k1_12_*`, `test_floor_off_pmin_zero_*`, parity FLOOR_BA06_* |
| S8 | `floor` / `floorEnergy` / `floorInit` slots reordered or counters not accumulated (M-F2) | Zone A K1.12 / K1.13 / initial-state tests, parity `floor_mismatch` |
| S9 | E_f sign / out-of-domain bound (S.52) | Zone A K1.13 `floorEnergy[0]` (Psi difference in the domain, W_f outside), parity E_f at 1e-10 |
| S10 | floored tangent with a bulk stiffness (regularisation) or the (S.51a) eps' term dropped (M-F3b) | Zone A K1.12 / K1.13 tangents, K1.14 (delta : C, a^e Phi), parity `dC_floor` and the tangent |
| S11 | q kept instead of eps_s under HAR (M-F6), pi_i or v altered by the projection (M-F7) | Zone A K1.14 (eps_s unchanged, q_f), K1.12 (pi_i, v) |
| S12 | `-pi0_auto` = the pre-round-3 apex rule, or the rule applied before the floor | `test_pi0_auto_is_the_unified_rule_k1_15`, parity `test_pi0_auto_and_the_floored_initial_state_through_the_shell_equal_o2` |
| S13 | initial state not projected / not counted | `test_initial_state_is_projected_counted_and_the_stress_replaced_ba06_and_har` |
| S14 | sendSelf / recvSelf drops the floor block (counters reset on restore) | `test_database_roundtrip_carries_the_floor_counters_and_energy` |
| S15 | DM04 mapping: (1+nu) <-> (1-nu), f(e) wrong, e_ref ignored | `test_dm04_mapping_gives_the_printed_g_and_k_and_the_same_response` |
| S16 | provider route off for HAR (LogStrain falls back to inv(D0)) | `test_rigid_rotation_is_objective_under_the_har_energy`, `test_provider_identity_under_har_*` |
| S17 | floor inside the local Newton (M-F9) or hidden stiffness: the limit load moves > 2 % | `test_g2_floor_sensitivity.py` (slow), parity FLOOR_HAR_FPf_x4 |
