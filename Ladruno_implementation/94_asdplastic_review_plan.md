---
title: ADR 94 — ASDPlasticMaterial3D implementation review
project: Ladruno
status: complete (2026-09-07) — verdict #804, fix wave #806/#809/#815 merged
priority: high
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR 94 — ASDPlasticMaterial3D implementation review (plan)

> Number 94 allocated 2026-09-06 (93 is held by the open draft #801). Bare "ADR 94"
> is legal — no other doc carries it. Verify with `ls Ladruno_implementation/94_*`
> before citing.

## 1. What

A warrant-first (ADR-87) review of the whole `SRC/material/nD/ASDPlasticMaterial3D/`
tree: jaabell's template framework (`ASDPlasticMaterial3D.h`, 4,165 lines; 7
integrators, 5 tangent operators, the `utuple_storage` state machine, the OPS
parser) **plus** everything the fork layered on it (Hoek–Brown, StiffSoil
shear/cap, MohrCoulombTensionCutoff, `strict_convergence`, the sentinel
refusal, the ResponseType labels). 58 files, ~12k lines, 46 registered
YF×PF×EL×IV specializations, class tag 10000.

**In scope:** correctness of the state machine (trial/commit/revert/copy),
the return maps and their tangents, component derivatives, the parser
contract, the host-element / parallel contract, test coverage, and drift
against jaabell's `ASDP` branch.

**Out of scope (separate WPs if the review warrants them):** rewriting an
integrator (e.g. a true closest-point BE, see H6), **anything OpenSeesSP/MP**
(owner decision D3, 2026-09-06: parallel is not a target for ASDP; H3 is
recorded, not worked), plane-strain support,
and the upstream port of framework fixes (goes through
[[upstream_pr_campaign]]).

**Deliverable of this ADR:** a severity-ranked verdict
(`reviews/adr94_verdict.md`) in which every finding is either CONFIRMED with
a pinned Zone-A test or REFUTED with the reproducer that failed to show it,
plus a fix plan split into fork PRs and jaabell-bound PRs. Fixes ship in
their own `wp/` branches, not in this one.

## 2. Why

- ASDP is the fork's geomaterial workhorse: Cerro Lindo (ADR-0005 M3/M5) runs on
  it, ADR-84 built MCTC on it, the finite-strain plan (ADR-17 D3) uses it as
  the *independent return-map oracle*, and Hoek–Brown + StiffSoil were added
  for rock/stiff-clay work.
- Only one slice has ever been reviewed: the MCTC/`Backward_Euler` path
  (ADR-84 P0–P3, ADR-86b review). ADR-84 already found **three** upstream
  silent-accept defects in that slice alone; the other six integrators, the
  numerical tangents, and every component except MC/MCTC have never been
  read critically.
- A first read-through for this plan surfaced a further set of candidate
  defects (§4), one of which — a class-wide `static` tangent — would make
  every multi-Gauss-point element receive the wrong consistent tangent and
  is invisible to every existing test.
- Our Hoek–Brown copy is one commit *behind* jaabell's `ASDP` branch (§4 H10);
  the upstream campaign row 1.5 assumed the opposite.
- The ADR-84 §6 phasing left P2(b–e) and P4 open; this review decides them.

## 3. Where — inventory (verified 2026-09-06, HEAD `76668332f`)

| Layer | What exists | Registered / reachable | Tested |
|---|---|---|---|
| Elasticity | `LinearIsotropic3D_EL`, `DuncanChang_EL`, `StiffSoil_EL` | LinearIsotropic, StiffSoil; **DuncanChang commented out** of the generator | LinearIsotropic only |
| Yield functions | VonMises, DruckerPrager, MohrCoulomb, RoundedMohrCoulomb, TensionCutoff, HoekBrown, StiffSoilShear, StiffSoilCap, MohrCoulombTensionCutoff | all but **TensionCutoff (commented out)**; RoundedMC not in the generator list | MC, MCTC only |
| Flow rules | VonMises, DruckerPrager, MohrCoulomb, ConstantDilatancy, HoekBrown, StiffSoilShear/Cap, MCTC | all but **ConstantDilatancy (commented out)** | MC, MCTC only |
| Hardening | Linear (scalar/tensor), Null (scalar/tensor), ArmstrongFrederick; Exponential commented out; StiffSoil-specific | — | none directly |
| Integrators | `Forward_Euler`, `Forward_Euler_Subincrement`, `Backward_Euler` (default), `Backward_Euler_LineSearch`, `Modified_Euler_Error_Control`, `Runge_Kutta_45_Error_Control`, `Runge_Kutta_45_Error_Control_old` | all 7 parse-reachable; `Full_Backward_Euler`, `*_Crisfield`, `Multistep_*` exist in the enum/setter but have no dispatch case and no parser token | **BE only** |
| Tangents | Elastic, Continuum, Secant (default), Numerical_Algorithmic_First/SecondOrder | all | Elastic/Continuum/Secant (ADR-84 P3, MCTC only) |
| Parser options | `f_absolute_tol`, `stress_absolute_tol`, `n_max_iterations`, `strict_convergence`, `rk45_dT_min`, `rk45_niter_max`, `return_to_yield_surface`, `integration_method`, `tangent_type` | — | `strict_convergence` |
| Standalone C++ tests | `test_HoekBrown.cpp`, `PlasticFlowDirections/test_pf.cpp`, `YieldFunctions/test_00_VonMises_YF.cpp` | **not built by any CMakeLists** | never run in CI |
| Zone-A pytest | `test_asdplastic_mctc.py`, `test_asdplastic_response_tags.py`, `test_adr84_p2a_strict_convergence.py`, `test_adr84_p3_confined_corner.py` | — | MC/MCTC/VM on BE; hosts `stdBrick`, `TenNodeTetrahedron`, `tri6n` |
| Fork markers | 41 `// Ladruno` lines in 9 files; `LEDGER_vanilla_files.md` rows 300, 337 (**UNMARKED** HB/StiffSoil block), 490–500 | — | — |
| Upstream | `jaabell/master` (2026-05-15) has nothing we lack; `jaabell/ASDP` has `60d9b9b23` "more changes to HB" (2026-05-12) that we lack | — | — |

Docs that already bear on this tree: [[84_ladruno_mc_tension_cutoff_adr]] (§6, §6b, §6c, §9),
[[reviews/adr86b_verdict]], `LEDGER_quirks.md` §"ASDPlasticMaterial3D" (five entries: NaN
heap, FullGeneral N=0, plastic-path tuning, BE silent accept, tangent-override, zero-free-DOF
driver), [[17_finite_strain_validation_plan]] §D3, [[upstream_pr_campaign]] rows 0.4 and 1.5.

## 4. Candidate findings from the read-through (the H-list)

These come from reading, not measuring. Each one is a hypothesis the review must
CONFIRM (pin with a test) or REFUTE (record the reproducer that failed). Severity
is provisional. Line numbers are for `ASDPlasticMaterial3D.h` at `76668332f`
unless a file is named.

| # | Claim | Where | Why it matters | How to test |
|---|---|---|---|---|
| **H1** | `Stiffness`, `dsigma`, `depsilon_elpl`, `intersection_stress/strain` are **class-wide `static` members** (4116–4120, defined 4150+), shared by every instance of a specialization. `getTangent()` (698) copies the static. Host elements call `setTrialStrain` for all GPs in `update()` and `getTangent()` for all GPs in a *separate* loop (`Brick.cpp:1069` vs `:1201`; `LadrunoBrick.cpp:1340`; `TenNodeTetrahedron.cpp:1355`), so **every GP receives the last GP's tangent**. | 698, 4116–4120 | Wrong consistent tangent in any element with a strain gradient; degrades global Newton; invisible in every existing test because all drivers are single-element homogeneous strain. Also makes ASDP permanently unsafe for ADR-75b threaded assembly. | Two-element (or one bent element) model: record per-GP `getTangent()` vs an FD tangent from each GP's own response; `ops.testIter()` count vs a per-instance-member build. Oracle: numpy DP with two different strain states. |
| **H2** | `getClassType()` returns `name.c_str()` of a **local `std::string`** — dangling pointer, UB. | 196–201 | Works by SSO luck; any recorder/`Print` path that reads it after the call is UB. | Trivial; fix is one line. Confirm by reading. |
| **H3** | `sendSelf`/`recvSelf` are stubs that print "not implemented" **and return 0 (success)**. The integration options live in per-tag `static std::map`s that no channel ever ships. | 1223–1238, 4100–4111 | OpenSeesMP/SP ship an unconfigured material silently. | **Out of scope (D3, 2026-09-06).** Recorded as a known limitation in the verdict; no np2 run, no fix. Revisit only if an MP run ever needs ASDP. |
| **H4** | `revertToLastCommit()` body is commented out (no-op); `revertToStart()` returns −1. | 752–778 | After a cutback (adaptive stepping, `DisplacementControl` retry, `LadrunoBeginAugment` sweep, arc-length) `getStress()` reports the *failed trial*; `revertToStart` breaks `reset`/`revertToStart` analyses. | Force a Newton failure + cutback and read `getStress` between revert and the next trial; call `ops.reset()` after a run. |
| **H5** | The "f decreased ⇒ elastic step" exit that ADR-84 P2a gated in BE exists **unguarded at seven more sites**: `compute_local_stress` (508), FE (1423), FE_sub (1599), BE_LS (2406), RK45_old (2667), ME (3088), RK45 (3435). `strict_convergence` only reaches BE. | listed | Same silent-perpetuation of `f>0` from an inadmissible commit on every non-default integrator and inside the numerical tangents. | Re-run `test_adr84_p2a_strict_convergence::test_flag_on_refuses…` parametrized over `integration_method`. |
| **H6** | `Backward_Euler` is a **cutting-plane** (Ortiz–Simo 1986) algorithm, not a closest-point projection: `n`, `m`, `H` are re-evaluated at the current iterate and the stress is corrected incrementally (2244–2320); IVs are advanced with `h` at the *updated* stress. The "Continuum" tangent (453–465) is the continuum operator, not the consistent tangent of this map, and `Numerical_Algorithmic_*` differentiate a *third* map (`compute_local_stress`, ADR-84 P4). | 2034–2360, 438–482 | First-order accuracy; no material-level quadratic convergence; no consistent tangent exists for any hardening YF. Not wrong, but mis-named and undocumented, and it caps what any tangent option can deliver. | Convergence-order study on VM+linear hardening (closed form) and DP: error vs Δε at 1st vs 2nd order; Newton iteration counts with free DOFs. |
| **H7** | BE's `dLambda + deltaLambda < 0` fallback (≈2300) returns 0 with `Stiffness = Eelastic` but leaves `TrialStress` at the *partially corrected* iterate, not the elastic predictor. BE's Continuum tangent reads static `depsilon_elpl`, which BE never sets (stale from another GP or integrator). | ≈2300, 456 | Committed state is neither elastic nor on the surface; non-associated PFs that use `depsilon` get a foreign strain increment. | Construct a path that triggers the branch (print), compare commit to predictor; check `depsilon_elpl` provenance under Continuum. |
| **H8** | `Backward_Euler_LineSearch`: hardcoded `max_iter = 30`, `tol_rel = 1e-8`; the "line search" tests a **linear prediction** `Phi + dPhi·dl`, never re-evaluates `Phi`, so with a Newton direction it always accepts α=1; ignores `strict_convergence`; the split loop reuses `Eelastic` at commit for all substeps. `Forward_Euler_Subincrement` reuses `n_max_iterations` as the **substep count**. | 2359–2600, 1622 | A misnamed algorithm the user might select believing it is robust; one option meaning two things. | Read + a stiff MC path comparing BE vs BE_LS iteration/refusal behaviour. |
| **H9** | ME and RK45 drift checks are **empty `if` blocks** (3246, 3767); the only correction is `return_to_yield_surface`. `Runge_Kutta_45_Error_Control_old` is still parse-reachable. Dead enum values (`Full_Backward_Euler`, `*_Crisfield`, `Multistep_*`) are accepted by the setter but have no dispatch case → runtime −1 if ever selected programmatically. | listed | Explicit integrators can commit `f ≫ tol` silently; a dead API surface. | Drift measurement on DP triaxial with ME/RK45 at coarse steps, `return_to_yield_surface Disabled`. |
| **H10** | **Hoek–Brown drift.** `jaabell/ASDP` `60d9b9b23` (2026-05-12) rewrote `HoekBrown_YF` as a composite `max(f_shear, f_tension)` with `pow(max(arg,0), a)`; ours (`e65e89203`, 2026-05-10) is the previous if/else version. The fork is one HB commit behind, and [[upstream_pr_campaign]] row 1.5 assumed we were ahead. `DruckerPrager_YF::CHECK_APEX_REGION` returns `false` with an `// Implement!!!` comment while declaring `yf_has_apex = true`. | `YieldFunctions/HoekBrown_YF.h`; `DruckerPrager_YF.h:120–133` | HB tension behaviour differs between the two trees; DP apex handling is a stub. | `git diff jaabell/ASDP -- …/HoekBrown_YF.h`; decide port direction with José; DP apex probe. |
| **H11** | **Coverage.** HB, StiffSoil (3 combos), DP, VM (beyond a P2a inertness check), RoundedMC, DuncanChang: no Zone-A test. FE, FE_sub, BE_LS, ME, RK45 and both `Numerical_Algorithmic_*` tangents: no test at all. The three standalone `.cpp` tests are not built. | `tests/`, CMake | The review cannot certify what it cannot run. | R4 matrix (§5). |
| **H12** | **Diagnostics/perf.** 161 `cout` vs 9 `opserr`; `commitState` prints a line per material per step whenever `GLOBAL_INT_max_iter > 0`; `INT_OPT_*[ASDP_TAG]` `std::map::operator[]` lookups on every GP call (default-inserting). | 727–748, throughout | Invisible in OpenSeesPy captured output, interleaved under MPI, stdout flood on real models; measurable per-GP overhead. | Count lines on a 1k-element run; profile one `setTrialStrain`. |
| **H13** | **Parser contract.** Unknown `Begin_Integration_Options` token: no `else` after the if-chain, so a typo (`strict_convergance 1`) is silently skipped and its value is consumed as the next name. Unknown model parameter: verify `setParameterByName_impl` fallthrough. Unset parameters default to 0 (ADR-84 P2(e) still deferred). `getCopy(type)` refuses everything but 3D. | `OPS_AllASDPlasticMaterial3Ds.cpp:358–470`, `utuple_storage.h:189` | Silent misconfiguration is the most expensive failure mode this fork has recorded (ADR-84 §9.1). | Typo decks; assert loud failure. |
| **H14** | `first_step`/`InitialP0`/`stress_set_externally` handshake: `getCopy` does not copy `first_step`; `InitialP0` is written into `CommitStress` on the first `setTrialStrain` of *every* copy; interaction with `setParameter stress`/`updateParameter` ordering unverified. | 229–246, 780–830 | Initial-stress workflows (geostatic, staged) depend on it. | `InitialP0` + `setParameter stress` before/after first analyze; compare. |
| **H15** | BE and BE_LS use `E(σ_commit)` for the whole return; ME/RK45 re-evaluate `E` per stage. For stress-dependent elasticity (`StiffSoil_EL`, `DuncanChang_EL`) the implicit maps are inconsistent with the elastic law. | 2055, 2372 | StiffSoil on the default integrator is a different material than StiffSoil on RK45. | StiffSoil triaxial BE vs RK45 at fine steps; document or fix. |

> [!question] H1 is the one to run first. If CONFIRMED it re-reads several Cerro Lindo
> "material stall" observations (ADR-84 §9.1 measured the corner tangent under a
> single-element driver — which cannot see H1 at all).

## 5. How — review method (phases)

Every phase produces a committed artefact; nothing is "known" until it is in a file.

**R0 — Freeze and baseline (½ session).**
- Build `OpenSees OpenSeesPy` in this worktree; open every probe with `ops.ladrunoBuild()`.
- Dump `nDMaterial ASDPlasticMaterial3D 999 list` → `_adr94_inventory.md` (the 46 registered strings, which are reachable, which are dead).
- Run the four existing ASDP test files; record pass counts and wall time as the baseline.
- Confirm the `LEDGER_vanilla_files.md` row-337 UNMARKED block and add the missing `// Ladruno` markers (bookkeeping, no behaviour change).

**R1 — Verify the H-list (1 session).**
- One reproducer per H, in `tests/test_adr94_hlist.py` (Zone-A, `t0m`), each test named after its H and asserting the *defect* (so a fix flips it, per the sentinel rule in `LEDGER_quirks.md`).
- H1 and H6 get numpy oracles under `Ladruno_implementation/adr94_oracle/` (`python3.12`): two-state DP consistent tangent (H1); VM+linear-hardening closed-form vs cutting-plane order (H6).
- Output: `_adr94_hlist_verdicts.md` — CONFIRMED / REFUTED / severity / blast radius per H.

**R2 — Red/blue review of the framework core (1 session).**
- Same format as `_adr92_p1_redblue/`: lanes `cpp` (state machine, statics, copy/revert/serialize, parser), `numerics` (integrators, tangents, apex/corner handling, tolerance semantics), `tests_process` (coverage, ledgers, markers, upstream provenance).
- Red argues "unsafe to rely on"; blue argues "adequate as shipped, opt-in fixes only". The verdict doc adjudicates with R1 evidence, not opinion.

**R3 — Component audit (1–2 sessions).**
- Per YF/PF: `f`, `∂f/∂σ`, `hardening` vs finite differences on a stress cloud that includes the Lode edges (θ = ±30°), the apex, hydrostatic states, and `J2 → 0`; non-finite guards; `yf_has_apex` honesty (DP stub).
- Per hardening law: sign/units of `h`, tensor-vs-scalar policy, AF saturation.
- Per EL: SPD of `E(σ)` over the admissible range (StiffSoil, DuncanChang).
- Hoek–Brown: reconcile with `jaabell/ASDP` `60d9b9b23` (port his composite or argue ours), rebuild `test_HoekBrown.cpp` into a CMake `ctest` or a pytest twin.
- Delivery: a component table in the verdict + `tests/test_adr94_components.py` (material-point drivers with **free z-DOFs** so tangents are observable — the ADR-84 lesson).

**R4 — Integrator × tangent matrix (1–2 sessions).**
- {VM, DP, MC, MCTC, HB} × {FE, FE_sub, BE, BE_LS, ME, RK45} × {Elastic, Continuum, Secant, NumAlg1, NumAlg2} on three paths: triaxial compression to the surface, simple shear, tension to a cutoff.
- Metrics: `|f|` at every commit, drift vs Δε (convergence order), `ops.testIter()` with free DOFs, refusal behaviour under `strict_convergence`, wall time per step.
- Output: the **recommended-configuration table** for the user guide (which integrator/tangent per YF, and which combinations are refused or unsupported). This table *is* the verification manifest ADR-87 asks for.

**R5 — Host-element contract (½ session, gated on H4).**
- `LADRUNO_MATERIAL_REFUSED` propagation from every integrator, not only BE (today the sentinel exists only at two BE sites).
- Revert semantics under `LadrunoBeginAugment`/`EndAugment` and adaptive stepping.
- No MP work (D3). H3 gets one sentence in the verdict's scope section.

**R6 — Verdict and fix plan (½ session).**
- `reviews/adr94_verdict.md` in the `adr86b_verdict.md` shape: default-inertness, fail-loud, conventions, warrant, scope, fix list, test results.
- Fix list split into: (a) fork WPs (opt-in or default-inert unless proven safe), (b) jaabell-bound framework fixes (H1, H2, H4, H5, H13 — no AI traces, fresh-branch ports per the campaign rules), (c) documentation-only (H6, H8, H15 if not fixed).
- Ledger rows, README active-plan line, and `banner_features.txt` untouched (no new feature ships from this ADR).

## 6. Acceptance (for the review itself)

- Every H in §4 has a verdict backed by a committed reproducer or oracle; "read the code" is not a verdict.
- Every CONFIRMED defect has a pinned Zone-A test that FAILS on `76668332f` and names the H.
- The R4 table exists and every cell is measured (or marked "refused"/"crashed" with the log line).
- Default behaviour of every registered material on `Backward_Euler` + `Secant` is byte-identical before and after any bookkeeping change made under this ADR (marker comments, ledger rows).
- Zone-A green on the review branch; the four pre-existing ASDP test files unchanged in pass count.

## 7. Branch / PR

- Rename this worktree's branch to `wp/94-asdplastic-review` (ADR-87 D9) and open a **draft** PR on day one; commit R0–R6 artefacts continuously.
- Fixes do **not** land on this branch. Each fix gets `wp/94x-<slug>` with its own mutation gate, so the review PR stays a review.
- Check `git worktree list` and `gh pr list` before starting a fix WP (the #753/#754 duplicate-lane lesson).

## 8. Risks and traps (already recorded, restated so nobody rediscovers them)

- `OPS_AllASDPlasticMaterial3Ds.cpp` is a `/bigobj` TU that instantiates all 46 specializations: every header change is a long rebuild — batch C++ edits, one build per batch.
- `system("FullGeneral")` hard-crashes fully-prescribed drivers (N = 0). Use `UmfPack`.
- A fully-prescribed material-point driver cannot see a wrong tangent. Leave DOFs free and gate on `ops.testIter()`.
- The NaN-heap trap: a standalone probe passes, pytest's churned heap fails. Run the battery, not the probe.
- `stdBrick`/`BrickUP`/`QuadUP` swallow material refusals; use `LadrunoBrick` or `TenNodeTetrahedron` for any refusal gate.
- Framework fixes that change default behaviour for VM/DP users are upstream-facing; coordinate with José before shipping them as fork defaults (ADR-84 §3 "bit-identical strategy").
- Number collision: 94 is claimed by this file; a second `94_*` must not appear.

## 9. Decisions needed from the owner

> [!question] **D1 — Scope.** Review + verdict + fix *plan* (this ADR) with fixes as separate WPs. Confirm, or fold the top-severity fixes (H1, H2) into the review branch.

> [!question] **D2 — H6.** If the cutting-plane classification is confirmed: document and leave (recommended for this ADR) or open ADR 95 for a true closest-point BE with a consistent tangent for VM/DP/MC.

> **D3 — DECIDED 2026-09-06 (owner): SP/MP is not a target.** H3 is recorded as a limitation; no refusal, no implementation, no np2 gate.

> [!question] **D4 — H10.** Port jaabell's HB composite (his tree is newer) or send him ours. Needs a direct exchange; do not guess.

## 10. Orchestration — who runs what

The owner's session (Fable) orchestrates with the `Agent` tool: phases run in
order, agents inside a phase run in parallel, all in the background. The
orchestrator never reads a transcript — only the file each agent was told to
write. Twelve agents total; one build.

| Phase | Agent | Model · effort | Lane / inputs | Writes | Returns |
|---|---|---|---|---|---|
| R0 | 1 × general-purpose | sonnet · medium | build `OpenSees OpenSeesPy` (WMI-launched `build.bat`, foreground wait via a log `Monitor`, never background its own build); `list` dump; run the 4 ASDP test files; add row-337 markers | `_adr94_inventory.md`, marker commits | ≤20 lines: build hash, pyd mtime, pass counts |
| R1-A | 1 × general-purpose | **opus · high** | H1, H6, H7, H8 — the numerics lane; gets the H rows verbatim + line ranges, not "read the header" | `adr94_oracle/*.py`, `tests/test_adr94_hlist.py` (its H's), its section of `_adr94_hlist_verdicts.md` | ≤30 lines: verdict per H + the one number that decides each |
| R1-B | 1 × general-purpose | sonnet · medium | H2, H4, H5, H9, H12, H13, H14, H15 — mechanical reproducers (pytest only, no C++) | same test file (its H's), its verdict section | ≤30 lines |
| R1-C | 1 × general-purpose | sonnet · medium | H10 (`git diff jaabell/ASDP`, DP apex probe); H3 gets a one-line "out of scope" entry, no run | verdict sections + `_adr94_hb_drift.md` | ≤20 lines |
| R2 | 3 red + 1 blue, general-purpose | red-numerics **opus · high**; red-cpp, red-tests, blue **sonnet · medium** | R1 verdict file as the only input; blue covers all three lanes in one agent | `_adr94_redblue/{red1_cpp,red2_numerics,red3_tests,blue}.md` | ≤15 lines each |
| R3 | 2 × general-purpose | harness **sonnet · medium**; apex/Lode/HB reconciliation **opus · high** | generic free-DOF FD material-point harness, then every registered YF/PF/EL; opus decides HB port direction with the diff | `tests/test_adr94_components.py`, component table section | ≤30 lines |
| R4 | 1 × general-purpose | sonnet · medium | parametrized matrix (5 YF × 6 integrators × 5 tangents × 3 paths); long-running, background | `tests/test_adr94_matrix.py`, `_adr94_matrix.md` (the table) | the table only |
| R5 | 1 × general-purpose | sonnet · medium | sentinel propagation per integrator; revert under augment/adaptive stepping | verdict section | ≤20 lines |
| R6 | orchestrator | fable · high | reads the six `_adr94_*` files, not the transcripts | `reviews/adr94_verdict.md`, ledger rows, README line | — |

Model rationale: opus where mechanics judgment decides the answer (H1/H6/H7/H8,
the numerics red lane, apex/Lode/HB); sonnet where a harness or a matrix is
written once and run many times; haiku nowhere (Zone-A tests need care).

Token discipline (binding for every agent prompt):

- One build, in R0. R1–R5 edit **no C++**. If a finding needs a source change to
  *prove* it (H1's per-instance member), one agent builds it once on a scratch
  branch at the end of R1 — not before, not twice.
- Every agent writes its output to the named file and returns a capped summary;
  the orchestrator forbids file dumps in the prompt. Agents get file:line ranges
  and the H rows, never the whole 4k-line header.
- No full Zone-A sweep; run only the named test files (memory: sessions ran too
  long doing the 27-minute sweep by default).
- `run_in_background: true` everywhere; wait on notifications, never poll. A build
  agent that backgrounds its own build is never re-woken (recorded trap) — R0's
  prompt says so explicitly.
- Phase gates: R2 waits for R1's verdict file; R5 waits for the H4 verdict; R6
  waits for all. Nothing else blocks, so R3 and R4 start as soon as R0's build
  exists.

## Implementation log

- **2026-09-06** — plan drafted; D3 decided by owner (SP/MP out of scope). Merged as #802 (`220fbac6f`).
- **2026-09-07** — R0–R6 executed on `wp/94-asdplastic-review-exec` (#804), 13 agents, one build
  (`52314165a`). R0 inventory: 46 specializations, baseline 22 passed. R1: H1/H7/H8 (opus),
  H2/H4/H5/H9/H12–15 (sonnet), H10 (sonnet). R2: three red lanes + one blue. R3a: standalone
  g++ FD harness over every registered component — worked without linking OpenSees; found the
  shear-slot convention split (B5) after a first "VM derivative bug" reading was withdrawn.
  R3b (apex/HB opus lane) dropped as redundant with R1-C + red-numerics. R4: 450-cell matrix;
  DP/HB rows not warrant-grade (harness sign convention disputed, see verdict §2). R5: sentinel
  at 2 of 15 sites; hosts drop bare −1. R6 verdict: [[reviews/adr94_verdict]] — 15/16
  confirmed, H6 accuracy REFUTED on VM, H1 downgraded from "results" to "cost + threading"
  after blue's two-cube measurement. Full battery 69 passed / 4 skipped / 0 failed in 34 s.
- **Surprises:** `capfd` cannot see the `.pyd`'s `cout` (child-process capture needed);
  `printA('-ret')` is empty except under `FullGeneral`; the R0 agent parked on its own build
  monitor (the recorded trap) but was woken by Monitor events; `stdBrick` swallows the sentinel
  too, so ADR-84 P2a is a no-op on the default hex; `meanStress()` is tension-positive and the
  DP comment is wrong.
- **Row 337 marker debt** closed in the same PR (comment-only markers on the HB/StiffSoil
  integration hunks).
- **2026-09-07, fix wave** (owner: "we are allowed to change this vanilla material" — the
  jaabell-bound column of the verdict collapsed into fork WPs; "merge, continue the orchestration").
  Four PRs, all merged the same day: **#806** `wp/94d-hb-port` (jaabell's composite Hoek–Brown,
  tension plateau 587 → 245.0 kPa = textbook σt, compression bit-identical); **#809**
  `wp/94a-fail-loud` (sentinel at all 15 failure sites, strict mode on every integrator, loud parser
  + required-parameter check, Eigen `setZero` at 6 sites, `Backward_Euler_LineSearch` and
  `Runge_Kutta_45_Error_Control_old` refused; 78 passed); **#815** `wp/94c-numerics` carrying
  `wp/94b-statics` (#813 closed as superseded after a syntax slip in a test edit): per-instance
  `Stiffness`/buffers, real `revertToLastCommit`/`revertToStart`, `getCopy` copies `first_step`
  (two-cube iteration penalty +62 % → −12 %); ONE shear-slot convention (Voigt) at every YF/PF and
  every contraction site — VonMises simple shear vs closed-form radial return **1.5e-13**, VM/DP FD
  gradient error 3.5e-1 / 9.7e-1 → 1e-8; live apex return with DP apex methods (hydrostatic tension
  40/40 finite); opt-in `f_relative_tol` with per-YF `strength_scale()`; 89 passed.
- **Results changes stated plainly:** DP paths change by design (2.7e-2 rel on the review deck);
  Hoek–Brown moves 2.5e-3 rel between two admissible states (cutting-plane path dependence, verdict
  M3); VM/MC/MCTC shear-free paths differ only at FP re-association level (≤ 2e-9 rel).
- **Left open after the wave:** `stdBrick` still swallows every material return code (by design,
  pinned); no consistent tangent exists for the cutting-plane map (D2, a closest-point rewrite would
  be ADR 95); `StiffSoilShear` NaNs on step 1 (M9, untested combo); `RK45_old` final NaN guard still
  calls `exit(-1)` (integrator is parser-refused); `tests/test_adr94_matrix.py` regenerates the
  tracked `_adr94_matrix.md` on every run and its HB oracle predates the composite port; banner text
  amended after the last build (cosmetic).
- **Orchestration lessons** (also in `LEDGER_quirks.md`): MSVC green ≠ GCC green (temporaries into
  non-const refs); cross-platform float pins need global-tolerance-size bounds (≥ 1e-6 rel), 1e-9
  failed twice; a bit-identity gate cannot certify a static-state/tangent fix; agents that stage
  edits as scripts survive session restarts (three were orphaned), builds launched via WMI survive
  too, `Start-Process` ones do not; run a test file once before pushing a "trivial" edit to it.

## See also

[[reviews/adr94_verdict]] · [[84_ladruno_mc_tension_cutoff_adr]] · [[reviews/adr86b_verdict]] ·
[[17_finite_strain_validation_plan]] · [[upstream_pr_campaign]] · [[87_ladruno_depth_with_width_adr]] ·
[[LEDGER_vanilla_files]] · [[LEDGER_quirks]] (the ASDPlasticMaterial3D entries, now twelve) ·
[[LEDGER_implementations]] (ADR-84 and ADR-94 rows)
