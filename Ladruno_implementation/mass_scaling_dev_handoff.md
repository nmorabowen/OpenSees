---
title: "Mass scaling + HRZ lumping — developer handoff guide"
project: Ladruno
status: >
  SHIPPED to `ladruno` — the explicit-integrator SMS axis is COMPLETE. HRZ lumping
  (ADR 35, `-lump hrz`) + selective mass scaling, LUMPED and CONSISTENT (Olovsson),
  on ALL THREE explicit integrators: CentralDifference (33007/33008), ExplicitBathe
  (33009/33010), ExplicitBatheLNVD (33011/33012). ADR 36 (lumped), ADR 38 (consistent),
  ADR 37 validation. PRs #295/#303/#306/#308/#311/#314 (CD lumped+validation),
  #320 (CD consistent), #324 (ExplicitBathe + LNVD families; subsumed #322). All merged.
  V4 (consistent-mass energy-recorder KE closure) DONE. T-MPI: the LUMPED path is
  VALIDATED parallel-correct (np=1 vs np=2 bit-identical; the distributed/MPI diagonal solver
  sums shared-node ΔM at solve time) and the stale "not reduced" warnings were corrected;
  the CONSISTENT (Olovsson) path now has a DISTRIBUTED PCG (V5) — VALIDATED under OpenSeesMP
  `system MPIDiagonal` (MPI np=2 == serial reference; PR #335). SMS is now also reachable from
  the legacy Tcl `integrator` parser (PR #340). The SMS axis is COMPLETE for every realistic
  config (serial + OpenSeesMP, lumped + consistent, openseespy + Tcl). OpenSeesSP
  (DistributedDiagonalSOE) was INVESTIGATED and DEFERRED — no local validation path + no
  demonstrated explicit-SP demand (see §0). NOT CI-gated (single-process CI). See §0.
tags:
  - integrator
  - explicit
  - mass-scaling
  - handoff
---

# Mass scaling + HRZ lumping — handoff to the next developer

Companion to the ADRs ([[35_ladruno_hrz_lumped_mass_adr]], [[36_ladruno_selective_mass_scaling_adr]],
[[37_ladruno_mass_scaling_validation_plan]]). This consolidates what shipped, where the
code lives, the non-obvious decisions, the validated behavior, and what is left.

---

## 0. Current state (2026-06-21) — handoff for the next session

> **ADR-52 addendum (2026-07-02) — read this box before the table below.** The six
> `ExplicitBathe*` classes listed in this document were **collapsed into ONE flag class**
> (ADR-52 W1-E2, PR #419): the primary command is now
> `integrator ExplicitBathe $p <-lnvd <alpha>> <-sms $dtTarget> <-consistent>` (classTag
> 33000). The five other spellings/tags (`ExplicitBatheLNVD` 33002, `ExplicitBatheSMS`
> 33009, `ExplicitBatheSMSConsistent` 33010, `ExplicitBatheLNVDSMS` 33011,
> `ExplicitBatheLNVDSMSConsistent` 33012) remain **deprecated-but-recognized aliases**
> for one release; both brokers route all six tags through `ExplicitBathe::makeForBroker`.
> The CentralDifference family kept its classes (33003/33007/33008). The 2026-07-01
> explicit review then landed three fix PRs on top: **#468** (betaKinit/betaKcomm damped
> SMS sizing, `revertToLastStep()`, NaN-capable circuit breaker), **#472** (`update()`
> misuse guard, unified `-sms` Diagonal sizing default, sub-step-2 solve return codes,
> serialization superset), and the kernel-guard + diagnostics batch PR (consistent-builder
> ndf clamp, lumped Σndf pre-walk, PRE-/POST-SCALING dt_cr warning reword, `-recompute` N
> consume, PCG-breakdown warning, DGGEV complex-pair skip, prevKE reset, getVel guards).
> Shipped log: [[52_ladruno_integrator_strengthening_adr]]. Everything below this box is
> the pre-collapse per-class picture — still accurate on mechanism, superseded on naming.

**The explicit-integrator SMS axis is COMPLETE and merged to `ladruno`.** Selective mass
scaling exists in BOTH flavors — **lumped** (conventional/DT2MS, ADR 36) and **consistent**
(Olovsson, ADR 38) — on ALL THREE Ladruno explicit integrators. Six integrator classTags:

| Integrator | Lumped | Consistent (Olovsson) | shipped via |
|---|---|---|---|
| `CentralDifference[Ladruno]` | `CentralDifferenceSMS` **33007** | `CentralDifferenceSMSConsistent` **33008** | #295/#303/#306/#308/#311/#314 ; #320 |
| `ExplicitBathe` (Noh-Bathe) | `ExplicitBatheSMS` **33009** | `ExplicitBatheSMSConsistent` **33010** | #324 (subsumed #322) |
| `ExplicitBatheLNVD` (+FLAC) | `ExplicitBatheLNVDSMS` **33011** | `ExplicitBatheLNVDSMSConsistent` **33012** | #324 |

All on `ladruno` (confirmed: `git show origin/ladruno:SRC/classTags.h` has all 6 tags;
the 8 EB/LNVD family files are in the tree). ADR-37 validation (Tier 0/1/2) complete for
the CD lumped path; the consistent + EB/LNVD families are validated by their own Zone-A
suites (see §4).

**The consistent (Olovsson) architecture — the one big idea (ADR 38).** Lumped scaling
adds an isotropic diagonal `(s−1)·mₐ` to each throttling element's nodes, which shifts ALL
global frequencies incl. f1. Consistent scaling injects the **centroidal** scaling mass
`M̄ₑ = β[diag(mₐ) − m mᵀ/Mₑ]` whose **row sums are zero** ⇒ `M̄·t_rigid = 0` ⇒ rigid
translation gets no added inertia ⇒ **f1 preserved** while only the relative/deformation
modes (which govern the explicit step) are loaded. Measured (transient FFT, oracle Case-A
bar): **f1 −0.17% (consistent) vs −53.4% (lumped)** at the same dtTarget — identical across
all three integrators.

`M̃ = M_lump + ΣM̄ₑ` is SPD but **non-diagonal**, so it can't live in `Node` mass and
`system Diagonal` can't invert it trivially. The solve is the key mechanism, shared by all
three consistent integrators:
- A protected **no-op `virtual refineAccel(Vector&)` hook** was added to each explicit base
  (`CentralDifferenceLadruno`, `ExplicitBathe`, `ExplicitBatheLNVD`), called at every site
  that consumes a diagonal `a = M_lump⁻¹r` solve (CD: starter + update; the Noh-Bathe pair:
  both sub-steps `A_tpdt`/`A_tdt`). **Default no-op ⇒ the base + lumped paths are
  byte-identical** (a guarantee, not something to prove; base regressions confirm it).
- The consistent subclass overrides it to: read the **factored `DiagonalSOE`** (post-solve
  `getDiagonalA()` = `1/mass`, the Jacobi preconditioner AND the way to recover the RHS as
  `r = a / Ainv`), then run a **matrix-free Jacobi-preconditioned CG** (`Ladruno::consistentPCG`)
  to replace `a` with `M̃⁻¹r`. Converges in 3–21 iters (M̄ is a small perturbation of the
  dominant lumped diagonal). It NEVER mutates `Node` mass ⇒ no inject/restore lifecycle.
- All sizing (`dt_e`, self-report skip, betaK-damped `s = T²+2Tc`, MP-node exclusion, cap)
  is the ONE shared kernel `LadrunoMassScaling.h` (`buildMassScaling` / `applyMassScaling`
  for lumped; `buildMassScalingConsistent` / `consistentMatVec` / `consistentPCG` +
  `struct ConsistentBlock` for consistent). Oracle: `Ladruno_implementation/mass_scaling_consistent/`.

**Two adversarial-review findings folded (don't re-discover):**
- **CONSISTENT excludes BOTH slave AND master MP nodes** (lumped excludes only slaves).
  Under the default `Transformation` handler an element touching an MP *master/retained*
  node has a transformed `FE_Element::getID()` LARGER than the element `n` → pairing it with
  the n×n `M̄` was an OOB read. Fix: exclude slave+master + a defensive `getID().Size()==n`
  guard (also keeps `M̄` off every ADR-30 projector equation). See [[LEDGER_quirks]].
- **LNVD factories validate `alpha ∈ [0,1)`** to match the parent `ExplicitBatheLNVD` contract.

**Also fixed earlier (#306):** an upstream openseespy bug — `PythonStream::err_out` fed the
message AS the printf format to `PySys_FormatStderr`, eating literal `%` in `opserr`. Fix:
`PySys_FormatStderr("%s", msg)`. Upstreamable (CWE-134).

**What's NEXT:**
- **V4 — consistent-mass energy closure. DONE (this session).** The `EnergyBalanceRecorder`
  builds KE from node/element `getMass()` (the lumped diagonal), so it could not see the
  cross-node `M̄ₑ` that the consistent path keeps only inside the integrator (matrix-free PCG
  operand). The recorder holds only a `Domain*` (no integrator handle), so the fix is a
  process-global **`Ladruno::MassScalingEnergyRegistry`** (`SRC/analysis/integrator/
  LadrunoMassScalingEnergy.{h,cpp}`): the active consistent integrator `publish()`es its
  per-element node-major `M̄ₑ` (keyed by element tag — `ConsistentBlock` gained an `eleTag`)
  at the end of `domainChanged` and `clear()`s on teardown (owner-guarded so a stale dtor
  can't wipe a newer publisher); the shared `EnergyBalanceKernel.h` queries it per element
  and adds the missing `½vᵀM̄ₑv` (vel is already node-major — no equation-number work).
  Empty for the lumped path + every base integrator (`active()==false`) ⇒ byte-identical
  there, no double count. Wired on all three consistent integrators. Validated by
  `tests/test_massScaling_consistent_energy.py` (Zone-A 4/4): an analytic `½vᵀM̃v` oracle
  (recorder KE == lumped+`M̄` to ~machine precision, all 3 integrators, `M̄`≈77% of KE on a
  deformation-rich IC) + CD mechanical-energy conservation (drift <5% vs the pre-V4
  reconstruction `½vᵀM_lump v + IE`, which oscillates with the omitted `M̄`). See [[LEDGER_quirks]].
- **T-MPI — parallel shared-node ΔM reduction. LUMPED path: VALIDATED parallel-correct;
  CONSISTENT path: distributed PCG IMPLEMENTED + VALIDATED (V5). The parallel SMS axis is
  COMPLETE for OpenSeesMP.** The earlier worry
  ("a partition-boundary node gets only rank-local Δmₑ … shared-node masses desync") was
  **WRONG for the lumped path.** Why: in a parallel build `system Diagonal` auto-swaps to
  `DistributedDiagonalSOE` (OpenSeesSP) and `system MPIDiagonal` → `MPIDiagonalSOE`
  (OpenSeesMP); **both solvers SUM the shared-DOF diagonal across ranks at solve time**
  (`DistributedDiagonalSolver`: gather→P0→`+=`→broadcast; `MPIDiagonalSolver`:
  `intersectionsAB` accumulates A on the first solve). `CentralDifferenceLadruno` reads M
  *through* `formTangent → SOE → solve` (the **reduced** diagonal), NOT raw `Node::getMass()`.
  Each element lives wholly on one rank ⇒ per-element `dt_e`/`s` are correct rank-locally;
  SMS injects `(s−1)·m` into each rank's local node copy via `setMass`; the solve-time
  reduction then sums all ranks' contributions to exactly the right physical total.
  **Empirical proof:** `Ladruno_implementation/mass_scaling_mpi/` — a 1D fixed-free bar with
  a fine zone straddling the central node, run at `np=1` (whole bar, all-local injection) vs
  `np=2` (split at the shared node, where elements 10/rank0 and 11/rank1 BOTH inject ΔM into
  shared node 11). Tip-disp history matches **bit-identical** (`max |abs diff| = 0.000e+00`,
  150 steps, 36.6% added mass on the fine zone). Run: `OpenSeesPyMP` + `system MPIDiagonal`
  (the `pyd` lives in `dist/openseesmp/`; see `run.ps1`).
  **Stale warnings CORRECTED** in the three lumped integrators (`CentralDifferenceSMS`,
  `ExplicitBatheSMS`, `ExplicitBatheLNVDSMS` `.cpp` limitation (3) + `CentralDifferenceSMS.h`
  scope comment): they used to say "parallel shared/boundary nodes are not mass-reduced
  across ranks", now they state the truth (lumped IS reduced via the distributed/MPI diagonal
  solver; consistent is not).
- **CONSISTENT (Olovsson) parallel — V5, DONE.** The serial `consistentPCG`/`consistentMatVec`
  ARE rank-local (local `res^z`/`p^Ap`, no shared-DOF `M̄` exchange), so a **distributed PCG**
  was added: `consistentParPCG`/`consistentParMatVec` + `ConsistentParOps` in
  `LadrunoMassScaling.h`, driven from the shared `LadrunoConsistentRefine.h` (one body for all
  3 consistent integrators, serial + MPI). **The one idea** is the weight `wᵢ=1/multiplicityᵢ`:
  the matvec applies the GLOBAL (replicated) lumped diagonal WEIGHTED + off-diagonal `M̄ₑ` in
  FULL then `assembleSharedSum`s shared DOFs across ranks (diagonal → full-once, off-diagonal
  accumulates); inner products use the same `w` + `globalReduceSum` (all-reduce). Every CG
  control scalar is global ⇒ identical iter count on all ranks ⇒ collectives lockstep (no
  deadlock); `w≡1`/no-op assemble/identity reduce ⇒ collapses to the serial PCG at `np=1`.
  **Architecture gotcha (see [[LEDGER_quirks]]):** the integrators are in the shared
  `OpenSeesLIB`, which can NOT `#ifdef _PARALLEL_INTERPRETERS` nor reference `MPIDiagonalSOE`
  (both exist only in the MP executables). Dispatch goes through 4 new `LinearSOE` base virtuals
  (serial no-op defaults; `MPIDiagonalSOE` overrides). **Validated** (`mass_scaling_mpi/
  run_consistent.ps1`): MPI `np=2` == serial `DiagonalSOE`+`consistentPCG` gold reference
  (max abs diff 0 to recorder precision) AND `np=1`==`np=2` AND consistent (1.746e-3) ≠ lumped
  (2.35e-3); serial Zone-A 36/36 green. **Caveat: the Zone-A CI gate is single-process**, so
  the 2-rank `mpiexec` tests (lumped + consistent) cannot be gated in CI — local-only validation.
  **Tcl-parser gap FIXED (PR #340):** `integrator CentralDifferenceSMS …` used to error in the
  *Tcl* `OpenSees.exe`/`OpenSeesMP.exe` ("No Integrator type exists") — SMS was registered only
  in the interpreter/openseespy layer (`OpenSeesCommands.cpp`), not the legacy
  `SRC/tcl/commands.cpp` `specifyIntegrator()` dispatch, despite the splash banner advertising it.
  Now all 6 are wired into the Tcl parser (mirrors the `CentralDifferenceLadruno` branch);
  smoke-tested via `OpenSees.exe Ladruno_implementation/mass_scaling_mpi/sms_tcl_smoke.tcl`.

**Merge mechanics lesson (this session).** #320/#322/#324 were a `--base ladruno` stack.
After #320 squash-merged, #322 went CONFLICTING (its branch carried #320's pre-squash
commits). Resolved by the CLAUDE.md-endorsed recovery: rebase/merge the **leaf** (#324)
onto current `ladruno` (one conflict resolution — all conflicts were empty-`theirs`
"keep-ours" in the shared registration/ledger files), merge #324 carrying BOTH the EB and
LNVD families, and **close #322 as superseded** (its content shipped via #324). For future
stacks: merge the bottom first, then immediately rebase the leaf; or just land one leaf PR.

---

## 1. Files & how they relate

```
SRC/analysis/integrator/
  CriticalTimeStep.{h,cpp}    SHARED kernel — lumpElementMass, elementLambdaMax,
                              elementCriticalDt (undamped 2/w), computeCriticalTimeStep
                              (global dt_cr query, damped+undamped). SMS and the dt_cr
                              query lump/eigensolve IDENTICALLY through this.
  LadrunoMassLumping.h        HRZ kernel (ADR 35): Ladruno::hrzLump + OpenSees-free
                              hrzLumpRaw (g++-verified). CTSLumping::HRZ mode.
  LadrunoMassScaling.h        SMS core (ADR 36): buildMassScaling (sizing + the guards) +
                              applyMassScaling (commit/restore nodal ΔM). HEADER-ONLY.
                              ConsistentBlock now also carries eleTag (V4 energy conduit).
  LadrunoMassScalingEnergy.{h,cpp}  V4 energy conduit: process-global owner-guarded
                              MassScalingEnergyRegistry (eleTag -> node-major M̄ₑ) so the
                              EnergyBalanceRecorder can add ½vᵀM̄ₑv (consistent KE closure).
  CentralDifferenceSMS.{h,cpp}  The integrator. Subclass of CentralDifferenceLadruno via a
                              protected classTag ctor. domainChanged injects; dtor restores.

tests/
  test_massScaling_validation.py     8 tests (ADR 37): T-ACC, T-MODAL, T-HRZ (Tier 0);
                                     T-CAP, T-ENERGY, T-SELFREP, T-CONSTR (Tier 1);
                                     T-BETAK (betaK sizing). zone_a.
  test_massScaling_validation_zoneb.py  T-ZONEB (Tier 2, gmsh refined 3D bar). zone_b.
  test_centralDifferenceSMS_integrator.py  SMS-1..4 (construct/run, tiny-element stability
                                     win, no-op below target, dtor-restore). zone_a.
  test_hrz_lumped_mass.py            6 HRZ tests. tests/_hrz_verify/hrz_standalone.cpp.
```

The `-lump hrz` option is wired on all 3 explicit integrators (CDL / ExplicitBathe /
LNVD); it only changes the dt_cr **estimate** — the explicit run uses `system Diagonal`
regardless (diagonal-of-consistent). HRZ's distinct value (mass conservation on
rotational DOFs) is kernel-tested, not a bigger safe step for a Diagonal run.

---

## 2. HRZ lumping (ADR 35)

Per-direction α_d = m_d / S_d (S_d = sum of the diagonal mass in direction d); rotational
DOFs scaled by the mean translational α; u-p / shell non-positive DOFs pass through
(NOT a whole-element fallback). Reduces to row-sum on a regular element. The point: row-sum
zeros rotational mass on consistent-mass beams/shells → the dt_cr eigensolve rejects the
indefinite pairs → grossly inflated (unsafe) dt_cr. HRZ keeps rotational mass positive.

> **T-HRZ pins this:** on a consistent-mass beam, `-lump rowsum` reports a dt_cr ~3× the
> true diagonal-run limit (running at 0.9× it diverges); `-lump hrz` ≈ `-lump diagonal`
> (both safe). The run mass is diagonal-of-consistent either way; rowsum's *estimate* is
> the hazard.

---

## 3. The SMS integrator (ADR 36)

`integrator CentralDifferenceSMS $dtTarget <-maxAddedMass f> <-lump rowsum|diagonal|hrz> <-tangent> <-verbose>`

For every element whose per-element stable step `dt_e < dtTarget`, scale its mass by `s`
so `dt_e → dtTarget`, implemented by ADDING fictitious nodal mass `(s−1)·m` to the
element's nodes (additive nodal mass — the canonical M every consumer reads: the leap-frog
M⁻¹ AND the EnergyBalanceRecorder KE). The increment is recorded per node so it is
re-baselined on a re-entrant `domainChanged` and **restored on integrator destruction** (M5
— never a permanent Domain mutation; corrupts a later gravity→modal→SMS→post workflow
otherwise).

### 3a. Sizing — undamped → betaK-damped (the #314 closed form)

Undamped: `s = (dtTarget/dt_e)²`. But stiffness-proportional (betaK) Rayleigh damping
shrinks the explicit step (ξ = betaK·ω/2 grows with ω), so the undamped scale UNDER-scales
→ the element is still unstable at dtTarget. The betaK-damped step at scale `s` is

```
dt_d(s) = (2/ω₀)·(√(s + c²) − c),    c = betaK/dt_e  (= ½·betaK·ω₀,  ω₀ = 2/dt_e)
```

which inverts in **closed form** (no iteration):

```
s = T² + 2·T·c,    T = dtTarget/dt_e
```

- Reduces to the undamped `T²` when `betaK = 0` → **no-damping models byte-identical**.
- The SKIP test uses the **damped** step: skip iff `dt_d(s=1) ≥ dtTarget` (an element whose
  *undamped* step clears dtTarget but whose *damped* step does not is now correctly scaled).
- **alphaM (mass-proportional) is intentionally EXCLUDED**: ξ_α = alphaM/2ω *decreases*
  with ω so it does not reduce the governing high-frequency step, and folding it in across
  scales is non-monotonic (`dt_d(s)` would asymptote, possibly never reaching dtTarget).
- Only `betaK = ele->getRayleighDampingFactors()(1)` enters — mirroring the betaK term of
  `computeCriticalTimeStep`'s damped estimate, so SMS sizes consistently with what
  `criticalTimeStep()` reports.

### 3b. The guards (skip + report), in `buildMassScaling`, in order

1. **Self-report** (`getExplicitCriticalTimeStep() > 0`): an element with a MASS-INDEPENDENT
   bound (bipenalty coupling: 2√(mₚ/k)) cannot be helped by adding mass → skip, count
   `nSelfReport`. (This is how RBE2/RBE3 couplings are handled — they are ELEMENTS, not
   MP constraints.)
2. **betaK-damped skip** (3a): skip if already stable incl. damping.
3. **Constraint exclusion (#311):** skip any sub-target element touching an MP-constraint
   SLAVE node (`Domain::getMPs()->getNodeConstrained()` — equalDOF / rigidDiaphragm /
   rigidLink / generic MP). A slave DOF is eliminated by the handler and the injected mass
   is redistributed to the retained node via `Tᵀ M T`, so the dt boost would NOT land →
   silent mass mis-distribution + instability. Excluded elements are counted (`nConstrained`)
   and reported as STILL GOVERNING at their un-scaled (damped) `dt_e`. **SP/fix nodes are
   NOT excluded** (mass on a removed DOF is inert; the motivating SSI/pile case refines at
   fixed supports). The retained (master) node is also not excluded (consistent with SMS's
   existing constraint-blind dt_e approximation).
4. **Partial-inject guard:** stage each element's increments locally; commit only if the
   node-by-node DOF walk maps (`pos == n`). A non-node-major / null-node element aborts with
   NOTHING committed (vs. a half-injected, un-rolled-back state), counted `nMismatch`.

### 3c. The `-maxAddedMass` cap (the SMS-CAP-DEAD fix)

`frac = addedMass / modelMass`. The denominator MUST be the **element** translational
lumped mass summed over elements — nodal `getMass()` is **0** on `-rho` element models, so
the original nodal denominator gave a dead 0% cap. Over the cap → a WARNING (proceeds; the
inertia shifts frequencies). The cap is a warning, not an abort.

### 3d. Rejected options

`-cflAbort` / `-recompute` are REJECTED at parse (MF-1): their inherited path re-runs the
element-mass eigensolve, which cannot see the nodal augmentation, and would spuriously
abort a run that is in fact stable.

---

## 4. Validation coverage (ADR 37) — the 8 + 1 tests

Each test has a control that FAILS on a broken / no-op / blind SMS (the Tier-0 review
caught a vacuous first cut; every test is now a genuine change-detector).

| Test | Proves | Non-vacuity |
|---|---|---|
| T-ACC | supra-stable step: plain CD diverges, SMS holds; modest-scaling fidelity vs fine-dt ref | no-scaling floor control |
| T-MODAL | SMS injects the analytic (s−1)·m (eigen <2%) + f1 falls monotonically with dtTarget | analytic eigen oracle |
| T-HRZ | rowsum's dt_cr is unsafe on a consistent beam; hrz≈diagonal | run at 0.9× rowsum diverges |
| T-CAP | cap fires with the REAL element-mass-denominator frac; the `%`-text survives | 5000% cap control must NOT warn |
| T-ENERGY | KE+IE closure under SMS + the uplift = analytic injected nodal KE within 15% | a nodal-mass-blind recorder ties the levels |
| T-SELFREP | bipenalty coupling skipped while a normal truss is scaled (`1/2`) | nodeMass injection + count |
| T-CONSTR | MP-slave element EXCLUDED (slave stays massless), free truss still scaled | targeted control |
| T-BETAK | betaK injects MORE than the undamped formula, matching `(T²+2Tc−1)/(T²−1)` | a betaK-blind SMS ties them |
| T-ZONEB | SMS on a real refined 3D bar: necessity + Δf₁<1% + period/peak fidelity (zone_b) | CD@bulk diverges |

T-ZONEB measured: 88 hexes, 16 fine scaled, 22.4% added mass, Δf₁ 0.028%, period dP 0.44%,
peak ratio 1.000.

---

## 5. Quirks & gotchas (see [[LEDGER_quirks]] for the full entries)

- **The explicit run uses `system Diagonal` regardless of `-lump`.** `-lump` drives only the
  dt_cr ESTIMATE. So `-lump hrz` does not give a bigger safe step for the Diagonal run.
- **SMS over-scales.** It sizes against the conservative per-element pencil `dt_e`, well
  below the true global stable step → it adds more mass than strictly needed → bigger
  frequency shift. Accuracy is good only at modest scaling / when the tiny elements are a
  small, low-participation fraction of the model (the T-ZONEB regime: fine zone at the clamp,
  Δf₁ ≈ 0, despite ~22% added mass).
- **openseespy `opserr` ate literal `%`** (fixed #306). When asserting on warning TEXT under
  openseespy, remember pre-fix binaries dropped `%`; CPython's `PySys_FormatStderr` still
  caps output at ~1000 bytes (per-token emission makes that essentially never bite).
- **`capfd`, not `capsys`.** opserr → StandardStream → `std::cerr` (fd 2). Use pytest `capfd`
  (fd-level) to capture warnings; `capsys` (Python sys.stderr) misses them.
- **betaK-damped sizing excludes alphaM** (§3a) — by design, not an omission.
- **Constraint guard is MP-slave only**, not SP, not the retained node, not RBE2/RBE3
  elements (those go through the self-report skip).

---

## 6. Recipes

- **Build** (only if C++ changes; tests-only need NO build): from the worktree ROOT via the
  PowerShell tool — `Set-Location '<root>'; cmd /c "Ladruno_scripts\build.bat OpenSeesPy"`.
  Bash and PowerShell tools SHARE cwd; a stray Bash `cd` breaks the build path. pyd lands at
  `dist/bin/opensees.pyd`.
- **Run tests** (no build for new tests): `pythoncore-3.12-64\python.exe` with
  `os.add_dll_directory(DIST)` + `sys.path.insert(0, DIST)` + `sys.path.insert(0, 'tests')`,
  DIST = `<worktree>\dist\bin`, then `pytest.main([...])`. (Zone-A: `-m zone_a`; the Zone-B
  T-ZONEB runs `-m zone_b`, nightly self-hosted — NOT in the PR gate.)
- **CI gates:** `ci/check_classtags.py`, `ci/check_manifest.py`. The PR CI = build +
  `pytest -m zone_a` + gates; Zone-B is nightly self-hosted.
- **Branch hygiene:** the fork auto-merges fast (squash). One logical PR per branch off
  `ladruno`; verify `gh pr view <n> --json state` is not MERGED before any follow-up push
  (stranded-commit trap). If a PR goes CONFLICTING, `git merge origin/ladruno` and re-push;
  watch for inherited G9 manifest-row debt + Linux-only build breaks a sibling PR left on
  `ladruno` HEAD.

---

## 7. Where to start the next increment

- **Consistent/Olovsson scaling (recommended next):** sequential, Zone-A-testable, real
  value. Inject a consistent ΔM matrix (not a diagonal) sized to preserve the element's
  frequencies; pair with `T-CONSISTENT` (frequency-preservation vs lumped at the same
  %added-mass). `applyMassScaling` currently writes diagonal-only via `Node::setMass`; the
  consistent path needs the off-diagonal nodal mass route.
- **MPI shared-node reduction: DONE — lumped AND consistent.** Lumped is parallel-correct via
  the distributed/MPI diagonal solver's solve-time shared-DOF sum; consistent got a distributed
  CG (V5: `consistentParPCG`/`consistentParMatVec` with global all-reduced inner products + a
  shared-DOF `M̄` assemble, weight `wᵢ=1/multiplicityᵢ`) — both validated 2-rank `mpiexec`,
  adversarially reviewed clean (PR #335). See §0 + `mass_scaling_mpi/`.
- **Tcl-parser wiring: DONE (PR #340).** All 6 SMS integrators reachable from the legacy Tcl
  `integrator` parser, not just openseespy.
- **OpenSeesSP (DistributedDiagonalSOE) consistent path: DEFERRED (decision 2026-06-21).** See
  the ADR-38 log entry — codeable (the integrator/helper layer is already SP-agnostic; just add
  the 4 `LinearSOE` overrides to `DistributedDiagonalSOE` via Channel gather→P0→sum→broadcast,
  noting that solver leaves `A` as *mass* so `getScalingDiagonalA` needs a cached inverse), but
  NOT responsibly shippable now: no local explicit-SP validation harness, no demonstrated demand
  (SP examples are all `Mumps`+implicit), and it would put untested distributed numerics in
  vanilla core. Revisit only with a real `mpiexec OpenSeesSP` harness + a concrete use case.
  **The mass-scaling axis is otherwise COMPLETE.**
