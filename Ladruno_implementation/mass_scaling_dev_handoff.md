---
title: "Mass scaling + HRZ lumping — developer handoff guide"
project: Ladruno
status: >
  SHIPPED to `ladruno`. HRZ mass-conserving lumping (ADR 35, `-lump hrz`) + selective
  mass scaling `CentralDifferenceSMS` (ADR 36, INTEGRATOR classTag 33007). Validation
  COMPLETE (ADR 37: Tier 0 fidelity #303, Tier 1 robustness #306, Tier 2 motivating-case
  T-ZONEB #308). v1.1 features: constraint-exclusion guard #311, betaK-damped sizing #314.
  NEXT = T-MPI (parallel shared-node ΔM reduction — CI-testability gap) or T-CONSISTENT
  (consistent/Olovsson anisotropic scaling). See §0.
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

## 0. Current state (2026-06-20) — handoff for the next session

**The feature is shipped and validated.** Both ADR-35 (HRZ) and ADR-36 (SMS) are on
`ladruno`, the ADR-37 validation plan is complete (Tier 0/1/2), and the two highest-value
v1.1 correctness features (constraint-exclusion guard, betaK-damped sizing) are in.

| Increment | What | PR |
|---|---|---|
| HRZ + SMS | `-lump hrz`, `CentralDifferenceSMS` (classTag 33007) | #295 |
| Tier 0 | fidelity: reference-match, (s−1)·m modal magnitude, HRZ-rowsum-unsafe | #303 |
| Tier 1 | robustness: cap fires, energy closure, self-report skip, constraint disclosure | #306 |
| Tier 2 | T-ZONEB: SMS on a gmsh-meshed refined 3D bar (motivating case) | #308 |
| v1.1 | constraint-exclusion guard (MP slave nodes) | #311 |
| v1.1 | betaK-damped sizing (closed form) | #314 |

**Also fixed along the way (#306):** an upstream openseespy bug — `PythonStream::err_out`
fed the message AS the printf format to `PySys_FormatStderr`, so any literal `%` in an
`opserr` message was silently eaten (the SMS cap warning was illegible). Fix:
`PySys_FormatStderr("%s", msg)`. Upstreamable (CWE-134). See [[LEDGER_quirks]].

**What's NEXT (both substantial, neither started):**
- **T-MPI — parallel shared-node ΔM reduction.** v1 is sequential / partition-interior
  only; a partition-boundary node gets only rank-local Δmₑ and the only MPI reduction is
  the scalar dt, so shared-node masses desync across ranks. The fix is `MPI_Allreduce` of
  the per-shared-node injected ΔM. **Caveat: the Zone-A CI gate is single-process**, so
  T-MPI (a 2-rank `mpiexec` test) cannot be gated in CI — validate locally only. First do a
  read-only investigation of how OpenSeesMP assembles shared-node mass (it may already sum
  `setMass` contributions, shrinking or removing the need).
- **T-CONSISTENT — consistent/Olovsson anisotropic scaling.** v1 injects an isotropic
  diagonal Δm. Consistent (anisotropic) scaling preserves frequencies better at the same
  %added-mass — the "v2" selling point. Larger change (inject a consistent ΔM matrix), but
  sequential and fully Zone-A-testable.

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
- **MPI shared-node reduction:** start read-only — determine whether OpenSeesMP already sums
  shared-node `setMass` contributions across ranks (Subdomain / parallel assembly). If yes,
  the "limitation" is overstated and the fix is small; if no, `MPI_Allreduce` the per-shared
  -node ΔM in `domainChanged`. Remember T-MPI cannot gate in the single-process CI.
