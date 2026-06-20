---
title: Mass scaling + HRZ lumping — validation plan (accuracy/fidelity)
project: Ladruno
status: ready-to-implement
priority: high
owner: nmora
tags:
  - validation
  - integrator
  - explicit
  - mass-scaling
---

# Mass scaling + HRZ lumping — validation plan

> **Validation plan (post-merge).** ADR 35 (HRZ, [[35_ladruno_hrz_lumped_mass_adr]])
> and ADR 36 (CentralDifferenceSMS, [[36_ladruno_selective_mass_scaling_adr]]) shipped
> in **PR #295**. The shipped tests prove the features **run, stabilize, and clean up**
> — but mass scaling's whole purpose is a **speed↔accuracy trade**, and only the speed
> side is validated. This plan closes the **fidelity** gap and converts the review's
> warned/probed behaviors into asserted tests. Follow-up PR off `ladruno` (not an amend
> to #295).

## What is already covered (do not re-do)

- HRZ kernel math: g++ standalone, 18 invariant checks (conservation, reduce-to-rowsum,
  direction-dependent α, rotational mean, guard). `tests/_hrz_verify/hrz_standalone.cpp`.
- HRZ integration: `-lump hrz` parses/runs on all 3 explicit integrators; HRZ≠diagonal
  regression guard; HRZ==rowsum on a lumped truss. `tests/test_hrz_lumped_mass.py` (6).
- SMS: constructs/runs; **tiny-element stability win** (plain CD diverges at the bulk
  step, SMS stays bounded); no-op below target; **dtor-restore**.
  `tests/test_centralDifferenceSMS_integrator.py` (4).

## The gap (why this plan)

The SMS stability test asserts the scaled run stays **bounded**, not that it **matches
the unscaled reference**. ADR 36's own Test 2 (reference-match) and Test 4
(frequency-shift) — the *accuracy* checks — were never implemented. Likewise the
Review-3 fixes (cap, energy closure, self-report skip, constraint guard) are confirmed
by a manual probe or a one-time warning, not by a regression test.

## Conventions

- **Zone-A, numpy-free** where possible (pure-`math` RMS/eigen helpers, mirroring
  `tests/test_centralDifferenceLadruno_integrator.py`) so the tests run in the lean CI.
  `ops.eigen` and `ops.nodeMass`/`ops.nodeEigenvector` are built-ins (no numpy needed).
- Bootstrap: `from _testbed import ops` (dist/bin on `sys.path` + `os.add_dll_directory`).
- Each test maps to an ADR test id; status starts `planned`.

---

## Tier 0 — fidelity (the actual gap) — Zone-A

### T-ACC — SMS reference-match  *(ADR 36 Test 2)*
- **Model:** fixed–free chain, one short/stiff/light element (tiny `dt_e`) + soft bulk,
  seeded with a small initial displacement. Choose a **modest** `dtTarget` so the added
  mass is a few % (NOT the ~2500× of the stability test) — i.e. the *useful* regime.
- **Procedure:** (a) unscaled `CentralDifferenceLadruno` at the true small Δt = reference
  history; (b) `CentralDifferenceSMS dtTarget` at the bulk Δt. Sample the free-node
  displacement at common times.
- **Assert:** `RMS(u_sms − u_ref)/RMS(u_ref) < 5%` (tighten to ~2% if achievable).
- **Why it matters:** this is the test that proves the scaled answer is *usable*, not
  merely bounded. **The single most important missing test.**

### T-MODAL — frequency-shift tradeoff  *(ADR 36 Test 4, Zone-A surrogate)*
- **Model:** a multi-DOF lumped-mass shear column (masses + `Elastic`/`Truss` springs),
  one sub-target stiff story. (The SSI soil-column version is Zone-B, deferred.)
- **Procedure:** `ops.eigen(n)` → record `f₁..f₃`. Apply SMS (run one `domainChanged` via
  a priming `analyze`), re-`eigen`. Sweep `dtTarget` over several values.
- **Assert:** Δf₁ < 1% at the modest `dtTarget`; and Δf₁ **monotonically grows** as
  `dtTarget` is pushed (the honest tradeoff curve — documents how hard a user can push).

### T-HRZ-END2END — HRZ where row-sum is genuinely wrong
- **Model:** a beam/shell transient (consistent mass, rotational DOFs) where row-sum
  zeros the rotational mass → garbage/inflated `dt_cr`.
- **Assert:** `-lump hrz` yields a finite, physically-correct `dt_cr` and a stable run at
  `dt < dt_cr,hrz`, while `-lump rowsum` mis-estimates it. Stronger than the current
  "finite/positive" + "≠ diagonal" assertions.

---

## Tier 1 — robustness (review findings → asserted) — Zone-A

### T-CAP — `-maxAddedMass` cap actually fires  *(locks SMS-CAP-DEAD)*
- Drive a model where the added-mass fraction exceeds a small `-maxAddedMass`; **capture
  stderr** (`capsys`/redirect) and assert the WARNING text appears and `frac` is the
  real (element-mass-denominator) value — not the dead `0%`. (We only probed this.)

### T-ENERGY — energy closure under scaling  *(locks the M2 claim)*
- Run SMS with an `EnergyBalanceRecorder`; assert the energy balance closes (KE computed
  against the **scaled** nodal mass), i.e. the residual stays within the recorder's
  tolerance over the run.

### T-SELFREP — self-report elements skipped  *(locks SMS-SELFREPORT)*
- Include a `LadrunoKinematicCoupling` (RBE2, self-reporting, mass-independent bound)
  below `dtTarget`; assert it is **not** scaled (`nSelfReport>0` warning emitted), and a
  pure-eigensolve element in the same model **is** scaled.

### T-CONSTR — constrained-node behavior  *(the deferred v1 guard)*
- An `equalDOF`/RBE2-tied tiny element; assert the v1 limitation warning fires and the
  documented behavior holds (dt under-delivered + which elements remain governing
  reported). Becomes a stricter test once the constraint-exclusion guard is built.

---

## Tier 2 — motivating case + deferred features (later)

### T-ZONEB — realistic refined 3D mesh  *(Zone-B, gmsh)*
- A gmsh-meshed 3D model with a **refinement zone** of tiny elements (the actual SSI /
  contact / pile use case). Assert SMS runs stably at the bulk Δt with bounded error vs a
  fine-Δt reference, and report the achieved %added-mass + Δf₁. This is "does it work on
  a *real* model," and the honest proof the feature delivers its motivating value.

### T-MPI — parallel shared-node correctness  *(gated on the MPI reduction feature)*
- OpenSeesMP, 2-rank partition with a boundary node on a scaled element; assert the
  shared-node ΔM is reduced consistently across ranks. v1 only warns; this test lands
  with the deferred MPI feature.

### T-CONSISTENT — consistent/Olovsson scaling  *(gated on the v2 feature)*
- Frequency-preservation of consistent (anisotropic) scaling vs lumped at the same
  %added-mass — the v2 selling point.

---

## Acceptance / exit criteria

- **Tier 0 green** → the accuracy claim is earned; the project memory's "tested for
  stability not accuracy" caveat can be lifted.
- **Tier 1 green** → every Review-3 fix has a regression test (not just a probe/warning).
- Tier 2 is opportunistic / feature-gated.
- Update `Ladruno_implementation/testbed/manifest.yaml` (the `CentralDifferenceSMS` row's
  test coverage) and flip the `project_ladruno_mass_scaling` memory caveat on Tier-0 pass.

## Implementation log

- **2026-06-20 — Tier 0 landed** (`tests/test_massScaling_validation.py`), two adversarial
  review rounds. Round 1 found the first cut **vacuous** (a zero-injection no-op passed the
  "accuracy" gate better than the real feature); reworked. Round 2 confirmed the 5 criticals
  dead and flagged loose bands; tightened (floor-control band, analytic magnitude check,
  honest HRZ estimate-vs-run, eigen mode guard). Tests:
  - **T-ACC** — (A) control leg: plain CD **diverges** at a supra-stable step (dt=0.011,
    above the global limit ~0.0105) while SMS holds → fails on a broken/no-op SMS; (B)
    fidelity at modest scaling within engineering tolerance, with a no-scaling
    discretization-floor control so a no-op fails the lower leg.
  - **T-MODAL** — SMS injects the **analytic `(s-1)·m`** mass (SMS-independent eigen oracle,
    <2%) + monotone dtTarget sweep lowering f1.
  - **T-HRZ** — on a consistent-mass beam, rowsum yields indefinite rotational lumped mass →
    eigensolve rejects those pairs → grossly inflated (unsafe) dt_cr; hrz≈diagonal (both
    safe); trusting rowsum's number destabilizes a run.
  - **Findings (for the report):** SMS **over-scales** (sizes against the conservative element
    pencil, not the global step); at `dtTarget` above the bulk element step the bulk elements
    scale too; the explicit run uses `system Diagonal` regardless of `-lump` (so `-lump`
    drives only the dt_cr ESTIMATE).
  - **Tier 1 deferred to a follow-up session** (T-CAP/ENERGY/SELFREP/CONSTR).

- **2026-06-20 — Tier 1 landed** (`tests/test_massScaling_validation.py`, appended below Tier 0),
  one adversarial-review round (no criticals/vacuity; strengthened loose checks per the review).
  Four robustness tests, each with a control that fails on a broken/no-op SMS:
  - **T-CAP** — the `-maxAddedMass` cap WARNING fires with the REAL element-mass-denominator
    fraction (kills SMS-CAP-DEAD's dead 0%); asserts the `%`-straddling `"% exceeds -maxAddedMass
    cap"` survives (also guards the PythonStream fix below); positive control: a 5000% cap +
    `-verbose` still runs + reports a sub-cap fraction but does NOT trip the cap.
  - **T-ENERGY** — undamped free-vibration KE+IE closure under SMS (recorder mass == integrator
    scaled mass, drift <5%) + a PINNED magnitude leg: the SMS conserved-energy uplift over the
    unscaled run equals the analytic injected nodal KE `½v₀²·Σ(injected mₓ)` within 15% (a recorder
    blind to nodal mass predicts 0 uplift → fails; a mis-scaled mass → wrong magnitude). The
    EnergyBalanceKernel sums BOTH element AND nodal mass, and SMS injects at nodes, so closure holds.
  - **T-SELFREP** — a bipenalty `LadrunoKinematicCoupling` (mass-independent self-reported bound
    < dtTarget) is SKIPPED (`nSelfReport` warning + slave node stays massless) while a normal truss
    in the same model IS scaled (injected nodeMass>0, `scaled 1/2`); guarded by a truss-only
    `dt_e < dtTarget` pre-check so a platform drift can't masquerade as a skip regression.
  - **T-CONSTR** — the one-time v1 constrained-node limitation is disclosed, and the documented
    "constrained nodes are NOT excluded" behavior is pinned (constrained node still injected) — a
    change-detector that flips when the exclusion guard is built.
  - **REAL BUG surfaced + fixed (`SRC/interpreter/PythonStream.h`):** `PythonStream::err_out` passed
    the already-formatted message AS the format to `PySys_FormatStderr(msg.c_str())`, so any literal
    `%` in an openseespy `opserr` message was eaten as a bogus printf conversion — the SMS cap
    warning (`"... % exceeds -maxAddedMass cap 5%"`) was illegible to the very Python audience the
    SMS-CAP-DEAD fix serves. Fix: `PySys_FormatStderr("%s", msg.c_str())`. Tcl `StandardStream`
    (`cerr << s`) was unaffected. Upstreamable (CWE-134); ledgered in `LEDGER_vanilla_files.md`
    (upstreamable-bugfixes) + `LEDGER_quirks.md` (incl. the residual ~1000-byte `PySys_FormatStderr`
    ceiling caveat). All 7 validation tests + a 65-test integrator/energy/coupling swath green.
  - **Tier 0 acceptance earned:** the [[project_ladruno_mass_scaling]] "tested for stability not
    accuracy" caveat can be lifted (Tier 0 reference-match + Tier 1 closure now in place).

- **2026-06-20 — Tier 2 / T-ZONEB landed** (`tests/test_massScaling_validation_zoneb.py`, Zone-B,
  raw-gmsh). The motivating-case proof on a real meshed 3D bar: a 1x1x10 hex prism (2x2 across)
  with a localized FINE band (dz=0.25, 2x finer than the coarse dz=0.5) near the clamped end — the
  SSI/contact/pile pattern. One test, `test_zoneb_refined_bar_sms_fidelity`:
  - **NECESSITY** — plain `CentralDifferenceLadruno` at the bulk dt DIVERGES (tiny elements above
    their step) while `CentralDifferenceSMS` at the same dt stays bounded (a no-op SMS diverges
    here too → the leg is non-vacuous).
  - **FIDELITY (modal)** — the ~22% fictitious mass sits in a LOW-participation zone, so f1 shifts
    **0.03%** (eigen). The honest selling point: a lot of added mass, almost no frequency cost.
  - **FIDELITY (response)** — free-end axial period (dP 0.44%) and peak amplitude (ratio 1.00) track
    the fine-dt reference. Period/peak deliberately ISOLATE the scaling fidelity from bulk-dt
    time-integration dispersion (a raw full-history RMS folds in ~12% dispersion, which is the point
    of running at the big step, not a scaling error).
  - **Report** — headline numbers (88 hexes, 16 scaled, %added, df1, dP, peak ratio) surfaced in the
    assert/print. **Measured: 22.4% added, df1 0.028%, dP 0.44%, peak 1.000.**
  - Runs in the **nightly Zone-B job** (`pytest -m zone_b`, self-hosted), NOT the PR Zone-A gate, so
    it's validated locally (Windows) + shares the gmsh hex-order assumption with the green nightly
    bend test. **T-MPI / T-CONSISTENT remain feature-gated** (need the MPI-reduction / consistent-
    scaling features first). With T-ZONEB the validation plan is COMPLETE for the v1 feature set.
