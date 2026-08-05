---
title: Ledger — OpenSees quirks & gotchas
project: Ladruno
tags:
  - ledger
  - quirks
  - gotchas
---

# Ledger — OpenSees quirks & gotchas we learned

Surprising, undocumented, or bug-prone behaviours of upstream OpenSees that
cost us time. Recording them here so we (and future us) stop re-discovering
them. This is observation-only — fixes we actually applied are tracked in
[[LEDGER_vanilla_files]] / [[LEDGER_implementations]].

## Conventions

- **One section per quirk.** Title = the symptom you'd search for.
- State: *what bites*, *why*, *workaround/status*, and the *date* learned.
- If a quirk drove a code change, cross-link the ledger row / PR.
- Deep build/toolchain quirks may live in
  [[../Ladruno_internal/01_compilation_journal]]; link rather than duplicate.

## Quirks

### A new element registered ONLY in `OpenSeesElementCommands.cpp` (Python functionMap) is SILENTLY unreachable from the Tcl `OpenSees.exe` — the `element` command has a SEPARATE dispatch and the miss looks like a DLL-load failure, not a wiring bug
- **Bites:** you add a modern C++ element (`void* OPS_Foo()`), register it in `SRC/interpreter/OpenSeesElementCommands.cpp`'s `functionMap` (which serves **openseespy** + the interpreter runtime), smoke-test it in Python — `ops.element('Foo', ...)` works — and ship. But `OpenSees.exe` (classic Tcl) prints `checking library: OPS_Foo` then `ERROR -- element of type Foo not known` and fails. It reads like a missing DLL / build-staleness problem; it is actually a **second, independent registration** you never did. (Real instance: `BezierTri6`/`BezierTet10` were functionMap-only for weeks; only caught by a Tcl smoke test after a clean rebuild "should have" fixed it — the rebuild couldn't, because the source was never wired for Tcl.)
- **Why:** two element dispatchers coexist. (1) `OpenSeesElementCommands.cpp` `functionMap` → Python/interpreter. (2) The classic Tcl `element` command lives in `SRC/element/TclElementCommands.cpp` as a big `strcmp(argv[1], ...)` if/else chain; on a miss it calls `OPS_GetElementType()` (`elementAPI_TCL.cpp`), which searches only the **legacy C `eleObj`/function-pointer** registry (NOT the modern functionMap), and on a further miss falls through to `getLibraryFunction()` — hence the misleading `checking library: OPS_Foo` message before "not known". Modern C++ `OPS_*` factories are never in that legacy registry, so a Tcl-side element MUST be added to the `TclElementCommands.cpp` chain explicitly.
- **Extra trap:** don't fix it with a fresh `else if` — that chain sits at MSVC's **C1061 "blocks nested too deeply"** limit (one-branch-per-element already broke the Windows build once). The fork routes all its elements through a single `else if` guarded by `strncmp(argv[1],"Ladruno"/"Bezier",...)` that dispatches a **flat factory table** (`ladrunoElementTable`); add new fork elements as a one-line table row + an `extern void *OPS_Foo(void);` decl, never a new branch.
- **Workaround/status (2026-07-08):** `BezierTri6`/`BezierTet10` wired into the Tcl table + guard broadened to match `Bezier`/`bezier`; both verified constructing in `OpenSees.exe` and openseespy. Rule going forward: **every new element needs BOTH registrations** (functionMap for Python, `TclElementCommands.cpp` table for Tcl) — and a Tcl smoke test, not just a Python one. See [[LEDGER_vanilla_files]] `TclElementCommands.cpp` rows.

### `algorithm Linear -factorOnce` cuts explicit CDL wall ~18% bit-identically — but a mid-run `domainChanged` leaves the skipped tangent STALE (never combine with re-emission / element removal / anything that resizes the SOE)
- **Bites (the win):** the standard explicit recipe (`system Diagonal` + `algorithm Linear` + `CentralDifferenceLadruno`) reassembles the **constant** lumped-mass effective tangent into the SOE every step — measured 16.2% of explicit wall (2500-step LadrunoBrick+J2 wave run). Vanilla `Linear` already has `-factorOnce` (forms once, then skips — `Linear.cpp:110`): measured **79.9 → 65.6 s (−17.9%), md5-bit-identical** displacement field.
- **Bites (the trap):** `factorOnce` has **no domainChanged reset** — nothing in `Linear` or the analysis re-arms it. After a mid-run `domainChanged` (ADR-60 contact re-emission fires one from `Domain::commit()`; ADR-51 removal; SMS mass injection re-handles) the SOE is re-sized/zeroed but the tangent re-form is skipped → solve against a stale/empty matrix. The invalidators are exactly the events that change **M** or resize the SOE: `domainChanged`, lumped-SMS injection, `updateParameter` on density. **Material/geometric nonlinearity is NOT an invalidator** — CDL's assembled matrix is pure mass (`formEleTangent` = `zeroTangent()+addMtoTang()`, `CentralDifferenceLadruno.cpp:227-232`); all nonlinearity lives in the residual, which `-factorOnce` never touches (the −17.9% A/B was a *yielding* J2 run, bit-identical). Δt changes do NOT invalidate the matrix either (no dt in pure M) — but they independently require `revertToLastStep` for the leap-frog reseed (separate quirk above).
- **Why:** `factorOnce` was designed for static linear re-analysis; the flag latches 1→2 on first form and only a convergence *failure* path (ModifiedNewton's, not Linear's) ever resets it.
- **Workaround/status (2026-07-06, 40b lane-D drill-down — SUPERSEDED same day):** the win is now **built into CentralDifferenceLadruno as a default-ON constant-mass tangent cache** (ADR-67 P-NEW-1, adversarially gated): invalidated in `domainChanged`/`setConstraintProjector`/`revertToLastStep`, so re-emission / removal / staged / SMS models get it safely too; opt out with `integrator CentralDifferenceLadruno ... -noMassCache`. `-factorOnce` on CDL decks is therefore obsolete (and still carries this trap — prefer the integrator cache). The ONE event the cache cannot see: `updateParameter` on a DENSITY (no domainChanged fires) — use `-noMassCache` there; a one-time runtime note fires when the deck declares parameters. See [[40b_phase0_dominance_report]] §lane-D addendum (which also records the 2nd-constitutive-pass finding, ~22% of explicit wall).

### MUMPS Error −9 ("work array too small") is RANK-COUNT-dependent and non-monotonic — and ICNTL14 strongly perturbs wall time, so hold it uniform across any np sweep
- **Bites:** parallel implicit runs (`system Mumps`) that work at np=1/2/8 can **fail at step 0 with INFO(1)=−9 at np=4** (or any specific rank count) — the workspace under-prediction depends on how the partition shapes the distributed factorization, not on problem size alone. Worse for benchmarking: the `-ICNTL14` value needed to survive also changes performance dramatically (measured: np=8 wall 90.5 s at `-ICNTL14 200` vs 33.0 s at `-ICNTL14 2000` on the same model), so a sweep with per-np ICNTL14 values produces meaningless comparisons.
- **Why:** ICNTL(14) is MUMPS's % workspace relaxation over its analysis-phase estimate; the estimate quality varies with the partition/ordering interaction. OpenSees defaults it to 20 (`OpenSeesCommands.cpp:4252`); the −9 handler prints the "make ICNTL14 larger" hint (`MumpsParallelSolver.cpp:168`). Default-20 failed at np=4 on a plain 5488-hex/18.9k-DOF block; 200 did NOT rescue it; 2000 did.
- **Workaround/status (2026-07-06, ADR-40 rank-3 MUMPS scaling measurement):** for any np sweep set `-ICNTL14` high (≥200, be ready for 2000) and **uniform across all rank counts**; verify every config actually completed (a DNF leaves stale output/h5 files from prior runs — delete artifacts before a fresh sweep). Measured sweep + table in [[40b_phase0_dominance_report]] §MUMPS addendum.

### `criticalTimeStep()` under SMS returns the PRE-scaling element pencil — the post-scaling effective limit is report-only with NO getter
- **Bites:** any driver/script that queries `criticalTimeStep()` on a mass-scaled run (CentralDifferenceSMS/SMSConsistent, ExplicitBathe `-sms`) expecting the stable-Δt after scaling. It gets the un-augmented sliver limit instead — `dt = safety*criticalTimeStep()` collapses Δt to the pre-scaling value and silently defeats SMS (conservative, so it "works", just wastes the entire scaling benefit).
- **Why:** `CentralDifferenceLadruno::getCriticalTimeStep()` (`CentralDifferenceLadruno.cpp:742–746`) and `ExplicitBathe::getCriticalTimeStep()` (`ExplicitBathe.cpp:1435–1441`) return the raw element-pencil value, which cannot see nodal injected mass (scaling writes `Node::setMass`; the pencil reads `ele->getMass()`). The post-scaling limit exists (`setSMSEffectiveLimit` ← `minDtSelfReport`, P3 #475) but is consumed ONLY by the `newStep()` "[PRE-SCALING estimate]" report — protected setter, no getter, no command plumbing.
- **Workaround/status (2026-07-02, ADR-65 adversarial gate):** treat the query as a conservative lower bound, or don't query under SMS at all (SMS runs are sized to `dtTarget` at construction — just use `dtTarget`). An SMS-aware adaptive driver needs a new `getSMSEffectiveLimit()` accessor first. See [[65_ladruno_explicit_dt_strategies_adr]] §Route B.

### Changing Δt between explicit `analyze()` calls without `revertToLastStep()` mis-centers the leap-frog by `(Δt_new−Δt_old)/2·a` per change
- **Bites:** any variable-Δt driver loop over `CentralDifferenceLadruno` (and the SMS subclasses). `newStep(dt_new)` after a step committed at `dt_old` advances `Vhalf += Δt_new·a_n` (`CentralDifferenceLadruno.cpp:574`) with no previous-Δt memory — but `v_{n−1/2}` is staggered for `Δt_old`, so the correct advance is `((Δt_old+Δt_new)/2)·a_n`. Each un-reseeded Δt change injects a systematic velocity error; small per change, compounding during ramps. No warning, nothing fails.
- **Why:** the leap-frog kernel is uniform-Δt by design; only the failure path was ever exercised with Δt changes, and that path happens to be correct because `revertToLastStep()` re-arms `firstStep` and rebuilds `v_{−1/2}` from committed state at the new Δt (`:420–436`, valid standalone per the `:413–419` comment).
- **Workaround/status (2026-07-02, ADR-65 adversarial gate):** call `revertToLastStep()` after the last commit before ANY `newStep` at a different Δt — grow or shrink, not just on failure. See [[65_ladruno_explicit_dt_strategies_adr]] §Route B (stagger reseed).

### Tangent-tracked `criticalTimeStep()` (`-recompute`/`-tangent`) OVER-reports the safe explicit step by √(E0/E_tan) for any reloadable material — unsafe as an adaptive-growth target
- **Bites:** any adaptive-Δt driver that grows the step toward `safety·criticalTimeStep()` on a model with plasticity or damage-with-reload. As the material softens/yields the tangent stiffness drops, so `criticalTimeStep()` returns a *larger* Δt_cr (Δt_cr = 2/ω_max ∝ √(m/k_tan)). Grow into it and the run blows up the instant an elastic-reload wave hits the softened element — because reload happens at the UNDAMAGED modulus E0, not the current tangent.
- **Why:** the pencil sees only the current (tangent) stiffness; a reloadable material can stiffen back to E0 in one step. Measured over-report = √(E0/E_tan): exactly 7.07× on Steel01 b=0.02 (=√(1/b)) and 12.5× on a Concrete01 softening bar (ADR-65 P0 oracle, `adr65_headroom_oracle/`, 2026-07-03). The tangent tracking itself is CORRECT and stays positive on the softening branch (0 bad samples / 2200 steps) — it's the *interpretation as a safe step* that's wrong.
- **Workaround/status (2026-07-03, ADR-65 Route B oracle):** an adaptive driver must clamp Δt growth to the UNDAMAGED/elastic Δt_cr, not the instantaneous tangent — which makes growth headroom ~1.00× (nil) for reloadable nonlinearity, i.e. adaptive *growth* buys nothing there. The safe use of the tangent query is the SHRINK direction (Δt_cr dropping mid-run) + permanent non-reloadable stiffness loss (erosion). See [[65_ladruno_explicit_dt_strategies_adr]] §Route B P0 oracle RUN.

### A held-load augmentation loop that re-enters `Domain::commit()` SILENTLY corrupts the recorder stream + `commitTag` — suppress them, don't try to undo them
- **Bites:** any "within-step" contact-ALM driver that drives `‖ḡ‖→augTol` at a fixed load by re-`analyze(1)`-ing under a zero-increment `LoadControl` (the shipped `analyze_augmented` proc, ADR-41 C2.2/C4/D1). Each held re-solve commits through `Domain::commit()`, which (a) fires EVERY recorder `theRecorders[i]->record(commitTag, currentTime)` and (b) bumps `commitTag++` — so N augmentation passes inject N SPURIOUS duplicate recorder samples at the SAME pseudo-time and jump `commitTag` by N. A model with a Node/Element recorder gets a corrupted output file (extra rows) the user never asked for, and downstream tooling keyed on `commitTag` mis-indexes.
- **Why:** `Domain::commit()` is the integrator-agnostic choke point — it is the SOLE recorder-firing path under a static `analyze(1)` (StaticAnalysis has no separate `theDomain->record()`; the record happens INSIDE commit). The contact Uzawa λ update is wired there too (`theContactDomain->commit()`), so you cannot skip the whole commit — you specifically need "commit the contact multiplier, do NOT record / do NOT advance commitTag". You also can't UNDO a recorder sample after the fact: `Recorder::record` may have already flushed to a file/stream, so a "snapshot + restore commitTag" protocol does not un-write the row — the only correct fix is to SUPPRESS the firing, not reverse it.
- **Workaround/status (2026-06-24, ADR-41 D1):** a `bool Domain::contactAugmenting` flag (default false; `set/isContactAugmenting`); `Domain::commit()` wraps ONLY the recorder loop + `commitTag++` in `if (!contactAugmenting) { … }` (everything else — node/element commit, `theContactDomain->commit()`, `committedTime=currentTime` — runs unchanged, so flag-off is byte-identical to stock). The `analyze_augmented` proc brackets its held-load sweep with `ladrunoBeginAugment`/`ladrunoEndAugment`, in a `try/finally` so an exception can't leave the flag stuck ON (recorders muted). Two things that make this SAFE-as-shipped but worth knowing: (1) `committedTime=currentTime` runs during augmentation but is a NO-OP because a zero-increment `LoadControl` holds `currentTime` constant — if a future caller ran a TIME-ADVANCING integrator under the flag, committedTime WOULD silently track forward with no recorder sample (the proc owns the `LoadControl(0.0)` switch, so this is unreachable today but is a latent trap if someone sets the flag by hand). (2) The held re-solves do NOT double-count the load because `LoadControl::newStep` does `applyLoadDomain(committedLambda + 0)` and `Domain::applyLoad` ASSIGNS the factor (zeroes then re-applies) — re-applying the same factor N times is idempotent, not Nλ. See [[48_ladruno_contact_capstone_adr]] contract #3 + the D1 row; vanilla edit in [[LEDGER_vanilla_files]].

### Adopting LadrunoProjection in a MULTI-SOLVE explicit integrator: `buildMass` can't ride the algorithm's solve — form a one-shot mass tangent in `newStep`
- **Bites:** making a Noh–Bathe integrator (`ExplicitBathe`/`ExplicitBatheLNVD`) a `LadrunoProjectionConsumer` (ADR-30 P5). The projector's `buildMass(theLinSOE)` reads `diag(M)` from the **Diagonal SOE** — but `DiagonalDirectSolver` overwrites `A[i]←1/A[i]` on its factor pass, so the mass must be read **before** the first solve. In `CentralDifferenceLadruno` the starter solve lives INSIDE `newStep`, so buildMass slots between its `formTangent` and `solve`. In `ExplicitBathe`/LNVD the FIRST solve is the **solution algorithm's** (after `newStep` returns) and the second is manual inside `update()` — the integrator never sees a "between formTangent and solve" window it controls.
- **Why:** the mass diagonal is transient: live as `m` only until the first `solve()` factors it to `1/m`, and the Bathe `update()` reuses that factored `1/m` (DiagonalDirectSolver factors once). By the time any integrator hook runs after the algorithm's solve, `diag(A)` is already `1/m`.
- **Workaround/status (2026-06-21):** in `newStep`, gate on `theProjector && !massBuilt` and do a **one-shot** `this->formTangent(CURRENT_TANGENT)` (re-assembles raw `M` into the SOE) → `theProjector->buildMass(theLinSOE)` → `massBuilt=true` → `project(*A_t)` (committed a0). The solution algorithm re-forms `M` before its own solve, so the extra assembly is harmless (one mass assembly on the first step of each stage). Gate on `!massBuilt` (re-armed by `domainChanged`/`setConstraintProjector`), **not** `firstStep` — `ExplicitBathe::domainChanged` does NOT reset `firstStep`, so a firstStep gate would skip the rebuild after staged construction. Also: the 4 SMS/Consistent subclasses (33009–33012) **inherit** projection for free because they chain `domainChanged()` to the base and don't override `newStep/update/commit` — only the two bases need editing. See [[30_ladruno_explicit_constraint_projection_adr]] / `projection_handler_handoff.md` §P5.

### A failed `OPS_GetIntInput`/`OPS_GetDoubleInput` consumes the arg in openseespy but NOT in the classic Tcl exe — blind `OPS_ResetCurrentInputArg(-1)` overshoots in Tcl
- **Bites:** any greedy "read ints until the next flag" parse scan that un-gets the
  flag with `OPS_ResetCurrentInputArg(-1)` after a failed numeric read. Symptom:
  the command parses fine from openseespy but is a parse error from the classic
  Tcl exe **only when the greedy option is followed by another flag** — e.g.
  `recorder ladruno f -G energy 1 -T nsteps 1` died with *"option -T nsteps
  requires a number-of-steps argument"* while `… -T nsteps 1 -G energy 1` ran.
- **Why:** the two API stacks disagree on cursor consumption at parse failure.
  openseespy (`PythonModule.cpp::getInt`) calls `incrCurrentArg()` *before* the
  type check, so the failed token is consumed and `-1` restores it — net zero.
  The classic Tcl exe (`SRC/api/elementAPI_TCL.cpp::OPS_GetIntInput`, ~line 224)
  returns `-1` *without* advancing, so the `-1` rewind steps back onto the
  previous already-consumed token and the whole tail parse shifts by one
  (a local `numdata` countdown then starves the line's last token). String reads
  (`OPS_GetString`/`OPS_GetStringFromAll`) consume in both stacks — only the
  numeric readers diverge.
- **Workaround/status (2026-06-11):** never blind-rewind — capture
  `int oldn = OPS_GetNumRemainingInputArgs();` before the numeric read and rewind
  only `if (OPS_GetNumRemainingInputArgs() < oldn)`. Applied to the
  LadrunoRecorder `-G` scan, both LadrunoMonitorRecorder scans, and the
  LadrunoBrick `-hourglass` coeff probe (LadrunoRCConcrete's `readList` already
  did it right). Guarded by the `TCL FLAG ORDER` regression gate
  (`flag_order_model.tcl`); the openseespy ordering was already covered by
  `energy_model.py`.

### Cloning the ASDConcrete3D plastic-damage spine: two silent-but-fatal requirements (`/E` measure + E-consistent backbone `q`)
- **Bites:** re-implementing the `ASDConcrete3D` update (e.g. `LadrunoRCKernel.h`). Two omissions each yield a material that *compiles, runs, and even passes a β-ratio test* yet is physically wrong:
  1. **The equivalent-strain measures must be divided by E.** `equivalentTensile/CompressiveStrainMeasure` return `lublinerCriterion(...) / E` (`ASDConcrete3DMaterial.cpp:2509,2522`) — the Lubliner criterion is a STRESS; `/E` converts it to the *strain* abscissa that indexes the hardening backbone. Omit it and the abscissa is ~E× too large → the lookup lands in deep softening → `dt̄/dc̄≈1` → nominal stress collapses to ~0 (looks like near-zero stiffness). A β-*ratio* test CANNOT catch this (the `/E` scaling cancels in the ratio); only an ABSOLUTE-stress test (σ=E·ε) does.
  2. **The effective-stress backbone `q` is E-consistent BY CONSTRUCTION — it is NOT `q=y/(1−d)` on the raw user points.** `HardeningLaw` c-tor + `adjust()` (`ASDConcrete3DMaterial.cpp:869-939,1134-1180`) prepend `(0,0,0)`, force `d[0]=0`, force the first segment elastic (`y[1]=E·x[1]`), **cap every secant slope at E**, enforce monotone plastic strain + non-decreasing damage, and ONLY THEN set `q=y/(1−d)`. Compute `q=y/(1−d)` on the raw points and any segment with effective slope > E gives `q>E·x` ⇒ `dc_plastic=1−q/(E·Δx) < 0` ⇒ NEGATIVE damage ⇒ the *elastic* stress is amplified (we saw exactly 8/7× too stiff). The tension branch can pass by luck if its slope already equals E.
- **Why:** the hardening abscissa is a strain; `q` is the undamaged (plastic-only) effective stress whose elastic branch must have slope exactly E so the damage split `d=1−y/q` starts at 0.
- **Workaround/status:** both replicated in `LadrunoRCKernel.h` (`/E` at the measure call site; `buildBackbone()` mirrors the c-tor+`adjust()`); proven by the reduce-to-ASDConcrete3D Zone-A gate (tension **and** compression byte-match). Also: a single OpenSees brick with ALL 24 DOFs prescribed (homogeneous-strain probe) **segfaults the `Transformation` constraint handler** (0 free equations) — use `constraints Penalty 1e14 1e14` instead. Learned 2026-06-03 building [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] (PR #155).

### Crack-band materials read element size via `ops_TheActiveElement->getCharacteristicLength()` — and the base default is wrong for high-order elements
- **Bites:** `ASDConcrete3DMaterial` (and any crack-band/smeared-crack material) auto-
  regularizes its softening branch by the element characteristic length `lch`,
  read **once** on the first `setTrialStrain` via the global
  `ops_TheActiveElement->getCharacteristicLength()`
  (`ASDConcrete3DMaterial.cpp:1614`, guarded by `regularization_done`,
  `ASDConcrete3DMaterial.h:386`). If your element doesn't override
  `getCharacteristicLength()`, it inherits `Element::getCharacteristicLength()`
  (`SRC/element/Element.cpp:682`), which returns the **minimum inter-node
  distance**. On a quadratic element (BezierTri6, BezierTet10, TenNodeTetrahedron,
  the quad/brick "...N" siblings) the closest node pair is corner-to-mid-edge,
  i.e. ≈ **½** the true edge length → `lch` under-estimated ~2× → fracture energy
  smeared over too small a band → **over-softening / spurious snap-back**.
- **Why:** the global is set by the framework, not the element — `Domain::addElement`
  (`Domain.cpp:447`) and, crucially, the `Domain::update()` loop
  (`Domain.cpp:2263`) sets `ops_TheActiveElement = theEle` immediately before
  `theEle->update()`. So as long as the element pushes strain **only inside
  `update()`** (both Bézier elements do), the active pointer is correct when the
  one-time regularization fires — no wrong-element window. The value is just
  geometrically poor for the min-distance default on high-order elements.
- **Also:** `getCharacteristicLength()` returns a **single element scalar**, read
  once per material instance — so a true per-Gauss-point `lch = sqrt(detJ·w)` is
  *not* expressible through this seam; every GP material on the element gets the
  same value. Pick one representative element size.
- **Workaround/status (2026-06-01):** override `getCharacteristicLength()` with an
  element-size equivalent from the integrated area/volume — BezierTri6 returns
  `sqrt(2·A)` (right-isosceles-triangle edge), BezierTet10 returns `cbrt(6·V)`
  (right-tetrahedron leg); both recover the true edge length on right-simplex
  meshes and err in the safe (under-estimating) direction on curved Bézier edges.
  LadrunoBrick is exempt — its 8 corner nodes give the correct min-edge already,
  and its only strain-push site is `update()` (verified read-only assembly).
  Done in `BezierTri6.cpp` / `BezierTet10.cpp` `getCharacteristicLength()`.

### A new element needs registering in THREE dispatch sites, not one — the standalone Tcl `OpenSees.exe` uses a different table than OpenSeesPy
- **Bites:** you add a new `element` (e.g. `LadrunoDispBeamColumn`), wire `classTags.h` + `FEM_ObjectBrokerAllClasses.cpp` + the `functionMap` in `SRC/interpreter/OpenSeesElementCommands.cpp`, it **builds and links clean**, works in OpenSeesPy — but the standalone `OpenSees.exe` (Tcl) reports `ERROR -- element of type X not known`.
- **Why:** the Tcl `element` command is `TclModelBuilder_addElement` in `SRC/element/TclElementCommands.cpp`, which dispatches through its OWN `ladrunoElementTable` / built-in tables — NOT the `OpenSeesElementCommands.cpp` `functionMap` (that serves the runtime/Python `OPS_Element()` path) and NOT `getLibraryFunction` (that does `LoadLibrary("Name.dll")`, which returns -1 for a built-in, so it never finds an in-exe symbol). The two error strings are nearly identical (`"element of type X not known"` in both `TclElementCommands.cpp:2195` and `runtime/commands/modeling/element.cpp:491`), so the message doesn't tell you which path failed.
- **Workaround/status (2026-06-16):** register a new Ladruno element in **all** of: (1) `SRC/classTags.h`, (2) `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` (include + `case`), (3) `SRC/interpreter/OpenSeesElementCommands.cpp` `functionMap` (OpenSeesPy), (4) `SRC/element/TclElementCommands.cpp` — both the `extern void *OPS_X(void);` block (~line 108) and the `ladrunoElementTable` (~line 582). Grep the WORKING sibling across `SRC/` to enumerate sites: `grep -rln "OPS_LadrunoIMKBeam" SRC/ | grep -v ladrunoIMKBeam/` lists exactly the registration files to mirror. See [[32_ladruno_dispbeamcolumn_regularization_adr]].

### zeroLength ignores stiffness-proportional Rayleigh unless `-doRayleigh 1`
- **Bites:** a `zeroLength` / `zeroLengthSection` element contributes **zero**
  stiffness-proportional Rayleigh damping (`betaK`, `betaKinit`/`betaK0`,
  `betaKcomm`/`betaKc`) by default. You set `rayleigh 0 0 0.0159 0`, expect
  ζ≈0.05, and measure ζ≈0. Mass-proportional `alphaM` (which lives on the node,
  not the element) works regardless, which masks the problem.
- **Why:** the element carries an internal `doRayleigh` flag, default **0**, that
  gates whether `getDamp()`/`getResistingForceIncInertia()` include the element's
  stiffness term. The `-doRayleigh` option flips it: `element zeroLength … -dir 1
  -doRayleigh 1`. Most other elements default the flag on; zeroLength does not.
- **Workaround/status (2026-06-01):** pass `-doRayleigh 1` whenever you want
  stiffness-proportional Rayleigh on a zeroLength; or (better) model the damping
  physically with a `Viscous`/`ViscousDamper` uniaxial material on the DOF — that
  enters R(u̇) directly and is explicit-safe. Pinned by
  `tests/test_damping_channels.py::test_zeroLength_doRayleigh_default_off`
  (ζ≈0 with default) and `::test_betaK0_realises_target_zeta` (ζ=0.05 with flag).
  Full map: [[12_damping_channels]].

### `modalDampingQ` (force-only modal damping) applies damping with the WRONG SIGN
- **Bites:** `modalDampingQ ζ` on a free-vibration SDOF/MDOF *amplifies* the
  response instead of damping it — measured ζ comes out ≈ **−ζ_target**.
  `modalDamping ζ` (the matrix form) is correct (+ζ).
- **Why (evidence, not yet line-pinned):** the sign error is **Δt-independent** —
  refining the step (steps/period 100→800) converges Q to −0.04996…→−0.04998 and
  the matrix form to +0.05000. A lag/explicit-stability artifact would vanish on
  refinement; this doesn't, so it's a **structural sign inconsistency**, not a
  numerical one. The only live force routine is the M-weighted
  `IncrementalIntegrator::addModalDampingForce` (`IncrementalIntegrator.cpp:502`,
  `setB` at :556); the two earlier variants (:303, :349) are commented out. Both
  `modalDamping` and `modalDampingQ` add that SAME `−Cv` force in `formUnbalance`
  (`TransientIntegrator.cpp:135`); the ONLY difference is `modalDamping` also adds
  `+c·C` to the tangent (`addModalDampingMatrix` :563, `formTangent` :88). So the
  matrix term is somehow compensating a force that, alone, has the wrong net sign.
- **PURE SIGN INVERSION — confirmed:** `modalDampingQ(-0.05)` damps correctly at
  **+0.05003**. So the force-only path injects `+Cv` energy instead of dissipating
  `-Cv`. The residual sign convention itself is fine (`FE_Element::addRIncInertia
  ToResidual` does `theResidual += -1.0·getResistingForceIncInertia`,
  `FE_Element.cpp:517`; the modal `-Cv` from `addModalDampingForce`/`setB` matches
  it). So the inversion is NOT in the force expression — it's that a velocity-
  proportional force applied through the *residual only* (no `+c·C` in the tangent,
  the term `modalDamping` adds and `modalDampingQ` omits) is integrated by Newmark
  with the opposite effective sign. i.e. `modalDampingQ` as a standalone force-only
  mode appears to have never worked; only `modalDamping` (matrix + force together)
  is correct.
- **Workaround/status (2026-06-01):** **use `modalDamping`, never `modalDampingQ`.**
  Reproduced under both `Newton` and `Linear`, Δt-independent. Pinned as a `strict`
  xfail: `tests/test_damping_channels.py::test_modalDampingQ_force_only_matches_matrix`.
  Genuine upstream bug. **DECISION (2026-06-01): DOCUMENT-ONLY** — not patched, not reported
  upstream; `modalDamping` (matrix) covers the use case. The `strict` xfail auto-detects any
  future upstream fix (would start XPASSing). A naive "flip the force sign" would break
  `modalDamping` (shares the same force), so any fix must special-case the `inclMatrix==false`
  path. Full map: [[12_damping_channels]].

### Assumed-strain hourglass: the dev-projection vs reduced-shear interaction is nu-coupled
- **Bites:** building the LadrunoBrick `physical` hourglass straight from Belytschko
  eq 8.7.26 (pointwise-*isochoric* assumed strain: 2/3,-1/3 dev-projection on the
  normal hourglass strains + mode-subset reduced shear) gave an element that was
  ~75% **too soft** in bending and got *worse* with nu (ratio 1.73 @nu=0 -> 2.59
  @nu=0.499 vs the analytic 1.0). Patch test + rank were exact, so the bug hid
  from the usual gates.
- **Why:** the dev-projection IS the proper B-bar mean-dilatation treatment for
  the hourglass normals (the algebra collapses to the same 2/3,-1/3). On its own
  (with full shear) it's fine; combined with the **reduced** assumed shear it
  removes too much energy -> over-soft, and the error scales with lambda(nu).
  Dropping the dev-projection (FULL compatible normal strains + reduced shear)
  gives a **correct shear-locking cure** — matches OpenSees `SSPbrick` to 3 digits
  and converges (0.94->1.005, nx=2..32) at nu=0 — but then VOLUMETRIC-locks at
  nu->0.5. There is **no single static projection** correct across nu with this
  shear field; a general all-nu element needs the coupled SSP/ASQBI operator
  (Belytschko sec 8.7.8 explicitly: 3D assumed-strain structure "not fully developed").
- **The validating oracle:** patch + rank CANNOT validate an assumed-strain
  element (gamma-orthogonality makes any variant pass). Use a **bending-convergence
  benchmark** and cross-check against `SSPbrick` (a proven OpenSees assumed-strain
  hex, ~1.0 across all nu). `tests/test_ladrunoBrick_bending.py`.
- **Status (2026-06-01):** shipped `physical` as the FULL-normals + reduced-shear
  **shear-locking cure** (verified vs SSPbrick); documented that near-incompressible
  needs `-formulation bbar`. A coupled general-nu operator is future work.
- **The definitive difference vs SSPbrick (read its source).** `SSPbrick.cpp` is an
  **EAS element = bbar + statically-condensed enhanced strain**. (1) Volumetric:
  its constant `Bnot` uses `dNmod` = mean-dilatation (B-bar) modified gradients
  (`SSPbrick.cpp:1254,1266`). (2) Shear/bending: 9 internal **enhanced-strain
  modes** `Fe`, condensed out — `interior = FCF − K_uα·K_αα⁻¹·K_αu`
  (`SSPbrick.cpp:1968`), then `Kstab = Mbenᵀ·interior·Mben`. The **static
  condensation** is why SSP works for ALL nu: the internal modes *adapt to C*. My
  `physical` is a single FIXED assumed-strain B (no internal DOFs/condensation) →
  can cure shear OR volumetric, never both across nu. **Upshot: a general-nu
  "physical" = our reserved `eas` formulation (v2), and SSPbrick is the production
  blueprint (bbar constant part + condensed EAS).** See `SSPbrick.cpp:1053`
  (`GetStab`), `:1243` (G/gamma), `:1647` (enhanced-strain block), `:1968`
  (condensation).
- **SHIPPED (2026-06-01): `LadrunoBrick -formulation eas` is the SSPbrick port.**
  Confirmed while porting: SSPbrick condenses the enhanced modes with the
  **initial** tangent (`GetStab` is called once in `setDomain`), so `Bnot`/`Kstab`
  are **constant** — there is **no per-step α internal state**, contrary to the
  general textbook EAS picture. That collapses the "heavy bit" (no
  `commitState`/`sendSelf` of α): the operators are deterministic from geometry +
  C(0), so the parallel receive side just rebuilds them in `setDomain` and
  `sendSelf` ships nothing extra. Validation gate: for a linear-elastic material
  the assembled `eas` stiffness is *identical* to `SSPbrick`, so the bending-
  benchmark tip matches SSPbrick to ~1e-6 across ν∈{0,0.3,0.45,0.499} (where
  `physical` vol-locks). One caveat: SSPbrick itself sends `Bnot`/`Kstab`/`J[]`
  over `sendSelf` (its null-ctor sets `mInitialize=false` → skips `GetStab` on
  recv); LadrunoBrick instead always rebuilds in `setDomain` — simpler, same
  result, smaller stream.

### `BbarBrick` has no `update()` — a bare `eleResponse("stresses")` reads the *predictor* (u=0) state
- **Bites:** after a static/linear solve, `ops.eleResponse(tag, "stresses")` on an
  upstream `bbarBrick` returns **all zeros** even though the displacements are
  correct. A regression that compares element stresses against `bbarBrick` then
  fails with `ours=<real> vs ref=0`. (`stdBrick` does **not** show this — it reads
  back real stresses.)
- **Why:** `Brick` overrides `Element::update()` to push the material trial strain
  every step, so its committed `getStress()` reflects the solved displacement.
  `BbarBrick` has **no `update()`** — it calls `setTrialStrain` lazily, only inside
  `formResidAndTangent`. The `"stresses"` response (responseID 3) just returns
  `materialPointers[i]->getStress()` *without* recomputing, so it reflects whatever
  the last `formResidAndTangent` saw — which, after a linear step, is the predictor
  state (u=0 → zero strain → zero stress).
- **Workaround:** read `"forces"` (responseID 1 → `getResistingForce` →
  `formResidAndTangent` → `setTrialStrain` at the committed u) **before** reading
  `"stresses"`. Then both lazy- and eager-strain elements report the solved state.
  `LadrunoBrick` implements `update()` (like `Brick`), so it is order-insensitive —
  this only matters when comparing against `bbarBrick`.
- **How it surfaced:** `tests/test_ladrunoBrick_element.py::test_bbar_matches_bbarBrick`
  (LadrunoBrick `bbar`↔`bbarBrick` regression, PR [#65](https://github.com/nmorabowen/OpenSees/pull/65)).
  Displacements matched to 1e-9 (kernel-equivalent); only the stress readback timing
  differed. Date learned: 2026-06-01.

### `OPS_GetString` returns the literal `"Invalid String Input!"` for a Python NUMERIC arg — never use it to peek an optional numeric token
- **Bites:** a factory that peeks an optional trailing numeric arg by calling
  `OPS_GetString()` and parsing it (e.g. `strtod`) works under Tcl (args arrive as
  strings, `"0.05"` parses) but **silently drops the value under openseespy**. For a
  Python number (`ops.element(..., '-hourglass','viscous', 0.05)`), `0.05` is a
  `PyFloat`, not a `PyUnicode`. `PythonModule::getString()` (`PyUnicode_Check` fails)
  returns NULL → `OPS_GetString` (`OpenSeesCommands.cpp`) substitutes the literal
  string `"Invalid String Input!"`. Your `strtod` then fails, you push the token
  back, the option loop re-reads it (NULL again → unknown-option WARNING), and the
  coefficient is **consumed-and-discarded** — the element silently falls back to the
  default coeff. Tcl path is unaffected, so it hides from Tcl-only testing.
- **Why:** `OPS_GetString` is for *string* options; `OPS_GetDoubleInput` /
  `OPS_GetIntInput` are the number-aware readers (their Python backends accept
  `PyFloat`/`PyLong`/`PyBool`). 
- **Fix / idiom for an OPTIONAL trailing numeric arg that may instead be the next
  flag:** read it with `OPS_GetDoubleInput(&n1,&tmp)`; on success use it, on failure
  `OPS_ResetCurrentInputArg(-1)` to un-get it for the option loop. `OPS_GetDoubleInput`
  advances the arg cursor by one even when the conversion fails, so the single
  `-1` reset un-gets exactly that token — works on **both** Tcl and Python.
- **How it surfaced:** adversarial review of LadrunoBrick PR [#75](https://github.com/nmorabowen/OpenSees/pull/75)
  (eas + viscous). The first `-hourglass <type> -lumped` fix used `OPS_GetString`+`strtod`;
  caught before merge by adding `test_hourglass_coefficient_reaches_kernel` (a numeric
  Python-float coeff must change the response). Date learned: 2026-06-01.

### `Truss` (and most elements) default to a LUMPED mass — `-lump diagonal` ≡ `-lump rowsum` on them
- **Bites:** the ExplicitBathe / CriticalTimeStep `-lump <rowsum|diagonal>` option
  only changes `dt_cr` when the element's `getMass()` returns a **consistent**
  (non-diagonal) matrix. The `Truss` element defaults to `cMass=0` (lumped:
  `0.5*rho*L` on the diagonal, zero off-diagonals). For a diagonal matrix,
  "diagonal-of-consistent" (take `M_ii`) and "row-sum" (sum each row onto the
  diagonal) are identical, so both lumpings return the same `dt_cr` (= `le/c` for
  a 2-node bar). A test that asserts `dt_diag < dt_rowsum` on a default lumped bar
  is **vacuous and fails** (the two tie).
- **Why:** `Truss::getMass()` returns lumped mass unless built with `-cMass 1`;
  the consistent branch is `rho*A*Le/6 * [[2,1],[1,2]]`. Only then do the lumpings
  diverge: diagonal-of-consistent = `rho*A*Le/3`, row-sum = `rho*A*Le/2`, giving
  `dt_diag/dt_rowsum = sqrt((Le/3)/(Le/2)) = sqrt(2/3) ≈ 0.816` (per-element pencil).
- **How it surfaced:** `tests/test_explicitBathe_integrator.py::test_eb_lump_diagonal`
  failed on Zone-A (Ubuntu) CI with `diagonal=0.02000 vs rowsum=0.02000` (a tie at
  `le/c`). The Linux result is **correct**; the source `_verify_explicit.py` test 11
  had the same faulty premise (built a lumped bar yet asserted the inequality).
- **Workaround/status:** build the bar with `-cMass 1` (consistent mass) for any
  test that means to exercise the diagonal-vs-rowsum *difference*. Fixed in the
  Zone-A port (PR #40) and mirrored back into `Ladruno_scripts/_verify_explicit.py`.
  *(Learned 2026-05-31, Zone-A ExplicitBathe battery.)*

### MPCORecorder `exit(-1)` kills the kernel
- **Bites:** ~25 raw libc `exit(-1)` calls in upstream `MPCORecorder.cpp` hard-kill
  the Jupyter kernel on any recorder error and leave the `.mpco` unflushed.
- **Why:** upstream uses `exit()` for error handling; no exception path. The file
  has been byte-identical to `master` and frozen ~4.5 years.
- **Status:** diagnose-only as of 2026-05-18. A local `ladruno` patch is the only
  fix; deferred. The Ladruno recorder rewrite (`ladruno`) avoids the pattern.

### `patch_banner.py` path went stale after workspace consolidation
- **Bites:** `patch_banner.py` computed `TARGET = ROOT/OpenSees/SRC/...`, assuming an
  inner `OpenSees/` dir. After the 2026-05-30 consolidation `SRC/` sits directly
  under the repo root, so the script silently targeted a non-existent path.
- **Why:** the script predated moving the workspace into the repo root.
- **Status:** **fixed** — paths recomputed to `ROOT/SRC/...`; the script now also
  patches the feature list into both `tclMain.cpp` and `PythonModule.cpp`.

### MPCO is nodal-blind without an element-result parity gate
- **Bites:** a recorder can pass all *nodal* parity checks yet silently drop or
  mis-tag *element/material* results (e.g. Lagrange quad/tri `setResponse` tags).
- **Why:** nodal and element result paths are independent; testing one doesn't
  cover the other.
- **Status:** addressed by the element-result parity gate ([#14](https://github.com/nmorabowen/OpenSees/pull/14)) and the upstream `setResponse` fix ([#7](https://github.com/nmorabowen/OpenSees/pull/7)).

### Damping shrinks the explicit critical time step — and βK is a trap
- **Bites:** adding damping *reduces* `dt_cr` in both coupled and explicit modes;
  stiffness-proportional damping (`βK`) is especially punishing under explicit
  integration.
- **Why:** βK inflates the effective high-frequency content that sets `dt_cr`.
- **Status:** use mass-proportional (`αM`) damping in explicit; captured in the
  CentralDifferenceLadruno plan (`project_robust_central_difference`).

### `Node::getMass()` returns ONLY the nodal `mass` command value — element-density mass is invisible there
- **Bites:** any code that reasons about a node's mass by reading `Node::getMass()` —
  e.g. the ADR-39 B1 SOFT=1 contact penalty, which needs the gap-mode effective mass
  `m_eff` from the mass the explicit integrator actually inverts. For a model whose mass
  comes from ELEMENT density (a solid `LadrunoBrick`/truss with `-rho`, no per-node
  `mass`), `Node::getMass()` returns a **zeroed** matrix (`Node.cpp:1214` — `mass==0 ⇒
  return a zero matrix`). The first B1 cut sized `m_eff` from `Node::getMass()` and so saw
  `m=0` for every solid-body contact node → silently fell back to the stiff base kn → the
  exact divergence SOFT exists to prevent. Tests that set `ops.mass(...)` directly never
  catch it (found by the B1 adversarial code gate, MAJOR).
- **Why:** OpenSees keeps nodal mass (the `mass` command, stored on the `Node`) and
  element mass (each element's `getMass()`) as **separate** contributions. The global mass
  is assembled from BOTH: the integrator's `formNodTangent → DOF_Group::addMtoTang →
  Node::getMass()` AND `formEleTangent → Element::addMtoTang() → Element::getMass()`. There
  is no `Node::addMass` that folds element mass back onto the node. So `Node::getMass()` is
  the nodal-`mass`-only piece, never the assembled total.
- **Workaround/status (2026-06-24):** to get the mass the explicit solve inverts (the
  assembled global diagonal, for `system Diagonal`), reconstruct it: `m[d] = nodal
  mass(d) + Σ_elements diag(M_e) at that node's translational DOFs`. B1's handler does this
  once per `handle()` (`ladrunoBuildNodalMass` in `LadrunoContactHandler.cpp`) and caches it
  on `LadrunoContactDomain` for the stateless adapter. Matches `diag(global M)` for `system
  Diagonal` (the `DiagonalSOE` default extracts `M(i,i)`, no row-sum); a row-sum/lumped
  distributed SOE (OpenSeesMP `MPIDiagonal`) would differ — soft contact is serial-only
  today. Translation-first DOF order (`[u | θ]`) makes a node's first 3 element DOFs
  translational.

### An explicit SOFT/Courant penalty `k ∝ 1/dt²` keeps the contact event at a CONSTANT step count — energy error converges in SOFSCL, not dt
- **Bites:** validating the ADR-39 B1 SOFT=1 penalty's energy balance. The instinct is to
  refine `dt` and watch the impact restitution `e → 1` (energy conservation). It does NOT
  improve — `1−e` is flat across `dt`. The naive read is "the penalty leaks energy."
- **Why:** SOFT sizes `k_soft = SOFSCL·4·m_eff/dt²`, so the contact period
  `T_contact = 2π√(m_eff/k_soft) ∝ dt`. The steps spanning a contact event,
  `T_contact/dt = π/√SOFSCL`, are **independent of dt** — refining dt shrinks the period
  and the step in lockstep, so the contact is always resolved by the same ~`π/√SOFSCL`
  steps. The discrete one-sided-contact engagement error (EITHER sign at coarse
  resolution — the chatter D2 viscous damps) is therefore a function of **SOFSCL**, not dt.
- **Workaround/status (2026-06-24):** test energy convergence by refining **SOFSCL** (not
  dt): `proto_b1_soft_penalty.py` T3/T4 sweep SOFSCL∈{0.1, 0.025, 0.00625} and show
  `|1−e|`, `|ΔKE/KE₀|` → 0. Frame it as "bounded & SOFSCL-convergent (not a formulation
  leak)" — NOT "conservative at the shipped SOFSCL=0.1" (there the bounded error is ~1–2%,
  sign-indefinite; that's the chatter, by design left for `-visc`). SOFSCL is the
  accuracy/stability knob: smaller = stiffer + better-resolved + less penetration.

### SOFT=2 (segment-based explicit penalty): the per-node Courant bound is necessary-NOT-sufficient — inter-node coupling in the assembled K_c raises ω_max·dt ~2×
- **Bites:** sizing the ADR-39 B2 SOFT=2 segment-based penalty (`contact … -mortar -soft <SOFSCL>`).
  Each slave node I gets a per-node Courant penalty `k_soft,I = SOFSCL·4·m_eff,I/dt²` so that node's
  contact mode `ω_I·dt = 2√SOFSCL` (≤ 2 for SOFSCL ≤ 1, exactly like the B1 NTS rule). The instinct
  is therefore "SOFSCL up to ~1 is stable." It is NOT for the segment lane: the ASSEMBLED contact
  stiffness `K_c = Σ_I k_soft,I·B_Iᵀ B_I` couples the slave + master facet nodes that multiple `p_I`
  springs SHARE (via the dense D/M mortar matrices), so the max eigenvalue of `M⁻¹K_c` is larger than
  any single `ω_I²`. On the oracle's non-matched facet pair the coupled `ω_max·dt = 1.19` at SOFSCL=0.1
  (vs the per-node `2√0.1 = 0.63`) — a ~1.9× amplification — and a SOFSCL sweep crosses the
  central-difference limit `ω_max·dt = 2` near **SOFSCL ≈ 0.3**, not 1.0.
- **Why:** the B1 NTS gap operator `B = [n | −Nᵢ n]` couples one slave node to one segment's nodes; at
  a shared segment node the multiple pairs' springs add, but the coupling is mild. The B2 segment
  operator `B_I = [D,−M]/a_I` is DENSE (every slave node of a facet couples to every master node it
  overlaps), so the assembled Gram `Σ k_I B_Iᵀ B_I` concentrates much more on the shared modes. The
  per-node sizing bounds each node's OWN mode, not the coupled spectrum.
- **Workaround/status (2026-06-24):** keep the **default SOFSCL = 0.10** (comfortable margin — coupled
  `ω·dt ≈ 1.2`). The command surface warns when `-mortar -soft SOFSCL > 0.25` (the per-node `> 1`
  warning understates the segment lane); a denser / more-coupled mesh can push the amplification
  higher, so a user raising SOFSCL must verify a dt margin (or drop dt). Quantified by
  `proto_b2_soft2_segment.py` T4 (the assembled-`K_c` eigenbound). This is the segment-lane analogue
  of B1's per-pair-conservatism caveat ([[#402](https://github.com/nmorabowen/OpenSees/pull/402)]).
  Also inherited from B1 (not B2-specific): `ladrunoBuildNodalMass` reconstructs `m_eff` from the
  DIAGONAL-of-consistent element mass; under `CentralDifferenceLadruno -lump rowsum|hrz` the integrator
  inverts a different lumping, modestly eroding the SOFSCL margin — diagonal/HRZ is fine, rowsum drifts.

### `partition` needs the METIS 5 API (Patch 9)
- **Bites:** the openseespy `partition` command fails / mislinks without the
  METIS 5 API path.
- **Why:** upstream assumes an older METIS; the MP build links METIS 5.
- **Status:** guarded with a clear error message in `OpenSeesMiscCommands.cpp`
  ([#1](https://github.com/nmorabowen/OpenSees/pull/1)); `OPS_HAVE_METIS5` follow-up noted.

### HDF5 won't create a group under a dataset — the "writeSections" red herring
- **Bites:** `ladruno` emitted ~3 non-fatal `HDF5-DIAG` errors to stderr for
  every custom-rule element (Lobatto beam, MVLEM, BezierTri6). The error
  signature (`invalid location identifier` / `dset_id is not a dataset ID`) was
  misattributed to `writeSections` for a while.
- **Why:** `writeModelElements` made `MODEL/ELEMENTS/<name>` the CONNECTIVITY
  **dataset**, then tried to `H5Gcreate` a `QUADRATURE` child **group** under it.
  HDF5 datasets cannot parent groups, so the handle was invalid and the
  `GP_PARAM`/`GP_WEIGHT` writes under it failed — silently losing the quadrature
  group (data was otherwise intact; GP_X attr carried the coords). The opserr
  buffering made the errors *look* like they came from `writeSections`, which is
  clean. Reliably-ordered `fprintf(stderr)` pinned them to `writeModelElements`.
- **Status:** **fixed.** [#16](https://github.com/nmorabowen/OpenSees/pull/16)
  removed the broken block as a stopgap (GP_X-only); [#18](https://github.com/nmorabowen/OpenSees/pull/18)
  did the real schema-v1 fix — `<name>` is now a GROUP holding `CONNECTIVITY` +
  `QUADRATURE`/{`GP_PARAM`,`GP_WEIGHT`}. See [[LEDGER_implementations]] Ladruno row.

### `Response`/`Information` stores a Matrix and a Vector separately — `getData()` only returns the Vector
- **Bites:** an element `setResponse` that returns a `Matrix` (via `setMatrix`,
  e.g. BezierTri6 `integrationPoints` → `Matrix(nGP,2)`) is invisible to code
  reading `response->getInformation().getData()` (that returns the **Vector**
  slot, which is empty/unrelated). The recorder's legacy `getCustomGaussPointLocations`
  read the Vector and silently missed the matrix-typed barycentric GP coords.
- **Why:** `Information` is a tagged union (`theType` ∈ {…, VectorType, MatrixType}
  with separate `theVector`/`theMatrix` pointers); `getData()` hardcodes the Vector.
- **Status:** **handled** — check `info.theType == MatrixType` and read `*info.theMatrix`
  for multi-dim parametric rules ([#18](https://github.com/nmorabowen/OpenSees/pull/18)).

### Gauss-point ordering is per-element, NOT a standard tensor order
- **Bites:** there is no global "Gauss point #k → natural coords" convention across
  OpenSees solid elements. Deriving GP coordinates from a `(geometry, rule)` table and
  assuming the usual lexicographic tensor order silently mis-maps results (stress at
  GP `k` paired with the wrong location) — invisible to nodal parity and to any
  value-only check (the frozen MPCO recorder never writes GP coords; STKO holds a
  hardcoded per-rule table reader-side).
- **Why:** each element hardcodes its own integration-point loop. Verified from sources:
  `FourNodeQuad` (4-pt 2×2) walks **counter-clockwise** `(−,−),(+,−),(+,+),(−,+)`
  (FourNodeQuad.cpp:298-305) — *not* the lexicographic `(−,−),(+,−),(−,+),(+,+)`;
  `Brick`/stdBrick (8-pt) walks nested `for i{for j{for k}}` ⇒ x-outer..z-inner
  **lexicographic** (Brick.cpp:536-542); `NineNodeQuad` (9-pt) walks 4 CCW corners →
  4 CCW edge-mids → center (NineNodeQuad.cpp:127-144), a serendipity-style order, not
  a 3×3 tensor sweep. So even two "quad GL" elements need not share GP order.
- **Status:** for `Ladruno`, the standard-quadrature `GP_PARAM[k]` MUST equal the
  element's own k-th GP natural coords (so it pairs with result `gauss_id k`); the
  table is verified against each canonical element's source, and the belt-and-suspenders
  `GLOBAL_GP_COORDS` round-trip oracle (`x(GP_PARAM[k])` vs the C++-computed global GP)
  catches any ordering/basis mismatch at write time. Learned 2026-05-30 building the
  standard-rule QUADRATURE table.

### `FourNodeTetrahedron` has 1 Gauss point, not 4 (the recorder comment lies)
- **Bites:** the recorder's `getGeometryAndIntRuleByClassTag` comment reads
  "4-node tetrahedron with **4 gp**" and maps it to `Tetrahedron_GaussLegendre_1`; an
  implementer trusting that would write 4 GP coordinates for a 1-GP element.
- **Why:** `FourNodeTetrahedron.cpp` defines `sg[]={0.25}` (a single abscissa) and its
  Gauss loop is collapsed: `i = j = k = 0; // Just one Gauss point in a tet`
  (FourNodeTetrahedron.cpp:226,574). A linear tet integrates exactly with one centroid
  point `(¼,¼,¼)` in barycentric coords; the "4 gp" comment is stale/wrong.
- **Status:** table uses **1 GP** for the tet rule. Authoritative GP count comes from
  the element source / the OutputDescriptor response walk, never the recorder comment.
  Learned 2026-05-30.

### `openseessp` is an unbuilt subsystem (Python has no SP)
- **Bites:** expecting an `import openseessp` analogous to `openseesmp`.
- **Why:** the SP parallel engine exists only on the Tcl side; the Python engine
  never had an SP build. Its build trace was fully removed from the fork.
- **Status:** by design — use `openseesmp`. Rationale in `02_openseespymp.md`.

### A `fix`/`sp` added mid-analysis snaps the node to the REFERENCE frame, not the current deformed shape
- **Bites:** in a staged analysis (deform under stage 1, hold load, then constrain
  part of the model), calling `fix nodeTag dof` on a node that has already displaced
  by `d` drives that DOF back toward `u = 0` on the next `analyze`, dragging the node
  to its **original undeformed location** and dumping spurious strains/forces into the
  attached elements. People expect the new support to "catch" the structure at its
  current deformed shape; it does the opposite.
- **Why — the conceptual trap:** an SP constraint prescribes the **total value of the
  displacement DOF**, `u = value` (`fix` ⇒ `u = 0`), and `u` is *always* measured from
  the original mesh at t=0. There is **only one displacement frame and it never
  re-zeros** — not at a stage boundary, not ever. Since `position = X_ref + u`, the
  statements "fix the deformation to zero" and "send the node back to its original
  position" are *identical* (`u=0 ⟺ position = X_ref`). The constraint is an algebraic
  equation on absolute `u`, not an incremental/ratchet condition on the change-from-now,
  so adding it later does **not** rebase `u` to the current state. "Constraints fix
  deformation, not position" is true but misleading — it's deformation *measured from
  the undeformed configuration*.
- **Source mechanics (this build):** the current displacement at constraint-add time
  *is* captured — `SP_Constraint::setDomain()` records `initialValue = U(dofNumber)`
  (SP_Constraint.cpp:380) — and all three handlers are wired to subtract it (Penalty
  `resid = alpha*(constraint - (nodeDisp - initialValue))`, PenaltySP_FE.cpp:139;
  Lagrange LagrangeSP_FE.cpp:143; Transformation under `#define TRANSF_INCREMENTAL_SP`
  in TransformationDOF_Group.h:44 → TransformationDOF_Group.cpp:1055). BUT that
  "stay-in-place" path only fires when `retZeroInitValue == false`, and
  `getInitialValue()` returns `0` whenever it's `true` (SP_Constraint.cpp:317).
  **`fix` and `sp` both default `retZeroInitValue = true`**, and the `sp -subtractInit`
  flag *also* just sets it `true` (OpenSeesPatternCommands.cpp:1065) — so the
  incremental/stay-in-place branch is compiled in but **not cleanly reachable from
  script**. Default behavior across all handlers = snap to reference.
- **Workaround:** to install a support that holds the *current* deformed position with
  zero initial force, prescribe the current displacement explicitly rather than `fix`:
  `d = ops.nodeDisp(n, dof); ops.sp(n, dof, d, '-const')` (needs an active pattern or
  `-pattern N`). To genuinely return the DOF to its t=0 position, `fix` is correct and
  the forces are physical. In dynamics, any sudden BC change also injects an impulse;
  ramp the prescribed value via a timeSeries. Learned 2026-05-31.

### …but `equalDOF`/MP constraints are the OPPOSITE — added mid-stage they PRESERVE the offset (no snap)
- **The asymmetry (this is the surprising part):** unlike SP, an MP constraint
  (`equalDOF`, `rigidLink`, `rigidDiaphragm`) added after a node has displaced does
  **not** snap the constrained node onto the retained one. It ties their *future
  increments* together while **preserving the relative offset that existed at tie-time**,
  with **zero initial constraint force**. This is the "install at the current deformed
  state" behavior you'd *wish* `fix` had — and for MP it's the default, no flag needed.
- **Why:** `MP_Constraint::setDomain()` captures BOTH nodes' current disps at add-time
  — `Uc0` (constrained), `Ur0` (retained) (MP_Constraint.cpp:294-313) — and every
  MP-capable handler enforces the relation on the **offset-removed** displacements,
  *unconditionally* (no `retZeroInitValue` equivalent): Penalty/Lagrange build the
  residual from `Uc - Uc0` and `Ur - Ur0` (PenaltyMP_FE.cpp:230-238 → equilibrium is
  `(Uc - Uc0) = C·(Ur - Ur0)`, not `Uc = C·Ur`); the Transformation handler under
  `TRANSF_INCREMENTAL_MP` transforms only the **increment** (`modUnbalance -=
  modTrialDispOld`, TransformationDOF_Group.cpp:525) and applies it via `incrTrialDisp`
  (line 560), so the standing offset is never overwritten. At tie-time `Uc=Uc0`,
  `Ur=Ur0` ⇒ constraint satisfied with zero force and zero jump.
- **Net rule:** SP (`fix`/`sp`) defaults to enforcing the **absolute** value → snaps to
  reference; MP (`equalDOF`/rigid) defaults to enforcing the **increment** → preserves
  the current offset. Same "capture init disp at add-time" machinery underneath,
  **opposite defaults** (MP was hardened for staged construction; SP kept its legacy
  absolute-value default and never exposed the incremental toggle to a command flag).
  Caveat: holds for the MP-capable handlers (Transformation / Penalty / Lagrange); the
  Plain handler isn't the one to use for nontrivial MP staged ties. Learned 2026-05-31.
  Full write-up with source trail: [[constraints_reference_position]].

## Quirk: bundled `OTHER/LAPACK` ships `dsygvx` but NOT `dsygv` (Linux link break)

The fork's in-tree reference LAPACK (`OTHER/LAPACK/`, statically linked by the
Unix/Ubuntu Zone-A build) is a *curated subset* of LAPACK, not the full library.
It provides the **expert/range** symmetric-definite generalized eigensolver
`dsygvx.f` (line ~147 of `OTHER/LAPACK/CMakeLists.txt`) but does **not** provide
the simple driver `dsygv.f`. It also ships `dsyevx` but not `dsyev`, `dggev` but
the `x`-suffixed expert drivers are the ones actually bundled for the symmetric
path.

Consequence: fork code that calls `dsygv_` compiles per-TU fine and links fine
on **Windows** (full MKL resolves `DSYGV`), but fails the **Linux** link with
`undefined reference to 'dsygv_'`. This silently blocked the entire ladruno
Zone-A CI Build step (and thus every downstream runtime test battery) after
`CriticalTimeStep.cpp` (CentralDifferenceLadruno, PR #22) introduced a `dsygv_`
call — even though the recorder/MPCO work was unrelated.

**Rule for fork code:** for a symmetric-definite generalized eigenproblem, call
the bundled **`dsygvx_`** (expert driver), not `dsygv_`. To get just the largest
eigenvalue use `RANGE='I'` with `IL=IU=N`; the value lands in `W[0]`. Mirror the
existing, proven usage in
`SRC/system_of_eqn/eigenSOE/SymmGeneralizedEigenSolver.cpp`. Adding `dsygv.f` to
the bundle instead is the wrong fix: it cascades (dsygv calls `dsyev`, also not
bundled) and edits the upstream LAPACK source list.
*(Fixed in CriticalTimeStep.cpp via `fix/criticaltimestep-dsygvx-link`.)*

## OPS_Recorder is compiled sequentially (no `_PARALLEL_PROCESSING`) for ALL targets

`SRC/recorder/LadrunoRecorder.cpp` (and the frozen `MPCORecorder.cpp`) live in the shared
`OPS_Recorder` static lib, which CMake compiles **once with the sequential define
set** (no `-D_PARALLEL_PROCESSING`, no `-D_MUMPS`) and links into OpenSeesPy,
OpenSeesMP **and** openseesmp alike. Consequence: any `#ifdef _PARALLEL_PROCESSING`
code in a recorder is **dead** — it compiles out in every artifact, including the
MP ones. The frozen `MPCORecorder`'s only `_PARALLEL_PROCESSING` block (the
per-step stage-stamp `MPI_Allreduce`) is therefore never active in this build; its
per-partition output works purely via `sendSelf`/`recvSelf` + `send_self_count`
(always compiled, no MPI). A recorder that needs its MPI rank cannot call
`MPI_Comm_rank` — read the launcher's `PMI_RANK`/`PMI_SIZE` env instead (Intel MPI
and MS-MPI both export them per rank; verified `mpiexec -n N` sets them).

Also: the openseespy parallel model here is `_PARALLEL_INTERPRETERS` (one
independent interpreter per rank, manual `getPID`/`getNP` partition — see
`example_mpi_paralleltruss_explicit.py`), NOT `_PARALLEL_PROCESSING`
(PartitionedDomain + `sendSelf` broadcast). The `partition` command (auto domain
decomposition, the path that *would* exercise `sendSelf`) is **METIS-4-blocked**
(`OPS_partition` returns an error; needs METIS 5 / `OPS_HAVE_METIS5`), so the
broadcast path is not runtime-testable in this build.

### Ladruno recorder nested result-group names need the `<display>` parent pre-created
- **Bites:** writing an element result whose `schema.name` is nested
  (`"<display>/<bucket>"`, e.g. `stress/204-FourNodeQuad[201:0:0]`) under a parent
  that does NOT already contain the intermediate `<display>` group → the
  `H5Gcreate` in `createResultGroup` returns an invalid handle and the datasets
  underneath silently fail (HDF5-DIAG noise / missing data).
- **Why:** the recorder's group-creation property list (`h_group_proplist`) is a
  GCPL with `CRT_ORDER_TRACKED|INDEXED` but carries **no create-intermediate-group
  LINK property**, and `createResultGroup` passes `lcpl=H5P_DEFAULT`. `H5Gcreate`
  only creates the *final* path component; any intermediate must already exist.
- **Status/workaround:** the time-series `StreamingSink` path is safe because the
  recorder pre-creates `ON_ELEMENTS/<display>` at init. The `EnvelopeSink` is
  self-contained, so it pre-creates each prefix segment of `m_name` inside
  `writeEnvelope` (the per-flush `H5Ldelete` only deletes the leaf link, so the
  intermediate persists). This is why element envelopes were latently broken until
  PR #45 — only flat node/domain names had ever been written. Learned 2026-05-31.

### `decode_columns` in `ladruno_format.py` must flatten COLUMN_MAP arrays — the recorder writes them 2-D `[k×1]`
- **Bites:** `int(mult[i])` (and the other per-block scalars) in `decode_columns`
  throws `TypeError: only 0-dimensional arrays can be converted to Python scalars`
  under numpy 2.x when reading a **real recorder** `.ladruno`.
- **Why:** the C++ writer stores each per-block COLUMN_MAP array via
  `createAndWrite(vec, k, 1)` → **2-D `[k×1]`**, so `arr[i]` is a `(1,)`-array, not
  a scalar. `make_synthetic.py` writes them 1-D `[k]`, which masked it; and the
  element-**parity** gate keys results by flat column index and never calls
  `decode_columns`, so no test exercised it on real 2-D output until PR #45's
  element-envelope checker.
- **Status:** fixed (PR #45) — `decode_columns` `.reshape(-1)`s GAUSS_ID/SECTION_TAG/
  FIBER_ID/NUM_COMP/MULTIPLICITY (LEVELS is consumed via `np.atleast_1d`, left as-is).
  Learned 2026-05-31.

### h5py reads of a freshly-written `.mpco`/`.ladruno` HANG on HDF5 file locking
- **Bites:** a venv-python checker calling `h5py.File(path, "r")` on a file the
  build-python just wrote (in the same gate run) **hangs indefinitely** — no error,
  just blocks (e.g. `parity_check.py` stuck forever).
- **Why:** HDF5's default file locking; the writer's lock/superblock state isn't
  cleared promptly on a synced/Temp FS, so the reader blocks acquiring the lock.
- **Workaround:** set `HDF5_USE_FILE_LOCKING=FALSE` in the checker's environment
  (`$env:HDF5_USE_FILE_LOCKING="FALSE"`) → opens instantly (parity 80/80). Apply to
  ALL `.ladruno`/`.mpco` read steps after a recorder run. Learned 2026-05-31.

### Parallel build OOM (`cl.exe` C1060 "out of heap") on the giant template TUs under RAM pressure
- **Bites:** with low free RAM (~1–2 GB of 28), `cmake --build ... -j8`/`-j16` dies
  with `fatal error C1060: compiler is out of heap space` on the huge template TUs
  (`OPS_AllASDPlasticMaterial3Ds.cpp`, `MPCORecorder.cpp`, `LadrunoRecorder.cpp`),
  and even ordinary TUs get OS-killed (ninja `FAILED: [code=2]` with no compiler
  diagnostic = the process was killed, not a code error).
- **Workaround:** compile the monsters **serially first** —
  `ninja -j1 CMakeFiles/OPS_Material.dir/SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp.obj`
  (and the two MPCO recorder objs) — then `cmake --build build\build\Release
  --target OpenSeesPy -j2` for the rest (ninja resumes cached objs). Don't assume a
  `code=2` with no error text is a code bug; check free RAM first. Learned 2026-05-31.

### `getCommitTag()` is a GLOBAL monotonic counter — it does NOT reset on `wipe()`
- **Bites:** any per-step recorder/series that uses `Domain::getCommitTag()` as its
  step axis (the analysis monitor's `STEP`, the profiler per-step series). Across
  several `analyze()` runs in one interpreter session — even with `wipe()` between
  them — the commitTag keeps climbing (run 1 → 0..199, run 2 → 200..399, ...). It is
  NOT a within-analysis 0-based step index.
- **Implication:** don't assert absolute step values or compare step arrays across
  runs by value; compare the *stride* (`np.diff(step) == every`) instead. For a live
  viewer, treat `STEP` as a monotonic id, not "step N of this analysis."
- Learned 2026-05-31 building the analysis monitor (`08_analysis_monitor.md`).

### `recorder ladruno` is NOT wired into the classic Tcl `OpenSees.exe` (was; now fixed)
- **Bites:** `recorder ladruno ...` works from OpenSeesPy/openseesmp and the
  interpreter-based Tcl (`TclWrapper`→`OPS_Recorder`, the shared map in
  `OpenSeesOutputCommands.cpp`), but the **classic** Tcl `OpenSees.exe`
  (`commands.cpp`→`addRecorder`→`TclAddRecorder` in `TclRecorderCommands.cpp`)
  hardcodes its recorder dispatch in a *separate* file that the rename PR never
  touched — it had `mpco`/`vtkhdf`/`gmsh`/`EnergyBalance` but no `ladruno`. So
  `recorder ladruno` raised "recorder type ladruno is unknown" only in `OpenSees.exe`.
- **Fix:** added the `else if (strcmp(argv[1],"ladruno")==0)` branch (+ extern
  `OPS_LadrunoRecorder`) mirroring the `mpco` block in `TclRecorderCommands.cpp`.
  Lesson: there are TWO recorder-dispatch tables (shared `OPS_Recorder` map for
  Py/interpreter-Tcl; hardcoded `TclAddRecorder` for classic Tcl) — wire new
  recorders into BOTH. Learned 2026-05-31.

### Ladruno recorder `modesOfVibration` (eigen) output writes no data — modal `DATA` group/dataset collision
- **Bites:** `recorder ladruno -N modesOfVibration` after `ops.eigen(n)` creates the
  `MODEL_STAGE[*]/RESULTS/ON_NODES/MODES_OF_VIBRATION(U)` group with a valid schema
  (ID, COMPONENTS) but `DATA` stays empty `(0, nNodes, nComp)` and no `MODE_k`
  datasets appear; HDF5-DIAG "can't synchronously write data / Write failed" errors
  fire. Happens for ANY model (reproduced with a bare elasticBeamColumn portal —
  NOT the known fiber-section `writeSections` noise), under both `ops.record()` and
  an `analyze()` step.
- **Why:** `LadrunoRecorder::recordModeChannel` (LadrunoRecorder.cpp ~1692) calls
  `ch.sink->begin()`, which (StreamingSink) creates `.../MODES_OF_VIBRATION(U)/DATA`
  as a **chunked dataset** (normal time-series layout), then tries to
  `h5::group::create` a **group** at `DATA/STEP_<step>` with `MODE_k` datasets
  under it (the MPCO modal layout). A group cannot be created beneath an existing
  dataset → the HDF5 calls fail, no modal data is written.
- **Status:** **FIXED 2026-05-31.** `recordModeChannel` no longer calls the StreamingSink
  `begin()`. It now owns the modal init, mirroring frozen `ResultRecorderModesOfVibration::record`:
  once per stage (idempotent via `H5Lexists`) it creates the result group
  (`h5::group::createResultGroup`) + `ID` dataset + `DATA` **group**, then per step writes
  `DATA/STEP_<step>/MODE_<k>` datasets with MODE/LAMBDA/OMEGA/FREQUENCY/PERIOD attrs. The
  validator `ladruno_format.py::_check_data_shape` was taught the modal layout (DATA = group
  of STEP_<step> groups of MODE_<k> datasets). **Modal eigenvectors now match frozen mpco to
  1e-12**; the EIGEN gate is promoted into the counted regression battery and `eigen_check.py`
  does a real modal value-parity diff vs the ref `.mpco`. Files: `SRC/recorder/LadrunoRecorder.cpp`,
  `ladruno_format.py`, `eigen_check.py`, `run_regression.bat`. Found + fixed 2026-05-31 by the
  new eigen coverage gate (the test scheme catching, then confirming the fix of, a real bug).

### Energy-balance check script `energy_check.py` was stale vs the chunked `DATA` layout (fixed)
- **Bites:** `energy_check.py` failed with `TypeError: Only 1D arrays allowed for
  fancy indexing` — its `_read_result` iterated `grp["DATA"]` as a per-step *group*,
  but the recorder writes `ON_DOMAIN/ON_REGIONS energyBalance DATA` as a chunked
  `[T×nrows×ncomp]` **dataset** (the standardized streaming layout; recorder output
  is correct).
- **Fix:** read via `lf.iter_step_slices(grp)` (the canonical slicer the parity
  checks already use; tolerates both chunked and legacy `DATA/STEP_k` layouts).
  After the fix the energy kernel matches the EnergyBalance text sidecar to ~5e-9.
  Learned 2026-05-31.

### `CrdTransf::getLocalAxes` base default zeros the axes — wiring `localAxes` on an element whose transform doesn't override it emits a *degenerate* frame
- **Bites:** when extending the Ladruno `"localAxes"` (response id 30) coverage to
  more beam elements, `DispBeamColumn2dInt` *looks* trivially wireable — it owns a
  `crdTransf` and the copy-paste pattern compiles. But its transform is
  `LinearCrdTransf2dInt`, which does **not** implement `getLocalAxes`, so it falls
  through to the base `CrdTransf::getLocalAxes` (`SRC/coordTransformation/CrdTransf.cpp:126`)
  that simply `Zero()`s xAxis/yAxis/zAxis and returns 0.
- **Why it matters:** the recorder's `writeModelLocalAxes` records *any* element that
  answers `"localAxes"`. A non-null response carrying an all-zeros frame would be
  written/quaternion-converted as a **degenerate** orientation — strictly worse than
  the current behaviour, where a silent (no-response) element falls back to a clean
  identity quaternion.
- **Status/workaround:** `DispBeamColumn2dInt` is **deliberately left unwired**. To
  wire it correctly, first give `LinearCrdTransf2dInt` a real `getLocalAxes` (mirror
  `LinearCrdTransf2d::getLocalAxes`), then add the id-30 response. The standard-transform
  beams (Elastic/Force/Disp/Mixed/GradientInelastic Beam(Column)2d/3d) are all safe —
  their LinearCrdTransf2d/3d, PDelta*, and Corot* transforms all override `getLocalAxes`.
  Learned 2026-06-03.

### LadrunoRecorder shell per-layer stress: the verb is a DOTTED token, and `material.fiber.*` needed a fix
- **Bites (1) — syntax:** the `-E` element verb is a single **dot-joined** token that the
  recorder splits on `'.'` (`recorder('ladruno', f, '-E', 'section.fiber.stress')`), NOT
  space-separated args. Passing `'-E','section','fiber','stress'` stores only `["section"]`
  (the rest are parsed as later options), which then **segfaults** (see bite 3).
- **Bites (2) — the real bug:** for a layered shell, the *obvious* per-layer verb
  `material.fiber.stress` silently emitted **no element bucket**. Root cause: the request
  builder only set `do_all_fibers` for `fiber` under `do_all_sections`, and the
  `do_all_materials` path had no fiber-index expansion — so the fiber id was never
  substituted and `setResponse` returned null for every element. `section.fiber.stress`
  worked because the recorder swaps section->material for shells and runs the section
  fiber-expansion. The section-level read itself was always fine
  (`LayeredShellFiberSection`/`MembranePlateFiberSection` answer `fiber <i> <resp>` and tag
  `FiberOutput`; `eleResponse(tag,'material','1','fiber','1','stress')` works directly).
- **Fix:** extend the `fiber` trigger to `do_all_sections || do_all_materials` and give the
  `do_all_materials` branch the same per-(gp,layer) expansion as `do_all_sections` (driven
  by the shared `elem_ngauss_nfiber_info` discovery table). Now `material.fiber.stress`
  emits the per-layer bucket **byte-identical** to `section.fiber.stress` (verified: 4 GP ×
  3 layers × 5 comp = 60 cols, maxdiff 0.0). Regression gate `SHELL LAYER STRESS`
  (`shell_layer_model.py`/`shell_layer_check.py`) in `run_regression.bat`.
- **Bites (3) — latent crash (FIXED):** a bare `-E section` / `-E material` (the keyword
  with NO sub-verb) segfaulted — `request_mod` was `["section",""]` (argc=2), so the element
  was queried as `["section"/"material", <id>]` and forwarded a **zero-length** arg list
  (`&argv[2]`, argc-2==0) to its section's `setResponse`, which derefs `argv[0]`. Affected
  any shell or beam. **Fix:** guard both non-fiber `setResponse` call sites in
  `initElementSources` — only call when `argc > (int)<keyword>_id_placeholder_index + 1`
  (i.e. at least one sub-verb token follows the id). A bare verb now emits no bucket instead
  of crashing. Regression gate `SHELL BARE VERB` (`shell_bare_verb_model.py` /
  `shell_bare_verb_check.py`). Learned + fixed 2026-06-03.
- **Known cosmetic limitation:** the shell per-layer COLUMN_MAP records `fiber_id=-1`,
  `section_tag=-1`, and `UnknownStress` component names (layer identity is flattened into a
  running `gauss_id`). This is IDENTICAL on the already-working `section.fiber.stress` path
  (not introduced here) — a metadata refinement for later, the stress *values* are correct.
  (Adversarial-review note: a reviewer claimed this *collapses* layers into one
  `normalize_element` key and loses data — REFUTED: the 3 layers map to distinct `gauss_id`
  0..11, so the parity dict has 4·12·5=240 distinct entries, no loss.)

### LadrunoRecorder `sendSelf`/`recvSelf` config must stay in LOCKSTEP — two fields were missing (FIXED) + the rule
- **Bites:** any recorder config field set by the OPS_ parser must be (a) `ser.put_*`'d in
  `sendSelf` AND (b) `de.get_*`'d in `recvSelf` **in the same order**, or OpenSeesMP worker
  ranks (built via `recvSelf`) silently diverge from P0 (built via the parser). Two fields
  were configured but NOT transmitted: `m_data->envelope_mode` (`-envelope`) and
  `m_data->info.store_data_f32` (`-precision f32`). Consequence: with `-envelope` or
  `-precision f32` under MP, P0 wrote ENVELOPES/f32 while every worker wrote full
  time-series/f64 → `.part-N` files with a *different schema* than `.part-0`, breaking the
  apeGmsh stitch-on-read. **Fix:** append both to `sendSelf` (after the elemental-results
  block) and `recvSelf` (after the same block), same order. **Rule for the next agent: when
  you add ANY field to the recorder config, grep `sendSelf`/`recvSelf` and add it to both.**
  (Found by the 2026-06-03 adversarial review. The MP round-trip itself was verified by
  construction — symmetric put/get — not by a live 2-rank run, since the worktree had only
  the OpenSeesPy build; the `mp_parallel` gate needs OpenSeesPyMP.)

### LadrunoRecorder smaller robustness fixes from the adversarial review (2026-06-03)
- **Mixed-dimension node OOB (FIXED):** `writeModelNodes` latched `ndim` from the *first*
  node then read `crds[1]`/`crds[2]` unchecked for every node — a node carrying fewer coords
  than `ndim` read past its `Vector`. Upstream MPCO guards this; the port had dropped it.
  Fix: clamp each read to `crds.Size()` (pad 0.0), mirroring the GLOBAL_GP path.
- **`eo_response` leak (FIXED):** in `initElementSources`, a `CompositeResponse` built but
  then rejected because `eo_stream.error_code != OK` (e.g. `ERROR_CODE_GENERIC`) was owned by
  nobody. Fix: `else if (eo_response) delete eo_response;` after the registration block.
- **Recorder `exit(-1)` aborting the analysis (FIXED):** the `OutputDescriptorStream`
  tag/attr parser and `mapElements` in `Ladruno_ElementResults.h` called `exit(-1)` on an
  element output-tag nesting they didn't expect (invalid parent for SectionOutput/FiberOutput,
  invalid tag at level, empty item-list at a walked level) or on two same-classTag elements
  with differing `getNumExternalNodes()` — a *recorder* killing the whole run. **Fix:** the 7
  stream sites now `error_code = ERROR_CODE_GENERIC; return -1;` (the offending element's bucket
  is dropped by the existing `error_code != OK` gate — same mechanism `ensureItemsOfUniformType`
  already used at line ~1036); `mapElements` now `continue;`s past the inconsistent element.
  No live `exit(-1)` remains (two pre-existing commented-out ones at ~770/~1033 untouched).
  Happy path unchanged (full recorder regression green). A runtime trigger needs a custom
  element that emits an unsupported tag nesting (not reachable from standard openseespy), so
  this is verified by the no-regression run + reuse of the already-proven GENERIC drop-path.
- **Still open:** `domainChanged()`/`restart()`/`setDomain()` are no-ops, so the only
  source-rebuild trigger is the `hasDomainChanged()` stamp inside `record()` (cached
  `Element*`/`Response*` can dangle if a model edit doesn't move the stamp across a record).
### uri `viscous` hourglass: explicit run goes silently unstable at large eps + dt (CDL doesn't trap it)
- **Bites:** `LadrunoBrick -hourglass viscous` adds a velocity-proportional damping
  force but NO hourglass stiffness. Under CentralDifferenceLadruno the viscous term
  has its own explicit stability bound; at `eps≈0.5` with `dt = 0.1·le/c` the
  hourglass mode blew up to `nodeDisp ~ 1e+99` — yet `analyze()` still returned 0
  (CDL does not check for NaN/Inf), so the run *looks* like it completed. Smaller
  eps tolerates the larger dt; large eps needs a smaller dt.
- **Fix / rule:** for viscous-hourglass explicit runs use a conservative step
  (`dt ≈ 0.02–0.03·le/c`) and modest `eps` (≤0.1). Don't trust a clean `analyze()`
  return alone — assert `isfinite(nodeDisp)`. (Element-level: the viscous tangent
  is rank-deficient ⇒ statics is singular; it is explicit-only by construction.)
  Learned 2026-06-01.

### Viscous dissipation reported via the discrete work integral `∫f·du` is a DIAGNOSTIC, not an exact energy balance
- **Bites:** `LadrunoBrick::hourglassEnergy()` for uri-viscous returns a committed
  accumulator `hgDissipated += c_visc·Σ q̇·Δq` (work against the FB rate damper).
  For LIGHT damping this tracks the true dissipated energy well (≈82% of the
  hourglass KE recovered over a long run); for STRONG damping the per-step velocity
  collapse in the leapfrog stagger makes `f·Δu` UNDER-count, so "more damping ⇒
  more reported dissipation" is FALSE as measured (it is non-monotone in eps).
- **Rule:** treat it as a monotone, energy-bounded spurious-mode diagnostic
  (GLSTAT-style), and write tests around the robust properties — non-decreasing,
  positive under hourglass excitation, `0 < E ≤ KE_imparted` across eps, exactly 0
  for a rigid/constant-strain velocity (γ⟂linear) — NOT exact energy convergence.
  Learned 2026-06-01.

### F-bar element tangent is GENERALLY UNSYMMETRIC (needs an unsymmetric solver)
- **Bites:** `LadrunoBrick -geom finite -formulation bbar` (F-bar, dSNPO §15.1). The
  consistent tangent (eq 15.10) is `K = ∫Gᵀa G dv + ∫Gᵀq(G₀−G) dv`; the second
  (F-bar coupling) term is **not symmetric** in general — the book says so
  explicitly (after eq 15.10): "the additional stiffness term … is generally
  unsymmetric and, therefore, requires an unsymmetric solver." So `system BandSPD`
  / `ProfileSPD` (symmetric storage) will silently use only the upper triangle and
  **converge to the wrong answer or stall**. Use `system FullGeneral` (or any
  unsymmetric solver). The symmetric storage path will *not* error — it just drops
  the lower triangle, so this is a silent-wrong-answer trap.
- **Why:** F̄ = (J₀/J)^(1/3)F couples every Gauss point's stress to the element
  centroid dilatation; the linearization mixes the GP gradient G with the centroid
  gradient G₀, and `q⊗(G₀−G)` is a non-symmetric outer product. (The plain
  `-formulation std` finite tangent stays symmetric — the coupling term is absent.)
- **Also:** the eq 15.11 coupling tensor is `q = (1/3)a:(I⊗I) − (2/3)σ⊗I` with `a`
  the **full** spatial tangent (= `c − σδ`, the same modulus as the standard term),
  **not** the material modulus `c` alone. A first-principles spatial shortcut
  (`dσ̄ = c̄:sym(L̄)`) drops the `−(2/3)σ⊗I` and gives a subtly wrong tangent that
  still "looks right" under a crude FD check — the element FD-tangent test against an
  *analytic* material tangent is what catches it. Learned 2026-06-02 (F-bar impl).
### Wrapping a KINEMATIC-hardening material in `LogStrain` (finite strain) is NOT objective under large rotation — backstress doesn't co-rotate
- **Bites:** `nDMaterial LogStrain $t $j2` over a combined-hardening `LadrunoJ2`
  (or any backstress material) gives finite-strain J2 that is **exact for the
  isotropic part but loses frame-indifference for the kinematic part under a large
  rotation**: superpose a finite rotation `Q` on a plastically-loaded state and the
  Cauchy stress does NOT come back as `Q σ Qᵀ` (principal stresses change).
- **Why:** the log-strain (MATISU) wrapper co-rotates only the elastic state it owns,
  `Bᵉᵗʳ = F_Δ Bᵉ_n F_Δᵀ`. Isotropic yield sees only `‖s‖,ε̄ᵖ` (rotation-invariant), so
  it is objective. The **backstress `α` lives inside the UNCHANGED small-strain inner
  in a FIXED frame** — the wrapper never rotates it — so `‖M‖=‖s−α‖` is not
  rotation-invariant once `α≠0`. It is a framework limit, not a bug: the direct
  Box-14.4 chain shows the identical behaviour. This is exactly the
  kinematic-hardening-at-finite-strain case dSNPO defers to **§14.11**.
- **Rule:** the simple `LogStrain`-wrap of a backstress material is correct only for
  **no / small rotation** (or pure isotropic hardening). For exact large-rotation
  combined hardening you need a finite-strain-NATIVE material that co-rotates `α`
  every step (push `α` forward by the incremental rotation) — which is why the J2
  return map was extracted into the OpenSees-free `LadrunoJ2Kernel.h`: a future
  `FiniteStrainNDMaterial`-subclass J2 can reuse the verified map and add the `α`
  push-forward. Pin the boundary with a **strict xfail** so v2 flips it green
  (`tests/test_ladrunoJ2_finite.py`). Learned 2026-06-02.

### LadrunoJ2 consistent-tangent denominator `h` assumes NON-softening hardening
- **Bites:** the analytic consistent tangent in `LadrunoJ2Kernel.h` divides by
  `h = dtheta + (2/3)·sig_y'(pbar) − n:Mp` (commented `= −df/ddG > 0`) to form
  `betaNN` and `betaMpN`. For standard hardening (`Hiso,Qinf,Cₖ,γₖ ≥ 0`) `h>0`
  always. But the material accepts arbitrary user params: a **negative `Hiso`
  (linear softening) or `Qinf<0`** makes `sig_y'` negative and can drive `h→0` or
  `h<0` ⇒ an `inf`/`NaN` or a sign-flipped (non-physical) tangent that poisons the
  global Newton with no diagnostic. The local scalar-Newton residual `dR` has the
  same exposure.
- **Why it's not "fixed":** this is **pre-existing** behaviour inherited verbatim
  from the original `integrate()` (the 2026-06-02 kernel extraction was deliberately
  bit-identical, so it was preserved, not introduced). The model is designed for
  hardening/perfect plasticity; softening is out of its intended envelope.
- **Rule:** do not feed `LadrunoJ2` softening hardening params. If softening is ever
  wanted, the right fix is **parameter validation at construction** (reject/warn on
  `Hiso<0`/`Qinf<0` that can violate `h>0`) plus a kernel guard
  (`if (fabs(h) < eps·stressScale)` → fall back to the elastic-predictor tangent,
  mirroring the existing `‖M‖→0` `normFloor` treatment) — a separate PR, not a
  bit-identical-extraction change. Surfaced by the PR #97 adversarial review,
  2026-06-02.

### `ZeroLengthSection` requires `-ndf 3` (2D) / `-ndf 6` (3D) — silently absent otherwise
- **Bites:** building a `zeroLengthSection` in a reduced-DOF model (e.g. an axial
  SDOF on `-ndf 2`) prints *"ZeroLengthSection::setDomain() -- element only works
  for 3 (2d) or 6 (3d) dof per node"* ([ZeroLengthSection.cpp:247]) and then the
  element is **not added** — but `analyze()` still runs, on a model with no spring,
  so the response looks like an undamped/zero-stiffness free body (constant disp,
  no oscillation). Plain `ZeroLength`/`TwoNodeLink` have no such restriction.
- **Why:** ZeroLengthSection maps the full section response set (P, Vy, Mz, …) onto
  the element DOFs and assumes the complete 3-/6-dof node layout.
- **Rule:** use `-ndf 3`/`-ndf 6` and fix the unused DOFs; never `-ndf 2`. Caught
  by `tests/test_spring_damping_claims.py`. Learned 2026-06-02.

### Prescribing ALL of a node's DOFs via `sp` with the Transformation handler ⇒ 0 free equations ⇒ process terminates
- **Bites:** a static test that imposes both DOFs of the only free node via `ops.sp`
  (with the other node fully `fix`ed) leaves the system with **zero unknowns**.
  Under `constraints('Transformation')` the solve does not return an error code —
  it **terminates the process** mid-`analyze()` (no Python traceback, exit 0),
  which under pytest aborts the whole run with no summary (looks like a hang).
- **Why:** the Transformation handler condenses out the constrained DOFs; with none
  left the assembled system is degenerate and the path hits a hard exit rather than
  a graceful failure (same family as the MPCO `exit(-1)` kernel-kill pattern).
- **Workaround:** use `constraints('Penalty', 1e14, 1e14)` for fully-prescribed
  configurations (DOFs are retained and penalised, so ≥1 equation remains), and
  read element force via `eleForce` rather than `nodeReaction`. Learned 2026-06-02.

## `setStrain` (testUniaxialMaterial) COMMITS — central FD tangent is invalid for plasticity

The Python `setStrain(eps)` command (`SRC/interpreter/OpenSeesCommandsPython.cpp`,
`ops_setStrain`) calls **both** `setTrialStrain(eps)` **and** `commitState()`. So
the `_testbed.fem_checks.uniaxial_tangent_fd` central difference — which probes
`strain0-d`, `strain0`, `strain0+d` expecting all three from one committed state —
straddles an **elastic unload** on the minus side for any path-dependent (plastic)
material, returning ~`(E + E_alg)/2` instead of the consistent tangent (caught on
LadrunoUniaxialJ2: analytic 176.7 vs the broken FD 597 ≈ (1000+176)/2). The helper
is only valid for nonlinear-elastic materials.

- **Rule:** to FD a *plasticity* consistent tangent, probe each point as an
  INDEPENDENT one-step return from a FRESH material (`wipe` → redefine →
  `testUniaxialMaterial` → `setStrain`), so all probes are one-step-from-zero; the
  central difference of that one-step map IS the algorithmic tangent (O(d²)). See
  `test_consistent_tangent_fd` (V6) in `tests/test_ladrunoUniaxialJ2_material.py`.
  Also note `testUniaxialMaterial(tag)` returns the SAME object (not a copy), so it
  does not reset committed state — you must rebuild the material to get a clean one.
  Learned 2026-06-02 (LadrunoUniaxialJ2 polish).

## Lemaitre damage tangent: the `dotT` weight cancels the tensor↔engineering shear factor

When assembling `∂D/∂ε` for the Lemaitre coupled-damage consistent tangent
(`LadrunoJ2Kernel.h::returnMapDamaged`), the `∂(p̄)/∂ε` term goes through
`∂‖M‖/∂M = wt·M/‖M‖` where `wt = (1,1,1,2,2,2)` are the `dotT` factor-2 weights on
the off-diagonal (shear) pairs. The strain variable being differentiated is then
converted tensor→engineering by `w = (1,1,1,½,½,½)`. **`wt_j · w_j = 1` for every
component** (normal: 1·1; shear: 2·½), so the net engineering derivative is just
`(2G/h)·n_j` with **NO shear half-factor**. The natural-but-wrong instinct is to
carry the `w_j=½` on shear (mirroring the strain-mapping elsewhere in the file);
that silently halves the shear columns of `∂p̄/∂ε` and breaks the tangent ONLY on
shear-coupled (3D-mixed) states — uniaxial paths pass clean, so a uniaxial-only FD
check misses it. Caught by the 3D-mixed FD case in `tests/ladruno_damage_check.cpp`
(error was a fixed 1.5e-5 independent of the FD step ⇒ a real bug, not truncation).

- **Rule:** for a `dotT`-norm gradient differentiated w.r.t. engineering strain, the
  `dotT` weight and the tensor→engineering factor cancel — use the bare component.
  Always FD-check the consistent tangent on a **shear-inclusive** state, not just
  uniaxial. Learned 2026-06-02 (Lemaitre damage, [[15_lemaitre_ductile_damage_adr]]).

## LadrunoIMKBeam: "elastic end" ≠ "released end" (pin); fake a release with a tiny Elastic hinge

A `LadrunoIMKBeam` end with **no hinge material** in that slot is **elastic** — it
still carries moment with the full elastic rotational stiffness (`4EI/L`,
far-end-fixed). That is NOT a structural **release** (pin), where the end moment is
identically zero and the rotation is free (the far end condenses out, dropping the
active end stiffness to `3EI/L`). The element has no `-release` flag (deferred, see
[[14_ladruno_imk_beam]] §8).

- **Workaround:** put a near-zero-stiffness `Elastic` uniaxial in that end/axis
  slot (e.g. `-matZj`). The series flexibility `1/k → ∞` zeroes the end moment.
  Use `k ≈ 1e-5·(4EI/L)` of the bending axis (~1e-5 of the elastic rotational
  stiffness): verified to give `k_i = 3EI/L` to ~5 figures and residual
  `M_j/M_i ≈ 7e-6`. Stay above the element `ktFloor` guard (`1e-8·4EI/L`); a
  factor in `[1e-6, 1e-4]` is the sweet spot (smaller hurts conditioning, larger
  leaves a non-negligible residual moment). Learned 2026-06-02 (LadrunoIMKBeam).

## LadrunoBrick `-formulation eas` is a misnomer (it's SSPbrick) → rename to `ssp` RECOMMENDED

The `-formulation eas` single-point element is **never** a Simo–Rifai enhanced-assumed-strain
element. It is a **verbatim port of `UWelements/SSPbrick::GetStab`** (McGann/Arduino) — B̄ +
a statically-condensed assumed-strain hourglass `Kstab` baked once at `setDomain` on the
*initial* tangent (no per-step α). Proven by source (`LadrunoBrick.cpp:1887`) and by
byte-identity to `SSPbrick` (6 figs elastic, 4 figs plastic, 0.2% with damage — see the
Lemaitre validation §4.6 + `tests/test_ladrunoBrick_element.py::test_eas_matches_sspbrick_distorted_hex`).

- **RECOMMENDED rename (NOT yet in source):** call it `-formulation ssp`; keep `eas` as a
  deprecated alias (warn → ssp) so the name `eas` is freed for a **true** Simo–Rifai EAS element.
  Suggested impl: enum `Formulation::SSP` with `EAS = SSP` back-compat alias; parser branch in
  `OPS_LadrunoBrick.cpp`; `formulationName→"ssp"`. The owner will land the C++ rename separately
  (docs-only PR keeps the current `eas` name). Tests can use a runtime fallback (`ssp` if built,
  else `eas`) to stay green across the rename.
- **Consequence (why it's not a bug):** being a single-point element, `eas`/`ssp` over-predicts load
  in a plastic/damage stress gradient (Jensen on the concave σ(ε); centroid under-samples) — but it
  shares this *exactly* with the validated `SSPbrick`, and it converges to `bbar` under refinement
  (gap 17.4%→3.7% to h=0.5). Use `bbar` for gradient/fracture fields, `eas`/`ssp` for smooth + cost.

## Shared append-point files conflict on stale branches — append at END + `merge=union`

Every new fork feature touches the same hotspot files (`SRC/classTags.h`, the
broker `FEM_ObjectBrokerAllClasses.cpp`, the `OpenSees*Commands.cpp` registries,
`*/CMakeLists.txt`, the banner `banner_features.txt`/`tclMain.cpp`/`PythonModule.cpp`,
and the `LEDGER_*.md` / testbed `manifest.yaml` bookkeeping). When a feature branch
falls behind `ladruno` (e.g. #127 was 48 commits stale), these are exactly where
merge conflicts land — git auto-merges *additions at different lines*, but two edits
to **adjacent** lines (a new row inserted next to a row another PR also edited) do
NOT auto-merge.

Prevention (all three, not just one):
- **Reconcile with latest `ladruno` right before merging** (rebase or merge-from-base
  — equivalent under squash-merge; rebase = linear + force-push, merge = no force-push).
  This fixes staleness, but NOT contention.
- **Append new entries at the END of a list/section, never interleaved.** A `classTags.h`
  tag appended after the last one auto-merges; a `LEDGER_*.md` row inserted *mid-table*
  next to a row another PR edits will conflict (the #127 case — its LadrunoJ2Finite row
  sat above the UniaxialJ2 row that the Lemaitre PR was editing).
- **`merge=union` driver** (set in `.gitattributes`) for the append-only logs
  (`LEDGER_*.md`, `banner_features.txt`) — git keeps BOTH sides instead of conflicting.
  NEVER apply `union` to source code (it interleaves → garbage). Learned 2026-06-02
  (the #127 rebase past the Lemaitre-damage merge; [[16_finite_native_j2_adr]]).
## Corot solid wrapper: external dead loads must NOT pass through `globalizeForce`; and corot IS objective for kinematic hardening

Two corot-seam gotchas, from the finite-strain trifecta deep review (2026-06-02,
[[10_solid_corotational_adr]]):

- **External dead loads must stay in the global frame.** `SolidTransformationCorot::
  globalizeForce` pushes the core force forward by R (`f_global = R f_d − …`). That is
  correct ONLY for the *internal* force (`∫ Bᵀσ`, self-equilibrated). A fixed-direction
  body/self-weight load (`-b`, `eleLoad -selfWeight`) is a GLOBAL-frame quantity — if it
  is folded into the core force before `globalizeForce`, corot rotates gravity WITH the
  element (wrong, non-conservative; was the COROT-1 bug). **Fix/pattern:** `LadrunoBrick`
  accumulates the body load in a separate `bodyForce` vector and adds it back AFTER
  `globalizeForce`/`globalizeStiff` (also keeps the spurious body-load term out of the
  corot geometric stiffness). Behavior-neutral under `-geom linear` (identity globalize);
  the `-geom finite` path was already correct (assembles in the spatial frame, no
  globalize). Any new fold-then-globalize site must keep external loads out of the core.

- **Corot is objective for KINEMATIC-hardening materials (unlike the LogStrain finite
  path).** A natural worry is that corot shares the dSNPO §14.11 backstress-frame
  non-objectivity. It does NOT: corot feeds the material `u_d = Rᵀ x_rel − X_rel`
  (REFERENCE frame) with reference-config gradients, so the small-strain material — and
  its backstress α — live in a FIXED reference frame across commits (the element's R
  rotates; the material's frame does not). Since `polar(Q·H) = Q·polar(H)` exactly for
  rigid Q, `u_d` is rigid-rotation-invariant ⇒ identical deformational-strain history ⇒
  exact objectivity (verified, `test_corot_kinematic_hardening_objectivity`). The
  LogStrain path differs precisely because `bᵉ_tr = f_Δ bᵉ_n f_Δᵀ` co-rotates the stress
  `s` into the current frame while α stays fixed — THAT is §14.11. Lesson: "the element's
  R rotates between commits" does NOT imply "the material's frame rotates."

## Finite-strain elastoplastic bending/necking BVPs need KrylovNewton (plain Newton diverges); F-bar needs an unsymmetric solver

From the finite-strain validation Phase P1 (2026-06-02, [[18_finite_strain_validation_report]]).

- **`LogStrain(LadrunoJ2)` + `LadrunoBrick -geom finite` bending *into plasticity* does
  NOT converge under `algorithm Newton` — nor under `NewtonLineSearch`.** The residual
  grows from the very first increment (a `NormDispIncr`/`EnergyIncr` norm that climbs, not
  shrinks). It is not a tangent bug (the consistent tangent passes the FD gate on
  homogeneous states): bending+plasticity on a low-order hex is just a stiff, badly-scaled
  Newton basin for a full step. **`algorithm KrylovNewton` (+ `test EnergyIncr 1e-6`) is
  robust** and converges quadratically-ish; the necking bar (C1) and the isochoric-J2
  locking cantilever (B3) both rely on it. Homogeneous single-element states and elastic
  bending converge fine under plain Newton — the divergence is specific to *inhomogeneous
  plastic* finite-strain BVPs.
- **A 1-element-wide cross-section bends too poorly for stable plastic Newton.** A 1×1×nz
  column under transverse displacement control diverges even with KrylovNewton; a ≥2×2
  cross-section is needed. (Elastic load-control on the 1-wide column is fine — it just
  locks.)
- **F-bar (`-formulation bbar -geom finite`) has an UNSYMMETRIC tangent** (dSNPO eq 15.10)
  ⇒ use `system FullGeneral` (dense) or, much faster for meshed studies, `system UmfPack`
  (unsymmetric sparse). A symmetric solver silently mis-solves. `UmfPack` made the 128–576
  hex necking runs tractable where `FullGeneral` would be dense-O(N³) per step.
- **Plastic finite-strain stress paths are path-dependent** (obvious, but it bites tests):
  a sub-stepped element solve does NOT equal a single-step return-map oracle for
  *non-proportional* loading (simple shear, equibiaxial). Drive ONE increment ref→F when
  comparing to a one-step oracle, or step the oracle incrementally over the same F_k.

## Explicit `-geom finite`: `criticalTimeStep()` is reference-config (must margin dt), and the EnergyBalance recorder reports IE with a flipped sign

From the finite-strain validation Phase P4 (Taylor-bar impact, 2026-06-02,
[[18_finite_strain_validation_report]] §7; `tests/test_finite_strain_P4_explicit.py`).

- **`ops.criticalTimeStep()` does NOT shrink as elements compress.** On the Taylor
  bar the cylinder shortened ~33 % and the impact face mushroomed >2×, yet
  `criticalTimeStep()` was *bit-identical* before and after (ratio 1.000). It is
  computed from the **reference** configuration characteristic length (review
  GEOM-2). So an explicit `-geom finite` run is only conservatively safe until
  strong compression; past that the *true* stable dt is smaller than reported.
  **Carry a safety factor < 1** — the Taylor bar uses `dt = 0.3·dt_cr` (0.5 is
  stable for the early/short transit but risks instability through full
  mushrooming). A future improvement would update dt_cr from the current config.
- **`EnergyBalance` recorder reports IE (internal energy) with a flipped SIGN for
  the finite-strain element.** On the Taylor bar `KE0=2.34e5`, `KE_final=1.0e4`
  (4.3 %, the rest absorbed plastically), and `IE_final=−2.36e5` — the MAGNITUDE
  equals the absorbed kinetic energy (≈ KE0−KE_final, within ~5 %) but the sign is
  negative, so the recorder's `RES`/`ERR%` columns read ~100 % (spurious). The KE
  column is correct (it's the validated getMass aliasing-fix path,
  `test_energyBalanceRecorder.py`); only IE's sign is off for the
  `LogStrain`/`LadrunoBrick -geom finite` path. **Work around it by comparing
  `|IE|` to the kinetic-energy change**; do not trust `ERR%` for finite-strain
  elements until the IE-increment sign convention is reconciled (likely the
  recorder integrates fᵀΔu with the internal-force sign opposite to what the
  finite element returns). Candidate follow-up: audit
  `EnergyBalanceRecorder.cpp` internal-energy accumulation vs `LadrunoBrick`
  `getResistingForce` sign under `-geom finite`.

- **`ASDConcrete3D` confines emergently, but there is NO dilation-angle input.**
  Measured (RC-3D Gate 2, `Ladruno_implementation/rc3d_gates/gate2_concrete_confinement.py`):
  a single brick under constant lateral pressure `p` + axial displacement control
  develops a confined peak `fcc` within **~5 % of Mander** for `p/fc ∈ [0, 0.20]`
  (unconfined recovers `fc` exactly), and the peak strain grows with `p` — so
  confinement is a REAL emergent property of the Lubliner triaxial surface; do
  **not** pre-inflate `fc` à la Mander in a 3D solid (that double-counts). BUT the
  *amount* of confinement is governed by the **`Kc` triaxial-meridian parameter +
  the compression hardening backbone**, NOT a dilation angle — `ASDConcrete3D`
  exposes no dilatancy/flow-rule input (grep the header for `dilatan` → nothing).
  So: validate `fcc(p)` against test data / Mander before trusting confined-member
  results; the lever to tune is `Kc` + the `-Ce/-Cs/-Cd` curve. **Backbone calibration
  gotcha:** the first compression point must be the *elastic limit* (`σ = E·ε`, so
  `Cd = 0` there); putting the first point past the elastic line makes the model run
  elastic up to that strain and the unconfined peak overshoots `fc` (≈2× in an early
  Gate-2 draft). **Solver:** confined softening needs `KrylovNewton` (or the blessed
  `Ladruno_scripts/ladruno_solve.py` adaptive driver) — plain Newton fixed-step
  diverges past the peak.

- **openseespy parsers must peek a maybe-numeric arg with `OPS_GetStringFromAll`,
  never `OPS_GetString`.** openseespy passes TYPED args; `OPS_GetString()` returns
  the sentinel `"Invalid String Input!"` when the current arg is an int or float,
  so any parser that peeks a position which could be a number (a positional count,
  or a flag value that might be `auto`/numeric like `-kt`) blows up — while string
  args at the same slot pass, making the failure look maddeningly selective. Use
  `char buf[N]; OPS_GetStringFromAll(buf, N);` — it stringifies any arg (`%d` for
  int, `%.20f` for double → exact `atof` round-trip) AND advances the cursor, then
  `atoi`/`atof`/`strcmp`. Tcl is all-strings so it never reproduces there. Bit us on
  `LadrunoEmbeddedRebar` (`-host` vs positional `nHost`, and `-kt auto` vs numeric
  `-kt`) — PRs #175→#177; the bug was masked in #175 because that build was broken
  (see the next quirk's CI note) so Zone-A pytest never ran.

- **ladruno auto-merge gates ONLY on the classTag+manifest fast check — NOT the
  Zone-A (Ubuntu) job at all (neither the build nor the pytest).** A PR that does
  not even COMPILE can merge (PR #175 did: a `getInterpolationWeights` override used
  `numberNodes`, a per-method `static const` local in `LadrunoBrick`, not a member).
  A broken ladruno HEAD then makes EVERY later PR's Zone-A red, and since the build
  dies the pytest phase never runs — masking test bugs until someone fixes the
  compile. After pushing C++ to a fork PR, **watch the Zone-A job**
  (`gh pr checks <n> --watch`): a fast (~1-2 min) fail = compile error, a slow
  (~5-6 min) fail = test failure. Don't trust a green fast-gate.

- **`LadrunoEmbeddedNode` is WIDE but only the U+`g0` core is VALIDATED — and
  `getInitialStiff` aliases the D9 tangent.** The element exposes five flag-gated capabilities
  (U · UP · UR · D9 · enforcement), but the [[23_ladruno_embedded_node_adr|ADR §14]] re-scope
  declares **only the U translational tie + `g0` stress-free birth + penalty/AL/bipenalty** as
  the *validated, world-class* core ([[27_ladruno_embedded_node_validation_plan]]). **Do not
  cite UR/UP/D9/`-corot` as validated** — UR is `½curl(u)` SPIN (not moment transfer; rigid
  spin on CST/TET4), UP is niche poromechanics, D9 is interface/contact-flavored (uncoupled
  friction only approximate). Their Zone-A *mechanics* tests prove they run, **not** that
  they're validated. **The one real latent bug:** `getInitialStiff()` aliases
  `getTangentStiff()` → `formTransTraction()` → `setTrialStrain()`, so in **D9 mode** the
  "initial" stiffness is **state-dependent** and **mutates material state during a query**.
  Harmless for the U core (`matMode 0` → `K_u·I`, exact/state-independent) but a real bug that
  **gates D9 promotion** — fix it to use each direction's *initial* tangent with no side
  effect. Also: `sendSelf`/`recvSelf` has **no version field** despite the format changing every
  phase (hdr→29 in #214) — add one (retroactively; pre-#214 DBs already incompatible). 2026-06-07.

- **`Ladruno_scripts\build.bat` takes ONE target argument, not a list.** It reads only
  `%1` (`set "MODE=%1"` → `set "TARGETS=%MODE%"`), so `build.bat OpenSees OpenSeesSP
  OpenSeesMP` builds **only `OpenSees`** and silently ignores `%2 %3 …` — exit code 0, no
  warning. (The `~/.claude/CLAUDE.md` example showing a multi-target list is misleading.)
  To build several targets either run it once per target, or run it with **no arguments**
  (`build.bat` alone builds all five: OpenSees, OpenSeesSP, OpenSeesMP, OpenSeesPy,
  OpenSeesPyMP — incremental via Ninja, so cheap after the first). The Python test module
  is `OpenSeesPy` → `dist\bin\opensees.pyd`; the Tcl exes are `OpenSees/SP/MP.exe`. Symptom
  of the trap: after a "successful" multi-arg build, `dist\bin` has `opensees.pyd` but no
  `OpenSees*.exe`. 2026-06-07.

- **Anisotropic embedded coupling (`LadrunoEmbeddedRebar`) needs a CO-ROTATED bar
  axis under large host rotation; isotropic node ties (`ASDEmbeddedNodeElement`) do
  not.** The frozen reference `dir` is the *only* true large-rotation defect: the gap
  `g` and the host weights `N_i(ξ)` are already frame-objective, but the axial/
  transverse split `s = g·dir`, `g_t = g − s·dir` taken against a FROZEN `dir`
  registers spurious axial slip under pure rigid rotation and yields a non-objective
  traction. Fix (ADR 20 §10.5, `-corot`): recompute `dir` each step as the secant of
  two embed points (embed point + a point B along the bar) from CURRENT host node
  positions. This is why `ASDEmbeddedNodeElement` recomputes geometry from REFERENCE
  coords yet stays objective — its `iK·BᵀB` penalty is isotropic, so there is no axis
  to go stale. (v1 omits the `∂dir/∂u` consistent-tangent term — EICR practice: exact
  for explicit, converges under step-halving for implicit.)

### LadrunoRecorder `-precision f32` is ignored in `-envelope` mode — STORED_PRECISION now honest (FIXED)
- **Bites:** `-precision f32` only changes the dtype of the streaming per-step DATA
  datasets (`StreamingSink::createTimeSeries3d`, `Ladruno_Sinks.cpp` — `H5T_IEEE_F32LE`).
  In `-envelope` mode there are no streaming DATA datasets; the only result datasets are
  the EnvelopeSink MIN/MAX/ABSMAX, which are **always f64**. But `initialize()` stamped
  `INFO/STORED_PRECISION` purely from the `store_data_f32` flag → an `-envelope -precision f32`
  file claimed `f32` while every dataset in it was f64. A reader trusting the attribute to
  pick its diff tolerance would be misled.
- **Fix:** `STORED_PRECISION` is now `f32` only when `store_data_f32 && !envelope_mode`
  (it must describe what is actually on disk); a one-time warning is emitted if `-precision f32`
  is combined with `-envelope`. (Honoring f32 *inside* the envelope datasets is a separate,
  judgment-dependent enhancement — not done; the label-honesty fix is unambiguous.) 2026-06-03.

### LadrunoRecorder `domainChanged()`/`restart()`/`setDomain()` being no-ops is INTENTIONAL — do not "fix"
- **Why it looks wrong:** an adversarial review flagged that these lifecycle hooks are inert,
  so cached `Element*`/`Response*` could dangle after a model edit. **Verified NON-issue:**
  the *only* source-rebuild trigger is the `domain->hasDomainChanged()` stamp checked inside
  `record()` (the `rebuild_model` block) — and **every** structural edit (`addElement`/
  `removeElement`/etc.) bumps the domain's geometry tag, so the stamp moves and the rebuild
  fires, re-acquiring fresh pointers and (re)writing the MODEL_STAGE. This is the exact frozen
  `MPCORecorder::record()` pattern (the code comment says so). Forcing a rebuild in
  `domainChanged()` would be redundant with the stamp check and risk breaking the multi-stage
  logic. **Leave them as no-ops.** (The only genuinely-real lifecycle gap is minor: the `-T`
  frequency gate can stall if `commitTag` regresses across a second `analyze()` after
  `wipeAnalysis` — a defensive guard, not yet added.) 2026-06-03.

### A node-embedding ROTATION tie needs the host's `∂N/∂x`, not its weights `N_i` (ADR 23 Phase 2 / UR)
- **Why it bites:** the translational (U) and pressure (UP) ties only need the host
  shape-function WEIGHTS `N_i(ξ)` (`getInterpolationWeights`, ADR 20). The rotation
  (UR) tie ties the constrained node's rotations to the host CONTINUUM rotation
  `θ = ½ curl(u) = skew(∇u)`, which is built from the host displacement GRADIENT — so
  it needs `∂N_i/∂x` (cartesian shape derivatives), a DIFFERENT host query that
  weights cannot supply. Hence the new vanilla `Element::getInterpolationGradients(ξ,dNdx)`
  (default −1; overridden on `LadrunoBrick` via `shp3d`, `BezierTet10` via
  `computeJacobian`). The translational rows of the UR `B`-matrix still use `N_i`; only
  the rotation rows use `∂N/∂x`.
- **Volume host ⇒ PURE `skew(∇u)`, not ASD's mixed convention.** `ASDEmbeddedNodeElement`
  embeds into a planar tri/tet *surface*, so it builds a 2D local frame and uses the
  surface SLOPE (factor 1) for the two bending rotations + `½ curl` (factor ½) only for
  the drilling — a mixed convention forced by the missing out-of-plane derivative.
  `LadrunoEmbeddedNode` embeds into a 3D VOLUME host (hex/tet) where all 9 gradient
  components are available, so it uses the dimensionally-clean **pure continuum rotation**
  `θ = ½ Σ_i (∇N_i × u_i)` (½ on all three, no local frame, frame-objective) — de Souza
  Neto §3. The host operator block is `½·skew(∇N_i)`; the gradient virtual returns global
  cartesian `∂N/∂x` directly, so NO `R`-rotation of the block is needed.
- **UR is mesh-limited (UR-4):** on a CST (3-node tri) / TET4 (4-node tet) host `∂N/∂x`
  is element-CONSTANT ⇒ the UR constraint collapses to a single element-wide RIGID-SPIN
  tie (no intra-element rotational gradient). Moment-critical embeds (anchors, headed
  studs) need a higher-order host (`BezierTet10`) where `∂N/∂x` varies with ξ. Document,
  don't silently sell as exact. 2026-06-04.

### The node-embed element needs `-corot` ONLY for the D9 MATERIAL frame, never for the penalty tie (ADR 23 Phase 2b v2)
- **The split:** the isotropic/penalty U/UP/UR tie is already frame-objective (`D=K·I` has no
  preferred axis; gap + weights both transform with the host), so it needs **no** co-rotation —
  unlike the anisotropic `LadrunoEmbeddedRebar`, whose frozen bar `dir` goes stale. But the D9
  **material interface** reintroduces a preferred axis (the `-normal`/tangent frame carrying
  per-direction uniaxials), so a directional contact normal DOES go stale under host rotation.
  `-corot` co-rotates that frame with the host CONTINUUM rotation `θ_host = skew(∇u)|_ξ`,
  `frameCur = R(θ_host)·frame` — **reusing the UR `∂N/∂x`/`rotOper` machinery verbatim** (not the
  rebar's secant-to-point-B trick; the node element has a normal+tangents frame, not a bar axis,
  so the natural rotation source is the host continuum spin, factored into `hostContinuumRotation`).
  3D = Rodrigues exp-map of the axial vector `θ_host`; 2D = the single drilling planar rotation.
- **Mechanically distinct from UR:** UR ties the cNode's rotation DOFs to `θ_host` (a constraint
  on a DOF); `-corot` rotates the *material frame* used by the translational interface (no rotation
  DOF needed — material mode runs at ndf=ndm). They share only `θ_host`/`gradN`; `-corot` is
  material-mode-only (parse-time reject otherwise) and independent of `-rot`.
- **The dropped `∂e_d/∂u` tangent term is EXACT, not approximate, when the host DOFs are
  prescribed.** `-corot` inherits the rebar's dropped consistent-tangent term (R7/D9-5: residual
  exact, tangent inexact ⇒ frame-objective for explicit, step-halving for implicit, may slow
  Newton on stiff-normal large-per-step-rotation contact). But note `frameCur` depends on
  `θ_host`, which is a function of the HOST translations **only** (`∂frameCur/∂u_cNode = 0`). So in
  a Zone-A mechanics test where the host is prescribed (`sp`/`fix`), the cNode tangent is *exact*
  and Newton converges quadratically — the inexactness only bites when the host continuum is free
  and spinning per-step. 2026-06-07.

### LadrunoEmbeddedNode v1 dropped the parent `m_U0` offset-capture → absolute tie yanks on staged addition (FIXED, ADR 23 Phase 2c)
- **The bug:** v1 computed every gap as a pure TRIAL-DISPLACEMENT difference
  (`g = u_c − Σ N_i u_host`, kernel `LadrunoEmbedded::computeGap`; likewise `g_p`, `g_r`), so
  the penalty enforced an ABSOLUTE tie `u_c = Σ N_i u_host`. An element added MID-ANALYSIS to a
  host that has already deformed (staged construction) activates with `g = −N·u_host ≠ 0` and the
  penalty **yanks the slave** by the full accumulated host displacement — a spurious force spike.
  The parent `ASDEmbeddedNodeElement` (and `equalDOF_Mixed`) already capture this offset
  (`m_U0` snapshot at `setDomain`, `getGlobalDisplacements()` returns `U − m_U0`); the fork's v1
  port silently dropped it.
- **The fix:** at `setDomain` capture each ACTIVE gap ONCE (`g0`/`gp0`/`gr0`, guarded by
  `g0Computed`) and drive ALL traction from the RELATIVE gap `(g − g0)`. Subtract the offset
  **inside** `computeGap`/`computeGapP`/`computeGapR` (NOT at each call site) so every consumer
  is covered in one place.
- **Force-free ≠ stress-free — the trap.** Zeroing only the penalty force is NOT enough in the
  D9 material mode: the gap also drives `matDir[d]->setTrialStrain(g·e_d)`. If the offset is an
  additive force correction, the material still sees the ABSOLUTE gap and is born PRE-STRAINED
  (a cohesive law partway up its backbone, a gap material already closed, bond pre-slipped) —
  force-corrected but NOT stress-free. Subtracting `g0` inside the gap (so the material's strain
  ORIGIN shifts) is what makes it genuinely stress-free. This is why "shift the canonical gap"
  beats "subtract at each consumer."
- **Default ON; no-op when undeformed.** Capture is ON by default (restores parent behavior);
  when the element is added at the undeformed state `g0 = 0` ⇒ byte-identical to the absolute
  tie, so the whole v1 battery is unaffected. `-absolute` (alias `-noInitGap`) opts out (legacy
  tie / a deliberate snap-to-host). `g0Computed` is serialized so `recvSelf` restores the
  captured offset instead of re-capturing. UR is linearized ⇒ `gr0` subtraction is exact for
  small inter-stage rotation, approximate for large. 2026-06-07.

### Per-DOF-class bipenalty: a translational `m_p` CANNOT bound the rotation mode (ADR 23 M1/ES-1)
- **Why it bites:** the bipenalty mass penalty `m_p` (lumped on the slave's translational
  DOFs) bounds the explicit `dt_cr` of the TRANSLATIONAL coupling only. The rotation tie's
  penalty `K_r` has different units (moment/rotation), so a translational-only `m_p` leaves
  the rotation mode UNBOUNDED in explicit (`dt_cr → 0`). Fix: give the rotation class its
  OWN inertia `I_p = K_r·(dt/2)²` (the SAME `-dtcr`/`-wcap` budget formula but with `K_r`),
  lumped on the slave's ROTATION DOFs. Then `dt_r = 2√(I_p/K_r) = dt_u` and the `lch²` in
  `K_r` cancels (it's also in `I_p`), so the rotation mode self-bounds at the SAME `dt`.
  `getExplicitCriticalTimeStep` reports the MIN over active DOF classes. (The same pattern
  generalizes to a pressure class if pressure bipenalty is ever added — pressure is
  implicit-recommended for now.) 2026-06-04.

### D9 interface material returns FORCE, not stress — no `bondScale` (unlike the rebar)
- **Why it bites:** `LadrunoEmbeddedRebar`'s axial slot drives its bond material with the
  axial SLIP and the material works in STRESS units (τ–s), so the element multiplies by
  `bondScale = perimeter·L_trib` to get a nodal force. `LadrunoEmbeddedNode`'s D9 interface
  materials are driven by the displacement GAP `g·e_d` (metres) and are expected to RETURN
  FORCE directly (`stress()` in N) — so there is **NO bondScale converter**: `t_d =
  mat_d->getStress()`, `k_d = mat_d->getTangent()` go straight into `t=Σ t_d e_d`,
  `D=Σ k_d e_d⊗e_d`. Pick/define the uniaxial accordingly (an Elastic of "stiffness" K is
  a penalty of force-per-metre K; a cohesive law's peak is a force, not a traction). The
  D9 grammar confines `-mat*` to the TRANSLATIONAL normal/tangent directions, so the
  rotation-unit (M1/ES-1) problem never arises in a material slot. 2026-06-04.
- **v1 uses the REFERENCE frame; `-corot` frame co-rotation is DEFERRED.** A material on a
  specific direction (esp. a unilateral-contact normal) ideally co-rotates with the host;
  v1 keeps the frame fixed (valid for small-rotation interfaces). The v2 corot would reuse
  the Phase-2 continuum-rotation (`skew(∇u)` from host gradients) to rotate the frame — and
  would re-introduce the rebar's dropped `∂e_d/∂u` consistent-tangent caveat (frame-objective
  for explicit, converges under step-halving for implicit, may slow Newton for stiff-normal
  large-rotation contact). Large-rotation RIGOROUS contact is the separate `LadrunoContact`.

### `UniformExcitation` writes to `theDof` with NO per-node ndf bounds check (mixed-ndf footgun)
- **Why it bites:** `UniformExcitation::applyLoad` calls `theNode->setR(theDof, 0, fact)` for every
  node in the domain, guarding only on `ndm` (`theDof < 1/2/3`), NEVER on the node's actual
  `numberDOF` (`SRC/domain/pattern/UniformExcitation.cpp:303-365`, writes at `:318/:323/:335`).
  In a MIXED-ndf model a single `UniformExcitation` hits ALL nodes, so exciting e.g. `dof 2`
  is fine on the ndf>=3 nodes but writes OUT OF BOUNDS on any ndf=2 (plane-solid) node sharing
  the domain. No warning, no skip. Apply ground motion only to a node set that actually owns
  the dof, or fix the affected nodes out of the excited direction. See
  [[ndf_and_mixed_models_guide]] §6. 2026-06-07.

### Explicit zero-mass rotational DOF on an ndf=6 node ⇒ `dt_cr` silently overestimated
- **Why it bites:** lumped element mass (beam lumped, `ASDShellQ4` rotational mass is
  EXPLICITLY omitted, `ASDShellQ4.cpp:1152`) and translational-only nodal `-mass` leave ZERO
  mass on the rotational dofs of an ndf=6 node ⇒ singular `M`. `CriticalTimeStep` does NOT
  warn — its DGGEV path FILTERS near-massless eigenpairs via a relative beta threshold
  (`betaTol = 1e-12*max|beta|`, `SRC/analysis/integrator/CriticalTimeStep.cpp:165-170`), so
  the zero-mass mode is dropped from omega_max and `dt_cr = 2/omega_max` comes back
  UNCONSERVATIVELY LARGE ⇒ the explicit run can go unstable with no diagnostic. In a mixed-ndf
  explicit model the binding constraint is MASS not stiffness: give every ACTIVE dof (incl.
  rotations on ndf=6 nodes) nonzero mass, use consistent mass, or restrain the massless dofs.
  Relevant to [[central_difference_ladruno_guide]]. See [[ndf_and_mixed_models_guide]] §7.
  2026-06-07.

### A no-op `setRayleighDampingFactors` WITHOUT a `getDamp` override ⇒ hard crash in implicit transient
- **Why it bites:** a pure-coupling/penalty element often overrides
  `setRayleighDampingFactors(...)` to a no-op `return 0` to refuse Rayleigh damping (so a
  `betaK` can't spuriously shrink its explicit `dt_cr`). But the base `Element::getDamp()`
  (`SRC/element/Element.cpp:211`) does `if (index==-1) this->setRayleighDampingFactors(...)`
  then `theMatrices[index]->Zero()` — and it is the BASE `setRayleighDampingFactors` that
  lazily allocates `theMatrices[index]` and sets `index>=0`. A no-op override never
  allocates, so `index` stays at its ctor default −1 and `theMatrices[-1]` is an
  out-of-bounds dereference → **hard crash**. `FE_Element::addCtoTang(fact)` calls
  `getDamp()` whenever `fact!=0`, and the Newmark/HHT velocity coefficient
  `c2=γ/(βΔt)` is ALWAYS nonzero ⇒ getDamp fires in EVERY implicit transient step. The
  residual damping-force path `addD_Force` calls it too. `getMass`,
  `getResistingForceIncInertia`, and `getRayleighDampingForces` share the same
  `theMatrices[index]` landmine — the first two are usually already overridden (so safe);
  **`getDamp` (and `getRayleighDampingForces`) are the ones people forget.**
- **Why quasi-static tests miss it:** `LoadControl`/`Static` never form a C-tangent;
  `CentralDifference` (explicit) dodges it when the model has no Rayleigh. Only an implicit
  transient (Newmark/HHT) — or any transient with Rayleigh — triggers it.
- **Fix:** override `getDamp()` and `getRayleighDampingForces()` to return an element-owned
  ZEROED `Matrix`/`Vector` (sized `nDOF`, allocated alongside the mass matrix), bypassing the
  base index path. `D≡0` is physically correct for a pure coupling; mass/inertia still come
  from `getMass` + bipenalty. Confirmed: `LadrunoDistributingCoupling` (RBE3, 33011) crashed
  a Newmark transient (exit 5) before the override, passes after (regression test added).
  **`LadrunoEmbeddedNode` (33006) + `LadrunoEmbeddedRebar` (33005) had the SAME latent bug**
  (no-op setRayleigh, no getDamp override) — **FIXED 2026-06-09** with the identical
  element-owned zeroed `C0`/`dampF` pattern + a `Newmark 0.5 0.25` `test_transient_newmark_smoke`
  regression in each Zone-A battery (empirically reproduced: pre-fix the smoke test segfaults
  `0xC0000005`, post-fix 77/77 pass). See [[LEDGER_implementations]] rows 33005/33006, PR #220.
  2026-06-07 (RBE3) / 2026-06-09 (embedded).

### `LadrunoArcLength -stabilize` (33004): what viscous regularization can and cannot pass
Measured on the live build (2026-06-16) while building the ADR-31 rung-4 seam. Four
non-obvious behaviours, all relevant to anyone wiring `-stabilize` into a driver:
- **Pure monotone softening is NOT passable by `-stabilize`.** In stabilize mode the
  integrator IS load control (mutually exclusive with the arc-length quadratic), and a
  softening branch has no equilibrium above the peak. Stabilized load control stalls at
  the strength peak *exactly* like plain load control (softening Concrete02 truss: both
  stop at `λ=30.0, ε=−0.002`). Softening is rung-3's job (switch to displacement control),
  NOT rung-4's. `-stabilize` only helps **snap-through-style limit points** where a
  continuing branch exists.
- **`-adaptStab` PREVENTS crossing a hard limit point.** It rescales `cVisc` each commit
  to hold `dissipVisc/Estrain0 ≈ fTarget`, which keeps the viscous force too weak to push
  through. On the von Mises snap-through truss, `-stabilize -adaptStab` (any `f`) stalls at
  the limit; only `-stabilize` *without* adaptStab and at an **elevated** `f` (≈1e-3…1e-2)
  crosses — and then by a slow **diffusive crawl** (ran to the 2000-step budget, cumulative
  ratio ~3e3), i.e. the "dynamic jump across the unstable branch" ADR-31 flags as R-LOG-MASK,
  not a traced path. Crossing is also **non-monotonic in `f`**: `f`=1e-3,1e-2 cross but
  `f`=5e-2,1e-1,5e-1 stall again (too much damping re-freezes the step).
- **`-cVisc N` is silently overwritten by the first-commit calibration.** `commit()` runs
  the `fTarget` calibration under `if (!cCalibrated)` with NO guard for a user-supplied
  `cVisc`, so `-cVisc 1`, `10`, `100` all produce the identical run (ratio 0.4495). The code
  comment "one-shot calibration … (skipped if -cVisc given)" is aspirational — the skip is
  not implemented. If you need an explicit coefficient today, the `-cVisc` path does not
  deliver it. (Out of scope to fix under ADR-31; flagged for a future ADR-20 follow-up.)
- **openseespy turns a command's stderr WARNING into a raised `OpenSeesError`.** The
  `scaleCVisc(factor)` guard `factor>0` writes a warning to `opserr` and returns −1, but in
  the Python module any `opserr` output during a command raises `opensees.OpenSeesError`
  rather than returning the −1. A driver must pre-validate `factor>0` itself (never rely on
  catching the −1) or wrap the call in try/except. (Tcl sees the −1 return; Python does not.)
  Verified by `torture_stabilize.py` + `test_robust_battery.py::test_stabilize_*` (4 cases).
- **RC-shell Phase 2a interlock — the MCFT `v_ci,max` formula is UNIT-DEPENDENT (SI: N, mm).**
  `LadrunoRCConcrete -interlock` bounds the crack-plane shear at the Vecchio–Collins limit
  `v_ci,max = 0.18·√fc' / (0.31 + 24·w/(a_g+16))` with crack width `w = eps_n·s_theta`. The
  numeric constants (0.18, 0.31, 24, 16) are empirical in **MPa and mm** — `√fc'` is `√(MPa)`,
  `w` and `a_g` are in **mm**. The rest of the kernel is unit-agnostic, but this one law is not:
  use it on an N–mm–tonne–s model (fc' in MPa = N/mm², lengths in mm) or rescale the constants.
  `s_theta` defaults to `-crackSpacing`, else `lch`, else 1.0; `a_g` (`-agg`) defaults to 16 mm.
  NB Phase 2a **CLIPS the smeared (damage-reduced) crack-plane shear `τ_sm=m_σ·sig_ip` to ±v_ci,max**
  (a bound), it does NOT substitute bare-elastic `G·γ`; below the cap the stress is unchanged.
- **Fixed-crack interlock only engages under NON-PROPORTIONAL loading.** The crack normal freezes to
  the principal-tensile direction at cracking; under a *proportional* path stress/strain stay coaxial,
  so the crack-plane shear is ~0 and `v_ci,max` never binds — interlock looks inert. It engages only
  once the principal direction ROTATES off the frozen normal (tension-then-shear, any non-radial path).
  TEST consequence: an off-axis interlock test MUST be two-stage (freeze oblique, then rotate shear onto
  it); a single proportional ramp gives `τ_nt≈0` and tests nothing.
- **rc_shell_ref.py oracle uses RAW `(x,y,q)` backbones; the C++ kernel ADJUSTS them
  (`buildBackbone` E-consistency).** So oracle-vs-C++ **absolute** stress only matches for
  quantities INDEPENDENT of the backbone q-adjustment. The Phase-1 β gate dodges this by being
  a RATIO; the Phase-2a interlock cap test dodges it by asserting only the **CAPPED** crack-plane shear
  `|τ_nt| == v_ci,max` (backbone-free); sub-cap `σ_xy` is the damaged smeared shear and differs
  C++-vs-oracle. To compare absolute normal stresses oracle-vs-C++ you must first port
  `buildBackbone`'s adjust() into the numpy Backbone. (Two-stage driver `_path_tension_then_shear`:
  tension on dof-1∝X, shear on dof-2∝X — disjoint DOFs let one SP per DOF realize the path under Penalty;
  pass `gamma0>0` to freeze an OBLIQUE crack for the off-axis rotation test.)

### `timeSeries Path` returns 0 BEYOND its last time node — float-accumulated pseudo-time overshoots and collapses prescribed strains
- **Bites:** a multi-stage prescribed-strain cyclic driver built from `timeSeries('Path', ...)` + `LoadControl(1/nper)` over N stages. The intended end pseudo-time is `N`, but `nper×N` accumulations of `1/nper` (e.g. `1/80`, not exact in binary) land at `N + epsilon`. `Path` returns **0** outside `[t_first, t_last]`, so at the FINAL step every `sp`-prescribed DOF drops to 0 → the element snaps to ~zero strain. For the RC cyclic interlock this looked like a phantom `-3.18` crack-shear spike (crack "closed" at `en≈0` ⇒ `v_ci,max` jumps to its max `0.18√fc/0.31`) — a TEST artifact, not a kernel bug.
- **Why:** `PathSeries::getFactor(t)` returns 0 for `t > t_last` (and `< t_first`). The classic mis-diagnosis is "the cap is wrong"; the real tell is reading the GP strain at the offending step (it's ~0, not the held value).
- **Fix (robust):** pad the Path with an extra HOLD node beyond the analysis end — `times=[0..N, N+1]`, `values=[...,last, last]` — so any overshoot interpolates between two equal final values. (`-useLast` as a trailing openseespy arg did NOT take effect in this build; the pad is reliable.) Learned 2026-06-16 building the Phase-2b cyclic shear driver `_path_stages`.

### RC fixed-crack cyclic PINCHING is a panel phenomenon, not a material-point one
- **Bites:** expecting a pinched (waisted) `τ_nt`–`γ_nt` hysteresis loop from a single material point / one element under homogeneous cyclic shear at constant normal strain. You get a FAT loop instead.
- **Why:** with the crack frozen and `en` (hence `v_ci,max`) constant, the only nonlinearity is the `±v_ci,max` clamp; the elastic interlock band has width `2·v_ci,max/G ≈ 5e-4` in slip — sub-step at normal resolution — so the loop is essentially a `±v_ci,max` rectangle (maximum dissipation, zero pinch). A pinched waist requires the crack to OPEN/CLOSE during the cycle (`v_ci,max` low near slip-reversal, high at slip-peaks), which needs the principal direction to ROTATE relative to the fixed crack — i.e. a real panel / non-homogeneous stress field. So Phase-2b.1 material-point tests assert the MECHANISM (reversal re-cap, unload=G, closure cap-recovery, energy>0), and the pinching-shape + hysteretic-energy acceptance is a panel/experiment (Tran–Wallace) gate deferred to 2b.2. Learned 2026-06-16.

### IMPL-EX in a STATIC analysis: `ops_Dt` (load-factor pseudo-time) is erratic — guard the extrapolation time-factor
- **Bites:** porting the ASDConcrete3D IMPL-EX recipe (`tf = dtime_n / dtime_n_commit * alpha`, `dtime_n = ops_Dt`) into a material and running it under a STATIC `DisplacementControl`/`LoadControl` analysis (e.g. a quasi-static cyclic wall). The IMPL-EX extrapolation `x_ext = x_n + tf·(x_n − x_{n-1})` detonates — damage jumps to garbage, the first step diverges immediately — even though the same material is fine in dynamics.
- **Why:** in a STATIC analysis the domain "time" IS the load factor λ, not a physical Δt. With `DisplacementControl`, λ's increment is whatever satisfies the controlled DOF (can be ~1e5 for a stiff structure), and `loadConst('-time',0.0)` (the standard gravity-then-pushover idiom) RESETS λ to 0 mid-run. So `ops_Dt = λ_n − λ_{n-1}` is huge / tiny / negative across steps, and `tf = ops_Dt/ops_Dt_commit` becomes a wild multiplier on the threshold extrapolation. ASDConcrete3D ships no clamp because it is typically exercised in dynamics (or with `-dtime`/user-defined Δt) where Δt is smooth.
- **Fix (proven):** clamp the time-factor in `implexTimeFactor()`: fall back to `alpha` whenever `!commitDone` or `dtime_n<=0` or `dtime_n_commit<=0` or the ratio is non-finite; clamp the result to `[_, 2·alpha]` so a single pseudo-time spike cannot blow up the extrapolation. For static IMPL-EX this effectively makes `tf≈alpha` (uniform extrapolation) — correct, since static load steps are meant to be uniform. Learned 2026-06-17 pulling Phase-4 IMPL-EX forward for `LadrunoRCConcrete`.

### Build gotcha: `cp` (git-bash) / `Copy-Item` to the main checkout can silently skip the rebuild
- **Bites:** editing fork source in the WORKTREE, copying to the main checkout, running `build.bat`, and testing — but the binary still shows OLD behavior (a new `-flag` is silently ignored, a new response returns empty). Two distinct traps stack: (1) `cp "src" "C:\…\dst"` from the Bash tool can mis-resolve the Windows backslash destination and write nowhere useful (the file you think you copied is unchanged in the main checkout — `grep -c <newsymbol>` there returns 0); (2) even after a correct `Copy-Item`, ninja compares mtimes and **`Copy-Item` PRESERVES the source's (older) mtime**, so if the worktree file was edited before the last build's `.obj`, ninja sees the object as newer and SKIPS recompiling — `opensees.pyd`'s timestamp never advances.
- **Tells:** the built `.pyd` LastWriteTime does not change after a "successful" build; `grep -c "<your new symbol>" <main-checkout-source>` returns 0; a parser silently ignores your new token (the RC parser treats unknown tokens as no-ops — forward-compat — so a non-compiled flag fails OPEN, not loud).
- **Fix (proven):** copy via `Copy-Item` (reliable on Windows), then explicitly bump the destination mtime `(Get-Item $dst).LastWriteTime = Get-Date` (and the `.cpp` that `#include`s a changed header) BEFORE `build.bat`, and confirm the build log shows `Building CXX object …<file>.cpp.obj` + `copying OpenSeesPy.dll -> opensees.pyd` and that the `.pyd` timestamp advanced. Run `build.bat` via the PowerShell tool, not the Bash heredoc (the latter captured only the banner here). Learned 2026-06-17.

### Cyclic softening RC shell wall (ASDShellQ4 + LadrunoRCConcrete) walls on Newton convergence past first crack
- **Bites:** building a meshed reinforced RC shear wall — `ASDShellQ4` grid on a `LayeredShell`
  (`LadrunoRCConcrete` concrete layers + `PlateRebar(Steel02)` web steel) — and pushing it cyclically
  under `DisplacementControl` to demonstrate panel-scale pinching. It assembles fine and runs the first
  small-drift cycle (real `V`–`δ` hysteresis with dissipation), then `analyze` returns `-3`
  (`NormDispIncr` stalls with a large residual `deltaR`) at the next, larger amplitude — even with
  `-implex`, `KrylovNewton`, and `NewtonLineSearch`/`ModifiedNewton`/`Broyden` fallbacks.
- **Why:** the softening plastic-damage + crack-localization makes the global tangent indefinite/ill-
  conditioned at the load-redistribution events (crack formation, interlock cap engaging across a row of
  elements); a load-/displacement-controlled Newton has no way around the limit/snap-back points. IMPL-EX
  helps the MATERIAL tangent stay SPD-secant but does not fix the STRUCTURAL snap-through. This is the
  textbook reason squat-wall validation is hard, not a bug in the material (the material-point + single-
  shell-element cyclic gates all pass).
- **Fix (the deferred path, not yet done):** drive with an arc-length / indirect-displacement control
  (`LadrunoIndirectControl` / `LadrunoArcLength` are built for exactly this snap-back), or dynamic
  relaxation / quasi-static transient with mass; finer substeps; possibly an IMPL-EX-error step-cut.
  Harness scaffold: `tests/_testbed/rc_wall_harness.py`. Learned 2026-06-17 building the Phase-2b.2c.4
  Tran–Wallace pinching validation for [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]].

### Quasi-static EXPLICIT is the cyclic-softening tool; the monotonic solvers are NOT — and ASDShellQ4 has no per-element dt_cr
- **Bites:** trying to validate a CYCLIC softening RC wall (RC shell stack 2b.2c). Implicit Newton diverges at the first reversal (the plastic-damage consistent tangent is indefinite on the softening branch). Reaching for the fork's softening solvers (`LadrunoArcLength`, `LadrunoIndirectControl`, `LadrunoDynamicRelaxation`, the `robust_drive` rung ladder) does NOT help — they are all MONOTONIC tools (they trace ONE equilibrium path through a limit point / snap-back); a load reversal is not a single monotonic path, so arc-length/indirect-control/DR have nothing to follow.
- **Fix (proven):** drive the wall as a QUASI-STATIC EXPLICIT problem — `integrator CentralDifferenceLadruno` (`-cflAbort -lump diagonal`) + `system Diagonal` + `algorithm Linear`. Explicit forms NO stiffness tangent, so the indefinite-tangent stall simply does not exist; reversals + softening integrate through. A single ASDShellQ4 panel AND a 4×3 gmsh-meshed squat wall both completed the full ±drift schedule this way (where implicit diverged immediately). This is how LS-DYNA et al. run cyclic RC walls.
- **Sub-gotchas for explicit + ASDShellQ4 (each cost real time):**
  1. **Need ELEMENT mass via the material `-rho`** (e.g. `-rho 2.4e-9` tonne/mm³). Nodal `ops.mass(...)` is NOT enough — the dt_cr eigensolve (and a clean mass matrix) needs element/section mass; with only nodal mass you get "no element produced a finite, positive estimate."
  2. **ASDShellQ4 supplies NO per-element critical-time-step estimate** — `ops.criticalTimeStep()` returns **-1** ("no element produced a finite estimate"). `CentralDifferenceLadruno`'s `CriticalTimeStep` eigensolve gets nothing usable from the shell. FALL BACK to a manual wave-speed bound: `dt ≈ frac · h / sqrt(E/rho)` (frac ~0.2; with E=30000 MPa, rho=2.4e-9, h=500 mm → c≈3.5e6 mm/s, dt_cr≈1.4e-4, dt≈2.8e-5). Querying `criticalTimeStep()` and blindly using it gives `dt = 0.2·(-1) < 0` → instant blow-up.
  3. **NO `equalDOF` / rigid ties** — `dt_cr` and the CD stability bound IGNORE constraints (documented in [[central_difference_ladruno_guide]]), so a rigid top via `equalDOF` makes the stability estimate a lie and the run detonates. Prescribe the rigid-top drift by putting the SAME `sp` value on every top node instead.
  4. **Mass-proportional damping only** — stiffness-proportional (`betaK`) Rayleigh collapses `dt_cr` ~quadratically.
  5. **Quasi-static = loading period ≫ structure period** (use a smooth cosine drift on a dt-based time grid, segment_time = steps_per_seg·dt with steps_per_seg ~ a few thousand). NOTE step count for quasi-static is ~independent of mass scaling (T and dt scale together).
- **Panel-scale finding:** at the SINGLE-element scale the cyclic interlock terms are nearly inert (concrete-damage-dominated; cyclic≈monotone). At the MESH scale the cyclic friction-slip + X-crack wear become clearly load-bearing — a 4×3 squat wall dissipated ~28% LESS hysteretic energy and ~11% lower peak shear with `-cyclic -xcrack` vs the monotone `-interlock` bound. The pinched-loop *waist* sharpens with the mesh (principal rotation); a quantitative experiment (Tran–Wallace) match is a further calibration step. Tests: `tests/test_ladrunoRCConcrete_wall.py` (Zone-B). Learned 2026-06-17.

### Tension-stiffening rank-1 floor: the equibiaxial self-consistency coefficient is 0.5, not 1
- **Bites:** the RC-shell tension-stiffening floor (RC stack Phase 3a, `-tensStiff`) pins the principal tensile stress to `σ_ts(ε1)` by injecting `Δ·(p1⊗p1)` and measuring `n^Tσn`. For a true unit eigenvector the self-consistency coefficient `a²+b²+2(ab)² = (p1x²+p1y²)² = 1`, so one injection pins the normal exactly. In the DEGENERATE (equibiaxial) branch the natural fallback `(a,b,ab)=(0.5,0.5,0)` is **rank-2 (isotropic), not rank-1**, so its coefficient is `0.5²+0.5² = 0.5` — injecting `Δ·(0.5,0.5,0)` raises the mean normal by only `0.5·Δ`, so the floor reaches **halfway** to `σ_ts`, never `σ_ts`.
- **Fix (proven):** in the degenerate branch use distinct injection vs measurement vectors — inject `g=σ_ts−mean` to BOTH in-plane normals with full weight (`ts_inj=(1,1,0)`), measure the mean (`ts_meas=(0.5,0.5,0)`); then `ts_meas·ts_inj=1` and each in-plane normal reaches `σ_ts` exactly (verified to 1e-16 in the standalone g++ equibiaxial gate). General rule: the floor's self-consistency requires `ts_meas·ts_inj==1`, which holds for a rank-1 `p⊗p` but NOT for the rank-2 isotropic `½(I)`. A pure uniaxial gate (the original T1) never exercises the degenerate branch, so this only surfaces under an EQUIBIAXIAL test. Learned 2026-06-18 (caught by a 3-agent adversarial review of [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a).

### Tension-stiffening pinning tangent must span ALL 6 columns of the floored rows, not just in-plane
- **Bites:** the tension-stiffening floor pins `n^Tσn` to `σ_ts(ε1)` independent of the bare stress. The bare in-plane normal stress depends on **out-of-plane** strain too (eps_zz via the elastic λ coupling), so pinning it removes that dependence. A consistent tangent that subtracts the bare sensitivity `d(n^Tσn)/dε` only over the in-plane columns `{0,1,3}` leaves the **eps_zz column** (and any out-of-plane shear columns) at the now-pinned-away bare value → a forward-difference tangent check shows ~0.98 rel-error at `D[*][2]`.
- **Fix (proven):** compute `row[c] = ts_meas·Dtan[in-plane rows][c]` for ALL c in 0..5 and subtract over all 6 columns; the `dσ_ts/dε1·de1/dε` add-back is membrane-only (nonzero only on `{0,1,3}`). Note the BASELINE W_B secant tangent is itself only approximate in the softening regime (it omits `d(dt_bar)/dε`), so a global FD-tangent gate "fails" for both on and off; the right TS tangent checks are (a) the PINNED direction `D[0][0]==dσ_ts/dε1` (exact, since `σ0=σ_ts(ε1)` is smooth) and (b) TS does not DEGRADE the worst FD rel-error vs the off baseline. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a (`tests/_testbed/rc_tensstiff_gpp.cpp`).

### Tension stiffening is monotonic-scope (live-ε1 floor re-inflates on unload); TS+interlock mixes live-p1 with frozen-crack
- **Bites:** `-tensStiff` floors the stress to `σ_ts(ε1)` using the LIVE membrane principal strain `ε1` (no `ε1max` memory). Because `σ_ts` *decreases* with `ε1`, on UNLOADING the floor **re-inflates** (tracks `σ_ts(live ε1)` back UP) — a load→unload analysis shows the floored stress rising as strain drops. This is correct on a monotone loading branch but is NOT a hysteretic cyclic-tension model. Separately, TS uses the LIVE principal axis `p1` while the fixed-crack `-interlock` uses the FROZEN crack normal; once principal axes rotate, the TS normal-stress injection leaks a small shear `Δ·sin θ_rel cos θ_rel` onto the frozen crack plane that the interlock then bounds.
- **Fix / scope:** use `-tensStiff` for MONOTONIC / pushover analyses; combined TS+interlock is validated for PROPORTIONAL (non-rotating) loading only. The cyclic upgrade (deferred) is an `ε1max`-envelope floor + secant unload + flooring along the frozen crack plane when cracked. Documented in the kernel comment + [[LadrunoRCConcrete_guide]] §4.7. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a.

### Crack-band regularization of a TABULATED backbone must re-apply adjust() after scaling
- **Bites:** cloning ASDConcrete3D's `HardeningLaw::regularize` (RC stack Phase 3b, `-autoRegularization`). The post-peak strain is rescaled to hit the target fracture energy `g_reg=G_f0·(lch_ref/lch)`. ASDConcrete3D calls `adjust()` INSIDE the scaling loop and once after — easy to omit because the **fracture energy comes out correct without it** (g uses only x,y). But the stress-update reads `q`/damage, and stretching x while preserving the plastic-to-inelastic ratio can drive the tail plastic strain `x−q/E` BACKWARD (non-physical) → wrong `q`/`dt_bar`/`dc_bar` (~6% on steep-mid-softening-damage backbones). Energy-objectivity gates (and the numpy/g++ energy gates) CANNOT catch it.
- **Fix (proven):** mirror the reference — call `adjustBackbone()` (E-cap + monotone plastic strain + non-decreasing damage + `q=y/(1-d)`) after the q/d update each iteration AND once after the loop. Gate it with a steep-damage backbone asserting the post-peak plastic strain is monotone non-decreasing (the energy gate alone is blind). Also: `adjustBackbone` derives `d` from the stored `q`, so the FIRST point must have `q[0]=0` (⇒ d0=0), matching what `buildBackbone` always produces — a hand-built oracle backbone with `q[0]=1e-12` makes adjust read `d0=1` and corrupts everything. Learned 2026-06-18 (3-agent review of [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3b).

### nDMaterial that regularizes in place: getCopy must PROPAGATE the latch, and a loud-fail must NOT latch
- **Bites:** Phase 3b regularizes the backbone IN PLACE (mutates P.ht/P.hc) on a one-time `regularizationDone` latch at first setTrialStrain. Two lifecycle traps: (1) `getCopy` that reconstructs from `P` while RESETTING `regularizationDone=false` will, for a copy made AFTER the source regularized, carry the already-stretched backbone + an un-latched flag ⇒ it regularizes AGAIN (double-scaled, silently). ASDConcrete3D dodges this by using its COPY-CTOR (copies `regularization_done`). (2) A "loud failure" (no lch available) that sets `regularizationDone=true` before returning −1 disarms itself: a `-1` from setTrialStrain is a RECOVERABLE convergence signal, so the algorithm cuts the step and RETRIES → the retry finds it "done" and silently proceeds on the UN-regularized backbone — the exact silent-fallback the guard forbids.
- **Fix (proven):** `getCopy` must copy `regularizationDone`/`regLch` to the new instance (mirror the reference copy-ctor) so a copy-of-a-used-instance is safe. The loud-fail branch must NOT set `regularizationDone` (so it re-fails every step and can never silently continue); gate the message to once with a separate transient `regWarned` flag. In the normal host flow getCopy runs before any setTrialStrain so (1) is dodged anyway, but the safety must be explicit not incidental. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3b.

### Structural crack-band objectivity (Bažant bar): localize via a thinner SECTION, and the un-regularized FINE mesh snap-back-diverges
- **Bites:** building the *structural* mesh-objectivity gate for `-autoRegularization` (RC stack Phase 3b-struct, `tests/test_ladrunoRCConcrete_meshobj.py`) — a row of `ASDShellQ4` driven to full tensile softening, the same specimen meshed coarse/medium/fine, asserting the total dissipated energy is mesh-objective. Two traps. (1) **How to localize the crack to one band:** a *weaker material* (lower `ft`) in the band changes the backbone/fracture energy, so the band-dissipation factor no longer cancels cleanly between meshes. Localize instead with a **10%-thinner SECTION** (one `LayeredShell` of reduced thickness on the left element column) — material, backbone and regularization stay identical everywhere; only the cross-section triggers the band, so the thickness factor cancels and the energy comparison is clean. (2) **The un-regularized control is numerically pathological at fine mesh:** without regularization the one-element-wide band is BRITTLE, and as the mesh refines the band's softening tail shortens until the series elastic-unloading outruns it and the response SNAP-BACKS → static displacement control DIVERGES at the finest mesh (this is the disease, not a bug).
- **Fix / recipe (proven, Zone-A, no kernel change):** square `ASDShellQ4` with all transverse + out-of-plane DOF locked (a clean 1-D X-chain), right edge prescribed in displacement control; read the bar force as `Nxx·width` at a *normal* (full-thickness) rightmost element (series bar ⇒ axial force uniform); integrate `∫F dδ` to `F→0` for the dissipated energy. NO `-lch` — each element latches its OWN `getCharacteristicLength()` (= edge length for `ASDShellQ4`), so refining the mesh is exactly what exercises the regularization. The mesh-DEPENDENCE negative control compares the two COARSER meshes (where reg-OFF still converges and shows the spurious ~0.5 energy halving); the reg-ON sequence runs all three meshes cleanly (the rescaled longer tail removes the snap-back). Result: reg-ON energy ratios 1.000/1.000/1.000 across a 4× refine vs reg-OFF fine/coarse ≈ 0.50. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3b-struct.

### The generic LogStrainNDMaterial wrapper is UNSOUND for a stiffness-degrading (damage) inner — lift damage materials with a native FiniteStrainNDMaterial subclass
- **Bites:** lifting a small-strain plastic-DAMAGE material (LadrunoRCConcrete) to finite strain. The obvious route — `nDMaterial LogStrain $ftag $rcTag` — looks free (it worked for LadrunoJ2). It is NOT correct for a damage material. `LogStrainNDMaterial::setTrialF` commits the elastic left Cauchy–Green by recovering the elastic strain through the inner's INITIAL tangent: `εᵉ_{n+1} = C₀⁻¹:τ` (LogStrainNDMaterial.cpp ~195, comment literally "v1 assumes a linear-elastic inner law", `τ = Dᵉ:εᵉ`). For a damage inner `τ = (1−d)·Dᵉ:εᵉ`, so the recovered `εᵉ` and the committed `bᵉ` are shrunk by `(1−d)` ⇒ the trial elastic state drifts as damage grows (negligible at small strain, corrupts the finite-strain trajectory once `d` is appreciable). Silent — the wrapper still runs and looks plausible.
- **Fix (proven):** a native `FiniteStrainNDMaterial` subclass (LadrunoRCFiniteStrain, classTag 33018). Crucially, the RC plastic-damage spine carries NO tensorial plastic strain (its inelasticity is the scalar damage thresholds xt/xc + the effective-stress recompose; returnMap3D is a pure function of the TOTAL strain), so the multiplicative split is trivial: the elastic left Cauchy–Green is just the TOTAL `B = F Fᵀ`, recomputed from the element's F each step — there is NO `bᵉ` to track or mis-commit, and no `(1−d)` recovery. Reuse `LogStrainKernel.h` only for the kinematics (`hencky_voigt`, `assemble_material`, `spatial_tangent_full`). VALIDATION that catches the seam exactly: the finite GP Cauchy must equal `(small-strain material at ε=½ln(F Fᵀ)) / J` — a same-binary cross-check needing no numpy oracle (`tests/test_ladrunoRCFiniteStrain.py::test_finite_rc_patch_matches_smallstrain_oracle`). The RC secant tangent makes a strict K==FD test fail in the damage range, so gate the spatial-tangent seam in the ELASTIC regime only. Build gotcha: a fully-prescribed (all-DOF) finite probe SEGFAULTS under `constraints Transformation` (0 free DOFs) — use `Lagrange` (or `Penalty`); reconfirms the existing RC quirk. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 4b.
### Condensed-frame inner Newton: return BEST-EFFORT on maxIter, do not abort the global step
- **Bites:** an element with an internal iteration condensed before `CrdTransf` (the LadrunoDispBeamColumn3d embedded hinge `solveHingeJump`/`solveHingeJumpBiaxial`, and the coupled `LadrunoCohesiveHingeBiaxial`) is tempting to make "correct" by returning a non-zero (failure) code from `update()` when the inner Newton hits `maxIter` — so the global solver cuts the step instead of condensing a stiffness inconsistent with the unconverged internal state (ADR 33 §invariant literally says "cut the step"). **Returning `-1` makes softening analyses FRAGILE:** the inner Newton transiently misses at the cohesive activation kink (the coupled law's elliptical onset `r→1`) for trial `v` the GLOBAL Newton probes mid-iteration; aborting the whole step there fails analyses the global Newton would have recovered (a fixed-increment DisplacementControl/LoadControl run with no substepping just dies). Caught 4 previously-green coupled tests at load factor ~0.03 (the onset step).
- **Fix / convention:** return **0 with the last iterate + a loud `opserr` warning** — the standard OpenSees internal-iteration convention (`ForceBeamColumn` et al. do exactly this). The global `NormDispIncr` test is the real arbiter: an inconsistent internal state keeps the global residual up, so a truly-bad step cannot be silently accepted, and the CONVERGED global state always carries a converged inner state. Truly honoring the strict "cut the step" intent needs a robust inner solve that does not miss at the onset (a consistent off-radial / line-searched tangent) — deferred. Learned 2026-06-18 (adversarial review of [[33_ladruno_dispbeamcolumn3d_hinge_adr]] PR-3b + [[34_ladruno_cohesive_hinge_biaxial_adr]] PR-4a).

### Eigenvalue-floored 2x2 condensation: scale the floor by the DIAGONAL only, not the off-diagonal
- **Bites:** the biaxial hinge condenses a coupled 2x2 `K_aa` via an eigenvalue-floored inverse (`ladrunoFlooredInv2x2`), with the floor magnitude `1e-8*(|bulk_zz|+|bulk_yy|+|Czz|+|Cyy|)`. Folding the off-diagonal coupling `|Kaa_zy|` into that sum (to "scale with the full matrix") **over-floors** the modes near the elliptical onset — the floored inverse is then too damped, the radial inner Newton oscillates and hits `maxIter`, and 4 coupled tests break. The diagonal-only scale is the validated choice; the off-diagonal already enters the eigenvalues themselves, so it must NOT also inflate the floor. Learned 2026-06-18 (same review).

### CDPM2 unilateral crack-closure is AUTOMATIC if you re-split the converged effective stress every step
- **Bites:** the LadrunoConcrete3D damage nominal stress is `σ = (1−ω_t)σ̄_t + (1−ω_c)σ̄_c` (Eq.1) with `σ̄_t`/`σ̄_c` the Macaulay (positive/negative-eigenvalue) parts of the EFFECTIVE stress. The tempting "optimization" — split once at first cracking and reuse the FROZEN tension/compression projectors — looks fine on monotonic loading but **silently kills unilateral recovery**: when a previously-tensile principal direction goes into compression (crack closing), the frozen projector keeps multiplying that stress by `(1−ω_t)` (≈0 near peak), so a cracked-then-closed material carries almost no compressive stress — the cyclic value proposition dies exactly where it is sold (ADR §4.3 BLOCKING; the same trap as ASDConcrete3D's IMPL-EX frozen `PT`).
- **Fix (proven, zero extra state):** recompute the spectral split from the **converged effective stress every step**. Then a principal that flips negative is routed into the `(1−ω_c)` channel and is no longer touched by `ω_t` — full stiffness recovers automatically, with NO `s_rec·g_close` knob and NO stored projector. It is cheap (one symmetric 3×3 eig you already do for the return map) and tier-independent (freeze only the damage-driver extrapolation in IMPL-EX, never the recomposition split). Oracle gate `run_p2c_gate` DT3: a tension→compression reversal with `ω_t→1` recovers `nominal/effective → 1.000` in early compression. Partial recovery (debris-retained `s_rec<1`) is then a small additive refinement on top, not the mechanism. NB the per-channel softening law must be driven by the EXTREME effective principal of its sign (`max⟨σ̄_i⟩₊` / `max⟨−σ̄_i⟩₊`) so it reduces to the validated uniaxial drivers; do not feed it `E·κ_d` (correct for tension, gives the wrong peak `f_t≠f_c` for compression). Learned 2026-06-18, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] P2c.

### Spectral/eigenprojection gradients need the Voigt double-contraction weight; per-component micro-FD does NOT
- **Bites:** building the LadrunoConcrete3D analytic damaged tangent (P2e), the chain `∂(scalar of σ̄)/∂ε = (∂scalar/∂σ̄) : C_eff` is computed two different-looking ways, and they need DIFFERENT shear handling. An **analytic tensor gradient** like `∂λ_max/∂σ̄ = E_max` (the eigenprojection `nₐ⊗nₐ`) is a true tensor; contracting it with the stress increment is a tensor double contraction `A:B = ΣᵢAᵢBᵢ + 2Σ_offdiag` ⇒ in Voigt `{00,11,22,01,12,02}` you must apply the weight `[1,1,1,2,2,2]` to the off-diagonals BEFORE `Cᵀ_eff @ g`. Omitting the ×2 looks fine on coaxial (no-shear) paths and silently gives a WRONG tangent the moment shear is present (a 1e-2 rel FD-mismatch that vanishes in uniaxial tests — easy to miss). By contrast a gradient obtained by **per-Voigt-component micro-FD** (perturb `σ̄[k]` directly, e.g. `∂ε̃/∂σ̄`) is ALREADY the per-component partial, so it chains as a plain `Cᵀ_eff @ g` with NO extra weight — mixing the two conventions is the trap.
- **Why:** perturbing the single Voigt slot `σ̄_01` sets BOTH tensor entries `σ̄[0,1]=σ̄[1,0]`, so the per-component partial of an eigenvalue is `∂λ/∂σ̄_01 = 2 n₀n₁ = 2·E_max[0,1]` — i.e. the per-component partial already equals `W·mat_to_voigt(E_max)`. So: use `W=[1,1,1,2,2,2]` when you START from an analytic tensor gradient (`E_max`, `E_min`, `Δε_p/‖Δε_p‖`); use NO weight when you START from a per-component micro-FD. Verified by FD: shear-state tangent went 1.3e-2 → 1.3e-10 the instant the eigenprojection terms got the weight. Also: the `σ̄_lat=0` Macaulay kink (uniaxial-STRESS compression) is a genuine non-differentiable point — the analytic tangent is a valid subgradient (agrees on the loaded axial component, the central diff crosses the kink on the ~zero-stress lateral directions); gate the analytic tangent at smooth states (eigenvalues bounded from 0). Learned 2026-06-19, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] P2e (the C++ analytic tangent port must replicate this exactly).

### Selective mass scaling: the cap denominator must be ELEMENT mass, and the injection MUST be restored on teardown
- **Bites:** CentralDifferenceSMS ([[36_ladruno_selective_mass_scaling_adr]]) inflates throttling elements' mass by ADDING `(s_e−1)·m` to their NODES (`Node::setMass`). Two non-obvious traps caught by the milestone-3 adversarial review (NEEDS-REWORK, 2026-06-19): (1) **the `-maxAddedMass` cap was DEAD** because the %added-mass denominator summed `Node::getMass()`, which is **zero on a standard element-`-rho` model** (element mass is not stored as nodal mass) → `frac≡0` → the cap never fires and verbose lies "0% added". The denominator MUST be the sum of element translational lumped mass (`Σ_elements mdiag_trans`), not nodal mass. (2) **The injected fictitious mass is a persistent mutation of `Domain` node state**; if not restored it corrupts every later stage on the same model (the normal gravity→modal→SMS→post workflow gets wrong frequencies/base shear with NO warning). Fix: keep an integrator-owned `appliedDomain` + per-node ΔM record, restore (subtract) at the START of every `domainChanged` (re-baseline, so re-entrant calls don't compound) AND in the destructor (so switching integrators mid-session leaves baseline masses).
- **Why:** these are silent — all 21 green Zone-A tests passed with the dead cap and no-restore. Also from the same review: (a) self-reporting elements (bipenalty/kinematic couplings, `getExplicitCriticalTimeStep()>0`) carry a MASS-INDEPENDENT bound, so adding mass can't help — SKIP them (don't report "scaled"); (b) size against the per-element step from `elementCriticalDt` (self-report aware), and STAGE per-element injection locally, committing only if the node-by-node DOF walk maps (`pos==n`) so a non-node-major element can't leave a half-injected, un-rolled-back state; (c) `-cflAbort`/`-recompute` are REJECTED with SMS because their inherited path re-runs the element-mass eigensolve (which can't see the nodal augmentation) and would spuriously hard-abort a stable run. v1 is sequential / partition-interior only (parallel shared-node ΔM is not MPI-reduced) and does NOT yet exclude constrained nodes — both are one-time `domainChanged` warnings.

### The CDPM2 damage ω-solve needs a physical stress FLOOR, or a numerical residual flips it 0↔1
- **Bites:** in a damaged dual-projector model the per-channel softening ω is solved from the EXTREME effective principal of that sign (`sig_t_drive = max⟨σ̄_i⟩₊`, `sig_c_drive = max⟨−σ̄_i⟩₊`). In a uniaxial-STRESS state the lateral effective stress is only driven to ~0 by an inner lateral Newton, leaving a residual of order the Newton tolerance (~1e-10 MPa). If you gate the ω-solve on `sig_t_drive > 0` (any positive), that residual's SIGN decides whether `ω_t` activates: a +1e-10 MPa lateral "tension" + an already-accumulated tensile history (CDPM2 Eq.43 grows `κ_dt` even in compression) makes `f·exp(−κ_dt1/ε_f)→0` so `ω_t→1` in PURE COMPRESSION. It is mask-hidden (the compressive axial principal routes through the `(1−ω_c)` channel so `σ11` is unaffected) but it poisons the tensile damage state and is sign-of-noise fragile. **Fix:** floor the drive — solve ω only when `sig_t_drive > 1e-6·ft` (mirror `ω_c`/`fc`); a 1e-10 MPa residual is not a tensile state. Apply at EVERY site that recomputes ω (the stress update AND the analytic-tangent recompute, or the tangent decouples from the update). Gate `max(ω_t)≈0` over pure monotonic compression (and mirror). Caught by a 4-dim adversarial review, 2026-06-19, LadrunoConcrete3D P2e review-fix. **Distinct, deferred** issue (DT5 diagnostic, P2f): even with the floor, a REAL compression-then-tension reload still pre-damages tension to ~0 because literal CDPM2 Eq.43 (`κ̇_dt=ε̃̇`, no `(1−α_c)`) accumulates `κ_dt` in compression — that is the cyclic T/C-coupling temper (`β_c` Eq.50 + the `α_t`-weighting question), not the floor.

### Verify a ported analytic tangent against a NUMERICAL diff of the same-binary stress, not a cross-platform fixture of the oracle analytic tangent
- **Bites:** porting the LadrunoConcrete3D analytic DAMAGED tangent (P3b) to the C++ kernel, the obvious gate is "dump the oracle `damaged_tangent_analytic` 6x6 to the fixture and diff the C++ analytic tangent against it." But the oracle analytic tangent contains isolated micro-FD scalar gradients (`_dscalar_dsig`, h=1e-6) AND its `C_eff` is the oracle numerical consistent tangent (central diff, rel_step 1e-6) — both AMPLIFY last-ULP libm differences ~1e6x across platforms (g++ libstdc++ vs CPython/numpy), so a cross-platform numeric diff is noisy at ~1e-5..1e-3 and forces a weak tolerance (the same reason the committed-fixture rot-guard uses 1e-3).
- **Fix (proven, P3b):** make the tangent gate SELF-CONTAINED + same-binary — compute the C++ analytic damaged tangent (`returnMap` doTangent) AND a NUMERICAL central difference of the SAME C++ damaged nominal stress (the operator the global Newton consumes), assert analytic==numerical (rel ~1e-7, no cross-platform term) — exactly how the oracle `run_p2e_gate` checks analytic-vs-numerical. RIGOROUS because the C++ damage STRESS is already pinned to the oracle by the B3 fixture (~1e-14): a numerical tangent of an oracle-pinned stress certifies the analytic formula with zero cross-platform FD noise. The cross-platform oracle-analytic diff is a loose belt-and-suspenders check only (~1e-10..2e-7 here), never the gate floor. General rule for porting ANY analytic tangent: gate it against a numerical diff of the same-binary stress; pin the stress to the spec separately. Learned 2026-06-19, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] P3b.

### Uniaxial single-element tests leave the shear/off-diagonal tangent conversions at value 0 — probe the nDMaterial directly with `NDTest`
- **Bites:** an `nDMaterial` wrapper over a tensor-convention kernel must convert engineering↔tensor shear at the OpenSees boundary (for LadrunoConcrete3D: strain ×0.5 in / ×2 out; stress unscaled; tangent shear **COLUMNS** ×0.5 because the kernel tangent is `dσ_ij=2G dε_ij`). A single-`stdBrick` battery driven by `DisplacementControl` on **one** node-dof produces a pure uniaxial-stress state where **every engineering shear strain is identically zero**, so all three shear conversions are exercised only at the value 0 (`0.5·0 == 2.0·0 == 0`). A wrong factor (dropped `0.5`, or `2.0` instead of `0.5`) then leaves **every gate green** — the axial asserts never touch the shear block, and the global Newton merely converges a bit slower. The companion numpy-oracle test bypasses the wrapper entirely, and the g++ byte-check calls the kernel in TENSOR convention — so neither catches a wrapper-side shear-factor bug. (Found by the LadrunoConcrete3D adversarial review, 2026-06-19: MAJOR false-green coverage gap; the convention itself was verified mathematically CORRECT — this was a *test* gap, not a numeric bug.)
- **Fix (proven):** probe the material OBJECT directly with arbitrary 6-strain vectors **including shear** via the `NDTest` facility — `ops.NDTest("SetStrain", tag, e0..e5)` → `setTrialStrain`, `ops.NDTest("GetStress", tag)` / `"GetStrain"` / `"GetTangentStiffness", tag)` (36, row-major) → `getStress`/`getStrain`/`getTangent`. No element/solve/constraints in the way. Assert against the numpy oracle: elastic simple shear `σ_xy = G·γ_xy` (others 0); a multiaxial+shear strain's full 6-stress vs `elastic_C @ ε_tensor`; the 6×6 tangent vs `elastic_C` with shear **columns** halved (`C0[i][j]` for `j<3`, `0.5·C0[i][j]` for `j≥3`). The material is created with a bare `ops.model` + `ops.nDMaterial` (no element needed — `OPS_getNDMaterial(tag)` looks it up). General rule for ANY new nD wrapper: don't trust uniaxial element gates to verify the strain/stress/tangent conventions — add a direct `NDTest` shear/multiaxial leg vs the oracle. Learned 2026-06-19, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] post-ship adversarial review.
- **BUT `NDTest` is 3D-ONLY — it cannot reach the reduced views.** `OPS_NDSetStrain` (`SRC/interpreter/OpenSeesNDTestCommands.cpp`) always reads a fixed **6-component** strain and calls `setTrialStrain` on the **parser-built** material — which for a `dim`-mode wrapper is the 3D prototype. A `PlaneStrain`/`PlaneStress`/`AxiSymmetric`/`PlateFiber` view only exists once an ELEMENT calls `getCopy(type)`, so the reduced views must be gated through an element: `quad ... PlaneStrain|PlaneStress` (`FourNodeQuad`), `bbarQuad` (`ConstantPressureVolumeQuad`, requests `getCopy("AxiSymmetric2D")`), or a shell+fiber section (PlateFiber). The quad forwards `eleResponse(ele,"material",gp,"tangent")` to the material's `getTangent`, so you can read the wrapper's **condensed** reduced tangent for a rigorous check. For a peer-less material (CDPM2 has no upstream OpenSees equivalent) verify reduced views SELF-REFERENTIALLY against the shipped 3D material under the matching out-of-plane constraint: PlaneStrain in-plane stress == the 3D `NDTest` ε_22=0 slice; PlaneStress in-plane stress == a test-side ε_22-Newton on the 3D material driving σ_22→0 (elastic + a nonlinear damaging strain-path REPLAY with `CommitState` each step); pin both elastic tangents to the CLOSED-FORM plane-strain/plane-stress isotropic moduli. Learned 2026-06-19, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] Phase-2 reduced views (#299).

### Static condensation (σ_22=0) COMMUTES with the tensor-tangent shear-column halving — condense first, then halve
- **Bites:** a `dim`-mode `nDMaterial` whose kernel returns the tangent in the TENSOR convention (`dσ_ij=2G dε_ij`, e.g. LadrunoConcrete3D) needs BOTH (a) static condensation of the out-of-plane `33` dof for PlaneStress/PlateFiber (`D'[I][J]=D[I][J]−D[I][2]D[2][J]/D[2][2]`) and (b) halving the shear COLUMNS to hand the element `dσ/d(engineering strain)`. Order matters if you get it wrong. It is SAFE because the out-of-plane normal is index 2 (its column is unscaled, `s_2=1`), so scaling column `J` by `s_J` and then condensing gives `s_J·D'[I][J]` — identical to condensing then scaling. So: run `condenseTangent()` on the raw tensor `Dtan6` in `setTrialStrain` (the same convention the σ_22=0 nested Newton uses, where the `d22=Dtan6[2][2]` pivot is also a normal-component derivative, unscaled), then let `getTangent` halve the shear columns of the already-condensed matrix. The σ_22 condensation Newton and the rank-1 update both live in tensor convention; only the final element-facing map applies the ×0.5. (Verified: PlaneStress elastic tangent == the closed-form `E/(1−ν²)·[[1,ν,0],[ν,1,0],[0,0,(1−ν)/2]]`.) Learned 2026-06-19, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] reduced views (#299).

### Dual-damage IMPL-EX: the degraded-elastic secant `D_dam(ω):C0` is SPD only on SINGLE-SIGN principal states — and the time-ratio MUST be clamped
- **Bites:** a Tier-2 IMPL-EX scheme over a DUAL (tension `ω_t` / compression `ω_c`) damage model (LadrunoConcrete3D). The textbook IMPL-EX promise is "freeze the internal variables ⇒ the algorithmic tangent is a constant degraded-elastic secant `(1−ω̃)·C0`, symmetric SPD." That is true for SCALAR damage. For DUAL damage the secant is `D_dam(ω̃):C0` where `D_dam` is the spectral derivative of the per-principal damaged stress — and `D_dam` is symmetric but `D_dam @ C0` does **NOT commute**, so its symmetric part goes **INDEFINITE** on any **MIXED-SIGN, high-ω direction-contrast** state (a tensile-damaged principal `ω_t>~0.97` beside an undamaged/compressive one — utterly routine: a tension crack carrying lateral compression). Measured `λ_min(sym)≈−5e2` for `σ̄=[1,−2,−2]`, `ω_t=0.99`. So **"unconditional SPD across the snap-back" is FALSE** — it holds only where all principals share a sign (all-tensile / all-compressive), where `D_dam` collapses to one positive scalar. Genuine dual-damage *consistency* (explicit→implicit as Δt→0) and unconditional SPD are **mutually exclusive**; don't claim both. A uniaxial-strain gate never sees it (all principals ≥0 ⇒ SPD ⇒ false green). The secant IS still the exact `d(σ_rep)/d(Δε)` (FD-verifiable) and far better-conditioned than Tier-1, and committed physics stays exact — so it's a real robustness win, just CONDITIONAL. **Second bite:** the extrapolation time-ratio `r=Δt/Δt_n` must be **clamped to `[0, r_max≈2]`**. Unclamped, an adaptive step-GROWTH (small step then a big one — what a step-cutting solver does) drives `r` large ⇒ the linear extrapolation of the *bounded* `ω̃` overshoots past `[0,1)` and clamps to ~1 ⇒ the reported stress **collapses to ~0** with a near-singular/indefinite tangent; and `Δε̃ᵖ=r·Δεᵖ_n` injects an **unbounded spurious plastic strain** (huge sign-flipped stress). A negative `Δt` (only `Δt_n>0` was guarded) runs damage **backward** (spurious healing) — floor `r` at 0. ASDConcrete3D guards the same ratio (there via error-control step reduction). **Gate the conditional truth, not the wish:** assert SPD on single-sign; FD-verify the secant == `d(σ_rep)/d(Δε)`; PIN the mixed-sign indefiniteness (assert `λ_min<0` there, like the P2e Macaulay-kink "valid subgradient" gate); and exercise NON-uniform/negative Δt. Also: measure IMPL-EX `O(Δt)` convergence on the SMOOTH softening region (fixed-strain or L2 norm, fitted order across ≥3 levels) — the GLOBAL-max overstress is dominated by an irreducible one-step lag at the damage-ONSET C0 kink that does NOT refine away, so a global-max single-ratio-at-one-N gate is brittle/misleading. Learned 2026-06-20, [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]] P3 IMPL-EX adversarial review (#301 oracle → #304 fixes; 6 findings, 0 refuted). Carries to the C++ Tier-2 port.
- **C++ PORT (#309) caveat — IMPL-EX wants ~UNIFORM pseudo-time; `-implex` + `DisplacementControl` past a limit point DIVERGES.** The wrapper reads the extrapolation time `dt = ops_Dt` (the ASDConcrete3D pattern). Under `LoadControl`/transient the pseudo-time is monotone-uniform ⇒ `r=dt/dt_n≈1` and IMPL-EX is well-behaved (verified: `-implex` LoadControl compression converges 300/300, matches Tier-1 to ~4e-11). Under `DisplacementControl`, `ops_Dt` is the LOAD-FACTOR (λ) increment, which is NON-uniform/non-monotone near a limit point — so `r` spikes and, even clamped to 2, a 2× extrapolation of a steep damage-onset `dwt` saturates `ω̃→1` ⇒ near-singular secant ⇒ the global Newton diverges. Tier-2's headline use is softening, which needs DisplacementControl/arc-length — so this is a real usability boundary, NOT just a test artifact: for **softening + DisplacementControl use Tier-1 + an unsymmetric solver** (or Tier-3 explicit / transient). `-implex` is for implicit TRANSIENT dynamics + uniform quasi-static. Element tests therefore use LoadControl (robust) + NDTest (strain-controlled, dodges the global limit point: drive tension softening strain-by-strain, the axial `C[0][0]` is <0 for Tier-1 but >0 for `-implex` — the SPD payoff, convention-robust on the normal diagonal). Future improvement: derive the IMPL-EX `dt` from a monotone control parameter instead of λ. Learned 2026-06-20 (#309).

### A time-varying prescribed SP (support motion) DOES work under explicit (CentralDifference / ExplicitBathe) — it is NOT "an implicit/static construct"
- **Bites:** the Tier-3 explicit demo for LadrunoConcrete3D (#328) concluded that "a ramped prescribed-SP is an implicit/static construct that fails under CentralDifference" and fell back to **free dynamics** (initial face velocity), which can only show "runs through cracking, stays bounded" — NOT a clean quasi-static peak+softening backbone (the elastic strain-to-peak `ε0=ft/E~1e-4` is a blink under free dynamics). That conclusion is **wrong**. A time-varying single-point constraint imposed through the `Transformation` handler is exactly **support-motion input** and integrates fine under both fork explicit steppers.
- **Recipe (proven, #333):** `ops.sp(node, dof, rate)` referenced by a `pattern Plain` + `timeSeries Linear` (so the imposed displacement is `rate·t`) + `constraints Transformation` + `system Diagonal` + `algorithm Linear` + element mass via the material's `-rho`. Step with `dt ≈ 0.3–0.4 · le/c_d` (`c_d=sqrt((K+4G/3)/rho)`). Keep loading SLOW — hundreds of wave transits over the run — so the gauss stress tracks the quasi-static backbone, and add a small **mass-proportional** `rayleigh(alphaM,0,0,0)` (stays tangent-free) to settle the lateral (Poisson) ringing the cracking onset excites. The result: the nominal stress peaks at `~ft` then softens with `ω_t→~1`, tracking the numpy oracle backbone — under BOTH `CentralDifferenceLadruno` and `ExplicitBathe`, with NO unsymmetric solver and NO step-cutting.
- **Why the original confusion:** free dynamics + a heavy element makes the driven face translate ballistically (no restoring control over the rate), so it never produces a controlled strain ramp; that is a property of *free* dynamics, not of prescribed SP under explicit. Verified 2026-06-20 (`tests/test_ladrunoConcrete3D_element.py` "Tier-3 EXPLICIT-DYNAMICS" section, #333; full file 30/30 on a fresh OpenSeesPy build). [[31_ladruno_concrete3d_adr|LadrunoConcrete3D]].

### `system Diagonal` silently drops the off-diagonal of a condensed mass; `FullGeneral`/`BandGeneral` swallow a singular mass and report success
- **Bites:** the explicit recipe is `Transformation` (or any handler that builds a condensed `TᵀMT`) + `system Diagonal`. When a constraint's condensed mass has genuine off-diagonal coupling — e.g. a `rigidLink -beam`/`rigidDiaphragm` slave whose offset translational mass produces a translation↔rotation transport term `m·d` on the retained master — `DiagonalSOE::addA` (default `lumpDiagonal=false`) keeps only `diag(A)` and **silently discards the coupling**, with NO warning. The run completes cleanly (all steps, all-finite) but the inertia is wrong: the coupling-induced response is *identically zero* and the coupled mode is detuned (a verified case: rz period 1.54 s → 1.40 s, the induced `uy` 0.075 → 0.0). This is why `Transformation`-in-explicit is not a substitute for ADR-30's projection handler.
- **Second trap (same family):** a free DOF with **zero lumped mass** makes the assembled `M` singular, and the SOE layer's response is *inconsistent and untrustworthy* — `system Diagonal` aborts at solve time with an opaque `DiagonalDirectSolver aii = 0` (analyze → −2), but `system FullGeneral` **and** `system BandGeneral` let the LAPACK factorization fail (`matrix singular U(i,i)=0`) and then **return success (rc = 0) with a non-physical result** — silently wrong. Worse, the garbage can numerically coincide with a well-posed control answer, so "differs from a good run" is not a safe detector. Lesson for any new explicit handler: detect massless DOFs up front (handle()-time scan) with a named, actionable error; never rely on the SOE to catch them.
- **Why it stays hidden:** both behaviors pass every smoke test that only checks "did the analysis run and stay finite". Empirically pinned 2026-06-19 (ADR-30 Phase 0 / Gate-0, build 8c8edb2); tests `tests/test_adr30_projection_p0.py`. Mechanism confirmed against `TransformationFE::getTangent` (forms dense `TᵀMT` with off-diagonals) → `DiagonalSOE::addA` (reads only `m(i,i)`).

### The biaxial cohesive damage D(r) is mix-INDEPENDENT at the linear mix — the "frozen-mix tangent" was already exact
- **Bites:** ADR 34 / the handoff flagged `LadrunoCohesiveHingeBiaxial`'s 2x2 tangent as "frozen-mix, exact only on radial paths" and listed a "consistent off-radial tangent (carry ∂mode-mix/∂α)" as the HIGHEST-VALUE next item. Implementing that term and FD-checking it reveals it is **identically zero** under the default LINEAR mix: the exact-reduction calibration `Kpen_i = ratio·Mc_i²/(2 Gf_i)` forces `c_i = Mc_i²/Kpen_i = 2 Gf_i/ratio`, so `S = (2/ratio)·Gf_lin` AND `Esoft = Gf_lin·(1−1/ratio)` are BOTH proportional to the linear `Gf_lin` → the softening rate `A = S/Esoft = 2/(ratio−1)` is a **constant independent of the mode mix** → damage `D = 1 − T_eff(r)/(S r)` depends on `r` ALONE. Verified against the compiled material: at fixed `r`, damage spread across pure-z / 45° / pure-y is `0` (EXP) / `1e-16` (LINEAR). So the v1 tangent was ALREADY the exact consistent off-radial Jacobian; do NOT spend a session "fixing" a non-existent frozen-mix error.
- **Fix / consequence (PR-4b):** the `∂mode-mix/∂α` term only becomes nonzero with a NONLINEAR `Gf_mix` — implemented as Benzeggagh-Kenane `Gf_mix = Gf_z + (Gf_y−Gf_z) w_y^η` (`-bk η`, default 1 = linear, bit-identical → the 88 stay green). With `η≠1`, `S` is no longer ∝ `Gf_mix`, `D` is mix-dependent, and the consistent mode-mix tangent term is LIVE (FD-gated in `tests/test_ladrunoCohesiveHingeBiaxial_material.py`). SECOND subtlety: even with B-K, at the near-rigid DEFAULT `ratio=1000` the term is tiny (~1e-4 of the tangent) because `D ≈ 1 − 1/r` (mix-independent secant geometry) dominates and the mix only enters via `A = O(1/ratio)`; it is numerically significant only for softer penalties (~1e-2 at `ratio=10`, where the FD discrimination test runs). THIRD: the real weak-axis-dominant indirect-control robustness gap is the non-smooth KINK at the elliptical onset `r=1` (a SOLVER / inner-line-search matter), NOT tangent completeness — keep the best-effort inner Newton + prescribe-rotations workaround. Learned 2026-06-19, [[34_ladruno_cohesive_hinge_biaxial_adr]] PR-4b.

### openseespy `opserr` messages with a literal `%` are silently mangled — `PySys_FormatStderr` treated the message as a format string
- **Bites:** any `opserr << "... % ..."` text printed under **openseespy** lost everything from the `%` onward (or garbled it): `PythonStream::err_out` (`SRC/interpreter/PythonStream.h`) passed the already-formatted message AS the format arg — `PySys_FormatStderr(msg.c_str())` — so a literal `%` became a bogus printf conversion. The Tcl `StandardStream` path (`cerr << s`) was UNAFFECTED, so the bug is invisible from `OpenSees.exe`. Surfaced validating the SMS `-maxAddedMass` cap warning: `"added mass 160.68% of model mass"` printed as `"added mass 160.68"` and `"% exceeds -maxAddedMass cap 5%"` vanished entirely — making the SMS-CAP-DEAD fix illegible to Python users (the very audience).
- **Fix (proven):** `PySys_FormatStderr("%s", msg.c_str())` — pass the message as a `%s` ARGUMENT, never as the format. One-line, strictly more correct, fixes EVERY `%`-containing `opserr` message across openseespy (not just SMS). Marked `// Ladruno`, ledgered as an upstreamable bugfix (CWE-134 uncontrolled format string; the line was the pristine upstream form — worth a PR to OpenSeesFramework). `err_out` is called PER-TOKEN (each `operator<<` fragment), so a multi-part `opserr << ... << x << ...` is many separate calls — the old bug was per-fragment (only `%`-bearing fragments mangled), matching the symptom. **Residual ceiling (NOT fixed by this — it's CPython):** `PySys_FormatStderr` still caps output at ~1000 bytes internally, so a SINGLE string-literal fragment longer than ~1000 bytes can truncate; the `%s` form is strictly safer on length but not unlimited (don't misattribute a long-message truncation to this line). General rule: when a test asserts on `opserr` warning TEXT under openseespy, remember pre-fix binaries drop `%`; and never feed user/contentful strings to a `*printf`-family format parameter. Learned 2026-06-20, mass-scaling validation (T-CAP), [[37_ladruno_mass_scaling_validation_plan]].

### Selective mass scaling EXCLUDES elements at MP-constrained slave nodes (they still govern) — but NOT SP/fixed nodes
- **Bites:** `CentralDifferenceSMS` sizes each sub-target element against its BARE element pencil dt_e and injects `(s-1)m` onto that element's nodes. If a node is the CONSTRAINED (slave) side of an MP_Constraint (`equalDOF`/`rigidDiaphragm`/`rigidLink`/generic MP), the handler eliminates that DOF and redistributes its mass to the retained node via `T^T M T` — so the injected fictitious mass does NOT land where the element pencil assumed, the dt boost silently fails to materialize, and mass is mis-distributed. (NB RBE2/RBE3 are Ladruno ELEMENTS, not MP_Constraints — they don't appear in `getMPs()`; a bipenalty coupling is handled by the SEPARATE self-report skip, [[project_ladruno_mass_scaling]].) v1 just warned; **v1.1 (ADR-36) EXCLUDES** any sub-target element touching an MP slave node: skip it, count it (`nConstrained`), and report it STILL GOVERNS at its un-scaled `dt_e` (so the user honestly learns dtTarget is not delivered for it — lower dt or remove the constraint). The exclusion is keyed on `Domain::getMPs()->getNodeConstrained()` only.
- **Do NOT exclude on SP/fixed nodes:** mass on a fully-removed DOF is inert (the handler drops it), and the MOTIVATING case (a refinement zone at a fixed support / pile tip) puts tiny elements right on SP-fixed nodes — excluding those would defeat SMS exactly where it's needed. A partially-fixed node still scales its FREE DOFs fine. So only MP *slave* nodes trigger exclusion, never SP.
- **Consequence for the user:** excluding a constrained element means it keeps its tiny dt_e and GOVERNS, so running at dtTarget would make it unstable — the warning says so (`still GOVERN at dt_e=... < dtTarget`). This is the honest trade vs v1's silent mass-corruption: SMS refuses to fake a stable big step it cannot deliver through a constraint. Gated by `tests/test_massScaling_validation.py::test_constrained_element_excluded_free_scaled` (constrained truss excluded + slave stays massless, free truss still scaled `1/2`). Learned 2026-06-20, mass-scaling v1.1 constraint-exclusion guard, [[36_ladruno_selective_mass_scaling_adr]] / [[project_ladruno_mass_scaling]].
### LadrunoProjection (ADR-30): every TIED DOF needs lumped mass, and the transport Ccr is frozen small-rotation
- **Bites (massless tied DOF):** unlike `constraints Transformation` (which ELIMINATES slave DOFs), `LadrunoProjection` keeps every constrained/retained DOF in the equation set and reads its diagonal mass for the projection weight `LᵀML`. So a tied DOF with zero lumped mass is refused at the first solve (`buildMass`, named error). This bites `rigidLink -beam` and 3D `rigidDiaphragm`, which tie the slave's **perpendicular rotation** (rz) — physical models usually give those nodes zero rotational mass. Fix: add a small rotational mass (~0.01–0.1% of the node's translational mass) to every tied rotational DOF, OR use Penalty/Transformation. (Relaxing this = ADR-30 P4 SOE-cooperative elimination.) Learned 2026-06-20, ADR-30 P2 (Gate-C).
- **Bites (frozen Ccr):** the transport lever-arm coefficients (`Ccr`) are captured ONCE at constraint construction (small-rotation). For an explicit run that accumulates large rotation (|θ| ≳ 0.1 rad) the projection enforces the tie with stale geometry → silently wrong. This is the SAME limitation `constraints Transformation` carries (both read the same frozen `Ccr` from the `MP_Constraint`), so it is not a regression — but there is no runtime guard yet (a per-step lever-arm staleness check is deferred). For finite-rotation rigid offsets use the RBE2 element route (`LadrunoKinematicCoupling`) or Transformation. Learned 2026-06-20, ADR-30 P2 (Gate-C).
- **Minor (near-singular `LᵀML`):** the per-group rank check (`Matrix::Solve`/DGESV) catches an EXACT zero pivot (massless direction, exact redundancy) but not a merely ill-conditioned `LᵀML` (e.g. a tied DOF with mass ~1e-12). Acceptable for v1 (ADR O1); a condition-number gate is a deferred hardening.

### `wipe` / `Domain::clearAll()` did NOT clear EQ_Constraints (upstream bug) — leaks across models
- **Bites:** `EQ_Constraint` (the `equationConstraint` command, a later upstream addition) was never wired into `Domain::clearAll()` — it clears `theSPs/thePCs/theMPs` but omitted `theEQs`. So `ops.wipe()` (and any `clearAll`) LEAVES equation constraints in the domain; the next model silently inherits them. This is invisible with the stock handlers (none of them iterate `getEQs()` in the common path), so it sat latent. It surfaces the moment a handler reads `getEQs()` — `LadrunoProjectionHandler` (ADR-30 P3) does, and a stale EQ from a prior model then mis-assembles a constraint group (wrong groups / partition-guard refusal / wrong projection) in the NEXT analysis. Found 2026-06-20 by the ADR-30 P3 full-suite regression: the EQ test poisoned every subsequent projection test ONLY in combined runs (passed in isolation) — the classic test-ordering signature of leaked global domain state.
- **Fix:** add `theEQs->clearAll();` to `Domain::clearAll()` (one line, mirrors `theMPs->clearAll()`). Upstreamable. Ledgered in LEDGER_vanilla. **General lesson:** a "passes alone, fails in a combined pytest run" ordering failure = leaked global OpenSees state; suspect a container `wipe`/`clearAll` doesn't clear (here EQ_Constraints). Run new constraint/handler features in a COMBINED suite, not just in isolation, to catch it.
### Prescribed motion under a Plain-style handler is the HANDLER's `applyLoad()` job; a plain SP supplies disp ONLY (vel/accel=0); a prescribed constraint-MASTER would silently force slaves to 0
- **Bites (the silent one):** under `LadrunoProjectionHandler` (ADR-30 P4b) a prescribed-motion DOF (non-homogeneous `SP_Constraint` or `imposedMotion`) is SP-excluded (eqn=−1) just like a `fix`. If that DOF is also a constraint MASTER (e.g. `equalDOF`/`rigidLink` retained node), the master's accel is taken as 0 (column dropped → slaves routed to the "SP-fixed-master" fixed set → slaves' accel forced to **0**) — a SILENTLY WRONG answer (the slaves should follow `C·a_prescribed`). P4b refuses this with a named "PRESCRIBED MASTER" error rather than run wrong; driving slaves from a prescribed master (the literal "overwrite a before projecting" known-RHS projection) is deferred to P4c (Tier 2). A prescribed DOF that is also a constraint SLAVE = overconstraint → "PRESCRIBED SLAVE" refusal.
- **Mechanism / how to enforce it:** the prescribed DISPLACEMENT is NOT set by the integrator (its `setResponse` scatter skips eqn<0 DOFs) and NOT by `ImposedMotionSP::applyConstraint` (which sets node vel/accel only — "disp is the responsibility of the constraint handler"). It is the HANDLER's `applyLoad()` override: `node->setTrialDisp(sp->getValue()+sp->getInitialValue(), dof)` each step, called from `AnalysisModel::updateDomain(time,dt)` AFTER `Domain::applyLoad` (so an `ImposedMotionSP`'s node vel/accel already set survive). `DOF_Group::setNode*` skipping eqn<0 is exactly what lets the imposed value survive the integrator's full-vector scatter — the free DOFs then feel the support through `F_int`/damping, integrator UNCHANGED. **A plain (non-`imposedMotion`) `SP_Constraint` supplies only the displacement; its prescribed vel/accel stay 0** — the SAME limitation `constraints Transformation` carries (`TransformationDOF_Group::enforceSPs` also only `setTrialDisp(getValue())`), NOT a regression; use `imposedMotion`+`groundMotion` when prescribed vel/accel matter. Learned 2026-06-20, ADR-30 P4b, [[30_ladruno_explicit_constraint_projection_adr]] / [[project_explicit_constraint_projection]].
### A prescribed-motion MASTER must drive its slaves KINEMATICALLY, not by acceleration projection (the displacement tie won't hold otherwise)
- **Bites (the abandoned approach):** ADR §2.4 says drive a prescribed master's slaves by "overwriting a on those DOFs with the prescribed acceleration before projecting the rest" — the M-orthogonal known-RHS projection `a_c = C_f a_f + C_p a_p`. Acceleration-exact, but the DISPLACEMENT tie `u_c = C u_master` DRIFTS and **does not converge**: a prescribed master's disp is externally imposed (an `imposedMotion` GroundMotion integrates accel→disp at its OWN internal `dtInt`; a constant SP holds disp fixed with `a_p=0`), while the slave is leap-frog integrated at the analysis `dt`. Two different integrators across the tie → slave disp ≠ master disp. Measured rel ~1.5e-3, and the error GREW as `dt` moved away from `dtInt` (smallest when `dt==dtInt`). A constant prescribed master is worse: `a_p=0` ⇒ the slave never leaves its IC. Every OTHER constraint in the handler is machine-exact precisely because both tied DOFs leap-frog the SAME projected accel — a prescribed master breaks that co-integration.
- **Fix (P4c, shipped):** a slave driven PURELY by prescribed master(s) is itself fully determined ⇒ **KINEMATIC imposition** — exclude it from the equation set (`eqn=-1`, like its masters) and set `u_c=ΣC_k u_{m_k}+delta`, `v_c=ΣC_k v_{m_k}`, `a_c=ΣC_k a_{m_k}` directly each step in the handler's `applyLoad()` (AFTER the masters' own disp is imposed). Exact, zero drift, like Transformation's elimination. A slave tied to BOTH a free and a prescribed master (MIXED) is refused (its disp tie can't ride the free-DOF projection); a group left with zero free equations is refused (clean error, not a segfault — ADR §5.2). **General lesson:** in an explicit scheme a kinematic tie to an EXTERNALLY-imposed DOF must be imposed kinematically (set u/v/a), never enforced through an acceleration-only projection — the integrators won't match.
- **Companion quirk (`Transformation`+CDL, found validating P4c):** when a node is BOTH a `rigidLink` master AND carries an `imposedMotion` SP on one of its DOFs, `Transformation`+`CentralDifferenceLadruno` SILENTLY DROPS the rigid tie on the OTHER (free) DOFs of that master (a horizontal `rigidLink -bar` master with `imposedMotion` on `ux` → the slave's `uy` stops tracking). `LadrunoProjection` holds it correctly, so Transformation is NOT a trustworthy reference there — `tests/test_adr30_projection_p4c.py::TC3` self-verifies the exact ties instead. Learned 2026-06-20, ADR-30 P4c, [[30_ladruno_explicit_constraint_projection_adr]] / [[project_explicit_constraint_projection]].
### SMS sizes against the betaK-DAMPED explicit step (closed form s=T^2+2Tc), NOT the undamped 2/w — and excludes alphaM
- **Bites:** with stiffness-proportional (`betaK`) Rayleigh damping the central-difference stable step SHRINKS (`xi = betaK*w/2` grows with `w`), so SMS sizing against the UNDAMPED element step `dt_e = 2/w_max` UNDER-scales — the scaled element is still unstable at `dtTarget` and the run diverges despite SMS reporting success. v1.1 sizes against the damped step: with `c = betaK/dt_e` (= 0.5*betaK*w_max), the betaK-damped step at mass-scale `s` is `dt_d(s) = (2/w0)(sqrt(s+c^2) - c)`, which inverts in CLOSED FORM to **`s = T^2 + 2*T*c`, `T = dtTarget/dt_e`** (no bisection). Reduces to the undamped `T^2` when `betaK=0` (no-damping models byte-identical). The SKIP test also moved from undamped `dt_e>=dtTarget` to damped `dt_d(1)>=dtTarget` — an element whose undamped step clears `dtTarget` but whose damped step does not is now correctly scaled.
- **alphaM (mass-proportional) is intentionally EXCLUDED from the sizing:** it does not reduce the high-frequency step that governs explicit stability (`xi_alpha = alphaM/(2w)` DECREASES with `w`), and folding it in across scales is non-monotonic (`dt_d(s)` would asymptote and may never reach `dtTarget`). Only `betaK = getRayleighDampingFactors()(1)` enters — mirroring the betaK term of `computeCriticalTimeStep`'s damped estimate, so SMS sizes consistently with what `criticalTimeStep()` reports. Gated by `tests/test_massScaling_validation.py::test_betaK_damped_sizing` (betaK injects MORE mass than undamped, matching the closed-form ratio; a betaK-blind SMS ties them). Learned 2026-06-20, mass-scaling v1.1 betaK-damped sizing, [[36_ladruno_selective_mass_scaling_adr]] / [[project_ladruno_mass_scaling]].
### The explicit accel solve `a=M⁻¹r` is driven by the ALGORITHM, not the integrator — inject a non-diagonal mass via a post-solve hook + the FACTORED DiagonalSOE
- **Bites:** wiring CONSISTENT (Olovsson) mass scaling (ADR 38, `CentralDifferenceSMSConsistent`) means solving `M̃ a = r` with a NON-diagonal `M̃ = M_lump + ΣM_bar_e`. The instinct is "override the integrator's solve" — but in OpenSees the per-step acceleration solve is performed by the **`Linear` SolutionAlgorithm** (`theSOE->solve()` → `integrator->update(a)`), NOT inside the integrator. The integrator only forms M (`formEleTangent`/`addMtoTang`) and consumes the solved `a`. So a non-diagonal mass solve must be injected at the TWO sites that consume the diagonal result: the `newStep()` first-step starter (which DOES solve inline) and `update()`.
- **The clean trick (no stored state):** the leap-frog uses `system Diagonal`, whose `DiagonalDirectSolver` factors A IN PLACE — **post-solve `getDiagonalA()` holds `A[i] = 1/mass_i`**. That single array is BOTH the Jacobi preconditioner AND the way to recover the RHS: the incoming `a = M_lump⁻¹ r`, so `r = M_lump .* a = a(i)/A[i]`. A protected base no-op hook `refineAccel(Vector&)` (default no-op ⇒ the lumped path stays **byte-identical**, not something to prove) is overridden to recover `r`, run a matrix-free Jacobi-preconditioned CG (`consistentMatVec` gathers/scatters element DOFs via `FE_Element::getID()`), and replace `a` with `M̃⁻¹ r`. CG converges in **3–21 iters** because `M_bar` is a small perturbation of the dominant lumped diagonal. NO `Node::setMass` mutation — the cross-node coupling can't live in per-node mass, so there's also no inject/restore lifecycle. **ExplicitBathe variant (33010):** the Noh-Bathe scheme takes TWO diagonal solves per step (`A_tpdt` from the external `Linear` solve consumed in `update()`, `A_tdt` from an INLINE second solve inside `update()`); hook `refineAccel` at both. The inline second solve does NOT re-`formTangent`, but `DiagonalDirectSolver` only factors when `isAfactored==false` (the second `solve()` takes the "just solve" `X=A·B` branch), so `getDiagonalA()` is still `1/mass` at both sites — `r = a/Ainv` recovery is valid for both. The LUMPED ExplicitBatheSMS (33009) needs NO hook (ExplicitBathe assembles only mass on the RHS, so nodal injection is seen directly).
- **The clean trick (no stored state):** the leap-frog uses `system Diagonal`, whose `DiagonalDirectSolver` factors A IN PLACE — **post-solve `getDiagonalA()` holds `A[i] = 1/mass_i`**. That single array is BOTH the Jacobi preconditioner AND the way to recover the RHS: the incoming `a = M_lump⁻¹ r`, so `r = M_lump .* a = a(i)/A[i]`. A protected base no-op hook `refineAccel(Vector&)` (default no-op ⇒ the lumped path stays **byte-identical**, not something to prove) is overridden to recover `r`, run a matrix-free Jacobi-preconditioned CG (`consistentMatVec` gathers/scatters element DOFs via `FE_Element::getID()`), and replace `a` with `M̃⁻¹ r`. CG converges in **3–21 iters** because `M_bar` is a small perturbation of the dominant lumped diagonal. NO `Node::setMass` mutation — the cross-node coupling can't live in per-node mass, so there's also no inject/restore lifecycle.
- **General lessons:** (1) to change WHAT a leap-frog explicit integrator inverts without writing a new SOE/solver, hook the post-solve accel and reuse the factored `DiagonalSOE` as preconditioner; (2) consistent (Olovsson) scaling preserves f1 (`−0.17%` measured via transient FFT) where lumped craters it (`−53%`) at the same dtTarget, because `M_bar`'s row sums are zero (rigid translation gets no added inertia). Learned 2026-06-20, [[38_ladruno_consistent_mass_scaling_adr]] / [[project_ladruno_mass_scaling]].
- **Consistent SMS needs a STRICTER constraint exclusion than lumped (slave AND master MP nodes), or the matvec reads M_bar out of bounds.** The lumped SMS excludes only MP *slave* nodes (`getNodeConstrained`) — fine, because it injects into `Node` mass and the constraint handler then transforms it (`TᵀMT`). The CONSISTENT path is matrix-free in equation space via `FE_Element::getID()`, so it bypasses the handler. Under the **default `TransformationConstraintHandler`** an element attached to an MP-constrained node gets a TRANSFORMED FE id — different basis AND different size: a *retained/master* node's `TransformationDOF_Group` absorbs the slave's DOFs, so `getID().Size()` can EXCEED the element's own `n`. Pairing that oversized id with the original-basis `n×n` M_bar both scatters in the wrong basis and reads `Mbar(a,b)` out of bounds (no bounds check in release `Matrix::operator()`). Fix (ADR 38): exclude any scaled element touching EITHER a slave OR a master MP node (`getNodeConstrained` + `getNodeRetained`) AND a defensive `if (getID().Size()!=n) skip`. Side benefit: with no M_bar on any constrained equation, the ADR-30 projector's lumped-diagonal mass stays exact (no consistent/projector mass mismatch). Caught by a 45-min adversarial review (the OOB fires on the exact SSI/diaphragm models the feature targets); regression `tests/test_centralDifferenceSMSConsistent_integrator.py::test_consistent_excludes_constraint_touching_elements`. **General lesson: a matrix-free element-DOF kernel that uses `FE_Element::getID()` MUST assume the id can be transformed (resized + rebased) by the constraint handler — never assume `getID().Size()==getMass().noRows()`.** Learned 2026-06-20, [[38_ladruno_consistent_mass_scaling_adr]].
- **Consistent SMS energy: the `EnergyBalanceRecorder` KE did NOT include the M_bar term — FIXED by V4 (global registry conduit).** KE is `½vᵀMv` summed over element+nodal `getMass()` (`EnergyBalanceKernel.h`); the consistent scaling mass `M_bar` lives only in the integrator's matrix-free blocks, never in `Node`/element mass, so the recorder under-reported KE by `½vᵀ(ΣM_bar)v` and the balance did NOT close for the consistent integrators (it DOES for the lumped siblings, whose nodal injection is visible to the recorder). **The recorder holds only a `Domain*` — no integrator handle — so the fix is a process-global `Ladruno::MassScalingEnergyRegistry` (`LadrunoMassScalingEnergy.{h,cpp}`):** the active consistent integrator `publish()`es its per-element node-major `M̄ₑ` keyed by element tag (`ConsistentBlock` gained `eleTag`) at the end of `domainChanged` and `clear()`s on dtor/recvSelf (owner-guarded — a stale teardown can't wipe a newer publisher); the shared kernel queries it per element and adds `½vᵀM̄ₑv`. The `M̄ₑ` is stored node-major, the exact order `addElementEnergy` already gathers velocities, so NO equation-number bookkeeping. EMPTY for the lumped path + every base integrator (`active()==false`) ⇒ recorded KE byte-identical there, no double count. **General lesson: when a recorder needs live state the active integrator owns but the recorder can't reach (only `Domain*`), a process-global owner-guarded registry the integrator publishes to in `domainChanged` / clears in its dtor is the minimal seam (single active integrator per process; NOT parallel-shared-node aware — that's a separate concern).** Learned 2026-06-20, V4 (this PR), [[38_ladruno_consistent_mass_scaling_adr]].
- **Energy-balance closure test design: a `RES`/`ERR%` gate is the WRONG metric for the consistent path; use `KE+IE` conservation with a VELOCITY IC, or an instantaneous analytic-KE oracle.** Two traps when validating V4: (1) **the closure RES sign depends on the IC.** `RES = ULW − (KE+IE+DW)`; for free vibration (`ULW=DW=0`) `RES = −(KE+IE)`. A **displacement** IC (`KE₀=0`) gives `IE = SE(t)−SE₀` and `KE+IE → 0` so `RES→0` (closes); a **velocity** IC seeds `KE₀=E₀` as an *initial condition* the balance never accounts as work, so `RES = −E₀` PERMANENTLY (`ERR%≈100%`) — by design, not a bug. The existing lumped T-ENERGY test exploits this: velocity IC, assert `KE+IE` stays at the constant level `E₀` (drift), NOT `RES→0`. (2) **`ERR%` does NOT converge with `dt` for the consistent path** because consistent scaling pushes the *deformation* modes to MARGINAL stability (`dt_e→dtTarget`), and the recorder's instantaneous `½vᵀM̃v` + trapezoidal `IE` can't track central-difference's *modified* energy for a mode at the stability edge — `ERR%` actually grew (14%→30%) as `dt` shrank. The robust V4 test is the **instantaneous analytic oracle**: seed an ALTERNATING (deformation-rich) velocity field so `M̄` carries a large majority of the KE (≈77% here — `M̄`'s zero row sums mean a SMOOTH/rigid field loads it almost not at all), then per step assert recorder `KE == ½vᵀM_lump v + ½Σ_e β_e(m_a/2)(v_a−v_b)²` (truss closed form, `dt_e=L/c`, `β=(dtTarget/dt_e)²−1`) to ~machine precision. Dissipation- and stability-independent ⇒ works for the Noh-Bathe families too. Learned 2026-06-20, V4 `tests/test_massScaling_consistent_energy.py`, [[38_ladruno_consistent_mass_scaling_adr]].
### LUMPED SMS is parallel-correct WITHOUT an explicit ΔM reduction — the distributed/MPI diagonal solver already sums shared-node mass across ranks. CONSISTENT (Olovsson) is NOT (rank-local PCG). Tcl can't reach SMS.
- **The non-obvious truth (validated bit-identical):** the instinct is that lumped SMS injecting `(s−1)·m` via `Node::setMass` on a partition-boundary node only adds rank-local mass (so shared nodes "desync"). FALSE. A parallel build auto-swaps the explicit `system Diagonal` → **`DistributedDiagonalSOE`** (OpenSeesSP, `SRC/tcl/commands.cpp:3213` under `_PARALLEL_PROCESSING`) and `system MPIDiagonal` → **`MPIDiagonalSOE`** (OpenSeesMP). BOTH solvers **sum the shared-DOF diagonal across ranks at solve time** (`DistributedDiagonalSolver`: gather→P0→`*vectShared+=otherShared`→broadcast; `MPIDiagonalSolver::intersectionsAB`: accumulates `A` on the first solve, reduces only `B` thereafter). `CentralDifferenceLadruno` reads M through `formTangent→theLinSOE->solve()` — i.e. the REDUCED diagonal, not raw `Node::getMass()`. Each element lives wholly on one rank ⇒ per-element `dt_e`/`s` are correct rank-locally; each rank injects its own elements' ΔM into its local node copy; the solve-time sum reconstructs the exact physical total. `dtTarget` is a user input (same on all ranks) ⇒ no dt desync.
- **Proof:** `Ladruno_implementation/mass_scaling_mpi/` — 1D fixed-free bar, fine zone straddling the central shared node (elements 10/rank0 + 11/rank1 BOTH inject into shared node 11). `np=1` (whole bar) vs `np=2` (split) tip-disp matches `max |abs diff| = 0.000e+00` over 150 steps (36.6% added mass). The three lumped integrators' "(3) … not mass-reduced across ranks" warnings + `CentralDifferenceSMS.h` scope comment were CORRECTED to state the truth (2026-06-21).
- **CONSISTENT (Olovsson) is NOW parallel-correct too (V5, 2026-06-21):** the serial `consistentPCG`/`consistentMatVec` ARE rank-local (local `res^z`/`p^Ap`, no shared-DOF `M̄` exchange) — so a real **distributed PCG** was added (`consistentParPCG`/`consistentParMatVec` + `LadrunoConsistentRefine.h`). See the dedicated V5 quirk below. The consistent-variant warnings were corrected to state the parallel path.
- **Tcl reaching SMS (FIXED):** `integrator CentralDifferenceSMS …` USED to error "No Integrator type exists" in the *Tcl* `OpenSees.exe`/`OpenSeesMP.exe` — SMS was registered only in the interpreter/openseespy layer (`SRC/interpreter/OpenSeesCommands.cpp`), NOT the legacy `SRC/tcl/commands.cpp` `specifyIntegrator()` parser, even though the Tcl splash banner advertised it. NOW WIRED: all 6 SMS integrators have `else if` branches in `specifyIntegrator()` (→ the same `OPS_*()` factories, null-guarded), mirroring the existing `CentralDifferenceLadruno` Tcl wiring. The OPS_ arg-parsing plumbing is already set up for the Tcl `integrator` command (CentralDifferenceLadruno/ExplicitBathe already use it), so the SMS factories read their args identically. Smoke-tested via `OpenSees.exe Ladruno_implementation/mass_scaling_mpi/sms_tcl_smoke.tcl` (all 6 build + step). **CI caveat: the 2-rank `mpiexec` test is single-process-CI-ungated — local validation only.** Learned 2026-06-21, T-MPI lumped validation + SMS Tcl wiring, [[36_ladruno_selective_mass_scaling_adr]].

### A shared-OpenSeesLIB class CANNOT `#ifdef _PARALLEL_INTERPRETERS` or reference `MPIDiagonalSOE` — that define + that SOE only exist in the MP executables. Use a `LinearSOE` base virtual.
- **The architecture trap (ADR-38 V5):** the consistent SMS integrators live in `OPS_Analysis`→`OpenSeesLIB`, which is compiled ONCE with NO parallel define and linked into BOTH the serial and MP binaries. `_PARALLEL_INTERPRETERS` is defined ONLY on the two special interpreter TUs (`OpenSeesCommands.cpp`/`OpenSeesMiscCommands.cpp` → `OPS_InterpPyCmds_MP`; see `CMakeLists.txt:~634`), and `MPIDiagonalSOE.cpp` is compiled PER-MP-TARGET (`CMakeLists.txt:~860/930/1031`), never into the shared lib. So in a shared-lib TU: (a) `#ifdef _PARALLEL_INTERPRETERS` is ALWAYS false (would compile the MP path out of every build); (b) `dynamic_cast<MPIDiagonalSOE*>` / any direct reference fails to LINK in the serial binary (the typeinfo/vtable isn't there).
- **The fix:** put the capability on the `LinearSOE` base as virtuals defaulting to a serial no-op (`isDistributedDiagonal()`, `getScalingDiagonalA()`, `assembleSharedSum(Vector&)`, `globalReduceSum(double)`); `MPIDiagonalSOE` overrides them. The integrator dispatches polymorphically through `LinearSOE*` at RUNTIME — works in every build, no MPI symbol in the shared lib. This is the idiomatic OpenSees pattern (LinearSOE already has many such optional virtuals). Cost: touching `LinearSOE.h` forces a wide recompile (it is included nearly everywhere) — a one-time hit. Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]].

### Distributed CG for a sum-assembled (additive) operator: one weight `wᵢ=1/multiplicityᵢ` fixes BOTH the replicated-diagonal double-count AND the dot products — no DOF ownership needed.
- **Setup:** the MPI diagonal SOE's `getDiagonalA()` is the GLOBAL (cross-rank-summed) lumped diagonal, so its value at a shared DOF is REPLICATED on every sharing rank. The off-diagonal `M̄ₑ` blocks each live on one rank. To matvec `y=M̃x` by local-compute + `assembleSharedSum` (sum shared entries across ranks): apply the diagonal WEIGHTED by `wᵢ` and the off-diagonal in FULL — after the assemble, the shared diagonal sums to `mult·(1/mult)·full = full` once, while the off-diagonal accumulates over all element-owning ranks. The SAME `wᵢ` makes global inner products count each shared DOF once: `⟨x,y⟩ = allreduce(Σᵢ wᵢ xᵢ yᵢ)`. Compute `multᵢ = assembleSharedSum(ones)` once (1 on purely-local DOFs ⇒ `w=1` there). With `w≡1` + no-op assemble + identity reduce it collapses EXACTLY to the serial PCG, so `np=1` reproduces serial and the serial path is untouched. **Deadlock-avoidance:** every convergence/participation decision (residual norm, `gActive`, `gBad`) MUST be a GLOBAL reduction so all ranks run the identical iteration count and call the per-iter collectives in lockstep; a rank with zero local scaled elements still participates if ANY rank scales (its shared DOFs feed a neighbour's `M̄`). Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]] V5.

### `build.bat` honored only the FIRST target arg (`set MODE=%1`) — `build.bat OpenSeesPy OpenSeesPyMP` silently built ONLY serial. And a parallel A/B compare can falsely PASS on TWO diverged runs.
- **The build trap (FIXED):** `Ladruno_scripts/build.bat` took `set "MODE=%1"`, so extra target args were dropped — the MP `.pyd` was never rebuilt and a STALE pre-edit binary ran (old warnings + the un-scaled fallback → the run diverged). Fixed to route the whole non-`clean`/`rebuild` arg list (`set MODE=%*`). Symptom to recognize: the run prints OLD `opserr` warning text you already changed, or the built `.pyd` mtime is older than your edits — always confirm the binary is fresh and the build log shows `Step 4: Building targets: <your target>` + your changed TUs recompiling.
- **The validation trap (FIXED):** a stale-binary divergence made BOTH `np=1` and `np=2` overflow to the SAME garbage (`-2.1e+179`), so a tip-disp A/B comparator saw `diff=0` and FALSELY PASSED. A parallel correctness compare MUST reject non-finite / unphysically large output BEFORE comparing (`compare*.py` now guard `|disp|>1.0` and `isfinite`). Cross-check against an INDEPENDENT serial reference (serial `DiagonalSOE`+`consistentPCG`), not only `np=1`-vs-`np=2` of the same new code. Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]] V5.

### `MPIDiagonalSOE::assembleSharedSum` reuses structures built by the FIRST `solve()` — calling it pre-solve SILENTLY drops the cross-rank sum (now tripwired).
- The neighbour exchange (`myActualNeighborsBsToSend`/`myNeighborsSizes`/posloc) and the factored GLOBAL `getScalingDiagonalA()` (= `1/mass` summed across ranks) are produced by `MPIDiagonalSolver::solve()`'s first (`notSet`) pass. The distributed consistent PCG's correctness therefore depends on the implicit invariant **first `solve()` precedes first `refineAccel()` on every rank** — true for all 3 shipped consistent integrators (CDL starter + ExplicitBathe/LNVD both sub-steps all do `solve()`→`refineAccel`). A new explicit integrator using the consistent path MUST preserve that order. A one-time tripwire warning now fires if `assembleSharedSum` runs with neighbours but un-built buffers (was a silent no-op → wrong answer). This was the unanimous residual flag from the 4-lens adversarial review (which found ZERO actual bugs). Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]] V5.
### LadrunoConcrete3D `-eta` (Duvaut–Lions): `dt≤0` falls back to INVISCID (β=1), NOT the elastic β→0 limit; and the oracle's deep-compression apex chaos forces stress-controlled gate states
- **dt-fallback (deliberate):** the relaxation factor is `β = dt/(η+dt)` only when `η>0 AND dt>0`; otherwise `β=1` (pure inviscid). Mathematically `dt→0` is the elastic limit (`β→0`, frozen trial), but operationally a missing/zero time increment (static, no pseudo-time) must NOT silently turn the material elastic — so `_dl_beta`/the kernel gate to inviscid. Consequence: like `-implex`, `-eta` is **inert without a positive `ops_Dt`** (transient, or uniform LoadControl/pseudo-time); it shares the IMPL-EX `dt=ops_Dt` caveat (`DisplacementControl` λ-increment is non-uniform near limit points). The elastic limit is reached only via large `η` (β→0 with finite dt), which is the correct knob. The wrapper applies `-eta` in the Tier-1 path only; under `-implex` the implicit solve stays inviscid (combined mode deferred, warned).
- **Gate-construction trap (the expensive one):** building a plastic committed state by uniaxial-**STRAIN** compression drives the numpy ORACLE into its known deep-compression apex regime (`kp<0`, the return is spuriously judged ELASTIC, `gap=σ̄_tr−σ̄_inv=0`). A viscous test built that way passes **tautologically** — the relaxation has nothing to relax (trial==inviscid). The eta gates (oracle PV5b/PV6, g++ NETA non-tautology check) therefore build the plastic state via the uniaxial-**STRESS** driver `drive_damaged_unified` (replay its solved `eps_lat` as a diagonal 6-strain path), which stays on the compressive meridian and yields a genuine plastic `gap`. Same family as the damaged-oracle "use STRESS-controlled, not confined-STRAIN" rule (#287). The g++ B6 block asserts `max viscous-inviscid gap > 1e-3` so an all-elastic (tautological) fixture can't pass. Learned 2026-06-20, `-eta` C++ port + wrapper, [[31_ladruno_concrete3d_adr]].
### LadrunoConcrete3D P2f: the real CDPM2 beta_c (Eq.50) makes compression MUCH more ductile; omega_c heals on unload (implicit re-solve); C2's Gc integral needs an analytic exp tail
- **beta_c is NOT ~1 in monotonic compression (it's ~0.058).** The P2b/P2c monotonic slice dropped beta_c (=1); the real beta_c = ft*qh2*sqrt(2/3)/(rho_bar*sqrt(1+2Df^2)) ≈ ft/(fc*sqrt(1+2Df^2)) ≈ 0.058 at the uniaxial-compression peak (ft/fc=0.1, Df=1). It scales the PLASTIC-strain part kappa_dc1 (Eq.48), suppressing it ~17×, so compression is MARKEDLY more ductile — the post-peak nominal stress differs by ~23 MPa (fc=30) vs beta_c=1 (real beta_c softens to ~0 only by ~-0.3 strain, vs ~-0.05 for beta_c=1). This is the faithful CDPM2 behaviour (user decision 2026-06-20: restore beta_c always). **Gotcha for gate design:** at DEEP strain omega_c saturates to ~1 for BOTH beta_c values, so the wc GAP is tiny there — use the post-peak STRESS gap (tens of MPa, unambiguous) as the non-tautology metric, not the wc gap. Also the damage-ONSET region is step-size sensitive (a coarse step overshoots the sharp onset, over-accumulating ‖Δεp‖ → more early damage); use a resolved path (~5000 steps over -0.3) for compression-damage gates.
- **omega_c HEALS on elastic unload (the cyclic gap).** The oracle solves omega_c IMPLICITLY against the CURRENT effective stress every step ((1-wc)*D = fc*exp(-eps_i/eps_fc), eps_i carries wc). On an elastic UNLOAD the effective stress D drops with the history (kdc1/kdc2) frozen, so the solved wc DECREASES — the material spuriously heals and the reload modulus returns to ~E. beta_c (the damage RATE factor) does NOT fix this; a cyclic-correct omega_c must be driven by the MONOTONE history (omega_c <- max over the path), a separate fix touching every driver + the committed state + the C++ kernel (the P2f monotone-omega_c slice). Reported as the run_p2f_gate F4 diagnostic, not gated.
- **C2's "by construction" Gc integral truncates with the ductile beta_c — add the analytic exp tail.** The omega_c solve enforces |sig_c_nom| = fc*exp(-eps_i/eps_fc) EXACTLY, so int_0^inf |sig| d eps_i = fc*eps_fc = Gc/lch. Over a FIXED compression strain the more-ductile beta_c response doesn't reach eps_i->inf at small lch (large eps_fc), truncating the integral (4.58 vs 5.0 at lch=50). FIX: add the analytic tail int_{eps_i,last}^inf = eps_fc*|sig_c_last| to the trapezoidal path integral — exact (5.000 at all lch), no need for an absurd strain range. **The analytic damaged tangent's ∂beta_c/∂eps** is a composite micro-FD THROUGH the return map (beta_c depends on BOTH rho_bar(sig_bar) AND qh2(kp)) — the C++ kernel `damagedTangent` MIRRORS that composite FD (re-runs `returnMapTensor` per component; needs a forward-decl since returnMapTensor is defined after damagedTangent) rather than deriving ∂kp/∂eps analytically. **beta_c is oracle+C++ in ONE PR (#321), NOT split like eta** — it is always-on (no default-off), so the regenerated g++ fixture's compression cases diverge from a beta_c=1 kernel and the byte-check FORCES oracle+kernel together (like #249); no wrapper/serialization change. Learned 2026-06-20, P2f beta_c oracle + C++ port (#321), [[31_ladruno_concrete3d_adr]].
### LadrunoConcrete3D P2g: damage healed on unload because omega was solved against the LIVE drive — fix = drive omega with the running MAX (monotone). The kappa histories were already monotone
- **Root cause:** the dual damage solved `omega` from `(1-omega)*D - f*exp(-(kd1+omega*kd2)/eps_f) = 0` with `D` = the **LIVE** extreme effective principal stress (`sig_t_drive`/`sig_c_drive`). The inelastic histories `kdt1/kdt2/kdc1/kdc2` ALREADY accumulate only when loading (monotone), so the live `D` was the **only** non-monotone input. On an elastic UNLOAD `D` drops with the histories frozen, `F(0)=D−f*exp(−kd1/eps_f)` goes ≤0, and the bracketed solve relaxes `omega` back toward 0 — the crack spuriously HEALED and the reload modulus returned to ~E (the #321 F4 diagnostic). Wrong: CDPM2 states `omega = omega(kappa_d)`, a function of the MONOTONE history only.
- **Fix (P2g, #325):** track the running MAX of each channel's drive stress (`sigt_max`/`sigc_max`, two new committed scalars) and solve `omega` against THAT. `omega` is then monotone-nondecreasing; on unload the max (and histories) are frozen ⇒ `omega` frozen ⇒ the nominal stress unloads along the degraded **secant** `(1-omega)*sig_bar`. On ANY monotonic path `max == live`, so the change is **byte-identical** to the pre-P2g drivers (DT1/DT2 reduce-to-P2a/P2b, all P2e/P2f gates, every monotonic g++ case — verified). Mirror in BOTH uniaxial reference drivers AND `drive_damaged_unified`/`damaged_step_tensor`/the C++ `damagedUpdate` so the whole damage subsystem is uniformly no-heal.
- **Tangent on unload = the SPD secant.** The analytic damaged tangent's `−sig⊗∂omega/∂eps` rank-update must VANISH on an unloading channel (frozen drive ⇒ `∂omega/∂eps=0`), leaving `D_dam:C_eff`. Gate the eigenprojection gradient `dD*_deps` on a per-channel "advancing the max" flag (`Dt >= committed sigt_max`); during loading `max==live` so the gradient is unchanged (P2e gates hold), on unload it is zeroed. The unload secant is **SPD** (`lambda_min>0`) — the well-conditioned branch, in contrast to the INDEFINITE Tier-1 loading tangent (gate TD2). Self-verified analytic==numerical-FD (~1e-11) at an unload state.
- **Verification trap (discriminating fixture):** a monotonic committed state has `sigt_max == live`, so it canNOT distinguish a correct kernel from one that ignores `sigt_max`. The g++ byte-check needs a CYCLIC committed state — build tension into softening, THEN elastically unload so `sigt_max > live drive`, then probe a further unload step (`dmg_cyclic_unload`). A live-drive kernel heals there and diverges; the monotone kernel byte-matches (`nom_sig_err=0`). Also serialize the two new scalars (`LC3D_NDATA` +2, send/recv) — else a parallel/database round-trip resets the drive max and re-heals. Gated by oracle `run_p2g_gate` G1-G6 + pytest `test_p2g_monotone_damage_gate` + element `test_cyclic_no_heal_unload`. Learned 2026-06-20, P2g monotone-omega no-heal (#325), [[31_ladruno_concrete3d_adr]].
### LadrunoConcrete3D P2h: compression->tension damage temper (-ctTemper) — the plastic-strain's OWN positive part is NOT a valid tensile shield (dilatancy); project onto the tensile-STRESS frame
- **The problem (DT5):** literal CDPM2 (Eq.43, `κ̇_dt = ε̃̇`, NO `(1-α_c)`) accumulates the TENSILE damage history during COMPRESSION too (the equiv strain `ε̃` grows in compression), so a compression excursion pre-damages a subsequent tension reload to ~0. User decision (2026-06-20): offer BOTH temper modes behind a flag `-ctTemper {none|alphat|proj}` (default `none` = literal CDPM2, byte-identical). A tensile weight `w_t` scales the kdt1/kdt2 accumulation: `alphat` → `w_t=1-α_c`; `proj` → tensile-stress-projected plastic-strain fraction.
- **alphat is the clean one:** `w_t=1-α_c` leaves BOTH monotonic backbones EXACT (α_c=0 in pure tension ⇒ w_t=1; the compression backbone rides the kdc channel untouched) and removes ONLY the cross-coupling. Restores tension-after-compression to ~ft. RECOMMENDED.
- **proj GOTCHA — do NOT use the plastic strain's own `<Δε_p>+`:** in compression the DILATANT non-associated flow (volumetric `m_v=Df·m0/(√3 fc)>0`) makes the LATERAL plastic strains POSITIVE, so `‖<Δε_p>+‖/‖Δε_p‖` stays LARGE in compression ⇒ no shield (the first proj attempt did nothing: tac peak stayed 0). The correct "tensile-plastic-strain projection" projects `Δε_p` onto the POSITIVE effective-STRESS principal directions (`w_t = ‖P+ Δε_p‖/‖Δε_p‖`, `P+` = eigenprojection where `σ̄_a>0`): compression has no tensile-stress directions ⇒ `w_t=0` (full shield). proj lightly softens the monotonic tension backbone (the loaded axial carries <100% of `‖Δε_p‖`).
- **Tangent + verification:** `∂w_t/∂ε` = `-∂α_c/∂ε` (alphat, analytic, reuses `dac_deps`) / composite micro-FD through the return map (proj). Oracle `run_p2h_gate` H0-H4 (H1 alphat restores + byte-identical backbone; H2 proj restores; H4 analytic==numerical tangent both ~6e-9). **g++ fixture-case GOTCHA:** a compression→tension committed state at the tension-softening tail has `wt=1` (degenerate tangent, the ORACLE fails its own analytic==numerical there too — NOT a C++ bug) and the unload overshoots into tension; use a **BIAXIAL** loading state (axial tension + lateral compression, partial α_c, 0<wt<1) for the alphat byte-check and PURE TENSION for proj — both smooth + discriminating (alphat wt 0.78 vs none 0.82). DMG fixture line gains the ctTemper int; wrapper serializes it (`LC3D_NDATA` +1). No new classTag/banner. Learned 2026-06-20, P2h ctTemper oracle+C+++wrapper (#327), [[31_ladruno_concrete3d_adr]].
### LadrunoConcrete3D P5b: the confined-fiber (BeamFiber) view condenses {11,22,12} vs the hoop on the EFFECTIVE lateral stress; the consistent tangent MIXES the damaged + effective operators
- **The view (§4.6 "Mander by mechanism"):** a `DIM_BEAMFIBER` mode (`getCopy("BeamFiber")`, order-3 retained {00,01,02} = axial + 2 transverse shears) consumed by the stock `NDFiberSection3d`. The lateral block {11,22,12} is NOT the `σ_22=0` single-index condensation the plane views use — it is a 3-unknown nested Newton (kernel `driveConfinedFiber`) balancing the passive hoop: `sigEff_11+σ_hoop(ε_11)=0`, `sigEff_22+σ_hoop(ε_22)=0`, `sigEff_12=0`, with `σ_hoop(ε_lat)=min(K·ε_lat,fy)` tension-only. `hoopK=0` ⇒ free reduction (plain BeamFiber). So the wrapper keeps `condense=false` and a separate `confined=true` flag selects the dedicated kernel call.
- **Balance on EFFECTIVE, not nominal (matches the shipped P5a oracle):** the lateral Newton targets the EFFECTIVE (undamaged) lateral stress + hoop, NOT the nominal damaged stress. This is a deliberate match to the shipped P5a oracle `drive_confined_fiber` (so the g++ `CFIB` byte-check passes without reopening #341). At/near the compression peak the damage is ~0 (effective-plasticity is monotonic, peak=onset κp=1, damage starts AT peak) so nominal==effective there ⇒ the Mander match (evaluated at peak) is unaffected by the choice. Post-peak the hoop reacts against the undamaged effective lateral (mildly optimistic). A nominal-balance variant is a future P5c.
- **The condensed consistent tangent is a MIX:** `C_RR = Cdam_RR − Cdam_RL · (Ceff_LL + H)⁻¹ · Ceff_LR`. The output rows use the DAMAGED tangent `Cdam` (the nominal stress the element consumes); the constraint linearization `∂g_L/∂ε` uses the EFFECTIVE tangent `Ceff_LL` (the balance is on effective stress) + the hoop stiffness `H=hoopStiffness` on the lateral-normal diagonal. A pure `Cdam` condensation would be wrong (the constraint is on effective). Reduces to a plain symmetric static condensation where ω→0 (`Cdam→Ceff`). `getInitialTangent` condenses the elastic C over {1,2,4} (hoop slack at ε=0 ⇒ H=0).
- **Byte-check structure:** the oracle gained a single-step primitive `confined_step` (an exact extract of the `drive_confined_fiber` per-step body; gate F5 asserts looping it == the path driver to 0.0, so the shipped F1–F4 stay untouched). The `CFIB` fixture dumps the committed kernel State (6-Voigt diagonal; `epl_prev` is NOT dumped — `damagedUpdate` derives it via `plasticStrain6(in.sigEff,in.eps)` == the oracle's tracked `epl_prev`) + one axial increment + hoopK/fy; C++ `driveConfinedFiber` reproduces nominal axial stress + p_conf + the converged lateral strain to ~1e-11 (reduce-to-free / confined / hoop-yield). No new classTag; send/recv +2 (hoopK,hoopFy). Learned 2026-06-21, P5b confined-fiber C++ view, [[31_ladruno_concrete3d_adr]].
### ADR-30 P6 — LtML condition gate + frozen-Ccr staleness; and a Windows stale-object build trap
- **The projection's singularity guard must be a CONDITION-NUMBER gate, not an exact-pivot test-solve.** `LadrunoConstraintProjector::buildMass` used to probe `LtML = LᵀML` with a single `Matrix::Solve(e0,x0)`, which only fails on an EXACT zero pivot. A *near*-singular `LtML` (a barely-dependent retained direction, hugely disparate tied masses, a near-redundant constraint) passes the solve, then `project()` amplifies round-off into a garbage acceleration — silently wrong. P6 estimates `cond = λmax/λmin` of the SPD `LtML` via a self-contained cyclic-**Jacobi** symmetric eigensolve (`ladrunoSymEigJacobi`, ~50 lines, no LAPACK — deliberately, to sidestep the bundled-LAPACK `dsygv_`-missing gap, [[project_zonea_link_blocker]]); **refuse** above `1e12` (also catches exact-singular `λmin≤0`), **warn** above `1e8`. Jacobi `t`-form root `θ=(aqq−app)/(2apq)`, `t=sign(θ)/(|θ|+√(θ²+1))` (smaller root) is exact even for a 1e13:1 eigenvalue spread.
- **To TEST the cond gate you need a group with nRet≥2** — a single `equalDOF` makes `LtML` a 1×1 scalar (cond≡1, never trips). Use a **multi-master `equationConstraint`** (`equationConstraint(cN,cD,1.0, r1N,r1D,-0.5, r2N,r2D,-0.5)` ⇒ `u_c=0.5u_{r1}+0.5u_{r2}`, two retained DOFs in ONE group). A retained-mass ratio of `1e13` gives `LtML≈[[1e13,0.25],[0.25,1.25]]`, `cond≈8e12` → refused; `1e10` → `cond≈8e9` (warn, runs). A tiny ROTATIONAL mass on a rigidLink-beam master does NOT ill-condition `LtML` (the slave's translational mass props the retained `rz` up through the lever arm — that is the projection working, not a defect).
- **Frozen-Ccr staleness: the lever arm that goes stale is the master-ROTATION → slave-TRANSLATION cross term, NOT any rotational tie.** A direct `equalDOF` on `rz` (rotation→rotation, coeff 1) is exact under any rotation and must NOT be flagged; a rigidLink-beam / rigidDiaphragm couples master `rz` into slave `ux,uy` (the offset lever arm) and DOES drift past ~0.1 rad. `flagRotMonitor` requires `masterDof` rotational AND `slaveDof` translational (rotational test by OpenSees convention from `ndm=node->getCrds().Size()`, `ndf`: `(ndm==2&&dof≥2)||(ndm==3&&ndf≥6&&dof≥3)`; translational = `dof<ndm`). The warn-once latch + per-step drift read live in the handler's `applyLoad()` (already the per-step hook).
- **NEVER use `ops.logFile` to capture `opserr` in a pytest — it POLLUTES every later test in the shared openseespy process. Use pytest `capfd`.** `ops.logFile(path,"-noEcho")` calls `opserr.setFile(path,...)` and there is NO command to restore `opserr` to the console — so once any test redirects `opserr` to a file, EVERY subsequent test in the same process (openseespy is ONE process per pytest session) that reads `opserr` from stderr gets an EMPTY string (`assert 'X' in ''`). The P6 staleness test first used `logFile` (it passed in isolation AND in the projection-only battery, since nothing captured `opserr` after it) and on the FULL Zone-A run it silently broke 7 downstream mass-scaling / consistent-PCG tests that capture `opserr` via `capfd` — **green locally, red on CI**. **Fix: capture with pytest's `capfd` fixture** (`out = capfd.readouterr(); txt = out.err + out.out`) — file-descriptor level, per-test, non-polluting, does NOT redirect `opserr` away (it is also how the mass-scaling tests capture their `-verbose`/warning reports). **General lesson: any global OpenSees stream/state mutation in a test (logFile, a left-open recorder, defaultUnits) leaks across the whole shared-process suite — before trusting a new test's green, run it TOGETHER with the broader integrator/mass-scaling set, not just its own file.** Learned 2026-06-21, P6 (#337 → #338 test-fix follow-up).
- **WINDOWS STALE-OBJECT BUILD TRAP: `build.bat OpenSees OpenSeesPy` can link a STALE object for a just-edited `.cpp` into the `.pyd`.** After editing `LadrunoConstraintProjector.cpp`, a combined `OpenSees OpenSeesPy` build (exit 0) produced an `opensees.pyd` that still ran the OLD `buildMass` (cond gate absent → an ill-conditioned model that should refuse ran clean). Re-running `build.bat OpenSeesPy` ALONE (after touching the file again) recompiled it and the gate appeared. Symptom: a behavior change you just implemented is absent despite a green build. **Fix/avoid:** when a single-file edit's behavior is missing after a multi-target build, rebuild the **single** consumer target (`OpenSeesPy`) to force the recompile, or wipe the build dir. CI (fresh g++ from scratch) is unaffected — this is a Windows incremental-build quirk only. Learned 2026-06-21, P6 (#334-followup).

### ADR-41 C1 mortar kernel: the constant-pressure patch test does NOT catch the things you'd expect — it reduces to `Σφ=1` and passes even on broken geometry; the clip's real job is exactness, and silent area-bias hides under it
- **Bites:** trusting the partition-of-unity (`Σ_K M_IK == Σ_J D_IJ`) and constant-pressure patch tests as a proof of mortar correctness. Both reduce **per-Gauss-point** to `Σ_K φ_K^m == 1` (master shape functions sum to 1 — *everywhere*, including extrapolated/out-of-bounds points), so they hold to ~1e-16 on tilted masters, distorted quads, slivers, OOB-projecting GPs, **even a single non-convex slave facet**. They are necessary, not sufficient. The C1 adversarial review (#369 follow-up) confirmed this on every pathological config.
- **Why:** the constant field `c` gives `M·(c·1)=c·(ΣM rows)` and `D·(c·1)=c·(ΣD rows)`; equal row sums ⇒ exact transfer, independent of the integration MEASURE. So a WRONG area Jacobian `J` (it cancels between D and M), a non-convex facet integrated over its convex hull, or a warped-facet measure error are all **invisible** to the patch test. What the patch test *cannot* hide, and what the clip actually buys, is reproduction of a **non-affine** master field: the integrand `N_I^s φ_K^m` has kinks at master-facet edges, and a clip-blind quadrature is ~40% biased there (oracle T5), converging to the clipped value only under heavy refinement.
- **Two silent-wrong-answer traps the gates miss (now guarded/documented):** (1) a **non-convex** (concave/bow-tie) slave facet projects to a non-convex polygon → Sutherland-Hodgman integrates it over its **convex hull** (e.g. a concave quad of true area 0.40 integrates to 1.80, 4.5×) AND the C++ `clipPolygon` truncates at the `MAXV` clamp, diverging ~30% from the numpy oracle. FIX: `isConvex2` guard in `integratePair` (and the oracle's `clip_subtris`) **refuses** the pair (`status=-1`). (2) a **warped (non-planar)** slave facet biases the `J = A_aux/|n_s·n0|` flat-facet area ratio (~0.7% at a 0.3·edge out-of-plane lift) — *not* refused (deformed meshes legitimately warp), but DOCUMENTED in the header + oracle T10; the exact fix is a per-sub-triangle slave Jacobian `|g1×g2|`, deferred to **C2/ADR-47**.
- **Status (2026-06-23):** guards shipped in `LadrunoMortarKernel.h` + mirrored in `proto_c1_mortar.py` (oracle 30/30); `inverseIsomap2D` also gained a convergence flag (skip the GP rather than integrate a wrong root). See [[48_ladruno_contact_capstone_adr]] C1 row, [[_adr41_c1_design]], #369.

### The `LadrunoContactBucketSort` superset contract holds for POINT slaves only — querying by a slave-FACET centroid silently drops overlaps (coarse-slave/fine-master)
- **Bites:** the ADR-41 C2.1 mortar pairing. The NTS broad phase (`LadrunoContactBucketSort::Grid`) was designed + proven (`proto_bucket_sort.py`) for slave **NODES** (points): `candidates(point)` returns the segments registered in the point's bucket ±1 neighbour, a guaranteed SUPERSET of near pairs. The mortar lane pairs slave **FACETS** with master facets; querying the grid by the slave-facet **centroid** is NOT a superset for an EXTENDED facet — the cell size is the **median MASTER-segment diagonal**, so a slave facet wider than ~3 master buckets loses its outer overlaps. Measured at the C2.1 code gate: a 2×2 slave facet over an 8×8 master mesh settled at `-1.78e-4` vs the full-coverage `-1.0e-4` — **~44% of the contact area silently dropped** (the kernel clip can only reject candidates it's GIVEN; it can't recover a master the broad phase never offered). Symptom: a coarse-slave/fine-master mortar interface transmits too little contact force, with NO error.
- **Why:** `Grid` registers each master facet in every bucket its corner-span touches (`LadrunoContactBucketSort.h` span-registration), and `lmax = lmaxFrac·median(masterDiag)`. The query is point-based; the slave facet's spatial extent is never considered.
- **Workaround/status (2026-06-23):** C2.1 ships **brute-force** mortar pairing (every master facet a candidate per slave facet; the clip rejects non-overlaps) — correctness-first, exactly as NTS shipped brute force at P2b-1 before the P2.5 bucket sort. A slave-aware mortar broad phase (query by the slave facet's corner-span, or size the cell from `max(masterDiag, slaveDiag)`) is the deferred "mortar P2.5" optimization. General lesson: a broad-phase "superset" proof is tied to the QUERY primitive (point vs extended) — re-validate it before reusing a point-proven grid for facet/segment queries. Found by the C2.1 adversarial code gate (#374). See [[48_ladruno_contact_capstone_adr]] C2 row.

### Per-facet mortar adapters reading a RUNNING global gap give an order-dependent (non-Newton) residual at shared slave nodes
- **Bites:** the ADR-41 C2.2 Uzawa ALM. The handoff ([[_adr41_c2_design]] §crux) recommended the augmented pressure `p_I = min(0, λ_I + epsN·ḡ_I^global)` read the GLOBAL weighted gap (summed over every facet a slave node touches) "lagged by one iteration." But `LadrunoContactFE` is PER-FACET and `NewtonRaphson::solveCurrentStep` forms the residual sweep **facet-by-facet within one `formUnbalance`**. If each adapter accumulates its facet's `g̃` into a Domain-side running sum and then reads it back, the FIRST facet to touch a shared node sees only its own contribution, the SECOND sees two, etc. — so a shared node gets a DIFFERENT pressure in each facet's force assembly. The residual is then a function of facet **evaluation order**, not just the displacement ⇒ not a clean `R(u)` ⇒ Newton's tangent is inconsistent and the solve goes **singular** (`FullGenLinLapackSolver U(i,i)=0`). Symptom: the matched (1-facet/node) C2.1 test still converged, but the non-matched + coarse-slave/fine-master cases (shared nodes) failed to converge after the C2.2 change.
- **Why:** the oracle's "lag" (`proto_c2_alm.py` T7c) was a **frozen** per-sweep value — every node read the SAME gap for the whole residual sweep, refreshed once per Newton iteration. Reproducing a frozen-per-sweep value in a per-facet adapter needs sweep-boundary detection (no per-residual-sweep hook exists on the contact engine), and seeding it to zero kills the first sweep's contact (→ singular tangent).
- **Resolution (shipped C2.2):** keep the force/tangent on each facet's **LOCAL** gap `p_I = min(0, λ_I + epsN·ḡ_I^facet)` — deterministic, exactly the C2.1 penalty — plus the per-GLOBAL-node multiplier `λ_I` (which assembles globally for free: `Σ_facets D_KI^facet λ_I = D_KI^global λ_I`, λ shared). The GLOBAL weighted gap IS still accumulated on the Domain (idempotent delta updates keyed `(contactTag, nodeTag, feTag)`), but ONLY for the once-per-`commit()` Uzawa update `λ_I ← min(0, λ_I + epsN·ḡ_I^global)` and the `ladrunoMortarPenetration` query — never read back into the same sweep's force. At convergence the penalty term → 0 and the consistent global `λ` carries the load, so the result is variationally consistent + epsN-independent. Pinned oracle-first: `proto_c2_alm.py` T8 (local-gap force + global-λ Uzawa → epsN-independent penetration, 28/28). General lesson: a "lag-by-one-iteration" scheme is only well-posed if the lagged quantity is FROZEN across the whole residual sweep; a value mutated *during* the sweep makes the residual order-dependent. Found while transcribing C2.2 (the non-matched battery regression caught it). See [[48_ladruno_contact_capstone_adr]] C2 row, #375.

### Mortar tangential SLIP must come from DISPLACEMENTS, not positions — the closest-point projection makes the weighted relative POSITION purely normal
- **Bites:** the ADR-41 C3.1 mortar friction. The natural-looking weighted relative position
  `r_I = Σ_J D_IJ x_s,J − Σ_K M_IK x_m,K` (= `∫N_I(x_s − x_m(ξ̄)) dΓ`) is **purely NORMAL**: `n·r_I = g̃_I`
  (the weighted normal gap) and its TANGENTIAL part is ≈ 0, because the closest-point projection `ξ̄`
  places the master point directly "under" the slave point (`x_s − x_m(ξ̄) ∥ n` by construction). So a
  friction slip built from positions is ~0 even when the slave has slid a finite tangential distance —
  the return map sees `gTeff ≈ 0`, stays in STICK, and assembles ZERO friction force (symptom: a driven
  block accelerates at the frictionless `a = Q/m`, friction silently inert, to 1e-13). Verified by a
  stderr probe: `gTeff=(2.8e-17, …)` at a step where the slave had displaced `x=1e-3`.
- **Why:** mortar inherits the NTS lesson — the ADR-39 `SEGMENT` path's `segmentActive` ALREADY documents
  "the closest-point projection makes (x_s − x̄) ∥ n, so POSITIONS carry NO tangential information; the slip
  is the slave DISPLACEMENT minus the interpolated master DISPLACEMENT at the projection: `d = u_s − Σ N_i u_i`."
  The C3.1 first draft re-made the position mistake the NTS path had already solved.
- **Fix (shipped C3.1):** build the slip from DISPLACEMENTS — `r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K`
  (`u = getTrialDisp()`), tangential part `/a_I`, minus the engagement origin `gT0_I`. This is the `D/M`-
  weighted generalisation of the NTS `u_s − Σ N_i u_i`. The normal gap still uses positions (it IS the
  normal projection); only the tangential slip switches to displacements. General lesson: in any
  closest-point-projected contact, the normal gap is a POSITION quantity and the tangential slip is a
  DISPLACEMENT quantity — they are not interchangeable. Found while bringing up C3.1 (the driven-block
  gate caught it). See [[_adr41_c3_design]] §mechanics step 1, #377.

### Mortar friction committed slip is last-writer-wins (order-dependent) at SHARED slave nodes — fenced to matched/explicit C3.1, must be guarded before non-matched friction
- **Bites:** ADR-41 C3.1 mortar friction at a slave node shared by ≥2 (slave-facet, master-facet) pairs.
  The per-global-node committed slip `st.gpTtrial` (`LadrunoContactFE::addMortarFriction`) is a plain
  OVERWRITE: each facet visiting the node computes its OWN LOCAL `gbarT` (its own clip/projection) and the
  LAST facet evaluated in the residual sweep wins the committed slip. The *force* is still deterministic
  (every facet reads the same read-only committed `gpT`, so `R(u)` is clean — no singular solve), but the
  committed plastic slip carried to the next step depends on FE-tag ordering. The normal gap dodged this
  with an idempotent delta-accumulator keyed `(c,node,feTag)` (`accumulateMortarGap`); the friction slip has
  no equivalent because the slip is a return-map OUTPUT, not a linear accumulation.
- **Why it's fenced (for now):** C3.1 ships matched-facet + explicit (CDL) only — one facet per node, so the
  race never fires (the battery is matched). It is within the design's accepted "standard-basis LOCAL
  approximation at shared nodes" ([[_adr41_c3_design]]). But it is UNGUARDED and untested for non-matched
  friction. **Before C4 / non-matched frictional meshes:** add a shared-node friction regression + either a
  per-(node,feTag) slip reconciliation or an explicit area-weighted blend. Found by the C3.1 adversarial gate
  (MAJOR-1, #377). **C3.3 update (#379):** the tangential multiplier `λ_T` (committed from `lambdaTtrial`,
  written per-facet last-writer-wins) and `gpT` BOTH inherit this — unlike the normal `λ_N` Uzawa, which
  augments from the order-INDEPENDENT global accumulator `gtGlobal/aGlobal`. So the per-node `λ_T`/`gpT`
  reconciliation is the same single fix for the whole friction state. Still fenced to matched/explicit; the
  C3.3 gate (MINOR-1) re-confirmed it is inherited, not introduced.
- **C4 update (#381) — RESOLVED for the TIE path; STILL FENCED for FRICTION.** C4 mesh-tying hits shared
  slave nodes immediately (non-matching meshes are the whole point), so the pre-req had to be discharged
  before relying on it. The tie state (`λ_tie`, the full 3-vec relative displacement `r_I`) does NOT inherit
  the bug, because `r_I = Σ D u_s − Σ M u_m` is a **LINEAR accumulation** — not a return-map OUTPUT — so it
  uses the SAME order-independent global accumulator `λ_N` already uses (`accumulateMortarTie` delta-update
  keyed `(c,node,feTag)` → `rtGlobal`, Uzawa'd in `commit()` no-clamp). The FORCE reads each facet's LOCAL
  `r` (deterministic R(u), the C2.2 rule); the GLOBAL `r` feeds only the commit Uzawa + the `‖r‖` query.
  Pinned by `test_adr41_mortar_c4_1::test_c4_1_shared_node_order_independent` (a slave node shared by 2 tie
  facets ⇒ either facet order gives a BIT-identical converged solution) + oracle T6. **The FRICTION
  last-writer-wins (gpT/λ_T) is STILL fenced** — C4 is mutually exclusive with friction, so the tie never
  touches the friction slip; non-matching FRICTIONAL meshes still need the per-(node,feTag) slip
  reconciliation (an area-weighted blend, since the slip is a return-map output, not linear).

### Mortar friction gT0/engaged are captured in getResidual and NOT reverted — latent until the C3.2 implicit tangent
- **Bites:** ADR-41 C3.2 (NOT C3.1). `revertToLastCommit` drops only `gpTtrial=gpT` for mortar slots
  (`LadrunoContactDomain.cpp`); the engagement origin `gT0`/`engaged` (set once in `addMortarFriction`) are
  never reverted. A rejected Newton step that FIRST-engages a node latches `gT0` from the rejected trial
  config; the retry keeps that stale origin (`engaged` stays true) ⇒ a spurious stick offset. Identical to
  the shipped NTS SEGMENT behavior (which also doesn't revert `engaged`), so NOT a C3.1 regression, and
  **unreachable under C3.1's explicit-only path** (CDL never reverts mid-step). It goes live when the C3.2
  friction tangent lands and an implicit Newton step is rejected. **RESOLVED in C3.2 (#378):** `gT0`/`engaged`
  are double-buffered (`gT0committed`/`engagedCommitted`), promoted in `commit()` and restored in
  `revertToLastCommit()`. Found by the C3.1 gate (MAJOR-2, #377), fixed in C3.2.
- **2026-07-02 (contact-review P2): the SAME fix was finally BACK-PORTED to the NTS lane it was copied
  from.** The paragraph above ("identical to the shipped NTS SEGMENT behavior") documented the shared flaw
  but only excused it as not-a-C3.1-regression — and NTS friction is a first-class IMPLICIT path since
  P3.5 (#361), so the "unreachable under explicit" shield never applied there. NTS `FrictionState` now
  carries the same `gT0committed`/`engagedCommitted` double-buffer (commit promotes / revert restores).
  LESSON: when a gate fixes a defect on a COPIED lane, grep for the source lane the copy came from — a
  ledger sentence acknowledging "the sibling has it too" is a fix obligation, not an absolution. Same PR:
  `LadrunoContactDomain::revertToStart()` (hooked from `Domain::revertToStart` — `ops.reset()`) drops ALL
  contact path state (friction slip/origins, ALM λ_N/λ_T/λ_tie, edge signs, re-emit anchors/fp/trigger,
  NTS-force + nodal-mass caches; NormalField σ+frozen sign KEPT — reference-geometry datum, re-derivable
  identical). Pre-fix, a re-run after `ops.reset()` started from the previous run's committed slip gpT ⇒ a
  large spurious backward stick force at first contact — silently different from a fresh model. Gates:
  `tests/test_contact_review_p2_lifecycle.py` (failed-step retry ≡ never-failed reference bit-tight;
  reset re-run ≡ first run bit-tight — both FAIL pre-fix).

### B3 (P2b-2c): the LadrunoContact handler does NOT enforce non-zero SP (imposed displacement) — NOW WARNS
- **Bites:** any model that drives a `LadrunoContact` analysis by imposed nodal displacement
  (`ops.sp(node, dof, value)` with a non-zero value in a `pattern`) — e.g. a displacement-controlled
  indenter. **Verified:** a free-x node with `sp(node,1,0.05)` ends at `ux=0.000` under
  `constraints LadrunoContact`, vs `ux=0.050` under `constraints Transformation`. The contact handler is
  Plain-style — it REMOVES the constrained DOF (ID=-1) but never applies the non-homogeneous value (only
  Transformation/Penalty/Lagrange do). So imposed-displacement DRIVING silently does nothing.
- **NOT a regression:** stock `PlainHandler` does the IDENTICAL thing (`ux=0.000`) — it is the documented
  PlainHandler limitation. The contact handler MUST be Plain-style (it injects the contact FE adapters in
  `handle()`, which Transformation/Penalty can't do), so it inherits this. **Fix shipped:** the handler now
  emits the same non-homogeneous-SP WARNING as `PlainHandler` (it was SILENT before — worse than stock
  Plain), pointing at `DisplacementControl`. So a displacement-driven contact model fails LOUD, not zero.
- **How to drive instead:** external LOAD (`ops.load`, `LoadControl`), the `DisplacementControl` integrator
  on a FREE DOF (it augments the system with the load factor — works with the Plain/contact handler, no SP
  needed), or a fixed-geometry incompatibility (pre-set node positions, like `block_on_block`'s 1e-8
  pre-penetration). FULL imposed-SP support would need a Transformation-style contact handler (no current
  consumer; deferred). Found while building the Hertz benchmark (capstone B3 gate 2).

### Contact-review P3 (2026-07-02): τmax-only friction was an UNBOUNDED bond; validation choke points; mortar isomap scale fix
- **τmax-only (HIGH-2):** `-tauMax` without `-mu`/`-cohesion` is the unified cone `min(μN+c, τmax) =
  min(0, τmax) = 0` — a ZERO cone radius (free slip). The kernel's `cap≤0` branch returned the RAW
  ELASTIC traction ("byte-identity is the guard's job") with a consistent stick tangent, so the config
  silently became an UNBOUNDED elastic bond — the FE guards treat `tauMax>0` as a friction REQUEST and
  the handler auto-defaults `kt`, so it fell exactly between the tested branches (Tresca is tested as
  μ=0 WITH cohesion). FIX at two layers: the kernel `cap≤0` branch now FREE-SLIPS (zero traction, slip
  absorbs the motion, zero tangent block — `frictionReturnMap` + `frictionTangentBlock` +
  the numpy mirror in `proto_a1_friction.py`), and the command surface REFUSES `-tauMax`/`-edgeTauMax`
  without a `-mu`/`-cohesion` (the sanctioned shear-capped BOND is `-cohesion`, optionally + `-tauMax`);
  on the NTS lane `-tauMax` is refused outright as mortar-only (it was silently INERT there — the
  positional-μ NTS cone has no shear cap; adversarial-gate MINOR-1, fail-loud like `-geomtan`).
  cap>0 branches byte-unchanged (166k-case gate fuzz bitwise-identical + oracle byte-guards pass on
  BOTH pre/post kernels). Oracle: `contact_prototypes/proto_friction_validation.cpp` (12 checks,
  6 FAIL pre-fix). GATE-SURFACED FOLLOW-UP (pre-existing, not fixed here): mortar AUTO-orientation
  silently degenerates for COINCIDENT conforming facets — `orientDir = scen − mcen = 0` and the
  facet-normal flip test never fires, leaving the sign to winding luck; pass `-outward` for coincident
  interfaces (the shipped c2_1 tests always do). Candidate warn/refuse in the review PR-5 batch.
- **Validation choke points (previously silent):** duplicate CONTACT tags now refused across
  contact/mortar/contactPlane (`contactTagInUse` — the tag is the leading key of every path-state
  store: a duplicate ALIASED friction slots last-writer-wins and ping-ponged the re-emit fingerprint
  every handle); `kt<0` refused (mirrors the kn guard); faceted-surface connectivity must be a whole
  number of segments (`size % nps == 0` — the trailing partial segment was silently DROPPED);
  a missing node or a 2D (`-ndm 2`) node in a contact surface now SKIPS the whole contact LOUDLY at
  handle() (`ladrunoSurfaceNodesOk` — previously silent dead pairs, and a 2D node read `getCrds()(2)`
  OUT OF BOUNDS, unchecked in release).
- **Mortar `inverseIsomap2D` (the PR-1 gate follow-up):** absolute `tolR=1e-13` on a LENGTH-unit
  residual (aux-plane UV ~ facet size h, noise floor ~eps·h) ⇒ GPs silently SKIPped for h ≳ 2000
  (378/400 dead at h=5e4) ⇒ `-mortar` contact and `LadrunoTie -mortar` quietly lost their integration
  at mm-unit-building scales. FIX = the PR-1 parametric-step escape (`|dξ|+|dη| < 1e-10` after the
  update). In-solver gate: `tests/test_contact_review_p3_validation.py` mortar press at h=1 AND h=5e4
  settle to the SAME penalty prediction (pre-fix: free fall at h=5e4).

### Contact-review P4 (2026-07-02): the contact handler silently DROPPED every equationConstraint — LadrunoTie + contact = dead ties
- **Bites:** any model combining `constraints LadrunoContact` with `equationConstraint` — notably the
  fork's own **LadrunoTie (ADR-62)**, which emits ordinary EQ rows. The handler REPLICATES PlainHandler's
  DOF/FE loop (so the contact-FE start tag is knowable), but upstream PlainHandler gained an
  EQ_Constraint block AFTER the replica was written (enforce trivial-identity ⇒ DOF −4; loudly
  warn-and-ignore non-trivial); the replica had NEITHER ⇒ ties ran completely dead with NO diagnostic —
  the exact silent-zero class the handler's non-homogeneous-SP warning guards (review HIGH-3).
- **FIX (upstream parity, no more):** the EQ block is ported — trivial-identity EQs are ENFORCED
  (DOF mark −4, the numberer's EQ code) and non-trivial ones are warned NOT-ENFORCED, pointing at
  `constraints LadrunoProjection`. **RULE: contact + non-trivial ties in ONE analysis remains
  unsupported** (the handlers are mutually exclusive; actual EQ enforcement is the projection
  handler's job) — it now says so instead of silently producing a wrong answer. Gate:
  `tests/test_contact_review_p4_eq_parity.py` (capfd warning assert — FAILS pre-fix on silence).

### Contact-review P5 (2026-07-02): per-contact NormalField sign re-key + the hygiene batch
- **Shared-master `-smoothNormal` sign (MED, the headline):** the ADR-63 nodal-normal field was keyed
  by MASTER SURFACE tag, but its frozen global outward SIGN is a per-CONTACT datum (parity with the
  faceted per-contact `orientDir`). A SECOND contact sharing the master (a two-sided baffle: slaves on
  BOTH faces of one declared plate) inherited the FIRST contact's frozen sign — `signCaptured` skipped
  the whole vote, so even its own explicit `-outward` was IGNORED — and its slaves were read as
  penetrating from the wrong side and EJECTED through the plate (silent pass-through). FIX =
  `theNormalFields` keyed by CONTACT tag (each contact votes + freezes its own sign; the topological
  σ/sharedEdge cache duplicates per contact on a shared master — cheap, fingerprint-cached). Gate:
  `tests/test_contact_review_p5_percontact.py::test_p5_two_sided_baffle_per_contact_sign`.
- **One-time warnings were PROCESS-lifetime:** every contact warning latch was a function-local
  `static bool warned` ⇒ after the first warning anywhere, every later model in the same interpreter
  (wipe/new model) stayed SILENT, and a second contact with the same defect was silenced by the first.
  FIX = per-(contactTag, topic) latches on the engine (`LadrunoContactDomain::warnOnce`, reset
  naturally on wipe since `Domain::clearAll` deletes the engine; kept across `revertToStart` — a reset
  re-run repeats the same configuration). The rigid-plane adapter now carries its contactPlane tag for
  the latch key.
- **`-visc` was NOT statics-inert:** the "v≡0 in statics ⇒ byte-identical" claim ignored that trial
  velocities are STATE — a static stage AFTER a transient one (or `setNodeVel`) keeps the last
  committed velocities, so the dashpot injected a constant unphysical force with NO matching tangent
  (`StaticIntegrator::formEleTangent` never calls `addCtoTang`) and silently shifted the static
  equilibrium. FIX = the dashpot is DISABLED under a `StaticIntegrator` (dynamic_cast gate in
  `LadrunoContactFE::viscousActive`, warn once per contact); genuinely-zero velocities made the old
  force exactly 0.0, so the gate is byte-identical for a pure static run (the shipped D2 static test).
- **Edge-edge first-capture sign could coin-flip:** the E2 body-fixed sign is captured ONCE from
  `sign(n·orientDir)` — with the reference nearly ⟂ n the capture was a numerical coin flip that then
  BOUND the pair for its whole life. FIX = the H2-style conditioning guard (mirrors
  `normalOriented`'s `|proj| < 1e-12·|ref|` refusal) DEFERS the capture (the eval is treated as
  no-contact; retried at the next, better-conditioned config) — and the defer is LOUD (gate lens-A):
  on a symmetric rig the reference can stay ⟂ n forever (e.g. an axis-aligned `-outward` exactly in
  the crossing plane), so a one-time per-contact warning names the pair inert and points at
  `-outward`. Also closes the degenerate-normal hole where `edgeGap`'s failure path left `n`
  uninitialized but the caller kept assembling.
- **`LadrunoEdgeKernel` TAU_LEN was an absolute floor (the PR-1/PR-3 unit-trap class):** `a ≤ 1e-9` on
  a SQUARED length ⇒ every edge shorter than ~3e-5 length units was silently DEGENERATE (edge-edge
  contact dead at small length units). FIX = RELATIVE floor `TAU_LEN_REL = 1e-12` (squared-length
  ratio vs the LONGER edge of the pair; both-zero ⇒ degenerate). Oracle mirrors
  (`proto_e0_closest_point.py` / `proto_e2_penalty.py`) updated in lockstep, all green.
- **Tri/quad shared-edge guard 2× factor:** `onSharedInteriorEdge`'s parametric `edgeTol` (0.1) was
  calibrated on the QUAD span of 2 (5% of span) but applied verbatim to the TRI span of 1 (10%) — the
  P2.1 ownership guard fired twice as aggressively on tri meshes. FIX = tri branch uses `edgeTol/2`
  (same fraction-of-span; quad byte-unchanged — the calibrated P2.1 ridge gates are quad meshes).
- **Coincident/staggered zero-gap mortar AUTO orientation (the PR-3 gate follow-up):** the kernel
  signs the facet normal by flipping it toward `orientDir = scen − mcen`, so the hazard is the
  NORMAL COMPONENT of that vector vanishing — a COINCIDENT conforming pair (`orientDir ≈ 0`) AND a
  zero-gap STAGGERED non-conforming pair (`orientDir` tangential, O(h) — the canonical mortar
  interface; gate lens-B) both leave the flip test to roundoff ⇒ the mortar sign is WINDING LUCK
  (an unlucky winding = attractive contact, no diagnostic). Now WARNS once per contact when
  `|n̂_m·orientDir| ≤ 1e-6·h` at a proximate pair (`|orientDir| ≤ 2h`; far coplanar pairs stay
  silent); `-tie` exempt — its full 3-vec bond is sign-independent and coincident is its common
  case. Pass `-outward` for zero-gap interfaces.

### B3: NTS contact force is NOT in nodeReaction — use the `ladrunoContactForce` query
- **Bites:** reading per-node contact pressure. The NTS contact traction is assembled by an injected
  `LadrunoContactFE` adapter (an FE_Element with no backing Domain Element), so it does NOT contribute to
  `Node::addReactionForce` ⇒ `ops.nodeReaction(slave, 3)` returns 0 for the contact force (only real
  elements + nodal loads/inertia accumulate into reactions). **Fix shipped (B3):** the SEGMENT adapter
  reports its `tn = kn·<−gap>₊` into a Domain snapshot (`set/getNtsForce`, cleared each handle in
  `frictionGCBegin`); query `ladrunoContactForce slaveNodeTag` returns the Σ over the node's pairs. Pure
  side-channel (no resid/tang effect).

### B3: the geometric ∂n/∂u normal tangent is SYMMETRIC but INDEFINITE — a convergence-basin tradeoff
- **Bites:** `contact … -geomtan` on large, soft-master contact patches. `K_geom = kn·gN·H` with the gap
  `gN < 0` in contact, so the geometric block SUBTRACTS from the main `kn·BᵀB` (PSD) ⇒ the contact tangent
  can be indefinite (still symmetric — it's the Hessian of the scalar gap — so ProfileSPD factors it, but
  Newton can leave the convergence basin far from the solution). **Observed:** a 749-slave-node fixed-sphere
  Hertz patch on a soft deformable master DIVERGES with `-geomtan` but CONVERGES without it; a single
  warped-quad slide (few DOFs) converges FASTER with it (quadratic, 4 vs 7 iters). So the geometric tangent
  improves LOCAL (near-solution) convergence but can shrink the GLOBAL basin on big soft patches.
- **Why it's GATED off-default:** exactly this tradeoff. `-geomtan` is an opt-in refinement (like
  `-consistanttan`), default OFF ⇒ the robust `kn·BᵀB` + byte-identity. Turn it on for curved /
  large-sliding interfaces where quadratic Newton is wanted and the patch is well-conditioned; pair with a
  line search / load stepping for robustness on large soft patches. Found by the B3 Hertz study.

### B3: penalty NTS contact bootstrap — a curved indenter starts contact at ONE point
- **Bites:** driving a curved (sphere/cylinder) indenter into a half-space by force. The first contact is a
  single point; the not-yet-contacting indenter material/columns are unsupported ⇒ free-fall under load ⇒
  Newton diverges at step 1 (seen across LoadControl / DisplacementControl / weak-spring-stabilized free
  indenters AND one-shot full pre-penetration of a fixed rigid sphere). The shipped `block_on_block` test
  converges only because its interface is FLAT (all slaves engage at once) at a tiny 1e-8 pre-penetration.
- **What works for a Hertz patch:** a FIXED rigid sphere (slaves pinned at the sphere surface ⇒ no free
  body) with a MODEST approach δ + moderate penalty kn + a looser convergence tol (1e-10 stalls on the
  stiff patch), and/or ramping the indentation gently. Robust quantitative 3D Hertz remains sensitive — it
  motivates pairing B3 with displacement control or D1 within-step augmentation. Found by the B3 gate 2.

### GeneralizedAlpha's tangent is INCONSISTENT with its own residual for αM≠1 (matters for DDM/Newton)
- **Bites:** anything that relies on `GeneralizedAlpha::formEleTangent` being the true ∂R/∂U — DDM
  sensitivity (ADR-52 W3-I2 `LadrunoGeneralizedAlpha`) most sharply, and Newton convergence rate generally.
- **The quirk:** `GeneralizedAlpha::update()` calls `setResponse(*Ualpha, *Ualphadot, *Udotdot)` — it sets
  the model acceleration to the **full-step `Udotdot`**, not the αM-intermediate `Ualphadotdot`. So the
  PRIMAL dynamic residual the elements assemble (`getResistingForceIncInertia`) is
  `R = F(t+αF·dt) − P(Ualpha) − C·Ualphadot − M·Udotdot`, whose consistent Jacobian has **M-coef `c3`**.
  But `formEleTangent` emits `αF·K + αF·c2·C + **αM·c3·M**` — M-coef `αM·c3`. For αM≠1 the tangent ≠
  ∂R/∂U. (Strictly this is also a non-textbook generalized-α: inertia "should" act at `Ualphadotdot`. HHT
  is FINE — it has no αM, M acts at the full step, tangent `c3·M` is consistent.)
- **Why it usually goes unnoticed:** for a LINEAR system Newton still converges to the residual-defined
  fixed point regardless of the (inexact) tangent — just not quadratically. Only when you DIFFERENTIATE the
  converged solution (DDM) does the tangent inconsistency bite: a tangent-reuse DDM gives a biased gradient
  (the W3-I2 PR2 FD-vs-DDM oracle failed ~2e-3 at αM=0.9 before the fix).
- **Fix used by `LadrunoGeneralizedAlpha` (#415):** build the sensitivity RESIDUAL primal-consistent (M at
  `Udotdot`, no αM — like Newmark) AND **re-form the sensitivity-solve tangent with `c3·M`** (a
  `sensTangentFlag` branch in `formEleTangent`/`formNodTangent` + a `formTangent(CURRENT_TANGENT)` in
  `computeSensitivities`) instead of reusing the inconsistent factored primal tangent. Primal path
  untouched ⇒ byte-identical. The K/C terms (αF) need no such fix — K acts at `Ualpha`, C at `Ualphadot`,
  consistent with the tangent. Found by the W3-I2 PR2 adversarial review + the Zone-A FD oracle.

## Collapsing N sibling classes into one without breaking serialization (ADR-52 W1-E2)

**Pattern (reusable):** when you fold a family of integrator/element/material *classTags* into ONE
class selected by flags, the trap is `sendSelf`/`recvSelf` and the object brokers — a saved DB or an
MPI rank reconstructs an object by its **classTag**, via `FEM_ObjectBroker::getNewX(classTag)`, then
calls `recvSelf` on it. If you delete the retired tags or point them at a default-constructed unified
object, the flag state is lost and the reconstructed object behaves wrong (silently — no error).

**What works (ExplicitBathe 6→1, #419):**
1. **Keep all N retired `#define`s** in `classTags.h` (annotate DEPRECATED; never free/reuse — a
   future feature grabbing 33009 would alias a saved model). The collapse frees ZERO tags.
2. **Derive the classTag from the flag combo** in the ctor (`tagForFlags(lnvd,sms,consistent)`),
   passing it to the base `TransientIntegrator(classTag)`. So `getClassTag()` of a `-lnvd -sms` object
   IS the legacy `ExplicitBatheLNVDSMS` tag ⇒ `sendSelf` writes the right tag with no extra field.
3. **Route every retired tag through one static factory** `X::makeForBroker(classTag)` that decodes
   the flags from the tag (`flagsForTag`) and constructs the unified object pre-set. Both brokers
   (`FEM_ObjectBrokerAllClasses.cpp` + `runtime/.../TclPackageClassBroker.cpp`) use a 6-way
   case-fallthrough → one `return X::makeForBroker(classTag);`.
4. **Send a FIXED-SIZE param superset** in `sendSelf`/`recvSelf` (union of all variants' payloads).
   The flags do NOT travel in the payload — they are implied by the classTag the broker already used.
   recvSelf only fills the numerics; it must reset transient bookkeeping (the SMS injected-map,
   energy-registry, warn-once flags).
- **Validity caveat:** the flag↔tag map must be a true bijection over the *valid* combos. ExplicitBathe
  has 6 valid combos (`consistent` implies `sms`), matching the 6 legacy tags exactly — so `tagForFlags`
  defends `if (!sms) consistent=false`. If your flag space has combos with NO legacy tag, you need a NEW
  tag for them (and a manifest row), not a silent collapse.
- **Architecture gotcha that fed this:** the "6 siblings" were really a **2-base × 3-mode lattice**
  (`ExplicitBathe`/`ExplicitBatheLNVD` bases, each ×{none,sms-lumped,sms-consistent}). The two bases'
  `newStep`/`update`/`commit` were **byte-identical**; only LNVD's `formUnbalance` override + the SMS
  `domainChanged` injection differed. Diff the candidate classes' hot paths BEFORE assuming a merge is
  behavior-preserving — here it was, so the byte-identity tests (assert on disp/vel/accel, not stderr)
  hold. Keep each retired OPS_ parser verbatim (exact historical positional grammar) one release; each
  just constructs the unified class with fixed flags.

## `Domain::getElements()` is a SHARED singleton iterator — never iterate the Domain from inside an element callback (ADR-58 P2-S2)

`Domain::getElements()` returns `*theEleIter` after calling `theEleIter->reset()` — **one shared
`SingleDomEleIter`**, not a fresh object (same for `getNodes()` etc.). `Domain::commit()` and
`Domain::update()` walk elements through that single iterator (`while ((e = theEleIter()) != 0)
e->commitState()/update()`). So if an element's `commitState()`/`update()`/`getResistingForce()`
calls **any** Domain method that re-iterates elements — most notably `Domain::calculateNodalReactions()`
(it does `getElements()` + `addResistingForceToNodalReaction` on every element) — the nested
`reset()` **rewinds the iterator the outer loop is using**, so the outer loop terminates early and
**silently skips `commitState()` on every element after the caller**. No crash, no warning — just
some elements never commit (subtly wrong results). This is NOT theoretical: it bit the rigid-body
moment gather, which needed the toe-spring reaction inside `commitState`.
- **Fix pattern:** never iterate the Domain from an element callback. If you need another element's
  force, **cache its tag at `setDomain`** (called from `Domain::addElement`, OUTSIDE any iteration —
  safe) and re-resolve with `Domain::getElement(tag)` + read `getResistingForce()` directly. Cache
  TAGS not `Element*` so a removed element is skipped (`getElement` returns 0), not deref'd after free.
  (`Node`/`Element` do not store incident elements, so there is no per-node shortcut — the setDomain
  scan is the way.)

## Transformation handler is INCREMENTAL (`TRANSF_INCREMENTAL_MP`); an element's `update()` setTrialDisp on a constrained slave SURVIVES (ADR-58 P2-S2)

`TransformationDOF_Group::setNodeDisp` is compiled with `TRANSF_INCREMENTAL_MP` (defined in the
header), so it imposes a constrained (MP-eliminated) node as `slave += T·δu_retained` —
**incremental**, off the retained node's per-step increment, homogeneous (no constant term; `Uc0/Ur0`
are never read on this path). Two consequences exploited for the rigid body's finite-rotation slaving:
1. A linear/time-varying MP can carry only what `T·u_retained` spans — it **cannot** place a slave at
   a nonlinear `(R−I)d⁰` offset whose source (the body orientation `q`) is not a retained DOF. A
   `setTrialDisp` done *inside* a custom MP's `getConstraint()` is **overwritten** by the subsequent
   `incrTrialDisp(T·δu)` (the MP_Joint3D "length-correction" is likewise a slow cross-step nudge only).
2. BUT a `setTrialDisp` done in the **owning Element's `update()`** SURVIVES into residual formation:
   per step the order is `newStep` → `applyLoad`/enforceSPs → `Domain::update()` (element `update()`
   imposes the absolute slave position) → `formUnbalance` (reads it). `setResponse` in the post-solve
   `update(U)` does the handler's incremental write, then `Domain::update()` re-imposes — the element
   always wins the **last write before the next residual**. Verified: slaves track `(R−I)d⁰` to 5e-16.
   This is why the rigid body imposes slaves directly (mechanism "C3") instead of via a custom MP.
   Caveat: the committed slave then lags `q` by one step (the element's `update()` uses the pre-commit
   `qTrial`; `q` advances later in `commitState`) — a standard explicit half-step offset, and the slave
   **velocity/accel are not re-imposed**, so velocity-dependent elements on a slave see an inconsistent
   `v` (keep incident elements displacement-only, or also impose `v_i = v_R + ω×(R·d⁰)`).

## A kinematic mesh-tie (LadrunoTie, ADR-62) inherits ALL of the projection handler's requirements; the GENERATOR must refuse-and-hand-off where they aren't met

`LadrunoTie` emits ordinary `EQ_Constraint`s (`u_s = Σ N_i u_{m,i}`) for the **shipped**
`LadrunoProjectionHandler` (ADR-30) to enforce. So the tie is only as usable as that handler:
`system Diagonal` only, explicit only, no MP-chains / double-constraints, lumped mass on every tied
DOF, ICs on the constraint manifold, partition-interior. The generator front-loads the detectable
violations as **named refusals at model-build** (clearer than a mid-analysis singular solve):
- **node-disjoint + one-facet-per-slave (BLOCKER-1).** A node that is both a slave and a master
  facet node, or a slave listed twice, is refused — those are exactly the MP-chain / double-constraint
  topologies the handler refuses at `handle()`. Collocation (one facet per slave) guarantees disjoint
  slave sets, so no slave appears in two constraints.
- **massed tied DOF (BLOCKER-2).** The projection keeps slave DOFs in the equation set (does NOT
  eliminate them), so a massless tied DOF makes `(LᵀML)` singular. **Subtlety: at the generator's
  model-build time the nodal `mass()` is usually still zero — solid nodes get their mass from element
  `-rho`, ASSEMBLED later.** So the generator's mass check must read `Element::getMass()` (formable
  from rho+geometry at build time, the same call `consistentMassGuard` makes at `handle()`) and mark a
  node "massed" if it has nodal mass OR belongs to an element with a nonzero mass diagonal. Consequence:
  **define `LadrunoTie` AFTER the elements/masses** — emitting it before the mass-bearing elements exist
  trips a false BLOCKER-2 refusal.
- **conforming-at-interface ICs (OQ-3).** The displacement tie is trivially on-manifold for a fresh
  model (all `u=0`), so the IC concern is purely geometric: the slave's reference coords must lie on
  the master surface or the projection/weights are meaningless. The generator refuses a slave whose
  closest-point projection lands farther than `-tol * facet-size` off the surface (default 1e-6).
  v1 does **not** snap ICs (`-projectICs` was declined) — non-conforming means different mesh
  *resolutions* on a *shared* surface, where slave nodes are already on a master facet.

Also: the generator drops near-zero shape weights (`|N_i| < 1e-12`) before emitting, so a slave that
projects onto a facet corner/edge ties only to the master nodes that actually carry a share (a corner
collocation degenerates to a clean `u_s = u_{m,corner}`, i.e. an `equalDOF`), avoiding spurious group
connectivity in the handler's connected-component grouping.

## Integral-mortar ties (LadrunoTie `-mortar`, ADR-62 P2): global D⁻¹ pre-inversion DODGES the "needs handler chain support" wall — but P is dense (one big handler group)

The ADR-62 plan assumed integral-mortar ties would "couple slave nodes ⇒ MP-chains ⇒ need the
projection handler's deferred chain support." **That premise is avoidable.** The mortar constraint
`D u_s = M u_m` (D = slave-interface consistent mass, `D_IJ=∫N_I^s N_J^s dΓ`; M = `∫N_I^s φ_K^m dΓ`)
is condensed **once at model-build** by pre-inverting D over the WHOLE interface: `u_s = P u_m`,
`P = D⁻¹M`. Each row of P then ties a slave to **master nodes only** (no slave appears as a retainer),
so it emits as an ordinary `EQ_Constraint` the SHIPPED `LadrunoProjectionHandler` already accepts
(verified: it allows dense rows + master-only retainers + many slaves sharing masters; only a DOF that
is both retained-master and constrained-slave, or constrained twice, is a "chain"/"double"). **⇒ P2
needed NO handler change and NO kernel change** — it reuses `LadrunoMortarKernel::integratePair` verbatim
and only adds a setup-time generator + one `Matrix::Solve` (DGESV). Confirmed in `proto_p2_mortar_tie.py`
(P·1=1, linear-completeness patch, master-only rows) + `tests/test_ladrunoTie_mortar.py` (genuinely
non-matching solid patch test, 6/6).

CONSEQUENCES / gotchas for the standard-basis condensation:
- **P is DENSE.** `D⁻¹` of a sparse SPD mass is full, so every slave couples to every master in the
  connected interface ⇒ the handler builds ONE large group (all interface DOFs) and factorizes
  `(LᵀML)` over all master DOFs once per `domainChanged`. Fine for typical tie interfaces; for a HUGE
  interface this is the cost a **dual/biorthogonal basis** would remove (diagonal D ⇒ sparse P ⇒ small
  local groups) — that's the deferred P2.1 optimization. (Row-sum LUMPING D is NOT a shortcut: it keeps
  partition-of-unity but BREAKS linear completeness ⇒ fails the constant-stress patch.)
- **DGESV only flags an EXACT zero pivot**, not near-singular/ill-conditioned D. So the generator must
  guard coverage BEFORE the solve: compute each slave node's FULL tributary area `fullCov[I]=∫N_I^s` over
  the whole slave surface via a **self-clip** (`integratePair(npsS,Xs,npsS,Xs,...)` — a facet clipped
  against itself = the full facet, reusing the kernel), then refuse if the master-overlapped `cover[I]` is
  `< (1−1e-3)·fullCov[I]` (the slave protrudes past the master — a partial/extrapolated bond). A SINGLE
  slave facet half-overlapping the master does NOT give any node `cover≈0` (its shape fn spans the
  overlap), so the cover≈0 test alone misses protrusion — the cover/fullCov RATIO is what catches it.
  Belt-and-braces: a **post-solve partition-of-unity check** (`|Σ_k P_Ik − 1| < 1e-6`) catches any
  ill-conditioned solve that slipped through (P·1=1 is algebraically exact, so drift ⇒ bad D).
- **Reference coords, not trial.** Feed `integratePair` the as-built `getCrds()` (NOT X+u): a mesh-tie
  freezes the bond at the reference config (same as P1).
- **refDir (mortar normal orientation)** defaults to the average MASTER facet normal; it only orients n
  (used for the gap g̃ — magnitude only — and the aux plane), so its sign is irrelevant to D/M/P. If the
  master normals cancel (folded/curved surface) the generator refuses and asks for `-outward ox oy oz`.

## `LadrunoContactBucketSort::Grid::runawayGuardFired()` is NOT a clean "node ran away" signal (ADR-60 R8)

The broad-phase grid's runaway guard clamps the centroid bbox to the `[clipPct, 100−clipPct]`
percentiles and sets `guardFired_` whenever that clamp **moves a bound**. With the shipped `clipPct=1.0`
that is the 1/99 percentile, so for ANY mesh with **>100 segment-centroids** the tails are clipped *by
design* and `guardFired_` is true — it does not mean a node diverged. So do **not** auto-surface
`runawayGuardFired()` as a warning (it would fire on every normal large model). ADR-60 R8 deliberately
leaves it a debug-only accessor. The real instability safety on the finite-sliding re-emit deformed feed
is `clipPct=0` (clip disabled so a genuinely-diverging node can't collapse the grid and silently drop
pairs) — a *behavior*, not a warning. Note also: with `clipPct=0` the guard never even computes
(`clip()` early-returns before the `guardFired_` test), so on the re-emit feed it is always false anyway.

## NTS contact slide-off-the-surface is already safe — no explicit detection needed (ADR-60 R5)

A slave that slides clean off the master (out of every segment's parametric domain) needs **no** special
"slide-off" code: `LadrunoContactProjection::evalSegment` gates on penetrating-AND-in-bounds, so an
out-of-bounds slave yields zero force; `project()` returns the out-of-bounds parametric coords
**UNclamped** (no edge-clamp that could hold spurious traction); the migration trigger + friction-slot GC
drop the now-stale adapter within a bounded window; and D4 fresh-slot re-engagement keeps any later
re-pairing traction-continuous. Empirically (`Ladruno_scripts/_probe_r5_slideoff.py` →
`tests/test_adr60_reemit_p4_slideoff.py`): a frictional slave flung off a finite strip's end departs with
force→0, falls freely, and retains its tangential velocity. CAVEAT for diagnostics: the `ladrunoContactForce`
(B3) snapshot is only refreshed on a re-handle (cleared in `frictionGCBegin`), so on the **frozen**
non-`-reemit` path it can report a STALE force after the slave has left contact; `-reemit` clears it each
re-emit so the readout is live there.

### `ShellMITC4` / `ASDShellQ4` `getMass()` NEGLECT rotational inertia (zero on DOFs 4,5,6)
- **Bites:** any code that reasons about a shell node's mass PER DOF — e.g. the ADR-62 P3 shell mesh-tie,
  which keeps tied slave DOFs in the explicit equation set and so needs nonzero mass on **every** tied DOF
  (a massless tied DOF ⇒ singular projection). A shell tie that defaults to tying rotations (DOFs 4,5,6)
  will find them massless.
- **Why:** both stock shells lump **translational** mass only. `ShellMITC4::formInertiaTerms` says verbatim
  "translational mass only // rotational inertia terms are neglected"; `ASDShellQ4::getMass()` says
  "Rotational mass neglected for the moment" and only fills the translational diagonals (`index+q`, q=0..2).
  So `getMass()` returns exactly 0.0 on the rotational diagonals.
- **How to apply:** don't assume a shell's rotational DOFs carry mass. LadrunoTie's BLOCKER-2 was made
  **per-DOF** for exactly this: it names-refuses a tied rotational DOF with no mass and tells the user to add
  nodal rotary mass (`mass $node 0 0 0 mrx mry mrz`) or drop the rotations (`-dof 3 1 2 3`). The per-DOF
  element scan maps global DOF `d` of the node at external position `k` to mass diagonal
  `(Σ preceding-node getNumberDOF()) + d − 1` (standard consecutive-per-node layout). Related:
  [[project_explicit_constraint_projection]] (the handler requires lumped mass on every tied DOF).

### To sparsify a mortar transfer `P=D⁻¹M`, use the DUAL basis — row-sum LUMPING D breaks the patch
- **Bites:** anyone trying to make the mortar slave mass `D` diagonal (so `P` is local/sparse instead of
  dense) — the ADR-62 P2.1 `LadrunoTie -mortar -dual`.
- **Why:** the obvious "diagonalize D" = row-sum LUMPING (`D_II ← Σ_J D_IJ`) keeps partition of unity
  (`ΣP=1`) but is NOT biorthogonal, so it **breaks linear completeness** ⇒ FAILS the constant-stress
  patch (oracle measured 0.28 error vs 2e-15 for dual). The CORRECT diagonalizer is the biorthogonal
  DUAL basis (Wohlmuth): replace the slave TEST functions `N_I` with `ψ_I = Aᵉ·N`, `Aᵉ = diag(cᵉ)(Dᵉ)⁻¹`,
  `cᵉ_a = Σ_b Dᵉ_ab = ∫N_a` ⇒ `Dᵉ_dual = AᵉDᵉ = diag(cᵉ)` EXACTLY (any facet) while preserving the patch.
- **How to apply:** the dual transform is a per-slave-FACET `npsS×npsS` solve applied to that facet's own
  `Mᵉ` — built from the SAME `LadrunoMortarKernel::integratePair` `Dᵉ/Mᵉ` ⇒ **no kernel change, no handler
  change** (same trick as P2's global-D⁻¹ dodge). `P(I,k)=Y(I,k)=(Dᵉ⁻¹Mᵉ)(I,k)`; PoU holds exactly because
  `Dᵉ` and `Mᵉ` integrate over the SAME overlap (`Dᵉ·1 = Mᵉ·1 = cᵉ`). Keep the default (standard dense P)
  byte-identical — the dual path is a separate post-guard block, opt-in `-dual`.

### Mixed-DOF EQ rows (w←(w,θ)) are ALREADY legal end-to-end — but a mortar basis change needs kernel help
- **Bites:** anyone extending a tie/constraint transfer beyond same-DOF rows — the ADR-62 P3.1
  `-hermite` Hermite w–θ shell-edge transfer, and any future shell-to-solid (`θ=½∇×u`) coupling.
- **Why (the good news, verified in source):** `EQ_Constraint` stores per-retained `(node, dof, coef)`
  TRIPLES (`rCoef_i·rDOF_i(rNode_i)`, EQ_Constraint.h:37-45) and `LadrunoProjectionHandler::buildGroups`
  creates one union-find vertex per `(node, dof)` walking arbitrary retained-DOF lists — there is NO
  retained-dof==constrained-dof assumption anywhere. So a slave-translation row that references master
  ROTATION DOFs is an ordinary EQ row: no constraint-class, kernel, or handler change (P3.1 shipped as a
  pure emission-level transform on the collocation path).
- **The catch:** this dodge is COLLOCATION-only for basis changes. `LadrunoMortarKernel::integratePair`
  returns only ACCUMULATED `D/M/g̃` (`PairResult` has no per-GP hook), so putting a different master basis
  (e.g. Hermite functions of `(w, θ_t)`) inside the weak-form `M` integral requires a kernel extension —
  that is why `-hermite -mortar` is a named refusal (mortar-Hermite = deferred P3.1b).
- **Bonus oracle gotcha:** when testing a Hermite w-row against slope-inconsistent (Mindlin) data, do NOT
  sample slaves at edge midpoints — the shear error enters through `H2+H4 ∝ ξ−3ξ²+2ξ³`, which has a root
  at exactly ξ=½, silently hiding the Kirchhoff-assumption error (proto_p3_1_hermite_tie.py T6).

## ADR-63 #4a — averaged nodal-normal smoothing (NTS)

**The auto global-sign vote uses the master-surface centroid over UNIQUE nodes — NOT the flat `mTags`
connectivity.** `LadrunoContactSurface::getNodeTags()` returns the flat per-segment connectivity, in which
edge/ridge nodes shared by K segments appear K times. Averaging that flat list double-counts the
high-valence ridge nodes and biases the "master centroid" toward the ridge; for a `-smoothNormal` master
whose slave cloud sits near that biased centroid plane the auto seed `slave_centroid − master_centroid` can
flip sign → the whole nodal-normal field points INWARD → `gap = n_smooth·(x_s−x̄) > 0` reads "not
penetrating" → silent pass-through (looks exactly like the R3 bug the feature is meant to fix). Caught on the
convex-ridge gate during P1 bring-up (a 90°-ish tent: flat-list centroid z=0.5 ≈ slave z=0.5 → seed_z≈−1e-3 →
G=−1). Fix: dedupe master node tags (a `std::set<int>`) before averaging (`LadrunoContactHandler.cpp`,
ADR-63 field-build block). The auto sign is still fragile when the slave cloud straddles the master centroid
plane — use `-outward` for such masters (the documented escape). The per-segment Newell area-normal vote
(`Σ σ_s·newellAreaNormal(s)`) is NOT affected (it sums over SEGMENTS, each once).

**ADR-63 #4a — sharp-ridge facet ownership (a smoothed-normal limitation).** A slave sitting AT a
sharp convex ridge has its closest-point projection land on the SHARED EDGE of the *adjacent* facet
(barycentric ≈ 0, marginally in-bounds), so that neighbor reads a large SPURIOUS penetration and
injects a big ejecting force. This is a pre-existing NTS facet-ownership issue: the FACETED path only
*incidentally* prunes it (at a 90° ridge the neighbor normal is ~⟂ the per-pair `orientDir`, so
`normalOriented`'s perpendicular fail-safe kills the pair) — the SMOOTHED path (which does not use the
per-facet `orientDir` for the normal) has no such prune. A blunt interior-margin gate in
`evalSegmentSmooth` was tried and REVERTED: it removed a spuriously-*helpful* neighbor force and
destabilized other geometries (the spurious force helps at some ridge angles, hurts at others). Proper
fix = closest-facet / interior facet ownership at a shared edge (a P2 item, related to ADR-57 #4b
edge-handoff). For P1: the R3 SIGN fix is validated by the quad convex-ridge gate; tri-3 chain coverage
uses a FLAT patch to avoid the pathology; a slave pressed onto a facet INTERIOR (away from ridges) is
unaffected.

**ADR-63 P2.1 — the sharp-ridge ownership fix is GAP-AWARE, not a near-edge reject (2026-07-01).** The
above limitation is RESOLVED, but the obvious fix is a TRAP. A per-segment shared-edge mask
(`segmentSharedEdges`, topological, cached with σ) is threaded into `evalSegmentSmooth`; but rejecting
ANY projection that lands within a parametric band of a SHARED edge reproduces EXACTLY the reverted
blunt fix and causes PASS-THROUGH: a frictionless slave slides UP-slope to the ridge apex (the smoothed
normal is more vertical than the facet normal, so a facet-perpendicular load tilts up-ridge), and at the
apex BOTH facets project onto the shared edge ⇒ both rejected ⇒ the slave is driven straight through
(regressed the P1 ridge gate to min_d=−338). The spurious double-activation is LOAD-BEARING at the apex.
The working rule is GAP-AWARE: reject a shared-edge projection ONLY when `−gap > edgeGapFrac·(|g1|+|g2|)`
(`edgeGapFrac=0.05`) — the spurious non-owner reads a penetration ∝ its LATERAL distance from the ridge
(large), while the true owner AND a genuine at-apex contact read a SMALL gap (kept). `edgeGapFrac`/
`edgeTol=0.1` are dimensionless (gap-vs-facet-size / parametric) ⇒ robust across ridge angle + scale.
FREE/boundary edges are never rejected ⇒ R5 slide-off untouched. Residual: near the apex the non-owner
is still kept (a harmless small force ⇒ mild over-stiffness at the exact ridge, NOT pass-through); full
closest-facet selection under sliding = ADR-57 #4b, a P2.3 concern.

**ADR-63 P2.1 — an ill-posed R3 test can pass for the WRONG reason (2026-07-01).** The P1 R3 gate
(`test_p1_smoothnormal_holds_the_ridge_facet`) originally ran the slave LATERALLY FREE and only
"passed" because the P1 spurious ridge ejection (the very bug) pinned `min_d ≥ 0` early. A frictionless
slave under a load with a lateral component on a CONVEX ridge has NO equilibrium — it slides up-slope and
launches — so a laterally-free rig cannot test "held". Once P2.1 removes the ejection the free slave
correctly slides off (min_d≪0). Lesson: to gate a SIGN/normal claim, constrain the confounding DOF —
FIX the slave's x,y (only z free) so the test states the well-posed thing (smooth stays repulsive,
min_d≈−1e-3; faceted flips, min_d≈−318). A "held" assertion over a frictionless free body on a curved
master is suspect.

**ADR-63 build — a deep worktree path overflows the cl.exe command line (2026-07-01).** Building
OpenSeesPy inside `…\.claude\worktrees\<name>\` fails at the include-heavy TUs (OpenSeesPy / SparsePython)
with `CreateProcess failed. The parameter is incorrect.` / `ninja: fatal: CreateProcess` — the long
absolute path, repeated across ~50 `-I` flags, blows the Windows ~32 KB command-line limit. Fix (in
`build.bat` configure): `-DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON` pushes include/object/library lists into
`.rsp` files. Harmless on short paths. NB: adding this flag to an existing build tree triggers a full
recompile (every compile rule's command hash changes).

**ADR-63 P2.2/P2.3 — friction composes with `-smoothNormal` for FREE, but the AUTO sign still needs
`-outward` for curved masters (2026-07-01).** (1) `LadrunoContactFE::segmentActive` builds the friction
slip via `LadrunoFrictionKernel::tangentPart(drel, n, gTvec)` with the SAME `n` the gap operator uses — the
smoothed normal when `useSmoothNormal` — so friction is projected against `n_smooth` with NO new code; on a
flat master (n_smooth==n_facet) a frictional slide is byte-identical smooth-vs-faceted. `-reemit`/
`-smoothNormal`/`-mu` compose (parser refuses only vs `-mortar`), and `-reemit -smoothNormal -mu -outward`
sustains a frictional crossing over a convex curved master (the ADR-60 "exposed combo" closed) — WITHOUT
`-reemit` it still passes through (smoothing fixes the SIGN, re-emit fixes the SEARCH; orthogonal). (2)
TRAP: the ADR-60 R3 `-outward` caveat is lifted only when the global sign vote is WELL-CONDITIONED. A single
slave *starting to the side* of a curved arc votes a near-horizontal seed (slave − master-centroid ≈ ⟂ the
up-field) whose tiny z-component can flip the sign INWARD ⇒ the smoothed field points inward ⇒ pass-through
even with `-smoothNormal` (the pre-existing F2/F3/F5 low-confidence warning fires). So `-smoothNormal` is NOT
a blanket lift of `-outward` for curved masters — it lifts it only for slave clouds sitting OVER the master
(seed ∥ field). Always pass `-outward` for edge-grazing / side-approaching slaves on a curve. **[SUPERSEDED
by P2.5 below — the auto sign is now a robust per-slave majority vote that holds for over-the-master edge/
side-approaching clouds without `-outward`; only a genuinely two-sided or multi-shell cloud still needs it.]**
(3) The P2.1
gap-aware guard's near-apex over-stiffness under SLIDING is a mild quality effect (a small extra bump as a
block crests a SHARP ridge at speed; negligible on realistic shallow arcs, maxpen ~0.01; never diverges),
not a pass-through — full single-owner selection is ADR-57 #4b.

**ADR-63 P2.4 — the frozen-field smoothed tangent CONVERGES implicitly; the dropped `∂n_smooth/∂u` is
`O(kn·gN)` ⇒ sub-dominant on a penalty contact (2026-07-01).** The Q-IMPLICIT-NEWTON tripwire (does Newton
converge with the SUPPRESSED `∂n_smooth/∂u`, i.e. the frozen-field symmetric `kn·BᵀB`, on a genuinely
curved implicit master?) resolved to **outcome (a): converges** — 2 Newton iterations per step, INDEPENDENT
of load-step coarseness (even a single step dragging a slave across the whole facet = maximal within-step
normal rotation still converges in 2). The reason is structural: the dropped consistent-tangent term is
`kn·gN·∂²gN/∂u²`, i.e. scaled by the penalty PENETRATION `gN ≈ press/kn`, which is small in any well-posed
penalty contact ⇒ it never dominates the kept `kn·BᵀB` (this is the SAME reason the shipped faceted default
drops its own B3 block by default). The ONLY regime where smoothed Newton iterations climb (a swept
soft-penalty misuse, `gN ≳ 15%` of the facet) is exactly where the FACETED `-geomtan` consistent tangent
ALSO diverges (and even the seat step fails) — a penalty-NTS-breakdown, NOT a smoothed-normal defect, so
P3's full `∂n_smooth/∂u` would not rescue it. **⇒ P3 (`-consistentNormalSmooth`) stays a genuinely-optional,
evidence-deferred follow-up, not a required item.** Rig gotchas (why the test looks the way it does): (a) a
lone frictionless slave on a convex master has ZERO lateral contact stiffness at the crest (vertical n ⇒
`kn·nx²=0`) and no lateral equilibrium off-crest ⇒ the seat solve is singular/runaway — use a weak lateral
spring (ks≪kn) + seat at the SYMMETRIC centre; (b) a single combined-load `(Fx,0,-P)` under
DisplacementControl is DEGENERATE for a frictionless slave (one load factor scales both the press and the
drive ⇒ equilibrium pins to a single slope) — decouple the constant press via `loadConst` + a separate
lateral drive pattern; (c) DisplacementControl is the ONLY implicit displacement driver for a
`constraints LadrunoContact` model — the handler REFUSES a non-homogeneous (imposed-displacement) SP; (d)
there is NO validated static+`-reemit` path (every ADR-60 reemit test is explicit/CDL), so a FIXED master
(constant nodal-normal field, no re-handle needed) sidesteps it for the implicit rig.

## The `ladruno_opensees.pth` boot module pins ONE worktree's pyd — a fresh build in ANOTHER worktree is silently ignored (ADR-66 P5.1)

`Ladruno_scripts/wire_venv_pth.py` writes `_ladruno_opensees_boot.py` into the py-3.12
site-packages with the generating checkout's `dist\bin` HARD-CODED, `sys.path.insert(0)`-ed at
interpreter startup, and — the sharp edge — an EAGER `import opensees` (for the
`openseespy` aliasing), so `opensees` is already in `sys.modules` before any test bootstrap or
`PYTHONPATH` entry can win. **Symptom:** you build a NEW element in worktree B, the build exits 0,
the pyd timestamp is fresh — and pytest says `element type X is unknown`, because the import came
from worktree A (check `opensees.__file__` FIRST when a freshly-built symbol is "unknown").
**Bypass without touching the other session's wiring:** set `PMI_RANK=1` in the child env (the boot's
MPI guard skips the eager import + aliasing) and `sys.path.insert(0, <your dist\bin>)` +
`os.add_dll_directory` + PATH-prepend in a small driver BEFORE importing pytest
(the P5.1 `run_gates.py` pattern). Re-running `wire_venv_pth.py` re-pins instead, but stomps the
sibling session.

## Solid-shell patch tests on a 1-element-thick mesh: the interior-node patch MUST use a traction-consistent field (ADR-66 P5.1)

Every node of a one-element-thick patch lies ON the free top/bottom faces. A full affine gradient
carries `sigma·e_z != 0` there, so with no applied face tractions the TRUE solution legitimately
deviates from the affine field — a plain-displacement std brick "fails" this exactly like the ANS
element does (~40% at the interior node; replica-verified). This is an ILL-POSED TEST, not element
failure. **Fix:** choose the patch gradient with `eps_13 = eps_23 = 0` and
`eps_33 = −lam(eps_11+eps_22)/(lam+2mu)` (so `sigma·e_z = 0`); the interior node then lands on the
affine field to machine precision (1e-16 in the numpy replica; 1e-6 through the Penalty solve) for
ans and std alike, and the `E33` channel is still exercised (`eps_33 != 0`). Full-traction GP-level
exactness belongs to the FULLY-PRESCRIBED single-element patches. Corollary for reviewers: a
solid-shell "patch test failure" report must state the face-traction handling before it counts.

## `getResistingForceIncInertia` MUST snapshot the shared static `resid` before calling `getRayleighDampingForces()` — else stiffness-proportional Rayleigh silently drops element inertia (ADR-66 P5.1)

An element that builds its residual in a **static class-member** `Vector resid` (the OpenSees
norm) and does `formInternal(); formInertia(); resid += getRayleighDampingForces();` has a hidden
re-entrancy bug: `Element::getRayleighDampingForces()` (Element.cpp:347/349) calls
`this->getTangentStiff()` (when `betaK != 0`) or `this->getInitialStiff()` (when `betaK0 != 0`),
which re-enter the element's own form routine and **`resid.Zero()`** it. The `resid +=` then adds
the damping force to a resid that has been wiped back to `f_int` (or, for the first `betaK0` call
before `Ki` is cached, to **zero**), so the returned unbalance is missing `M·a` — and Newton still
CONVERGES (the Newmark tangent keeps its `c3·M` term) to a wrong dynamics solution, or (as observed
for this element) fails to converge outright. **Fix = the LadrunoBrick donor pattern
(LadrunoBrick.cpp:688-700): `static Vector res(24); res = resid;` BEFORE the Rayleigh call, then
accumulate into `res`.** SYMPTOM: a transient with `rayleigh 0 <betaK> 0 0` gives quantitatively
wrong periods/accelerations (or diverges) while a `betaK=0` run is fine. GATE: compare a `betaK=0`
transient to a tiny-`betaK` one — light damping must change the peak <5%; the bug collapses the
tiny-`betaK` run to quasi-static (or non-convergence). This is INVISIBLE to every static/patch gate
— a dynamic Rayleigh regression is mandatory for any new element with mass. (Verified by
reverting the fix + rebuilding: the gate fails `analyze -3` on the buggy binary, passes on the fixed.)

**2026-07-11 recurrence (ADR-70 P4a adversarial gate): the WHOLE plane family had it.**
`LadrunoQuad`/`LadrunoCST`/`LadrunoLST` under `-geom finite` and `LadrunoCSTPair` all used the
unsnapshotted pattern; their `getTangentStiff()` → `formFinite`/`formPair(1)` refills the shared
static `P` (small-strain paths are safe — those fill only `K`). Extra sting on the pair: the
re-entry also wiped `−Q`, so UniformExcitation ground-motion loads were dropped too. All four fixed
with the donor snapshot in the same PR; parametrized regression gate
`tests/test_ladrunoplane_dynamics.py::test_dynamic_rayleigh_preserves_inertia[pair|cst|quad|lst]`.
LESSON: the quirk was already in this ledger and the P4a author even reasoned about it in a code
comment — and still missed that the clobber happens INSIDE `getRayleighDampingForces()`. When a
quirk names a pattern, grep for the pattern (`P += this->getRayleighDampingForces`), don't reason
about the instance.
**`Vector::pNorm(0)` is NaN-BLIND — never use it for a divergence/NaN check (2026-07-02).** `pNorm(0)`
implements max via `value = (fabs(data) > value) ? fabs(data) : value`; every comparison against NaN is
FALSE, so NaN entries are silently SKIPPED and the returned max is never NaN. The explicit integrators'
circuit breakers (`A_max = U.pNorm(0); if (A_max != A_max || A_max == inf)`) therefore only ever fired on
±Inf: an all-NaN acceleration (NaN material state, poisoned IC, `inf − inf` residual) sailed through and
was COMMITTED into the node state. Fixed by `vectorIsFinite()` (`CriticalTimeStep.{h,cpp}`, std::isfinite
scan) in `CentralDifferenceLadruno`/`ExplicitBathe`/`LadrunoDynamicRelaxation::update()`. Related trap when
WRITING the test: a *seeded* Inf displacement degenerates to NaN before reaching the breaker
(`inf + dt·(−inf) = NaN`), so the honest Inf-path fixture is a genuinely divergent run driven to its first
overflow, and the honest NaN fixture is a NaN committed displacement (`Fint = k·NaN`). See
`tests/test_explicit_nan_breaker.py`.

**betaKinit/betaKcomm MUST be summed with betaK for any explicit stable-step bound (2026-07-02).**
`C = αM·M + βK·K + βK0·K_init + βKc·K_commit` — all three β slots are stiffness-proportional and shrink the
explicit stable step identically (ξ = β·ω/2 grows with ω), and `rayleigh 0 0 βKinit 0` is the MOST COMMON
form in practice (chosen precisely to avoid the current tangent). Reading `getRayleighDampingFactors()(1)`
alone made the SMS damped sizing AND the damped dt_cr estimate blind to it → under-scaled → unstable at
dtTarget for exactly the users the betaK feature targets. Use `stiffnessRayleighBeta(ele)`
(`CriticalTimeStep.{h,cpp}`): per-slot clamp at 0, then sum. Exact at the initial state (K == K0 == Kc);
conservative under softening — the right side to err on for a stability bound.

**The `-divergence` KE proxy FALSE-TRIPS on free vibration at velocity troughs (2026-07-02; FIXED review-P3 #475 -- baseline is now the running MAX of KE).**
The breaker compares per-step `ke/prevKE` against the factor, with `prevKE` updated every step it is
positive. In plain free vibration the velocity passes through ~0 every half period; the step nearest the
zero leaves `prevKE ≈ ε`, and the next steps' quadratic KE regrowth off that near-zero floor produces an
unbounded ratio — phase luck decides whether a given trough exceeds the factor (observed: an SDOF at
dt=0.005, ω=10 tripping factor 10 at its second trough). Workaround in tests/models: excite a
CIRCULAR-motion state (equal springs x+y, quadrature seed) whose total KE is constant, or keep
`-divergence` for monotonic-divergence detection only. Real fix (PR-3 diagnostics batch candidate):
compare against a running MAX of KE, not the previous step.
**ADR-63 P2.5 — the AUTO outward sign is a per-slave MAJORITY vote; a LOCAL closest-point vote beats the
aggregate-normal·global-chord coin-flip (2026-07-01).** The frozen global sign on the auto (no-`-outward`)
path used to be `sign(Σ_a σ_a n_a · (slaveCentroid − masterCentroid))` — an AGGREGATE normal (which nearly
cancels on a curved/domed master) dotted against a GLOBAL chord (which goes ~tangent to the field when the
slave cloud grazes the master edge-on) ⇒ `vote·seed ≈ 0` is a coin-flip and a tiny wrong-signed component
flipped the WHOLE field inward → silent pass-through even with `-smoothNormal` (F2/F3/F5). FIX =
`LadrunoContactProjection::voteSignRobust`: each slave projects onto its NEAREST facet (closest-point,
clamped) to get that slave's LOCAL coherent unit normal `n̂ = σ_{s*}·newell̂(s*)`, then votes
`w = n̂·(slave − surfaceCentroid)` (surfaceCentroid = slot-average of the facet nodes, an INTERIOR
reference for an open convex patch); the surface takes the DISTANCE-WEIGHTED majority `sign(Σ w)`. TWO
adversarial-forced choices: (i) the LOCAL normal (not the aggregate `Σσn`) supplies the lateral
component the aggregate lacked at an edge-grazing slave — the actual F2/F3/F5 fix; (ii) the INTERIOR
CENTROID reference (not the local footpoint separation) keeps a single slave seeded slightly
PENETRATING voting outward — a footpoint separation points inward for such a slave and with one slave
there's no majority to protect it (the P1 sign gate does exactly this; adversarial F2). `|w|`-weighting
then lets a clearly-separated majority dominate. The vote runs on the REFERENCE coords of BOTH master
and slaves (adversarial F1 — the DEFORMED master vs reference slave mix mis-signs on restart / mid-run
recapture; `setNormalField` takes `refSegCoords`, the DEFORMED `segCoords` still drives the per-handle
field). RESIDUAL (LOW): a non-convex open patch whose centroid is not interior ⇒ pass `-outward`.
KEY POINTS: (1) it decides only the
ONE frozen sign (D2/F1) — still captured once, still frozen; (2) `-outward` given ⇒ the aggregate
`sign(vote·outward)` path is UNCHANGED (byte-behavior preserved) — the robust vote runs only on the auto
path; (3) `nVoted==0` (no slave projected) ⇒ fall back to the aggregate seed vote; (4) a genuinely two-sided
cloud yields margin≈0 ⇒ the SAME `conf<0.1` handler warning fires (the ambiguity is DETECTED, recommend
`-outward`); (5) a disconnected multi-shell master is still refused at `propagateOrientation` — a
per-connected-component vote (run `voteSignRobust` per component vs its own nearest slaves) is the
Q-MULTISHELL follow-up; (6) the slave coords fed to the vote are the REFERENCE coords (config-independent —
captured once); (7) RESIDUAL: the degenerate-BLEND fallback still orients by the aggregate seed (review
GAP-2), so a degenerate blend AND an edge-grazing cloud together can still drop a pair (fails safe) — pass
`-outward` for that compound corner. `-smoothNormal` OFF stays byte-identical; no classTag; no vanilla touch.

**An ABSOLUTE tolerance on a DIMENSIONAL residual silently killed all contact away from the origin
(2026-07-02, contact-review fix PR-1).** `LadrunoContactProjection::project()` converged on
`|R| < 1e-12` where `R = d·g_α` has units **length²**: its floating-point noise floor is `~eps·|X|·|g|`,
so for coordinates far from the origin the test was UNREACHABLE — e.g. a plain mm-unit building
(h~500 mm facets at x~5e4 mm, noise ~3e-10) failed **200/200** projections; every pair evaluated
inactive and contact SILENTLY vanished (slave free-falls, no warning). Never caught because every
contact gate ran at unit scale near the origin, where the products happen to be exact in binary.
FIX = a scale-free **parametric-step escape** `|dξ|+|dη| < 1e-8` (parent coords are O(1)) checked AFTER
the Newton update. HONEST behavior contract (the adversarial gate REFUTED a stronger "bit-preserving"
claim): NO previously-converging input is ever LOST (0/5M trials), and on a FLAT facet the escape exits
one iteration early with the IDENTICAL (ξ,η); on a WARPED facet Gauss-Newton contracts only linearly, so
~19% of converging warped cases exit with a footpoint within ~tolP parametric of the residual-converged
answer (measured drift ≤ ~1e-9; gap error second-order ⇒ physically nil; the full 199-test contact+tie
battery, incl. its exact-`==` byte-identity gates, is unchanged) — in-bounds classification can flip ONLY
for a footpoint within the 1e-9 parent-boundary slack (~1e-10·h of slave positions). Oracle:
`contact_prototypes/proto_projection_offset.cpp` (13/21 checks fail on the pre-fix header); in-solver:
`tests/test_contact_projection_offset.py`. THREE general lessons: (1) any tolerance compared against a
quantity with length units must be RELATIVE to a local metric (the `detK` degeneracy guard had the twin
disease — length⁴ vs a length² floor — now the pure angle test `detK < 1e-14·K00·K11`, i.e. sin²θ);
(2) the SAME absolute test is too LOOSE at micro scale: for h≲1e-5 it passes at the INITIAL GUESS
(R ~ |d||g| < 1e-12 before any iteration), so micro-facets get centre-footpoint "convergence" — gap on a
flat facet is footpoint-independent so contact stays functional, just parametrically sloppy (documented,
unchanged; tightening it would break byte-identity for nothing); (3) test meshes at unit scale near the
origin CANNOT catch dimensional-tolerance bugs — put one offset/scaled case in every geometric oracle.
Same review pass: the bucket-grid cap arithmetic used 32-bit `long` (LLP64 Windows) — a diverging
`clipPct=0` re-emit feed with per-axis cell counts in the [~5e4, 2e9] window wrapped the product, exited
the cap loop early, and `grid_.assign()` could request a ruinous allocation (`bad_alloc` kills the
process mid-run). FIX = per-axis pre-clamp to 5000.0 IN DOUBLE before the int cast (also removes the
double→int UB; the total is capped at min(nSeg,5000) cells anyway) + `long long` product arithmetic.
**OPEN FOLLOW-UP (found by the PR-1 adversarial gate, pre-existing, NOT fixed here):** the mortar
back-map `LadrunoMortarKernel::inverseIsomap2D` has the SAME disease — `tolR = 1e-13` ABSOLUTE on a
LENGTH-unit residual (aux-plane UV ~ facet size h, noise floor ~eps·h): measured 0/400 GPs dead at
h≤500 but **241/400 dead at h=5e3 and 378/400 at h=5e4**, and the caller silently SKIPs a failed GP ⇒
`-mortar` contact and `LadrunoTie -mortar` quietly lose most of their integration for facet edges
≳ 2000 length units. Same fix pattern (parametric-step escape or eps-relative tolR) + its own oracle —
file with the review-fix PR-3 hardening batch. Same-theme, likely benign: `isConvex2` tol=1e-12
(length²) and `dedupe` tol=1e-12 (length) in the same header.
**Collapsing a command family into flags? DIFF THE DEFAULTS, not just the grammar (2026-07-02).**
The ADR-52 W1-E2 collapse kept every deprecated alias parsing byte-compatibly, but the UNIFIED command
carried its own `-lump` default (RowSum, upstream-compatible for the bare dt_cr estimate) while the
alias impl defaulted Diagonal ("matches the system Diagonal run") — so `ExplicitBathe p -sms dt` and
`ExplicitBatheSMS p dt` sized the SAME model's scaling with DIFFERENT lumping: 3.73 vs 31.88 injected
mass (8.5x) on a consistent-mass beam, silently. Per-combo byte-identity tests passed because the test
elements had rowsum == diagonal (Truss — the [[project_zonea_link_blocker]] CDL battery caught that
equivalence once before). Fixed (review-P2): `-sms` without an explicit `-lump` flips the sizing
default to Diagonal. LESSON: when merging commands, enumerate every DEFAULT each retired parser had and
prove the merged parser reproduces them per mode — and put a rotational-DOF (rowsum≠diagonal) model in
the byte-identity battery.

## Shell-to-solid plane-section tie (ADR-64 `-shellSolid`): the paid-for gotchas

The P4 rung ties an ndf-6 shell EDGE (master polyline) to an ndf-3 solid FACE (slaves) with
rigid plane-section arm rows `u_s = Σ N_j (u_j + θ_j×d)` (three mixed EQ rows per solid face
node; drilling drops out via the 1e-12 filter). What we paid to learn:

- **Cross-ndf EQ rows are `Transformation`-INCOMPATIBLE** — the handler's factorization goes
  singular (`U(i,i)=0`). Enforce with `Lagrange` (static) or `LadrunoProjection` (explicit)
  only. Documented-unsupported, no investigation planned (ADR-64 OQ-6).
- **Honest kinematic limits, GATED not hidden (ADR-64 OQ-2):** the 3-translation rigid arm
  suppresses through-thickness Poisson stretch at the seam (St-Venant boundary layer, misfit
  = ν·ε·|z| ∝ ν·t, EXACT at ν=0 — oracle T6), and Timoshenko shear warping leaves an O(γ·t)
  local misfit (analytic max γc/(3√3)) that does NOT shrink with in-plane refinement while
  resultant transfer stays exact (oracle T7). Model at ν=0/thin seams or accept the layer.
- **Explicit rotary-mass rule (CORRECTED from the ADR plan):** the generator has NO rotary
  precondition (slaves are solid translations, `-rho` covers BLOCKER-2; static Lagrange
  needs nothing) — but under `LadrunoProjection` the master edge's θx/θy still need a small
  nodal rotary mass (`mass $n 0 0 0 mr mr mr`): the projector keeps every group DOF in the
  explicit equation set, so a zero-mass DOF has no equation of motion (`rigidLink -beam`
  has the identical requirement). The "tied solid supplies the edge's rotary inertia via
  CᵀM_ccC" claim in the ADR draft was measured WRONG for the shipped projector. ALSO: the
  shell element must LUMP mass on the tied DOFs — `ShellMITC4::getMass()` is CONSISTENT
  (off-diagonal) and the handler's element guard refuses it; `ASDShellQ4` lumps
  translational mass and works.
- **`mass` sizes by the BUILDER's ndf, not the node's** (`OPS_addNodalMass` uses
  `OPS_GetNDF()`): in a mixed ndf-3/ndf-6 model, `ops.mass(shellNode, 6 values)` after a
  `model('basic','-ndf',3)` builds a 3×3 and dies with `Node::setMass - incompatible
  matrices` even though the node was created with the per-node `-ndf 6` override. Re-issue
  `ops.model('basic','-ndm',3,'-ndf',6)` before creating/massing the shell side.
- **Collocation force-transfer needs NESTED grids for an exact force patch:** the kinematic
  rows are exact on ANY grid (oracle T1–T5), but the transmitted interface FORCES only match
  the receiving side's consistent nodal pattern when the slave grid nests in the master's
  (2:1 mid-splits etc. — the same property P1's patch test used silently). A non-nested
  split (e.g. 0.4 vs 0.5) redistributes interface forces at the few-% level and fails a
  tight stress patch. Consistent transfer on non-nested grids = the mortar variant
  (deferred: needs the plane-section basis inside the M integral, a kernel per-GP hook).
- **Testing against a stale `dist/bin` binary: EQ_Constraints survive `wipe()`.** The main
  checkout's old pyd predates the ADR-30 `theEQs->clearAll()` fix in `Domain::clearAll`, so
  a second model built in the same process inherits the previous model's EQ rows → duplicate
  (linearly dependent) Lagrange multiplier rows → `U(i,i)=0`. Current source is fixed; the
  trap only bites pre-flight scripts run against an outdated binary — one model per process
  there.
**A `pos < n` bounded staged walk can "pass" its own post-check — pre-walk the totals (2026-07-02).**
The lumped SMS injection walked element nodes with `for (...; pos < n)` and then rejected non-node-major
layouts with `if (pos != n)`. For an element whose nodes' TOTAL ndf EXCEEDS its mass size n, the bounded
loop stops MID-NODE at exactly pos == n, the post-check sees n, and the misaligned mass is silently
committed. The consistent sibling was immune only because its first pass summed ndf UNBOUNDED before
comparing. Same family of trap next door: the consistent M̄ builder indexed `mdiag[base[a]+d]` for
d < min-ndm without clamping by each node's OWN ndf — an ndf < ndm node spills the write into the next
node's block (off the end of Mbar/mdiag on the last node). Fixed (review-P3): pre-walk Σndf and reject
≠ n up front; clamp `ndmOf[a] = min(ndm_a, ndf_a)`. LESSON: a walk over a matrix whose layout a DIFFERENT
object (the nodes) defines must validate the FULL mapping before mutating anything — loop bounds are not
validation. Related: an iterative-solver "did not converge" warning must key on the RESIDUAL, not on
iters == maxIt — the consistent PCG's pAp≤0 SPD-breakdown guard exits EARLY (iters < maxIt, resid > tol)
and the old `iters >= maxIt` condition swallowed it silently (also fixed review-P3).

### Isotropic sqrt(A) `getCharacteristicLength` regularizes a mesh-objectivity band ONLY on in-plane-SQUARE elements (dx == dy)
- **Bites:** any crack-band/energy mesh-objectivity study (or user model) with `LadrunoSolidShell` (and any element whose lch is the isotropic sqrt of the midsurface/element area) meshed with in-plane rectangles. The crack band localizes in ONE element column, so the physical band width is the element size ALONG the band normal (dx) — but the material regularizes with lch = sqrt(dx*dy). For dx != dy the dissipated energy is off by sqrt(dy/dx), and a "refinement" that changes the aspect ratio reads as spurious energy drift even with `-autoRegularization` working perfectly.
- **Why:** the scalar lch has no direction; sqrt(A) == dx only when dx == dy. The through-thickness projection is already excluded by design (ADR 66 D6), but the in-plane anisotropy is not.
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** keep localization-band meshes in-plane square (the G5 gate enforces dx == dy at all three densities: spread 0.9% across a 4x size range, vs the fixed-lch control at ~4x energy error). A directional `lch(n)` API is the ADR 66 O4 backlog item (shared with ADR 19's sqrt(2)-strut residual).

### A hairline-weakened localization band (ft_band ~ 0.95 ft) drives BULK GPs into the Concrete3D return-map apex — seed bands with >= 20% strength margin
- **Bites:** weakened-band localization tests/models with `LadrunoConcrete3D`: with the band only 5% weaker, the bulk sits at ~95% of its own tensile onset at band peak, and the global Newton's trial excursions through the localization transition push bulk GPs past onset into the deep-tension apex regime — the run drowns in "return map did not converge -> step-cut" (the kernel's safe fallback), grinding to a crawl without ever failing outright.
- **Why:** Newton iterates are not monotone: mid-iteration trial strains overshoot the converged state by far more than the 5% margin; the Concrete3D apex regime is exactly where the return map is trajectory-fragile (documented handoff §6 gap).
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** seed localization with ft_band = 0.8*ft (dissipation is Gf-governed, so energy gates are unaffected by the seed strength); the churn disappears entirely.

### Mass scaling multiplies the support-motion START SHOCK by sqrt(s) — a Linear sp ramp that is harmless unscaled can crack elements at the moving support under SMS
- **Bites:** quasi-static explicit runs driven by prescribed support motion (`sp` + timeSeries, the #333 recipe) under `CentralDifferenceSMS`. A `Linear` series applies a velocity STEP v at t=0; the wave it launches carries sigma ~ rho'*c'*v = sqrt(rho'*E)*v — and mass scaling inflates rho' by s = (dtTarget/dt_e)^2, so the shock stress grows by sqrt(s) = the dtTarget factor. Measured (ADR 66 G9c): a rate whose unscaled shock is a trivial 0.35 MPa hit 3.5 MPa > ft at the 10x target and cracked the pulled-face concrete element outright (omega_t -> 0.97) before any real loading happened.
- **Why:** impedance rho*c = sqrt(rho*E); uniform scaling multiplies rho by s and leaves E alone. The shock rides the SCALED impedance while the "quasi-static" rate was budgeted against the unscaled one.
- **Workaround/status (2026-07-06, ADR 66 G9):** drive support motion with a SMOOTHSTEP displacement protocol (lam_end * u^2(3-2u), zero start/end velocity — a ~60-point Path series), never a raw Linear ramp, whenever mass scaling is on; budget the ramp against the SCALED wave speed c' = c/sqrt(s). Same discipline applies to any velocity IC under SMS.

### A consistent-mass element under the explicit recipe (`system Diagonal`) runs 8/27 mass-deficient — and the CDL `-lump` flag lumps only the PENCIL, not the runtime
- **Bites:** any explicit run of an element whose `getMass()` returns the CONSISTENT matrix (LadrunoSolidShell pre-G9; any vanilla brick/solid without a lumping option). `system Diagonal` keeps just the raw diagonal — for a trilinear solid that is `rho*int N_a^2 dV` = 1/27 of the mass per node, 8/27 total — so the dynamics run 3.375x light (frequencies 1.84x high) with NO warning. Worse, `integrator CentralDifferenceLadruno ... -lump hrz` does NOT fix it: the `-lump` flag selects the model for the `criticalTimeStep()` PENCIL only (`CentralDifferenceLadruno.cpp:338/475` pass it to `computeCriticalTimeStep`); the runtime mass is whatever `getMass()` feeds the SOE. The combination is the trap: the hrz pencil reports the full-mass thickness CFL `t/c` while the run evolves on the light diagonal — the step is over-reported by sqrt(27/8) = 1.84x and the run BLOWS UP at 0.9x the reported dt_cr (measured on the P5.1 solid-shell: true boundary ~0.54x = exactly the consistent-eig value).
- **Why:** lumping is an ELEMENT property in this architecture (LadrunoBrick `-lumped`, truss default-diagonal); the integrator flag exists so the pencil can MATCH whatever the element does, not to impose it. The CDL default `-lump diagonal` (diagonal-of-consistent) is self-consistent with a consistent-mass element under `system Diagonal` — but then the physics is the 8/27-deficient one.
- **Workaround/status (2026-07-06, ADR 66 G9):** give the element a row-sum/HRZ `-lumped` flag (LadrunoSolidShell has one as of G9; LadrunoBrick always had it) and pair it with `-lump hrz`/`rowsum` on the integrator — then pencil == runtime == full mass (solid-shell gates: period == 4L/c to 5%, pencil == t/c, bisection-validated). RULE: under `system Diagonal`, every element in the model should provide a genuinely lumped `getMass()`.

### Multi-element Concrete3D softening under implicit Newton cut-crawls even where single elements pass — use `-implex` (uniform LoadControl) for band/localization runs
- **Bites:** implicit displacement-ramped runs of a MULTI-element `LadrunoConcrete3D` specimen through localization: plain Newton + step-cutting (the recipe that walks a SINGLE element through its limit point) converges only at micro-steps once several elements carry the indefinite softening tangent simultaneously — the G5 coarse band (4 elements!) took ~2400 micro-steps / 13 min wall; refinement makes it worse.
- **Why:** the indefinite/non-symmetric tangent cluster around the band throws the global Newton into cut/recover oscillation; no single step fails permanently, so nothing surfaces except wall time.
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** put the MATERIALS in `-implex` and drive with CONSTANT-dlam LoadControl (the uniform-pseudo-time regime IMPL-EX wants; kinematic sp ramp): the SPD-ish secant lets Newton track the full localization at the planned step size (13 min -> 16 s on the G5 coarse mesh; committed states stay implicit-exact). The ADR 66 risk register lists exactly this toolbox row; the `-implex`+DisplacementControl limit-point trap (Concrete3D ledger) does NOT apply because LoadControl dlam is uniform.
**`LysmerTriangle` stage-3 `getResistingForce()` MUTATES state on every call — any recorder that
reads element forces perturbs it (2026-07-05, ADR-69).** At stage 3 ("preserve elastic spring
forces after gravity") `getResistingForce()` executes `internalForces -= springForces` on EVERY
invocation — it is not idempotent. The EnergyBalanceRecorder (v1 AND v2) calls it once per record
per element, so each record subtracts `springForces` again from the member the residual path also
serves (rebuilt only at the next `getResistingForceIncInertia`). Consequences: (a) stage-3 Lysmer
energy readings are untrustworthy; (b) anything else querying element forces between residual
formations (nodal reactions, other recorders) compounds it. The ADR-69 leak publisher deliberately
RECOMPUTES `R_inj = getDamp()*v_gnd` in `commitState` instead of reading the member, so E_inject is
immune. Upstream-origin behavior — left unfixed (vanilla change budget); avoid stage 3 + per-step
force recorders in the same model, or accept the drift.
**Tcl `eleLoad` SILENTLY ACCEPTS unknown `-type` flags (returns TCL_OK, no warning) — a
no-op that looks like success (2026-07-06, ADR-69 P0.5).** The eleLoad handler's tail falls
through to `return 0` when no `-type` branch matches, so `eleLoad -type -fooBarBazNotALoad`
"succeeds". Any deck relying on a loader that is not actually wired (LysmerVelocityLoader was
exactly this for 15 years) runs unloaded with zero diagnostics. Left unfixed (changing the
return could break decks); when a load seems dead, FIRST verify the `-type` string exists in
TclModelBuilder.cpp before debugging the physics.

**Stage-0 `LysmerTriangle` under implicit Newmark realizes only ~HALF the dashpot energy
(DW = 0.50*ULW measured; 2026-07-06, ADR-69 P0.5 F2).** The element's
`getResistingForceIncInertia` uses `0*v_node + gnd_velocity` — the node-velocity damping
force C*v NEVER enters the element residual; under Newmark damping then acts only through the
a1*C term in the effective tangent, giving an energy-inconsistent solve (the recorder's
DW = int v'Cv dt books the full ideal dashpot power and RES exposes the ~0.5*W gap). EXPLICIT
integrators assemble the damping force from getDamp() directly and are consistent. For
implicit absorbing runs prefer ASDAbsorbingBoundary; for Lysmer prefer explicit. NOT a
recorder bug — the recorder is the instrument that surfaced it.

**`UniformExcitation` input work HIDES INSIDE the EnergyBalance recorder's IE column (IE = -DW
exactly, RES accidentally closed; 2026-07-06, ADR-69 P1.6).** Elements that implement
`addInertiaLoadToUnbalance` (FourNodeQuad etc.) store `-M*ug''` in their element load vector
`Q`, and `getResistingForce` returns `K*u - Q` — so the recorder's IE integral
(`int F.v dt`) silently accumulates MINUS the seismic input work. The balance then "closes"
with IE the negative mirror of the genuine absorbed/damped energy and ULW = 0. Consequence
for closure gates: NEVER drive an energy-balance validation model with UniformExcitation —
the input-work pollution drowns whatever leak the gate is trying to isolate (use initial
velocities: no patterns, Q = 0, IE = pure strain energy). Not a recorder bug per se — a
consequence of OpenSees folding element loads into the resisting force.

**`ASDAbsorbingBoundary2D/3D::addInertiaLoadToUnbalance` is a deliberate NO-OP ("we don't
need this!") — free-field columns are NEVER driven by `UniformExcitation` (2026-07-06,
ADR-69 P1.6).** The FF masses (`addMff`) receive no `-M*ug''` effective load, so under
uniform excitation the FF column rides rigidly in relative coordinates: zero strain, zero
`addRffToSoil` transfer, dead lateral boundary. The intended input path for ASD boundaries
is the BOTTOM compliant base (`"B" -fx/-fy` time series); lateral elements take no
time-series args at all (parser rejects them for non-bottom). Also note the lateral
`addClk` dashpot writes only SOIL rows (one-way coupling): the FF column is UNDAMPED unless
element Rayleigh (`addCff`, alphaM) is set — an undamped FF column rings forever.

**openseesmp flat-per-rank nodal-term summation is convention-dependent — the split-mass idiom
sums CORRECTLY, only full-mirror emits double-count (2026-07-06, ADR-69 P2, measured).** The
upstream MPI example declares `mass(4, m, m)` for a shared node on BOTH ranks; the parallel
diagonal assembly SUMS duplicate contributions, so the assembled system has `2m` and each
rank's EnergyBalance recorder books only its own share — the naive cross-rank sum of the nodal
columns equals the serial (2m) reference EXACTLY (gate `energy_v2/p2_mpi_owned_nodes.py`, G2).
The "nodal terms multiply-counted on shared boundary nodes" hazard (ADR-69) applies only to
FULL-MIRROR conventions: PartitionedDomain-style external-node mirrors, or an emitter writing
the full nodal mass/load on every touching rank (which also changes the assembled physics
unless the solver dedups). Consequence: don't "fix" per-rank energy sums blindly — first
determine which convention the model uses; `-ownedNodes <regionTag>` is the dedup tool for
mirror conventions and a no-op burden otherwise. Also note per-rank output files
(`stem.part-<rank>.ext`) are auto-suffixed under a detected MPI launcher since ADR-69 P2 —
ranks racing a single recorder file was the previous (corrupting) behavior.

**Modal-damping energy is published only by integrators using the BASE
`IncrementalIntegrator::commit()` (2026-07-06, ADR-69 P2).** 35 integrators override
`commit()` without chaining (HHT family, `*_TP` explicit): they still APPLY modal forces (via
`TransientIntegrator::formUnbalance` / their own `addModalDampingForce` calls) but never reach
the publish site, so their modal dissipation stays in RES and no `E_modal` column appears
(declare-on-first-publish prevents a silent zero column). Newmark does not override commit and
is fully covered. If you need E_modal under HHT: either chain the override to the base commit
(vanilla edit, ledger it) or accept the documented RES drift.

**Recorders NEVER receive domainChanged() — any recorder caching pointers is a
use-after-free waiting for `remove element`/`remove node` (2026-07-06, ADR-69 P2.1).**
`Domain::removeElement` calls `domainChange()` (sets a flag) but `Domain::record()` invokes
recorders directly with no invalidation hook, and the analysis propagates domain changes only
to its own components (handler/numberer/integrator/algorithm). Worse, a recorder CANNOT poll
`Domain::hasDomainChanged()` — it is STATEFUL (consumes the flag, increments currentGeoTag,
resets graph-built flags) and belongs to the Analysis; calling it from a recorder would eat
the analysis's own change detection. `getDomainChangeFlag()` is a pure read but is usually
already consumed by the time record() runs. The working patterns (both in-tree): (1)
re-resolve entities BY TAG on every emit (LadrunoMonitorRecorder, #489); (2) structural
re-validation per record — compare `getNumElements()/getNumNodes()` (O(1)) against cached
sentinels + verify cached tag→pointer identity before any virtual call through a cached
object, rebuild on mismatch (EnergyBalanceRecorder P2.1). Key ALL membership maps by TAG,
never by pointer — a freed pointer can be REUSED by a new allocation and silently inherit the
old binning. Note the count sentinel alone is blind to remove-then-readd-same-tag (counts
restore); the tag→pointer identity check is what catches it.

### `-implex` looks INERT under ASDShellQ4 in fully-prescribed rigs — it is not: ASDShellQ4 reports the POST-COMMIT state, and IMPL-EX re-integrates implicitly at commit

- **Bites:** A/B-ing `-implex` on/off (LadrunoRCConcrete) with a fully prescribed
  (every-DOF `sp`) single-element bending rig under `ASDShellQ4` + `LayeredShell`:
  the recorded responses are BIT-IDENTICAL (rel 4.7e-16 over 60 softening steps) —
  it looks like the flag is dropped somewhere in the section copy chain. It isn't:
  the same rig under `ShellMITC4` (same section) or `LadrunoSolidShell` (3D view)
  shows the expected ~3.5% IMPL-EX extrapolation difference.
- **Why:** LadrunoRCConcrete's IMPL-EX `commitState` re-integrates implicitly to
  advance the TRUE thresholds (the ASDConcrete3D pattern), so the post-commit
  material state is the implicit one. ASDShellQ4's reported element forces reflect
  the post-commit section state; in a rig with NO free DOFs the extrapolated trial
  stresses are never consulted by any equilibrium iteration, so the recorded curve
  collapses exactly onto the implicit run. ShellMITC4 / LadrunoSolidShell report
  the converged TRIAL state, where the extrapolation lives.
- **Workaround/status (measured 2026-07-07, ADR-66 G7):** not a defect on either
  side, but three consequences: (1) never "verify implex engaged" with a prescribed
  probe under ASDShellQ4 — it is structurally invisible there; (2) cross-element
  parity benchmarks must run implicit-vs-implicit (as G7 does) or accept a
  reporting-path asymmetry masquerading as element deviation (~2% here); (3) on
  free-DOF problems implex under ASDShellQ4 IS active (0.7% path shift measured on
  the G7b rig) — the wall-harness usage is fine. Pinned discriminatingly by
  `tests/test_ladrunoSolidShell_flexure.py::test_implex_reporting_paths`.

### `section LayeredShell` bending is exact only on its own midpoint rule — a predictable Σ h³/12 stiffness deficit vs the continuum (≈2% at 5 uniform layers)

- **Bites:** elastic cylindrical bending of a LayeredShell with n uniform layers
  undershoots `E·t³/12(1−ν²)` by `Σ E_i·h_i³/12` (each layer is ONE fiber at its
  centroid: the midpoint rule loses the layer's self-inertia). Compared against
  LadrunoSolidShell — whose `-nz` gauss/lobatto rule integrates z² EXACTLY — this
  reads as "the solid-shell is too stiff". It is the layered quadrature, on both
  the concrete AND the displaced-by-rebar bookkeeping.
- **Workaround/status (measured 2026-07-07, ADR-66 G7):** predict it (the G7
  elastic anchor asserts the layered arm to 1e-4 against the midpoint-rule closed
  form, deficit 2.07% at 5 core layers) or halve the layer thickness (error ∝ h²;
  0.5% at the G7 production layering of 3+10+3).

### `equationConstraint` refuses zero coefficients ("WARNING invalid rcoef inputs") — skip the zero-arm terms when emitting plane-section rows

- **Bites:** emitting `u_i − u0 − θ·z_i = 0` rows for a node ON the reference
  plane (z_i = 0) fails parse: the EQ parser hard-rejects `coef == 0.0`.
- **Workaround (2026-07-07, ADR-66 G7):** drop the θ term when `|z_i| < tol`
  (the row degenerates to `u_i = u0`, which is exactly right) — the same filter
  the LadrunoTie shell-solid tests use for near-zero shape weights.
## `wipe()` does NOT recreate the Domain — new domain-level state MUST be reset in `Domain::clearAll()` (ADR 46 P1)

`ops.wipe()` calls `Domain::clearAll()` on the SAME Domain object; nothing is
reconstructed. Any new domain-level member you add therefore leaks across model
generations unless you reset it in `clearAll()` yourself (the ADR-30 `theEQs` and
ADR-39 contact-engine cleanups are the same lesson). ADR 46 P1 hit this twice in
one PR: (1) the new Rayleigh-factor domain copy survived wipe → the NEXT model
reported phantom damping; (2) the UPSTREAM `theEigenvalues`/`theEigenvalueSetTime`
also survive wipe — latent forever because nothing read the spectrum across a wipe
until `complexEigen` did (a stale spectrum from the previous model silently answers
for the new one). Both now reset in `clearAll()` (`// Ladruno ADR46`).

Related trap: `Domain::getEigenvalues()` **exit(-1)s** (kernel-killer, MPCO class)
when the spectrum was never set — any graceful-failure caller must probe via the
additive `Domain::getNumEigenvalues()` (ADR 46) BEFORE calling it.

Related trap: `region ... -rayleigh` and per-element `-rayleigh` write factors
straight onto elements via `MeshRegion::setRayleighDampingFactors` / element
`setRayleighDampingFactors` and NEVER touch `Domain::setRayleighDampingFactors` —
the domain-level copy reflects only the last GLOBAL `rayleigh` call. `complexEigen`
scans elements and warns on mismatch (D1 policy: warn, never silently absorb).

ADR 46 P2 pre-warning (from the P1 Opus gate): once Route B projects a full
`ΦᵀCΦ` and emits `ψ=Φz`, REPEATED eigenvalues need the basis M-orthonormalized
WITHIN each repeated eigenspace (ARPACK does not guarantee it); the P1 closed form
is immune (uses eigenvalues only, never Φ). **RESOLVED BY DESIGN at P2:** the
assembler projects the FULL `M̃=ΦᵀMΦ` (never assumes `I`) and synthesizes
`K̃=sym(M̃·diag(ω²))` (exact: `KΦ=MΦΛ` column-wise, and `M̃` commutes with `Λ`
within a repeated eigenspace), so the QZ pencil is exact in ANY normalization and
for repeated eigenvalues — no re-orthonormalization step exists to get wrong.

Related trap (ADR 46 P2, OpenSees classic): `getDamp()`/`getMass()` on Elements
AND Nodes habitually return the SAME per-class/per-size static scratch matrix —
never hold references to both at once; deep-copy the first before calling the
second (`const Matrix Ce(el->getDamp());` then consume `el->getMass()` directly).

Route-B validation of the `-doRayleigh` contract (ADR 46 P2 test discovery):
**`Truss` also defaults `-doRayleigh` to 0** (the [[project_damping_channels]]
"OFF for the rest" family, not just zeroLength) — a default Truss feels NO global
`rayleigh` betaK in a transient, and the assembled complexEigen (correctly!)
reports zeta=0 for such a model while `-closedForm` claims bK·w/2. When a test or
model expects stiffness-Rayleigh on trusses, build them with `-doRayleigh 1`. The
divergence between the two routes on default elements is the P1-vs-P2 contract
difference, not a bug: Route B = the C the transient feels; closed form = the
global factors as if every element carried them.
is immune (uses eigenvalues only, never Φ).

### `rigidLink beam` across an ndf mismatch (6-DOF master, 3-DOF slave) WARNS and silently adds NOTHING — the model solves DISCONNECTED

- **Bites:** the natural first attempt at a shell↔solid seam — `rigidLink('beam',
  shellNode, solidNode)` per matched node column — prints
  `RigidBeam - mismatch in numDOF between constrained Node ... and Retained node ...`
  to opserr and returns WITHOUT adding any constraint and WITHOUT raising. The
  analysis then runs to completion with the two meshes disconnected: under a
  moment the shell rides freely as a hinge and the solid stays unloaded — easy
  to read as "the connection is just flexible".
- **Why:** vanilla `RigidBeam` builds a square ndf×ndf constraint matrix and
  hard-requires matching DOF counts; the failure path is a warning + early
  return, and the interpreter command does not convert it into an error.
- **Workaround/status (measured 2026-07-07, ADR-66 G10):** the seam recipe is
  `LadrunoTie -shellSolid` (exact plane-section rows, works on non-matched
  grids, dt_cr-neutral under the explicit projector). The disconnect trap is
  pinned by `tests/test_ladrunoSolidShell_seam.py::test_rigidlink_cross_ndf_silently_disconnects`
  so a future rigidLink behavior change surfaces loudly.

### `gmsh.initialize()` REPLACES the process's native Win32 PATH — later child processes fail DLL/exe resolution battery-wide

- **Bites:** any pytest battery (or long-lived python session) that runs gmsh
  and LATER spawns a child by bare name or runs a freshly-built MinGW exe.
  Symptoms observed (2026-07-07): `test_hrz_standalone_kernel` exe dies at load
  with `0xC0000135` (STATUS_DLL_NOT_FOUND — libstdc++ not on the inherited
  PATH) and `test_cpp_kernel_matches_oracle_dump` raises
  `FileNotFoundError: WinError 2` on bare `"g++"`. Both PASS in isolation and
  fail in the full battery — the classic order-dependence smell.
- **Why (proven with a ctypes probe):** `gmsh.initialize()` calls
  `SetEnvironmentVariable("PATH", <system dirs + C:\ + sys.prefix>)` — a ~315
  char stub replacing the real ~2 kB PATH — and `finalize()` does NOT restore
  it. Python's `os.environ` is a STARTUP SNAPSHOT, so python-side code (and
  `shutil.which`) still sees the good PATH; only native consumers break:
  `CreateProcess` bare-name lookup and the DLL loader of child processes both
  read the LIVE native block. Import alone is harmless — the nuke fires at
  `initialize()`. The worst variant was a zone_b module running a module-level
  gmsh mesh: pytest COLLECTION imports every module (deselection happens
  after), so the corruption hit even `-m zone_a` batteries that never run a
  single gmsh test.
- **Fixes shipped:** (1) never call gmsh at module level in a test file (the
  zoneb mesh is now lazy via `_mesh_once()`); (2) `tests/conftest.py` autouse
  `_native_path_resync` re-syncs the native PATH from `os.environ` before
  every test on Windows; (3) the two g++ checker tests are hermetic anyway —
  absolute `shutil.which("g++")` path, `env=os.environ.copy()` on every
  `subprocess.run`, exe outputs under `tmp_path`.
- **General rule:** in any long-lived process that touched gmsh, pass
  `env=os.environ.copy()` and absolute executable paths to `subprocess`.

### `MumpsParallelSolver` accepts an `mpi_comm` constructor argument and silently ignores it — always binds to `MPI_COMM_WORLD`

- **Bites:** ADR 43 P3 (FEAST-over-sub-communicator, D2 spike, 2026-07-07). Anyone
  assuming `MumpsParallelSOE`/`MumpsParallelSolver` can be pointed at an
  `MPI_Comm_split` sub-communicator today by passing a comm handle — it compiles and
  runs, but silently uses `MPI_COMM_WORLD` for the factorization anyway (wrong ranks
  participate, or it deadlocks/hangs on a partial-world sub-comm).
- **Why:** `MumpsParallelSolver::MumpsParallelSolver(int mpi_comm, int ICNTL7, int
  ICNTL14)` (`MumpsParallelSolver.cpp:54-64`) takes the parameter but never stores it —
  no member is set. `initializeMumps()` (`:93-105`) hardcodes
  `id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD)` on the Intel-MPI path (this fork's
  Windows/oneAPI build), `0` (MUMPS's own WORLD) under `_OPENMPI`; the rank/size probe
  two lines later also reads `MPI_COMM_WORLD` directly. Dead parameter, not a config
  toggle.
- **Fix (ADR 43 P3a, not yet built):** store the passed communicator on the solver,
  `MPI_Comm_c2f` *it* (not WORLD) into `id.comm_fortran`, and use it for the
  `MPI_Comm_rank/size` probe too; thread it through `MumpsParallelSOE`'s constructor,
  which today never passes one either.

### A fully-prescribed model (zero free DOFs) under `constraints Transformation` FATALLY exits the process — `FullGenLinSOE::getX - vectX == 0`

- **Bites:** any prescribed-displacement rig that pins/`sp()`s EVERY DOF of a small
  patch model (single-element material-response probes are the classic case: all x
  fixed for uniaxial strain, bottom y fixed, top y driven by `ops.sp` in a pattern).
  Under `constraints Transformation` the sp-handled DOFs are condensed OUT, the
  equation count hits 0, and `FullGenLinSOE::getX()` hits a raw
  `opserr << "FATAL ..."; exit()` — killing the whole Python kernel/pytest run with
  no traceback (surfaced 2026-07-07 while building `tests/test_planestrain_sigma_zz.py`,
  PR #525). Other SOEs have sibling zero-size exits; this is not FullGeneral-specific.
- **Why:** Transformation removes constrained DOFs from the numbered system; a model
  where every DOF is fixed or sp-prescribed leaves size-0 vectors that the SOE treats
  as an allocation failure, and OpenSees's error path is `exit`, not a recoverable
  analysis error.
- **Fix:** for fully-prescribed rigs use `constraints("Penalty", 1e15, 1e15)` — the
  prescribed DOFs then STAY in the system (size > 0) and the penalty violation at
  1e15 vs typical stiffness is ~1e-10 relative, invisible to material-response
  checks. Alternatively leave at least one genuinely free DOF in the model.
- **Fix (SHIPPED, ADR 43 P3a):** `setCommunicator(MPI_Comm)` on the solver stores the
  comm, `MPI_Comm_c2f`s *it* (not WORLD) into `id.comm_fortran`, uses it for the
  rank/size probe, tears down any live MUMPS instance, and clears the SOE's
  `factored` flag. Test hook: `system('Mumps', '-commSplit', color)` (collective —
  every rank must call it). Gate: `feast_d2_spike/p3a_commsplit_gate.py` (4 ranks,
  2 concurrent groups vs serial oracles). **Residual subtlety:** `MPI_Channel`
  hardcodes WORLD/tag-0, so only the MUMPS factor/solve is comm-isolated — the SOE's
  B/X exchange rides WORLD envelopes, safe today via disjoint (src,dst) pairs +
  MPI non-overtaking with phase ordering; true envelope isolation is a P3c-MPI item.

### MKL `dfeast_srci` SHRINKS its in/out `m0` argument in place — don't feed the shrunken value back into your own subspace bookkeeping

- **Bites:** any FEAST RCI driver that passes its subspace-size variable by pointer
  into `dfeast_srci` and afterwards compares the found count `m` against that same
  variable (saturation / auto-enlarge logic). Surfaced building ADR 43 P3c-serial:
  a band holding 5 modes with a seeded `m0 = 15` came back with `m0` REWRITTEN to 5,
  so `m == m0` looked like subspace saturation, the enlargement re-ran from the
  shrunken value, plateaued (`1.5*5+8 → shrunk to 5 → again`), and after
  `maxEnlarge` retries the solver REFUSED a perfectly complete result.
- **Why:** the RCI's `m0` is in/out session state — when FEAST's reduced Rayleigh–Ritz
  detects the band holds fewer modes than the subspace, it adapts the working
  subspace DOWN and records that in `m0`. That is convergence bookkeeping, not a
  saturation verdict. (The packed driver `dfeast_scsrgv` does NOT exhibit this —
  its `m0` comes back unchanged, which is why the P1/P2 battery never saw it.)
- **Fix:** give the RCI session its OWN `m0` copy (`int m0Rci = m0;` …
  `dfeast_srci(..., &m0Rci, ...)`) and keep the caller's requested `m0` for the
  `m >= m0` saturation compare. Post-shrink, the live right-hand sides at
  `ijob=11` are the first `m0Rci` columns of `workc` — solve with `m0Rci` RHS,
  not the original `m0`. See `runFeastRci` in `FeastEigenSolver.cpp`.

### A REPLICATED MKL FEAST loop across MPI ranks DEADLOCKS unless MKL is pinned single-threaded

- **Bites:** any design that runs an MKL routine (FEAST `dfeast_srci`/`dfeast_scsrgv`,
  or any threaded LAPACK) REPLICATED on every MPI rank in lockstep, with
  collective work between the rank-local MKL calls. This is exactly the ADR 43
  P3c-MPI (L3-only) model: every rank runs the same `dfeast_srci` outer loop and
  they cooperate on ONE distributed MUMPS inner solve at each `ijob=10/11`.
- **Why:** multi-threaded MKL LAPACK is **not bitwise-reproducible across
  processes** (dynamic work scheduling, thread count, CPU affinity all vary the
  reduction order → last-ULP differences). Two failure modes, both silent until
  they aren't: (1) the reduced m0×m0 eig inside `dfeast_srci` rounds differently
  per rank → the refinement-loop **count** diverges → one rank returns `ijob=0`
  (done) while another issues another collective inner solve → **hard deadlock**;
  (2) the stochastic auto-seed (`dfeast_scsrgv` with `fpm[13]=2`, taken when
  `-m0` is unset) rounds `mEst` differently per rank → different `m0` → different
  `nrhs` at the first solve → `MPI_Bcast` length mismatch → MPI abort. MKL often
  runs small dense solves single-threaded (below an internal size threshold), so
  it **usually works and intermittently hangs** — the worst profile. Caught by
  the P3c-MPI adversarial gate (the differential test masked it by pinning
  `MKL_NUM_THREADS=1` in its env; the shipped code enforced nothing).
- **Fix (what actually shipped, after a false start):** the obvious fix —
  `mkl_set_num_threads_local(1)` around the replicated span — **SEGFAULTS at
  runtime in this MP MKL DLL layout** (c0000005 on the first call; the global
  `mkl_set_num_threads` too). Bisected: guard disabled ⇒ distributed solve runs
  clean; guard enabled ⇒ crash. The symbol links (import lib has it) but the
  call dies — never exercised in the serial build because the guard is inactive
  there. So do NOT reach for MKL thread control here. Instead, attack the
  **catastrophic** desync directly: broadcast the stochastic auto-seed `m0` from
  rank 0 (`LadrunoFeastInnerSolve::agreeInt` → `MPI_Bcast`, serial default =
  identity) right after the estimate, so the FIRST distributed solve's block
  width (hence its `MPI_Bcast` length) is identical on every rank. The *other*
  desync source — the reduced m0×m0 eig inside `dfeast_srci` rounding
  differently under threading — does **not** bite in practice: m0 is small
  (tens), below MKL's LAPACK threading threshold, so it runs single-threaded and
  stays bit-identical; the enlargement decisions key off the collective solve's
  `m`/`info` (identical on all ranks post-broadcast). **PROVEN**: the P3c-MPI
  gate passes at `MKL_NUM_THREADS=4` (multi-threaded) with only the broadcast,
  no thread pin. For pathologically large bands (m0 in the hundreds, where MKL
  might thread the reduced eig) set `MKL_NUM_THREADS=1` — documented, not
  enforced (the enforcement API is the thing that crashes). See `agreeInt` +
  the m0-broadcast in `FeastEigenSolver::solve()`.

### MUMPS SYM=2 wants the LOWER triangle; PARDISO mtype −2 wants the UPPER

- When porting the block-real `(zM−K)` solve from serial PARDISO
  (`LadrunoBlockZKernel`, mtype −2, **upper** triangle CSR) to distributed MUMPS
  (`LadrunoDistBlockZKernel`, SYM=2), the stored triangle FLIPS: MUMPS symmetric
  expects entries with global **row ≥ col** (lower), matching how OpenSees's own
  `MumpsSOE`/`MumpsParallelSOE` assemble symmetric matrices (they store
  `row > vertexTag`). For the 2n block `[[aM−K,−bM],[−bM,−(aM−K)]]` the lower
  triangle is: `(i,j) j≤i` = aM−K; `(n+i, j)` for **all** j = the −bM block
  (row n+i ≥ n > j, wholly lower, supplied in full, its transpose block NOT
  supplied); `(n+i, n+j) j≤i` = −(aM−K). Verified transpose-consistent by the
  P3c-MPI adversarial gate.

### ADR-44 LadrunoModalResponse — modal-transient (P1a) + FRF/SSD (P2) gotchas

- **Element-level stiffness-proportional (`betaK`) Rayleigh damping ≠ assembled
  `a1·K` (Truss).** A direct `Newmark` run with `rayleigh 0 0 0 a1` (or `betaKcurr`
  / `betaKinit` — all three identical) on a linear Truss chain differs from the
  EXACT classical-modal solution of `M ü + a1 K u̇ + K u = −MR ü_g` by **several
  percent, dt-INVARIANT** (measured 4.4% on a 2-DOF chain; does NOT shrink as
  dt→0, so it is NOT Newmark truncation). The exact modal answer was pinned three
  ways (numpy modal-superposition 1e-18, numpy full-matrix Newmark of (M,a1K,K)
  ~9e-5, and `modalResponseHistory` itself). Mass-proportional (`alphaM`) Rayleigh
  DOES match OpenSees-Newmark to truncation (~1e-4). **Consequence for validation:**
  never oracle an exact modal-damping feature against an OpenSees direct-`Newmark`
  run that uses `betaK` — use `alphaM`-only, `modalDamping`, or an explicit
  full-matrix Newmark reference. (Root cause: element `getRayleighDampingForces`
  builds the damping force from a per-element stiffness that is not the same as the
  globally-assembled `K` for the truss basic-system mapping — not chased further.)

- **`Path` timeSeries `getFactor(t)` returns 0 at exactly the record end.** For a
  transient whose last station time `t = nsteps·dt` coincides with the final
  abscissa, floating-point round-off of `nsteps·dt` lands just past the last
  sample ⇒ `getFactor` returns 0, not the last value. Symptom: the FINAL committed
  station is slightly off (all earlier stations exact). Fix in models/tests: pad
  the record with ≥1 trailing sample so the analysis window is strictly interior.
  (This is inherent to any getFactor-sampling integrator, incl. UniformExcitation.)

- **`eigen` on a tiny model needs `-fullGenLapack`.** The default ARPACK path
  requires `NEV < N` (NCV must be > NEV and ≤ N); a 1- or 2-DOF model with
  `eigen 1`/`eigen 2` fails `_saupd info = -3`. Use `eigen -fullGenLapack N`
  (dense) for small verification models.

- **Recorder text precision.** Node recorders default to ~6 significant figures;
  a modal-transient vs analytic assert at 1e-8 is dominated by that rounding.
  Add `-precision 15`/`16` to the recorder (already noted for ADR-46; re-bitten here).

- **`responseSpectrumAnalysis -scale` is a dead knob (stock).** The Petracca
  `ResponseSpectrumAnalysis::solveMode` computes `u = V*Vscale*MPF*Sa/λ` and never
  multiplies by `m_scale`, so the `-scale` factor the command parses is silently
  ignored. Not our bug (upstream); the ADR-44 P1b `-combine` path matches this
  behavior (no scale) for consistency rather than silently diverging. Flagged
  during the P1b wiring.

- **FRF sign convention (P2) is `e^{+iΩt}` — pin it or the phase silently flips.**
  `frequencyResponse`/`steadyStateDynamics` use `H_a(Ω)=1/(ω_a²−Ω²+iΩd_a)` with a
  `+iΩd_a` imaginary part, i.e. the `e^{+iΩt}` time convention → the response LAGS
  90° at resonance (`u=+i/(ωd)` for the mass-normalized SDOF, `angle=+90°` because
  the base-accel modal load carries the extra `−Γ`). The opposite (`e^{−iΩt}`,
  `−iΩd_a`) magnitude is IDENTICAL and only the phase sign differs, so a flipped
  convention passes every |·|/RMS test and is invisible until someone reads phase —
  the classic frequency-domain bug. Pinned two ways: the resonance-phase assert in
  `test_sdof_frf_closed_form`, and the end-to-end match to the direct complex solve
  `(K−Ω²M+iΩC)⁻¹(−MR)` (same convention on both sides) in `test_mdof_frf_vs_direct`
  and `modal_response_p2_spike/frf_oracle.py`. Magnitude gates alone would NOT catch
  a sign flip.

- **P2 frequencies are Hz, not rad/s.** `-freq fmin fmax nf` bounds and the output
  `f` column are Hz; internally `Ω=2πf`. So the FRF magnitude peaks at `f_a=ω_a/2π`,
  not at `ω_a`. (A `-radps` option was deliberately NOT added in v1 to keep the knob
  count down; documented in the guide.) Reusing the P1a `LadrunoModalDamping.coeff`
  gives `d_a=2ξ_aω_a` with `ω_a` in rad/s, consistent with the internal `Ω`.

- **P2 reuses P1a's normalization by construction — do NOT re-derive it.** The FRF
  recovery weight is `ψ_a(node,dof)·(−Γ_a)` with `ψ_a=eigenvector·Vscale` and `Γ_a`
  the `modalProperties` participation factor — byte-for-byte the P1a per-mode
  recovery ingredients, so the frequency-domain result is the analytic steady state
  of the exact same modal ODE P1a integrates in time. This is WHY the numpy oracle
  can use clean mass-normalized modes (`m_a=1`, `ψ=φ`, `Γ=φᵀMR`) while the C++ uses
  OpenSees' Vscale/Γ: P1a already proved the two normalizations agree physically
  (it matched direct Newmark), so P2 inherits that instead of re-validating it.

- **An UNDAMPED mode sampled exactly at `Ω=ω_a` gives an infinite FRF.** With
  `d_a=0` the denominator `ω_a²−Ω²+iΩd_a` is real and hits 0 at resonance → `inf`.
  Harmless with any real damping (imag part ≠ 0). **CORRECTION (ADR-44 P2 review):**
  the `-biased` grid clusters points in a ±5% window around each in-band `ω_a/2π`
  AND its `k=0` sample lands **EXACTLY ON** `ω_a/2π` (`f = fa + halfw·(k/NCLUST)`,
  k=0 ⇒ f=fa — `LadrunoModalResponse.cpp:632-634`). So a `-biased` sweep with
  **zero** damping DOES manufacture this singularity (inf/NaN row at every in-band
  mode) — the parser now emits a one-line WARNING when `-biased` is combined with an
  exactly-zero `d_a`. A `-lin` grid landing exactly on an undamped resonance would
  do the same. The same singularity appears
  at `Ω=0` (the `f=0` sample of a `fmin=0` sweep) if a RIGID mode is retained
  (`ω_a=0` → denom `0+0i`) — an unrestrained base-excited structure genuinely has an
  undefined DC steady displacement. Both are honest physical singularities (NaN/inf
  row), documented not guarded; fully-supported base-excitation models never carry a
  rigid mode.

## `OPS_GetStringFromAll` buffer contract differs between interpreters (classic-Tcl reads garbage)

**Symptom** (ADR-71 P1, 2026-07-10): `element LadrunoUP ...` in classic Tcl fails
with "got 0 leading integers" while the identical command works in openseespy.

**Cause**: the two backends implement `OPS_GetStringFromAll(buf, len)` with
DIFFERENT buffer semantics. The modern interpreter
(`OpenSeesCommands.cpp:1185`) copies the token into `buf` and returns it; the
classic-Tcl shim (`elementAPI_TCL.cpp:494`) just does `return OPS_GetString();`
— **it never touches `buf`**. Any parser that calls the function and then reads
`buf` (a common family idiom) parses garbage under classic Tcl and silently
"works" everywhere else.

**Rule**: always consume the RETURN VALUE, never the buffer. If you need the
token in a local buffer, normalize:
`const char* s = OPS_GetStringFromAll(buf,len); if (s && s != buf) strncpy(...)`
(see `upGetTok()` in `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp`). Fix the vanilla
shim only via an ADR — several upstream parsers may depend on the current
behavior.

## Symmetric-storage solvers silently DROP one u-p coupling block on LadrunoUP models — ProfileSPD (the no-`system`-command DEFAULT) returns rc = 0 with garbage pore pressures

**Symptom** (ADR-71 P1, 2026-07-11): a LadrunoUP transient run without a
`system` command (or with `system ProfileSPD` / any symmetric-storage SOE)
"succeeds" — `analyze()` returns 0 every step, displacements look plausible —
but the pore-pressure field is wrong by ~87 orders of magnitude. Measured on
the identical Terzaghi 1×10 Q4 column script: `system UmfPack` gives
p ∈ [0, 5] kPa (matches the Terzaghi series); `system ProfileSPD` gives
p ~ 1e88–1e89 with every `analyze()` still returning 0. Nothing fails loudly.

**Cause**: the honest-p contract (ADR-71 §3.2 ⟨FW-F1⟩) makes the effective
transient tangent UNSYMMETRIC — −Q lives in `getTangentStiff()` (u-rows) and
+Qᵀ in `getDamp()` (p-rows), so c₁K + c₂C never has matching off-diagonal
pairs. Symmetric-profile assembly stores only upper-triangle-in-profile
entries: one Q block is silently discarded at `addA()` time, and the solver
then factors and solves the mutilated (still well-conditioned-looking) system
cleanly. No framework hook lets an element reject an SOE.

**Rule**: every LadrunoUP model MUST name a general solver — `system UmfPack`
/ `SuperLU` / `FullGeneral` / `BandGeneral` (serial), MUMPS with `SYM=0` in
the MPI targets. The parser prints a one-line notice at element creation; the
Zone-B battery pins the divergence
(`tests/test_ladruno_up_element_analytic.py::test_wrong_solver_divergence_profilespd_vs_umfpack`).

### `ops.printA('-ret')` returns the raw OpenSees Matrix buffer, which is COLUMN-major — read row-major, an unsymmetric tangent looks transposed and FD-vs-tangent tests false-fail at exactly the asymmetry magnitude
- **Bites:** any test/tool that reshapes `printA('-ret')` into `(neq, neq)` C-order and compares against an oracle or FD residual. On symmetric tangents the bug is invisible; on LadrunoUP's unsymmetric [K,−Q;0,H] the −Q block appears in the p-row/u-col slot and the check fails at |Q|/|K| (measured 2.6e-5 — small enough to chase as a "tolerance problem" for hours).
- **Why:** `OpenSeesCommands.cpp:2590` hands back `&A(0,0)` flat; OpenSees `Matrix` storage is column-major (Fortran order).
- **Workaround:** reshape Fortran-order or transpose after reshape (`np.array(ret).reshape(neq, neq, order='F')`). Pinned in `tests/test_ladruno_up_element_equiv.py` helper (ADR-71 P1, 2026-07-11).

### Staged NONZERO pressure `sp` added mid-analysis under `constraints('Transformation')` converges cleanly to a WRONG steady state on LadrunoUP (interior p wildly off; Penalty/Lagrange correct)
- **Bites:** the ADR-71 §3.2 initialization recipe — run a stage, then add `sp p=<head>` and continue. Under Transformation the model converges (rc=0) to interior p ≈ −73 for a top head of +1 on a sealed 4-element column (measured); the identical `sp` present from step 1 is handled correctly by all three handlers, and Penalty/Lagrange are correct in both sequences.
- **Why (suspected):** Transformation condenses constrained DOFs at analysis-setup time; a mid-analysis `sp` after `wipeAnalysis` re-setup interacts with the committed-but-unconstrained p state; exact mechanism not chased (P4 revisit alongside the gravity/hydrostatic init recipe).
- **Rule for the guide:** staged prescribed-head sequences use `constraints('Penalty', ...)` (or Lagrange); mirrors the existing fully-prescribed-rig Transformation trap. Repro pinned in `tests/test_ladruno_up_element_analytic.py` (ADR-71 P1, 2026-07-11).

### LadrunoUP `-dynSeepage on` (the default) DIVERGES under Δt-refinement in quasi-static consolidation runs — the ü-term feeds integrator noise into the seepage source
- **Bites:** ZS84-class consolidation column, Newmark γ=0.6: with `-dynSeepage off` the error converges 1.3e-3 → 7e-4 as Δt shrinks 0.08 → 0.005; with the default `on` it GROWS 1.8e-2 → 8.7e-1. Smaller Δt is WORSE: trial accelerations of numerically-damped compressible-wave modes are noise, and f_seep integrates them.
- **Why:** the dynamic-seepage drive (b − ü) is physically right for genuine dynamics (B5-class, P4-gated) but quasi-static consolidation has no meaningful ü — the term is pure noise amplification there.
- **Rule (AMENDED at P4):** the default is now **`off`** — the B5 Simon gate measured the failure in genuine dynamics too (wandering post-front p ≈ 1.7–2.0 vs β = 0.973; unbounded shallow-station growth). `-dynSeepage on` is an explicit research opt-in (ADR §12 log 2026-07-13). Companion: `-stab auto` adds ~10% spurious ringing on wave problems — wave runs use `-stab off`. Measured in `tests/test_ladruno_up_element_analytic.py` sweep (ADR-71 P1, 2026-07-11).

## T6 quirks: constant-mean B-bar/F-bar is RANK-DEFICIENT; nodal-lumped corner masses are ZERO (ADR 70 P3)

Two hard-won T6 (quadratic triangle) facts from building `LadrunoLST`:

1. **Constant element-mean dilatation (B-bar / centroid-sampled F-bar) loses
   rank on the T6.** The two quadratic CONFORMAL displacement modes (Re/Im of
   z²: u=(x²−y², −2xy)-type) have identically zero deviatoric strain, and
   their linear dilatation has zero element mean (mean over the 3-interior-
   point rule = centroid value for a linear integrand). Any constant-mean
   averaging therefore assigns them ZERO energy: a free element shows **5**
   zero-energy modes (3 RBM + 2 spurious; stacked-B̄ rank 7 of the required 9).
   The Q4 is immune only because z² fields are not in its bilinear space —
   "J varies so F-bar can average" is NOT sufficient for rank. Caught by the
   locked single-free-element `assert_zero_energy` T0 gate; confirmed by numpy
   rank and compiled `eigen`. Cure DECIDED by the ADR-70 P4 spike (2026-07-11):
   disjoint 2-triangle patch macro-element — patch-constant J̄ for T3 (dSNPO
   §15.1.9), patch-P1 projected dilatation for T6 (see the [[#Volumetric-
   projection traps on triangles (ADR 70 P4 spike)|P4 quirk entry]]). Until it
   ships, triangles are `std` and near-incompressible plane problems use the
   quad `bbar`.

2. **SixNodeTri-style plain N-lumping gives EXACTLY ZERO corner masses on the
   T6** — ∫N_corner over the 3-interior-point rule integrates to 0, and goes
   NEGATIVE on distorted elements (adversarial gate verified both). Unusable
   for explicit dynamics (diagonal M⁻¹). `LadrunoLST::getMass` therefore
   deliberately DIVERGES from upstream and uses **HRZ lumping** (∫ρN²_a dV
   rescaled to the exact total — strictly positive), the same call ADR-72 made
   for the H20 (row-sum corners were −M/8 there). The reduce-to-`SixNodeTri`
   gate is static, so the anchor is unaffected. Documented in `LadrunoLST.h`;
   gated by a free-free finite-eigenvalue test.

Also: upstream `SixNodeTri` registers in the shared functionMap as **`tri6n`**
(not "SixNodeTri"), and its `setPressureLoadAtNodes` carries a side-61
`dx61 = x4-x6` copy-paste typo that breaks closed-contour equilibrium of the
consistent pressure load — `LadrunoLST` fixes it (`x1-x6`) and pins the fix
with a zero-net-resultant test.
## `randomResponse` PSD convention is ONE-SIDED in Hz — the factor bugs have unmistakable signatures (ADR 44 P3)

The `-inputPSD` time series is sampled at **f in Hz** and read as the **one-sided**
PSD `G(f)` of the base acceleration (`σ_üg² = ∫₀^∞ G df` — the wind / floor-vibration
/ equipment-spec convention). Against the random-vibration-textbook **two-sided
rad/s** PSD `S(Ω)`: `G(f) = 4π·S(Ω=2πf)`, and the white-noise SDOF anchor becomes
`σ_x² = G0/(8ξω³)` (NOT the textbook `πS0/(2ξω³)`). If a future edit scrambles the
convention, the Monte-Carlo gate reads it immediately: a one-sided/two-sided mixup
shows as a ~41 % (√2) RMS error, an Hz/rad mixup as ~150 % (√2π) — both pinned in
`modal_response_p3_spike/psd_rms_oracle.py` (0.6 % agreement when correct). Related
trap in the same spike: a synthetic realization `Σ√(2G·df)·cos(2πf_k t+φ_k)` has
EXACT variance only over a full period `T = 1/df` — validate over whole periods or
the input-variance check itself wobbles.

## Staleness guards cannot see a stale `DomainModalProperties` that reproduces the SAME spectrum — write guard tests with unique stiffnesses (ADR 44 P3)

`DomainModalProperties` survives `wipe()` (the [[#`wipe()` does NOT recreate the
Domain — new domain-level state MUST be reset in `Domain::clearAll()` (ADR 46 P1)|
clearAll leak]] family). The P1a/P2/P3 staleness guards compare eigenvalue count +
element-wise values between the Domain and the snapshot — so a
`wipe(); rebuild-IDENTICAL-model; eigen` sequence leaves a stale-but-equal snapshot
the guard legitimately CANNOT distinguish (same spectrum ⇒ same Γ/Vscale up to sign
⇒ numerically the same answer, so it is also harmless). The trap is in TESTS: a
`guard_no_modalproperties` pytest that rebuilds the same `m,k` as any earlier test
in the file will NOT raise. Give guard-test models a stiffness unique within the
file (`test_ladrunoRandomResponse.py` uses k=512 for exactly this reason).

## Volumetric-projection traps on triangles (ADR 70 P4 spike)

Three traps from the P4 design spike (all numpy-pinned + adversarially
twin-verified in `Ladruno_implementation/adr70_p4_spike/`):

1. **A quadrature L2 projection is the IDENTITY whenever #GP = dim(projection
   space) at unisolvent points.** On the T6's 3-interior-point rule, projecting
   the dilatation onto P1 (3 modes, 3 samples) interpolates — b̃ ≡ b at the GPs
   to ~4e-15, straight OR curved (weights and detJ cancel algebraically). An
   "element-local P1-projected dilatation" formulation flag would silently ship
   `std`. Same trap one level up: patch-P2 over a 2-T6 patch (6 GPs = dim P2)
   is also the identity. Volumetric relief on the triangle REQUIRES coupling
   beyond the element's own quadrature — there is no element-local escape.

2. **Enlarging the constant-mean averaging region can never fix the T6.** The
   P3 conformal modes re-center: u = a(z−z_p)² about the REGION centroid z_p
   has ε_dev ≡ 0 pointwise and zero region-mean dilatation, for any region. So
   dSNPO §15.1.9 F-bar-Patch (one J̄ per patch) cures the T3 pair (exactly 3
   RBM) but leaves the T6 pair with 5 zero-energy modes. The T6 cure must see
   the LINEAR variation of the dilatation (patch-P1 does; rank exactly 3,
   robust distorted + conforming-curved, inf-sup β_h ≈ 0.46 plateau on straight
   structured pair-meshes — unstructured/curved inf-sup unproven).

3. **F-bar-Patch cross-element tangent blocks are STRESS-PROPORTIONAL** (dSNPO
   eqs. 15.37–15.38): K^(es) vanishes at F = I / zero prestress, where the
   macro tangent is symmetric. An FD consistent-tangent gate run at zero
   stress CANNOT catch a wrong/missing/symmetrized cross block — the gate must
   run at finite stress, and the pair element must declare an unsymmetric
   tangent. Also: a shared edge whose mid-nodes don't exactly coincide silently
   cracks the patch and re-opens spurious modes (bit the adversarial twin's own
   curved test — assert conformity).

### `InitialStateAnalysis on|off` (openseespy/interpreter path) heap-corrupts the process (0xc0000374) on the NEXT model operation — UPSTREAM dangling-parameter bug, any element — **FIXED by backport (see below)**
- **Bites:** `ops.InitialStateAnalysis("off")` appears to work (prints its notice, `revertToStart` runs, displacements AND honest-p pressures zero correctly) — then the next `ops.wipe()` / model build crashes with a Windows heap-corruption fault. Found while gating the ADR-71 init-sequencing recipe (P4).
- **Why (upstream, verified in source):** `OPS_InitialStateAnalysis` (SRC/interpreter/OpenSeesMiscCommands.cpp:1352-1354 and :1367-1369) does `theDomain->addParameter(theP); delete theP;` — but `Domain::addParameter` STORES the pointer in the parameter container (Domain.cpp:897). The container now holds freed memory; the wipe-time parameter cleanup double-frees. Not element-specific and not a Ladruno defect — the MSVC/ucrt heap checker is just the first to notice.
- **FIXED (backport of upstream `191c67c2d`, [#568](https://github.com/nmorabowen/OpenSees/pull/568), 2026-07-13):** all three handler copies (interpreter, classic Tcl, unbuilt `runtime/`) now use a stack `InitialStateParameter` + `setDomain` — nothing is registered in the domain, so wipe/exit are clean AND repeat calls actually toggle the flag. NOTE the naive fix (deleting the `delete theP;` lines) would have been WRONG: the fixed tag 0 means the second call hits `Domain::addParameter`'s already-exists early return *before* `setDomain`, so `off`-after-`on` never flipped `ops_InitialStateAnalysis` — post-ISA `ops.reset()` then silently kept PM4Sand-family committed state (`revertToStart` no-ops while the flag is stuck true) and *compounded* `InitialStateAnalysisWrapper`'s ε₀ (`mEpsilon_o += mStrain` fires again → ~3× initial stresses after one reset). Regression: `tests/test_initial_state_analysis_lifetime.py`. Historical workaround (`ops.reset()` for the zeroing step, pinned in `tests/test_ladruno_up_init_recorders.py`) remains valid and is still the right recipe for honest-p sequencing — but is only SAFE with materials that don't consume `ops_InitialStateAnalysis`; pre-fix builds must not combine ISA + reset with the UW soil materials.
## ASDConcrete3D's HardeningLawStorage is a process-global store-if-absent registry keyed by MATERIAL TAG — it survives `ops.wipe()` (ADR 72 P1)

`ASDConcrete3DMaterial.cpp::HardeningLawStorage::store` (static singleton;
`if (item == nullptr) item = make_shared(hl)`) latches the FIRST hardening law
ever constructed for a given material tag, for the life of the process;
`recover(tag, type)` hands that original law to every later material with the
same tag. `ops.wipe()` does not clear it. Consequence in a shared pytest
process (Zone-A runs the whole battery in ONE process): a test file that
builds `ASDConcrete3D` with a common tag (e.g. 1) and *different parameters*
than another file silently poisons the latter — the ADR-72 P1 battery's
advisory test (tag 1, toy Gf) drove `test_ladrunoBrick_asdconcrete.py`'s
mesh-objectivity dissipation ratio from ≈4 to 7.03 with zero diagnostic,
alphabetical test order deciding the victim. **Rule: every test file gives its
ASDConcrete3D materials a file-unique tag** (ADR-72 uses 337218). Same class
of static-registry risk: `CrackPlanesStorage` (same file). Upstream ASDEA
code — fix-on-touch only.

## `printModel('-JSON')` without `-file` prints NOTHING (upstream interpreter quirk, ADR 72 P2)

`OPS_printModel` (SRC/interpreter/OpenSeesCommands.cpp:2927-2930) treats
`-JSON` as flag-only: the argument loop sets `flag = OPS_PRINT_PRINTMODEL_JSON`
and moves on, but `theDomain->Print(...)` is reached ONLY inside the
filename branch (`print <file>` / `-file <file>`). Bare `ops.printModel()`
prints the whole domain (classic format, opserr); bare
`ops.printModel("-JSON")` silently emits nothing at all — a capfd-based
"assert JSON contains X" test then fails with empty capture (bit the ADR-72
P2 S8 acceptance test). **Rule: to assert on model JSON, always use
`ops.printModel("-JSON", "-file", tmpfile)` and read the file** (OVERWRITE
mode, so the file is complete + fresh). Classic-format assertions via capfd
remain fine. Upstream behavior — fix-on-touch only.

### FourNodeQuad-family element `rho` is NOT serialized and is uninitialized by the broker ctor → garbage mass on `database`/`restore` (and OpenSeesMP) → non-deterministic, indefinite M on transient restart — **UPSTREAM bug, FIXED (see [[LEDGER_vanilla_files]])**

The continuum plane elements `FourNodeQuad`, `Tri31`, `NineNodeQuad`,
`EightNodeQuad`, `SixNodeTri` (and `FourNodeQuad3d`) carry an element-level
density member `rho` and build their mass matrix as
`if (rho == 0) rhoi = theMaterial[i]->getRho(); else rhoi = rho;`. Upstream
`sendSelf`/`recvSelf` pack `[tag, thickness, b0, b1, pressure, alphaM, betaK,
betaK0, betaKc, ...]` but **never `rho`**, and the no-arg constructor the
`FEM_ObjectBroker` uses on restore initialized `thickness`/`pressure` but **left
`rho` uninitialized** (the one exception, `FourNodeQuad3d`, happens to zero it).

Consequence: after `ops.database('File',p); save; wipe; restore` (or any
OpenSeesMP inter-process send), the reconstructed element's `rho` is garbage
heap memory. Garbage ≠ 0, so the `if (rho==0)` gate FALLS THROUGH and the mass
matrix is built from the garbage value instead of the material density —
**non-deterministic across fresh processes** (uninitialized read; measured
denormalized/`-2.5e-297` and `-2.5e6` generalized eigenvalues), an **indefinite
M**, and a transient (Newmark) restart that diverges from the uninterrupted run
(measured rel. diffs 0.5 … 5.6 across identical fresh runs).

Why it hid so well:
- The **committed nodal disp/vel/accel round-trip bit-exact** (FE_Datastore stores
  them fine), so a plain nodeDisp `database_roundtrip` PASSES — the corruption
  only surfaces once the CONTINUED transient uses the mass matrix.
- The **material `rho` is restored correctly**; it's the element's *shadow* `rho`
  that's garbage and *overrides* the material via the `if (rho==0)` gate.
- With **no LoadPattern / no excitation** the restart is bit-exact (u≡0, so a
  wrong M is never exercised) — which masks the bug and made it look
  LoadPattern-triggered. It is not; it's mass-triggered.

Diagnostic recipe: build a small dynamic plane-strain model WITH a load, run a
few Newmark steps, `save→wipe→restore`, continue, and compare to the
uninterrupted run (or `ops.eigen` before-vs-after — negative eigenvalues = smoking
gun). See `tests/test_quad_tri_rho_db_restart.py`.

Fix (strictly additive, this PR): append `rho` to each element's send/recv data
Vector (+1 slot) AND add `rho(0.0)` to each blank constructor. Vanilla bug —
every mainline OpenSees user doing a database or parallel dynamic restart with
these elements is affected; upstream-PR candidate.

**General lesson (this is the 2nd fork bug of this exact shape):** any scalar
member an element/material reads at analysis time MUST be either serialized in
`sendSelf`/`recvSelf` OR initialized in the no-arg broker constructor — ideally
both. A `save→wipe→restore→continue` (not just `→compare committed disp`)
round-trip test is the cheap guard; the plain nodeDisp round-trip is blind to
any state that only feeds future steps (mass, damping, committed internal vars).
### LadrunoPorousOverlay pattern: the TimeSeries/load factor is IGNORED by design — the overlay owns its force amplitudes
- **Bites:** attaching a `timeSeries` to `pattern LadrunoPorousOverlay ...` (or expecting `loadConst`-style factor scaling) does nothing: the injected nodal forces are always `+Q·p_committed` at full amplitude — the pore-pressure field, not a factored load, is the amplitude. Silent expectation mismatch if you try to "ramp" the overlay.
- **Why:** the overlay is a domain ENGINE riding the LoadPattern plumbing (H5DRM precedent); its `applyLoad(time)` ignores `time` and any series. A one-time informational notice prints if a series was assigned. Note the python/Tcl surface cannot even attach a series to it structurally (verified 1.E-ii, 2026-07-14).
- **Workaround/status:** ramp the SOLID loads (they live in ordinary patterns); stage the fluid via `-pInit` / `-staticMode`. ADR-73 §4.1; P1 battery gates the bit-exactness of "factor changes nothing".

### One overlay per water body: two `LadrunoPorousOverlay`s must not share elements — and a shared water table needs ONE overlay, not one per soil layer
- **Bites:** modeling a layered deposit as one overlay per layer disconnects the fluid: each overlay owns an independent p-field with its own drained set — no cross-layer flow, wrong consolidation. Conversely two overlays CLAIMING the same element double-count the fluid.
- **Why:** the p-field lives per-overlay (own CSR system); continuity exists only inside one region. Layered properties belong to `-layer` blocks INSIDE one overlay.
- **Workaround/status:** element overlap across overlays is a snapshot FATAL ⟨A-13⟩ (P1 battery-gated); the one-overlay-per-water-body rule is a modeling discipline the P4 guide owns. 2026-07-14.

### FileDatastore silently CLOBBERS same-type objects stored with the same (dbTag, commitTag, SIZE) — pack multi-part sendSelf payloads under distinct dbTags
- **Bites:** a class whose `sendSelf` sends TWO Vectors (or two IDs) on its one dbTag works fine — until a model size where the second object's length coincides with the first's; then the later write overwrites the earlier entry and `recvSelf` restores garbage (measured: LadrunoPorousOverlay configs where `nPI + 6·nLay + 2·nRN == 38` corrupted the scalar block → huge bogus counts → abort/exit-127 on restore; other coincidences restored silently-wrong state). Non-deterministic-LOOKING because it is config-size-dependent.
- **Why:** `FileDatastore` files entries per (type, size) and keys by (dbTag, commitTag) within that file — same type + same size + same tags = same slot.
- **Workaround/status (2026-07-14, ADR-73 P1):** the upstream matDbTag idiom — grab a second `theChannel.getDbTag()` for the payload object and transmit it inside the first (fixed-size) block. Overlay fixed this way; DB round-trip battery-gated. Audit any future multi-send class for size coincidences.

### RESOLVED (root cause = the quad/tri rho row above, fixed #577): transient DB restart looked "LoadPattern-triggered, non-deterministically corrupt" — the real trigger was mass-driven response over garbage element rho
- **Bites (historical, 2026-07-14 / ADR-73 P1):** `database File` + `save`/`wipe`/`restore` of a TRANSIENT model did not reproduce the continuous run, non-deterministically across fresh processes (rel 0…5.6; restored-pencil eigenvalues went NEGATIVE), while every python-visible restored quantity checked out (per-node committed d/v/a, SP sets, pattern state) — which makes this class of bug look like *your* serialization bug. It is not.
- **Why (corrected by the #577 root-cause session):** quad/tri elements never serialized their ELEMENT-LEVEL `rho` and the broker ctor left it uninitialized — heap garbage overrode the correctly-restored material density via the `if (rho==0)` gate, corrupting the restored mass matrix. The P1-era "any LoadPattern in the stream triggers it" attribution was an ARTIFACT: the no-pattern control was nearly motionless, so the corrupt M went unexercised. Lesson: a "pattern-correlated" restart failure can be a mass-path bug that only excitation reveals.
- **Still open, separate:** FileDatastore does NOT restore a sibling `Plain` pattern's NODAL loads (base reaction collapses post-restore — keep restart gates self-contained); a `-pInit` list overlay crashed on restore in one 1.E-ii observation (unre-verified since the rho fix). `ops.logFile()` remains the only way to capture opserr from python.
- **Workaround/status (2026-07-15):** FIXED by #577 (rho serialized + zero-init, all six quad/tri elements). Deformable transient DB restart is SUPPORTED again and hard-gated: `tests/test_ladruno_overlay_framework.py` 6c (promoted from xfail per its XPASS contract — 6/6 fresh-process restarts must reproduce; the multi-fresh-process structure is the regression net for uninitialized-serialization bugs, since a single run can land on benign heap garbage) + `tests/test_quad_tri_rho_db_restart.py` (zone_a). The two fork-side hardenings found while chasing this stay in force: the FileDatastore same-size clobber (dbTag2 idiom — dedicated row) and lazy `uSnapshot_` re-derivation (restore-ordering hygiene).

### `LadrunoStaggeredAnalyze` overrides `-subcycle` while driving — the driver syncs the fluid every step by construction
- **Bites:** an overlay built with `-subcycle N>1` (or `-subcycle auto`) accumulates the Δu window across N commits under plain `analyze`, but under the ADR-73 P2 driver `LadrunoStaggeredAnalyze` the fluid is advanced (and committed) at EVERY driver step — the configured N is ignored. If you expect the same subcycled fluid cadence you configured, you don't get it while driving.
- **Why:** the driver's whole point is the iterated per-step fixed-point solve; a multi-commit accumulation window is incompatible with re-solving the solid against the current-step fluid iterate. The latched `onDomainCommit` (SM_MARCH branch) does `commitFluid` + counter reset only — no window accumulation, no extra `advanceTrial`.
- **Workaround/status (2026-07-17, ADR-73 P2):** intended, not a bug. A one-time advisory prints ("LadrunoStaggeredAnalyze overrides -subcycle while driving") when the driven overlay has `subcycleN>1` or `-subcycle auto`. A pending window at driver entry is caught up with an early fs1 sync (`catchUpPendingWindow`) so the first driver advance doesn't pair a multi-commit Δu with the single driver dt. Counters are zeroed at each latched commit, so a post-driver plain `analyze` restarts its window cleanly.

### `updateMaterialStage` CANNOT reach a `LoadPattern` (overlay moduli must be re-set via the `parameter ... loadPattern` route after a stage flip)
- **Bites:** you flip a soil constitutive stage mid-analysis with `updateMaterialStage $matTag $stage` (PDMY/PM4Sand elastic→plastic), expecting the `LadrunoPorousOverlay`'s fixed-stress `L` factor (which uses the drained skeleton moduli) to follow the new stage. It does not — the overlay keeps its stage-0 moduli, so `L` is stale (stable, but convergence-degraded — the fixed-stress split's L is a preconditioner-like term, not a physics error).
- **Why:** `MaterialStageParameter` registers only the FIRST accepting ELEMENT in its domain scan and never scans load patterns (the ADR-71 sibling-broadcast trap, family-documented). A `LoadPattern` subclass like the overlay is simply not on the path `updateMaterialStage` walks.
- **Workaround/status (2026-07-17, ADR-73 P2):** the transport contract is explicit — after a stage flip the USER re-sets overlay moduli through the EXISTING parameter route: `parameter $p loadPattern $overlayTag E $newE` (or `nu`, `layerE $i`, `layerNu $i`), which marks the overlay `moduliDirty_` and lazily rebuilds `aS_`/`aL_` at the next fluid use. A flip without a re-set keeps stage-0 `L`. The PDMY staged-liquefaction battery + the P4 guide pin the recipe.

### Overlay `-layer` moduli overrides use a `> 0` sentinel — `layerNu 0` is unreachable via the parameter route (rejected loudly)
- **Bites:** trying to set a per-layer Poisson ratio of exactly 0.0 through `parameter $p loadPattern $tag layerNu $i` + `updateParameter` fails with a warning, even though nu = 0 is a physically legal value (and the global `nu` accepts it).
- **Why:** the `Layer` struct encodes "inherit the overlay-global value" as `nu <= 0` (P1 sentinel, serialized that way); a stored layer nu of 0.0 would be silently re-interpreted as "unset" by `resolveCellModuli`, turning the update into a no-op. The P2 panel (robustness-7) flagged the silent path; the fix rejects it loudly instead.
- **Workaround/status (2026-07-17, ADR-73 P2):** set the GLOBAL `nu` to 0 (parameter id `nu`) and leave the layer inheriting, or use a tiny positive value. Changing the sentinel to an explicit per-field override flag would touch the serialized layer payload — deferred until a real user needs layered nu = 0.

### `LadrunoStaggeredAnalyze` failure surface in openseespy is split: parse/analysis-setup fatals raise, run-time aborts return a negative int
- **Bites:** `ops.LadrunoStaggeredAnalyze(...)` raises `OpenSeesError` for a static-analysis-active or no-transient-analysis misuse (detected in the command body), but returns a plain negative int (no exception) for run-time aborts — empty driven set, bad args caught in the core, solve/fluid failures, maxIter (−1/−2/−3/−6/−7). `assert rc == 0` catches the second class only.
- **Why:** the run path deliberately mirrors classic `analyze` (negative return, "failed, returned: N error flag" print) so scripted retry logic works the same for both commands; the command-body fatals happen before a result exists, so the wrapper turns them into exceptions (`Py_ops_...` NULL-return convention).
- **Workaround/status (2026-07-17, ADR-73 P2):** by design; check the integer return like you would for `analyze`, and wrap in try/except only for setup misuse. Panel 2.D robustness-9 recorded the split; the P2 battery handles both surfaces.

### Staggered-overlay twin models: plain quad `b1 b2` = body FORCE/VOLUME, `LadrunoUP -body` = ACCELERATION — copying the same number silently unloads (or double-loads) the solid
- **Bites:** building the ADR-73 staggered twin of a monolithic `LadrunoUP` model (plain `quad` + `LadrunoPorousOverlay`), you copy the monolithic `-body 0 $bY` value into the quad's trailing `b1 b2` args. The staggered solid then carries `bY` per unit VOLUME instead of `rho_mix*bY`, and with a hydrostatic overlay `+Q·p` force the net solid load can cancel to ~zero — the P2 battery measured settlement 1e-9 vs 3.5e-4 (u-trace rel diff exactly 1.0) before the fix. Nothing errors; the fluid side looks perfect.
- **Why:** `FourNodeQuad` applies `b` directly (`P(ia) -= dvol*shp*b`, FourNodeQuad.cpp:900 — force density, rho only builds mass), while `LadrunoUP`/upstream `quadUP` scale `-body` by the mixture density (acceleration semantics). Same flag name, different units.
- **Workaround/status (2026-07-18, ADR-73 P2):** staggered twin recipe — quad gets `b2 = rho_mix * bY_accel` (full mixture weight as force/volume; the overlay's `+Q·p` supplies the pore-pressure part of effective stress), overlay `-fluidBody` keeps the acceleration form (`f_seep` scales by `rhoF` internally, matching `LadrunoUP -fluidBody`). Pinned by battery gate (d)(ii) (`tests/test_ladruno_overlay_driver.py`); the P4 guide inherits the recipe.

### Honest-p `LadrunoUP` + upstream `CentralDifference` = Richardson-unstable pore pressure (CD leapfrogs a pure diffusion operator)
- **Bites:** running a monolithic `LadrunoUP` (honest-p) model under the upstream explicit `CentralDifference` integrator produces unstable / garbage pore pressure at ANY dt -- p oscillates and grows even when the displacement field looks plausible for a while. No warning fires; the run may limp along before blowing up. INCUBATION TRAP (P3 measured): at small dt the blow-up incubates for ~1300 steps (dt=2e-4; ~200 at dt=1e-3) while tracking the reference to ~1e-4 -- a short validation march "passes" and the production run dies later. Never certify this combination from a few-hundred-step check.
- **Why:** honest-p carries NO p-mass (the fluid row is first-order: `S p' + H p = ...`), so central difference's second-order leapfrog applied to that row is the Richardson / DuFort-Frankel-class explicit scheme for a pure diffusion operator -- unconditionally unstable in the Richardson form. This is structural (the integrator discretizes a row with no inertia), not a tuning problem: no dt, damping, or mass-scaling choice fixes it.
- **Workaround/status (2026-07-18, ADR-73 P3):** use the overlay explicit lane instead (plain solid elements + `LadrunoPorousOverlay -fsL zero` under `CentralDifferenceLadruno` at dt <= 0.5x the undrained pencil -- the fluid solve stays implicit SPD at commit), or keep `LadrunoUP` under an implicit integrator (Newmark/HHT). Pinned by the P3 battery's expected-bad gate (e): the CD+LadrunoUP run must NOT track the reference -- if that gate ever starts passing, this row is stale.

### Overlay energy accounting: the `+Q.p` coupling work rides ADR-69's ULW (external-load-work) channel -- closure holds, attribution is merged
- **Bites:** reading an ADR-69 `EnergyBalanceRecorder` breakdown on an overlay run, there is no separate "pore-coupling work" channel and the external-load work looks inflated -- you cannot tell coupling work from genuine external load work in the per-channel numbers.
- **Why:** the overlay injects its forces through `Node::addUnbalancedLoad`, and the ADR-69 kernel's ULW = integral of v^T P_ext dt reads `Node::getUnbalancedLoad` -- so the `+Q.p` forces are INSIDE the external-work channel by construction (verified at P3 pin 3.A). The closure residual ERR therefore stays within the ADR-69 bound (measured by battery gate (g)); only the ATTRIBUTION is merged, not the balance.
- **Workaround/status (2026-07-18, ADR-73 P3):** documented, not silently absent -- energy closure on overlay runs is trustworthy; per-channel attribution of coupling work is not separable. P4 may split it into its own channel; no recorder code change at P3.

### `-fsL zero` under a quasi-static implicit fs1 march = the naive drained split -- diverges in ~4 steps at soil coupling (explicit-lane-only setting)
- **Bites:** an overlay built with `-fsL zero` (the ADR-73 P3 explicit-lane setting) but driven by a quasi-static / implicit fs1 march (plain `analyze` under Newmark at consolidation time scales) diverges within ~4 steps, ~10 orders of magnitude, at realistic soil coupling strength tau = (alpha^2/K_dr)/storage ~ 1e3 (measured, ADR-73 SS3.2) -- loudly, not wrongly.
- **Why:** L = 0 removes the fixed-stress relaxation entirely, so the split IS the naive drained split the ADR bans for the implicit lane. The setting exists ONLY for the explicit lane, where stability is governed by dt <= 0.5x the discrete undrained pencil (E7.2: the L=0 implicit-fluid boundary = exactly 1.000x the pencil) and no iteration/L is needed.
- **Workaround/status (2026-07-18, ADR-73 P3):** by design. The parser prints a one-time loud advisory at `-fsL zero`; `LadrunoStaggeredAnalyze` refuses FSL_ZERO overlays with a loud fatal (iterating with L = 0 is the same drained split). Use `-fsL classic|oedometric` for implicit/driver lanes; reserve `zero` for `CentralDifferenceLadruno`/explicit runs at dt <= 0.5x the (overlay-aware) `criticalTimeStep` pencil.

### SMS mass scaling + LadrunoPorousOverlay = certified-stable-but-actually-UNSTABLE (UNSUPPORTED until ADR-73 P3b)
- **Bites:** you run `CentralDifferenceSMS`/selective mass scaling (lumped OR consistent Olovsson) on a model carrying a `LadrunoPorousOverlay` and pick `dtTarget` trusting the SMS report. SMS sizes every element against the DRAINED per-element pencil (`elementCriticalDt` receives no overlay augmentation), but overlay cells' true explicit limit is the UNDRAINED pencil — AT LEAST `sqrt(1+Kf/(n*M_oed))` smaller (the material formula is a LOWER bound on the discrete per-element factor: measured ~26x on the e72 soft soil vs the 21.4x material value — mode-shape excess, ADR-73 §12 P3 item 1). Mass scaling raises both pencils by the same `sqrt(s)`, so after sizing the overlay cells are still factor-x short of `dtTarget`: the run is certified stable and then blows up.
- **Why:** the ADR-73 P3 advisory augmentation lives in the shared `computeCriticalTimeStep` scan (the `criticalTimeStep()` report IS overlay-aware, even under SMS), but SMS sizing calls the per-element `elementCriticalDt` kernel directly and no caller passes the `Kadd` seam yet (panel 3.D advisory-critic finding H).
- **Workaround/status (AMENDED 2026-07-19, ADR-73 P3b — FIXED for the LUMPED builder, warning retired):** BOTH builders now pass the per-element undrained augmentation through the `elementCriticalDt` Kadd seam, so SMS sizing prices the UNDRAINED pencil (the blanket UNSUPPORTED warning is RETIRED; a one-time INFO line `SMS sizing priced the UNDRAINED pencil for N overlay-owned elements (ADR-73 P3b)` prints instead; residual drained-priced configs are exactly `getUndrainedAugmentation`'s loud refusals — singular S_e, exotic-ndf size mismatch, failed snapshot). **Lumped `CentralDifferenceSMS` is the composability deliverable**: the certified dtTarget march is battery-HARD-gated stable (4000 steps, gate (d)). **Consistent/Olovsson `CentralDifferenceSMSConsistent` remains limited — see its own row below** (under-delivers on the coupling mode; loud warning ships). ADR §12 P3b item 5 governs.

### Consistent/Olovsson SMS + LadrunoPorousOverlay: correct undrained PRICING but measured UNDER-DELIVERY on the coupling mode (use lumped SMS)
- **Bites:** `CentralDifferenceSMSConsistent` on an overlay model, post-P3b: the sizing report correctly prices the undrained pencil (INFO line prints), yet the certified `dtTarget` march can still diverge — measured uniform ~×1.83/step growth from step 1 at dtTarget = 3× the unscaled pencil on the e72 column (battery gate (d) EXPECTED-LIMITED record), while the LUMPED builder's certified march on the same model is stable 4000+ steps.
- **Why (ADR §12 P3b item 5, panel-checked):** the Olovsson centroid-preserving M̄ blocks add inertia only to the NON-RIGID element modes (that is their design — element mass distribution preserved). The overlay's undrained volumetric coupling mode carries a large rigid-translation component that stays UNSCALED, so the coupled frequency scales by less than √s and the certified step over-promises for overlay-owned cells. This is a scheme interaction, not a wiring bug — the lumped builder injects real nodal mass and delivers.
- **Workaround/status (2026-07-19, ADR-73 P3b):** use lumped `CentralDifferenceSMS` with overlays, or size dt from the overlay-aware `criticalTimeStep()` report. The consistent builder keeps the undrained pricing (honest report) and prints a loud one-time warning when it scales overlay-owned elements ("the Olovsson centroid-preserving blocks under-scale the undrained COUPLING mode").

### Overlay explicit lane on a ZERO-drainage undamped column: no asymptotically stable dt — the "stable at <= 0.5x pencil" pin is HORIZON-relative there (secular energy pumping)
- **Bites:** an `-fsL zero` overlay under CentralDifferenceLadruno on a pathological column with effectively zero drainage (k-bar ~ 1e-11) and NO damping shows slow secular energy growth at ANY dt: measured steps-to-blowup 48k @ 0.4x pencil, 30k @ 0.5x, 7.8k @ 1.0x, 855 @ 3.0x (growth per unit time ~ O(dt)). A short march looks stable; a long enough march always diverges. The drained control (same solid, no overlay) stays bounded — the pumping is the frozen-force staggered coupling, not the solid integrator.
- **Why:** the ADR-73 explicit lane holds the +Q.p force FROZEN over each step (fs1 operator split). With zero drainage and zero damping there is no dissipation channel to absorb the O(dt) splitting-energy input, so it accumulates secularly. This is exactly the ADR §3.4 A-7 hedge landing: the ZPC-1988 stability proof covers THEIR scheme, not the frozen-force variant — the E7.2 toy "stable" boundary was itself horizon-limited on this pathology. Any physical drainage or damping absorbs the input (the ZS84 two-leg gate marches 6.5k steps clean and dt-converges; a fast-draining column is clean).
- **Workaround/status (2026-07-18, ADR-73 P3):** real soils drain and real models have damping — the lane is production-fit there (battery gates (a)/(f) green). For near-zero-k undamped models: add any small alphaM Rayleigh or use the implicit lane. The frozen "L=0 stable at <= 0.5x undrained pencil" pin is re-stated as horizon-relative on this pathology (0.5x buys ~4x the blowup horizon of 1.0x, ~35x of 3.0x); ADR §12 P3 entry records the measurement; battery gate (b) documents itself as a fixed-horizon envelope and prints the secular-pumping record block.
- **AMENDED 2026-07-19 (ADR-73 P3b, measured — gate (g)):** this row is now **IMPLICIT-`-fsL zero`-lane-specific**. Under the P3b fully-explicit fluid (`-fluidUpdate explicit`, load-application-time advance), the SAME pathology legs are **BOUNDED to their caps**: 60k steps @ 0.4× pencil, 45k @ 0.5×, 15k @ 1.0× (vs implicit 48344/30018/7761). The lumped forward p-step does not feed the secular channel the implicit-at-commit variant does (direction consistent with the toy: e74 explicit boundary 1.32× > e72 implicit 1.000×). Near-zero-k undamped models are therefore BETTER served by `-fluidUpdate explicit`; the horizon-relative caveat stays for the implicit lane. ADR §12 P3b item 2.

### `remove loadPattern` does NOT fire domainChange (unless the pattern owned SP_Constraints) — cached pattern pointers dangle silently
- **Bites:** any object caching a `LoadPattern*` across steps (a recorder result source, an engine seam, a driver) keeps a freed pointer after `remove loadPattern $tag`: the interpreter command DELETES the pattern object (`OpenSeesMiscCommands.cpp` remove path), but `Domain::removeLoadPattern` calls `domainChange()` only when the pattern carried SP_Constraints — a LadrunoPorousOverlay owns none, so NO domain-change stamp bumps and NO recorder/model rebuild fires. The ADR-73 P4 panel measured the consequence: an early `-overlay` recorder build cached the pattern pointer and wrote subnormal garbage (9.9e-312) into every post-removal row — silent use-after-free, allocator-dependent whether it corrupts or segfaults.
- **Why:** element/node removal goes through paths that mark the domain changed; load-pattern removal is only conditionally marked. Anything keyed off `hasDomainChanged()` (LadrunoRecorder writeModel rebuild, MODEL_STAGE rollover) will NOT observe a pattern removal.
- **Workaround/status (2026-07-18, ADR-73 P4):** never cache a `LoadPattern*` across commits — re-resolve by tag each use (`domain->getLoadPattern(tag)` + classTag check) and act loudly/zero-fill when absent. `OverlayPressureSource::evaluate` and the Monitor overlay path both do this now (battery gate (i) pins it: post-removal rows identically 0.0, file readable). Audit any future pattern-consuming seam for the same hole.

### Commit-hook "next-newStep" force refresh is TWO commits stale under the CDL family — fatal for explicit-in-time load values, fine for implicit
- **Bites:** a LoadPattern that refreshes its injected nodal forces at `Domain::commit` (the ADR-73 overlay hook idiom) and expects the solid's NEXT step to move under the refreshed value: under `CentralDifferenceLadruno` (and its SMS subclasses) the displacement advance of step k+1 uses `Aprev` — the acceleration formed during step k's SOLVE, i.e., with the forces applied at step k's newStep, i.e., the values committed at step k−1. The refreshed value first moves u at step k+2. For the ADR-73 explicit-fluid lane this one extra lag step was UNCONDITIONALLY destabilizing (vertical p-checkerboard, growth ∝ dt, no stable Δt — toy-lag replica reproduced the divergence to 7 significant digits); the implicit-fluid lanes tolerate it (dissipation absorbs an O(Δt) lag; P1–P3 shipped and gated with it).
- **Why:** CDL's leapfrog does advance-then-updateDomain: `newStep` moves u with the stored `Aprev`, THEN applies loads and solves for the NEW acceleration. Equilibrium pairing is correct for TIME-series loads (evaluated at the right time); it is one step stale for loads whose VALUE is refreshed at commit.
- **Workaround/status (2026-07-19, ADR-73 P3b):** value-refreshed injected forces that must pair tightly with the current step belong at LOAD-APPLICATION time (inside `applyLoad`, after newStep's trial set, reading TRIAL state), not at the commit hook — the P3b explicit lane does exactly this (advance-from-committed on the trial window Δu, idempotent under re-applies). Audit any future commit-hook force refresher against this pairing before assuming ZPC-class stability results transfer. ADR §12 P3b item 1.

### Explicit-fluid lane + `-subcycle N>1`: the UNDRAINED CFL binds the SYNC interval N*dt — the implicit lane's large-N freedom does NOT transfer
- **Bites:** carrying the E7.3a intuition ("all N <= 50 stable, error ~ N^1.2") from the implicit `-fsL zero` lane to `-fluidUpdate explicit`: N=4 at 0.4x the undrained pencil (sync interval 1.6x) diverges in ~400 steps on the e72 column — where the SAME dt at N=1 is bounded for 60k steps. Measured and toy-twinned (C++ step 395 / toy step 408).
- **Why:** with the fluid explicit, the coupled staggered stability is set by the frozen-force interval = the sync interval N*dt (the undrained stiffening acts once per sync); the implicit-at-commit fluid absorbed exactly that stiffening, which is where E7.3a's freedom came from. Subcycling also has no purpose on the explicit lane — the fluid step is an axpy, there is no solve to amortize.
- **Workaround/status (2026-07-19, ADR-73 P3b):** `-subcycle auto` under `-fluidUpdate explicit` resolves N=1 with a notice; manual N>1 prints a loud one-time sync-CFL warning (keep N*dt within the N=1 margin). Battery s5 pins the expected-diverge demo, the bounded sync-0.4x leg, and the auto->1 notice. ADR §12 P3b item 9.
### `analyze()` returns rc=0 on a NaN-poisoned system — "rc==0" does NOT mean "numbers"
- **Bites:** trusting the analyze return code as a health signal. `FullGeneral` + `algorithm Linear` solved a stiffness matrix full of NaN (degenerate-eas 1/det ≈ 1e17 blowup) and reported SUCCESS; `nodeDisp` was NaN with no error printed at any level.
- **Why:** LAPACK `dgesv` propagates NaN without setting its info flag, and no layer above it (SOE / algorithm / analysis) checks the solution for finiteness.
- **Consequence:** element-level degeneracy/finiteness guards are the ONLY defense against silent-NaN results in linear analyses. Tests asserting "the analysis must fail" must either construct a genuinely SINGULAR system (zeroed row → dgesv info>0 → rc!=0) or assert output finiteness explicitly — never rely on NaN tripping the solver.
- **Status (2026-07-20):** found while probing the eas degeneracy-guard axis-collapse hole; the guard fix makes eas elements refuse loudly (their zeroed block → singular SOE → honest rc!=0), but the generic solver blind spot remains (upstream-class, unfixed).

### `ID::insert` is a SORTED-SET insert, not an append — so every sparsity structure built through it is ascending by construction (and that is what makes a binary search in `addA` legal)
- **Bites (in the good direction, and in the bad one).** Good: you want to replace an `addA` linear scan with a binary search and need to know whether the column/row indices are ordered. They are — `ID::insert` (`SRC/matrix/ID.cpp`) is a **binary-search insertion into a sorted array with dedup** (`if (x == dataMiddle) return 1; // already there`), so any structure filled by it comes out **strictly ascending and unique**. `UmfpackGenLinSOE::setSize` builds each CSC column that way (`ID col; col.insert(diag); col.insert(row)...` → `Ai.push_back(col(i))`), and the PARDISO SOE does the CSR equivalent. Bad: it is easy to read `insert` as "append" and conclude the opposite, or to assume ordering that a *different* fill path does not provide — `PARDISOSymLinSOE` (the dead pair, see its own row) appends adjacency **in ID order without sorting** and would hand PARDISO a non-ascending CSR.
- **Why it matters:** binary search on a non-ascending column does not fail loudly — it returns the wrong slot or reports "absent", so an element's stiffness lands in the wrong entry or is silently discarded. A converged, plausible, wrong answer.
- **Workaround/status:** "by construction" is not "enforced" — the construction can change under you. Both live SOEs now **CHECK** the invariant at the end of `setSize` and refuse (`return -1`) with a loud message: `PARDISOGenLinSOE` (ADR-75 P1f, [#636](https://github.com/nmorabowen/OpenSees/pull/636)) and `UmfpackGenLinSOE` (ADR-75 P1g). The check is O(nnz) once per `setSize` against an assembly loop that runs for the whole analysis — free. General rule for this fork: when a fast path depends on an invariant someone else establishes, **check it where you depend on it**, don't cite the place that happens to establish it today. *2026-07-25 (ADR-75 P1f/P1g).*

### Any leftmost-pivot quicksort fed a graph-derived list is deterministically O(n²) — MapOfTaggedObjects iterates ASCENDING, so the input is already sorted
- **Bites:** `MPIDiagonalSOE::setSize` sorted its DOF list with a hand-rolled `q_sort` (leftmost pivot). The DOF graph is `std::map`-backed, its iterator hands tags over ascending ⇒ the sort's textbook worst case, per rank, every domainChanged: 501 s at the 2.0 M np8 rung, ~N^1.9 — masqueraded for years as "setup cost". The same trap arms ANY hand-rolled pivot sort downstream of a `MapOfTaggedObjects`/`ArrayOfTaggedObjects` iteration (both are ascending for dense tags).
- **Why:** classic quicksort worst case = sorted input + first-element pivot; map iteration guarantees sorted input.
- **Workaround/status:** FIXED in-place with `std::sort` ([#593](https://github.com/nmorabowen/OpenSees/pull/593), output provably identical). Rule: never hand-roll a sort on tag/equation lists; grep found no other `q_sort` copies. *2026-07-22.*

### `ID::getLocation` (or a full constraint/domain iteration) inside a per-node / per-element / per-DOF-group loop is a per-rank quadratic — the "scan-in-loop" family
- **Bites:** five independent instances cost real wall: `TransformationConstraintHandler::handle()` (N^1.94, 14.4 s at 2.0 M np8 — element classification scanning the SP list per element-node), the `TransformationDOF_Group` SP-only ctor (swept EVERY domain SP per constrained node — invisible until the handler fix landed), `PlainHandler` per-node `getMPs()` AND `getEQs()` sweeps, and the `-4` fixup full-MP sweep in `DOF_Numberer` + `PlainNumberer` (both variants each). All MP/SP-count-driven: zero cost on unconstrained decks, quadratic on slab meshes (constrained area ∝ N) and tie-heavy decks.
- **Why:** `ID::getLocation` is a linear scan (`ID.cpp`); constraint iterators restart from scratch each call. O(outer) × O(scan) with both ∝ N/P.
- **Workaround/status:** ALL FIXED with one-pass hash/multimap indexes, order-preserving, byte-identity-gated ([#595](https://github.com/nmorabowen/OpenSees/pull/595), [#598](https://github.com/nmorabowen/OpenSees/pull/598); the parallel numberer's own copies in #592). Audit method that found them: `dc.*` profiler brackets + a fixed-np rung sweep (exponent), then a fixed-V np-sweep (per-rank vs global discrimination: per-rank quadratics FALL ~1/P², global-serial ones are np-invariant and mimic an Amdahl fraction in scaling studies). **Fix one, re-measure — the ctor sweep was invisible behind the handler scans.** *2026-07-22.*

### `ParallelNumberer` silently FUSES node-less DOF groups (Lagrange multipliers) across ranks — `getRef()` returns −1 for all of them
- **Bites:** under MP with the Lagrange handler, every rank's multiplier DOF groups share ref=−1; the stock gather-merge dedups by ref ⇒ all of them collapse into ONE merged vertex ⇒ silently wrong numbering. Never observed in production only because the fork's MP lanes use Transformation/Plain handlers.
- **Workaround/status:** `LadrunoParallelNumberer` hard-errors on ref<0 with a message naming the gap; stock still fuses. Use Transformation/Penalty handlers under MP. *2026-07-22 (found in the ADR-74 N2 review).*

### `ParallelPlain` numbers vertex 0 LAST (the tag-0 quirk) and `ParallelNumberer::numberDOF(ID&)` never worked at all
- **Bites:** (a) stock `ParallelPlain` checks "already ordered" against a zero-filled ID, so merged tag 0 always reads as present and gets pushed to the end — a valid but surprising permutation (bandwidth outlier on the first node). (b) The `numberDOF(ID& lastDOFs)` variant (SP/DomainDecomposition lane) has its numbering call commented out upstream and a mismatched recv layout — any caller gets garbage start-DOFs.
- **Workaround/status:** `LadrunoParallelPlain` fixes (a) (G1b-gated: valid bijection, differs from stock by design); the Ladruno override hard-errors on (b). *2026-07-22.*

### Under `system MPIDiagonal` the global equation numbering exists only TRANSIENTLY — setSize rewrites every DOF id to rank-local 0..n−1
- **Bites:** any oracle/debug dump of equation ids taken AFTER analysis setup shows rank-LOCAL ids (shared boundary nodes legitimately disagree across ranks); naive cross-rank identity checks fail on correct runs. MUMPS SOEs do NOT localize — dumps keep globals.
- **Why:** `MPIDiagonalSOE::setSize` builds its shared-DOF exchange from the globals, then deliberately compacts per rank ("renumber DOFs 0 through size") and has FE elements re-cache.
- **Workaround/status:** by design, no physics defect. Numbering oracles must be two-deck: `system Mumps` for strict global identity, MPIDiagonal for end-state identity (`tests/test_adr74_numberer_1.py`). *2026-07-22 (ADR-74 N0's first catch).*

### MUMPS: symbolic analysis is NOT in `setSize` (runs at first solve), and error −13 in "substitution" is plain out-of-memory
- **Bites:** (a) profiling `dc.setSize` on the implicit lane and expecting the MUMPS analysis there — it's in the first `solveCurrentStep` instead (`MumpsParallelSolver::setSize` only sets `needsSetSize`). (b) A 1.0 M-node (3.1 M-eq) LU does not fit a 64 GB box; the failure surfaces as `Error -13 returned in substitution dmumps()` per rank, easily misread as a numerical bug.
- **Workaround/status:** attribution: setSize = graph build + triplet fill (measured linear, [#604](https://github.com/nmorabowen/OpenSees/pull/604)); budget first-solve separately. Local implicit ceiling ≈ the 0.5 M rung. *2026-07-22.*

### `Graph::getVertexPtr` cost is set by the STORAGE the graph was built on — `ArrayOfTaggedObjects` is O(1) for dense tags, `MapOfTaggedObjects` is O(log n)
- **Bites:** an RCM BFS (or any per-vertex loop calling `getVertexPtr`) on a `MapOfTaggedObjects`-backed graph pays a `std::map::find` per read — ~1.6×10⁹ `find`s ≈ 25-45 min at 19 M nodes, the whole post-T0 numberer residual. It reads like "RCM is just slow"; it's actually the storage choice underneath `getVertexPtr`. `AnalysisModel::getDOFGroupGraph` builds on a `MapOfTaggedObjects`, so every parallel-numberer merge inherited it.
- **Why:** `Graph::getVertexPtr` → `TaggedObjectStorage::getComponentPtr`. `ArrayOfTaggedObjects` returns `theComponents[tag]` directly (O(1)) when the tag sits at its own index — true for the dense 0..N-1 tags graphs use; `MapOfTaggedObjects` is a red-black tree (O(log n)) always.
- **Workaround/status:** the T1 lever — `LadrunoParallelNumberer` builds its merged graph on **owned `ArrayOfTaggedObjects` storage** + a `tag → Vertex*` mirror, so RCM's BFS reads are O(1) and edge inserts skip `Graph::addEdge`'s two lookups entirely (`Vertex::addEdge` direct). numberDOF 15.5 → 6.4 s at 2.0 M ([#594](https://github.com/nmorabowen/OpenSees/pull/594)). Bit-identity holds because the adjacency `ID::insert` is a sorted set (insertion-order-canonical) and dense tags iterate ascending in BOTH storages. Rule: if you build a Graph you will read per-vertex, build it on array storage and confirm tags are dense. *2026-07-22.*

### `domainChanged()` re-runs the ENTIRE setup path (handle + numberDOF + setSize) on EVERY domain change — it is a K× cost, not a one-time one — and `eigen()` routes through it too
- **Bites:** treating "setup" as paid once at step 1. It is re-paid on every `domainChange`: apeGmsh emits one per **stage** (gravity → dynamic, staged construction, SSI decks — recorder MODEL_STAGE splitting, apeGmsh PR #633), ADR-51 element removal bumps the stamp per removal event, ADR-55 contact re-discovery likewise, and `DirectIntegrationAnalysis::eigen()` (`:323`) calls `domainChanged()` — so modal/FEAST runs pay it too. A 20-event removal history at 19 M re-crosses the hour line **even post-fix** if any setup term were still super-linear; a progressive-collapse (AEM) run at 10 M+ with hundreds of events multiplies whatever residual remains.
- **Why:** `domainChanged` unconditionally re-forms DOF groups, re-numbers, and re-sizes the SOE — there is no "incremental re-setup" path; a domain that changed by one element pays the full O(model) again.
- **Workaround/status:** the reason the ADR-74 fixes had to make EVERY setup term linear, not just fast-at-K=1: the numberer (T0/T1), setSize (#593), handle (#595), and the MP/EQ `-4` sweeps (#598) are all now ~linear, so K× is bounded. The per-run benchmark hides this — a single-`domainChange` deck (the G3 plane wave) is the K=1 best case and must be flagged as such; staged/removal/contact decks are where K bites. `LadrunoParallelPlain` (no RCM pass) further cuts the per-event cost on the explicit + MUMPS lanes. *2026-07-22.*

### `mpiexec` returns exit code 0 even when every rank died in Tcl-init or hit a script `error` — the process rc is NOT the MP-harness failure signal
- **Bites:** an MP test/sweep harness that gates on `$LASTEXITCODE` / subprocess rc passes green while the run produced nothing — a missing `TCL_LIBRARY`, a bad `source`, a deck `error`, or a per-rank abort all leave `mpiexec` reporting success. Sibling to the serial `analyze() rc=0 on NaN` quirk above, but worse because it hides TOTAL failure, not just a bad answer.
- **Why:** the fork's Intel-MPI `mpiexec` propagates the launcher's exit status, not the ranks' Tcl interpreter status; OpenSees does not `MPI_Abort` with a nonzero code on a Tcl error.
- **Workaround/status:** MP harnesses must assert on **artifact existence + content** (dump-file count, expected line count, a sentinel `puts` like `TIEGATE_DONE`/`ANALYZE_MS`), never on rc. Build-tree exes additionally need `TCL_LIBRARY` exported (the packaged `openseesmp.sh` sets it; a raw `mpiexec … OpenSeesMP.exe` does not) or every rank dies in init — silently, rc=0. *2026-07-22 (banked across the ADR-74 rung/tie/checkpoint harnesses).*

### `SRC/system_of_eqn/linearSOE/pardiso/` holds TWO symmetric PARDISO implementations — one of them is dead and unbuilt
- **Bites:** an agent asked to "add symmetric PARDISO" finds `PARDISOSymLinSOE.{h,cpp}` + `PARDISOSymLinSolver.{h,cpp}` already sitting in the directory and wires *those* up. They are a 2019 contributed prototype (M. Salehi, same author as the `Gen` pair) and are **not listed in `pardiso/CMakeLists.txt`** — only the `Gen` pair is compiled. They carry every defect ADR-75 P1a had to fix in the `Gen` pair and one more of their own: their `setSize` appends adjacency entries **in ID order without sorting**, relying on the adjacency already being ascending, so they would hand PARDISO a CSR whose columns are not guaranteed ascending.
- **Why:** the live symmetric path is `PARDISOGenLinSOE` with `matType != 0` (ADR-75 P1d, [#630](https://github.com/nmorabowen/OpenSees/pull/630)) — one class covering unsym + SPD + symmetric-indefinite, so the hardened factorization-reuse/`mtype`-derivation logic is shared rather than duplicated. The `Sym` pair was left untouched rather than deleted, since it is upstream-contributed and touching it widens the vanilla footprint for no gain.
- **Workaround/status:** use `system Pardiso -matrixType 1|2`. Do **not** wire up `PARDISOSymLin*`; if it ever gets built, its unsorted-column fill is the first thing to fix. Note the ascending-column requirement is now *checked* at the end of `PARDISOGenLinSOE::setSize`, so a future mistake here fails loudly instead of returning a plausible wrong answer. *2026-07-25 (ADR-75 P1d adversarial review).*

### Threaded MKL PARDISO is NOT byte-reproducible run-to-run — and one run per thread count cannot detect it
- **Bites:** anyone building a byte-identical CI gate, an oracle comparison, or an A/B that judges a code change by comparing outputs, while `MKL_NUM_THREADS > 1`. At 4 threads a 14³ Lane-B model returns **two distinct tip displacements across 10 runs of the SAME binary** (`2.13985383446834687` vs `…732`, ~1 ULP, a 5/5 split). At 1 thread it is 10/10 identical. It is **size-dependent**: an 8³ model is reproducible even at 4 threads.
- **Why it was missed:** ADR-75 P1 concluded "bit-identical at every thread count ⇒ threading introduces no FP drift ⇒ the determinism concern is Lane-3-only" from **one run per thread count**. That design cannot distinguish *deterministic* from *the same value came up twice* — with a 50/50 split, a single-sample check passes half the time. The claim was published in `RESULTS_p1_pardiso.md` and the ADR and stood for a day; both now carry corrections.
- **How it surfaced:** an A/B of the P1f `addA` change reported "ux DIFFERS" at 4 threads. The natural reading is "my change broke exactness". The correct next step was **not** to debug the change but to ask whether each binary reproduces *itself* — the OLD binary showed the identical 5/5 split, proving the variation was MKL's and the change was exact (confirmed at 1 thread, 10/10 identical between old and new).
- **Workaround/status:** pin `MKL_NUM_THREADS=1` for any byte-identical gate. PARDISO's own CNR control is `iparm[33]`, and it **is** available to this fork — Intel only forbids it when `iparm[1]=3` (parallel METIS) and we set `iparm[1]=2`. Not currently exposed as an option; wire it if a threaded reproducible mode is ever needed. Drift is last-bit, so this is a *reproducibility* problem, not an accuracy one — every configuration still matches UmfPack to `0.0` on a single run. **General lesson: "deterministic" is a claim about a DISTRIBUTION; never conclude it from n=1.** *2026-07-25 (ADR-75 P1f).*

### A PARDISO CGS *success* (`iparm[3]`) leaves the stored factors STALE — the next phase-33-only solve then answers the PREVIOUS matrix
- **Bites:** anyone adding `iparm[3]` (preconditioned CGS/CG) on top of a factorization-reuse gate. The natural reading of "CGS replaces the computation of LU" is that PARDISO updated something; it did not. On a CGS **win** the handle still holds the L/U of an *older* A — CGS merely used it as a preconditioner while iterating against the current A, so the returned x is correct but the factors are one tangent behind. A later phase-33 call (a second RHS, a re-solve without re-forming the tangent) is then a solve of the **previous** matrix that returns `error = 0`. Silent wrong answer, not a crash. The inverse case is equally easy to get backwards: on CGS **failure** under phase 23 PARDISO *does* refactor, so there the factors ARE current.
- **Why:** `theSOE->factored` answers "has A been re-assembled since the last solve", which is a different question from "do the stored factors correspond to A". With a purely direct solver the two coincide, which is why the distinction never had to exist before. `PARDISOGenLinSolver` now tracks them separately as `haveFactors` (is there a preconditioner at all) and `factorsCurrent` (does it match A); the phase-33 shortcut requires **both** `factored` and `factorsCurrent`.
- **Also:** the automatic direct fallback is documented for **phase 23 only**. Driving CGS from phase 33 turns a failed iteration into `error = -4` with no factorization — so "call phase 22 then phase 33" is *not* a valid decomposition once `iparm[3] != 0`. And `iparm[3]` must be left at 0 for the first numeric pass of a pattern: there is no previous factorization to precondition with. *2026-07-25 (ADR-75 P1e).*
### Returning `const Matrix &` from `getTangentStiff()` makes essentially EVERY OpenSees element non-re-entrant **by construction** — the `static Matrix` scratch buffer is the idiom, not a bug
- **Bites:** any plan that says "de-`static` the element kernels" as a bounded task before threading the element loop. This is not a handful of oversights: `SRC/element/` + `SRC/material/` hold **~5,600** function-/file-scope `static Matrix|Vector|ID` declarations across **587 distinct files** (1,686 Matrix, 3,535 Vector, 372 ID; re-counted 2026-07-25), plus **711** class-level `static Matrix|Vector|ID` members declared in `SRC/element/*.h`. `Element` itself owns instance-shared pools — `Element::theMatrices`, `theVectors1`, `theVectors2`, `numMatrices` (`Element.cpp:49-52`). Worked examples: `LadrunoBrick::getTangentStiff` (`:453`) / `getResistingForce` (`:682`) run on ~12 statics (`stiffJK`, `dd`, `BJ`, `BJtran`, `BK`, `BJtranD`, `shp`, `Shape` — `:534-544`, `:691`, `:718-721`); `ForceBeamColumn2d::update()` (`:1235`) runs its **entire interior Newton** on 14 statics declared inside the function (`dv` :1248, `vin` :1255, `vr` :1268, `f` :1269, `I` :1271, `dSe` :1281, `dvToDo` :1282, `dvTrial` :1283, `SeTrial` :1284, `kvTrial` :1285, `Ss`/`dSs`/`dvs`/`fb` :1342-1345).
- **Why:** the `Element` interface returns `const Matrix &`, so the conventional implementation needs storage outliving the call. A function-scope `static` is the cheapest way to get it, and upstream uses it everywhere. Counts above are **upper bounds on the hazard** — many sit in `sendSelf`/`recvSelf` and never run on a threaded path — but they size the audit honestly.
- **Failure mode if threaded anyway:** two elements in `getTangent` concurrently overwrite the same buffer, and the caller then does `addA(*eleTangent, ...)` on a buffer another thread is mid-write on ⇒ **silent, thread-count-dependent wrong answers**, no crash.
- **Workaround/status:** per-classTag **opt-in allowlist** (default empty) + ThreadSanitizer + bit-identical acceptance, per [[75b_ladruno_threaded_assembly_adr]] §5/§7. Prefer **per-element** buffers over `thread_local` statics: `thread_local` fixes the race but yields only `nthreads` buffers, which forecloses the exact gather assembly the determinism policy depends on (§4.2). *2026-07-25 (ADR-75b L3-0).*

### `ops_TheActiveElement` is a mutable GLOBAL written per-element inside `Domain::update()`'s loop and read by materials — a `static`-only re-entrancy audit misses it entirely
- **Bites:** you audit element/material kernels for `static` scratch, conclude the `update` loop is clean, thread it, and get plausible-but-wrong regularized softening. The hazard is not a `static` — it is a file-scope global: `Element *ops_TheActiveElement` (`SRC/element/Element.cpp:47`, `extern` in `SRC/G3Globals.h:45`), assigned **per element inside the loop** at `Domain.cpp:2401` (and in `Element`'s ctor `Element.cpp:65`, `Domain.cpp:461`, `OpenSeesCommands.cpp:2865`, `LadrunoDispBeamColumn2d.cpp:528`, `3d.cpp:651`).
- **Why:** it is a deliberate fork idiom — the "Phase-3b lch latch". Materials read it on first `setTrialStrain` to fetch a regularization characteristic length: `LadrunoJ2.cpp:353`, `LadrunoConcrete3D.cpp:347`, `ASDConcrete3DMaterial.cpp:1614`, `LadrunoRCConcrete.cpp:330`, `LadrunoRCFiniteStrain.cpp`, plus documented reliance in `BezierTet10.cpp`, `BezierTri6.cpp`, `LadrunoUP.cpp`.
- **Failure mode:** an element's material latches **another element's** characteristic length ⇒ a converged, plausible, wrong softening response. Silent.
- **Workaround/status:** must become `thread_local` (or be plumbed down the call chain) before `Domain::update()` is threaded — a prerequisite of [[75b_ladruno_threaded_assembly_adr]] L3-1, and it touches vanilla files (`Element.cpp`, `G3Globals.h`) so it owes [[LEDGER_vanilla_files]] rows when it lands. Note there are also several unrelated *definitions* of the same symbol in non-linked TUs (`NeesDataTest.cpp:45`, `TestDataOutput{Database,File,Stream}Handler.cpp:49`) — don't mistake those for the live one. *2026-07-25 (ADR-75b L3-0).*

### Profiler: a named deep scope (`elem.tangent`) times the WHOLE loop **including the `addA` scatter**, while its `elem_by_type` rows time ONLY the kernel call — reading the scope as "element cost" overstates the threadable fraction
- **Bites:** you read `formTangent` (or its `elem.tangent` child) as the element fraction, apply a ">40% element work ⇒ thread it" gate, and over-predict the win. `OPS_PROFILE_SCOPE_DEEP_NAMED(_ops_elemTan, "elem.tangent")` (`IncrementalIntegrator.cpp:109`) wraps the **entire** element loop — `getTangent` *and* `theSOE->addA` — whereas the per-classTag `elem_by_type` bucket, filled by `OPS_PROFILE_FE_ELEM_SCOPE` (`:117`), covers only the `getTangent` call. Threading buys the second, not the first.
- **So:** `scatter = scope_wall − Σ(elem_by_type wall)` is the non-threadable remainder, and it is not small. Measured (ADR-75b L3-0, `lane3/RESULTS_l3a_update_scope.md`): the `addA` share of the tangent loop is **19.1%** on Lane B under UmfPack, **9.5%** under PARDISO — and on Lane A's cheap `forceBeamColumn` tangents the **scatter (44.1 ms) EXCEEDS the kernel (35.4 ms)**, i.e. that assembly loop is scatter-bound, not element-bound.
- **Also:** a loop scope appears at **several places in the tree** — `elem.update` shows up under `newStep`, under `solveCurrentStep/update`, and (Lane A) directly under `solveCurrentStep` from `DisplacementControl`. Summing only one site undercounts by ~2× on Lane A and ~2× on Lane D. Same trap ADR-40b hit with the hidden second `soe.factor`.
- **Workaround/status:** use `Ladruno_files/testbed/perf/lane3/parse_lane3.py`, which sums every site and reports kernel-vs-scatter per loop. *2026-07-25 (ADR-75b L3-0).*

### `zeroA()` clears only the VALUES — the CSR pattern (`colA`/`rowStartA`) is built once in `setSize()` and survives every assembly, so OpenSees already satisfies "freeze the sparsity"
- **Bites (in the good direction):** a plan that budgets real work for building and maintaining a frozen sparsity graph before threaded assembly (the expensive half of the Kratos atomic-scatter pattern). It is already paid. In `PARDISOGenLinSOE`: `colA`/`rowStartA` are filled once from the DOF graph (`:221-297`); `zeroA()` zeros `A[]` and clears `factored` and touches neither; `addA` (`:370-478`) then locates each target by a **read-only** linear search over the frozen row and does exactly one read-modify-write, `A[k] += m(i,j)`.
- **Consequences:** (1) the entire implicit assembly race is *one* `+=` on a shared `double` at an index computed from immutable data — an atomic on one statement per SOE, no coloring, no thread-private matrices; (2) but the inner search makes the scatter **O(idSize² × rowlen)**, a pre-existing *serial* inefficiency worth fixing independently of threading (this is the "check is O(nnz) against an O(nnz × rowlen) loop" shape); (3) anything that can resize the SOE mid-run (`domainChanged` from ADR-51 element removal / ADR-60 contact re-emission) invalidates the freeze and must force a threaded phase **off**, not race with it — same family as the `-factorOnce` staleness caveat.
- **Workaround/status:** recorded as settled evidence in [[75b_ladruno_threaded_assembly_adr]] §2.2/§4.1. *2026-07-25 (ADR-75b L3-0).*

### `linearSolve` is NOT the solver cost when the integrator solves too — sum `soe.factor`/`soe.trisolve` over EVERY site
- **Bites.** You read the `linearSolve` phase as "time in the solver" and understate it. Under `DisplacementControl` the integrator calls `setB(phat); solve()` for the reference displacement `dUhat` in **both** `update()` and `newStep()` — a full extra factorization of the same `K`, booked under those phases, not under `linearSolve`. Measured on lane A (ADR-75b L3-0): `soe.factor` summed over all sites = **505.7 ms > `linearSolve` = 449.4 ms**, so true solver work is **7.86% of step, not the 5.65% `linearSolve` reports**. ADR-40b found the same thing on lane E at far greater severity (59% of step in factorization, two-thirds of it outside `linearSolve`) — and a later report still walked into it, which is why this is its own row.
- **Workaround/status:** to get true solver cost, sum `soe.factor + soe.trisolve` over **every** site in the tree, not the `linearSolve` phase. Affected integrators are any that solve inside `update`/`newStep` — `DisplacementControl` confirmed; `ArcLength`/`MinUnbalDispNorm` share the shape. *2026-07-25 (ADR-75b L3-0 adversarial review pass 2).*
- ✅ **The PARDISO half of this row is RESOLVED** — see the next row. It used to read "you try the same cross-check on a PARDISO run, get `soe.factor = 0.00%`, and conclude no factorization cost", because the scopes lived in `UmfpackGenLinSolver.cpp` only. PARDISO has had its own brackets since 2026-07-27.

### A profile recorded BEFORE 2026-07-27 shows `soe.factor = 0.00%` on a PARDISO run — that means "not instrumented", never "free"
- **Bites:** you re-read one of the archived ADR-75/75b profiles (e.g. `Ladruno_files/testbed/perf/lane3/l3a_laneB_pard*.json`, dated 2026-07-25) and conclude PARDISO's factorization is free, or that `linearSolve` is entirely triangular solve. Both wrong. Until [#667](https://github.com/nmorabowen/OpenSees/pull/667) (`62768d1f1`) `PARDISOGenLinSolver`/`PARDISOGenLinSOE` carried **zero** `OPS_PROFILE` scopes, so on a `system Pardiso` run the whole phase 11/22/23/33 cycle collapsed into one opaque `linearSolve` blob and every `soe.*` counter read zero.
- **Now:** 7 brackets exist — `soe.symbolic` (phase 11), `soe.factor` (22), `soe.trisolve` (33) using **UmfPack's exact names** so a cross-solver profile lines the phases up; `soe.cgs` (23) as its OWN bracket; `dc.s.fill`/`dc.s.verify` on the `setSize` CSR build; and a DEEP-gated `soe.addA`. CI: `tests/test_pardiso_solver.py::test_profiler_brackets_present`.
- **Two traps in reading the new brackets.** (1) **`soe.cgs` is not a subset of `soe.factor`.** When `-krylov`'s CGS gives up, PARDISO refactorizes *inside* the phase-23 call (Intel's automatic fallback), so that factorization bills to `soe.cgs` — a `-krylov` run's `soe.factor` therefore **understates** total factorization work. Use `iparm[19]` / `-stats` win-rate to interpret it. (2) **`soe.addA` is DEEP-gated**, so a coarse `profiler start` run shows it absent; that is the gating working, not a missing scatter cost.
- ✅ **The run-attribute half is FIXED too** (ADR-75 P1i, same day): `threads` was `threads_.size()` — profiler-*registered* threads, 1 on any single-threaded command layer regardless of `MKL_NUM_THREADS` — and `nElem`/`nNode` were promised by a comment and filled by nobody. **But any profile recorded before 2026-07-27 still carries `threads=1`**, so treat the attribute as unreliable on archived files and trust the harness log instead. `nnz` is still a uniform 0 by decision. *2026-07-27.*

### Two configs benchmarked in ONE process: whichever must GROW the heap pays first-touch page faults, so config order confounds the comparison at large N
- **Bites:** you loop configs inside one process to save build time, and read the difference as a property of the configs. Measured on ADR-75 P1j at 136,080 DOF: `-matrixType 2` took **310 s / 289 s running first** and **190 s running second** — a **1.58× swing from position alone**, with assembly −34% and solve −40%. It made symmetric assembly look **37% SLOWER** than unsymmetric at n=35 while being *faster* at n=15, i.e. it manufactured a plausible "the advantage inverts with size" story. Measured warm, symmetric assembly is faster at n=35 too (64.7 s vs 69.2 s) — **the artifact pointed the exact opposite way from the truth.**
- **Why the asymmetry:** unsymmetric CSR is ~2× the entries, so symmetric-after-unsymmetric is the only ordering that fits inside memory the process has already faulted in. Unsymmetric is position-indifferent (294 s first vs ~290 s second) because nothing before it ever allocated more. Rule of thumb: the penalty lands on whichever config **grows** the heap, which at fixed N is whichever runs first.
- **Workaround:** one configuration per process (what the P1h sweep did, which is why P1h was unaffected), or rotate the order and report position. **And prefer a WITHIN-RUN ratio to an absolute wall as the headline metric** — P1j's `fac/(fac+tri)` read 93.68 / 93.81 / 93.50% across both orderings and that 1.58× wall swing. It was chosen to be robust to this box's ±30% background-load noise and turned out to be robust to a much larger effect nobody anticipated. *2026-07-28 (ADR-75 P1j).*

### The Profiler is a process-global singleton and `ops.wipe()` does NOT reset it — a multi-run script without `profiler reset` reports the SUM of every run so far
- **Bites:** you write one script that loops over sizes / configs / repeats in a single process, `start`/`stop`/`report` around each, and get a table that is **clean, monotone and completely plausible** — and entirely wrong. Every report after the first contains all preceding runs. Measured on the ADR-75 P1j size sweep: `soe.factor` call counts came out **44, 88, 132, 176 … 936** (should be ~44 every run) and `step_ms` was a running total — 1952 s reported against a 348 s measured wall. Only run #1 was correct. Cost: one full 53-minute sweep, and the bogus trend was *in the same direction as the real one*, which is what makes it dangerous.
- **Why:** `ops_profiler::theProfiler()` is a process-global; `wipe()` tears down the Domain and touches nothing in the profiler. ADR-40c hit the same thing (its `soe.factor` call-count proof was written off as "inconclusive — in-process profiler accumulation") but it was never written down as a quirk, so it was rediscovered.
- **Workaround:** call `ops.profiler("reset")` immediately before every `profiler start` in any multi-run script. **And do not trust the profiler to police itself** — accumulation is invisible to every check that only reads the profiler, because every internal number stays self-consistent. The only reliable guard compares the profiler's `step` total against a clock it cannot influence: record `time.perf_counter()` around the same loop and assert agreement (`p1j_size_trend.py` writes `wall_by_run`; `p1j_rollup.py` refuses to print a table if any row is >25% off). One process per run also works and needs no guard.
- **Not applicable to** the single-run harnesses (`laneB_model.py`, the P1h sweep) — those fork a process per run, which is why P1h was unaffected. *2026-07-27 (ADR-75 P1j).*

### `nSteps=0` in a profile is NOT a bug — it means the run had no `-perStep`
- **Bites:** you see `nSteps=0` next to a healthy 15-step rollup, add it to a bug list, and go looking for a counter that was never broken. (Done — an early ADR-75 P1h draft published it as a defect alongside the two real `threads`/`nElem` bugs, and it had to be retracted.)
- **Why:** `nSteps` is derived from the per-step **series**, exactly like `dt_min`/`dt_max` — `buildMeta()` only fills it under `if (config_.perStep)`. A coarse run has no series to count, so 0 is the correct answer. The rollup's `root/step` scope still carries the true `calls` count if you need it.
- **Workaround/status:** by design. Read step count from `root/step` `calls`, or run with `-perStep`. Distinguish this from the genuinely-broken attributes in the row above, which were wrong *regardless* of flags. *2026-07-27 (ADR-75 P1h/P1i).*

### ADR-40b's lane-D "formTangent = 12.9 s (16.2%)" is STALE — ADR-67 P-NEW-1's constant-mass tangent cache is on by default and removed it
- **Bites:** re-running lane D expecting the ADR-40b explicit phase mix, seeing `formTangent ≈ 0.00%`, and hunting for a broken model or a lost scope. Nothing is broken: `CentralDifferenceLadruno` now ships `massCache = true` by default (`:89`, `:154`, `formTangent` override at `:259`) — the ADR-67 P-NEW-1 constant-mass tangent cache, i.e. "`-factorOnce` behaviour with safe invalidation", which is exactly the fix ADR-40b's Finding-3 item 1 recommended. Measured 2026-07-25: lane D `formTangent` **1.3 ms of a 29.5 s step** (was 12.9 s of 79.9 s), the rest of the mix shifting to `formUnbalance` 34.0% / `newStep` 32.0% / `update` 29.9%.
- **Workaround/status:** general rule for this fork — **a phase baseline in a dated report may have been optimized away by a later ADR.** Re-measure before quoting; check `LEDGER_implementations` for a shipped fix on that path first. *2026-07-25 (ADR-75b L3-0).*

### The Tcl `system Mumps` ladder DROPS unknown options SILENTLY — so a cluster BLR run measures full-rank and calls it a BLR result
- **Bites:** you follow the ADR-75 handoff's item 1 ("~10 min, no code"), run `system Mumps -ICNTL14 200 -stats` vs `... -BLR 1e-8 -stats` on the cluster, compare `INFOG(21)`, see no movement and no stats output, and conclude **"BLR is not the memory lever at our scale either."** Every token was discarded. `SRC/tcl/commands.cpp`'s `system Mumps` branch had **only** `-ICNTL14`/`-ICNTL7`/`-matrixType`; its terminal `else currentArg++` arm skipped everything else **without a warning**, and it constructed `MumpsParallelSolver(icntl7, icntl14)` — the 2-arg form, so `ICNTL35`/`CNTL7`/`printStats` could not be set even in principle. P2/P2b wired `-BLR`/`-stats` into `SRC/interpreter/OpenSeesCommands.cpp` only.
- **Why this instance was worse than P1d's and P1e's.** Those were *availability* gaps: `system Pardiso` simply did not exist in Tcl, so a Tcl user got a loud error. This one is a **silent-wrong-configuration** gap, and it sat on the **only interpreter the cluster uses** — `OpenSeesMP` is Tcl-driven (apeGmsh-emitted partitioned decks auto-emit `system Mumps`), while the Python `openseesmp` module that *did* have the wiring has never been built there. So the one measurement ADR-75 still owed was, by construction, unobtainable and would have returned a **confident false negative**.
- **The general rule, restated with the sharper edge:** the two `system` if-ladders drift (banked at P1d, re-banked at P1e) — **and a ladder that ignores unknown options silently converts that drift into fabricated measurements.** When a "just run it" item depends on an option, *prove the option is parsed on the interpreter you will actually run*, not on the one where it was implemented. Check the invariant where you depend on it.
- **Workaround/status:** ✅ fixed 2026-07-25 (ADR-75 P2h) — `-BLR`/`-stats`/`-ICNTL35`/`-CNTL7` wired into the Tcl ladder mirroring the Python one, the 5-arg ctor used on both parallel arms, a missing-value bounds check added (the ladder also read `argv[currentArg+1]` past the end for a trailing value-taking option), and the silent skip replaced by `Mumps Warning: unrecognized option '<opt>' -- ignored`. Unknown options remain non-fatal. **Compile-verified on all three `#ifdef` arms; NOT yet executed** — a Tcl+MUMPS binary needs the deferred cluster rebuild, so the runtime parity smoke is owed at the item-1 run. *2026-07-25 (ADR-75 P2h).*

### esmeralda: SLURM is at `/opt/slurm/bin` and is NOT on `PATH` — this has now cost two different things
- **Bites (1), the big one:** `which sbatch` / `which sinfo` return nothing, which reads as "there is no scheduler / the cluster is down". It is not. `esmeralda` had **33 days of uptime** while an entire ADR-75 session recorded the cluster as down and deferred both cluster-gated items on that basis. `/opt/slurm/bin/sinfo` answers instantly: 18 nodes, 32 cores / 60 GB each.
- **Bites (2):** a `#SBATCH` script that calls bare `srun` dies with **`rc=127` (command not found)** and an **empty log**. Two consecutive sweep submissions (jobs 144449, 144451) produced zero results this way. It is easy to misread as an MPI/launcher problem, because the *previous* failure on the same script genuinely was one (see the next row).
- **Workaround/status:** call `/opt/slurm/bin/{sbatch,srun,squeue,sinfo,sacct}` by absolute path from scripts, or prepend it to `PATH` at the top of every batch script. *2026-07-26 (ADR-75 P2h).*

### `mpirun` works inside a 1-node SLURM allocation and dies instantly across 2 nodes
- **Bites:** the identical sweep script runs fine on `--nodes=1` and, on `--nodes=2`, every launch fails immediately with `[[...]] FORCE-TERMINATE AT (null):1 - error plm_slurm_module.c(471)` / `An internal error has occurred in ORTE`. OpenMPI's ORTE SLURM launcher, not the model.
- **Workaround/status:** use `/opt/slurm/bin/srun --cpu-bind=cores --mpi=pmix_v3 <wrapper> deck.tcl` for anything multi-node — the launcher `02_esmeralda_linux_build_guide.md` §7 already documents. `srun --mpi=list` confirms `pmix_v3` is available. *2026-07-26 (ADR-75 P2h).*

### A sweep script that pipes a run straight into `grep` can hide the reason it failed — including from itself
- **Bites:** the ADR-75 P2h sweep piped each run's output into `grep -E "P2H_RESULT|...|ERROR|rror"`. When `srun` failed with `srun: fatal: ...` / `command not found`, **none of the filter's patterns matched**, so the job produced a clean-looking log, printed its `..._DONE` banner and exited 0 — with zero results and zero explanation. Two takes were burned before the cause was visible. Compounding it: `mpiexec`/`mpirun`/`srun` wrappers **return 0 even when every rank died**, so `rc` proved nothing either.
- **Workaround/status:** write each run's FULL output to its own file, then `grep` the **file**; and gate on an **artifact** (`if grep -q P2H_RESULT "$LOG"`), dumping the log head when the artifact is missing. This is the banked "never let a bench script grep away its own log" rule — it was violated in the very sweep meant to honour it, which is why it is re-banked here with the concrete symptom. *2026-07-26 (ADR-75 P2h).*

### Never rebuild the OpenSees binary while a sweep is running
- **Bites:** a multi-mode sweep re-execs the binary once per mode (`openseesmp.sh` → `exec`). Rebuilding mid-job silently swaps the executable between modes, so the A/B comparison spans two different binaries and nothing in the output says so.
- **Workaround/status:** hold rebuilds until `squeue` is clear, or build to a distinct path and point the wrapper at it explicitly. *2026-07-26 (ADR-75 P2h).*
### MUMPS `ICNTL(7)=7` (auto ordering) picks PORD and PORD dies on a DENSE-as-CSR pencil
- **Bites:** you hand MUMPS a matrix that is *stored* sparse but *structurally* dense — which is exactly what a Craig-Bampton reduction produces — and the ANALYSIS phase aborts with `Error in function orderMinPriority / no valid number of stages in multisector (#stages = N)`, having burned enormous memory (15.6 GiB on Building 1A) before any factorization happened. Because it dies in analysis, it reads like an out-of-memory or a corrupt-matrix bug, and because it is threshold-driven it looks *rank-count*-dependent: ADR-1000 saw it at np=2 while np=4 and np=6 passed, and the resulting diagnosis ("MUMPS chose *parallel* analysis; force `ICNTL(28)=1`") was wrong.
- **Why:** `ICNTL(7)=7` lets MUMPS choose the ordering, and for this pattern it chooses PORD (`ICNTL(7)=4`), whose multisector minimum-priority ordering cannot form valid stages on a dense pattern. `ICNTL(28)` (sequential vs parallel analysis) is irrelevant — measured: `ICNTL(28)=1` with `ICNTL(7)=7` still fails identically. Rank count is also irrelevant: with a dense order-12000 pencil the failure reproduces at **both np=2 and np=4**, while a *sparse* order-64000 pencil at the same rank counts is fine. What varies with np on a real model is the size/density of the rank-local pencil, i.e. whether you cross the threshold — not the process count itself.
- **Workaround/status:** pin `ICNTL(7)=0` (AMD — no multisector concept, always compiled into MUMPS). Done at **both** `LadrunoCMSMumps.cpp` factorization sites: the distributed one (`MPI_COMM_WORLD`, `ICNTL(18)=3`) and, more importantly, the serial `MumpsSPD` one (`MPI_COMM_SELF`) that the rank-local CB pencils actually go through — the serial site fails the same way on one rank. Reproducer + regression guard: `testDistributedMumpsAtScale` / `testSerialMumpsAtScale` in `tests/ladruno_cms_mumps_check.cpp`, sized by `LADRUNO_CMS_CHECK_DENSE_ORDER` (defaults are CI-cheap; set it to 12000 to reproduce the failure on the pre-fix code). *2026-07-26 (ADR-1000 Part 0).*

### A distributed MPI check that guards its collective leg on `size == N` silently tests NOTHING at other sizes
- **Bites:** ADR-1000's P3 plan prescribed "run the checks at np=2, they must now pass" as the acceptance gate for a 2-rank MUMPS fix. Both distributed legs (`testDistributedMumps`, `checkDistributedFourRankFlow`) opened with `if (size != 4) return;`, so at np=2 they returned instantly and the check printed `passed` having exercised no distributed code at all. A green np=2 run would have "validated" the fix while proving nothing.
- **Workaround/status:** when a check's coverage depends on `MPI_Comm_size`, either make the fixture size-generic or make the skip **loud** (print what was skipped and why). The tiny 2x2 distributed fixture was size-generic all along — only the guard was wrong. Note the second-order trap: a fixture can be the right *size* and still be the wrong *shape* — the 2x2 one runs at np=2 but is far too small for MUMPS to make an ordering decision, so it could never have caught the bug the gate was written for. *2026-07-26 (ADR-1000 Part 0).*
### `algorithm Newton -initial` re-assembles AND re-factorizes the initial tangent on EVERY iteration — `ModifiedNewton -initial` is the same algorithm at 1/2 the cost per step ONLY when the initial tangent is genuinely state-independent
- **Bites:** anyone reaching for `Newton -initial` as *the* robust fallback (standard practice in pushover work). `NewtonRaphson::solveCurrentStep` puts `formTangent` **inside** the iteration `do`-loop (`NewtonRaphson::solveCurrentStep`, the `do`-loop) and special-cases only `INITIAL_THEN_CURRENT_TANGENT`; `INITIAL_TANGENT` falls through to the generic branch and is re-formed unconditionally. Each re-formation is a full `zeroA()` + `FE_Element` assembly loop **and** it clears the SOE `factored` flag, forcing the solver to redo its numeric factorization. Measured by the TIMs reporter on a 38 984-DOF `SSPbrickUP`/`ShellMITC4` model with MKL PARDISO/8 threads: `Newton -initial` **1.504 s/step** vs `ModifiedNewton -initial` **0.936 s/step** (**1.61x**) for an IDENTICAL converged settlement in an IDENTICAL iteration count (2.00 iters/step both). Unit costs: assembly 325 ms, numeric factorization 255 ms, triangular substitution 32 ms. Same ratio reproduced on an earlier UMFPACK build (6.56 vs 3.88 s/step).
- **Why:** `ModifiedNewton::solveCurrentStep` forms the tangent **once before** its `do`-loop. **Conditionally** the same iteration: when the assembled initial tangent really is state-independent, `Newton -initial` is simply the expensive spelling of `ModifiedNewton -initial`. When it is NOT — corot3d, contact, any `-damping` element, or any transient run with `betaK`/`betaK0` != 0 — they are genuinely different algorithms (Newton re-linearizes on the current configuration each iteration; ModifiedNewton freezes at step start), and iteration counts and convergence can differ. ⚠ **And the usual "`K_f` only sets the path, never the equilibrium" reassurance is FALSE under a displacement- or energy-based convergence test**: `NormDispIncr`/`RelativeNormDispIncr`/`EnergyIncr` accept on `dU = K_f^-1 R`, so a different `K_f` accepts at a *different point* — a stiffer-than-appropriate `K_f` gives small `dU` with non-small `R`, the classic false convergence. It also changes whether a step fits inside `maxIter`, which in any adaptive-substepping driver changes the load path. The path-only claim holds only for a force-residual test, conditional on convergence. Every sparse solver in the fork already gates its numeric phase on `theSOE->factored` (`PARDISOGenLinSolver.cpp:378`, `UmfpackGenLinSolver`, `MumpsSolver`); the flag is cleared *because* the matrix is re-assembled, so the saving is there for the taking and nothing upstream takes it.
- **The trap inside the trap — "static ⇒ the matrix is invariant" is NOT sound.** `addKiToTang()` reaches `Element::getInitialStiff()`, and in OpenSees that name is a convention, not a contract. **Confirmed configuration-dependent:** `CorotCrdTransf3d::getInitialGlobalStiffMatrix` (`:1452`) triple-products through the member `T`, recomputed by `compTransfMatrixBasicGlobal()` at the tail of `update()` (`:549`); `LadrunoContactFE::addKiToTang` (`:1590`) mirrors `addKtToTang` including the augmented active mask; `updateMaterialStage` changes `getInitialTangent()` **without** calling `Domain::domainChange()` (`MaterialStageParameter` never does), so it is invisible to the topology stamp. **Confirmed clean** (don't re-derive): `CorotCrdTransf2d::getInitialGlobalStiffMatrix` uses only `cosTheta`/`sinTheta`/`L` from `initialize()` (`:318`), never the current `cosAlpha`/`sinAlpha`; `PDeltaCrdTransf3d` builds `T_bl` from `L` alone; `CorotTruss::getInitialStiff` fixes `R` in `setDomain` (`:512`). So invariance is an **element/transformation** property, not an integrator-family one.
- **Transient extra:** under `Newmark` the `INITIAL_TANGENT` matrix is `c1*Ki + c2*C + c3*M` (`Newmark.cpp:295-298`) and `C = alphaM*M + betaK*Kt + ...` (`Element.cpp:222`), so with `betaK != 0` the matrix genuinely DOES change between iterations and the re-formation is legitimate work. Related surprise worth knowing: with `rayleigh $a $b 0 0` the `-initial` flag is very nearly **inert** — in the reporter's model (`Newmark(0.5,0.25)`, `dt=1e-3`, `zeta=0.5`, T=0.5/0.2 s) `c1=1`, `c2=2000`, `b=0.02274`, giving `1.0*Ki + 45.5*Kt + 4.02e6*M` under INITIAL vs `46.5*Kt + 4.02e6*M` under CURRENT. `-initial` exchanges ~**2%** of the stiffness content there; the rest re-enters through the damping. ⚠ **That 2% does NOT generalise** — it needs BOTH a large zeta and a small dt. Same periods at zeta=0.05, dt=0.01: `b=0.00227`, `c2=200`, `c2*b=0.455`, so `-initial` exchanges ~**69%** and the flag is doing very nearly what the user expects. The lesson is to *compute `c2*betaK` for your own model* before assuming `-initial` is either effective or inert; at high damping ratios with fine time steps it is close to inert, and you pay full price for it.
- **Workaround/status (2026-07-25, ADR-76 from a TIMs project issue report):** use **`algorithm ModifiedNewton -initial`** (one assembly + one factorization per step), or **`-initial -factoronce`** for one per *analysis* — that pair only started working in this fork's `OPS_ModifiedNewton` as of ADR-76 R4 (see next quirk). An engine-side tangent-version counter was designed and then **WITHDRAWN** after adversarial review — see [[76_ladruno_tangent_reuse_adr]] §4 for why (its invalidator set was incomplete, and the flagship measurement is on a model where it must never fire). Appendix A of that ADR is the superseded design, kept only as a record; do not implement from it.

### `OPS_ModifiedNewton` read exactly ONE option — `algorithm ModifiedNewton -initial -factoronce` silently dropped `-factoronce` (fixed ADR-76); and `factorOnce` still has no `domainChanged` reset
- **Bites:** `algorithm ModifiedNewton -initial -factoronce` — the natural spelling of "initial stiffness, assembled and factorized once", and for a static analysis arguably the cheapest robust algorithm the framework offers. Upstream `OPS_ModifiedNewton` opened with `if (OPS_GetNumRemainingInputArgs() > 0) { const char* type = OPS_GetString(); ... }` — a single read feeding an `else if` chain — so the *first* option won and every later one was silently ignored. No warning; the deck runs, just not the algorithm you asked for. Whichever option you put first is the one you get.
- **Why:** an `if` where the sibling factory `OPS_NewtonRaphsonAlgorithm()` (same file family) has always used `while (OPS_GetNumRemainingInputArgs() > 0)`. Pure oversight, invisible because the single-option spellings all work.
- **Fixed (2026-07-25, ADR-76 R4):** `if` → `while` in `OPS_ModifiedNewton`. **Adversarial review then found four defects in that one-line cut, all fixed in the same commit:** (1) `-factorOnce` (camelCase — what `Linear.cpp:67`, `ExpressNewton.cpp:78` and both `commands.cpp` sites accept, and what THIS ledger writes) was still being silently dropped by the very fix meant to stop options being dropped; now accepted, with `-Initial`/`-Secant` for NewtonRaphson parity. (2) a failed `-hall` factor read used to `return 0`, and openseespy's `OPS_Algorithm` used to silently discard a null and report SUCCESS, leaving the PREVIOUS algorithm in force; the parser now warns and keeps the 0.1/0.9 defaults, AND the choke point itself now errors on any null factory result (`OpenSeesCommands.cpp:2060` — see the `OpenSeesCommands.cpp` ADR-76 row in [[LEDGER_vanilla_files]]). (3) `-secant`/`-initial` now reset `iFactor`/`cFactor` as `OPS_NewtonRaphsonAlgorithm` does — unnecessary under an `if`, required once a later option can override an earlier one. (4) **unknown tokens now WARN** instead of vanishing. Deliberately preserved: `-hall`'s trailing-factor read still gates on `OPS_GetNumRemainingInputArgs() == 2` (so `-hall a b` must END the command for the factors to be read — the identical quirk `NewtonRaphson` carries; widening it is a separate decision). `OPS_ModifiedNewton` is the single choke point for BOTH interpreters — classic Tcl `specifyAlgorithm` (`SRC/tcl/commands.cpp:4467`) and openseespy (`SRC/interpreter/OpenSeesCommands.cpp:2017`) both call it. **Note also:** `SRC/runtime/commands/analysis/algorithm.cpp` parses `algorithm` a *third* time and drops `factorOnce` entirely — but that whole tree appears in no `add_subdirectory()` (`SRC/CMakeLists.txt:8-30`) and is dead in this build. Don't "fix" it and don't trust it when grepping for parser behaviour.
- **Trap the fix makes newly reachable:** `factorOnce` has **no `domainChanged` reset**. `SolutionAlgorithm::domainChanged()` is virtual but `ModifiedNewton` does not override it, and the 1→2 latch only ever resets on a convergence *failure* (`ModifiedNewton::solveCurrentStep`, the `result == -2` branch). After a mid-run domain change the SOE is re-sized/zeroed but the tangent re-form is skipped ⇒ solve against a stale matrix. Never combine `-factoronce` with element removal (ADR-51), contact re-emission (ADR-60), or staged construction. This is the SAME trap already recorded above for `algorithm Linear -factorOnce`; R4 just adds a second spelling that can reach it. The designed fix is the `LadrunoModifiedNewton` fork class of [[76_ladruno_tangent_reuse_adr]] §4.5 (re-arm on `domainChanged()` + every negative return; covers only the topology subset of invalidators, so this warning survives it); unstarted.

### `NDMaterial::getInitialTangent()` DEFAULTS to `getTangent()` — so on most solid models `algorithm Newton -initial` is silently full Newton, at full cost, and "initial stiffness" is a naming convention rather than a contract
- **Bites:** anyone using `-initial` (or `ModifiedNewton -initial`, or `-intialThenCurrent`, or Hall) as a robustness fallback on a **solid/continuum** model and wondering why it neither behaves nor converges like initial-stiffness iteration. `SRC/material/nD/NDMaterial.h:64` is literally `virtual const Matrix &getInitialTangent(void) {return this->getTangent();};` — the base-class default for the *initial* tangent is *the current* tangent. Any nD material that does not override it makes `addKiToTang()` and `addKtToTang()` assemble the **same matrix**, so `-initial` degenerates to `-current` with no diagnostic. You pay the full assembly + factorization for an iteration strategy you are not getting.
- **Why the asymmetry with frames:** `UniaxialMaterial::getInitialTangent` (`SRC/material/uniaxial/UniaxialMaterial.h:68`) is **pure virtual** (`= 0`) — uniaxial materials are *forced* to implement it. That is why fiber/truss/frame models generally do get a real initial stiffness and solid models generally do not. The two base classes made opposite choices and nothing documents it.
- **It is not only the materials.** ~41 element `getInitialStiff()` implementations in `SRC/element/` are a straight `return getTangentStiff();` — including the whole `SSPquad`/`SSPbrick`/`SSPquadUP`/`SSPbrickUP` family. `Joint2D::getInitialStiff` (`SRC/element/joint/Joint2D.cpp:1135`) calls `theSprings[i]->getTangent()`, not `getInitialTangent()`; `Joint3D.cpp:573` likewise. The `zeroLength` contact family and the `UWelements` contact set branch on the current active set (`ZeroLengthContact3D::getInitialStiff:360` literally calls `formResidAndTangent(tang_flag=1)`). **Correction (fact-check):** `ASDShellQ4`/`ASDShellT3` were wrongly grouped here — they pass `OPT_LHS_IS_INITIAL` and genuinely select `getInitialTangent()` on that flag, i.e. they have a real initial path. Materials that *do* override can still return a state-updated member: `stressDensity` recomputes `initialTangent` from the current trial mean stress; `ManzariDafalias` returns `mCe`, which the integrator rewrites on every `setTrialStrain`; `PM4Sand`/`PM4Silt` rebuild `mCe` in `commitState()`.
- **Separate latent bug found alongside:** `CorotCrdTransf3d`'s transformation matrix `T` is declared **`static`** (`SRC/coordTransformation/CorotCrdTransf3d.h:135`), i.e. shared across every instance, and `DispBeamColumn3d::getInitialStiff` never calls `crdTransf->update()`. So it triple-products through whatever `T` the *last element to call `update()`* left behind — cross-element aliasing, not merely configuration dependence. Worth reporting upstream; not fork-introduced.
- **Workaround/status (2026-07-25, ADR-76 adversarial review):** treat "does this element/material actually have an initial tangent?" as a **per-class question you must check**, never as a property of the `-initial` flag. Practically: on a solid model, if `-initial` changes neither the iteration count nor the convergence character versus the default, that is the expected result, not a bug in your deck. This is also why [[76_ladruno_tangent_reuse_adr]] Appendix A.4 makes `Element::isInitialStiffInvariant()` default **false** — an earlier draft defaulted it true and would have silently opted in most of `SRC/element/`.

### A first-pivot singularity on `system BandGeneral` / `FullGeneral` / `BandSPD` returned SUCCESS — `analyze` gave 0 and the "displacement" was the load vector (fixed; but the assumption "BandGeneral fails loudly on a singular matrix" was never true)
- **Bites:** anyone relying on `analyze` returning nonzero to detect a singular tangent — i.e. every scripted driver, every adaptive-substepping loop, every automated study. Three LAPACK-backed solvers mapped a positive LAPACK `info` with `return -info+1;`, which C parses as `(-info)+1`. LAPACK sets `info = i` when `U(i,i)` is exactly zero, so `info == 1` — a singularity at the **first** pivot, exactly what an all-zero assembled `A` produces — returned **0 == SUCCESS**. A WARNING went to stderr, but the return code lied, and in a long run that warning just scrolls past.
- **Why it is a silent WRONG ANSWER, not a missing error:** the solvers copy `B` into `X` **before** calling LAPACK (`BandGenLinLapackSolver.cpp:116-118`), and DGBSV/DGESV/DPBSV do **not** compute `X` when `info > 0`. So `X` is left equal to `B`, the caller (e.g. `ModifiedNewton::solveCurrentStep`, gating on `theSOE->solve() < 0`) accepts it, and the analysis proceeds using the **load vector as a displacement increment**. Compounding it, `theSOE->factored = true` sits *after* the error return (`:166` post-fix), so the flag stays false and every subsequent iteration silently repeats the identical non-solve. **Measured** (3 free DOF, `uniaxialMaterial Elastic` with E=0 so A is identically zero, load 7.0, `algorithm Linear`): `analyze` returned **0** with `nodeDisp == 7.0` on BandGeneral, FullGeneral **and** BandSPD.
- **⚠ It fires ONLY at `info == 1`, and that is why it survived for decades.** `(-info)+1` is already negative for every `info >= 2`, so those cases were detected correctly. `info == 1` means the singularity is at the **first pivot**. Measured on the same 4-node chain: a free node numbered **last** gives `info = 3` -> old return `-2` -> **correctly detected**; the same free node numbered **first** gives `info = 1` -> old return `0` -> **bug**. So the "obvious" singular model (a forgotten unfixed node) reproduces this only when its zero row happens to land on the first pivot. Anyone spot-checking with a casual singular deck will most likely conclude there is no bug.
- **`ProfileSPD` was always correct** — `ProfileSPDLinDirectSolver` hard-returns `-2` on a non-positive first diagonal. That asymmetry is why a model can "work" on one `system` and fail on another for reasons that have nothing to do with the solver's numerics: **if you need a singular tangent to fail loudly, `ProfileSPD` was the only one of these four that did.**
- **⚠ Trap when writing a repro — a minimal test can be TOO minimal.** `BandSPDLinLapackSolver.cpp:98` and `ProfileSPDLinDirectSolver.cpp:150` carry 1x1 special-case guards (`singular 1x1 system` / `singular 1x1 (|aii| < 1e-15)`) that intercept a single-DOF singular system **before** the LAPACK call. Our first probe used 1 DOF and reported BandSPD as CLEAN; at n=3 it fails with the other two. Use n >= 2 for anything in this area.
- **Workaround/status (2026-07-25, spun out of the ADR-76 `factorOnce`/`domainChanged` audit):** **FIXED** — `return -info;` in all three (always negative in that branch, and it keeps the failing pivot index in the magnitude; the idiom already used at `SymBandEigenSolver.cpp:228` / `FullGenEigenSolver.cpp:188`). No caller reads the magnitude — I found no comparison against a specific value anywhere in `SRC/`. **Correction (fact-check):** an earlier wording of this row said "every call site gates on `< 0`", which is FALSE — dozens of call sites ignore the return code entirely (`ArcLength.cpp:244,1053,1064`, `DisplacementControl.cpp:321,1092,1102`, `MinUnbalDispNorm`, `EQPath`, `LadrunoArcLength.cpp:379`, `Newmark.cpp:1047`, `SensitivityAlgorithm.cpp:133`, the accelerator family, ...). The load-bearing half holds, so the fix is safe — but it also means the BUG was worse than first documented: on an arc-length or displacement-control step the singular-solve return was never inspected at all. Regression deck: `Ladruno_implementation/lapack_singular_regression/`. **Expect fallout:** models that previously limped past a singular step now stop. That is correct, but it will surface latent singularities (free nodes under `constraints Plain`, zero-stiffness materials, fully-released members) that were being masked. See [[LEDGER_vanilla_files]]. **CI (2026-07-26):** the regression deck runs in the Zone-A job, gated on its terminal marker (never the exit code — see the parse-error quirk below).
- **The predicted fallout landed in our own suite first (2026-07-26):** the first post-merge Zone-A run on ladruno failed **7 zone_a tests that had been relying on the limp-past**. Six were static rigs whose FIRST tangent is legitimately singular — a node/facet held by NOTHING until contact engages (`test_adr39_contact_p4_soft`, `p5_soft2`, `test_adr41_viscous_d2`, `soft2_visc`, `test_contact_review_p5_percontact`) or a fully-softened cohesive hinge with tangent exactly 0 (`test_ladrunoCohesiveHinge_material`) — and the SWALLOWED singular solve (X = B) was the very thing that seated the contact / limped past exhaustion. Fixed by grounding the free DOFs with a spring ~7 orders softer than the penalty (analytic asserts hold to ~1e-10) or stopping the push at 0.999·kappaf. The 7th (`test_adr30_projection_p0` massless-DOF) ASSERTED the old swallow as its premise — its own assertion message said "§2.5 premise drifted"; rewritten to assert the refusal and guard against the defect regressing. **Moral, twice over:** (a) a model that "converges" on its first static step while held only by contact is probably riding a swallowed singular solve — post-fix such rigs need seating or grounding; (b) the fallout was invisible for a day because the ADR-76 session was STRANDED with no PR and the rewrite PR auto-merged before Zone-A finished — auto-merge + non-required checks means CI gates nothing.
- **Two dormant threaded solvers NOT fixed (recorded so a resurrection inherits the warning, not the bug):** `ThreadedSuperLU.cpp:120` is a true same-shape copy — it forwards **pdgstrf**'s positive zero-pivot `info` via `return info;`, so callers gating on `< 0` see success on a singular factor. `BandSPDLinThreadSolver` is worse in a different way: its factorization runs in worker threads and lands `info` in `TCB_BandSPDLinThreadSolver.info` (`:158`), which is **never checked** after the `threadsDone` wait — a singular factor sails straight into `dpbtrs_` (whose raw `info` is forwarded at `:189`, though `dpbtrs` itself only reports negative/illegal-argument codes). Neither file is compiled by any CMake target (only the legacy `SRC/Makefile` `PROGRAMMING_MODE=THREADS` build references ThreadedSuperLU), so both defects are consequence-free today — fix them the day either threaded solver is brought back into a build.

### `OpenSees.exe` exits 0 after a Tcl **parse** error — a deck that never ran looks like a deck that passed
- **Bites:** any CI/harness that gates on the process exit status. A brace mismatch (or any Tcl syntax error) aborts the script with `missing close-brace` on stderr, the deck's own `exit 1` on the failure path is never reached, and the process still exits **0**. Hit live while extending `Ladruno_implementation/lapack_singular_regression/` — the run printed its banner, printed nothing else, and reported success. A harness would have recorded a green run for a deck that executed zero assertions.
- **Why:** the Tcl interpreter reports the error and returns; `tclMain` does not translate a script error into a nonzero process status. Same family as the two exit-status traps already recorded here (`analyze()` returning 0 on a NaN field; `mpiexec` returning 0 when every rank died).
- **Workaround/status (2026-07-25):** never gate solely on the exit code. Gate on a **positive terminal marker** the deck prints only on the success path (e.g. grep for a final `=== ... all checks passed ===` line), or count the expected number of PASS lines. Both the ADR-76 smoke and the LAPACK regression print such a marker; the checker should require it, not merely tolerate its absence.

### `require(call(...), "..." + message)` reports an EMPTY diagnostic — argument evaluation order is unspecified
- **Bites:** every CMS standalone check uses `require(bool, std::string)` and the natural idiom `require(doThing(..., message) == 0, "doThing failed: " + message)`. C++ does not specify the order in which function arguments are evaluated, so MSVC builds the message string **before** running the call — capturing `message` while it is still empty. The failure prints `FAIL: doThing failed:` with nothing after the colon, which reads like the callee returned no diagnostic and sends you looking in the wrong place. Cost the P3d work a full debug cycle: the real message was `invalid distributed hierarchy input on at least one rank`, which points straight at the cause.
- **Why it is not a correctness bug:** the *condition* is still evaluated correctly, so nothing passes that should fail. Only the diagnostic is lost — and only on the failure path, which is exactly when you need it.
- **Workaround/status:** call first, store the bool, then `require(ok, "... " + message)`. Swept 2026-07-26: `REQUIRE_CALL(status, text)` evaluates the call first and only then builds the diagnostic; 11 sites converted across `assembly` (4), `lanczos` (4), `mumps` (2) and `topology` (1). Sites that already computed the bool into a variable were left alone -- they never had the problem. **Use `REQUIRE_CALL` for any new `require(someCall(...), "..." + message)`.** *2026-07-26 (ADR-1000 P3d).*

### `DistributedHierarchyInput::fine` must equal the MPI rank — it is not a free label
- **Bites:** you try to express "rank r owns partition p" by setting `input.fine = p`, and `solveDistributedHierarchy` rejects the whole collective with `invalid distributed hierarchy input on at least one rank`. The validation is `input.fine != rank` (`LadrunoCMSHierarchy.cpp:1293`).
- **Consequence for testing:** a rank/partition permutation can only be expressed by **moving the data** (which subdomain's equations/stiffness/mass a rank carries), keeping `fine = rank`. That is also the physically honest formulation — it is what a different partitioner would hand you.
- **Workaround/status:** by design, documented here so the next agent does not read the rejection as a bug. *2026-07-26 (ADR-1000 P3d).*

### The Tcl `nDMaterial` command is a hand-written strcmp ladder — the fork's own nD materials are openseespy-ONLY
- **Bites:** `nDMaterial LadrunoJ2 ...` from a `.tcl` deck under `OpenSees`/`OpenSeesMP` fails, while the identical call works from openseespy. `TclModelBuilderNDMaterialCommand` (`SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp`) dispatches by an explicit `strcmp(argv[1], ...)` chain with **no generic `OPS_*` fallback**, and **no `Ladruno*` nD material appears in it** — the fork's materials were only ever registered on the Python side (`SRC/interpreter/OpenSeesNDMaterialCommands.cpp`).
- **Contrast with elements, which are fine.** `SRC/element/TclElementCommands.cpp` *does* carry a dispatch TABLE with the fork's entries (`{"LadrunoBrick", "ladrunoBrick", OPS_LadrunoBrick}` at :596, `LadrunoBrick20` at :603, …), so `element LadrunoBrick` works from Tcl. **Elements were wired, materials were not** — do not assume that because one fork object reaches Tcl, its neighbours do.
- **Same family as the `system Mumps` gap** (ADR-75 P2h) and the `system Pardiso` gap (P1d) and the SMS-integrator gap (#340): this fork keeps discovering that the Tcl and Python command surfaces are independent and drift. This instance is the widest one found so far — it is *every* fork nD material, not one verb.
- **Consequence, concretely:** it blocked ADR-75b's G-L3 gate from using the fork's own `LadrunoJ2`, forcing vanilla `stdBrick`+`J2Plasticity` as a proxy. That happened to be survivable (the gate failed by ~42×, far outside the proxy's error bar) but it would not be survivable for a measurement with a narrower margin.
- **Workaround/status:** ✅ **FIXED 2026-07-26** — all **10 fork factories / 15 names** (including every alias: `LogStrain`, `LogStrain2D`, `LadrunoCohesiveHingeBiaxial`, `LadrunoJ2`, `LadrunoJ2Finite`, `InitDefGrad`/`StagedDefGrad` ×4, `StagedStrain` ×2, `LadrunoRCConcrete`, `LadrunoRCFiniteStrain`, `LadrunoConcrete3D`) wired into the Tcl ladder, mirroring the Python map one-for-one. Purely additive. Verified by compile + negative control, by a real link on the cluster, and at RUNTIME (`nDMaterial LadrunoJ2` + `element LadrunoBrick` converge from a `.tcl` deck under `OpenSeesMP`). **The lesson survives the fix:** the two command surfaces are independent and *silently* diverge — when a fork object works from one interpreter, that is no evidence at all about the other. *2026-07-26 (ADR-75b G-L3).*

### `recvSelf` into a LIVE element can flip "construction-fixed" inputs — any per-instance cache keyed on the mutable-only signature must be invalidated there
- **Bites:** the mass-cache contract (`SRC/element/LadrunoMassCache.h`) lets elements omit construction-fixed inputs (`massType`, `quadz`, `formulation`, `cMass`) from the guard signature — correct for every normal flow, because broker-built receives start with an empty cache. But `recvSelf` on an *already-live* element (database `restore`, object reuse) rewrites exactly those members from the stream; if rho and coords happen to match, the next `getMass()` is a clean-guard HIT serving the pre-recv matrix — a consistent mass on a now-lumped element, Gauss-rule mass on a now-Lobatto element. The A/B on/off byte-equality battery can never see it (it never re-receives into a live object). Found by the ADR-77 review wave in ALL six cache-bearing elements, the LadrunoBrick donor included — a pattern-level gap, not a transcription error; the tell was that SolidShell's `recvSelf` already invalidated `Ki` for this exact reason while leaving `massCache` alone.
- **Rule:** whatever a cache's guards exempt as immutable, `recvSelf` must `invalidate()` — the two claims ("fixed at construction" and "recvSelf may rewrite it") are both true and only compatible if the cache dies on receive.
- **Workaround/status:** ✅ FIXED — `massCache.invalidate()` (brick: `delete Mi`) added at the end of all six `recvSelf`s. *2026-07-27 (ADR-77 review wave).*

### The bare `this->getMass();` side-effect idiom still lives in the cache-LESS `LadrunoCST`/`LadrunoCSTPair` — four pre-planted traps for the next mass-cache extension
- **Bites:** `LadrunoCST.cpp` (`addInertiaLoadToUnbalance`, `getResistingForceIncInertia`) and `LadrunoCSTPair.cpp` (same two) call `this->getMass();` purely to refill class-static `K`, then read `K(i,i)` directly — byte-for-byte the idiom that silently corrupted Quad/LST inertia when the G2 cache landed there (a hit skips the formation; `K` holds the last *tangent*). Correct today only because these two elements have no cache; the ledger row for `LadrunoMassCache` lists them as "NOT applied (trivial formation)", which invites exactly the future extension that would fire all four sites at once.
- **Workaround/status:** in-source `// Ladruno (ADR-77 review wave): DO NOT add the G2 cache here without first rewriting this to consume getMass()'s return` breadcrumbs at all four sites. If the cache is ever extended: fix the callers FIRST (the Quad/LST pattern, PR #650), and remember the Zone-A dynamics/Rayleigh battery — not A/B equality — is the oracle that catches it. *2026-07-27 (ADR-77 review wave).*

### `-noMassCache` is per-instance and NOT serialized — the A/B escape hatch is effective on the building rank only
- **Bites:** under SP/MP/DDM, broker-built remote copies of a cache-bearing element default-construct the cache **enabled** (the flag rides no sendSelf vector — deliberate T7/brick policy). A user running `-noMassCache` as a workaround or as the off-arm of a parallel A/B gets the cache silently re-enabled on every non-building rank, so a "cache-off" parallel comparison is not what it claims to be. Results remain bit-identical by the guard construction (G-BYTE), so this is a diagnostics honesty gap, not a correctness one — but it is precisely the parallel context where an escape hatch matters most.
- **Workaround/status:** documented in the `LadrunoMassCache.h` contract (review wave); serializing the flag would cost a stream-format change on six elements for a diagnostic — declined. Treat serial runs as the authoritative A/B arena. *2026-07-27 (ADR-77 review wave).*

## The splash banner breaks `subprocess.run(text=True)` test harnesses on cp1252 consoles (ADR 78)

- **Bites:** tests that spawn a python child with `capture_output=True, text=True`
  and import the engine (`test_ladruno_up_element_th.py` winding gate,
  `test_ladruno_up_mp_smoke.py` serial roundtrip) die with
  `TypeError: argument of type 'NoneType' is not iterable` on a cp1252-locale
  Windows box: the child prints the splash banner, whose UTF-8 art/feature text
  contains bytes cp1252 cannot map (e.g. the superscript minus in `F0⁻¹`,
  byte 0x81, present since the StagedDefGrad banner line), the reader thread
  raises `UnicodeDecodeError`, and `proc.stdout` comes back `None`. Looks like a
  physics regression; is an encoding trap. Predates ADR 78 (verified: the byte
  is in the committed `banner_features.txt`).
- **Workaround/status:** set `LADRUNO_OPENSEES_QUIET=1` in the child env (both
  tests pass then), or pass `encoding="utf-8", errors="replace"` to
  `subprocess.run`. Proper fix (harden the harnesses) spun off as its own task.
  Keep NEW banner-feature lines ASCII-only. *2026-07-28 (ADR 78).*

## u-p storage coupling under corot: velocity contraction is chord-poisoned; incremental-coupling-only is a pump (ADR 78)

- **Bites (two ways, both measured):** (1) contracting the damp p-row `(R̄Q)ᵀ`
  against integrator velocities of a ROTATING body picks up the chord's
  apparent volumetric rate `2(1−cosΔθ)/Δt` per step — one-signed (a systematic
  dilation in the current frame; never averages out),
  amplified by `Q̄ ≈ K_f/n` undrained (order-1 spurious p at bearing-mechanism
  rotations); no velocity-linear operator can remove it. (2) Fixing ONLY the
  coupling to the incremental `QᵀΔu_d/Δt` while leaving `S·ṗ` on the Newmark
  velocity breaks the skew-symmetry of the discrete coupling pair and PUMPS the
  structural ringing (consolidation column grew a ±100·q p-oscillation at
  Δt=0.02 where linear decays).
- **Workaround/status:** the WHOLE p-row rate block goes incremental under
  corot — `QᵀΔu_d/Δt + (S+αH̃)Δp/Δt` (the Book's GN22/GN11 pairing);
  first-order convergent (measured 2.3e-2→4.2e-3 over Δt 0.16→0.02) and
  rigid-motion exact. Full analysis: ADR 78 §3.3. *2026-07-28 (ADR 78).*
### `Element::setResponse` opens its OWN `ElementOutput` tag — chain to it AFTER `output.endTag()`, never before
- **Bites:** `Element::setResponse` is not a stub. It serves `force`/`globalForce`/`dampingForce`/`dynamicForce`/`inertialForce` and is the only way an element gets those for free — but it *also* calls `output.tag("ElementOutput")` + `attr(...)` itself. Delegate to it from inside your own open tag and the XML/MPCO stream gets a **nested duplicate** `ElementOutput`, which silently corrupts the column metadata rather than failing.
- **Rule:** `output.endTag(); if (theResponse == 0) return this->Element::setResponse(argv, argc, output);` — end your tag first, then delegate. `LadrunoDispBeamColumn2d` had this right from the start; every other fork element (its own 3d twin included) simply did not chain at all, so `inertialForce` / `dampingForce` / `dynamicForce` were unavailable fork-wide until the token sweep. Mirror it in `getResponse` with `return this->Element::getResponse(responseID, eleInfo);` instead of `return -1` — the base IDs are 111111/222222/333333/444444, so they cannot collide with an element's own.
- *2026-07-28 (recorder-token consistency sweep).*

### An element that overrides `setRayleighDampingFactors` to ignore Rayleigh damping makes 11 `Element` methods dereference `theMatrices[-1]`
- **Bites:** the `if (index == -1) this->setRayleighDampingFactors(...)` self-heal at the top of `getRayleighDampingForces`, `getResistingForceIncInertia`, the five `*Sensitivity` getters and `getGeometricTangentStiff` calls the **virtual**. An override that just stores factors and returns (the honest thing for a pure-penalty tie, and what six fork elements do — `LadrunoUP`, `LadrunoRigidBody`, `LadrunoEmbeddedNode`, `LadrunoEmbeddedRebar`, `LadrunoKinematicCoupling`, `LadrunoDistributingCoupling`) never reaches the base allocator, so `index` stays `-1` and the very next line indexes the static scratch pool at −1. **Access violation, not a wrong number** — the whole interpreter dies.
- **Not a fork problem.** Upstream's own `Subdomain` (which IS an `Element` subclass) overrides the virtual as a pure forwarder to `Domain::setRayleighDampingFactors` and likewise never sets `Element::index` — vanilla OpenSees carries the identical crash.
- **Why it stayed hidden:** nothing *reached* those paths on those elements. The moment the token sweep made them chain to `Element::setResponse`, `eleResponse(tag, "dampingForce")` routed straight into it and killed pytest on the first run.
- **Also note `Element::getRayleighDampingForces` is NOT virtual.** `LadrunoUP` and `LadrunoRigidBody` declare their own — that is *shadowing*, not overriding, so anything calling it from inside `Element` (e.g. `Element::getResponse(222222)`) gets the BASE version. Both elements therefore answer `dampingForce` themselves rather than letting the chain do it; for `LadrunoUP` the base version would additionally compose `betaK*getTangentStiff()` and inject the −Q/H coupling blocks into what it calls "Rayleigh damping", which its class contract explicitly forbids.
- **Workaround/status:** ✅ FIXED in `SRC/element/Element.cpp` — qualified to `this->Element::setRayleighDampingFactors(...)` at all 11 sites (see [[LEDGER_vanilla_files]], also filed as upstreamable). **If you add an element that overrides `setRayleighDampingFactors`, this is now safe — but if you shadow `getRayleighDampingForces`, answer `dampingForce` in your own `setResponse` or the recorder will report the base's composition, not yours.** *2026-07-28 (recorder-token consistency sweep).*

### Recorder-token aliases must never change the emitted `ResponseType` labels
- **Bites:** the obvious way to add an alias is to canonicalize the token and let one branch serve several spellings — which is right — but if the alias is allowed to influence what the branch *emits* (the `output.tag("ResponseType", ...)` labels, the response ID, or the vector width), then a recorder's column headers start depending on how the user spelled the token in their script. MPCO/STKO readers key on those labels; `stress` and `stresses` producing different metadata for the same numbers is worse than `stress` not working at all.
- **Rule:** `LadrunoResp::is()` (`SRC/element/LadrunoResponseTokens.h`) decides only WHICH branch runs. Everything the branch emits stays keyed to the branch. Matching is symmetric, so an element cannot accidentally test against a non-canonical spelling, and an unregistered token falls back to an exact `strcmp` — never worse than upstream.
- **Corollary for new elements:** canonicalization is GLOBAL across the table. Before giving a new element a short token (`k`, `dir`, `L`, `gap`, …), check the table — those are already claimed by the penalty/tie and rigid-body families. *2026-07-28 (recorder-token consistency sweep).*

### A token that answers to the wrong quantity can permanently SHADOW a later branch in the same `strcmp` ladder
- **Bites:** `LadrunoEmbeddedNode::setResponse` listed `localForce` twice — once as a second spelling of the global-component tie traction (branch 1) and once, 16 branches later, for the genuine D9 local/interface-frame force (branch 17). First match wins, so branch 17's `localForce` was **dead code from the day it was written**, and anyone recording `localForce` silently got global components. `LadrunoIMKBeam`/`2d` had the same lie without the shadow: `force`, `globalForce` and `localForce` all returned `getResistingForce()`.
- **Tell:** a spelling appearing in two branches of one ladder is always a bug — either a copy/paste or a name that means two things. `grep -c '"localForce"'` per file is a cheap audit.
- **Workaround/status:** ✅ FIXED — `localForce` now resolves to the local-frame branch on `LadrunoEmbeddedNode`, is withdrawn from `LadrunoEmbeddedRebar` (no local frame exists there), and is a REAL element-frame force on both IMK beams (basic-`q` mapping, as DispBeamColumn). **Behaviour change for anyone who recorded `localForce` on those elements** — noted in [[LEDGER_implementations]]. *2026-07-28 (recorder-token consistency sweep).*

### A singular tangent can pass on Windows and fail on Linux — a green local run is not proof the model is well-posed
- **Bites:** `tests/test_ladruno_response_tokens.py::test_distributing_coupling_tokens` tied a `LadrunoDistributingCoupling` reference node to two COLLINEAR hex corners, leaving the reference rotation about that line unconstrained. Same source, same `system BandGeneral` + `numberer RCM`: the **Windows/MKL build returned `analyze == 0`**, Zone-A on **Ubuntu failed the step** with `BandGenLinLapackSolver::solve() -factorization failed, matrix singular U(i,i) = 0`. The LAPACK behind `BandGeneral` differs between the two (MKL vs the CI's reference/OpenBLAS) and so does the pivot at which it calls a matrix singular.
- **Consequence:** a model that is *actually* degenerate can develop and pass locally for as long as you only run it on Windows, then fail the moment CI builds it. The green local run proved nothing; the Linux failure was the correct answer. Same family as the ADR-76 gate ("a singular matrix must not report SUCCESS") — that one is about the reported status of a known-singular solve, this one is about the two platforms disagreeing on whether the solve is singular at all.
- **Rule for new test models:** when an element PRINTS a well-posedness warning, treat it as an error in a test. This element said `reference rotation about axis (1, 0, 0) is unconstrained (degenerate independent set)` on the very run that "passed" — the diagnostic was right and the assertion was wrong. For `LadrunoDistributingCoupling` / `LadrunoKinematicCoupling` specifically, the independent node set must SPAN (non-collinear for a rotation-carrying reference node); check `nKept` == the number of rotation axes you expect.
- *2026-07-28 (recorder-token consistency sweep).*

### `timeSeries Path -time/-values -useLast` silently dropped the flag — Path-driven sp snapped to zero at the FINAL analysis step
- **Bites:** any Python/interpreter model driving a `Plain` pattern's `sp` with `timeSeries("Path", tag, "-time", ..., "-values", ..., "-useLast")`. The interpreter parse (`OPS_PathSeries`, `SRC/domain/pattern/PathSeries.cpp`) READ `-useLast` but then constructed the `-time` variant as `PathTimeSeries(tag, path, time, factor)` — dropping the flag (5th ctor arg defaults false). Because the DOMAIN time accumulates dt in floating point, the final step's pseudo-time overshoots the last path point by a few ulps, `getFactor` takes the beyond-the-end branch, and with `useLast == false` returns **factor 0**: the constrained DOFs snap back to zero exactly at the last step.
- **Tell:** a monotone Path-driven quantity that grows correctly for N−1 steps and collapses at step N (the ADR-79 P2 undrained gate saw p: 6.16e5 → 5.2e4 at the final step, top displacement 0.0975 → 0.0000). Invisible when the final state happens to satisfy the assertion anyway — the ADR-78 corot gate-3 test asserts p ≈ 0 under rigid rotation, so its final-step snap-back to u = 0 ALSO read p ≈ 0 and passed.
- **Workaround/status:** ✅ FIXED (ADR-79 P2 PR) — the parsed `useLast` is forwarded to the `PathTimeSeries` ctor (`// Ladruno` mark; [[LEDGER_vanilla_files]] row). Belt-and-braces for test authors: give the Path an extra terminal point beyond the last analysis time so no route depends on the beyond-the-end branch. *2026-07-28 (ADR-79 P2).*

### An INSTALLED Ladruno hijacks `import opensees` in every venv it has wired — `sys.path.insert` cannot win
- **Bites:** any script that bootstraps a *worktree* build with the standard
  `sys.path.insert(0, "<worktree>/dist/bin"); import opensees`. The Ladruno
  installer writes `ladruno_opensees.pth` into the venv's `site-packages`,
  which imports `_ladruno_opensees_boot` — and that module prepends
  `C:\Program Files\Ladruno\OpenSees\bin` to `sys.path` **and runs
  `import opensees`** (to alias `openseespy` onto it) at INTERPRETER STARTUP.
  By the time line 1 of your script executes, `sys.modules["opensees"]` is
  already the installed build; the `sys.path.insert` is a no-op.
- **Tell:** the feature you just built is "missing". Measured this session:
  `ERROR LadrunoUP 1 -- -geom 'corot' not supported: only 'linear' is accepted
  (ADR 71 §2.4)` on a worktree whose own `dist/bin` build accepts it fine. The
  banner is the giveaway — it prints the *installed* build's commit, not the
  worktree's.
- **Workaround:** run campaign/testbed scripts with an interpreter that has no
  Ladruno `.pth` (the base `C:\Users\<u>\AppData\Local\Programs\Python\Python312\
  python.exe`), and **assert which engine loaded** — compare
  `os.path.dirname(opensees.__file__)` against the intended `dist/bin` and
  `raise SystemExit` on mismatch. `Ladruno_files/testbed/hypo_bearing/
  bearing_backbone.py` carries that guard as the reference pattern. Note the
  venvs that DO have the `.pth` (e.g. `opensees_env`) are still the right ones
  for apeGmsh mesh work — just not for running a worktree build.
  *2026-07-28 (ADR-79 bearing campaign).*

### A converged push increment is a property of the MESH, not of the problem — recalibrate it after any refinement
- **Bites:** reusing a step size proven on one discretization. ADR-79 P3 proved
  2.5 mm push increments on 1 m hexes; on the graded bearing mesh (0.5 m hexes
  at the surface, 4.5 B clearance, 2816 UP-H8) the same 2.5 mm increment does
  not converge **for any** geometry method — `-geom linear` included, which is
  what proves it is the model/mesh and not the geometry lane.
- **What it is NOT** (all measured, so don't re-derive them): not a tolerance
  artifact — Newton diverges (‖Δu‖ ≈ 0.2 m) and KrylovNewton stalls at
  ‖Δu‖ ≈ 3e-5 for *every* test tried, `NormDispIncr` 1e-8…1e-4,
  `RelativeNormDispIncr`, `EnergyIncr`; not a `dt` effect — 2.5 / 25 / 250 s
  all fail identically at the same displacement increment; not the u-p
  stabilization or formulation — `-stab auto 0.10/0.25/0.50`, `-stab off` and
  `-formulation bbar` all fail identically at 0.25 mm.
- **What it IS:** a FIRST-LOADING shock off the `updateMaterialStage 1` PDMY
  state. From the gravity state, 0.05 mm converges 12/12 while 0.10 mm fails —
  but once a few small increments have been taken, the SAME model happily
  accepts 0.4 mm (8x). The threshold is transient, so a fixed ladder is the
  wrong tool.
- **Rule:** drive displacement-controlled pushes with a 2-point linear ramp
  (`timeSeries Path -time 0 T -values u0 u1`) so the increment is carried by
  `dt`, then ADAPT it: halve on failure, grow back after a run of successes,
  and truncate honestly at a floor. This replaces P3's fixed `dt/10` fallback
  and self-tunes across the backbone. Also note KrylovNewton is the *primary*
  algorithm on this problem class, not a fallback — plain Newton diverged at
  every increment tested. *2026-07-28 (ADR-79 bearing campaign).*
### A relative-to-first-residual convergence test with a `count >= 1` guard can NEVER pass when the loop starts converged — and eas starts converged on every assembly
- **Bites:** `LadrunoBrick::formEAStrue` / `LadrunoQuad::formEAStrue` (the Simo-Rifai inner Newton on the enhanced parameters). The element's `update()` solves alpha, but `formResidAndTangent` calls `formEAStrue` AGAIN on the same trial state — so the inner loop routinely *enters* with `r0` at the fp noise floor (e.g. 8.4e-15). The old test `count >= 1 && r <= tolRel*r0 + tolAbs` (a) refused to declare convergence before taking one Newton step, and (b) that forced step solves a pure-roundoff RHS and lands at the material-evaluation noise floor — which for large-unit models sits ABOVE the fixed `tolAbs = 1e-12` (units: force×length², so the floor's reachability scales with E·V; the observed run stalled at 1.4e-12 vs a 1.0e-12 threshold). Result: 12 wasted material sweeps + a "did not converge" warning **per element per global iteration** — a 6-step probe wrote 37 MB of log; the queued 2000-step run would have written ~12 GB. It also inflated the eas per-step cost (~12x the inner-loop material work).
- **Tell:** warnings quoting `||r||` a few ×1e-12 against `r0` at 1e-14/1e-15 — the "residual" is smaller than anything physical; the test is comparing noise to noise.
- **Rule:** an iterative-solve floor on a residual with physical units must be SCALED. The fix (this PR) tests from `count == 0`, adds a noise floor `tolNoise * rScale` with `rScale = Σ_GP ||M^T σ dV||` (the magnitude of the terms whose cancellation produces the residual — below ~1e-12 of that, the residual is roundoff, not signal), and breaks silently on stagnation (`r >= rPrev && r <= r0`: a Newton step that fails to reduce the residual while no worse than entry means roundoff dominates). Same family as the buildEAStrue scale-invariant degeneracy check: absolute thresholds on dimensional quantities are always wrong at some unit system.
- *2026-07-29 (eas inner-Newton warning spam).*
### A worktree's `dist/bin` can be MONTHS behind its branch, and an "engine guard" that checks the module's PATH will happily wave it through
- **Bites:** any harness that pins itself to a worktree build. The `hypo_bearing`
  runner already carries a guard (scoping finding 5) against the installed
  Ladruno's site `.pth` pre-importing `opensees` — but that guard asserts
  `os.path.dirname(ops.__file__) == dist/bin`, i.e. *where* the module loaded
  from. It says nothing about *what is in it*. A fresh worktree checked out at
  a branch with ADR-78/79 merged had a `dist/bin/opensees.pyd` from an earlier
  build, which passed the location guard and then refused the feature under
  test: `-geom 'corot' not supported: only 'linear' is accepted (the axis is
  reserved ... ADR 71 §2.4)` — an ADR-71-era binary answering for an
  ADR-79-era branch.
- **Tell:** a capability error naming an ADR *older* than your branch, on a
  worktree you never built in. `git log` looks right; the `.pyd` mtime predates
  the feature commits.
- **Rule:** a location guard is necessary but not sufficient — assert
  CAPABILITY. Cheapest form is to construct the thing under test (one element
  with the flags the campaign needs) before committing hours to a run, which is
  also what catches a stale build in a *shared* checkout that another agent
  rebuilt on a different branch. Rebuild the worktree (`build.bat OpenSeesPy`)
  and verify `dist/bin/opensees.pyd` mtime moved. Cross-ref
  [[ladruno-build-in-worktree-not-shared-checkout]]. *2026-07-30 (ADR-79 locking leg).*
### Two long parallel runs writing incremental CSVs: re-running one leg's name TRUNCATES the live file under it, and the survivor keeps writing at its old offset
- **Bites:** any campaign whose legs stream results to `f = open(path, "w")` +
  `flush()` per step and run for hours in parallel — the `hypo_bearing` legs, and
  the same idiom in the perf testbeds. Starting a short smoke of a leg whose
  full run is ALREADY in flight reopens the same path with `"w"`. The smoke
  truncates the file to zero; the live process still holds its own descriptor at
  (say) byte 12000, so its next write lands there and the OS zero-fills the gap.
  Measured result: a 15 753-byte CSV that was 14 589 NUL bytes, holding the
  smoke's 11 rows, then padding, then the real run's tail — the live leg's
  entire early backbone (through s/B = 3.37%, ~80 min of compute, including the
  1% and 2.5% checkpoints) simply gone.
- **Tell:** NUL bytes in a CSV; a first data row whose settlement is *larger*
  than rows further down; a file far bigger than its line count justifies. It
  does NOT crash, and the live process's end-of-run summary still prints correct
  numbers from its in-memory rows — so the loss is silent unless you read the file.
- **Rule:** two guards, both cheap. (1) A capped/smoke run must write a
  DIFFERENT filename (`backbone_<leg>__smoke.csv`), so a smoke can never address
  a real leg's output. (2) Refuse to open an existing output file modified within
  ~180 s — that means another process is actively appending — with an explicit
  env override for the case where you know it is dead. Note Windows does not
  block the second open, so nothing but your own guard prevents this.
  *2026-07-30 (ADR-79 locking leg).*
### A frictional verification model can be 95 % mobilised by its own gravity state, because the elastic K0 is what sets the initial stress ratio
- **Bites:** any limit-analysis / bearing-capacity check that picks a "nice
  easy friction angle" to validate the machinery against a closed-form answer.
  Under 1-D gravity or surcharge loading the initial stress ratio is the
  ELASTIC `K0 = nu/(1-nu)`, not `1 - sin(phi)`, so the mobilisation of the
  yield surface before anything is loaded is
  `m = (1-K0) / (sqrt(3) * alpha * (1+2*K0))` and depends only on nu and the
  surface — never on the load magnitude. At PDMY's moduli (`K = 1.5e5`,
  `G = 5.5e4` => nu = 0.3366, K0 = 0.507) that mobilises 19.1 deg. Validating
  against a phi_txc = 20 deg cone therefore starts with the WHOLE DOMAIN at
  **m = 0.950 of yield**, and the "verification" measures its own initial
  condition: 44.5 % of elements at m > 0.99 after 2 mm of footing settlement,
  the yielded zone already touching the roller boundary, and no convergence
  past s/B = 0.0019.
- **Tell:** a validation leg that dies almost immediately with a large fraction
  of the mesh yielded and the plastic zone at the boundary, while the SAME
  machinery runs fine on the (steeper, supposedly harder) real surface.
- **Rule:** compute `m` of the initial state and print it before the run; treat
  `m > 0.8` as void. The fix is free — a collapse load of an
  elastic-perfectly-plastic body does not depend on its elastic constants, so
  raise nu for the verification legs only (nu = 0.45 gives K0 = 0.818,
  m = 0.268) and MEASURE the independence with a nu-pair rather than asserting
  it. *2026-07-30 (ADR-79 collapse study).*
### `nDMaterial DruckerPrager` cannot have pressure-dependent moduli AND plasticity, and `updateMaterialStage` does not reach it
- **Bites:** anyone reaching for `DruckerPrager` as a cheap perfectly plastic
  soil next to PDMY. `mElastFlag` is `0 = elastic+no param update, 1 =
  elastic+param update, 2 = elastoplastic` (default), and `updateElasticParam`
  only rescales `K, G` by `sqrt(1 + p/p_atm)` when `mElastFlag == 1` — which is
  an ELASTIC state. So the PDMY-style `G ~ sqrt(p)` and plasticity are mutually
  exclusive; a plastic DruckerPrager has constant moduli, full stop. Separately,
  `setParameter` returns -1 for `"updateMaterialStage"` (it answers to
  `"materialState"` instead), so the usual `ops.updateMaterialStage(...)` staged
  gravity idiom silently does nothing here.
- **Rule:** for a collapse load this does not matter (a limit load is
  independent of the elastic constants), but the SETTLEMENT at which it arrives
  is not — so any `q` quoted at a fixed `s/B` criterion is affected, and that
  must be stated. For the gravity stage, check whether the elastic K0 state is
  admissible (previous entry) and just solve it plastically in one step rather
  than hunting for a stage flip that is not wired up. Perfect plasticity needs
  `Kinf = Ko = delta1 = delta2 = H = theta = 0`; `mHprime = (1-theta)*H` so
  `H = 0` alone suffices. The conversion from a measured cone is
  `rho = sqrt(2)*alpha` and the sqrt(J2) cohesion intercept is `sigma_y/sqrt(3)`
  (yield is `||s|| + rho*I1 - sqrt(2/3)*sigma_y`); the tension cutoff `mTo` is
  placed exactly at the cone apex, so it adds nothing.
  *2026-07-30 (ADR-79 collapse study).*
### Consolidating with a `Linear` time series and then adding a deviatoric pattern applies the deviator AT FULL AMPLITUDE on the first increment
- **Bites:** every single-element stress-path probe built as "ramp the
  confinement to lambda = 1, then add a second pattern and keep stepping". Both
  patterns scale with the SAME load factor, so the first post-consolidation
  increment (lambda = 1.0025) evaluates the deviatoric pattern at 1.0025x its
  full amplitude, not at 0.0025x. The probe then reports failure at step 0 with
  ZERO measured deviator — `cone_probe.py` avoided this with `Path` series keyed
  to pseudo-time and it is easy to miss when porting the same probe to a
  different material.
- **Tell:** every path of a stress probe returns the consolidation state
  unchanged (`sqrt(J2) = 0`, `alpha = 0`) and "failure at load factor ~1.00".
- **Rule:** `ops.loadConst("-time", 0.0)` after the consolidation stage, so the
  held pattern stops scaling and the new pattern ramps from zero. Same rule as
  the push stage of any displacement-controlled runner.
  *2026-07-30 (ADR-79 collapse study).*
### "Is the plastic zone clear of the boundary?" cannot be tested as a fraction of the domain extent on a GRADED mesh
- **Bites:** any field post-processing that asks whether a mechanism is
  contained. The obvious test — is the outermost yielded element's centroid
  within 90 % of the half-width? — is silently wrong when the mesh grades
  outward, because the far-field element is enormous. On the `hypo_bearing`
  benchmark the x-columns run `..., 4.00, 5.83, 8.45` in a 10 m half-domain:
  the LAST element is 3.1 m wide and its centroid sits at 8.45 m, so an element
  physically TOUCHING the roller passes a "< 9 m" test. Measured: the collapse
  leg reported "contained" while 16 of the 352 elements in the
  boundary-adjacent column were at m > 0.99.
- **Tell:** the reported extent equals one of the coarse outer centroids
  exactly, and equals it for several different legs — because it is a mesh
  coordinate, not a mechanism coordinate.
- **Rule:** test **element-column / row membership**, not distance. Build
  `np.unique(np.round(np.abs(centroid[:,0]), 4))` and compare against
  `cols.max()`; report `n_yielded / n_in_column` so the reader sees how much of
  the boundary layer is engaged rather than a yes/no. Same trap applies to any
  "did the wave reach the absorbing boundary" or "is the damage localized"
  check on a graded or adaptively refined mesh.
  *2026-07-30 (ADR-79 collapse study).*
### `Vector::operator()` is UNCHECKED outside `_G3DEBUG` — a mis-sized `static Vector` in `getResponse` is a silent heap overrun, not a caught index error
- **Bites:** any `getResponse` that fills a fixed-size scratch `Vector` in a
  Gauss-point loop. Found in vanilla `TenNodeTetrahedron::getResponse`, where
  `static Vector stresses(6)` — copy-pasted from `FourNodeTetrahedron`, which
  has ONE Gauss point so 6 is right there — is written by BOTH the `stresses`
  and `strains` branches looping over all `NumGaussPoints=4` points × 6
  components = **24 doubles into a 6-double block**. 18 doubles / 144 bytes past
  the end, on every recorder step. `Brick` gets it right (`stresses(48)` for
  8 GP), so the pattern is fine — only the size was left behind when the loop
  bound was edited.
- **Why it is invisible:** `Vector::operator()` bounds-checks ONLY under
  `_G3DEBUG` (`SRC/matrix/Vector.h`), which release builds do not define — the
  checked accessor is `operator[]`, which nobody uses in element code. And the
  buffer is `static`, so it is heap-allocated once on first call and the SAME
  144 bytes past it are stomped forever after. The crash therefore surfaces at
  whatever the allocator placed next — often a later `free`/alloc far from the
  write — so the backtrace does not point at the element.
- **Second-order damage even without a crash:** `Information::setVector` does
  `*theVector = newVector`, and `Vector::operator=` REALLOCATES on size
  mismatch. So the `ElementResponse`'s advertised `Vector(6*nGP)` gets silently
  shrunk to the scratch size, while `ElementRecorder` already sized its columns
  from the advertised size at setup (`ElementRecorder.cpp`) and then copies
  `eleData.Size()` per element — the column layout desynchronises across
  elements. Garbage output, no diagnostic.
- **Measured A/B** (one tet10, uniform-strain patch, `eleResponse(ele,'stresses')`):
  pre-fix build returns **6** values, post-fix returns **24** — so the cheap,
  deterministic tell is a response list that is a WHOLE FRACTION of the
  advertised length, one GP block instead of nGP. Both builds report
  `σxx = 10000.0 = E·ε` exactly, i.e. the physics was never wrong — only the
  buffer. Note the 1-element control did NOT crash: the overrun is real on
  every call but whether it segfaults depends on what the allocator put after
  the block, so a small repro proving "no crash" proves nothing. Trust the
  length, not the absence of a crash.
- **Rule:** size the scratch from the same expression `setResponse` advertises
  (`6*NumGaussPoints`, never a literal), and treat any `static Vector`/`Matrix`
  in a response path as a place to check the loop bound against the declared
  size. When a new element is derived by copy-paste from one with a different
  Gauss-point count, audit every fixed size in the file, not just the loops.
  Cross-ref [[ladruno-adr79-geom-hypo]] (tet10 was the recorder in use).
  *2026-08-04 (tet10 recorder segfault).*
### A numpy ORACLE can silently lag a reviewed C++ fix — a red zone_a gate may indict the reference, not the shipped code
- **Bites:** any subsystem gated by a hand-written numpy oracle that mirrors a
  C++ kernel (`concrete3d_ref.py`↔`LadrunoConcrete3DKernel.h`, and the same
  pattern in the J2/logstrain/up/hypo `*_reference.py` pairs). A fix applied to
  the kernel during PR review does NOT propagate to the oracle, and nothing in
  CI notices: the kernel-vs-oracle test compares them on a COMMITTED FIXTURE of
  paths that may never exercise the diverging branch.
- **The case:** `test_p2i_multiaxial_apportioning_gate` asserts
  `I3_pure_compression_wt < 1e-9` (no spurious compression→tension damage) and
  measured **0.997**. It looks like a formulation bug in the tensile damage
  gate. It is not — the gate `sig_t_drive = E*et if max(w) > 1e-6*ft else 0`
  is CORRECT and opens legitimately, because the effective stress handed to it
  genuinely IS tensile: under uniaxial-STRAIN compression the hardening Newton
  overshoots to `rho<0`, the apex branch teleports to the hydrostatic-TENSION
  vertex, and a trial with max principal **−23.76** returns **[+2.94,+2.94,+2.94]**
  with `conv=True`. `f==0` holds at the apex BY CONSTRUCTION, so the oracle's
  `(converged or apex) and |f_indep|<tol` reports success for a sign-flipped,
  inadmissible state. `kp` also jumps 0.034→0.182 in one step.
- **The kernel was already right.** PR #249's adversarial review added to
  `returnMapHardening` an admissibility test (`dlam>=0 && kp>=kp_n`), refusal to
  report converged, and a fallback to the ELASTIC PREDICTOR so the caller cuts
  the step — with the comment "this deliberately diverges from the numpy
  oracle's (equally-arbitrary) apex teleport — the kernel is the safe reference
  here." The oracle received only the HONEST-f-recompute half of that fix and
  kept the literal pre-fix expression the C++ comment calls out as lying. So
  **the shipped material never had the bug**; only the reference did.
- **It shipped red and stayed red.** The gate produces the byte-identical
  0.9971183764898133 at `c349e8763` (#336), the very commit that introduced it —
  #336 landed AFTER #249, so the gate was authored against an oracle that
  already lagged the kernel. It has never passed.
- **Rule:** when a zone_a gate goes red on a subsystem that has BOTH a numpy
  oracle and a C++ kernel, diff the two implementations of the disputed branch
  BEFORE touching the assertion — the oracle is as likely to be stale as the
  code. And when a PR-review fix lands in a kernel, port it to the oracle in the
  SAME PR, because the fixture-based agreement test will not catch the drift.
  Do not "fix" a gate by weakening its assertion until you have established
  which side is wrong. Cross-ref [[ladruno-adr79-geom-hypo]].
  *2026-08-04 (P2i apex-teleport hunt).*
### quad8 mortar dual masses are NEGATIVE at corners (−A/12) by construction — sign-aware guards, do not "fix"
- **Bites:** the serendipity quad8 corner shape functions integrate to **−A/12** over the facet (mids +A/3; the total is still A), so the ADR-62 P2.1 dual condensation `Aᵉ=diag(∫N)(Dᵉ)⁻¹` legitimately produces NEGATIVE diagonal dual masses at every quad8 corner node. The sign cancels exactly in `P = Mdual/Ddual` (partition of unity `P·1=1` holds algebraically), but any guard, recorder, or future reader that assumes `Ddual > 0` — the shipped `Ddual[I] <= 1e-300` "uncovered node" refusal did exactly this — false-refuses every valid quad8 tie. Same trap wherever a rowsum `∫N_I ≥ 0` assumption hides: the pre-ADR-78 coverage ratio (`cover/fullCov ≥ 1−1e-3` FLIPS for negative rowsums) and the `|gap|/cover` normalization.
- **Rule:** for serendipity bases, node-wise rowsum measures are SIGNED; guards must be sign-free (areas, L1 integrals) or sign-aware (`|Ddual|`). ADR-78 D3 unified the mortar-tie guards on per-facet AREA coverage + per-facet `∫|g_N|/area`.
- **Workaround/status:** ✅ shipped that way (ADR-78). *2026-08-04.*

### tri6 SLAVE facets are structurally incompatible with the dual mortar — corner ∫N = 0, refused by name
- **Bites:** on the reference triangle the tri6 CORNER shape functions integrate to exactly ZERO (the three midsides carry the whole area). The dual scaling divides by the per-facet corner rowsum ⇒ division by zero, and no tolerance rescues it — it is structural, not conditioning. The rowsum-based coverage machinery is equally meaningless there. quad8 corners (−A/12 ≠ 0) are fine.
- **Rule:** `LadrunoTie -mortar` refuses `npsS == 6` in BOTH bases (ADR-78 D2, mirroring apeGmsh ADR 0086 v1); tri6 MASTER facets are fully supported. Remedy in the message: swap master/slave, or put quad8/hex20 faces on the slave side. Revisit only if a tet10-interface user materializes.
- **Workaround/status:** ✅ named refusal shipped (ADR-78). *2026-08-04.*

### ADR-78 unified the mortar-tie guards for ALL facet orders — linear decks with cancelling gap fields now refuse (by design)
- **Bites:** the pre-ADR-78 conforming-gap guard tested the per-node SIGNED weighted gap `|Σ∫N_I g_N|/cover`, which a gap field that cancels inside a node's support could slip through (a warp/antisymmetric offset reading as "conforming"). The ADR-78 per-facet `∫|g_N|/area` L1 guard has no cancellation blind spot and — per the OQ-2 "unify" sign-off — applies to tri3/quad4 decks too. A linear deck that previously built its tie may now refuse with the conforming-gap message. The emitted P (weights) is byte-identical for linear inputs; only refusal behaviour changed.
- **Rule:** if a formerly-working linear tie now refuses on the gap guard, the geometry genuinely is off the master surface somewhere — fix the interface or consciously relax `-tol`.
- **Workaround/status:** ✅ intentional behaviour change, documented here + ADR-78 D3/BLOCKER-3; regression-tested (`test_refuse_accordion_gap_L1`). Two adversarial-gate amendments same day: (a) the threshold scale is `0.5·sqrt(areaFull)` — a bare `sqrt(A)` was 2× LOOSER than the shipped per-node tributary scale (review MINOR); (b) the per-facet area-coverage sum counts MULTIPLICITY, so a self-overlapping / doubly-listed MASTER surface could exactly mask an uncovered slave strip — now a named refusal (coincident master-master overlap, mean-gap-gated so curved masters stay legal; review MAJOR, `test_refuse_duplicate_master_facets`). A linear deck with duplicated master facets that previously "worked" now refuses. *2026-08-04.*

### A non-homogeneous `sp` makes `enforceSPs` PRE-UPDATE the driven element layer, so the first constitutive evaluation of every increment is over-strained by L/h — harmless for an elastic law, ×28 iterations for a plastic one
- **Bites:** a static push driven by a prescribed displacement on a face (`sp` in a `Plain` pattern + `Linear` series, `constraints Transformation`, `integrator LoadControl`) costs an order of magnitude more solver work **if the elements at that face carry a path-dependent material**. Measured on a 3D solid, small-strain `LadrunoJ2`, `H/E` = 1 %, all four runs pinned to `KrylovNewton`: driven-face elements **elastic** → 30 increments / 0 cutbacks / **149** Newton iterations; **plastic** → 324 / **52** / **4 255** (**×28.6**). Same converged answer to **0.00 %** (920.0 vs 920.1 kN at a stated internal rotation) — pure conditioning. Replace the prescribed displacement with a **traction** on the same face and the plastic material is FREE (30 / 0 / 145, ×0.95) — *even though the cover then grossly plastically hinges*, while the expensive displacement-driven run's covers **never yield at all** (max σ_vM 295 MPa vs `f_y` 379.5, 0.00 % plastic at every step). So "the BC singularity yields and wrecks the tangent" is **refuted**: a fully hinging cover converges in 4.6 iterations/increment, a non-yielding one takes 13.4 with 35 cutbacks.
- **Why:** `LoadControl::newStep()` calls only `applyLoadDomain(λ)` (`SRC/analysis/integrator/LoadControl.cpp:130`) — **statics has no predictor**. That reaches `TransformationConstraintHandler::applyLoad()` → `enforceSPs()` (`SRC/analysis/handler/TransformationConstraintHandler.cpp:496-524`), which writes the full prescribed value into each constrained node (`TransformationDOF_Group.cpp:1071`, `setTrialDisp`) **and then calls `theEle->updateElement()` on every element touching a constrained node, at lines 518-521**. So the first constitutive evaluation happens *inside `newStep`*, with the driven face advanced by the whole increment and every interior node still at the last converged position: that layer sees **Δδ/h_element** instead of the physical **Δδ/L_model**. The return map yields it spuriously, its consistent tangent collapses toward `2G·H/(3G+H)` ≈ 1 % of elastic, the predictor built on it overshoots, the layer unloads elastically next iteration, and the active set oscillates. Trial states are recomputed from the committed state each iteration, so **nothing survives into committed state** — which is exactly why the answer is right and only the iteration count betrays it. The pre-update is not gratuitous: the `sp`'d DOF is genuinely eliminated (`setID(dof,-1)`, `TransformationDOF_Group.cpp:1034`; numberers assign equations only to `-2`) and `TransformationFE::getResidual` is only `Tᵗ·R` from element state (`TransformationFE.cpp:391-394`), so there is no column for a `K·Δu_prescribed` term — pre-updating is how the prescribed motion reaches the RHS at all.
- **Diagnosis trap:** raising `sig0` out of reach "fixes" it in 2 iterations, which looks like it exonerates the constraint and indicts yielding. It does neither — the mechanism is *yielding caused by the constraint enforcement*, a conjunction, and that control removes only the yielding half. The discriminating test is the **traction** drive, which keeps the material able to yield and removes the constraint.
- **Contrast worth knowing:** `DisplacementControl::newStep` already has the right ordering — `formTangent → solve → deltaU → incrDisp(deltaU) → applyLoadDomain(λ) → updateDomain()` (`SRC/analysis/integrator/DisplacementControl.cpp` ~35-98), i.e. the interior is advanced by a tangent-consistent predictor *before* the SPs are enforced. Abaqus avoids the whole thing with a default-ON 100 % extrapolation predictor; Kratos ships `use_old_stiffness_in_first_iteration` (assemble iteration 1 at the converged state, inject `b −= K·Δu_D`) and **enables it in its own J2 / plastic-damage / fatigue tests while disabling it for elastic ones**.
- **Workaround/status (2026-08-04):** ⚠ **open — no code change.** Practical workaround: give the load-introduction volumes an **elastic** material (measured neutral to 0.00 % at equal internal rotation; it is what elastic load blocks / rigid platens achieve in commercial practice, though **no code's manual actually recommends it** — searched). Candidate fix is a static **predictor** (`LadrunoLoadControl -extrapolate`), NOT a new SP handler; do **not** "fix" it by adding `updateDomain()` to `LoadControl::newStep` (the driven layer is already updated, so that buys an extra full constitutive sweep per increment for zero behaviour change). Evidence, cross-code comparison and validation gates: [[80_sp_prescribed_displacement_findings]]. Related fork SP fix: the `sp -subtractInit` no-op row in [[LEDGER_vanilla_files]].

### `AutoConstraintHandler::applyLoad()` omits the `updateElement()` loop that Transformation has — with a non-homogeneous `sp` that is a SILENT WRONG ANSWER, not an error
- **Bites:** `constraints Auto` looks like a drop-in replacement for `Transformation`, and for homogeneous BCs it is. With a **non-homogeneous** `sp` (a prescribed, load-factor-scaled displacement) an increment can be committed with only the boundary layer moved and equilibrium never re-checked — no warning, no error, plausible-looking output.
- **Why:** `AutoConstraintHandler::applyLoad()` (`SRC/analysis/handler/AutoConstraintHandler.cpp:573-584`) runs the two `enforceSPs` loops and **stops**, where `TransformationConstraintHandler::enforceSPs` follows them with `theEle->updateElement()` over every element touching a constrained node (`TransformationConstraintHandler.cpp:518-521`). Without that loop the first residual of the increment is formed from **stale** element state, so `dU ≈ 0` comes out of the first solve — and `CTestNormDispIncr` has **no minimum-iteration guard** (`SRC/convergenceTest/CTestNormDispIncr.cpp:160-175`: `if (norm <= tol) return currentIter;`, with only a `currentIter == 0` "start() never invoked" check), so it reports convergence at iteration 1. Neither `NewtonRaphson.cpp` nor `KrylovNewton.cpp` guards against it either.
- **CONFIRMED and FIXED (ADR-80 gate G4, 2026-08-04).** No longer an inference. Measured A/B on the pre-fix binary, two equal trusses driven to δ=3.0 with the midpoint free: `Transformation` gives **u_mid = 1.5 in 2 iterations**, `Auto` gives **u_mid = 0.0 in 1 iteration** with `analyze()` returning **0**. Identical signature on two stacked `stdBrick`s. The **single Newton iteration is the diagnostic tell** — it means the first residual carried no information about the prescribed motion.
- **The inference UNDERSTATED the blast radius: `EnergyIncr` is fooled too**, not just `CTestNormDispIncr`. Energy is `dU·R`, so anything built on `dU` collapses with it. Measured: `NormDispIncr` → 0.0 (wrong), `EnergyIncr` → 0.0 (wrong), `NormUnbalance` → 1.5 (right).
- **Test-dependence IS the mechanism, and it is why this hid for so long.** `NormUnbalance` looks at the large *stale* residual, refuses to converge at iteration 1, and the very next `LoadControl::update` → `updateDomain()` silently repairs the state — so **the same deck is right or wrong purely by choice of convergence test.** If you are auditing an old `constraints Auto` result, the question to ask is not "did it converge" but "what did `test` say".
- **Contrast — `Plain` is also wrong here but LOUD:** it prints `non-homogeneos constraint for node N homogeneous constraint assumed` and gives 0.0. That is documented by-design behaviour (Plain cannot do non-homogeneous SPs), not this bug. `Auto` was wrong and **silent** — that is the whole defect.
- **Workaround/status:** ✅ **FIXED** — `theFEs` membership + the `updateElement()` loop in `applyLoad()`, mirroring the Transformation handler; see the `AutoConstraintHandler` row in [[LEDGER_vanilla_files]]. Gated by `tests/test_auto_handler_sp_update.py` (verified to FAIL on the pre-fix binary). Homogeneous-only (`fix`-only) models are numerically unaffected. Still do not reach for `constraints Auto` as a workaround for the `sp` × plasticity *conditioning* quirk above — that one is unrelated and remains open. See [[80_ladruno_sp_imposition_strengthening_adr]] S2 and [[80_sp_prescribed_displacement_findings]] §4 candidate B.

### The 16 `*_cpp.py` kernel self-checks FAIL under plain `pytest` and PASS under `pytest -s` — it's stdout capture breaking `subprocess`, not a regression
- **Bites:** a full Zone-A sweep reports `16 failed, 1742 passed`, and every failure is a C++ kernel self-check (`test_hypo_kernel_cpp`, `test_ladrunoJ2_*_cpp`, `test_logstrain*_cpp`, `test_ladrunoConcrete3D_material`, `test_ladrunoCMS_core_cpp`, `test_ladruno_up_kernel_cpp`, `test_ladrunoRCConcrete_{reg,tensstiff}_cpp`, …). It reads as if a source change just broke every material kernel at once. It did not.
- **Why:** those tests shell out to `g++` via `subprocess.run(..., capture_output=True)`. Under pytest's default capture the parent's stdio handles are not real console handles, and Windows fails the spawn with **`OSError: [WinError 50] The request is not supported`** at `subprocess.py:1416` — *before* any C++ is compiled or any OpenSees code runs. `g++` itself is fine (`/c/msys64/mingw64/bin/g++`, MSYS2 rev8 15.2.0) and the same `subprocess.run` succeeds from a plain `python3.12 -c`.
- **Tell:** the failure is an `OSError` at spawn time, never an assertion or a numeric mismatch. A real kernel regression fails on *values*.
- **Workaround/status (2026-08-04):** run those tests with **`pytest -s`** (capture disabled) — `16 passed` on the identical binary that "failed" 16. So the honest Zone-A total is **1758 passed, 0 real failures**. Verify a suspected kernel regression with `-s` **before** believing it. Two further modules (`test_ladruno_overlay_{driver,physics}.py`) cannot even be *collected* here because `matplotlib` is absent — also unrelated; `--ignore` them or install matplotlib.

### A cross-check patch built on RECTANGLES cannot detect a quadrature mismatch — patch-exactness and PoU are both quadrature-independent, so the rule hides until the facets stop being parallelograms
- **Bites:** the ADR-78 → apeGmsh mortar contract (R3) told the porter to match the fork's 12-point degree-6 Dunavant rule "so the cross-check matches to 1e-12, not 1e-6", and R7 defined that cross-check on uniform grids over [0,1]². apeGmsh's port (ADR 0086 S1, PR #898) used a completely different rule — a Duffy-collapsed 5×5 tensor rule, 25 points — and **reproduced the fork's `P` to 5.0e-13 anyway**. R3 had no teeth: a porter could ignore it entirely and still pass R7.
- **Why:** on an affine facet (rectangle, parallelogram, any triangle) the isoparametric pullback of the `N·φ` product is a POLYNOMIAL of degree ≤ 6, so *every* rule exact to degree 6 returns the same integral to round-off — measured 1.9e-15 fork-side, 4.5e-15 apeGmsh-side. Make the facet a non-parallelogram and the pullback goes RATIONAL: no finite rule is exact, and the same two rules diverge by ~1e-7 in `P` (fork: 1.4e-07 on a trapezoid, 3.5e-07 on a stronger skew; apeGmsh: 8.4e-08 / 5.9e-08). Four to five orders above the stated tolerance. The linear **patch test** cannot expose this either — it is exact under both rules on both patch types (`D` and `M` share Gauss points, so the error cancels in `P·1`); PoU likewise. The *only* signal is `P` itself, on a non-affine patch.
- **Rule:** any two-implementation agreement protocol for an integral-mortar / segment-to-segment kernel MUST include at least one non-parallelogram facet, or it is silently only testing the algebra, not the quadrature. Sanity-check the reverse too: a rule of INSUFFICIENT degree *is* caught on affine facets (the degree-4 rule shifts `P_dual` by 4.7e-04 on the same patch), so an affine-only protocol distinguishes "wrong degree" but not "different adequate rule".
- **Related:** the same protocol quoted its reference row to 9 decimals while asking for a 1e-12 row-by-row comparison (caps any row check at ~5e-10), and keyed that row to node INDICES that were described as row-major but are actually the mesh helper's `nid()` **creation** order. Publish cross-check references as a machine-readable file at round-trip float precision, keyed by COORDINATES, not indices — `kinematic_tie_validation/mortar_crosscheck_reference.json` now does.
- **Workaround/status:** ✅ fixed — contract revision 2 adds R7 **case B** (the same 2×2 quad8 / 3×3 quad4 topology bilinearly mapped onto a trapezoid; mesh lines stay straight and midsides stay at midpoints, so R1/R2 are untouched) with published reference `P`, and the oracle gates it with **T12** (case B assembly + patch exactness) and **T13** (Dunavant-12 vs Duffy-5×5: agree ≤1e-12 affine, disagree ≥1e-9 non-affine). *2026-08-04.*
