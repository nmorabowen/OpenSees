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
### SMS sizes against the betaK-DAMPED explicit step (closed form s=T^2+2Tc), NOT the undamped 2/w — and excludes alphaM
- **Bites:** with stiffness-proportional (`betaK`) Rayleigh damping the central-difference stable step SHRINKS (`xi = betaK*w/2` grows with `w`), so SMS sizing against the UNDAMPED element step `dt_e = 2/w_max` UNDER-scales — the scaled element is still unstable at `dtTarget` and the run diverges despite SMS reporting success. v1.1 sizes against the damped step: with `c = betaK/dt_e` (= 0.5*betaK*w_max), the betaK-damped step at mass-scale `s` is `dt_d(s) = (2/w0)(sqrt(s+c^2) - c)`, which inverts in CLOSED FORM to **`s = T^2 + 2*T*c`, `T = dtTarget/dt_e`** (no bisection). Reduces to the undamped `T^2` when `betaK=0` (no-damping models byte-identical). The SKIP test also moved from undamped `dt_e>=dtTarget` to damped `dt_d(1)>=dtTarget` — an element whose undamped step clears `dtTarget` but whose damped step does not is now correctly scaled.
- **alphaM (mass-proportional) is intentionally EXCLUDED from the sizing:** it does not reduce the high-frequency step that governs explicit stability (`xi_alpha = alphaM/(2w)` DECREASES with `w`), and folding it in across scales is non-monotonic (`dt_d(s)` would asymptote and may never reach `dtTarget`). Only `betaK = getRayleighDampingFactors()(1)` enters — mirroring the betaK term of `computeCriticalTimeStep`'s damped estimate, so SMS sizes consistently with what `criticalTimeStep()` reports. Gated by `tests/test_massScaling_validation.py::test_betaK_damped_sizing` (betaK injects MORE mass than undamped, matching the closed-form ratio; a betaK-blind SMS ties them). Learned 2026-06-20, mass-scaling v1.1 betaK-damped sizing, [[36_ladruno_selective_mass_scaling_adr]] / [[project_ladruno_mass_scaling]].
### LadrunoConcrete3D `-eta` (Duvaut–Lions): `dt≤0` falls back to INVISCID (β=1), NOT the elastic β→0 limit; and the oracle's deep-compression apex chaos forces stress-controlled gate states
- **dt-fallback (deliberate):** the relaxation factor is `β = dt/(η+dt)` only when `η>0 AND dt>0`; otherwise `β=1` (pure inviscid). Mathematically `dt→0` is the elastic limit (`β→0`, frozen trial), but operationally a missing/zero time increment (static, no pseudo-time) must NOT silently turn the material elastic — so `_dl_beta`/the kernel gate to inviscid. Consequence: like `-implex`, `-eta` is **inert without a positive `ops_Dt`** (transient, or uniform LoadControl/pseudo-time); it shares the IMPL-EX `dt=ops_Dt` caveat (`DisplacementControl` λ-increment is non-uniform near limit points). The elastic limit is reached only via large `η` (β→0 with finite dt), which is the correct knob. The wrapper applies `-eta` in the Tier-1 path only; under `-implex` the implicit solve stays inviscid (combined mode deferred, warned).
- **Gate-construction trap (the expensive one):** building a plastic committed state by uniaxial-**STRAIN** compression drives the numpy ORACLE into its known deep-compression apex regime (`kp<0`, the return is spuriously judged ELASTIC, `gap=σ̄_tr−σ̄_inv=0`). A viscous test built that way passes **tautologically** — the relaxation has nothing to relax (trial==inviscid). The eta gates (oracle PV5b/PV6, g++ NETA non-tautology check) therefore build the plastic state via the uniaxial-**STRESS** driver `drive_damaged_unified` (replay its solved `eps_lat` as a diagonal 6-strain path), which stays on the compressive meridian and yields a genuine plastic `gap`. Same family as the damaged-oracle "use STRESS-controlled, not confined-STRAIN" rule (#287). The g++ B6 block asserts `max viscous-inviscid gap > 1e-3` so an all-elastic (tautological) fixture can't pass. Learned 2026-06-20, `-eta` C++ port + wrapper, [[31_ladruno_concrete3d_adr]].
### LadrunoConcrete3D P2f: the real CDPM2 beta_c (Eq.50) makes compression MUCH more ductile; omega_c heals on unload (implicit re-solve); C2's Gc integral needs an analytic exp tail
- **beta_c is NOT ~1 in monotonic compression (it's ~0.058).** The P2b/P2c monotonic slice dropped beta_c (=1); the real beta_c = ft*qh2*sqrt(2/3)/(rho_bar*sqrt(1+2Df^2)) ≈ ft/(fc*sqrt(1+2Df^2)) ≈ 0.058 at the uniaxial-compression peak (ft/fc=0.1, Df=1). It scales the PLASTIC-strain part kappa_dc1 (Eq.48), suppressing it ~17×, so compression is MARKEDLY more ductile — the post-peak nominal stress differs by ~23 MPa (fc=30) vs beta_c=1 (real beta_c softens to ~0 only by ~-0.3 strain, vs ~-0.05 for beta_c=1). This is the faithful CDPM2 behaviour (user decision 2026-06-20: restore beta_c always). **Gotcha for gate design:** at DEEP strain omega_c saturates to ~1 for BOTH beta_c values, so the wc GAP is tiny there — use the post-peak STRESS gap (tens of MPa, unambiguous) as the non-tautology metric, not the wc gap. Also the damage-ONSET region is step-size sensitive (a coarse step overshoots the sharp onset, over-accumulating ‖Δεp‖ → more early damage); use a resolved path (~5000 steps over -0.3) for compression-damage gates.
- **omega_c HEALS on elastic unload (the cyclic gap).** The oracle solves omega_c IMPLICITLY against the CURRENT effective stress every step ((1-wc)*D = fc*exp(-eps_i/eps_fc), eps_i carries wc). On an elastic UNLOAD the effective stress D drops with the history (kdc1/kdc2) frozen, so the solved wc DECREASES — the material spuriously heals and the reload modulus returns to ~E. beta_c (the damage RATE factor) does NOT fix this; a cyclic-correct omega_c must be driven by the MONOTONE history (omega_c <- max over the path), a separate fix touching every driver + the committed state + the C++ kernel (the P2f monotone-omega_c slice). Reported as the run_p2f_gate F4 diagnostic, not gated.
- **C2's "by construction" Gc integral truncates with the ductile beta_c — add the analytic exp tail.** The omega_c solve enforces |sig_c_nom| = fc*exp(-eps_i/eps_fc) EXACTLY, so int_0^inf |sig| d eps_i = fc*eps_fc = Gc/lch. Over a FIXED compression strain the more-ductile beta_c response doesn't reach eps_i->inf at small lch (large eps_fc), truncating the integral (4.58 vs 5.0 at lch=50). FIX: add the analytic tail int_{eps_i,last}^inf = eps_fc*|sig_c_last| to the trapezoidal path integral — exact (5.000 at all lch), no need for an absurd strain range. **The analytic damaged tangent's ∂beta_c/∂eps** is a composite micro-FD THROUGH the return map (beta_c depends on BOTH rho_bar(sig_bar) AND qh2(kp)) — the C++ kernel `damagedTangent` MIRRORS that composite FD (re-runs `returnMapTensor` per component; needs a forward-decl since returnMapTensor is defined after damagedTangent) rather than deriving ∂kp/∂eps analytically. **beta_c is oracle+C++ in ONE PR (#321), NOT split like eta** — it is always-on (no default-off), so the regenerated g++ fixture's compression cases diverge from a beta_c=1 kernel and the byte-check FORCES oracle+kernel together (like #249); no wrapper/serialization change. Learned 2026-06-20, P2f beta_c oracle + C++ port (#321), [[31_ladruno_concrete3d_adr]].
