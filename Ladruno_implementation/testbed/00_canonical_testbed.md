# Canonical Test Bed — Ladruno OpenSees fork

> Source of truth for how every extension to this fork is tested. Settled
> after a 5-lens design review (2026-05-30). The design's *shape* (two zones,
> tiered FEM checks) was confirmed; this document fixes the **enforcement and
> the exact numbers** the review found missing.

The meta-finding that drove every decision: **discipline-only fails, and this
repo proves it.** `SRC/classTags.h` has already drifted 205 lines from
`DEVELOPER/core/classTags.h`, carries collision hacks (`ModElasticBeam3d=41234`,
`E_SFI_MVLEM_3D=259259`), and ships `exit(-1)` inside `FourNodeQuad::recvSelf`.
A checklist nobody's CI runs is how all three got in. So every rule below that
*can* be a machine check **is** one.

---

## 0. The prime directive

Extend OpenSees heavily, but keep each change **additive and localized** so
upstream `master` could later accept it. The test bed serves that directive:
a feature's tests must (a) prove correctness to *you*, and (b) travel with the
PR and run on the maintainer's CI with **zero extra setup**.

---

## 1. Two zones

| | **Zone A — upstreamable** | **Zone B — fork-local** |
|---|---|---|
| Lives in | `tests/` (pytest), `tests/tcl/` + `EXAMPLES/verification/` (Tcl) | `Ladruno_files/testbed/` |
| Deps | OpenSees / OpenSeesPy **only** | apeGmsh, gmsh, LS-DYNA refs, numpy/scipy, perf baselines |
| Travels with | the PR to master | **never** offered to master |
| Runs | every push to `ladruno` **and** upstream CI | nightly, self-hosted only |
| Marker | `@pytest.mark.zone_a` | `@pytest.mark.zone_b` |

**Rule:** the upstream PR carries only Zone A. Anything needing a mesher, a
commercial code, or a GPU stays in Zone B and is excluded from the PR diff.

---

## 2. Three axes → tiers

### Axis 1 — C++ architecture debt (is the diff upstreamable?)

Not a runtime test. A **diff-shape discipline backed by two automated gates.**

Per-feature contract for a new Element / Material / Integrator:

- [ ] classTag added in the **ladruno private band ≥ 33000** (see §5) —
      `ci/check_classtags.py` enforces uniqueness + cross-file identity.
- [ ] registered in `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`
      (`getNewElement`/`getNewUniaxialMaterial`/… switch).
- [ ] full virtual interface: `getTangentStiff` · `getInitialStiff` ·
      `getResistingForce` · `commitState` · `revertToLastCommit` ·
      `revertToStart` · `sendSelf` · `recvSelf` · `setResponse` · `getResponse`.
- [ ] **`recvSelf` returns `-1` on error — never `exit()`** (the
      `FourNodeQuad.cpp:1178,1189` anti-pattern kills the Python/Jupyter kernel).
- [ ] **both** dispatch surfaces wired (they are separate, hand-maintained files):
      `SRC/element/TclElementCommands.cpp` (parser) **and**
      `SRC/interpreter/OpenSeesElementCommands.cpp` (`functionMap.insert`).
      Wire one, forget the other → "works in pytest, missing in Tcl."
- [ ] CMake: `target_sources(OPS_Element PRIVATE new.cpp PUBLIC new.h)` in the
      element's dir.
- [ ] builds clean on the matrix: Ubuntu + Mac (upstream) **and** Windows oneAPI.

**Diff-shape budget.** Count touched *upstream core* files. A clean additive
element is `classTags.h` + `FEM_ObjectBrokerAllClasses.cpp` + the two command
files + a `CMakeLists.txt`. Editing `Domain.cpp` / `AnalysisModel` / a solver
is a red flag → needs an ADR before it can be upstreamed. Tag every unavoidable
core edit `// LADRUNO: <reason>`. `grep -rn "LADRUNO:" SRC` is your exact
rebase surface against a moving master.

**Why two gates and not just review:** the two defect classes below are
*invisible* to every runtime test and every human reviewer to date —

1. **classTag collision / drift** → `ci/check_classtags.py`
2. **interface completeness + serialization** → the Zone-A Python round-trip
   (`tests/_testbed/roundtrip.py`): `database -save/-restore` + `setResponse`
   smoke. Catches a forgotten `recvSelf` or a missing response branch without
   new build infra. *(A single vendored `doctest` C++ target proving
   bit-identical `sendSelf`/`recvSelf` is deferred until the first element whose
   serialization you genuinely distrust lands — see ADR note at bottom.)*

### Axis 2 — FEM correctness (does the physics work?)

Tiers, cheapest first. **Tolerances are fixed here, not left to the author** —
a patch test at `1e-3` silently hides a B-bar bug; a rank test at machine-eps
*false-fails* because zero-energy eigenvalues come back at ~`1e-9·λmax`.

| Tier | What | Zone | Tolerance (locked) |
|---|---|---|---|
| **T0** | single-element rank / patch / rigid-body | A | see below |
| **T0m** | material-point driver (for new `nDMaterial`/`uniaxialMaterial`) | A | tangent FD rtol `1e-6` |
| **T1** | closed-form (cantilever, SS beam, SDOF period, frame eigen) | A | `isclose` rel_tol `1e-9` |
| **T2a** | fixed-mesh benchmark decks (Cook, MacNeal–Harder, PinchedCylinder) | A | published targets, §4 |
| **T2b** | apeGmsh convergence sweeps (varying density → rate) | B | observed rate `p > 1.8` |
| **T3** | cross-code / complex-geometry (apeGmsh drive, LS-DYNA oracle) | B | case-specific |

**T0 — locked numbers:**
- **Patch test:** prescribe a linear field `u = a + bx + cy (+ dz)` on the
  boundary nodes of a ≥2-element **distorted** patch, interior node free. Assert
  every Gauss-point stress equals the closed-form constant stress at
  **rtol ≤ 1e-9**, atol floor `1e-9·E·|strain|`. Read via
  `eleResponse(tag, 'material', gp, 'stress')`.
- **Rank / zero-energy:** eigen the free element K (unit mass on every DOF →
  eigenvalues of K). Classify `λ_i < 1e-8·λmax` as zero; assert
  `count == n_rbm` (**3 in 2D, 6 in 3D**); fail on `n_rbm+1` (spurious
  hourglass / mechanism). The `1e-8` cut is **relative** and
  conditioning-aware — not machine-eps.
- **Rigid-body:** prescribe pure translation + small pure rotation; assert
  `‖getResistingForce‖ < 1e-8·E·A`.

**T0m — material point:** drive the material directly
(`setStrain`/`getStress`/`getTangent` via a single-element or the material
driver), assert a known stress–strain path **and** finite-difference the
consistent tangent against `getTangent` to **rtol ~1e-6**
(`eleResponse(1,'material',1,'tangent')` for an `nDMaterial`).

**Integrator (e.g. `ExplicitBatheLNVD`) — 4-part battery** (Zone A where
dep-free; this scheme is **Noh–Bathe two-sub-step + FLAC local-non-viscous
damping for dynamic relaxation**, *not* a generic transient integrator, so one
"SDOF period" check is worthless):
1. **Stability bracket:** `dt = 0.95·dt_cr` stays bounded over 1000 steps **and**
   `dt = 1.05·dt_cr` diverges (`EB_NB_STABILITY_FACTOR = 2.0`, ≈2× the
   central-difference `2/ω` limit).
2. **Period elongation:** `(T_num − T_exact)/T_exact ~ O((dt/T)²)`, measured at
   `dt/T = 0.01` and `0.1`.
3. **Energy drift:** undamped free vibration, `|E(t) − E0|/E0 < 1e-3` over 50
   periods — read via the **EnergyBalanceRecorder (classTag 26)**. Pure
   self-consistency, no external dep → **Zone A**.
4. **LNVD consistency:** a 2-DOF chain relaxes to the **exact** linear-static
   solution, unbalanced-force norm `< 1e-8` — proves local damping does not bias
   equilibrium.

### Axis 3 — performance (no regression) — Zone B only

Hosted-runner variance makes absolute timing useless as a gate, so this never
touches upstream-facing CI.

- **Metric:** state-determination **ns/call** — loop one element's
  `getResistingForce` + `getTangentStiff` `1e5–1e6×`. **Not** wall-clock of
  `analyze()`.
- **Environment:** the pinned self-hosted Windows oneAPI box, **MKL threads = 1**,
  warmup discarded, **median-of-7**.
- **Gate:** baseline JSON keyed by `(machine-id, classTag, op)`; **warn at +10%,
  hard-fail at +25%** of the median. Re-baseline only via an explicit
  `perf: re-baseline <reason>` commit editing the JSON — never CI-auto-updated.
- **Complexity shape:** a **ratio** test (assembly time across 1k/2k/4k elements
  grows ~linearly by ratio) catches accidental O(n²). Ratios survive machine
  variance, so this one *may* run near CI; absolute timings may not.

---

## 3. The merge gate (feature-class-dependent)

No merge to `ladruno` until the required tier is green **and** both automated
Axis-1 gates pass.

| Feature class | Required tiers | + automated gates |
|---|---|---|
| Material (`nDMaterial`/`uniaxialMaterial`) | **T0m + T1** | classtag + manifest + roundtrip |
| Continuum / shell element | **T0 + T1 + ≥1 T2 rate test** | classtag + manifest + roundtrip |
| Integrator | **4-part dynamics battery** | classtag + manifest |
| Derivative (section / recorder / transform) | **T0 + T1 floor** | classtag + manifest |

The convergence-rate test (assert observed rate `> 1.8`) is *the* signature that
separates a correct element from a plausible-but-wrong one. It is cheap enough to
gate on. **Use a SMOOTH problem for the rate test** — a manufactured solution with
the exact field on the whole boundary and a consistent body force (rate ~2 for a
quadratic element). **Do not use Cook's membrane for the rate gate:** its
clamped-corner stress singularity caps the displacement rate at ~1.3 (measured for
BezierTri6, 2026-05-30), so Cook is an *accuracy* benchmark (compare a fixed mesh
to a reference), not a rate test. Worked example:
`tests/test_bezierTri6_element.py::test_mms_convergence_rate`.

---

## 4. Concrete specifics the scaffold hard-codes

**Tolerances** — patch `1e-9` · rank `λ<1e-8·λmax` → count 3 (2D)/6 (3D) ·
rigid-body `<1e-8·E·A` · T1 `isclose` rel_tol `1e-9` · material tangent FD
rtol `1e-6` · energy drift `<1e-3`/50 periods · perf warn `+10%`/fail `+25%`.

**Benchmarks — published targets (do not invent):**
- **Smooth MMS (rate gate):** unit square, exact field on the full boundary via
  the element's interpolation, consistent body force → observed displacement rate
  `> 1.8` (BezierTri6 measures ~2.07). This — not Cook — is the rate test.
- **Cook's membrane (accuracy only):** N = {2,4,8,16,32}/side, tip → **23.96**.
  Singularity-limited (rate ~1.3); compare a fixed mesh to the reference, do NOT
  assert a rate on it.
- **MacNeal–Harder:** normalized tip ∈ **[0.95, 1.05]** (rectangular mesh);
  *record* (baseline, not hard-pass) the trapezoidal/parallelogram locking ratios.
- **Scordelis–Lo** 0.3024 · **hemisphere** radial 0.094 · **pinched cylinder**
  `164.24·P/(E·t)` (already in `EXAMPLES/verification/PinchedCylinder.tcl`, 5% tol).
- **1D elastic-bar wave** (Noh–Bathe's own 2013 validation): front arrival < 2%.
- **Bézier / high-order elements:** patch test at **elevated polynomial order**
  (a complete p-th-order element must reproduce a p-th-order field exactly).

**apeGmsh T3 harness rules (encode in the template, else silent false-greens):**
- `apeSees(fem)` does **not** auto-emit loads / masses / fix / sp — **re-declare**
  on the bridge (`ops.fix`, `ops.mass`, `p.load`, `p.sp`). MP constraints **do**
  auto-emit (the asymmetry is the trap → a converging *zero-load* model "passes"
  by matching ~0).
- **`assert_equilibrium(Σreactions + Σapplied ≈ 0)`** is a **universal guard** in
  every T1/T2/T3 case — one line that kills the entire load-omission class.
- Mesh step order: `generate → split_higher_order_lines(frame_pg, 'split')`
  *(only if a 2nd-order continuum forces frame lines to Line3; **incompatible**
  with HingeRadau/HingeMidpoint — re-mesh nonlinear frames at order 1)*
  `→ crack() → renumber(dim=extract_dim, base=1) → get_fem_data(dim=extract_dim)`.
- Beam internal forces → **native HDF5 path only**
  (`ops.h5()` → `OpenSeesModel.from_h5(fem_root='/model')` → `Results.from_native`).
  This build writes **no beam vecxz to MPCO**; `line_force` over MPCO is garbage.
- Zone-B determinism: `pytest.importorskip('gmsh'/'apeGmsh')`,
  `APEGMSH_SKIP_VIEWER=1`, recorder `file=` → `tmp_path` (**never** the
  OneDrive-synced worktree — `.mpco` on a synced FS aborts the kernel at
  `ops.wipe()`), assert on PG-scoped scalars (tip disp, max von Mises, reaction
  sum, period) — never raw gmsh node-id-keyed fields; pin the gmsh version.

**LS-DYNA (deferred — D6):** not an oracle for anything shipping today
(truss/zeroLength/SDOF all have tighter, license-free closed-form answers; a
number-to-number match would false-fail on legitimate scheme differences and
risks a shared-bug false positive). First non-redundant use: a **fully-integrated
brick** (hourglass off, bulk-viscosity `q1=q2=0`) compared only at
`dt ≤ 0.1·dt_cr`. Cite LS-DYNA Theory Manual §central-difference / §bulk
viscosity / §energy balance in the test docstring; never reproduce text.

---

## 5. classTag policy (the live G2 example)

`ELE_TAG_BezierTri6` was originally `272` — **one above** upstream's highest
(`PML3DVISCOUS = 271`). Upstream's next element would take 272 → **silent
collision**, wrong constructor on `recvSelf`/DB-restore, no compile/link/​
serial-runtime symptom. It was the first remediation: moved `272 → 33000`
(2026-05-30), the first real exercise of `ci/check_classtags.py`.

**Policy:** every ladruno-authored classTag uses the **private band ≥ 33000**,
allocated **from the floor inclusive** (`BezierTri6 = 33000`, the next ladruno
tag → `33001`, …). The band is per-family — an `ELE_TAG` and a `MAT_TAG` may both
start at 33000 since the broker dispatches them through separate switches. 33000
is a high, round offset upstream's low sequential allocation (~272 today) is
very unlikely to reach; it is a strong convention, not a hard guarantee (upstream
itself uses `ELE_TAG_ExternalElement = 99990`), so the checker — not the band — is
the real safety net.
`ci/check_classtags.py` (a) asserts no two `#define`s in a TAG family share a
value, (b) asserts the shared symbols in `SRC/classTags.h` and
`DEVELOPER/core/classTags.h` agree, and (c) **warns** on any ladruno tag (marked
`// ... (Ladruno ...)`) below 33000. At upstreaming time, renumber the band tag
to whatever low integer the maintainer assigns — a one-line, greppable change.

> **Migration (done):** `BezierTri6` moved `272 → 33000` on 2026-05-30. All code
> referenced the *symbol* `ELE_TAG_BezierTri6`, so only the `#define`, one header
> comment, and two docs changed — no broker/dispatch edits (the fork's
> `getNewElement` is a stub; see manifest). The element's T0/T1/T2a tests remain
> the pending part of its `remediation` row.

---

## 6. CI map

- **`.github/workflows/build_cmake.yml`** (upstream, untouched): Ubuntu/Mac/​
  win-oneAPI build + `pytest -v tests/` on master. Leave it alone.
- **`.github/workflows/ladruno.yml`** (new, guarded `if: github.repository ==
  'nmorabowen/OpenSees'`): on every push to `ladruno` —
  `build-{ubuntu,mac,win-oneapi}` (clone the existing oneAPI job, don't rebuild
  it) → `pytest -m "zone_a" tests/` → `tclsh runVerificationSuite.tcl &&
  python ci/check_tcl_results.py` → `python ci/check_classtags.py` →
  `python ci/check_manifest.py`.
- **Nightly self-hosted** `[self-hosted, windows, ladruno-perf]`:
  `pytest -m "zone_b" tests/` + perf micro-bench + LS-DYNA cross-checks. Not a
  required check.

---

## 7. Per-feature manifest

The enforced ledger lives at `Ladruno_implementation/testbed/manifest.yaml`.
`ci/check_manifest.py` diffs `SRC/classTags.h` against the base ref, extracts
new ladruno `#define`s, and asserts each has a manifest row **and** a
`tests/test_<feature>.py` (or an explicit `tcl: WAIVED(<reason>)`). One row per
feature; required fields documented in that file's header.

---

## ADR notes (decisions of record, 2026-05-30)

- **D1** pytest always; **Tcl required only on a Tcl command surface**, else
  `WAIVED(reason)`. Both read one shared params table (`tests/params/<f>.json`)
  so they cannot diverge. Doubled cost → ~1.2×.
- **D2** merge gate is **feature-class-dependent** (§3).
- **D3** perf: median-of-7, MKL=1, warn +10 / fail +25, self-hosted; ratio
  shape test may run near CI.
- **D4** Python round-trip + setResponse smoke **now**; single `doctest` C++
  target **later**, only when a distrusted `sendSelf` lands.
- **D5** T2 split: **T2a** fixed-mesh decks → Zone A; **T2b** apeGmsh sweeps →
  Zone B.
- **D6** LS-DYNA deferred to the brick/shell roadmap; energy-closure check is
  Zone A.
