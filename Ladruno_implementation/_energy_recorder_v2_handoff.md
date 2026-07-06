---
title: "EnergyBalanceRecorder v2 energy channels (ADR-69) — session handoff"
project: Ladruno
type: handoff / next-session pointer
related:
  - "[[69_ladruno_energy_recorder_channels_adr]]"
updated: 2026-07-06
---

# EnergyBalanceRecorder v2 energy channels (ADR-69) — session handoff

Next-session pointer. Full detail (source-verified mechanism traces, gate
numbers, adversarial-review findings) lives in the ADR's **Implementation
log**: [[69_ladruno_energy_recorder_channels_adr]]. This is the short version.

## Status

| Piece | State | PR |
|---|---|---|
| **P0.5** — sign-pinning (Lysmer injection leaks into `IE` w/ resisting-force sign; outcome (a) confirmed empirically) | ✅ done, no code | — |
| **P1** — `EnergyChannelRegistry`, `-v2` recorder layout, `ExplicitBathe -lnvd` publisher, `LysmerTriangle` Tcl publisher + `eleLoad -lysmerVelocityLoader` (Tcl) | ✅ merged on `ladruno` | [#488](https://github.com/nmorabowen/OpenSees/pull/488) |
| **P1.5** — `LysmerTriangle` + loader openseespy dispatch, `ASDAbsorbingBoundary2D/3D` bottom compliant-base publisher | ✅ merged on `ladruno` | [#491](https://github.com/nmorabowen/OpenSees/pull/491) |
| **P1.6** — ASD **lateral** free-field injection (`addRffToSoil`) publish | ⬜ not started | — |
| **P2** — modal-damping publisher, hourglass energy via `getResponse`, count-where-owned MPI gate | ⬜ not started | — |
| **P3** — mass-scaling fidelity advisory column | ⬜ not started | — |

The recorder today (`recorder EnergyBalance ... -v2`): legacy 6-column output
is byte-identical when the flag is omitted. With `-v2`: `KE_ele KE_nod IE
DW_ele DW_nod ULW [E_inject] [E_lnvd] RES ERR%`. `E_inject` covers **Lysmer**
(full) and **ASD bottom compliant-base** (full) — **not** ASD's lateral
free-field boundary (documented gap, recorder warns when ASD elements are
present). `E_lnvd` covers `ExplicitBathe -lnvd` FLAC dissipation. Modal
damping still pollutes `RES` (P2).

## Where to work

- **Worktree:** `C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\strange-hawking-15f131`
  on branch `guppi/strange-hawking-15f131` (already merged twice to `ladruno`;
  keep using this branch/worktree for P1.6+ unless told otherwise — `git merge
  origin/ladruno` before starting new work, don't rebase).
- **MUMPS junctions already set up** (`mumps-install`/`mumps-archive`/`mumps-src`
  → PowerShell `New-Item -ItemType Junction` from the main checkout) — don't
  recreate them, don't delete them.
- **Build:** `$env:LADRUNO_OPENSEES_QUIET = "1"; cmd /c "Ladruno_scripts\build.bat OpenSeesPy"`
  (or `OpenSees` for the Tcl exe) via the **PowerShell tool**, `run_in_background: true`
  — cold build is ~20-25 min, incremental is faster. Always check the background
  output for `error C` before trusting a run.
- **Run tests against the worktree's own pyd** (NOT the installed/main-checkout
  one — a boot `.pth` will otherwise silently hijack `import opensees` to a
  stale build): `Ladruno_scripts/_run_energy_tests.py` is an UNTRACKED
  (machine-local, not committed — same policy as the ADR-60 dev helpers)
  wrapper that does `python -S` + manual `sys.path`/`os.add_dll_directory`
  wiring + asserts the loaded `__file__` path. Recreate it if missing:
  ```python
  import os, sys
  WORKTREE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  DIST_BIN = os.path.join(WORKTREE, "dist", "bin")
  sys.path.insert(0, r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
  sys.path.insert(0, DIST_BIN)
  os.add_dll_directory(DIST_BIN)
  import opensees
  assert DIST_BIN.lower() in os.path.normpath(opensees.__file__).lower()
  import pytest, sys
  sys.exit(pytest.main(sys.argv[1:] or [os.path.join(WORKTREE, "tests", "test_energyBalanceRecorder.py"), "-x"]))
  ```
  Run via `"C:/Users/nmora/AppData/Local/Python/pythoncore-3.12-64/python.exe" -S Ladruno_scripts/_run_energy_tests.py [files...]`.
- **Full-battery ordering matters, isolated runs don't catch it.** Channel
  declarations in `Ladruno::EnergyChannelRegistry` are **process-sticky**
  (monotone, no reset). A test that hard-codes a column count or indexes a
  channel from the back will pass alone and fail in the real CI ordering if
  an earlier-sorting test file declares a channel first. **This bit twice
  already** (once in P1, once in P1.5's own new test). When adding a new
  `-v2` test: locate a channel by its **guaranteed write-order position**
  (front-anchored if it's unconditionally declared by your own model, e.g.
  `E_inject` is always col 7 whenever `chInject` is true; back-anchored only
  for the genuinely-always-last channel, `E_lnvd`/`RES`/`ERR%`) — never by an
  assumed total column count. Before trusting a new test, run it via
  `_run_energy_tests.py tests/test_adr52_w1e2_explicitbathe_collapse.py tests/test_energyBalanceRecorder.py`
  (that ordering declares `LNVD_WORK` first and is what caught both bugs).
- **Gate evidence lives in `Ladruno_implementation/energy_v2/`**: `p05_*`
  (sign-pin), `p1_ricker_closure.py` (Lysmer, Tcl), `p15_asd_ricker_closure.py`
  (ASD, openseespy) — all runnable standalone against the worktree build.

## P1.6 starting point (ASD lateral free-field injection)

`addRffToSoil()` (2D: `ASDAbsorbingBoundary2D.cpp` ~line 1416; 3D: `~line
2291`) transfers the free-field column's response into the soil domain —
physically a source term like `addBaseActions`, but unlike it, its magnitude
depends on the free-field column's own internal kinematic state (built up
over the whole run via `addKff`/`addMff`/`addCff`), not just a time-series
lookup. That's **why P1.5 scoped this out**: recomputing it "for free" at
commit time is not obviously safe the way `addBaseActions` was — need to
verify whether calling it a second time (outside the normal
`getResistingForce`/`getResistingForceIncInertia` call sites) reads any
state that could have already advanced, before assuming the same "recompute
into a scratch Vector, dot with getVelocity()" pattern applies unchanged.
Start by reading `addRffToSoil` end-to-end (both dims) the way P1.5 read
`addBaseActions`, then re-run the mutual-exclusion argument in reverse (does
publishing this term double-count anything `addRff`/`addBaseActions`
already covers?).

## P2 starting points

- **Modal damping publisher:** site is `TransientIntegrator.cpp`
  (`modalDampingFactors`, application inside the solve) — not yet
  read end-to-end this line; start there.
- **Hourglass via `getResponse`:** preferred zero-vanilla route identified in
  the ADR (§ Mechanism C) — the recorder resolves `Response*` objects once in
  `buildCache()`, fork elements (`LadrunoBrick`) add an `"energy"` response.
  Not started.
- **Count-where-owned MPI gate:** needs a `Node` owner/rank query that
  doesn't exist yet (`Node.cpp` has no `processID` field) — the harder,
  more invasive item; likely wants its own scoping pass before code.
