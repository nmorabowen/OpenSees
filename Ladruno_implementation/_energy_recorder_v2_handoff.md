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
| **P1.6** — ASD **lateral** free-field injection (`addRffToSoil`) publish | ✅ merged on `ladruno` | [#496](https://github.com/nmorabowen/OpenSees/pull/496) |
| **P2** — `E_modal` publisher, `E_hg` hourglass pull, `-ownedNodes` MPI gate + per-rank files | ✅ done (this session's PR) | — |
| **P3** — mass-scaling fidelity advisory column | ⬜ not started | — |

The recorder today (`recorder EnergyBalance ... -v2`): legacy 6-column output
is byte-identical when the flag is omitted. With `-v2`: `KE_ele KE_nod IE
DW_ele DW_nod ULW [E_inject] [E_lnvd] [E_modal] [E_hg] RES ERR%` (parse by
the echoed column NAMES, never position). `E_inject` covers **Lysmer** and
**ASD** in full. `E_lnvd` covers `ExplicitBathe -lnvd`. `E_modal` covers
modal damping for integrators on the base `IncrementalIntegrator::commit()`
(Newmark; HHT-family overrides = documented gap). `E_hg` is the
recorder-pulled hourglass attribution (`IE_display = IE_raw − L − E_hg`,
RES-invariant; LadrunoBrick today, any element exposing `hourglassEnergy`
tomorrow). `-ownedNodes <regionTag>` filters nodal terms for MPI dedup;
file outputs auto-suffix `part-<rank>` under a launcher. KEY P2 finding:
flat-per-rank split-mass models sum nodal terms CORRECTLY naively — the
multiply-count hazard is full-mirror-convention-only (quirks ledger).

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

## P1.6 outcome (ASD lateral free-field injection — DONE)

The P1.5 worry was **refuted by source verification**: `addRffToSoil` is a
pure, stateless function of `getTrialDisp()` (no accumulated internal state
— `addKff`/`addMff`/`addCff` build the FF column's *dynamics*, not a cached
state the transfer reads), so the recompute-at-commit pattern applied
unchanged as a `BND_BOTTOM ? addBaseActions : addRffToSoil` branch in the
existing publisher. Full mechanism trace + gate numbers in the ADR's
implementation log. Two NEW gate-design quirks ledgered: UniformExcitation
input work hides in the recorder's IE via element load vectors (`F = K*u −
Q`), and ASD's `addInertiaLoadToUnbalance` no-op means FF columns are never
driven by uniform excitation (drive closure gates with initial velocities +
alphaM Rayleigh instead — the lateral FF column is otherwise undamped).
Gate evidence: `energy_v2/p16_asd_lateral_closure.py`.

## P2 outcome (all three items — DONE)

- **E_modal:** published from `IncrementalIntegrator::commit()` (rate
  captured per-iterate in `addModalDampingForce`, trapezoid at commit,
  declare-on-first-publish). COVERAGE: base-commit integrators only
  (Newmark ✓; HHT family / `*_TP` override commit without chaining —
  documented gap, quirks-ledgered).
- **E_hg:** the element side had ALREADY SHIPPED (`LadrunoBrick`
  `hourglassEnergy` response, responseID 8) — P2 is recorder-pull only
  (probe cache in `buildCache`, `chHourglass` from the CURRENT model at
  initialize, no registry). Rebucket uniform with E_inject. Gate needs
  `-precision 12` (near-cancelling identity vs 6-digit file quantization).
  Fast-follow candidate: add the same response to `LadrunoQuad` (ssp Kstab
  exists there); the recorder picks it up by name automatically.
- **MPI:** `-ownedNodes <regionTag>` + auto `part-<rank>` file suffixing;
  Node-owner-field and in-recorder-MPI options REJECTED (flat-per-rank has
  no ownership authority; OPS_Recorder lib compiles without MPI). MEASURED:
  split-mass idiom sums nodal terms correctly naively; the hazard is
  full-mirror conventions only (see quirks). Gate:
  `energy_v2/p2_mpi_owned_nodes.py` (needs `dist/openseesmp`, build target
  `OpenSeesPyMP`).

## P3 starting point (mass-scaling fidelity advisory column)

The last open ADR-69 item. ADR § "Proposed column schema" sketches an
advisory column for mass-scaling KE fidelity (the consistent-SMS energy
conduit already closes the BALANCE — this is about surfacing the added-mass
kinetic-energy fraction as a quality metric, LS-DYNA GLSTAT-style). Read the
ADR schema section + `LadrunoMassScalingEnergy.h` (registry already feeds
the kernel's KE via `elementScalingKE`) before scoping.
