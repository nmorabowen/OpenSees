---
title: Robust-solve driver — session handoff (cold-resume map)
project: Ladruno
status: in-progress
owner: nmora
updated: 2026-06-16
tags:
  - handoff
  - solver
  - robust-solve
  - integrator
---

# Robust-solve driver — session handoff

Cold-resume map for the robust nonlinear-solution driver track (ADR-31). Read
this first on resume. Design + risks: [[31_ladruno_robust_solve_driver_adr]].

## TL;DR — where we are

- **Layer-0 is MERGED** to `ladruno` (PR #242): the pure-Python `robust_drive()`
  spine (adaptive stepping + algorithm ladder + **rung-3 load→displacement
  constraint switch** + cheap peak detector + JSONL decision log) and the pytest
  acceptance gate (**8 green** on the live build).
- **Layer-1.5 C++ getters are now BUILT + VERIFIED** (2026-06-16). Compiled the
  `86af9af37` getters into the robust worktree's `dist\bin\opensees.pyd` and
  locked them into the battery with a new fixture (`torture_stabilize.py`) + 4
  pytest cases — **12 green** total. Still on `guppi/robust-solve-driver`, NOT
  pushed; the build + tests are uncommitted working-tree changes (see below).
- **A build blocker was fixed in this worktree:** the robust branch's `build.bat`
  hardcoded a stale `PYEXE` (`C:\Users\nmora\…`), so CMake auto-picked Python 3.14
  and configure failed. Ported the machine-agnostic `PYEXE` probe block (the one on
  the beam branch's `18cb1baba`) into `Ladruno_scripts\build.bat` here.
- **rung-4 is still NOT wired** into `robust_drive.py` (`stabilize=True` still
  raises). What changed is the *design*: measuring the seam revealed `-stabilize`
  is a narrow last-resort, not a limit-point hero (see "rung-4 reality check").

## Branch / worktree topology (GET THIS RIGHT before touching git)

| Path | Branch | Contents |
|---|---|---|
| `C:\Users\nmb\Documents\Github\OpenSees` (main worktree) | `guppi/ladruno-dispbeamcolumn` | **THE BEAM AGENT'S branch** — a *different* concurrent agent's uncommitted work (ladrunoDispBeamColumn element, ADR-32). **Do NOT disrupt.** |
| `C:\Users\nmb\Documents\Github\OpenSees-robust` (worktree) | `guppi/robust-solve-driver` | THIS track. The Layer-1.5 draft lives here. |

> [!danger] One shared repo, two agents
> There is a SINGLE git worktree per branch and `git checkout` switches HEAD
> globally. A previous `git checkout` in the main worktree switched the beam
> agent's HEAD out from under it (recovered, nothing lost). **Rule: never run
> `git checkout` in the main worktree while the beam agent is active.** Do all
> robust-solve git work via `git -C C:\Users\nmb\Documents\Github\OpenSees-robust`
> (or future agents: `git worktree add` your own).

## What is merged vs. draft

- Merged via #242 (now on `ladruno`): `robust_drive.py`, `robust_solve_tests/`
  (`torture_softening.py`, `torture_snapthrough.py`, `conftest.py`,
  `test_robust_battery.py`), ADR-31.
- Built + tested, committed at `86af9af37`, **NOT pushed**:
  - `SRC/analysis/integrator/LadrunoArcLength.{h,cpp}` — `getStabilizationDissipatedEnergy`,
    `getReferenceStrainEnergy`, `getStabilizationDissipationRatio`, `scaleCVisc`.
  - `SRC/analysis/integrator/LadrunoDynamicRelaxation.{h,cpp}` — `getResidualNorm`,
    `getKineticEnergy`.
  - `SRC/interpreter/OpenSeesCommands.cpp` — `ladrunoArcLength` subcommands
    `dissipationRatio | dissipatedEnergy | referenceEnergy | scaleCVisc`.
  - No `sendSelf/recvSelf` change was needed (read-only on existing state).
- **Uncommitted working-tree changes in this worktree** (2026-06-16, to commit next):
  - `Ladruno_scripts/build.bat` — machine-agnostic `PYEXE` fix (the blocker).
  - `Ladruno_scripts/robust_solve_tests/torture_stabilize.py` — new fixture.
  - `Ladruno_scripts/robust_solve_tests/test_robust_battery.py` — +4 `test_stabilize_*`.
  - `Ladruno_implementation/LEDGER_quirks.md` — the `-stabilize` reality-check quirk.
  - this handoff.

## rung-4 reality check (measured 2026-06-16 — READ before wiring rung-4)

Building the seam let us actually *measure* what `-stabilize` does. The ADR's
optimistic framing ("arm STABILIZE at a limit point") needs tightening — full
detail in [[LEDGER_quirks]] "what viscous regularization can and cannot pass":

- **Pure softening is unpassable by `-stabilize`** (it IS load control; no
  equilibrium above the peak). It stalls at the peak exactly like plain load
  control. Softening → rung-3 (already shipped), never rung-4.
- **`-adaptStab` prevents crossing** a hard limit (it weakens `c` to hold
  `fTarget`). Crossing a snap-through needs `-stabilize` *without* adaptStab at an
  **elevated** `f` (≈1e-3…1e-2), and even then it's a **diffusive crawl / dynamic
  jump** (R-LOG-MASK), not a traced path. Non-monotonic in `f` (too much damping
  re-freezes the step).
- **`-cVisc` is overwritten** by the first-commit calibration (latent ADR-20 bug).
- **openseespy raises `OpenSeesError`** on any command stderr warning, so the
  driver must pre-validate `scaleCVisc` factor > 0 (cannot catch the −1 return).

**Implication for rung-4 wiring:** rung-4 is a *narrow last resort* — appropriate
only when rung-3 is unavailable (no single control DOF: distributed buckling,
contact chatter) AND a continuing branch exists. It must run **without** adaptStab
at an elevated `f`, gate hard on `dissipationRatio` (R-DIFFUSION), and stamp
`verdict="regularized"`/`"unverified"` (never `"equilibrium"`) unless the mandatory
c-reduction drift check passes. This is a genuine design decision — confirm scope
with the track owner before implementing (the optimistic "stabilize beats the
limit point" reflex would be wrong for the current battery, where rung-3 wins).

## Resume recipe (next session)

```
# 1. operate in the robust-solve worktree (NOT the main one)
cd C:\Users\nmb\Documents\Github\OpenSees-robust

# 2. the build is already done; the seam pyd is at this worktree's dist\bin.
#    To REBUILD after a C++ edit (build machine must be free — no concurrent
#    beam-agent oneAPI build): set PYEXE explicitly to dodge a stray py3.14,
#    then build only the Python target and let it refresh dist\bin\opensees.pyd.
set PYEXE=C:\Users\nmb\AppData\Local\Programs\Python\Python312\python.exe
Ladruno_scripts\build.bat OpenSeesPy

# 3. run the battery against THIS worktree's dist (not the main one)
set LADRUNO_DIST=C:\Users\nmb\Documents\Github\OpenSees-robust\dist\bin
set LADRUNO_OPENSEES_QUIET=1
python -m pytest Ladruno_scripts\robust_solve_tests -q     # 12 green

# 4. commit the working-tree changes, then push Layer-1.5 to its OWN PR
#    (base ladruno). Rebase on origin/ladruno first if it has advanced.
```

## Immediate next steps (priority order)

1. **Commit the uncommitted working-tree changes** (build.bat fix, the new
   fixture + 4 tests, the quirk, this handoff) and push Layer-1.5 to its own PR
   (base `ladruno`). The getters are built + 12-green; nothing blocks the PR.
2. **Decide rung-4 scope with the track owner** (see "rung-4 reality check"):
   the optimistic "arm STABILIZE at the limit point" reflex is wrong for the
   current battery (rung-3 wins on both fixtures). If proceeding, rung-4 wiring in
   `robust_drive.py` = replace the `stabilize=True` refusal with: arm `-stabilize`
   **without** adaptStab at an elevated `f`, pre-validate `scaleCVisc` factor>0,
   read `dissipationRatio` each step, ramp `c` down once past the limit, the
   mandatory c-reduction verification pass (R-DIFFUSION), and a hard
   `verdict="regularized"/"unverified"` stamp. Add the **L0 energy partition**
   post-processing (`SW := W_stab`, `RES_true := RES − SW`) — ADR-31 "Energy
   accounting". A snap-through (not softening) fixture is the right regression
   vehicle, but note it crosses only by a diffusive crawl.
3. **Rung-5**: the `ladrunoDR` runtime command (4-file Python/Tcl wrapper,
   deferred per R-SCOPE) + the implicit→DR handoff primitive (R-HANDOFF), then
   wire the dynamics fall-through.
4. **Battery**: buckling-column fixture (rung-3 limit-point, no build), then the
   slack-cable snap (needs rung-5).

## Key facts — do NOT re-derive

- **Gate vs audit** (the load-bearing decision): the dissipation signal comes
  from the **integrator getter**, never the energy recorder — `f_v`/`M*` never
  touch the domain the recorder sweeps. The recorder is the *audit* and sources
  `W_stab` from the same getter. (ADR-31 "Energy accounting".)
- **`RES = W_stab`**: on a stabilized run the recorder's existing closure
  residual already equals the artificial dissipation; the fix is a *partition*
  (one scalar from the getter), not a new measurement.
- **No `sendSelf` change** to the shipped 33004 class was needed (open-Q#2
  resolved safely).
- **Snap-through ≠ clean killer**: adaptive LOAD control can *dynamically jump* a
  snap-through's unstable branch (it "succeeds" while skipping the true path) —
  softening (no far branch) is the unambiguous load-control killer. The driver
  must flag branch-jumps (R-LOG-MASK).
- **rung-4 is REFUSED until built**: `robust_drive(stabilize=True)` raises by
  design (R-OBS — fail loud, not diffuse silent).

## Build / run environment (this machine)

- Full toolchain installed (VS2022 + Intel oneAPI + CMake/Ninja/Conan/Py3.12)
  via `Ladruno_scripts\install_prereqs.ps1` (PR #238). Build: `build.bat`.
- Python: `C:\Users\nmb\AppData\Local\Programs\Python\Python312\python.exe`.
- Tests resolve the built module from `LADRUNO_DIST` (default
  `C:\Users\nmb\Documents\Github\OpenSees\dist\bin` — the MAIN worktree's dist;
  set it to the robust worktree's dist after building there). Banner off with
  `LADRUNO_OPENSEES_QUIET=1`.
