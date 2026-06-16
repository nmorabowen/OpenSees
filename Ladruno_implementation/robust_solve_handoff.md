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
- **Layer-1.5 C++ getters are DRAFTED but NOT YET BUILT** — committed on
  `guppi/robust-solve-driver` (commit `86af9af37`), not pushed, not compiled.
  They unlock rung-4 stabilization. Building them is the next step.

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
- Draft, unmerged, **UNBUILT** (commit `86af9af37`):
  - `SRC/analysis/integrator/LadrunoArcLength.{h,cpp}` — `getStabilizationDissipatedEnergy`,
    `getReferenceStrainEnergy`, `getStabilizationDissipationRatio`, `scaleCVisc`.
  - `SRC/analysis/integrator/LadrunoDynamicRelaxation.{h,cpp}` — `getResidualNorm`,
    `getKineticEnergy`.
  - `SRC/interpreter/OpenSeesCommands.cpp` — `ladrunoArcLength` subcommands
    `dissipationRatio | dissipatedEnergy | referenceEnergy | scaleCVisc`.
  - No `sendSelf/recvSelf` change was needed (read-only on existing state).

## Resume recipe (next session)

```
# 1. operate in the robust-solve worktree (NOT the main one)
cd C:\Users\nmb\Documents\Github\OpenSees-robust
git fetch origin
git rebase origin/ladruno        # Layer-0 merged; leaves only 86af9af37 ahead

# 2. BUILD (only when the build machine is free — the beam agent's heavy
#    oneAPI build must not be running; two concurrent builds thrash).
#    The worktree has no warm build/ tree, so first build is cold (~MUMPS+
#    conan+full compile). Faster alternative: an incremental build tree.
Ladruno_scripts\build.bat OpenSeesPy        # or full; then refresh dist\bin\opensees.pyd

# 3. VERIFY the getters compiled + add the rung-4 stabilization regression
#    (a softening model run under LadrunoArcLength -stabilize; assert
#     ladrunoArcLength('dissipationRatio') < fTarget, and scaleCVisc ramps it).
python -m pytest Ladruno_scripts\robust_solve_tests -q

# 4. push Layer-1.5 to its OWN PR (base ladruno) once green
```

## Immediate next steps (priority order)

1. **Build + verify the Layer-1.5 getters** (above). The only thing blocking
   rung-4.
2. **Rung-4 wiring in `robust_drive.py`**: replace the `stabilize=True` refusal
   with the real reflex — arm `-stabilize`, read `dissipationRatio` each step,
   `scaleCVisc` to ramp down, and the mandatory c-reduction verification pass
   (R-DIFFUSION). Add the **L0 energy partition** post-processing
   (`SW := W_stab`, `RES_true := RES − SW`) — ADR-31 "Energy accounting".
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
