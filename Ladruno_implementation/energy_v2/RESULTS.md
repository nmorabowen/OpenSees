# ADR-68 P0.5 — LysmerTriangle / EnergyBalanceRecorder sign-pinning (2026-07-06)

Diagnostic-only experiment against the **installed baseline** build (no source
changes, nothing compiled). Driver: `p05_lysmer_sign_pin.py` (this dir); full
console capture in `p05_run_output.txt`; per-run Tcl decks, console logs, and
recorder outputs (`*_energy.txt`, `*_vels.txt`) alongside.

## Binaries used

- **Interpreter**: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`
  run with `-S` (the boot `.pth` currently hijacks `import opensees` to a STALE
  worktree pyd: `...\worktrees\vigilant-solomon-8c93bd\dist\bin\opensees.pyd`,
  Jun 22, ADR-39 contact branch — NOT the install).
- **Engine actually exercised**: `C:\Program Files\Ladruno\OpenSees\bin\OpenSees.exe`
  (Jun 25, version 3.8.0, commit `bab6995cc`) — the canonical installer build.
  The matching installed pyd `C:\Program Files\Ladruno\OpenSees\bin\opensees.pyd`
  imports fine under `-S`, **but Python could not be used for the model** (below),
  so the experiment runs Tcl decks through the installed exe from a Python harness.

## Why Tcl, and why the loader cannot be driven at all

Grep evidence (paths relative to `C:\Users\nmora\Github\OpenSees_Compile\OpenSees\SRC`,
identical in the installed commit `bab6995cc`):

1. **`LysmerTriangle` is Tcl-only.** Registered in
   `element/TclElementCommands.cpp` and `runtime/commands/modeling/element.hpp`,
   but ABSENT from the Python map in `interpreter/OpenSeesElementCommands.cpp`
   (247 `functionMap.insert` entries, no Lysmer). Syntax
   (`element/absorbentBoundaries/LysmerTriangle.cpp`, `OPS_LysmerTriangle`):
   `element LysmerTriangle eleTag iNode jNode kNode rho Vp Vs <eleLength> <stage>`
   — default `stage 0` = pure damping. (Parser quirk: passing only ONE optional
   arg makes the stage read fail silently; pass none or both.)

2. **`LysmerVelocityLoader` has NO interpreter invocation — Tcl or Python.**
   `grep -rn "LysmerVelocityLoader" SRC` hits only: `classTags.h:651`
   (`LOAD_TAG_LysmerVelocityLoader 17`), `domain/load/LysmerVelocityLoader.{h,cpp}`
   (the class), `element/absorbentBoundaries/LysmerTriangle.cpp:426` (the
   `addLoad` consumer), and Makefile/CMake object lists. The Tcl `eleLoad` types
   (`modelbuilder/tcl/TclModelBuilder.cpp`) are: beamUniform, beamPoint, BrickW,
   surfaceLoad, selfWeight, shellThermal, ThermalWrapper, beamThermal, beamTemp.
   The Python set (`interpreter/OpenSeesPatternCommands.cpp`) adds IGAFollowerLoad.
   No Lysmer anywhere; `LysmerTriangle::setParameter` exposes only
   stage/rho/Vp/Vs (no gnd_velocity back door). The loader is C++-only dead code.

3. **Runtime proof (runs B, B2).** `eleLoad -ele 101 102 -type -lysmerVelocityLoader 1`
   (and `-LysmerVelocityLoader`, `-lysmer`, even `-fooBarBazNotALoad`) all return
   TCL_OK **silently** — the classic eleLoad handler's tail is
   `// if get here we have successfully created the load` → `return TCL_OK;` for
   any unmatched `-type`. Run B2 drives the model with ONLY the loader eleLoad
   under a Ricker pattern and measures `max |vx| = 0.000e+00` over all nodes and
   600 steps: **silent no-op confirmed**.

### Surrogate used to pin the sign empirically

The suspected leak mechanism is *"work of an elemental load rides
`getResistingForce()` into IE with resisting-force sign, while ULW (from
`Node::getUnbalancedLoad`) never sees it"*. That exact pathway IS reachable via
`stdBrick` body forces: `element stdBrick ... 1 1.0 0.0 0.0` +
`eleLoad -ele 1..10 -type -BrickW` under a Path(Ricker) pattern injects x-body
force whose equivalent nodal work is computed independently in the harness
(trilinear tributary volumes × recorded nodal velocities). Mechanism-identical
to the Lysmer loader term `internalForces = C·v_gnd` (LysmerTriangle.cpp:500-513,
note the `0*v(i) + gnd_velocity(i)` — nodal velocity is deliberately zeroed;
stage-0 `getResistingForce()` returns that stale member).

## Model (all runs)

1×1×10 m soil column, 10 `stdBrick`, `ElasticIsotropic` E=2e8, nu=0.25,
rho=2000 (Vs=200, Vp=346.41); 1-D x-polarized shear (y,z fixed everywhere,
x free); base quad covered by 2 `LysmerTriangle` (stage 0) with same rho/Vp/Vs;
`recorder EnergyBalance -file ... -time`; Newmark 0.5/0.25, dt=1e-3, 600 steps
(0.6 s), algorithm Linear, UmfPack, constraints Plain, numberer RCM. Ricker:
f0=10 Hz, t0=0.15 s. All `analyze` returns 0.

## Results

### Run C — surrogate injection (elemental BrickW Ricker drive, Lysmer base)

Independent injection work `W = Σ∫f·v dt = 7.055e-7 J`. Finals (t=0.6 s):

| KE | IE | DW | ULW | RES | ERR% |
|---|---|---|---|---|---|
| 3.7e-12 | **−7.055e-7** | 3.529e-7 | 0.0 | +3.525e-7 | +15.8 |

- **IE_end / W = −0.9999** → the elemental-load work leaks into IE with
  **negative** sign, magnitude exactly W_inject.
- ULW_end = 0 exactly → ULW is blind to elemental loads.
- RES_end / W = +0.4997, DW_end / W = 0.5003 (see finding F2).

### Run A2 — control, NODAL Ricker pulse on top nodes (clean ULW path)

Finals: KE=8.9e-7, IE=1.7e-6 (both ≈0 ✓), ULW=7.082e-2, **DW=3.545e-2 =
0.5006·ULW**, RES=+3.536e-2 = +0.499·ULW, ERR=+27.4%.
Nodal-load work IS booked in ULW (contrast with run C), IE stays ≥0 and returns
to ~0 ✓ — but half the input energy is never booked as damping work (F2).

### Run A1 — control, initial velocity vx=0.1 on top nodes, no loads

KE 2.984 → 0.042 (98.6% absorbed — radiation works), IE ≥ 0 throughout
(max 2.04, final 0.27), ULW=0, DW_end=1.717, RES −2.03…−3.63. RES≠0 here is
expected from the kernel's design (initial KE exists before any integrated work
— RES(t) = −(KE+IE+DW) in a load-free run) plus F2.

## VERDICT — outcome (a), with a superposed second defect

**Leak sign pinned: (a). leak = −W_inject** (IE_end/W = −0.9999). The recorder's
IE integrates `getResistingForce()`-signed power, and the injection enters the
global residual as an applied force `−F_stale` (FE_Element::addRIncInertiaToResidual
adds resisting force with −1), so IE books exactly minus the injected work while
ULW books nothing. In an energy-consistent absorber this alone drives RES → 0:
**the balance accidentally closes while IE trends strongly negative** — the
dangerous outcome (a): a Lysmer compliant-base run would look "closed" (small
RES) while IE is silently polluted by −W_inject.

The textbook RES≈0 signature is masked in THIS implicit harness by an
independent finding:

- **F2 — Lysmer stage-0 absorption under implicit Newmark books only half the
  absorbed energy.** Both A2 and C show DW_end = 0.500×(input work) with
  KE,IE→0: the discrete trajectory loses the other half without any booked
  force work (RES = +0.5·W_in, ERR ≈ +16…27%). Consistent with the element
  never contributing `C·v_node` to the residual (`0*v` at
  LysmerTriangle.cpp:500-502; damping acts only through the c2·C term in the
  Newmark LHS) → the solve itself is energy-inconsistent, and the recorder
  correctly exposes it. Hypothesis, not proven: the ½ factor is the classic
  "damping on velocity increments only" signature.
- **F3 — Tcl `eleLoad` silently accepts ANY unknown `-type`** (returns TCL_OK,
  no warning, creates nothing). Foot-gun independent of ADR-68.

## Implications for ADR-68

1. The sign-pin needed for P0.5 is **(a)**: fix must ensure injection work is
   booked as external work (ULW-like channel) or excluded from IE — NOT relied
   upon to cancel, since it currently *hides* inside a seemingly-closed RES.
2. The velocity loader is unreachable from any interpreter in the baseline —
   any fix + test for the loader path must first add an `eleLoad` hook (and a
   Python element registration if the test is to be openseespy-based).
3. Stage-0 Lysmer + implicit Newmark is energy-inconsistent independent of the
   recorder (F2) — relevant to how ADR-68 words its "RES closes" acceptance
   criteria for absorbing-boundary runs.
