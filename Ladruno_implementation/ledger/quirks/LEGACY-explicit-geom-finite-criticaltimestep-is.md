---
wp: LEGACY
title: "Explicit -geom finite: criticalTimeStep() is reference-config (must margin dt), and the EnergyBalance recorder reports IE with a flipped sign"
legacy_seq: 64
---
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
