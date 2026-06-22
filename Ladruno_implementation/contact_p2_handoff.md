# ADR-39 ContactDomain — handoff (resume: P2a code gate → P2b)

> Handoff for the next session. Everything is durable in git + `_adr39_loop_state.md`
> (the live driver — read its iteration log 1–11 for the full trail). Gated specs:
> `_adr39_p2_design.md`. Parent: `39_ladruno_contact_domain_adr.md`.

## Status (2026-06-22)

- **P1 MERGED** (#345, ladruno squash cfd35f458): ContactDomain skeleton + handler
  + lifecycle hooks + parser. ZERO-FORCE. 7/7 green.
- **P2a in PR #346** (branch `guppi/contact-p2a`, base ladruno): **rigid analytical
  plane, penalty NTS, frictionless — FIRST LOAD-BEARING contact.** 3/3 green (10/10
  with P1 regression). Mechanics verified EXACT (static P/kn rel-1e-6; explicit
  impact pen=v·√(m/kn)=0.002 exact, restitution e=1.0; inclined off-axis). Verify
  #346 merged (`gh pr view 346 --json state`) before stacking P2b; if it stranded,
  recover per [[feedback_stranded_commits_after_automerge]]. **Branch FRESH off the
  updated ladruno** for P2b (stacked-PR pitfall).

## FIRST next-session task: the P2a CODE gate (user will run)

Run an adversarial code review of the P2a C++ (the user deferred it to fold into
the P2b gate). Focus the reviewers on:
- The `getTangent → integrator->formEleTangent` routing + the virtual `addKtToTang`/
  `addMtoTang`/`zeroTangent` overrides (the explicit-mass-pollution fix). Is it
  correct for ALL integrators (statics, Newmark, HHT, the SMS/Bathe family)? Does
  `c1·K_c` match the consistent dynamic tangent?
- Sign/active-set: `getResidual = −kn⟨−g⟩₊ n`, tangent `kn n⊗n` — never attracts;
  active check consistent between getResidual and getTangent within a Newton iter.
- The **element-free teardown crash** (known limitation, see below): decide the
  real fix — make `FE_Element::theMatrices/theVectors` protected (1-line vanilla,
  ledger) so the adapter can guard-allocate, OR document the ≥1-element requirement.
- Lifecycle: bound adapter holds `Node*` resolved at handle(); survives
  destroy/rebuild across domainChanged + ops.wipe (P1 lifecycle hooks).
- The negligible z-anchor truss in the tests — is it truly zero-perturbation, or
  should the gate prefer a cleaner test model?

## P2a known limitation (carry into the gate)

A truly **element-free** contact model (single mass on a rigid plane, NO structural
elements) **crashes on teardown**: the base `FE_Element` shared static scratch
(`theMatrices`/`theVectors`) is **private** and only allocated by an element-backed
FE ctor, so a contact-only model leaves it null → the base `~FE_Element` null-derefs
it when the last FE (a contact adapter) is destroyed (numFEs→0). The P2a tests
sidestep it with a negligible z-anchor truss (tiny EA, no z-motion → zero physics
effect). Realistic meshed-body contact always has elements, so this is off the v1
path — but the gate should pick the fix.

## Then: P2b — deformable mechanics (the design-gate BLOCKERs live here)

`_adr39_p2_design.md` has the full spec. Faceted-master projection rung → 2
deformable `LadrunoBrick` blocks + Hertz + `-kn auto` + SOFT floor + the **∂n/∂u
tangent**. Honor the gate BLOCKERs verbatim:
- **Derive outward normal from the master ELEMENT centroid** (`n·(x̄−x_elem)>0`);
  winding-flip → identical force (gate). A flipped normal silently passes contact through.
- **Bounded projection Newton**: cap 10 + `|detK|<ε‖g1‖‖g2‖` guard + non-convergence
  SENTINEL (scan skips). Do NOT copy SimpleContact3D's unbounded `while` (:600-635).
- **Tangent = `kn BᵀB` + ∂n/∂u block** (symmetric for frictionless), drop only
  `O(g·κ)`. FD-on-flat-ROTATED must PASS (∂n/∂u≠0 for a flat seg under rotation —
  the whole point); FD-on-curved carries an O(g·κ) residual that shrinks w/ refinement.
- Concave-corner tie-break; getID superset ENFORCED (immutable surface + assert);
  implicit RE-PROJECTS per Newton iter (explicit freezes per step); implicit
  anti-chatter (freeze-set-per-step or NewtonLineSearch).
- **Oracle**: hand-placed `ZeroLengthContactASDimplex` pre-defined pair, rel 1e-6
  (THE deformable cross-check); Hertz = refinement-convergence benchmark.
- Reusable math → header-only OpenSees-free `SRC/domain/contact/LadrunoContactKernel.h`
  (mirror `LadrunoJ2Kernel.h`). `ContactMaterial3D` has NO penalty tangent to lift.
- `-kn auto = f_si·K·A²/V` (26.14a): K/A/V have NO Element API → cache at setDomain
  from `material->getInitialTangent()(0,0)` of a master GP + V from nodal coords;
  retain the master element pointer.

## Then: P2.5 bucket sort · P3 IMPL-EX Coulomb friction (SHIP) · P3.5 implicit
friction · P4 SOFT · P5 segment-based · P6 tied (rides ADR-30 projection).

## Hard-won C++ lessons (do not rediscover)

1. **A backing-element-less FE_Element (myEle==0) MUST route `getTangent` through
   `integrator->formEleTangent(this)`** and override the virtual `addKtToTang`/
   `addMtoTang`/`zeroTangent` — returning K directly from `getTangent` corrupts the
   EXPLICIT mass matrix (CDL re-forms the tangent each step and assembles getTangent
   into M → M_xx≈kn → a≈0 → node coasts through).
2. Custom FE connectivity (PenaltySP_FE pattern): `FE_Element(tag, 1, ndof)` +
   `myDOF_Groups(0)=node->getDOF_GroupPtr()->getTag()`; base `setID()` fills `myID`
   with the node's first `numDOF` equation numbers.
3. Sign: `getTangent`=+K, `getResidual`=force toward the constraint (PenaltySP).
4. Element-free model → base shared scratch null → teardown crash (above).
5. A new classTag needs a `Ladruno_implementation/testbed/manifest.yaml` row in the
   SAME commit (the G9 CI gate), not just the LEDGER_implementations table — P1
   failed CI on this. (P2a added no classTag, so it was fine.)

## Environment / build / test (verified)

- Run python: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe` (numpy+pytest).
- DIST (worktree-local): `<worktree>\dist\bin`. Bootstrap:
  `os.add_dll_directory(D); sys.path.insert(0,D); import opensees`.
- Build (PowerShell tool, NOT Bash): `cmd /c "Ladruno_scripts\build.bat OpenSeesPy"`
  (incremental ~min; `OpenSees OpenSeesSP OpenSeesMP` for PR-parity). Worktree build,
  no warm cache shared with the compile-root checkout.
- Run a test file (forward slashes in the path string — `\U` breaks the escape):
  `python -c "import os,sys; D='<wt>/dist/bin'; os.add_dll_directory(D); sys.path.insert(0,D); import pytest; sys.exit(pytest.main(['tests/test_adr39_contact_p2a.py','-v']))"`
- pytest crashes hide the cause (faulthandler). Reproduce hard crashes with a DIRECT
  `python -u` script + `sys.stderr.write(...,flush)` to find the exact failing step.
- New authored files: add to `Ladruno_scripts/stamp_headers.py` GLOBS + run it.
  Ledger new files (LEDGER_implementations) + vanilla edits (LEDGER_vanilla_files).
- Command classifier may have transient outages; if Bash/PowerShell get gated, retry
  or do read-only/file work.

## Loop discipline (keep it)
Per phase: design → adversarial gate at critical junctions (P2/P3 = full multi-agent)
→ code → build → test → code gate → ledger → commit → PR (base ladruno). Keep
`_adr39_loop_state.md` current every iteration; commit at each phase boundary.
