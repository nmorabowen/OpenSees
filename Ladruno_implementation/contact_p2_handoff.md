# ADR-39 ContactDomain — handoff (resume at P2a)

> Handoff for the next session. Everything is durable in git + the loop-state doc.
> Read this, then `_adr39_loop_state.md` (live driver), then `_adr39_p2_design.md`
> (the gated P2 spec). Parent: `39_ladruno_contact_domain_adr.md`.

## Status (2026-06-21)

- **P1 shipped** in PR [#345](https://github.com/nmorabowen/OpenSees/pull/345)
  (base `ladruno`): ContactDomain skeleton + `LadrunoContactHandler` +
  `LadrunoContactFE` + Domain lifecycle hooks + parser. **Zero-force** (plumbing
  only). 7/7 Zone-A green (`tests/test_adr39_contact_p1.py`). Opt-in, byte-identical
  to stock when unused. Verify #345 merged before building on it
  (`gh pr view 345 --json state`) — this fork auto-merges; if it stranded, recover
  per [[feedback_stranded_commits_after_automerge]].
- **P2 designed + gated** (4-agent gate wp3cr60mf, SALVAGEABLE; all fixes folded
  into `_adr39_p2_design.md`). **No P2 code yet.**

## Resume = code **P2a** (rigid/inclined analytical plane, `-kn $val`)

The gate split P2: **P2a first** (isolates penalty + B + assembly + active-set;
NO projection kernel, NO auto-kn, NO master mass). P2a adapter = per slave node,
connectivity `{slave}` only, `B = nᵀ`, `K_c = kₙ n⊗n` on the slave.

**P2a steps:**
1. Study the **custom-`FE_Element`-subtype** precedent before coding the adapter
   connectivity: `SRC/analysis/handler/{PenaltySP_FE,TransformationFE,LagrangeSP_FE}.cpp`
   — how a constraint-handler FE_Element sets its own `myDOF_Groups` + `myID` and
   returns a real tangent/residual. (P1's adapter was empty-connectivity; P2a needs
   the slave node's DOF_Group bound at handle()-time.)
2. `LadrunoContactDomain`: add a **rigid-plane** contact def (point `p0` + outward
   normal `n̂` + `kn` + slave-surface tag) + per-slave pair state (`gₙ0`, active flag).
3. Handler: build one adapter per slave node of the contact's slave surface, passing
   the slave node's DOF_Group + the plane + the ContactDomain ptr (for state).
4. `LadrunoContactFE`: real connectivity = slave DOFs; `getResidual` = `tₙ=kₙ⟨−gₙ⟩₊`,
   `gₙ = n̂·(x_s − p0)`, force `−tₙ n̂` on slave; `getTangent` = bare `kₙ n̂⊗n̂`
   (Newmark scales by c1 — NO internal pre-multiply).
5. Parser: extend `contact` (or a `contactPlane` cmd) to define the rigid plane.
6. Build `OpenSeesPy`, test `tests/test_adr39_contact_p2a.py`.

**P2a acceptance (gate-mandated):** penetration `g=P/kₙ` (rel 1e-8, explicit +
implicit); inclined plane `g=P(n̂·d̂)/kₙ`; release/reopen → exact F=0; sign
(penetration→restoring, NEVER attract — use the Macaulay bracket `⟨−gₙ⟩₊`, never
bare `−kₙgₙ`); explicit 500-step stable; implicit converges with the **anti-chatter
rule** (freeze active set per step, or `algorithm NewtonLineSearch`); E_contact
load–unload → 0.

## P2b (after P2a) — the deformable mechanics, where the BLOCKERs live

Faceted-master projection rung → 2 deformable `LadrunoBrick` blocks + Hertz +
`-kn auto` + SOFT floor + the **∂n/∂u tangent**. Honor the gate BLOCKERs:
- **Derive the outward normal from the master ELEMENT centroid**, never trust node
  winding (`n·(x̄ − x_elem_centroid) > 0`). Gate: winding-flip → identical force.
- **Bounded projection Newton**: cap 10 iters + `|detK|<ε‖g1‖‖g2‖` guard +
  non-convergence **sentinel** (scan skips, never assembles garbage). Do NOT copy
  SimpleContact3D's unbounded `while` (`:600-635`).
- **Tangent = `kₙ BᵀB` + ∂n/∂u block** (symmetric for frictionless), drop only
  `O(gₙ·κ)`. FD-on-flat-ROTATED must PASS; FD-on-curved carries an `O(gₙ·κ)`
  residual that shrinks with refinement.
- Concave-corner tie-break (max-penetration in-bounds, else edge/vertex avg normal).
- `getID` superset ENFORCED (immutable surface + subset assertion in getResidual).
- Implicit RE-PROJECTS per Newton iteration; explicit freezes per step.
- **Oracle**: hand-placed `ZeroLengthContactASDimplex` pre-defined pair, rel 1e-6
  (THE deformable cross-check); Hertz = refinement-convergence benchmark.
- Reusable math → header-only OpenSees-free `SRC/domain/contact/LadrunoContactKernel.h`
  (mirror `LadrunoJ2Kernel.h`). `ContactMaterial3D` has NO penalty tangent to lift.
- `-kn auto` = `f_si·K·A²/V` (26.14a) + SOFT floor (26.15): K/A/V have **no Element
  API** — cache at setDomain from `material->getInitialTangent()(0,0)` of a master
  GP + V from nodal coords; retain the master element pointer.

## Later: P2.5 bucket sort · P3 IMPL-EX Coulomb friction (SHIP) · P3.5 implicit
friction · P4 SOFT · P5 segment-based · P6 tied (rides ADR-30 projection).

## Environment / build / test (verified, don't rediscover)

- Run python (numpy+pytest): `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`.
- DIST (worktree-local, built here): `<worktree>\dist\bin`. Bootstrap:
  `os.add_dll_directory(DIST); sys.path.insert(0, DIST); import opensees`.
- Build (run via `cmd /c`, PowerShell tool, NOT Bash): `Ladruno_scripts\build.bat OpenSeesPy`
  (incremental ~min after the first cold build; full `OpenSees OpenSeesSP OpenSeesMP`
  for a PR-parity build). build.bat roots at its own dir → worktree build, no warm
  cache shared with the compile-root checkout.
- Run a test file:
  `python -c "import os,sys; D='<wt>/dist/bin'; os.add_dll_directory(D); sys.path.insert(0,D); import pytest; sys.exit(pytest.main(['tests/test_adr39_contact_p2a.py','-v']))"`
  (forward slashes in the path string — backslashes break the `\U` escape).
- New authored files: add to `Ladruno_scripts/stamp_headers.py` GLOBS + run it.
  Ledger every new file (LEDGER_implementations) + every vanilla edit (LEDGER_vanilla_files).
- Command classifier had a transient outage this session — if Bash/PowerShell get
  gated ("cannot determine safety"), it's infra; retry, or do read-only/file work.

## Loop discipline (keep it)
Per phase: design → **adversarial gate at critical junctions** (P2/P3 = full
multi-agent; mechanical ports = light) → code → build → test → code gate → ledger
→ commit. Keep `_adr39_loop_state.md` current every iteration. Commit at each
phase boundary so the loop survives context resets.
