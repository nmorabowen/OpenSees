# ADR-39 ContactDomain — handoff (RESUME: P3 IMPL-EX Coulomb friction — CODE the C++)

> **⇒ RESUME POINT (2026-06-22): P3 friction — design+oracle+gate DONE, write the C++.**
> Branch `guppi/contact-p3-friction` (off ladruno, has P2.5) already carries the P3 design
> commit. The SPEC is `_adr39_p3_design.md` (gate-HARDENED — read it, the `[GATE]` folds are
> law). The math reference is `contact_prototypes/proto_p3_friction.py` (6/6). Driver
> `_adr39_loop_state.md` iter 20 has the full trail.
>
> **What to build (v1 = EXPLICIT ship, FORCE-only — no friction tangent, no IMPL-EX):**
> 1. Kernel (`LadrunoContactKernel.h`): clean Coulomb DIRECT return map returning the
>    NEGATED applied friction force `t_fric=−tT` (‖tT‖=μN at slip) — the sign lives in ONE
>    place. Physically-scaled near-zero-slip guard (stick if ‖tT*‖<tol·μN). IMPL-EX fns kept
>    for P3.5 only.
> 2. Engine (`LadrunoContactDomain`): per-pair state `{gpT[3], gT0[3], engaged}` keyed
>    `(contactTag, slaveTag, GLOBAL segIndex)`; real `commit()`(gpT=gpT_trial)/`revert()`; a
>    live-key-set GC each handle(); cross-contact shared-slave warn.
> 3. Adapter (`LadrunoContactFE`): kt/mu + LAZY engine re-fetch via `Domain*` (wipe deletes the
>    engine); capture `gT0` at first activation; `mu≤0` SHORT-CIRCUITS before any slot touch
>    (byte-identical P2b); friction force added to the residual mirroring the normal block.
> 4. Handler: thread `ct.kt, ct.mu` + the live key-set to the engine for GC.
> 5. Tests `tests/test_adr39_contact_p3_friction.py`: **INCLINE through REAL CDL+addB reading
>    nodal accel from the SOE (NOT the oracle's drive−fric) ⇒ a=g(sinθ−μcosθ) — THE gate that
>    catches the sign bug**; + stick; slip-caps-μN; energy dissipation; `mu=0` frictionless
>    byte-identical to P2b; commit/revert; wipe-then-reanalyze.
>
> **Discipline reminders:** branch off ladruno; **wait for Zone-A GREEN before merging** (this
> repo's `gh pr merge --auto` merges IMMEDIATELY, before Linux CI — #358 merged pre-Zone-A this
> session, was fine but risky); `set LADRUNO_OPENSEES_QUIET=1` before any build (a wired `.pth`
> makes python print the banner → poisons CMake's Python probe); copy `mumps-archive/
> mumps_src.tar.gz` into a fresh worktree before a cold build.
>
> P3.5 (DEFERRED): implicit frictional Newton + the non-symmetric consistent tangent
> (`∂tT/∂gT=(μN·kt/‖tT*‖)(I−n̂⊗n̂)`, `∂tT/∂gN=−μ·kn·n̂`; oracle has it) + IMPL-EX + committed-N
> cap. P2b-2c (DEFERRED): ∂n/∂u normal tangent + Hertz + SOFT Courant floor for auto-kn.

## (historical) earlier resume point was P2b-2b `-kn auto`

> Read this, then `_adr39_loop_state.md` (the LIVE driver — iteration log 1–17 has the
> full trail), then `_adr39_p2b_design.md` (gated P2b spec). Parent ADR
> `39_ladruno_contact_domain_adr.md`. Everything below is durable in git + those docs.

## Status (2026-06-22)

Phased explicit-first NTS penalty contact, each phase design→(gate)→code→build→test→
code-gate→PR. **SHIPPED to `ladruno`:**
- **P1** (#345) — ContactDomain skeleton + handler + lifecycle hooks + parser. Zero-force.
- **P2a** (#346 + gate-fix #350) — rigid analytical plane, penalty NTS, frictionless.
  FIRST load-bearing. Code-gate PASS (B1 `FE_Element` subtype-ctor scratch null-deref fixed
  at source + kn/surface-kind/ndf guards + addKiToTang).
- **P2b-1** (#354) — faceted node-to-segment vs a FIXED master segment. NEW header-only
  `SRC/domain/contact/LadrunoContactKernel.h` (shape, bounded projection Newton, direction-
  oriented winding-immune normal, gap, penalty) + `LadrunoContactFE` SEGMENT mode +
  one-adapter-per-(slave,segment) handler + `contact … -outward`. Oracle-first
  (`contact_prototypes/proto_p2b_nts.py` 7/7). Code-gate PASS.
- **P2b-2a** (PR #355, validation-only, no C++) — proved the P2b-1 code ALREADY handles a
  DEFORMABLE master (`tests/test_adr39_contact_p2b2.py` 2/2: slave-on-deformable-brick +
  block-on-block, analytic-exact). ⇒ the deferred ∂n/∂u block is a CONVERGENCE REFINEMENT,
  not a correctness gap.
- **P2b-2b `-kn auto`** (#357, merge c49de8773, Zone-A green) — auto-size the penalty per
  (slave,segment) from the owning solid's `getInitialStiff()`: kₙ=f_si·mean(nᵀK_block n),
  f_si=0.10 (element-stiffness form of LS-DYNA 26.14a), GENERIC via base-`Element` virtuals.
  Parser literal `auto`. Oracle `proto_p2b2b_autokn.py` 6/6 + tests 4/4 (incl. ABSOLUTE
  P/pen==oracle + no-owning-solid skip). Code-gate PASS (folded parser stray-token error).
- **P2.5 bucket-sort BROAD PHASE** (#358, merge f891287bf, Zone-A green) — NEW header-only
  `SRC/domain/contact/LadrunoContactBucketSort.h` (§26.11 grid: median-diag cell, runaway
  percentile clip, cap≤min(nSeg,5000), segment-SPAN registration, 27-neighbour candidates w/
  sparse-set de-dup) replaces the brute-force O(nSlave·nSeg) pairing. `contact … -cell <frac>`
  (huge ⇒ 1 bucket = brute force). Oracle `proto_bucket_sort.py` 4/4; gate=verify==brute force
  (`test_adr39_contact_p25_bucketsort.py` bitwise-identical). Code-gate PASS (250k-probe superset
  re-derivation). SCOPE: grid at handle() from reference coords ⇒ ==brute force for static/small-
  motion; large sliding needs epoch re-emit (follow-on).

Test battery: 29/29 (3 P2.5 + 4 P2b-2b + 2 P2b-2a + 8 P2b-1 + 5 P2a + 7 P1). HANDLER_TAG 33002; ELE_TAG 33015
still RESERVED for a future contact element (P2b rides the handler, no new classTag).

## RESUME = P2b-2b: `-kn auto` (the next real feature)

Auto-size the penalty: `kn = f_si·K·A²/V` (LS-DYNA 26.14a, f_si≈0.10) + SOFT Courant floor
`max(½·SOFSCL·m*/Δt_c², …)` (26.15). The **design problem to solve first:**
- K (master stiffness), A (segment area), V (master element volume) have **no Element API**
  (`Element.h` exposes only `getCharacteristicLength`). The master `LadrunoContactSurface`
  stores NODE TAGS ONLY — there is **no link to the parent `LadrunoBrick`**. So step 1 is a
  data-model/parser decision: **auto-detect the owning solid element from the face nodes**,
  OR add a user-supplied master element tag to `contactSurface -master`/`contact`.
- Once linked: K = `masterEle->...materialPointers[gp]->getInitialTangent()(0,0)` (LadrunoBrick
  exposes the GP material array); V = `getCharacteristicLength()`³ or from nodal coords; A from
  the segment node coords. **Cache at setDomain** (not reachable at getResidual — no per-step
  Element access), retain the master element pointer on the surface/contact.
- Parser: `contact … -kn auto` (the value path `-kn $val` stays). Reject auto for the rigid
  plane (P2a) — undefined there (no deformable master).
- Tests: a deformable block where `-kn auto` lands a physically-reasonable kn (penetration
  small vs element size, no lock/blow-up); compare to a hand-tuned kn.

## Then P2b-2c: ∂n/∂u consistent tangent + Hertz (convergence refinement)

The honest NTS tangent `kn·BᵀB + ∂n/∂u` block (drop only O(gn·κ)); the oracle
`proto_p2b_nts.py` already has the FD ground truth (∂n/∂u≈3.6%, symmetric). Gate =
**FD-on-ROTATED** (flat PASS, curved O(gn·κ) shrinks with refinement) + Hertz
`p(r)=p₀√(1−(r/a)²)` refinement benchmark. LOW priority: the residual is exact, the main
term converges (P2b-2a proved it), and the explicit-ship path never forms the tangent.

## Later: P2.5 bucket sort (drop-in broad phase) · P3 IMPL-EX Coulomb friction (SHIP) ·
P3.5 implicit friction · P4 SOFT · P5 segment-based · P6 tied.

## Environment / build / test (verified, don't rediscover)

- Run python (numpy+pytest): `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`.
  Module is `import opensees` (the built `.pyd`), NOT the `openseespy` package.
- DIST (worktree-local): `<worktree>\dist\bin`. Bootstrap: `os.add_dll_directory(DIST);
  sys.path.insert(0,DIST); import opensees`. Tests: `from _testbed import ops`.
- Build (`cmd /c` via PowerShell tool, NOT Bash): `Ladruno_scripts\build.bat OpenSeesPy`.
  **A FRESH worktree is a COLD build — copy `mumps-archive/mumps_src.tar.gz` in first** (from
  compile-root `OpenSees/mumps-archive/`) or build.bat fails on the MUMPS download.
- New authored files → add to `Ladruno_scripts/stamp_headers.py` GLOBS + run it. Ledger every
  new file (LEDGER_implementations) + vanilla edit (LEDGER_vanilla_files) in the SAME PR.
- Branch fresh off `ladruno` per phase (stacked-PR pitfall). Re-check `gh pr view <n> --json
  state` IMMEDIATELY before any follow-up push — this fork auto-merges within minutes.

## Loop discipline (keep it)
Per phase: design → adversarial gate at critical junctions (P2/P3 full multi-agent; mechanical
ports light) → code → build → test → code gate → ledger → commit/PR. Keep `_adr39_loop_state.md`
current every iteration; commit at each phase boundary so the loop survives context resets.
Oracle-first (numpy reference before C++) has caught real issues cheaply every phase.
