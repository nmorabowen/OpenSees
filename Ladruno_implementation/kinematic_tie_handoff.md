---
title: "LadrunoTie (ADR-62) — handoff: P1 collocation + P2 integral-mortar SHIPPED; P2.1 + P3 deferred"
project: Ladruno
type: handoff
status: P0 (numpy) + P1 (collocation, PR #454) + P2 (integral-mortar, PR #455) MERGED to ladruno. Deferred = P2.1 dual-basis sparsification + P3 shell/rotational ties.
related:
  - "[[62_ladruno_kinematic_mesh_tie_adr]]"            # the spec (P1/P2 marked SHIPPED)
  - "[[30_ladruno_explicit_constraint_projection_adr]]" # the enforcement handler (SHIPPED, reused unchanged)
  - "[[41_ladruno_mortar_alm_contact_adr]]"            # the mortar D/M kernel reused by P2
  - "[[39_ladruno_contact_domain_adr]]"                # bucket-sort + projection geometry
  - "[[projection_handler_handoff]]"                   # handler API + limits
tags: [adr, constraints, explicit-dynamics, mesh-tie, projection, kinematic, mortar, handoff]
---

# LadrunoTie — handoff (P1 + P2 shipped; P2.1 / P3 next)

## TL;DR — what exists now (all on `ladruno`)

A non-conforming explicit **mesh-tie**, enforced KINEMATICALLY (not by penalty): the generator
emits ordinary `EQ_Constraint`s and the **shipped** `LadrunoProjectionHandler` (ADR-30) enforces
them — exact, `dt_cr`-neutral, momentum-clean, no fictitious mass. ONE command, two modes:

```
# P1 — node-to-segment COLLOCATION (each slave node -> one master facet)
LadrunoTie -slaveNodes <ns> s1..  -masterFacets <npsM> <nf> m..  [-dof nd d..] [-tol frac]

# P2 — integral MORTAR (slave SURFACE -> master surface, weak form)
LadrunoTie -mortar -slaveFacets <npsS> <nf> s..  -masterFacets <npsM> <nf> m..
           [-dof nd d..] [-tol frac] [-outward ox oy oz]
```

Both register in openseespy and the modern Tcl interpreter (NOT the classic `OpenSees.exe` shell —
like `equationConstraint`, they live in the `DL_Interpreter`/`TclWrapper` path). **No new class tag**
(emits `EQ_Constraint`, enforced by handler 33001). **No vanilla source edits** beyond the banner.

- Code: `SRC/domain/constraints/LadrunoTie.{h,cpp}` — `LadrunoTie::generate` (P1) + `::generateMortar` (P2) + `OPS_LadrunoTie`.
- Oracles: `kinematic_tie_validation/proto_p{0,1}_kinematic_tie*.py` (P1) + `proto_p2_mortar_tie.py` (P2, 13/13).
- Tests: `tests/test_ladrunoTie_patch.py` (8, P1) + `tests/test_ladrunoTie_mortar.py` (7, P2). Full suite 15/15.
- PRs: ADR/oracles #449, P1 #454, P2 #455 — all merged.

## The architecture (so the next agent doesn't re-derive it)

A mesh-tie is **static** (a permanent, non-sliding bond), so geometry is computed ONCE at model-build
and frozen into constraints. Both modes condense to a per-slave-node transfer `u_s = P u_m` and emit
one `EQ_Constraint` per slave node **per translational DOF** with the P row as coefficients (Ccr = P
row, NO negation — the ModelBuilder `EQ_Constraint` ctor stores the coef vector directly):

- **P1 collocation:** `P_sk = N_k(ξ̄_s)` — the master facet shape functions at the slave node's
  closest-point projection. Reuses `LadrunoContactProjection` (project + shape) + `LadrunoContactBucketSort`.
- **P2 mortar:** assemble global `D_IJ=∫N_I^s N_J^s dΓ`, `M_IK=∫N_I^s φ_K^m dΓ` over the clipped
  overlap (reusing `LadrunoMortarKernel::integratePair` VERBATIM), then `P = D⁻¹M` via one
  `Matrix::Solve` (DGESV). **KEY INSIGHT: pre-inverting D globally makes every P row tie to MASTER
  nodes only ⇒ no MP-chains ⇒ the shipped handler accepts it (dense, master-only rows) with NO handler
  change and NO kernel change.** The ADR's original "P2 needs deferred MP-chain support" premise was
  WRONG — global D⁻¹ is what dodges it.

Both inherit the handler's requirements (and the generator front-loads them as named refusals):
`system Diagonal` + explicit (for the projection path; Lagrange/Penalty work for static), every tied
slave DOF needs LUMPED mass, node-disjoint chain-free surfaces, ICs on the manifold (conforming-at-
interface geometry — NO snapping).

## GOTCHAS already paid for (don't rediscover)

1. **Explicit projection needs LUMPED mass on tied DOFs.** `stdBrick` assembles CONSISTENT mass →
   the handler's R5 guard REFUSES it on a tied DOF. Use **`SSPbrick`** for explicit tie tests/models.
   (`stdBrick` is fine for the static Lagrange patch tests — no mass there.)
2. **BLOCKER-2 mass check scans ELEMENT mass.** At model-build the nodal `mass()` is still 0 (solids
   get mass from element `-rho`, assembled later), so the generator scans `Element::getMass()` and
   marks a node massed if any incident element has a nonzero mass diagonal. ⇒ **Define `LadrunoTie`
   AFTER the elements/masses**, or it false-refuses.
3. **`dt_cr` kernel reads ELEMENT mass** (use `-rho`), not nodal `mass`. And `criticalTimeStep()`
   needs `integrator CentralDifferenceLadruno -cfl` + one safe `analyze(1, dt0)` before the query.
4. **Mortar coverage guard must use the per-node FULL tributary ratio, not the interface average.**
   `cover[I]` (overlap integral) is compared to `fullCov[I]` (= ∫N_I over the WHOLE slave facet,
   computed by a SELF-CLIP `integratePair(npsS,Xs,npsS,Xs,...)`). A single slave facet half-overlapping
   the master gives NO node `cover≈0` (its shape fn spans the overlap) — only the `cover/fullCov` RATIO
   catches a slave protruding past the master. DGESV flags only an EXACT-zero pivot, so a post-solve
   partition-of-unity check (`|Σ P_Ik − 1| < 1e-6`) backstops a near-singular D.
5. **Standard-basis P is DENSE** (every slave couples to every master through `D⁻¹`) ⇒ ONE large
   handler group (`(LᵀML)` over all master DOFs). Correct + fine for typical interfaces; the cost is
   the P2.1 motivation. **Row-sum LUMPING D is NOT a shortcut** — it keeps partition-of-unity but
   BREAKS linear completeness ⇒ fails the constant-stress patch.
6. **Build/run env:** set `LADRUNO_OPENSEES_QUIET=1` when building (else the splash banner pollutes
   CMake's FindPython version string on a reconfigure → "Invalid character escape '\U'" configure
   fail). To run tests against THIS worktree's `.pyd`, set `PMI_SIZE=1` (makes the site `.pth` boot
   module skip its pre-import of another worktree's opensees) + `sys.path.insert(0, dist/bin)` +
   `os.add_dll_directory(dist/bin)`.

## Build + run recipe (this machine)

- Python 3.12: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe` (numpy + pytest).
- Build (PowerShell tool, NOT Bash): `set LADRUNO_OPENSEES_QUIET=1 && Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP OpenSeesPy`. ~30 min cold; incremental after. ZLIB-not-found ⇒ CMake-4.3 shadow, use conan cmake by full path ([[project_build_env_cmake43_conan_zlib]]). The `.pyd` lands in `dist/bin`.
- Run a test:
  ```
  $env:LADRUNO_OPENSEES_QUIET="1"; $env:PMI_SIZE="1"
  python -c "import os,sys; D=r'...\dist\bin'; os.add_dll_directory(D); sys.path.insert(0,D); import pytest; sys.exit(pytest.main(['tests/test_ladrunoTie_mortar.py','-v']))"
  ```
- Stamp headers after adding source files: `python Ladruno_scripts/stamp_headers.py` (LadrunoTie.* already in GLOBS).

## Deferred backlog (pick the next session's target)

### P2.1 — dual / biorthogonal basis (sparsify the mortar transfer) — OPTIMIZATION
The standard-basis condensation gives a DENSE `P` (one big interface group; `(LᵀML)` over all master
DOFs). A **dual (biorthogonal) Lagrange-multiplier basis** (Wohlmuth 2000) makes `D` DIAGONAL ⇒
`P = D⁻¹M` is LOCAL (each slave ties only to masters under its own support) ⇒ sparse rows, small
local handler groups, cheap projection. Requires constructing the dual shape functions on each slave
facet (a per-facet linear transform of the standard basis) — a NEW kernel path
(`LadrunoMortarKernel` currently emits STANDARD-basis `D`/`M`). Preserves linear completeness (that's
the whole point of dual vs lumped). Value: only at LARGE interfaces (100s+ of interface nodes); the
dense path is correct and fine below that. Oracle-first: extend `proto_p2_mortar_tie.py` with a
`dual_D` (diagonal) variant and show same patch test + sparse P. Low correctness risk, medium effort.

### P3 — shell / rotational ties (ndf 6) — CAPABILITY (likely higher user value)
Today the generator ties **translational DOFs only** (default `-dof 1..ndm`). A shell-to-shell or
shell-to-solid tie must also tie the **rotational DOFs** (4,5,6 for an ndf-6 shell node). Open issues
to resolve (worth an oracle + a short design note, NOT a full ADR):
- **Rotational transfer weights:** does the rotation field interpolate with the same `N_i`/`P` as the
  translation field (so `θ_s = Σ P_sk θ_{m,k}`), or does a shell-to-solid tie need the solid's
  translation field to DRIVE the shell rotation (a `rigidLink -beam`-style `θ = ½∇×u` coupling)?
  Shell-to-shell is the clean first case (same ndf, same P on rotations).
- **Rotational lumped mass on tied DOFs:** the projection handler needs nonzero mass on every tied
  DOF; shell rotational mass is small but present (mirror the `rigidLink -beam` / 3D-diaphragm rule the
  handler already documents). BLOCKER-2's element-mass scan must pick up the shell's rotational mass
  diagonal — verify `ShellMITC4`/`ASDShellQ4` `getMass()` has nonzero rotational entries, else refuse.
- **Drilling DOF** (6th, in-plane rotation) is often near-zero-stiffness — decide whether to tie it
  or leave it free (likely leave free / tie only the two bending rotations, configurable via `-dof`).
- Scope first cut: shell-to-shell, same ndf, `θ_s = Σ P_sk θ_{m,k}` (reuse the existing P), gated by a
  rotational-mass check. Test: a constant-curvature (bending) patch must cross the tie exactly.

### Also noted (out of scope for a tie, but adjacent)
Finite-sliding re-emission (the ADR-60 hook) is for a tie that must survive LARGE interface rotation;
a tie doesn't slide, so this is only relevant if the frozen-`Ccr` small-rotation limit (~0.1 rad) is
exceeded. Not a planned increment.

## Bookkeeping state (all current on `ladruno`)
- `LEDGER_implementations.md` — LadrunoTie row says "shipped (P1 + P2 mortar)", PRs #454/#455.
- `LEDGER_vanilla_files.md` — interpreter-wiring + banner-regen rows for #454/#455.
- `LEDGER_quirks.md` — two entries: the handler-requirements refusal contract, and the mortar
  global-D⁻¹/dense-P/self-clip-coverage gotcha.
- `banner_features.txt` — "LadrunoTie — kinematic ... (collocation + integral-mortar; ...)".
- `classTags.h` — unchanged (no tag). Memory: [[project_adr61_bipenalty_shelved]] carries the full
  LadrunoTie record (P1 + P2).
