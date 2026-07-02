---
title: "LadrunoTie (ADR-62) — handoff: P1 collocation + P2 integral-mortar + P2.1 dual + P3 shell/rotational + P3.1 Hermite SHIPPED; shell-to-solid deferred"
project: Ladruno
type: handoff
status: P0 (numpy) + P1 (collocation, PR #454) + P2 (integral-mortar, PR #455) + P2.1 (dual/sparse mortar, PR #462) + P3 (shell/rotational ndf-6, PR #459) + P2.1×P3 composition test (PR #464) + P3.1 (Hermite w–θ edge transfer, -hermite) MERGED to ladruno. Deferred = mortar-Hermite (P3.1b, needs a kernel per-GP hook) + shell-to-solid ties.
related:
  - "[[62_ladruno_kinematic_mesh_tie_adr]]"            # the spec (P1/P2 marked SHIPPED)
  - "[[30_ladruno_explicit_constraint_projection_adr]]" # the enforcement handler (SHIPPED, reused unchanged)
  - "[[41_ladruno_mortar_alm_contact_adr]]"            # the mortar D/M kernel reused by P2
  - "[[39_ladruno_contact_domain_adr]]"                # bucket-sort + projection geometry
  - "[[projection_handler_handoff]]"                   # handler API + limits
tags: [adr, constraints, explicit-dynamics, mesh-tie, projection, kinematic, mortar, handoff]
---

# LadrunoTie — handoff (P1 + P2 + P2.1 + P3 shipped; P3.1 / shell-to-solid next)

## TL;DR — what exists now (all on `ladruno`)

A non-conforming explicit **mesh-tie**, enforced KINEMATICALLY (not by penalty): the generator
emits ordinary `EQ_Constraint`s and the **shipped** `LadrunoProjectionHandler` (ADR-30) enforces
them — exact, `dt_cr`-neutral, momentum-clean, no fictitious mass. Ties translational DOFs (solids)
AND the rotational DOFs of ndf-6 shells (**P3**: default `-dof` → 1..6 for a 3D ndf-6 node,
`θ_s = Σ P_sk θ_{m,k}` with the same P; `-dof`-selectable). ONE command, two modes:

```
# P1 — node-to-segment COLLOCATION (each slave node -> one master facet)
# -hermite = P3.1 rotation-consistent cubic Hermite w–θ transfer (ndf-6 shell EDGE ties)
LadrunoTie -slaveNodes <ns> s1..  -masterFacets <npsM> <nf> m..  [-dof nd d..] [-tol frac] [-hermite]

# P2 — integral MORTAR (slave SURFACE -> master surface, weak form); -dual = P2.1 sparse basis
LadrunoTie -mortar -slaveFacets <npsS> <nf> s..  -masterFacets <npsM> <nf> m..
           [-dof nd d..] [-tol frac] [-outward ox oy oz] [-dual]
```

Both register in openseespy and the modern Tcl interpreter (NOT the classic `OpenSees.exe` shell —
like `equationConstraint`, they live in the `DL_Interpreter`/`TclWrapper` path). **No new class tag**
(emits `EQ_Constraint`, enforced by handler 33001). **No vanilla source edits** beyond the banner.

- Code: `SRC/domain/constraints/LadrunoTie.{h,cpp}` — `LadrunoTie::generate` (P1) + `::generateMortar` (P2; `dual` param = P2.1) + `OPS_LadrunoTie`; P3 adds `ltDefaultDofs`/`ltScanMassedDOFs`/`ltCheckTiedDofMass` helpers (per-DOF mass) + a shell-to-solid master-DOF guard in both generators.
- Oracles: `kinematic_tie_validation/proto_p{0,1}_kinematic_tie*.py` (P1) + `proto_p2_mortar_tie.py` (P2, 13/13) + `proto_p2_1_dual_mortar.py` (P2.1, 12/12) + `proto_p3_rotational_tie.py` (P3, 8/8) + `proto_p3_1_hermite_tie.py` (P3.1, 13/13).
- Tests: `tests/test_ladrunoTie_patch.py` (8, P1) + `tests/test_ladrunoTie_mortar.py` (11 = 7 P2 + 4 P2.1 dual) + `tests/test_ladrunoTie_shell.py` (13 = 6 P3 + 1 P2.1×P3 composition + 6 P3.1 hermite). Full suite 32/32.
- PRs: ADR/oracles #449, P1 #454, P2 #455, P3 shell/rotational #459, P2.1 dual #462, P2.1×P3 composition test #464, P3.1 hermite #467 — all merged.

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

## Shipped increments + deferred backlog

### P2.1 — dual / biorthogonal basis (`-mortar -dual`) — SHIPPED (PR #462)
Sparsifies the mortar transfer: a **dual (biorthogonal) basis** (Wohlmuth 2000) makes the slave mass
`D` DIAGONAL ⇒ `P = D_dual⁻¹M` is LOCAL (each slave ties only to masters under its own facet support)
⇒ sparse rows, small local handler groups. **KEY IMPLEMENTATION INSIGHT (avoided the feared "new kernel
path"):** the dual basis ψ=Aᵉ·N is a per-slave-FACET linear transform, `Aᵉ=diag(∫N)·(Dᵉ)⁻¹`, so
`Dᵉ_dual = Aᵉ Dᵉ = diag(∫N)` EXACTLY — built from the SAME `integratePair` Dᵉ/Mᵉ via a per-facet
npsS×npsS solve at the `LadrunoTie` level ⇒ NO kernel change, NO handler change (the SAME dodge as P2's
global-D⁻¹). `P(I,k)=Y(I,k)=(Dᵉ⁻¹Mᵉ)(I,k)`; PoU holds exactly (Dᵉ·1=Mᵉ·1 over the same overlap). The
naive alternative — row-sum LUMPING D — keeps PoU but BREAKS linear completeness (fails the patch; oracle
T6 = 0.28 err). Default OFF = the standard dense P byte-identical (dual is a separate post-guard block).
Oracle `proto_p2_1_dual_mortar.py` (12/12: D diagonal, sparse, dual≠lumped, dual≠dense) + 4 FE tests
(dual patch/split/explicit + `-dual`-requires-`-mortar`). The RUNTIME win is in the handler (sparse P ⇒
small groups); the setup does a 2nd integration pass (opt-in, one-time). Value: LARGE interfaces (100s+
interface nodes) — the dense path is correct and fine below that.
**Composition with P3 verified:** `-mortar -dual` on an ndf-6 SHELL tie
(`test_mortar_dual_rigid_rotation_crosses_tie` in `test_ladrunoTie_shell.py`) crosses a rigid rotation
exactly on all 6 DOFs — as expected, since `-dual` only changes how P is computed while the
rotational-DOF defaulting / per-DOF mass guard / shell-to-solid guard live in the shared emission path.

### P3 — shell / rotational ties (ndf 6) — SHIPPED (PR #459)
The generator now ties the **rotational DOFs** (4,5,6) of an ndf-6 shell node with the SAME per-slave
transfer `P` (`θ_s = Σ P_sk θ_{m,k}`), in both collocation and mortar modes. Resolutions of the open
issues that were listed here:
- **Rotational transfer weights:** shell-to-shell reuses the SAME `N_i` / `D⁻¹M` row as the
  translations (no new P math). Shell-to-**solid** (which WOULD need the `θ=½∇×u` curl coupling) is
  OUT of v1 scope — guarded with a named master-DOF-existence refusal in both generators, deferred.
- **Rotational lumped mass:** **surprise — `ShellMITC4` AND `ASDShellQ4` `getMass()` NEGLECT rotary
  inertia** (zero on DOFs 4,5,6; verified in source — see `LEDGER_quirks.md`). So BLOCKER-2 was made
  **PER-DOF**: it names-refuses a tied rotational DOF with no mass and tells the user to add nodal
  rotary mass (`mass $n 0 0 0 mrx mry mrz`) or tie translations only (`-dof 3 1 2 3`). The default
  `-dof` becomes 1..6 for a 3D ndf-6 node (`ltDefaultDofs`).
- **Drilling DOF (OQ-P3, user-decided):** tie ALL 3 rotations (4,5,6) by default, `-dof`-selectable.
  Rationale: the drilling rotation is the shell LOCAL-normal rotation, a mix of global 4/5/6 unless the
  shell is in a global plane — so "skip DOF 6" is only correct for a global-Z-normal shell; the mass
  guard already refuses a massless drilling DOF cleanly.
- **Honest limit (P3.1):** `P` is linearly complete ⇒ rotations + CYLINDRICAL bending with the
  curvature axis ∥ the interface cross EXACTLY (the FE patch test in `test_ladrunoTie_shell.py` uses
  the stock `ShellMITC4`). A curvature varying ALONG the interface leaves an O(h²) residual on the
  quadratic `w` (never on θ) — the Hermite w–θ transfer is the deferred P3.1 rung.
- Oracle `proto_p3_rotational_tie.py` (8/8) + `tests/test_ladrunoTie_shell.py` (6/6). Handler
  UNCHANGED (rotational EQ rows are just more master-only union-find vertices).

### P3.1 — Hermite (rotation-consistent) w–θ shell transfer (`-hermite`) — SHIPPED
Makes a curvature varying ALONG the interface cross exactly (the P3 honest limit): with `-hermite`
(collocation mode, ndf-6 shell EDGE/butt-joint ties) each slave's shell-NORMAL translation and
interface-TANGENT slope are reconstructed with the CUBIC HERMITE basis over the master edge nodes'
`(w, θ_t)` pairs; in-plane translations + non-slope rotations keep the linear weights:
```
u_s = (I−n⊗n)·Σ N_k u_k + n·[Σ Hw_k (n·u_k) + Hp_k (e·θ_k)]
θ_s = (I−e⊗e)·Σ N_k θ_k + e·[Σ dHw_k (n·u_k) + dHp_k (e·θ_k)]      e = t×n, slope = θ·e
```
(from `u = θ×r ⇒ ∇w = n×θ ⇒ dw/ds = θ·(t×n)`; the transfer is invariant under n→−n). Exact through
CUBIC w (O(h⁴) beyond — strict dominance over the linear P's O(h²)); every P3-exact case still crosses.
**OQ answered (user signed off collocation-only): it stays a LadrunoTie-LEVEL transform** —
`EQ_Constraint` natively stores per-retained `(node,dof,coef)` triples (EQ_Constraint.h:37-45) and
`LadrunoProjectionHandler::buildGroups` builds one union-find vertex per (node,dof) with NO
retained-dof==constrained-dof assumption ⇒ mixed `w←(w,θ)` rows are ordinary EQ rows, NO kernel/handler
change. **Mortar-Hermite is the exception** (`integratePair` returns only ACCUMULATED D/M — no per-GP
hook — so the Hermite master basis can't enter the M integral without a kernel extension) ⇒ deferred
(P3.1b); `-hermite -mortar` is a named refusal. Other named refusals: full 6-DOF tie required (the rows
couple u↔θ), ndf-6 on BOTH sides, slave must project onto a master facet EDGE (or corner → identity
rows); interior projections refused (surface-overlap collocation keeps the standard tie). HONEST
LIMIT: Kirchhoff on the tie line (slope=rotation) — Mindlin/shear-deformable data carries a bounded
O(γ·h) tie error (γ = shear angle), vanishing thin + on refinement. Oracle
`proto_p3_1_hermite_tie.py` (13/13: rigid modes, quadratic/cubic exact, O(h⁴)-vs-O(h²) rates,
3D-inclined global-DOF rows, γ-bound; gotcha: the shear-error test must NOT sample midpoints —
`H2+H4 ∝ ξ−3ξ²+2ξ³` has a root at ξ=½). FE tests: headline along-interface bending crosses exactly
at the non-conforming mid-edge node + linear-tie contrast (κ/4 vs exact κ/8) + aligned-bending
no-regression + 3 refusals. **Shell-to-solid** (the `θ=½∇×u` coupling) remains the open sibling rung.

### Also noted (out of scope for a tie, but adjacent)
Finite-sliding re-emission (the ADR-60 hook) is for a tie that must survive LARGE interface rotation;
a tie doesn't slide, so this is only relevant if the frozen-`Ccr` small-rotation limit (~0.1 rad) is
exceeded. Not a planned increment.

## Bookkeeping state (all current on `ladruno`)
- `LEDGER_implementations.md` — LadrunoTie row says "shipped (P1 + P2 mortar + P2.1 dual + P3 shell/rotational)", PRs #454/#455/#459/#462.
- `LEDGER_vanilla_files.md` — interpreter-wiring + banner-regen rows for #454/#455/#459/#462.
- `LEDGER_quirks.md` — four entries: the handler-requirements refusal contract, the mortar
  global-D⁻¹/dense-P/self-clip-coverage gotcha, the shell `getMass()`-neglects-rotary-inertia gotcha (P3),
  and the dual-basis-vs-lumping sparsifier note (P2.1).
- `banner_features.txt` — "LadrunoTie — kinematic ... (collocation + integral-mortar, -dual sparse biorthogonal basis; solid + ndf-6 shell rotational DOFs; ...)".
- `classTags.h` — unchanged (no tag). Memory: [[project_adr61_bipenalty_shelved]] carries the full
  LadrunoTie record (P1 + P2 + P2.1 + P3).
