---
title: "apeGmsh adoption guide — the 2026-09-07 fork batch (TIMs proposed-model asks)"
project: Ladruno
type: guide
audience: apeGmsh team (bridge `apeSees`, Ladruno recorder/reader, Substep controller, profiler)
status: current as of ladruno `a240b9183` (2026-09-07)
owner: nmora
related:
  - "[[96_ladruno_contact_passenger_dof_adr]]"
  - "[[ladruno_apegmsh_contract]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[ladruno_solver_flag_guide]]"
  - "[[LadrunoKinematicCoupling_guide]]"
  - "[[ndf_and_mixed_models_guide]]"
  - "[[zerolength_and_link_springs_guide]]"
tags: [apegmsh, adoption, tims, contact, zerolength, sanisand, pardiso, recorder]
---

# apeGmsh adoption guide — the 2026-09-07 fork batch

Eight fork PRs merged on 2026-09-07 in answer to the TIMs proposed-model request
(`_tims_proposed_model_requests_2026-09-07.md`): #805, #808, #810, #811, #812, #814, #820,
#821. This guide lists, per feature, **what changed in the fork, what apeGmsh must do to use
it, and how to verify the adoption**. Minimum fork build: `ops.ladrunoBuild()` must be a
descendant of `a240b9183`. Everything below is serial; nothing changes the MP contract.

Items are ordered by what they unblock on the TIMs ladder (rung 1 spring, rung 2 contact,
post-processing, profiling), then the contract clarifications.

---

## 1. Interfaces on u-p nodes: `zeroLength` and contact accept ndf >= 3 (ADR 96, #808)

**Unblocks:** TIMs rung 1 (node-to-node no-tension spring) and rung 2 (contact) on the
coupled model; apeGmsh slice **A10** (`interface()` in 3D was blocked on ZeroLength's ndf
guard) is unblocked.

**What changed.** The pressure DOF (or any DOF past the third) rides as a *passenger*: never
read, never written, never coupled. This is **not** ADR 47 deferral 9 — no gap flow, no
pressure penetration.

- **`zeroLength` in 3D** now accepts any node pair with **both ndf >= 3** outside the vanilla
  `(3,3)` / `(6,6)` table: `(3,4)`, `(4,3)`, `(4,4)`, `(3,6)`, `(6,4)`. The element keeps its
  6-slot translational core and scatters it into an `ndf1 + ndf2` element; each node's
  first three DOFs are the translations.
  - Only translational `-dir 1 2 3` exist on such a pair. A rotational `-dir` (4–6) is
    **refused** with `"... (passenger mode, ADR-96): only translational -dir 1..3 exist there ...
    element disabled"`; the element stays in the domain and contributes nothing.
  - The vanilla refusal for any other mismatch (`(2,3)` etc.) keeps its wording
    `"... have differing dof at ends for ZeroLength ..."` but **no longer crashes** at the
    `element` command (vanilla returned half-initialised and the post-add `update()`
    dereferenced a NULL; measured). If apeGmsh modelled that crash as a guard, the guard is
    now a warning plus an inert element.
  - **Reader impact:** the element `force` response is **element-sized** (`ndf1 + ndf2`
    slots), node 2's translations at `ndf1 .. ndf1+2`, every passenger slot identically
    `0.0`. A reader that assumes a 6-vector on every zeroLength must size it by the two
    nodes' ndf. `deformation`, `material` and `basicForce` responses stay the 6-slot core.
- **Contact (NTS, mortar, edge-edge, rigid plane), 3-D lane:** every node on either side may
  have ndf >= 3 (u-p ndf 4, beam/shell ndf 6); the adapter's ID map takes each node's first
  three equations. The handler's six `ndf != 3` guards are now `ndf < 3`. **The 2-D lane is
  unchanged**: `ndf == ndm == 2` exactly (its apeGmsh adoption record still holds).

**What apeGmsh does.**
- Lift the A10 block: `interface()` in 3D may emit `element zeroLength <tag> <skinNode>
  <soilNode> -mat <ENT> -dir 3` (or `-dir 1 2 3` with three materials) between an ndf-3 skin
  node and an ndf-4 `LadrunoUP` node. Always pass `-dir` explicitly; never a rotational one
  on a mixed pair.
- Contact surfaces may list u-p nodes as slaves or masters unchanged; nothing to add to the
  emitted `contactSurface` / `contact` lines.
- If apeGmsh has an ndf-equality pre-check on zeroLength pairs, relax it to
  `ndm == 3 and min(ndf) >= 3` (plus the vanilla pairs) and keep refusing 2-D mismatches.

**Verify on your side** (the fork's G2/G3 gates, `tests/test_adr96_passenger_dof.py`):
- Spring: ENT zeroLength across a `(3,4)` pair, DOF 4 of the ndf-4 node imposed by `sp` at a
  non-zero value → same displacement as a `(3,3)` pair (never read) and reaction exactly `0`
  on DOF 4 (never written). Use `constraints Transformation`: `Plain` silently homogenises a
  non-zero `sp`.
- Contact: a `LadrunoUP` column under an ndf-3 platen in NTS, static drained → summed normal
  traction = applied load (fork measured 6.6e-12 rel), all slaves in compression, pore
  pressure identical to the `equalDOF 1 2 3` twin (fork measured 2.6e-24).

Docs: `96_ladruno_contact_passenger_dof_adr.md`, `ndf_and_mixed_models_guide.md` (rows for
ZeroLength and contact), `zerolength_and_link_springs_guide.md` §2.1.

---

## 2. Material responses `psi` and `yieldDistance` on `LadrunoSANISAND` (#805) and named components on every fork response (#820) — apeGmsh slice A12

**Unblocks:** the mechanism block's per-Gauss-point state without a post-processing step.

**What changed.**
- Two new read-only scalar responses (committed state):
  - `psi` (alias `stateParameter`) = `e - e_c(p')`, `p' = p + p_residual` floored at `1e-10`
    — the model's own `GetPSI`, the psi behind `M^b`/`M^d`. Default `p_r = 0` means plain
    `e - e_c(p)` from `state[24]` and the mean stress.
  - `yieldDistance` (alias `yieldFunction`) = `|s - p' alpha| - sqrt(2/3) m p'` on the
    committed pair: negative inside the cone, `~mTolF` (1e-7 abs, default) on it, never
    positive after a converged return. Fork measured `|f| / (sqrt(2/3) m p') ~ 1e-5` on the
    surface at ~2.5 kPa.
- All seven fork responses now emit `ResponseType` component names (XML header and
  `.ladruno` HDF5 `COMP_NAMES`), with the strings agreed with the apeGmsh session:

| response (token) | slots | component names |
|---|---|---|
| `psi` | 1 | `psi` |
| `yieldDistance` | 1 | `yieldDistance` |
| `implexError` | 1 | `implexError` |
| `avgImplexError` | 1 | `avgImplexError` (process-wide mean; every GP reports the same) |
| `substeps` | 2 | `substeps_me`, `substeps_capHit` |
| `implexDetail` | 6 | `implexDetail_total`, `_dev`, `_vol`, `_clampFired`, `_clampCount`, `_f` |
| `implexRefusals` | 4 | `implexRefusals_total`, `_signChange`, `_control`, `_companion` (process-wide ledger) |

  Vanilla names are unchanged and still unnamed (C1..Cn): `stress` 6, `strain` 6, `state` 26
  (void ratio `[24]`, dGamma `[25]`), `alpha` 6, `fabric` 6, `alpha_in` 6, `estrains` 6,
  `plasticstrains` 6. Vanilla `ManzariDafalias` answers none of the fork names (empty).
  Plane-strain wrapper: `stress` is the 3-vector; `psi` / `yieldDistance` are still computed
  from the full internal 6-vectors.

**What apeGmsh does (A12).** Register the tokens above in the Gauss-level canonicaliser
(`_KIND_TO_ROOT` / `_CONTINUUM_SCALAR_TOKENS` / `RESPONSE_CATALOG`) as per-GP scalars and
named small vectors; with #820 the names arrive as real `COMP_NAMES`, so the `C1..Cn` fallback
is no longer the path for them. The Ladruno recorder passes `-E` tokens verbatim, so the
emitter needs nothing new: `recorder ladruno ... -E psi yieldDistance implexRefusals`.

**Verify:** `ops.eleResponse(ele, 'material', gp, 'psi')` returns one float; a recorded
`.ladruno` file lists `psi` under `COMP_NAMES`; on the fork's confine-first cube
`yieldDistance` is `-sqrt(2/3) m p'` exactly before yield (alpha = 0).

Docs: `LadrunoSANISAND_implex_guide.md` §6, §6.1.

---

## 3. `LadrunoKinematicCoupling` refuses an ambiguous slave ndf without `-dof` (#814)

**apeGmsh already guards this (A1, #1100).** The fork now refuses it for every caller: a
slave whose ndf is neither `ndm` nor `ndm + nrot` (an ndf-4 u-p node in 3D; an ndf-3 node in
2D is fine) with the **default** component list is refused at the parser:

```
WARNING LadrunoKinematicCoupling <tag>: slave node <n> has ndf = 4, which is neither 3
(translations) nor 6 (translations + rotations) in 3D; the default component list would tie
node DOF 4 (a pressure or other passenger DOF on a u-p node) to a master rotation. Pass -dof
explicitly (e.g. -dof 1 2 3) -- REFUSED
```

**What apeGmsh does.** Keep emitting `-dof 1 2 3` (or the explicit list) for any u-p slave;
the emitted command for ndf-3 / ndf-6 slaves is unchanged. If the bridge surfaces fork
warnings, match on `REFUSED` and `-dof`.

Docs: `LadrunoKinematicCoupling_guide.md` §1.2 and the option table.

---

## 4. Desktop factorisation statistics: `system Pardiso -stats` (#821), and why not MUMPS (#810)

**Unblocks:** PM-01 D26 (memory and fill in every desktop leg's log) and the apeGmsh profiler
slice **A8** (`s.profile(...)`).

**What changed.** `system Pardiso -stats` (alias `-pardisoStats`) existed since ADR-75 P1d but
printed once per sparsity pattern, in MB, with Mflops always 0. It now prints the MUMPS-shaped
block **after every numeric factorisation** (phase 22, including refactorisations), through
`opserr`:

```
PARDISO stats: n=<n> nnz(A)=<nnz> matrixType=<mtype> threads=<nthreads>
  factor entries iparm(18)  = <nnz in L+U>
  peak memory KB iparm(15)  = <peak during symbolic>
  perm memory KB iparm(16)  = <permanent>
  fact memory KB iparm(17)  = <numerical factorization + solve>
  factor Mflops  iparm(19)  = <mflops>
```

Labels are the 1-based MKL names. The nnz and Mflops sentinels (`iparm[17]`, `iparm[18]` set
to -1) are armed only with `-stats`, so runs without the flag are byte-identical to before.
**The serial `MumpsSolver` is never compiled in this fork** (`_MUMPS` is defined only for the
parallel targets): `system Mumps` on `OpenSees.exe` / the desktop pyd answers *unknown system
type*. MUMPS statistics exist only on `OpenSeesMP` / `openseesmp` rank 0 (`system Mumps -stats`).

**What apeGmsh does.** In the profiler, emit `system Pardiso -stats` (with the deck's
`-matrixType` if used), capture stderr, and parse the block: the first line's `n=`, `nnz(A)=`,
`matrixType=`, `threads=`, then five `label = value` lines. Expect one block per factorisation
(so per step under `Newton` with a fresh tangent), and store the per-stage max of the three
memory lines. Do not emit `system Mumps` on a desktop leg.

**Verify:** `tests/test_pardiso_stats.py` in the fork asserts the exact labels; a 54-DOF brick
gives `factor entries iparm(18) = 1836`.

Docs: `ladruno_solver_flag_guide.md` §`-stats` (PARDISO subsection).

---

## 5. Contract clarifications that affect apeGmsh code paths (no fork feature)

- **Sealed static u-p system is SILENT, not loud** (ADR-71 §3.2 reworded, #814). Every serial
  general solver factorises the structurally singular all-impervious static system through
  round-off and returns `rc = 0` with an arbitrary pressure level. The guard is apeGmsh's build
  gate **G4** (A2, #1102: a fixed pressure DOF in every hydraulically connected region) — keep
  it mandatory for static u-p decks; do not expect the fork to refuse.
- **`Results.from_ladruno` is shipped** on the apeGmsh side (`apeGmsh/results/Results.py:491`);
  the fork's `ladruno_apegmsh_contract.md` now says so (#814).
- **Stage switching is process-wide.** `mElastFlag` is a `static` on `ManzariDafalias`: one
  `updateMaterialStage` (command or `parameter`) through ANY SANISAND tag flips every SANISAND
  instance in the process. `s.update_parameter` users cannot hold one layer elastic while
  another is plastic; stage all SANISAND materials together (PM-01 §17.4).
- **`setParameter` contract on `LadrunoSANISAND`** (for A4 `s.update_parameter`): the three fork
  names `implexError`, `avgImplexError` (read-only) and `implexDt` (writable, >= 0) match on
  `argv[0]` only — an appended material tag is ignored, harmlessly. Every base name
  (`updateMaterialStage`, `materialState`, `poissonRatio` — the K0 trick —, `refShearModulus`,
  `voidRatio`, `IntegrationScheme`, `Jacobian`, ...) needs `argv[1] == materialTag`. Constructor
  options (`-Presidual`, `-Pmin`, `-implex*`, `-maxSubsteps`, `-honorTolR`, `-tanType`) are not
  parameters. With `-implex` on, the IMPL-EX history initialises **at** the stage flip: set
  `implexDt` after it, not before.
- **Signals for the A7 Substep controller** (no parameter needed, all through responses):
  `yieldDistance` per GP as the "on the surface" flag (`>= -tol * sqrt(2/3) m p'`); the per-step
  delta of `implexRefusals[0]` — a refused step returns `analyze` rc `= -33086`
  (`LADRUNO_MATERIAL_REFUSED`, propagated by exact value) and is the retry-with-smaller-step
  signal; `substeps[1]` cap-hit and `implexDetail[3]` clamp-fired as quality warnings.
- **Rigid footing driver is geometrically linear** (F6 scoped, #812): `LadrunoKinematicCoupling`
  builds its gap operator once from the reference lever arms; the reference-point moment
  transfer TIMs post-process is exact under that linearisation only. TIMs' G4 measures whether
  a corotational transport is needed; nothing to adopt now.
- **ADR 93 candidate II.1** (elastic-only floor `p_r,e`, #811) is the fallback of record, **not
  built**. If it is built, `LadrunoSANISAND` gains a `-Pelastic <p>` constructor flag (default
  0, byte-identical); the emitter would pass it through like `-Presidual`. Nothing to do now.

---

## 6. Checklist for the adoption PR(s) on the apeGmsh side

1. Bridge: require `ops.ladrunoBuild()` descendant of `a240b9183` for the features above (the
   bridge's existing fail-loud build check).
2. A10: `interface()` 3D emits `zeroLength ... -dir ...` across ndf-3/ndf-4 pairs; pre-check
   relaxed; reader sizes `force` by `ndf1 + ndf2`.
3. A12: canonicaliser registers the seven fork tokens with the component names in §2.
4. Profiler (A8): emits `system Pardiso -stats`, parses the §4 block, never `system Mumps`
   on desktop.
5. Emitter: `-dof` always explicit on u-p slaves (already A1); nothing else changes.
6. Docs: point `internal_docs/contact_2d_adoption.md`'s 3-D sibling (or the interface doc) at
   ADR 96 and record that the 2-D lane is unchanged.

Fork-side tests to mirror: `tests/test_adr96_passenger_dof.py`,
`tests/test_ladruno_sanisand_responses.py`, `tests/test_ladruno_sanisand_responsetype.py`,
`tests/test_pardiso_stats.py`, `tests/test_ladrunoKinematicCoupling_element.py` (the three
`up_slave` tests).
