---
title: "ADR 214 — Native per-mode RSA output: responseSpectrumPerMode (WP-216, program link P1)"
project: Ladruno
type: ADR
status: "PR #1032 open (WP-216)"
priority: medium
owner: nmora
related:
  - "[[213_ladruno_modal_basis_completeness_adr]] (mode sets, LadrunoModeKinds, D3b Sa-period override)"
  - "WP-207 (`responseSpectrumAnalysis -increment`), WP-212 (mass guard)"
  - "Program issue #1024 (dynamic analysis), chain M link P1"
  - "apeWorkbench opensees-dynamic-analysis, branch work/claude/rsa-gap-study: studies/rsa-gaps/apegmsh-adoption-guide.md §1.1, §5"
tags: [adr, response-spectrum, performance, reactions, wp-216]
updated: 2026-10-27
---

# ADR 214 — Native per-mode RSA output

## 1. Context

- **Evidence from a real tower.** A shell-wall tower (ASDShellQ4 walls,
  ElasticTimoshenko frames, LadrunoMembrane floors, lateral lumped mass,
  rigid diaphragm, ~52k nodes) runs RSA as a Python loop over modes:
  `responseSpectrumAnalysis -mode k -increment`, `reactions` + `nodeReaction`,
  pier `eleForce`, CQC in Python. With 60 modes and 2 directions RSA took
  232 s of a 9-min run: 102 s in `responseSpectrumAnalysis` (0.84 s/call),
  81 s in `reactions` (0.67 s/mode), 29 s in pier `eleForce`. Per-mode
  readback, not mass completeness, is the main RSA cost.
- **Where the time goes (read from the source, then measured, §4).** One
  `-increment` call is two full domain sweeps: `revertDomainToLastCommit`
  reverts every node and element, re-applies the loads and runs
  `Domain::update()` (every element's state determination), then the mode
  is imposed on every node and `updateDomain()` runs `Domain::update()`
  again. `reactions` is a third sweep: `Domain::calculateNodalReactions`
  re-applies the loads and calls `getResistingForce()` on every element. The
  quantities wanted per mode need a few hundred elements of tens of
  thousands: those attached to a support, and the pier groups.

## 2. Decisions

### D1 — a separate Ladruno command, not an RSA option

```tcl
responseSpectrumPerMode (-ts $tag | -Tn $T.. -Sa $Sa.. | -fn $f.. -Sa $Sa..)
    -dir $d1 <$d2 ..> <-scale $s> <-modes $m1 $m2 ..>
    <-reactions <-about $x $y <$z>>>
    <-reactionSet $name $node.. <-about ..>> ...
    <-eleGroup $name -ele $ele.. <-cut $node..> <-about ..>> ...
    <-inertia $name $node.. <-about ..>> ...
```

- `responseSpectrumAnalysis` has one positional direction and leaves state
  behind (a commit, or an `-increment` trial); the per-mode output wants
  several directions in one call and no state. A separate command keeps the
  stock RSA syntax untouched.
- The spectrum options mean what they mean in `responseSpectrumAnalysis`
  (`-ts` is its `$tsTag`; `-Tn`/`-fn`/`-Sa` lists and `-scale` are the same).
  The command builds a `ResponseSpectrumAnalysis` and calls its `check()`
  (modalProperties of this mode set, else refused), its WP-212 `checkMass()`
  once per requested direction (no free or no participating mass is refused,
  as `-mode k -increment` refuses it; the cumulative-mass report and the 90 %
  warning print the same) and `getSa()` (same interpolation, `-scale` applied
  once, WP-202). These three members are made public; that is the whole
  vanilla touch.
- Every number list (`-Tn`, `-fn`, `-Sa`, `-dir`, `-modes`, node and element
  lists) is read either expanded (`{*}$list`, Python `*list`) or, in classic
  Tcl, as one unexpanded list argument (`-Tn $Tn`), as
  `responseSpectrumAnalysis` reads `-Tn`/`-Sa` (review of #1032).
- `-modes` defaults to every mode of the current set; order is kept.
- `-about` applies to the group defined just before it and adds r × F about
  the point (undeformed coordinates) to the moment components. Without it the
  moments are the plain sum of nodal moments (what summing `nodeReaction` /
  `eleForce` gives).
- Wired once in `SRC/interpreter/LadrunoCommandTable.h` (Python, DL Tcl,
  classic Tcl bridge).

### D2 — the quantity is the increment over the committed state

Every output is the value at trial = committed + u_k,d minus the value at
the committed state. With an unloaded, undeformed committed state (the usual
RSA) it is exactly what the loop reads. After a gravity step + `loadConst`
it is the modal part alone (the loop then has to subtract the gravity
reactions itself; the gate does). The domain is left at its committed state
(one `revertDomainToLastCommit` at the start, the touched objects reverted at
the end), so a call never leaves an `-increment` trial behind.

### D3 — same amplitude as `-mode k -increment`, override per ADR 213 D3b

u_k,d = φ_k s_k Γ_k,d Sa(T*) / λ_k, written with the same expressions in the
same order as `ResponseSpectrumAnalysis::solveMode` (T = 1/(√λ/2π);
T* = the `LadrunoModeKinds` Sa-period override when set, else T). So the
trial displacements are bit-identical to the loop's. The command never
combines modes; it returns, per mode, what a combination needs:

| key | per | meaning |
|---|---|---|
| `period` | mode | 2π/√λ (Rayleigh period for residual / Ritz modes) |
| `saPeriod` | mode | the period Sa was read at (T*, D3b) |
| `omega` | mode | the CQC ω: 2π/T* for an overridden mode, √λ otherwise |
| `term` | mode | combination term, 1-based: overridden modes sharing one T* share one term (D3b: sum their responses algebraically into one rigid term before CQC/SRSS) |
| `kind` | mode | 0 eigen, 1 residual, 2 Ritz (`LadrunoModeKinds`) |
| `Sa`, `gamma` | dir × mode | Sa(T*) with -scale; Γ_k,d |

### D4 — update only what a quantity needs, in the `-increment` sequence

- Elements touched: those with a node carrying an SP (when reactions are
  requested; found by one pass over the connectivity, no state
  determination) and the `-eleGroup` elements. Nodes touched: their nodes.
- Per (direction, mode), on the touched objects only: nodes and elements
  `revertToLastCommit`, elements `update()` at the committed state, nodes
  `setTrialDisp(committed + u)` (solveMode), elements `update()` and
  `getResistingForce()`. That is the loop's sequence restricted to the
  touched objects. The first `update()` is not optional: ASDShellQ4 condenses
  its internal DOFs from the previous trial, so updating straight from a
  revert lands 1e-12 away from the loop on cancelling sums (measured on the
  gate model: 2.6e-12 on a core-wall reaction set).
- No load re-application: element loads (Q in `getResistingForce`) are the
  same at both states and cancel in the increment.
- Serial. The touched set is small; the threaded `Domain::update()` relies on
  a per-element-type audit that a hand-picked subset would have to repeat.

### D5 — reactions on SP-constrained DOFs only

A reaction is summed only on a DOF that carries an SP. At a free DOF of a
partially fixed node (a 2D floor node fixed in UY, a diaphragm master fixed
out of plane) `nodeReaction` returns the element force sum there, which for
a mode shape is not zero: it is the mode's out-of-balance (inertia) force,
λ M φ. Adding it to a base reaction double-counts inertia. At a fully fixed
node this is exactly `nodeReaction`. `-reactions` sums every SP node;
`-reactionSet` nodes must carry an SP.

### D6 — inertia resultants from the nodal mass

`-inertia` sums f = M_n φ_k s_k Γ_k,d Sa(T*) (= λ_k M_n u_k,d) over the
named nodes, M_n = `Node::getMass()`. No element work. Element mass matrices
(ρ on shells, consistent-mass frames) are not included: on such models use
reactions or element groups for storey forces. For an eigenmode the total
base reaction equals minus the total inertia (ιᵀKφ = 0 ⇒ R = −λ ιᵀMφ); the
gate checks it to 1e-6. For residual and Ritz vectors Kψ ≠ λMψ, so the
balance holds only after summation over the complete set, not per vector.

### D7 — output: a dict of flat lists, [dir][mode][component]

`OPS_SetDoubleDictListOutput`: a Python dict / a Tcl dict.

- `shape` = [ndir, nmode, ncomp], ncomp = 6 in 3D (Fx Fy Fz Mx My Mz),
  3 in 2D (Fx Fy Mz) — the ndf of `modalProperties`.
- `dirs`, `modes` (1-based, as given), the D3 table keys.
- `reactions`, `reaction:<name>`, `ele:<name>`, `inertia:<name>`: flat
  ndir·nmode·ncomp, row-major [dir][mode][component];
  `numpy.reshape(shape)` recovers the array.
- `numSupports`, `numUpdated` (how many elements a mode updates).
- Signs: reactions and element-group resultants follow `nodeReaction` and
  `eleForce` (the element's resisting force at its nodes).

### D8 — refusals

Refused before anything is touched: no AnalysisModel; OpenSeesSP; no
`modalProperties` or stale modal properties (`check()`); a direction with no
free or no participating mass (`checkMass()`, per direction); no spectrum, both
`-ts` and lists, unequal or non-increasing lists, negative Sa, unknown time
series; no `-dir`, a direction out of 1..ndf or repeated; nothing requested;
a mode out of range, repeated, or with λ ≤ 0 (the loop would divide by it);
an unknown node or element; a `-reactionSet` node without an SP; an empty
group; an element listed twice; a `-cut` node on no element of its group;
a group name used twice; `-about` with no group before it or with the wrong
number of coordinates; an unknown option.

### D9 — a TransformationFE overflow fixed here (#1028)

The gate model crashed (0xc0000374 at the second or third `wipe`) under
Transformation only. `TransformationFE::dofData` holds one DOF count per
node but was sized `MAX_NUM_DOF` = 64, so a LadrunoMembrane floor with ~100
nodes wrote past it on every assembly. It is fixed in this PR, because the
gate cannot run without it: the buffer and its WP-192 thread-scratch copy
are grown to the largest element node count. Any Transformation model with
a membrane of more than 64 nodes was exposed, including the tower that
motivated P1.

### Rejected

- **Linearity shortcut** (one update per mode, scale by Γ_d Sa_d for the
  other directions): halves the element work but is exact only for linear
  elements; a `-pDelta` preload or corotational element would make the
  output depend on the shortcut. Not worth a second code path at the
  measured cost (§4).
- **`calculateNodalReactions` per mode**: the whole-domain sweep is the cost
  being removed.
- **Combining inside the command**: CQC/SRSS per quantity is cheap in the
  caller once the per-mode arrays exist, and keeps signed / directional
  combinations (adoption guide A5) out of the fork.

## 3. Verification (`tests/test_wp216_rsa_per_mode.py`)

Oracle: the stock per-mode loop on the same model (`-mode k -increment`,
`reactions` + `nodeReaction` masked to the SP DOFs, `eleForce`,
λ_k M · IncrDisp), every quantity minus its committed value. Model: the
synthetic tower of `tests/_testbed/rsa_tower.py` at 3 storeys (1001 nodes,
membrane floors of ~100 nodes), 6 modes.

**Tolerance: 1e-12 of the term scale, per mode and component.** The error
of each (direction, mode) is the largest, over the components, of
|new − loop| divided by that component's Σ|term|: forces against forces,
moments against moments (review of #1032; it first divided by the largest
scale of any component, which let forces through ~10–20× looser). The terms
are the ones each sum is formed from: for reactions, the element force blocks
at the support nodes (what `calculateNodalReactions` adds; several shells meet
at a base node and largely cancel, so the nodal reactions themselves
understate the floor: a weak-direction base shear of 1.7 against nodal terms
of 160 missed 1e-12 by 1.8×); for element groups, the blocks at the cut; for
inertia, the nodal forces. A component whose terms are all exactly zero must
match exactly. Relative to the result itself a 1e-12 gate is below the loop's own
reproducibility: ASDShellQ4's AGQI update starts from the `KQU`/`Q_residual`
of the element's last `calculateAll` call, which `revertToLastCommit` does
not restore, so the stock loop read with and without its `eleForce` calls,
or in reverse mode order, differs at round-off (ledger quirk WP-216). On a
cancelling sum (the vertical reactions of the core under overturning) that
is ~1e-12 of the result. Measured on the residual set, new − loop and
loop − loop(reversed) agree quantity by quantity, e.g. 3.2e-12 against
3.4e-12 of the term scale on the pier at storey 2 in the 0.75 ms pseudo-mode
(which follows modes of much larger amplitude). The residual-set check
therefore holds the command to max(1e-12, 4 × the loop's own spread),
the spread being the larger of two history changes (reverse mode order;
no `eleForce` reads). One sample of round-off is a noisy estimate of its
size: 2 × a single sample failed once by 15 %.

| Check | Result |
|---|---|
| (a) 6 modes × 2 dirs, Transformation / Penalty / Lagrange; total, core set, 3 pier groups, 2 inertia groups, `-about` | ≤ 1e-12 of term scale (per node ~3e-16) |
| (a) `-modes 5 2`, `-ts` Path series, `-scale 0.7` | ≤ 1e-12 |
| (a) gravity preload + `loadConst` (oracle subtracts the gravity state) | ≤ 1e-9 |
| eigenmodes: base reaction = − total inertia (force, moment about a point) | < 1e-6 |
| (b) `residualModes -dir 1 2 -period 0.05` set (10 modes) | ≤ max(1e-12, 4 × loop spread); spread < 1e-10; `kind`, `saPeriod` = 0.05, one shared `term`, `omega` = 2π/0.05 |
| (b) `ritzVectors 6 -dir 1 2` set | ≤ 1e-12; `kind` = 2 |
| (c) 23 bad calls refused; the next good call is bit-identical; stale modal properties refused; λ ≤ 0 refused | pass |
| the loop itself: with vs without `eleForce` reads | not bit-equal; < 1e-12 of term scale |
| a call after an `-increment` trial: domain left committed, two calls bit-identical | pass |
| TransformationFE, 104-node membrane (#1028) | Transformation = Plain (1e-10), 4 rounds, no crash |
| classic Tcl bridge = Python (2D), unexpanded `-Tn $Tn -Sa $Sa`, node list `$top` | 1e-12 |
| a direction without mass (`-dir 2`, alone or with `-dir 1`): refused in Python and Tcl, as stock `-mode k -increment` | pass |

## 4. Profile and gain

`tests/bench_wp216_rsa_per_mode.py` on the synthetic tower
(`tests/_testbed/rsa_tower.py`; no project data). The machine is the nmb
desktop: Windows, MSVC Release build, serial element loop, `system Pardiso`.
Both directions are run. The old loop is the one the tower project ran:
per (direction, mode) `-mode k -increment`, `reactions`, `nodeReaction` on
the supports, `eleForce` on the piers of storeys 1 and nst/2, and the
inertia of two storey groups in Python (one `nodeDisp` per node, as there;
the right read is IncrDisp, which costs more, so the old timing is a lower
bound). The new call requests the same quantities in one command.

| tower | nodes / elements / supports | modes × dirs | old, ms per (mode, dir): total (rsa / reactions / reads / inertia) | new, ms per (mode, dir) | elements updated per mode | gain |
|---|---|---|---|---|---|---|
| 10 storeys, h = 0.4 m | 18 529 / 17 830 / 239 | 25 × 2 | **334.5** (135.8 / 149.5 / 1.9 / 47.3) | **10.1** | 332 | 33× |
| 45 storeys, h = 0.5 m | 50 190 / 48 555 / 195 | 60 × 2 | **921.5** (385.0 / 409.9 / 1.6 / 124.9) | **18.3** | 268 | 50× |

- **Profile → target.** In the old loop, `responseSpectrumAnalysis` and
  `reactions` take 86 % of the time, and both grow with the element count:
  two whole-domain updates per RSA call, one whole-domain
  `getResistingForce` sweep per `reactions`. The readback calls themselves
  (`nodeReaction`, `eleForce`) cost under 1 %. So the target was to take every
  whole-domain sweep off the per-mode path, leaving a cost proportional to
  the touched elements and the group nodes. Making the sweeps faster
  (threads) would have kept them and attacked the wrong term.
- **New-call breakdown (10 storeys, per (mode, dir)).** Reactions only:
  4.6 ms. Reactions + piers: 6.2 ms. Inertia only (two groups, ~33k node
  visits, no element work): 3.8 ms. Each call also pays one
  `revertDomainToLastCommit` (a whole-domain update) and one connectivity
  pass, spread over its pairs. What is left grows with the group sizes, not
  with the model.
- **Tower scale.** The tower that motivated P1 spent 232 s on 120 pairs at
  ~52k nodes. At the 45-storey rate (18 ms per pair) the same readback is
  about 2 s.

## 5. Consequences

- apeGmsh (A5) can read per-mode reactions, pier forces and storey inertia
  in one call per spectrum and combine in numpy.
- `nodeDisp` after `-increment` is the committed displacement (ledger quirk
  WP-216): per-mode Python loops that form inertia from `nodeDisp` read zero.
- Element mass is not in `-inertia` (D6).

## 6. Follow-ups

- Side items found on the same tower, filed separately: #1029 (UmfPack out
  of memory reported as a bare status −1), #1030 (`eigen -genBandArpack` with
  Pardiso factorizes twice), #1031 (keep FEAST off default paths).
- Re-run the tower project's RSA with the #1028 fix (it uses membrane floors
  under Transformation).
