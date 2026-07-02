# LadrunoTie shell-to-solid — plane-section rigid-arm coupling of an ndf-6 shell edge to an ndf-3 solid face

> ADR-64. Status: **SHIPPED (PR #478, 2026-07-02).** OQ-1 (operator b-B) and OQ-2
> (frozen-arm / gated-limits v1 contract) signed off; OQ-3 warn-if-omitted, OQ-4..7
> confirmed as recommended. P4.0 oracle `proto_p4_shell_solid_tie.py` 22/22 (T1–T8, both
> honest-limit scalings land exactly on the predicted ν·t and γ·t lines; T7 misfit within
> 0.3% of the analytic γc/(3√3) bound); `.pyd` pre-flight green (rigid arm-engaged +
> constant-moment states enforce under Lagrange via hand-emitted `equationConstraint`
> rows). P4.1+P4.2 `tests/test_ladrunoTie_shellsolid.py` 18/18 (membrane/rigid/section/
> HEADLINE bending/hinge-contrast + 8 refusals + explicit dt_cr-neutral/zero-violation/
> momentum kick–coast); full LadrunoTie battery 50/50 + ADR-30 projection battery 49/49
> no-regression. **One MEASURED correction vs the plan: Decision 1 win 4 — see the
> struck-through note there.** This is the "P4" rung of the LadrunoTie family (ADR-62):
> the one remaining seam where the transfer is a genuinely NEW operator, not a reuse of
> the P row. See `kinematic_tie_handoff.md` (shell-to-solid section) for the probe record.
> Family: ADR-62 (LadrunoTie — the generator this extends, P1/P2/P2.1/P3/P3.1 SHIPPED) ·
> ADR-30 (LadrunoProjectionHandler — the explicit enforcement, reused unchanged) ·
> ADR-41 (mortar D/M) · ADR-39 (ContactDomain projection geometry — deliberately NOT
> reused here, see Where) · ADR-61 (bipenalty — SHELVED, stays shelved). Next free ADR
> slot is 64 (60 = finite-sliding re-emission, 61 = bipenalty-shelved, 62 = LadrunoTie,
> 63 = nodal-normal smoothing). Author: N. Mora-Bowen (Ladruno).

---

## What

Tie an **ndf-6 shell EDGE to an ndf-3 solid FACE**. The solid face carries no rotational
DOF, so the shell rotation must couple to the solid **translation field** — a straight
per-DOF transfer `P` cannot express it. Both LadrunoTie generators currently
**name-REFUSE** exactly this case ("a shell-to-solid tie needs a rotation-from-translation
coupling", `SRC/domain/constraints/LadrunoTie.cpp:531-548` in `generate`,
`:863-876` in `generateMortar`, plus the Hermite corner guard at `:223-229`). **This ADR
replaces that refusal** with a decided operator: the **linearized rigid plane-section map**
(Abaqus shell-to-solid coupling, TG §6.6.2 / §6.6.6), emitted as ordinary mixed-triple
`EQ_Constraint`s through the existing `ltEmitMixedRow` (`LadrunoTie.cpp:154-186`) — the
same shape P3.1's Hermite rows already emit — and enforced by the shipped handlers
(`Lagrange` static, `LadrunoProjection` explicit). No new enforcement code, no new
classTag, no contact-machinery reuse: the master is a short 1D edge polyline, closed-form
point-to-segment projection suffices.

**Feasibility gate — PASSED (record).** Before any design effort, a go/no-go
ndf-coexistence probe was run against the built `.pyd` (scratchpad, recorded in
`kinematic_tie_handoff.md`). All four legs green:

1. **Per-node mixed ndf works** — `node(tag, x, y, z, '-ndf', 6)` creates an ndf-6 node
   inside a default-ndf-3 `model('basic','-ndm',3,'-ndf',3)` (the `BasicBuilder`
   per-node `-ndf` override).
2. **Mixed models assemble + run** — `SSPbrick` (ndf-3) and `ShellMITC4` (ndf-6) coexist,
   assemble, and `analyze(1)` returns 0.
3. **Cross-ndf `EQ_Constraint` rows construct AND enforce through a real solve** under
   `Lagrange`: same-DOF (`shell uz = solid uz`) to ~1e-16, AND a curl-like row
   (`shell θ_y = solid uz-gradient`) reproduced exactly.
4. **`Transformation` goes SINGULAR on cross-ndf rows** (`U(i,i)=0` factorization fail);
   `Lagrange` (static) and `LadrunoProjection` (explicit) enforce cleanly.

⇒ **The remaining difficulty is the operator / work-conjugacy design, not the plumbing.**
That design is decided below.

---

## Why

Shell-to-solid transitions are the standard economy move of 3D structural modeling: a
solid region where the stress state is genuinely triaxial (a joint, a base block, an
anchorage) handing off to shells where plane-section kinematics hold (walls, slabs,
tanks). Every mature code ships this coupling (Abaqus `*COUPLING`-based shell-to-solid,
TG §6.6.6; LS-DYNA `*CONSTRAINED_SHELL_TO_SOLID`). The Ladruno fork now has the entire
tie stack — collocation, mortar, dual basis, ndf-6 rotational rows, Hermite edge transfer
— and this is the last named refusal in the family. Closing it closes the whole seam.

Why NOT the "obvious" curl route the P3 refusal messages themselves point at
(`θ = ½∇×u`): investigated and **rejected** — see Decision 1. The point of this ADR is
precisely that the refusal text guessed the wrong operator; the production answer is the
Abaqus-style plane-section rigid arm, which stays inside the existing surface-only,
emission-only architecture and inherits every guarantee the family already validated.

The moment-transfer question is not academic: a translations-only tie (`-dof 3 1 2 3`,
available today) passes a membrane patch but leaves the shell edge a **hinge** — bending
cannot cross. The headline of this rung is a **constant-MOMENT patch crossing between
dissimilar element types** (oracle gate T3, FE test), which no currently-shippable mode
can pass.

---

## Where

### New code — emission-level only, all in the existing generator files

- **`SRC/domain/constraints/LadrunoTie.{h,cpp}`** — everything lands here. NO gather
  machinery, NO handler change, NO kernel change, NO contact-header change:
  - New flag **`-shellSolid`** (sibling of `-hermite` / `-dual`). **Collocation-only**:
    `-shellSolid -mortar` is refused with a named error (the mortar/weak variant needs
    the plane-section basis inside the M integral — a kernel per-GP hook — deferred with
    mortar-Hermite P3.1b).
  - New command token **`-masterEdge <nseg> a1 b1 a2 b2 ..`** — the master is a **1D edge
    polyline** (segment node pairs), not facets; plus optional **`-thickness <t>`** for
    the arm-length guard.
  - New **`LadrunoTie::generateShellSolid(...)`** — sibling of `generate` /
    `generateMortar`, reusing `generate`'s guard skeleton (domain lookups, per-DOF mass
    scan via `ltCheckTiedDofMass`, named refusals) but with **dofs = {1,2,3}** on the
    slave side (solid translations).
  - New static helper **`ltProjectEdge`** — closed-form point-to-segment projection +
    linear `N_A/N_B`. It does **NOT** touch the facet-only
    `LadrunoContactProjection.h`, and does **NOT** need `LadrunoContactBucketSort.h` —
    edge lists are short; brute-force over segments is the right tool.
  - Emission: **3 translation rows per solid slave node**, each a mixed `nk×6` triple row
    through **`ltEmitMixedRow` verbatim** (`LadrunoTie.cpp:154`), including its 1e-12
    near-zero coefficient filter (`:161`) — which is load-bearing here (Decision 2).
- **The two existing shell-to-solid refusals are KEPT but reworded** — they still guard
  the plain shell-to-shell modes (`generate` `:531-548`, `generateMortar` `:863-876`),
  and now point the user at `-shellSolid` instead of at a nonexistent curl coupling.
- **Oracle** `Ladruno_implementation/kinematic_tie_validation/proto_p4_shell_solid_tie.py`
  (numpy-only, mirrors `proto_p3_1_hermite_tie.py`) + FE tests
  `tests/test_ladrunoTie_shellsolid.py`.

### Modify vanilla — NONE expected

The generator emits standard `EQ_Constraint`s via the public Domain API, exactly like
P1–P3.1. Banner regeneration (`tclMain.cpp` / `PythonModule.cpp` FEATURES blocks via
`patch_banner.py`) is the only vanilla touch → one `LEDGER_vanilla_files` row.

### classTags

**None** — still plain `EQ_Constraint`, still enforced by
`HANDLER_TAG_LadrunoProjectionHandler 33001` (explicit) or `Lagrange` (static).
`classTags.h` unchanged.

---

## How

### The constraint (per solid face node, captured once)

Direction is **inverted from every prior LadrunoTie mode**: the **solid face nodes are
the SLAVES**, the **shell edge nodes (A, B per segment) are the MASTERS**. Each solid
face node `s` projects onto the master edge polyline at parametric `ξ` (closest point,
`ltProjectEdge`) with a **frozen arm** `d = x_s − x̄(ξ)`. The tie is the linearized rigid
plane-section map — three translation rows:

```
u_{s,i} = Σ_{j=A,B} N_j(ξ) · [ u_{j,i} + (θ_j × d)_i ]        i = 1..3
```

i.e. per master node `j`, per slave DOF `i`, coefficients `N_j` on the translation
columns and `N_j·(e_i · (I×d))` on the rotation columns — exactly the `nk×6` mixed-triple
shape `ltEmitMixedRow` already consumes. This IS the through-thickness plane-section
("normals remain straight") kinematics of Abaqus shell-to-solid coupling (TG §6.6.2 /
§6.6.6; consulted via the `abaqus-theory-contact-loading` skill), in its
kinematic/RBE2-flavored form.

### Decision 1 — the operator: rigid plane-section arms (b-B), not curl (a), not gradient rows (b-A)

Three candidates were weighed. Record all three; the losers are load-bearing context.

**(a) Curl synthesis `θ = ½∇×u` — REJECTED.** The route the old refusal text named.
Three independent disqualifiers:
1. It needs the **full 3D displacement gradient**: the in-face `∂/∂n` terms live on the
   solid face, but the `∂/∂a` term (a = the axis running INTO the shell surface) needs
   volume nodes OFF the face ⇒ breaks the surface-only architecture. A volume-node
   neighborhood gather + MLS gradient reconstruction would be new machinery with nothing
   to reuse (the facet-only `LadrunoContactProjection.h` / `LadrunoContactBucketSort.h`
   don't help).
2. The **continuum spin ½∇×u is not the shell section rotation**: they differ by ½γ under
   transverse shear, and that error is baked into EVERY state, not just refinement-limited.
3. The synthesized moment **reacts on interior volume nodes with a mesh-dependent lever
   arm** — the transferred couple changes with the volume mesh. Not a defensible
   production tie.

**(b-A) Shell = slave, through-thickness gradient rows — REJECTED as primary (kept as
documented fallback).** The direction the probe's leg-3 curl-like row demonstrated:
`θ_⊥ = n×(u⁺−u⁻)/t` from a pair of face nodes bracketing the shell mid-surface.
Structural deficiency: the through-thickness node line carries **no information about
drilling** (θ·n). Emitting all 3 rotation rows clamps drilling to zero ⇒ **FAILS the
rigid-body gate**; leaving it free is only expressible as an EQ row when n happens to be
a global axis (rows constrain global DOFs). It also puts the shell's rotational DOFs on
the slave side ⇒ per-DOF BLOCKER-2 demands nodal rotary mass on every tied shell node
(shells' `getMass()` neglects rotary inertia — the P3 gotcha), and it ties only the
bracketing pair, leaving interior face-layer nodes unpinned.

**(b-B) Solid face nodes = SLAVES, shell edge = MASTER, rigid plane-section arms —
RECOMMENDED / DECIDED.** The formula above. Why it wins, point by point:

1. **Drilling resolves itself.** The arm `d` runs ~along the shell normal, so the θ·n
   column of `θ×d` is ~zero and is **dropped by the existing 1e-12 filter in
   `ltEmitMixedRow`** (`LadrunoTie.cpp:161`) ⇒ shell drilling is auto-free (the Abaqus
   behavior), no rigid-body failure, no local-frame gymnastics.
2. **Rigid-body + constant-moment exact by construction** — the rows ARE
   `dx_s = dx_m + dθ×d`; all 6 rigid modes and a constant transferred couple satisfy them
   identically.
3. **Work-conjugate / momentum-clean.** Both handlers apply `f_master = Cᵀλ`. The
   translation columns satisfy partition of unity ⇒ linear momentum transfer exact. The
   θ columns hand the master node the moment `N_j·(d×λ)` ⇒ the transferred couple is
   **statically exact to the frozen-arm (small-rotation) approximation** — the identical
   status `rigidLink -beam` has always had. The projection handler's frozen-lever-arm
   staleness monitor (`flagRotMonitor`, 0.1 rad) auto-arms for these rows for free.
4. **No GENERATOR-level rotary-mass precondition.** The slaves are solid translations —
   always massed via `-rho`, so the P3 per-DOF BLOCKER-2 refusal never fires for this
   mode and static Lagrange needs no rotary mass anywhere. **MEASURED CORRECTION
   (P4.2): the original claim that "the tied solid supplies the shell edge's rotary
   inertia inside the `(LᵀML)` group" is WRONG for the shipped projector** — the
   handler keeps every group DOF in the explicit equation set, so the master edge's
   θx/θy still need a small nodal rotary mass under `LadrunoProjection`
   (`LadrunoConstraintProjector::buildMass` refuses a massless group DOF, and
   correctly so: a zero-mass DOF has no explicit equation of motion). This is generic
   to EVERY rotational tie under handler 33001 (`rigidLink -beam` has the identical
   requirement), not a b-B penalty — b-A would need it too, plus slave-side rotary
   mass. Additionally the shell element's `getMass()` must be LUMPED on the tied
   DOFs: `ShellMITC4` assembles a CONSISTENT (off-diagonal) mass and is refused by
   the handler's element guard; `ASDShellQ4` lumps translational mass and works.
5. **Closes the whole seam** — every solid face node in the through-thickness patch is
   constrained (no unpinned interior layers), one mode, one flag.

**Honest limit — gate it, don't hide it.** The 3-translation rigid arm **suppresses
through-thickness Poisson stretch at the seam**. Abaqus frees it via the SLIDER-style
internal length-change DOF; our EQ rows constrain **global** DOFs and cannot free a skew
direction. This is the kinematic/RBE2 variant of the coupling: a **local St-Venant
boundary layer**, exact at ν=0, error ∝ ν·t·strain, not growing with the model. v1 ships
it suppressed + documented (oracle gate T6 measures it); the free-stretch variant is a
deferred rung, because "constrain a non-global direction" is an architecture change.

### Decision 2 — emission shape: reuse `ltEmitMixedRow` verbatim, dofs = {1,2,3}

Each slave row is `nk` master nodes × 6 DOFs of coefficients — precisely the P3.1 Hermite
signature. `ltEmitMixedRow` (`LadrunoTie.cpp:154-186`) is reused with no change; its
near-zero drop (`:161`) is what silently frees drilling (Decision 1, win 1). Because the
slave DOFs are the solid translations, `ltCheckTiedDofMass` runs on dofs 1..3 only ⇒ the
solid's `-rho` consistent-lumped mass satisfies BLOCKER-2 out of the box — **the P3
rotary-mass requirement DISAPPEARS in this mode.**

### Decision 3 — master geometry: `-masterEdge` polyline, closed-form projection, no broad-phase

The master is a 1D polyline of shell-edge segments (`-masterEdge <nseg> a1 b1 a2 b2 ..`),
not a facet surface, so the ADR-39 bucket-sort + `LadrunoContactProjection::project` are
deliberately NOT reused: `ltProjectEdge` does point-to-segment closest-point in closed
form, brute-force over segments (edge lists are short). Curved/kinked edges are handled
as straight-segment polylines with **per-node frozen arms** in v1 (each slave keeps the
arm off its own closest point — no smoothing rung needed for a tie). Out-of-polyline
projections (ξ outside [0,1] on every segment beyond tol) are refused, named.

### Decision 4 — guards (all named refusals, generator-side)

- Master edge nodes must have **ndf ≥ 6**; slave face nodes **ndf ≥ 3**.
- `-shellSolid -mortar` refused (Decision 6).
- Out-of-bounds edge projection refused (Decision 3).
- **`-thickness <t>` arm guard**: `|d| > (0.5+tol)·t` refused — a slave that far off the
  mid-surface is not this shell's through-thickness material and would get a spurious
  long lever arm. If `-thickness` is omitted: see OQ-3 (warn vs require).
- Coincident-node degeneracy (`|d| ≈ 0`, slave on the edge): emit identity translation
  rows (`u_s = Σ N_j u_j`, θ columns vanish naturally).
- BLOCKER-1 topology (node-disjoint, one row set per slave) inherited from `generate`'s
  skeleton unchanged.

### Decision 5 — handlers: Lagrange (static) + LadrunoProjection (explicit); Transformation documented-unsupported

The probe (leg 4) showed `Transformation` factorization goes singular (`U(i,i)=0`) on
cross-ndf rows. Decision: **no investigation this rung** — record it as a documented
limitation. Static verification uses `Lagrange`; production explicit uses the shipped
`LadrunoProjection`. Test-authoring consequence of the handler's `classifyDerivedSlaves`
mixed-master rule: a test that prescribes shell-edge motion must prescribe/fix **EVERY
referenced master DOF** (sp + homogeneous fix is OK; a free+prescribed mix on one row's
masters refuses).

### Decision 6 — collocation-only in v1; mortar/weak variant deferred

The integral-mortar form of this coupling needs the plane-section basis
`[N_j, N_j(·×d)]` evaluated per Gauss point inside `integratePair`'s M integral — the
same kernel per-GP hook mortar-Hermite (P3.1b) needs. Both defer together; `-shellSolid`
is collocation-only and says so in its refusal of `-mortar`.

---

## Design-gate BLOCKERs

**BLOCKER-1 (inherited) — chain-free, node-disjoint topology.** Solid face nodes are
slaves, shell edge nodes are masters; no node on both lists; each slave gets exactly one
row set off its closest-point segment. The `generate` guard skeleton already enforces
this; `generateShellSolid` keeps it verbatim.

**BLOCKER-2 (inherited, SATISFIED by construction at the GENERATOR) — mass on every tied
DOF.** Slave DOFs are solid translations; `-rho` mass is always there. `ltCheckTiedDofMass`
still runs (defense-in-depth: a genuinely massless solid node is still refused, named).
The P3 generator-level rotary-mass refusal does not apply. **Explicit caveat (measured,
P4.2 — see the Decision 1 win 4 correction): under `LadrunoProjection` the master edge's
θx/θy still need a small nodal rotary mass (every group DOF stays in the explicit
equation set; `rigidLink -beam` precedent), and the shell element must lump its mass
(`ASDShellQ4` yes, `ShellMITC4` consistent → refused).** Static Lagrange needs neither.

**BLOCKER-3 — the frozen arm is a small-rotation operator.** `d` is captured at setup;
under finite interface rotation the transferred couple lags. Same status and same
mitigation as `rigidLink -beam` + the handler's `flagRotMonitor` (0.1 rad) — which
auto-arms on these rows. Gate: the oracle T4 verifies static exactness at the frozen
configuration; the ADR-60 re-emission hook is the (out-of-scope) escape for a tie that
must survive large rotation.

**BLOCKER-4 — the honest kinematic limits must be GATED, not hidden.** Two of them:
(a) through-thickness Poisson stretch suppressed at the seam (error ∝ ν·t·strain, exact
at ν=0) — oracle T6; (b) Timoshenko shear makes the true through-thickness line non-plane
⇒ O(γ·t) local misfit that does NOT shrink with in-plane refinement, while resultant
transfer stays exact — oracle T7. Ship with both measured and documented; if either gate
shows anything other than the predicted scaling, STOP — the operator is wrong, not the
tolerance.

**BLOCKER-5 — cross-ndf rows are Transformation-incompatible.** Probe leg 4. Named in the
docs and the quirks ledger; tests pin `Lagrange`/`LadrunoProjection`. No fix attempted
this rung.

---

## Phased implementation plan (oracle-first, mirroring the fork's method)

- **P4.0 — sign-off + oracle.** Resolve the Open Questions below. Then
  `kinematic_tie_validation/proto_p4_shell_solid_tie.py` (numpy-only, mirrors
  `proto_p3_1_hermite_tie.py`): build the constraint rows exactly as
  `generateShellSolid` will, assemble a reduced model, gates:
  - **T1 — rigid-body, 6 modes** incl. drilling-falls-out-free: verify the θ·n column is
    dropped by the 1e-12 filter (not clamped) and all 6 modes are strain-free.
  - **T2 — constant-stress patch** (membrane, non-conforming through-thickness layout).
  - **T3 — constant-MOMENT patch (the headline)**: moment crossing between dissimilar
    element types — shell side `θ = κx`, `w = ½κx²` ↔ solid side `u_a = −κxz`
    linear-in-thickness, ν = 0. Exact by construction; the gate proves it.
  - **T4 — work-conjugacy / static equivalence**: `Cᵀf` reproduces the exact edge
    resultants; the θ columns contribute zero net force (pure couple); reduced stiffness
    `K_red` symmetric to 1e-14.
  - **T5 — 3D orientation invariance** (rotated global frame, skew edge).
  - **T6 — honest limit, Poisson**: thickness-stretch suppression misfit ∝ ν·t, → 0 at
    ν = 0.
  - **T7 — honest limit, shear**: Timoshenko plane-section misfit O(γ·t), NOT shrinking
    with in-plane refinement; resultant transfer stays exact.
  - **T8 — contrast**: θ columns zeroed (= translations-only tie) passes T2 but FAILS T3
    (hinge) — proves the θ columns are the load-bearing delta of this rung.
  Oracle green ⇒ a **5-minute `.pyd` pre-flight**: hand-emit one mirrored plane-section
  row set via `equationConstraint` on the built binary (the probe's leg-3 pattern, b-B
  direction) to confirm the row direction enforces before any C++ is written.
- **P4.1 — `generateShellSolid` + parser + build + static tests.** The C++ surface of
  Where; `Ladruno_scripts\build.bat`; `tests/test_ladrunoTie_shellsolid.py` under
  `Lagrange`: static membrane patch, **bending-crosses-the-tie** (T3's FE twin),
  rigid-rotation, translations-only contrast (hinge), + the named refusals
  (`-mortar` combo, ndf guards, out-of-bounds projection, `-thickness` arm guard).
- **P4.2 — explicit.** Same models on `CentralDifferenceLadruno` + `LadrunoProjection` +
  `system Diagonal`: momentum-conservation gate + `dt_cr`-neutrality gate (the family
  thesis), and the tied-solid-supplies-rotary-inertia claim observed (group solve
  well-posed with no nodal rotary mass added on the tied group).
- **P4.3 — bookkeeping + ONE PR.** Ledgers + banner + handoff + this ADR flipped to
  SHIPPED; base `ladruno`, watch `gh pr checks`.

---

## Adversarial gate decision

Per the gate-when rule: **full gate on the OPERATOR (done at design time — it is this
ADR), targeted verification on the C++ delta.**

- The novel math is the operator choice, and it received the adversarial treatment
  up-front: three candidates, two rejected with named structural disqualifiers (curl:
  volume-gradient gather + spin≠section-rotation + interior lever arm; gradient rows:
  drilling clamp fails rigid-body + rotary-mass precondition + unpinned layers), the
  winner cross-checked against the Abaqus TG formulation and carrying its two honest
  limits as explicit gates (T6/T7) rather than hopes.
- The C++ delta is then a mechanical sibling of shipped code: `generate`'s guard skeleton
  + `ltEmitMixedRow` verbatim + one closed-form projection helper. That profile matches
  P3's "targeted pass, not a full loop" precedent — one focused adversarial read of
  `generateShellSolid` (index maps, arm sign, filter interaction, refusal reachability)
  after the oracle is green.
- Escalate to a full gate only if P4.0's oracle contradicts a predicted scaling (BLOCKER-4)
  or the pre-flight row direction fails on the real solver.

---

## Ledger / classTag bookkeeping

- `LEDGER_implementations.md` — **extend the existing LadrunoTie row**: shell-to-solid
  `-shellSolid` plane-section coupling, ADR-62 P4 / ADR-64, status per phase, PR.
- `LEDGER_vanilla_files.md` — one banner-regen row (`tclMain.cpp` + `PythonModule.cpp`
  FEATURES blocks via `patch_banner.py`).
- `LEDGER_quirks.md` — sub-points under a new shell-to-solid section:
  (1) cross-ndf EQ rows are **Transformation-incompatible** (`U(i,i)=0`) — use
  `Lagrange`/`LadrunoProjection`; (2) the plane-section tie **suppresses through-thickness
  Poisson stretch** at the seam (St-Venant boundary layer, ∝ ν·t, exact at ν=0) — honest
  limit, not a bug; (3) **explicit rotary-mass rule (CORRECTED from the plan)**: the
  master edge θx/θy need a small nodal rotary mass under the projector (generic
  rotational-tie rule, `rigidLink -beam` precedent) and the shell element must LUMP mass
  (`ASDShellQ4` yes, `ShellMITC4` consistent → refused); the GENERATOR precondition is
  gone and static Lagrange needs nothing; (4) the `mass` command sizes its matrix by the
  BUILDER's ndf, so mixed-ndf models must re-issue `model -ndf 6` before massing shell
  nodes; (5) the collocation force patch is exact only for NESTED interfaces (2:1 etc.) —
  non-nested grids redistribute interface forces at the few-% level (kinematics stay
  exact on any grid; a consistent transfer on non-nested grids is mortar's territory,
  deferred with the `-shellSolid -mortar` variant).
- `banner_features.txt` — extend the LadrunoTie line with `-shellSolid`; rerun
  `patch_banner.py`.
- `kinematic_tie_handoff.md` — flip the shell-to-solid section from "NEXT RUNG" to
  shipped-state on P4.3; cross-link this ADR.
- `classTags.h` — **unchanged** (still plain `EQ_Constraint`). `stamp_headers.py` — not
  needed: no new SOURCE files (`generateShellSolid` lands in the already-stamped
  `LadrunoTie.{h,cpp}`; the oracle/tests are py).

---

## Open questions — need sign-off before coding

- **OQ-1 (the operator — HEADLINE).** Confirm (b) through-thickness coupling over
  (a) curl synthesis, and within (b) the **b-B direction** (solid nodes = slaves, shell
  edge = master, rigid plane-section arms) over b-A (shell = slave, gradient rows;
  retained as documented fallback only). I recommend b-B — Decision 1 carries the case.
- **OQ-2 (work-conjugacy scope — HEADLINE).** Accept as v1: "statically-exact couple,
  frozen small-rotation arm, thickness-stretch suppressed (RBE2-style)", with the ν
  (T6) and γ (T7) limits **GATED and documented, not engineered away**. The free-stretch
  (Abaqus SLIDER) variant is a deferred rung — it needs a "constrain a non-global
  direction" capability, an architecture change.
- **OQ-3 (`-thickness`).** Warn-if-omitted (skip the arm guard) or require it? I lean
  warn: the guard is a footgun-catcher, not a correctness precondition.
- **OQ-4 (drilling).** Confirm drilling-left-free (the θ·n column dropped by the filter)
  is the intended contract — matching Abaqus — rather than adding an optional drilling
  clamp.
- **OQ-5 (mortar/weak variant).** Confirm deferred (with mortar-Hermite P3.1b, shared
  kernel per-GP hook); `-shellSolid -mortar` refuses in v1.
- **OQ-6 (Transformation).** Confirm documented-unsupported for cross-ndf rows, no
  investigation this rung.
- **OQ-7 (curved/kinked edges).** Confirm v1 = straight-segment polyline with per-node
  frozen arms (no smoothing rung); a tie doesn't slide, so ADR-63-style smoothing is not
  pulled in.
