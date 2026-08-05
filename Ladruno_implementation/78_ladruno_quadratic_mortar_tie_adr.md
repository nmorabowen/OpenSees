# LadrunoTie quadratic mortar — quad8/tri6 (serendipity) facets for the integral-mortar mesh-tie

> ADR-78. Status: **IMPLEMENTED + ADVERSARIALLY GATED (2026-08-04) — Q0–Q3
> shipped same day; Q4 (apeGmsh cross-check) CLOSED 2026-08-04, apeGmsh PR #898,
> case-A agreement 5.0e-13 — see Q4 for the contract revision it forced.**
> Oracle `proto_p2_2_quad8_mortar.py` 38/38; FE
> `tests/test_ladrunoTie_mortar_quad8.py` 13/13; all shipped tie suites (40) +
> contact suites (91) green unchanged. Gate: 2 independent refute-oriented
> reviews + 2 killed mutations; 3 MAJORs found and fixed (master self-overlap
> false-accept, negative-mass false-pass, unverified quadratic-master basis) —
> see "Adversarial gate — RUN". Signed off same day. OQ-1 refuse
> tri6 slaves in both bases, OQ-3 degree-6 rule for quadratic pairs only, OQ-4
> new ADR-78, all as recommended; **OQ-2 resolved as "unify"** — the sign-free
> guards (per-facet area coverage, per-facet L1 gap, total-area zero-overlap,
> sign-aware dual mass) replace the rowsum guards for ALL facet orders, one code
> path (see the revised D3/BLOCKER-3). Extends ADR-62 P2/P2.1 (the
> `-mortar` / `-dual` mesh-tie, SHIPPED) and ADR-41 C1 (`LadrunoMortarKernel`).
> Trigger: apeGmsh ADR 0086 D2 resolved as **"both"** — the apeGmsh numpy mortar
> kernel ships there, and extending THIS kernel to quadratic facets is the named
> fork-side follow-up so the two kernels can cross-check long-term
> (apeGmsh repo, `src/apeGmsh/opensees/architecture/decisions/0086-mortar-tie-composed.md`).
> Family: ADR-62 (kinematic mesh-tie — the command + generator this extends) ·
> ADR-41 (mortar kernel) · ADR-30 (LadrunoProjectionHandler — enforcement, untouched) ·
> ADR-35 (HRZ lumped mass — the hex20 mass precondition) · ADR-72 (second-order
> bricks — where the quad8 faces come from). Next free ADR slot is 78 (77 =
> implicit transient). Author: N. Mora-Bowen (Ladruno).

---

## What

Extend the integral-mortar mesh-tie path — `LadrunoTie -mortar [-dual]`, the
`LadrunoMortarKernel` beneath it, and the `ltReadFacets` parser — from tri3/quad4
facets to **quadratic serendipity facets**:

* **quad8** accepted as **slave and master** (the driving case: a hex20 element
  face is quad8);
* **tri6** accepted as **master only**; a tri6 **slave** surface is a **named
  refusal** (the dual-basis degeneracy, D2 below);
* tri3/quad4 behaviour stays **byte-identical** (D3).

The collocation mode, `-shellSolid`, `-hermite`, and the whole contact runtime
(`LadrunoContactFE` / `LadrunoContactHandler` mortar contact) stay **linear-only**
— out of scope here, with named pointers where a user could plausibly try.

## Why

ADR-62's mortar tie exists to couple non-matching meshes, but its real payoff is
coupling meshes of **mixed element order** — hex20 ribs on hex8 plates — and that
is exactly the case it cannot express today: `ltReadFacets` rejects `nps ∉ {3,4}`
(`LadrunoTie.cpp:1191`), `generateMortar` re-checks (`:592-595`), and the kernel
arrays are fixed at `X[4][3]` / `D[4][4]` / `M[4][4]`
(`LadrunoMortarKernel.h:83-90`). apeGmsh ADR 0086 measured the cost of the
workarounds: collocating hex20 onto hex8 over-constrains the rib root
(a = 0.195 vs 0.230 conformal, K +7.7 %); going all-hex20-fine costs +53 % solve
time. apeGmsh is shipping a quad8-capable numpy port of this kernel; the fork
side must catch up so (a) native decks can express the same tie and (b) the two
independent implementations cross-check each other on the same patch geometry —
the whole point of the "both" resolution.

## Where

* `SRC/domain/contact/LadrunoMortarKernel.h` — the generalization (header-only,
  OpenSees-free, oracle-testable — the ADR-41 discipline is preserved).
* `SRC/domain/constraints/LadrunoTie.cpp` — `generateMortar` (guards + dual
  condensation + local arrays), `ltReadFacets`, the collocation `generate()`
  error message.
* `Ladruno_implementation/kinematic_tie_validation/proto_p2_2_quad8_mortar.py` —
  the new numpy oracle (mirrors `proto_p2_mortar_tie.py` / `proto_p2_1_dual_mortar.py`).
* `tests/test_ladrunoTie_mortar.py` — FE additions (hex20-on-hex8 patch, refusals).
* **Modify vanilla: NONE.** **classTags: none** (same `EQ_Constraint` emission,
  same handler).

Blast radius outside the tie: `PairResult` grows (`D[8][8]`, `M[8][8]`, `g[8]`)
⇒ `LadrunoContactFE.cpp:349` recompiles against the bigger struct but its code
and behaviour are unchanged (it still calls with nps 3/4, and array-of-array
parameters decay to `(*)[3]`, so no signature break). The contact handler's
quadratic-facet skip (`LadrunoContactHandler.cpp:559-561`) stays as-is: mortar
**contact** at runtime needs the linearization stub filled and is not a tie
concern.

## How — the mechanics

### D1 — geometry from the corner polygon; quadratic shape functions at the mapped Gauss points

The serendipity node ordering is corners-first (gmsh / OpenSees 20-node-brick
face convention): quad8 = corners 0-3 then midsides 4-7 (4 between 0-1, … 7
between 3-0); tri6 = corners 0-2 then midsides 3-5. Let `ncorner = nps ≤ 4 ? nps
: nps/2`.

All **geometry** — `auxPlane`, `facetNormal`, the Sutherland-Hodgman clip,
`inverseIsomap2D`, `projectFull` — runs on the **corner sub-facet** (first
`ncorner` nodes), exactly as today. Only the **basis evaluation** changes: at
each mapped Gauss point `(ξ,η)` the kernel evaluates the full serendipity basis
(new `shapeFull(nps, ξ, η, N[8])`; quad8 corners ¼(1+ξξᵢ)(1+ηηᵢ)(ξξᵢ+ηηᵢ−1),
mids ½(1−ξ²)(1+ηηᵢ) / ½(1−η²)(1+ξξᵢ); tri6 corners λ(2λ−1), mids 4λₐλᵦ) and
accumulates `D`, `M`, `g̃` over all `nps` rows/columns. `LadrunoContactProjection`
is **untouched** — `projectFull` is called with the corner facet, and the kernel
evaluates the quadratic master `φ` itself at the returned `(ξ,η)`.

Why this is exact, not an approximation, in scope: when the midside nodes sit at
the edge midpoints, the serendipity isoparametric map **reduces identically to
the corner map** (the quadratic terms cancel — review-verified numerically to
~1e-15 **including warped corner sets**), so the corner-based `(ξ,η)` IS the
quadratic facet's parametrization. The **guard** checks midside-at-midpoint only
(`|X_mid − (X_a+X_b)/2| ≤ tol·edge`, a 3D distance, so a midside lifted off its
straight edge is caught); corner coplanarity is deliberately NOT checked — a
warped quad8 degrades exactly like the shipped warped quad4 (the documented
~0.7 % constant-J area bias), no new wrongness. A violating facet is a **named
refusal, not a silent skip** —
`integratePair` today returns −1 for both "empty overlap" (fine, skip) and
"degenerate" (silent); quadratic-geometry violations get a distinct **status −2**
which `generateMortar` converts to a hard error naming the facet and the remedy
(this tie is the ADR-0086 v1 scope: flat, straight-edged, coincident interfaces —
curved quadratic interfaces are future work, refused loudly until then).

### D2 — tri6 slave surfaces: refused, named

On the reference triangle the tri6 **corner** shape functions integrate to
**exactly zero** (∫N_corner = 0; the midsides carry the whole area A: 3 × A/3).
The dual condensation scales by the per-facet row-sum `cᵉ_a = ∫ₑ N_a`
(`Aᵉ = diag(cᵉ)(Dᵉ)⁻¹`) — for tri6 corners that is a division by zero, and no
tolerance rescues it: it is structural. The rowsum-based coverage machinery
(`cover[]`, the near-zero-row backstop) is equally meaningless at tri6 corners.
quad8 corners integrate to **−A/12 ≠ 0** and are fine (D3 handles the sign).

Decision: refuse a tri6 **slave** surface in both bases (not just `-dual`) with
a named error stating the degeneracy and the remedies — swap master/slave (tri6
master is fine), or use quad8/hex20 on the slave side. Allowing tri6 slaves in
the standard basis only would fork the guard machinery for a case with no
driving user (tet10 interfaces are deferred fork-wide); revisit when one exists.
This mirrors apeGmsh ADR 0086 v1 exactly.

### D3 — sign-free guards, unified across all facet orders (OQ-2: "unify")

quad8 corner row-sums are **negative** (−A/12), which breaks three shipped
guards that implicitly assume `∫N_I ≥ 0`. Per the OQ-2 sign-off the sign-free
replacements apply to **every** facet order — one code path, and the L1 gap
check is strictly stronger than the signed nodal average it replaces (no
cancellation blind spot). This is a **behaviour change for shipped linear
decks**: marginal geometry that previously slipped the per-node guards (e.g. an
antisymmetric gap averaging to ~0) now refuses, by design. P itself is
unchanged for linear inputs (guards never touch the weights):

1. **Zero-overlap detection** (`maxCover ≤ 0` today): use total integrated
   overlap **area** (Σ `pr.area`), which is basis-independent and positive.
2. **Coverage / protrusion** (per-node `cover/fullCov ≥ 1−1e-3` today — the
   inequality *flips* for negative row-sums): per-**facet area** coverage
   instead — Σ_fm `pr.area` vs the self-clip `prSelf.area`, refuse below
   `1−1e-3`. Facet-complete coverage ⇒ node-complete tributary integrals —
   equivalent detection **given a non-self-overlapping master surface**, which
   item 5 makes a checked precondition (the area sum counts multiplicity).
3. **Conforming-gap** (per-node `|gap|/cover` today — divides by ~0 at quad8
   corners): the kernel additionally accumulates `gapL1 = ∫|g_N| dΓ` per pair;
   the guard becomes per-facet `gapL1/areaCov ≤ tolFrac·facetScale` with
   `facetScale = 0.5·sqrt(areaFull)` — the 0.5 restores the shipped per-node
   tributary scale (`sqrt(A/4)` for a quad; a bare `sqrt(A)` would have been 2×
   LOOSER than the old guard — review finding). The L1 numerator has no
   cancellation blind spot (an accordion surface whose shared-node signed gaps
   cancel to ~1e-20 reads d/2 in L1 — oracle T11).
4. **Dual mass refusal** (`Ddual[I] ≤ 1e-300` today — would refuse every
   legitimate negative corner): sign-aware AND relative,
   `|Ddual[I]| ≤ 1e-12·max_J|Ddual[J]|` refuses; a negative corner dual mass is
   **correct**, the sign cancels in `P = Mdual/Ddual`. (Relative, because a
   near-cancelled dual mass produces huge P rows that the partition-of-unity
   backstop provably cannot catch — dual rowsums are identically 1.)
5. **Master self-overlap refusal** (review fix MAJOR-1, new): before assembly,
   any pair of master facets with a finite **coincident** mutual overlap (clip
   area > dust, mean gap within the conforming tolerance) is a named refusal —
   a doubly-listed master facet double-counts item 2's area sum and can exactly
   mask an uncovered slave strip (probed: the old per-node guard refused that
   input; the unified guard without this check accepted it). The mean-gap gate
   keeps legitimately curved master surfaces (large-gap shadow overlaps on the
   aux plane) out of the refusal; edge-adjacent facets clip to slivers and are
   skipped.

The two backstops that are already basis-valid stay for all orders: the
post-solve **partition-of-unity** check (`P·1 = D⁻¹D·1 = 1` holds algebraically
for any basis, any sign) and `De.Solve` singularity refusal.

### D4 — quadrature

`N_I^s φ_K^m` products of serendipity bases reach polynomial degree 6 on the
mapped sub-triangles; the shipped degree-4 6-point rule under-integrates them.
`triQuadrature` gains the 12-point degree-6 rule; `generateMortar` passes
order 6 whenever either facet is quadratic (linear pairs keep order 4 —
byte-identical). `MAXGP` 6 → 12. Setup-time cost only. Note the patch test is
quadrature-*independent* (D and M share Gauss points, and the serendipity basis
reproduces linears at midpoint configuration, so `P·(linear)` is exact even
under inexact quadrature) — the higher rule buys D-conditioning and honest
`area`/`gapL1` guards, and the oracle quantifies order-4 vs order-6 (Q0 gate).

### D5 — parser and mode validation

`ltReadFacets` accepts `nps ∈ {3,4,6,8}`. The per-mode guards then own the
policy:

* `generateMortar`: `npsS ∈ {3,4,8}` (tri6-slave named error per D2),
  `npsM ∈ {3,4,6,8}`;
* collocation `generate()` (`LadrunoTie.cpp:348-351`): still refuses `nps > 4`,
  message extended to point at `-mortar` ("collocation onto a quadratic facet
  needs the quadratic closest-point projection — use -mortar, which handles
  quad8/tri6 masters");
* `-shellSolid` untouched (edge polyline, no facets).

## Design-gate BLOCKERs

**BLOCKER-1 — hex20 slave-node lumped mass (the inherited ADR-62 BLOCKER-2,
sharpened).** The projection handler requires positive lumped mass on every tied
slave DOF. The sharp edge: **row-sum lumping of a serendipity element gives
zero/NEGATIVE corner masses** — a hex20 slave surface only satisfies the check
under **HRZ lumping (ADR-35, `LadrunoBrick20 -lumped`)**. Review finding
(MAJOR-2): the shipped check tested `fabs(m) > 0`, so a NEGATIVE assembled
diagonal *passed* and would have poisoned the projection handler (whose own R5
guard scans off-diagonals only). Fix shipped with this ADR: `ltScanMassedDOFs`
now accumulates the SIGNED per-DOF element-mass sum, `ltCheckTiedDofMass`
refuses a negative total with the HRZ remedy by name (all tie modes — a
negative tied diagonal is broken for the handler everywhere), and the mortar
path adds the serendipity hint on any mass refusal of a quadratic slave face.

**BLOCKER-2 — no silent zero-force path for invalid quadratic geometry.** A
curved or midside-displaced quad8 facet must not quietly integrate as its corner
shadow (wrong measure) or quietly skip (lost bond). The −2 status + hard error
(D1) is the gate; the oracle includes a deliberately-curved facet that must be
refused, not mis-integrated.

**BLOCKER-3 — linear regression: P byte-identical, guards honestly changed.**
Every shipped tri3/quad4 mortar test (`tests/test_ladrunoTie_mortar.py`) must
pass unchanged, and the emitted P (weights) must be byte-identical for linear
inputs (order stays 4 for linear pairs; the unified guards never touch the
weights, only refusal behaviour). The guard unification (D3) is a deliberate
behaviour change for linear decks — documented in the ledger, exercised by a
test (a gap field that cancels in the nodal average but fails the L1 check must
now refuse).

## Phased plan (oracle-first, each PR-able)

* **Q0 — the oracle, no build.** `proto_p2_2_quad8_mortar.py`: (T1) serendipity
  shape sanity — Kronecker delta at nodes, Σ N = 1, ∫N_corner = −A/12,
  ∫N_mid = A/3 (quad8) and ∫N_corner = 0 (tri6 — the D2 demonstration);
  (T2) D/M assembly on a **quad8-on-quad4 non-matching patch** (2×2 quad8 slave
  against 3×3 quad4 master, coincident plane) — partition of unity
  Σ_K M_IK = Σ_J D_IJ to 1e-12; (T3) standard P: linear patch exact to ~1e-12;
  (T4) dual condensation: D_dual diagonal with negative corners, P·1 = 1, linear
  patch exact — vs (T5) row-sum lumping breaking the patch (the P2.1 contrast,
  re-proved for quad8); (T6) order-4 vs order-6 quadrature deltas on D
  conditioning and the guards; (T7) the curved-facet refusal case.
* **Q1 — the kernel.** `LadrunoMortarKernel.h` generalization per D1/D4:
  `shapeFull`, `ncorner`, arrays → 8, `gapL1`, status −2, degree-6 rule.
  Verified against Q0 numbers (the oracle's D/M/P are the kernel's expected
  outputs on the same coordinates).
* **Q2 — generator + parser.** `generateMortar` guards per D2/D3, dual sign
  handling, local arrays, `ltReadFacets` + message edits per D5.
* **Q3 — FE tests (build).** Additions to `tests/test_ladrunoTie_mortar.py`:
  hex20-block-on-hex8-block linear patch with `-mortar -dual` (exact, ν = 0 K
  check); tri6-slave named refusal; collocation-quad8 named refusal; curved-quad8
  named refusal; the shipped linear suite green unchanged (BLOCKER-3).
* **Q4 — the cross-check. CLOSED 2026-08-04** (apeGmsh ADR 0086 S1 merged, PR
  #898). Both kernels ran the same quad8-on-quad4 patch: `‖P_dual‖_F` agrees to
  |diff| = **5.0e-13** and the reference corner row to every quoted digit
  (row sum 1.000000000000). Recorded in both repos — the oracle's Q4 block +
  `kinematic_tie_validation/mortar_crosscheck_reference.json` here,
  `test_q4_crosscheck_matches_fork_oracle` there. The contract
  ([[78_apegmsh_mortar_crosscheck_requirements]]) went to **revision 2** on the
  strength of apeGmsh's report, which found four defects in it. The load-bearing
  one: **apeGmsh passed while using a different quadrature rule** (Duffy 5×5, not
  the R3 Dunavant-12), because the Q4 patch is affine and any two degree-6-exact
  rules agree there to ~1e-15 — R3 was unenforceable by its own cross-check. Q4
  is therefore now a **two-case** protocol: case A (the affine patch above) plus
  case B, the same topology bilinearly mapped onto a trapezoid, where the
  pullback goes rational and the two rules part at ~1e-7. Oracle gates T12
  (case B) and T13 (the discrimination measurement); case B's reference `P` is
  published and awaits a second implementation. See [[LEDGER_quirks]].

## Adversarial gate — RUN (2026-08-04), findings fixed

Q1/Q2 is novel math on a shipped, gated kernel — it got the full treatment:
oracle first, then **two independent refute-oriented review agents** (mortar
math; robustness/API) on the C++ delta, plus **mutation testing** of the FE
suite. Results:

* **Mutations (both killed):** (A) a corrupted quad8 midside shape function →
  the 3 quad8 correctness tests fail loudly (surfaces as a 0.33 mean-gap
  refusal), linear tests untouched; (C) reverting `|Ddual|` to the signed
  pre-ADR-78 guard → exactly the two `-dual` quad8 tests fail with the false
  "zero dual mass" corner refusal. The suite discriminates.
* **MAJOR-1 (fixed, D3 item 5):** overlapping/duplicated master facets defeated
  the per-facet area-coverage sum — a false-accept REGRESSION vs the old
  per-node guard (probed with an exactly-masking duplicated master half). Fixed
  with the master self-overlap refusal + regression test.
* **MAJOR-2 (fixed, BLOCKER-1):** negative lumped mass passed `fabs(m) > 0`;
  the HRZ hint was unreachable in the scenario it named. Fixed with the signed
  mass scan + negative-mass named refusal (all tie modes).
* **MAJOR-3 (fixed):** the C++ tri6/quad8 MASTER basis was verified by nothing
  — a cyclic tri6-midside permutation preserves Σφ=1, passes every internal
  gate, and leaves a 0.15 patch error. Fixed with oracle T10 (tri6-master,
  quad8-master, tri3-on-quad8 patches, ≤3e-15) + the FE spring-rig transfer
  tests (`test_tri6_master_linear_field_exact`,
  `test_quad8_master_linear_field_exact`) + an explicit-`LadrunoProjection`
  quad8 test (no test had ever met the handler the emitted message names).
* **MINOR (fixed):** gap-guard scale was 2× looser than shipped (`0.5·sqrt(A)`
  restores parity — D3 item 3); `|Ddual|` threshold shipped absolute instead of
  the signed-off relative form (D3 item 4); the oracle had a dead quadrature
  branch and lacked the kernel's step-converged projection escape; the D1
  "coplanarity guarded" wording was overstated (guard checks midsides only —
  reviewer verified the corner-map identity holds for warped corners, so the
  behaviour matches the shipped warped-quad4 path); collocation silently
  dropped a stray `-slaveFacets` (now a named refusal); the quad8+`-hermite`
  advice loop was circular (message now states no quadratic Hermite exists).
* **Accepted (documented, not fixed):** a degenerate (collapsed-edge) slave
  facet is exempt from the coverage guard via its failed self-clip — it still
  ends in a named D-singular/dual refusal downstream (inherited); the
  `LadrunoContactFE` mortar-contact route would swallow a −2 but is unreachable
  (contact surfaces are nps 3/4 by construction — cross-subsystem invariant,
  noted); refusal tests use bare `pytest.raises` (opserr text is not carried in
  the Python exception); a consistent-mass hex20 passes BLOCKER-1 with positive
  diagonals and is later refused by the handler's R5 with its generic message.

### Post-merge amendment (2026-08-04) — the gate's own "FE 13/13" was pivot luck

`test_quad8_split_bar_equivalence_dual` went red on `ladruno` from the merge of
PR #694 onward and was the only Zone-A failure. **Not a kernel defect** — the
split-bar rig was the one test in the file that did *not* skip `fix(t,0,1,1)` on
the tied facet nodes, so the 8 tied nodes carried an SP on y/z that the tie's own
EQ rows already impose (the masters' y/z are fixed and `Σ P_m = 1`, so
`EQ_row = SP_slave − Σ P_m·SP_master` identically). Under `constraints Lagrange`
that is 16 linearly dependent multiplier rows.

Measured on the merged tip: KKT 176×176, **rank 160 → deficiency exactly 16**
(= 8 tied nodes × 2 redundantly-fixed DOFs), `cond = 1.1e20`; with the redundant
`fix` removed, 160×160 at full rank, `cond = 5.9e6`. The displacement block was
never affected — `|u − U/2| = 4.5e-16` and `σ_xx = 10` exactly in *both* — only
the multipliers were indeterminate. So solvability rode entirely on pivot order:
`numberer Plain` on an MKL desktop build solved (hence the gate's honest-looking
"FE 13/13"), while `RCM` and `AMD` locally, and CI, all failed with
`U(i,i) = 0`. The suite's other 12 tests never tripped it because each already
skips `fix` on `SLAVE_IFACE`, and the projection/transformation handlers condense
the slave DOF rather than multiply it, absorbing the redundancy outright.

Fix: drop the redundant `fix` on the tied facet and parametrize the test over
`numberer` `Plain`/`RCM`, so an over-constrained rig fails loudly rather than
intermittently. **Assertions unchanged** — they were correct throughout, and they
pass tighter on the well-posed system (`4.16e-17`). Lesson recorded in
`LEDGER_quirks.md`; note that the "diff the oracle against the kernel before
touching a red gate" rule had nothing to diff here, because the oracle scopes
itself to kernel math while the disputed branch was enforcement (C++-only).

## Ledger / bookkeeping

* `LEDGER_implementations.md` — extend the LadrunoTie row: quadratic (quad8/tri6)
  mortar facets, ADR-78, PR #.
* `LEDGER_quirks.md` — three entries: (a) *quad8 dual mortar masses are negative
  at corners (−A/12) by construction; sign-aware guards, do not "fix"*; (b) *tri6
  slave surfaces are refused (corner ∫N = 0 ⇒ dual scaling and rowsum coverage
  degenerate); tri6 masters fine*; (c) *ADR-78 unified the mortar-tie guards
  (per-facet area coverage + L1 gap) for ALL facet orders — linear decks with
  cancelling gap fields that previously passed the signed nodal average now
  refuse, by design.*
* `classTags.h` — no change.

## Open questions — RESOLVED at sign-off (2026-08-04)

* **OQ-1 (tri6 slave policy, D2).** → **Refuse in both bases** (as recommended).
* **OQ-2 (guard scoping, D3).** → **Unify all orders** on the sign-free guards
  (against the recommendation, accepting the linear-deck behaviour change for
  one code path and the strictly stronger L1 gap check).
* **OQ-3 (quadrature, D4).** → **Degree-6 for quadratic pairs only** (as
  recommended); linear pairs keep order 4.
* **OQ-4 (ADR home).** → **New ADR-78** (as recommended); ADR-62 gets a
  one-line pointer.
