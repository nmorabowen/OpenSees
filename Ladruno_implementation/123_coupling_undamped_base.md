# WP-123 — One shared "undamped element" base for the four coupling/embedded elements

Revision 1. Not adversarially reviewed: a behaviour-preserving refactor proven bit-identical and
mutation-tested (below). Per `feedback_adversarial_gate_when`, the full gate is for novel math, core or
vanilla code, or weak coverage; none applies (the vanilla `Element.cpp` is only edited *temporarily* in a
mutation row).

Status: **built and verified locally; draft PR #858.** Next: Zone-A by `workflow_dispatch` on the draft;
mark ready once it is green. The owner merges.

Scoped 2026-09-25. Branch `wp/123-coupling-damp-helper`, cut from `ladruno` @ `bc5c33453` (not stacked on
#855 / #857). Refactor candidate 1 of WP-120 (#855) R3.

## Problem

`LadrunoDistributingCoupling` (RBE3, 33011), `LadrunoKinematicCoupling` (RBE2, 33012),
`LadrunoEmbeddedNode` (33006) and `LadrunoEmbeddedRebar` (33005) are pure penalty couplings with no
physical Rayleigh damping. Each carried its own copy of:

- a no-op `setRayleighDampingFactors`, which refuses the factors so a `betaK` cannot shrink the reported
  explicit `dt_cr` (`CriticalTimeStep` reads the element's factors; ADR 28 §5, ADR 20 §10.6);
- a zero `getDamp` returning an element-owned `C0` (nDOF × nDOF);
- a zero `getRayleighDampingForces` returning an element-owned `dampF`;
- the `C0`/`dampF` members, their allocations, deletes and recv-side re-allocation.

The same defect was fixed in them one after another: the no-op override without a `getDamp` override
hard-crashed every implicit transient. The fix landed in RBE3 (#219, 2026-06-08), then in both embedded
elements (#220, 2026-06-09); RBE2 was born with the override the same day (#221). LEDGER_quirks: "A no-op
`setRayleighDampingFactors` WITHOUT a `getDamp` override".

Reading the code for this WP turned up two further facts:

1. **The `getRayleighDampingForces` copies were dead.** `Element::getRayleighDampingForces` is not virtual
   (`Element.h:207`). Nothing in the four elements calls its own copy, and external callers hold an
   `Element*`, so they reach the base. The copies only shadowed it.
2. **The crash they guarded against was fixed at the root on 2026-07-28.** `Element.cpp` now self-heals
   with the qualified `this->Element::setRayleighDampingFactors(...)` at all 11 sites (LEDGER_quirks
   "makes 11 `Element` methods dereference `theMatrices[-1]`"). The base `getDamp` therefore already
   answers zero for these elements. The per-element `getDamp` is defence-in-depth against that vanilla
   edit being lost in an upstream sync.

## Shape

1. **Baseline first.** Build unchanged `ladruno` (all 5 targets). Record `wp123_undamped/fingerprint.py`:
   56 series, 7,504 exact float reprs. Per element: Newmark and HHT, 40 steps each, with
   `rayleigh(0.5, 1e-3, 1e-3, 1e-3)` and with `rayleigh(0,0,0,0)`, disp + vel of every free DOF, then
   `dampingForce` and `force`. Plus `CentralDifferenceLadruno -cfl` with and without `betaK`:
   self-reported `dtcr`, 40 steps of disp, and `criticalTimeStep()`.
   *Accept:* the new test and the four element batteries pass on the baseline. **Met:** 16/16 new;
   178/178 across the batteries + response tokens + new.
2. **`LadrunoUndampedElement`** (`SRC/element/ladrunoEmbeddedRebar/LadrunoUndampedElement.h`,
   header-only, stamped, in that directory's CMake list). It lives next to `LadrunoEmbeddedKernel`, which all
   four already include; the directory is already in GLOBS and on the global include path. It is
   `: public Element`, with `setRayleighDampingFactors(...) { return 0; }` and `getDamp()` returning a
   per-instance `Matrix` resized to `getNumDOF()` and zeroed. No `getRayleighDampingForces`: the header
   says why.
3. **The four elements** derive from it. Their three declarations, `C0`/`dampF` members, allocations,
   deletes and definitions are removed; each keeps a one-line pointer comment. All four still chain
   `this->Element::commitState()` (L4 passes; `Element` is now an indirect base).
4. **Bit-identity.** Full 5-target rebuild; 0 compiler errors or warnings.
   *Accept:* fingerprint identical to the baseline. **Met: 56/56 series, 7,504/7,504 values identical.**
   178/178 batteries; 10/10 `test_massScaling_validation` + `test_wp103_getstringfromall_tcl`.
5. **New test** `tests/test_ladruno_undamped_couplings.py` (zone_a, 20 cases, < 1 s). The old smoke tests
   only assert "Newmark does not crash"; this pins the contract for all four elements:
   - a damped Newmark/HHT run is **bit-identical** to the undamped run. The models have no nodal masses, so
     the factors could act only through the element;
   - `dampingForce` is exactly zero while moving, and must be non-empty (no vacuous pass);
   - the self-reported `dtcr` and the integrator's `criticalTimeStep()` do not depend on `betaK`;
   - a free vibration under Newmark average acceleration with `algorithm Linear` does not decay. Added after
     mutation row C; see Results.
6. **Mutation rows** (`wp123_undamped/mutation_rows.py`: rebuild `opensees.pyd` only, the five test files
   as separate processes, sources restored in a `finally`), then a final full rebuild. Results below.
7. **Ledgers and guide.** `ladruno-new-element` guide: "Ignoring Rayleigh? derive from
   `LadrunoUndampedElement`"; the old item told authors to override `getRayleighDampingForces`, which does
   nothing. LEDGER_quirks: status + correction on the #219/#220 entry. `LEDGER_implementations` row.

## Results

### Mutation rows

`wp123_undamped/mutation_rows.py` (rebuilds `opensees.pyd` only per row; the five test files run as separate
processes; sources restored in a `finally`). "Batteries" means the four existing element suites
(147 cases).

| Row | Mutation (one edit to the shared code) | New test | Batteries |
|---|---|---|---|
| A | none (final code, full 5-target rebuild) | 20/20 pass | 147/147 pass |
| B | base **stores** the factors (`setRayleighDampingFactors` forwards to `Element`) | **4 fail**: `test_damping_force_is_zero_while_moving` ×4 | pass (**blind**) |
| C | `getDamp` returns a **nonzero** entry (C(0,0) = 10³) | **4 fail**: `test_free_vibration_does_not_decay` ×4 | pass (**blind**) |
| D | no `getDamp` override **and** `Element.cpp`'s `getDamp` self-heal un-qualified (the #219 shape) | **crash** 0xC0000005 | **crash** 0xC0000005 ×4 |
| E | no `getDamp` override, `Element.cpp` fix intact | 20/20 pass: **equivalent mutant** | pass |

What the rows established:
- **The new test closes two holes the batteries could not see (B, C).**
- **Row C first slipped through.** The first version of the test (16 cases) passed it: comparing damped vs
  undamped *in the same build* cannot see a C that is nonzero regardless of the factors. Test 4 was added:
  Newmark average acceleration conserves energy exactly, so a free vibration must not decay.
  - Under Newton it caught only 1/4. These elements put no D·v in their residual, so a spurious C only
    pollutes the tangent, and Newton iterates it away. `EmbeddedRebar` failed only because its tiny
    bipenalty mass made Newton stall.
  - With `algorithm Linear` (one solve per step with the element's own tangent, exact for these linear ties)
    it caught 4/4.
- **Row D reproduces the 2026-06 incident.** The existing smoke tests (added by #219/#220) were already its
  guard.
- **Row E is behaviourally equivalent today.** The base `getDamp` answers zero since the 2026-07-28
  `Element.cpp` fix; the override only matters if that vanilla edit is lost in an upstream sync (D ≠ E).
- **Test 3 (`dt_cr` independent of `betaK`) did not fire under B.** In these models the reported step is the
  elements' self-reported bipenalty bound, which does not read the factors, so the design comment "refuse the
  factors so a βK can't shrink dt_cr" is not observable here. The test is kept as a regression pin on that
  bound, and the doc does not claim more.

## Rejected approaches

- **Free functions in `LadrunoEmbeddedKernel` (WP-120's wording).** A namespace cannot provide virtual
  overrides. Each element would still declare and define the three methods and call the helper, so the
  duplication that let #219/#220 happen stays. A base class makes forgetting impossible.
- **Delete the overrides and rely on the 2026-07-28 `Element.cpp` fix.** It is less code, but the
  elements' safety would then depend on a fork edit to a vanilla file surviving every upstream sync, and
  the four no-op `setRayleighDampingFactors` would stay duplicated. Row E measures this option: it is
  behaviourally equivalent today.
- **Keep a `getRayleighDampingForces` in the base.** It would still be a shadow (not virtual). Keeping it
  would teach the next author that the override does something.
- **Include `LadrunoUP` / `LadrunoRigidBody`.** They also refuse Rayleigh factors, but they answer
  `dampingForce` themselves with their own semantics (LEDGER_quirks "makes 11 `Element` methods") and
  have no replicated-bug history. WP-120 did not name them.

## Open questions

- None blocking. Whether to upstream the 2026-07-28 `Element.cpp` qualification (already recorded as
  upstreamable in `LEDGER_vanilla_files`) is independent of this WP.
