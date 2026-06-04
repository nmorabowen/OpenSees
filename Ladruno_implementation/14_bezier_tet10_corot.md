# BezierTet10 — `-geom corot` (LOCKED plan)

Status: **LOCKED** (2026-06-02), post adversarial review (5 lenses, 30 findings,
18 confirmed / 6 refuted). Scope = corotational geometry method only. **Finite
strain is OUT OF SCOPE** (separate task: per-element UL assembly + the
`FiniteStrain`/`LogStrain` material adaptors).

Goal: give `BezierTet10` a `-geom corot` option by consuming the existing
[`SolidTransformation`](../SRC/element/solidTransformation/SolidTransformation.h)
/ `SolidTransformationCorot` machinery, mirroring how `LadrunoBrick` (8-node hex)
already consumes it. `BezierTet10` stays `NEN=10, NELD=30, NSTRESS=6`.

---

## 0. What the review SETTLED (refuted findings — do NOT re-litigate)

- **`MAX_NODES 8→10` is the only transformation change and is safe.** Every loop
  in `SolidTransformationCorot.cpp` runs on the runtime `n`/`nn`/`ndof`; the
  `MAX_NODES`/`MAX_DOF` macros only *size* fixed arrays (`Xrel/xrel`, `Gamma`,
  the `Rf` scratch) and bound-check `update()`. The polar, spin map Γ, and EICR
  force/stiffness are node-count-agnostic. `SolidTransformationCorot` has exactly
  one consumer today (LadrunoBrick, n=8), which is unaffected.
- **Rigid-body zero-stress identity holds for the 10-node cloud** (it is an
  algebraic identity of the polar fit, independent of node distribution).
- **Unweighted Procrustes "bias" / deferred-geometric-term "worse for quadratic"**
  — refuted as defects: accurate code facts, no demonstrated correctness problem.
  The element's lumped volume share is uniform `Ve/10`, so unweighted nodal
  Procrustes is defensible here. (Still covered by the FD-tangent gate, §5.)
- **Null-ctor `theGeom`** — mirror LadrunoBrick: the null ctor pre-sets
  `theGeom = new SolidTransformationLinear()`; `recvSelf` rebuilds it.

---

## 1. Shared machinery (one line)

`SRC/element/solidTransformation/SolidTransformationCorot.h` — enum
`MAX_NODES = 8` → `10` (auto-bumps `MAX_DOF` to 30). No logic change.

---

## 2. Element wiring in `BezierTet10` (mirror LadrunoBrick)

- **Member / lifecycle**: add `SolidTransformation *theGeom`, include
  `SolidTransformation.h`. Full ctor gains a **trailing default param**
  `int geomMethodID = SolidTransformation::METHOD_LINEAR` (preserves every
  existing call site) → `theGeom = SolidTransformation::create(geomMethodID)`
  with `SolidTransformationLinear` fallback. Null ctor sets
  `theGeom = new SolidTransformationLinear()`. Destructor deletes it.
- **`computeLocalDisp()`**: build `refCrds(10,3)`, `curCrds(10,3)`,
  `uGlobal(30)`; `theGeom->update(NEN, refCrds, curCrds)` then
  `localizeDisp(uGlobal, uCore)`; return `uCore(30)`. (Buffers sized `NELD`, not 24.)
- **`update()`**: replace the direct `getTrialDisp` loop with
  `uCore = computeLocalDisp()`; build strain from `uCore`. `-geom linear` ⇒
  identity ⇒ bit-identical to today.
- **`getResistingForce()`** — see §3 for the load-frame contract.
- **`getTangentStiff()`** — see §3 (fCore) and §4 (`getInitialStiff`).
- **`sendSelf`/`recvSelf`** — see §4 (buffer is full, must widen).
- **`commitState`/`revertToLastCommit`/`revertToStart`**: corot is stateless;
  base no-ops suffice (still forward the calls to `theGeom` for safety).

---

## 3. Load-frame contract (the one genuine DESIGN FORK — DECIDED)

The reviewers split on how body force / pressure / Q interact with the rotation.
**Decision (diverges from LadrunoBrick, intentionally — and is more correct):**

> **`fCore` = pure internal force `Bᵀσ` only.** Body force, the `+z` pressure
> hack, and `Q` (applied/inertia load) are all applied in the **GLOBAL frame
> AFTER `globalizeForce`** — they are fixed-direction loads and must NOT be
> co-rotated by `R`.

Rationale: the corotational geometric stiffness `K_geo` arises from the
configuration-dependence of `R` acting on the **internal** force only. Fixed
external loads are `R`-independent, so they correctly contribute **zero** to
`K_geo`. LadrunoBrick folds body force into the rotated path (giving gravity a
follower character + a spurious body-force `K_geo` term) — a small, validated
effect, but Tet10 takes the cleaner contract. Bonus: with `fCore = Bᵀσ` in both
paths, f/K consistency is **automatic**.

Concrete code shape:

- Factor a private helper **`computeCoreInternalForce(Vector &fInt)`** = `Bᵀσ`
  over the GP loop, **core frame, no external loads**. Call it from BOTH
  `getResistingForce` and `getTangentStiff` so the `fCore` fed to
  `globalizeStiff` is byte-identical to what `globalizeForce` rotates.
- `getResistingForce()`: `computeCoreInternalForce(fInt)` →
  `theGeom->globalizeForce(fInt, P)` → then in the **global frame** subtract
  body force (`b`/`appliedB`), the `+z` pressure term (move it OUT of the GP
  loop), and `Q` (already global at `BezierTet10.cpp:650`).
- `getTangentStiff()`: assemble `K` core-frame, `computeCoreInternalForce(fInt)`,
  `theGeom->globalizeStiff(K, fInt, K)`.
- **Frame refresh**: both `getResistingForce` and `getTangentStiff` refresh the
  frame (`theGeom->update(NEN, ref, cur)`) at the top before reading
  `getStress()`/`getTangent()`, so f and K assemble against one fresh `R` even
  outside the normal `update()`-first integrator contract.
- **v1 guard**: REJECT `-geom corot` when `pressure != 0` (the `+z` volume hack
  is unvalidated under rotation) with a clear parse-time message.
- **Mass**: lumped/consistent mass stays in the **global frame, unchanged** —
  corot does not de-rotate mass. Document.

---

## 4. `getInitialStiff` + serialization (confirmed plumbing fixes)

- **`getInitialStiff()`** (was under-specified). Before `globalizeStiff`, build
  `refC(10,3)` and `curC = refC` from node reference coords and call
  `theGeom->update(NEN, refC, refC)` to **pin `R = I`**, THEN
  `globalizeStiff(K, zeroForce, K)`. Without the update-to-reference,
  `globalizeStiff` would use a STALE current-frame `R`. The `Ki` memo is safe
  because `R=I` is deterministic on every entry — but the refresh MUST run on the
  first-compute path.
- **`sendSelf`/`recvSelf`**: the `iData` ID buffer is **FULL (15/15)**. Widen
  `ID iData(15)` → `iData(16)` in both; `sendSelf` writes
  `iData(15) = theGeom->getMethodID()`; `recvSelf` rebuilds
  `theGeom = SolidTransformation::create(iData(15))` with `Linear` fallback for
  unknown/0. **Do NOT use the `+1000` packing trick** (no packed slot exists;
  a dedicated index is correct). `dData` stays at 5.
- **Broker registration (discovered during build, NOT in the original plan).**
  The serialization round-trip test surfaced a **pre-existing latent defect**:
  `BezierTet10` (and its sibling `BezierTri6`) were never registered in
  `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`, so database/MPI
  `recvSelf` failed with *"no Element type exists for class tag 33001"*. Added
  `#include` + `case ELE_TAG_BezierTri6/BezierTet10: return new …()` (vanilla
  file → `// Ladruno` hooks + `LEDGER_vanilla_files.md` row). This is exactly the
  class of bug the round-trip test was added to catch.

---

## 5. OPS parser (`OPS_BezierTet10.cpp`)

- Add an explicit `else if (strcmp(opt,"-geom")==0 || strcmp(opt,"-geometry")==0)`
  branch: check remaining args, read the token, map `linear`→`METHOD_LINEAR`,
  `corot`/`corotational`→`METHOD_COROT`; **explicit `return 0` with a clear
  message** for `finite` and any unrecognized value (do NOT fall through to a
  warn-only path). Pass `geomMethodID` to the ctor.
- **Note**: Tet10's formulation axis is a **bool `-bbar`**, not a `std|bbar`
  enum. Both default (std) and `-bbar` are allowed under corot. (`-pressure`
  rejected under corot in v1, §3.)

---

## 6. Test matrix (Zone-A pytest; the FD-tangent gate is load-bearing)

Mirror `tests/test_ladrunoBrick_corot.py`. Required:

1. **Rigid-body rotation → zero stress/force** (std AND `-bbar`).
2. **`-geom linear` bit-identical to no `-geom`** (~1e-12) — seam proof.
3. **Small-strain patch**: corot == linear at infinitesimal disp.
4. **FD-tangent gate** (THE arbiter of the pure-`Bᵀσ` Option-A assembly):
   `eleResponse('stiff')` vs central-differenced `eleForce`; gap → 0 ~linearly as
   strain→0 (gates `R k_d Rᵀ` + `G1+G1ᵀ`), with a documented loose bound at finite
   strain. PLUS a direct cross-check that `fCore` inside `getTangentStiff` equals
   `getResistingForce`'s internal force before the global load subtractions.
   (mirror `test_ladrunoBrick_corot.py:240,263`)
5. **`getInitialStiff` R=I guard**: `update()` at a deformed config, then
   `getInitialStiff()` == reference-config initial stiffness.
6. **Large-rotation cantilever** — REPLACE the invalid hex-vs-tet pointwise
   oracle with a **mesh-refinement convergence gate**: corot Tet10 tip
   displacement converges (monotone/Cauchy) toward an analytic large-rotation
   elastica reference as `nz` increases (or a numpy corot-beam oracle at loose
   tol). Energy/monotonicity/convergence, NOT pointwise cross-element match.
7. **bbar under corot**: (a) rigid rotation + `-bbar` → zero force; (b)
   near-incompressible (ν≈0.4999) `-bbar`+corot stays strictly more flexible than
   std+corot under rotation + small deviatoric strain (volumetric split survives
   the rotation strip). (mirror `test_ladrunoBrick_element.py:213`)
8. **Serialization round-trip** with `-geom corot` via
   `tests/_testbed/roundtrip.py` `database_roundtrip()` (build+analyze, save,
   wipe, restore, assert committed disp) — catches a forgotten geom-id pack that
   would silently revert corot→linear on a parallel worker (no MPI needed).
9. **Existing battery unchanged**: `tests/test_bezierTet10_element.py` passes
   verbatim with the no-arg ctor (run before/after the wiring change).
10. **(low/optional)** explicit-dynamics smoke: corot Tet10 under central
    difference in free rigid rotation stays stress-free. Document: corot Tet10 v1
    is verified quasi-statically; explicit use rides the unchanged lumped-mass +
    corot resisting-force path; recommend **αM-only** damping in explicit runs
    (βK uses the rotated/K_geo tangent and shrinks dt_cr — see LEDGER_quirks).

---

## 7. Ledgers / banner

- `LEDGER_implementations.md`: row for BezierTet10 corot (**no new classTag** —
  reuses 33001).
- `banner_features.txt` line if shipped → `patch_banner.py`.
- The `MAX_NODES` edit touches only the fork-authored `SolidTransformationCorot.h`
  (already header-stamped) ⇒ **no `LEDGER_vanilla` row**.
- No new files unless the helper is split out; `BezierTet10.*` already stamped.

---

## 8. Effort / sequencing

Transformation change ~1 line. Element wiring ~1 day (validated pattern). The
only real engineering was §3 (load-frame fork — DECIDED) and §4 (getInitialStiff
refresh + buffer widen). `finite` remains the separate, larger follow-up.
