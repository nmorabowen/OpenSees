---
wp: WP-99
title: "stdBrick/BbarBrick/BrickUP swallow the material's refusal — a THIRD silent accept, in a vanilla ELEMENT (found ADR-84 §6c finding 2; status updated ADR-86b; ro…"
legacy_seq: 302
---
### `stdBrick`/`BbarBrick`/`BrickUP` swallow the material's refusal — a THIRD silent accept, in a vanilla ELEMENT (found ADR-84 §6c finding 2; status updated ADR-86b; **roster corrected WP-99**)

- **TITLE AND ROSTER CORRECTED 2026-09-14 (WP-99 / F7).** This entry was headed
  "`stdBrick`/`BrickUP`/`QuadUP`" and that third name was **wrong**:
  `FourNodeQuadUP::update()` does `ret += theMaterial[i]->setTrialStrain(eps)`
  (`FourNodeQuadUP.cpp:419`) and therefore **PROPAGATES** any nonzero code.
  `BbarBrick` (`:951`), `SSPbrick` (`:445`), `SSPquad` (`:426`) and
  `LadrunoSolidShell` (`:670`) are the names that belong on the list instead, and
  `stdBrick` **is** `Brick` (`TclBrickCommand.cpp:210`). Full audited
  classification in the `Domain::commit()` entry at the end of this ledger.
- **Bites:** any material that refuses a trial strain (a strict-mode ASDP
  rejection, a `ManzariDafalias`/`LadrunoSANISAND` substep-cap refusal, ...)
  hosted in `stdBrick` (= `Brick`), `BbarBrick`, `BrickUP`, `SSPbrick`,
  `SSPquad` or `LadrunoSolidShell`. `Brick::update()` writes
  `success = ...->setTrialStrain(strain);` and then **unconditionally**
  `return 0;` — the code is assigned and never read. `BbarBrick`/`BrickUP` call
  `setTrialStrain` inside a *void* `formResidAndTangent`, so there is no return
  path for the code at all. So the material refuses, prints its `opserr` line,
  returns a failure code — and the analysis reports success regardless of what
  that code was.
- **Why:** these are vanilla elements, written before any Ladruno material
  needed a fail-loud contract; nobody expected `setTrialStrain` to return
  anything worth checking. Found while measuring ADR-84 P2a's
  `strict_convergence` gate (§6c finding 2) — `TenNodeTetrahedron` already
  accumulates and returns the sum (a pre-existing fork fix, TIMs report item
  8), which is why the ADR-84 contract tests use the tet host instead of
  `stdBrick`.
- **Deliberately NOT fixed for `stdBrick`:** `return success` is an
  unconditional behaviour change for every `stdBrick` + every material, which
  is precisely the blast radius an opt-in refusal contract exists to avoid.
  Pinned by `tests/test_adr84_p2a_strict_convergence.py::test_stdbrick_swallows_the_refusal`
  so that fixing `Brick.cpp` shows up as a loud, informative test failure
  rather than a mystery elsewhere.
- **Workaround/status (2026-09-05, ADR-86b review-fix):** **`LadrunoBrick` now
  propagates the sentinel (`LADRUNO_MATERIAL_REFUSED`) on ALL FIVE `update()`
  paths, including `updateHypo` and `formEAStrue`** (ADR-86b's original repair
  covered four of five; the review pass confirmed `updateHypo`/`formEAStrue`
  are also sentinel-aligned, not the blanket `< 0` an earlier ledger row
  mistakenly claimed — see `LEDGER_implementations.md`'s ADR-86b row). **`stdBrick`,
  `BbarBrick`, `BrickUP`, `SSPbrick`, `SSPquad` and `LadrunoSolidShell` still
  swallow the refusal, unchanged** — this defect is not fixed on any vanilla
  element, only worked around by using an element that forwards the code
  (`LadrunoBrick`, `LadrunoBrick20`, `LadrunoQuad`/`CST`/`LST`,
  `BezierTet10`/`Tri6`, `FourNodeQuad`, `FourNodeQuadUP`) for every gate that
  needs the return code to mean something.
