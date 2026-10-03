---
wp: WP-99
title: "Element refusal roster — who acts on a material's setTrialStrain return code (full audit, WP-99)"
legacy_seq: 447
---
### Element refusal roster — who acts on a material's `setTrialStrain` return code (full audit, WP-99)

- **Bites:** a material returns a failure code from `setTrialStrain` and nothing
  happens — or it happens on one element of your model and not the one beside
  it. There is no engine-wide contract here at all: each element's `update()`
  decides on its own, and roughly half of them throw the code away.
- **Why:** `setTrialStrain` predates any fail-loud convention. An element either
  accumulates the codes into the value its `update()` returns (so the algorithm
  sees a failed state determination and the step is cut), or calls
  `setTrialStrain` from a `void form*` routine / assigns the result to a variable
  it never reads / has no `update()` override at all (so `Element::update()`
  returns 0 and the refusal is invisible).
- **The audit (2026-09-14, on `9c2f964ea`+WP-99).** Every `.cpp` under
  `SRC/element` that mentions both `setTrialStrain` and `NDMaterial`, classified
  by what its own `update()` does. **52 elements: 26 FORWARD, 1 SENTINEL-only,
  25 DISCARD.** Reproduce with
  `grep -rln setTrialStrain SRC/element --include=*.cpp` and read each
  `update()`.
  **If you script that grep, do not key the definition on `)` followed by `{`.**
  Six of these files write the opening brace *below a comment line*
  (`int\nSSPquad::update(void)\n// this function updates ...\n{`), and three more
  give the class a different name from the file (`Nine_Four_Node_QuadUP.cpp`
  defines `NineFourNodeQuadUP`). A first pass of this audit reported all nine as
  "no `update()` override" — the verdicts happened to stay right, because those
  six return 0 anyway, but the stated *reason* was wrong for six rows until
  re-verification caught it.
  - **FORWARD** — any nonzero code reaches the return of `update()`, so ANY
    material refusal cuts the step.
  - **SENTINEL** — only `LADRUNO_MATERIAL_REFUSED` cuts the step; every other
    nonzero code is ignored. Deliberate, per ADR-33/34 (ASDConcrete3D's negative
    "best-state" codes must not fail a step). `LadrunoBrick` is the only one.
  - **DISCARD** — the code cannot reach the analysis at all, either because
    the element has no `update()` override (so `Element::update()` returns 0),
    or because its `update()` calls `setTrialStrain` and returns 0 regardless,
    or because `setTrialStrain` is reached only from a `void form*` routine.
    The table says which.

  **THIS TABLE IS THE ONE AUTHORITATIVE COPY.** The `opserr` strings in
  `LadrunoSANISAND.cpp` / `ManzariDafalias.cpp` and the guides name EXAMPLES and
  point here — a second closed copy of this list is exactly how the previous,
  wrong one survived in three documents at once.

| element | verdict | evidence | first `setTrialStrain` |
|---|---|---|---|
| `BBarFourNodeQuadUP` | **FORWARD** | `update()`@327: `ret += ...setTrialStrain(...)`, `return ret` | `:369` |
| `BezierTet10` | **FORWARD** | `update()`@370: `ret += ...setTrialStrain(...)`, `return ret` | `:408` |
| `BezierTri6` | **FORWARD** | `update()`@393: `ret += ...setTrialStrain(...)`, `return ret` | `:461` |
| `ConstantPressureVolumeQuad` | **FORWARD** | `update()`@362: `success += ...setTrialStrain(...)`, `return success` | `:502` |
| `E_SFI` | **FORWARD** | `update()`@599: `errCode1 += ...setTrialStrain(...)`, `return errCode1` | `:617` |
| `E_SFI_MVLEM_3D` | **FORWARD** | `update()`@781: `errCode += ...setTrialStrain(...)`, `return errCode` | `:798` |
| `EightNodeQuad` | **FORWARD** | `update()`@386: `ret += ...setTrialStrain(...)`, `return ret` | `:437` |
| `FourNodeQuad` | **FORWARD** | `update()`@576: `ret += ...setTrialStrain(...)`, `return ret` | `:615` |
| `FourNodeQuad3d` | **FORWARD** | `update()`@384: `ret += ...setTrialStrain(...)`, `return ret` | `:424` |
| `FourNodeQuadUP` | **FORWARD** | `update()`@359: `ret += ...setTrialStrain(...)`, `return ret` | `:419` |
| `FourNodeQuadWithSensitivity` | **FORWARD** | `update()`@346: `ret += ...setTrialStrain(...)`, `return ret` | `:385` |
| `LadrunoBrick20` | **FORWARD** | `update()`@945: `ret += ...setTrialStrain(...)`, `return ret` | `:972` |
| `LadrunoCST` | **FORWARD** | `update()`@212: `ret += ...setTrialStrain(...)`, `return ret` | `:234` |
| `LadrunoLST` | **FORWARD** | `update()`@252: `ret += ...setTrialStrain(...)`, `return ret` | `:271` |
| `LadrunoQuad` | **FORWARD** | `update()`@710: `ret += ...setTrialStrain(...)`, `return ret` | `:534` |
| `LadrunoUP` | **FORWARD** | `update()`@856: `ret += ...setTrialStrain(...)`, `return ret` | `:916` |
| `Nine_Four_Node_QuadUP` | **FORWARD** | `update()`@507: `ret += ...setTrialStrain(...)`, `return ret` | `:572` |
| `Nine_Four_Node_QuadUPOld` | **FORWARD** | `update()`@235: `ret += ...setTrialStrain(...)`, `return ret` | `:266` |
| `NineNodeQuad` | **FORWARD** | `update()`@392: `ret += ...setTrialStrain(...)`, `return ret` | `:446` |
| `SFI_MVLEM` | **FORWARD** | `update()`@767: `errCode1 += ...setTrialStrain(...)`, `return errCode1` | `:785` |
| `SFI_MVLEM_3D` | **FORWARD** | `update()`@893: `errCode += ...setTrialStrain(...)`, `return errCode` | `:911` |
| `SixNodeTri` | **FORWARD** | `update()`@352: `ret += ...setTrialStrain(...)`, `return ret` | `:397` |
| `TenNodeTetrahedron` | **FORWARD** | `update()`@1029: `success += ...setTrialStrain(...)`, `return success` | `:1213` |
| `Tri31` | **FORWARD** | `update()`@550: `ret += ...setTrialStrain(...)`, `return ret` | `:586` |
| `Twenty_Eight_Node_BrickUP` | **FORWARD** | `update()`@819: `ret += ...setTrialStrain(...)`, `return ret` | `:983` |
| `Twenty_Node_Brick` | **FORWARD** | `update()`@427: `ret += ...setTrialStrain(...)`, `return ret` | `:509` |
| `LadrunoBrick` | **SENTINEL** | `update()`@985 tests `== LADRUNO_MATERIAL_REFUSED` | `:1034` |
| `AC3D8HexWithSensitivity` | **DISCARD** | `update()`@267 calls it and drops the code (`return 0`) | `:289` |
| `BbarBrick` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:951` |
| `BBarBrickUP` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:1021` |
| `BbarBrickWithSensitivity` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:965` |
| `BeamContact2D` | **DISCARD** | `update()`@386 calls it and drops the code (`return 0`) | `:466` |
| `BeamContact2Dp` | **DISCARD** | `update()`@375 calls it and drops the code (`return 0`) | `:466` |
| `BeamContact3D` | **DISCARD** | `update()`@580 calls it and drops the code (`return 0`) | `:694` |
| `BeamContact3Dp` | **DISCARD** | `update()`@454 calls it and drops the code (`return 0`) | `:563` |
| `Brick` | **DISCARD** | `update()`@912 calls it and drops the code (`return 0`) | `:1069` |
| `BrickUP` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:1069` |
| `EmbeddedEPBeamInterface` | **DISCARD** | `update()`@688 calls it and drops the code (`return 0`) | `:755` |
| `EnhancedQuad` | **DISCARD** | `update()`@1289 does not call it (called from a void `form*` routine); returns 0 | `:1077` |
| `FourNodeTetrahedron` | **DISCARD** | `update()`@973 calls it and drops the code (`return 0`) | `:1144` |
| `IGAKLShell` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:3165` |
| `IGAKLShell_BendingStrip` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:2445` |
| `LadrunoDispBeamColumn3d` | **DISCARD** | `update()`@646 does not call it (called from a void `form*` routine); returns 0/err/solveHingeJump(v, L)/solveHingeJumpBiaxial(v, L) | `:832` |
| `LadrunoSolidShell` | **DISCARD** | `update()`@274 does not call it (called from a void `form*` routine); returns 0 | `:670` |
| `NineNodeMixedQuad` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:1003` |
| `SimpleContact2D` | **DISCARD** | `update()`@368 calls it and drops the code (`return 0`) | `:432` |
| `SimpleContact3D` | **DISCARD** | `update()`@470 calls it and drops the code (`return 0`) | `:552` |
| `SSPbrick` | **DISCARD** | `update()`@402 calls it and drops the code (`return 0`) | `:445` |
| `SSPbrickUP` | **DISCARD** | `update()`@369 calls it and drops the code (`return 0`) | `:412` |
| `SSPquad` | **DISCARD** | `update()`@404 calls it and drops the code (`return 0`) | `:426` |
| `SSPquadUP` | **DISCARD** | `update()`@360 calls it and drops the code (`return 0`) | `:382` |
| `ZeroLengthND` | **DISCARD** | no `update()` override at all -> `Element::update()` returns 0 | `:383` |

- **The u-p family is the one to notice.** Every vanilla `*QuadUP` / `*BrickUP`
  element that has its own `update()` FORWARDS (`FourNodeQuadUP`,
  `BBarFourNodeQuadUP`, `Nine_Four_Node_QuadUP`, `Twenty_Eight_Node_BrickUP`),
  while `BrickUP` and `BBarBrickUP` (no `update()` override) and `SSPquadUP` /
  `SSPbrickUP` (an `update()` that calls `setTrialStrain` and returns 0
  regardless) DISCARD. "the UP family swallows refusals" was stated in four
  fork documents and is wrong for half of them — and u-p is SANISAND's canonical
  host, so it is the half that matters.
- **Two fork edits are already in the FORWARD column** and are easy to mistake
  for vanilla behaviour: `TenNodeTetrahedron` (`success +=`, the TIMs report
  item 8 fix) and `LadrunoUP`.
- **AT COMMIT TIME THE TABLE IS IRRELEVANT: nothing propagates.**
  `Domain::commit()` is `elePtr->commitState();` with the return value dropped,
  for every element in the table. See the entry
  "`Domain::commit()` discards element commit returns" above for what WP-99 does
  about it (a material declares the refusal out of band and `Domain::commit()`
  aborts).
- **Workaround/status (2026-09-14):** none of the vanilla DISCARD elements is
  fixed — `return success` on `Brick` is an unconditional behaviour change for
  every material (`tests/test_adr84_p2a_strict_convergence.py::test_stdbrick_swallows_the_refusal`
  pins it). Pick a FORWARD element for any gate whose meaning depends on a
  refusal being seen.
