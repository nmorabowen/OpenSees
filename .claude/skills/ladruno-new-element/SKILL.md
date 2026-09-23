---
name: ladruno-new-element
description: >
  Checklist for adding or changing an ELEMENT in the Ladruno OpenSees fork (SRC/element/,
  fork-authored classes such as LadrunoBrick/Quad/CST/LST/UP/SolidShell, Bezier*, couplings).
  Use before writing or modifying getResistingForceIncInertia, getTangentStiff, getMass,
  getDamp/Rayleigh, sendSelf/recvSelf, setResponse/getResponse, an OPS_ parser, or when
  registering a new element classTag. Every item points to the LEDGER_quirks entry that
  explains it.
---

# New or changed element — checklist

Read this before adding or changing an element under `SRC/element/`. Each item names the
`Ladruno_implementation/LEDGER_quirks.md` heading to grep for; read that entry when the item
applies. Items marked **[lint]** are enforced by `python ci/check_quirk_patterns.py`.

## Registration (a new element)

- [ ] classTag in `SRC/classTags.h` (fork band), recorded in `LEDGER_implementations.md`;
      `python ci/check_classtags.py` clean.
- [ ] Register in ALL dispatch sites: `classTags.h`, `FEM_ObjectBrokerAllClasses.cpp`, the
      `functionMap` in `OpenSeesElementCommands.cpp` (Python), AND the Tcl table in
      `TclElementCommands.cpp`. Quirks: "A new element needs registering in THREE dispatch sites".
- [ ] Manifest row + a test (`python ci/check_manifest.py`); header stamp (add the files to
      `GLOBS` in `Ladruno_scripts/stamp_headers.py`, then run it); banner line via
      `banner_features.txt` → `patch_banner.py`. See `AGENTS.md`.

## Dynamics, mass and damping

- [ ] **[lint]** In `getResistingForceIncInertia`, snapshot the shared static residual into a
      LOCAL before adding `getRayleighDampingForces()`: betaK Rayleigh calls `getTangentStiff()`,
      which may refill the static and silently drop inertia (and `-Q`). Quirks: "MUST snapshot
      the shared static `resid`" — read its 2026-07-11 recurrence note.
- [ ] Any element with mass gets a **dynamic Rayleigh regression test** (betaK ≠ 0, transient).
      Same entry: "a dynamic Rayleigh regression is mandatory".
- [ ] Ignoring Rayleigh? Override `getDamp` AND `getRayleighDampingForces` too. Quirks: "A no-op
      `setRayleighDampingFactors` WITHOUT a `getDamp` override" and "makes 11 `Element` methods
      dereference `theMatrices[-1]`".
- [ ] Quadratic/serendipity/T6: nodal-lumped corner masses can be zero or negative; use HRZ.
      Quirks: "T6 quirks", "runs 8/27 mass-deficient".
- [ ] `rho` and every construction input are serialized in `sendSelf`/`recvSelf` and
      zero-initialized in the broker ctor. Quirks: "element `rho` is NOT serialized",
      "FileDatastore silently CLOBBERS", "`recvSelf` into a LIVE element".

## State and re-entrancy

- [ ] `static Matrix`/`static Vector` scratch returned by reference is non-re-entrant (matters
      under the OpenMP element loop). Quirks: "Returning `const Matrix &` from `getTangentStiff()`".
- [ ] Materials read `ops_TheActiveElement` (e.g. crack-band `lch`); it is a global written per
      element. Quirks: "`ops_TheActiveElement` is a mutable GLOBAL", "Crack-band materials read
      element size" (the base `getCharacteristicLength` is wrong for high-order elements).
- [ ] Never iterate the Domain from inside an element callback. Quirks: "`Domain::getElements()`
      is a SHARED singleton iterator".
- [ ] Size every `static Vector` in `getResponse` exactly: `Vector::operator()` is unchecked in
      release. Quirks: "`Vector::operator()` is UNCHECKED".
- [ ] `setResponse`: chain to `Element::setResponse` AFTER `output.endTag()`. Quirks:
      "`Element::setResponse` opens its OWN `ElementOutput` tag".
- [ ] Refuse degenerate geometry with a scale-free metric that also catches axis COLLAPSE, and
      test the refusal directly. Quirks: "A Jacobian degeneracy guard must normalize by the
      LARGEST column".

## Parser (`OPS_<Element>`)

- [ ] Never peek optional numeric tokens with `OPS_GetString` (returns the literal
      `"Invalid String Input!"`); a failed `OPS_GetIntInput` consumes the arg in Python but not
      in Tcl. Quirks: "`OPS_GetString` returns the literal", "A failed `OPS_GetIntInput`".
- [ ] Don't use `OPS_GetNDM()` as a dimension oracle. Quirks: "`OPS_GetNDM()` is NOT a safe".

## Validation

- [ ] Surface loads on quadratic/Bezier elements: corner weights can be negative, and the
      base-reaction identity cannot see the error. Quirks: "uniform surface pressure has NEGATIVE
      corner weights", "Lagrange-consistent surface loads on control DOFs".
- [ ] Break each new gate on purpose once (revert the fix, confirm it fails). Quirks: "A test can
      be GREEN because of the very bug".

Found a new trap? Add it to `LEDGER_quirks.md`, then add one line here pointing to it. If the
trap has a greppable pattern, add a rule to `ci/check_quirk_patterns.py` instead.
