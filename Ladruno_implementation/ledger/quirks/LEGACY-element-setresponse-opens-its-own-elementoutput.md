---
wp: LEGACY
title: "Element::setResponse opens its OWN ElementOutput tag — chain to it AFTER output.endTag(), never before"
legacy_seq: 244
---
### `Element::setResponse` opens its OWN `ElementOutput` tag — chain to it AFTER `output.endTag()`, never before
- **Bites:** `Element::setResponse` is not a stub. It serves `force`/`globalForce`/`dampingForce`/`dynamicForce`/`inertialForce` and is the only way an element gets those for free — but it *also* calls `output.tag("ElementOutput")` + `attr(...)` itself. Delegate to it from inside your own open tag and the XML/MPCO stream gets a **nested duplicate** `ElementOutput`, which silently corrupts the column metadata rather than failing.
- **Rule:** `output.endTag(); if (theResponse == 0) return this->Element::setResponse(argv, argc, output);` — end your tag first, then delegate. `LadrunoDispBeamColumn2d` had this right from the start; every other fork element (its own 3d twin included) simply did not chain at all, so `inertialForce` / `dampingForce` / `dynamicForce` were unavailable fork-wide until the token sweep. Mirror it in `getResponse` with `return this->Element::getResponse(responseID, eleInfo);` instead of `return -1` — the base IDs are 111111/222222/333333/444444, so they cannot collide with an element's own.
- *2026-07-28 (recorder-token consistency sweep).*
