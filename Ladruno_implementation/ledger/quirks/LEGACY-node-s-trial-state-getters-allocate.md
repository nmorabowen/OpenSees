---
wp: LEGACY
title: "Node's trial-state GETTERS allocate"
legacy_seq: 453
---
### `Node`'s trial-state GETTERS allocate

`Node::getTrialDisp()` (`Node.cpp:590`), `getTrialVel()`, `getTrialAccel()` and
friends lazily call `createDisp()` / `createVel()` / `createAccel()`, which
`new` a `4*numberDOF` array and build four `Vector`s over it. A read-only-looking
`const Vector &d = theNodes[a]->getTrialDisp();` in an element's `update()` is
therefore a **write** the first time it runs on a node. Two elements sharing a
fresh node both allocate: last-writer-wins, the other allocation leaks, and a
`const Vector&` already handed out points into the freed buffer.

Exactly three `create*` functions exist and all lazy getters funnel through them,
so a serial pre-pass touching those three getters on every node closes the hazard
completely — which is what `Domain::ladrunoThreadedUpdate()` does before entering
its parallel region.
