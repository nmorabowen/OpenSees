---
wp: LEGACY
title: "FourNodeQuad-family element rho is NOT serialized and is uninitialized by the broker ctor → garbage mass on database/restore (and OpenSeesMP) → non-determinist…"
legacy_seq: 181
---
### FourNodeQuad-family element `rho` is NOT serialized and is uninitialized by the broker ctor → garbage mass on `database`/`restore` (and OpenSeesMP) → non-deterministic, indefinite M on transient restart — **UPSTREAM bug, FIXED (see [[LEDGER_vanilla_files]])**

The continuum plane elements `FourNodeQuad`, `Tri31`, `NineNodeQuad`,
`EightNodeQuad`, `SixNodeTri` (and `FourNodeQuad3d`) carry an element-level
density member `rho` and build their mass matrix as
`if (rho == 0) rhoi = theMaterial[i]->getRho(); else rhoi = rho;`. Upstream
`sendSelf`/`recvSelf` pack `[tag, thickness, b0, b1, pressure, alphaM, betaK,
betaK0, betaKc, ...]` but **never `rho`**, and the no-arg constructor the
`FEM_ObjectBroker` uses on restore initialized `thickness`/`pressure` but **left
`rho` uninitialized** (the one exception, `FourNodeQuad3d`, happens to zero it).

Consequence: after `ops.database('File',p); save; wipe; restore` (or any
OpenSeesMP inter-process send), the reconstructed element's `rho` is garbage
heap memory. Garbage ≠ 0, so the `if (rho==0)` gate FALLS THROUGH and the mass
matrix is built from the garbage value instead of the material density —
**non-deterministic across fresh processes** (uninitialized read; measured
denormalized/`-2.5e-297` and `-2.5e6` generalized eigenvalues), an **indefinite
M**, and a transient (Newmark) restart that diverges from the uninterrupted run
(measured rel. diffs 0.5 … 5.6 across identical fresh runs).

Why it hid so well:
- The **committed nodal disp/vel/accel round-trip bit-exact** (FE_Datastore stores
  them fine), so a plain nodeDisp `database_roundtrip` PASSES — the corruption
  only surfaces once the CONTINUED transient uses the mass matrix.
- The **material `rho` is restored correctly**; it's the element's *shadow* `rho`
  that's garbage and *overrides* the material via the `if (rho==0)` gate.
- With **no LoadPattern / no excitation** the restart is bit-exact (u≡0, so a
  wrong M is never exercised) — which masks the bug and made it look
  LoadPattern-triggered. It is not; it's mass-triggered.

Diagnostic recipe: build a small dynamic plane-strain model WITH a load, run a
few Newmark steps, `save→wipe→restore`, continue, and compare to the
uninterrupted run (or `ops.eigen` before-vs-after — negative eigenvalues = smoking
gun). See `tests/test_quad_tri_rho_db_restart.py`.

Fix (strictly additive, this PR): append `rho` to each element's send/recv data
Vector (+1 slot) AND add `rho(0.0)` to each blank constructor. Vanilla bug —
every mainline OpenSees user doing a database or parallel dynamic restart with
these elements is affected; upstream-PR candidate.

**General lesson (this is the 2nd fork bug of this exact shape):** any scalar
member an element/material reads at analysis time MUST be either serialized in
`sendSelf`/`recvSelf` OR initialized in the no-arg broker constructor — ideally
both. A `save→wipe→restore→continue` (not just `→compare committed disp`)
round-trip test is the cheap guard; the plain nodeDisp round-trip is blind to
any state that only feeds future steps (mass, damping, committed internal vars).
