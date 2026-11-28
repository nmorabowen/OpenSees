---
wp: WP-262
title: "A UniaxialMaterial without DDM methods contributes d(sigma)/d(theta) = 0 to a fiber beam, silently"
pr: "#1172"
date: 2026-11-28
status: "open"
---
### DDM sensitivity of a fiber-section beam ignores a fork fiber material: `sensNodeDisp` misses its whole contribution, no warning

- **Bites:** a DDM run (`parameter` + `sensitivityAlgorithm`) on a `dispBeamColumn`,
  `forceBeamColumn` or Ladruno beam-column with a `Fiber` section whose fibers are a
  material that does not override the sensitivity virtuals. `FiberSection2d/3d` call
  `getStressSensitivity(g, true)` and `commitSensitivity(...)` on every fiber. The
  `UniaxialMaterial` base returns `0.0` for every `get*Sensitivity` and `-1` from
  `commitSensitivity`, and nothing checks that return. So the material's
  contribution is zero. If the parameter belongs to that material,
  `sensNodeDisp` is zero. If it belongs to another one, the response sensitivity
  is wrong, because the stubbed fibers carry no history sensitivity. The run
  reports success either way.
- **Why:** the DDM interface is opt-in per class. A `setParameter` that accepts the
  name (as `LadrunoUniaxialJ2` always did, for `updateParameter` / staging) makes
  the parameter look supported. `activateParameter` storing an id that no
  sensitivity method reads made it look supported a second time.
- **Fix (WP-262), for `LadrunoUniaxialJ2`:** the return map moved into the templated
  `LadrunoUniaxialJ2Kernel.h`. Its double instantiation is bit-identical to the old
  in-class code, checked against a verbatim transcription in
  `tests/ladruno_uniaxialj2_dual_check.cpp`. `getStressSensitivity`,
  `getTangentSensitivity`, `getInitialTangentSensitivity` and `commitSensitivity`
  are each one `Dual<1>` pass. The pass seeds the active parameter and the stored
  history sensitivities (`SHV`, 12 doubles per gradient: εᵖ, p̄, D, ΔD, X_k), plus
  the strain gradient at commit. The local Newton is closed with the IFT step.
  Gates:
  - kernel DDM through history vs FD of the whole path, 1e-9 – 1e-12, for E, sig0,
    Qinf, b and Hiso, with damage and IMPL-EX;
  - dropping the history sensitivities fails by orders of magnitude;
  - end to end: `tests/test_ladrunoUniaxialJ2_ddm_fiber.py`, a cyclic fiber
    cantilever where `sensNodeDisp` equals the FD, with Steel01 as the harness
    control.
- **Still stubbed** (they return 0): every other fork uniaxial (`LadrunoCohesiveHinge`,
  `LadrunoBondSlip`, `LadrunoRebarBuckling`, …) and every fork nD material. Not
  differentiated in `LadrunoUniaxialJ2`: the damage constants and `lchRef`.
  `SHV` is not sent by `sendSelf`, so DDM across OpenSeesSP/MP repartitioning or a
  `database` save/restore is unsupported. Before trusting a DDM run, check that every
  material in the path implements the sensitivity methods.
- **`sensitivityAlgorithm -computeByCommand` is wrong for path-dependent materials**
  (measured error 1.0 for Steel01 and for `LadrunoUniaxialJ2` alike). It computes
  the sensitivities after `commit`, when the trial and committed states already
  coincide. Path-dependent DDM needs `-computeAtEachStep`.
- **A second trap found by the WP-262 review: the active element is stale in the
  DDM pass.** `ops_TheActiveElement` is set only inside `Domain::update`. Any
  sensitivity method that re-reads it (for a crack-band `lch`, say) gets the LAST
  element updated. Latch whatever the trial state used. `LadrunoUniaxialJ2`
  latches `lchScaleTrial`, with a truss-chain regression in
  `tests/test_ladrunoUniaxialJ2_ddm_fiber.py`.

#### Follow-up (WP-267)

`LadrunoCohesiveHinge` and `LadrunoBondSlip` now implement the four sensitivity
methods by the same recipe (templated kernels `LadrunoCohesiveHingeKernel.h`,
`LadrunoBondSlipKernel.h`; PR #1176). `LadrunoBondSlip` also gained the
`setParameter` / `updateParameter` it never had. **Still stubbed** (they return 0):
`LadrunoRebarBuckling` (WP-264 owns the file this wave), every other fork uniaxial,
and every fork nD material, `LadrunoCohesiveHingeBiaxial` included. Not
differentiated in `LadrunoCohesiveHinge`: the parent beam's axial force `P_axial`
(a P-M hinge in a beam misses dP/dθ). `SHV` is not sent by `sendSelf` in either.
