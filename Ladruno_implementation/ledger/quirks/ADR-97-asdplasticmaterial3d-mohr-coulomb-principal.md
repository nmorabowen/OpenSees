---
wp: ADR-97
title: "ASDPlasticMaterial3D — Mohr-Coulomb principal-space return (ADR-97 P2)"
legacy_seq: 416
---
## ASDPlasticMaterial3D — Mohr-Coulomb principal-space return (ADR-97 P2)

**`Backward_Euler` reproduces the exact Mohr-Coulomb return ONLY through the
header's own finite difference.** `MohrCoulomb_YF::df_dsigma_ij` and
`MohrCoulomb_PF` both branch on `MC_ds`: `> 0` central-differences their own
`f` / `g` over the six Voigt slots, `== 0` uses an ANALYTIC Lode-angle
`c1/c2/c3` expression. On a proportional face-return deck (E 30000, nu 0.25,
phi 30, psi 10, c 10; trial `sigma = [-10,-40,-100]`) `Backward_Euler` lands on
the exact closest point to **2.1e-14** with `MC_ds = 1e-4`, 1.3e-12 with 1e-6,
1.0e-10 with 1e-8 — and **2.9e-1 away, step-size dependent**, with `MC_ds = 0`.
Mohr-Coulomb's flow direction is CONSTANT inside a sextant (the surface is
piecewise linear), so the cutting plane and the closest point are provably the
same point there; `f` is exactly linear in principal stress, so a central
difference of it is the EXACT 6D gradient and its accuracy *improves* with a
larger step. Two independent references therefore agree against the analytic
coefficients. **Every Mohr-Coulomb deck in `tests/` passes `MC_ds 0.0`**, i.e.
runs on the analytic branch; the ADR-84 MCTC battery is largely insulated
because `special_return` does not use that gradient. Pinned in both directions
by `tests/test_adr97_p2_principal.py::
test_gate4_backward_euler_agrees_only_through_its_own_finite_difference`.
Fixing it changes `Backward_Euler`, so ADR-97 D1 defers it to its own PR.

**An oedometric Mohr-Coulomb deck at `nu = 0.25` with `phi = 30` NEVER YIELDS.**
The elastic lateral-stress ratio `K0 = nu/(1-nu) = 1/3` coincides EXACTLY with
the Mohr-Coulomb compression meridian `(1-sin phi)/(1+sin phi) = 1/3`, so the
uniaxial-strain stress path runs PARALLEL to the yield surface and `f` is
identically `-c cos(phi)` at every stress level, however hard you press. A
tangent or iteration gate built on that rig is silently vacuous: it measures the
ELASTIC operator and passes. Assert the state is plastic before measuring
anything, and pick `nu < 0.25` for `phi = 30` (the P2 gates use 0.15). The same
coincidence exists for any `nu = (1-sin phi)/2`.

**`yf_tolerance()` is not a usable admissibility tolerance for a stress
reassembled from a spectral decomposition.** It defaults to the ABSOLUTE
`f_absolute_tol = 1e-6` (`f_relative_tol` defaults to 0, ADR-94 M5). Recomputing
the header's `f` from `Q diag(y) Q^T` on the ADR-84 MCTC deck (kPa,
`|sigma| ~ 5.4e3`) gives 3.4e-6 — 6e-10 RELATIVE, i.e. round-off — and a check
written against `yf_tolerance()` refuses the step. The round-off is amplified
because an EDGE return lands exactly on a corner, where the Lode angle is ill
conditioned (`dtheta/dJ3 ~ 1/cos(3 theta)` diverges and `dA/dtheta` is not
stationary). Check admissibility where the return was COMPUTED — in principal
space, against all three surfaces of the sextant, which is exact and perfectly
conditioned — and keep the invariant-form check only as a loose guard, at a
tolerance relative to `max(strength_scale, |sigma|)`.

**A load-driven cube rig cannot reach a Mohr-Coulomb FACE state.** Both rigs in
`adr97_oracle/fd_tangent_driver.py` land on `s1 == s2`: `uniaxial` fixes x and y
on every node, so the state is axisymmetric by construction, and `full` (top face
free) passes its limit point under any lateral load before it yields in a
three-distinct-principal state. Measuring a tangent in the face region needs a
kinematically over-determined rig — ADR-97 P2 uses a FREE-NODE rig: the
homogeneous strain field prescribed on seven of the eight nodes, node 7 left
free. Seven prescribed nodes mean no limit point at any stress level, so any
state can be reached, and the 3x3 assembled block at the free node is compared
with a central difference of its own reaction (the same two-rig scheme, since
`setNodeDisp` does not trigger `Domain::update`).
**A `template class` explicit instantiation is the WRONG shape for a
`g++ -fsyntax-only` pre-flight of `ASDPlasticMaterial3D`.** It instantiates
EVERY member of the specialization, including members the real build never
touches because their only call site sits under an `if constexpr` -- e.g.
`cp_apex_return`, which calls `yf.apex_stress()` and therefore fails to compile
for any yield function without an apex (`VonMises_YF`,
`MohrCoulombTensionCutoff_YF`). The result is a page of errors about code that is
correct and unreachable. Instantiate the MEMBERS instead:

    #define private public
    #include ".../AllASDPlasticMaterial3Ds.h"
    typedef ASDPlasticMaterial3D<LinearIsotropic3D_EL, ...> MCMC_t;
    template int MCMC_t::Closest_Point(const VoigtVector&);
    static_assert(MCMC_t::supportsClosestPoint(), "...");

`if constexpr` then discards the unreachable branches exactly as it does in the
real build, and the `static_assert`s turn the support matrix itself into a
compile-time gate (P2 pins all six mixed pairings that way, so a widened family
trait fails at pre-flight instead of at run time). The include set comes from
`adr97_scripts/mk_incs.py` and the build tree's `build/build/Release/build.ninja`
-- note the doubled `build/build`, which `mk_incs.py`'s own usage line does not
say.
