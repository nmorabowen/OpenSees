# WP-116 — Stamp LadrunoDispBeamColumn so the quirk lint scans it

Revision 1. Follow-up to WP-115 (#850); its adversarial review raised the gap.

Status: **built; draft PR #851** (2026-09-23).

Scoped 2026-09-23. Branch `wp/116-dispbeam-stamp`, cut from `ladruno` @ `79e062367` (not stacked on
#850). The quirk lint `ci/check_quirk_patterns.py` lands with #850; until then it is run here from
WP-115's copy with `--root`.

## Problem

`ci/check_quirk_patterns.py` scans only fork-authored sources, identified by the
`LADRUNO-HEADER-START` stamp. `LadrunoDispBeamColumn2d/3d` (ELE_TAG 33013/33014) were fork classes
without the stamp and absent from `Ladruno_scripts/stamp_headers.py` `GLOBS`, so the L1 Rayleigh check
never saw them. The WP-115 reviewer called their Rayleigh sites safe, unverified.

## Shape

1. **Stamp.** Files are fork-authored: added by fork commits `07666f7ff` (2D) and `bc5a11f19` (3D),
   classTags 33013/33014, `LEDGER_implementations` row "LadrunoDispBeamColumn2d / 3d". They live in
   `SRC/element/ladrunoDispBeamColumn/` (not `dispBeamColumn/`). Added to `GLOBS`; four files stamped;
   `stamp_headers.py --check` clean (227 files). The parser `OPS_LadrunoDispBeamColumn.cpp` was already
   stamped.
2. **Lint.** L1 then reports five sites: `getResistingForceIncInertia` twice per class (rho ≠ 0 and
   rho = 0 branches) accumulating Rayleigh into the class-shared `static Vector P`, and 2D `getResponse`
   id 12 (`dampingForces`) doing `P.Zero(); P.addVector(Rayleigh)`.
3. **Trace (independent of the WP-115 review).** Using the lint's own C++ function parser on both
   classes: `P` is written only by `getResistingForce`, `getResistingForceIncInertia`, `getResponse` and
   `getResistingForceSensitivity`. From `Element::getRayleighDampingForces()` the class methods reached
   are `getMass`, `getTangentStiff`, `getInitialStiff`, `getBasicStiff`; none writes `P` or calls
   `getResistingForce*`. The only callback from sections/materials into the element is the lch
   channel (`ops_TheActiveElement->getCharacteristicLength()`), which returns a double
   (`current_section_lch` or the base fallback) and writes nothing. Verdict: all five were **safe**.
4. **Conversion anyway** (so safety stops depending on that invariant): `getResistingForceIncInertia`
   accumulates into a function-local `static Vector res(6|12)` seeded by `res = getResistingForce()`,
   same operations in the same order, `return res`; `getResponse` id 12 gets its own `static Vector damp`.
   L1 is clean for both files afterwards.
5. **Test** `tests/test_rayleigh_inertia_dispbeam.py` (zone_a, 34 cases). No transient Rayleigh test
   existed. Differential against `elasticBeamColumn` (elastic section + 3 Legendre IPs is the elastic
   beam): {alphaM, betaK, betaK0} × {lumped element mass, nodal mass} × {step, UniformExcitation} plus
   consistent `-cMass` step legs, 2D and 3D. The lumped/nodal oracle carries nodal masses and the
   `-cMass` oracle runs step loads only, because vanilla `ElasticBeam2d` with element mass subtracts
   the ground-motion Q twice (LEDGER_quirks, found by WP-115). Plus a closed-form check that
   `dampingForces` equals `betaK·K_e·v_e` from the analytical Euler–Bernoulli stiffness.

## Results

(filled in below after the build and the break-on-purpose runs)

## Rejected approaches

- **Waive the five sites.** They are safe today, but a waiver keeps safety dependent on nobody ever
  making `getTangentStiff` touch `P`. The conversion is bit-identical, so it costs nothing.
- **Stack on #850** to get the lint file. Stacked PRs strand their merges; cut from `ladruno` and run the
  lint with `--root` instead.
- **`elasticBeamColumn -mass` as the ground-motion oracle.** Vanilla `ElasticBeam2d` double-counts Q.

## Open questions

- None specific to this WP.
