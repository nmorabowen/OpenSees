---
wp: WP-227
title: "GP_PARAM alone places nothing: a per-GP response with no GaussPoint tags records as one flat C1..Cn row (GAUSS_ID -1) (WP-227)"
pr: "TBD"
date: 2026-11-07
---
### GP_PARAM alone places nothing: a per-GP response with no `GaussPoint` tags records as one flat C1..Cn row (GAUSS_ID -1) (WP-227)
- **Bites:** the Ladruno and MPCO recorders learn the GP split of a response
  ONLY from the `GaussPoint` / `NdMaterialOutput` / `ResponseType` tags the
  element emits in `setResponse`. An element that returns
  `new ElementResponse(this, id, Vector(nComp * nGP))` with no tags records one
  element-level column block named `C1..Cn` with `GAUSS_ID = -1`
  (`workaroundForSizeInconsistency`). The group's `QUADRATURE/GP_PARAM` is then
  correct but pairs with nothing, so a viewer still cannot place the values.
  LadrunoLST / SolidShell / UP `stress` had this until WP-227. LadrunoQuad /
  LadrunoCST `stress` / `strain` still do (only their `stressPlaneStrain` is
  tagged).
- **Rule:** a per-GP response opens one `GaussPoint` (`number` = gp+1, in the
  element's own GP order) per point, with an `NdMaterialOutput` and one
  `ResponseType` per component, before returning the response (LadrunoBrick
  `stress` is the pattern). Tags change only metadata: the values and their
  order stay the same.
- **Also:** one class tag that carries several topologies (LadrunoUP) cannot be
  described by the class-tag table. Its arm is geometry Custom + custom rule; the
  element answers `basisInfo`, and the writer derives the geometry per bucket.
  Never route such a class through a standard arm by node count: its GP order
  differs per topology (UP's H8 is not Brick's order).
