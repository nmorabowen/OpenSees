# WP-124 — Shared Element-contract helpers for the continuum element shells

Revision 0 (scoping). Refactor candidate 2 of WP-120 (#855) R3.

Status: **scoping; draft PR.** No production code changed yet. Implementation waits for the inventory
(below) and for WP-123's local builds to finish (one build machine).

Scoped 2026-09-25. Branch `wp/124-continuum-shell-helpers`, cut from `ladruno` @ `bc5c33453` (not stacked
on #855 / #857 / #858).

## Problem

WP-120 R3 family 1: eight fork continuum elements (LadrunoQuad, LadrunoCST, LadrunoLST, LadrunoCSTPair,
LadrunoBrick, LadrunoBrick20, BezierTri6, BezierTet10) share an `Element`-contract shell (1,002 exact /
2,365 renamed-identifier shared lines). Eight PRs had to fix the same defect in several copies:

| PR | Defect | Copies |
|---|---|---|
| #562 | Rayleigh P-clobber in `getResistingForceIncInertia` | Quad, CST, LST, CSTPair |
| #228 | `getInitialStiff` not returning the cached `*Ki` | Quad, CST |
| #670 | mass cache not invalidated in `recvSelf` | 8 (incl. SolidShell) |
| #683 | EAS inner-Newton warning spam | Brick, Quad |
| #588 | EAS degeneracy guard blind to axis collapse | Brick, Quad |
| #709 | unsymmetric plastic tangent symmetrised | BezierTet10, BezierTri6 |
| #224 | `setParameter` not forwarded to GP materials | BezierTet10, BezierTri6 |
| #852 | ground-motion inertia sign | BezierTet10, BezierTri6 |

## Shape (draft — to be fixed after the inventory)

1. **Inventory** (read-only): the six shared parts × eight elements, the exact operation order per copy,
   and which copies can share one helper **bit-identically**.
2. **Helpers**: free functions (the elements have different bases and DOF layouts), in a header such as
   `SRC/element/LadrunoElementShell.h`.
3. **Staging**: plane four first (they already share `LadrunoFiniteStrain2DKernel.h`), then Bezier, then
   Brick/Brick20. One commit per stage; each stage proven bit-identical (baseline fingerprint vs refactor)
   with a test that fails on a one-line mutation of the shared helper.

## Rejected approaches

- **A common base class.** The eight derive from `Element` with different node counts, DOF layouts and
  mass/stiffness caches; a base class would force one layout or a template hierarchy. That is a bigger
  change than the defects justify.

## Open questions

- Which parts cannot be unified without changing results (operation order, mass source): from the inventory.
