---
wp: WP-195
title: "Caching an element tangent: key it on the raw bits of its inputs, because a material parameter update never reaches the element"
date: 2026-10-14
---
### Caching an element stiffness: a material `updateParameter` never reaches the element, so key the cache on the raw bits of the material tangent

- **What bites.** The obvious invalidation hooks for an element-level stiffness cache are
  the element's own `revertToStart`, `setDomain`, `setParameter` and `updateParameter`.
  The last one is a trap: `parameter 1 element 5 E` makes `Brick::setParameter` forward
  to each Gauss-point material, and the MATERIAL registers itself with the `Parameter`
  (`param.addObject(id, this)`). `updateParameter 1 <value>` then calls the material's
  `updateParameter` directly; `Brick::updateParameter` is never called. An element cache
  invalidated only by its own hooks returns the old stiffness after a sensitivity /
  parameter update -- silently, with a converged (wrong) answer. Likewise `setNodeCoord`
  changes the geometry without touching the element, and a staged material
  (`updateMaterialStage`, a plastic point leaving the elastic range) changes the tangent
  without telling anyone.
- **What works (WP-195, stdBrick).** Never decide "elastic" by class or flag. At every
  tangent call fetch the material tangents anyway (cheap for an elastic material), copy
  their raw bits and the nodal coordinates' bits, and `memcmp` them with the bits that
  built the cached stiffness. Equal bits => the cached matrix is exactly what the full
  loop nest would compute from those inputs (same code, same order) => bit-identical by
  construction. Different bits => rebuild with the same loop nest, re-cache. The element's
  own hooks still invalidate (cheap belt and braces), but correctness does not depend on
  them. Compare BITS, not `==`: `+0.0 == -0.0` and both feed signed zeros into the
  products, so a value-equal compare can return a matrix that differs in the sign of a
  zero.
- **Corollary.** A dirty-flag cache (a "tangent changed" flag set by the material) needs
  every path that can change the tangent to set it -- parameter updates, stage changes,
  `revertToStart`, `recvSelf` -- and ADR-40 lists dirty-flag memoization in the J2
  material as an anti-goal for exactly that reason. The bit compare costs 288 doubles per
  element per tangent call and needs no cooperation from the material.
- **Status.** Tested by `tests/test_wp195_brick_tangent_cache.py` (both parameter routes,
  `setNodeCoord`, `reset`, a yielding J2 block; cache on == cache off, full-field md5).
