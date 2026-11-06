---
wp: LEGACY
title: "Node::getMass() returns ONLY the nodal mass command value — element-density mass is invisible there"
legacy_seq: 25
---
### `Node::getMass()` returns ONLY the nodal `mass` command value — element-density mass is invisible there
- **Bites:** any code that reasons about a node's mass by reading `Node::getMass()` —
  e.g. the ADR-39 B1 SOFT=1 contact penalty, which needs the gap-mode effective mass
  `m_eff` from the mass the explicit integrator actually inverts. For a model whose mass
  comes from ELEMENT density (a solid `LadrunoBrick`/truss with `-rho`, no per-node
  `mass`), `Node::getMass()` returns a **zeroed** matrix (`Node.cpp:1214` — `mass==0 ⇒
  return a zero matrix`). The first B1 cut sized `m_eff` from `Node::getMass()` and so saw
  `m=0` for every solid-body contact node → silently fell back to the stiff base kn → the
  exact divergence SOFT exists to prevent. Tests that set `ops.mass(...)` directly never
  catch it (found by the B1 adversarial code gate, MAJOR).
- **Why:** OpenSees keeps nodal mass (the `mass` command, stored on the `Node`) and
  element mass (each element's `getMass()`) as **separate** contributions. The global mass
  is assembled from BOTH: the integrator's `formNodTangent → DOF_Group::addMtoTang →
  Node::getMass()` AND `formEleTangent → Element::addMtoTang() → Element::getMass()`. There
  is no `Node::addMass` that folds element mass back onto the node. So `Node::getMass()` is
  the nodal-`mass`-only piece, never the assembled total.
- **Workaround/status (2026-06-24):** to get the mass the explicit solve inverts (the
  assembled global diagonal, for `system Diagonal`), reconstruct it: `m[d] = nodal
  mass(d) + Σ_elements diag(M_e) at that node's translational DOFs`. B1's handler does this
  once per `handle()` (`ladrunoBuildNodalMass` in `LadrunoContactHandler.cpp`) and caches it
  on `LadrunoContactDomain` for the stateless adapter. Matches `diag(global M)` for `system
  Diagonal` (the `DiagonalSOE` default extracts `M(i,i)`, no row-sum); a row-sum/lumped
  distributed SOE (OpenSeesMP `MPIDiagonal`) would differ — soft contact is serial-only
  today. Translation-first DOF order (`[u | θ]`) makes a node's first 3 element DOFs
  translational.
