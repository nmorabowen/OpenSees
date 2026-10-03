---
wp: LEGACY
title: "A CONSTANT perturbed-pivot count from the elastic stage onward is a mesh-topology symptom — check for element-less (orphan) nodes before blaming the element or…"
legacy_seq: 472
---
### A CONSTANT perturbed-pivot count from the elastic stage onward is a mesh-topology symptom — check for element-less (orphan) nodes before blaming the element or the material
- **Bites:** Pardiso (or MUMPS) reports the same number of perturbed pivots on every factorisation, starting at the first elastic K0 solve and never changing through the plastic stages. It looks like an element rank defect, but an element defect would change with the material state and would show up as near-null elastic modes of the element stiffness. WP-114 could not reproduce the TIMs report's 1 816 pivots on any Bézier deck (0 pivots across every bisection, no near-null modes). The one time the harness did produce such pivots, the cause was tri6 mid-edge nodes that no element referenced (their DOFs have exactly zero stiffness), which is the likely cause in a mesh export too.
- **Why:** a node with no element has an all-zero row and column (plus any fix/sp). A direct solver has to perturb each of its free DOFs, and that count depends only on the topology, so it is constant from the first step.
- **Census recipe (one line, after the model is built):** `used = {n for e in ops.getEleTags() for n in ops.eleNodes(e)}; orphans = [n for n in ops.getNodeTags() if n not in used]`. `2 × len(orphans)` (ndf=2, minus fixed DOFs) is the pivot count to expect. Remove the orphans or fix their DOFs.
- **Status (2026-09-18):** diagnostic only, WP-114 (`tests/wp114/pivot_bisection.py`). Not chased further; the reporter's mesh export is the suspect.
