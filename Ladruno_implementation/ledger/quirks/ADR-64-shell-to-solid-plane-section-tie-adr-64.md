---
wp: ADR-64
title: "Shell-to-solid plane-section tie (ADR-64 -shellSolid): the paid-for gotchas"
legacy_seq: 149
---
## Shell-to-solid plane-section tie (ADR-64 `-shellSolid`): the paid-for gotchas

The P4 rung ties an ndf-6 shell EDGE (master polyline) to an ndf-3 solid FACE (slaves) with
rigid plane-section arm rows `u_s = Σ N_j (u_j + θ_j×d)` (three mixed EQ rows per solid face
node; drilling drops out via the 1e-12 filter). What we paid to learn:

- **Cross-ndf EQ rows are `Transformation`-INCOMPATIBLE** — the handler's factorization goes
  singular (`U(i,i)=0`). Enforce with `Lagrange` (static) or `LadrunoProjection` (explicit)
  only. Documented-unsupported, no investigation planned (ADR-64 OQ-6).
- **Honest kinematic limits, GATED not hidden (ADR-64 OQ-2):** the 3-translation rigid arm
  suppresses through-thickness Poisson stretch at the seam (St-Venant boundary layer, misfit
  = ν·ε·|z| ∝ ν·t, EXACT at ν=0 — oracle T6), and Timoshenko shear warping leaves an O(γ·t)
  local misfit (analytic max γc/(3√3)) that does NOT shrink with in-plane refinement while
  resultant transfer stays exact (oracle T7). Model at ν=0/thin seams or accept the layer.
- **Explicit rotary-mass rule (CORRECTED from the ADR plan):** the generator has NO rotary
  precondition (slaves are solid translations, `-rho` covers BLOCKER-2; static Lagrange
  needs nothing) — but under `LadrunoProjection` the master edge's θx/θy still need a small
  nodal rotary mass (`mass $n 0 0 0 mr mr mr`): the projector keeps every group DOF in the
  explicit equation set, so a zero-mass DOF has no equation of motion (`rigidLink -beam`
  has the identical requirement). The "tied solid supplies the edge's rotary inertia via
  CᵀM_ccC" claim in the ADR draft was measured WRONG for the shipped projector. ALSO: the
  shell element must LUMP mass on the tied DOFs — `ShellMITC4::getMass()` is CONSISTENT
  (off-diagonal) and the handler's element guard refuses it; `ASDShellQ4` lumps
  translational mass and works.
- **`mass` sizes by the BUILDER's ndf, not the node's** (`OPS_addNodalMass` uses
  `OPS_GetNDF()`): in a mixed ndf-3/ndf-6 model, `ops.mass(shellNode, 6 values)` after a
  `model('basic','-ndf',3)` builds a 3×3 and dies with `Node::setMass - incompatible
  matrices` even though the node was created with the per-node `-ndf 6` override. Re-issue
  `ops.model('basic','-ndm',3,'-ndf',6)` before creating/massing the shell side.
- **Collocation force-transfer needs NESTED grids for an exact force patch:** the kinematic
  rows are exact on ANY grid (oracle T1–T5), but the transmitted interface FORCES only match
  the receiving side's consistent nodal pattern when the slave grid nests in the master's
  (2:1 mid-splits etc. — the same property P1's patch test used silently). A non-nested
  split (e.g. 0.4 vs 0.5) redistributes interface forces at the few-% level and fails a
  tight stress patch. Consistent transfer on non-nested grids = the mortar variant
  (deferred: needs the plane-section basis inside the M integral, a kernel per-GP hook).
- **Testing against a stale `dist/bin` binary: EQ_Constraints survive `wipe()`.** The main
  checkout's old pyd predates the ADR-30 `theEQs->clearAll()` fix in `Domain::clearAll`, so
  a second model built in the same process inherits the previous model's EQ rows → duplicate
  (linearly dependent) Lagrange multiplier rows → `U(i,i)=0`. Current source is fixed; the
  trap only bites pre-flight scripts run against an outdated binary — one model per process
  there.
**A `pos < n` bounded staged walk can "pass" its own post-check — pre-walk the totals (2026-07-02).**
The lumped SMS injection walked element nodes with `for (...; pos < n)` and then rejected non-node-major
layouts with `if (pos != n)`. For an element whose nodes' TOTAL ndf EXCEEDS its mass size n, the bounded
loop stops MID-NODE at exactly pos == n, the post-check sees n, and the misaligned mass is silently
committed. The consistent sibling was immune only because its first pass summed ndf UNBOUNDED before
comparing. Same family of trap next door: the consistent M̄ builder indexed `mdiag[base[a]+d]` for
d < min-ndm without clamping by each node's OWN ndf — an ndf < ndm node spills the write into the next
node's block (off the end of Mbar/mdiag on the last node). Fixed (review-P3): pre-walk Σndf and reject
≠ n up front; clamp `ndmOf[a] = min(ndm_a, ndf_a)`. LESSON: a walk over a matrix whose layout a DIFFERENT
object (the nodes) defines must validate the FULL mapping before mutating anything — loop bounds are not
validation. Related: an iterative-solver "did not converge" warning must key on the RESIDUAL, not on
iters == maxIt — the consistent PCG's pAp≤0 SPD-breakdown guard exits EARLY (iters < maxIt, resid > tol)
and the old `iters >= maxIt` condition swallowed it silently (also fixed review-P3).
