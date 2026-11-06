---
wp: ADR-63
title: "ADR-63 #4a — averaged nodal-normal smoothing (NTS)"
legacy_seq: 145
---
## ADR-63 #4a — averaged nodal-normal smoothing (NTS)

**The auto global-sign vote uses the master-surface centroid over UNIQUE nodes — NOT the flat `mTags`
connectivity.** `LadrunoContactSurface::getNodeTags()` returns the flat per-segment connectivity, in which
edge/ridge nodes shared by K segments appear K times. Averaging that flat list double-counts the
high-valence ridge nodes and biases the "master centroid" toward the ridge; for a `-smoothNormal` master
whose slave cloud sits near that biased centroid plane the auto seed `slave_centroid − master_centroid` can
flip sign → the whole nodal-normal field points INWARD → `gap = n_smooth·(x_s−x̄) > 0` reads "not
penetrating" → silent pass-through (looks exactly like the R3 bug the feature is meant to fix). Caught on the
convex-ridge gate during P1 bring-up (a 90°-ish tent: flat-list centroid z=0.5 ≈ slave z=0.5 → seed_z≈−1e-3 →
G=−1). Fix: dedupe master node tags (a `std::set<int>`) before averaging (`LadrunoContactHandler.cpp`,
ADR-63 field-build block). The auto sign is still fragile when the slave cloud straddles the master centroid
plane — use `-outward` for such masters (the documented escape). The per-segment Newell area-normal vote
(`Σ σ_s·newellAreaNormal(s)`) is NOT affected (it sums over SEGMENTS, each once).

**ADR-63 #4a — sharp-ridge facet ownership (a smoothed-normal limitation).** A slave sitting AT a
sharp convex ridge has its closest-point projection land on the SHARED EDGE of the *adjacent* facet
(barycentric ≈ 0, marginally in-bounds), so that neighbor reads a large SPURIOUS penetration and
injects a big ejecting force. This is a pre-existing NTS facet-ownership issue: the FACETED path only
*incidentally* prunes it (at a 90° ridge the neighbor normal is ~⟂ the per-pair `orientDir`, so
`normalOriented`'s perpendicular fail-safe kills the pair) — the SMOOTHED path (which does not use the
per-facet `orientDir` for the normal) has no such prune. A blunt interior-margin gate in
`evalSegmentSmooth` was tried and REVERTED: it removed a spuriously-*helpful* neighbor force and
destabilized other geometries (the spurious force helps at some ridge angles, hurts at others). Proper
fix = closest-facet / interior facet ownership at a shared edge (a P2 item, related to ADR-57 #4b
edge-handoff). For P1: the R3 SIGN fix is validated by the quad convex-ridge gate; tri-3 chain coverage
uses a FLAT patch to avoid the pathology; a slave pressed onto a facet INTERIOR (away from ridges) is
unaffected.

**ADR-63 P2.1 — the sharp-ridge ownership fix is GAP-AWARE, not a near-edge reject (2026-07-01).** The
above limitation is RESOLVED, but the obvious fix is a TRAP. A per-segment shared-edge mask
(`segmentSharedEdges`, topological, cached with σ) is threaded into `evalSegmentSmooth`; but rejecting
ANY projection that lands within a parametric band of a SHARED edge reproduces EXACTLY the reverted
blunt fix and causes PASS-THROUGH: a frictionless slave slides UP-slope to the ridge apex (the smoothed
normal is more vertical than the facet normal, so a facet-perpendicular load tilts up-ridge), and at the
apex BOTH facets project onto the shared edge ⇒ both rejected ⇒ the slave is driven straight through
(regressed the P1 ridge gate to min_d=−338). The spurious double-activation is LOAD-BEARING at the apex.
The working rule is GAP-AWARE: reject a shared-edge projection ONLY when `−gap > edgeGapFrac·(|g1|+|g2|)`
(`edgeGapFrac=0.05`) — the spurious non-owner reads a penetration ∝ its LATERAL distance from the ridge
(large), while the true owner AND a genuine at-apex contact read a SMALL gap (kept). `edgeGapFrac`/
`edgeTol=0.1` are dimensionless (gap-vs-facet-size / parametric) ⇒ robust across ridge angle + scale.
FREE/boundary edges are never rejected ⇒ R5 slide-off untouched. Residual: near the apex the non-owner
is still kept (a harmless small force ⇒ mild over-stiffness at the exact ridge, NOT pass-through); full
closest-facet selection under sliding = ADR-57 #4b, a P2.3 concern.

**ADR-63 P2.1 — an ill-posed R3 test can pass for the WRONG reason (2026-07-01).** The P1 R3 gate
(`test_p1_smoothnormal_holds_the_ridge_facet`) originally ran the slave LATERALLY FREE and only
"passed" because the P1 spurious ridge ejection (the very bug) pinned `min_d ≥ 0` early. A frictionless
slave under a load with a lateral component on a CONVEX ridge has NO equilibrium — it slides up-slope and
launches — so a laterally-free rig cannot test "held". Once P2.1 removes the ejection the free slave
correctly slides off (min_d≪0). Lesson: to gate a SIGN/normal claim, constrain the confounding DOF —
FIX the slave's x,y (only z free) so the test states the well-posed thing (smooth stays repulsive,
min_d≈−1e-3; faceted flips, min_d≈−318). A "held" assertion over a frictionless free body on a curved
master is suspect.

**ADR-63 build — a deep worktree path overflows the cl.exe command line (2026-07-01).** Building
OpenSeesPy inside `…\.claude\worktrees\<name>\` fails at the include-heavy TUs (OpenSeesPy / SparsePython)
with `CreateProcess failed. The parameter is incorrect.` / `ninja: fatal: CreateProcess` — the long
absolute path, repeated across ~50 `-I` flags, blows the Windows ~32 KB command-line limit. Fix (in
`build.bat` configure): `-DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON` pushes include/object/library lists into
`.rsp` files. Harmless on short paths. NB: adding this flag to an existing build tree triggers a full
recompile (every compile rule's command hash changes).

**ADR-63 P2.2/P2.3 — friction composes with `-smoothNormal` for FREE, but the AUTO sign still needs
`-outward` for curved masters (2026-07-01).** (1) `LadrunoContactFE::segmentActive` builds the friction
slip via `LadrunoFrictionKernel::tangentPart(drel, n, gTvec)` with the SAME `n` the gap operator uses — the
smoothed normal when `useSmoothNormal` — so friction is projected against `n_smooth` with NO new code; on a
flat master (n_smooth==n_facet) a frictional slide is byte-identical smooth-vs-faceted. `-reemit`/
`-smoothNormal`/`-mu` compose (parser refuses only vs `-mortar`), and `-reemit -smoothNormal -mu -outward`
sustains a frictional crossing over a convex curved master (the ADR-60 "exposed combo" closed) — WITHOUT
`-reemit` it still passes through (smoothing fixes the SIGN, re-emit fixes the SEARCH; orthogonal). (2)
TRAP: the ADR-60 R3 `-outward` caveat is lifted only when the global sign vote is WELL-CONDITIONED. A single
slave *starting to the side* of a curved arc votes a near-horizontal seed (slave − master-centroid ≈ ⟂ the
up-field) whose tiny z-component can flip the sign INWARD ⇒ the smoothed field points inward ⇒ pass-through
even with `-smoothNormal` (the pre-existing F2/F3/F5 low-confidence warning fires). So `-smoothNormal` is NOT
a blanket lift of `-outward` for curved masters — it lifts it only for slave clouds sitting OVER the master
(seed ∥ field). Always pass `-outward` for edge-grazing / side-approaching slaves on a curve. **[SUPERSEDED
by P2.5 below — the auto sign is now a robust per-slave majority vote that holds for over-the-master edge/
side-approaching clouds without `-outward`; only a genuinely two-sided or multi-shell cloud still needs it.]**
(3) The P2.1
gap-aware guard's near-apex over-stiffness under SLIDING is a mild quality effect (a small extra bump as a
block crests a SHARP ridge at speed; negligible on realistic shallow arcs, maxpen ~0.01; never diverges),
not a pass-through — full single-owner selection is ADR-57 #4b.

**ADR-63 P2.4 — the frozen-field smoothed tangent CONVERGES implicitly; the dropped `∂n_smooth/∂u` is
`O(kn·gN)` ⇒ sub-dominant on a penalty contact (2026-07-01).** The Q-IMPLICIT-NEWTON tripwire (does Newton
converge with the SUPPRESSED `∂n_smooth/∂u`, i.e. the frozen-field symmetric `kn·BᵀB`, on a genuinely
curved implicit master?) resolved to **outcome (a): converges** — 2 Newton iterations per step, INDEPENDENT
of load-step coarseness (even a single step dragging a slave across the whole facet = maximal within-step
normal rotation still converges in 2). The reason is structural: the dropped consistent-tangent term is
`kn·gN·∂²gN/∂u²`, i.e. scaled by the penalty PENETRATION `gN ≈ press/kn`, which is small in any well-posed
penalty contact ⇒ it never dominates the kept `kn·BᵀB` (this is the SAME reason the shipped faceted default
drops its own B3 block by default). The ONLY regime where smoothed Newton iterations climb (a swept
soft-penalty misuse, `gN ≳ 15%` of the facet) is exactly where the FACETED `-geomtan` consistent tangent
ALSO diverges (and even the seat step fails) — a penalty-NTS-breakdown, NOT a smoothed-normal defect, so
P3's full `∂n_smooth/∂u` would not rescue it. **⇒ P3 (`-consistentNormalSmooth`) stays a genuinely-optional,
evidence-deferred follow-up, not a required item.** Rig gotchas (why the test looks the way it does): (a) a
lone frictionless slave on a convex master has ZERO lateral contact stiffness at the crest (vertical n ⇒
`kn·nx²=0`) and no lateral equilibrium off-crest ⇒ the seat solve is singular/runaway — use a weak lateral
spring (ks≪kn) + seat at the SYMMETRIC centre; (b) a single combined-load `(Fx,0,-P)` under
DisplacementControl is DEGENERATE for a frictionless slave (one load factor scales both the press and the
drive ⇒ equilibrium pins to a single slope) — decouple the constant press via `loadConst` + a separate
lateral drive pattern; (c) DisplacementControl is the ONLY implicit displacement driver for a
`constraints LadrunoContact` model — the handler REFUSES a non-homogeneous (imposed-displacement) SP; (d)
there is NO validated static+`-reemit` path (every ADR-60 reemit test is explicit/CDL), so a FIXED master
(constant nodal-normal field, no re-handle needed) sidesteps it for the implicit rig.
