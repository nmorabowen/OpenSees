# WP-232i — B-spline basis, boundary treatment, MMS paths (round 1: the freeze candidate)

Status: **proposed** by the 232i C++ worker, 2026-11-16; binding once the orchestrator says
"frozen". Authority: ADR 232 (§4, §5.2, §5.3 + 232f amendments, §5.4, §6, §8 P3, §9 phase 3) and
the "232i" brief (K1–K4, gate table). Anchors quoted from `origin/wp/232i-mpm-bspline-mms` @
`d93ce4dd4`; line numbers drift. `tests/mpm_reference.py` and `tests/mpm_implicit_reference.py`
were not read. Decisions are numbered (B*n*) so a review can cite them. Conventions are
`tests/mpm_interface.md` §1 (tension positive, row-major F, Voigt with engineering shear strain).

## 1. The basis

- **B1 — the 1D factors.** Node-centred uniform B-splines on the existing vertex lattice
  `x_i = x0 + i h` (knots at `x_i ± h/2` for degree 2, at the nodes for degree 3). With
  `r = x_p − x_i`, `ξ = r/h`, and `dS` = d/dx_p (the sign convention of interface §4.1):

  | basis | region | `S` | `dS` |
  |---|---|---|---|
  | `bspline2` | `\|ξ\| < 1/2` | `3/4 − ξ²` | `−2ξ/h` |
  | | `1/2 ≤ \|ξ\| < 3/2` | `½(3/2 − \|ξ\|)²` | `−sign(ξ)(3/2 − \|ξ\|)/h` |
  | `bspline3` | `\|ξ\| < 1` | `2/3 − ξ² + \|ξ\|³/2` | `(−2ξ + (3/2)ξ\|ξ\|)/h` |
  | | `1 ≤ \|ξ\| < 2` | `(2 − \|ξ\|)³/6` | `−sign(ξ)(2 − \|ξ\|)²/(2h)` |
  | both | beyond | 0 | 0 |

  Cardinal values (G1): quadratic `S(0) = 3/4`, `S(±h/2) = 1/2`, `S(±h) = 1/8`, `S(±3h/2) = 0`
  (the brief's "S(±h/2) = 1/2 and 1/8" is read as these two values); cubic `S(0) = 2/3`,
  `S(±h) = 1/6`. Both are C¹ at every knot (cubic C²), and `S`, `dS` are 0 at the support end,
  so a region boundary needs no tie rule. Tensor product exactly as interface §4.3. `lp` is
  ignored (as for `linear`). Source: de Boor, *A Practical Guide to Splines* (cardinal B-splines);
  the MPM use: Steffen, Kirby & Berzins (2008).
- **B2 — support rule.** On axis d, `c = floor((x_p − x0_d)/h)`, candidates `c−1 … c+2`, entries
  with `S = dS = 0` dropped — the existing GIMP scan, unchanged:
  `for (int i = c - 1; i <= c + 2; i++) {` … `if (S == 0.0 && dS == 0.0) continue;`
  (`LadrunoMPMCore.h:200`, `:204`). Quadratic keeps 3 nodes per axis (`round(ξ)−1 … +1`, a subset
  of the candidates), cubic 4; `AxisSupport` already holds `int idx[4]` (`:171`).
- **B3 — one ghost layer suffices, admissibility.** A particle is admissible iff its centre lies in
  the closed real box on every axis — the `linear` rule (`admissibleAxis`, `:182-183`), `lp` not
  checked. From `x_p ≥ x0`, node `−2` is at `≥ 2h`: outside the quadratic support (`3h/2`) and on
  the cubic support's edge (`S = dS = 0`). So `Grid::ghost = 1` (`:122`) holds both degrees, node
  ids are unchanged (interface §2), and partition of unity, `Σ∇S = 0` and linear reproduction
  `Σ S_i x_i = x_p` hold exactly on the extended lattice at every admissible position. A B-spline
  grid needs `n_d ≥ 2` on every axis (refused at `mpm grid`), so a ghost's image (B5) is never a
  node of the same face.
- **B4 — plumbing.** `enum Basis { BASIS_GIMP = 0, BASIS_LINEAR = 1, BASIS_BSPLINE2 = 2,
  BASIS_BSPLINE3 = 3 }` (`:74`); `bsplineFactor(int degree, double r, double h, double &S, double
  &dS)` beside `gimpFactor` (`:81`); `axisSupport` (`:187`) dispatches; `particleSupport` (`:212`)
  is untouched. One `basisName(Basis)` helper replaces every two-way ternary. Every site that
  switches on the basis today:

  | Site | Today | 232i |
  |---|---|---|
  | `OPS_LadrunoMPM.cpp:176-180` (`mpm grid -basis`) | `gimp`, `linear` | + `bspline2`, `bspline3`; `n_d ≥ 2` |
  | `OPS_LadrunoMPM.cpp:430-437` (`debug shape1d`) | always `gimpFactor` | `[-basis b]`, default `gimp` (unchanged) |
  | `OPS_LadrunoMPM.cpp:438-447` (`debug shape`) | extended-lattice S | unchanged; `-folded` adds B5's S′ |
  | `LadrunoMPMPattern.cpp:330-331` (Print), `:860-862` (C1 compare), `:995` (C2 meta) | `== GIMP ? "gimp" : "linear"` — would print `linear` for a B-spline | `basisName()` |
  | `LadrunoMPMIO.cpp:399-400` (C1-BASIS read) | `gimp\|linear` | **unchanged** (K2) |
  | `LadrunoMPMIO.cpp:858` (C2 `GRID/BASIS` write) | writes `m.basis` | per K3 (§5) |
  | `LadrunoMPMImplicit.cpp:118` (trial admissibility), `:530` (supports) | via `admissibleAxis`, `particleSupport` | unchanged call; B5 folding after `:530` |
  | `LadrunoMPMImplicit.cpp:1096`, `:1142` (restart `put`/`get` of `par.basis`) | numeric | round-trips 2 and 3; `get` refuses a code outside 0–3 |
  | `ladruno_mpm_harness.py:64` (`basis` key) | passes the string | unchanged |

## 2. K4 — the boundary treatment: mirrored ghosts by images, on B-spline grids

- **B5 — decision.** On a B-spline grid, every ghost node `g` behind a face that carries `fix`
  is **slaved to its image**: per axis, `i < 0` behind a fixed min face maps to `−i`, `i > n`
  behind a fixed max face to `2n − i`; `R_g = Π_faces diag(s_a)` with `s_a = −1` for a component
  the face fixes (odd) and `+1` for a free one (even). Then `v_g = R_g v_ι(g)`. This is the method
  of images — the body plus its reflection across a symmetry (roller) plane — written as the
  constraint `u_g = R u_ι` and reduced exactly (`Tᵀ(·)T`), so it is a Galerkin system, not a
  patch. Equivalent per-component basis: `S′_{ι,a} = S_ι + s_a S_g`, `∇S′_{ι,a} = ∇S_ι + s_a ∇S_g`.
- **Why.** Zeroed ghosts (A5) force `v_{−1} = 0` where a field vanishing on the face needs
  `v_{−1} = −v_1`: for `v = kx` the quadratic weight of node −1 gives an O(1) error in `∂v/∂x`
  over a layer `h/2` wide, an O(h) displacement error everywhere, so p ≈ 1 by construction and G2
  could not pass. With images, any field odd in the fixed components and even in the free ones
  about the face is reproduced as on the infinite grid: the MMS rollers (`u_x ∝ sin 2πX` is odd
  about X = 0 and X = 1; `u_y` does not depend on X) are reproduced exactly near every face. For
  quadratic, `S′_{1,n} = N_1 − N_{−1} = ξ` on `ξ = x/h ∈ [0, ½)` — non-negative, and `S′_1 x_1 = x`.
  A clamped face (all components odd) imposes `∂²v/∂x² = 0` there: consistent, O(h²) locally.
  The open-knot / clamped basis (Steffen et al. 2008; Gan et al. 2018) was **not chosen**: its
  control points are not the vertex lattice (linear reproduction holds at Greville abscissae, not
  at `x0 + ih`), so node ids, face sets, mode-A node coordinates and Contract-2 grid metadata would
  all change. A5 zeroing kept: rejected, it fails G2 by construction.
- **B6 — mode F.** After P2G, for every touched slaved `g`: `m_ι += m_g`, `p_ι += R p_g`,
  `f_ι += R (fint_g + fext_g)` (before damping). Interface §6.4–6.5 then runs on the own nodes
  only (active test on folded masses; `g` is active iff `ι(g)` is); the face-layer fixes zero
  their components as today (`applyConditions`, `:789`). Then `v_g^{n+1} = R v_ι`, `a_g = R a_ι`.
  The second projection folds `p*` the same way, conditions, then `v*_g = R v*_ι`. G2P, `L`, `x`
  sum over the unfolded support with slaved values — identical to summing `S′`. The lumped mass of
  `ι` is `Σ m_p (S_ι + S_g)` for every component (the image's mass, always positive).
- **`fix` faces / `wall` faces / free faces.** `fix`: B5–B6. Free face: ghosts stay free nodes
  (A5, unchanged). `wall` (unilateral contact, so no image): the A5 operator on face and ghost
  layers, unchanged — first order at a wall, documented, gated nowhere. A face with both `fix`
  and `wall` is refused at `mpm analyze` on a B-spline grid.
- **GIMP and linear: nothing changes.** The fold is compiled in and run only for `BASIS_BSPLINE*`.
  Reason: interface §2's proof — the uGIMP support reaches the ghost layer only for an
  overhanging particle (`x_p − x0 < lp`), the linear hat never — and changing them would move
  frozen slice-1 and 232f values (the overhang cases). R1's GIMP legs therefore run under A5.

## 3. Mode A (particle elements)

- **B7 — node set.** Support at the step-start position: quadratic ≤ 3 per axis (≤ 9 / 27 nodes),
  cubic ≤ 4 (≤ 16 / 64); folding only merges. Ghosts behind fixed faces are **never Domain
  nodes** (folded away), so no MP constraint exists. `MAX_DOF = 81` (`LadrunoMPMPoint.cpp:52`,
  checked at `:178`) holds quadratic 3D exactly; cubic raises it to 192 (the Matrix pool too).
  Quadratic 3D supports are always 27 nodes (GIMP 8–27), so M0b's late-step memory
  (37.6–39.4 kB) is B-spline's typical one — recorded, not gated.
- **B8 — per-component S.** The element's `S`, `dS` (`LadrunoMPMPoint.h:184`, scalar per node)
  become `S[nen·nd]`, `dS[nen·nd·nd]`, used at every site (`F_rel`, G2P `:381`, FLIP `:398-401`,
  force `:479`, `K`): `F_rel,id = δ_id + Σ_a u_{a,i} ∂S′_{a,i}/∂x_d`; `x_i += Σ_a S′_{a,i} u_{a,i}`;
  `f_{a,i} = V Σ_d σ_id g′_{a,i,d} − S′_{a,i} m_p b_i`; `K = Tᵀ K T`. Lumped mass
  `M_aa = W_a m_p I` with the unsigned weight `W_ι = S_ι + S_g` (`:626`, today `S[k] * mp`) — the
  mode-F folded mass. GIMP / linear pass equal components (bitwise the shipped arithmetic).
- **B9 — order with the empty-node treatment.** Fold first; the cutoff tests folded masses;
  `renormalize` (`LadrunoMPMImplicit.cpp:563-576`) uses the unsigned `W` for its sums and scales
  every component alike: `S′ ← S′/ΣW`, `∇S′ ← (∇S′ − S′ Σ∇W)/ΣW`. With nothing dropped `ΣW = 1`,
  `Σ∇W = 0`: the identity; for GIMP and linear `W = S`: today's rule.
- **B10 — refusals.** A `gridvel` box containing a ghost behind a fixed face is refused on a
  B-spline grid. A prescribed B-spline node velocity is a control-point value, not a material
  velocity: the guide says so.
- **Why not an MP with factor −1.** `PlainHandler` ignores it: "constraint matrix not identity,
  ignoring constraint" (`PlainHandler.cpp:198`) — a warning and a wrong answer under the default
  handler. Admissible under Transformation, Penalty or Lagrange, but that would force a handler
  choice on every B-spline roller and churn MPs every step (W5: removal cleans nothing).

## 4. MMS paths: particle body force and initial velocity

- **B11 — `mpm particleBodyForce <ndm·N doubles> | -clear`.** Body force per unit **mass**
  (acceleration units, like `mpm gravity`), one flat list in storage order (= `mpm particleData`
  order). No MMS formula in C++. **Not scaled** by the host series λ: the caller owns its time
  dependence. Mode F: `fext_i += S_ip m_p b_p` beside gravity (`LadrunoMPMCore.h:533`), held in
  `Engine` so the g++ check reaches it. Mode A: `mpmBodyForce(p)` (`LadrunoMPMImplicit.cpp:99-110`,
  today `lambda * (m g + pforce)`) adds `m_p b_p` unscaled.
- **Lifecycle.** Set between `mpm analyze` calls, any number of times; length must be `ndm·N` and
  every value finite (else refused); held until replaced. `mpm analyze` refuses a stale length
  (particles added since). Cleared by `mpm reset`; gone with the host on `wipe` and
  `remove loadPattern`. Mode-A restart sends it with the host (as `pforce`, `:1110`). A cut-and-
  retried step keeps it constant over its sub-steps. Recommended use: mode F sets `b(X_p, t_n)`,
  mode A `b(X_p, t_{n+1})` (λ_{n+1} in quasi-static), `X_p` the initial position; any timing error
  is inside the brief's Δt validity check.
- **B12 — `mpm particleVelocity <ndm·N doubles>`.** An initial condition, before the first step
  (as `mpm velocity`, refused after it); kept by the `mpm reset` snapshot. `mpm velocity … sine`
  cannot serve: it zeroes the other components (`LadrunoMPMPattern.cpp:580`).
- **Tcl.** Both commands read through the shared elementAPI cursor, so classic Tcl accepts them
  (`mpm particleBodyForce {*}$b`); the per-step pattern is a Python idiom and the MMS harness is
  Python only. The guide says so.

## 5. Kernel facts, contracts, syntax

- **Kernel.** The strain feed `eps_hat` is the sum of Hughes–Winget increments in the unrotated
  frame, `G = (F₁ − Fₙ)·F_mid⁻¹`, `Δε̂ = R_midᵀ sym(G) R_mid` (interface §8.2); the material gets
  the **total** `eps_hat` (§9); `sig` is **Cauchy**, `R σ̂ Rᵀ` (§3, §6.3: "σ Cauchy"), volume
  `V = J V0`. Plane strain embeds `F` with `F33 = 1`, `ε̂_zz = 0` (§8.3); `σ_zz` is not carried.
  For the MMS, `F` is diagonal, `R = I`, and per step `Δε̂_ii = 2(F₁ − Fₙ)/(F₁ + Fₙ)`, which sums to
  `ln F_ii` + O(Δt²). `ElasticIsotropic` gives `σ = λ tr(ln F) I + 2μ ln F` (Cauchy), so
  `P = J σ F⁻ᵀ` and `b = a − (1/ρ₀) ∇_X·P` with `(∇_X·P)_i = ∂P_ii/∂X_i` (P diagonal). Round 2
  derives it by hand; `test_mpm_mms_reference.py` checks it by centred differences.
- **K2.** C1 unchanged. `mpm particles -file` into a B-spline grid is refused by the existing
  compare: the file's `GRID/BASIS` is `gimp|linear` (`LadrunoMPMIO.cpp:399`), so "the set's
  GRID/BASIS '…' differs from the deck's '…'" (`LadrunoMPMPattern.cpp:861-862`) fires once
  `basisName()` names the deck's basis. Contract 1 v2 is a follow-up for T.
- **K3 (open; costs).** (a) *apeGmsh tolerates an unknown value:* write `bspline2|bspline3`;
  amend `232a_mpm_contracts.md` §2 C2-GRID to that set for C2 only (no minor bump: an added value
  that readers pass through); split `check_grid`'s value set per contract
  (`ladruno_format.py:1197`, today `gimp|linear` for both); a C2 round-trip test. Fork only, about a
  C++ helper and 15 lines of Python. (b) *It refuses:* writing the value would make every
  B-spline `.ladruno` unopenable; `mpm output` is refused loudly on a B-spline grid until T lands a
  reader amendment; G7's `.ladruno` leg becomes "refusal tested". No G2–G5 leg reads `.ladruno`
  (the harness reads `mpm particleData`). Writing `gimp` instead is a false record: excluded.
- **Syntax.** `mpm grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp|linear|bspline2|bspline3]`;
  `mpm info basis` → 0 gimp, 1 linear, 2 bspline2, 3 bspline3; `mpm debug shape <x..> [-lp <..>]
  [-folded]` (folded: per node, id then per component `S′` and `∇S′`);
  `mpm debug shape1d <h> <lp> <r> [-basis b]`.

## 6. Proposed amendments and objections (for the orchestrator)

1. **G7 wording.** `mpm info summary` carries no basis (it returns `[mass, momentum.., max |v|,
   KE]`, `OPS_LadrunoMPM.cpp:535`); the round trip is read through `mpm info basis` instead.
2. **G1 under K4.** B5 changes DOFs, not `S`, so G1's rows are checked on the extended-lattice
   `S` exactly as written, plus a G1-K4 row on the folded `S′` within 2h of each face: odd
   component `Σ S′ (x_ι − x_face)_n = (x_p − x_face)_n`; even component `Σ S′ = 1`,
   `Σ S′ x_t = x_{p,t}`; `S′ ≥ 0` — at the G1 tolerances.
3. **Where e is measured.** At `t_end` = one period the exact `u` is 0, so `e` is a drift norm.
   The SCALE mutant (×1.5 on the G2P position increment) then scales e by about 1.5 at small
   amplitude and is expected to survive G2 and die on G3 — consistent with G6's "(or G3)" and the
   ADR §8 mutation note. Proposed: also report `e(T/4)` (peak displacement), not gated.
4. **K1.** Cubic falls out of the same path (degree parameter, `idx[4]`, one ghost layer, B7's
   `MAX_DOF`); delivered in both modes unless the orchestrator says otherwise; G5 never delays
   the quadratic gates.
